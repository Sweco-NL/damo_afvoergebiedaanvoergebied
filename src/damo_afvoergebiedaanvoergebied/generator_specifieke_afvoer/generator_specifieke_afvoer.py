import logging
import time
from pathlib import Path

import folium
import geopandas as gpd
import pandas as pd
import imod
import numpy as np
import pandas as pd
import pyflwdir
import rioxarray
import xarray
import networkx as nx
import rasterstats
import xarray as xr
from geocube.api.core import make_geocube
from pydantic import ConfigDict
from tqdm import tqdm

from ..generator_basis import GeneratorBasis
from ..utils.folium_map import generate_folium_map
from ..utils.network_functions import sum_edge_node_values_through_network


class GeneratorSpecifiekeAfvoer(GeneratorBasis):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str = "1_resultaat"
    waterschap: str = None

    read_results: bool = False
    write_results: bool = False
    crs: int = 28992

    water_lines: list[str] = ["hydroobject"]
    
    required_results: list[str] = [
        "hydroobject",
        "hydroobject_processed_0", 
        "hydroobject_processed_1", 
        "overige_watergang_processed_4", 
        "potential_culverts_5", 
        "outflow_nodes",
        "edges", 
        "nodes",
    ]

    hydroobject: gpd.GeoDataFrame = None
    hydroobject_processed_0: gpd.GeoDataFrame = None
    hydroobject_processed_1: gpd.GeoDataFrame = None
    
    overige_watergang: gpd.GeoDataFrame = None
    overige_watergang_processed_4: gpd.GeoDataFrame = None
    potential_culverts_5: gpd.GeoDataFrame = None
    afvoergebied_gdf: gpd.GeoDataFrame = None

    snapping_distance: float = 0.05
    use_specifieke_afvoer: float = 1.0

    outflow_nodes: gpd.GeoDataFrame = None
    splitsing: gpd.GeoDataFrame = None
    split_nodes: gpd.GeoDataFrame = None

    file_name_specifieke_afvoer: str = None
    specifieke_afvoer: xarray.Dataset = None
    specifieke_afvoer_filled: xarray.Dataset = None
    
    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None
    
    folium_map: folium.Map = None
    folium_html_path: str = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.use_processed_hydroobject(force_preprocess=False)
        # if self.edges is None:
        self.create_graph_from_network(water_lines=self.water_lines)
        self.analyse_netwerk_add_information_to_nodes_edges(min_difference_angle=20.0)


    def generate_distribution_splits_downstream(self):
        self.split_nodes = self.nodes[
            (self.nodes["no_downstream_edges"] > 1) &
            (self.nodes["no_upstream_edges"] > 0)
        ].copy()
        self.split_nodes["downstream_splits_dist"] = self.split_nodes["no_downstream_edges"].apply(
            lambda x: ",".join([str(1/x) for i in range(x)])
        ).values.tolist()

        self.split_nodes = self.split_nodes.sjoin(
            self.splitsing, 
            how="left", 
            predicate="dwithin",
            distance=5.0,
            lsuffix=None, 
            rsuffix='x'
        )

        def replace_downstream_split_dist(row, max_angle_difference=20.0):
            if pd.isna(row["downstream_splits_dist_x"]) or \
                pd.isna(row["downstream_angles_x"]):
                return row

            old_angles = [float(i) for i in row["downstream_angles"].split(",")]
            old_dist = [float(i) for i in row["downstream_splits_dist"].split(",")]

            new_angles = [float(i) for i in row["downstream_angles_x"].split(",")]
            new_dist = [float(i) for i in row["downstream_splits_dist_x"].split(",")]

            diff_angles_1 = [abs(a1 - a2) % 360.0 for a1, a2 in zip(old_angles, new_angles)]
            diff_angles_2 = [abs(a1 - a2) % 360.0 for a1, a2 in zip(old_angles, new_angles[::-1])]

            if sum(diff_angles_1) < sum(diff_angles_2):
                row["downstream_splits_dist"] = ",".join([str(i) for i in new_dist])
            else:
                row["downstream_splits_dist"] = ",".join([str(i) for i in new_dist[::-1]])
            return row

        self.split_nodes = self.split_nodes.apply(lambda x: replace_downstream_split_dist(x), axis=1)
        self.split_nodes = self.split_nodes.drop(columns=[
            "upstream_angles_x",
            "downstream_angles_x",
            "downstream_splits_dist_x",
        ])

        downstream_splits = self.split_nodes[["downstream_edges", "downstream_splits_dist"]].copy()
        downstream_splits["downstream_edges"] = downstream_splits["downstream_edges"].apply(lambda x: x.split(","))
        downstream_splits["downstream_splits_dist"] = downstream_splits["downstream_splits_dist"].apply(lambda x: [float(i) if i!='' else None for i in x.split(",")])
        downstream_splits = downstream_splits.explode(["downstream_edges", "downstream_splits_dist"])
        downstream_splits = downstream_splits[downstream_splits.downstream_splits_dist.notnull()]
        downstream_splits.rename(columns={"downstream_edges": "downstream_edge"}, inplace=True)
        self.edges = self.edges.drop(
            columns=["downstream_splits_dist", "downstream_edge"],
            errors="ignore"
        ).merge(
            downstream_splits,
            left_on="code", 
            right_on="downstream_edge", 
            how="left"
        )
        self.edges["downstream_splits_dist"] = self.edges["downstream_splits_dist"].fillna(1.0)


    def read_specifieke_afvoer(self, file_name_specifieke_afvoer: str = None):
        """Read specific discharge data from a file."""
        if file_name_specifieke_afvoer is None:
            return None
        logging.info("   x read topographical data as input")
        self.file_name_specifieke_afvoer = file_name_specifieke_afvoer
        self.specifieke_afvoer = xr.open_dataset(Path(self.path, self.dir_basisdata, file_name_specifieke_afvoer))
        self.specifieke_afvoer = self.specifieke_afvoer.rename_vars({"__xarray_dataarray_variable__": "specifieke_afvoer"})
        if self.specifieke_afvoer.rio.crs is None:
            self.specifieke_afvoer = self.specifieke_afvoer.rio.write_crs(28992)
        # interpolate
        self.specifieke_afvoer = self.specifieke_afvoer.interpolate_na(dim="x", method="nearest")
        return self.specifieke_afvoer
    

    def add_specifieke_afvoer_to_discharge_units(self, level_discharge_units=0, use_specifieke_afvoer=0):
        """Add specific discharge to discharge units"""
        logging.info(
            f"   x add specific discharge to discharge units"
        )
        if 0 <= level_discharge_units <= 4:
            discharge_units_file_path = Path(self.dir_results, f"afvoergebied_{level_discharge_units}_gdf.gpkg")
            if discharge_units_file_path.exists():
                self.afvoergebied_gdf = gpd.read_file(discharge_units_file_path)
                area_afvoergebied = self.afvoergebied_gdf.geometry.area.sum()/10000.0
                logging.info(
                    f"     - afvoergebied level{level_discharge_units} [{len(self.afvoergebied_gdf)}] with area {round(area_afvoergebied, 2)}[ha]"
                )

        if self.afvoergebied_gdf is None or use_specifieke_afvoer<0.0:
            logging.info("     - no drainage units or specific discharge defined")
            return self.afvoergebied_gdf

        if use_specifieke_afvoer == 0:
            # TODO: link this to specific discharge area discharge units
            logging.info("     - add distributed specific discharge to afvoergebied [l/s]")
            afvoergebied_geojson = rasterstats.zonal_stats(
                self.afvoergebied_gdf,
                self.specifieke_afvoer["specifieke_afvoer"].data[0],
                nodata=np.nan,
                affine=self.specifieke_afvoer.rio.transform(),
                stats=["sum", "mean"],
                geojson_out=True
            )
            afvoergebied_gdf = gpd.GeoDataFrame.from_features(afvoergebied_geojson)
            self.afvoergebied_gdf = afvoergebied_gdf.rename(
                columns={"sum": "specifieke_afvoer", "mean": "mean_specifieke_afvoer"}
            )
            self.afvoergebied_gdf["specifieke_afvoer"] = self.afvoergebied_gdf["specifieke_afvoer"]*25*25/1000/24/3600
            self.afvoergebied_gdf["specifieke_afvoer"] = self.afvoergebied_gdf["specifieke_afvoer"].fillna(0.0)
            
        elif use_specifieke_afvoer > 0:
            logging.info(f"     - add homogenic specific discharge {round(use_specifieke_afvoer, 2)} [l/ha/s] to drainage units")
            self.afvoergebied_gdf["specifieke_afvoer"] = self.afvoergebied_gdf.geometry.area / 10000.0 * use_specifieke_afvoer

        total_specifieke_afvoer = self.afvoergebied_gdf["specifieke_afvoer"].sum() / 1000.0
        logging.info(f"     - total specific discharge: {round(total_specifieke_afvoer, 2)} [m3/s]")
        return self.afvoergebied_gdf


    def add_specifieke_afvoer_to_edges(self):
        """Specify specific discharge to edges and nodes."""
        logging.info("   x add specific discharge to edges and nodes")
        # TODO: link this to specific discharge discharge units
        if self.afvoergebied_gdf is None:
            self.edges["specifieke_afvoer"] = self.edges.geometry.length / 1000.0
            self.nodes["specifieke_afvoer"] = 0.0
            total_specifieke_afvoer = self.edges["specifieke_afvoer"].sum()
            logging.info(
                f"     - no afvoergebied, total length edges used: {total_specifieke_afvoer} [km]"
            )
        else:
            specifieke_afvoer = self.edges[["code"]].merge(
                self.afvoergebied_gdf[["code", "specifieke_afvoer"]], 
                how="left", 
                left_on="code", 
                right_on="code"
            )
            self.edges = self.edges.merge(
                specifieke_afvoer,
                how="left", 
                left_on="code", 
                right_on="code"
            )
            self.edges["specifieke_afvoer"] = self.edges["specifieke_afvoer"].fillna(0.0)
            self.nodes["specifieke_afvoer"] = 0.0
            total_specifieke_afvoer = self.edges["specifieke_afvoer"].sum() / 1000.0
            logging.info(
                f"     - total discharge added to edges: {total_specifieke_afvoer} [m3/s]"
            )
        return self.edges, self.nodes
    

    def sum_specifieke_afvoer_through_network(self):
        """Sum specific discharge through the network."""
        logging.info("   x sum specific discharge through network")
        edges = self.edges.copy()
        nodes = self.nodes.copy()
        
        edges, nodes = sum_edge_node_values_through_network(
            edges=edges, 
            nodes=nodes, 
            edges_id_column="code",
            nodes_id_column="nodeID",
            column_to_sum="specifieke_afvoer", 
            sum_column="total_specifieke_afvoer"
        )
        edges["log10_total_specifieke_afvoer"] = edges["total_specifieke_afvoer"].apply(lambda x: np.log10(x) if x>0.0 else np.nan)
        
        self.edges = edges.copy()
        self.nodes = nodes.copy()
        
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "edges",
                "nodes",
                "split_nodes",
            ])
        return self.edges, self.nodes


    def generate_folium_map(self, html_file_name:str="", **kwargs):
        if html_file_name == '':
            html_file_name = self.name + "_specifieke_afvoer"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            **kwargs
        )
        return self.folium_map

