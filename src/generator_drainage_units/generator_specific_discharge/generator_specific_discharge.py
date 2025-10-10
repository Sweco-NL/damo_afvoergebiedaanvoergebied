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
from geocube.api.core import make_geocube
from pydantic import ConfigDict
from tqdm import tqdm

from ..generator_basis import GeneratorBasis
from ..utils.folium_map import generate_folium_map
from ..utils.network_functions import sum_edge_node_values_through_network


class GeneratorSpecificDischarge(GeneratorBasis):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str = "1_resultaat"
    waterschap: str = None

    read_results: bool = False
    write_results: bool = False
    crs: int = 28992

    water_lines: list[str] = ["hydroobjecten"]
    
    required_results: list[str] = [
        "hydroobjecten",
        "hydroobjecten_processed_0", 
        "hydroobjecten_processed_1", 
        "overige_watergangen_processed_4", 
        "potential_culverts_5", 
        "outflow_nodes",
        "edges", 
        "nodes",
    ]

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    hydroobjecten_processed_1: gpd.GeoDataFrame = None
    
    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed_4: gpd.GeoDataFrame = None
    potential_culverts_5: gpd.GeoDataFrame = None
    drainage_units_gdf: gpd.GeoDataFrame = None

    snapping_distance: float = 0.05
    use_specific_discharge: float = 1.0

    outflow_nodes: gpd.GeoDataFrame = None
    defined_split_nodes: gpd.GeoDataFrame = None
    split_nodes: gpd.GeoDataFrame = None

    specific_discharge_file_name: str = None
    specific_discharge: xarray.Dataset = None
    specific_discharge_filled: xarray.Dataset = None
    
    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None
    
    folium_map: folium.Map = None
    folium_html_path: str = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.use_processed_hydroobjecten(force_preprocess=False)
        # if self.edges is None:
        self.create_graph_from_network(water_lines=self.water_lines)
        self.analyse_netwerk_add_information_to_nodes_edges(min_difference_angle=20.0)


    def generate_distribution_splits_downstream(self):
        self.nodes["downstream_splits_dist"] = \
            self.nodes["no_downstream_edges"].apply(
                lambda x: ",".join([str(1/x) for i in range(x)])
            ).values.tolist()
        downstream_splits = self.nodes[["downstream_edges", "downstream_splits_dist"]].copy()
        downstream_splits["downstream_edges"] = downstream_splits["downstream_edges"].apply(lambda x: x.split(","))
        downstream_splits["downstream_splits_dist"] = downstream_splits["downstream_splits_dist"].apply(lambda x: [float(i) if i!='' else None for i in x.split(",")])
        downstream_splits = downstream_splits.explode(["downstream_edges", "downstream_splits_dist"])
        downstream_splits = downstream_splits[downstream_splits.downstream_splits_dist.notnull()]
        downstream_splits.rename(columns={"downstream_edges": "downstream_edge"}, inplace=True)
        self.edges = self.edges.merge(
            downstream_splits, 
            left_on="code", 
            right_on="downstream_edge", 
            how="left"
        )


    def read_specific_discharge(self, specific_discharge_file_name: str = None):
        """Read specific discharge data from a file."""
        if specific_discharge_file_name is None:
            return None
        logging.info("   x read topographical data as input")
        self.specific_discharge_file_name = specific_discharge_file_name
        self.specific_discharge = rioxarray.open_rasterio(Path(self.path, self.dir_basisdata, specific_discharge_file_name))
        self.specific_discharge.name = "specific_discharge"
        if self.specific_discharge.rio.crs is None:
            self.specific_discharge = self.specific_discharge.rio.write_crs(28992)
        # interpolate
        self.specific_discharge = self.specific_discharge.interpolate_na(dim="x", method="nearest")
        return self.specific_discharge
    

    def add_specific_discharge_to_discharge_units(self, use_specific_discharge=0, level_discharge_units=0):
        """Add specific discharge to discharge units"""
        logging.info(
            f"   x add specific discharge to discharge units"
        )
        if 0 <= level_discharge_units <= 3:
            discharge_units_file_path = Path(self.dir_results, f"drainage_units_{level_discharge_units}_gdf.gpkg")
            if discharge_units_file_path.exists():
                self.drainage_units_gdf = gpd.read_file(discharge_units_file_path)
                area_drainage_units = self.drainage_units_gdf.geometry.area.sum()/10000.0
                logging.info(
                    f"     - drainage_units level{level_discharge_units} [{len(self.drainage_units_gdf)}] with area {round(area_drainage_units, 2)}[ha]"
                )

        if self.drainage_units_gdf is None or use_specific_discharge<0.0:
            logging.info("     - no drainage units or specific discharge defined")
            return self.drainage_units_gdf

        if use_specific_discharge == 0:
            # TODO: link this to specific discharge area discharge units
            logging.info("     - add distributed specific discharge to drainage_units [l/s]")
            drainage_units_geojson = rasterstats.zonal_stats(
                self.drainage_units_gdf,
                self.specific_discharge.data[0],
                nodata=self.specific_discharge.rio.nodata,
                affine=self.specific_discharge.rio.transform(),
                stats=["sum"],
                geojson_out=True
            )
            drainage_units_gdf = gpd.GeoDataFrame.from_features(drainage_units_geojson)
            self.drainage_units_gdf = drainage_units_gdf.rename(columns={"sum": "specific_discharge"})
            self.drainage_units_gdf["specific_discharge"] = self.drainage_units_gdf["specific_discharge"]*25*25/1000/24/3600
            self.drainage_units_gdf["specific_discharge"] = self.drainage_units_gdf["specific_discharge"].fillna(0.0)
            
        elif use_specific_discharge > 0:
            logging.info(f"     - add homogenic specific discharge {round(use_specific_discharge, 2)} [l/ha/s] to drainage units")
            self.drainage_units_gdf["specific_discharge"] = self.drainage_units_gdf.geometry.area / 10000.0 * use_specific_discharge

        total_specific_discharge = self.drainage_units_gdf["specific_discharge"].sum() / 1000.0
        logging.info(f"     - total specific discharge: {round(total_specific_discharge, 2)} [m3/s]")
        return self.drainage_units_gdf


    def add_specific_discharge_to_edges(self):
        """Specify specific discharge to edges and nodes."""
        logging.info("   x add specific discharge to edges and nodes")
        # TODO: link this to specific discharge discharge units
        if self.drainage_units_gdf is None:
            self.edges["specific_discharge"] = self.edges.geometry.length / 1000.0
            self.nodes["specific_discharge"] = 0.0
            total_specific_discharge = self.edges["specific_discharge"].sum()
            logging.info(
                f"     - no drainage_units, total length edges used: {total_specific_discharge} [km]"
            )
        else:
            specific_discharge = self.edges[["code"]].merge(
                self.drainage_units_gdf[["code", "specific_discharge"]], 
                how="left", 
                left_on="code", 
                right_on="code"
            )
            self.edges = self.edges.merge(
                specific_discharge,
                how="left", 
                left_on="code", 
                right_on="code"
            )
            self.edges["specific_discharge"] = self.edges["specific_discharge"].fillna(0.0)
            self.nodes["specific_discharge"] = 0.0
            total_specific_discharge = self.edges["specific_discharge"].sum() / 1000.0
            logging.info(
                f"     - total discharge added to edges: {total_specific_discharge} [m3/s]"
            )
        return self.edges, self.nodes
    

    def sum_specific_discharge_through_network(self):
        """Sum specific discharge through the network."""
        logging.info("   x sum specific discharge through network")
        edges = self.edges.copy()
        nodes = self.nodes.copy()
        
        edges, nodes = sum_edge_node_values_through_network(
            edges=edges, 
            nodes=nodes, 
            edges_id_column="code",
            nodes_id_column="nodeID",
            column_to_sum="specific_discharge", 
            sum_column="total_specific_discharge"
        )
        edges["log10_total_specific_discharge"] = edges["total_specific_discharge"].apply(lambda x: np.log10(x) if x>0.0 else np.nan)
        self.edges = edges.copy()
        self.nodes = nodes.copy()
        
        self.split_nodes = self.nodes[
            (self.nodes["no_downstream_edges"] > 1) &
            (self.nodes["no_upstream_edges"] > 0)
        ].copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "edges",
                "nodes",
                "split_nodes",
            ])
        return self.edges, self.nodes


    def generate_folium_map(self, html_file_name:str="", **kwargs):
        if html_file_name == '':
            html_file_name = self.name + "_specific_discharge"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            **kwargs
        )
        return self.folium_map

