import logging
import os
import random
import webbrowser
from pathlib import Path

import folium
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
import ast

from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import find_nodes_edges_for_direction
from ..utils.general_functions import (
    define_list_upstream_downstream_edges_ids,
    remove_holes_from_polygons,
    calculate_angle,
    find_closest_edge
)


class GeneratorNetworkLumping(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basis_data: str = "0_basisdata"
    dir_inter_results: str = "1_tussenresultaat"
    dir_results: str = "2_resultaat"

    direction: str = "upstream"
    read_results: bool = False
    write_results: bool = False

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_extra: gpd.GeoDataFrame = None
    rivieren: gpd.GeoDataFrame = None
    afwateringseenheden: gpd.GeoDataFrame = None
    
    inflow_outflow_points: gpd.GeoDataFrame = None
    inflow_outflow_splits: gpd.GeoDataFrame = None

    afwateringseenheden_0: gpd.GeoDataFrame = None
    afwateringseenheden_1: gpd.GeoDataFrame = None

    inflow_outflow_edges: gpd.GeoDataFrame = None
    inflow_outflow_nodes: gpd.GeoDataFrame = None
    inflow_outflow_areas_0: gpd.GeoDataFrame = None
    inflow_outflow_areas_1: gpd.GeoDataFrame = None
    inflow_outflow_areas_2: gpd.GeoDataFrame = None
    inflow_outflow_splits_0: gpd.GeoDataFrame = None
    inflow_outflow_splits_1: gpd.GeoDataFrame = None
    inflow_outflow_splits_2: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    network_positions: dict = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: str = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.check_case_path_directory(path=self.path)
            self.read_data_from_case()


    def check_case_path_directory(self, path: Path):
        """Checks if case directory exists and if required directory structure exists

        Parameters
        ----------
        path : Path
            path to case directory. name of directory is used as case name.
            self.path and self.name are set

        Raises ValueErrors in case directory and 0_basisdata directory not exist
        """
        if not path.exists() and path.is_dir():
            raise ValueError(
                f"provided path [{path}] does not exist or is not a directory"
            )
        self.path = path
        self.name = self.path.name
        logging.info(f' ### Case "{self.name.capitalize()}" ###')

        # check if directories 0_basisdata and 1_tussenresultaat exist
        if not Path(self.path, self.dir_basis_data).exists():
            raise ValueError(f"provided path [{path}] exists but without a 0_basisdata")
        for folder in [self.dir_inter_results, self.dir_results]:
            if not Path(self.path, folder).exists():
                Path(self.path, folder).mkdir(parents=True, exist_ok=True)


    def read_data_from_case(self, path: Path = None, read_results: bool = None):
        """Read data from case: including basis data and intermediate results

        Parameters
        ----------
        path : Path, optional
            Path to the case directory including directories 0_basisdata and
            1_tussenresultaat. Directory name is used as name for the case,
            by default None
        read_results : bool, optional
            if True, it reads already all resulst from, by default None
        """
        if path is not None and path.exists():
            self.check_case_path_directory(path=path)
        logging.info(f"   x read basisdata")
        basisdata_gpkgs = [
            Path(self.path, self.dir_basis_data, f + ".gpkg")
            for f in [
                "hydroobjecten",
                "hydroobjecten_extra",
                "rivieren",
                "krw_lijn",
                "krw_vlak",
                "krw_deelgebieden",
                "afwateringseenheden",
                "inflow_outflow_points",
                "inflow_outflow_splits",
            ]
        ]
        if isinstance(read_results, bool):
            self.read_results = read_results
        baseresults_gpkgs = (
            [
                Path(self.path, self.dir_inter_results, f + ".gpkg")
                for f in [
                    "inflow_outflow_edges",
                    "inflow_outflow_nodes",
                    "inflow_outflow_areas_0",
                    "inflow_outflow_areas_1",
                    "inflow_outflow_areas_2",
                ]
            ]
            if self.read_results
            else []
        )
        for list_gpkgs in [basisdata_gpkgs, baseresults_gpkgs]:
            for x in list_gpkgs:
                if x.is_file():
                    if hasattr(self, x.stem):
                        logging.debug(f"    - get dataset {x.stem}")
                        gdf = gpd.read_file(x, layer=x.stem)
                        if "CODE" in gdf.columns:
                            gdf = gdf.rename(columns={"CODE": "code"})
                        setattr(self, x.stem, gdf)


    def create_graph_from_network(
        self, water_lines=["rivieren", "hydroobjecten", "hydroobjecten_extra"]
    ):
        """_summary_

        _extended_summary_

        Parameters
        ----------
        water_lines : list, optional
            _description_, by default ["rivieren", "hydroobjecten", "hydroobjecten_extra"]
        """
        if water_lines is None:
            water_lines = ["hydroobjecten"]
        logging.info("  x create network graph")
        self.inflow_outflow_edges = None
        for water_line in water_lines:
            gdf_water_line = getattr(self, water_line)
            if gdf_water_line is None:
                continue
            if self.inflow_outflow_edges is None:
                self.inflow_outflow_edges = gdf_water_line.explode()
            else:
                self.inflow_outflow_edges = pd.concat([
                    self.inflow_outflow_edges, 
                    gdf_water_line.explode()
                ])
        self.nodes, self.edges, self.graph = create_graph_from_edges(
            self.inflow_outflow_edges
        )
        self.network_positions = {n: [n[0], n[1]] for n in list(self.graph.nodes)}


    def find_upstream_downstream_nodes_edges(
        self, direction: str = "upstream", no_inflow_outflow_points: int = None
    ):
        """_summary_

        _extended_summary_

        Parameters
        ----------
        direction : str, optional
            _description_, by default "upstream"
        """
        if direction not in ["upstream", "downstream"]:
            raise ValueError(f" x direction needs to be 'upstream' or 'downstream'")
        self.direction = direction
        logging.info(
            f"  x find {direction} nodes and edges for {len(self.inflow_outflow_points)} outflow locations"
        )

        if no_inflow_outflow_points is not None:
            self.inflow_outflow_points = self.nodes.sample(n=no_inflow_outflow_points)

        self.inflow_outflow_points["representative_node"] = (
            self.inflow_outflow_points.geometry.apply(
                lambda x: self.nodes.geometry.distance(x).idxmin()
            )
        )

        # split_points for inflow_outflow. check if which version 2, 1 or 0 needs to be used.
        inflow_outflow_splits = None
        for i_splits, splits in enumerate([
            self.inflow_outflow_splits_2, 
            self.inflow_outflow_splits_1, 
            self.inflow_outflow_splits_0, 
            self.inflow_outflow_splits
        ]):
            if splits is not None:
                splits[["upstream_edge", "downstream_edge"]] = splits[["upstream_edge", "downstream_edge"]].replace("", None)
                for direction in ["upstream", "downstream"]:
                    if f"{direction}_edge" not in splits.columns:
                        if f"selected_{direction}_edge" in splits.columns: 
                            splits[f"{direction}_edge"] = splits[f"selected_{direction}_edge"]
                        else:
                            splits[f"{direction}_edge"] = None
                    else:
                        if f"selected_{direction}_edge" not in splits.columns: 
                            splits[f"selected_{direction}_edge"] = splits[f"{direction}_edge"]
                inflow_outflow_splits = splits
                break

        self.inflow_outflow_nodes, self.inflow_outflow_edges = find_nodes_edges_for_direction(
            nodes=self.nodes,
            edges=self.edges,
            node_ids=self.inflow_outflow_points["representative_node"].to_numpy(),
            border_node_ids=self.inflow_outflow_points["representative_node"].to_numpy(),
            direction=self.direction,
            split_points=inflow_outflow_splits
        )
        self.inflow_outflow_nodes = define_list_upstream_downstream_edges_ids(
            self.inflow_outflow_nodes.nodeID.unique(),
            self.inflow_outflow_nodes,
            self.inflow_outflow_edges,
        )
        return self.inflow_outflow_nodes


    def detect_split_points(self):
        """Detect all split points where the basins of two or more outflow/inflow points are connecting

        Returns
        -------
        self.inflow_outflow_splits_1: gpd.GeoDataFrame
            gdf with splitpoints, nodeid, downstream_edges_ids, upstream_edges_ids, etc.
        """
        logging.info("  x search for split points based on the basins of outflow/inflow points")

        inflow_outflow_nodes = self.inflow_outflow_points.representative_node.values
        if self.direction == "downstream":
            search_direction = "upstream"
            opposite_direction = "downstream"
            node_search = "node_end"
        else:
            search_direction = "downstream"
            opposite_direction = "upstream"
            node_search = "node_start"

        inflow_outflow_edges = None

        inflow_outflow_nodes = self.inflow_outflow_nodes[self.inflow_outflow_nodes.no_downstream_edges>1]

        for i_node, node in enumerate(inflow_outflow_nodes.nodeID.values):
            if i_node%50 == 0:
                logging.debug(f"    - detect points: {i_node}/{len(inflow_outflow_nodes)}")
            upstream_edges = self.inflow_outflow_edges[
                self.inflow_outflow_edges[node_search] == node
            ].copy()

            inflow_outflow_points_columns = [
                f"{self.direction}_node_{n}" for n in self.inflow_outflow_points.representative_node.values
            ]

            for col in inflow_outflow_points_columns:
                upstream_edges[col] = upstream_edges[col].replace(False, np.nan).copy()
            
            upstream_edges = upstream_edges.dropna(
                subset=inflow_outflow_points_columns, how="all"
            )

            # checked if all columns are equal
            edges = upstream_edges.drop_duplicates(subset=inflow_outflow_points_columns)
            if len(edges) > 1:
                if inflow_outflow_edges is None:
                    inflow_outflow_edges = edges.copy()
                else:
                    inflow_outflow_edges = pd.concat([inflow_outflow_edges, edges])

        for col in inflow_outflow_points_columns:
            inflow_outflow_edges[col] = (
                inflow_outflow_edges[col].replace(np.nan, False).copy()
            )

        self.inflow_outflow_splits_0 = define_list_upstream_downstream_edges_ids(
            inflow_outflow_edges[node_search].unique(),
            self.inflow_outflow_nodes,
            self.inflow_outflow_edges,
        )

        for edge in [f'{search_direction}_edge', f'{opposite_direction}_edge']:
            self.inflow_outflow_splits_0[edge] = self.inflow_outflow_splits_0.apply(
                lambda x: None if len(x[edge+'s'].split(','))>1 else x[edge+'s'], 
                axis=1
            )

        if self.inflow_outflow_splits is None:
            self.inflow_outflow_splits_1 = self.inflow_outflow_splits_0.copy()
        else:
            inflow_outflow_splits = self.inflow_outflow_splits.copy()
            for edge in [f'{search_direction}_edge', f'{opposite_direction}_edge']:
                    inflow_outflow_splits[edge] = inflow_outflow_splits.apply(
                    lambda x: x[edge] if len(x[edge+'s'].split(','))>1 else x[edge+'s'], 
                    axis=1
                )
            self.inflow_outflow_splits_1 = pd.concat([
                inflow_outflow_splits, 
                self.inflow_outflow_splits_0
            ]).reset_index(drop=True).drop_duplicates(subset="nodeID", keep='first')
        
        if self.inflow_outflow_splits is None:
            logging.debug(f"   - no. of splits as input: {0}")
        else:
            logging.debug(f"   - no. of splits as input: {len(self.inflow_outflow_splits)}")
        logging.debug(f"   - no. of splits found in network: {len(self.inflow_outflow_splits_0)}")
        logging.debug(f"   - no. of splits in total: {len(self.inflow_outflow_splits_1)}")
        
        return self.inflow_outflow_splits_1


    def export_detected_split_points(self):
        if self.inflow_outflow_splits_1 is None:
            logging.info("  x Splitsingen nog niet gevonden. gebruik functie .detect_split_points()")
        else:
            base_dir = Path(self.path, self.dir_basis_data)
            file_detected_points = "inflow_outflow_splits_detected.gpkg"
            logging.info(f"  x Split points found: saved as {self.dir_basis_data}/{file_detected_points}")
            detected_inflow_outflow_splits = self.inflow_outflow_splits_1.drop(
                columns=["selected_upstream_edge", "selected_downstream_edge"],
                errors='ignore'
            )
            detected_inflow_outflow_splits.to_file(Path(base_dir, file_detected_points))

    def calculate_angles_of_edges_at_splitpoints(self):
        self.inflow_outflow_splits_0 = self.inflow_outflow_splits.copy()

        # Initialize angle columns as empty strings
        self.inflow_outflow_splits_0['upstream_angles'] = ""
        self.inflow_outflow_splits_0['downstream_angles'] = ""

        for index, row in self.inflow_outflow_splits_0.iterrows():
            # Convert the string representation of downstream_edges to a list
            upstream_edges_str = row['upstream_edges']  # Corrected to 'upstream_edges'
            upstream_edges = ast.literal_eval(f"[{upstream_edges_str}]")
            
            downstream_edges_str = row['downstream_edges']
            downstream_edges = ast.literal_eval(f"[{downstream_edges_str}]")

            # Initialize lists to hold angles
            upstream_angles = []
            downstream_angles = []

            # Calculate angles for upstream edges
            for edge_id in upstream_edges:
                matching_lines = self.inflow_outflow_edges.loc[self.inflow_outflow_edges['code'] == str(edge_id), 'geometry']
                line = matching_lines.values[0]
                angle = calculate_angle(line, 'upstream').round(2)
                upstream_angles.append(angle)

            # Join and assign as a string
            if upstream_angles:
                self.inflow_outflow_splits_0.at[index, 'upstream_angles'] = ', '.join(map(str, upstream_angles))

            # Calculate angles for downstream edges
            for edge_id in downstream_edges:
                matching_lines = self.inflow_outflow_edges.loc[self.inflow_outflow_edges['code'] == str(edge_id), 'geometry']
                line = matching_lines.values[0]
                angle = calculate_angle(line, 'downstream').round(2)
                downstream_angles.append(angle)

            # Join and assign as a string
            if downstream_angles:
                self.inflow_outflow_splits_0.at[index, 'downstream_angles'] = ', '.join(map(str, downstream_angles))

        return self.inflow_outflow_splits_0
    
    
    def select_directions_for_splits_based_on_angle(self):
        self.inflow_outflow_splits_1 = self.inflow_outflow_splits_0.copy()

        for index, row in self.inflow_outflow_splits_1.iterrows():
            if self.direction == 'upstream':
                # Get upstream angles and edges
                upstream_angles = row['upstream_angles']
                print(upstream_angles)
                upstream_edges = row['upstream_edges']
                downstream_angles_str = row['downstream_angles']
                downstream_edges_str = row['downstream_edges']
                # Convert strings to lists
                downstream_edges_list = downstream_edges_str.split(',')
                downstream_angles_list = downstream_angles_str.split(',')
                
                # Assuming there's a reference angle, e.g., the first angle in the list
                reference_angle = upstream_angles
                
                selected_edge = find_closest_edge(reference_angle, downstream_angles_list, downstream_edges_list)
                
                # Update the selected columns
                self.inflow_outflow_splits_1.at[index, 'selected_downstream_edge'] = selected_edge

            elif self.direction == 'downstream':
                # Get downstream angles and edges
                downstream_angles = row['downstream_angles']
                downstream_edges = row['downstream_edges']
                upstream_angles_str = row['upstream_angles']
                upstream_edges_str = row['upstream_edges']
                upstream_edges_list = upstream_edges_str.split(',')
                upstream_angles_list = upstream_angles_str.split(',')
                
                # Assuming there's a reference angle, e.g., the first angle in the list
                reference_angle = downstream_angles
                
                selected_edge = find_closest_edge(reference_angle, upstream_angles_list, upstream_edges_list)
                
                # Update the selected columns
                self.inflow_outflow_splits_1.at[index, 'selected_upstream_edge'] = selected_edge

        return self.inflow_outflow_splits_1

    def select_directions_for_splits(self, fillna_with_random=False):
        """_summary_: This function can be used to define the random direction at split points. 
        
        inflow_outflow_splits_0 is the GeoDataFrame with the detected points, 
        without decisions. These decisions can be made manually in GIS, 
        by adding the correct hydroobject in the column {search_direction}_edges, 
        and saving the file as inflow_outflow_splits_1. 
        
        This function fills the entire column when there is no inflow_outflow_splits_1
        GeoDataFrame. When it is present, the function only fills empty instances in the column.
        """
        # check whether to use inflow_outflow_splits_1 or inflow_outflow_splits_0
        if self.inflow_outflow_splits_1 is not None and not self.inflow_outflow_splits_1.empty:
            self.inflow_outflow_splits_2 = self.inflow_outflow_splits_1.copy()
        elif self.inflow_outflow_splits_0 is not None and not self.inflow_outflow_splits_0.empty:
            self.inflow_outflow_splits_2 = self.inflow_outflow_splits_0.copy()
        else:
            logging.debug(" x no splits found: no direction for splits selected")
            return None

        logging.info(" x search for direction in splits")
        for search_direction in ["upstream", "downstream"]:
            no_splits_known = len(self.inflow_outflow_splits_2[
                ~self.inflow_outflow_splits_2[f'{search_direction}_edge'].isna()
            ])
            logging.debug(f"   - known {search_direction} direction at splits: {no_splits_known}/{len(self.inflow_outflow_splits_2)}")
            self.inflow_outflow_splits_2[f'selected_{search_direction}_edge'] = (
                self.inflow_outflow_splits_2.apply(
                    lambda x: random.choice(x[f"{search_direction}_edges"].split(','))
                    if x[f'{search_direction}_edge'] is None
                    else x[f'{search_direction}_edge'],
                    axis=1,
                )
            )
            logging.debug(f"   - randomly choosen {search_direction} direction at splits: {
                len(self.inflow_outflow_splits_2) - no_splits_known}/{len(self.inflow_outflow_splits_2)}")
        return self.inflow_outflow_splits_2


    def assign_drainage_units_to_outflow_points_based_on_id(self):
        self.afwateringseenheden["gridcode"] = (
            self.afwateringseenheden["gridcode"].round(0).astype("Int64").astype(str)
        )
        self.inflow_outflow_edges["code"] = self.inflow_outflow_edges["code"].astype(str)

        upstream_downstream_columns = [
            f"{self.direction}_node_{node}"
            for node in self.inflow_outflow_points["representative_node"].tolist()
        ]
        self.afwateringseenheden_0 = self.afwateringseenheden.merge(
            self.inflow_outflow_edges[
                ["code"] + [f"{column}" for column in upstream_downstream_columns]
            ],
            how="left",
            left_on="gridcode",
            right_on="code",
        )
        self.afwateringseenheden_0 = self.afwateringseenheden_0.drop(columns=["code"])
        self.afwateringseenheden_0[upstream_downstream_columns] = self.afwateringseenheden_0[upstream_downstream_columns].fillna(False)


    def assign_drainage_units_to_outflow_points_based_on_length_hydroobject(self):
        if self.afwateringseenheden is None:
            return None
        self.afwateringseenheden["unique_id"] = self.afwateringseenheden.index
        self.afwateringseenheden["savedgeom"] = self.afwateringseenheden.geometry

        joined = gpd.sjoin(
            self.inflow_outflow_edges.rename(columns={"code": "code_hydroobject"}),
            self.afwateringseenheden,
            how="inner",
            predicate="intersects",
        )

        joined["intersection_length"] = joined.apply(
            lambda row: row.geometry.intersection(row.savedgeom).length, axis=1
        )

        merged = self.afwateringseenheden.merge(
            joined[["unique_id", "code_hydroobject", "intersection_length"]],
            on="unique_id",
            how="inner",
        )


        # Select the rows from the merged GeoDataFrame that correspond to those indices
        max_intersections = merged.groupby("unique_id")["intersection_length"].idxmax()
        result = merged.loc[max_intersections]
        result = result.drop(columns=["savedgeom"]).reset_index(drop=True)

        upstream_columns = [
            f"{self.direction}_node_{node}"
            for node in self.inflow_outflow_points["representative_node"].tolist()
        ]
        self.afwateringseenheden_1 = result.merge(
            self.inflow_outflow_edges[
                ["code"] + [f"{column}" for column in upstream_columns]
            ],
            how="left",
            left_on="code_hydroobject",
            right_on="code",
        ).reset_index(drop=True)
        self.afwateringseenheden_1 = self.afwateringseenheden_1.loc[
            :, ~self.afwateringseenheden_1.columns.duplicated()
        ].copy()

        for col in upstream_columns:
            self.afwateringseenheden_1[col] = self.afwateringseenheden_1[col].fillna(False)
        return self.afwateringseenheden_1


    def dissolve_assigned_drainage_units(self):
        if self.afwateringseenheden_1 is None:
            return None
        self.inflow_outflow_areas_0 = None
        for node in list(self.inflow_outflow_points["representative_node"].unique()):
            filtered_areas = self.afwateringseenheden_1[
                self.afwateringseenheden_1[f"{self.direction}_node_{node}"]
            ]
            # Step 2: Dissolve the filtered geometries
            dissolved_areas = filtered_areas[["geometry"]].dissolve().explode()
            dissolved_areas["inflow_outflow_point"] = node
            dissolved_areas["area"] = dissolved_areas.geometry.area

            if self.inflow_outflow_areas_0 is None:
                self.inflow_outflow_areas_0 = dissolved_areas.reset_index(drop=True)
            else:
                self.inflow_outflow_areas_0 = pd.concat([
                    self.inflow_outflow_areas_0, 
                    dissolved_areas
                ]).reset_index(drop=True)

        self.inflow_outflow_areas_0.geometry = remove_holes_from_polygons(
            self.inflow_outflow_areas_0.geometry, min_area=50
        )
        self.inflow_outflow_areas_0.geometry = self.inflow_outflow_areas_0.geometry.buffer(0.1)
        return self.inflow_outflow_areas_0


    def export_results_all(
        self,
        html_file_name: str = None,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
    ):
        """Export results to geopackages and folium html"""
        self.export_results_to_gpkg()
        self.export_results_to_html_file(
            html_file_name=html_file_name,
            width_edges=width_edges,
            opacity_edges=opacity_edges,
        )


    def export_results_to_gpkg(self):
        """Export results to geopackages in folder 1_tussenresultaat"""
        results_dir = Path(self.path, self.dir_inter_results)
        logging.info(f"  x export results")
        for layer in [
            "inflow_outflow_points",
            "inflow_outflow_edges",
            "inflow_outflow_nodes",
            "inflow_outflow_areas_0",
            "inflow_outflow_areas_1",
            "inflow_outflow_areas_2",
            "afwateringseenheden_0",
            "afwateringseenheden_1",
        ]:
            result = getattr(self, layer)
            if result is not None:
                logging.debug(f"   - {layer}")
                result.to_file(Path(results_dir, f"{layer}.gpkg"))


    def export_results_to_html_file(
        self,
        html_file_name: str = None,
        include_areas: bool = True,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
        open_html: bool = False
    ):
        """Export results to folium html file

        Parameters
        ----------
        html_file_name : str, optional
            filename folium html, by default None
        width_edges : float, optional
            width (meters) of edges in folium html, by default 10.0
        opacity_edges : float, optional
            opacity of edges in folium html, by default 0.5
        """
        logging.info(f'   x saving html file')

        nodes_selection = self.inflow_outflow_points.representative_node.to_numpy()
        no_nodes = len(self.inflow_outflow_points) + 1
        nodes_colors = plt.get_cmap("hsv", no_nodes)
        i_nodes_colors = np.arange(start=0, stop=no_nodes-1)
        np.random.shuffle(i_nodes_colors)
        nodes_colors = [nodes_colors(i) for i in i_nodes_colors]
        nodes_4326 = self.inflow_outflow_nodes.to_crs(4326)

        m = folium.Map(
            location=[nodes_4326.geometry.y.mean(), nodes_4326.geometry.x.mean()],
            tiles=None,
            zoom_start=12,
        )

        fg = folium.FeatureGroup(name=f"Watergangen", control=True).add_to(m)

        folium.GeoJson(
            self.inflow_outflow_edges.buffer(width_edges / 2.0),
            color="grey",
            weight=5,
            z_index=0,
            opacity=0.25,
        ).add_to(fg)

        folium.GeoJson(
            self.inflow_outflow_nodes,
            marker=folium.Circle(
                radius=width_edges,
                fill_color="darkgrey",
                fill_opacity=0.5,
                color="darkgrey",
                weight=1,
                z_index=1,
            ),
        ).add_to(fg)

        if self.afwateringseenheden is not None and include_areas:
            folium.GeoJson(
                self.afwateringseenheden[['geometry']].explode(ignore_index=True),
                fill_opacity=0.0,
                color="grey",
                weight=0.5,
                z_index=10,
                name="Afwateringseenheden",
            ).add_to(m)

        for inflow_outflow_areas in [
            self.inflow_outflow_areas_2, 
            self.inflow_outflow_areas_1, 
            self.inflow_outflow_areas_0
        ]:
            if inflow_outflow_areas is not None:
                break

        inflow_outflow = "instroom" if self.direction == "downstream" else "uitstroom"

        fg = folium.FeatureGroup(
            name=f"{inflow_outflow.capitalize()}punten", 
            control=True, 
            show=True,
            z_index=0
        ).add_to(m)
        for i, node_selection in enumerate(nodes_selection):
            inflow_outflow = "instroom" if self.direction == "downstream" else "uitstroom"
            c = matplotlib.colors.rgb2hex(nodes_colors[i])
        
            folium.GeoJson(
                self.inflow_outflow_nodes.iloc[[node_selection]],
                marker=folium.Circle(
                    radius=width_edges * 7.5,
                    fill_color=c,
                    fill_opacity=1.0,
                    color=c,
                    opacity=1.0,
                    weight=4,
                    z_index=10,
                ),
            ).add_to(fg)

        for i, node_selection in enumerate(nodes_selection):
            c = matplotlib.colors.rgb2hex(nodes_colors[i])
            fg = folium.FeatureGroup(
                name=f"Gebied {inflow_outflow}punt {node_selection}", 
                control=True, 
                show=True,
                z_index=2
            ).add_to(m)

            sel_inflow_outflow_edges = self.inflow_outflow_edges[
                self.inflow_outflow_edges[f"{self.direction}_node_{node_selection}"]
            ].copy()

            sel_inflow_outflow_edges.geometry = sel_inflow_outflow_edges.buffer(width_edges / 2.0)
            if len(sel_inflow_outflow_edges) > 0:
                folium.GeoJson(
                    sel_inflow_outflow_edges[["geometry"]],
                    color=c,
                    weight=5,
                    z_index=2,
                    opacity=opacity_edges,
                    fill_opacity=opacity_edges,
                ).add_to(fg)

            if inflow_outflow_areas is not None and include_areas:
                inflow_outflow_areas_node = inflow_outflow_areas[inflow_outflow_areas[f"inflow_outflow_point"]==node_selection].copy()
                folium.GeoJson(
                    inflow_outflow_areas_node,
                    color=c,
                    weight=1,
                    z_index=2,
                    opacity=opacity_edges * 0.5,
                    fill_opacity=opacity_edges * 0.5,
                ).add_to(fg)

            folium.GeoJson(
                self.inflow_outflow_nodes.iloc[[node_selection]],
                marker=folium.Circle(
                    radius=width_edges * 7.5,
                    fill_color=c,
                    fill_opacity=1.0,
                    color=c,
                    opacity=1.0,
                    weight=4,
                    z_index=0,
                ),
            ).add_to(fg)
        
        # Voorgedefinieerde splits
        if self.inflow_outflow_splits is not None:
            folium.GeoJson(
                self.inflow_outflow_splits.loc[
                    self.inflow_outflow_splits[
                        ["upstream_edge", "downstream_edge"]
                    ].dropna(how='all').index
                ],
                marker=folium.Circle(
                    radius=width_edges * 7.5,
                    fill_color="black",
                    fill_opacity=0.1,
                    color="black",
                    weight=3,
                    z_index=10,
                ),
                name=f"Voorgedefinieerde splitsingen",
                show=True,
            ).add_to(m)

        # Voorgedefinieerde splits
        if self.inflow_outflow_splits_1 is not None:
            search_direction = "upstream" if self.direction == "downstream" else "downstream"
            folium.GeoJson(
                self.inflow_outflow_splits_1.loc[
                    self.inflow_outflow_splits_1[f"{search_direction}_edge"].isna()
                ],
                marker=folium.Circle(
                    radius=width_edges * 7.5,
                    fill_color="red",
                    fill_opacity=0.1,
                    color="red",
                    weight=3,
                    z_index=1,
                ),
                name=f"Extra gevonden splitsingen",
                show=True,
            ).add_to(m)

        folium.TileLayer(
            "openstreetmap", 
            name="Open Street Map", 
            show=False
        ).add_to(m)
        folium.TileLayer(
            "cartodbpositron", 
            name="Light Background", 
            show=True
        ).add_to(m)

        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m

        if html_file_name is None:
            html_file_name = self.name
        
        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f' x html file saved: {html_file_name}.html')

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m
        
