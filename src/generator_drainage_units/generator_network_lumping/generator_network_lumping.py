import ast
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

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.folium_utils import add_basemaps_to_folium_map
from ..utils.general_functions import (
    calculate_angles_of_edges_at_nodes,
    define_list_upstream_downstream_edges_ids,
    find_edge_smallest_angle_difference,
    remove_holes_from_polygons,
    remove_z_dims,
)
from ..utils.network_functions import find_nodes_edges_for_direction


class GeneratorNetworkLumping(GeneratorBasis):
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

    inflow_outflow_edges: gpd.GeoDataFrame = None
    inflow_outflow_nodes: gpd.GeoDataFrame = None
    inflow_outflow_areas: gpd.GeoDataFrame = None

    inflow_outflow_splits_0: gpd.GeoDataFrame = None
    inflow_outflow_splits_1: gpd.GeoDataFrame = None
    inflow_outflow_splits_2: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    network_positions: dict = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: str = None


    def create_graph_from_network(
        self, water_lines=["rivieren", "hydroobjecten", "hydroobjecten_extra"]
    ):
        if water_lines is None:
            water_lines = ["hydroobjecten"]
        logging.info("   x create network graph")
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
        if direction not in ["upstream", "downstream"]:
            raise ValueError(f" x direction needs to be 'upstream' or 'downstream'")
        self.direction = direction

        if self.inflow_outflow_points is None:
            logging.info(
                f"   x no outflow locations available"
            )
            return None

        logging.info(
            f"   x find {direction} nodes and edges for {len(self.inflow_outflow_points)} outflow locations"
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
        for i_splits, splits in enumerate(
            [
                self.inflow_outflow_splits_2,
                self.inflow_outflow_splits_1,
                self.inflow_outflow_splits_0,
                self.inflow_outflow_splits,
            ]
        ):
            if splits is None:
                continue
            splits[["upstream_edge", "downstream_edge"]] = splits[
                ["upstream_edge", "downstream_edge"]
            ].replace("", None)
            for direction in ["upstream", "downstream"]:
                if f"selected_{direction}_edge" not in splits.columns:
                    splits[f"selected_{direction}_edge"] = splits[f"{direction}_edge"]
                else:
                    splits[f"selected_{direction}_edge"] = None
                inflow_outflow_splits = splits
            break
        
        self.inflow_outflow_nodes, self.inflow_outflow_edges = (
            find_nodes_edges_for_direction(
                nodes=self.nodes,
                edges=self.edges,
                node_ids=self.inflow_outflow_points["representative_node"].to_numpy(),
                border_node_ids=self.inflow_outflow_points[
                    "representative_node"
                ].to_numpy(),
                direction=self.direction,
                split_points=inflow_outflow_splits,
            )
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
        logging.info(
            "   x search for split points based on the basins of outflow/inflow points"
        )

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
        
        inflow_outflow_nodes = define_list_upstream_downstream_edges_ids(
            self.inflow_outflow_nodes.nodeID.unique(),
            self.inflow_outflow_nodes,
            self.inflow_outflow_edges,
        )
        inflow_outflow_nodes = inflow_outflow_nodes[
            inflow_outflow_nodes.no_downstream_edges > 1
        ]

        for i_node, node in enumerate(inflow_outflow_nodes.nodeID.values):
            if i_node % 50 == 0:
                logging.debug(
                    f"    - detect points: {i_node}/{len(inflow_outflow_nodes)}"
                )
            upstream_edges = self.inflow_outflow_edges[
                self.inflow_outflow_edges[node_search] == node
            ].copy()

            inflow_outflow_points_columns = [
                f"{self.direction}_node_{n}"
                for n in self.inflow_outflow_points.representative_node.values
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

        for edge in [f"{search_direction}_edge", f"{opposite_direction}_edge"]:
            self.inflow_outflow_splits_0[edge] = self.inflow_outflow_splits_0.apply(
                lambda x: None if len(x[edge + "s"].split(",")) > 1 else x[edge + "s"],
                axis=1,
            )

        if self.inflow_outflow_splits is None:
            self.inflow_outflow_splits_1 = self.inflow_outflow_splits_0.copy()
        else:
            inflow_outflow_splits = self.inflow_outflow_splits.copy()
            inflow_outflow_splits = define_list_upstream_downstream_edges_ids(
                inflow_outflow_splits.nodeID.unique(),
                inflow_outflow_splits,
                self.inflow_outflow_edges,
            )
            for edge in [f"{search_direction}_edge", f"{opposite_direction}_edge"]:
                inflow_outflow_splits[edge] = inflow_outflow_splits.apply(
                    lambda x: x[edge]
                    if len(x[edge + "s"].split(",")) > 1
                    else x[edge + "s"],
                    axis=1,
                )
            self.inflow_outflow_splits_1 = (
                pd.concat([inflow_outflow_splits, self.inflow_outflow_splits_0])
                .reset_index(drop=True)
                .drop_duplicates(subset="nodeID", keep="first")
            )

        if self.inflow_outflow_splits is None:
            logging.debug(f"    - no. of splits as input: {0}")
        else:
            logging.debug(
                f"    - no. of splits as input: {len(self.inflow_outflow_splits)}"
            )
        logging.debug(
            f"    - no. of splits found in network: {len(self.inflow_outflow_splits_0)}"
        )
        logging.debug(
            f"    - no. of splits in total: {len(self.inflow_outflow_splits_1)}"
        )

        return self.inflow_outflow_splits_1

    def export_detected_split_points(self):
        if self.inflow_outflow_splits_1 is None:
            logging.info(
                "   x splitsingen nog niet gevonden. gebruik functie .detect_split_points()"
            )
        else:
            results_dir = Path(self.path, self.dir_results)
            file_detected_points = "inflow_outflow_splits_detected.gpkg"
            logging.info(
                f"   x split points found: saved as {self.dir_results}/{file_detected_points}"
            )
            detected_inflow_outflow_splits = self.inflow_outflow_splits_1.drop(
                columns=["selected_upstream_edge", "selected_downstream_edge"],
                errors="ignore",
            )[
                [
                    "nodeID",
                    "downstream_edges",
                    "no_downstream_edges",
                    "downstream_edge",
                    "upstream_edges",
                    "no_upstream_edges",
                    "upstream_edge",
                    "geometry",
                ]
            ]
            detected_inflow_outflow_splits.to_file(Path(results_dir, file_detected_points))


    def calculate_angles_of_edges_at_nodes(self):
        self.inflow_outflow_nodes, self.inflow_outflow_edges = (
            calculate_angles_of_edges_at_nodes(
                nodes=self.inflow_outflow_nodes,
                edges=self.inflow_outflow_edges
            )
        )
        return self.inflow_outflow_nodes


    def select_directions_for_splits_based_on_angle(self):
        self.inflow_outflow_splits_1 = self.inflow_outflow_splits_0.copy()

        for index, row in self.inflow_outflow_splits_1.iterrows():
            if self.direction == "upstream":
                # Get upstream angles and edges
                upstream_angles = row["upstream_angles"]
                print(upstream_angles)
                upstream_edges = row["upstream_edges"]
                downstream_angles_str = row["downstream_angles"]
                downstream_edges_str = row["downstream_edges"]
                # Convert strings to lists
                downstream_edges_list = downstream_edges_str.split(",")
                downstream_angles_list = downstream_angles_str.split(",")

                # Assuming there's a reference angle, e.g., the first angle in the list
                reference_angle = upstream_angles

                selected_edge = find_edge_smallest_angle_difference(
                    reference_angle, downstream_angles_list, downstream_edges_list
                )

                # Update the selected columns
                self.inflow_outflow_splits_1.at[index, "selected_downstream_edge"] = (
                    selected_edge
                )

            elif self.direction == "downstream":
                # Get downstream angles and edges
                downstream_angles = row["downstream_angles"]
                downstream_edges = row["downstream_edges"]
                upstream_angles_str = row["upstream_angles"]
                upstream_edges_str = row["upstream_edges"]
                upstream_edges_list = upstream_edges_str.split(",")
                upstream_angles_list = upstream_angles_str.split(",")

                # Assuming there's a reference angle, e.g., the first angle in the list
                reference_angle = downstream_angles

                selected_edge = find_edge_smallest_angle_difference(
                    reference_angle, upstream_angles_list, upstream_edges_list
                )

                # Update the selected columns
                self.inflow_outflow_splits_1.at[index, "selected_upstream_edge"] = (
                    selected_edge
                )

        return self.inflow_outflow_splits_1

    def select_directions_for_splits(self, fillna_with_random=False):
        # check whether to use inflow_outflow_splits_1 or inflow_outflow_splits_0
        if (
            self.inflow_outflow_splits_1 is not None
            and not self.inflow_outflow_splits_1.empty
        ):
            self.inflow_outflow_splits_2 = self.inflow_outflow_splits_1.copy()
        elif (
            self.inflow_outflow_splits_0 is not None
            and not self.inflow_outflow_splits_0.empty
        ):
            self.inflow_outflow_splits_2 = self.inflow_outflow_splits_0.copy()
        else:
            logging.debug(" x no splits found: no direction for splits selected")
            return None

        logging.info("   x search for direction in splits")
        for search_direction in ["upstream", "downstream"]:
            no_splits_known = len(
                self.inflow_outflow_splits_2[
                    ~self.inflow_outflow_splits_2[f"{search_direction}_edge"].isna()
                ]
            )
            logging.debug(
                f"    - known {search_direction} direction at splits: {no_splits_known}/{len(self.inflow_outflow_splits_2)}"
            )
            self.inflow_outflow_splits_2[f"selected_{search_direction}_edge"] = (
                self.inflow_outflow_splits_2.apply(
                    lambda x: random.choice(x[f"{search_direction}_edges"].split(","))
                    if x[f"{search_direction}_edge"] is None
                    else x[f"{search_direction}_edge"],
                    axis=1,
                )
            )
            logging.debug(
                f"     - randomly choosen {search_direction} direction at splits: {
                    len(self.inflow_outflow_splits_2) - no_splits_known}/{len(self.inflow_outflow_splits_2)
                }"
            )
        return self.inflow_outflow_splits_2

    def assign_drainage_units_to_outflow_points_based_on_id(self):
        self.inflow_outflow_edges["code"] = self.inflow_outflow_edges["code"].astype(
            str
        )

        upstream_downstream_columns = [
            f"{self.direction}_node_{node}"
            for node in self.inflow_outflow_points["representative_node"]
            .unique()
            .tolist()
        ]
        self.inflow_outflow_areas = self.afwateringseenheden.merge(
            self.inflow_outflow_edges[
                ["code"] + [f"{column}" for column in upstream_downstream_columns]
            ],
            how="left",
            on="code",
        )
        self.inflow_outflow_areas[upstream_downstream_columns] = (
            self.inflow_outflow_areas[upstream_downstream_columns].fillna(False)
        )


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
        result = result.drop(columns=["code", "savedgeom"]).reset_index(drop=True)

        upstream_columns = [
            f"{self.direction}_node_{node}"
            for node in self.inflow_outflow_points["representative_node"].tolist()
        ]
        self.inflow_outflow_areas = result.merge(
            self.inflow_outflow_edges[
                ["code"] + [f"{column}" for column in upstream_columns]
            ],
            how="left",
            left_on="code_hydroobject",
            right_on="code",
        ).reset_index(drop=True)
        self.inflow_outflow_areas = self.inflow_outflow_areas.loc[
            :, ~self.inflow_outflow_areas.columns.duplicated()
        ].copy()

        for col in upstream_columns:
            self.inflow_outflow_areas[col] = self.inflow_outflow_areas[col].fillna(
                False
            )
        return self.inflow_outflow_areas


    def dissolve_assigned_drainage_units(self):
        if self.inflow_outflow_areas is None:
            return None
        inflow_outflow_areas = None
        for node in list(self.inflow_outflow_points["representative_node"].unique()):
            filtered_areas = self.inflow_outflow_areas[
                self.inflow_outflow_areas[f"{self.direction}_node_{node}"]
            ].copy()
            # Step 2: Dissolve the filtered geometries
            dissolved_areas = filtered_areas[["geometry"]].dissolve().explode()
            dissolved_areas["inflow_outflow_point"] = node
            dissolved_areas["area"] = dissolved_areas.geometry.area

            if inflow_outflow_areas is None:
                inflow_outflow_areas = dissolved_areas.reset_index(drop=True)
            else:
                inflow_outflow_areas = pd.concat(
                    [inflow_outflow_areas, dissolved_areas]
                ).reset_index(drop=True)

        self.inflow_outflow_areas = inflow_outflow_areas.copy()
        self.inflow_outflow_areas.geometry = remove_holes_from_polygons(
            inflow_outflow_areas.geometry, min_area=50
        )
        self.inflow_outflow_areas.geometry = self.inflow_outflow_areas.geometry.buffer(
            0.1
        ).buffer(-0.1)
        return self.inflow_outflow_areas


    def export_results_all(
        self,
        html_file_name: str = None,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
    ):
        """Export results to geopackages and folium html"""
        self.export_results_to_gpkg()
        self.generate_folium_map(
            html_file_name=html_file_name,
            width_edges=width_edges,
            opacity_edges=opacity_edges,
        )

    def export_results_to_gpkg(self):
        """Export results to geopackages in folder 1_tussenresultaat"""
        results_dir = Path(self.path, self.dir_results)
        logging.info(f"  x export results")
        for layer in [
            "inflow_outflow_points",
            "inflow_outflow_edges",
            "inflow_outflow_nodes",
            "inflow_outflow_areas",
        ]:
            result = getattr(self, layer)
            if result is None:
                logging.debug(f"    - {layer} not available")
            else:
                logging.debug(f"    - {layer} ({len(result)})")
                result.to_file(Path(results_dir, f"{layer}.gpkg"))


    def generate_folium_map(
        self,
        html_file_name: str = None,
        include_areas: bool = True,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
        open_html: bool = False,
        base_map: str = "OpenStreetMap",
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
        logging.info(f"   x saving html file")

        nodes_selection = self.inflow_outflow_points.representative_node.to_numpy()
        no_nodes = len(self.inflow_outflow_points) + 1
        nodes_colors = plt.get_cmap("hsv", no_nodes)
        i_nodes_colors = np.arange(start=0, stop=no_nodes - 1)
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
                self.afwateringseenheden[["geometry"]].explode(ignore_index=True),
                fill_opacity=0.0,
                color="grey",
                weight=0.5,
                z_index=10,
                name="Afwateringseenheden",
            ).add_to(m)

        inflow_outflow = "instroom" if self.direction == "downstream" else "uitstroom"

        fg = folium.FeatureGroup(
            name=f"{inflow_outflow.capitalize()}punten",
            control=True,
            show=True,
            z_index=0,
        ).add_to(m)
        for i, node_selection in enumerate(nodes_selection):
            inflow_outflow = (
                "instroom" if self.direction == "downstream" else "uitstroom"
            )
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
                z_index=2,
            ).add_to(m)

            sel_inflow_outflow_edges = self.inflow_outflow_edges[
                self.inflow_outflow_edges[f"{self.direction}_node_{node_selection}"]
            ].copy()

            sel_inflow_outflow_edges.geometry = sel_inflow_outflow_edges.buffer(
                width_edges / 2.0
            )
            if len(sel_inflow_outflow_edges) > 0:
                folium.GeoJson(
                    sel_inflow_outflow_edges[["geometry"]],
                    color=c,
                    weight=5,
                    z_index=2,
                    opacity=opacity_edges,
                    fill_opacity=opacity_edges,
                ).add_to(fg)

            if self.inflow_outflow_areas is not None and include_areas:
                inflow_outflow_areas_node = self.inflow_outflow_areas[
                    self.inflow_outflow_areas[f"inflow_outflow_point"] == node_selection
                ].copy()
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
                    self.inflow_outflow_splits[["upstream_edge", "downstream_edge"]]
                    .dropna(how="all")
                    .index
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
            search_direction = (
                "upstream" if self.direction == "downstream" else "downstream"
            )
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

        m = add_basemaps_to_folium_map(m=m, base_map=base_map)
        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m

        if html_file_name is None:
            html_file_name = self.name

        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m
