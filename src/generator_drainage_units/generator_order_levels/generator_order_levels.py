import logging
from pathlib import Path
import webbrowser
import datetime

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import (
    define_list_upstream_downstream_edges_ids,
    calculate_angles_of_edges_at_nodes,
    select_downstream_upstream_edges,
    find_node_edge_ids_in_directed_graph,
)
from ..utils.folium_utils import (
    add_labels_to_points_lines_polygons,
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
)


class GeneratorOrderLevels(GeneratorBasis):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    path: Path = None
    name: str = None
    base_dir: Path = None
    waterschap: str = None
    range_order_code_min: int = None
    range_order_code_max: int = None

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed: gpd.GeoDataFrame = None

    # overige_watergangen: gpd.GeoDataFrame = None
    # overige_watergangen_processed: gpd.GeoDataFrame = None

    rws_wateren: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    dead_end_edges: gpd.GeoDataFrame = None
    dead_start_edges: gpd.GeoDataFrame = None
    outflow_edges: gpd.GeoDataFrame = None

    dead_end_nodes: gpd.GeoDataFrame = None
    dead_start_nodes: gpd.GeoDataFrame = None

    outflow_nodes: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None
    network_positions: dict = None

    folium_map: folium.Map = None
    folium_html_path: Path = None

    def create_graph_from_network(
        self, water_lines=["rivieren", "hydroobjecten", "hydroobjecten_extra"]
    ):
        if water_lines is None:
            water_lines = ["hydroobjecten"]
        edges = None
        for water_line in water_lines:
            gdf_water_line = getattr(self, water_line)
            if gdf_water_line is None:
                continue
            if edges is None:
                edges = gdf_water_line.explode()
            else:
                edges = pd.concat([edges, gdf_water_line.explode()])
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        self.network_positions = {n: [n[0], n[1]] for n in list(self.graph.nodes)}
        logging.info(
            f"   x create network graph ({len(self.edges)} edges, {len(self.nodes)} nodes)"
        )
        return self.nodes, self.edges, self.graph

    def define_list_upstream_downstream_edges_ids(self):
        self.nodes = define_list_upstream_downstream_edges_ids(
            node_ids=self.nodes.nodeID.values, nodes=self.nodes, edges=self.edges
        )

    def calculate_angles_of_edges_at_nodes(self):
        logging.info("   x calculate angles of edges to nodes")
        self.nodes, self.edges = calculate_angles_of_edges_at_nodes(
            nodes=self.nodes, edges=self.edges
        )
        return self.nodes

    def select_downstream_upstream_edges(self, min_difference_angle: str = 20.0):
        logging.info("   x find downstream upstream edges")
        self.nodes = select_downstream_upstream_edges(self.nodes, min_difference_angle)
        return self.nodes

    def find_end_points_hydroobjects(self, buffer_width=0.5, direction="upstream"):
        # Copy hydroobject data to new variable 'hydroobjects' and make dataframes with start and end nodes
        logging.info("   x find start and end nodes hydrobojects")

        dead_end_edges = self.edges[
            ~self.edges.node_end.isin(self.edges.node_start.values)
        ].copy()
        dead_end_edges["nodeID"] = dead_end_edges["node_end"]
        dead_end_edges = dead_end_edges[["code", "nodeID", "geometry"]].rename(
            columns={"code": "edge_code"}
        )
        dead_end_nodes = dead_end_edges.copy()
        dead_end_nodes["geometry"] = dead_end_nodes["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )

        dead_start_edges = self.edges[
            ~self.edges.node_start.isin(self.edges.node_end.values)
        ].copy()
        dead_start_edges["nodeID"] = dead_start_edges["node_start"]
        dead_start_edges = dead_start_edges[["code", "nodeID", "geometry"]].rename(
            columns={"code": "edge_code"}
        )
        dead_start_nodes = dead_start_edges.copy()
        dead_start_nodes["geometry"] = dead_start_nodes["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )

        self.dead_end_edges = dead_end_edges.copy()

        self.dead_end_nodes = dead_end_nodes.copy()
        self.dead_start_edges = dead_start_edges.copy()
        self.dead_start_nodes = dead_start_nodes.copy()
        return self.dead_end_nodes, self.dead_start_nodes

    def generate_rws_code_for_all_outflow_points(self, buffer_rws=10.0):
        logging.info("   x generating order code for all outflow points")
        rws_wateren = self.rws_wateren.copy()
        rws_wateren.geometry = rws_wateren.geometry.buffer(buffer_rws)

        outflow_edges = (
            self.dead_end_edges.sjoin(rws_wateren[["geometry", "rws_code"]])
            .drop(columns="index_right")
            .reset_index(drop=True)
        )

        outflow_edges_all_waters = None
        for rws_code in outflow_edges.rws_code.unique():
            outflows = outflow_edges[outflow_edges.rws_code == rws_code].reset_index(
                drop=True
            )
            outflows["rws_code_no"] = self.range_order_code_min + outflows.index
            outflows["order_code"] = outflows.apply(
                lambda x: f"{x.rws_code}.{str(x.rws_code_no)}", axis=1
            )
            if outflows["rws_code_no"].max() > self.range_order_code_max:
                logging_message = f" XXX aantal uitstroompunten op RWS-water ({rws_code}) hoger dan range order_code waterschap"
                logging.debug(logging_message)

            if outflow_edges_all_waters is None:
                outflow_edges_all_waters = outflows.copy()
            else:
                outflow_edges_all_waters = pd.concat(
                    [outflow_edges_all_waters, outflows]
                )
            logging_message = (
                f"    - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            )
            logging.debug(logging_message)

        logging_message = f"    - total no. outflow_point on outside waters for {self.name}: {len(outflow_edges_all_waters)}"
        logging.debug(logging_message)

        self.outflow_edges = outflow_edges_all_waters.rename(
            columns={"nodeID": "node_end"}
        ).reset_index(drop=True)
        self.outflow_edges["order_no"] = 2
        self.outflow_nodes = self.outflow_edges.copy()
        self.outflow_nodes["geometry"] = self.outflow_nodes["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )

        return self.outflow_edges

    def generate_orde_level_for_hydroobjects(self):
        logging.info(f"   x generate order levels for hydroobjects")
        # self.outflow_edges = self.outflow_edges[self.outflow_edges["order_no"]==2].copy()

        edges_left = self.edges.copy()  # .drop(columns=["order_edge_no"]).copy()
        nodes_left = self.nodes.copy()  # .drop(columns=["order_node_no"]).copy()

        order_no = 2

        edges_all_orders = None
        nodes_all_orders = None
        outflow_edges_orders = self.outflow_edges[
            self.outflow_edges["order_no"] == 2
        ].copy()
        new_outflow_edges_order = outflow_edges_orders.copy()
        new_outflow_edges_order = new_outflow_edges_order.rename(
            columns={"nodeID": "node_end"}
        )

        while not new_outflow_edges_order.empty and order_no < 100:
            logging_message = (
                f"    - order {order_no}: {len(new_outflow_edges_order)} outflow edges"
            )
            logging.debug(logging_message)
            outflow_edges = new_outflow_edges_order.copy()
            new_outflow_edges_order = None
            outflow_edges_edges = None
            outflow_edges_nodes = None

            for i_edge, outflow_edge in outflow_edges.reset_index(drop=True).iterrows():
                # logging.debug(f"      * {i_node+1}/{len(outflow_edges)}")
                outflow_edge_nodes_id, outflow_edge_edges_code, new_outflow_edges = (
                    find_node_edge_ids_in_directed_graph(
                        from_node_ids=edges_left.node_start.to_numpy(),
                        to_node_ids=edges_left.node_end.to_numpy(),
                        edge_ids=edges_left.code.to_numpy(),
                        outflow_edge_ids=[outflow_edge["edge_code"]],
                        search_node_ids=nodes_left.nodeID.to_numpy(),
                        search_edge_ids=edges_left.code.to_numpy(),
                        border_node_ids=None,
                        direction="upstream",
                        split_points=nodes_left,
                        order_first=True,
                        new_order_at_splits=True,
                    )
                )

                outflow_edge_edges = pd.DataFrame(
                    data={
                        "code": outflow_edge_edges_code[0],
                        "rws_code": outflow_edge["rws_code"],
                        "rws_code_no": outflow_edge["rws_code_no"],
                        "rws_order_code": outflow_edge["order_code"],
                        "order_no": order_no,
                        "outflow_edge": outflow_edge["edge_code"],
                        "order_edge_no": range(len(outflow_edge_edges_code[0])),
                    }
                )
                outflow_edge_edges = edges_left.merge(
                    outflow_edge_edges, how="right", on="code"
                )
                outflow_edge_edges = outflow_edge_edges.sort_values(
                    "order_no"
                ).drop_duplicates(subset="code", keep=False)

                outflow_edge_nodes = pd.DataFrame(
                    data={
                        "nodeID": outflow_edge_nodes_id[0],
                        "rws_code": outflow_edge["rws_code"],
                        "rws_code_no": outflow_edge["rws_code_no"],
                        "rws_order_code": outflow_edge["order_code"],
                        "order_no": order_no,
                        "outflow_edge": outflow_edge["edge_code"],
                        "order_node_no": range(len(outflow_edge_nodes_id[0])),
                    }
                )
                outflow_edge_nodes = self.nodes.merge(
                    outflow_edge_nodes, how="right", on="nodeID"
                )
                outflow_edge_nodes = outflow_edge_nodes.sort_values(
                    "order_no"
                ).drop_duplicates(subset="nodeID", keep=False)

                if outflow_edges_edges is None:
                    outflow_edges_edges = outflow_edge_edges.copy()
                else:
                    outflow_edges_edges = pd.concat(
                        [outflow_edges_edges, outflow_edge_edges.copy()]
                    )
                if outflow_edges_nodes is None:
                    outflow_edges_nodes = outflow_edge_nodes.copy()
                else:
                    outflow_edges_nodes = pd.concat(
                        [outflow_edges_nodes, outflow_edge_nodes.copy()]
                    )

                new_outflow_edges = self.edges[
                    self.edges["code"].isin(new_outflow_edges)
                ][["code", "geometry", "node_end"]]
                new_outflow_edges = new_outflow_edges.rename(
                    columns={"code": "edge_code"}
                )
                new_outflow_edges["rws_code"] = outflow_edge["rws_code"]
                new_outflow_edges["rws_code_no"] = outflow_edge["rws_code_no"]
                new_outflow_edges["order_code"] = outflow_edge["order_code"]
                new_outflow_edges["order_no"] = order_no + 1

                if new_outflow_edges_order is not None:
                    new_outflow_edges = new_outflow_edges[
                        ~new_outflow_edges.edge_code.isin(
                            list(new_outflow_edges_order.edge_code.unique())
                        )
                    ]
                if outflow_edges_orders is not None:
                    new_outflow_edges = new_outflow_edges[
                        ~new_outflow_edges.edge_code.isin(
                            list(outflow_edges_orders.edge_code.unique())
                        )
                    ]

                if new_outflow_edges_order is None:
                    new_outflow_edges_order = new_outflow_edges.copy()
                else:
                    new_outflow_edges_order = pd.concat(
                        [new_outflow_edges_order, new_outflow_edges]
                    )

            edges_order = outflow_edges_edges.copy()
            edges_left = edges_left[~edges_left["code"].isin(edges_order.code)].copy()

            nodes_order = outflow_edges_nodes.copy()
            nodes_left = nodes_left[
                (nodes_left["nodeID"].isin(edges_left.node_start))
                | (nodes_left["nodeID"].isin(edges_left.node_end))
            ].copy()

            if edges_all_orders is None:
                edges_all_orders = edges_order.copy()
            else:
                edges_all_orders = pd.concat([edges_all_orders, edges_order])
            if nodes_all_orders is None:
                nodes_all_orders = nodes_order.copy()
            else:
                nodes_all_orders = pd.concat([nodes_all_orders, nodes_order])

            if outflow_edges_orders is None:
                outflow_edges_orders = new_outflow_edges_order.copy()
            else:
                outflow_edges_orders = pd.concat(
                    [outflow_edges_orders, new_outflow_edges_order]
                )

            order_no = order_no + 1

        edges_left["rws_code"] = ""
        edges_left["rws_code_no"] = -999
        edges_left["order_no"] = -999
        edges_left["outflow_edge"] = -999
        edges_left["order_edge_no"] = -999

        nodes_left["rws_code"] = ""
        nodes_left["rws_code_no"] = -999
        nodes_left["order_no"] = -999
        nodes_left["outflow_edge"] = -999
        nodes_left["order_node_no"] = -999

        edges_all_orders = edges_all_orders.sort_values("order_no").drop_duplicates(
            subset="code", keep="first"
        )
        nodes_all_orders = nodes_all_orders.sort_values("order_no").drop_duplicates(
            subset="nodeID", keep="first"
        )
        outflow_edges_orders = outflow_edges_orders.sort_values(
            "order_no"
        ).drop_duplicates(subset=["node_end", "edge_code"], keep="first")

        self.edges = pd.concat([edges_all_orders, edges_left]).reset_index(drop=True)
        self.nodes = pd.concat([nodes_all_orders, nodes_left]).reset_index(drop=True)
        self.outflow_edges = outflow_edges_orders.copy()
        self.outflow_nodes = self.outflow_edges.copy()
        self.outflow_nodes["geometry"] = self.outflow_nodes["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )
        self.outflow_nodes = (
            self.outflow_nodes.groupby("node_end")
            .agg(
                {
                    col: list if col == "edge_code" else "first"
                    for col in self.outflow_nodes.columns
                }
            )
            .rename(columns={"edge_code": "edge_codes"})
        )
        self.outflow_nodes = gpd.GeoDataFrame(
            self.outflow_nodes, geometry="geometry", crs=self.outflow_edges.crs
        ).reset_index(drop=True)

        len_edges_with_order = len(self.edges[self.edges.order_no > 0])
        len_edges_without_order = len(self.edges[self.edges.order_no < 0])
        logging_message = f"     - order levels generated: {len_edges_with_order} edges - {len_edges_without_order} left"
        logging.info(logging_message)

    def generate_order_code_for_edges(self, order_for_each_edge=False):
        logging.info(f"   x generate order code for edges")
        edges = (
            self.edges[
                [
                    "code",
                    "node_start",
                    "node_end",
                    "rws_code",
                    "rws_code_no",
                    "order_no",
                    "outflow_edge",
                    "order_edge_no",
                    "geometry",
                ]
            ]
            .copy()
            .sort_values(
                ["rws_code", "rws_code_no", "order_no", "outflow_edge", "order_edge_no"]
            )
        )

        # edges = edges[edges["rws_code_no"]==728]

        logging.info(f"     - generate order code: preparation")
        edges_orders = None
        for order_no in [
            order_no for order_no in edges.order_no.unique() if order_no > 0
        ]:
            edges_order = edges[edges.order_no == order_no].copy()
            edges_next_order = edges[edges.order_no == order_no + 1].copy()
            edges_order = edges_order.merge(
                edges_next_order[["code", "node_end"]].rename(
                    columns={"code": "edge_codes"}
                ),
                how="left",
                left_on="node_start",
                right_on="node_end",
            )
            sort_columns = [
                "rws_code",
                "rws_code_no",
                "order_no",
                "outflow_edge",
                "order_edge_no",
            ]
            edges_order = (
                edges_order.sort_values(sort_columns)
                .groupby("node_start")
                .agg(
                    {
                        k: list if k == "edge_codes" else "first"
                        for k in edges_order.columns
                    }
                )
                .reset_index(drop=True)
                .sort_values(sort_columns)
            )

            if edges_orders is None:
                edges_orders = edges_order.copy()
            else:
                edges_orders = pd.concat([edges_orders, edges_order]).copy()

        edges_orders["edge_codes"] = edges_orders["edge_codes"].apply(
            lambda x: [s for s in x if isinstance(s, str)]
        )
        edges_orders["no_edge_codes"] = edges_orders["edge_codes"].apply(
            lambda x: len(x)
        )

        edges_orders = gpd.GeoDataFrame(
            edges_orders, geometry="geometry", crs=28992
        ).reset_index(drop=True)

        edges_outflow_edges = None
        edges_outflow_edges_all = None

        logging.info(f"     - generate order code: per order")
        for order_no in [
            order_no for order_no in edges.order_no.unique() if order_no > 0
        ]:
            order_outflow_edges = self.outflow_edges[
                (self.outflow_edges.order_no == order_no)
            ]
            logging_message = (
                f"    - order {order_no}: {len(order_outflow_edges)} outflow edges"
            )
            logging.debug(logging_message)

            if edges_outflow_edges is not None:
                new_outflow_edge = edges_outflow_edges[
                    [
                        "code",
                        "edge_codes",
                        "order_no",
                        "no_edge_codes",
                        "order_code_no",
                        "order_code",
                    ]
                ].copy()
                new_outflow_edge = new_outflow_edge.rename(
                    columns={"edge_codes": "edge_code"}
                )
                new_outflow_edge = new_outflow_edge.explode(column="edge_code")
                new_outflow_edge = new_outflow_edge.loc[
                    (new_outflow_edge["no_edge_codes"] > 0)
                    & (new_outflow_edge["order_no"] > 0)
                ]
                order_outflow_edges = order_outflow_edges.drop(
                    columns="order_code"
                ).merge(
                    new_outflow_edge[["edge_code", "order_code"]],
                    how="left",
                    on="edge_code",
                )
                loc_order_code_nan = order_outflow_edges["order_code"].isna()

                order_outflow_edges["order_code_no"] = -999
                order_outflow_edges.loc[~loc_order_code_nan, "order_code_no"] = (
                    order_outflow_edges.loc[~loc_order_code_nan, "order_code"]
                    .str[-3:]
                    .astype(int)
                )
                order_outflow_edges.loc[~loc_order_code_nan, "order_code_no"] = (
                    order_outflow_edges.loc[~loc_order_code_nan, "order_code_no"] + 1
                )
                order_outflow_edges.loc[~loc_order_code_nan, "order_code"] = (
                    order_outflow_edges.loc[~loc_order_code_nan, "order_code"].str[:-4]
                )
                edges_duplicates = order_outflow_edges[
                    ["node_end", "order_code_no"]
                ].duplicated()

                while edges_duplicates.sum() > 0:
                    order_outflow_edges["order_code_no"] = order_outflow_edges[
                        "order_code_no"
                    ] + edges_duplicates.astype(int)
                    edges_duplicates = order_outflow_edges[
                        ["node_end", "order_code_no"]
                    ].duplicated()

                order_outflow_edges["order_code"] = (
                    order_outflow_edges["order_code"]
                    + "."
                    + order_outflow_edges["order_code_no"].astype(str).str.zfill(3)
                )

            edges_outflow_edges = edges_orders[
                edges_orders["order_no"] == order_no
            ].copy()

            edges_outflow_edges["order_code_no"] = 0

            edges_outflow_edges["order_code_no"] = edges_outflow_edges["no_edge_codes"]
            if order_for_each_edge:
                edges_outflow_edges["order_code_no"] = (
                    edges_outflow_edges["order_code_no"] + 1
                )
            else:
                edges_outflow_edges.loc[
                    edges_outflow_edges["order_code_no"] > 0, "order_code_no"
                ] = (
                    edges_outflow_edges.loc[
                        edges_outflow_edges["order_code_no"] > 0, "order_code_no"
                    ]
                    + 1
                )
            edges_outflow_edges.loc[
                edges_outflow_edges["code"].isin(order_outflow_edges["edge_code"]),
                "order_code_no",
            ] = 1
            edges_outflow_edges["order_code"] = None

            for i_outflow_edge, outflow_edge in order_outflow_edges.iterrows():
                # print(f"{i_outflow_edge} / {len(order_outflow_edges)}", end="\r")
                loc_edges_outflow_edge = (
                    edges_outflow_edges.outflow_edge == outflow_edge.edge_code
                )
                if order_for_each_edge:
                    edges_outflow_edges.loc[loc_edges_outflow_edge, "order_code_no"] = (
                        edges_outflow_edges.loc[loc_edges_outflow_edge, "order_code_no"]
                        .shift(1)
                        .fillna(1)
                        .cumsum()
                    )
                else:
                    edges_outflow_edges.loc[loc_edges_outflow_edge, "order_code_no"] = (
                        edges_outflow_edges.loc[loc_edges_outflow_edge, "order_code_no"]
                        .shift(1)
                        .fillna(0)
                        .cumsum()
                    )
                    edges_outflow_edges.loc[
                        edges_outflow_edges.code==edges_outflow_edges.outflow_edge, 
                        "order_code_no"
                    ] = 1
                edges_outflow_edges.loc[loc_edges_outflow_edge, "order_code"] = (
                    outflow_edge["order_code"]
                )
            
            edges_outflow_edges["order_code"] = (
                edges_outflow_edges["order_code"]
                + "."
                + edges_outflow_edges["order_code_no"].astype(str).str.zfill(3)
            )

            if edges_outflow_edges_all is None:
                edges_outflow_edges_all = edges_outflow_edges.copy()
            else:
                edges_outflow_edges_all = pd.concat(
                    [edges_outflow_edges_all, edges_outflow_edges]
                ).copy()

        edges_outflow_edges_all = (
            gpd.GeoDataFrame(
                edges_outflow_edges_all, geometry="geometry", crs=self.edges.crs
            )
            .reset_index(drop=True)
            .drop(columns=["edge_codes"])
        )

        self.edges = self.edges.merge(
            edges_outflow_edges_all[
                ["code", "no_edge_codes", "order_code_no", "order_code"]
            ],
            how="left",
            on="code",
        )
        self.edges["no_edge_codes"] = self.edges["no_edge_codes"].fillna(-999)
        self.edges["order_code_no"] = self.edges["order_code_no"].fillna(-999)
        self.edges["order_code"] = self.edges["order_code"].fillna("")
        return self.edges

    def export_results_to_gpkg(self):
        """Export results to geopackages in folder 1_tussenresultaat"""
        results_dir = Path(self.path, self.dir_inter_results)
        logging.info(f"   x export results")
        for layer in [
            "dead_end_edges",
            "dead_start_edges",
            "outflow_edges",
            "outflow_nodes",
            "edges",
            "nodes",
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
        order_labels: bool = False,
    ):
        # Make figure
        outflow_nodes_4326 = self.outflow_nodes.to_crs(4326)

        m = folium.Map(
            location=[
                outflow_nodes_4326.geometry.y.mean(),
                outflow_nodes_4326.geometry.x.mean(),
            ],
            zoom_start=12,
            tiles=None,
        )

        folium.GeoJson(
            self.rws_wateren.geometry,
            name="RWS_Wateren",
            z_index=0,
        ).add_to(m)

        if "order_no" in self.edges.columns:
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.edges[self.edges["order_no"] > 1][
                    ["code", "order_no", "geometry"]
                ],
                layer_name="Orde-nummer watergangen",
                control=True,
                lines=True,
                line_color_column="order_no",
                line_color_cmap="hsv",
                label=False,
                line_weight=5,
                z_index=1,
            )

            if order_labels and "order_no" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"Orde-nummer watergangen (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["order_no"] > 1][
                        ["code", "order_no", "geometry"]
                    ],
                    column="order_no",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )

            if order_labels and "order_code" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"Orde-code watergangen (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                self.edges["order_code"].fill = ""

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["order_no"] > 1][
                        ["code", "order_code", "geometry"]
                    ],
                    column="order_code",
                    label_fontsize=8,
                    fg=fg,
                )

        folium.GeoJson(
            self.hydroobjecten.geometry,
            name="Watergangen",
            color="blue",
            fill_color="blue",
            zoom_on_click=True,
            show=False,
            z_index=2,
        ).add_to(m)

        fg = folium.FeatureGroup(
            name=f"Uitstroompunten RWS-water", control=True
        ).add_to(m)

        folium.GeoJson(
            self.outflow_nodes[self.outflow_nodes["order_no"] == 2],
            name="Uitstroompunten RWS-wateren",
            marker=folium.Circle(
                radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=3,
        ).add_to(fg)

        add_labels_to_points_lines_polygons(
            gdf=self.outflow_nodes[self.outflow_nodes["order_no"] == 2],
            column="order_code",
            label_fontsize=8,
            fg=fg,
        )

        folium.GeoJson(
            self.outflow_nodes[self.outflow_nodes["order_no"] > 2],
            name="Overige uitstroompunten",
            marker=folium.Circle(
                radius=15,
                fill_color="orange",
                fill_opacity=0.8,
                color="orange",
                weight=3,
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=2,
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
