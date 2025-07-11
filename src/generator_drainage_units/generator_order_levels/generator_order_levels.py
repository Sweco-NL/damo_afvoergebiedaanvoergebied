import logging
import webbrowser
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.folium_utils import (
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
    add_labels_to_points_lines_polygons,
)
from ..utils.network_functions import (
    calculate_angles_of_edges_at_nodes,
    define_list_upstream_downstream_edges_ids,
    find_node_edge_ids_in_directed_graph,
    select_downstream_upstream_edges,
)


class GeneratorOrderLevels(GeneratorBasis):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str | None = "1_resultaat"

    waterschap: str = None
    range_order_code_min: int = None
    range_order_code_max: int = None

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    hydroobjecten_processed_1: gpd.GeoDataFrame = None

    rws_wateren: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    required_results: list[str] = [
        "hydroobjecten_processed_0", 
        "rws_wateren",
        "overige_watergangen_processed_3", 
        "outflow_nodes_overige_watergangen",
    ]

    outflow_edges: gpd.GeoDataFrame = None
    outflow_nodes: gpd.GeoDataFrame = None

    outflow_nodes_overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed_3: gpd.GeoDataFrame = None
    overige_watergangen_processed_4: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: Path = None

    def create_graph_from_network(self, water_lines=["hydroobjecten"], processed="processed"):
        """Turns a linestring layer containing waterlines into a graph of edges and nodes. 

        Parameters
        ----------
        water_lines : list, optional
            List of waterline files names used to create graph, must refer to geopackages containing linestrings, by default ["hydroobjecten"]

        Returns
        -------
        self.nodes: gpd.GeoDataFrame
            Geodataframe containing nodes between waterlines
        self.edges: gpd.GeoDataFrame
            Geodataframe containing edges (waterlines)
        self.graph: nx.DiGraph
            Networkx graph containing the edges and nodes
        """
        
        edges = None
        for water_line in water_lines:
            gdf_water_line = getattr(self, water_line)
            for i in range(10):
                if not hasattr(self, f"{water_line}_{processed}_{i}"):
                    break
                gdf_water_line_processed = getattr(self, f"{water_line}_{processed}_{i}")
                if gdf_water_line_processed is None:
                    break
                else:
                    gdf_water_line = gdf_water_line_processed.copy()
            if gdf_water_line is None:
                continue
            if edges is None:
                edges = gdf_water_line.explode()
            else:
                edges = pd.concat([edges, gdf_water_line.explode()])
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        logging.info(
            f"   x create network graph ({len(self.edges)} edges, {len(self.nodes)} nodes)"
        )
        return self.nodes, self.edges, self.graph

    def define_list_upstream_downstream_edges_ids(self):
        """Get the upstream and downstream edges for each node. 

        Returns
        -------
        self.nodes: gpd.GeoDataFrame
            Geodataframe containing nodes between waterlines, including upstream and downstream edges
        """
        self.nodes = define_list_upstream_downstream_edges_ids(
            node_ids=self.nodes.nodeID.values, nodes=self.nodes, edges=self.edges
        )
        return self.nodes

    def calculate_angles_of_edges_at_nodes(self):
        """Calculates the angles of the upstream and downstream edges for each node. 

        Returns
        -------
        self.nodes: gpd.GeoDataFrame
            Geodataframe containing nodes between waterlines, including upstream and downstream edges and their angles
        """
        logging.info("   x calculate angles of edges to nodes")
        self.nodes, self.edges = calculate_angles_of_edges_at_nodes(
            nodes=self.nodes, edges=self.edges
        )
        return self.nodes

    def select_downstream_upstream_edges(self, min_difference_angle: str = 20.0):
        """select the upstream or downstream edge that represents the main channel, based on the smallest angle. When the angle of both edges is too large, no edge is selected.

        Parameters
        ----------
        min_difference_angle : str, optional
            minimum , by default 20.0

        Returns
        -------
        self.nodes: gpd.GeoDataFrame
            Geodataframe containing nodes between waterlines, including the selected upstream and downstream angles
        """
        logging.info("   x find downstream upstream edges")
        self.nodes = select_downstream_upstream_edges(
            self.nodes, min_difference_angle=min_difference_angle
        )
        return self.nodes

    def generate_rws_code_for_all_outflow_points(self, buffer_rws_water=10.0):
        """Generates an RWS code for al outflow points into rws water bodies. These are the points where the water flows out of the management area of the water board and therefore the start of the orde codes of the edges.

        Parameters
        ----------
        buffer_rws_water : float, optional
            buffers around the RWS water polygons, ensures that outflow points intersect with the RWS water, by default 10.0

        Returns
        -------
        self.outflow_edges: gpd.GeoDataFrame
            Geodataframe containing the outflow edges into the RWS waters
        """
        # Copy hydroobject data to new variable 'hydroobjects' and make dataframes with start and end nodes
        logging.info("   x find start and end nodes hydroobjects")

        dead_end_edges = self.edges[
            ~self.edges.node_end.isin(self.edges.node_start.values)
        ].copy()
        dead_end_edges["nodeID"] = dead_end_edges["node_end"]
        dead_end_edges = dead_end_edges[["code", "nodeID", "geometry"]].rename(
            columns={"code": "edge_code"}
        )

        logging.info("   x generating order code for all outflow points")
        rws_wateren = self.rws_wateren.copy()
        rws_wateren.geometry = rws_wateren.geometry.buffer(buffer_rws_water)

        outflow_edges = (
            dead_end_edges.sjoin(rws_wateren[["geometry", "rws_code"]])
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
                logging.info(logging_message)

            if outflow_edges_all_waters is None:
                outflow_edges_all_waters = outflows.copy()
            else:
                outflow_edges_all_waters = pd.concat(
                    [outflow_edges_all_waters, outflows]
                )
            logging_message = (
                f"     - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            )
            logging.info(logging_message)

        logging_message = f"     - total no. outflow_point on outside waters for {self.name}: {len(outflow_edges_all_waters)}"
        logging.info(logging_message)

        self.outflow_edges = outflow_edges_all_waters.rename(
            columns={"nodeID": "node_end"}
        ).reset_index(drop=True)
        self.outflow_edges["order_no"] = 2
        self.outflow_nodes = self.outflow_edges.copy()
        self.outflow_nodes["geometry"] = self.outflow_nodes["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )
        return self.outflow_edges


    def generate_order_level_for_hydroobjects(self):
        """Generates the order level of the hydroobjects. A hydroobject will get the same orde as the downstream edge. 
        When a hydroobject is split at a node, the edge with the larger angle difference will get the order number downstream +1.
        The order level will keeping increasing until the complete network has an order level.
        """
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
                f"     - order {order_no}: {len(new_outflow_edges_order)} outflow edges"
            )
            logging.info(logging_message)
            outflow_edges = new_outflow_edges_order.copy()
            new_outflow_edges_order = None
            outflow_edges_edges = None
            outflow_edges_nodes = None

            for i_edge, outflow_edge in outflow_edges.reset_index(drop=True).iterrows():
                # logging.info(f"      * {i_node+1}/{len(outflow_edges)}")
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
                outflow_edge_edges = edges_left.drop(
                    columns=[
                        "rws_code",
                        "rws_code_no",
                        "rws_order_code",
                        "order_no",
                        "outflow_edge",
                        "order_edge_no",
                    ],
                    errors="ignore",
                ).merge(outflow_edge_edges, how="right", on="code")
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
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "outflow_edges",
                "outflow_nodes",
            ])


    def generate_order_code_for_hydroobjects(self, order_for_each_edge=False):
        """
        Generates a code for each hydroobject based on the order level. Hydroobjects with a higher order level will get a longer code, that shows in to which hydroobject with a lower code it flows.
        """
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

        logging.info(f"     - generate order code: preparation")
        edges_orders = None
        for order_no in [
            order_no for order_no in edges.order_no.unique() if order_no > 0
        ]:
            # get edges from this order
            edges_order = edges[edges.order_no == order_no].copy()
            # get edges from next order
            edges_next_order = edges[edges.order_no == order_no + 1].copy()
            # find where edges from next order are entering this order
            edges_order = edges_order.merge(
                edges_next_order[["code", "node_end"]].rename(
                    columns={"code": "edge_codes", "node_end": "node_end2"}
                ),
                how="left",
                left_on="node_start",
                right_on="node_end2",
            ).drop(columns="node_end2")
            # sort and get inflow edges from next order
            sort_columns = [
                "rws_code",
                "rws_code_no",
                "order_no",
                "outflow_edge",
                "order_edge_no",
            ]
            edges_order = (
                edges_order.sort_values(sort_columns)
                .groupby(["node_start", "node_end"])
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
                self.outflow_edges.order_no == order_no
            ]
            logging_message = (
                f"     - order {order_no}: {len(order_outflow_edges)} outflow edges"
            )
            logging.info(logging_message)

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
                        edges_outflow_edges.code == edges_outflow_edges.outflow_edge,
                        "order_code_no",
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

        self.edges = self.edges.drop(
            columns=["no_edge_codes", "order_code_no", "order_code"], errors="ignore"
        ).merge(
            edges_outflow_edges_all[
                ["code", "no_edge_codes", "order_code_no", "order_code"]
            ],
            how="left",
            on="code",
        )
        self.edges["no_edge_codes"] = self.edges["no_edge_codes"].fillna(-999)
        self.edges["order_code_no"] = self.edges["order_code_no"].fillna(-999)
        self.edges["order_code"] = self.edges["order_code"].fillna("")

        for direction in ["upstream", "downstream"]:
            self.nodes[f"{direction}_order_no"] = self.nodes.apply(
                lambda x: list(
                    self.edges.loc[
                        self.edges["code"].isin(x[f"{direction}_edges"]), "order_no"
                    ].astype(str).values
                ),
                axis=1,
            )
            self.nodes[f"{direction}_order_code"] = self.nodes.apply(
                lambda x: list(
                    self.edges.loc[
                        self.edges["code"].isin(x[f"{direction}_edges"]), "order_code"
                    ].values
                ),
                axis=1,
            )
        
        self.hydroobjecten_processed_1 = self.edges.copy()
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "hydroobjecten_processed_1",
                "edges",
                "nodes",
            ])
        return self.hydroobjecten_processed_1


    def generate_order_no_order_code_for_other_waterlines(self):
        """Generates order level for overige watergangen, based on their outflow point in the hydroobjects. 

        Returns
        -------
        self.outflow_nodes_overige_watergangen: gpd.GeoDataFrame
            Geodataframe containing the outflow points of the overige watergang in the hydroobjects
        self.overige_watergangen_processed_4: gpd.GeoDataFrame
            Geodataframe containing the processed overige watergangen, including the order levels and order codes
        """
        logging.info(f"   x generate order code for overige watergangen")

        def string_to_list(string, sep=",", type=int):
            return [type(i.strip()[1:-1]) if i!="" else -1 for i in string[1:-1].split(sep)]

        # check if values are strings and change into lists
        for direction in ["upstream", "downstream"]:
            column = f"{direction}_order_no"
            if column in self.nodes.columns and isinstance(self.nodes[column].to_numpy()[0], str):
                self.nodes[f"{direction}_order_no"] = self.nodes[f"{direction}_order_no"].apply(lambda x: string_to_list(x, sep=",", type=int))
        if isinstance(self.nodes[f"downstream_edges"].values[0], str):
            self.nodes[f"downstream_edges"] = self.nodes[f"downstream_edges"].apply(lambda x: string_to_list(x, sep=",", type=str))
        if isinstance(self.nodes[f"downstream_order_code"].values[0], str):
            self.nodes[f"downstream_order_code"] = self.nodes[f"downstream_order_code"].apply(lambda x: string_to_list(x, sep=",", type=str))
        
        if self.outflow_nodes_overige_watergangen is None:
            return None, None
        
        outflow_nodes_overige_watergangen = self.outflow_nodes_overige_watergangen[
            ["nodeID", "geometry"]
        ].sjoin(
            self.nodes[
                [
                    "downstream_edges",
                    "selected_downstream_edge",
                    "downstream_order_no",
                    "downstream_order_code",
                    "geometry",
                ]
            ],
            how="left",
        )
        outflow_nodes = outflow_nodes_overige_watergangen[
            [
                "nodeID",
                "downstream_edges",
                "selected_downstream_edge",
                "downstream_order_no",
                "downstream_order_code",
            ]
        ].copy()

        outflow_nodes = (
            outflow_nodes.explode(
                ["downstream_edges", "downstream_order_no", "downstream_order_code"]
            )
            .sort_values(["nodeID", "downstream_order_no", "downstream_order_code"])
            .drop_duplicates()
        )
        outflow_nodes = outflow_nodes.groupby("nodeID").first().reset_index()

        self.outflow_nodes_overige_watergangen = (
            self.outflow_nodes_overige_watergangen.drop(
                columns=["downstream_edges", "downstream_order_no", "downstream_order_code"],
                errors="ignore",
            ).merge(
                outflow_nodes[
                    ["nodeID", "downstream_edges", "downstream_order_no", "downstream_order_code"]
                ],
                how="left",
                on="nodeID",
            )
        )

        edges = self.overige_watergangen_processed_3.merge(
            self.outflow_nodes_overige_watergangen[
                ["nodeID", "downstream_edges", "downstream_order_no", "downstream_order_code"]
            ],
            how="left",
            left_on="outflow_node",
            right_on="nodeID",
        )
        edges = edges.sort_values("downstream_order_code").drop_duplicates(subset="geometry", keep="first")
        
        edges["order_no"] = edges["downstream_order_no"].astype(int) + 1
        edges["order_code_no"] = (
            edges.groupby("downstream_order_code").cumcount().fillna(-1000).astype(int)
            + 1
        )
        edges["order_code"] = (
            edges["downstream_order_code"]
            + "-X"
            + edges["order_code_no"].astype(str).str.zfill(4)
        )
        self.overige_watergangen_processed_4 = edges.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "outflow_edges",
                "outflow_nodes",
                "outflow_nodes_overige_watergangen",
                "overige_watergangen_processed_4",
            ])
        return (
            self.outflow_nodes_overige_watergangen,
            self.overige_watergangen_processed_4,
        )


    def generate_folium_map(
        self,
        html_file_name: str = None,
        include_areas: bool = True,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
        open_html: bool = False,
        base_map: str = "Light Mode",
        order_labels: bool = True,
    ):
        """Generate folium maps to display the results of the order level generator

        Parameters
        ----------
        html_file_name : str, optional
            _description_, by default None
        include_areas : bool, optional
            _description_, by default True
        width_edges : float, optional
            _description_, by default 10.0
        opacity_edges : float, optional
            _description_, by default 0.5
        open_html : bool, optional
            _description_, by default False
        base_map : str, optional
            _description_, by default "Light Mode"
        order_labels : bool, optional
            _description_, by default True

        Returns
        -------
        _type_
            _description_
        """
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

        folium.GeoJson(
            self.hydroobjecten.geometry,
            name="AB-Watergangen",
            color="black",
            fill_color="black",
            zoom_on_click=True,
            show=True,
            z_index=2,
        ).add_to(m)

        if "order_no" in self.edges.columns:
            edges = self.edges[self.edges["order_no"] > 1][
                ["code", "order_no", "order_code", "geometry"]
            ].sort_values("order_no", ascending=False)

            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges,
                layer_name="AB-watergangen - Orde-nummer",
                control=True,
                lines=True,
                line_color_column="order_no",
                line_color_cmap="hsv_r",
                label=False,
                line_weight=5,
                z_index=1,
            )

            edges_labels = edges.drop_duplicates(subset="order_code", keep="first")

            if order_labels and "order_no" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"AB-watergangen - Orde-nummer (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=edges_labels,
                    column="order_no",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )

            if order_labels and "order_code" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"AB-watergangen - Orde-code (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=edges_labels,
                    column="order_code",
                    label_fontsize=8,
                    center=True,
                    fg=fg,
                )

        fg = folium.FeatureGroup(
            name=f"Uitstroompunten in RWS-water", control=True
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

        if self.outflow_nodes_overige_watergangen is not None:
            outflow_nodes_overige_watergangen = (
                self.outflow_nodes_overige_watergangen.sjoin(
                    self.nodes.drop(columns=["x", "y", "nodeID"]), how="left"
                )
            )
            folium.GeoJson(
                outflow_nodes_overige_watergangen.geometry,
                name="C-watergangen - Uitstroompunten",
                marker=folium.Circle(
                    radius=3,
                    fill_color="orange",
                    fill_opacity=0.8,
                    color="black",
                    weight=3,
                ),
                show=False,
                highlight_function=lambda x: {"fillOpacity": 0.8},
                zoom_on_click=True,
                z_index=2,
            ).add_to(m)

        if (
            self.overige_watergangen_processed_4 is not None
            and "outflow_node" in self.overige_watergangen_processed_4.columns
        ):
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.overige_watergangen_processed_4[
                    ["outflow_node", "geometry"]
                ],
                layer_name=f"C-Watergangen - Gegroepeerd per uitstroompunt",
                control=True,
                lines=True,
                line_color_column="outflow_node",
                line_color_cmap=None,
                show=False,
                z_index=1,
            )

            if order_labels and "order_code" in self.overige_watergangen_processed_4:
                fg = folium.FeatureGroup(
                    name=f"C-Watergangen - Orde-codering",
                    control=True,
                    show=False,
                    z_index=2,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.overige_watergangen_processed_4,
                    column="order_code",
                    label_fontsize=7,
                    label_decimals=0,
                    center=True,
                    fg=fg,
                )

        m = add_basemaps_to_folium_map(m=m, base_map=base_map)

        folium.LayerControl(collapsed=False).add_to(m)
        m.add_child(folium.plugins.MeasureControl())

        self.folium_map = m

        if html_file_name is None:
            html_file_name = self.name + "_order_code"

        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m
