import logging
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.folium_map import generate_folium_map

from ..assets import waterschappen_order_codes
from ..utils.network_functions import (
    calculate_angles_of_edges_at_nodes,
    define_list_upstream_downstream_edges_ids,
    find_node_edge_ids_in_directed_graph,
    select_downstream_upstream_edges_angle,
)


class GeneratorGebiedsOrde(GeneratorBasis):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str | None = "1_resultaat"

    waterschap: str = None

    hydroobject: gpd.GeoDataFrame = None
    hydroobject_processed_0: gpd.GeoDataFrame = None
    hydroobject_processed_1: gpd.GeoDataFrame = None

    snapping_distance: float = 0.05

    rws_water: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    required_results: list[str] = [
        "hydroobject_processed_0", 
        "rws_water",
        "overige_watergang_processed_3", 
        "outflow_nodes_overig",
        "nodes_hydro",
        "edges_hydro",
    ]

    outflow_edges_hydro: gpd.GeoDataFrame = None
    outflow_nodes_hydro: gpd.GeoDataFrame = None

    outflow_nodes_overig: gpd.GeoDataFrame = None
    overige_watergang: gpd.GeoDataFrame = None
    overige_watergang_processed_3: gpd.GeoDataFrame = None

    edges_hydro: gpd.GeoDataFrame = None
    nodes_hydro: gpd.GeoDataFrame = None
    graph_hydro: nx.DiGraph = None

    edges_overig: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: Path = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.use_processed_hydroobject(force_preprocess=False)


    def read_outflow_nodes_hydro_with_rws_code(self, outflow_nodes_hydro=None, buffer_outflow_nodes_hydro=50.0):
        """Uses the outflow_nodes_hydro from the basisdata directory to search for closest outflow points of 
        waterways. Only waterways outflow points within a certain distance (buffer_outflow_nodes_hydro) are used.

        Args:
            outflow_nodes_hydro (gpd.GeoDataFrame): 
                outflow_nodes_hydro locations with columns: rws_code, rws_code_no, order_code, order_no
            buffer_outflow_nodes_hydro (float, optional): 
                maximum distance from node to waterway outflow point. Defaults to 50.0.

        Returns:
            gpd.GeoDataFrame: outflow_edges_hydro 
            gpd.GeoDataFrame: outflow_nodes_hydro
        """
        if outflow_nodes_hydro is not None:
            self.outflow_nodes_hydro = outflow_nodes_hydro

        logging.info("   x find final end nodes hydroobjects (dead ends)")
        # get dead end nodes
        dead_end_edges_hydro = self.edges_hydro[
            ~self.edges_hydro.node_end.isin(self.edges_hydro.node_start.values)
        ].copy()
        dead_end_edges_hydro = dead_end_edges_hydro.rename(columns={"code": "edge_code"})

        dead_end_nodes = dead_end_edges_hydro.copy()
        dead_end_nodes.geometry = dead_end_nodes.geometry.apply(lambda x: Point(x.coords[-1]))

        # For each point, find the closest endpoint
        def closest_point(base_point, points):
            distances = points.geometry.distance(base_point)
            points.loc[:, "distance"] = distances
            return points.loc[distances.idxmin()]

        logging.info("   x get dead ends (nodes/edges_hydro) closest to outflow_nodes_hydro")
        dead_end_nodes["distance"] = 0.0
        self.outflow_nodes_hydro[["edge_code", "node_end", "distance", "geometry"]] = self.outflow_nodes_hydro.geometry.apply(
            lambda x: closest_point(x, dead_end_nodes[["edge_code", "node_end", "distance", "geometry"]])
        )
        if self.outflow_nodes_hydro["distance"].max() > buffer_outflow_nodes_hydro:
            no_outflow_nodes_hydro = self.outflow_nodes_hydro.loc[
                self.outflow_nodes_hydro["distance"]>buffer_outflow_nodes_hydro, 
                "order_code"
            ].values
            logging.info(f"   x no dead ends close to outflow_nodes_hydro {no_outflow_nodes_hydro}")
            self.outflow_nodes_hydro = self.outflow_nodes_hydro[self.outflow_nodes_hydro["distance"] <= buffer_outflow_nodes_hydro]

        self.outflow_edges_hydro = self.edges_hydro[["node_end", "geometry"]].merge(
            self.outflow_nodes_hydro.drop(columns=["geometry"]),
        )
        return self.outflow_edges_hydro, self.outflow_nodes_hydro


    def generate_rws_code_for_outflow_points(self, search_range_outflow_nodes_hydro=50.0):
        """Generates an RWS code for al outflow points into rws water bodies. '
        These are the points where the water flows out of the management area of the water board and therefore the start of the orde codes of the edges_hydro.

        Parameters
        ----------
        search_range_outflow_nodes_hydro : float, optional
            buffers around the RWS water polygons, ensures that outflow points intersect with the RWS water, by default 10.0

        Returns
        -------
        self.outflow_edges_hydro: gpd.GeoDataFrame
            Geodataframe containing the outflow edges_hydro into the RWS waters
        """
        logging.info("   x find final end nodes hydroobjects (dead ends)")
        dead_end_edges_hydro = self.edges_hydro[
            ~self.edges_hydro.node_end.isin(self.edges_hydro.node_start.values)
        ].copy()
        dead_end_edges_hydro = dead_end_edges_hydro[["code", "node_end", "geometry"]]
        dead_end_edges_hydro = dead_end_edges_hydro.rename(columns={"code": "edge_code"})

        logging.info("   x generating order code for all outflow edges_hydro")
        rws_water_buffer = self.rws_water[["geometry", "rws_code"]].copy()
        rws_water_buffer.geometry = rws_water_buffer.buffer(search_range_outflow_nodes_hydro)

        outflow_edges_hydro = (
            dead_end_edges_hydro.sjoin(rws_water_buffer)
            .drop(columns="index_right")
            .reset_index(drop=True)
        )

        # Check waterschap name is right and then get rws_order_code range (min/max)
        def get_rws_order_code_no(waterschap):
            if waterschap is None or waterschap not in waterschappen_order_codes["Waterschap"].values:
                logging.error("  x Waterschap is not set, cannot generate outflow points. use waterschap=")
                for waterschap in waterschappen_order_codes["Waterschap"].values:
                    logging.error(f"    - {waterschap}")
                raise ValueError("Waterschap is not set or not valid.")
            
            rws_order_code_min = waterschappen_order_codes.set_index("Waterschap").loc[
                waterschap, "order_code_min"
            ]
            rws_order_code_max = waterschappen_order_codes.set_index("Waterschap").loc[
                waterschap, "order_code_max"
            ]
            return rws_order_code_min, rws_order_code_max

        rws_order_code_min, rws_order_code_max = get_rws_order_code_no(self.waterschap)
        outflow_edges_hydro_all_waters = None
        for rws_code in outflow_edges_hydro.rws_code.unique():
            outflows = outflow_edges_hydro[outflow_edges_hydro.rws_code==rws_code].reset_index(
                drop=True
            ).copy()
            outflows["rws_code_no"] = rws_order_code_min + outflows.index
            outflows["order_code"] = outflows.apply(
                lambda x: f"{x.rws_code}.{str(x.rws_code_no).zfill(3)}", axis=1
            )
            if outflows["rws_code_no"].max() > rws_order_code_max:
                logging_message = f" XXX aantal outflow_nodes_hydroen op RWS-water ({rws_code}) hoger dan range order_code waterschap"
                logging.info(logging_message)

            if outflow_edges_hydro_all_waters is None:
                outflow_edges_hydro_all_waters = outflows.copy()
            else:
                outflow_edges_hydro_all_waters = pd.concat([outflow_edges_hydro_all_waters, outflows])
            logging_message = (
                f"     - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            )
            logging.info(logging_message)

        logging.info(
            f"     - total no. outflow_point on outside waters for {self.name}: {len(outflow_edges_hydro_all_waters)}"
        )

        # set outflow_edges_hydro and set order_no to 2        
        self.outflow_edges_hydro = outflow_edges_hydro_all_waters.reset_index(drop=True)
        self.outflow_edges_hydro["order_no"] = 2
        self.outflow_edges_hydro["distance"] = 0.0

        # get outflow nodes from outflow edges_hydro
        self.outflow_nodes_hydro = self.outflow_edges_hydro.copy()
        self.outflow_nodes_hydro.geometry = self.outflow_nodes_hydro.geometry.apply(
            lambda x: Point(x.coords[-1])
        )
        return self.outflow_edges_hydro, self.outflow_nodes_hydro


    def generate_order_level_for_hydroobjects(self, max_order_no: int = 1000):
        """Generates the order level of the hydroobjects. A hydroobject will get the same orde as the downstream edge. 
        When a hydroobject is split at a node, the edge with the larger angle difference will get the order number downstream +1.
        The order level will keeping increasing until the complete network has an order level.
        """
        logging.info(f"   x generate order levels for hydroobjects")
        outflow_edges_hydro_orders = self.outflow_edges_hydro.copy()

        edges_hydro_left = self.edges_hydro.drop(
            columns=[
                "rws_code",
                "rws_code_no",
                "rws_order_code",
                "order_no",
                "outflow_edge",
                "order_edge_no",
            ],
            errors="ignore",
        ).copy()
        nodes_left = self.nodes_hydro.copy()

        order_no = 2

        edges_hydro_all_orders = None
        nodes_all_orders = None
        
        while order_no <= max_order_no:
            outflow_edges_hydro_order = outflow_edges_hydro_orders[
                outflow_edges_hydro_orders["order_no"] == order_no
            ].copy()

            if outflow_edges_hydro_order.empty:
                break
            logging.info(
                f"     - order {order_no}: {len(outflow_edges_hydro_order)} outflow edges_hydro"
            )
            outflow_edges_hydro_edges_hydro = None
            outflow_edges_hydro_nodes = None

            new_outflow_edges_hydro_order = None

            for i_edge, outflow_edge in outflow_edges_hydro_order.reset_index(drop=True).iterrows():

                # logging.info(f"      * {i_node+1}/{len(outflow_edges_hydro)}")
                outflow_edge_nodes_id, outflow_edge_edges_hydro_code, new_outflow_edges_hydro = (
                    find_node_edge_ids_in_directed_graph(
                        from_node_ids=edges_hydro_left.node_start.to_numpy(),
                        to_node_ids=edges_hydro_left.node_end.to_numpy(),
                        edge_ids=edges_hydro_left.code.to_numpy(),
                        outflow_edge_ids=[outflow_edge["edge_code"]],
                        search_node_ids=nodes_left.nodeID.to_numpy(),
                        search_edge_ids=edges_hydro_left.code.to_numpy(),
                        border_node_ids=None,
                        direction="upstream",
                        split_points=nodes_left,
                        order_first=True,
                        new_order_at_splits=True,
                    )
                )

                outflow_edge_edges_hydro = pd.DataFrame(
                    data={
                        "code": outflow_edge_edges_hydro_code[0],
                        "rws_code": outflow_edge["rws_code"],
                        "rws_code_no": outflow_edge["rws_code_no"],
                        "rws_order_code": outflow_edge["order_code"],
                        "order_no": order_no,
                        "outflow_edge": outflow_edge["edge_code"],
                        "order_edge_no": range(len(outflow_edge_edges_hydro_code[0])),
                    }
                )
                outflow_edge_edges_hydro = edges_hydro_left.merge(
                    outflow_edge_edges_hydro, 
                    how="right", 
                    on="code"
                )
                outflow_edge_edges_hydro = outflow_edge_edges_hydro.sort_values(
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
                outflow_edge_nodes = self.nodes_hydro.drop(
                    columns=[
                        "rws_code",
                        "rws_code_no",
                        "rws_order_code",
                        "order_no",
                        "outflow_edge",
                        "order_node_no",
                    ], errors="ignore"
                ).merge(
                    outflow_edge_nodes, how="right", on="nodeID"
                )
                outflow_edge_nodes = (
                    outflow_edge_nodes
                    .sort_values("order_no")
                    .drop_duplicates(subset="nodeID", keep=False)
                )

                if outflow_edges_hydro_edges_hydro is None:
                    outflow_edges_hydro_edges_hydro = outflow_edge_edges_hydro.copy()
                else:
                    outflow_edges_hydro_edges_hydro = pd.concat([
                        outflow_edges_hydro_edges_hydro, 
                        outflow_edge_edges_hydro.copy()
                    ])
                if outflow_edges_hydro_nodes is None:
                    outflow_edges_hydro_nodes = outflow_edge_nodes.copy()
                else:
                    outflow_edges_hydro_nodes = pd.concat([
                        outflow_edges_hydro_nodes, 
                        outflow_edge_nodes.copy()
                    ])

                new_outflow_edges_hydro = self.edges_hydro[
                    self.edges_hydro["code"].isin(new_outflow_edges_hydro)
                ][["code", "geometry", "node_end"]]
                new_outflow_edges_hydro = new_outflow_edges_hydro.rename(
                    columns={"code": "edge_code"}
                )
                new_outflow_edges_hydro["rws_code"] = outflow_edge["rws_code"]
                new_outflow_edges_hydro["rws_code_no"] = outflow_edge["rws_code_no"]
                new_outflow_edges_hydro["order_code"] = outflow_edge["order_code"]
                new_outflow_edges_hydro["order_no"] = order_no + 1

                if new_outflow_edges_hydro_order is None:
                    new_outflow_edges_hydro_order = new_outflow_edges_hydro.copy()
                else:
                    new_outflow_edges_hydro_order = pd.concat(
                        [
                            new_outflow_edges_hydro_order, 
                            new_outflow_edges_hydro,
                        ]
                    )
            
            # collect all edges_hydro with associated order code
            edges_hydro_order = outflow_edges_hydro_edges_hydro.copy()
            if edges_hydro_all_orders is None:
                edges_hydro_all_orders = edges_hydro_order.copy()
            else:
                edges_hydro_all_orders = pd.concat([edges_hydro_all_orders, edges_hydro_order])

            # collect all nodes with associated order code
            nodes_order = outflow_edges_hydro_nodes.copy()
            if nodes_all_orders is None:
                nodes_all_orders = nodes_order.copy()
            else:
                nodes_all_orders = pd.concat([nodes_all_orders, nodes_order])

            # find all edges_hydro not yet linked to an order
            edges_hydro_left = edges_hydro_left[
                ~edges_hydro_left["code"].isin(edges_hydro_order.code)
            ].copy()
            nodes_left = nodes_left[
                (nodes_left["nodeID"].isin(edges_hydro_left.node_start))
                | (nodes_left["nodeID"].isin(edges_hydro_left.node_end))
            ].copy()
            
            # filter 
            new_outflow_edges_hydro_order = new_outflow_edges_hydro_order[
                ~new_outflow_edges_hydro_order.edge_code.isin(edges_hydro_all_orders.code)
            ]
            new_outflow_edges_hydro_order = new_outflow_edges_hydro_order.drop_duplicates(
                subset="edge_code", 
                keep="first"
            )

            # get all outflow edges_hydro for all orders.
            if outflow_edges_hydro_orders is None:
                outflow_edges_hydro_orders = new_outflow_edges_hydro_order.copy()
            else:
                outflow_edges_hydro_orders = pd.concat(
                    [
                        outflow_edges_hydro_orders, 
                        new_outflow_edges_hydro_order,
                    ]
                )

            order_no = order_no + 1

        edges_hydro_left["rws_code"] = ""
        edges_hydro_left["rws_code_no"] = -999
        edges_hydro_left["order_no"] = -999
        edges_hydro_left["outflow_edge"] = -999
        edges_hydro_left["order_edge_no"] = -999

        nodes_left["rws_code"] = ""
        nodes_left["rws_code_no"] = -999
        nodes_left["order_no"] = -999
        nodes_left["outflow_edge"] = -999
        nodes_left["order_node_no"] = -999

        edges_hydro_all_orders = edges_hydro_all_orders.sort_values("order_no").drop_duplicates(
            subset="code", keep="first"
        )
        nodes_all_orders = nodes_all_orders.sort_values("order_no").drop_duplicates(
            subset="nodeID", keep="first"
        )
        outflow_edges_hydro_orders = outflow_edges_hydro_orders.sort_values(
            "order_no"
        ).drop_duplicates(subset=["node_end", "edge_code"], keep="first")

        self.edges_hydro = pd.concat([edges_hydro_all_orders, edges_hydro_left]).reset_index(drop=True)
        self.edges_hydro["source"] = "hydroobject"
        self.nodes_hydro = pd.concat([nodes_all_orders, nodes_left]).reset_index(drop=True)
        
        self.outflow_edges_hydro = outflow_edges_hydro_orders.copy()
        self.outflow_nodes_hydro = self.outflow_edges_hydro.copy()
        self.outflow_nodes_hydro["geometry"] = self.outflow_nodes_hydro["geometry"].apply(
            lambda x: Point(x.coords[-1])
        )
        self.outflow_nodes_hydro = (
            self.outflow_nodes_hydro.groupby("node_end")
            .agg({
                col: list if col == "edge_code" else "first"
                for col in self.outflow_nodes_hydro.columns
            })
            .rename(columns={"edge_code": "edge_codes"})
        )
        self.outflow_nodes_hydro = gpd.GeoDataFrame(
            self.outflow_nodes_hydro, geometry="geometry", crs=self.outflow_edges_hydro.crs
        ).reset_index(drop=True)

        len_edges_hydro_with_order = len(self.edges_hydro[self.edges_hydro.order_no > 0])
        len_edges_hydro_without_order = len(self.edges_hydro[self.edges_hydro.order_no < 0])
        logging_message = f"     - order levels generated: {len_edges_hydro_with_order} edges_hydro - {len_edges_hydro_without_order} left"
        logging.info(logging_message)
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "edges_hydro",
                "nodes_hydro",
                "outflow_edges_hydro",
                "outflow_nodes_hydro",
            ])


    def generate_order_code_for_hydroobjects(self, order_for_each_edge=False):
        """
        Generates a code for each hydroobject based on the order level. Hydroobjects with a higher order level will get a longer code, that shows in to which hydroobject with a lower code it flows.
        """
        logging.info(f"   x generate order code for edges_hydro")
        edges_hydro = (
            self.edges_hydro[
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
        edges_hydro_orders = None
        for order_no in [
            order_no for order_no in edges_hydro.order_no.unique() if order_no > 0
        ]:
            # get edges_hydro from this order and next order
            edges_hydro_order = edges_hydro[edges_hydro.order_no == order_no].copy()
            edges_hydro_next_order = edges_hydro[edges_hydro.order_no == order_no + 1].copy()
            edges_hydro_next_order = edges_hydro_next_order[["code", "node_end"]].rename(
                columns={"code": "edge_codes", "node_end": "node_end2"}
            )

            # find where edges_hydro from next order are entering this order
            edges_hydro_order = edges_hydro_order.merge(
                edges_hydro_next_order,
                how="left",
                left_on="node_start",
                right_on="node_end2",
            ).drop(columns="node_end2")

            # sort and get inflow edges_hydro from next order (aggregate to list)
            sort_columns = [
                "rws_code",
                "rws_code_no",
                "order_no",
                "outflow_edge",
                "order_edge_no",
            ]
            edges_hydro_order = (
                edges_hydro_order.sort_values(sort_columns)
                .groupby(["node_start", "node_end"])
                .agg(
                    {
                        k: list if k == "edge_codes" else "first"
                        for k in edges_hydro_order.columns
                    }
                )
                .reset_index(drop=True)
                .sort_values(sort_columns)
            )

            if edges_hydro_orders is None:
                edges_hydro_orders = edges_hydro_order.copy()
            else:
                edges_hydro_orders = pd.concat([edges_hydro_orders, edges_hydro_order]).copy()
        
        # remove nan from list in case of no inflow_edge from next order
        edges_hydro_orders["edge_codes"] = edges_hydro_orders["edge_codes"].apply(
            lambda x: [s for s in x if isinstance(s, str)]
        )
        # get number of inflow_edges_hydro from next order
        edges_hydro_orders["no_edge_codes"] = edges_hydro_orders["edge_codes"].apply(
            lambda x: len(x)
        )
        edges_hydro_orders = gpd.GeoDataFrame(
            edges_hydro_orders, geometry="geometry", crs=28992
        ).reset_index(drop=True)

        # Generate order codes
        edges_hydro_outflow_edges_hydro = None
        edges_hydro_outflow_edges_hydro_all = None

        logging.info(f"     - generate order code: per order")

        for order_no in sorted([
            order_no for order_no in edges_hydro.order_no.unique() if order_no > 0
        ]):
            order_outflow_edges_hydro = self.outflow_edges_hydro[
                self.outflow_edges_hydro.order_no == order_no
            ]
            logging_message = (
                f"     - order {order_no}: {len(order_outflow_edges_hydro)} outflow edges_hydro"
            )
            logging.info(logging_message)

            if edges_hydro_outflow_edges_hydro is not None:
                # get from existing edges_hydro (lower order) the new outflow edges_hydro
                new_outflow_edge = edges_hydro_outflow_edges_hydro.loc[
                    (edges_hydro_outflow_edges_hydro["no_edge_codes"] > 0) & 
                    (edges_hydro_outflow_edges_hydro["order_no"] > 0)
                ].copy()
                new_outflow_edge = new_outflow_edge[[
                    "code", "edge_codes", "order_no", "no_edge_codes", "order_code_no", "order_code"
                ]]
                new_outflow_edge = new_outflow_edge.rename(
                    columns={"edge_codes": "edge_code"}
                )
                new_outflow_edge = new_outflow_edge.explode(column="edge_code")
                
                # check for new outflow_edges_hydro for system:
                old_order_outflow_edges_hydro = order_outflow_edges_hydro[~order_outflow_edges_hydro["distance"].isna()].copy()
                old_order_outflow_edges_hydro["order_code_no"] = 1
                
                # outflow_edges_hydro to lower order:
                new_order_outflow_edges_hydro = order_outflow_edges_hydro[order_outflow_edges_hydro["distance"].isna()].copy()
                new_order_outflow_edges_hydro = new_order_outflow_edges_hydro.drop(
                    columns="order_code"
                ).merge(
                    new_outflow_edge[["edge_code", "order_code", "order_code_no"]],
                    how="left",
                    on="edge_code",
                )

                loc_order_code_nan = new_order_outflow_edges_hydro["order_code"].isna()
                new_order_outflow_edges_hydro.loc[loc_order_code_nan, "order_code_no"] = -999
                new_order_outflow_edges_hydro.loc[~loc_order_code_nan, "order_code_no"] = (
                    new_order_outflow_edges_hydro.loc[~loc_order_code_nan, "order_code_no"] + 1
                )
                new_order_outflow_edges_hydro.loc[~loc_order_code_nan, "order_code"] = (
                    new_order_outflow_edges_hydro.loc[~loc_order_code_nan, "order_code"].str[:-4]
                )

                edges_hydro_duplicates = new_order_outflow_edges_hydro[
                    ["node_end", "order_code_no"]
                ].duplicated()

                while edges_hydro_duplicates.sum() > 0:
                    new_order_outflow_edges_hydro["order_code_no"] = new_order_outflow_edges_hydro[
                        "order_code_no"
                    ] + edges_hydro_duplicates.astype(int)
                    edges_hydro_duplicates = new_order_outflow_edges_hydro[
                        ["node_end", "order_code_no"]
                    ].duplicated()

                new_order_outflow_edges_hydro["order_code"] = (
                    new_order_outflow_edges_hydro["order_code"]
                    + "."
                    + new_order_outflow_edges_hydro["order_code_no"].astype(str).str.zfill(3)
                )

                order_outflow_edges_hydro = pd.concat([
                    old_order_outflow_edges_hydro,
                    new_order_outflow_edges_hydro
                ])

            # Get edges_hydro from order x
            edges_hydro_outflow_edges_hydro = edges_hydro_orders[
                edges_hydro_orders["order_no"] == order_no
            ].copy()
            edges_hydro_outflow_edges_hydro["order_code_no"] = edges_hydro_outflow_edges_hydro["no_edge_codes"]
            
            if order_for_each_edge:
                # in case of new ordercode for each edge give all a number 1 (+ no_inflow_edges_hydro)
                edges_hydro_outflow_edges_hydro["order_code_no"] = (
                    edges_hydro_outflow_edges_hydro["order_code_no"] + 1
                )
            else:
                # in case of not a new ordercode for each edge 
                # give all a number 0 and a 1 at location of inflow_edges_hydro (+number of inflow_edges_hydro)
                edges_hydro_outflow_edges_hydro.loc[
                    edges_hydro_outflow_edges_hydro["order_code_no"] > 0, "order_code_no"
                ] = edges_hydro_outflow_edges_hydro.loc[
                    edges_hydro_outflow_edges_hydro["order_code_no"] > 0, "order_code_no"
                ] + 1
            
            # give all outflow_edges_hydro should a number 1
            # edges_hydro_outflow_edges_hydro.loc[
            #     edges_hydro_outflow_edges_hydro["code"].isin(order_outflow_edges_hydro["edge_code"]),
            #     "order_code_no",
            # ] = 1
            edges_hydro_outflow_edges_hydro["order_code"] = None

            # Loop through all outflow_edges_hydro
            for i_outflow_edge, outflow_edge in order_outflow_edges_hydro.iterrows():
                # get all edges_hydro with this outflow_edge
                loc_edges_hydro_outflow_edge = (
                    edges_hydro_outflow_edges_hydro.outflow_edge == outflow_edge.edge_code
                )
                if order_for_each_edge:
                    # shift 1, fillna with 1 and cumsum to get order_code_no
                    edges_hydro_outflow_edges_hydro.loc[loc_edges_hydro_outflow_edge, "order_code_no"] = (
                        edges_hydro_outflow_edges_hydro.loc[loc_edges_hydro_outflow_edge, "order_code_no"]
                        .shift(1)
                        .fillna(1)
                        .cumsum()
                    )
                else:
                    # shift 1, fillna with 0 and cumsum to get order_code_no
                    edges_hydro_outflow_edges_hydro.loc[loc_edges_hydro_outflow_edge, "order_code_no"] = (
                        edges_hydro_outflow_edges_hydro.loc[loc_edges_hydro_outflow_edge, "order_code_no"]
                        .shift(1)
                        .fillna(1)
                        .cumsum()
                    )
                    # edges_hydro_outflow_edges_hydro.loc[
                    #     edges_hydro_outflow_edges_hydro.code == edges_hydro_outflow_edges_hydro.outflow_edge,
                    #     "order_code_no",
                    # ] = 1
                
                # give all edges_hydro the order code of the outflow_edge
                edges_hydro_outflow_edges_hydro.loc[loc_edges_hydro_outflow_edge, "order_code"] = (
                    outflow_edge["order_code"]
                )

            # add order_code_no to order_code - str with zfill(3)
            edges_hydro_outflow_edges_hydro["order_code"] = (
                edges_hydro_outflow_edges_hydro["order_code"]
                + "."
                + edges_hydro_outflow_edges_hydro["order_code_no"].astype(str).str.zfill(3)
            )
            if edges_hydro_outflow_edges_hydro_all is None:
                edges_hydro_outflow_edges_hydro_all = edges_hydro_outflow_edges_hydro.copy()
            else:
                edges_hydro_outflow_edges_hydro_all = pd.concat(
                    [edges_hydro_outflow_edges_hydro_all, edges_hydro_outflow_edges_hydro]
                ).copy()

        edges_hydro_outflow_edges_hydro_all = (
            gpd.GeoDataFrame(
                edges_hydro_outflow_edges_hydro_all, geometry="geometry", crs=self.edges_hydro.crs
            )
            .reset_index(drop=True)
            .drop(columns=["edge_codes"])
        )

        self.edges_hydro = self.edges_hydro.drop(
            columns=["no_edge_codes", "order_code_no", "order_code"], errors="ignore"
        ).merge(
            edges_hydro_outflow_edges_hydro_all[
                ["code", "no_edge_codes", "order_code_no", "order_code"]
            ],
            how="left",
            on="code",
        )
        self.edges_hydro["no_edge_codes"] = self.edges_hydro["no_edge_codes"].fillna(-999)
        self.edges_hydro["order_code_no"] = self.edges_hydro["order_code_no"].fillna(-999)
        self.edges_hydro["order_code"] = self.edges_hydro["order_code"].fillna("")

        for direction in ["upstream", "downstream"]:
            self.nodes_hydro[f"{direction}_order_no"] = self.nodes_hydro.apply(
                lambda x: ",".join(
                    self.edges_hydro.loc[
                        self.edges_hydro["code"].isin(x[f"{direction}_edges"].split(",")), "order_no"
                    ].astype(str).values
                ),
                axis=1,
            )
            self.nodes_hydro[f"{direction}_order_code"] = self.nodes_hydro.apply(
                lambda x: ",".join(
                    self.edges_hydro.loc[
                        self.edges_hydro["code"].isin(x[f"{direction}_edges"].split(",")), "order_code"
                    ].astype(str).values
                ),
                axis=1,
            )
        
        self.hydroobject_processed_1 = self.edges_hydro.copy()
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "hydroobject_processed_1",
                "edges_hydro",
                "nodes_hydro",
            ])
        return self.hydroobject_processed_1


    def generate_order_no_order_code_for_other_waterlines(self):
        """Generates order level for overige watergangen, based on their outflow point in the hydroobjects. 

        Returns
        -------
        self.outflow_nodes_overig: gpd.GeoDataFrame
            Geodataframe containing the outflow points of the overige watergang in the hydroobjects
        self.edges_overig: gpd.GeoDataFrame
            Geodataframe containing the processed overige watergangen, including the order levels and order codes
        """
        logging.info(f"   x generate order code for overige watergangen")

        def string_to_list(string, sep=",", str_type=int):
            if str_type==int:
                return [str_type(s) if s!="" else -1 for s in string.split(sep)]
            else:
                return [str_type(s) if s!="" else "" for s in string.split(sep)]

        def list_to_string(lst, sep=","):
            return ",".join([str(i) for i in lst])

        # check if values are strings and change into lists
        if self.outflow_nodes_overig is None:
            return None, None
        
        logging.info(f"     - coupling overige watergangen to outflow_nodes_hydro")
        
        outflow_nodes_overig = self.outflow_nodes_overig[
            ["nodeID", "geometry"]
        ].sjoin(
            self.nodes_hydro[
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
        outflow_nodes_hydro = outflow_nodes_overig[
            [
                "nodeID",
                "downstream_edges",
                "selected_downstream_edge",
                "downstream_order_no",
                "downstream_order_code",
            ]
        ].copy()

        for col in ["downstream_edges", "downstream_order_no", "downstream_order_code"]:
            outflow_nodes_hydro[col] = outflow_nodes_hydro[col].fillna("")
            outflow_nodes_hydro[col] = outflow_nodes_hydro[col].apply(
                lambda x: string_to_list(
                    x, sep=",", str_type=int if col == "downstream_order_no" else str
                )
            )


        outflow_nodes_hydro["no_downstream_edges"] = outflow_nodes_hydro["downstream_edges"].apply(lambda x: len(x))
        outflow_nodes_hydro["no_downstream_order_no"] = outflow_nodes_hydro["downstream_order_no"].apply(lambda x: len(x))
        outflow_nodes_hydro["no_downstream_order_code"] = outflow_nodes_hydro["downstream_order_code"].apply(lambda x: len(x))
        outflow_nodes_hydro["equal"] = outflow_nodes_hydro.apply(
            lambda x: (x.no_downstream_edges == x.no_downstream_order_no) and (x.no_downstream_edges == x.no_downstream_order_code),
            axis=1
        )

        outflow_nodes_hydro = (
            outflow_nodes_hydro.explode(
                ["downstream_edges", "downstream_order_no", "downstream_order_code"]
            )
            .sort_values(["nodeID", "downstream_order_no", "downstream_order_code"])
            .drop_duplicates()
        )
        outflow_nodes_hydro = outflow_nodes_hydro.groupby("nodeID").first().reset_index()

        # filter out outflow_nodes_hydro without downstream_order_no
        logging.info(f"     - filter overige watergangen without downstream order no")
        outflow_nodes_overig = outflow_nodes_overig[
            ~outflow_nodes_overig["downstream_order_no"].isnull()
        ]

        self.outflow_nodes_overig = (
            outflow_nodes_overig.drop(
                columns=["downstream_edges", "downstream_order_no", "downstream_order_code"],
                errors="ignore",
            ).merge(
                outflow_nodes_hydro[
                    ["nodeID", "downstream_edges", "downstream_order_no", "downstream_order_code"]
                ],
                how="left",
                on="nodeID",
            )
        ).drop(columns="index_right")

        if self.overige_watergang_processed_3 is not None:
            edges_overig = self.overige_watergang_processed_3.merge(
                self.outflow_nodes_overig[
                    ["nodeID", "downstream_edges", "downstream_order_no", "downstream_order_code"]
                ],
                how="left",
                left_on="outflow_nodes_hydro",
                right_on="nodeID",
            )
            edges_overig = edges_overig.sort_values("downstream_order_code").drop_duplicates(subset="geometry", keep="first")

            # filter out waterways with outflow_nodes_hydro without downstream_order_no
            edges_overig = edges_overig[edges_overig["downstream_order_no"] > 0]

            # use order_no of downstream principal waterways
            edges_overig["order_no"] = edges_overig["downstream_order_no"].astype(int) + 1
            edges_overig["order_code_no"] = (
                edges_overig.groupby("downstream_order_code").cumcount().fillna(-1000).astype(int)
                + 1
            )
            edges_overig["order_code"] = (
                edges_overig["downstream_order_code"]
                + "-X"
                + edges_overig["order_code_no"].astype(str).str.zfill(4)
            )
            self.edges_overig = edges_overig.copy()
            self.edges_overig["source"] = "overige_watergang"
            logging.info(f"     - overige watergangen with order code: {len(self.edges_overig)}")
        
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "edges_hydro",
                "nodes_hydro",
                "outflow_edges_hydro",
                "outflow_nodes_hydro",
                "outflow_nodes_overig",
                "edges_overig",
            ])
        return (
            self.outflow_nodes_overig,
            self.edges_overig,
        )


    def generate_folium_map(self, html_file_name:str="", **kwargs):
        if html_file_name == "":
            html_file_name = self.name + "_order_code"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            **kwargs
        )
        return self.folium_map
