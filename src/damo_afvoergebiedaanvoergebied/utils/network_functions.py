import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import logging
from ..utils.general_functions import calculate_angle_difference, calculate_angle


def find_predecessors_graph(
    from_node_ids: np.array,
    to_node_ids: np.array,
    edge_ids: np.array,
    edge_id: int,
    border_node_ids: np.array = None,
    pred_nodes: list = [],
    pred_edges: list = [],
):
    """
    Find predecessors within graph for specified node_id.

    Note: recursive function!
    """
    node_id = to_node_ids[np.where(edge_ids == edge_id)]
    from_node = from_node_ids[np.where(edge_ids == edge_id)]
    pred_node = from_node_ids[np.where(to_node_ids == from_node)]
    pred_edge = edge_ids[np.where(to_node_ids == from_node)]

    for i in range(pred_node.shape[0]):
        p = pred_node[i]
        e = pred_edge[i]
        if e not in pred_edges:
            pred_edges = pred_edges + [e]
        if p not in pred_nodes:
            pred_nodes = pred_nodes + [int(p)]
            if border_node_ids is None or p not in border_node_ids:
                pred_nodes, pred_edges = find_predecessors_graph(
                    from_node_ids,
                    to_node_ids,
                    edge_ids,
                    e,
                    border_node_ids,
                    pred_nodes,
                    pred_edges,
                )
    return pred_nodes, pred_edges


def find_predecessors_graph_with_splits(
    from_node_ids,
    to_node_ids,
    edge_ids,
    edge_id,
    border_node_ids=None,
    split_node_edge_ids=None,
    split_node_edge_ids2=None,
    pred_nodes=[],
    pred_edges=[],
    new_outflow_edges=[],
):
    """
    Find predecessors within graph for specified node_id.

    Note: recursive function!
    """
    node_id = to_node_ids[np.where(edge_ids == edge_id)][0]
    from_node = from_node_ids[np.where(edge_ids == edge_id)][0]
    pred_node = from_node_ids[np.where(to_node_ids == from_node)]
    pred_edge = edge_ids[np.where(to_node_ids == from_node)]

    for i in range(pred_node.shape[0]):
        p = pred_node[i]
        e = pred_edge[i]

        # if edge_id in ["WL_88", "WL_86", "WL_284"]: #"WL_269", "WL_106618-0", "WL_271-8"]:
        # logging.debug("xxxxxxxxxxxx")
        # logging.debug(edge_id)
        # logging.debug(from_node)
        # logging.debug(pred_node)
        # logging.debug(pred_edge)
        # logging.debug(split_node_edge_ids2)
        # logging.debug("xxxxxxxxxxxx")

        if (
            split_node_edge_ids2 is not None
            and from_node in split_node_edge_ids2
            and split_node_edge_ids2[from_node] != e
        ):
            # logging.debug(e)
            new_outflow_edges = new_outflow_edges + [e]
            continue

        if e in pred_edges:
            continue

        pred_edges = pred_edges + [e]
        pred_nodes = pred_nodes + [int(p)]

        if border_node_ids is not None and p in border_node_ids:
            continue

        if (
            split_node_edge_ids is None
            or p not in split_node_edge_ids
            or split_node_edge_ids[p] == e
        ):
            pred_nodes, pred_edges, new_outflow_edges = (
                find_predecessors_graph_with_splits(
                    from_node_ids,
                    to_node_ids,
                    edge_ids,
                    e,
                    border_node_ids,
                    split_node_edge_ids,
                    split_node_edge_ids2,
                    pred_nodes,
                    pred_edges,
                    new_outflow_edges,
                )
            )
    return pred_nodes, pred_edges, new_outflow_edges


def accumulate_values_graph(
    from_node_ids,
    to_node_ids,
    node_ids,
    edge_ids,
    values_node_ids,
    values_nodes,
    values_edge_ids,
    values_edges,
    border_node_ids=None,
    itself=False,
    direction="upstream",
    decimals=None,
):
    """Calculate for all node_ids the accumulated values of all predecessors with values."""
    len_node_ids = np.shape(node_ids)[0]
    results_nodes = np.zeros(np.shape(node_ids))
    results_edges = np.zeros(np.shape(node_ids))
    logging.info(f"accumulate values using graph for {len(node_ids)} node(s)")
    for i in range(node_ids.shape[0]):
        print(f" * {i+1}/{len_node_ids} ({(i+1)/len(node_ids):.2%})", end="\r")
        node_id = node_ids[i]
        if direction == "upstream":
            pred_nodes, pred_edges = find_predecessors_graph(
                from_node_ids,
                to_node_ids,
                edge_ids,
                node_id,
                border_node_ids,
                np.array([]),
            )
        else:
            pred_nodes, pred_edges = find_predecessors_graph(
                to_node_ids,
                from_node_ids,
                edge_ids,
                node_id,
                border_node_ids,
                np.array([]),
            )
        if itself:
            pred_nodes = np.append(pred_nodes, node_id)
        pred_nodes_sum = np.sum(
            values_nodes[np.searchsorted(values_node_ids, pred_nodes)]
        )
        pred_edges_sum = np.sum(
            values_edges[np.searchsorted(values_edge_ids, pred_edges)]
        )
        if decimals is None:
            results_nodes[i] = pred_nodes_sum
            results_edges[i] = pred_edges_sum
        else:
            results_nodes[i] = np.round(pred_nodes_sum, decimals=decimals)
            results_edges[i] = np.round(pred_edges_sum, decimals=decimals)
    return results_nodes, results_edges


def find_node_edge_ids_in_directed_graph(
    from_node_ids,
    to_node_ids,
    edge_ids,
    outflow_edge_ids,
    search_node_ids,
    search_edge_ids,
    border_node_ids=None,
    direction="upstream",
    split_points=None,
    set_logging=False,
    order_first=False,
    new_order_at_splits=False,
):
    results_nodes = [
        [
            int(to_node_ids[np.where(edge_ids == e)][0]),
            int(from_node_ids[np.where(edge_ids == e)][0]),
        ]
        for e in outflow_edge_ids
    ]
    results_edges = [[e] if e in search_edge_ids else [] for e in outflow_edge_ids]
    if set_logging:
        logging.debug(
            f"     - find {direction} nodes/edges for {len(outflow_edge_ids)}/{len(search_edge_ids)} nodes"
        )
    search_direction = "upstream" if direction == "downstream" else "downstream"
    opposite_direction = "downstream" if direction == "downstream" else "upstream"
    if isinstance(outflow_edge_ids, list):
        outflow_edge_ids = np.array(outflow_edge_ids)

    new_outflow_edges = []
    if split_points is None:
        for i in range(outflow_edge_ids.shape[0]):
            # print(f" * {i+1}/{len_outflow_edge_ids} ({(i+1)/len(outflow_edge_ids):.2%})")
            edge_id = outflow_edge_ids[i]
            if direction == "upstream":
                pred_nodes, pred_edges = find_predecessors_graph(
                    from_node_ids, to_node_ids, edge_ids, edge_id, border_node_ids, []
                )
            else:
                pred_nodes, pred_edges = find_predecessors_graph(
                    to_node_ids, from_node_ids, edge_ids, edge_id, border_node_ids, []
                )
            results_nodes[i] = results_nodes[i] + [
                p for p in pred_nodes if p in search_node_ids
            ]
            results_edges[i] = results_edges[i] + [
                p for p in pred_edges if p in search_edge_ids
            ]
    else:
        split_node_edge_ids = split_points.set_index("nodeID")[
            f"selected_{search_direction}_edge"
        ].to_dict()
        split_node_edge_ids2 = split_points.set_index("nodeID")[
            f"selected_{opposite_direction}_edge"
        ].to_dict()

        if not new_order_at_splits:
            split_node_edge_ids = {
                k: v for k, v in split_node_edge_ids.items() if v not in [None, ""]
            }
            split_node_edge_ids2 = {
                k: v for k, v in split_node_edge_ids2.items() if v not in [None, ""]
            }

        for i in range(outflow_edge_ids.shape[0]):
            # print(f" * {i+1}/{len_outflow_nodes_ids} ({(i+1)/len_outflow_nodes_ids:.2%})", end="\r")
            edge_id = outflow_edge_ids[i]
            if direction == "upstream":
                pred_nodes, pred_edges, new_outflow_edges = (
                    find_predecessors_graph_with_splits(
                        from_node_ids,
                        to_node_ids,
                        edge_ids,
                        edge_id,
                        border_node_ids,
                        split_node_edge_ids if not order_first else None,
                        split_node_edge_ids2,
                        [],
                        [],
                    )
                )
            else:
                pred_nodes, pred_edges, new_outflow_edges = (
                    find_predecessors_graph_with_splits(
                        to_node_ids,
                        from_node_ids,
                        edge_ids,
                        edge_id,
                        border_node_ids,
                        split_node_edge_ids2,
                        split_node_edge_ids if not order_first else None,
                        [],
                        [],
                    )
                )
            results_nodes[i] = results_nodes[i] + [
                p for p in pred_nodes if p in search_node_ids
            ]
            results_edges[i] = results_edges[i] + [
                p for p in pred_edges if p in search_edge_ids
            ]
    return results_nodes, results_edges, new_outflow_edges


def find_nodes_edges_for_direction(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    outflow_edge_ids: list,
    border_node_ids: list = None,
    direction: str = "upstream",
    split_points: gpd.GeoDataFrame = None,
    order_first: bool = False,
):
    nodes_direction, edges_direction, new_outflow_edges = (
        find_node_edge_ids_in_directed_graph(
            from_node_ids=edges.node_start.to_numpy(),
            to_node_ids=edges.node_end.to_numpy(),
            edge_ids=edges.code.to_numpy(),
            outflow_edge_ids=outflow_edge_ids,
            search_node_ids=nodes.nodeID.to_numpy(),
            search_edge_ids=edges.code.to_numpy(),
            border_node_ids=border_node_ids,
            direction=direction,
            split_points=split_points,
            order_first=order_first,
        )
    )
    for edge_id, node_direction, edge_direction in zip(
        outflow_edge_ids, nodes_direction, edges_direction
    ):
        nodes[f"{direction}_edge_{edge_id}"] = False
        nodes.loc[
            nodes["nodeID"].isin(node_direction), f"{direction}_edge_{edge_id}"
        ] = True
        edges[f"{direction}_edge_{edge_id}"] = False
        edges.loc[edges["code"].isin(edge_direction), f"{direction}_edge_{edge_id}"] = (
            True
        )
    return nodes, edges, new_outflow_edges


def calculate_angles_of_edges_at_nodes(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    nodes_id_column: str = "nodeID",
):
    """Calculates the angles of the upstream and downstream edges for each node. 

    Returns
    -------
    self.nodes: gpd.GeoDataFrame
        Geodataframe containing nodes between waterlines, including upstream and downstream edges and their angles
    """
    edges["upstream_angle"] = edges["geometry"].apply(
        lambda x: calculate_angle(x, "upstream").round(2)
    )
    edges["downstream_angle"] = edges["geometry"].apply(
        lambda x: calculate_angle(x, "downstream").round(2)
    )

    def group_angles(x):
        if isinstance(x, float):
            if np.isnan(x):
                return ""
            else:
                return str(x)
        elif (isinstance(x[0], float) and not np.isnan(x[0])):
            return ",".join([str(a) for a in x])
        else:
            return ""

    for direction, opp_direction in zip(
        ["upstream", "downstream"], ["downstream", "upstream"]
    ):
        node_end = "node_end" if direction == "upstream" else "node_start"
        temp = nodes.merge(
            edges[[node_end, f"{opp_direction}_angle"]].rename(
                columns={node_end: nodes_id_column}
            ),
            how="left",
            on=nodes_id_column,
        )
        temp = temp.groupby(nodes_id_column)
        temp = temp.agg({f"{opp_direction}_angle": list})
        nodes[f"{direction}_angles"] = temp
        nodes[f"{direction}_angles"] = nodes[f"{direction}_angles"].apply(
            lambda x: group_angles(x)
        )
    return nodes, edges


def calculate_discharges_of_edges_at_nodes(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    nodes_id_column: str = "nodeID",
):
    """Calculates the angles of the upstream and downstream edges for each node. 

    Returns
    -------
    self.nodes: gpd.GeoDataFrame
        Geodataframe containing nodes between waterlines, including upstream and downstream edges and their angles
    """
    if "total_specifieke_afvoer" not in edges.columns:
        return nodes, edges
    
    for direction, opp_direction in zip(
        ["upstream", "downstream"], ["downstream", "upstream"]
    ):
        node_end = "node_end" if direction == "upstream" else "node_start"
        temp = nodes.merge(
            edges[[node_end, "total_specifieke_afvoer"]].rename(
                columns={node_end: nodes_id_column, "total_specifieke_afvoer": f"{direction}_discharge"}
            ),
            how="left",
            on=nodes_id_column,
        )
        temp = temp.groupby(nodes_id_column).agg({f"{direction}_discharge": list})
        temp[f"{direction}_discharge"] = temp[f"{direction}_discharge"].apply(
            lambda x: ",".join([f"{a:.3f}" for a in x if ~np.isnan(a)])
        )
        nodes[f"{direction}_discharges"] = temp[f"{direction}_discharge"]
    return nodes, edges


def select_downstream_upstream_edges_angle(nodes, min_difference_angle: float = 20.0):
    """select the upstream or downstream edge that represents the main channel, based on the smallest angle. When the angle of both edges is too large, no edge is selected.

    Parameters
    ----------
    min_difference_angle : str, optional
        minimum , by default 20.0

    Returns
    -------
    gpd.GeoDataFrame: self.nodes
        Geodataframe containing nodes between waterlines, including the selected upstream and downstream angles
    """
    def select_downstream_upstream_edges_angle_per_node(x, min_difference_angle: float = 20.0):
        upstream_edges = [
            a for a in x["upstream_edges"].split(",") if a != ""
        ]
        downstream_edges = [
            a for a in x["downstream_edges"].split(",") if a != ""
        ]
        upstream_angles = [
            float(a) for a in x["upstream_angles"].split(",") if a != ""
        ]
        downstream_angles = [
            float(a) for a in x["downstream_angles"].split(",") if a != ""
        ]

        angle_differences = [
            [round(abs(au - ad), 2) for ad in downstream_angles]
            for au in upstream_angles
        ]

        smallest_angle1 = None
        smallest_angle2 = None
        selected_upstream_edge = None
        selected_downstream_edge = None
        iteration = 0

        for upstream_edge, upstream_angle_differences in zip(
            upstream_edges, angle_differences
        ):
            for downstream_edge, angle_difference in zip(
                downstream_edges, upstream_angle_differences
            ):
                iteration = iteration + 1
                if smallest_angle1 is None or angle_difference < smallest_angle1:
                    smallest_angle2 = smallest_angle1
                    smallest_angle1 = angle_difference
                    selected_upstream_edge = upstream_edge
                    selected_downstream_edge = downstream_edge
                elif smallest_angle2 is None or angle_difference < smallest_angle2:
                    smallest_angle2 = angle_difference

        x["selected_upstream_edge"] = None
        x["selected_downstream_edge"] = None
        if (
            smallest_angle2 is None
            or smallest_angle1 < smallest_angle2 - min_difference_angle
        ):
            x["selected_upstream_edge"] = selected_upstream_edge
            x["selected_downstream_edge"] = selected_downstream_edge

        if x["no_downstream_edges"] == 1:
            x["selected_downstream_edge"] = downstream_edges[0]
        if x["no_upstream_edges"] == 1:
            x["selected_upstream_edge"] = upstream_edges[0]

        return x

    nodes = nodes.apply(
        lambda x: select_downstream_upstream_edges_angle_per_node(x, min_difference_angle),
        axis=1,
    )
    return nodes


def select_downstream_upstream_edges_discharge(nodes, min_difference_discharge_factor: float = 2.0):
    """select the upstream or downstream edge that represents the main channel, based on the smallest angle. When the angle of both edges is too large, no edge is selected.

    Parameters
    ----------
    min_difference_angle : str, optional
        minimum , by default 20.0

    Returns
    -------
    gpd.GeoDataFrame: self.nodes
        Geodataframe containing nodes between waterlines, including the selected upstream and downstream angles
    """
    def select_downstream_upstream_edges_discharge_per_node(x, min_difference_discharge_factor: float = 2.0):
        upstream_edges = [
            a for a in x["upstream_edges"].split(",") if a != ""
        ]
        downstream_edges = [
            a for a in x["downstream_edges"].split(",") if a != ""
        ]
        upstream_discharges = [
            float(a) for a in x["upstream_discharges"].split(",") if a != ""
        ]
        downstream_discharges = [
            float(a) for a in x["downstream_discharges"].split(",") if a != ""
        ]

        x["selected_upstream_edge"] = None
        x["selected_downstream_edge"] = None

        if len(upstream_discharges) == 1:
            x["selected_upstream_edge"] = upstream_edges[0]
        elif len(upstream_discharges) > 1:
            max_upstream_discharge = max(upstream_discharges)
            upstream_discharges_sort = sorted(upstream_discharges, reverse=True)
            if upstream_discharges_sort[0] >= upstream_discharges_sort[1] * min_difference_discharge_factor:
                index_max = upstream_discharges.index(max_upstream_discharge)
                x["selected_upstream_edge"] = upstream_edges[index_max]

        if len(downstream_discharges) == 1:
            x["selected_downstream_edge"] = downstream_edges[0]
        elif len(downstream_discharges) > 1:
            max_downstream_discharge = max(downstream_discharges)
            downstream_discharges_sort = sorted(downstream_discharges, reverse=True)
            if downstream_discharges_sort[0] >= downstream_discharges_sort[1] * min_difference_discharge_factor:
                index_max = downstream_discharges.index(max_downstream_discharge)
                x["selected_downstream_edge"] = downstream_edges[index_max]

        return x

    nodes = nodes.apply(
        lambda x: select_downstream_upstream_edges_discharge_per_node(x, min_difference_discharge_factor),
        axis=1,
    )
    return nodes


def define_list_upstream_downstream_edges_ids(
    node_ids: np.array,
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    nodes_id_column: str = "nodeID",
    edges_id_column: str = "code",
):
    """Get the upstream and downstream edges for each node. 

    Returns
    -------
    self.nodes: gpd.GeoDataFrame
        Geodataframe containing nodes between waterlines, including upstream and downstream edges
    """
    logging.info("   x find connected edges for nodes")
    nodes_sel = nodes[nodes.nodeID.isin(node_ids)].copy()
    nodes_sel.index = nodes_sel[nodes_id_column].values

    if any(edges[edges_id_column].apply(lambda x: x is None)):
        none_code_edges = edges[edges[edges_id_column].isna()]
        logging.info(f"     - edges found without code: {len(none_code_edges)}")
        edges = edges[edges[edges_id_column].notna()]

    for direction in ["upstream", "downstream"]:
        node_end = "node_end" if direction == "upstream" else "node_start"
        direction_edges = nodes_sel.merge(
            edges[[node_end, "code"]],
            how="left",
            left_on=nodes_id_column,
            right_on=node_end,
        )
        nodes_sel[f"{direction}_edges"] = direction_edges.groupby(nodes_id_column).agg(
            {edges_id_column: list}
        )

        def len_edges(x):
            if ~(isinstance(x[0], float) and np.isnan(x[0])):
                return len(x)
            else: 
                return 0

        nodes_sel[f"no_{direction}_edges"] = nodes_sel[f"{direction}_edges"].apply(
            lambda x: len_edges(x)
        )

        def list_to_str(x):
            if (~(isinstance(x[0], float) and np.isnan(x[0])) and x[0] is not None):
                return ",".join(x) 
            else:
                return  ",".join([])

        nodes_sel[f"{direction}_edges"] = nodes_sel[f"{direction}_edges"].apply(
            lambda x: list_to_str(x)
        )
    nodes_sel = nodes_sel.reset_index(drop=True)
    return nodes_sel


def get_start_edges_nodes(edges, nodes, direction="downstream"):
    """    Function to get the starting edges of the graph.
    This function identifies nodes that do not have"""

    if direction == "downstream":
        node_end = "node_end"
        node_start = "node_start"
    else:
        node_end = "node_start"
        node_start = "node_end"

    # any incoming edges (i.e., nodes that are not the end of any edge)
    nodes_end_ids = edges[node_end].unique()
    start_nodes = nodes[~nodes.index.isin(nodes_end_ids)]
    other_nodes = nodes[nodes.index.isin(nodes_end_ids)]

    # Identify edges that start from the identified start nodes
    nodes_start_ids = start_nodes.index
    start_edges = edges[edges[node_start].isin(nodes_start_ids)]
    other_edges = edges[~edges[node_start].isin(nodes_start_ids)]
    return start_nodes, other_nodes, start_edges, other_edges


def sum_edge_node_values_through_network(
    edges, 
    nodes, 
    edges_nodes: str = "edges", 
    edges_id_column: str = "code", 
    nodes_id_column: str = "nodeID", 
    direction: str = "downstream", 
    column_to_sum: str = "specifieke_afvoer", 
    sum_column: str = "total_specifieke_afvoer",
):
    nodes = nodes.set_index(nodes_id_column)
    edges = edges.set_index(edges_id_column)

    edges[sum_column] = 0.0
    if sum_column not in nodes.columns:
        nodes[sum_column] = 0.0

    other_edges = edges.copy()
    other_nodes = nodes.copy()
    start_edges = edges.copy()

    logging.info(f"   x sum of '{column_to_sum}' '{direction}' through the network (total {len(other_edges)} edges)...")
    iteration = 0
    while not other_edges.empty and not start_edges.empty:
        iteration += 1
        # find all start edges and nodes
        start_nodes, other_nodes, start_edges, other_edges = get_start_edges_nodes(
            other_edges, 
            other_nodes, 
            direction=direction, 
        )
        logging.info(f"     - Found {len(start_edges)} start edges ({len(other_edges)} left)")

        # add the total values of start nodes to start edges
        start_edges = start_edges.drop(
            columns=[sum_column], 
        ).merge(
            start_nodes[sum_column],
            left_on="node_start",
            right_index=True,
            how="left"
        )

        # fillna and multiply with downstream distribution factor
        start_edges[sum_column] = start_edges[sum_column].fillna(0.0)
        start_edges = start_edges[~start_edges.index.duplicated()]

        if "downstream_splits_dist" in start_edges.columns:
            start_edges[sum_column] = (
                start_edges[sum_column] * start_edges["downstream_splits_dist"]
            )
        start_edges[column_to_sum] += start_edges[sum_column]

        edges = edges[~edges.index.duplicated(keep='first')]
        edges.loc[start_edges.index, sum_column] += start_edges[column_to_sum].values

        specifieke_afvoer_start_nodes = start_edges[["node_end", column_to_sum]].groupby("node_end").sum()
        nodes.loc[specifieke_afvoer_start_nodes.index, sum_column] += specifieke_afvoer_start_nodes[column_to_sum].values
        other_nodes = nodes.loc[other_nodes.index]

    edges[sum_column] = edges[sum_column].astype(float)
    nodes[sum_column] = nodes[sum_column].astype(float)
    nodes = nodes.reset_index()
    edges = edges.reset_index()

    return edges, nodes

