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
        #     logging.debug("xxxxxxxxxxxx")
        #     logging.debug(edge_id)
        #     logging.debug(from_node)
        #     logging.debug(pred_node)
        #     logging.debug(pred_edge)
        #     logging.debug(split_node_edge_ids2)
        #     logging.debug(split_node_edge_ids2[from_node])
        #     logging.debug("xxxxxxxxxxxx")

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
    len_outflow_edge_ids = np.shape(outflow_edge_ids)[0]
    results_nodes = [
        [int(to_node_ids[np.where(edge_ids == e)][0])] for e in outflow_edge_ids
    ]
    results_edges = [[e] if e in search_edge_ids else [] for e in outflow_edge_ids]
    if set_logging:
        logging.debug(
            f"    - find {direction} nodes/edges for {len(outflow_edge_ids)}/{len(search_edge_ids)} nodes"
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
            # print(f" * {i+1}/{len_outflow_node_ids} ({(i+1)/len_outflow_node_ids:.2%})", end="\r")
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
    edges["upstream_angle"] = edges["geometry"].apply(
        lambda x: calculate_angle(x, "upstream").round(2)
    )
    edges["downstream_angle"] = edges["geometry"].apply(
        lambda x: calculate_angle(x, "downstream").round(2)
    )
    for direction, opp_direction in zip(
        ["upstream", "downstream"], ["downstream", "upstream"]
    ):
        node_end = "node_end" if direction == "upstream" else "node_start"
        nodes[f"{direction}_angles"] = (
            nodes.merge(
                edges[[node_end, f"{opp_direction}_angle"]].rename(
                    columns={node_end: nodes_id_column}
                ),
                how="left",
                on=nodes_id_column,
            )
            .groupby(nodes_id_column)
            .agg({f"{opp_direction}_angle": list})
        )
        nodes[f"{direction}_angles"] = nodes[f"{direction}_angles"].apply(
            lambda x: ",".join([str(a) for a in x])
            if ~(isinstance(x[0], float) and np.isnan(x[0]))
            else ",".join([])
        )
    return nodes, edges


def select_downstream_upstream_edges(nodes, min_difference_angle: str = 20.0):
    def select_downstream_upstream_edges_per_node(x, min_difference_angle: str = 20.0):
        upstream_edges = x["upstream_edges"] = [
            a for a in x["upstream_edges"].split(",") if a != ""
        ]
        downstream_edges = x["downstream_edges"] = [
            a for a in x["downstream_edges"].split(",") if a != ""
        ]
        upstream_angles = x["upstream_angles"] = [
            float(a) for a in x["upstream_angles"].split(",") if a != ""
        ]
        downstream_angles = x["downstream_angles"] = [
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
            x["selected_downstream_edge"] = x["downstream_edges"][0]
        if x["no_upstream_edges"] == 1:
            x["selected_upstream_edge"] = x["upstream_edges"][0]

        return x

    nodes = nodes.apply(
        lambda x: select_downstream_upstream_edges_per_node(x, min_difference_angle),
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
    logging.info("   x find connected edges for nodes")
    nodes_sel = nodes[nodes.nodeID.isin(node_ids)].copy()
    nodes_sel.index = nodes_sel[nodes_id_column].values
    for direction in ["upstream", "downstream"]:
        node_end = "node_end" if direction == "upstream" else "node_start"
        direction_edges = nodes_sel.merge(
            edges[[node_end, "code"]],
            how="left",
            left_on=nodes_id_column,
            right_on=node_end,
        )
        nodes_sel[f"{direction}_edges"] = (
            direction_edges
            .groupby(nodes_id_column)
            .agg({edges_id_column: list})
        )
        nodes_sel[f"no_{direction}_edges"] = nodes_sel[f"{direction}_edges"].apply(
            lambda x: len(x)
            if ~(isinstance(x[0], float) and np.isnan(x[0]))
            else 0
        )
        nodes_sel[f"{direction}_edges"] = nodes_sel[f"{direction}_edges"].apply(
            lambda x: ",".join(x)
            if ~(isinstance(x[0], float) and np.isnan(x[0]))
            else ",".join([])
        )
    nodes_sel = nodes_sel.reset_index(drop=True)
    return nodes_sel
