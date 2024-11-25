import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import logging


def find_predecessors_graph(
    from_node_ids, to_node_ids, node_id, border_node_ids=None, preds=np.array([])
):
    """
    Find predecessors within graph for specified node_id.

    Note: recursive function!
    """
    pred = from_node_ids[np.where(to_node_ids == node_id)]
    for i in range(pred.shape[0]):
        p = pred[i]
        if p not in preds:
            preds = np.append(preds, p)
            if border_node_ids is None or p not in border_node_ids:
                preds = find_predecessors_graph(
                    from_node_ids, to_node_ids, p, border_node_ids, preds
                )
    return preds


def find_predecessors_graph_with_splits(
    from_node_ids, to_node_ids, edge_ids, node_id, border_node_ids=None, split_node_edge_ids=None, preds=np.array([])
):
    """
    Find predecessors within graph for specified node_id.

    Note: recursive function!
    """
    pred = from_node_ids[np.where(to_node_ids == node_id)]
    edge = edge_ids[np.where(to_node_ids == node_id)]

    for i in range(pred.shape[0]):
        p = pred[i]
        e = edge[i]
        if p not in preds:
            preds = np.append(preds, p)
            if border_node_ids is None or p not in border_node_ids:
                if split_node_edge_ids is None or p not in split_node_edge_ids or split_node_edge_ids[p] == e:
                    preds = find_predecessors_graph_with_splits(
                        from_node_ids, to_node_ids, edge_ids, p, border_node_ids, split_node_edge_ids, preds
                    )
    return preds


def accumulate_values_graph(
    from_node_ids,
    to_node_ids,
    node_ids,
    values_node_ids,
    values,
    border_node_ids=None,
    itself=False,
    direction="upstream",
    decimals=None,
):
    """Calculate for all node_ids the accumulated values of all predecessors with values."""
    len_node_ids = np.shape(node_ids)[0]
    results = np.zeros(np.shape(node_ids))
    logging.info(f"accumulate values using graph for {len(node_ids)} node(s)")
    for i in range(node_ids.shape[0]):
        print(f" * {i+1}/{len_node_ids} ({(i+1)/len(node_ids):.2%})", end="\r")
        node_id = node_ids[i]
        if direction == "upstream":
            pred = find_predecessors_graph(
                from_node_ids, to_node_ids, node_id, border_node_ids, np.array([])
            )
        else:
            pred = find_predecessors_graph(
                to_node_ids, from_node_ids, node_id, border_node_ids, np.array([])
            )
        if itself:
            pred = np.append(pred, node_id)
        pred_sum = np.sum(values[np.searchsorted(values_node_ids, pred)])
        if decimals is None:
            results[i] = pred_sum
        else:
            results[i] = np.round(pred_sum, decimals=decimals)
    return results


def find_node_ids_in_directed_graph(
    from_node_ids,
    to_node_ids,
    edge_ids,
    node_ids,
    search_node_ids,
    border_node_ids=None,
    direction="upstream",
    split_points=None,
):
    """Find node ids from a certain list whether they are present between predecessors.

    Parameters
    ----------
    from_node_ids : _type_
        _description_
    to_node_ids : _type_
        _description_
    node_ids : _type_
        _description_
    search_node_ids : _type_
        _description_
    border_node_ids : _type_, optional
        _description_, by default None
    direction : str, optional
        _description_, by default "upstream"

    Returns
    -------
    List of Lists
        for each node_id in node_ids a list of all upstream or downstream located node_ids
    """
    len_node_ids = np.shape(node_ids)[0]
    results = []
    logging.debug(
        f"   - find {direction} nodes/edges for {len(node_ids)}/{len(search_node_ids)} nodes"
    )
    search_direction = "upstream" if direction == "downstream" else "downstream"
    
    if split_points is None:
        for i in range(node_ids.shape[0]):
            print(f" * {i+1}/{len_node_ids} ({(i+1)/len(node_ids):.2%})", end="\r")
            node_id = node_ids[i]
            if direction == "upstream":
                pred = find_predecessors_graph(
                    from_node_ids, to_node_ids, node_id, border_node_ids, np.array([])
                )
            else:
                pred = find_predecessors_graph(
                    to_node_ids, from_node_ids, node_id, border_node_ids, np.array([])
                )
            results += [[p for p in pred if p in search_node_ids]]
        
    else:
        split_node_edge_ids = split_points.set_index("nodeID")[
            f"selected_{search_direction}_edge"
        ].to_dict()
        split_node_edge_ids = {k:v for k, v in split_node_edge_ids.items() if v not in [None, '']}


        for i in range(node_ids.shape[0]):
            print(f" * {i+1}/{len_node_ids} ({(i+1)/len(node_ids):.2%})", end="\r")
            node_id = node_ids[i]
            if direction == "upstream":
                pred = find_predecessors_graph_with_splits(
                    from_node_ids, to_node_ids, edge_ids, node_id, 
                    border_node_ids, split_node_edge_ids, np.array([])
                )
            else:
                pred = find_predecessors_graph_with_splits(
                    to_node_ids, from_node_ids, edge_ids, node_id, 
                    border_node_ids, split_node_edge_ids, np.array([])
                )
            results += [[p for p in pred if p in search_node_ids]]

    
    return results


def find_nodes_edges_for_direction(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    node_ids: list,
    border_node_ids: list = None,
    direction: str = "upstream",
    split_points: gpd.GeoDataFrame = None
):
    """Find nodes edges upstream or downstream (direction)

    Parameters
    ----------
    nodes : gpd.GeoDataFrame
        nodes of network
    edges : gpd.GeoDataFrame
        edges of network
    node_ids : list
        node_ids to start searching from
    border_node_ids : list, optional
        list of border_node_ids, by default None
    direction : str, optional
        search direction "upstream" or "downstream", by default "upstream"
    split_points: gpd.GeoDataFrame
        selected split_points with column selected_downstream_edge and selected_upstream_edge

    Returns
    -------
    _type_
        _description_
    """
    nodes_direction = find_node_ids_in_directed_graph(
        from_node_ids=edges.node_start.to_numpy(),
        to_node_ids=edges.node_end.to_numpy(),
        edge_ids=edges.code.to_numpy(),
        node_ids=node_ids,
        search_node_ids=nodes.nodeID.to_numpy(),
        border_node_ids=border_node_ids,
        direction=direction,
        split_points=split_points,
    )

    for node_id, node_direction in zip(node_ids, nodes_direction):
        nodes[f"{direction}_node_{node_id}"] = False
        nodes.loc[[int(n) for n in node_direction], f"{direction}_node_{node_id}"] = (
            True
        )
        edges[f"{direction}_node_{node_id}"] = False
        if direction == "upstream":
            edges.loc[
                edges.node_end.isin(node_direction + [node_id])
                & edges.node_start.isin(node_direction),
                f"{direction}_node_{node_id}",
            ] = True
        else:
            edges.loc[
                edges.node_start.isin(node_direction + [node_id])
                & edges.node_end.isin(node_direction),
                f"{direction}_node_{node_id}",
            ] = True
    return nodes, edges
