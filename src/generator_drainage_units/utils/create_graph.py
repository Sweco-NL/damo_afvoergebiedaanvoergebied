import pandas as pd
import geopandas as gpd
import networkx as nx
import momepy
from shapely.geometry import Point
import logging
import logging


def create_graph_from_edges(edges: gpd.GeoDataFrame):
    edges["geom_length"] = edges.geometry.length
    G = momepy.gdf_to_nx(
        edges,
        approach="primal",
        directed=True,
        length="geom_length",
    )
    nodes, edges = momepy.nx_to_gdf(G)
    return nodes, edges, G


def generate_nodes_from_edges(
    edges: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Generate start/end nodes from edges and update node information in edges GeoDataFrame.
    Return updated edges geodataframe and nodes geodataframe

    Parameters
    ----------
    edges : gpd.GeoDataFrame
        Line feature dataset containing edges

    Returns
    -------
    Tuple containing GeoDataFrame with edges and GeoDataFrame with nodes
    """
    edges["edge_no"] = range(len(edges))
    edges.index = edges["edge_no"].values

    # Generate nodes from edges and include extra information in edges
    edges[["from_node", "to_node"]] = [
        [g.coords[0], g.coords[-1]] for g in edges.geometry
    ]  # generate endpoints
    _nodes = pd.unique(
        edges["from_node"].tolist() + edges["to_node"].tolist()
    )  # get unique nodes
    indexer = dict(zip(_nodes, range(len(_nodes))))
    nodes = gpd.GeoDataFrame(
        data={"node_no": [indexer[x] for x in _nodes]},
        index=[indexer[x] for x in _nodes],
        geometry=[Point(x) for x in _nodes],
        crs=edges.crs,
    )
    edges[["from_node", "to_node"]] = edges[["from_node", "to_node"]].map(
        indexer.get
    )  # get node id instead of coords
    return edges, nodes


def create_graph_based_on_nodes_edges(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    directional_graph: bool = True,
    add_edge_length_as_weight: bool = False,
    print_logmessage: bool = True,
) -> nx.Graph | nx.DiGraph:
    """
    create networkx graph based on geographic nodes and edges.
    default a directional graph.
    """
    if directional_graph:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()
    if nodes is not None:
        for i, node in nodes.iterrows():
            graph.add_node(node.node_no, pos=(node.geometry.x, node.geometry.y))
    if edges is not None:
        for i, edge in edges.iterrows():
            if add_edge_length_as_weight:
                graph.add_edge(
                    edge.from_node, edge.to_node, weight=edge.geometry.length
                )
            else:
                graph.add_edge(edge.from_node, edge.to_node)
    if print_logmessage:
        logging.info(
            f" - create network graph from nodes ({len(nodes)}x) and edges ({len(edges)}x)"
        )
    return graph


def add_basin_code_from_network_to_nodes_and_edges(
    graph: nx.DiGraph,
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
):
    """add basin (subgraph) code to nodes and edges"""
    subgraphs = list(nx.weakly_connected_components(graph))
    if nodes is None or edges is None:
        return None, None
    nodes["basin"] = -1
    edges["basin"] = -1
    for i, subgraph in enumerate(subgraphs):
        node_ids = list(subgraph)
        edges.loc[
            edges["from_node"].isin(node_ids) & edges["to_node"].isin(node_ids),
            "basin",
        ] = i + 1
        nodes.loc[nodes["node_no"].isin(list(subgraph)), "basin"] = i + 1
    loggin.info(
        f" - define numbers Ribasim-Basins ({len(subgraphs)}x) and join edges/nodes"
    )
    return nodes, edges
