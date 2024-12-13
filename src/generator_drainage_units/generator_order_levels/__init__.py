from pathlib import Path
import geopandas as gpd

from .generator_order_levels import GeneratorOrderLevels
from ..utils.general_functions import (
    define_list_upstream_downstream_edges_ids,
    calculate_angles_of_edges_at_nodes,
)


def run_generator_order_levels(
    path: Path,
    waterschap: str,
    range_orde_code_min: int,
    range_orde_code_max: int,
    water_lines: list[str] = None,
    read_results: bool = True,
    write_results: bool = True,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorOrderLevels:
    order = GeneratorOrderLevels(
        path=path,
        waterschap=waterschap,
        range_orde_code_min=range_orde_code_min,
        range_orde_code_max=range_orde_code_max,
        read_results=read_results,
        write_results=write_results,
    )
    order.create_graph_from_network(
        water_lines=water_lines
    )
    order.create_graph_from_network(water_lines=water_lines)
    order.nodes = define_list_upstream_downstream_edges_ids(
        node_ids=order.nodes.nodeID.values, nodes=order.nodes, edges=order.edges
    )
    order.nodes, order.edges = calculate_angles_of_edges_at_nodes(
        nodes=order.nodes, edges=order.edges
    )
    order.find_end_points_hydroobjects()
    order.generate_rws_code_for_all_outflow_points()
    order.generate_order_levels_to_outflow_nodes_edges(max_order=None)

    if create_html_map:
        order.generate_folium_map(
            open_html=open_html
        )
    return order
