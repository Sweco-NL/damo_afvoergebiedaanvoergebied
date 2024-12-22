from pathlib import Path
import geopandas as gpd

from .generator_order_levels import GeneratorOrderLevels


def run_generator_order_levels(
    path: Path,
    waterschap: str,
    range_order_code_min: int,
    range_order_code_max: int,
    buffer_rws: float = 10.0,
    generate_order_no: float = True,
    generate_order_code: float = True,
    order_for_each_edge: bool = True,
    water_lines: list[str] = None,
    read_results: bool = True,
    write_results: bool = True,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorOrderLevels:
    order = GeneratorOrderLevels(
        path=path,
        waterschap=waterschap,
        range_order_code_min=range_order_code_min,
        range_order_code_max=range_order_code_max,
        read_results=read_results,
        write_results=write_results,
    )

    if generate_order_no:
        order.create_graph_from_network(water_lines=water_lines)

        order.define_list_upstream_downstream_edges_ids()

        order.calculate_angles_of_edges_at_nodes()

        order.select_downstream_upstream_edges(min_difference_angle=20.0)

        order.find_end_points_hydroobjects()

        order.generate_rws_code_for_all_outflow_points(buffer_rws=buffer_rws)

        order.generate_orde_level_for_hydroobjects()

    if generate_order_code:
        order.generate_order_code_for_edges(order_for_each_edge=order_for_each_edge)

    if write_results:
        order.export_results_to_gpkg()

    if create_html_map:
        order.generate_folium_map(open_html=open_html)
    return order
