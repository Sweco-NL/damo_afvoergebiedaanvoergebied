from pathlib import Path
import geopandas as gpd
import logging
import time
from .generator_order_levels import GeneratorOrderLevels


def run_generator_order_levels(
    path: Path,
    dir_basisdata: str = "0_basisdata",
    dir_results: str = "1_resultaat",
    waterschap: str = None,
    generate_new_outflow_nodes: bool = False,
    search_range_outflow_nodes: float = 50.0,
    generate_order_no: bool = True,
    generate_order_code: bool = True,
    generate_order_code_sub_waterlines: bool = False,
    order_for_each_edge: bool = True,
    water_lines: list[str] = None,
    read_results: bool = True,
    write_results: bool = True,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorOrderLevels:
    start_time = time.time()
    order = GeneratorOrderLevels(
        path=path,
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        read_results=read_results,
        write_results=write_results,
    )
    order.read_required_data_from_case()

    if generate_order_no:
        order.create_graph_from_network(water_lines=water_lines)
        order.define_list_upstream_downstream_edges_ids()
        order.calculate_angles_of_edges_at_nodes()
        order.select_downstream_upstream_edges(min_difference_angle=20.0)

        if not generate_new_outflow_nodes and order.outflow_nodes is not None:
            order.read_outflow_nodes_with_rws_code(
                buffer_outflow_nodes=search_range_outflow_nodes
            )
        else:
            order.generate_rws_code_for_outflow_points(
                search_range_outflow_nodes=search_range_outflow_nodes
            )

        order.generate_order_level_for_hydroobjects()

    if generate_order_code:
        order.generate_order_code_for_hydroobjects(
            order_for_each_edge=order_for_each_edge
        )

    if generate_order_code_sub_waterlines:
        order.generate_order_no_order_code_for_other_waterlines()

    if create_html_map:
        order.generate_folium_map(open_html=open_html)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return order
