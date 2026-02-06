from pathlib import Path
import geopandas as gpd
import logging
import time
from .generator_gebiedsorde import GeneratorGebiedsOrde


def run_generator_gebiedsorde(
    path: Path,
    dir_basisdata: str = "0_basisdata",
    dir_results: str = "1_resultaat",
    waterschap: str = None,
    generate_new_outflow_nodes_hydro: bool = False,
    search_range_outflow_nodes_hydro: float = 50.0,
    generate_order_no: bool = True,
    max_order_no: int = 100,
    generate_order_code: bool = True,
    generate_order_code_overig: bool = False,
    order_for_each_edge: bool = True,
    read_results: bool = True,
    write_results: bool = True,
    create_html_map: bool = False,
    open_html_map: bool = False,
    snapping_distance: float = 0.05,
) -> GeneratorGebiedsOrde:

    start_time = time.time()
    order = GeneratorGebiedsOrde(
        path=path,
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        read_results=read_results,
        write_results=write_results,
        snapping_distance=snapping_distance,
    )

    if generate_order_no:
        if order.edges is None:
            order.create_graph_from_network()

        if not generate_new_outflow_nodes_hydro and order.outflow_nodes_hydro is not None:
            order.read_outflow_nodes_hydro_with_rws_code(
                buffer_outflow_nodes_hydro=search_range_outflow_nodes_hydro
            )
        else:
            order.generate_rws_code_for_outflow_points(
                search_range_outflow_nodes_hydro=search_range_outflow_nodes_hydro
            )

        order.generate_order_level_for_hydroobjects(max_order_no=max_order_no)

    if generate_order_code:
        order.generate_order_code_for_hydroobjects(
            order_for_each_edge=order_for_each_edge
        )

    if generate_order_code_overig:
        order.generate_order_no_order_code_for_other_waterlines()

    if create_html_map:
        order.generate_folium_map(open_html_map=open_html_map)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return order
