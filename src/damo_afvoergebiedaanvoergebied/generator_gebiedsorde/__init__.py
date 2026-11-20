from pathlib import Path
import geopandas as gpd
import logging
import time
from .generator_gebiedsorde import GeneratorGebiedsOrde


def run_generator_gebiedsorde(
    path: Path,
    generate_order_no: bool = True,
    generate_order_code: bool = True,
    generate_order_code_overig: bool = True,
    create_html_map: bool = False,
    open_html_map: bool = False,
    generate_new_outflow_nodes: bool = False,
    search_range_outflow_nodes: float = 50.0, 
    buffer_outflow_nodes: float = 50.0,
    max_order_no: int = 1000,
    order_for_each_edge: bool = True,
    html_file_name: str = ""
) -> GeneratorGebiedsOrde:

    start_time = time.time()
    order = GeneratorGebiedsOrde(path=path)

    if generate_order_no:
        logging.info("generate_order_no")
        order.generate_outflow_nodes(
            generate_new_outflow_nodes=generate_new_outflow_nodes,
            search_range_outflow_nodes=search_range_outflow_nodes, 
            buffer_outflow_nodes=buffer_outflow_nodes
        )
        order.generate_order_level_for_hydroobjects(
            max_order_no=max_order_no
        )

    if generate_order_code:
        logging.info("generate_order_code")
        order.generate_order_code_for_hydroobjects(
            order_for_each_edge=order_for_each_edge
        )

    if generate_order_code_overig:
        logging.info("generate_order_code_overig")
        order.generate_order_no_order_code_for_other_waterlines()

    if create_html_map:
        logging.info("generate_folium_map")
        order.generate_folium_map(
            html_file_name=html_file_name,
            open_html_map=open_html_map
        )

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return order
