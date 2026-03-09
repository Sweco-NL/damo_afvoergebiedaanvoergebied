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
) -> GeneratorGebiedsOrde:

    start_time = time.time()
    order = GeneratorGebiedsOrde(path=path)

    if generate_order_no:
        order.generate_outflow_nodes()
        order.generate_order_level_for_hydroobjects()

    if generate_order_code:
        order.generate_order_code_for_hydroobjects()

    if generate_order_code_overig:
        order.generate_order_no_order_code_for_other_waterlines()

    if create_html_map:
        order.generate_folium_map(open_html_map=open_html_map)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return order
