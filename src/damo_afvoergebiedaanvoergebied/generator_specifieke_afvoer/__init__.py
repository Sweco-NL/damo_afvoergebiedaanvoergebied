import logging
import time
from pathlib import Path

from .generator_specifieke_afvoer import GeneratorSpecifiekeAfvoer


def run_generator_specifieke_afvoer(
    path: Path,
    generate_values: bool = True,
    sum_values: bool = True,
    create_html_map: bool = False,
    open_html_map: bool = False,
) -> GeneratorSpecifiekeAfvoer:
    start_time = time.time()
    afvoer = GeneratorSpecifiekeAfvoer(path=path)

    if generate_values:
        afvoer.generate_distribution_splits_downstream()

        afvoer.read_specifieke_afvoer()
        afvoer.add_specifieke_afvoer_to_discharge_units()
        afvoer.add_specifieke_afvoer_to_edges()
        afvoer.fill_specifieke_afvoer_water_supply_structures()
    
    if sum_values:
        afvoer.sum_specifieke_afvoer_through_network()
    
    if create_html_map:
        afvoer.generate_folium_map(open_html_map=open_html_map)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return afvoer


