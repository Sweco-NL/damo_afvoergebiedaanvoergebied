import logging
import time
from pathlib import Path

from .generator_specifieke_afvoer import GeneratorSpecifiekeAfvoer


def run_generator_specifieke_afvoer(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    waterschap: str = None,
    water_lines: list[str] = ["hydroobject"],
    generate_specifieke_afvoer: bool = True,
    sum_specifieke_afvoer: bool = True,
    level_discharge_units: int = 0, # level discharge units 0/1/2/3
    use_specifieke_afvoer: float = 1.0,
    file_name_specifieke_afvoer: str = None,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorSpecifiekeAfvoer:
    start_time = time.time()
    afvoer = GeneratorSpecifiekeAfvoer(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        water_lines=water_lines,
        read_results=read_results, 
        write_results=write_results,
    )

    if generate_specifieke_afvoer:
        afvoer.generate_distribution_splits_downstream()

        afvoer.read_specifieke_afvoer(
            file_name_specifieke_afvoer=file_name_specifieke_afvoer
        )
        afvoer.add_specifieke_afvoer_to_discharge_units(
            use_specifieke_afvoer=use_specifieke_afvoer,
            level_discharge_units=level_discharge_units
        )
        afvoer.add_specifieke_afvoer_to_edges()
        afvoer.fill_specifieke_afvoer_water_supply_structures()
    
    if sum_specifieke_afvoer:
        afvoer.sum_specifieke_afvoer_through_network()
    
    if create_html_map:
        afvoer.generate_folium_map(open_html=open_html)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return afvoer


