import logging
import time
from pathlib import Path

from .generator_specifieke_afvoer import GeneratorSpecifiekeAfvoer


def run_generator_specifieke_afvoer(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    waterschap: str = None,
    water_lines: list[str] = ["hydroobjecten"],
    level_discharge_units: int = 0, # level discharge units 0/1/2/3
    use_specific_discharge: float = 1.0,
    specific_discharge_file_name: str = None,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorSpecifiekeAfvoer:
    start_time = time.time()
    discharge = GeneratorSpecifiekeAfvoer(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        water_lines=water_lines,
        read_results=read_results, 
        write_results=write_results,
    )

    discharge.generate_distribution_splits_downstream()

    discharge.read_specific_discharge(
        specific_discharge_file_name=specific_discharge_file_name
    )
    discharge.add_specific_discharge_to_discharge_units(
        use_specific_discharge=use_specific_discharge,
        level_discharge_units=level_discharge_units
    )
    discharge.add_specific_discharge_to_edges()
    discharge.sum_specific_discharge_through_network()

    if create_html_map:
        discharge.generate_folium_map(open_html=open_html)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return discharge

