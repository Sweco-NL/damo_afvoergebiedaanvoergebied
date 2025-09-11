import logging
import time
from pathlib import Path

from .generator_specific_discharge import GeneratorSpecificDischarge


def run_generator_specific_discharge(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    waterschap: str = None,
    water_lines: list[str] = ["hydroobjecten"],
    specific_discharge_file_name: str = None,
    calculate_specific_discharge: bool = True,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
    open_html: bool = False,
) -> GeneratorSpecificDischarge:
    start_time = time.time()
    discharge = GeneratorSpecificDischarge(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        water_lines=water_lines,
        read_results=read_results, 
        write_results=write_results,
    )

    if specific_discharge_file_name is not None:
        discharge.read_specific_discharge(ghg_file_name=specific_discharge_file_name)

    if calculate_specific_discharge:
        discharge.add_specific_discharge_to_edges()
        discharge.generate_distribution_splits_downstream()
        discharge.sum_specific_discharge_through_network()

    if create_html_map:
        discharge.generate_folium_map(open_html=open_html)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return discharge

