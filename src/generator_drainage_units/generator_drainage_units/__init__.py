
import logging
import time
from pathlib import Path
import xarray as xr

from .generator_drainage_units import GeneratorDrainageUnits


def run_generator_drainage_units(
    path: Path,
    dir_basis_data: Path,
    dir_inter_results: Path = None,
    dir_results: Path = None,
    ghg_file_name: str = None,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
    water_lines: list[str] = ["hydroobjecten"]
) -> GeneratorDrainageUnits:
    """Run Generator Culvert Locations (Duikergenerator)

    Parameters
    ----------
    path : Path
        Path to the case directory including directories 0_basisdata and
        1_tussenresultaat. Directory name is used as name for the case,
        by default None
    read_results : bool, optional
        option (True/False) to read previous results from gpkg, by default None
    write_results : bool, optional
        option (True/False) to write results to case folder in gpkg, by default None

    Returns
    -------
    GeneratorCulvertLocations
        An instance of the GeneratorCulvertLocations including basisdata and results
    """
    start_time = time.time()
    gdu = GeneratorDrainageUnits(
        path=path, 
        read_results=read_results, 
        write_results=write_results
    )
    gdu.ghg = xr.open_dataset(Path(gdu.path, gdu.dir_basis_data, ghg_file_name))["__xarray_dataarray_variable__"][0]
    gdu.ghg.name = "GHG_2000-2010_L1"

    # create map
    if create_html_map:
        gdu.generate_folium_map()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return gdu
