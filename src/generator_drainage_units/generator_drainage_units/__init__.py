import logging
import time
from pathlib import Path
import rioxarray

from .generator_drainage_units import GeneratorDrainageUnits


def run_generator_drainage_units(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    ghg_file_name: str = None,
    preprocess: bool = True,
    process: bool = True,
    resolution: float = 2.0,
    depth_waterways: float = 1.0,
    buffer_waterways: float = 2.5,
    smooth_distance: float = 25.0,
    iterations: int = 2000,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
) -> GeneratorDrainageUnits:
    """Run Generator Culvert Locations (Duikergenerator)

    Parameters
    ----------
    path : Path
        Path to the case directory. Directory name is used as name for the case
    dir_basisdata : str | pathlib.Path
        String representing subfolder with basisdata
    dir_results : str | pathlib.Path
        String representing subfolder with results
    read_results : bool, optional
        option (True/False) to read previous results from gpkg, by default None
    write_results : bool, optional
        option (True/False) to write results to case folder in gpkg, by default None

    Returns
    -------
    GeneratorDrainageUnits
        An instance of the GeneratorDrainageUnits including basisdata and results
    """
    start_time = time.time()
    gdu = GeneratorDrainageUnits(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        read_results=read_results, 
        write_results=write_results,
    )
    if ghg_file_name is not None:
        gdu.read_ghg(ghg_file_name=ghg_file_name)

        if preprocess:
            gdu.preprocess_ghg(
                resolution=resolution, 
                depth_waterways=depth_waterways,
                buffer_waterways=buffer_waterways,
                smooth_distance=smooth_distance,
            )
        if process:
            gdu.generate_drainage_units(iterations=iterations)

    # create map
    if create_html_map:
        gdu.generate_folium_map()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return gdu
