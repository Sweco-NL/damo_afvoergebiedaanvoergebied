import logging
import time
from pathlib import Path

from .generator_drainage_units import GeneratorDrainageUnits


def run_generator_drainage_units(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    method: str = "pyflwdir",
    ghg_file_name: str = None,
    preprocess: bool = False,
    process: bool = False,
    postprocess: bool = False,
    resolution: float = 2.0,
    depth_waterways: float = 1.0,
    buffer_waterways: float = 2.5,
    smooth_distance: float = 25.0,
    iterations: int = 2000,
    iteration_group: int = 100,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
) -> GeneratorDrainageUnits:
    """_summary_

    _extended_summary_

    Parameters
    ----------
    path : Path
        Path to the case directory. Directory name is used as name for the case
    dir_basisdata : str | pathlib.Path
        String representing subfolder with basisdata
    dir_results : str | pathlib.Path
        String representing subfolder with results
    ghg_file_name : str, optional
        _description_, by default None
    preprocess : bool, optional
        _description_, by default True
    process : bool, optional
        _description_, by default True
    postprocess : bool, optional
        _description_, by default True
    resolution : float, optional
        _description_, by default 2.0
    depth_waterways : float, optional
        _description_, by default 1.0
    buffer_waterways : float, optional
        _description_, by default 2.5
    smooth_distance : float, optional
        _description_, by default 25.0
    iterations : int, optional
        _description_, by default 2000

    Returns
    -------
    _type_
        _description_

    Yields
    ------
    GeneratorDrainageUnits
        _description_
    """
    start_time = time.time()
    gdu = GeneratorDrainageUnits(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        read_results=read_results, 
        write_results=write_results,
        method=method,
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
            gdu.generate_drainage_units(
                iterations=iterations,
                iteration_group=iteration_group
            )

        if postprocess:
            gdu.aggregate_drainage_units()

    if create_html_map:
        gdu.generate_folium_map()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return gdu
