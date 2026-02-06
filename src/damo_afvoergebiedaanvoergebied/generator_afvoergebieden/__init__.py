import logging
import time
from pathlib import Path

from .generator_afvoergebieden import GeneratorAfvoergebieden


def run_generator_afvoergebieden(
    path: Path,
    dir_basisdata: Path,
    dir_results: Path = None,
    waterschap: str = None,
    method: str = "pyflwdir",
    ghg_file_name: str = None,
    preprocess: bool = False,
    process: bool = False,
    postprocess1: bool = False,
    postprocess2: bool = False,
    resolution: float = 2.0,
    depth_waterways: float = 1.0,
    buffer_waterways: float = 2.5,
    smooth_distance: float = 25.0,
    afvoergebied_cmap: str = "Pastel2",
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
) -> GeneratorAfvoergebieden:
    start_time = time.time()
    gdu = GeneratorAfvoergebieden(
        path=path, 
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
        waterschap=waterschap,
        read_results=read_results, 
        write_results=write_results,
        method=method,
    )
    if ghg_file_name is not None:
        gdu.read_ghg(ghg_file_name=ghg_file_name)

        if preprocess:
            if gdu.edges is None:
                gdu.create_graph_from_network()
            
            gdu.preprocess_ghg(
                resolution=resolution, 
                depth_waterways=depth_waterways,
                buffer_waterways=buffer_waterways,
                smooth_distance=smooth_distance,
            )
        if process:
            gdu.generate_afvoergebied()

        if postprocess1:
            gdu.aggregate_afvoergebied_tot_2()

        if gdu.afvoergebied_2_gdf is not None and postprocess2:
            gdu.aggregate_afvoergebied_tot_4()

    if create_html_map:
        gdu.generate_folium_map()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return gdu
