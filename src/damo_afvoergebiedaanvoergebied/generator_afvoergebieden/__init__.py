import logging
import time
from pathlib import Path

from .generator_afvoergebieden import GeneratorAfvoergebieden


def run_generator_afvoergebieden(
    path: Path,
    preprocess: bool = False,
    process: bool = False,
    postprocess1: bool = False,
    postprocess2: bool = False,
    create_html_map: bool = False,
    open_html_map: bool = False,
) -> GeneratorAfvoergebieden:
    start_time = time.time()
    gdu = GeneratorAfvoergebieden(path=path)
    gdu.read_topo()

    if preprocess:
        if gdu.edges is None:
            gdu.create_graph_from_network()
        gdu.preprocess_topo()

    if process:
        gdu.generate_afvoergebied()

    if postprocess1:
        gdu.aggregate_afvoergebied_tot_stroomgebied_1()

    if gdu.afvoergebied_2_gdf is not None and postprocess2:
        gdu.aggregate_afvoergebied_tot_stroomgebied_2()

    if create_html_map:
        gdu.generate_folium_map(open_html_map=open_html_map)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return gdu
