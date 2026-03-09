import logging
import time
from pathlib import Path

from .generator_duikers import GeneratorDuikers


def run_generator_duikers(
    path: Path,
    create_html_map: bool = False,
    open_html_map: bool = False,
) -> GeneratorDuikers:
    """Run Generator Culvert Locations (Duikergenerator)

    Parameters
    ----------
    path : Path
        Path to the case directory. Directory name is used as name for the case
    settings : Dict
        Dict with all settings and locations and names files
    read_results : bool, optional
        option (True/False) to read previous results from gpkg, by default None
    write_results : bool, optional
        option (True/False) to write results to case folder in gpkg, by default None

    Returns
    -------
    GeneratorDuikers
        An instance of the GeneratorDuikers including basisdata and results
    """
    start_time = time.time()
    culvert = GeneratorDuikers(path=path)

    # generate all vertices every 10 meters
    culvert.generate_vertices_along_waterlines()

    # generate all potential culverts with max lenght = 40
    culvert.find_potential_culvert_locations()

    # check intersections culvert with objects and first filter
    culvert.check_intersections_potential_culverts()

    # assing scores to potential culverts
    culvert.assign_scores_to_potential_culverts()

    # select best culverts
    culvert.select_correct_score_based_on_score_and_length()

    # post processing culverts
    culvert.post_process_potential_culverts()

    # split hydroobjects at endpoints culverts
    culvert.splits_hydroobject_by_endpoints_of_culverts_and_combine()

    # check if culverts are in the right direction
    culvert.check_culverts_direction()

    culvert.combine_culvert_with_line()

    culvert.splits_hydroobject_by_endpoind_of_culverts_and_combine_2()

    culvert.get_shortest_path_from_overige_watergang_to_hydroobjects()

    # create map
    if create_html_map:
        culvert.generate_folium_map(open_html_map=open_html_map)

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return culvert
