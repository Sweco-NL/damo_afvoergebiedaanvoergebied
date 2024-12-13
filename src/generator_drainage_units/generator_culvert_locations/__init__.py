import logging
import time
from pathlib import Path

from .generator_culverts_locations import GeneratorCulvertLocations


def run_generator_culvert_locations(
    path: Path,
    distance_vertices: float = 10.0,
    max_culvert_length: float = 40.0,
    read_results: bool = False,
    write_results: bool = False,
    create_html_map: bool = False,
) -> GeneratorCulvertLocations:
    """Run Generator Culvert Locations (Duikergenerator)

    Parameters
    ----------
    path : Path
        Path to the case directory including directories 0_basisdata and
        1_tussenresultaat. Directory name is used as name for the case,
        by default None
    distance_vertices : float, optional
        distacne between vertices, by default 10.0
    max_culvert_length : int, optional
        maximum culvert length: in case of larger distance between points, connections are not made, by default 40
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
    culverts_generator = GeneratorCulvertLocations(
        path=path,
        read_results=read_results, 
        write_results=write_results
    )
    
    # generate all vertices every 10 meters
    culverts_generator.generate_vertices_along_waterlines(
        distance_vertices=distance_vertices, write_results=write_results
    )

    # generate all potential culverts with max lenght = 40
    culverts_generator.find_potential_culvert_locations(
        max_culvert_length=max_culvert_length, write_results=write_results
    )

    # check intersections culvert with objects and first filter
    culverts_generator.check_intersections_potential_culverts()

    # assing scores to potential culverts
    culverts_generator.assign_scores_to_potential_culverts()

    # select best culverts
    culverts_generator.select_correct_score_based_on_score_and_length()

    # post processing culverts
    culverts_generator.post_process_potential_culverts()

    # split hydroobjects at endpoints culverts
    culverts_generator.splits_hydroobjecten_by_endpoints_of_culverts_and_combine()

    # check if culverts are in the right direction
    culverts_generator.check_culverts_direction()

    culverts_generator.combine_culvert_with_line()

    # check if culverts are in the right direction
    culverts_generator.check_culverts_direction()

    culverts_generator.combine_culvert_with_line()

    culverts_generator.splits_hydroobjecten_by_endpoind_of_culverts_and_combine_2()

    # create map
    if create_html_map:
        culverts_generator.generate_folium_map()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return culverts_generator
