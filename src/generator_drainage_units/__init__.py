import logging
import time
from pathlib import Path

from .generator_culvert_locations import GeneratorCulvertLocations
from .generator_order_levels import GeneratorOrderLevels


def run_generator_culvert_locations(
    path: Path,
    distance_vertices: float = 10.0,
    max_culvert_length: float = 40.0,
    read_results: bool = False,
    write_results: bool = False,
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
        read_results=read_results, write_results=write_results
    )

    # read basis data from folder 0_basisdata
    culverts_generator.read_data_from_case(path=path)

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

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return culverts_generator
