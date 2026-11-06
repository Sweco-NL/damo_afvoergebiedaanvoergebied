import logging
import time
from pathlib import Path

from .generator_duikers import GeneratorDuikers


def _format_duration(duration_seconds: float) -> str:
    minutes, seconds = divmod(int(duration_seconds), 60)
    return f"{minutes:02d}:{seconds:02d}"


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
    start_time = time.perf_counter()
    culvert = GeneratorDuikers(path=path)

    # generate all vertices every 10 meters
    step_start_time = time.perf_counter()
    culvert.generate_vertices_along_waterlines()
    logging.info(
        f"   x Generated vertices along waterlines in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )  
    # generate all potential culverts with max lenght = 40
    step_start_time = time.perf_counter()
    culvert.find_potential_culvert_locations()
    logging.info(
        f"   x Found potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # check intersections culvert with objects and first filter
    step_start_time = time.perf_counter()
    culvert.check_intersections_potential_culverts()
    logging.info(
        f"   x Checked intersections potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # assing scores to potential culverts
    step_start_time = time.perf_counter()
    culvert.assign_scores_to_potential_culverts()
    logging.info(
        f"   x Assigned scores to potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # Calculate lengths and angles
    step_start_time = time.perf_counter()
    culvert.calculate_length_and_angles_potential_culverts_2a()
    logging.info(
        "   x Calculated lengths and angles in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )

    # Filter by length and score group
    step_start_time = time.perf_counter()
    culvert.filter_potential_culverts_2b()
    logging.info(
        "   x Filtered potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )

    # Select the best score within the remaining groups
    step_start_time = time.perf_counter()
    culvert.select_best_potential_culverts_2c()
    logging.info(
        "   x Selected best potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )

    # post processing culverts
    step_start_time = time.perf_counter()
    culvert.post_process_potential_culverts()
    logging.info(
        f"   x Post processed potential culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # split hydroobjects at endpoints culverts
    step_start_time = time.perf_counter()
    culvert.splits_hydroobject_by_endpoints_of_culverts_and_combine()
    logging.info(
        f"   x Split hydroobjects at endpoints of culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # check if culverts are in the right direction
    step_start_time = time.perf_counter()
    culvert.check_culverts_direction()
    logging.info(
        f"   x Checked culverts direction in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    step_start_time = time.perf_counter()
    culvert.combine_culvert_with_line()
    logging.info(
        f"   x Combined culverts with lines in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    step_start_time = time.perf_counter()
    culvert.splits_hydroobject_by_endpoind_of_culverts_and_combine_2()
    logging.info(
        f"   x Split hydroobjects at endpoints of culverts in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    step_start_time = time.perf_counter()
    culvert.get_shortest_path_from_overige_watergang_to_hydroobjects()
    logging.info(
        f"   x Get shortest path from overige watergang to hydroobjects in "
        f"{_format_duration(time.perf_counter() - step_start_time)}"
    )
    # create map
    if create_html_map:
        step_start_time = time.perf_counter()
        culvert.generate_folium_map(open_html_map=open_html_map)
        logging.info(
            f"   x Generated map in "
            f"{_format_duration(time.perf_counter() - step_start_time)}"
        )

    logging.info(
        f"   x Case finished in "
        f"{_format_duration(time.perf_counter() - start_time)}"
    )
    return culvert
