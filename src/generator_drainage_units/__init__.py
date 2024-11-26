import logging
import time
from pathlib import Path

from .generator_network_lumping import GeneratorNetworkLumping
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

    culverts_generator.select_correct_score_based_on_score_and_length()

    culverts_generator.post_process_potential_culverts()

    culverts_generator.splits_hydroobjecten_by_endpoints_of_culverts_and_combine()

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return culverts_generator



def run_generator_network_lumping(
    path: Path,
    direction: str = "upstream",
    water_lines: list[str] = None,
    include_areas: bool = True,
    no_inflow_outflow_points: int = None,
    detect_split_points: bool = False,
    write_results: bool = False,
    html_file_name: str = None,
    width_edges: float = 10.0,
    opacity_edges: float = 0.5,
):
    network = GeneratorNetworkLumping()
    network.read_data_from_case(path=path)
    network.create_graph_from_network(water_lines=water_lines)

    network.find_upstream_downstream_nodes_edges(
        direction=direction,
        no_inflow_outflow_points=no_inflow_outflow_points,
    )
    network.calculate_angles_of_edges_at_splitpoints()
    network.select_directions_for_splits_based_on_angle()

    # Find upstream nodes again?
    #network.find_upstream_downstream_nodes_edges(direction=network.direction)
    network.assign_drainage_units_to_outflow_points_based_on_length_hydroobject()
    network.dissolve_assigned_drainage_units()

    if detect_split_points:
        network.detect_split_points()
        network.export_detected_split_points()

    if write_results:
        network.export_results_to_gpkg()
        network.export_results_to_html_file(
            html_file_name=html_file_name,
            width_edges=width_edges,
            opacity_edges=opacity_edges,
        )
    return network


def run_network_lumping_with_random_selection_splits(
    network: GeneratorNetworkLumping,
    include_areas: bool = True,
    write_html: bool = False
):
    network.select_directions_for_splits()
    network.find_upstream_downstream_nodes_edges(direction=network.direction)
    network.assign_drainage_units_to_outflow_points_based_on_length_hydroobject()
    network.dissolve_assigned_drainage_units()
    network.export_results_to_html_file(html_file_name=f"{network.name}_random_selection_splits")
    return network
