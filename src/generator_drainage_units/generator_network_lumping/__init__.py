from pathlib import Path

from .generator_network_lumping import GeneratorNetworkLumping


def run_generator_network_lumping(
    path: Path,
    dir_basis_data: str = "0_basisdata",
    dir_inter_results: str = "1_tussenresultaat",
    dir_results: str = "2_resultaat",
    direction: str = "upstream",
    water_lines: list[str] = None,
    include_areas: bool = True,
    no_inflow_outflow_points: int = None,
    detect_split_points: bool = False,
    read_results: bool = True,
    write_results: bool = False,
    html_file_name: str = None,
    width_edges: float = 10.0,
    opacity_edges: float = 0.5,
):
    network = GeneratorNetworkLumping(
        path=path,
        dir_basis_data=dir_basis_data,
        dir_inter_results=dir_inter_results,
        dir_results=dir_results,
    )
    network.create_graph_from_network(water_lines=water_lines)

    network.find_upstream_downstream_nodes_edges(
        direction=direction,
        no_inflow_outflow_points=no_inflow_outflow_points,
    )
    network.calculate_angles_of_edges_at_nodes()

    if include_areas:
        network.assign_drainage_units_to_outflow_points_based_on_length_hydroobject()
        network.dissolve_assigned_drainage_units()

    if detect_split_points:
        network.detect_split_points()
        network.export_detected_split_points()

    if write_results:
        network.export_results_to_gpkg()
        network.generate_folium_map(
            html_file_name=html_file_name,
            width_edges=width_edges,
            opacity_edges=opacity_edges,
        )
    return network

