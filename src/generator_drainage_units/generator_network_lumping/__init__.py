from pathlib import Path

from .generator_network_lumping import GeneratorNetworkLumping


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
    write_html: bool = False,
):
    network.select_directions_for_splits()
    network.find_upstream_downstream_nodes_edges(direction=network.direction)
    network.assign_drainage_units_to_outflow_points_based_on_length_hydroobject()
    network.dissolve_assigned_drainage_units()
    network.export_results_to_html_file(
        html_file_name=f"{network.name}_random_selection_splits"
    )
    return network
