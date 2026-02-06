import logging
import time
from pathlib import Path
import geopandas as gpd

from .generator_network_lumping import GeneratorNetworkLumping


def run_generator_network_lumping(
    path: Path,
    dir_basisdata: str = "0_basisdata",
    dir_results: str = "1_resultaat",
    direction: str = "upstream",
    water_lines: list[str] = None,
    afvoergebied_from_results: str = "afvoergebied_1_gdf.gpkg",
    include_areas_based_on_id: bool = False,
    include_areas_based_on_length: bool = False,
    no_inflow_outflow_points: int = None,
    detect_split_points: bool = False,
    smooth_area: bool = False,
    html_include_units: bool = True,
    read_results: bool = True,
    write_results: bool = False,
    html_file_name: str = None,
    width_edges: float = 10.0,
    opacity_edges: float = 0.5,
    create_html_map: bool = False,
):
    start_time = time.time()
    network = GeneratorNetworkLumping(
        path=path,
        dir_basisdata=dir_basisdata,
        dir_results=dir_results,
    )
    network.create_graph_from_network()

    network.find_upstream_downstream_nodes_edges(
        direction=direction,
        no_inflow_outflow_points=no_inflow_outflow_points,
    )
    if network.afwateringseenheden is None:
        path_afvoergebied = Path(network.dir_results, afvoergebied_from_results)
        if path_afvoergebied.exists():
            network.afwateringseenheden = gpd.read_file(path_afvoergebied)

    network.calculate_angles_of_edges_at_nodes()
    network.select_downstream_upstream_edges_angle(min_difference_angle=20.0)

    if include_areas_based_on_id:
        network.assign_afvoergebied_to_outflow_points_based_on_id()
        network.dissolve_assigned_afvoergebied(smooth_area=smooth_area)
    elif include_areas_based_on_length:
        network.assign_afvoergebied_to_outflow_points_based_on_length_hydroobject()
        network.dissolve_assigned_afvoergebied(smooth_area=smooth_area)

    if detect_split_points:
        network.detect_split_points()
        network.export_detected_split_points()

    if write_results:
        network.export_results_to_gpkg_or_nc(
            list_layers=[
                "inflow_outflow_points",
                "inflow_outflow_edges",
                "inflow_outflow_nodes_hydro",
                "inflow_outflow_areas",
            ]
        )
    if create_html_map:
        network.generate_folium_map(
            html_file_name=html_file_name,
            width_edges=width_edges,
            opacity_edges=opacity_edges,
            html_include_units=html_include_units,
            include_areas=include_areas_based_on_id or include_areas_based_on_length,
        )

    logging.info(f"   x Case finished in {round(time.time()-start_time, 3)} seconds")
    return network
