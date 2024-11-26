from pathlib import Path

from .generator_order_levels import GeneratorOrderLevels


def run_generator_order_levels(
    path: Path,
    waterschap: str,
    range_orde_code_min: int,
    range_orde_code_max: int,
    read_results: bool = True,
    write_results: bool = True,
    create_html_map: bool = False,
) -> GeneratorOrderLevels:
    order = GeneratorOrderLevels(
        waterschap=waterschap,
        range_orde_code_min=range_orde_code_min,
        range_orde_code_max=range_orde_code_max,
        read_results=read_results,
        write_results=write_results,
    )
    order.read_data_from_case(path=path)

    order.find_end_points_hydroobjects()
    order.generate_rws_code_for_all_outflow_points()

    if create_html_map:
        order.generate_folium_map()
    return order
