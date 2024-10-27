from pathlib import Path
from .generator_culvert_locations import GeneratorCulvertLocations


def run_generator_culvert_locations(
    path: Path,
    distance_vertices: float = 10.0,
    max_culvert_length: float = 40.0,
    read_results: bool = False,
    write_results: bool = False,
):
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
    return culverts_generator
