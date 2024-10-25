from pathlib import Path
from .generator_culvert_locations import GeneratorCulvertLocations


def run_generator_culvert_locations(
    path: Path, 
):
    culverts_generator = GeneratorCulvertLocations()

    # read basis data from folder 0_basisdata
    culverts_generator.read_basisdata_from_case(path=path)

