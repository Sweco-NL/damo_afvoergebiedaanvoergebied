from pathlib import Path
import logging
from dotenv import dotenv_values
import time

from generator_drainage_units import run_generator_culvert_locations


logging.basicConfig(level=logging.DEBUG)

config = dotenv_values(".env")
base_dir = config["BASE_DIR"]
# case_name = "vallei_en_veluwe"
case_name = "geerestein"
# case_name = "hattemerbroek"
# case_name = "pangelerbeek"

case_path = Path(base_dir, case_name)

cg = run_generator_culvert_locations(
    path=case_path,
    distance_vertices=10,
    max_culvert_length=40,
    read_results=True,
    write_results=True,
)

logging.info("Analyse afgerond")
