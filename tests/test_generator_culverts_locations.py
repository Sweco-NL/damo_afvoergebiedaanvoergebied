from pathlib import Path
import logging
from dotenv import dotenv_values

from generator_drainage_units import run_generator_culvert_locations

logging.basicConfig(level=logging.INFO)

config = dotenv_values(".env")
base_dir = config["BASE_DIR"]
# case_names = ["vallei_en_veluwe"]
case_names = ["geerestein", "hattemerbroek", "pangelerbeek"]
    
for case_name in case_names:
    case_path = Path(base_dir, case_name)

    run_generator_culvert_locations(
        path=case_path,
        distance_vertices=10,
        max_culvert_length=40,
        write_results=True
    )

logging.info('ready')
