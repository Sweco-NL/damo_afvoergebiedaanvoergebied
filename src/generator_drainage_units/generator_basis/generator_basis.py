import logging
import warnings
from pathlib import Path

import folium
import geopandas as gpd
import rioxarray as rio
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from shapely.geometry import LineString, Point


class GeneratorBasis(BaseModel):
    """Basis class for all generators"""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None

    dir_basisdata: str | Path = "0_basisdata"
    dir_inter_results: str | Path | None = "1_tussenresultaat"
    dir_results: str | Path | None = "2_resultaat"

    read_results: bool = False
    write_results: bool = False

    folium_map: folium.Map = None

    results: list = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.check_case_path_directory(path=self.path)
            self.read_data_from_case()
            self.use_processed_hydroobjecten()

    def check_case_path_directory(self, path: Path):
        """Checks if case directory exists and if required directory structure exists

        Parameters
        ----------
        path : Path
            path to case directory. name of directory is used as case name.
            self.path and self.name are set

        Raises ValueErrors in case directory and 0_basisdata directory not exist
        """
        if not path.exists() and path.is_dir():
            raise ValueError(
                f"provided path [{path}] does not exist or is not a directory"
            )
        self.path = path
        self.name = self.path.name
        logging.info(f' ### Case "{self.name.capitalize()}" ###')

        # check if directories 0_basisdata and 1_tussenresultaat and 2_resultaat exist
        if isinstance(self.dir_basisdata, str):
            self.dir_basisdata = Path(self.path, self.dir_basisdata)
        if not isinstance(self.dir_basisdata, Path) or not self.dir_basisdata.exists():
            raise ValueError(
                f"provided [{self.dir_basisdata}] is not a path or does not exist"
            )

        if self.dir_inter_results is not None:
            if isinstance(self.dir_inter_results, str):
                self.dir_inter_results = Path(self.path, self.dir_inter_results)
            if isinstance(self.dir_inter_results, Path):
                if not self.dir_inter_results.exists():
                    self.dir_inter_results.mkdir(parents=True, exist_ok=True)

        if self.dir_results is not None:
            if isinstance(self.dir_results, str):
                self.dir_results = Path(self.path, self.dir_results)
            if isinstance(self.dir_results, Path):
                if not self.dir_results.exists():
                    self.dir_results.mkdir(parents=True, exist_ok=True)

        logging.debug(f"    - dir basisdata    = {self.dir_basisdata}")
        logging.debug(f"    - dir interresults = {self.dir_inter_results}")
        logging.debug(f"    - dir results      = {self.dir_results}")

    def read_data_from_case(self, path: Path = None, read_results: bool = None):
        """Read data from case: including basis data and intermediate results

        Parameters
        ----------
        path : Path, optional
            Path to the case directory including directories 0_basisdata and
            1_tussenresultaat. Directory name is used as name for the case,
            by default None
        read_results : bool, optional
            if True, it reads already all resulst from, by default None
        """
        if path is not None and path.exists():
            self.check_case_path_directory(path=path)

        def read_attributes_from_folder(path_dir: Path):
            for f in path_dir.glob("**/*"):
                if hasattr(self, f.stem):
                    logging.debug(f"    - get dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        setattr(self, f.stem, rio.open_dataset(f))

        logging.info(f"   x read basisdata")
        if self.dir_basisdata is not None and self.dir_basisdata.exists():
            read_attributes_from_folder(self.dir_basisdata)

        if self.read_results:
            logging.info(f"   x read results")
            if self.dir_inter_results is not None and self.dir_inter_results.exists():
                read_attributes_from_folder(self.dir_inter_results)
            if self.dir_results is not None and self.dir_results.exists():
                read_attributes_from_folder(self.dir_results)

    def use_processed_hydroobjecten(self, processed_file="processed"):
        for watergang in ["hydroobjecten", "overige_watergangen"]:
            if getattr(self, watergang, None) is None:
                logging.debug(f"    - attribute {watergang} does not exist")
                continue

            watergang_processed_file = None
            for dir_results in [self.dir_inter_results, self.dir_results]:
                if dir_results is None:
                    continue
                files_in_dir = [f for f in dir_results.glob("**/*")]
                for f in files_in_dir:
                    if f.stem == f"{watergang}_{processed_file}":
                        watergang_processed_file = f

            if watergang_processed_file is not None:
                logging.debug(f"    - get dataset processed {watergang}")
                setattr(self, watergang, gpd.read_file(watergang_processed_file))
