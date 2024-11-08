import logging
from pathlib import Path

from pydantic import BaseModel, ConfigDict
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

from ..utils.general_functions import shorten_line_two_vertices, line_to_vertices


class GeneratorOrderLevels(BaseModel):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None

    hydroobjecten: gpd.GeoDataFrame = None
    rws_wateren: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    # water_line_pnts: gpd.GeoDataFrame = None
    results_0: gpd.GeoDataFrame = None
    results_1: gpd.GeoDataFrame = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.check_case_path_directory(path=self.path)
            self.read_data_from_case()

    def check_case_path_directory(self, path: Path):
        """Checks if case directory exists and if required directory structure exists

        Parameters
        ----------
        path : Path
            path to case directory. name of directory is used as case name.
            self.path and self.name are set

        Raises ValueErrors in case directory and 0_basisdata directory do not exist
        """
        if not path.exists() and path.is_dir():
            raise ValueError(
                f"provided path [{path}] does not exist or is not a directory"
            )
        self.path = path
        self.name = self.path.name
        logging.info(f' ### Case "{self.name.capitalize()}" ###')
        # check if directories 0_basisdata and 1_tussenresultaat exist
        if not Path(self.path, "0_basisdata").exists():
            raise ValueError(f"provided path [{path}] exists but without a 0_basisdata")
        for folder in ["1_tussenresultaat", "2_resultaat"]:
            if not Path(self.path, folder).exists():
                Path(self.path, folder).mkdir(parents=True, exist_ok=True)

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
        logging.info(f"   x read basisdata")
        basisdata_gpkgs = [
            Path(self.path, "0_basisdata", f + ".gpkg")
            for f in [
                "hydroobjecten",
                "rws_wateren",
                "overige_watergangen"
            ]
        ]
        if isinstance(read_results, bool):
            self.read_results = read_results
        baseresults_gpkgs = (
            [
                Path(self.path, "1_tussenresultaat", f + ".gpkg")
                for f in [
                    "water_line_pnts",
                    "results_0",
                    "results_1",
                ]
            ]
            if self.read_results
            else []
        )
        for list_gpkgs in [basisdata_gpkgs, baseresults_gpkgs]:
            for x in list_gpkgs:
                if x.is_file():
                    if hasattr(self, x.stem):
                        logging.debug(f"    - get dataset {x.stem}")
                        setattr(self, x.stem, gpd.read_file(x, layer=x.stem))

