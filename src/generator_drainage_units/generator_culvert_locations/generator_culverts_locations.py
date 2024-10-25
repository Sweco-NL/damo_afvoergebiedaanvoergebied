from pydantic import BaseModel, ConfigDict
import pandas as pd
import geopandas as gpd
from pathlib import Path
import logging
from shapely.geometry import Point, LineString, MultiLineString, Polygon

from ..utils.general_functions import line_to_vertices


class GeneratorCulvertLocations(BaseModel):
    """"Module to guess (best-guess) the locations of culverts 
    based on existing water network, other water bodies (c-watergangen), 
    roads and level areas (peilgebieden)."""
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None
    hydroobjecten: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None
    bebouwing: gpd.GeoDataFrame = None
    hydroobjecten: gpd.GeoDataFrame = None
    keringen: gpd.GeoDataFrame = None
    nwb: gpd.GeoDataFrame = None
    peilgebieden: gpd.GeoDataFrame = None
    snelwegen: gpd.GeoDataFrame = None
    spoorwegen: gpd.GeoDataFrame = None

    water_line_pnts: gpd.GeoDataFrame = None
    potential_culverts: gpd.GeoDataFrame = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.check_case_path_directory(path=self.path)
            self.read_basisdata_from_case()

    def check_case_path_directory(self, path: Path):
        if not path.exists():
            raise ValueError(f"provided path [{path}] does not exist")
        self.path = path
        self.name = self.path.name
        logging.info(f"Case: {self.name.capitalize()} initiated")
        # check if directories 0_basisdata and 1_tussenresultaat exist
        if not Path(self.path, "0_basisdata").exists():
            raise ValueError(f"provided path exists but does not include a 0_basisdata")
        for folder in ["1_tussenresultaat"]:
            if not Path(self.path, folder).exists():
                Path(self.path, folder).mkdir(parents=True, exist_ok=True)

    def read_basisdata_from_case(self, path: Path = None):
        """Read basisdata from case: 
        using: get data from subfolder 0_basisdata within path of case"""
        if path is not None and path.exists():
            self.check_case_path_directory(path=path)
        logging.info(f" x read basisdata")
        list_input_gpkgs = Path(self.path, "0_basisdata").glob("**/*")
        for x in list_input_gpkgs:
            if x.is_file():
                if hasattr(self, x.stem):
                    logging.debug(f" - get dataset {x.stem}")
                    setattr(self, x.stem, gpd.read_file(x, layer=x.stem))

    def generate_vertices_along_waterlines(
        self, 
        distance_vertices=10,
        waterlines=["hydroobjecten", "overige_watergangen"],
        write_results=False
    ):
        waterlines = pd.concat([getattr(self, l) for l in waterlines])
        logging.info(f" x generate vertices for {len(waterlines)} waterlines")
        self.water_line_pnts = line_to_vertices(waterlines, distance=distance_vertices)
        if write_results:
            dir_results = Path(self.path, "1_tussenresultaat")
            self.water_line_pnts.to_file(Path(dir_results, "water_line_pnts.gpkg"), layer="water_line_pnts")
        return self.water_line_pnts
    
