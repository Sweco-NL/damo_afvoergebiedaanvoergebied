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
    
    def find_potential_culvert_locations(self, water_line_pnts=None, max_culvert_length=40, write_results=False):
        if water_line_pnts is None:
            if self.water_line_pnts is None:
                self.generate_water_line_pnts_along_waterlines()
            water_line_pnts = self.water_line_pnts.copy()

        logging.info(f" x find potential culvert locations")
        water_line_pnts["unique_id"] = water_line_pnts.index

        # Filter for end points
        end_pnts = water_line_pnts[water_line_pnts['line_type'] == 'dangling']
        end_pnts = end_pnts.rename(columns={"CODE": "dangling_CODE", "unique_id": "dangling_id"})
        end_pnts_orig_geometry = end_pnts.geometry

        # end points have buffer as geometry, perform spatial join with points
        logging.debug(f" - spatial join vertices")
        end_pnts["geometry"] = end_pnts.buffer(max_culvert_length)
        end_pnts = gpd.sjoin(
            end_pnts, 
            water_line_pnts[["geometry", "unique_id", "CODE"]], 
            how='left', 
            predicate='intersects'
        ).drop(columns="index_right")

        # remove connections with same hydroobject
        end_pnts = end_pnts[end_pnts['dangling_CODE'] != end_pnts['CODE']]

        # make water_line_pnts id into unique_id again, 
        # restore dangling geometry to points, 
        # merge with original water_line_pnts, 
        logging.debug(f" - add data to potential culverts")
        end_pnts["geometry"] = end_pnts_orig_geometry
        end_pnts = end_pnts.rename(columns={"geometry": "geometry2"})

        end_pnts = pd.merge(
            water_line_pnts,
            end_pnts[["unique_id", "dangling_id", "dangling_CODE", 'geometry2']], 
            on='unique_id', 
            how='inner'
        )

        #create lines
        potential_culverts = end_pnts.copy()
        potential_culverts.dropna(subset=['dangling_id'], inplace=True)
        potential_culverts["geometry"] = potential_culverts.apply(
            lambda x: LineString([x["geometry"], x["geometry2"]]), 
            axis=1
        )
        self.potential_culverts = potential_culverts.drop(columns="geometry2")

        if write_results:
            dir_results = Path(self.path, "1_tussenresultaat")
            self.potential_culverts.to_file(Path(dir_results, "potential_culverts.gpkg"), layer="potential_culverts")

        logging.debug(f" - {len(self.potential_culverts)} potential culverts generated")
        return self.potential_culverts
        

