import logging
from pathlib import Path

from pydantic import BaseModel, ConfigDict
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

from ..utils.general_functions import line_to_vertices


class GeneratorCulvertLocations(BaseModel):
    """ "Module to guess (best-guess) the locations of culverts
    based on existing water network, other water bodies (c-watergangen),
    roads and level areas (peilgebieden)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None

    hydroobjecten: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None
    bebouwing: gpd.GeoDataFrame = None
    keringen: gpd.GeoDataFrame = None
    nwb: gpd.GeoDataFrame = None
    peilgebieden: gpd.GeoDataFrame = None
    snelwegen: gpd.GeoDataFrame = None
    spoorwegen: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    water_line_pnts: gpd.GeoDataFrame = None
    potential_culverts: gpd.GeoDataFrame = None

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

        Raises ValueErrors in case directory and 0_basisdata directory not exist
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
        for folder in ["1_tussenresultaat"]:
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
                "overige_watergangen",
                "bebouwing",
                "keringen",
                "nwb",
                "peilgebieden",
                "snelwegen",
                "spoorwegen",
            ]
        ]
        if isinstance(read_results, bool):
            self.read_results = read_results
        baseresults_gpkgs = (
            [
                Path(self.path, "1_tussenresultaat", f + ".gpkg")
                for f in ["water_line_pnts", "potential_culverts"]
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


    def generate_vertices_along_waterlines(
        self,
        distance_vertices: float = 10.0,
        waterlines: list[str] = ["hydroobjecten", "overige_watergangen"],
        read_results: bool = None,
        write_results: bool = None,
    ) -> gpd.GeoDataFrame:
        """Generate vertices (water_line_pnts) along waterlines every x meters

        Parameters
        ----------
        distance_vertices : float, optional
            distacne between vertices, by default 10.0
        waterlines : list[str], optional
            list of attributes to be included, by default ["hydroobjecten", "overige_watergangen"]
        read_results : bool, optional
            option (True/False) to read previous results from gpkg, by default None
        write_results : bool, optional
            option (True/False) to write results to case folder in gpkg, by default None

        Returns
        -------
        vertices (water_line_pnts): gpd.GeoDataFrame
            All points within a geodataframe
        """
        if isinstance(read_results, bool):
            self.read_results = read_results
        if isinstance(write_results, bool):
            self.write_results = write_results

        waterlines = pd.concat([getattr(self, l) for l in waterlines])
        if self.read_results and self.water_line_pnts is not None:
            logging.info(
                f"   x {len(self.water_line_pnts)} vertices for {len(waterlines)} waterlines already generated"
            )
            return self.water_line_pnts

        logging.info(f"   x generate vertices for {len(waterlines)} waterlines")
        self.water_line_pnts = line_to_vertices(waterlines, distance=distance_vertices)
        if self.write_results:
            dir_results = Path(self.path, "1_tussenresultaat")
            self.water_line_pnts.to_file(
                Path(dir_results, "water_line_pnts.gpkg"), layer="water_line_pnts"
            )
        return self.water_line_pnts

    def find_potential_culvert_locations(
        self,
        water_line_pnts=None,
        max_culvert_length=40,
        read_results=None,
        write_results=None,
    ) -> gpd.GeoDataFrame:
        """Find potential culvert locations based on water_line_pnts.
        THe connections between two points from different waterlines
        with a maximum distance of x m (max_culvert_length)

        Parameters
        ----------
        water_line_pnts : _type_, optional
            points every x m along the waterlines; includes a water_line_id, by default None
        max_culvert_length : int, optional
            maximum culvert length: in case of larger distance between points, connections are not made, by default 40
        read_results : bool, optional
            option (True/False) to read previous results from gpkg, by default None
        write_results : bool, optional
            option (True/False) to write results to case folder in gpkg, by default None

        Returns
        -------
        potential_culverts: gpd.GeoDataFrame
            Locations potential culverts (between different water_lines)
        """
        # check read_results and write_results
        if isinstance(read_results, bool):
            self.read_results = read_results
        if isinstance(write_results, bool):
            self.write_results = write_results

        # check if water_line_pnts are given. if not known, use function to generate.
        if water_line_pnts is None:
            if self.water_line_pnts is None:
                self.generate_water_line_pnts_along_waterlines()
        else:
            self.water_line_pnts = water_line_pnts

        # check if potential culverts already exists and should be read
        if self.read_results and self.potential_culverts is not None:
            logging.info(
                f"   x {len(self.potential_culverts)} potential culverts already generated"
            )
            return self.potential_culverts

        logging.info("   x find potential culvert locations")
        self.water_line_pnts["unique_id"] = self.water_line_pnts.index

        # Filter for end points
        end_pnts = self.water_line_pnts[self.water_line_pnts["line_type"] == "dangling"]
        end_pnts = end_pnts.rename(
            columns={"CODE": "dangling_CODE", "unique_id": "dangling_id"}
        )
        end_pnts_orig_geometry = end_pnts.geometry

        # end points have buffer as geometry, perform spatial join with points
        logging.debug("    - spatial join vertices")
        end_pnts["geometry"] = end_pnts.buffer(max_culvert_length)
        end_pnts = gpd.sjoin(
            end_pnts,
            self.water_line_pnts[["geometry", "unique_id", "CODE"]],
            how="left",
            predicate="intersects",
        ).drop(columns="index_right")

        # remove connections with same hydroobject
        end_pnts = end_pnts[end_pnts["dangling_CODE"] != end_pnts["CODE"]]

        # make water_line_pnts id into unique_id again,
        # restore dangling geometry to points,
        # merge with original water_line_pnts,
        logging.debug("    - add data to potential culverts")
        end_pnts["geometry"] = end_pnts_orig_geometry
        end_pnts = end_pnts.rename(columns={"geometry": "geometry2"})

        end_pnts = pd.merge(
            self.water_line_pnts,
            end_pnts[["unique_id", "dangling_id", "dangling_CODE", "geometry2"]],
            on="unique_id",
            how="inner",
        )

        # create potential culverts (linestrings)
        potential_culverts = end_pnts.copy()
        potential_culverts.dropna(subset=["dangling_id"], inplace=True)
        potential_culverts["geometry"] = potential_culverts.apply(
            lambda x: LineString([x["geometry"], x["geometry2"]]), axis=1
        )
        self.potential_culverts = potential_culverts.drop(columns="geometry2")

        if write_results:
            dir_results = Path(self.path, "1_tussenresultaat")
            self.potential_culverts.to_file(
                Path(dir_results, "potential_culverts.gpkg"), layer="potential_culverts"
            )

        logging.debug(
            f"    - {len(self.potential_culverts)} potential culverts generated"
        )
        return self.potential_culverts
