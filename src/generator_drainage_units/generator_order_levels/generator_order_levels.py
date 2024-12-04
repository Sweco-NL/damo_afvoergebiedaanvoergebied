import logging
from pathlib import Path

from pydantic import BaseModel, ConfigDict
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point
import folium
from folium.features import DivIcon

from ..utils.general_functions import shorten_line_two_vertices, line_to_vertices
from ..utils.folium_utils import add_labels_to_points_lines_polygons


class GeneratorOrderLevels(BaseModel):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None
    waterschap: str = None
    range_orde_code_min: int = None
    range_orde_code_max: int = None

    hydroobjecten: gpd.GeoDataFrame = None
    rws_wateren: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    # water_line_pnts: gpd.GeoDataFrame = None
    dead_end_nodes: gpd.GeoDataFrame = None
    outflow_nodes_all: gpd.GeoDataFrame = None

    folium_map: folium.Map = None


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
            for f in ["hydroobjecten", "rws_wateren", "overige_watergangen"]
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


    def find_end_points_hydroobjects(self, buffer_width=0.5):
        # Copy hydroobject data to new variable 'hydroobjects' and make dataframes with start and end nodes
        hydroobjects = self.hydroobjecten[["CODE", "geometry"]].explode()
        hydroobjects["start_node"] = hydroobjects.geometry.apply(
            lambda x: Point(x.coords[0][0:2])
        )
        hydroobjects["end_node"] = hydroobjects.geometry.apply(
            lambda x: Point(x.coords[-1][0:2])
        )

        start_nodes = hydroobjects[["start_node"]].rename(
            columns={"start_node": "geometry"}
        )
        end_nodes = hydroobjects[["end_node"]].rename(columns={"end_node": "geometry"})

        logging.info(f"   - Start and end nodes defined")

        # Determine dead-end end nodes (using buffer to deal with inaccuracies in geometry hydroobjects) and make new dataframe
        start_nodes["geometry"] = start_nodes.buffer(buffer_width)
        end_nodes["geometry"] = end_nodes.buffer(buffer_width)

        # Determine where start and end nodes overlap
        dead_end_nodes = start_nodes.sjoin(end_nodes, how="right")
        dead_end_nodes = hydroobjects.loc[
            dead_end_nodes[dead_end_nodes.index_left.isna()].index, ["CODE", "end_node"]
        ].copy()
        dead_end_nodes = dead_end_nodes.rename(columns={"end_node": "geometry"})

        dead_end_nodes["geometry"] = dead_end_nodes["geometry"].buffer(buffer_width)
        sjoin = dead_end_nodes.sjoin(hydroobjects[["CODE", "geometry"]])
        no_dead_end_node_ids = sjoin[sjoin["CODE_left"] != sjoin["CODE_right"]].index

        self.dead_end_nodes = dead_end_nodes[
            ~dead_end_nodes.index.isin(no_dead_end_node_ids)
        ].copy()
        self.dead_end_nodes.geometry = self.dead_end_nodes.geometry.centroid

        logging.info(f"   - Dead end nodes defined")

        return self.dead_end_nodes


    def generate_rws_code_for_all_outflow_points(self, buffer_rws=10.0):
        logging.info(f"   - Generating order code for all outflow points")
        rws_wateren = self.rws_wateren.copy()
        rws_wateren.geometry = rws_wateren.geometry.buffer(buffer_rws)
        outflow_nodes = (
            self.dead_end_nodes.sjoin(rws_wateren[["geometry", "rws_code"]])
            .drop(columns="index_right")
            .reset_index(drop=True)
        )

        outflow_nodes = outflow_nodes.rename(columns={"CODE": "hydroobject_code"})
        outflow_nodes_all = None
        for rws_code in outflow_nodes.rws_code.unique():
            outflows = outflow_nodes[outflow_nodes.rws_code == rws_code].reset_index(
                drop=True
            )
            outflows["rws_code_no"] = self.range_orde_code_min + outflows.index
            outflows["orde_code"] = outflows.apply(
                lambda x: f"{x.rws_code}.{str(x.rws_code_no)}", axis=1
            )
            if outflows["rws_code_no"].max() > self.range_orde_code_max:
                logging.info(
                    f" XXX aantal uitstroompunten op RWS-water ({rws_code}) hoger dan range order_code waterschap"
                )
            if outflow_nodes_all is None:
                outflow_nodes_all = outflows
            else:
                outflow_nodes_all = pd.concat([outflow_nodes_all, outflows])
            logging.info(
                f"   - RWS-water ({rws_code}): {len(outflows)} uitstroompunten"
            )

        logging.info(
            f"   - Totaal aantal uitstroompunten op RWS-wateren voor {self.name}: {len(outflow_nodes_all)}"
        )
        self.outflow_nodes_all = outflow_nodes_all.reset_index(drop=True)
        return self.outflow_nodes_all


    def generate_folium_map(self):
        # Make figure
        outflow_nodes_4326 = self.outflow_nodes_all.to_crs(4326)

        m = folium.Map(
            location=[
                outflow_nodes_4326.geometry.y.mean(),
                outflow_nodes_4326.geometry.x.mean(),
            ],
            zoom_start=12,
        )

        folium.GeoJson(
            self.rws_wateren.geometry,
            name="RWS_Wateren",
            z_index=0,
        ).add_to(m)

        folium.GeoJson(
            self.hydroobjecten.geometry,  # .buffer(10),
            name="Watergangen",
            color="blue",
            fill_color="blue",
            zoom_on_click=True,
            z_index=1,
        ).add_to(m)

        folium.GeoJson(
            self.dead_end_nodes,
            name="Outflow points",
            marker=folium.Circle(
                radius=10,
                fill_color="orange",
                fill_opacity=0.4,
                color="orange",
                weight=3,
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=2,
        ).add_to(m)

        fg = folium.FeatureGroup(
            name=f"Uitstroompunten RWS-water", control=True
        ).add_to(m)

        folium.GeoJson(
            self.outflow_nodes_all,
            name="Uitstroompunten RWS-wateren",
            marker=folium.Circle(
                radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=3,
        ).add_to(fg)

        add_labels_to_points_lines_polygons(
            gdf=self.outflow_nodes_all, column="orde_code", label_fontsize=8, fg=fg
        )

        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m
        return m
