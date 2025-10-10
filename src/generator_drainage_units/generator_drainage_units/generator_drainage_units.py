import logging
import time
from pathlib import Path

import folium
import geopandas as gpd
import imod
import numpy as np
import pandas as pd
import pyflwdir
import rioxarray
import xarray
from geocube.api.core import make_geocube
from pydantic import ConfigDict
from rasterio.enums import Resampling
from scipy.ndimage import distance_transform_edt
from tqdm import tqdm

from ..generator_basis import GeneratorBasis
import pyflwdir
from ..utils.pyflwdir import run_pyflwdir
from ..utils.folium_map import generate_folium_map


def change_flow_direction_d8_to_d16(dir_d8):
    """Change flow direction from D8 to D16."""
    dir_d16 = dir_d8.copy()
    for i, j in zip([1, 2, 4, 8, 16, 32, 64, 128], [7, 9, 11, 13, 15, 1, 3, 5]):
        dir_d16.data[dir_d8 == i] = j
    return dir_d16


def get_resolution_2d_array(dataarray, x: str = 'x', y: str = 'y', decimals=2):
    # Access the coordinates, e.g., 'x' and 'y' for a 2D grid
    x_coords = dataarray.coords[x].values
    y_coords = dataarray.coords[y].values

    # Calculate the resolution
    x_resolution = round(abs(x_coords[1] - x_coords[0]), decimals)
    y_resolution = round(abs(y_coords[1] - y_coords[0]), decimals)
    return x_resolution, y_resolution


class GeneratorDrainageUnits(GeneratorBasis):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str = "1_resultaat"
    waterschap: str = None

    method: str = "pyflwdir"

    read_results: bool = False
    write_results: bool = False
    crs: int = 28992
    
    required_results: list[str] = [
        "hydroobjecten",
        "hydroobjecten_processed_0", 
        "overige_watergangen", 
        "outflow_nodes_overige_watergangen",
        "overige_watergangen_processed_3", 
        "overige_watergangen_processed_4", 
        "potential_culverts_5", 
        "outflow_nodes",
        "edges", 
        "nodes",
        "all_waterways_0",
        "drainage_units_0",
    ]

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    
    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed_3: gpd.GeoDataFrame = None
    overige_watergangen_processed_4: gpd.GeoDataFrame = None
    outflow_nodes_overige_watergangen: gpd.GeoDataFrame = None
    potential_culverts_5: gpd.GeoDataFrame = None

    outflow_nodes: gpd.GeoDataFrame = None

    ghg_file_name: str = None
    ghg: xarray.Dataset = None
    ghg_filled: xarray.Dataset = None
    ghg_fills: xarray.Dataset = None

    new_resolution: float = 5.0*5.0

    all_waterways_0: gpd.GeoDataFrame = None
    all_waterways_0_buffer: gpd.GeoDataFrame = None
    all_waterways_1: gpd.GeoDataFrame = None

    ghg_waterways: xarray.Dataset = None
    ghg_waterways_distance: xarray.Dataset = None

    ghg_processed: xarray.Dataset = None
    ghg_processed_adapt: xarray.Dataset = None

    flw: pyflwdir.FlwdirRaster = None

    flow_direction: xarray.Dataset = None
    flow_direction_d8: xarray.Dataset = None
    flow_direction_d8_ind: xarray.Dataset = None
    flow_direction_d8_fills: xarray.Dataset = None
    flow_direction_d16: xarray.Dataset = None
    flow_direction_d16_ind: xarray.Dataset = None
    flow_direction_d16_fills: xarray.Dataset = None
    
    drainage_units_clean: gpd.GeoDataFrame = None
    
    drainage_units_0: xarray.Dataset = None
    drainage_units_0_gdf: gpd.GeoDataFrame = None
    drainage_units_1: xarray.Dataset = None
    drainage_units_1_gdf: gpd.GeoDataFrame = None
    drainage_units_2: xarray.Dataset = None
    drainage_units_2_gdf: gpd.GeoDataFrame = None
    drainage_units_3: xarray.Dataset = None
    drainage_units_3_gdf: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    
    folium_map: folium.Map = None
    folium_html_path: str = None


    def read_ghg(self, ghg_file_name: str):
        logging.info("   x read topographical data as input")
        self.ghg_file_name = ghg_file_name
        self.ghg = rioxarray.open_rasterio(Path(self.path, self.dir_basisdata, ghg_file_name))
        self.ghg.name = "ghg"
        if self.ghg.rio.crs is None:
            self.ghg = self.ghg.rio.write_crs(28992)
        return self.ghg
    

    def preprocess_ghg(self, resolution=2.0, depth_waterways=1.0, buffer_waterways=None, smooth_distance=25.0):
        logging.info("   x preprocessing GHG data")
        logging.info("     - resampling data to new resolution")

        self.new_resolution = resolution
        if buffer_waterways is None:
            buffer_waterways = 1.5 * resolution

        old_resolution = 25.0
        upscale_factor = old_resolution / resolution
        new_width = int(np.ceil(self.ghg.rio.width * upscale_factor))
        new_height = int(np.ceil(self.ghg.rio.height * upscale_factor))

        self.ghg_processed = self.ghg.rio.reproject(
            self.ghg.rio.crs,
            shape=(new_height, new_width),
            resampling=Resampling.bilinear,
        )
        self.ghg_processed.name = "ghg_processed"

        logging.info("     - fill holes")
        self.ghg_filled = self.ghg_processed.copy()
        self.flow_direction_d8 = self.ghg_processed.copy()
        self.flow_direction_d8.name = "flow_direction_d8"
        self.flow_direction_d8 = self.flow_direction_d8.astype(np.int32)

        # use pyflwdir to fill depressions in the GHG data and get D8-direction
        self.ghg_filled.data[0], self.flow_direction_d8.data[0] = pyflwdir.fill_depressions(
            self.ghg_filled.data[0],
            nodata=self.ghg_filled._FillValue,
        )
        # find out where holes are filled
        self.ghg_fills = self.ghg_filled - self.ghg
        # get D8-direction and corresponding D16-direction
        self.flow_direction_d8_fills = self.flow_direction_d8.where(self.ghg_fills>0.001).fillna(-1).astype(np.int32)
        self.flow_direction_d16_fills = change_flow_direction_d8_to_d16(self.flow_direction_d8_fills)

        logging.info("     - select waterways and add depth at waterways")
        # combine all waterways and filter on order_no
        edges = self.edges[["code", "order_no", "geometry"]].reset_index(drop=True)
        # select only with order_no
        edges = edges[edges["order_no"]>0]
        self.all_waterways_0 = edges[["code", "geometry"]].reset_index(drop=True)
        # do the same for the other waterways and combine
        if self.overige_watergangen_processed_4 is not None:
            self.all_waterways_0 = pd.concat([
                self.all_waterways_0,
                self.overige_watergangen_processed_4[["code", "geometry"]]
            ]).reset_index(drop=True)
        # add depth
        self.all_waterways_0["depth_waterways"] = depth_waterways

        logging.info("     - give each waterway an unique id")
        self.all_waterways_0["drainage_unit_id"] = self.all_waterways_0.index
        self.all_waterways_0["color_id"] = np.random.shuffle(np.arange(len(self.all_waterways_0)))
        self.all_waterways_0_buffer = self.all_waterways_0.copy()
        self.all_waterways_0_buffer.geometry = self.all_waterways_0.geometry.buffer(buffer_waterways, cap_style="flat")
        
        ghg_waterways = make_geocube(
            vector_data=self.all_waterways_0_buffer,
            measurements=["depth_waterways"],
            like=self.ghg_processed,
        )["depth_waterways"].fillna(0.0)

        logging.info("     - calculate distance to waterways to add depth")
        ghg_waterways_distance = distance_transform_edt(
            ghg_waterways==0.0, 
            sampling=2.0,
        )
        ghg_waterways_distance = xarray.DataArray(
            ghg_waterways_distance, 
            dims=ghg_waterways.dims, 
            coords=ghg_waterways.coords
        )
        ghg_waterways = depth_waterways * 0.5**(ghg_waterways_distance / smooth_distance)
        ghg_waterways.data[ghg_waterways.data<0.001] = 0.0

        self.ghg_processed.data = self.ghg_processed.data - ghg_waterways.data
        self.ghg_processed.data[self.ghg_processed.data<-900.0] = -999.99
        
        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "all_waterways_0", 
                    "all_waterways_0_buffer", 
                    "ghg_processed"
                ]
            )
        return self.ghg_processed


    def generate_drainage_units(self, iterations=2000, iteration_group=100, flow_method="d8"):
        logging.info("   x generate drainage units for each waterway")
        # create raster with unique id of each waterway
        logging.info("     - give each waterway an unique id")
        if self.all_waterways_0_buffer is None or self.ghg_processed is None:
            raise ValueError("   x run preprocess_ghg to preprocess the data")

        self.drainage_units_0 = make_geocube(
            vector_data=self.all_waterways_0_buffer,
            measurements=["drainage_unit_id"],
            like=self.ghg_processed,
        )["drainage_unit_id"].fillna(-1)
        self.drainage_units_0.name = "drainage_units_0"
        self.drainage_units_0.attrs["_FillValue"] = -1
        
        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "drainage_units_0", 
                ]
            )

        # PYFLWDIR
        if self.method == "pyflwdir":
            logging.info("   x run pyflwdir")
            self.drainage_units_0, self.ghg_processed_adapt, flow_direction = run_pyflwdir(
                dem=self.ghg_processed,
                waterways=self.drainage_units_0,
                iterations=iterations,
                iteration_group=iteration_group,
                flow_method=flow_method,
                flow_direction_d16_fills=self.flow_direction_d16_fills
            )
            if flow_method == "d8":
                # self.flow_direction_d8 = self.ghg_processed.copy()
                self.flow_direction_d8.data = flow_direction.reshape(
                    self.ghg_processed.data.shape
                )
                self.flow_direction_d8.name = "flow_direction_d8"
                self.flow_direction_d8.attrs["_FillValue"] = -1
            elif flow_method == "d16":
                self.flow_direction_d16 = flow_direction.copy()
                self.flow_direction_d16.name = "flow_direction_d16"
        elif self.method == "pcraster":
            raise ValueError("method pcraster not yet installed")
        else:
            raise ValueError("method wrong")

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "flow_direction_d16",
                    "drainage_units_0"
                ]
            )
        return self.drainage_units_0


    def aggregate_drainage_units(self):
        self.drainage_units_0_gdf = None
        self.drainage_units_1 = None
        self.drainage_units_1_gdf = None
        self.drainage_units_2 = None
        self.drainage_units_2_gdf = None

        logging.info(f"     - polygonize rasters to polygons")
        if self.drainage_units_0.dims == ("y", "x"):
            gdf = imod.prepare.polygonize(self.drainage_units_0)
        else:
            gdf = imod.prepare.polygonize(self.drainage_units_0[0])

        gdf = gdf.rename(columns={"value": "drainage_unit_id"})
        gdf["drainage_unit_id"] = gdf["drainage_unit_id"].astype(int)
        gdf = gdf.dissolve(by="drainage_unit_id", aggfunc="first").reset_index()
        gdf = gdf.set_crs(self.hydroobjecten.crs)
        random_color_id = np.random.randint(0, 25, size=len(gdf))
        gdf["color_id"] = random_color_id

        gdf["drainage_unit_area"] = gdf["geometry"].area
        gdf = gdf.explode().reset_index(drop=True)
        gdf["part_count"] = gdf[["drainage_unit_id"]].groupby("drainage_unit_id").transform("count").reset_index()
        gdf = gdf[gdf["drainage_unit_id"] > -1]
        gdf = gdf.dissolve(by="drainage_unit_id").reset_index()

        self.drainage_units_0_gdf = gdf.copy()
        self.drainage_units_clean = gdf.copy()
        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "drainage_units_clean"
                ]
            )

        logging.info("   x aggregation/lumping of drainage units: aggregate 'overige watergangen'")
        logging.info("     - define new drainage_unit_ids for all 'overige watergangen'")

        self.edges["downstream_edges"] = self.edges["code"]
        all_waterways_1 = self.all_waterways_0.merge(
            pd.concat([
                self.edges[["code", "downstream_edges", "order_code"]], 
                self.overige_watergangen_processed_4[["code", "downstream_edges", "downstream_order_code"]]
            ]),
            how="left",
            on="code"
        )

        all_waterways_1["order_code"] = all_waterways_1["order_code"].fillna(all_waterways_1["downstream_order_code"])
        all_waterways_1 = all_waterways_1[all_waterways_1["downstream_edges"] != ''].copy()
        all_waterways_1 = all_waterways_1[all_waterways_1["order_code"] != ''].copy()
        self.all_waterways_1 = all_waterways_1.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "all_waterways_1"
                ]
            )
        
        logging.info(f"     - aggregate sub drainage units: replace {len(all_waterways_1)} drainage_unit_ids")
        drainage_units_0_gdf = self.drainage_units_0_gdf.merge(
            self.all_waterways_1[
                ["drainage_unit_id", "code", "downstream_edges", "order_code", "downstream_order_code"]
            ],
            how="left",
            on="drainage_unit_id"
        )
        self.drainage_units_0_gdf = drainage_units_0_gdf.copy()

        self.drainage_units_0_gdf["downstream_order_code"] = self.drainage_units_0_gdf["downstream_order_code"].fillna("")
        self.drainage_units_0_gdf = self.drainage_units_0_gdf.sort_values("downstream_order_code", ascending=True)
        self.drainage_units_1_gdf = (
            self.drainage_units_0_gdf
            .dissolve(by="downstream_edges")
            .drop(columns="downstream_order_code")
            .reset_index()
        )
        random_color_id = np.random.randint(0, 25, size=len(self.drainage_units_1_gdf))
        self.drainage_units_1_gdf["color_id"] = random_color_id

        self.drainage_units_2_gdf = self.drainage_units_1_gdf.dissolve(by="order_code").reset_index()
        random_color_id = np.random.randint(0, 25, size=len(self.drainage_units_2_gdf))
        self.drainage_units_2_gdf["color_id"] = random_color_id
        self.drainage_units_2_gdf["order_code_no"] = self.drainage_units_2_gdf.order_code.str[:6]

        self.drainage_units_3_gdf = self.drainage_units_2_gdf.dissolve("order_code_no").reset_index()
        random_color_id = np.random.randint(0, 25, size=len(self.drainage_units_3_gdf))
        self.drainage_units_3_gdf["color_id"] = random_color_id

        # rasterize gdfs
        def dataarray_from_gdf(raster, gdf, raster_name):
            raster.data = imod.prepare.rasterize(
                gdf.reset_index(drop=True), 
                column="color_id",
                like=raster,
                fill=-1,
            )
            raster.name = raster_name
            return raster
        
        if self.drainage_units_0.dims != ("y", "x"):
            drainage_units_0 = self.drainage_units_0[0]
        else:
            drainage_units_0 = self.drainage_units_0

        self.drainage_units_0 = dataarray_from_gdf(
            drainage_units_0, 
            self.drainage_units_0_gdf, 
            "drainage_units_0"
        )
        self.drainage_units_1 = self.drainage_units_0.copy()
        self.drainage_units_1 = dataarray_from_gdf(
            self.drainage_units_1, 
            self.drainage_units_1_gdf, 
            "drainage_units_1"
        )
        self.drainage_units_2 = self.drainage_units_0.copy()
        self.drainage_units_2 = dataarray_from_gdf(
            self.drainage_units_2, 
            self.drainage_units_2_gdf, 
            "drainage_units_2"
        )

        self.drainage_units_3 = self.drainage_units_0.copy()
        self.drainage_units_3 = dataarray_from_gdf(
            self.drainage_units_3, 
            self.drainage_units_3_gdf, 
            "drainage_units_3"
        )

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "drainage_units_0_gdf",
                    "drainage_units_1_gdf",
                    "drainage_units_2_gdf",
                    "drainage_units_3_gdf",
                    "drainage_units_1",
                    "drainage_units_2",
                    "drainage_units_3",
                ]
            )
        return self.drainage_units_2_gdf


    def generate_folium_map(self, html_file_name:str="", **kwargs):
        if html_file_name == '':
            html_file_name = self.name + "_drainage_units"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            **kwargs
        )
        return self.folium_map
