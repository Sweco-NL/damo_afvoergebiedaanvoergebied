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
import xarray as xr
import networkx as nx
from geocube.api.core import make_geocube
from pydantic import ConfigDict
from rasterio.enums import Resampling
from scipy.ndimage import distance_transform_edt
from tqdm import tqdm

from ..generator_basis import GeneratorBasis
import pyflwdir
from ..utils.pyflwdir import run_pyflwdir
from ..utils.folium_map import generate_folium_map


def change_stroomrichting_d8_to_d16(dir_d8):
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


def dataarray_from_gdf(raster, gdf, raster_name):
    raster.data = imod.prepare.rasterize(
        gdf.reset_index(drop=True), 
        column="color_id",
        like=raster,
        fill=-1,
    )
    raster.name = raster_name
    return raster


class GeneratorAfvoergebieden(GeneratorBasis):
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
        "hydroobject",
        "hydroobject_processed_0", 
        "overige_watergang", 
        "outflow_nodes_overige_watergang",
        "overige_watergang_processed_3", 
        "overige_watergang_processed_4", 
        "outflow_nodes",
        "alle_watergangen_0",
        "afvoergebied_0",
        "afvoergebied_0_gdf",
        "afvoergebied_1",
        "afvoergebied_1_gdf",
        "afvoergebied_2",
        "afvoergebied_2_gdf",
        "edges",
        "nodes",
    ]

    hydroobject: gpd.GeoDataFrame = None
    hydroobject_processed_0: gpd.GeoDataFrame = None
    hydroobject_processed_1: gpd.GeoDataFrame = None
    
    overige_watergang: gpd.GeoDataFrame = None
    overige_watergang_processed_3: gpd.GeoDataFrame = None
    overige_watergang_processed_4: gpd.GeoDataFrame = None
    outflow_nodes_overige_watergang: gpd.GeoDataFrame = None

    outflow_nodes: gpd.GeoDataFrame = None

    ghg_file_name: str = None
    ghg: xr.Dataset = None
    ghg_filled: xr.Dataset = None
    ghg_fills: xr.Dataset = None

    new_resolution: float = 5.0*5.0

    alle_watergangen_0: gpd.GeoDataFrame = None
    alle_watergangen_0_buffer: gpd.GeoDataFrame = None
    alle_watergangen_1: gpd.GeoDataFrame = None

    ghg_waterways: xr.Dataset = None
    ghg_waterways_distance: xr.Dataset = None

    ghg_processed: xr.Dataset = None

    flw: pyflwdir.FlwdirRaster = None

    afvoergebied_0: xr.Dataset = None
    afvoergebied_0_gdf: gpd.GeoDataFrame = None
    afvoergebied_1: xr.Dataset = None
    afvoergebied_1_gdf: gpd.GeoDataFrame = None
    afvoergebied_2: xr.Dataset = None
    afvoergebied_2_gdf: gpd.GeoDataFrame = None
    afvoergebied_3: xr.Dataset = None
    afvoergebied_3_gdf: gpd.GeoDataFrame = None
    afvoergebied_4: xr.Dataset = None
    afvoergebied_4_gdf: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None

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

        logging.info("     - select waterways and add depth at waterways")
        # combine all waterways and filter on order_no
        if "order_no" in self.edges.columns:
            edges = self.edges[["code", "order_no", "geometry"]].reset_index(drop=True)
            edges = edges[edges["order_no"]>0]
        else:
            edges = self.edges[["code", "geometry"]].reset_index(drop=True)

        self.alle_watergangen_0 = edges[["code", "geometry"]].reset_index(drop=True)

        # # do the same for the other waterways and combine
        if self.overige_watergang_processed_4 is not None:
            self.alle_watergangen_0 = pd.concat([
                self.alle_watergangen_0,
                self.overige_watergang_processed_4[["code", "geometry"]]
            ]).reset_index(drop=True)
        # add depth
        self.alle_watergangen_0["depth_waterways"] = depth_waterways

        logging.info("     - give each waterway an unique id")
        self.alle_watergangen_0["drainage_unit_id"] = self.alle_watergangen_0.index
        self.alle_watergangen_0["color_id"] = np.random.shuffle(np.arange(len(self.alle_watergangen_0)))
        self.alle_watergangen_0_buffer = self.alle_watergangen_0.copy()
        self.alle_watergangen_0_buffer.geometry = self.alle_watergangen_0.geometry.buffer(buffer_waterways, cap_style="flat")
        
        ghg_waterways = make_geocube(
            vector_data=self.alle_watergangen_0_buffer,
            measurements=["depth_waterways"],
            like=self.ghg_processed,
        )["depth_waterways"].fillna(0.0)

        logging.info("     - calculate distance to waterways to add depth")
        ghg_waterways_distance = distance_transform_edt(
            ghg_waterways==0.0, 
            sampling=2.0,
        )
        ghg_waterways_distance = xr.DataArray(
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
                    "alle_watergangen_0", 
                    "alle_watergangen_0_buffer", 
                    "ghg_processed"
                ]
            )
        return self.ghg_processed


    def generate_afvoergebied(self):
        logging.info("   x generate drainage units for each waterway")
        # create raster with unique id of each waterway
        logging.info("     - give each waterway an unique id")
        if self.alle_watergangen_0_buffer is None or self.ghg_processed is None:
            raise ValueError("   x run preprocess_ghg to preprocess the data")

        self.afvoergebied_0 = make_geocube(
            vector_data=self.alle_watergangen_0_buffer,
            measurements=["drainage_unit_id"],
            like=self.ghg_processed,
        )["drainage_unit_id"].fillna(-1)
        self.afvoergebied_0.name = "afvoergebied_0"
        self.afvoergebied_0.attrs["_FillValue"] = -1
        
        # PYFLWDIR
        if self.method == "pyflwdir":
            logging.info("     - run pyflwdir to create drainage units (basins)")
            flw_pyflwdir = pyflwdir.from_dem(
                data=self.ghg_processed.data[0],
                nodata=self.ghg_processed._FillValue,
                transform=self.ghg_processed.rio.transform(),
                latlon=False,
            )
            self.flw = flw_pyflwdir
            
            idx = self.afvoergebied_0.astype(np.int32).data.flatten()
            idxs = np.where(idx > 0)[0]
            ids = idx[idx > 0]

            self.afvoergebied_0 = xr.DataArray(
                self.flw.basins(
                    idxs=idxs,
                    ids=ids
                ),
                dims=("y", "x"),
                coords={
                    "y": self.ghg_processed.y, 
                    "x": self.ghg_processed.x
                },
                name="afvoergebied_0",
            )

        elif self.method == "pcraster":
            raise ValueError("method pcraster not yet installed")
        else:
            raise ValueError("method wrong")

        self.afvoergebied_0_gdf = None

        logging.info(f"     - polygonize rasters to polygons")
        if self.afvoergebied_0.dims == ("y", "x"):
            gdf = imod.prepare.polygonize(self.afvoergebied_0)
        else:
            gdf = imod.prepare.polygonize(self.afvoergebied_0[0])

        gdf = gdf.rename(columns={"value": "drainage_unit_id"})
        gdf["drainage_unit_id"] = gdf["drainage_unit_id"].astype(int)
        gdf = gdf.dissolve(by="drainage_unit_id", aggfunc="first").reset_index()
        gdf = gdf.set_crs(self.hydroobject.crs)
        random_color_id = np.random.randint(0, 25, size=len(gdf))
        gdf["color_id"] = random_color_id

        gdf["drainage_unit_area"] = gdf["geometry"].area
        gdf = gdf.explode().reset_index(drop=True)
        gdf["part_count"] = gdf[["drainage_unit_id"]].groupby("drainage_unit_id").transform("count").reset_index()
        gdf = gdf[gdf["drainage_unit_id"] > -1]
        gdf = gdf.dissolve(by="drainage_unit_id").reset_index()

        gdf = gdf.merge(
            self.alle_watergangen_0[["code", "drainage_unit_id"]],
            how="left",
            on="drainage_unit_id",
        )

        self.afvoergebied_0_gdf = gdf.copy()

        # if self.afvoergebied_0.dims != ("y", "x"):
        #     afvoergebied_0 = self.afvoergebied_0[0]
        # else:
        #     afvoergebied_0 = self.afvoergebied_0

        # self.afvoergebied_0 = dataarray_from_gdf(
        #     afvoergebied_0, 
        #     self.afvoergebied_0_gdf, 
        #     "afvoergebied_0"
        # )

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "afvoergebied_0",
                    "afvoergebied_0_gdf",
                ]
            )
        return self.afvoergebied_0


    def aggregate_afvoergebied_tot_2(self):
        self.afvoergebied_1 = None
        self.afvoergebied_1_gdf = None
        self.afvoergebied_2 = None
        self.afvoergebied_2_gdf = None

        logging.info("   x aggregation/lumping of drainage units: aggregate 'overige watergangen'")
        logging.info("     - define new drainage_unit_ids for all 'overige watergangen'")

        self.edges["downstream_edges"] = self.edges["code"]

        alle_watergangen_1 = self.alle_watergangen_0.merge(
            pd.concat([
                self.edges[["code", "downstream_edges", "order_code"]], 
                self.overige_watergang_processed_4[["code", "downstream_edges", "downstream_order_code"]]
            ]),
            how="left",
            on="code"
        )

        alle_watergangen_1["order_code"] = alle_watergangen_1["order_code"].fillna(alle_watergangen_1["downstream_order_code"])
        alle_watergangen_1 = alle_watergangen_1[alle_watergangen_1["downstream_edges"] != ''].copy()
        alle_watergangen_1 = alle_watergangen_1[alle_watergangen_1["order_code"] != ''].copy()
        self.alle_watergangen_1 = alle_watergangen_1.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "alle_watergangen_1"
                ]
            )
        
        logging.info(f"     - aggregate sub drainage units: replace {len(alle_watergangen_1)} drainage_unit_ids")
        self.afvoergebied_1_gdf = self.afvoergebied_0_gdf.merge(
            self.alle_watergangen_1.drop(columns=["color_id", "code", "geometry"]),
            how="left",
            on="drainage_unit_id"
        )
        
        self.afvoergebied_1_gdf["downstream_order_code"] = self.afvoergebied_1_gdf["downstream_order_code"].fillna("")
        self.afvoergebied_1_gdf = self.afvoergebied_1_gdf.sort_values("downstream_order_code", ascending=True)
        self.afvoergebied_2_gdf = (
            self.afvoergebied_1_gdf
            .dissolve(by="downstream_edges")
            .drop(columns="downstream_order_code")
            .reset_index()
        )

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "afvoergebied_1_gdf",
                    "afvoergebied_2_gdf",
                ]
            )

        self.afvoergebied_1 = self.afvoergebied_0.copy()
        self.afvoergebied_1 = dataarray_from_gdf(
            self.afvoergebied_1, 
            self.afvoergebied_1_gdf, 
            "afvoergebied_1"
        )
        self.afvoergebied_2 = self.afvoergebied_0.copy()
        self.afvoergebied_2 = dataarray_from_gdf(
            self.afvoergebied_2, 
            self.afvoergebied_2_gdf, 
            "afvoergebied_2"
        )

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "afvoergebied_1",
                    "afvoergebied_2",
                ]
            )
        return self.afvoergebied_2_gdf


    def aggregate_afvoergebied_tot_4(self):
        self.afvoergebied_3 = None
        self.afvoergebied_3_gdf = None
        self.afvoergebied_4 = None
        self.afvoergebied_4_gdf = None

        logging.info("   x aggregation/lumping of drainage units to basins")
        logging.info("     - aggregate sub drainage units to basins (level 3 and 4)")
        self.afvoergebied_3_gdf = self.afvoergebied_2_gdf.dissolve(by="order_code").reset_index()
        random_color_id = np.random.randint(0, 25, size=len(self.afvoergebied_3_gdf))
        self.afvoergebied_3_gdf["color_id"] = random_color_id
        self.afvoergebied_3_gdf["order_code_no"] = self.afvoergebied_3_gdf.order_code.str[:6]

        self.afvoergebied_4_gdf = self.afvoergebied_3_gdf.dissolve("order_code_no").reset_index()
        random_color_id = np.random.randint(0, 25, size=len(self.afvoergebied_4_gdf))
        self.afvoergebied_4_gdf["color_id"] = random_color_id

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "afvoergebied_3_gdf",
                    "afvoergebied_4_gdf",
                ]
            )

        self.afvoergebied_3 = self.afvoergebied_0.copy()
        self.afvoergebied_3 = dataarray_from_gdf(
            self.afvoergebied_3, 
            self.afvoergebied_3_gdf, 
            "afvoergebied_3"
        )

        self.afvoergebied_4 = self.afvoergebied_0.copy()
        self.afvoergebied_4 = dataarray_from_gdf(
            self.afvoergebied_4, 
            self.afvoergebied_4_gdf, 
            "afvoergebied_4"
        )

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "afvoergebied_3",
                    "afvoergebied_4",
                ]
            )
        return self.afvoergebied_4_gdf


    def generate_folium_map(self, html_file_name:str="", **kwargs):
        if html_file_name == '':
            html_file_name = self.name + "_afvoergebied"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            **kwargs
        )
        return self.folium_map
