import logging
import time
import webbrowser
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
from ..utils.folium_utils import (
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
    add_graduated_raster_to_map,
    add_labels_to_points_lines_polygons,
)


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

    method: str = "pyflwdir"

    read_results: bool = False
    write_results: bool = False
    crs: int = 28992
    
    required_results: list[str] = [
        "hydroobjecten",
        "hydroobjecten_processed_0", 
        "overige_watergangen_processed_4", 
        "potential_culverts_5", 
        "outflow_nodes",
        "edges", 
    ]

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    
    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed_4: gpd.GeoDataFrame = None
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

        # add depth at location hydroobjects and other waterways (m)
        logging.info("     - add depth at waterways")
        if self.overige_watergangen_processed_4 is None:
            self.all_waterways_0 = self.edges[["code", "geometry"]].reset_index(drop=True)
        else:
            self.all_waterways_0 = pd.concat([
                self.edges[["code", "geometry"]],
                self.overige_watergangen_processed_4[["code", "geometry"]]
            ]).reset_index(drop=True)
        
        logging.info("     - give each waterway an unique id")
        self.all_waterways_0["depth_waterways"] = depth_waterways
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
        # gdf = gdf[gdf["drainage_unit_id"] >= 0]
        gdf = gdf.set_crs(self.hydroobjecten.crs)
        random_color_id = np.random.randint(0, 25, size=len(gdf))
        gdf["color_id"] = random_color_id

        gdf["drainage_unit_area"] = gdf["geometry"].area
        gdf = gdf.explode().reset_index(drop=True)
        gdf["part_count"] = gdf[["drainage_unit_id"]].groupby("drainage_unit_id").transform("count").reset_index()
        gdf = gdf[gdf["drainage_unit_id"] > -1]

        # area_lim = self.new_resolution * self.new_resolution * 1.5
        # gdf = gdf.loc[(gdf.geometry.area>=area_lim) | (gdf.part_count<2)]
        # gdf.geometry = remove_holes_from_polygons(gdf.geometry, min_area=area_lim)
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
        drainage_units_0_gdf = gpd.sjoin(
            self.drainage_units_0_gdf,
            self.all_waterways_1.drop(columns=["color_id", "drainage_unit_id"]),
            how="inner",
            predicate="intersects"
        )
        drainage_units_0_gdf = drainage_units_0_gdf.merge(
            self.all_waterways_1[["geometry", "downstream_edges"]].rename(columns={"geometry": "waterway_geometry"}), 
            how="left", 
            on="downstream_edges", 
        )
        drainage_units_0_gdf['overlap_length'] = (
            drainage_units_0_gdf
            .geometry
            .intersection(
                drainage_units_0_gdf.waterway_geometry
            ).length
        )
        drainage_units_0_gdf = drainage_units_0_gdf.loc[
            drainage_units_0_gdf.groupby('drainage_unit_id')['overlap_length'].idxmax()
        ]
        self.drainage_units_0_gdf = (
            drainage_units_0_gdf
            .reset_index(drop=True)
            .drop(columns=["index_right", "waterway_geometry", "overlap_length"])
        )
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
        
        self.drainage_units_0 = dataarray_from_gdf(
            self.drainage_units_0, 
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


    def generate_folium_map(
        self,
        base_map="Light Mode",
        html_file_name=None,
        open_html=False,
        order_labels=False,
        zmin=None,
        zmax=None,
        dx=0.0,
        dy=0.0,
    ):
        # Make figure
        m = folium.Map(
            location=[
                self.hydroobjecten.geometry.centroid.to_crs(4326).y.mean(),
                self.hydroobjecten.geometry.centroid.to_crs(4326).x.mean(),
            ],
            zoom_start=12,
            tiles=None,
        )

        if self.outflow_nodes is not None:
            fg = folium.FeatureGroup(
                name=f"Uitstroompunten in RWS-wateren", 
                control=True,
                z_index=0
            ).add_to(m)

            folium.GeoJson(
                self.outflow_nodes[self.outflow_nodes["order_no"] == 2],
                name="Uitstroompunten RWS-wateren",
                marker=folium.Circle(
                    radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
                ),
                highlight_function=lambda x: {"fillOpacity": 0.8},
                zoom_on_click=True,
            ).add_to(fg)

            add_labels_to_points_lines_polygons(
                gdf=self.outflow_nodes[self.outflow_nodes["order_no"] == 2],
                column="order_code",
                label_fontsize=8,
                fg=fg,
            )

        folium.GeoJson(
            self.hydroobjecten.geometry,
            name="AB-Watergangen",
            color="blue",
            weight=4,
            zoom_on_click=True,
            show=True,
            z_index=2,
        ).add_to(m)

        if "order_no" in self.edges.columns:
            edges = self.edges[self.edges["order_no"] > 1][
                ["code", "order_no", "order_code", "geometry"]
            ].sort_values("order_no", ascending=False)

            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges,
                layer_name="A/B Watergangen: orde-nummers",
                control=True,
                show=False,
                lines=True,
                line_color_column="order_no",
                line_color_cmap="hsv_r",
                label=False,
                line_weight=5,
                z_index=1,
            )

            if order_labels and "order_code" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"A/B Watergangen: orde-code (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["order_no"] > 1][
                        ["code", "order_code", "geometry"]
                    ],
                    column="order_code",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )

        if self.overige_watergangen is not None:
            folium.GeoJson(
                self.overige_watergangen.geometry,
                name="C-Watergangen - Zonder duikers",
                color="#0287c3",
                weight=2,
                z_index=0,
                show=False,
            ).add_to(m)

        if self.potential_culverts_5 is not None:
            folium.GeoJson(
                self.potential_culverts_5.geometry,
                name="C-Watergangen - Gevonden Duikers",
                color="red",
                weight=3,
                z_index=1,
                show=False,
            ).add_to(m)

        if self.overige_watergangen_processed_4 is not None:
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.overige_watergangen_processed_4[
                    ["outflow_node", "geometry"]
                ],
                layer_name=f"C-Watergangen - gegroepeerd",
                control=True,
                lines=True,
                line_color_column="outflow_node",
                line_color_cmap=None,
                show=True,
                z_index=3,
            )

        if self.ghg is not None:
            add_graduated_raster_to_map(
                m=m,
                raster=self.ghg,
                layer_name="GHG",
                unit="m NAP",
                control=True,
                vmin=self.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
                vmax=self.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
                legend=False,
                opacity=1.00,
                show=False,
                cmap='Spectral_r',
                dx=dx,
                dy=dy,
            )

        if self.ghg_processed is not None:
            add_graduated_raster_to_map(
                m=m,
                raster=self.ghg_processed,
                layer_name="GHG-processed",
                unit="m NAP",
                control=True,
                vmin=self.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
                vmax=self.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
                legend=False,
                opacity=1.00,
                show=False,
                cmap='Spectral_r',
                dx=dx,
                dy=dy,
            )
        if self.ghg_processed_adapt is not None:
            add_graduated_raster_to_map(
                m=m,
                raster=self.ghg_processed_adapt,
                layer_name="GHG-processed-adapt",
                unit="m NAP",
                control=True,
                vmin=self.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
                vmax=self.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
                legend=False,
                opacity=1.00,
                show=False,
                cmap='Spectral_r',
                dx=dx,
                dy=dy,
            )

        if self.drainage_units_0 is not None:
            drainage_units_0 = self.drainage_units_0.where(self.drainage_units_0 > -1.0)
            drainage_units_0 = drainage_units_0.rio.write_crs(self.crs)
            add_graduated_raster_to_map(
                m=m,
                raster=drainage_units_0,
                layer_name="Afwateringseenheden (A/B/C watergangen)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_0.data.max()),
                legend=False,
                opacity=1.0,
                show=True,
                dx=dx,
                dy=dy,
            )

        if self.drainage_units_1 is not None:
            drainage_units_1 = self.drainage_units_1.where(self.drainage_units_1 > -1.0)
            drainage_units_1 = drainage_units_1.rio.write_crs(self.crs)
            add_graduated_raster_to_map(
                m=m,
                raster=drainage_units_1,
                layer_name="Afwateringseenheden (A/B watergangen)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_1.data.max()),
                legend=False,
                opacity=1.0,
                show=False,
                dx=dx,
                dy=dy,
            )
        
        if self.drainage_units_2 is not None:
            drainage_units_2 = self.drainage_units_2.where(self.drainage_units_2 > -1.0)
            drainage_units_2 = drainage_units_2.rio.write_crs(self.crs)
            add_graduated_raster_to_map(
                m=m,
                raster=drainage_units_2,
                layer_name="Afwateringseenheden (orde-code)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_2.data.max()),
                legend=False,
                opacity=1.0,
                show=False,
                dx=dx,
                dy=dy,
            )

        if self.drainage_units_3 is not None:
            drainage_units_3 = self.drainage_units_3.where(self.drainage_units_3 > -1.0)
            drainage_units_3 = drainage_units_3.rio.write_crs(self.crs)
            add_graduated_raster_to_map(
                m=m,
                raster=drainage_units_3,
                layer_name="Afwateringseenheden (stroomgebied)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_3.data.max()),
                legend=False,
                opacity=1.0,
                show=False,
                dx=dx,
                dy=dy,
            )

        m = add_basemaps_to_folium_map(m=m, base_map=base_map)
        folium.LayerControl(collapsed=False).add_to(m)
        m.add_child(folium.plugins.MeasureControl())

        self.folium_map = m
        if html_file_name is None:
            html_file_name = self.name + "_drainage_units"

        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m

