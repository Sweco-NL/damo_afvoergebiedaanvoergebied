import logging
import webbrowser
from pathlib import Path
import time
from tqdm import tqdm

import folium
import geopandas as gpd
import numpy as np
import pandas as pd
import pyflwdir
import rioxarray
import xarray
from geocube.api.core import make_geocube
from numba import njit
from pydantic import ConfigDict
from rasterio.enums import Resampling
from scipy.ndimage import distance_transform_edt

from ..generator_basis import GeneratorBasis
from ..utils.folium_utils import (
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
    add_graduated_raster_to_map,
    add_labels_to_points_lines_polygons,
)


class GeneratorDrainageUnits(GeneratorBasis):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basisdata: str = "0_basisdata"
    dir_results: str = "1_resultaat"

    read_results: bool = False
    write_results: bool = False
    crs: int = 28992
    
    required_results: list[str] = [
        "hydroobjecten",
        "hydroobjecten_processed_0", 
        "overige_watergangen",
        "overige_watergangen_processed_4", 
        "outflow_nodes",
        "edges", 
    ]

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    
    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed_4: gpd.GeoDataFrame = None

    outflow_nodes: gpd.GeoDataFrame = None

    ghg_file_name: str = None
    ghg: xarray.Dataset = None

    all_waterways_0: gpd.GeoDataFrame = None
    all_waterways_1: gpd.GeoDataFrame = None

    ghg_waterways: xarray.Dataset = None
    ghg_waterways_distance: xarray.Dataset = None

    ghg_processed: xarray.Dataset = None

    flw: pyflwdir.FlwdirRaster = None
    
    drainage_units_0: xarray.Dataset = None
    drainage_units_1: xarray.Dataset = None
    
    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    
    folium_map: folium.Map = None
    folium_html_path: str = None


    def read_ghg(self, ghg_file_name: str):
        """Read GHG file with different name then GHG

        Read GHG file with different name then GHG

        Parameters
        ----------
        ghg_file_name : str
            Name of GHG netcdf file ('GHG.nc')

        Returns
        -------
        GHG xarray data array
            GHG DataArray
        """
        self.ghg_file_name = ghg_file_name
        self.ghg = rioxarray.open_rasterio(Path(self.path, self.dir_basisdata, ghg_file_name))
        self.ghg.name = "GHG_2000-2010_L1"
        return self.ghg


    def preprocess_ghg(self, resolution=2.0, depth_waterways=1.0, buffer_waterways=2.5, smooth_distance=25.0):
        # resample to new resolution (m)
        logging.info("   x preprocessing GHG data")
        logging.info("     - resampling data to new resolution")
        old_resolution = 25.0
        upscale_factor = old_resolution / resolution
        new_width = int(np.ceil(self.ghg.rio.width * upscale_factor))
        new_height = int(np.ceil(self.ghg.rio.height * upscale_factor))

        ghg_processed = self.ghg.rio.reproject(
            self.ghg.rio.crs,
            shape=(new_height, new_width),
            resampling=Resampling.bilinear,
        )

        # add depth at location hydroobjects and other waterways (m)
        logging.info("     - add depth at waterways")
        self.all_waterways_0 = pd.concat([
            self.hydroobjecten_processed_0[["code", "geometry"]],
            self.overige_watergangen_processed_4[["code", "geometry"]]
        ]).reset_index(drop=True).reset_index(drop=True)
        
        self.all_waterways_0["depth_waterways"] = depth_waterways
        self.all_waterways_0["drainage_unit_id"] = self.all_waterways_0.index
        self.all_waterways_0.geometry = self.all_waterways_0.geometry.buffer(buffer_waterways)

        ghg_waterways = make_geocube(
            vector_data=self.all_waterways_0,
            measurements=["depth_waterways"],
            like=ghg_processed,
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

        ghg_processed.data = ghg_processed.data - ghg_waterways.data
        ghg_processed.data[ghg_processed.data<-900.0] = -999.99

        self.ghg_processed = ghg_processed.copy()
        if self.write_results:
            self.all_waterways_0.to_file(
                Path(self.dir_results, "all_waterways_0.gpkg"), 
                layer="all_waterways_0"
            )

            # netcdf_file_path = Path(self.dir_results, "ghg_processed.nc")
            # encoding = {
            #     self.ghg_file_name.replace(".nc", "").replace(".NC", ""): {
            #         'dtype': 'float32',
            #         'zlib': True,
            #         'complevel': 9,
            #     },
            # }
            # self.ghg_processed.to_netcdf(
            #     netcdf_file_path, 
            #     encoding=encoding
            # )
        return ghg_processed


    def generate_drainage_units(self, iterations=2000):
        logging.info("   x generate drainage units for each waterway")
        # create raster with unique id of each waterway
        logging.info("     - give each waterway an unique id")
        if self.all_waterways_0 is None or self.ghg_processed is None:
            raise ValueError("   x run preprocess_ghg to preprocess the data")
        
        self.drainage_units_0 = make_geocube(
            vector_data=self.all_waterways_0,
            measurements=["drainage_unit_id"],
            like=self.ghg_processed,
        )["drainage_unit_id"].fillna(-1)

        if self.write_results:
            self.drainage_units_0.name = 'drainage_units_0'
            netcdf_file_path = Path(self.dir_results, "drainage_units_0_test.nc")
            encoding = {
                'drainage_units_0': {
                    'dtype': 'float32',
                    'zlib': True,
                    'complevel': 9,
                },
            }
            self.drainage_units_0.to_netcdf(
                netcdf_file_path, 
                encoding=encoding
            )

        # create pyflwdir object
        logging.info("     - create pyflwdir object to calculate downstream direction")
        flw = pyflwdir.from_dem(
            data=self.ghg_processed.data[0],
            nodata=self.ghg_processed._FillValue,
            transform=self.ghg_processed.rio.transform(),
            latlon=False,
        )

        # get upstream values
        def get_upstream_values(
            flw_mask: np.ndarray, 
            flw_idxs_ds: np.ndarray, 
            drainage_units_flat: np.ndarray, 
            iterations: int
        ):
            for i in range(iterations):
                time_start = time.time()
                upstream_values = drainage_units_flat.copy()
                upstream_values[flw_mask] = upstream_values[flw_idxs_ds[flw_mask]]
                
                new_filled_cells = (drainage_units_flat == -1.0) & (upstream_values != -1.0)
                drainage_units_flat = np.where(new_filled_cells, upstream_values, drainage_units_flat)

                number_new_filled_cells = new_filled_cells.sum()

                if number_new_filled_cells == 0:
                    print("     * break at iteration: ", i)
                    break
                print1 = f"  * iteration: {i}/{iterations}"
                print2 = f" | number new cells: {number_new_filled_cells}"
                print3 = f"({round(time.time()-time_start, 2)} seconds)"
                print(print1 + print2 + print3, end="\r")
            return drainage_units_flat
        
        logging.info("     - get upstream area of each waterway")
        drainage_units_flat_basis = flw._check_data(self.drainage_units_0.data, "data")
        drainage_units_flat_new = get_upstream_values(
            flw.mask, 
            flw.idxs_ds, 
            drainage_units_flat_basis, 
            iterations
        )
        self.drainage_units_0.data = drainage_units_flat_new.reshape(
            self.drainage_units_0.data.shape
        )
        self.drainage_units_0.name = "drainage_units_0"
        if self.write_results:
            netcdf_file_path = Path(self.dir_results, "drainage_units_0.nc")
            encoding = {
                'drainage_units_0': {
                    'dtype': 'float32',
                    'zlib': True,
                    'complevel': 9,
                },
            }
            self.drainage_units_0.to_netcdf(
                netcdf_file_path, 
                encoding=encoding
            )
        return self.drainage_units_0


    def aggregate_drainage_units(self):
        logging.info("   x aggregation/lumping of drainage units: aggregate 'overige watergangen'")
        logging.info("     - define new drainage_unit_ids for all 'overige watergangen'")

        all_waterways_1 = self.all_waterways_0.merge(
            self.edges[["code", "order_code"]],
            how="left",
            on="code"
        ).merge(
            self.overige_watergangen_processed_4[["code", "downstream_order_code"]],
            how="left",
            on="code"
        )
        all_waterways_1 = all_waterways_1.merge(
            all_waterways_1[["drainage_unit_id", "order_code"]].rename(
                columns={"drainage_unit_id": "new_drainage_unit_id"}
            ).dropna(subset="order_code"),
            how="left",
            left_on="downstream_order_code",
            right_on="order_code",
            suffixes=["", "_r"]
        ).drop(columns=["order_code_r"])
        
        all_waterways_1 = all_waterways_1.drop_duplicates(
            subset=["geometry"], 
            keep="first"
        )
        all_waterways_1 = all_waterways_1.dropna(subset="new_drainage_unit_id")
        all_waterways_1 = all_waterways_1[all_waterways_1['order_code'].isna()]
        all_waterways_1["new_drainage_unit_id"] = all_waterways_1["new_drainage_unit_id"].astype(int)

        self.all_waterways_1 = all_waterways_1.copy()
        if self.write_results and self.all_waterways_1 is not None:
            self.all_waterways_1.to_file(
                Path(self.dir_results, "all_waterways_1.gpkg"), 
                layer="all_waterways_1"
            )

        # translations_drainage_unit_id = all_waterways_1.groupby("new_drainage_unit_id").agg({"drainage_unit_id": list}).reset_index()
        translations_drainage_unit_id = all_waterways_1[["drainage_unit_id", "new_drainage_unit_id"]]

        logging.info(f"     - make the translations from {len(all_waterways_1)} to {len(translations_drainage_unit_id)}")
        drainage_units_1 = self.drainage_units_0.copy()
        drainage_units_1_flat = drainage_units_1[0].data.flatten()

        @njit
        def replace_values_in_array(array, old_values, new_values):
            perc_array = 0
            index = 0
            for old_value, new_value in zip(old_values, new_values):
                array[array==old_value] = new_value
                index += 1
                if index > perc_array:
                    print(index, "/", len(old_values))
                    perc_array += len(old_values)/20.0
            return array
        
        drainage_units_1_flat_new = replace_values_in_array(
            drainage_units_1_flat, 
            translations_drainage_unit_id["drainage_unit_id"].values,
            translations_drainage_unit_id["new_drainage_unit_id"].values,
        )
        drainage_units_1[0].data = np.reshape(drainage_units_1_flat_new, drainage_units_1[0].data.shape)
        
        self.drainage_units_1 = drainage_units_1.copy()
        self.drainage_units_1.name = "drainage_units_1"

        if self.write_results:
            netcdf_file_path = Path(self.dir_results, "drainage_units_1.nc")
            encoding = {
                'drainage_units_1': {
                    'dtype': 'float32',
                    'zlib': True,
                    'complevel': 9,
                },
            }
            self.drainage_units_1.to_netcdf(
                netcdf_file_path, 
                encoding=encoding
            )
        return self.drainage_units_1


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
            add_categorized_lines_to_map(
                m=m,
                # feature_group=fg,
                lines_gdf=self.edges[self.edges["order_no"] > 1][
                    ["code", "order_no", "geometry"]
                ].sort_values("order_no", ascending=False),
                layer_name="A/B Watergangen: orde-nummers",
                control=True,
                show=False,
                lines=True,
                line_color_column="order_no",
                line_color_cmap=None,
                label=False,
                line_weight=4,
                z_index=1,
            )

            if order_labels and "order_no" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"A/B Watergangen: orde-nummers (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["order_no"] > 1][
                        ["code", "order_no", "geometry"]
                    ],
                    column="order_no",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )
        else:
            folium.GeoJson(
                self.hydroobjecten_processed_0.geometry,
                name="A/B Watergangen",
                color="blue",
                fill_color="blue",
                zoom_on_click=True,
                show=False,
                z_index=2,
            ).add_to(m)

        if self.overige_watergangen_processed_4 is not None:
            folium.GeoJson(
                self.overige_watergangen_processed_4.geometry,
                name="C-watergangen",
                color="blue",
                weight=1,
                zoom_on_click=True,
                show=True,
                z_index=3,
            ).add_to(m)

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
                show=False,
                z_index=3,
            )

        if self.drainage_units_0 is not None:
            drainage_units_0 = self.drainage_units_0.where(self.drainage_units_0 > -1.0)
            drainage_units_0 = drainage_units_0.rio.write_crs(self.crs)
            add_graduated_raster_to_map(
                m=m,
                raster=drainage_units_0,
                layer_name="Afwateringseenheden (per watergangsdeel)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_0.data.max()),
                legend=False,
                opacity=0.75,
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
                layer_name="Afwateringseenheden (hydroobjecten)",
                unit="unique id",
                control=True,
                vmin=0,
                vmax=int(self.drainage_units_1.data.max()),
                legend=False,
                opacity=0.75,
                show=True,
                dx=dx,
                dy=dy,
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
                legend=True,
                opacity=0.75,
                show=False,
                cmap='Spectral_r',
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

