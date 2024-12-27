import ast
import logging
import os
import random
import webbrowser
from pathlib import Path
import xarray

import folium
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.folium_utils import add_basemaps_to_folium_map, add_categorized_lines_to_map, add_labels_to_points_lines_polygons
from ..utils.network_functions import find_nodes_edges_for_direction


class GeneratorDrainageUnits(GeneratorBasis):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    dir_basis_data: str = "0_basisdata"
    dir_inter_results: str = "1_tussenresultaat"
    dir_results: str = "2_resultaat"

    direction: str = "upstream"
    read_results: bool = False
    write_results: bool = False

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed: gpd.GeoDataFrame = None

    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed: gpd.GeoDataFrame = None

    overige_watergangen_connections: gpd.GeoDataFrame = None

    ghg: xarray.Dataset = None
    afwateringseenheden: gpd.GeoDataFrame = None
    outflow_nodes: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    network_positions: dict = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: str = None

    def generate_folium_map(
        self, 
        base_map="Light Mode", 
        html_file_name=None, 
        open_html=False, 
        order_labels=False
    ):
        # Make figure
        outflow_nodes_4326 = self.outflow_nodes.to_crs(4326)

        m = folium.Map(
            location=[
                outflow_nodes_4326.geometry.y.mean(),
                outflow_nodes_4326.geometry.x.mean(),
            ],
            zoom_start=12,
            tiles=None,
        )

        if "order_no" in self.edges.columns:
            add_categorized_lines_to_map(
                m=m,
                # feature_group=fg,
                lines_gdf=self.edges[self.edges["order_no"] > 1][
                    ["code", "order_no", "geometry"]
                ],
                layer_name="Orde-nummer watergangen",
                control=True,
                lines=True,
                line_color_column="order_no",
                line_color_cmap="hsv",
                label=False,
                line_weight=5,
                z_index=1,
            )

            if order_labels and "order_no" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"Orde-nummer watergangen (labels)",
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

            if order_labels and "order_code" in self.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"Orde-code watergangen (labels)",
                    control=True,
                    show=False,
                ).add_to(m)

                self.edges["order_code"].fill = ""

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["order_no"] > 1][
                        ["code", "order_code", "geometry"]
                    ],
                    column="order_code",
                    label_fontsize=8,
                    fg=fg,
                )
        else:
            folium.GeoJson(
                self.hydroobjecten.geometry,
                name="Watergangen",
                color="blue",
                fill_color="blue",
                zoom_on_click=True,
                show=False,
                z_index=2,
            ).add_to(m)

        fg = folium.FeatureGroup(
            name=f"Uitstroompunten RWS-water", control=True
        ).add_to(m)

        folium.GeoJson(
            self.outflow_nodes[self.outflow_nodes["order_no"]==2],
            name="Uitstroompunten RWS-wateren",
            marker=folium.Circle(
                radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=3,
        ).add_to(fg)

        add_labels_to_points_lines_polygons(
            gdf=self.outflow_nodes[self.outflow_nodes["order_no"]==2],
            column="order_code", 
            label_fontsize=8, 
            fg=fg
        )
        m = add_basemaps_to_folium_map(m=m, base_map=base_map)

        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m
        if html_file_name is None:
            html_file_name = self.name

        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m
