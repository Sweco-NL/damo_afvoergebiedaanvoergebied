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
from ..utils.folium_utils import add_basemaps_to_folium_map
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

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    network_positions: dict = None
    graph: nx.DiGraph = None

    folium_map: folium.Map = None
    folium_html_path: str = None

    def generate_folium_map(self, base_map="OpenStreetMap"):
        # Make figure
        outflow_nodes_4326 = self.outflow_nodes_all.to_crs(4326)

        m = folium.Map(
            location=[
                outflow_nodes_4326.geometry.y.mean(),
                outflow_nodes_4326.geometry.x.mean(),
            ],
            zoom_start=12,
            tiles=None,
        )

        folium.GeoJson(
            self.rws_wateren.geometry,
            name="RWS_Wateren",
            z_index=0,
        ).add_to(m)

        if "orde_nr" in self.edges.columns:
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.edges[self.edges["orde_nr"] > 1],
                layer_name="Orde watergangen",
                control=True,
                lines=True,
                line_color_column="orde_nr",
                label=True,
                label_column="orde_nr",
                label_decimals=0,
                line_weight=5,
                z_index=1,
            )
        folium.GeoJson(
            self.hydroobjecten.geometry,  # .buffer(10),
            name="Watergangen",
            color="blue",
            fill_color="blue",
            zoom_on_click=True,
            show=False,
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
        m = add_basemaps_to_folium_map(m=m, base_map=base_map)

        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m
        return m
