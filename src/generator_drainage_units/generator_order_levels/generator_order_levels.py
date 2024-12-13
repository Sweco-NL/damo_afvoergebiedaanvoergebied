import logging
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import pandas as pd
from pydantic import BaseModel, ConfigDict
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.general_functions import shorten_line_two_vertices, line_to_vertices
from ..utils.folium_utils import (
    add_labels_to_points_lines_polygons,
    add_basemaps_to_folium_map,
    add_lines_to_map,
    add_categorized_lines_to_map
)


class GeneratorOrderLevels(GeneratorBasis):
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
    hydroobjecten_processed: gpd.GeoDataFrame = None

    overige_watergangen: gpd.GeoDataFrame = None
    overige_watergangen_processed: gpd.GeoDataFrame = None

    rws_wateren: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    dead_end_nodes: gpd.GeoDataFrame = None
    outflow_nodes_all: gpd.GeoDataFrame = None
    outflow_nodes: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None
    network_positions: dict = None

    folium_map: folium.Map = None


    def create_graph_from_network(
        self, water_lines=["rivieren", "hydroobjecten", "hydroobjecten_extra"]
    ):
        if water_lines is None:
            water_lines = ["hydroobjecten"]
        logging.info("   x create network graph")
        edges = None
        for water_line in water_lines:
            gdf_water_line = getattr(self, water_line)
            if gdf_water_line is None:
                continue
            if edges is None:
                edges = gdf_water_line.explode()
            else:
                edges = pd.concat(
                    [edges, gdf_water_line.explode()]
                )
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        self.network_positions = {n: [n[0], n[1]] for n in list(self.graph.nodes)}
        return self.nodes, self.edges, self.graph


    def find_end_points_hydroobjects(self, buffer_width=0.5):
        # Copy hydroobject data to new variable 'hydroobjects' and make dataframes with start and end nodes
        hydroobjects = self.hydroobjecten[["code", "geometry"]].explode()
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

        logging.info(f"     - start and end nodes defined")

        # Determine dead-end end nodes (using buffer to deal with inaccuracies in geometry hydroobjects) and make new dataframe
        start_nodes["geometry"] = start_nodes.buffer(buffer_width)
        end_nodes["geometry"] = end_nodes.buffer(buffer_width)

        # Determine where start and end nodes overlap
        dead_end_nodes = start_nodes.sjoin(end_nodes, how="right")
        dead_end_nodes = hydroobjects.loc[
            dead_end_nodes[dead_end_nodes.index_left.isna()].index, ["code", "end_node"]
        ].copy()
        dead_end_nodes = dead_end_nodes.rename(columns={"end_node": "geometry"})

        dead_end_nodes["geometry"] = dead_end_nodes["geometry"].buffer(buffer_width)
        sjoin = dead_end_nodes.sjoin(hydroobjects[["code", "geometry"]])
        no_dead_end_node_ids = sjoin[sjoin["code_left"] != sjoin["code_right"]].index

        self.dead_end_nodes = dead_end_nodes[
            ~dead_end_nodes.index.isin(no_dead_end_node_ids)
        ].copy()
        self.dead_end_nodes.geometry = self.dead_end_nodes.geometry.centroid

        logging.info(f"     - dead end nodes defined")

        return self.dead_end_nodes

    def generate_rws_code_for_all_outflow_points(self, buffer_rws=10.0):
        logging.info(f"     - generating order code for all outflow points")
        rws_wateren = self.rws_wateren.copy()
        rws_wateren.geometry = rws_wateren.geometry.buffer(buffer_rws)
        outflow_nodes = (
            self.dead_end_nodes.sjoin(rws_wateren[["geometry", "rws_code"]])
            .drop(columns="index_right")
            .reset_index(drop=True)
        )

        outflow_nodes = outflow_nodes.rename(columns={"code": "hydroobject_code"})
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
                f"     - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            )

        logging.info(
            f"     - total no. outflow_point on RWS-waters for {self.name}: {len(outflow_nodes_all)}"
        )
        self.outflow_nodes_all = outflow_nodes_all.reset_index(drop=True)
        return self.outflow_nodes_all


    def generate_order_levels_to_outflow_nodes_edges(self, max_order: int = 10):
        logging.info(f"   x Generating order levels for all edges")
        self.edges['orde_nr'] = 0

        self.outflow_nodes = self.outflow_nodes_all.copy()
        self.outflow_nodes = self.outflow_nodes.rename(columns={'hydroobject_code':'hydroobject_code_in'})
        self.outflow_nodes['hydroobject_code_out'] = None
        self.outflow_nodes["orde_nr"] = 1
        for orde_nr in range(1,max_order):
            print(orde_nr)
            outflow_nodes_order = self.outflow_nodes[self.outflow_nodes['orde_nr'] == orde_nr].copy()
            outflow_nodes_order['geometry'] = outflow_nodes_order.buffer(0.0001).to_crs(28992)
            
            end_points = self.edges.copy()
            end_points['geometry'] = end_points['geometry'].apply(lambda geom: Point(geom.coords[-1]))
            sel_edges = gpd.sjoin(
                outflow_nodes_order,
                end_points[['code', 'geometry']],
                how="left",
                predicate="intersects",
            )
            codes = list(sel_edges['code'].values)

            while codes:
                # Find hydroobjects in gdf_hydroobjects whose endpoints match the start points of the current hydroobjects
                start_points = self.edges[self.edges['code'].isin(codes)]['geometry'].apply(lambda x: x.coords[0]).tolist()
                next_edges = self.edges[self.edges['geometry'].apply(lambda x: x.coords[-1] in start_points)].copy()

                # Add a column 'next_edge' to store the hydroobject it ends at
                next_edges['next_edge'] = next_edges['geometry'].apply(
                    lambda x: [code for code, geom in zip(self.edges['code'], self.edges['geometry']) if geom.coords[0] == x.coords[-1]]
                )
                self.edges.loc[self.edges['code'].isin(codes), 'orde_nr'] = orde_nr + 1
                
                # Identify duplicated 'next_edge' values
                next_edges_1 = next_edges[~next_edges.duplicated('next_edge', keep=False)]
                next_edges_2 = next_edges[next_edges.duplicated('next_edge', keep=False)]
                
                # Create a new outflow_nodes based on multiple_next_hydro
                for next_edge in next_edges_2['next_edge'].explode().unique():
                    # Get hydroobjects with the same next_hydro
                    next_edge_same = next_edges_2[next_edges_2['next_edge'].explode() == next_edge]
                    
                    # Aggregate the hydroobject codes
                    # hydroobjects_codes = ','.join(hydroobjects_with_same_next_hydro['CODE'].astype(str))
                    edges_codes = list(next_edge_same['code'].astype(str).values)
                    
                    # Create a new outflow_node
                    new_outflow_node = {
                        'geometry': Point(next_edge_same.iloc[0]['geometry'].coords[-1]),
                        'hydroobject_code_out': next_edge,
                        'hydroobject_code_in': edges_codes,
                        'orde_nr': int(orde_nr + 1)
                    }
                    new_outflow_node = gpd.GeoDataFrame(
                        pd.DataFrame([new_outflow_node]), 
                        geometry='geometry', 
                        crs=self.outflow_nodes.crs
                    )
                    # Append the new outflow_node to outflow_nodes
                    self.outflow_nodes = pd.concat([self.outflow_nodes, new_outflow_node], ignore_index=True)

                codes = next_edges_1['code'].tolist()
        return self.edges, self.outflow_nodes


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

        if 'orde_nr' in self.edges.columns:
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.edges[self.edges['orde_nr']>1],
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
