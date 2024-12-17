import logging
from pathlib import Path
import webbrowser

import folium
import geopandas as gpd
import networkx as nx
import pandas as pd
from pydantic import ConfigDict
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import (
    define_list_upstream_downstream_edges_ids,
    calculate_angles_of_edges_at_nodes,
    select_downstream_upstream_edges
)
from ..utils.folium_utils import (
    add_labels_to_points_lines_polygons,
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
)


class GeneratorOrderLevels(GeneratorBasis):
    """Module to generate partial networks and order levels for all water bodies,
    based on ..."""

    path: Path = None
    name: str = None
    base_dir: Path = None
    waterschap: str = None
    range_orde_code_min: int = None
    range_orde_code_max: int = None

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed: gpd.GeoDataFrame = None

    # overige_watergangen: gpd.GeoDataFrame = None
    # overige_watergangen_processed: gpd.GeoDataFrame = None

    rws_wateren: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    dead_end_nodes: gpd.GeoDataFrame = None
    dead_start_nodes: gpd.GeoDataFrame = None

    outflow_nodes_all: gpd.GeoDataFrame = None
    outflow_nodes: gpd.GeoDataFrame = None

    edges: gpd.GeoDataFrame = None
    nodes: gpd.GeoDataFrame = None
    graph: nx.DiGraph = None
    network_positions: dict = None

    folium_map: folium.Map = None
    folium_html_path: Path = None

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
                edges = pd.concat([edges, gdf_water_line.explode()])
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        self.network_positions = {n: [n[0], n[1]] for n in list(self.graph.nodes)}
        return self.nodes, self.edges, self.graph


    def define_list_upstream_downstream_edges_ids(self):
        self.nodes = define_list_upstream_downstream_edges_ids(
            node_ids=self.nodes.nodeID.values, 
            nodes=self.nodes, 
            edges=self.edges
        )


    def calculate_angles_of_edges_at_nodes(self):
        logging.info("   x calculate angles of edges to nodes")
        self.nodes, self.edges = calculate_angles_of_edges_at_nodes(
            nodes=self.nodes, edges=self.edges
        )
        return self.nodes


    def select_downstream_upstream_edges(self, min_difference_angle: str = 20.0):
        logging.info("   x find downstream upstream edges")
        self.nodes = select_downstream_upstream_edges(self.nodes)
        return self.nodes


    def find_end_points_hydroobjects(self, buffer_width=0.5, direction="upstream"):
        # Copy hydroobject data to new variable 'hydroobjects' and make dataframes with start and end nodes
        logging.info("   x find start and end nodes hydrobojects")

        dead_end_nodes = self.edges[~self.edges.node_end.isin(self.edges.node_start.values)].copy()
        dead_end_nodes["geometry"] = dead_end_nodes["geometry"].apply(lambda x: Point(x.coords[-1]))
        dead_end_nodes["nodeID"] = dead_end_nodes["node_end"]
        dead_end_nodes = dead_end_nodes[["code", "nodeID", "geometry"]].rename(columns={"code": "hydroobject_code"})

        dead_start_nodes = self.edges[~self.edges.node_start.isin(self.edges.node_end.values)].copy()
        dead_start_nodes["geometry"] = dead_start_nodes["geometry"].apply(lambda x: Point(x.coords[0]))
        dead_start_nodes["nodeID"] = dead_start_nodes["node_start"]
        dead_start_nodes = dead_start_nodes[["code", "nodeID", "geometry"]].rename(columns={"code": "hydroobject_code"})
        
        self.dead_end_nodes = dead_end_nodes.copy()
        self.dead_start_nodes = dead_start_nodes.copy()
        return self.dead_end_nodes, self.dead_start_nodes


    def generate_rws_code_for_all_outflow_points(self, buffer_rws=10.0):
        logging.info("   x generating order code for all outflow points")
        rws_wateren = self.rws_wateren.copy()
        rws_wateren.geometry = rws_wateren.geometry.buffer(buffer_rws)
        outflow_nodes = (
            self.dead_end_nodes.sjoin(rws_wateren[["geometry", "rws_code"]])
            .drop(columns="index_right")
            .reset_index(drop=True)
        )

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
                logging_message = f" XXX aantal uitstroompunten op RWS-water ({rws_code}) hoger dan range order_code waterschap"
                logging.debug(logging_message)

            if outflow_nodes_all is None:
                outflow_nodes_all = outflows
            else:
                outflow_nodes_all = pd.concat([outflow_nodes_all, outflows])
            logging_message = f"    - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            logging.debug(logging_message)

        logging_message = f"    - total no. outflow_point on RWS-waters for {self.name}: {len(outflow_nodes_all)}"
        logging.debug(logging_message)
        self.outflow_nodes_all = outflow_nodes_all.reset_index(drop=True)
        self.outflow_nodes_all["orde_nr"] = 2
        return self.outflow_nodes_all


    def export_results_to_gpkg(self):
        """Export results to geopackages in folder 1_tussenresultaat"""
        results_dir = Path(self.path, self.dir_inter_results)
        logging.info(f"  x export results")
        for layer in [
            "outflow_nodes_all",
            "edges",
            "nodes",
        ]:
            result = getattr(self, layer)
            if result is None:
                logging.debug(f"    - {layer} not available")
            else:
                logging.debug(f"    - {layer} ({len(result)})")
                result.to_file(Path(results_dir, f"{layer}.gpkg"))
    

    def generate_folium_map(
        self,
        html_file_name: str = None,
        include_areas: bool = True,
        width_edges: float = 10.0,
        opacity_edges: float = 0.5,
        open_html: bool = False,
        base_map: str = "OpenStreetMap",
        order_labels: bool = False
    ):
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
                lines_gdf=self.edges[self.edges["orde_nr"] > 1][
                    ["code", "orde_nr", "geometry"]
                ],
                layer_name="Orde-nummer watergangen",
                control=True,
                lines=True,
                line_color_column="orde_nr",
                line_color_cmap="hsv",
                label=False,
                line_weight=5,
                z_index=1,
            )

            if order_labels:
                fg = folium.FeatureGroup(
                    name=f"Orde-nummer watergangen (labels)", 
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.edges[self.edges["orde_nr"] > 1][
                        ["code", "orde_nr", "geometry"]
                    ], 
                    column="orde_nr", 
                    label_fontsize=8, 
                    label_decimals=0,
                    fg=fg
                )
        
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
            self.outflow_nodes_all[self.outflow_nodes_all["orde_nr"]==2],
            name="Uitstroompunten RWS-wateren",
            marker=folium.Circle(
                radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=3,
        ).add_to(fg)

        add_labels_to_points_lines_polygons(
            gdf=self.outflow_nodes_all[self.outflow_nodes_all["orde_nr"]==2], 
            column="orde_code", 
            label_fontsize=8, fg=fg
        )
        
        folium.GeoJson(
            self.outflow_nodes_all[self.outflow_nodes_all["orde_nr"]>2],
            name="Overige uitstroompunten",
            marker=folium.Circle(
                radius=15,
                fill_color="orange",
                fill_opacity=0.8,
                color="orange",
                weight=3,
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=2,
        ).add_to(m)

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
