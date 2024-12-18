import logging
from pathlib import Path
import webbrowser

import folium
import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import (
    define_list_upstream_downstream_edges_ids,
    calculate_angles_of_edges_at_nodes,
    select_downstream_upstream_edges,
    find_node_edge_ids_in_directed_graph
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
    range_order_code_min: int = None
    range_order_code_max: int = None

    hydroobjecten: gpd.GeoDataFrame = None
    hydroobjecten_processed: gpd.GeoDataFrame = None

    # overige_watergangen: gpd.GeoDataFrame = None
    # overige_watergangen_processed: gpd.GeoDataFrame = None

    rws_wateren: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    dead_end_nodes: gpd.GeoDataFrame = None
    dead_start_nodes: gpd.GeoDataFrame = None

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
        logging.info(f"   x create network graph ({len(self.edges)} edges, {len(self.nodes)} nodes)")
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
        dead_end_nodes = dead_end_nodes.drop_duplicates("nodeID", keep="first")

        dead_start_nodes = self.edges[~self.edges.node_start.isin(self.edges.node_end.values)].copy()
        dead_start_nodes["geometry"] = dead_start_nodes["geometry"].apply(lambda x: Point(x.coords[0]))
        dead_start_nodes["nodeID"] = dead_start_nodes["node_start"]
        dead_start_nodes = dead_start_nodes[["code", "nodeID", "geometry"]].rename(columns={"code": "hydroobject_code"})
        dead_start_nodes = dead_start_nodes.drop_duplicates("nodeID", keep="first")
        
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

        outflow_nodes_all_waters = None
        for rws_code in outflow_nodes.rws_code.unique():
            outflows = outflow_nodes[outflow_nodes.rws_code == rws_code].reset_index(
                drop=True
            )
            outflows["rws_code_no"] = self.range_order_code_min + outflows.index
            outflows["order_code"] = outflows.apply(
                lambda x: f"{x.rws_code}.{str(x.rws_code_no)}", axis=1
            )
            if outflows["rws_code_no"].max() > self.range_order_code_max:
                logging_message = f" XXX aantal uitstroompunten op RWS-water ({rws_code}) hoger dan range order_code waterschap"
                logging.debug(logging_message)

            if outflow_nodes_all_waters is None:
                outflow_nodes_all_waters = outflows.copy()
            else:
                outflow_nodes_all_waters = pd.concat([outflow_nodes_all_waters, outflows])
            logging_message = f"    - RWS-water ({rws_code}): {len(outflows)} outflow_points"
            logging.debug(logging_message)

        logging_message = f"    - total no. outflow_point on outside waters for {self.name}: {len(outflow_nodes_all_waters)}"
        logging.debug(logging_message)
        self.outflow_nodes = outflow_nodes_all_waters.reset_index(drop=True)
        self.outflow_nodes["order_no"] = 2

        return self.outflow_nodes


    def generate_orde_level_for_hydroobjects(self):
        logging.info(f"   x generate order levels for hydroobjects")
        self.outflow_nodes = self.outflow_nodes[self.outflow_nodes["order_no"]==2].copy()

        edges_left = self.edges.copy() #.drop(columns=["order_edge_no"]).copy()
        nodes_left = self.nodes.copy() #.drop(columns=["order_node_no"]).copy()
        order_no = 2
        
        edges_all_orders = None
        nodes_all_orders = None
        outflow_nodes_orders = self.outflow_nodes[self.outflow_nodes["order_no"]==2].copy()
        new_outflow_nodes_order = outflow_nodes_orders.copy()
        
        while not new_outflow_nodes_order.empty and order_no<100:
            logging.debug(f"    - order {order_no}: {len(new_outflow_nodes_order)} outflow nodes")
            outflow_nodes = new_outflow_nodes_order.copy()
            new_outflow_nodes_order = None
            outflow_nodes_edges = None
            outflow_nodes_nodes = None

            for i_node, outflow_node in outflow_nodes.reset_index(drop=True).iterrows():
                # logging.debug(f"      * {i_node+1}/{len(outflow_nodes)}")
                outflow_node_nodes_id, outflow_node_edges_code, new_outflow_nodes = find_node_edge_ids_in_directed_graph(
                    from_node_ids=edges_left.node_start.to_numpy(),
                    to_node_ids=edges_left.node_end.to_numpy(),
                    edge_ids=edges_left.code.to_numpy(),
                    node_ids=[outflow_node["nodeID"]],
                    search_node_ids=nodes_left.nodeID.to_numpy(),
                    search_edge_ids=edges_left.code.to_numpy(),
                    border_node_ids=None,
                    direction="upstream",
                    split_points=nodes_left,
                    order_first=True
                )
                outflow_node_edges = pd.DataFrame(data={
                    "code": outflow_node_edges_code[0],
                    "rws_code": outflow_node["rws_code"],
                    "rws_code_no": outflow_node["rws_code_no"],
                    "rws_order_code": outflow_node["order_code"],
                    "order_no": order_no,
                    "outflow_node": outflow_node["nodeID"],
                    "order_edge_no": range(len(outflow_node_edges_code[0])), 
                })
                outflow_node_edges = edges_left.merge(
                    outflow_node_edges, how="right", on="code"
                )
                outflow_node_edges = (
                    outflow_node_edges
                    .sort_values("order_no")
                    .drop_duplicates(subset="code", keep=False)
                )

                outflow_node_nodes = pd.DataFrame(data={
                    "nodeID": outflow_node_nodes_id[0],
                    "rws_code": outflow_node["rws_code"],
                    "rws_code_no": outflow_node["rws_code_no"],
                    "rws_order_code": outflow_node["order_code"],
                    "order_no": order_no,
                    "outflow_node": outflow_node["nodeID"],
                    "order_node_no": range(len(outflow_node_nodes_id[0])), 
                })
                outflow_node_nodes = self.nodes.merge(
                    outflow_node_nodes, how="right", on="nodeID"
                )
                outflow_node_nodes = (
                    outflow_node_nodes
                    .sort_values("order_no")
                    .drop_duplicates(subset="nodeID", keep=False)
                )

                if outflow_nodes_edges is None:
                    outflow_nodes_edges = outflow_node_edges.copy()
                else:
                    outflow_nodes_edges = pd.concat([outflow_nodes_edges, outflow_node_edges.copy()])
                if outflow_nodes_nodes is None:
                    outflow_nodes_nodes = outflow_node_nodes.copy()
                else:
                    outflow_nodes_nodes = pd.concat([outflow_nodes_nodes, outflow_node_nodes.copy()])
                
                new_outflow_nodes = self.nodes[self.nodes["nodeID"].isin(new_outflow_nodes)][["nodeID", "geometry"]]
                new_outflow_nodes["rws_code"] = outflow_node["rws_code"]
                new_outflow_nodes["rws_code_no"] = outflow_node["rws_code_no"]
                new_outflow_nodes["order_code"] = outflow_node["order_code"]
                new_outflow_nodes["order_no"] = order_no + 1

                if new_outflow_nodes_order is None:
                    new_outflow_nodes_order = new_outflow_nodes.copy()
                else:
                    new_outflow_nodes_order = pd.concat([new_outflow_nodes_order, new_outflow_nodes])

            edges_order = outflow_nodes_edges[outflow_nodes_edges["outflow_node"] >= 0]
            edges_left = edges_left[~edges_left.code.isin(edges_order.code)].copy()

            nodes_order = outflow_nodes_nodes[outflow_nodes_nodes["outflow_node"] >= 0]
            nodes_left = nodes_left[nodes_left.nodeID.isin(edges_left.node_start)].copy()

            if edges_all_orders is None:
                edges_all_orders = edges_order.copy()
            else:
                edges_all_orders = pd.concat([edges_all_orders, edges_order])
            if nodes_all_orders is None:
                nodes_all_orders = nodes_order.copy()
            else:
                nodes_all_orders = pd.concat([nodes_all_orders, nodes_order])
            
            if outflow_nodes_orders is None:
                outflow_nodes_orders = new_outflow_nodes_order.copy()
            else:
                outflow_nodes_orders = pd.concat([outflow_nodes_orders, new_outflow_nodes_order])
            
            order_no = order_no + 1
        
        edges_left["rws_code"] = ""
        edges_left["rws_code_no"] = -999
        edges_left["order_no"] = -999
        edges_left["outflow_node"] = -999
        edges_left["order_edge_no"] = -999

        nodes_left["rws_code"] = ""
        nodes_left["rws_code_no"] = -999
        nodes_left["order_no"] = -999
        nodes_left["outflow_node"] = -999
        nodes_left["order_node_no"] = -999
        
        edges_all_orders = edges_all_orders.sort_values("order_no").drop_duplicates(subset="code", keep="first")
        nodes_all_orders = nodes_all_orders.sort_values("order_no").drop_duplicates(subset="nodeID", keep="first")
        outflow_nodes_orders = outflow_nodes_orders.sort_values("order_no").drop_duplicates(subset="nodeID", keep="first")

        self.edges = pd.concat([edges_all_orders, edges_left]).reset_index(drop=True)
        self.nodes = pd.concat([nodes_all_orders, nodes_left]).reset_index(drop=True)
        self.outflow_nodes = outflow_nodes_orders.copy()

        len_edges_with_order = len(self.edges[self.edges.order_no > 0])
        len_edges_without_order = len(self.edges[self.edges.order_no < 0])
        logging.info(f"     - order levels generated for {len_edges_with_order} edges - {len_edges_without_order} left")


    def export_results_to_gpkg(self):
        """Export results to geopackages in folder 1_tussenresultaat"""
        results_dir = Path(self.path, self.dir_inter_results)
        logging.info(f"   x export results")
        for layer in [
            "dead_end_nodes",
            "dead_start_nodes",
            "outflow_nodes",
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
        outflow_nodes_4326 = self.outflow_nodes.to_crs(4326)

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

        if "order_no" in self.edges.columns:
            add_categorized_lines_to_map(
                m=m,
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

            if order_labels:
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
            label_fontsize=8, fg=fg
        )
        
        folium.GeoJson(
            self.outflow_nodes[self.outflow_nodes["order_no"]>2],
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
