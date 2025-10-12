import logging
import webbrowser
from pathlib import Path
import folium
from .folium_utils import (
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
    add_labels_to_points_lines_polygons,
    add_graduated_raster_to_map,
    add_lines_to_map,
)
logging.basicConfig(level=logging.INFO)


def generate_folium_map(
    generator, 
    all_culverts=False,
    order_labels=True,
    drop_duplicate_codes=True,
    show_other_waterways_culverts=False,
    drainage_units_cmap="Pastel2",
    drainage_units_opacity=0.5,
    specific_discharge_cmap="turbo",
    html_file_name=None, 
    base_map="Light Mode", 
    save_html=True,
    open_html=False, 
    zoom_start=11,
    zmin=None,
    zmax=None,
    dx=0.0,
    dy=0.0,
):
    def is_attribute_not_none(obj, attribute):
        return hasattr(obj, attribute) and getattr(obj, attribute) is not None

    if not is_attribute_not_none(generator, "hydroobjecten"):
        raise ValueError("Generator does not have hydroobjecten")

    # Calculate the extent (bounding box) of your GeoDataFrame
    hydro_4326 = generator.hydroobjecten.to_crs(4326).copy()
    bounds = hydro_4326.total_bounds  # returns (minx, miny, maxx, maxy)
    
    # Center the map around the mean coordinates of the bounds
    center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]
    m = folium.Map(
        location=center,
        zoom_start=zoom_start,
        tiles=None,
    )
    logging.info('   x creating map layers')

    if is_attribute_not_none(generator, "rws_wateren"):
        logging.info('     - rws_wateren')
        folium.GeoJson(
            generator.rws_wateren.geometry,
            name="RWS_Wateren",
            z_index=0,
        ).add_to(m)

    if is_attribute_not_none(generator, "hydroobjecten"):
        logging.info('     - hydroobjecten')
        add_lines_to_map(
            m=m,
            lines_gdf=generator.hydroobjecten,
            layer_name="AB-Watergangen",
            control=True,
            lines=True,
            label=False,
            line_weight=2.5,
            z_index=1,
        )

    # if is_attribute_not_none(generator, "nodes"):
    #     folium.GeoJson(
    #         generator.nodes,
    #         name="nodes",
    #         marker=folium.Circle(
    #             radius=2.5, fill_color="black", fill_opacity=0.5, color="black", width=1,
    #         ),
    #         highlight_function=lambda x: {"fillOpacity": 0.8},
    #         z_index=1,
    #         control=False,
    #     ).add_to(m)

    if is_attribute_not_none(generator, "edges"):
        if "order_no" in generator.edges.columns:
            if "order_edge_no" in generator.edges.columns:
                edges = generator.edges[
                    generator.edges["order_no"] > 1
                ].sort_values(["order_no", "order_edge_no"], ascending=[False, True])
            else:
                edges = generator.edges[
                    generator.edges["order_no"] > 1
                ].sort_values(["order_no"], ascending=[False])
            edges_left = generator.edges[generator.edges["order_no"] < 2].copy()
            edges_labels = edges.copy()

            if "order_code" in edges_labels.columns and drop_duplicate_codes:
                edges_labels = edges_labels.drop_duplicates(subset="order_code", keep="first")

            logging.info('     - edges - order_no')
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges[["geometry", "order_no"]],
                layer_name=f"AB-watergangen - Orde-nummer ({len(edges)})",
                control=True,
                lines=True,
                line_color_column="order_no",
                line_color_cmap="hsv_r",
                label=False,
                line_weight=3,
                z_index=2,
                show=True,
            )

            logging.info('     - edges - no order_no')
            add_lines_to_map(
                m=m,
                lines_gdf=edges_left,
                layer_name=f"AB-Watergangen - Geen Orde-nummer ({len(edges_left)})",
                line_color="black",
                control=True,
                lines=True,
                label=False,
                show=False,
                line_weight=3,
                z_index=0,
            )

            if order_labels and "order_no" in generator.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"AB-watergangen - Orde-nummer (labels)",
                    control=True,
                    show=False,
                    z_index=2,
                ).add_to(m)

                logging.info('     - edges - order_no (labels)')
                add_labels_to_points_lines_polygons(
                    gdf=edges_labels,
                    column="order_no",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )

            if order_labels and "order_code" in generator.edges.columns:
                fg = folium.FeatureGroup(
                    name=f"AB-watergangen - Orde-code (labels)",
                    control=True,
                    show=False,
                    z_index=2,
                ).add_to(m)

                logging.info('     - edges - order_code (labels)')
                add_labels_to_points_lines_polygons(
                    gdf=edges_labels,
                    column="order_code",
                    label_fontsize=8,
                    center=True,
                    fg=fg,
                )
    
    if is_attribute_not_none(generator, "outflow_nodes"):
        fg = folium.FeatureGroup(
            name=f"Uitstroompunten in RWS-water", control=True
        ).add_to(m)

        logging.info('     - rws_water: outflow nodes')
        if "distance" in generator.outflow_nodes.columns:
            outflow_nodes = generator.outflow_nodes.dropna(subset="distance")
        else:
            outflow_nodes = generator.outflow_nodes.copy()
        folium.GeoJson(
            outflow_nodes,
            name="Uitstroompunten RWS-wateren",
            marker=folium.Circle(
                radius=25, fill_color="red", fill_opacity=0.4, color="red", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            zoom_on_click=True,
            z_index=3,
        ).add_to(fg)

        add_labels_to_points_lines_polygons(
            gdf=outflow_nodes,
            column="order_code",
            label_fontsize=8,
            fg=fg,
        )

    if is_attribute_not_none(generator, "overige_watergangen"):
        logging.info('     - other waterways - without culverts')
        add_lines_to_map(
            m=m,
            lines_gdf=generator.overige_watergangen[["geometry"]],
            layer_name=f"C-Watergangen - Zonder duikers ({len(generator.overige_watergangen)})",
            line_color="lightblue",
            control=True,
            lines=True,
            label=False,
            line_weight=2.5,
            z_index=0,
            show=show_other_waterways_culverts,
        )

    for i in range(5,0,-1):
        if is_attribute_not_none(generator, f"potential_culverts_{i}"):
            logging.info(f'     - other waterways - potential culverts ({i})')
            add_lines_to_map(
                m=m,
                lines_gdf=getattr(generator, f"potential_culverts_{i}")[["geometry"]],
                layer_name=f"C-Watergangen - Gevonden Duikers ({i})",
                line_color="red",
                control=True,
                lines=True,
                label=False,
                line_weight=2.5,
                z_index=1,
                show=show_other_waterways_culverts,
            )
            if not all_culverts:
                break
            show_other_waterways_culverts = False

    if is_attribute_not_none(generator, f"outflow_nodes_overige_watergangen"):
        logging.info(f'     - other waterways - outflow nodes')
        folium.GeoJson(
            generator.outflow_nodes_overige_watergangen.geometry,
            name="C-Watergangen - Uitstroompunten",
            marker=folium.Circle(
                radius=3,
                fill_color="orange",
                fill_opacity=1.0,
                color="orange",
                weight=1,
                z_index=3,
            ),
            show=False,
        ).add_to(m)

    if is_attribute_not_none(generator, f"overige_watergangen_processed_4"):
        overige_watergangen_processed = generator.overige_watergangen_processed_4.copy()
    elif is_attribute_not_none(generator, f"overige_watergangen_processed_3"):
        overige_watergangen_processed = generator.overige_watergangen_processed_3.copy()
    else:
        overige_watergangen_processed = None

    if overige_watergangen_processed is not None:
        logging.info(f'     - other waterways - aggregated per outflow node')
        add_categorized_lines_to_map(
            m=m,
            lines_gdf=overige_watergangen_processed,
            layer_name=f"C-Watergangen - Gegroepeerd per uitstroompunt",
            control=True,
            lines=True,
            line_color_column="outflow_node",
            line_color_cmap=None,
            show=False,
            z_index=2,
        )
        
        if order_labels and "order_code" in overige_watergangen_processed.columns:
            logging.info('     - other waterways - order code (labels)')
            fg = folium.FeatureGroup(
                name=f"C-watergangen - Orde-code (labels)",
                control=True,
                show=False,
                z_index=2,
            ).add_to(m)

            add_labels_to_points_lines_polygons(
                gdf=overige_watergangen_processed[["geometry", "order_code"]],
                column="order_code",
                label_fontsize=8,
                label_decimals=1,
                fg=fg,
            )

    if is_attribute_not_none(generator, f"ghg"):
        logging.info(f'     - raster topography or groundwaterlevel')
        add_graduated_raster_to_map(
            m=m,
            raster=generator.ghg,
            layer_name="GHG",
            unit="m NAP",
            control=True,
            vmin=generator.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
            vmax=generator.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
            legend=False,
            opacity=1.00,
            show=False,
            cmap='Spectral_r',
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"ghg_processed"):
        logging.info(f'     - raster topography or groundwaterlevel (processed)')
        add_graduated_raster_to_map(
            m=m,
            raster=generator.ghg_processed,
            layer_name="GHG-processed",
            unit="m NAP",
            control=True,
            vmin=generator.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
            vmax=generator.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
            legend=False,
            opacity=1.00,
            show=False,
            cmap='Spectral_r',
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"ghg_processed_adapt"):
        logging.info(f'     - raster topography or groundwaterlevel (adapted)')
        add_graduated_raster_to_map(
            m=m,
            raster=generator.ghg_processed_adapt,
            layer_name="GHG-processed-adapt",
            unit="m NAP",
            control=True,
            vmin=generator.ghg.STATISTICS_MINIMUM if zmin is None else zmin,
            vmax=generator.ghg.STATISTICS_MAXIMUM if zmax is None else zmax,
            legend=False,
            opacity=1.00,
            show=False,
            cmap='Spectral_r',
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"drainage_units_0"):
        logging.info(f'     - drainage units (raster, level 0)')
        drainage_units_0 = generator.drainage_units_0.where(generator.drainage_units_0 > -1.0)
        drainage_units_0 = drainage_units_0.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=drainage_units_0,
            layer_name="Afwateringseenheden (A/B/C watergangen)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.drainage_units_0.data.max()),
            cmap=drainage_units_cmap,
            legend=False,
            opacity=drainage_units_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"drainage_units_1"):
        logging.info(f'     - drainage units (raster, level 1)')
        drainage_units_1 = generator.drainage_units_1.where(generator.drainage_units_1 > -1.0)
        drainage_units_1 = drainage_units_1.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=drainage_units_1,
            layer_name="Afwateringseenheden (A/B watergangen)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.drainage_units_1.data.max()),
            cmap=drainage_units_cmap,
            legend=False,
            opacity=drainage_units_opacity,
            show=True,
            dx=dx,
            dy=dy,
        )
    
    if is_attribute_not_none(generator, f"drainage_units_2"):
        logging.info(f'     - drainage units (raster, level 2)')
        drainage_units_2 = generator.drainage_units_2.where(generator.drainage_units_2 > -1.0)
        drainage_units_2 = drainage_units_2.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=drainage_units_2,
            layer_name="Afwateringseenheden (orde-code)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.drainage_units_2.data.max()),
            cmap=drainage_units_cmap,
            legend=False,
            opacity=drainage_units_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"drainage_units_3"):
        logging.info(f'     - drainage units (raster, level 3)')
        drainage_units_3 = generator.drainage_units_3.where(generator.drainage_units_3 > -1.0)
        drainage_units_3 = drainage_units_3.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=drainage_units_3,
            layer_name="Afwateringseenheden (stroomgebied)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.drainage_units_3.data.max()),
            cmap=drainage_units_cmap,
            legend=False,
            opacity=drainage_units_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, "edges"):
        if "total_specific_discharge" in generator.edges.columns:
            edges = generator.edges[
                generator.edges["total_specific_discharge"] > 0.0
            ].sort_values("total_specific_discharge", ascending=True)
            
            logging.info('     - edges - total specific discharge')
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges[["geometry", "total_specific_discharge"]],
                layer_name="AB-watergangen - Total specific discharge",
                control=True,
                lines=True,
                line_color_column="total_specific_discharge",
                line_color_cmap=specific_discharge_cmap,
                label=False,
                line_weight=3,
                z_index=2,
                show=False,
            )

            logging.info('     - edges - total specific discharge (log10)')
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges[["geometry", "log10_total_specific_discharge"]],
                layer_name="AB-watergangen - Total specific discharge (log10)",
                control=True,
                lines=True,
                line_color_column="log10_total_specific_discharge",
                line_color_cmap=specific_discharge_cmap,
                label=False,
                line_weight=3,
                z_index=2,
                show=True,
            )

            logging.info('     - edges - total specific discharge (labels)')
            fg = folium.FeatureGroup(
                name=f"AB-watergangen - Total specific discharge (labels)",
                control=True,
                show=False,
                z_index=2,
            ).add_to(m)

            add_labels_to_points_lines_polygons(
                gdf=edges[["geometry", "total_specific_discharge"]],
                column="total_specific_discharge",
                label_fontsize=8,
                label_decimals=1,
                fg=fg,
            )

            edges_left = generator.edges[
                generator.edges["total_specific_discharge"] <= 0.0
            ]
            add_lines_to_map(
                m=m,
                lines_gdf=edges_left[["geometry"]],
                layer_name=f"AB-Watergangen - Geen specifieke afvoer ({len(edges_left)})",
                line_color="black",
                control=True,
                lines=True,
                label=False,
                show=False,
                line_weight=2.5,
                z_index=0,
            )
        
    if is_attribute_not_none(generator, "nodes"):
        logging.info('     - split nodes')
        split_points_downstreams = generator.nodes[
            (generator.nodes["no_downstream_edges"] > 1) &
            (generator.nodes["no_upstream_edges"] > 0)
        ].copy()
        folium.GeoJson(
            split_points_downstreams,
            name=f"Splitsingspunten benedenstrooms ({len(split_points_downstreams)})",
            marker=folium.Circle(
                radius=25, fill_color="purple", fill_opacity=0.4, color="purple", weight=3
            ),
            highlight_function=lambda x: {"fillOpacity": 0.8},
            z_index=3,
            show=False,
        ).add_to(m)

    m = add_basemaps_to_folium_map(m=m, base_map=base_map)

    folium.LayerControl(collapsed=False).add_to(m)
    m.add_child(folium.plugins.MeasureControl())

    generator.folium_map = m
    if save_html:
        if html_file_name is None:
            html_file_name = generator.name

        generator.folium_html_path = Path(generator.path, f"{html_file_name}.html")
        m.save(generator.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")
        if open_html:
            webbrowser.open(Path(generator.path, f"{html_file_name}.html"))
    return m