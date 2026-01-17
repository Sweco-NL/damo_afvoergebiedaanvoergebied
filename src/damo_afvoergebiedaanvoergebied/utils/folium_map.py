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
    afvoergebied_cmap="Pastel2",
    afvoergebied_opacity=0.5,
    specifieke_afvoer_cmap="turbo",
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

    if not is_attribute_not_none(generator, "hydroobject"):
        raise ValueError("Generator does not have hydroobject")

    # Calculate the extent (bounding box) of your GeoDataFrame
    hydro_4326 = generator.hydroobject.to_crs(4326).copy()
    bounds = hydro_4326.total_bounds  # returns (minx, miny, maxx, maxy)
    
    # Center the map around the mean coordinates of the bounds
    center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]
    m = folium.Map(
        location=center,
        zoom_start=zoom_start,
        tiles=None,
    )
    logging.info('   x creating map layers')

    if is_attribute_not_none(generator, "rws_water"):
        logging.info('     - rws_water')
        folium.GeoJson(
            generator.rws_water.geometry,
            name="rws_water",
            z_index=0,
        ).add_to(m)

    if is_attribute_not_none(generator, "hydroobject"):
        logging.info('     - hydroobject')
        add_lines_to_map(
            m=m,
            lines_gdf=generator.hydroobject,
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
            name=f"outflow_nodesen in RWS-water", control=True
        ).add_to(m)

        logging.info('     - rws_water: outflow nodes')
        if "distance" in generator.outflow_nodes.columns:
            outflow_nodes = generator.outflow_nodes.dropna(subset="distance")
        else:
            outflow_nodes = generator.outflow_nodes.copy()
        folium.GeoJson(
            outflow_nodes,
            name="outflow_nodesen RWS-wateren",
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

    if is_attribute_not_none(generator, "overige_watergang"):
        logging.info('     - other waterways - without culverts')
        add_lines_to_map(
            m=m,
            lines_gdf=generator.overige_watergang[["geometry"]],
            layer_name=f"C-Watergangen - Zonder duikers ({len(generator.overige_watergang)})",
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

    if is_attribute_not_none(generator, f"outflow_nodes_overige_watergang"):
        logging.info(f'     - other waterways - outflow nodes')
        folium.GeoJson(
            generator.outflow_nodes_overige_watergang.geometry,
            name="C-Watergangen - outflow_nodesen",
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

    if is_attribute_not_none(generator, f"overige_watergang_processed_4"):
        overige_watergang_processed = generator.overige_watergang_processed_4.copy()
    elif is_attribute_not_none(generator, f"overige_watergang_processed_3"):
        overige_watergang_processed = generator.overige_watergang_processed_3.copy()
    else:
        overige_watergang_processed = None

    if overige_watergang_processed is not None:
        logging.info(f'     - other waterways - aggregated per outflow node')
        add_categorized_lines_to_map(
            m=m,
            lines_gdf=overige_watergang_processed,
            layer_name=f"C-Watergangen - Gegroepeerd per outflow_nodes",
            control=True,
            lines=True,
            line_color_column="outflow_nodes",
            line_color_cmap=None,
            show=False,
            z_index=2,
        )
        
        if order_labels and "order_code" in overige_watergang_processed.columns:
            logging.info('     - other waterways - order code (labels)')
            fg = folium.FeatureGroup(
                name=f"C-watergangen - Orde-code (labels)",
                control=True,
                show=False,
                z_index=2,
            ).add_to(m)

            add_labels_to_points_lines_polygons(
                gdf=overige_watergang_processed[["geometry", "order_code"]],
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

    if is_attribute_not_none(generator, f"afvoergebied_0"):
        logging.info(f'     - drainage units (raster, level 0)')
        afvoergebied_0 = generator.afvoergebied_0.where(generator.afvoergebied_0 > -1.0)
        afvoergebied_0 = afvoergebied_0.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=afvoergebied_0,
            layer_name="Afwateringseenheden - basis",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.afvoergebied_0.data.max()),
            cmap=afvoergebied_cmap,
            legend=False,
            opacity=afvoergebied_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"afvoergebied_1"):
        logging.info(f'     - drainage units (raster, level 1)')
        afvoergebied_1 = generator.afvoergebied_1.where(generator.afvoergebied_1 > -1.0)
        afvoergebied_1 = afvoergebied_1.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=afvoergebied_1,
            layer_name="Afwateringseenheden (A/B/C watergangen)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.afvoergebied_1.data.max()),
            cmap=afvoergebied_cmap,
            legend=False,
            opacity=afvoergebied_opacity,
            show=True,
            dx=dx,
            dy=dy,
        )
    
    if is_attribute_not_none(generator, f"afvoergebied_2"):
        logging.info(f'     - drainage units (raster, level 2)')
        afvoergebied_2 = generator.afvoergebied_2.where(generator.afvoergebied_2 > -1.0)
        afvoergebied_2 = afvoergebied_2.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=afvoergebied_2,
            layer_name="Afwateringseenheden (A/B watergangen)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.afvoergebied_2.data.max()),
            cmap=afvoergebied_cmap,
            legend=False,
            opacity=afvoergebied_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"afvoergebied_3"):
        logging.info(f'     - drainage units (raster, level 3)')
        afvoergebied_3 = generator.afvoergebied_3.where(generator.afvoergebied_3 > -1.0)
        afvoergebied_3 = afvoergebied_3.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=afvoergebied_3,
            layer_name="Afwateringseenheden (orde-code)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.afvoergebied_3.data.max()),
            cmap=afvoergebied_cmap,
            legend=False,
            opacity=afvoergebied_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, f"afvoergebied_4"):
        logging.info(f'     - drainage units (raster, level 4)')
        afvoergebied_4 = generator.afvoergebied_4.where(generator.afvoergebied_4 > -1.0)
        afvoergebied_4 = afvoergebied_4.rio.write_crs(generator.crs)
        add_graduated_raster_to_map(
            m=m,
            raster=afvoergebied_4,
            layer_name="Afwateringseenheden (stroomgebied)",
            unit="unique id",
            control=True,
            vmin=0,
            vmax=int(generator.afvoergebied_4.data.max()),
            cmap=afvoergebied_cmap,
            legend=False,
            opacity=afvoergebied_opacity,
            show=False,
            dx=dx,
            dy=dy,
        )

    if is_attribute_not_none(generator, "edges"):
        if "total_specifieke_afvoer" in generator.edges.columns:
            edges = generator.edges[
                generator.edges["total_specifieke_afvoer"] > 0.0
            ].sort_values("total_specifieke_afvoer", ascending=True)
            
            logging.info('     - edges - total specific discharge')
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges[["geometry", "total_specifieke_afvoer"]],
                layer_name="AB-watergangen - Total specific discharge",
                control=True,
                lines=True,
                line_color_column="total_specifieke_afvoer",
                line_color_cmap=specifieke_afvoer_cmap,
                label=False,
                line_weight=3,
                z_index=2,
                show=False,
            )

            logging.info('     - edges - total specific discharge (log10)')
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=edges[["geometry", "log10_total_specifieke_afvoer"]],
                layer_name="AB-watergangen - Total specific discharge (log10)",
                control=True,
                lines=True,
                line_color_column="log10_total_specifieke_afvoer",
                line_color_cmap=specifieke_afvoer_cmap,
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
                gdf=edges[["geometry", "total_specifieke_afvoer"]],
                column="total_specifieke_afvoer",
                label_fontsize=8,
                label_decimals=1,
                fg=fg,
            )

            edges_left = generator.edges[
                generator.edges["total_specifieke_afvoer"] <= 0.0
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
        
    if is_attribute_not_none(generator, "split_nodes"):
        logging.info('     - split nodes')
        folium.GeoJson(
            generator.split_nodes,
            name=f"Splitsingspunten benedenstrooms ({len(generator.split_nodes)})",
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