import folium
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point, Polygon


def add_labels_to_points_lines_polygons(
    gdf: gpd.GeoDataFrame,
    column: str,
    label_fontsize: int = 14,
    label_unit: str = "",
    label_decimals: int = 2,
    show: bool = True,
    center=True,
    fg=None,
    fgs=None,
):
    gdf = gdf.to_crs(4326).copy()

    if column not in gdf.columns:
        return

    for element_id, element in gdf.iterrows():
        if center:
            html_style1 = f'<div style="font-size: {label_fontsize}pt; color: black">'
        else:
            html_style1 = f'<div style="font-size: {label_fontsize}pt; color: black">'
        if isinstance(element.geometry, Polygon):
            point = element.geometry.representative_point()
        elif isinstance(element.geometry, LineString):
            point = element.geometry.interpolate(0.5, normalized=True)
        elif isinstance(element.geometry, Point):
            point = element.geometry
        else:
            raise ValueError(" * GeoDataFrame does not have the right geometry")

        label_value = element[column]
        if isinstance(label_value, float) and np.isnan(label_value):
            return

        if (
            isinstance(label_value, float) or isinstance(label_value, int)
        ) and label_decimals is not None:
            if label_unit == "%":
                label_str = f"{float(label_value):0.{label_decimals}%}"
            else:
                label_str = f"{float(label_value):0.{label_decimals}f}"
        else:
            label_str = f"{label_value}"
        html_style2 = f"<b>{label_str}{label_unit}</b></div>"
        if center:
            icon = folium.DivIcon(
                icon_size=(200, 50),
                icon_anchor=(-10, 15),
                html=html_style1 + html_style2,
            )
        else:
            icon = folium.DivIcon(
                icon_size=(200, 50),
                icon_anchor=(-10, 20),
                html=html_style1 + html_style2,
            )
        _label = folium.Marker(location=[point.y, point.x], icon=icon, show=show)
        if fgs is not None:
            _label.add_to(fgs)
        else:
            _label.add_to(fg)


def add_basemaps_to_folium_map(m: folium.Map, base_map="OpenStreetMap"):
    m.tiles = None
    basemaps = ["ESRI Luchtfoto", "Dark Mode", "Light Mode", "OpenStreetMap"]
    basemap_types = [
        {
            "tiles": "cartodbpositron",
            "name": "Light Mode",
            "attr": None,
            "control": True,
            "maxNativeZoom": 20,
            "maxZoom": 20,
        },
        {
            "tiles": "openstreetmap",
            "name": "OpenStreetMap",
            "attr": None,
            "control": True,
            "maxNativeZoom": 19,
            "maxZoom": 19,
            "show": True,
        },
        {
            "tiles": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
            "attr": "Esri",
            "name": "ESRI Luchtfoto",
            "control": True,
            "maxNativeZoom": 21,
            "maxZoom": 21,
            "show": True,
        },
        {
            "tiles": "cartodbdark_matter",
            "name": "Dark Mode",
            "attr": None,
            "control": True,
            "maxNativeZoom": 20,
            "maxZoom": 20,
            "show": True,
        },
        {
            "tiles": "Stamen Toner",
            "name": "Stamen Toner",
            "attr": None,
            "control": True,
            "maxNativeZoom": 17,
            "maxZoom": 17,
            "show": True,
        },
    ]

    for bm in basemaps:
        basemap = [o for o in basemap_types if o["name"] == bm][0]
        folium.TileLayer(
            tiles=basemap["tiles"],
            name=basemap["name"],
            attr=basemap["attr"],
            control=basemap["control"],
            maxNativeZoom=basemap["maxNativeZoom"],
            maxZoom=basemap["maxZoom"],
            show=True if basemap["name"] == base_map else False,
        ).add_to(m)
    return m
