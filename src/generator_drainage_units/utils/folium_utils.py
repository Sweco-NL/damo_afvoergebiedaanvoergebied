from shapely.geometry import Polygon, LineString, Point
import folium
import geopandas as gpd


def add_labels_to_points_lines_polygons(
        gdf: gpd.GeoDataFrame, column: str,
        label_fontsize: int = 14, label_unit: str = '',
        label_decimals: int = 2, show: bool = True,
        center=True, fg=None, fgs=None
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
            raise ValueError(' * GeoDataFrame does not have the right geometry')
        label_value = element[column]
        if (isinstance(label_value, float) or isinstance(label_value, int)) and label_decimals is not None:
            if label_unit == '%':
                label_str = f'{float(label_value):0.{label_decimals}%}'
            else:
                label_str = f'{float(label_value):0.{label_decimals}f}'
        else:
            label_str = f'{label_value}'
        html_style2 = f'<b>{label_str}{label_unit}</b></div>'
        if center:
            icon = folium.DivIcon(icon_size=(200, 50), icon_anchor=(-10, 15), html=html_style1 + html_style2)
        else:
            icon = folium.DivIcon(icon_size=(200, 50), icon_anchor=(-10, 20), html=html_style1 + html_style2)
        _label = folium.Marker(
            location=[point.y, point.x],
            icon=icon,
            show=show
        )
        if fgs is not None:
            _label.add_to(fgs)
        else:
            _label.add_to(fg)