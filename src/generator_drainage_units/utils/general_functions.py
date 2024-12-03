import geopandas as gpd
import pandas as pd
from shapely import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import split, linemerge
import numpy as np
from typing import TypeVar


def shorten_line_two_vertices(line, offset):
    # Ensure the geometry is a LineString
    if isinstance(line, LineString):
        # Calculate the length of the line
        length = line.length

        # Create new points for the shortened line
        start_point = line.interpolate(offset)
        end_point = line.interpolate(length - offset)

        # Return the new shortened line
        return LineString([start_point, end_point])
    else:
        return line  # Return the original line if not a LineString


def line_to_vertices(gdf, distance=10.0, keep_columns=["CODE", "WaterLineType"]):
    all_points = []

    for _, row in gdf.iterrows():
        geometry = row["geometry"]

        # Check if the geometry is a LineString or MultiLineString
        if isinstance(geometry, LineString):
            lines = [geometry]
        elif isinstance(geometry, MultiLineString):
            lines = list(geometry.geoms)
        else:
            continue  # Skip if not a recognized geometry typ

        for line in lines:
            # Add start and end points to the list with 'dangling' label
            start_point_x = {
                "geometry": Point(line.coords[0][0], line.coords[0][1]),
                "line_type": "dangling",
                "distance_from_start": 0,
            }
            end_point_x = {
                "geometry": Point(line.coords[-1][0], line.coords[-1][1]),
                "line_type": "dangling",
                "distance_from_start": line.length,
            }
            for point in [start_point_x, end_point_x]:
                for col in keep_columns:
                    point[col] = row[col]
                all_points.append(point)

            # Generate intermediate points at specified intervals
            num_points = int(line.length / distance)
            for i in range(1, num_points):
                new_distance = line.length / num_points
                # Calculate the point at a distance along the line
                point = line.interpolate(i * new_distance)
                point_2d = Point(point.x, point.y)
                # Calculate the cumulative distance from the starting point
                cumulative_distance = i * new_distance
                point = {
                    "geometry": point_2d,
                    "line_type": "other",
                    "distance_from_start": cumulative_distance,
                }
                for col in keep_columns:
                    point[col] = row[col]
                all_points.append(point)

    # Create a new GeoDataFrame for the vertices
    points_gdf = gpd.GeoDataFrame(all_points, crs=gdf.crs)

    return points_gdf


def split_waterways_by_endpoints(hydroobjects, relevant_culverts):
    end_points = relevant_culverts["geometry"].apply(
        lambda line: Point(line.coords[-1])
    )
    end_points_gdf = (
        relevant_culverts.copy()
    )  # Make a copy to retain all the other data

    # Replace the geometry column with the new end_points
    end_points_gdf.set_geometry(end_points, inplace=True)

    # Example of how to ensure your points are 2D
    end_points_gdf["geometry"] = end_points_gdf["geometry"].apply(
        lambda geom: Point(geom.x, geom.y)
        if geom.geom_type == "Point" and geom.has_z
        else geom
    )

    # Step 1: Convert MultiLineStrings to LineStrings
    hydroobjects["geometry"] = hydroobjects["geometry"].apply(
        lambda geom: [LineString(part.coords) for part in geom.geoms]
        if isinstance(geom, MultiLineString)
        else geom
    )

    # Since the previous step might convert some geometries to lists, we need to flatten them
    hydroobjects = hydroobjects.explode("geometry").reset_index(drop=True)

    # Step 2: Ensure that all LineStrings are 2D
    hydroobjects["geometry"] = hydroobjects["geometry"].apply(
        lambda geom: LineString([(x, y) for x, y, *_ in geom.coords])
        if geom.has_z
        else geom
    )

    hydroobjects = gpd.GeoDataFrame(hydroobjects, geometry='geometry', crs=28992)

    # Create lists to store start and end points
    start_points = []
    end_points = []

    # Iterate through the lines and extract start and end points
    for line in hydroobjects.geometry:
        start_points.append(Point(line.coords[0]))  # Start point
        end_points.append(Point(line.coords[-1]))

    # If you want to keep them as separate rows instead of columns, you can do this:
    start_points_gdf = gpd.GeoDataFrame(geometry=start_points, crs="EPSG:28992")
    end_points_gdf_ho = gpd.GeoDataFrame(geometry=end_points, crs="EPSG:28992")

    # Combine them into a single GeoDataFrame if desired
    combined_points_gdf = gpd.GeoDataFrame(
        pd.concat([start_points_gdf, end_points_gdf_ho], ignore_index=True),
        crs="EPSG:28992",
    )
    combined_end_points_gdf = gpd.GeoDataFrame(
        pd.concat([combined_points_gdf, end_points_gdf], ignore_index=True),
        crs="EPSG:28992",
    )
    combined_end_points_gdf = combined_end_points_gdf.drop_duplicates(subset='geometry')

    def split_lines_at_points(gdf_lines, gdf_points):

        #Create list of coordinates of linestring in gdf_lines
        gdf_lines['coords'] = gdf_lines.apply(lambda x: list(x.geometry.coords), axis=1)
        
        #Create segments with coords
        gdf_lines['segments'] = gdf_lines.apply(lambda x: [LineString(x['coords'][i-1 : i + 1]) for i in range(1, len(x['coords']))], axis=1)
        
        #Create gdf from segments of linestring
        gdf_segments = gpd.GeoDataFrame(
            gdf_lines[["CODE", "segments"]].explode("segments").rename(columns={'segments': "geometry"}), 
            geometry="geometry", 
            crs=gdf_lines.crs
        ).reset_index(drop=True)
        
        #Determine which segment point of points_gdf are located on
        gdf_segments = gdf_segments.merge(
            gdf_points[["dangling_id", 'geometry']].sjoin_nearest(
                gdf_segments.drop(columns="CODE")
            ).groupby('index_right').agg({'dangling_id': list, 'geometry': lambda x: list(x.apply(lambda geom: (geom.x, geom.y)))}).reset_index(), 
            how='outer', 
            left_index=True, 
            right_on="index_right",
            suffixes=("", "_2")
        ).reset_index(drop=True)
        
        def split_linestring(linestring, coordinates):
            # Convert coordinates to Points
            points = [Point(coord) for coord in coordinates]
            
            # Sort points by their position on the LineString
            points = sorted(points, key=lambda point: linestring.project(point))
            
            segments = []
            prev_point = Point(linestring.coords[0])
            
            for point in points:
                segment = LineString([prev_point, point])
                segments.append(segment)
                prev_point = point
            
            # Add the last segment
            segments.append(LineString([prev_point, Point(linestring.coords[-1])]))
            return segments

        def split_and_explode(row):
            if not row['geometry_2'] or isinstance(row['geometry_2'], float):
                return [row['geometry']]
            else:
                return split_linestring(row['geometry'], row['geometry_2'])

        gdf_segments["geometry"] = gdf_segments.apply(split_and_explode, axis=1)
        
        # Explode the list of geometries into separate rows
        gdf_segments = gdf_segments.explode("geometry").reset_index(drop=True)

        #Determine which segments are located between the same points of points_gdf
        gdf_segments["index_right_diff"] = gdf_segments["index_right"].diff(1)
        gdf_segments["index_right_cumsum"] = gdf_segments["index_right_diff"].map({0.0: 1, 1.0: 0})
        gdf_segments.loc[0, "index_right_cumsum"] = 0.0
        gdf_segments["index_right_cumsum"] = gdf_segments["index_right_cumsum"].cumsum().astype(int)
        
        #Group segments that are between points, creating a list of geometries in geometry
        gdf_segments = gdf_segments.groupby(['index_right_cumsum', 'CODE']).agg({
            "CODE": 'first', 
            "geometry": list, 
            'dangling_id': 'first',
        }).reset_index(drop=True)
        
        #Remove segments with a length of 0 (created points that are already between 2 linestrings of gdf_lines)
        gdf_segments["geometry_len"] = gdf_segments.apply(lambda x: sum([y.length for y in x["geometry"]]), axis=1)
        gdf_segments = gdf_segments[gdf_segments["geometry_len"]>0.0]

        #Merge the geometries of the segments that are listed in geometry
        gdf_segments["geometry"] = gdf_segments.apply(lambda x: linemerge(x["geometry"]), axis=1)

        #Create gdf from 
        gdf_segments = gpd.GeoDataFrame(gdf_segments.reset_index(drop=True), geometry='geometry', crs = gdf_lines.crs)

        # Function to add suffixes to duplicate CODEs
        def add_suffixes(df, column):
            counts = df[column].value_counts()
            suffixes = {code: 0 for code in counts.index if counts[code] > 1}
            
            def suffix_code(code):
                if code in suffixes:
                    suffix = suffixes[code]
                    suffixes[code] += 1
                    return f"{code}-{suffix}"
                return code
    
            df[column] = df[column].apply(suffix_code)
            return df
        # Apply the function to your GeoDataFrame
        gdf_segments = add_suffixes(gdf_segments, 'CODE')

        return gdf_segments

    result = split_lines_at_points(hydroobjects, combined_end_points_gdf)
    return result

def _remove_holes(geom, min_area):
    def p(p: Polygon, min_area) -> Polygon:
        holes = [i for i in p.interiors if not Polygon(i).area > min_area]
        return Polygon(shell=p.exterior, holes=holes)

    def mp(mp: MultiPolygon, min_area) -> MultiPolygon:
        return MultiPolygon([p(i, min_area) for i in mp.geoms])

    if isinstance(geom, Polygon):
        return p(geom, min_area)
    elif isinstance(geom, MultiPolygon):
        return mp(geom, min_area)
    else:
        return geom


_Geom = TypeVar("_Geom", Polygon, MultiPolygon, gpd.GeoSeries, gpd.GeoDataFrame)


def remove_holes_from_polygons(geom: _Geom, min_area: float) -> _Geom:
    """Remove all holes from a geometry that satisfy the filter function."""
    if isinstance(geom, gpd.GeoSeries):
        return geom.apply(_remove_holes, min_area=min_area)
    elif isinstance(geom, gpd.GeoDataFrame):
        geom = geom.copy()
        geom["geometry"] = remove_holes_from_polygons(
            geom["geometry"], min_area=min_area
        )
        return geom
    return _remove_holes(geom, min_area=min_area)


def remove_holes_from_basin_areas(basin_areas: gpd.GeoDataFrame, min_area: float):
    print(f" - remove holes within basin areas with less than {min_area/10000.0:.2f}ha")
    return remove_holes_from_polygons(geom=basin_areas, min_area=min_area)


def define_list_upstream_downstream_edges_ids(
    node_ids, nodes: gpd.GeoDataFrame, edges: gpd.GeoDataFrame
):
    nodes_sel = nodes[nodes.nodeID.isin(node_ids)].copy()
    for direction, node in zip(["upstream", "downstream"], ["node_end", "node_start"]):
        nodes_sel[f"{direction}_edges"] = nodes_sel.apply(
            lambda x: ",".join(list(edges[edges[node] == x.nodeID].code.values)), axis=1
        )
        nodes_sel[f"no_{direction}_edges"] = nodes_sel.apply(
            lambda x: len(x[f"{direction}_edges"].split(",")), axis=1
        )
    nodes_sel = nodes_sel.reset_index(drop=True)
    return nodes_sel
