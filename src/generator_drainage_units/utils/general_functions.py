import geopandas as gpd
import pandas as pd
from shapely import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import split
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

    def split_lines_at_points(gdf_lines, gdf_points):
        new_lines = []

        # Iterate over each line in the GeoDataFrame
        for index, line in gdf_lines.iterrows():
            # Convert the LineString to a series of coordinates
            coords = list(line.geometry.coords)

            # Prepare to track segments
            segments = []
            segment_start = 0

            # Iterate through the coordinates to create segments
            for i in range(1, len(coords)):
                segment = LineString(
                    coords[segment_start : i + 1]
                )  # Create segment from start to current point
                segments.append(segment)
                segment_start = i  # Move the start to the current point

            # Now check for points on each segment
            for segment in segments:
                for _, point in gdf_points.iterrows():
                    if (
                        segment.distance(point.geometry) < 1e-8
                    ):  # Check if the point is very close to the segment
                        # Add point's coordinate twice
                        new_coords = list(segment.coords)
                        new_coords.append(
                            (point.geometry.x, point.geometry.y)
                        )  # Add point
                        new_coords.append(
                            (point.geometry.x, point.geometry.y)
                        )  # Add point again

                        # Create new segments including the added point
                        # Sort coordinates so points are in correct order
                        new_coords.sort(key=lambda coord: segment.project(Point(coord)))

                        # Split the segment at the double point
                        for j in range(len(new_coords) - 1):
                            new_segment = LineString(
                                new_coords[j : j + 2]
                            )  # Create a new segment
                            new_lines.append(new_segment)

                        # Adjust to start new segment from next coordinate
                        break  # Exit after adding point to avoid duplicate processing
                else:  # This else runs if no break occurs
                    new_lines.append(
                        segment
                    )  # Append the original segment if no point was added

        # Convert the new lines back into a GeoDataFrame
        gdf_new_lines = gpd.GeoDataFrame(geometry=new_lines)
        return gdf_new_lines

    segments = split_lines_at_points(hydroobjects, end_points_gdf)

    def merge_segments(segments_gdf, full_endpoints_gdf):
        merged_lines = []
        current_line_coords = []

        # Get all the end points for easy checking
        endpoint_coords = set(
            zip(full_endpoints_gdf.geometry.x, full_endpoints_gdf.geometry.y)
        )

        for index, row in segments_gdf.iterrows():
            segment = row.geometry
            segment_start = segment.coords[0]
            segment_end = segment.coords[-1]

            # If we are starting a new segment, reset the current line
            if len(current_line_coords) == 0:
                current_line_coords.append(segment_start)

            # Append current segment coordinates
            current_line_coords.extend(
                list(segment.coords)[1:]
            )  # Avoid duplicating the starting point

            # Check if the segment ends at an endpoint from full_endpoints_gdf
            if segment_end in endpoint_coords:
                current_line_coords.append(
                    segment_end
                )  # Ensure the endpoint is included

                # Create a LineString and add to merged lines
                merged_lines.append(LineString(current_line_coords))
                current_line_coords = []  # Reset for the next set
            elif segment_start in endpoint_coords:
                # If it starts at an endpoint, we need to finalize before that
                if current_line_coords:
                    current_line_coords.append(segment_start)
                    merged_lines.append(LineString(current_line_coords))
                    current_line_coords = []

        # If there's any remaining line not finalized, add it
        if current_line_coords:
            merged_lines.append(LineString(current_line_coords))

        # Create a GeoDataFrame from the merged segments
        merged_gdf = gpd.GeoDataFrame(geometry=merged_lines, crs=segments_gdf.crs)
        return merged_gdf

    result = merge_segments(segments, combined_end_points_gdf)

    def clean_doubled_lines(gdf_segments):
        cleaned_segments = []

        for i in range(len(gdf_segments)):
            segment = gdf_segments.geometry.iloc[i]
            coords = list(segment.coords)

            # Check if the segment has more than one point
            if len(coords) > 1:
                # If the first and last coordinates are the same
                if coords[0] == coords[-1]:
                    # Check if there are more than two points (to have an in-between point)
                    if len(coords) > 2:
                        # Keep the first point and the in-between point
                        new_segment = LineString([coords[0], coords[1]])
                        cleaned_segments.append(new_segment)
                    # If there are only two points, we simply skip this segment (remove it)
                else:
                    # If the segment does not have the same start and end, keep it
                    cleaned_segments.append(segment)

        # Create a GeoDataFrame from the cleaned segments
        cleaned_gdf = gpd.GeoDataFrame(geometry=cleaned_segments)
        return cleaned_gdf

    cleaned = clean_doubled_lines(result)

    def merge_lines_based_on_endpoints(segments_gdf, end_points_gdf):
        merged_lines = []

        # Create a set of all endpoint coordinates in end_points_gdf for quick lookup
        endpoint_coords = set(zip(end_points_gdf.geometry.x, end_points_gdf.geometry.y))

        # Iterate over each segment in the GeoDataFrame
        for index, line in segments_gdf.iterrows():
            current_line = line.geometry

            # Check if the endpoint of the current line is in the endpoint coordinates
            endpoint = (
                current_line.coords[-1][0],
                current_line.coords[-1][1],
            )  # Get the endpoint as a tuple

            if endpoint not in endpoint_coords:
                # Find the next line with the same start point as the current line's endpoint
                for next_index, next_line in segments_gdf.iterrows():
                    if (
                        line.geometry.coords[-1] == next_line.geometry.coords[0]
                    ):  # Check if the start of the next line matches
                        # Merge the two lines
                        current_line = gpd.GeoSeries(
                            [current_line, next_line.geometry]
                        ).union_all()
                        break  # Exit the inner loop after merging

            # Append the (possibly merged) line to the list
            merged_lines.append(current_line)

        # Create a new GeoDataFrame from the merged lines
        merged_gdf = gpd.GeoDataFrame(geometry=merged_lines)

        # Remove duplicate geometries
        merged_gdf = merged_gdf[~merged_gdf.geometry.duplicated()]

        return merged_gdf

    final_result = merge_lines_based_on_endpoints(cleaned, combined_end_points_gdf)

    def remove_covered_lines(merged_gdf):
        # Create a GeoDataFrame to store the lines that are not covered
        filtered_gdf = merged_gdf.copy()

        # Iterate through each line in the GeoDataFrame
        for idx, line in merged_gdf.iterrows():
            # Check if this line is covered by any other lines
            if any(
                filtered_gdf.geometry.notna()
                & (filtered_gdf.geometry != line.geometry)
                & line.geometry.within(filtered_gdf.geometry)
            ):
                # Remove the covered line
                filtered_gdf = filtered_gdf[filtered_gdf.geometry != line.geometry]

        return filtered_gdf

    final_results_cleaned = remove_covered_lines(final_result)

    return final_results_cleaned


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
    return nodes_sel.reset_index(drop=True)
