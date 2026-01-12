import geopandas as gpd
import pandas as pd
from shapely.geometry import (
    Point,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
    base,
)
from shapely.ops import split, linemerge, nearest_points
import numpy as np
from typing import TypeVar
import time
import logging
import warnings
import shapely


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


def line_to_vertices(gdf, distance=10.0, keep_columns=["code", "WaterLineType"]):
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
                "line_type": "dangling_start",
                "distance_from_start": 0,
            }
            end_point_x = {
                "geometry": Point(line.coords[-1][0], line.coords[-1][1]),
                "line_type": "dangling_end",
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
    end_points_gdf = relevant_culverts.copy()
    end_points_gdf.set_geometry(end_points, inplace=True)
    end_points_gdf = remove_z_dims(end_points_gdf)

    # Step 1: Convert MultiLineStrings to LineStrings
    hydroobjects = hydroobjects.explode("geometry").reset_index(drop=True)
    hydroobjects = remove_z_dims(hydroobjects)
    # hydroobjects = gpd.GeoDataFrame(hydroobjects, geometry="geometry", crs=28992)

    # Create 'NAAM' and 'status' if not present
    if "status" not in hydroobjects.columns:
        hydroobjects["status"] = "unchanged"

    if "NAAM" not in hydroobjects.columns:
        hydroobjects["NAAM"] = "unknown"

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
    combined_end_points_gdf = combined_end_points_gdf.drop_duplicates(subset="geometry")

    def split_lines_at_points(gdf_lines, gdf_points):
        # Create list of coordinates of linestring in gdf_lines
        gdf_lines["coords"] = gdf_lines.apply(lambda x: list(x.geometry.coords), axis=1)

        # Create segments with coords
        gdf_lines["segments"] = gdf_lines.apply(
            lambda x: [
                LineString(x["coords"][i - 1 : i + 1])
                for i in range(1, len(x["coords"]))
            ],
            axis=1,
        )

        # Create gdf from segments of linestring
        gdf_segments = gpd.GeoDataFrame(
            gdf_lines[["code", "segments", "NAAM", "status"]]
            .explode("segments")
            .rename(columns={"segments": "geometry"}),
            geometry="geometry",
            crs=gdf_lines.crs,
        ).reset_index(drop=True)

        # Determine which segment point of points_gdf are located on
        gdf_segments = gdf_segments.merge(
            gdf_points[["geometry"]]
            .sjoin_nearest(gdf_segments.drop(columns="code"))
            .groupby("index_right")
            .agg(
                {
                    "geometry": lambda x: list(x.apply(lambda geom: (geom.x, geom.y))),
                }
            )
            .reset_index(),
            how="outer",
            left_index=True,
            right_on="index_right",
            suffixes=("", "_2"),
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
            if not row["geometry_2"] or isinstance(row["geometry_2"], float):
                return [row["geometry"]]
            else:
                return split_linestring(row["geometry"], row["geometry_2"])

        gdf_segments["geometry"] = gdf_segments.apply(split_and_explode, axis=1)

        # Explode the list of geometries into separate rows
        gdf_segments = gdf_segments.explode("geometry").reset_index(drop=True)

        # Determine which segments are located between the same points of points_gdf
        gdf_segments["index_right_diff"] = gdf_segments["index_right"].diff(1)
        gdf_segments["index_right_cumsum"] = gdf_segments["index_right_diff"].map(
            {0.0: 1, 1.0: 0}
        )
        gdf_segments.loc[0, "index_right_cumsum"] = 0.0
        gdf_segments["index_right_cumsum"] = (
            gdf_segments["index_right_cumsum"].cumsum().astype(int)
        )

        # Group segments that are between points, creating a list of geometries in geometry
        gdf_segments = (
            gdf_segments.groupby(["index_right_cumsum", "code"])
            .agg(
                {
                    "code": "first",
                    "geometry": list,
                    "NAAM": "first",
                    "status": "first",
                }
            )
            .reset_index(drop=True)
        )

        # Remove segments with a length of 0 (created points that are already between 2 linestrings of gdf_lines)
        gdf_segments["geometry_len"] = gdf_segments.apply(
            lambda x: sum([y.length for y in x["geometry"]]), axis=1
        )
        gdf_segments = gdf_segments[gdf_segments["geometry_len"] > 0.0]

        # Merge the geometries of the segments that are listed in geometry
        gdf_segments["geometry"] = gdf_segments.apply(
            lambda x: linemerge(x["geometry"]), axis=1
        )

        # Create gdf from
        gdf_segments = gpd.GeoDataFrame(
            gdf_segments.reset_index(drop=True), geometry="geometry", crs=gdf_lines.crs
        )

        # Function to add suffixes to duplicate codes
        def add_suffixes(df, column):
            counts = df[column].value_counts()
            suffixes = {code: 0 for code in counts.index if counts[code] > 1}

            def suffix_code(code):
                if code in suffixes:
                    suffix = suffixes[code]
                    suffixes[code] += 1
                    df.loc[df[column] == code, "status"] = "split"
                    return f"{code}-{suffix}"
                return code

            df[column] = df[column].apply(suffix_code)
            return df

        # Apply the function to your GeoDataFrame
        gdf_segments = add_suffixes(gdf_segments, "code")

        return gdf_segments

    result = split_lines_at_points(hydroobjects, combined_end_points_gdf)
    return result


# Function to check if line endpoints need to be flipped
def check_and_flip(line, start_points, end_points):
    start_point = line.coords[0]
    end_point = line.coords[-1]

    # Check if start and end points are in the correct positions
    if start_point in end_points and end_point in start_points:
        return line, False
    else:
        return LineString(line.coords[::-1]), True


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


def calculate_angle(line, direction):
    if direction == "upstream":
        coords = list(line.coords)
        p1, p2 = coords[0], coords[1]  # First segment
    elif direction == "downstream":
        coords = list(line.coords)
        p1, p2 = coords[-2], coords[-1]  # Last segment
    else:
        raise ValueError("Direction must be 'upstream' or 'downstream'")

    angle = np.arctan2(p2[0] - p1[0], p2[1] - p1[1])  # Angle in radians
    angle_degrees = np.degrees(angle)
    return angle_degrees


def calculate_angle_reverse(line):
    coords = list(line.coords)
    p1, p2 = coords[1], coords[0]

    angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # Angle in radians
    angle_degrees = np.degrees(angle)
    angle_degrees = angle_degrees % 360

    return angle_degrees


def calculate_angle_start(line):
    coords = list(line.coords)
    p1, p2 = coords[0], coords[1]

    angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # Angle in radians
    angle_degrees = np.degrees(angle)
    angle_degrees = angle_degrees % 360

    return angle_degrees


def calculate_angle_end(line):
    coords = list(line.coords)
    p1, p2 = coords[-2], coords[-1]

    angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # Angle in radians
    angle_degrees = np.degrees(angle)
    angle_degrees = angle_degrees % 360

    return angle_degrees


def calculate_angle_difference(angle1, angle2):
    diff = abs(angle1 % 360 - angle2 % 360)
    if diff > 180:
        diff = 360 - diff
    return diff


def find_edge_smallest_angle_difference(reference_angle, angles, edge_codes):
    if reference_angle is None:
        return [None for a in angles], None
    reference_angle = float(reference_angle)
    angle_differences = np.array(
        [calculate_angle_difference(angle, reference_angle) for angle in angles]
    )
    min_index = np.argmin(angle_differences)
    return angle_differences, edge_codes[min_index]


def remove_z_dims(_gdf: gpd.GeoDataFrame):
    _gdf.geometry = [
        (Point(g.coords[0][:2]) if len(g.coords[0]) > 2 else Point(g.coords[0]))
        if isinstance(g, Point)
        else (
            (LineString([c[:2] if len(c) > 2 else c for c in g.coords]))
            if isinstance(g, LineString)
            else Polygon([c[:2] if len(c) > 2 else c for c in g.exterior.coords])
        )
        for g in _gdf.geometry.values
    ]
    return _gdf


def get_endpoints_from_lines(lines: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Extract all unique endpoints of line features from vector data

    Args:
        lines (gpd.GeoDataFrame): GeoDataFrame containing line features

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing all unique endpoints from
        line features
    """
    lines[["startpoint", "endpoint"]] = lines["geometry"].apply(
        lambda x: pd.Series([x.coords[0], x.coords[-1]])
    )
    endpoints = pd.unique(lines[["startpoint", "endpoint"]].values.ravel("K"))
    endpoints = gpd.GeoDataFrame({"coordinates": endpoints})
    endpoints["starting_lines"] = endpoints["coordinates"].apply(
        lambda x: lines["code"][lines["startpoint"] == x].values
    )
    endpoints["ending_lines"] = endpoints["coordinates"].apply(
        lambda x: lines["code"][lines["endpoint"] == x].values
    )
    endpoints["starting_line_count"] = endpoints.apply(
        lambda x: len(list(x["starting_lines"])), axis=1
    )
    endpoints["ending_line_count"] = endpoints.apply(
        lambda x: len(list(x["ending_lines"])), axis=1
    )
    endpoints["connected_line_count"] = endpoints.apply(
        lambda x: x["starting_line_count"] + x["ending_line_count"], axis=1
    )
    endpoints_geometry = endpoints.coordinates.apply(lambda x: Point(x))
    endpoints = endpoints.set_geometry(endpoints_geometry)
    return endpoints


def snap_unconnected_endpoints_to_endpoint_or_line(
    hydroobject, snapping_distance=0.05
):
    hydroobject = hydroobject.explode()

    endpoints = get_endpoints_from_lines(hydroobject)

    endpoints["ID"] = endpoints.index
    endpoints = endpoints[
        (endpoints["starting_line_count"] == 0) | (endpoints["ending_line_count"] == 0)
    ]
    endpoints = endpoints.rename(columns={"coordinates": "geometry"})
    endpoints = gpd.GeoDataFrame(endpoints, geometry="geometry", crs=28992)

    # Create buffer and copy of original_geometry
    original_geometry = endpoints[["ID", "geometry"]].copy()
    endpoints["geometry"] = endpoints.buffer(snapping_distance)

    endpoints = endpoints.to_crs(28992)
    original_geometry = original_geometry.to_crs(28992)
    # Perform a spatial join to find endpoints within the buffer
    joined_points = gpd.sjoin(
        original_geometry, endpoints, how="inner", predicate="intersects"
    )
    # Filter out the endpoints with the same ID
    unconnected_endpoints_points = joined_points[
        joined_points["ID_left"] != joined_points["ID_right"]
    ]

    merged_df = unconnected_endpoints_points.merge(
        unconnected_endpoints_points,
        left_on="ID_left",
        right_on="ID_right",
        suffixes=("_left", "_right"),
    )

    # Explode the DataFrame to handle list elements individually
    exploded_df = merged_df.explode("ending_lines_left")

    # Filtering out rows where 'ending_lines_left' is nan
    exploded_df = exploded_df[exploded_df["ending_lines_left"].notna()]

    point_df = (
        exploded_df[["ending_lines_left", "geometry_left"]]
        .dropna()
        .reset_index(drop=True)
    )
    point_df.columns = ["code", "geometry_left"]

    def replace_last_coordinate(row, point_df):
        if row["code"] in point_df["code"].values:
            new_point = point_df.loc[
                point_df["code"] == row["code"], "geometry_left"
            ].values[0]
            line = row["geometry"]
            new_coords = list(line.coords[:-1]) + [new_point.coords[0]]
            return LineString(new_coords), "snapped"
        return row["geometry"], "unchanged"

    hydroobject[["geometry", "status"]] = hydroobject.apply(
        lambda row: replace_last_coordinate(row, point_df), axis=1, result_type="expand"
    )

    joined_lines = gpd.sjoin(
        endpoints, hydroobject, how="left", predicate="intersects"
    )

    joined_lines = joined_lines[
        joined_lines.apply(
            lambda x: x["code"] not in x["starting_lines"]
            and x["code"] not in x["ending_lines"],
            axis=1,
        )
    ]

    # Merge the DataFrames on the 'ID' column to bring the correct geometry back
    joined_lines = joined_lines.merge(
        original_geometry[["ID", "geometry"]], on="ID", suffixes=("", "_original")
    )

    # Replace the geometry column in joined_lines with the one from original_geometry
    joined_lines["geometry"] = joined_lines["geometry_original"]

    # Drop the now redundant 'geometry_original' column
    joined_lines.drop(columns=["geometry_original"], inplace=True)

    joined_lines["line_key"] = joined_lines.apply(
        lambda row: str(
            sorted(
                row["starting_lines"] if isinstance(row["starting_lines"], list) else []
            )
        )
        + "-"
        + str(
            sorted(row["ending_lines"] if isinstance(row["ending_lines"], list) else [])
        ),
        axis=1,
    )
    unconnected_endpoints_points["point_key"] = unconnected_endpoints_points.apply(
        lambda row: str(
            sorted(
                row["starting_lines"] if isinstance(row["starting_lines"], list) else []
            )
        )
        + "-"
        + str(
            sorted(row["ending_lines"] if isinstance(row["ending_lines"], list) else [])
        ),
        axis=1,
    )
    filtered_joined_lines = joined_lines[
        ~joined_lines["line_key"].isin(unconnected_endpoints_points["point_key"])
    ]

    def project_point_onto_line(point, line):
        if point is None or line is None:
            return None  # Handle None values gracefully
        # Use shapely's nearest_points to find the closest point on the line
        projected_point = nearest_points(point, line)[1]
        return projected_point

    # Create a new column for the projected geometries
    filtered_joined_lines["projected_geometry"] = filtered_joined_lines.apply(
        lambda row: project_point_onto_line(
            row["geometry"]
            if isinstance(row["geometry"], (shapely.geometry.base.BaseGeometry, list))
            else None,
            hydroobject.loc[hydroobject["code"] == row["code"], "geometry"].values[
                0
            ]
            if not pd.isna(row["code"])
            and not hydroobject.loc[
                hydroobject["code"] == row["code"], "geometry"
            ].empty
            else None,
        ),
        axis=1,
    )

    # Explode the lists in starting_lines and ending_lines
    filtered_joined_lines_exploded = filtered_joined_lines.explode(
        "starting_lines"
    ).explode("ending_lines")

    # Create a new DataFrame from filtered_joined_lines_exploded
    point_df_line = filtered_joined_lines_exploded[
        ["starting_lines", "ending_lines", "projected_geometry"]
    ]

    def replace_coordinate_with_projection(row, point_df):
        if (
            row["code"] in point_df["starting_lines"].values
            or row["code"] in point_df["ending_lines"].values
        ):
            new_point = point_df.loc[
                (point_df["starting_lines"] == row["code"])
                | (point_df["ending_lines"] == row["code"]),
                "projected_geometry",
            ].values[0]
            line = row["geometry"]
            if row["code"] in point_df["starting_lines"].values:
                # Replace the first coordinate
                new_coords = [new_point.coords[0]] + list(line.coords[1:])
            elif row["code"] in point_df["ending_lines"].values:
                # Replace the last coordinate
                new_coords = list(line.coords[:-1]) + [new_point.coords[0]]
            return LineString(new_coords), "snapped"
        return row["geometry"], row["status"]

    hydroobject[["geometry", "status"]] = hydroobject.apply(
        lambda row: replace_coordinate_with_projection(row, point_df_line),
        axis=1,
        result_type="expand",
    )

    return hydroobject


def check_duplicate_codes(gdf, column):
    temp_column = gdf[column].copy()
    duplicates = gdf[column][gdf[column].duplicated(keep=False)]

    suffix_dict = {}
    for idx, value in duplicates.items():
        if value not in suffix_dict:
            suffix_dict[value] = 0
        else:
            suffix_dict[value] += 1
        temp_column[idx] = f"{value}-{suffix_dict[value]}"

    gdf[column] = temp_column
    return gdf
