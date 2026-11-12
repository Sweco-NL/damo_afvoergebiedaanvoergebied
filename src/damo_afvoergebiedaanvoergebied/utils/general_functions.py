import geopandas as gpd
import pandas as pd
from shapely.geometry import (
    Point,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
    base,
    MultiPoint,
)
from shapely.ops import split, linemerge, nearest_points, substring
import numpy as np
from typing import TypeVar
import time
import logging
import warnings
import shapely
from collections import defaultdict


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


def line_to_vertices(gdf, distance=10.0, keep_columns=["code", "waterline_type"]):
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
    """
    Split hydroobjects at their own endpoints and at the final endpoints
    of relevant culverts.

    The original point coordinates are retained as the canonical split
    coordinates. This is important because later network construction
    functions require neighbouring geometries to have exactly identical
    endpoint coordinates.
    """

    # ------------------------------------------------------------------
    # Prepare input data
    # ------------------------------------------------------------------

    lines = hydroobjects.copy()
    culverts = relevant_culverts.copy()

    # Explode multipart geometries before calling remove_z_dims().
    # This is necessary because the current remove_z_dims() function does
    # not support MultiLineString geometries.
    lines = (
        lines
        .explode(index_parts=False)
        .reset_index(drop=True)
    )

    culverts = (
        culverts
        .explode(index_parts=False)
        .reset_index(drop=True)
    )

    # Now all MultiLineStrings have become LineStrings.
    lines = remove_z_dims(lines)
    culverts = remove_z_dims(culverts)

    if lines.crs != culverts.crs:
        culverts = culverts.to_crs(lines.crs)

    # Keep only valid LineStrings.
    lines = lines.loc[
        lines.geometry.notna()
        & lines.geometry.geom_type.eq("LineString")
    ].copy()

    culverts = culverts.loc[
        culverts.geometry.notna()
        & culverts.geometry.geom_type.eq("LineString")
    ].copy()

    if lines.empty:
        return lines

    # Keep the columns expected by later functions.
    # if "status" not in lines.columns:
    #     lines["status"] = "unchanged"

    # if "NAAM" not in lines.columns:
    #     lines["NAAM"] = "unknown"

    lines = lines.reset_index(drop=True)
    lines.index = pd.RangeIndex(len(lines))
    lines.index.name = None

    # ------------------------------------------------------------------
    # Create split points
    #
    # The actual coordinates of these points are retained. Do not replace
    # them with the result of line.interpolate(), because the old function
    # used the original point coordinates as network node coordinates.
    # ------------------------------------------------------------------

    points = []

    # Hydroobject start and end points
    for line in lines.geometry:
        points.append(Point(line.coords[0]))
        points.append(Point(line.coords[-1]))

    # Final endpoint of each relevant culvert
    for line in culverts.geometry:
        points.append(Point(line.coords[-1]))

    split_points = gpd.GeoDataFrame(
        {"geometry": points},
        crs=lines.crs,
    )

    split_points = (
        split_points
        .drop_duplicates(subset="geometry")
        .reset_index(drop=True)
    )

    if split_points.empty:
        return lines

    # ------------------------------------------------------------------
    # Match points to complete hydroobjects
    # ------------------------------------------------------------------

    line_geometries = lines[["geometry"]].copy()
    line_geometries.index = pd.RangeIndex(len(line_geometries))
    line_geometries.index.name = None

    matches = gpd.sjoin_nearest(
        split_points,
        line_geometries,
        how="inner",
        distance_col="_distance_to_line",
    )

    if matches.empty:
        return lines

    # ------------------------------------------------------------------
    # Determine split distances
    #
    # The distance is calculated on the hydroobject, but the point stored
    # with it remains the original point from split_points.
    # ------------------------------------------------------------------

    split_locations = {}

    for source_index, point in zip(
        matches["index_right"],
        matches.geometry,
    ):
        line = lines.loc[source_index, "geometry"]

        distance = line.project(point)
        line_length = line.length

        # Snap numerical noise to the actual line endpoints.
        tolerance = max(line_length * 1e-12, 1e-9)

        if distance <= tolerance:
            distance = 0.0

        elif line_length - distance <= tolerance:
            distance = line_length

        split_locations.setdefault(source_index, []).append(
            (distance, point)
        )

    # ------------------------------------------------------------------
    # Split hydroobjects
    # ------------------------------------------------------------------

    output_rows = []

    for source_index, row in lines.iterrows():
        line = row.geometry
        line_length = line.length

        locations = split_locations.get(source_index, [])

        # Sort by position along the hydroobject.
        locations = sorted(
            locations,
            key=lambda item: item[0],
        )

        # Remove duplicate or numerically identical positions.
        unique_locations = []

        tolerance = max(line_length * 1e-12, 1e-9)

        for distance, point in locations:
            if not unique_locations:
                unique_locations.append((distance, point))
                continue

            previous_distance = unique_locations[-1][0]

            if abs(distance - previous_distance) <= tolerance:
                # Keep the first point as the canonical coordinate.
                continue

            unique_locations.append((distance, point))

        # Only interior points cause a split.
        internal_locations = [
            (distance, point)
            for distance, point in unique_locations
            if 0.0 < distance < line_length
        ]

        if not internal_locations:
            output_rows.append(row.to_dict())
            continue

        # Distances used by substring().
        boundaries = [
            0.0,
            *[
                distance
                for distance, point in internal_locations
            ],
            line_length,
        ]

        # Canonical coordinates for every boundary.
        #
        # Interior boundaries use the original split point coordinates.
        # The first and last boundary use the original hydroobject
        # coordinates.
        boundary_points = {
            0.0: Point(line.coords[0]),
            line_length: Point(line.coords[-1]),
        }

        for distance, point in internal_locations:
            boundary_points[distance] = point

        # Create segments.
        for start_distance, end_distance in zip(
            boundaries[:-1],
            boundaries[1:],
        ):
            if end_distance <= start_distance:
                continue

            segment = substring(
                line,
                start_distance,
                end_distance,
            )

            if segment.is_empty or segment.length == 0:
                continue

            start_point = boundary_points[start_distance]
            end_point = boundary_points[end_distance]

            # Explicitly overwrite the endpoints generated by substring().
            #
            # This is the crucial step: the end of one segment and the
            # start of the next segment receive exactly the same coordinate.
            coordinates = list(segment.coords)

            if len(coordinates) < 2:
                continue

            coordinates[0] = start_point.coords[0]
            coordinates[-1] = end_point.coords[0]

            segment = LineString(coordinates)

            values = row.to_dict()
            values["geometry"] = segment
            output_rows.append(values)

    # ------------------------------------------------------------------
    # Create result
    # ------------------------------------------------------------------

    result = gpd.GeoDataFrame(
        output_rows,
        geometry="geometry",
        crs=lines.crs,
    ).reset_index(drop=True)

    result["geometry_len"] = result.geometry.length

    # ------------------------------------------------------------------
    # Add suffixes to duplicate codes
    # ------------------------------------------------------------------

    if "code" in result.columns:
        result["code"] = result["code"].astype(str)

        duplicate_codes = result["code"].duplicated(
            keep=False
        )

        result.loc[duplicate_codes, "status"] = "split"

        occurrence = result.groupby(
            "code",
            sort=False,
        ).cumcount()

        result.loc[duplicate_codes, "code"] = (
            result.loc[duplicate_codes, "code"]
            + "-"
            + occurrence.loc[duplicate_codes].astype(str)
        )

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
    ].copy()
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
