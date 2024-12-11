import geopandas as gpd
import pandas as pd
from shapely import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import split, linemerge
import numpy as np
from typing import TypeVar
import time
import logging
import warnings


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

    hydroobjects = gpd.GeoDataFrame(hydroobjects, geometry="geometry", crs=28992)

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
            gdf_lines[["CODE", "segments"]]
            .explode("segments")
            .rename(columns={"segments": "geometry"}),
            geometry="geometry",
            crs=gdf_lines.crs,
        ).reset_index(drop=True)

        # Determine which segment point of points_gdf are located on
        gdf_segments = gdf_segments.merge(
            gdf_points[["dangling_id", "geometry"]]
            .sjoin_nearest(gdf_segments.drop(columns="CODE"))
            .groupby("index_right")
            .agg(
                {
                    "dangling_id": list,
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
            gdf_segments.groupby(["index_right_cumsum", "CODE"])
            .agg(
                {
                    "CODE": "first",
                    "geometry": list,
                    "dangling_id": "first",
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
        gdf_segments = add_suffixes(gdf_segments, "CODE")

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
            lambda x: ",".join(
                list(edges.loc[edges[node] == x["nodeID"], "code"].values) if len(edges.loc[edges[node] == x["nodeID"]])>0 else []
            ),
            axis=1,
        )
        nodes_sel[f"no_{direction}_edges"] = nodes_sel.apply(
            lambda x: len(x[f"{direction}_edges"].split(",")), axis=1
        )
    nodes_sel = nodes_sel.reset_index(drop=True)
    return nodes_sel


def calculate_angle(line, direction):
    if direction == "downstream":
        # Get the first segment for downstream
        coords = list(line.coords)
        p1, p2 = coords[0], coords[1]  # First segment
    elif direction == "upstream":
        # Get the last segment for upstream
        coords = list(line.coords)
        p1, p2 = coords[-2], coords[-1]  # Last segment
    else:
        raise ValueError("Direction must be 'upstream' or 'downstream'")

    # Calculate the angle relative to the north (0 degrees)
    angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # Angle in radians
    angle_degrees = np.degrees(angle)
    angle_degrees = angle_degrees % 360  # Normalize to 0-360 degrees

    return angle_degrees


def find_closest_edge(reference_angle, angles, edge_codes):
    angles = np.array(
        angles, dtype=float
    )  # Convert to numpy array for easier calculations
    edge_codes = np.array(edge_codes)
    reference_angle = float(reference_angle)

    # Calculate the angle differences
    angle_differences = np.abs(angles - reference_angle)

    # Find the index of the minimum angle difference
    min_index = np.argmin(angle_differences)

    return edge_codes[min_index]


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

def connect_lines_by_endpoints(
    split_endpoints: gpd.GeoDataFrame, lines: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Connects boundary lines to other lines, based on instructions for each line endpoint. Connections are
    created by inserting the endpoints into their target lines. The target line features are
    split afterwards in order to create nodes at the intersection of the connected linestrings.

    Args:
        split_endpoints (gpd.GeoDataFrame): Dataframe containing line endpoints and instructions
        lines (list): Vector Dataset containing line features

    Returns:
        gpd.GeoDataFrame: line feature dataframe
    """
    split_endpoints["lines_to_add_point"] = split_endpoints.apply(
        lambda x: np.append(x["starting_lines"], x["ending_lines"]), axis=1
    )
    split_endpoints["points_to_add"] = split_endpoints.apply(
        lambda x: np.array([Point(p) for p in x["points_to_add"]]), axis=1
    )
    split_endpoints["lines_to_add_point"] = split_endpoints.apply(
        lambda x: x["unconnected_lines"]
        if x["points_to_add"].size == 0
        else x["lines_to_add_point"],
        axis=1,
    )
    split_endpoints["points_to_add"] = split_endpoints.apply(
        lambda x: [Point(p) for p in set(x["points_to_add"])]
        if list(x["lines_to_add_point"]) != list(x["unconnected_lines"])
        else [x["coordinates"]],
        axis=1,
    )
    connections_to_create = split_endpoints[["lines_to_add_point", "points_to_add"]]
    connections_to_create["inserted"] = False

    # A dataframe is created to store the resulting linestrings from the splits
    split_lines = gpd.GeoDataFrame(columns=lines.columns)
    split_lines["preprocessing_split"] = None
    split_lines["new_code"] = split_lines["code"]

    for ind, split_action in split_endpoints.iterrows():
        logging.info(split_action["lines_to_add_point"])
        line = lines[lines["code"].isin(split_action["lines_to_add_point"])].iloc[0]
        linestring = line["geometry"]
        nodes_to_add = []
        for node in split_action["points_to_add"]:
            new_linestring, index_new_point = add_point_to_linestring(node, linestring)

            if index_new_point == 0 and linestring.coords[0] in list(
                connections_to_create.loc[
                    connections_to_create["inserted"] == True, "points_to_add"
                ].values
            ):
                continue

            elif index_new_point == len(linestring.coords) - 1 and linestring.coords[
                -1
            ] in list(
                connections_to_create.loc[
                    connections_to_create["inserted"] == True, "points_to_add"
                ].values
            ):
                continue

            new_linestring = new_linestring.simplify(tolerance=0.0)
            if linestring == new_linestring:
                connections_to_create = connections_to_create[
                    [
                        False if split_action["lines_to_add_point"][0] in c else True
                        for c in connections_to_create["lines_to_add_point"].values
                    ]
                ]
                continue

            linestring = new_linestring

            connections_to_create["inserted"][
                connections_to_create["points_to_add"] == node
            ] = True

            nodes_to_add += [node]

        # The modified line will be divided in seperate linestrings
        split_indices = [
            list(linestring.coords).index(node.coords[0]) for node in nodes_to_add
        ]
        split_linestrings = split_linestring_by_indices(linestring, split_indices)

        # Append split linestrings to collection of new lines
        for k, split_linestring in enumerate(split_linestrings):
            snip_line = line.copy()
            snip_line["geometry"] = split_linestring
            snip_line["preprocessing_split"] = "split"
            snip_line["new_code"] = f'{snip_line["code"]}-{k}'
            split_lines = pd.concat(
                [split_lines, snip_line.to_frame().T], axis=0, join="inner"
            )

    # Remove lines that have been divided from original geodataframe, and append resulting lines
    unedited_lines = lines[
        ~lines["code"].isin(
            [x for c in connections_to_create["lines_to_add_point"].values for x in c]
        )
    ]
    lines = pd.concat([unedited_lines, split_lines], axis=0, join="outer").reset_index(
        drop=True
    )

    # Remove excessive split lines
    lines.geometry = lines.geometry.simplify(tolerance=0.0)
    lines = lines.drop_duplicates(subset="geometry", keep="first")
    lines = lines[lines.geometry.length > 0]
    return lines


def connect_endpoints_by_buffer(
    lines: gpd.GeoDataFrame, buffer_distance: float = 0.5
) -> gpd.GeoDataFrame:
    """
    Connects boundary line endpoints within a vector dataset to neighbouring lines that pass
    within a specified buffer distance with respect to the the boundary endpoints. The boundary
    endpoints are inserted into the passing linestrings

    Args:
        lines (gpd.GeoDataFrame): Line vector data
        buffer distance (float): Buffer distance for connecting line boundary endpoints, expressed
        in the distance unit of vector line dataset

    Returns:
        gpd.Geo: list of resulting linestrings
    """
    warnings.filterwarnings("ignore")
    start_time = time.time()
    iterations = 0
    unconnected_endpoints_count = 0
    finished = False
    lines["preprocessing_split"] = None

    logging.info(
        f"Detect unconnected endpoints nearby linestrings, buffer distance: {buffer_distance}m"
    )

    while not finished:
        endpoints = get_endpoints_from_lines(lines)

        boundary_endpoints = gpd.GeoDataFrame(
            endpoints[
                (endpoints["starting_line_count"] == 0)
                | (endpoints["ending_line_count"] == 0)
            ]
        )
        lines["buffer_geometry"] = lines.geometry.buffer(
            buffer_distance, join_style="round"
        )

        boundary_endpoints["overlaying_line_buffers"] = list(
            map(
                lambda x: lines[lines.buffer_geometry.contains(x)].code.tolist(),
                boundary_endpoints.geometry,
            )
        )

        # check on unconnected endpoints
        boundary_endpoints["unconnected_lines"] = boundary_endpoints.apply(
            lambda x: np.array(
                [
                    b
                    for b in x["overlaying_line_buffers"]
                    if (b not in x["starting_lines"] and b not in x["ending_lines"])
                ]
            ),
            axis=1,
        )
        boundary_endpoints["startpoint_unconnected_lines"] = boundary_endpoints.apply(
            lambda x: lines[lines.code.isin(x["unconnected_lines"])].endpoint.values,
            axis=1,
        )
        boundary_endpoints["endpoint_unconnected_lines"] = boundary_endpoints.apply(
            lambda x: lines[lines.code.isin(x["unconnected_lines"])].startpoint.values,
            axis=1,
        )
        boundary_endpoints["unconnected_lines_present"] = boundary_endpoints.apply(
            lambda x: len(x["unconnected_lines"]) > 0, axis=1
        )

        unconnected_endpoints = boundary_endpoints[
            boundary_endpoints["unconnected_lines_present"] == True
        ].reset_index(drop=True)

        if unconnected_endpoints.empty:
            unconnected_endpoints["unconnected_lines_present"] = None
            unconnected_endpoints["points_to_add"] = None
        else:
            unconnected_endpoints["points_to_add"] = unconnected_endpoints.apply(
                lambda x: np.array(
                    [
                        p
                        for p in np.append(
                            x["startpoint_unconnected_lines"],
                            x["endpoint_unconnected_lines"],
                        )
                        if x["coordinates"].buffer(buffer_distance).covers(Point(p))
                    ]
                ),
                axis=1,
            )

        previous_unconnected_endpoints_count = unconnected_endpoints_count
        unconnected_endpoints_count = len(unconnected_endpoints)
        if iterations == 0:
            unconnected_endpoints_count_total = unconnected_endpoints_count
        logging.info(
            f"{unconnected_endpoints_count} unconnected endpoints nearby intersecting lines (iteration {iterations})"
        )
        if (
            unconnected_endpoints_count != 0
            and unconnected_endpoints_count != previous_unconnected_endpoints_count
        ):
            lines = connect_lines_by_endpoints(unconnected_endpoints, lines)
            iterations += 1
        else:
            lines = lines.drop(["startpoint", "endpoint", "buffer_geometry"], axis=1)
            finished = True

    end_time = time.time()
    passed_time = end_time - start_time
    logging.info(
        f"Summary:\n\n\
          Detected unconnected endpoints nearby intersecting lines: {unconnected_endpoints_count_total} \n\
          Connected endpoints: {unconnected_endpoints_count_total-unconnected_endpoints_count} \n\
          Remaining unconnected endpoints: {unconnected_endpoints_count}\n\
          Iterations: {iterations}"
    )
    logging.info(f"Finished within {passed_time}")
    return lines

def split_linestring_by_indices(linestring: LineString, split_indices: list) -> list:
    """
    Divides a linestring into multiple linestrings based on a list that contains
    the indices of the split points within the linestring

    Args:
        linestring (LineString): Linestring
        split_indices (list): List of indices of split nodes within the linestring

    Returns:
        list: list of resulting linestrings
    """
    split_linestrings = []
    split_indices = sorted(
        list(set([0] + split_indices + [len(linestring.coords) - 1]))
    )
    for i in range(len(split_indices) - 1):
        split_linestrings.append(
            LineString(linestring.coords[split_indices[i] : split_indices[i + 1] + 1])
        )

    return split_linestrings

def add_point_to_linestring(point: Point, linestring: LineString) -> LineString:
    """
    Inserts point into a linestring, placing the point next to its
    nearest neighboring point in a way that minimizes the total length
    of the linestring.

    Args:
        point (Point): point
        linestring (LineString): linestring

    Returns:
        LineString: resulting linestring
    """
    distances = [point.distance(Point(line_point)) for line_point in linestring.coords]
    index_nearest_neighbour = distances.index(min(distances))
    modified_linestring1 = LineString(
        list(linestring.coords)[: index_nearest_neighbour + 1]
        + [point.coords[0]]
        + list(linestring.coords)[index_nearest_neighbour + 1 :]
    )
    modified_linestring2 = LineString(
        list(linestring.coords)[:index_nearest_neighbour]
        + [point.coords[0]]
        + list(linestring.coords)[index_nearest_neighbour:]
    )
    modified_linestring = (
        modified_linestring1
        if modified_linestring1.length < modified_linestring2.length
        else modified_linestring2
    )
    return modified_linestring, index_nearest_neighbour

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
