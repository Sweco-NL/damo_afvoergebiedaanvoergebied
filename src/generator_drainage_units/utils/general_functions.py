import geopandas as gpd
from shapely import Point, LineString, MultiLineString


def line_to_vertices(gdf, distance):
    all_points = []

    for _, row in gdf.iterrows():
        geometry = row['geometry']
        
        # Check if the geometry is a LineString or MultiLineString
        if isinstance(geometry, LineString):
            lines = [geometry]
        elif isinstance(geometry, MultiLineString):
            lines = list(geometry.geoms)
        else:
            continue  # Skip if not a recognized geometry typ
        
        for line in lines:
        # Generate start and end points
            start_point = Point(line.coords[0][0], line.coords[0][1])
            end_point = Point(line.coords[-1][0], line.coords[-1][1])
            distance_from_start_end = line.length  
            distance_from_start_start = 0      
            # Add start and end points to the list with 'dangling' label
            all_points.append({'geometry': start_point, 'line_type': 'dangling','distance_from_start': distance_from_start_start, 'CODE':row['CODE']})
            all_points.append({'geometry': end_point, 'line_type': 'dangling','distance_from_start': distance_from_start_end, 'CODE':row['CODE']})

            # Generate intermediate points at specified intervals
            num_points = int(line.length / distance)
            for i in range(1, num_points):
                new_distance = line.length / num_points
                # Calculate the point at a distance along the line
                point = line.interpolate(i * new_distance)
                point_2d = Point(point.x, point.y)
                # Calculate the cumulative distance from the starting point
                cumulative_distance = i * new_distance
                all_points.append({'geometry': point_2d, 'line_type': 'other', 'distance_from_start': cumulative_distance, 'CODE':row['CODE']})
    
    # Create a new GeoDataFrame for the vertices
    points_gdf = gpd.GeoDataFrame(all_points, crs=gdf.crs)

    return points_gdf