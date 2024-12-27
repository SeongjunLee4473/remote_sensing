#---------------------------------------------------------------------------------------------------#
# Import libraries
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import folium
#---------------------------------------------------------------------------------------------------#
# Function to plot points based on ptNo
def plot_points(csv_file_path, ptNo):
    '''
    This function takes a CSV file path and a point number (ptNo) as inputs,
    and plots the corresponding points on a folium map.

    CSV structure: The CSV file should contain the following columns:
    - ptNo: Point number (example: 3008B40)
    - ptNm: Point name (example: 대청댐1)
    - lon: Longitude (decimal degree, example: 127.123456)
    - lat: Latitude (decimal degree, example: 35.123456)

    Last modified: 2024-12-27 / Rename from plot_point to plot_points
    Created on: 2024-12-25
    Created by: Seongjun Lee
    '''
    # Load CSV data
    df = pd.read_csv(csv_file_path)

    # Check if ptNo is four digits
    if len(str(ptNo)) == 4:
        point = df[df['ptNo'].astype(str).str[:4] == str(ptNo)]
    else:
        point = df[df['ptNo'] == ptNo]
    
    if point.empty:
        print(f"ptNo {ptNo} not found.")
        return

    # Create a folium map centered at the average of points
    center_lat = point['lat'].mean()
    center_lon = point['lon'].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=12, tiles='OpenStreetMap',
                    attr='Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors')
    
    for _, row in point.iterrows():
        folium.Marker(
            location=[row['lat'], row['lon']],
            popup=f"{row['ptNm']} (ptNo: {row['ptNo']})",
            icon=folium.Icon(color='blue')
        ).add_to(m)
    
    # Display the map directly
    return m

'''
# Example usage
csv_file_path = '/home/seongjun/water_quality/obs_data/obs_coordinate_2024.csv'
point_code = '3008B40'  # Change to the ptNo you want to plot (If you type 4 digits, all points with the same 4 digits will be plotted)
plot_points(csv_file_path, point_code)
'''
#---------------------------------------------------------------------------------------------------#
# Function to load and prepare data
def points_to_gdf(csv_file_path):
    """
    This function loads a CSV file containing point coordinates (longitude and latitude),
    creates a GeoDataFrame (Point) from these coordinates.

    CSV structure: The CSV file should contain the following columns:
    - ptNo: Point number (example: 3008B40)
    - ptNm: Point name (example: 대청댐1)
    - lon: Longitude (decimal degree, example: 127.123456)
    - lat: Latitude (decimal degree, example: 35.123456)

    Last modified: 2024-12-27 / Rename from create_point_gdf to points_to_gdf
    Created on: 2024-12-27
    Created by: Seongjun Lee
    """
    # Load the filtered GeoDataFrame from the CSV file
    df = pd.read_csv(csv_file_path)

    # Create a GeoDataFrame from the DataFrame using lon and lat
    point_gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df['lon'], df['lat']),
        crs="EPSG:4326"
    )

    # Drop the original lon and lat columns
    point_gdf.drop(columns=['lon', 'lat'], inplace=True)

    return point_gdf

'''
# Example usage:
csv_file_path = 'path/to/obs_data.csv'
obs_geometry_gdf = points_to_gdf(csv_file_path)
'''
#---------------------------------------------------------------------------------------------------#
# Load raster data
def read_raster(raster_file_path):
    '''
    This function loads raster data from a specified file and returns the raster data
    and its affine transform.

    Last modified: 2024-12-27 / Rename from load_raster to read_raster
    Created on: 2024-12-27
    Created by: Seongjun Lee
    '''
    # Load Pekel Data
    dataset = rasterio.open(raster_file_path)

    # Get the affine transform and data from the dataset
    raster_transform = dataset.transform
    raster_data = dataset.read(1)

    return raster_data, raster_transform

'''
# Example usage:
raster_file_path = 'path/to/raster.tif'
raster_data, raster_transform = read_raster(raster_file_path)
'''
#---------------------------------------------------------------------------------------------------#
# Function to extract pixel values around points
def extract_pixels_near_points(gdf, raster_data, transform, radius_deg, min_pixel_value=None, max_pixel_value=None):
    """
    This function extracts pixel values near points within a specified radius and pixel value range.

    Parameters:
    - gdf: GeoDataFrame (must have the following columns: ptNo, geometry)
    - raster_data: Raster data
    - transform: Raster transform
    - radius_deg: Radius in degrees
    - min_pixel_value: Minimum pixel value
    - max_pixel_value: Maximum pixel value

    Returns:
    Final output is a DataFrame with the following columns:
    - ptNo: Point number
    - point: Point coordinates
    - pixel_value: Pixel value
    - pixel_lon: Pixel longitude
    - pixel_lat: Pixel latitude
    """
    # Initialize an empty list to store results
    results = []

    # Calculate the pixel radius from the degree radius
    pixel_radius_x = int(np.ceil(radius_deg / abs(transform[0])))  # Width in pixels
    pixel_radius_y = int(np.ceil(radius_deg / abs(transform[4])))  # Height in pixels

    # Iterate over each point in the GeoDataFrame
    for idx, row in gdf.iterrows():
        point = row['geometry']
        pt_no = row['ptNo']
        coords = (point.x, point.y)

        # Convert geographic coordinates to pixel coordinates
        pixel_x, pixel_y = ~transform * coords
        pixel_x = int(round(pixel_x))
        pixel_y = int(round(pixel_y))

        # Loop over the square of pixels around the center point
        for dx in range(-pixel_radius_x, pixel_radius_x + 1):
            for dy in range(-pixel_radius_y, pixel_radius_y + 1):
                # Convert pixel indices back to geographic coordinates
                x = pixel_x + dx
                y = pixel_y + dy

                # Check bounds to avoid indexing errors
                if 0 <= x < raster_data.shape[1] and 0 <= y < raster_data.shape[0]:
                    pixel_value = raster_data[y, x]  # Note: raster data is row-major (y, x)

                    # Extract pixel values within specified range if min/max values are provided
                    if (min_pixel_value is None or pixel_value >= min_pixel_value) and \
                       (max_pixel_value is None or pixel_value <= max_pixel_value):
                        lon, lat = transform * (x, y)
                        distance_deg = np.sqrt((lon - coords[0])**2 + (lat - coords[1])**2)

                        # Check if the geographic distance is within the radius
                        if distance_deg <= radius_deg:
                            results.append({
                                'ptNo': pt_no,
                                'point': coords,
                                'pixel_value': pixel_value,
                                'pixel_lon': lon,
                                'pixel_lat': lat
                            })

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)

    return results_df

'''
# Example usage
radius_deg = 0.001
min_pixel_value = 30
max_pixel_value = 100
pekel_results = extract_pixels_near_points(obs_geometry_gdf, pekel_data, pekel_transform, radius_deg, min_pixel_value, max_pixel_value)
'''
#---------------------------------------------------------------------------------------------------#