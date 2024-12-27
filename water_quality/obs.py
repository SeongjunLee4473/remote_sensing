import pandas as pd
import geopandas as gpd
import folium
import rasterio
#---------------------------------------------------------------------------------------------------#
# Function to plot points based on ptNo
'''
This function takes a CSV file path and a point number (ptNo) as inputs,
and plots the corresponding points on a folium map.

CSV structure: The CSV file should contain the following columns:
- ptNo: Point number (example: 3008B40)
- ptNm: Point name (example: 대청댐1)
- lon: Longitude (decimal degree, example: 127.123456)
- lat: Latitude (decimal degree, example: 35.123456)

Last modified: 2024-12-25
Created on: 2024-12-25
Created by: Seongjun Lee
'''
def plot_point(csv_file_path, ptNo):
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
plot_point(csv_file_path, point_code)
'''
#---------------------------------------------------------------------------------------------------#
# Function to load and prepare data
"""
This function loads a CSV file containing point coordinates (longitude and latitude),
creates a GeoDataFrame (Point) from these coordinates.

CSV structure: The CSV file should contain the following columns:
- ptNo: Point number (example: 3008B40)
- ptNm: Point name (example: 대청댐1)
- lon: Longitude (decimal degree, example: 127.123456)
- lat: Latitude (decimal degree, example: 35.123456)

Created on: 2024-12-27
Created by: Seongjun Lee
"""
def create_point_gdf(csv_file_path):

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
obs_geometry_gdf = create_point_gdf(csv_file_path)
'''
#---------------------------------------------------------------------------------------------------#
# Load raster data
"""
This function loads raster data from a specified file and returns the raster data
and its affine transform.

Created on: 2024-12-27
Created by: Seongjun Lee
"""
def load_raster(raster_file_path):
    # Load Pekel Data
    dataset = rasterio.open(raster_file_path)

    # Get the affine transform and data from the dataset
    raster_transform = dataset.transform
    raster_data = dataset.read(1)

    return raster_data, raster_transform

'''
# Example usage:
raster_file_path = 'path/to/raster.tif'
raster_data, raster_transform = load_raster(raster_file_path)
'''
#---------------------------------------------------------------------------------------------------#