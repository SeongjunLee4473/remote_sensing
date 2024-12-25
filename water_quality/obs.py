import pandas as pd
import folium
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