#---------------------------------------------------------------------------------------------------#
# Import libraries
import geopandas as gpd
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Read vector data
def read(vector_file_path):
    """
    This function reads a vector file and returns a GeoDataFrame.

    Parameters:
    - file_path: Path to the vector file (e.g., Shapefile)

    Returns:
    - GeoDataFrame containing the vector data.
    
    Raises:
    - ValueError if the file does not exist or is not a valid vector file.
    """
    if not os.path.isfile(vector_file_path):
        raise ValueError(f"Vector file '{vector_file_path}' does not exist.")
    
    if not vector_file_path.endswith('.shp'):
        raise ValueError("The file must be a .shp file.")
    
    return gpd.read_file(vector_file_path)

'''
# Example usage:
vector_file_path = 'path/to/vector.shp'
vector_data = read(vector_file_path)
'''
#---------------------------------------------------------------------------------------------------#
# Check the coordinate reference system (CRS) of a vector file
def check_crs(shp_file):
    """
    This function checks the coordinate reference system (CRS) of a vector file.

    Parameters:
    - shp_file: Path to the input vector file (e.g., Shapefile) or GeoDataFrame

    Returns:
    - The EPSG code of the CRS if defined, otherwise a message indicating that the CRS is not defined.
    """
    # Check if the input is a GeoDataFrame or a file path
    if isinstance(shp_file, gpd.GeoDataFrame):
        vector_gdf = shp_file  # Use the provided GeoDataFrame directly
    elif isinstance(shp_file, str) and shp_file.endswith('.shp'):
        vector_gdf = gpd.read_file(shp_file)  # Load the vector data from the file
    else:
        raise ValueError("Input must be a GeoDataFrame or a valid .shp file path.")

    if vector_gdf.crs is None:
        return "The vector file does not have a defined coordinate reference system (CRS)."
    else:
        return f"EPSG:{vector_gdf.crs.to_epsg()}"
    
'''
# Example usage:
shp_file = 'path/to/vector.shp'
crs_info = check_crs(shp_file)
print(crs_info)
'''
#---------------------------------------------------------------------------------------------------#
def reproject(input_path, output_path, epsg_code):
    """
    This function reprojects a vector file to a specified EPSG coordinate system.

    Parameters:
    - input_path: Path to the input vector file (e.g., Shapefile)
    - output_path: Path to save the reprojected vector file
    - epsg_code: EPSG code of the desired coordinate system (example: 3857 -> WGS 84 Pseudo-Mercator)

    Last modified: None / None
    Created on: 2024-12-27
    Created by: Seongjun Lee
    """
    # Check if the input is a GeoDataFrame or a file path
    if isinstance(input_path, gpd.GeoDataFrame):
        vector_gdf = input_path  # Use the provided GeoDataFrame directly
    elif isinstance(input_path, str) and input_path.endswith('.shp'):
        vector_gdf = gpd.read_file(input_path)  # Load the vector data from the file
    else:
        raise ValueError("Input must be a GeoDataFrame or a valid .shp file path.")

    # Check if the CRS is defined
    if vector_gdf.crs is None:
        raise ValueError("The input vector file does not have a defined coordinate reference system (CRS).")

    # Check if the input CRS is the same as the desired EPSG code
    if vector_gdf.crs.to_epsg() == epsg_code:
        print("The input and output CRS are the same. No reprojection needed.\n"
              f"Current CRS: {vector_gdf.crs.to_string()}")
        return  # Exit the function if no reprojection is needed

    # Ensure the output path ends with .shp
    if not output_path.lower().endswith('.shp'):
        raise ValueError("The output must be a .shp file.")

    # Reproject the GeoDataFrame to the target EPSG code
    vector_gdf = vector_gdf.to_crs(epsg=epsg_code)

    # Save the reprojected vector data to the output path
    vector_gdf.to_file(output_path)

    return vector_gdf

'''
# Example usage:
input_path = 'path/to/input.shp'
output_path = 'path/to/output.shp'
epsg_code = 4326
reprojected_gdf = reproject(input_path, output_path, epsg_code)
'''
#---------------------------------------------------------------------------------------------------#
def plot(vector_file,
                bounds=None, 
                title=None, 
                label=None, 
                cmap=None,
                tick_interval=None,
                grid=True):
    '''
    This function plots vector data.

    Last modified: None / None
    Created on: 2024-12-27
    Created by: Seongjun Lee
    '''
    import geopandas as gpd

    # Check if vector_file is a valid path or a GeoDataFrame
    if isinstance(vector_file, str):
        if vector_file.endswith('.shp'):
            if not os.path.isfile(vector_file):
                raise ValueError(f"Vector file '{vector_file}' does not exist.")
            vector_data = gpd.read_file(vector_file)  # Read if it's a .shp file
        else:
            raise ValueError("vector_file must be a .shp file path.")
    elif isinstance(vector_file, gpd.GeoDataFrame):
        vector_data = vector_file  # If it's already a GeoDataFrame
    else:
        raise ValueError("vector_file must be a valid file path or a GeoDataFrame.")

    # Check if the label exists in the DataFrame
    if label is not None and label not in vector_data.columns:
        raise ValueError(f"Label '{label}' not found in the vector data columns: {vector_data.columns.tolist()}")

    # Check CRS
    if vector_data.crs is None:
        raise ValueError("The CRS of the vector data is not defined.")
    elif vector_data.crs.to_string() != 'EPSG:4326':
        print(f"Warning: This code is optimized for the WGS 84 coordinate system (EPSG:4326).\n"
              f"Please reproject to EPSG:4326 for accurate plotting.\n"
              f"Current CRS: {vector_data.crs.to_string()}.")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    vector_data.plot(ax=ax, cmap=cmap if cmap is not None else 'viridis')

    # Set the limits based on the provided bounds
    if bounds is not None:
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(bounds[2], bounds[3]) 

    # X/Y Ticks
    if tick_interval is not None:
        x_ticks = np.arange(bounds[0], bounds[1] + tick_interval, tick_interval)
        y_ticks = np.arange(bounds[2], bounds[3] + tick_interval, tick_interval)
        plt.xticks(x_ticks, fontsize=13, rotation=45)
        plt.yticks(y_ticks, fontsize=13)
    else:
        plt.xticks(fontsize=13, rotation=45)
        plt.yticks(fontsize=13)

    # Title
    plt.title(title if title is not None else 'Vector Data Visualization', fontsize=20, fontweight='bold')

    # Colorbar (if applicable)
    if label is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap if cmap is not None else 'viridis', 
                                   norm=plt.Normalize(vmin=vector_data[label].min(), 
                                                      vmax=vector_data[label].max()))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(label, fontsize=15)

    # Grid
    if grid:
        ax.grid(True, linestyle='--', color='gray', alpha=0.5)

    # Show the plot
    plt.show()

'''
# Example usage:
1. If you want to directly load vector file
vector_file_path = 'path/to/vector.shp'
plot(vector_file=vector_file_path)

2. If you already read vector data as GeoDataFrame
plot(vector_file=vector_data)
'''
#---------------------------------------------------------------------------------------------------#