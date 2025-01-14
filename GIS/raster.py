#---------------------------------------------------------------------------------------------------#
# Import libraries
import os
import stat
from tqdm.auto import tqdm
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from osgeo import gdal, osr
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Check the validity of TIFF files
def check_tif_files(input_dir):
    '''
    This function checks the validity of TIFF files in the specified directory.

    It collects all .tif files from the input directory and its subdirectories, 
    and verifies each file to determine if it is valid or invalid. A file is 
    considered valid if it exists, is not empty, and can be opened by GDAL 
    without errors. 
    
    This function explores subdirectories of the input directory.

    The function returns two lists: one containing valid files 
    and another containing invalid files.

    Parameters:
    input_dir (str): The directory path where the TIFF files are located.

    Returns:
    tuple: A tuple containing two lists:
        - valid_files (list): A list of valid TIFF file paths.
        - invalid_files (list): A list of invalid TIFF file paths.
    
    Last modified: None / None
    Created on: 2025-01-08
    Created by: Seongjun Lee
    '''
    valid_files = []
    invalid_files = []
    raster_files = []

    # Enable GDAL exceptions
    gdal.UseExceptions()
    
    # Collect all .tif files from the input directory and its subdirectories
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                raster_files.append(os.path.join(root, file))
    
    for file in raster_files:
        # Check if the file exists and is not empty
        if not os.path.exists(file) or os.path.getsize(file) == 0:
            invalid_files.append(file)
            continue

        try:
            ds = gdal.Open(file)
            if ds is None or ds.RasterCount == 0:
                invalid_files.append(file)
            else:
                try:
                    ds.ReadAsArray()
                except Exception as e:
                    invalid_files.append(file)
                else:
                    valid_files.append(file)
                finally:
                    ds = None
        except RuntimeError as e:  # Handle GDAL errors as RuntimeError
            invalid_files.append(file)
    
    return valid_files, invalid_files

'''
# Example usage
raster_dir = '/data/HLS/Korea/L30/new'
valid_tif_list, invalid_tif_list = check_tif_files(raster_dir)
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Remove invalid TIFFfiles
def remove_invalid_tif(invalid_tif_list):
    '''
    This function removes invalid files from the filesystem.

    Parameters:
    invalid_files (list): A list of file paths to be removed.

    Returns:
    None

    Last modified: None / None
    Created on: 2025-01-08
    Created by: Seongjun Lee
    '''
    for file in tqdm(invalid_tif_list, desc="Removing files"):
        print(f"Attempting to remove: {file}")  # Debugging output
        
        # Change file permissions to allow deletion
        try:
            os.chmod(file, stat.S_IWUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)  # Read and write permissions
        except Exception as e:
            print(f"Error changing permissions for {file}: {e}")
            continue  # Skip to the next file if permission change fails

        try:
            if os.path.exists(file):
                os.remove(file)
                print(f"Removed: {file}")
            else:
                print(f"File not found: {file}")
        except Exception as e:
            print(f"Error removing {file}: {e}")

'''
# Example usage
_, invalid_tif_list = check_tif_files('Your_Path')
remove_invalid_tif(invalid_tif_list)
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Read raster data
def read(raster_file_path):
    '''
    This function loads raster data from a specified file and returns the raster data
    and its affine transform.

    Last modified: 2024-12-27 / Add raster_crs (To check the coordinate reference system)
    Created on: 2024-12-27
    Created by: Seongjun Lee
    '''
    # Load Pekel Data
    dataset = rasterio.open(raster_file_path)

    # Get the affine transform and data from the dataset
    raster_transform = dataset.transform
    raster_data = dataset.read(1)
    
    # Get the CRS from the dataset
    raster_crs = dataset.crs

    return raster_data, raster_transform, raster_crs

'''
# Example usage:
raster_file_path = 'path/to/raster.tif'
raster_data, raster_transform = read_raster(raster_file_path)
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Get lon_1d and lat_1d from the given TIF file path or data and transform.
def get_lon_lat(input_data, transform=None, dem=1):
    """
    Calculate lon_1d and lat_1d from the given TIF file path or data and transform.

    Parameters:
    input_data (str or np.ndarray): The file path to the TIF file or the water percentage data.
    transform (Affine, optional): The affine transformation if input_data is not a file path.
    dem (int, optional): The dimension of the output. 1 for 1D (lon_1d, lat_1d), 2 for 2D (lon_2d, lat_2d).

    Returns:
    tuple: lon_1d and lat_1d arrays or lon_2d and lat_2d arrays based on dem.
    """
    if isinstance(input_data, str) and input_data.endswith('.tif'):
        # Load Pekel data from TIF file
        water_percentage, water_transform, _ = raster.read(input_data)
    elif transform is not None:
        # Use provided data and transform
        water_percentage = input_data
        water_transform = transform
    else:
        raise ValueError("Invalid input: Provide a .tif file path or data with transform.")

    # Extract parameters from the Affine transformation
    pixel_width = water_transform.a
    top_left_x = water_transform.c
    pixel_height = water_transform.e  # This is negative
    top_left_y = water_transform.f

    ncols = water_percentage.shape[1]  # Number of columns
    nrows = water_percentage.shape[0]   # Number of rows

    if dem == 1:
        lon_1d = top_left_x + np.arange(ncols) * pixel_width
        lat_1d = top_left_y + np.arange(nrows) * pixel_height  # Note: pixel_height is negative
        return lon_1d, lat_1d
    elif dem == 2:
        lon_2d, lat_2d = np.meshgrid(
            top_left_x + np.arange(ncols) * pixel_width,
            top_left_y + np.arange(nrows) * pixel_height
        )
        return lon_2d, lat_2d
    else:
        raise ValueError("Invalid dem value: must be 1 or 2.")
    
'''
# Example usage with .tif file
# 1. Calculate lon_1d and lat_1d
tif_path = os.path.join(cpuserver_data_FP, 'HLS/Korea/L30/20130411/B02_merged_20130411.tif')
lon_1d, lat_1d = get_lon_lat(tif_path)

# 2. Calculate lon_2d and lat_2d
lon_2d, lat_2d = get_lon_lat(tif_path, dem=2)
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Function to reproject raster to a desired EPSG coordinate system
def reproject(input_path, output_path, epsg_code):
    """
    This function reprojects a raster file to a specified EPSG coordinate system.

    Parameters:
    - input_path: Path to the input raster file
    - output_path: Path to save the reprojected raster file
    - epsg_code: EPSG code of the desired coordinate system (example: 3857 -> WGS 84 Pseudo-Mercator)

    Last modified: None / None
    Created on: 2024-12-27
    Created by: Seongjun Lee
    """
    # Enable exception handling
    gdal.UseExceptions()

    # Open the input raster file
    dataset = gdal.Open(input_path)
    if dataset is None:
        raise FileNotFoundError(f"Input file {input_path} not found.")

    # Get the current projection
    source_srs = osr.SpatialReference()
    source_srs.ImportFromWkt(dataset.GetProjection())

    # Create the target spatial reference
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(epsg_code)

    # Create a coordinate transformation
    transform = osr.CoordinateTransformation(source_srs, target_srs)

    # Reproject the raster
    gdal.Warp(output_path, dataset, dstSRS=target_srs.ExportToWkt())

    # Close the dataset
    dataset = None

'''
# Example usage
input_path = '/YOUR_PATH/EPSG4326.tif'
output_path = '/YOUR_PATH/EPSG3857.tif'
epsg_code = 3857
reproject_raster(input_path, output_path, epsg_code)
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
def plot(raster_file,
                transform=None,
                bounds=None, 
                title=None, 
                label=None, 
                cmap=None, 
                vmin=None, 
                vmax=None, 
                tick_interval=None,
                grid=True):
    '''
    This function plots raster data.

    Last modified: None / None
    Created on: 2024-12-27
    Created by: Seongjun Lee
    '''
    # Read raster data
    if isinstance(raster_file, np.ndarray):  # Check if raster_file is a numpy array
        raster_data = raster_file
        if transform is None:
            raise ValueError("Transform must be provided when raster_data is a numpy array.")
    elif transform is None:  # Use raster file path
        raster_data, transform, raster_crs = read(raster_file)
        # Check if the coordinate reference system is EPSG:4326
        if raster_crs.to_string() != 'EPSG:4326':
            raise ValueError(f"This code is based on the WGS 84 coordinate system, "
                             f"please reproject to EPSG:4326 and try again. \n"
                             f"Current CRS: {raster_crs.to_string()}.")
    elif transform is not None and raster_file.endswith('.tif'):
        raise ValueError("You have provided both a raster file path and transform. \n"
                         "Please input either just the raster file path without transform or "
                         "provide the already loaded raster data along with its transform.")
    else:
        raster_data, transform = raster_file, transform

    # Plot
    plt.figure(figsize=(12, 8), dpi=600)

    # Set default values if not provided
    cmap = cmap if cmap is not None else 'viridis'
    vmin = vmin if vmin is not None else np.nanmin(raster_data)
    vmax = vmax if vmax is not None else np.nanmax(raster_data)

    # Display the raster data
    plt.imshow(raster_data, extent=(
        transform[2], transform[2] + transform[0] * raster_data.shape[1],
        transform[5] + transform[4] * raster_data.shape[0], transform[5]),
        cmap=cmap, vmin=vmin, vmax=vmax, interpolation='none')

    # Set the limits based on the provided extent
    if bounds is not None:
        plt.xlim(bounds[0], bounds[1])
        plt.ylim(bounds[2], bounds[3]) 

    # Title
    plt.title(title if title is not None else 'Raster Data Visualization', fontsize=20, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(label=label if label is not None else 'Value', orientation='vertical')
    cbar.ax.set_ylabel(label if label is not None else 'Value', fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    # X/Y Ticks
    if tick_interval is not None:
        x_ticks = np.arange(bounds[0], bounds[1] + tick_interval, tick_interval)
        y_ticks = np.arange(bounds[2], bounds[3] + tick_interval, tick_interval)
        plt.xticks(x_ticks, fontsize=13, rotation=45)
        plt.yticks(y_ticks, fontsize=13)
    else:
        plt.xticks(fontsize=13, rotation=45)
        plt.yticks(fontsize=13)

    # Grid
    if grid:
        plt.grid(True, linestyle='--', color='gray', alpha=0.5)

    # Show the plot
    plt.show()

'''
# Example usage
1. If you want to directly load raster file
raster_file_path = '/home/seongjun/water_quality/obs_data/pekel_data/pekel_2024.tif'
plot_raster(raster_file=raster_file_path,
            bounds=[127.0, 127.5, 35.0, 35.5],
            title='Pekel Data',
            label='Percentage of Water Occurence (%)',
            cmap='Blues',
            vmin=0,
            vmax=100,
            tick_interval=1)

2. If you already have raster array and transform array
plot_raster(raster_file=pekel_data,
            transform=pekel_transform,
            bounds=[127.0, 127.5, 35.0, 35.5],
            title='Pekel Data',
            label='Percentage of Water Occurence (%)',
            cmap='Blues',
            vmin=0,
            vmax=100,
            tick_interval=1)
'''
#---------------------------------------------------------------------------------------------------#