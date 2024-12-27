#---------------------------------------------------------------------------------------------------#
# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from osgeo import gdal, osr
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
# Function to reproject raster to a desired EPSG coordinate system
def reproject(input_path, output_path, epsg_code):
    """
    This function reprojects a raster file to a specified EPSG coordinate system.

    Parameters:
    - input_path: Path to the input raster file
    - output_path: Path to save the reprojected raster file
    - epsg_code: EPSG code of the desired coordinate system (example: 3857 -> WGS 84 Pseudo-Mercator)
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
        raster_data, transform, raster_crs = read_raster(raster_file)
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
    plt.figure(figsize=(12, 8))

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