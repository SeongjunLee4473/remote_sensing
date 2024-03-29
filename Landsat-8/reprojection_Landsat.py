# reprojection_landsat.py

import os
from osgeo import gdal
gdal.UseExceptions()

def reprojection(input_directory, output_directory):
    # Get all TIF files in the input directory
    tif_files = [file for file in os.listdir(input_directory) if file.endswith('.TIF')]

    # Create the output directory if it does not exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Perform geometric correction and output for each TIF file
    for tif_file in tif_files:
        input_file = os.path.join(input_directory, tif_file)
        output_file = os.path.join(output_directory, tif_file)

        # Open the image file
        input_dataset = gdal.Open(input_file)

        # Create a geometrically corrected image
        gdal.Warp(output_file, input_dataset, dstSRS='EPSG:4326')  # EPSG:4326은 WGS84 좌표계를 의미합니다.

        # Clean up dataset and memory
        input_dataset = None

    print("Geometric correction is complete.")