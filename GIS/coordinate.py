#---------------------------------------------------------------------------------------------------#
# Import libraries
import pandas as pd
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Convert DMS to DD
def dms_to_dd(dms):
    """Convert a DMS string to decimal degrees."""
    # Check if the input is already in decimal format
    if isinstance(dms, str) and '°' in dms:
        dms_parts = dms.replace('°', ' ').replace("'", ' ').replace('"', '').split()
        degrees = float(dms_parts[0])
        minutes = float(dms_parts[1])
        seconds = float(dms_parts[2])
        return degrees + (minutes / 60) + (seconds / 3600)
    else:
        return float(dms)  # Return as float if already in decimal format

'''
# Usage example
dms_example = "37° 43' 55.0\""
decimal_degrees = dms_to_dd(dms_example)
print(f"The decimal degrees: {decimal_degrees}")
'''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# Convert DMS to DD in a CSV file
def convert_dms_to_dd(file_path, save=False, save_path=None):
    """
    Convert DMS (Degrees, Minutes, Seconds) coordinates in a CSV file to DD (decimal degrees).

    Parameters:
    - file_path (str): The path to the input CSV file containing DMS coordinates.
    - save (bool): If True, save the converted DataFrame to a new CSV file.
    - save_path (str): The path to save the converted CSV file if save is True.

    Returns:
    - pd.DataFrame: A DataFrame containing the original data with additional decimal degree columns for longitude and latitude.

    Note:
    - The input CSV must have columns named 'lon' and 'lat' for longitude and latitude coordinates, respectively.
    """
    import pandas as pd
    
    # Read the CSV file
    df = pd.read_csv(file_path, encoding='utf-8')
    
    # Create new columns for the decimal coordinates
    df['lon'] = df['lon'].apply(dms_to_dd).astype(float)
    df['lat'] = df['lat'].apply(dms_to_dd).astype(float)

    if save:
        df.to_csv(save_path, encoding='utf-8-sig', index=False)
        print(f"'{save_path}' is saved")

    return df

'''
# Apply the function to the uploaded file
file_path = '/home/seongjun/Downloads/obs_coords_2025.csv'
save_path = '/home/seongjun/Downloads/obs_coords_2025_new.csv'
converted_data = convert_dms_to_dd(file_path, save=True, save_path=save_path)
'''