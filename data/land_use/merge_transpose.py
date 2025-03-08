# Transpose latitude and longitude of land use dataset (mixed up originally)
# Merge layers of all rainfed crops and irrigation crops

import xarray as xr
import numpy as np
import os
from pathlib import Path

def merge_netcdf_variables(input_file, output_file):
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Open dataset
    ds = xr.open_dataset(input_file).transpose()
    
    # Calculate sums of odd and even variables (15,17,...,29 and 16,18,...,30)
    odd_vars = ["PFT" + str(i) for i in range(15, 30, 2)]
    even_vars = ["PFT" + str(i) for i in range(16, 31, 2)]
    
    # Create new variables by summing
    var33 = sum(ds[var] for var in odd_vars)
    var34 = sum(ds[var] for var in even_vars)
    
    # Create new dataset with variables 1-14, 31-33, and the new sums
    keep_vars = ["PFT" + str(i) for i in list(range(0, 15)) + list(range(31, 33))]
    new_ds = xr.Dataset()
    
    # Copy kept variables
    for var in keep_vars:
        new_ds[var] = ds[var]
    
    # Add new summed variables
    new_ds['Rainfed'] = var33
    new_ds['Irrigated'] = var34
    
    # Preserve coordinates and attributes
    new_ds.attrs = ds.attrs
    
    # Save and close
    new_ds.to_netcdf(output_file)
    ds.close()
    new_ds.close()

def process_directory(input_dir, output_dir):
    # Convert to Path objects for easier manipulation
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Process all .nc files in the directory and its subdirectories
    for input_file in input_path.rglob('*.nc'):
        # Skip macOS metadata files
        if input_file.name.startswith('._'):
            print(f"Skipping macOS metadata file: {input_file}")
            continue
            
        # Skip .DS_Store files
        if input_file.name == '.DS_Store':
            continue
        
        # Create corresponding output path
        relative_path = input_file.relative_to(input_path)
        output_file = output_path / relative_path
        
        print(f"\nProcessing: {input_file} -> {output_file}")
        try:
            merge_netcdf_variables(str(input_file), str(output_file))
            print(f"Successfully processed {input_file}")
        except Exception as e:
            print(f"Error processing {input_file}: {str(e)}")
            continue

if __name__ == "__main__":
    process_directory("original", "merged")
