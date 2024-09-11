import h5py
import os
import numpy as np

# Constants
INPUT_FILE = "/media/ireas/Data/v6/merged_h5/all_8+j_1529106e.h5"
OUTPUT_PATH = '/media/ireas/Data/v6/split_h5/'
SPLIT_RATIO = 0.1  # e.g., 10% split


# Ensure output directory exists
os.makedirs(OUTPUT_PATH, exist_ok=True)


# Split filenames based on the ratio
base_name = os.path.basename(INPUT_FILE)
name, ext = os.path.splitext(base_name)
file1_name = f"{name}_{int(SPLIT_RATIO*100)}split{ext}"
file2_name = f"{name}_{int((1-SPLIT_RATIO)*100)}split{ext}"


# File paths
file1_path = os.path.join(OUTPUT_PATH, file1_name)
file2_path = os.path.join(OUTPUT_PATH, file2_name)


def split_data(input_group, output_group1, output_group2, split_ratio):
    for key in input_group.keys():
        item = input_group[key]
        print(f"  splitting key {key}")
        
        if isinstance(item, h5py.Dataset):
            data = item[:]
            total_length = len(data)
            indices = np.arange(total_length)
            
            # Determine step size
            step_size = int(1 / split_ratio)
            
            # Select every i-th entry
            selected_indices = indices[::step_size]
            remaining_indices = np.setdiff1d(indices, selected_indices)
            
            output_group1.create_dataset(key, data=data[selected_indices])
            output_group2.create_dataset(key, data=data[remaining_indices])
        
        elif isinstance(item, h5py.Group):
            grp1 = output_group1.create_group(key)
            grp2 = output_group2.create_group(key)
            split_data(item, grp1, grp2, split_ratio)


# Load input file
print(f"Spltting input file at {INPUT_FILE}")

with h5py.File(INPUT_FILE, 'r') as input_file:
    # Create output files
    with h5py.File(file1_path, 'w') as file1, h5py.File(file2_path, 'w') as file2:
        split_data(input_file, file1, file2, SPLIT_RATIO)