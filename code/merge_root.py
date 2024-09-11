import os
import sys
import uproot
import numpy as np
import awkward as ak
from timeit import default_timer as timer


# constants
INPUT_FOLDER = "all_5+j_truncated_50_ttHWW_amplified/"
INPUT_PATH = "/media/ireas/Data/v6/matched/"
OUTPUT_PATH = "/media/ireas/Data/v6/merged_root/"


# main method 
def main():
	# access input
	input_files = []
	input_directory = INPUT_PATH+INPUT_FOLDER
	print(f"Accessing input files in {input_directory}")
	for input_destination in os.listdir(input_directory):
		print(f" > {input_destination}")
		input_files.append(input_destination)
	print("")
	print(f"{len(input_files)} files found!")

	
	# List to store the data from each file
	all_data_dict = {}


	# Loop over each file and read the data
	for (i, file) in enumerate(input_files):
		with uproot.open(INPUT_PATH+INPUT_FOLDER+file) as root_file:
			new_data_dict = root_file["matched"].arrays(library="np")
			
			# check for empty files
			if(len(new_data_dict.keys())==0):
				continue
			
			number_of_events = len(new_data_dict["eventNumber"])
			print(f" ({i+1}/{len(input_files)})> {file}: including {number_of_events} events")

    		# Combine keys from both dictionaries
			all_keys = set(all_data_dict.keys()).union(set(new_data_dict.keys()))

			for key in all_keys:
				# If key is in both dictionaries, append the values
				if key in all_data_dict and key in new_data_dict:
					all_data_dict[key] = np.append(all_data_dict[key], new_data_dict[key])
				
				# If key is only in dict_a
				elif key in all_data_dict:
					all_data_dict[key] = all_data_dict[key]
				
				# If key is only in dict_b
				elif key in new_data_dict:
					all_data_dict[key] = new_data_dict[key]
				

	number_of_events = len(all_data_dict["eventNumber"])
	

	# Write the merged data to a new ROOT file
	print(f"Recreating Output File {input_directory} with {number_of_events} events")
	with uproot.recreate(f"{OUTPUT_PATH}{INPUT_FOLDER[:-1]}_{number_of_events}e_merged.root") as merged_root_file:
		merged_root_file["matched"] = all_data_dict


if __name__ == '__main__':
	main()
