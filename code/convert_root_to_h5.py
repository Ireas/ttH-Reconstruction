import sys
import uproot # root in python
import h5py # h5 in python
import awkward # working with irregular arrays
import numpy as np

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# converts .root file to a .h5 file at same destination and with same name
# todo: implement creating nested dataset/groups 
# todo: improve speed (mutlithreading, non-uniform data, ...)
# todo: take input from argument as print function
# todo: make print root file own script

ROOT_FILE_DESTINATION = "../output/v1/output_data.root"
H5_FILE_DESTINATION = "../output/v1/output_data.h5"


def main():
	root_file = uproot.open(ROOT_FILE_DESTINATION)
	convert_to_h5(root_file)


def convert_to_h5(root_file):
	h5_file = h5py.File(H5_FILE_DESTINATION, 'w')

	maximum_number_of_events = root_file['other']['max_number_of_events'].array()[0]
	maximum_number_of_jets = root_file['other']['max_number_of_jets'].array()[0]
	mask_is_set = False	

	source_group = h5_file.create_group("SOURCE")
	mask = source_group.create_dataset("MASK", (maximum_number_of_events,maximum_number_of_jets), dtype=bool)
	
	print("group: SOURCE")
	
	for branch in root_file['source'].keys():
		print("dataset: " + branch)
		content = root_file['source'][branch].array()

		dt = h5py.special_dtype(vlen=np.dtype(float)) # use variable length float array because jet number varies per event
		dataset = source_group.create_dataset(branch, (len(content),), dtype=dt)
		for i, entry in enumerate(content): # add each event one-by-one
			dataset[i] = np.array(entry, dtype=float)
			if not mask_is_set:
				for j in range(len(entry)):	
					mask[i][j] = True		
		
		mask_is_set = True
		print("successfully written: " + branch + " to SOURCE")
	
	print("group: INDICIES")
	
	indicies_group = h5_file.create_group("INDICIES")
	for branch in root_file['indicies'].keys():
		print("dataset: " + branch)
		content = root_file['indicies'][branch].array()
		indicies_group.create_dataset(branch, data=content)
		print("successfully written: " + branch + " to INDICIES")
		
	print("Indicies successful")
	
	h5_file.close()


if __name__ == '__main__':
	main()
