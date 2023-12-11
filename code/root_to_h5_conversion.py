import uproot # root in python
import h5py # h5 in python
import awkward # working with irregular arrays
import numpy as np

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# prints out the .root file at fixed destination and converts it to a flat .h5 file
# todo: implement creating nested dataset/groups 
# todo: improve speed (mutlithreading, non-uniform data, ...)

ROOT_FILE_DESTINATION = "../output/v1/output_data.root"
H5_FILE_DESTINATION = "../output/v1/output_data.h5"



def main():
	root_file = uproot.open(ROOT_FILE_DESTINATION)
	print_root_file(root_file)
	convert_to_h5(root_file)



def print_root_file(root_file):
	print("=====  ROOT FILE STRUCTURE  =====")
	for tree in root_file.keys():
		print("TTree: " + tree)
		for branch in root_file[tree].keys():
			content = root_file[tree][branch].array()
			print("\tTBranch: " + branch + ": " + str(content))



def convert_to_h5(root_file):
	h5_file = h5py.File(H5_FILE_DESTINATION, 'w')

	for tree in root_file.keys():
		maximum_number_of_jets = max(root_file[tree]['number_of_jets'].array())
		for branch in root_file[tree].keys():
			content = root_file[tree][branch].array()
			
			if type(content[0])==awkward.highlevel.Array: # manual treatment of non-uniform arrays
				dt = h5py.special_dtype(vlen=np.dtype(float)) # use variable length float array because jet number varies per event
				dataset = h5_file.create_dataset(branch, (len(content),), dtype=dt)
				for i, entries in enumerate(content): # add each event one-by-one
					dataset[i] = np.array(entries, dtype=float)
			else:
				h5_file.create_dataset(branch, data=content)
			print("successfully written: " + branch)

	h5_file.close()



if __name__ == '__main__':
	main()
