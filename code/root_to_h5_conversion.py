import uproot # root in python
import h5py # h5 in python
import awkward # working with irregular arrays
import numpy as np

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# prints out the .root file at fixed destination and converts it to a flat .h5 file
# todo: implement creating nested dataset/groups 
# todo: improve speed (mutlithreading, non-uniform data, ...)
# todo: take input from argument as print function
# todo: make print root file own script

ROOT_FILE_DESTINATION = "../output/v1/output_data.root"
H5_FILE_DESTINATION = "../output/v1/output_data.h5"



def main():
	root_file = uproot.open(ROOT_FILE_DESTINATION)
	#print_root_file(root_file)
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

	maximum_number_of_jets = root_file['other']['max_number_of_jets'].array()[0]
	
	source_group = h5_file.create_group("SOURCE")
	for branch in root_file['source'].keys():
		content = root_file['source'][branch].array()
		padded_content = np.zeros((len(content),maximum_number_of_jets), dtype=float)
		for i, subarray in enumerate(content):
			padded_content[i, :len(subarray)] = subarray
	
		source_group.create_dataset(branch, data=padded_content)
		print("successfully written: " + branch)
		
	
	h5_file.close()

#	for tree in root_file.keys():
#		for branch in root_file[tree].keys():
#			content = root_file[tree][branch].array()
#			
#			if type(content[0])==awkward.highlevel.Array: # manual treatment of non-uniform arrays
#				dt = h5py.special_dtype(vlen=np.dtype(float)) # use variable length float array because jet number varies per event
#				dataset = h5_file.create_dataset(branch, (len(content),), dtype=dt)
#				for i, entries in enumerate(content): # add each event one-by-one
#					dataset[i] = np.array(entries, dtype=float)
#			else:
#				h5_file.create_dataset(branch, data=content)




if __name__ == '__main__':
	main()
