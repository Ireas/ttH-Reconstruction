import sys
import uproot # root in python
import h5py # h5 in python
import awkward # working with irregular arrays
import numpy as np

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# converts .root file to a .h5 file at same destination and with same name
# todo: improve speed (mutlithreading, non-uniform data, ...)
# todo: take input from argument as print function
# todo: add timestamps to output estimated time

ROOT_FILE_DESTINATION = "../output/rdataframes_output.root"
H5_FILE_DESTINATION = "../output/output.h5"


def main():
	root_file = uproot.open(ROOT_FILE_DESTINATION)
	convert_to_h5(root_file)



def convert_to_h5(root_file):
	h5_file = h5py.File(H5_FILE_DESTINATION, 'w')

	number_of_events = len(root_file['matched/nJets'].array())
	max_number_of_jets = max(root_file['matched/nJets'].array())
	print(number_of_events)
	print(max_number_of_jets)

	input_group = h5_file.create_group("INPUTS")
	source_group = input_group.create_group("SOURCE")
	mask = source_group.create_dataset("MASK", (number_of_events,max_number_of_jets), dtype=bool)
	jet_match_mask = source_group.create_dataset("JET_MATCH_MASK", (number_of_events,max_number_of_jets), dtype=int)
	jet_pt = source_group.create_dataset("JET_PT", (number_of_events,max_number_of_jets), dtype=float)
	jet_eta = source_group.create_dataset("JET_ETA", (number_of_events,max_number_of_jets), dtype=float)
	jet_phi = source_group.create_dataset("JET_PHI", (number_of_events,max_number_of_jets), dtype=float)
	jet_e = source_group.create_dataset("JET_E", (number_of_events,max_number_of_jets), dtype=float)

		
	print("group: SOURCE")
	for i in range(10):
		number_of_jets = root_file['matched/nJets'].array()[i]
		for j in range(number_of_jets):
			# important that both indicies are used simultanously, otherwise data is set in copied array and is discarded!
			mask[i,j] = True
			jet_pt[i,j] = root_file['matched/jet_pt_NOSYS'].array()[i,j]
			jet_eta[i,j] = root_file['matched/jet_eta'].array()[i,j]
			jet_phi[i,j] = root_file['matched/jet_phi'].array()[i,j]
			jet_e[i,j] = root_file['matched/jet_e_NOSYS'].array()[i,j]
			jet_match_mask[i,j] = root_file['matched/jet_match_mask'].array()[i,j]

	
	print("group: TARGET")
	target_group = h5_file.create_group("TARGETS")
	t1_group = target_group.create_group("t1")
	b1 = t1_group.create_dataset("b", (number_of_events), dtype=int)
	q1_1 = t1_group.create_dataset("q1", (number_of_events), dtype=int)
	q1_2 = t1_group.create_dataset("q2", (number_of_events), dtype=int)
	
	t2_group = target_group.create_group("t2")
	b2 = t2_group.create_dataset("b", (number_of_events), dtype=int)
	q2_1 = t2_group.create_dataset("q1", (number_of_events), dtype=int)
	q2_2 = t2_group.create_dataset("q2", (number_of_events), dtype=int)
	
	
	h5_file.close()
	exit()
	
#for branch in root_file['source'].keys():
#		print("dataset: " + branch)
#		content = root_file['source'][branch].array()
#
#		dt = h5py.special_dtype(vlen=np.dtype(float)) # use variable length float array because jet number varies per event
#		dataset = source_group.create_dataset(branch, (len(content),), dtype=dt)
#		for i, entry in enumerate(content): # add each event one-by-one
#			dataset[i] = np.array(entry, dtype=float)
#			if not mask_is_set:
#				for j in range(len(entry)):	
#					mask[i][j] = True		
#		
#		mask_is_set = True
#		print("successfully written: " + branch + " to SOURCE")
	


if __name__ == '__main__':
	main()
