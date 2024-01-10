import sys
import uproot # root in python
import h5py # h5 in python
import awkward # working with irregular arrays
import numpy as np
import time

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# converts .root file to a .h5 file at same destination and with same name
# TODO: improve speed (mutlithreading, non-uniform data, ...)
# TODO: CORRECT INDEXING FOR TARGETS

TESTING_MODE = True


def main():
	if(len(sys.argv)<2):
		print("Error: no file for conversion is given, exiting")
		exit()
	if(len(sys.argv)<3):
		print("Error: no file for conversion is given, exiting")
		exit()
	root_file = uproot.open(sys.argv[1])
	convert_to_h5(root_file, sys.argv[2])



def convert_to_h5(root_file, destination):
	h5_file = h5py.File(destination, 'w')
	
	start = time.time()

	number_of_events = len(root_file['matched/nJets'].array())
	max_number_of_jets = max(root_file['matched/nJets'].array())
	if TESTING_MODE and number_of_events>1000:
		number_of_events = 1000

	
	print("creating group 'INPUTS'")
	input_group = h5_file.create_group("INPUTS")
	source_group = input_group.create_group("Source")
	mask = source_group.create_dataset("MASK", (number_of_events,max_number_of_jets), dtype=bool)
	jet_e = source_group.create_dataset("mass", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_pt = source_group.create_dataset("pt", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_eta = source_group.create_dataset("eta", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_phi = source_group.create_dataset("phi", (number_of_events,max_number_of_jets), dtype=np.float32)
	#jet_match_mask = source_group.create_dataset("JET_MATCH_MASK", (number_of_events,max_number_of_jets), dtype=int)
	btag = source_group.create_dataset("btag", (number_of_events,max_number_of_jets), dtype=int)
		
	start_group = time.time()
	step = 1/number_of_events
	ratio = 0.0
	displayed_percentage = 0

	for i in range(number_of_events):
		number_of_jets = root_file['matched/nJets'].array()[i]
		for j in range(number_of_jets):
			# important that both indicies are used simultanously, otherwise data is set in copied array and is discarded!
			mask[i,j] = True
			jet_e[i,j] = root_file['matched/jet_e_NOSYS'].array()[i,j]
			jet_pt[i,j] = root_file['matched/jet_pt_NOSYS'].array()[i,j]
			jet_eta[i,j] = root_file['matched/jet_eta'].array()[i,j]
			jet_phi[i,j] = root_file['matched/jet_phi'].array()[i,j]
			#jet_match_mask[i,j] = root_file['matched/jet_match_mask'].array()[i,j]
			btag[i,j] = 0

		ratio+= step
		if ratio>=0.1:
			displayed_percentage+= 1
			ratio-= 0.1
			partial = time.time()
			print("  >", 10*displayed_percentage, "% (estimated remaining time: ", round((partial-start_group)/displayed_percentage*(10-displayed_percentage)), "s)")
	
	
	print()	
	print("creating group 'TARGET'")
	target_group = h5_file.create_group("TARGETS")
	t1_group = target_group.create_group("t1")
	b1 = t1_group.create_dataset("b", (number_of_events), dtype=int)
	q1_1 = t1_group.create_dataset("q1", (number_of_events), dtype=int)
	q1_2 = t1_group.create_dataset("q2", (number_of_events), dtype=int)
	
	t2_group = target_group.create_group("t2")
	b2 = t2_group.create_dataset("b", (number_of_events), dtype=int)
	q2_1 = t2_group.create_dataset("q1", (number_of_events), dtype=int)
	q2_2 = t2_group.create_dataset("q2", (number_of_events), dtype=int)
	
	start_group = time.time()
	ratio = 0.0
	displayed_percentage = 0
	# dummy variable	
	for i in range(number_of_events):
		b1[i] = 0
		q1_1[i] = 1
		q1_2[i] = 2
		b2[i] = 3
		q2_1[i] = 4
		q2_2[i] = 5

		ratio+= step
		if ratio>=0.1:
			displayed_percentage+= 1
			ratio-= 0.1
			partial = time.time()
			print("  >", 10*displayed_percentage, "% (estimated remaining time: ", round((partial-start_group)/displayed_percentage*(10-displayed_percentage)), "s)")
	
	h5_file.close()
	end = time.time()
	print()
	print("Time needed: ", round(end-start), "s for ", number_of_events, "events")
		


if __name__ == '__main__':
	main()
