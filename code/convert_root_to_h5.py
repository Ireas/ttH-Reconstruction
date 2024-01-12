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


MAX_NUMBER_OF_EVENTS = 1000 # set to -1 to use all events


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
	# open file
	h5_file = h5py.File(destination, 'w')
	

	# start timer 
	timer_start_conversion = time.time()


	# get event numbers
	number_of_jets = root_file['matched/nJets'].array()
	number_of_events = len(number_of_jets)
	max_number_of_jets = max(number_of_jets)
	if MAX_NUMBER_OF_EVENTS>-1 and number_of_events>MAX_NUMBER_OF_EVENTS:
		number_of_events = MAX_NUMBER_OF_EVENTS


	# create INPUTS group
	print("creating group 'INPUTS'")
	input_group = h5_file.create_group("INPUTS")
	source_group = input_group.create_group("Source")
	mask = source_group.create_dataset("MASK", (number_of_events,max_number_of_jets), dtype=bool)
	jet_e = source_group.create_dataset("energy", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_pt = source_group.create_dataset("pt", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_eta = source_group.create_dataset("eta", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_phi = source_group.create_dataset("phi", (number_of_events,max_number_of_jets), dtype=np.float32)
	jet_match_mask = source_group.create_dataset("match_mask", (number_of_events,max_number_of_jets), dtype=int)
	

	# timer for input
	timer_start_input = time.time()
	step = 1/number_of_events
	ratio = 0.0
	displayed_percentage = 0


	# prepare root_files for fast acces
	root_energy = root_file['matched/jet_e_NOSYS'].array()
	root_pt = root_file['matched/jet_pt_NOSYS'].array()
	root_eta = root_file['matched/jet_eta'].array()
	root_phi = root_file['matched/jet_phi'].array()
	root_match_mask = root_file['matched/jet_match_mask'].array()
	


	# loop for INPUTS group
	for i in range(number_of_events):
		# fill dataset
		for j in range(number_of_jets[i]):
			# important that both indicies are used simultanously, otherwise data is set in copied array and is discarded!
			mask[i,j] = True
			jet_e[i,j] = root_energy[i,j]
			jet_pt[i,j] = root_pt[i,j]
			jet_eta[i,j] = root_eta[i,j]
			jet_phi[i,j] = root_phi[i,j]
			jet_match_mask[i,j] = root_match_mask[i,j]


		# show progress
		ratio+= step
		if ratio>=0.1:
			displayed_percentage+= 1
			ratio-= 0.1
			timer_partial = time.time()
			print("  >", 10*displayed_percentage, "% (estimated remaining time: ", round((timer_partial-timer_start_input)/displayed_percentage*(10-displayed_percentage)), "s)")
	
	
	# create TARGET group
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
	

	# timer for target
	start_group = time.time()
	ratio = 0.0
	displayed_percentage = 0
		

	# prepare root_files for fast access
	indicies = root_file["matched/indicies"].array()


	# loop for TARGET group
	for i in range(number_of_events):
		# fill datasets
		b1[i] = indicies[i][0]
		q1_1[i] = indicies[i][1]
		q1_2[i] = indicies[i][2]
		b2[i] = indicies[i][3]
		q2_1[i] = indicies[i][4]
		q2_2[i] = indicies[i][5]


		# show progress
		ratio+= step
		if ratio>=0.1:
			displayed_percentage+= 1
			ratio-= 0.1
			partial = time.time()
			print("  >", 10*displayed_percentage, "% (estimated remaining time: ", round((partial-start_group)/displayed_percentage*(10-displayed_percentage)), "s)")
	

	# print time for conversion
	timer_end_conversion = time.time()
	print()
	print("Time needed: ", round(timer_end_conversion-timer_start_conversion), "s for ", number_of_events, "events")

		
	# close file
	h5_file.close()


if __name__ == '__main__':
	main()
