import sys
import uproot # root in python
import h5py # h5 in python
import numpy as np
from timeit import default_timer as timer


# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# converts given .root file to a .h5 file
#
# usage: python3 convert_root_to_h5.py [INPUT_ROOT_FILE] [OUTPUT_H5_DESTINATION]


MAX_NUMBER_OF_EVENTS = -1 # set to -1 to use all events



# main method 
def main():
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for conversion was given, exiting")
		exit()
	if(len(sys.argv)<3):
		print("Error: no .h5 file destination was given, exiting")
		exit()

	# open .root file
	root_file = uproot.open(sys.argv[1])
	
	# create .h5 file
	h5_file = h5py.File(sys.argv[2], 'w')
	
	# fill .h5 file with .root information
	fill_h5_from_root(h5_file, root_file)
	
	# close .h5 file, .root cannot be closed
	h5_file.close()



def fill_h5_from_root(h5_file, root_file):
	# start conversion timer 
	timer_start_conversion = timer()


	# get number of jets and events, apply limit on number of events
	number_of_jets = root_file['matched/number_of_jets'].array()
	number_of_events = len(number_of_jets)
	max_number_of_jets = max(number_of_jets)
	if MAX_NUMBER_OF_EVENTS>-1 and number_of_events>MAX_NUMBER_OF_EVENTS:
		number_of_events = MAX_NUMBER_OF_EVENTS


	# create INPUTS group for SPANet (reconstruction information e.g. jets, missing transverse energy, btagging, mask,...)
	print("creating group 'INPUTS'")
	input_group = h5_file.create_group("INPUTS")


	# create SOURCE subgroup for SPANet (source information per event)
	source_group = input_group.create_group("Source")


	# set fixed out array length, use MASK to mask which jets are actually in an event
	spanet_source_dimension = (number_of_events, max_number_of_jets)
	mask = source_group.create_dataset("MASK", spanet_source_dimension, dtype=bool)
	
	
	# add custom information per event for SPANet to use in training
	# here branches can be customized if desired
	jet_e = source_group.create_dataset("energy", spanet_source_dimension, dtype=np.float32)
	jet_pt = source_group.create_dataset("pt", spanet_source_dimension, dtype=np.float32)
	jet_eta = source_group.create_dataset("eta", spanet_source_dimension, dtype=np.float32)
	jet_phi = source_group.create_dataset("phi", spanet_source_dimension, dtype=np.float32)
	jet_match_mask = source_group.create_dataset("match_mask", spanet_source_dimension, dtype=int) # integer in bit representation yields matched particle
	
	
	# prepare root_files for fast access
	root_energy = root_file['matched/jet_e_NOSYS'].array()
	root_pt = root_file['matched/jet_pt_NOSYS'].array()
	root_eta = root_file['matched/jet_eta'].array()
	root_phi = root_file['matched/jet_phi'].array()
	root_match_mask = root_file['matched/jet_final_match_mask'].array()


	# start timer for INPUT steps and prepare variables for printing
	step = 1/number_of_events
	ratio = 0.0
	displayed_ratio= 0
	timer_start_input = timer()


	# loop for INPUTS group
	for i in range(number_of_events):
		for j in range(number_of_jets[i]):
			# for each event fill dataset per jet j in that specific event i
			# it is important that both indicies are used simultanously e.g. [i,j] instead of [i][j]
			# otherwise data is set in copied array and is discarded!
			mask[i,j] = True # mask is true for every jet that is filled


			# here custom branches must be filled
			jet_e[i,j] = root_energy[i,j]
			jet_pt[i,j] = root_pt[i,j]
			jet_eta[i,j] = root_eta[i,j]
			jet_phi[i,j] = root_phi[i,j]
			jet_match_mask[i,j] = root_match_mask[i,j]


		# print progress every 10%
		ratio+= step #1/number_of_events is dont
		if ratio>=0.1: 
			displayed_ratio+= 1
			ratio-= 0.1
			timer_partial = timer()
			print("  >", 10*displayed_ratio, "% (estimated remaining time: ", round((timer_partial-timer_start_input)/displayed_ratio*(10-displayed_ratio)), "s)")
	print()	
	
	
	# create TARGET group for SPANET (truth information about matching e.g. jet indicies of truth object)
	# here different outputs can be defined, if desired
	print("creating group 'TARGET'")
	target_group = h5_file.create_group("TARGETS")

	# set fixed out array length with one entry per event
	spanet_target_dimension = (number_of_events)

	# create t1 subgroup for truth infromation of first t
	t1_group = target_group.create_group("t1")
	b1 = t1_group.create_dataset("b", spanet_target_dimension, dtype=int)
	q1_1 = t1_group.create_dataset("q1", spanet_target_dimension, dtype=int)
	q1_2 = t1_group.create_dataset("q2", spanet_target_dimension, dtype=int)
	
	# create t2 subgroup for truth infromation of second t
	t2_group = target_group.create_group("t2")
	b2 = t2_group.create_dataset("b", spanet_target_dimension, dtype=int)
	q2_1 = t2_group.create_dataset("q1", spanet_target_dimension, dtype=int)
	q2_2 = t2_group.create_dataset("q2", spanet_target_dimension, dtype=int)
	
	# create HW subgroup for truth infromation of the hadronic decaying W from H
	HW_group = target_group.create_group("HW")
	HW_1 = HW_group.create_dataset("q1", spanet_target_dimension, dtype=int)
	HW_2 = HW_group.create_dataset("q2", spanet_target_dimension, dtype=int)

	# prepare root_files for fast access
	indicies = root_file["matched/jet_to_object_indicies_fixed"].array()

	
	# timer for target and reset variables for printing
	ratio = 0.0
	displayed_ratio = 0
	start_group = timer()
		

	# loop for TARGET group
	for i in range(number_of_events):
		# fill datasets with custom data
		# important is fixed order within root file structured set in c++ code
		b1[i] = indicies[i][0]
		q1_1[i] = indicies[i][1]
		q1_2[i] = indicies[i][2]
		b2[i] = indicies[i][3]
		q2_1[i] = indicies[i][4]
		q2_2[i] = indicies[i][5]
		HW_1[i] = indicies[i][6]
		HW_2[i] = indicies[i][7]


		# print progress
		ratio+= step
		if ratio>=0.1:
			displayed_ratio+= 1
			ratio-= 0.1
			partial = timer()
			print("  >", 10*displayed_ratio, "% (estimated remaining time: ", round((partial-start_group)/displayed_ratio*(10-displayed_ratio)), "s)")
	

	# print time for conversion
	timer_end_conversion = timer()
	print()
	print("Time needed: ", round(timer_end_conversion-timer_start_conversion), "s for ", number_of_events, "events")



if __name__ == '__main__':
	main()
