import os
import os.path
import sys
import uproot # root in python
import h5py # h5 in python
import numpy as np
from timeit import default_timer as timer


# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
CURRENT_FOLDER = "all_8+j/"

INPUT_PATH = "/media/ireas/Data/v6/matched/"
OUTPUT_PATH = "/media/ireas/Data/v6/converted/"


# main method 
def main():
	directory = os.fsencode(INPUT_PATH+CURRENT_FOLDER)
    
	for (i,file) in enumerate(os.listdir(directory)):
		filename = os.fsdecode(file)

		target_file = OUTPUT_PATH+CURRENT_FOLDER+"converted"+filename[7:-5]+".h5"

		print(f" ({i+1}/{len(os.listdir(directory))}) - converting {INPUT_PATH+CURRENT_FOLDER+filename} -> converted{filename[7:-5]}.h5")
		if os.path.exists(target_file):
			print(" > target file exits already, skipping...")
			print()
			continue
		

		# fill h5 file with root file 
		with uproot.open(INPUT_PATH+CURRENT_FOLDER+filename) as root_file:			
			# check if root_file has any entries
			if not check_root_integrety(root_file):
				print(" > no entries, exiting....")
				print()
				continue
			
			# create new h5 file
			h5_file = h5py.File(target_file, 'w')
			fill_h5_from_root(h5_file, root_file)
			h5_file.close()
			print()


def check_root_integrety(root_file):
	if (len(root_file["matched"].keys())==0):
		return False

	return True


def fill_h5_from_root(h5_file, root_file):
	# get number of jets and events, apply limit on number of events
	number_of_jets = root_file['matched/number_of_jets'].array()
	number_of_events = len(number_of_jets)
	max_number_of_jets = max(number_of_jets)

	print(" > using " +str(number_of_events)+ " events!")
	print()


	# create INPUTS group for SPANet (reconstruction information e.g. jets, missing transverse energy, btagging, mask,...)
	print(" > creating 'INPUTS'")
	input_group = h5_file.create_group("INPUTS")


	# create SOURCE subgroup for SPANet (source information per event)
	source_group = input_group.create_group("Source")
	met_group = input_group.create_group("Met")


	# set fixed out array length, use MASK to mask which jets are actually in an event
	spanet_source_dimension = (number_of_events, max_number_of_jets)
	spanet_met_dimension = (number_of_events)
	mask = source_group.create_dataset("MASK", spanet_source_dimension, dtype=bool)
	
	
	# add custom information per event for SPANet to use in training
	# here branches can be customized if desired
	jet_e = source_group.create_dataset("energy", spanet_source_dimension, dtype=np.float32)
	jet_pt = source_group.create_dataset("pt", spanet_source_dimension, dtype=np.float32)
	jet_eta = source_group.create_dataset("eta", spanet_source_dimension, dtype=np.float32)
	jet_phi = source_group.create_dataset("phi", spanet_source_dimension, dtype=np.float32)
	jet_btag = source_group.create_dataset("btag", spanet_source_dimension, dtype=bool)
	
	# global
	global_met_value = met_group.create_dataset("value", spanet_met_dimension, dtype=np.float32)
	global_met_phi = met_group.create_dataset("phi", spanet_met_dimension, dtype=np.float32)

	
	
	# prepare root_files for fast access
	root_energy = root_file['matched/jet_e_NOSYS'].array()
	root_pt = root_file['matched/jet_pt_NOSYS'].array()
	root_eta = root_file['matched/jet_eta'].array()
	root_phi = root_file['matched/jet_phi'].array()
	root_btag = root_file['matched/jet_btag_85wp'].array()

	root_met_value = root_file['matched/reco_met_value'].array()
	root_met_phi = root_file['matched/reco_met_phi'].array()


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
			jet_btag[i,j] = root_btag[i,j]
		
		# add gloabl met variables
		global_met_value[i] = root_met_value[i]
		global_met_phi[i] = root_met_phi[i]
		


		# print progress every 10%
		ratio+= step #1/number_of_events is dont
		if ratio>=0.1: 
			displayed_ratio+= 1
			ratio-= 0.1
			timer_partial = timer()
			print(" >>", 10*displayed_ratio, "% (estimated remaining time: ", round((timer_partial-timer_start_input)/displayed_ratio*(10-displayed_ratio)), "s)")
	print()	
	
	
	# create TARGET group for SPANET (truth information about matching e.g. jet indicies of truth object)
	# here different outputs can be defined, if desired
	print(" > creating 'TARGET'")
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
			print("  >>", 10*displayed_ratio, "% (estimated remaining time: ", round((partial-start_group)/displayed_ratio*(10-displayed_ratio)), "s)")
	

	# create OTHER group for validation
	print()
	print(" > creating 'OTHER'")
	other_group = h5_file.create_group("OTHER")
	spanet_other_dimension = (number_of_events)


	# here branches can be customized if desired
	event_number = other_group.create_dataset("eventNumber", spanet_other_dimension, dtype=np.intc)
	mc_channel_number = other_group.create_dataset("mcChannelNumber", spanet_other_dimension, dtype=np.intc)
	classification_event_completion = other_group.create_dataset("classification_event_completion", spanet_other_dimension, dtype=np.intc)
	classification_true_higgs_decay = other_group.create_dataset("classification_true_higgs_decay", spanet_other_dimension, dtype=np.intc)
	classification_true_t1_decay = other_group.create_dataset("classification_true_t1_decay", spanet_other_dimension, dtype=np.intc)
	classification_true_t2_decay = other_group.create_dataset("classification_true_t2_decay", spanet_other_dimension, dtype=np.intc)
	classification_onshell_whad = other_group.create_dataset("classification_onshell_whad", spanet_other_dimension, dtype=np.intc)
	
	
	# prepare root_files for fast access
	root_event_number = root_file['matched/eventNumber'].array()
	root_mc_channel_number = root_file['matched/mcChannelNumber'].array()
	root_classification_event_completion = root_file['matched/classification_event_completion'].array()
	root_classification_true_higgs_decay = root_file['matched/classification_true_higgs_decay'].array()
	root_classification_true_t1_decay = root_file['matched/classification_true_t1_decay'].array()
	root_classification_true_t2_decay = root_file['matched/classification_true_t2_decay'].array()
	root_cclassification_onshell_whad = root_file['matched/classification_onshell_whad'].array()



	# timer for target and reset variables for printing
	ratio = 0.0
	displayed_ratio = 0
	start_group = timer()

	# loop for OTHER group
	for i in range(number_of_events):
		# fill datasets with custom data
		# important is fixed order within root file structured set in c++ code
		event_number[i] = root_event_number[i]
		mc_channel_number[i] = root_mc_channel_number[i]
		classification_event_completion[i] = root_classification_event_completion[i]
		classification_true_higgs_decay[i] = root_classification_true_higgs_decay[i]
		classification_true_t1_decay[i] = root_classification_true_t1_decay[i]
		classification_true_t2_decay[i] = root_classification_true_t2_decay[i]
		classification_onshell_whad[i] = root_cclassification_onshell_whad[i]

		# print progress
		ratio+= step
		if ratio>=0.1:
			displayed_ratio+= 1
			ratio-= 0.1
			partial = timer()
			print("  >>", 10*displayed_ratio, "% (estimated remaining time: ", round((partial-start_group)/displayed_ratio*(10-displayed_ratio)), "s)")
	print()



if __name__ == '__main__':
	main()
