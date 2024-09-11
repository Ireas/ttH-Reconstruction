import os
import sys
import h5py # h5 in python
import numpy as np
from timeit import default_timer as timer


# constants
MAX_EVENTS_PER_FILE = -1 #-1 to use all events
SUCCESSFUL_ONLY = False

INPUT_FOLDER = "all_5+j_truncated_50_ttHWW_amplified/"
INPUT_PATH = "/media/ireas/Data/v6/converted/"
OUTPUT_PATH = "/media/ireas/Data/v6/merged_h5/"


successful_event_history = {}

# main method 
def main():
	# access input
	input_files = []

	input_directory = INPUT_PATH+INPUT_FOLDER
	print(f"Accessing input files in {input_directory}")
	for input_destination in os.listdir(input_directory):
		print(f" > {input_destination}")
		input_file = h5py.File(input_directory+input_destination, 'r')
		input_files.append(input_file)
	print("")
	print(f"{len(input_files)} files found!")


	# create output
	print("Create Output")
	output_destination = OUTPUT_PATH + INPUT_FOLDER[:-1] + ".h5"
	print(f" > {output_destination}")
	output_file = h5py.File(output_destination, 'w')
	print("")


	# fill with data
	number_of_events = fill_output(output_file, input_files)	
	

	# close all input files
	for input_file in input_files:
		input_file.close()

	if(SUCCESSFUL_ONLY):
		os.rename(output_destination, output_destination[:-3]+"_"+str(number_of_events)+"e+_SuccessfulOnly.h5")
	else:
		os.rename(output_destination, output_destination[:-3]+"_"+str(number_of_events)+"e.h5")

def get_success_history(input_file):
	# access input file 
	events_in_input_file = input_file['INPUTS']['Source']['MASK'][()].shape[0] 
	history = np.array([])

	ind_t1q1 = input_file['TARGETS']['t1']['q1'][()]
	ind_t1q2 = input_file['TARGETS']['t1']['q2'][()]
	ind_t1b = input_file['TARGETS']['t1']['b'][()]
	ind_t2q1 = input_file['TARGETS']['t2']['q1'][()]
	ind_t2q2 = input_file['TARGETS']['t2']['q2'][()]
	ind_t2b = input_file['TARGETS']['t2']['b'][()]
	ind_HWq1 = input_file['TARGETS']['HW']['q1'][()]
	ind_HWq2 = input_file['TARGETS']['HW']['q2'][()]

	for i in range(events_in_input_file):
		# empty event
		if ( ind_t1q1[i]==0 and ind_t1q2[i]==0 and ind_t1b[i]==0 and ind_t2q1[i]==0 and ind_t2q2[i]==0 and ind_t2b[i]==0 and ind_HWq1[i]==0 and ind_HWq2[i]==0): 
			history = np.append(history, [False])
			continue
	
		# and invalid assignments
		if ( ind_t1q1[i]==-1 or ind_t1q2[i]==-1 or ind_t1b[i]==-1 or ind_t2q1[i]==-1 or ind_t2q2[i]==-1 or ind_t2b[i]==-1 or ind_HWq1[i]==-1 or ind_HWq2[i]==-1): 
			history = np.append(history, [False])
			continue

		history = np.append(history, [True])


	return history


def fill_output(output_file, input_files):
	# determine length of files first
	events_in_output_file = 0
	max_njets_in_output_file  = 0

	print("Count Events")
	for input_file in input_files:	
		# get number of events in file	
		events_in_input_file = input_file['INPUTS']['Source']['MASK'][()].shape[0] 
		jets_in_input_file = input_file['INPUTS']['Source']['MASK'][()].shape[1] 
		
		# trim if needed
		number_of_events = events_in_input_file
		if(SUCCESSFUL_ONLY):
			successful_event_history[str(input_file)] = get_success_history(input_file)
			number_of_events = int(successful_event_history[str(input_file)].sum())


		if MAX_EVENTS_PER_FILE>0 and number_of_events>MAX_EVENTS_PER_FILE:
			number_of_events = MAX_EVENTS_PER_FILE
		
		max_njets_in_output_file = jets_in_input_file if jets_in_input_file>max_njets_in_output_file else max_njets_in_output_file
		
		# print
		print(f" > {number_of_events} events \t| {os.path.splitext(os.path.basename(input_file.filename))[0]}")
		
		# add to overall length of output
		events_in_output_file+= number_of_events

	
	# prepare timer
	t_start = timer()

	# print skipped events
	print()		
	print(f"Including {events_in_output_file} events")
	
	# Create INPUTS group
	output_inputs_dimension = (events_in_output_file, max_njets_in_output_file)
	output_met_dimension = (events_in_output_file)

	print()		
	print("Fill group INPUTS")
	input_group = output_file.create_group("INPUTS")
	source_group = input_group.create_group("Source")
	met_group = input_group.create_group("Met")
	
	mask = source_group.create_dataset("MASK", output_inputs_dimension, dtype=bool)	
	jet_e = source_group.create_dataset("energy", output_inputs_dimension, dtype=np.float32)	
	jet_pt = source_group.create_dataset("pt", output_inputs_dimension, dtype=np.float32)	
	jet_eta = source_group.create_dataset("eta", output_inputs_dimension, dtype=np.float32)	
	jet_phi = source_group.create_dataset("phi", output_inputs_dimension, dtype=np.float32)	
	jet_btag = source_group.create_dataset("btag", output_inputs_dimension, dtype=bool)	
	
	met_value = met_group.create_dataset("value", output_met_dimension, dtype=np.float32)	
	met_phi = met_group.create_dataset("phi", output_met_dimension, dtype=np.float32)	

	
	# move start index when using multiple files	
	start_index = 0

	for input_file in input_files:
		number_of_events = input_file['INPUTS']['Source']['MASK'][()].shape[0]
		event_limit = input_file['INPUTS']['Source']['MASK'][()].shape[0]

		if(SUCCESSFUL_ONLY):
			event_limit = int(successful_event_history[str(input_file)].sum())
		
		event_limit = MAX_EVENTS_PER_FILE if (MAX_EVENTS_PER_FILE>0 and event_limit>MAX_EVENTS_PER_FILE) else event_limit
		jet_limit = input_file['INPUTS']['Source']['MASK'][()].shape[1] 

		t_new = timer()
		if start_index==0:
			print(f" > ETA ---\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")
		else:
			eta = round( (events_in_output_file-start_index) * (t_new-t_start)/start_index )
			print(f" > ETA: {eta}s\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")

		# access
		input_mask = input_file['INPUTS']['Source']['MASK'][()]
		input_energy = input_file['INPUTS']['Source']['energy'][()]
		input_pt = input_file['INPUTS']['Source']['pt'][()]
		input_eta = input_file['INPUTS']['Source']['eta'][()]
		input_phi = input_file['INPUTS']['Source']['phi'][()]
		input_btag = input_file['INPUTS']['Source']['btag'][()]
		input_met_value = input_file['INPUTS']['Met']['value'][()]
		input_met_phi = input_file['INPUTS']['Met']['phi'][()]

		if(SUCCESSFUL_ONLY):
			success_history = successful_event_history[str(input_file)]
		success_events = 0


		# fill
		for i in range(number_of_events):
			
			for j in range(jet_limit):
				mask[start_index+success_events,j] = input_mask[i,j] 
				jet_e[start_index+success_events,j] = input_energy[i,j] 
				jet_pt[start_index+success_events,j] = input_pt[i,j] 
				jet_eta[start_index+success_events,j] = input_eta[i,j] 
				jet_phi[start_index+success_events,j] = input_phi[i,j] 
				jet_btag[start_index+success_events,j] = input_btag[i,j]
			
			# global variabls
			met_value[start_index+success_events] = input_met_value[i]
			met_phi[start_index+success_events] = input_met_phi[i]

			# count success event
			if(SUCCESSFUL_ONLY):
				if(success_history[i]):
					success_events+= 1	
				else:
					continue
			else:
				success_events+= 1
			
			# manual break
			if(success_events>=event_limit):
				break
				
		start_index+= success_events



	# reset timer
	t_start = timer()
	

	# Create TARGETS group
	output_targets_dimension = (events_in_output_file)
	
	print()		
	print("Fill group TARGETS")
	target_group = output_file.create_group("TARGETS")
	t1_group = target_group.create_group("t1")
	b1 = t1_group.create_dataset("b", output_targets_dimension, dtype=int)
	q1_1 = t1_group.create_dataset("q1", output_targets_dimension, dtype=int)
	q1_2 = t1_group.create_dataset("q2", output_targets_dimension, dtype=int)
	
	t2_group = target_group.create_group("t2")
	b2 = t2_group.create_dataset("b", output_targets_dimension, dtype=int)
	q2_1 = t2_group.create_dataset("q1", output_targets_dimension, dtype=int)
	q2_2 = t2_group.create_dataset("q2", output_targets_dimension, dtype=int)
	
	HW_group = target_group.create_group("HW")
	HW_1 = HW_group.create_dataset("q1", output_targets_dimension, dtype=int)
	HW_2 = HW_group.create_dataset("q2", output_targets_dimension, dtype=int)
	
	# move start index when using multiple files	
	start_index = 0
		
	for input_file in input_files:
		number_of_events = input_file['INPUTS']['Source']['MASK'][()].shape[0]
		event_limit = input_file['INPUTS']['Source']['MASK'][()].shape[0]
		
		if(SUCCESSFUL_ONLY):
			event_limit = int(successful_event_history[str(input_file)].sum())
		
		event_limit = MAX_EVENTS_PER_FILE if (MAX_EVENTS_PER_FILE>0 and event_limit>MAX_EVENTS_PER_FILE) else event_limit
		
		t_new = timer()
		if start_index==0:
			print(f" > ETA ---\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")
		else:
			eta = round( (events_in_output_file-start_index) * (t_new-t_start)/start_index )
			print(f" > ETA: {eta}s\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")

		# access
		input_t1_b = input_file['TARGETS']['t1']['b'][()]
		input_t1_q1 = input_file['TARGETS']['t1']['q1'][()]
		input_t1_q2 = input_file['TARGETS']['t1']['q2'][()]
		input_t2_b = input_file['TARGETS']['t2']['b'][()]
		input_t2_q1 = input_file['TARGETS']['t2']['q1'][()]
		input_t2_q2 = input_file['TARGETS']['t2']['q2'][()]
		input_HW_q1 = input_file['TARGETS']['HW']['q1'][()]
		input_HW_q2 = input_file['TARGETS']['HW']['q2'][()]
		
		if(SUCCESSFUL_ONLY):
			success_history = successful_event_history[str(input_file)]
		success_events = 0

		# fill
		for i in range(number_of_events):
			b1[start_index+success_events] = input_t1_b[i] 
			q1_1[start_index+success_events] = input_t1_q1[i] 
			q1_2[start_index+success_events] = input_t1_q2[i] 
			b2[start_index+success_events] = input_t2_b[i] 
			q2_1[start_index+success_events] = input_t2_q1[i] 
			q2_2[start_index+success_events] = input_t2_q2[i] 
			HW_1[start_index+success_events] = input_HW_q1[i] 
			HW_2[start_index+success_events] = input_HW_q2[i] 
			

			
			if(SUCCESSFUL_ONLY):
				if(success_history[i]):
					success_events+= 1	
				else:
					continue
			else:
				success_events+= 1

			if(success_events>=event_limit):
				break
			
		start_index+= success_events
		
	# reset timer
	t_start = timer()
	

	# Create OTHER group
	print()		
	print("Fill group OTHER")
	output_other_dimension = (events_in_output_file)

	other_group = output_file.create_group("OTHER")
	event_number = other_group.create_dataset("eventNumber", output_other_dimension, dtype=np.intc)
	mc_channel_number = other_group.create_dataset("mcChannelNumber", output_other_dimension, dtype=np.intc)
	classification_event_channel = other_group.create_dataset("classification_event_channel", output_other_dimension, dtype=np.intc)
	classification_true_higgs_decay = other_group.create_dataset("classification_true_higgs_decay", output_other_dimension, dtype=np.intc)
	classification_true_t1_decay = other_group.create_dataset("classification_true_t1_decay", output_other_dimension, dtype=np.intc)
	classification_true_t2_decay = other_group.create_dataset("classification_true_t2_decay", output_other_dimension, dtype=np.intc)
	classification_onshell_whad = other_group.create_dataset("classification_onshell_whad", output_other_dimension, dtype=np.intc)
	
	# move start index when using multiple files	
	start_index = 0
	
	for input_file in input_files:
		number_of_events = input_file['INPUTS']['Source']['MASK'][()].shape[0]
		event_limit = input_file['INPUTS']['Source']['MASK'][()].shape[0]

		if(SUCCESSFUL_ONLY):
			event_limit = int(successful_event_history[str(input_file)].sum())
		
		event_limit = MAX_EVENTS_PER_FILE if (MAX_EVENTS_PER_FILE>0 and event_limit>MAX_EVENTS_PER_FILE) else event_limit
		
		t_new = timer()
		if start_index==0:
			print(f" > ETA ---\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")
		else:
			eta = round( (events_in_output_file-start_index) * (t_new-t_start)/start_index )
			print(f" > ETA: {eta}s\t| currently processing {event_limit} events from {os.path.splitext(os.path.basename(input_file.filename))[0]}")
		
		# access
		input_event_number = input_file['OTHER']['eventNumber'][()]
		input_mc_channel_number = input_file['OTHER']['mcChannelNumber'][()]
		input_classification_event_channel = input_file['OTHER']['classification_event_channel'][()]
		input_classification_true_higgs_decay = input_file['OTHER']['classification_true_higgs_decay'][()]
		input_classification_true_t1_decay = input_file['OTHER']['classification_true_t1_decay'][()]
		input_classification_true_t2_decay = input_file['OTHER']['classification_true_t2_decay'][()]
		input_classification_onshell_whad = input_file['OTHER']['classification_onshell_whad'][()]
				
		if(SUCCESSFUL_ONLY):
			success_history = successful_event_history[str(input_file)]
		success_events = 0

		# fill
		for i in range(number_of_events):
			event_number[start_index+success_events] = input_event_number[i] 
			mc_channel_number[start_index+success_events] = input_mc_channel_number[i] 
			classification_event_channel[start_index+success_events] = input_classification_event_channel[i] 
			classification_true_higgs_decay[start_index+success_events] = input_classification_true_higgs_decay[i] 
			classification_true_t1_decay[start_index+success_events] = input_classification_true_t1_decay[i] 
			classification_true_t2_decay[start_index+success_events] = input_classification_true_t2_decay[i] 
			classification_onshell_whad[start_index+success_events] = input_classification_onshell_whad[i] 
						
			if(SUCCESSFUL_ONLY):
				if(success_history[i]):
					success_events+= 1	
				else:
					continue
			else:
				success_events+= 1

			if(success_events>=event_limit):
				break

		start_index+= success_events


	return events_in_output_file



if __name__ == '__main__':
	main()
