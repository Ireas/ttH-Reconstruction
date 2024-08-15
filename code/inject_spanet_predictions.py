import sys
import shutil
import uproot # root in python
import h5py # h5 in python
import numpy as np



# ==========  INJECT SPANet PREDICTION  ============
# ==================================================
# injects SPANet prediction into root file as new prediction tree 


# CONSTANTS
MATCHED_ROOT_FILE = "/media/ireas/Data/v5/merged_root/all_8+j_1530766e_merged.root"
CONVERTED_H5_FILE = "/media/ireas/Data/v5/merged_h5/all_8+j_1530766e.h5"
SPANET_PREDICTION_H5_FILE = "/media/ireas/Data/v5/predicted/prediction_all_8+j_1530766e.h5"
INJECTED_OUTPUT_ROOT_FILE = "/media/ireas/Data/v5/injected/all_8+j_1530766e_injected.root"


def main():
	print(f"copying from {MATCHED_ROOT_FILE} to {INJECTED_OUTPUT_ROOT_FILE}")
	shutil.copy(MATCHED_ROOT_FILE, INJECTED_OUTPUT_ROOT_FILE)

	# open spanet prediction file
	print(f"creating lookup table {SPANET_PREDICTION_H5_FILE}")
	spanet_prediction_file = h5py.File(SPANET_PREDICTION_H5_FILE, 'r')
	converted_file = h5py.File(CONVERTED_H5_FILE, 'r')
	lookup_table = create_lookup_table(converted_file, spanet_prediction_file)

	# fill .h5 file with .root information
	print(f"filling lookup table into {INJECTED_OUTPUT_ROOT_FILE}")
	with uproot.update(INJECTED_OUTPUT_ROOT_FILE) as injected_root_file:
		with uproot.open(MATCHED_ROOT_FILE) as matched_root_file:
			inject_prediction(injected_root_file, matched_root_file, lookup_table)

	# close .h5 file
	converted_file.close()
	spanet_prediction_file.close()


def create_lookup_table(converted_file, spanet_prediction_h5_file):
	# define lookup dictionary for fast access, ordered by eventnumber, mcchannelnumber
	lookup_table = {}

	# access data
	event_numbers = np.array(converted_file["OTHER"]["eventNumber"][()], dtype=np.intc)
	mc_channel_numbers = np.array(converted_file["OTHER"]["mcChannelNumber"][()], dtype=np.intc)

	pred_t1_q1 = np.array(spanet_prediction_h5_file["TARGETS"]["t1"]["q1"][()], dtype=np.intc)
	pred_t1_q2 = np.array(spanet_prediction_h5_file["TARGETS"]["t1"]["q2"][()], dtype=np.intc)
	pred_t1_b = np.array(spanet_prediction_h5_file["TARGETS"]["t1"]["b"][()], dtype=np.intc)
	pred_t2_q1 = np.array(spanet_prediction_h5_file["TARGETS"]["t2"]["q1"][()], dtype=np.intc)
	pred_t2_q2 = np.array(spanet_prediction_h5_file["TARGETS"]["t2"]["q2"][()], dtype=np.intc)
	pred_t2_b = np.array(spanet_prediction_h5_file["TARGETS"]["t2"]["b"][()], dtype=np.intc)
	pred_HW_q1 = np.array(spanet_prediction_h5_file["TARGETS"]["HW"]["q1"][()], dtype=np.intc)
	pred_HW_q2 = np.array(spanet_prediction_h5_file["TARGETS"]["HW"]["q2"][()], dtype=np.intc)

	t1_assignment_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t1']['assignment_probability'][()] , dtype=np.float32)
	t1_detection_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t1']['detection_probability'][()] , dtype=np.float32)
	t1_marginal_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t1']['marginal_probability'][()] , dtype=np.float32)
	t2_assignment_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t2']['assignment_probability'][()] , dtype=np.float32)
	t2_detection_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t2']['detection_probability'][()] , dtype=np.float32)
	t2_marginal_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['t2']['marginal_probability'][()] , dtype=np.float32)
	HW_assignment_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['HW']['assignment_probability'][()] , dtype=np.float32)
	HW_detection_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['HW']['detection_probability'][()] , dtype=np.float32)
	HW_marginal_probabilities = np.array(spanet_prediction_h5_file['TARGETS']['HW']['marginal_probability'][()] , dtype=np.float32)

	# loop trough entries and add them accordendly
	for (
		mc_channel_number, event_number,
		t1_q1, t1_q2, t1_b, t2_q1, t2_q2, t2_b, HW_q1, HW_q2,
		t1_assignment_probability, t1_detection_probability, t1_marginal_probability,
		t2_assignment_probability, t2_detection_probability, t2_marginal_probability,
		HW_assignment_probability, HW_detection_probability, HW_marginal_probability,
		) in zip(
			mc_channel_numbers, event_numbers,
			pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
			t1_assignment_probabilities, t1_detection_probabilities, t1_marginal_probabilities,
			t2_assignment_probabilities, t2_detection_probabilities, t2_marginal_probabilities,
			HW_assignment_probabilities, HW_detection_probabilities, HW_marginal_probabilities,
		):
		
		lookup_table[(mc_channel_number, event_number)] = (
			t1_q1,
			t1_q2,
			t1_b,
			t2_q1,
			t2_q2,
			t2_b,
			HW_q1,
			HW_q2,
			t1_assignment_probability,
			t1_detection_probability,
			t1_marginal_probability,
			t2_assignment_probability,
			t2_detection_probability,
			t2_marginal_probability,
			HW_assignment_probability,
			HW_detection_probability,
			HW_marginal_probability
		)

	return lookup_table


def inject_prediction(injected_root_file, matched_root_file, lookup_table):
	# select correct entry by eventNumber and mcChannelNumber from lookup table!
	mc_channel_numbers = np.array( matched_root_file["matched/mcChannelNumber"].array(), dtype=np.intc )
	event_numbers = np.array( matched_root_file["matched/eventNumber"].array(), dtype=np.intc )
	
	# sanity check
	if (len(mc_channel_numbers)!=len(event_numbers)):
		print(f"Error: length are not equal! ({len(mc_channel_numbers)},{len(event_numbers)}), exiting")
		return
	
	# variable containers
	t1_q1 = np.array([], dtype=np.intc)
	t1_q2 = np.array([], dtype=np.intc)
	t1_b = np.array([], dtype=np.intc)
	t2_q1 = np.array([], dtype=np.intc)
	t2_q2 = np.array([], dtype=np.intc)
	t2_b = np.array([], dtype=np.intc)
	HW_q1 = np.array([], dtype=np.intc)
	HW_q2 = np.array([], dtype=np.intc)
	t1_detection_probabilities = np.array([], dtype=np.float32)
	t1_assignment_probabilities = np.array([], dtype=np.float32)
	t1_marginal_probabilities = np.array([], dtype=np.float32)
	t2_detection_probabilities = np.array([], dtype=np.float32)
	t2_assignment_probabilities = np.array([], dtype=np.float32)
	t2_marginal_probabilities = np.array([], dtype=np.float32)
	HW_detection_probabilities = np.array([], dtype=np.float32)
	HW_assignment_probabilities = np.array([], dtype=np.float32)
	HW_marginal_probabilities = np.array([], dtype=np.float32)

	# loop over all events and order them correctly
	for (mc_channel_number,event_number) in zip(mc_channel_numbers, event_numbers):
		# get entries
		entry = lookup_table[(mc_channel_number,event_number)]

		t1_q1 = np.append( t1_q1, entry[0] )
		t1_q2 = np.append( t1_q2, entry[1] )
		t1_b = np.append( t1_b, entry[2] )
		t2_q1 = np.append( t2_q1, entry[3] )
		t2_q2 = np.append( t2_q2, entry[4] )
		t2_b = np.append( t2_b, entry[5] )
		HW_q1 = np.append( HW_q1, entry[6] )
		HW_q2 = np.append( HW_q2, entry[7] )
		
		t1_detection_probabilities = np.append( t1_detection_probabilities, entry[8] )
		t1_assignment_probabilities = np.append( t1_assignment_probabilities, entry[9] )
		t1_marginal_probabilities = np.append( t1_marginal_probabilities, entry[10] )
		t2_detection_probabilities = np.append( t2_detection_probabilities, entry[11] )
		t2_assignment_probabilities = np.append( t2_assignment_probabilities, entry[12] )
		t2_marginal_probabilities = np.append( t2_marginal_probabilities, entry[13] )
		HW_detection_probabilities = np.append( HW_detection_probabilities, entry[14] )
		HW_assignment_probabilities = np.append( HW_assignment_probabilities, entry[15] )
		HW_marginal_probabilities = np.append( HW_marginal_probabilities, entry[16] )


	# inject them properly
	injected_root_file["spanet"] = {
		# mc sample ordering
		"mcChannelNumber": mc_channel_numbers,
		"eventNumber": event_numbers,
		# spanet jet assignment prediction
		"t1_q1": t1_q1,
		"t1_q2": t1_q2,
		"t1_b": t1_b,
		"t2_q1": t2_q1,
		"t2_q2": t2_q2,
		"t2_b": t2_b,
		"HW_q1": HW_q1,
		"HW_q2": HW_q2,
		# prediction confidences
		"t1_detection_probability": t1_detection_probabilities,
		"t1_assignment_probability": t1_assignment_probabilities,
		"t1_marginal_probability": t1_marginal_probabilities,
		"t2_detection_probability": t2_detection_probabilities,
		"t2_assignment_probability": t2_assignment_probabilities,
		"t2_marginal_probability": t2_marginal_probabilities,
		"HW_detection_probability": HW_detection_probabilities,
		"HW_assignment_probability": HW_assignment_probabilities,
		"HW_marginal_probability": HW_marginal_probabilities,
	}
	

if __name__ == '__main__':
	main()
