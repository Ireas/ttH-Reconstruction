import sys
import shutil
import uproot # root in python
import h5py # h5 in python
import numpy as np



# ==========  INJECT SPANet PREDICTION  ============
# ==================================================
# injects SPANet prediction into root file as new prediction tree 


# CONSTANTS

MATCHED_ROOT_FILE = "/media/ireas/Data/data/temp2.root"
CONVERTED_H5_FILE = "/media/ireas/Data/data/temp.h5"
SPANET_PREDICTION_H5_FILE = "/media/ireas/Data/data/temp_predicted.h5"
INJECTED_OUTPUT_ROOT_FILE = "/media/ireas/Data/data/temp_injected.root"



def main():
	print(f"copying from {MATCHED_ROOT_FILE} to {INJECTED_OUTPUT_ROOT_FILE}")
	shutil.copy(MATCHED_ROOT_FILE, INJECTED_OUTPUT_ROOT_FILE)

	# open spanet prediction file
	spanet_prediction_file = h5py.File(SPANET_PREDICTION_H5_FILE, 'r')
	converted_file = h5py.File(CONVERTED_H5_FILE, 'r')
	

	# fill .h5 file with .root information
	with uproot.update(INJECTED_OUTPUT_ROOT_FILE) as injected_root_file:
		with uproot.open(MATCHED_ROOT_FILE) as matched_root_file:
			inject_prediction(injected_root_file, matched_root_file, spanet_prediction_file, converted_file)

	# close .h5 file
	converted_file.close()
	spanet_prediction_file.close()


def inject_prediction(injected_root_file, matched_root_file, spanet_prediction_h5_file, converted_file):
	# select correct entry by eventNumber and mcChannelNumber from lookup table!
	event_numbers_root = np.array( matched_root_file["prepared/eventNumber"].array(), dtype=np.intc )
	event_numbers_pred = np.array(converted_file["OTHER"]["eventNumber"][()], dtype=np.intc)

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


	# sanity check
	for (i, (event_number_root, event_number_pred)) in enumerate(zip(event_numbers_root, event_numbers_pred)):
		if(event_number_root!=event_number_pred):
			print(f"  Warning: event numbers dont match at event {i}")
			return


	container_t1_q1 = []
	container_t1_q2 = []
	container_t1_b = []
	container_t2_q1 = []
	container_t2_q2 = []
	container_t2_b = []
	container_HW_q1 = []
	container_HW_q2 = []
	container_t1_assignment_probability = []
	container_t2_assignment_probability = []
	container_HW_assignment_probability = []
	container_t1_detection_probability = []
	container_t2_detection_probability = []
	container_HW_detection_probability = []
	container_t1_marginal_probability = []
	container_t2_marginal_probability = []
	container_HW_marginal_probability = []
	
	# loop trough entries and add them accordendly
	for (
		t1_q1, t1_q2, t1_b, t2_q1, t2_q2, t2_b, HW_q1, HW_q2,
		t1_assignment_probability, t1_detection_probability, t1_marginal_probability,
		t2_assignment_probability, t2_detection_probability, t2_marginal_probability,
		HW_assignment_probability, HW_detection_probability, HW_marginal_probability,
		) in zip(
			pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
			t1_assignment_probabilities, t1_detection_probabilities, t1_marginal_probabilities,
			t2_assignment_probabilities, t2_detection_probabilities, t2_marginal_probabilities,
			HW_assignment_probabilities, HW_detection_probabilities, HW_marginal_probabilities,
		):	
		
		container_t1_q1.append(t1_q1)
		container_t1_q2.append(t1_q2)
		container_t1_b.append(t1_b)
		container_t2_q1.append(t2_q1)
		container_t2_q2.append(t2_q2)
		container_t2_b.append(t2_b)
		container_HW_q1.append(HW_q1)
		container_HW_q2.append(HW_q2)
		container_t1_assignment_probability.append(t1_assignment_probability)
		container_t2_assignment_probability.append(t2_assignment_probability)
		container_HW_assignment_probability.append(HW_assignment_probability)
		container_t1_detection_probability.append(t1_detection_probability)
		container_t2_detection_probability.append(t2_detection_probability)
		container_HW_detection_probability.append(HW_detection_probability)
		container_t1_marginal_probability.append(t1_marginal_probability)
		container_t2_marginal_probability.append(t2_marginal_probability)
		container_HW_marginal_probability.append(HW_marginal_probability)

	container_t1_q1 = np.array( container_t1_q1, dtype=np.intc)
	container_t1_q2 = np.array( container_t1_q2, dtype=np.intc)
	container_t1_b = np.array( container_t1_b, dtype=np.intc)
	container_t2_q1 = np.array( container_t2_q1, dtype=np.intc)
	container_t2_q2 = np.array( container_t2_q2, dtype=np.intc)
	container_t2_b = np.array( container_t2_b, dtype=np.intc)
	container_HW_q1 = np.array( container_HW_q1, dtype=np.intc)
	container_HW_q2 = np.array( container_HW_q2, dtype=np.intc)
	container_t1_assignment_probability = np.array( container_t1_assignment_probability, dtype=np.float32)
	container_t2_assignment_probability = np.array( container_t2_assignment_probability, dtype=np.float32)
	container_HW_assignment_probability = np.array( container_HW_assignment_probability, dtype=np.float32)
	container_t1_detection_probability = np.array( container_t1_detection_probability, dtype=np.float32)
	container_t2_detection_probability = np.array( container_t2_detection_probability, dtype=np.float32)
	container_HW_detection_probability = np.array( container_HW_detection_probability, dtype=np.float32)
	container_t1_marginal_probability = np.array( container_t1_marginal_probability, dtype=np.float32)
	container_t2_marginal_probability = np.array( container_t2_marginal_probability, dtype=np.float32)
	container_HW_marginal_probability = np.array( container_HW_marginal_probability, dtype=np.float32)
	

	# inject them properly
	injected_root_file["spanet"] = {
		# mc sample ordering
		"eventNumber": event_numbers_pred,
		# spanet jet assignment prediction
		"t1_q1": container_t1_q1,
		"t1_q2": container_t1_q2,
		"t1_b": container_t1_b,
		"t2_q1": container_t2_q1,
		"t2_q2": container_t2_q2,
		"t2_b": container_t2_b,
		"HW_q1": container_HW_q1,
		"HW_q2": container_HW_q2,
		# prediction confidences
		"t1_assignment_probability": 	container_t1_assignment_probability,
		"t2_assignment_probability": 	container_t2_assignment_probability,
		"HW_assignment_probability": 	container_HW_assignment_probability,
		"t1_detection_probability": 	container_t1_detection_probability,
		"t2_detection_probability": 	container_t2_detection_probability,
		"HW_detection_probability": 	container_HW_detection_probability,
		"t1_marginal_probability": 	container_t1_marginal_probability,
		"t2_marginal_probability": 	container_t2_marginal_probability,
		"HW_marginal_probability": 	container_HW_marginal_probability,
	}
	

if __name__ == '__main__':
	main()
