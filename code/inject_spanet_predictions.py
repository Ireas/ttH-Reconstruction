import sys
import uproot # root in python
import h5py # h5 in python
import numpy as np


# ==========  INJECT SPANet PREDICTION  ============
# ==================================================
# injects SPANet prediction into root file as new prediction tree 
#
# arguments [H5 SPANET PREDICTION] [h5 ORIGINAL TRUTH MATCHED FILE]


ROOT_OUTPUT_DESTINATION = "/media/ireas/Data/v4/injected/injection.root"

def main():
	# validate arguments
	if(len(sys.argv)<1):
		print("Error: no .h5 SPANet prediction was given, exiting")
		exit()


	# open spanet prediction file
	spanet_prediction = h5py.File(sys.argv[1], 'r')

	# fill .h5 file with .root information
	with uproot.update(ROOT_OUTPUT_DESTINATION) as root_file:
		inject_prediction(root_file, spanet_prediction)


	# close .h5 file
	spanet_prediction.close()



def inject_prediction(root_file, spanet_prediction):
	
	root_file["spanet"] = {
		"t1_q1": np.array(spanet_prediction['TARGETS']['t1']['q1'][()], dtype=np.intc),
		"t1_q2": np.array(spanet_prediction['TARGETS']['t1']['q2'][()], dtype=np.intc),
		"t1_b" : np.array(spanet_prediction['TARGETS']['t1']['b'][()] , dtype=np.intc),
		"t1_assignment_probability" : np.array(spanet_prediction['TARGETS']['t1']['assignment_probability'][()] , dtype=np.float32),
		"t1_detection_probability" : np.array(spanet_prediction['TARGETS']['t1']['detection_probability'][()] , dtype=np.float32),
		"t1_marginal_probability" : np.array(spanet_prediction['TARGETS']['t1']['marginal_probability'][()] , dtype=np.float32),
		"t2_q1": np.array(spanet_prediction['TARGETS']['t2']['q1'][()], dtype=np.intc),
		"t2_q2": np.array(spanet_prediction['TARGETS']['t2']['q2'][()], dtype=np.intc),
		"t2_b" : np.array(spanet_prediction['TARGETS']['t2']['b'][()] , dtype=np.intc),
		"t2_assignment_probability" : np.array(spanet_prediction['TARGETS']['t2']['assignment_probability'][()] , dtype=np.float32),
		"t2_detection_probability" : np.array(spanet_prediction['TARGETS']['t2']['detection_probability'][()] , dtype=np.float32),
		"t2_marginal_probability" : np.array(spanet_prediction['TARGETS']['t2']['marginal_probability'][()] , dtype=np.float32),
		"HW_q1": np.array(spanet_prediction['TARGETS']['HW']['q1'][()], dtype=np.intc),
		"HW_q2": np.array(spanet_prediction['TARGETS']['HW']['q2'][()], dtype=np.intc),
		"HW_assignment_probability" : np.array(spanet_prediction['TARGETS']['HW']['assignment_probability'][()] , dtype=np.float32),
		"HW_detection_probability" : np.array(spanet_prediction['TARGETS']['HW']['detection_probability'][()] , dtype=np.float32),
		"HW_marginal_probability" : np.array(spanet_prediction['TARGETS']['HW']['marginal_probability'][()] , dtype=np.float32),
	}
	

if __name__ == '__main__':
	main()
