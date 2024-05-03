import sys
import uproot # root in python
import h5py # h5 in python
import numpy as np


# ==========  INJECT SPANet PREDICTION  ============
# ==================================================
# injects SPANet prediction into root file as new prediction tree 
#
# arguments [ROOT FILE] [H5 SPANET PREDICTION]



def main():
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for injection was given, exiting")
		exit()
	if(len(sys.argv)<3):
		print("Error: no .h5 SPANet prediction was given, exiting")
		exit()


	# open .root file
	root_file = uproot.update(sys.argv[1])
	spanet_prediction = h5py.File(sys.argv[2], 'r')
	

	# fill .h5 file with .root information
	inject_prediction(root_file, spanet_prediction)

	
	# close .h5 file, .root cannot be closed
	spanet_prediction.close()
	root_file.close()



def inject_prediction(root_file, spanet_prediction):
	
	root_file["prediction"] = {
		"t1_q1": np.array(spanet_prediction['TARGETS']['t1']['q1'][()]),
		"t1_q2": np.array(spanet_prediction['TARGETS']['t1']['q2'][()]),
		"t1_b" : np.array(spanet_prediction['TARGETS']['t1']['b'][()]),
		"t2_q1": np.array(spanet_prediction['TARGETS']['t2']['q1'][()]),
		"t2_q2": np.array(spanet_prediction['TARGETS']['t2']['q2'][()]),
		"t2_b" : np.array(spanet_prediction['TARGETS']['t2']['b'][()]),
		"HW_q1": np.array(spanet_prediction['TARGETS']['HW']['q1'][()]),
		"HW_q2": np.array(spanet_prediction['TARGETS']['HW']['q2'][()])
	}



if __name__ == '__main__':
	main()
