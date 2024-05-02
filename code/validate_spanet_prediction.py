import sys
import h5py
import numpy as np


# TODO: Implement partial top reconstructions?

CONFIDENCE_PROBABILITY = "marginal"
CONFIDENCE_THRESHOLD = 0.0



def main():
	if(len(sys.argv)<3):
		print("ERROR: Wrong Format")
		print("use \"python compare.py PREDICTION TRUTH\"")
		print("Exiting...")
		return

	# open files
	file_predi = h5py.File(sys.argv[1], 'r')
	file_truth = h5py.File(sys.argv[2], 'r')
	
	# transform to dictionaries for easy access
	dict_predi, dict_truth = files_to_dictionaries(file_predi, file_truth)

	number_of_events = len(dict_truth["HW_q1"])


	# compare resonance particles
	print("=========================")
	compare_HW(dict_predi, dict_truth, number_of_events)
	print("=========================")
	compare_t(dict_predi, dict_truth, number_of_events)
	print("=========================")
	compare_event(dict_predi, dict_truth, number_of_events)
	print("=========================")

	# close files
	file_predi.close()
	file_truth.close()


def compare_event(dict_predi, dict_truth, number_of_events):
	n_success = 0
	n_fail = 0
	n_unconfident = 0
	n_impossible = 0

	for event_idx in range(0, number_of_events):
		# access SPANet confidence
		t1_confidence = dict_predi["t1_p_"+CONFIDENCE_PROBABILITY][event_idx]	
		t2_confidence = dict_predi["t2_p_"+CONFIDENCE_PROBABILITY][event_idx]	
		HW_confidence = dict_predi["HW_p_"+CONFIDENCE_PROBABILITY][event_idx]	

		# access t1 partons
		t1_q1_p = dict_predi["t1_q1"][event_idx]
		t1_q2_p = dict_predi["t1_q2"][event_idx]
		t1_b_p  = dict_predi["t1_b"][event_idx]
		t1_q1_t = dict_truth["t1_q1"][event_idx]
		t1_q2_t = dict_truth["t1_q2"][event_idx]
		t1_b_t  = dict_truth["t1_b"][event_idx]
		
		# access t2 partons
		t2_q1_p = dict_predi["t2_q1"][event_idx]
		t2_q2_p = dict_predi["t2_q2"][event_idx]
		t2_b_p  = dict_predi["t2_b"][event_idx]
		t2_q1_t = dict_truth["t2_q1"][event_idx]
		t2_q2_t = dict_truth["t2_q2"][event_idx]
		t2_b_t  = dict_truth["t2_b"][event_idx]
		
		# access HW partons
		HW_q1_p = dict_predi["HW_q1"][event_idx]
		HW_q2_p = dict_predi["HW_q2"][event_idx]
		HW_q1_t = dict_truth["HW_q1"][event_idx]
		HW_q2_t = dict_truth["HW_q2"][event_idx]


		# check skip connection 	
		if t1_q1_t==-1 or t1_q2_t==-1 or t1_b_t==-1 or t2_q1_t==-1 or t2_q2_t==-1 or t2_b_t==-1 or HW_q1_t==-1 or HW_q2_t==-1:
			n_impossible+=1	
			continue
		if t1_confidence<CONFIDENCE_THRESHOLD or t2_confidence<CONFIDENCE_THRESHOLD or HW_confidence<CONFIDENCE_THRESHOLD:
			n_unconfident+= 1
			continue 



		# assigment and permutation check
		t_success = False
		HW_success = False
		## no perm
		## t1:(q1,q2) perm 
		## t2:(q1,q2) perm 
		## t1:(q1,q2) + t2:(q1,q2) perm 
		## (t1,t2) perm 
		## (t1,t2) + t1:(q1,q2) perm 
		## (t1,t2) + t2:(q1,q2) perm 
		## (t1,t2) + t1:(q1,q2) + t2:(q1,q2) perm 
		if 	(t1_q1_p==t1_q1_t and t1_q2_p==t1_q2_t and t2_q1_p==t2_q1_t and t2_q2_p==t2_q2_t and t1_b_p==t1_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t1_q2_t and t1_q2_p==t1_q1_t and t2_q1_p==t2_q1_t and t2_q2_p==t2_q2_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t1_q1_t and t1_q2_p==t1_q2_t and t2_q1_p==t2_q2_t and t2_q2_p==t2_q1_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t1_q2_t and t1_q2_p==t1_q1_t and t2_q1_p==t2_q2_t and t2_q2_p==t2_q1_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t2_q1_t and t1_q2_p==t2_q2_t and t2_q1_p==t1_q1_t and t2_q2_p==t1_q2_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q2_t and t1_q2_p==t2_q1_t and t2_q1_p==t1_q1_t and t2_q2_p==t1_q2_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q1_t and t1_q2_p==t2_q2_t and t2_q1_p==t1_q2_t and t2_q2_p==t1_q1_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q2_t and t1_q2_p==t2_q1_t and t2_q1_p==t1_q2_t and t2_q2_p==t1_q1_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t):
			t_success = True

		## no perm
		## HW:(q1,q2) perm 
		if (HW_q1_p==HW_q1_t and HW_q2_p==HW_q2_t) or (HW_q1_p==HW_q2_t and HW_q2_p==HW_q1_t):
			HW_success = True


		if t_success and HW_success:
			n_success+= 1
		else:
			n_fail+= 1
		

	# print evaluation
	print("> event")
	print("  -Total: ", number_of_events)
	print("  -Fails: ", n_fail)
	print("  -Successes: ", n_success)
	if n_success+n_fail>0:
		print("  -Success Ratio: ", np.round(n_success/(n_success+n_fail),4) )
	print("  -Unconfident", n_unconfident)
	print("  -Impossible", n_impossible)


def compare_t(dict_predi, dict_truth, number_of_events):
	n_success = 0
	n_fail = 0
	n_unconfident = 0
	n_impossible = 0
	
	for event_idx in range(0, number_of_events):
		# access SPANet confidence
		t1_confidence = dict_predi["t1_p_"+CONFIDENCE_PROBABILITY][event_idx]	
		t2_confidence = dict_predi["t2_p_"+CONFIDENCE_PROBABILITY][event_idx]	
		
		# access t1 partons
		t1_q1_p = dict_predi["t1_q1"][event_idx]
		t1_q2_p = dict_predi["t1_q2"][event_idx]
		t1_b_p  = dict_predi["t1_b"][event_idx]
		t1_q1_t = dict_truth["t1_q1"][event_idx]
		t1_q2_t = dict_truth["t1_q2"][event_idx]
		t1_b_t  = dict_truth["t1_b"][event_idx]
		
		# access t2 partons
		t2_q1_p = dict_predi["t2_q1"][event_idx]
		t2_q2_p = dict_predi["t2_q2"][event_idx]
		t2_b_p  = dict_predi["t2_b"][event_idx]
		t2_q1_t = dict_truth["t2_q1"][event_idx]
		t2_q2_t = dict_truth["t2_q2"][event_idx]
		t2_b_t  = dict_truth["t2_b"][event_idx]
		

		# check skip connection 
		if t1_q1_t==-1 or t1_q2_t==-1 or t1_b_t==-1 or t2_q1_t==-1 or t2_q2_t==-1 or t2_b_t==-1:
			n_impossible+=1	
			continue
		if t1_confidence<CONFIDENCE_THRESHOLD or t2_confidence<CONFIDENCE_THRESHOLD:
			n_unconfident+= 1
			continue 


		# assigment and permutation check
		## no perm
		## t1:(q1,q2) perm 
		## t2:(q1,q2) perm 
		## t1:(q1,q2) + t2:(q1,q2) perm 
		## (t1,t2) perm 
		## (t1,t2) + t1:(q1,q2) perm 
		## (t1,t2) + t2:(q1,q2) perm 
		## (t1,t2) + t1:(q1,q2) + t2:(q1,q2) perm 
		if 	(t1_q1_p==t1_q1_t and t1_q2_p==t1_q2_t and t2_q1_p==t2_q1_t and t2_q2_p==t2_q2_t and t1_b_p==t1_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t1_q2_t and t1_q2_p==t1_q1_t and t2_q1_p==t2_q1_t and t2_q2_p==t2_q2_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t1_q1_t and t1_q2_p==t1_q2_t and t2_q1_p==t2_q2_t and t2_q2_p==t2_q1_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t1_q2_t and t1_q2_p==t1_q1_t and t2_q1_p==t2_q2_t and t2_q2_p==t2_q1_t and t1_b_p==t1_b_t and t2_b_p==t2_b_t) or \
		 	(t1_q1_p==t2_q1_t and t1_q2_p==t2_q2_t and t2_q1_p==t1_q1_t and t2_q2_p==t1_q2_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q2_t and t1_q2_p==t2_q1_t and t2_q1_p==t1_q1_t and t2_q2_p==t1_q2_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q1_t and t1_q2_p==t2_q2_t and t2_q1_p==t1_q2_t and t2_q2_p==t1_q1_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t) or \
		 	(t1_q1_p==t2_q2_t and t1_q2_p==t2_q1_t and t2_q1_p==t1_q2_t and t2_q2_p==t1_q1_t and t1_b_p==t2_b_t and t2_b_p==t1_b_t):
			n_success+= 1
		else:
			n_fail+= 1

	# print evaluation
	print("> t1 & t2")
	print("  -Total: ", number_of_events)
	print("  -Fails: ", n_fail)
	print("  -Successes: ", n_success)
	if n_success+n_fail>0:
		print("  -Success Ratio: ", np.round(n_success/(n_success+n_fail),4) )
	print("  -Unconfident", n_unconfident)
	print("  -Impossible", n_impossible)



def compare_HW(dict_predi, dict_truth, number_of_events):
	n_success = 0
	n_fail = 0
	n_unconfident = 0
	n_impossible = 0
	
	for event_idx in range(0, number_of_events):
		# accessSPANet confidence
		HW_confidence = dict_predi["HW_p_"+CONFIDENCE_PROBABILITY][event_idx]	
		
		# access HW partons 
		HW_q1_p = dict_predi["HW_q1"][event_idx] 
		HW_q2_p = dict_predi["HW_q2"][event_idx]
		HW_q1_t = dict_truth["HW_q1"][event_idx]
		HW_q2_t = dict_truth["HW_q2"][event_idx]
	

		# check skip connection 
		if HW_q1_t==-1 or HW_q2_t==-1:
			n_impossible+=1	
			continue
		if HW_confidence<CONFIDENCE_THRESHOLD:
			n_unconfident+= 1
			continue 


		# assigment and permutation check
		## no perm
		## HW:(q1,q2) perm 
		if (HW_q1_p==HW_q1_t and HW_q2_p==HW_q2_t) or (HW_q1_p==HW_q2_t and HW_q2_p==HW_q1_t):
			n_success+= 1
		else:
			n_fail+= 1

	# print evaluation
	print("> HW")
	print("  -Total: ", number_of_events)
	print("  -Fails: ", n_fail)
	print("  -Successes: ", n_success)
	if n_success+n_fail>0:
		print("  -Success Ratio: ", np.round(n_success/(n_success+n_fail),4) )
	print("  -Unconfident", n_unconfident)
	print("  -Impossible", n_impossible)


def files_to_dictionaries(file_predi, file_truth):
	## prediction
	dict_predi = {}
	dict_predi["t1_b"]  = file_predi["TARGETS"]["t1"]["b"]
	dict_predi["t1_q1"] = file_predi["TARGETS"]["t1"]["q1"]
	dict_predi["t1_q2"] = file_predi["TARGETS"]["t1"]["q2"]
	dict_predi["t2_b"]  = file_predi["TARGETS"]["t2"]["b"]
	dict_predi["t2_q1"] = file_predi["TARGETS"]["t2"]["q1"]
	dict_predi["t2_q2"] = file_predi["TARGETS"]["t2"]["q2"]
	dict_predi["HW_q1"] = file_predi["TARGETS"]["HW"]["q1"]
	dict_predi["HW_q2"] = file_predi["TARGETS"]["HW"]["q2"]
	dict_predi["t1_p_assigment"] = file_predi["TARGETS"]["t1"]["assignment_probability"]
	dict_predi["t1_p_detection"] = file_predi["TARGETS"]["t1"]["detection_probability"]
	dict_predi["t1_p_marginal"]  = file_predi["TARGETS"]["t1"]["marginal_probability"]
	dict_predi["t2_p_assigment"] = file_predi["TARGETS"]["t2"]["assignment_probability"]
	dict_predi["t2_p_detection"] = file_predi["TARGETS"]["t2"]["detection_probability"]
	dict_predi["t2_p_marginal"]  = file_predi["TARGETS"]["t2"]["marginal_probability"]
	dict_predi["HW_p_assigment"] = file_predi["TARGETS"]["HW"]["assignment_probability"]
	dict_predi["HW_p_detection"] = file_predi["TARGETS"]["HW"]["detection_probability"]
	dict_predi["HW_p_marginal"]  = file_predi["TARGETS"]["HW"]["marginal_probability"]
	
	## true
	dict_truth = {}
	dict_truth["t1_b"]  = file_truth["TARGETS"]["t1"]["b"]
	dict_truth["t1_q1"] = file_truth["TARGETS"]["t1"]["q1"]
	dict_truth["t1_q2"] = file_truth["TARGETS"]["t1"]["q2"]
	dict_truth["t2_b"]  = file_truth["TARGETS"]["t2"]["b"]
	dict_truth["t2_q1"] = file_truth["TARGETS"]["t2"]["q1"]
	dict_truth["t2_q2"] = file_truth["TARGETS"]["t2"]["q2"]
	dict_truth["HW_q1"] = file_truth["TARGETS"]["HW"]["q1"]
	dict_truth["HW_q2"] = file_truth["TARGETS"]["HW"]["q2"]


	return dict_predi, dict_truth
	

if __name__=="__main__":
	main()
