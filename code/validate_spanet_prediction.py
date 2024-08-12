import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt


# constants
CONFIDENCE_PROBABILITY = "assignment"
CONFIDENCE_THRESHOLD = 0.0

FULL_VALID_EVENTS_ONLY = False
SHOW_PLOTS = False

PREDICTION_FILE = "/media/ireas/Data/v5/predicted/prediction_all_8+j_1530766e.h5"
TRUTH_FILE = "/media/ireas/Data/v5/merged_h5/all_8+j_1530766e.h5"


def main():
	# open files
    file_pred = h5py.File(PREDICTION_FILE, 'r')
    file_true = h5py.File(TRUTH_FILE, 'r')
	
	# transform to dictionaries for easy access
    dict_pred, dict_true = files_to_dictionaries(file_pred, file_true)
    
    validate_all_events(dict_pred, dict_true, FULL_VALID_EVENTS_ONLY)
    
    # plot as baptiste 
    #plot_t1(dict_pred, dict_true, "assignment", FULL_VALID_EVENTS_ONLY, True)
    #plot_t1(dict_pred, dict_true, "assignment", FULL_VALID_EVENTS_ONLY, False)
    #plot_t1(dict_pred, dict_true, "detection", FULL_VALID_EVENTS_ONLY, True)
    #plot_t1(dict_pred, dict_true, "detection", FULL_VALID_EVENTS_ONLY, False)
    #plot_t1(dict_pred, dict_true, "marginal", FULL_VALID_EVENTS_ONLY, True)
    #plot_t1(dict_pred, dict_true, "marginal", FULL_VALID_EVENTS_ONLY, False)


    # close files
    file_pred.close()
    file_true.close()

    # get information
	#number_of_events = len(dict_true["t1_q1"])


def plot_t1(dict_pred, dict_true, variable_name, use_full_valid_events_only, normalize): 
    # transform    
    confidence_correct = np.array([])
    confidence_swapped = np.array([])
    confidence_failed = np.array([])
    confidence_impossible = np.array([])
   

    for (
            pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
            true_t1_q1, true_t1_q2, true_t1_b, true_t2_q1, true_t2_q2, true_t2_b, true_HW_q1, true_HW_q2,
            pred_confidence_t1, pred_confidence_t2, pred_confidence_HW
        ) in zip(
            dict_pred["t1_q1"], dict_pred["t1_q2"], dict_pred["t1_b"], dict_pred["t2_q1"], dict_pred["t2_q2"], dict_pred["t2_b"], dict_pred["HW_q1"], dict_pred["HW_q2"],
            dict_true["t1_q1"], dict_true["t1_q2"], dict_true["t1_b"], dict_true["t2_q1"], dict_true["t2_q2"], dict_true["t2_b"], dict_true["HW_q1"], dict_true["HW_q2"],
            dict_pred["t1_prediction_"+variable_name], dict_pred["t2_prediction_"+variable_name], dict_pred["HW_prediction_"+variable_name]
        ):
            # confidence threshold
            if pred_confidence_t1<CONFIDENCE_THRESHOLD:
                continue
            if pred_confidence_t2<CONFIDENCE_THRESHOLD:
                continue
            if pred_confidence_HW<CONFIDENCE_THRESHOLD:
                continue

            # check which resonance are possible
            possible_t1W = (true_t1_q1!=-1 and true_t1_q2!=-1)
            possible_t1 = (possible_t1W and true_t1_b!=-1)
            possible_t2W = (true_t2_q1!=-1 and true_t2_q2!=-1)
            possible_t2 = (possible_t2W and true_t2_b!=-1)
            possible_ttbar = (possible_t1 and possible_t2)
            possible_HW = (true_HW_q1!=-1 and true_HW_q2!=-1)
            possible_event = (possible_ttbar and possible_HW)

            # skip events that are not fully assigned
            if use_full_valid_events_only and not possible_event:
                continue
            
            if possible_t1:
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b):
                    confidence_correct = np.append(confidence_correct, [pred_confidence_t1])
                # permutate q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b):
                    confidence_swapped = np.append(confidence_swapped, [pred_confidence_t1])
                else:
                    confidence_failed = np.append(confidence_failed, [pred_confidence_t1])
            else: 
                confidence_impossible = np.append(confidence_impossible, [pred_confidence_t1])

 
    bins = np.linspace(0,1,25)
    if normalize:
        # plot histogram
        plt.figure()

	    # weights
        weights_correct = np.ones_like(confidence_correct)/len(confidence_correct)
        weights_swapped = np.ones_like(confidence_swapped)/len(confidence_swapped)
        weights_failed = np.ones_like(confidence_failed)/len(confidence_failed)
        weights_impossible = np.ones_like(confidence_impossible)/len(confidence_impossible)

	    # histograms
        plt.hist(confidence_correct, bins, label="correct", fill=False, histtype="step", weights=weights_correct)
        plt.hist(confidence_swapped, bins, label="swapped", fill=False, histtype="step", weights=weights_swapped)
        plt.hist(confidence_failed, bins, label="failed", fill=False, histtype="step", weights=weights_failed)
        plt.hist(confidence_impossible, bins, label="impossible", fill=False, histtype="step", weights=weights_impossible)

	    # labels
        plt.xlabel(f"SPANet Probability ({variable_name})")
        plt.ylabel("Events Normalized")
        plt.xlim([0,1])
        plt.legend()

        # output
        plt.savefig(f"spanet_probability_{variable_name}_normed.png")
    else:
        # plot histogram
        plt.figure()

	    # histograms
        plt.hist(confidence_correct, bins, label="correct", fill=False, histtype="step")
        plt.hist(confidence_swapped, bins, label="swapped", fill=False, histtype="step")
        plt.hist(confidence_failed, bins, label="failed", fill=False, histtype="step")
        plt.hist(confidence_impossible, bins, label="impossible", fill=False, histtype="step")

	    # labels
        plt.xlabel(f"SPANet Probability ({variable_name})")
        plt.ylabel("Events")
        plt.xlim([0,1])
        plt.yscale("log")
        plt.legend()

        # output
        plt.savefig(f"spanet_probability_{variable_name}.png")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.clf()



def validate_all_events(dict_pred, dict_true, use_full_valid_events_only): 
    # define counters for event counting
    n_t1W_correct = 0
    n_t1W_false = 0
    n_t1_correct = 0
    n_t1_false = 0
    n_t2W_correct = 0
    n_t2W_false = 0
    n_t2_correct = 0
    n_t2_false = 0
    n_HW_correct = 0
    n_HW_false = 0
    n_ttbar_correct = 0
    n_ttbar_false = 0
    n_event_correct = 0
    n_event_false = 0

    for (
            pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
            true_t1_q1, true_t1_q2, true_t1_b, true_t2_q1, true_t2_q2, true_t2_b, true_HW_q1, true_HW_q2,
            pred_confidence_t1, pred_confidence_t2, pred_confidence_HW
        ) in zip(
            dict_pred["t1_q1"], dict_pred["t1_q2"], dict_pred["t1_b"], dict_pred["t2_q1"], dict_pred["t2_q2"], dict_pred["t2_b"], dict_pred["HW_q1"], dict_pred["HW_q2"],
            dict_true["t1_q1"], dict_true["t1_q2"], dict_true["t1_b"], dict_true["t2_q1"], dict_true["t2_q2"], dict_true["t2_b"], dict_true["HW_q1"], dict_true["HW_q2"],
            dict_pred["t1_prediction_"+CONFIDENCE_PROBABILITY], dict_pred["t2_prediction_"+CONFIDENCE_PROBABILITY], dict_pred["HW_prediction_"+CONFIDENCE_PROBABILITY]
        ):
            # confidence threshold
            if pred_confidence_t1<CONFIDENCE_THRESHOLD:
                continue
            if pred_confidence_t2<CONFIDENCE_THRESHOLD:
                continue
            if pred_confidence_HW<CONFIDENCE_THRESHOLD:
                continue

            # check which resonance are possible
            possible_t1W = (true_t1_q1!=-1 and true_t1_q2!=-1)
            possible_t1 = (possible_t1W and true_t1_b!=-1)
            possible_t2W = (true_t2_q1!=-1 and true_t2_q2!=-1)
            possible_t2 = (possible_t2W and true_t2_b!=-1)
            possible_ttbar = (possible_t1 and possible_t2)
            possible_HW = (true_HW_q1!=-1 and true_HW_q2!=-1)
            possible_event = (possible_ttbar and possible_HW)

            # skip events that are not fully assigned
            if use_full_valid_events_only and not possible_event:
                continue

            # check only possible particle
            if possible_t1W:
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2):
                    n_t1W_correct+= 1
                # permutate q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1):
                    n_t1W_correct+= 1
                else:
                    n_t1W_false+= 1
            

            if possible_t1:
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b):
                    n_t1_correct+= 1
                # permutate q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b):
                    n_t1_correct+= 1
                else:
                    n_t1_false+= 1
            
            
            if possible_t2W:
                if(pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2):
                    n_t2W_correct+= 1
                # permutate q1,q2
                elif(pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1):
                    n_t2W_correct+= 1
                else:
                    n_t2W_false+= 1
            
            
            if possible_t2:
                if(pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    n_t2_correct+= 1
                # permutate q1,q2
                elif(pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    n_t2_correct+= 1
                else:
                    n_t2_false+= 1
           

            if possible_HW:
                if(pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    n_HW_correct+= 1
                # permutate q1,q2
                elif(pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    n_HW_correct+= 1
                else:
                    n_HW_false+= 1
            

            if possible_ttbar:
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    n_ttbar_correct+= 1
                # permutate t1:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    n_ttbar_correct+= 1
                # permutate t2:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    n_ttbar_correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    n_ttbar_correct+= 1
                else:
                    n_ttbar_false+= 1
           

            if possible_event:
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    n_event_correct+= 1
                # permutate t1:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    n_event_correct+= 1
                # permutate t2:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    n_event_correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    n_event_correct+= 1
                # permutate HW:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    n_event_correct+= 1
                # permutate t1:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    n_event_correct+= 1
                # permutate t2:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    n_event_correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    n_event_correct+= 1
                else:
                    n_event_false+= 1

    print()
    print("==========  OUTPUT  ==========")
    print("Particle: total possible |  correct,false  |  ratio_correct")
    print()
    print(f"t1W:  {n_t1W_correct+n_t1W_false}  |  {n_t1W_correct}, {n_t1W_false}  ->  {np.round(n_t1W_correct/float(n_t1W_correct+n_t1W_false), 4)}")
    print(f"t1:  {n_t1_correct+n_t1_false}  |  {n_t2_correct}, {n_t1_false}  ->   {np.round(n_t1_correct/float(n_t1_correct+n_t1_false), 4)}")
    print(f"t2W:  {n_t2W_correct+n_t2W_false}  |  {n_t2W_correct}, {n_t2W_false}  ->  {np.round(n_t2W_correct/float(n_t2W_correct+n_t2W_false), 4)}")
    print(f"t2:  {n_t2_correct+n_t2_false}  |  {n_t2_correct}, {n_t2_false}  ->   {np.round(n_t2_correct/float(n_t2_correct+n_t2_false), 4)}")
    print(f"ttbar:  {n_ttbar_correct+n_ttbar_false}  |  {n_ttbar_correct}, {n_ttbar_false}  ->  {np.round(n_ttbar_correct/float(n_ttbar_correct+n_ttbar_false), 4)}")
    print(f"HW:  {n_HW_correct+n_HW_false}  |  {n_HW_correct}, {n_HW_false}  ->  {np.round(n_HW_correct/float(n_HW_correct+n_HW_false), 4)}")
    print()
    print(f"Event:  {n_event_correct+n_event_false}  |  {n_event_correct}, {n_event_false}  ->  {np.round(n_event_correct/float(n_event_correct+n_event_false), 4)}")
    print()


def files_to_dictionaries(file_pred, file_true):
    ## prediction
    dict_pred = {}
    dict_pred["t1_b"]  = file_pred["TARGETS"]["t1"]["b"]
    dict_pred["t1_q1"] = file_pred["TARGETS"]["t1"]["q1"]
    dict_pred["t1_q2"] = file_pred["TARGETS"]["t1"]["q2"]
    dict_pred["t2_b"]  = file_pred["TARGETS"]["t2"]["b"]
    dict_pred["t2_q1"] = file_pred["TARGETS"]["t2"]["q1"]
    dict_pred["t2_q2"] = file_pred["TARGETS"]["t2"]["q2"]
    dict_pred["HW_q1"] = file_pred["TARGETS"]["HW"]["q1"]
    dict_pred["HW_q2"] = file_pred["TARGETS"]["HW"]["q2"]
    
    ## resonance prediction probabilites
    dict_pred["t1_prediction_assignment"] = file_pred["TARGETS"]["t1"]["assignment_probability"]
    dict_pred["t1_prediction_detection"] = file_pred["TARGETS"]["t1"]["detection_probability"]
    dict_pred["t1_prediction_marginal"]  = file_pred["TARGETS"]["t1"]["marginal_probability"]
    dict_pred["t2_prediction_assignment"] = file_pred["TARGETS"]["t2"]["assignment_probability"]
    dict_pred["t2_prediction_detection"] = file_pred["TARGETS"]["t2"]["detection_probability"]
    dict_pred["t2_prediction_marginal"]  = file_pred["TARGETS"]["t2"]["marginal_probability"]
    dict_pred["HW_prediction_assignment"] = file_pred["TARGETS"]["HW"]["assignment_probability"]
    dict_pred["HW_prediction_detection"] = file_pred["TARGETS"]["HW"]["detection_probability"]
    dict_pred["HW_prediction_marginal"]  = file_pred["TARGETS"]["HW"]["marginal_probability"]
    
    ## true
    dict_true = {}
    dict_true["t1_b"]  = file_true["TARGETS"]["t1"]["b"]
    dict_true["t1_q1"] = file_true["TARGETS"]["t1"]["q1"]
    dict_true["t1_q2"] = file_true["TARGETS"]["t1"]["q2"]
    dict_true["t2_b"]  = file_true["TARGETS"]["t2"]["b"]
    dict_true["t2_q1"] = file_true["TARGETS"]["t2"]["q1"]
    dict_true["t2_q2"] = file_true["TARGETS"]["t2"]["q2"]
    dict_true["HW_q1"] = file_true["TARGETS"]["HW"]["q1"]
    dict_true["HW_q2"] = file_true["TARGETS"]["HW"]["q2"]
    
    
    return dict_pred, dict_true
	

if __name__=="__main__":
	main()
