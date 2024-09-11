import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt


# Options
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


# Constants
CONFIDENCE_PROBABILITY = "assignment"

FULL_VALID_EVENTS_ONLY = False
SHOW_PLOTS = False

PREDICTION_FILE = "/media/ireas/Data/v6/predicted/all_5+j_truncated_20_ttHWW_amplified_predicted.h5"
TRUTH_FILE = "/media/ireas/Data/v6/merged_h5/all_5+j_truncated_20_ttHWW_amplified_2503624e.h5"


# Functions
def main():
	# open files
    file_pred = h5py.File(PREDICTION_FILE, 'r')
    file_true = h5py.File(TRUTH_FILE, 'r')
	
	# transform to dictionaries for easy access
    dict_pred, dict_true = files_to_dictionaries(file_pred, file_true)
    
    #validate_all_events(dict_pred, dict_true, FULL_VALID_EVENTS_ONLY)
    
    # plot as baptiste 
    plot_t1(dict_pred, dict_true, "assignment", FULL_VALID_EVENTS_ONLY, True)
    plot_t1(dict_pred, dict_true, "detection", FULL_VALID_EVENTS_ONLY, True)
    plot_t1(dict_pred, dict_true, "marginal", FULL_VALID_EVENTS_ONLY, True)
    plot_t2(dict_pred, dict_true, "assignment", FULL_VALID_EVENTS_ONLY, True)
    plot_t2(dict_pred, dict_true, "detection", FULL_VALID_EVENTS_ONLY, True)
    plot_t2(dict_pred, dict_true, "marginal", FULL_VALID_EVENTS_ONLY, True)

    validate_all_events(dict_pred, dict_true, False)

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
            pred_confidence_t1, pred_confidence_t2, pred_confidence_HW, event_status
        ) in zip(
            dict_pred["t1_q1"], dict_pred["t1_q2"], dict_pred["t1_b"], dict_pred["t2_q1"], dict_pred["t2_q2"], dict_pred["t2_b"], dict_pred["HW_q1"], dict_pred["HW_q2"],
            dict_true["t1_q1"], dict_true["t1_q2"], dict_true["t1_b"], dict_true["t2_q1"], dict_true["t2_q2"], dict_true["t2_b"], dict_true["HW_q1"], dict_true["HW_q2"],
            dict_pred["t1_prediction_"+variable_name], dict_pred["t2_prediction_"+variable_name], dict_pred["HW_prediction_"+variable_name], dict_true["classification_event_completion"]
        ):
            
            if event_status!=1:
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
        plt.title("SPANet Prediction")
        plt.xlabel(f"$t_1$ {variable_name[0].upper() + variable_name[1:]} Probability")
        plt.ylabel("Relative Event Yield")
        plt.xlim([0,1])
        plt.legend()

        # output
        plt.savefig(f"/home/ireas/git_repos/master/plots/baptiste_plot/spanet_probability_{variable_name}_complete_normed.png")
    else:
        # plot histogram
        plt.figure()

	    # histograms
        plt.hist(confidence_correct, bins, label="correct", fill=False, histtype="step")
        plt.hist(confidence_swapped, bins, label="swapped", fill=False, histtype="step")
        plt.hist(confidence_failed, bins, label="failed", fill=False, histtype="step")
        plt.hist(confidence_impossible, bins, label="impossible", fill=False, histtype="step")

	    # labels
        plt.title("SPANet Prediction")
        plt.xlabel(f"$t_1$ {variable_name[0].upper() + variable_name[1:]} Probability")
        plt.ylabel("Event Yield")
        plt.xlim([0,1])
        plt.yscale("log")
        plt.legend()

        # output
        plt.savefig(f"/home/ireas/git_repos/master/plots/baptiste_plot/spanet_probability_{variable_name}_complete.png")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.clf()


def plot_t2(dict_pred, dict_true, variable_name, use_full_valid_events_only, normalize): 
    # transform    
    confidence_correct = np.array([])
    confidence_swapped = np.array([])
    confidence_failed = np.array([])
    confidence_impossible = np.array([])
   

    for (
            pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
            true_t1_q1, true_t1_q2, true_t1_b, true_t2_q1, true_t2_q2, true_t2_b, true_HW_q1, true_HW_q2,
            pred_confidence_t1, pred_confidence_t2, pred_confidence_HW, event_status
        ) in zip(
            dict_pred["t1_q1"], dict_pred["t1_q2"], dict_pred["t1_b"], dict_pred["t2_q1"], dict_pred["t2_q2"], dict_pred["t2_b"], dict_pred["HW_q1"], dict_pred["HW_q2"],
            dict_true["t1_q1"], dict_true["t1_q2"], dict_true["t1_b"], dict_true["t2_q1"], dict_true["t2_q2"], dict_true["t2_b"], dict_true["HW_q1"], dict_true["HW_q2"],
            dict_pred["t1_prediction_"+variable_name], dict_pred["t2_prediction_"+variable_name], dict_pred["HW_prediction_"+variable_name], dict_true["classification_event_completion"]
        ):
            
            if event_status!=1:
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
            
            if possible_t2:
                if(pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    confidence_correct = np.append(confidence_correct, [pred_confidence_t2])
                # permutate q1,q2
                elif(pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    confidence_swapped = np.append(confidence_swapped, [pred_confidence_t2])
                else:
                    confidence_failed = np.append(confidence_failed, [pred_confidence_t2])
            else: 
                confidence_impossible = np.append(confidence_impossible, [pred_confidence_t2])

 
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
        plt.title("SPANet Prediction")
        plt.xlabel(f"$t_2$ {variable_name[0].upper() + variable_name[1:]} Probability")
        plt.ylabel("Relative Event Yield")
        plt.xlim([0,1])
        plt.legend()

        # output
        plt.savefig(f"/home/ireas/git_repos/master/plots/baptiste_plot/spanet_probability_{variable_name}_complete_t2_normed.png")
    else:
        # plot histogram
        plt.figure()

	    # histograms
        plt.hist(confidence_correct, bins, label="correct", fill=False, histtype="step")
        plt.hist(confidence_swapped, bins, label="swapped", fill=False, histtype="step")
        plt.hist(confidence_failed, bins, label="failed", fill=False, histtype="step")
        plt.hist(confidence_impossible, bins, label="impossible", fill=False, histtype="step")

	    # labels
        plt.title("SPANet Prediction")
        plt.xlabel(f"$t_2$ {variable_name[0].upper() + variable_name[1:]} Probability")
        plt.ylabel("Event Yield")
        plt.xlim([0,1])
        plt.yscale("log")
        plt.legend()

        # output
        plt.savefig(f"/home/ireas/git_repos/master/plots/baptiste_plot/spanet_probability_{variable_name}_complete_t2.png")

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
    
    dict_true["classification_event_completion"] = file_true["OTHER"]["classification_event_completion"]
    
    return dict_pred, dict_true
	

if __name__=="__main__":
	main()
