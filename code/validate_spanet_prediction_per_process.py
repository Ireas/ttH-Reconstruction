import h5py
import numpy as np
import matplotlib.pyplot as plt


# constants
PREDICTION_FILE = "/media/ireas/Data/v5/predicted/prediction_all_8+j_1530766e.h5"
TRUTH_FILE = "/media/ireas/Data/v5/merged_h5/all_8+j_1530766e.h5"

FULL_VALID_EVENTS_ONLY = False

# try to make a number sheme (2-gidits=(process,particle)->correct percentage)

def main():
	# open files
    file_pred = h5py.File(PREDICTION_FILE, 'r')
    file_true = h5py.File(TRUTH_FILE, 'r')
	
	# transform to dictionaries for easy access
    dict_pred, dict_true = files_to_dictionaries(file_pred, file_true)
    
    # print
    digit_dict = make_digit_dict(dict_pred, dict_true)

    for (key,value) in digit_dict.items():
        print(f"{key}:  {int(value[1])} / {int(value[0])} = {np.round(value[2],3)}")
    

    # close files
    file_pred.close()
    file_true.close()


def make_digit_dict(dict_pred, dict_true):
    # hard-coded dictionary
    digit_dict = {}

    # complete and onshell had events = fully ttH available + hadronic onshell
    digit_dict[11] = np.array([0,0,0]) #event
    digit_dict[12] = np.array([0,0,0]) #HW
    digit_dict[13] = np.array([0,0,0]) #ttbar
    digit_dict[14] = np.array([0,0,0]) #t1
    digit_dict[15] = np.array([0,0,0]) #t2
    digit_dict[16] = np.array([0,0,0]) #t1W
    digit_dict[17] = np.array([0,0,0]) #t2W

    # complete and offshell had events = fully ttH available + hadronoc offshell
    digit_dict[21] = np.array([0,0,0]) #event
    digit_dict[22] = np.array([0,0,0]) #HW
    digit_dict[23] = np.array([0,0,0]) #ttbar
    digit_dict[24] = np.array([0,0,0]) #t1
    digit_dict[25] = np.array([0,0,0]) #t2
    digit_dict[26] = np.array([0,0,0]) #t1W
    digit_dict[27] = np.array([0,0,0]) #t2W

    # incomplete events = fully ttH available, but reco jets missing
    digit_dict[31] = np.array([0,0,0]) #event
    digit_dict[32] = np.array([0,0,0]) #HW
    digit_dict[33] = np.array([0,0,0]) #ttbar
    digit_dict[34] = np.array([0,0,0]) #t1
    digit_dict[35] = np.array([0,0,0]) #t2
    digit_dict[36] = np.array([0,0,0]) #t1W
    digit_dict[37] = np.array([0,0,0]) #t2W

    # ttH->bb events
    digit_dict[41] = np.array([0,0,0]) #event
    digit_dict[42] = np.array([0,0,0]) #HW
    digit_dict[43] = np.array([0,0,0]) #ttbar
    digit_dict[44] = np.array([0,0,0]) #t1
    digit_dict[45] = np.array([0,0,0]) #t2
    digit_dict[46] = np.array([0,0,0]) #t1W
    digit_dict[47] = np.array([0,0,0]) #t2W

    # ttH->tautau events
    digit_dict[51] = np.array([0,0,0]) #event
    digit_dict[52] = np.array([0,0,0]) #HW
    digit_dict[53] = np.array([0,0,0]) #ttbar
    digit_dict[54] = np.array([0,0,0]) #t1
    digit_dict[55] = np.array([0,0,0]) #t2
    digit_dict[56] = np.array([0,0,0]) #t1W
    digit_dict[57] = np.array([0,0,0]) #t2W

    # ttbar events
    digit_dict[61] = np.array([0,0,0]) #event
    digit_dict[62] = np.array([0,0,0]) #HW
    digit_dict[63] = np.array([0,0,0]) #ttbar
    digit_dict[64] = np.array([0,0,0]) #t1
    digit_dict[65] = np.array([0,0,0]) #t2
    digit_dict[66] = np.array([0,0,0]) #t1W
    digit_dict[67] = np.array([0,0,0]) #t2W

    # ttH->qqqq events
    digit_dict[71] = np.array([0,0,0]) #event
    digit_dict[72] = np.array([0,0,0]) #HW
    digit_dict[73] = np.array([0,0,0]) #ttbar
    digit_dict[74] = np.array([0,0,0]) #t1
    digit_dict[75] = np.array([0,0,0]) #t2
    digit_dict[76] = np.array([0,0,0]) #t1W
    digit_dict[77] = np.array([0,0,0]) #t2W

    # other events
    digit_dict[81] = np.array([0,0,0]) #event
    digit_dict[82] = np.array([0,0,0]) #HW
    digit_dict[83] = np.array([0,0,0]) #ttbar
    digit_dict[84] = np.array([0,0,0]) #t1
    digit_dict[85] = np.array([0,0,0]) #t2
    digit_dict[86] = np.array([0,0,0]) #t1W
    digit_dict[87] = np.array([0,0,0]) #t2W

    testcounter = 0 

    for (
            pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2,
            true_t1_q1, true_t1_q2, true_t1_b, true_t2_q1, true_t2_q2, true_t2_b, true_HW_q1, true_HW_q2,
            classification_event_completion, classification_true_higgs_decay, classification_onshell_whad, classification_true_t1_decay, classification_true_t2_decay
        ) in zip(
            dict_pred["t1_q1"], dict_pred["t1_q2"], dict_pred["t1_b"], dict_pred["t2_q1"], dict_pred["t2_q2"], dict_pred["t2_b"], dict_pred["HW_q1"], dict_pred["HW_q2"],
            dict_true["t1_q1"], dict_true["t1_q2"], dict_true["t1_b"], dict_true["t2_q1"], dict_true["t2_q2"], dict_true["t2_b"], dict_true["HW_q1"], dict_true["HW_q2"],
            dict_true["classification_event_completion"], dict_true["classification_true_higgs_decay"], dict_true["classification_onshell_whad"], dict_true["classification_true_t1_decay"], dict_true["classification_true_t2_decay"]
        ):
            # check which resonance are possible
            possible_t1W = False #(true_t1_q1!=-1 and true_t1_q2!=-1)
            possible_t1 = classification_true_t1_decay==1
            possible_t2W = False #(true_t2_q1!=-1 and true_t2_q2!=-1)
            possible_t2 = classification_true_t2_decay==1
            possible_ttbar = (possible_t1 and possible_t2)
            possible_HW = classification_true_higgs_decay==10 or classification_true_higgs_decay==11 or classification_true_higgs_decay==12 or classification_true_higgs_decay==13 or classification_true_higgs_decay==14 
            possible_event = classification_event_completion==1

            # skip events that are not fully assigned
            ##if use_full_valid_events_only and not possible_event:
            #    continue

            classification_number = -1

            # check for ttH->WW->qqlv signature
            if classification_event_completion==1:
                # check for onshell
                if classification_onshell_whad==1 or classification_onshell_whad==2 :
                    classification_number = 10 # complete and onshell whad
                else:
                    classification_number = 20 # complete and offshell whad
            
            # check for missing reco jets
            elif classification_event_completion==-2:
                classification_number = 30 # correct but jets missing
            
            # check for ttH->bb
            elif classification_true_higgs_decay==2:
                classification_number = 40
            
            # check for ttH->tatau
            elif classification_true_higgs_decay==3:
                classification_number = 50
                
            # check for ttbar
            elif classification_true_higgs_decay==-1:
                if possible_HW:
                    print(f"Warning: particle category: {classification_event_completion}, {classification_true_higgs_decay}, {classification_onshell_whad} => {true_HW_q1},{true_HW_q2}")
                classification_number = 60

            # check for ttH->qqqq
            elif classification_true_higgs_decay==10:
                classification_number = 70
            
            else:
                classification_number = 80

                # check for wierd stuff
                #if(classification_true_higgs_decay!=0 or classification_onshell_whad!=0):
                #    print(f"Warning: No particle category for this: {classification_event_completion}, {classification_true_higgs_decay}, {classification_onshell_whad}, skipping")
                #    continue

            # check only possible particle
            if possible_t1W:
                # get current values
                total = digit_dict[classification_number+6][0]
                correct = digit_dict[classification_number+6][1]
                
                # correct assignment
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2):
                    correct+= 1
                # permutate q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1):
                    correct+= 1

                total+=1
                # update current values
                digit_dict[classification_number+6] = np.array([total, correct, correct/total])

            if possible_t1:
                # get current values
                total = digit_dict[classification_number+4][0]
                correct = digit_dict[classification_number+4][1]
                
                # correct assignment
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b):
                    correct+= 1
                # permutate q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b):
                    correct+= 1

                total+=1
                # update current values
                digit_dict[classification_number+4] = np.array([total, correct, correct/total])
            
            
            if possible_t2W:
                # get current values
                total = digit_dict[classification_number+7][0]
                correct = digit_dict[classification_number+7][1]

                # correct assignment
                if(pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2):
                    correct+= 1
                # permutate q1,q2
                elif(pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1):
                    correct+= 1
                
                total+=1
                # update current values
                digit_dict[classification_number+7] = np.array([total, correct, correct/total])
            
            
            if possible_t2:
                # get current values
                total = digit_dict[classification_number+5][0]
                correct = digit_dict[classification_number+5][1]

                # correct assignment
                if(pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    correct+= 1
                # permutate q1,q2
                elif(pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    correct+= 1
                
                total+=1
                # update current values
                digit_dict[classification_number+5] = np.array([total, correct, correct/total])            

            if possible_ttbar:
                # get current values
                total = digit_dict[classification_number+3][0]
                correct = digit_dict[classification_number+3][1]

                # correct assignment
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    correct+= 1
                # permutate t1:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b):
                    correct+= 1
                # permutate t2:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b):
                    correct+= 1
                
                total+=1
                # update current values
                digit_dict[classification_number+3] = np.array([total, correct, correct/total])
           

            if possible_HW:
                # get current values
                total = digit_dict[classification_number+2][0]
                correct = digit_dict[classification_number+2][1]

                # correct assignment
                if(pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    correct+= 1
                # permutate q1,q2
                elif(pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    correct+= 1
                
                total+=1
                # update current values
                digit_dict[classification_number+2] = np.array([total, correct, correct/total])

            if possible_event:
                # get current values
                total = digit_dict[classification_number+1][0]
                correct = digit_dict[classification_number+1][1]

                # correct assignment
                if(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    correct+= 1
                # permutate t1:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    correct+= 1
                # permutate t2:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q1 and pred_HW_q2==true_HW_q2):
                    correct+= 1
                # permutate HW:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    correct+= 1
                # permutate t1:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q1 and pred_t2_q2==true_t2_q2 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    correct+= 1
                # permutate t2:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q1 and pred_t1_q2==true_t1_q2 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    correct+= 1
                # permutate t1:q1,q2 and t2:q1,q2 and HW:q1,q2
                elif(pred_t1_q1==true_t1_q2 and pred_t1_q2==true_t1_q1 and pred_t1_b==true_t1_b and pred_t2_q1==true_t2_q2 and pred_t2_q2==true_t2_q1 and pred_t2_b==true_t2_b and pred_HW_q1==true_HW_q2 and pred_HW_q2==true_HW_q1):
                    correct+= 1
                
                total+=1
                # update current values
                digit_dict[classification_number+1] = np.array([total, correct, correct/total])

    print(testcounter)    
    return digit_dict


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

    # event classifiers
    dict_true["classification_event_completion"] = file_true["OTHER"]["classification_event_completion"]
    dict_true["classification_true_higgs_decay"] = file_true["OTHER"]["classification_true_higgs_decay"]
    dict_true["classification_true_t1_decay"] = file_true["OTHER"]["classification_true_t1_decay"]
    dict_true["classification_true_t2_decay"] = file_true["OTHER"]["classification_true_t2_decay"]
    dict_true["classification_onshell_whad"] = file_true["OTHER"]["classification_onshell_whad"]
    

    return dict_pred, dict_true
	

if __name__=="__main__":
	main()
