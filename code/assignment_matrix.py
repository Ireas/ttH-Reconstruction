import numpy as np
import uproot
import h5py
import matplotlib.pyplot as plt


# OPTIONS
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


# CONSTANTS
INPUT_ROOT_FILE = "/media/ireas/Data/v6/injected/all_5+j_truncated_50_ttHWW_amplified_1675331e_injected.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/assignment_matrix/"

PREDICTION_FILE = "/media/ireas/Data/v6/predicted/all_5+j_truncated_50_ttHWW_amplified_1675331e_predicted.h5"
TRUTH_FILE = "/media/ireas/Data/v6/merged_h5/all_5+j_truncated_50_ttHWW_amplified_1675331e.h5"


MATRIX_LABELS = [
    r"$t_1 q_1$",
    r"$t_1 q_2$",
    r"$t_1 b$",
    r"$t_2 q_1$",
    r"$t_2 q_2$",
    r"$t_2 b$",
    r"$W_\text{had} q_1$",
    r"$W_\text{had} q_2$",
]



def main():
    # open files
    pred_file = h5py.File(PREDICTION_FILE, 'r')
    true_file = h5py.File(TRUTH_FILE, 'r')
    
    # access predicted assignments
    pred_t1q1 = np.array( pred_file['TARGETS']['t1']['q1'][()] )
    pred_t1q2 = np.array( pred_file['TARGETS']['t1']['q2'][()] )
    pred_t1b = np.array( pred_file['TARGETS']['t1']['b'][()] )
    pred_t2q1 = np.array( pred_file['TARGETS']['t2']['q1'][()] )
    pred_t2q2 = np.array( pred_file['TARGETS']['t2']['q2'][()] )
    pred_t2b = np.array( pred_file['TARGETS']['t2']['b'][()] )
    pred_HWq1 = np.array( pred_file['TARGETS']['HW']['q1'][()] )
    pred_HWq2 = np.array( pred_file['TARGETS']['HW']['q2'][()] )
    
    # access true assignments
    true_t1q1 = np.array( true_file['TARGETS']['t1']['q1'][()] )
    true_t1q2 = np.array( true_file['TARGETS']['t1']['q2'][()] )
    true_t1b = np.array( true_file['TARGETS']['t1']['b'][()] )
    true_t2q1 = np.array( true_file['TARGETS']['t2']['q1'][()] )
    true_t2q2 = np.array( true_file['TARGETS']['t2']['q2'][()] )
    true_t2b = np.array( true_file['TARGETS']['t2']['b'][()] )
    true_HWq1 = np.array( true_file['TARGETS']['HW']['q1'][()] )
    true_HWq2 = np.array( true_file['TARGETS']['HW']['q2'][()] )

    event_channel = False#np.array( true_file["OTHER"]["classification_event_channel"][()])
    classification_hw = np.array( true_file['OTHER']['classification_true_higgs_decay'][()] )
    classification_t1 = np.array( true_file['OTHER']['classification_true_t1_decay'][()] )
    classification_t2 = np.array( true_file['OTHER']['classification_true_t2_decay'][()] )

    # create arrays with fixed positions (t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2, invalid)
    row_t1q1 = calculate_row(pred_t1q1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_t1q2 = calculate_row(pred_t1q2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_t1b = calculate_row(pred_t1b, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_t2q1 = calculate_row(pred_t2q1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_t2q2 = calculate_row(pred_t2q2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_t2b = calculate_row(pred_t2b, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_HWq1 = calculate_row(pred_HWq1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)
    row_HWq2 = calculate_row(pred_HWq2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2)


    # print output
    print()
    print("==========  OUTPUT  ==========")
    print("pred\t|| t1q1\t| t1q2 \t| t1b \t| t2q1 \t| t2q2 \t| t1b \t| HWq1 \t| HWq2 \t|| invalid")
    print("----------------------------------------------------------------------------------")
    print_row("t1q1", row_t1q1)
    print_row("t1q2", row_t1q2)
    print_row("t1b ", row_t1b)
    print_row("t2q1", row_t2q1)
    print_row("t2q2", row_t2q2)
    print_row("t2b ", row_t2b)
    print_row("HWq1", row_HWq1)
    print_row("HWq2", row_HWq2)
    print("----------------------------------------------------------------------------------")
    print_row("t1q1", row_t1q1, False)
    print_row("t1q2", row_t1q2, False)
    print_row("t1b ", row_t1b, False)
    print_row("t2q1", row_t2q1, False)
    print_row("t2q2", row_t2q2, False)
    print_row("t2b ", row_t2b, False)
    print_row("HWq1", row_HWq1, False)
    print_row("HWq2", row_HWq2, False)
    print("----------------------------------------------------------------------------------")


    # plot
    matrix = np.array([row_t1q1, row_t1q2, row_t1b, row_t2q1, row_t2q2, row_t2b, row_HWq1, row_HWq2])[:,:-1]


    matrix_norm_row = np.array([
        row_t1q1[:-1]/matrix[:,0].sum(), 
        row_t1q2[:-1]/matrix[:,1].sum(), 
        row_t1b[:-1]/matrix[:,2].sum(), 
        row_t2q1[:-1]/row_t2q1[:-1].sum(), 
        row_t2q2[:-1]/row_t2q2[:-1].sum(), 
        row_t2b[:-1]/row_t2b[:-1].sum(),
        row_HWq1[:-1]/row_HWq1[:-1].sum(),
        row_HWq2[:-1]/row_HWq2[:-1].sum()
    ])

    plt.figure(figsize=(10,8))
    plt.pcolormesh(
        np.arange(-0.5, matrix_norm_row.shape[1]),
        np.arange(-0.5, matrix_norm_row.shape[0]),
        matrix_norm_row,
    )

    for i in range(matrix_norm_row.shape[0]):
        for j in range(matrix_norm_row.shape[1]):
            if matrix_norm_row[i,j]>0:
                plt.text(j, i, np.round(matrix_norm_row[i,j], 2), ha='center', va='center')

    plt.xlabel("True Labels")
    plt.xticks(np.arange(len(MATRIX_LABELS)),labels=MATRIX_LABELS)
    plt.ylabel("Predicted Labels")
    plt.yticks(np.arange(len(MATRIX_LABELS)),labels=MATRIX_LABELS)

    plt.savefig(f"{PLOT_DESTINATION}assignment_matrix.png")
    #plt.show()



def calculate_row(pred_indicies, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, classification_hw, classification_t1, classification_t2):
    # create array with fixed positions (t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2, invalid)
    fixed_row = np.array([0,0,0,0,0,0,0,0,0], dtype=np.intc)

    # fill array
    for (class_hw, class_t1, class_t2, pred_index, t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2) in zip(classification_hw, classification_t1, classification_t2, pred_indicies, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2):
        # event is invalid
        if (t1q1==0 and t1q2==0 and t1b==0 and t2q1==0 and t2q2==0 and t2b==0 and HWq1==0 and HWq2==0):
            continue
    
        # use only events where all particles are included
        if(t1q1==-1 or t1q2==-1 or t1b==-1 or t2q1==-1 or t2q2==-1 or t2b==-1 or HWq1==-1 or HWq2==-1):
            continue
            
        #if not ( (class_hw==11 or class_hw==12 or class_hw==13 or class_hw==14) and class_t1==1 and class_t2==1):
        #    continue 
        # ttHWW events only
        #if event_channel<10:
        #    continue
        
        # predicted particle is t1q1
        elif pred_index==t1q1:
            fixed_row[0]+= 1
        # predicted particle is t1q2
        elif pred_index==t1q2:
            fixed_row[1]+= 1
        # predicted particle is t1b
        elif pred_index==t1b:
            fixed_row[2]+= 1
        # predicted particle is t2q1
        elif pred_index==t2q1:
            fixed_row[3]+= 1
        # predicted particle is t2q1
        elif pred_index==t2q2:
            fixed_row[4]+= 1
        # predicted particle is t2b
        elif pred_index==t2b:
            fixed_row[5]+= 1
        # predicted particle is HWq1
        elif pred_index==HWq1:
            fixed_row[6]+= 1
        # predicted particle is HWq2
        elif pred_index==HWq2:
            fixed_row[7]+= 1
        # predicted particle is invalid
        else:
            fixed_row[8]+= 1

    return fixed_row



def print_row(title, row, relative=True):
    if relative:
        print(f"{title}\t|| {np.round(row[0]/row.sum(),2)} \t| {np.round(row[1]/row.sum(),2)} \t| {np.round(row[2]/row.sum(),2)} \t| {np.round(row[3]/row.sum(),2)} \t| {np.round(row[4]/row.sum(),2)} \t| {np.round(row[5]/row.sum(),2)} \t| {np.round(row[6]/row.sum(),2)} \t| {np.round(row[7]/row.sum(),2)} \t|| {np.round(row[8]/row.sum(),2)}")
    else:
        print(f"{title}\t|| {row[0]} \t| {row[1]} \t| {row[2]} \t| {row[3]} \t| {row[4]} \t| {row[5]} \t| {row[6]} \t| {row[7]} \t|| {row[8]}")


if __name__=="__main__":
    main()