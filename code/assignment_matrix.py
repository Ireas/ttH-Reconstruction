import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt


#plt.rc('font', size=20)          # controls default text sizes
#plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
#plt.rc('legend', fontsize=14)    # legend fontsize
#plt.rc('figure', titlesize=14)  # fontsize of the figure title

FILE_PREDICTION = "/media/ireas/Data/v5/predicted/prediction_all_8+j_1530766e.h5"
FILE_TRUTH = "/media/ireas/Data/v5/merged_h5/all_8+j_1530766e.h5"

def main():
    # open files
    pred_file = h5py.File(FILE_PREDICTION, 'r')
    true_file = h5py.File(FILE_TRUTH, 'r')
    
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

    event_status = np.array( true_file["OTHER"]["classification_event_completion"][()])
    t1_decays = np.array( true_file["OTHER"]["classification_true_t1_decay"][()])
    t2_decays = np.array( true_file["OTHER"]["classification_true_t2_decay"][()])
    higgs_decays = np.array( true_file["OTHER"]["classification_true_higgs_decay"][()])


    # create arrays with fixed positions (t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2, invalid)
    row_t1q1 = calculate_row(pred_t1q1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_t1q2 = calculate_row(pred_t1q2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_t1b = calculate_row(pred_t1b, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_t2q1 = calculate_row(pred_t2q1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_t2q2 = calculate_row(pred_t2q2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_t2b = calculate_row(pred_t2b, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_HWq1 = calculate_row(pred_HWq1, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)
    row_HWq2 = calculate_row(pred_HWq2, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays)


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

    plt.figure(figsize=(8,5))
    plt.pcolormesh(
        np.arange(-0.5, matrix_norm_row.shape[1]),
        np.arange(-0.5, matrix_norm_row.shape[0]),
        matrix_norm_row,
    )

    for i in range(matrix_norm_row.shape[0]):
        for j in range(matrix_norm_row.shape[1]):
            if matrix_norm_row[i,j]>0:
                plt.text(j, i, np.round(matrix_norm_row[i,j], 2), ha='center', va='center')

    plt.xlabel("true labels")
    plt.xticks(
        [0,1,2,3,4,5,6,7],
        labels=[r"$t_{1,q1}$",r"$t_{1,q2}$",r"$t_{1,b}$",r"$t_{2,q1}$",r"$t_{2,q2}$",r"$t_{2,b}$",r"$W_{\text{had},q1}$",r"$W_{\text{had},q2}$"]
    )
    plt.ylabel("predicted labels")
    plt.yticks(
        [0,1,2,3,4,5,6,7],
        labels=[r"$t_{1,q1}$",r"$t_{1,q2}$",r"$t_{1,b}$",r"$t_{2,q1}$",r"$t_{2,q2}$",r"$t_{2,b}$",r"$W_{\text{had},q1}$",r"$W_{\text{had},q2}$"]
    )

    plt.savefig("/home/ireas/git_repos/master/plots/assignment_matrix/assignment_ttbar_only.png")
    #plt.show()



def calculate_row(pred_indicies, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays):
    # create array with fixed positions (t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2, invalid)
    fixed_row = np.array([0,0,0,0,0,0,0,0,0], dtype=np.intc)

    # fill array
    for (pred_index, t1q1, t1q2, t1b, t2q1, t2q2, t2b, HWq1, HWq2, status, t1_decay, t2_decay, HW_decay) in zip(pred_indicies, true_t1q1, true_t1q2, true_t1b, true_t2q1, true_t2q2, true_t2b, true_HWq1, true_HWq2, event_status, t1_decays, t2_decays, higgs_decays):
        # event is invalid
        if (t1q1==0 and t1q2==0 and t1b==0 and t2q1==0 and t2q2==0 and t2b==0 and HWq1==0 and HWq2==0):
            continue
        
        # skip everything with H for ttbar matrix
        if (HW_decay!=-1):
            continue


        # only ttH->WW->qqlv?0
        #if status!=1:
        #    continue

        # ignore incomplete events?
        #if (t1q1==-1 or t1q2==-1 or t1b==-1 or t2q1==-1 or t2q2==-1 or t2b==-1 or HWq1==-1 or HWq2==-1):
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