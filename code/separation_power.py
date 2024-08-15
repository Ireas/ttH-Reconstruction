import os
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt


# CONSTANTS
INJECTED_ROOT_FILE = "/media/ireas/Data/v5/injected/all_8+j_1530766e_injected.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/separation_power/"

CLASSIFIER_VARIABLE_NAME = "matched/classification_event_completion"

SHOW_PLOTS = False

DISTRIBUTION_VARIABLE_NAMES = [
#    "matched/number_of_jets",
    "matched/number_of_bjets",
#    "matched/number_of_leptons",
    "t1_energy",
    "t2_energy",
    "HW_energy",
    "spanet/t1_marginal_probability",
    "spanet/t1_detection_probability",
    "spanet/t1_assignment_probability",
    "spanet/t2_marginal_probability",
    "spanet/t2_detection_probability",
    "spanet/t2_assignment_probability",
    "spanet/HW_marginal_probability",
    "spanet/HW_detection_probability",
    "spanet/HW_assignment_probability",
]

DISTRIBUTION_BINS = [
#    np.arange(8,15,1),
    np.arange(2,6,1),
    np.linspace(0,2e6,21),
    np.linspace(0,2e6,21),
    np.linspace(0,2e6,21),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
    np.linspace(0,1,11),
]


def main():
    # augment the dataset with custom variables
    print(f"augmenting {INJECTED_ROOT_FILE}")
    augmentations = augment_variables()
    
    print(f"calculating and plotting")
    # loop over all targeted variables
    for (variable, bins) in zip(DISTRIBUTION_VARIABLE_NAMES, DISTRIBUTION_BINS):
        # access root file
        print(f" > {variable}")
        print(f" >> binning split into {len(bins)-1} bins from {np.round(min(bins),3)} to {np.round(max(bins),3)}")
        
        # get signal and background
        signal, background = get_distributions(variable, augmentations)

        if (len(signal)==0 or len(background)==0 ):
            print(" >> signal", len(signal), signal)
            print(" >> background", len(background), background)
            print(" >> empty entries, skipping")
            continue
        
        # print true distribution
        print(f" >> distributed reaches from {np.round(min(min(signal),min(background)),3)} to {np.round(max(max(signal),max(background)),3)}")
        
        # print separation power value
        print(f" >> separation power equals {np.round(calculate_separation_power(signal, background, bins), 5)}")

        # show separation plot
        plot_separation(signal, background, variable, bins)
        print(f" >> plots produced")


def augment_variables():
    # open root file
    with uproot.open(INJECTED_ROOT_FILE) as root_file:     
        # define dict and keys
        t1_energies = np.array([], dtype=np.float32)
        t2_energies = np.array([], dtype=np.float32)
        HW_energies = np.array([], dtype=np.float32)

        # get needed varaiables
        pred_t1_q1 = np.array( root_file["spanet/t1_q1"].array() )
        pred_t1_q2 = np.array( root_file["spanet/t1_q2"].array() )
        pred_t1_b = np.array( root_file["spanet/t1_b"].array() )
        pred_t2_q1 = np.array( root_file["spanet/t2_q1"].array() )
        pred_t2_q2 = np.array( root_file["spanet/t2_q2"].array() )
        pred_t2_b = np.array( root_file["spanet/t2_b"].array() )
        pred_HW_q1 = np.array( root_file["spanet/HW_q1"].array() )
        pred_HW_q2 = np.array( root_file["spanet/HW_q2"].array() )
        all_jet_energies = root_file["matched/lvecs_jets.fCoordinates.fE"].array()

        # fill augmentations
        for (t1_q1, t1_q2, t1_b, t2_q1, t2_q2, t2_b, HW_q1, HW_q2, jet_energies) in zip(pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2, all_jet_energies):
            # sanity check
            if (len(jet_energies)<HW_q1+1) or (len(jet_energies)<HW_q2+1):
                print(f"WARNING: {len(jet_energies)}, {HW_q1}, {HW_q2}")
                continue
            
            # combine jets to resonance particle
            t1_energies = np.append(t1_energies, jet_energies[t1_q1] + jet_energies[t1_q2] + jet_energies[t1_b])
            t2_energies = np.append(t2_energies, jet_energies[t1_q1] + jet_energies[t2_q2] + jet_energies[t2_b])
            HW_energies = np.append(HW_energies, jet_energies[HW_q1] + jet_energies[HW_q2])

    # return augmentations and close root file
    augmentation = {
        "t1_energy" : t1_energies, 
        "t2_energy" : t2_energies, 
        "HW_energy" : HW_energies, 
    }

    return augmentation


def get_distributions(variable, augmentations):
    signal = np.array([])
    background = np.array([])

    with uproot.open(INJECTED_ROOT_FILE) as root_file:			
        classifier = np.array( root_file[CLASSIFIER_VARIABLE_NAME].array() )
        
        # get distribution from augments or rootfile directly 
        distribution = np.array([])
        if variable in augmentations:
            distribution = augmentations[variable]
        else:
            distribution = np.array( root_file[variable].array() )

        # filter by event type
        for (classification, value) in zip(classifier, distribution):
            if classification== 1:
                signal = np.append(signal, [value])
            else:
                background = np.append(background, [value])

    return signal, background


def transform_to_histogram(signal, background, bins, normalize=False):
    # get binning parameters
    min_value = np.min([np.min(signal), np.min(background)])
    max_value = np.max([np.max(signal), np.max(background)])

    # calculate distribution with same binning
    hist_signal, _ = np.histogram(signal, bins)
    hist_background, _ = np.histogram(background, bins)

    # normalize histogram to 1 if needed
    if normalize:
        hist_signal = hist_signal/hist_signal.sum()
        hist_background = hist_background/hist_background.sum()
    
    return hist_signal, hist_background


def calculate_separation_power(signal, background, bins):
    # create even histograms
    hist_signal, hist_background = transform_to_histogram(signal, background, bins, True)

    # calculate separation power
    separation_power = 0

    for (bin_signal, bin_background) in zip(hist_signal, hist_background):
        # check for empty bins
        if(bin_signal+bin_background==0):
            continue
        
        separation_power+= (bin_signal-bin_background)**2 / (bin_signal+bin_background)

    # correctly scale separation power for double counted events
    separation_power/= 2

    return separation_power


def plot_separation(signal, background, variable, bins):
    # create plot
    plt.figure()

    # needed weights to scale histogram to 1
    weights_signal = np.ones_like(signal)/float(len(signal))
    weights_background = np.ones_like(background)/float(len(background))

    # create histogram
    plt.hist(
        [signal, background],
        bins,
        weights=[weights_signal, weights_background], # normalize
        label=["Signal","Background"],
        histtype='step',
        stacked=False,
        fill=False
    )
    
    # beautificate
    plt.title(f"{variable.split('/')[-1]}")
    plt.legend()
    plt.ylabel("relative yield")
    
    min_value = np.min([np.min(bins), np.min(bins)])
    max_value = np.max([np.max(bins), np.max(bins)])
    plt.xlim([min_value,max_value])

    plt.savefig(PLOT_DESTINATION+f"separation_{variable.split('/')[-1]}.png")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()

if __name__=="__main__":
    main()