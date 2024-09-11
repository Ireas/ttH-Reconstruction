import os
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt


#plt.rc('font', size=20)          # controls default text sizes
#plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=11)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title


# CONSTANTS
INPUT_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_20-1436656e_weighted.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/separation_power/"
CLASSIFIER_VARIABLE_NAME = "neutrino_weighting/classification_event_completion"

SHOW_PLOTS = False
SKIP_AUGMENTING = True


PLOT_TUPLES = [
    (
        "neutrino_weighting/reco_met_value", # variable name
        np.linspace(0,120e3,13), # binning
        "Event Property", # title
        r"Missing Transverse Energy in [MeV]", #x-axis
        False,
    ),
    (
        "neutrino_weighting/reco_lepton_e", # variable name
        np.linspace(0,160e3,17), # binning
        "Lepton Property", # title
        r"Lepton Energy [MeV]", #x-axis
        False,
    ),
#    (
#        "short_NW_weight", # variable name
#        np.linspace(0,1,11), # binning
#        "NW Weights (no solution excluded)", # title
#        r"NW Weight", #x-axis
#        False,
#    ),
    (
        "neutrino_weighting/NW_weight", # variable name
        np.linspace(-1,1,21), # binning
        "NW Weights ", # title
        r"NW Weight", #x-axis
        False,
    ),
#    (
#        "short_NW_weight_sum", # variable name
#        np.linspace(0,1000,11), # binning
#        "Sum of NW Weights (no solution excluded)", # title
#        r"Sum of NW Weight", #x-axis
#        False,
#    ),
#    (
#        "short_NW_wlep_mass", # variable name
#        np.linspace(0,50,11), # binning
#        "NW Prediction (no solution excluded)", # title
#        r"Mass $m_{W,\text{lep}}$ [GeV]", #x-axis
#        False,
#    ),
    (
        "neutrino_weighting/spanet_t1_assignment_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Assignment Probability $t_1$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_t1_detection_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Detection Probability $t_1$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_t1_marginal_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Margin Probability $t_1$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_t2_assignment_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Assignment Probability $t_2$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_t2_detection_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Detection Probability $t_2$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_t2_marginal_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Margin Probability $t_2$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_HW_assignment_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Assignment Probability $W_H$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_HW_detection_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Detection Probability $W_H$", #x-axis
        False,
    ),
    (
        "neutrino_weighting/spanet_HW_marginal_probability", # variable name
        np.linspace(0,1,21), # binning
        "SPANet Prediction", # title
        r"Margin Probability $W_H$", #x-axis
        False,
    ),
#    (
#        "t1_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        "SPANet Prediction", # title
#        r"Sum of Jet Energies $t_1$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "t2_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        "SPANet Prediction", # title
#        r"Sum of Jet Energies $t_2$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "HW_energy", # variable name
#        np.linspace(0,1e3,21), # binning
#        "SPANet Prediction", # title
#        r"Sum of Jet Energies $W_\text{had.}$ [GeV]", #x-axis
#        False,
#    ),
    (
        "neutrino_weighting/number_of_jets", # variable name
        np.arange(5,12,1), # binning
        "Event Property", # title
        r"Number of Jets", #x-axis
        True,
    ),
    (
        "neutrino_weighting/number_of_bjets", # variable name
        np.arange(2,7,1), # binning
        "Event Property", # title
        r"Number of $b$-Jets", #x-axis
        True,
    ),
]


def main():
    # augment the dataset with custom variables
    augmentations = {}
    if not SKIP_AUGMENTING:
        print(f"augmenting {INPUT_ROOT_FILE}")
        augmentations = augment_variables()
    

    print(f"calculating and plotting")
    # loop over all targeted variables
    for target_tuple in PLOT_TUPLES:
        # access truple
        variable = target_tuple[0]
        bins = target_tuple[1]
        
        print(f" > {variable}")
        

        # ity check
        if len(bins-1)==0:
            print(f" >> Error: bins are empty, skipping")
            continue
    

        # access root file
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
        plot_separation(signal, background, target_tuple)
        print(f" >> plots produced")
        print()


def augment_variables():
    # open root file
    with uproot.open(INPUT_ROOT_FILE) as root_file:     
        #### NW shortend
        # define dict and keys
        short_NW_weight = []
        short_NW_nu_eta = []
        short_NW_nu_phi = []
        short_NW_wlep_mass = []

        # get needed varaiables
        NW_weight = np.array( root_file["neutrino_weighting/NW_weight"].array() )
        NW_nu_eta = np.array( root_file["neutrino_weighting/NW_nu_eta"].array() )
        NW_nu_phi = np.array( root_file["neutrino_weighting/NW_nu_phi"].array() )
        NW_wlep_mass = np.array( root_file["neutrino_weighting/NW_wlep_mass"].array() )
       
        # fill augmentations
        for (weight, nu_eta, nu_phi, wlep_mass) in zip(NW_weight, NW_nu_eta, NW_nu_phi, NW_wlep_mass):
            # check prediction
            if weight<0:
                continue
            
            # fill shorted stuff
            short_NW_weight.append(weight)
            short_NW_nu_eta.append(nu_eta)
            short_NW_nu_phi.append(nu_phi)
            short_NW_wlep_mass.append(wlep_mass) 

        short_NW_weight = np.array(short_NW_weight)
        short_NW_nu_eta = np.array(short_NW_nu_eta)
        short_NW_nu_phi = np.array(short_NW_nu_phi)
        short_NW_wlep_mass = np.array(short_NW_wlep_mass)

        ### t1,t2,HW energies
        # define dict and keys
        t1_energies = []
        t2_energies = []
        HW_energies = []
        

        # get needed varaiables
        pred_t1_q1 = np.array( root_file["neutrino_weighting/spanet_t1_q1"].array() )
        pred_t1_q2 = np.array( root_file["neutrino_weighting/spanet_t1_q2"].array() )
        pred_t1_b = np.array( root_file["neutrino_weighting/spanet_t1_b"].array() )
        pred_t2_q1 = np.array( root_file["neutrino_weighting/spanet_t2_q1"].array() )
        pred_t2_q2 = np.array( root_file["neutrino_weighting/spanet_t2_q2"].array() )
        pred_t2_b = np.array( root_file["neutrino_weighting/spanet_t2_b"].array() )
        pred_HW_q1 = np.array( root_file["neutrino_weighting/spanet_HW_q1"].array() )
        pred_HW_q2 = np.array( root_file["neutrino_weighting/spanet_HW_q2"].array() )
        all_jet_energies = root_file["neutrino_weighting/jet_e_NOSYS"].array()

        # fill augmentations
        for (t1_q1, t1_q2, t1_b, t2_q1, t2_q2, t2_b, HW_q1, HW_q2, jet_energies) in zip(pred_t1_q1, pred_t1_q2, pred_t1_b, pred_t2_q1, pred_t2_q2, pred_t2_b, pred_HW_q1, pred_HW_q2, all_jet_energies):
            # sanity check
            if (len(jet_energies)<HW_q1+1) or (len(jet_energies)<HW_q2+1):
                print(f"WARNING: {len(jet_energies)}, {HW_q1}, {HW_q2}")
                continue
            
            # combine jets to resonance particle
            t1_energies.append( jet_energies[t1_q1] + jet_energies[t1_q2] + jet_energies[t1_b] )
            t2_energies.append( jet_energies[t2_q1] + jet_energies[t2_q2] + jet_energies[t2_b] )
            HW_energies.append( jet_energies[HW_q1] + jet_energies[HW_q2] )


    t1_energies = np.array(t1_energies)
    t2_energies = np.array(t2_energies)
    HW_energies = np.array(HW_energies)

    # return augmentations and close root file
    augmentation = {
        "t1_energy" : t1_energies/1e3, #scale to MeV 
        "t2_energy" : t2_energies/1e3, #scale to MeV
        "HW_energy" : HW_energies/1e3, #scale to MeV
        "short_NW_weight" : short_NW_weight,
        "short_NW_nu_eta" : short_NW_nu_eta,
        "short_NW_nu_phi" : short_NW_nu_phi,
        "short_NW_wlep_mass" : short_NW_wlep_mass/1e3, #scale to MeV 
    }

    return augmentation


def get_distributions(variable, augmentations):
    signal = []
    background = []

    with uproot.open(INPUT_ROOT_FILE) as root_file:			
        classifier = np.array( root_file[CLASSIFIER_VARIABLE_NAME].array() )
        
        # get distribution from augments or rootfile directly 
        distribution = np.array([])
        if variable in augmentations:
            distribution = augmentations[variable]
        else:
            distribution = np.array( root_file[variable].array() )

        # filter by event type
        for (classification, value) in zip(classifier, distribution):
            if classification==1 or classification==-2:
                signal.append(value)
            else:
                background.append(value)

    signal = np.array(signal)
    background = np.array(background)

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


def plot_separation(signal, background, target_tuple):
    # access target tuple
    variable = target_tuple[0]
    bins = target_tuple[1]
    title = target_tuple[2]
    x_label = target_tuple[3]
    is_int_value = target_tuple[4]

    
    # get correct values
    min_value = np.min([np.min(bins), np.min(bins)])
    max_value = np.max([np.max(bins), np.max(bins)])
    
    # needed weights to scale histogram to 1
    weights_signal = np.ones_like(signal)/float(len(signal))
    weights_background = np.ones_like(background)/float(len(background))

    
    # create plot
    plt.figure()

    # int value centered bins
    if is_int_value:
        bins = bins - 0.5
    

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
    plt.title(title)    
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel("Relative Event Yield")

    
    plt.xlim([min_value,max_value])

    # nice int ticks 
    if is_int_value:
        plt.xlim([min_value-0.5,max_value-0.5])
        plt.xticks(range(int(min_value), int(max_value)))

    # save figure
    plt.savefig(PLOT_DESTINATION+f"separation_{variable.split('/')[-1]}.png")
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()

if __name__=="__main__":
    main()