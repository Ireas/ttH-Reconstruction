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
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title


# CONSTANTS
INPUT_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_20-1436656e_weighted.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/more_plots/"

CLASSIFIER_VARIABLE_NAME = "neutrino_weighting/classification_event_completion"


SKIP_AUGMENTING = False


PLOT_TUPLES = [
#    (
#        "short_NW_weight", # variable name
#        np.linspace(0,1,11), # binning
#        "NW Weights (no solution excluded)", #title
#        r"Weight", #x-axis
#        False,
#    ),
#    (
#        "short_NW_wlep_mass", # variable name
#        np.linspace(0,50,11), # binning
#        "NW Prediction (no solution excluded)", #title
#        r"Mass $m_{W,\text{lep}}$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_assignment_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Assignment Probability $t_1$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_detection_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Detection Probability $t_1$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_marginal_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Margin Probability $t_1$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_assignment_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Assignment Probability $t_2$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_detection_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Detection Probability $t_2$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_marginal_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Margin Probability $t_2$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_assignment_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Assignment Probability $W_H$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_detection_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Detection Probability $W_H$", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_marginal_probability", # variable name
#        np.linspace(0,1,11), # binning
#        "SPANet Prediction", #title
#        r"Margin Probability $W_H$", #x-axis
#        False,
#    ),
#    (
#        "t1_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        "SPANet Prediction", #title
#        r"Energy $t_1$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "t2_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        "SPANet Prediction", #title
#        r"Energy $t_2$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "HW_energy", # variable name
#        np.linspace(0,1e3,21), # binning
#        "SPANet Prediction", #title
#        r"Energy $H_W$ [GeV]", #x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/number_of_jets", # variable name
#        np.arange(8,14,1), # binning
#        "Event Property", #title
#        r"Number of Jets", #x-axis
#        True,
#    ),
#    (
#        "neutrino_weighting/number_of_bjets", # variable name
#        np.arange(2,7,1), # binning
#        "Event Property", #title
#        r"Number of $b$-Jets", #x-axis
#        True,
#    ),
]

def plot_2d_hist(name, dist1, bin1, dist2, bin2, normalize="not", plot_options=None):
    # Compute the 2D histogram
    hist, xedges, yedges = np.histogram2d(dist1, dist2, bins=[bin1, bin2])
        
    # Normalize the histogram based on the provided option
    if normalize == "hist":
        hist = hist / np.sum(hist)
    elif normalize == "row":
        hist = hist / hist.sum(axis=1, keepdims=True)
    elif normalize == "col":
        hist = hist / hist.sum(axis=0, keepdims=True)
    elif normalize == "not":
        pass  # Do nothing, no normalization
    else:
        print("Warning: Invalid normalization option. No normalization applied.")
    
    # Create the plot
    fig, ax = plt.subplots()
    
    # Plot the 2D histogram
    im = ax.imshow(hist.T, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto')
    
    # Check bin size condition
    if len(bin1) < 10 and len(bin2) < 10:
        # Add labels to each bin
        for i in range(len(xedges) - 1):
            for j in range(len(yedges) - 1):
                text = ax.text((xedges[i] + xedges[i + 1]) / 2, 
                               (yedges[j] + yedges[j + 1]) / 2, 
                               f'{hist[i, j]:.2f}', 
                               ha='center', va='center', color='white')
    else:
        # Plot colorbar if bin count is large
        plt.colorbar(im, ax=ax)
    
    # Apply plot options
    if plot_options:
        if 'title' in plot_options:
            ax.set_title(plot_options['title'])
        if 'xlabel' in plot_options:
            ax.set_xlabel(plot_options['xlabel'])
        if 'ylabel' in plot_options:
            ax.set_ylabel(plot_options['ylabel'])
    
    # Set limits
    ax.set_xlim([xedges[0], xedges[-1]])
    ax.set_ylim([yedges[0], yedges[-1]])
    
    plt.savefig(f"{PLOT_DESTINATION}{name}.png")
    plt.clf()


def main():
    with uproot.open(INPUT_ROOT_FILE) as root_file:			
        nw_weights = np.array( root_file["neutrino_weighting/NW_weight"].array() )
        met_values = np.array( root_file["neutrino_weighting/reco_met_value"].array() )/1e3
        classification_event = np.array( root_file["neutrino_weighting/classification_event_channel"].array() )


        spanet_probabilities_t1_assignment = np.array( root_file["neutrino_weighting/spanet_t1_assignment_probability"].array() )
        spanet_probabilities_t1_detection = np.array( root_file["neutrino_weighting/spanet_t1_detection_probability"].array() )
        spanet_probabilities_t1_marginal = np.array( root_file["neutrino_weighting/spanet_t1_marginal_probability"].array() )
        spanet_probabilities_t2_assignment = np.array( root_file["neutrino_weighting/spanet_t2_assignment_probability"].array() )
        spanet_probabilities_t2_detection = np.array( root_file["neutrino_weighting/spanet_t2_detection_probability"].array() )
        spanet_probabilities_t2_marginal = np.array( root_file["neutrino_weighting/spanet_t2_marginal_probability"].array() )
        spanet_probabilities_HW_assignment = np.array( root_file["neutrino_weighting/spanet_HW_assignment_probability"].array() )
        spanet_probabilities_HW_detection = np.array( root_file["neutrino_weighting/spanet_HW_detection_probability"].array() )
        spanet_probabilities_HW_marginal = np.array( root_file["neutrino_weighting/spanet_HW_marginal_probability"].array() )


        short_t1_ass = []
        short_t2_ass = []
        short_nw_weight = []
        #short_t1 = np.array([])
        #short_t2 = np.array([])
        #short_HW = np.array([])
            
        for (event_channel, t1_assignment, t2_assignment, nw_weight) in zip(classification_event, spanet_probabilities_t1_assignment, spanet_probabilities_t2_assignment, nw_weights):
            if event_channel!=1:
                continue
            
            short_t1_ass.append(t1_assignment)
            short_t2_ass.append(t2_assignment)
            short_nw_weight.append(nw_weight)
            #short_t1 = np.append(short_t1, t1)
            #short_t2 = np.append(short_t2, t2)
            #short_HW = np.append(short_HW, higgs)


        # test plot
        plot_2d_hist(
            "ttbar_only_t1_assignment_vs_nw_weight", 
            short_t1_ass,
            np.linspace(0,1,11),
            short_nw_weight,
            np.linspace(0,1,11),
            normalize = "col",
            plot_options = {
                "title" : r"$t\bar{t} only",
                "xlabel" : r"$t_1$ Assignment Probability",
                "ylabel" : r"Neutrino Weighting Weight",
            }
        )
        plot_2d_hist(
            "ttbar_only_t2_assignment_vs_nw_weight", 
            short_t2_ass,
            np.linspace(0,1,11),
            short_nw_weight,
            np.linspace(0,1,11),
            normalize = "col",
            plot_options = {
                "title" : r"$t\bar{t} only",
                "xlabel" : r"$t_2$ Assignment Probability",
                "ylabel" : r"Neutrino Weighting Weight",
            }
        )
        #plot_2d_hist(
        #    "t1_assignment_vs_t2_decay", 
        #    short_t1_ass,
        #    np.linspace(0,1,6),
        #    short_t2,
        #    np.linspace(-1,4,6),
        #    normalize = True,
        #    plot_options = {
        #        "title" : r"$t\bar{t}(H\rightarrow WW)$ only",
        #        "xlabel" : r"$t_1$ Assignment Probability",
        #        "ylabel" : r"Classifier t1 Decay",
        #    }
        #)
        #plot_2d_hist(
        #    "t1_assignment_vs_H_decay", 
        #    short_t1_ass,
        #    np.linspace(0,1,6),
        #    short_HW,
        #    np.linspace(10,15,7),
        #    normalize = True,
        #    plot_options = {
        #        "title" : r"$t\bar{t}(H\rightarrow WW)$ only",
        #        "xlabel" : r"$t_1$ Assignment Probability",
        #        "ylabel" : r"Classifier Higgs Decay",
        #    }
        #)        

        # test plot
#        plot_2d_hist(
#            "nw_weight_vs_met_value", 
#            nw_weights,
#            np.linspace(0,1,5),
#            met_values,
#            np.linspace(0,100,10),
#            normalize = True,
#            plot_options = {
#                "xlabel" : "NW Weight",
#                "ylabel" : r"Missing $E_T$ [GeV]",
#            }
#        )
#        plot_2d_hist(
#            "nw_weight_vs_lepton_energies", 
#            nw_weights,
#            np.linspace(0,1,5),
#            lepton_energies,
#            np.linspace(0,100,10),
#            normalize = True,
#            plot_options = {
#                "xlabel" : "NW Weight",
#                "ylabel" : r"Lepton $E$ [GeV]",
#            }
#        )
#        plot_2d_hist(
#            "nw_weight_vs_lepton_pts", 
#            nw_weights,
#            np.linspace(0,1,5),
#            lepton_pts,
#            np.linspace(0,100,10),
#            normalize = True,
#            plot_options = {
#                "xlabel" : "NW Weight",
#                "ylabel" : r"Lepton $p_t$ [GeV]",
#            }
#        )
#        plot_2d_hist(
#            "t1_vs_t2", 
#            classification_t1,
#            np.linspace(-1,4,6),
#            classification_t2,
#            np.linspace(-1,4,6),
#            normalize = True,
#            plot_options = {
#                "xlabel" : r"classifier $t_1$",
#                "ylabel" : r"classifier $t_2$",
#            }
#        )
#        plot_2d_hist(
#            "spanet_detection_t1_vs_t2", 
#            spanet_probabilities_t1_detection,
#            np.linspace(0,1,20),
#            spanet_probabilities_t2_detection,
#            np.linspace(0,1,20),
#            normalize = True,
#            plot_options = {
#                "xlabel" : r"SPA-Net $t_1$ Detection Probabiity",
#                "ylabel" : r"SPA-Net $t_2$ Detection Probabiity",
#            }
#        )
#        plot_2d_hist(
#            "spanet_detection_t1_vs_HW", 
#            spanet_probabilities_t1_detection,
#            np.linspace(0,1,20),
#            spanet_probabilities_HW_detection,
#            np.linspace(0,1,20),
#            normalize = True,
#            plot_options = {
#                "xlabel" : r"SPA-Net $t_1$ Detection Probabiity",
#                "ylabel" : r"SPA-Net $W_\text{had}$ Detection Probabiity",
#            }
#        )
#        plot_2d_hist(
#            "spanet_assignment_t1_vs_t2", 
#            spanet_probabilities_t1_assignment,
#            np.linspace(0,1,20),
#            spanet_probabilities_t2_assignment,
#            np.linspace(0,1,20),
#            normalize = True,
#            plot_options = {
#                "xlabel" : r"SPA-Net $t_1$ Assignment Probabiity",
#                "ylabel" : r"SPA-Net $t_2$ Assignment Probabiity",
#            }
#        )
#        plot_2d_hist(
#            "spanet_assignment_t1_vs_HW", 
#            spanet_probabilities_t1_assignment,
#            np.linspace(0,1,20),
#            spanet_probabilities_HW_assignment,
#            np.linspace(0,1,20),
#            normalize = True,
#            plot_options = {
#                "xlabel" : r"SPA-Net $t_1$ Assignment Probabiity",
#                "ylabel" : r"SPA-Net $W_\text{had}$ Assignment Probabiity",
#            }
#        )
    
    exit()
    


def augment_variables():
    # open root file
    with uproot.open(INJECTED_ROOT_FILE) as root_file:     
        #### NW shortend
        # define dict and keys
        short_NW_weight = np.array([], dtype=np.float32)
        short_NW_nu_eta = np.array([], dtype=np.float32)
        short_NW_nu_phi = np.array([], dtype=np.float32)
        short_NW_wlep_mass = np.array([], dtype=np.float32)

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
            short_NW_weight = np.append(short_NW_weight, weight)
            short_NW_nu_eta = np.append(short_NW_nu_eta, nu_eta)
            short_NW_nu_phi = np.append(short_NW_nu_phi, nu_phi)
            short_NW_wlep_mass = np.append(short_NW_wlep_mass, wlep_mass) 


        ### t1,t2,HW energies
        # define dict and keys
        t1_energies = np.array([], dtype=np.float32)
        t2_energies = np.array([], dtype=np.float32)
        HW_energies = np.array([], dtype=np.float32)
        

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
            t1_energies = np.append(t1_energies, jet_energies[t1_q1] + jet_energies[t1_q2] + jet_energies[t1_b])
            t2_energies = np.append(t2_energies, jet_energies[t2_q1] + jet_energies[t2_q2] + jet_energies[t2_b])
            HW_energies = np.append(HW_energies, jet_energies[HW_q1] + jet_energies[HW_q2])


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
            if classification==1 or classification==-2:
                signal = np.append(signal, value)
            else:
                background = np.append(background, value)

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