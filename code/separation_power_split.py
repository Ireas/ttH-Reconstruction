import os
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt


plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


# CONSTANTS
INJECTED_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_50_ttHWW_amplified_1675331e_weighted.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/separation_power/"
CLASSIFIER_VARIABLE_NAME = "neutrino_weighting/classification_event_channel"
SHOW_PLOTS = False
SKIP_AUGMENTING = False


PLOT_TUPLES = [
    (
        "neutrino_weighting/NW_H_mass", # variable name
        np.arange(105e3,145e3,1e3), # binning
        r"Mass of $H$-boson $m_H$ [MeV]", # x-axis
        False,
    ),
    (
        "short_NW_weight", # variable name
        np.linspace(0,1,11), # binning
        r"NW Weight", # x-axis
        False,
    ),
    (
        "neutrino_weighting/NW_weight", # variable name
        np.linspace(-1,1,21), # binning
        r"NW Weight", # x-axis
        False,
    ),
#    (
#        "short_NW_wlep_mass", # variable name
#        np.linspace(0,50,11), # binning
#        r"Mass $m_{W,\text{lep}}$ [GeV]", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_assignment_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Assignment Probability $t_1$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_detection_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Detection Probability $t_1$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t1_marginal_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Margin Probability $t_1$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_assignment_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Assignment Probability $t_2$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_detection_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Detection Probability $t_2$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_t2_marginal_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Margin Probability $t_2$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_assignment_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Assignment Probability $W_H$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_detection_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Detection Probability $W_H$", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/spanet_HW_marginal_probability", # variable name
#        np.linspace(0,1,21), # binning
#        r"Margin Probability $W_H$", # x-axis
#        False,
#    ),
#    (
#        "t1_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        r"Sum of Jet Energies $t_1$ [GeV]", # x-axis
#        False,
#    ),
#    (
#        "t2_energy", # variable name
#        np.linspace(0,2e3,21), # binning
#        r"Sum of Jet Energies $t_2$ [GeV]", # x-axis
#        False,
#    ),
#    (
#        "neutrino_weighting/number_of_jets", # variable name
#        np.arange(5,12,1), # binning
#        r"Number of Jets", # x-axis
#        True,
#    ),
#    (
#        "neutrino_weighting/number_of_bjets", # variable name
#        np.arange(2,7,1), # binning
#        r"Number of $b$-Jets", # x-axis
#        True,
#    ),
#    (
#        "reco_met_value_gev", # variable name
#        np.linspace(0,100,11), # binning
#        r"Missing Transverse Energy [GeV]", # x-axis
#        False,
#    ),
#    (
#        "reco_lepton_e_gev", # variable name
#        np.linspace(0,160,11), # binning
#        r"Energy of the Lepton [GeV]", # x-axis
#        False
#    ),
#    (
#        "reco_lepton_pt_gev", # variable name
#        np.linspace(0,120,11), # binning
#        r"Transverse Momentum of the Lepton [GeV]", # x-axis
#        False
#    ),
]


def main():
    # augment the dataset with custom variables
    augmentations = {}
    if not SKIP_AUGMENTING:
        print(f"augmenting {INJECTED_ROOT_FILE}")
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
        signal, split_background = get_distributions(variable, augmentations)

        background = np.array([])
        background = np.append(background, split_background[0])
        background = np.append(background, split_background[1])
        background = np.append(background, split_background[2])
        background = np.append(background, split_background[3])
        background = np.append(background, split_background[4])
        background = np.append(background, split_background[5])
        background = np.append(background, split_background[6])
        background = np.append(background, split_background[7])
        background = np.append(background, split_background[8])
        background = np.array(background)

        if (len(signal)==0):
            print(" >> signal", len(signal), signal)
            print(" >> empty entries, skipping")
            continue
        

        # print true distribution
        print(f" >> distributed reaches from {np.round(min(min(signal), min(background)),3)} to {np.round(max(max(signal),max(background)),3)}")
        

        # show separation plot
        plot_separation(signal, split_background, target_tuple)
        print(f" >> plots produced")
        print()


def augment_variables():
    # open root file
    with uproot.open(INJECTED_ROOT_FILE) as root_file:     
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
            
            # combine jets to resonance particle if possible
            if not (t1_q1==t1_q2 or t1_q1==t1_b or t1_q2==t1_b or len(jet_energies)<t1_q1+1 or len(jet_energies)<t1_q2+1 or len(jet_energies)<t1_b+1):
                t1_energies.append( jet_energies[t1_q1] + jet_energies[t1_q2] + jet_energies[t1_b] )
            if not (t2_q1==t2_q2 or t2_q1==t2_b or t2_q2==t2_b or len(jet_energies)<t2_q1+1 or len(jet_energies)<t2_q2+1 or len(jet_energies)<t2_b+1):
                t2_energies.append( jet_energies[t2_q1] + jet_energies[t2_q2] + jet_energies[t2_b] )            
            if not (HW_q1==HW_q2 or len(jet_energies)<HW_q1+1 or len(jet_energies)<HW_q2+1):
                HW_energies.append( jet_energies[HW_q1] + jet_energies[HW_q2] )


        met_value = np.array( root_file["neutrino_weighting/reco_met_value"].array() )
        lepton_e = np.array( root_file["neutrino_weighting/reco_lepton_pt"].array() )
        lepton_pt = np.array( root_file["neutrino_weighting/reco_lepton_e"].array() )

    # return augmentations and close root file
    augmentation = {
        "t1_energy" : np.array(t1_energies)/1e3, #scale to MeV 
        "t2_energy" : np.array(t2_energies)/1e3, #scale to MeV
        "HW_energy" : np.array(HW_energies)/1e3, #scale to MeV
        "short_NW_weight" : np.array(short_NW_weight),
        "short_NW_nu_eta" : np.array(short_NW_nu_eta),
        "short_NW_nu_phi" : np.array(short_NW_nu_phi),
        "short_NW_wlep_mass" : np.array(short_NW_wlep_mass)/1e3, #scale to MeV 
        "reco_met_value_gev" : met_value/1e3,
        "reco_lepton_pt_gev" : lepton_e/1e3,
        "reco_lepton_e_gev" : lepton_pt/1e3
    }

    return augmentation


def get_distributions(variable, augmentations):
    signal = []
    background_ttbar = []
    background_ttHbb = []
    background_ttHcc = []
    background_ttHtautau = []
    background_ttHZZ = []
    background_ttHyy = []
    background_ttHWWqqqq = []
    background_ttHWWlvlv = []
    background_other = []
    
    
    with uproot.open(INJECTED_ROOT_FILE) as root_file:			
        event_channels = np.array( root_file[CLASSIFIER_VARIABLE_NAME].array() )
        
        # get distribution from augments or rootfile directly 
        distribution = np.array([])
        
        if variable in augmentations:
            distribution = augmentations[variable]
        else:
            distribution = np.array( root_file[variable].array() )

        # filter by event type
        for (event_channel, value) in zip(event_channels, distribution):
            if event_channel==11:
                signal.append(value)
            else:
                if event_channel==1:
                    background_ttbar.append(value)
                elif event_channel==2:
                    background_ttHbb.append(value)
                elif event_channel==3:
                    background_ttHcc.append(value)
                elif event_channel==4:
                    background_ttHtautau.append(value)
                elif event_channel==5:
                    background_ttHZZ.append(value)
                elif event_channel==6:
                    background_ttHyy.append(value)
                elif event_channel==10:
                    background_ttHWWqqqq.append(value)
                elif event_channel==12:
                    background_ttHWWlvlv.append(value)
                else:
                    background_other.append(value)

    return signal, [background_ttbar, background_ttHbb, background_ttHcc, background_ttHtautau, background_ttHZZ, background_ttHyy, background_ttHWWqqqq, background_ttHWWlvlv, background_other]



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


def plot_separation(signal, split_background, target_tuple):
    # access backgrounds
    ttbar_background = split_background[0]
    ttHbb_background = split_background[1]
    ttHcc_background = split_background[2]
    ttHtautau_background = split_background[3]
    #ttHZZ_background = split_background[4]
    #ttHyy_background = split_background[5]
    ttHWWqqqq_background = split_background[6]
    ttHWWlvlv_background = split_background[7]
    other_background = split_background[8] + split_background[4] + split_background[5]
    
    ttHother = split_background[2] + split_background[3]
    ttHWW = split_background[6] + split_background[7]

    # access target tuple
    variable = target_tuple[0]
    bins = target_tuple[1]
    x_label = target_tuple[2]
    is_int_value = target_tuple[3]

    
    # get correct values
    min_value = np.min([np.min(bins), np.min(bins)])
    max_value = np.max([np.max(bins), np.max(bins)])
    
    # needed weights to scale histogram to 1
    weights_signal = np.ones_like(signal)/float(len(signal))
    weights_ttbar_background = np.ones_like(ttbar_background)/float(len(ttbar_background))
    weights_ttHbb_background = np.ones_like(ttHbb_background)/float(len(ttHbb_background))
    weights_ttHcc_background = np.ones_like(ttHcc_background)/float(len(ttHcc_background))
    weights_ttHtautau_background = np.ones_like(ttHtautau_background)/float(len(ttHtautau_background))
    #weights_ttHZZ_background = np.ones_like(ttHZZ_background)/float(len(ttHZZ_background))
    #weights_ttHyy_background = np.ones_like(ttHyy_background)/float(len(ttHyy_background))
    weights_ttHWWqqqq_background = np.ones_like(ttHWWqqqq_background)/float(len(ttHWWqqqq_background))
    weights_ttHWWlvlv_background = np.ones_like(ttHWWlvlv_background)/float(len(ttHWWlvlv_background))
    weights_other_background = np.ones_like(other_background)/float(len(other_background))

    weights_ttHother_background = np.ones_like(ttHother)/float(len(ttHother))
    weights_ttHWW_background = np.ones_like(ttHWW)/float(len(ttHWW))


    # int value centered bins
    if is_int_value:
        bins = bins - 0.5


    # create plot
    plt.figure(figsize=(8,6))


    # create histogram
    plt.hist(
        [# splt samples
            other_background,
            #ttHyy_background,
            #ttHZZ_background,
            ttbar_background,
            ttHtautau_background,
            ttHcc_background,
            ttHbb_background,
            ttHWWlvlv_background,
            ttHWWqqqq_background,
            signal, 
        ],
        bins, # same bins for everyone
        weights=[# normalize
            weights_other_background,
            #weights_ttHyy_background,
            #weights_ttHZZ_background,
            weights_ttbar_background,
            weights_ttHtautau_background,
            weights_ttHcc_background,
            weights_ttHbb_background,
            weights_ttHWWlvlv_background,
            weights_ttHWWqqqq_background,
            weights_signal, 
        ],
        label=[ # legend
            r"other",
            #r"$t\bar{t}(H\rightarrow \gamma\gamma)$",
            #r"$t\bar{t}(H\rightarrow ZZ)$",
            r"$t\bar{t}$",
            r"$t\bar{t}(H\rightarrow\tau\tau)$",
            r"$t\bar{t}(H\rightarrow cc)$",
            r"$t\bar{t}(H\rightarrow bb)$",
            r"$t\bar{t}(H\rightarrow WW)_\text{dilep}$",
            r"$t\bar{t}(H\rightarrow WW)_\text{had}$",
            r"$t\bar{t}(H\rightarrow WW)_\text{semilep}$",
        ],
        histtype='step',
        stacked=False,
        fill=False
    )
    
    # beautificate
    plt.legend(loc="upper right")
    plt.xlabel(x_label)
    plt.ylabel("Relative Event Yield")
    plt.xlim([min_value,max_value])

    # nice int ticks 
    if is_int_value:
        plt.xlim([min_value-0.5,max_value-0.5])
        plt.xticks(range(int(min_value), int(max_value)))

    # save figure
    plt.savefig(PLOT_DESTINATION+f"separation_split_{variable.split('/')[-1]}.png")
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()



    
    # create plot
    plt.figure(figsize=(8,6))
    

    # create histogram
    plt.hist(
        [# splt samples
            other_background,
            #ttHyy_background,
            #ttHZZ_background,
            ttbar_background,
            ttHbb_background,
            ttHother,
            ttHWW,
            signal, 
        ],
        bins, # same bins for everyone
        weights=[# normalize
            weights_other_background,
            #weights_ttHyy_background,
            #weights_ttHZZ_background,
            weights_ttbar_background,
            weights_ttHbb_background,
            weights_ttHother_background,
            weights_ttHWW_background,
            weights_signal, 
        ],
        label=[ # legend
            r"other",
            #r"$t\bar{t}(H\rightarrow \gamma\gamma)$",
            #r"$t\bar{t}(H\rightarrow ZZ)$",
            r"$t\bar{t}$",
            r"$t\bar{t}(H\rightarrow bb)$",
            r"$t\bar{t}(H\rightarrow\text{other})$",
            r"$t\bar{t}(H\rightarrow WW)$",
            r"$t\bar{t}(H\rightarrow WW)_\text{semilep}$",
        ],
        histtype='step',
        stacked=False,
        fill=False
    )
    
    # beautificate
    plt.legend(loc="upper right")
    plt.xlabel(x_label)
    plt.ylabel("Relative Event Yield")
    plt.xlim([min_value,max_value])

    # nice int ticks 
    if is_int_value:
        plt.xlim([min_value-0.5, max_value-0.5])
        plt.xticks(range(int(min_value), int(max_value)))

    # save figure
    plt.savefig(PLOT_DESTINATION+f"separation_split2_{variable.split('/')[-1]}.png")
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()



if __name__=="__main__":
    main()