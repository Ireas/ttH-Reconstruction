import uproot
import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS
#INJECTED_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_8+j_1530766e_weighted.root"
#TREE = "neutrino_weighting"

INJECTED_ROOT_FILE = "/media/ireas/Data/v6/injected/all_8+j_1530766e_injected.root"

def main():
    # define region
    signal_region = []
    control_region_Hbb = []
    rest_region = []

    n_events = 0
    n_events_signal = 0


    # fill regions
    with uproot.open(INJECTED_ROOT_FILE) as root_file:			
        # classifications
        classifications_event_channel = np.array( root_file["matched/classification_event_channel"].array() )
    
        # spanet prediction
        probabilities_t1_marginal = np.array( root_file["spanet/t1_marginal_probability"].array() )
        probabilities_t2_marginal = np.array( root_file["spanet/t2_marginal_probability"].array() )
        probabilities_HW_marginal = np.array( root_file["spanet/HW_marginal_probability"].array() )

        # jet multiplicities
        b_jet_multiplicities = np.array( root_file["matched/number_of_bjets"].array() )

        # kinematics
        met_values_gev = np.array( root_file["matched/reco_met_value"].array() )/1e3
        
        sm_weights = np.array( root_file["matched/SM_event_weight"].array() )
    

        # loop and count
        for(
            sm_weight,
            channel,
            probability_t1_marginal, probability_t2_marginal, probability_HW_marginal,
            b_jet_multiplicity, met_value_gev
        ) in zip(
            sm_weights,
            classifications_event_channel,
            probabilities_t1_marginal, probabilities_t2_marginal, probabilities_HW_marginal,
            b_jet_multiplicities, met_values_gev,
        ):
            # transform to uniform event_var
            
            # signal region selection
            if (probability_t1_marginal>0.3) and (probability_t2_marginal>0.3) and (probability_HW_marginal>0.3):
                signal_region.append( (channel, sm_weight) )
            # control region for tt(H->bb)
            elif (b_jet_multiplicity>2) and (met_value_gev<20):
                control_region_Hbb.append( (channel, sm_weight) )
            # rest region dump
            else:
                rest_region.append( (channel, sm_weight) )
            
            
            # counters
            n_events+= sm_weight
            if channel==11:
                n_events_signal+= sm_weight


    # calculate sm weighted events
    sr_ttHWWlvqq = 0
    sr_ttHWWqqqq = 0
    sr_ttHWWlvlv = 0
    sr_ttHbb = 0
    sr_ttHtautau = 0
    sr_ttbar = 0
    sr_other = 0
    sr_all = 0
    cr_ttHWWlvqq = 0
    cr_ttHWWqqqq = 0
    cr_ttHWWlvlv = 0
    cr_ttHbb = 0
    cr_ttHtautau = 0
    cr_ttbar = 0
    cr_other = 0
    cr_all = 0
    rr_ttHWWlvqq = 0
    rr_ttHWWqqqq = 0
    rr_ttHWWlvlv = 0
    rr_ttHbb = 0
    rr_ttHtautau = 0
    rr_ttbar = 0
    rr_other = 0
    rr_all = 0

    for (channel, sm_weight) in signal_region:
        sr_all+= sm_weight
        
        if channel==11:
            sr_ttHWWlvqq+= sm_weight
        elif channel==10:
            sr_ttHWWqqqq+= sm_weight
        elif channel==12:
            sr_ttHWWlvlv+= sm_weight
        elif channel==2:
            sr_ttHbb+= sm_weight
        elif channel==4:
            sr_ttHtautau+= sm_weight
        elif channel==1:
            sr_ttbar+= sm_weight
        else:
            sr_other+= sm_weight

    for (channel, sm_weight) in control_region_Hbb:
        cr_all+= sm_weight

        if channel==11:
            cr_ttHWWlvqq+= sm_weight
        elif channel==10:
            cr_ttHWWqqqq+= sm_weight
        elif channel==12:
            cr_ttHWWlvlv+= sm_weight
        elif channel==2:
            cr_ttHbb+= sm_weight
        elif channel==4:
            cr_ttHtautau+= sm_weight
        elif channel==1:
            cr_ttbar+= sm_weight
        else:
            cr_other+= sm_weight

    for (channel, sm_weight) in rest_region:
        rr_all+= sm_weight
        
        if channel==11:
            rr_ttHWWlvqq+= sm_weight
        elif channel==10:
            rr_ttHWWqqqq+= sm_weight
        elif channel==12:
            rr_ttHWWlvlv+= sm_weight
        elif channel==2:
            rr_ttHbb+= sm_weight
        elif channel==4:
            rr_ttHtautau+= sm_weight
        elif channel==1:
            rr_ttbar+= sm_weight
        else:
            rr_other+= sm_weight



    # print regions
    print(f"Signal Region: {np.round(sr_all, 2)} ({np.round( 100*sr_all/n_events ,2)}%) events")    
    print(f" > tt(H->WW->qqlv): {np.round(sr_ttHWWlvqq, 2)} ({np.round( 100*sr_ttHWWlvqq/sr_all, 2)}%) ==> ({np.round( 100*sr_ttHWWlvqq/n_events_signal, 2)}%) of all signal")
    print(f" > tt(H->WW->qqqq): {np.round(sr_ttHWWqqqq, 2)} ({np.round( 100*sr_ttHWWqqqq/sr_all, 2)}%)")
    print(f" > tt(H->WW->lvlv): {np.round(sr_ttHWWlvlv, 2)} ({np.round( 100*sr_ttHWWlvlv/sr_all, 2)}%)")
    print(f" > tt(H->bb):       {np.round(sr_ttHbb, 2)} ({np.round( 100*sr_ttHbb/sr_all, 2)}%)")
    print(f" > tt(H->tautau):   {np.round(sr_ttHtautau, 2)} ({np.round( 100*sr_ttHtautau/sr_all, 2)}%)")
    print(f" > ttbar:           {np.round(sr_ttbar, 2)} ({np.round( 100*sr_ttbar/sr_all, 2)}%)")
    print(f" > other:           {np.round(sr_other, 2)} ({np.round( 100*sr_other/sr_all, 2)}%)")
    print(f"")
    print(f"Control Region tt(H->bb): {np.round(cr_all, 2)} ({np.round( 100*cr_all/n_events ,2)}%) events")    
    print(f" > tt(H->WW->qqlv): {np.round(cr_ttHWWlvqq, 2)} ({np.round( 100*cr_ttHWWlvqq/cr_all, 2)}%) ==> ({np.round( 100*cr_ttHWWlvqq/n_events_signal, 2)}%) of all signal")
    print(f" > tt(H->WW->qqqq): {np.round(cr_ttHWWqqqq, 2)} ({np.round( 100*cr_ttHWWqqqq/cr_all, 2)}%)")
    print(f" > tt(H->WW->lvlv): {np.round(cr_ttHWWlvlv, 2)} ({np.round( 100*cr_ttHWWlvlv/cr_all, 2)}%)")
    print(f" > tt(H->bb):       {np.round(cr_ttHbb, 2)} ({np.round( 100*cr_ttHbb/cr_all, 2)}%)")
    print(f" > tt(H->tautau):   {np.round(cr_ttHtautau, 2)} ({np.round( 100*cr_ttHtautau/cr_all, 2)}%)")
    print(f" > ttbar:           {np.round(cr_ttbar, 2)} ({np.round( 100*cr_ttbar/cr_all, 2)}%)")
    print(f" > other:           {np.round(cr_other, 2)} ({np.round( 100*cr_other/cr_all, 2)}%)")
    print(f"")
    print(f"Rest Region: {np.round(rr_all, 2)} ({np.round( 100*rr_all/n_events ,2)}%) events")    
    print(f" > tt(H->WW->qqlv): {np.round(rr_ttHWWlvqq, 2)} ({np.round( 100*rr_ttHWWlvqq/rr_all, 2)}%) ==> ({np.round( 100*rr_ttHWWlvqq/n_events_signal, 2)}%) of all signal")
    print(f" > tt(H->WW->qqqq): {np.round(rr_ttHWWqqqq, 2)} ({np.round( 100*rr_ttHWWqqqq/rr_all, 2)}%)")
    print(f" > tt(H->WW->lvlv): {np.round(rr_ttHWWlvlv, 2)} ({np.round( 100*rr_ttHWWlvlv/rr_all, 2)}%)")
    print(f" > tt(H->bb):       {np.round(rr_ttHbb, 2)} ({np.round( 100*rr_ttHbb/rr_all, 2)}%)")
    print(f" > tt(H->tautau):   {np.round(rr_ttHtautau, 2)} ({np.round( 100*rr_ttHtautau/rr_all, 2)}%)")
    print(f" > ttbar:           {np.round(rr_ttbar, 2)} ({np.round( 100*rr_ttbar/rr_all, 2)}%)")
    print(f" > other:           {np.round(rr_other, 2)} ({np.round( 100*rr_other/rr_all, 2)}%)")



def transform_to_event_var(event_completion, true_t1_decay, true_t2_decay, true_higgs_decay):
    # tt(H->WW->qqlv)
    if event_completion==1 or event_completion==-2:
        return 1

    # tt(H->WW->qqqq)
    if true_higgs_decay==10:
        return 2

    # tt(H->bb)
    if true_higgs_decay==2:
        return 3

    # other    
    return 0




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