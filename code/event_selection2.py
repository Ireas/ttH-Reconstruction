import os
import uproot
import numpy as np
import matplotlib.pyplot as plt


# Options
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


# Constants
INPUT_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_50_ttHWW_amplified_1675331e_weighted_on_reco_without_spanet.root"
PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/regions/"

REGION_LABELS = [
    r"SR",
    r"CR$_{t\bar{t}}$",
    r"CR$_{H->bb}$",
]

CHANNEL_LABELS = {
    0  : r"undefined",
    1  :r"$t\bar{t}$ (scaled 1:100)",
    2  :r"$t\bar{t}(H\rightarrow bb)$",
    3  :r"$t\bar{t}(H\rightarrow cc)$",
    4  :r"$t\bar{t}(H\rightarrow\tau\tau)$",
    5  :r"$t\bar{t}(H\rightarrow\gamma\gamma)$",
    6  :r"$t\bar{t}(H\rightarrow ZZ)$",
    7  :r"undefined",
    8  :r"undefined",
    9  :r"$t\bar{t}(H\rightarrow other)$",
    10 :r"$t\bar{t}(H\rightarrow WW\rightarrow qqqq)$",
    11 :r"$t\bar{t}(H\rightarrow WW\rightarrow qqlv)$",
    12 :r"$t\bar{t}(H\rightarrow WW\rightarrow lvlv)$",
    21 :r"$t\bar{t}(H\rightarrow\text{other})$",
    22 :r"$t\bar{t}(H\rightarrow WW)$",
}

# Functions
def main():
    with uproot.open(INPUT_ROOT_FILE) as root_file:
        TR, SR, CR_ttbar, CR_Hbb = create_regions(root_file)

        # calculate metrics
        SR_background_rejection = calculate_background_rejection(TR, SR)
        SR_signal_acceptance = calculate_channel_acceptance(TR, SR, [11])
        SR_signal_purity, SR_signal_yield = calculate_channel_purity(SR, [11])
        
        CR_ttbar_acceptance = calculate_channel_acceptance(TR, CR_ttbar, [1])
        CR_ttbar_signal_contamination, CR_ttbar_signal_yield = calculate_channel_purity(CR_ttbar, [11])

        CR_Hbb_acceptance = calculate_channel_acceptance(TR, CR_Hbb, [2])
        CR_Hbb_signal_contamination, CR_Hbb_signal_yield = calculate_channel_purity(CR_Hbb, [11])

        # print values
        print(f" SR background rejection: {np.round(SR_background_rejection, 3)}")
        print(f" SR signal acceptance: {np.round(SR_signal_acceptance, 3)}")
        print(f" SR signal purity: {np.round(SR_signal_purity, 3)} @ {np.round(SR_signal_yield,1)} events")
        print(f"")
        print(f" CR ttbar acceptance: {np.round(CR_ttbar_acceptance, 3)}")
        print(f" CR ttbar signal contamination: {np.round(CR_ttbar_signal_contamination, 3)} @ {np.round(CR_ttbar_signal_yield,1)} events")
        print(f"")
        print(f" CR Hbb acceptance: {np.round(CR_Hbb_acceptance, 3)}")
        print(f" CR Hbb signal contamination: {np.round(CR_Hbb_signal_contamination, 3)}  @ {np.round(CR_Hbb_signal_yield,1)} events")

        # plot regions
        plot_regions(TR, SR, CR_ttbar, CR_Hbb)


def create_regions(root_file):
    # create containers
    TR, SR, CR_ttbar, CR_Hbb, RR = {},{},{},{},{}

    # access variables
    sm_weights = np.array( root_file["neutrino_weighting/SM_event_weight"].array() )
    event_channels = np.array(root_file["neutrino_weighting/classification_event_channel"].array())

    jet_multiplicities = np.array(root_file["neutrino_weighting/number_of_jets"].array())
    bjet_multiplicities = np.array(root_file["neutrino_weighting/number_of_bjets"].array())



    # loop trough all events and fill regions
    for(
        sm_weight, event_channel,
        njets, nbjets
    ) in zip(
        sm_weights, event_channels,
        jet_multiplicities, bjet_multiplicities
    ):
        # TR (all events)
        TR[event_channel] = TR.get(event_channel, 0) + sm_weight


        # SR
        if(njets>=6 and nbjets==2):
            SR[event_channel] = SR.get(event_channel, 0) + sm_weight
        # CR_ttbar
        elif(njets<6 and nbjets==2):
            CR_ttbar[event_channel] = CR_ttbar.get(event_channel, 0) + sm_weight
        # CR_Hbb
        elif(nbjets>2):
            CR_Hbb[event_channel] = CR_Hbb.get(event_channel, 0) + sm_weight
        # RR (should never reach this)
        else:
            RR[event_channel] = RR.get(event_channel, 0) + sm_weight



    # check if rest region has any values
    if( len(RR.items())!=0 ):
        print("Warning: RR is not empty!")
        for (key,value) in RR.items():
            print(f"  >{key}: {value}")

    return TR, SR, CR_ttbar, CR_Hbb


def calculate_background_rejection(TR, SR):
    # estimate total background
    total_number_of_background_events = 0
    for (key, value) in TR.items():
        if(key!=11):
            total_number_of_background_events+= value
    
    if(total_number_of_background_events==0):
        print("Warning: no background events in TR!")
        return 1
    
    
    # estimate background in SR
    number_of_background_events_in_signal_region = 0
    for (key, value) in SR.items():
        if(key!=11):
            number_of_background_events_in_signal_region+= value
    
    # calculate background rejection
    background_rejection = 1 - number_of_background_events_in_signal_region / total_number_of_background_events

    return background_rejection


def calculate_channel_acceptance(TR, region, accepted_channels):
    total_number_of_events = 0
    for (key, value) in TR.items():
        if(key in accepted_channels):
            total_number_of_events+= value
    
    if(total_number_of_events==0):
        print("Warning: no channel events in TR!")
        return 1
    
    
    # estimate background in SR
    number_of_events_in_region = 0
    for (key, value) in region.items():
        if(key in accepted_channels):
            number_of_events_in_region+= value

    # calculate acceptance
    acceptance = number_of_events_in_region / total_number_of_events

    return acceptance



def calculate_channel_purity(region, accepted_channels):
    number_of_pure_events = 0
    number_of_impure_events = 0

    # count events
    for (key, value) in region.items():
        if(key in accepted_channels):
            number_of_pure_events+= value
        else:
            number_of_impure_events+= value

    if(number_of_pure_events+number_of_impure_events==0):
        print("Warning: region is empty!")
        return 0

    # calculate purity
    purity = number_of_pure_events /(number_of_pure_events + number_of_impure_events)
    
    return purity, number_of_pure_events



def plot_regions(TR, SR, CR_ttbar, CR_Hbb):
    # prepare data
    b_values = np.zeros(len(REGION_LABELS))
    
    weights = {}
    for channel in TR.keys():
        weights[channel] = np.array( [ SR.get(channel,0), CR_ttbar.get(channel,0), CR_Hbb.get(channel,0) ] )

    # combine events for clarity
    weights[21] = weights[9] + weights[6] + weights[5] + weights[3]
    weights[22] = weights[10] + weights[12]

    # scale ttbar for clarity
    weights[1] = weights[1]/1e2 


    # create a figure
    plt.figure(figsize=(10,8))

    
    for (channel) in [11,22,2,4,21,1]:
        plt.bar(REGION_LABELS, weights[channel], label=CHANNEL_LABELS[channel], bottom = b_values)
        b_values+= weights[channel]


    # beautificate
    #plt.yscale('log')
    plt.legend()

    # save the figure
    plt.savefig(f"{PLOT_DESTINATION}regions.png")
    plt.close()



if __name__ == "__main__":
    main()