import os
import uproot
import numpy as np
import matplotlib.pyplot as plt


# Constants
INPUT_ROOT_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_20-1436656e_weighted.root"
OUTPUT_PATH = "/home/ireas/git_repos/master/txts/"


EVENT_SELECTIONS = [
    (
        "Event Selection v1", 
        [
            "Signal Region: njets>=8, bjets==2, t_probability>0.85, HW_probability>0.85", 
            "Control Region ttbar: njets<8, bjets==2, t_probability<0.85, HW_probability<0.85",
            "Control Region Hbb: bjets>=3",
        ], 
        lambda root_file: my_regions(root_file)
    ),
]

REGION_NAMES = {
    2 : "Signal Region",
    3 : "Control Region ttbar",
    4 : "Control Region tt(H->bb)",
    5 : "Rest Region (should be empty)",
}



def main():
    with uproot.open(INPUT_ROOT_FILE) as root_file:
        # loop troughe each selection
        for (title, description, function) in EVENT_SELECTIONS:
            print(f"> {title}")
            txt_file = os.path.join(OUTPUT_PATH, f"{title}.txt")
            results = function(root_file)
            WriteToFile(txt_file, results, description)
            exit()


def WriteToFile(txt_file, results, Description):
    with open(txt_file, 'w') as f:
        for desc in Description:
            f.write(desc + "\n")

        f.write(f"\n")
        f.write(f"Number of Events: {np.round(20*results[0], 2)} \n")
        f.write(f"Number of Signal Events: {np.round(20*results[1], 2)} ({np.round(20*results[1]/results[0 ], 4)})\n\n\n")

        
        for (i,region) in enumerate(results):
            # skip number of events as no dictionary
            if(i<2):
                continue

            # region
            f.write(f"{REGION_NAMES[i]}\n")
            f.write(f"  (0)  total:           {np.round(20*region.get(0,0), 2)} (total yield: {np.round(100*region.get(0,0)/results[0], 2)}%)\n") 
            f.write(f"  (11) tt(H->WW->qqlv): {np.round(20*region.get(11,0), 2)} ({np.round(region.get(11,0)/region.get(0,1), 4)})  ->  {np.round(100*region.get(11,0)/results[1], 2)}% of all signal | {np.round(100*region.get(11,0)/region.get(0,1), 2)}% of that region\n") 
            f.write(f"  (10) tt(H->WW->qqqq): {np.round(20*region.get(10,0), 2)} ({np.round(region.get(10,0)/region.get(0,1), 4)})\n")
            f.write(f"  (12) tt(H->WW->lvlv): {np.round(20*region.get(12,0), 2)} ({np.round(region.get(12,0)/region.get(0,1), 4)})\n")
            f.write(f"  (2)  tt(H->bb):       {np.round(20*region.get(2,0), 2)} ({np.round(region.get(2,0)/region.get(0,1), 4)})\n")
            f.write(f"  (3)  tt(H->cc):       {np.round(20*region.get(3,0), 2)} ({np.round(region.get(3,0)/region.get(0,1), 4)})\n")
            f.write(f"  (4)  tt(H->tautau):   {np.round(20*region.get(4,0), 2)} ({np.round(region.get(4,0)/region.get(0,1), 4)})\n")
            f.write(f"  (1)  ttbar:           {np.round(20*region.get(1,0), 2)} ({np.round(region.get(1,0)/region.get(0,1), 4)})\n")
            f.write(f"  (9)  other:           {np.round(20*region.get(9,0), 2)} ({np.round(region.get(9,0)/region.get(0,1), 4)})\n")
            f.write(f"\n")


def my_regions(root_file):
    # prepare regions
    SR, CR_ttbar, CR_Hbb, RR = {}, {}, {}, {}
    
    # access variables from the root file
    njets = np.array(root_file["neutrino_weighting/number_of_jets"].array())
    nbjets = np.array(root_file["neutrino_weighting/number_of_bjets"].array())
    event_channels = np.array(root_file["neutrino_weighting/classification_event_channel"].array())
    sm_weights = np.array( root_file["neutrino_weighting/SM_event_weight"].array() )
    missing_energies = np.array( root_file["neutrino_weighting/reco_met_value"].array() )/1e3

    jet_energies_vecs = root_file["neutrino_weighting/jet_e_NOSYS"].array()
    prediction_t1_q1 = np.array( root_file["neutrino_weighting/spanet_t1_q1"].array(), dtype=np.intc )
    prediction_t1_q2 = np.array( root_file["neutrino_weighting/spanet_t1_q2"].array(), dtype=np.intc )
    prediction_t1_b = np.array( root_file["neutrino_weighting/spanet_t1_b"].array(), dtype=np.intc )
    prediction_t2_q1 = np.array( root_file["neutrino_weighting/spanet_t2_q1"].array(), dtype=np.intc )
    prediction_t2_q2 = np.array( root_file["neutrino_weighting/spanet_t2_q2"].array(), dtype=np.intc )
    prediction_t2_b = np.array( root_file["neutrino_weighting/spanet_t2_b"].array(), dtype=np.intc )
    prediction_HW_q1 = np.array( root_file["neutrino_weighting/spanet_HW_q1"].array(), dtype=np.intc )
    prediction_HW_q2 = np.array( root_file["neutrino_weighting/spanet_HW_q2"].array(), dtype=np.intc )

    
    t1_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_t1_assignment_probability"].array() )
    t2_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_t2_assignment_probability"].array() )
    HW_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_HW_assignment_probability"].array() )

    number_of_events = 0
    number_of_signal_events = 0

    # loop through the variables
    for (
        event_channel, sm_weight, njet, nbjet, met,
        t1_q1, t1_q2, t1_b, t2_q1, t2_q2, t2_b, HW_q1, HW_q2,
        t1_assignment_probability, t2_assignment_probability, HW_assignment_probability 
    ) in zip(
        event_channels, sm_weights, njets, nbjets, missing_energies,
        prediction_t1_q1, prediction_t1_q2, prediction_t1_b, prediction_t2_q1, prediction_t2_q2, prediction_t2_b, prediction_HW_q1, prediction_HW_q2,
        t1_assignment_probabilities, t2_assignment_probabilities, HW_assignment_probabilities 
    ):

        number_of_events+= sm_weight
        if(event_channel==11):
            number_of_signal_events+= sm_weight

        # SR
        if njet>=8 and nbjet==2 and t1_assignment_probability>=0.85 and HW_assignment_probability>=0.85:
            SR[0] = SR.get(0, 0) + sm_weight
            SR[event_channel] = SR.get(event_channel, 0) + sm_weight
        # CR ttbar
        elif (njet<8 and nbjet==2) or (njet>=8 and nbjet==2 and (t1_assignment_probability<0.85 or HW_assignment_probability<0.85) ):
            CR_ttbar[0] = CR_ttbar.get(0, 0) + sm_weight
            CR_ttbar[event_channel] = CR_ttbar.get(event_channel, 0) + sm_weight
        # CR Hbb
        elif nbjet>=3:
            CR_Hbb[0] = CR_Hbb.get(0, 0) + sm_weight
            CR_Hbb[event_channel] = CR_Hbb.get(event_channel, 0) + sm_weight
        # Rest, should be empty
        else:
            print(njet, nbjet)
            RR[0] = RR.get(0, 0) + sm_weight
            RR[event_channel] = RR.get(event_channel, 0) + sm_weight
    
    return number_of_events, number_of_signal_events, SR, CR_ttbar, CR_Hbb, RR



def preselection(root_file):
    # limits
    min_jet_limit = 7
    exact_bjet_limit = 3
    upper_met_limit = 50

    # Access variables from the root file
    event_channels = np.array(root_file["neutrino_weighting/classification_event_channel"].array())
    sm_weights = np.array( root_file["neutrino_weighting/SM_event_weight"].array() )

    njets = np.array(root_file["neutrino_weighting/number_of_jets"].array())
    nbjets = np.array(root_file["neutrino_weighting/number_of_bjets"].array())
    missing_energies = np.array( root_file["neutrino_weighting/reco_met_value"].array() )/1e3

    spanet_t1_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_t1_assignment_probability"].array() )
    spanet_t2_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_t2_assignment_probability"].array() )
    spanet_HW_assignment_probabilities = np.array( root_file["neutrino_weighting/spanet_HW_assignment_probability"].array() )
    spanet_t1_detection_probabilities = np.array( root_file["neutrino_weighting/spanet_t1_detection_probability"].array() )
    spanet_t2_detection_probabilities = np.array( root_file["neutrino_weighting/spanet_t2_detection_probability"].array() )
    spanet_HW_detection_probabilities = np.array( root_file["neutrino_weighting/spanet_HW_detection_probability"].array() )

    nw_weights = np.array( root_file["neutrino_weighting/NW_weight"].array() )


    # define dictionaries (0=all and event_channel)
    unselected = {}
    selected_jets = {}
    selected_bjets = {}
    selected_nw = {}
    selected_spanet = {}


    # loop through preselection and collect stuff
    for (
        event_channel, sm_weight, njet, nbjet, met, nw_weight, 
        t1_assignment, t2_assignment, HW_assignment,
        t1_detection, t2_detection, HW_detection
    ) in zip(
            event_channels, sm_weights, njets, nbjets, missing_energies, nw_weights, 
            spanet_t1_assignment_probabilities, spanet_t2_assignment_probabilities, spanet_HW_assignment_probabilities,
            spanet_t1_detection_probabilities, spanet_t2_detection_probabilities, spanet_HW_detection_probabilities
    ):
        
        # all events
        unselected[0] = unselected.get(0, 0) + sm_weight
        unselected[event_channel] = unselected.get(event_channel, 0) + sm_weight


        # enough jet events
        if not (njet>min_jet_limit):
            continue
        selected_jets[0] = selected_jets.get(0, 0) + sm_weight
        selected_jets[event_channel] = selected_jets.get(event_channel, 0) + sm_weight


        # enough bjet events
        if not (nbjet<=exact_bjet_limit):
            continue
        selected_bjets[0] = selected_bjets.get(0, 0) + sm_weight
        selected_bjets[event_channel] = selected_bjets.get(event_channel, 0) + sm_weight


        # nw
        if not (nw_weight>=0.0):
            continue
        selected_nw[0] = selected_nw.get(0, 0) + sm_weight
        selected_nw[event_channel] = selected_nw.get(event_channel, 0) + sm_weight


        # spanet
        if not (HW_assignment>0.85 and t1_assignment>0.85):
            continue
        selected_spanet[0] = selected_spanet.get(0, 0) + sm_weight
        selected_spanet[event_channel] = selected_spanet.get(event_channel, 0) + sm_weight


    # print results
    print(f"No Selection")
    print(f"  > total:         {np.round(20*unselected.get(0,0), 2)} events ({ np.round(100* unselected.get(0,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttbar:         {np.round(20*unselected.get(1,0), 2)} events ({ np.round(100* unselected.get(1,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->bb:       {np.round(20*unselected.get(2,0), 2)} events ({ np.round(100* unselected.get(2,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->cc:       {np.round(20*unselected.get(3,0), 2)} events ({ np.round(100* unselected.get(3,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->tautau:   {np.round(20*unselected.get(4,0), 2)} events ({ np.round(100* unselected.get(4,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->ZZ:       {np.round(20*unselected.get(5,0), 2)} events ({ np.round(100* unselected.get(5,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->yy:       {np.round(20*unselected.get(6,0), 2)} events ({ np.round(100* unselected.get(6,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->other:    {np.round(20*unselected.get(9,0), 2)} events ({ np.round(100* unselected.get(9,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqqq: {np.round(20*unselected.get(10,0), 2)} events ({ np.round(100* unselected.get(10,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->WW->lvlv: {np.round(20*unselected.get(12,0), 2)} events ({ np.round(100* unselected.get(12,0)/unselected.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqlv: {np.round(20*unselected.get(11,0), 2)} events ({ np.round(100* unselected.get(11,0)/unselected.get(0,1), 2) }%)")
    print(f"")

    
    print(f"Jet Selection (>{min_jet_limit} jets)")
    print(f"  > total:         {np.round(20*selected_jets.get(0,0), 2)} events ({ np.round(100* selected_jets.get(0,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttbar:         {np.round(20*selected_jets.get(1,0), 2)} events ({ np.round(100* selected_jets.get(1,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->bb:       {np.round(20*selected_jets.get(2,0), 2)} events ({ np.round(100* selected_jets.get(2,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->cc:       {np.round(20*selected_jets.get(3,0), 2)} events ({ np.round(100* selected_jets.get(3,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->tautau:   {np.round(20*selected_jets.get(4,0), 2)} events ({ np.round(100* selected_jets.get(4,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->ZZ:       {np.round(20*selected_jets.get(5,0), 2)} events ({ np.round(100* selected_jets.get(5,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->yy:       {np.round(20*selected_jets.get(6,0), 2)} events ({ np.round(100* selected_jets.get(6,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->other:    {np.round(20*selected_jets.get(9,0), 2)} events ({ np.round(100* selected_jets.get(9,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqqq: {np.round(20*selected_jets.get(10,0), 2)} events ({ np.round(100* selected_jets.get(10,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->lvlv: {np.round(20*selected_jets.get(12,0), 2)} events ({ np.round(100* selected_jets.get(12,0)/selected_jets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqlv: {np.round(20*selected_jets.get(11,0), 2)} events ({ np.round(100* selected_jets.get(11,0)/selected_jets.get(0,1), 2) }%)")
    print(f"")
    

    print(f"b-Jet Selection (={exact_bjet_limit} bjets)")
    print(f"  > total:         {np.round(20*selected_bjets.get(0,0), 2)} events ({ np.round(100* selected_bjets.get(0,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttbar:         {np.round(20*selected_bjets.get(1,0), 2)} events ({ np.round(100* selected_bjets.get(1,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->bb:       {np.round(20*selected_bjets.get(2,0), 2)} events ({ np.round(100* selected_bjets.get(2,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->cc:       {np.round(20*selected_bjets.get(3,0), 2)} events ({ np.round(100* selected_bjets.get(3,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->tautau:   {np.round(20*selected_bjets.get(4,0), 2)} events ({ np.round(100* selected_bjets.get(4,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->ZZ:       {np.round(20*selected_bjets.get(5,0), 2)} events ({ np.round(100* selected_bjets.get(5,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->yy:       {np.round(20*selected_bjets.get(6,0), 2)} events ({ np.round(100* selected_bjets.get(6,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->other:    {np.round(20*selected_bjets.get(9,0), 2)} events ({ np.round(100* selected_bjets.get(9,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqqq: {np.round(20*selected_bjets.get(10,0), 2)} events ({ np.round(100* selected_bjets.get(10,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->lvlv: {np.round(20*selected_bjets.get(12,0), 2)} events ({ np.round(100* selected_bjets.get(12,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqlv: {np.round(20*selected_bjets.get(11,0), 2)} events ({ np.round(100* selected_bjets.get(11,0)/selected_bjets.get(0,1), 2) }%)")
    print(f"")


    print(f"NW Selection")
    print(f"  > total:         {np.round(20*selected_nw.get(0,0), 2)} events ({ np.round(100* selected_nw.get(0,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttbar:         {np.round(20*selected_nw.get(1,0), 2)} events ({ np.round(100* selected_nw.get(1,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->bb:       {np.round(20*selected_nw.get(2,0), 2)} events ({ np.round(100* selected_nw.get(2,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->cc:       {np.round(20*selected_nw.get(3,0), 2)} events ({ np.round(100* selected_nw.get(3,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->tautau:   {np.round(20*selected_nw.get(4,0), 2)} events ({ np.round(100* selected_nw.get(4,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->ZZ:       {np.round(20*selected_nw.get(5,0), 2)} events ({ np.round(100* selected_nw.get(5,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->yy:       {np.round(20*selected_nw.get(6,0), 2)} events ({ np.round(100* selected_nw.get(6,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->other:    {np.round(20*selected_nw.get(9,0), 2)} events ({ np.round(100* selected_nw.get(9,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqqq: {np.round(20*selected_nw.get(10,0), 2)} events ({ np.round(100* selected_nw.get(10,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->WW->lvlv: {np.round(20*selected_nw.get(12,0), 2)} events ({ np.round(100* selected_nw.get(12,0)/selected_nw.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqlv: {np.round(20*selected_nw.get(11,0), 2)} events ({ np.round(100* selected_nw.get(11,0)/selected_nw.get(0,1), 2) }%)")
    print(f"")


    print(f"SPANet Selection")
    print(f"  > total:         {np.round(20*selected_spanet.get(0,0), 2)} events ({ np.round(100* selected_spanet.get(0,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttbar:         {np.round(20*selected_spanet.get(1,0), 2)} events ({ np.round(100* selected_spanet.get(1,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->bb:       {np.round(20*selected_spanet.get(2,0), 2)} events ({ np.round(100* selected_spanet.get(2,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->cc:       {np.round(20*selected_spanet.get(3,0), 2)} events ({ np.round(100* selected_spanet.get(3,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->tautau:   {np.round(20*selected_spanet.get(4,0), 2)} events ({ np.round(100* selected_spanet.get(4,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->ZZ:       {np.round(20*selected_spanet.get(5,0), 2)} events ({ np.round(100* selected_spanet.get(5,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->yy:       {np.round(20*selected_spanet.get(6,0), 2)} events ({ np.round(100* selected_spanet.get(6,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->other:    {np.round(20*selected_spanet.get(9,0), 2)} events ({ np.round(100* selected_spanet.get(9,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqqq: {np.round(20*selected_spanet.get(10,0), 2)} events ({ np.round(100* selected_spanet.get(10,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->WW->lvlv: {np.round(20*selected_spanet.get(12,0), 2)} events ({ np.round(100* selected_spanet.get(12,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"  > ttH->WW->qqlv: {np.round(20*selected_spanet.get(11,0), 2)} events ({ np.round(100* selected_spanet.get(11,0)/selected_spanet.get(0,1), 2) }%)")
    print(f"")

    return 0



def plot_event_distribution(data_list, output_filename):
    # Unpack the input list
    total_events = data_list[0]       # Total number of events (not used for plot)
    signal_events = data_list[1]      # Number of signal events (not used for plot)
    regions = data_list[2:]           # List of SR and CRs (dictionaries)

    # Get unique event channels
    event_channels = set()
    for region in regions:
        event_channels.update(region.keys())
    event_channels = list(event_channels)  # Convert to list for consistent ordering

    # Create an array to store the event counts for each region and channel
    event_counts = []
    for region in regions:
        counts = [region.get(channel, 0) for channel in event_channels]  # Get count for each channel, default to 0 if not present
        event_counts.append(counts)

    # Convert event_counts to a NumPy array for easier manipulation
    event_counts = np.array(event_counts)

    # Generate the plot
    plt.figure(figsize=(10, 6))
    
    # Create stacked bar plot
    bottom = np.zeros(len(regions))  # Initialize bottom for stacking
    for i, channel in enumerate(event_channels):
        plt.bar(range(len(regions)), event_counts[:, i], bottom=bottom, label=channel)
        bottom += event_counts[:, i]  # Update bottom for next stack

    # Add labels and title
    plt.xlabel('Regions (SR, CR_1, CR_2, etc.)')
    plt.ylabel('Number of Events')
    plt.yscale('log')  # Set y-scale to logarithmic
    plt.xticks(range(len(regions)), [f'Region {i+1}' for i in range(len(regions))], rotation=45)
    plt.legend(title='Event Channels')

    # Save the figure
    plt.tight_layout()
    plt.savefig(output_filename)



if __name__ == "__main__":
    main()