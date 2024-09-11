import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


WHAD_MASS_PDG = 80.4335
WHAD_MASS_MAX_DEVIATION = 1 # maximum valid diviation from PDG mass for on-shell definition

PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/validate_nw/"
WEIGHTED_FILE = "/media/ireas/Data/v6/weighted/all_5+j_truncated_20_1436656e_weighted_on_truth_125_2.root"


plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize


TARGETS = [
	(
		r"delta_mw", # x-data name
		r"delta_nueta", # y-data name
		[np.linspace(-0.5,0.5,11), np.linspace(-0.1,0.1,11)], # binning
		r"Difference in Estimation $\Delta m_{W^*}$ [GeV]", # x-label
		r"Difference in Estimation $\Delta \eta_\nu$", # y-label
		True, # normalise
	),
]


def main():
	# count events
	count_events()

	# create dictionary
	root_dictionary = create_dictionary()
	
	# loop trough targets
	for target in TARGETS:
		plot_2d(root_dictionary, target)


def count_events():
	with uproot.open(WEIGHTED_FILE) as root_file:
		total_counts = {}
		success_counts = {}
		good_success_counts = {}

		NW_weights = root_file["neutrino_weighting/NW_weight"].array()
		event_channels = root_file["neutrino_weighting/classification_event_channel"].array()

		for (weight, channel) in zip(NW_weights, event_channels):
			total_counts[channel] = total_counts.get(channel,0) + 1
			if(weight>=0):
				success_counts[channel] = success_counts.get(channel,0) + 1
			if(weight>=0.5):
				good_success_counts[channel] = good_success_counts.get(channel,0) + 1

		print_channel("ttbar             ", 1 , total_counts, success_counts, good_success_counts)
		print_channel("ttHbb             ", 2,  total_counts, success_counts, good_success_counts)
		print_channel("ttHcc             ", 3,  total_counts, success_counts, good_success_counts)
		print_channel("ttHtautau         ", 4,  total_counts, success_counts, good_success_counts)
		print_channel("ttHZZ             ", 5,  total_counts, success_counts, good_success_counts)
		print_channel("ttHyy             ", 6,  total_counts, success_counts, good_success_counts)
		print_channel("other             ", 9 , total_counts, success_counts, good_success_counts)
		print_channel("ttHWW semileptonic", 11, total_counts, success_counts, good_success_counts)
		print_channel("ttHWW dileptonic  ", 12, total_counts, success_counts, good_success_counts)
		print_channel("ttHWW hadronic    ", 10, total_counts, success_counts, good_success_counts)

def print_channel(name, channel, total, success, good):
	print(f"{name}: {success.get(channel,0)} / {total.get(channel,0)} = {np.round(success.get(channel,0)/total.get(channel,1), 2)} | {good.get(channel,0)} / {total.get(channel,0)} = {np.round(good.get(channel,0)/total.get(channel,1), 2)}")



def create_dictionary():	
	# create empty dictionary
	root_dictionary = {}

	with uproot.open(WEIGHTED_FILE) as root_file:
		# fill dictionary
		root_dictionary["NW_H_mass"] = root_file["neutrino_weighting/NW_H_mass"].array()
		NW_weights = root_file["neutrino_weighting/NW_weight"].array()

		NW_wlep_mass = root_file["neutrino_weighting/NW_wlep_mass"].array()
		NW_nu_eta = root_file["neutrino_weighting/NW_nu_eta"].array()
		true_wlep_m = root_file["neutrino_weighting/true_wlep_m"].array()
		true_neutrino_eta = root_file["neutrino_weighting/true_neutrino_eta"].array()
		delta_mw = []
		delta_nueta = []
		
		for (weight, nw_mw, nw_nueta, true_mv, true_nueta) in zip(NW_weights, NW_wlep_mass, NW_nu_eta, true_wlep_m, true_neutrino_eta):
			if(weight<0.05):
				continue

			delta_mw.append(nw_mw-true_mv)
			delta_nueta.append(nw_nueta-true_nueta)
		
		root_dictionary["delta_mw"] = np.array(delta_mw)/1e3 #convert to GeV
		root_dictionary["delta_nueta"] = np.array(delta_nueta)
		
		root_dictionary["NW_weights"] = NW_weights


		#root_dictionary["NW_on_true_weights"] = root_file["neutrino_weighting/NW_on_true_weight"].array()
		#root_dictionary["NW_on_true_H_smear"] = root_file["neutrino_weighting/NW_on_true_H_smear"].array()
		#root_dictionary["NW_on_true_wlep_masses"] = root_file["neutrino_weighting/NW_on_true_wlep_mass"].array()

		#root_dictionary["NW_on_reco_weights"] = root_file["neutrino_weighting/NW_on_reco_weight"].array()
		#root_dictionary["NW_on_reco_H_smear"] = root_file["neutrino_weighting/NW_on_reco_H_smear"].array()
		#root_dictionary["NW_on_reco_wlep_masses"] = root_file["neutrino_weighting/NW_on_reco_wlep_mass"].array()

		#root_dictionary["reco_lepton_lvecs"] = root_file["neutrino_weighting/reco_lepton_lvec"].array()
		#root_dictionary["reco_whad_lvecs"] = root_file["neutrino_weighting/reco_whad_lvec"].array()

		#root_dictionary["true_lepton_lvecs"] = root_file["neutrino_weighting/true_lepton_lvec"].array()
		#root_dictionary["true_whad_lvecs"] = root_file["neutrino_weighting/true_whad_lvec"].array()
		#root_dictionary["true_wlep_lvecs"] = root_file["neutrino_weighting/true_wlep_lvec"].array()

	# return dictionary
	return root_dictionary



def plot_2d(root_dictionary, tuple):
    # access tuple
    x_data = np.array(root_dictionary[tuple[0]])
    y_data = np.array(root_dictionary[tuple[1]])
    binning = tuple[2]
    xlabel = tuple[3]
    ylabel = tuple[4]
    normalize = tuple[5]  # boolean value to determine if normalization is needed

    # create plot
    plt.figure(figsize=(12, 8))

    # check if normalization is needed
    if normalize:
        # Calculate the weights to normalize the data
        weights = np.ones_like(x_data) / len(x_data)
    else:
        weights = None  # No normalization, use default weights

    # plot data with the optional weights
    hist, xedges, yedges, im = plt.hist2d(x_data, y_data, bins=binning, weights=weights)

    # add colorbar
    plt.colorbar(im)

    # beautificate    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # save fig
    plt.savefig(f"{PLOT_DESTINATION}2d_{tuple[0]}_{tuple[1]}.png")




if __name__ == '__main__':
	main()

