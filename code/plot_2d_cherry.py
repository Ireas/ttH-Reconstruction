import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


# ==========  VALIDATE NEUTRINO WEIGHTING  =========
# ==================================================
WHAD_MASS_PDG = 80.4335
WHAD_MASS_MAX_DEVIATION = 10 # maximum valid diviation from PDG mass for on-shell definition

CHERRY_INDEX = 2


def create_dictionary(root_file_directory):	
	# create empty dictionary
	root_dictionary = {}

	# open .root file
	root_file = uproot.open(root_file_directory)
	
	# fill dictionary
	root_dictionary["NW_vec_weights"] = root_file["neutrino_weighting/NW_on_true_vec_weights"].array()
	root_dictionary["NW_vec_wlep_masses"] = root_file["neutrino_weighting/NW_on_true_vec_wlep_masses"].array()
	root_dictionary["NW_vec_nu_etas"] = root_file["neutrino_weighting/NW_on_true_vec_nu_etas"].array()
	root_dictionary["true_wlep_lvecs"] = root_file["neutrino_weighting/true_wlep_lvec"].array()
	root_dictionary["true_nu_lvecs"] = root_file["neutrino_weighting/true_neutrino_lvec"].array()

	# close .root file afterwards
	root_file.close

	# return dictionary
	return root_dictionary


def plot_2d(root_dictionary, index):
	# access
	weights = np.array(root_dictionary["NW_vec_weights"][index])
	wlep_masses = np.array(root_dictionary["NW_vec_wlep_masses"][index])
	nu_etas = np.array(root_dictionary["NW_vec_nu_etas"][index])

	true_nu_eta = root_dictionary["true_nu_lvecs"][index].fCoordinates.tolist().get("fCoordinates.fEta")
	true_wlep_mass= root_dictionary["true_wlep_lvecs"][index].fCoordinates.tolist().get("fCoordinates.fM")


	# get highest values
	highest_nw_weight_wlep_mass = wlep_masses[weights==max(weights)][0]
	highest_nw_weight_nu_eta = nu_etas[weights==max(weights)][0]
	

	# print values
	print(true_nu_eta)
	print(true_wlep_mass)
	
	print(highest_nw_weight_wlep_mass)
	print(highest_nw_weight_nu_eta)


	# replace invalid solutions (plot them as white)
	weights[weights<=0] = np.nan


	# plot
	plt.figure()
	hb = plt.hist2d(wlep_masses/1e3, nu_etas, bins=[251,121], weights=weights, range=[[0,50],[-3,3]])
	
	cb = plt.colorbar(hb[3])
	cb.set_label('Weight')

	
	plt.scatter(highest_nw_weight_wlep_mass/1e3, highest_nw_weight_nu_eta, color='red', label='Best Solution', marker='*', s=90)
	plt.scatter(true_wlep_mass/1e3, true_nu_eta, color='red', label='True Values', marker='X', s=90)
	
	plt.xlabel("Sampled $W_{lep}$ mass M")
	plt.ylabel("Sampled $\\eta_{\\nu}$")

	plt.legend()
	plt.savefig("cherry_" + str(index) + ".png")
	plt.clf()



if __name__ == '__main__':
	## INPUT 
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for validation was given, exiting")
		exit()
	

	# create dictionary
	root_dictionary = create_dictionary(sys.argv[1])
	
	# call methods
	for i in range(0,100,1):
		plot_2d(root_dictionary, i)

