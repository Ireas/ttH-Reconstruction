import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


# ==========  VALIDATE NEUTRINO WEIGHTING  =========
# ==================================================
WHAD_MASS_PDG = 80.4335
WHAD_MASS_MAX_DEVIATION = 10 # maximum valid diviation from PDG mass for on-shell definition

CHERRY_INDICIES = [13]

PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/"

def create_dictionary(root_file_directory):	
	# create empty dictionary
	root_dictionary = {}

	# open .root file
	root_file = uproot.open(root_file_directory)
	
	# fill dictionary
	root_dictionary["NW_vec_weights"] = root_file["neutrino_weighting/NW_full_weight"].array()
	root_dictionary["NW_vec_wlep_masses"] = root_file["neutrino_weighting/NW_full_wlep_mass"].array()
	root_dictionary["NW_vec_nu_etas"] = root_file["neutrino_weighting/NW_full_nu_eta"].array()

	root_dictionary["true_wlep_lvecs"] = root_file["neutrino_weighting/true_wlep_lvec"].array()
	root_dictionary["true_nu_lvecs"] = root_file["neutrino_weighting/true_neutrino_lvec"].array()

	root_dictionary["NW_weights"] = root_file["neutrino_weighting/NW_weight"].array()
	root_dictionary["NW_nu_pxs"] = root_file["neutrino_weighting/NW_nu_px"].array()
	root_dictionary["NW_nu_pys"] = root_file["neutrino_weighting/NW_nu_py"].array()	


	# close .root file afterwards
	root_file.close()

	# return dictionary
	return root_dictionary


def plot_cherry(root_dictionary, index):
	# access
	weights = np.array(root_dictionary["NW_vec_weights"][index])
	wlep_masses = np.array(root_dictionary["NW_vec_wlep_masses"][index])
	nu_etas = np.array(root_dictionary["NW_vec_nu_etas"][index])

	true_nu_eta = root_dictionary["true_nu_lvecs"][index].fCoordinates.tolist().get("fCoordinates.fEta")
	true_wlep_mass= root_dictionary["true_wlep_lvecs"][index].fCoordinates.tolist().get("fCoordinates.fM")


	# get highest values
	highest_nw_weight_wlep_mass = wlep_masses[weights==max(weights)][0]
	highest_nw_weight_nu_eta = nu_etas[weights==max(weights)][0]
	



	# replace invalid solutions (plot them as white)
	weights[weights<=0] = np.nan


	# plot
	plt.figure()
	hb = plt.hist2d(wlep_masses/1e3, nu_etas, bins=[501,601], weights=weights, range=[[0,50],[-3,3]])
	
	cb = plt.colorbar(hb[3])
	cb.set_label('Weight')

	
	plt.scatter(highest_nw_weight_wlep_mass/1e3, highest_nw_weight_nu_eta, color='red', label='Best Solution', marker='*', s=90)
	plt.scatter(true_wlep_mass/1e3, true_nu_eta, color='red', label='True Values', marker='X', s=90)
	
	plt.xlabel("Sampled $W_{lep}$ mass M")
	plt.ylabel("Sampled $\\eta_{\\nu}$")

	plt.legend()
	plt.savefig(PLOT_DESTINATION+"cherries/cherry_" + str(index) + ".png")


def heatmap(root_dictionary):
	# access
	weights_vecs = np.array(root_dictionary["NW_vec_weights"])
	wlep_masses_vecs = np.array(root_dictionary["NW_vec_wlep_masses"])
	nu_etas_vecs = np.array(root_dictionary["NW_vec_nu_etas"])
	true_nu_lvecs = root_dictionary["true_nu_lvecs"]
	true_wlep_lvecs = root_dictionary["true_wlep_lvecs"]
	weights = root_dictionary["NW_weights"]


	# plotting stuff
	heat_delta_wlep_mass = np.array([])
	heat_delta_nu_eta = np.array([])



	# event loop
	for (weights, wlep_masses, nu_etas, true_nu_lvec, true_wlep_lvec, weight) in zip(weights_vecs, wlep_masses_vecs, nu_etas_vecs, true_nu_lvecs, true_wlep_lvecs, weights):
		if weight<0:
			continue

		# get true values
		true_wlep_mass = true_wlep_lvec.fCoordinates.tolist().get("fCoordinates.fM")
		true_nu_eta = true_nu_lvec.fCoordinates.tolist().get("fCoordinates.fEta")

		# get best estimate
		highest_nw_weight_wlep_mass = wlep_masses[weights==max(weights)][0]
		highest_nw_weight_nu_eta = nu_etas[weights==max(weights)][0]

		# fill heat data
		heat_delta_wlep_mass = np.append( heat_delta_wlep_mass, [highest_nw_weight_wlep_mass-true_wlep_mass] )
		heat_delta_nu_eta = np.append( heat_delta_nu_eta, [highest_nw_weight_nu_eta-true_nu_eta] )


	# plot
	plt.figure(figsize=(7,5))
	hb = plt.hist2d(heat_delta_wlep_mass/1e3, heat_delta_nu_eta , bins=[np.linspace(-50,50,25), np.linspace(-6,3,10)])
	cb = plt.colorbar(hb[3])
	cb.set_label('Number of Events')
	plt.xlabel("Difference $\\Delta(M_\\text{est} - M_\\text{true}$) for leptonic $W$-boson GeV")
	plt.ylabel("Difference $\\Delta(\\eta_\\text{est} - \\eta_\\text{true}$) for $\\nu$")
	plt.savefig(PLOT_DESTINATION+"heatmap.png")



def heatmap_nu(root_dictionary):
	# access
	true_nu_lvecs = root_dictionary["true_nu_lvecs"]
	weights = root_dictionary["NW_weights"]
	estimated_nu_pxs = root_dictionary["NW_nu_pxs"]
	estimated_nu_pys = root_dictionary["NW_nu_pys"]


	# plotting stuff
	heat_delta_px = np.array([])
	heat_delta_py = np.array([])



	# event loop
	for (true_nu_lvec, estimated_nu_px, estimated_nu_py, weight) in zip(true_nu_lvecs, estimated_nu_pxs, estimated_nu_pys, weights):
		if weight<0:
			continue
		
		# get true values
		true_nu_pt = true_nu_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		true_nu_phi = true_nu_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		
		true_nu_px = true_nu_pt * np.cos( true_nu_phi )
		true_nu_py = true_nu_pt * np.sin( true_nu_phi )

		# fill heat data
		heat_delta_px = np.append( heat_delta_px, [estimated_nu_px-true_nu_px] )
		heat_delta_py = np.append( heat_delta_py, [estimated_nu_py-true_nu_py] )


	# plot
	plt.figure(figsize=(7,5))
	hb = plt.hist2d(heat_delta_px, heat_delta_py , bins=[np.linspace(-300,300,25), np.linspace(-300,300,25)])
	cb = plt.colorbar(hb[3])
	cb.set_label('Number of Events')
	plt.xlabel("Difference $\\Delta(p_\\text{x,est} - p_\\text{x,true}$) for neutrino in Mev")
	plt.ylabel("Difference $\\Delta(p_\\text{y,est} - p_\\text{y,true}$) for neutrino in Mev")
	plt.savefig(PLOT_DESTINATION+"heatmap_nu.png")



if __name__ == '__main__':
	## INPUT 
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for validation was given, exiting")
		exit()
	

	# create dictionary
	root_dictionary = create_dictionary(sys.argv[1])
	
	# generate Delta Heatmap
	#heatmap(root_dictionary)

	#heatmap_nu(root_dictionary)

	# generate 2D distribution cherries
	for cherry_index in CHERRY_INDICIES:
		plot_cherry(root_dictionary, cherry_index)
	
