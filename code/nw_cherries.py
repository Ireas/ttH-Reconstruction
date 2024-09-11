import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


# CONSTATNTS
CHERRY_INDICIES = [14]#np.arange(2,3,1)

PLOT_DESTINATION = "/home/ireas/git_repos/master/plots/NW_cherries/"
ROOT_INPUT_FILE = "/media/ireas/Data/v6/weighted/ttHWW_full_only_2000eta_2000m_50e.root"


NW_MASS_RANGE = [0,50] # Mass in GeV
NW_ETA_RANGE = [-3,3]
CHERRY_ETA_BINS = 2000
CHERRY_MASS_BINS = 2000


plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

COLORMAP = 'magma'
HIGHLIGHT_COLOR = 'red'


def create_dictionary(root_file_directory):	
	# create empty dictionary
	root_dictionary = {}

	# open .root file
	root_file = uproot.open(root_file_directory)
	
	# fill dictionary
	root_dictionary["NW_vec_weights"] = root_file["neutrino_weighting/NW_full_weight"].array()
	root_dictionary["NW_vec_wlep_masses"] = root_file["neutrino_weighting/NW_full_wlep_mass"].array()
	root_dictionary["NW_vec_nu_etas"] = root_file["neutrino_weighting/NW_full_nu_eta"].array()

	root_dictionary["true_wlep_m"] = root_file["neutrino_weighting/true_wlep_m"].array()
	root_dictionary["true_nu_eta"] = root_file["neutrino_weighting/true_neutrino_eta"].array()

	#root_dictionary["NW_weights"] = root_file["neutrino_weighting/NW_weight"].array()
	#root_dictionary["NW_nu_pxs"] = root_file["neutrino_weighting/NW_nu_px"].array()
	#root_dictionary["NW_nu_pys"] = root_file["neutrino_weighting/NW_nu_py"].array()

	# close .root file afterwards
	root_file.close()

	# return dictionary
	return root_dictionary


def plot_cherry(root_dictionary, index):
	global roi_delta_wlep_mass, roi_delta_nu_eta, delta_wlep_mass, delta_nu_eta
	# access values
	# neutrino weighting output
	weights = np.array(root_dictionary["NW_vec_weights"][index])
	wlep_masses = np.array(root_dictionary["NW_vec_wlep_masses"][index])
	nu_etas = np.array(root_dictionary["NW_vec_nu_etas"][index])

	# true values
	true_nu_eta = root_dictionary["true_nu_eta"][index]
	true_wlep_mass= root_dictionary["true_wlep_m"][index]

	# highest rated values
	highest_nw_weight_wlep_mass = wlep_masses[weights==max(weights)][0]
	highest_nw_weight_nu_eta = nu_etas[weights==max(weights)][0]
	
	# skip cherry if no solution is found anywhere
	if max(weights)<0:
		return


	# replace invalid solutions (plot them as white)
	weights[weights<=0] = np.nan
	

	# create plot
	plt.figure(figsize=(10,8))

	# plot 2d cherry distribution
	hb = plt.hist2d(wlep_masses/1e3, nu_etas, bins=[CHERRY_MASS_BINS+1,CHERRY_ETA_BINS+1], weights=weights, range=[NW_MASS_RANGE,NW_ETA_RANGE])
	cb = plt.colorbar(hb[3])
	cb.set_label(r"NW Weight")

	# mark best and true solution
	plt.scatter(highest_nw_weight_wlep_mass/1e3, highest_nw_weight_nu_eta, color=HIGHLIGHT_COLOR, label='Best Solution RoI', marker='*', s=120)
	plt.scatter(true_wlep_mass/1e3, true_nu_eta, color=HIGHLIGHT_COLOR, label='True Values', marker='X', s=120)
	
	# customise
	plt.xlabel(r"Sampled $m_{W*}$ [GeV]")
	plt.ylabel(r"Sampled $\eta_\nu$")
	plt.legend()

	# save file, dont show
	plt.savefig(f"{PLOT_DESTINATION}cherry_{index}.png")
	plt.close()


def heatmap_nw(root_dictionary):
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
	hb = plt.hist2d(heat_delta_wlep_mass/1e3, heat_delta_nu_eta, bins=[np.linspace(-30,30,20), np.linspace(-2,2,20)], cmap=COLORMAP)
	cb = plt.colorbar(hb[3])
	cb.set_label('Number of Events')
	plt.xlabel("Difference $\\Delta(M_\\text{est} - M_\\text{true}$) for leptonic $W$-boson GeV")
	plt.ylabel("Difference $\\Delta(\\eta_\\text{est} - \\eta_\\text{true}$) for $\\nu$")
	plt.savefig(PLOT_DESTINATION+"heatmap-nw.png")



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
	hb = plt.hist2d(heat_delta_px/1e3, heat_delta_py/1e3, bins=[np.linspace(-50,50,20), np.linspace(-50,50,20)], cmap=COLORMAP)
	cb = plt.colorbar(hb[3])
	cb.set_label('Number of Events')
	plt.xlabel("Difference $\\Delta(p_\\text{est} - p_\\text{true})_x$ for $\\nu$ in GeV")
	plt.ylabel("Difference $\\Delta(p_\\text{est} - p_\\text{true})_y$ for $\\nu$ in GeV")
	plt.savefig(PLOT_DESTINATION+"heatmap-nu-pt.png")



if __name__ == '__main__':
	# create dictionary
	print("create dictionary")
	root_dictionary = create_dictionary(ROOT_INPUT_FILE)
	
	# generate Delta Heatmap
	#heatmap_nw(root_dictionary)
	#heatmap_nu(root_dictionary)

	# generate 2D distribution cherries
	print("plot cherries")
	for cherry_index in CHERRY_INDICIES:
		print(f" cherry: {cherry_index} / {len(CHERRY_INDICIES)}")
		plot_cherry(root_dictionary, cherry_index)

