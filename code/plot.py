import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt



OUTPUT_DESTINATION = "../output/"
SHOW_PLOTS = False

COLORS= [
    'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
    'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
    'deepskyblue', 'darkorange', 'forestgreen', 'firebrick', 'mediumslateblue',
    'sienna', 'lightcoral', 'dimgrey', 'darkkhaki', 'mediumturquoise'
]




def main():
	# access given root file
	assert len(sys.argv)==2, "Error: root file must be given as only argument"
	root_file = uproot.open(sys.argv[1])

	generate_plots_success(root_file)
	generate_plots_mass(root_file)



def generate_plots_success(root_file):
	"""
	DESCRIPTION:
	Generates success plots for event reconstruction from given .root file.
	"""

	### PREPARE INPUT
	# read input from .root file
	successful_reconstruction = root_file['matched/successful_reconstruction'].array()
	

	# determine max number of jets
	number_of_jets = root_file['matched/number_of_jets'].array()
	max_number_of_jets = max(number_of_jets)
	
	
	# create empty dictionaries for each jet multiplicity
	dict_success = {}
	dict_unsuccess = {}
	categories = []
	for i in range(max_number_of_jets+1):
		dict_success[i] = 0
		dict_unsuccess[i] = 0
		categories+= [str(i)]



	### PROCESS INPUTS
	# fill dictionaries with success
	for success,jets in zip(successful_reconstruction, number_of_jets):
		if success:
			dict_success[jets]+= 1;
		else:
			dict_unsuccess[jets]+= 1;


	# calculate ratio
	ratio = []
	for success, unsuccess in zip(dict_success.values(), dict_unsuccess.values()):
		if unsuccess==0:
			ratio.append(0)
		else:
			ratio.append(success/unsuccess)



	### GENERATE PLOTS
	# plot successful events 
	plot_grouped_bar(
		data1 = dict_success.values(),
		data2 = dict_unsuccess.values(),
		label1 = "successful",
		label2 = "unsuccessful",
		categories = categories,
		title = "Successful Reconstruction by Number of Jets",
		xlabel = "Number of Jets",
		ylabel = "Number of Events",
		file_name = "successful_events"
	)
	

	# plot successful event ratio
	plot_bar(
		data = ratio,
		categories = categories,
		custom_xlims = [6,16],
		custom_ylims = [0,1],
		title = "Success Ratio by Number of Jets",
		xlabel = "Number of Jets",
		ylabel = "Fraction of Events",
		file_name="success_ratio"
	)



def generate_plots_mass(root_file):
	"""
	DESCRIPTION:
	Generates mass plots for W-boson and top quark from given .root file.
	"""

	### PREPARE INPUT
	# read data from .root file
	number_of_jets = root_file['matched/number_of_jets'].array()

	
	# convert MeV to GeV
	masses_W_from_t_truth = root_file['matched/truth_W_from_t_m'].array()*1e-3
	masses_W_from_t_reconstructed = root_file['matched/reconstructed_W_from_t_m'].array()*1e-3
	masses_W_from_tbar_truth = root_file['matched/truth_W_from_tbar_m'].array()*1e-3
	masses_W_from_tbar_reconstructed = root_file['matched/reconstructed_W_from_tbar_m'].array()*1e-3
	masses_t_truth = root_file['matched/truth_t_m'].array()*1e-3
	masses_t_reconstructed = root_file['matched/reconstructed_t_m'].array()*1e-3
	masses_tbar_truth = root_file['matched/truth_tbar_m'].array()*1e-3
	masses_tbar_reconstructed = root_file['matched/reconstructed_tbar_m'].array()*1e-3


	# create arrays for saving mass differences
	mass_difference_W = np.array([])	
	mass_difference_t = np.array([])	


	# create dictionaries for saving reco mass by jet multiplicity
	dict_m_w = {}
	dict_m_t = {}
	
	for i in range(max(number_of_jets)+1):
		dict_m_w[i] = np.array([])
		dict_m_t[i] = np.array([])
	

	# loop over .root file, skip unsuccessful reconstructions
	for m_truth, m_reco, n_jets in zip(masses_W_from_t_truth, masses_W_from_t_reconstructed, number_of_jets):
		if m_reco<1:
			continue
		mass_difference_W = np.append(mass_difference_W, m_reco - m_truth)
		dict_m_w[n_jets] = np.append(dict_m_w[n_jets], m_reco)
	
	for m_truth, m_reco, n_jets in zip(masses_W_from_tbar_truth, masses_W_from_tbar_reconstructed, number_of_jets):
		if m_reco<1:
			continue
		mass_difference_W = np.append(mass_difference_W, m_reco - m_truth)
		dict_m_w[n_jets] = np.append(dict_m_w[n_jets], m_reco)

	for m_truth, m_reco, n_jets in zip(masses_t_truth, masses_t_reconstructed, number_of_jets):
		if m_reco<1:
			continue
		mass_difference_t = np.append(mass_difference_t, m_reco - m_truth)
		dict_m_t[n_jets] = np.append(dict_m_t[n_jets], m_reco)
	
	for m_truth, m_reco, n_jets in zip(masses_tbar_truth, masses_tbar_reconstructed, number_of_jets):
		if m_reco<1:
			continue
		mass_difference_t = np.append(mass_difference_t, m_reco - m_truth)
		dict_m_t[n_jets] = np.append(dict_m_t[n_jets], m_reco)


	## PROCESS INPUT
	# create compact dictionaries that combine higher jet multiplicities	
	dict_m_w_compact = {
		"nJets=6" : np.array([]),
		"nJets=7" : np.array([]),
		"nJets=8" : np.array([]),
		"nJets>8" : np.array([])
	}

	dict_m_t_compact = {
		"nJets=6" : np.array([]),
		"nJets=7" : np.array([]),
		"nJets=8" : np.array([]),
		"nJets>8" : np.array([])
	}

	# fill compact dictionaries
	for i in range(max(number_of_jets)):
		if i==6:
			dict_m_w_compact["nJets=6"] = np.append(dict_m_w_compact["nJets=6"], dict_m_w[i])
			dict_m_t_compact["nJets=6"] = np.append(dict_m_t_compact["nJets=6"], dict_m_t[i])
		elif i==7:
			dict_m_w_compact["nJets=7"] = np.append(dict_m_w_compact["nJets=7"], dict_m_w[i])
			dict_m_t_compact["nJets=7"] = np.append(dict_m_t_compact["nJets=7"], dict_m_t[i])
		elif i==8:
			dict_m_w_compact["nJets=8"] = np.append(dict_m_w_compact["nJets=8"], dict_m_w[i])
			dict_m_t_compact["nJets=8"] = np.append(dict_m_t_compact["nJets=8"], dict_m_t[i])
		else:
			dict_m_w_compact["nJets>8"] = np.append(dict_m_w_compact["nJets>8"], dict_m_w[i])
			dict_m_t_compact["nJets>8"] = np.append(dict_m_t_compact["nJets>8"], dict_m_t[i])
	

	### GENERATE PLOTS 	
	# plot mass difference between truth and reco
	plot_histogram(
		data = mass_difference_W,
		bins = np.arange(-50,140,10), 
		title = "$\Delta M$ for W-bosons from (Anti-)Top", 
		xlabel = "$\Delta M$ in GeV", 
		ylabel = "Number of Events", 
		file_name = "delta_m_w"
	)
	
	plot_histogram(
		data = mass_difference_t,
		bins = np.arange(-80,200,10), 
		title = "$\Delta M$ for (Anti-)Top", 
		xlabel = "$\Delta M$ in GeV", 
		ylabel = "Number of Events", 
		file_name = "delta_m_t"
	)
	

	# plot reconstructed mass in respect to jet multiplicity
	plot_multidata_histogram(
		data_dict = dict_m_w_compact,
		fixed_bins = np.arange(20,220,20),
		stacked = True,
		title = "Reconstructed W-boson Mass by Jets",
		xlabel = "Mass in GeV",
		ylabel = "Number of Events",
		file_name = "multidata_m_w_by_jets"
	)

	plot_multidata_histogram(
		data_dict = dict_m_t_compact,
		fixed_bins = np.arange(100,340,20),
		stacked = True,
		title = "Reconstructed t Quark Mass by Jets",
		xlabel = "Mass in GeV",
		ylabel = "Number of Events",
		file_name = "multidata_m_t_by_jets"
	)



def plot_bar(data, categories, custom_xlims=None, custom_ylims=None, title='Bar Plot', xlabel='x-axis', ylabel='y-axis', width=0.35, file_name=None):	
	"""
	DESCRIPTION
	Plots simple bar plot and applies given labels.
	"""	

	# calculate the position for each group of bars
	bar_width = width
	bar_positions = np.arange(len(categories))
	
	# create bar plot
	plt.bar(bar_positions, data, color=COLORS[0], alpha=0.7, width=bar_width)
	
	# customize plot
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks(bar_positions, categories)
	if custom_xlims:
		plt.xlim([custom_xlims[0]-0.5, custom_xlims[1]+0.5])
	if custom_ylims:
		plt.ylim([custom_ylims[0], custom_ylims[1]])
    
	# save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")

    # show plot
	if SHOW_PLOTS:
		plt.show()

	# clean up
	plt.clf()



def plot_grouped_bar(data1, data2, categories, title='Grouped Bar Plot', xlabel='x-axis', ylabel='y-axis', label1='label1', label2='label2', width=0.35, file_name=None):	
	"""
	DESCRIPTION
	Plots grouped bar plot for two datasets side-by-side.
	"""	

	# calculate the position for each group of bars
	bar_width = width
	bar_positions1 = np.arange(len(categories))
	bar_positions2 = bar_positions1 + bar_width
	
	# create grouped bar plot
	plt.bar(bar_positions1, data1, color=COLORS[0], alpha=0.7, width=bar_width, label=label1)
	plt.bar(bar_positions2, data2, color=COLORS[1], alpha=0.7, width=bar_width, label=label2)
	
	# customize plot
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks(bar_positions1 + bar_width / 2, categories)
	plt.legend()
    
	# save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")

    # show plot
	if SHOW_PLOTS:
		plt.show()

	# clean up
	plt.clf()



def plot_histogram(data, bins=[], title='Histogram', xlabel="x-axis", ylabel="y-axis", yscale="linear", file_name=None):	
	"""
	DESCRIPTION
	Plots simple histrogram from data array. Bins must be array-like.
	"""	
	
	# creates histogram
	data = np.clip(data, bins[0], bins[-1])
	plt.hist(data, bins=bins, range=(bins[0], bins[-1]), color=COLORS[0], alpha=0.7, edgecolor="black")


    # customize plot
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.yscale(yscale)
	plt.xlim([bins[0],bins[-1]])
		
	# save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")
	
	# show plot
	if SHOW_PLOTS:
		plt.show()

	# clean up
	plt.clf()



def plot_multidata_histogram(data_dict, fixed_bins, stacked=True, title="Stacked Histogram", xlabel="x-axis", ylabel="y-label", yscale="linear",  file_name=None):
	"""
	DESCRIPTION
	Plots dictionary containing multiple dataset as (stacked) histrogram. Bins should be fixed manually. Dictionary is expected in the following format.
		{ [LABEL] : [DATA_ARRAY] }
	"""    

	# format input to arrays
	data_array = []
	labels_array = []
	colors_array = []
	for i, (key, value) in enumerate(data_dict.items()):
		data_array+= [value]
		labels_array+= [key]
		colors_array+= [COLORS[i]]

	# plot array of data
	plt.hist(data_array, bins=fixed_bins, label=labels_array, color=colors_array, alpha=0.7, edgecolor="black", stacked=stacked)
	
	# customize_plot
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.yscale(yscale)
	plt.xlim([fixed_bins[0],fixed_bins[-1]])
	plt.legend()

	# save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")

	# shot plot
	if SHOW_PLOTS:
		plt.show()

	# clean up
	plt.clf()



if __name__=='__main__':
	main()	
