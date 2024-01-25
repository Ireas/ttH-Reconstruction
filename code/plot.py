import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

import module_plot_mass_differences as plot_delta
import module_plot_higgs_decay_modes as plot_higgs
import module_plot_jet_multiplicity as plot_njets


OUTPUT_DESTINATION = "../output/plots/"
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
	
	# mass difference
	print("Verbose: verify mass difference")	
	plot_delta.verify(root_file, OUTPUT_DESTINATION)

	# higgs decay mode
	print("Verbose: verify higgs decay modes")	
	plot_higgs.verify(root_file, OUTPUT_DESTINATION)
	
	# jet multiplicity
	print("Verbose: verify jet multiplicities")	
	plot_njets.verify(root_file, OUTPUT_DESTINATION)
	
	# old
	generate_plots_success(root_file)
	#generate_plots_mass(root_file)



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



if __name__=='__main__':
	main()	
