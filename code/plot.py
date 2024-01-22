import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt


OUTPUT_DESTINATION = "../output/"
SHOW_PLOTS = False

def main():
	# access given root file
	assert len(sys.argv)==2, "Error: root file must be given as only argument"
	root_file = uproot.open(sys.argv[1])


	# ===========  MASS PLOTS

	## read from .root file
	masses_W_from_t_truth = root_file['matched/truth_W_from_t_m'].array()*1e6
	masses_W_from_t_reconstructed = root_file['matched/reconstructed_W_from_t_m'].array()*1e6
	masses_W_from_tbar_truth = root_file['matched/truth_W_from_tbar_m'].array()*1e6
	masses_W_from_tbar_reconstructed = root_file['matched/reconstructed_W_from_tbar_m'].array()*1e6
	masses_t_truth = root_file['matched/truth_t_m'].array()*1e6
	masses_t_reconstructed = root_file['matched/reconstructed_t_m'].array()*1e6
	masses_tbar_truth = root_file['matched/truth_tbar_m'].array()*1e6
	masses_tbar_reconstructed = root_file['matched/reconstructed_tbar_m'].array()*1e6


	## calculate mass_difference for each successful reconstruction
	mass_difference_W = np.array([])	
	successfull_reconstructed_t_mass = np.array([])	
	successfull_reconstructed_tbar_mass = np.array([])	
	successfull_reconstructed_w_from_t_mass = np.array([])	
	successfull_reconstructed_w_from_tbar_mass = np.array([])	

	for m_truth, m_reco in zip(masses_W_from_t_truth, masses_W_from_t_reconstructed):
		if m_reco<1:
			continue
		mass_difference_W = np.append(mass_difference_W, m_reco - m_truth)
		successfull_reconstructed_w_from_t_mass = np.append(successfull_reconstructed_w_from_t_mass, m_reco)
	
	for m_truth, m_reco in zip(masses_W_from_tbar_truth, masses_W_from_tbar_reconstructed):
		if m_reco<1:
			continue
		mass_difference_W = np.append(mass_difference_W, m_reco - m_truth)
		successfull_reconstructed_w_from_tbar_mass = np.append(successfull_reconstructed_w_from_tbar_mass, m_reco)

	mass_difference_t = np.array([])	
	for m_truth, m_reco in zip(masses_t_truth, masses_t_reconstructed):
		if m_reco<1:
			continue
		mass_difference_t = np.append(mass_difference_t, m_reco - m_truth)
		successfull_reconstructed_t_mass = np.append(successfull_reconstructed_t_mass, m_reco)
	
	for m_truth, m_reco in zip(masses_tbar_truth, masses_tbar_reconstructed):
		if m_reco<1:
			continue
		mass_difference_t = np.append(mass_difference_t, m_reco - m_truth)
		successfull_reconstructed_tbar_mass = np.append(successfull_reconstructed_tbar_mass, m_reco)


	## transform to MeV and GeV for plotting
	mass_difference_W_mev = mass_difference_W*1e-6
	mass_difference_W_gev = mass_difference_W*1e-9
	mass_difference_t_mev = mass_difference_t*1e-6
	mass_difference_t_gev = mass_difference_t*1e-9
	successfull_reconstructed_t_mass_gev = successfull_reconstructed_t_mass*1e-9
	successfull_reconstructed_tbar_mass_gev = successfull_reconstructed_tbar_mass*1e-9
	successfull_reconstructed_w_from_t_mass_gev = successfull_reconstructed_w_from_t_mass*1e-9
	successfull_reconstructed_w_from_tbar_mass_gev = successfull_reconstructed_w_from_tbar_mass*1e-9


	
	## plot and save histrogram
	plot_histogram_and_save(
		mass_difference_W_gev, 
		bins=30, 
		custom_range=(-50,150),
		plot_title="$\Delta M$ for W-bosons from Top and Antitop Reconstruction (successful only)", 
		x_label="$\Delta M$ in GeV", 
		y_label='Number of Events', 
		file_name='delta_m_w',
		show_plot = SHOW_PLOTS
	)
	
	plot_histogram_and_save(
		mass_difference_t_gev, 
		bins=30, 
		custom_range=(-100,200),
		plot_title="$\Delta M$ for Top and Antitop Reconstruction (successful only)", 
		x_label="$\Delta M$ in GeV", 
		y_label='Number of Events', 
		file_name='delta_m_t',
		show_plot = SHOW_PLOTS
	)
	
	plot_histogram_and_save(
		successfull_reconstructed_t_mass_gev, 
		bins=30, 
		custom_range=(100,300),
		plot_title="Reconstructed $M_t$ for Top Quark (successful only)", 
		x_label="$M_t$ in GeV", 
		y_label='Number of Events', 
		file_name='m_t',
		show_plot = SHOW_PLOTS
	)
	
	plot_histogram_and_save(
		successfull_reconstructed_tbar_mass_gev, 
		bins=30, 
		custom_range=(100,300),
		plot_title="Reconstructed $M_t$ for Antitop Quark (successful only)", 
		x_label="$M_t$ in GeV", 
		y_label='Number of Events', 
		file_name='m_tbar',
		show_plot = SHOW_PLOTS
	)
	
	plot_histogram_and_save(
		successfull_reconstructed_w_from_t_mass_gev, 
		bins=30, 
		custom_range=(30,180),
		plot_title="Reconstructed $M_W$ for W-Boson from Top Quark (successful only)", 
		x_label="$M_W$ in GeV", 
		y_label='Number of Events', 
		file_name='m_W_from_t',
		show_plot = SHOW_PLOTS
	)
	
	plot_histogram_and_save(
		successfull_reconstructed_w_from_tbar_mass_gev, 
		bins=30, 
		custom_range=(30,180),
		plot_title="Reconstructed $M_W$ for W-Boson from Antitop Quark (successful only)", 
		x_label="$M_W$ in GeV", 
		y_label='Number of Events', 
		file_name='m_W_from_tbar',
		show_plot = SHOW_PLOTS
	)


	# ==========  SUCCESS PLOTS
	successful_reconstruction = root_file['matched/successful_reconstruction'].array()
	number_of_jets = root_file['matched/number_of_jets'].array()
	max_number_of_jets = max(number_of_jets)
	
	# prepare dictionaries
	dict_success = {}
	dict_unsuccess = {}
	categories = []
	for i in range(max_number_of_jets+1):
		dict_success[i] = 0
		dict_unsuccess[i] = 0
		categories+= [str(i)]

	# fill dictionaries
	for success,jets in zip(successful_reconstruction, number_of_jets):
		if success:
			dict_success[jets]+= 1;
		else:
			dict_unsuccess[jets]+= 1;

	plot_grouped_bar_and_save(
		data1 = dict_success.values(),
		data2 = dict_unsuccess.values(),
		labels = categories,
		title = "Successful Reconstruction by Number of Jets",
		xlabel = "Number of Jets",
		ylabel = "Number of Events",
		label1 = "successful",
		label2 = "unsuccessful",
		file_name="success",
		show_plot = SHOW_PLOTS
	)
	
	ratio = []
	for success, unsuccess in zip(dict_success.values(), dict_unsuccess.values()):
		if unsuccess==0:
			ratio.append(0)
		else:
			ratio.append(success/unsuccess)

	plot_bar_and_save(
		data = ratio,
		labels = categories,
		custom_xlims = [6,16],
		custom_ylims = [0,1],
		title = "Success Ratio by Number of Jets",
		xlabel = "Number of Jets",
		ylabel = "Fraction of Events",
		file_name="success_ratio",
		show_plot = SHOW_PLOTS
	)





def plot_histogram_and_save(data, bins=10, custom_range=None, plot_title='Histogram', x_label='X-axis', y_label='Frequency', file_name=None, show_plot=False):	
	# create histogram
	if isinstance(bins, list) and len(bins)>=2:
		data = np.clip(data, bins[0], bins[-1])
		plt.hist(data, bins=bins, range=(bins[0],bins[-1]), edgecolor='black')
	elif custom_range:
		plt.hist(data, bins=bins, range=custom_range, edgecolor='black')
	else:
		plt.hist(data, bins=bins, edgecolor='black')

    # Customize plot options
	plt.title(plot_title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
		
	# Save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")
	
	# Show plot
	if show_plot:
		plt.show()
	plt.clf()


def plot_grouped_bar_and_save(data1, data2, labels, color1='skyblue', color2='orange', title='Grouped Bar Plot', xlabel='Categories', ylabel='Values', label1='label1', label2='label2', width=0.35, rotation=0, file_name=None, show_plot=False):	
	# Calculate the position for each group of bars
	bar_width = width
	bar_positions1 = np.arange(len(labels))
	bar_positions2 = bar_positions1 + bar_width
	
	# Create grouped bar plot
	plt.bar(bar_positions1, data1, color=color1, width=bar_width, label=label1)
	plt.bar(bar_positions2, data2, color=color2, width=bar_width, label=label2)
	
	# Customize plot options
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks(bar_positions1 + bar_width / 2, labels, rotation=rotation)
	plt.legend()
    
	# Save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")

    # Show plot
	if show_plot:
		plt.show()
	plt.clf()

def plot_bar_and_save(data, labels, custom_xlims=None, custom_ylims=None, color='skyblue', title='Bar Plot', xlabel='Categories', ylabel='Values', label=None, width=0.35, rotation=0, file_name=None, show_plot=False):	
	# Calculate the position for each group of bars
	bar_width = width
	bar_positions = np.arange(len(labels))
	
	# Create grouped bar plot
	plt.bar(bar_positions, data, color=color, width=bar_width, label=label)
	
	# Customize plot options
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks(bar_positions, labels, rotation=rotation)
	if custom_xlims:
		plt.xlim([custom_xlims[0]-0.5, custom_xlims[1]+0.5])
	if custom_ylims:
		plt.ylim([custom_ylims[0], custom_ylims[1]])
	if label:
		plt.legend()
    
	# Save plot
	if file_name:
		plt.savefig(OUTPUT_DESTINATION+file_name+".png")

    # Show plot
	if show_plot:
		plt.show()
	plt.clf()


if __name__=='__main__':
	main()	
