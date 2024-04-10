import numpy as np
import module_bar_plots as bars


OUTPUT_DESTINATION = "./"
LABELS = { 
	-1: r"Other",
	6: r"6 Jets",
	7: r"7 Jets",
	8: r"8 Jets",
}


def plot_success(root_file):
	# read from root_file
	success_ratings = root_file["matched/successful_reconstruction"].array()
	jet_multiplicities = root_file["matched/number_of_jets"].array()
	

	# preprocess data
	success_by_jets = {}
	unsuccess_by_jets = {}
	ratio_by_jets = {}
	total = 0

	for i in range(max(jet_multiplicities)+1):
		success_by_jets[i] = 0
		unsuccess_by_jets[i] = 0
		ratio_by_jets[i] = 0

	for (jet_multiplicity, success) in zip(jet_multiplicities, success_ratings):
		if success:
			success_by_jets[jet_multiplicity]+= 1
		else:	
			unsuccess_by_jets[jet_multiplicity]+= 1
		total+= 1
	
	if total==0:
		print("Warning: No events successful, exiting")
		return	

	for (jet_multiplicity, success_multiplicity) in success_by_jets.items():
		ratio_by_jets[jet_multiplicity] = success_multiplicity/total


		
	# create and plot containers
	source = bars.BarPlotSource(
		data = ratio_by_jets.values(),
	)

	options = bars.BarPlotOptions(
		categories = ratio_by_jets.keys(),
		title = "Success Ratio by Jet Multiplicity",
		x_label = r"Number of Jets",
		y_label = r"Ratio of Events",
		file_destination= OUTPUT_DESTINATION + "success_ratio.png"
	)

	bars.plot_single_dataset(source, options)

	
	sources = []
	sources.append( bars.BarPlotSource( data=success_by_jets.values(), label="success" ) )
	sources.append( bars.BarPlotSource( data=unsuccess_by_jets.values(), label="unsuccess" ) )
	
	options = bars.BarPlotOptions(
		categories = success_by_jets.keys(),
		title = "Successes by Jet Multiplicity",
		x_label = r"Number of Jets",
		y_label = r"Number of Events",
		file_destination= OUTPUT_DESTINATION + "success_by_jets.png"
	)

	bars.plot_multiple_datasets(sources, options)



def plot(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_success(root_file)
