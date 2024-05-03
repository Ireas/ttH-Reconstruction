import numpy as np
import module_bar_plots as bars


OUTPUT_DESTINATION = "./"
LABELS = { 
	-1: r"Other",
	8: r"8 Jets",
	9: r"9 Jets",
	10: r"10 Jets",
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

	for (jet_multiplicity, success_multiplicity, unsuccess_multiplicity) in zip(success_by_jets.keys(), success_by_jets.values(), unsuccess_by_jets.values()):
		if (success_multiplicity + unsuccess_multiplicity)>0:
			ratio_by_jets[jet_multiplicity] = success_multiplicity/(success_multiplicity + unsuccess_multiplicity)
		else:
			ratio_by_jets[jet_multiplicity] = 0


	ratio_by_jets_filtered = {a:b for (a,b) in ratio_by_jets.items() if (a>=8 and a<=16)}
	success_by_jets_filtered = {a:b for (a,b) in success_by_jets.items() if (a>=8 and a<=16)}
	unsuccess_by_jets_filtered = {a:b for (a,b) in unsuccess_by_jets.items() if (a>=8 and a<=16)}
		
	# create and plot containers
	ratio_err = []
	for (ratio, succ_event, unsucc_event) in zip(ratio_by_jets_filtered.values(), success_by_jets_filtered.values(), unsuccess_by_jets_filtered.values()):
		ratio_err.append( ratio*(1-ratio)/(succ_event+unsucc_event) ) # binomial error

	ratio_err = np.sqrt(np.array(ratio_err))

	source = bars.BarPlotSource(
		data = ratio_by_jets_filtered.values(),
		err = ratio_err
	)

	options = bars.BarPlotOptions(
		categories = ratio_by_jets_filtered.keys(),
		title = "Success Ratio by Jet Multiplicity",
		x_label = r"Number of Jets",
		y_label = r"Ratio of Events",
		file_destination= OUTPUT_DESTINATION + "success_ratio.png"
	)

	bars.plot_single_dataset(source, options)



	sources = []
	sources.append(
		bars.BarPlotSource(
			data = success_by_jets_filtered.values(),
			err = np.sqrt( np.array( list(success_by_jets_filtered.values())) ),
			label = "success" 
		)
	)
	sources.append(
		bars.BarPlotSource(
			data = unsuccess_by_jets_filtered.values(),
			err = np.sqrt( np.array( list(unsuccess_by_jets_filtered.values())) ),
			label = "failure"
		)
	)
	
	options = bars.BarPlotOptions(
		categories = success_by_jets_filtered.keys(),
		title = "Successes by Jet Multiplicity",
		x_label = r"Number of Jets",
		y_label = r"Number of Events",
		y_scale = "log",
		file_destination= OUTPUT_DESTINATION + "success_by_jets.png"
	)

	bars.plot_multiple_datasets(sources, options)


def plot_success_all_hww(root_file):
	# read from root_file
	success_ratings = root_file["matched/successful_reconstruction"].array()
	jet_multiplicities = root_file["matched/number_of_jets"].array()
	root_tree = root_file["matched"]

	r_h_dm = root_tree["higgs_decay_mode_custom"].array()
	r_h_d1_d1_pdgids = root_tree["higgs_decay1_decay1_filtered_pdgid"].array()
	r_h_d1_d2_pdgids = root_tree["higgs_decay1_decay2_filtered_pdgid"].array()
	r_h_d2_d1_pdgids = root_tree["higgs_decay2_decay1_filtered_pdgid"].array()
	r_h_d2_d2_pdgids = root_tree["higgs_decay2_decay2_filtered_pdgid"].array()

	# preprocess data
	ratio_success_jets_all = {}
	ratio_success_jets_hww = {}
	total_all = 0
	total_hww = 0

	for i in range(max(jet_multiplicities)+1):
		ratio_success_jets_all[i] = 0
		ratio_success_jets_hww[i] = 0


	for (jet_multiplicity, success, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(jet_multiplicities, success_ratings, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
		if not success:
			continue

		ratio_success_jets_all[jet_multiplicity]+= 1
		total_all+= 1

      # ww-decay
		if h_dm!=5:
			continue

        # invalid
		if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
			continue
        
        # full hadronic
		if ( (pdgid_11>0 and pdgid_11<9) and (pdgid_12>0 and pdgid_12<9) and (pdgid_21>0 and pdgid_21<9) and (pdgid_22>0 and pdgid_22<9) ):
			continue

        # full leptonic
		if ( (pdgid_11>10 and pdgid_11<19) and (pdgid_12>10 and pdgid_12<19) and (pdgid_21>10 and pdgid_21<19) and (pdgid_22>10 and pdgid_22<19) ):
			continue

		ratio_success_jets_hww[jet_multiplicity]+= 1
		total_hww+= 1

	
	print(ratio_success_jets_all.items())
	print(ratio_success_jets_hww.items())
	print(total_all)
	print(total_hww)

	for (key,value) in ratio_success_jets_all.items():
		ratio_success_jets_all[key]/= total_all

	for (key,value) in ratio_success_jets_hww.items():
		ratio_success_jets_hww[key]/= total_hww

	print(ratio_success_jets_all.items())
	print(ratio_success_jets_hww.items())	


	ratio_success_jets_all_fltered = {a:b for (a,b) in ratio_success_jets_all.items() if (a>6 and a<10)}
	ratio_success_jets_hww_fltered = {a:b for (a,b) in ratio_success_jets_hww.items() if (a>6 and a<10)}
	
	print(ratio_success_jets_all_fltered.items())
	print(ratio_success_jets_hww_fltered.items())

	# create and plot containers
	sources = []
	sources.append( bars.BarPlotSource( data=ratio_success_jets_all_fltered.values(), label=r"all events" ) )
	sources.append( bars.BarPlotSource( data=ratio_success_jets_hww_fltered.values(), label=r"semi-leptonic $H\rightarrow WW$" ) )
	
	options = bars.BarPlotOptions(
		categories = ratio_success_jets_all_fltered.keys(),
		title = "Successes Ratios by Jet Multiplicity",
		x_label = r"Number of Jets",
		y_label = r"Ratio of Events",
		y_lim = 1,
		file_destination= OUTPUT_DESTINATION + "success_ratio_event_type.png"
	)

	bars.plot_multiple_datasets(sources, options)



def plot(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_success(root_file)
	#plot_success_all_hww(root_file)
