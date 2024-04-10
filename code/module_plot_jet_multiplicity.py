import numpy as np
import module_histograms as histogram


OUTPUT_DESTINATION = "./"
LABELS = { 
	-1: r"Other",
	6: r"6 Jets",
	7: r"7 Jets",
	8: r"8 Jets",
}



def plot_number_of_jets(root_file):
	# read from root_file
	m_reco_w_from_t = root_file["matched/reconstructed_W_from_t_m"].array()*1e-3
	m_reco_w_from_tbar = root_file["matched/reconstructed_W_from_tbar_m"].array()*1e-3
	jet_multiplicities= root_file["matched/number_of_jets"].array()
	

	# fill reconstructed masses by number of jets
	recos_by_jet_multiplicity = {}
	
	for i, jet_multiplicity in enumerate(jet_multiplicities):
		if not jet_multiplicity in recos_by_jet_multiplicity:
			recos_by_jet_multiplicity[jet_multiplicity] = []

		if m_reco_w_from_t[i]>0: 
			recos_by_jet_multiplicity[jet_multiplicity].append(m_reco_w_from_t[i]) 
		if m_reco_w_from_tbar[i]>0: 
			recos_by_jet_multiplicity[jet_multiplicity].append(m_reco_w_from_tbar[i]) 


	# combine multiple jet multiplicities
	recos_by_jet_multiplicity_reduced = {
		-1 : [],
		6 : [],
		7 : [],
		8 : [],
	}

	for (jet_multiplicity, masses) in recos_by_jet_multiplicity.items():
		if not jet_multiplicity in recos_by_jet_multiplicity_reduced:
			recos_by_jet_multiplicity_reduced[-1].extend(masses)
		else:
			recos_by_jet_multiplicity_reduced[jet_multiplicity].extend(masses)

		
	# create and plot containers
	for (jet_multiplicity, masses) in recos_by_jet_multiplicity_reduced.items():
		source = histogram.HistogramSource(
			data = masses,
			label = "yield = "+str(len(masses))
		)

		options = histogram.HistogramOptions(
			bins = np.arange(20,200,10),
			title = "Reconstructed Mass of W-Boson for " + LABELS[jet_multiplicity] + " (normalized)",
			x_label = r"$M_\text{reco}$ in GeV",
			y_label = r"Fraction of Events",
			normalize = True,
			file_destination= OUTPUT_DESTINATION + "mass_njets_" + str(jet_multiplicity) +".png"
		)

		histogram.plot_single_dataset(source, options)



def plot(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_number_of_jets(root_file)
