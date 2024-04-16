import numpy as np
import module_histograms as histogram


OUTPUT_DESTINATION = "./"
LABELS = { 
	-1: r"Other",
	6: r"6 Jets",
	7: r"7 Jets",
	8: r"8 Jets",
}


def plot_mass_difference(root_file):
	# read from root_file
	m_reco_w_from_t = root_file["matched/reconstructed_W_from_t_m"].array()*1e-3
	m_reco_w_from_tbar = root_file["matched/reconstructed_W_from_tbar_m"].array()*1e-3
	m_truth_w_from_t = root_file["matched/truth_W_from_t_m"].array()*1e-3
	m_truth_w_from_tbar = root_file["matched/truth_W_from_tbar_m"].array()*1e-3

	m_reco_t = root_file["matched/reconstructed_t_m"].array()*1e-3
	m_reco_tbar = root_file["matched/reconstructed_tbar_m"].array()*1e-3
	m_truth_t = root_file["matched/truth_t_m"].array()*1e-3
	m_truth_tbar = root_file["matched/truth_tbar_m"].array()*1e-3
	

	# fill reconstructed masses by number of jets
	mass_difference_w = []

	for m_reco, m_truth in zip(m_reco_w_from_t, m_truth_w_from_t):
		if m_reco<1:
			continue
		mass_difference_w.append( m_reco-m_truth )
	
	for m_reco, m_truth in zip(m_reco_w_from_tbar, m_truth_w_from_tbar):
		if m_reco<1:
			continue
		mass_difference_w.append( m_reco-m_truth )
	
		
	# create and plot containers
	source = histogram.HistogramSource(
		data = mass_difference_w,
		#label = "yield = "+str(len(mass_difference_w))
	)

	options = histogram.HistogramOptions(
		bins = np.arange(-50, 120, 10),
		title = r"$\Delta M$ of W-Boson",
		x_label = r"$\Delta M$ in GeV",
		y_label = r"Number of Events",
		normalize = False,
		wip = True,
		file_destination= OUTPUT_DESTINATION + "delta_mass_w.png"
	)

	histogram.plot_single_dataset(source, options)



	# fill reconstructed masses by number of jets
	mass_difference_t = []

	for m_reco, m_truth in zip(m_reco_t, m_truth_t):
		if m_reco<1:
			continue
		mass_difference_t.append( m_reco-m_truth )
	
	for m_reco, m_truth in zip(m_reco_tbar, m_truth_tbar):
		if m_reco<1:
			continue
		mass_difference_t.append( m_reco-m_truth )
	
		
	# create and plot containers
	source = histogram.HistogramSource(
		data = mass_difference_t
		#label = "yield = "+str(len(mass_difference_w))
	)

	options = histogram.HistogramOptions(
		bins = np.arange(-70, 160, 10),
		title = r"$\Delta M$ of top quark",
		x_label = r"$\Delta M$ in GeV",
		y_label = r"Number of Events",
		normalize = False,
		wip = True,
		file_destination= OUTPUT_DESTINATION + "delta_mass_t.png"
	)

	histogram.plot_single_dataset(source, options)




def plot(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_mass_difference(root_file)
