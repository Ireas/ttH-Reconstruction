import numpy as np
import module_histograms as histogram


OUTPUT_DESTINATION = "./"


def plot_pt_difference(root_file):
	# read from root_file
	reco_w_from_t_pt = root_file["matched/reconstructed_W_from_t_pt"].array()*1e-3
	reco_w_from_tbar_pt = root_file["matched/reconstructed_W_from_tbar_pt"].array()*1e-3
	truth_w_from_t_pt = root_file["matched/truth_W_from_t_pt"].array()*1e-3
	truth_w_from_tbar_pt = root_file["matched/truth_W_from_tbar_pt"].array()*1e-3
	
	reco_w_from_t_m = root_file["matched/reconstructed_W_from_t_m"].array()*1e-3
	reco_w_from_tbar_m = root_file["matched/reconstructed_W_from_tbar_m"].array()*1e-3
	

	# fill reconstructed masses by number of jets
	pt_difference_w = []

	for m_reco, pt_reco, pt_truth in zip(reco_w_from_t_m, reco_w_from_t_pt, truth_w_from_t_pt):
		if m_reco<1e-3:
			continue
		pt_difference_w.append( pt_reco-pt_truth )
	
	for m_reco, pt_reco, pt_truth in zip(reco_w_from_tbar_m, reco_w_from_tbar_pt, truth_w_from_tbar_pt):
		if m_reco<1e-3:
			continue
		pt_difference_w.append( pt_reco-pt_truth )
	
		
	# create and plot containers
	source = histogram.HistogramSource(
		data = pt_difference_w,
		#label = "yield = "+str(len(pt_difference_w))
	)

	options = histogram.HistogramOptions(
		bins = np.arange(-60, 120, 10),
		title = "$\Delta p_t$ of W-Boson (Reconstructed - Truth)",
		x_label = r"$\Delta p_t$ in GeV",
		y_label = r"Number of Events",
		normalize = False,
		file_destination= OUTPUT_DESTINATION + "delta_pt_w.png"
	)

	histogram.plot_single_dataset(source, options)



def plot(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_pt_difference(root_file)
