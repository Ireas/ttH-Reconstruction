import numpy as np
import module_plot as plot


OUTPUT_DESTINATION = "./"
LABELS = { 
	-1: r"Other",
	6: r"6 Jets",
	7: r"7 Jets",
	8: r"8 Jets",
}



class ContainerMassW:
	def __init__(self, data):
		self.data = data

	def plot_histogram(self):
		source = plot.HistogramSource(
			data = self.data,
			label = "yield = "+str(len(self.data))
		)

		options = plot.HistogramOptions(
			bins = np.arange(-50, 120, 10),
			title = "Mass Difference of Reconstructed W-Boson",
			x_label = r"$\Delta M_\text{reco}$ in GeV",
			y_label = r"Number of Events",
			normalize = False,
			file_destination= OUTPUT_DESTINATION + "mass_w_delta.png"
		)

		plot.histogram(source, options)
	


def plot_mass_difference(root_file):
	# read from root_file
	m_reco_w_from_t = root_file["matched/reconstructed_W_from_t_m"].array()*1e-3
	m_reco_w_from_tbar = root_file["matched/reconstructed_W_from_tbar_m"].array()*1e-3
	m_truth_w_from_t = root_file["matched/truth_W_from_t_m"].array()*1e-3
	m_truth_w_from_tbar = root_file["matched/truth_W_from_tbar_m"].array()*1e-3
	

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
	container = ContainerMassW(mass_difference_w)
	container.plot_histogram()



def verify(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_mass_difference(root_file)
