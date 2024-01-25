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
	def __init__(self, data, jet_multiplicity):
		self.data = data
		self.jet_multiplicity = jet_multiplicity

	def plot_histogram(self):
		source = plot.HistogramSource(
			data = self.data,
			label = "yield = "+str(len(self.data))
		)

		options = plot.HistogramOptions(
			bins = np.arange(20,200,10),
			title = "Reconstructed Mass of W-Boson for " + LABELS[self.jet_multiplicity] + " (normalized)",
			x_label = r"$M_\text{reco}$ in GeV",
			y_label = r"Fraction of Events",
			normalize = True,
			file_destination= OUTPUT_DESTINATION + "mass_njets_" + str(self.jet_multiplicity) +".png"
		)

		plot.histogram(source, options)
	


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
		container = ContainerMassW(masses, jet_multiplicity)
		container.plot_histogram()



def verify(root_file, output_destination):
	global OUTPUT_DESTINATION 
	OUTPUT_DESTINATION = output_destination
	plot_number_of_jets(root_file)
