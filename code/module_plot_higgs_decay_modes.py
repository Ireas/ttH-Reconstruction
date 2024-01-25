import numpy as np

import module_plot as plot


OUTPUT_DESTINATION = "./"

LABELS= {
	0 : r"$bb$",
	1 : r"$ee$",
	2 : r"$\mu\mu$",
	3 : r"$\tau\tau$",
	4 : r"$\gamma\gamma$",
	5 : r"$WW$",
	6 : r"$ZZ$",
	8 : r"Other",
}



class ContainerMassWDecayMode:
	def __init__(self, data, decay_mode_id):
		self.data = data
		self.decay_mode_id = decay_mode_id

	def plot_histogram(self):
		source = plot.HistogramSource(
			data = self.data, 
			label = "yield = "+str(len(self.data))
		)

		options = plot.HistogramOptions(
			bins = np.arange(20,200,10), 
			title = "Reconstructed Mass of W-Boson for " + LABELS[self.decay_mode_id] + " H Decay Mode (normalized)",
			x_label = r"$M_\text{reco}$ in GeV",
			y_label = r"Fraction of Events",
			normalize = True,
			file_destination = OUTPUT_DESTINATION + "mass_higgs_dm_" + str(self.decay_mode_id) +".png"
		)
		
		plot.histogram(source, options)
	


def check_higgs_decay_modes(root_file):
	m_reco_w_from_t = root_file["matched/reconstructed_W_from_t_m"].array()*1e-3
	m_reco_w_from_tbar = root_file["matched/reconstructed_W_from_tbar_m"].array()*1e-3

	dm_ids = root_file["matched/higgs_decay_mode_custom"].array()
	
	recos_by_dm_id = {}

	for i, dm_id in enumerate(dm_ids):
		if dm_id==-1:
			continue

		if not dm_id in recos_by_dm_id:
			recos_by_dm_id[dm_id] = []
		

		if m_reco_w_from_t[i]>0: 
			recos_by_dm_id[dm_id].append(m_reco_w_from_t[i]) 
		if m_reco_w_from_tbar[i]>0: 
			recos_by_dm_id[dm_id].append(m_reco_w_from_tbar[i]) 


	containers = []
	for (dm_id, masses) in recos_by_dm_id.items():
		container = ContainerMassWDecayMode(masses, dm_id)
		container.plot_histogram()
		containers.append(container)	



def verify(root_file, output_destination):
	global OUTPUT_DESTINATION
	OUTPUT_DESTINATION = output_destination
	check_higgs_decay_modes(root_file)
