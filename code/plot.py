import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt


# Example usage:
def main():
	assert len(sys.argv)==2, "Error: root file must be given as only argument"
	root_file = uproot.open(sys.argv[1])

	masses_t_truth = root_file['matched/truth_Tth_MC_W_from_t_m'].array()*1e6
	masses_t_reconstructed = root_file['matched/reconstructed_W_from_t_m'].array()*1e6
	masses_tbar_truth = root_file['matched/truth_Tth_MC_W_from_tbar_m'].array()*1e6
	masses_tbar_reconstructed = root_file['matched/reconstructed_W_from_tbar_m'].array()*1e6

	mass_difference = np.array([])	
	for m_truth, m_reco in zip(masses_t_truth, masses_t_reconstructed):
		if m_reco<1:
			continue
		mass_difference = np.append(mass_difference, m_reco - m_truth)
	
	for m_truth, m_reco in zip(masses_tbar_truth, masses_tbar_reconstructed):
		if m_reco<1:
			continue
		mass_difference = np.append(mass_difference, m_reco - m_truth)


	mass_difference_mev = mass_difference*1e-6
	mass_difference_gev = mass_difference*1e-9

	
	plot_histogram_and_save(
		mass_difference_gev, 
		bins=[-50,0,50,100,150,200,250,300], 
		plot_title="$\Delta M$ for W-boson reconstruction", 
		x_label="$\Delta M$ in GeV", 
		y_label='Number of Events', 
		save_path='delta_m_w.pdf'
	)


def plot_histogram_and_save(data, bins=10, lower_boundary=-1, upper_boundary=-1, plot_title='Histogram', x_label='X-axis', y_label='Frequency', save_path='histogram.pdf'):
	"""
	Create a histogram with common options, label axes, and save as a PDF.
	
	Parameters:
	- data (array-like): Input data for the histogram.
	- bins (int or array-like, optional): Specification of the bins (default is 10).
	- plot_title (str): Title of the plot (default is 'Histogram').
	- x_label (str): Label for the X-axis (default is 'X-axis').
	- y_label (str): Label for the Y-axis (default is 'Frequency').
	- save_path (str): Path to save the plot as a PDF (default is 'histogram.pdf').
	"""
	
	# Create a histogram
	if isinstance(bins, list) and len(bins)>=2:
		data = np.clip(data, bins[0], bins[-1])
		plt.hist(data, bins=bins, range=(bins[0],bins[-1]), edgecolor='black')
	else:
		plt.hist(data, bins=bins, edgecolor='black')

	# Set plot title and labels
	plt.title(plot_title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
		
	# Save the plot as a PDF file
	plt.savefig(save_path)
	
	# Show the plot (optional)
	plt.show()



if __name__=='__main__':
	main()	
