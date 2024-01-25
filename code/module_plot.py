import numpy as np
import matplotlib.pyplot as plt

COLOR = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
    '#5254a3', '#637939', '#8c6d31', '#843c39', '#7b4173'
]

def plot_histogram(data, bins, normalize=False, title="Histogram", x_label="x-axis", y_label="y-axis", data_label=None, y_scale="linear", y_range=None, grid=False, file_name=None, show_plot=False):
	"""
	Plot a histogram.

	Parameters:
	    data (array-like): Input data for the histogram.
	    bins (int or array-like): Specification of histogram bins.
	    normalize (bool, optional): If True, normalize the histogram to form a probability density. Default is False.
	    title (str, optional): Title of the plot. Default is "Histogram".
	    x_label (str, optional): Label for the x-axis. Default is "x-axis".
	    y_label (str, optional): Label for the y-axis. Default is "y-axis".
	    data_label (str, optional): Label for the histogram data. Default is None.
	    y_scale (str, optional): Scale for the y-axis. Default is "linear".
	    y_range (tuple, optional): Tuple specifying the range of the y-axis.
	    grid (bool, optional): If True, display grid lines. Default is False.
	    file_name (str, optional): If provided, save the plot to the specified file.
	    show_plot (bool, optional): If True, display the plot. If False, close the plot. Default is False.
	"""

	# Preprocessing
	if data is None or bins is None:
		raise ValueError("Mandatory parameters 'data' and 'bins' cannot be None.")
	data = np.clip(data, bins[0], bins[-1])

	# Plotting
	plt.figure()

	plt.hist(data, bins=bins, density=normalize, label=data_label, edgecolor='black', alpha=0.8, color=COLOR[0])

	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)

	plt.xlim(bins[0], bins[-1])
	plt.yscale(y_scale)
	
	if y_range:
		plt.ylim(y_range)
	if grid:
		plt.grid(True)
	if data_label:
		plt.legend()

	# Output
	if file_name:
		plt.savefig(file_name)
	if show_plot:
		plt.show()
	else:
		plt.close()
