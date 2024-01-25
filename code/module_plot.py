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


class HistogramSource:
	def __init__(self, data, label=None):
		self.data = data
		self.label = label


class HistogramOptions:
	def __init__(self, bins, title="title", x_label="x-axis", y_label="y-axis", y_scale="linear", y_range=None, normalize=False, grid=False, file_destination=None, show_plot=False):
		self.bins = bins
		self.title = title
		self.x_label = x_label
		self.y_label = y_label
		self.y_scale = y_scale
		self.y_range = y_range
		self.normalize = normalize		
		self.grid = grid
		self.file_destination = file_destination
		self.show_plot = show_plot


def histogram(source, options):
	# Preprocessing
	if source==None or not isinstance(source, HistogramSource):
		raise ValueError("Error: Source not properly set.")
	if options==None or not isinstance(options, HistogramOptions):
		raise ValueError("Error: Options not properly set.")

	data = np.clip(source.data, options.bins[0], options.bins[-1])

	# Plotting
	plt.figure()

	plt.hist(data, bins=options.bins, density=options.normalize, label=source.label, edgecolor='black', alpha=0.8, color=COLOR[0])

	
	# Customize
	plt.title(options.title)
	plt.xlabel(options.x_label)
	plt.ylabel(options.y_label)
	plt.yscale(options.y_scale)
	
	plt.xlim(options.bins[0], options.bins[-1])
	
	if options.y_range:
		plt.ylim(options.y_range)
	if options.grid:
		plt.grid(True)
	if source.label:
		plt.legend()

	# Output
	if options.file_destination:
		plt.savefig(options.file_destination)
	if options.show_plot:
		plt.show()
	else:
		plt.close()
