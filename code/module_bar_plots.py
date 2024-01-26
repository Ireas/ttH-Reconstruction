import numpy as np
import matplotlib.pyplot as plt

COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
    '#5254a3', '#637939', '#8c6d31', '#843c39', '#7b4173'
]

class BarPlotSource:
	def __init__(self, data, label=None):
		self.data = data
		self.label = label
	

class BarPlotOptions:
	def __init__(self, categories, title="title", x_label="x-axis", y_label="y-axis", width=0.35, file_destination=None, show_plot=False):
		self.categories = categories
		self.title = title
		self.x_label = x_label
		self.y_label = y_label
		self.width = width
		self.file_destination = file_destination
		self.show_plot = show_plot


def plot_single_dataset(source, options):
	# determine bar positions
	bar_positions = np.arange(len(options.categories))

	# plot bars
	plt.bar(bar_positions, source.data, label=source.label, color=COLORS[0], alpha=0.7, width=options.width)

	# customize plot
	plt.title(options.title)
	plt.xlabel(options.x_label)
	plt.ylabel(options.y_label)
	plt.xticks(bar_positions, options.categories)

	if source.label:
		plt.legend()
	
	if options.file_destination:
		plt.savefig(options.file_destination)
	if options.show_plot:
		plt.show()
	else:
		plt.close()


def plot_multiple_datasets(sources, options):
	# determine bar positions
	bar_positions_center = np.arange(len(options.categories))

	# plot bars
	for (i, source) in enumerate(sources):
		bar_positions = np.array([float(x) for x in bar_positions_center])
		bar_positions+= i*options.width		
		plt.bar(bar_positions, source.data, label=source.label, color=COLORS[i], alpha=0.7, width=options.width)

	# customize plot
	plt.title(options.title)
	plt.xlabel(options.x_label)
	plt.ylabel(options.y_label)
	plt.xticks(bar_positions - options.width/2 * (len(sources)-1), options.categories)

	if source.label:
		plt.legend()
	
	if options.file_destination:
		plt.savefig(options.file_destination)
	if options.show_plot:
		plt.show()
	else:
		plt.close()
	

