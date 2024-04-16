import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size=14)

COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
    '#5254a3', '#637939', '#8c6d31', '#843c39', '#7b4173'
]

class BarPlotSource:
	def __init__(self, data, err=None, label=None):
		self.data = data
		self.err = err
		self.label = label
	

class BarPlotOptions:
	def __init__(self, categories, title="title", x_label="x-axis", y_label="y-axis", y_scale="", wip=False, width=0.35, file_destination=None, show_plot=False):
		self.categories = categories
		self.title = title
		self.x_label = x_label
		self.y_label = y_label
		self.y_scale = y_scale
		self.width = width
		self.wip = wip
		self.file_destination = file_destination
		self.show_plot = show_plot


def plot_single_dataset(source, options):
	# determine bar positions
	bar_positions = np.arange(len(options.categories))

	# plot bars
	plt.bar(bar_positions, source.data, yerr=source.err, label=source.label, color=COLORS[0], alpha=0.7, width=options.width, capsize=2)

	# customize plot
	plt.title(options.title)
	plt.xlabel(options.x_label)
	plt.ylabel(options.y_label)
	plt.xticks(bar_positions, options.categories)

	if source.label:
		plt.legend()


	max_y = max(source.data)
	
	atlas_y = 0
	if options.y_scale=="log":
		atlas_y = 5*max_y
		plt.ylim([1, 10*max_y])
	else:
		atlas_y = 1.15*max_y
		plt.ylim([0, 1.2*max_y])

	if options.wip:
		plt.text(0, atlas_y, r"\textbf{\textit{ATLAS}} Simulation Work in Progress", ha='left', va='top', fontsize=14, color='black', usetex=True)
	
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
		plt.bar(bar_positions, source.data, yerr=source.err, label=source.label, color=COLORS[i], alpha=0.7, width=options.width, capsize=2)

	# customize plot
	plt.title(options.title)
	plt.xlabel(options.x_label)
	plt.ylabel(options.y_label)
	plt.xticks(bar_positions - options.width/2 * (len(sources)-1), options.categories)
	if options.y_scale:
		plt.yscale(options.y_scale)

	if source.label:
		plt.legend()
	
	max_y = 0
	for source in sources:
		if max(source.data)>max_y:
			max_y = max(source.data)
	
	atlas_y = 0
	if options.y_scale=="log":
		atlas_y = 5*max_y
		plt.ylim([1, 10*max_y])
	else:
		atlas_y = 1.1*max_y
		plt.ylim([0, 1.2*max_y])

	if options.wip:
		plt.text(0, atlas_y, r"\textbf{\textit{ATLAS}} Simulation Work in Progress", ha='left', va='top', fontsize=14, color='black', usetex=True)

	if options.file_destination:
		plt.savefig(options.file_destination)
	if options.show_plot:
		plt.show()
	else:
		plt.close()
	

