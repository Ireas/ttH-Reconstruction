import sys
import numpy as np
import uproot

import module_plot_mass_differences as plot_delta_m
import module_plot_pt_differences as plot_delta_pt
import module_plot_higgs_decay_modes as plot_higgs
import module_plot_jet_multiplicity as plot_njets
import module_plot_success as plot_success


OUTPUT_DESTINATION = "../output/plots/"


def main():
	# access given root file
	assert len(sys.argv)==2, "Error: root file must be given as only argument"
	root_file = uproot.open(sys.argv[1])
	
	
	# jet multiplicity
	#print("Verbose: plotting jet multiplicities")	
	#plot_njets.plot(root_file, OUTPUT_DESTINATION)
	
	# success
	print("Verbose: plotting success")	
	plot_success.plot(root_file, OUTPUT_DESTINATION)
	
	# mass difference
	print("Verbose: plotting mass difference")	
	plot_delta_m.plot(root_file, OUTPUT_DESTINATION)
	
	# mass difference
	#print("Verbose: plotting pt difference")	
	#plot_delta_pt.plot(root_file, OUTPUT_DESTINATION)

	# higgs decay mode
	print("Verbose: plotting higgs decay modes")	
	plot_higgs.plot(root_file, OUTPUT_DESTINATION)



if __name__=='__main__':
	main()	
