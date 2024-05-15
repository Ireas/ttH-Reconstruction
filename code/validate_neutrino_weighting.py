import sys
import uproot # root in python
import numpy as np


# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
#
# usage: python3 convert_root_to_h5.py [ROOT_FILE_NEUTRINO_WEIGHTING]



# main method 
def main():
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for validation was given, exiting")
		exit()

	# open .root file
	root_file = uproot.open(sys.argv[1])
	
	
	# fill .h5 file with .root information


	# get number of jets and events, apply limit on number of events
	true_lvecs_neutrino = root_file['neutrino_weighting/lvec_neutrino_true'].array()
	nw_eta_nu = root_file['neutrino_weighting/estimated_eta_nu'].array()
		


	m = 0
	#for (true,esti) in zip(true_lvecs_neutrino, esti_lvecs_neutrino):
	for true in true_lvecs_neutrino:
		if(m>20):
			exit()
		
	
		#true_pt  = true.fCoordinates.tolist().get('fCoordinates.fPt')
		true_eta = true.fCoordinates.tolist().get('fCoordinates.fEta')
		#true_phi = true.fCoordinates.tolist().get('fCoordinates.fPhi')
		#true_m   = true.fCoordinates.tolist().get('fCoordinates.fM')
		
		#esti_pt  = esti.fCoordinates.tolist().get('fCoordinates.fPt')
		#esti_eta = esti.fCoordinates.tolist().get('fCoordinates.fEta')
		#esti_phi = esti.fCoordinates.tolist().get('fCoordinates.fPhi')
		#esti_m   = esti.fCoordinates.tolist().get('fCoordinates.fM')

		current_eta_nu = nw_eta_nu[m] 

		if(true_eta==0):
			continue

		print("==========")
		#print("pt: ", true_pt, " | ", esti_pt )
		print("eta: ", np.round(true_eta,3), " | ", np.round(current_eta_nu,3) )
		#print("phi: ", true_phi, " | ", esti_phi )
		#print("m: ", true_m, " | ", esti_m )
		m+= 1





if __name__ == '__main__':
	main()
