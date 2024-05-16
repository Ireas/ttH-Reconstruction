import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


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
	nw_weight = root_file['neutrino_weighting/neutrino_weight'].array()
	
	val_py =  root_file['neutrino_weighting/val_px'].array()
	val_py = root_file['neutrino_weighting/val_py'].array()
	met_value= root_file['neutrino_weighting/missing_energy_value'].array()
	met_phi = root_file['neutrino_weighting/missing_energy_phi'].array()
		


	d_eta = np.array([])
	d_px = np.array([])
	d_py = np.array([])


	#for (true,esti) in zip(true_lvecs_neutrino, esti_lvecs_neutrino):
	for true_lvec, weight, estimated_eta, px, py, met_v, met_p in zip(true_lvecs_neutrino, nw_weight, nw_eta_nu, val_py, val_py, met_value, met_phi):
		
		#true_pt  = true.fCoordinates.tolist().get('fCoordinates.fPt')
		true_eta = true_lvec.fCoordinates.tolist().get('fCoordinates.fEta')
		#true_phi = true.fCoordinates.tolist().get('fCoordinates.fPhi')
		#true_m   = true.fCoordinates.tolist().get('fCoordinates.fM')
		
		#esti_pt  = esti.fCoordinates.tolist().get('fCoordinates.fPt')
		#esti_eta = esti.fCoordinates.tolist().get('fCoordinates.fEta')
		#esti_phi = esti.fCoordinates.tolist().get('fCoordinates.fPhi')
		#esti_m   = esti.fCoordinates.tolist().get('fCoordinates.fM')
		if (true_eta==0):
			continue

		if (weight<-1):
			continue	
		

		
		print("==========")
		print("eta: ", np.round(true_eta,3), " | ", np.round(estimated_eta,3)," | ", weight )

		true_px = met_v * np.cos(met_p)
		true_py = met_v * np.sin(met_p)

	
		d_eta = np.append(d_eta, [estimated_eta-true_eta])
		d_px = np.append(d_px, [px-true_px])
		d_py = np.append(d_py, [py-true_py])
		



	plt.title("delta eta")
	plt.hist(d_eta, np.arange(-3,3,0.5))
	plt.savefig("delta_eta.png")
	plt.show()
	plt.title("delta px")
	plt.hist(d_px, np.arange(-1e5, 1e5, 2e4))
	plt.savefig("delta_px.png")
	plt.show()
	plt.title("delta py")
	plt.hist(d_py, np.arange(-2e4, 2e4, 5e3))
	plt.savefig("delta_py.png")
	plt.show()






if __name__ == '__main__':
	main()
