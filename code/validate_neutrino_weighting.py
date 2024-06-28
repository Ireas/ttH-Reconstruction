import sys
import uproot # root in python
import numpy as np
import matplotlib.pyplot as plt


# ==========  VALIDATE NEUTRINO WEIGHTING  =========
# ==================================================
WHAD_MASS_PDG = 80.4335
WHAD_MASS_MAX_DEVIATION = 10 # maximum valid diviation from PDG mass for on-shell definition


def create_dictionary(root_file_directory):	
	# create empty dictionary
	root_dictionary = {}

	# open .root file
	root_file = uproot.open(root_file_directory)
	
	# fill dictionary
	root_dictionary["NW_on_true_weights"] = root_file["neutrino_weighting/NW_on_true_weight"].array()
	root_dictionary["NW_on_reco_weights"] = root_file["neutrino_weighting/NW_on_reco_weight"].array()
	
	root_dictionary["NW_on_true_wlep_masses"] = root_file["neutrino_weighting/NW_on_true_wlep_mass"].array()
	root_dictionary["NW_on_reco_wlep_masses"] = root_file["neutrino_weighting/NW_on_reco_wlep_mass"].array()
	
	root_dictionary["true_lepton_lvecs"] = root_file["neutrino_weighting/true_lepton_lvec"].array()
	root_dictionary["reco_lepton_lvecs"] = root_file["neutrino_weighting/reco_lepton_lvec"].array()
	
	root_dictionary["true_whad_lvecs"] = root_file["neutrino_weighting/true_whad_lvec"].array()
	root_dictionary["reco_whad_lvecs"] = root_file["neutrino_weighting/reco_whad_lvec"].array()
	
	root_dictionary["true_wlep_lvecs"] = root_file["neutrino_weighting/true_wlep_lvec"].array()
	
	root_dictionary["event_abs_lepton_pdgids"] = root_file["neutrino_weighting/signature_abs_lepton_pdgid"].array()
	root_dictionary["classifier_whad_possible"] = root_file["neutrino_weighting/classifier_whad_possible"].array()


	# close .root file afterwards
	root_file.close

	# return dictionary
	return root_dictionary


def val_reco_vs_truth_lepton_success(root_dictionary):
	# access
	NW_on_reco_weights = root_dictionary["NW_on_reco_weights"]
	reco_lepton_lvecs = root_dictionary["reco_lepton_lvecs"]
	true_abs_lepton_pdgids = root_dictionary["event_abs_lepton_pdgids"]

	# prepare
	plot_lepton_success_pts = np.array([])
	plot_lepton_success_etas = np.array([])
	plot_lepton_success_phis = np.array([])
	plot_lepton_fail_pts = np.array([])
	plot_lepton_fail_etas = np.array([])
	plot_lepton_fail_phis = np.array([])
	
	n_electron_success = 0
	n_electron_fail = 0
	n_muon_success = 0
	n_muon_fail = 0

	n_electron_total = 0
	n_muon_total = 0
	n_total = 0
		

	# fill arrays 
	for (reco_lepton_lvec, lepton_pdgid, weight) in zip(reco_lepton_lvecs, true_abs_lepton_pdgids, NW_on_reco_weights):
		reco_lepton_pt = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		reco_lepton_eta = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		reco_lepton_phi = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		
		if(weight<0):
			if(lepton_pdgid==11):
				n_electron_fail+= 1
				n_electron_total+= 1
			if(lepton_pdgid==13):
				n_muon_fail+= 1
				n_muon_total+= 1
			
			plot_lepton_fail_pts = np.append(plot_lepton_fail_pts, [reco_lepton_pt])
			plot_lepton_fail_etas = np.append(plot_lepton_fail_etas, [reco_lepton_eta])
			plot_lepton_fail_phis = np.append(plot_lepton_fail_phis, [reco_lepton_phi])
		else:
			if(lepton_pdgid==11):
				n_electron_success+= 1
				n_electron_total+= 1
			if(lepton_pdgid==13):
				n_muon_success+= 1
				n_muon_total+= 1
			
			plot_lepton_success_pts = np.append(plot_lepton_success_pts, [reco_lepton_pt])
			plot_lepton_success_etas = np.append(plot_lepton_success_etas, [reco_lepton_eta])
			plot_lepton_success_phis = np.append(plot_lepton_success_phis, [reco_lepton_phi])
		
		n_total+= 1


	# plot
	plt.xlabel(r"Lepton $p_t$ [GeV]")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_lepton_success_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"NW success")
	plt.hist(plot_lepton_fail_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"NW failed")
	plt.legend()
	plt.savefig("out_success_lepton_pt.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\eta$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_lepton_success_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"NW success")
	plt.hist(plot_lepton_fail_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"NW failed")
	plt.legend()
	plt.savefig("out_success_lepton_eta.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\phi$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_lepton_success_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"NW success")
	plt.hist(plot_lepton_fail_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"NW failed")
	plt.legend()
	plt.savefig("out_success_lepton_phi.png")
	plt.clf()
	
	plt.xlabel(r"Lepton Flavour")
	plt.ylabel(r"Relative Number of Events")
	plt.bar( 
		["electron", "muon"], 
		[n_electron_success/n_electron_total, n_muon_success/n_muon_total], 
		label="success",
		bottom = 0
	)
	plt.bar( 
		["electron", "muon"], 
		[n_electron_fail/n_electron_total, n_muon_fail/n_muon_total], 
		label="fail",
		bottom = [n_electron_success/n_electron_total, n_muon_success/n_muon_total]
	)
	plt.legend()
	plt.savefig("out_success_lepton_flavour.png")
	plt.clf()


def val_reco_vs_truth_lepton(root_dictionary):
	# access
	true_lepton_lvecs = root_dictionary["true_lepton_lvecs"]
	reco_lepton_lvecs = root_dictionary["reco_lepton_lvecs"]
	true_abs_lepton_pdgids = root_dictionary["event_abs_lepton_pdgids"]

	# prepare
	plot_true_electron_pts = np.array([])
	plot_true_electron_etas = np.array([])
	plot_true_electron_phis = np.array([])
	plot_true_electron_energies = np.array([])
	plot_reco_electron_pts = np.array([])
	plot_reco_electron_etas = np.array([])
	plot_reco_electron_phis = np.array([])
	plot_reco_electron_energies = np.array([])
	
	plot_true_muon_pts = np.array([])
	plot_true_muon_etas = np.array([])
	plot_true_muon_phis = np.array([])
	plot_true_muon_energies = np.array([])
	plot_reco_muon_pts = np.array([])
	plot_reco_muon_etas = np.array([])
	plot_reco_muon_phis = np.array([])
	plot_reco_muon_energies = np.array([])
	
	
	# fill arrays 
	for (true_lepton_lvec, reco_lepton_lvec, lepton_pdgid) in zip(true_lepton_lvecs, reco_lepton_lvecs, true_abs_lepton_pdgids):
		# true lepton
		true_lepton_pt = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		true_lepton_eta = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		true_lepton_phi = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		true_lepton_energy = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fE")	
		reco_lepton_pt = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		reco_lepton_eta = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		reco_lepton_phi = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		reco_lepton_energy = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fE")

		if(lepton_pdgid==11):	
			plot_true_electron_pts = np.append(plot_true_electron_pts, [true_lepton_pt])
			plot_true_electron_etas = np.append(plot_true_electron_etas, [true_lepton_eta])
			plot_true_electron_phis = np.append(plot_true_electron_phis, [true_lepton_phi])
			plot_true_electron_energies = np.append(plot_true_electron_energies, [true_lepton_energy])
			plot_reco_electron_pts = np.append(plot_reco_electron_pts, [reco_lepton_pt])
			plot_reco_electron_etas = np.append(plot_reco_electron_etas, [reco_lepton_eta])
			plot_reco_electron_phis = np.append(plot_reco_electron_phis, [reco_lepton_phi])
			plot_reco_electron_energies = np.append(plot_reco_electron_energies, [reco_lepton_energy])
		
		if(lepton_pdgid==13):	
			plot_true_muon_pts = np.append(plot_true_muon_pts, [true_lepton_pt])
			plot_true_muon_etas = np.append(plot_true_muon_etas, [true_lepton_eta])
			plot_true_muon_phis = np.append(plot_true_muon_phis, [true_lepton_phi])
			plot_true_muon_energies = np.append(plot_true_muon_energies, [true_lepton_energy])
			plot_reco_muon_pts = np.append(plot_reco_muon_pts, [reco_lepton_pt])
			plot_reco_muon_etas = np.append(plot_reco_muon_etas, [reco_lepton_eta])
			plot_reco_muon_phis = np.append(plot_reco_muon_phis, [reco_lepton_phi])
			plot_reco_muon_energies = np.append(plot_reco_muon_energies, [reco_lepton_energy])
		


	# plot
	plt.xlabel(r"Lepton $p_t$ [GeV]")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_electron_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"true value (electron)")
	plt.hist(plot_reco_electron_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"reco value (electron)")
	plt.hist(plot_true_muon_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"true value (muon)")
	plt.hist(plot_reco_muon_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"reco value (muon)")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_flavour_pt.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\eta$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_electron_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"true value (electron)")
	plt.hist(plot_reco_electron_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"reco value (electron)")
	plt.hist(plot_true_muon_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"true value (muon)")
	plt.hist(plot_reco_muon_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"reco value (muon)")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_flavour_eta.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\phi$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_electron_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"true value (electron)")
	plt.hist(plot_reco_electron_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"reco value (electron)")
	plt.hist(plot_true_muon_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"true value (muon)")
	plt.hist(plot_reco_muon_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"reco value (muon)")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_flavour_phi.png")
	plt.clf()



	
def val_reco_vs_truth_kinematics(root_dictionary):
	# access
	true_lepton_lvecs = root_dictionary["true_lepton_lvecs"]
	reco_lepton_lvecs = root_dictionary["reco_lepton_lvecs"]
	true_whad_lvecs = root_dictionary["true_whad_lvecs"]
	reco_whad_lvecs = root_dictionary["reco_whad_lvecs"]
	
	whad_possibles = root_dictionary["classifier_whad_possible"]

	# prepare
	plot_true_lepton_pts = np.array([])
	plot_true_lepton_etas = np.array([])
	plot_true_lepton_phis = np.array([])
	plot_true_lepton_energies = np.array([])
	plot_reco_lepton_pts = np.array([])
	plot_reco_lepton_etas = np.array([])
	plot_reco_lepton_phis = np.array([])
	plot_reco_lepton_energies = np.array([])
	
	plot_true_whad_pts = np.array([])
	plot_true_whad_etas = np.array([])
	plot_true_whad_phis = np.array([])
	plot_true_whad_energies = np.array([])
	plot_reco_whad_pts = np.array([])
	plot_reco_whad_etas = np.array([])
	plot_reco_whad_phis = np.array([])
	plot_reco_whad_energies = np.array([])

	
	# fill arrays 
	for (true_lepton_lvec, reco_lepton_lvec, true_whad_lvec, reco_whad_lvec, whad_possible) in zip(true_lepton_lvecs, reco_lepton_lvecs, true_whad_lvecs, reco_whad_lvecs, whad_possibles):
		# true lepton
		true_lepton_pt = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		true_lepton_eta = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		true_lepton_phi = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		true_lepton_energy = true_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fE")	

		plot_true_lepton_pts = np.append(plot_true_lepton_pts, [true_lepton_pt])
		plot_true_lepton_etas = np.append(plot_true_lepton_etas, [true_lepton_eta])
		plot_true_lepton_phis = np.append(plot_true_lepton_phis, [true_lepton_phi])
		plot_true_lepton_energies = np.append(plot_true_lepton_energies, [true_lepton_energy])
		
		# reco lepton	
		reco_lepton_pt = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		reco_lepton_eta = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		reco_lepton_phi = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		reco_lepton_energy = reco_lepton_lvec.fCoordinates.tolist().get("fCoordinates.fE")
		
		plot_reco_lepton_pts = np.append(plot_reco_lepton_pts, [reco_lepton_pt])
		plot_reco_lepton_etas = np.append(plot_reco_lepton_etas, [reco_lepton_eta])
		plot_reco_lepton_phis = np.append(plot_reco_lepton_phis, [reco_lepton_phi])
		plot_reco_lepton_energies = np.append(plot_reco_lepton_energies, [reco_lepton_energy])

		if (whad_possible<0):
			continue

		# true whad
		true_whad_pt = true_whad_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		true_whad_eta = true_whad_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		true_whad_phi = true_whad_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		true_whad_energy = true_whad_lvec.fCoordinates.tolist().get("fCoordinates.fE")	
		
		plot_true_whad_pts = np.append(plot_true_whad_pts, [true_whad_pt])
		plot_true_whad_etas = np.append(plot_true_whad_etas, [true_whad_eta])
		plot_true_whad_phis = np.append(plot_true_whad_phis, [true_whad_phi])
		plot_true_whad_energies = np.append(plot_true_whad_energies, [true_whad_energy])
		
		# reco whad	
		reco_whad_pt = reco_whad_lvec.fCoordinates.tolist().get("fCoordinates.fPt")
		reco_whad_eta = reco_whad_lvec.fCoordinates.tolist().get("fCoordinates.fEta")
		reco_whad_phi = reco_whad_lvec.fCoordinates.tolist().get("fCoordinates.fPhi")
		reco_whad_energy = reco_whad_lvec.fCoordinates.tolist().get("fCoordinates.fE")
		
		plot_reco_whad_pts = np.append(plot_reco_whad_pts, [reco_whad_pt])
		plot_reco_whad_etas = np.append(plot_reco_whad_etas, [reco_whad_eta])
		plot_reco_whad_phis = np.append(plot_reco_whad_phis, [reco_whad_phi])
		plot_reco_whad_energies = np.append(plot_reco_whad_energies, [reco_whad_energy])


	# plot
	plt.xlabel(r"Lepton $p_t$ [GeV]")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_lepton_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_lepton_pts/1e3, np.linspace(0,200,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_pt.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\eta$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_lepton_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_lepton_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_eta.png")
	plt.clf()
	
	plt.xlabel(r"Lepton $\phi$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_lepton_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_lepton_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_lepton_phi.png")
	plt.clf()

	# energies not supported by truth objects by default	
	#plt.xlabel(r"Lepton $E$ [GeV]")
	#plt.ylabel(r"Number of Events")
	#plt.hist(plot_true_lepton_energies, np.linspace(0,400,20), histtype="step", fill=False, label=r"true value")
	#plt.hist(plot_reco_lepton_energies, np.linspace(0,400,20), histtype="step", fill=False, label=r"reco value")
	#plt.legend()
	#plt.savefig("out_reco_vs_truth_lepton_energy.png")
	#plt.clf()

	
	plt.xlabel(r"Hadronic W-Boson $p_t$ [GeV]")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_whad_pts/1e3, np.linspace(0,400,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_whad_pts/1e3, np.linspace(0,400,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_whad_pt.png")
	plt.clf()
	
	plt.xlabel(r"Hadronic W-Boson $\eta$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_whad_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_whad_etas, np.linspace(-3,3,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_whad_eta.png")
	plt.clf()
	
	plt.xlabel(r"Hadronic W-Boson $\phi$")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_whad_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"true value")
	plt.hist(plot_reco_whad_phis, np.linspace(-4,4,20), histtype="step", fill=False, label=r"reco value")
	plt.legend()
	plt.savefig("out_reco_vs_truth_whad_phi.png")
	plt.clf()

	# energies not supported by truth objects by default	
	#plt.xlabel(r"Hadronic W-Boson $E$ [GeV]")
	#plt.ylabel(r"Number of Events")
	#plt.hist(plot_true_whad_energies, np.linspace(0,400,20), histtype="step", fill=False, label=r"true value")
	#plt.hist(plot_reco_whad_energies, np.linspace(0,400,20), histtype="step", fill=False, label=r"reco value")
	#plt.legend()
	#plt.savefig("out_reco_vs_truth_whad_energy.png")
	#plt.clf()



def val_weights(root_dictionary):
	# access
	weights_on_true = root_dictionary["NW_on_true_weights"]
	weights_on_reco = root_dictionary["NW_on_reco_weights"]
	true_abs_lepton_pdgids = root_dictionary["event_abs_lepton_pdgids"]


	# prepare
	plot_weights_on_true = np.array([])
	plot_weights_on_reco = np.array([])

	plot_weights_on_reco_electron = np.array([])
	plot_weights_on_reco_muon = np.array([])
	plot_weights_on_reco_tau = np.array([])
	plot_weights_on_reco_invalid = np.array([])


	# fill arrays 
	for (weight_on_true, weight_on_reco, lepton_pdgid) in zip(weights_on_true, weights_on_reco, true_abs_lepton_pdgids):
		plot_weights_on_true = np.append(plot_weights_on_true, [weight_on_true])
		plot_weights_on_reco = np.append(plot_weights_on_reco, [weight_on_reco])
		
		# check reco distribution by lepton flavour
		if lepton_pdgid==11:
			plot_weights_on_reco_electron = np.append(plot_weights_on_reco_electron, [weight_on_reco])
		elif lepton_pdgid==13:
			plot_weights_on_reco_muon = np.append(plot_weights_on_reco_muon, [weight_on_reco])
		elif lepton_pdgid==15:
			plot_weights_on_reco_tau = np.append(plot_weights_on_reco_tau, [weight_on_reco])
		else:
			plot_weights_on_reco_invalid = np.append(plot_weights_on_reco_invalid, [weight_on_reco])
	

	# plot
	plt.xlabel(r"NW weight (only NW using reco)")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_weights_on_reco_electron, np.linspace(-1,1,20), histtype="step", fill=False, label=f"electron ({len(plot_weights_on_reco_electron)} events)")
	plt.hist(plot_weights_on_reco_muon, np.linspace(-1,1,20), histtype="step", fill=False, label=f"muon ({len(plot_weights_on_reco_muon)} events)")
	plt.hist(plot_weights_on_reco_tau, np.linspace(-1,1,20), histtype="step", fill=False, label=f"tau ({len(plot_weights_on_reco_tau)} events)")
	plt.hist(plot_weights_on_reco_invalid, np.linspace(-1,1,20), histtype="step", fill=False, label=f"invalid ({len(plot_weights_on_reco_invalid)} events)")
	plt.legend()
	plt.xlim([-1,1])
	plt.savefig("out_weights_reco_flavour.png")
	plt.clf()



def val_wlep_mass(root_dictionary):
	# access
	true_wlep_lvecs = root_dictionary["true_wlep_lvecs"]
	wlep_masses_on_true= root_dictionary["NW_on_true_wlep_masses"]
	wlep_masses_on_reco = root_dictionary["NW_on_reco_wlep_masses"]
	
	true_abs_lepton_pdgids = root_dictionary["event_abs_lepton_pdgids"]
	

	# preprocess informaton 
	plot_true_wlep_masses = np.array([])
	plot_wlep_masses_on_true = np.array([])
	plot_wlep_masses_on_reco = np.array([])
	
	plot_wlep_masses_on_reco_electron = np.array([])
	plot_wlep_masses_on_reco_muon = np.array([])
	plot_wlep_masses_on_reco_tau = np.array([])
	plot_wlep_masses_on_reco_invalid = np.array([])

	for (true_wlep_lvec, wlep_mass_on_true, wlep_mass_on_reco, lepton_pdgid) in zip(true_wlep_lvecs, wlep_masses_on_true, wlep_masses_on_reco, true_abs_lepton_pdgids):
		true_wlep_mass = true_wlep_lvec.fCoordinates.tolist().get("fCoordinates.fM")

		plot_true_wlep_masses = np.append(plot_true_wlep_masses, [true_wlep_mass])
		plot_wlep_masses_on_true = np.append(plot_wlep_masses_on_true, [wlep_mass_on_true])
		plot_wlep_masses_on_reco = np.append(plot_wlep_masses_on_reco, [wlep_mass_on_reco])
			
		# check reco distribution by lepton flavour
		if lepton_pdgid==11:
			plot_wlep_masses_on_reco_electron = np.append(plot_wlep_masses_on_reco_electron, [wlep_mass_on_reco])
		elif lepton_pdgid==13:
			plot_wlep_masses_on_reco_muon = np.append(plot_wlep_masses_on_reco_muon, [wlep_mass_on_reco])
		elif lepton_pdgid==15:
			plot_wlep_masses_on_reco_tau = np.append(plot_wlep_masses_on_reco_tau, [wlep_mass_on_reco])
		else:
			plot_wlep_masses_on_reco_invalid = np.append(plot_wlep_masses_on_reco_invalid, [wlep_mass_on_reco])


	# plot true vs NW on truth vs NW on reco
	plt.xlabel(r"Mass leptonic $W$-boson [GeV]")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_true_wlep_masses/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=r"true mass")
	plt.hist(plot_wlep_masses_on_true/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=r"NW on truth")
	plt.hist(plot_wlep_masses_on_reco/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=r"NW on reco")
	plt.legend()
	plt.xlim([0,50])
	plt.savefig("out_wlep_mass.png")
	plt.clf()

	# plot NW on reco for lepton flavours
	plt.xlabel(r"Mass leptonic $W$-boson [GeV] (only NW on reco)")
	plt.ylabel(r"Number of Events")
	plt.hist(plot_wlep_masses_on_reco_electron/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=f"electron flavour ({len(plot_wlep_masses_on_reco_electron)} events)")
	plt.hist(plot_wlep_masses_on_reco_muon/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=f"muon flavour ({len(plot_wlep_masses_on_reco_muon)} events)")
	plt.hist(plot_wlep_masses_on_reco_tau/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=f"tau flavour ({len(plot_wlep_masses_on_reco_tau)} events)")
	plt.hist(plot_wlep_masses_on_reco_invalid/1e3, np.linspace(0,50,20), histtype="step", fill=False, label=f"invalid lepton ({len(plot_wlep_masses_on_reco_invalid)} events)")
	plt.legend()
	plt.xlim([0,50])
	plt.savefig("out_wlep_mass_reco_flavour.png")
	plt.clf()
	


if __name__ == '__main__':
	## INPUT 
	# validate arguments
	if(len(sys.argv)<2):
		print("Error: no .root file for validation was given, exiting")
		exit()
	

	# create dictionary
	root_dictionary = create_dictionary(sys.argv[1])
	
	# call methods
	val_weights(root_dictionary)
	val_wlep_mass(root_dictionary)
	val_reco_vs_truth_kinematics(root_dictionary)
	val_reco_vs_truth_lepton(root_dictionary)
	val_reco_vs_truth_lepton_success(root_dictionary)

