#include <string>
#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"


// ==========  DELTA R MATCHING  ==========
// ========================================
// matches jets with truth level particles using the delta R variable
//
// todo: add timestep/completion prints
// todo: add information about invariant mass to improve matching (?)


// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 0.4;
const std::string INPUT_PATH = "/home/ireas/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/master/output/v1/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
bool fill_with_best_indicies(
	int* indicies,
	std::vector<Float_t> jet_eta,
	std::vector<Float_t> jet_phi,
	Float_t wboson1_decay1_eta,
	Float_t wboson1_decay1_phi,
	Float_t wboson1_decay2_eta,
	Float_t wboson1_decay2_phi,
	Float_t wboson2_decay1_eta,
	Float_t wboson2_decay1_phi,
	Float_t wboson2_decay2_eta,
	Float_t wboson2_decay2_phi,
	Float_t bquark1_eta,
	Float_t bquark1_phi,
	Float_t bquark2_eta,
	Float_t bquark2_phi
);

int find_best_match(
	std::unordered_set<int>* unavailable_indicies,
	std::vector<Float_t> jet_eta,
	std::vector<Float_t> jet_phi,
	Float_t truth_eta,
	Float_t truth_phi
);

float calc_delta_R(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);


// ==========  MAIN  ==========
// ============================
int main(){
	// ========== FILES  
	// ================
	

	// collect all truth/reco trees
	TChain reco_chain("reco");
	TChain truth_chain("truth");
	for(int i=0; i<(int)(sizeof(INPUT_FILE_NAMES)/sizeof(INPUT_FILE_NAMES[0])); i++){
		truth_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
		reco_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
	}
	

	// ==========  SETUP
	// =================

	// global
	ULong64_t current_event_number = 0;
	std::unordered_map<ULong64_t, int> truth_event_indicies; // unordered set for fast finding matching eventNumbers
	
	
	// reco: jets
	std::vector<Float_t>* jet_e = nullptr;
	std::vector<Float_t>* jet_pt = nullptr;
	std::vector<Float_t>* jet_eta = nullptr;
	std::vector<Float_t>* jet_phi = nullptr;

	reco_chain.SetBranchAddress("jet_e_NOSYS", &jet_e);
	reco_chain.SetBranchAddress("jet_pt_NOSYS", &jet_pt);
	reco_chain.SetBranchAddress("jet_eta", &jet_eta);
	reco_chain.SetBranchAddress("jet_phi", &jet_phi);
	
	// truth: wboson1 - decay1
	Float_t wboson1_decay1_m = 0;
	Float_t wboson1_decay1_pt = 0;
	Float_t wboson1_decay1_phi = 0;
	Float_t wboson1_decay1_eta = 0;
	Int_t wboson1_decay1_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_m", &wboson1_decay1_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_pt", &wboson1_decay1_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_eta", &wboson1_decay1_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_phi", &wboson1_decay1_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_pdgId", &wboson1_decay1_id);
	
	// truth: wboson1 - decay2
	Float_t wboson1_decay2_m = 0;
	Float_t wboson1_decay2_pt = 0;
	Float_t wboson1_decay2_phi = 0;
	Float_t wboson1_decay2_eta = 0;
	Int_t wboson1_decay2_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_m", &wboson1_decay2_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_pt", &wboson1_decay2_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_eta", &wboson1_decay2_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_phi", &wboson1_decay2_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_pdgId", &wboson1_decay2_id);
	
	// truth: wboson1 - decay1
	Float_t wboson2_decay1_m = 0;
	Float_t wboson2_decay1_pt = 0;
	Float_t wboson2_decay1_phi = 0;
	Float_t wboson2_decay1_eta = 0;
	Int_t wboson2_decay1_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_m", &wboson2_decay1_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_pt", &wboson2_decay1_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_eta", &wboson2_decay1_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_phi", &wboson2_decay1_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_pdgId", &wboson2_decay1_id);
	
	// truth: wboson1 decay2
	Float_t wboson2_decay2_m = 0;
	Float_t wboson2_decay2_pt = 0;
	Float_t wboson2_decay2_phi = 0;
	Float_t wboson2_decay2_eta = 0;
	Int_t wboson2_decay2_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_m", &wboson2_decay2_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_pt", &wboson2_decay2_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_eta", &wboson2_decay2_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_phi", &wboson2_decay2_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_pdgId", &wboson2_decay2_id);
	
	// truth: bquark1
	Float_t bquark1_m = 0;
	Float_t bquark1_pt = 0;
	Float_t bquark1_phi = 0;
	Float_t bquark1_eta = 0;

	truth_chain.SetBranchAddress("Tth_MC_b_from_t_m", &bquark1_m);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_pt", &bquark1_pt);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_eta", &bquark1_eta);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_phi", &bquark1_phi);
	
	// truth: bquark2
	Float_t bquark2_m = 0;
	Float_t bquark2_pt = 0;
	Float_t bquark2_phi = 0;
	Float_t bquark2_eta = 0;

	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_m", &bquark2_m);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_pt", &bquark2_pt);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_eta", &bquark2_eta);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_phi", &bquark2_phi);
	

	// initiate histograms (name, title, bins, min_val, max_val) to save
	TH1F hist_deltaMass_Wbosons = *(new TH1F("Delta Mass W-Boson", "Delta Mass (W-Boson) in MeV", 100, -50e3, 50e3));
	TH1F hist_deltaMass_bquark = *(new TH1F("Delta Mass b-Quark", "Delta Mass (B-quark) in MeV", 100, -30e3, 30e3));
	TH1F hist_wboson1_ids = *(new TH1F("WBoson1_decay_ids", "Partle Data Group IDs for Wboson1", 10, 0, 10));
	TH1F hist_wboson2_ids = *(new TH1F("WBoson2_decay_ids", "Partle Data Group IDs for Wboson2", 10, 0, 10));
	TH1F hist_bquark1_deltaR = *(new TH1F("DeltaR_bquark1", "Delta R for bquark1", 12, 0, 6));
	TH1F hist_bquark2_deltaR = *(new TH1F("DeltaR_bquark2", "Delta R for bquark2", 12, 0, 6));


	// jet indicies for matching during training
	TTree tree("jet_indicies", "jet indicies for wbosons and bquarks");
	
	bool reconstruction_was_successful = false;
	int wboson1_decay1_jet_index = -1;
	int wboson1_decay2_jet_index = -1;
	int wboson2_decay1_jet_index = -1;
	int wboson2_decay2_jet_index = -1;
	int bquark1_jet_index = -1;
	int bquark2_jet_index = -1;
	
	tree.Branch("reconstruction_was_successful", &reconstruction_was_successful);
	tree.Branch("wboson1_decay1_jet_index", &wboson1_decay1_jet_index);
	tree.Branch("wboson1_decay2_jet_index", &wboson1_decay2_jet_index);
	tree.Branch("wboson2_decay1_jet_index", &wboson2_decay1_jet_index);
	tree.Branch("wboson2_decay2_jet_index", &wboson2_decay2_jet_index);
	tree.Branch("bquark1_jet_index", &bquark1_jet_index);
	tree.Branch("bquark2_jet_index", &bquark2_jet_index);



	// ==========  EVENT LOOP
	// ======================

	// truth loop: save all entries in hashmap for fast and labeled access
	truth_chain.SetBranchAddress("eventNumber", &current_event_number);
	for(int i=0; i<truth_chain.GetEntries(); i++){
		truth_chain.GetEntry(i);	
		if(truth_event_indicies.count(current_event_number)!=0){
			std::cout << "Warning: multiple events with same event number found (" << current_event_number << "), skipping" << std::endl;
			continue;
		}

		truth_event_indicies.insert({current_event_number,i}); // map eventNumber to index for reco loop
	}	
	
	// reco loop: match jets to objects via delta R
	reco_chain.SetBranchAddress("eventNumber", &current_event_number);
	for(int i=0; i<reco_chain.GetEntries(); i++){
		reco_chain.GetEntry(i);
		if(truth_event_indicies.count(current_event_number)==0){
			std::cout << "Warning: event number not found (" << current_event_number << "), skipping" << std::endl;
			continue;
		}
		
		truth_chain.GetEntry(truth_event_indicies[current_event_number]); // get corresponding truth event
		truth_event_indicies.erase(current_event_number); // remove event from mapi (is not really important i guess)
		
		// get best matching jet indicies
		int best_indicies[6];
		reconstruction_was_successful = fill_with_best_indicies(
			best_indicies, 
			(*jet_eta), (*jet_phi), 
			wboson1_decay1_eta, wboson1_decay1_phi,
			wboson1_decay2_eta, wboson1_decay2_phi,
			wboson2_decay1_eta, wboson2_decay1_phi,
			wboson2_decay2_eta, wboson2_decay2_phi,
			bquark1_eta, bquark1_phi,
			bquark2_eta, bquark2_phi
		);
		wboson1_decay1_jet_index = best_indicies[0];
		wboson1_decay2_jet_index = best_indicies[1];
		wboson2_decay1_jet_index = best_indicies[2];
		wboson2_decay2_jet_index = best_indicies[3];
		bquark1_jet_index = best_indicies[4];
		bquark2_jet_index = best_indicies[5];

		// fill indicies
		tree.Fill();
	}


	// ==========  OUTPUT
	// ==================
	TFile* output_file = new TFile( (OUTPUT_PATH+"output_data.root").c_str(), "RECREATE" );
	output_file->cd();
	tree.Write();
	output_file->Close();	
	
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================


bool fill_with_best_indicies(
	int* indicies,
	std::vector<Float_t> jet_eta,
	std::vector<Float_t> jet_phi,
	Float_t wboson1_decay1_eta,
	Float_t wboson1_decay1_phi,
	Float_t wboson1_decay2_eta,
	Float_t wboson1_decay2_phi,
	Float_t wboson2_decay1_eta,
	Float_t wboson2_decay1_phi,
	Float_t wboson2_decay2_eta,
	Float_t wboson2_decay2_phi,
	Float_t bquark1_eta,
	Float_t bquark1_phi,
	Float_t bquark2_eta,
	Float_t bquark2_phi
)
{		
	// define needed variables
	int best_index = -1;
	std::unordered_set<int> unavailable_indicies;

	indicies[0] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, wboson1_decay1_eta, wboson1_decay1_phi);
	indicies[1] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, wboson1_decay2_eta, wboson1_decay2_phi);
	indicies[2] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, wboson2_decay1_eta, wboson2_decay1_phi);	
	indicies[3] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, wboson2_decay2_eta, wboson2_decay2_phi);
	indicies[4] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, bquark1_eta, bquark1_phi);
	indicies[5] = find_best_match(&unavailable_indicies, jet_eta, jet_phi, bquark2_eta, bquark2_phi);

	bool success = true;
	if(jet_eta.size()<6)
		success = false;
	else{
		for(int i=0; i<6; i++){
			if(indicies[i]==-1)
				success = false;
				break;
		}
	}

	return success;
}

int find_best_match(
	std::unordered_set<int>* unavailable_indicies,
	std::vector<Float_t> jet_eta,
	std::vector<Float_t> jet_phi,
	Float_t truth_eta,
	Float_t truth_phi
)
{
	int number_of_jets = jet_eta.size();
	float best_delta_R = 999;
	int best_delta_R_index = -1;

	for(int current_index=0; current_index<number_of_jets; current_index++){
		if((*unavailable_indicies).count(current_index)>0)
			continue;

		float delta_R = calc_delta_R(truth_eta, truth_phi, jet_eta[current_index], jet_phi[current_index]);
		if( (delta_R<best_delta_R) && (delta_R<=DELTA_R_THRESHOLD) ){ // must be below threshold
			best_delta_R = delta_R;
			best_delta_R_index = current_index;
		}
	}

	(*unavailable_indicies).insert(best_delta_R_index);
	return best_delta_R_index;
}

float calc_delta_R(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2){
	return sqrt( (eta2-eta1)*(eta2-eta1) + (phi2-phi1)*(phi2-phi1) );
}



// ==========  OLD CODE  ==========
// ================================
//		// find two jets to match 80GeV (W-Boson) as close as possible (simple algorithm)
//		for(int a=0; a<nJets; a++){
//			for(int b=a+1; b<nJets; b++){
//				float invariantMass = sqrt( ((*jet_e)[a]+(*jet_e)[b])*((*jet_e)[a]+(*jet_e)[b]) - ((*jet_pt)[a]+(*jet_pt)[b])*((*jet_pt)[a]+(*jet_pt)[b]) ); // m_0^2 = M^2 = E^2 - p^2
//				//std::cout << "  (" << a << "," << b << "): "<< invariantMass-80e3 << " " << recoMassW1-80e3 << std::endl;
//					
//				// Wboson1
//				if( abs(invariantMass-80e3) < abs(recoMassW1-80e3) ){
//					indexW1a = a;
//					indexW1b = b;
//					recoMassW2 = recoMassW1;				
//					recoMassW1 = invariantMass;
//				}	
//				// Wboson2
//				else if( (abs(invariantMass-80e3) < abs(recoMassW2-80e3)) && (a!=indexW1a) && (a!=indexW1b) && (b!=indexW1a) && (b!=indexW1b)){ // second best mass with different jets
//					indexW2a = a;
//					indexW2b = b;
//					recoMassW2 = invariantMass;	
//				}	
//			} 
//		}
