#include <string>
#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include <Math/Vector4D.h> // for PtEtaPhiEVector
using namespace ROOT::Math;
#include <chrono> // for time measurements during the program


// ==========  DELTA R MATCHING  ==========
// ========================================
// matches jets with truth level particles using the delta R variable
//
// todo: add information about invariant mass to improve matching (?)



// ==========  CONSTANTS  ==========
// =================================
const bool VERBOSE = true;
const float DELTA_R_THRESHOLD = 0.4;
const std::string INPUT_PATH = "/home/ireas/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/master/output/v1/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};


// ==========  FUNCTION DECLARATION  ==========
// ============================================
// description: fills int-array with best matching indicies of jets for the truth objects, also checks whether reconstruction is successful
void fill_with_best_indicies(
	bool* reconstruction_successful,
	int* indicies,
	int number_of_jets,
	PtEtaPhiEVector* jet_lvectors,
	PtEtaPhiMVector wboson1_decay1_lvector,
	PtEtaPhiMVector wboson1_decay2_lvector,
	PtEtaPhiMVector wboson2_decay1_lvector,
	PtEtaPhiMVector wboson2_decay2_lvector,
	PtEtaPhiMVector bquark1_lvector,
	PtEtaPhiMVector bquark2_lvector
);

// description: finds best matching jet for one truth object via deltaR matching, ignores previously matched jets
int find_best_match(
	std::unordered_set<int>* unavailable_indicies,
	int number_of_jets,
	PtEtaPhiEVector* jet_lvectors,
	PtEtaPhiMVector truth_object_lvector
);

// description: calculates deltaR of two objects, ROOT::MATH::LorentzVector does not have deltaR() function
float calc_delta_R(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);


// ==========  MAIN  ==========
// ============================
int main(){
	// ========== CHAIN EVENTS
	// =======================
	TChain reco_chain("reco");
	TChain truth_chain("truth");
	for(int i=0; i<(int)(sizeof(INPUT_FILE_NAMES)/sizeof(INPUT_FILE_NAMES[0])); i++){
		truth_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
		reco_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
	}
	

	// ==========  SETUP
	// =================

	// global variables
	ULong64_t current_event_number = 0;
	std::unordered_map<ULong64_t, int> truth_event_indicies; // unordered set for fast finding matching eventNumbers
	
	
	// reco: jets
	std::vector<Float_t>* jet_pt = nullptr;
	std::vector<Float_t>* jet_eta = nullptr;
	std::vector<Float_t>* jet_phi = nullptr;
	std::vector<Float_t>* jet_e = nullptr;

	reco_chain.SetBranchAddress("jet_pt_NOSYS", &jet_pt);
	reco_chain.SetBranchAddress("jet_eta", &jet_eta);
	reco_chain.SetBranchAddress("jet_phi", &jet_phi);
	reco_chain.SetBranchAddress("jet_e_NOSYS", &jet_e);
	
	// truth: wboson1 - decay1
	PtEtaPhiMVector wboson1_decay1_lvector(0,0,0,0);
	Float_t wboson1_decay1_pt = 0;
	Float_t wboson1_decay1_eta = 0;
	Float_t wboson1_decay1_phi = 0;
	Float_t wboson1_decay1_m = 0;
	Int_t wboson1_decay1_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_pt", &wboson1_decay1_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_eta", &wboson1_decay1_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_phi", &wboson1_decay1_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_m", &wboson1_decay1_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_t_pdgId", &wboson1_decay1_id);
	
	// truth: wboson1 - decay2
	PtEtaPhiMVector wboson1_decay2_lvector(0,0,0,0);
	Float_t wboson1_decay2_pt = 0;
	Float_t wboson1_decay2_eta = 0;
	Float_t wboson1_decay2_phi = 0;
	Float_t wboson1_decay2_m = 0;
	Int_t wboson1_decay2_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_pt", &wboson1_decay2_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_eta", &wboson1_decay2_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_phi", &wboson1_decay2_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_m", &wboson1_decay2_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_t_pdgId", &wboson1_decay2_id);
	
	// truth: wboson1 - decay1
	PtEtaPhiMVector wboson2_decay1_lvector(0,0,0,0);
	Float_t wboson2_decay1_pt = 0;
	Float_t wboson2_decay1_eta = 0;
	Float_t wboson2_decay1_phi = 0;
	Float_t wboson2_decay1_m = 0;
	Int_t wboson2_decay1_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_pt", &wboson2_decay1_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_eta", &wboson2_decay1_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_phi", &wboson2_decay1_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_m", &wboson2_decay1_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_pdgId", &wboson2_decay1_id);
	
	// truth: wboson1 decay2
	PtEtaPhiMVector wboson2_decay2_lvector(0,0,0,0);
	Float_t wboson2_decay2_pt = 0;
	Float_t wboson2_decay2_eta = 0;
	Float_t wboson2_decay2_phi = 0;
	Float_t wboson2_decay2_m = 0;
	Int_t wboson2_decay2_id = 0;

	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_pt", &wboson2_decay2_pt);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_eta", &wboson2_decay2_eta);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_phi", &wboson2_decay2_phi);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_m", &wboson2_decay2_m);
	truth_chain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_pdgId", &wboson2_decay2_id);
	
	// truth: bquark1
	PtEtaPhiMVector bquark1_lvector(0,0,0,0);
	Float_t bquark1_pt = 0;
	Float_t bquark1_eta = 0;
	Float_t bquark1_phi = 0;
	Float_t bquark1_m = 0;

	truth_chain.SetBranchAddress("Tth_MC_b_from_t_pt", &bquark1_pt);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_eta", &bquark1_eta);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_phi", &bquark1_phi);
	truth_chain.SetBranchAddress("Tth_MC_b_from_t_m", &bquark1_m);
	
	// truth: bquark2
	PtEtaPhiMVector bquark2_lvector(0,0,0,0);
	Float_t bquark2_pt = 0;
	Float_t bquark2_eta = 0;
	Float_t bquark2_phi = 0;
	Float_t bquark2_m = 0;

	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_pt", &bquark2_pt);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_eta", &bquark2_eta);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_phi", &bquark2_phi);
	truth_chain.SetBranchAddress("Tth_MC_b_from_tbar_m", &bquark2_m);

	// jet indicies for matching during training
	TTree tree("jet_indicies", "jet indicies for wbosons and bquarks");
	
	bool reconstruction_successful = false;
	int wboson1_decay1_jet_index = -1;
	int wboson1_decay2_jet_index = -1;
	int wboson2_decay1_jet_index = -1;
	int wboson2_decay2_jet_index = -1;
	int bquark1_jet_index = -1;
	int bquark2_jet_index = -1;
	
	tree.Branch("reconstruction_successful", &reconstruction_successful);
	tree.Branch("wboson1_decay1_jet_index", &wboson1_decay1_jet_index);
	tree.Branch("wboson1_decay2_jet_index", &wboson1_decay2_jet_index);
	tree.Branch("wboson2_decay1_jet_index", &wboson2_decay1_jet_index);
	tree.Branch("wboson2_decay2_jet_index", &wboson2_decay2_jet_index);
	tree.Branch("bquark1_jet_index", &bquark1_jet_index);
	tree.Branch("bquark2_jet_index", &bquark2_jet_index);

	// ==========  EVENT LOOPS
	// =======================

	// truth loop: save all entries in hashmap for fast and labeled access
	std::chrono::steady_clock::time_point timestamp_begin_truth_chain = std::chrono::steady_clock::now();
	truth_chain.SetBranchAddress("eventNumber", &current_event_number);
	for(int i=0; i<truth_chain.GetEntries(); i++){
		truth_chain.GetEntry(i);	
		if(truth_event_indicies.count(current_event_number)!=0){
			std::cout << "Warning: multiple events with same event number found (" << current_event_number << "), skipping" << std::endl;
			continue;
		}

		truth_event_indicies.insert({current_event_number,i}); // map eventNumber to index for reco loop
	}
	std::chrono::steady_clock::time_point timestamp_end_truth_chain = std::chrono::steady_clock::now();
	if(VERBOSE)
		std::cout << "> time needed for truth chain: " << std::chrono::duration_cast<std::chrono::seconds>(timestamp_end_truth_chain - timestamp_begin_truth_chain).count() << "s" << std::endl;

	
	// reco loop: match jets to objects via delta R
	std::chrono::steady_clock::time_point timestamp_begin_reco_chain = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point timestamp_last_checkpoint = std::chrono::steady_clock::now();
	float percentage = 0.0;
	float percentage_step = 1.0;
	float percentage_target = 0.1;
	if(reco_chain.GetEntries()!=0)
		percentage_step = 1.0/reco_chain.GetEntries();

	reco_chain.SetBranchAddress("eventNumber", &current_event_number);
	for(int i=0; i<reco_chain.GetEntries(); i++){
		// print estimated time
		percentage+= percentage_step;
		if(percentage<1.0 && percentage>=percentage_target){
			std::chrono::steady_clock::time_point timestamp_checkpoint = std::chrono::steady_clock::now();
			int duration = std::chrono::duration_cast<std::chrono::seconds>(timestamp_checkpoint- timestamp_last_checkpoint).count();			
			if(VERBOSE)
				std::cout << "> time for " << static_cast<int>(100*percentage_target) << "%: " << duration << "s (est. remaining time: " << static_cast<int>(duration*10*(1.0-percentage_target)) << "s)"<< std::endl;
			timestamp_last_checkpoint = timestamp_checkpoint;
			percentage_target+= 0.1;
		}
		
		// get entries
		reco_chain.GetEntry(i);
		if(truth_event_indicies.count(current_event_number)==0){
			std::cout << "Warning: event number not found (" << current_event_number << "), skipping" << std::endl;
			continue;
		}
	
		// get corresponding truth event from map
		truth_chain.GetEntry(truth_event_indicies[current_event_number]);
		truth_event_indicies.erase(current_event_number); // remove event from map (is not really important i guess)


		// put reco jets into lorentzvectors	
		int number_of_jets = jet_pt->size();
		PtEtaPhiEVector jet_lvectors[number_of_jets];
		for(int j=0; j<number_of_jets; j++){
			jet_lvectors[j].SetCoordinates( (*jet_pt)[j], (*jet_eta)[j], (*jet_phi)[j], (*jet_e)[j] );
		}
		

		// put truth objects into lorentzvectors
		wboson1_decay1_lvector.SetCoordinates( wboson1_decay1_pt, wboson1_decay1_eta, wboson1_decay1_phi, wboson1_decay1_m );
		wboson1_decay2_lvector.SetCoordinates( wboson1_decay2_pt, wboson1_decay2_eta, wboson1_decay2_phi, wboson1_decay2_m );
		wboson2_decay1_lvector.SetCoordinates( wboson2_decay1_pt, wboson2_decay1_eta, wboson2_decay1_phi, wboson2_decay1_m );
		wboson2_decay2_lvector.SetCoordinates( wboson2_decay2_pt, wboson2_decay2_eta, wboson2_decay2_phi, wboson2_decay2_m );
		bquark1_lvector.SetCoordinates( bquark1_pt, bquark1_eta, bquark1_phi, bquark1_m );
		bquark2_lvector.SetCoordinates( bquark2_pt, bquark2_eta, bquark2_phi, bquark2_m );
		

		// get best matching jet indicies
		int best_indicies[6];
		
		fill_with_best_indicies(
			&reconstruction_successful,
			best_indicies, 
			number_of_jets,
			jet_lvectors,
			wboson1_decay1_lvector,
			wboson1_decay2_lvector,
			wboson2_decay1_lvector,
			wboson2_decay2_lvector,
			bquark1_lvector,
			bquark2_lvector
		);

		wboson1_decay1_jet_index = best_indicies[0];
		wboson1_decay2_jet_index = best_indicies[1];
		wboson2_decay1_jet_index = best_indicies[2];
		wboson2_decay2_jet_index = best_indicies[3];
		bquark1_jet_index = best_indicies[4];
		bquark2_jet_index = best_indicies[5];

		// fill data to tree
		tree.Fill();
	}

	std::chrono::steady_clock::time_point timestamp_end_reco_chain = std::chrono::steady_clock::now();
	if(VERBOSE){
		std::cout << "> time needed for reco chain: " << std::chrono::duration_cast<std::chrono::seconds>(timestamp_end_reco_chain - timestamp_begin_reco_chain).count() << "s" << std::endl;
		std::cout << "> time needed to complete: " << std::chrono::duration_cast<std::chrono::seconds>(timestamp_end_reco_chain - timestamp_begin_truth_chain).count() << "s" << std::endl;
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


void fill_with_best_indicies(
	bool* reconstruction_successful,
	int* indicies,
	int number_of_jets,
	PtEtaPhiEVector* jet_lvectors,
	PtEtaPhiMVector wboson1_decay1_lvector,
	PtEtaPhiMVector wboson1_decay2_lvector,
	PtEtaPhiMVector wboson2_decay1_lvector,
	PtEtaPhiMVector wboson2_decay2_lvector,
	PtEtaPhiMVector bquark1_lvector,
	PtEtaPhiMVector bquark2_lvector
)
{		
	// set of unavalable indicies (jets that have already been matched)
	std::unordered_set<int> unavailable_indicies;


	// find best match for each truth object
	indicies[0] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, wboson1_decay1_lvector);
	indicies[1] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, wboson1_decay2_lvector);
	indicies[2] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, wboson2_decay1_lvector);
	indicies[3] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, wboson2_decay2_lvector);
	indicies[4] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, bquark1_lvector);
	indicies[5] = find_best_match(&unavailable_indicies, number_of_jets, jet_lvectors, bquark2_lvector);


	// check whether matching was successful or not
	if(number_of_jets<6) 
		*reconstruction_successful = false; // unsuccessful, not enough jets for all truth objects
	
	for(int i=0; i<6; i++){ 
		if(indicies[i]==-1)
			*reconstruction_successful = false; // unsuccessful, not every truth object has a matching jet index
	}

	*reconstruction_successful = true; // successful matching
}

int find_best_match(
	std::unordered_set<int>* unavailable_indicies,
	int number_of_jets,
	PtEtaPhiEVector* jet_lvectors,
	PtEtaPhiMVector truth_object_lvector
)
{
	float best_delta_R = 999;
	int best_delta_R_index = -1;

	for(int current_index=0; current_index<number_of_jets; current_index++){
		if(unavailable_indicies->count(current_index)>0) // if unavailable, skip
			continue;

		float delta_R = calc_delta_R(truth_object_lvector.Eta(), truth_object_lvector.Phi(), jet_lvectors[current_index].Eta(), jet_lvectors[current_index].Phi() );
		if( (delta_R<best_delta_R) && (delta_R<=DELTA_R_THRESHOLD) ){ // must be better than previous match and below threshold
			best_delta_R = delta_R;
			best_delta_R_index = current_index;
		}
	}

	unavailable_indicies->insert(best_delta_R_index); // index was matched, add to unavailable indicies for next matching proccess
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
				// WARNing FROM THIS IS TRANSVERSE MASS NOT INVARIANT MASS!!!
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
