#include <string>
#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"


int main(){
	// define input/output filepaths
	std::string pathInput = "/home/ireas/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
	std::string pathOutput = "/home/ireas/master/output/v1/";
	
	// put into array for easier access later
	const char* fileNames[1] = {
		"user.ravinab.35392295._000001.output.root"
	};

	// collect all truth/reco trees in array for accessing events
	TChain recoChain("reco");
	TChain truthChain("truth");
	for(int i=0; i<(int)(sizeof(fileNames)/sizeof(fileNames[0])); i++){
		truthChain.Add( (pathInput+fileNames[i]).c_str() ); 
		recoChain.Add( (pathInput+fileNames[i]).c_str() ); 
	}
	

	// initiate variables to fill
	ULong64_t currentEventNumber = 0;
	std::vector<float_t>* jet_e = nullptr;
	std::vector<float_t>* jet_pt = nullptr;
	std::vector<float_t>* jet_phi = nullptr;
	std::vector<float_t>* jet_eta = nullptr;
	
	float_t wBoson1_m = 0;
	float_t wBoson2_m = 0;
	
	float_t bquark1_m = 0;
	float_t bquark1_phi = 0;
	float_t bquark1_eta = 0;
	float_t bquark2_m = 0;
	float_t bquark2_phi = 0;
	float_t bquark2_eta = 0;
	
	Int_t wboson1_decay1_id = 0;
	Int_t wboson1_decay2_id = 0;
	Int_t wboson2_decay1_id = 0;
	Int_t wboson2_decay2_id = 0;

	// set branch addresses to data
	truthChain.SetBranchAddress("eventNumber", &currentEventNumber);
	
	recoChain.SetBranchAddress("jet_e_NOSYS", &jet_e);
	recoChain.SetBranchAddress("jet_pt_NOSYS", &jet_pt);
	recoChain.SetBranchAddress("jet_phi", &jet_phi);
	recoChain.SetBranchAddress("jet_eta", &jet_eta);
	
	truthChain.SetBranchAddress("Tth_MC_W_from_t_m", &wBoson1_m);
	truthChain.SetBranchAddress("Tth_MC_W_from_tbar_m", &wBoson2_m);

	truthChain.SetBranchAddress("Tth_MC_Wdecay1_from_t_pdgId", &wboson1_decay1_id);
	truthChain.SetBranchAddress("Tth_MC_Wdecay2_from_t_pdgId", &wboson1_decay2_id);
	truthChain.SetBranchAddress("Tth_MC_Wdecay1_from_tbar_pdgId", &wboson2_decay1_id);
	truthChain.SetBranchAddress("Tth_MC_Wdecay2_from_tbar_pdgId", &wboson2_decay2_id);
	
	truthChain.SetBranchAddress("Tth_MC_b_from_t_m", &bquark1_m);
	truthChain.SetBranchAddress("Tth_MC_b_from_t_phi", &bquark1_phi);
	truthChain.SetBranchAddress("Tth_MC_b_from_t_eta", &bquark1_eta);
	truthChain.SetBranchAddress("Tth_MC_b_from_tbar_m", &bquark2_m);
	truthChain.SetBranchAddress("Tth_MC_b_from_tbar_phi", &bquark2_phi);
	truthChain.SetBranchAddress("Tth_MC_b_from_tbar_eta", &bquark2_eta);
	
	// unordered set for fast finding matching eventNumbers
	std::unordered_map<ULong64_t, int> truthEventIndicies;

	// initiate histograms (name, title, bins, min_val, max_val) to save
	TH1F hist_deltaMass_Wbosons = *(new TH1F("Delta Mass W-Boson", "Delta Mass (W-Boson) in MeV", 100, -50e3, 50e3));
	TH1F hist_deltaMass_bquark = *(new TH1F("Delta Mass b-Quark", "Delta Mass (B-quark) in MeV", 100, -30e3, 30e3));
	TH1F hist_wboson1_ids = *(new TH1F("WBoson1_decay_ids", "Partle Data Group IDs for Wboson1", 10, 0, 10));
	TH1F hist_wboson2_ids = *(new TH1F("WBoson2_decay_ids", "Partle Data Group IDs for Wboson2", 10, 0, 10));
	TH1F hist_bquark1_deltaR = *(new TH1F("DeltaR_bquark1", "Delta R for bquark1", 12, 0, 6));
	TH1F hist_bquark2_deltaR = *(new TH1F("DeltaR_bquark2", "Delta R for bquark2", 12, 0, 6));


	// ttree
	TTree t1("name", "description");
	UInt_t run_number = 0;
	recoChain.SetBranchAddress("runNumber", &run_number);
	t1.Branch("name", &run_number);

	// event loop: get truth entries
	for(int i=0; i<truthChain.GetEntries(); i++){
		truthChain.GetEntry(i);	
		if(truthEventIndicies.count(currentEventNumber)!=0){
			std::cout << "Warning: multiple detected! " << currentEventNumber << std::endl;
			continue;
		}

		truthEventIndicies.insert({currentEventNumber,i}); // map eventNumber to index for reco loop
	}	
	
	// event loop: reconstruct W-Boson/b-Quark mass
	recoChain.SetBranchAddress("eventNumber", &currentEventNumber);

	for(int i=0; i<recoChain.GetEntries(); i++){
		recoChain.GetEntry(i);
		if(truthEventIndicies.count(currentEventNumber)==0){
			std::cout << "Warning: event number not found! " << currentEventNumber << std::endl;
			continue;
		}
		
		truthChain.GetEntry(truthEventIndicies[currentEventNumber]); // get corresponding truth event
		truthEventIndicies.erase(currentEventNumber); // remove event from map
		
		// ttree
		t1.Fill();
		
		// define indicies to gurantee unique match of jet to object
		int nJets = (*jet_e).size();
		int indexW1a = -1;
		int indexW1b = -1;
		int indexW2a = -1;
		int indexW2b = -1;
		int indexb1 = -1;
		int indexb2 = -1;

		float recoMassW1 = -9e5;
		float recoMassW2 = -9e5;
		float recoMassb1 = -9e5;
		float recoMassb2 = -9e5;

		float reco_b1_eta = 0;
		float reco_b1_phi = 0;
		float reco_b2_eta = 0;
		float reco_b2_phi = 0;
	
		// find two jets to match 80GeV (W-Boson) as close as possible (simple algorithm)
		for(int a=0; a<nJets; a++){
			for(int b=a+1; b<nJets; b++){
				float invariantMass = sqrt( ((*jet_e)[a]+(*jet_e)[b])*((*jet_e)[a]+(*jet_e)[b]) - ((*jet_pt)[a]+(*jet_pt)[b])*((*jet_pt)[a]+(*jet_pt)[b]) ); // m_0^2 = M^2 = E^2 - p^2
				//std::cout << "  (" << a << "," << b << "): "<< invariantMass-80e3 << " " << recoMassW1-80e3 << std::endl;
					
				// Wboson1
				if( abs(invariantMass-80e3) < abs(recoMassW1-80e3) ){
					indexW1a = a;
					indexW1b = b;
					recoMassW2 = recoMassW1;				
					recoMassW1 = invariantMass;
				}	
				// Wboson2
				else if( (abs(invariantMass-80e3) < abs(recoMassW2-80e3)) && (a!=indexW1a) && (a!=indexW1b) && (b!=indexW1a) && (b!=indexW1b)){ // second best mass with different jets
					indexW2a = a;
					indexW2b = b;
					recoMassW2 = invariantMass;	
				}	
			} 
		}
		
		// find one jet to match 4.12Gev (b-Quark) as close as possible
		for(int a=0; a<nJets; a++){
			float invariantMass = sqrt( (*jet_e)[a]*(*jet_e)[a] - (*jet_pt)[a]*(*jet_pt)[a] ); // m_0^2 = M^2 = E^2 - p^2
			
			if( (a==indexW1a) ||  (a==indexW1b) || (a==indexW2a) || (a==indexW2b) )
				continue;
			
			// bquark1
			if( abs(invariantMass-4.12e3) < abs(recoMassb1-4.12e3)  ){
				indexb1 = a;
				recoMassb1 = invariantMass;
				reco_b1_eta = (*jet_eta)[indexb2];
				reco_b1_phi = (*jet_phi)[indexb2];
			}
			// bquark2
			else if( a!=indexb1 && (abs(invariantMass-4.12e3) < abs(recoMassb1-4.12e3)) ){
				indexb2 = a;
				recoMassb2 = invariantMass;
				reco_b2_eta = (*jet_eta)[indexb2];
				reco_b2_phi = (*jet_phi)[indexb2];
			}
		} 

		// safe mass difference between truth and reconstructed (matched) mass
		hist_deltaMass_Wbosons.Fill( recoMassW1-wBoson1_m );
		hist_deltaMass_Wbosons.Fill( recoMassW2-wBoson2_m);
		hist_deltaMass_bquark.Fill( recoMassb1-bquark1_m);
		hist_deltaMass_bquark.Fill( recoMassb2-bquark2_m );
		hist_wboson1_ids.Fill( wboson1_decay1_id );
		hist_wboson1_ids.Fill( wboson1_decay2_id );
		hist_wboson2_ids.Fill( wboson2_decay1_id );
		hist_wboson2_ids.Fill( wboson2_decay2_id );
		hist_bquark1_deltaR.Fill( sqrt( (reco_b1_eta-bquark1_eta)*(reco_b1_eta-bquark1_eta) + (reco_b1_phi-bquark1_phi)*(reco_b1_phi-bquark1_phi) ) );
		hist_bquark2_deltaR.Fill( sqrt( (reco_b2_eta-bquark2_eta)*(reco_b2_eta-bquark2_eta) + (reco_b2_phi-bquark2_phi)*(reco_b2_phi-bquark2_phi) ) );
//		std::cout << "truth: " << wBoson1Mass*1e-3 << " | " << wBoson2Mass*1e-3 << std::endl;
//		std::cout << "reco:  " << recoMassW1*1e-3 << " | " << recoMassW2*1e-3 << std::endl;
//		std::cout << nJets << ": (" << indexW1a << "," << indexW1b << ") / (" << indexW2a << "," << indexW2b << ")" << std::endl << std::endl;
	}


	// print histogram into output_data.root file
	TFile* fileOutput = new TFile( (pathOutput+"output_data.root").c_str(), "RECREATE" );
	fileOutput->cd();
	hist_deltaMass_Wbosons.Write("DeltaMass Wbosons");
	hist_deltaMass_bquark.Write("DeltaMass bquark");
	hist_wboson1_ids.Write("Wboson1 decay ids");
	hist_wboson2_ids.Write("Wboson2 decay ids");
	hist_bquark1_deltaR.Write("DeltaR bquark1");
	hist_bquark2_deltaR.Write("DeltaR bquark2");
	t1.Write();
	fileOutput->Close();	
	
	return 0;
}
