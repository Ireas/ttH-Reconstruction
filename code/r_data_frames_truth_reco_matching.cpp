#include<string>
#include<iostream>
#include<vector>
#include"TROOT.h"
#include"TFile.h"
#include"TChain.h"
#include"TTree.h"
#include"TH1F.h"
#include"ROOT/RDataFrame.hxx"
#include<Math/Vector4D.h> // for PtEtaPhiEVector
using namespace ROOT;
using namespace ROOT::Math;

#include<chrono> // for time measurements during the program

#include<plog/Log.h> // c++ logger + init headers
#include<plog/Formatters/TxtFormatter.h>
#include<plog/Initializers/ConsoleInitializer.h>



// ==========  DELTA R MATCHING  ==========
// ========================================
// uses ROOTs modern RDataFrames to simplify the event loop and branchaddress allocation
//
// todo: add information about invariant mass to improve matching (?)




// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 0.4;

const std::string INPUT_PATH = "/home/ireas/git_repos/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/git_repos/master/output/v1/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};

const plog::Severity LOG_LEVEL = plog::Severity::debug; // set logger output severity filer


// TODO: implement bit_definitions for truth objects
enum BIT_ID{
	undefined = 0,
	b_from_t = 1,
	b_from_tbar = 2
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
float calc_delta_R(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2){
	return sqrt( (eta2-eta1)*(eta2-eta1) + (phi2-phi1)*(phi2-phi1) );
}


PtEtaPhiMVector generate_lvec_from_m(Float_t obj_pt, Float_t obj_eta, Float_t obj_phi, Float_t obj_m){
	PtEtaPhiMVector lvec(obj_pt, obj_eta, obj_phi, obj_m);
	return lvec;
}


PtEtaPhiEVector generate_lvec(Float_t obj_pt, Float_t obj_eta, Float_t obj_phi, Float_t obj_e){
	PtEtaPhiEVector lvec(obj_pt, obj_eta, obj_phi, obj_e);
	return lvec;
}


std::vector<PtEtaPhiEVector> generate_lvecs(std::vector<Float_t> obj_pt, std::vector<Float_t> obj_eta, std::vector<Float_t> obj_phi, std::vector<Float_t> obj_e){
	std::vector<PtEtaPhiEVector> lorentz_vectors;
	for(int i=0; i<obj_pt.size(); i++){
		PtEtaPhiEVector lorentz_vector(obj_pt[i], obj_eta[i], obj_phi[i], obj_e[i]);
		lorentz_vectors.push_back(lorentz_vector);
	}
	return lorentz_vectors;
}


std::vector<int> generate_bit_id_mask(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj){
	std::vector<int> bit_id_mask;
	
	float best_delta_R = 999;
	int index_best = 0;
	bool found_best_match = false;

	for(int i=0; i<jets.size(); i++){
		bit_id_mask.push_back(0);

		float current_delta_R = calc_delta_R(truth_obj.Eta(), truth_obj.Phi(), jets[i].Eta(), jets[i].Phi());
		if(current_delta_R<DELTA_R_THRESHOLD && current_delta_R<best_delta_R){
			best_delta_R = current_delta_R;
			index_best = 0;
		}
	}

	if(found_best_match){
		bit_id_mask[index_best] = 1;
	}

	return bit_id_mask;
}



// ==========  MAIN  ==========
// ============================
int main(){
	// ==========  SETUP
	// =================
	// setup plog (logger)
	static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
	plog::init(LOG_LEVEL, &consoleAppender);
	

	// setup ROOT
	PLOG_DEBUG << "start: setup ROOT";
	//ROOT::EnableImplicitMT()

	
	// setup TChain
	PLOG_DEBUG << "start: setup TChain";
	TChain reco_chain("reco");
	TChain truth_chain("truth");
	for(int i=0; i<(int)(sizeof(INPUT_FILE_NAMES)/sizeof(INPUT_FILE_NAMES[0])); i++){
		truth_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
		reco_chain.Add( (INPUT_PATH+INPUT_FILE_NAMES[i]).c_str() ); 
	}
	reco_chain.AddFriend(&truth_chain);


	// setup RDataFrame
	PLOG_DEBUG << "start: setup RDataFrame";
	RDataFrame data_frame(reco_chain);
	auto a = data_frame.Range(10);
	auto b = a.Define("jet_lvecs", generate_lvecs, {"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"});
	auto c = b.Define("jet_lvec_decay1", generate_lvec_from_m, {"truth.Tth_MC_Higgs_decay1_pt", "truth.Tth_MC_Higgs_decay1_eta", "truth.Tth_MC_Higgs_decay1_phi", "truth.Tth_MC_Higgs_decay1_m"});
	auto f = c.Define("bit_id_mask", generate_bit_id_mask, {"jet_lvecs", "jet_lvec_decay1"});


	PLOG_DEBUG << "start: accessing RDataFrame";
	f.Display("jet_lvecs")->Print();
	f.Display("jet_lvec_decay1")->Print();
	f.Display("bit_id_mask")->Print();
	
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================
