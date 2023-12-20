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
// TODO:improve matching function



// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 0.4;

const std::string INPUT_PATH = "/home/ireas/git_repos/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/git_repos/master/output/v1/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};

const plog::Severity LOG_LEVEL = plog::Severity::debug; // set logger output severity filer


enum TRUTH_OBJ_MATCH_VALUES{
	unmatched = 0,
	b_from_t = 1,
	b_from_tbar = 2,
	Wdecay1_from_t = 3,
	Wdecay2_from_t = 4,
	Wdecay1_from_tbar = 6,
	Wdecay2_from_tbar = 7,
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
// calculates delta_R for 2 obj
float delta_R(float eta1, float phi1, float eta2, float phi2);
float delta_R(PtEtaPhiEVector jet, PtEtaPhiMVector obj);


// generate lorentz vector for truth object
PtEtaPhiMVector generate_lorentz_vector_m(Float_t pt, Float_t eta, Float_t phi, Float_t m);
PtEtaPhiEVector generate_lorentz_vector_e(Float_t pt, Float_t eta, Float_t phi, Float_t e);


// generate lorentz vectors for jets
std::vector<PtEtaPhiEVector> generate_lvecs(std::vector<Float_t> obj_pt, std::vector<Float_t> obj_eta, std::vector<Float_t> obj_phi, std::vector<Float_t> obj_e);


// generate jet bit_id_mask
std::vector<int> match_jets_to_objs(
	std::vector<PtEtaPhiEVector> jets, 
	PtEtaPhiMVector obj_b_from_t, 
	PtEtaPhiMVector obj_b_from_tbar,
	PtEtaPhiMVector obj_Wdecay1_from_t,
	PtEtaPhiMVector obj_Wdecay2_from_t,
	PtEtaPhiMVector obj_Wdecay1_from_tbar,
	PtEtaPhiMVector obj_Wdecay2_from_tbar
);
int find_best_jet(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj);


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


	// for speedy testing
	auto loop_manager = data_frame.Range(10);


	// generate jet lorentz vectors
	loop_manager = loop_manager.Define("jet_lvecs", generate_lvecs, {"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"});
	

	// generate truth obj lorentz vectors
	loop_manager = loop_manager.Define(
		"lvec_b_from_t", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_b_from_t_pt", "truth.Tth_MC_b_from_t_eta", "truth.Tth_MC_b_from_t_phi", "truth.Tth_MC_b_from_t_m"}
	);
	loop_manager = loop_manager.Define(
		"lvec_b_from_tbar", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_b_from_tbar_pt", "truth.Tth_MC_b_from_tbar_eta", "truth.Tth_MC_b_from_tbar_phi", "truth.Tth_MC_b_from_tbar_m"}
	);
	loop_manager = loop_manager.Define(
		"lvec_Wdecay1_from_t", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_Wdecay1_from_t_pt", "truth.Tth_MC_Wdecay1_from_t_eta", "truth.Tth_MC_Wdecay1_from_t_phi", "truth.Tth_MC_Wdecay1_from_t_m"}
	);
	loop_manager = loop_manager.Define(
		"lvec_Wdecay2_from_t", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_Wdecay2_from_t_pt", "truth.Tth_MC_Wdecay2_from_t_eta", "truth.Tth_MC_Wdecay2_from_t_phi", "truth.Tth_MC_Wdecay2_from_t_m"}
	);
	loop_manager = loop_manager.Define(
		"lvec_Wdecay1_from_tbar", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_Wdecay1_from_tbar_pt", "truth.Tth_MC_Wdecay1_from_tbar_eta", "truth.Tth_MC_Wdecay1_from_tbar_phi", "truth.Tth_MC_Wdecay1_from_tbar_m"}
	);
	loop_manager = loop_manager.Define(
		"lvec_Wdecay2_from_tbar", 
		generate_lorentz_vector_m, 
		{"truth.Tth_MC_Wdecay2_from_tbar_pt", "truth.Tth_MC_Wdecay2_from_tbar_eta", "truth.Tth_MC_Wdecay2_from_tbar_phi", "truth.Tth_MC_Wdecay2_from_tbar_m"}
	);


	// matching jets and objs
	loop_manager = loop_manager.Define(
		"jet_match_mask", 
		match_jets_to_objs, 
		{"jet_lvecs", "lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar"}
		);


	// accessing stuff
	PLOG_DEBUG << "start: accessing RDataFrame";
	loop_manager.Display("jet_lvecs")->Print();
	loop_manager.Display("jet_match_mask")->Print();
	

	// finishing
	PLOG_DEBUG << "start: writing to disk";
	loop_manager.Snapshot("matched", OUTPUT_PATH+"rdataframes_output.root");
	PLOG_DEBUG << "start: finished successfully";
	
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================
float delta_R(float eta1, float phi1, float eta2, float phi2){ 
	return sqrt( (eta2-eta1)*(eta2-eta1) + (phi2-phi1)*(phi2-phi1) ); 
}


float delta_R(PtEtaPhiEVector jet, PtEtaPhiMVector obj){ 
	return delta_R(jet.Eta(), jet.Phi(), obj.Eta(), obj.Phi());
}



PtEtaPhiMVector generate_lorentz_vector_m(Float_t pt, Float_t eta, Float_t phi, Float_t m){
	PtEtaPhiMVector lvec(pt, eta, phi, m);
	return lvec;
}


PtEtaPhiEVector generate_lorentz_vector_e(Float_t pt, Float_t eta, Float_t phi, Float_t e){
	PtEtaPhiEVector lvec(pt, eta, phi, e);
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





std::vector<int> match_jets_to_objs(
	std::vector<PtEtaPhiEVector> jets, 
	PtEtaPhiMVector obj_b_from_t, 
	PtEtaPhiMVector obj_b_from_tbar,
	PtEtaPhiMVector obj_Wdecay1_from_t,
	PtEtaPhiMVector obj_Wdecay2_from_t,
	PtEtaPhiMVector obj_Wdecay1_from_tbar,
	PtEtaPhiMVector obj_Wdecay2_from_tbar
){
	std::vector<int> jet_match_masks(jets.size());	

	int current_index;
	
	//match b_from_t
	current_index = find_best_jet(jets, obj_b_from_t);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_t;


	//match b_from_tbar
	current_index = find_best_jet(jets, obj_b_from_tbar);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_tbar;
	

	//match Wdecay1_from_t
	current_index = find_best_jet(jets, obj_Wdecay1_from_t);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_t;
	

	//match Wdecay2_from_t
	current_index = find_best_jet(jets, obj_Wdecay2_from_t);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_t;
	

	//match Wdecay1_from_tbar
	current_index = find_best_jet(jets, obj_Wdecay1_from_t);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_tbar;
	

	//match Wdecay2_from_tbar
	current_index = find_best_jet(jets, obj_Wdecay2_from_tbar);
	if(current_index!=-1)
		jet_match_masks[current_index]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_tbar;

	
	return jet_match_masks;
}

int find_best_jet(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj){
	int best_index = -1;
	float best_delta_R = 999;

	for(int current_index=0; current_index<jets.size(); current_index++){
		float current_delta_R = delta_R(jets[current_index], truth_obj);
		if(current_delta_R>DELTA_R_THRESHOLD)
			continue;
		if(current_delta_R>best_delta_R)
			continue;

		best_delta_R = current_delta_R;
		best_index = current_index;
	}

	return best_index;
}


bool check_match(std::vector<int> jet_match_masks, int truth_obj_bit_value){
	for(int i=0; i<jet_match_masks.size(); i++){
		if( (jet_match_masks[i] & 1<<truth_obj_bit_value)!=0 ){
			return true;
		}
	}

	return false;
}


int generate_truth_match_value(){return 3;}
