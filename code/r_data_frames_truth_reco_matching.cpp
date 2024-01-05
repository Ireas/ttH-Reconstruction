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
const float DELTA_R_THRESHOLD = 1;

const std::string INPUT_PATH = "/home/ireas/git_repos/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/git_repos/master/output/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};

const plog::Severity LOG_LEVEL = plog::Severity::debug; // set logger output severity filer

// bit-shift amounts for matching
enum TRUTH_OBJ_MATCH_VALUES{
	b_from_t = 0,
	b_from_tbar = 1,
	Wdecay1_from_t = 2,
	Wdecay2_from_t = 3,
	Wdecay1_from_tbar = 4,
	Wdecay2_from_tbar = 5,
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

int find_best_jet_index(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj);
std::vector<int> find_all_jet_indicies_within_threshold(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj);


// combine jets to parent object
PtEtaPhiEVector combine_tbar(std::vector<PtEtaPhiEVector> jets, std::vector<int> jet_match_masks);
PtEtaPhiEVector combine_t(std::vector<PtEtaPhiEVector> jets, std::vector<int> jet_match_masks);


float mass_from_lvec_in_gev(PtEtaPhiEVector jet);



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
	//ROOT::EnableImplicitMT() // commented out for using fixxed range

	
	// setup TChain
	PLOG_DEBUG << "start: setup TChain";
	TChain reco_chain("reco");
	TChain truth_chain("truth");

	for (auto input_file_name : INPUT_FILE_NAMES) {
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);

		truth_chain.Add(file_path.c_str());
		reco_chain.Add(file_path.c_str());
	}

	// index TruthChain to kick out un-matched events - Then declare friends
	truth_chain.BuildIndex("mcChannelNumber", "eventNumber");  // Just for security, use DSID too
	reco_chain.AddFriend(&truth_chain);

	// setup RDataFrame
	PLOG_DEBUG << "start: setup RDataFrame";
	auto data_frame = RDataFrame(reco_chain);


	// for speedy testing
	auto loop_manager = data_frame.Range(1000);

	// First check that event numbers agree
	auto nTotalEvents = loop_manager.Count();
	auto nMismatchedEvents = loop_manager.Filter(
			"mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber"
	).Count();

	if (nMismatchedEvents.GetValue() > 0) {
		PLOG_ERROR << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue()
		           << " mismatched events!";
		exit(1);
	}

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


	// matching jets and objs (fixed object order)
	loop_manager = loop_manager.Define(
		"jet_match_mask", 
		match_jets_to_objs, 
		{"jet_lvecs", "lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar"}
	);

	
	// combine jets to top and anti-top quarks
	loop_manager = loop_manager.Define(
		"t_combined",
		combine_t,
		{"jet_lvecs", "jet_match_mask"}
	);
	
	loop_manager = loop_manager.Define(
		"tbar_combined",
		combine_tbar,
		{"jet_lvecs", "jet_match_mask"}
	);


	// get mass from top and anti_top quarks
	loop_manager = loop_manager.Define(
		"t_combined_m_GeV",
		mass_from_lvec_in_gev,
		{"t_combined"}
	);
	
	loop_manager = loop_manager.Define(
		"tbar_combined_m_GeV",
		mass_from_lvec_in_gev,
		{"tbar_combined"}
	);




	// accessing and printing information
	PLOG_DEBUG << "start: accessing RDataFrame";
	loop_manager.Display("jet_match_mask")->Print();
	loop_manager.Display("t_combined")->Print();
	loop_manager.Display("tbar_combined")->Print();
	

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
	std::vector<int> jet_indicies_within_threshold;
	
	
	//match b_from_t
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_b_from_t);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_t;
		}
	}
	


	//match b_from_tbar
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_b_from_tbar);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_tbar;
		}
	}
	

	//match Wdecay1_from_t
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_Wdecay1_from_t);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_t;
		}
	}
	

	//match Wdecay2_from_t
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_Wdecay2_from_t);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_t;
		}
	}
	

	//match Wdecay1_from_tbar
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_Wdecay1_from_t);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_tbar;
		}
	}
	

	//match Wdecay2_from_tbar
	jet_indicies_within_threshold = find_all_jet_indicies_within_threshold(jets, obj_Wdecay2_from_tbar);
	if(!jet_indicies_within_threshold.empty()){
		for(int i=0; i<jet_indicies_within_threshold.size(); i++){
			jet_match_masks[jet_indicies_within_threshold[i]]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_tbar;
		}
	}

	
	return jet_match_masks;
}


int find_best_jet_index(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj){
	int best_jet_index = -1;
	float best_delta_R = 999;

	for(int i=0; i<jets.size(); i++){
		float current_delta_R = delta_R(jets[i], truth_obj);
		if(current_delta_R<=DELTA_R_THRESHOLD && current_delta_R<best_delta_R){
			best_delta_R = current_delta_R;
			best_jet_index = i;
		}
	}

	return best_jet_index;
}



std::vector<int> find_all_jet_indicies_within_threshold(std::vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj){
	std::vector<int> all_jet_indicies;

	for(int i=0; i<jets.size(); i++){
		float current_delta_R = delta_R(jets[i], truth_obj);
		if(current_delta_R<=DELTA_R_THRESHOLD){
			all_jet_indicies.push_back(i);
		}
	}

	return all_jet_indicies;
}



int find_match(std::vector<int> jet_match_masks, int truth_obj_bit_value, std::vector<int> ignored_match_mask_indicies = {}){
	for(int i=0; i<jet_match_masks.size(); i++){
		if( !ignored_match_mask_indicies.empty() ){
			bool skip_jet_mask = false;
			
			for(int j=0; j<ignored_match_mask_indicies.size(); j++){
				if(ignored_match_mask_indicies[j]==i){
					skip_jet_mask = true;
					break;
				}
			}
			
			if(skip_jet_mask)
				continue;
		}

	

		if( (jet_match_masks[i] & 1<<truth_obj_bit_value)!=0 ){
			return i;
		}
	}

	return -1;
}


PtEtaPhiEVector combine_t(std::vector<PtEtaPhiEVector> jets, std::vector<int> jet_match_masks){
	if(jets.size()!=jet_match_masks.size()){
		PLOG_ERROR << "ERROR: jet size does not match jet_match_masks size";
		return PtEtaPhiEVector();
	}
	

	std::vector<int> matched_indicies;

	int index_b = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::b_from_t, matched_indicies);
	matched_indicies.push_back(index_b);
	if(index_b==-1)
		return PtEtaPhiEVector();
	
	int index_Wdecay1 = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_t, matched_indicies);
	matched_indicies.push_back(index_Wdecay1);
	if(index_Wdecay1==-1)
		return PtEtaPhiEVector();
	
	int index_Wdecay2 = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_t, matched_indicies);
	matched_indicies.push_back(index_Wdecay2);
	if(index_Wdecay2==-1)
		return PtEtaPhiEVector();


	return jets[index_b] + jets[index_Wdecay1] + jets[index_Wdecay2];
}


PtEtaPhiEVector combine_tbar(std::vector<PtEtaPhiEVector> jets, std::vector<int> jet_match_masks){
	if(jets.size()!=jet_match_masks.size()){
		PLOG_ERROR << "ERROR: jet size does not match jet_match_masks size";
		return PtEtaPhiEVector();
	}
		

	std::vector<int> matched_indicies;

	int index_b = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::b_from_tbar, matched_indicies);
	matched_indicies.push_back(index_b);
	if(index_b==-1)
		return PtEtaPhiEVector();
	
	int index_Wdecay1 = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_tbar, matched_indicies);
	matched_indicies.push_back(index_Wdecay1);
	if(index_Wdecay1==-1)
		return PtEtaPhiEVector();
	
	int index_Wdecay2 = find_match(jet_match_masks, TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_tbar, matched_indicies);
	matched_indicies.push_back(index_Wdecay2);
	if(index_Wdecay2==-1)
		return PtEtaPhiEVector();


	return jets[index_b] + jets[index_Wdecay1] + jets[index_Wdecay2];
}

float mass_from_lvec_in_gev(PtEtaPhiEVector jet){
	if(jet.M()<1e-3)
		return -999-3;
	return jet.M()*1e-3;
}
