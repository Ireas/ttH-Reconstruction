#include<string>
#include<iostream>
#include<vector>
using namespace std;
#include"TROOT.h"
#include"TFile.h"
#include"TLorentzVector.h"
#include"TChain.h"
#include"TTree.h"
#include"TH1F.h"
#include"ROOT/RDataFrame.hxx"
#include<Math/VectorUtil.h> // for DeltaR
#include<Math/Vector4D.h> // for PtEtaPhiEVector and PtEtaPhiMVector
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;

#include<chrono> // for time measurements during the program

#include<plog/Log.h> // c++ logger + init headers
#include<plog/Formatters/TxtFormatter.h>
#include<plog/Initializers/ConsoleInitializer.h>



// ==========  DELTA R MATCHING  ==========
// ========================================
// uses ROOTs modern RDataFrames to simplify the event loop and branchaddress allocation



// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 0.4;


const int MAX_NUMBER_OF_EVENTS = 1e6;//1e6; // set to 0 for no limit
const plog::Severity LOG_LEVEL = plog::Severity::verbose; // set logger output severity filter


const string INPUT_PATH = "/home/ireas/git_repos/master/input/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const string OUTPUT_PATH = "/home/ireas/git_repos/master/output/";


const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};


const initializer_list<string> OUTPUT_COLOUMN_NAMES = { // put into array for easier access
	"jet_potential_match_mask",
	"jet_final_match_mask",
	"jet_to_object_indicies_fixed",
	"jet_lvecs",
	"jet_e_NOSYS",
	"jet_pt_NOSYS",
	"jet_eta",
	"jet_phi",
	"number_of_jets",
	"reconstructed_t_lvec",
	"reconstructed_tbar_lvec",
	"reconstructed_W_from_t_lvec",
	"reconstructed_W_from_tbar_lvec",
	"reconstructed_HW_lvec",
	"reconstructed_t_m",
	"reconstructed_tbar_m",
	"reconstructed_W_from_t_m",
	"reconstructed_W_from_tbar_m",
	"truth_t_m",
	"truth_tbar_m",
	"truth_W_from_t_m",
	"truth_W_from_tbar_m",
	"reconstructed_t_pt",
	"reconstructed_tbar_pt",
	"reconstructed_W_from_t_pt",
	"reconstructed_W_from_tbar_pt",
	"truth_t_pt",
	"truth_tbar_pt",
	"truth_W_from_t_pt",
	"truth_W_from_tbar_pt",
	"successful_reconstruction",
	"higgs_decay_mode_custom",
	"higgs_decay_decay_mode",
	"higgs_decay1_decay1_filtered_pdgid",
	"higgs_decay1_decay2_filtered_pdgid",
	"higgs_decay2_decay1_filtered_pdgid",
	"higgs_decay2_decay2_filtered_pdgid",
	"reco_mc_channel_number",
	"reco_event_number",
	"truth_mc_channel_number",
	"truth_event_number",
};


// indicies and bit-shift amounts for different truth objects
int NUMBER_OF_TRUTH_OBJECTS = 8;
enum TRUTH_PARTONS{
	b_from_t = 0,
	b_from_tbar = 1,
	Wdecay1_from_t = 2,
	Wdecay2_from_t = 3,
	Wdecay1_from_tbar = 4,
	Wdecay2_from_tbar = 5,
	Wdecay1_from_H = 6,
	Wdecay2_from_H = 7,
};

enum HIGGS_DECAY_MODE{
	undefined = -1,
	b_b = 0,
	e_e = 1,
	mu_mu = 2,
	tau_tau = 3,
	y_y = 4,
	w_w = 5,
	z_z = 6,
	other = 7,
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
// generate lorentz vector for truth object
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass);
PtEtaPhiMVector GenerateLorentzVectorMHiggsDecision(Float_t pt1, Float_t eta1, Float_t phi1, Float_t mass1, Int_t pdgId1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t mass2, Int_t pdgId2);
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy);


// generate vectors of lorentz vectors for easier access
vector<PtEtaPhiMVector> GenerateTruthLvecs(
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromH, 
	PtEtaPhiMVector truthLvecWdecay2FromH 
);

vector<PtEtaPhiEVector> GenerateJetLvecs(
	vector<Float_t> pts, 
	vector<Float_t> etas, 
	vector<Float_t> phis, 
	vector<Float_t> energies
);


// generate jet potential match masks (all truth objects within delta R range)
vector<int> GenerateJetPotentialMatchMasks(vector<PtEtaPhiEVector> jetLvecs, vector<PtEtaPhiMVector> truthLvecs); 
int GenerateJetPotentialMatchMask(PtEtaPhiEVector jetLvec, vector<PtEtaPhiMVector> truthLvecs); 


// generate jet final match masks (only closest truth object to jets)	
vector<int> GenerateJetFinalMatchMasks(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetPotentialMatchMasks, vector<PtEtaPhiMVector> truthLvecs);
int ClosestTruthMatchValueToJet(PtEtaPhiEVector jetLvec, int jetPotentialMatchMask, vector<PtEtaPhiMVector> truthLvecs, vector<int> unavailableTruthValues);


// recollect the jet index for each matched truth object in fixxed order
vector<int> CollectJetToObjectIndiciesFixed(vector<int> jetFinalMatchMasks);


// reconstruct matched jets to parent objects
PtEtaPhiEVector ReconstructT(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks);
PtEtaPhiEVector ReconstructTBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks);
PtEtaPhiEVector ReconstructWFromT(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks);
PtEtaPhiEVector ReconstructWFromTBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks);
PtEtaPhiEVector ReconstructHW(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks);


// higgs
int GenerateHiggsDecayModeCustom(int higgsDecay1PdgId, int higgsDecay2PdgId);
vector<int> GenerateHiggsDecayDecayMode(int higgsDecayMode, int higgsDecay11, int higgsDecay12, int higgsDecay21, int higgsDecay22);
int GetFilteredPdgIDs(Int_t pdgId);


int GetNumberOfJets(vector<PtEtaPhiEVector> jetLvecs);
bool CheckReconstruction(PtEtaPhiEVector tLvec, PtEtaPhiEVector tBarLvev, PtEtaPhiEVector HLvec);

float RenameFloat(float target);
int RenameInt(unsigned int target);
int RenameInt2(unsigned long long target);
float GetMass(PtEtaPhiEVector lvec);
float GetPt(PtEtaPhiEVector lvec);



// ==========  MAIN  ==========
// ============================
int main(){
	// ==========  SETUP
	// =================
	// setup plog (logger)
	static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
	plog::init(LOG_LEVEL, &consoleAppender);
	

	// setup ROOT
	PLOG_DEBUG << "Start: setup ROOT";
	//ROOT::EnableImplicitMT();// use multithreading when no range is set 

	
	// setup TChain
	PLOG_DEBUG << "Start: setup TChain";
	TChain rRecoChain("reco");
	TChain rTruthChain("truth");

	for(auto input_file_name:INPUT_FILE_NAMES){
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);

		rTruthChain.Add(file_path.c_str());
		rRecoChain.Add(file_path.c_str());
	}


	// index TruthChain to kick out un-matched events - then declare friends
	rTruthChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rRecoChain.AddFriend(&rTruthChain);


	// setup RDataFrame
	PLOG_DEBUG << "Start: setup RDataFrame";
	auto rDataFrame = RDataFrame(rRecoChain);
	auto rLoopManager = rDataFrame.Range(MAX_NUMBER_OF_EVENTS); // limit input for testing


	// check if chains are matched properly
	auto nTotalEvents = rLoopManager.Count();
	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0){ 
		PLOG_ERROR << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!";
		exit(1);
	}


	// generate jet lorentz vectors
	rLoopManager = rLoopManager.Define(
		"jet_lvecs", 
		GenerateJetLvecs, 
		{"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"number_of_jets", 
		GetNumberOfJets, 
		{"jet_lvecs"}
	);


	// apply jet filter
	auto rLoopManagerFiltered = rLoopManager.Filter("number_of_jets>=8");


	// generate truth obj lorentz vectors
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_b_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_b_from_t_pt", "truth.Tth_MC_b_from_t_eta", "truth.Tth_MC_b_from_t_phi", "truth.Tth_MC_b_from_t_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_b_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_b_from_tbar_pt", "truth.Tth_MC_b_from_tbar_eta", "truth.Tth_MC_b_from_tbar_phi", "truth.Tth_MC_b_from_tbar_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay1_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay1_from_t_pt", "truth.Tth_MC_Wdecay1_from_t_eta", "truth.Tth_MC_Wdecay1_from_t_phi", "truth.Tth_MC_Wdecay1_from_t_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay2_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay2_from_t_pt", "truth.Tth_MC_Wdecay2_from_t_eta", "truth.Tth_MC_Wdecay2_from_t_phi", "truth.Tth_MC_Wdecay2_from_t_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay1_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay1_from_tbar_pt", "truth.Tth_MC_Wdecay1_from_tbar_eta", "truth.Tth_MC_Wdecay1_from_tbar_phi", "truth.Tth_MC_Wdecay1_from_tbar_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay2_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay2_from_tbar_pt", "truth.Tth_MC_Wdecay2_from_tbar_eta", "truth.Tth_MC_Wdecay2_from_tbar_phi", "truth.Tth_MC_Wdecay2_from_tbar_m"}
	);

	// find hadronic H decay and generate turth obj lorentz vectors
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay1_from_H", 
		GenerateLorentzVectorMHiggsDecision, 
		{"truth.Tth_MC_Higgs_decay1_from_decay1_pt", "truth.Tth_MC_Higgs_decay1_from_decay1_eta", "truth.Tth_MC_Higgs_decay1_from_decay1_phi", "truth.Tth_MC_Higgs_decay1_from_decay1_m", "truth.Tth_MC_Higgs_decay1_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay1_from_decay2_pt", "truth.Tth_MC_Higgs_decay1_from_decay2_eta", "truth.Tth_MC_Higgs_decay1_from_decay2_phi", "truth.Tth_MC_Higgs_decay1_from_decay2_m","truth.Tth_MC_Higgs_decay1_from_decay2_pdgId"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"lvec_Wdecay2_from_H", 
		GenerateLorentzVectorMHiggsDecision, 
		{"truth.Tth_MC_Higgs_decay2_from_decay1_pt", "truth.Tth_MC_Higgs_decay2_from_decay1_eta", "truth.Tth_MC_Higgs_decay2_from_decay1_phi", "truth.Tth_MC_Higgs_decay2_from_decay1_m","truth.Tth_MC_Higgs_decay2_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay2_pt", "truth.Tth_MC_Higgs_decay2_from_decay2_eta", "truth.Tth_MC_Higgs_decay2_from_decay2_phi", "truth.Tth_MC_Higgs_decay2_from_decay2_m", "truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);


	// generate truth obj vector for each event for easier access
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_lvecs",
		GenerateTruthLvecs,
		{"lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar", "lvec_Wdecay1_from_H", "lvec_Wdecay2_from_H"}
	);


	// match jets to objs
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"jet_potential_match_mask", 
		GenerateJetPotentialMatchMasks, 
		{"jet_lvecs", "truth_lvecs"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"jet_final_match_mask", 
		GenerateJetFinalMatchMasks, 
		{"jet_lvecs", "jet_potential_match_mask", "truth_lvecs"}
	);
	

	// collect jet indicies for SPANet training
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"jet_to_object_indicies_fixed",
		CollectJetToObjectIndiciesFixed,
		{"jet_final_match_mask"}
	);

	
	// reconstruct objects from matched jets
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_t_lvec",
		ReconstructT,
		{"jet_lvecs", "jet_final_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_tbar_lvec",
		ReconstructTBar,
		{"jet_lvecs", "jet_final_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_t_lvec",
		ReconstructWFromT,
		{"jet_lvecs", "jet_final_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_tbar_lvec",
		ReconstructWFromTBar,
		{"jet_lvecs", "jet_final_match_mask"}
	);

	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_HW_lvec",
		ReconstructHW,
		{"jet_lvecs", "jet_final_match_mask"}
	);

	
	// get mass from reconstructed objects
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_t_m",
		GetMass,
		{"reconstructed_t_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_tbar_m",
		GetMass,
		{"reconstructed_tbar_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_t_m",
		GetMass,
		{"reconstructed_W_from_t_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_tbar_m",
		GetMass,
		{"reconstructed_W_from_tbar_lvec"}
	);

	// get pt from reconstruced object
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_t_pt",
		GetPt,
		{"reconstructed_t_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_tbar_pt",
		GetPt,
		{"reconstructed_tbar_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_t_pt",
		GetPt,
		{"reconstructed_W_from_t_lvec"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_tbar_pt",
		GetPt,
		{"reconstructed_W_from_tbar_lvec"}
	);


	 


	// rename truth trees
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_t_m",
		RenameFloat,
		{"truth.Tth_MC_t_afterFSR_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_tbar_m",
		RenameFloat,
		{"truth.Tth_MC_tbar_afterFSR_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_W_from_t_m",
		RenameFloat,
		{"truth.Tth_MC_W_from_t_m"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_W_from_tbar_m",
		RenameFloat,
		{"truth.Tth_MC_W_from_tbar_m"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_t_pt",
		RenameFloat,
		{"truth.Tth_MC_t_afterFSR_pt"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_tbar_pt",
		RenameFloat,
		{"truth.Tth_MC_tbar_afterFSR_pt"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_W_from_t_pt",
		RenameFloat,
		{"truth.Tth_MC_W_from_t_pt"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_W_from_tbar_pt",
		RenameFloat,
		{"truth.Tth_MC_W_from_tbar_pt"}
	);


	// save eventNumber
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reco_mc_channel_number",
		RenameInt,
		{"mcChannelNumber"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_mc_channel_number",
		RenameInt,
		{"truth.mcChannelNumber"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"truth_event_number",
		RenameInt2,
		{"truth.eventNumber"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reco_event_number",
		RenameInt2,
		{"truth.eventNumber"}
	);


	// check if event is fully reconstructed
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"successful_reconstruction",
		CheckReconstruction,
		{"reconstructed_t_lvec", "reconstructed_tbar_lvec", "reconstructed_HW_lvec"}
	);


	// get higgs information
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay_mode_custom",
		GenerateHiggsDecayModeCustom,
		{"truth.Tth_MC_Higgs_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_pdgId"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay_decay_mode",
		GenerateHiggsDecayDecayMode,
		{"higgs_decay_mode_custom", "truth.Tth_MC_Higgs_decay1_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay1_from_decay2_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);

	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay1_decay1_filtered_pdgid",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay1_from_decay1_pdgId"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay1_decay2_filtered_pdgid",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay2_from_decay1_pdgId"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay2_decay1_filtered_pdgid",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay1_from_decay2_pdgId"}
	);
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"higgs_decay2_decay2_filtered_pdgid",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);


	// save snapshot to disk
	PLOG_DEBUG << "Start: saving snapshot";
	rLoopManagerFiltered.Snapshot(
		"matched", 
		OUTPUT_PATH+"rdataframes_output.root",
		OUTPUT_COLOUMN_NAMES
	);
		
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================

// ==========  GENERATE LORENTZ VECTORS
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass){return PtEtaPhiMVector(pt,eta,phi,mass);}
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy){return PtEtaPhiEVector(pt,eta,phi,energy);}
PtEtaPhiMVector GenerateLorentzVectorMHiggsDecision(Float_t pt1, Float_t eta1, Float_t phi1, Float_t mass1, Int_t pdgId1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t mass2, Int_t pdgId2){
	if(abs(pdgId1)>=1 && abs(pdgId1)<=8){
		return PtEtaPhiMVector(pt1,eta1,phi1,mass1);
	}

	if(abs(pdgId2)>=1 && abs(pdgId2)<=8){
		return PtEtaPhiMVector(pt2,eta2,phi2,mass2);
	}
	
	
	//PLOG_WARNING << "Warning: found event with no hadronic H decay! (" << pdgId1 << "," << pdgId2 << ")";
	return PtEtaPhiMVector(0,0,0,0);
}



vector<PtEtaPhiMVector> GenerateTruthLvecs(
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromH, 
	PtEtaPhiMVector truthLvecWdecay2FromH 
){
	vector<PtEtaPhiMVector> truthLvecs(NUMBER_OF_TRUTH_OBJECTS);
	
	truthLvecs[TRUTH_PARTONS::b_from_t] = truthLvecBFromT;
	truthLvecs[TRUTH_PARTONS::b_from_tbar] = truthLvecBFromTbar;
	truthLvecs[TRUTH_PARTONS::Wdecay1_from_t] = truthLvecWdecay1FromT;
	truthLvecs[TRUTH_PARTONS::Wdecay2_from_t] = truthLvecWdecay2FromT;
	truthLvecs[TRUTH_PARTONS::Wdecay1_from_tbar] = truthLvecWdecay1FromTbar;
	truthLvecs[TRUTH_PARTONS::Wdecay2_from_tbar] = truthLvecWdecay2FromTbar;
	truthLvecs[TRUTH_PARTONS::Wdecay1_from_H] = truthLvecWdecay1FromH;
	truthLvecs[TRUTH_PARTONS::Wdecay2_from_H] = truthLvecWdecay2FromH;

	return truthLvecs;
}


vector<PtEtaPhiEVector> GenerateJetLvecs(
	vector<Float_t> pts,
	vector<Float_t> etas, 
	vector<Float_t> phis, 
	vector<Float_t> energies
){
	vector<PtEtaPhiEVector> jetLvecs;
	for(int i=0; i<pts.size(); i++){
		PtEtaPhiEVector jetLvec(pts[i], etas[i], phis[i], energies[i]);
		jetLvecs.push_back(jetLvec);
	}
	return jetLvecs;
}



// ==========  GENERATE POTENTIAL MATCH MASKS AND FINAL MATCH MASKS
vector<int> GenerateJetPotentialMatchMasks(vector<PtEtaPhiEVector> jetLvecs, vector<PtEtaPhiMVector> truthLvecs){
	vector<int> jetPotentialMatchMasks(jetLvecs.size());	
	
	for(int i=0; i<jetPotentialMatchMasks.size(); i++)
		jetPotentialMatchMasks[i] = GenerateJetPotentialMatchMask(jetLvecs[i], truthLvecs);

	return jetPotentialMatchMasks;
}


int GenerateJetPotentialMatchMask(PtEtaPhiEVector jetLvec, vector<PtEtaPhiMVector> truthLvecs){ 
	int jetPotentialMatchMask = 0;
	
	// iterate all truth objects
	for(int i=0; i<truthLvecs.size(); i++){ // i equals the corresponding TRUTH_MATCH_VALUE
		float currentDeltaR = DeltaR(jetLvec, truthLvecs[i]);
		if(currentDeltaR<=DELTA_R_THRESHOLD)
			jetPotentialMatchMask+= 1<<i;
	}

	return jetPotentialMatchMask;
}


vector<int> GenerateJetFinalMatchMasks(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetPotentialMatchMasks, vector<PtEtaPhiMVector> truthLvecs){
	vector<int> jetFinalMatchMasks(jetLvecs.size());	
	vector<int> unavailableTruthValues;	
	
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		int closestTruthMatchValue = ClosestTruthMatchValueToJet(jetLvecs[i], jetPotentialMatchMasks[i], truthLvecs, unavailableTruthValues);

		if(closestTruthMatchValue!=-1){ // closest available truth object is matched
			jetFinalMatchMasks[i] = 1<<closestTruthMatchValue;
			unavailableTruthValues.push_back(closestTruthMatchValue);
		}	
		else{
			jetFinalMatchMasks[i] = 0; // no truth object has been matched successfully, unmatched jet
		}	

	}

	return jetFinalMatchMasks;
}


int ClosestTruthMatchValueToJet(PtEtaPhiEVector jetLvec, int jetPotentialMatchMask, vector<PtEtaPhiMVector> truthLvecs, vector<int> unavailableTruthValues){ 
	float bestValue = -1;
	int bestTruthMatchValue = -1;
		
	// check all truth objects for closest deltaR match
	for(int i=0; i<truthLvecs.size(); i++){// i equals the corresponding TRUTH_MATCH_VALUE 
		if( (jetPotentialMatchMask & 1<<i)==0 ) // jet is not potentially matched to truth, skipping 
			continue;

		// take best 
		float currentValue = DeltaR(jetLvec, truthLvecs[i]);

		if( currentValue <= bestValue || bestValue<0 ){
			bool truthValueIsAvailable = true;
			// check if new best truth match value is available
			for(int j=0; j<unavailableTruthValues.size(); j++){
				if(i==unavailableTruthValues[j]){
					truthValueIsAvailable = false;
					break;
				}
			}
			
			// if new best truth match value is available update values
			if(truthValueIsAvailable){
				bestValue = currentValue;
				bestTruthMatchValue = i;
			}
		}
	}
		
	return bestTruthMatchValue;
}


// ==========  RECONSTRUCT OBJECTS FROM JETS
PtEtaPhiEVector ReconstructT(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks){
	int jetIndexB = -1;	
	int jetIndexWdecay1 = -1;	
	int jetIndexWdecay2 = -1;	

	//  matched jets if possible
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		// match b quark
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::b_from_t)!=0 ){
			jetIndexB = i;
			continue;
		}
		
		// match Wdecay1
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay1_from_t)!=0 ){
			jetIndexWdecay1 = i;
			continue;
		}
		
		// match Wdecay2
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay2_from_t)!=0 ){
			jetIndexWdecay2 = i;
			continue;
		}
	}

	// check if any object is not matched
	if(jetIndexB==-1 || jetIndexWdecay1==-1 || jetIndexWdecay2==-1)
		return PtEtaPhiEVector(0,0,0,0);
	
	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}
	
PtEtaPhiEVector ReconstructTBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks){
	int jetIndexB = -1;	
	int jetIndexWdecay1 = -1;	
	int jetIndexWdecay2 = -1;	

	//  matched jets if possible
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		// match b quark
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::b_from_tbar)!=0 ){
			jetIndexB = i;
			continue;
		}
		
		// match Wdecay1
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay1_from_tbar)!=0 ){
			jetIndexWdecay1 = i;
			continue;
		}
		
		// match Wdecay2
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay2_from_tbar)!=0 ){
			jetIndexWdecay2 = i;
			continue;
		}
	}

	// check if any object is not matched
	if(jetIndexB==-1 || jetIndexWdecay1==-1 || jetIndexWdecay2==-1)
		return PtEtaPhiEVector(0,0,0,0);
	
	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}


PtEtaPhiEVector ReconstructWFromT(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks){
	int jetIndexWdecay1 = -1;	
	int jetIndexWdecay2 = -1;	

	//  matched jets if possible
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		// match Wdecay1
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay1_from_t)!=0 ){
			jetIndexWdecay1 = i;
			continue;
		}
		
		// match Wdecay2
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay2_from_t)!=0 ){
			jetIndexWdecay2 = i;
			continue;
		}
	}

	// check if any object is not matched
	if(jetIndexWdecay1==-1 || jetIndexWdecay2==-1)
		return PtEtaPhiEVector(0,0,0,0);
	
	return jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}

PtEtaPhiEVector ReconstructWFromTBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks){
	int jetIndexWdecay1 = -1;	
	int jetIndexWdecay2 = -1;	

	//  matched jets if possible
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		// match Wdecay1
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay1_from_tbar)!=0 ){
			jetIndexWdecay1 = i;
			continue;
		}
		
		// match Wdecay2
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay2_from_tbar)!=0 ){
			jetIndexWdecay2 = i;
			continue;
		}
	}

	// check if any object is not matched
	if(jetIndexWdecay1==-1 || jetIndexWdecay2==-1)
		return PtEtaPhiEVector(0,0,0,0);
	
	return jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}


PtEtaPhiEVector ReconstructHW(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetFinalMatchMasks){
	int jetIndexWdecay1 = -1;	
	int jetIndexWdecay2 = -1;	

	//  matched jets if possible
	for(int i=0; i<jetFinalMatchMasks.size(); i++){
		// match Wdecay1
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay1_from_H)!=0 ){
			jetIndexWdecay1 = i;
			continue;
		}
		
		// match Wdecay2
		if( (jetFinalMatchMasks[i] & 1<<TRUTH_PARTONS::Wdecay2_from_H)!=0 ){
			jetIndexWdecay2 = i;
			continue;
		}
	}

	// check if any object is not matched
	if(jetIndexWdecay1==-1 || jetIndexWdecay2==-1)
		return PtEtaPhiEVector(0,0,0,0);
	
	return jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}


// ==========  COLLECT JET INDICIES FROM MATCHED OBJECTS
vector<int> CollectJetToObjectIndiciesFixed(vector<int> jetFinalMatchMasks){
	vector<int> jetToObjectIndicies(NUMBER_OF_TRUTH_OBJECTS);
	
	for(int i=0; i<NUMBER_OF_TRUTH_OBJECTS; i++){// i equals the corresponding TRUTH_PARTON value
		jetToObjectIndicies[i] = -1;
		for(int j=0; j<jetFinalMatchMasks.size(); j++){
			if( (jetFinalMatchMasks[j] & 1<<i)!=0 ){
				jetToObjectIndicies[i] = j;
				break;
			}
		}
	}
	
	return jetToObjectIndicies;
}


// ==========  HIGGS
int GenerateHiggsDecayModeCustom(int higgsDecay1PdgId, int higgsDecay2PdgId){
	int pdgid1 = -1;
	int pdgid2 = -1;
	
	higgsDecay1PdgId = abs(higgsDecay1PdgId);
	higgsDecay2PdgId = abs(higgsDecay2PdgId);

	if(higgsDecay1PdgId!=higgsDecay1PdgId)
		return HIGGS_DECAY_MODE::undefined;

	if(higgsDecay1PdgId>999)
		return HIGGS_DECAY_MODE::undefined;
	
	if(higgsDecay1PdgId==5)
		return HIGGS_DECAY_MODE::b_b;
	
	if(higgsDecay1PdgId==11)
		return HIGGS_DECAY_MODE::e_e;
	
	if(higgsDecay1PdgId==13)
		return HIGGS_DECAY_MODE::mu_mu;
	
	if(higgsDecay1PdgId==15)
		return HIGGS_DECAY_MODE::tau_tau;
	
	if(higgsDecay1PdgId==22)
		return HIGGS_DECAY_MODE::y_y;
	
	if(higgsDecay1PdgId==23)
		return HIGGS_DECAY_MODE::z_z;
	
	if(higgsDecay1PdgId==24)
		return HIGGS_DECAY_MODE::w_w;
	
	return HIGGS_DECAY_MODE::other;
}

vector<int> GenerateHiggsDecayDecayMode(int higgsDecayMode, Int_t higgsDecay11, Int_t higgsDecay12, Int_t higgsDecay21, Int_t higgsDecay22){
	if(higgsDecayMode!=HIGGS_DECAY_MODE::w_w && higgsDecayMode!=HIGGS_DECAY_MODE::z_z){
		return vector<int>{-1};
	}

	int indexHiggsDecay11 = -1;
	int indexHiggsDecay12 = -1;
	int indexHiggsDecay21 = -1;
	int indexHiggsDecay22 = -1;

	if(abs(higgsDecay11)<1000)
		indexHiggsDecay11 = abs(higgsDecay11);
	if(abs(higgsDecay12)<1000)
		indexHiggsDecay12 = abs(higgsDecay12);
	if(abs(higgsDecay21)<1000)
		indexHiggsDecay21 = abs(higgsDecay21);
	if(abs(higgsDecay22)<1000)
		indexHiggsDecay22 = abs(higgsDecay22);

	if(indexHiggsDecay11<0)
		indexHiggsDecay11 = -1;
	if(indexHiggsDecay12<0)
		indexHiggsDecay12 = -1;
	if(indexHiggsDecay21<0)
		indexHiggsDecay21 = -1;
	if(indexHiggsDecay22<0)
		indexHiggsDecay22 = -1;

	return vector<int>{indexHiggsDecay11, indexHiggsDecay12, indexHiggsDecay21, indexHiggsDecay22};
}

int GetFilteredPdgIDs(Int_t pdgId){
	int newPdgId = abs(pdgId);
	
	if(abs(newPdgId)>1000 || newPdgId==0){
		return -1;
	}

	return abs(newPdgId);
}



// ==========  HELPER FUNCTONS
int GetNumberOfJets(vector<PtEtaPhiEVector> jetLvecs){return jetLvecs.size();}

bool CheckReconstruction(PtEtaPhiEVector tLvec, PtEtaPhiEVector tBarLvec, PtEtaPhiEVector HLvec){
	if(tLvec.M()>0 && tBarLvec.M()>0 && HLvec.M()>0){
		return true;
	}
	return false;
}

float RenameFloat(float target){return target;}
int RenameInt(unsigned int target){return target;}
int RenameInt2(unsigned long long target){return target;}

float GetMass(PtEtaPhiEVector lvec){return lvec.M();}
float GetPt(PtEtaPhiEVector lvec){return lvec.Pt();}
