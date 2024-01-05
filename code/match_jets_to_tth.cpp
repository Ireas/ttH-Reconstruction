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
// TODO:No multiple calculations by using lazy actions
// TODO:improve matching function




// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 1;
const plog::Severity LOG_LEVEL = plog::Severity::debug; // set logger output severity filer


const std::string INPUT_PATH = "/home/ireas/git_repos/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const std::string OUTPUT_PATH = "/home/ireas/git_repos/master/output/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};


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
// calculates deltaR for jet and obj
float DeltaR(PtEtaPhiEVector jetLvec, PtEtaPhiMVector truthLvec);


// generate lorentz vector for truth object
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass);
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy);

// generate lorentz vectors for jets
vector<PtEtaPhiEVector> GenerateLvecs(vector<Float_t> pts, vector<Float_t> etas, vector<Float_t> phis, vector<Float_t> energies);


// generate jet match mask using bit ids
vector<int> GenerateJetMatchMasks(
	vector<PtEtaPhiEVector> jetLvecs, 
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar 
);

// find closets deltaR jet for object
int FindBestJetIndex(vector<PtEtaPhiEVector> jetLvecs, PtEtaPhiMVector truthLvec, vector<int> unavailableJetIndicies);

// find matched jet for object
int FindFirstMatchIndex(vector<int> jetMatchMasks, int truthMatchValue, vector<int> unavailableJetIndicies);


// combine jets to parent object
PtEtaPhiEVector CombineTop(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);
PtEtaPhiEVector CombineTopBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);


// return mass in gev for plotting
float MassInGeV(PtEtaPhiEVector jet);



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
	auto rLoopManager = rDataFrame.Range(100000); // limit input for testing


	// check if chains are matched properly
	auto nTotalEvents = rLoopManager.Count();
	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0){ 
		PLOG_ERROR << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!";
		exit(1);
	}


	// generate jet lorentz vectors
	rLoopManager = rLoopManager.Define("jet_lvecs", GenerateLvecs, {"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"});
	

	// generate truth obj lorentz vectors
	rLoopManager = rLoopManager.Define(
		"lvec_b_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_b_from_t_pt", "truth.Tth_MC_b_from_t_eta", "truth.Tth_MC_b_from_t_phi", "truth.Tth_MC_b_from_t_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_b_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_b_from_tbar_pt", "truth.Tth_MC_b_from_tbar_eta", "truth.Tth_MC_b_from_tbar_phi", "truth.Tth_MC_b_from_tbar_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay1_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay1_from_t_pt", "truth.Tth_MC_Wdecay1_from_t_eta", "truth.Tth_MC_Wdecay1_from_t_phi", "truth.Tth_MC_Wdecay1_from_t_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay2_from_t", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay2_from_t_pt", "truth.Tth_MC_Wdecay2_from_t_eta", "truth.Tth_MC_Wdecay2_from_t_phi", "truth.Tth_MC_Wdecay2_from_t_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay1_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay1_from_tbar_pt", "truth.Tth_MC_Wdecay1_from_tbar_eta", "truth.Tth_MC_Wdecay1_from_tbar_phi", "truth.Tth_MC_Wdecay1_from_tbar_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay2_from_tbar", 
		GenerateLorentzVectorM, 
		{"truth.Tth_MC_Wdecay2_from_tbar_pt", "truth.Tth_MC_Wdecay2_from_tbar_eta", "truth.Tth_MC_Wdecay2_from_tbar_phi", "truth.Tth_MC_Wdecay2_from_tbar_m"}
	);


	// matching jets and objs (fixed object order)
	rLoopManager = rLoopManager.Define(
		"jet_match_mask", 
		GenerateJetMatchMasks, 
		{"jet_lvecs", "lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar"}
	);

	
	// combine jets to top and anti-top quarks
	rLoopManager = rLoopManager.Define(
		"t_combined",
		CombineTop,
		{"jet_lvecs", "jet_match_mask"}
	);
	
	rLoopManager = rLoopManager.Define(
		"tbar_combined",
		CombineTopBar,
		{"jet_lvecs", "jet_match_mask"}
	);


	// get mass from top and anti_top quarks
	rLoopManager = rLoopManager.Define(
		"t_combined_m_GeV",
		MassInGeV,
		{"t_combined"}
	);
	
	rLoopManager = rLoopManager.Define(
		"tbar_combined_m_GeV",
		MassInGeV,
		{"tbar_combined"}
	);




	// accessing and printing information
	PLOG_DEBUG << "Start: accessing RDataFrame";
	rLoopManager.Display("jet_match_mask")->Print();
	rLoopManager.Display("t_combined")->Print();
	rLoopManager.Display("tbar_combined")->Print();
	

	// finishing
	PLOG_DEBUG << "Start: saving snapshot";
	rLoopManager.Snapshot("matched", OUTPUT_PATH+"rdataframes_output.root");
	
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================
float DeltaR(PtEtaPhiEVector jetLvec, PtEtaPhiMVector truthLvec){
	return sqrt( (jetLvec.eta()-truthLvec.eta())*(jetLvec.eta()-truthLvec.eta()) + (jetLvec.phi()-truthLvec.phi())*(jetLvec.phi()-truthLvec.phi()) );
}

PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass){return PtEtaPhiMVector(pt,eta,phi,mass);}
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy){return PtEtaPhiEVector(pt,eta,phi,energy);}

vector<PtEtaPhiEVector> GenerateLvecs(vector<Float_t> pts, vector<Float_t> etas, vector<Float_t> phis, vector<Float_t> energies){
	vector<PtEtaPhiEVector> lvecs;
	for(int i=0; i<pts.size(); i++){
		PtEtaPhiEVector lvec(pts[i], etas[i], phis[i], energies[i]);
		lvecs.push_back(lvec);
	}
	return lvecs;
}



vector<int> GenerateJetMatchMasks(
	vector<PtEtaPhiEVector> jetLvecs, 
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar 
){
	vector<int> jetMatchMasks(jetLvecs.size());	
	vector<int> unavailableJetIndicies;
	int bestJetIndex;

	
	// b from t
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecBFromT, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_t;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	
	
	// b from tbar
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecBFromTbar, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::b_from_tbar;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	
	
	// Wdecay1 from t
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecWdecay1FromT, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_t;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	
	
	// Wdecay2 from t
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecWdecay2FromT, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_t;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	
	
	// Wdecay1 from tbar
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecWdecay1FromTbar, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_tbar;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	
	
	// Wdecay2 from tbar
	bestJetIndex = FindBestJetIndex(jetLvecs, truthLvecWdecay2FromT, unavailableJetIndicies);
	if(bestJetIndex!=-1){
		jetMatchMasks[bestJetIndex]+= 1<<TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_tbar;
		unavailableJetIndicies.push_back(bestJetIndex);
	}	


	return jetMatchMasks;
}


int FindBestJetIndex(vector<PtEtaPhiEVector> jetLvecs, PtEtaPhiMVector truthLvec, vector<int> unavailableJetIndicies={}){
	int bestJetIndex = -1;
	float bestDeltaR = 999;

	for(int i=0; i<jetLvecs.size(); i++){
		// check if index is available
		bool indexIsAvailable = true;
		for(int j=0; j<unavailableJetIndicies.size(); j++){
			if(i==unavailableJetIndicies[j]){
				indexIsAvailable = false;
				break;
			}
		}
		
		// skip if index is not available
		if(!indexIsAvailable)
			continue;

		// calculate deltaR if index is available
		float currentDeltaR = DeltaR(jetLvecs[i], truthLvec);
		if(currentDeltaR<=DELTA_R_THRESHOLD && currentDeltaR<bestDeltaR){
			bestDeltaR = currentDeltaR;
			bestJetIndex = i;
		}
	}

	return bestJetIndex;
}


int FindFirstMatchIndex(vector<int> jetMatchMasks, int truthMatchValue, vector<int> unavailableJetIndicies={}){
	for(int i=0; i<jetMatchMasks.size(); i++){
		// check if index is available
		bool indexIsAvailable = true;
		for(int j=0; j<unavailableJetIndicies.size(); j++){
			if(i==unavailableJetIndicies[j]){
				indexIsAvailable = false;
				break;
			}
		}
		
		// skip if index is not available
		if(!indexIsAvailable)
			continue;	

		// check if jetMatchMask has truthMatchValue
		if( (jetMatchMasks[i] & 1<<truthMatchValue)!=0 ){
			return i;
		}
	}

	return -1;
}


PtEtaPhiEVector CombineTop(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "ERROR: jetLvecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;
	
	// b from t
	int jetIndexB = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::b_from_t, unavailableJetIndicies);
	if(jetIndexB==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexB);

	// Wdecay1 from t
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_t, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from t
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_t, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}


PtEtaPhiEVector CombineTopBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "ERROR: jetLVecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;
	
	// b from tbar
	int jetIndexB = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::b_from_tbar, unavailableJetIndicies);
	if(jetIndexB==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexB);

	// Wdecay1 from tbar
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::Wdecay1_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from tbar
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_OBJ_MATCH_VALUES::Wdecay2_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}



float MassInGeV(PtEtaPhiEVector jetLvec){ // simple function for testing output
	if(jetLvec.M()*1e-3<1)
		return -999;
	return jetLvec.M()*1e-3;
}


// =========  OLD CODE  ======================
// ===========================================
//vector<int> FindAllJetIndicies(vector<PtEtaPhiEVector> jets, PtEtaPhiMVector truth_obj){
//	vector<int> all_jet_indicies;
//
//	for(int i=0; i<jets.size(); i++){
//		float current_delta_R = DeltaR(jets[i], truth_obj);
//		if(current_delta_R<=DELTA_R_THRESHOLD){
//			all_jet_indicies.push_back(i);
//		}
//	}
//
//	return all_jet_indicies;
//}
