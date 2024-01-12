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
// TODO: Verify no multiple calculations by using lazy actions




// ==========  CONSTANTS  ==========
// =================================
const float DELTA_R_THRESHOLD = 0.4;
const plog::Severity LOG_LEVEL = plog::Severity::verbose; // set logger output severity filter


const string INPUT_PATH = "/home/ireas/git_repos/master/data/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const string OUTPUT_PATH = "/home/ireas/git_repos/master/output/";

const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};


const initializer_list<string> OUTPUT_COLOUMN_NAMES = { // put into array for easier access
	"jet_match_mask",
	"indicies",
	"jet_e_NOSYS",
	"jet_pt_NOSYS",
	"jet_eta",
	"jet_phi",
	"nJets",
	"jet_lvecs",
	"reconstructed_t",
	"reconstructed_tbar",
	"reconstructed_W_from_t",
	"reconstructed_W_from_tbar",
	"successful_reconstruction",
};

// bit-shift amounts for matching
enum TRUTH_MATCH_VALUES{
	b_from_t = 0,
	b_from_tbar = 1,
	Wdecay1_from_t = 2,
	Wdecay2_from_t = 3,
	Wdecay1_from_tbar = 4,
	Wdecay2_from_tbar = 5,
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
// generate lorentz vector for truth object
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass);
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy);

// generate lorentz vectors for jets
vector<PtEtaPhiEVector> GenerateLvecs(vector<Float_t> pts, vector<Float_t> etas, vector<Float_t> phis, vector<Float_t> energies);
int GetNumberOfJets(vector<PtEtaPhiEVector> jetLvecs);

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
int GenerateJetMatchMask(
	PtEtaPhiEVector jetLvec, 
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar 
);
	

// find matched jet for object
int FindFirstMatchIndex(vector<int> jetMatchMasks, int truthMatchValue, vector<int> unavailableJetIndicies);


// combine jets to parent object
PtEtaPhiEVector CombineW(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);
PtEtaPhiEVector CombineWBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);
PtEtaPhiEVector CombineTop(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);
PtEtaPhiEVector CombineTopBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);


// return mass in gev for plotting
float MassInGeV(PtEtaPhiEVector jet);
float MassDifference(PtEtaPhiEVector reconstructedWLvec, float truthMass);

float EvaluateReconstruction(PtEtaPhiEVector reconstructedTLvec, PtEtaPhiEVector reconstructedTbarLvec, PtEtaPhiEVector reconstructedWfromTLvec, PtEtaPhiEVector reconstructedWfromTbarLvec);
int total = 0;	
int nSuccessfulReconstruction = 0;
int nSuccessfulReconstructionT = 0;
int nSuccessfulReconstructionTbar = 0;
int nSuccessfulReconstructionWfromT = 0;
int nSuccessfulReconstructionWfromTbar = 0;


// Get indicies of matched truth
vector<int> GenerateIndiciesFixed(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks);



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
	rLoopManager = rLoopManager.Define(
		"jet_lvecs", 
		GenerateLvecs, 
		{"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"nJets", 
		GetNumberOfJets, 
		{"jet_lvecs"}
	);
	
	auto rLoopManagerFiltered = rLoopManager.Filter("nJets>=6");


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


	// matching jets and objs (fixed object order)
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"jet_match_mask", 
		GenerateJetMatchMasks, 
		{"jet_lvecs", "lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar"}
	);

	
	// combine jets to top and anti-top quarks
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_t",
		CombineTop,
		{"jet_lvecs", "jet_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_tbar",
		CombineTopBar,
		{"jet_lvecs", "jet_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_t",
		CombineTop,
		{"jet_lvecs", "jet_match_mask"}
	);
	
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"reconstructed_W_from_tbar",
		CombineTopBar,
		{"jet_lvecs", "jet_match_mask"}
	);


	// get indicies
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"indicies",
		GenerateIndiciesFixed,
		{"jet_lvecs", "jet_match_mask"}
	);


	// accessing and printing information
	PLOG_DEBUG << "Start: accessing RDataFrame";
	rLoopManagerFiltered.Display("jet_match_mask")->Print();
		
	
	// evaluate t and tbar reconstruction
	rLoopManagerFiltered = rLoopManagerFiltered.Define(
		"successful_reconstruction",
		EvaluateReconstruction,
		{"reconstructed_t", "reconstructed_tbar", "reconstructed_W_from_t", "reconstructed_W_from_tbar"}
	);
	

	// finishing
	PLOG_DEBUG << "Start: saving snapshot";
	rLoopManagerFiltered.Snapshot(
		"matched", 
		OUTPUT_PATH+"rdataframes_output.root",
		OUTPUT_COLOUMN_NAMES
);
	

	// print evaluation after lazy action is done
	float evalRatio = ((float)nSuccessfulReconstruction/total);
	float evalRatioT = ((float)nSuccessfulReconstructionT/total);
	float evalRatioTbar = ((float)nSuccessfulReconstructionTbar/total);
	float evalRatioWfromT = ((float)nSuccessfulReconstructionWfromT/total);
	float evalRatioWfromTbar = ((float)nSuccessfulReconstructionWfromTbar/total);

	PLOG_VERBOSE << "Result: Total: " << total;
	PLOG_VERBOSE << "Result: Successful Reconstruction: " << nSuccessfulReconstruction << " (" << evalRatio << ")";
	PLOG_VERBOSE << "Result: Successful Reconstruction t: " << nSuccessfulReconstructionT << " (" << evalRatioT << ")";
	PLOG_VERBOSE << "Result: Successful Reconstruction tbar: " << nSuccessfulReconstructionTbar << " (" << evalRatioTbar << ")";
	PLOG_VERBOSE << "Result: Successful Reconstruction W form t: " << nSuccessfulReconstructionWfromT << " (" << evalRatioWfromT << ")";
	PLOG_VERBOSE << "Result: Successful Reconstruction W form tbar: " << nSuccessfulReconstructionWfromTbar << " (" << evalRatioWfromTbar << ")";
		
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================
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

int GetNumberOfJets(vector<PtEtaPhiEVector> jetLvecs){
	return jetLvecs.size();
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
	int bestJetIndex;

	
	for(int i=0; i<jetMatchMasks.size(); i++){
		jetMatchMasks[i] = GenerateJetMatchMask(
			jetLvecs[i], 
			truthLvecBFromT, 
			truthLvecBFromTbar, 
			truthLvecWdecay1FromT, 
			truthLvecWdecay2FromT, 
			truthLvecWdecay1FromTbar, 
			truthLvecWdecay2FromTbar
		);
	}

	return jetMatchMasks;
}



int GenerateJetMatchMask(
	PtEtaPhiEVector jetLvec, 
	PtEtaPhiMVector truthLvecBFromT, 
	PtEtaPhiMVector truthLvecBFromTbar, 
	PtEtaPhiMVector truthLvecWdecay1FromT, 
	PtEtaPhiMVector truthLvecWdecay2FromT, 
	PtEtaPhiMVector truthLvecWdecay1FromTbar, 
	PtEtaPhiMVector truthLvecWdecay2FromTbar 
){	
	int jetMatchMask = 0;
	float currentDeltaR = 0;

	// b from t
	currentDeltaR = DeltaR(jetLvec, truthLvecBFromT);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::b_from_t;	
	}
	
	// b from tbar
	currentDeltaR = DeltaR(jetLvec, truthLvecBFromTbar);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::b_from_tbar;	
	}
	
	// Wdecay1 from t
	currentDeltaR = DeltaR(jetLvec, truthLvecWdecay1FromT);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::Wdecay1_from_t;	
	}
	
	// Wdecay2 from t
	currentDeltaR = DeltaR(jetLvec, truthLvecWdecay2FromT);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::Wdecay2_from_t;	
	}
	
	// Wdecay1 from tbar
	currentDeltaR = DeltaR(jetLvec, truthLvecWdecay1FromTbar);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::Wdecay1_from_tbar;	
	}
	
	// Wdecay2 from tbar
	currentDeltaR = DeltaR(jetLvec, truthLvecWdecay2FromTbar);
	if(currentDeltaR<=DELTA_R_THRESHOLD){
		jetMatchMask+= 1<<TRUTH_MATCH_VALUES::Wdecay2_from_tbar;	
	}


	return jetMatchMask;
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



PtEtaPhiEVector CombineW(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "ERROR: jetLvecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;

	// Wdecay1 from t
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_t, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from t
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_t, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}



PtEtaPhiEVector CombineWBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "jetLVecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;
	

	// Wdecay1 from tbar
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from tbar
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}



PtEtaPhiEVector CombineTop(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "ERROR: jetLvecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;
	
	// b from t
	int jetIndexB = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::b_from_t, unavailableJetIndicies);
	if(jetIndexB==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexB);

	// Wdecay1 from t
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_t, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from t
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_t, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}



PtEtaPhiEVector CombineTopBar(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	if(jetLvecs.size()!=jetMatchMasks.size()){
		PLOG_ERROR << "jetLVecs size does not match jetMatchMasks size";
		exit(1);
	}

	vector<int> unavailableJetIndicies;
	
	// b from tbar
	int jetIndexB = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::b_from_tbar, unavailableJetIndicies);
	if(jetIndexB==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexB);

	// Wdecay1 from tbar
	int jetIndexWdecay1 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay1==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay1);
	
	// Wdecay2 from tbar
	int jetIndexWdecay2 = FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_tbar, unavailableJetIndicies);
	if(jetIndexWdecay2==-1){
		return PtEtaPhiEVector(0,0,0,0);
	}	
	unavailableJetIndicies.push_back(jetIndexWdecay2);


	return jetLvecs[jetIndexB] + jetLvecs[jetIndexWdecay1] + jetLvecs[jetIndexWdecay2];
}



float MassInGeV(PtEtaPhiEVector jetLvec){ // simple function for testing output
	if(jetLvec.M()<1){
		return -999;
	}
	return jetLvec.M()*1e-3;
}


float MassDifference(PtEtaPhiEVector reconstructedWLvec, float truthMass){
	if(reconstructedWLvec.M()<1){
		return -999e3;
	}
	return reconstructedWLvec.M() - truthMass;
}



float EvaluateReconstruction(PtEtaPhiEVector reconstructedTLvec, PtEtaPhiEVector reconstructedTbarLvec, PtEtaPhiEVector reconstructedWfromTLvec, PtEtaPhiEVector reconstructedWfromTbarLvec){
	bool success = false;

	total++;
	if(reconstructedTLvec.M()>1)
		nSuccessfulReconstructionT++;
	if(reconstructedTbarLvec.M()>1)
		nSuccessfulReconstructionTbar++;
	if(reconstructedWfromTLvec.M()>1)
		nSuccessfulReconstructionWfromT++;
	if(reconstructedWfromTbarLvec.M()>1)
		nSuccessfulReconstructionWfromTbar++;

	if(reconstructedTLvec.M()>1 && reconstructedTbarLvec.M()>1 && reconstructedWfromTLvec.M()>1 && reconstructedWfromTbarLvec.M()>1){
		nSuccessfulReconstruction++;
		success = true;
	}

	return success;	
}






// Indicies stuff
vector<int> GenerateIndiciesFixed(vector<PtEtaPhiEVector> jetLvecs, vector<int> jetMatchMasks){
	vector<int> objectIndicies;
	

	// b from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::b_from_t, objectIndicies) );
 
	// Wdecay1 from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_t, objectIndicies) );
 	
	// Wdecay2 from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_t, objectIndicies) );


	// b from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::b_from_tbar, objectIndicies) );
 
	// Wdecay1 from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay1_from_tbar, objectIndicies) );
 	
	// Wdecay2 from tbar
	objectIndicies.push_back( FindFirstMatchIndex(jetMatchMasks, TRUTH_MATCH_VALUES::Wdecay2_from_tbar, objectIndicies) );

	
	return objectIndicies;
}
