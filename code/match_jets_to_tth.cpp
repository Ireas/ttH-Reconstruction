#include<string>
#include<iostream>
#include<vector>
using namespace std;
#include"TROOT.h"
#include"TFile.h"
#include"TLorentzVector.h"
#include"TChain.h"
#include"TTree.h"
#include"ROOT/RDataFrame.hxx"
#include<Math/VectorUtil.h> // for DeltaR
#include<Math/Vector4D.h> // for PtEtaPhiEVector and PtEtaPhiMVector
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;

#include<chrono> // for time measurements during the program



// ==========  DELTA R MATCHING  ==========
// ========================================
// uses ROOTs modern RDataFrames to simplify the event loop and branchaddress allocationso



// ==========  CONSTANTS  ==========
// ================================= 
const int MAX_NUMBER_OF_EVENTS = 1e4; // set to 0 for no limit

const float PDG_MASS_WBOSON = 80369.2; // mass in MeV
const float THRESHOLD_ONSHELL_DEFINITION = 1e3; // maximum deviation from DPG mass in MeV to classify as onshell

const float THRESHOLD_DELTA_R = 0.4; // maximum reco jet deviation from the truth for matching





const string INPUT_PATH = "/home/ireas/git_repos/master/samples/input/v3/";
const char* INPUT_FILE_NAMES[] = { // put into array for easier access
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000001.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000002.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000003.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000004.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000001.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000002.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000003.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000004.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000005.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000006.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000007.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000008.output.root",
//	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043366._000001.output.root",
//	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043370._000001.output.root",
//	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043373._000001.output.root",
//	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043367._000001.output.root",
//	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043367._000002.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043371._000001.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000001.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000002.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000003.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000004.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043368._000001.output.root",
//	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043368._000002.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000001.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000002.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000003.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000004.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000001.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000002.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000003.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000004.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000005.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000006.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000007.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000008.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000009.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000010.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000011.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000012.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000001.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000002.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000003.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000004.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000005.output.root",
//	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000006.output.root",
};

const string OUTPUT_PATH = "/home/ireas/git_repos/master/samples/matched/";
const string OUTPUT_FILE = "_matched.root";
const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
	// logistics
	"mcChannelNumber",
	"eventNumber",
	// event global information
	"number_of_jets",
	"reco_met_value", 
	"reco_met_phi",
	// jet particle information
	"jet_DL1dv01_FixedCutBEff_85_select",
	"lvecs_jets",
	"jet_e_NOSYS",
	"jet_pt_NOSYS",
	"jet_eta",
	"jet_phi",
	"jet_final_match_mask",
	"jet_to_object_indicies_fixed",
	// true event signatures
	"signature_higgs_decay",
	"signature_onshell_whad",
	"signature_abs_lepton_pdgid",
	// reco event classifier
	"classifier_event_status",
	"higgs_decay_mode_custom",
	"higgs_decay_decay_mode",
	// lepton particle information
	"true_lepton_pt",
	"true_lepton_eta",
	"true_lepton_phi",
	"true_lepton_m",
	"true_lepton_lvec",
	"reco_lepton_pt",
	"reco_lepton_eta",
	"reco_lepton_phi",
	"reco_lepton_e",
	"reco_lepton_lvec",
	// neutrino particle information
	"true_neutrino_lvec",
	// whad particle information
	"true_whad_lvec",
	// wlep particle information
	"true_wlep_lvec",
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


enum CLASSIFIER_EVENT_STATUS{
	impossible = -1,
	possible = 0,
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


float ExtractRecoInformationLepton(int classifierLeptonFlavour, vector<float> info_electron, vector<char> pass_electron_selection, vector<float> info_muon, vector<char> pass_muon_selection){
	// lepton is electron
	if(classifierLeptonFlavour==1)
	{
		for(int i=0; i<info_electron.size(); i++)
		{
			if( (bool)pass_electron_selection[i] )
			{
				return info_electron[i];
			}
		}

		// no electron passes selection
		return -2;
	}
		

	// lepton is muon
	if(classifierLeptonFlavour==2)
	{
		for(int i=0; i<info_muon.size(); i++)
		{
			if( (bool)pass_muon_selection[i] )
			{
				return info_muon[i];
			}
		}
		
		// no muon passes selection
		return -2;
	}

	// no lepton detected
	return -1;
}



vector<int> CombineHiggsDecayPDGIDs(int pdgIdHDecay11, int pdgIdHDecay12, int pdgIdHDecay21, int pdgIdHDecay22){
	return vector<int>{pdgIdHDecay11, pdgIdHDecay12, pdgIdHDecay21, pdgIdHDecay22};
}
vector<PtEtaPhiMVector> CombineHiggsDecayLorentzVectors(PtEtaPhiMVector lvecHdecay11, PtEtaPhiMVector lvecHdecay12, PtEtaPhiMVector lvecHdecay21, PtEtaPhiMVector lvecHdecay22){
	return vector<PtEtaPhiMVector>{lvecHdecay11, lvecHdecay12, lvecHdecay21, lvecHdecay22};
}



int ClassifyRecoLeptonFlavour(char passElectronChar, char passMuonChar)
{
	bool passElectron = (bool)passElectronChar;
	bool passMuon = (bool)passMuonChar;

	// both
	if( passElectron && passMuon )
	{
		std::cout << "?" << std::endl;
		return 3;
	}
	// electron
	if(passElectron)
	{
		return 1;
	}
	// muon
	if(passMuon)
	{
		return 2;
	}

	return 0;
}

int ClassifyEventStatus(vector<int> fixedAssignment){
	for(int i=0; i<NUMBER_OF_TRUTH_OBJECTS; i++){
		if(fixedAssignment[i]==-1){
			return CLASSIFIER_EVENT_STATUS::impossible;
		}
	}

	return CLASSIFIER_EVENT_STATUS::possible;
}

PtEtaPhiMVector GenerateLorentzVectorNeutrinoTrue(vector<PtEtaPhiMVector> lvecsHdecay, vector<int> pdgIdsHdecay){
	for(int i=0; i<pdgIdsHdecay.size(); i++){
		if(abs(pdgIdsHdecay[i])==12 || abs(pdgIdsHdecay[i])==14 || abs(pdgIdsHdecay[i])==16){
			return lvecsHdecay[i];
		}
	}
	std::cout << "IMPOSSIBLE: " << pdgIdsHdecay[0] << "," << pdgIdsHdecay[1] << "," << pdgIdsHdecay[2] << "," << pdgIdsHdecay[3] << std::endl;
 	return PtEtaPhiMVector(0,0,0,0);
}

PtEtaPhiMVector GenerateLorentzVectorLeptonTrue(vector<PtEtaPhiMVector> lvecsHdecay, vector<int> pdgIdsHdecay){
	for(int i=0; i<pdgIdsHdecay.size(); i++){
		if(abs(pdgIdsHdecay[i])==11 || abs(pdgIdsHdecay[i])==13 || abs(pdgIdsHdecay[i])==15){
			return lvecsHdecay[i];
		}
	}
	return PtEtaPhiMVector(0,0,0,0);
}

PtEtaPhiMVector CombineTwoPtEtaPhiM(PtEtaPhiMVector v1, PtEtaPhiMVector v2){return v1+v2;}




int SignateAbsLeptonPdgid(int pdgid_11, int pdgid_12, int pdgid_21, int pdgid_22){
	// electron
	if( abs(pdgid_11)==11 || abs(pdgid_12)==11 || abs(pdgid_21)==11 || abs(pdgid_22)==11 )
	{
		if( abs(pdgid_11)==12 || abs(pdgid_12)==12 || abs(pdgid_21)==12 || abs(pdgid_22)==12 )
		{
			return 11;
		}
		return -1; //mismatched flavours lepton <-> neutrino 
	}

	// muon
	if( abs(pdgid_11)==13 || abs(pdgid_12)==13 || abs(pdgid_21)==13 || abs(pdgid_22)==13 )
	{
		if( abs(pdgid_11)==14 || abs(pdgid_12)==14 || abs(pdgid_21)==14 || abs(pdgid_22)==14 )
		{
			return 13;
		}
		return -1; //mismatched flavours lepton <-> neutrino 
	}

	// tau
	if( abs(pdgid_11)==15 || abs(pdgid_12)==15 || abs(pdgid_21)==15 || abs(pdgid_22)==15 )
	{
		if( abs(pdgid_11)==16 || abs(pdgid_12)==16 || abs(pdgid_21)==16 || abs(pdgid_22)==16 )
		{
			return 15;
		}
		return -1; //mismatched flavours lepton <-> neutrino 
	}

	return 0; // no lepton found	
}



int SignateHiggsDecay(int pdgid_h1, int pdgid_h2, int pdgid_w11, int pdgid_w12, int pdgid_w21, int pdgid_w22)
{	
	// check for H->WW signature
	if( abs(pdgid_h1)!=24 || abs(pdgid_h2)!=24 )
	{
		return -3;
	}

	// check for H->WW->lvqq signature 
	if( abs(pdgid_w11)==11 || abs(pdgid_w11)==13 || abs(pdgid_w11)==15 )
	{
		// condition: partner must be neutrino
		if( abs(pdgid_w12)!=abs(pdgid_w11)+1 )
		{
			return -1; // mismatched lepton-neutrino flavours
		}

		// condition: other pair must be hadronic 
		if( abs(pdgid_w21)<1 || abs(pdgid_w21)>6 )
		{
			return -2; // no hadronic decays found
		}
		if( abs(pdgid_w22)<1 || abs(pdgid_w22)>6 )
		{
			return -2; // no hadronic decays found
		}

		// all conditions fulfilled
		return 1;
	}

	// check for H->WW->vlqq signature 
	if( abs(pdgid_w12)==11 || abs(pdgid_w12)==13 || abs(pdgid_w12)==15 )
	{
		// condition: partner must be neutrino
		if( abs(pdgid_w11)!=abs(pdgid_w12)+1 )
		{
			return -1; // mismatched lepton-neutrino flavours
		}

		// condition: other pair must be hadronic 
		if( abs(pdgid_w21)<1 || abs(pdgid_w21)>6 )
		{
			return -2; // no hadronic decays found
		}
		if( abs(pdgid_w22)<1 || abs(pdgid_w22)>6 )
		{
			return -2; // no hadronic decays found
		}

		// all conditions fulfilled
		return 2;
	}

	// check for H->WW->qqlv signature 
	if( abs(pdgid_w21)==11 || abs(pdgid_w21)==13 || abs(pdgid_w21)==15 )
	{
		// condition: partner must be neutrino
		if( abs(pdgid_w22)!=abs(pdgid_w21)+1 )
		{
			return -1; // mismatched lepton-neutrino flavours
		}

		// condition: other pair must be hadronic 
		if( abs(pdgid_w11)<1 || abs(pdgid_w11)>6 )
		{
			return -2; // no hadronic decays found
		}
		if( abs(pdgid_w12)<1 || abs(pdgid_w12)>6 )
		{
			return -2; // no hadronic decays found
		}

		// all conditions fulfilled
		return 3;
	}

	// check for H->WW->qqvl signature 
	if( abs(pdgid_w22)==11 || abs(pdgid_w22)==13 || abs(pdgid_w22)==15 )
	{
		// condition: partner must be neutrino
		if( abs(pdgid_w21)!=abs(pdgid_w22)+1 )
		{
			return -1; // mismatched lepton-neutrino flavours
		}

		// condition: other pair must be hadronic 
		if( abs(pdgid_w11)<1 || abs(pdgid_w11)>6 )
		{
			return -2; // no hadronic decays found
		}
		if( abs(pdgid_w12)<1 || abs(pdgid_w12)>6 )
		{
			return -2; // no hadronic decays found
		}

		// all conditions fulfilled
		return 4;
	}

	// H->WW->qqqq signature 
	return 0;
}


int SignateOnshellHadronicWboson(int signature_higgs_decay, float m_true_w1, float m_true_w2){
	// w1 decays hadronically
	if(signature_higgs_decay==3 || signature_higgs_decay==4){
		if( abs(m_true_w1-PDG_MASS_WBOSON)<THRESHOLD_ONSHELL_DEFINITION )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	// w2 decays hadronically
	if(signature_higgs_decay==1 || signature_higgs_decay==2){
		if( abs(m_true_w2-PDG_MASS_WBOSON)<THRESHOLD_ONSHELL_DEFINITION )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	return -1;
}


float ExtractTrueInformationLeptonFromHiggs(int signatureHiggsDecay, float info_w11, float info_w12, float info_w21, float info_w22)
{
	// H->WW->lvqq signature
	if(signatureHiggsDecay==1)
	{
		return info_w11;
	}
	
	// H->WW->vlqq signature
	if(signatureHiggsDecay==2)
	{
		return info_w12;
	}
	
	// H->WW->qqlv signature
	if(signatureHiggsDecay==3)
	{
		return info_w21;
	}

	// H->WW->qqvl signature
	if(signatureHiggsDecay==4)
	{
		return info_w22;
	}

	// no lepton from Higgs decays
	return -1;
}


PtEtaPhiMVector CombineTrueWHadFromHiggs(int signatureHiggsDecay, PtEtaPhiMVector lvec_w11, PtEtaPhiMVector lvec_w12, PtEtaPhiMVector lvec_w21, PtEtaPhiMVector lvec_w22)
{
	// H->WW->lvqq or H->WW->vlqq signature
	if(signatureHiggsDecay==1 || signatureHiggsDecay==2)
	{
		return lvec_w21 + lvec_w22;
	}

	// H->WW->lvqq or H->WW->vlqq signature
	if(signatureHiggsDecay==3 || signatureHiggsDecay==4)
	{
		return lvec_w11 + lvec_w12;
	}

	// no fully hadronic decay from Higgs decays
	return PtEtaPhiMVector(-999,-999,-999,-999);
}

PtEtaPhiMVector CombineTrueWLepFromHiggs(int signatureHiggsDecay, PtEtaPhiMVector lvec_w11, PtEtaPhiMVector lvec_w12, PtEtaPhiMVector lvec_w21, PtEtaPhiMVector lvec_w22)
{
	// H->WW->lvqq or H->WW->vlqq signature
	if(signatureHiggsDecay==1 || signatureHiggsDecay==2)
	{
		return lvec_w11 + lvec_w12;
	}

	// H->WW->lvqq or H->WW->vlqq signature
	if(signatureHiggsDecay==3 || signatureHiggsDecay==4)
	{
		return lvec_w21 + lvec_w22;
	}

	// no fully hadronic decay from Higgs decays
	return PtEtaPhiMVector(-999,-999,-999,-999);
}




// ==========  MAIN  ==========
// ============================
int main(){
	// ==========  SETUP
	// =================	
	// setup TChain
	cout << "Start: setup TChain" << endl;
	TChain rRecoChain("reco");
	TChain rTruthChain("truth");

	for(auto input_file_name:INPUT_FILE_NAMES){
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);
		rRecoChain.Add(file_path.c_str());
		rTruthChain.Add(file_path.c_str());
	}


	// index TruthChain to kick out un-matched events - then declare friends
	rTruthChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rRecoChain.AddFriend(&rTruthChain);


	// setup RDataFrame
	cout << "Start: setup RDataFrame" << endl;
	auto rDataFrame = RDataFrame(rRecoChain);
	auto rLoopManager = rDataFrame.Range(0); // limit input for testing


	// check if chains are matched properly
	auto nTotalEvents = rLoopManager.Count();

	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0){ 
		cout << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!" << endl;
		exit(1);
	}


	// generate jet lorentz vectors
	rLoopManager = rLoopManager.Define(
		"lvecs_jets", 
		GenerateJetLvecs, 
		{"jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"number_of_jets", 
		GetNumberOfJets, 
		{"lvecs_jets"}
	);


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

	// find hadronic H decay and generate turth obj lorentz vectors
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay1_from_H", 
		GenerateLorentzVectorMHiggsDecision, 
		{"truth.Tth_MC_Higgs_decay1_from_decay1_pt", "truth.Tth_MC_Higgs_decay1_from_decay1_eta", "truth.Tth_MC_Higgs_decay1_from_decay1_phi", "truth.Tth_MC_Higgs_decay1_from_decay1_m", "truth.Tth_MC_Higgs_decay1_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay1_from_decay2_pt", "truth.Tth_MC_Higgs_decay1_from_decay2_eta", "truth.Tth_MC_Higgs_decay1_from_decay2_phi", "truth.Tth_MC_Higgs_decay1_from_decay2_m","truth.Tth_MC_Higgs_decay1_from_decay2_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Wdecay2_from_H", 
		GenerateLorentzVectorMHiggsDecision, 
		{"truth.Tth_MC_Higgs_decay2_from_decay1_pt", "truth.Tth_MC_Higgs_decay2_from_decay1_eta", "truth.Tth_MC_Higgs_decay2_from_decay1_phi", "truth.Tth_MC_Higgs_decay2_from_decay1_m","truth.Tth_MC_Higgs_decay2_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay2_pt", "truth.Tth_MC_Higgs_decay2_from_decay2_eta", "truth.Tth_MC_Higgs_decay2_from_decay2_phi", "truth.Tth_MC_Higgs_decay2_from_decay2_m", "truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);


	// generate truth obj vector for each event for easier access
	rLoopManager = rLoopManager.Define(
		"truth_lvecs",
		GenerateTruthLvecs,
		{"lvec_b_from_t", "lvec_b_from_tbar", "lvec_Wdecay1_from_t", "lvec_Wdecay2_from_t", "lvec_Wdecay1_from_tbar", "lvec_Wdecay2_from_tbar", "lvec_Wdecay1_from_H", "lvec_Wdecay2_from_H"}
	);



	// TRUTH 

	// signatures
	rLoopManager = rLoopManager.Define(
		"signature_higgs_decay",
		SignateHiggsDecay,
		{"Tth_MC_Higgs_decay1_pdgId", "Tth_MC_Higgs_decay2_pdgId", "Tth_MC_Higgs_decay1_from_decay1_pdgId", "Tth_MC_Higgs_decay2_from_decay1_pdgId", "Tth_MC_Higgs_decay1_from_decay2_pdgId", "Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);
	
	rLoopManager = rLoopManager.Define(
		"signature_onshell_whad",
		SignateOnshellHadronicWboson,
		{"signature_higgs_decay", "Tth_MC_Higgs_decay1_m", "Tth_MC_Higgs_decay2_m"}
	);
	
	// Higgs decays decay products
	rLoopManager = rLoopManager.Define(
		"true_w11_lvec",
		GenerateLorentzVectorM,
		{"Tth_MC_Higgs_decay1_from_decay1_pt", "Tth_MC_Higgs_decay1_from_decay1_eta", "Tth_MC_Higgs_decay1_from_decay1_phi", "Tth_MC_Higgs_decay1_from_decay1_m"}
	);
	rLoopManager = rLoopManager.Define(
		"true_w12_lvec",
		GenerateLorentzVectorM,
		{"Tth_MC_Higgs_decay2_from_decay1_pt", "Tth_MC_Higgs_decay2_from_decay1_eta", "Tth_MC_Higgs_decay2_from_decay1_phi", "Tth_MC_Higgs_decay2_from_decay1_m"}
	);
	rLoopManager = rLoopManager.Define(
		"true_w21_lvec",
		GenerateLorentzVectorM,
		{"Tth_MC_Higgs_decay1_from_decay2_pt", "Tth_MC_Higgs_decay1_from_decay2_eta", "Tth_MC_Higgs_decay1_from_decay2_phi", "Tth_MC_Higgs_decay1_from_decay2_m"}
	);
	rLoopManager = rLoopManager.Define(
		"true_w22_lvec",
		GenerateLorentzVectorM,
		{"Tth_MC_Higgs_decay2_from_decay2_pt", "Tth_MC_Higgs_decay2_from_decay2_eta", "Tth_MC_Higgs_decay2_from_decay2_phi", "Tth_MC_Higgs_decay2_from_decay2_m"}
	);

	// Lepton from Higgs decays
	rLoopManager = rLoopManager.Define(
		"true_lepton_pt",
		ExtractTrueInformationLeptonFromHiggs,
		{"signature_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_pt", "Tth_MC_Higgs_decay2_from_decay1_pt", "Tth_MC_Higgs_decay1_from_decay2_pt", "Tth_MC_Higgs_decay2_from_decay2_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_eta",
		ExtractTrueInformationLeptonFromHiggs,
		{"signature_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_eta", "Tth_MC_Higgs_decay2_from_decay1_eta", "Tth_MC_Higgs_decay1_from_decay2_eta", "Tth_MC_Higgs_decay2_from_decay2_eta"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_phi",
		ExtractTrueInformationLeptonFromHiggs,
		{"signature_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_phi", "Tth_MC_Higgs_decay2_from_decay1_phi", "Tth_MC_Higgs_decay1_from_decay2_phi", "Tth_MC_Higgs_decay2_from_decay2_phi"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_m",
		ExtractTrueInformationLeptonFromHiggs,
		{"signature_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_m", "Tth_MC_Higgs_decay2_from_decay1_m", "Tth_MC_Higgs_decay1_from_decay2_m", "Tth_MC_Higgs_decay2_from_decay2_m"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_lvec",
		GenerateLorentzVectorM,
		{"true_lepton_pt", "true_lepton_eta", "true_lepton_phi", "true_lepton_m"}
	);

	// Combein true W decays from Higgs
	rLoopManager = rLoopManager.Define(
		"true_whad_lvec",
		CombineTrueWHadFromHiggs,
		{"signature_higgs_decay", "true_w11_lvec", "true_w12_lvec", "true_w21_lvec", "true_w22_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_wlep_lvec",
		CombineTrueWLepFromHiggs,
		{"signature_higgs_decay", "true_w11_lvec", "true_w12_lvec", "true_w21_lvec", "true_w22_lvec"}
	);




	// match jets to objs
	rLoopManager = rLoopManager.Define(
		"jet_potential_match_mask", 
		GenerateJetPotentialMatchMasks, 
		{"lvecs_jets", "truth_lvecs"}
	);
	rLoopManager = rLoopManager.Define(
		"jet_final_match_mask", 
		GenerateJetFinalMatchMasks, 
		{"lvecs_jets", "jet_potential_match_mask", "truth_lvecs"}
	);
	

	// collect jet indicies for SPANet training
	rLoopManager = rLoopManager.Define(
		"jet_to_object_indicies_fixed",
		CollectJetToObjectIndiciesFixed,
		{"jet_final_match_mask"}
	);

	
	// reconstruct objects from matched jets
	rLoopManager = rLoopManager.Define(
		"reconstructed_t_lvec",
		ReconstructT,
		{"lvecs_jets", "jet_final_match_mask"}
	);
	
	rLoopManager = rLoopManager.Define(
		"reconstructed_tbar_lvec",
		ReconstructTBar,
		{"lvecs_jets", "jet_final_match_mask"}
	);
	
	rLoopManager = rLoopManager.Define(
		"reconstructed_W_from_t_lvec",
		ReconstructWFromT,
		{"lvecs_jets", "jet_final_match_mask"}
	);
	
	rLoopManager = rLoopManager.Define(
		"reconstructed_W_from_tbar_lvec",
		ReconstructWFromTBar,
		{"lvecs_jets", "jet_final_match_mask"}
	);

	rLoopManager = rLoopManager.Define(
		"reconstructed_HW_lvec",
		ReconstructHW,
		{"lvecs_jets", "jet_final_match_mask"}
	); 



	


	// rename truth trees
	rLoopManager = rLoopManager.Define(
		"truth_t_m",
		RenameFloat,
		{"truth.Tth_MC_t_afterFSR_m"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_tbar_m",
		RenameFloat,
		{"truth.Tth_MC_tbar_afterFSR_m"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_W_from_t_m",
		RenameFloat,
		{"truth.Tth_MC_W_from_t_m"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_W_from_tbar_m",
		RenameFloat,
		{"truth.Tth_MC_W_from_tbar_m"}
	);
	
	rLoopManager = rLoopManager.Define(
		"truth_t_pt",
		RenameFloat,
		{"truth.Tth_MC_t_afterFSR_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_tbar_pt",
		RenameFloat,
		{"truth.Tth_MC_tbar_afterFSR_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_W_from_t_pt",
		RenameFloat,
		{"truth.Tth_MC_W_from_t_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"truth_W_from_tbar_pt",
		RenameFloat,
		{"truth.Tth_MC_W_from_tbar_pt"}
	);


	rLoopManager = rLoopManager.Define(
		"reco_met_value",
		RenameFloat,
		{"met_met_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_met_phi",
		RenameFloat,
		{"met_phi_NOSYS"}
	);


	// get higgs information
	rLoopManager = rLoopManager.Define(
		"higgs_decay_mode_custom",
		GenerateHiggsDecayModeCustom,
		{"truth.Tth_MC_Higgs_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_pdgId"}
	);
	
	rLoopManager = rLoopManager.Define(
		"higgs_decay_decay_mode",
		GenerateHiggsDecayDecayMode,
		{"higgs_decay_mode_custom", "truth.Tth_MC_Higgs_decay1_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay1_pdgId", "truth.Tth_MC_Higgs_decay1_from_decay2_pdgId", "truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);

	
	//>> filter pdgIds (exclude invalid ids and take absolute value)
	rLoopManager = rLoopManager.Define(
		"pdgid_Hdecay11_filtered",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay1_from_decay1_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"pdgid_Hdecay12_filtered",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay2_from_decay1_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"pdgid_Hdecay21_filtered",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay1_from_decay2_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"pdgid_Hdecay22_filtered",
		GetFilteredPdgIDs,
		{"truth.Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);

	rLoopManager = rLoopManager.Define(
		"pdgids_Hdecay_ordered_filtered",
		CombineHiggsDecayPDGIDs,
		{"pdgid_Hdecay11_filtered", "pdgid_Hdecay12_filtered", "pdgid_Hdecay21_filtered", "pdgid_Hdecay22_filtered"}
	);

	
	rLoopManager = rLoopManager.Define(
		"lvec_Hdecay11",
		GenerateLorentzVectorM,
		{"truth.Tth_MC_Higgs_decay1_from_decay1_pt", "truth.Tth_MC_Higgs_decay1_from_decay1_eta", "truth.Tth_MC_Higgs_decay1_from_decay1_phi", "truth.Tth_MC_Higgs_decay1_from_decay1_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Hdecay12",
		GenerateLorentzVectorM,
		{"truth.Tth_MC_Higgs_decay2_from_decay1_pt", "truth.Tth_MC_Higgs_decay2_from_decay1_eta", "truth.Tth_MC_Higgs_decay2_from_decay1_phi", "truth.Tth_MC_Higgs_decay2_from_decay1_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Hdecay21",
		GenerateLorentzVectorM,
		{"truth.Tth_MC_Higgs_decay1_from_decay2_pt", "truth.Tth_MC_Higgs_decay1_from_decay2_eta", "truth.Tth_MC_Higgs_decay1_from_decay2_phi", "truth.Tth_MC_Higgs_decay1_from_decay2_m"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Hdecay22",
		GenerateLorentzVectorM,
		{"truth.Tth_MC_Higgs_decay2_from_decay2_pt", "truth.Tth_MC_Higgs_decay2_from_decay2_eta", "truth.Tth_MC_Higgs_decay2_from_decay2_phi", "truth.Tth_MC_Higgs_decay2_from_decay2_m"}
	);

	rLoopManager = rLoopManager.Define(
		"lvecs_Hdecay_ordered",
		CombineHiggsDecayLorentzVectors,
		{"lvec_Hdecay11", "lvec_Hdecay12", "lvec_Hdecay21", "lvec_Hdecay22"}
	);





	//>> Generate Neutrino Information
	rLoopManager = rLoopManager.Define(
		"true_neutrino_lvec",
		GenerateLorentzVectorNeutrinoTrue,
		{"lvecs_Hdecay_ordered", "pdgids_Hdecay_ordered_filtered"}
	);

	//>> Generate Lepton Information
	rLoopManager = rLoopManager.Define(
		"classifier_lepton_flavour",
		ClassifyRecoLeptonFlavour,
		{"pass_ejets_NOSYS", "pass_mujets_NOSYS"}
	);


	rLoopManager = rLoopManager.Define(
		"signature_abs_lepton_pdgid",
		SignateAbsLeptonPdgid,
		{"Tth_MC_Higgs_decay1_from_decay1_pdgId", "Tth_MC_Higgs_decay2_from_decay1_pdgId", "Tth_MC_Higgs_decay1_from_decay2_pdgId", "Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);





	






	//>> get correct lepton information
	rLoopManager = rLoopManager.Define(
		"reco_lepton_pt",
		ExtractRecoInformationLepton,
		{"classifier_lepton_flavour", "el_pt_NOSYS", "el_select_loose_NOSYS", "mu_pt_NOSYS", "mu_select_tight_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_lepton_eta",
		ExtractRecoInformationLepton,
		{"classifier_lepton_flavour", "el_eta", "el_select_loose_NOSYS", "mu_eta", "mu_select_tight_NOSYS"}
	);

	rLoopManager = rLoopManager.Define(
		"reco_lepton_phi",
		ExtractRecoInformationLepton,
		{"classifier_lepton_flavour", "el_phi", "el_select_loose_NOSYS", "mu_phi", "mu_select_tight_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_lepton_e",
		ExtractRecoInformationLepton,
		{"classifier_lepton_flavour", "el_e_NOSYS", "el_select_loose_NOSYS", "mu_e_NOSYS", "mu_select_tight_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_lepton_lvec",
		GenerateLorentzVectorE,
		{"reco_lepton_pt", "reco_lepton_eta", "reco_lepton_phi", "reco_lepton_e"}
	);






	rLoopManager = rLoopManager.Define(
		"classifier_event_status",
		ClassifyEventStatus,
		{"jet_to_object_indicies_fixed"}
	);




	// apply filter	and limit if needed
	auto rLoopManagerFiltered = rLoopManager.Filter( 
		"(number_of_jets>=8) && (signature_higgs_decay>0) && (signature_onshell_whad==1) && (signature_abs_lepton_pdgid==11 || signature_abs_lepton_pdgid==13)" 
	).Range(MAX_NUMBER_OF_EVENTS);

	cout << rLoopManagerFiltered.Count().GetValue() << " events passed the selection!" << endl;

	// save snapshot to disk
	cout << "Start: saving snapshot" << endl;
	rLoopManagerFiltered.Snapshot(
		"matched", 
		OUTPUT_PATH+OUTPUT_FILE,
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
		if(currentDeltaR<=THRESHOLD_DELTA_R)
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
			if( (jetFinalMatchMasks[j] & 1<<i)!=0 ){ //final match is 1 for current TRUTH_PARTON
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
