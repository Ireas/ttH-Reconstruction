#include<filesystem>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
using namespace std;
#include"TROOT.h"
#include"TFile.h"
#include"TLorentzVector.h"
#include"TChain.h"
#include"TTree.h"
#include"TH1.h"
#include"TKey.h"
#include"TClass.h"
#include"ROOT/RDataFrame.hxx"
#include<Math/VectorUtil.h> // for DeltaR
#include<Math/Vector4D.h> // for PtEtaPhiEVector and PtEtaPhiMVector
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;

#include<chrono> // for time measurements during the program



// ==========  CONSTANTS  ==========
// =================================
// Event Filter String
const string FILTER = "(number_of_jets>=0)";// && classification_true_higgs_decay==-1 && classification_true_t1_decay==1 && classification_true_t2_decay==1";
//&& ( signature_higgs_decay<0 || classification_onshell_whad==1 )

// Paths
const string INPUT_PATH = "/media/ireas/Data/download/";
const string OUTPUT_PATH = "/media/ireas/Data/v6/matched/all_8+j/";
	
// branches
const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
	// logistics
	"eventNumber",
	
	// global event information
	"number_of_jets",
	"number_of_bjets",
	"number_of_leptons",

	// reco values
	"jet_pt_NOSYS", // jets
	"jet_eta",
	"jet_phi",
	"jet_e_NOSYS",
	"jet_btag_continous", // b-tagging
	"jet_btag_85wp",
	"jet_btag_77wp",
	"jet_btag_70wp",
	"jet_btag_60wp",
	"reco_lepton_pt", // lepton
	"reco_lepton_eta",
	"reco_lepton_phi",
	"reco_lepton_e",
	"reco_met_value", // met
	"reco_met_phi",
};



// ==========  FUNCTION DECLARATION  ==========
// ============================================
// get number of jets
int GetNumberOfJets(vector<PtEtaPhiEVector> jetLvecs){return jetLvecs.size();}
int GetNumberOfBJets(vector<char> jet_btagging)
{
	int btags = 0;
	for(char tag:jet_btagging)
	{
		if((int)tag==1)
		{
			btags++;
		}
	}

	return btags;
}
int GetNumberOfLeptons(char el_select_loose_NOSYS, char muon_select_loose_NOSYS)
{
	int leptons = 0;

	// check for electron 
	if((int)el_select_loose_NOSYS==1)
	{
		leptons++;
	}

	// check for muons 
	if((int)muon_select_loose_NOSYS==1)
	{
		leptons++;
	}

	return leptons;
}

// rename variables for new tree
float RenameFloat(float target){return target;}
int RenameInt(int target){return target;}
vector<int> RenameVectorInt(vector<int> target){return target;}

// convert char to int
int ConvertCharToInt(char target){return (int)target;}
vector<int> ConvertVectorCharToVectorInt(vector<char> target){
	vector<int> new_target;
	for(char element:target)
	{
		new_target.push_back((int)element);
	}
	return new_target;
}

// generate lorentz vector for truth object
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy){return PtEtaPhiEVector(pt,eta,phi,energy);}
float ExtractPt(PtEtaPhiEVector lvec){return lvec.Pt();}
float ExtractEta(PtEtaPhiEVector lvec){return lvec.Eta();}
float ExtractPhi(PtEtaPhiEVector lvec){return lvec.Phi();}
float ExtractE(PtEtaPhiEVector lvec){return lvec.E();}

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

int ClassifyRecoLeptonFlavour(char passElectronChar, char passMuonChar)
{
	bool passElectron = (bool)passElectronChar;
	bool passMuon = (bool)passMuonChar;

	// both
	if( passElectron && passMuon )
	{
		std::cout << "help me ?" << std::endl;
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

PtEtaPhiMVector CombineTwoPtEtaPhiM(PtEtaPhiMVector v1, PtEtaPhiMVector v2){return v1+v2;}







// Function to split a string by a delimiter and return the parts as a vector
std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to extract DSID
std::string extractDSID(const std::string &filename) {
    std::vector<std::string> parts = split(filename, '.');
    if (parts.size() > 2) {
        return parts[2]; // DSID is the third element (index 2)
    } else {
        return ""; // Return an empty string if the format is not as expected
    }
}


// ==========  MAIN  ==========
// ===========================
const int LAST_ENTRY_IN_GRP15 = 186;
const int LAST_ENTRY_IN_GRP16 = 924;
const int LAST_ENTRY_IN_GRP17 = 1043;
const int LAST_ENTRY_IN_GRP18 = 904;

std::map<int,int> GROUP_NUMBERS = {
	{15, 40043329},
	{16, 40043331},
	{17, 40043332},
	{18, 40043333}
};


int Prepare(string input_file);


string GetCorrectFilePath(int group, int entry)
{
	string file_path = "/media/ireas/Data/data/grp"; 
	file_path.append(to_string(group));  
	file_path.append("/user.chscheul.");
	file_path.append( to_string(GROUP_NUMBERS[group]) );
	file_path.append("._");
	if(entry<10)
	{
		file_path.append("00000" + to_string(entry) );
	}
	else if(entry<100)
	{
		file_path.append("0000" + to_string(entry) );
	}
	else if(entry<1000)
	{
		file_path.append("000" + to_string(entry) );
	}
	else if(entry<10000)
	{
		file_path.append("00" + to_string(entry) );
	}
	else if(entry<100000)
	{
		file_path.append("0" + to_string(entry) );
	}
	else
	{
		file_path.append("" + to_string(entry) );
	}
	file_path.append(".output.root");
	return file_path;
}

int main(int argc, char** argv)
{
	std::vector<string> valid_files;

	for(int i=0; i<LAST_ENTRY_IN_GRP15; i++)
	{
		string file_path = GetCorrectFilePath(15,i+1);
		std::cout << file_path << std::endl;
		if( std::filesystem::exists(file_path) )
		{ // file exists
			TFile *file = TFile::Open(file_path.c_str());
    			if (!file || file->IsZombie()) {
        		std::cout << "  Error: file not found!" << std::endl << std::endl;
        		continue;
    		}
    		// Get the TTree object
    		TTree *tree = (TTree*)file->Get("reco"); // Replace "tree_name" with the actual tree name in the file
    		if (!tree) {
    		    std::cout << "  Error: reco tree not found!" << std::endl << std::endl;
    		    continue;
    		}
    		
    		TObjArray *subBranches = tree->GetListOfBranches();
			if(subBranches && subBranches->GetEntries() > 0)
			{
				valid_files.push_back(file_path);
    		} 
			else 
			{
        		std::cout << "'reco' branch does not have any sub-branches." << std::endl;
				continue;
    		}

		}
		else
		{ // file is missing
			std::cout << "  Error: file is missing!" << std::endl << std::endl;
			continue;
		}
        
		// success
		std::cout << "  Perfect!" << std::endl << std::endl;
	}



	// Create a TChain to combine the trees
    TChain chain("reco"); // Replace "tree_name" with the actual name of the TTree

    // Add all input files to the chain
    for (const auto& file : valid_files) {
        chain.Add(file.c_str());
    }

    // Create the output file
	string outputFile = "/media/ireas/Data/data/temp.root";
    TFile outFile(outputFile.c_str(), "RECREATE");

    // Clone the structure of the input TTree into the output TTree
    TTree* mergedTree = chain.CloneTree(-1, "fast");

    // Write the merged TTree into the output file
    mergedTree->Write();

    // Close the output file
    outFile.Close();

    std::cout << "Merged " << valid_files.size() << " files into " << outputFile << std::endl;

	Prepare(outputFile);


	return 0;
}





int Prepare(string input_file)
{
	// SETUP
	// ==============================	
	// setup TChain
	cout << " > setup TChain" << endl;
	TChain rRecoChain("reco");

	// link input file
	cout << " > preparing " << input_file << endl;
	auto file_path = std::string();
	file_path.append(input_file);
	rRecoChain.Add(file_path.c_str());

	// index TruthChain to kick out un-matched events - then declare friends
	rRecoChain.BuildIndex("eventNumber");  // just for security, use DSID too

	// setup RDataFrame
	cout << " > setup RDataFrame" << endl;
	auto rDataFrame = RDataFrame(rRecoChain);
	auto rLoopManager = rDataFrame.Range(0); // no limit on input

	// check if chains are matched properly
	auto nTotalEvents = rLoopManager.Count();

	// JETS
	// ==============================
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
	rLoopManager = rLoopManager.Define(
		"number_of_bjets",
		GetNumberOfBJets, 
		{"jet_DL1dv01_FixedCutBEff_85_select"}
	);
	rLoopManager = rLoopManager.Define(
		"number_of_leptons",
		GetNumberOfLeptons, 
		{"pass_ejets_NOSYS", "pass_mujets_NOSYS"}
	);



	// Lepton
	// ==============================
	// identification
	rLoopManager = rLoopManager.Define(
		"classifier_lepton_flavour",
		ClassifyRecoLeptonFlavour,
		{"pass_ejets_NOSYS", "pass_mujets_NOSYS"}
	);

	// get correct lepton information
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




	// RENAMING
	// ==============================
	// b-tagging
	rLoopManager = rLoopManager.Define(
		"jet_btag_continous",
		RenameVectorInt,
		{"jet_DL1dv01_Continuous_quantile"}
	);
	rLoopManager = rLoopManager.Define(
		"jet_btag_60wp",
		ConvertVectorCharToVectorInt,
		{"jet_DL1dv01_FixedCutBEff_60_select"}
	);
	rLoopManager = rLoopManager.Define(
		"jet_btag_70wp",
		ConvertVectorCharToVectorInt,
		{"jet_DL1dv01_FixedCutBEff_70_select"}
	);
	rLoopManager = rLoopManager.Define(
		"jet_btag_77wp",
		ConvertVectorCharToVectorInt,
		{"jet_DL1dv01_FixedCutBEff_77_select"}
	);
	rLoopManager = rLoopManager.Define(
		"jet_btag_85wp",
		ConvertVectorCharToVectorInt,
		{"jet_DL1dv01_FixedCutBEff_85_select"}
	);


	// reco met values
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


	// FINALISE
	// ==============================
	// apply filter	and limit if needed
	auto rLoopManagerFiltered = rLoopManager.Filter(FILTER);
	cout << " > " << rLoopManagerFiltered.Count().GetValue() << " events passed the following pre-selection: " << FILTER << endl;


	// save snapshot to disk
	cout << " > saving snapshot" << endl;

	string outputFile = "/media/ireas/Data/data/temp2.root";
	rLoopManagerFiltered.Snapshot(
		"prepared", 
		outputFile,
		OUTPUT_COLOUMN_NAMES
	);
		
	return 0;
}
