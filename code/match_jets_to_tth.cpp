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
const int MAX_NUMBER_OF_EVENTS = 0; // set to 0 for no limit
const float PDG_MASS_WBOSON = 80.3692e3; // mass in MeV
const float THRESHOLD_DELTA_R = 0.4; // maximum reco jet deviation from the truth for matching
const float THRESHOLD_ONSHELL_DEFINITION = 1e3; // maximum deviation from DPG mass in MeV to classify as onshell

const float LUMINOSITY = 139 * 1e3; //convert fm⁻1 to pb^-1, value from Baptiste full ATLAS run 2 set
const int TRUNCATE = 50; // Take only every i-th event to reduce file size (0 to disable truncating)

// Event Filter String
const string FILTER = "(count_var==1) && (number_of_jets>=5) && (classification_event_channel>=0) && (number_of_matches>0)";// && classification_true_higgs_decay==-1 && classification_true_t1_decay==1 && classification_true_t2_decay==1";
const bool BREAK_AFTER_FIRST_FILE = false;
//&& ( signature_higgs_decay<0 || classification_onshell_whad==1 )


// Paths
const string INPUT_PATH = "/media/ireas/Data/download/";
const string OUTPUT_PATH = "/media/ireas/Data/v6/matched/all_5+j_truncated_50_ttHWW_amplified/";

// File path
std::string XSEC_PATH = "/media/ireas/Data/PMGxsecDB_mc16.txt";
	


const char* INPUT_FILE_NAMES[] = { // put into array for easier access	
	// PowhegPythia-ttH (125 GeV, allhad) 
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000001.output.root", 
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000002.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000003.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043369._000004.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000001.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000002.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000003.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000004.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000005.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000006.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000007.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043372._000008.output.root",
	"user.chscheul.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043366._000001.output.root",
	// PowhegPythia-ttH (125 GeV, semilep) 
	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043370._000001.output.root",
	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043373._000001.output.root",
	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043367._000001.output.root",
	"user.chscheul.346344.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043367._000002.output.root",
	// PowhegPythia-ttH (125 GeV, dilep) 
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043371._000001.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000001.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000002.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000003.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043374._000004.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043368._000001.output.root",
	"user.chscheul.346345.PhPy8EG.DAOD_PHYS.e7148_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043368._000002.output.root",
	// PowhegPythia-ttbar (nonallhad) 
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000001.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000002.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000003.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043337._000004.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000001.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000003.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000004.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000005.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000006.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000007.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000008.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000009.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000010.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000011.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043340._000012.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000001.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000002.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000003.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000004.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000005.output.root",
	"user.chscheul.410470.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043334._000006.output.root",
	// PowhegPythia-ttbar (allhad) 
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000001.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000004.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000005.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000006.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000008.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000009.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000010.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000013.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043338._000014.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043341._000001.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043341._000002.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043341._000003.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000001.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000002.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000003.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000004.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000005.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000006.output.root",
	"user.chscheul.410471.PhPy8EG.DAOD_PHYS.e6337_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043335._000007.output.root",
	// PowhegPythia-ttbar (dilep) 
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000001.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000002.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000003.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000004.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000005.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000006.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043339._000008.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000001.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000002.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000004.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000005.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000006.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000013.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000015.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000016.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000017.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000018.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000019.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000020.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000021.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000022.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000024.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000025.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000026.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000027.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000028.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000029.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000030.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000031.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000032.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000033.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000034.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000035.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000036.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000038.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000040.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000042.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000043.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000044.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000045.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000046.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000047.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000048.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043342._000049.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000001.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000002.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000003.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000004.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000005.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000007.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000008.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000009.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000010.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000011.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000012.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000013.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000014.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000015.output.root",
	"user.chscheul.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p6026.ttHWW-240620-v1_output/user.chscheul.40043336._000016.output.root",
	// ??? 
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000001.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000002.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000003.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000004.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000005.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000006.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000007.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000008.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000009.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000010.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000011.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000012.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000013.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000014.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000015.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000016.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000017.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000018.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000019.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000020.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043352._000021.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000001.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000002.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000003.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000004.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000005.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000006.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000007.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000008.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000009.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000013.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000014.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000015.output.root",
	"user.chscheul.411316.PowhegHerwig7EvtGen.DAOD_PHYS.e7765_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043355._000016.output.root",
	// ??? 
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000002.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000003.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000004.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000005.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000006.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000007.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000008.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000009.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000010.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043358._000011.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000001.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000002.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000003.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000004.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000005.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000006.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000007.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000008.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000009.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000010.output.root",
	"user.chscheul.700121.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043362._000011.output.root",
	// ??? 
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000001.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000002.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000003.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000004.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000005.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000006.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000007.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000008.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000009.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000010.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000011.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043359._000012.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000001.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000002.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000003.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000004.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000006.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000007.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000008.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000009.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000010.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000011.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000012.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000013.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000014.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000015.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000016.output.root",
	"user.chscheul.700122.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043363._000021.output.root",
	// ??? 
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000001.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000002.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000003.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000004.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000005.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000006.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000007.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043360._000008.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000004.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000005.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000006.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000007.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000008.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000009.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000011.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000012.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000013.output.root",
	"user.chscheul.700123.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043364._000014.output.root",
	// ??? 
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000001.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000002.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000003.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000004.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000005.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000006.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000007.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000008.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000009.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000010.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000011.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13144_p6026.ttHWW-240620-v1_output/user.chscheul.40043361._000012.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000002.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000003.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000004.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000005.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000006.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000007.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000009.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000010.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000011.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000012.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000013.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000014.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000015.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000016.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000017.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000018.output.root",
	"user.chscheul.700124.Sh.DAOD_PHYS.e8253_s3681_r13145_p6026.ttHWW-240620-v1_output/user.chscheul.40043365._000019.output.root",
};


const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
	// logistics
	"mcChannelNumber",
	"eventNumber",

	// weights
	"SM_event_xsecs",
	"SM_event_weight",
	
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
	"jet_final_match_mask",
	"reco_lepton_pt", // lepton
	"reco_lepton_eta",
	"reco_lepton_phi",
	"reco_lepton_e",
	"reco_met_value", // met
	"reco_met_phi",
	"reco_whad_pt",
	"reco_whad_eta",
	"reco_whad_phi",
	"reco_whad_e",

	// true values
	"true_lepton_pt", // lepton
	"true_lepton_eta",
	"true_lepton_phi",
	"true_lepton_m",
	"true_neutrino_pt", // neutrino
	"true_neutrino_eta",
	"true_neutrino_phi",
	"true_neutrino_m",
	"true_whad_pt", // hadronic w from H
	"true_whad_eta",
	"true_whad_phi",
	"true_whad_m",
	"true_wlep_pt", // leptonic W from H
	"true_wlep_eta",
	"true_wlep_phi",
	"true_wlep_m",

	// stuff
	"jet_to_object_indicies_fixed",

	// classificiation
	"classification_event_channel",
	"classification_true_t1_decay",
	"classification_true_t2_decay",
	"classification_true_higgs_decay",
	"classification_event_completion",
	"classification_onshell_whad",

	// true event signatures
	"signature_abs_lepton_pdgid",

	// evaluation
	"successful_matches",
};

// Dictionaries to store dataset_number: crossSection_pb and dataset_number: kFactor
int currentDSID = -1;
std::map<int, double> xSecs;
std::map<int, double> kFactors;
std::map<int, double> genFiltersEff;
std::map<int, double> sumOfWeights;

std::map<int, int> successfulMatchHistory;
std::map<int, int> successfulPossibleMatchHistory;
std::map<int, int> successfulMatchHistorySignal;
std::map<int, int> successfulPossibleMatchHistorySignal;


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


// ==========  FUNCTION DECLARATION  ==========
// ============================================


// Utilities
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
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass){return PtEtaPhiMVector(pt,eta,phi,mass);}
PtEtaPhiEVector GenerateLorentzVectorE(Float_t pt, Float_t eta, Float_t phi, Float_t energy){return PtEtaPhiEVector(pt,eta,phi,energy);}
PtEtaPhiMVector GenerateLorentzVectorMHiggsDecision(Float_t pt1, Float_t eta1, Float_t phi1, Float_t mass1, Int_t pdgId1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t mass2, Int_t pdgId2);
float ExtractPt(PtEtaPhiMVector lvec){return lvec.Pt();}
float ExtractEta(PtEtaPhiMVector lvec){return lvec.Eta();}
float ExtractPhi(PtEtaPhiMVector lvec){return lvec.Phi();}
float ExtractM(PtEtaPhiMVector lvec){return lvec.M();}

float ExtractPtE(PtEtaPhiEVector lvec){return lvec.Pt();}
float ExtractEtaE(PtEtaPhiEVector lvec){return lvec.Eta();}
float ExtractPhiE(PtEtaPhiEVector lvec){return lvec.Phi();}
float ExtractEE(PtEtaPhiEVector lvec){return lvec.E();}

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


int GetFilteredPdgIDs(Int_t pdgId);





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



vector<int> CombineHiggsDecayPDGIDs(int pdgIdHDecay11, int pdgIdHDecay12, int pdgIdHDecay21, int pdgIdHDecay22){return vector<int>{pdgIdHDecay11, pdgIdHDecay12, pdgIdHDecay21, pdgIdHDecay22};}
vector<PtEtaPhiMVector> CombineHiggsDecayLorentzVectors(PtEtaPhiMVector lvecHdecay11, PtEtaPhiMVector lvecHdecay12, PtEtaPhiMVector lvecHdecay21, PtEtaPhiMVector lvecHdecay22){return vector<PtEtaPhiMVector>{lvecHdecay11, lvecHdecay12, lvecHdecay21, lvecHdecay22};}


int ClassifyTrueTopDecay(int pdgID_Wq1, int pdgID_Wq2)
{
	// check if any particle is invalid
	if(abs(pdgID_Wq1)==0 || abs(pdgID_Wq2)==0)
	{
		return -1;
	}

	// check for hadronic decay qq
	if( (abs(pdgID_Wq1)>0 && abs(pdgID_Wq1)<7) && (abs(pdgID_Wq2)>0 && abs(pdgID_Wq2)<7) )
	{
		return 1;
	}
	
	// check for leptonic decay lv
	if( (abs(pdgID_Wq1)==11 && abs(pdgID_Wq2)==12) || (abs(pdgID_Wq1)==13 && abs(pdgID_Wq2)==14) || (abs(pdgID_Wq1)==15 && abs(pdgID_Wq2)==16) )
	{
		return 2;
	}
	// check for leptonic decay vl
	if( (abs(pdgID_Wq1)==12 && abs(pdgID_Wq2)==11) || (abs(pdgID_Wq1)==14 && abs(pdgID_Wq2)==13) || (abs(pdgID_Wq1)==16 && abs(pdgID_Wq2)==15) )
	{
		return 3;
	}

	// default case, should never reach this point
	return 0;
}

int ClassifyTrueHiggsDecay(int pdgID_H1, int pdgID_H2, int pdgID_H11, int pdgID_H12, int pdgID_H21, int pdgID_H22)
{
	// check if Higgs doesnt exists
	if(abs(pdgID_H1)==0 || abs(pdgID_H2)==0 || abs(pdgID_H2)>100 || abs(pdgID_H1)>100)
	{
		return -1;
	}

	// check for bb
	if(abs(pdgID_H1)==5 || abs(pdgID_H2)==5)
	{
		return 2;
	}
	// check for tautau
	if(abs(pdgID_H1)==15 || abs(pdgID_H2)==15)
	{
		return 3;
	}

	// check for WW
	if(abs(pdgID_H1)==24 || abs(pdgID_H2)==24)
	{
		// check if all particles are valid objects
		if( abs(pdgID_H11)==0 || abs(pdgID_H12)==0 || abs(pdgID_H21)==0 || abs(pdgID_H2)==0)
		{
			return -2;
		}

		// check for qqqq
		if( (abs(pdgID_H11)>0 && abs(pdgID_H11)<7) && (abs(pdgID_H12)>0 && abs(pdgID_H12)<6) && (abs(pdgID_H21)>0 && abs(pdgID_H21)<6) && (abs(pdgID_H22)>0 && abs(pdgID_H22)<6) )
		{
			return 10;
		}

		// check for qqlv
		if( (abs(pdgID_H11)>0 && abs(pdgID_H11)<7) && (abs(pdgID_H12)>0 && abs(pdgID_H12)<6) && ( (abs(pdgID_H21)==11 && abs(pdgID_H22)==12) || (abs(pdgID_H21)==13 && abs(pdgID_H22)==14) || (abs(pdgID_H21)==15 && abs(pdgID_H22)==16) ) )
		{
			return 11;
		}
		// check for qqvl
		if( (abs(pdgID_H11)>0 && abs(pdgID_H11)<7) && (abs(pdgID_H12)>0 && abs(pdgID_H12)<6) && ( (abs(pdgID_H21)==12 && abs(pdgID_H22)==11) || (abs(pdgID_H21)==14 && abs(pdgID_H22)==13) || (abs(pdgID_H21)==16 && abs(pdgID_H22)==15) ) )
		{
			return 12;
		}
		// check for lvqq
		if( ( (abs(pdgID_H11)==11 && abs(pdgID_H12)==12) || (abs(pdgID_H11)==13 && abs(pdgID_H12)==14) || (abs(pdgID_H11)==15 && abs(pdgID_H12)==16) ) && (abs(pdgID_H21)>0 && abs(pdgID_H21)<7) && (abs(pdgID_H22)>0 && abs(pdgID_H22)<6) )
		{
			return 13;
		}
		// check for vlqq
		if( ( (abs(pdgID_H11)==12 && abs(pdgID_H12)==11) || (abs(pdgID_H11)==14 && abs(pdgID_H12)==13) || (abs(pdgID_H11)==16 && abs(pdgID_H12)==15) ) && (abs(pdgID_H21)>0 && abs(pdgID_H21)<7) && (abs(pdgID_H22)>0 && abs(pdgID_H22)<6) )
		{
			return 14;
		}
	}

	// default case for stuff like H->gammagamma, H->cc
	return 0;
}

int ClassifyEventCompletion(int classifier_true_t1_decay, int classifier_true_t2_decay, int classifier_true_higgs_decay, vector<int> jet_to_object_indicies_fixed){
	// check for allhad ttbar and semileptonic H->WW
	if(classifier_true_t1_decay!=1 || classifier_true_t2_decay!=1 || (classifier_true_higgs_decay!=11 && classifier_true_higgs_decay!=12 && classifier_true_higgs_decay!=13 && classifier_true_higgs_decay!=14) )
	{
		return -1;
	}

	// check if any object has no jet
	for(int index:jet_to_object_indicies_fixed)
	{
		if(index==-1)
		{
			return -2;
		}
	}

	return 1;
}

int ClassifyEventChannel(
	int pgdid_t1_W_q1, int pgdid_t1_W_q2, int pgdid_t2_W_q1, int pgdid_t2_W_q2, // decay partcles from W from t 
	int pgdid_H_d1, int pgdid_H_d2, // decay particles from H
	int pgdid_H_d1_d1, int pgdid_H_d1_d2, int pgdid_H_d2_d1, int pgdid_H_d2_d2 // decay particles from decay from H
)
{
	// check proper t decay
	if( abs(pgdid_t1_W_q1)>100 || abs(pgdid_t1_W_q2)>100 || abs(pgdid_t2_W_q1)>100 || abs(pgdid_t2_W_q2)>100 )
	{
		return -1; // invalid
	}

	// check if H does not exits
	if( abs(pgdid_H_d1)>100 || abs(pgdid_H_d2)>100 )
	{
		return 1; // ttbar
	}

	// check for H->bb
	if( abs(pgdid_H_d1)==5 && abs(pgdid_H_d2)==5 )
	{
		return 2; // tt(H->bb)
	}
	// check for H->cc
	if( abs(pgdid_H_d1)==4 && abs(pgdid_H_d2)==4 )
	{
		return 3; // tt(H->cc)
	}
	// check for H->tautau
	if( abs(pgdid_H_d1)==15 && abs(pgdid_H_d2)==15 )
	{
		return 4; // tt(H->tautau)
	}
	// check for H->ZZ
	if( abs(pgdid_H_d1)==23 && abs(pgdid_H_d2)==23 )
	{
		return 5; // tt(H->ZZ)
	}
	// check for H->yy
	if( abs(pgdid_H_d1)==22 && abs(pgdid_H_d2)==22 )
	{
		return 6; // tt(H->yy)
	}

	// check for others not H->WW decays
	if( !(abs(pgdid_H_d1)==24 && abs(pgdid_H_d2)==24) )
	{
		return 9; // tt(H->OTHER)
	}


	// check for H->WW
	if( abs(pgdid_H_d1)==24 && abs(pgdid_H_d2)==24 )
	{
		//check for full hadronic decays
		if( 
			(abs(pgdid_H_d1_d1)>0 && abs(pgdid_H_d1_d1)<7) && 
			(abs(pgdid_H_d1_d2)>0 && abs(pgdid_H_d1_d2)<7) && 
			(abs(pgdid_H_d2_d1)>0 && abs(pgdid_H_d2_d1)<7) && 
			(abs(pgdid_H_d2_d2)>0 && abs(pgdid_H_d2_d2)<7) 
		)
		{
			return 10; // tt(H->qqqq)
		}
		
		//check for semileptonic decays
		if( 
			(abs(pgdid_H_d1_d1)>10 && abs(pgdid_H_d1_d1)<17) && 
			(abs(pgdid_H_d1_d2)>10 && abs(pgdid_H_d1_d2)<17) &&
			(abs(pgdid_H_d2_d1)>0 && abs(pgdid_H_d2_d1)<7) && 
			(abs(pgdid_H_d2_d2)>0 && abs(pgdid_H_d2_d2)<7)
		)
		{
			return 11; // tt(H->lvqq) including tau as l
		}
		if( 
			(abs(pgdid_H_d1_d1)>0 && abs(pgdid_H_d1_d1)<7) && 
			(abs(pgdid_H_d1_d2)>0 && abs(pgdid_H_d1_d2)<7) &&
			(abs(pgdid_H_d2_d1)>10 && abs(pgdid_H_d2_d1)<17) && 
			(abs(pgdid_H_d2_d2)>10 && abs(pgdid_H_d2_d2)<17)
		)
		{
			return 11; // tt(H->qqlv) including tau as l
		}

		// check for dileptonic decays
		if( 
			(abs(pgdid_H_d1_d1)>10 && abs(pgdid_H_d1_d1)<17) && 
			(abs(pgdid_H_d1_d2)>10 && abs(pgdid_H_d1_d2)<17) && 
			(abs(pgdid_H_d2_d1)>10 && abs(pgdid_H_d2_d1)<17) && 
			(abs(pgdid_H_d2_d2)>10 && abs(pgdid_H_d2_d2)<17) 
		)
		{
			return 12; // tt(H->lvlv)
		}

		// sanity check, should never reach this place
		else
		{
			return -2; // invalid
		}
		
	}

	// default, should never reach this place
	return 0;
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

PtEtaPhiMVector GenerateLorentzVectorNeutrinoTrue(vector<PtEtaPhiMVector> lvecsHdecay, vector<int> pdgIdsHdecay){
	for(int i=0; i<pdgIdsHdecay.size(); i++){
		if(abs(pdgIdsHdecay[i])==12 || abs(pdgIdsHdecay[i])==14 || abs(pdgIdsHdecay[i])==16){
			return lvecsHdecay[i];
		}
	}
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


PtEtaPhiEVector GenerateWhadLvec(vector<PtEtaPhiEVector> lvecs, vector<int> indicies_fixxed)
{
	if(indicies_fixxed[6]==-1 || indicies_fixxed[7]==-1 || indicies_fixxed[6]==indicies_fixxed[7])
	{
		return PtEtaPhiEVector(0,0,0,0);
	}

	return lvecs[indicies_fixxed[6]] + lvecs[indicies_fixxed[7]];
}



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


int ClassifyTrueOnshellWHad(int classification_true_higgs_decay, float m_true_w1, float m_true_w2){
	// w1 decays hadronically
	if(classification_true_higgs_decay==11 || classification_true_higgs_decay==12){
		if( abs(m_true_w1-PDG_MASS_WBOSON)<THRESHOLD_ONSHELL_DEFINITION )
		{
			return 1;
		}
		else
		{
			return -1;
		}
	}

	// w2 decays hadronically
	if(classification_true_higgs_decay==13 || classification_true_higgs_decay==14){
		if( abs(m_true_w2-PDG_MASS_WBOSON)<THRESHOLD_ONSHELL_DEFINITION )
		{
			return 2;
		}
		else
		{
			return -2;
		}
	}

	return 0;
}


float ExtractTrueInformationLeptonFromHiggs(int classification_true_higgs_decay, float info_w11, float info_w12, float info_w21, float info_w22)
{
	
	// H->WW->qqlv signature
	if(classification_true_higgs_decay==11)
	{
		return info_w21;
	}
	// H->WW->qqvl signature
	if(classification_true_higgs_decay==12)
	{
		return info_w22;
	}
	// H->WW->lvqq signature
	if(classification_true_higgs_decay==13)
	{
		return info_w11;
	}
	// H->WW->vlqq signature	
	if(classification_true_higgs_decay==14)
	{
		return info_w12;
	}

	// no lepton from Higgs decays
	return -1;
}


PtEtaPhiMVector CombineTrueWHadFromHiggs(int classification_true_higgs_decay, PtEtaPhiMVector lvec_w11, PtEtaPhiMVector lvec_w12, PtEtaPhiMVector lvec_w21, PtEtaPhiMVector lvec_w22)
{
	// H->WW->qqlv or H->WW->qqvl signature
	if(classification_true_higgs_decay==11 || classification_true_higgs_decay==12)
	{
		return lvec_w11 + lvec_w12;
	}

	// H->WW->lvqq or H->WW->vlqq signature
	if(classification_true_higgs_decay==13 || classification_true_higgs_decay==14)
	{
		return lvec_w21 + lvec_w22;
	}

	// no fully hadronic decay from Higgs decays
	return PtEtaPhiMVector(-999,-999,-999,-999);
}

PtEtaPhiMVector CombineTrueWLepFromHiggs(int classification_true_higgs_decay, PtEtaPhiMVector lvec_w11, PtEtaPhiMVector lvec_w12, PtEtaPhiMVector lvec_w21, PtEtaPhiMVector lvec_w22)
{
	// H->WW->qqlv or H->WW->qqvl signature
	if(classification_true_higgs_decay==11 || classification_true_higgs_decay==12)
	{
		return lvec_w21 + lvec_w22;
	}

	// H->WW->lvqq or H->WW->vlqq signature
	if(classification_true_higgs_decay==13 || classification_true_higgs_decay==14)
	{
		return lvec_w11 + lvec_w12;
	}

	// no fully hadronic decay from Higgs decays
	return PtEtaPhiMVector(-999,-999,-999,-999);
}



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


int FillXsecMaps()
{
    // Open file
    std::ifstream file(XSEC_PATH);
    if (!file.is_open())
	{
        std::cerr << "Unable to open file: " << XSEC_PATH << std::endl;
        return 1;
    }


	// Fill File
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int dataset_number;
        double crossSection_pb, kFactor, genFiltEff;

        // Columns we are interested in: (1) dataset_number, (3) crossSection_pb, (5) kFactor
        std::string dataset_str, physics_short, genFiltEff_str, kFactor_str, relUncertUP, relUncertDOWN, generator_name, etag;

        if (!(iss >> dataset_number >> physics_short >> crossSection_pb >> genFiltEff >> kFactor >> relUncertUP >> relUncertDOWN >> generator_name >> etag)) {
            continue; // Skip the line if it can't be parsed
        }

        // Fill the dictionaries
        xSecs[dataset_number] = crossSection_pb;
		genFiltersEff[dataset_number] = genFiltEff;
        kFactors[dataset_number] = kFactor;
    }


    // Close the file
    file.close();

    return 0;
}

int FillSumOfWeights(string fileName, int DSID)
{
	string filePath = INPUT_PATH + fileName;

    // Open the ROOT file
    TFile *file = TFile::Open(filePath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return 1;
    }

	// Variables to store the histogram once found
    TH1 *hist = nullptr;
    std::string histName;

    // Iterate over the keys in the ROOT file to find the matching histogram
    TIter next(file->GetListOfKeys());
    TKey *key;
    
	while((key = (TKey*)next()))
	{
        // Check if the object is a histogram (TH1F)
        if(key->GetClassName() == std::string("TH1F"))
		{
            histName = key->GetName();

            // Check if the histogram name matches the desired pattern
            if(histName.find("CutBookkeeper_") == 0)
			{
                //std::cout << "Found matching histogram: " << histName << std::endl;

                // Get the histogram
                hist = (TH1*)file->Get(histName.c_str());
                break;
            }
        }
    }

    // If a matching histogram was found, retrieve and print the content of bin 2
    if(hist) 
	{
		// get correct content
        double binContent = hist->GetBinContent(2);

		// save sum of weights of same DSID
		if(sumOfWeights.find(DSID) != sumOfWeights.end())
		{
			// DSID already exists
        	sumOfWeights[DSID] += binContent;
    	}
		else
		{
        	// Key does not exist, create it with value 'x'
        	sumOfWeights[DSID] = binContent;
    	}
    } 
	else {
        std::cerr << "No matching histogram found in " << fileName << std::endl;
		exit(1);
    }


    // Clean up
    file->Close();
	return 0;
}

int CountMatches(vector<int> jetToObjectIndiciesFixed)
{
	int counter = 0;
	for(int i=0; i<jetToObjectIndiciesFixed.size(); i++)
	{
		if(jetToObjectIndiciesFixed[i]!=-1)
		{
			counter++;
		}
	}
	return counter;
}


double GenerateSMxSec()
{
	return xSecs[currentDSID] * genFiltersEff[currentDSID] * kFactors[currentDSID];
}
double GenerateSMWeights(float mc_weight)
{
	return mc_weight * LUMINOSITY * xSecs[currentDSID] * genFiltersEff[currentDSID] * kFactors[currentDSID] / sumOfWeights[currentDSID];
}

int eventCounter = 0;

int CountEvents(int event_channel){
	// take every signal event
	if(event_channel>=10)
	{
		return 1;
	}

	eventCounter++;
	// ignore event truntuation
	if(TRUNCATE<=0)
	{
		return 1;
	}
	// count events, truncate to every x element
	else if(eventCounter==TRUNCATE)
	{
		eventCounter = 0; // reset counter
		return 1;
	}
	
	return -1;
}

int EvalutateMatching(
	int pgdid_t1_W_q1, int pgdid_t1_W_q2, int pgdid_t2_W_q1, int pgdid_t2_W_q2, 
	int pgdid_H_d1, int pgdid_H_d2, 
	int pgdid_H_d1_d1, int pgdid_H_d1_d2, int pgdid_H_d2_d1, int pgdid_H_d2_d2,
	vector<int> final_match_mask, int event_channel
){
	if(event_channel<0)
	{
		return -2;
	}

	vector<int> successful_matches = {-1,-1,-1,-1,-1};

	// check t1
	if( !( abs(pgdid_t1_W_q1)>100 || abs(pgdid_t1_W_q2)>100 ) )
	{
		if(
			final_match_mask[TRUTH_PARTONS::b_from_t]!=-1 &&
			final_match_mask[TRUTH_PARTONS::Wdecay1_from_t]!=-1 &&
			final_match_mask[TRUTH_PARTONS::Wdecay2_from_t]!=-1
		)
		{
			successful_matches[3] = 1;
		}
		else
		{
			successful_matches[3] = 0;
		}
	}

	// check t2
	if( !( abs(pgdid_t2_W_q1)>100 || abs(pgdid_t2_W_q2)>100 ) )
	{
		if(
			final_match_mask[TRUTH_PARTONS::b_from_tbar]!=-1 &&
			final_match_mask[TRUTH_PARTONS::Wdecay1_from_tbar]!=-1 &&
			final_match_mask[TRUTH_PARTONS::Wdecay2_from_tbar]!=-1
		)
		{
			successful_matches[4] = 1;
		}
		else
		{
			successful_matches[4] = 0;
		}
	}


	// check ttbar
	if(( successful_matches[3]!=-1 || successful_matches[4]!=-1 ) )
	{
		if(successful_matches[3]==1 && successful_matches[4]==1)
		{
			successful_matches[2] = 1;
		}
		else
		{
			successful_matches[2] = 0;
		}
	}

	if( !( abs(pgdid_H_d1)>100 || abs(pgdid_t2_W_q2)>100 || abs(pgdid_H_d1_d1)>100 || abs(pgdid_H_d2_d1)>100 || abs(pgdid_H_d1_d2)>100 || abs(pgdid_H_d2_d2)>100 ) )
	{
		if(
			final_match_mask[TRUTH_PARTONS::Wdecay1_from_H]!=-1 &&
			final_match_mask[TRUTH_PARTONS::Wdecay2_from_H]!=-1
		)
		{
			successful_matches[1] = 1;
		}
		else
		{
			successful_matches[1] = 0;
		}
	}

	// check event
	if( successful_matches[1]!=-1 && successful_matches[2]!=-1 )
	{
		if(successful_matches[1]==1 && successful_matches[2]==1)
		{
			successful_matches[0] = 1;
		}
		else
		{
			successful_matches[0] = 0;
		}
	}

	// save entries for later, only which are possible
	if(successful_matches[0]!=-1) // ttH full event
	{
		successfulMatchHistory[0]+= successful_matches[0];
		successfulPossibleMatchHistory[0]+= 1;
	}
	if(successful_matches[1]!=-1) // H hadronic W
	{
		successfulMatchHistory[1]+= successful_matches[1];
		successfulPossibleMatchHistory[1]+= 1;
	}
	if(successful_matches[2]!=-1) // ttbar
	{
		successfulMatchHistory[2]+= successful_matches[2];
		successfulPossibleMatchHistory[2]+= 1;
	}
	if(successful_matches[3]!=-1) // t1
	{
		successfulMatchHistory[3]+= successful_matches[3];
		successfulPossibleMatchHistory[3]+= 1;
	}
	if(successful_matches[4]!=-1) // t2
	{
		successfulMatchHistory[4]+= successful_matches[4];
		successfulPossibleMatchHistory[4]+= 1;
	}

	// second signal only check
	if(event_channel==11)
	{	
		if(successful_matches[0]!=-1) // ttH full event
		{
			successfulMatchHistorySignal[0]+= successful_matches[0];
			successfulPossibleMatchHistorySignal[0]+= 1;
		}
		if(successful_matches[1]!=-1) // H hadronic W
		{
			successfulMatchHistorySignal[1]+= successful_matches[1];
			successfulPossibleMatchHistorySignal[1]+= 1;
		}
		if(successful_matches[2]!=-1) // ttbar
		{
			successfulMatchHistorySignal[2]+= successful_matches[2];
			successfulPossibleMatchHistorySignal[2]+= 1;
		}
		if(successful_matches[3]!=-1) // t1
		{
			successfulMatchHistorySignal[3]+= successful_matches[3];
			successfulPossibleMatchHistorySignal[3]+= 1;
		}
		if(successful_matches[4]!=-1) // t2
		{
			successfulMatchHistorySignal[4]+= successful_matches[4];
			successfulPossibleMatchHistorySignal[4]+= 1;
		}
	}


	return successful_matches[0];
}



// ==========  MAIN  ==========
// ===========================
int match(string input_file);


int main(int argc, char** argv)
{
	// prepare evaluation
	successfulMatchHistory[0] = 0;
	successfulMatchHistory[1] = 0;
	successfulMatchHistory[2] = 0;
	successfulMatchHistory[3] = 0;
	successfulMatchHistory[4] = 0;
	successfulPossibleMatchHistory[0] = 0;
	successfulPossibleMatchHistory[1] = 0;
	successfulPossibleMatchHistory[2] = 0;
	successfulPossibleMatchHistory[3] = 0;
	successfulPossibleMatchHistory[4] = 0;
	successfulMatchHistorySignal[0] = 0;
	successfulMatchHistorySignal[1] = 0;
	successfulMatchHistorySignal[2] = 0;
	successfulMatchHistorySignal[3] = 0;
	successfulMatchHistorySignal[4] = 0;
	successfulPossibleMatchHistorySignal[0] = 0;
	successfulPossibleMatchHistorySignal[1] = 0;
	successfulPossibleMatchHistorySignal[2] = 0;
	successfulPossibleMatchHistorySignal[3] = 0;
	successfulPossibleMatchHistorySignal[4] = 0;


	// fill dictionaries beforehand for properly weighted MC events
	std::cout << "Calculate proper SM weights" << std::endl;
	FillXsecMaps();


	// loop all samples for sum of weights 
	for(int i=0; i<sizeof(INPUT_FILE_NAMES)/sizeof(char*); i++)
	{
    	// gets DSID and estimates weights 
		currentDSID = std::stoi( extractDSID(INPUT_FILE_NAMES[i]) );	
		FillSumOfWeights(INPUT_FILE_NAMES[i], currentDSID);
	}

	for(const auto& entry : sumOfWeights)
	{
        std::cout << " > summed SM weights for " << entry.first << ": " << entry.second << std::endl;
    }


	std::cout << std::endl << "Loop trough all input files..." << std::endl << std::endl;
	// loop all samples for matching
	for(int i=0; i<sizeof(INPUT_FILE_NAMES)/sizeof(char*); i++)
	{
    	// get current DSID    
		currentDSID = std::stoi( extractDSID(INPUT_FILE_NAMES[i]) );		
		std::cout << INPUT_FILE_NAMES[i] << " -> "  << currentDSID  << std::endl;

		std::cout << "(" << i+1 << "/" << sizeof(INPUT_FILE_NAMES)/sizeof(char*) << ") - Processing: " << INPUT_FILE_NAMES[i] << std::endl;
		match(INPUT_FILE_NAMES[i]);
		std::cout << std::endl;

		if(i==0 && BREAK_AFTER_FIRST_FILE){
			break;
		}
	}


	// print matching evaluation
	std::cout << std::endl << "Overall Truth-Matching Evaluation" << std::endl;
	std::cout << " > all valid events" << std::endl;
	for(int i=0; i<5; i++)
	{
		if(successfulPossibleMatchHistory[i]>0)
		{
			std::cout << "   " << i << ": " << successfulMatchHistory[i] << "/" << successfulPossibleMatchHistory[i] << " = " << (float)successfulMatchHistory[i]/successfulPossibleMatchHistory[i] << std::endl;
		}
		else
		{
			std::cout << "   " <<  i << ": " << successfulMatchHistory[i] << "/" << successfulPossibleMatchHistory[i] << std::endl; 
		}
	}

	std::cout << std::endl << " > valid signal events only" << std::endl;
	for(int i=0; i<5; i++)
	{
		if(successfulPossibleMatchHistory[i]>0)
		{
			std::cout << "   " << i << ": " << successfulMatchHistorySignal[i] << "/" << successfulPossibleMatchHistorySignal[i] << " = " << (float)successfulMatchHistorySignal[i]/successfulPossibleMatchHistorySignal[i] << std::endl;
		}
		else
		{
			std::cout << "   " <<  i << ": " << successfulMatchHistorySignal[i] << "/" << successfulPossibleMatchHistorySignal[i] << std::endl; 
		}
	}

	return 0;
}


int match(string input_file)
{
	// SETUP
	// ==============================	
	// setup TChain
	cout << " > setup TChain" << endl;
	TChain rRecoChain("reco");
	TChain rTruthChain("truth");

	// link input file
	auto file_path = std::string();
	file_path.append(INPUT_PATH).append(input_file);
	rRecoChain.Add(file_path.c_str());
	rTruthChain.Add(file_path.c_str());

	// index TruthChain to kick out un-matched events - then declare friends
	rTruthChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rRecoChain.AddFriend(&rTruthChain);

	// setup RDataFrame
	cout << " > setup RDataFrame" << endl;
	auto rDataFrame = RDataFrame(rRecoChain);
	auto rLoopManager = rDataFrame.Range(0); // no limit on input

	// check if chains are matched properly
	auto nTotalEvents = rLoopManager.Count();
	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0){ 
		cout << " > Warning: there are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events! skipping..." << endl;
		return -1;
	}




	// GENERATE SM WEIGHTS
	// ==============================
	// get proper SM weights
	rLoopManager = rLoopManager.Define(
		"SM_event_xsecs", 
		GenerateSMxSec, 
		{}
	);
	rLoopManager = rLoopManager.Define(
		"SM_event_weight", 
		GenerateSMWeights, 
		{"weight_mc_NOSYS"}
	);


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
	rLoopManager = rLoopManager.Define(
		"jet_to_object_indicies_fixed",
		CollectJetToObjectIndiciesFixed,
		{"jet_final_match_mask"}
	);

	rLoopManager = rLoopManager.Define(
		"reco_whad_lvec",
		GenerateWhadLvec,
		{"lvecs_jets", "jet_to_object_indicies_fixed"}
	);

	rLoopManager = rLoopManager.Define(
		"reco_whad_pt",
		ExtractPtE,
		{"reco_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_whad_eta",
		ExtractEtaE,
		{"reco_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_whad_phi",
		ExtractPhiE,
		{"reco_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_whad_e",
		ExtractEE,
		{"reco_whad_lvec"}
	);


	// classification
	rLoopManager = rLoopManager.Define(
		"classification_true_t1_decay",
		ClassifyTrueTopDecay,
		{"Tth_MC_Wdecay1_from_t_pdgId", "Tth_MC_Wdecay2_from_t_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"classification_true_t2_decay",
		ClassifyTrueTopDecay,
		{"Tth_MC_Wdecay1_from_tbar_pdgId", "Tth_MC_Wdecay2_from_tbar_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"classification_true_higgs_decay",
		ClassifyTrueHiggsDecay,
		{"Tth_MC_Higgs_decay1_pdgId", "Tth_MC_Higgs_decay2_pdgId", "Tth_MC_Higgs_decay1_from_decay1_pdgId", "Tth_MC_Higgs_decay2_from_decay1_pdgId", "Tth_MC_Higgs_decay1_from_decay2_pdgId", "Tth_MC_Higgs_decay2_from_decay2_pdgId"}
	);
	rLoopManager = rLoopManager.Define(
		"classification_event_completion",
		ClassifyEventCompletion,
		{"classification_true_t1_decay", "classification_true_t2_decay", "classification_true_higgs_decay", "jet_to_object_indicies_fixed"}
	);
	rLoopManager = rLoopManager.Define(
		"classification_event_channel",
		ClassifyEventChannel,
		{
			"Tth_MC_Wdecay1_from_t_pdgId", "Tth_MC_Wdecay2_from_t_pdgId", "Tth_MC_Wdecay1_from_tbar_pdgId", "Tth_MC_Wdecay2_from_tbar_pdgId",
			"Tth_MC_Higgs_decay1_pdgId", "Tth_MC_Higgs_decay2_pdgId", 
			"Tth_MC_Higgs_decay1_from_decay1_pdgId", "Tth_MC_Higgs_decay2_from_decay1_pdgId", "Tth_MC_Higgs_decay1_from_decay2_pdgId", "Tth_MC_Higgs_decay2_from_decay2_pdgId"
		}
	);



	rLoopManager = rLoopManager.Define(
		"classification_onshell_whad",
		ClassifyTrueOnshellWHad,
		{"classification_true_higgs_decay", "Tth_MC_Higgs_decay1_m", "Tth_MC_Higgs_decay2_m"}
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
		{"classification_true_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_pt", "Tth_MC_Higgs_decay2_from_decay1_pt", "Tth_MC_Higgs_decay1_from_decay2_pt", "Tth_MC_Higgs_decay2_from_decay2_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_eta",
		ExtractTrueInformationLeptonFromHiggs,
		{"classification_true_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_eta", "Tth_MC_Higgs_decay2_from_decay1_eta", "Tth_MC_Higgs_decay1_from_decay2_eta", "Tth_MC_Higgs_decay2_from_decay2_eta"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_phi",
		ExtractTrueInformationLeptonFromHiggs,
		{"classification_true_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_phi", "Tth_MC_Higgs_decay2_from_decay1_phi", "Tth_MC_Higgs_decay1_from_decay2_phi", "Tth_MC_Higgs_decay2_from_decay2_phi"}
	);
	rLoopManager = rLoopManager.Define(
		"true_lepton_m",
		ExtractTrueInformationLeptonFromHiggs,
		{"classification_true_higgs_decay", "Tth_MC_Higgs_decay1_from_decay1_m", "Tth_MC_Higgs_decay2_from_decay1_m", "Tth_MC_Higgs_decay1_from_decay2_m", "Tth_MC_Higgs_decay2_from_decay2_m"}
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
		{"classification_true_higgs_decay", "true_w11_lvec", "true_w12_lvec", "true_w21_lvec", "true_w22_lvec"}
	);
	
	rLoopManager = rLoopManager.Define(
		"true_whad_pt",
		ExtractPt,
		{"true_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_whad_eta",
		ExtractEta,
		{"true_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_whad_phi",
		ExtractPhi,
		{"true_whad_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_whad_m",
		ExtractM,
		{"true_whad_lvec"}
	);
	
	rLoopManager = rLoopManager.Define(
		"true_wlep_lvec",
		CombineTrueWLepFromHiggs,
		{"classification_true_higgs_decay", "true_w11_lvec", "true_w12_lvec", "true_w21_lvec", "true_w22_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_wlep_pt",
		ExtractPt,
		{"true_wlep_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_wlep_eta",
		ExtractEta,
		{"true_wlep_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_wlep_phi",
		ExtractPhi,
		{"true_wlep_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_wlep_m",
		ExtractM,
		{"true_wlep_lvec"}
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
	rLoopManager = rLoopManager.Define(
		"true_neutrino_pt",
		ExtractPt,
		{"true_neutrino_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_neutrino_eta",
		ExtractEta,
		{"true_neutrino_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_neutrino_phi",
		ExtractPhi,
		{"true_neutrino_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_neutrino_m",
		ExtractM,
		{"true_neutrino_lvec"}
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
		"number_of_matches",
		CountMatches,
		{"jet_to_object_indicies_fixed"}
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

	// truth top kinematics
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
	



	// TRUNCTUATE
	// ==============================
	rLoopManager = rLoopManager.Define(
		"count_var", 
		CountEvents, 
		{"classification_event_channel"}
	);


	// EVALUATION
	// ==============================
	// count successful matches
	rLoopManager = rLoopManager.Define(
		"successful_matches",
		EvalutateMatching,
		{
			"Tth_MC_Wdecay1_from_t_pdgId", "Tth_MC_Wdecay2_from_t_pdgId", "Tth_MC_Wdecay1_from_tbar_pdgId", "Tth_MC_Wdecay2_from_tbar_pdgId",
			"Tth_MC_Higgs_decay1_pdgId", "Tth_MC_Higgs_decay2_pdgId", 
			"Tth_MC_Higgs_decay1_from_decay1_pdgId", "Tth_MC_Higgs_decay2_from_decay1_pdgId", "Tth_MC_Higgs_decay1_from_decay2_pdgId", "Tth_MC_Higgs_decay2_from_decay2_pdgId",
			"jet_final_match_mask", "classification_event_channel"
		}
	);


	// FINALISE
	// ==============================
	// apply filter	and limit if needed
	auto rLoopManagerFiltered = rLoopManager.Filter(FILTER).Range(MAX_NUMBER_OF_EVENTS);
	cout << " > " << rLoopManagerFiltered.Count().GetValue() << " events passed the following selection: " << FILTER << endl;


	// save snapshot to disk
	cout << " > saving snapshot" << endl;
	string output_file_name = OUTPUT_PATH + "matched_" + input_file.substr(14,6) +  "_" + input_file.substr(input_file.find("r13"),6) + "_" + input_file.substr(input_file.size()-18,6) + "_" + std::to_string( rLoopManagerFiltered.Count().GetValue() ) + "e.root";
	rLoopManagerFiltered.Snapshot(
		"matched", 
		output_file_name,
		OUTPUT_COLOUMN_NAMES
	);
		
	return 0;
}



// =========  FUNCTION DEFINITION  ===========
// ===========================================

// ==========  GENERATE LORENTZ VECTORS
PtEtaPhiMVector GenerateLorentzVectorMHiggsDecision(Float_t pt1, Float_t eta1, Float_t phi1, Float_t mass1, Int_t pdgId1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t mass2, Int_t pdgId2){
	if(abs(pdgId1)>=1 && abs(pdgId1)<=8){
		return PtEtaPhiMVector(pt1,eta1,phi1,mass1);
	}

	if(abs(pdgId2)>=1 && abs(pdgId2)<=8){
		return PtEtaPhiMVector(pt2,eta2,phi2,mass2);
	}
	
	
	return PtEtaPhiMVector(0,0,0,0);
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

		if(closestTruthMatchValue!=-1)
		{ // closest available truth object is matched
			jetFinalMatchMasks[i] = 1<<closestTruthMatchValue;
			unavailableTruthValues.push_back(closestTruthMatchValue);
		}	
		else{
			jetFinalMatchMasks[i] = -1; // no truth object has been matched successfully, unmatched jet
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
			// skip invalid entries
			if(jetFinalMatchMasks[j]<0)
			{
				continue;
			}

			if( (jetFinalMatchMasks[j] & 1<<i)!=0 ){ //final match is 1 for current TRUTH_PARTON
				jetToObjectIndicies[i] = j;
				break;
			}
		}
	}

	return jetToObjectIndicies;
}





// ==========  HIGGS

int GetFilteredPdgIDs(Int_t pdgId){
	int newPdgId = abs(pdgId);
	
	if(abs(newPdgId)>100 || newPdgId==0){
		return -1;
	}

	return abs(newPdgId);
}

