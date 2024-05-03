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
#include<Math/VectorUtil.h>
#include<Math/Vector4D.h> // for PtEtaPhiEVector and PtEtaPhiMVector
#include<TMath.h>
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;



///==========  CONSTANTS  ==========///
//>> assumptions
const double MASS_HIGGS = 125; //SM Higgs mass
const double SIGMA_X = 0.5; //resolutions missing momentum
const double SIGMA_Y = 0.5; //resolutions missing momentum

//>> sampling neutrino eta
const double SAMPLE_ETA_NU_MIN	= -3;
const double SAMPLE_ETA_NU_MAX	= 3;
const double SAMPLE_ETA_NU_STEP	= 0.1;

//>> sampling mass leptonic W-boson
const double SAMPLE_MASS_WLEP_MIN 	= 0;
const double SAMPLE_MASS_WLEP_MAX 	= 45;
const double SAMPLE_MASS_WLEP_STEP 	= 1; 


//>> input/output directories/files
const string INPUT_PATH = "/home/ireas/git_repos/master/input/v1/user.ravinab.346343.PhPy8EG.DAOD_PHYS.e7148_s3681_r13144_p5855.20231104-v0_output/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"user.ravinab.35392295._000001.output.root"
};

const string OUTPUT_PATH = "/home/ireas/git_repos/master/output/";
const string OUTPUT_FILE = "final.root";
const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
	"lvec_Whad",
};



///==========  DECLARATION  ==========///
// using the final state information (lepton, hadronic W), solve for leptonic W
int ReconstructLeptonicWBoson(PtEtaPhiEVector lepton, PtEtaPhiEVector parton_w_had, float missing_energy_value, float missing_energy_phi);

// for a given assumption (leptonic W mass, neutrino eta, Higgs Mass) calculate neutrino solutions
vector<PtEtaPhiEVector> EstimateNeutrinoSolutions(PtEtaPhiEVector lvec_lepton, PtEtaPhiEVector lvec_Whad, double m_Wlep, double eta_nu, double m_H);

// calculate neutrino weight
double CalculateNeutrinoWeight(PtEtaPhiEVector lvec_neutrino, double missing_energy_x, double missing_energy_y);

// get hadronic lorentz vector from prediced indicies
PtEtaPhiEVector GenerateLorentzVectorWHad(vector<PtEtaPhiEVector> lvecs_jets, int index_q1, int index_q2);



///==========  MAIN BODY  ==========///
int main(){
	//>> setup TChain
	TChain rRecoChain("reco");
	TChain rTruthChain("truth");
	TChain rPredictionChain("prediction");

	for(auto input_file_name:INPUT_FILE_NAMES){
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);
		rRecoChain.Add(file_path.c_str());
		rTruthChain.Add(file_path.c_str());
		rPredictionChain.Add(file_path.c_str());
	}


	//>> build befriended rDataFrame
	rPredictionChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rTruthChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rRecoChain.AddFriend(&rPredictionChain);
	rRecoChain.AddFriend(&rTruthChain);

	auto rDataFrame = RDataFrame(rRecoChain);
	auto rLoopManager = rDataFrame.Range(0);

	//auto nTotalEvents = rLoopManager.Count();
	//auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != truth.mcChannelNumber || eventNumber != truth.eventNumber").Count();
	//if(nMismatchedEvents.GetValue()>0){ 
	//	PLOG_ERROR << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!";
	//	exit(1);
	//}


	//>> prepare inputs
	rLoopManager = rLoopManager.Define(
		"lvec_Whad", 
		GenerateLorentzVectorWHad, 
		{"jet_lvecs", "prediction/HW_q1", "prediction/HW_q2"}
	);


	//>> reconstruct event
	rLoopManager = rLoopManager.Define(
		"test_reco_int", 
		ReconstructLeptonicWBoson, 
		{"lvec_lepton", "lvec_Whad", "missing_energy_value", "missing_energy_phi"}
	);


	//>> save snapshot
	rLoopManagerFiltered.Snapshot(
		"final", 
		OUTPUT_PATH+OUTPUT_FILE,
		OUTPUT_COLOUMN_NAMES
	);

	return 0;
}



///==========  DEFINITION  ==========///
int ReconstructLeptonicWBoson(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_Whad,
	float missing_energy_value,
	float missing_energy_phi
)
{	
	///----------  calculate intermediate variables  ----------///
	double missing_energy_x = missing_energy_value * cos( missing_energy_phi ) / 1e3;
	double missing_energy_y = missing_energy_value * sin( missing_energy_phi ) / 1e3;
	
	
	///----------  grid sampling for unknown parameters  ----------///
	//>> sampling done in two loops
	//>> start at constant MIN value
	//>> add STEP value each iteration
	//>> stop when exceeding MAX value
	for(double sample_mass_wlep=SAMPLE_MASS_WLEP_MIN; sample_mass_wlep<=SAMPLE_MASS_WLEP_MAX; sample_mass_wlep+=SAMPLE_MASS_WLEP_STEP){
		for(double sample_eta_nu=SAMPLE_ETA_NU_MIN; sample_eta_nu<=SAMPLE_ETA_NU_MAX; sample_eta_nu+=SAMPLE_ETA_NU_STEP){
			//>> sampled variables:
			//>> sample_mass_wlep = mass of the leptonic W boson from H (assumed off-shell)
			//>> sample_eta_nu = eta of the unknown neutrino from leptonic W boson

			vector<PtEtaPhiEVector> neutrino_solutions = EstimateNeutrinoSolutions(lvec_lepton, lvec_Whad, sample_mass_wlep, sample_eta_nu, MASS_HIGGS);
			cout << "(" << sample_mass_wlep << "," << sample_eta_nu << ") yields " << neutrino_solutions.size() << " solutions:" << endl;
			
			for(PtEtaPhiEVector neutrino_solution : neutrino_solutions){
				double neutrino_weight = CalculateNeutrinoWeight(neutrino_solution, missing_energy_x, missing_energy_y);
				cout << "  > " << neutrino_weight << endl;
			}
		}
	}

	return 1;
}


vector<PtEtaPhiEVector> EstimateNeutrinoSolutions(
	PtEtaPhiEVector lvec_lepton, 
	PtEtaPhiEVector lvec_Whad,
    double m_Wlep,	//sampled 
	double eta_nu,	//sampled
	double m_H		//assumed
)
{
	///-----  save lvec properties for faster access  -----///	
	double E_lepton  = lvec_lepton.E();
	double px_lepton = lvec_lepton.Px();
	double py_lepton = lvec_lepton.Py();
	double pz_lepton = lvec_lepton.Pz();
	double E_Whad  = lvec_Whad.E();
	double px_Whad = lvec_Whad.Px();
	double py_Whad = lvec_Whad.Py();
	double pz_Whad = lvec_Whad.Pz();
	double m_Whad = lvec_Whad.M();
	double m_Wlep_squared = m_Wlep*m_Wlep;

	double cosh_eta_nu = cosh(eta_nu);
	double sinh_eta_nu = sinh(eta_nu);

	///-----  calculate intermediate variables  -----///	
	double E_lepton_prime  = E_lepton * cosh_eta_nu - pz_lepton * sinh_eta_nu; //TODO: what is this and following?
	double E_Whad_prime  = E_Whad * cosh_eta_nu - pz_Whad * sinh_eta_nu;

	double A 	= (py_lepton * E_Whad_prime - py_Whad * E_lepton_prime) / (px_Whad * E_lepton_prime - px_lepton * E_Whad_prime);
	double B 	= (E_lepton_prime * (m_H * m_H - m_Wlep_squared - m_Whad * m_Whad - 2. * lvec_lepton.Dot(lvec_Whad)) - E_Whad_prime * m_Wlep_squared) / (2. * (px_lepton * E_Whad_prime - px_Whad * E_lepton_prime));
	double par1 = (px_lepton * A + py_lepton) / E_lepton_prime;
	double C    = A * A + 1. - par1 * par1;
  	double par2 = (m_Wlep_squared / 2. + px_lepton * B) / E_lepton_prime;
  	double D    = 2. * (A * B - par2 * par1);
	double F    = B * B - par2 * par2;

	double det = D * D - 4. * C * F; // determinante yields number of solutions


	///-----  calculate solutions for neutrinos  -----///
  	vector<PtEtaPhiEVector> sol;
	
  	//>>> 1 solution case
  	if(det==0.){ //TODO: is this working properly? zero float comparison problem
    	double py1 = -D / (2. * C); //TODO: What happens here
    	double px1 = A * py1 + B;
    	double pT2_1 = px1 * px1 + py1 * py1;
    	double pz1 = sqrt(pT2_1) * sinh_eta_nu;

    	PtEtaPhiEVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    	if (!TMath::IsNaN(a1.E()) )
    		sol.push_back(a1);
  	}
  	//>>> 2 solution case
	else if(det>0){
    	double tmp   = sqrt(det) / (2. * C);
    	double py1   = -D / (2. * C) + tmp;
    	double py2   = -D / (2. * C) - tmp;
    	double px1   = A * py1 + B;
    	double px2   = A * py2 + B;
    	double pT2_1 = px1 * px1 + py1 * py1;
    	double pT2_2 = px2 * px2 + py2 * py2;
    	double pz1   = sqrt(pT2_1) * sinh_eta_nu;
    	double pz2   = sqrt(pT2_2) * sinh_eta_nu;

    	PtEtaPhiEVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    	PtEtaPhiEVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    	if (!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
    		sol.push_back(a1);
    		sol.push_back(a2);
    	}
  	}

  	return sol;
}


double CalculateNeutrinoWeight(PtEtaPhiEVector lvec_neutrino, double missing_energy_x, double missing_energy_y){
	double weight_x = exp( -pow( (lvec_neutrino.Px()-missing_energy_x), 2) / pow(SIGMA_X, 2));
	double weight_y = exp( -pow( (lvec_neutrino.Py()-missing_energy_y), 2) / pow(SIGMA_Y, 2));
	return weight_x * weight_y;           
}
