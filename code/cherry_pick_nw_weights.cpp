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
const double MASS_HIGGS = 125e3; //SM Higgs mass in MeV
const double SIGMA = 10e3; //resolutions missing transverse energy in MeV

//>> sampling neutrino eta
const double SAMPLE_ETA_NU_MIN = -3;
const double SAMPLE_ETA_NU_MAX = 3;
const double SAMPLE_ETA_NU_STEP	= 0.05; //120bins

//>> sampling mass leptonic W-boson in MeV
const double SAMPLE_MASS_WLEP_MIN = 0;
const double SAMPLE_MASS_WLEP_MAX = 50e3;
const double SAMPLE_MASS_WLEP_STEP = 2e2; //250bins



//>> input/output directories/files
const string INPUT_PATH = "/home/ireas/git_repos/master/samples/injected/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"_injection.root"
};

const string OUTPUT_PATH = "/home/ireas/git_repos/master/samples/output/";
const string OUTPUT_FILE = "_cherry.root";
const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
	// logistics
	"mcChannelNumber",
	"eventNumber",
	// event global information
	"number_of_jets",
	"reco_met_value", 
	"reco_met_phi",
	// true event signatures
	"signature_higgs_decay",
	"signature_onshell_whad",
	"signature_abs_lepton_pdgid",
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
	// neutrino weighting output
	"NW_on_true_vec_weights",
	"NW_on_true_vec_nu_etas",
	"NW_on_true_vec_wlep_masses",
};




///==========  FUNCTIONS  ==========///    
PtEtaPhiEVector ConvertLorentzVectorMToE(PtEtaPhiMVector lvec){return PtEtaPhiEVector{lvec.Pt(), lvec.Eta(), lvec.Phi(), lvec.E()};}


int ClassifyWhadPossible(vector<PtEtaPhiEVector> lvec_jets, int index1, int index2)
{
	if(index1==-1 || index2==-1)
		return -1;

	if(index1>lvec_jets.size() || index2>lvec_jets.size())
		return -2;

	if(index1==index2)
		return -3;
	
	return 0;
}

PtEtaPhiEVector CombineJetsFromIndicies(vector<PtEtaPhiEVector> lvec_jets, int index1, int index2)
{
	// check for invalid indicies
	if(index1==-1 || index2==-1 || index1>lvec_jets.size() || index2>lvec_jets.size() || index1==index2)
	{
		return PtEtaPhiEVector{-999,-999,-999,-999};
	}

	return lvec_jets[index1] + lvec_jets[index2];
}


// extract truth information from neutrino
float ExtractEnergy(PtEtaPhiMVector lvec){return lvec.E();}
float ExtractPhi(PtEtaPhiMVector lvec){return lvec.Phi();}


// temporary extraction from output float of NW algorithm
vector<float> ExtractWeights(vector<vector<float>> neutrino_weight_output){return neutrino_weight_output[0];}
vector<float> ExtractWlepMasses(vector<vector<float>> neutrino_weight_output){return neutrino_weight_output[1];}
vector<float> ExtractNuEtas(vector<vector<float>> neutrino_weight_output){return neutrino_weight_output[2];}



///==========  NEUTRINO WEIGHTING  ==========///    
std::vector<TLorentzVector> solveForNeutrinoEta(
	TLorentzVector* lvec_lepton, 
    TLorentzVector* lvec_whad, 
    double nu_eta, 
    double higgs_mass, 
    double wlep_mass
){
	double nu_cosh = cosh(nu_eta);
	double nu_sinh = sinh(nu_eta);
	double wlep_mass_squared   = wlep_mass*wlep_mass;
	double Whadmass = lvec_whad->M();
	double Elprime  = lvec_lepton->E() * nu_cosh - lvec_lepton->Pz() * nu_sinh;
	double Ebprime  = lvec_whad->E()   * nu_cosh - lvec_whad->Pz()   * nu_sinh;

	double A = (lvec_lepton->Py() * Ebprime - lvec_whad->Py() * Elprime) / (lvec_whad->Px() * Elprime - lvec_lepton->Px() * Ebprime);
	double B = (Elprime * (higgs_mass * higgs_mass - wlep_mass_squared - Whadmass * Whadmass - 2. * lvec_lepton->Dot(*lvec_whad)) - Ebprime * wlep_mass_squared) / (2. * (lvec_lepton->Px() * Ebprime - lvec_whad->Px() * Elprime));
	
	double par1 = (lvec_lepton->Px() * A + lvec_lepton->Py()) / Elprime;
	double C = A * A + 1. - par1 * par1;
	double par2 = (wlep_mass_squared / 2. + lvec_lepton->Px() * B) / Elprime;
	double D = 2. * (A * B - par2 * par1);
	double F = B * B - par2 * par2;
	double det = D * D - 4. * C * F;

  	
	std::vector<TLorentzVector> sol;

	///-- 0 solutions case --///
	if(det<0){
		return std::move(sol);
	}

	///-- Only one real solution case --///
  	if(det==0.){
    	double py1 = -D / (2. * C);
	    double px1 = A * py1 + B;
	    double pT2_1 = px1 * px1 + py1 * py1;
	    double pz1 = sqrt(pT2_1) * nu_sinh;
	
	    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    	if(!TMath::IsNaN(a1.E())){
			sol.push_back(a1);
		}

		return std::move(sol);
  }

	///-- 2 solutions case --///
	if(det>0){
		double tmp   = sqrt(det) / (2. * C);
		double py1   = -D / (2. * C) + tmp;
		double py2   = -D / (2. * C) - tmp;
    	double px1   = A * py1 + B;
    	double px2   = A * py2 + B;
    	double pT2_1 = px1 * px1 + py1 * py1;
    	double pT2_2 = px2 * px2 + py2 * py2;
    	double pz1   = sqrt(pT2_1) * nu_sinh;
    	double pz2   = sqrt(pT2_2) * nu_sinh;
    
		TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    	TLorentzVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    	if(!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
			sol.push_back(a1);
			sol.push_back(a2);
    	}

    	return std::move(sol);
	}
  
  
  ///-- Should never reach this point --///
  return std::move(sol);
}     


vector<vector<float>> NeutrinoWeighting(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi,
	float higgs_mass_smear
){ 

	// define weights
	float best_weight = -1;
	float sum_weight = -1;


	float missing_energy_x = met_value * cos(met_phi);
	float missing_energy_y = met_value * sin(met_phi);
	
	TLorentzVector* particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  	TLorentzVector* particle_reco_nu = new TLorentzVector(0.,0.,0.,0.);





	// convert PtEtaPhiEVector to legacy TLorentzVector
	TLorentzVector* legacy_lvec_lepton = new TLorentzVector();
  	TLorentzVector* legacy_lvec_whad = new TLorentzVector();
  	legacy_lvec_lepton->SetPtEtaPhiM(lvec_lepton.pt(), lvec_lepton.eta(), lvec_lepton.phi(), lvec_lepton.mass());
  	legacy_lvec_whad->SetPtEtaPhiM(lvec_whad.pt(), lvec_whad.eta(), lvec_whad.phi(), lvec_whad.mass());   


	// vectors for plotting cherry-picked	
	vector<float> vec_weights;	
	vector<float> vec_wlep_mass;	
	vector<float> vec_nu_eta;	


	// loop trough sampled points (wlep_m, nu_eta)
	for(double wlep_sampled_mass=SAMPLE_MASS_WLEP_MIN; wlep_sampled_mass<=SAMPLE_MASS_WLEP_MAX; wlep_sampled_mass+=SAMPLE_MASS_WLEP_STEP)
	{
    	for(double nu_sampled_eta=SAMPLE_ETA_NU_MIN; nu_sampled_eta<=SAMPLE_ETA_NU_MAX; nu_sampled_eta+=SAMPLE_ETA_NU_STEP)
		{
			// solve for neutrinos
			std::vector<TLorentzVector> neutrinos;
          	neutrinos = solveForNeutrinoEta(legacy_lvec_lepton , legacy_lvec_whad, nu_sampled_eta, MASS_HIGGS+higgs_mass_smear, wlep_sampled_mass);
          
          	double temp_weight_ex = 0;
          	double temp_weight_ey = 0;
          	double temp_weight = 0;
			double best_weight = -1;

			// check all possible neutrino solutions for best 
          	for(auto neutrino : neutrinos){
	
            	TLorentzVector temp_w_boson = *(legacy_lvec_lepton) + neutrino;

              	//TLorentzVector* particle_nu = new TLorentzVector();
              	temp_weight_ex = exp( -1 * pow( (neutrino.Px() - missing_energy_x) ,2) / pow(SIGMA,2) );
              	temp_weight_ey = exp( -1 * pow( (neutrino.Py() - missing_energy_y) ,2) / pow(SIGMA,2) );
				temp_weight = temp_weight_ex*temp_weight_ey;

				if(temp_weight>best_weight)
				{
					best_weight = temp_weight;
				}
			}
			
			vec_weights.push_back(best_weight);	
			vec_wlep_mass.push_back(wlep_sampled_mass);	
			vec_nu_eta.push_back(nu_sampled_eta);
		}
	}
	
	return vector<vector<float>>{vec_weights, vec_wlep_mass, vec_nu_eta};
}


vector<vector<float>> NeutrinoWeightingWrapper(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi
){ 
	// NW without any Higgs Smearing
	vector<vector<float>> neutrino_weighting_output = NeutrinoWeighting(lvec_lepton, lvec_whad, met_value, met_phi, 0);
	return neutrino_weighting_output;
}




///==========  MAIN BODY  ==========///
int main(){
	///----------  Preparation  ----------///
	//>> setup TChain
	TChain rMatchedChain("matched");
	TChain rPredictionChain("prediction");

	for(auto input_file_name:INPUT_FILE_NAMES){
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);
		rMatchedChain.Add(file_path.c_str());
		rPredictionChain.Add(file_path.c_str());
	}

	//>> build befriended rDataFrames
	rPredictionChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rMatchedChain.AddFriend(&rPredictionChain);

	auto rDataFrame = RDataFrame(rMatchedChain);
	auto rLoopManager = rDataFrame.Range(0);

	auto nTotalEvents = rLoopManager.Count();
	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != prediction.mcChannelNumber || eventNumber != prediction.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0){ 
		cout << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!" << endl;
		return -1;
	}



	///----------  Calculation  ----------///

	//>> prepare reco information
	rLoopManager = rLoopManager.Define(
		"true_met_value", 
		ExtractEnergy, 
		{"true_neutrino_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_met_phi", 
		ExtractPhi, 
		{"true_neutrino_lvec"}
	);

	rLoopManager = rLoopManager.Define(
		"true_lepton_lvec_converted", 
		ConvertLorentzVectorMToE, 
		{"true_lepton_lvec"}
	);
	rLoopManager = rLoopManager.Define(
		"true_whad_lvec_converted", 
		ConvertLorentzVectorMToE, 
		{"true_whad_lvec"}
	);
	

	//>> estimate using truth information 
	rLoopManager = rLoopManager.Define(
		"prediction_on_true", 
		NeutrinoWeightingWrapper, 
		{"true_lepton_lvec_converted", "true_whad_lvec_converted", "true_met_value", "true_met_phi"}
	);

	rLoopManager = rLoopManager.Define(
		"NW_on_true_vec_weights", 
		ExtractWeights, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_vec_wlep_masses", 
		ExtractWlepMasses, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_vec_nu_etas", 
		ExtractNuEtas, 
		{"prediction_on_true"}
	);
	
	///----------  Output  ----------///
	//>> save snapshot
	rLoopManager.Snapshot(
		"neutrino_weighting", 
		OUTPUT_PATH+OUTPUT_FILE,
		OUTPUT_COLOUMN_NAMES
	);

	return 0;
}
