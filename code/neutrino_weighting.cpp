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
const double SAMPLE_ETA_NU_STEP	= 0.1;

//>> sampling mass leptonic W-boson	in MeV
const double SAMPLE_MASS_WLEP_MIN = 0;
const double SAMPLE_MASS_WLEP_MAX = 50e3;
const double SAMPLE_MASS_WLEP_STEP = 2e2; 


//>> input/output directories/files
const string INPUT_PATH = "/home/ireas/git_repos/master/samples/injected/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"_injection.root"
};

const string OUTPUT_PATH = "/home/ireas/git_repos/master/samples/output/";
const string OUTPUT_FILE = "_output.root";
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
	// reco event classifier
	"classifier_whad_possible",
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
	"reco_whad_lvec",
	// wlep particle information
	"true_wlep_lvec",
	// neutrino weighting output
	"NW_on_true_weight",
	"NW_on_true_nu_eta",
	"NW_on_true_nu_phi",
	"NW_on_true_nu_px",
	"NW_on_true_nu_py",
	"NW_on_true_wlep_mass",
	"NW_on_reco_weight",
	"NW_on_reco_nu_eta",
	"NW_on_reco_nu_phi",
	"NW_on_reco_nu_px",
	"NW_on_reco_nu_py",
	"NW_on_reco_wlep_mass",
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
float ExtractWeight(vector<float> neutrino_weight_output){return neutrino_weight_output[0];}
float ExtractEstimatedEtaNu(vector<float> neutrino_weight_output){return neutrino_weight_output[1];}
float ExtractEstimatedMassWLep(vector<float> neutrino_weight_output){return neutrino_weight_output[2];}
float ExtractEstimatedPx(vector<float> neutrino_weight_output){return neutrino_weight_output[3];}
float ExtractEstimatedPy(vector<float> neutrino_weight_output){return neutrino_weight_output[4];}
float ExtractEstimatedPhiNu(vector<float> neutrino_weight_output){return neutrino_weight_output[5];}



///==========  NEUTRINO WEIGHTING  ==========///    
std::vector<TLorentzVector> solveForNeutrinoEta(
	TLorentzVector* lepton, 
    TLorentzVector* Whad, 
    double nu_eta, 
    double mHiggs, 
    double mWlep
){
	double nu_cosh = cosh(nu_eta);
	double nu_sinh = sinh(nu_eta);
	double Wmass2   = mWlep*mWlep;
	double Whadmass = Whad->M();
	double Elprime  = lepton->E() * nu_cosh - lepton->Pz() * nu_sinh;
	double Ebprime  = Whad->E()   * nu_cosh - Whad->Pz()   * nu_sinh;

	double A = (lepton->Py() * Ebprime - Whad->Py() * Elprime) / (Whad->Px() * Elprime - lepton->Px() * Ebprime);
	double B = (Elprime * (mHiggs * mHiggs - Wmass2 - Whadmass * Whadmass - 2. * lepton->Dot(*Whad)) - Ebprime * Wmass2) / (2. * (lepton->Px() * Ebprime - Whad->Px() * Elprime));
	
	double par1 = (lepton->Px() * A + lepton->Py()) / Elprime;
	double C = A * A + 1. - par1 * par1;
	double par2 = (Wmass2 / 2. + lepton->Px() * B) / Elprime;
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


vector<float> NeutrinoWeightingAlgorithm(
	PtEtaPhiEVector lepton,
	PtEtaPhiEVector parton_w_had,
	float met_value,
	float met_phi
){ 

	// define weights
	float best_weight = -1;
	float sum_weight = -1;


	float missing_energy_x = met_value * cos(met_phi);
	float missing_energy_y = met_value * sin(met_phi);
	
	float val_best_px = 0.0;
	float val_best_py = 0.0;  
	
	TLorentzVector* particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  	TLorentzVector* particle_reco_nu = new TLorentzVector(0.,0.,0.,0.);


	// convert PtEtaPhiEVector to legacy TLorentzVector
	TLorentzVector* lepton_lv = new TLorentzVector();
  	TLorentzVector* parton_w_had_lv = new TLorentzVector();
  	lepton_lv->SetPtEtaPhiM(lepton.pt(), lepton.eta(), lepton.phi(), lepton.mass());
  	parton_w_had_lv->SetPtEtaPhiM(parton_w_had.pt(), parton_w_had.eta(), parton_w_had.phi(), parton_w_had.mass());   

	// loop trough sampled points (wlep_m, nu_eta)
	for(double w_mass=SAMPLE_MASS_WLEP_MIN; w_mass<=SAMPLE_MASS_WLEP_MAX; w_mass+=SAMPLE_MASS_WLEP_STEP)
	{
    	for(double nu_eta=SAMPLE_ETA_NU_MIN; nu_eta<=SAMPLE_ETA_NU_MAX; nu_eta+=SAMPLE_ETA_NU_STEP)
		{
			// solve for neutrinos
			std::vector<TLorentzVector> neutrinos;
          	neutrinos = solveForNeutrinoEta(lepton_lv , parton_w_had_lv, nu_eta, MASS_HIGGS, w_mass);
          
          	double temp_weight_ex = 0;
          	double temp_weight_ey = 0;
          	double temp_weight = 0;

			// check all possible neutrino solutions for best 
          	for(auto neutrino : neutrinos){
	
            	TLorentzVector temp_w_boson = *(lepton_lv) + neutrino;

              	//TLorentzVector* particle_nu = new TLorentzVector();
              	temp_weight_ex = exp( -1 * pow( (neutrino.Px() - missing_energy_x) ,2) / pow(SIGMA,2) );
              	temp_weight_ey = exp( -1 * pow( (neutrino.Py() - missing_energy_y) ,2) / pow(SIGMA,2) );
				temp_weight = temp_weight_ex*temp_weight_ey;     
				
				if(temp_weight>best_weight){
					best_weight = temp_weight;
					val_best_px = neutrino.Px();
					val_best_py = neutrino.Py();
                	particle_reco_w_lep->SetPtEtaPhiM(temp_w_boson.Pt(), temp_w_boson.Eta(), temp_w_boson.Phi(), temp_w_boson.M());
                	particle_reco_nu->SetPtEtaPhiM(neutrino.Pt(), neutrino.Eta(), neutrino.Phi(), neutrino.M()); 
              	}
			}
		}
	}
	
	// count the events for overview
	float best_phi_nu = particle_reco_nu->Phi();
	float best_eta_nu = particle_reco_nu->Eta();
	float best_mass_Wlep = particle_reco_w_lep->M();

	return vector<float>{best_weight, best_eta_nu, best_mass_Wlep, val_best_px, val_best_py, best_phi_nu};
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
		NeutrinoWeightingAlgorithm, 
		{"true_lepton_lvec_converted", "true_whad_lvec_converted", "true_met_value", "true_met_phi"}
	);

	rLoopManager = rLoopManager.Define(
		"NW_on_true_weight", 
		ExtractWeight, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_nu_eta", 
		ExtractEstimatedEtaNu, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_wlep_mass", 
		ExtractEstimatedMassWLep, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_nu_px", 
		ExtractEstimatedPx, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_nu_py", 
		ExtractEstimatedPy, 
		{"prediction_on_true"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_true_nu_phi", 
		ExtractEstimatedPhiNu, 
		{"prediction_on_true"}
	);


	//>> prepare reco information
	rLoopManager = rLoopManager.Define(
		"classifier_whad_possible", 
		ClassifyWhadPossible, 
		{"lvecs_jets", "prediction.HW_q1", "prediction.HW_q2"}
	);

	rLoopManager = rLoopManager.Define(
		"reco_whad_lvec", 
		CombineJetsFromIndicies, 
		{"lvecs_jets", "prediction.HW_q1", "prediction.HW_q2"}
	);


	//>> estimate using reco information 
	rLoopManager = rLoopManager.Define(
		"prediction_on_reco", 
		NeutrinoWeightingAlgorithm, 
		{"reco_lepton_lvec", "reco_whad_lvec", "reco_met_value", "reco_met_phi"}
	);

	rLoopManager = rLoopManager.Define(
		"NW_on_reco_weight", 
		ExtractWeight, 
		{"prediction_on_reco"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_reco_nu_eta", 
		ExtractEstimatedEtaNu, 
		{"prediction_on_reco"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_reco_wlep_mass", 
		ExtractEstimatedMassWLep, 
		{"prediction_on_reco"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_reco_nu_px", 
		ExtractEstimatedPx, 
		{"prediction_on_reco"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_reco_nu_py", 
		ExtractEstimatedPy, 
		{"prediction_on_reco"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_on_reco_nu_phi", 
		ExtractEstimatedPhiNu, 
		{"prediction_on_reco"}
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