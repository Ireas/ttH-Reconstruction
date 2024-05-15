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
const double SIGMA_X = 0.5; //resolutions missing transverse energy in x
const double SIGMA_Y = 0.5; //resolutions missing transverse energy in y

//>> sampling neutrino eta
const double SAMPLE_ETA_NU_MIN	= -3;
const double SAMPLE_ETA_NU_MAX	= 3;
const double SAMPLE_ETA_NU_STEP	= 0.02;

//>> sampling mass leptonic W-boson	
const double SAMPLE_MASS_WLEP_MIN 	= 0;
const double SAMPLE_MASS_WLEP_MAX 	= 45;
const double SAMPLE_MASS_WLEP_STEP 	= 1; 


//>> input/output directories/files
const string INPUT_PATH = "/home/ireas/git_repos/master/samples/injected/";
const char* INPUT_FILE_NAMES[1] = { // put into array for easier access
	"_injection.root"
};

const string OUTPUT_PATH = "/home/ireas/git_repos/master/samples/output/";
const string OUTPUT_FILE = "_output.root";
const initializer_list<string> OUTPUT_COLOUMN_NAMES = {
//	"lvec_Whad_true",
//	"lvec_Whad_predicted",
	"lvec_neutrino_true",
//	"lvec_neutrino_estimated",
//	"lvec_lepton",
	"legacy_prediction",
	"neutrino_weight",
	"estimated_eta_nu",
	"estimated_mass_Wlep",
};


PtEtaPhiEVector TempConv(PtEtaPhiMVector t){return PtEtaPhiEVector(t.Pt(), t.Eta(), t.Phi(), t.E());}



///==========  DECLARATION  ==========///
// using the final state information (lepton, hadronic W), solve for leptonic W
vector<float> ReconstructLeptonicWBoson(PtEtaPhiEVector lvec_lepton, PtEtaPhiEVector lvec_Whad, float missing_energy_value, float missing_energy_phi);

// for a given assumption (leptonic W mass, neutrino eta, Higgs Mass) calculate neutrino solutions
vector<TLorentzVector> EstimateNeutrinoSolutions(PtEtaPhiEVector lvec_lepton, PtEtaPhiEVector lvec_Whad, float m_Wlep, float eta_nu, float m_H);

// calculate neutrino weight
double CalculateNeutrinoWeight(PtEtaPhiEVector lvec_neutrino, double missing_energy_x, double missing_energy_y);

// get hadronic lorentz vector from prediced indicies
PtEtaPhiEVector GenerateLorentzVectorWHad(vector<PtEtaPhiEVector> lvecs_jets, int index_q1, int index_q2){
	return lvecs_jets[index_q1] + lvecs_jets[index_q2];
}


    
vector<float> OLDWReconstruction(
	PtEtaPhiEVector lepton,
	PtEtaPhiEVector parton_w_had,
	PtEtaPhiEVector nu
);



///==========  DEFINITION  ==========///
PtEtaPhiEVector GenerateLorentzVectorE(float pt, float eta, float phi, float energy){return PtEtaPhiEVector(pt,eta,phi,energy);}
PtEtaPhiMVector GenerateLorentzVectorM(Float_t pt, Float_t eta, Float_t phi, Float_t mass){return PtEtaPhiMVector(pt,eta,phi,mass);}

float ExtractWeight(vector<float> neutrino_weight_output){return neutrino_weight_output[0];}
float ExtractEstimatedEtaNu(vector<float> neutrino_weight_output){return neutrino_weight_output[1];}
float ExtractEstimatedMassWLep(vector<float> neutrino_weight_output){return neutrino_weight_output[2];}



///==========  MAIN BODY  ==========///
int main(){
	///----------  preparation  ----------///
	//>> setup TChain
	TChain rMatchedChain("matched");
	TChain rPredictionChain("prediction");

	for(auto input_file_name:INPUT_FILE_NAMES){
		auto file_path = std::string();
		file_path.append(INPUT_PATH).append(input_file_name);
		rMatchedChain.Add(file_path.c_str());
		rPredictionChain.Add(file_path.c_str());
	}

	//>> build befriended rDataFrameQ
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


	///----------  calculation  ----------///
//
//	//>> prepare inputs
//	rLoopManager = rLoopManager.Define(
//		"lvec_Whad_predicted", 
//		GenerateLorentzVectorWHad, 
//		{"lvecs_jets", "prediction.HW_q1", "prediction.HW_q2"}
//	);
//	rLoopManager = rLoopManager.Define(
//		"lvec_lepton", 
//		GenerateLorentzVectorE, 
//		{"lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"}
//	);
//
//
//	//>> neutrino weighting algorithm
//	rLoopManager = rLoopManager.Define(
//		"neutrino_weighting_output_vector", 
//		ReconstructLeptonicWBoson, 
//		{"lvec_lepton", "lvec_Whad_predicted", "missing_energy_value", "missing_energy_phi"}
//	);
//
//


	rLoopManager = rLoopManager.Define(
		"lvec_neutrino", 
		TempConv, 
		{"lvec_neutrino_true"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_lepton", 
		TempConv, 
		{"lvec_lepton_true"}
	);
	rLoopManager = rLoopManager.Define(
		"lvec_Whad", 
		TempConv, 
		{"lvec_Whad_true"}
	);

	rLoopManager = rLoopManager.Define(
		"legacy_prediction", 
		OLDWReconstruction, 
		{"lvec_lepton", "lvec_Whad", "lvec_neutrino"}
	);
	//>> extract estimations
	rLoopManager = rLoopManager.Define(
		"neutrino_weight", 
		ExtractWeight, 
		{"legacy_prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"estimated_eta_nu", 
		ExtractEstimatedEtaNu, 
		{"legacy_prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"estimated_mass_Wlep", 
		ExtractEstimatedMassWLep, 
		{"legacy_prediction"}
	);


	///----------  output  ----------///
	//>> save snapshot
	rLoopManager.Snapshot(
		"neutrino_weighting", 
		OUTPUT_PATH+OUTPUT_FILE,
		OUTPUT_COLOUMN_NAMES
	);

	return 0;
}



///==========  DEFINITION  ==========///
vector<float> ReconstructLeptonicWBoson(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_Whad,
	float missing_energy_value,
	float missing_energy_phi
)
{	
	float best_weight = -999;
	float best_eta_nu = -999;
	float best_mass_Wlep = -999;
	
	///----------  calculate intermediate variables  ----------///
	float missing_energy_x = missing_energy_value * cos( missing_energy_phi ) / 1e3;
	float missing_energy_y = missing_energy_value * sin( missing_energy_phi ) / 1e3;
	
	
	///----------  grid sampling for unknown parameters  ----------///
	//>> sampling done in two loops
	//>> start at constant MIN value
	//>> add STEP value each iteration
	//>> stop when exceeding MAX value
	for(float sample_mass_wlep=SAMPLE_MASS_WLEP_MIN; sample_mass_wlep<=SAMPLE_MASS_WLEP_MAX; sample_mass_wlep+=SAMPLE_MASS_WLEP_STEP){
		for(float sample_eta_nu=SAMPLE_ETA_NU_MIN; sample_eta_nu<=SAMPLE_ETA_NU_MAX; sample_eta_nu+=SAMPLE_ETA_NU_STEP){
			//>> sampled variables:
			//>> sample_mass_wlep = mass of the leptonic W boson from H (assumed off-shell)
			//>> sample_eta_nu = eta of the unknown neutrino from leptonic W boson

			vector<TLorentzVector> neutrino_solutions_legacy = EstimateNeutrinoSolutions(lvec_lepton, lvec_Whad, sample_mass_wlep, sample_eta_nu, MASS_HIGGS);
			
			//>> convert legacy class TLorentzVector to new PtEtaPhiEVector
			vector<PtEtaPhiEVector> neutrino_solutions = {};
			for(TLorentzVector neutrino_solution : neutrino_solutions_legacy){
				neutrino_solutions.push_back( PtEtaPhiEVector(neutrino_solution.Pt(), neutrino_solution.Eta(), neutrino_solution.Phi(), neutrino_solution.E() ) );
			}

			for(PtEtaPhiEVector neutrino_solution : neutrino_solutions){
				float neutrino_weight = CalculateNeutrinoWeight(neutrino_solution, missing_energy_x, missing_energy_y);
				
				if(neutrino_weight>best_weight){
					best_weight = neutrino_weight;
					best_eta_nu = sample_eta_nu;
					best_mass_Wlep = sample_mass_wlep;
				}
			}
		}
	}

	return vector<float>{best_weight, best_eta_nu, best_mass_Wlep};
}


vector<TLorentzVector> EstimateNeutrinoSolutions(
	PtEtaPhiEVector lvec_lepton, 
	PtEtaPhiEVector lvec_Whad,
    float m_Wlep,	//sampled 
	float eta_nu,	//sampled
	float m_H		//assumed
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
  	vector<TLorentzVector> sol;
	
  	//>>> 1 solution case
  	if(det==0.){ //TODO: is this working properly? zero float comparison problem
    	double py1 = -D / (2. * C); //TODO: What happens here
    	double px1 = A * py1 + B;
    	double pT2_1 = px1 * px1 + py1 * py1;
    	double pz1 = sqrt(pT2_1) * sinh_eta_nu;

    	TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    	if(!TMath::IsNaN(a1.E())){
			sol.push_back(a1);

		}
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

    	TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    	TLorentzVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    	if(!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
    		sol.push_back(a1);
    		sol.push_back(a2);
    	}
  	}

  	return std::move(sol);
}


double CalculateNeutrinoWeight(PtEtaPhiEVector lvec_neutrino, double missing_energy_x, double missing_energy_y){
	double weight_x = exp( -pow( (lvec_neutrino.Px()-missing_energy_x), 2) / pow(SIGMA_X, 2));
	double weight_y = exp( -pow( (lvec_neutrino.Py()-missing_energy_y), 2) / pow(SIGMA_Y, 2));
	return weight_x * weight_y;           
}



double OLDCalculateNeutrinoWeight(
	PtEtaPhiEVector lvec_neutrino,
	double missing_energy_x, 
	double missing_energy_y
){
	double weight_x = exp( -pow( (lvec_neutrino.Px()-missing_energy_x), 2) / pow(SIGMA_X, 2));
	double weight_y = exp( -pow( (lvec_neutrino.Py()-missing_energy_y), 2) / pow(SIGMA_Y, 2));
	return weight_x * weight_y;           
}



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
  double C    = A * A + 1. - par1 * par1;
  double par2 = (Wmass2 / 2. + lepton->Px() * B) / Elprime;
  double D    = 2. * (A * B - par2 * par1);
  double F    = B * B - par2 * par2;
  double det  = D * D - 4. * C * F;

  std::vector<TLorentzVector> sol;

  ///-- 0 solutions case --///
  if (det < 0.0){
    return std::move(sol);
  }

  ///-- Only one real solution case --///
  if (det == 0.) {
    double py1 = -D / (2. * C);
    double px1 = A * py1 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pz1 = sqrt(pT2_1) * nu_sinh;

    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    if (!TMath::IsNaN(a1.E()) )
      sol.push_back(a1);
    return std::move(sol);
  }

  ///-- 2 solutions case --///
  if(det > 0){
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

    if (!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
      sol.push_back(a1);
      sol.push_back(a2);
    }
    return std::move(sol);
  }
  
  ///-- Should never reach this point --///
  return std::move(sol);
}


      
      
vector<float> OLDWReconstruction(
	PtEtaPhiEVector lepton,
	PtEtaPhiEVector parton_w_had,
	PtEtaPhiEVector nu
){
  float m_best_NW_weight = -1000;
  TLorentzVector* particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector* particle_reco_nu= new TLorentzVector(0.,0.,0.,0.);
  //Implementing only the linear mass range variation, for now
  std::vector<double> wmass_points;
  int wmass_linear_nbins = 50;
  double wmass_linear_min   = 0;
  double wmass_linear_max   = 50;
  double wmass_stepsize = (wmass_linear_max - wmass_linear_min)/wmass_linear_nbins;
  for (int i=0; i< wmass_linear_nbins; i++){
      wmass_points.push_back(wmass_linear_min + i*wmass_stepsize);
  }

  // Nu eta linear
  std::vector<double> nueta_linear_points;
  int nueta_linear_nbins = 50;
  double nueta_linear_min   = -3.0;
  double nueta_linear_max   = +3.0;
  double nueta_stepsize = (nueta_linear_max - nueta_linear_min)/nueta_linear_nbins;
  for (int i=0; i< nueta_linear_nbins; i++){
      nueta_linear_points.push_back(nueta_linear_min + i*nueta_stepsize + 0.001);
  }


  TLorentzVector* lepton_lv= new TLorentzVector();
  TLorentzVector* parton_w_had_lv= new TLorentzVector();
  lepton_lv->SetPtEtaPhiM(lepton.pt(), lepton.eta(), lepton.phi(), lepton.mass());
  parton_w_had_lv->SetPtEtaPhiM(parton_w_had.pt(), parton_w_had.eta(), parton_w_had.phi(), parton_w_had.mass());   

  double weight = -99.;
  double weight_for_hist = -99;
  double temp_weight_ex = 0;
  double temp_weight_ey = 0;
  double temp_weight =0;

  for(double w_mass : wmass_points){
      for (double nu_eta : nueta_linear_points){
          std::vector<TLorentzVector> neutrinos;
          neutrinos = solveForNeutrinoEta(lepton_lv , parton_w_had_lv, nu_eta, 125., w_mass);
          
          weight_for_hist = -99;
          temp_weight_ex = 0;
          temp_weight_ey = 0;
          temp_weight =0;

          for (auto neutrino : neutrinos){

              TLorentzVector temp_w_boson = *(lepton_lv) + neutrino;
              TLorentzVector* particle_nu = new TLorentzVector();
              particle_nu->SetPtEtaPhiE(nu.pt(), nu.eta(), nu.phi(), nu.E());
              temp_weight_ex = exp(-pow((neutrino.Px() - particle_nu->Px()),2) / 2*pow(0.5,2));
              temp_weight_ey = exp(-pow((neutrino.Py() - particle_nu->Py()),2) / 2*pow(0.5,2));

              temp_weight    = temp_weight_ex*temp_weight_ey;                 
                      
              if (temp_weight > weight){
                weight = temp_weight;
                particle_reco_w_lep->SetPtEtaPhiM(temp_w_boson.Pt(), temp_w_boson.Eta(), temp_w_boson.Phi(), temp_w_boson.M());
                particle_reco_nu->SetPtEtaPhiM(neutrino.Pt(), neutrino.Eta(), neutrino.Phi(), neutrino.M()); 
                m_best_NW_weight=weight;
              }  

              if (temp_weight > weight_for_hist){
                weight_for_hist = temp_weight;
              }  
              }

              if (weight_for_hist == -99){
                weight_for_hist = 0.0;
              }
              
              
          }

}
	float best_weight = m_best_NW_weight;
	float best_eta_nu = particle_reco_nu->Eta();
	float best_mass_Wlep = particle_reco_w_lep->M();
return vector<float>{best_weight, best_eta_nu, best_mass_Wlep};
}