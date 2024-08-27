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
const double MASS_HIGGS = 125e3; //SM Higgs mass in MeV -> switch to pdg mass on data
const double SIGMA = 10e3; //resolutions missing transverse energy in MeV

//>> sampling neutrino eta
const double SAMPLE_ETA_NU_MIN = -3;
const double SAMPLE_ETA_NU_MAX = 3;
const double SAMPLE_ETA_NU_STEP	= 6e-2; //100 bins

//>> sampling mass leptonic W-boson in MeV
const double SAMPLE_MASS_WLEP_MIN = 0;
const double SAMPLE_MASS_WLEP_MAX = 50e3;
const double SAMPLE_MASS_WLEP_STEP = 50e1; //100 bins

//>> sampling mass higgs boson for smearing in MeV
const double SAMPLE_MASS_HIGGS_SMEAR_MIN = -1e3;
const double SAMPLE_MASS_HIGGS_SMEAR_MAX = 1e3;
const double SAMPLE_MASS_HIGGS_SMEAR_STEP = 5e2; // 5 steps 

//>> fine sampling neutrino eta in roi
const double SAMPLE_FINE_ETA_NU_MIN = -0.8;
const double SAMPLE_FINE_ETA_NU_MAX = +0.8;
const double SAMPLE_FINE_ETA_NU_STEP = 4e-2; //40 bins

//>> sampling mass higgs boson for smearing in MeV
const double SAMPLE_FINE_MASS_WLEP_MIN = -4000;
const double SAMPLE_FINE_MASS_WLEP_MAX = 4000;
const double SAMPLE_FINE_MASS_WLEP_STEP = 200; //40 bins


//>> input/output directories/files
const string INPUT_FILE = "/media/ireas/Data/v6/injected/all_8+j_1530557e_injected.root";
const string OUTPUT_FILE = "/media/ireas/Data/v6/weighted/all_8+j_1530557e_merged.root";

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
	"jet_btag_60wp", // btags
	"jet_btag_70wp",
	"jet_btag_77wp",
	"jet_btag_85wp",
	"jet_btag_continous",
	"jet_final_match_mask",
	"reco_lepton_pt", // lepton
	"reco_lepton_eta",
	"reco_lepton_phi",
	"reco_lepton_e",
	"reco_met_value", // met
	"reco_met_phi",

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

	// classificiation
	"classification_true_t1_decay",
	"classification_true_t2_decay",
	"classification_true_higgs_decay",
	"classification_event_completion",
	"classification_onshell_whad",

	// true event signatures
	"signature_abs_lepton_pdgid",

	// reco event classifier
	"higgs_decay_mode_custom",
	"higgs_decay_decay_mode",

	// NW results
	"NW_weight",
	"NW_weight_sum",
	"NW_H_smear",
	"NW_nu_eta",
	"NW_nu_phi",
	"NW_nu_px",
	"NW_nu_py",
	"NW_wlep_mass",

	// NW results full
	//"NW_full_weight",
	//"NW_full_nu_eta",
	//"NW_full_wlep_mass",

	// SPANet stuff
	"spanet.t1_q1",
	"spanet.t1_q2",
	"spanet.t1_b",
	"spanet.t2_q1",
	"spanet.t2_q2",
	"spanet.t2_b",
	"spanet.HW_q1",
	"spanet.HW_q2",
	"spanet.t1_assignment_probability",
	"spanet.t1_detection_probability",
	"spanet.t1_marginal_probability",
	"spanet.t2_assignment_probability",
	"spanet.t2_detection_probability",
	"spanet.t2_marginal_probability",
	"spanet.HW_assignment_probability",
	"spanet.HW_detection_probability",
	"spanet.HW_marginal_probability",
};




///==========  FUNCTIONS  ==========///    
// convert truth to reco lorentzvector
PtEtaPhiEVector ConvertLorentzVectorMToE(PtEtaPhiMVector lvec){return PtEtaPhiEVector(lvec.Pt(), lvec.Eta(), lvec.Phi(), lvec.E());}

PtEtaPhiMVector GenerateLorentzVectorM(Double_t pt, Float_t eta, Float_t phi, Float_t mass){return PtEtaPhiMVector(pt,eta,phi,mass);}
PtEtaPhiEVector GeneratePtEtaPhiEVector(Float_t pt, Float_t eta, Float_t phi, Float_t e){return PtEtaPhiEVector(pt,eta,phi,e);}


vector<PtEtaPhiEVector> GenerateJetLvecs(int number_of_jets, ROOT::RVec<Double_t> pts, ROOT::RVec<Double_t> etas, ROOT::RVec<Double_t> phis, ROOT::RVec<Double_t> energies)
{
	vector<PtEtaPhiEVector> jetLvecs;

	for(int i=0; i<number_of_jets; i++)
	{
		jetLvecs.push_back( PtEtaPhiEVector(pts[i],etas[i],phis[i],energies[i]) );
	}

	return jetLvecs;
}

// combine jets from given indicies
PtEtaPhiEVector CombineJetsFromIndicies(vector<PtEtaPhiEVector> lvec_jets, int index1, int index2)
{
	// check for invalid indicies
	if(index1==-1 || index2==-1 || index1>lvec_jets.size() || index2>lvec_jets.size() || index1==index2)
	{
		return PtEtaPhiEVector{-999,-999,-999,-999};
	}

	return lvec_jets[index1] + lvec_jets[index2];
}


// extract information from lvec
float ExtractPT(PtEtaPhiMVector lvec){return lvec.Pt();}
float ExtractPhi(PtEtaPhiMVector lvec){return lvec.Phi();}


// extraction from output float of NW algorithm
float ExtractWeight(vector<float> NWOutput){return NWOutput[0];}
float ExtractWeightSum(vector<float> NWOutput){return NWOutput[7];}
float ExtractHiggsSmear(vector<float> NWOutput){return NWOutput[1];}
float ExtractEstimatedEtaNu(vector<float> NWOutput){return NWOutput[2];}
float ExtractEstimatedMassWLep(vector<float> NWOutput){return NWOutput[3];}
float ExtractEstimatedPx(vector<float> NWOutput){return NWOutput[4];}
float ExtractEstimatedPy(vector<float> NWOutput){return NWOutput[5];}
float ExtractEstimatedPhiNu(vector<float> NWOutput){return NWOutput[6];}

vector<float> ExtractWeightFull(vector<vector<float>> NWOutputFull){return NWOutputFull[0];}
vector<float> ExtractWlepMassesFull(vector<vector<float>> NWOutputFull){return NWOutputFull[1];}
vector<float> ExtractNuEtaFull(vector<vector<float>> NWOutputFull){return NWOutputFull[2];}


///==========  NEUTRINO WEIGHTING  ==========///    
std::vector<TLorentzVector> solveForNeutrinos(
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


vector<float> NeutrinoWeighting(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi,
	float higgs_mass_smear
){ 

	// define weights
	float best_weight = -1;
	float sum_all_weights = -1;


	float missing_energy_x = met_value * cos(met_phi);
	float missing_energy_y = met_value * sin(met_phi);
	
	float val_best_px = 0.0;
	float val_best_py = 0.0;  
	
	TLorentzVector* particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  	TLorentzVector* particle_reco_nu = new TLorentzVector(0.,0.,0.,0.);


	// convert PtEtaPhiEVector to legacy TLorentzVector
	TLorentzVector* legacy_lvec_lepton = new TLorentzVector();
  	TLorentzVector* legacy_lvec_whad = new TLorentzVector();
  	legacy_lvec_lepton->SetPtEtaPhiM(lvec_lepton.pt(), lvec_lepton.eta(), lvec_lepton.phi(), lvec_lepton.mass());
  	legacy_lvec_whad->SetPtEtaPhiM(lvec_whad.pt(), lvec_whad.eta(), lvec_whad.phi(), lvec_whad.mass());   

	// loop trough sampled points (wlep_m, nu_eta)
	for(double wlep_sampled_mass=SAMPLE_MASS_WLEP_MIN; wlep_sampled_mass<=SAMPLE_MASS_WLEP_MAX; wlep_sampled_mass+=SAMPLE_MASS_WLEP_STEP)
	{
    	for(double nu_sampled_eta=SAMPLE_ETA_NU_MIN; nu_sampled_eta<=SAMPLE_ETA_NU_MAX; nu_sampled_eta+=SAMPLE_ETA_NU_STEP)
		{
			// solve for neutrinos
			std::vector<TLorentzVector> neutrinos;
          	neutrinos = solveForNeutrinos(legacy_lvec_lepton , legacy_lvec_whad, nu_sampled_eta, MASS_HIGGS+higgs_mass_smear, wlep_sampled_mass);
          
          	double temp_weight_ex = 0;
          	double temp_weight_ey = 0;
          	double temp_weight = 0;

			// check all possible neutrino solutions for best 
          	for(auto neutrino : neutrinos)
			{
            	TLorentzVector temp_w_boson = *(legacy_lvec_lepton) + neutrino;

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

			if(best_weight>-1)
			{
				sum_all_weights+= best_weight;
			}
		}
	}
	
	// count the events for overview
	float best_phi_nu = particle_reco_nu->Phi();
	float best_eta_nu = particle_reco_nu->Eta();
	float best_mass_Wlep = particle_reco_w_lep->M();

	return vector<float>{best_weight, higgs_mass_smear, best_eta_nu, best_mass_Wlep, val_best_px, val_best_py, best_phi_nu, sum_all_weights};
}


vector<float> SampleROI(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi,
	float higgs_mass_smear,
	vector<float> neutrino_weighting_output
){ 
	// check if no solution has been found before, then no finetuning needed
	if( ExtractWeight(neutrino_weighting_output)>=0 )
	{
		return neutrino_weighting_output;
	}


	// extract information
	float best_weight = ExtractWeight(neutrino_weighting_output);
	float estimated_eta_nu = ExtractEstimatedEtaNu(neutrino_weighting_output);
	float estimated_mass_wlep = ExtractEstimatedMassWLep(neutrino_weighting_output);


	// new containers for neutrino weighting outputs
	float best_eta_nu = ExtractEstimatedEtaNu(neutrino_weighting_output);
	float best_mass_Wlep = ExtractEstimatedMassWLep(neutrino_weighting_output);
	float val_best_px = ExtractEstimatedPx(neutrino_weighting_output);
	float val_best_py = ExtractEstimatedPy(neutrino_weighting_output);
	float best_phi_nu = ExtractEstimatedPhiNu(neutrino_weighting_output);

	// empty objects to fill in neutrino weighting
	TLorentzVector* particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  	TLorentzVector* particle_reco_nu = new TLorentzVector(0.,0.,0.,0.);
	float missing_energy_x = met_value * cos(met_phi);
	float missing_energy_y = met_value * sin(met_phi);

	// convert PtEtaPhiEVector to legacy TLorentzVector
	TLorentzVector* legacy_lvec_lepton = new TLorentzVector();
  	TLorentzVector* legacy_lvec_whad = new TLorentzVector();
  	legacy_lvec_lepton->SetPtEtaPhiM(lvec_lepton.pt(), lvec_lepton.eta(), lvec_lepton.phi(), lvec_lepton.mass());
  	legacy_lvec_whad->SetPtEtaPhiM(lvec_whad.pt(), lvec_whad.eta(), lvec_whad.phi(), lvec_whad.mass());   


	// recheck ROI with fine sampling steps
	for(double wlep_sampled_mass=estimated_mass_wlep+SAMPLE_FINE_MASS_WLEP_MIN; wlep_sampled_mass<=estimated_mass_wlep+SAMPLE_FINE_MASS_WLEP_MAX; wlep_sampled_mass+=SAMPLE_FINE_MASS_WLEP_STEP)
	{
    	for(double nu_sampled_eta=estimated_eta_nu+SAMPLE_FINE_ETA_NU_MIN; nu_sampled_eta<=estimated_eta_nu+SAMPLE_FINE_ETA_NU_MAX; nu_sampled_eta+=SAMPLE_FINE_ETA_NU_STEP)
		{
			// solve for neutrinos
			std::vector<TLorentzVector> neutrinos;
          	neutrinos = solveForNeutrinos(legacy_lvec_lepton , legacy_lvec_whad, nu_sampled_eta, MASS_HIGGS+higgs_mass_smear, wlep_sampled_mass);
          
          	double temp_weight_ex = 0;
          	double temp_weight_ey = 0;
          	double temp_weight = 0;

			// check all possible neutrino solutions for best 
          	for(auto neutrino : neutrinos)
			{
            	TLorentzVector temp_w_boson = *(legacy_lvec_lepton) + neutrino;

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


	// get best values for output
	best_phi_nu = particle_reco_nu->Phi();
	best_eta_nu = particle_reco_nu->Eta();
	best_mass_Wlep = particle_reco_w_lep->M();
	

	return vector<float>{best_weight, higgs_mass_smear, best_eta_nu, best_mass_Wlep, val_best_px, val_best_py, best_phi_nu};
}

int current_event_index = 0;

vector<float> NeutrinoWeightingWrapper(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi
){ 
	// verbose output for me
	cout << " > current_event_index = " << current_event_index << endl;
	current_event_index++;


	// NW without any Higgs Smearing
	vector<float> neutrino_weighting_output = NeutrinoWeighting(lvec_lepton, lvec_whad, met_value, met_phi, 0);

	return neutrino_weighting_output;
	// ignore everything related to optimises to reduce time


//	if( ExtractWeight(neutrino_weighting_output)>=0 )
//	{
//		// sample ROI for better results
//		vector<float> finetuned_neutrino_weighting_output = SampleROI(lvec_lepton, lvec_whad, met_value, met_phi, 0, neutrino_weighting_output);
//		return finetuned_neutrino_weighting_output;
//	}
//
//
//	// if no valid solution has been found, try with smearing higgs mass
//	for(float higgsSmear=SAMPLE_MASS_HIGGS_SMEAR_MIN; higgsSmear<=SAMPLE_MASS_HIGGS_SMEAR_MAX; higgsSmear+= SAMPLE_MASS_HIGGS_SMEAR_STEP)
//	{
//		neutrino_weighting_output = NeutrinoWeighting(lvec_lepton, lvec_whad, met_value, met_phi, higgsSmear);
//
//		if( ExtractWeight(neutrino_weighting_output)>=0 )
//		{
//			// sample ROI for better results
//			vector<float> finetuned_neutrino_weighting_output = SampleROI(lvec_lepton, lvec_whad, met_value, met_phi, 0, neutrino_weighting_output);
//			return finetuned_neutrino_weighting_output;
//		}
//	}
//
//	// no solution found, no need for sample roi
//	return neutrino_weighting_output;
}



vector<vector<float>> NeutrinoWeightingFull(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi,
	float higgs_mass_smear
){ 

	// define weights
	float missing_energy_x = met_value * cos(met_phi);
	float missing_energy_y = met_value * sin(met_phi);
	


	// convert PtEtaPhiEVector to legacy TLorentzVector
	TLorentzVector* particle_reco_w_lep = new TLorentzVector();
  	TLorentzVector* particle_reco_nu = new TLorentzVector();

	TLorentzVector* legacy_lvec_lepton = new TLorentzVector();
  	TLorentzVector* legacy_lvec_whad = new TLorentzVector();
  	legacy_lvec_lepton->SetPtEtaPhiM(lvec_lepton.pt(), lvec_lepton.eta(), lvec_lepton.phi(), lvec_lepton.mass());
  	legacy_lvec_whad->SetPtEtaPhiM(lvec_whad.pt(), lvec_whad.eta(), lvec_whad.phi(), lvec_whad.mass());   


	// vectors for plotting cherries	
	vector<float> vec_weights;	
	vector<float> vec_wlep_mass;	
	vector<float> vec_nu_eta;	


	// loop trough sampled points (wlep_m, nu_eta)
	for(double wlep_sampled_mass=SAMPLE_MASS_WLEP_MIN; wlep_sampled_mass<=SAMPLE_MASS_WLEP_MAX; wlep_sampled_mass+=SAMPLE_MASS_WLEP_STEP)
	{
    	for(double nu_sampled_eta=SAMPLE_ETA_NU_MIN; nu_sampled_eta<=SAMPLE_ETA_NU_MAX; nu_sampled_eta+=SAMPLE_ETA_NU_STEP)
		{
			// solve for neutrinos
			std::vector<TLorentzVector> neutrinos = solveForNeutrinos(legacy_lvec_lepton , legacy_lvec_whad, nu_sampled_eta, MASS_HIGGS+higgs_mass_smear, wlep_sampled_mass);
          
			double best_weight = -1;
          	double temp_weight = 0;
          	double temp_weight_ex = 0;
          	double temp_weight_ey = 0;

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


vector<vector<float>> NeutrinoWeightingWrapperFull(
	PtEtaPhiEVector lvec_lepton,
	PtEtaPhiEVector lvec_whad,
	float met_value,
	float met_phi
){ 
	// NW without any Higgs Smearing
	vector<vector<float>> neutrino_weighting_output = NeutrinoWeightingFull(lvec_lepton, lvec_whad, met_value, met_phi, 0);
	return neutrino_weighting_output;
}


// utilities
float RenameFloat(float target){return target;}


///==========  MAIN BODY  ==========///
int main(){
	///----------  Preparation  ----------///
	//>> setup TChain
	cout << " > setup TChain" << endl;
	TChain rMatchedChain("matched");
	TChain rPredictionChain("spanet");

	// append input file to chains
	rMatchedChain.Add(INPUT_FILE.c_str());
	rPredictionChain.Add(INPUT_FILE.c_str());


	//>> build befriended rDataFrames
	rPredictionChain.BuildIndex("mcChannelNumber", "eventNumber");  // just for security, use DSID too
	rMatchedChain.AddFriend(&rPredictionChain);

	auto rDataFrame = RDataFrame(rMatchedChain);
	auto rLoopManager = rDataFrame.Range(0);
	
	auto nTotalEvents = rLoopManager.Count();
	auto nMismatchedEvents = rLoopManager.Filter("mcChannelNumber != spanet.mcChannelNumber || eventNumber != spanet.eventNumber").Count();
	if(nMismatchedEvents.GetValue()>0)
	{ 
		cout << "There are " << nMismatchedEvents.GetValue() << " / " << nTotalEvents.GetValue() << " mismatched events!" << endl;
		return -1;
	}



	///----------  Calculation  ----------///

	cout << " > define branches" << endl;
	//>> truth conversions
	rLoopManager = rLoopManager.Define(
		"met_value_for_truth", 
		RenameFloat, 
		{"true_neutrino_pt"}
	);
	rLoopManager = rLoopManager.Define(
		"met_phi_for_truth", 
		RenameFloat, 
		{"true_neutrino_pt"}
	);

//	rLoopManager = rLoopManager.Define(
//		"lepton_lvec_for_truth", 
//		ConvertLorentzVectorMToE, 
//		{"true_lepton_lvec"}
//	);
//	rLoopManager = rLoopManager.Define(
//		"whad_lvec_for_truth", 
//		ConvertLorentzVectorMToE, 
//		{"true_whad_lvec"}
//	);

	rLoopManager = rLoopManager.Define(
		"lvecs_jets", 
		GenerateJetLvecs, 
		{"number_of_jets", "jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_whad_lvec", 
		CombineJetsFromIndicies, 
		{"lvecs_jets", "spanet.HW_q1", "spanet.HW_q2"}
	);
	rLoopManager = rLoopManager.Define(
		"reco_lepton_lvec", 
		GeneratePtEtaPhiEVector, 
		{"reco_lepton_pt", "reco_lepton_eta", "reco_lepton_phi", "reco_lepton_e"}
	);


	

	//>> NW information
	rLoopManager = rLoopManager.Define(
		"prediction", 
		NeutrinoWeightingWrapper, 
	//	{"lepton_lvec_for_truth", "whad_lvec_for_truth", "met_value_for_truth", "met_phi_for_truth"}
		{"reco_lepton_lvec", "reco_whad_lvec", "reco_met_value", "reco_met_phi"}
	);

	rLoopManager = rLoopManager.Define(
		"NW_weight", 
		ExtractWeight, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_weight_sum", 
		ExtractWeightSum, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_H_smear", 
		ExtractHiggsSmear, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_nu_eta", 
		ExtractEstimatedEtaNu, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_wlep_mass", 
		ExtractEstimatedMassWLep, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_nu_px", 
		ExtractEstimatedPx, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_nu_py", 
		ExtractEstimatedPy, 
		{"prediction"}
	);
	rLoopManager = rLoopManager.Define(
		"NW_nu_phi", 
		ExtractEstimatedPhiNu, 
		{"prediction"}
	);




	//>> full NW information for 2D plots
	//rLoopManager = rLoopManager.Define(
	//	"prediction_full",
	//	NeutrinoWeightingWrapperFull,
	////	{"lepton_lvec_for_truth", "whad_lvec_for_truth", "met_value_for_truth", "met_phi_for_truth"}
	//	{"reco_lepton_lvec", "reco_whad_lvec", "reco_met_value", "reco_met_phi"}
	//);
//
	//rLoopManager = rLoopManager.Define(
	//	"NW_full_weight",
	//	ExtractWeightFull,
	//	{"prediction_full"}
	//);
	//rLoopManager = rLoopManager.Define(
	//	"NW_full_wlep_mass",
	//	ExtractWlepMassesFull,
	//	{"prediction_full"}
	//);
	//rLoopManager = rLoopManager.Define(
	//	"NW_full_nu_eta",
	//	ExtractNuEtaFull,
	//	{"prediction_full"}
	//);
	
	


	///----------  Output  ----------///
	//>> save snapshot
	cout << " > saving snapshot to " << OUTPUT_FILE << endl;
	rLoopManager.Snapshot(
		"neutrino_weighting", 
		OUTPUT_FILE,
		OUTPUT_COLOUMN_NAMES
	);

	return 0;
}