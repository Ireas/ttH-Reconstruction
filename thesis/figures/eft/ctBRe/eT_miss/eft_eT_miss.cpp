#include "TFile.h"
#include <string>


int eft_eT_miss(){

	// Define input/output filepath
	std::string pathInput = "/mnt/c9e9c9e1-af0b-4d2a-854e-cf170c26e797/OutputNtuplesFull/ttZ_3L_v4.1_noZ/";
	std::string pathOutput = "/home/ireas/Desktop/eft_analysis/ctBRe/eT_miss/";

	int index_default = 0;
	int index_c1 = 44;  //p0.3
	int index_c2 = 2;   //m0.3
	int index_c3 = 55;  //p0.9
	int index_c4 = 920; //0.9



	// Collect all trees into Chains
	TChain c_smeft("nominal"); //Data Chain
	c_smeft.Add( (pathInput+"exported_508772.MGPy8_ttll_SMEFTsim_reweighted.ttZ_tL_v4.e8379_a875_r9364_p4514.root").c_str() ); 
	c_smeft.Add( (pathInput+"exported_508772.MGPy8_ttll_SMEFTsim_reweighted.ttZ_tL_v4.e8379_a875_r10201_p4514.root").c_str() ); 
	c_smeft.Add( (pathInput+"exported_508772.MGPy8_ttll_SMEFTsim_reweighted.ttZ_tL_v4.e8379_a875_r10724_p4514.root").c_str() ); 
	//c_smeft.Add( (pathInput+"exported_508773.MGPy8_ttll_SMEFTsim_reweighted_prop.ttZ_tL_v4.e8379_a875_r9364_p4514.root").c_str() ); 
	//c_smeft.Add( (pathInput+"exported_508773.MGPy8_ttll_SMEFTsim_reweighted_prop.ttZ_tL_v4.e8379_a875_r10201_p4514.root").c_str() ); 
	//c_smeft.Add( (pathInput+"exported_508773.MGPy8_ttll_SMEFTsim_reweighted_prop.ttZ_tL_v4.e8379_a875_r10724_p4514.root").c_str() ); 



	// Initiate Variables
	bool HasOSSFPair = false;
	bool tight_1lep = false;
	bool tight_2lep = false;
	bool tight_3lep = false;
	Float_t pT_1lep = 0.0;
	Float_t pT_2lep = 0.0;
	Float_t pT_3lep = 0.0;
	Float_t min_mll = 0.0;
	Int_t nJets = 0;	
	Int_t nEl = 0;
	Int_t nMu = 0;
	
	Float_t MyModel_ttZ = 0.0;
	Float_t MyModel_ttW = 0.0;
	Float_t MyModel_Fakes = 0.0;
	
	Float_t eT_miss = 0.0;
	
	Double_t XSecWeight = 0.0;
	Float_t weight_mc = 0.0;
	Float_t weight_pileup = 0.0;
	Float_t weight_leptonSF = 0.0;
	Float_t weight_bTagSF_DL1r_Continuous = 0.0;
	Float_t weight_jvt = 0.0;
	Float_t weight_year = 0.0;
	Float_t weight_photonSF = 0.0;
	Float_t integrated_luminosity = 138965.2; //Constant
	
	Double_t weight_nominal = 0.0;
	Double_t weight_default = 0.0; 
	Double_t weight_c1 = 0.0;
	Double_t weight_c2 = 0.0;
	Double_t weight_c3 = 0.0;
	Double_t weight_c4 = 0.0;
	


	// Initiate MC  Histogram
	std::vector<TH1F*> hist_array_ttZ; 
	TH1F* hist_ttZ_default = new TH1F("Region00d", "ttZ_default", 3,0,140);
	TH1F* hist_ttZ_c1 = new TH1F("Region00c1", "ttZ_c1", 3,0,140);
	TH1F* hist_ttZ_c2 = new TH1F("Region00c2", "ttZ_c2", 3,0,140);
	TH1F* hist_ttZ_c3 = new TH1F("Region00c3", "ttZ_c3", 3,0,140);
	TH1F* hist_ttZ_c4 = new TH1F("Region00c4", "ttZ_c4", 3,0,140);
	hist_array_ttZ.push_back(hist_ttZ_default);
	hist_array_ttZ.push_back(hist_ttZ_c1);
	hist_array_ttZ.push_back(hist_ttZ_c2);
	hist_array_ttZ.push_back(hist_ttZ_c3);
	hist_array_ttZ.push_back(hist_ttZ_c4);

	std::vector<TH1F*> hist_array_ttW; 
	TH1F* hist_ttW_default = new TH1F("Region01d", "ttW_default", 3,0,140);
	TH1F* hist_ttW_c1 = new TH1F("Region01c1", "ttW_c1", 3,0,140);
	TH1F* hist_ttW_c2 = new TH1F("Region01c2", "ttW_c2", 3,0,140);
	TH1F* hist_ttW_c3 = new TH1F("Region01c3", "ttW_c3", 3,0,140);
	TH1F* hist_ttW_c4 = new TH1F("Region01c4", "ttW_c4", 3,0,140);
	hist_array_ttW.push_back(hist_ttW_default);
	hist_array_ttW.push_back(hist_ttW_c1);
	hist_array_ttW.push_back(hist_ttW_c2);
	hist_array_ttW.push_back(hist_ttW_c3);
	hist_array_ttW.push_back(hist_ttW_c4);
	
	std::vector<TH1F*> hist_array_Fakes; 
	TH1F* hist_Fakes_default = new TH1F("Region02d", "Fakes_default", 3,0,140);
	TH1F* hist_Fakes_c1 = new TH1F("Region02c1", "Fakes_c1", 3,0,140);
	TH1F* hist_Fakes_c2 = new TH1F("Region02c2", "Fakes_c2", 3,0,140);
	TH1F* hist_Fakes_c3 = new TH1F("Region02c3", "Fakes_c3", 3,0,140);
	TH1F* hist_Fakes_c4 = new TH1F("Region02c4", "Fakes_c4", 3,0,140);
	hist_array_Fakes.push_back(hist_Fakes_default);
	hist_array_Fakes.push_back(hist_Fakes_c1);
	hist_array_Fakes.push_back(hist_Fakes_c2);
	hist_array_Fakes.push_back(hist_Fakes_c3);
	hist_array_Fakes.push_back(hist_Fakes_c4);
	
	//Set Branch Adress to MC Variable
	c_smeft.SetBranchAddress("HasOSSFPair", &HasOSSFPair);
	c_smeft.SetBranchAddress("tight_1lep", &tight_1lep);
	c_smeft.SetBranchAddress("tight_2lep", &tight_2lep);
	c_smeft.SetBranchAddress("tight_3lep", &tight_3lep);
	c_smeft.SetBranchAddress("pT_1lep", &pT_1lep);
	c_smeft.SetBranchAddress("pT_2lep", &pT_2lep);
	c_smeft.SetBranchAddress("pT_3lep", &pT_3lep);
	c_smeft.SetBranchAddress("min_mll", &min_mll);
	c_smeft.SetBranchAddress("nJets", &nJets);
	c_smeft.SetBranchAddress("nEl", &nEl);
	c_smeft.SetBranchAddress("nMu", &nMu);

	c_smeft.SetBranchAddress("MyModel_ttZ", &MyModel_ttZ);
	c_smeft.SetBranchAddress("MyModel_ttW", &MyModel_ttW);
	c_smeft.SetBranchAddress("MyModel_Fakes", &MyModel_Fakes);

	c_smeft.SetBranchAddress("eT_miss", &eT_miss);

	c_smeft.SetBranchAddress("XSecWeight", &XSecWeight);
	c_smeft.SetBranchAddress("weight_mc", &weight_mc);
	c_smeft.SetBranchAddress("weight_pileup", &weight_pileup);
	c_smeft.SetBranchAddress("weight_leptonSF", &weight_leptonSF);
	c_smeft.SetBranchAddress("weight_bTagSF_DL1r_Continuous", &weight_bTagSF_DL1r_Continuous);
	c_smeft.SetBranchAddress("weight_jvt", &weight_jvt);
	c_smeft.SetBranchAddress("weight_year", &weight_year);
	c_smeft.SetBranchAddress("weight_photonSF", &weight_photonSF);

	std::vector<float>* mc_generator_weights = new std::vector<float>; 
	c_smeft.SetBranchAddress("mc_generator_weights", &mc_generator_weights);
	
	//Event Loop
	for(int i=0; i<c_smeft.GetEntries(); i++){
		c_smeft.GetEntry(i);
	
		//Standard Selection
		if( !(nJets>=3) )
			continue;
		if( !HasOSSFPair )
			continue;
		if( !tight_1lep )
			continue;
		if( !tight_2lep )
			continue;
		if( !tight_3lep )
			continue;
		if( !(pT_1lep>27) )
			continue;
		if( !(pT_2lep>20) )
			continue;
		if( !(pT_3lep>15) )
			continue;
		if( !((nEl+nMu)==3) )
			continue;
		if( !(min_mll>10) )
			continue;
		 

		
		//Calculate Weights
		weight_nominal = XSecWeight * weight_mc * weight_pileup * weight_leptonSF * weight_bTagSF_DL1r_Continuous * weight_jvt * weight_year * weight_photonSF * integrated_luminosity; 
		weight_default = weight_nominal * mc_generator_weights->at(index_default) / weight_mc;
		weight_c1 = weight_nominal * mc_generator_weights->at(index_c1) / weight_mc;
		weight_c2 = weight_nominal * mc_generator_weights->at(index_c2) / weight_mc;
		weight_c3 = weight_nominal * mc_generator_weights->at(index_c3) / weight_mc;
		weight_c4 = weight_nominal * mc_generator_weights->at(index_c4) / weight_mc;

	
		//Fill the Histogram
		if( (MyModel_ttW<0.4) && (MyModel_Fakes<0.4) ){		
			(hist_array_ttZ[0])->Fill(eT_miss, weight_default);
			(hist_array_ttZ[1])->Fill(eT_miss, weight_c1);
			(hist_array_ttZ[2])->Fill(eT_miss, weight_c2);
			(hist_array_ttZ[3])->Fill(eT_miss, weight_c3);
			(hist_array_ttZ[4])->Fill(eT_miss, weight_c4);
		}
		if( (MyModel_ttW>=0.4) ){
			(hist_array_ttW[0])->Fill(eT_miss, weight_default);
			(hist_array_ttW[1])->Fill(eT_miss, weight_c1);
			(hist_array_ttW[2])->Fill(eT_miss, weight_c2);
			(hist_array_ttW[3])->Fill(eT_miss, weight_c3);
			(hist_array_ttW[4])->Fill(eT_miss, weight_c4);
		}
		if( (MyModel_ttW<0.4) && (MyModel_Fakes>=0.4) ){		
			(hist_array_Fakes[0])->Fill(eT_miss, weight_default);
			(hist_array_Fakes[1])->Fill(eT_miss, weight_c1);
			(hist_array_Fakes[2])->Fill(eT_miss, weight_c2);
			(hist_array_Fakes[3])->Fill(eT_miss, weight_c3);
			(hist_array_Fakes[4])->Fill(eT_miss, weight_c4);
		}

	}
	
	// Print histograms into output.root
	TFile* fileOutput = new TFile( (pathOutput+"output.root").c_str(), "RECREATE" );
	fileOutput->cd();
		(hist_array_ttZ[0])->Write("ttZ_default");	
		(hist_array_ttZ[1])->Write("ttZ_c1");	
		(hist_array_ttZ[2])->Write("ttZ_c2");	
		(hist_array_ttZ[3])->Write("ttZ_c3");	
		(hist_array_ttZ[4])->Write("ttZ_c4");	
		(hist_array_ttW[0])->Write("ttW_default");	
		(hist_array_ttW[1])->Write("ttW_c1");	
		(hist_array_ttW[2])->Write("ttW_c2");	
		(hist_array_ttW[3])->Write("ttW_c3");	
		(hist_array_ttW[4])->Write("ttW_c4");	
		(hist_array_Fakes[0])->Write("Fakes_default");	
		(hist_array_Fakes[1])->Write("Fakes_c1");	
		(hist_array_Fakes[2])->Write("Fakes_c2");	
		(hist_array_Fakes[3])->Write("Fakes_c3");	
		(hist_array_Fakes[4])->Write("Fakes_c4");	
	fileOutput->Close();
	

	//Normalise
	for(int k=0; k<5; k++){
		Float_t factorA = (hist_array_ttZ[k])->Integral();
		if(factorA!=0){
			(hist_array_ttZ[k])->Scale(1/factorA);
		}

		Float_t factorB = (hist_array_ttW[k])->Integral();
		if(factorB!=0){
			(hist_array_ttW[k])->Scale(1/factorB);
		}
		
		Float_t factorC = (hist_array_Fakes[k])->Integral();
		if(factorC!=0){
			(hist_array_Fakes[k])->Scale(1/factorC);
		}
	}
	// Print nominalised histograms into output_nominalised.root
	fileOutput = new TFile( (pathOutput+"output_nominalised.root").c_str(), "RECREATE" );
	fileOutput->cd();
		(hist_array_ttZ[0])->Write("ttZ_default");	
		(hist_array_ttZ[1])->Write("ttZ_c1");	
		(hist_array_ttZ[2])->Write("ttZ_c2");	
		(hist_array_ttZ[3])->Write("ttZ_c3");	
		(hist_array_ttZ[4])->Write("ttZ_c4");	
		(hist_array_ttW[0])->Write("ttW_default");	
		(hist_array_ttW[1])->Write("ttW_c1");	
		(hist_array_ttW[2])->Write("ttW_c2");	
		(hist_array_ttW[3])->Write("ttW_c3");	
		(hist_array_ttW[4])->Write("ttW_c4");	
		(hist_array_Fakes[0])->Write("Fakes_default");	
		(hist_array_Fakes[1])->Write("Fakes_c1");	
		(hist_array_Fakes[2])->Write("Fakes_c2");	
		(hist_array_Fakes[3])->Write("Fakes_c3");	
		(hist_array_Fakes[4])->Write("Fakes_c4");	
	fileOutput->Close();
	

	// Clean up
	c_smeft.ResetBranchAddresses();
	
	delete(mc_generator_weights);
	delete(hist_array_ttZ[0]);
	delete(hist_array_ttZ[1]);
	delete(hist_array_ttZ[2]);
	delete(hist_array_ttZ[3]);
	delete(hist_array_ttZ[4]);
	delete(hist_array_ttW[0]);
	delete(hist_array_ttW[1]);
	delete(hist_array_ttW[2]);
	delete(hist_array_ttW[3]);
	delete(hist_array_ttW[4]);
	delete(hist_array_Fakes[0]);
	delete(hist_array_Fakes[1]);
	delete(hist_array_Fakes[2]);
	delete(hist_array_Fakes[3]);
	delete(hist_array_Fakes[4]);
	delete(fileOutput);

	return 0;
}
