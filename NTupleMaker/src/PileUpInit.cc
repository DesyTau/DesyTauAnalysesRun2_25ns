#include "DesyTauAnalyses/NTupleMaker/interface/PileUpInit.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h" 

PileUpSyst* PileUpInit::initPU (Config cfg ){


	std::string runPeriod = cfg.get<std::string>("RunPeriod");
	
	std::string RootFilePath = "/nfs/dust/cms/user/bottav/CMSSW_7_4_7_patch1/src/DesyTauAnalyses/NTupleMaker/test/PileUp/";
	std::string RootFileName = "PileUpDistributions_"+runPeriod;

	//const std::string sampleName = cfg.get<std::string>("SampleName");      
	
	//TH1::AddDirectory(kFALSE);
	//TDirectory* dir = gDirectory;
	TFile * inputRootFile = new TFile(TString(RootFilePath)+TString(RootFileName)+".root","read");
	//dir->cd();

	// configuration of PUsyst
		// 1) configure PU central
 	PileUp * PUcentral = new PileUp();
	TH1D * h_data_central = new TH1D(); h_data_central = (TH1D*)(inputRootFile->Get("data_central")->Clone());
	TH1D * h_MC_central   = new TH1D(); h_MC_central = (TH1D*)(inputRootFile->Get("Dummy_MC")->Clone());  // --------> fix to get w/ sampleName
	PUcentral->set_h_data(h_data_central);
	PUcentral->set_h_MC(h_MC_central);

		// 2) get shifted histograms
	TH1D * h_data_up = new TH1D();   h_data_up = (TH1D*)(inputRootFile->Get("data_systUp")->Clone());
	TH1D * h_data_down = new TH1D(); h_data_down = (TH1D*)(inputRootFile->Get("data_systDown")->Clone());	

		// 3) create PUsyst 
	PileUpSyst * PUsyst = new PileUpSyst(PUcentral);
	PUsyst->set_shifted_data(h_data_up, "up");
	PUsyst->set_shifted_data(h_data_down, "down");

	inputRootFile->Close();
	return PUsyst;

} 


// runPeriod = latest
// SampleName = Dummy_MC

