#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TKey.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Spring15Tree.h"

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#define pi 	3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 		   0.10565837
#define tauMass 		   1.77682
#define pionMass 		   0.1396

int main(int argc, char * argv[]){

  // first argument - config file for analysis
  // second argument - file
  
  using namespace std;

  string cmsswBase = (getenv ("CMSSW_BASE"));
  
  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load( pathToCrystalLib );
  if (openSuccessful !=0 ) {
    cout<<pathToCrystalLib<<" does not exists. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit( -1 );
  }

  // **** configuration

  Config cfg(argv[1]);
  const string inName = argv[2];
  
  // No Data
  const bool isData = cfg.get<bool>("isData");
  if (isData) exit(0);

  // which SF
  const bool ApplyPUweight   = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepIdIsoSF = cfg.get<bool>("ApplyLepIdIsoSF"); 
  const bool ApplyTriggerSF  = cfg.get<bool>("ApplyTriggerSF");
  const bool ApplyTrackSF    = cfg.get<bool>("ApplyTrackSF");

  // PU reweighting - initialization
  PileUp * PUofficial = new PileUp();
  if(ApplyPUweight){
    const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
    const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
    
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(pileUpInDataFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(pileUpInMCFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  

  // Lepton Trigger and Id+Iso SF
  ScaleFactor * SF_lepIdIso = new ScaleFactor();
  ScaleFactor * SF_lepTrigger = new ScaleFactor();

  if(ApplyLepIdIsoSF){
    const string idIsoEffFile = cfg.get<string>("idIsoEffFile");
    SF_lepIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(idIsoEffFile));
  }
  
  if(ApplyTriggerSF){  
    const string trigEffFile = cfg.get<string>("trigEffFile");
    SF_lepTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(trigEffFile));
  }
  
  // tau ID SF and Track SF
  TString workspace_filename = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v4.root";
  TFile *f_workspace = new TFile(workspace_filename,"read");
  if (f_workspace->IsZombie()) {std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");

  // output fileName with histograms
  TFile* ifile = new TFile( TString(inName) ,"read");
  
  TString outName(inName);
  outName.ReplaceAll(".root","_updated.root");
  TFile *ofile = new TFile( outName ,"recreate");
  ofile->cd("");
  
  TIter next(ifile->GetListOfKeys());
  TKey *key;	
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree")) continue;
    TTree *itree = (TTree*)key->ReadObj();
    itree->LoadTree(0);
    itree->SetBranchStatus("*", 1);

    Float_t pt_1_in = 0.;
    Float_t eta_1_in = 0.;
    Float_t npu_in = 0.;
    Float_t m_1_in = 0.;
    Int_t gen_match_2_in = 0;
    Float_t pu_weight_in = 0.;
    Float_t trigweight_1_in = 0.;
    Float_t idisoweight_1_in = 0.;
    Double_t trkeffweight_1_in = 0.;
    Float_t trigweight_2_in = 0.;
    Float_t idisoweight_2_in = 0.;

    itree->SetBranchAddress("pt_1",&pt_1_in);
    itree->SetBranchAddress("eta_1",&eta_1_in);
    itree->SetBranchAddress("npu",&npu_in);
    itree->SetBranchAddress("m_1",&m_1_in);
    itree->SetBranchAddress("gen_match_2",&gen_match_2_in);
    itree->SetBranchAddress("pu_weight", &pu_weight_in);
    itree->SetBranchAddress("trigweight_1",&trigweight_1_in);
    itree->SetBranchAddress("idisoweight_1",&idisoweight_1_in);
    itree->SetBranchAddress("trkeffweight_1",&trkeffweight_1_in);
    itree->SetBranchAddress("trigweight_2",&trigweight_2_in);  
    itree->SetBranchAddress("idisoweight_2",&idisoweight_2_in);
    
    TTree *otree = itree->GetTree()->CloneTree(0);
    Float_t pu_weight_out = 0.;
    Float_t trigweight_1_out = 0.;
    Float_t idisoweight_1_out = 0.;
    Double_t trkeffweight_1_out = 0.;
    Float_t trigweight_2_out = 0.;
    Float_t idisoweight_2_out = 0.;
    Float_t effweight_out = 0.;
    
    otree->SetBranchAddress("pu_weight", &pu_weight_out);
    otree->SetBranchAddress("trigweight_1",&trigweight_1_out);
    otree->SetBranchAddress("idisoweight_1",&idisoweight_1_out);
    otree->SetBranchAddress("trkeffweight_1",&trkeffweight_1_out);
    otree->SetBranchAddress("trigweight_2",&trigweight_2_out);  
    otree->SetBranchAddress("idisoweight_2",&idisoweight_2_out);
    otree->SetBranchAddress("effweight",&effweight_out);
    
    Long64_t numberOfEntries = itree->GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;

    Int_t nEvents = 0;
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {       
      itree->GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%100000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
      
      pu_weight_out = pu_weight_in;
      trigweight_1_out = trigweight_1_in;
      idisoweight_1_out = idisoweight_1_in;
      trigweight_2_out = trigweight_2_in;   
      idisoweight_2_out = idisoweight_2_in;
      trkeffweight_1_out = trkeffweight_1_in;
      
      // PU weight
      if (!isData && ApplyPUweight) 
	pu_weight_out = float(PUofficial->get_PUweight(double(npu_in)));
      
      // Scale Factor trigger and Id+Iso
      if (!isData && ApplyTriggerSF)
	trigweight_1_out = SF_lepTrigger->get_EfficiencyData(pt_1_in, eta_1_in);
      
      if (!isData && ApplyLepIdIsoSF)    
	idisoweight_1_out = SF_lepIdIso->get_ScaleFactor(pt_1_in, eta_1_in);
      
      // tracking efficiency weight
      if (!isData && ApplyTrackSF) {      
	if(m_1_in > 0.1 && m_1_in < 0.11) {
	  w->var("m_eta")->setVal(eta_1_in); 
	  trkeffweight_1_out = w->function("m_trk_ratio")->getVal();
	}
	else if( m_1_in > 0. && m_1_in < 0.001){
	  w->var("e_eta")->setVal(pt_1_in); 
	  w->var("e_pt")->setVal(eta_1_in);
	  trkeffweight_1_out = w->function("e_trk_ratio")->getVal();
	}
      }
      
      // tauID weight
      if (!isData && gen_match_2_in==5)
	idisoweight_2_out = 0.83;
      
      effweight_out = (trigweight_1_out)*(idisoweight_1_out)*(trigweight_2_out)*(idisoweight_2_out);
      
      otree->Fill();
    } // end of file processing (loop over events in one file)
    
    std::cout << std::endl;
    std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
    std::cout << std::endl;
    
    ofile->cd("");
    otree->Write(); 
  }
  ofile->Close(); 
}
