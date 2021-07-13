#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string> 
#include <map>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/src/Config.cc"
#include "TRandom.h"
#include "TRandom3.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TSystem.h"

using namespace std;

const float PionMass = 0.13957;

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : analysis_macro [config] [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(1);
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// ************* Configuration ********************* //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool applyHiggsPtWeight = cfg.get<bool>("ApplyHiggsPtWeight");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");

  // topological cuts
  const float dRMuonsCut     = cfg.get<float>("dRMuonsCut");
  const bool sameSign        = cfg.get<bool>("SameSignMuons");

  const bool applyOSdimuonVeto = cfg.get<bool>("ApplyOppositeSignDimuonVeto");
  const float maxPtSumOSdimuons   = cfg.get<float>("MaxPtOppositeSignDimuon");

  // track selection
  const float ptSumCut        = cfg.get<float>("ptSumCut");
  const float dRIso           = cfg.get<float>("dRIso");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile            = cfg.get<string>("jsonFile");

  // trigger
  const bool applyTriggerMatch          = cfg.get<bool>("ApplyTriggerMatch");
  const string dimuonTriggerName        = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName     = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName1     = cfg.get<string>("MuonLowPtFilterName1");
  const string muonLowPtFilterName2     = cfg.get<string>("MuonLowPtFilterName2");
  const string dimuonDzFilterName       = cfg.get<string>("DimuonDzFilterName");
  const string dimuonSameSignFilterName = cfg.get<string>("DimuonSameSignFilterName");
  
  // trigger matching
  const float DRTrigMatch          = cfg.get<float>("DRTrigMatch"); 
  const float effDzSS              = cfg.get<float>("effDzSS");
  const float trkIsoSF             = cfg.get<float>("TrackIsolationSF");
  const unsigned int numberOfMuons = cfg.get<unsigned int>("NumberOfMuons");

  TString DiMuonTriggerName(dimuonTriggerName);
  TString MuonHighPtFilterName(muonHighPtFilterName);
  TString MuonLowPtFilterName1(muonLowPtFilterName1);
  TString MuonLowPtFilterName2(muonLowPtFilterName2);
  TString DiMuonDzFilterName(dimuonDzFilterName);
  TString DiMuonSameSignFilterName(dimuonSameSignFilterName);

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile   = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  const string Muon17TriggerFile = cfg.get<string>("MuonHighPtTriggerEff");
  const string Muon8TriggerFile  = cfg.get<string>("MuonLowPtTriggerEff");
  const string CorrectionWSFileName = cfg.get<string>("CorrectionWorkspaceFileName");

  // Higgs pt reweighting
  const string higgsPtFileName = cfg.get<string>("HiggsPtFileName");
  TString HiggsPtFileName(higgsPtFileName);
  const bool isVH = cfg.get<bool>("IsVH");

  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  std::ifstream fileList(argv[2]);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// ************* NTuple Branches ********************* //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // event info
  ULong64_t event_nr;
  unsigned int event_run;
  unsigned int event_luminosityblock;

  // tracks 
  UInt_t track_count;
  int track_ID[1000];
  float track_px[1000];
  float track_py[1000];
  float track_pz[1000];
  float track_pt[1000];
  float track_eta[1000];
  float track_phi[1000];
  float track_charge[1000];
  float track_mass[1000];
  float track_dxy[1000];
  float track_dxyerr[1000];
  float track_dz[1000];
  float track_dzerr[1000];
  bool track_highPurity[1000];
  
  // muons
  UInt_t muon_count;
  UInt_t muon_nMuonStations[1000];
  UInt_t muon_nMuonHits[1000];
  UInt_t muon_nPixelHits[1000];
  UInt_t muon_nTrackerHits[1000];
  float muon_px[1000];
  float muon_py[1000];
  float muon_pz[1000];
  float muon_pt[1000];
  float muon_eta[1000];
  float muon_phi[1000];
  float muon_pterror[1000];
  float muon_chi2[1000];
  float muon_ndof[1000];
  float muon_charge[1000];
  float muon_dxy[1000];
  float muon_dxyerr[1000];
  float muon_dz[1000];
  float muon_dzerr[1000];
  float muon_chargedHadIso[1000];
  float muon_neutralHadIso[1000];
  float muon_photonIso[1000];
  float muon_puIso[1000];
  bool muon_isPF[1000];
  bool muon_isGlobal[1000];
  bool muon_isTracker[1000];
  bool muon_isTight[1000];
  bool muon_isLoose[1000];
  bool muon_isMedium[1000];
  bool muon_isICHEP[1000];

  UInt_t genparticles_count;
  Float_t genparticles_e[1000];
  Float_t genparticles_px[1000];
  Float_t genparticles_py[1000];
  Float_t genparticles_pz[1000];
  Int_t genparticles_pdgid[1000];
  Int_t genparticles_status[1000];
  UInt_t genparticles_info[1000];
  Int_t genparticles_fromHardProcess[1000];

  float genweight;

  float metx;
  float mety;
  //  float met;
  //  float metphi;
  
  // Trigger
  unsigned int trigobject_count;
  float trigobject_px[1000];
  float trigobject_py[1000];
  float trigobject_pz[1000];
  float trigobject_pt[1000];
  float  trigobject_eta[1000];
  float trigobject_phi[1000];
  bool trigobject_filters[1000][200];

  // Gen jets
  UInt_t genjets_count;
  Float_t genjets_e[100];
  Float_t genjets_px[100];
  Float_t genjets_py[100];
  Float_t genjets_pz[100];
  Float_t genjets_pt[100];
  Float_t genjets_eta[100];
  Float_t genjets_phi[100];
  Int_t genjets_pdgid[100];
  Int_t genjets_status[100];

  // Jets
  unsigned int pfjet_count;
  float pfjet_e[200];
  float pfjet_px[200];
  float pfjet_py[200];
  float pfjet_pz[200];
  float pfjet_pt[200];
  float pfjet_eta[200];
  float pfjet_phi[200];
  int pfjet_flavour[200];
  float pfjet_btag[200][10];
  Bool_t pfjet_pu_jet_fullId_loose[100];
  Bool_t pfjet_pu_jet_fullId_medium[100];
  Bool_t pfjet_pu_jet_fullId_tight[100];

  float numtruepileupinteractions;

  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  

  std::string rootFileName(argv[2]);  
  std::string chainName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;
  if (TStrName.Contains("Signal")) {
    std::cout << "=============================" << std::endl;
    std::cout << "=== Running on Signal MC ====" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
  }

  TString FullName = TStrName;      
  
  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// ************* HISTOGRAMS & TREES ********************* //////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////
  ////////// Variables for Dataset /////////
  //////////////////////////////////////////
  
  float MuLMuT_Mass, MuLMuT_DR, MuLMuT_DPhi, MuLTrk_Mass, MuLTrk_Pt, // 5
        MuLTrk_DR, MuLTrkMET_DPhi, MuTTrk_Mass, MuTTrk_Pt, MuTTrk_DR, // 5
        MuTTrkMET_DPhi, MuLTrkMuTTrk_Mass, MuLTrkMuTTrk_Pt, MET_Pt, MuLTrkMuTTrkMET_Mass; // 5
  float Eventweight; // 1
  float Eventweight_TrkIso_Up,Eventweight_TrkIso_Down;
  
  //////////////////////////////////////////
  /////////////// Histograms ///////////////
  //////////////////////////////////////////
  
  /////////////////// Counters /////////////////////
  TH1D * counter_InputEventsH=new TH1D("counter_InputEventsH","",1,0.,2.);
  TH1D * counter_MuonSizeGTE2H=new TH1D("counter_MuonSizeGTE2H","",1,0.,2.);
  TH1D * counter_MuonKinematicsH=new TH1D("counter_MuonKinematicsH","",1,0.,2.);         
  TH1D * counter_FinalEventsH=new TH1D("counter_FinalEventsH","",1,0.,2.);         
  
  /////// Muons, Trigger and PU ////////////////////
  TH1D * muonCountH = new TH1D("muonCountH","",11,-0.5,10.5);
  TH1D * nGoodMuonsH = new TH1D("nGoodMuonsH","",11,-0.5,10.5);
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);
  TH1D * MuTrkWeightH = new TH1D("MuTrkWeightH","",80,0.8,1.2);
  TH1D * HiggsPtWeightH = new TH1D("HiggsPtWeightH","",250,0,5);
  
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);
  
  //////// Kinematics Muons ////////////////////
  TH1D * ptLeadingMuH = new TH1D("ptLeadingMuH","",400,0,400);
  TH1D * ptTrailingMuH = new TH1D("ptTrailingMuH","",400,0,400);
  TH1D * etaLeadingMuH = new TH1D("etaLeadingMuH","",48,-2.4,2.4);
  TH1D * etaTrailingMuH = new TH1D("etaTrailingMuH","",48,-2.4,2.4);
  TH1D * dimuonMassH = new TH1D("dimuonMassH","",500,0,500);
  
  //////// Track Multiplicity ////////////////////
  TH1D * nTracksLeadingMuH = new TH1D("nTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nTracksTrailingMuH = new TH1D("nTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nTracksTrkLeadingMuH = new TH1D("nTracksTrkLeadingMuH","",21,-0.5,20.5);
  TH1D * nTracksTrkTrailingMuH = new TH1D("nTracksTrkTrailingMuH","",21,-0.5,20.5);
  
  TH1D * nSoftTracksLeadingMuH = new TH1D("nSoftTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksTrailingMuH = new TH1D("nSoftTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksTrkLeadingMuH = new TH1D("nSoftTracksTrkLeadingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksTrkTrailingMuH = new TH1D("nSoftTracksTrkTrailingMuH","",21,-0.5,20.5);
  
  ////////////////// Signal Region //////////////////
  TTree * tree_Sel = new TTree("tree_Sel","tree_Sel");
  tree_Sel->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_Sel->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_Sel->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_Sel->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_Sel->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_Sel->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_Sel->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_Sel->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_Sel->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_Sel->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_Sel->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_Sel->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_Sel->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_Sel->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_Sel->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_Sel->Branch("Eventweight",&Eventweight,"Eventweight");
  tree_Sel->Branch("Eventweight_TrkIso_Up",&Eventweight_TrkIso_Up,"Eventweight_TrkIso_Up");
  tree_Sel->Branch("Eventweight_TrkIso_Down",&Eventweight_TrkIso_Down,"Eventweight_TrkIso_Down");
  
  
  ////// Control Regions for Bkgd Estimation ////////
  TTree * tree_NoSR = new TTree("tree_NoSR","tree_NoSR");
  tree_NoSR->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_NoSR->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_NoSR->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_NoSR->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_NoSR->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_NoSR->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_NoSR->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_NoSR->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_NoSR->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_NoSR->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_NoSR->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_NoSR->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_NoSR->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_NoSR->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_NoSR->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_NoSR->Branch("Eventweight",&Eventweight,"Eventweight");
  
  TTree * tree_SemiIso = new TTree("tree_SemiIso","tree_SemiIso");
  tree_SemiIso->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_SemiIso->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_SemiIso->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_SemiIso->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_SemiIso->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_SemiIso->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_SemiIso->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_SemiIso->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_SemiIso->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_SemiIso->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_SemiIso->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_SemiIso->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_SemiIso->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_SemiIso->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_SemiIso->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_SemiIso->Branch("Eventweight",&Eventweight,"Eventweight");
  
  TTree * tree_LooseIso = new TTree("tree_LooseIso","tree_LooseIso");
  tree_LooseIso->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_LooseIso->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_LooseIso->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_LooseIso->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_LooseIso->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_LooseIso->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_LooseIso->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_LooseIso->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_LooseIso->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_LooseIso->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_LooseIso->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_LooseIso->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_LooseIso->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_LooseIso->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_LooseIso->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_LooseIso->Branch("Eventweight",&Eventweight,"Eventweight");
  
  TTree * tree_LooseSemiIso = new TTree("tree_LooseSemiIso","tree_LooseSemiIso");
  tree_LooseSemiIso->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_LooseSemiIso->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_LooseSemiIso->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_LooseSemiIso->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_LooseSemiIso->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_LooseSemiIso->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_LooseSemiIso->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_LooseSemiIso->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_LooseSemiIso->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_LooseSemiIso->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_LooseSemiIso->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_LooseSemiIso->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_LooseSemiIso->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_LooseSemiIso->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_LooseSemiIso->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_LooseSemiIso->Branch("Eventweight",&Eventweight,"Eventweight");
  
  TTree * tree_LeadingSemiIso = new TTree("tree_LeadingSemiIso","tree_LeadingSemiIso");
  tree_LeadingSemiIso->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_LeadingSemiIso->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_LeadingSemiIso->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_LeadingSemiIso->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_LeadingSemiIso->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_LeadingSemiIso->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_LeadingSemiIso->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_LeadingSemiIso->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_LeadingSemiIso->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_LeadingSemiIso->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_LeadingSemiIso->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_LeadingSemiIso->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_LeadingSemiIso->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_LeadingSemiIso->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_LeadingSemiIso->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_LeadingSemiIso->Branch("Eventweight",&Eventweight,"Eventweight");
  
  TTree * tree_LeadingLooseIso = new TTree("tree_LeadingLooseIso","tree_LeadingLooseIso");
  tree_LeadingLooseIso->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_LeadingLooseIso->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_LeadingLooseIso->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_LeadingLooseIso->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_LeadingLooseIso->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_LeadingLooseIso->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_LeadingLooseIso->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_LeadingLooseIso->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_LeadingLooseIso->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_LeadingLooseIso->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_LeadingLooseIso->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_LeadingLooseIso->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_LeadingLooseIso->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_LeadingLooseIso->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_LeadingLooseIso->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_LeadingLooseIso->Branch("Eventweight",&Eventweight,"Eventweight");
  
  
  ////// Generator Level Distributions ////////
  TTree * tree_GenLev = new TTree("tree_GenLev","tree_GenLev");
  tree_GenLev->Branch("MuLMuT_Mass",&MuLMuT_Mass,"MuLMuT_Mass");
  tree_GenLev->Branch("MuLMuT_DR",&MuLMuT_DR,"MuLMuT_DR");
  tree_GenLev->Branch("MuLMuT_DPhi",&MuLMuT_DPhi,"MuLMuT_DPhi");
  tree_GenLev->Branch("MuLTrk_Mass",&MuLTrk_Mass,"MuLTrk_Mass");
  tree_GenLev->Branch("MuLTrk_Pt",&MuLTrk_Pt,"MuLTrk_Pt");
  tree_GenLev->Branch("MuLTrk_DR",&MuLTrk_DR,"MuLTrk_DR");
  tree_GenLev->Branch("MuLTrkMET_DPhi",&MuLTrkMET_DPhi,"MuLTrkMET_DPhi");
  tree_GenLev->Branch("MuTTrk_Mass",&MuTTrk_Mass,"MuTTrk_Mass");
  tree_GenLev->Branch("MuTTrk_Pt",&MuTTrk_Pt,"MuTTrk_Pt");
  tree_GenLev->Branch("MuTTrk_DR",&MuTTrk_DR,"MuTTrk_DR");
  tree_GenLev->Branch("MuTTrkMET_DPhi",&MuTTrkMET_DPhi,"MuTTrkMET_DPhi");
  tree_GenLev->Branch("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass,"MuLTrkMuTTrk_Mass");
  tree_GenLev->Branch("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt,"MuLTrkMuTTrk_Pt");
  tree_GenLev->Branch("MET_Pt",&MET_Pt,"MET_Pt");
  tree_GenLev->Branch("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass,"MuLTrkMuTTrkMET_Mass");
  tree_GenLev->Branch("Eventweight",&Eventweight,"Eventweight");

  ////// Higgs ////////// 
  float higgsPt;
  float higgsSMPt;
  bool isWplus;
  bool isWminus;
  bool isZ;
  float genWeight;

  TTree * higgsTree = new TTree("higgsTree","");
  higgsTree->Branch("HiggsPt",&higgsPt,"HiggsPt/F");
  higgsTree->Branch("isWplus",&isWplus,"isWplus/O");
  higgsTree->Branch("isWminus",&isWminus,"isWminus/O");
  higgsTree->Branch("isZ",&isZ,"isZ/O");
  higgsTree->Branch("genweight",&genWeight,"genweight/F");

  TTree * higgsSMTree = new TTree("higgsSMTree","");
  higgsSMTree->Branch("HiggsSMPt",&higgsSMPt,"HiggsSMPt/F");
  higgsSMTree->Branch("genweight",&genWeight,"genweight/F");
  
  ////// Others //////
  TH1D * triggerH = new TH1D("triggerH","",2,-0.5,1.5);
  TH1D * triggerSigH = new TH1D("triggerSigH","",2,-0.5,1.5);
  TH1D * nmuonsH = new TH1D("nmuonsH","",10,-0.5,9.5);
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// ************* DATA & MC Corrections ********************* //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Run-lumi selector
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
  std::vector<Period> periods;  
  if (isData) { // read the good runs 
    std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
    if (inputFileStream.fail() ) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
    }
    
    for(std::string s; std::getline(inputFileStream, s); ) {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Higgs reweighting
  TString fullpath_HiggsPtFile = TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+HiggsPtFileName;
  TFile * higgsPtFile = NULL;
  TH1D * higgsPtH = NULL;
  TH1D * higgsPt_WPlusH = NULL;
  TH1D * higgsPt_WMinusH = NULL;
  TH1D * higgsPt_ZH = NULL;
  if (applyHiggsPtWeight) { 
    std::cout << "ApplyHiggsPtWeight = " << applyHiggsPtWeight << std::endl;
    higgsPtFile = new TFile(fullpath_HiggsPtFile);
    if (higgsPtFile->IsZombie()) {
      std::cout << fullpath_HiggsPtFile << "  not found" << std::endl;
      exit(-1);
    }
    if (isVH) {
      std::cout << "IsVH = " << isVH << std::endl;
      higgsPt_WPlusH = (TH1D*)higgsPtFile->Get("kfactor_WplusH");
      higgsPt_WMinusH = (TH1D*)higgsPtFile->Get("kfactor_WminusH");
      higgsPt_ZH = (TH1D*)higgsPtFile->Get("kfactor_ZH");
      if (higgsPt_WPlusH==NULL) {
	std::cout << "histogram kfactor_WplusH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
      if (higgsPt_WMinusH==NULL) {
	std::cout << "histogram kfactor_WminusH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
      if (higgsPt_ZH==NULL) {
	std::cout << "histogram kfactor_ZH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
    }
    else {
      higgsPtH = (TH1D*)higgsPtFile->Get("kfactor");
      if (higgsPtH==NULL) {
	std::cout << "histogram kfactor is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);	
      }
    }
  }
  //  std::cout << "Higgs Pt histogram : " << higgsPtH << std::endl;

  // Trigger efficiencies
  ScaleFactor * SF_muon17 = new ScaleFactor();
  SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));
  ScaleFactor * SF_muon8 = new ScaleFactor();
  SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));

  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load( pathToCrystalLib );
  if (openSuccessful !=0 ) {
    cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit( -1 );
  }

  // Correction workspace
  TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/"+TString(CorrectionWSFileName);
  TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
  if (correctionWorkSpaceFile->IsZombie()) {
    std::cout << correctionsWorkspaceFileName << " does not exist " << std::endl;
  }
  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");
  if (correctionWS==NULL) {
    std::cout << "correction workspace not found " << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// ************* Processing NTuples Files ********************* //////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TString filen;
  int iFiles = 0;
  int events = 0;
  while (fileList >> filen) { // Starting Filelist loop
   iFiles++;
   cout << "file " << iFiles << " : " << filen << endl;
   
   TFile * file_ = TFile::Open(TString(filen));
   if (file_==NULL) continue;

   // sum of weights (needed for normalization)
   TTree * _inittree = (TTree*)file_->Get(TString(initNtupleName));
   if (_inittree!=NULL) {
     Float_t Genweight;
     if (!isData)
       _inittree->SetBranchAddress("genweight",&Genweight);
     Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
     std::cout << "Number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
     for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
       _inittree->GetEntry(iEntry);
       if (isData)
	 histWeightsH->Fill(0.,1.);
       else
	 histWeightsH->Fill(0.,Genweight);
     }
   }

   TTree * tree_ = (TTree*)file_->Get(TString(chainName));
   
   if (tree_==NULL) continue;

   tree_->SetMaxVirtualSize(3000000);
   // event info
   tree_->SetBranchAddress("event_nr", &event_nr);
   tree_->SetBranchAddress("event_run", &event_run);
   tree_->SetBranchAddress("event_luminosityblock", &event_luminosityblock);

   // Muons
   tree_->SetBranchAddress("muon_count", &muon_count);
   tree_->SetBranchAddress("muon_nMuonStations", muon_nMuonStations);
   tree_->SetBranchAddress("muon_nMuonHits", muon_nMuonHits);
   tree_->SetBranchAddress("muon_nPixelHits", muon_nPixelHits);
   tree_->SetBranchAddress("muon_nTrackerHits", muon_nTrackerHits);
   tree_->SetBranchAddress("muon_px", muon_px);
   tree_->SetBranchAddress("muon_py", muon_py);
   tree_->SetBranchAddress("muon_pz", muon_pz);
   tree_->SetBranchAddress("muon_pt", muon_pt);
   tree_->SetBranchAddress("muon_eta", muon_eta);
   tree_->SetBranchAddress("muon_phi", muon_phi);
   tree_->SetBranchAddress("muon_pterror", muon_pterror);
   tree_->SetBranchAddress("muon_chi2", muon_chi2);
   tree_->SetBranchAddress("muon_ndof", muon_ndof);
   tree_->SetBranchAddress("muon_charge", muon_charge);
   tree_->SetBranchAddress("muon_dxy", muon_dxy);
   tree_->SetBranchAddress("muon_dxyerr", muon_dxyerr);
   tree_->SetBranchAddress("muon_dz", muon_dz);
   tree_->SetBranchAddress("muon_dzerr", muon_dzerr);
   tree_->SetBranchAddress("muon_chargedHadIso", muon_chargedHadIso);
   tree_->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso);
   tree_->SetBranchAddress("muon_photonIso", muon_photonIso);
   tree_->SetBranchAddress("muon_puIso", muon_puIso);
   tree_->SetBranchAddress("muon_isMedium", muon_isMedium);
   tree_->SetBranchAddress("muon_isICHEP", muon_isICHEP);

   // MET
   tree_->SetBranchAddress("pfmetcorr_ex", &metx);
   tree_->SetBranchAddress("pfmetcorr_ey", &mety);
 

   // Tracks
   tree_->SetBranchAddress("track_count", &track_count);
   tree_->SetBranchAddress("track_ID", track_ID);
   tree_->SetBranchAddress("track_px", track_px);
   tree_->SetBranchAddress("track_py", track_py);
   tree_->SetBranchAddress("track_pz", track_pz);
   tree_->SetBranchAddress("track_pt", track_pt);
   tree_->SetBranchAddress("track_eta", track_eta);
   tree_->SetBranchAddress("track_phi", track_phi);
   tree_->SetBranchAddress("track_mass", track_mass);
   tree_->SetBranchAddress("track_charge", track_charge);
   tree_->SetBranchAddress("track_dxy", track_dxy);
   tree_->SetBranchAddress("track_dxyerr", track_dxyerr);
   tree_->SetBranchAddress("track_dz", track_dz);
   tree_->SetBranchAddress("track_dzerr",track_dzerr);
   tree_->SetBranchAddress("track_highPurity", track_highPurity);

   // trigger objects
   tree_->SetBranchAddress("trigobject_count", &trigobject_count);
   tree_->SetBranchAddress("trigobject_px", trigobject_px);
   tree_->SetBranchAddress("trigobject_py", trigobject_py);
   tree_->SetBranchAddress("trigobject_pz", trigobject_pz);
   tree_->SetBranchAddress("trigobject_pt", trigobject_pt);
   tree_->SetBranchAddress("trigobject_eta", trigobject_eta);
   tree_->SetBranchAddress("trigobject_phi", trigobject_phi);
   tree_->SetBranchAddress("trigobject_filters",trigobject_filters);

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   //   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);


   if (!isData) {
     tree_->SetBranchAddress("genweight",&genweight);
     tree_->SetBranchAddress("genparticles_count", &genparticles_count);
     tree_->SetBranchAddress("genparticles_e", genparticles_e);
     tree_->SetBranchAddress("genparticles_px", genparticles_px);
     tree_->SetBranchAddress("genparticles_py", genparticles_py);
     tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
     tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
     tree_->SetBranchAddress("genparticles_status", genparticles_status);
     tree_->SetBranchAddress("genparticles_info", genparticles_info);
     tree_->SetBranchAddress("genparticles_fromHardProcess", genparticles_fromHardProcess);
     tree_->SetBranchAddress("genparticles_status", genparticles_status);
     tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);
   }   

   // jets
   tree_->SetBranchAddress("pfjet_count", &pfjet_count);
   tree_->SetBranchAddress("pfjet_e", pfjet_e);
   tree_->SetBranchAddress("pfjet_px", pfjet_px);
   tree_->SetBranchAddress("pfjet_py", pfjet_py);
   tree_->SetBranchAddress("pfjet_pz", pfjet_pz);
   tree_->SetBranchAddress("pfjet_pt", pfjet_pt);
   tree_->SetBranchAddress("pfjet_eta", pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi", pfjet_phi);
   tree_->SetBranchAddress("pfjet_flavour", pfjet_flavour);
   tree_->SetBranchAddress("pfjet_btag", pfjet_btag);
   //   tree_->SetBranchAddress("pfjet_pu_jet_fullId_loose", pfjet_pu_jet_fullId_loose);
   //   tree_->SetBranchAddress("pfjet_pu_jet_fullId_medium", pfjet_pu_jet_fullId_medium);
   //   tree_->SetBranchAddress("pfjet_pu_jet_fullId_tight", pfjet_pu_jet_fullId_tight);
   
   // genjets
   if (!isData) {
     tree_->SetBranchAddress("genjets_count",&genjets_count);
     tree_->SetBranchAddress("genjets_e",genjets_e);
     tree_->SetBranchAddress("genjets_px",genjets_px);
     tree_->SetBranchAddress("genjets_py",genjets_py);
     tree_->SetBranchAddress("genjets_pz",genjets_pz);
     tree_->SetBranchAddress("genjets_pt",genjets_pt);
     tree_->SetBranchAddress("genjets_eta",genjets_eta);
     tree_->SetBranchAddress("genjets_phi",genjets_phi);
     tree_->SetBranchAddress("genjets_pdgid",genjets_pdgid);
     tree_->SetBranchAddress("genjets_status",genjets_status);
   }


  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
   TString sample(argv[2]);
   bool is4tau = sample.Contains("ToAA_AToTauTau");
  
   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events in Main Tree = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {  /////// Starting Event Loop ........ ////////
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     if (!isData) {
       weight *= genweight;
     }

     // ********************
     //   gen-level study
     // ********************
     
     genWeight = genweight;

     std::vector<unsigned int> posPion; posPion.clear();
     std::vector<unsigned int> negPion; negPion.clear();
     std::vector<unsigned int> posMuon; posMuon.clear();
     std::vector<unsigned int> negMuon; negMuon.clear();

     bool HiggsFound = false;
     unsigned int higgsIndex = 0;
     bool HiggsSMFound = false;
     unsigned int higgsSMIndex = 0;
     unsigned int nmuons = 0;

     isWplus = false;
     isWminus = false;
     isZ = false;

     if (!isData) {
       
       for (unsigned int iP=0; iP<genparticles_count; ++iP) 
       {
		   
	       if (genparticles_pdgid[iP]==13||genparticles_pdgid[iP]==-13) nmuons++;
	       if (genparticles_status[iP]==1&&genparticles_info[iP]==12) 
	       {
	           if (genparticles_pdgid[iP]==13) negMuon.push_back(iP);
	           if (genparticles_pdgid[iP]==-13) posMuon.push_back(iP);
	           if (genparticles_pdgid[iP]==211) posPion.push_back(iP);
	           if (genparticles_pdgid[iP]==-211) negPion.push_back(iP);
	       }
	       if (genparticles_pdgid[iP]==35) 
	       {
	           higgsIndex = iP;
	           HiggsFound = true;
	       }
	       if (genparticles_pdgid[iP]==25) 
	       {
	       higgsSMIndex = iP;
	       HiggsSMFound = true;
	       }
	       if (genparticles_pdgid[iP]==24) isWplus = true;
	       if (genparticles_pdgid[iP]==-24) isWminus = true;
	       if (genparticles_pdgid[iP]==23) isZ = true;
       }
       
       if (posMuon.size()==2&&negMuon.size()==0&&negPion.size()==2&&posPion.size()==0) 
       {
		   
		   unsigned int indexLeading = posMuon.at(0);
		   unsigned int indexTrailing = posMuon.at(1);
		   
		   
		   TLorentzVector LmuonLV; LmuonLV.SetXYZT(genparticles_px[indexLeading],
			                                       genparticles_py[indexLeading],
			                                       genparticles_pz[indexLeading],
			                                       genparticles_e[indexLeading]);
			                                       
		   TLorentzVector TmuonLV; TmuonLV.SetXYZT(genparticles_px[indexTrailing],
			                                       genparticles_py[indexTrailing],
			                                       genparticles_pz[indexTrailing],
			                                       genparticles_e[indexTrailing]);
		   
		   if(TmuonLV.Pt()>LmuonLV.Pt())
		   {
			  TLorentzVector muonLV=LmuonLV;
			  LmuonLV=TmuonLV;
			  TmuonLV=muonLV;
		   }
		   
		   float dRmuons = deltaR(LmuonLV.Eta(),LmuonLV.Phi(),
				                  TmuonLV.Eta(),TmuonLV.Phi());
				                  
		   if(dRmuons<1.5) continue;
		   
		   unsigned int LpionIndex = negPion.at(0);
		   unsigned int TpionIndex = negPion.at(1);
		   
		   TLorentzVector LpionLV; LpionLV.SetXYZT(genparticles_px[LpionIndex],
						                           genparticles_py[LpionIndex],
						                           genparticles_pz[LpionIndex],
						                           genparticles_e[LpionIndex]);
						                           
		   TLorentzVector TpionLV; TpionLV.SetXYZT(genparticles_px[TpionIndex],
						                           genparticles_py[TpionIndex],
						                           genparticles_pz[TpionIndex],
						                           genparticles_e[TpionIndex]);			                           
						                           
	      
	       float dRL = deltaR(LpionLV.Eta(),LpionLV.Phi(),LmuonLV.Eta(),LmuonLV.Phi());
	       float dRT = deltaR(TpionLV.Eta(),TpionLV.Phi(),LmuonLV.Eta(),LmuonLV.Phi());
	       
	       if(dRT<dRL)
	       {
			   TLorentzVector pionLV=LpionLV;
			   LpionLV=TpionLV;
			   TpionLV=pionLV;
		   }
				             
	       TLorentzVector DiMuonLV = LmuonLV + TmuonLV;
	       TLorentzVector Met4LV; Met4LV.SetXYZM(metx,mety,0,0);
	       TLorentzVector TrkLMuonLV = LmuonLV + LpionLV;
	       TLorentzVector TrkTMuonLV = TmuonLV + TpionLV;
	       TLorentzVector VisibleLV = TrkLMuonLV + TrkTMuonLV;
	       TLorentzVector VisibleandMETLV = VisibleLV + Met4LV;
	       
	       ////**** Setting the values of the variables *****/////
           MuLMuT_Mass  = DiMuonLV.M();
           MuLMuT_DR  = dRmuons;
           MuLMuT_DPhi = dPhiFrom2P(LmuonLV.Px(),LmuonLV.Py(),TmuonLV.Px(),TmuonLV.Py());
     
           MuLTrk_Mass = TrkLMuonLV.M();
           MuLTrk_Pt = TrkLMuonLV.Pt();
           MuLTrk_DR = deltaR(LmuonLV.Eta(),LmuonLV.Phi(),LpionLV.Eta(),LpionLV.Phi());
           MuLTrkMET_DPhi = dPhiFrom2P(TrkLMuonLV.Px(),TrkLMuonLV.Py(),Met4LV.Px(),Met4LV.Py());
      
           MuTTrk_Mass = TrkTMuonLV.M();
           MuTTrk_Pt = TrkTMuonLV.Pt();
           MuTTrk_DR = deltaR(TmuonLV.Eta(),TmuonLV.Phi(),TpionLV.Eta(),TpionLV.Phi());
           MuTTrkMET_DPhi = dPhiFrom2P(TrkTMuonLV.Px(),TrkTMuonLV.Py(),Met4LV.Px(),Met4LV.Py());
     
           MuLTrkMuTTrk_Mass = VisibleLV.M();
           MuLTrkMuTTrk_Pt = VisibleLV.Pt();
           MuLTrkMuTTrkMET_Mass = VisibleandMETLV.M();
           MET_Pt = Met4LV.Pt();
     
           Eventweight = 1.;
           
           tree_GenLev->Fill();
       }
       
       
       if (posMuon.size()==0&&negMuon.size()==2&&negPion.size()==0&&posPion.size()==2) 
       {
		   
		   unsigned int indexLeading = negMuon.at(0);
		   unsigned int indexTrailing = negMuon.at(1);
		   
		   
		   TLorentzVector LmuonLV; LmuonLV.SetXYZT(genparticles_px[indexLeading],
			                                       genparticles_py[indexLeading],
			                                       genparticles_pz[indexLeading],
			                                       genparticles_e[indexLeading]);
			                                       
		   TLorentzVector TmuonLV; TmuonLV.SetXYZT(genparticles_px[indexTrailing],
			                                       genparticles_py[indexTrailing],
			                                       genparticles_pz[indexTrailing],
			                                       genparticles_e[indexTrailing]);
		   
		   if(TmuonLV.Pt()>LmuonLV.Pt())
		   {
			  TLorentzVector muonLV=LmuonLV;
			  LmuonLV=TmuonLV;
			  TmuonLV=muonLV;
		   }
		   
		   float dRmuons = deltaR(LmuonLV.Eta(),LmuonLV.Phi(),
				                  TmuonLV.Eta(),TmuonLV.Phi());
				                  
		   if(dRmuons<1.5) continue;
		   
		   unsigned int LpionIndex = posPion.at(0);
		   unsigned int TpionIndex = posPion.at(1);
		   
		   TLorentzVector LpionLV; LpionLV.SetXYZT(genparticles_px[LpionIndex],
						                           genparticles_py[LpionIndex],
						                           genparticles_pz[LpionIndex],
						                           genparticles_e[LpionIndex]);
						                           
		   TLorentzVector TpionLV; TpionLV.SetXYZT(genparticles_px[TpionIndex],
						                           genparticles_py[TpionIndex],
						                           genparticles_pz[TpionIndex],
						                           genparticles_e[TpionIndex]);			                           
						                           
	      
	       float dRL = deltaR(LpionLV.Eta(),LpionLV.Phi(),LmuonLV.Eta(),LmuonLV.Phi());
	       float dRT = deltaR(TpionLV.Eta(),TpionLV.Phi(),LmuonLV.Eta(),LmuonLV.Phi());
	       
	       if(dRT<dRL)
	       {
			   TLorentzVector pionLV=LpionLV;
			   LpionLV=TpionLV;
			   TpionLV=pionLV;
		   }
				             
	       TLorentzVector DiMuonLV = LmuonLV + TmuonLV;
	       TLorentzVector Met4LV; Met4LV.SetXYZM(metx,mety,0,0);
	       TLorentzVector TrkLMuonLV = LmuonLV + LpionLV;
	       TLorentzVector TrkTMuonLV = TmuonLV + TpionLV;
	       TLorentzVector VisibleLV = TrkLMuonLV + TrkTMuonLV;
	       TLorentzVector VisibleandMETLV = VisibleLV + Met4LV;
	       
	       ////**** Setting the values of the variables *****/////
           MuLMuT_Mass  = DiMuonLV.M();
           MuLMuT_DR  = dRmuons;
           MuLMuT_DPhi = dPhiFrom2P(LmuonLV.Px(),LmuonLV.Py(),TmuonLV.Px(),TmuonLV.Py());
     
           MuLTrk_Mass = TrkLMuonLV.M();
           MuLTrk_Pt = TrkLMuonLV.Pt();
           MuLTrk_DR = deltaR(LmuonLV.Eta(),LmuonLV.Phi(),LpionLV.Eta(),LpionLV.Phi());
           MuLTrkMET_DPhi = dPhiFrom2P(TrkLMuonLV.Px(),TrkLMuonLV.Py(),Met4LV.Px(),Met4LV.Py());
      
           MuTTrk_Mass = TrkTMuonLV.M();
           MuTTrk_Pt = TrkTMuonLV.Pt();
           MuTTrk_DR = deltaR(TmuonLV.Eta(),TmuonLV.Phi(),TpionLV.Eta(),TpionLV.Phi());
           MuTTrkMET_DPhi = dPhiFrom2P(TrkTMuonLV.Px(),TrkTMuonLV.Py(),Met4LV.Px(),Met4LV.Py());
     
           MuLTrkMuTTrk_Mass = VisibleLV.M();
           MuLTrkMuTTrk_Pt = VisibleLV.Pt();
           MuLTrkMuTTrkMET_Mass = VisibleandMETLV.M();
           MET_Pt = Met4LV.Pt();
     
           Eventweight = 1.;
           
           tree_GenLev->Fill();
       }
       
       
     }
     
     
     if (HiggsFound) 
     {
       TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[higgsIndex],
					       genparticles_py[higgsIndex],
					       genparticles_pz[higgsIndex],
					       genparticles_e[higgsIndex]);
       higgsPt = higgsLV.Pt();
       higgsTree->Fill();
       
       if (applyHiggsPtWeight&&is4tau) 
       {
	   double HiggsPtForWeighting = higgsPt;
	   if (higgsPt>500) HiggsPtForWeighting = 499;
	   double higgsPtWeight = 1;
	   if (isVH) {
	     if (isWplus)
	       higgsPtWeight = higgsPt_WPlusH->GetBinContent(higgsPt_WPlusH->FindBin(HiggsPtForWeighting));
	     if (isWminus)
	       higgsPtWeight = higgsPt_WMinusH->GetBinContent(higgsPt_WMinusH->FindBin(HiggsPtForWeighting));
	     if (isZ)
	       higgsPtWeight = higgsPt_ZH->GetBinContent(higgsPt_ZH->FindBin(HiggsPtForWeighting));
	   }
	   else {
	     higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
	   }
	   weight *= higgsPtWeight;
	   HiggsPtWeightH->Fill(higgsPtWeight);
	   //	   std::cout << "HiggsPt weight (Pythia) = " << higgsPtWeight << std::endl;
       }
       
     }

     if (HiggsSMFound) 
     {
       TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[higgsSMIndex],
                                               genparticles_py[higgsSMIndex],
                                               genparticles_pz[higgsSMIndex],
                                               genparticles_e[higgsSMIndex]);
       higgsSMPt = higgsLV.Pt();

       if (applyHiggsPtWeight&&!is4tau)
	 {
           double HiggsPtForWeighting = higgsSMPt;
           if (higgsSMPt>500) HiggsPtForWeighting = 499;
           double higgsPtWeight = 1;
	   higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
	   weight *= higgsPtWeight;
	   HiggsPtWeightH->Fill(higgsPtWeight);
	   //	   std::cout << "HiggsPt weight (Madgraph) = " << higgsPtWeight << std::endl;
	 }

       higgsSMTree->Fill();
     }
     
     counter_InputEventsH->Fill(1.0,weight);


     // ***********************
     // Data taking corrections
     // ***********************
     
     if (isData) {
	  if (applyGoodRunSelection) {
	  bool lumi = false;
	  int n=event_run;
	  int lum = event_luminosityblock;
	  
	  std::string num = std::to_string(n);
	  std::string lnum = std::to_string(lum);
	  for(const auto& a : periods)
	    {
	      if ( num.c_str() ==  a.name ) {
		//	      std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
		//std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
		
		for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		  
		  //   cout<<b->lower<<"  "<<b->bigger<<endl;
		  if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
		}
		auto last = std::prev(a.ranges.end());
		// std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
		if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	      }
	    }
	  if (!lumi) continue;
	  }
     }

     float puweight = 1;
     if (!isData) {
       puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
       //       std::cout << "n(true interactions) = " << numtruepileupinteractions << "   :  PU weight = " << puweight << std::endl; 
     }
     puWeightH->Fill(puweight,1.0);
     weight *= puweight;

     nmuonsH->Fill(double(nmuons));

     
     // ***********************
     //   HLT Filter library
     // ***********************
    
     unsigned int nMu8Leg   = 0;
     unsigned int nMu17Leg  = 0;
     unsigned int nDZFilter = 0;
     unsigned int nSSFilter = 0;
     bool isMu8Leg = false;
     bool isMu17Leg = false;
     bool isDZFilter = false;
     bool isSSFilter = false;

     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) 
     {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       if (HLTFilter==MuonHighPtFilterName) 
       {
	       nMu17Leg = i;
	       isMu17Leg = true;
       }
       if (HLTFilter==MuonLowPtFilterName1||HLTFilter==MuonLowPtFilterName2) 
       {
	       nMu8Leg = i;
	       isMu8Leg = true;
       }
       if (HLTFilter==DiMuonDzFilterName) 
       {
	       nDZFilter = i;
	       isDZFilter = true;
       }
       if (HLTFilter==DiMuonSameSignFilterName) 
       {
	       nSSFilter = i;
	       isSSFilter = true;
       }
     }
     
     if (!isMu17Leg) 
     {
       cout << "Filter " << MuonHighPtFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isMu8Leg) 
     {
       cout << "Filters " << MuonLowPtFilterName1 
	    << " or " << MuonLowPtFilterName2
	    << " not found " << endl;
       exit(-1);
     }
     if (!isDZFilter) 
     {
       cout << "Filter " << DiMuonDzFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isSSFilter) 
     {
       cout << "Filter " << DiMuonSameSignFilterName << " not found " << endl;
       exit(-1);
     }


     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();

     for(UInt_t i=0;i<muon_count;i++)
     {
		 
         bool muonID = muon_isMedium[i]; // MC          
	 
         if (isData) 
         {
	         if (event_run >= 278820 && muon_isMedium[i]) muonID = true; // Run2016G-H
	         if (event_run < 278820 && muon_isICHEP[i]) muonID = true; // Run2016B-F
		   
         }

         if (!muonID) continue;
         if(fabs(muon_dxy[i])>dxyMuonCut) continue;
         if(fabs(muon_dz[i])>dzMuonCut) continue;
         if(muon_pt[i]<ptMuonLowCut) continue;
         if(fabs(muon_eta[i])>etaMuonLowCut) continue;

         rand.SetSeed((int)((muon_eta[i]+2.41)*100000));
         double rannum  = rand.Rndm();
         
         muons.push_back(i);
     } // end loop over good muons

     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<2) continue; // quit event if number of good muons < 2
     
     counter_MuonSizeGTE2H->Fill(1.,weight);

     bool vetoOSdimuons = false;
     if (applyOSdimuonVeto) {
       for (unsigned int im=0; im<muons.size()-1; ++im) {
	 for (unsigned int jm=im+1; jm<muons.size(); ++jm) {
	   unsigned int index1 = muons[im];
	   unsigned int index2 = muons[jm];
	   float charge = muon_charge[index1] *  muon_charge[index2];
	   if (charge>0.0) continue;
	   TLorentzVector muon1; muon1.SetXYZM(muon_px[index1],
					       muon_py[index1],
					       muon_pz[index1],
					       MuMass); 
	   TLorentzVector muon2; muon2.SetXYZM(muon_px[index2],
					       muon_py[index2],
					       muon_pz[index2],
					       MuMass); 
	   TLorentzVector dimuons_os = muon1+muon2;
	   float ptSum = dimuons_os.Pt();
	   if (ptSum>maxPtSumOSdimuons) vetoOSdimuons = true;

	 }
       }
     }

     if (vetoOSdimuons) continue;
     
     // **********************
     // selecting good tracks
     // **********************
     vector<unsigned int> tracks; tracks.clear();

     for (unsigned int iTrk=0; iTrk<track_count; ++iTrk)
      {
          if (fabs(track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
          if (fabs(track_dxy[iTrk])>dxyTrkLooseCut) continue;
          if (fabs(track_dz[iTrk])>dzTrkLooseCut) continue;
          if (fabs(track_eta[iTrk])>etaTrkCut) continue;
          if (fabs(track_pt[iTrk])<ptTrkLooseCut) continue;


          tracks.push_back(iTrk);

      } // end loop over good tracks

     // *************************
     // selection of dimuon pair
     // *************************
     
     float maxPtSum = -1;
     int iLeading = -1;
     int iTrailing = -1;

     for (unsigned int i1=0; i1<muons.size()-1; ++i1) // first loop over muons
     {
          int index1 = muons.at(i1);
           
          for (unsigned int i2=i1+1; i2<muons.size(); ++i2) // second loop over muons
          {
	           int index2 = muons.at(i2);
	           float ptSum = muon_pt[index1] + muon_pt[index2];
	           float charge = muon_charge[index1] *  muon_charge[index2];
	           bool chargeSelection = charge<0;
	           if (sameSign) chargeSelection = charge>0;
	           if (!chargeSelection) continue;
	           
	           float dRmuons = deltaR(muon_eta[index1],muon_phi[index1],
				                      muon_eta[index2],muon_phi[index2]);

	           if (dRmuons<dRMuonsCut) continue;
	           
	           bool mu1MatchMu17 = false;
	           bool mu1MatchMu8  = false;
	           bool mu1MatchDz   = false;
	           bool mu1MatchSS   = false;
	 
	           for (unsigned int iT=0; iT<trigobject_count; ++iT) 
	           {
	                float dRtrig = deltaR(muon_eta[index1],muon_phi[index1],
				                          trigobject_eta[iT],trigobject_phi[iT]);
	   
	                if (dRtrig>DRTrigMatch) continue;
	                if (trigobject_filters[iT][nMu17Leg]) mu1MatchMu17 = true;
	                if (trigobject_filters[iT][nMu8Leg]) mu1MatchMu8 = true;
	                if (trigobject_filters[iT][nDZFilter]) mu1MatchDz = true;
	                if (trigobject_filters[iT][nSSFilter]) mu1MatchSS = true;
	           }
	           
	           bool mu2MatchMu17 = false;
	           bool mu2MatchMu8  = false;
	           bool mu2MatchDz   = false;
	           bool mu2MatchSS   = false;
	 
	           for (unsigned int iT=0; iT<trigobject_count; ++iT) 
	           {
	                float dRtrig = deltaR(muon_eta[index2],muon_phi[index2],
				                          trigobject_eta[iT],trigobject_phi[iT]);
	   
	                if (dRtrig>DRTrigMatch) continue;
	                if (trigobject_filters[iT][nMu17Leg]) mu2MatchMu17 = true;
	                if (trigobject_filters[iT][nMu8Leg]) mu2MatchMu8 = true;
	                if (trigobject_filters[iT][nDZFilter]) mu2MatchDz = true;
	                if (trigobject_filters[iT][nSSFilter]) mu2MatchSS = true;
	           }

	           bool mu1PtHigh = muon_pt[index1]>ptMuonHighCut && fabs(muon_eta[index1])<etaMuonHighCut;
	           bool mu1PtLow  = muon_pt[index1]>ptMuonLowCut && fabs(muon_eta[index1])<etaMuonLowCut;
	           bool mu2PtHigh = muon_pt[index2]>ptMuonHighCut && fabs(muon_eta[index2])<etaMuonHighCut;
	           bool mu2PtLow  = muon_pt[index2]>ptMuonLowCut && fabs(muon_eta[index2])<etaMuonLowCut;

	           // trigger condition
	           bool isTriggerMatched = true;
	           if (applyTriggerMatch) 
	           {
	               isTriggerMatched = (mu1MatchMu17&&mu2MatchMu8&&mu1PtHigh&&mu2PtLow)||(mu1MatchMu8&&mu2MatchMu17&&mu1PtLow&&mu2PtHigh);
	               if (isData) 
	               {
	                   isTriggerMatched = isTriggerMatched && mu1MatchSS && mu2MatchSS;
	                   
	                   if (event_run<=274442||(event_run>=280919&&event_run<=284044)) isTriggerMatched = isTriggerMatched && mu1MatchDz && mu2MatchDz; // when dZ filter is present
	               }
	           }
	           else 
	           {
	               isTriggerMatched = (mu1PtHigh&&mu2PtLow) || (mu1PtLow&&mu2PtHigh);
	           }
	           
	           if (!isTriggerMatched) continue;
	           
	           if (ptSum>maxPtSum) // choose the mair with maximum Sum(pT)
	           { 
	               maxPtSum = ptSum;
	               
	               if (muon_pt[index1]>muon_pt[index2]) 
	               {
	                   iLeading = index1;
	                   iTrailing = index2;
	               }
	               else 
	               {
	                   iLeading = index2;
	                   iTrailing = index1;
	               }
	           }
	           
          } // end of first loop over muons
          
     } // end of second loop over muons


     if (iLeading<0) continue;
     if (iTrailing<0) continue;

     // *************************************************
     //   Trigger efficiency and muon tracking efficiencies
     // *************************************************
     
     double triggerWeight = 1;
     double effTrkLeading = 1;
     double effTrkTrailing = 1;

     if (!isData)
     {   
         double ptLeading = muon_pt[iLeading];
         double etaLeading = muon_eta[iLeading];
         double ptTrailing = muon_pt[iTrailing];
         double etaTrailing = muon_eta[iTrailing];

         double effMu17dataTrailing = SF_muon17->get_EfficiencyData(ptTrailing,etaTrailing);
         double effMu8dataTrailing = SF_muon8->get_EfficiencyData(ptTrailing,etaTrailing);
         double effMu17dataLeading = SF_muon17->get_EfficiencyData(ptLeading,etaLeading);
         double effMu8dataLeading = SF_muon8->get_EfficiencyData(ptLeading,etaLeading);

         double effMu17MCTrailing = SF_muon17->get_EfficiencyMC(ptTrailing,etaTrailing);
         double effMu8MCTrailing = SF_muon8->get_EfficiencyMC(ptTrailing,etaTrailing);
         double effMu17MCLeading = SF_muon17->get_EfficiencyMC(ptLeading,etaLeading);
         double effMu8MCLeading = SF_muon8->get_EfficiencyMC(ptLeading,etaLeading);

         double trigWeightData = effMu17dataLeading*effMu8dataTrailing + effMu17dataTrailing*effMu8dataLeading - effMu17dataLeading*effMu17dataTrailing;
         double trigWeightMC = effMu17MCLeading*effMu8MCTrailing + effMu17MCTrailing*effMu8MCLeading - effMu17MCLeading*effMu17MCTrailing;

         if (applyTriggerMatch) 
         {
	 
	         if (trigWeightMC>0) triggerWeight = trigWeightData/trigWeightMC;
         }
         else
         {
	         triggerWeight = trigWeightData;
	     }
	     
         triggerWeight *= effDzSS;

	 correctionWS->var("m_pt")->setVal(ptLeading);
	 correctionWS->var("m_eta")->setVal(etaLeading);
	 effTrkLeading = correctionWS->function("m_trk_ratio")->getVal();
	 correctionWS->var("m_pt")->setVal(ptTrailing);
	 correctionWS->var("m_eta")->setVal(etaTrailing);
	 effTrkTrailing = correctionWS->function("m_trk_ratio")->getVal();

     }
     
     triggerWeightH->Fill(triggerWeight,1.0);
     weight = weight*triggerWeight*effTrkLeading*effTrkTrailing;
     MuTrkWeightH->Fill(effTrkLeading*effTrkTrailing,1.0);

     // *****************************
     // Histos after dimuon selection
     // *****************************
     
     TLorentzVector LeadingMuon4; LeadingMuon4.SetXYZM(muon_px[iLeading],
						       muon_py[iLeading],
						       muon_pz[iLeading],
						       MuMass);
     
     TLorentzVector TrailingMuon4; TrailingMuon4.SetXYZM(muon_px[iTrailing],
							 muon_py[iTrailing],
							 muon_pz[iTrailing],
							 MuMass);
     TLorentzVector diMuon4 = LeadingMuon4 + TrailingMuon4;
     
     TLorentzVector Met4; Met4.SetXYZM(metx,mety,0,0);
     
     float dimuonMass = diMuon4.M();
     float dimuonDR  = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),TrailingMuon4.Eta(),TrailingMuon4.Phi());
     float dimuonDPhi  = dPhiFrom2P(LeadingMuon4.Px(),LeadingMuon4.Py(),TrailingMuon4.Px(),TrailingMuon4.Py());
     float metPt = Met4.Pt();


     // filling histograms (muon kinematics)
     dimuonMassH->Fill(dimuonMass,weight);
     ptLeadingMuH->Fill(muon_pt[iLeading],weight);
     ptTrailingMuH->Fill(muon_pt[iTrailing],weight);
     etaLeadingMuH->Fill(muon_eta[iLeading],weight);
     etaTrailingMuH->Fill(muon_eta[iTrailing],weight);
     counter_MuonKinematicsH->Fill(1.0,weight);                  
 
     /////////////////////////////////////////////////////////////////////////
     //////// ****** Selecting a Track around each muon ********* ////////////
     /////////////////////////////////////////////////////////////////////////
     
     ///////// Track around leading Muon //////////
     float maxPtSumMuLTrk = ptSumCut;
     int iTrkLeading = -1;

     for (unsigned int iTrk=0; iTrk<tracks.size(); ++iTrk)
     {
          if (fabs(track_dxy[iTrk])>dxyTrkCut) continue;
          if (fabs(track_dz[iTrk])>dzTrkCut) continue;
          if (track_pt[iTrk]<ptTrkCut) continue;
       
          int index = tracks.at(iTrk);

          TLorentzVector trk4; trk4.SetXYZM(track_px[index],
					    track_py[index],
					    track_pz[index],
					    track_mass[index]);

          TLorentzVector leadingMuDiff = LeadingMuon4 - trk4;
          TLorentzVector trailingMuDiff = TrailingMuon4 - trk4;

          TLorentzVector MuLTrack4 = LeadingMuon4 + trk4;

          if (leadingMuDiff.P()<0.1 || trailingMuDiff.P()<0.1) continue; // track is not any of the muons

	      float ptSum = MuLTrack4.Pt();
	      float drTrkMu = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				                 trk4.Eta(),trk4.Phi());
	      float qTrkLeadingMu = track_charge[index]*muon_charge[iLeading];
	      bool chargeSelection = qTrkLeadingMu<0;

	      if (drTrkMu<dRMuonsCut && chargeSelection &&  ptSum>maxPtSumMuLTrk)
		{
	          iTrkLeading=index;
		  maxPtSumMuLTrk=ptSum;
		}


     }

     if (iTrkLeading<0) continue;

     TLorentzVector TrkLeadingMuon4; TrkLeadingMuon4.SetXYZM(track_px[iTrkLeading],
							     track_py[iTrkLeading],
							     track_pz[iTrkLeading],
							     track_mass[iTrkLeading]);

     ///////// Track around trailing Muon //////////
     float maxPtSumMuTTrk = ptSumCut;
     int iTrkTrailing = -1;

     for (unsigned int iTrk=0; iTrk<tracks.size(); ++iTrk)
     { 
          if (fabs(track_dxy[iTrk])>dxyTrkCut) continue;
          if (fabs(track_dz[iTrk])>dzTrkCut) continue;
          if (track_pt[iTrk]<ptTrkCut) continue;
       
          int index = tracks.at(iTrk);

          TLorentzVector trk4; trk4.SetXYZM(track_px[index],
					    track_py[index],
					    track_pz[index],
					    track_mass[index]);

          TLorentzVector leadingMuDiff = LeadingMuon4 - trk4;
          TLorentzVector trailingMuDiff = TrailingMuon4 - trk4;
          TLorentzVector trkleadingMuDiff = TrkLeadingMuon4 - trk4;

          TLorentzVector MuTTrack4 = TrailingMuon4 + trk4;

          if (leadingMuDiff.P()<0.1 || trailingMuDiff.P()<0.1 || trkleadingMuDiff.P()<0.1) continue; // track is not any of the muons nor the previous track

	      float ptSum = MuTTrack4.Pt();
	      float drTrkMu = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				                 trk4.Eta(),trk4.Phi());
	      float qTrkTrailingMu = track_charge[index]*muon_charge[iTrailing];
	      bool chargeSelection = qTrkTrailingMu<0;

	      if (drTrkMu<dRMuonsCut && chargeSelection && ptSum>maxPtSumMuTTrk)
		{
	          iTrkTrailing=index;
		  maxPtSumMuTTrk=ptSum;
		}


     }

     if (iTrkTrailing<0) continue;
     
     
     
     // ************************************
     // Kinematics after 4 object selection
     // ************************************

     TLorentzVector TrkTrailingMuon4; TrkTrailingMuon4.SetXYZM(track_px[iTrkTrailing],
						                                       track_py[iTrkTrailing],
						                                       track_pz[iTrkTrailing],
						                                       track_mass[iTrkTrailing]);
						                         
						                         
	 TLorentzVector DiTauLeading4 = LeadingMuon4 + TrkLeadingMuon4;
	 TLorentzVector DiTauTrailing4 = TrailingMuon4 + TrkTrailingMuon4;
	
	 TLorentzVector Visible4 = DiTauLeading4 + DiTauTrailing4;
	 TLorentzVector VisibleandMET4 = Visible4 + Met4;
	
	
	 float muonLtrkMass = DiTauLeading4.M();
 	 float muonLtrkPt = DiTauLeading4.Pt();
	 float muonLtrkDR = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),TrkLeadingMuon4.Eta(),TrkLeadingMuon4.Phi());
	 float muonLtrkmetDPhi = dPhiFrom2P(DiTauLeading4.Px(),DiTauLeading4.Py(),Met4.Px(),Met4.Py());
 	 
	 float muonTtrkMass = DiTauTrailing4.M();
	 float muonTtrkPt = DiTauTrailing4.Pt();
	 float muonTtrkDR = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),TrkTrailingMuon4.Eta(),TrkTrailingMuon4.Phi());
	 float muonTtrkmetDPhi = dPhiFrom2P(DiTauTrailing4.Px(),DiTauTrailing4.Py(),Met4.Px(),Met4.Py());
 	 
	 float visibleMass = Visible4.M();
	 float visiblePt = Visible4.Pt();
	 float visibleandMetMass = VisibleandMET4.M();
	 
     
     
     // ************************************
     // Counting Tracks around each object
     // ************************************
        
     std::vector<unsigned int> trksLeadingMuon; trksLeadingMuon.clear(); // all tracks
     std::vector<unsigned int> trksTrailingMuon; trksTrailingMuon.clear(); // all tracks
     std::vector<unsigned int> trksTrkLeadingMuon; trksTrkLeadingMuon.clear(); // all tracks
     std::vector<unsigned int> trksTrkTrailingMuon; trksTrkTrailingMuon.clear(); // all tracks
     
     std::vector<unsigned int> softtrksLeadingMuon; softtrksLeadingMuon.clear(); // all tracks
     std::vector<unsigned int> softtrksTrailingMuon; softtrksTrailingMuon.clear(); // all tracks
     std::vector<unsigned int> softtrksTrkLeadingMuon; softtrksTrkLeadingMuon.clear(); // all tracks
     std::vector<unsigned int> softtrksTrkTrailingMuon; softtrksTrkTrailingMuon.clear(); // all tracks

     for (unsigned int iTrk=0; iTrk<track_count; ++iTrk) // loop over tracks
     {
          if (fabs(track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
          if (fabs(track_dxy[iTrk])>dxyTrkLooseCut) continue;
          if (fabs(track_dz[iTrk])>dzTrkLooseCut) continue;
          if (fabs(track_eta[iTrk])>etaTrkCut) continue;
          if (fabs(track_pt[iTrk])<ptTrkLooseCut) continue;

          TLorentzVector trk4; trk4.SetXYZM(track_px[iTrk],
					                        track_py[iTrk],
					                        track_pz[iTrk],
					                        track_mass[iTrk]);

         TLorentzVector leadingMuDiff = LeadingMuon4 - trk4;
         TLorentzVector trailingMuDiff = TrailingMuon4 - trk4;
         TLorentzVector trkleadingMuDiff = TrkLeadingMuon4 - trk4;
         TLorentzVector trktrailingMuDiff = TrkTrailingMuon4 - trk4;

         if (leadingMuDiff.P()<0.1 || trailingMuDiff.P()<0.1 || trkleadingMuDiff.P()<0.1 || trktrailingMuDiff.P()<0.1) continue; // track is not any of the selected objects
         {

	   float drtrkMuL = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				   track_eta[iTrk],track_phi[iTrk]);
	   float drtrkMuT = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				   track_eta[iTrk],track_phi[iTrk]);
	   float drtrkTrkMuL = deltaR(TrkLeadingMuon4.Eta(),TrkLeadingMuon4.Phi(),
				      track_eta[iTrk],track_phi[iTrk]);
	   float drtrkTrkMuT = deltaR(TrkTrailingMuon4.Eta(),TrkTrailingMuon4.Phi(),
				      track_eta[iTrk],track_phi[iTrk]);
	   
	   
	   if (drtrkMuL < dRIso)
	     { 
	       trksLeadingMuon.push_back(iTrk);
	       
	       if (fabs(track_pt[iTrk])<ptTrkCut) softtrksLeadingMuon.push_back(iTrk);
	     }   
	   if (drtrkMuT < dRIso)
	     {
	       trksTrailingMuon.push_back(iTrk);
	       
	       if (fabs(track_pt[iTrk])<ptTrkCut) softtrksTrailingMuon.push_back(iTrk);
	     }
	   if (drtrkTrkMuL < dRIso) 
	     { 
	       trksTrkLeadingMuon.push_back(iTrk);
	       
	       if (fabs(track_pt[iTrk])<ptTrkCut) softtrksTrkLeadingMuon.push_back(iTrk);
	     }
	   if (drtrkTrkMuT < dRIso) 
	     {
	       trksTrkTrailingMuon.push_back(iTrk);
	       
	       if (fabs(track_pt[iTrk])<ptTrkCut) softtrksTrkTrailingMuon.push_back(iTrk);
	     } 
	   
	   
         }

     } // end loop over tracks

     nTracksLeadingMuH->Fill(float(trksLeadingMuon.size()),weight);
     nTracksTrailingMuH->Fill(float(trksTrailingMuon.size()),weight);
     nTracksTrkLeadingMuH->Fill(float(trksTrkLeadingMuon.size()),weight);
     nTracksTrkTrailingMuH->Fill(float(trksTrkTrailingMuon.size()),weight);
     
     nSoftTracksLeadingMuH->Fill(float(softtrksLeadingMuon.size()),weight);
     nSoftTracksTrailingMuH->Fill(float(softtrksTrailingMuon.size()),weight);
     nSoftTracksTrkLeadingMuH->Fill(float(softtrksTrkLeadingMuon.size()),weight);
     nSoftTracksTrkTrailingMuH->Fill(float(softtrksTrkTrailingMuon.size()),weight);
     
     
     
     //////////////////////////////////////// ************************************** ////////////////////////////////////////////
     ////////////////////////////////////////  Defining Sidebands and Signal Region  ////////////////////////////////////////////
     //////////////////////////////////////// ************************************** ////////////////////////////////////////////


     // Signal Regions
     bool signalRegion = trksLeadingMuon.size()==0&&trksTrkLeadingMuon.size()==0&&trksTrailingMuon.size()==0&&trksTrkTrailingMuon.size()==0;


     // Control Regions for Background Modeling
     bool controlNoSR = !signalRegion;
     
     bool controlSemiIso = ( trksLeadingMuon.size()==0&&trksTrkLeadingMuon.size()==0&&trksTrailingMuon.size()!=0&&trksTrkTrailingMuon.size()!=0 )
                           ||
                           ( trksLeadingMuon.size()!=0&&trksTrkLeadingMuon.size()!=0&&trksTrailingMuon.size()==0&&trksTrkTrailingMuon.size()==0 );
                           

     bool controlLooseIso = ( trksLeadingMuon.size()!=0&&softtrksLeadingMuon.size()==trksLeadingMuon.size() )
                            &&
                            ( trksTrkLeadingMuon.size()!=0&&softtrksTrkLeadingMuon.size()==trksTrkLeadingMuon.size() )
                            &&
                            ( trksTrailingMuon.size()!=0&&softtrksTrailingMuon.size()==trksTrailingMuon.size() )
                            &&
                            ( trksTrkTrailingMuon.size()!=0&&softtrksTrkTrailingMuon.size()==trksTrkTrailingMuon.size() );
                            
                            
     bool controlLooseSemiIso = ( ( trksLeadingMuon.size()==0&&trksTrkLeadingMuon.size()==0 )
                                  &&
                                  ( trksTrailingMuon.size()!=0&&softtrksTrailingMuon.size()==trksTrailingMuon.size() )
                                  &&
                                  ( trksTrkTrailingMuon.size()!=0&&softtrksTrkTrailingMuon.size()==trksTrkTrailingMuon.size() )
                                )
                                ||
                                ( ( trksTrailingMuon.size()==0&&trksTrkTrailingMuon.size()==0 )
                                  &&
                                  ( trksLeadingMuon.size()!=0&&softtrksLeadingMuon.size()==trksLeadingMuon.size() )
                                  &&
                                  ( trksTrkLeadingMuon.size()!=0&&softtrksTrkLeadingMuon.size()==trksTrkLeadingMuon.size() )
                                );
     
                               
     bool controlLeadingSemiIso = ( trksLeadingMuon.size()==0&&trksTrkLeadingMuon.size()==0&&trksTrailingMuon.size()!=0&&trksTrkTrailingMuon.size()!=0 );
                          
                              
     bool controlLeadingLooseIso = ( trksLeadingMuon.size()==0&&trksTrkLeadingMuon.size()==0
                                   &&
                                   ( trksTrailingMuon.size()!=0&&softtrksTrailingMuon.size()==trksTrailingMuon.size() )
                                   &&
                                   ( trksTrkTrailingMuon.size()!=0&&softtrksTrkTrailingMuon.size()==trksTrkTrailingMuon.size() )
                                   );
         
         
     
     /////////////////////////****************///////////////////////////////
     //////******** Filling histograms by Regions & Categories ********//////
     /////////////////////////***************///////////////////////////////
     
     ////**** Setting the values of the variables *****/////
     MuLMuT_Mass  = dimuonMass;
     MuLMuT_DR  = dimuonDR;
     MuLMuT_DPhi = dimuonDPhi;
     
     MuLTrk_Mass = muonLtrkMass;
     MuLTrk_Pt = muonLtrkPt;
     MuLTrk_DR = muonLtrkDR;
     MuLTrkMET_DPhi = muonLtrkmetDPhi;
     
     MuTTrk_Mass = muonTtrkMass;
     MuTTrk_Pt = muonTtrkPt;
     MuTTrk_DR = muonTtrkDR;
     MuTTrkMET_DPhi = muonTtrkmetDPhi;
     
     MuLTrkMuTTrk_Mass = visibleMass;
     MuLTrkMuTTrk_Pt = visiblePt;
     MuLTrkMuTTrkMET_Mass = visibleandMetMass;
     MET_Pt = metPt;
     
     Eventweight = weight;
     
     
     ////**** Filling Tress *****/////
     if(controlNoSR) tree_NoSR->Fill();
	 
	        
     if(controlSemiIso) tree_SemiIso->Fill();
     

     if(controlLooseIso) tree_LooseIso->Fill();
     
     
     if(controlLooseSemiIso) tree_LooseSemiIso->Fill();
    
     
     if(controlLeadingSemiIso) tree_LeadingSemiIso->Fill();
     
     
     if(controlLeadingLooseIso) tree_LeadingLooseIso->Fill();
     
        
     
     if (signalRegion)
     {
       //event info
       //std::cout << event_run <<":"<< event_luminosityblock <<":"<< event_nr << std::endl;
       
       //////////////////////////////////////////////////////////////////////////
       /////////////////////// ** Signal Corrections ** /////////////////////////
       //////////////////////////////////////////////////////////////////////////
       if (!isData)
       {
       /////////////////////// ** Reweighting TrackIso ** /////////////////////////
	 // Will use flat SF and uncertainty
	 /*
       double scaleFIsoL = -1.45573e-03*TrkLeadingMuon4.Pt()+9.37871e-01;
       double scaleFIsoT = -1.45573e-03*TrkTrailingMuon4.Pt()+9.37871e-01;
       
       /////////////////////// ** Uncertainty Studies ** /////////////////////////
       // Track Iso uncertainty studies
       double derivativeL[2];
       derivativeL[0]=1;derivativeL[1]=TrkLeadingMuon4.Pt();
       double derivativeT[2];
       derivativeT[0]=1;derivativeT[1]=TrkTrailingMuon4.Pt();
       
       double ParError[2];
       ParError[0]=7.38139e-02;ParError[1]=4.11611e-03;
       double CorrMatrix[2][2];
       CorrMatrix[0][0]=1.0;CorrMatrix[1][1]=1.0;CorrMatrix[1][0]=-0.919;CorrMatrix[0][1]=-0.919;
       
       double sumdeltaL=0,sumdeltaT=0;
       
       for(int npar=0;npar<2;npar++)
       {
       for(int mpar=0;mpar<2;mpar++)
       {
        sumdeltaL=sumdeltaL+derivativeL[npar]*derivativeL[mpar]*ParError[npar]*ParError[mpar]*CorrMatrix[npar][mpar];
        sumdeltaT=sumdeltaT+derivativeT[npar]*derivativeT[mpar]*ParError[npar]*ParError[mpar]*CorrMatrix[npar][mpar];
       }
       } 
       double DeltascaleFIsoL = TMath::Sqrt(sumdeltaL);
       double DeltascaleFIsoT = TMath::Sqrt(sumdeltaT);

       
       double scaleFIsoLUp = scaleFIsoL+DeltascaleFIsoL;
       double scaleFIsoLDown = scaleFIsoL-DeltascaleFIsoL;
       double scaleFIsoTUp = scaleFIsoT+DeltascaleFIsoT;
       double scaleFIsoTDown = scaleFIsoT-DeltascaleFIsoT;

       double weightUp = weight*scaleFIsoLUp*scaleFIsoTUp;
       double weightDown = weight*scaleFIsoLDown*scaleFIsoTDown;
       
	 */       
       weight *= trkIsoSF*trkIsoSF;
       
       ////**** Setting the values of the variables *****/////
       Eventweight = weight;
       //       Eventweight_TrkIso_Up=weightUp;
       //       Eventweight_TrkIso_Down=weightDown;
       Eventweight_TrkIso_Up=weight;
       Eventweight_TrkIso_Down=weight;
       
       }

       //////////////////////////////////////////////////////////////////////////
       /////////////////////// ** Filling Histogrmas ** /////////////////////////
       //////////////////////////////////////////////////////////////////////////
       tree_Sel->Fill();
        
       counter_FinalEventsH->Fill(1.,weight);
      

     }


     
     
   } /////// End Event Loop ........ ////////
   
   delete tree_;
   file_->Close();
   delete file_;
   
  } // End Filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
} 

 

