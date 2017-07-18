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
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/src/Config.cc"
#include "TRandom.h"
#include "TRandom3.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/QCDModel.h"

using namespace std;

// obsolete parameterization
double MuLegEfficiency(double pt, double eta, double ptThres, double etaThres) {

  double absEta = fabs(eta);

  if (absEta>etaThres) return 0;

  double effEtaLt0p9 = 0.932;
  double effEta0p9to1p2 = 0.922;
  double effEtaGt1p2 = 0.950;

  double eff = 1.0;
  double ptThresLow = ptThres-1.0;
  double ptThresHigh = ptThres+1.0;
  if (ptThres>30) {
    ptThresLow = ptThres-3.0;
    ptThresHigh = ptThres+3.0;
    effEtaLt0p9 = 0.919;
    effEta0p9to1p2 = 0.841;
    effEtaGt1p2 = 0.866;
  }
  else if (ptThres>16) {
    effEtaLt0p9 =  0.931;
    effEta0p9to1p2 = 0.919;
    effEtaGt1p2 = 0.926;
  }

  if (pt<ptThresLow) {
    eff = 0;
  }
  else if (pt>=ptThresLow&&pt<ptThresHigh) {
    if (absEta<0.9) 
      eff = 0.5*effEtaLt0p9;
    else if (absEta>=0.9&&absEta<1.2)
      eff = 0.5*effEta0p9to1p2;
    else
      eff = 0.5*effEtaGt1p2;
  }
  else {
    if (absEta<0.9)
      eff = effEtaLt0p9;
    else if (absEta>=0.9&&absEta<1.2)
      eff = effEta0p9to1p2;
    else
      eff = effEtaGt1p2;
  }
  
  return eff;
  
}

const float PionMass = 0.13957;

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : mutrk_mass [config_file] [file_list]" << std::endl;
    exit(1);
  }


  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool debug  = cfg.get<bool>("Debug");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");

  // topological cuts
  const float dRMuonsCut   = cfg.get<float>("dRMuonsCut");
  const bool sameSign      = cfg.get<bool>("SameSignMuons");

  // track selection
  const float dRIsoMuon       = cfg.get<float>("dRIsoMuon");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");

  // trigger
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const string dimuonTriggerName = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName1 = cfg.get<string>("MuonLowPtFilterName1");
  const string muonLowPtFilterName2 = cfg.get<string>("MuonLowPtFilterName2");
  const string dimuonDzFilterName = cfg.get<string>("DimuonDzFilterName");
  const string dimuonSameSignFilterName = cfg.get<string>("DimuonSameSignFilterName");
  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
  const float effDzSS        = cfg.get<float>("effDzSS");
  const unsigned int numberOfMuons = cfg.get<unsigned int>("NumberOfMuons");

  TString DiMuonTriggerName(dimuonTriggerName);
  TString MuonHighPtFilterName(muonHighPtFilterName);
  TString MuonLowPtFilterName1(muonLowPtFilterName1);
  TString MuonLowPtFilterName2(muonLowPtFilterName2);
  TString DiMuonDzFilterName(dimuonDzFilterName);
  TString DiMuonSameSignFilterName(dimuonSameSignFilterName);

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
  const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");


  // ********** end of configuration *******************

  std::ifstream fileList(argv[2]);

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


  float genweight;

  float metx;
  float mety;
  float met;
  float metphi;

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
  
   // Trigger
  unsigned int trigobject_count;
  float trigobject_px[1000];
  float trigobject_py[1000];
  float trigobject_pz[1000];
  float trigobject_pt[1000];
  float  trigobject_eta[1000];
  float trigobject_phi[1000];
  bool trigobject_filters[1000][200];

  float numtruepileupinteractions;

  //unsigned int iLeadingPosTrig = 0;
  //vector<bool> trigobject_filter; trigobject_filter.clear();

  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();

  std::string rootFileName(argv[2]);
  
  std::string chainName("makeroottree/AC1B");
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
  
   // Muons
  TH1D * muonCountH = new TH1D("muonCountH","",11,-0.5,10.5);

  TH1D * nGoodMuonsH = new TH1D("nGoodMuonsH","",11,-0.5,10.5);
  TH1D * nGoodIsoMuonsH = new TH1D("nGoodIsoMuonsH","",11,-0.5,10.5);
  TH1D * nGoodLooseIsoMuonsH = new TH1D("nGoodLooseIsoMuonsH","",11,-0.5,10.5);
  TH1D * nGoodSbMuonsH = new TH1D("nGoodSbMuonsH","",11,-0.5,10.5);

  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);

  // histograms a selection

  TH1D * ptMuH = new TH1D("ptMuH","",400,0,400);
  TH1D * etaMuH = new TH1D("etaMuH","",48,-2.4,2.4);
  TH1D * dxyMuH = new TH1D("dxyMuH","",200,-0.5,0.5);
  TH1D * dzMuH = new TH1D("dzMuH","",200,-1,1);

  TH1D * nTracksMuH = new TH1D("nTracksMuH","",21,-0.5,20.5);
  TH1D * nSignalTracksMuH = new TH1D("nSignalTracksMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksMuH = new TH1D("nSoftTracksMuH","",21,-0.5,20.5);
  
  // isolated muon-tracks
  TH1D * ptIsoTrackH = new TH1D("ptIsoTrackH","",100,0,100);
  TH1D * etaIsoTrackH = new TH1D("etaIsoTrackH","",48,-2.4,2.4);
  TH1D * dxyIsoTrackH = new TH1D("dxyIsoTrackH","",200,-0.5,0.5);
  TH1D * dzIsoTrackH = new TH1D("dzIsoTrackH","",200,-1,1);

  TH1D * ptIsoMuH = new TH1D("ptIsoMuH","",400,0,400);
  TH1D * etaIsoMuH = new TH1D("etaIsoMuH","",48,-2.4,2.4);
  TH1D * dxyIsoMuH = new TH1D("dxyIsoMuH","",200,-0.5,0.5);
  TH1D * dzIsoMuH = new TH1D("dzIsoMuH","",200,-1,1);

  TH1D * deltaRMuTrkIsoH = new TH1D("deltaRMuTrkIsoH","",100,0,1.0);
  
  // muon-tracks from Bkgd CR "LooseIso"
  TH1D * ptLooseIsoTrackH = new TH1D("ptLooseIsoTrackH","",100,0,100);
  TH1D * etaLooseIsoTrackH = new TH1D("etaLooseIsoTrackH","",48,-2.4,2.4);
  TH1D * dxyLooseIsoTrackH = new TH1D("dxyLooseIsoTrackH","",200,-0.5,0.5);
  TH1D * dzLooseIsoTrackH = new TH1D("dzLooseIsoTrackH","",200,-1,1);

  TH1D * ptLooseIsoMuH = new TH1D("ptLooseIsoMuH","",400,0,400);
  TH1D * etaLooseIsoMuH = new TH1D("etaLooseIsoMuH","",48,-2.4,2.4);
  TH1D * dxyLooseIsoMuH = new TH1D("dxyLooseIsoMuH","",200,-0.5,0.5);
  TH1D * dzLooseIsoMuH = new TH1D("dzLooseIsoMuH","",200,-1,1);

  TH1D * deltaRMuTrkLooseIsoH = new TH1D("deltaRMuTrkLooseIsoH","",100,0,1.0);

  // muon-tracks from Bkgd CR "LooseIso" and signal region "Iso"
  // adopted name Sb
  TH1D * ptSbTrackH = new TH1D("ptSbTrackH","",100,0,100);
  TH1D * etaSbTrackH = new TH1D("etaSbTrackH","",48,-2.4,2.4);
  TH1D * dxySbTrackH = new TH1D("dxySbTrackH","",200,-0.5,0.5);
  TH1D * dzSbTrackH = new TH1D("dzSbTrackH","",200,-1,1);

  TH1D * ptSbMuH = new TH1D("ptSbMuH","",400,0,400);
  TH1D * etaSbMuH = new TH1D("etaSbMuH","",48,-2.4,2.4);
  TH1D * dxySbMuH = new TH1D("dxySbMuH","",200,-0.5,0.5);
  TH1D * dzSbMuH = new TH1D("dzSbMuH","",200,-1,1);

  TH1D * deltaRMuTrkSbH = new TH1D("deltaRMuTrkSbH","",100,0,1.0);

  int nPartonMomBins = 3;
  float partonMomBins[4] = {0,50,100,150};

  TString partonFlavor[4] = {"uds","g","c","b"};
  TString muonPartonNetCharge[2] = {"opposite","same"};
  TString partonMomRange[3] = {"Lt50","50to100","Gt100"};

  TH1D * PartonMomBinsH = new TH1D("PartonMomBinsH","",nPartonMomBins,partonMomBins);
  TH1D * MuonPartonNetChargeH = new TH1D("MuonPartonNetChargeH","",2,0,2);
  TH1D * PartonFlavorH = new TH1D("PartonFlavorH","",4,0,4);
  for (int iB=0; iB<nPartonMomBins; ++iB)
    PartonMomBinsH->GetXaxis()->SetBinLabel(iB+1,partonMomRange[iB]);
  for (int iB=0; iB<4; ++iB)
    PartonFlavorH->GetXaxis()->SetBinLabel(iB+1,partonFlavor[iB]);
  for (int iB=0; iB<2; ++iB)
    MuonPartonNetChargeH->GetXaxis()->SetBinLabel(iB+1,muonPartonNetCharge[iB]);


  // Signal region
  TH1D * InvMassHighMuIsoH = new TH1D("InvMassHighMuIsoH","",20,0.,20.);
  TH1D * ModelInvMassHighMuIsoH = new TH1D("ModelInvMassHighMuIsoH","",20,0.,20.);
  TH1D * InvMassHighMuIsoAllH = new TH1D("InvMassHighMuIsoAllH","",20,0.,20.);
  TH1D * InvMassHighMuIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuIsoH = new TH1D("InvMassLowMuIsoH","",20,0.,20.);
  TH1D * ModelInvMassLowMuIsoH = new TH1D("ModelInvMassLowMuIsoH","",20,0.,20.);
  TH1D * InvMassLowMuIsoAllH = new TH1D("InvMassLowMuIsoAllH","",20,0.,20.);
  TH1D * InvMassLowMuIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuIsoFlavorChargeMomH[4][2][10];

  // LooseIso region
  TH1D * InvMassHighMuLooseIsoH = new TH1D("InvMassHighMuLooseIsoH","",20,0.,20.);
  TH1D * ModelInvMassHighMuLooseIsoH = new TH1D("ModelInvMassHighMuLooseIsoH","",20,0.,20.);
  TH1D * InvMassHighMuLooseIsoAllH = new TH1D("InvMassHighMuLooseIsoAllH","",20,0.,20.);
  TH1D * InvMassHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuLooseIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuLooseIsoH = new TH1D("InvMassLowMuLooseIsoH","",20,0.,20.);
  TH1D * ModelInvMassLowMuLooseIsoH = new TH1D("ModelInvMassLowMuLooseIsoH","",20,0.,20.);
  TH1D * InvMassLowMuLooseIsoAllH = new TH1D("InvMassLowMuLooseIsoAllH","",20,0.,20.);
  TH1D * InvMassLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuLooseIsoFlavorChargeMomH[4][2][10];

  // Sb region (Signal + Background region)
  TH1D * InvMassHighMuSbH = new TH1D("InvMassHighMuSbH","",20,0.,20.);
  TH1D * ModelInvMassHighMuSbH = new TH1D("ModelInvMassHighMuSbH","",20,0.,20.);
  TH1D * InvMassHighMuSbAllH = new TH1D("InvMassHighMuSbAllH","",20,0.,20.);
  TH1D * InvMassHighMuSbFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuSbFlavorChargeH[4][2];
  TH1D * PartonMomHighMuSbFlavorChargeH[4][2];
  TH1D * InvMassHighMuSbFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuSbH = new TH1D("InvMassLowMuSbH","",20,0.,20.);
  TH1D * ModelInvMassLowMuSbH = new TH1D("ModelInvMassLowMuSbH","",20,0.,20.);
  TH1D * InvMassLowMuSbAllH = new TH1D("InvMassLowMuSbAllH","",20,0.,20.);
  TH1D * InvMassLowMuSbFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuSbFlavorChargeH[4][2];
  TH1D * PartonMomLowMuSbFlavorChargeH[4][2];
  TH1D * InvMassLowMuSbFlavorChargeMomH[4][2][10];

  // momentum
  TH1D * PartonMomFlavorH[4];
  TH1D * PartonMomHighMuFlavorChargeH[4][2];
  TH1D * PartonMomLowMuFlavorChargeH[4][2];

  TH1D * PartonMultiplicityH = new TH1D("PartonMultiplicityH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityH  = new TH1D("PFJetMultiplicityH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityH = new TH1D("GenJetMultiplicityH","",20,-0.5,19.5);

  TH1D * PartonMultiplicityMuH = new TH1D("PartonMultiplicityMuH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuH  = new TH1D("PFJetMultiplicityMuH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuH = new TH1D("GenJetMultiplicityMuH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuH = new TH1D("deltaRPartonMuH","",100,0.,1.);

  TH1D * PartonMultiplicityMuIsoH = new TH1D("PartonMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuIsoH  = new TH1D("PFJetMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuIsoH = new TH1D("GenJetMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuIsoH = new TH1D("deltaRPartonMuIsoH","",100,0.,1.);

  TH1D * PartonMultiplicityMuLooseIsoH = new TH1D("PartonMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuLooseIsoH  = new TH1D("PFJetMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuLooseIsoH = new TH1D("GenJetMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuLooseIsoH = new TH1D("deltaRPartonMuLooseIsoH","",100,0.,1.);

  TH1D * PartonMultiplicityMuSbH = new TH1D("PartonMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuSbH  = new TH1D("PFJetMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuSbH = new TH1D("GenJetMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuSbH = new TH1D("deltaRPartonMuSbH","",100,0.,1.);
  
  TH1D * InvMassH = new TH1D("InvMassH","",20,0.,20.);
  TH2D * InvMass2DH = new TH2D("InvMass2DH","",20,0.,20.,20,0.,20.);

  // Correlation Plots
  TH1D * InvMassTrackPlusMuon1D_ControlH = new TH1D("InvMassTrackPlusMuon1D_ControlH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlH = new TH2D("InvMassTrackPlusMuon2D_ControlH","",20,0.,20.,20,0.,20.);

  TH1D * InvMassTrackPlusMuon1D_ControlXH = new TH1D("InvMassTrackPlusMuon1D_ControlXH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlXH = new TH2D("InvMassTrackPlusMuon2D_ControlXH","",20,0.,20.,20,0.,20.);
   
  TH1D * InvMassTrackPlusMuon1D_ControlYH = new TH1D("InvMassTrackPlusMuon1D_ControlYH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlYH = new TH2D("InvMassTrackPlusMuon2D_ControlYH","",20,0.,20.,20,0.,20.);

  for (int iF=0; iF<4; ++iF) {

    PartonMomFlavorH[iF] = new TH1D("partonMomFlavor_"+partonFlavor[iF],"",500,0.,5000.);
    

    for (int iQ=0; iQ<2; ++iQ) {

      PartonMomHighMuFlavorChargeH[iF][iQ] = new TH1D("partonMomHighMu_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuFlavorChargeH[iF][iQ] = new TH1D("partonMomLowMu_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      
      InvMassHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      InvMassLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      PartonMomHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      InvMassHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      InvMassLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      PartonMomHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      InvMassHighMuSbFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      InvMassLowMuSbFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassHighMuSbFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      ModelInvMassLowMuSbFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      PartonMomHighMuSbFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuSbFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      for (int iM=0; iM<nPartonMomBins; ++iM) {

	InvMassLowMuIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);
	InvMassLowMuLooseIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);
	InvMassLowMuSbFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);

	InvMassHighMuIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);
	InvMassHighMuLooseIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);
	InvMassHighMuSbFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",20,0.,20.);

      }      
      
    }
  }

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

  // QCD Model
  TString fileNameQCDModel = TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/QCDModel.root");
  QCDModel * qcdModel = new QCDModel(fileNameQCDModel);


  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Trigger efficiencies

  ScaleFactor * SF_muon17 = new ScaleFactor();
  SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));
  ScaleFactor * SF_muon8 = new ScaleFactor();
  SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));

  TString filen;
  int iFiles = 0;
  int events = 0;
  while (fileList >> filen) {
   iFiles++;
   cout << "file " << iFiles << " : " << filen << endl;
   
   TFile * file_ = TFile::Open(TString(filen));
   if (file_==NULL) continue;

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
   tree_->SetBranchAddress("pfmet_ex", &metx);
   tree_->SetBranchAddress("pfmet_ey", &mety);
   tree_->SetBranchAddress("pfmet_pt", &met);
   tree_->SetBranchAddress("pfmet_phi",&metphi);
 

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

   // genjets
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

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   //   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

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
   }   

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     if (!isData) {
       weight *= genweight;
     }


     histWeightsH->Fill(1.0,weight);

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

     if (debug) {
       std::cout << std::endl;
       std::cout << "+++++++++++++++++++++++++++++" << std::endl;
       std::cout << std::endl;
       std::cout << "genjets : " << genjets_count << std::endl;
     }
     int partons = 0;
     int pfjets  = 0;
     int genjets = 0;
     if (!isData) {
       for (unsigned int igen=0; igen<genjets_count; ++igen) {
	 int pdgId = genjets_pdgid[igen];
	 TLorentzVector genjetLV; genjetLV.SetXYZT(genjets_px[igen],
						   genjets_py[igen],
						   genjets_pz[igen],
						   genjets_e[igen]);
       
	 if (genjetLV.Pt()>10&&fabs(genjetLV.Eta())<3.0) {
	   genjets++;
	   if (debug)
	     printf("  flavor = %3i   pT = %7.2f   eta = %5.2f   phi = %5.2f   status = %3i\n",
		    genjets_pdgid[igen],genjetLV.Pt(),genjetLV.Eta(),genjetLV.Phi(),genjets_status[igen]);
	 }
       }
     }
     if (debug) {
       std::cout << std::endl;
       std::cout << "jets -> " << std::endl;
     }
     std::vector<int> partonPdgId; partonPdgId.clear();
     std::vector<TLorentzVector> partonLV; partonLV.clear();
     for (unsigned ijet=0; ijet<pfjet_count; ++ijet) {
       if (pfjet_pt[ijet]>10&&fabs(pfjet_eta[ijet])<3.0) {
	 pfjets++;
	 if (pfjet_flavour[ijet]!=0&&!isData) { 
	   partons++;
	   int absPdgId = TMath::Abs(pfjet_flavour[ijet]);
	   int iflav = 0;
	   if (absPdgId==21) iflav = 1;
	   if (absPdgId==4)  iflav = 2;
	   if (absPdgId==5)  iflav = 3;
	   TLorentzVector jetLV; jetLV.SetXYZT(pfjet_px[ijet],
					       pfjet_py[ijet],
					       pfjet_pz[ijet],
					       pfjet_e[ijet]);
	   TLorentzVector partLV = jetLV;
	   float dRMin = 0.5;
	   for (unsigned int igen=0; igen<genjets_count; ++igen) {
	     TLorentzVector genjetLV; genjetLV.SetXYZT(genjets_px[igen],
						       genjets_py[igen],
						       genjets_pz[igen],
						       genjets_e[igen]);
	     float dRJets = deltaR(jetLV.Eta(),jetLV.Phi(),
				   genjetLV.Eta(),genjetLV.Phi());
	     if (dRJets<dRMin) {
	       dRMin = dRJets;
	       partLV = genjetLV;
	     }	     
	   }
	   partonPdgId.push_back(pfjet_flavour[ijet]);
	   partonLV.push_back(partLV);
	   PartonMomFlavorH[iflav]->Fill(partLV.P(),weight);
	 }
	 if (debug)
	   printf("  flavor = %3i   pT = %7.2f   eta = %5.2f   phi = %5.2f\n",
		  pfjet_flavour[ijet],pfjet_pt[ijet],pfjet_eta[ijet],pfjet_phi[ijet]);
       
       }
     }
    
     // filling histograms 
     PartonMultiplicityH->Fill(float(partons),weight);
     PFJetMultiplicityH->Fill(float(pfjets),weight);
     GenJetMultiplicityH->Fill(float(genjets),weight);

     // finding HLT filters in the HLT Filter library
     unsigned int nMu8Leg   = 0;
     unsigned int nMu17Leg  = 0;
     bool isMu8Leg = false;
     bool isMu17Leg = false;

     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       if (HLTFilter==MuonHighPtFilterName) {
	 nMu17Leg = i;
	 isMu17Leg = true;
       }
       if (HLTFilter==MuonLowPtFilterName1||HLTFilter==MuonLowPtFilterName2) {
	 nMu8Leg = i;
	 isMu8Leg = true;
       }
     }
     
     if (!isMu17Leg) {
       cout << "Filter " << MuonHighPtFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isMu8Leg) {
       cout << "Filters " << MuonLowPtFilterName1 
	    << " or " << MuonLowPtFilterName2
	    << " not found " << endl;
       exit(-1);
     }


     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       if (isData) {
	 if (event_run >= 278820 && muon_isMedium[i]) muonID = true; // Run2016G-H
	 if (event_run < 278820  && muon_isICHEP[i]) muonID = true; // Run2016B-F
       }
       if (!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       if(muon_pt[i]<ptMuonLowCut) continue;
       if(fabs(muon_eta[i])>etaMuonLowCut) continue;
       //      cout << "muon pt = " << muon_pt[i] << endl;
       muons.push_back(i);
     }
     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<1) continue; // quit event if number of good muons < 1

     int nIsoMuons = 0;
     int nLooseIsoMuons = 0;
     int nSbMuons = 0;

     if (debug) {
       std::cout << std::endl;
       std::cout << "muons -> " << std::endl;
     }

     for (unsigned int imu=0; imu<muons.size(); ++imu) {
       bool muMatchMu8  = false;
       bool muMatchMu17 = false;
       unsigned int index = muons.at(imu);
       for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	 float dRtrig = deltaR(muon_eta[index],muon_phi[index],
			       trigobject_eta[iT],trigobject_phi[iT]);
	 if (dRtrig>DRTrigMatch) continue;
	 if (trigobject_filters[iT][nMu17Leg])
	   muMatchMu17 = true;
	 if (trigobject_filters[iT][nMu8Leg])
	   muMatchMu8 = true;
       }
       if (!applyTriggerMatch) {
	 muMatchMu17 = true;
	 muMatchMu8 = true;
       }
       

       bool muHighPassed = muon_pt[index]>ptMuonHighCut && fabs(muon_eta[index])<etaMuonHighCut && muMatchMu17;
       bool muLowPassed  = muon_pt[index]>ptMuonLowCut  && fabs(muon_eta[index])<etaMuonLowCut  && muMatchMu8;

       // trigger weight for given muon
       float trigWeightMu17 = 1;
       float trigWeightMu8  = 1;
       if (!isData) { // trigger efficiency here
	   double ptMuon = muon_pt[index];
	   double etaMuon = muon_eta[index];

	   double effMu8data = SF_muon8->get_EfficiencyData(ptMuon,etaMuon); 
	   double effMu8MC   = SF_muon8->get_EfficiencyMC(ptMuon,etaMuon); 
	   double effMu17data = SF_muon17->get_EfficiencyData(ptMuon,etaMuon); 
	   double effMu17MC   = SF_muon17->get_EfficiencyMC(ptMuon,etaMuon); 
	   
	   if (applyTriggerMatch) {
	     if (effMu8MC>0)
	       trigWeightMu8 = effMu8data/effMu8MC;
	     if (effMu17MC>0)
	       trigWeightMu17 = effMu17data/effMu17MC;
	   }
	   else {
	     trigWeightMu8  = effMu8data;
	     trigWeightMu17 = effMu17data;
	   }
       }	   
       // Muon
       TLorentzVector Muon4; Muon4.SetXYZM(muon_px[index],
					   muon_py[index],
					   muon_pz[index],
					   MuMass);
       
       // determine flavour of jet
       float dRmin = 1.0;
       int flavour = 0;
       float qnet  = 0.0;
       int pdgId = 0;
       bool matchedParton = false;
       TLorentzVector matchedPartonLV;
       for (unsigned int ip=0; ip<partonPdgId.size(); ++ip) {
	 TLorentzVector partLV = partonLV.at(ip);
	 float drJetMuon = deltaR(muon_eta[index],muon_phi[index],
				  partLV.Eta(),partLV.Phi());
	 if (drJetMuon<dRmin) {
	   dRmin = drJetMuon;
	   int absFlav = TMath::Abs(partonPdgId.at(ip));
	   pdgId = partonPdgId.at(ip);
	   flavour = 0;
	   if (absFlav==21) 
	     flavour = 1;
	   if (absFlav==4)
	     flavour = 2;
	   if (absFlav==5)
	     flavour = 3;
	   qnet = float(muon_charge[index])*float(pdgId);
	   matchedParton = true;
	   matchedPartonLV = partLV;
	 }
       }

       int net = 0;
       if (qnet>0.0) net = 1;
       if (flavour==1) {
	 if (muon_charge[index]<0) net = 0;
	 if (muon_charge[index]>0) net = 1;
       }

       // counting tracks around each muon
       std::vector<unsigned int> trkMu; trkMu.clear(); // all tracks
       std::vector<unsigned int> trkSignalMu; trkSignalMu.clear(); // signal tracks
       std::vector<unsigned int> trkSoftMu; trkSoftMu.clear(); // tracks control region

       TLorentzVector trackLV; trackLV.SetXYZM(0,0,0,PionMass);
     
       for (unsigned int iTrk=0; iTrk<track_count; ++iTrk) {
	 if (fabs(track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
	 if (fabs(track_dxy[iTrk])>dxyTrkLooseCut) continue;
	 if (fabs(track_dz[iTrk])>dzTrkLooseCut) continue;
	 if (fabs(track_eta[iTrk])>etaTrkCut) continue;
	 if (fabs(track_pt[iTrk])<ptTrkLooseCut) continue;
       
	 TLorentzVector trk4; trk4.SetXYZM(track_px[iTrk],
					   track_py[iTrk],
					   track_pz[iTrk],
					   track_mass[iTrk]);
	 
	 TLorentzVector MuDiff = Muon4 - trk4;
	 if (MuDiff.P()>0.1) { // track is not leading muon
	   float drTrkMu = deltaR(muon_eta[index],muon_phi[index],
				  track_eta[iTrk],   track_phi[iTrk]);
	   float qTrkMu = track_charge[iTrk]*muon_charge[index];
	   if (drTrkMu<dRIsoMuon){
	     trkMu.push_back(iTrk);
	     if (track_pt[iTrk]>ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	       trkSoftMu.push_back(iTrk);
	   }
	   if (drTrkMu<dRIsoMuon && qTrkMu<0 
	       && fabs(track_dxy[iTrk])<dxyTrkCut 
	       && fabs(track_dz[iTrk])<dzTrkCut 
	       && track_pt[iTrk]>ptTrkCut) {
	     trkSignalMu.push_back(iTrk);
	     trackLV = trk4;
	   }
	 }
       }

       TLorentzVector muonTrkLV = trackLV + Muon4;
       float muonTrkMass = muonTrkLV.M();
       float deltaRMuonTrk = 1.1;
       unsigned int indexTrk = 0;
       if (trkSignalMu.size()==1) { 
	 indexTrk = trkSignalMu.at(0);
	 deltaRMuonTrk = deltaR(Muon4.Eta(),Muon4.Phi(),
				trackLV.Eta(),trackLV.Phi());
       }

       if (debug) {
	 printf("  pT = %6.2f   eta = %5.2f   phi = %5.2f   jetFlavor = %3i  nTrk = %2i  nSigTrk = %2i\n",
		Muon4.Pt(),Muon4.Eta(),Muon4.Phi(),pdgId,int(trkMu.size()),int(trkSignalMu.size()));
       }

       bool sigMu = trkSignalMu.size()==1 && trkMu.size()==1;

       bool bkgdMu = 
	 (trkSignalMu.size()==1 && trkSoftMu.size()==1 && trkMu.size()==2) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==2 && trkMu.size()==3) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==3 && trkMu.size()==4);

       bool sbMu = sigMu || bkgdMu;
       
       float weightHigh = weight * trigWeightMu17;
       float weightLow  = weight * trigWeightMu8;

       PartonMultiplicityMuH->Fill(float(partons),weightLow);
       PFJetMultiplicityMuH->Fill(float(pfjets),weightLow);
       GenJetMultiplicityMuH->Fill(float(genjets),weightLow);

       ptMuH->Fill(muon_pt[index],weightLow);
       etaMuH->Fill(muon_eta[index],weightLow);
       dxyMuH->Fill(muon_dxy[index],weightLow);
       dzMuH->Fill(muon_dz[index],weightLow);
       
       nTracksMuH->Fill(float(trkMu.size()),weightLow);
       nSoftTracksMuH->Fill(float(trkSoftMu.size()),weightLow);
       nSignalTracksMuH->Fill(float(trkSoftMu.size()),weightLow);

       int partonMomBin = 0;
       if (matchedParton) {
	 partonMomBin = binNumber(TMath::Min(float(matchedPartonLV.P()),float(partonMomBins[nPartonMomBins]-0.1)),nPartonMomBins,partonMomBins);
	 deltaRPartonMuH->Fill(dRmin,weightLow);
	 if (muHighPassed) PartonMomHighMuFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightHigh);
	 if (muLowPassed) PartonMomLowMuFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightLow);
	 if (muHighPassed) {
	   for (int iM=0; iM<20; ++iM) {
	     double mass = double(iM) + double(0.5);
	     int muType = 0;
	     int ireg = 0;
	     double pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
	     ModelInvMassHighMuIsoH->Fill(mass,weightHigh*pdf);
	     ModelInvMassHighMuIsoFlavorChargeH[flavour][net]->Fill(mass,weightHigh*pdf);
	     ireg = 1;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassHighMuLooseIsoH->Fill(mass,weightHigh*pdf);
             ModelInvMassHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(mass,weightHigh*pdf);
	     ireg = 2;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassHighMuSbH->Fill(mass,weightHigh*pdf);
             ModelInvMassHighMuSbFlavorChargeH[flavour][net]->Fill(mass,weightHigh*pdf);
	   }
	 }
	 if (muLowPassed) {
	   for (int iM=0; iM<20; ++iM) {
	     double mass = double(iM) + double(0.5);
	     int muType = 1;
	     int ireg = 0;
	     double pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
	     ModelInvMassLowMuIsoH->Fill(mass,weightHigh*pdf);
	     ModelInvMassLowMuIsoFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	     ireg = 1;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassLowMuLooseIsoH->Fill(mass,weightHigh*pdf);
             ModelInvMassLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(mass,weightLow*pdf);
	     ireg = 2;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassLowMuSbH->Fill(mass,weightHigh*pdf);
             ModelInvMassLowMuSbFlavorChargeH[flavour][net]->Fill(mass,weightLow*pdf);
	   }
	 }
       }

       if (sigMu) {
	 
	 nIsoMuons++;
	 ptIsoTrackH->Fill(trackLV.Pt(),weightLow);
	 etaIsoTrackH->Fill(trackLV.Eta(),weightLow);
	 dxyIsoTrackH->Fill(track_dxy[indexTrk],weightLow);
	 dzIsoTrackH->Fill(track_dz[indexTrk],weightLow);
	 ptIsoMuH->Fill(Muon4.Pt(),weightLow);
	 etaIsoMuH->Fill(Muon4.Eta(),weightLow);
	 dxyIsoMuH->Fill(muon_dxy[index],weightLow);
	 dzIsoMuH->Fill(muon_dz[index],weightLow);
	 PartonMultiplicityMuIsoH->Fill(float(partons),weightLow);
	 PFJetMultiplicityMuIsoH->Fill(float(pfjets),weightLow);
	 GenJetMultiplicityMuIsoH->Fill(float(genjets),weightLow);
	 deltaRMuTrkIsoH->Fill(deltaRMuonTrk,weightLow);


	 if (muHighPassed) InvMassHighMuIsoAllH->Fill(muonTrkMass,weightHigh);
	 if (muLowPassed)  InvMassLowMuIsoAllH->Fill(muonTrkMass,weightLow);
	 if (matchedParton) {
	   deltaRPartonMuIsoH->Fill(dRmin,weightLow);
	   if (muHighPassed) { 
	     InvMassHighMuIsoH->Fill(muonTrkMass,weightHigh);
	     InvMassHighMuIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightHigh);
	     PartonMomHighMuIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightHigh);
	     InvMassHighMuIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightHigh);
	   }
	   if (muLowPassed) { 
	     InvMassLowMuIsoH->Fill(muonTrkMass,weightLow);
	     InvMassLowMuIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightLow);
             PartonMomLowMuIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightLow);
             InvMassLowMuIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightLow);
	   }
	 }
       }

       if (bkgdMu) {

	 nLooseIsoMuons++;
	 ptLooseIsoTrackH->Fill(trackLV.Pt(),weightLow);
         etaLooseIsoTrackH->Fill(trackLV.Eta(),weightLow);
         dxyLooseIsoTrackH->Fill(track_dxy[indexTrk],weightLow);
         dzLooseIsoTrackH->Fill(track_dz[indexTrk],weightLow);
         ptLooseIsoMuH->Fill(Muon4.Pt(),weightLow);
         etaLooseIsoMuH->Fill(Muon4.Eta(),weightLow);
         dxyLooseIsoMuH->Fill(muon_dxy[index],weightLow);
         dzLooseIsoMuH->Fill(muon_dz[index],weightLow);
         PartonMultiplicityMuLooseIsoH->Fill(float(partons),weightLow);
         PFJetMultiplicityMuLooseIsoH->Fill(float(pfjets),weightLow);
         GenJetMultiplicityMuLooseIsoH->Fill(float(genjets),weightLow);
	 deltaRMuTrkLooseIsoH->Fill(deltaRMuonTrk,weightLow);

	 if (muHighPassed) InvMassHighMuLooseIsoAllH->Fill(muonTrkMass,weightHigh);
         if (muLowPassed)  InvMassLowMuLooseIsoAllH->Fill(muonTrkMass,weightLow);

	 
	 if (matchedParton) {
           deltaRPartonMuLooseIsoH->Fill(dRmin,weightLow);
	   if (muHighPassed) {
             InvMassHighMuLooseIsoH->Fill(muonTrkMass,weightHigh);
             InvMassHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightHigh);
             PartonMomHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightHigh);
             InvMassHighMuLooseIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightHigh);
           }
           if (muLowPassed) {
             InvMassLowMuLooseIsoH->Fill(muonTrkMass,weightLow);
             InvMassLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightLow);
             PartonMomLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightLow);
             InvMassLowMuLooseIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightLow);
           }

         }
       }

       if (sbMu) {

	 nSbMuons++;
	 ptSbTrackH->Fill(trackLV.Pt(),weightLow);
         etaSbTrackH->Fill(trackLV.Eta(),weightLow);
         dxySbTrackH->Fill(track_dxy[indexTrk],weightLow);
         dzSbTrackH->Fill(track_dz[indexTrk],weightLow);
         ptSbMuH->Fill(Muon4.Pt(),weightLow);
         etaSbMuH->Fill(Muon4.Eta(),weightLow);
         dxySbMuH->Fill(muon_dxy[index],weightLow);
         dzSbMuH->Fill(muon_dz[index],weightLow);
         PartonMultiplicityMuSbH->Fill(float(partons),weightLow);
         PFJetMultiplicityMuSbH->Fill(float(pfjets),weightLow);
         GenJetMultiplicityMuSbH->Fill(float(genjets),weightLow);
	 deltaRMuTrkSbH->Fill(deltaRMuonTrk,weightLow);

         if (muHighPassed) InvMassHighMuSbAllH->Fill(muonTrkMass,weightHigh);
         if (muLowPassed)  InvMassLowMuSbAllH->Fill(muonTrkMass,weightLow);
	 if (matchedParton) {
           deltaRPartonMuSbH->Fill(dRmin,weightLow);
           if (muHighPassed) {
             InvMassHighMuSbH->Fill(muonTrkMass,weightHigh);
             InvMassHighMuSbFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightHigh);
             PartonMomHighMuSbFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightHigh);
             InvMassHighMuSbFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightHigh);
           }
           if (muLowPassed) {
             InvMassLowMuSbH->Fill(muonTrkMass,weightLow);
             InvMassLowMuSbFlavorChargeH[flavour][net]->Fill(muonTrkMass,weightLow);
             PartonMomLowMuSbFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weightLow);
             InvMassLowMuSbFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weightLow);
           }
         }
       }

     }

     nGoodIsoMuonsH->Fill(float(nIsoMuons),weight);
     nGoodLooseIsoMuonsH->Fill(float(nLooseIsoMuons),weight);
     nGoodSbMuonsH->Fill(float(nSbMuons),weight);

   }
   delete tree_;
   file_->Close();
   delete file_;
   
  }// filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
}// int main loop 

 

