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

using namespace std;
double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));
  return Eta;

}

double dPhiFrom2P(double Px1, double Py1,
                  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);

  double cosDPhi = prod/(mod1*mod2);

  return TMath::ACos(cosDPhi);

}

double deltaR(double Eta1, double Phi1,
              double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

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


const float MuMass = 0.105658367;
const float PionMass = 0.13957;

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : mutrk_mass [config_file] [file_list]" << std::endl;
    exit(1);
  }


  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");

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
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);

  // histograms a selection

  TH1D * ptMuH = new TH1D("ptTrailingMuH","",400,0,400);
  TH1D * etaMuH = new TH1D("etaMuH","",48,-2.4,2.4);

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
  float partonMomBins[4] = {0,50,100,1000000};

  TString partonFlavor[4] = {"uds","g","c","b"};
  TString muonPartonNetCharge[2] = {"same","opposite"};
  TString partonMomRange[3] = {"Lt50","50to100","Gt100"};

  // Signal region
  TH1D * InvMassHighMuIsoH = new TH1D("InvMassHighMuIsoH","",20,0.,20.);
  TH1D * InvMassHighMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuIsoH = new TH1D("InvMassLowMuIsoH","",20,0.,20.);
  TH1D * InvMassLowMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuIsoFlavorChargeMomH[4][2][10];

  // LooseIso region
  TH1D * InvMassHighMuLooseIsoH = new TH1D("InvMassHighMuLooseIsoH","",20,0.,20.);
  TH1D * InvMassHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuLooseIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuLooseIsoH = new TH1D("InvMassLowMuLooseIsoH","",20,0.,20.);
  TH1D * InvMassLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuLooseIsoFlavorChargeMomH[4][2][10];

  // Sb region
  TH1D * InvMassHighMuSbH = new TH1D("InvMassHighMuSbH","",20,0.,20.);
  TH1D * InvMassHighMuSbFlavorChargeH[4][2];
  TH1D * PartonMomHighMuSbFlavorChargeH[4][2];
  TH1D * InvMassHighMuSbFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuSbH = new TH1D("InvMassLowMuSbH","",20,0.,20.);
  TH1D * InvMassLowMuSbFlavorChargeH[4][2];
  TH1D * PartonMomLowMuSbFlavorChargeH[4][2];
  TH1D * InvMassLowMuSbFlavorChargeMomH[4][2][10];

  TH1D * partonMomFlavor[4];
  TH1D * partonMomHighMuMatchedFlavorCharge[4][2];
  TH1D * partonMomHighMuFlavorCharge[4][2];
  TH1D * partonMomLowMuMatchedFlavorCharge[4][2];
  TH1D * partonMomLowMuFlavorCharge[4][2];

  for (int iF=0; iF<4; ++iF) {
    for (int iQ=0; iQ<2; ++iQ) {

      InvMassHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      InvMassLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",20,0.,20.);
      PartonMomHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      
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

   tree_->SetBranchAddress("pfjet_count", &pfjet_count);
   tree_->SetBranchAddress("pfjet_e", pfjet_e);
   tree_->SetBranchAddress("pfjet_px", pfjet_px);
   tree_->SetBranchAddress("pfjet_py", pfjet_py);
   tree_->SetBranchAddress("pfjet_pz", pfjet_pz);
   tree_->SetBranchAddress("pfjet_pt", pfjet_pt);
   tree_->SetBranchAddress("pfjet_eta", pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi", pfjet_phi);
   tree_->SetBranchAddress("pfjet_flavour", pfjet_flavour);



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
     

       bool muHighPassed = muon_pt[index]>ptMuonHighCut && fabs(muon_eta[index])<etaMuonHighCut;
       bool muLowPassed  = muon_pt[index]>ptMuonLowCut  && fabs(muon_eta[index])<etaMuonLowCut;

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
       float dRmin = 0.5;
       int flavour = 0;
       float qnet    = 0.0;
       for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	 float drJetMuon = deltaR(muon_eta[index],muon_phi[index],
				  pfjet_eta[jet],pfjet_phi[jet]);
	 if (drJetMuon<dRmin) {
	   dRmin = drJetMuon;
	   int absFlav = abs(pfjet_flavour[jet]);
	   if (absFlav<=3) 
	     flavour = 0;
	   else if (absFlav==4)
	     flavour = 1;
	   else if (absFlav==5)
	     flavour = 2;
	   else if (absFlav==21)
	     flavour = 3;
	   else 
	     flavour = 0;
	   qnet = float(muon_charge[index])*float(pfjet_flavour[jet]);
	 }
       }
       int net = 0;
       if (qnet>0.0) net = 1;

     
       // counting tracks around each muon
       std::vector<unsigned int> trkMu; trkMu.clear(); // all tracks
       std::vector<unsigned int> trkSignalMu; trkSignalMu.clear(); // signal tracks
       std::vector<unsigned int> trkSoftMu; trkSoftMu.clear(); // tracks control region
     
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
	   }
	 }
       }
       
       bool sigMu = trkSignalMu.size()==1 && trkMu.size()==1;

       bool bkgdMu = sigMu || 
       (trkSignalMu.size()==1 && trkSoftMu.size()==1 && trkMu.size()==2) ||
       (trkSignalMu.size()==1 && trkSoftMu.size()==2 && trkMu.size()==3);


     }
     
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

 

