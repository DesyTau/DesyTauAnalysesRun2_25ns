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

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"


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
    std::cout << "Usage of the program : Hto4TausAnalysis [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(1);
  }


  // **** configuration
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

  //  const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
  //  const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");
  const string Muon17TriggerFile = cfg.get<string>("MuonHighPtTriggerEff");
  const string Muon8TriggerFile = cfg.get<string>("MuonLowPtTriggerEff");

  // Higgs pt reweighting
  const string higgsPtFileName = cfg.get<string>("HiggsPtFileName");
  TString HiggsPtFileName(higgsPtFileName);


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
  
   // Muons
  TH1D * muonCountH = new TH1D("muonCountH","",11,-0.5,10.5);
  TH1D * nGoodMuonsH = new TH1D("nGoodMuonsH","",11,-0.5,10.5);
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);

  // histograms after dimuon selection
  TH1D * ptLeadingMuH = new TH1D("ptLeadingMuH","",400,0,400);
  TH1D * ptTrailingMuH = new TH1D("ptTrailingMuH","",400,0,400);
  TH1D * etaLeadingMuH = new TH1D("etaLeadingMuH","",48,-2.4,2.4);
  TH1D * etaTrailingMuH = new TH1D("etaTrailingMuH","",48,-2.4,2.4);
  TH1D * dimuonMassH = new TH1D("dimuonMassH","",500,0,500);

  TH1D * nTracksLeadingMuH = new TH1D("nTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nTracksTrailingMuH = new TH1D("nTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nSigTracksLeadingMuH = new TH1D("nSigTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nSigTracksTrailingMuH = new TH1D("nSigTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksLeadingMuH = new TH1D("nSoftTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksTrailingMuH = new TH1D("nSoftTracksTrailingMuH","",21,-0.5,20.5);
  
  // isolated muons
  TH1D * ptTrackLeadingMuH = new TH1D("ptTrackLeadingMuH","",100,0,100);
  TH1D * etaTrackLeadingMuH = new TH1D("etaTrackLeadingMuH","",48,-2.4,2.4);
  TH1D * dxyTrackLeadingMuH = new TH1D("dxyTrackLeadingMuH","",200,-0.5,0.5);
  TH1D * dzTrackLeadingMuH = new TH1D("dzTrackLeadingMuH","",200,-1,1);

  TH1D * ptTrackTrailingMuH = new TH1D("ptTrackTrailingMuH","",100,0,100);
  TH1D * etaTrackTrailingMuH = new TH1D("etaTrackTrailingMuH","",48,-2.4,2.4);
  TH1D * dxyTrackTrailingMuH = new TH1D("dxyTrackTrailingMuH","",200,-0.5,0.5);
  TH1D * dzTrackTrailingMuH = new TH1D("dzTrackTrailingMuH","",200,-1,1);

  TH1D * ptTrackN1H = new TH1D("ptTrackN1H","",100,0,100);
  TH1D * etaTrackN1H = new TH1D("etaTrackN1H","",48,-2.4,2.4);
  TH1D * dxyTrackN1H = new TH1D("dxyTrackN1H","",200,-0.5,0.5);
  TH1D * dzTrackN1H = new TH1D("dzTrackN1H","",200,-1,1);

  TH1D * ptTrackH = new TH1D("ptTrackH","",100,0,100);
  TH1D * etaTrackH = new TH1D("etaTrackH","",48,-2.4,2.4);
  TH1D * dxyTrackH = new TH1D("dxyTrackH","",200,-0.5,0.5);
  TH1D * dzTrackH = new TH1D("dzTrackH","",200,-1,1);

  // Signal region
  TH1D * InvMassLeadingH = new TH1D("InvMassLeadingH","",20,0.,20.);
  TH1D * InvMassTrailingH = new TH1D("InvMassTrailingH","",20,0.,20.);
  TH1D * InvMassH = new TH1D("InvMassH","",20,0.,20.);
  TH2D * InvMass2DH = new TH2D("InvMass2DH","",20,0.,20.,20,0.,20.);

  TH1D * MetSelH = new TH1D("MetH","",400,0.,400.);
  TH1D * mTtotSelH = new TH1D("mTtotSelH","",400,0.,400.);

  TH1D * ptLeadingMuSelH = new TH1D("ptLeadingMuSelH","",100,0,100);
  TH1D * ptLeadingMuTrkSelH = new TH1D("ptLeadingMuTrkSelH","",100,0,100);
  TH1D * ptTrailingMuSelH = new TH1D("ptTrailingMuSelH","",100,0,100);
  TH1D * ptTrailingMuTrkSelH = new TH1D("ptTrailingMuTrkSelH","",100,0,100);

  TH1D * ptMuSelH = new TH1D("ptMuSelH","",100,0,100);
  TH1D * ptTrkSelH = new TH1D("ptTrkSelH","",100,0,100);

  TH1D * etaLeadingMuSelH = new TH1D("etaLeadingMuSelH","",48,-2.4,2.4);
  TH1D * etaLeadingMuTrkSelH = new TH1D("etaLeadingMuTrkSelH","",48,-2.4,2.4);
  TH1D * etaTrailingMuSelH = new TH1D("etaTrailingMuSelH","",48,-2.4,2.4);
  TH1D * etaTrailingMuTrkSelH = new TH1D("etaTrailingMuTrkSelH","",48,-2.4,2.4);

  TH1D * etaMuSelH = new TH1D("etaMuSelH","",48,-2.4,2.4);
  TH1D * etaTrkSelH = new TH1D("etaTrkSelH","",48,-2.4,2.4);

  TH1D * dRMuTrkSelH = new TH1D("dRMuTrkSelH","",50,0.,0.5);

  TH1D * dimuonMassSelH = new TH1D("dimuonMassSelH","",500,0,500);
  TH1D * invMass2Mu2TrkSelH = new TH1D("invMass2Mu2TrkSelH","",500,0,500);

  // Counters
  TH1D * counter_InputEventsH=new TH1D("counter_InputEventsH","",1,0.,2.);
  TH1D * counter_MuonSizeGTE2H=new TH1D("counter_MuonSizeGTE2H","",1,0.,2.);
  TH1D * counter_MuonKinematicsH=new TH1D("counter_MuonKinematicsH","",1,0.,2.);         
  TH1D * counter_nMuTrackSigH=new TH1D("counter_nMuTrackSigH","",1,0.,2.);
  TH1D * counter_FinalEventsH=new TH1D("counter_FinalEventsH","",1,0.,2.);         
  TH1D * counter_ControlEventsH=new TH1D("counter_ControlEventsH","",1,0.,2.);         
  TH1D * counter_ControlXEventsH=new TH1D("counter_ControlXEventsH","",1,0.,2.);         
  TH1D * counter_ControlYEventsH=new TH1D("counter_ControlYEventsH","",1,0.,2.);         
  

  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);
  TH1D * histWeightsSingleMuH = new TH1D("histWeightsSingleMuH","",1,0.,2.);
  TH1D * histWeightsTripleMuH = new TH1D("histWeightsTripleMuH","",1,0.,2.);
  TH1D * histWeightsDoubleMuSSH = new TH1D("histWeightsDoubleMuSSH","",1,0.,2.);
  TH1D * histWeightsAllTriggersH = new TH1D("histWeightsAllTriggersH","",1,0.,2.);

  // Background studies
  // N23 
  TH1D * InvMassN23leadingH = new TH1D("InvMassN23leadingH","",20,0.,20.);
  TH1D * InvMassN23trailingH = new TH1D("InvMassN23trailingH","",20,0.,20.);
  TH1D * InvMassN23H = new TH1D("InvMassN23H","",20,0.,20.);
  TH1D * InvMassN45H = new TH1D("InvMassN45H","",20,0.,20.);

  TH1D * ptMuN23H = new TH1D("ptMuN23H","",100,0,100);
  TH1D * ptTrkN23H = new TH1D("ptTrkN23H","",100,0,100);

  TH1D * etaMuN23H = new TH1D("etaMuN23H","",48,-2.4,2.4);
  TH1D * etaTrkN23H = new TH1D("etaTrkN23H","",48,-2.4,2.4);

  TH1D * dRMuTrkN23H = new TH1D("dRMuTrkN23H","",50,0.,0.5);

  // N1trk, N23trk 

  TH1D * InvMassHardestNtrk23leadingH = new TH1D("InvMassHardestNtrk23leadingH","",20,0.,20.);
  TH1D * InvMassHardestNtrk23trailingH = new TH1D("InvMassHardestNtrk23trailingH","",20,0.,20.);
  TH1D * InvMassHardestNtrk23H = new TH1D("InvMassHardestNtrk23H","",20,0.,20.);

  TH1D * InvMassSoftestNtrk23leadingH = new TH1D("InvMassSoftestNtrk23leadingH","",20,0.,20.);
  TH1D * InvMassSoftestNtrk23trailingH = new TH1D("InvMassSoftestNtrk23trailingH","",20,0.,20.);
  TH1D * InvMassSoftestNtrk23H = new TH1D("InvMassSoftestNtrk23H","",20,0.,20.);

  TH1D * InvMassHardestNtrk1leadingH = new TH1D("InvMassHardestNtrk1leadingH","",20,0.,20.);
  TH1D * InvMassHardestNtrk1trailingH = new TH1D("InvMassHardestNtrk1trailingH","",20,0.,20.);
  TH1D * InvMassHardestNtrk1H = new TH1D("InvMassHardestNtrk1H","",20,0.,20.);

  TH1D * InvMassSoftestNtrk1leadingH = new TH1D("InvMassSoftestNtrk1leadingH","",20,0.,20.);
  TH1D * InvMassSoftestNtrk1trailingH = new TH1D("InvMassSoftestNtrk1trailingH","",20,0.,20.);
  TH1D * InvMassSoftestNtrk1H = new TH1D("InvMassSoftestNtrk1H","",20,0.,20.);

  // Correlation Plots
  TH1D * InvMassTrackPlusMuon1D_ControlH = new TH1D("InvMassTrackPlusMuon1D_ControlH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlH = new TH2D("InvMassTrackPlusMuon2D_ControlH","",20,0.,20.,20,0.,20.);

  TH1D * InvMassTrackPlusMuon1D_ControlXH = new TH1D("InvMassTrackPlusMuon1D_ControlXH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlXH = new TH2D("InvMassTrackPlusMuon2D_ControlXH","",20,0.,20.,20,0.,20.);
   
  TH1D * InvMassTrackPlusMuon1D_ControlYH = new TH1D("InvMassTrackPlusMuon1D_ControlYH","",20,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlYH = new TH2D("InvMassTrackPlusMuon2D_ControlYH","",20,0.,20.,20,0.,20.);
   
  // Monte Carlo information
  TH1D * deltaRMuonPionH = new TH1D("deltaRMuonPionH","",200,0,2);
  TH1D * pionPtH = new TH1D("pionPtH","",100,0,100);

  float higgsPt;

  TTree * higgsTree = new TTree("higgsTree","");
  higgsTree->Branch("HiggsPt",&higgsPt,"HiggsPt/F");

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
  TFile * higgsPtFile = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+HiggsPtFileName);
  TH1D * higgsPtH = (TH1D*)higgsPtFile->Get("kfactor");
  //  std::cout << "Higgs Pt histogram : " << higgsPtH << std::endl;

  // Trigger efficiencies
  ScaleFactor * SF_muon17 = new ScaleFactor();
  SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));
  ScaleFactor * SF_muon8 = new ScaleFactor();
  SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));

  // Correction workspace
  TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_5.root";
  TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  TString filen;
  int iFiles = 0;
  int events = 0;
  while (fileList >> filen) {
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
   tree_->SetBranchAddress("pfmetcorr_pt", &met);
   tree_->SetBranchAddress("pfmetcorr_phi",&metphi);
 

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

     std::vector<unsigned int> posPion; posPion.clear();
     std::vector<unsigned int> negPion; negPion.clear();
     std::vector<unsigned int> posMuon; posMuon.clear();
     std::vector<unsigned int> negMuon; negMuon.clear();

     bool HiggsFound = false;
     unsigned int higgsIndex = 0;
     if (!isData) {
       //       std::cout << "Generated particles = " << genparticles_count << std::endl;
       for (unsigned int iP=0; iP<genparticles_count; ++iP) {
	 if (genparticles_status[iP]==1&&genparticles_info[iP]==12) {
	   if (genparticles_pdgid[iP]==13) negMuon.push_back(iP);
	   if (genparticles_pdgid[iP]==-13) posMuon.push_back(iP);
	   if (genparticles_pdgid[iP]==211) posPion.push_back(iP);
	   if (genparticles_pdgid[iP]==-211) negPion.push_back(iP);
	 }
	 if (genparticles_pdgid[iP]==35) {
	   higgsIndex = iP;
	   HiggsFound = true;
	   TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[iP],
						   genparticles_py[iP],
						   genparticles_pz[iP],
						   genparticles_e[iP]);
	   
	 }
       }
       //       if (posMuon.size()==2||negMuon.size()==2) {
       //	 std::cout << "H->aa->4tau : " << std::endl;
       //	 std::cout << "Number of mu-   : " << negMuon.size() << std::endl;
       //	 std::cout << "Number of mu+   : " << posMuon.size() << std::endl;
       //	 std::cout << "Number of pion- : " << negPion.size() << std::endl;
       //	 std::cout << "Number of pion+ : " << posPion.size() << std::endl;
       //	 std::cout << std::endl;
       //       }
       if (posMuon.size()==2&&negMuon.size()==0&&negPion.size()==2&&posPion.size()==0) {
	 for (unsigned int iPion=0; iPion<negPion.size(); ++iPion) {
	   unsigned int pionIndex = negPion.at(iPion);
	   TLorentzVector pionLV; pionLV.SetXYZT(genparticles_px[pionIndex],
						 genparticles_py[pionIndex],
						 genparticles_pz[pionIndex],
						 genparticles_e[pionIndex]);
	   pionPtH->Fill(pionLV.Pt(),weight);
	   float dRMuonPion = 1e+8;
	   for (unsigned int iMuon=0; iMuon<posMuon.size(); ++iMuon) {
	     unsigned int muonIndex = posMuon.at(iMuon);
	     TLorentzVector muonLV; muonLV.SetXYZT(genparticles_px[muonIndex],
						   genparticles_py[muonIndex],
						   genparticles_pz[muonIndex],
						   genparticles_e[muonIndex]);
	     float dRx = deltaR(pionLV.Eta(),pionLV.Phi(),
				muonLV.Eta(),muonLV.Phi());
	     if (dRx<dRMuonPion) dRMuonPion = dRx;
	   }
	   deltaRMuonPionH->Fill(dRMuonPion,weight);
	 }
       }
       if (posMuon.size()==0&&negMuon.size()==2&&negPion.size()==0&&posPion.size()==2) {
	 for (unsigned int iPion=0; iPion<posPion.size(); ++iPion) {
	   unsigned int pionIndex = posPion.at(iPion);
	   TLorentzVector pionLV; pionLV.SetXYZT(genparticles_px[pionIndex],
						 genparticles_py[pionIndex],
						 genparticles_pz[pionIndex],
						 genparticles_e[pionIndex]);
	   pionPtH->Fill(pionLV.Pt(),weight);
	   float dRMuonPion = 1e+8;
	   for (unsigned int iMuon=0; iMuon<negMuon.size(); ++iMuon) {
	     unsigned int muonIndex = negMuon.at(iMuon);
	     TLorentzVector muonLV; muonLV.SetXYZT(genparticles_px[muonIndex],
						   genparticles_py[muonIndex],
						   genparticles_pz[muonIndex],
						   genparticles_e[muonIndex]);
	     float dRx = deltaR(pionLV.Eta(),pionLV.Phi(),
				muonLV.Eta(),muonLV.Phi());
	     if (dRx<dRMuonPion) dRMuonPion = dRx;
	   }
	   deltaRMuonPionH->Fill(dRMuonPion,weight);
	 }
       }
     }
     if (HiggsFound) {
       TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[higgsIndex],
					       genparticles_py[higgsIndex],
					       genparticles_pz[higgsIndex],
					       genparticles_e[higgsIndex]);
       higgsPt = higgsLV.Pt();
       higgsTree->Fill();
       if (applyHiggsPtWeight) {
	 double HiggsPtForWeighting = higgsPt;
	 if (higgsPt>500) HiggsPtForWeighting = 499;
	 double higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
	 weight *= higgsPtWeight;
	 //	 std::cout << "Higgs : " << genparticles_pdgid[higgsIndex] 
	 //		   << "    pT = " << higgsLV.Pt()
	 //		   << "    eta = " << higgsLV.Eta() << " weight = " << higgsPtWeight << std::endl;

       }
     }

     //     histWeightsH->Fill(1.0,weight);
     counter_InputEventsH->Fill(1.0,weight);

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

     // checking if dimuon trigger bit is ON
     //     bool isDimuonTrigger = false;
     //     for (std::map<string,int>::iterator it=hltriggerresults->begin(); it!=hltriggerresults->end(); ++it) {
     //       TString trigName(it->first);
     //       if (trigName.Contains(DiMuonTriggerName)) {
     //  if (it->second==1)
     //    isDimuonTrigger = true;
     //       }
     //     }
     //     if (!isDimuonTrigger) continue;

     //     unsigned int ntrig = hltriggerresults->size();
     //     std::cout << "ntrig = " << ntrig << std::endl;
     //     for (std::map<string,int>::iterator it=hltriggerresults->begin(); it!=hltriggerresults->end(); ++it) 
       //       std::cout << it->first << "  :  "  << it->second << std::endl;
       //     std::cout << std::endl;

     //     unsigned int npres = hltriggerprescales->size();
     //     std::cout << "npres = " << npres << std::endl;
     //     for (std::map<string,int>::iterator it=hltriggerprescales->begin(); it!=hltriggerprescales->end(); ++it) 
     //       std::cout << it->first << "  :  "  << it->second << std::endl;
     //     std::cout << std::endl;
     
     // finding HLT filters in the HLT Filter library
     unsigned int nMu8Leg   = 0;
     unsigned int nMu17Leg  = 0;
     unsigned int nDZFilter = 0;
     unsigned int nSSFilter = 0;
     bool isMu8Leg = false;
     bool isMu17Leg = false;
     bool isDZFilter = false;
     bool isSSFilter = false;

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
       if (HLTFilter==DiMuonDzFilterName) {
	 nDZFilter = i;
	 isDZFilter = true;
       }
       if (HLTFilter==DiMuonSameSignFilterName) {
	 nSSFilter = i;
	 isSSFilter = true;
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
     if (!isDZFilter) {
       cout << "Filter " << DiMuonDzFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isSSFilter) {
       cout << "Filter " << DiMuonSameSignFilterName << " not found " << endl;
       exit(-1);
     }

     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     vector<bool> passSingleMu; passSingleMu.clear();
     vector<bool> passMu17; passMu17.clear();
     vector<bool> passMu8;  passMu8.clear();
     vector<bool> passMu12; passMu12.clear();
     vector<bool> passMu10; passMu10.clear();
     vector<bool> passMu5;  passMu5.clear();
     vector<bool> passMu45; passMu45.clear();
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
       rand.SetSeed((int)((muon_eta[i]+2.41)*100000));
       double rannum  = rand.Rndm();
       // old staff for trigger efficiency studies
       double effMu17 = MuLegEfficiency(muon_pt[i],muon_eta[i],17.0,2.4);
       double effMu8  = MuLegEfficiency(muon_pt[i],muon_eta[i],8.0,2.4);
       double effMu12 = MuLegEfficiency(muon_pt[i],muon_eta[i],12.0,2.4);
       double effMu10 = MuLegEfficiency(muon_pt[i],muon_eta[i],10.0,2.4);
       double effMu5  = MuLegEfficiency(muon_pt[i],muon_eta[i],5.0,2.4);
       double effMu45 = MuLegEfficiency(muon_pt[i],muon_eta[i],45.0,2.1);
       bool isMu17fired = effMu17 > rannum;
       bool isMu8fired  = effMu8  > rannum;
       bool isMu12fired = effMu12 > rannum;
       bool isMu10fired = effMu10 > rannum;
       bool isMu5fired  = effMu5  > rannum;
       bool isMu45fired = effMu45 > rannum;
       passMu17.push_back(isMu17fired);
       passMu12.push_back(isMu12fired);
       passMu10.push_back(isMu10fired);
       passMu8.push_back(isMu8fired);
       passMu5.push_back(isMu5fired);
       passMu45.push_back(isMu45fired);
       //
       muons.push_back(i);
     }
     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<2) continue; // quit event if number of good muons < 2
     counter_MuonSizeGTE2H->Fill(1.,weight);

     // ************************************************
     // ********** Some trigger studies ****************
     // *** TripleMu triggers vs. DoubleMuSS trigger ***
     // ************************************************

     // checking 4-muon final states
     if (muons.size()>=numberOfMuons) {
       bool isSingleMuon = false;
       for (unsigned int im=0; im<muons.size(); ++im) {
	 if (passMu45.at(im)) {
	   isSingleMuon = true;
	   break;
	 }
       }
       bool isDoubleMuSS = false;
       for (unsigned int im1=0; im1<muons.size()-1; ++im1) {
	 for (unsigned int im2=im1+1; im2<muons.size(); ++im2) {
	   int mu1 = muons.at(im1);
	   int mu2 = muons.at(im2);
	   if (muon_charge[mu1]*muon_charge[mu2]>0.0) {
	     if ( (passMu17.at(im1)&&passMu8.at(im2)) || (passMu17.at(im2)&&passMu8.at(im1))) {
	       isDoubleMuSS = true;
	       break;
	     }
	   }
	 }
	 if (isDoubleMuSS) break;
       }
       // 012
       bool triple012 = (passMu12.at(0)&&passMu10.at(1)&&passMu5.at(2));
       bool triple021 = (passMu12.at(0)&&passMu10.at(2)&&passMu5.at(1));
       bool triple102 = (passMu12.at(1)&&passMu10.at(0)&&passMu5.at(2));
       bool triple120 = (passMu12.at(1)&&passMu10.at(2)&&passMu5.at(0));
       bool triple210 = (passMu12.at(2)&&passMu10.at(1)&&passMu5.at(0));
       bool triple201 = (passMu12.at(2)&&passMu10.at(0)&&passMu5.at(1));
       bool Triple012 = triple012 || triple021 || triple102 || triple120 || triple210 || triple201;
       // 013
       bool Triple013 = false;
       if (muons.size()>3) {
	 bool triple013 = (passMu12.at(0)&&passMu10.at(1)&&passMu5.at(3));
	 bool triple031 = (passMu12.at(0)&&passMu10.at(3)&&passMu5.at(1));
	 bool triple103 = (passMu12.at(1)&&passMu10.at(0)&&passMu5.at(3));
	 bool triple130 = (passMu12.at(1)&&passMu10.at(3)&&passMu5.at(0));
	 bool triple310 = (passMu12.at(3)&&passMu10.at(1)&&passMu5.at(0));
	 bool triple301 = (passMu12.at(3)&&passMu10.at(0)&&passMu5.at(1));
	 Triple013 = triple013 || triple031 || triple103 || triple130 || triple310 || triple301;
       }
       // 123
       bool Triple123 = false;
       if (muons.size()>3) {
	 bool triple123 = (passMu12.at(1)&&passMu10.at(2)&&passMu5.at(3));
	 bool triple132 = (passMu12.at(1)&&passMu10.at(3)&&passMu5.at(2));
	 bool triple213 = (passMu12.at(2)&&passMu10.at(1)&&passMu5.at(3));
	 bool triple231 = (passMu12.at(2)&&passMu10.at(3)&&passMu5.at(1));
	 bool triple312 = (passMu12.at(3)&&passMu10.at(1)&&passMu5.at(2));
	 bool triple321 = (passMu12.at(3)&&passMu10.at(2)&&passMu5.at(1));
	 Triple123 = triple123 || triple132 || triple213 || triple231 || triple312 || triple321;
       }
       // TripleMu
       bool isTripleMuon = Triple012 || Triple013 || Triple123;
       if (isSingleMuon) histWeightsSingleMuH->Fill(1.0,weight);
       if (isDoubleMuSS) histWeightsDoubleMuSSH->Fill(1.0,weight);
       if (isTripleMuon) histWeightsTripleMuH->Fill(1.0,weight);
       bool isAllTriggers = isTripleMuon || isDoubleMuSS || isSingleMuon;
       if (isAllTriggers) histWeightsAllTriggersH->Fill(1.0,weight);

     }

     // *************************
     // selection of dimuon pairs 
     // *************************

     float maxPtSum = -1;
     int iLeading = -1;
     int iTrailing = -1;
     for (unsigned int i1=0; i1<muons.size()-1; ++i1) {
       int index1 = muons.at(i1);
       for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {
	 int index2 = muons.at(i2);
	 float ptSum = muon_pt[index1] + muon_pt[index2];
	 float charge = muon_charge[index1] *  muon_charge[index2];
	 bool chargeSelection = charge<0;
	 if (sameSign)
	   chargeSelection = charge>0;
	 if (!chargeSelection) continue;
	 float dRmuons = deltaR(muon_eta[index1],muon_phi[index1],
				muon_eta[index2],muon_phi[index2]);
	 
	 if (dRmuons<dRMuonsCut) continue;
	 bool mu1MatchMu17 = false;
	 bool mu1MatchMu8  = false;
	 bool mu1MatchDz   = false;
	 bool mu1MatchSS   = false;
	 for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	   float dRtrig = deltaR(muon_eta[index1],muon_phi[index1],
				 trigobject_eta[iT],trigobject_phi[iT]);
	   if (dRtrig>DRTrigMatch) continue;
	   if (trigobject_filters[iT][nMu17Leg])
	     mu1MatchMu17 = true;
	   if (trigobject_filters[iT][nMu8Leg])
	     mu1MatchMu8 = true;
	   if (trigobject_filters[iT][nDZFilter])
	     mu1MatchDz = true;
	   if (trigobject_filters[iT][nSSFilter])
	     mu1MatchSS = true;
	 }
	 bool mu2MatchMu17 = false;
	 bool mu2MatchMu8  = false;
	 bool mu2MatchDz   = false;
	 bool mu2MatchSS   = false;
	 for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	   float dRtrig = deltaR(muon_eta[index2],muon_phi[index2],
				 trigobject_eta[iT],trigobject_phi[iT]);
	   if (dRtrig>DRTrigMatch) continue;
	   if (trigobject_filters[iT][nMu17Leg])
	     mu2MatchMu17 = true;
	   if (trigobject_filters[iT][nMu8Leg])
	     mu2MatchMu8 = true;
	   if (trigobject_filters[iT][nDZFilter])
	     mu2MatchDz = true;
	   if (trigobject_filters[iT][nSSFilter])
	     mu2MatchSS = true;
	 }

	 bool mu1PtHigh = muon_pt[index1]>ptMuonHighCut && fabs(muon_eta[index1])<etaMuonHighCut;
	 bool mu1PtLow  = muon_pt[index1]>ptMuonLowCut && fabs(muon_eta[index1])<etaMuonLowCut;
	 bool mu2PtHigh = muon_pt[index2]>ptMuonHighCut && fabs(muon_eta[index2])<etaMuonHighCut;
	 bool mu2PtLow  = muon_pt[index2]>ptMuonLowCut && fabs(muon_eta[index2])<etaMuonLowCut;
	 
	 // trigger condition
	 bool isTriggerMatched = true;
	 if (applyTriggerMatch) {
	   isTriggerMatched = (mu1MatchMu17&&mu2MatchMu8&&mu1PtHigh&&mu2PtLow)||(mu1MatchMu8&&mu2MatchMu17&&mu1PtLow&&mu2PtHigh);
	   if (isData) {
	     isTriggerMatched = isTriggerMatched && mu1MatchSS && mu2MatchSS;
	     if (event_run<=274442||event_run>=280919) // when dZ filter is present
	       isTriggerMatched = isTriggerMatched && mu1MatchDz && mu2MatchDz;
	   }
	 }
	 else {
	   isTriggerMatched = (mu1PtHigh&&mu2PtLow) || (mu1PtLow&&mu2PtHigh);
	 }
	 if (!isTriggerMatched) continue;
	 if (ptSum>maxPtSum) { // choose the mair with maximum Sum(pT)
	   maxPtSum = ptSum;
	   if (muon_pt[index1]>muon_pt[index2]) {
	     iLeading = index1;
	     iTrailing = index2;
	   }
	   else {
	     iLeading = index2;
	     iTrailing = index1;
	   }
	 }
       }
     }


     if (iLeading<0) continue;
     if (iTrailing<0) continue;

     double triggerWeight = 1;
     double idLeadingWeight = 1;
     double idTrailingWeight = 1;
     if (!isData) { // trigger efficiency here
       double ptLeading = muon_pt[iLeading];
       double etaLeading = muon_eta[iLeading];
       double ptTrailing = muon_pt[iTrailing];
       double etaTrailing = muon_eta[iTrailing];

       if (ptLeading<10.0) ptLeading = 10.01;
       if (ptTrailing<10.0) ptTrailing = 10.01;

       correctionWS->var("m_pt")->setVal(ptLeading);
       correctionWS->var("m_eta")->setVal(etaLeading);
       double idLeadingW  = correctionWS->function("m_id_ratio")->getVal();
       double trkLeadingW = correctionWS->function("m_trk_ratio")->getVal();
       idLeadingWeight = idLeadingW * trkLeadingW;
       correctionWS->var("m_pt")->setVal(ptTrailing);
       correctionWS->var("m_eta")->setVal(etaTrailing);
       double idTrailingW = correctionWS->function("m_id_ratio")->getVal();
       double trkTrailingW = correctionWS->function("m_trk_ratio")->getVal();
       idTrailingWeight = idTrailingW * trkTrailingW;

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

       if (applyTriggerMatch) {
	 if (trigWeightMC>0)
	   triggerWeight = trigWeightData/trigWeightMC;
       }
       else
	 triggerWeight = trigWeightData;

       triggerWeight *= effDzSS;

       /*       
       std::cout << "pT(leading) = " << ptLeading
		 << "   eta(leading) = " << etaLeading
		 << "   pT(trailing) = " << ptTrailing
		 << "   eta(trailing) = " << etaTrailing << std::endl;
       std::cout << "IdW(leading) = " << idLeadingW
		 << "   TrkW(leading) = " << trkLeadingW 
		 << "   IdW(trailing) = " << idTrailingW
		 << "   TrkW(trailing) = " << trkTrailingW << std::endl;
       std::cout << "Trigger weight = " << triggerWeight << std::endl;
       */

     }
     triggerWeightH->Fill(triggerWeight,1.0);
     weight *= triggerWeight;
     weight *= idLeadingWeight;
     weight *= idTrailingWeight;

     // dimuon selection passed 
     TLorentzVector LeadingMuon4; LeadingMuon4.SetXYZM(muon_px[iLeading],
						       muon_py[iLeading],
						       muon_pz[iLeading],
						       MuMass);
     
     TLorentzVector TrailingMuon4; TrailingMuon4.SetXYZM(muon_px[iTrailing],
							 muon_py[iTrailing],
							 muon_pz[iTrailing],
							 MuMass);
     
     TLorentzVector diMuon4 = LeadingMuon4 + TrailingMuon4;
     
     float dimuonMass = diMuon4.M();

     // trigger weight (obsolete)
     //     if (!isData) {
     //       float effMu17Leading  = MuLegEfficiency(LeadingMuon4.Pt(),LeadingMuon4.Eta(),17.0,2.4);
     //       float effMu8Leading   = MuLegEfficiency(LeadingMuon4.Pt(),LeadingMuon4.Eta(), 8.0,2.4);
     //       float effMu17Trailing = MuLegEfficiency(TrailingMuon4.Pt(),TrailingMuon4.Eta(),17.0,2.4);
     //       float effMu8Trailing  = MuLegEfficiency(TrailingMuon4.Pt(),TrailingMuon4.Eta(), 8.0,2.4);
     //       float trigWeight = 0.935*(effMu17Leading*effMu8Trailing+effMu8Leading*effMu17Trailing-effMu17Leading*effMu17Trailing);
     //       std::cout << "trig weight = " << trigWeight << std::endl;
     //       weight *= trigWeight; 
     //     }

     // filling histograms (muon kinematics)
     dimuonMassH->Fill(dimuonMass,weight);
     ptLeadingMuH->Fill(muon_pt[iLeading],weight);
     ptTrailingMuH->Fill(muon_pt[iTrailing],weight);
     etaLeadingMuH->Fill(muon_eta[iLeading],weight);
     etaTrailingMuH->Fill(muon_eta[iTrailing],weight);
     counter_MuonKinematicsH->Fill(1.0,weight);                  

     TLorentzVector Met4; Met4.SetXYZM(metx,mety,0,0);


     // counting tracks around each muon
     std::vector<unsigned int> trkLeadingMu; trkLeadingMu.clear(); // all tracks
     std::vector<unsigned int> trkTrailingMu; trkTrailingMu.clear(); // all tracks
     std::vector<unsigned int> trkSigLeadingMu; trkSigLeadingMu.clear(); // signal tracks
     std::vector<unsigned int> trkSigTrailingMu; trkSigTrailingMu.clear(); // signal tracks
     std::vector<double> trkSigDRLeadingMu; trkSigDRLeadingMu.clear();
     std::vector<double> trkSigDRTrailingMu; trkSigDRTrailingMu.clear();
     unsigned int hardestTrkLeading = 0; // index of hardest track around leading mu
     unsigned int hardestTrkTrailing = 0; // index of hardest track around trailing mu
     unsigned int softestTrkLeading = 0; // index of softest track around leading mu
     unsigned int softestTrkTrailing = 0; // index of softest track around trailing mu
     

     float ptHardestLeading = 0;
     float ptHardestTrailing = 0;
     float ptSoftestLeading = 1e+10;
     float ptSoftestTrailing = 1e+10;
     std::vector<unsigned int> Soft_trkLeadingMu; Soft_trkLeadingMu.clear(); // Lead Muon tracks control region
     std::vector<unsigned int> Soft_trkTrailingMu; Soft_trkTrailingMu.clear(); // Trail Muon tracks control region

     
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
       
       TLorentzVector leadingMuDiff = LeadingMuon4 - trk4;
       if (leadingMuDiff.P()>0.1) { // track is not leading muon
	 float drTrkMu = deltaR(muon_eta[iLeading],muon_phi[iLeading],
				track_eta[iTrk],   track_phi[iTrk]);
	 float qTrkLeadingMu = track_charge[iTrk]*muon_charge[iLeading];
	 if (drTrkMu<dRIsoMuon){
	   trkLeadingMu.push_back(iTrk);
	   if (track_pt[iTrk]>ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	     Soft_trkLeadingMu.push_back(iTrk);
	 }
	 if (drTrkMu<dRIsoMuon && qTrkLeadingMu<0 && fabs(track_dxy[iTrk])<dxyTrkCut && fabs(track_dz[iTrk])<dzTrkCut && track_pt[iTrk]>ptTrkCut) {
	   trkSigLeadingMu.push_back(iTrk);
	   trkSigDRLeadingMu.push_back(drTrkMu);
	   if (track_pt[iTrk]>ptHardestLeading) {
	     ptHardestLeading = track_pt[iTrk];
	     hardestTrkLeading = iTrk;
	   }
	   if (track_pt[iTrk]<ptSoftestLeading) {
	     ptSoftestLeading = track_pt[iTrk];
	     softestTrkLeading = iTrk;
	   }
	 }
       }
       
       TLorentzVector trailingMuDiff = TrailingMuon4 - trk4;
       if (trailingMuDiff.P()>0.1) { // track is not trailing muon
	 float drTrkMu = deltaR(muon_eta[iTrailing],muon_phi[iTrailing],
				track_eta[iTrk],track_phi[iTrk]);
	 float qTrkTrailingMu = track_charge[iTrk]*muon_charge[iTrailing];         
	 if (drTrkMu<dRIsoMuon){
	   trkTrailingMu.push_back(iTrk);
	   if (track_pt[iTrk] > ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	     Soft_trkTrailingMu.push_back(iTrk);
	 }
	 if (drTrkMu<dRIsoMuon && qTrkTrailingMu<0 && fabs(track_dxy[iTrk])<dxyTrkCut && fabs(track_dz[iTrk])<dzTrkCut && track_pt[iTrk]>ptTrkCut) {
	   trkSigTrailingMu.push_back(iTrk);
	   trkSigDRTrailingMu.push_back(drTrkMu);
	   if (track_pt[iTrk]>ptHardestTrailing) {
	     ptHardestTrailing = track_pt[iTrk];
	     hardestTrkTrailing = iTrk;
	   }
	   if (track_pt[iTrk]<ptSoftestTrailing) {
	     ptSoftestTrailing = track_pt[iTrk];
	     softestTrkTrailing = iTrk;
	   }
	   
	 }
       }
       
     }
     
     nTracksLeadingMuH->Fill(float(trkLeadingMu.size()),weight);
     nTracksTrailingMuH->Fill(float(trkTrailingMu.size()),weight);
     nSigTracksLeadingMuH->Fill(float(trkSigLeadingMu.size()),weight);
     nSigTracksTrailingMuH->Fill(float(trkSigTrailingMu.size()),weight);
     nSoftTracksLeadingMuH->Fill(float(Soft_trkLeadingMu.size()),weight);
     nSoftTracksTrailingMuH->Fill(float(Soft_trkTrailingMu.size()),weight);

     float chargedIsoLeading = muon_chargedHadIso[iLeading];
     for (unsigned int iTrk=0; iTrk<trkSigLeadingMu.size(); ++iTrk) {
       unsigned int indexTrack = trkSigLeadingMu.at(iTrk);       
       if (trkSigDRLeadingMu.at(iTrk)<0.4) {
	 TLorentzVector trk4; trk4.SetXYZM(track_px[indexTrack],
					   track_py[indexTrack],
					   track_pz[indexTrack],
					   track_mass[indexTrack]);
	 chargedIsoLeading -= trk4.Pt();
       }
     }
     if (chargedIsoLeading<0.0) chargedIsoLeading = 0;
     float neutralIsoLeading = muon_neutralHadIso[iLeading] + muon_photonIso[iLeading] - 0.5*muon_puIso[iLeading];
     if (neutralIsoLeading<0) neutralIsoLeading = 0;
     float isoLeading = (chargedIsoLeading+neutralIsoLeading)/muon_pt[iLeading];
     
     float chargedIsoTrailing = 0;
     for (unsigned int iTrk=0; iTrk<trkSigTrailingMu.size(); ++iTrk) {
       unsigned int indexTrack = trkSigTrailingMu.at(iTrk);       
       if (trkSigDRTrailingMu.at(iTrk)<0.4) {
	 TLorentzVector trk4; trk4.SetXYZM(track_px[indexTrack],
					   track_py[indexTrack],
					   track_pz[indexTrack],
					   track_mass[indexTrack]);
	 chargedIsoTrailing -= trk4.Pt();
       }
     }
     if (chargedIsoTrailing<0.0) chargedIsoTrailing = 0;
     float neutralIsoTrailing = muon_neutralHadIso[iTrailing] + muon_photonIso[iTrailing] - 0.5*muon_puIso[iTrailing];
     if (neutralIsoTrailing<0) neutralIsoTrailing = 0;
     float isoTrailing = (chargedIsoTrailing+neutralIsoTrailing)/muon_pt[iTrailing];
     
     // definition of signal muon+track
     bool signalLeadingMu  = trkLeadingMu.size()==1 && trkSigLeadingMu.size()==1;
     bool signalTrailingMu = trkTrailingMu.size()==1 && trkSigTrailingMu.size()==1;

     double weightTrkLeading = 1;
     double weightTrkTrailing = 1;
     if (trkSigLeadingMu.size()>0&&!isData) {
       unsigned int iTrkLeading = trkSigLeadingMu.at(0);
       int absPdgId = TMath::Abs(track_ID[iTrkLeading]);
       if (absPdgId==11) {
	 correctionWS->var("e_pt")->setVal(track_pt[iTrkLeading]);
	 correctionWS->var("e_eta")->setVal(track_eta[iTrkLeading]);
	 weightTrkLeading = correctionWS->function("e_trk_ratio")->getVal();
       }
       else {
	 correctionWS->var("m_pt")->setVal(track_pt[iTrkLeading]);
         correctionWS->var("m_eta")->setVal(track_eta[iTrkLeading]);
         weightTrkLeading = correctionWS->function("m_trk_ratio")->getVal();
       }
       //       std::cout << "track around leading : pT = " << track_pt[iTrkLeading]
       //		 << "  eta = " << track_eta[iTrkLeading] 
       //		 << "  Id = " << track_ID[iTrkLeading]
       //		 << "  weight(trk) = " << weightTrkLeading << std::endl;
	
     }
     weight *= weightTrkLeading;

     if (trkSigTrailingMu.size()>0&&!isData) {
       unsigned int iTrkTrailing = trkSigTrailingMu.at(0);
       int absPdgId = TMath::Abs(track_ID[iTrkTrailing]);
       if (absPdgId==11) {
         correctionWS->var("e_pt")->setVal(track_pt[iTrkTrailing]);
         correctionWS->var("e_eta")->setVal(track_eta[iTrkTrailing]);
         weightTrkTrailing = correctionWS->function("e_trk_ratio")->getVal();
       }
       else {
         correctionWS->var("m_pt")->setVal(track_pt[iTrkTrailing]);
         correctionWS->var("m_eta")->setVal(track_eta[iTrkTrailing]);
         weightTrkTrailing = correctionWS->function("m_trk_ratio")->getVal();
       }
       //       std::cout << "track around trailing : pT = " << track_pt[iTrkTrailing]
       //                 << "  eta = " << track_eta[iTrkTrailing]
       //                 << "  Id = " << track_ID[iTrkTrailing]
       //                 << "  weight(trk) = " << weightTrkTrailing << std::endl;
     }
     //     std::cout << std::endl;
     weight *= weightTrkTrailing;


     if (trkLeadingMu.size()==1) {

       unsigned int iTrkLeading = trkLeadingMu.at(0);

       ptTrackLeadingMuH->Fill(track_pt[iTrkLeading],weight);
       etaTrackLeadingMuH->Fill(track_eta[iTrkLeading],weight);
       dxyTrackLeadingMuH->Fill(track_dxy[iTrkLeading],weight);
       dzTrackLeadingMuH->Fill(track_dz[iTrkLeading],weight);

       ptTrackN1H->Fill(track_pt[iTrkLeading],weight);
       etaTrackN1H->Fill(track_eta[iTrkLeading],weight);
       dxyTrackN1H->Fill(track_dxy[iTrkLeading],weight);
       dzTrackN1H->Fill(track_dz[iTrkLeading],weight);
       
     }

     if (trkTrailingMu.size()==1) {

       unsigned int iTrkTrailing = trkTrailingMu.at(0);

       ptTrackTrailingMuH->Fill(track_pt[iTrkTrailing],weight);
       etaTrackTrailingMuH->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackTrailingMuH->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackTrailingMuH->Fill(track_dz[iTrkTrailing],weight);

       ptTrackN1H->Fill(track_pt[iTrkTrailing],weight);
       etaTrackN1H->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackN1H->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackN1H->Fill(track_dz[iTrkTrailing],weight);
       
     }

     if (signalLeadingMu) {
       unsigned int iTrkLeading = trkLeadingMu.at(0);

       ptTrackH->Fill(track_pt[iTrkLeading],weight);
       etaTrackH->Fill(track_eta[iTrkLeading],weight);
       dxyTrackH->Fill(track_dxy[iTrkLeading],weight);
       dzTrackH->Fill(track_dz[iTrkLeading],weight);

       counter_nMuTrackSigH->Fill(1.0,weight);                 

     }

     if (signalTrailingMu) {
       unsigned int iTrkTrailing = trkTrailingMu.at(0);
       
       ptTrackH->Fill(track_pt[iTrkTrailing],weight);
       etaTrackH->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackH->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackH->Fill(track_dz[iTrkTrailing],weight);

       counter_nMuTrackSigH->Fill(1.0,weight);                 

     }
     
     // sideband N23 and N45
     bool isN23leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
     bool isN45leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==4||trkTrailingMu.size()==5);
     bool isN23trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==2||trkLeadingMu.size()==3); 
     bool isN45trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==4||trkLeadingMu.size()==5); 

     
     // sidebands Ntrk23
     bool isNtrk23leading  = trkSigLeadingMu.size()>0 && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
     bool isNtrk23trailing = trkSigTrailingMu.size()>0 && (trkLeadingMu.size()==2||trkLeadingMu.size()==3);
     
     // sidebands Ntrk1
     bool isNtrk1leading  = trkSigLeadingMu.size()>0 && trkTrailingMu.size()==1;
     bool isNtrk1trailing = trkSigTrailingMu.size()>0 && trkLeadingMu.size()==1;
     
     // signal region
     bool signalRegion = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkSigTrailingMu.size()==1&&trkTrailingMu.size()==1);
     
     // sidebands N23 (see definition of this sideband in HIG-14-019)
     if (isN23leading) {
       int iTrk = trkSigLeadingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
       float deltaRMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN23H->Fill(mass,weight);
     }


     if (isN45leading) {
       int iTrk = trkSigLeadingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
       float deltaRMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN45H->Fill(mass,weight);
     }
     
     if (isN23trailing) {
       int iTrk = trkSigTrailingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
       float deltaRMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN23H->Fill(mass,weight);
     }      
     
     if (isN45trailing) {
       int iTrk = trkSigTrailingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
       float deltaRMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN45H->Fill(mass,weight);
     }      
     

     // sidebands Ntrk23 (see definition of this sideband in HIG-14-019)
     if (isNtrk23leading) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkLeading],
							   track_py[hardestTrkLeading],
							   track_pz[hardestTrkLeading],
							   track_mass[hardestTrkLeading]);
       TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk23leadingH->Fill(mass,weight);
       InvMassHardestNtrk23H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkLeading],
							   track_py[softestTrkLeading],
							   track_pz[softestTrkLeading],
							   track_mass[softestTrkLeading]);
       TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk23leadingH->Fill(mass,weight);
       InvMassSoftestNtrk23H->Fill(mass,weight);
     }      
     if (isNtrk23trailing) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkTrailing],
							   track_py[hardestTrkTrailing],
							   track_pz[hardestTrkTrailing],
							   track_mass[hardestTrkTrailing]);
       TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk23trailingH->Fill(mass,weight);
       InvMassHardestNtrk23H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkTrailing],
							   track_py[softestTrkTrailing],
							   track_pz[softestTrkTrailing],
							   track_mass[softestTrkTrailing]);
       TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk23trailingH->Fill(mass,weight);
       InvMassSoftestNtrk23H->Fill(mass,weight);
     }      
      
     // sidebands Ntrk1 (see definition of this sideband in HIG-14-019)
     if (isNtrk1leading) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkLeading],
							   track_py[hardestTrkLeading],
							   track_pz[hardestTrkLeading],
							   track_mass[hardestTrkLeading]);
       TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk1leadingH->Fill(mass,weight);
       InvMassHardestNtrk1H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkLeading],
							   track_py[softestTrkLeading],
							   track_pz[softestTrkLeading],
							   track_mass[softestTrkLeading]);
       TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk1leadingH->Fill(mass,weight);
       InvMassSoftestNtrk1H->Fill(mass,weight);
     }      
     if (isNtrk1trailing) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkTrailing],
							   track_py[hardestTrkTrailing],
							   track_pz[hardestTrkTrailing],
							   track_mass[hardestTrkTrailing]);
       TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk1trailingH->Fill(mass,weight);
       InvMassHardestNtrk1H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkTrailing],
							   track_py[softestTrkTrailing],
							   track_pz[softestTrkTrailing],
							   track_mass[softestTrkTrailing]);
       TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk1trailingH->Fill(mass,weight);
       InvMassSoftestNtrk1H->Fill(mass,weight);
     }      

     // **************
     // signal region
     // *************
     if (signalRegion) {
       counter_FinalEventsH->Fill(1.,weight);
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector MuonTrackLeading4 = LeadingMuon4 + TrackLeading4;
       float massTrkMuLeading = MuonTrackLeading4.M();
       
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector MuonTrackTrailing4 = TrailingMuon4 + TrackTrailing4;
       float massTrkMuTrailing = MuonTrackTrailing4.M();
       
       TLorentzVector Visible4 = MuonTrackTrailing4 + MuonTrackLeading4;
       TLorentzVector Total4 = Visible4 + Met4;

       float dRLeadingMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				     TrackLeading4.Eta(),TrackLeading4.Phi());


       float dRTrailingMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				      TrackTrailing4.Eta(),TrackTrailing4.Phi());



       MetSelH->Fill(Met4.Pt(),weight);
       mTtotSelH->Fill(Total4.M(),weight);
       dimuonMassSelH->Fill(diMuon4.M(),weight);
       invMass2Mu2TrkSelH->Fill(Visible4.M(),weight);

       ptLeadingMuSelH->Fill(LeadingMuon4.Pt(),weight);
       ptTrailingMuSelH->Fill(TrailingMuon4.Pt(),weight);
       ptLeadingMuTrkSelH->Fill(TrackLeading4.Pt(),weight);
       ptTrailingMuTrkSelH->Fill(TrackTrailing4.Pt(),weight);

       ptMuSelH->Fill(LeadingMuon4.Pt(),weight);
       ptMuSelH->Fill(TrailingMuon4.Pt(),weight);
       ptTrkSelH->Fill(TrackLeading4.Pt(),weight);
       ptTrkSelH->Fill(TrackTrailing4.Pt(),weight);

       etaLeadingMuSelH->Fill(LeadingMuon4.Eta(),weight);
       etaTrailingMuSelH->Fill(TrailingMuon4.Eta(),weight);
       etaLeadingMuTrkSelH->Fill(TrackLeading4.Eta(),weight);
       etaTrailingMuTrkSelH->Fill(TrackTrailing4.Eta(),weight);

       etaMuSelH->Fill(LeadingMuon4.Eta(),weight);
       etaMuSelH->Fill(TrailingMuon4.Eta(),weight);
       
       etaTrkSelH->Fill(TrackLeading4.Eta(),weight);
       etaTrkSelH->Fill(TrackTrailing4.Eta(),weight);

       dRMuTrkSelH->Fill(dRLeadingMuTrk,weight);
       dRMuTrkSelH->Fill(dRTrailingMuTrk,weight);

       InvMassLeadingH->Fill(massTrkMuLeading,weight);
       InvMassTrailingH->Fill(massTrkMuTrailing,weight);
       InvMassH->Fill(massTrkMuLeading,weight);
       InvMassH->Fill(massTrkMuTrailing,weight);
       InvMass2DH->Fill(massTrkMuLeading,massTrkMuTrailing,weight);
       
     }

     // ******************************************************************
     // Control region to study mass correlations (RegionA in HIG-14-019)
     // *****************************************************************
     
     bool bkgdLeadingMu = 
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3);
     
     bool bkgdTrailingMu = 
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3);

     bool bkgdXLeadingMu = 
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==3 && trkLeadingMu.size()==4);
     
     bool bkgdXTrailingMu = 
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==3 && trkTrailingMu.size()==4);
     
     bool ControlAll = (signalLeadingMu&&bkgdTrailingMu) || 
       (signalTrailingMu&&bkgdLeadingMu) || (bkgdLeadingMu&&bkgdTrailingMu);
     
     // * Now we can use this boolean to select bkg sideband
     // * where correlation coefficients are computed
     // * It is sufficient to use only one boolean - ControlAll
     if(ControlAll){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();
       
       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlH->Fill(masslow, masshigh, weight);

       counter_ControlEventsH->Fill(1.0,weight);

     }

     // ********** ControlX ****************************
     bool ControlAllX = bkgdXLeadingMu&&bkgdXTrailingMu;
     if(ControlAllX){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();
       
       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlXH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlXH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlXH->Fill(masslow, masshigh, weight);

       counter_ControlXEventsH->Fill(1.0,weight);

     }
     
     // ********* ControlY *********************
     bool ControlAllY = (signalLeadingMu&&bkgdXTrailingMu) ||
       (signalTrailingMu&&bkgdXLeadingMu) || (bkgdXLeadingMu&&bkgdXTrailingMu);
     if(ControlAllY){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();
       
       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlYH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlYH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlYH->Fill(masslow, masshigh, weight);

       counter_ControlYEventsH->Fill(1.0,weight);

     }
     
   } // icand loop
   
   delete tree_;
   file_->Close();
   delete file_;
   
  }// filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
}// int main loop 

 

