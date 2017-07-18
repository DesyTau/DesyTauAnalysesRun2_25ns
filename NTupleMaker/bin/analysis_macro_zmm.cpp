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
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "TSystem.h"

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

bool looseId(
	     float energy,
	     float eta,
	     float energycorr,
	     float chargedhadronicenergy,
	     float neutralhadronicenergy,
	     float neutralemenergy,
	     float chargedemenergy,
	     float muonenergy,
	     unsigned int chargedmulti,
	     unsigned int neutralmulti
	     ) {
  bool looseJetID = false;
  energy *= energycorr; // uncorrected energy must be used                                                                                                             
  float chf = chargedhadronicenergy/energy;
  float nhf = neutralhadronicenergy/energy;
  float nem = neutralemenergy/energy;
  float elf = chargedemenergy/energy;
  float muf = muonenergy/energy;
  float chm = chargedmulti;
  float nm  = neutralmulti;
  float npr = chargedmulti + neutralmulti;

  if (fabs(eta)<=2.7)
    {
      looseJetID = (nhf < 0.99 && nem < 0.99 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0 && elf < 0.99));
    }
  else if (fabs(eta)<=3.0)
    {
      looseJetID = (nem < 0.9 && nm > 2);
    }
  else
    {
      looseJetID = nem < 0.9 && nm > 10;
    }

  return looseJetID;


}

float topPtWeight(float pt1,
                  float pt2,
		  bool run1) {
    
  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter

  if (run1) {
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
    a = 0.156;    // Run1 a parameter
    b = -0.00137;  // Run1 b parameter
  }
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
  
  return TMath::Sqrt(w1*w2);
    
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
  const bool isTOP  = cfg.get<bool>("IsTOP");
  const bool isDY   = cfg.get<bool>("IsDY");
  const bool isW    = cfg.get<bool>("IsW");
  const bool isZMM  = cfg.get<bool>("IsZMM");

  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");

  
  // kinematic cuts on muons
  const float ptMuonCut      = cfg.get<float>("ptMuonLowCut");
  const float ptMuonLooseCut = cfg.get<float>("ptMuonLooseCut");
  const float etaMuonCut  = cfg.get<float>("etaMuonHighCut");
  const float dxyMuonCut  = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut   = cfg.get<float>("dzMuonCut");
  const float isoMuonCut  = cfg.get<float>("isoMuonCut");
  const float isoMuonLooseCut  = cfg.get<float>("isoMuonLooseCut");

  // topological cuts
  const float dPhiMuonTrkCut   = cfg.get<float>("dPhiMuonTrkCut");

  // track selection
  const float dRIsoMuon       = cfg.get<float>("dRIsoMuon");
  const bool sameSign      = cfg.get<bool>("SameSignMuons");

  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  // jets
  const float jetEtaCut = cfg.get<float>("jetEtaCut"); 
  const float jetPtCut  = cfg.get<float>("jetPtCut");
  const float deltaRJetLeptonCut = cfg.get<float>("deltaRJetLeptonCut");

  // trigger
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const string singleMuonTriggerName1 = cfg.get<string>("SingleMuonTriggerName1");
  const string singleMuonFilterName1 = cfg.get<string>("SingleMuonFilterName1");
  const string singleMuonTriggerName2 = cfg.get<string>("SingleMuonTriggerName2");
  const string singleMuonFilterName2 = cfg.get<string>("SingleMuonFilterName2");

  TString SingleMuonTriggerName1(singleMuonTriggerName1);
  TString SingleMuonFilterName1(singleMuonFilterName1);
  TString SingleMuonTriggerName2(singleMuonTriggerName2);
  TString SingleMuonFilterName2(singleMuonFilterName2);

  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  const string MuonTriggerEffFile = cfg.get<string>("MuonTriggerEffFile");
  const string MuonIdEffFile     = cfg.get<string>("MuonIdEffFile");

  const string recoilFileName = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

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

  Int_t genparticles_fromHardProcess[1000];
  Int_t genparticles_fromHardProcessBeforeFSR[1000];
  Int_t genparticles_isDecayedLeptonHadron[1000];
  Int_t genparticles_isDirectHadronDecayProduct[1000];
  Int_t genparticles_isDirectHardProcessTauDecayProduct[1000];
  Int_t genparticles_isDirectPromptTauDecayProduct[1000];
  Int_t genparticles_isDirectTauDecayProduct[1000];
  Int_t genparticles_isFirstCopy[1000];
  Int_t genparticles_isHardProcess[1000];
  Int_t genparticles_isHardProcessTauDecayProduct[1000];
  Int_t genparticles_isLastCopy[1000];
  Int_t genparticles_isLastCopyBeforeFSR[1000];
  Int_t genparticles_isPrompt[1000];
  Int_t genparticles_isPromptTauDecayProduct[1000];
  Int_t genparticles_isTauDecayProduct[1000];

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

  // pat jets                                                                                                                                                                                                      
  UInt_t pfjet_count;
  Float_t pfjet_e[1000];
  Float_t pfjet_px[1000];
  Float_t pfjet_py[1000];
  Float_t pfjet_pz[1000];
  Float_t pfjet_pt[1000];
  Float_t pfjet_eta[1000];
  Float_t pfjet_phi[1000];
  Float_t pfjet_energycorr[1000];

  Float_t pfjet_neutralhadronicenergy[1000];
  Float_t pfjet_chargedhadronicenergy[1000];
  Float_t pfjet_neutralemenergy[1000];
  Float_t pfjet_chargedemenergy[1000];
  Float_t pfjet_muonenergy[1000];
  Float_t pfjet_chargedmuonenergy[1000];
  UInt_t pfjet_chargedmulti[1000];
  UInt_t pfjet_neutralmulti[1000];
  UInt_t pfjet_chargedhadronmulti[1000];

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
  
  // weights
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * muTrigWeightH = new TH1D("muTrigWeightH","",100,0,2);
  TH1D * muIdWeightH = new TH1D("muIdWeightH","",100,0,2);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);
  //TH1D * histWeightsH = (TH1D*)file->Get("makeroottree/histWeightsH");

  // histograms after final selection
  TH1D * ptMuonH         = new TH1D("ptMuonH","",100,0,500);
  TH1D * etaMuonH = new TH1D("etaMuonH","",25,-2.5,2.5);
  TH1D * ptTrackH = new TH1D("ptTrackH","",400,0,400);
  TH1D * ptTrackLowH = new TH1D("ptTrackLowH","",100,0,50);
  TH1D * etaTrackH = new TH1D("etaTrackH","",25,-2.5,2.5);
  TH1D * mTH         = new TH1D("mTH","",500,0,500);
  TH1D * invMassMu1Mu2H         = new TH1D("invMassMu1Mu2H","",500,0,500);
  TH1D * metH         = new TH1D("metH","",500,0,500);
  TH1D * metH_2p5to5  = new TH1D("metH_2p5to5","",500,0,500);
  TH1D * metH_5to7p5  = new TH1D("metH_5to7p5","",500,0,500);
  TH1D * metH_7p5to10 = new TH1D("metH_7p5to10","",500,0,500);
  TH1D * metH_10to15  = new TH1D("metH_10to15","",500,0,500);
  TH1D * metH_15to20  = new TH1D("metH_15to20","",500,0,500);
  TH1D * metH_20toInf = new TH1D("metH_20toInf","",500,0,500);
  TH1D * dPhiMuonMetH = new TH1D("dPhiMuonMetH","",32,0,3.15);
  TH1D * dPhiMuonJetH = new TH1D("dPhiMuonJetH","",32,0,3.15);
  TH1D * nJetH = new TH1D("nJetH","",10,0,10);

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load( pathToCrystalLib );
  if (openSuccessful !=0 ) {
    cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit( -1 );
  }

  //Lepton Scale Factors from KIT - from synch ntuple producer
  TString filename_kitLeptonSF = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_5.root";
  TFile f_kitLeptonSF(filename_kitLeptonSF);
  if (f_kitLeptonSF.IsZombie()) {std::cout << " workspace file " << filename_kitLeptonSF << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w_kitLeptonSF = (RooWorkspace*)f_kitLeptonSF.Get("w");
  f_kitLeptonSF.Close();

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
  ScaleFactor * SF_muonId = new ScaleFactor();
  SF_muonId->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdEffFile));
  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerEffFile));

  // Z-pt weights from correction WS                    
  TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_5.root";
  TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  // Recoil correctot
  RecoilCorrector recoilMetCorrector(RecoilFileName);

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
     tree_->SetBranchAddress("genparticles_fromHardProcess", genparticles_fromHardProcess);
     tree_->SetBranchAddress("genparticles_isDirectHardProcessTauDecayProduct", genparticles_isDirectHardProcessTauDecayProduct);
   }   

   tree_->SetBranchAddress("pfjet_count",&pfjet_count);
   tree_->SetBranchAddress("pfjet_e",pfjet_e);
   tree_->SetBranchAddress("pfjet_px",pfjet_px);
   tree_->SetBranchAddress("pfjet_py",pfjet_py);
   tree_->SetBranchAddress("pfjet_pz",pfjet_pz);
   tree_->SetBranchAddress("pfjet_pt",pfjet_pt);
   tree_->SetBranchAddress("pfjet_eta",pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi",pfjet_phi);
   tree_->SetBranchAddress("pfjet_energycorr",pfjet_energycorr);
   tree_->SetBranchAddress("pfjet_neutralhadronicenergy",pfjet_neutralhadronicenergy);
   tree_->SetBranchAddress("pfjet_chargedhadronicenergy",pfjet_chargedhadronicenergy);
   tree_->SetBranchAddress("pfjet_neutralemenergy",pfjet_neutralemenergy);
   tree_->SetBranchAddress("pfjet_chargedemenergy",pfjet_chargedemenergy);
   tree_->SetBranchAddress("pfjet_muonenergy",pfjet_muonenergy);
   tree_->SetBranchAddress("pfjet_chargedmulti",pfjet_chargedmulti);
   tree_->SetBranchAddress("pfjet_neutralmulti",pfjet_neutralmulti);
   tree_->SetBranchAddress("pfjet_chargedhadronmulti",pfjet_chargedhadronmulti);

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Event Loop
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     if (!isData) weight *= genweight;

     histWeightsH->Fill(1.0,weight);

     float topPt = -1;
     float antitopPt = -1;
     TLorentzVector genBosonLV;genBosonLV.SetXYZT(0,0,0,0);
     TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);

     if (!isData) {
       for (unsigned int igen=0; igen<genparticles_count; ++igen) {
	 TLorentzVector genLV; genLV.SetXYZT(genparticles_px[igen],
					     genparticles_py[igen],
					     genparticles_pz[igen],
					     genparticles_e[igen]);
	 
	 if (genparticles_pdgid[igen]==6){
	   topPt = TMath::Sqrt(genparticles_px[igen]*genparticles_px[igen]+
			       genparticles_py[igen]*genparticles_py[igen]);
	 }
	 if (genparticles_pdgid[igen]==-6){
	   antitopPt = TMath::Sqrt(genparticles_px[igen]*genparticles_px[igen]+
				   genparticles_py[igen]*genparticles_py[igen]);
	 }
	 bool fromHardProcessFinalState = genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1;
	 bool isMuon = false;
	 bool isElectron = false;
	 bool isNeutrino = false;
	 bool isDirectHardProcessTauDecayProduct = genparticles_isDirectHardProcessTauDecayProduct[igen];

	 if (abs(genparticles_pdgid[igen])==11) isElectron = true;
	 if (abs(genparticles_pdgid[igen])==13) isMuon = true;
	 if (abs(genparticles_pdgid[igen])==12||
	     abs(genparticles_pdgid[igen])==14||
	     abs(genparticles_pdgid[igen])==16) isNeutrino = true;

	 bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
	 bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);

	 if (isBoson)        genBosonLV += genLV;
	 if (isVisibleBoson) genVisBosonLV += genLV;
       }
     }
     
     
     if (isDY) { // applying Z pt mass weights
       float zptmassweight = 1;
       float bosonMass = genBosonLV.M();
       float bosonPt = genBosonLV.Pt();
       if (bosonMass>50.0) {
	 float bosonMassX = bosonMass;
	 float bosonPtX = bosonPt;
	 if (bosonMassX>1000.) bosonMassX = 1000.;
	 if (bosonPtX<1.)      bosonPtX = 1.;
	 if (bosonPtX>1000.)   bosonPtX = 1000.;
	 //  zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
	 //  histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
	 correctionWS->var("z_gen_pt")->setVal(bosonPtX);
	 correctionWS->var("z_gen_mass")->setVal(bosonMassX);
	 zptmassweight = correctionWS->function("zpt_weight_nom")->getVal();
       }
       weight *= zptmassweight;
     }

     if (isTOP) { // applying top pT weights
       float topptweight = 1;
       if (topPt>0&&antitopPt>0) topptweight = topPtWeight(topPt,antitopPt,true);
       weight *= topptweight;
     }

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
     }
     puWeightH->Fill(puweight,1.0);
     weight *= puweight;

     // finding HLT filters in the HLT Filter library
     unsigned int singleMuonFilter1   = 0;
     bool isSingleMuonFilter1 = false;
     unsigned int singleMuonFilter2   = 0;
     bool isSingleMuonFilter2 = false;

     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       if (HLTFilter==SingleMuonFilterName1) {
	 singleMuonFilter1 = i;
	 isSingleMuonFilter1 = true;
       }
       if (HLTFilter==SingleMuonFilterName2) {
	 singleMuonFilter2 = i;
	 isSingleMuonFilter2 = true;
       }
     }
     
     if (!isSingleMuonFilter1) {
       cout << "Filter " << SingleMuonFilterName1 << " not found " << endl;
       exit(-1);
     }
     if (!isSingleMuonFilter2) {
       cout << "Filter " << SingleMuonFilterName2 << " not found " << endl;
       exit(-1);
     }

     // Check if it is ztt or zll (for dy sample only)
     if(isDY){
       bool isDYtoMM = false;
       int nPromptElectrons = 0;
       int nPromptMuons     = 0;
       for (unsigned int igen=0; igen<genparticles_count; ++igen) {
	 if(abs(genparticles_pdgid[igen])==11&&genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1){
	   nPromptElectrons += 1;
	 }
	 
	 if(abs(genparticles_pdgid[igen])==13&&genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1){
	   nPromptMuons += 1;
	 }
       }

       if(nPromptMuons == 2) isDYtoMM = true;
       if(isZMM && !isDYtoMM) continue;
       if(!isZMM && isDYtoMM) continue;
     }


     // ********************
     // selecting 2 good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     vector<unsigned int> looseMuons; looseMuons.clear();
     int index1Mu = -1 ;
     int index2Mu = -1 ;
     double leadingMuonPt = 0;
     double relIsoMu1 = 0;
     double relIsoMu2 = 0;
     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       if (isData) {
	 if (event_run >= 278820 && muon_isMedium[i]) muonID = true; // Run2016G-H
	 if (event_run < 278820  && muon_isICHEP[i]) muonID = true; // Run2016B-F
       }
       if(!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       if(muon_pt[i]<ptMuonLooseCut) continue;
       if(fabs(muon_eta[i])>etaMuonCut) continue;
       float neutralHadIsoMu = muon_neutralHadIso[i];
       float photonIsoMu = muon_photonIso[i];
       float chargedHadIsoMu = muon_chargedHadIso[i];
       float puIsoMu = muon_puIso[i];
       float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
       neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
       float absIsoMu = chargedHadIsoMu + neutralIsoMu;
       float relIsoMu = absIsoMu/muon_pt[i];
       if (relIsoMu>isoMuonCut)	continue;
       index2Mu = i;
       relIsoMu2 = relIsoMu;
       if (muon_pt[i]>ptMuonCut && muon_pt[i]>leadingMuonPt){
	 index2Mu = index1Mu;
	 index1Mu = i;
	 leadingMuonPt = muon_pt[i];
	 relIsoMu2 = relIsoMu1;
	 relIsoMu1 = relIsoMu;
       }
     }

     if(index1Mu < 0 || index2Mu < 0) continue;

     TLorentzVector muon1LV; muon1LV.SetXYZM(muon_px[index1Mu],
					    muon_py[index1Mu],
					    muon_pz[index1Mu],
					    MuMass);

     TLorentzVector muon2LV; muon2LV.SetXYZM(muon_px[index2Mu],
					    muon_py[index2Mu],
					    muon_pz[index2Mu],
					    MuMass);


     bool muMatchedTrigger = false;
     for (unsigned int iT=0; iT<trigobject_count; ++iT) {
       float dRtrig = deltaR(muon_eta[index1Mu],muon_phi[index1Mu],trigobject_eta[iT],trigobject_phi[iT]);
       if (dRtrig>DRTrigMatch) continue;
       if (trigobject_filters[iT][singleMuonFilter1] || trigobject_filters[iT][singleMuonFilter2]) muMatchedTrigger = true;
     }

     if (applyTriggerMatch && !muMatchedTrigger) continue;

     double invMassMu1Mu2 = (muon1LV + muon2LV).M();

     // jets 
     int njets = 0;
     double leadingJetPt = 0;
     TLorentzVector jetLV; 
     for (unsigned int ijet=0; ijet<pfjet_count; ++ijet) {
       float absJetEta = fabs(pfjet_eta[ijet]);
       float jetEta = pfjet_eta[ijet];
       float jetPt = pfjet_pt[ijet];
       if (absJetEta>jetEtaCut) continue;
       if (jetPt<jetPtCut) continue;
       bool jetId = looseId(pfjet_e[ijet],
			    pfjet_eta[ijet],
			    pfjet_energycorr[ijet],
			    pfjet_chargedhadronicenergy[ijet],
			    pfjet_neutralhadronicenergy[ijet],
			    pfjet_neutralemenergy[ijet],
			    pfjet_chargedemenergy[ijet],
			    pfjet_muonenergy[ijet],
			    pfjet_chargedmulti[ijet],
			    pfjet_neutralmulti[ijet]);
       if (!jetId) continue;

       float dRMuJet = deltaR(pfjet_eta[ijet],pfjet_phi[ijet],muon1LV.Eta(),muon1LV.Phi());
       if (dRMuJet<deltaRJetLeptonCut) continue;
       dRMuJet = deltaR(pfjet_eta[ijet],pfjet_phi[ijet],muon2LV.Eta(),muon2LV.Phi());
       if (dRMuJet<deltaRJetLeptonCut) continue;
       
       if(pfjet_pt[ijet] > leadingJetPt){
	 jetLV.SetPxPyPzE(pfjet_px[ijet],pfjet_py[ijet],pfjet_pz[ijet],pfjet_e[ijet]);
	 leadingJetPt = pfjet_pt[ijet];
       }
       

       njets++;
     }
     
     // recoil corrections
     int njets_forrecoil = njets;
     if (isW) njets_forrecoil = njets + 1;

     if (isDY||isW) { // applying recoil corrections
       float metcorr_x = metx;
       float metcorr_y = mety;
       float bosonPx = genBosonLV.Px();
       float bosonPy = genBosonLV.Py();
       float lepPx = genVisBosonLV.Px();
       float lepPy = genVisBosonLV.Py();
       recoilMetCorrector.CorrectByMeanResolution(metx,mety,bosonPx,bosonPy,lepPx,lepPy,njets_forrecoil,metcorr_x,metcorr_y);
       metx = metcorr_x;
       mety = metcorr_y;
       met = TMath::Sqrt(metx*metx+mety*mety);
       metphi = TMath::ATan2(mety,metx);
     }

     TLorentzVector Met4; Met4.SetXYZM(metx,mety,0,0);

     float dPhiMuJet = dPhiFrom2P(muon1LV.Px(),muon1LV.Py(),jetLV.Px(),jetLV.Py());
     float dPhiMuMet = dPhiFrom2P(muon1LV.Px(),muon1LV.Py(),metx,mety);
     float mT = TMath::Sqrt(2*met*muon1LV.Pt()*(1-TMath::Cos(dPhiMuMet)));

     //if(mT > 30) continue;
     
     // muon id/iso/trigger weights here
     double muTrigWeight = 1;
     double muIdWeight1  = 1;
     double muIdWeight2  = 1;
     if (!isData) { 

       // from valeria
       //muIdWeight1 = SF_muonId->get_ScaleFactor(muon1LV.Pt(),muon1LV.Eta()); 
       muIdWeight2 = SF_muonId->get_ScaleFactor(muon2LV.Pt(),muon2LV.Eta()); 

       // from KIT
       w_kitLeptonSF -> var("m_eta") -> setVal(muon1LV.Eta());
       w_kitLeptonSF -> var("m_pt")  -> setVal(muon1LV.Pt());
       w_kitLeptonSF -> var("m_iso") -> setVal(relIsoMu1);
       muIdWeight1 = (w_kitLeptonSF->function("m_id_ratio")->getVal())*(w_kitLeptonSF->function("m_iso_binned_ratio")->getVal());
       w_kitLeptonSF -> var("m_eta") -> setVal(muon2LV.Eta());
       w_kitLeptonSF -> var("m_pt")  -> setVal(muon2LV.Pt());
       w_kitLeptonSF -> var("m_iso") -> setVal(relIsoMu2);
       muIdWeight2 = (w_kitLeptonSF->function("m_id_ratio")->getVal())*(w_kitLeptonSF->function("m_iso_binned_ratio")->getVal());
       
       double trigWeightData = SF_muonTrigger->get_EfficiencyData(muon1LV.Pt(),muon1LV.Eta());
       double trigWeightMC   = SF_muonTrigger->get_EfficiencyMC(muon1LV.Pt(),muon1LV.Eta());

       if (applyTriggerMatch) {
	 if (trigWeightMC>0){
	   //muTrigWeight = trigWeightData/trigWeightMC;
	   // From KIT
	   w_kitLeptonSF -> var("m_eta") -> setVal(muon1LV.Eta());
	   w_kitLeptonSF -> var("m_pt")  -> setVal(muon1LV.Pt());
	   w_kitLeptonSF -> var("m_iso") -> setVal(relIsoMu1);
	   muTrigWeight = w_kitLeptonSF->function("m_trgOR4_binned_ratio")->getVal();
	 }
       }
       else
	 muTrigWeight = trigWeightData;
     }

     muTrigWeightH->Fill(muTrigWeight,1.0);
     weight *= muTrigWeight;
     muIdWeightH->Fill(muIdWeight1*muIdWeight2,1.0);
     weight *= muIdWeight1*muIdWeight2;
     ptMuonH->Fill(muon1LV.Pt(),weight);
     etaMuonH->Fill(muon1LV.Eta(),weight);
     mTH->Fill(mT,weight);
     invMassMu1Mu2H->Fill(invMassMu1Mu2,weight);
     metH->Fill(met,weight);
     dPhiMuonMetH->Fill(dPhiMuMet,weight);
     dPhiMuonJetH->Fill(dPhiMuJet,weight);
     nJetH->Fill(njets,weight);
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

 

