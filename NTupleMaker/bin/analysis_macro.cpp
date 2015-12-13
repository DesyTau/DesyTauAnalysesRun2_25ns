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
#include "H2to2H1to4Taus/HTauTau2MuonAnalysis/interface/Config.h"
#include "H2to2H1to4Taus/HTauTau2MuonAnalysis/src/Config.cc"
//#include "H2_2H1_4Taus/HTauTau2MuonAnalysis/interface/Config.h"

#include "TRandom.h"
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

  // trigger
  const string dimuonTriggerName = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName = cfg.get<string>("MuonLowPtFilterName");
  const string dimuonDzFilterName = cfg.get<string>("DimuonDzFilterName");
  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 

  TString DiMuonTriggerName(dimuonTriggerName);
  TString MuonHighPtFilterName(muonHighPtFilterName);
  TString MuonLowPtFilterName(muonLowPtFilterName);
  TString DiMuonDzFilterName(dimuonDzFilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  std::ifstream fileList(argv[2]);

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
 
  
   // Trigger
  unsigned int trigobject_count;
  float trigobject_px[1000];
  float trigobject_py[1000];
  float trigobject_pz[1000];
  float trigobject_pt[1000];
  float  trigobject_eta[1000];
  float trigobject_phi[1000];
  bool trigobject_filters[1000][50];

  //unsigned int iLeadingPosTrig = 0;
  //vector<bool> trigobject_filter; trigobject_filter.clear();

  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();

  std::string rootFileName(argv[2]);
  
  std::string chainName("NTupleProducer/AC1B");
  //std::string chainNameGen("NTupleProducer/AC1B");
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
  TH1F * muonCountH = new TH1F("muonCountH","",11,-0.5,10.5);
  TH1F * nGoodMuonsH = new TH1F("nGoodMuonsH","",11,-0.5,10.5);

  // histograms after dimuon selection
  TH1F * ptLeadingMuH = new TH1F("ptLeadingMuH","",50,0,100);
  TH1F * ptTrailingMuH = new TH1F("ptTrailingMuH","",50,0,100);
  TH1F * etaLeadingMuH = new TH1F("etaLeadingMuH","",50,-2.5,2.5);
  TH1F * etaTrailingMuH = new TH1F("etaTrailingMuH","",50,-2.5,2.5);
  TH1F * dimuonMassH = new TH1F("dimuonMassH","",500,0,500);
  TH1F * nTracksLeadingMuH = new TH1F("nTracksLeadingMuH","",21,-0.5,20.5);
  TH1F * nTracksTrailingMuH = new TH1F("nTracksTrailingMuH","",21,-0.5,20.5);
  TH1F * nSigTracksLeadingMuH = new TH1F("nSigTracksLeadingMuH","",21,-0.5,20.5);
  TH1F * nSigTracksTrailingMuH = new TH1F("nSigTracksTrailingMuH","",21,-0.5,20.5);

  TH1F * ptTrackLeadingMuH = new TH1F("ptTrackLeadingMuH","",100,0,100);
  TH1F * etaTrackLeadingMuH = new TH1F("etaTrackLeadingMuH","",50,-2.5,2.5);
  TH1F * dxyTrackLeadingMuH = new TH1F("dxyTrackLeadingMuH","",200,-0.5,0.5);
  TH1F * dzTrackLeadingMuH = new TH1F("dzTrackLeadingMuH","",200,-1,1);

  TH1F * ptTrackTrailingMuH = new TH1F("ptTrackTrailingMuH","",100,0,100);
  TH1F * etaTrackTrailingMuH = new TH1F("etaTrackTrailingMuH","",50,-2.5,2.5);
  TH1F * dxyTrackTrailingMuH = new TH1F("dxyTrackTrailingMuH","",200,-0.5,0.5);
  TH1F * dzTrackTrailingMuH = new TH1F("dzTrackTrailingMuH","",200,-1,1);

  TH1F * ptLeadingMuSelH = new TH1F("ptLeadingMuSelH","",50,0,100);
  TH1F * ptTrailingMuSelH = new TH1F("ptTrailingMuSelH","",50,0,100);
  TH1F * etaLeadingMuSelH = new TH1F("etaLeadingMuSelH","",50,-2.5,2.5);
  TH1F * etaTrailingMuSelH = new TH1F("etaTrailingMuSelH","",50,-2.5,2.5);
  TH1F * dimuonMassSelH = new TH1F("dimuonMassSelH","",500,0,500);
  TH1F * invMass2Mu2TrkSelH = new TH1F("invMass2Mu2TrkSelH","",500,0,500);

  TH1F * massLeadingMuTrkH = new TH1F("massLeadingMuTrkH","",500,0,50);
  TH1F * massTrailingMuTrkH = new TH1F("massTrailingMuTrkH","",500,0,50);
  TH1F * massMuTrkH = new TH1F("massMuTrkH","",500,0,50);
  TH2F * m1m2SelH = new TH2F("m1m2SelH","",20,0,20,20,0,20);

  // Counters
  TH1F * counter_InputEventsH=new TH1F("counter_InputEventsH","",1,0.,2.);
  TH1F * counterMuonIDH=new TH1F("counterMuonIDH","",1,0.,2.);
  TH1F * counter_MuonIPH=new TH1F("counter_MuonIPH","",1,0.,2.);
  TH1F * counter_MuonEtaH=new TH1F("counter_MuonEtaH","",1,0.,2.);
  TH1F * counter_MuonSizeGTE2H=new TH1F("counter_MuonSizeGTE2H","",1,0.,2.);
  TH1F * counter_MuonSameSignH=new TH1F("counter_MuonSameSignH","",1,0.,2.);
  TH1F * counter_dRMuonH=new TH1F("counter_dRMuonH","",1,0.,2.);
  TH1F * counter_Mu1_Mu2TrigH=new TH1F("counter_Mu1_Mu2TrigH","",1,0.,2.);
  TH1F * counter_MaxMuon_ptSumH=new TH1F("counter_Muon_ptSumH","",1,0.,2.);
  TH1F * counter_MuonKinematicsH=new TH1F("counter_MuonKinematicsH","",1,0.,2.);         
  TH1F * counter_nTracksH=new TH1F("counter_nTracksH","",1,0.,2.);
  TH1F * counter_OneTrackOnlyH=new TH1F("counter_OneTrackOnlyH","",1,0.,2.);                  
  TH1F * counter_ChargeReqLeadMuonTrkH=new TH1F("counter_ChargeReqLeadMuonTrkH","",1,0.,2.);         
  TH1F * counter_ChargeReqTrailMuonTrkH=new TH1F("counter_ChargeReqTrailMuonTrkH","",1,0.,2.);         
  TH1F * counter_TightTrkIPReqH=new TH1F("counter_TightTrkIPReqH","",1,0.,2.);         
  TH1F * counter_TightTrkptReqH=new TH1F("counter_TightTrkptReqH","",1,0.,2.);         
  TH1F * counter_FinalEventsH=new TH1F("counter_FinalEventsH","",1,0.,2.);         

   // Background studies
  TH1F * InvMassN23leadingH = new TH1F("InvMassN23leadingH","",10,0.,10.);
  TH1F * InvMassN23trailingH = new TH1F("InvMassN23trailingH","",10,0.,10.);
  TH1F * InvMassN23H = new TH1F("InvMassN23H","",10,0.,10.);

  TH1F * InvMassHardestNtrk23leadingH = new TH1F("InvMassHardestNtrk23leadingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk23trailingH = new TH1F("InvMassHardestNtrk23trailingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk23H = new TH1F("InvMassHardestNtrk23H","",10,0.,10.);

  TH1F * InvMassSoftestNtrk23leadingH = new TH1F("InvMassSoftestNtrk23leadingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk23trailingH = new TH1F("InvMassSoftestNtrk23trailingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk23H = new TH1F("InvMassSoftestNtrk23H","",10,0.,10.);

  TH1F * InvMassHardestNtrk1leadingH = new TH1F("InvMassHardestNtrk1leadingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk1trailingH = new TH1F("InvMassHardestNtrk1trailingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk1H = new TH1F("InvMassHardestNtrk1H","",10,0.,10.);

  TH1F * InvMassSoftestNtrk1leadingH = new TH1F("InvMassSoftestNtrk1leadingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk1trailingH = new TH1F("InvMassSoftestNtrk1trailingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk1H = new TH1F("InvMassSoftestNtrk1H","",10,0.,10.);

  TH1F * InvMassLeadingH = new TH1F("InvMassLeadingH","",10,0.,10.);
  TH1F * InvMassTrailingH = new TH1F("InvMassTrailingH","",10,0.,10.);
  TH1F * InvMassH = new TH1F("InvMassH","",10,0.,10.);
  TH2F * InvMass2DH = new TH2F("InvMass2DH","",10,0.,10.,10,0.,10.);

  // Correlation Plots
  TH1F * InvMassTrackPlusMuon1D_ControlH = new TH1F("InvMassTrackPlusMuon1D_ControlH","",100,0.,10.); 
  TH2F * InvMassTrackPlusMuon2D_ControlLeadH = new TH2F("InvMassTrackPlusMuon2D_ControlLeadH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_ControlTrailH = new TH2F("InvMassTrackPlusMuon2D_ControlTrailH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_ControlBothH = new TH2F("InvMassTrackPlusMuon2D_ControlBothH","",100,0.,10.,100,0.,10.);
   

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
   tree_->SetBranchAddress("muon_chargedHadIso", &muon_chargedHadIso);
   tree_->SetBranchAddress("muon_neutralHadIso", &muon_neutralHadIso);
   tree_->SetBranchAddress("muon_photonIso", &muon_photonIso);
   tree_->SetBranchAddress("muon_puIso", &muon_puIso);
   tree_->SetBranchAddress("muon_isPF", muon_isPF);
   tree_->SetBranchAddress("muon_isGlobal", muon_isGlobal);
   tree_->SetBranchAddress("muon_isTracker", muon_isTracker);
   tree_->SetBranchAddress("muon_isTight", muon_isTight);
   tree_->SetBranchAddress("muon_isLoose", muon_isLoose);

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
   //   tree_->SetBranchAddress("run_hltnames", &run_hltnames);
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   //   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);
   //     tree_->SetBranchAddress("hltriggerresultsV", &hltriggerresultsV_);
   
   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%1000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     counter_InputEventsH->Fill(1.);
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
     bool isMu8Leg = false;
     bool isMu17Leg = false;
     bool isDZFilter = false;
     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       if (HLTFilter==MuonHighPtFilterName) {
   nMu17Leg = i;
   isMu17Leg = true;
  }
       if (HLTFilter==MuonLowPtFilterName) {
   nMu8Leg = i;
   isMu8Leg = true;
  }
       if (HLTFilter==DiMuonDzFilterName) {
   nDZFilter = i;
   isDZFilter = true;
  }
      }

      if (!isMu17Leg) {
  cout << "Filter " << MuonHighPtFilterName << " not found " << endl;
  exit(-1);
      }
      if (!isMu8Leg) {
  cout << "Filter " << MuonLowPtFilterName << " not found " << endl;
  exit(-1);
      }
      if (!isDZFilter) {
  cout << "Filter " << DiMuonDzFilterName << " not found " << endl;
  exit(-1);
      }
      //      std::cout << std::endl;

      muonCountH->Fill(float(muon_count),weight);

      // selecting good muons
      vector<unsigned int> muons; muons.clear();
      for(UInt_t i=0;i<muon_count;i++){
  if(!muon_isLoose[i]) continue;
  if(!muon_isPF[i]) continue;
           counterMuonIDH->Fill(1.);
  if(fabs(muon_dxy[i])>dxyMuonCut) continue;
  if(fabs(muon_dz[i])>dzMuonCut) continue;
           counter_MuonIPH->Fill(1.);
  if(muon_pt[i]<ptMuonLowCut) continue;
  if(fabs(muon_eta[i])>etaMuonLowCut) continue;
           counter_MuonEtaH->Fill(1.);
  //      cout << "muon pt = " << muon_pt[i] << endl;
  muons.push_back(i);
      }

      nGoodMuonsH->Fill(float(muons.size()),weight);
      
      if (muons.size()<2) continue; // quit event if number of good muons < 2
           counter_MuonSizeGTE2H->Fill(1.);
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
            counter_MuonSameSignH->Fill(1.);
    float dRmuons = deltaR(muon_eta[index1],muon_phi[index1],
         muon_eta[index2],muon_phi[index2]);

    if (dRmuons<dRMuonsCut) continue;
            counter_dRMuonH->Fill(1.);    
    bool mu1MatchMu17 = false;
    bool mu1MatchMu8  = false;
    bool mu1MatchDz   = false;
    for (unsigned int iT=0; iT<trigobject_count; ++iT) {
      float dRtrig = deltaR(muon_eta[index1],muon_phi[index1],
          trigobject_eta[iT],trigobject_phi[iT]);
      if (dRtrig>DRTrigMatch) continue;
      if (trigobject_filters[iT][nMu17Leg] && muon_pt[index1]>ptMuonHighCut && fabs(muon_eta[index1])<etaMuonHighCut)
        mu1MatchMu17 = true;
      if (trigobject_filters[iT][nMu8Leg] && muon_pt[index1]>ptMuonLowCut && fabs(muon_eta[index1])<etaMuonLowCut)
        mu1MatchMu8 = true;
      if (trigobject_filters[iT][nDZFilter])
        mu1MatchDz = true;
    }
    bool mu2MatchMu17 = false;
    bool mu2MatchMu8  = false;
    bool mu2MatchDz   = false;
    for (unsigned int iT=0; iT<trigobject_count; ++iT) {
      float dRtrig = deltaR(muon_eta[index2],muon_phi[index2],
          trigobject_eta[iT],trigobject_phi[iT]);
      if (dRtrig>DRTrigMatch) continue;
      if (trigobject_filters[iT][nMu17Leg] && muon_pt[index2]>ptMuonHighCut && fabs(muon_eta[index2])<etaMuonHighCut)
        mu2MatchMu17 = true;
      if (trigobject_filters[iT][nMu8Leg] && muon_pt[index2]>ptMuonLowCut && fabs(muon_eta[index2])<etaMuonLowCut)
        mu2MatchMu8 = true;
      if (trigobject_filters[iT][nDZFilter])
        mu2MatchDz = true;
    }
    
    // trigger condition
    bool isTriggerMatched = ((mu1MatchMu17&&mu2MatchMu8)||(mu1MatchMu8&&mu2MatchMu17)) && mu1MatchDz && mu2MatchDz;
    if (!isTriggerMatched) continue;
              counter_Mu1_Mu2TrigH->Fill(1.);                 
    if (ptSum>maxPtSum) {
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
          counter_MaxMuon_ptSumH->Fill(1.);                 

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

      // filling histograms (muon kinematics)
      dimuonMassH->Fill(dimuonMass,weight);
      ptLeadingMuH->Fill(muon_pt[iLeading],weight);
      ptTrailingMuH->Fill(muon_pt[iTrailing],weight);
      etaLeadingMuH->Fill(muon_eta[iLeading],weight);
      etaTrailingMuH->Fill(muon_eta[iTrailing],weight);
      counter_MuonKinematicsH->Fill(1.);                  

      // counting tracks around each muon
      std::vector<unsigned int> trkLeadingMu; trkLeadingMu.clear(); // all tracks
      std::vector<unsigned int> trkTrailingMu; trkTrailingMu.clear(); // all tracks
      std::vector<unsigned int> trkSigLeadingMu; trkSigLeadingMu.clear(); // signal tracks
      std::vector<unsigned int> trkSigTrailingMu; trkSigTrailingMu.clear(); // signal tracks
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

      // defining sidebands and signal region

      // sideband N23
      bool isN23leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
      bool isN23trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==2||trkLeadingMu.size()==3); 

      // sidebands Ntrk23
      bool isNtrk23leading  = trkSigLeadingMu.size()>0 && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
      bool isNtrk23trailing = trkSigTrailingMu.size()>0 && (trkLeadingMu.size()==2||trkLeadingMu.size()==3);

      // sidebands Ntrk1
      bool isNtrk1leading  = trkSigLeadingMu.size()>0 && trkTrailingMu.size()==1;
      bool isNtrk1trailing = trkSigTrailingMu.size()>0 && trkLeadingMu.size()==1;

      // signal region
      bool signalRegion = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkSigTrailingMu.size()==1&&trkTrailingMu.size()==1);
      
      counter_nTracksH->Fill(1.);                 

      // sidebands N23
      if (isN23leading) {
	int iTrk = trkSigLeadingMu[0];
	TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					      track_py[iTrk],
					      track_pz[iTrk],
					      track_mass[iTrk]);
	TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
	float mass = TrackPlusMuon4.M();
	InvMassN23leadingH->Fill(mass,weight);
	InvMassN23H->Fill(mass,weight);
      }      
      if (isN23trailing) {
	int iTrk = trkSigTrailingMu[0];
	TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					      track_py[iTrk],
					      track_pz[iTrk],
					      track_mass[iTrk]);
	TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
	float mass = TrackPlusMuon4.M();
	InvMassN23trailingH->Fill(mass,weight);
	InvMassN23H->Fill(mass,weight);
      }      

      // sidebands Ntrk23
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
      
      // sidebands Ntrk1
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
      
      // signal region
      if (signalRegion) {
	counter_FinalEventsH->Fill(1.,weight);
	int iTrkLeading = trkSigTrailingMu[0];
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
	
	InvMassLeadingH->Fill(massTrkMuLeading,weight);
	InvMassTrailingH->Fill(massTrkMuTrailing,weight);
	InvMassH->Fill(massTrkMuLeading,weight);
	InvMassH->Fill(massTrkMuTrailing,weight);
	InvMass2DH->Fill(massTrkMuLeading,massTrkMuTrailing,weight);
      }
      
      // Correlation distributions
      // ****
      // *** Amithabh, you should add conditions that there are
      // *** no additional tracks apart from signal track and soft tracks
      // *** I modified code below
      // ****
      // Old code
      //      bool LeadMuonControl_2Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1) && trkSigTrailingMu.size()==1;
      //      bool TrailMuonControl_2Tracks = (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1) && trkSigLeadingMu.size()==1;
      //      bool LeadMuonControl_3Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2) && trkSigTrailingMu.size()==1;
      //      bool TrailMuonControl_3Tracks = (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2) && trkSigLeadingMu.size()==1;
      //      bool BothMuonControl_2Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1) && (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1) ;
      //      bool BothMuonControl_3Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2) && (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2) ;
      //      bool ControlAll = LeadMuonControl_2Tracks || TrailMuonControl_2Tracks || LeadMuonControl_3Tracks || TrailMuonControl_3Tracks || BothMuonControl_2Tracks || BothMuonControl_3Tracks;
      // New code
      bool signalLeadingMu = trkSigLeadingMu.size()==1 && trkLeadingMu.size()==1;

      bool signalTrailingMu = trkSigTrailingMu.size()==1 && trkTrailingMu.size()==1;

      bool bkgdLeadingMu = 
	(trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
	(trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3);

      bool bkgdTrailingMu = 
	(trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
	(trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3);

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

 

