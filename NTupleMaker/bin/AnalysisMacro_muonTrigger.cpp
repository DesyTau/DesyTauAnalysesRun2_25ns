#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  std::cout << "What the fuck" << std::endl;

  // Data
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection"); 

  // kinematic cuts on muons
  const float ptMuonCut      = cfg.get<float>("ptMuonCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float ptMuonTagCut   = cfg.get<float>("ptMuonTagCut");

  const float etaMuonCut = cfg.get<float>("etaMuonCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonTagCut  = cfg.get<float>("etaMuonTagCut");

  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZleptonsCut   = cfg.get<float>("dZleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");
  int chargeTagMuon          = cfg.get<int>("chargeTagMuon"); // -1/0/+1 (neg,all,pos)

  const string l1seed        = cfg.get<string>("L1Seed");
  const float ptL1TauCut     = cfg.get<float>("ptL1TauCut");
  const float etaL1TauCut    = cfg.get<float>("etaL1TauCut");
  const float detaL1TauCut   = cfg.get<float>("detaL1TauCut");
  const bool matchL1Tau      = cfg.get<bool>("matchL1Tau");


  // tau
  const float ptTauCut  = cfg.get<float>("ptTauCut");
  const float etaTauCut = cfg.get<float>("etaTauCut");
  const float deltaRMuTauCut = cfg.get<float>("deltaRMuTauCut");

  // triggers
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");

  const string hltIsoSingleMu   = cfg.get<string>("HLTIsoSingleMu");  
  const string hltIsoMuProbe       = cfg.get<string>("HLTIsoMuProbe");
  const string hltMu50          = cfg.get<string>("HLTMu50");  

  // HLT filters
  const string hltIsoSingleMuFilter = cfg.get<string>("HLTIsoSingleMuFilter");
  const string hltIsoMuProbeFilter     = cfg.get<string>("HLTIsoMuProbeFilter");
  const string hltIsoMuProbeFilter2    = cfg.get<string>("HLTIsoMuProbeFilter2");
  const string hltMu50Filter  = cfg.get<string>("HLTMu50Filter");

  const unsigned int runRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int runRangeMax = cfg.get<unsigned int>("RunRangeMax");

  const string jsonFile = cfg.get<string>("jsonFile");
  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);

  TString L1Seed(l1seed);
  TString HLTIsoSingleMu(hltIsoSingleMu);
  TString HLTIsoMuProbe(hltIsoMuProbe);
  TString HLTMu50(hltMu50);

  TString HLTIsoSingleMuFilter(hltIsoSingleMuFilter);
  TString HLTIsoMuProbeFilter(hltIsoMuProbeFilter);
  TString HLTIsoMuProbeFilter2(hltIsoMuProbeFilter2);
  TString HLTMu50Filter(hltMu50Filter);
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // **** end of configuration

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

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  bool isHLTIsoSingleMu = 0;

  unsigned int eventNumber = 0;
  unsigned int runNumber = 0;
  unsigned int lumiBlock = 0;

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");
  /*
  TTree * treePrescales = new TTree("Prescales","Prescales");
  treePrescales->Branch("Run",&runNumber,"Run/i");
  treePrescales->Branch("Event",&eventNumber,"Event/i");
  treePrescales->Branch("Lumi",&lumiBlock,"Run/i");

  treePrescales->Branch("HLTDoubleMu",&isHLTDoubleMu,"HLTDoubleMu/O");
  treePrescales->Branch("HLTDoubleMuDZ",&isHLTDoubleMuDZ,"HLTDoubleMuDZ/O");
  treePrescales->Branch("HLTDoubleMuSameSignDZ",&isHLTDoubleMuSameSignDZ,"HLTDoubleMuSameSignDZ/O");
  treePrescales->Branch("HLTIsoSingleMu",&isHLTIsoSingleMu,"HLTIsoSingleMu/O");
  treePrescales->Branch("HLTIsoMu18",&isHLTIsoMu18,"HLTIsoMu18/O");
  treePrescales->Branch("HLTIsoMu20",&isHLTIsoMu20,"HLTIsoMu20/O");
  treePrescales->Branch("HLTIsoMu22",&isHLTIsoMu22,"HLTIsoMu22/O");
  treePrescales->Branch("HLTMu27",&isHLTMu27,"HLTMu27/O");
  treePrescales->Branch("HLTMu50",&isHLTMu50,"HLTMu50/O");

  treePrescales->Branch("PrescaleHLTDoubleMu",&presHLTDoubleMu,"PrescaleHLTDoubleMu/I");
  treePrescales->Branch("PrescaleHLTDoubleMuDZ",&presHLTDoubleMuDZ,"PrescaleHLTDoubleMuDZ/I");
  treePrescales->Branch("PrescaleHLTDoubleMuSameSignDZ",&presHLTDoubleMuSameSignDZ,"PrescaleHLTDoubleMuSameSignDZ/I");
  treePrescales->Branch("PrescaleHLTIsoSingleMu",&presHLTIsoSingleMu,"PrescaleHLTIsoSingleMu/I");
  treePrescales->Branch("PrescaleHLTIsoMu18",&presHLTIsoMu18,"PrescaleHLTIsoMu18/I");
  treePrescales->Branch("PrescaleHLTIsoMu20",&presHLTIsoMu20,"PrescaleHLTIsoMu20/I");
  treePrescales->Branch("PrescaleHLTIsoMu22",&presHLTIsoMu22,"PrescaleHLTIsoMu22/I");
  treePrescales->Branch("PrescaleHLTMu27",&presHLTMu27,"PrescaleHLTMu27/I");
  treePrescales->Branch("PrescaleHLTMu50",&presHLTMu50,"PrescaleHLTMu50/I");
  */

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);

  int nEtaBins = 4;
  float etaBins[5] = {0, 0.9, 1.2, 2.1, 2.4};

  int nEtaFineBins = 14;
  float etaFineBins[15] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};

  int nPtBins = 18;
  float ptBins[19] = {10, 12, 14, 16, 18, 
		      20, 22, 24, 26, 28,
		      30, 32, 34, 36, 38,
		      40, 45, 50, 100}; 

  int nPtBins50 = 9;
  float ptBins50[10] = {20,25,30,35,40,45,50,55,60,1000};

  int nDzBins = 6;
  float dzBins[7] = {-0.001,0.05,0.1,0.15,0.2,0.3,0.5};

  int nDPhiBins = 20;
  float dPhiBins[21];
  for (int ip=0; ip<=20; ++ip)
    dPhiBins[ip] = 0.05*TMath::Pi()*float(ip);

  int nDRBins = 10;
  float dRBins[11] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

  TH1F * ZMassIsoMuLegPassH = new TH1F("ZMassIsoMuLegPassH","",60,60,120);
  TH1F * ZMassIsoMuLegFailH = new TH1F("ZMassIsoMuLegFailH","",60,60,120);

  TH1F * ZMassMu50LegPassH = new TH1F("ZMassMu50LegPassH","",60,60,120);
  TH1F * ZMassMu50LegFailH = new TH1F("ZMassMu50FailH","",60,60,120);


  TString EtaBins[4] = {"EtaLt0p9","Eta0p9to1p2","Eta1p2to2p1","EtaGt2p1"};
  TString EtaFineBins[14];
  for (int iB=0; iB<nEtaFineBins; ++iB) {
    char BinChar[4];  
    if (iB<10)
      sprintf(BinChar,"%1i",iB);
    else 
      sprintf(BinChar,"%2i",iB);
    EtaFineBins[iB] = "EtaBin" + TString(BinChar);
  }

  //  TString PtBins[18] = {"Pt5to7","Pt7to9","Pt9to11",
  //			"Pt11to13","Pt13to15","Pt15to17","Pt17to19",
  //			"Pt19to21","Pt21to23","Pt23to25","Pt25to27",
  //			"Pt27to30","Pt30to40","Pt40to50","Pt50to60",
  //			"Pt60to80","Pt80to100","PtGt100"};

  TString PtBins[18] = {"Pt10to12","Pt12to14","Pt14to16","Pt16to18","Pt18to20",
			"Pt20to22","Pt22to24","Pt24to26","Pt26to28","Pt28to30",
			"Pt30to32","Pt32to34","Pt34to36","Pt36to38","Pt38to40",
			"Pt40to45","Pt45to50","PtGt50"};

  TString DzBins[6] = {"Dz0to0p5","Dz0p5to1p0","Dz1p0to1p5",
		       "Dz1p5to2p0","Dz2p0to3p0","DzGt3p0"};
  TString PtBins50[9] = {"Pt20to25","Pt25to30","Pt30to35",
			 "Pt35to40","Pt40to45","Pt45to50",
			 "Pt50to55","Pt55to60","Pt60toInf"};

  TString DRBins[10] = {"dRLt0p5","dR0p5to1p0","dR1p0to1p5","dR1p5to2p0","dR2p0to2p5",
			"dR2p5to3p0","dR3p0to3p5","dR3p5to4p0","dR4p0to4p5","dRGt4p5"};
  TString DPhiBins[20] = {"dPhi1","dPhi2","dPhi3","dPhi4","dPhi5",
			  "dPhi6","dPhi7","dPhi8","dPhi9","dPhi10",
			  "dPhi11","dPhi12","dPhi13","dPhi14","dPhi15",
			  "dPhi16","dPhi17","dPhi18","dPhi19","dPhi20"};

  TH1F * EtaBinsH  = new TH1F("EtaBinsH","",nEtaBins,etaBins);
  TH1F * EtaFineBinsH  = new TH1F("EtaFineBinsH","",nEtaFineBins,etaFineBins);
  TH1F * PtBinsH   = new TH1F("PtBinsH","",nPtBins,ptBins);
  TH1F * DzBinsH   = new TH1F("DzBinsH","",nDzBins,dzBins); 
  TH1F * PtBins50H = new TH1F("PtBins50H","",nPtBins50,ptBins50);


  TH1F * DPhiBinsH = new TH1F("DPhiBinsH","",nDPhiBins,dPhiBins); 
  TH1F * DRBinsH   = new TH1F("DRBinsH","",nDRBins,dRBins);

  for (int iB=0; iB<nEtaFineBins; ++iB)
    EtaFineBinsH->GetXaxis()->SetBinLabel(iB+1,EtaFineBins[iB]);

  for (int iB=0; iB<nEtaBins; ++iB)
    EtaBinsH->GetXaxis()->SetBinLabel(iB+1,EtaBins[iB]);

  for (int iB=0; iB<nDPhiBins; ++iB)
    DPhiBinsH->GetXaxis()->SetBinLabel(iB+1,DPhiBins[iB]);

  for (int iB=0; iB<nPtBins; ++iB)
    PtBinsH->GetXaxis()->SetBinLabel(iB+1,PtBins[iB]);

  for (int iB=0; iB<nPtBins50; ++iB)
    PtBins50H->GetXaxis()->SetBinLabel(iB+1,PtBins50[iB]);

  for (int iB=0; iB<nDzBins; ++iB)
    DzBinsH->GetXaxis()->SetBinLabel(iB+1,DzBins[iB]);

  // (Pt,Eta)

  TH1F * ZMassIsoMuLegPtEtaPassH[4][18];
  TH1F * ZMassIsoMuLegPtEtaFailH[4][18];

  TH1F * ZMassLowPtLegPtEtaPassH[4][18];
  TH1F * ZMassLowPtLegPtEtaFailH[4][18];

  TH1F * ZMassIsoMuLegEtaPassH[4];
  TH1F * ZMassIsoMuLegEtaFailH[4];

  TH1F * ZMassIsoMuLegEtaFinePassH[14];
  TH1F * ZMassIsoMuLegEtaFineFailH[14];

  TH1F * ZMassIsoMuDPhiMuMuPassH[10];
  TH1F * ZMassIsoMuDPhiMuMuFailH[10];

  // IsoMu standard
  TH1F * ZMassIsoMuLegPtPassH[18];
  TH1F * ZMassIsoMuLegPtFailH[18];

  // IsoMu relative to reco tau
  TH1F * ZMassIsoMuToTauLegPtPassH[18];
  TH1F * ZMassIsoMuToTauLegPtFailH[18];

  // IsoMu relative to L1Seed
  TH1F * ZMassIsoMuToL1SeedLegPtPassH[18];
  TH1F * ZMassIsoMuToL1SeedLegPtFailH[18];

  // L1Seed standard
  TH1F * ZMassL1SeedLegPtPassH[18];
  TH1F * ZMassL1SeedLegPtFailH[18];

  // L1Seed to reco tau
  TH1F * ZMassL1SeedToTauLegPtPassH[18];
  TH1F * ZMassL1SeedToTauLegPtFailH[18];

  // L1object 
  TH1F * ZMassL1ObjectPtPassH[18];
  TH1F * ZMassL1ObjectPtFailH[18];

  // L1Seed to l1 object
  TH1F * ZMassL1SeedToL1ObjectLegPtPassH[18];
  TH1F * ZMassL1SeedToL1ObjectLegPtFailH[18];

  TH1F * ZMassL1SeedToL1ObjectNotTagFiredLegPtPassH[18];
  TH1F * ZMassL1SeedToL1ObjectNotTagFiredLegPtFailH[18];

  TH1F * ZMassL1SeedToL1ObjectTagFiredLegPtPassH[18];
  TH1F * ZMassL1SeedToL1ObjectTagFiredLegPtFailH[18];

  TH1F * ZMassL1SeedToL1ObjectOrTagFiredLegPtPassH[18];
  TH1F * ZMassL1SeedToL1ObjectOrTagFiredLegPtFailH[18];

  TH1F * ZMassTagL1SeedToL1ObjectLegPtPassH[18];
  TH1F * ZMassTagL1SeedToL1ObjectLegPtFailH[18];

  TH1F * ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtPassH[18];
  TH1F * ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtFailH[18];

  TH1F * ZMassTagL1SeedToL1ObjectProbeFiredLegPtPassH[18];
  TH1F * ZMassTagL1SeedToL1ObjectProbeFiredLegPtFailH[18];

  TH1F * ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtPassH[18];
  TH1F * ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtFailH[18];

  // dPhi dependece
  TH1F * ZMassL1SeedToL1ObjectLegDPhiPassH[20];
  TH1F * ZMassL1SeedToL1ObjectLegDPhiFailH[20];
  

  // Mu50 efficiency
  TH1F * ZMassMu50LegPtPassH[9];
  TH1F * ZMassMu50LegPtFailH[9];

  TH1F * ZMassMu50LegEtaPassH[4];
  TH1F * ZMassMu50LegEtaFailH[4];

  TH1F * ZMassMu50LegPtEtaPassH[4][9];
  TH1F * ZMassMu50LegPtEtaFailH[4][9];

  // dPhi
  for (int iDPhi=0; iDPhi<nDPhiBins; ++iDPhi) {
    ZMassL1SeedToL1ObjectLegDPhiPassH[iDPhi] = new TH1F("ZMassL1SeedToL1ObjectLeg_"+DPhiBins[iDPhi]+"_PassH","",60,60,120);
    ZMassL1SeedToL1ObjectLegDPhiFailH[iDPhi] = new TH1F("ZMassL1SeedToL1ObjectLeg_"+DPhiBins[iDPhi]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaFineBins; ++iEta) {
    ZMassIsoMuLegEtaFinePassH[iEta] = new TH1F("ZMassIsoMuLeg_"+EtaFineBins[iEta]+"_PassH","",60,60,120);
    ZMassIsoMuLegEtaFineFailH[iEta] = new TH1F("ZMassIsoMuLeg_"+EtaFineBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    ZMassIsoMuLegEtaPassH[iEta] = new TH1F("ZMassIsoMuLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassIsoMuLegEtaFailH[iEta] = new TH1F("ZMassIsoMuLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
    ZMassMu50LegEtaPassH[iEta]   = new TH1F("ZMassMu50Leg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassMu50LegEtaFailH[iEta]   = new TH1F("ZMassMu50Leg_"+EtaBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iPt=0; iPt<nPtBins; ++iPt) {

      ZMassIsoMuLegPtPassH[iPt] = new TH1F("ZMassIsoMuLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassIsoMuLegPtFailH[iPt] = new TH1F("ZMassIsoMuLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassIsoMuToTauLegPtPassH[iPt] = new TH1F("ZMassIsoMuToTauLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassIsoMuToTauLegPtFailH[iPt] = new TH1F("ZMassIsoMuToTauLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassIsoMuToL1SeedLegPtPassH[iPt] = new TH1F("ZMassIsoMuToL1SeedLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassIsoMuToL1SeedLegPtFailH[iPt] = new TH1F("ZMassIsoMuToL1SeedLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassL1SeedLegPtPassH[iPt] = new TH1F("ZMassL1SeedLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedLegPtFailH[iPt] = new TH1F("ZMassL1SeedLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassL1SeedToTauLegPtPassH[iPt] = new TH1F("ZMassL1SeedToTauLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedToTauLegPtFailH[iPt] = new TH1F("ZMassL1SeedToTauLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      //

      ZMassL1ObjectPtPassH[iPt] = new TH1F("ZMassL1Object_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1ObjectPtFailH[iPt] = new TH1F("ZMassL1Object_"+PtBins[iPt]+"_FailH","",60,60,120);

      // Probe muon

      ZMassL1SeedToL1ObjectLegPtPassH[iPt] = new TH1F("ZMassL1SeedToL1ObjectLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedToL1ObjectLegPtFailH[iPt] = new TH1F("ZMassL1SeedToL1ObjectLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassL1SeedToL1ObjectNotTagFiredLegPtPassH[iPt] = new TH1F("ZMassL1SeedToL1ObjectNotTagFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedToL1ObjectNotTagFiredLegPtFailH[iPt] = new TH1F("ZMassL1SeedToL1ObjectNotTagFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassL1SeedToL1ObjectTagFiredLegPtPassH[iPt] = new TH1F("ZMassL1SeedToL1ObjectTagFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedToL1ObjectTagFiredLegPtFailH[iPt] = new TH1F("ZMassL1SeedToL1ObjectTagFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassL1SeedToL1ObjectOrTagFiredLegPtPassH[iPt] = new TH1F("ZMassL1SeedToL1ObjectOrTagFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassL1SeedToL1ObjectOrTagFiredLegPtFailH[iPt] = new TH1F("ZMassL1SeedToL1ObjectOrTagFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      // Tag muon

      ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtPassH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectNotProbeFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtFailH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectNotProbeFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassTagL1SeedToL1ObjectProbeFiredLegPtPassH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectProbeFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassTagL1SeedToL1ObjectProbeFiredLegPtFailH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectProbeFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtPassH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectOrProbeFiredLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtFailH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectOrProbeFiredLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

      ZMassTagL1SeedToL1ObjectLegPtPassH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassTagL1SeedToL1ObjectLegPtFailH[iPt] = new TH1F("ZMassTagL1SeedToL1ObjectLeg_"+PtBins[iPt]+"_FailH","",60,60,120);

  }

  for (int iPt=0; iPt<nPtBins50; ++iPt) {
      ZMassMu50LegPtPassH[iPt] = new TH1F("ZMassMu50Leg_"+PtBins50[iPt]+"_PassH","",60,60,120);
      ZMassMu50LegPtFailH[iPt] = new TH1F("ZMassMu50Leg_"+PtBins50[iPt]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {

    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassIsoMuLegPtEtaPassH[iEta][iPt] = new TH1F("ZMassIsoMuLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassIsoMuLegPtEtaFailH[iEta][iPt] = new TH1F("ZMassIsoMuLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
    }

    for (int iPt=0; iPt<nPtBins50; ++iPt) {
      ZMassMu50LegPtEtaPassH[iEta][iPt] = new TH1F("ZMassMu50Leg_"+EtaBins[iEta]+"_"+PtBins50[iPt]+"_PassH","",60,60,120);
      ZMassMu50LegPtEtaFailH[iEta][iPt] = new TH1F("ZMassMu50Leg_"+EtaBins[iEta]+"_"+PtBins50[iPt]+"_FailH","",60,60,120);
    }

  }

  unsigned int iRun;
  unsigned int iEvent;
  //  TTree * eventTree = new TTree("eventTree","eventTree");
  //  eventTree->Branch("Run",&iRun,"Run/i");
  //  eventTree->Branch("Event",&iEvent,"Event/i");

  int nFiles = 0;
  int nEvents = 0;

  int selEvents = 0;
  int selEventsHLTDoubleMu = 0;
  int selEventsHLTDoubleMuDZ = 0;
  int selEventsHLTDoubleMuSameSignDZ = 0;

  int selPairs = 0;
  int selPairsHLTDoubleMu = 0;
  int selPairsHLTDoubleMuDZ = 0;
  int selPairsHLTDoubleMuSameSignDZ = 0;

  unsigned int minRun = 99999999;
  unsigned int maxRun = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
  std::vector<unsigned int> allGoodRuns; allGoodRuns.clear();


  // PileUp
  PileUp * PUofficial = new PileUp();
  TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUDataFile,"read");
  TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUMCFile, "read");
  TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
  TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get("pileup");
  //  PUofficial->set_h_data(PUOfficial_data);
  //  PUofficial->set_h_MC(PUOfficial_mc);

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("initroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      if (analysisTree.event_run>maxRun)
	maxRun = analysisTree.event_run;

      if (analysisTree.event_run<minRun)
	minRun = analysisTree.event_run;


      bool isNewRun = true;
      if (allRuns.size()>0) {
	for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	  if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }

      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);

      weightsH->Fill(0.0,weight);

      if (isData) {
	if (applyGoodRunSelection) {
	  bool lumi = false;
	  int n=analysisTree.event_run;
	  int lum = analysisTree.event_luminosityblock;
	  
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

      if (analysisTree.event_run<runRangeMin) continue;
      if (analysisTree.event_run>runRangeMax) continue;

      isNewRun = true;
      if (allGoodRuns.size()>0) {
	for (unsigned int iR=0; iR<allGoodRuns.size(); ++iR) {
	  if (analysisTree.event_run==allGoodRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }

      if (isNewRun) 
	allGoodRuns.push_back(analysisTree.event_run);


      if (!isData) {
	//	float puWeight =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	float puWeight =  float(PUOfficial_data->GetBinContent(PUOfficial_data->GetXaxis()->FindBin(analysisTree.numtruepileupinteractions)));
	weight *= puWeight;
      }

      //      std::cout << "Triggers" << std::endl;

      // triggers
      isHLTIsoSingleMu = false;
      
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTIsoSingleMu)) {
	  if (it->second==1)
            isHLTIsoSingleMu = true;
	}
      }
      
      //      runNumber = analysisTree.event_run;
      //      eventNumber = analysisTree.event_nr;
      //      lumiBlock = analysisTree.event_luminosityblock;
      //      treePrescales->Fill();

      if (!applyTrigger) {
	isHLTIsoSingleMu = true;
      }

      unsigned int nIsoSingleMuFilter = 0;
      bool isIsoSingleMuFilter = false;
      unsigned int nIsoMuProbeFilter = 0;
      bool isIsoMuProbeFilter = false;
      unsigned int nIsoMuProbeFilter2 = 0;
      bool isIsoMuProbeFilter2 = false;
      unsigned int nMu50Filter = 0;
      bool isMu50Filter = false;

      unsigned int nL1Seed = 0;
      bool isL1Seed = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==HLTIsoSingleMuFilter) {
	  nIsoSingleMuFilter = i;
	  isIsoSingleMuFilter = true;
	}
	if (HLTFilter==HLTIsoMuProbeFilter) {
	  nIsoMuProbeFilter = i;
	  isIsoMuProbeFilter = true;
	}
	if (HLTFilter==HLTIsoMuProbeFilter2) {
	  nIsoMuProbeFilter2  = i;
	  isIsoMuProbeFilter2 = true;
	}
	if (HLTFilter==HLTMu50Filter) {
	  nMu50Filter = i;
	  isMu50Filter = true;
	}
	if (HLTFilter==L1Seed) {
	  nL1Seed = i;
	  isL1Seed = true;
	}
      }

      if (!isIsoSingleMuFilter)
	std::cout << "Filter " << HLTIsoSingleMuFilter << "  does not exist" << std::endl;
      if (!isIsoMuProbeFilter)
	std::cout << "Filter " << HLTIsoMuProbeFilter << "  does not exist" << std::endl;
      if (!isIsoMuProbeFilter2)
	std::cout << "Filter " << HLTIsoMuProbeFilter2 << "  does not exist" << std::endl;
      if (!isMu50Filter)
	std::cout << "Filter " << HLTMu50Filter << "  does not exist" << std::endl;
      //      if (!isL1Seed)
      //	std::cout << "L1Seed " << L1Seed << "  does not exist" << std::endl;

      /*
      int numL1muons = 0;
      for (unsigned int imu=0; imu<analysisTree.l1muon_count; ++imu) {
        TLorentzVector l1muonLV; l1muonLV.SetXYZM(analysisTree.l1muon_px[imu],
                                                  analysisTree.l1muon_py[imu],
                                                  analysisTree.l1muon_pz[imu],
                                                  muonMass);
        if (l1muonLV.Pt()>20&&fabs(l1muonLV.Eta())<2.1) {
	  numL1muons++;
	}
      }

      if (numL1muons>1) {

      }
      */

      //


      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isTight[im]) continue;
	muons.push_back(im);
      }

      if (muons.size()<2) continue;

      bool isPairSelected = false;
      bool isPairSelectedHLTDoubleMu = false;
      bool isPairSelectedHLTDoubleMuDZ = false;
      bool isPairSelectedHLTDoubleMuSameSignDZ = false;

      // selecting muon pair
      for (unsigned int im1=0; im1<muons.size(); ++im1) {
	//	  std::cout << "Muon " << im << std::endl;
	int  mu1Index = muons[im1];
	bool mu1MatchIsoSingleMu = false;
	bool mu1MatchL1Seed = false;
	bool mu1MatchL1object = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nIsoSingleMuFilter]) // HLT_IsoMu filter
	    mu1MatchIsoSingleMu = true;
	  if (analysisTree.trigobject_filters[iT][nL1Seed]) // L1Seed Mu+Tau
	    mu1MatchL1Seed = true;
	}

	bool mu1IsoSingleMu = 
	  mu1MatchIsoSingleMu && 
	  analysisTree.muon_pt[mu1Index]>ptMuonTagCut && 
	  fabs(analysisTree.muon_eta[mu1Index])<etaMuonTagCut;

	float q1 = analysisTree.muon_charge[mu1Index];
	
	for (unsigned int im2=0; im2<muons.size(); ++im2) {

	  if (im1==im2) continue;

	  int  mu2Index = muons[im2];

	  float q2 = analysisTree.muon_charge[mu2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index]);

	  bool mu2MatchMu50   = false;
	  bool mu2MatchIsoMuProbe1 = false;
	  bool mu2MatchIsoMuProbe2 = false;
	  bool mu2MatchL1object = false;
	  bool mu2MatchL1Seed   = false; 
	  bool tauIsFound       = false;

	  for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) {
	    if (analysisTree.tau_pt[itau]<ptTauCut) continue;
	    if (fabs(analysisTree.tau_eta[itau])>etaTauCut) continue;
	    if (analysisTree.tau_decayModeFinding[itau]<0.5) continue;
	    if (analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[itau]<0.5) continue;
	    if (analysisTree.tau_againstElectronLooseMVA6[itau]<0.5) continue;
	    if (analysisTree.tau_againstMuonTight3[itau]<0.5) continue;
	    float deltaRmu1tau = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
					analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
	    if (deltaRmu1tau<deltaRMuTauCut) continue;
	    float deltaRmu2tau = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
					analysisTree.tau_eta[itau],analysisTree.tau_phi[itau]);
	    if (deltaRmu2tau<deltaRMuTauCut) continue;
	    tauIsFound = true;

	  }


	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMu50Filter]) // HLT_Mu50 filter
	      mu2MatchMu50 = true;
	    if (isIsoMuProbeFilter) {
	      if (analysisTree.trigobject_filters[iT][nIsoMuProbeFilter]) // HLT_IsoMu filter
		mu2MatchIsoMuProbe1 = true;
	    }
	    if (isIsoMuProbeFilter2) {
	      if (analysisTree.trigobject_filters[iT][nIsoMuProbeFilter2]) // HLT_IsoMu filter
		mu2MatchIsoMuProbe2 = true;
	    }
	    if (analysisTree.trigobject_filters[iT][nL1Seed]) // L1 Seed
	      mu2MatchL1Seed = true;
	  }
	  bool mu2MatchIsoMuProbe = mu2MatchIsoMuProbe1 || mu2MatchIsoMuProbe2;
	  bool l1muonFound = false;
	  for (unsigned int imu=0; imu<analysisTree.l1muon_count; ++imu) {                                                                                                                  
	    TLorentzVector l1muonLV; l1muonLV.SetXYZM(analysisTree.l1muon_px[imu],                                                                                                          
						      analysisTree.l1muon_py[imu],                                                                                                          
						      analysisTree.l1muon_pz[imu],                                                                                                          
						      muonMass);                                                                                                                            
	    if (analysisTree.l1muon_iso[imu]>0 && l1muonLV.Pt()>20 && fabs(l1muonLV.Eta())<2.1) {
	      float deltaRl1mu = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
					l1muonLV.Eta(),l1muonLV.Phi());
	      if (deltaRl1mu<DRTrigMatch) mu2MatchL1object = true;
	      deltaRl1mu = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				  l1muonLV.Eta(),l1muonLV.Phi());
	      if (deltaRl1mu<DRTrigMatch) mu1MatchL1object = true;

	    }                                                                                                                                                                               
	  }

	  if (fabs(analysisTree.muon_eta[mu2Index])<ptMuonHighCut)


	  bool mu2Mu50 = mu2MatchMu50;
	  bool mu2IsoMuProbe = mu2MatchIsoMuProbe;

	  float dZ = fabs(analysisTree.muon_dz[mu1Index]-analysisTree.muon_dz[mu2Index]);

	  TLorentzVector mu1lv; mu1lv.SetXYZM(analysisTree.muon_px[mu1Index],
					      analysisTree.muon_py[mu1Index],
					      analysisTree.muon_pz[mu1Index],
					      muonMass);
	  TLorentzVector mu2lv; mu2lv.SetXYZM(analysisTree.muon_px[mu2Index],
                                              analysisTree.muon_py[mu2Index],
                                              analysisTree.muon_pz[mu2Index],
                                              muonMass);

	  float mass = (mu1lv+mu2lv).M();

	  if (mass<50.0) continue;

	  float absIso1 = analysisTree.muon_chargedHadIso[mu1Index];
	  float neutralIso1 = 
	    analysisTree.muon_neutralHadIso[mu1Index] +
            analysisTree.muon_photonIso[mu1Index] -
            0.5*analysisTree.muon_puIso[mu1Index];
          neutralIso1 = TMath::Max(float(0),neutralIso1);
          absIso1 += neutralIso1;
	  float relIso1 = absIso1/analysisTree.muon_pt[mu1Index];

	  float absIso2 = analysisTree.muon_chargedHadIso[mu2Index];
	  float neutralIso2 = 
	    analysisTree.muon_neutralHadIso[mu2Index] +
            analysisTree.muon_photonIso[mu2Index] -
            0.5*analysisTree.muon_puIso[mu2Index];
          neutralIso2 = TMath::Max(float(0),neutralIso2);
          absIso2 += neutralIso2;
	  float relIso2 = absIso2/analysisTree.muon_pt[mu2Index];

	  float mu1RelIso = relIso1;
	  float mu2RelIso = relIso2;

	  if (analysisTree.muon_pt[mu2Index]>analysisTree.muon_pt[mu1Index]) {
	    float temp = relIso1;
	    relIso1 = relIso2;
	    relIso2 = temp;
	  }

	  bool dirIso = (relIso2<isoMuonCut) && (relIso1<isoMuonCut);

	  float dPhiMuMu = dPhiFrom2P(mu1lv.Px(),mu1lv.Py(),
				      mu2lv.Px(),mu2lv.Py());
	  
	  float dPhiMuMuX = TMath::Max(float(0.001),TMath::Min(dPhiMuMu,float(TMath::Pi()-0.001)));
	  int dPhiMuMuBin = binNumber(dPhiMuMuX,nDPhiBins,dPhiBins);


	  bool foundL1Seed = false;

          for (unsigned int iT=0; iT<analysisTree.l1tau_count; ++iT) {

	    TLorentzVector l1LV; l1LV.SetXYZM(analysisTree.l1tau_px[iT],
					      analysisTree.l1tau_py[iT],
					      analysisTree.l1tau_pz[iT],
					      pionMass);
	    
	    if (analysisTree.l1tau_bx[iT]!=0) continue;
            if (l1LV.Pt()<ptL1TauCut) continue;
	    if (fabs(l1LV.Eta())>etaL1TauCut) continue;
            float dRl1tau = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
                                   l1LV.Eta(),l1LV.Phi());
	    if (dRl1tau<DRTrigMatch) continue;
	    dRl1tau = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
			     l1LV.Eta(),l1LV.Phi());
	    if (dRl1tau<DRTrigMatch) continue;
	    /*  
	    std::cout << "Yep : Muon2 eta = " << analysisTree.muon_eta[mu2Index] 
		      << "   phi = " << analysisTree.muon_phi[mu2Index] 
		      << "   pt = " << analysisTree.muon_pt[mu2Index] << std::endl;
	    std::cout << "   l1tau eta = " << analysisTree.l1isotau_eta[iT]
		      << "   phi = " << analysisTree.l1isotau_phi[iT]
		      << "   pt = "  << analysisTree.l1isotau_pt[iT] << std::endl;
	    */
            foundL1Seed = true;
          }
	  bool printEvent = analysisTree.event_run==297487 && analysisTree.event_nr==77934257;
	  printEvent = printEvent || (analysisTree.event_run==297488 && analysisTree.event_nr==270262596);
	  printEvent = printEvent || (analysisTree.event_run==297488 && analysisTree.event_nr==269862269);
	  
	  if (printEvent) {
	    std::cout << "Run = " << analysisTree.event_run 
		      << "     Lumi = " << analysisTree.event_luminosityblock 
		      << "     Event = " << analysisTree.event_nr << std::endl;
	    printf("muon1 : pt = %6.2f   eta = %5.2f   phi = %5.2f   iso = %5.2f\n",
		   mu1lv.Pt(),
		   mu1lv.Eta(),
		   mu1lv.Phi(),
		   relIso1);

	    printf("muon2 : pt = %6.2f   eta = %5.2f   phi = %5.2f   iso = %5.2f   HLTMatch = %1i    L1Match = %1i\n",
		   mu2lv.Pt(),
		   mu2lv.Eta(),
		   mu2lv.Phi(),
		   relIso2,
		   mu2MatchIsoMuProbe,
		   mu2MatchL1Seed);

	    printf("mass = %6.2f\n",mass);

	    // Check of L1 Seeds
	    for (unsigned int imu=0; imu<analysisTree.l1muon_count; ++imu) {
	      TLorentzVector l1muonLV; l1muonLV.SetXYZM(analysisTree.l1muon_px[imu],
							analysisTree.l1muon_py[imu],
							analysisTree.l1muon_pz[imu],
							muonMass);
	      if (l1muonLV.Pt()>20&&fabs(l1muonLV.Eta())<2.1) {
		printf("l1muon : pt = %6.2f   eta = %5.2f   phi = %5.2f   iso = %2i\n",
		       l1muonLV.Pt(),
		       l1muonLV.Eta(),
		       l1muonLV.Phi(),
		       analysisTree.l1muon_iso[imu]);
	      }
	    }
	    
	    for (unsigned int itau=0; itau<analysisTree.l1tau_count; ++itau) {
	      TLorentzVector l1tauLV; l1tauLV.SetXYZM(analysisTree.l1tau_px[itau],
						      analysisTree.l1tau_py[itau],
						      analysisTree.l1tau_pz[itau],
						      pionMass);
	      if (l1tauLV.Pt()>24&&fabs(l1tauLV.Eta())<2.1) {
		printf("l1tau  : pt = %6.2f   eta = %5.2f   phi = %5.2f   iso = %2i\n",
		       l1tauLV.Pt(),
		       l1tauLV.Eta(),
		       l1tauLV.Phi(),
		       analysisTree.l1tau_iso[itau]);
	      }
	    }
	    
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      if (analysisTree.trigobject_filters[iT][nL1Seed]) 
		printf("L1Seed object : pt = %6.2f   eta = %5.2f   phi = %5.2f\n",
		       analysisTree.trigobject_pt[iT],
		       analysisTree.trigobject_eta[iT],
		       analysisTree.trigobject_phi[iT]);
	    }
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      if (analysisTree.trigobject_filters[iT][nIsoMuProbeFilter]) 
		printf("IsoMu20 object : pt = %6.2f   eta = %5.2f   phi = %5.2f\n",
		       analysisTree.trigobject_pt[iT],
		       analysisTree.trigobject_eta[iT],
		       analysisTree.trigobject_phi[iT]);
	    }
	    std::cout << "deltaPhi(mu1,mu2) = " << dPhiMuMu << std::endl;
	    std::cout << "L1 tau found = " << foundL1Seed << std::endl;
	    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	    std::cout << std::endl;
	    std::cout << std::endl;
	  }
	 
	  if (!matchL1Tau) foundL1Seed = true;

	  float mu2AbsEta = fabs(analysisTree.muon_eta[mu2Index]);
	  float mu2Pt = TMath::Max(float(5.01),TMath::Min(float(analysisTree.muon_pt[mu2Index]),float(99.9)));
	  int etaBin  = binNumber(mu2AbsEta,nEtaBins,etaBins);
	  float etaFine = analysisTree.muon_eta[mu2Index];
	  if (etaFine<-2.4) etaFine = -2.39;
	  if (etaFine>2.4) etaFine = 2.39;
	  int etaFineBin = binNumber(etaFine,nEtaFineBins,etaFineBins);
	  int ptBin   = binNumber(mu2Pt,nPtBins,ptBins);
	  int ptBin50 = binNumber(mu2Pt,nPtBins50,ptBins50);

	  float mu1AbsEta = fabs(analysisTree.muon_eta[mu1Index]);
	  float mu1Pt = TMath::Max(float(5.01),TMath::Min(float(analysisTree.muon_pt[mu1Index]),float(99.9)));
	  int ptBinMu1   = binNumber(mu1Pt,nPtBins,ptBins);

	  bool chargeTagPassed = true;
	  if (chargeTagMuon<0 && analysisTree.muon_charge[mu1Index]>0) chargeTagPassed = false;
	  if (chargeTagMuon>0 && analysisTree.muon_charge[mu1Index]<0) chargeTagPassed = false;
	  
	  if (!chargeTagPassed) continue;

	  // L1 object matching
	  if ((mu2AbsEta<etaMuonHighCut) && (mu2RelIso<isoMuonCut)) {
	    if (mu2MatchL1object) ZMassL1ObjectPtPassH[ptBin]->Fill(mass,weight);
	    else ZMassL1ObjectPtFailH[ptBin]->Fill(mass,weight);
	      
	  }

	  // tau presence
	  if (isHLTIsoSingleMu && mu1IsoSingleMu && tauIsFound && dRmumu>dRleptonsCut) {
	    if ( (mu2AbsEta<etaMuonHighCut) && (mu2RelIso<isoMuonCut) ) { // isolated muon path

	      if (mu2MatchL1Seed) ZMassL1SeedToTauLegPtPassH[ptBin]->Fill(mass,weight);
              else ZMassL1SeedToTauLegPtFailH[ptBin]->Fill(mass,weight);

	      if (mu2MatchIsoMuProbe) ZMassIsoMuToTauLegPtPassH[ptBin]->Fill(mass,weight);
	      else ZMassIsoMuToTauLegPtFailH[ptBin]->Fill(mass,weight);

	    }
	  }
	  // l1 seed
	  if (isHLTIsoSingleMu && mu1IsoSingleMu && mu2MatchL1Seed && dRmumu>dRleptonsCut) {
            if ( (mu2AbsEta<etaMuonHighCut) && (mu2RelIso<isoMuonCut) ) { // isolated muon path

              if (mu2MatchIsoMuProbe) ZMassIsoMuToL1SeedLegPtPassH[ptBin]->Fill(mass,weight);
              else ZMassIsoMuToL1SeedLegPtFailH[ptBin]->Fill(mass,weight);

            }
          }
	  // l1 tau presence and l1 muon presence 
          if (isHLTIsoSingleMu && mu1IsoSingleMu && mu2MatchL1object && foundL1Seed && dRmumu>dRleptonsCut) { // Single muon selection
	    if ( (mu2AbsEta<etaMuonHighCut) && (mu2RelIso<isoMuonCut) ) { // isolated muon path

	      if (mu2MatchL1Seed) { 
		ZMassL1SeedToL1ObjectLegPtPassH[ptBin]->Fill(mass,weight);
		ZMassL1SeedToL1ObjectLegDPhiPassH[dPhiMuMuBin]->Fill(mass,weight);
	      }
	      else {
		ZMassL1SeedToL1ObjectLegPtFailH[ptBin]->Fill(mass,weight);
		ZMassL1SeedToL1ObjectLegDPhiFailH[dPhiMuMuBin]->Fill(mass,weight);
	      }
	      if (mu2MatchL1Seed&&!mu1MatchL1Seed) ZMassL1SeedToL1ObjectNotTagFiredLegPtPassH[ptBin]->Fill(mass,weight);
              else ZMassL1SeedToL1ObjectNotTagFiredLegPtFailH[ptBin]->Fill(mass,weight);

	      if (mu2MatchL1Seed&&mu1MatchL1Seed) ZMassL1SeedToL1ObjectTagFiredLegPtPassH[ptBin]->Fill(mass,weight);
              else ZMassL1SeedToL1ObjectTagFiredLegPtFailH[ptBin]->Fill(mass,weight);

	      if (mu2MatchL1Seed||mu1MatchL1Seed) ZMassL1SeedToL1ObjectOrTagFiredLegPtPassH[ptBin]->Fill(mass,weight);
              else ZMassL1SeedToL1ObjectOrTagFiredLegPtFailH[ptBin]->Fill(mass,weight);

	    }
	  }

	  // l1 tau presence and l1 muon presence (tag muon)
	  if (isHLTIsoSingleMu && mu1IsoSingleMu && mu1MatchL1object && foundL1Seed && dRmumu>dRleptonsCut) { // Single muon selection
            if ( (mu1AbsEta<etaMuonHighCut) && (mu1RelIso<isoMuonCut) ) { // isolated muon path

              if (mu1MatchL1Seed) ZMassTagL1SeedToL1ObjectLegPtPassH[ptBinMu1]->Fill(mass,weight);
              else ZMassTagL1SeedToL1ObjectLegPtFailH[ptBinMu1]->Fill(mass,weight);

              if (mu1MatchL1Seed&&!mu2MatchL1Seed) ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtPassH[ptBinMu1]->Fill(mass,weight);
              else ZMassTagL1SeedToL1ObjectNotProbeFiredLegPtFailH[ptBinMu1]->Fill(mass,weight);

              if (mu1MatchL1Seed&&mu2MatchL1Seed) ZMassTagL1SeedToL1ObjectProbeFiredLegPtPassH[ptBinMu1]->Fill(mass,weight);
              else ZMassTagL1SeedToL1ObjectProbeFiredLegPtFailH[ptBinMu1]->Fill(mass,weight);

              if (mu1MatchL1Seed||mu2MatchL1Seed) ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtPassH[ptBinMu1]->Fill(mass,weight);
              else ZMassTagL1SeedToL1ObjectOrProbeFiredLegPtFailH[ptBinMu1]->Fill(mass,weight);


            }
          }

	  // l1 tau presence
	  if (isHLTIsoSingleMu && mu1IsoSingleMu && foundL1Seed && dRmumu>dRleptonsCut) { // Single muon selection

	    isPairSelected = true;
	    selPairs++;

	    if ( (mu2AbsEta<etaMuonHighCut) && (mu2RelIso<isoMuonCut) ) { // isolated muon path

	      if (mu2MatchL1Seed) ZMassL1SeedLegPtPassH[ptBin]->Fill(mass,weight);
	      else ZMassL1SeedLegPtFailH[ptBin]->Fill(mass,weight);

	      if (mu2MatchIsoMuProbe) {
		ZMassIsoMuLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		ZMassIsoMuLegPtPassH[ptBin]->Fill(mass,weight);
	      }
	      else {
		ZMassIsoMuLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);
		ZMassIsoMuLegPtFailH[ptBin]->Fill(mass,weight);
	      }
	      
	      if (mu2Pt>ptMuonHighCut) {
		if (mu2MatchIsoMuProbe) {
		  ZMassIsoMuLegEtaPassH[etaBin]->Fill(mass,weight);
		  ZMassIsoMuLegEtaFinePassH[etaFineBin]->Fill(mass,weight);
		  ZMassIsoMuLegPassH->Fill(mass,weight);
		} 
		else {
		  ZMassIsoMuLegEtaFailH[etaBin]->Fill(mass,weight);
		  ZMassIsoMuLegEtaFineFailH[etaFineBin]->Fill(mass,weight);
		  ZMassIsoMuLegFailH->Fill(mass,weight);
		}
	      }


	    }
	      
	    if(mu2AbsEta<2.4) { // Mu50 path
	      
	      if (mu2MatchMu50) {
		ZMassMu50LegPtEtaPassH[etaBin][ptBin50]->Fill(mass,weight);
		ZMassMu50LegPtPassH[ptBin50]->Fill(mass,weight);
	      }
	      else {
		ZMassMu50LegPtEtaFailH[etaBin][ptBin50]->Fill(mass,weight);
		ZMassMu50LegPtFailH[ptBin50]->Fill(mass,weight);
	      }
	      
	      if (analysisTree.muon_pt[mu2Index]>55) {
		if (mu2MatchIsoMuProbe) {
		  ZMassMu50LegEtaPassH[etaBin]->Fill(mass,weight);
		  ZMassMu50LegPassH->Fill(mass,weight);
		}
                  else {
                    ZMassMu50LegEtaFailH[etaBin]->Fill(mass,weight);
		    ZMassMu50LegFailH->Fill(mass,weight);
		  }
	      }
	    }
	  }
	}
      }

      if (isPairSelected)
	selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << "Total number of selected pairs  = " << selPairs << std::endl;
  std::cout << std::endl;
  std::cout << "Run range " << minRun << ":" << maxRun << std::endl;
  std::cout << std::endl;
  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "runs      : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  std::sort (allGoodRuns.begin(), allGoodRuns.end(), myobject);
  std::cout << "good runs : ";
  for (unsigned int iR=0; iR<allGoodRuns.size(); ++iR)
    std::cout << " " << allGoodRuns.at(iR);
  std::cout << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}




//  LocalWords:  puMCFile
