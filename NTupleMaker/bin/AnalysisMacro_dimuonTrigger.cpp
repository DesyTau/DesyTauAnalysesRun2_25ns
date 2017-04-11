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

  // Data
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection"); 

  // kinematic cuts on muons
  const float ptMuonCut      = cfg.get<float>("ptMuonCut");
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float ptMuonTagCut   = cfg.get<float>("ptMuonTagCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonTagCut  = cfg.get<float>("etaMuonTagCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  const bool applyMuonIso    = cfg.get<bool>("ApplyMuonIso");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZleptonsCut   = cfg.get<float>("dZleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");
  int chargeTagMuon          = cfg.get<int>("chargeTagMuon"); // -1/0/+1 (neg,all,pos)

  // triggers
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");

  const string hltDoubleMu           = cfg.get<string>("HLTDoubleMu");
  const string hltDoubleMuDZ         = cfg.get<string>("HLTDoubleMuDZ");
  const string hltDoubleMuSameSignDZ = cfg.get<string>("HLTDoubleMuSameSignDZ");

  const string hltIsoSingleMu   = cfg.get<string>("HLTIsoSingleMu");  
  const string hltIsoMu18       = cfg.get<string>("HLTIsoMu18");  
  const string hltIsoMu20       = cfg.get<string>("HLTIsoMu20");  
  const string hltIsoMu22       = cfg.get<string>("HLTIsoMu22");  
  const string hltMu27          = cfg.get<string>("HLTMu27");  
  const string hltMu45Eta2p1    = cfg.get<string>("HLTMu45Eta2p1");  


  // HLT filters
  const string hltHighPtLeg   = cfg.get<string>("HLTHighPtLeg");
  const string hltLowPtLeg1    = cfg.get<string>("HLTLowPtLeg1");
  const string hltLowPtLeg2    = cfg.get<string>("HLTLowPtLeg2");
  const string dzFilter       = cfg.get<string>("dzFilter"); 
  const string sameSignFilter = cfg.get<string>("sameSignFilter");

  const string hltIsoSingleMuFilter = cfg.get<string>("HLTIsoSingleMuFilter");
  const string hltIsoMu18Filter     = cfg.get<string>("HLTIsoMu18Filter");
  const string hltIsoMu20Filter     = cfg.get<string>("HLTIsoMu20Filter");
  const string hltIsoMu22Filter     = cfg.get<string>("HLTIsoMu22Filter");
  const string hltMu27Filter        = cfg.get<string>("HLTMu27Filter");
  const string hltMu45Eta2p1Filter  = cfg.get<string>("HLTMu45Eta2p1Filter");

  const unsigned int runRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int runRangeMax = cfg.get<unsigned int>("RunRangeMax");

  const string jsonFile = cfg.get<string>("jsonFile");
  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);

  // convesrsion from string to TString
  TString HLTDoubleMu(hltDoubleMu);
  TString HLTDoubleMuDZ(hltDoubleMuDZ);
  TString HLTDoubleMuSameSignDZ(hltDoubleMuSameSignDZ);

  TString HLTIsoSingleMu(hltIsoSingleMu);
  TString HLTIsoMu18(hltIsoMu18);
  TString HLTIsoMu20(hltIsoMu20);
  TString HLTIsoMu22(hltIsoMu22);
  TString HLTMu27(hltMu27);
  TString HLTMu45Eta2p1(hltMu45Eta2p1);

  TString HLTHighPtLeg(hltHighPtLeg);
  TString HLTLowPtLeg1(hltLowPtLeg1);
  TString HLTLowPtLeg2(hltLowPtLeg2);

  TString DZFilter(dzFilter);
  TString SameSignFilter(sameSignFilter);

  TString HLTIsoSingleMuFilter(hltIsoSingleMuFilter);
  TString HLTIsoMu18Filter(hltIsoMu18Filter);
  TString HLTIsoMu20Filter(hltIsoMu20Filter);
  TString HLTIsoMu22Filter(hltIsoMu22Filter);
  TString HLTMu27Filter(hltMu27Filter);
  TString HLTMu45Eta2p1Filter(hltMu45Eta2p1Filter);
  
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

  int presHLTDoubleMu = 0;
  int presHLTDoubleMuDZ = 0;
  int presHLTDoubleMuSameSignDZ = 0;
  int presHLTIsoSingleMu = 0;
  int presHLTIsoMu18 = 0;
  int presHLTIsoMu20 = 0;
  int presHLTIsoMu22 = 0;
  int presHLTMu27 = 0;
  int presHLTMu45Eta2p1 = 0;

  bool isHLTDoubleMu = 0;
  bool isHLTDoubleMuDZ = 0;
  bool isHLTDoubleMuSameSignDZ = 0;
  bool isHLTIsoSingleMu = 0;
  bool isHLTIsoMu18 = 0;
  bool isHLTIsoMu20 = 0;
  bool isHLTIsoMu22 = 0;
  bool isHLTMu27 = 0;
  bool isHLTMu45Eta2p1 = 0;

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
  treePrescales->Branch("HLTMu45Eta2p1",&isHLTMu45Eta2p1,"HLTMu45Eta2p1/O");

  treePrescales->Branch("PrescaleHLTDoubleMu",&presHLTDoubleMu,"PrescaleHLTDoubleMu/I");
  treePrescales->Branch("PrescaleHLTDoubleMuDZ",&presHLTDoubleMuDZ,"PrescaleHLTDoubleMuDZ/I");
  treePrescales->Branch("PrescaleHLTDoubleMuSameSignDZ",&presHLTDoubleMuSameSignDZ,"PrescaleHLTDoubleMuSameSignDZ/I");
  treePrescales->Branch("PrescaleHLTIsoSingleMu",&presHLTIsoSingleMu,"PrescaleHLTIsoSingleMu/I");
  treePrescales->Branch("PrescaleHLTIsoMu18",&presHLTIsoMu18,"PrescaleHLTIsoMu18/I");
  treePrescales->Branch("PrescaleHLTIsoMu20",&presHLTIsoMu20,"PrescaleHLTIsoMu20/I");
  treePrescales->Branch("PrescaleHLTIsoMu22",&presHLTIsoMu22,"PrescaleHLTIsoMu22/I");
  treePrescales->Branch("PrescaleHLTMu27",&presHLTMu27,"PrescaleHLTMu27/I");
  treePrescales->Branch("PrescaleHLTMu45Eta2p1",&presHLTMu45Eta2p1,"PrescaleHLTMu45Eta2p1/I");
  */

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);

  // prescales
  TH1F * prescaleHLTDoubleMuH = new TH1F("prescaleHLTDoubleMuH","",21,-0.5,20.5);
  TH1F * prescaleHLTDoubleMuDZH = new TH1F("prescaleHLTDoubleMuDZH","",21,-0.5,20.5);
  TH1F * prescaleHLTDoubleMuSameSignDZH = new TH1F("prescaleHLTDoubleMuSameSignDZH","",21,-0.5,20.5);
  TH1F * prescaleHLTIsoSingleMuH = new TH1F("prescaleHLTIsoSingleMuH","",21,-0.5,20.5);
  TH1F * prescaleHLTIsoMu18H = new TH1F("prescaleHLTIsoMu18H","",101,-0.5,100.5);
  TH1F * prescaleHLTIsoMu20H = new TH1F("prescaleHLTIsoMu20H","",101,-0.5,100.5);
  TH1F * prescaleHLTIsoMu22H = new TH1F("prescaleHLTIsoMu22H","",101,-0.5,100.5);
  TH1F * prescaleHLTMu27H = new TH1F("prescaleHLTMu27H","",101,-0.5,100.5);
  TH1F * prescaleHLTMu45Eta2p1H = new TH1F("prescaleHLTMu45Eta2p1H","",101,-0.5,100.5);

  // J/Psi ->
  TH1F * JPsiMassDZFilterPassH =  new TH1F("JPsiMassDZFilterPassH","",200,2,4);
  TH1F * JPsiMassDZFilterFailH =  new TH1F("JPsiMassDZFilterFailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz0to1PassH =  new TH1F("JPsiMassDZFilterDz0to1PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz0to1FailH =  new TH1F("JPsiMassDZFilterDz0to1FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz1to2PassH =  new TH1F("JPsiMassDZFilterDz1to2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz1to2FailH =  new TH1F("JPsiMassDZFilterDz1to2FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDzGt2PassH =  new TH1F("JPsiMassDZFilterDzGt2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDzGt2FailH =  new TH1F("JPsiMassDZFilterDzGt2FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterPassH =  new TH1F("JPsiMassSameSignFilterPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterFailH =  new TH1F("JPsiMassSameSignFilterFailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDR0to0p15PassH =  new TH1F("JPsiMassSameSignFilterDR0to0p15PassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDR0to0p15FailH =  new TH1F("JPsiMassSameSignFilterDR0to0p15FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDRGt0p15PassH =  new TH1F("JPsiMassSameSignFilterDRGt0p15PassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDRGt0p15FailH =  new TH1F("JPsiMassSameSignFilterDRGt0p15FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDirIsoPassH =  new TH1F("JPsiMassSameSignFilterDirIsoPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDirIsoFailH =  new TH1F("JPsiMassSameSignFilterDirIsoFailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterInvIsoPassH =  new TH1F("JPsiMassSameSignFilterInvIsoPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterInvIsoFailH =  new TH1F("JPsiMassSameSignFilterInvIsoFailH","",200,2,4);

  TH1F * muIsoLeadJPsiHLTDoubleMuH   = new TH1F("muIsoLeadJPsiHLTDoubleMuH","",200,0,2);
  TH1F * muIsoLeadJPsiHLTDoubleMuDZH = new TH1F("muIsoLeadJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * muIsoTrailJPsiHLTDoubleMuH   = new TH1F("muIsoTrailJPsiHLTDoubleMuH","",200,0,2);
  TH1F * muIsoTrailJPsiHLTDoubleMuDZH = new TH1F("muIsoTrailJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * dRmumuJPsiHLTDoubleMuH   = new TH1F("dRmumuJPsiHLTDoubleMuH","",200,0,2);
  TH1F * dRmumuJPsiHLTDoubleMuDZH = new TH1F("dRmumuJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * dZmumuJPsiHLTDoubleMuH   = new TH1F("dZmumuJPsiHLTDoubleMuH","",100,0,1);
  TH1F * dZmumuJPsiHLTDoubleMuDZH = new TH1F("dZmumuJPsiHLTDoubleMuDZH","",100,0,1);

  // Z ->
  TH1F * ZMassDZFilterPassH =  new TH1F("ZMassDZFilterPassH","",60,60,120);
  TH1F * ZMassDZFilterFailH =  new TH1F("ZMassDZFilterFailH","",60,60,120);

  TH1F * ZMassDZFilterDz0to1PassH =  new TH1F("ZMassDZFilterDz0to1PassH","",60,60,120);
  TH1F * ZMassDZFilterDz0to1FailH =  new TH1F("ZMassDZFilterDz0to1FailH","",60,60,120);

  TH1F * ZMassDZFilterDz1to2PassH =  new TH1F("ZMassDZFilterDz1to2PassH","",60,60,120);
  TH1F * ZMassDZFilterDz1to2FailH =  new TH1F("ZMassDZFilterDz1to2FailH","",60,60,120);

  TH1F * ZMassDZFilterDzGt2PassH =  new TH1F("ZMassDZFilterDzGt2PassH","",60,60,120);
  TH1F * ZMassDZFilterDzGt2FailH =  new TH1F("ZMassDZFilterDzGt2FailH","",60,60,120);

  TH1F * ZMassSameSignFilterPassH =  new TH1F("ZMassSameSignFilterPassH","",60,60,120);
  TH1F * ZMassSameSignFilterFailH =  new TH1F("ZMassSameSignFilterFailH","",60,60,120);

  int nEtaBins = 4;
  float etaBins[5] = {-0.001, 0.9, 1.2, 2.1, 2.4};
  int nPtBins = 18;
  float ptBins[19] = { 5,  7,  9, 11, 13, 
		      15, 17, 19, 21, 23, 
		      25, 27, 30, 40, 50,
		      60, 80, 100, 200};

  int nPtBins45 = 9;
  float ptBins45[10] = {21,26,31,36,41,46,51,56,70,1000};

  int nDzBins = 6;
  float dzBins[7] = {-0.001,0.05,0.1,0.15,0.2,0.3,0.5};


  TH1F * ZMassHighPtLegPassH = new TH1F("ZMassHighPtLegPassH","",60,60,120);
  TH1F * ZMassHighPtLegFailH = new TH1F("ZMassHighPtLegFailH","",60,60,120);
  TH1F * ZMassLowPtLegPassH = new TH1F("ZMassLowPtLegPassH","",60,60,120);
  TH1F * ZMassLowPtLegFailH = new TH1F("ZMassLowPtLegFailH","",60,60,120);


  TString EtaBins[4] = {"EtaLt0p9","Eta0p9to1p2","Eta1p2to2p1","EtaGt2p1"};
  TString PtBins[18] = {"Pt5to7","Pt7to9","Pt9to11",
			"Pt11to13","Pt13to15","Pt15to17","Pt17to19",
			"Pt19to21","Pt21to23","Pt23to25","Pt25to27",
			"Pt27to30","Pt30to40","Pt40to50","Pt50to60",
			"Pt60to80","Pt80to100","PtGt100"};
  TString DzBins[6] = {"Dz0to0p5","Dz0p5to1p0","Dz1p0to1p5",
		       "Dz1p5to2p0","Dz2p0to3p0","DzGt3p0"};
  TString PtBins45[9] = {"Pt21to26","Pt26to31","Pt31to36",
			 "Pt36to41","Pt41to46","Pt46to51",
			 "Pt51to56","Pt56to70","Pt70toInf"};

  TH1F * EtaBinsH = new TH1F("EtaBinsH","",nEtaBins,etaBins);
  TH1F * PtBinsH  = new TH1F("PtBinsH","",nPtBins,ptBins);
  TH1F * DzBinsH  = new TH1F("DzBinsH","",nDzBins,dzBins); 
  TH1F * PtBins45H = new TH1F("PtBins45H","",nPtBins45,ptBins45);

  for (int iB=0; iB<nEtaBins; ++iB)
    EtaBinsH->GetXaxis()->SetBinLabel(iB+1,EtaBins[iB]);

  for (int iB=0; iB<nPtBins; ++iB)
    PtBinsH->GetXaxis()->SetBinLabel(iB+1,PtBins[iB]);

  for (int iB=0; iB<nPtBins45; ++iB)
    PtBins45H->GetXaxis()->SetBinLabel(iB+1,PtBins45[iB]);

  for (int iB=0; iB<nDzBins; ++iB)
    DzBinsH->GetXaxis()->SetBinLabel(iB+1,DzBins[iB]);

  // (Pt,Eta)

  TH1F * ZMassHighPtLegPtEtaPassH[4][18];
  TH1F * ZMassHighPtLegPtEtaFailH[4][18];

  TH1F * ZMassLowPtLegPtEtaPassH[4][18];
  TH1F * ZMassLowPtLegPtEtaFailH[4][18];

  // Eta dependence

  TH1F * ZMassHighPtLegEtaPassH[4];
  TH1F * ZMassHighPtLegEtaFailH[4];

  TH1F * ZMassLowPtLegEtaPassH[4];
  TH1F * ZMassLowPtLegEtaFailH[4];

  // Pt dependence

  TH1F * ZMassHighPtLegPtPassH[18];
  TH1F * ZMassHighPtLegPtFailH[18];

  TH1F * ZMassLowPtLegPtPassH[18];
  TH1F * ZMassLowPtLegPtFailH[18];
  TH1F * JPsiMassDzPassH[6];
  TH1F * JPsiMassDzFailH[6];

  TH1F * ZMassMu45LegPtPassH[9];
  TH1F * ZMassMu45LegPtFailH[9];

  TH1F * ZMassMu45LegEtaPassH[4];
  TH1F * ZMassMu45LegEtaFailH[4];

  TH1F * ZMassMu45LegPtEtaPassH[4][9];
  TH1F * ZMassMu45LegPtEtaFailH[4][9];

  for (int iDz=0; iDz<nDzBins; ++iDz) {
    JPsiMassDzPassH[iDz] = new TH1F("JPsiMass_"+DzBins[iDz]+"_PassH","",200,2,4);
    JPsiMassDzFailH[iDz] = new TH1F("JPsiMass_"+DzBins[iDz]+"_FailH","",200,2,4);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    ZMassHighPtLegEtaPassH[iEta] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassHighPtLegEtaFailH[iEta] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
    ZMassLowPtLegEtaPassH[iEta]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassLowPtLegEtaFailH[iEta]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
    ZMassMu45LegEtaPassH[iEta]   = new TH1F("ZMassMu45Leg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassMu45LegEtaFailH[iEta]   = new TH1F("ZMassMu45Leg_"+EtaBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassHighPtLegPtPassH[iPt] = new TH1F("ZMassHighPtLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassHighPtLegPtFailH[iPt] = new TH1F("ZMassHighPtLeg_"+PtBins[iPt]+"_FailH","",60,60,120);
      ZMassLowPtLegPtPassH[iPt]  = new TH1F("ZMassLowPtLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassLowPtLegPtFailH[iPt]  = new TH1F("ZMassLowPtLeg_"+PtBins[iPt]+"_FailH","",60,60,120);
  }

  for (int iPt=0; iPt<nPtBins45; ++iPt) {
      ZMassMu45LegPtPassH[iPt] = new TH1F("ZMassMu45Leg_"+PtBins45[iPt]+"_PassH","",60,60,120);
      ZMassMu45LegPtFailH[iPt] = new TH1F("ZMassMu45Leg_"+PtBins45[iPt]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {

    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassHighPtLegPtEtaPassH[iEta][iPt] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassHighPtLegPtEtaFailH[iEta][iPt] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
      ZMassLowPtLegPtEtaPassH[iEta][iPt]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassLowPtLegPtEtaFailH[iEta][iPt]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
    }

    for (int iPt=0; iPt<nPtBins45; ++iPt) {
      ZMassMu45LegPtEtaPassH[iEta][iPt] = new TH1F("ZMassMu45Leg_"+EtaBins[iEta]+"_"+PtBins45[iPt]+"_PassH","",60,60,120);
      ZMassMu45LegPtEtaFailH[iEta][iPt] = new TH1F("ZMassMu45Leg_"+EtaBins[iEta]+"_"+PtBins45[iPt]+"_FailH","",60,60,120);
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
  PUofficial->set_h_data(PUOfficial_data);
  PUofficial->set_h_MC(PUOfficial_mc);

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
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
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
	if (analysisTree.event_run<runRangeMin) continue;
	if (analysisTree.event_run>runRangeMax) continue;
      }

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

      weightsH->Fill(0.0,weight);

      if (!isData) {
	float puWeight =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	weight *= puWeight;
      }
 

      // triggers
      isHLTDoubleMu = false;
      isHLTDoubleMuDZ = false;
      isHLTDoubleMuSameSignDZ = false;
      isHLTIsoSingleMu = false;
      isHLTIsoMu18 = false;
      isHLTIsoMu20 = false;
      isHLTIsoMu22 = false;
      isHLTMu27 = false;
      isHLTMu45Eta2p1 = false;
      
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTDoubleMu)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
	    isHLTDoubleMu = true;
	}
	if (trigName.Contains(HLTDoubleMuDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTDoubleMuDZ = true;
	}
	if (trigName.Contains(HLTDoubleMuSameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTDoubleMuSameSignDZ = true;
	}
	if (trigName.Contains(HLTIsoSingleMu)) {
	  if (it->second==1)
            isHLTIsoSingleMu = true;
	}
	if (trigName.Contains(HLTIsoMu18)) {
	  if (it->second==1)
            isHLTIsoMu18 = true;
	}
	if (trigName.Contains(HLTIsoMu20)) {
	  if (it->second==1)
            isHLTIsoMu20 = true;
	}
	if (trigName.Contains(HLTIsoMu22)) {
	  if (it->second==1)
            isHLTIsoMu22 = true;
	}
	if (trigName.Contains(HLTMu27)) {
	  if (it->second==1)
            isHLTMu27 = true;
	}
	if (trigName.Contains(HLTMu45Eta2p1)) {
	  if (it->second==1)
	    isHLTMu45Eta2p1 = true;
	}

      }
      //      std::cout << "Filters found " << std::endl;

      presHLTDoubleMu = 0;
      presHLTDoubleMuDZ = 0;
      presHLTDoubleMuSameSignDZ = 0;
      presHLTIsoSingleMu = 0;
      presHLTIsoMu18 = 0;
      presHLTIsoMu20 = 0;
      presHLTIsoMu22 = 0;
      presHLTMu27 = 0;
      presHLTMu45Eta2p1 = 0;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerprescales->begin(); it!=analysisTree.hltriggerprescales->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTDoubleMu)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMu = it->second;
	}
	if (trigName.Contains(HLTDoubleMuDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMuDZ = it->second;
	}
	if (trigName.Contains(HLTDoubleMuSameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMuSameSignDZ = it->second;
	}
	if (trigName.Contains(HLTIsoSingleMu)) {
	  presHLTIsoSingleMu = it->second;
	}
	if (trigName.Contains(HLTIsoMu18)) {
	  presHLTIsoMu18 = it->second;
	  //	  if (presHLTIsoMu18>100)
	  //	    cout << "HLTIsoMu18 prescale : " << presHLTIsoMu18 << std::endl;
	}
	if (trigName.Contains(HLTIsoMu20)) {
	  presHLTIsoMu20 = it->second;
	  // if (presHLTIsoMu20>100)
	  //   cout << "HLTIsoMu20 prescale : " << presHLTIsoMu20 << std::endl;
	}
	if (trigName.Contains(HLTIsoMu22)) {
	  presHLTIsoMu22 = it->second;
	  // if (presHLTIsoMu22>100)
	  //   cout << "HLTIsoMu22 prescale : " << presHLTIsoMu22 << std::endl;
	}
	if (trigName.Contains(HLTMu27)) {
	  presHLTMu27 = it->second;
	  // if (presHLTMu27>100)
	  //   cout << "HLTMu27 prescale : " << presHLTMu27 << std::endl;
	}
	if (trigName.Contains(HLTMu45Eta2p1)) {
	  presHLTMu45Eta2p1 = it->second;
	}
      }
      
      runNumber = analysisTree.event_run;
      eventNumber = analysisTree.event_nr;
      lumiBlock = analysisTree.event_luminosityblock;
      //      treePrescales->Fill();

      if (!applyTrigger) {
	isHLTDoubleMu = true;
	isHLTDoubleMuDZ = true;
	isHLTDoubleMuSameSignDZ = true;
	isHLTIsoSingleMu = true;
      }

      unsigned int nHighPtLeg = 0;
      bool isHighPtLeg = false;
      unsigned int nLowPtLeg = 0;
      bool isLowPtLeg = false;
      unsigned int nDZFilter = 0;
      bool isDZFilter = false;
      unsigned int nSameSignFilter = 0;
      bool isSameSignFilter = false;

      unsigned int nIsoSingleMuFilter = 0;
      bool isIsoSingleMuFilter;
      unsigned int nIsoMu20Filter = 0;
      bool isIsoMu20Filter;
      unsigned int nIsoMu18Filter = 0;
      bool isIsoMu18Filter;
      unsigned int nIsoMu22Filter = 0;
      bool isIsoMu22Filter;
      unsigned int nMu27Filter = 0;
      bool isMu27Filter;
      unsigned int nMu45Eta2p1Filter = 0;
      bool isMu45Eta2p1Filter;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==HLTHighPtLeg) {
	  nHighPtLeg = i;
	  isHighPtLeg = true;
	}
	if (HLTFilter==HLTLowPtLeg1||HLTFilter==HLTLowPtLeg2) {
	  nLowPtLeg = i;
	  isLowPtLeg = true;
	}
	if (HLTFilter==DZFilter) {
	  nDZFilter = i;
	  isDZFilter = true;
	}
	if (HLTFilter==SameSignFilter) {
	  nSameSignFilter = i;
	  isSameSignFilter = true;
	}
	if (HLTFilter==HLTIsoSingleMuFilter) {
	  nIsoSingleMuFilter = i;
	  isIsoSingleMuFilter = true;
	}
	if (HLTFilter==HLTIsoMu20Filter) {
	  nIsoMu20Filter = i;
	  isIsoMu20Filter = true;
	}
	if (HLTFilter==HLTIsoMu18Filter) {
	  nIsoMu18Filter = i;
	  isIsoMu18Filter = true;
	}
	if (HLTFilter==HLTIsoMu22Filter) {
	  nIsoMu22Filter = i;
	  isIsoMu22Filter = true;
	}
	if (HLTFilter==HLTMu27Filter) {
	  nMu27Filter = i;
	  isMu27Filter = true;
	}
	if (HLTFilter==HLTMu45Eta2p1Filter) {
	  nMu45Eta2p1Filter = i;
	  isMu45Eta2p1Filter = true;
	}
      }
      if (!isHighPtLeg) {
	std::cout << "HLT filter " << HLTHighPtLeg << " not found" << std::endl;
	exit(-1);
      }
      if (!isLowPtLeg) {
	std::cout << "HLT filter " << HLTLowPtLeg1
		  << " or " << HLTLowPtLeg2
		  << " not found" << std::endl;
	exit(-1);
      }
      if (!isDZFilter) {
	std::cout << "HLT filter " << DZFilter << " not found" << std::endl;
	exit(-1);
      }
      if (!isSameSignFilter) {
	std::cout << "HLT filter " << SameSignFilter << " not found" << std::endl;
	exit(-1);
      }
      if (!isIsoSingleMuFilter) {
	std::cout << "HLT filter " << HLTIsoSingleMuFilter << " not found" << std::endl;
        exit(-1);
      }

      // vertex cuts

      prescaleHLTDoubleMuH->Fill(float(presHLTDoubleMu),weight);
      prescaleHLTDoubleMuDZH->Fill(float(presHLTDoubleMuDZ),weight);
      prescaleHLTDoubleMuSameSignDZH->Fill(float(presHLTDoubleMuSameSignDZ),weight);
      prescaleHLTIsoSingleMuH->Fill(float(presHLTIsoSingleMu),weight);
      prescaleHLTIsoMu18H->Fill(float(presHLTIsoMu18),weight);
      prescaleHLTIsoMu20H->Fill(float(presHLTIsoMu20),weight);
      prescaleHLTIsoMu22H->Fill(float(presHLTIsoMu22),weight);
      prescaleHLTMu27H->Fill(float(presHLTMu27),weight);
      prescaleHLTMu45Eta2p1H->Fill(float(presHLTMu45Eta2p1),weight);
      
      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isICHEP[im]) continue;
	float absIso = analysisTree.muon_chargedHadIso[im];
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut&&applyMuonIso) continue;
	muons.push_back(im);
      }

      if (muons.size()<2) continue;

      bool isPairSelected = false;
      bool isPairSelectedHLTDoubleMu = false;
      bool isPairSelectedHLTDoubleMuDZ = false;
      bool isPairSelectedHLTDoubleMuSameSignDZ = false;

      // selecting muon pair
      for (unsigned int im1=0; im1<muons.size()-1; ++im1) {
	//	  std::cout << "Muon " << im << std::endl;
	int  mu1Index = muons[im1];
	bool mu1MatchHighPt = false;
	bool mu1MatchLowPt  = false;
	bool mu1MatchMu45   = false;
	bool mu1MatchDz     = false;
	bool mu1MatchSS     = false;
	bool mu1MatchIsoSingleMu = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nHighPtLeg]) // Muon17 Leg
	    mu1MatchHighPt = true;
	  if (analysisTree.trigobject_filters[iT][nLowPtLeg]) // Muon8 Leg
	    mu1MatchLowPt = true;
	  if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
	    mu1MatchDz = true;
	  if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
	    mu1MatchSS = true;
	  if (analysisTree.trigobject_filters[iT][nIsoSingleMuFilter]) // HLT_IsoMu filter
	    mu1MatchIsoSingleMu = true;
	  if (analysisTree.trigobject_filters[iT][nMu45Eta2p1Filter]) // HLT_Mu45Eta2p1 filter
	    mu1MatchMu45 = true;
	}

	bool mu1HighPt = mu1MatchHighPt && analysisTree.muon_pt[mu1Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonHighCut;
	bool mu1LowPt  = mu1MatchLowPt && analysisTree.muon_pt[mu1Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonLowCut;
	bool mu1IsoSingleMu = mu1MatchIsoSingleMu && analysisTree.muon_pt[mu1Index]>ptMuonTagCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonTagCut;
	bool mu1Mu45 = mu1MatchMu45 && analysisTree.muon_pt[mu1Index]>46 && fabs(analysisTree.muon_eta[mu1Index])<2.1;

	float q1 = analysisTree.muon_charge[mu1Index];
	
	for (unsigned int im2=im1+1; im2<muons.size(); ++im2) {

	  int  mu2Index = muons[im2];

	  float q2 = analysisTree.muon_charge[mu2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index]);

	  bool mu2MatchHighPt = false;
	  bool mu2MatchLowPt  = false;
	  bool mu2MatchDz   = false;
	  bool mu2MatchSS   = false;
	  bool mu2MatchMu45   = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nHighPtLeg]) // Muon17 Leg
              mu2MatchHighPt = true;
	    if (analysisTree.trigobject_filters[iT][nLowPtLeg]) // Muon8 Leg
              mu2MatchLowPt = true;
	    if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
              mu2MatchDz = true;
	    if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
              mu2MatchSS = true;
	    if (analysisTree.trigobject_filters[iT][nMu45Eta2p1Filter]) // HLT_Mu45Eta2p1 filter
	      mu2MatchMu45 = true;

	  }
	  bool mu2HighPt = mu2MatchHighPt && analysisTree.muon_pt[mu2Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonHighCut;
	  bool mu2LowPt  = mu2MatchLowPt && analysisTree.muon_pt[mu2Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonLowCut;

	  bool mu2Mu45 = mu2MatchMu45 && analysisTree.muon_pt[mu2Index]>46 && fabs(analysisTree.muon_eta[mu2Index])<2.1;

	  bool triggerMatch = (mu1HighPt&&mu2LowPt) || (mu1LowPt&&mu2HighPt);
	  bool triggerMatchDz = triggerMatch && mu1MatchDz && mu2MatchDz;
	  bool triggerMatchSS = triggerMatchDz && mu1MatchSS && mu2MatchSS;

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

	  float absIso1 = analysisTree.muon_chargedHadIso[mu1Index];
	  float relIso1 = absIso1/analysisTree.muon_pt[mu1Index];

	  float absIso2 = analysisTree.muon_chargedHadIso[mu2Index];
	  float relIso2 = absIso2/analysisTree.muon_pt[mu2Index];

	  float mu1RelIso = relIso1;

	  if (analysisTree.muon_pt[mu2Index]>analysisTree.muon_pt[mu1Index]) {
	    float temp = relIso1;
	    relIso1 = relIso2;
	    relIso2 = temp;
	  }

	  bool dirIso = (relIso2<isoMuonCut) && (relIso1<isoMuonCut);

	  isPairSelected = true;
	  selPairs++;

	  float DZ = TMath::Min(float(dZ),float(0.5));
	  int DZBin  = binNumber(DZ,nDzBins,dzBins);
	  
	  //	  std::cout << "dZ = " << dZ << "  ->  " << DZBin << std::endl;

	  if (isHLTDoubleMu && triggerMatch) { // pass HLT_HighPt_LowPt
	    isPairSelectedHLTDoubleMu = true;
	    selPairsHLTDoubleMu++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTDoubleMuH->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTDoubleMuH->Fill(relIso2,weight);
	      dRmumuJPsiHLTDoubleMuH->Fill(dRmumu,weight);
	      dZmumuJPsiHLTDoubleMuH->Fill(dZ,weight);
	    }
	    if (triggerMatchDz) { // pass HLT_HighPt_LowPt_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterPassH->Fill(mass,weight);
		if (dRmumu>dRleptonsCut) ZMassDZFilterPassH->Fill(mass,weight);
	      }
	      JPsiMassDzPassH[DZBin]->Fill(mass,weight);
	      if (dZ<0.1) {
		JPsiMassDZFilterDz0to1PassH->Fill(mass,weight);
		if (dRmumu>dRleptonsCut) ZMassDZFilterDz0to1PassH->Fill(mass,weight);
	      }
	      else if (dZ<0.2) {
		JPsiMassDZFilterDz1to2PassH->Fill(mass,weight); 
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz1to2PassH->Fill(mass,weight);
	      }
	      else {
		JPsiMassDZFilterDzGt2PassH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDzGt2PassH->Fill(mass,weight);
	      }
	    }
	    else { // fail HLT_HighPt_LowPt_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterFailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterFailH->Fill(mass,weight);
              }
	      JPsiMassDzFailH[DZBin]->Fill(mass,weight);
	      if (dZ<0.1) {
                JPsiMassDZFilterDz0to1FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz0to1FailH->Fill(mass,weight);
              }
              else if (dZ<0.2) {
                JPsiMassDZFilterDz1to2FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz1to2FailH->Fill(mass,weight);
              }
              else {
                JPsiMassDZFilterDzGt2FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDzGt2FailH->Fill(mass,weight);
              }
	    }
	  }

	  if (isHLTDoubleMuDZ && triggerMatchDz) { // pass HLT_HighPt_LowPt_DZ
	    isPairSelectedHLTDoubleMuDZ = true;
            selPairsHLTDoubleMuDZ++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTDoubleMuDZH->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTDoubleMuDZH->Fill(relIso2,weight);
	      dRmumuJPsiHLTDoubleMuDZH->Fill(dRmumu,weight);
	      dZmumuJPsiHLTDoubleMuDZH->Fill(dZ,weight);
	    }
	    if (triggerMatchSS) { // pass HLT_HighPt_LowPt_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassSameSignFilterPassH->Fill(mass,weight);
		if (dirIso)
		  JPsiMassSameSignFilterDirIsoPassH->Fill(mass,weight);
		else 
		  JPsiMassSameSignFilterInvIsoPassH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassSameSignFilterPassH->Fill(mass,weight);
		if (dRmumu<0.15)
		  JPsiMassSameSignFilterDR0to0p15PassH->Fill(mass,weight);
		else
		  JPsiMassSameSignFilterDRGt0p15PassH->Fill(mass,weight);
              }
            }
            else { // fail HLT_HighPt_LowPt_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassSameSignFilterFailH->Fill(mass,weight);
		if (dirIso)
                  JPsiMassSameSignFilterDirIsoFailH->Fill(mass,weight);
                else
                  JPsiMassSameSignFilterInvIsoFailH->Fill(mass,weight); 
                if (dRmumu>dRleptonsCut) ZMassSameSignFilterFailH->Fill(mass,weight);
		if (dRmumu<0.15)
		  JPsiMassSameSignFilterDR0to0p15FailH->Fill(mass,weight);
		else
		  JPsiMassSameSignFilterDRGt0p15FailH->Fill(mass,weight);
              }
            }
	  }

	  if (isHLTDoubleMuSameSignDZ && triggerMatchSS) {
	    isPairSelectedHLTDoubleMuSameSignDZ = true;
            selPairsHLTDoubleMuSameSignDZ++;
	  }

	  if (isHLTIsoSingleMu && mu1IsoSingleMu && mu1RelIso<isoMuonCut && dZ<dZleptonsCut && dRmumu>dRleptonsCut) { // Single muon selection

	    float mu2AbsEta = fabs(analysisTree.muon_eta[mu2Index]);
	    float mu2Pt = TMath::Max(float(5.01),TMath::Min(float(analysisTree.muon_pt[mu2Index]),float(99.9)));
	    int etaBin = binNumber(mu2AbsEta,nEtaBins,etaBins);
	    int ptBin  = binNumber(mu2Pt,nPtBins,ptBins);
	    int ptBin45 = binNumber(mu2Pt,nPtBins45,ptBins45);

	    bool chargeTagPassed = true;
	    if (chargeTagMuon<0 && analysisTree.muon_charge[mu1Index]>0) chargeTagPassed = false;
	    if (chargeTagMuon>0 && analysisTree.muon_charge[mu1Index]<0) chargeTagPassed = false;

	    if (chargeTagPassed) {

	      if (mu2MatchLowPt) { // LowPt Leg
		if (analysisTree.muon_pt[mu2Index]>ptMuonLowCut&&mu2AbsEta<etaMuonLowCut) 
		  ZMassLowPtLegPassH->Fill(mass,weight);
		if (analysisTree.muon_pt[mu2Index]>ptMuonLowCut)
		  ZMassLowPtLegEtaPassH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegPtPassH[ptBin]->Fill(mass,weight);
	      }
	      else {
		if (analysisTree.muon_pt[mu2Index]>ptMuonLowCut&&mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegFailH->Fill(mass,weight);
		if (analysisTree.muon_pt[mu2Index]>ptMuonLowCut)
		  ZMassLowPtLegEtaFailH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegPtFailH[ptBin]->Fill(mass,weight);
	      }
	     
	      if (mu2MatchHighPt) { // HighPt Leg
		if (analysisTree.muon_pt[mu2Index]>ptMuonHighCut&&mu2AbsEta<etaMuonHighCut)  
		  ZMassHighPtLegPassH->Fill(mass,weight);
		if (analysisTree.muon_pt[mu2Index]>ptMuonHighCut) 
		  ZMassHighPtLegEtaPassH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegPtPassH[ptBin]->Fill(mass,weight);
	      }
	      else {
		if (analysisTree.muon_pt[mu2Index]>ptMuonHighCut&&mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegFailH->Fill(mass,weight);
		if (analysisTree.muon_pt[mu2Index]>ptMuonHighCut) 
		  ZMassHighPtLegEtaFailH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegPtFailH[ptBin]->Fill(mass,weight);
	      }

	      if (mu2AbsEta<2.4) { // bins in (eta,pt)

		if (mu2MatchLowPt)
		  ZMassLowPtLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassLowPtLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

		if (mu2MatchHighPt)
		  ZMassHighPtLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassHighPtLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

	      }

	      if (mu2MatchMu45) { // Mu45 Leg
		if (analysisTree.muon_pt[mu2Index]>46) 
		  ZMassMu45LegEtaPassH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<2.1)
		  ZMassMu45LegPtPassH[ptBin45]->Fill(mass,weight);
	      }
	      else {
		if (analysisTree.muon_pt[mu2Index]>46) 
		  ZMassMu45LegEtaFailH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<2.1)
		  ZMassMu45LegPtFailH[ptBin45]->Fill(mass,weight);
	      }


	      if (mu2AbsEta<2.1) { // Mu45 Leg bins in (eta,pt)

		if (mu2MatchMu45)
                  ZMassMu45LegPtEtaPassH[etaBin][ptBin45]->Fill(mass,weight);
                else
                  ZMassMu45LegPtEtaFailH[etaBin][ptBin45]->Fill(mass,weight);

	      }

	    }

	  }

	}
      }
    
      if (isPairSelected)
	selEvents++;
      if (isPairSelectedHLTDoubleMu)
	selEventsHLTDoubleMu++;
      if (isPairSelectedHLTDoubleMuDZ)
	selEventsHLTDoubleMuDZ++;
      if (isPairSelectedHLTDoubleMuSameSignDZ)
	selEventsHLTDoubleMuSameSignDZ++;

      
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
  std::cout << "Total number of selected events (HLT_Dimuon) = " << selEventsHLTDoubleMu << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon) = " << selPairsHLTDoubleMu  << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Dimuon_DZ) = " << selEventsHLTDoubleMuDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon_DZ) = " << selPairsHLTDoubleMuDZ << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Dimuon_SameSign_DZ) = " << selEventsHLTDoubleMuSameSignDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon_SameSign_DZ) = " << selPairsHLTDoubleMuSameSignDZ << std::endl;
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



