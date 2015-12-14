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

  const string hltMu17Mu8           = cfg.get<string>("HLTMu17Mu8");
  const string hltMu17Mu8DZ         = cfg.get<string>("HLTMu17Mu8DZ");
  const string hltMu17Mu8SameSignDZ = cfg.get<string>("HLTMu17Mu8SameSignDZ");

  const string hltIsoMu24Eta2p1     = cfg.get<string>("HLTIsoMu24Eta2p1");  
  const string hltMu24              = cfg.get<string>("HLTMu24");  
  const string hltMu24Eta2p1        = cfg.get<string>("HLTMu24Eta2p1");  
  const string hltMu27              = cfg.get<string>("HLTMu27");  
  const string hltMu34              = cfg.get<string>("HLTMu34");  
  const string hltMu45Eta2p1        = cfg.get<string>("HLTMu45Eta2p1");  


  // HLT filters
  const string hltMu17Leg     = cfg.get<string>("HLTMu17Leg");
  const string hltMu8Leg      = cfg.get<string>("HLTMu8Leg");
  const string dzFilter       = cfg.get<string>("dzFilter"); 
  const string sameSignFilter = cfg.get<string>("sameSignFilter");

  const string hltIsoMu24Eta2p1Filter = cfg.get<string>("HLTIsoMu24Eta2p1Filter");
  const string hltMu24Filter          = cfg.get<string>("HLTMu24Filter");
  const string hltMu24Eta2p1Filter    = cfg.get<string>("HLTMu24Eta2p1Filter");
  const string hltMu27Filter          = cfg.get<string>("HLTMu27Filter");
  const string hltMu34Filter          = cfg.get<string>("HLTMu34Filter");
  const string hltMu45Eta2p1Filter    = cfg.get<string>("HLTMu45Eta2p1Filter");

  // convesrsion from string to TString
  TString HLTMu17Mu8(hltMu17Mu8);
  TString HLTMu17Mu8DZ(hltMu17Mu8DZ);
  TString HLTMu17Mu8SameSignDZ(hltMu17Mu8SameSignDZ);

  TString HLTIsoMu24Eta2p1(hltIsoMu24Eta2p1);
  TString HLTMu24(hltMu24);
  TString HLTMu24Eta2p1(hltMu24Eta2p1);
  TString HLTMu27(hltMu27);
  TString HLTMu34(hltMu34);
  TString HLTMu45Eta2p1(hltMu45Eta2p1);

  TString HLTMu17Leg(hltMu17Leg);
  TString HLTMu8Leg(hltMu8Leg);
  TString DZFilter(dzFilter);
  TString SameSignFilter(sameSignFilter);

  TString HLTIsoMu24Eta2p1Filter(hltIsoMu24Eta2p1Filter);
  TString HLTMu24Filter(hltMu24Filter);
  TString HLTMu24Eta2p1Filter(hltMu24Eta2p1Filter);
  TString HLTMu27Filter(hltMu27Filter);
  TString HLTMu34Filter(hltMu34Filter);
  TString HLTMu45Eta2p1Filter(hltMu45Eta2p1Filter);
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");


  // **** end of configuration

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Run-lumi selector

  std::vector<Period> periods;
    
  std::fstream inputFileStream("temp_muon", std::ios::in);
  for(std::string s; std::getline(inputFileStream, s); )
    {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  int presHLTMu17Mu8 = 0;
  int presHLTMu17Mu8DZ = 0;
  int presHLTMu17Mu8SameSignDZ = 0;
  int presHLTIsoMu24Eta2p1 = 0;
  int presHLTMu24 = 0;
  int presHLTMu24Eta2p1 = 0;
  int presHLTMu27 = 0;
  int presHLTMu34 = 0;
  int presHLTMu45Eta2p1 = 0;
  unsigned int eventNumber = 0;
  unsigned int runNumber = 0;
  unsigned int lumiBlock = 0;

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");
  TTree * treePrescales = new TTree("Prescales","Prescales");
  treePrescales->Branch("Run",&runNumber,"Run/i");
  treePrescales->Branch("Event",&eventNumber,"Event/i");
  treePrescales->Branch("Lumi",&lumiBlock,"Run/i");
  treePrescales->Branch("HLTMu17Mu8",&presHLTMu17Mu8,"HLTMu17Mu8/I");
  treePrescales->Branch("HLTMu17Mu8DZ",&presHLTMu17Mu8DZ,"HLTMu17Mu8DZ/I");
  treePrescales->Branch("HLTMu17Mu8SameSignDZ",&presHLTMu17Mu8SameSignDZ,"HLTMu17Mu8SameSignDZ/I");
  treePrescales->Branch("HLTIsoMu24Eta2p1",&presHLTIsoMu24Eta2p1,"HLTIsoMu24Eta2p1/I");
  treePrescales->Branch("HLTMu24Eta2p1",&presHLTMu24Eta2p1,"HLTMu24Eta2p1/I");
  treePrescales->Branch("HLTMu24",&presHLTMu24,"HLTMu24/I");
  treePrescales->Branch("HLTMu27",&presHLTMu27,"HLTMu27/I");
  treePrescales->Branch("HLTMu34",&presHLTMu34,"HLTMu34/I");
  treePrescales->Branch("HLTMu45Eta2p1",&presHLTMu45Eta2p1,"HLTMu45Eta2p1/I");

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);

  // prescales
  TH1F * prescaleHLTMu17Mu8H = new TH1F("prescaleHLTMu17Mu8H","",21,-0.5,20.5);
  TH1F * prescaleHLTMu17Mu8DZH = new TH1F("prescaleHLTMu17Mu8DZH","",21,-0.5,20.5);
  TH1F * prescaleHLTMu17Mu8SameSignDZH = new TH1F("prescaleHLTMu17Mu8SameSignDZH","",21,-0.5,20.5);
  TH1F * prescaleHLTIsoMu24Eta2p1H = new TH1F("prescaleHLTIsoMu24Eta2p1H","",21,-0.5,20.5);
  TH1F * prescaleHLTMu24H = new TH1F("prescaleHLTMu24H","",101,-0.5,100.5);
  TH1F * prescaleHLTMu24Eta2p1H = new TH1F("prescaleHLTMu24Eta2p1H","",101,-0.5,100.5);
  TH1F * prescaleHLTMu27H = new TH1F("prescaleHLTMu27H","",101,-0.5,100.5);
  TH1F * prescaleHLTMu34H = new TH1F("prescaleHLTMu34H","",101,-0.5,100.5);
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

  TH1F * muIsoLeadJPsiHLTMu17Mu8H   = new TH1F("muIsoLeadJPsiHLTMu17Mu8H","",200,0,2);
  TH1F * muIsoLeadJPsiHLTMu17Mu8DZH = new TH1F("muIsoLeadJPsiHLTMu17Mu8DZH","",200,0,2);

  TH1F * muIsoTrailJPsiHLTMu17Mu8H   = new TH1F("muIsoTrailJPsiHLTMu17Mu8H","",200,0,2);
  TH1F * muIsoTrailJPsiHLTMu17Mu8DZH = new TH1F("muIsoTrailJPsiHLTMu17Mu8DZH","",200,0,2);

  TH1F * dRmumuJPsiHLTMu17Mu8H   = new TH1F("dRmumuJPsiHLTMu17Mu8H","",200,0,2);
  TH1F * dRmumuJPsiHLTMu17Mu8DZH = new TH1F("dRmumuJPsiHLTMu17Mu8DZH","",200,0,2);

  TH1F * dZmumuJPsiHLTMu17Mu8H   = new TH1F("dZmumuJPsiHLTMu17Mu8H","",100,0,1);
  TH1F * dZmumuJPsiHLTMu17Mu8DZH = new TH1F("dZmumuJPsiHLTMu17Mu8DZH","",100,0,1);

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

  int nEtaBins = 2;
  float etaBins[3] = {-0.1, 1.48, 2.4};
  int nPtBins = 6;
  float ptBins[7] = {5, 10, 15, 20, 25, 30, 40};

  TH1F * ZMassMu17LegPassH = new TH1F("ZMassMu17LegPassH","",60,60,120);
  TH1F * ZMassMu17LegFailH = new TH1F("ZMassMu17LegFailH","",60,60,120);
  TH1F * ZMassMu8LegPassH = new TH1F("ZMassMu8LegPassH","",60,60,120);
  TH1F * ZMassMu8LegFailH = new TH1F("ZMassMu8LegFailH","",60,60,120);

  TH1F * ZMassMu17LegBarrelPassH = new TH1F("ZMassMu17LegBarrelPassH","",60,60,120);
  TH1F * ZMassMu17LegBarrelFailH = new TH1F("ZMassMu17LegBarrelFailH","",60,60,120);
  TH1F * ZMassMu8LegBarrelPassH = new TH1F("ZMassMu8LegBarrelPassH","",60,60,120);
  TH1F * ZMassMu8LegBarrelFailH = new TH1F("ZMassMu8LegBarrelFailH","",60,60,120);

  TH1F * ZMassMu17LegEndcapPassH = new TH1F("ZMassMu17LegEndcapPassH","",60,60,120);
  TH1F * ZMassMu17LegEndcapFailH = new TH1F("ZMassMu17LegEndcapFailH","",60,60,120);
  TH1F * ZMassMu8LegEndcapPassH = new TH1F("ZMassMu8LegEndcapPassH","",60,60,120);
  TH1F * ZMassMu8LegEndcapFailH = new TH1F("ZMassMu8LegEndcapFailH","",60,60,120);



  TString EtaBins[2] = {"Barrel","Endcap"};
  TString PtBins[6]  = {"Pt5to10","Pt10to15","Pt15to20","Pt20to25","Pt25to30","Pt30toInf"};

  TH1F * ZMassMu17LegPtEtaPassH[2][6];
  TH1F * ZMassMu17LegPtEtaFailH[2][6];

  TH1F * ZMassMu8LegPtEtaPassH[2][6];
  TH1F * ZMassMu8LegPtEtaFailH[2][6];

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassMu17LegPtEtaPassH[iEta][iPt] = new TH1F("ZMassMu17Leg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassMu17LegPtEtaFailH[iEta][iPt] = new TH1F("ZMassMu17Leg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
      ZMassMu8LegPtEtaPassH[iEta][iPt]  = new TH1F("ZMassMu8Leg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassMu8LegPtEtaFailH[iEta][iPt]  = new TH1F("ZMassMu8Leg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
    }
  }

  unsigned int iRun;
  unsigned int iEvent;
  TTree * eventTree = new TTree("eventTree","eventTree");
  eventTree->Branch("Run",&iRun,"Run/i");
  eventTree->Branch("Event",&iEvent,"Event/i");

  int nFiles = 0;
  int nEvents = 0;

  int selEvents = 0;
  int selEventsHLTMu17Mu8 = 0;
  int selEventsHLTMu17Mu8DZ = 0;
  int selEventsHLTMu17Mu8SameSignDZ = 0;

  int selPairs = 0;
  int selPairsHLTMu17Mu8 = 0;
  int selPairsHLTMu17Mu8DZ = 0;
  int selPairsHLTMu17Mu8SameSignDZ = 0;

  unsigned int minRun = 99999999;
  unsigned int maxRun = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
  std::vector<unsigned int> allGoodRuns; allGoodRuns.clear();

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

      // triggers
      bool isHLTMu17Mu8 = false;
      bool isHLTMu17Mu8DZ = false;
      bool isHLTMu17Mu8SameSignDZ = false;
      bool isHLTIsoMu24Eta2p1 = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTMu17Mu8)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
	    isHLTMu17Mu8 = true;
	}
	if (trigName.Contains(HLTMu17Mu8DZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTMu17Mu8DZ = true;
	}
	if (trigName.Contains(HLTMu17Mu8SameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTMu17Mu8SameSignDZ = true;
	}
	if (trigName.Contains(HLTIsoMu24Eta2p1)) {
	  if (it->second==1)
            isHLTIsoMu24Eta2p1 = true;
	}
      }

      presHLTMu17Mu8 = 0;
      presHLTMu17Mu8DZ = 0;
      presHLTMu17Mu8SameSignDZ = 0;
      presHLTIsoMu24Eta2p1 = 0;
      presHLTMu24 = 0;
      presHLTMu24Eta2p1 = 0;
      presHLTMu27 = 0;
      presHLTMu34 = 0;
      presHLTMu45Eta2p1 = 0;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerprescales->begin(); it!=analysisTree.hltriggerprescales->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTMu17Mu8)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTMu17Mu8 = it->second;
	}
	if (trigName.Contains(HLTMu17Mu8DZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTMu17Mu8DZ = it->second;
	}
	if (trigName.Contains(HLTMu17Mu8SameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTMu17Mu8SameSignDZ = it->second;
	}
	if (trigName.Contains(HLTIsoMu24Eta2p1)) {
	  presHLTIsoMu24Eta2p1 = it->second;
	}
	if (trigName.Contains(HLTMu24)) {
	  presHLTMu24 = it->second;
	  //	  if (presHLTMu24>100)
	  //	    cout << "HLTMu24 prescale : " << presHLTMu24 << std::endl;
	}
	if (trigName.Contains(HLTMu24Eta2p1)) {
	  presHLTMu24Eta2p1 = it->second;
	  // if (presHLTMu24Eta2p1>100)
	  //   cout << "HLTMu24Eta2p1 prescale : " << presHLTMu24Eta2p1 << std::endl;
	}
	if (trigName.Contains(HLTMu27)) {
	  presHLTMu27 = it->second;
	  // if (presHLTMu27>100)
	  //   cout << "HLTMu27 prescale : " << presHLTMu27 << std::endl;
	}
	if (trigName.Contains(HLTMu34)) {
	  presHLTMu34 = it->second;
	  // if (presHLTMu34>100)
	  //   cout << "HLTMu34 prescale : " << presHLTMu34 << std::endl;
	}
	if (trigName.Contains(HLTMu45Eta2p1)) {
	  presHLTMu45Eta2p1 = it->second;
	}
      }
      
      runNumber = analysisTree.event_run;
      eventNumber = analysisTree.event_nr;
      lumiBlock = analysisTree.event_luminosityblock;
      treePrescales->Fill();

      if (!applyTrigger) {
	isHLTMu17Mu8 = true;
	isHLTMu17Mu8DZ = true;
	isHLTMu17Mu8SameSignDZ = true;
	isHLTIsoMu24Eta2p1 = true;
      }

      unsigned int nMu17Leg = 0;
      bool isMu17Leg = false;
      unsigned int nMu8Leg = 0;
      bool isMu8Leg = false;
      unsigned int nDZFilter = 0;
      bool isDZFilter = false;
      unsigned int nSameSignFilter = 0;
      bool isSameSignFilter = false;

      unsigned int nIsoMu24Eta2p1Filter = 0;
      bool isIsoMu24Eta2p1Filter;
      unsigned int nMu24Eta2p1Filter = 0;
      bool isMu24Eta2p1Filter;
      unsigned int nMu24Filter = 0;
      bool isMu24Filter;
      unsigned int nMu27Filter = 0;
      bool isMu27Filter;
      unsigned int nMu34Filter = 0;
      bool isMu34Filter;
      unsigned int nMu45Eta2p1Filter = 0;
      bool isMu45Eta2p1Filter;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==HLTMu17Leg) {
	  nMu17Leg = i;
	  isMu17Leg = true;
	}
	if (HLTFilter==HLTMu8Leg) {
	  nMu8Leg = i;
	  isMu8Leg = true;
	}
	if (HLTFilter==DZFilter) {
	  nDZFilter = i;
	  isDZFilter = true;
	}
	if (HLTFilter==SameSignFilter) {
	  nSameSignFilter = i;
	  isSameSignFilter = true;
	}
	if (HLTFilter==HLTIsoMu24Eta2p1Filter) {
	  nIsoMu24Eta2p1Filter = i;
	  isIsoMu24Eta2p1Filter = true;
	}
	if (HLTFilter==HLTMu24Eta2p1Filter) {
	  nMu24Eta2p1Filter = i;
	  isMu24Eta2p1Filter = true;
	}
	if (HLTFilter==HLTMu24Filter) {
	  nMu24Filter = i;
	  isMu24Filter = true;
	}
	if (HLTFilter==HLTMu27Filter) {
	  nMu27Filter = i;
	  isMu27Filter = true;
	}
	if (HLTFilter==HLTMu34Filter) {
	  nMu34Filter = i;
	  isMu34Filter = true;
	}
	if (HLTFilter==HLTMu45Eta2p1Filter) {
	  nMu45Eta2p1Filter = i;
	  isMu45Eta2p1Filter = true;
	}
      }
      if (!isMu17Leg) {
	std::cout << "HLT filter " << HLTMu17Leg << " not found" << std::endl;
	exit(-1);
      }
      if (!isMu8Leg) {
	std::cout << "HLT filter " << HLTMu8Leg << " not found" << std::endl;
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

      // vertex cuts
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;

      prescaleHLTMu17Mu8H->Fill(float(presHLTMu17Mu8),weight);
      prescaleHLTMu17Mu8DZH->Fill(float(presHLTMu17Mu8DZ),weight);
      prescaleHLTMu17Mu8SameSignDZH->Fill(float(presHLTMu17Mu8SameSignDZ),weight);
      prescaleHLTIsoMu24Eta2p1H->Fill(float(presHLTIsoMu24Eta2p1),weight);
      prescaleHLTMu24H->Fill(float(presHLTMu24),weight);
      prescaleHLTMu24Eta2p1H->Fill(float(presHLTMu24Eta2p1),weight);
      prescaleHLTMu27H->Fill(float(presHLTMu27),weight);
      prescaleHLTMu34H->Fill(float(presHLTMu34),weight);
      prescaleHLTMu45Eta2p1H->Fill(float(presHLTMu45Eta2p1),weight);
      
      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	float absIso = analysisTree.muon_chargedHadIso[im];
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut&&applyMuonIso) continue;
	muons.push_back(im);
      }

      if (muons.size()<2) continue;

      bool isPairSelected = false;
      bool isPairSelectedHLTMu17Mu8 = false;
      bool isPairSelectedHLTMu17Mu8DZ = false;
      bool isPairSelectedHLTMu17Mu8SameSignDZ = false;

      // selecting muon pair
      for (unsigned int im1=0; im1<muons.size()-1; ++im1) {
	//	  std::cout << "Muon " << im << std::endl;
	int  mu1Index = muons[im1];
	bool mu1MatchMu17 = false;
	bool mu1MatchMu8  = false;
	bool mu1MatchDz   = false;
	bool mu1MatchSS   = false;
	bool mu1MatchIsoMu24Eta2p1 = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nMu17Leg]) // Muon17 Leg
	    mu1MatchMu17 = true;
	  if (analysisTree.trigobject_filters[iT][nMu8Leg]) // Muon8 Leg
	    mu1MatchMu8 = true;
	  if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
	    mu1MatchDz = true;
	  if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
	    mu1MatchSS = true;
	  if (analysisTree.trigobject_filters[iT][nIsoMu24Eta2p1Filter]) // HLT_IsoMu24 filter
	    mu1MatchIsoMu24Eta2p1 = true;
	}

	bool mu1Mu17 = mu1MatchMu17 && analysisTree.muon_pt[mu1Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonHighCut;
	bool mu1Mu8  = mu1MatchMu8 && analysisTree.muon_pt[mu1Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonLowCut;
	bool mu1IsoMu24Eta2p1 = mu1MatchIsoMu24Eta2p1 && analysisTree.muon_pt[mu1Index]>ptMuonTagCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonTagCut;

	float q1 = analysisTree.muon_charge[mu1Index];
	
	for (unsigned int im2=im1+1; im2<muons.size(); ++im2) {

	  int  mu2Index = muons[im2];

	  float q2 = analysisTree.muon_charge[mu2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index]);

	  bool mu2MatchMu17 = false;
	  bool mu2MatchMu8  = false;
	  bool mu2MatchDz   = false;
	  bool mu2MatchSS   = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMu17Leg]) // Muon17 Leg
              mu2MatchMu17 = true;
	    if (analysisTree.trigobject_filters[iT][nMu8Leg]) // Muon8 Leg
              mu2MatchMu8 = true;
	    if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
              mu2MatchDz = true;
	    if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
              mu2MatchSS = true;

	  }
	  bool mu2Mu17 = mu2MatchMu17 && analysisTree.muon_pt[mu2Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonHighCut;
	  bool mu2Mu8  = mu2MatchMu8 && analysisTree.muon_pt[mu2Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonLowCut;

	  bool triggerMatch = (mu1Mu17&&mu2Mu8) || (mu1Mu8&&mu2Mu17);
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
	  
	  if (isHLTMu17Mu8 && triggerMatch) { // pass HLT_Mu17_Mu8
	    isPairSelectedHLTMu17Mu8 = true;
	    selPairsHLTMu17Mu8++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTMu17Mu8H->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTMu17Mu8H->Fill(relIso2,weight);
	      dRmumuJPsiHLTMu17Mu8H->Fill(dRmumu,weight);
	      dZmumuJPsiHLTMu17Mu8H->Fill(dZ,weight);
	    }
	    if (triggerMatchDz) { // pass HLT_Mu17_Mu8_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterPassH->Fill(mass,weight);
		if (dRmumu>dRleptonsCut) ZMassDZFilterPassH->Fill(mass,weight);
	      }
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
	    else { // fail HLT_Mu17_Mu8_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterFailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterFailH->Fill(mass,weight);
              }
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

	  if (isHLTMu17Mu8DZ && triggerMatchDz) { // pass HLT_Mu17_Mu8_DZ
	    isPairSelectedHLTMu17Mu8DZ = true;
            selPairsHLTMu17Mu8DZ++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTMu17Mu8DZH->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTMu17Mu8DZH->Fill(relIso2,weight);
	      dRmumuJPsiHLTMu17Mu8DZH->Fill(dRmumu,weight);
	      dZmumuJPsiHLTMu17Mu8DZH->Fill(dZ,weight);
	    }
	    if (triggerMatchSS) { // pass HLT_Mu17_Mu8_SameSign_DZ
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
            else { // fail HLT_Mu17_Mu8_SameSign_DZ
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

	  if (isHLTMu17Mu8SameSignDZ && triggerMatchSS) {
	    isPairSelectedHLTMu17Mu8SameSignDZ = true;
            selPairsHLTMu17Mu8SameSignDZ++;
	  }

	  if (isHLTIsoMu24Eta2p1 && mu1IsoMu24Eta2p1 && mu1RelIso<isoMuonCut && dZ<dZleptonsCut && dRmumu>dRleptonsCut) { // Single muon selection

	    float mu2AbsEta = fabs(analysisTree.muon_eta[mu2Index]);
	    float mu2Pt = TMath::Max(float(5.),TMath::Min(float(analysisTree.muon_pt[mu2Index]),float(39.9)));
	    int etaBin = binNumber(mu2AbsEta,nEtaBins,etaBins);
	    int ptBin  = binNumber(mu2Pt,nPtBins,ptBins);
	    bool chargeTagPassed = true;
	    if (chargeTagMuon<0 && analysisTree.muon_charge[mu1Index]>0) chargeTagPassed = false;
	    if (chargeTagMuon>0 && analysisTree.muon_charge[mu1Index]<0) chargeTagPassed = false;

	    if (chargeTagPassed) {

	      if (analysisTree.muon_pt[mu2Index]>ptMuonLowCut) { // Mu8 Leg
		if (mu2MatchMu8) {
		  ZMassMu8LegPassH->Fill(mass,weight);
		  if (mu2AbsEta<1.48) 
		    ZMassMu8LegBarrelPassH->Fill(mass,weight);
		  else if (mu2AbsEta<2.4)
		    ZMassMu8LegEndcapPassH->Fill(mass,weight);
		}
		else {
		  ZMassMu8LegFailH->Fill(mass,weight);
		  if (mu2AbsEta<1.48) 
		    ZMassMu8LegBarrelFailH->Fill(mass,weight);
		  else if (mu2AbsEta<2.4)
		    ZMassMu8LegEndcapFailH->Fill(mass,weight);
		}
	      }
	     
	      if (analysisTree.muon_pt[mu2Index]>ptMuonHighCut) { // Mu17 Leg
		if (mu2MatchMu17) {
		  ZMassMu17LegPassH->Fill(mass,weight);
		  if (mu2AbsEta<1.48) 
		    ZMassMu17LegBarrelPassH->Fill(mass,weight);
		  else if (mu2AbsEta<2.4)
		    ZMassMu17LegEndcapPassH->Fill(mass,weight);
		}
		else {
		  ZMassMu17LegFailH->Fill(mass,weight);
		  if (mu2AbsEta<1.48) 
		    ZMassMu17LegBarrelFailH->Fill(mass,weight);
		  else if (mu2AbsEta<2.4)
		    ZMassMu17LegEndcapFailH->Fill(mass,weight);
		}
	      }

	      if (mu2AbsEta<2.4) { // bins in eta,pt

		if (mu2MatchMu8)
		  ZMassMu8LegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassMu8LegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

		if (mu2MatchMu17)
		  ZMassMu17LegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassMu17LegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

	      }

	    }

	  }

	}
      }
    
      if (isPairSelected)
	selEvents++;
      if (isPairSelectedHLTMu17Mu8)
	selEventsHLTMu17Mu8++;
      if (isPairSelectedHLTMu17Mu8DZ)
	selEventsHLTMu17Mu8DZ++;
      if (isPairSelectedHLTMu17Mu8SameSignDZ)
	selEventsHLTMu17Mu8SameSignDZ++;

      
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
  std::cout << "Total number of selected events (HLT_Mu17_Mu8) = " << selEventsHLTMu17Mu8 << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8) = " << selPairsHLTMu17Mu8  << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Mu17_Mu8_DZ) = " << selEventsHLTMu17Mu8DZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8_DZ) = " << selPairsHLTMu17Mu8DZ << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Mu17_Mu8_SameSign_DZ) = " << selEventsHLTMu17Mu8SameSignDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8_SameSign_DZ) = " << selPairsHLTMu17Mu8SameSignDZ << std::endl;
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



