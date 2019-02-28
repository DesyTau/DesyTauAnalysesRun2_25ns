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
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"


int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // Data
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection"); 

  // kinematic cuts on elecs
  const float ptElecCut      = cfg.get<float>("ptElecCut");
  const float ptElecHighCut  = cfg.get<float>("ptElecHighCut");
  const float ptElecTagCut   = cfg.get<float>("ptElecTagCut");

  const float etaElecCut = cfg.get<float>("etaElecCut");
  const float etaElecHighCut = cfg.get<float>("etaElecHighCut");
  const float etaElecTagCut  = cfg.get<float>("etaElecTagCut");

  const float dxyElecCut     = cfg.get<float>("dxyElecCut");
  const float dzElecCut      = cfg.get<float>("dzElecCut");
  const float isoElecCut     = cfg.get<float>("isoElecCut");
  const int  electronIdType  = cfg.get<int>("ElectronIdType");
  const bool applyRhoCorrectedIso = cfg.get<bool>("ApplyRhoCorrectedIso");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZleptonsCut   = cfg.get<float>("dZleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");
  const float ptL1TauCut     = cfg.get<float>("ptL1TauCut");
  const float etaL1TauCut    = cfg.get<float>("etaL1TauCut");
  const float detaL1TauCut   = cfg.get<float>("detaL1TauCut");
  const bool matchL1Tau      = cfg.get<bool>("matchL1Tau");

  const bool matchL1Elec     = cfg.get<bool>("matchL1Elec");
  const float ptL1ElecCut     = cfg.get<float>("ptL1ElecCut");
  const float etaL1ElecCut    = cfg.get<float>("etaL1ElecCut");


  // triggers
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");

  const string hltSingleEle  = cfg.get<string>("HLTSingleEle");  

  // HLT filters
  const string hltSingleEleFilter = cfg.get<string>("HLTSingleEleFilter");
  const string hltEleProbeFilter  = cfg.get<string>("HLTEleProbeFilter");
  const string hltEleProbeFilter2 = cfg.get<string>("HLTEleProbeFilter2");

  const unsigned int runRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int runRangeMax = cfg.get<unsigned int>("RunRangeMax");

  const string jsonFile = cfg.get<string>("jsonFile");
  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");
  const string puMCHist = cfg.get<string>("PileUpMCHist");

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);
  TString PUMCHist(puMCHist);

  TString HLTSingleEle(hltSingleEle);

  TString HLTSingleEleFilter(hltSingleEleFilter);
  TString HLTEleProbeFilter(hltEleProbeFilter);
  TString HLTEleProbeFilter2(hltEleProbeFilter2);
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  const float isoCone = 0.3;

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

  bool isHLTSingleEle = 0;

  unsigned int eventNumber = 0;
  unsigned int runNumber = 0;
  unsigned int lumiBlock = 0;

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);

  int nEtaBins = 5;
  float etaBins[6] = {0, 1.0, 1.479, 1.653, 2.1, 2.5};


  int nEtaFineBins  = 10;
  float etaFineBins[11] = {-2.5, -2.1, -1.57, -1.44, -0.8, 0., 0.8, 1.44, 1.57, 2.1, 2.5};
  
  int nPtBins = 18;
  float ptBins[19] = {10, 12, 14, 16, 18,
                      20, 22, 24, 26, 28,
                      30, 32, 34, 36, 38,
                      40, 45, 50, 100};
  
  //  int nPtBins = 7;
  //  float ptBins[8] = {10, 15, 20, 25, 30, 35, 40, 100}; 


  TH1F * ZMassSingleEleLegPassH = new TH1F("ZMassSingleEleLegPassH","",60,60,120);
  TH1F * ZMassSingleEleLegFailH = new TH1F("ZMassSingleEleLegFailH","",60,60,120);


  TString EtaBins[5] = {"EtaLt1p0",
                        "Eta1p0to1p48",
                        "Eta1p48to1p65",
                        "Eta1p65to2p1",
                        "EtaGt2p1"};

  TString EtaFineBins[10];
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

  //  TString PtBins[7] = {"Pt10to15","Pt15to20","Pt20to25","Pt25to30","Pt30to35","Pt35to40","PtGt40"};
  TString PtBins[18] = {"Pt10to12","Pt12to14","Pt14to16","Pt16to18","Pt18to20",
			"Pt20to22","Pt22to24","Pt24to26","Pt26to28","Pt28to30",
			"Pt30to32","Pt32to34","Pt34to36","Pt36to38","Pt38to40",
			"Pt40to45","Pt45to50","PtGt50"};


  TH1F * EtaBinsH = new TH1F("EtaBinsH","",nEtaBins,etaBins);
  TH1F * EtaFineBinsH = new TH1F("EtaFineBinsH","",nEtaFineBins,etaFineBins);
  TH1F * PtBinsH  = new TH1F("PtBinsH","",nPtBins,ptBins);

  for (int iB=0; iB<nEtaFineBins; ++iB)
    EtaFineBinsH->GetXaxis()->SetBinLabel(iB+1,EtaFineBins[iB]);

  for (int iB=0; iB<nEtaBins; ++iB)
    EtaBinsH->GetXaxis()->SetBinLabel(iB+1,EtaBins[iB]);

  for (int iB=0; iB<nPtBins; ++iB)
    PtBinsH->GetXaxis()->SetBinLabel(iB+1,PtBins[iB]);

  int nIsoBins = 4;
  float isoBins[5] = {-0.01,0.10,0.2,0.5,1.0};
  TString IsoBins[4] = {"IsoLt0p1",
			"Iso0p1to0p2",
			"Iso0p2to0p5",
			"IsoGt0p5"};

  TH1F * isoBinsH = new TH1F("isoBinsH","",nIsoBins,isoBins);
  for (int iB=0; iB<nIsoBins; ++iB)
    isoBinsH->GetXaxis()->SetBinLabel(iB+1,IsoBins[iB]);

  // (Pt,Eta)

  TH1F * ZMassSingleEleLegPtEtaPassH[5][18];
  TH1F * ZMassSingleEleLegPtEtaFailH[5][18];

  TH1F * ZMassSingleEleLegIsoPtEtaPassH[4][5][18];
  TH1F * ZMassSingleEleLegIsoPtEtaFailH[4][5][18];

  TH1F * ZMassSingleEleLegEtaPassH[5];
  TH1F * ZMassSingleEleLegEtaFailH[5];

  TH1F * ZMassSingleEleLegEtaFinePassH[10];
  TH1F * ZMassSingleEleLegEtaFineFailH[10];

  TH1F * ZMassSingleEleLegPtPassH[18];
  TH1F * ZMassSingleEleLegPtFailH[18];


  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    ZMassSingleEleLegEtaPassH[iEta] = new TH1F("ZMassSingleEleLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassSingleEleLegEtaFailH[iEta] = new TH1F("ZMassSingleEleLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaFineBins; ++iEta) {
    ZMassSingleEleLegEtaFinePassH[iEta] = new TH1F("ZMassSingleEleLeg_"+EtaFineBins[iEta]+"_PassH","",60,60,120);
    ZMassSingleEleLegEtaFineFailH[iEta] = new TH1F("ZMassSingleEleLeg_"+EtaFineBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassSingleEleLegPtPassH[iPt] = new TH1F("ZMassSingleEleLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassSingleEleLegPtFailH[iPt] = new TH1F("ZMassSingleEleLeg_"+PtBins[iPt]+"_FailH","",60,60,120);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {

    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassSingleEleLegPtEtaPassH[iEta][iPt] = new TH1F("ZMassSingleEleLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassSingleEleLegPtEtaFailH[iEta][iPt] = new TH1F("ZMassSingleEleLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
      for (int iIso=0; iIso<nIsoBins; ++iIso) {
	ZMassSingleEleLegIsoPtEtaPassH[iIso][iEta][iPt] = new TH1F("ZMassSingleEleLeg_"+IsoBins[iIso]+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
	ZMassSingleEleLegIsoPtEtaFailH[iIso][iEta][iPt] = new TH1F("ZMassSingleEleLeg_"+IsoBins[iIso]+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
	
      }
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
  TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get(PUMCHist);
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
	float puWeight =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	//	float puWeight =  float(PUOfficial_data->GetBinContent(PUOfficial_data->GetXaxis()->FindBin(analysisTree.numtruepileupinteractions)));
	weight *= puWeight;
	//	std::cout << "PU weight = " << puWeight << std::endl;
      }

      //      std::cout << "Triggers" << std::endl;

      // triggers
      isHLTSingleEle = false;
      
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTSingleEle)) {
	  if (it->second==1)
            isHLTSingleEle = true;
	}
      }
      
      //      runNumber = analysisTree.event_run;
      //      eventNumber = analysisTree.event_nr;
      //      lumiBlock = analysisTree.event_luminosityblock;
      //      treePrescales->Fill();

      if (!applyTrigger) {
	isHLTSingleEle = true;
      }

      unsigned int nSingleEleFilter = 0;
      bool isSingleEleFilter = false;
      unsigned int nEleProbeFilter = 0;
      bool isEleProbeFilter = false;
      unsigned int nEleProbeFilter2 = 0;
      bool isEleProbeFilter2 = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==HLTSingleEleFilter) {
	  nSingleEleFilter = i;
	  isSingleEleFilter = true;
	}
	if (HLTFilter==HLTEleProbeFilter) {
	  nEleProbeFilter = i;
	  isEleProbeFilter = true;
	}
	if (HLTFilter==HLTEleProbeFilter2) {
	  nEleProbeFilter2 = i;
	  isEleProbeFilter2 = true;
	}
      }

      if (!isSingleEleFilter)
	std::cout << "Filter " << HLTSingleEleFilter << "  does not exist" << std::endl;
      //      if (!isEleProbeFilter)
	//	std::cout << "Filter " << HLTEleProbeFilter << "  does not exist" << std::endl;
      //      if (!isEleProbeFilter2)
	//	std::cout << "Filter " << HLTEleProbeFilter2 << "  does not exist" << std::endl;

      // elec selection
      vector<int> elecs; elecs.clear();
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {
	if (analysisTree.electron_pt[im]<ptElecCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElecCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElecCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElecCut) continue;
	bool electronId = true;
	if (electronIdType==1)
          electronId = analysisTree.electron_mva_wp80_Iso_Fall17_v1[im]>0.5;
        else if (electronIdType==2)
          electronId = analysisTree.electron_mva_wp90_Iso_Fall17_v1[im]>0.5;
        else if (electronIdType==3)
          electronId = analysisTree.electron_mva_wp80_general_Spring16_v1[im]>0.5;
        else if (electronIdType==4)
          electronId = analysisTree.electron_mva_wp90_general_Spring16_v1[im]>0.5;
	electronId = electronId && analysisTree.electron_pass_conversion[im] > 0.5 &&
	  analysisTree.electron_nmissinginnerhits[im] <= 1;
	if (!electronId) continue;
	elecs.push_back(im);
      }

      if (elecs.size()<2) continue;

      //      std::cout << "nember of electrons : " << elecs.size() << std::endl;

      bool isPairSelected = false;

      // selecting elec pair
      for (unsigned int im1=0; im1<elecs.size(); ++im1) {

	int  ele1Index = elecs[im1];
	bool ele1MatchSingleEle = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.electron_eta[ele1Index],analysisTree.electron_phi[ele1Index],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nSingleEleFilter]) 
	    ele1MatchSingleEle = true;
	}
	bool foundL1ElecSeed = false;
	for (unsigned int il1=0; il1<analysisTree.l1egamma_count; ++il1) {
	  if (analysisTree.l1egamma_bx[il1]!=0) continue;
	  TLorentzVector l1egammaLV; l1egammaLV.SetXYZM(analysisTree.l1egamma_px[il1],
						       analysisTree.l1egamma_py[il1],
						       analysisTree.l1egamma_pz[il1],
						       0.1);
	  if (l1egammaLV.Pt()<ptL1ElecCut) continue;
	  if (fabs(l1egammaLV.Eta())>etaL1ElecCut) continue;
	  float dRmatch = deltaR(analysisTree.electron_eta[ele1Index],analysisTree.electron_phi[ele1Index],
				 l1egammaLV.Eta(),l1egammaLV.Phi());
	  if (dRmatch<DRTrigMatch) 
	    foundL1ElecSeed = true;
	    
	}

	//	if (ele1MatchSingleEle)	std::cout << "ele1 match : " << ele1MatchSingleEle << std::endl;

	bool ele1SingleEle = 
	  ele1MatchSingleEle && 
	  analysisTree.electron_pt[ele1Index]>ptElecTagCut && 
	  fabs(analysisTree.electron_eta[ele1Index])<etaElecTagCut;

	/*
	if (ele1MatchSingleEle&&!foundL1ElecSeed) {
	  std::cout << "Elec   pt = " << analysisTree.electron_pt[ele1Index] 
		    << "   eta = " << analysisTree.electron_eta[ele1Index]
		    << "   L1Seed = " << foundL1ElecSeed
		    << "   HLTmatch = " << ele1MatchSingleEle << std::endl;
	}
	*/

	float q1 = analysisTree.electron_charge[ele1Index];
	
	for (unsigned int im2=0; im2<elecs.size(); ++im2) {

	  if (im1==im2) continue;

	  int  ele2Index = elecs[im2];

	  float q2 = analysisTree.electron_charge[ele2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(analysisTree.electron_eta[ele1Index],analysisTree.electron_phi[ele1Index],
				analysisTree.electron_eta[ele2Index],analysisTree.electron_phi[ele2Index]);

	  bool ele2MatchSingleEleProbe1 = false;
	  bool ele2MatchSingleEleProbe2 = false;

	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.electron_eta[ele2Index],analysisTree.electron_phi[ele2Index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (isEleProbeFilter) {
	      if (analysisTree.trigobject_filters[iT][nEleProbeFilter]) // HLT_SingleEle filter
		ele2MatchSingleEleProbe1 = true;
	    }
	    if (isEleProbeFilter2) { 
	      if (analysisTree.trigobject_filters[iT][nEleProbeFilter2]) // HLT_SingleEle filter
		ele2MatchSingleEleProbe2 = true;
	    }
	  }

	  bool ele2MatchSingleEleProbe = ele2MatchSingleEleProbe1 || ele2MatchSingleEleProbe2;

	  float dZ = fabs(analysisTree.electron_dz[ele1Index]-analysisTree.electron_dz[ele2Index]);

	  TLorentzVector ele1lv; ele1lv.SetXYZM(analysisTree.electron_px[ele1Index],
						analysisTree.electron_py[ele1Index],
						analysisTree.electron_pz[ele1Index],
						electronMass);
	  TLorentzVector ele2lv; ele2lv.SetXYZM(analysisTree.electron_px[ele2Index],
						analysisTree.electron_py[ele2Index],
						analysisTree.electron_pz[ele2Index],
						electronMass);

	  float mass = (ele1lv+ele2lv).M();

	  float absIso1 = 0;
	  if (applyRhoCorrectedIso) {
	    absIso1 = abs_Iso_et(ele1Index, &analysisTree, isoCone);
	  } 
	  else {
	   absIso1 = analysisTree.electron_r03_sumChargedHadronPt[ele1Index];
	   float neutralIso1 = 
	     analysisTree.electron_r03_sumNeutralHadronEt[ele1Index] +
	     analysisTree.electron_r03_sumPhotonEt[ele1Index] -
	     0.5*analysisTree.electron_r03_sumPUPt[ele1Index];
	   neutralIso1 = TMath::Max(float(0),neutralIso1);
	   absIso1 += neutralIso1;
	  }	  
	  float relIso1 = absIso1/analysisTree.electron_pt[ele1Index];

	  float absIso2 = 0;
	  if (applyRhoCorrectedIso) {
            absIso2 = abs_Iso_et(ele2Index, &analysisTree, isoCone);
          }
          else {
	    absIso2 = analysisTree.electron_r03_sumChargedHadronPt[ele2Index];
	    float neutralIso2 =
	      analysisTree.electron_r03_sumNeutralHadronEt[ele2Index] +
	      analysisTree.electron_r03_sumPhotonEt[ele2Index] -
	      0.5*analysisTree.electron_r03_sumPUPt[ele2Index];
	    neutralIso2 = TMath::Max(float(0),neutralIso2);
	    absIso2 += neutralIso2;
          }
	  float relIso2 = absIso2/analysisTree.electron_pt[ele2Index];

	  float ele1RelIso = relIso1;
	  float ele2RelIso = relIso2;

	  if (analysisTree.electron_pt[ele2Index]>analysisTree.electron_pt[ele1Index]) {
	    float temp = relIso1;
	    relIso1 = relIso2;
	    relIso2 = temp;
	  }

	  bool dirIso = (relIso2<isoElecCut) && (relIso1<isoElecCut);

	  bool foundL1TauSeed = false;
	  /*
	  std::cout << "mu1 :  pT = " << analysisTree.electron_pt[ele1Index]
		    << "   eta = " << analysisTree.electron_eta[ele1Index]
		    << "   phi = " << analysisTree.electron_phi[ele1Index] << std::endl;
	  std::cout << "mu2 :  pT = " << analysisTree.electron_pt[ele2Index]
		    << "   eta = " << analysisTree.electron_eta[ele2Index]
		    << "   phi = " << analysisTree.electron_phi[ele2Index] << std::endl;
	  std::cout << "l1 taus : " << analysisTree.l1isotau_count << std::endl;
	  */
	  for (unsigned int iT=0; iT<analysisTree.l1tau_count; ++iT) {
	    /*
	    std::cout << iT << "   pT = " << analysisTree.l1isotau_pt[iT]
		      << "   eta = " << analysisTree.l1isotau_eta[iT]
		      << "   phi = " << analysisTree.l1isotau_phi[iT]
		      << "   iso = " << analysisTree.l1isotau_iso[iT] << std::endl;
	    */
	    TLorentzVector l1LV; l1LV.SetXYZM(analysisTree.l1tau_px[iT],
                                              analysisTree.l1tau_py[iT],
                                              analysisTree.l1tau_pz[iT],
                                              pionMass);

	    if (analysisTree.l1tau_bx[iT]!=0) continue;
	    if (analysisTree.l1tau_iso[iT]==0) continue;
	    if (l1LV.Pt()<ptL1TauCut) continue;
	    if (fabs(l1LV.Eta())>etaL1TauCut) continue;
	    float dRl1tau = deltaR(analysisTree.electron_eta[ele2Index],analysisTree.electron_phi[ele2Index],
				   l1LV.Eta(),l1LV.Phi());
	    if (dRl1tau<DRTrigMatch) continue;
	    float dEtal1tau = fabs(l1LV.Eta()-analysisTree.electron_eta[ele2Index]);
	    if (dEtal1tau<detaL1TauCut) continue;
	    foundL1TauSeed = true;
	  }
	  //	  std::cout << "L1 found = " << foundL1Seed << std::endl;
	  //	  std::cout << std::endl;

	  //	  std::cout << "Ok1" << std::endl;
	  bool foundL1Seed = true;
	  if (matchL1Tau) foundL1Seed = foundL1TauSeed;
	  if (matchL1Elec) foundL1Seed = foundL1ElecSeed;
	  //	  std::cout << "foundL1Seed : " << << std::endl;

	  // **************************************
	  // *********** Tag-&-Probe **************
	  // **************************************
	  if (isHLTSingleEle && ele1SingleEle && foundL1Seed && dRmumu>dRleptonsCut) { 

	    float ele2AbsEta = fabs(analysisTree.electron_eta[ele2Index]);
	    float ele2Pt = TMath::Max(float(5.01),TMath::Min(float(analysisTree.electron_pt[ele2Index]),float(99.9)));
	    float etaFine = analysisTree.electron_eta[ele2Index];
	    if (etaFine<-2.5) etaFine = -2.49;
	    if (etaFine>2.5) etaFine = 2.49;
	    int etaFineBin = binNumber(etaFine,nEtaFineBins,etaFineBins);
	    int etaBin  = binNumber(ele2AbsEta,nEtaBins,etaBins);
	    int ptBin   = binNumber(ele2Pt,nPtBins,ptBins);
	    int isoBin  = binNumber(ele2RelIso,nIsoBins,isoBins);

	    isPairSelected = true;
	    selPairs++;

	    if (ele2AbsEta<etaElecHighCut) { 

	      if (ele2MatchSingleEleProbe) {
		ZMassSingleEleLegIsoPtEtaPassH[isoBin][etaBin][ptBin]->Fill(mass,weight);
	      }
	      else {
		ZMassSingleEleLegIsoPtEtaFailH[isoBin][etaBin][ptBin]->Fill(mass,weight);
	      }

	      if (ele2RelIso<isoElecCut) { // isolated electron path

		if (ele2MatchSingleEleProbe) {
		  ZMassSingleEleLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		  ZMassSingleEleLegPtPassH[ptBin]->Fill(mass,weight);
		}
		else {
		  ZMassSingleEleLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);
		  ZMassSingleEleLegPtFailH[ptBin]->Fill(mass,weight);
		}
		
		if (ele2Pt>ptElecHighCut) {
		  if (ele2MatchSingleEleProbe) {
		    ZMassSingleEleLegEtaPassH[etaBin]->Fill(mass,weight);
		    ZMassSingleEleLegEtaFinePassH[etaFineBin]->Fill(mass,weight);
		    ZMassSingleEleLegPassH->Fill(mass,weight);
		  } 
		  else {
		    ZMassSingleEleLegEtaFailH[etaBin]->Fill(mass,weight);
		    ZMassSingleEleLegEtaFineFailH[etaFineBin]->Fill(mass,weight);
		    ZMassSingleEleLegFailH->Fill(mass,weight);
		  }
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
