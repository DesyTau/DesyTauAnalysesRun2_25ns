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
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "Riostream.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "TGraphAsymmErrors.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;


  // **** configuration
  Config cfg(argv[1]);
  string Channel="mutau";

  // kinematic cuts on electrons
  bool fillplots= true;
  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
  const bool applyPUreweighting_vertices = cfg.get<bool>("ApplyPUreweighting_vertices");
  const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");
 
  const bool InvertTauIso = cfg.get<bool>("InvertTauIso");
  const bool InvertMuIso = cfg.get<bool>("InvertMuIso");
  const bool InvertMET = cfg.get<bool>("InvertMET");
  const double ptElectronLowCut   = cfg.get<double>("ptElectronLowCut");
  const double ptElectronHighCut  = cfg.get<double>("ptElectronHighCut");
  const double etaElectronCut     = cfg.get<double>("etaElectronCut");
  const double dxyElectronCut     = cfg.get<double>("dxyElectronCut");
  const double dzElectronCut      = cfg.get<double>("dzElectronCut");
  const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCut");
  const double isoElectronHighCut = cfg.get<double>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");
  // kinematic cuts on muons
  const double ptMuonLowCuT   = cfg.get<double>("ptMuonLowCut");
  const double ptMuonHighCut  = cfg.get<double>("ptMuonHighCut");
  const double etaMuonCut     = cfg.get<double>("etaMuonCut");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonLowCut  = cfg.get<double>("isoMuonLowCut");
  const double isoMuonHighCut = cfg.get<double>("isoMuonHighCut");
  const double isoMuonHighCutQCD = cfg.get<double>("isoMuonHighCutQCD");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
 
  const double ptMuonLowCutmtau = cfg.get<double>("ptMuonLowCutmtau"); 
  const double etaMuonCutmtau = cfg.get<double>("etaMuonCutmtau"); 


  string TrigLeg  ;
   if (!isData) TrigLeg  = cfg.get<string>("Mu17LegMC");
   if (isData) TrigLeg  = cfg.get<string>("Mu18LegData");
  const string Mu17Tau20MuLegA  = cfg.get<string>("Mu17Tau20MuLegA");
  const string Mu17Tau20MuLegB  = cfg.get<string>("Mu17Tau20MuLegB");
  const string Mu17Tau20TauLegA  = cfg.get<string>("Mu17Tau20TauLegA");
  const string Mu17Tau20TauLegB  = cfg.get<string>("Mu17Tau20TauLegB");

  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");
  
  
  // topological cuts
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double dZetaCut       = cfg.get<double>("dZetaCut");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
 
  // tau
  const double taupt    = cfg.get<double>("taupt");
  const double taueta    = cfg.get<double>("taueta");
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double   decayModeFindingNewDMs  = cfg.get<double>("decayModeFindingNewDMs");
  const double   againstElectronVLooseMVA5  = cfg.get<double>("againstElectronVLooseMVA5");
  const double   againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double   vertexz =  cfg.get<double>("vertexz");
  const double   byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  
  const string dataBaseDir = cfg.get<string>("DataBaseDir");

  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
  
  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  // lepton scale factors
  const string muonSfDataBarrel = cfg.get<string>("MuonSfDataBarrel");
  const string muonSfDataEndcap = cfg.get<string>("MuonSfDataEndcap");
  const string muonSfMcBarrel = cfg.get<string>("MuonSfMcBarrel");
  const string muonSfMcEndcap = cfg.get<string>("MuonSfMcEndcap");
  
  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

  
  // Run-lumi selector
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
	  //std::fstream inputFileStream("temp", std::ios::in);
           periods.push_back(Period());
           std::stringstream ss(s);
           ss >> periods.back();
          }
  }

  TString MainTrigger(TrigLeg);
  TString Muon17Tau20MuLegA (Mu17Tau20MuLegA );
  TString Muon17Tau20MuLegB (Mu17Tau20MuLegB );
  TString Muon17Tau20TauLegA (Mu17Tau20TauLegA );
  TString Muon17Tau20TauLegB (Mu17Tau20TauLegB );



  const double Lumi   = cfg.get<double>("Lumi");
  const double bTag   = cfg.get<double>("bTag");
  const double metcut = cfg.get<double>("metcut");
  
  xs=1;fact=1;fact2=1;
  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;
        
  ifstream ifs("xsecs");
  string line;

  while(std::getline(ifs, line)) // read one line from ifs
    {
		
      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt,st1,st2;st1="stau2_1";st2="stau5_2";
      iss >> dt >> xs >> fact >> fact2;
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      //if (dt==argv[2]) {
      //if (std::string::npos != dt.find(argv[2])) {
      if (  dt == argv[2]) {
	XSec= xs*fact*fact2;
	cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
      }
      if (isData) XSec=1.;
    }

  if (XSec<0&& !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}

  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;

  char ff[100];

	
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  // reading vertex weights
  TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertDataFileName);
  TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertMcFileName);

  TH1D * vertexDataH = (TH1D*)fileDataNVert->Get(TString(vertHistName));
  TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get(TString(vertHistName));

  float normVertexData = vertexDataH->GetSumOfWeights();
  float normVertexMc   = vertexMcH->GetSumOfWeights();

  vertexDataH->Scale(1/normVertexData);
  vertexMcH->Scale(1/normVertexMc);


 
  PileUp * PUofficial = new PileUp();

  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Nov17.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring15_PU25_Startup.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);


  TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataBarrel);  // mu SF barrel data
  TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataEndcap); // mu SF endcap data
  TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcBarrel);  // mu SF barrel MC
  TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcEndcap); // mu SF endcap MC 
  
  TGraphAsymmErrors *hEffBarrelData = (TGraphAsymmErrors*)f10->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapData = (TGraphAsymmErrors*)f11->Get("ZMassEndcap");
  TGraphAsymmErrors *hEffBarrelMC = (TGraphAsymmErrors*)f12->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapMC = (TGraphAsymmErrors*)f13->Get("ZMassEndcap");
  
  double * dataEffBarrel = new double[10];
  double * dataEffEndcap = new double[10];
  double * mcEffBarrel = new double[10];
  double * mcEffEndcap = new double[10];
  
  dataEffBarrel = hEffBarrelData->GetY();
  dataEffEndcap = hEffEndcapData->GetY();
  mcEffBarrel = hEffBarrelMC->GetY();
  mcEffEndcap = hEffEndcapMC->GetY();
  double Weight=0;
  int nTotalFiles = 0;
  int nominator=-1;int denominator = -1;
  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  
  TString era=argv[3];
  TString invMuStr,invTauStr,invMETStr;
  if(InvertMuIso) invMuStr = "_InvMuIso_";
  if(InvertTauIso) invTauStr = "_InvTauIso_";
  if(InvertMET) invMETStr = "_InvMET_";

  TString TStrName(rootFileName+invMuStr+invTauStr+invMETStr+"_"+Region+"_"+Sign);
  std::cout <<" The filename will be "<<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file ;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  TH1::SetDefaultSumw2(true);
  //file->mkdir(Channel.c_str());
  //file->cd();
  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
  TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);
  TH1D * NumberOfVerticesH = new TH1D("NumberOfVerticesH","",51,-0.5,50.5);
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * histWeightsSkimmedH = new TH1D("histWeightsSkimmedH","",1,-0.5,0.5);
  TH1D * hMT = new TH1D("hMT","",20,0,200);
  TH1D * hRatioSum = new TH1D("hRatioSum","",10,0,1);
  TH1D * hDPhi = new TH1D("hDPhi","",70,0,3.5);
  TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);

  int nPtBins = 7;
  float ptBins[8] = {20,25,30,40,60,100,150,200};

  
  int nEtaBins = 3;
  float etaBins[4] = {0,0.9,1.2,2.4}; 
  
  TString PtBins[8] = {"Pt20to25",
		       "Pt25to30",
		       "Pt30to40",
		       "Pt40to60",
		       "Pt60to100",
		       "Pt100to150",
		       "Pt150to200",
		        "PtGt200"};
  

  TString EtaBins[3] = {"EtaLt0p9",
			"Eta0p9to1p2",
		        "EtaGt1p2"};



  TH1D * FakeRatePtIncLoose[3];
  TH1D * FakeRatePtIncTight[3];
  //TH1D * FakeRatePt[3][7];
  //TH1D * FakeRateNV[3][7];
  //TH1D * FakeRateEta[3][7];

  
  TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  
  
  	for (int iEta=0; iEta<nEtaBins; ++iEta) {
      FakeRatePtIncLoose[iEta] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta],"",nPtBins,ptBins);
      FakeRatePtIncTight[iEta] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta],"",nPtBins,ptBins);
    //for (int iPt=0; iPt<nPtBins; ++iPt) {
    //  FakeRatePt[iEta][iPt] = new TH1D("FakeRatePt"+EtaBins[iEta]+PtBins[iPt],"",100,0,1000);
    //  FakeRateNV[iEta][iPt] = new TH1D("FakeRateNV"+EtaBins[iEta]+PtBins[iPt],"",50,0,50);
    //  FakeRateEta[iEta][iPt] = new TH1D("FakeRateEta"+EtaBins[iEta]+PtBins[iPt],"",80,-4,4);
   // }
  
  }

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;
  int iCut=0;
      bool lumi=false;
      bool isLowIsoMu=false;
      bool isHighIsoMu = false;
      bool isLowIsoTau=false;
      bool isHighIsoTau = false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
 
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
    int NE = int(histoInputEvents->GetEntries());
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
    std::cout << "      number of input events         = " << NE << std::endl;

    TTree * _inittree = NULL;
    _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree==NULL) continue;
    Float_t genweight;
    if (!isData)
      _inittree->SetBranchAddress("genweight",&genweight);
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
    std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
      _inittree->GetEntry(iEntry);
      if (isData)
	histWeightsH->Fill(0.,1.);
      else
	histWeightsH->Fill(0.,genweight);
    }

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);



    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 
     
      float weight = 1.;
      float puweight = 1.;
      analysisTree.GetEntry(iEntry);
      nEvents++;
    
      iCut = 0;
      
    //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;
     
      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
     //nominator=0;denominator=0;
 

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  

      //isData= false;
      bool lumi=false;
      isLowIsoMu=false;
      isHighIsoMu = false;
      isLowIsoTau=false;
      isHighIsoTau = false;


      if (!isData )  { 

	  weight *= analysisTree.genweight;
	  lumi=true;
        //cout<<"  weight from init "<<genweights<< "  "<<analysisTree.genweight<<"  "<<weight<<endl;
      
      
	if (applyPUreweighting_vertices) {
	  int binNvert = vertexDataH->FindBin(analysisTree.primvertex_count);
	  float_t dataNvert = vertexDataH->GetBinContent(binNvert);
	  float_t mcNvert = vertexMcH->GetBinContent(binNvert);
	  if (mcNvert < 1e-10){mcNvert=1e-10;}
	  float_t vertWeight = dataNvert/mcNvert;
	  weight *= vertWeight;
	  //	  cout << "NVert = " << analysisTree.primvertex_count << "   weight = " << vertWeight << endl;
	}

        if (applyPUreweighting_official) {
	  nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
	  double Ninteractions = analysisTree.numtruepileupinteractions;
	  double PUweight = PUofficial->get_PUweight(Ninteractions);
	  weight *= float(PUweight);
	  	  PUweightsOfficialH->Fill(PUweight);
	  //	  cout << PUweight << endl;
        }
      
      
      
      
      }


      if (isData)  {
	XSec = 1.;

	//if (analysisTree.event_run<RunRangeMin || analysisTree.event_run>RunRangeMax)
	//      cout<< "  event_run "<<analysisTree.event_run<<"  Min  "<<RunRangeMin<<" Max "<<RunRangeMax<<endl;
	if (analysisTree.event_run<RunRangeMin) continue;
	if (analysisTree.event_run>RunRangeMax) continue;
      

	if (analysisTree.event_run<RunMin)
	  RunMin = analysisTree.event_run;
      
	if (analysisTree.event_run>RunMax)
	  RunMax = analysisTree.event_run;

	//std::cout << " Run : " << analysisTree.event_run << std::endl;

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

   
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;

	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {
        
	    if ( num.c_str() ==  a.name ) {
	      //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	           //std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	       //   std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }
	
	  }
    
      }
      //lumi=true;
      if (!lumi) continue;
      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();
      mu_index=-1;
      tau_index=-1;
      el_index=-1;
      
      //if(fillplots)
      iCut++;
	
      bool trigAccept = false;
      
      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;
      

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}
      

      
      }


      if (!isMainTrigger) {
	std::cout << "HLT filter for Mu20 " << MainTrigger << " not found" << std::endl;
	return(-1);
      }
       ////////// first make a jets collection from pfjet candidates

      	if (analysisTree.pfjet_count<1) continue;

	      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
		float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	
		if (absJetEta > etaJetCut) continue;
                if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
		bool  looseJetID = looseJetiD(analysisTree,jet);
	          if (!looseJetID) continue;	  
	
		JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
		JetsMV.push_back(JetsV);
	      }

		sort(JetsMV.begin(), JetsMV.end(),ComparePt);
		if (JetsMV.size()==0) continue;

		//////make sure that has fired the trigger
		//
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nMainTrigger]) { // Mu17 Leg
	    double dRtrig = deltaR(JetsMV.at(0).Eta(),JetsMV.at(0).Phi(),
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>deltaRTrigMatch) continue;
	 //cout<< " "<<analysisTree.muon_pt[mu_index]<<"   "<<analysisTree.trigobject_pt[iT]<<endl;
	    //if (!isData && analysisTree.trigobject_filters[iT][nMainTrigger] && analysisTree.trigobject_pt[iT]>18) 
	      isMainTrigger = true;
	    
	  }
	}
	
	if (!isMainTrigger) continue;

      iCut++;
	//////find tau candidates
      double isoTauMin = 999;
      bool tau_iso = false;
      bool isTight = false;
      bool isLoose = false;
      unsigned tau_loose=-1;
      vector<int> tau; tau.clear();
      for (unsigned  int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptMuonLowCutmtau || fabs(analysisTree.tau_eta[it])> etaMuonCutmtau) continue;
	if (analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	/*if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
	if (analysisTree.tau_againstElectronVLooseMVA5[it]<againstElectronVLooseMVA5) continue;
	if (analysisTree.tau_againstMuonTight3[it]<againstMuonTight3) continue;
*/
	double  tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];

	if (tauIso < byCombinedIsolationDeltaBetaCorrRaw3Hits && fabs(analysisTree.tau_leadchargedhadrcand_dz[it])< leadchargedhadrcand_dz && analysisTree.tau_againstElectronVLooseMVA5[it]>againstElectronVLooseMVA5 
			&& analysisTree.tau_againstMuonTight3[it]>againstMuonTight3) isTight = true;

	if (tauIso > byCombinedIsolationDeltaBetaCorrRaw3Hits && tauIso < 3*byCombinedIsolationDeltaBetaCorrRaw3Hits) {isLoose = true; tau_loose = int(it);}

	if (isTight){

		
	if (tauIso<isoTauMin ) {
	  isoTauMin  = tauIso;
	  tau_iso=true;
	  tau_index = int(it);
	  tau.push_back(tau_index);
	  TauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_index], analysisTree.tau_eta[tau_index], analysisTree.tau_phi[tau_index], tauMass);
	  TauMV.push_back(TauV);

	}

	if (tauIso==isoTauMin && it != tau_index) {
	  analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ? tau_index = it : tau_index = tau_index;
	}

	}

      }

	//if (!isLoose && !isTight) continue;
	if (!isLoose) continue;


	double dPhi=0;double MT=0 ; double RatioSums=0;

	if (isLoose && int(tau_loose) !=-1 && JetsMV.size()!=0){

      double met = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
	
      double SumPtJetTau = sqrt ( (analysisTree.tau_px[tau_loose]+JetsMV.at(0).Px())*(analysisTree.tau_px[tau_loose]+JetsMV.at(0).Px()) + 
	      		 (analysisTree.tau_py[tau_loose]+JetsMV.at(0).Py())*(analysisTree.tau_py[tau_loose]+JetsMV.at(0).Py()));

      RatioSums = (SumPtJetTau - met )/ (SumPtJetTau + met );
      	

      TLorentzVector MetV; 
      MetV.SetPx(analysisTree.pfmet_ex); 
      MetV.SetPy(analysisTree.pfmet_ey);
      
      TLorentzVector tauV;  tauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_loose], analysisTree.tau_eta[tau_loose], analysisTree.tau_phi[tau_loose], tauMass);

      
       dPhi = dPhiFrom2P( JetsMV.at(0).Px(), JetsMV.at(0).Py(), analysisTree.tau_px[tau_loose], analysisTree.tau_py[tau_loose] );


	hRatioSum->Fill(RatioSums,weight);
	hDPhi->Fill(dPhi, weight);
	}


       if (dPhi < 2 || RatioSums > 0.25) continue;

      //if (tau.size()==0 || !tau_iso ) continue;

      //cout<< " Lets check  "<<mu_index <<"  "<<int(tau_index) <<"  "<<isTight<<"  "<<isLoose<<endl;
      //cout<<"  "<<endl;
	if (isLoose) denominator++;
	if (isTight) nominator++;

	      if (int(tau_loose)!=-1){
	      float ptProbe = TMath::Min(float(analysisTree.tau_pt[tau_loose]),float(ptBins[nPtBins]-0.1));
	      float absEtaProbe = fabs(analysisTree.tau_eta[tau_loose]);
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      if (ptBin<0) continue;
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      if (etaBin<0) continue;

//cout<< "filling here  "<<analysisTree.tau_pt[tau_loose]<<"  "<<ptBin<<"  "<<etaBin<<"  "<<weight<<endl;
		//FakeRatePt[etaBin][ptBin]->Fill(analysisTree.tau_pt[tau_loose],weight);
		
		if (int(tau_loose)!=-1) FakeRatePtIncLoose[etaBin]->Fill(analysisTree.tau_pt[tau_loose],weight);
		if (int(tau_index)!=-1) FakeRatePtIncTight[etaBin]->Fill(analysisTree.tau_pt[tau_index],weight);

		//FakeRateEta[etaBin][ptBin]->Fill(analysisTree.tau_eta[tau_loose],weight);
		//FakeRateNV[etaBin][ptBin]->Fill(analysisTree.tau_vertexz[tau_loose],weight);

	      }
	      

      	      // double q = analysisTree.tau_charge[tau_index] * analysisTree.muon_charge[mu_index];
  
        // if (q>0 && Sign=="OS" ) continue;
        // if (q<0 && Sign=="SS" ) continue;


      iCut++;


      selEvents++;
    } // end of file processing (loop over events in one file)
  
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;
	double fr =  double(nominator/denominator) ;
  if (denominator>0) cout <<" the tau fake rate "<<fr<<endl;
    cout <<" the tau fake rate "<<nominator<<"  "<<denominator<<endl;
  cout<<" Total events  "<<nEvents<<"  Will use weight  "<<histWeightsH->GetSumOfWeights()<<" Norm Factor "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<endl;
  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;


  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  file->cd();
  hxsec->Fill(XSec);
 // hxsec->Write();
  file->Write();
  file->Close();
   
  delete file;
  
}



