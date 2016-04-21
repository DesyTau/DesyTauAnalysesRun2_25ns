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
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist



  // **** configuration
  Config cfg(argv[1]);
  string Channel="eltau";

  // kinematic cuts on electrons
  bool fillplots= false;
  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
  const bool InvertTauIso = cfg.get<bool>("InvertTauIso");
  const bool InvertLeptonIso = cfg.get<bool>("InvertLeptonIso");
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

  const double ptTauLowCut = cfg.get<double>("ptTauLowCut"); 
  const double etaTauCut = cfg.get<double>("etaTauCut"); 

  const string dataBaseDir = cfg.get<string>("DataBaseDir");

  string TrigLeg  ;
  if (!isData) TrigLeg  = cfg.get<string>("El22LegMC");
  if (isData) TrigLeg  = cfg.get<string>("El23LegData");
  const float singleElectronTriggerPtCut = cfg.get<float>("SingleElectronTriggerPtCut");
  const float singleElectronTriggerEtaCut = cfg.get<float>("SingleElectronTriggerEtaCut");

  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");

  const string SingleElectronTriggerFile  = cfg.get<string>("SingleElectronTriggEff");

  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");


  // topological cuts
  const double dRleptonsCuteltau   = cfg.get<double>("dRleptonsCuteltau");
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

  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");
  const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");

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



  const double Lumi   = cfg.get<double>("Lumi");
  const double bTag   = cfg.get<double>("bTag");
  const double metcut = cfg.get<double>("metcut");

  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("No cut after PU");
  CutList.push_back("el");
  CutList.push_back("tau");
  CutList.push_back("Trigger");
  CutList.push_back("DiElectron-Veto");
  CutList.push_back("3rd lepV");
  CutList.push_back("Lepton SF");
  CutList.push_back("TauFakeRate");
  CutList.push_back("topPtRwgt");
  CutList.push_back("Cleaned jets");
  CutList.push_back("b-Veto");
  CutList.push_back("MET gt 50");
  CutList.push_back("MET gt 100");
  CutList.push_back("MET gt 150");
  CutList.push_back("MET gt 200");
  CutList.push_back("InvMass gt 100");
  CutList.push_back("dR gt 3");
  CutList.push_back("Jets lt 3");

  int CutNumb = int(CutList.size());
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
      //datasetName = dt.c_str();
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
      ChiMass=0.0;
    }

  if (XSec<0&& !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}
	
  xsecs=XSec;

  std::vector<unsigned int> allRuns; allRuns.clear();

  cout<<" ChiMass is "<<ChiMass<<"  "<<mIntermediate<<endl;
  bool doThirdLeptVeto=true;
  bool doDiElVeto=true;

  //CutList[CutNumb]=CutListt[CutNumb];
  char ff[100];


  sprintf(ff,"%s/%s",argv[3],argv[2]);


  // reweighting with vertices

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


  // Lepton Scale Factors 

  TH1D * ElSF_IdIso_El1H = new TH1D("ElIdIsoSF_El1H", "ElIdIsoSF_El1", 100, 0.5,1.5);

  ScaleFactor * SF_elIdIso; 
  if (applyLeptonSF) {
    SF_elIdIso = new ScaleFactor();
    SF_elIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
  }

	
  ScaleFactor * SF_electronTrigger = new ScaleFactor();
  SF_electronTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(SingleElectronTriggerFile));

  cout<<" Will try to initialize the TFR now.... "<<endl;
  ScaleFactor * SF_TFR; 
  bool applyTFR = true;
  if (applyTFR) {
    SF_TFR = new ScaleFactor();
    SF_TFR->init_ScaleFactorb(TString(cmsswBase)+"/src/"+TString(TauFakeRateFile),applyTFR);
  }



  double Weight=0;
  int nTotalFiles = 0;
  int iCut=0;
  double CFCounter[CutNumb];
  double statUnc[CutNumb];
  int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; i++){
    CFCounter[i] = 0;
    CFCounter_[i] = 0;
    iCFCounter[i] = 0;
    statUnc[i] =0;
  }
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
  if(InvertLeptonIso) invMuStr = "_InvMuIso_";
  if(InvertTauIso) invTauStr = "_InvTauIso_";
  if(InvertMET) invMETStr = "_InvMET_";

  TString TStrName(rootFileName+invMuStr+invTauStr+invMETStr+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  regionName = invMuStr+invTauStr+invMETStr+"_"+Region.c_str()+"_"+Sign.c_str();
  std::cout <<" The filename will be "<<TStrName <<std::endl;  
  // output fileName with histograms
  TFile * file;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoElons = 0;
  bool lumi=false;
  bool isLowIsoEl=false;
  bool isHighIsoEl = false;
  bool isLowIsoTau=false;
  bool isHighIsoTau = false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  //string treename = rootFileName+"_tree.root";

  SetupTree(); 
  SetupHists(CutNumb); 
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
	///Start reading the ntuple
	//
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
    int NE = (int)histoInputEvents->GetEntries();
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
      //else
	//histWeightsH->Fill(0.,genweight);
    }

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);


    float genweights=1;
    float topPt = -1;
    float antitopPt = -1;
    bool isZTT = false;
    if(!isData) 
      {

	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");

	genweightsTree->SetBranchAddress("genweight",&genweights);
	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) { 
	  genweightsTree->GetEntry(iEntryInit);
	  histWeightsH->Fill(0.,genweights);
	}


      }


    //	    if (std::string::npos != rootFileName.find("TTJetsLO") || std::string::npos != rootFileName.find("TTPow"))


    //numberOfEntries = 1000;

     //numberOfEntries = 1000;
    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1;
      Float_t puweight = 1;
      Float_t topptweight = 1;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      iCut = 0;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
	/////////////////////
      	//Good vertex
	////////////////////
	//
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  

      //isData= false;
      bool lumi=false;
      isLowIsoEl=false;
      isHighIsoEl = false;
      isLowIsoTau=false;
      isHighIsoTau = false;

      Float_t genweights;
      float topPt = 0;
      float antitopPt = 0;
      bool isZTT = false;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;
      if(!isData) 
	{
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

	    if (analysisTree.genparticles_pdgid[igen]==6)
	      topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    if (analysisTree.genparticles_pdgid[igen]==-6)
	      antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);

	  }    

	  if (!isData ) {
	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;
	  }


	  lumi=true;
	  //cout<<"  weight from init "<<genweights<< "  "<<analysisTree.genweight<<"  "<<weight<<endl;
	}


      if (isData)  {
	XSec = 1.;
	histRuns->Fill(analysisTree.event_run);
	///////////////according to dimuons
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;

	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {

	    if ( num.c_str() ==  a.name ) {
	      //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	      //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }

	  }
	//lumi=true;
	if (!lumi) continue;
	//if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;

      }

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




      if (!lumi) continue;
      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();
      mu_index=-1;
      tau_index=-1;
      el_index=-1;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;




      double MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      METV.SetPx(analysisTree.pfmet_ex);	      
      METV.SetPy(analysisTree.pfmet_ey);


      for (unsigned int ijj = 0; ijj<analysisTree.pfjet_count; ++ijj) {
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[ijj], analysisTree.pfjet_py[ijj], analysisTree.pfjet_pz[ijj], analysisTree.pfjet_e[ijj]);
	JetsMV.push_back(JetsV);
      } 


      for (unsigned int imm = 0; imm<analysisTree.muon_count; ++imm) {
	MuV.SetPtEtaPhiM(analysisTree.muon_pt[imm], analysisTree.muon_eta[imm], analysisTree.muon_phi[imm], muonMass);
	MuMV.push_back(MuV);
	//	mu_index=0;
      }

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
	ElMV.push_back(ElV);
	//	el_index=0;
      }

      for (unsigned int itt = 0; itt<analysisTree.tau_count; ++itt) {
	TauV.SetPtEtaPhiM(analysisTree.tau_pt[itt], analysisTree.tau_eta[itt], analysisTree.tau_phi[itt], tauMass);
	TauMV.push_back(TauV);
	//	tau_index=0;
      }



      if (!isData) 
	{
	  if (applyPUreweighting)	 {
	    puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	    weight *=puweight;      
	    pu_weight = puweight;
	  }

	}

      if (fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      bool trigAccept = false;

      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;



      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
//		std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}

      }


      if (!isMainTrigger) {
	std::cout << "HLT filter for SingleEle " << MainTrigger << " not found" << std::endl;
	return(-1);
      }
      /////now clear the Mu.El.Jets again to fill them again after cleaning
      MuMV.clear();
      ElMV.clear();
      TauMV.clear();
      LeptMV.clear();

      double isoElMin = 9999;
      bool el_iso=false;
      vector<int> electrons; electrons.clear();
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {
	if (analysisTree.electron_pt[im]<ptElectronHighCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut) continue;


	bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[im];
	if (!electronMvaId&&applyElectronId) continue;
	if (!analysisTree.electron_pass_conversion[im]&&applyElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[im]>1&&applyElectronId) continue;

	double absIso= analysisTree.electron_r03_sumChargedHadronPt[im]
	  + max(analysisTree.electron_r03_sumNeutralHadronEt[im] + analysisTree.electron_r03_sumPhotonEt[im]
		- 0.5 * analysisTree.electron_r03_sumPUPt[im],0.0);


	double relIso = absIso/analysisTree.electron_pt[im];

	if (relIso<isoElectronLowCut) continue;

	if (double(relIso)<double(isoElMin)) {
	  isoElMin  = relIso;
	  el_index = int(im);
	  el_iso=true;
	  electrons.push_back(im);
	}

	if (relIso == isoElMin && im != el_index) {
	  if (analysisTree.electron_pt[im] > analysisTree.electron_pt[el_index] ) el_index = int(im) ;
	}

      }
      if (electrons.size()==0 || !el_iso ) continue;
	  
      ElV.SetPtEtaPhiM(analysisTree.electron_pt[el_index], analysisTree.electron_eta[el_index], analysisTree.electron_phi[el_index], electronMass);
      ElMV.push_back(ElV);
      LeptMV.push_back(ElV);

      double absIso= analysisTree.electron_r03_sumChargedHadronPt[el_index]
	+ max(analysisTree.electron_r03_sumNeutralHadronEt[el_index] + analysisTree.electron_r03_sumPhotonEt[el_index]
	      - 0.5 * analysisTree.electron_r03_sumPUPt[el_index],0.0);


      double relIso = absIso/analysisTree.electron_pt[el_index];
	
      el_relIso[0]=relIso;
      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      if (LeptMV.size() == 0 ) continue; 
/*
      if (relIso>isoMuonHighCut && !InvertLeptonIso) continue;

      if (relIso>isoMuonHighCutQCD ) { isHighIsoEl=true ;isLowIsoEl=false;}
      else    { isHighIsoEl = false;isLowIsoEl=true;}


      if (InvertLeptonIso && !isHighIsoEl) continue;
      if (!InvertLeptonIso && isHighIsoEl) continue;
      if (InvertLeptonIso && isLowIsoEl) continue;
*/

      //cout<<"  Iso check  "<<relIso<<" InvertLeptonIso "<<InvertLeptonIso<<" isHighIsoEl "<<isHighIsoEl<<" isLowIsoEl "<<isLowIsoEl<<" cutQCD "<<isoMuonHighCutQCD<<endl;


      /////////////////////////////////
      // Apply trigger SF
      // ////////////////////////////
      double ptEl1 = (double)analysisTree.electron_pt[el_index];
      double etaEl1 = (double)analysisTree.electron_eta[el_index];
      float trigweight=1.;


      float El22EffData = (float)SF_electronTrigger->get_EfficiencyData(double(ptEl1),double(etaEl1));
      float El22EffMC   = (float)SF_electronTrigger->get_EfficiencyMC(double(ptEl1),double(etaEl1));

      if (!isData) {
	if (El22EffMC>1e-6)
	  trigweight = El22EffData / El22EffMC;
	weight *= trigweight;
	trig_weight = trigweight;
      }


      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      double isoTauMin = 999;
      bool tau_iso = false;
      vector<int> tau; tau.clear();


      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauLowCut || fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;

	if (analysisTree.tau_againstElectronVLooseMVA5[it]<againstElectronVLooseMVA5) continue;
	if (analysisTree.tau_againstMuonTight3[it]<againstMuonTight3) continue;

	
	////if (!InvertTauIso && analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it] > byCombinedIsolationDeltaBetaCorrRaw3Hits ) continue;
	if (!InvertTauIso && analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[it] < 0.5 ) continue;

	ta_IsoFlag=analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[it];
	double  tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];

	if (tauIso<isoTauMin ) {
	  //      cout<<"  there was a chenge  "<<tauIso<<"  "<<isoTauMin<<" it "<<it<<" tau_index "<<tau_index<<"  "<<analysisTree.tau_count<<endl;
	  isoTauMin  = tauIso;
	  tau_iso=true;
	  tau_index = (int)it;
	  tau.push_back(tau_index);
	  TauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_index], analysisTree.tau_eta[tau_index], analysisTree.tau_phi[tau_index], tauMass);
	  TauMV.push_back(TauV);

	}

	if (tauIso==isoTauMin && it != tau_index) {
	  //analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ? tau_index = it : tau_index = tau_index;
	  if (analysisTree.tau_pt[it] > analysisTree.tau_pt[tau_index] ) tau_index = (int)it ;
	}
      }
      if (tau.size()==0 || !tau_iso ) continue;

 
	ta_relIso[0]=isoTauMin;

////////////////////change to new tau inverted definition 
      if ( abs(analysisTree.tau_charge[tau_index]) != 1) continue;
      if ( abs(analysisTree.electron_charge[el_index]) != 1) continue;
      double q = analysisTree.tau_charge[tau_index] * analysisTree.electron_charge[el_index];

      event_sign  = q;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      /////////////////Make sure that the selected electron fired the trigger and it is within a Dr < 0.5

      bool isdRLeptonMatched = false;
      for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	if (analysisTree.trigobject_filters[iT][nMainTrigger]) { // Mu17 Leg
	  double dRtrig = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
				 analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);

	  if (analysisTree.trigobject_pt[iT]>singleElectronTriggerPtCut && fabs( analysisTree.trigobject_eta[iT])<singleElectronTriggerEtaCut  && dRtrig<deltaRTrigMatch)
	    isdRLeptonMatched = true;

	}
      }

      event_leptonDrTrigger = isdRLeptonMatched;

      if (!isdRLeptonMatched) continue;

      double dR = deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			 analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index]);

      if (dR<dRleptonsCuteltau) continue;


      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



      //Set this flag if there is an opposite-charge muon pair in the event with muons separated by DR>0.15 and both passing the loose selection: 

      bool DiElVeto=false;

      if (doDiElVeto){

	if (electrons.size()>1){
	  for (unsigned  int imv = 0; imv<analysisTree.electron_count; ++imv) {
	    if ( imv != el_index ){
	      double absIso= analysisTree.electron_r03_sumChargedHadronPt[imv]
		+ max(analysisTree.electron_r03_sumNeutralHadronEt[imv] + analysisTree.electron_r03_sumPhotonEt[imv]
		      - 0.5 * analysisTree.electron_r03_sumPUPt[imv],0.0);


	      double relIso = absIso/analysisTree.electron_pt[imv];


	      double dRr = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
				  analysisTree.electron_eta[imv],analysisTree.electron_phi[imv]);
	      bool OSCharge = false;  
	      if ( imv != el_index  && analysisTree.electron_charge[imv] != analysisTree.electron_charge[el_index] ) OSCharge=true;
	      if ( analysisTree.electron_charge[imv] * analysisTree.electron_charge[el_index] < 0) OSCharge=true;

	      if ( analysisTree.electron_charge[imv] != analysisTree.electron_charge[el_index] &&  analysisTree.electron_cutId_veto_Spring15[imv]  
		   &&  analysisTree.electron_pt[imv]> 15 &&  fabs(analysisTree.electron_eta[imv])< 2.5 && fabs(analysisTree.electron_dxy[imv])<0.045 
		   && fabs(analysisTree.electron_dz[imv] < 0.2 && relIso< 0.3 &&   OSCharge)) //removed from last recipe
		DiElVeto=true;
	    }
	  }
	}
      }
	event_secondLeptonVeto = DiElVeto;
      if (DiElVeto) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



      bool ThirdLeptVeto=false;

      if (doThirdLeptVeto){
	if (analysisTree.electron_count>0) {
	  for (unsigned int iev = 0; iev<analysisTree.electron_count; ++iev) {

	    double IsoWithEA = analysisTree.electron_r03_sumChargedHadronPt[iev] 
	      + max(analysisTree.electron_r03_sumNeutralHadronEt[iev] + analysisTree.electron_r03_sumPhotonEt[iev]
		    - 0.5 * analysisTree.electron_r03_sumPUPt[iev], 0.0) ;

	    double relIsoV = IsoWithEA/analysisTree.electron_pt[iev];


	    bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[iev], analysisTree.electron_superclusterEta[iev], analysisTree.electron_mva_id_nontrigPhys14[iev]);


	    if ( iev != el_index && analysisTree.electron_pt[iev] > 10 &&  fabs(analysisTree.electron_eta[iev]) < 2.5 && fabs(analysisTree.electron_dxy[iev])<0.045
		 && fabs(analysisTree.electron_dz[iev]) < 0.2 && relIsoV< 0.3 && electronMvaId && analysisTree.electron_pass_conversion[iev] 
		 && analysisTree.electron_nmissinginnerhits[iev] <=1) ThirdLeptVeto=true;

	  }
	}


	if (analysisTree.muon_count>0){
	  for (unsigned int imvv = 0; imvv<analysisTree.muon_count; ++imvv) {

	    //       if ( imvv != mu_index  && analysisTree.muon_charge[imvv] != analysisTree.muon_charge[mu_index] ){

	    double absIso= analysisTree.muon_r03_sumChargedHadronPt[imvv]
	      + max(analysisTree.muon_r03_sumNeutralHadronEt[imvv] + analysisTree.muon_r03_sumPhotonEt[imvv]
		    - 0.5 * analysisTree.muon_r03_sumPUPt[imvv],0.0);

	    double relIso = absIso/analysisTree.muon_pt[imvv];


	    if ( imvv != mu_index &&  analysisTree.muon_isMedium[imvv] &&  analysisTree.muon_pt[imvv]> 10 &&  fabs(analysisTree.muon_eta[imvv])< 2.4 && fabs(analysisTree.muon_dxy[imvv])<0.045 
		 && fabs(analysisTree.muon_dz[imvv] < 0.2 && relIso< 0.3 && analysisTree.muon_isMedium[imvv]) ) ThirdLeptVeto=true;
	  }
	}
      }

	event_thirdLeptonVeto = ThirdLeptVeto;
      if (ThirdLeptVeto) continue;


      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      if (!isData && applyLeptonSF) {

	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
	double ptEl1 = (double)analysisTree.electron_pt[el_index];
	double etaEl1 = (double)analysisTree.electron_eta[el_index];
	double IdIsoSF_el1 = SF_elIdIso->get_ScaleFactor(ptEl1, etaEl1);

	ElSF_IdIso_El1H->Fill(IdIsoSF_el1);
	//	  if (ptEl1<20||ptMu2<20) {
	//	    std::cout << "mu 1 ->  pt = " << ptEl1 << "   eta = " << etaEl1 << std::endl;
	//	    std::cout << "eff data mu 1 = " << SF_elIdIso->get_EfficiencyData(ptEl1, etaEl1)<< " |  eff mc mu 1 = " << SF_elIdIso->get_EfficiencyMC(ptEl1, etaEl1)<<std::endl;
	//	    std::cout << "mu 2 ->  pt = " << ptMu2 << "   eta = " << etaMu2 << std::endl;
	//	    std::cout << "eff data mu 2 = " << SF_elIdIso->get_EfficiencyData(ptMu2, etaMu2)<< " |  eff mc mu 2 = " << SF_elIdIso->get_EfficiencyMC(ptMu2, etaMu2)<<std::endl;
	//	    std::cout << "SF mu1 = " << IdIsoSF_el1 << std::endl;
	//	    std::cout << "SF mu2 = " << IdIsoSF_mu2 << std::endl;
	//	    
	//	    std::cout << " mass = " << massSel << std::endl;
	//	    std::cout << std::endl;
	//	  }
	weight *= IdIsoSF_el1;
	LSF_weight = IdIsoSF_el1;
      }
      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
 
///////////////Check if the selected tau has a gen matched - if not, apply Tau Fake Rate
      TLorentzVector genTauV;  
	bool isTauMatched = false;
      for (unsigned int gt = 0 ; gt < analysisTree.gentau_count;++gt){

      genTauV.SetXYZT(0.,0.,0.,0.);

      genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);
	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
				genTauV.Eta(), genTauV.Phi());
	if (Drr < 0.3) isTauMatched = true;

      }
      if (!isData && applyTFR && !isTauMatched) {


	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	

	double ptTau1 = (double)analysisTree.tau_pt[tau_index];
	double etaTau1 = (double)analysisTree.tau_eta[tau_index];
	double TFRSF_mu1 = SF_TFR->get_ScaleFactor(ptTau1, etaTau1);

	//MuSF_IdIso_Mu1H->Fill(TFRSF_mu1);
	weight *= TFRSF_mu1;
	TFR_weight  = TFRSF_mu1;
	//cout<<"  "<<TFRSF_mu1<<"  for  eta  "<<etaTau1<<  " pT  "<< ptTau1<<endl;
      }

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

	///////////apply TopPtWeight
      if (!isData && ( string::npos != filen.find("TTJets")  || string::npos != filen.find("TTPowHeg")) ) 
	{

	  if (topPt>0.&&antitopPt>0.) {
	    float topptweight = topPtWeight(topPt,antitopPt);
	    //  cout<<"  "<<topPt<<"  "<<antitopPt<<endl;
	    weight *= topptweight;
	    top_weight = topptweight;
	  }
	}

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      ///////////////////////////////////////////////////////////
      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)analysisTree.muon_count;

      //cout<<" here ============================> "<<iEntry<<"  "<<mu_count<<"  "<<(int)analysisTree.muon_count<<"  "<<analysisTree.muon_count<<endl;
      for (unsigned int im=0;im<analysisTree.muon_count;im++){
	mu_px[im]=analysisTree.muon_px[im];
	mu_py[im]=analysisTree.muon_py[im];
	mu_pz[im]=analysisTree.muon_pz[im];
	mu_eta[im]=analysisTree.muon_eta[im];
	mu_pt[im]=analysisTree.muon_pt[im];
	mu_phi[im]=analysisTree.muon_phi[im];
	mu_charge[im]=analysisTree.muon_charge[im];
	mu_dxy[im]=analysisTree.muon_dxy[im];
	mu_dz[im]=analysisTree.muon_dz[im];
      }
				

      el_count=(int)analysisTree.electron_count;
      for (unsigned int ie=0;ie<analysisTree.electron_count;ie++){
	el_px[ie]=analysisTree.electron_px[ie];
	el_py[ie]=analysisTree.electron_py[ie];
	el_pz[ie]=analysisTree.electron_pz[ie];
	el_eta[ie]=analysisTree.electron_eta[ie];
	el_pt[ie]=analysisTree.electron_pt[ie];
	el_phi[ie]=analysisTree.electron_phi[ie];
	el_charge[ie]=analysisTree.electron_charge[ie];
	el_dxy[ie]=analysisTree.electron_dxy[ie];
	el_dz[ie]=analysisTree.electron_dz[ie];
      }

				
      ta_count=(int)analysisTree.tau_count;
      for (unsigned int it=0;it<analysisTree.tau_count;it++){
	ta_px[it]=analysisTree.tau_px[it];
	ta_py[it]=analysisTree.tau_py[it];
	ta_pz[it]=analysisTree.tau_pz[it];
	ta_eta[it]=analysisTree.tau_eta[it];
	ta_pt[it]=analysisTree.tau_pt[it];
	ta_phi[it]=analysisTree.tau_phi[it];
	ta_charge[it]=analysisTree.tau_charge[it];
	ta_dxy[it]=analysisTree.tau_dxy[it];
	ta_dz[it]=analysisTree.tau_dz[it];
      }
      jet_count=(int)analysisTree.pfjet_count;
      for (unsigned int jj=0;jj<analysisTree.pfjet_count;jj++){

	jet_e[jj] = analysisTree.pfjet_e[jj];
	jet_px[jj] = analysisTree.pfjet_px[jj];
	jet_py[jj] = analysisTree.pfjet_py[jj];
	jet_pz[jj] = analysisTree.pfjet_pz[jj];
	jet_pt[jj] = analysisTree.pfjet_pt[jj];
	jet_eta[jj] = analysisTree.pfjet_eta[jj];
	jet_phi[jj] = analysisTree.pfjet_phi[jj];
	jet_flavour[jj] = analysisTree.pfjet_flavour[jj];
	jet_btag[jj] = analysisTree.pfjet_btag[jj][0];
      }
      met_ex = analysisTree.pfmet_ex;
      met_ey = analysisTree.pfmet_ey;
      met_ez = analysisTree.pfmet_ez;
      met_pt = analysisTree.pfmet_pt;
      met_phi = analysisTree.pfmet_phi;

      all_weight = weight;


      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;

      //      continue;

      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      bool btagged= false;


      JetsMV.clear();
      	/////////////////////////////////////////////
	///////////////// Clean jets from selected eletrons
	/////////////////////////////////////////////
      float jetEtaCut = 2.4;
      float DRmax = 0.5;
      int countjets = 0;
      bool dRmuJet = false;
      bool dRtauJet = false;
      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);

	if (absJetEta > etaJetCut) continue;
	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;



	bool isPFJetId = false ; 
	isPFJetId =looseJetiD(analysisTree,jet);

	//				cout<<" 1- jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  iEntry "<<iEntry<<endl;
	if (!isPFJetId) continue;
	jet_isLoose[jet] = isPFJetId;
	//				cout<<"  jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  "<<iEntry<<endl;
	//for (unsigned int lep=0;LeptMV.size();lep++){
	//double Dr=(LeptMV.at(lep).Eta(),LeptMV.at(lep).Phi(),
	double Dr=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Dr  < DRmax)  continue;

	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Drr  < DRmax)  continue;


	if (analysisTree.pfjet_btag[jet][0]  > bTag) btagged = true;
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
	JetsMV.push_back(JetsV);
	countjets++;
      }

      T->Fill();
	//cout<<"   Will fill the tree now..."<<endl;
      continue;
      ///////////////////////////////////////
      //////////////////////////////////////////

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      if (btagged) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      double ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      if (InvertMET && ETmiss > 50.) continue;
      if (!InvertMET && ETmiss < 50.) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      if (InvertMET && ETmiss > 100.) continue;
      if (!InvertMET && ETmiss < 100.) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;

      iCFCounter[iCut]++;
      iCut++;

      if (InvertMET && ETmiss > 150.) continue;
      if (!InvertMET && ETmiss < 150.) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      if (InvertMET && ETmiss > 200.) continue;
      if (!InvertMET && ETmiss < 200.) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      //if (dRmuJet || dRtauJet || countjets >2) continue;

      // pt Scalar

      //cout<<"  "<<mu_index<<"  "<<tau_index<<"   "<<MuMV.at(mu_index).M()<<"  "<<TauMV.at(tau_index).M()<<endl;

      TLorentzVector muVc ;  muVc.SetPtEtaPhiM(analysisTree.muon_pt[mu_index], analysisTree.muon_eta[mu_index], analysisTree.muon_phi[mu_index], muonMass);
      TLorentzVector tauVc;  tauVc.SetPtEtaPhiM(analysisTree.tau_pt[tau_index], analysisTree.tau_eta[tau_index], analysisTree.tau_phi[tau_index], tauMass);

      TLorentzVector diL = muVc + tauVc;
      if ( diL.M() <100 ) continue;
      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      double MTv = mT(diL,METV);
      /* if (ETmiss < 100) continue;
      //if (MTv>80 && MTv<120) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      if (ETmiss < 120) continue;
      FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      // topological cut
      //if (DZeta<dZetaCut) continue;
      */


      double dRr = deltaR(muVc.Eta(), muVc.Phi(), tauVc.Eta(), tauVc.Phi());

      if (dRr>3) continue;

      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      if (countjets >2) continue;
      if(fillplots)
	FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;



      selEvents++;
    } // end of file processing (loop over events in one file)



    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;


  cout<<" Total events  "<<nEvents<<"  Will use weight  "<<histWeightsH->GetSumOfWeights()<<" Norm Factor for a Lumi of "<<Lumi<<"/pb is "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<endl;
  cout<<" First content "<<CFCounter[0]<<endl;
  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;

  //write out cutflow
  ofstream tfile;
  // TString outname = argv[argc-1];
  TString outname=argv[2];
  TString textfilename = "cutflow_"+outname+"_"+Channel+"_"+argv[3]+".txt";
  //  tfile.open(textfilename);
  cout<<" iCFCounter is the pure event count, CFCounter is normalized to "<<Lumi<<" /pb and CFCounter_ is the sum of weights for a given bin"<<endl;
  //  tfile << "########################################" << endl;
  for(int ci = 0; ci < CutNumb; ci++)
    {
      // tfile << CutList[ci]<<"\t & \t"
      //	    << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
      CFCounter_[iCut] = iCFCounter[iCut];
      CutFlowUnW->SetBinContent(1+ci,0);
      CutFlow->SetBinContent(1+ci,0);
      CutFlowUnW->SetBinContent(1+ci,float(CFCounter[ci]) );
      //CFCounter_->Fill();


      CFCounter[ci] *= float(XSec*Lumi/( histWeightsH->GetSumOfWeights()));

      CutFlow->SetBinContent(1+ci,float(CFCounter[ci]));

      cout << " i "<<ci<<" iCFCounter "<<iCFCounter[ci]<<" CFCounter  "<<CFCounter[ci]<<" CFCounter_ " <<CFCounter_[ci]<<"  "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<" UnWBinContent "<<CutFlowUnW->GetBinContent(1+ci)<<" Norm to Lumi Content "<<CutFlow->GetBinContent(1+ci)<<endl;   
      if (iCFCounter[ci] <0.2) statUnc[ci] =0;
      //else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
      else statUnc[ci] = sqrt(CFCounter[ci]);
    }


  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;



  file->cd(Channel.c_str());
  //WriteTree();
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histRuns->Write();
  CutFlowUnW->Write();
  CutFlow->Write();
  ElSF_IdIso_El1H->Write();
  file->Write();
  file->Close();


  delete file;

}



