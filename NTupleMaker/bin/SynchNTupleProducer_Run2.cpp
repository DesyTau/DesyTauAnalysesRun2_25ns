#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TSystem.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys_WIP.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"
#include "HTT-utilities/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"

#define pi 	3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 		   0.105658
#define tauMass 		   1.77682
#define pionMass 		   0.1396


void FillMuTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso);
void FillETau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso);
void FillTau(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex);
void initializeCPvar(Synch17Tree *otree);
void SaveRECOVertices(const AC1B * analysisTree,Synch17Tree *otree, const bool isData);
void initializeGenTree(Synch17GenTree *gentree);
void FillGenTree(const AC1B * analysisTree, Synch17GenTree *gentree);

/*Notifications Merijn
-to my experience currently macro only works for mu-tau (unless udpated in the meantime)
-added a histogram and leafs in the otree to keep track of the pdg codes of the hadronically decaying tau, and ctr to keep track how often 2nd tau does NOT decay hadronically
-we call an improved CP function, in which a number of bugs are solved. It takes the decay mode also as argument, to construct certain kinematic variables
-Merijn suspects that second reco tau in mu-tau is currently not necessarily hadronically decaying, since RECO decay products substatial fraction 
-see functionsCP.h for updates on the CP calculation
*/

int main(int argc, char * argv[]){
// first argument - config file for analysis
// second argument - file list (MUST BE IN THE SAME DIRECTORY OF THE EXECUTABLE)
// third argument - channel ("et" or "mt")
// third argument - index of first file to run on (optional, ignored if only one file is used)
// fourth argument - index of last file to run on (optional, ignored if only one file is used)

  using namespace std;
  gErrorIgnoreLevel = kFatal;
  string cmsswBase = (getenv("CMSSW_BASE"));
  
  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load(pathToCrystalLib);
  if (openSuccessful != 0) {
    cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit(-1);
  }

  if(argc < 4){
    std::cout << "RUN ERROR: wrong number of arguments"<< std::endl;
    std::cout << "Please run the code in the following way:"<< std::endl;
    std::cout << "SynchNtupleProducer_2017 NameOfTheConfigurationFile FileList Channel" << std::endl;
    std::cout << "example: SynchNtupleProducer_2017 analysisMacroSynch_lept_mt_DATA_SingleMuon.conf DATA_SingleMuon mt" << std::endl;
    exit(-1);
  }

  // **** configuration analysis  
  Config cfg(argv[1]);

  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];
  TString ch = argv[3];
  std::string lep;

  if ((ch != "mt" ) && (ch != "et")) {
  	std::cout << " Channel " << ch << " is not a valid choice. Please try again with 'mt' or 'et'. Exiting. " << std::endl;
    exit(0);
  }
  
  if (ch == "mt") lep = "Muon"; 
  if (ch == "et") lep = "Electron";

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase) + "/src/" + TString(json_name)).Data(), json);
  }

  const int era = cfg.get<int>("era");
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters  = cfg.get<bool>("ApplyMetFilters");

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");

  //lep eff
  const string idIsoEffFile = cfg.get<string>("idIsoEffFile");
  const string singleLepTrigEffFile = cfg.get<string>("singleLepTrigEffFile");
  const string xTrigLepLegEffFile   = cfg.get<string>("xTrigLepLegEffFile");

  const string idIsoEffFile_antiiso = cfg.get<string>("idIsoEffFile_antiiso");
  const string singleLepTrigEffFile_antiiso = cfg.get<string>("singleLepTrigEffFile_antiiso");
  const string xTrigLepLegEffFile_antiiso   = cfg.get<string>("xTrigLepLegEffFile_antiiso");

  // tau trigger efficiency
  TauTriggerSFs2017 *tauTriggerSF = new TauTriggerSFs2017(cmsswBase + "/src/HTT-utilities/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root", "tight");
  	
  //svfit
  const string svFitPtResFile = TString(TString(cmsswBase) + "/src/" + TString(cfg.get<string>("svFitPtResFile"))).Data();

  //zptweight file 
  const string ZptweightFile = cfg.get<string>("ZptweightFile");

  //b-tag scale factors
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  TString pathToBtagScaleFactors = (TString) cmsswBase + "/src/" + BtagSfFile;
  if( ApplyBTagScaling && gSystem->AccessPathName(pathToBtagScaleFactors) ){
    cout<<pathToBtagScaleFactors<<" not found. Please check."<<endl;
    exit(-1);
  }//cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv"
  
  BTagCalibration calib("csvv2", (string) pathToBtagScaleFactors );
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM, "central");
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM, "central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM, "central");
  if(ApplyBTagScaling){
    reader_B.load(calib, BTagEntry::FLAV_B, "comb");
    reader_C.load(calib, BTagEntry::FLAV_C, "comb");
    reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  }
  
  const string TaggingEfficienciesFile = cfg.get<string>("BtagMCeffFile");
  TString pathToTaggingEfficiencies = (TString) cmsswBase + "/src/" + TaggingEfficienciesFile;
  if (ApplyBTagScaling && gSystem->AccessPathName(pathToTaggingEfficiencies)){
    cout<<pathToTaggingEfficiencies<<" not found. Please check."<<endl;
    exit(-1);
  }
  
  TFile *fileTagging  = new TFile(pathToTaggingEfficiencies);
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  
  if(ApplyBTagScaling){
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  }
  TRandom3 *rand = new TRandom3();
  const struct btag_scaling_inputs inputs_btag_scaling_medium = {reader_B, reader_C, reader_Light, tagEff_B, tagEff_C, tagEff_Light, rand};

  // MET Recoil Corrections
  const bool ApplyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  const bool isDY = infiles.find("DY") == infiles.rfind("/")+1;
  const bool isWJets = (infiles.find("WJets") == infiles.rfind("/")+1) || (infiles.find("W1Jets") == infiles.rfind("/")+1) || (infiles.find("W2Jets") == infiles.rfind("/")+1) || (infiles.find("W3Jets") == infiles.rfind("/")+1) || (infiles.find("W4Jets") == infiles.rfind("/")+1) || (infiles.find("EWK") == infiles.rfind("/")+1);
  const bool isVBForGGHiggs = (infiles.find("VBFHTo")== infiles.rfind("/")+1) || (infiles.find("GluGluHTo")== infiles.rfind("/")+1);
  const bool isEWKZ =  infiles.find("EWKZ") == infiles.rfind("/")+1;
  const bool isMG = infiles.find("madgraph") != string::npos;
  const bool isMSSMsignal =  (infiles.find("SUSYGluGluToHToTauTau")== infiles.rfind("/")+1) || (infiles.find("SUSYGluGluToBBHToTauTau")== infiles.rfind("/")+1);
  
  RecoilCorrector *recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
  
  if(!isData && ApplyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal) ){
    TString RecoilFilePath = cfg.get<string>("RecoilFilePath");
    std::cout << RecoilFilePath << std::endl;
    recoilPFMetCorrector = new RecoilCorrector( RecoilFilePath);
  }
  
  // Read in HLT filter
  vector<string> filterSingleLep;
  vector<string> filterXtriggerLepLeg;
  vector<string> filterXtriggerTauLeg;
  
  filterSingleLep.push_back(cfg.get<string>("filterSingleLep1"));
  filterSingleLep.push_back(cfg.get<string>("filterSingleLep2"));
  filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg1"));
  filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg2"));
  filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg1"));
  filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg2"));
  
  cout<<"Number of single lepton trigger legs = "<<filterSingleLep.size()<<endl;
  cout<<"Number of X trigger legs (lep leg)   = "<<filterXtriggerLepLeg.size()<<endl;
  cout<<"Number of X trigger legs (tau leg)   = "<<filterXtriggerTauLeg.size()<<endl;

  // tau cuts
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut       = cfg.get<float>("dzTauCut");

  // tau energy scale corrections
  const float shift_tes_1prong = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  
  // for lep->tau fakes
  const float shift_tes_lepfake_1prong = cfg.get<float>("TauEnergyScaleShift_LepFake_OneProng");
  const float shift_tes_lepfake_1p1p0  = cfg.get<float>("TauEnergyScaleShift_LepFake_OneProngOnePi0");
  const float shift_tes_lepfake_3prong = cfg.get<float>("TauEnergyScaleShift_LepFake_ThreeProng");
  
  // for tau Id efficiency
  const float tau_id_sf = cfg.get<float>("TauIdSF");
  
  // pair selection
  const float dRleptonsCut = cfg.get<float>("dRleptonsCut");

  // extra electron veto
  const float ptVetoElectronCut  = cfg.get<float>("ptVetoElectronCut");  
  const float etaVetoElectronCut = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut = cfg.get<float>("dxyVetoElectronCut");  
  const float dzVetoElectronCut  = cfg.get<float>("dzVetoElectronCut"); 
  const bool applyVetoElectronId = cfg.get<bool>("applyVetoElectronId");
  const float isoVetoElectronCut = cfg.get<float>("isoVetoElectronCut");  
  
  // extra muon veto
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");  
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut  = cfg.get<float>("dxyVetoMuonCut");  
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut"); 
  const bool  applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
  const float isoVetoMuonCut  = cfg.get<float>("isoVetoMuonCut");

 // lepton cuts
  const float ptLeptonLowCut   = cfg.get<float>("pt" + lep + "LowCut");
  const float ptLeptonHighCut  = cfg.get<float>("pt" + lep + "HighCut");
  const float etaLeptonCut     = cfg.get<float>("eta" + lep + "Cut");
  const float dxyLeptonCut     = cfg.get<float>("dxy" + lep + "Cut");
  const float dzLeptonCut      = cfg.get<float>("dz" + lep + "Cut");

  cout<<"dxyLeptonCut "<<dxyLeptonCut<<endl;
  cout<<"dzLeptonCut "<<dzLeptonCut<<endl;
  cout<<"dzTauCut "<<dzTauCut<<endl;

  const bool  ApplyLeptonId    = cfg.get<bool>("Apply" + lep + "Id");

  const float dRTrigMatch = cfg.get<float>("dRTrigMatch");
  const float dRiso = cfg.get<float>("dRiso");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");

  // check overlap
  const bool checkOverlap = cfg.get<bool>("CheckOverlap");
  const bool debug = cfg.get<bool>("debug");
  
  // **** end of configuration analysis

  //file list creation
  int ifile = 0;
  int jfile = -1;

  if (argc > 4)
    ifile = atoi(argv[4]);
  if (argc > 5)
    jfile = atoi(argv[5]);
  
  // create input files list
  std::vector<std::string> fileList;  
  int NumberOfFiles = 0;
  if (infiles.find(".root") != std::string::npos){
    ifile = 0;
    jfile = 1;
    fileList.push_back(infiles);
  }
  else {
    ifstream input;
    std::string infile;  
    input.open(infiles);

    while(true){
      input>>infile;
      if(!input.eof()){
	       if (infile.length() > 0){
	          fileList.push_back(infile);
	          NumberOfFiles += 1 ;
	         }
         }
      else
	     break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

  if(NumberOfFiles < jfile) jfile = NumberOfFiles;

  for (int iF = ifile; iF < jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }

  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");

  // PU reweighting - initialization
  PileUp *PUofficial = new PileUp();
  if(ApplyPUweight){
    TFile *filePUdistribution_data = new TFile(TString(cmsswBase) + "/src/" + TString(pileUpInDataFile), "read");
    TFile *filePUdistribution_MC = new TFile (TString(cmsswBase) + "/src/" + TString(pileUpInMCFile), "read");
    TH1D *PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    //std::cout << filePUdistribution_data << std::endl;
    //std::cout << filePUdistribution_MC << std::endl;
    
    TH1D *PU_mc;
    if (era == 2017) {
      PU_mc = (TH1D *)filePUdistribution_MC->Get(TString(pileUpforMC));
    } else {
      PU_mc = (TH1D *)filePUdistribution_MC->Get(TString("pileup"));
    }
    if (PU_mc == NULL) {
      std::cout << "Histogram " << pileUpforMC << " is not present in pileup file" << std::endl;
      exit(-1);
    }
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  

  // Lepton Scale Factors
  ScaleFactor *SF_lepIdIso = new ScaleFactor();
  ScaleFactor *SF_lepIdIso_antiiso = new ScaleFactor();
  ScaleFactor *SF_SingleLepTrigger = new ScaleFactor();
  ScaleFactor *SF_SingleLepTrigger_antiiso = new ScaleFactor();
  ScaleFactor *SF_XTriggerLepLeg   = new ScaleFactor();
  ScaleFactor *SF_XTriggerLepLeg_antiiso   = new ScaleFactor();

  if(ApplyLepSF){
    SF_lepIdIso->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(idIsoEffFile));
    SF_lepIdIso_antiiso->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(idIsoEffFile_antiiso));
    SF_SingleLepTrigger->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(singleLepTrigEffFile));
    SF_SingleLepTrigger_antiiso->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(singleLepTrigEffFile_antiiso));
    SF_XTriggerLepLeg->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(xTrigLepLegEffFile));
    SF_XTriggerLepLeg_antiiso->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(xTrigLepLegEffFile_antiiso));
  }
  // For tau leg of cross-trigger a different implementation is used

  // Workspace containing tracking efficiency weights 
  //TString workspace_filename = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_3.root";
  TString workspace_filename = TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/htt_scalefactors_v17_1.root";
  TFile *f_workspace = new TFile(workspace_filename, "read");
  if (f_workspace->IsZombie()) {
    std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
     exit(1);
   }
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");

  // Zpt reweighting for LO DY samples 
  TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
  //TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+"DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016_BtoH.root","read");
  //Merijn: the file will point now for 2017  instead to DesyTauAnalyses/NTupleMaker/data/zpt_weights_2017.root
  TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

  // lepton to tau fake init
  LepTauFakeRate *leptauFR = new LepTauFakeRate();
  leptauFR->Init();

  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_" + ch + "_Sync.root";
    
  std::cout <<rootFileName <<std::endl;  

  TFile *file = new TFile(rootFileName, "recreate");
  file->cd("");

  TH1D *inputEventsH = new TH1D("inputEventsH", "", 1, -0.5, 0.5);
  TH1D *nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5, 0.5);
  
  TTree *tree = new TTree("TauCheck", "TauCheck");
  TTree *gtree = new TTree("GenTauCheck", "GenTauCheck");

  //Merijn added a histogram to spot the pdg codes of the decaying hadronic tau
  TH1F *ConstitsPDG = new TH1F("ConstitsPDG", "ConstitsPDG", 500, -250, 250);
  int nonpionphotonctr = 0;

  Synch17Tree *otree = new Synch17Tree(tree);
  initializeCPvar(otree);
  
  //Synch17GenTree *gentree = new Synch17GenTree(gtree);
  
  Synch17GenTree *gentree = new Synch17GenTree(gtree);
  Synch17GenTree *gentreeForGoodRecoEvtsOnly = new Synch17GenTree(tree);
    
  int nTotalFiles = 0;
  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  vector<unsigned int> runList; runList.clear();
  vector<unsigned int> eventList; eventList.clear();

  int nonOverlap = 0;

  std::ifstream fileEvents("overlap.txt");
  int Run, Event, Lumi;
  if (checkOverlap) {
    std::cout << "Non-overlapping events ->" << std::endl;
    while (fileEvents >> Run >> Event >> Lumi) {
      runList.push_back(Run);
      eventList.push_back(Event);
      std::cout << Run << ":" << Event << std::endl;
    }
    std::cout << std::endl;
  }
  std::ofstream fileOutput("overlap.out");


  //svFit
  TH1::AddDirectory(false);  
  TFile *inputFile_visPtResolution = new TFile(svFitPtResFile.data());

  //Systematics init
  TauScaleSys *tauScaleSys = 0;
  TauOneProngScaleSys *tauOneProngScaleSys = 0;
  TauOneProngOnePi0ScaleSys *tauOneProngOnePi0ScaleSys = 0;
  TauThreeProngScaleSys *tauThreeProngScaleSys = 0;

  LepTauFakeOneProngScaleSys *lepTauFakeOneProngScaleSys = 0;
  LepTauFakeOneProngOnePi0ScaleSys *lepTauFakeOneProngOnePi0ScaleSys = 0;
  LepTauFakeThreeProngScaleSys  *lepTauFakeThreeProngScaleSys = 0;

  ZPtWeightSys* zPtWeightSys = 0;
  TopPtWeightSys* topPtWeightSys = 0;
  std::vector<JetEnergyScaleSys*> jetEnergyScaleSys;
  JESUncertainties * jecUncertainties = 0;

  if(!isData && ApplySystShift){
    tauScaleSys = new TauScaleSys(otree);
    tauScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauScaleSys->SetUseSVFit(ApplySVFit);
    zPtWeightSys = new ZPtWeightSys(otree);
    topPtWeightSys = new TopPtWeightSys(otree);

    if (cfg.get<bool>("splitJES")){
      JESUncertainties *jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Summer16_UncertaintySources_AK4PFchs.txt");
      std::vector<std::string> JESnames = jecUncertainties->getUncertNames();
      for (unsigned int i = 0; i < JESnames.size(); i++) std::cout << "i: "<< i << ", JESnames.at(i) : " << JESnames.at(i) << std::endl;
      for (unsigned int i = 0; i < JESnames.size(); i++){
        JetEnergyScaleSys *aJESobject = new JetEnergyScaleSys(otree, TString(JESnames.at(i)));
        aJESobject->SetConfig(&cfg);
        aJESobject->SetBtagScaling(&inputs_btag_scaling_medium);
        aJESobject->SetJESUncertainties(jecUncertainties);
        jetEnergyScaleSys.push_back(aJESobject);
      }	  
    }
    else { // use JEC uncertainty from analysis tree
      JetEnergyScaleSys *singleJES = new JetEnergyScaleSys(otree, TString("JES"));
      singleJES->SetConfig(&cfg);
      singleJES->SetBtagScaling(&inputs_btag_scaling_medium);
      singleJES->SetJESUncertainties(jecUncertainties);
      jetEnergyScaleSys.push_back(singleJES);
    }    
  }

  // list of met filters from config
  std::vector<TString> met_filters_list;
  for (unsigned int i = 1; i < (unsigned int) cfg.get<int>("num_met_filters") + 1; i++) {
    met_filters_list.push_back(cfg.get<string>("met_filter_" + std::to_string(i)));
  }
  int counter[20] = {0};

  ///////////////FILE LOOP///////////////

  for (int iF = ifile; iF < jfile; ++iF) {
    std::cout << "file " << iF + 1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    
    TFile *file_ = TFile::Open(fileList[iF].data());
    TTree *_tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));  
    if (_tree == NULL) continue;
    
    TH1D *histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents == NULL) continue;
    int NE = int(histoInputEvents->GetEntries());
    std::cout << "      number of input events    = " << NE << std::endl;
    for (int iE = 0; iE < NE; ++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    // set AC1B for JES systematicsf_
    if (!isData && ApplySystShift && jetEnergyScaleSys.size() > 0){
      for (unsigned int i = 0; i < jetEnergyScaleSys.size(); i++)
      	(jetEnergyScaleSys.at(i))->SetAC1B(&analysisTree);
    }

    Long64_t numberOfEntries = analysisTree.GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    ///////////////EVENT LOOP///////////////

    for (Long64_t iEntry = 0; iEntry < numberOfEntries; iEntry++) {
      counter[0]++;
      analysisTree.GetEntry(iEntry);
      nEvents++;
    
      if (isData)
        nWeightedEventsH->Fill(0., 1.);
      else {
      	nWeightedEventsH->Fill(0., analysisTree.genweight);
      	FillGenTree(&analysisTree,gentree);
      	gentree->Fill();
      }
      
      //Skip events not passing the MET filters, if applied
      if (ApplyMetFilters && !passedAllMetFilters(&analysisTree, met_filters_list)) continue;
      counter[1]++;
      
      // Check if all triggers are existent in each event and save index
      vector<int> nSingleLepTrig(filterSingleLep.size(), -1);
      vector<int> nXTrigLepLeg(filterXtriggerLepLeg.size(), -1);
      vector<int> nXTrigTauLeg(filterXtriggerTauLeg.size(), -1);
      
      if(ApplyTrigger){
      	vector<bool> checkFilterSingleLep(filterSingleLep.size(), false); 
      	vector<bool> checkFilterXTrigLepLeg(filterXtriggerLepLeg.size(), false); 
      	vector<bool> checkFilterXTrigTauLeg(filterXtriggerTauLeg.size(), false);
      	unsigned int nfilters = analysisTree.run_hltfilters->size();
      
      	for (unsigned int i = 0; i < nfilters; ++i) {
      	  TString HLTFilter(analysisTree.run_hltfilters->at(i));
      	  for(unsigned int i_trig = 0; i_trig < filterSingleLep.size(); i_trig++){
      	    if (HLTFilter == filterSingleLep.at(i_trig)){ nSingleLepTrig.at(i_trig) = i; checkFilterSingleLep.at(i_trig) = true;}
      	  }
      	  for(unsigned int i_trig = 0; i_trig < filterXtriggerLepLeg.size(); i_trig++){
      	    if (HLTFilter == filterXtriggerLepLeg.at(i_trig)){ nXTrigLepLeg.at(i_trig) = i; checkFilterXTrigLepLeg.at(i_trig) = true;}
      	  }
      	  for(unsigned int i_trig = 0; i_trig < filterXtriggerTauLeg.size(); i_trig++){
      	    if (HLTFilter == filterXtriggerTauLeg.at(i_trig)){ nXTrigTauLeg.at(i_trig) = i; checkFilterXTrigTauLeg.at(i_trig) = true;}
      	  }
      	}
      }
    
      if (nEvents % 10000 == 0) 
      	cout << "      processed " << nEvents << " events" << endl; 
    
      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
    
      bool overlapEvent = true;
      for (unsigned int iEvent = 0; iEvent < runList.size(); ++iEvent) {
      	if (runList.at(iEvent) == otree->run && eventList.at(iEvent) == otree->evt) {
      	  overlapEvent = false;	  
      	}
      }
    
      if (overlapEvent && checkOverlap) continue;
      nonOverlap++;
      counter[2]++;
    
      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;
    
      initializeGenTree(gentree);
    
       // weights
      if(ApplyPUweight) fill_weight(&analysisTree, otree, PUofficial, isData);
    
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;
    
      // tau selection
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it < analysisTree.tau_count; ++it) { 
        if (analysisTree.tau_pt[it] <= ptTauLowCut) continue;
        if (fabs(analysisTree.tau_eta[it]) >= etaTauCut) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it]) >= dzTauCut) continue;
        if (fabs(fabs(analysisTree.tau_charge[it]) - 1) > 0.001) continue;
    
      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] < 0.5) continue;
      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[it] < 0.5) continue;
      	if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] < 0.5) continue;
    
        if (analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) continue; //always true, cut applied in NTupleMaker
        if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
    
        taus.push_back(it);
      }
      counter[3]++;
    
      //lepton selection
      vector<int> leptons; leptons.clear();
      if(ch == "et"){
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
          bool electronMvaId = analysisTree.electron_mva_wp80_Iso_Fall17_v1[ie];
      	  // bool electronMvaId = analysisTree.electron_mva_wp80_general_Spring16_v1[ie];
    
          if (analysisTree.electron_pt[ie] <= ptLeptonLowCut) continue;
          if (fabs(analysisTree.electron_eta[ie]) >= etaLeptonCut) continue;
          if (fabs(analysisTree.electron_dxy[ie]) >= dxyLeptonCut) continue;
      	  if (fabs(analysisTree.electron_dz[ie]) >= dzLeptonCut) continue;
          if (!electronMvaId && ApplyLeptonId) continue;
    
      	  //Meirjn 2019 8 20: reinstated. They are mentioned in the legacy twiki
      	  if (!analysisTree.electron_pass_conversion[ie] && ApplyLeptonId) continue;
      	  if (analysisTree.electron_nmissinginnerhits[ie] > 1 && ApplyLeptonId) continue;
          leptons.push_back(ie);
        }
      }
    
      if(ch == "mt"){
        for (unsigned int im = 0; im < analysisTree.muon_count; ++im) {
          bool muonMediumId = isIdentifiedMediumMuon(im, &analysisTree, isData);	          
          if (analysisTree.muon_pt[im] <= ptLeptonLowCut) continue;
          if (fabs(analysisTree.muon_eta[im]) >= etaLeptonCut) continue;
          if (fabs(analysisTree.muon_dxy[im]) >= dxyLeptonCut) continue;
          if (fabs(analysisTree.muon_dz[im]) >= dzLeptonCut) continue;
          if (!muonMediumId && ApplyLeptonId) continue;
          leptons.push_back(im);
          }
      }
      counter[4]++;
    
    
    	if (leptons.size() == 0) continue;
      if (taus.size() == 0) continue;
      counter[5]++;
    
      // selecting electron and tau pair (OS or SS) or ditau pair;
      int leptonIndex = -1;
      int tauIndex = -1;
    
      float isoLepMin   = 1e+10;
      float isoTauMax   = -1;      
      float lep_pt_max  = -1;
      float tau_pt_max  = -1;
    
      //////////////LOOP on Taus/////////////
    
      for (unsigned int it = 0; it < taus.size(); ++it) {
        counter[6]++;
      	unsigned int tIndex = taus.at(it);
    
    	//////////////LOOP on Leptons or second Tau/////////////
    
        for (unsigned int il = 0; il < leptons.size(); ++il) {
          unsigned int lIndex = leptons.at(il);
    
          float lep_pt     = -9999.;
          float lep_pt_max = -9999.;
          float lep_eta    = -9999.;
          float lep_phi    = -9999.;
          float relIsoLep  = -9999.;
    
          if (ch == "mt")  relIsoLep = (abs_Iso_mt(lIndex, &analysisTree, dRiso) / analysisTree.muon_pt[lIndex] );
          if (ch == "et")  relIsoLep = (abs_Iso_et(lIndex, &analysisTree, dRiso) / analysisTree.electron_pt[lIndex] );
    
          if(ch == "mt"){
            lep_pt  = analysisTree.muon_pt[lIndex]; 
            lep_eta = analysisTree.muon_eta[lIndex]; 
            lep_phi = analysisTree.muon_phi[lIndex];
          }
          if(ch == "et"){         
            lep_pt  = analysisTree.electron_pt[lIndex];
            lep_eta = analysisTree.electron_eta[lIndex]; 
            lep_phi = analysisTree.electron_phi[lIndex];
          }
          counter[7]++;
    
          float sortIsoTau = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex];
          float dR = deltaR(analysisTree.tau_eta[tIndex], analysisTree.tau_phi[tIndex], lep_eta, lep_phi);
          if (dR < dRleptonsCut) continue;
    
          // change pair
          bool changePair =  false;
          if (relIsoLep < isoLepMin)
           changePair = true;
          else if (fabs(relIsoLep - isoLepMin) < 1.e-5)
           {
                 if (lep_pt > lep_pt_max)
                  changePair = true;
                 else if (fabs(lep_pt - lep_pt_max) < 1.e-5) 
                 {
        	              if (sortIsoTau > isoTauMax)
                          changePair = true;
        	              else if ((sortIsoTau - isoTauMax) < 1.e-5)
                        {
        	                    if (analysisTree.tau_pt[tIndex] > tau_pt_max)
                                changePair = true;
        	              }
                 }
            }
          counter[8]++;
    
          if (changePair) {
            isoLepMin = relIsoLep;
            lep_pt_max = lep_pt;
            tau_pt_max = analysisTree.tau_pt[tIndex];
            leptonIndex = lIndex;
            isoTauMax = sortIsoTau;
            tauIndex = tIndex;
          }
        } // lepton loop
      } // tau loop
    
      if (leptonIndex < 0) continue;
      if (tauIndex < 0) continue;
      counter[9]++;
    
      float lep_pt  = -9999.;
      float lep_eta = -9999.;
      float lep_phi = -9999.;
    
      if (ch == "mt") {
      	lep_pt =      analysisTree.muon_pt[leptonIndex]; 
      	lep_eta =     analysisTree.muon_eta[leptonIndex]; 
      	lep_phi =     analysisTree.muon_phi[leptonIndex];
      }
      if(ch == "et") {         
      	lep_pt  = analysisTree.electron_pt[leptonIndex];
      	lep_eta = analysisTree.electron_eta[leptonIndex]; 
      	lep_phi = analysisTree.electron_phi[leptonIndex];
      }
    
      bool isSingleLepTrig = false;
      vector<bool> isXTrigLepLeg(filterXtriggerLepLeg.size(), false);
      vector<bool> isXTrigTauLeg(filterXtriggerTauLeg.size(), false);
      bool isXTrig         = false;
      bool isXTrigLep      = true;
      bool isXTrigTau      = true;
    
      otree->trg_singlemuon = false;
      otree->trg_singleelectron = false;
      otree->singleLepTrigger = false;
      otree->ditauTrigger = false;
      otree->xTrigger = false;
      otree->xTriggerLep = false;
      otree->xTriggerTau = false;
    
      // setting weights to 1
      otree->trkeffweight = 1;
      otree->trigweight_1 = 1;
      otree->trigweight_2 = 1;
      otree->effweight = 1;
      otree->trigweight = 1;
      otree->idisoweight_1 = 1;
      otree->idisoweight_2 = 1;
    
      for (unsigned int iT = 0; iT < analysisTree.trigobject_count; ++iT) {
         float dRtrigLep = deltaR(lep_eta, lep_phi, analysisTree.trigobject_eta[iT], analysisTree.trigobject_phi[iT]);        
         float dRtrigTau = deltaR(analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex], analysisTree.trigobject_eta[iT], analysisTree.trigobject_phi[iT]);        
    
         if (dRtrigLep < dRTrigMatch){
           for(unsigned int i_trig = 0; i_trig < filterSingleLep.size(); i_trig++)
           {
              if (nSingleLepTrig.at(i_trig) == -1) continue;
              if (analysisTree.trigobject_filters[iT][nSingleLepTrig.at(i_trig)]) isSingleLepTrig = true;
            }
           for(unsigned int i_trig = 0; i_trig < filterXtriggerLepLeg.size(); i_trig++)
           {
              if (nXTrigLepLeg.at(i_trig) == -1) continue;
              if (analysisTree.trigobject_filters[iT][nXTrigLepLeg.at(i_trig)]) isXTrigLepLeg.at(i_trig) = true;
            } 
          }
         if (dRtrigTau < dRTrigMatch){ 
            for(unsigned int i_trig = 0; i_trig < filterXtriggerTauLeg.size(); i_trig++)
            {
               if (nXTrigTauLeg.at(i_trig) == -1) continue;
               if (analysisTree.trigobject_filters[iT][nXTrigTauLeg.at(i_trig)]) isXTrigTauLeg.at(i_trig) = true;
            }  
         }	  
      }
    
      for(unsigned int i_trig = 0; i_trig < filterXtriggerTauLeg.size(); i_trig++)
        isXTrigTau = isXTrigTau && isXTrigTauLeg.at(i_trig);
      for(unsigned int i_trig = 0; i_trig < filterXtriggerLepLeg.size(); i_trig++)
        isXTrigLep = isXTrigLep && isXTrigLepLeg.at(i_trig);
      isXTrig = isXTrigTau && isXTrigLep;
    
      if (ch == "mt")
       otree->trg_singlemuon = isSingleLepTrig;
      if (ch == "et")
       otree->trg_singleelectron = isSingleLepTrig;
    
      otree->singleLepTrigger = isSingleLepTrig;
      otree->xTriggerLep = isXTrigLep;
      otree->xTriggerTau = isXTrigTau;
      otree->xTrigger = isXTrig;
    
    
      //filling variables
      TLorentzVector leptonLV;
    
      // used for trigger weights
      double eff_data_trig_lt_tau = 1;
      double eff_mc_trig_lt_tau   = 1;
      double eff_data_trig_lt_l   = 1;
      double eff_mc_trig_lt_l     = 1;
      double eff_data_trig_L      = 1;
      double eff_mc_trig_L        = 1;
      double sf_trig_ditau_tau1   = 1;
      double sf_trig_ditau_tau2   = 1;
      // reset efficiency weights
    
      //all criterua passed, we fill vertices here;	
      SaveRECOVertices(&analysisTree,otree, isData);
    
      //Merijn: save here all gen information for the selected RECO events, gets stored for convenience in the taucheck tree ;-). Note that no selection on gen level is applied..     
      //Merijn 2019 4 3:note that a separate fill is not needed. We store in the otree now, which is Filled at the bottom! Filling here will make things out of synch..
      if (!isData)
        FillGenTree(&analysisTree, gentreeForGoodRecoEvtsOnly);       
    
      if(ch == "mt") {
      	FillMuTau(&analysisTree, otree, leptonIndex, dRiso);
        leptonLV.SetXYZM(analysisTree.muon_px[leptonIndex], analysisTree.muon_py[leptonIndex], analysisTree.muon_pz[leptonIndex], muonMass);
    
      	// tracking efficiency weight	
        if (!isData && ApplyLepSF) {
           w->var("m_eta")->setVal(analysisTree.muon_eta[leptonIndex]); 
           otree->trkeffweight = (double)(w->function("m_trk_ratio")->getVal());
        }
      } 
      else if(ch == "et"){
      	FillETau(&analysisTree, otree, leptonIndex, dRiso);
        leptonLV.SetXYZM(analysisTree.electron_px[leptonIndex], analysisTree.electron_py[leptonIndex], analysisTree.electron_pz[leptonIndex], electronMass);
    
        // tracking efficiency weight
        if (!isData && ApplyLepSF) {
      	  w->var("e_eta")->setVal(analysisTree.electron_eta[leptonIndex]); 
      	  w->var("e_pt")->setVal(analysisTree.electron_pt[leptonIndex]); 	
      	  otree->trkeffweight = (double)( w->function("e_reco_ratio")->getVal());
         }
      }
    
      if (!isData && ApplyLepSF) {
    	  if (ch == "mt") {
    	    eff_data_trig_lt_tau = tauTriggerSF->getMuTauEfficiencyData(analysisTree.tau_pt[tauIndex], analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex]);
    	    eff_mc_trig_lt_tau = tauTriggerSF->getMuTauEfficiencyMC(analysisTree.tau_pt[tauIndex], analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex]);
    	  }
    	  else {
    	    eff_data_trig_lt_tau = tauTriggerSF->getETauEfficiencyData(analysisTree.tau_pt[tauIndex], analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex]);
    	    eff_mc_trig_lt_tau = tauTriggerSF->getETauEfficiencyMC(analysisTree.tau_pt[tauIndex], analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex]);
    	  }
    
    	  eff_data_trig_lt_l = SF_XTriggerLepLeg->get_EfficiencyData(leptonLV.Pt(), leptonLV.Eta());
    	  eff_mc_trig_lt_l = SF_XTriggerLepLeg->get_EfficiencyMC(leptonLV.Pt(), leptonLV.Eta());
    	  eff_data_trig_L = SF_SingleLepTrigger->get_EfficiencyData(leptonLV.Pt(), leptonLV.Eta());
    	  eff_mc_trig_L = SF_SingleLepTrigger->get_EfficiencyData(leptonLV.Pt(), leptonLV.Eta());
    
    	  otree->idisoweight_1 = SF_lepIdIso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
    	  otree->idisoweight_antiiso_1 = SF_lepIdIso_antiiso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
    	  otree->trigweight_1 = SF_SingleLepTrigger->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
    	  otree->trigweight_antiiso_1 = SF_SingleLepTrigger_antiiso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
    
    	  double eff_data_trig = eff_data_trig_L + (eff_data_trig_lt_l - eff_data_trig_L) * eff_data_trig_lt_tau;
    	  double eff_mc_trig = eff_mc_trig_L + (eff_mc_trig_lt_l - eff_mc_trig_L) * eff_mc_trig_lt_tau;
    
    	  if (eff_data_trig > 1e-4 && eff_mc_trig > 1e-4)
    	    otree->trigweight = eff_data_trig / eff_mc_trig;
      }
      counter[10]++;
    
      FillTau(&analysisTree, otree, tauIndex);
    
      if (!isData && analysisTree.tau_genmatch[tauIndex] == 5) otree->idisoweight_2 = tau_id_sf; 
    
      otree->effweight = otree->idisoweight_1 * otree->trkeffweight * otree->idisoweight_2 * otree->trigweight;
      otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
    
      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
      //MET
      //Merijn 2019 6 20: overloaded the function, it takes the era as arugment now, to take pfmetcorr for 2016 and 2017..
      fillMET(ch, leptonIndex, tauIndex, &analysisTree, otree, era);
      
      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);
      
      // Zpt weight
      otree->zptweight = 1.;
      if (!isData && ((isDY && isMG ) || isEWKZ)){
        genV = genTools::genV(analysisTree); // gen Z boson ?
      	float bosonMass = genV.M();
      	float bosonPt = genV.Pt();

        //Merijn determine here some min and max values:
        double massxmin = h_zptweight->GetXaxis()->GetXmin();
        double massxmax = h_zptweight->GetXaxis()->GetXmax();

        double ptxmin = h_zptweight->GetYaxis()->GetXmin();
        double ptxmax = h_zptweight->GetYaxis()->GetXmax();

      	//Merijn 2019 6 13: adjust to T/M functions, to get boundaries right. Otherwise, for 2017 data we get few outliers that screw up the weight histogram dramatically.
      	Float_t zptmassweight = 1;
      	if (bosonMass > 50.0) {
          float bosonMassX = bosonMass;
          float bosonPtX = bosonPt;
          if (bosonMassX > massxmax) bosonMassX = massxmax - h_zptweight->GetXaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;//Merijn: if doesn't work, lower by 1/2 binwidth..
          if (bosonPtX < ptxmin)     bosonPtX = ptxmin + h_zptweight->GetYaxis()->GetBinWidth(1)*0.5;
          if (bosonPtX > ptxmax)     bosonPtX = ptxmax - h_zptweight->GetYaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;
          zptmassweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(bosonMassX), h_zptweight->GetYaxis()->FindBin(bosonPtX));
          }	
          otree->zptweight = zptmassweight;
      }
      
      // topPt weight
      otree->topptweight = 1.;
      int a_topPtWeight = cfg.get<int>("a_topPtWeight");
      int b_topPtWeight = cfg.get<int>("b_topPtWeight");
      if(!isData)
         // otree->topptweight = genTools::return_topPtWeight(analysisTree, a_topPtWeight, b_topPtWeight);
         otree->topptweight = genTools::topPtWeight(analysisTree, 1); // 1 is for Run1 - use this reweighting as recommended by HTT 17
      counter[11]++;
      
      // lepton tau fakerates
      otree->mutaufakeweight = 1.;
      otree->etaufakeweight = 1.;
      if (!isData){
      	if (ch == "et") {
      	  otree->etaufakeweight = leptauFR->get_fakerate("electron", "Tight", otree->eta_2, otree->gen_match_2);
      	  otree->mutaufakeweight = leptauFR->get_fakerate("muon", "Loose", otree->eta_2, otree->gen_match_2);
      	}
    	  else if (ch == "mt") {
      	  otree->etaufakeweight = leptauFR->get_fakerate("electron", "VLoose", otree->eta_2, otree->gen_match_2);
      	  otree->mutaufakeweight = leptauFR->get_fakerate("muon", "Tight", otree->eta_2, otree->gen_match_2);
      	}
      }
    
      ////////////////////////////////////////////////////////////
      // MET Recoil Corrections
      ////////////////////////////////////////////////////////////
      
      otree->njetshad = otree->njets;
      
      if (!isData && ApplyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal) ){
      	genV = genTools::genV(analysisTree);
      	genL = genTools::genL(analysisTree);
      	if(isWJets) otree->njetshad += 1;

        genTools::RecoilCorrections( *recoilPFMetCorrector, 1, // dummy parameter
          otree->met, otree->metphi,
          genV.Px(), genV.Py(),
          genL.Px(), genL.Py(),
          otree->njetshad,
          otree->met_rcmr, otree->metphi_rcmr
        );
        
        // overwriting with recoil-corrected values 
        otree->met = otree->met_rcmr;
        otree->metphi = otree->metphi_rcmr;   
      }
      
      //ditau sytem
      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
  				     analysisTree.tau_py[tauIndex],
  				     analysisTree.tau_pz[tauIndex],
  				     analysisTree.tau_mass[tauIndex]);
    
      // using PF MET
      TLorentzVector metLV; 
      metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi), otree->met*TMath::Sin(otree->metphi), 0,
                    TMath::Sqrt( otree->met*TMath::Sin(otree->metphi)*otree->met*TMath::Sin(otree->metphi) +
  			               otree->met*TMath::Cos(otree->metphi)*otree->met*TMath::Cos(otree->metphi)));
    
    
      // shift the tau energy scale by decay mode and propagate to the met. 
      if (!isData) {
      	bool isOneProng = false;
      	float shift_tes = 0.0;
      	if (otree->gen_match_2 >= 5){
      	  if (otree->tau_decay_mode_2 == 0){
            shift_tes = shift_tes_1prong; 
            isOneProng = true;
          }
      	  else if (otree->tau_decay_mode_2 == 1) shift_tes = shift_tes_1p1p0; 
      	  else if (otree->tau_decay_mode_2 == 10) shift_tes = shift_tes_3prong;
         }
      	else if (otree->gen_match_2 < 5) {
      	  if (otree->tau_decay_mode_2 == 0)      {shift_tes = shift_tes_lepfake_1prong; isOneProng = true;}
      	  else if (otree->tau_decay_mode_2 == 1)  shift_tes = shift_tes_lepfake_1p1p0; 
      	  else if (otree->tau_decay_mode_2 == 10) shift_tes = shift_tes_lepfake_3prong; 
      	}
      	correctTauES(tauLV, metLV, shift_tes, isOneProng);	    
      	otree->pt_2 = tauLV.Pt();
      	otree->m_2 = tauLV.M();
      	otree->met = metLV.Pt();
      	otree->metphi = metLV.Phi();
       }
    
      if (!isData) {
      	if(otree->gen_match_2 == 5 && tauLV.E() <= 400 && tauLV.E() >= 20){
      	  if (otree->tau_decay_mode_2 == 0) tauLV *= (1-0.03);
      	  else if (otree->tau_decay_mode_2 < 5) tauLV *= (1-0.02);
      	  else if (otree->tau_decay_mode_2 == 10)tauLV *= (1-0.01);
      	  otree->pt_2 = tauLV.Pt();
      	  otree->m_2 = tauLV.M();
      	}
      }
    
      TLorentzVector dileptonLV = leptonLV + tauLV;
      otree->m_vis = dileptonLV.M();
      otree->pt_tt = (dileptonLV+metLV).Pt();   
    
      // mt TOT
      float mtTOT = 2*(otree->pt_1)*metLV.Pt()*(1-cos(DeltaPhi(leptonLV,metLV)));
      mtTOT += 2*(otree->pt_2)*metLV.Pt()*(1-cos(DeltaPhi(tauLV,metLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(leptonLV,tauLV))); 
      otree->mt_tot = TMath::Sqrt(mtTOT);
    
      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;
    
      // dilepton veto
      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
      if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree);
    
  	  //extra lepton veto
      otree->extraelec_veto = extra_electron_veto(leptonIndex, ch, &cfg, &analysisTree);
      otree->extramuon_veto = extra_muon_veto(leptonIndex, ch, &cfg, &analysisTree, isData);
    
      counter[13]++;
    
      otree->mt_1 = mT(leptonLV, metLV);
      otree->mt_2 = mT(tauLV, metLV);
    
      // bisector of lepton and tau transverse momenta
    
      float leptonUnitX = leptonLV.Px() / leptonLV.Pt();
      float leptonUnitY = leptonLV.Py() / leptonLV.Pt();
    
      float tauUnitX = tauLV.Px() / tauLV.Pt();
      float tauUnitY = tauLV.Py() / tauLV.Pt();
    
      float zetaX = leptonUnitX + tauUnitX;
      float zetaY = leptonUnitY + tauUnitY;
    
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
    
      zetaX = zetaX / normZeta;
      zetaY = zetaY / normZeta;
    
      float vectorVisX = leptonLV.Px() + tauLV.Px();
      float vectorVisY = leptonLV.Py() + tauLV.Py();
    
      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;
      otree->pzetamiss = otree->met*TMath::Cos(otree->metphi)*zetaX + otree->met*TMath::Sin(otree->metphi)*zetaY;
      counter[14]++;
    
      // svfit variables
      otree->m_sv   = -10;//Merijn updated for the DNN
      otree->pt_sv  = -9999;
      otree->eta_sv = -9999;
      otree->phi_sv = -9999;
      otree->met_sv = -9999;
      otree->mt_sv = -9999;
    
      //calculate SV fit only for events passing baseline selection and mt cut
      // fill otree only for events passing baseline selection 
      // for synchronisation, take all events
      const bool Synch = cfg.get<bool>("Synch"); 
      bool passedBaselineSel = false;
    
  	  // for SM analysis
      if (ch == "mt") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5 && 
  		      otree->againstElectronVLooseMVA6_2>0.5 && otree->againstMuonTight3_2>0.5  &&
  		      otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);
      if (ch == "et") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5 && 
                            otree->againstMuonLoose3_2>0.5 && otree->againstElectronTightMVA6_2>0.5 && 
                            otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);
    
      if(otree->iso_1<0.35 && 
    	 otree->againstMuonLoose3_2 > 0.5 && otree->againstElectronTightMVA6_2 > 0.5 && 
    	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[15]++;
      if(otree->iso_1 < 0.35 && otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 > 0.5 && 
    	 otree->againstElectronTightMVA6_2 > 0.5 && 
    	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[16]++;
      if(otree->iso_1 < 0.35 && otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 > 0.5 && 
    	 otree->againstMuonLoose3_2 > 0.5 && 
    	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[17]++;
      if(otree->iso_1 < 0.35 && otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 > 0.5 && 
      	 otree->againstMuonLoose3_2 > 0.5 && otree->againstElectronTightMVA6_2 > 0.5 ) counter[18]++;
    
      //if (!Synch && !passedBaselineSel) continue;
    
      if (ApplySVFit && otree->njetspt20 > 0) svfit_variables(ch, &analysisTree, otree, &cfg, inputFile_visPtResolution);
    
    
      //addition Merijn: here we select the constituent of the tau with highest pT
      int ncomponents = analysisTree.tau_constituents_count[tauIndex];
      float maxPt = -1;
      int sign = -1;
      int pdgcode = -9999;
      if(analysisTree.tau_charge[tauIndex] > 0) sign = 1; 
    
      for(int i = 0; i < ncomponents; i++){  
        if((analysisTree.tau_constituents_pdgId[tauIndex][i]*sign) > 0){
          TLorentzVector lvector; 
          lvector.SetXYZT(analysisTree.tau_constituents_px[tauIndex][i],
      				      analysisTree.tau_constituents_py[tauIndex][i],
      				      analysisTree.tau_constituents_pz[tauIndex][i],
      				      analysisTree.tau_constituents_e[tauIndex][i]);
          double Pt = lvector.Pt();
          if(Pt > maxPt){
          	pdgcode = analysisTree.tau_constituents_pdgId[tauIndex][i];
          	maxPt = Pt;
          }
        }
      }
    
      otree->pdgcodetau2 = pdgcode; //Merijn tried here to assign to our tree. Not working yet so put in histogam..
      ConstitsPDG->Fill(pdgcode);
      if(abs(pdgcode) != 211 && abs(pdgcode) != 22) nonpionphotonctr++;
    
        // evaluate systematics for MC 
      if(!isData && ApplySystShift){
      zPtWeightSys->Eval(); 
      topPtWeightSys->Eval();
      for(unsigned int i = 0; i < jetEnergyScaleSys.size(); i++)
        (jetEnergyScaleSys.at(i))->Eval(); 
      if (ch == "mt") tauScaleSys->Eval(utils::MUTAU);
      else if (ch == "et") tauScaleSys->Eval(utils::ETAU);
      }
      counter[19]++;
    
    
      //CP calculation. Updates Merijn: placed calculation at end, when all kinematic corrections are performed. Removed statement to only do calculation for tt. 
      //Created the acott_Impr function, which takes ch as input as well. See the funcrtion in functionsCP.h to see my updates to the function itself      
      //Merijn 2019 1 10 debug: a major source of problems was that indices were innertwined from the beginning...
      //one should note that in et or mt case,
      acott_Impr(&analysisTree, otree, leptonIndex, tauIndex, ch);
      selEvents++;
    
      //Merijn 2019 1 10: perhaps this should be called before moving to next event..
      otree->Fill();
    } // event loop

    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  } // file loop
   
  std::cout << "COUNTERS" << std::endl;
  for(int iC = 0; iC < 20; iC++) std::cout << "Counter " << iC << ":    " << counter[iC] << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of input events    = " << int(inputEventsH->GetEntries()) << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  cout<<"nonpionphotonctr ="<<nonpionphotonctr<<endl;

  file->cd("");
  file->Write();
  ConstitsPDG->Write();

  // delete systematics objects

  if(tauScaleSys != 0){
    tauScaleSys->Write();
    delete tauScaleSys;
  }

  if(tauOneProngScaleSys != 0){
    tauOneProngScaleSys->Write();
    delete tauOneProngScaleSys;
  }

  if(tauOneProngOnePi0ScaleSys != 0){
    tauOneProngOnePi0ScaleSys->Write();
    delete tauOneProngOnePi0ScaleSys;
  }

  if(tauThreeProngScaleSys != 0){
    tauThreeProngScaleSys->Write();
    delete tauThreeProngScaleSys;
  }

  if(zPtWeightSys != 0){
    zPtWeightSys->Write();
    delete zPtWeightSys;
  }

  if(topPtWeightSys != 0){
    topPtWeightSys->Write();
    delete topPtWeightSys;
  }

  if(jetEnergyScaleSys.size() > 0){
    for (unsigned int i = 0; i < jetEnergyScaleSys.size(); i++){
      (jetEnergyScaleSys.at(i))->Write();
      delete jetEnergyScaleSys.at(i);
    }
  }

  if(lepTauFakeOneProngScaleSys != 0){
    lepTauFakeOneProngScaleSys->Write();
    delete lepTauFakeOneProngScaleSys;
  }

  if(lepTauFakeOneProngOnePi0ScaleSys != 0){
    lepTauFakeOneProngOnePi0ScaleSys->Write();
    delete lepTauFakeOneProngOnePi0ScaleSys;
  }

  if(lepTauFakeThreeProngScaleSys != 0){
    lepTauFakeThreeProngScaleSys->Write();
    delete lepTauFakeThreeProngScaleSys;
  }

  file->Close();
  delete file;

}


////FILLING FUNCTIONS//////

void SaveRECOVertices(const AC1B *analysisTree, Synch17Tree *otree, const bool isData){

  otree->RecoVertexX = analysisTree->primvertex_x;
  otree->RecoVertexY = analysisTree->primvertex_y;
  otree->RecoVertexZ = analysisTree->primvertex_z;


  if(!isData){
    for (unsigned int igen = 0; igen < analysisTree->genparticles_count; ++igen) {

  //here fill the generator vertices to have the gen information present in tree PER GOOD RECO EVENT
  //Note: we may want to add constraint that the W and Z are prompt. If we remove these, may get in trouble with a DY or W MC sample..

      if ((analysisTree->genparticles_pdgid[igen] == 23 || analysisTree->genparticles_pdgid[igen] == 24 ||
  	       analysisTree->genparticles_pdgid[igen] == 25 || analysisTree->genparticles_pdgid[igen] == 35 || analysisTree->genparticles_pdgid[igen] == 36) && 
           analysisTree->genparticles_isLastCopy[igen] == 1 && analysisTree->genparticles_isPrompt[igen] == 1) {
        otree->GenVertexX = analysisTree->genparticles_vx[igen];
        otree->GenVertexY = analysisTree->genparticles_vy[igen];
        otree->GenVertexZ = analysisTree->genparticles_vz[igen];
        break;
      }
    }
  }
  else {//if it is data, fill with something recognisable nonsensible
    otree->GenVertexX = -9999;
    otree->GenVertexY = -9999;
    otree->GenVertexZ = -9999;
  }
}

void initializeGenTree(Synch17GenTree *gentree){
  gentree->Higgs_pt=-9999;
  gentree->Higgs_eta=-9999;
  gentree->Higgs_phi=-9999;
  gentree->Higgs_mass=-9999;
  gentree->pt_1=-9999;
  gentree->eta_1=-9999;
  gentree->phi_1=-9999;
  gentree->pt_2=-9999;
  gentree->eta_2=-9999;
  gentree->phi_2=-9999;
  gentree->acotautau_00 = -9999;
  gentree->acotautau_10 = -9999;
  gentree->acotautau_01 = -9999;
  gentree->acotautau_11 = -9999;
  gentree->acotautau_02=-9999;
  gentree->acotautau_20=-9999;
  gentree->acotautau_12=-9999;
  gentree->acotautau_21=-9999;
  gentree->acotautau_22=-9999;

   //Merijn added the angle psi, currently for debugging purpose. Later may extend to 3-prong..
  gentree->acotautauPsi_00=-9999;
  gentree->acotautauPsi_01=-9999;
  gentree->acotautauPsi_10=-9999;
  gentree->acotautauPsi_11=-9999;

  //init the new vertex variables to something nonsensible also
  gentree->VertexX=-9999;
  gentree->VertexY=-9999;
  gentree->VertexZ=-99999;

  gentree->VxConstitTau1=-9999;
  gentree->VyConstitTau1=-9999;
  gentree->VzConstitTau1=-9999;
  
  gentree->VxConstitTau2=-9999;
  gentree->VyConstitTau2=-9999;
  gentree->VzConstitTau2=-9999;

  //Merijn add initialiser for
  gentree->chconst_1_pt=-9999;
  gentree->chconst_1_eta=-9999;
  gentree->chconst_1_phi=-9999;
  
  gentree->chconst_2_pt=-9999;
  gentree->chconst_2_eta=-9999;
  gentree->chconst_2_phi=-9999;
  gentree->alphaminus=-9999;

}
void initializeCPvar(Synch17Tree *otree){
  otree->acotautau_00=-9999;
  otree->acotautau_01=-9999;
  otree->acotautau_10=-9999;
  otree->acotautau_11=-9999;
  /*
  otree->acotautau_20=-9999;
  otree->acotautau_02=-9999;
  otree->acotautau_21=-9999;
  otree->acotautau_12=-9999;
  otree->acotautau_22=-9999;
  */

  //Merijn added the angle psi, currently for debugging purpose. Later may extend to 3-prong..
  otree->acotautauPsi_00=-9999;
  otree->acotautauPsi_01=-9999;
  otree->acotautauPsi_10=-9999;
  otree->acotautauPsi_11=-9999;
  
  otree->tau1DecayPlaneX=-9999;
  otree->tau1DecayPlaneY=-9999;
  otree->tau1DecayPlaneZ=-9999;
  otree->tau2DecayPlaneX=-9999;
  otree->tau2DecayPlaneY=-9999;
  otree->tau2DecayPlaneZ=-9999;

  otree->VxConstitTau1=-9999;
  otree->VyConstitTau1=-9999;
  otree->VzConstitTau1=-9999;
  
  otree->VxConstitTau2=-9999;
  otree->VyConstitTau2=-9999;
  otree->VzConstitTau2=-9999;

  //Merijn add initialiser for
  otree->chconst_1_pt=-9999;
  otree->chconst_1_eta=-9999;
  otree->chconst_1_phi=-9999;
  
  otree->chconst_2_pt=-9999;
  otree->chconst_2_eta=-9999;
  otree->chconst_2_phi=-9999;
  otree->alphaminus=-9999;

}

void FillGenTree(const AC1B *analysisTree, Synch17GenTree *gentree){
  int ntaus=analysisTree->gentau_count;
  int npart=analysisTree->genparticles_count;
  int leptonid=15;
  TLorentzVector Tau1,Tau2,Tau;
  TLorentzVector Lepton;
  TLorentzVector lvector;
  int tauHIndex=-1;
  int tauLIndex=-1;
  int LeadingtauIndex=-1;
  int TrailingtauIndex=-1;
  double taumaxpt=-1;
  
  for(int itau=0;itau<ntaus;itau++){
    if(analysisTree->gentau_isLastCopy[itau]==1&&analysisTree->gentau_isPrompt[itau]==1){
      if(analysisTree->gentau_visible_pt[itau]>=taumaxpt) {
	LeadingtauIndex=itau; 
	taumaxpt=analysisTree->gentau_visible_pt[itau];
      }
    }
  }

  taumaxpt=-1; 
  for(int itau=0;itau<ntaus;itau++){
    if(analysisTree->gentau_isLastCopy[itau]==1&&analysisTree->gentau_isPrompt[itau]==1&&itau!=LeadingtauIndex){
      if(analysisTree->gentau_visible_pt[itau]>=taumaxpt) {
	TrailingtauIndex=itau; 
	taumaxpt=analysisTree->gentau_visible_pt[itau];
      }
    }
  }

  TLorentzVector genTauVis1; genTauVis1.SetXYZT(0,0,0,0);
  TLorentzVector genTauVis2; genTauVis2.SetXYZT(0,0,0,0);
  gentree->decaymode_1 = -1;
  gentree->decaymode_2 = -1;
  if (LeadingtauIndex>-1) {
    genTauVis1.SetXYZT(analysisTree->gentau_visible_px[LeadingtauIndex],
		       analysisTree->gentau_visible_py[LeadingtauIndex],
		       analysisTree->gentau_visible_pz[LeadingtauIndex],
		       analysisTree->gentau_visible_e[LeadingtauIndex]);
    gentree->decaymode_1 = analysisTree->gentau_decayMode[LeadingtauIndex];
  }
  if (TrailingtauIndex>-1) {
    genTauVis2.SetXYZT(analysisTree->gentau_visible_px[TrailingtauIndex],
		       analysisTree->gentau_visible_py[TrailingtauIndex],
		       analysisTree->gentau_visible_pz[TrailingtauIndex],
		       analysisTree->gentau_visible_e[TrailingtauIndex]);
    gentree->decaymode_2 = analysisTree->gentau_decayMode[TrailingtauIndex];
  }
  gentree->pt_1 = genTauVis1.Pt();
  gentree->eta_1 = genTauVis1.Eta();
  gentree->phi_1 = genTauVis1.Phi();

  gentree->pt_2 = genTauVis2.Pt();
  gentree->eta_2 = genTauVis2.Eta();
  gentree->phi_2 = genTauVis2.Phi();

  double dR;
  const double dRcut=0.3;
  for(int ipart=0;ipart<npart;ipart++){
    if((abs(analysisTree->genparticles_pdgid[ipart])==25||
	abs(analysisTree->genparticles_pdgid[ipart])==35||
	abs(analysisTree->genparticles_pdgid[ipart])==36)&&
       analysisTree->genparticles_isLastCopy[ipart]==1){
      TLorentzVector Higgs;
      Higgs.SetPxPyPzE(analysisTree->genparticles_px[ipart],
		       analysisTree->genparticles_py[ipart],
		       analysisTree->genparticles_pz[ipart],
		       analysisTree->genparticles_e[ipart]);
      gentree->Higgs_pt=Higgs.Pt();
      gentree->Higgs_eta=Higgs.Eta();
      gentree->Higgs_phi=Higgs.Phi();
      gentree->Higgs_mass=Higgs.M();
    }
  }

  
  if (LeadingtauIndex>-1&&TrailingtauIndex>-1)
    gen_acott(analysisTree,gentree,LeadingtauIndex,TrailingtauIndex);

  gentree->a1polarization_1=gen_A1Polarization(analysisTree,LeadingtauIndex);
  gentree->a1polarization_2=gen_A1Polarization(analysisTree,TrailingtauIndex);

//here fill the generator vertices to have the information present in tree
//Note: we may want to add constraint that the W and Z are prompt. If we remove these, may get in trouble with a DY or W MC sample..


  for (unsigned int igen=0; igen<analysisTree->genparticles_count; ++igen) {
    if ((analysisTree->genparticles_pdgid[igen]==23||analysisTree->genparticles_pdgid[igen]==24||
	analysisTree->genparticles_pdgid[igen]==25||analysisTree->genparticles_pdgid[igen]==35||analysisTree->genparticles_pdgid[igen]==36)&&analysisTree->genparticles_isLastCopy[igen]==1&&analysisTree->genparticles_isPrompt[igen]==1) {
      gentree->VertexX=analysisTree->genparticles_vx[igen];
      gentree->VertexY=analysisTree->genparticles_vy[igen];
      gentree->VertexZ=analysisTree->genparticles_vz[igen];
      break;
    }
  }
}

//fill the otree with the muon variables in channel mutau
void FillMuTau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->muon_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->muon_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->muon_phi[leptonIndex];
  otree->m_1 = 		muonMass;
  otree->q_1 = -1;
  if (analysisTree->muon_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->muon_genmatch[leptonIndex];

  otree->iso_1 = abs_Iso_mt(leptonIndex, analysisTree, dRiso) / analysisTree->muon_pt[leptonIndex];

  otree->d0_1 = analysisTree->muon_dxy[leptonIndex];
  otree->dZ_1 = analysisTree->muon_dz[leptonIndex];
  //otree->d0err_1 = analysisTree->muon_dxyerr[leptonIndex];
  //otree->dZerr_1 = analysisTree->muon_dzerr[leptonIndex];
  

  otree->tau_decay_mode_1 = -9999; 
  // otree->tau_decay_mode_1=analysisTree->tau_decayMode[leptonIndex]; can;'t do since its a lepton not a tau index
 
  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -9999;
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->againstMuonLoose3_1 = -9999;
  otree->againstMuonTight3_1 = -9999;
  otree->againstElectronVLooseMVA6_1 = -9999;
  otree->againstElectronTightMVA6_1 = -9999;

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree-> byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = -9999;
  otree->chargedIsoPtSum_1 = -9999;
  otree->neutralIsoPtSum_1 = -9999;
  otree->puCorrPtSum_1 = -9999;

}

//fill the otree with the electron variables in channel etau
void FillETau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->electron_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->electron_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->electron_phi[leptonIndex];
  otree->m_1 = 		electronMass;
    otree->q_1 = -1;
  if (analysisTree->electron_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->electron_genmatch[leptonIndex];

  otree->iso_1 =  abs_Iso_et(leptonIndex, analysisTree, dRiso) /analysisTree->electron_pt[leptonIndex];

  otree->d0_1 = analysisTree->electron_dxy[leptonIndex];
  otree->dZ_1 = analysisTree->electron_dz[leptonIndex];
  //otree->d0err_1 = analysisTree->electron_dxyerr[leptonIndex];
  //otree->dZerr_1 = analysisTree->electron_dzerr[leptonIndex]; 

  otree->tau_decay_mode_1 = -9999;
  // otree->tau_decay_mode_1=analysisTree->tau_decayMode[leptonIndex];// can;t do since its a lepton index..

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -9999;
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->againstMuonLoose3_1 = -9999;
  otree->againstMuonTight3_1 = -9999;
  otree->againstElectronVLooseMVA6_1 = -9999;
  otree->againstElectronTightMVA6_1 = -9999;

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree-> byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = -9999;
  otree->chargedIsoPtSum_1 = -9999;
  otree->neutralIsoPtSum_1 = -9999;
  otree->puCorrPtSum_1 = -9999;
  

}

//fill the otree with the tau variables 
void FillTau(const AC1B *analysisTree, Synch17Tree *otree, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->eta_2 = analysisTree->tau_eta[tauIndex];
  otree->phi_2 = analysisTree->tau_phi[tauIndex];
  otree->q_2 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_2 = analysisTree->tau_genmatch[tauIndex];
  otree->mva_2 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->mva17_2= analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->d0_2 = analysisTree->tau_leadchargedhadrcand_dxy[tauIndex];
  otree->dZ_2 = analysisTree->tau_leadchargedhadrcand_dz[tauIndex];      
  otree->iso_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->m_2 = analysisTree->tau_mass[tauIndex];
  otree->tau_decay_mode_2 = analysisTree->tau_decayMode[tauIndex];

  otree->tau_pca2D_x_2 = analysisTree->tau_pca2D_x[tauIndex];
  otree->tau_pca2D_y_2 = analysisTree->tau_pca2D_y[tauIndex];
  otree->tau_pca2D_z_2 = analysisTree->tau_pca2D_z[tauIndex];
  otree->tau_pca3D_x_2 = analysisTree->tau_pca3D_x[tauIndex];
  otree->tau_pca3D_y_2 = analysisTree->tau_pca3D_y[tauIndex];
  otree->tau_pca3D_z_2 = analysisTree->tau_pca3D_z[tauIndex];
  otree->tau_SV_x_2 = analysisTree->tau_SV_x[tauIndex];
  otree->tau_SV_y_2 = analysisTree->tau_SV_y[tauIndex];
  otree->tau_SV_z_2 = analysisTree->tau_SV_z[tauIndex];
  otree->tau_SV_covxx_2 = analysisTree->tau_SV_cov[tauIndex][0];
  otree->tau_SV_covyx_2 = analysisTree->tau_SV_cov[tauIndex][1];
  otree->tau_SV_covzx_2 = analysisTree->tau_SV_cov[tauIndex][2];
  otree->tau_SV_covyy_2 = analysisTree->tau_SV_cov[tauIndex][3];
  otree->tau_SV_covzy_2 = analysisTree->tau_SV_cov[tauIndex][4];
  otree->tau_SV_covzz_2 = analysisTree->tau_SV_cov[tauIndex][5];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstMuonLoose3_2 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_2 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronVLooseMVA6_2 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronTightMVA6_2 = analysisTree->tau_againstElectronTightMVA6[tauIndex];

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree-> byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex];
  otree->chargedIsoPtSum_2 = analysisTree->tau_chargedIsoPtSum[tauIndex];
  otree->neutralIsoPtSum_2 = analysisTree->tau_neutralIsoPtSum[tauIndex];
  otree->puCorrPtSum_2 = analysisTree->tau_puCorrPtSum[tauIndex];
  //otree->isolationGammaCands_size_2 = analysisTree->tau_isolationGammaCands_size[tauIndex];
  //otree->signalGammaCands_size_2 = analysisTree->tau_signalGammaCands_size[tauIndex];

  otree->correction_againstElectronVLooseMVA6_2 = 1;
  //otree->correction_againstElectronLooseMVA6_2 = 1;
  //otree->correction_againstElectronMediumMVA6_2 = 1;
  otree->correction_againstElectronTightMVA6_2 = 1;
  //otree->correction_againstElectronVTightMVA6_2 = 1;
  otree->correction_againstMuonLoose3_2 = 1;
  otree->correction_againstMuonTight3_2 = 1;

  if (analysisTree->tau_genmatch[tauIndex]==1){
    if (abs(otree->eta_2)<1.460){
      otree->correction_againstElectronVLooseMVA6_2 = 1.09;
      //		otree->correction_againstElectronLooseMVA6_2 = 1.17;
      //		otree->correction_againstElectronMediumMVA6_2 = 1.40;
      otree->correction_againstElectronTightMVA6_2 = 1.80;
      //		otree->correction_againstElectronVTightMVA6_2 = 1.96;
    }else if (abs(otree->eta_2)>1.558){
      otree->correction_againstElectronVLooseMVA6_2 = 1.19;
      //		otree->correction_againstElectronLooseMVA6_2 = 1.25;
      //		otree->correction_againstElectronMediumMVA6_2 = 1.21;
      otree->correction_againstElectronTightMVA6_2 = 1.53;
      //		otree->correction_againstElectronVTightMVA6_2 = 1.66;
    }
  }
  
  if (analysisTree->tau_genmatch[tauIndex]==2){
    if ((0 < abs(otree->eta_2))&&(abs(otree->eta_2)<= 0.4)){
      otree->correction_againstMuonLoose3_2 = 1.06;
      otree->correction_againstMuonTight3_2 = 1.17;
    }else if ((0.4 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 0.8)){
      otree->correction_againstMuonLoose3_2 = 1.02;
      otree->correction_againstMuonTight3_2 = 1.14;
    }else if ((0.8 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 1.2)){
      otree->correction_againstMuonLoose3_2 = 1.10;
      otree->correction_againstMuonTight3_2 = 1.14;
    }else if ((1.2 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 1.7)){
      otree->correction_againstMuonLoose3_2 = 1.03;
      otree->correction_againstMuonTight3_2 = 0.93;
    }else if ((1.7 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 2.3)){
      otree->correction_againstMuonLoose3_2 = 1.94;
      otree->correction_againstMuonTight3_2 = 1.61;
    }
  }
  ///////////////////////////////////////////////////////NEW
  otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  if (analysisTree->tau_genmatch[tauIndex]==5){
    otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.88;
    //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.86;
    //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.84;
  }

}
