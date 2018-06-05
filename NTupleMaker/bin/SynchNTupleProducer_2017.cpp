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
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
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
void FillTau_leading(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex);
void initializeCPvar(Synch17Tree *otree);
//void fillTTbarUncWeights(const AC1B * analysisTree, Synch17Tree *otree, bool isData, bool includeTTbarUncWeights);

int main(int argc, char * argv[]){

  // first argument - config file for analysis
  // second argument - file list (MUST BE IN THE SAME DIRECTORY OF THE EXECUTABLE)
  // third argument - channel ("et" or "mt")
  // third argument - index of first file to run on (optional, ignored if only one file is used)
  // fourth argument - index of last file to run on (optional, ignored if only one file is used)

  using namespace std;

  gErrorIgnoreLevel = kFatal;
  //gDebug = 2;

  string cmsswBase = (getenv ("CMSSW_BASE"));
  
  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load( pathToCrystalLib );
  if (openSuccessful !=0 ) {
    cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit( -1 );
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
	std::cout << " Channel " << ch << " is not a valid choice. Please try again with 'mt' or 'et'.Exiting. " << std::endl;
    exit(0);
  }

  if (ch=="mt") lep = "Muon"; 
  if (ch=="et") lep = "Electron";
  if (ch=="tt") lep = "Tau"; 

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyLeptonSFfromKIT  = cfg.get<bool>("ApplyLeptonSFfromKIT"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters  = cfg.get<bool>("ApplyMetFilters");

  bool checkDuplicateMuons = false;
  if (isData && ch == "mt") checkDuplicateMuons = cfg.get<bool>("checkDuplicateMuons");

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
  	
  //svfit
  //const string svFitPtResFile = cfg.get<string>("svFitPtResFile");
  const string svFitPtResFile = TString(TString(cmsswBase)+"/src/"+TString(cfg.get<string>("svFitPtResFile"))).Data();

  //zptweight file 
  const string ZptweightFile = cfg.get<string>("ZptweightFile");

  //b-tag scale factors
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  TString pathToBtagScaleFactors = (TString) cmsswBase+"/src/"+BtagSfFile;
  //TString pathToBtagScaleFactors = (TString) cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv";
  if( ApplyBTagScaling && gSystem->AccessPathName(pathToBtagScaleFactors) ){
    cout<<pathToBtagScaleFactors<<" not found. Please check."<<endl;
    exit( -1 );
  }//cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv"
  BTagCalibration calib("csvv2", (string) pathToBtagScaleFactors );
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central");
  if( ApplyBTagScaling ){
    reader_B.load(calib,BTagEntry::FLAV_B,"comb");
    reader_C.load(calib,BTagEntry::FLAV_C,"comb");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");
  }
  const string TaggingEfficienciesFile = cfg.get<string>("BtagMCeffFile");
  TString pathToTaggingEfficiencies = (TString) cmsswBase+"/src/"+TaggingEfficienciesFile;
  if ( ApplyBTagScaling && gSystem->AccessPathName(pathToTaggingEfficiencies) ){
    cout<<pathToTaggingEfficiencies<<" not found. Please check."<<endl;
    exit( -1 );
  }
  TFile *fileTagging  = new TFile( pathToTaggingEfficiencies );
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  
  if( ApplyBTagScaling ){
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  }
  TRandom3 *rand = new TRandom3();

  const struct btag_scaling_inputs inputs_btag_scaling_medium = { reader_B, reader_C, reader_Light, tagEff_B, tagEff_C, tagEff_Light, rand };

  // MET Recoil Corrections
  const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  const bool isDY = infiles.find("DY") == infiles.rfind("/")+1;
  const bool isWJets = (infiles.find("WJets") == infiles.rfind("/")+1) || (infiles.find("W1Jets") == infiles.rfind("/")+1) || (infiles.find("W2Jets") == infiles.rfind("/")+1) || (infiles.find("W3Jets") == infiles.rfind("/")+1) || (infiles.find("W4Jets") == infiles.rfind("/")+1) || (infiles.find("EWK") == infiles.rfind("/")+1);
  const bool isVBForGGHiggs = (infiles.find("VBFHTo")== infiles.rfind("/")+1) || (infiles.find("GluGluHTo")== infiles.rfind("/")+1);
  const bool isEWKZ =  infiles.find("EWKZ") == infiles.rfind("/")+1;
  const bool isMG = infiles.find("madgraph") != string::npos;
  const bool isMSSMsignal =  (infiles.find("SUSYGluGluToHToTauTau")== infiles.rfind("/")+1) || (infiles.find("SUSYGluGluToBBHToTauTau")== infiles.rfind("/")+1);
  //const bool applyRecoilCorrections = isDY || isWJets;
  
  
  RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
  RecoilCorrector* recoilMvaMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilMvaMetCorrector));
  
  if(!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal) ){
    TString RecoilDir("HTT-utilities/RecoilCorrections/data/");
    
    TString RecoilFileName = RecoilDir; RecoilFileName += "TypeI-PFMet_Run2016BtoH.root";
    std::cout<<RecoilFileName<<std::endl;
    recoilPFMetCorrector = new RecoilCorrector( RecoilFileName);
        
  
    RecoilFileName = RecoilDir; RecoilFileName += "MvaMET_2016BCD.root";
    std::cout<<RecoilFileName<<std::endl;
    recoilMvaMetCorrector = new RecoilCorrector( RecoilFileName);
  }
  
  // Read in HLT filter
  vector<string> filterSingleLep;
  vector<string> filterXtriggerLepLeg;
  vector<string> filterXtriggerTauLeg;
  if(ch == "mt"){
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep1"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg1"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg2"));
  }else if(ch == "et"){
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep1"));
    
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg1"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg2"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg1"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg2"));
    
  }else if(ch == "tt"){
    filterSingleLep.push_back(cfg.get<string>("filterDiTau"));
  }
  cout<<"Number of single lepton triggers = "<<filterSingleLep.size()<<endl;
  cout<<"Number of X triggers (lep leg)   = "<<filterXtriggerLepLeg.size()<<endl;
  cout<<"Number of X triggers (tau leg)   = "<<filterXtriggerTauLeg.size()<<endl;
  
  // tau cuts
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut       = cfg.get<float>("dzTauCut");
  const bool  applyTauId     = cfg.get<bool>("ApplyTauId");
  if(ch=="tt"){
    ptTauLowCut    = cfg.get<float>("ptTau2LowCut");
    etaTauCut      = cfg.get<float>("etaTau2Cut");
    dzTauCut       = cfg.get<float>("dzTau2Cut");
  }


  // tau energy scale corrections
  const float shift_tes_1prong = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0 = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  // for lep->tau fakes
  const float shift_tes_lepfake_1prong = cfg.get<float>("TauEnergyScaleShift_LepFake_OneProng");
  const float shift_tes_lepfake_1p1p0 = cfg.get<float>("TauEnergyScaleShift_LepFake_OneProngOnePi0");
  const float shift_tes_lepfake_3prong = cfg.get<float>("TauEnergyScaleShift_LepFake_ThreeProng");

  // pair selection
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // extra electron veto
  const float ptVetoElectronCut  = cfg.get<float>("ptVetoElectronCut");  
  const float etaVetoElectronCut = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut = cfg.get<float>("dxyVetoElectronCut");  
  const float dzVetoElectronCut  = cfg.get<float>("dzVetoElectronCut"); 
  const bool applyVetoElectronId = cfg.get<bool>("applyVetoElectronId");
  const float isoVetoElectronCut = cfg.get<float>("isoVetoElectronCut");  
  
  // extra muon veto
  const float ptVetoMuonCut  = cfg.get<float>("ptVetoMuonCut");  
  const float etaVetoMuonCut = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut = cfg.get<float>("dxyVetoMuonCut");  
  const float dzVetoMuonCut  = cfg.get<float>("dzVetoMuonCut"); 
  const bool  applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
  const float isoVetoMuonCut = cfg.get<float>("isoVetoMuonCut");

 // lepton cuts
  const float ptLeptonLowCut   = cfg.get<float>("pt"+lep+"LowCut");
  const float ptLeptonHighCut  = cfg.get<float>("pt"+lep+"HighCut");
  const float etaLeptonCut     = cfg.get<float>("eta"+lep+"Cut");
  const float dxyLeptonCut     = cfg.get<float>("dxy"+lep+"Cut");
  const float dzLeptonCut      = cfg.get<float>("dz"+lep+"Cut");

  const bool  applyLeptonId    = cfg.get<bool>("Apply"+lep+"Id");

  //dilepton veto
  const float ptDiLeptonVeto     = cfg.get<float>("ptDi"+lep+"Veto");  
  const float etaDiLeptonVeto    = cfg.get<float>("etaDi"+lep+"Veto");
  const float dxyDiLeptonVeto    = cfg.get<float>("dxyDi"+lep+"Veto");  
  const float dzDiLeptonVeto     = cfg.get<float>("dzDi"+lep+"Veto"); 
  const bool applyDiLeptonVetoId = cfg.get<bool>("applyDi"+lep+"VetoId");
  const bool applyDiLeptonOS     = cfg.get<bool>("applyDi"+lep+"OS");
  const float isoDiLeptonVeto    = cfg.get<float>("isoDi"+lep+"Veto");
  const float drDiLeptonVeto     = cfg.get<float>("drDi"+lep+"Veto"); 

  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
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

  unsigned int lhc_run_era = 2;
  bool applyRun1topPtWeights = true;
  if (applyRun1topPtWeights) lhc_run_era =1;

  bool includeTTbarUncWeights = true;

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
  else{
    ifstream input;
    std::string infile;
    
    input.open(infiles);

    while(true){
      input>>infile;
      if(!input.eof()){
	if (infile.length() > 0){
	  fileList.push_back(infile);
	  NumberOfFiles +=1 ;
	}
      }
      else
	break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

  if(NumberOfFiles < jfile) jfile = NumberOfFiles;

  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }

  
  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");

  // PU reweighting - initialization
  PileUp * PUofficial = new PileUp();
  if(ApplyPUweight){
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(pileUpInDataFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(pileUpInMCFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(TString(pileUpforMC));
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  

  // Lepton Scale Factors
  // Lepton Id+Iso scale factor
  ScaleFactor * SF_lepIdIso = new ScaleFactor();
  ScaleFactor * SF_lepIdIso_antiiso = new ScaleFactor();
  ScaleFactor * SF_SingleLepTrigger = new ScaleFactor();
  ScaleFactor * SF_XTriggerLepLeg   = new ScaleFactor();
  ScaleFactor * SF_SingleLepTrigger_antiiso = new ScaleFactor();
  ScaleFactor * SF_XTriggerLepLeg_antiiso   = new ScaleFactor();

  if(ApplyLepSF){
    SF_lepIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(idIsoEffFile));
    SF_lepIdIso_antiiso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(idIsoEffFile_antiiso));
    SF_SingleLepTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(singleLepTrigEffFile));
    SF_XTriggerLepLeg->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(xTrigLepLegEffFile));
    SF_SingleLepTrigger_antiiso -> init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(singleLepTrigEffFile_antiiso));
    SF_XTriggerLepLeg_antiiso   -> init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(xTrigLepLegEffFile_antiiso));
  }
  // For tau leg of cross-trigger a different implementation is used
  TString filename_XTrigTauLegSF = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_sm_moriond_v1.root";
  TFile f_XTrigTauLegSF(filename_XTrigTauLegSF);
  if (f_XTrigTauLegSF.IsZombie()) {std::cout << " workspace file " << filename_XTrigTauLegSF << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w_XTrigTauLegSF = (RooWorkspace*)f_XTrigTauLegSF.Get("w");
  f_XTrigTauLegSF.Close();

  //Lepton Scale Factors from KIT - for MSSM 2016
  TString filename_kitLeptonSF = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_4.root";
  TFile f_kitLeptonSF(filename_kitLeptonSF);
  if (f_kitLeptonSF.IsZombie()) {std::cout << " workspace file " << filename_kitLeptonSF << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w_kitLeptonSF = (RooWorkspace*)f_kitLeptonSF.Get("w");
  f_kitLeptonSF.Close();


  // Workspace containing tracking efficiency weights 
  TString workspace_filename = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_3.root";
  TFile *f_workspace = new TFile(workspace_filename,"read");
  if (f_workspace->IsZombie()) {std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");
  //f.Close();

  // Zpt reweighting for LO DY samples 
  TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+ZptweightFile,"read");
  //TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+"DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016_BtoH.root","read");
  TH2D * h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

  // lepton to tau fake init
  LepTauFakeRate *leptauFR = new LepTauFakeRate();
  leptauFR->Init();

  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_" + ch + "_Sync.root";

  std::cout <<rootFileName <<std::endl;  

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5,0.5);
  
  TTree * tree = new TTree("TauCheck","TauCheck");

  Synch17Tree *otree = new Synch17Tree(tree);
  initializeCPvar(otree);

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
  TFile* inputFile_visPtResolution = new TFile(svFitPtResFile.data());

  //Systematics init
  TauScaleSys * tauScaleSys =0;
  TauOneProngScaleSys* tauOneProngScaleSys =0;
  TauOneProngOnePi0ScaleSys* tauOneProngOnePi0ScaleSys=0;
  TauThreeProngScaleSys* tauThreeProngScaleSys=0;

  //LepTauFakeScaleSys * lepTauFakeScaleSys = 0;
  LepTauFakeOneProngScaleSys* lepTauFakeOneProngScaleSys =0;
  LepTauFakeOneProngOnePi0ScaleSys* lepTauFakeOneProngOnePi0ScaleSys=0;
  LepTauFakeThreeProngScaleSys * lepTauFakeThreeProngScaleSys =0;

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
      JESUncertainties * jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Summer16_UncertaintySources_AK4PFchs.txt");
      std::vector<std::string> JESnames = jecUncertainties->getUncertNames();
      for (unsigned int i = 0; i<JESnames.size(); i++) std::cout << "i: "<< i << ", JESnames.at(i) : " << JESnames.at(i) << std::endl;
      for (unsigned int i = 0; i<JESnames.size(); i++){
        JetEnergyScaleSys * aJESobject = new JetEnergyScaleSys(otree, TString(JESnames.at(i)));
        aJESobject->SetConfig(&cfg);
        aJESobject->SetBtagScaling(&inputs_btag_scaling_medium);
        aJESobject->SetJESUncertainties(jecUncertainties);
        jetEnergyScaleSys.push_back(aJESobject);
      }	  
    }
    else { // use JEC uncertainty from analysis tree
      JetEnergyScaleSys * singleJES = new JetEnergyScaleSys(otree, TString("JES"));
      singleJES->SetConfig(&cfg);
      singleJES->SetBtagScaling(&inputs_btag_scaling_medium);
      singleJES->SetJESUncertainties(jecUncertainties);
      jetEnergyScaleSys.push_back(singleJES);
    }
    
  }
  
  // list of met filters
  std::vector<TString> met_filters_list ;
  met_filters_list.push_back("Flag_HBHENoiseFilter");
  met_filters_list.push_back("Flag_HBHENoiseIsoFilter");
  met_filters_list.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  met_filters_list.push_back("Flag_goodVertices");
  met_filters_list.push_back("Flag_eeBadScFilter");
  met_filters_list.push_back("Flag_globalTightHalo2016Filter");
  met_filters_list.push_back("Flag_BadPFMuonFilter");
  met_filters_list.push_back("Flag_BadChargedCandidateFilter");


  int counter[20];

  ///////////////FILE LOOP///////////////

  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
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

    AC1B analysisTree(_tree, isData);

	// set AC1B for JES systematics
	if (!isData && ApplySystShift && jetEnergyScaleSys.size() >0){
		for (unsigned int i=0; i<jetEnergyScaleSys.size(); i++)
			(jetEnergyScaleSys.at(i))->SetAC1B(&analysisTree);
	}
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
  ///////////////EVENT LOOP///////////////
    

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {       
      counter[0]++;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (isData)
	       nWeightedEventsH->Fill(0., 1.);
      else
	       nWeightedEventsH->Fill(0., analysisTree.genweight);

	  //Skip events not passing the MET filters, if applied
      if (ApplyMetFilters && !passedAllMetFilters(&analysisTree, met_filters_list, isData)) continue;
      counter[1]++;

      // Check if all triggers are existent in each event and save index
      vector<int> nSingleLepTrig(filterSingleLep.size(),-1);
      vector<int> nXTrigLepLeg(filterXtriggerLepLeg.size(),-1);
      vector<int> nXTrigTauLeg(filterXtriggerTauLeg.size(),-1);

      if(ApplyTrigger){

	vector<bool> checkFilterSingleLep(filterSingleLep.size(), false); 
	vector<bool> checkFilterXTrigLepLeg(filterXtriggerLepLeg.size(), false); 
	vector<bool> checkFilterXTrigTauLeg(filterXtriggerTauLeg.size(), false);
	unsigned int nfilters = analysisTree.run_hltfilters->size();

	for (unsigned int i=0; i<nfilters; ++i) {
	  TString HLTFilter(analysisTree.run_hltfilters->at(i));
	  for(unsigned int i_trig=0; i_trig<filterSingleLep.size(); i_trig++){
	    if (HLTFilter==filterSingleLep.at(i_trig)){ nSingleLepTrig.at(i_trig) = i; checkFilterSingleLep.at(i_trig) = true;}
	  }
	  for(unsigned int i_trig=0; i_trig<filterXtriggerLepLeg.size(); i_trig++){
	    if (HLTFilter==filterXtriggerLepLeg.at(i_trig)){ nXTrigLepLeg.at(i_trig) = i; checkFilterXTrigLepLeg.at(i_trig) = true; }
	  }
	  for(unsigned int i_trig=0; i_trig<filterXtriggerTauLeg.size(); i_trig++){
	    if (HLTFilter==filterXtriggerTauLeg.at(i_trig)){ nXTrigTauLeg.at(i_trig) = i; checkFilterXTrigTauLeg.at(i_trig) = true;}
	  }
	}
      }
      
      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
	  

      bool overlapEvent = true;
      for (unsigned int iEvent=0; iEvent<runList.size(); ++iEvent) {
      	if (runList.at(iEvent)==otree->run && eventList.at(iEvent)==otree->evt) {
      	  overlapEvent = false;	  
      	}
      }

      if (overlapEvent&&checkOverlap) continue;
      nonOverlap++;

      counter[2]++;

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;

      // weights
      if(ApplyPUweight) fill_weight(&analysisTree, otree, PUofficial, isData);
      
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;

      // tau selection

      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) { 
       
        if (analysisTree.tau_pt[it]<=ptTauLowCut) continue;
        if (fabs(analysisTree.tau_eta[it])>=etaTauCut) continue;

        if (fabs(fabs(analysisTree.tau_charge[it])-1)>0.001) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>=dzTauCut) continue;

        if (applyTauId && analysisTree.tau_decayModeFinding[it] < 0.5) continue;
        taus.push_back(it);
      }
	  vector<int> duplicatemuons; duplicatemuons.clear();
      counter[3]++;

      //selecting leptons 
      vector<int> leptons; leptons.clear();
      if(ch == "et"){

        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
			
          bool electronMvaId = analysisTree.electron_mva_wp80_general_Spring16_v1[ie];

          if (analysisTree.electron_pt[ie]<=ptLeptonLowCut) continue;

          if (fabs(analysisTree.electron_eta[ie])>=etaLeptonCut) continue;
          if (fabs(analysisTree.electron_dxy[ie])>=dxyLeptonCut) continue;

	  if (fabs(analysisTree.electron_dz[ie])>=dzLeptonCut) continue;
          if (!electronMvaId  &&  applyLeptonId) continue;
          if (!analysisTree.electron_pass_conversion[ie]  && applyLeptonId) continue;

          if (analysisTree.electron_nmissinginnerhits[ie]>1 && applyLeptonId) continue;

          leptons.push_back(ie);
        }

      }

      if(ch == "mt"){
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

	  bool muonMediumId = isIdentifiedMediumMuon(im, &analysisTree, isData);
  			
          if (analysisTree.muon_pt[im]<=ptLeptonLowCut) continue;
          if (fabs(analysisTree.muon_eta[im])>=etaLeptonCut) continue;

          if (fabs(analysisTree.muon_dxy[im])>=dxyLeptonCut) continue;
          if (fabs(analysisTree.muon_dz[im])>=dzLeptonCut) continue;

          if (!muonMediumId &&  applyLeptonId) continue;

          leptons.push_back(im);
        }
      }
      counter[4]++;

      if (leptons.size()==0) continue;
      if (taus.size()==0) continue;
      counter[5]++;//

      // selecting electron and tau pair (OS or SS);
      int leptonIndex = -1;
      int tauIndex = -1;
      
      float isoLepMin = 1e+10;
      float isoTauMax = -1;      
      float lep_pt_max =  -1;
      float tau_pt_max =  -1;


      bool isSingleLepTrigStored = false;
      bool isXTrigStored         = false;
      
      bool isSingleLepTrig = false;
      vector<bool> isXTrigLepLeg(filterXtriggerLepLeg.size(), false);
      vector<bool> isXTrigTauLeg(filterXtriggerTauLeg.size(), false);
      bool isXTrig         = false;
      otree->singleLepTrigger = false;
      otree->xTrigger = false;
      otree->trg_singlemuon = false;
      otree->trg_singleelectron = false;
      for (unsigned int it=0; it<taus.size(); ++it) {
	counter[6]++;

	unsigned int tIndex = taus.at(it);
	unsigned int Size = leptons.size();
	if(ch=="tt")Size = taus.size();


	for (unsigned int il=0; il<Size; ++il) {
	  unsigned int lIndex  = leptons.at(il);
	  if(ch=="tt") lIndex  = taus.at(il);
	  float relIsoLep  = -9999.;
	  if(ch=="mt")  relIsoLep = (abs_Iso_mt(lIndex, &analysisTree, dRiso) / analysisTree.muon_pt[lIndex] );
	  else if(ch=="et")  relIsoLep = (abs_Iso_et(lIndex, &analysisTree, dRiso) / analysisTree.electron_pt[lIndex] );
	  
	  float lep_pt     = -9999.;
	  float lep_pt_max = -9999.;
	  float lep_eta    = -9999.;
	  float lep_phi    = -9999.;
	  
	  if(ch=="mt"){
	    lep_pt =      analysisTree.muon_pt[lIndex]; 
	    lep_eta =     analysisTree.muon_eta[lIndex]; 
	    lep_phi =     analysisTree.muon_phi[lIndex];}
	  if(ch=="et"){         
	    lep_pt =      analysisTree.electron_pt[lIndex];
	    lep_eta =     analysisTree.electron_eta[lIndex]; 
	    lep_phi =     analysisTree.electron_phi[lIndex];}
	  if(ch=="tt"){         
	    lep_pt =      analysisTree.tau_pt[lIndex];
	    lep_eta =     analysisTree.tau_eta[lIndex]; 
	    lep_phi =     analysisTree.tau_phi[lIndex];
	    if(lIndex==tIndex)continue;
	  }
	  
	  if(ApplyTrigger){  //TO FIX
	    
	    isSingleLepTrig = false;
	    isXTrigLepLeg.assign(isXTrigLepLeg.size(), false);
	    isXTrigTauLeg.assign(isXTrigTauLeg.size(), false);
	    isXTrig         = false;
	    
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      float dRtrigLep = deltaR(lep_eta, lep_phi, analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);        
	      float dRtrigTau = deltaR(analysisTree.tau_eta[tIndex], analysisTree.tau_phi[tIndex], analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);        
	      
 	      if (dRtrigLep < deltaRTrigMatch){
		
		for(unsigned int i_trig = 0; i_trig<filterSingleLep.size(); i_trig++){
		  if(ch=="mt" && lep_pt<23) continue;
		  if (nSingleLepTrig.at(i_trig) == -1) continue;
		  if (analysisTree.trigobject_filters[iT][nSingleLepTrig.at(i_trig)]) isSingleLepTrig = true;
		  
		}
		for(unsigned int i_trig = 0; i_trig<filterXtriggerLepLeg.size(); i_trig++){
		  if(ch=="mt" && lep_pt>=23) continue;
		  if (nXTrigLepLeg.at(i_trig) == -1) continue;
		  if (analysisTree.trigobject_filters[iT][nXTrigLepLeg.at(i_trig)]) isXTrigLepLeg.at(i_trig) = true;
		}
		
	      }
	      
	      if (dRtrigTau < deltaRTrigMatch && ch=="mt"){
		
		for(unsigned int i_trig = 0; i_trig<filterXtriggerTauLeg.size(); i_trig++){
		  if (nXTrigTauLeg.at(i_trig) == -1) continue;
		  if (analysisTree.trigobject_filters[iT][nXTrigTauLeg.at(i_trig)]) isXTrigTauLeg.at(i_trig) = true;
		}
	      }
	      
	    }
	    
	    for(unsigned int i_trig = 0; i_trig<filterXtriggerTauLeg.size(); i_trig++){
	      if(isXTrigLepLeg.at(i_trig) && isXTrigTauLeg.at(i_trig)) isXTrig = true;
	    }
	    
	    if ( !(isSingleLepTrig || isXTrig) ) continue;
	    
	    if (isSingleLepTrig){
			if (ch=="mt") otree->trg_singlemuon = true;
			else if (ch=="et") otree->trg_singleelectron = true;
			
	    }
	    
	  }
        
	  
	  counter[7]++;
	  
          float absIsoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
          float relIsoTau = absIsoTau / analysisTree.tau_pt[tIndex];

          
          float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
                lep_eta, lep_phi);

          if (dR<dRleptonsCut) continue;

          // kinematic match
          if (lep_pt<=ptLeptonHighCut) continue;

          // change pair
//        cout<<iEntry<<". "<<lIndex<<", "<<tauIndex<<". "<<lep_pt<<" - "<<lep_pt_max;
          bool changePair =  false;
          if (relIsoLep<isoLepMin) {

            changePair = true;
          }
          else if (fabs(relIsoLep - isoLepMin) < 1.e-5) {
            if (lep_pt > lep_pt_max) {

              changePair = true;
            }     
            else if (fabs(lep_pt - lep_pt_max) < 1.e-5) {
              if (absIsoTau > isoTauMax) {   

                changePair = true;
              }
              else if ((absIsoTau - isoTauMax) < 1.e-5){
                if (analysisTree.tau_pt[tIndex] > tau_pt_max){

                  changePair = true;
                }
              }
            }
          }
	  counter[8]++;
	  
          if (changePair){
	    isoLepMin  = relIsoLep;
	    lep_pt_max=lep_pt;
	    tau_pt_max=analysisTree.tau_pt[tIndex];
	    leptonIndex = lIndex;
	    isoTauMax = absIsoTau;
	    tauIndex = tIndex;
	    isSingleLepTrigStored = isSingleLepTrig;
	    isXTrigStored = isXTrig;
          }
        }
      }
      
      if (leptonIndex<0) continue;
      if (tauIndex<0) continue;
      counter[9]++;

      //filling variables
      TLorentzVector leptonLV;

      // used for trigger weights
      w_XTrigTauLegSF->var("t_pt")  -> setVal(analysisTree.tau_pt[tauIndex]);
      w_XTrigTauLegSF->var("t_eta") -> setVal(analysisTree.tau_eta[tauIndex]);
      w_XTrigTauLegSF->var("t_dm")  -> setVal(analysisTree.tau_decayMode[tauIndex]);
      double sf_trig_t = 0;
      
      if(ch=="mt") {
      	FillMuTau(&analysisTree, otree, leptonIndex, dRiso);
      	
        leptonLV.SetXYZM(analysisTree.muon_px[leptonIndex],
			 analysisTree.muon_py[leptonIndex],
			 analysisTree.muon_pz[leptonIndex],
			 muonMass);
	
	// used for trigger weights: channel dependent settings
	sf_trig_t  = w_XTrigTauLegSF -> function("t_genuine_TightIso_mt_ratio")->getVal();


	
	// tracking efficiency weight	
        if (!isData && ApplyLepSF) {
	  w->var("m_eta")->setVal(analysisTree.muon_eta[leptonIndex]); 
	  otree->trkeffweight = (double)( w->function("m_trk_ratio")->getVal());
        }
      } 
      else if(ch=="et"){
	FillETau(&analysisTree, otree, leptonIndex, dRiso);
	
        leptonLV.SetXYZM(analysisTree.electron_px[leptonIndex],
			 analysisTree.electron_py[leptonIndex],
			 analysisTree.electron_pz[leptonIndex],
			 electronMass);
	
	// used for trigger weights: channel dependent settings
	sf_trig_t  = w_XTrigTauLegSF -> function("t_genuine_TightIso_et_ratio")->getVal();

	// tracking efficiency weight
        if (!isData && ApplyLepSF) {
	  w->var("e_eta")->setVal(analysisTree.electron_eta[leptonIndex]); 
	  w->var("e_pt")->setVal(analysisTree.electron_pt[leptonIndex]); 	
	  otree->trkeffweight = (double)( w->function("e_trk_ratio")->getVal());
	}
      }
      else if(ch=="tt"){
	FillTau(&analysisTree, otree, leptonIndex);
	
        leptonLV.SetXYZM(analysisTree.tau_px[leptonIndex],
			 analysisTree.tau_py[leptonIndex],
			 analysisTree.tau_pz[leptonIndex],
			 tauMass);
	/*TO FIX
	// used for trigger weights: channel dependent settings
	sf_trig_t  = w_XTrigTauLegSF -> function("t_genuine_TightIso_et_ratio")->getVal();

	// tracking efficiency weight
        if (!isData && ApplyLepSF) {
	  w->var("e_eta")->setVal(analysisTree.tau_eta[leptonIndex]); 
	  w->var("e_pt")->setVal(analysisTree.tau_pt[leptonIndex]); 	
	  otree->trkeffweight = (double)( w->function("e_trk_ratio")->getVal());
	  }*/
      }
      
      if (!isData && ApplyLepSF) {//TO FIX
	//std::cout<< iEntry <<std::endl;
	otree->idisoweight_1 = SF_lepIdIso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	//std::cout<< "SF_lepIdIso" <<std::endl;

	// calculation of trigger weights
	double scalefactor = 1;
	if(isSingleLepTrigStored && !isXTrigStored)  //std::cout<< "isSingleLepTrigStored && !isXTrigStored:" << std::endl , std::cout<< "  leptonLV.Pt():  " << leptonLV.Pt() << "   leptonLV.Eta(): " << leptonLV.Eta() << std::endl ,
	  scalefactor = SF_SingleLepTrigger -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	else if(isXTrigStored && !isSingleLepTrigStored) scalefactor = SF_XTriggerLepLeg   -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta())*sf_trig_t;
	else if(isXTrigStored && isSingleLepTrigStored){
	  cout<<"Both triggers returned true. This should not happen. Please check. Exiting"<<endl;
	  exit(-1);
	}  
	otree->trigweight_1 = scalefactor;
	otree->singleLepTrigger = isSingleLepTrigStored;
	otree->xTrigger = isXTrigStored;

	if (otree->iso_1>=0.15 && otree->iso_1<=0.3){  
	  otree->idisoweight_antiiso_1 = SF_lepIdIso_antiiso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	  //std::cout<< "SF_lepIdIso_antiiso" << std::endl;
	  scalefactor = 1;
	  if(isSingleLepTrigStored && !isXTrigStored)      scalefactor = SF_SingleLepTrigger_antiiso -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	  else if(isXTrigStored && !isSingleLepTrigStored) scalefactor = SF_XTriggerLepLeg_antiiso   -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta())*sf_trig_t;
	  else if(isXTrigStored && isSingleLepTrigStored){
	    cout<<"Both triggers returned true. This should not happen. Please check. Exiting"<<endl;
	    exit(-1);
	  }
	  otree->trigweight_antiiso_1 = scalefactor;
	}

	// use lepton SF from KIT workspace - for MSSM 
	// already binned in lepton isolation, no need for antiiso weights in separate branches.
	if (ApplyLeptonSFfromKIT){
	  if (ch == "mt"){
	  w_kitLeptonSF ->var("m_eta")->setVal(leptonLV.Eta());
	  w_kitLeptonSF ->var("m_pt")->setVal(leptonLV.Pt());
	  w_kitLeptonSF ->var("m_iso")->setVal(otree->iso_1);
	  otree->idisoweight_1 = (w_kitLeptonSF->function("m_id_ratio")->getVal())*(w_kitLeptonSF->function("m_iso_binned_ratio")->getVal());
	  otree->trigweight_1 = w_kitLeptonSF->function("m_trgOR4_binned_ratio")->getVal();
	  }
	  else if (ch == "et"){
	  w_kitLeptonSF ->var("e_eta")->setVal(analysisTree.electron_superclusterEta[leptonIndex]);//leptonLV.Eta());
	  w_kitLeptonSF ->var("e_pt")->setVal(leptonLV.Pt());
	  w_kitLeptonSF ->var("e_iso")->setVal(otree->iso_1);
	  otree->idisoweight_1 = (w_kitLeptonSF->function("e_id_ratio")->getVal())*(w_kitLeptonSF->function("e_iso_binned_ratio")->getVal());
	  otree->trigweight_1 = w_kitLeptonSF->function("e_trg_binned_ratio")->getVal();
	  }
	}

      }
      
      counter[10]++;

      if(ch=="mt"||ch=="et") FillTau(&analysisTree, otree, tauIndex);
      else if(ch=="tt") FillTau_leading(&analysisTree, otree, tauIndex);
	  
      // tauID weight
      if (!isData && analysisTree.tau_genmatch[tauIndex]==5) otree->idisoweight_2 = 1.; 

      otree->effweight = (otree->trigweight_1)*(otree->idisoweight_1)*(otree->trigweight_2)*(otree->idisoweight_2);

      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
      //MET
      fillMET(ch, leptonIndex, tauIndex, &analysisTree, otree);

      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);

      // Zpt weight
      otree->zptweight = 1.;
      if (!isData && ((isDY && isMG ) || isEWKZ) ) {
        genV = genTools::genV(analysisTree); // gen Z boson ?
        otree->zptweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(genV.M()),h_zptweight->GetYaxis()->FindBin(genV.Pt()));
      }

      // topPt weight
	  otree->topptweight =1.;
      if(!isData){
		otree->topptweight = genTools::topPtWeight(analysisTree, lhc_run_era);
      }

	  // weights for ttbar samples to estimate uncertianties
      //fillTTbarUncWeights(&analysisTree, otree, isData, includeTTbarUncWeights);

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
      if (!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal) ){
				genV = genTools::genV(analysisTree);
				genL = genTools::genL(analysisTree);
				if(isWJets) otree->njetshad += 1;
      }

      // PF MET
      genTools::RecoilCorrections( *recoilPFMetCorrector, (!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal)) * genTools::MeanResolution,
			                     otree->met, otree->metphi,
			                     genV.Px(), genV.Py(),
			                     genL.Px(), genL.Py(),
			                     otree->njetshad,
			                     otree->met_rcmr, otree->metphi_rcmr
			                     );

      // overwriting with recoil-corrected values 
      otree->met = otree->met_rcmr;
      otree->metphi = otree->metphi_rcmr;   


      //ditau sytem
      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
					     analysisTree.tau_py[tauIndex],
					     analysisTree.tau_pz[tauIndex],
					     analysisTree.tau_mass[tauIndex]);


      // using PF MET
      TLorentzVector metLV; metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi),
					  otree->met*TMath::Sin(otree->metphi),
					  0,
					  TMath::Sqrt( otree->met*TMath::Sin(otree->metphi)*otree->met*TMath::Sin(otree->metphi) +
				          otree->met*TMath::Cos(otree->metphi)*otree->met*TMath::Cos(otree->metphi)));

	  if (ch == "mt" && checkDuplicateMuons && duplicatemuons.size()!=0 && isData){
		// add to the MET LV the duplicate muons pt
		TLorentzVector duplicateMuonsLV;
		float duplicatemu_px = 0; 
		float duplicatemu_py = 0; 
		for (unsigned int i=0; i<duplicatemuons.size(); i++){
			duplicatemu_px += analysisTree.muon_px[duplicatemuons[i]];
			duplicatemu_py += analysisTree.muon_py[duplicatemuons[i]];
			}
	    // add duplicate muons pT to MET
		float tot_px = metLV.Px()+ duplicatemu_px;
		float tot_py = metLV.Py() + duplicatemu_py;
		metLV.SetPx(tot_px);
		metLV.SetPy(tot_py);
		// save to otree
		otree->met= metLV.Pt();
		otree->metphi = metLV.Phi();
	  }
	  counter[12]++;

	  // shift the tau energy scale by decay mode and propagate to the met. 
	  if (!isData) {
	    bool isOneProng = false;
	    float shift_tes = 0.0;
		if (otree->gen_match_2 >=5){
	    	if (otree->tau_decay_mode_2 == 0){      shift_tes=shift_tes_1prong; isOneProng=true; }
	    	else if (otree->tau_decay_mode_2 ==1){  shift_tes=shift_tes_1p1p0; }
	    	else if (otree->tau_decay_mode_2 ==10){ shift_tes=shift_tes_3prong; }
		}

		else if (otree->gen_match_2 <5) {
	    	if (otree->tau_decay_mode_2 == 0){      shift_tes=shift_tes_lepfake_1prong; isOneProng=true; }
	    	else if (otree->tau_decay_mode_2 ==1){  shift_tes=shift_tes_lepfake_1p1p0; }
	    	else if (otree->tau_decay_mode_2 ==10){ shift_tes=shift_tes_lepfake_3prong; }
		}

	    correctTauES(tauLV, metLV, shift_tes, isOneProng);	    
	    otree->pt_2 = tauLV.Pt();
	    otree->m_2 = tauLV.M();
	    otree->met = metLV.Pt();
	    otree->metphi = metLV.Phi();
	  }

      TLorentzVector dileptonLV = leptonLV + tauLV;

      // visible mass
      otree->m_vis = dileptonLV.M();

	  // ditau pt
      otree->pt_tt = (dileptonLV+metLV).Pt();   

	  // mt TOT
      float mtTOT = 2*(otree->pt_1)*metLV.Pt()*(1-cos(DeltaPhi(otree->phi_1,otree->metphi)));
      mtTOT += 2*(otree->pt_2)*metLV.Pt()*(1-cos(DeltaPhi(otree->phi_2,otree->metphi))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(otree->phi_1,otree->phi_2))); 
	  otree->mt_tot = TMath::Sqrt(mtTOT);

      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;

      // dilepton veto
      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
      if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree);

	  //extra letpn veto
      otree->extraelec_veto = extra_electron_veto(leptonIndex, ch, &cfg, &analysisTree);
      otree->extramuon_veto = extra_muon_veto(leptonIndex, ch, &cfg, &analysisTree, isData);

      counter[13]++;


      otree->mt_1=mT(leptonLV,metLV);
      otree->mt_2=mT(tauLV,metLV);


      // bisector of lepton and tau transverse momenta

      float leptonUnitX = leptonLV.Px()/leptonLV.Pt();
      float leptonUnitY = leptonLV.Py()/leptonLV.Pt();

      float tauUnitX = tauLV.Px()/tauLV.Pt();
      float tauUnitY = tauLV.Py()/tauLV.Pt();

      float zetaX = leptonUnitX + tauUnitX;
      float zetaY = leptonUnitY + tauUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorVisX = leptonLV.Px() + tauLV.Px();
      float vectorVisY = leptonLV.Py() + tauLV.Py();

      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;
      otree->pzetamiss = otree->met*TMath::Cos(otree->metphi)*zetaX + otree->met*TMath::Sin(otree->metphi)*zetaY;
      counter[14]++;

      // svfit variables
      otree->m_sv   = -9999;
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
      
      if (ch=="mt") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
			      otree->againstElectronVLooseMVA6_2>0.5 && otree->againstMuonTight3_2>0.5  &&
			      otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);
      if (ch=="et") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
                            otree->againstMuonLoose3_2>0.5 && otree->againstElectronTightMVA6_2>0.5 && 
                            otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);

      if(otree->iso_1<0.35 && 
	 otree->againstMuonLoose3_2>0.5 && otree->againstElectronTightMVA6_2>0.5 && 
	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[15]++;
      if(otree->iso_1<0.35 && otree->byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
	 otree->againstElectronTightMVA6_2>0.5 && 
	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[16]++;
      if(otree->iso_1<0.35 && otree->byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
	 otree->againstMuonLoose3_2>0.5 && 
	 otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0) counter[17]++;
      if(otree->iso_1<0.35 && otree->byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
	 otree->againstMuonLoose3_2>0.5 && otree->againstElectronTightMVA6_2>0.5 ) counter[18]++;

      if (!Synch && !passedBaselineSel) continue;

      if (ApplySVFit && otree->njetspt20>0) svfit_variables(ch, &analysisTree, otree, &cfg, inputFile_visPtResolution);

      otree->Fill();

	  // evaluate systematics for MC 
      if(!isData && ApplySystShift){
	zPtWeightSys->Eval(); 
	topPtWeightSys->Eval();
	for(unsigned int i=0; i<jetEnergyScaleSys.size(); i++)
	  (jetEnergyScaleSys.at(i)) ->Eval(); 
	if (ch=="mt") {
           tauScaleSys->Eval(utils::MUTAU);
	}
	 else if (ch=="et") {
           tauScaleSys->Eval(utils::ETAU);
	 }
      }
      counter[19]++;

      selEvents++;
      //std::cout << "*************SelEv " << selEvents << "************" << std::endl;

    } // end of file processing (loop over events in one file)

    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
   }
  std::cout << "COUNTERS" << std::endl;
  for(int iC=0;iC<20;iC++) std::cout << "Counter " << iC << ":    " << counter[iC] << std::endl;

  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();

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

  if(jetEnergyScaleSys.size() >0 ){
    for (unsigned int i=0; i<jetEnergyScaleSys.size(); i++){
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

//fill the otree with the muon variables in channel mutau
void FillMuTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso){
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

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -9999;
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->againstMuonLoose3_1 = -9999;
  otree->againstMuonTight3_1 = -9999;
  otree->againstElectronVLooseMVA6_1 = -9999;
  otree->againstElectronTightMVA6_1 = -9999;

  otree->byVLooseIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byLooseIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byMediumIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byVTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byVVTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree-> byIsolationMVArun2v1DBoldDMwLTraw_1 = -9999;
  otree->chargedIsoPtSum_1 = -9999;
  otree->neutralIsoPtSum_1 = -9999;
  otree->puCorrPtSum_1 = -9999;

}

//fill the otree with the electron variables in channel etau
void FillETau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso){
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

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -9999;
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->againstMuonLoose3_1 = -9999;
  otree->againstMuonTight3_1 = -9999;
  otree->againstElectronVLooseMVA6_1 = -9999;
  otree->againstElectronTightMVA6_1 = -9999;

  otree->byVLooseIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byLooseIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byMediumIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byVTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree->byVVTightIsolationMVArun2v1DBoldDMwLT_1 = -9999;
  otree-> byIsolationMVArun2v1DBoldDMwLTraw_1 = -9999;
  otree->chargedIsoPtSum_1 = -9999;
  otree->neutralIsoPtSum_1 = -9999;
  otree->puCorrPtSum_1 = -9999;
  

}

//fill the otree with the tau variables 
void FillTau_leading(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex){
  otree->pt_1 = analysisTree->tau_pt[tauIndex];
  otree->eta_1 = analysisTree->tau_eta[tauIndex];
  otree->phi_1 = analysisTree->tau_phi[tauIndex];
  otree->q_1 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_1 = analysisTree->tau_genmatch[tauIndex];
  otree->mva_1 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->d0_1 = analysisTree->tau_leadchargedhadrcand_dxy[tauIndex];
  otree->dZ_1 = analysisTree->tau_leadchargedhadrcand_dz[tauIndex];      
  otree->iso_1 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->m_1 = analysisTree->tau_mass[tauIndex];
  otree->tau_decay_mode_1 = analysisTree->tau_decayMode[tauIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstMuonLoose3_1 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_1 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronVLooseMVA6_1 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronTightMVA6_1 = analysisTree->tau_againstElectronTightMVA6[tauIndex];

  otree->byVLooseIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byVLooseIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byLooseIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byLooseIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byMediumIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byMediumIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byTightIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byVTightIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byVTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byVVTightIsolationMVArun2v1DBoldDMwLT_1 = analysisTree->tau_byVVTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree-> byIsolationMVArun2v1DBoldDMwLTraw_1 = analysisTree->tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex];
  otree->chargedIsoPtSum_1 = analysisTree->tau_chargedIsoPtSum[tauIndex];
  otree->neutralIsoPtSum_1 = analysisTree->tau_neutralIsoPtSum[tauIndex];
  otree->puCorrPtSum_1 = analysisTree->tau_puCorrPtSum[tauIndex];
  //otree->isolationGammaCands_size_1 = analysisTree->tau_isolationGammaCands_size[tauIndex];
  //otree->signalGammaCands_size_1 = analysisTree->tau_signalGammaCands_size[tauIndex];
}

//fill the otree with the tau variables 
void FillTau(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->eta_2 = analysisTree->tau_eta[tauIndex];
  otree->phi_2 = analysisTree->tau_phi[tauIndex];
  otree->q_2 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_2 = analysisTree->tau_genmatch[tauIndex];
  otree->mva_2 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->d0_2 = analysisTree->tau_leadchargedhadrcand_dxy[tauIndex];
  otree->dZ_2 = analysisTree->tau_leadchargedhadrcand_dz[tauIndex];      
  otree->iso_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->m_2 = analysisTree->tau_mass[tauIndex];
  otree->tau_decay_mode_2 = analysisTree->tau_decayMode[tauIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstMuonLoose3_2 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_2 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronVLooseMVA6_2 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronTightMVA6_2 = analysisTree->tau_againstElectronTightMVA6[tauIndex];

  otree->byVLooseIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byVLooseIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byLooseIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byLooseIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byMediumIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byMediumIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byVTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byVTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree->byVVTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byVVTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
  otree-> byIsolationMVArun2v1DBoldDMwLTraw_2 = analysisTree->tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex];
  otree->chargedIsoPtSum_2 = analysisTree->tau_chargedIsoPtSum[tauIndex];
  otree->neutralIsoPtSum_2 = analysisTree->tau_neutralIsoPtSum[tauIndex];
  otree->puCorrPtSum_2 = analysisTree->tau_puCorrPtSum[tauIndex];
  //otree->isolationGammaCands_size_2 = analysisTree->tau_isolationGammaCands_size[tauIndex];
  //otree->signalGammaCands_size_2 = analysisTree->tau_signalGammaCands_size[tauIndex];
}

void initializeCPvar(Synch17Tree *otree){
  otree->acotautau=-9999;

  otree->tau1DecayPlaneX=-9999;
  otree->tau1DecayPlaneY=-9999;
  otree->tau1DecayPlaneZ=-9999;
  otree->tau2DecayPlaneX=-9999;
  otree->tau2DecayPlaneY=-9999;
  otree->tau2DecayPlaneZ=-9999;
}

