
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

#include "DesyTauAnalyses/NTupleMaker/interface/Spring15Tree.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"

#define pi 	3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 		   0.105658
#define tauMass 		   1.77682
#define pionMass 		   0.1396

typedef std::vector<std::pair<int,int> > lumi_json;

struct compare_lumi { //accepts two pairs, return 1 if left.first < right.first or left.first = right.first e left.second < right.second
  bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
    if (left.first < right.first)
      return 1;
    else if (left.first > right.first)
      return 0;
    else
      return left.second < right.second;
  }
};

int read_json(std::string filename, lumi_json& json);
bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json);
bool isGoodLumi(int run, int lumi, const lumi_json& json);
void fill_weight(const AC1B * analysisTree, Spring15Tree *otree, PileUp *PUofficial, bool isData);
float abs_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso);
float rel_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso);
void FillMuTau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso);
void FillETau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso);
void FillTau(const AC1B * analysisTree, Spring15Tree *otree, int tauIndex);
bool dilepton_veto_mt(const Config *cfg, const AC1B *analysisTree);
bool dilepton_veto_et(const Config *cfg, const AC1B *analysisTree);
bool extra_electron_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree);
bool extra_muon_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree, bool isData);
void fillMET(TString ch, int leptonIndex, int tauIndex, const AC1B * analysisTree, Spring15Tree *otree);
void mt_calculation(Spring15Tree *otree);
void svfit_variables(TString ch, const AC1B *analysisTree, Spring15Tree *otree, const Config *cfg, TFile *inputFile_visPtResolution);
bool isICHEPmed(int Index, const AC1B * analysisTree);
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData);
void correctTauES(TLorentzVector& Tau, TLorentzVector& Met, float relative_shift, bool tau_is_one_prong);
bool passedSummer16VetoId(const AC1B * analysisTree, int index);
bool SafeRatio(double denominator);

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
    else if (ch=="et") lep = "Electron";

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");

  bool checkDuplicateMuons = false;
  if (isData && ch == "mt") checkDuplicateMuons = cfg.get<bool>("checkDuplicateMuons");

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");

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

  //b-tag scale factors
  TString pathToBtagScaleFactors = (TString) cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv";
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
  TString pathToTaggingEfficiencies = (TString) cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_ichep2016.root";
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
  const bool isVBForGGHiggs = (infiles.find("VBFHToTauTau")== infiles.rfind("/")+1) || (infiles.find("GluGluHToTauTau")== infiles.rfind("/")+1);
  const bool isEWKZ =  infiles.find("EWKZ") == infiles.rfind("/")+1;
  const bool isMG = infiles.find("madgraph") != string::npos;
  //const bool applyRecoilCorrections = isDY || isWJets;
  

  RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
  //RecoilCorrector* recoilPuppiMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPuppiMetCorrector));
  RecoilCorrector* recoilMvaMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilMvaMetCorrector));

  if(!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs) ){
    TString RecoilDir("HTT-utilities/RecoilCorrections/data/");

    TString RecoilFileName = RecoilDir; RecoilFileName += "TypeI-PFMet_Run2016BtoH.root";
    std::cout<<RecoilFileName<<std::endl;
    recoilPFMetCorrector = new RecoilCorrector( RecoilFileName);
        
    //RecoilFileName = RecoilDir; RecoilFileName += "recoilPuppiMet.root";
    //std::cout<<RecoilFileName<<std::endl;
    //recoilPuppiMetCorrector = new RecoilCorrector( RecoilFileName);

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
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep2"));
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep3"));
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep4"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg"));
  }
  else if(ch == "et"){
    filterSingleLep.push_back(cfg.get<string>("filterSingleLep1"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg1"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg1"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg2"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg2"));
    filterXtriggerLepLeg.push_back(cfg.get<string>("filterXtriggerLepLeg3"));
    filterXtriggerTauLeg.push_back(cfg.get<string>("filterXtriggerTauLeg3"));
  }
  cout<<"Number of single lepton triggers = "<<filterSingleLep.size()<<endl;
  cout<<"Number of X triggers (lep leg)   = "<<filterXtriggerLepLeg.size()<<endl;
  cout<<"Number of X triggers (tau leg)   = "<<filterXtriggerTauLeg.size()<<endl;
  
  // tau cuts
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut        = cfg.get<float>("dzTauCut");
  const bool applyTauId      = cfg.get<bool>("ApplyTauId");

  // tau energy scale corrections
  const float shift_tes_1prong = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0 = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong = cfg.get<float>("TauEnergyScaleShift_ThreeProng");

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
  const bool applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
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
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
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

  // Workspace containing tracking efficiency weights 
  TString workspace_filename = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_3.root";
  TFile *f_workspace = new TFile(workspace_filename,"read");
  if (f_workspace->IsZombie()) {std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl; exit(1);}
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");
  //f.Close();

  // Zpt reweighting for LO DY samples
  TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+"DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016_BtoH.root","read");
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

  Spring15Tree *otree = new Spring15Tree(tree);

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
  TauOneProngScaleSys* tauOneProngScaleSys =0;
  TauOneProngOnePi0ScaleSys* tauOneProngOnePi0ScaleSys=0;
  TauThreeProngScaleSys* tauThreeProngScaleSys=0;
  ZPtWeightSys* zPtWeightSys = 0;
  TopPtWeightSys* topPtWeightSys = 0;
  std::vector<JetEnergyScaleSys*> jetEnergyScaleSys;
  JESUncertainties * jecUncertainties = 0;
  LepTauFakeScaleSys * lepTauFakeScaleSys = 0;
  if(!isData && ApplySystShift){
    tauOneProngScaleSys = new TauOneProngScaleSys(otree);
    tauOneProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauOneProngOnePi0ScaleSys = new TauOneProngOnePi0ScaleSys(otree);
    tauOneProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngScaleSys = new TauThreeProngScaleSys(otree);
    tauThreeProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    zPtWeightSys = new ZPtWeightSys(otree);
    topPtWeightSys = new TopPtWeightSys(otree);
    lepTauFakeScaleSys = new LepTauFakeScaleSys(otree);
    lepTauFakeScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
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
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {       
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (isData)
	       nWeightedEventsH->Fill(0., 1.);
      else
	       nWeightedEventsH->Fill(0., analysisTree.genweight);

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
	// Check if filters are available
	/*
	for(unsigned int i_trig=0; i_trig<filterSingleLep.size(); i_trig++){
	  if (checkFilterSingleLep.at(i_trig) == false && filterSingleLep.at(i_trig) != "" ){ cout<<"Following filter is not available : "<<filterSingleLep.at(i_trig)<<". Please check!"<<endl; exit( -1 );}
	}
	for(unsigned int i_trig=0; i_trig<filterXtriggerLepLeg.size(); i_trig++){
	  if (checkFilterXTrigLepLeg.at(i_trig) == false && filterXtriggerLepLeg.at(i_trig) != ""){ cout<<"Following filter is not available : "<<filterXtriggerLepLeg.at(i_trig)<<". Please check!"<<endl; exit( -1 ); }
	}
	for(unsigned int i_trig=0; i_trig<filterXtriggerLepLeg.size(); i_trig++){
	  if (checkFilterXTrigTauLeg.at(i_trig) == false && filterXtriggerTauLeg.at(i_trig) != "" ){ cout<<"Following filter is not available : "<<filterXtriggerTauLeg.at(i_trig)<<". Please check!"<<endl; exit( -1 ); }
	}
	*/
      }
      
      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
	  
	  //if (otree->evt != 8993558 && otree->evt != 8563251 && otree->evt != 8765306 && otree->evt != 8776865 && otree->evt != 8816054 && otree->evt !=8763551)	
	  //continue;

      bool overlapEvent = true;
      for (unsigned int iEvent=0; iEvent<runList.size(); ++iEvent) {
      	if (runList.at(iEvent)==otree->run && eventList.at(iEvent)==otree->evt) {
      	  overlapEvent = false;	  
      	}
      }

      if (overlapEvent&&checkOverlap) continue;
      nonOverlap++;

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
        if (applyTauId &&
            //analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) 
            analysisTree.tau_decayModeFinding[it] < 0.5) continue;
        
        taus.push_back(it);
      }

	  otree->n_badmuons=0;
	  otree->n_duplicatemuons =0;
	  vector<int> duplicatemuons; duplicatemuons.clear();

      //selecting leptons 
      vector<int> leptons; leptons.clear();
      if(ch == "et"){
        
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
			
          bool electronMvaId = analysisTree.electron_mva_wp80_general_Spring16_v1[ie];//analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie];
    
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

		  // check for duplicates and bad muons
         // check for duplicate muons, save their IDs in a vector and add to the Met later.

		  if (analysisTree.muon_isDuplicate[im]) {otree->n_duplicatemuons +=1; duplicatemuons.push_back(im);}
		  if (analysisTree.muon_isBad[im]) otree->n_badmuons+=1;
//          bool muonMediumId = analysisTree.muon_isMedium[im];
          bool muonMediumId = isIdentifiedMediumMuon(im, &analysisTree, isData);
  			
          if (analysisTree.muon_pt[im]<=ptLeptonLowCut) continue;
          if (fabs(analysisTree.muon_eta[im])>=etaLeptonCut) continue;
          if (fabs(analysisTree.muon_dxy[im])>=dxyLeptonCut) continue;
          if (fabs(analysisTree.muon_dz[im])>=dzLeptonCut) continue;
          if (!muonMediumId &&  applyLeptonId) continue;

          leptons.push_back(im);
        }
      }

      if (leptons.size()==0) continue;
      if (taus.size()==0) continue;


      // selecting electron and tau pair (OS or SS);
      int leptonIndex = -1;
      int tauIndex = -1;
      
      float isoLepMin = 1e+10;
      float isoTauMin = 1e+10;      

      bool isSingleLepTrigStored = false;
      bool isXTrigStored         = false;
      
      bool isSingleLepTrig = false;
      vector<bool> isXTrigLepLeg(filterXtriggerLepLeg.size(), false);
      vector<bool> isXTrigTauLeg(filterXtriggerTauLeg.size(), false);
      bool isXTrig         = false;
      otree->singleLepTrigger = false;
      otree->xTrigger = false;

      for (unsigned int il=0; il<leptons.size(); ++il) {
        unsigned int lIndex  = leptons.at(il);

	float relIsoLep = rel_Iso(lIndex, ch, &analysisTree, dRiso);

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

	for (unsigned int it=0; it<taus.size(); ++it) {
	  unsigned int tIndex = taus.at(it);

	  if(ApplyTrigger){  

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
	  }
        


          float absIsoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
          float relIsoTau = absIsoTau / analysisTree.tau_pt[tIndex];

          if(ch=="mt"){lep_pt_max =  analysisTree.muon_pt[leptonIndex];}
          if(ch=="et"){lep_pt_max =  analysisTree.electron_pt[leptonIndex];}

          float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
                lep_eta, lep_phi);

          if (dR<dRleptonsCut) continue;

          // kinematic match
          if (lep_pt<=ptLeptonHighCut) continue;
          
          // change pair
//        cout<<iEntry<<". "<<lIndex<<", "<<tauIndex<<". "<<lep_pt<<" - "<<lep_pt_max;
          bool changePair =  false;
          if (relIsoLep<isoLepMin) {
//          cout<<"ele iso - ";

            changePair = true;
          }
          else if (fabs(relIsoLep - isoLepMin) < 1.e-5) {
            if (lep_pt > lep_pt_max) {
//            cout<<"ele pt - "; 

              changePair = true;
            }     
            else if (fabs(lep_pt - lep_pt_max) < 1.e-5) {
              if (absIsoTau > isoTauMin) {    // -->>> changed here, before was "<" 
//                cout<<"tau iso - "; 

                changePair = true;
              }
              else if ((absIsoTau - isoTauMin) < 1.e-5){
                if (analysisTree.tau_pt[tIndex] > analysisTree.tau_pt[tauIndex]) {
//                cout<<"tau pt"; 

                  changePair = true;
                }
              }
            }
          }
//        cout<<endl;

          if (changePair){
	    isoLepMin  = relIsoLep;
	    leptonIndex = lIndex;
	    isoTauMin = absIsoTau;
	    tauIndex = tIndex;
	    isSingleLepTrigStored = isSingleLepTrig;
	    isXTrigStored = isXTrig;
          }
        }
      }
      
      if (leptonIndex<0) continue;
      if (tauIndex<0) continue;

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
	  otree->trkeffweight_1 = (double)( w->function("m_trk_ratio")->getVal());
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
	  otree->trkeffweight_1 = (double)( w->function("e_trk_ratio")->getVal());
	}
      }
      
      if (!isData && ApplyLepSF) {
	
	otree->idisoweight_1 = SF_lepIdIso->get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	
	// calculation of trigger weights
	double scalefactor = 1;
	if(isSingleLepTrigStored && !isXTrigStored)      scalefactor = SF_SingleLepTrigger -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
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
	  
	  scalefactor = 1;
	  if(isSingleLepTrigStored && !isXTrigStored)      scalefactor = SF_SingleLepTrigger_antiiso -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta());
	  else if(isXTrigStored && !isSingleLepTrigStored) scalefactor = SF_XTriggerLepLeg_antiiso   -> get_ScaleFactor(leptonLV.Pt(), leptonLV.Eta())*sf_trig_t;
	  else if(isXTrigStored && isSingleLepTrigStored){
	    cout<<"Both triggers returned true. This should not happen. Please check. Exiting"<<endl;
	    exit(-1);
	  }
	  otree->trigweight_antiiso_1 = scalefactor;
	}
      }
      
      
      FillTau(&analysisTree, otree, tauIndex);
	  
      // tauID weight
      if (!isData && analysisTree.tau_genmatch[tauIndex]==5) otree->idisoweight_2 = 0.95; 

	  // taking tau ID weight from the workspace
	  //if (!isData) {
	  //w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex]); //tau pt
	  //w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex]); // tau eta
	  //w->var("t_dm")->setVal(analysisTree.tau_genmatch[tauIndex]); // tau decay mode
	  //double tauIDweight =  w->function("t_iso_mva_m_pt30_sf")->getVal();
	  //otree->idisoweight_2 = (double)( w->function("t_iso_mva_m_pt30_sf")->getVal());
	  //}

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
      if(!isData)
	    otree->topptweight = genTools::topPtWeight(analysisTree);


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
      if (!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs) ){
				genV = genTools::genV(analysisTree);
				genL = genTools::genL(analysisTree);
				if(isWJets) otree->njetshad += 1;
      }

      // MVA MET      
      genTools::RecoilCorrections( *recoilMvaMetCorrector, (!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs)) * genTools::MeanResolution,
			                     otree->mvamet, otree->mvametphi,
			                     genV.Px(), genV.Py(),
			                     genL.Px(), genL.Py(),
			                     otree->njetshad,
			                     otree->mvamet_rcmr, otree->mvametphi_rcmr
			                     );

      // overwriting with recoil-corrected values 
      otree->mvamet = otree->mvamet_rcmr;
      otree->mvametphi = otree->mvametphi_rcmr;            
      //otree->mt_rcmr_1 = calc::mt(otree->pt_1, otree->phi_1, otree->mvamet_rcmr, otree->mvametphi_rcmr);
      //otree->mt_rcmr_2 = calc::mt(otree->pt_2, otree->phi_2, otree->mvamet_rcmr, otree->mvametphi_rcmr);     
      //otree->pzetamiss_rcmr = calc::pzetamiss( zetaX, zetaY, otree->mvamet_rcmr, otree->mvametphi_rcmr);

      // PF MET
      genTools::RecoilCorrections( *recoilPFMetCorrector, (!isData && applyRecoilCorrections && (isDY || isWJets || isVBForGGHiggs)) * genTools::MeanResolution,
			                     otree->met, otree->metphi,
			                     genV.Px(), genV.Py(),
			                     genL.Px(), genL.Py(),
			                     otree->njetshad,
			                     otree->met_rcmr, otree->metphi_rcmr
			                     );

      // overwriting with recoil-corrected values 
      otree->met = otree->met_rcmr;
      otree->metphi = otree->metphi_rcmr;   

      // save genV and genL
      otree->genV_px = genV.Px();
      otree->genV_py = genV.Py();
      otree->genV_pz = genV.Pz();
      otree->genV_e  = genV.E();
      otree->genL_px = genL.Px();
      otree->genL_py = genL.Py();
      otree->genL_pz = genL.Pz();
      otree->genL_e  = genL.E();

      //end MET Recoil Corrections

      //ditau sytem
      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
					     analysisTree.tau_py[tauIndex],
					     analysisTree.tau_pz[tauIndex],
					     analysisTree.tau_mass[tauIndex]);

	  // using MVA MET
      /*TLorentzVector metLV; metLV.SetXYZT(otree->mvamet*TMath::Cos(otree->mvametphi),
                						 otree->mvamet*TMath::Sin(otree->mvametphi),
               							 0,
                						 TMath::Sqrt( otree->mvamet*TMath::Sin(otree->mvametphi)*otree->mvamet*TMath::Sin(otree->mvametphi) +
                         				 otree->mvamet*TMath::Cos(otree->mvametphi)*otree->mvamet*TMath::Cos(otree->mvametphi)));*/

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

	  // shift the tau energy scale by decay mode and propagate to the met. 
	  if (!isData) {
	    bool isOneProng = false;
	    float shift_tes = 0.0;
	    if (otree->tau_decay_mode_2 == 0){      shift_tes=shift_tes_1prong; isOneProng=true; }
	    else if (otree->tau_decay_mode_2 ==1){  shift_tes=shift_tes_1p1p0; }
	    else if (otree->tau_decay_mode_2 ==10){ shift_tes=shift_tes_3prong; }
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
	  float mtTOT = 2*(otree->pt_1)*metLV.Pt()*(1-cos(otree->phi_1 - otree->mvametphi));
	  mtTOT += 2*(otree->pt_2)*metLV.Pt()*(1-cos(otree->phi_2 - otree->mvametphi)); 
	  mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(otree->phi_1-otree->phi_2)); 
	  otree->mt_tot = TMath::Sqrt(mtTOT);

      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;

      // dilepton veto
      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
	  if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree);

	  //extra letpn veto
	  otree->extraelec_veto = extra_electron_veto(leptonIndex, ch, &cfg, &analysisTree);
      otree->extramuon_veto = extra_muon_veto(leptonIndex, ch, &cfg, &analysisTree, isData);


      // define MET covariance
	  /*TMatrixD covMET(2, 2);
      covMET[0][0] = otree->mvacov00;
      covMET[1][0] = otree->mvacov10;
      covMET[0][1] = otree->mvacov01;
      covMET[1][1] = otree->mvacov11;*/

      // mt calculation and filling. MVA MET is used. 
      mt_calculation(otree);

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
      //otree->pzetamiss = otree->mvamet*TMath::Cos(otree->mvametphi)*zetaX + otree->mvamet*TMath::Sin(otree->mvametphi)*zetaY;
      otree->pzetamiss = otree->met*TMath::Cos(otree->metphi)*zetaX + otree->met*TMath::Sin(otree->metphi)*zetaY;
      //otree->pfpzetamiss = analysisTree.pfmetcorr_ex*zetaX + analysisTree.pfmetcorr_ey*zetaY; // this is not recoil-corrected  
      otree->pfpzetamiss = calc::pzetamiss( zetaX, zetaY, otree->met, otree->metphi);   
      otree->puppipzetamiss = analysisTree.puppimet_ex*zetaX + analysisTree.puppimet_ey*zetaY;  // this is not recoil-corrected  

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
      if (ch=="mt") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
                            otree->againstElectronVLooseMVA6_2>0.5 && otree->againstMuonTight3_2>0.5  &&
                            otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);
      if (ch=="et") 
        passedBaselineSel = ( otree->iso_1<0.35 && otree->byMediumIsolationMVArun2v1DBoldDMwLT_2>0.5 && 
                            otree->againstMuonLoose3_2>0.5 && otree->againstElectronTightMVA6_2>0.5 && 
                            otree->dilepton_veto == 0 && otree->extraelec_veto == 0 && otree->extramuon_veto == 0);


      if (!Synch && !passedBaselineSel)
        continue;

      if (ApplySVFit && otree->njetspt20>0) svfit_variables(ch, &analysisTree, otree, &cfg, inputFile_visPtResolution);

      otree->Fill();

	  // evaluate systematics for MC 
      if(!isData && ApplySystShift){
       zPtWeightSys->Eval(); 
       topPtWeightSys->Eval();
       for (unsigned int i=0; i<jetEnergyScaleSys.size(); i++)
         (jetEnergyScaleSys.at(i)) ->Eval(); 
         if (ch=="mt") {
           tauOneProngScaleSys->Eval(utils::MUTAU);
           tauOneProngOnePi0ScaleSys->Eval(utils::MUTAU);
           tauThreeProngScaleSys->Eval(utils::MUTAU);
           lepTauFakeScaleSys->Eval(utils::MUTAU);
         }
	 else if (ch=="et") {
           tauOneProngScaleSys->Eval(utils::ETAU);
           tauOneProngOnePi0ScaleSys->Eval(utils::ETAU);
           tauThreeProngScaleSys->Eval(utils::ETAU);
           lepTauFakeScaleSys->Eval(utils::ETAU);
         }
      }

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
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();

  // delete systematics objects

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

  if(lepTauFakeScaleSys != 0){
    lepTauFakeScaleSys->Write();
    delete lepTauFakeScaleSys;
  }

  file->Close();
  delete file;
    

}




///////////////////////////////////////////////
//////////////FUNCTION DEFINITION//////////////
///////////////////////////////////////////////

bool isICHEPmed(int Index, const AC1B * analysisTree) {

  bool goodGlob = analysisTree->muon_isGlobal[Index] && 
                  analysisTree->muon_normChi2[Index] < 3 && 
                  analysisTree->muon_combQ_chi2LocalPosition[Index] < 12 && 
                  analysisTree->muon_combQ_trkKink[Index] < 20;

  bool isICHEPmedium  = analysisTree->muon_isLoose[Index] &&
                        analysisTree->muon_validFraction[Index] >0.49 &&
                        analysisTree->muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
  
  return isICHEPmedium;

}

int read_json(std::string filename, lumi_json& json){

  std::pair <int,int> lumi;

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);

  BOOST_FOREACH(boost::property_tree::ptree::value_type &json_run, pt.get_child("")){
    int irun = atoi(json_run.first.data());
    BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_ranges, json_run.second.get_child("")){
      int ilumi[2] = {};

      int count = 0;
      BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_boundaries, lumi_ranges.second.get_child("")){
	ilumi[count] = atoi(lumi_boundaries.second.data().data());
	count++;
      }
      
      for (;ilumi[0] <= ilumi[1]; ilumi[0]++){
	lumi = std::make_pair(irun, ilumi[0]);
	json.push_back(lumi);
      }
    }
  }

  sort( json.begin(), json.end(),  compare_lumi());
  json.erase( unique( json.begin(), json.end() ), json.end() );
  
  return 0;
}


bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json){
  static compare_lumi compare;
  static std::pair<int,int> oldlumi = lumi;
  static bool old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  
  if(lumi.first != oldlumi.first || lumi.second != oldlumi.second){
    oldlumi = lumi;
    old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  }

  return old;
}

//accepts run number, lumi and json, make pair of (run,lumi) and starts isGoodLumi 
bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}

//fill the otree with the weights
void fill_weight(const AC1B * analysisTree, Spring15Tree *otree, PileUp *PUofficial, bool isData){
  
  otree->mcweight = 1.;
  otree->pu_weight = 0;

  if(!isData)
    otree->mcweight = analysisTree->genweight;

  otree->pu_weight = float(PUofficial->get_PUweight(double(analysisTree->numtruepileupinteractions)));

  otree->trigweight_1 = 0;
  otree->trigweight_2 = 1;
  otree->idisoweight_1 = 0;
  otree->idisoweight_2 = 1;
  otree->trigweight_antiiso_1 = 0;
  otree->idisoweight_antiiso_1 = 0;
  otree->trkeffweight_1=1;
  otree->effweight = 0;
  otree->weight = 1;
  otree->lheHt = analysisTree->genparticles_lheHt;
  otree->gen_noutgoing = analysisTree->genparticles_noutgoing;
}


// select medium id or ICHEP medium id for different runs, in data
// use ICHEP medium ID fro runs up to 278808 (end of Run2016F)
// on MC, always apply medium ID
bool isIdentifiedMediumMuon(int Index, const AC1B * analysisTree, bool isData){
  bool isGoodMuon;
  if (isData){
	if (analysisTree->event_run<=278808) isGoodMuon= isICHEPmed(Index, analysisTree);
    else isGoodMuon=analysisTree->muon_isMedium[Index]; 
	}
  else isGoodMuon=analysisTree->muon_isMedium[Index]; 
  return isGoodMuon;
}

//compute the absolute isolation for a given lepton labeled by Index in channel ch
float abs_Iso (int Index, TString ch, const AC1B * analysisTree, float dRiso){

  float neutralHadIso = -9999.;
  float photonIso     = -9999.;
  float chargedHadIso = -9999.;
  float puIso         = -9999.;

  if(ch=="mt"){
    neutralHadIso = analysisTree->muon_neutralHadIso[Index];
    photonIso =     analysisTree->muon_photonIso[Index];
    chargedHadIso = analysisTree->muon_chargedHadIso[Index];
    puIso =         analysisTree->muon_puIso[Index];
    if (dRiso>0.29) {
      neutralHadIso =     analysisTree->muon_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r03_sumPUPt[Index];
    }
    if (dRiso>0.39) {
      neutralHadIso =     analysisTree->muon_r04_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r04_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r04_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r04_sumPUPt[Index];
    }
  }
  if(ch=="et"){
    neutralHadIso = analysisTree->electron_neutralHadIso[Index];
    photonIso =     analysisTree->electron_photonIso[Index];
    chargedHadIso = analysisTree->electron_chargedHadIso[Index];
    puIso =         analysisTree->electron_puIso[Index];
    if (dRiso>0.29) {
      neutralHadIso =     analysisTree->electron_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->electron_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->electron_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->electron_r03_sumPUPt[Index];
    }
  }

  float neutralIso = neutralHadIso + photonIso -0.5*puIso;
  neutralIso = TMath::Max(float(0), neutralIso);
  return(chargedHadIso + neutralIso);
}

//compute the relative isolation for a given lepton labeled by Index in channel ch
float rel_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso){
  if(ch=="mt")  return(abs_Iso(Index, ch, analysisTree, dRiso) / analysisTree->muon_pt[Index] );
  else if(ch=="et")   return(abs_Iso(Index, ch, analysisTree, dRiso) / analysisTree->electron_pt[Index] );
    else return(-1.);
}


////FILLING FUNCTIONS//////

//fill the otree with the muon variables in channel mutau
void FillMuTau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->muon_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->muon_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->muon_phi[leptonIndex];
  otree->m_1 = 		muonMass;
  otree->q_1 = -1;
  if (analysisTree->muon_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->muon_genmatch[leptonIndex];

	otree->iso_1 = abs_Iso(leptonIndex, "mt", analysisTree, dRiso) / analysisTree->muon_pt[leptonIndex];

	otree->d0_1 = analysisTree->muon_dxy[leptonIndex];
	otree->dZ_1 = analysisTree->muon_dz[leptonIndex];

}

//fill the otree with the electron variables in channel etau
void FillETau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->electron_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->electron_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->electron_phi[leptonIndex];
  otree->m_1 = 		electronMass;
    otree->q_1 = -1;
  if (analysisTree->electron_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->electron_genmatch[leptonIndex];

  otree->iso_1 =  abs_Iso(leptonIndex, "et", analysisTree, dRiso) /analysisTree->electron_pt[leptonIndex];
  otree->mva_1 = analysisTree->electron_mva_value_Spring16_v1[leptonIndex];//analysisTree->electron_mva_value_nontrig_Spring15_v1[leptonIndex];

	otree->d0_1 = analysisTree->electron_dxy[leptonIndex];
	otree->dZ_1 = analysisTree->electron_dz[leptonIndex];

}

//fill the otree with the tau variables 
void FillTau(const AC1B * analysisTree, Spring15Tree *otree, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->eta_2 = analysisTree->tau_eta[tauIndex];
  otree->phi_2 = analysisTree->tau_phi[tauIndex];
  otree->q_2 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_2 = analysisTree->tau_genmatch[tauIndex];
  //if (analysisTree->tau_charge[tauIndex]>0)
  //otree->q_2 = 1;
  //otree->mva_2 = log(0);
  otree->mva_2 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex]; //before, tau_byIsolationMVArun2v1DBoldDMwLTraw
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
  otree->isolationGammaCands_size_2 = analysisTree->tau_isolationGammaCands_size[tauIndex];
  otree->signalGammaCands_size_2 = analysisTree->tau_signalGammaCands_size[tauIndex];
}

//////DILEPTON FUNCTIONS

//returns the dilepton veto for mt channel
bool dilepton_veto_mt(const Config *cfg,const  AC1B *analysisTree){

  for (unsigned int im = 0; im<analysisTree->muon_count; ++im) {

    if (analysisTree->muon_pt[im]<=cfg->get<float>("ptDiMuonVeto")) continue;
    if (fabs(analysisTree->muon_eta[im])>=cfg->get<float>("etaDiMuonVeto")) continue;	
		
    if (fabs(analysisTree->muon_dxy[im])>=cfg->get<float>("dxyDiMuonVeto")) continue;
    if (fabs(analysisTree->muon_dz[im])>=cfg->get<float>("dzDiMuonVeto")) continue;

    float relIsoMu = rel_Iso(im, "mt", analysisTree, cfg->get<float>("dRiso"));
    if(relIsoMu >= cfg->get<float>("isoDiMuonVeto")) continue;
		
    //bool passedVetoId = isICHEPmed(im, analysisTree);
    //if (!passedVetoId && cfg->get<bool>("applyDiMuonVetoId")) continue;
    if ( !(analysisTree->muon_isGlobal[im] && analysisTree->muon_isTracker[im] && analysisTree->muon_isPF[im]) ) continue;
    
    for (unsigned int je = im+1; je<analysisTree->muon_count; ++je) {
      
      if (analysisTree->muon_pt[je]<=cfg->get<float>("ptDiMuonVeto")) continue;
      if (fabs(analysisTree->muon_eta[je])>=cfg->get<float>("etaDiMuonVeto")) continue;	
		  
      if (fabs(analysisTree->muon_dxy[je])>=cfg->get<float>("dxyDiMuonVeto")) continue;
      if (fabs(analysisTree->muon_dz[je])>=cfg->get<float>("dzDiMuonVeto")) continue;
		  
      if (analysisTree->muon_charge[im] * analysisTree->muon_charge[je] > 0. && cfg->get<bool>("applyDiMuonOS")) continue;

      float relIsoMu = 	rel_Iso(je, "mt", analysisTree, cfg->get<float>("dRiso"));
      if(relIsoMu >= cfg->get<float>("isoDiMuonVeto")) continue;	

      //passedVetoId = isICHEPmed(je, analysisTree);
      //if (!passedVetoId && cfg->get<bool>("applyDiMuonVetoId")) continue;
      if ( ! (analysisTree->muon_isGlobal[je] && analysisTree->muon_isTracker[je] && analysisTree->muon_isPF[je]) ) continue;
		  
      float dr = deltaR(analysisTree->muon_eta[im],analysisTree->muon_phi[im],analysisTree->muon_eta[je],analysisTree->muon_phi[je]);

      if(dr<=cfg->get<float>("drDiMuonVeto")) continue;

      return(1);
    }
  }
  return(0);
}


//returns the dilepton veto fot et channel
bool dilepton_veto_et(const Config *cfg,const  AC1B *analysisTree){

  for (unsigned int ie = 0; ie<analysisTree->electron_count; ++ie) {

    if (analysisTree->electron_pt[ie]<=cfg->get<float>("ptDiElectronVeto")) continue;
    if (fabs(analysisTree->electron_eta[ie])>=cfg->get<float>("etaDiElectronVeto")) continue;	
    if (fabs(analysisTree->electron_dxy[ie])>=cfg->get<float>("dxyDiElectronVeto")) continue;
    if (fabs(analysisTree->electron_dz[ie])>=cfg->get<float>("dzDiElectronVeto")) continue;

    float absIsoEle =   abs_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));
    float relIsoEle =   rel_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));
    if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;
    
    bool passedVetoId =  passedSummer16VetoId(analysisTree, ie); //analysisTree->electron_cutId_veto_Summer16[ie];
    if (!passedVetoId && cfg->get<bool>("applyDiElectronVetoId")) continue;
		
    for (unsigned int je = ie+1; je<analysisTree->electron_count; ++je) {

      if (analysisTree->electron_pt[je]<=cfg->get<float>("ptDiElectronVeto")) continue;
      if (fabs(analysisTree->electron_eta[je])>=cfg->get<float>("etaDiElectronVeto")) continue;	
      if (fabs(analysisTree->electron_dxy[je])>=cfg->get<float>("dxyDiElectronVeto")) continue;
      if (fabs(analysisTree->electron_dz[je])>=cfg->get<float>("dzDiElectronVeto")) continue;
		  
      float absIsoEle =  abs_Iso(je, "et", analysisTree, cfg->get<float>("dRiso"));
      float relIsoEle =  rel_Iso(je, "et", analysisTree, cfg->get<float>("dRiso"));
      if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;	

      passedVetoId =  analysisTree->electron_cutId_veto_Summer16[je];
      if (!passedVetoId && cfg->get<bool>("applyDiElectronVetoId")) continue;

      if (analysisTree->electron_charge[ie] * analysisTree->electron_charge[je] > 0. && cfg->get<bool>("applyDiElectronOS")) continue;
		  
      float dr = deltaR(analysisTree->electron_eta[ie],analysisTree->electron_phi[ie],
                        analysisTree->electron_eta[je],analysisTree->electron_phi[je]);

      if(dr<=cfg->get<float>("drDiElectronVeto")) continue;

      return(1);
    }
  }
  return(0);
}

//////EXTRA LEPTON VETO FUNCTIONS

//returns the extra electron veto
bool extra_electron_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree){

  for (unsigned int ie = 0; ie<analysisTree->electron_count; ++ie) {

    if (ch=="et") {if (int(ie)==leptonIndex) continue;}

    if (analysisTree->electron_pt[ie]<=cfg->get<float>("ptVetoElectronCut")) continue;
    if (fabs(analysisTree->electron_eta[ie])>=cfg->get<float>("etaVetoElectronCut")) continue;
    if (fabs(analysisTree->electron_dxy[ie])>=cfg->get<float>("dxyVetoElectronCut")) continue;
    if (fabs(analysisTree->electron_dz[ie])>=cfg->get<float>("dzVetoElectronCut")) continue;

    bool electronMvaId = analysisTree->electron_mva_wp90_general_Spring16_v1[ie]; //analysisTree->electron_mva_wp90_nontrig_Spring15_v1[ie];
    if (!electronMvaId && cfg->get<bool>("applyVetoElectronId")) continue;
    if (!analysisTree->electron_pass_conversion[ie] && cfg->get<bool>("applyVetoElectronId")) continue;
    if (analysisTree->electron_nmissinginnerhits[ie]>1 && cfg->get<bool>("applyVetoElectronId")) continue;

    float relIsoEle = rel_Iso(ie, "et", analysisTree, cfg->get<float>("dRisoExtraElecVeto"));
    if (relIsoEle>=cfg->get<float>("isoVetoElectronCut")) continue;

    return(1);		
  }
  return(0);
}			

//returns the extra muon veto
bool extra_muon_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree, bool isData){

  for (unsigned int im = 0; im<analysisTree->muon_count; ++im) {

    if (ch=="mt") {if (int(im)==leptonIndex) continue;}

    if (analysisTree->muon_pt[im]<cfg->get<float>("ptVetoMuonCut")) continue;
    if (fabs(analysisTree->muon_eta[im])>cfg->get<float>("etaVetoMuonCut")) continue;
    if (fabs(analysisTree->muon_dxy[im])>cfg->get<float>("dxyVetoMuonCut")) continue;
    if (fabs(analysisTree->muon_dz[im])>cfg->get<float>("dzVetoMuonCut")) continue;

    if (cfg->get<bool>("applyVetoMuonId") && !(isIdentifiedMediumMuon(im,analysisTree,isData)) ) continue;
    float relIsoMu = rel_Iso(im, "mt", analysisTree, cfg->get<float>("dRisoExtraMuonVeto"));
    if (relIsoMu>cfg->get<float>("isoVetoMuonCut")) continue;

    return(1);
  }
  return(0);
}


//////MET FUNCTIONS

//fill the otree with the met variables
void fillMET(TString ch, int leptonIndex, int tauIndex, const AC1B * analysisTree, Spring15Tree *otree){

   // pfmet variables
  
  otree->met = TMath::Sqrt(analysisTree->pfmetcorr_ex*analysisTree->pfmetcorr_ex + analysisTree->pfmetcorr_ey*analysisTree->pfmetcorr_ey);
  otree->metphi = TMath::ATan2(analysisTree->pfmetcorr_ey,analysisTree->pfmetcorr_ex);
  otree->metcov00 = analysisTree->pfmetcorr_sigxx;
  otree->metcov01 = analysisTree->pfmetcorr_sigxy;
  otree->metcov10 = analysisTree->pfmetcorr_sigyx;
  otree->metcov11 = analysisTree->pfmetcorr_sigyy;
  float met_x = analysisTree->pfmetcorr_ex;
  float met_y = analysisTree->pfmetcorr_ey;

  float met_x2 = met_x * met_x;
  float met_y2 = met_y * met_y;

  // puppimet variables
  otree->puppimet = TMath::Sqrt(analysisTree->puppimet_ex*analysisTree->puppimet_ex +
		 analysisTree->puppimet_ey*analysisTree->puppimet_ey);
  otree->puppimetphi = TMath::ATan2(analysisTree->puppimet_ey,analysisTree->puppimet_ex);

  // choosing mva met
  int mva_break = -1;
  if (ch=="mt") mva_break = 3;
  else if (ch=="et") mva_break = 2;

  unsigned int iMet = 0;
  float mvamet_x = 0;
  float mvamet_y = 0;
  otree->mvacov00 = 0.;
  otree->mvacov01 = 0.;
  otree->mvacov10 = 0.;
  otree->mvacov11 = 0.;

  for (; iMet<analysisTree->mvamet_count; ++iMet) {
    if ((int)analysisTree->mvamet_channel[iMet]==mva_break) break;
  }

  if (iMet>=analysisTree->mvamet_count){
    otree->mvamet = log(0);
    otree->mvametphi = log(0);
    otree->mvacov00 = log(0);
    otree->mvacov01 = log(0);
    otree->mvacov10 = log(0);
    otree->mvacov11 = log(0);
  }
  else {
    // choosing mva met
    unsigned int iMet = 0;
    for (; iMet<analysisTree->mvamet_count; ++iMet) {
      if (analysisTree->mvamet_channel[iMet]==mva_break){
	if( ((int)analysisTree->mvamet_lep1[iMet])==tauIndex && ((int)analysisTree->mvamet_lep2[iMet])==leptonIndex)
	  break;
      }
    }
    
    float mvamet_x = 0;
    float mvamet_y = 0;
    otree->mvacov00 = 0.;
    otree->mvacov01 = 0.;
    otree->mvacov10 = 0.;
    otree->mvacov11 = 0.;
		
    if(iMet < analysisTree->mvamet_count){
      mvamet_x = analysisTree->mvamet_ex[iMet];
      mvamet_y = analysisTree->mvamet_ey[iMet];
      otree->mvacov00 = analysisTree->mvamet_sigxx[iMet];
      otree->mvacov01 = analysisTree->mvamet_sigxy[iMet];
      otree->mvacov10 = analysisTree->mvamet_sigyx[iMet];
      otree->mvacov11 = analysisTree->mvamet_sigyy[iMet];
    }
    else{
      std::cout<<"MVA MET not found! :: "<< analysisTree->mvamet_count << std::endl;
      iMet = 0;
      //std::cout<<"tau = "<<tauIndex<<" mu = "<<muonIndex<<std::endl;
      //for (; iMet<analysisTree->mvamet_count; ++iMet) {
      //std::cout<<"imet = "<<analysisTree->mvamet_channel[iMet]<<" itau = "<<analysisTree->mvamet_lep1[iMet]<<" imu = "<<analysisTree->mvamet_lep2[iMet]<<std::endl;
      //}
    }
		
    float mvamet_x2 = mvamet_x * mvamet_x;
    float mvamet_y2 = mvamet_y * mvamet_y;
    otree->mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
    otree->mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
  }

}


//caltulates and fill the otree with the MT variables
void mt_calculation(Spring15Tree *otree){
  //otree->mt_1 = sqrt(2*otree->pt_1*otree->mvamet*(1.-cos(otree->phi_1-otree->mvametphi)));
  //otree->mt_2 = sqrt(2*otree->pt_2*otree->mvamet*(1.-cos(otree->phi_2-otree->mvametphi)));
  otree->mt_1 = sqrt(2*otree->pt_1*otree->met*(1.-cos(otree->phi_1-otree->metphi)));
  otree->mt_2 = sqrt(2*otree->pt_2*otree->met*(1.-cos(otree->phi_2-otree->metphi)));

  
  otree->pfmt_1 = sqrt(2*otree->pt_1*otree->met*(1.-cos(otree->phi_1-otree->metphi)));
  otree->pfmt_2 = sqrt(2*otree->pt_2*otree->met*(1.-cos(otree->phi_2-otree->metphi)));

  otree->puppimt_1 = sqrt(2*otree->pt_1*otree->puppimet*(1.-cos(otree->phi_1-otree->puppimetphi)));
  otree->puppimt_2 = sqrt(2*otree->pt_2*otree->puppimet*(1.-cos(otree->phi_2-otree->puppimetphi)));
}

///////////////////////////////////
// SV fit 
///////////////////////////////////

void svfit_variables(TString ch, const AC1B *analysisTree, Spring15Tree *otree, const Config *cfg, TFile * inputFile_visPtResolution){

  // define MET covariance
  TMatrixD covMET(2, 2);
  // using MVA MET
  /*covMET[0][0] = otree->mvacov00;
  covMET[1][0] = otree->mvacov10;
  covMET[0][1] = otree->mvacov01;
  covMET[1][1] = otree->mvacov11;*/

  // using PF MET
  covMET[0][0] = otree->metcov00;
  covMET[1][0] = otree->metcov10;
  covMET[0][1] = otree->metcov01;
  covMET[1][1] = otree->metcov11;

  svFitStandalone::kDecayType type_ = svFitStandalone::kUndefinedDecayType;
  if (ch == "mt")      type_ = svFitStandalone::kTauToMuDecay;
  else if (ch == "et") type_ = svFitStandalone::kTauToElecDecay;
  // define lepton four vectors
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(type_,
								  otree->pt_1,
								  otree->eta_1,
								  otree->phi_1,
								  otree->m_1)); 
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay,
								  otree->pt_2,
								  otree->eta_2,
								  otree->phi_2,
								  otree->m_2,
								  otree->tau_decay_mode_2));

  //SVfitStandaloneAlgorithm algo(measuredTauLeptons, otree->mvamet * cos(otree->mvametphi), otree->mvamet * sin(otree->mvametphi), covMET, 0); // using MVA MET
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, otree->met * cos(otree->metphi), otree->met * sin(otree->metphi), covMET, 0); // using PF MET
  //edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  //TH1::AddDirectory(false);  
  //TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  //algo.addLogM(false);  
  algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();

  otree->m_sv   = algo.mass();
  otree->pt_sv  = algo.pt();
  otree->eta_sv = algo.eta();
  otree->phi_sv = algo.phi();      
  otree->met_sv = algo.fittedMET().Rho();
  otree->mt_sv  = algo.transverseMass();
/*
  std::cout << " SVFit in main  - Inputs -  " << std::endl;
  std::cout << " Lep 1 : pt = " << otree->pt_1 << " eta : " << otree->eta_1 << " phi : " << otree->phi_1 << " M : " << otree->m_1 <<   " DM : " << otree->tau_decay_mode_1 << std::endl;
  std::cout << " Lep 2 : pt = " << otree->pt_2 << " eta : " << otree->eta_2 << " phi : " << otree->phi_2 << " M : " << otree->m_2 << " DM : " << otree->tau_decay_mode_2 << std::endl;
  std::cout << " MET : px = " <<  otree->met * cos(otree->metphi) << " py = " << otree->met * sin(otree->metphi) << std::endl;
  std::cout << " MET COV : 00 = " << otree->metcov00 << " 01 " << otree->metcov01  << " 10 " << otree->metcov10 << " 11 " << otree->metcov11 << std::endl;
  std::cout << " SVFit in main  - Output -  " << std::endl;
  std::cout << " Mass sv = " << otree->m_sv << "pt sv = " << otree->pt_sv << std::endl;  
  std::cout << " SVFit in main  - END -  " << std::endl;
*/
  //if ( algo.isValidSolution() ) {
  //  std::cout << "found mass = " << mass << std::endl;
  //} else {
  //  std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  //}

  //delete inputFile_visPtResolution;
}


/// shift tau energy scale and propagate it to the met. 
void correctTauES(TLorentzVector& Tau, TLorentzVector& Met, float relative_shift, bool tau_is_one_prong){
  Met.SetPx(Met.Px()- (Tau.Px()*relative_shift) ) ;
  Met.SetPy(Met.Py()- (Tau.Py()*relative_shift) );
  if (tau_is_one_prong){ // don't scale the mass
  	Tau.SetXYZM( Tau.Px()*(1+relative_shift), Tau.Py()*(1+relative_shift), Tau.Pz()*(1+relative_shift), Tau.M());
  }
  else {
    Tau.SetXYZM( Tau.Px()*(1+relative_shift), Tau.Py()*(1+relative_shift), Tau.Pz()*(1+relative_shift), Tau.M()*(1+relative_shift));
  }
}

// Implementation of cut-based veto ID Summer 16 (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Electron_ID_Working_Points_WP_de)
bool passedSummer16VetoId(const AC1B * analysisTree, int index){

  if(abs(analysisTree->electron_superclusterEta[index]) <= 1.479){

    if( analysisTree->electron_full5x5_sigmaietaieta[index] >= 0.0155 )                       return false;
    if( abs(analysisTree->electron_detaInSeed[index]) >= 0.00749 )                            return false;
    if( abs(analysisTree->electron_deltaphisuperclustertrack[index]) >= 0.228 )               return false;
    if( analysisTree->electron_he[index] >= 0.356 )                                           return false;
    if( analysisTree->electron_eaIsolation[index]/analysisTree->electron_pt[index] >= 0.175 ) return false;
    if( abs(analysisTree->electron_ooemoop[index]) >= 0.299 )                                 return false;
    if( analysisTree->electron_nmissinginnerhits[index] > 2 )                                 return false;
    if( analysisTree->electron_pass_conversion[index] == false )                              return false;
    
    return true;
  }
  else{

    if( analysisTree->electron_full5x5_sigmaietaieta[index] >= 0.037 )                        return false;
    if( abs(analysisTree->electron_detaInSeed[index]) >= 0.00895 )                            return false;
    if( abs(analysisTree->electron_deltaphisuperclustertrack[index]) >= 0.213 )               return false;
    if( analysisTree->electron_he[index] >= 0.211 )                                           return false;
    if( analysisTree->electron_eaIsolation[index]/analysisTree->electron_pt[index] >= 0.159 ) return false;
    if( abs(analysisTree->electron_ooemoop[index]) >= 0.15 )                                  return false;
    if( analysisTree->electron_nmissinginnerhits[index] > 3 )                                 return false;
    if( analysisTree->electron_pass_conversion[index] == false )                              return false;
    
    return true;
  }
}


bool SafeRatio(double denominator){

  if(denominator<0.001) return false;
  else                  return true;

}
