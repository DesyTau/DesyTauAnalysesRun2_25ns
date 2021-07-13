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

#include "TVector3.h"
#include "TMatrix.h"

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
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PuppiMETSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PFMETSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/BtagSys_WIP.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"
//#include "HTT-utilities/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h"

//#include "DesyTauAnalyses/NTupleMaker/interface/ImpactParameter.h"
#include "HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "RooFunctor.h"

#define pi   3.14159265358979312
#define d2r  1.74532925199432955e-02
#define r2d  57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 	 0.105658
#define tauMass 	 1.77682
#define pionMass 	 0.1396

#define expectedtauspinnerweights 5

void initializeGenTree(Synch17GenTree *gentree);
void FillTauTau(const AC1B * analysisTree, Synch17Tree *otree, int tau1Index, int tau2Index, float shift_ts_1, float shift_ts_2);
void FillVertices(const AC1B * analysisTree,Synch17Tree *otree, const bool isData, int leptonIndex, int tauIndex, TString channel);
void FillGenTree(const AC1B * analysisTree, Synch17GenTree *gentree);
float getEmbeddedWeight(const AC1B * analysisTree, RooWorkspace* WS);
float shift_tauES(const AC1B * analysisTree, unsigned int itau, 
		  float shift_tes_1prong, 
		  float shift_tes_1p1p0, 
		  float shift_tes_3prong,
		  float shift_tes_3prong1p0,
		  float shift_tes_lepfake_1prong_barrel, 
		  float shift_tes_lepfake_1p1p0_barrel,
		  float shift_tes_lepfake_1prong_endcap, 
		  float shift_tes_lepfake_1p1p0_endcap
		  ); 
bool accessTriggerInfo(const AC1B * analysisTree, TString HLTFilterName, unsigned int &nHLTFilter)
{
   bool isHLTFilter = false;
   
   for (unsigned int i=0; i<analysisTree->run_hltfilters->size(); ++i) {
      TString HLTFilter(analysisTree->run_hltfilters->at(i));
      if (HLTFilter==HLTFilterName) {
         nHLTFilter = i;
         isHLTFilter = true;
      }
   }
   return isHLTFilter;
}
bool triggerMatching(AC1B * analysisTree, Float_t eta, Float_t phi, bool isFilter, unsigned int nFilter, float deltaRTrigMatch = 0.5)
{
   bool trigMatch = false;
   if (!isFilter) return trigMatch;
   for (unsigned int iT=0; iT<analysisTree->trigobject_count; ++iT) {
      float dRtrig = deltaR(eta,phi,analysisTree->trigobject_eta[iT],analysisTree->trigobject_phi[iT]);
      if (dRtrig> deltaRTrigMatch) continue;
      if (analysisTree->trigobject_filters[iT][nFilter]) trigMatch = true;
      
   }
   
   return trigMatch;
}

int main(int argc, char * argv[]){
// first argument - config file for analysis
// second argument - file list (MUST BE IN THE SAME DIRECTORY OF THE EXECUTABLE)
// third argument - index of first file to run on (optional, ignored if only one file is used)
// forth argument - index of last file to run on (optional, ignored if only one file is used)

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

  if(argc < 3){
    std::cout << "RUN ERROR: wrong number of arguments"<< std::endl;
    std::cout << "Please run the code in the following way:"<< std::endl;
    std::cout << "SynchNTupleProducer_Run2 NameOfTheConfigurationFile FileList" << std::endl;
    std::cout << "example: SynchNTupleProducer_tt_Run2 analysisMacroSynch_tt_18_data.conf Tau_Run2018A" << std::endl;
    exit(-1);
  }

  // **** configuration analysis  
  Config cfg(argv[1]);

  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];
  TString ch("tt");
  std::string lep;

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase) + "/src/" + TString(json_name)).Data(), json);
  }

  const int era = cfg.get<int>("era");
  const bool Synch = cfg.get<bool>("Synch"); 
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyFastMTT     = cfg.get<bool>("ApplyFastMTT");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters  = cfg.get<bool>("ApplyMetFilters");
  const bool usePuppiMET      = cfg.get<bool>("UsePuppiMET");
  const bool ApplyBTagCP5Correction = cfg.get<bool>("ApplyBTagCP5Correction");

  // JER
  //  const string jer_resolution = cfg.get<string>("JER_Resolution");
  //  const string jer_scalefactor = cfg.get<string>("JER_ScaleFactor");

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");

  // tau trigger efficiency
  std::string channel;
  if (ch == "mt") channel = "mutau"; 
  if (ch == "et") channel = "etau";
  if (ch == "tt") channel = "tautau";

  std::string year_label;
  if (era == 2016) year_label = "2016";
  else if (era == 2017) year_label = "2017";
  else if (era == 2018) year_label = "2018";	
  else {std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n'; exit(-1);}

  //svfit
  const string svFitPtResFile = TString(TString(cmsswBase) + "/src/" + TString(cfg.get<string>("svFitPtResFile"))).Data();

  //zptweight file 
  const string ZptweightFile = cfg.get<string>("ZptweightFile");

  //b-tag scale factors
  const string BTagAlgorithm = cfg.get<string>("BTagAlgorithm");
  const string BtagSfFile = cmsswBase + "/src/" + cfg.get<string>("BtagSfFile");
  if( ApplyBTagScaling && gSystem->AccessPathName( (TString) BtagSfFile) ){
    cout<<BtagSfFile<<" not found. Please check."<<endl;
    exit(-1);
  }
  
  // JER
  std::unique_ptr<JME::JetResolution> m_resolution_from_file;
  std::unique_ptr<JME::JetResolutionScaleFactor> m_scale_factor_from_file;
  if (era==2016) {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Summer16_25nsV1_MC_SF_AK4PFchs.txt"));
  }
  else if (era==2017) {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Fall17_V3_MC_SF_AK4PFchs.txt"));
  }
  else {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Autumn18_V7b_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Autumn18_V7b_MC_SF_AK4PFchs.txt"));    
  }

  JME::JetResolution resolution = *m_resolution_from_file;
  JME::JetResolutionScaleFactor resolution_sf = *m_scale_factor_from_file;

  cout<<"using "<<BTagAlgorithm<<endl;
  BTagCalibration calib;
  BTagCalibrationReader reader_B;
  BTagCalibrationReader reader_C;
  BTagCalibrationReader reader_Light;
  if(ApplyBTagScaling){
    calib = BTagCalibration(BTagAlgorithm, BtagSfFile);
    reader_B = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
    reader_C = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
    reader_Light = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
    reader_B.load(calib, BTagEntry::FLAV_B, "comb");
    reader_C.load(calib, BTagEntry::FLAV_C, "comb");
    reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  }
    
  TString pathToTaggingEfficiencies = (TString) cmsswBase + "/src/" + cfg.get<string>("BtagMCeffFile");
  if (ApplyBTagScaling && gSystem->AccessPathName(pathToTaggingEfficiencies)){
    cout<<pathToTaggingEfficiencies<<" not found. Please check."<<endl;
    exit(-1);
  }
    
  TFile *fileTagging  = new TFile(pathToTaggingEfficiencies);
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  TH2F  *tagEff_B_nonCP5     = 0;
  TH2F  *tagEff_C_nonCP5     = 0;
  TH2F  *tagEff_Light_nonCP5 = 0;
  TRandom3 *rand = new TRandom3();

  if(ApplyBTagScaling){
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
    if (ApplyBTagCP5Correction) {
      TString pathToTaggingEfficiencies_nonCP5 = (TString) cmsswBase + "/src/" + cfg.get<string>("BtagMCeffFile_nonCP5");
      if (gSystem->AccessPathName(pathToTaggingEfficiencies_nonCP5)) {
        cout<<pathToTaggingEfficiencies_nonCP5<<" not found. Please check."<<endl;
        exit(-1);
      } 
      TFile *fileTagging_nonCP5  = new TFile(pathToTaggingEfficiencies_nonCP5);
      tagEff_B_nonCP5     = (TH2F*)fileTagging_nonCP5->Get("btag_eff_b");
      tagEff_C_nonCP5     = (TH2F*)fileTagging_nonCP5->Get("btag_eff_c");
      tagEff_Light_nonCP5 = (TH2F*)fileTagging_nonCP5->Get("btag_eff_oth");
    }
  }  
  const struct btag_scaling_inputs inputs_btag_scaling_medium = {reader_B, reader_C, reader_Light, tagEff_B, tagEff_C, tagEff_Light, tagEff_B_nonCP5, tagEff_C_nonCP5, tagEff_Light_nonCP5, rand};

  TFile * ff_file = TFile::Open(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_tt_lite_"+TString(year_label)+"_dR_corr.root");
  if (ff_file->IsZombie()) {
    cout << "File " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_tt_lite_" << TString(year_label) << "_dR_corr.root not found" << endl;
    cout << "Quitting... " << endl;
    exit(-1);
  }

  std::shared_ptr<RooWorkspace> ff_ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  ff_ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  fns_["ff_tt_medium_dmbins"] = std::shared_ptr<RooFunctor>(ff_ws_->function("ff_tt_medium_dmbins")->functor(ff_ws_->argSet("pt,dm,njets,os,met_var_qcd,dR")));

  // MET Recoil Corrections
  const bool isDY = (infiles.find("DY") != string::npos) || (infiles.find("EWKZ") != string::npos);//Corrections that should be applied on EWKZ are the same needed for DY
  const bool isWJets = (infiles.find("WJets") != string::npos) || (infiles.find("W1Jets") != string::npos) || (infiles.find("W2Jets") != string::npos) || (infiles.find("W3Jets") != string::npos) || (infiles.find("W4Jets") != string::npos) || (infiles.find("EWK") != string::npos);
  const bool isHiggs = (infiles.find("VBFHTo")!= string::npos) || (infiles.find("WminusHTo")!= string::npos) || (infiles.find("WplusHTo")!= string::npos) || (infiles.find("ZHTo")!= string::npos) || (infiles.find("GluGluHTo")!= string::npos); 
  const bool isEWKZ =  infiles.find("EWKZ") != string::npos;
  const bool isMG = infiles.find("madgraph") != string::npos;
  const bool isMSSMsignal =  (infiles.find("SUSYGluGluToHToTauTau")!= string::npos) || (infiles.find("SUSYGluGluToBBHToTauTau")!= string::npos);
  const bool isTauSpinner = infiles.find("Uncorr") != string::npos;
  const bool isTTbar = infiles.find("TT") != string::npos;

  bool applyTauSpinnerWeights = false;
  if(isTauSpinner) applyTauSpinnerWeights = true;
  const bool isEmbedded = infiles.find("Embed") != string::npos;

  const bool ApplyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections") && !isEmbedded && !isData && (isDY || isWJets || isHiggs || isMSSMsignal);
  kit::RecoilCorrector recoilCorrector(cfg.get<string>("RecoilFilePath"));
  kit::MEtSys MetSys(cfg.get<string>("RecoilSysFilePath"));
  
  // tau cuts
  const float ptTauCut    = cfg.get<float>("ptTauCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut       = cfg.get<float>("dzTauCut");

  // tau energy scale corrections
  const float shift_tes_1prong     = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0      = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong     = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  const float shift_tes_3prong1p0 = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0");

  const float shift_tes_1prong_e = cfg.get<float>("TauEnergyScaleShift_OneProng_Error");
  const float shift_tes_1p1p0_e  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_Error");
  const float shift_tes_3prong_e = cfg.get<float>("TauEnergyScaleShift_ThreeProng_Error");
  const float shift_tes_3prong1p0_e = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0_Error");
  
  std::string year_label_fes;
  if (era == 2016) year_label_fes = "2016Legacy";
  else if (era == 2017) year_label_fes = "2017ReReco";
  else if (era == 2018) year_label_fes = "2018ReReco";	
  else {std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n'; exit(-1);}

  // lep->tau FES correction and uncertainties
  TFile TauFES_file(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_"+year_label_fes+".root"); 
  TGraphAsymmErrors* FES_graph = (TGraphAsymmErrors*) TauFES_file.Get("fes");
	
  // apply non-zero only for DY MC in et channel, HPS decay modes 0 and 1	
  const float shift_tes_lepfake_1prong_barrel =  FES_graph->GetY()[0];
  const float shift_tes_lepfake_1p1p0_barrel  = FES_graph->GetY()[1];
  const float shift_tes_lepfake_1prong_endcap = FES_graph->GetY()[2];
  const float shift_tes_lepfake_1p1p0_endcap  = FES_graph->GetY()[3];
  
  // for et take up/down values from file, for mt set to 1% (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2#Corrections_to_be_applied_to_AN2)
  // with the definition as below it can be accidently applied to data/embedded, so make sure to not do this!
  const float shift_tes_lepfake_1prong_barrel_up = FES_graph->GetErrorYhigh(0);
  const float shift_tes_lepfake_1p1p0_barrel_up  = FES_graph->GetErrorYhigh(1);
  const float shift_tes_lepfake_1prong_endcap_up = FES_graph->GetErrorYhigh(2);
  const float shift_tes_lepfake_1p1p0_endcap_up  = FES_graph->GetErrorYhigh(3);
	
  const float shift_tes_lepfake_1prong_barrel_down = FES_graph->GetErrorYlow(0);
  const float shift_tes_lepfake_1p1p0_barrel_down  = FES_graph->GetErrorYlow(1);
  const float shift_tes_lepfake_1prong_endcap_down = FES_graph->GetErrorYlow(2);
  const float shift_tes_lepfake_1p1p0_endcap_down  = FES_graph->GetErrorYlow(3);

  // lep->tau FR, DeepTau WPs
  const string leptauFake_wpVsEle = cfg.get<string>("LeptauFake_wpVsEle");
  const string leptauFake_wpVsMu = cfg.get<string>("LeptauFake_wpVsMu");
  TString LeptauFake_wpVsEle(leptauFake_wpVsEle);
  TString LeptauFake_wpVsMu(leptauFake_wpVsMu);

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

  const float dRTrigMatch = cfg.get<float>("dRTrigMatch");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");
  
  // Read in HLT filter
  vector<string> filterDiTau;
  vector<string> filterDiTau_before_HPS = cfg.get<vector<string>>("filterDiTau_before_HPS");
  vector<string> filterDiTau_HPS = cfg.get<vector<string>>("filterDiTau_HPS"); 
 
    // correction workspace
  const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");

  //file list creation
  int ifile = 0;
  int jfile = -1;

  if (argc > 3)
    ifile = atoi(argv[3]);
  if (argc > 4)
    jfile = atoi(argv[4]);
  
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
  std::string initNtupleName("initroottree/AC1B");

  // PU reweighting - initialization
  PileUp *PUofficial = new PileUp();
  if(ApplyPUweight){
    TFile *filePUdistribution_data = new TFile(TString(cmsswBase) + "/src/" + TString(pileUpInDataFile), "read");
    TFile *filePUdistribution_MC = new TFile (TString(cmsswBase) + "/src/" + TString(pileUpInMCFile), "read");
    TH1D *PU_data = (TH1D *)filePUdistribution_data->Get("pileup");    
    TH1D *PU_mc = (TH1D *)filePUdistribution_MC->Get(TString(pileUpforMC));
    if (PU_mc == NULL) {
      std::cout << "Histogram " << pileUpforMC << " is not present in pileup file" << std::endl;
      exit(-1);
    }
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  

  // Workspace with corrections
  TString workspace_filename = TString(cmsswBase) + "/src/" + CorrectionWorkspaceFileName;
  cout << "Taking correction workspace from " << workspace_filename << endl;
  TFile *f_workspace = new TFile(workspace_filename, "read");
  if (f_workspace->IsZombie()) {
    std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
     exit(-1);
   }
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");

  // Zpt reweighting for LO DY samples 
  TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
  TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

  // lepton to tau fake init
  TFile muTauFRfile(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSmu_"+year_label_fes+".root"); 
  TH1F * SF_muTau_hist = (TH1F*) muTauFRfile.Get(LeptauFake_wpVsMu);
  TFile eTauFRfile(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSe_"+year_label_fes+".root"); 
  TH1F * SF_eTau_hist = (TH1F*) eTauFRfile.Get(LeptauFake_wpVsEle);

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
  Synch17Tree *otree = new Synch17Tree(tree,ch,false);
  Synch17GenTree *gentree = new Synch17GenTree(gtree);
  Synch17GenTree *gentreeForGoodRecoEvtsOnly = new Synch17GenTree(tree);
    
  int nTotalFiles = 0;
  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  vector<unsigned int> runList; runList.clear();
  vector<unsigned int> eventList; eventList.clear();

  //svFit
  TH1::AddDirectory(false);  
  TFile *inputFile_visPtResolution = new TFile(svFitPtResFile.data());

  std::cout << "inputFile_visPtResolution : " << std::endl;

  //Systematics init
  
  TauScaleSys *tauScaleSys = 0;
  TauOneProngScaleSys *tauOneProngScaleSys = 0;
  TauOneProngOnePi0ScaleSys *tauOneProngOnePi0ScaleSys = 0;
  TauThreeProngScaleSys *tauThreeProngScaleSys = 0;
  TauThreeProngOnePi0ScaleSys *tauThreeProngOnePi0ScaleSys = 0;

  ZPtWeightSys* zPtWeightSys = 0;
  TopPtWeightSys* topPtWeightSys = 0;
  BtagSys * btagSys = 0;
  std::vector<JetEnergyScaleSys*> jetEnergyScaleSys;
  JESUncertainties * jecUncertainties = 0;

  std::vector<TString> metSysNames = {"CMS_scale_met_unclustered_13TeV"};
  std::vector<TString> recoilSysNames = {"CMS_htt_boson_reso_met_13TeV",
					 "CMS_htt_boson_scale_met_13TeV"};

  std::vector<PFMETSys*> metSys;
  std::vector<PuppiMETSys*> puppiMetSys;

  if((!isData||isEmbedded) && ApplySystShift){
    //    tauScaleSys = new TauScaleSys(otree);
    //    tauScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    //    tauScaleSys->SetUseSVFit(ApplySVFit);

    tauOneProngScaleSys = new TauOneProngScaleSys(otree);
    tauOneProngScaleSys->SetScale(shift_tes_1prong,shift_tes_1prong_e);
    tauOneProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauOneProngScaleSys->SetUseSVFit(ApplySVFit);
    tauOneProngScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauOneProngScaleSys->SetUsePuppiMET(usePuppiMET);

    tauOneProngOnePi0ScaleSys = new TauOneProngOnePi0ScaleSys(otree);
    tauOneProngOnePi0ScaleSys->SetScale(shift_tes_1p1p0,shift_tes_1p1p0_e);
    tauOneProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauOneProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
    tauOneProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauOneProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);

    tauThreeProngScaleSys = new TauThreeProngScaleSys(otree);
    tauThreeProngScaleSys->SetScale(shift_tes_3prong,shift_tes_3prong_e);
    tauThreeProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngScaleSys->SetUseSVFit(ApplySVFit);
    tauThreeProngScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauThreeProngScaleSys->SetUsePuppiMET(usePuppiMET);

    tauThreeProngOnePi0ScaleSys = new TauThreeProngOnePi0ScaleSys(otree);
    tauThreeProngOnePi0ScaleSys->SetScale(shift_tes_3prong1p0,shift_tes_3prong1p0_e);
    tauThreeProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
    tauThreeProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauThreeProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);

    // systematics only for MC
    if (!isEmbedded) {
      if (!isDY && !isWJets && !isHiggs) {
	btagSys = new BtagSys(otree,TString("Btag"));
	btagSys->SetConfig(&cfg);
	btagSys->SetBtagScaling(&inputs_btag_scaling_medium);
      }
      if (ApplyRecoilCorrections) {
	if (usePuppiMET) {
	  for (unsigned int i = 0; i < recoilSysNames.size(); ++i) {
	    PuppiMETSys * puppiMetRecoilSys = new PuppiMETSys(otree,recoilSysNames[i]);
	    puppiMetRecoilSys->SetMEtSys(&MetSys);
	    puppiMetSys.push_back(puppiMetRecoilSys);
	  }
	}
      }
      else {
	for (unsigned int i = 0; i<metSysNames.size(); ++i) {
	  if (usePuppiMET)
	    puppiMetSys.push_back(new PuppiMETSys(otree,metSysNames[i]));
	  else
	    metSys.push_back(new PFMETSys(otree,metSysNames[i]));
	}
      }
      if (cfg.get<bool>("splitJES")){
	JESUncertainties *jecUncertainties;
	if (era==2016) 
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
	else if (era==2017)
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
	else 
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
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
      JetEnergyScaleSys * JERsys = new JetEnergyScaleSys(otree, TString("JER"));
      JERsys->SetConfig(&cfg);
      JERsys->SetBtagScaling(&inputs_btag_scaling_medium);
      JERsys->SetJESUncertainties(jecUncertainties);
      jetEnergyScaleSys.push_back(JERsys);
    }
  }

  // list of met filters from config
  std::vector<TString> met_filters_list;
  for (unsigned int i = 1; i < (unsigned int) cfg.get<int>("num_met_filters") + 1; i++) {
    met_filters_list.push_back(cfg.get<string>("met_filter_" + std::to_string(i)));
  }

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
    // set AC1B for JES Btag and MET systematics
    if ( !isData && !isEmbedded && ApplySystShift) {
      if (!isDY && !isWJets && !isHiggs)
	btagSys->SetAC1B(&analysisTree);
      for (unsigned int i = 0; i < jetEnergyScaleSys.size(); i++)
      	(jetEnergyScaleSys.at(i))->SetAC1B(&analysisTree);
      for (unsigned int i = 0; i < metSys.size(); i++)
	(metSys.at(i))->SetAC1B(&analysisTree);
      for (unsigned int i = 0; i < puppiMetSys.size(); ++i)
	(puppiMetSys.at(i))->SetAC1B(&analysisTree);
    }
    
    //    double * TSweight = new double[expectedtauspinnerweights];
    //    TTree  * _treeTauSpinnerWeights = NULL;
    
    //    if(applyTauSpinnerWeights&&era==2017){ 
    //      _treeTauSpinnerWeights = (TTree*)file_->Get(TString(TauSpinnerWeightTreeName));
    //      _treeTauSpinnerWeights->SetBranchAddress("TauSpinnerWeights",TSweight);		
    //    }  
    

    TTree * _inittree = NULL;
    _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree!=NULL) {
        Float_t genweight;
        if (!isData)
            _inittree->SetBranchAddress("genweight",&genweight);
        Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
        std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
        for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
            _inittree->GetEntry(iEntry);
            if (isData && !isEmbedded)
                nWeightedEventsH->Fill(0.,1.);
            else
                nWeightedEventsH->Fill(0.,genweight);
        }
    }
    delete _inittree;

    ///////////////EVENT LOOP///////////////
    Long64_t numberOfEntries = analysisTree.GetEntries();

    for (Long64_t iEntry = 0; iEntry < numberOfEntries; iEntry++) {
      analysisTree.GetEntry(iEntry);
      nEvents++;

      filterDiTau = filterDiTau_before_HPS;
      if (era == 2018) {
      	if(isData && !isEmbedded) {
	  if (analysisTree.event_run < 317509)
	    filterDiTau = filterDiTau_before_HPS;
	  else
	    filterDiTau = filterDiTau_HPS;
	}
	else { // use with HPS in its name for MC and Embedded
      	  filterDiTau = filterDiTau_HPS;
      	}
      }    

      //      std::cout << "OK!!!!!!!!" << std::endl;
      //Skip events not passing the MET filters, if applied
      bool passed_all_met_filters = passedAllMetFilters(&analysisTree, met_filters_list);
      if (ApplyMetFilters && !Synch && !passed_all_met_filters) continue;
      otree->passedAllMetFilters = passed_all_met_filters;
      
      if (nEvents % 10000 == 0) 
      	cout << "      processed " << nEvents << " events" << endl; 
    
      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
    
      if ((isData || isEmbedded) && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;
    
      initializeGenTree(gentree);
      if (!isData) {
      	FillGenTree(&analysisTree,gentree);
      	gentree->Fill();
      }

      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;


      // tau selection
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it < analysisTree.tau_count; ++it) {
	float shift_tes  = shift_tauES(&analysisTree,it,
				       shift_tes_1prong,
				       shift_tes_1p1p0,
				       shift_tes_3prong,
				       shift_tes_3prong1p0,
				       shift_tes_lepfake_1prong_barrel,
				       shift_tes_lepfake_1p1p0_barrel,
				       shift_tes_lepfake_1prong_endcap,
				       shift_tes_lepfake_1p1p0_endcap);
 
	float ptTau = (1.+shift_tes)*analysisTree.tau_pt[it];

        if (ptTau < ptTauCut) continue;
        if (fabs(analysisTree.tau_eta[it]) >= etaTauCut) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it]) >= dzTauCut) continue;
        if (fabs(fabs(analysisTree.tau_charge[it]) - 1) > 0.001) continue;

      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] < 0.5) continue;
      	if (analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[it] < 0.5) continue;
      	if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] < 0.5) continue;
        if (analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) continue; //always true, cut applied in NTupleMaker
        if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6 || analysisTree.tau_decayMode[it] ==7) continue;
	//	if (analysisTree.tau_MVADM2017v1[it] < 0) continue; //prevents storing events with unidentified mva DM for the tau (-1)
    
        taus.push_back(it);
      }
    
      if (taus.size() < 2) continue;
    
      // selecting ditau pair;
      int tau1Index = -1;
      int tau2Index = -1;
    
      float isoTausMax   = -10;
      float shift_t1_es = 0.;
      float shift_t2_es = 0.;

      //////////////LOOP on Taus/////////////
    
      for (unsigned int it = 0; it < taus.size()-1; ++it) {
      	unsigned int tIndex = taus.at(it);
	float shift_tes_1 = 0.;
	if (!isData || isEmbedded ) shift_tes_1 = shift_tauES(&analysisTree,tIndex,
							      shift_tes_1prong,
							      shift_tes_1p1p0,
							      shift_tes_3prong,
							      shift_tes_3prong1p0,
							      shift_tes_lepfake_1prong_barrel,
							      shift_tes_lepfake_1p1p0_barrel,
							      shift_tes_lepfake_1prong_endcap,
							      shift_tes_lepfake_1p1p0_endcap);
	
	float ptTau1 = (1.+shift_tes_1)*analysisTree.tau_pt[tIndex];
	float sortIsoTau1 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex];
    
    	//////////////LOOP on Leptons or second Tau/////////////
    
        for (unsigned int il = it+1; il < taus.size(); ++il) {
          unsigned int lIndex = taus.at(il);

	  float shift_tes_2 = 0.; 
	  if (!isData || isEmbedded) shift_tes_2 = shift_tauES(&analysisTree,lIndex,
							       shift_tes_1prong,
							       shift_tes_1p1p0,
							       shift_tes_3prong,
							       shift_tes_3prong1p0,
							       shift_tes_lepfake_1prong_barrel,
							       shift_tes_lepfake_1p1p0_barrel,
							       shift_tes_lepfake_1prong_endcap,
							       shift_tes_lepfake_1p1p0_endcap);
 
	  float ptTau2 = (1.+shift_tes_2)*analysisTree.tau_pt[lIndex];        
          float sortIsoTau2 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex];
          float dR = deltaR(analysisTree.tau_eta[tIndex], analysisTree.tau_phi[tIndex], 
			    analysisTree.tau_eta[lIndex], analysisTree.tau_phi[lIndex]);
          if (dR < dRleptonsCut) continue;
    
          // change pair
          bool changePair =  false;
	  float sortIsoTaus = sortIsoTau1+sortIsoTau2;
          if (sortIsoTaus > isoTausMax)
	    changePair = true;
    
          if (changePair) {
	    isoTausMax = sortIsoTaus;
            tau1Index = tIndex;
            tau2Index = lIndex;
	    shift_t1_es = shift_tes_1;
	    shift_t2_es = shift_tes_2;
	    if (ptTau2>ptTau1) {
	      tau1Index = lIndex;
	      tau2Index = tIndex;
	      shift_t1_es = shift_tes_2;
	      shift_t2_es = shift_tes_1;
	    }
          }
        } // lepton loop
      } // tau loop
    
      //      std::cout << "OK1" << std::endl;

      if (tau1Index < 0) continue;
      if (tau2Index < 0) continue;
    
      FillTauTau(&analysisTree,otree,tau1Index,tau2Index,shift_t1_es,shift_t2_es);

      // rejecting events with pT(lep),pT(tau)<cut
      if (otree->pt_2<ptTauCut) continue;
      if (otree->pt_1<ptTauCut) continue;

      ////////////////////////////////////////////////////////////
      // Trigger matching
      ////////////////////////////////////////////////////////////
    
      vector<bool> isTriggerDiTau; isTriggerDiTau.clear();
      vector<unsigned int> nTriggerDiTau; nTriggerDiTau.clear();
      
      for (unsigned int iF=0; iF<filterDiTau.size(); ++iF) {
	unsigned int nFilter = 0;
	bool isFilterPresent = accessTriggerInfo(&analysisTree,TString(filterDiTau.at(iF)),nFilter);
	isTriggerDiTau.push_back(isFilterPresent);
	nTriggerDiTau.push_back(nFilter);
	cout << filterDiTau.at(iF) << " : " << isFilterPresent << " -> " << nFilter << std::endl;
      }
      //      cout << endl;

      bool tau1Matched = false;
      bool tau2Matched = false;
      for (unsigned int iF=0; iF<filterDiTau.size(); ++iF) {
	bool match1 = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isTriggerDiTau.at(iF),nTriggerDiTau.at(iF));
	tau1Matched = tau1Matched || match1;
	bool match2 = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isTriggerDiTau.at(iF),nTriggerDiTau.at(iF));
	tau2Matched = tau2Matched || match2;
      }      
    
      otree->trg_singlemuon = false;
      otree->trg_singleelectron = false;
      otree->singleLepTrigger = false;
      otree->trg_doubletau = tau1Matched && tau2Matched;
      otree->trg_mutaucross = false;
      otree->trg_mutaucross_mu = false;
      otree->trg_mutaucross_tau = false;
      otree->trg_etaucross = false;
      otree->trg_etaucross_e = false;
      otree->trg_etaucross_tau = false;
    
      cout << "pt_1 = " << otree->pt_1 
	   << "  eta_1 = " << otree->eta_1 
	   << "  phi_1 = " << otree->phi_1
	   << "  trg_match_1 = " << tau1Matched << endl;
      cout << "pt_2 = " << otree->pt_2 
	   << "  eta_2 = " << otree->eta_2 
	   << "  phi_2 = " << otree->phi_2
	   << "  trg_match_2 = " << tau2Matched << endl;
	
      ////////////////////////////////////////////////////////////
      // Filling variables
      ////////////////////////////////////////////////////////////

      //      if (!trg_doubletau) continue;

      //all criterua passed, we fill vertices here;	
      FillVertices(&analysisTree, otree, isData, tau1Index, tau2Index, ch);
    
      if (!isData)
        FillGenTree(&analysisTree, gentreeForGoodRecoEvtsOnly);       
    
      
      // initialize JER (including data and embedded) 
      otree->apply_recoil = ApplyRecoilCorrections;
      jets::initializeJER(&analysisTree);

      if (!isData && !isEmbedded) { // JER smearing
	jets::associateRecoAndGenJets(&analysisTree, resolution);
	jets::smear_jets(&analysisTree,resolution,resolution_sf,true);
      }

      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
  
      // setting weights to 1
      otree->trkeffweight = 1;
      otree->trigweight_1 = 1;
      otree->trigweight_2 = 1;
      otree->idisoweight_1 = 1;
      otree->idisoweight_antiiso_1 = 1;
      otree->idisoweight_2 = 1;
      otree->idisoweight_antiiso_2 = 1;
      otree->trigweight = 1;
      otree->trigweightSingle = 1;
      otree->trigweightExcl = 1;
      otree->effweight = 1;
      otree->effweightSingle = 1;
      otree->effweightExcl = 1;
      otree->puweight = 1; 
      otree->mcweight = 1;
      otree->weight = 1;
      otree->weightSingle = 1;
      otree->weightExcl = 1;
      otree->trigweight_l_lt = 1;
      otree->trigweight_t_lt = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM0_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM1_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM2_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM11_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM0_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM1_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM2_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_mt_MVADM11_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM0_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM1_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM2_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM10_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM11_13TeVUp = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM0_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM1_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM2_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM10_13TeVDown = 1;
      otree->weight_CMS_eff_Xtrigger_et_MVADM11_13TeVDown = 1;
			
      otree->weight_CMS_eff_t_pTlow_MVADM0_13TeVUp = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM1_13TeVUp = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM2_13TeVUp = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM10_13TeVUp = 1;
      otree->weight_CMS_eff_t_pTlow_MVADM11_13TeVUp = 1;
      otree->weight_CMS_eff_t_pThigh_MVADM0_13TeVUp = 1;
      otree->weight_CMS_eff_t_pThigh_MVADM1_13TeVUp = 1;
      otree->weight_CMS_eff_t_pThigh_MVADM2_13TeVUp = 1;
      otree->weight_CMS_eff_t_pThigh_MVADM10_13TeVUp = 1; 
      otree->weight_CMS_eff_t_pThigh_MVADM11_13TeVUp = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM0_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM1_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM2_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM10_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pTlow_MVADM11_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pThigh_MVADM0_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pThigh_MVADM1_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pThigh_MVADM2_13TeVDown = 1; 
      otree->weight_CMS_eff_t_pThigh_MVADM10_13TeVDown = 1;
      otree->weight_CMS_eff_t_pThigh_MVADM11_13TeVDown = 1;
			
      otree->weight_mufake_corr = 1; 
      otree->weight_CMS_mufake_mt_MVADM0_13TeVUp = 1; 
      otree->weight_CMS_mufake_mt_MVADM1_13TeVUp = 1; 
      otree->weight_CMS_mufake_mt_MVADM2_13TeVUp = 1; 
      otree->weight_CMS_mufake_mt_MVADM10_13TeVUp = 1;
      otree->weight_CMS_mufake_mt_MVADM11_13TeVUp = 1;
      otree->weight_CMS_mufake_mt_MVADM0_13TeVDown = 1; 
      otree->weight_CMS_mufake_mt_MVADM1_13TeVDown = 1; 
      otree->weight_CMS_mufake_mt_MVADM2_13TeVDown = 1; 
      otree->weight_CMS_mufake_mt_MVADM10_13TeVDown = 1;
      otree->weight_CMS_mufake_mt_MVADM11_13TeVDown = 1;
	
      if (ApplyPUweight) 
        otree->puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      if(!isData || isEmbedded){
        otree->mcweight = analysisTree.genweight;
        otree->gen_noutgoing = analysisTree.genparticles_noutgoing;
	if (isEmbedded&&otree->mcweight>1.0)
	  otree->mcweight = 0.0;
      }
      otree->weight = otree->puweight * otree->mcweight;

      // applying trigger/ID SFs       
      if ((!isData || isEmbedded) && ApplyLepSF) {

	// first tau
	w->var("t_pt")->setVal(otree->pt_1);
	w->var("t_eta")->setVal(otree->eta_1);
	w->var("t_phi")->setVal(otree->phi_1);
	w->var("t_dm")->setVal(analysisTree.tau_decayMode[tau1Index]);
	
	otree->idisoweight_1 = w->function("t_deeptauid_dm_medium")->getVal();
	if (isEmbedded) 
	  otree->idisoweight_1 = w->function("t_deeptauid_dm_embed_medium")->getVal();

	otree->trigweight_1 = w->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
	if (isEmbedded)
	  otree->trigweight_1 = w->function("t_trg_mediumDeepTau_ditau_embed_ratio")->getVal();

	// second tau
	w->var("t_pt")->setVal(otree->pt_2);
	w->var("t_eta")->setVal(otree->eta_2);
	w->var("t_phi")->setVal(otree->phi_2);
	w->var("t_dm")->setVal(analysisTree.tau_decayMode[tau2Index]);
	
	otree->idisoweight_2 = w->function("t_deeptauid_dm_medium")->getVal();
	if (isEmbedded) 
	  otree->idisoweight_2 = w->function("t_deeptauid_dm_embed_medium")->getVal();

	otree->trigweight_2 = w->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
	if (isEmbedded)
	  otree->trigweight_2 = w->function("t_trg_mediumDeepTau_ditau_embed_ratio")->getVal();

	if (otree->gen_match_1==1||otree->gen_match_1==3)
	  otree->idisoweight_1 = SF_eTau_hist->GetBinContent(SF_eTau_hist->GetXaxis()->FindBin(abs(otree->eta_1)));

	if (otree->gen_match_1==2||otree->gen_match_1==4)
	  otree->idisoweight_1 = SF_muTau_hist->GetBinContent(SF_muTau_hist->GetXaxis()->FindBin(abs(otree->eta_1)));

	if (otree->gen_match_2==1||otree->gen_match_2==3)
	  otree->idisoweight_2 = SF_eTau_hist->GetBinContent(SF_eTau_hist->GetXaxis()->FindBin(abs(otree->eta_2)));

	if (otree->gen_match_2==2||otree->gen_match_2==4)
	  otree->idisoweight_2 = SF_muTau_hist->GetBinContent(SF_muTau_hist->GetXaxis()->FindBin(abs(otree->eta_2)));

      }

      otree->trigweight = otree->trigweight_1 * otree->trigweight_2;

      otree->effweight = otree->idisoweight_1 * otree->idisoweight_2 * otree->trigweight;

      otree->weight *= otree->effweight;
      
      //Theory uncertainties 
      
      otree->weight_CMS_scale_gg_13TeVUp   = 1.;
      otree->weight_CMS_scale_gg_13TeVDown = 1.;

      otree->weight_CMS_PS_ISR_ggH_13TeVUp   = 1.;
      otree->weight_CMS_PS_ISR_ggH_13TeVDown = 1.;
      otree->weight_CMS_PS_FSR_ggH_13TeVUp   = 1.;
      otree->weight_CMS_PS_FSR_ggH_13TeVDown = 1.;

      if(isHiggs){

	otree->weight_CMS_scale_gg_13TeVUp   = analysisTree.weightScale4;
	otree->weight_CMS_scale_gg_13TeVDown = analysisTree.weightScale8;

	otree->weight_CMS_PS_ISR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[6];
	otree->weight_CMS_PS_ISR_ggH_13TeVDown = analysisTree.gen_pythiaweights[8];
	otree->weight_CMS_PS_FSR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[7];
	otree->weight_CMS_PS_FSR_ggH_13TeVDown = analysisTree.gen_pythiaweights[9];
      }

      //Prefiring weights
      otree->prefiringweight = 1.0;
      otree->prefiringweightUp = 1.0;
      otree->prefiringweightDown = 1.0;
      if (era<2018) {
	if (!isData) {
	  otree->prefiringweight     = analysisTree.prefiringweight;
	  otree->prefiringweightUp   = analysisTree.prefiringweight;
	  otree->prefiringweightDown = analysisTree.prefiringweight;
	  if (otree->prefiringweight>0.1) {
	    otree->prefiringweightUp   = analysisTree.prefiringweightup/analysisTree.prefiringweight;
	    otree->prefiringweightDown = analysisTree.prefiringweightdown/analysisTree.prefiringweight;
	  }
	}
      }
      otree->weight *= otree->prefiringweight;

      // embedded weight
      otree->embweight = 1;
      if (isEmbedded) {
      	otree->embweight = getEmbeddedWeight(&analysisTree, w);
	if (otree->embweight>10.0) {
	  cout << "warning : embedding weight = " << otree->embweight << endl;
	  otree->embweight = 0.0;
	}
      }
      otree->weight *= otree->embweight;



      ////////////////////////////////////////////////////////////
      // Z pt weight
      ////////////////////////////////////////////////////////////      
      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);
      otree->zptweight = 1.;
      if (!isData && isDY){
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
      otree->weight *= otree->zptweight;
      
      ////////////////////////////////////////////////////////////
      // Top pt weight
      ////////////////////////////////////////////////////////////
      otree->topptweight = 1.;
      if(!isData || isEmbedded){
        float a_topPtWeight = cfg.get<float>("a_topPtWeight");
        float b_topPtWeight = cfg.get<float>("b_topPtWeight");
        float c_topPtWeight = cfg.get<float>("c_topPtWeight");
        float max_pt_topPtWeight = cfg.get<float>("max_pt_topPtWeight");
        otree->topptweight = genTools::topPtWeight_Run2(analysisTree, a_topPtWeight, b_topPtWeight, c_topPtWeight, max_pt_topPtWeight);
        // otree->topptweight = genTools::topPtWeight(analysisTree, 1); // 1 is for Run1 - use this reweighting as recommended by HTT 17
      }
      otree->weight *= otree->topptweight;

      ////////////////////////////////////////////////////////////
      // MET and Recoil Corrections
      ////////////////////////////////////////////////////////////     
      otree->met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
      otree->metphi = TMath::ATan2(analysisTree.pfmetcorr_ey,analysisTree.pfmetcorr_ex);
      otree->metcov00 = analysisTree.pfmetcorr_sigxx;
      otree->metcov01 = analysisTree.pfmetcorr_sigxy;
      otree->metcov10 = analysisTree.pfmetcorr_sigyx;
      otree->metcov11 = analysisTree.pfmetcorr_sigyy;

      otree->puppimet = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex + analysisTree.puppimet_ey*analysisTree.puppimet_ey);
      otree->puppimetphi = TMath::ATan2(analysisTree.puppimet_ey,analysisTree.puppimet_ex);
      otree->puppimetcov00 = analysisTree.puppimet_sigxx;
      otree->puppimetcov01 = analysisTree.puppimet_sigxy;
      otree->puppimetcov10 = analysisTree.puppimet_sigyx;
      otree->puppimetcov11 = analysisTree.puppimet_sigyy;

      otree->met_uncorr = otree->met;
      otree->metphi_uncorr = otree->metphi;
      if (usePuppiMET) {
	otree->met_uncorr = otree->puppimet;
	otree->metphi_uncorr = otree->puppimetphi;
      }
      otree->njetshad = otree->njets;
      if (isWJets) otree->njetshad += 1;

      if(ApplyRecoilCorrections){        
      	genV = genTools::genV(analysisTree);
      	genL = genTools::genL(analysisTree);
        genTools::KITRecoilCorrections( recoilCorrector, ApplyRecoilCorrections, // pass the value != 0 to apply corrections
          otree->met_uncorr, otree->metphi_uncorr,
          genV.Px(), genV.Py(),
          genL.Px(), genL.Py(),
          otree->njetshad,
          otree->met_rcmr, otree->metphi_rcmr
        );
        
        // overwriting with recoil-corrected values 
	if (usePuppiMET) {
	  otree->puppimet = otree->met_rcmr;
	  otree->puppimetphi = otree->metphi_rcmr;
	}   
	else {
	  otree->met = otree->met_rcmr;
          otree->metphi = otree->metphi_rcmr;
	}
      }
      
      //ditau sytem
      TLorentzVector tau1LV; tau1LV.SetPtEtaPhiM(otree->pt_1,
						 otree->eta_1,
						 otree->phi_1,
						 otree->m_1);

      TLorentzVector tau1_uncorrLV; tau1_uncorrLV.SetPtEtaPhiM(otree->pt_uncorr_1,
							       otree->eta_1,
							       otree->phi_1,
							       otree->m_1);

      TLorentzVector tau2LV; tau2LV.SetPtEtaPhiM(otree->pt_2,
						 otree->eta_2,
						 otree->phi_2,
						 otree->m_2);

      TLorentzVector tau2_uncorrLV; tau2_uncorrLV.SetPtEtaPhiM(otree->pt_uncorr_2,
							       otree->eta_2,
							       otree->phi_2,
							       otree->m_2);

    
      // MET
      TLorentzVector metLV, puppimetLV; 
      TLorentzVector puppimetUnclLV_Up, puppimetUnclLV_Down;
      TLorentzVector metUnclLV_Up, metUnclLV_Down;

      float met_x = otree->met*TMath::Cos(otree->metphi);
      float met_y = otree->met*TMath::Sin(otree->metphi);

      float puppimet_x = otree->puppimet*TMath::Cos(otree->puppimetphi);
      float puppimet_y = otree->puppimet*TMath::Sin(otree->puppimetphi);

      float puppimetUncl_Up = TMath::Sqrt(analysisTree.puppimet_ex_UnclusteredEnUp*analysisTree.puppimet_ex_UnclusteredEnUp+analysisTree.puppimet_ey_UnclusteredEnUp*analysisTree.puppimet_ey_UnclusteredEnUp);
      float puppimetUncl_Down = TMath::Sqrt(analysisTree.puppimet_ex_UnclusteredEnDown*analysisTree.puppimet_ex_UnclusteredEnDown+analysisTree.puppimet_ey_UnclusteredEnDown*analysisTree.puppimet_ey_UnclusteredEnDown);

      float metUncl_Up = TMath::Sqrt(analysisTree.pfmetcorr_ex_UnclusteredEnUp*analysisTree.pfmetcorr_ex_UnclusteredEnUp+analysisTree.pfmetcorr_ey_UnclusteredEnUp*analysisTree.pfmetcorr_ey_UnclusteredEnUp);
      float metUncl_Down = TMath::Sqrt(analysisTree.pfmetcorr_ex_UnclusteredEnDown*analysisTree.pfmetcorr_ex_UnclusteredEnDown+analysisTree.pfmetcorr_ey_UnclusteredEnDown*analysisTree.pfmetcorr_ey_UnclusteredEnDown);
      
      metLV.SetXYZT(met_x,met_y,0.,otree->met);
      puppimetLV.SetXYZT(puppimet_x,puppimet_y,0.,otree->puppimet);
      puppimetUnclLV_Up.SetXYZT(analysisTree.puppimet_ex_UnclusteredEnUp,
				analysisTree.puppimet_ey_UnclusteredEnUp,
				0.,puppimetUncl_Up);
      puppimetUnclLV_Down.SetXYZT(analysisTree.puppimet_ex_UnclusteredEnDown,
				  analysisTree.puppimet_ey_UnclusteredEnDown,
				  0.,puppimetUncl_Down);
      metUnclLV_Up.SetXYZT(analysisTree.pfmetcorr_ex_UnclusteredEnUp,
			   analysisTree.pfmetcorr_ey_UnclusteredEnUp,
			   0.,metUncl_Up);
      metUnclLV_Down.SetXYZT(analysisTree.pfmetcorr_ex_UnclusteredEnDown,
			     analysisTree.pfmetcorr_ey_UnclusteredEnDown,
			     0.,metUncl_Down);
        
      ////////////////////////////////////////////////////////////
      // Tau ES shift -> propagate to MET
      ////////////////////////////////////////////////////////////      
      if (!isData||isEmbedded) {
	TLorentzVector corr_LV = tau1_uncorrLV + tau2_uncorrLV - tau1LV - tau2LV;
	puppimetLV += corr_LV;
	puppimetUnclLV_Up += corr_LV;
	puppimetUnclLV_Down += corr_LV;
	otree->puppimet = puppimetLV.Pt();
	otree->puppimetphi = puppimetLV.Phi();
	otree->puppimet_ex_UnclusteredEnUp = puppimetUnclLV_Up.Px();
	otree->puppimet_ey_UnclusteredEnUp = puppimetUnclLV_Up.Py();
	otree->puppimet_ex_UnclusteredEnDown = puppimetUnclLV_Down.Px();
	otree->puppimet_ey_UnclusteredEnDown = puppimetUnclLV_Down.Py();
	
	metLV += corr_LV;
	metUnclLV_Up += corr_LV;
	metUnclLV_Down += corr_LV;
	otree->met = metLV.Pt();
	otree->metphi = metLV.Phi();
	otree->met_ex_UnclusteredEnUp = metUnclLV_Up.Px();
	otree->met_ey_UnclusteredEnUp = metUnclLV_Up.Py();
	otree->met_ex_UnclusteredEnDown = metUnclLV_Down.Px();
	otree->met_ey_UnclusteredEnDown = metUnclLV_Down.Py();
	
      }
          
      TLorentzVector dileptonLV = tau1LV + tau2LV;
      otree->m_vis = dileptonLV.M();
      otree->pt_tt = (dileptonLV+metLV).Pt();   
      if (usePuppiMET)
	otree->pt_tt = (dileptonLV+puppimetLV).Pt();
    
      // mt TOT
      TLorentzVector metxLV = metLV;
      if (usePuppiMET) metxLV = puppimetLV;
      float mtTOT = 2*(otree->pt_1)*metxLV.Pt()*(1-cos(DeltaPhi(tau1LV,metxLV)));
      mtTOT += 2*(otree->pt_2)*metxLV.Pt()*(1-cos(DeltaPhi(tau2LV,metxLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(tau1LV,tau2LV))); 
      otree->mt_tot = TMath::Sqrt(mtTOT);
    
      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;
    
      //extra lepton vetos
      otree->extraelec_veto = extra_electron_veto_tt(&cfg, &analysisTree, otree, era, isEmbedded);
      otree->extramuon_veto = extra_muon_veto_tt(&cfg, &analysisTree, otree, isData);
    
      otree->mt_1 = mT(tau1LV, metLV);
      otree->mt_2 = mT(tau2LV, metLV);
      otree->puppimt_1 = mT(tau1LV, puppimetLV);
      otree->puppimt_2 = mT(tau2LV, puppimetLV);
    
      bool isSRevent = true; //boolean used to compute SVFit variables only on SR events, it is set to true when running Synchronization to run SVFit on all events
      if(!Synch){
	isSRevent = (otree->trg_doubletau>0.5 && otree->byVVVLooseDeepTau2017v2p1VSjet_1>0.5 && otree->byVVVLooseDeepTau2017v2p1VSjet_2>0.5);
      }
      /*
	if (!isSRevent) { 
	cout << "                                        " << endl;
	cout << "========================================" << endl;
	cout << "        Event is not selected           " << endl;
	cout << "========================================" << endl;
	cout << "                                        " << endl;
	}
      */
      double met_var_qcd_1 = puppimetLV.Pt()*cos(DeltaPhi(tau1LV,puppimetLV))/tau1LV.Pt();
      auto args = std::vector<double>{otree->pt_1,
				      static_cast<double>(otree->tau_decay_mode_1),
				      static_cast<double>(otree->njets),
				      static_cast<double>(otree->os),
				      met_var_qcd_1,
				      otree->dr_tt};
      otree->ff_nom = fns_["ff_tt_medium_dmbins"]->eval(args.data());
      otree->ff_nom_sys = 0.15;

      cout << "Trigger weight = " << otree->trigweight << endl;
      cout << "id/Iso1 weight = " << otree->idisoweight_1 << endl;
      cout << "id/Iso2 weight = " << otree->idisoweight_2 << endl;
      cout << "mc weight      = " << otree->mcweight << endl;
      cout << "pu weight      = " << otree->puweight << endl;
      cout << "embed weight   = " << otree->embweight << endl;
      cout << "prefire weight = " << otree->prefiringweight << endl;
      cout << "total weight   = " << otree->weight << endl;
      cout << "trg_doubletau  = " << otree->trg_doubletau << endl;
      cout << "Fake factor = " << otree->ff_nom << endl;
      cout << endl;

      // initialize svfit and fastMTT variables
      otree->m_sv   = -10;
      otree->pt_sv  = -10;
      otree->eta_sv = -10;
      otree->phi_sv = -10;
      otree->met_sv = -10;
      otree->mt_sv = -10;
      otree->m_fast = -10;
      otree->mt_fast = -10;
      otree->pt_fast = -10;
      otree->phi_fast = -10;
      otree->eta_fast = -10;
      if (!isSRevent && ApplySystShift) continue;
      if ( (ApplySVFit||ApplyFastMTT) && isSRevent ) svfit_variables(ch, &analysisTree, otree, &cfg, inputFile_visPtResolution);
        
      // evaluate systematics for MC 
      if( !isData && !isEmbedded && ApplySystShift){
	if (!isDY && !isWJets && !isHiggs) {
	  btagSys->Eval();
	}
	for(unsigned int i = 0; i < jetEnergyScaleSys.size(); i++) {
	  //	  cout << endl;
	  //	  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	  //	  cout << endl;	  
	  (jetEnergyScaleSys.at(i))->Eval(); 
	  //	  cout << endl;
	  //	  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl; 
	  //	  cout << endl;
	}
	if (usePuppiMET) {
	  for(unsigned int i = 0; i < puppiMetSys.size(); ++i)
	    (puppiMetSys.at(i))->Eval();
	}
	else {
	  for(unsigned int i = 0; i < metSys.size(); ++i) 
	    (metSys.at(i))->Eval();
	}
      }
			
      if ((!isData||isEmbedded) && ApplySystShift) {

	tauOneProngScaleSys->Eval(utils::MUTAU);
	tauOneProngOnePi0ScaleSys->Eval(utils::MUTAU);
	tauThreeProngScaleSys->Eval(utils::MUTAU);
	tauThreeProngOnePi0ScaleSys->Eval(utils::MUTAU);

      }
      
      selEvents++;
      
      
      //      cout << "Once again puppimet central : " << otree->puppimet << endl;
      //Merijn 2019 1 10: perhaps this should be called before moving to next event..
      otree->Fill();
    } // event loop
    
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  } // file loop
   

  std::cout << std::endl;
  std::cout << "Total number of input events    = " << int(inputEventsH->GetEntries()) << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();
  // delete systematics objects

  if(tauScaleSys != 0){
    tauScaleSys->Write("",TObject::kOverwrite);
    delete tauScaleSys;
  }

  if(tauOneProngScaleSys != 0){
    tauOneProngScaleSys->Write("",TObject::kOverwrite);
    delete tauOneProngScaleSys;
  }

  if(tauOneProngOnePi0ScaleSys != 0){
    tauOneProngOnePi0ScaleSys->Write("",TObject::kOverwrite);
    delete tauOneProngOnePi0ScaleSys;
  }

  if(tauThreeProngScaleSys != 0){
    tauThreeProngScaleSys->Write("",TObject::kOverwrite);
    delete tauThreeProngScaleSys;
  }

  if(tauThreeProngOnePi0ScaleSys != 0){
    tauThreeProngOnePi0ScaleSys->Write("",TObject::kOverwrite);
    delete tauThreeProngOnePi0ScaleSys;
  }

  if (btagSys != 0) {
    btagSys->Write("",TObject::kOverwrite);
    delete btagSys;
  }

  if(jetEnergyScaleSys.size() > 0){
    for (unsigned int i = 0; i < jetEnergyScaleSys.size(); i++){
      (jetEnergyScaleSys.at(i))->Write("",TObject::kOverwrite);
      delete jetEnergyScaleSys.at(i);
    }
  }

  if (metSys.size() > 0){
    for (unsigned int i = 0; i < metSys.size(); i++ ) {
      (metSys.at(i))->Write("",TObject::kOverwrite);
      delete metSys.at(i);
    }
  }

  if (puppiMetSys.size() > 0){
    for (unsigned int i = 0; i < puppiMetSys.size(); i++ ) {
      (puppiMetSys.at(i))->Write("",TObject::kOverwrite);
      delete puppiMetSys.at(i);
    }
  }

  file->Close();
  delete file;

}


////FILLING FUNCTIONS//////

void FillVertices(const AC1B *analysisTree, Synch17Tree *otree, const bool isData, int leptonIndex, int tauIndex, TString channel){

  otree->RecoVertexX = analysisTree->primvertex_x;
  otree->RecoVertexY = analysisTree->primvertex_y;
  otree->RecoVertexZ = analysisTree->primvertex_z;
  
  bool is_refitted_PV_with_BS = true;
  TVector3 vertex_refitted_BS = get_refitted_PV_with_BS(analysisTree, leptonIndex, tauIndex, channel, is_refitted_PV_with_BS);
  otree->pvx = vertex_refitted_BS.X();
  otree->pvy = vertex_refitted_BS.Y();
  otree->pvz = vertex_refitted_BS.Z();
  otree->is_refitted_PV_with_BS = is_refitted_PV_with_BS;

  otree->pvx_bs = analysisTree->primvertexwithbs_x;
  otree->pvy_bs = analysisTree->primvertexwithbs_y;
  otree->pvz_bs = analysisTree->primvertexwithbs_z;

  otree->tau_SV_x_2 = analysisTree->tau_SV_x[tauIndex];
  otree->tau_SV_y_2 = analysisTree->tau_SV_y[tauIndex];
  otree->tau_SV_z_2 = analysisTree->tau_SV_z[tauIndex];

  otree->tau_SV_covxx_2 = analysisTree->tau_SV_cov[tauIndex][0];
  otree->tau_SV_covyx_2 = analysisTree->tau_SV_cov[tauIndex][1];
  otree->tau_SV_covzx_2 = analysisTree->tau_SV_cov[tauIndex][2];
  otree->tau_SV_covyy_2 = analysisTree->tau_SV_cov[tauIndex][3];
  otree->tau_SV_covzy_2 = analysisTree->tau_SV_cov[tauIndex][4];
  otree->tau_SV_covzz_2 = analysisTree->tau_SV_cov[tauIndex][5];

  otree->SVminusRefitV_x=otree->tau_SV_x_2-vertex_refitted_BS.X();
  otree->SVminusRefitV_y=otree->tau_SV_y_2-vertex_refitted_BS.Y();
  otree->SVminusRefitV_z=otree->tau_SV_z_2-vertex_refitted_BS.Z();

  otree->SVminusPV_x=otree->tau_SV_x_2-otree->RecoVertexX;
  otree->SVminusPV_y=otree->tau_SV_y_2-otree->RecoVertexY;
  otree->SVminusPV_z=otree->tau_SV_z_2-otree->RecoVertexZ;


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

float getEmbeddedWeight(const AC1B *analysisTree, RooWorkspace * wEm) {

  std::vector<TLorentzVector> taus; taus.clear();
  float emWeight = 1;
  for (unsigned int igentau = 0; igentau < analysisTree->gentau_count; ++igentau) {
    TLorentzVector tauLV; tauLV.SetXYZT(analysisTree->gentau_px[igentau], 
					analysisTree->gentau_py[igentau],
					analysisTree->gentau_pz[igentau],
					analysisTree->gentau_e[igentau]);
    if (analysisTree->gentau_isPrompt[igentau]&&analysisTree->gentau_isFirstCopy[igentau]) {
      taus.push_back(tauLV);
    }
  }

  //  std::cout << "n taus = " << taus.size() << "  :  wEm = " << wEm << std::endl;

  if (taus.size() == 2) {
    double gt1_pt  = taus[0].Pt();
    double gt1_eta = taus[0].Eta();
    double gt2_pt  = taus[1].Pt();
    double gt2_eta = taus[1].Eta();
    wEm->var("gt_pt")->setVal(gt1_pt);
    wEm->var("gt_eta")->setVal(gt1_eta);
    double id1_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt_pt")->setVal(gt2_pt);
    wEm->var("gt_eta")->setVal(gt2_eta);
    double id2_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt1_pt")->setVal(gt1_pt);
    wEm->var("gt2_pt")->setVal(gt2_pt);
    wEm->var("gt1_eta")->setVal(gt1_eta);
    wEm->var("gt2_eta")->setVal(gt2_eta);
    double trg_emb = wEm->function("m_sel_trg_ic_ratio")->getVal();
    emWeight = id1_embed * id2_embed * trg_emb;
  }

  return emWeight;

}

float shift_tauES(const AC1B * analysisTree, 
		  unsigned int itau,
		  float shift_tes_1prong,
		  float shift_tes_1p1p0,
		  float shift_tes_3prong,
		  float shift_tes_3prong1p0,
		  float shift_tes_lepfake_1prong_barrel,
		  float shift_tes_lepfake_1p1p0_barrel,
		  float shift_tes_lepfake_1prong_endcap,
		  float shift_tes_lepfake_1p1p0_endcap
		  ) {
  float shift_tes = 0.0;
  int gen_match = analysisTree->tau_genmatch[itau];
  int decay_mode = analysisTree->tau_decayMode[itau];
  if (gen_match >= 5){
    if (decay_mode == 0) shift_tes = shift_tes_1prong; 
    else if (decay_mode >= 1 && decay_mode <= 3) shift_tes = shift_tes_1p1p0; 
    else if (decay_mode == 10) shift_tes = shift_tes_3prong;
    else if (decay_mode == 12) shift_tes = shift_tes_3prong1p0;
  }
  else if (gen_match==1 || gen_match==3) {
    if (decay_mode == 0){
      if (fabs(analysisTree->tau_eta[itau]) < 1.5) // barrel definition taken from https://github.com/cms-tau-pog/TauIDSFs/blob/b4963b627df0ba85bce4aeaeacf401bc35686246/python/TauIDSFTool.py#L231
	shift_tes = shift_tes_lepfake_1prong_barrel; 
      else
	shift_tes = shift_tes_lepfake_1prong_endcap; 
    } 	 
    else if (decay_mode >= 1 && decay_mode <= 3){
      if (fabs(analysisTree->tau_eta[itau]) < 1.5)
	shift_tes = shift_tes_lepfake_1p1p0_barrel; 
      else
	shift_tes = shift_tes_lepfake_1p1p0_endcap; 
    } 
  }
  return shift_tes;
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

//fill the otree with the tau variables 
void FillTauTau(const AC1B *analysisTree, Synch17Tree *otree, int tau1Index, int tau2Index, float shift_ts_1, float shift_ts_2){

  otree->pt_1 = (1.0+shift_ts_1)*analysisTree->tau_pt[tau1Index];
  otree->pt_uncorr_1 = analysisTree->tau_pt[tau1Index];
  otree->eta_1 = analysisTree->tau_eta[tau1Index];
  otree->phi_1 = analysisTree->tau_phi[tau1Index];
  otree->q_1 = analysisTree->tau_charge[tau1Index];
  otree->gen_match_1 = analysisTree->tau_genmatch[tau1Index];
  otree->d0_1 = analysisTree->tau_leadchargedhadrcand_dxy[tau1Index];
  otree->dZ_1 = analysisTree->tau_leadchargedhadrcand_dz[tau1Index];      
  otree->iso_1 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau1Index];
  otree->tau_decay_mode_1 = analysisTree->tau_decayMode[tau1Index];
  if (otree->tau_decay_mode_1==0)
    otree->m_1 = analysisTree->tau_mass[tau1Index];
  else
    otree->m_1 = (1.0+shift_ts_1)*analysisTree->tau_mass[tau1Index];
  otree->dm_1 = analysisTree->tau_decayMode[tau1Index];
  otree->dmMVA_1 = analysisTree->tau_MVADM2017v1[tau1Index];

  otree->deepTauVsEleRaw_1                = analysisTree->tau_byDeepTau2017v2p1VSeraw[tau1Index];
  otree->deepTauVsJetRaw_1                = analysisTree->tau_byDeepTau2017v2p1VSjetraw[tau1Index];
  otree->deepTauVsMuRaw_1                 = analysisTree->tau_byDeepTau2017v2p1VSmuraw[tau1Index];
  otree->byLooseDeepTau2017v2p1VSe_1      = analysisTree->tau_byLooseDeepTau2017v2p1VSe[tau1Index];
  otree->byLooseDeepTau2017v2p1VSjet_1    = analysisTree->tau_byLooseDeepTau2017v2p1VSjet[tau1Index];
  otree->byLooseDeepTau2017v2p1VSmu_1     = analysisTree->tau_byLooseDeepTau2017v2p1VSmu[tau1Index];
  otree->byMediumDeepTau2017v2p1VSe_1     = analysisTree->tau_byMediumDeepTau2017v2p1VSe[tau1Index];
  otree->byMediumDeepTau2017v2p1VSjet_1   = analysisTree->tau_byMediumDeepTau2017v2p1VSjet[tau1Index];
  otree->byMediumDeepTau2017v2p1VSmu_1    = analysisTree->tau_byMediumDeepTau2017v2p1VSmu[tau1Index];
  otree->byTightDeepTau2017v2p1VSe_1      = analysisTree->tau_byTightDeepTau2017v2p1VSe[tau1Index];
  otree->byTightDeepTau2017v2p1VSjet_1    = analysisTree->tau_byTightDeepTau2017v2p1VSjet[tau1Index];
  otree->byTightDeepTau2017v2p1VSmu_1     = analysisTree->tau_byTightDeepTau2017v2p1VSmu[tau1Index];
  otree->byVLooseDeepTau2017v2p1VSe_1     = analysisTree->tau_byVLooseDeepTau2017v2p1VSe[tau1Index];
  otree->byVLooseDeepTau2017v2p1VSjet_1   = analysisTree->tau_byVLooseDeepTau2017v2p1VSjet[tau1Index];
  otree->byVLooseDeepTau2017v2p1VSmu_1    = analysisTree->tau_byVLooseDeepTau2017v2p1VSmu[tau1Index];
  otree->byVTightDeepTau2017v2p1VSe_1     = analysisTree->tau_byVTightDeepTau2017v2p1VSe[tau1Index];
  otree->byVTightDeepTau2017v2p1VSjet_1   = analysisTree->tau_byVTightDeepTau2017v2p1VSjet[tau1Index];
  otree->byVVLooseDeepTau2017v2p1VSe_1    = analysisTree->tau_byVVLooseDeepTau2017v2p1VSe[tau1Index];
  otree->byVVLooseDeepTau2017v2p1VSjet_1  = analysisTree->tau_byVVLooseDeepTau2017v2p1VSjet[tau1Index];
  otree->byVVTightDeepTau2017v2p1VSe_1    = analysisTree->tau_byVVTightDeepTau2017v2p1VSe[tau1Index];
  otree->byVVTightDeepTau2017v2p1VSjet_1  = analysisTree->tau_byVVTightDeepTau2017v2p1VSjet[tau1Index];
  otree->byVVVLooseDeepTau2017v2p1VSe_1   = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSe[tau1Index];
  otree->byVVVLooseDeepTau2017v2p1VSjet_1 = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSjet[tau1Index];

  otree->pt_2 = (1.0+shift_ts_2)*analysisTree->tau_pt[tau2Index];
  otree->pt_uncorr_2 = analysisTree->tau_pt[tau2Index];
  otree->eta_2 = analysisTree->tau_eta[tau2Index];
  otree->phi_2 = analysisTree->tau_phi[tau2Index];
  otree->q_2 = analysisTree->tau_charge[tau2Index];
  otree->gen_match_2 = analysisTree->tau_genmatch[tau2Index];
  otree->d0_2 = analysisTree->tau_leadchargedhadrcand_dxy[tau2Index];
  otree->dZ_2 = analysisTree->tau_leadchargedhadrcand_dz[tau2Index];      
  otree->iso_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau2Index];
  otree->tau_decay_mode_2 = analysisTree->tau_decayMode[tau2Index];
  if (otree->tau_decay_mode_2==0)
    otree->m_2 = analysisTree->tau_mass[tau2Index];
  else
    otree->m_2 = (1.0+shift_ts_2)*analysisTree->tau_mass[tau2Index];
  otree->dm_2 = analysisTree->tau_decayMode[tau2Index];
  otree->dmMVA_2 = analysisTree->tau_MVADM2017v1[tau2Index];

  otree->deepTauVsEleRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSeraw[tau2Index];
  otree->deepTauVsJetRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSjetraw[tau2Index];
  otree->deepTauVsMuRaw_2                 = analysisTree->tau_byDeepTau2017v2p1VSmuraw[tau2Index];
  otree->byLooseDeepTau2017v2p1VSe_2      = analysisTree->tau_byLooseDeepTau2017v2p1VSe[tau2Index];
  otree->byLooseDeepTau2017v2p1VSjet_2    = analysisTree->tau_byLooseDeepTau2017v2p1VSjet[tau2Index];
  otree->byLooseDeepTau2017v2p1VSmu_2     = analysisTree->tau_byLooseDeepTau2017v2p1VSmu[tau2Index];
  otree->byMediumDeepTau2017v2p1VSe_2     = analysisTree->tau_byMediumDeepTau2017v2p1VSe[tau2Index];
  otree->byMediumDeepTau2017v2p1VSjet_2   = analysisTree->tau_byMediumDeepTau2017v2p1VSjet[tau2Index];
  otree->byMediumDeepTau2017v2p1VSmu_2    = analysisTree->tau_byMediumDeepTau2017v2p1VSmu[tau2Index];
  otree->byTightDeepTau2017v2p1VSe_2      = analysisTree->tau_byTightDeepTau2017v2p1VSe[tau2Index];
  otree->byTightDeepTau2017v2p1VSjet_2    = analysisTree->tau_byTightDeepTau2017v2p1VSjet[tau2Index];
  otree->byTightDeepTau2017v2p1VSmu_2     = analysisTree->tau_byTightDeepTau2017v2p1VSmu[tau2Index];
  otree->byVLooseDeepTau2017v2p1VSe_2     = analysisTree->tau_byVLooseDeepTau2017v2p1VSe[tau2Index];
  otree->byVLooseDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVLooseDeepTau2017v2p1VSjet[tau2Index];
  otree->byVLooseDeepTau2017v2p1VSmu_2    = analysisTree->tau_byVLooseDeepTau2017v2p1VSmu[tau2Index];
  otree->byVTightDeepTau2017v2p1VSe_2     = analysisTree->tau_byVTightDeepTau2017v2p1VSe[tau2Index];
  otree->byVTightDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVTightDeepTau2017v2p1VSjet[tau2Index];
  otree->byVVLooseDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVLooseDeepTau2017v2p1VSe[tau2Index];
  otree->byVVLooseDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVLooseDeepTau2017v2p1VSjet[tau2Index];
  otree->byVVTightDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVTightDeepTau2017v2p1VSe[tau2Index];
  otree->byVVTightDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVTightDeepTau2017v2p1VSjet[tau2Index];
  otree->byVVVLooseDeepTau2017v2p1VSe_2   = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSe[tau2Index];
  otree->byVVVLooseDeepTau2017v2p1VSjet_2 = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSjet[tau2Index];

  otree->dr_tt = deltaR(otree->eta_1,otree->phi_1,otree->eta_2,otree->phi_2);

}
