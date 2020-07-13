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

void initializeCPvar(Synch17Tree *otree);
void initializeGenTree(Synch17GenTree *gentree);
void FillMuTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex, float dRiso);
void FillETau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso, int era, bool isEmbedded);
void FillTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex);
void FillVertices(const AC1B * analysisTree,Synch17Tree *otree, const bool isData, int leptonIndex, int tauIndex, TString channel);
void FillGenTree(const AC1B * analysisTree, Synch17GenTree *gentree);
float getEmbeddedWeight(const AC1B * analysisTree, RooWorkspace* WS);
float shift_tauES(const AC1B * analysisTree, unsigned int itau, 
									float shift_tes_1prong, float shift_tes_1p1p0, float shift_tes_3prong,
								  float shift_tes_lepfake_1prong_barrel, float shift_tes_lepfake_1p1p0_barrel,
								  float shift_tes_lepfake_1prong_endcap, float shift_tes_lepfake_1p1p0_endcap
								); 

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
    std::cout << "SynchNTupleProducer_Run2 NameOfTheConfigurationFile FileList Channel" << std::endl;
    std::cout << "example: SynchNTupleProducer_Run2 analysisMacroSynch_lept_mt_DATA_SingleMuon.conf DATA_SingleMuon mt" << std::endl;
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
  else if (ch == "et") lep = "Electron";

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase) + "/src/" + TString(json_name)).Data(), json);
  }

  const int era = cfg.get<int>("era");
  const bool Synch = cfg.get<bool>("Synch"); 
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyFastMTT     = cfg.get<bool>("ApplyFastMTT");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters  = cfg.get<bool>("ApplyMetFilters");
  const bool usePuppiMET      = cfg.get<bool>("UsePuppiMET");
  const bool ApplyIpCorrection = cfg.get<bool>("ApplyIpCorrection");
  const bool ApplyBTagCP5Correction = cfg.get<bool>("ApplyBTagCP5Correction");

  // JER
  //  const string jer_resolution = cfg.get<string>("JER_Resolution");
  //  const string jer_scalefactor = cfg.get<string>("JER_ScaleFactor");

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");

  const string ipCorrFileNameLepton   = cfg.get<string>("IpCorrFileNameLepton");
  const string ipCorrFileNameLeptonBS = cfg.get<string>("IpCorrFileNameLeptonPVBS");
  const string ipCorrFileNamePion     = cfg.get<string>("IpCorrFileNamePion");
  const string ipCorrFileNamePionBS   = cfg.get<string>("IpCorrFileNamePionPVBS");

  TString IpCorrFileNameLepton(ipCorrFileNameLepton);
  TString IpCorrFileNameLeptonBS(ipCorrFileNameLeptonBS);
  TString IpCorrFileNamePion(ipCorrFileNamePion);
  TString IpCorrFileNamePionBS(ipCorrFileNamePionBS);


  // tau trigger efficiency
  std::string channel;
  if (ch == "mt") channel = "mutau"; 
  if (ch == "et") channel = "etau";

  TauTriggerSFs2017 *tauTriggerSF = new TauTriggerSFs2017(cmsswBase + "/src/TauAnalysisTools/TauTriggerSFs/data/tauTriggerEfficiencies" + to_string(era) + ".root", channel, to_string(era), "tight", "MVAv2");
  std::string year_label;
  if (era == 2016) year_label = "2016Legacy";
  else if (era == 2017) year_label = "2017ReReco";
  else if (era == 2018) year_label = "2018ReReco";	
  else {std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n'; exit(-1);}
  TauIDSFTool *tauIDSF_medium = new TauIDSFTool(year_label, "DeepTau2017v2p1VSjet", "Medium", false);

  IpCorrection *ipLepton   = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameLepton);
  IpCorrection *ipLeptonBS = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameLeptonBS);  
  IpCorrection *ipPion     = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNamePion);
  IpCorrection *ipPionBS   = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNamePionBS);  

  std::map<TString,IpCorrection*> ipCorrectors = {
    {"ipTau1"   ,ipLepton},
    {"ipTau1BS" ,ipLeptonBS},
    {"ipTau2"     ,ipPion},
    {"ipTau2BS"   ,ipPionBS}
  };
  

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
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Autumn18_V7b_MC_SF_AK8PFchs.txt"));    
  }

  JME::JetResolution resolution = *m_resolution_from_file;
  JME::JetResolutionScaleFactor resolution_sf = *m_scale_factor_from_file;

  cout<<"using "<<BTagAlgorithm<<endl;
  BTagCalibration calib(BTagAlgorithm, BtagSfFile);
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM, "central",{"up","down"});
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM, "central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM, "central");
  if(ApplyBTagScaling){
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

  TFile * ff_file = TFile::Open(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_"+TString(year_label)+".root");
  if (ff_file->IsZombie()) {
    cout << "File " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_" << TString(year_label) << ".root not found" << endl;
    cout << "Quitting... " << endl;
    exit(-1);
  }

  FakeFactor* ff = (FakeFactor*)ff_file->Get("ff_comb");
  std::shared_ptr<RooWorkspace> ff_ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  ff_ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  fns_["ff_mt_medium_dmbins"] = std::shared_ptr<RooFunctor>(ff_ws_->function("ff_mt_medium_dmbins")->functor(ff_ws_->argSet("pt,dm,njets,m_pt,os,met,mt,m_iso,pass_single,mvis")));
  fns_["ff_mt_medium_mvadmbins"] = std::shared_ptr<RooFunctor>(ff_ws_->function("ff_mt_medium_mvadmbins")->functor(ff_ws_->argSet("pt,mvadm,njets,m_pt,os,met,mt,m_iso,pass_single,mvis")));



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
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut       = cfg.get<float>("dzTauCut");

  // tau energy scale corrections
  const float shift_tes_1prong     = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0      = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong     = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  const float shift_tes_3prong1pi0 = cfg.get<float>("TauEnergyScaleShift_ThreeProng");

  const float shift_tes_1prong_e = cfg.get<float>("TauEnergyScaleShift_OneProng_Error");
  const float shift_tes_1p1p0_e  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_Error");
  const float shift_tes_3prong_e = cfg.get<float>("TauEnergyScaleShift_ThreeProng_Error");
  const float shift_tes_3prong1p0_e = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0_Error");
  

	// lep->tau FES correction and uncertainties
	TFile TauFES_file(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_"+year_label+".root"); 
	TGraphAsymmErrors* FES_graph = (TGraphAsymmErrors*) TauFES_file.Get("fes");
	
	// apply non-zero only for DY MC in et channel, HPS decay modes 0 and 1	
	const float shift_tes_lepfake_1prong_barrel((ch == "et" && isDY) ? FES_graph->GetY()[0] - 1.0 : 0.0 );
	const float shift_tes_lepfake_1p1p0_barrel ((ch == "et" && isDY) ? FES_graph->GetY()[1] - 1.0 : 0.0 );
	const float shift_tes_lepfake_1prong_endcap((ch == "et" && isDY) ? FES_graph->GetY()[2] - 1.0 : 0.0 );
	const float shift_tes_lepfake_1p1p0_endcap ((ch == "et" && isDY) ? FES_graph->GetY()[3] - 1.0 : 0.0 );
	
	// for et take up/down values from file, for mt set to 1% (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2#Corrections_to_be_applied_to_AN2)
	// with the definition as below it can be accidently applied to data/embedded, so make sure to not do this!
	const float shift_tes_lepfake_1prong_barrel_up((ch == "et") ? FES_graph->GetErrorYhigh(0) : 0.01 );
	const float shift_tes_lepfake_1p1p0_barrel_up ((ch == "et") ? FES_graph->GetErrorYhigh(1) : 0.01 );
	const float shift_tes_lepfake_1prong_endcap_up((ch == "et") ? FES_graph->GetErrorYhigh(2) : 0.01 );
	const float shift_tes_lepfake_1p1p0_endcap_up ((ch == "et") ? FES_graph->GetErrorYhigh(3) : 0.01 );
	
	const float shift_tes_lepfake_1prong_barrel_down((ch == "et") ? FES_graph->GetErrorYlow(0) : 0.01 );
	const float shift_tes_lepfake_1p1p0_barrel_down ((ch == "et") ? FES_graph->GetErrorYlow(1) : 0.01 );
	const float shift_tes_lepfake_1prong_endcap_down((ch == "et") ? FES_graph->GetErrorYlow(2) : 0.01 );
	const float shift_tes_lepfake_1p1p0_endcap_down ((ch == "et") ? FES_graph->GetErrorYlow(3) : 0.01 );

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

 // lepton cuts
  const float ptLeptonLowCut   = cfg.get<float>("pt" + lep + "LowCut");
  const float ptLeptonHighCut  = cfg.get<float>("pt" + lep + "HighCut");
  const float etaLeptonCut     = cfg.get<float>("eta" + lep + "Cut");
  const float dxyLeptonCut     = cfg.get<float>("dxy" + lep + "Cut");
  const float dzLeptonCut      = cfg.get<float>("dz" + lep + "Cut");

  cout<<"dxyLeptonCut "<<dxyLeptonCut<<endl;
  cout<<"dzLeptonCut "<<dzLeptonCut<<endl;
  cout<<"dzTauCut "<<dzTauCut<<endl;

  // const bool  ApplyLeptonId    = cfg.get<bool>("Apply" + lep + "Id");

  const float dRTrigMatch = cfg.get<float>("dRTrigMatch");
  const float dRiso = cfg.get<float>("dRiso");
  const float lepIsoCut = cfg.get<float>("LeptonIsoCut");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");

  // check overlap
  const bool checkOverlap = cfg.get<bool>("CheckOverlap");

  // correction workspace
  const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");

  bool triggerEmbed2017 = false;
  float ptTriggerEmbed2017 = 40;
  float etaTriggerEmbed2017 = 1.479;

  if (era==2017&&isEmbedded)
    triggerEmbed2017 = true;

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
  std::string initNtupleName("initroottree/AC1B");
  std::string TauSpinnerWeightTreeName("icTauSpinnerProducer/TauSpinnerWeightTree");

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
  Synch17Tree *otree = new Synch17Tree(tree);
  initializeCPvar(otree);  
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

  std::cout << "inputFile_visPtResolution : " << std::endl;

  //Systematics init
  
  TauScaleSys *tauScaleSys = 0;
  TauOneProngScaleSys *tauOneProngScaleSys = 0;
  TauOneProngOnePi0ScaleSys *tauOneProngOnePi0ScaleSys = 0;
  TauThreeProngScaleSys *tauThreeProngScaleSys = 0;
  TauThreeProngOnePi0ScaleSys *tauThreeProngOnePi0ScaleSys = 0;
  MuonScaleSys *muonScaleSys = 0;
  ElectronScaleSys *electronScaleSys = 0;

  LepTauFakeOneProngScaleSys *lepTauFakeOneProngScaleSys = 0;
  LepTauFakeOneProngOnePi0ScaleSys *lepTauFakeOneProngOnePi0ScaleSys = 0;

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

    muonScaleSys = new MuonScaleSys(otree);
    muonScaleSys->SetUseSVFit(ApplySVFit);
    muonScaleSys->SetUseFastMTT(ApplyFastMTT);
    muonScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    muonScaleSys->SetUsePuppiMET(usePuppiMET);
    
    electronScaleSys = new ElectronScaleSys(otree);
    electronScaleSys->SetUseSVFit(ApplySVFit);
    electronScaleSys->SetUseFastMTT(ApplyFastMTT);
    electronScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    electronScaleSys->SetUsePuppiMET(usePuppiMET);

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
    tauThreeProngOnePi0ScaleSys->SetScale(shift_tes_3prong,shift_tes_3prong_e);
    tauThreeProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
    tauThreeProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauThreeProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);

    if (isDY) {

      lepTauFakeOneProngScaleSys = new LepTauFakeOneProngScaleSys(otree);
      lepTauFakeOneProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
      lepTauFakeOneProngScaleSys->SetUseSVFit(ApplySVFit);
      lepTauFakeOneProngScaleSys->SetUseFastMTT(ApplyFastMTT);
      lepTauFakeOneProngScaleSys->SetUsePuppiMET(usePuppiMET);
			lepTauFakeOneProngScaleSys->SetBarrelEdge(1.5);
			lepTauFakeOneProngScaleSys->SetScaleBarrelUp(shift_tes_lepfake_1prong_barrel, shift_tes_lepfake_1prong_barrel_up);
			lepTauFakeOneProngScaleSys->SetScaleBarrelDown(shift_tes_lepfake_1prong_barrel, shift_tes_lepfake_1prong_barrel_down);
			lepTauFakeOneProngScaleSys->SetScaleEndcapUp(shift_tes_lepfake_1prong_endcap, shift_tes_lepfake_1prong_endcap_up);
			lepTauFakeOneProngScaleSys->SetScaleEndcapDown(shift_tes_lepfake_1prong_endcap, shift_tes_lepfake_1prong_endcap_down);

      lepTauFakeOneProngOnePi0ScaleSys = new LepTauFakeOneProngOnePi0ScaleSys(otree);
      lepTauFakeOneProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
      lepTauFakeOneProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
      lepTauFakeOneProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
      lepTauFakeOneProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);
			lepTauFakeOneProngOnePi0ScaleSys->SetBarrelEdge(1.5);
			lepTauFakeOneProngOnePi0ScaleSys->SetScaleBarrelUp(shift_tes_lepfake_1p1p0_barrel, shift_tes_lepfake_1p1p0_barrel_up);
			lepTauFakeOneProngOnePi0ScaleSys->SetScaleBarrelDown(shift_tes_lepfake_1p1p0_barrel, shift_tes_lepfake_1p1p0_barrel_down);
			lepTauFakeOneProngOnePi0ScaleSys->SetScaleEndcapUp(shift_tes_lepfake_1p1p0_endcap, shift_tes_lepfake_1p1p0_endcap_up);
			lepTauFakeOneProngOnePi0ScaleSys->SetScaleEndcapDown(shift_tes_lepfake_1p1p0_endcap, shift_tes_lepfake_1p1p0_endcap_down);

    }

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
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
	else if (era==2017)
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
	else 
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
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
    
    double * TSweight = new double[expectedtauspinnerweights];
    TTree  * _treeTauSpinnerWeights = NULL;
    
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


    // Read in HLT filter
    vector<string> filterSingleLep = cfg.get<vector<string>>("filterSingleLep");
    vector<string> filterXtriggerLepLeg;
    vector<string> filterXtriggerTauLeg;
    if (era != 2018) {
      filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg");
      filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg");
    }
      
    cout<<"Number of single lepton trigger legs = "<<filterSingleLep.size()<<endl;
    cout<<"Number of X trigger legs (lep leg)   = "<<filterXtriggerLepLeg.size()<<endl;
    cout<<"Number of X trigger legs (tau leg)   = "<<filterXtriggerTauLeg.size()<<endl;

    ///////////////EVENT LOOP///////////////
    Long64_t numberOfEntries = analysisTree.GetEntries();

    for (Long64_t iEntry = 0; iEntry < numberOfEntries; iEntry++) {
      counter[0]++;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (era == 2018) {
      	if(isData && !isEmbedded)
      	  {
            if (ch == "mt")
            {
              if (analysisTree.event_run < 315974) { // muon filter of the mutau triggers changed
                filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_before_run315974");
                filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_before_run315974");
              }
              else if (analysisTree.event_run < 317509) // HPS algorithm was introduced
              {
                filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_run315974_to_HPS");
                filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_run315974_to_HPS");    
              }
              else
              {
                filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_after_HPS");
                filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_after_HPS");    
              }
            }
            else if (ch == "et")
            {
              if (analysisTree.event_run < 317509) // HPS algorithm was introduced
              {
                filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_before_HPS");
                filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_before_HPS");    
              }
              else
              {
                filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_after_HPS");
                filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_after_HPS");    
              }
        	  }
          }
      	else
        { // use with HPS in its name for MC and Embedded
      	  filterXtriggerLepLeg = cfg.get<vector<string>>("filterXtriggerLepLeg_after_HPS");
      	  filterXtriggerTauLeg = cfg.get<vector<string>>("filterXtriggerTauLeg_after_HPS");          
      	}
      }    

      //TO DO: fix tauspinner weight implementation after deprecated method is not in use for 2017
      if(!applyTauSpinnerWeights){
      	for(int tsindex = 0; tsindex < expectedtauspinnerweights; tsindex++) 
          TSweight[tsindex] = 1;

	otree->TauSpinnerWeightsEven = TSweight[0];
      	gentree->sm_htt125 = TSweight[0];
      	gentreeForGoodRecoEvtsOnly->sm_htt125 = TSweight[0];

      	otree->TauSpinnerWeightsMaxMix = TSweight[1];
      	gentree->mm_htt125 = TSweight[1];
      	gentreeForGoodRecoEvtsOnly->mm_htt125 = TSweight[1];

      	otree->TauSpinnerWeightsOdd = TSweight[2];
      	gentree->ps_htt125 = TSweight[2];
      	gentreeForGoodRecoEvtsOnly->ps_htt125 = TSweight[2];

      	otree->TauSpinnerWeightsMinusMaxMix = TSweight[3];
      	gentree->minusmm_htt125 = TSweight[3];
      	gentreeForGoodRecoEvtsOnly->minusmm_htt125 = TSweight[3];

      	otree->TauSpinnerWeightsMix0p375 = TSweight[4];
      	gentree->mix0p375_htt125 = TSweight[4];
      	gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = TSweight[4];
      
      }
      else{
	otree->TauSpinnerWeightsEven = analysisTree.TauSpinnerWeight[0];
      	gentree->sm_htt125 = analysisTree.TauSpinnerWeight[0];
      	gentreeForGoodRecoEvtsOnly->sm_htt125 = analysisTree.TauSpinnerWeight[0];

      	otree->TauSpinnerWeightsMaxMix = analysisTree.TauSpinnerWeight[1];
      	gentree->mm_htt125 = analysisTree.TauSpinnerWeight[1];
      	gentreeForGoodRecoEvtsOnly->mm_htt125 = analysisTree.TauSpinnerWeight[1];

      	otree->TauSpinnerWeightsOdd = analysisTree.TauSpinnerWeight[2];
      	gentree->ps_htt125 = analysisTree.TauSpinnerWeight[2];
      	gentreeForGoodRecoEvtsOnly->ps_htt125 = analysisTree.TauSpinnerWeight[2];

	if(analysisTree.TauSpinAngles_count>=5){
	  otree->TauSpinnerWeightsMinusMaxMix = analysisTree.TauSpinnerWeight[3];
	  gentree->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
	  gentreeForGoodRecoEvtsOnly->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
	  
	  otree->TauSpinnerWeightsMix0p375 = analysisTree.TauSpinnerWeight[4];
	  gentree->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
	  gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
	}
      }

      if (!isData){
      	FillGenTree(&analysisTree,gentree);
      	gentree->Fill();
      }

      //      std::cout << "OK!!!!!!!!" << std::endl;

      //Skip events not passing the MET filters, if applied
      bool passed_all_met_filters = passedAllMetFilters(&analysisTree, met_filters_list);
      if (ApplyMetFilters && !Synch && !passed_all_met_filters) continue;
      otree->passedAllMetFilters = passed_all_met_filters;
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
    
      //      if (overlapEvent && checkOverlap) continue;
      //      nonOverlap++;
      //      counter[2]++;
    
      if ((isData || isEmbedded) && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;
    
      initializeGenTree(gentree);
    
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;

      // embedded weight
      otree->embweight = 1;
      if (isEmbedded) {
      	otree->embweight = getEmbeddedWeight(&analysisTree, w);
	if (otree->embweight>10.0)
	  cout << "warning : embedding weight = " << otree->embweight << endl;
      }

      // tau selection
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it < analysisTree.tau_count; ++it) {
	float shift_tes  = shift_tauES(&analysisTree,it,
				       shift_tes_1prong,
				       shift_tes_1p1p0,
				       shift_tes_3prong,
				       shift_tes_lepfake_1prong_barrel,
				       shift_tes_lepfake_1p1p0_barrel,
				       shift_tes_lepfake_1prong_endcap,
				       shift_tes_lepfake_1p1p0_endcap);
 
	float ptTau = (1.+shift_tes)*analysisTree.tau_pt[it];

        if (ptTau <= ptTauLowCut) continue;
        if (fabs(analysisTree.tau_eta[it]) >= etaTauCut) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it]) >= dzTauCut) continue;
        if (fabs(fabs(analysisTree.tau_charge[it]) - 1) > 0.001) continue;

      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] < 0.5) continue;
      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[it] < 0.5) continue;
      	if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] < 0.5) continue;
        if (analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) continue; //always true, cut applied in NTupleMaker
        if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
	//	if (analysisTree.tau_MVADM2017v1[it] < 0) continue; //prevents storing events with unidentified mva DM for the tau (-1)
    
        taus.push_back(it);
      }
      counter[3]++;
    
      // as of 30 June 2020 electron ES in Embedded samples isn't corrected in BigNTuples, so it needs to be applied downstream
      float sf_eleES = 1.0;   

      //lepton selection
      vector<int> leptons; leptons.clear();
      if(ch == "et"){
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
          bool electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie];
          
          if (isEmbedded) sf_eleES = EmbedElectronES_SF(&analysisTree, era, ie);          
          if (sf_eleES*analysisTree.electron_pt[ie] <= ptLeptonLowCut) continue;
          if (fabs(analysisTree.electron_eta[ie]) >= etaLeptonCut) continue;
          if (fabs(analysisTree.electron_dxy[ie]) >= dxyLeptonCut) continue;
      	  if (fabs(analysisTree.electron_dz[ie]) >= dzLeptonCut) continue;
	  float relIsoLep = (abs_Iso_et(ie, &analysisTree, dRiso) / (sf_eleES*analysisTree.electron_pt[ie]) );
	  if (relIsoLep>lepIsoCut) continue;
	  if (!electronMvaId) continue;
      	  if (!analysisTree.electron_pass_conversion[ie]) continue;
      	  if (analysisTree.electron_nmissinginnerhits[ie] > 1) continue;
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
          if (!muonMediumId) continue;
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
	float shift_tes  = shift_tauES(&analysisTree,tIndex,
				       shift_tes_1prong,
				       shift_tes_1p1p0,
				       shift_tes_3prong,
				       shift_tes_lepfake_1prong_barrel,
				       shift_tes_lepfake_1p1p0_barrel,
				       shift_tes_lepfake_1prong_endcap,
				       shift_tes_lepfake_1p1p0_endcap);
 
	float ptTau = (1.+shift_tes)*analysisTree.tau_pt[tIndex];
    
    	//////////////LOOP on Leptons or second Tau/////////////
    
        for (unsigned int il = 0; il < leptons.size(); ++il) {
          unsigned int lIndex = leptons.at(il);
    
          float lep_pt     = -9999.;
          float lep_eta    = -9999.;
          float lep_phi    = -9999.;
	  float relIsoLep  = -9999.;
             
    
          if(ch == "mt"){
            relIsoLep = (abs_Iso_mt(lIndex, &analysisTree, dRiso) / analysisTree.muon_pt[lIndex] );
            lep_pt  = analysisTree.muon_pt[lIndex]; 
            lep_eta = analysisTree.muon_eta[lIndex]; 
            lep_phi = analysisTree.muon_phi[lIndex];
          }
          if(ch == "et"){    
            sf_eleES = 1.;
            if (isEmbedded) sf_eleES = EmbedElectronES_SF(&analysisTree, era, lIndex);  
            relIsoLep = (abs_Iso_et(lIndex, &analysisTree, dRiso) / (sf_eleES*analysisTree.electron_pt[lIndex]) );     
            lep_pt  = sf_eleES*analysisTree.electron_pt[lIndex];
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
		 else if (fabs(sortIsoTau - isoTauMax) < 1.e-5)
		   {
		     if (ptTau > tau_pt_max)
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
    
      //      std::cout << "OK1" << std::endl;

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
        sf_eleES = 1.;
        if (isEmbedded) sf_eleES = EmbedElectronES_SF(&analysisTree, era, leptonIndex);  
      	lep_pt  = sf_eleES*analysisTree.electron_pt[leptonIndex];
      	lep_eta = analysisTree.electron_eta[leptonIndex]; 
      	lep_phi = analysisTree.electron_phi[leptonIndex];
      }


      ////////////////////////////////////////////////////////////
      // Trigger matching
      ////////////////////////////////////////////////////////////
    
      vector<bool> isSingleLepLeg(filterSingleLep.size(), false);
      vector<bool> isXTrigLepLeg(filterXtriggerLepLeg.size(), false);
      vector<bool> isXTrigTauLeg(filterXtriggerTauLeg.size(), false);
      bool isSingleLepTrig = false;
      bool isXTrig         = false;
      bool isXTrigLep      = true;
      bool isXTrigTau      = true;
    
      otree->trg_singlemuon = false;
      otree->trg_singleelectron = false;
      otree->singleLepTrigger = false;
      otree->trg_doubletau = false;
      otree->trg_mutaucross = false;
      otree->trg_mutaucross_mu = false;
      otree->trg_mutaucross_tau = false;
      otree->trg_etaucross = false;
      otree->trg_etaucross_e = false;
      otree->trg_etaucross_tau = false;
        
      for (unsigned int iT = 0; iT < analysisTree.trigobject_count; ++iT) {
	
	/*
	for(unsigned int i_trig = 0; i_trig < filterXtriggerLepLeg.size(); i_trig++)
	  {
	    if (nXTrigLepLeg.at(i_trig) == -1) continue;
	    if (analysisTree.trigobject_filters[iT][nXTrigLepLeg.at(i_trig)]) 
	      cout << filterXtriggerLepLeg.at(i_trig) << " : " << analysisTree.trigobject_pt[iT] << endl;
	    
	  }       
	for(unsigned int i_trig = 0; i_trig < filterXtriggerTauLeg.size(); i_trig++)
	  {
	    if (nXTrigTauLeg.at(i_trig) == -1) continue;
	    if (analysisTree.trigobject_filters[iT][nXTrigTauLeg.at(i_trig)])
	      cout << filterXtriggerTauLeg.at(i_trig) << " : " << analysisTree.trigobject_pt[iT] << endl;
	  }
	*/
	float dRtrigLep = deltaR(lep_eta, lep_phi, analysisTree.trigobject_eta[iT], analysisTree.trigobject_phi[iT]);        
	float dRtrigTau = deltaR(analysisTree.tau_eta[tauIndex], analysisTree.tau_phi[tauIndex], analysisTree.trigobject_eta[iT], analysisTree.trigobject_phi[iT]);        
    

	if (dRtrigLep < dRTrigMatch){
	  for(unsigned int i_trig = 0; i_trig < filterSingleLep.size(); i_trig++)
	    {
              if (nSingleLepTrig.at(i_trig) == -1) continue;
              if (analysisTree.trigobject_filters[iT][nSingleLepTrig.at(i_trig)]) isSingleLepLeg.at(i_trig) = true;
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
      /*
      for (unsigned int i=0; i<nXTrigTauLeg.size(); ++i) {
	cout << "Tau leg : " << filterXtriggerTauLeg.at(i) << " : " << nXTrigTauLeg.at(i) << endl;
      }
      for (unsigned int i=0; i<nXTrigLepLeg.size(); ++i) {
	cout << "Lepton leg : " << filterXtriggerLepLeg.at(i) << " : " << nXTrigLepLeg.at(i) << endl;
      }
      */

      if (era == 2017 && ch == "et") {
        int id_SingleEGO = -1;
        int id_Single32 = -1;    
        for(unsigned int i_trig = 0; i_trig < filterSingleLep.size(); i_trig++){
          if(filterSingleLep.at(i_trig) == "hltEGL1SingleEGOrFilter")
            id_SingleEGO = i_trig;
          if(filterSingleLep.at(i_trig) == "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")
            id_Single32 = i_trig;
        }
        isSingleLepLeg[id_SingleEGO] = isSingleLepLeg[id_SingleEGO] && isSingleLepLeg[id_Single32];
        isSingleLepLeg[id_Single32] = isSingleLepLeg[id_SingleEGO] && isSingleLepLeg[id_Single32];
      }
       
      for(unsigned int i_trig = 0; i_trig < filterSingleLep.size(); i_trig++)
        isSingleLepTrig = isSingleLepTrig || isSingleLepLeg.at(i_trig);
      for(unsigned int i_trig = 0; i_trig < filterXtriggerTauLeg.size(); i_trig++)
        isXTrigTau = isXTrigTau && isXTrigTauLeg.at(i_trig);
      for(unsigned int i_trig = 0; i_trig < filterXtriggerLepLeg.size(); i_trig++)
        isXTrigLep = isXTrigLep && isXTrigLepLeg.at(i_trig);
      isXTrig = isXTrigTau && isXTrigLep;
    
      otree->singleLepTrigger = isSingleLepTrig;
      otree->xTrigger = isXTrig;
      otree->xTriggerTau = isXTrigTau;
      otree->xTriggerLep = isXTrigLep;
      
      if (ch == "mt")
      {
        otree->trg_singlemuon = isSingleLepTrig;
        otree->trg_mutaucross_mu = isXTrigLep;
        otree->trg_mutaucross_tau = isXTrigTau;
        otree->trg_mutaucross = isXTrig;        
	otree->singleLepTrigger = otree->trg_singlemuon;
      }
      if (ch == "et")
      {  
	otree->trg_singleelectron = isSingleLepTrig;
	otree->trg_etaucross_e = isXTrigLep;
	otree->trg_etaucross_tau = isXTrigTau;
	otree->trg_etaucross = isXTrig;        
	if (triggerEmbed2017) {
	  if (lep_pt<ptTriggerEmbed2017&&fabs(lep_eta)>etaTriggerEmbed2017) {
	    otree->trg_singleelectron = true;
			otree->trg_etaucross = true;
	    // otree->trg_etaucross_e = true;
	    // otree->trg_etaucross = otree->trg_etaucross_tau;
	  }
	}
	otree->singleLepTrigger = otree->trg_singleelectron;
      } 
      /*
      cout << "xtrig_lep = " << isXTrigLep << "  xtrig_tau = " << isXTrigTau << endl;
      cout << "lep pt = " << lep_pt << "  eta = " << lep_eta << endl;
      cout << "tau pt = " 
	   << analysisTree.tau_pt[tauIndex] << "  eta = " 
	   << analysisTree.tau_eta[tauIndex] << endl;
      
      cout << endl;
      */

    
    
      ////////////////////////////////////////////////////////////
      // Filling variables
      ////////////////////////////////////////////////////////////

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
      double eff_data_trig_lt_tauUp   = 1;
      double eff_mc_trig_lt_tauUp     = 1;
      double eff_data_trig_lt_tauDown = 1;
      double eff_mc_trig_lt_tauDown   = 1;
      // reset efficiency weights
    
      //all criterua passed, we fill vertices here;	
      FillVertices(&analysisTree, otree, isData, leptonIndex, tauIndex, ch);
    
      //Merijn: save here all gen information for the selected RECO events, gets stored for convenience in the taucheck tree ;-). Note that no selection on gen level is applied..     
      //Merijn 2019 4 3:note that a separate fill is not needed. We store in the otree now, which is Filled at the bottom! Filling here will make things out of synch..
      if (!isData)
        FillGenTree(&analysisTree, gentreeForGoodRecoEvtsOnly);       
    
      if(ch == "mt") {
      	FillMuTau(&analysisTree, otree, leptonIndex, tauIndex, dRiso);
        leptonLV.SetXYZM(analysisTree.muon_px[leptonIndex], analysisTree.muon_py[leptonIndex], analysisTree.muon_pz[leptonIndex], muonMass);
      } 
      else if(ch == "et"){
      	FillETau(&analysisTree, otree, leptonIndex, dRiso, era, isEmbedded);
        leptonLV.SetXYZM(analysisTree.electron_px[leptonIndex], analysisTree.electron_py[leptonIndex], analysisTree.electron_pz[leptonIndex], electronMass);
        if(isEmbedded) leptonLV *= sf_eleES;
      }
      
      FillTau(&analysisTree, otree, leptonIndex, tauIndex);

      // initialize JER (including data and embedded) 
      otree->apply_recoil = ApplyRecoilCorrections;
      jets::initializeJER(&analysisTree);

      if (!isData && !isEmbedded) { // JER smearing
	jets::associateRecoAndGenJets(&analysisTree, resolution);
	jets::smear_jets(&analysisTree,resolution,resolution_sf,true);
      }

      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
  
      ////////////////////////////////////////////////////////////
      // ID/Iso and Trigger Corrections
      ////////////////////////////////////////////////////////////

      // setting weights to 1
      otree->trkeffweight = 1;
      otree->trigweight_1 = 1;
      otree->trigweight_2 = 1;
      otree->idisoweight_1 = 1;
      otree->idisoweight_antiiso_1 = 1;
      otree->idisoweight_2 = 1;
      otree->idisoweight_antiiso_2 = 1;
      otree->trigweight = 1;
      otree->effweight = 1;
      otree->puweight = 1; 
      otree->mcweight = 1;
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

      if ((!isData || isEmbedded) && ApplyLepSF) {
      	TString suffix = "mc";
      	TString suffixRatio = "ratio";
      	if (isEmbedded) {suffix = "embed"; suffixRatio = "embed_ratio";}
	TString mvadm = TString::Itoa(analysisTree.tau_MVADM2017v1[tauIndex],10);
	if (analysisTree.tau_MVADM2017v1[tauIndex]<0.0)
	  mvadm = TString::Itoa(analysisTree.tau_decayMode[tauIndex],10);;
	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex]);
	w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex]);
	w->var("t_phi")->setVal(analysisTree.tau_phi[tauIndex]);
	//w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex]);
	
	if (analysisTree.tau_MVADM2017v1[tauIndex]<0.)
	  w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex]);
	else
	  w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex]);
	
	if (ch == "mt") {
	  w->var("m_pt")->setVal(leptonLV.Pt());
	  w->var("m_eta")->setVal(leptonLV.Eta());
	  eff_data_trig_lt_tau = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_data")->getVal();
	  eff_mc_trig_lt_tau = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_" + suffix)->getVal();
	  eff_data_trig_lt_tauUp = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_data_mvadm"+mvadm+"_up")->getVal();
	  eff_mc_trig_lt_tauUp = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_" + suffix + "_mvadm"+mvadm+"_up")->getVal();
	  eff_data_trig_lt_tauDown = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_data_mvadm"+mvadm+"_down")->getVal();
	  eff_mc_trig_lt_tauDown = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_" + suffix + "_mvadm"+mvadm+"_down")->getVal();
	  eff_data_trig_L = w->function("m_trg_ic_data")->getVal();
	  eff_mc_trig_L = w->function("m_trg_ic_" + suffix)->getVal();
	  if (era == 2016) {
	    eff_data_trig_lt_l = w->function("m_trg_19_ic_data")->getVal();
	    eff_mc_trig_lt_l = w->function("m_trg_19_ic_" + suffix)->getVal();
	  }
	  else {
	    eff_data_trig_lt_l = w->function("m_trg_20_ic_data")->getVal();
	    eff_mc_trig_lt_l = w->function("m_trg_20_ic_" + suffix)->getVal();
	  }
	  otree->idisoweight_1 = w->function("m_idiso_ic_" + suffixRatio)->getVal();
	  otree->idisoweight_antiiso_1 = w->function("m_idiso_ic_" + suffixRatio)->getVal();
	  otree->trkeffweight = w->function("m_trk_ratio")->getVal();
	}
	else if (ch == "et") {
	  w->var("e_pt")->setVal(leptonLV.Pt());
	  w->var("e_eta")->setVal(leptonLV.Eta());
	  eff_data_trig_L = w->function("e_trg_ic_data")->getVal();
	  eff_mc_trig_L = w->function("e_trg_ic_" + suffix)->getVal();
	  if (era > 2016) {
	    eff_data_trig_lt_tau = w->function("t_trg_ic_deeptau_medium_mvadm_etau_data")->getVal();
	    eff_mc_trig_lt_tau = w->function("t_trg_ic_deeptau_medium_mvadm_etau_" + suffix)->getVal();
	    eff_data_trig_lt_l = w->function("e_trg_24_ic_data")->getVal();
	    eff_mc_trig_lt_l = w->function("e_trg_24_ic_" + suffix)->getVal();
	    eff_data_trig_lt_tauUp = w->function("t_trg_ic_deeptau_medium_mvadm_etau_data_mvadm"+mvadm+"_up")->getVal();
	    eff_mc_trig_lt_tauUp = w->function("t_trg_ic_deeptau_medium_mvadm_etau_" + suffix + "_mvadm"+mvadm+"_up")->getVal();
	    eff_data_trig_lt_tauDown = w->function("t_trg_ic_deeptau_medium_mvadm_etau_data_mvadm"+mvadm+"_down")->getVal();
	    eff_mc_trig_lt_tauDown = w->function("t_trg_ic_deeptau_medium_mvadm_etau_" + suffix + "_mvadm"+mvadm+"_down")->getVal();
	    if (triggerEmbed2017) {
	      if (leptonLV.Pt()<ptTriggerEmbed2017&&fabs(leptonLV.Eta())>etaTriggerEmbed2017) {
		eff_mc_trig_L = 1.0;
		eff_mc_trig_lt_l = 1.0;
		eff_mc_trig_lt_tau = 1.0;
	      }
	    }
	  }
	  else {
	    eff_data_trig_lt_tau = 0;
	    eff_mc_trig_lt_tau = 0;
	    eff_data_trig_lt_l = 0;
	    eff_mc_trig_lt_l = 0;
	    eff_data_trig_lt_tauUp = 0;
	    eff_mc_trig_lt_tauUp = 0;
	    eff_data_trig_lt_tauDown = 0;
	    eff_mc_trig_lt_tauDown = 0;
	  }
	  otree->idisoweight_1 = w->function("e_idiso_ic_" + suffixRatio)->getVal();
	  otree->idisoweight_antiiso_1 = w->function("e_idiso_ic_" + suffixRatio)->getVal();
	  otree->trkeffweight = w->function("e_trk_" + suffixRatio)->getVal();
	}
	otree->trigweight_1 = 1;
	otree->trigweight_2 = 1;

	//	if (eff_mc_trig_L>0.1)
	otree->trigweight_1 = eff_data_trig_L/eff_mc_trig_L;
	//	if (eff_mc_trig_lt_l>0.1&&eff_mc_trig_lt_tau>0.1) 
	otree->trigweight_2 = (eff_data_trig_lt_l*eff_data_trig_lt_tau)/(eff_mc_trig_lt_l*eff_mc_trig_lt_tau);
	
	if (std::isnan(otree->trigweight_1))
	  otree->trigweight_1 = eff_data_trig_L;
	if (std::isnan(otree->trigweight_2))
	  otree->trigweight_2 = eff_data_trig_lt_l*eff_data_trig_lt_tau;

	/*
	if (leptonLV.Pt()>26.&&leptonLV.Pt()<30.&&fabs(leptonLV.Eta())>1.65
	    &&otree->trg_etaucross_e) {
	  
	  std::cout << "run : " << otree->run << "    event : " << otree->evt << std::endl;
	  std::cout << "trg_singleelectron : " << otree->trg_singleelectron << std::endl;
	  std::cout << "trg_etaucross_e    : " << otree->trg_etaucross_e << endl;
	  std::cout << "trg_etaucross_tau  : " << otree->trg_etaucross_tau << endl;
	  std::cout << "electron pt = " << leptonLV.Pt() << "   eta = " << leptonLV.Eta() << std::endl;
	  std::cout << "eff(trig_L,Data) = " << eff_data_trig_L 
		    << "    eff(trig_L,MC) = " << eff_mc_trig_L << std::endl;
	  std::cout << "eff(Data)/Eff(MC) = " << otree->trigweight_1 << std::endl;
	  std::cout << "eff(trig_l,Data) = " << eff_data_trig_lt_l 
		    << "    eff(trig_l,MC) = " << eff_mc_trig_lt_l << std::endl;

	  std::cout << "eff(trig_t,Data) = " << eff_data_trig_lt_tau 
		    << "    eff(trig_t,MC) = " << eff_mc_trig_lt_tau << std::endl;
	  if (isEmbedded)
	    std::cout << "SF(trig_L) = " << w->function("e_trg_ic_embed_ratio")->getVal();
	  else
	    std::cout << "SF(trig_L) = " << w->function("e_trg_ic_ratio")->getVal();
	  std::cout << std::endl;
	}
	*/

	double eff_data_trig = eff_data_trig_L + (eff_data_trig_lt_l - eff_data_trig_L) * eff_data_trig_lt_tau;
	double eff_mc_trig = eff_mc_trig_L + (eff_mc_trig_lt_l - eff_mc_trig_L) * eff_mc_trig_lt_tau;                                                            
	double eff_data_trigUp = eff_data_trig_L + (eff_data_trig_lt_l - eff_data_trig_L) * eff_data_trig_lt_tauUp;
	double eff_mc_trigUp = eff_mc_trig_L + (eff_mc_trig_lt_l - eff_mc_trig_L) * eff_mc_trig_lt_tauUp;                                                   
	double eff_data_trigDown = eff_data_trig_L + (eff_data_trig_lt_l - eff_data_trig_L) * eff_data_trig_lt_tauDown;
	double eff_mc_trigDown = eff_mc_trig_L + (eff_mc_trig_lt_l - eff_mc_trig_L) * eff_mc_trig_lt_tauDown;

	

	if (eff_data_trig > 1e-4 && eff_mc_trig > 1e-4){
	  otree->trigweight = eff_data_trig / eff_mc_trig;
	  double trigweightUp = (eff_data_trigUp / eff_mc_trigUp) / (eff_data_trig / eff_mc_trig);
	  double trigweightDown = (eff_data_trigDown / eff_mc_trigDown) / (eff_data_trig / eff_mc_trig);                
	  if (ch == "mt"){
	    if(mvadm=="0"){
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM0_13TeVUp = trigweightUp;
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM0_13TeVDown = trigweightDown;
	    }else if(mvadm=="1"){
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM1_13TeVUp = trigweightUp;
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM1_13TeVDown = trigweightDown;
	    }else if(mvadm=="2"){
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM2_13TeVUp = trigweightUp;
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM2_13TeVDown = trigweightDown;
	    }else if(mvadm=="10"){
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp = trigweightUp;
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown = trigweightDown;
	    }else if(mvadm=="11"){
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM11_13TeVUp = trigweightUp;
	      otree->weight_CMS_eff_Xtrigger_mt_MVADM11_13TeVDown = trigweightDown;
	    }
	  }
	}      
      }
      counter[10]++;
    
      //cout <<"TauID SF" <<endl;

      if ((!isData || isEmbedded) && analysisTree.tau_genmatch[tauIndex] <= 5) { 
      	TString suffix = "";
      	if (isEmbedded) suffix = "_embed";
	
	TString mvadm = TString::Itoa(analysisTree.tau_MVADM2017v1[tauIndex],10);
	if (analysisTree.tau_MVADM2017v1[tauIndex]<0.0)
	  mvadm = TString::Itoa(analysisTree.tau_decayMode[tauIndex],10);;
	double t_pt = analysisTree.tau_pt[tauIndex];
	
      	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex]);
	if (analysisTree.tau_MVADM2017v1[tauIndex]<0.0) 
	  w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex]);
	else
	  w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex]);


	//corrections for genuine taus
	if( analysisTree.tau_genmatch[tauIndex] == 5 ){
	  TString forEtau="";
	  if(ch == "et") forEtau = "_tightvsele";
	  double nominalID= w->function("t_deeptauid_mvadm"+suffix+"_medium"+forEtau)->getVal();
	  otree->idisoweight_2 = nominalID;
	  //cout <<nominalID <<endl;
	  double tauIDlowpTUp = w->function("t_deeptauid_mvadm"+suffix+"_medium"+forEtau+"_lowpt_mvadm"+mvadm+"_up")->getVal() / nominalID;
	  double tauIDhighpTUp = w->function("t_deeptauid_mvadm"+suffix+"_medium"+forEtau+"_highpt_mvadm"+mvadm+"_up")->getVal() / nominalID;
	  double tauIDlowpTDown = w->function("t_deeptauid_mvadm"+suffix+"_medium"+forEtau+"_lowpt_mvadm"+mvadm+"_down")->getVal() / nominalID;
	  double tauIDhighpTDown = w->function("t_deeptauid_mvadm"+suffix+"_medium"+forEtau+"_highpt_mvadm"+mvadm+"_down")->getVal() / nominalID;
	  
	  if(mvadm=="0"){
	    if(t_pt<40){
	      otree->weight_CMS_eff_t_pTlow_MVADM0_13TeVUp = tauIDlowpTUp;
	      otree->weight_CMS_eff_t_pTlow_MVADM0_13TeVDown = tauIDlowpTDown;
	    }else{
	      otree->weight_CMS_eff_t_pThigh_MVADM0_13TeVUp = tauIDhighpTUp;
	      otree->weight_CMS_eff_t_pThigh_MVADM0_13TeVDown = tauIDhighpTDown;
	    }
	  }else if(mvadm=="1"){
	    if(t_pt<40){
	      otree->weight_CMS_eff_t_pTlow_MVADM1_13TeVUp = tauIDlowpTUp;
	      otree->weight_CMS_eff_t_pTlow_MVADM1_13TeVDown = tauIDlowpTDown;
	    }else{
	      otree->weight_CMS_eff_t_pThigh_MVADM1_13TeVUp = tauIDhighpTUp;
	      otree->weight_CMS_eff_t_pThigh_MVADM1_13TeVDown = tauIDhighpTDown;
	    }
	  }else if(mvadm=="2"){
	    if(t_pt<40){
	      otree->weight_CMS_eff_t_pTlow_MVADM2_13TeVUp = tauIDlowpTUp;
	      otree->weight_CMS_eff_t_pTlow_MVADM2_13TeVDown = tauIDlowpTDown;
	    }else{
	      otree->weight_CMS_eff_t_pThigh_MVADM2_13TeVUp = tauIDhighpTUp;
	      otree->weight_CMS_eff_t_pThigh_MVADM2_13TeVDown = tauIDhighpTDown;
	    }
	  }else if(mvadm=="10"){
	    if(t_pt<40){
	      otree->weight_CMS_eff_t_pTlow_MVADM10_13TeVUp = tauIDlowpTUp;
	      otree->weight_CMS_eff_t_pTlow_MVADM10_13TeVDown = tauIDlowpTDown;
	    }else{
	      otree->weight_CMS_eff_t_pThigh_MVADM10_13TeVUp = tauIDhighpTUp;
	      otree->weight_CMS_eff_t_pThigh_MVADM10_13TeVDown = tauIDhighpTDown;
	    }
	  }else if(mvadm=="11"){
	    if(t_pt<40){
	      otree->weight_CMS_eff_t_pTlow_MVADM11_13TeVUp = tauIDlowpTUp;
	      otree->weight_CMS_eff_t_pTlow_MVADM11_13TeVDown = tauIDlowpTDown;
	    }else{
	      otree->weight_CMS_eff_t_pThigh_MVADM11_13TeVUp = tauIDhighpTUp;
	      otree->weight_CMS_eff_t_pThigh_MVADM11_13TeVDown = tauIDhighpTDown;
	    }
	  }
	}else if( analysisTree.tau_genmatch[tauIndex] == 2 || analysisTree.tau_genmatch[tauIndex] == 4 ){
	  //corrections for mu->tau fakes

	  double FRcorr = w->function("t_mufake_mt_mvadm")->getVal();
	  double FRcorrUp = w->function("t_mufake_mt_mvadm_mvadm"+mvadm+"_up")->getVal() / FRcorr;
	  double FRcorrDown = w->function("t_mufake_mt_mvadm_mvadm"+mvadm+"_down")->getVal() / FRcorr;

	  if(mvadm=="0"){
	    otree->weight_mufake_corr = FRcorr;
	    otree->weight_CMS_mufake_mt_MVADM0_13TeVUp = FRcorrUp;
	    otree->weight_CMS_mufake_mt_MVADM0_13TeVDown = FRcorrDown;
	  }else if(mvadm=="1"){
	    otree->weight_mufake_corr = FRcorr;
	    otree->weight_CMS_mufake_mt_MVADM1_13TeVUp = FRcorrUp;
	    otree->weight_CMS_mufake_mt_MVADM1_13TeVDown = FRcorrDown;
	  }else if(mvadm=="2"){
	    otree->weight_mufake_corr = FRcorr;
	    otree->weight_CMS_mufake_mt_MVADM2_13TeVUp = FRcorrUp;
	    otree->weight_CMS_mufake_mt_MVADM2_13TeVDown = FRcorrDown;
	  }else if(mvadm=="10"){
	    otree->weight_mufake_corr = FRcorr;
	    otree->weight_CMS_mufake_mt_MVADM10_13TeVUp = FRcorrUp;
	    otree->weight_CMS_mufake_mt_MVADM10_13TeVDown = FRcorrDown;
	  }else if(mvadm=="11"){
	    otree->weight_mufake_corr = FRcorr;
	    otree->weight_CMS_mufake_mt_MVADM11_13TeVUp = FRcorrUp;
	    otree->weight_CMS_mufake_mt_MVADM11_13TeVDown = FRcorrDown;
	  }
	}else if( analysisTree.tau_genmatch[tauIndex] == 1 || analysisTree.tau_genmatch[tauIndex] == 3 ){
	  //corrections for e->tau fakes
	  //Currently not present
	}
	  
      }
    
      //      cout << "\n======================" << endl;
      //      cout << "Trigger weight = " << otree->trigweight << endl;
      //      cout << "Trk eff weight = " << otree->trkeffweight << endl;
      //      cout << "id/Iso 1       = " << otree->idisoweight_1 << endl;
      //      cout << "id/Iso 2       = " << otree->idisoweight_2 << endl;
      //      cout << "======================\n" << endl;

      otree->effweight = otree->idisoweight_1 * otree->trkeffweight * otree->idisoweight_2 * otree->trigweight;
      otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
      
      
      //Theory uncertainties for CP analysis
      
      otree->weight_CMS_scale_gg_13TeVUp   = analysisTree.weightScale4;
      otree->weight_CMS_scale_gg_13TeVDown = analysisTree.weightScale8;

      otree->weight_CMS_PS_ISR_ggH_13TeVUp   = 1.;
      otree->weight_CMS_PS_ISR_ggH_13TeVDown = 1.;
      otree->weight_CMS_PS_FSR_ggH_13TeVUp   = 1.;
      otree->weight_CMS_PS_FSR_ggH_13TeVDown = 1.;

      if(isHiggs){
	otree->weight_CMS_PS_ISR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[6];
	otree->weight_CMS_PS_ISR_ggH_13TeVDown = analysisTree.gen_pythiaweights[8];
	otree->weight_CMS_PS_FSR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[7];
	otree->weight_CMS_PS_FSR_ggH_13TeVDown = analysisTree.gen_pythiaweights[9];
      }

      //Prefiring weights for CP analysis
      otree->prefiringweight     = analysisTree.prefiringweight;
      otree->prefiringweightUp   = analysisTree.prefiringweightup;
      otree->prefiringweightDown = analysisTree.prefiringweightdown;



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
      counter[11]++;

      ////////////////////////////////////////////////////////////
      // Lep->tau fake weight
      ////////////////////////////////////////////////////////////

      otree->mutaufakeweight = 1.;
      otree->etaufakeweight = 1.;
      /* //DEPREACTED
      if (!isData){
      
      	if (ch == "et") {
      	  otree->etaufakeweight = leptauFR->get_fakerate("electron", "Tight", otree->eta_2, otree->gen_match_2);
      	  otree->mutaufakeweight = leptauFR->get_fakerate("muon", "Loose", otree->eta_2, otree->gen_match_2);
      	}
    	  else if (ch == "mt") {
      	  otree->etaufakeweight = leptauFR->get_fakerate("electron", "VLoose", otree->eta_2, otree->gen_match_2);
      	  otree->mutaufakeweight = leptauFR->get_fakerate("muon", "Tight", otree->eta_2, otree->gen_match_2);
      	}
	}*/

      if(!isData){
	if(otree->gen_match_2==2||otree->gen_match_2==4){
	  TFile muTauFRfile(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSmu_"+year_label+".root"); 
	  TH1F *SFhist = (TH1F*) muTauFRfile.Get(LeptauFake_wpVsMu);
	  otree->mutaufakeweight = SFhist->GetBinContent(SFhist->GetXaxis()->FindBin(abs(otree->eta_2)));
	}else if(otree->gen_match_2==1||otree->gen_match_2==3){
	  TFile eTauFRfile(TString(cmsswBase)+"/src/TauPOG/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSe_"+year_label+".root"); 
	  TH1F *SFhist = (TH1F*) eTauFRfile.Get(LeptauFake_wpVsEle);
	  otree->etaufakeweight = SFhist->GetBinContent(SFhist->GetXaxis()->FindBin(abs(otree->eta_2)));
	}
      }
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

      otree->met_uncorr = otree->puppimet;
      otree->metphi_uncorr = otree->puppimetphi;
      otree->njetshad = otree->njets;
      if (isWJets) otree->njetshad += 1;

      if(ApplyRecoilCorrections){        
      	genV = genTools::genV(analysisTree);
      	genL = genTools::genL(analysisTree);

        genTools::KITRecoilCorrections( recoilCorrector, ApplyRecoilCorrections, // pass the value != 0 to apply corrections
          otree->puppimet, otree->puppimetphi,
          genV.Px(), genV.Py(),
          genL.Px(), genL.Py(),
          otree->njetshad,
          otree->met_rcmr, otree->metphi_rcmr
        );
        
        // overwriting with recoil-corrected values 
        otree->puppimet = otree->met_rcmr;
        otree->puppimetphi = otree->metphi_rcmr;   
      }
      
      //ditau sytem
      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
  				     analysisTree.tau_py[tauIndex],
  				     analysisTree.tau_pz[tauIndex],
  				     analysisTree.tau_mass[tauIndex]);
    
      // using PF MET
      TLorentzVector metLV, puppimetLV; 
      metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi), otree->met*TMath::Sin(otree->metphi), 0,
                    TMath::Sqrt( otree->met*TMath::Sin(otree->metphi)*otree->met*TMath::Sin(otree->metphi) +
  			               otree->met*TMath::Cos(otree->metphi)*otree->met*TMath::Cos(otree->metphi)));
      puppimetLV.SetXYZT(otree->puppimet*TMath::Cos(otree->puppimetphi), otree->puppimet*TMath::Sin(otree->puppimetphi), 0,
                    TMath::Sqrt( otree->puppimet*TMath::Sin(otree->puppimetphi)*otree->puppimet*TMath::Sin(otree->puppimetphi) +
  			               otree->puppimet*TMath::Cos(otree->puppimetphi)*otree->puppimet*TMath::Cos(otree->puppimetphi)));
    
    
      ////////////////////////////////////////////////////////////
      // Tau ES shift + propagate to MET
      ////////////////////////////////////////////////////////////
      
      if (!isData||isEmbedded) {
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
				// do only mu->tau FES for mt channel (but effectively it is 0, see its definition above) and e->tau FES for et channel
      	else if ((ch == "mt" && (otree->gen_match_2 == 2 || otree->gen_match_2 == 4)) || (ch == "et" && (otree->gen_match_2 == 1 || otree->gen_match_2 == 3))) {
      	  if (otree->tau_decay_mode_2 == 0){
						if (fabs(otree->eta_2) < 1.5)
							shift_tes = shift_tes_lepfake_1prong_barrel; 
						else
							shift_tes = shift_tes_lepfake_1prong_endcap; 
						isOneProng = true;
					}
      	  else if (otree->tau_decay_mode_2 == 1){
						if (fabs(otree->eta_2) < 1.5)
							shift_tes = shift_tes_lepfake_1p1p0_barrel; 
						else
							shift_tes = shift_tes_lepfake_1p1p0_endcap; 
					}  
      	}
				
	if (usePuppiMET) {
	  TLorentzVector tauLV_unclES_UP = tauLV;
	  TLorentzVector tauLV_unclES_DOWN = tauLV;
	  float puppiMET_Up = TMath::Sqrt(analysisTree.puppimet_ex_UnclusteredEnUp*analysisTree.puppimet_ex_UnclusteredEnUp+
					  analysisTree.puppimet_ey_UnclusteredEnUp*analysisTree.puppimet_ey_UnclusteredEnUp);
	  float puppiMET_Down = TMath::Sqrt(analysisTree.puppimet_ex_UnclusteredEnDown*analysisTree.puppimet_ex_UnclusteredEnDown+
					    analysisTree.puppimet_ey_UnclusteredEnDown*analysisTree.puppimet_ey_UnclusteredEnDown);
	  TLorentzVector puppiUncl_UP; puppiUncl_UP.SetXYZT(analysisTree.puppimet_ex_UnclusteredEnUp,
							    analysisTree.puppimet_ey_UnclusteredEnUp,
							    0,
							    puppiMET_Up); 
	  TLorentzVector puppiUncl_DOWN; puppiUncl_DOWN.SetXYZT(analysisTree.puppimet_ex_UnclusteredEnDown,
								analysisTree.puppimet_ey_UnclusteredEnDown,
								0,
								puppiMET_Down); 
	  //	  std::cout << "Central : " << puppimetLV.Pt()
	  //		    << "    Up : " << puppiUncl_UP.Pt()
	  //		    << "    Down : " << puppiUncl_DOWN.Pt() << std::endl;
	  correctTauES(tauLV, puppimetLV, shift_tes, isOneProng);
	  correctTauES(tauLV_unclES_UP,puppiUncl_UP,shift_tes, isOneProng);
	  correctTauES(tauLV_unclES_DOWN,puppiUncl_DOWN,shift_tes, isOneProng);
	  otree->pt_2 = tauLV.Pt();
	  otree->m_2 = tauLV.M();
	  otree->puppimet = puppimetLV.Pt();
	  otree->puppimetphi = puppimetLV.Phi();
	  otree->puppimet_ex_UnclusteredEnUp = puppiUncl_UP.Px();
	  otree->puppimet_ey_UnclusteredEnUp = puppiUncl_UP.Py();
	  otree->puppimet_ex_UnclusteredEnDown = puppiUncl_DOWN.Px();
	  otree->puppimet_ey_UnclusteredEnDown = puppiUncl_DOWN.Py();
	  //
	  //	  std::cout << "Corrected -> " << std::endl;
	  //	  std::cout << "Central : " << puppimetLV.Pt()
	  //                    << "    Up : " << puppiUncl_UP.Pt()
	  //                    << "    Down : " << puppiUncl_DOWN.Pt() << std::endl;
	  //
	}
	else {
	  correctTauES(tauLV,metLV,shift_tes, isOneProng);
	  otree->pt_2 = tauLV.Pt();
          otree->m_2 = tauLV.M();
	  otree->met = metLV.Pt();
	  otree->metphi = metLV.Phi();
	}
      }
    
      // rejecting events with pT(lep),pT(tau)<cut
      if (otree->pt_2<ptTauLowCut) continue;
      if (otree->pt_1<ptLeptonLowCut) continue;

      // if (!isData) {
      // 	if(otree->gen_match_2 == 5 && tauLV.E() <= 400 && tauLV.E() >= 20){
      // 	  if (otree->tau_decay_mode_2 == 0) tauLV *= (1-0.03);
      // 	  else if (otree->tau_decay_mode_2 < 5) tauLV *= (1-0.02);
      // 	  else if (otree->tau_decay_mode_2 == 10)tauLV *= (1-0.01);
      // 	  otree->pt_2 = tauLV.Pt();
      // 	  otree->m_2 = tauLV.M();
      // 	}
      // }
      
      ////////////////////////////////////////////////////////////
      // Filling variables (with corrected MET and tau momentum)
      ////////////////////////////////////////////////////////////

      if (otree->gen_match_2 == 5 && !isData)
    	  otree->tauvsjetweightMedium_2 = tauIDSF_medium->getSFvsPT(otree->pt_2);
      else 
        otree->tauvsjetweightMedium_2 = 1.;
      
      TLorentzVector dileptonLV = leptonLV + tauLV;
      otree->m_vis = dileptonLV.M();
      otree->pt_tt = (dileptonLV+metLV).Pt();   
      if (usePuppiMET)
	otree->pt_tt = (dileptonLV+puppimetLV).Pt();
    
      // mt TOT
      TLorentzVector metxLV = metLV;
      if (usePuppiMET) metxLV = puppimetLV;
      float mtTOT = 2*(otree->pt_1)*metxLV.Pt()*(1-cos(DeltaPhi(leptonLV,metxLV)));
      mtTOT += 2*(otree->pt_2)*metxLV.Pt()*(1-cos(DeltaPhi(tauLV,metxLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(leptonLV,tauLV))); 
      otree->mt_tot = TMath::Sqrt(mtTOT);
    
      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;
    
      // dilepton veto
      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
      if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree, era, isEmbedded);
    
  	  //extra lepton veto
      otree->extraelec_veto = extra_electron_veto(leptonIndex, ch, &cfg, &analysisTree, era, isEmbedded);
      otree->extramuon_veto = extra_muon_veto(leptonIndex, ch, &cfg, &analysisTree, isData);
    
      counter[13]++;
    
      otree->mt_1 = mT(leptonLV, metLV);
      otree->mt_2 = mT(tauLV, metLV);
      otree->puppimt_1 = mT(leptonLV, puppimetLV);
      otree->puppimt_2 = mT(tauLV, puppimetLV);
    
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
      otree->puppipzetamiss = otree->puppimet*TMath::Cos(otree->puppimetphi)*zetaX + otree->puppimet*TMath::Sin(otree->puppimetphi)*zetaY;
      counter[14]++;
    
      bool isSRevent = true; //boolean used to compute SVFit variables only on SR events, it is set to true when running Synchronization to run SVFit on all events
      if(!Synch){
	isSRevent = (otree->dilepton_veto<0.5 &&  otree->extramuon_veto<0.5 && otree->extraelec_veto<0.5 && otree->pt_1>19 && otree->pt_2>19 && otree->byVVVLooseDeepTau2017v2p1VSjet_2>0.5);
	if(usePuppiMET) isSRevent = isSRevent && otree->puppimt_1<60;
	else isSRevent = isSRevent && otree->mt_1<60;
  if(ch == "mt") isSRevent = isSRevent && (otree->trg_singlemuon>0.5 || otree->trg_mutaucross>0.5);
  if(ch == "et") isSRevent = isSRevent && (otree->trg_singleelectron>0.5 || otree->trg_etaucross>0.5);
      }
      /*
	if (!isSRevent) { 
	cout << "                                        " << endl;
	cout << "========================================" << endl;
	cout << "        Event is not selected           " << endl;
	cout << "========================================" << endl;
	cout << "                                        " << endl;
	}
	if (otree->byMediumDeepTau2017v2p1VSjet_2<0.5){
	
	auto args = std::vector<double>{otree->pt_2,
	static_cast<double>(otree->tau_decay_mode_2),
	static_cast<double>(otree->njets),
	otree->pt_1,
	static_cast<double>(otree->os),
	otree->puppimet,
	otree->puppimt_1,
	otree->iso_1,
	static_cast<double>(otree->trg_singlemuon),
	otree->m_vis};
	otree->ff_nom = fns_["ff_mt_medium_dmbins"]->eval(args.data());
	
	auto args_mva = std::vector<double>{otree->pt_2,
	otree->dmMVA_2,
	static_cast<double>(otree->njets),
	otree->pt_1,
	static_cast<double>(otree->os),
	otree->puppimet,
	otree->puppimt_1,
	otree->iso_1,
	static_cast<double>(otree->trg_singlemuon),
	otree->m_vis};
	otree->ff_mva = fns_["ff_mt_medium_mvadmbins"]->eval(args_mva.data());
	
	//		cout << "ff_nom : " << ff_nom << "   ff_mva : " << ff_mva << endl;
	
	}else { 
	otree->ff_nom = 1.;
	otree->ff_mva = 1.;
	}
	otree->ff_nom_sys = 0.15;
	otree->ff_mva_sys = 0.15;
      */

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

      // ***********************************
      // ** IPSignificance calibration ->
      // ***********************************
      otree->v_tracks = 0;
      std::vector<float> PV_covariance; PV_covariance.clear();
      std::vector<float> PVBS_covariance; PVBS_covariance.clear();
      // by default store non-refitted PV with BS constraint if refitted one is not found

      float vtx_x = analysisTree.primvertexwithbs_x; 
      float vtx_y = analysisTree.primvertexwithbs_y;
      float vtx_z = analysisTree.primvertexwithbs_z;

      float vtx_bs_x = analysisTree.primvertexwithbs_x; 
      float vtx_bs_y = analysisTree.primvertexwithbs_y;
      float vtx_bs_z = analysisTree.primvertexwithbs_z;

      for (int j = 0; j<6 ; ++j) {
	PV_covariance.push_back(analysisTree.primvertexwithbs_cov[j]);
	PVBS_covariance.push_back(analysisTree.primvertexwithbs_cov[j]);
      }

      for(unsigned int i = 0; i < analysisTree.refitvertexwithbs_count; i++)
        {
	  bool lep_match = false;
	  if (ch=="mt")
	    lep_match = leptonIndex == analysisTree.refitvertexwithbs_muIndex[i][0] || leptonIndex == analysisTree.refitvertexwithbs_muIndex[i][1];
	  else if (ch=="et")
	    lep_match = leptonIndex == analysisTree.refitvertexwithbs_eleIndex[i][0] || leptonIndex == analysisTree.refitvertexwithbs_eleIndex[i][1];

          if( lep_match &&
              (tauIndex == analysisTree.refitvertexwithbs_tauIndex[i][0] || tauIndex == analysisTree.refitvertexwithbs_tauIndex[i][1]))
            {
              otree->v_tracks = analysisTree.refitvertexwithbs_ntracks[i];
	      vtx_x = analysisTree.refitvertexwithbs_x[i];
	      vtx_y = analysisTree.refitvertexwithbs_y[i];
	      vtx_z = analysisTree.refitvertexwithbs_z[i];
	      for (int j=0; j<6; ++j) 
		PV_covariance[j] = analysisTree.refitvertexwithbs_cov[i][j];
            }
        }

      std::map<TString,IpCorrection*> ipCorrectorsNULL = {
	{"ipTau1",NULL},
	{"ipTau1BS",NULL},
	{"ipTau2",NULL},
	{"ipTau2BS",NULL}
      };

      ImpactParameter IP;
      std::map<TString,IpCorrection*> ipCorrectorsPass = ipCorrectorsNULL;
      if (ApplyIpCorrection && (!isData || isEmbedded)) ipCorrectorsPass = ipCorrectors;
      //      std::cout << "before..." << std::endl;
      acott_Impr(&analysisTree, otree, leptonIndex, tauIndex, ch, era, isEmbedded, ipCorrectorsPass);
      //      std::cout << "after..." << std::endl;

      otree->acotautau_00 = otree->acotautau_refitbs_00;
      otree->acotautau_01 = otree->acotautau_refitbs_01;
      
      // *****
      // RefitV + BS
      // *****
      TVector3 vertex(vtx_x,vtx_y,vtx_z);
      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1;
      TVector3 IP1;
      double ipsig1 = IP_significance_helix_lep(&analysisTree,leptonIndex, ch, vertex,PV_covariance,ipCov1,IP1);

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2;
      TVector3 IP2;
      double ipsig2 = IP_significance_helix_tauh(&analysisTree,tauIndex,vertex,PV_covariance,ipCov2,IP2);
      // ******
      // PV+BS
      // ******
      TVector3 vertex_bs(vtx_bs_x,vtx_bs_y,vtx_bs_z);
      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1_bs;
      TVector3 IP1_bs;
      double ipsig1_bs = IP_significance_helix_lep(&analysisTree,leptonIndex,ch,vertex_bs,PVBS_covariance,ipCov1_bs,IP1_bs);

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2_bs;
      TVector3 IP2_bs;
      double ipsig2_bs = IP_significance_helix_tauh(&analysisTree,tauIndex,vertex_bs,PVBS_covariance,ipCov2_bs,IP2_bs);

      /*
      cout << "ipsig1 = " << ipsig1 << " ipsig2 = " << ipsig2 << endl;
      cout << "IP1  : x = " << IP1.x() 
	   << "  y = " << IP1.y() 
	   << "  z = " << IP1.z() << std::endl;
      cout << "       x = " << otree->ipx_uncorr_1 
	   << "  y = " << otree->ipy_uncorr_1 
	   << "  z = " << otree->ipz_uncorr_1 << std::endl;

      cout << "IP2  : x = " << IP2.x() 
	   << "  y = " << IP2.y() 
	   << "  z = " << IP2.z() << std::endl;
      cout << "       x = " << otree->ipx_uncorr_2 
	   << "  y = " << otree->ipy_uncorr_2 
	   << "  z = " << otree->ipz_uncorr_2 << std::endl;
      */
      // Uncorrected values

      TLorentzVector ip1; ip1.SetXYZM(otree->ipx_uncorr_1,otree->ipy_uncorr_1,otree->ipz_uncorr_1,0.);
      otree->ipxy_uncorr_1 = ip1.Pt();
      otree->ipn_uncorr_1 = ip1.P();
      otree->drip_uncorr_1 = deltaR(otree->eta_1,otree->phi_1,ip1.Eta(),ip1.Phi());
      otree->detaip_uncorr_1 = ip1.Eta() - otree->eta_1; 
      TVector3 vectIP = TVector3(otree->ipx_uncorr_1,otree->ipy_uncorr_1,0.);
      TVector3 vectP  = TVector3(leptonLV.Px(),leptonLV.Py(),0.);
      otree->dphiip_uncorr_1 = TMath::ACos(vectIP*vectP/(vectIP.Mag()*vectP.Mag()));
      //      cout << "dphi = " << otree->dphiip_uncorr_1 << std::endl;
      //      cout << "refit = " << otree->isrefitBS << std::endl;

      TVector3 Ip1(otree->ipx_uncorr_1,otree->ipy_uncorr_1,otree->ipz_uncorr_1);
      otree->IP_signif_RefitV_with_BS_uncorr_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1);
      Ip1.SetXYZ(otree->ipx_bs_uncorr_1,otree->ipy_bs_uncorr_1,otree->ipz_bs_uncorr_1);
      otree->IP_signif_PV_with_BS_uncorr_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_bs);

      otree->ip_covxx_1 = ipCov1(0,0);
      otree->ip_covxy_1 = ipCov1(0,1);
      otree->ip_covxz_1 = ipCov1(0,2);
      otree->ip_covyy_1 = ipCov1(1,1);
      otree->ip_covyz_1 = ipCov1(1,2);
      otree->ip_covzz_1 = ipCov1(2,2);

      TLorentzVector ip2; ip2.SetXYZM(otree->ipx_uncorr_2,otree->ipy_uncorr_2,otree->ipz_uncorr_2,0.);
      otree->ipxy_uncorr_2 = ip2.Pt();
      otree->ipn_uncorr_2 = ip2.P();
      otree->drip_uncorr_2 = deltaR(otree->eta_2,otree->phi_2,ip2.Eta(),ip2.Phi());
      otree->detaip_uncorr_2 = ip2.Eta() - otree->eta_2;
      vectIP.SetX(otree->ipx_uncorr_2);
      vectIP.SetY(otree->ipy_uncorr_2);
      vectIP.SetZ(0.);
      vectP.SetX(tauLV.Px());
      vectP.SetY(tauLV.Py());
      vectP.SetZ(0.);
      otree->dphiip_uncorr_2 = TMath::ACos(vectIP*vectP/(vectIP.Mag()*vectP.Mag()));
      TVector3 Ip2(otree->ipx_uncorr_2,otree->ipy_uncorr_2,otree->ipz_uncorr_2);
      otree->IP_signif_RefitV_with_BS_uncorr_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2);
      Ip2.SetXYZ(otree->ipx_bs_uncorr_2,otree->ipy_bs_uncorr_2,otree->ipz_bs_uncorr_2);
      otree->IP_signif_PV_with_BS_uncorr_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_bs);

      otree->ip_covxx_2 = ipCov2(0,0);
      otree->ip_covxy_2 = ipCov2(0,1);
      otree->ip_covxz_2 = ipCov2(0,2);
      otree->ip_covyy_2 = ipCov2(1,1);
      otree->ip_covyz_2 = ipCov2(1,2);
      otree->ip_covzz_2 = ipCov2(2,2);

      // Corrected values 

      ip1.SetXYZM(otree->ipx_1,otree->ipy_1,otree->ipz_1,0.);
      otree->ipxy_1 = ip1.Pt();
      otree->ipn_1 = ip1.P();
      otree->drip_1 = deltaR(otree->eta_1,otree->phi_1,ip1.Eta(),ip1.Phi());
      otree->detaip_1 = ip1.Eta() - otree->eta_1; 
      vectIP.SetX(otree->ipx_1);
      vectIP.SetY(otree->ipy_1);
      vectIP.SetZ(0.);
      vectP.SetX(leptonLV.Px());
      vectP.SetY(leptonLV.Py());
      vectP.SetZ(0.);
      otree->dphiip_1 = TMath::ACos(vectIP*vectP/(vectIP.Mag()*vectP.Mag()));

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1_corr;
      if (otree->ipx_1<-9998)  {
	otree->IP_signif_RefitV_with_BS_1 = -9999;
      }
      else {
	ipCov1_corr =ipCov1;
	Ip1.SetXYZ(otree->ipx_1,otree->ipy_1,otree->ipz_1);
	if (ipCorrectorsPass["ipTau1"]!=NULL)
	  ipCov1_corr = ipCorrectorsPass["ipTau1"]->correctIpCov(ipCov1,otree->eta_1);      
	otree->IP_signif_RefitV_with_BS_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_corr);
      }

      if (otree->ipx_bs_1<-9998) {
	otree->IP_signif_PV_with_BS_1 = -9999;
      }
      else {
	Ip1.SetXYZ(otree->ipx_bs_1,otree->ipy_bs_1,otree->ipz_bs_1);
	ipCov1_corr = ipCov1_bs;
	if (ipCorrectorsPass["ipTaus1BS"]!=NULL)
	  ipCov1_corr = ipCorrectorsPass["ipTau1BS"]->correctIpCov(ipCov1_bs,otree->eta_1);      
	otree->IP_signif_PV_with_BS_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_corr);
      }

      ip2.SetXYZM(otree->ipx_2,otree->ipy_2,otree->ipz_2,0.);
      otree->ipxy_2 = ip2.Pt();
      otree->ipn_2 = ip2.P();
      otree->drip_2 = deltaR(otree->eta_2,otree->phi_2,ip2.Eta(),ip2.Phi());
      otree->detaip_2 = ip2.Eta() - otree->eta_2; 
      vectIP.SetX(otree->ipx_2);
      vectIP.SetY(otree->ipy_2);
      vectIP.SetZ(0.);
      vectP.SetX(tauLV.Px());
      vectP.SetY(tauLV.Py());
      vectP.SetZ(0.);
      otree->dphiip_2 = TMath::ACos(vectIP*vectP/(vectIP.Mag()*vectP.Mag()));

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2_corr;
      if (otree->ipx_2<-9998) {
	otree->IP_signif_RefitV_with_BS_2 = -9999;
      }
      else {
	Ip2.SetXYZ(otree->ipx_2,otree->ipy_2,otree->ipz_2);
	ipCov2_corr = ipCov2;
	if (ipCorrectorsPass["ipTau2"]!=NULL)
	  ipCov2_corr = ipCorrectorsPass["ipTau2"]->correctIpCov(ipCov2,otree->eta_2);      
	otree->IP_signif_RefitV_with_BS_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_corr);      
      }

      if (otree->ipx_bs_2<-9998) {
	otree->IP_signif_PV_with_BS_2 = -9999;
      }
      else {
	Ip2.SetXYZ(otree->ipx_bs_2,otree->ipy_bs_2,otree->ipz_bs_2);
	ipCov2_corr = ipCov2_bs;
	if (ipCorrectorsPass["ipTau2BS"]!=NULL)
	  ipCov2_corr = ipCorrectors["ipTau2BS"]->correctIpCov(ipCov2_bs,otree->eta_2);      
	otree->IP_signif_PV_with_BS_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_corr);      
      }

      otree->tauspinnerH = 0.;
      otree->tauspinnerA = 0.;
      otree->tauspinnerMaxMix = 0.;
      if (!std::isnan(otree->TauSpinnerWeightsEven)) 
	otree->tauspinnerH = otree->TauSpinnerWeightsEven;
      if (!std::isnan(otree->TauSpinnerWeightsOdd)) 
	otree->tauspinnerA = otree->TauSpinnerWeightsOdd;
      if (!std::isnan(otree->TauSpinnerWeightsMaxMix)) 
	otree->tauspinnerMaxMix = otree->TauSpinnerWeightsMaxMix;

      /*
      std::cout << "IPSig(1) = " << otree->IP_signif_RefitV_with_BS_uncorr_1
		<< "  :  " << otree->IP_signif_RefitV_with_BS_1 
		<< "  -> " << otree->IP_signif_PV_with_BS_1 << std::endl;

      std::cout << "IPSig(2) = " << otree->IP_signif_RefitV_with_BS_uncorr_2
		<< "  :  " << otree->IP_signif_RefitV_with_BS_2 
		<< "  -> " << otree->IP_signif_PV_with_BS_2 << std::endl;
      std::cout << std::endl;
      */

      otree->ip_sig_1 = otree->IP_signif_RefitV_with_BS_1;
      otree->ip_sig_2 = otree->IP_signif_RefitV_with_BS_2;
        
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

	if (ch == "mt") { 
	  muonScaleSys->Eval(utils::MUTAU);
	  tauOneProngScaleSys->Eval(utils::MUTAU);
	  tauOneProngOnePi0ScaleSys->Eval(utils::MUTAU);
	  tauThreeProngScaleSys->Eval(utils::MUTAU);
	  tauThreeProngOnePi0ScaleSys->Eval(utils::MUTAU);
	  if (isDY) {
	    lepTauFakeOneProngScaleSys->Eval(utils::MUTAU);
	    lepTauFakeOneProngOnePi0ScaleSys->Eval(utils::MUTAU);
	  }
	}
	else if (ch == "et") {
    electronScaleSys->SetElectronIndex(leptonIndex);
    electronScaleSys->SetIsEmbedded(isEmbedded);
    electronScaleSys->SetAC1B(&analysisTree);
    electronScaleSys->Eval(utils::ETAU);
    
	  tauOneProngScaleSys->Eval(utils::ETAU);
	  tauOneProngOnePi0ScaleSys->Eval(utils::ETAU);
	  tauThreeProngScaleSys->Eval(utils::ETAU);
	  tauThreeProngOnePi0ScaleSys->Eval(utils::ETAU);
	  if (isDY) {
	    lepTauFakeOneProngScaleSys->Eval(utils::ETAU);
	    lepTauFakeOneProngOnePi0ScaleSys->Eval(utils::ETAU);	    
	  }
	}
      }
      
      counter[19]++;  

      selEvents++;
      
      
      //      cout << "Once again puppimet central : " << otree->puppimet << endl;
      //Merijn 2019 1 10: perhaps this should be called before moving to next event..
      otree->Fill();
    } // event loop
    
    delete _treeTauSpinnerWeights;
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  } // file loop
   
  //
  //  std::cout << "COUNTERS" << std::endl;
  //  for(int iC = 0; iC < 20; iC++) std::cout << "Counter " << iC << ":    " << counter[iC] << std::endl;
  //
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

  if (muonScaleSys != 0) {
    muonScaleSys->Write("",TObject::kOverwrite);
    delete muonScaleSys;
  }
  
  if (electronScaleSys != 0) {
    electronScaleSys->Write("",TObject::kOverwrite);
    delete electronScaleSys;
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

  if(lepTauFakeOneProngScaleSys != 0){
    lepTauFakeOneProngScaleSys->Write("",TObject::kOverwrite);
    delete lepTauFakeOneProngScaleSys;
  }

  if(lepTauFakeOneProngOnePi0ScaleSys != 0){
    lepTauFakeOneProngOnePi0ScaleSys->Write("",TObject::kOverwrite);
    delete lepTauFakeOneProngOnePi0ScaleSys;
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
    else if (decay_mode >= 10) shift_tes = shift_tes_3prong;
  }
  else {
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
void FillMuTau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex, float dRiso){
	otree->pt_1  = 	analysisTree->muon_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->muon_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->muon_phi[leptonIndex];
  otree->m_1   =  muonMass;
  otree->q_1   =  analysisTree->muon_charge[leptonIndex];
  otree->gen_match_1 = analysisTree->muon_genmatch[leptonIndex];
  otree->iso_1 = abs_Iso_mt(leptonIndex, analysisTree, dRiso) / analysisTree->muon_pt[leptonIndex];
  otree->d0_1 = analysisTree->muon_dxy[leptonIndex];
  otree->dZ_1 = analysisTree->muon_dz[leptonIndex];
  otree->tau_decay_mode_1 = -9999; 
  otree->dm_1 = 16;
  otree->dmMVA_1 = 16;
  
  /*
  std::vector<float> PV_with_BS_cov_components = {};
  for(auto i:  analysisTree->primvertexwithbs_cov) PV_with_BS_cov_components.push_back(i);	
  TVector3 PV_with_BS (analysisTree->primvertexwithbs_x, analysisTree->primvertexwithbs_y, analysisTree->primvertexwithbs_z );
  TVector3 ip; 
  ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCovariance;
  otree->IP_signif_PV_with_BS_1 = IP_significance_helix_lep(analysisTree, leptonIndex, "mt", PV_with_BS, PV_with_BS_cov_components,ipCovariance,ip);
  */
  TLorentzVector muon_P4;
  muon_P4.SetXYZM(analysisTree->muon_px[leptonIndex], analysisTree->muon_py[leptonIndex], analysisTree->muon_pz[leptonIndex], muonMass);
  otree->chpt_1 = muon_P4.Pt();
  otree->cheta_1 = muon_P4.Eta();
  otree->chphi_1 = muon_P4.Phi();
  otree->chm_1 = muon_P4.M();
  
  otree->npt_1 = (muon_P4 - muon_P4).Pt();
  otree->neta_1 = (muon_P4 - muon_P4).Eta();
  otree->nphi_1 = (muon_P4 - muon_P4).Phi();
  otree->nm_1 = (muon_P4 - muon_P4).M();
    
  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -9999;
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = -9999;
  
  otree->deepTauVsEleRaw_1  = -9999;
  otree->deepTauVsJetRaw_1  = -9999;
  otree->deepTauVsMuRaw_1  = -9999;
  
  otree->againstMuonLoose3_1 = -9999;
  otree->againstMuonTight3_1 = -9999;
  otree->againstElectronLooseMVA6_1 = -9999;
  otree->againstElectronMediumMVA6_1 = -9999;
  otree->againstElectronTightMVA6_1 = -9999;
  otree->againstElectronVLooseMVA6_1 = -9999;
  otree->againstElectronVTightMVA6_1 = -9999;

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = -9999;
  otree->byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = -9999;
  otree->byIsolationMVA3oldDMwLTraw_1 = -9999;
  otree->chargedIsoPtSum_1 = -9999;
  otree->neutralIsoPtSum_1 = -9999;
  otree->puCorrPtSum_1 = -9999;
}

//fill the otree with the electron variables in channel etau
void FillETau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, float dRiso, int era, bool isEmbedded){
  
  float sf_eleES = 1.;
  if (isEmbedded) sf_eleES = EmbedElectronES_SF(analysisTree, era, leptonIndex);  

  otree->pt_1  = 	sf_eleES*analysisTree->electron_pt[leptonIndex];
  otree->pt_uncorr_1 =  analysisTree->electron_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->electron_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->electron_phi[leptonIndex];
  otree->m_1 = 		electronMass;
    otree->q_1 = -1;
  if (analysisTree->electron_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->electron_genmatch[leptonIndex];

  otree->iso_1 =  abs_Iso_et(leptonIndex, analysisTree, dRiso) / (sf_eleES*analysisTree->electron_pt[leptonIndex]);

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
void FillTau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->pt_uncorr_2 = analysisTree->tau_pt[tauIndex];
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
  otree->dm_2 = analysisTree->tau_decayMode[tauIndex];
  otree->dmMVA_2 = analysisTree->tau_MVADM2017v1[tauIndex];
  otree->DM = analysisTree->tau_decayModeFinding[tauIndex];

  /*
  std::vector<float> PV_with_BS_cov_components = {};
  for(auto i:  analysisTree->primvertexwithbs_cov) PV_with_BS_cov_components.push_back(i);	
  TVector3 PV_with_BS (analysisTree->primvertexwithbs_x, analysisTree->primvertexwithbs_y, analysisTree->primvertexwithbs_z );
  ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCovariance;
  TVector3 ip;
  otree->IP_signif_PV_with_BS_2 = IP_significance_helix_tauh(analysisTree, tauIndex, PV_with_BS, PV_with_BS_cov_components, ipCovariance, ip);
  */

  TLorentzVector constituents_P4 = charged_constituents_P4(analysisTree, tauIndex);
  TLorentzVector tau_P4;
  tau_P4.SetXYZM(analysisTree->tau_px[tauIndex],
                 analysisTree->tau_py[tauIndex],
                 analysisTree->tau_pz[tauIndex],
                 analysisTree->tau_mass[tauIndex]);
  otree->chpt_2 = constituents_P4.Pt();
  otree->cheta_2 = constituents_P4.Eta();
  otree->chphi_2 = constituents_P4.Phi();
  otree->chm_2 = constituents_P4.M();
  
  otree->npt_2 = (tau_P4 - constituents_P4).Pt();
  otree->neta_2 = (tau_P4 - constituents_P4).Eta();
  otree->nphi_2 = (tau_P4 - constituents_P4).Phi();
  otree->nm_2 = (tau_P4 - constituents_P4).M();
  
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

  otree->deepTauVsEleRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSeraw[tauIndex];
  otree->deepTauVsJetRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSjetraw[tauIndex];
  otree->deepTauVsMuRaw_2                 = analysisTree->tau_byDeepTau2017v2p1VSmuraw[tauIndex];
  otree->byLooseDeepTau2017v2p1VSe_2      = analysisTree->tau_byLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byLooseDeepTau2017v2p1VSjet_2    = analysisTree->tau_byLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byLooseDeepTau2017v2p1VSmu_2     = analysisTree->tau_byLooseDeepTau2017v2p1VSmu[tauIndex];
  otree->byMediumDeepTau2017v2p1VSe_2     = analysisTree->tau_byMediumDeepTau2017v2p1VSe[tauIndex];
  otree->byMediumDeepTau2017v2p1VSjet_2   = analysisTree->tau_byMediumDeepTau2017v2p1VSjet[tauIndex];
  otree->byMediumDeepTau2017v2p1VSmu_2    = analysisTree->tau_byMediumDeepTau2017v2p1VSmu[tauIndex];
  otree->byTightDeepTau2017v2p1VSe_2      = analysisTree->tau_byTightDeepTau2017v2p1VSe[tauIndex];
  otree->byTightDeepTau2017v2p1VSjet_2    = analysisTree->tau_byTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byTightDeepTau2017v2p1VSmu_2     = analysisTree->tau_byTightDeepTau2017v2p1VSmu[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSe_2     = analysisTree->tau_byVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSmu_2    = analysisTree->tau_byVLooseDeepTau2017v2p1VSmu[tauIndex];
  otree->byVTightDeepTau2017v2p1VSe_2     = analysisTree->tau_byVTightDeepTau2017v2p1VSe[tauIndex];
  otree->byVTightDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVLooseDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVVLooseDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVTightDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVTightDeepTau2017v2p1VSe[tauIndex];
  otree->byVVTightDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVVLooseDeepTau2017v2p1VSe_2   = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVVVLooseDeepTau2017v2p1VSjet_2 = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSjet[tauIndex];

  otree->MVADM2017v1DM0raw_2 = analysisTree->tau_MVADM2017v1DM0raw[tauIndex];
  otree->MVADM2017v1DM10raw_2 = analysisTree->tau_MVADM2017v1DM10raw[tauIndex];
  otree->MVADM2017v1DM11raw_2 = analysisTree->tau_MVADM2017v1DM11raw[tauIndex];
  otree->MVADM2017v1DM1raw_2 = analysisTree->tau_MVADM2017v1DM1raw[tauIndex];
  otree->MVADM2017v1DM2raw_2 = analysisTree->tau_MVADM2017v1DM2raw[tauIndex];
  otree->MVADM2017v1DMotherraw_2 = analysisTree->tau_MVADM2017v1DMotherraw[tauIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstMuonLoose3_2 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_2 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronLooseMVA6_2 = analysisTree->tau_againstElectronLooseMVA6[tauIndex];
  otree->againstElectronMediumMVA6_2 = analysisTree->tau_againstElectronMediumMVA6[tauIndex];
  otree->againstElectronTightMVA6_2 = analysisTree->tau_againstElectronTightMVA6[tauIndex];
  otree->againstElectronVLooseMVA6_2 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronVTightMVA6_2 = analysisTree->tau_againstElectronVTightMVA6[tauIndex];

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex];
  otree->byIsolationMVA3oldDMwLTraw_2 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex];
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
