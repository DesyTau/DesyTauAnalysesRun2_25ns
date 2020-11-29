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
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"

#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
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
void FillVertices(const AC1B * analysisTree,Synch17Tree *otree, const bool isData);
void FillGenTree(const AC1B * analysisTree, Synch17GenTree *gentree);
float getEmbeddedWeight(const AC1B * analysisTree, RooWorkspace* WS);
void FillElMu(const AC1B *analysisTree, Synch17Tree *otree, int electronIndex, float dRisoElectron, int muonIndex, float dRIsoMuon, int era, bool isEmbedded);

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

void selectMuonElePair(AC1B *analysisTree, vector<int> muons, vector<int> electrons, bool isMuonIsoR03, bool isElectronIsoR03, float dRleptonsCut, float ptMuonHighCut, float ptElectronHighCut, int &electronIndex, int &muonIndex, float &isoMuMin, float &isoEleMin, int era, bool isEmbedded){
         
   for (unsigned int im=0; im<muons.size(); ++im) {
      unsigned int mIndex  = muons.at(im);
      float neutralHadIsoMu = analysisTree->muon_neutralHadIso[mIndex];
      float photonIsoMu = analysisTree->muon_photonIso[mIndex];
      float chargedHadIsoMu = analysisTree->muon_chargedHadIso[mIndex];
      float puIsoMu = analysisTree->muon_puIso[mIndex];
      if (isMuonIsoR03) {
         neutralHadIsoMu = analysisTree->muon_r03_sumNeutralHadronEt[mIndex];
         photonIsoMu = analysisTree->muon_r03_sumPhotonEt[mIndex];
         chargedHadIsoMu = analysisTree->muon_r03_sumChargedHadronPt[mIndex];
         puIsoMu = analysisTree->muon_r03_sumPUPt[mIndex];
      }
      float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
      neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
      float absIsoMu = chargedHadIsoMu + neutralIsoMu;
      float relIsoMu = absIsoMu/analysisTree->muon_pt[mIndex];
      
      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
         unsigned int eIndex = electrons.at(ie);
         float dR = deltaR(analysisTree->electron_eta[eIndex],analysisTree->electron_phi[eIndex],
                           analysisTree->muon_eta[mIndex],analysisTree->muon_phi[mIndex]);
         
         if (dR<dRleptonsCut) continue;
         
         double ele_sf = 1.0;
         if (isEmbedded) ele_sf= EmbedElectronES_SF(analysisTree, era, electronIndex );

         bool trigMatch =
            (analysisTree->muon_pt[mIndex]>ptMuonHighCut) ||
            (analysisTree->electron_pt[eIndex] * ele_sf >ptElectronHighCut);
         if (!trigMatch) continue;
         
         float absIsoEle; 
         float relIsoEle;
         float rhoNeutral = analysisTree->rho;
         float  eA = getEffectiveArea( fabs(analysisTree->electron_superclusterEta[eIndex]) );
         absIsoEle = analysisTree->electron_r03_sumChargedHadronPt[eIndex] +
            TMath::Max(0.0f,analysisTree->electron_r03_sumNeutralHadronEt[eIndex]+analysisTree->electron_r03_sumPhotonEt[eIndex]-eA*rhoNeutral);
         relIsoEle = absIsoEle/(analysisTree->electron_pt[eIndex] *ele_sf);
            //}
         if (int(mIndex)!=muonIndex) {
            if (relIsoMu==isoMuMin) {
               if (analysisTree->muon_pt[mIndex]>analysisTree->muon_pt[muonIndex]) {
                  isoMuMin  = relIsoMu;
                  muonIndex = int(mIndex);
                  isoEleMin = relIsoEle;
                  electronIndex = int(eIndex);
               }
            }
            else if (relIsoMu<isoMuMin) {
               isoMuMin  = relIsoMu;
               muonIndex = int(mIndex);
               isoEleMin = relIsoEle;
               electronIndex = int(eIndex);
            }
         }
         else {
            if (relIsoEle==isoEleMin) {
               if (analysisTree->electron_pt[eIndex]>analysisTree->electron_pt[electronIndex]) {
                  isoEleMin = relIsoEle;
                  electronIndex = int(eIndex);
               }
            }
            else if (relIsoEle<isoEleMin) {
               isoEleMin = relIsoEle;
               electronIndex = int(eIndex);
            }
         }
      }
   }
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


// Synch ntuple producer in the e+mu channel

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

  if(argc < 3){
    std::cout << "RUN ERROR: wrong number of arguments"<< std::endl;
    std::cout << "Please run the code in the following way:"<< std::endl;
    std::cout << "SynchNTupleProducer_em_Run2 NameOfTheConfigurationFile FileList" << std::endl;
    std::cout << "example: SynchNTupleProducer_Run2 analysisMacroSynch_em_2018.conf DATA_SingleMuon" << std::endl;
    exit(-1);
  }

  // **** configuration analysis  
  Config cfg(argv[1]);

  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase) + "/src/" + TString(json_name)).Data(), json);
  }

  const int era = cfg.get<int>("era");
  const bool synch            = cfg.get<bool>("Synch");
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
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

  const string ipCorrFileNameMuon   = cfg.get<string>("IpCorrFileNameMuon");
  const string ipCorrFileNameMuonBS = cfg.get<string>("IpCorrFileNameMuonBS");
  const string ipCorrFileNameElec     = cfg.get<string>("IpCorrFileNameElec");
  const string ipCorrFileNameElecBS   = cfg.get<string>("IpCorrFileNameElecBS");

  TString IpCorrFileNameMuon(ipCorrFileNameMuon);
  TString IpCorrFileNameMuonBS(ipCorrFileNameMuonBS);
  TString IpCorrFileNameElec(ipCorrFileNameElec);
  TString IpCorrFileNameElecBS(ipCorrFileNameElecBS);

  std::string year_label;
  if (era == 2016) year_label = "2016Legacy";
  else if (era == 2017) year_label = "2017ReReco";
  else if (era == 2018) year_label = "2018ReReco";	
  else {std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n'; exit(-1);}

  IpCorrection *CorrectorIpMuon   = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameMuon);
  IpCorrection *CorrectorIpMuonBS = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameMuonBS);  
  IpCorrection *CorrectorIpElec   = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameElec);
  IpCorrection *CorrectorIpElecBS = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileNameElecBS);  

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

  std::cout << "btagging loaded" << std::endl;

  // MET Recoil Corrections
  const bool isDY = (infiles.find("DY") != string::npos) || (infiles.find("EWKZ") != string::npos);//Corrections that should be applied on EWKZ are the same needed for DY
  const bool isWJets = (infiles.find("WJets") != string::npos) || (infiles.find("W1Jets") != string::npos) || (infiles.find("W2Jets") != string::npos) || (infiles.find("W3Jets") != string::npos) || (infiles.find("W4Jets") != string::npos) || (infiles.find("EWK") != string::npos);
  const bool isHiggs = (infiles.find("VBFHTo")!= string::npos) || (infiles.find("WminusHTo")!= string::npos) || (infiles.find("WplusHTo")!= string::npos) || (infiles.find("ZHTo")!= string::npos) || (infiles.find("GluGluHTo")!= string::npos);
  const bool isEWKZ =  infiles.find("EWKZ") != string::npos;
  const bool isMG = infiles.find("madgraph") != string::npos;
  const bool isMSSMsignal =  (infiles.find("SUSYGluGluToHToTauTau")!= string::npos) || (infiles.find("SUSYGluGluToBBHToTauTau")!= string::npos);
  const bool isTTbar = infiles.find("TT") != string::npos;

  const bool isEmbedded = infiles.find("Embed") != string::npos;

  const bool ApplyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections") && !isEmbedded && !isData && (isDY || isWJets || isHiggs || isMSSMsignal);
  kit::RecoilCorrector PFMetRecoilCorrector(cfg.get<string>("PFMetRecoilFilePath"));
  kit::RecoilCorrector recoilCorrector(cfg.get<string>("RecoilFilePath"));
  kit::MEtSys MetSys(cfg.get<string>("RecoilSysFilePath"));

  std::cout << "recoil corrections" << std::endl;

  // kinematic cuts on electrons
  const float ptElectronLowCut    = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut   = cfg.get<float>("ptElectronHighCut");
  const float ptElectronSingleCut = cfg.get<float>("ptElectronSingleCut");
  const float etaElectronCut      = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut      = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut       = cfg.get<float>("dzElectronCut");
  const string lowPtLegElectron   = cfg.get<string>("LowPtLegElectron");
  const string highPtLegElectron  = cfg.get<string>("HighPtLegElectron");
  const vector<string> singleLegElectron = cfg.get<vector<string> >("SingleLegElectron");
  
  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  
  // kinematic cuts on muons
  const float ptMuonLowCut    = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut   = cfg.get<float>("ptMuonHighCut");
  const float ptMuonSingleCut = cfg.get<float>("ptMuonSingleCut");
  const float etaMuonCut      = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut      = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut       = cfg.get<float>("dzMuonCut");
  const string lowPtLegMuon   = cfg.get<string>("LowPtLegMuon");
  const string highPtLegMuon  = cfg.get<string>("HighPtLegMuon");
  const vector<string> singleLegMuon = cfg.get<vector<string> >("SingleLegMuon");
   
  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  
  // dR trigger, leptons
  const float dRTrigMatch = cfg.get<float>("dRTrigMatch");
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const bool isMuonIsoR03 = cfg.get<bool>("IsMuonIsoR03");
  const bool isElectronIsoR03 = cfg.get<bool>("IsElectronIsoR03");

  // dz filter
  const bool applyDzFilterMatch = cfg.get<bool>("ApplyDzFilterMatch");
  const string mu23ele12DzFilter = cfg.get<string>("Mu23Ele12DzFilter");
  const string mu8ele23DzFilter = cfg.get<string>("Mu8Ele23DzFilter");

  float dRIsoMuon = 0.4;
  if (isMuonIsoR03) dRIsoMuon = 0.3;
  float dRIsoElectron = 0.4;
  if (isElectronIsoR03) dRIsoElectron = 0.3;

  TString LowPtLegMuon(lowPtLegMuon);
  TString LowPtLegElectron(lowPtLegElectron);
  TString HighPtLegMuon(highPtLegMuon);
  TString HighPtLegElectron(highPtLegElectron);
  vector<TString> SingleLegMuon;
  vector<TString> SingleLegElectron;
  for (unsigned int i=0; i<singleLegMuon.size(); ++i) 
    SingleLegMuon.push_back(TString(singleLegMuon.at(i)));
  for (unsigned int i=0; i<singleLegElectron.size(); ++i) 
    SingleLegElectron.push_back(TString(singleLegElectron.at(i)));
  TString Mu23Ele12DzFilter(mu23ele12DzFilter);
  TString Mu8Ele23DzFilter(mu8ele23DzFilter);

  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");

  // correction workspace
  const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");
  const string CorrectionWorkspaceFileNameKIT = cfg.get<string>("CorrectionWorkspaceFileNameKIT");

  bool triggerEmbed2017 = false;
  float ptTriggerEmbed2017 = 40;
  float etaTriggerEmbed2017 = 1.479;

  if (era==2017&&isEmbedded)
    triggerEmbed2017 = true;

  // **** end of configuration analysis

  std::cout << "end of configuration " << std::endl;

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

  std::cout << "Number of files" << jfile << std::endl;

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

  std::cout << "PU loaded" << std::endl;

  // Workspace with corrections
  TString workspace_filename = TString(cmsswBase) + "/src/" + CorrectionWorkspaceFileName;
  TString workspace_filename_kit = TString(cmsswBase) + "/src/" + CorrectionWorkspaceFileNameKIT;
  cout << "Taking correction workspace from " << workspace_filename << endl;
  TFile *f_workspace = new TFile(workspace_filename, "read");
  TFile *f_workspace_kit = new TFile(workspace_filename_kit, "read");
  if (f_workspace->IsZombie()) {
    std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
     exit(-1);
   }
  if (f_workspace_kit->IsZombie()) {
    std::cout << " workspace file " << workspace_filename_kit << " not found. Please check. " << std::endl;
     exit(-1);
   }

  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");
  RooWorkspace *correctionWS = (RooWorkspace*)f_workspace_kit->Get("w");

  // Zpt reweighting for LO DY samples 
  TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
  TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

  std::cout << "ZPt weights" << std::endl;

  // load QCD =========================
  const string qcdFileName = cfg.get<string>("QCDFileName");
  TFile *fQCD = new TFile(TString(cmsswBase)+"/src/"+qcdFileName);

  TF1 *OS_SS_njetgt1 = (TF1*)fQCD->Get("OS_SS_transfer_factors_njetgt1");
  TF1 *OS_SS_njet1 = (TF1*)fQCD->Get("OS_SS_transfer_factors_njet1");
  TF1 *OS_SS_njet0 = (TF1*)fQCD->Get("OS_SS_transfer_factors_njet0");
  
  TGraph *OS_SS_njet0_Par0_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njet0_UP");
  TGraph *OS_SS_njet0_Par0_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njet0_DOWN");
  TGraph *OS_SS_njet0_Par1_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njet0_UP");
  TGraph *OS_SS_njet0_Par1_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njet0_DOWN");  
  TGraph *OS_SS_njet0_Par2_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njet0_UP");
  TGraph *OS_SS_njet0_Par2_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njet0_DOWN");  
  
  TGraph *OS_SS_njet1_Par0_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njet1_UP");
  TGraph *OS_SS_njet1_Par0_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njet1_DOWN");
  TGraph *OS_SS_njet1_Par1_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njet1_UP");
  TGraph *OS_SS_njet1_Par1_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njet1_DOWN");  
  TGraph *OS_SS_njet1_Par2_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njet1_UP");
  TGraph *OS_SS_njet1_Par2_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njet1_DOWN");  
  
  TGraph *OS_SS_njetgt1_Par0_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njetgt1_UP");
  TGraph *OS_SS_njetgt1_Par0_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par0_njetgt1_DOWN");
  TGraph *OS_SS_njetgt1_Par1_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njetgt1_UP");
  TGraph *OS_SS_njetgt1_Par1_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par1_njetgt1_DOWN");  
  TGraph *OS_SS_njetgt1_Par2_UP = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njetgt1_UP");
  TGraph *OS_SS_njetgt1_Par2_DOWN = (TGraph*)fQCD->Get("OS_SS_transfer_factors_Par2_njetgt1_DOWN");  
   
  TH2F *hNonClosureCorrection = (TH2F*)fQCD->Get("NonClosureCorrection");
  TH2F *hIsolationCorrection = (TH2F*)fQCD->Get("IsolationCorrection");

  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_em_Sync.root";
    
  std::cout <<rootFileName <<std::endl;  

  TFile *file = new TFile(rootFileName, "recreate");
  file->cd("");

  TH1D *inputEventsH = new TH1D("inputEventsH", "", 1, -0.5, 0.5);
  TH1D *nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5, 0.5);
  
  TTree *tree = new TTree("TauCheck", "TauCheck");
  TTree *gtree = new TTree("GenTauCheck", "GenTauCheck");
  Synch17Tree *otree = new Synch17Tree(tree,true);
  Synch17GenTree *gentree = new Synch17GenTree(gtree);
    
  int nTotalFiles = 0;
  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  //svFit
  TH1::AddDirectory(false);  
  TFile *inputFile_visPtResolution = new TFile(svFitPtResFile.data());

  std::cout << "inputFile_visPtResolution : " << std::endl;

  //Systematics init
  
  MuonScaleSys *muonScaleSys = 0;
  ElectronScaleSys *electronScaleSys = 0;

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

  ///////////////FILE LOOP///////////////

  for (int iF = ifile; iF < jfile; ++iF) {
    std::cout << "file " << iF + 1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    
    TFile *file_ = TFile::Open(fileList[iF].data());
    TTree *_tree = (TTree*)file_->Get(TString(ntupleName));  
    if (_tree == NULL) {
      std::cout << "TTree " << ntupleName << " is absent" << std::endl;
      continue;
    }
    
    TH1D *histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents == NULL) continue;
    int NE = int(histoInputEvents->GetEntries());
    std::cout << "      number of input events    = " << NE << std::endl;
    for (int iE = 0; iE < NE; ++iE)
      inputEventsH->Fill(0.);

    TTree * _inittree = (TTree*)file_->Get(TString(initNtupleName));
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

    //    AC1B analysisTree(_tree, isData);
    AC1B analysisTree(_tree);
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
    
    ///////////////EVENT LOOP///////////////
    Long64_t numberOfEntries = analysisTree.GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;    

    for (Long64_t iEntry = 0; iEntry < numberOfEntries; iEntry++) {

      if (nEvents % 10000 == 0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      analysisTree.GetEntry(iEntry);
      nEvents++;

      // filling generator tree
      if (!isData){
	initializeGenTree(gentree);
      	FillGenTree(&analysisTree,gentree);
      	gentree->Fill();
      }

      //Skip events not passing the MET filters, if applied
      bool passed_all_met_filters = passedAllMetFilters(&analysisTree, met_filters_list);
      if (ApplyMetFilters && !synch && !passed_all_met_filters) continue;
      otree->passedAllMetFilters = passed_all_met_filters;
      
      // accessing trigger info ====
      vector<bool> isSingleLegMuon; isSingleLegMuon.clear();
      vector<bool> isSingleLegElectron; isSingleLegElectron.clear();
      vector<unsigned int> nSingleLegMuon; nSingleLegMuon.clear();
      vector<unsigned int> nSingleLegElectron; nSingleLegElectron.clear();

      bool isMu23Ele12DzFilter = false;
      bool isMu8Ele23DzFilter = false;
      unsigned int nMu23Ele12DzFilter = 0;
      unsigned int nMu8Ele23DzFilter = 0;

      unsigned int nLowPtLegElectron = 0;
      unsigned int nHighPtLegElectron = 0;
      unsigned int nLowPtLegMuon = 0;
      unsigned int nHighPtLegMuon = 0;
      bool isLowPtLegElectron = accessTriggerInfo(&analysisTree,LowPtLegElectron,nLowPtLegElectron);
      bool isHighPtLegElectron = accessTriggerInfo(&analysisTree,HighPtLegElectron,nHighPtLegElectron);
      bool isLowPtLegMuon = accessTriggerInfo(&analysisTree,LowPtLegMuon,nLowPtLegMuon);
      bool isHighPtLegMuon = accessTriggerInfo(&analysisTree,HighPtLegMuon,nHighPtLegMuon);
      for (unsigned int i=0; i<SingleLegMuon.size(); ++i) {
	unsigned int nfilter = 0;
	bool isfilter = accessTriggerInfo(&analysisTree,SingleLegMuon.at(i),nfilter);
	isSingleLegMuon.push_back(isfilter);
	nSingleLegMuon.push_back(nfilter);
      }
      for (unsigned int i=0; i<SingleLegElectron.size(); ++i) {
	unsigned int nfilter = 0;
	bool isfilter = accessTriggerInfo(&analysisTree,SingleLegElectron.at(i),nfilter);
	isSingleLegElectron.push_back(isfilter);
	nSingleLegElectron.push_back(nfilter);
      }
      if (applyDzFilterMatch){
	isMu23Ele12DzFilter = accessTriggerInfo(&analysisTree,Mu23Ele12DzFilter,nMu23Ele12DzFilter);
	isMu8Ele23DzFilter = accessTriggerInfo(&analysisTree,Mu8Ele23DzFilter,nMu8Ele23DzFilter);
      }
    
      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
    
      // lumi filter
      if ((isData || isEmbedded) && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;
        
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

      // PU weight
      if (ApplyPUweight) 
        otree->puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));

      // generator weight
      if(!isData || isEmbedded){
        otree->mcweight = analysisTree.genweight;
        otree->gen_noutgoing = analysisTree.genparticles_noutgoing;
	if (isEmbedded&&otree->mcweight>1.0)
	  otree->mcweight = 0.0;
      }

      // selecting electrons      
      float sf_eleES = 1.0;   
      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	bool electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie];        
	if (isEmbedded) sf_eleES = EmbedElectronES_SF(&analysisTree, era, ie);          
	if (sf_eleES*analysisTree.electron_pt[ie] <= ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie]) >= etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie]) >= dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie]) >= dzElectronCut) continue;
	if (!electronMvaId) continue;
	if (!analysisTree.electron_pass_conversion[ie]) continue;
	if (analysisTree.electron_nmissinginnerhits[ie] > 1) continue;
	electrons.push_back(ie);
      }
    
      // selecting muons
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im < analysisTree.muon_count; ++im) {
	bool muonMediumId = isIdentifiedMediumMuon(im, &analysisTree, isData);	          
	if (analysisTree.muon_pt[im] <= ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im]) >= etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im]) >= dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im]) >= dzMuonCut) continue;
	if (!muonMediumId) continue;
	muons.push_back(im);
      }
    
      //      std::cout << "muons = " << muons.size() << std::endl;
      //      std::cout << "electrons = " << electrons.size() << std::endl;

      if (muons.size() == 0) continue;
      if (electrons.size() == 0) continue;
    
      int electronIndex = -1;
      int muonIndex = -1;
      float isoMuMin = 1e+10;
      float isoEleMin = 1e+10;
      bool isMuonIsoR03 = false;
      bool isElectronR03 = true;
      selectMuonElePair(&analysisTree, muons, electrons, isMuonIsoR03, isElectronIsoR03, dRleptonsCut, ptMuonHighCut, ptElectronHighCut, electronIndex, muonIndex, isoMuMin, isoEleMin, era, isEmbedded);
      if (electronIndex<0) continue;
      if (muonIndex<0) continue;
      
      // Filling ntuple with the electron and muon information
      FillElMu(&analysisTree,otree,electronIndex,dRIsoElectron,muonIndex,dRIsoMuon,era,isEmbedded);

      //      std::cout << "fill emu " << std::endl;
      //      std::cout << "e : pt = " << otree->pt_1 << "  eta = " << otree->eta_1 << "  phi = " << otree->phi_1 << std::endl;
      //      std::cout << "m : pt = " << otree->pt_2 << "  eta = " << otree->eta_2 << "  phi = " << otree->phi_2 << std::endl;

      //all criterua passed, we fill vertices here;	
      FillVertices(&analysisTree, otree, isData);
      
      //      std::cout << "fill vertices" << std::endl;

      /////////////////////////////
      // Trigger matching
      /////////////////////////////
      vector<bool> isSingleMuonMatch(isSingleLegMuon.size(), false);
      vector<bool> isSingleElectronMatch(isSingleLegElectron.size(), false);
      bool isHighPtMuonMatch = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isHighPtLegMuon,nHighPtLegMuon,dRTrigMatch);
      bool isLowPtMuonMatch  = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isLowPtLegMuon,nLowPtLegMuon,dRTrigMatch);

      bool isHighPtElectronMatch = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isHighPtLegElectron,nHighPtLegElectron,dRTrigMatch);;
      bool isLowPtElectronMatch  = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isLowPtLegElectron,nLowPtLegElectron,dRTrigMatch);

      for ( unsigned int i=0; i<isSingleLegMuon.size(); ++i) {
	isSingleMuonMatch.at(i) = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isSingleLegMuon.at(i),nSingleLegMuon.at(i),dRTrigMatch);
      }
      for ( unsigned int i=0; i<isSingleLegElectron.size(); ++i) {
	isSingleElectronMatch.at(i) = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isSingleLegElectron.at(i),nSingleLegElectron.at(i),dRTrigMatch);
      }
      if (era == 2017) {
        int id_SingleEGO = -1;
        int id_Single32 = -1;    
        for(unsigned int i_trig = 0; i_trig < SingleLegElectron.size(); i_trig++){
          if(SingleLegElectron.at(i_trig) == "hltEGL1SingleEGOrFilter")
            id_SingleEGO = i_trig;
          if(SingleLegElectron.at(i_trig) == "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")
            id_Single32 = i_trig;
        }
        isSingleElectronMatch[id_SingleEGO] = isSingleElectronMatch[id_SingleEGO] && isSingleElectronMatch[id_Single32];
        isSingleElectronMatch[id_Single32] =  isSingleElectronMatch[id_SingleEGO] && isSingleElectronMatch[id_Single32];
      }
      bool isHighPtMuonDZMatch = true;
      bool isLowPtMuonDZMatch = true;
      bool isHighPtElectronDZMatch = true;
      bool isLowPtElectronDZMatch = true;

      if (applyDzFilterMatch) {
	isLowPtMuonDZMatch  = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isMu8Ele23DzFilter,nMu8Ele23DzFilter,dRTrigMatch);
	isHighPtMuonDZMatch = triggerMatching(&analysisTree,otree->eta_2,otree->phi_2,isMu23Ele12DzFilter,nMu23Ele12DzFilter,dRTrigMatch); 
	isLowPtElectronDZMatch = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isMu23Ele12DzFilter,nMu23Ele12DzFilter,dRTrigMatch);
	isHighPtElectronDZMatch = triggerMatching(&analysisTree,otree->eta_1,otree->phi_1,isMu8Ele23DzFilter,nMu8Ele23DzFilter,dRTrigMatch);
      }

      otree->trg_singlemuon = false;
      otree->trg_singleelectron = false;
      otree->singleLepTrigger = false;

      // irrelevant for emu channel
      otree->trg_doubletau = false;
      otree->trg_mutaucross = false;
      otree->trg_mutaucross_mu = false;
      otree->trg_mutaucross_tau = false;
      otree->trg_etaucross = false;
      otree->trg_etaucross_e = false;
      otree->trg_etaucross_tau = false;

      // emu triggers (add them to tree)
      otree->trg_mulow = isLowPtMuonMatch;
      otree->trg_muhigh = isHighPtMuonMatch;
      otree->trg_elow = isLowPtElectronMatch;
      otree->trg_ehigh = isHighPtElectronMatch;
      otree->trg_muhigh_elow = isHighPtMuonMatch && isLowPtElectronMatch && isHighPtMuonDZMatch && isLowPtElectronDZMatch;
      otree->trg_ehigh_mulow = isHighPtElectronMatch && isLowPtMuonMatch && isLowPtMuonDZMatch && isHighPtElectronDZMatch;

      for(unsigned int i_trig = 0; i_trig < SingleLegElectron.size(); i_trig++)
        otree->trg_singleelectron = otree->trg_singleelectron || isSingleElectronMatch.at(i_trig);
      for(unsigned int i_trig = 0; i_trig < SingleLegMuon.size(); i_trig++)
        otree->trg_singlemuon = otree->trg_singlemuon || isSingleMuonMatch.at(i_trig);
      if (triggerEmbed2017) {
	if (otree->pt_1<ptTriggerEmbed2017&&fabs(otree->eta_1)>etaTriggerEmbed2017) {
	  otree->trg_singleelectron = true;
	}
      }    
      otree->singleLepTrigger = otree->trg_singleelectron || otree->trg_singlemuon;
    
      //      std::cout << "after trigger" << std::endl;

      // initialize JER (including data and embedded) 
      otree->apply_recoil = ApplyRecoilCorrections;
      jets::initializeJER(&analysisTree);

      //      std::cout << "initialising JER" << std::endl;

      if (!isData && !isEmbedded) { // JER smearing
	jets::associateRecoAndGenJets(&analysisTree, resolution);
	jets::smear_jets(&analysisTree,resolution,resolution_sf,true);
      }

      //      std::cout << "smeared jets" << std::endl;

      //counting jet
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);
  
      //      std::cout << "counted jets" << std::endl;

      ////////////////////////////////////////////////////////////
      // ID/Iso and Trigger Corrections
      ////////////////////////////////////////////////////////////

      // setting weights to 1
      otree->trkeffweight_1 = 1;
      otree->trkeffweight_2 = 1;
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

      otree->weight = 1;
      otree->weightSingle = 1;
      otree->weightExcl = 1;

      otree->trigweight_l_lt = 1;
      otree->trigweight_t_lt = 1;
	
      float eff_data_trig_m = 1;
      float eff_mc_trig_m = 1;
      float sf_trig_m = 1;
      float eff_data_trig_e = 1;
      float eff_mc_trig_e = 1;
      float sf_trig_e = 1;
      
      float eff_data_trig_mhigh = 1;
      float eff_mc_trig_mhigh = 1;
      float eff_data_trig_mlow = 1;
      float eff_mc_trig_mlow = 1;

      float eff_data_trig_ehigh = 1;
      float eff_mc_trig_ehigh = 1;
      float eff_data_trig_elow = 1;
      float eff_mc_trig_elow = 1;
       
      //      std::cout << "before weights " << std::endl;

      if ((!isData || isEmbedded) && ApplyLepSF) {
      	TString suffix = "mc";
      	TString suffixRatio = "ratio";
      	if (isEmbedded) {suffix = "embed"; suffixRatio = "embed_ratio";}

	correctionWS->var("e_pt")->setVal(otree->pt_1);
	correctionWS->var("e_eta")->setVal(otree->eta_1);
	correctionWS->var("e_iso")->setVal(otree->iso_1);
	correctionWS->var("m_pt")->setVal(otree->pt_2);
	correctionWS->var("m_eta")->setVal(otree->eta_2);
	correctionWS->var("m_iso")->setVal(otree->iso_2);

	/*
	float eff_data_trig_mhigh_kit = correctionWS->function("m_trg_23_binned_ic_data")->getVal();
	float eff_data_trig_mlow_kit = correctionWS->function("m_trg_8_binned_ic_data")->getVal();
	float eff_mc_trig_mhigh_kit = correctionWS->function("m_trg_23_binned_ic_"+suffix)->getVal();
	float eff_mc_trig_mlow_kit = correctionWS->function("m_trg_8_binned_ic_"+suffix)->getVal();

	float eff_data_trig_ehigh_kit = correctionWS->function("e_trg_23_binned_ic_data")->getVal();
	float eff_data_trig_elow_kit = correctionWS->function("e_trg_12_binned_ic_data")->getVal();
	float eff_mc_trig_ehigh_kit = correctionWS->function("e_trg_23_binned_ic_"+suffix)->getVal();
	float eff_mc_trig_elow_kit = correctionWS->function("e_trg_12_binned_ic_"+suffix)->getVal();
	*/
	// muon weights
	w->var("m_pt")->setVal(otree->pt_2);
	w->var("m_eta")->setVal(otree->eta_2);
	w->var("m_iso")->setVal(otree->iso_2);

	eff_data_trig_m = w->function("m_trg_ic_data")->getVal();
	eff_mc_trig_m = w->function("m_trg_ic_" + suffix)->getVal();
	sf_trig_m = w->function("m_trg_ic_" + suffixRatio)->getVal();

	eff_data_trig_mhigh = w->function("m_trg_23_binned_ic_data")->getVal();
	eff_data_trig_mlow = w->function("m_trg_8_binned_ic_data")->getVal();
	eff_mc_trig_mhigh = w->function("m_trg_23_binned_ic_"+suffix)->getVal();
	eff_mc_trig_mlow = w->function("m_trg_8_binned_ic_"+suffix)->getVal();

;
	otree->idisoweight_2 = w->function("m_idiso_ic_" + suffixRatio)->getVal();
	otree->idisoweight_antiiso_2 = w->function("m_idiso_ic_" + suffixRatio)->getVal();
	otree->trkeffweight_2 = w->function("m_trk_ratio")->getVal(); //  may be wrong

	// electron weights
	w->var("e_pt")->setVal(otree->pt_1);
	w->var("e_eta")->setVal(otree->eta_1);
	w->var("e_iso")->setVal(otree->iso_1);
	eff_data_trig_e = w->function("e_trg_ic_data")->getVal();
	eff_mc_trig_e = w->function("e_trg_ic_" + suffix)->getVal();
	sf_trig_e = w->function("e_trg_ic_" + suffixRatio)->getVal();

	eff_data_trig_ehigh = w->function("e_trg_23_binned_ic_data")->getVal();
	eff_data_trig_elow = w->function("e_trg_12_binned_ic_data")->getVal();
	eff_mc_trig_ehigh = w->function("e_trg_23_binned_ic_"+suffix)->getVal();
	eff_mc_trig_elow = w->function("e_trg_12_binned_ic_"+suffix)->getVal();
	/*
	std::cout << "eff_data_trig_ehigh : " << eff_data_trig_ehigh << "  " << eff_data_trig_ehigh_kit << std::endl;
	std::cout << "eff_mc_trig_ehigh   :  " << eff_mc_trig_ehigh << "  " << eff_mc_trig_ehigh_kit << std::endl;
	std::cout << "eff_data_trig_elow : " << eff_data_trig_elow << "  " << eff_data_trig_elow_kit << std::endl;
	std::cout << "eff_mc_trig_elow   :  " << eff_mc_trig_elow << "  " << eff_mc_trig_elow_kit << std::endl;

	std::cout << "eff_data_trig_mhigh : " << eff_data_trig_mhigh << "  " << eff_data_trig_mhigh_kit << std::endl;
	std::cout << "eff_mc_trig_mhigh   :  " << eff_mc_trig_mhigh << "  " << eff_mc_trig_mhigh_kit << std::endl;
	std::cout << "eff_data_trig_mlow : " << eff_data_trig_mlow << "  " << eff_data_trig_mlow_kit << std::endl;
	std::cout << "eff_mc_trig_mlow   :  " << eff_mc_trig_mlow << "  " << eff_mc_trig_mlow_kit << std::endl;
	*/
	if (triggerEmbed2017) {
	  if (otree->pt_1<ptTriggerEmbed2017&&fabs(otree->eta_1)>etaTriggerEmbed2017) {
	    eff_mc_trig_e = 1.0;
	    sf_trig_e = eff_data_trig_e;
	  }
	}

	otree->idisoweight_1 = w->function("e_idiso_ic_" + suffixRatio)->getVal();
	otree->idisoweight_antiiso_1 = w->function("e_idiso_ic_" + suffixRatio)->getVal();
	otree->trkeffweight_1 = w->function("e_trk_" + suffixRatio)->getVal();

	float isoweight_1_kit = 1.0;
	float isoweight_2_kit = 1.0;
	float trkeffweight_1_kit = 1.0;
	float trkeffweight_2_kit = 1.0;

	// scale factors (from KIT)
	if (era==2016){
	  if (isEmbedded) {
	    isoweight_1_kit = correctionWS->function("e_idiso_ratio_emb")->getVal();
	    isoweight_2_kit = correctionWS->function("m_idlooseiso_binned_ic_embed_ratio")->getVal();
	  }
	  else {
	    isoweight_1_kit = correctionWS->function("e_idiso_ratio")->getVal();
	    isoweight_2_kit = correctionWS->function("m_idlooseiso_binned_ic_ratio")->getVal();
	  }
	}
	else{
	  if (isEmbedded) {
	    isoweight_1_kit = correctionWS->function("e_id90_embed_kit_ratio")->getVal() * correctionWS->function("e_iso_binned_embed_kit_ratio")->getVal();
	    isoweight_2_kit = correctionWS->function("m_looseiso_binned_ic_embed_ratio")->getVal()*correctionWS->function("m_id_embed_kit_ratio")->getVal();
	  }
	  else {
	    isoweight_1_kit = correctionWS->function("e_id90_kit_ratio")->getVal() * correctionWS->function("e_iso_binned_kit_ratio")->getVal();
	    isoweight_2_kit = correctionWS->function("m_looseiso_binned_ic_ratio")->getVal()*correctionWS->function("m_id_kit_ratio")->getVal();
	  }
	}
	if (!isEmbedded){
	  correctionWS->var("e_pt")->setVal(otree->pt_1);
	  correctionWS->var("e_eta")->setVal(otree->eta_1);
	  if (era == 2018) trkeffweight_1_kit = correctionWS->function("e_trk_ratio")->getVal();
	  correctionWS->var("m_eta")->setVal(otree->eta_2);
	  correctionWS->var("m_pt")->setVal(otree->pt_2);
	  if (era==2016 || era==2018) 
	    trkeffweight_2_kit = correctionWS->function("m_trk_ratio")->getVal();
	}
	if (era == 2017) trkeffweight_1_kit = correctionWS->function("e_trk_ratio")->getVal();
	//	isoweight_1_kit *= trkeffweight_1_kit;
	//	isoweight_2_kit *= trkeffweight_2_kit;
	// KIT SF
	otree->idisoweight_1 = isoweight_1_kit;
	otree->idisoweight_2 = isoweight_2_kit;
	otree->trkeffweight_1 = trkeffweight_1_kit;
	otree->trkeffweight_2 = trkeffweight_2_kit;

	otree->trigweight_1 = sf_trig_e;
	otree->trigweight_2 = sf_trig_m;
	
	float eff_single_data = 1.0 - (1.0-eff_data_trig_e)*(1.0-eff_data_trig_m);
	float eff_single_mc   = 1.0 - (1.0-eff_mc_trig_e)  *(1.0-eff_mc_trig_m);
	if (eff_single_mc<1e-3||eff_single_data<1e-3) 
	  otree->trigweightSingle = 0.0;
	else
	  otree->trigweightSingle = eff_single_data/eff_single_mc;

	if (otree->pt_1<ptElectronSingleCut)
	  otree->trigweightSingle = sf_trig_m;
	if (otree->pt_2<ptMuonSingleCut)
	  otree->trigweightSingle = sf_trig_e;
	if (otree->pt_1<ptElectronSingleCut&&otree->pt_2<ptMuonSingleCut) {
	  otree->trigweightSingle = 0;
	}
	float eff_emu_data = 
	  eff_data_trig_mhigh*eff_data_trig_elow + 
	  eff_data_trig_mlow*eff_data_trig_ehigh -
	  eff_data_trig_mhigh*eff_data_trig_ehigh;
	float eff_emu_mc = 
	  eff_mc_trig_mhigh*eff_mc_trig_elow + 
	  eff_mc_trig_mlow*eff_mc_trig_ehigh -
	  eff_mc_trig_mhigh*eff_mc_trig_ehigh;

	if (eff_emu_mc<1e-3||eff_emu_data<1e-3)
	  otree->trigweight = 0.0;
	else
	  otree->trigweight = eff_emu_data/eff_emu_mc;

	float sf_trig_mhigh = 0;
	float sf_trig_mlow = 0;
	float sf_trig_ehigh = 0;
	float sf_trig_elow = 0;

	if (eff_data_trig_mhigh>1e-3&&eff_mc_trig_mhigh>1e-3)
	  sf_trig_mhigh = eff_data_trig_mhigh/eff_mc_trig_mhigh;
	if (eff_data_trig_mlow>1e-3&&eff_mc_trig_mlow>1e-3)
	  sf_trig_mlow = eff_data_trig_mlow/eff_mc_trig_mlow;

	if (eff_data_trig_ehigh>1e-3&&eff_mc_trig_ehigh>1e-3)
	  sf_trig_ehigh = eff_data_trig_ehigh/eff_mc_trig_ehigh;
	if (eff_data_trig_elow>1e-3&&eff_mc_trig_elow>1e-3)
	  sf_trig_elow = eff_data_trig_elow/eff_mc_trig_elow;

	if (otree->pt_1<ptElectronHighCut)
	  otree->trigweight = sf_trig_mhigh*sf_trig_elow;
	if (otree->pt_2<ptMuonHighCut)
	  otree->trigweight = sf_trig_ehigh*sf_trig_mlow;
	if (otree->pt_1<ptElectronHighCut&&otree->pt_2<ptMuonHighCut) {
	  otree->trigweight = 0;
	}

	/*
	std::cout << "pt_1 : " << otree->pt_1 << "  eta_1 = " << otree->eta_1 << "  iso_1 = " << otree->iso_1 << std::endl; 
	std::cout << "idiso_1 = " << otree->idisoweight_1 
		  << "  idiso_kit_1 = " << isoweight_1_kit << std::endl;
	std::cout << "trkeff_1 = " << otree->trkeffweight_1 
		  << "  trkeff_kit_1 = " << trkeffweight_1_kit << std::endl;

	std::cout << "pt_2 : " << otree->pt_2 << "  eta_2 = " << otree->eta_2 << "  iso_2 = " << otree->iso_2 << std::endl; 
	std::cout << "idiso_2 = " << otree->idisoweight_2 
		  << "  idiso_kit_2 = " << isoweight_2_kit << std::endl;
	std::cout << "trkeff_2 = " << otree->trkeffweight_2 
		  << "  trkeff_kit_2 = " << trkeffweight_2_kit << std::endl;	
	
	std::cout << "eff_mc_trig_e   = " << eff_mc_trig_e << std::endl;
	std::cout << "eff_mc_trig_m   = " << eff_mc_trig_m << std::endl;
	std::cout << "eff_data_trig_e = " << eff_data_trig_e << std::endl;
	std::cout << "eff_data_trig_m = " << eff_data_trig_m << std::endl;
	std::cout << "eff_single_data = " << eff_single_data << std::endl;
	std::cout << "eff_single_mc   = " << eff_single_mc << std::endl;
	std::cout << "trigweight = " << otree->trigweightSingle << std::endl;
	std::cout << "trigweight(e+mu) = " << otree->trigweightExcl << std::endl;

	std::cout << std::endl;
	*/
	float eff_emu_weight = otree->idisoweight_1 * otree->trkeffweight_1 * otree->idisoweight_2 * otree->trkeffweight_2;
	otree->effweight = eff_emu_weight * otree->trigweight;
	otree->effweightSingle = eff_emu_weight * otree->trigweightSingle;
	otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
	otree->weightSingle = otree->effweightSingle * otree->puweight * otree->mcweight; 
	otree->weight *= otree->embweight;
	otree->weightSingle *= otree->embweight;
      }
      
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

      // ************************
      // QCD background weights *
      // ************************
      if(otree->njets==0){
	otree->qcdweight_deltaR =OS_SS_njet0->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_up =  OS_SS_njet0_Par0_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_down =  OS_SS_njet0_Par0_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_up =  OS_SS_njet0_Par1_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_down =  OS_SS_njet0_Par1_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_up =  OS_SS_njet0_Par2_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_down =  OS_SS_njet0_Par2_DOWN->Eval(otree->dr_tt);
      }
      else if(otree->njets ==1) {
	otree->qcdweight_deltaR = OS_SS_njet1->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_up =  OS_SS_njet1_Par0_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_down =  OS_SS_njet1_Par0_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_up =  OS_SS_njet1_Par1_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_down =  OS_SS_njet1_Par1_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_up =  OS_SS_njet1_Par2_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_down =  OS_SS_njet1_Par2_DOWN->Eval(otree->dr_tt);
      }
      else {
	otree->qcdweight_deltaR = OS_SS_njetgt1->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_up =  OS_SS_njetgt1_Par0_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par0_down =  OS_SS_njetgt1_Par0_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_up =  OS_SS_njetgt1_Par1_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par1_down =  OS_SS_njetgt1_Par1_DOWN->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_up =  OS_SS_njetgt1_Par2_UP->Eval(otree->dr_tt);
	otree->qcdweight_deltaR_Par2_down =  OS_SS_njetgt1_Par2_DOWN->Eval(otree->dr_tt);
        
      }
      otree->qcdweight_nonclosure = hNonClosureCorrection->GetBinContent(hNonClosureCorrection->GetXaxis()->FindBin(otree->pt_2),hNonClosureCorrection->GetYaxis()->FindBin(otree->pt_1));
      otree->qcdweight_isolationcorrection = hIsolationCorrection->GetBinContent(hIsolationCorrection->GetXaxis()->FindBin(otree->pt_2),hIsolationCorrection->GetYaxis()->FindBin(otree->pt_1));
        
      otree->qcdweight=otree->qcdweight_deltaR*otree->qcdweight_nonclosure*otree->qcdweight_isolationcorrection;

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

      ////////////////////////////////////////////////////////////
      // MET and Recoil Corrections !  
      ////////////////////////////////////////////////////////////
      // !!!!!!!!!!! include electron pT correction !!!!!!!!!!!!!
    
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
	
        genTools::KITRecoilCorrections( PFMetRecoilCorrector, ApplyRecoilCorrections, // pass the value != 0 to apply corrections
          otree->met, otree->metphi,
          genV.Px(), genV.Py(),
          genL.Px(), genL.Py(),
          otree->njetshad,
          otree->met_rcmr, otree->metphi_rcmr
        );
	otree->met = otree->met_rcmr;
	otree->metphi = otree->metphi_rcmr;
      }
      
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
      // Filling variables (with corrected MET and electron pt)
      ////////////////////////////////////////////////////////////

      TLorentzVector muonLV; muonLV.SetPtEtaPhiM(otree->pt_2,
						 otree->eta_2,
						 otree->phi_2,
						 muonMass);

      TLorentzVector electronLV; electronLV.SetPtEtaPhiM(otree->pt_1,
							 otree->eta_1,
							 otree->phi_1,
							 electronMass);

      TLorentzVector metLV; metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi),
					  otree->met*TMath::Sin(otree->metphi),
					  otree->met,
					  otree->met);

      TLorentzVector puppimetLV; puppimetLV.SetXYZT(otree->puppimet*TMath::Cos(otree->puppimetphi),
						    otree->puppimet*TMath::Sin(otree->puppimetphi),
						    otree->puppimet,
						    otree->puppimet);

      TLorentzVector dileptonLV = muonLV + electronLV;
      otree->m_vis = dileptonLV.M();
      otree->pt_tt = (dileptonLV+metLV).Pt();   
      otree->pt_tt_puppi = (dileptonLV+puppimetLV).Pt();
    
      // mt TOT
      float mtTOT = 2*(otree->pt_1)*metLV.Pt()*(1-cos(DeltaPhi(electronLV,metLV)));
      mtTOT += 2*(otree->pt_2)*metLV.Pt()*(1-cos(DeltaPhi(muonLV,metLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(electronLV,muonLV))); 
      otree->mt_tot = TMath::Sqrt(mtTOT);

      mtTOT = 2*(otree->pt_1)*puppimetLV.Pt()*(1-cos(DeltaPhi(electronLV,puppimetLV)));
      mtTOT += 2*(otree->pt_2)*puppimetLV.Pt()*(1-cos(DeltaPhi(muonLV,puppimetLV))); 
      mtTOT += 2*(otree->pt_1)*(otree->pt_2)*(1-cos(DeltaPhi(electronLV,muonLV))); 
      otree->mt_tot_puppi = TMath::Sqrt(mtTOT);
    
      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;
    
      // dilepton veto
      //      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
      //      if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree, era, isEmbedded);
    
      //extra lepton veto
      TString chE("et");
      TString chMu("mt");
      otree->extraelec_veto = extra_electron_veto(electronIndex, chE, &cfg, &analysisTree, era, isEmbedded);
      otree->extramuon_veto = extra_muon_veto(muonIndex, chMu, &cfg, &analysisTree, isData);
    
      otree->mt_1 = mT(electronLV, metLV);
      otree->mt_2 = mT(muonLV, metLV);
      otree->puppimt_1 = mT(electronLV, puppimetLV);
      otree->puppimt_2 = mT(muonLV, puppimetLV);
    
      // bisector of lepton and tau transverse momenta
    
      float electronUnitX = electronLV.Px() / electronLV.Pt();
      float electronUnitY = electronLV.Py() / electronLV.Pt();
    
      float muonUnitX = muonLV.Px() / muonLV.Pt();
      float muonUnitY = muonLV.Py() / muonLV.Pt();
    
      float zetaX = electronUnitX + muonUnitX;
      float zetaY = electronUnitY + muonUnitY;
    
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
    
      zetaX = zetaX / normZeta;
      zetaY = zetaY / normZeta;
    
      float vectorVisX = electronLV.Px() + muonLV.Px();
      float vectorVisY = electronLV.Py() + muonLV.Py();
    
      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;
      otree->pzetamiss = otree->met*TMath::Cos(otree->metphi)*zetaX + otree->met*TMath::Sin(otree->metphi)*zetaY;
      otree->puppipzetamiss = otree->puppimet*TMath::Cos(otree->puppimetphi)*zetaX + otree->puppimet*TMath::Sin(otree->puppimetphi)*zetaY;
      otree->pzeta = otree->pzetamiss - 0.85*otree->pzetavis;
      otree->puppipzeta = otree->puppipzetamiss - 0.85*otree->pzetavis;

    
      bool isSRevent = true; //boolean used to compute SVFit variables only on SR events, it is set to true when running Synchronization to run SVFit on all events
      //      if(!Synch){
      //	isSRevent = (otree->dilepton_veto<0.5 &&  otree->extramuon_veto<0.5 && otree->extraelec_veto<0.5 && otree->pt_1>19 && otree->pt_2>19 && otree->byVVVLooseDeepTau2017v2p1VSjet_2>0.5);
      //	if(usePuppiMET) isSRevent = isSRevent && otree->puppimt_1<60;
      //	else isSRevent = isSRevent && otree->mt_1<60;
      //	if(ch == "mt") isSRevent = isSRevent && (otree->trg_singlemuon>0.5 || otree->trg_mutaucross>0.5);
      //	if(ch == "et") isSRevent = isSRevent && (otree->trg_singleelectron>0.5 || otree->trg_etaucross>0.5);
      //      }

      /*
	if (!isSRevent) { 
	cout << "                                        " << endl;
	cout << "========================================" << endl;
	cout << "        Event is not selected           " << endl;
	cout << "========================================" << endl;
	cout << "                                        " << endl;
	}
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
	  bool muon_match = muonIndex == analysisTree.refitvertexwithbs_muIndex[i][0] || muonIndex == analysisTree.refitvertexwithbs_muIndex[i][1];
	  bool electron_match = electronIndex == analysisTree.refitvertexwithbs_eleIndex[i][0] || electronIndex == analysisTree.refitvertexwithbs_eleIndex[i][1];
	  bool lep_match = muon_match && electron_match;

          if( lep_match )
            {
              otree->v_tracks = analysisTree.refitvertexwithbs_ntracks[i];
	      vtx_x = analysisTree.refitvertexwithbs_x[i];
	      vtx_y = analysisTree.refitvertexwithbs_y[i];
	      vtx_z = analysisTree.refitvertexwithbs_z[i];
	      for (int j=0; j<6; ++j) 
		PV_covariance[j] = analysisTree.refitvertexwithbs_cov[i][j];
            }
        }

      otree->pvx = vtx_x;
      otree->pvy = vtx_y;
      otree->pvz = vtx_z;
      TVector3 vertex_bs(vtx_bs_x,vtx_bs_y,vtx_bs_z);
      TVector3 vertex(vtx_x,vtx_y,vtx_z);
      
      TLorentzVector muonIP = ipVec_Lepton_emu(&analysisTree,muonIndex,false,vertex,era,isEmbedded);
      TLorentzVector muonIP_pvbs = ipVec_Lepton_emu(&analysisTree,muonIndex,false,vertex_bs,era,isEmbedded);

      TLorentzVector electronIP = ipVec_Lepton_emu(&analysisTree,electronIndex,true,vertex,era,isEmbedded);
      TLorentzVector electronIP_pvbs = ipVec_Lepton_emu(&analysisTree,electronIndex,true,vertex_bs,era,isEmbedded);
      
      otree->ipx_uncorr_1 = electronIP.X();
      otree->ipy_uncorr_1 = electronIP.Y();
      otree->ipz_uncorr_1 = electronIP.Z();

      otree->ipx_bs_uncorr_1 = electronIP_pvbs.X();
      otree->ipy_bs_uncorr_1 = electronIP_pvbs.Y();
      otree->ipz_bs_uncorr_1 = electronIP_pvbs.Z();

      otree->ipx_uncorr_2 = muonIP.X();
      otree->ipy_uncorr_2 = muonIP.Y();
      otree->ipz_uncorr_2 = muonIP.Z();

      otree->ipx_bs_uncorr_2 = muonIP_pvbs.X();
      otree->ipy_bs_uncorr_2 = muonIP_pvbs.Y();
      otree->ipz_bs_uncorr_2 = muonIP_pvbs.Z();

      otree->gen_ipx_1 = -9999;
      otree->gen_ipy_1 = -9999;
      otree->gen_ipz_1 = -9999;

      otree->gen_ipx_2 = -9999;
      otree->gen_ipy_2 = -9999;
      otree->gen_ipz_2 = -9999;

      otree->ipx_1 = otree->ipx_uncorr_1;
      otree->ipy_1 = otree->ipy_uncorr_1;
      otree->ipz_1 = otree->ipz_uncorr_1;
      
      otree->ipx_bs_1 = otree->ipx_bs_uncorr_1;
      otree->ipy_bs_1 = otree->ipy_bs_uncorr_1;
      otree->ipz_bs_1 = otree->ipz_bs_uncorr_1;
      
      otree->ipx_2 = otree->ipx_uncorr_2;
      otree->ipy_2 = otree->ipy_uncorr_2;
      otree->ipz_2 = otree->ipz_uncorr_2;
      
      otree->ipx_bs_2 = otree->ipx_bs_uncorr_2;
      otree->ipy_bs_2 = otree->ipy_bs_uncorr_2;
      otree->ipz_bs_2 = otree->ipz_bs_uncorr_2;

      if (!isData || isEmbedded) {

	TVector3 genVertex(otree->GenVertexX,
			   otree->GenVertexY,
			   otree->GenVertexZ);

	TVector3 ipGenMuon = genParticleIP(&analysisTree,genVertex,otree->gen_match_2,otree->eta_2,otree->phi_2);
	TVector3 ipGenElec = genParticleIP(&analysisTree,genVertex,otree->gen_match_1,otree->eta_1,otree->phi_1);
      
	otree->gen_ipx_1 = ipGenElec.X();
	otree->gen_ipy_1 = ipGenElec.Y();
	otree->gen_ipz_1 = ipGenElec.Z();

	otree->gen_ipx_2 = ipGenMuon.X();
	otree->gen_ipy_2 = ipGenMuon.Y();
	otree->gen_ipz_2 = ipGenMuon.Z();
	
	if (ApplyIpCorrection) {

	  TVector3 ipMuon_uncorr = muonIP.Vect();
	  TVector3 ipMuon_uncorr_pvbs = muonIP_pvbs.Vect(); 
	  TVector3 ipElec_uncorr = electronIP.Vect();
	  TVector3 ipElec_uncorr_pvbs = electronIP_pvbs.Vect(); 
	
	  TVector3 ipMuon      = CorrectorIpMuon->correctIp(ipMuon_uncorr,ipGenMuon,otree->eta_2);
	  TVector3 ipMuon_pvbs = CorrectorIpMuonBS->correctIp(ipMuon_uncorr_pvbs,ipGenMuon,otree->eta_2);

	  TVector3 ipElec      = CorrectorIpElec->correctIp(ipElec_uncorr,ipGenElec,otree->eta_1);
	  TVector3 ipElec_pvbs = CorrectorIpElecBS->correctIp(ipElec_uncorr_pvbs,ipGenElec,otree->eta_1);

	  otree->ipx_1 = ipElec.X();
	  otree->ipy_1 = ipElec.Y();
	  otree->ipz_1 = ipElec.Z();

	  otree->ipx_bs_1 = ipElec_pvbs.X();
	  otree->ipy_bs_1 = ipElec_pvbs.Y();
	  otree->ipz_bs_1 = ipElec_pvbs.Z();
	  
	  otree->ipx_2 = ipMuon.X();
	  otree->ipy_2 = ipMuon.Y();
	  otree->ipz_2 = ipMuon.Z();
	  
	  otree->ipx_bs_2 = ipMuon_pvbs.X();
	  otree->ipy_bs_2 = ipMuon_pvbs.Y();
	  otree->ipz_bs_2 = ipMuon_pvbs.Z();
	}
      }
      
      ImpactParameter IP;
      
      // *****
      // RefitV + BS
      // *****

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1;
      TVector3 IP1;
      IP_significance_helix_lep(&analysisTree,electronIndex, "et", vertex,PV_covariance,ipCov1,IP1);

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2;
      TVector3 IP2;
      IP_significance_helix_lep(&analysisTree, muonIndex, "mt", vertex,PV_covariance,ipCov2,IP2);

      // ******
      // PV+BS
      // ******

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1_bs;
      TVector3 IP1_bs;
      IP_significance_helix_lep(&analysisTree,electronIndex,"et",vertex_bs,PVBS_covariance,ipCov1_bs,IP1_bs);

      ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2_bs;
      TVector3 IP2_bs;
      IP_significance_helix_lep(&analysisTree,muonIndex,"mt",vertex_bs,PVBS_covariance,ipCov2_bs,IP2_bs);

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

      TVector3 Ip1(otree->ipx_uncorr_1,otree->ipy_uncorr_1,otree->ipz_uncorr_1);
      otree->IP_signif_RefitV_with_BS_uncorr_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1);
      otree->IP_signif_RefitV_with_BS_1 = otree->IP_signif_RefitV_with_BS_uncorr_1;

      Ip1.SetXYZ(otree->ipx_bs_uncorr_1,otree->ipy_bs_uncorr_1,otree->ipz_bs_uncorr_1);
      otree->IP_signif_PV_with_BS_uncorr_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_bs);
      otree->IP_signif_PV_with_BS_1 = otree->IP_signif_PV_with_BS_uncorr_1;

      otree->ip_covxx_1 = ipCov1(0,0);
      otree->ip_covxy_1 = ipCov1(0,1);
      otree->ip_covxz_1 = ipCov1(0,2);
      otree->ip_covyy_1 = ipCov1(1,1);
      otree->ip_covyz_1 = ipCov1(1,2);
      otree->ip_covzz_1 = ipCov1(2,2);

      TVector3 Ip2(otree->ipx_uncorr_2,otree->ipy_uncorr_2,otree->ipz_uncorr_2);
      otree->IP_signif_RefitV_with_BS_uncorr_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2);
      otree->IP_signif_RefitV_with_BS_2 = otree->IP_signif_RefitV_with_BS_uncorr_2;

      Ip2.SetXYZ(otree->ipx_bs_uncorr_2,otree->ipy_bs_uncorr_2,otree->ipz_bs_uncorr_2);
      otree->IP_signif_PV_with_BS_uncorr_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_bs);
      otree->IP_signif_PV_with_BS_2 = otree->IP_signif_PV_with_BS_uncorr_2;

      otree->ip_covxx_2 = ipCov2(0,0);
      otree->ip_covxy_2 = ipCov2(0,1);
      otree->ip_covxz_2 = ipCov2(0,2);
      otree->ip_covyy_2 = ipCov2(1,1);
      otree->ip_covyz_2 = ipCov2(1,2);
      otree->ip_covzz_2 = ipCov2(2,2);

      // ***********************
      // Corrected values 
      // ***********************

      if ((!isData || isEmbedded) && ApplyIpCorrection) {
	ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov1_corr;
	ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2_corr;

	ipCov1_corr =ipCov1;
	Ip1.SetXYZ(otree->ipx_1,otree->ipy_1,otree->ipz_1);
	ipCov1_corr = CorrectorIpElec->correctIpCov(ipCov1,otree->eta_1);      
	otree->IP_signif_RefitV_with_BS_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_corr);

	Ip1.SetXYZ(otree->ipx_bs_1,otree->ipy_bs_1,otree->ipz_bs_1);
	ipCov1_corr = ipCov1_bs;
	ipCov1_corr = CorrectorIpElecBS->correctIpCov(ipCov1_bs,otree->eta_1);      
	otree->IP_signif_PV_with_BS_1 = IP.CalculateIPSignificanceHelical(Ip1, ipCov1_corr);

	Ip2.SetXYZ(otree->ipx_2,otree->ipy_2,otree->ipz_2);
	ipCov2_corr = ipCov2;
	ipCov2_corr = CorrectorIpMuon->correctIpCov(ipCov2,otree->eta_2);      
	otree->IP_signif_RefitV_with_BS_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_corr);      

	Ip2.SetXYZ(otree->ipx_bs_2,otree->ipy_bs_2,otree->ipz_bs_2);
	ipCov2_corr = ipCov2_bs;
	ipCov2_corr = CorrectorIpMuonBS->correctIpCov(ipCov2_bs,otree->eta_2);      
	otree->IP_signif_PV_with_BS_2 = IP.CalculateIPSignificanceHelical(Ip2, ipCov2_corr);      
      }

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

	muonScaleSys->Eval(utils::EMU);
	electronScaleSys->SetElectronIndex(electronIndex);
	electronScaleSys->SetIsEmbedded(isEmbedded);
	electronScaleSys->SetAC1B(&analysisTree);
	electronScaleSys->Eval(utils::EMU);
    
      }

      selEvents++;
      
      //      cout << "Once again puppimet central : " << otree->puppimet << endl;
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

  if (muonScaleSys != 0) {
    muonScaleSys->Write("",TObject::kOverwrite);
    delete muonScaleSys;
  }
  
  if (electronScaleSys != 0) {
    electronScaleSys->Write("",TObject::kOverwrite);
    delete electronScaleSys;
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

//// FILLING FUNCTIONS //////

void FillVertices(const AC1B *analysisTree, Synch17Tree *otree, const bool isData){

  otree->RecoVertexX = analysisTree->primvertex_x;
  otree->RecoVertexY = analysisTree->primvertex_y;
  otree->RecoVertexZ = analysisTree->primvertex_z;  

  otree->pvx_bs = analysisTree->primvertexwithbs_x;
  otree->pvy_bs = analysisTree->primvertexwithbs_y;
  otree->pvz_bs = analysisTree->primvertexwithbs_z;

  if(!isData){
    for (unsigned int igen = 0; igen < analysisTree->genparticles_count; ++igen) {

  //here fill the generator vertices to have the gen information present in tree PER GOOD RECO EVENT
  //Note: we may want to add constraint that the W and Z are prompt. If we remove these, may get in trouble with a DY or W MC sample..

      if ( analysisTree->genparticles_pdgid[igen] == 23 || 
	   analysisTree->genparticles_pdgid[igen] == 24 ||
	   analysisTree->genparticles_pdgid[igen] == -24 ||
	   analysisTree->genparticles_pdgid[igen] == 25 || 
	   analysisTree->genparticles_pdgid[igen] == 35 || 
	   analysisTree->genparticles_pdgid[igen] == 36 ||
	   analysisTree->genparticles_pdgid[igen] == 6 ||
	   analysisTree->genparticles_pdgid[igen] == -6 ) {
        otree->GenVertexX = analysisTree->genparticles_vx[igen];
        otree->GenVertexY = analysisTree->genparticles_vy[igen];
        otree->GenVertexZ = analysisTree->genparticles_vz[igen];
        break;
      }
    }
  }
  else {//if it is data, fill with something recognisable nonsensible
    otree->GenVertexX = 0;
    otree->GenVertexY = 0;
    otree->GenVertexZ = 0;
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

  gentree->VertexX=-9999;
  gentree->VertexY=-9999;
  gentree->VertexZ=-99999;


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

//fill the otree with the electron/muon variables in channel emu
void FillElMu(const AC1B *analysisTree, Synch17Tree *otree, int electronIndex, float dRIsoElectron, int muonIndex, float dRIsoMuon, int era, bool isEmbedded){
  
  float sf_eleES = 1.;
  if (isEmbedded) sf_eleES = EmbedElectronES_SF(analysisTree, era, electronIndex);  

  otree->pt_1  = 	sf_eleES*analysisTree->electron_pt[electronIndex];
  otree->pt_uncorr_1 =  analysisTree->electron_pt[electronIndex];
  otree->eta_1 = 	analysisTree->electron_eta[electronIndex];
  otree->phi_1 = 	analysisTree->electron_phi[electronIndex];
  otree->m_1 = 		electronMass;
  otree->q_1 = -1;
  if (analysisTree->electron_charge[electronIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->electron_genmatch[electronIndex];

  otree->iso_1 =  abs_Iso_et(electronIndex, analysisTree, dRIsoElectron) / (sf_eleES*analysisTree->electron_pt[electronIndex]);

  otree->d0_1 = analysisTree->electron_dxy[electronIndex];
  otree->dZ_1 = analysisTree->electron_dz[electronIndex];

  otree->pt_2  = 	analysisTree->muon_pt[muonIndex];
  otree->eta_2 = 	analysisTree->muon_eta[muonIndex];
  otree->phi_2 = 	analysisTree->muon_phi[muonIndex];
  otree->m_2   =  muonMass;
  otree->q_2   =  analysisTree->muon_charge[muonIndex];
  otree->gen_match_2 = analysisTree->muon_genmatch[muonIndex];
  otree->iso_2 = abs_Iso_mt(muonIndex, analysisTree, dRIsoMuon) / analysisTree->muon_pt[muonIndex];
  otree->d0_2 = analysisTree->muon_dxy[muonIndex];
  otree->dZ_2 = analysisTree->muon_dz[muonIndex];
  otree->q_2 = -1;
  if (analysisTree->muon_charge[muonIndex]>0)
    otree->q_2 = 1;

  otree->dr_tt = deltaR(otree->eta_1,otree->phi_1,otree->eta_2,otree->phi_2);

}
