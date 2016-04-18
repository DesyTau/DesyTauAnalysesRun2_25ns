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
#include "TRandom3.h"

#include "TLorentzVector.h"

#include "TRandom.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "HTT-utilities/QCDModelingEMu/interface/QCDModelForEMu.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

float totalTransverseMass(TLorentzVector l1, 
			  TLorentzVector l2,
			  TLorentzVector l3) {

  TLorentzVector totalLV = l1 + l2 + l3;
  float    totalET = l1.Pt() +  l2.Pt() + l3.Pt();
  float     mTtot = TMath::Sqrt(totalET*totalET-totalLV.Pt()*totalLV.Pt());
  return    mTtot;

}

void computeDzeta(float metX,  float metY,
		  float zetaX, float zetaY,
		  float pzetavis,
		  float & pzetamiss,
		  float & dzeta) {

  pzetamiss = metX*zetaX + metY*zetaY;
  dzeta = pzetamiss - 0.85*pzetavis;

}


float topPtWeight(float pt1,
		  float pt2) {

  if (pt1>400) pt1 = 400;
  if (pt2>400) pt2 = 400;

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter
  
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  return TMath::Sqrt(w1*w2);

}

SVfitStandaloneAlgorithm SVFitMassComputation(svFitStandalone::MeasuredTauLepton svFitEle,
					      svFitStandalone::MeasuredTauLepton svFitMu,
					      double measuredMVAMETx,
					      double measuredMVAMETy,
					      TMatrixD covMVAMET,
					      TFile * inputFile_visPtResolution
					      ) {

  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitEle);
  measuredTauLeptons.push_back(svFitMu);
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMVAMETx, measuredMVAMETy, covMVAMET, 0);
  algo.addLogM(false);  
  algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();
  
  return algo;

}


//struct myclass {
//  bool operator() (int i,int j) { return (i<j);}
//} myobject, myobjectX;

//const float electronMass = 0;
//const float muonMass = 0.10565837;
//const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool applyInclusiveSelection = cfg.get<bool>("ApplyInclusiveSelection");
  const bool computeSVFitMass = cfg.get<bool>("ComputeSVFitMass");

  const bool isData = cfg.get<bool>("IsData");
  const bool isDY   = cfg.get<bool>("IsDY");
  const bool isW    = cfg.get<bool>("IsW");
  const bool isTOP    = cfg.get<bool>("IsTOP");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");
  const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
  const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");

  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
  const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");

  // jets
  const string bTagDiscriminator = cfg.get<string>("BTagDiscriminator");
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");
  const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
  const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");

  TString LowPtLegElectron(lowPtLegElectron);
  TString HighPtLegElectron(highPtLegElectron);
  
  TString LowPtLegMuon(lowPtLegMuon);
  TString HighPtLegMuon(highPtLegMuon);

  TString BTagDiscriminator(bTagDiscriminator);

  const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");

  const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
  const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");

  const string Electron17TriggerFile = cfg.get<string>("Electron17TriggerEff");
  const string Electron12TriggerFile = cfg.get<string>("Electron12TriggerEff");

  const float muonScale = cfg.get<float>("MuonScale");
  const float eleScaleBarrel = cfg.get<float>("EleScaleBarrel");
  const float eleScaleEndcap = cfg.get<float>("EleScaleEndcap");

  const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
  TString RecoilMvaFileName(recoilMvaFileName);

  const string metSysFileName   = cfg.get<string>("MetSysFileName");
  TString MetSysFileName(metSysFileName);

  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
  TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
  TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
  TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * histWeightsSkimH = new TH1D("histWeightsSkimH","",1,-0.5,0.5);
  TH1D * histWeightsTTH = new TH1D("histWeightsTTH","",1,-0.5,0.5);


  TTree * tree = new TTree("TauCheck","TauCheck");
   // Declaration of leaf types
  unsigned int minRun = 99999999;
  unsigned int maxRun = 0;
  Int_t           run;
  Int_t           lumi;
  Int_t           evt;
  Int_t           npv;
  Int_t           npu;
  Float_t         rho;
  Float_t         mcweight;
  Float_t         puweight;
  Float_t         trigweight_1;
  Float_t         trigweight_2;
  Float_t         trigweight;
  Float_t         idweight_1;
  Float_t         idweight_2;
  Float_t         isoweight_1;
  Float_t         isoweight_2;
  Float_t         effweight;
  Float_t         fakeweight;
  Float_t         embeddedWeight;
  Float_t         signalWeight;
  Float_t         topptweight;

  Float_t         qcdweight;
  Float_t         qcdweightup;
  Float_t         qcdweightdown;

  Float_t         qcdweight_nodzeta;
  Float_t         qcdweightup_nodzeta;
  Float_t         qcdweightdown_nodzeta;

  Float_t         weight;

  Float_t         m_vis;
  Float_t         m_vis_muUp;
  Float_t         m_vis_muDown;
  Float_t         m_vis_eUp;
  Float_t         m_vis_eDown;
  Float_t         m_vis_scaleUp;
  Float_t         m_vis_scaleDown;
  Float_t         m_vis_resoUp;
  Float_t         m_vis_resoDown;

  Float_t         mTtot;
  Float_t         mTtot_muUp;
  Float_t         mTtot_muDown;
  Float_t         mTtot_eUp;
  Float_t         mTtot_eDown;
  Float_t         mTtot_scaleUp;
  Float_t         mTtot_scaleDown;
  Float_t         mTtot_resoUp;
  Float_t         mTtot_resoDown;

  Float_t         m_sv;
  Float_t         mt_sv;
  Float_t         pt_sv;
  Float_t         eta_sv;
  Float_t         phi_sv;
  
  Float_t         m_sv_eUp;
  Float_t         m_sv_eDown;
  Float_t         m_sv_muUp;
  Float_t         m_sv_muDown;
  Float_t         m_sv_scaleUp;
  Float_t         m_sv_scaleDown;
  Float_t         m_sv_resoUp;
  Float_t         m_sv_resoDown;

  Float_t         mt_sv_eUp;
  Float_t         mt_sv_eDown;
  Float_t         mt_sv_muUp;
  Float_t         mt_sv_muDown;
  Float_t         mt_sv_scaleUp;
  Float_t         mt_sv_scaleDown;
  Float_t         mt_sv_resoUp;
  Float_t         mt_sv_resoDown;


  Float_t         pt_1;
  Float_t         pt_Up_1;
  Float_t         pt_Down_1;

  Float_t         phi_1;
  Float_t         eta_1;
  Float_t         m_1;
  Int_t           q_1;
  Float_t         iso_1;
  Float_t         mva_1;
  Float_t         d0_1;
  Float_t         dZ_1;
  Float_t         mt_1;

  Float_t         pt_2;
  Float_t         pt_Up_2;
  Float_t         pt_Down_2;

  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         m_2;
  Int_t           q_2;
  Float_t         iso_2;
  Float_t         d0_2;
  Float_t         dZ_2;
  Float_t         mva_2;
  Float_t         mt_2;

  Bool_t          os;
  Bool_t          dilepton_veto;
  Bool_t          extraelec_veto;
  Bool_t          extramuon_veto;

  Float_t         met;
  Float_t         metphi;
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;

  Float_t         mvamet;
  Float_t         mvametphi;
  Float_t         mvacov00;
  Float_t         mvacov01;
  Float_t         mvacov10;
  Float_t         mvacov11;

  Float_t         mvamet_uncorr;
  Float_t         mvametphi_uncorr;

  Float_t         mvamet_resoUp;
  Float_t         mvametphi_resoUp;

  Float_t         mvamet_resoDown;
  Float_t         mvametphi_resoDown;

  Float_t         mvamet_scaleUp;
  Float_t         mvametphi_scaleUp;

  Float_t         mvamet_scaleDown;
  Float_t         mvametphi_scaleDown;

  Float_t         genmet;
  Float_t         genmetphi;

  Float_t         pt_tt;
  Float_t         dr_tt;
  Float_t         dphi_tt;

  Float_t         pzetavis;
  Float_t         pzetamiss;
  Float_t         dzeta;

  Float_t         pzetamiss_mvamet; 
  Float_t         dzeta_mvamet;
  
  Float_t         pzetamiss_mvamet_uncorr; 
  Float_t         dzeta_mvamet_uncorr;
  
  Float_t         pzetamiss_mvamet_scaleUp; 
  Float_t         dzeta_mvamet_scaleUp;

  Float_t         pzetamiss_mvamet_scaleDown; 
  Float_t         dzeta_mvamet_scaleDown;

  Float_t         pzetamiss_mvamet_resoUp; 
  Float_t         dzeta_mvamet_resoUp;

  Float_t         pzetamiss_mvamet_resoDown; 
  Float_t         dzeta_mvamet_resoDown;

  Float_t         pzetamiss_genmet; 
  Float_t         dzeta_genmet;

  Float_t         mva_gf;

  Int_t           njets;
  Int_t           njetspt20;
  Float_t         jpt_1;
  Float_t         jeta_1;
  Float_t         jphi_1;
  Float_t         jptraw_1;
  Float_t         jptunc_1;
  Float_t         jmva_1;
  Float_t         jlrm_1;
  Int_t           jctm_1;
  Int_t           gen_match_1;

  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jptraw_2;
  Float_t         jptunc_2;
  Float_t         jmva_2;
  Float_t         jlrm_2;
  Int_t           jctm_2;
  Int_t           gen_match_2;

  Float_t         mjj;
  Float_t         jdeta;
  Int_t           njetingap;

  Int_t           nbtag;
  Float_t         bpt;
  Float_t         beta;
  Float_t         bphi;

  Float_t         nuPx;
  Float_t         nuPy;
  Float_t         nuPz;

  Float_t         lepPx;
  Float_t         lepPy;
  Float_t         lepPz;
  
  Float_t         bosonPx;
  Float_t         bosonPy;
  Float_t         bosonPz;
  Float_t         bosonMass;
  
  UInt_t          npartons;

  Bool_t isZLL;
  Bool_t isZMM;
  Bool_t isZEE;
  Bool_t isZTT;

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("evt", &evt, "evt/I");
  tree->Branch("npv", &npv, "npv/I");
  tree->Branch("npu", &npu, "npu/I");
  tree->Branch("rho", &rho, "rho/F");

  tree->Branch("isZLL",&isZLL,"isZLL/O");
  tree->Branch("isZEE",&isZEE,"isZEE/O");
  tree->Branch("isZMM",&isZMM,"isZMM/O");
  tree->Branch("isZTT",&isZTT,"isZTT/O");

  tree->Branch("mcweight", &mcweight, "mcweight/F");
  tree->Branch("puweight", &puweight, "puweight/F");
  tree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
  tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  tree->Branch("trigweight", &trigweight, "trigweight/F");
  tree->Branch("idweight_1", &idweight_1, "idweight_1/F");
  tree->Branch("idweight_2", &idweight_2, "idweight_2/F");
  tree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
  tree->Branch("isoweight_2", &isoweight_2, "isoweight_2/F");
  tree->Branch("effweight", &effweight, "effweight/F");
  tree->Branch("fakeweight", &fakeweight, "fakeweight/F");
  tree->Branch("embeddedWeight", &embeddedWeight, "embeddedWeight/F");
  tree->Branch("signalWeight", &signalWeight, "signalWeight/F");
  tree->Branch("topptweight", &topptweight, "topptweight/F");

  tree->Branch("qcdweight", &qcdweight, "qcdweight/F");
  tree->Branch("qcdweightup", &qcdweightup, "qcdweightup/F");
  tree->Branch("qcdweightdown", &qcdweightdown, "qcdweightdown/F");

  tree->Branch("qcdweight_nodzeta", &qcdweight_nodzeta, "qcdweight_nodzeta/F");
  tree->Branch("qcdweightup_nodzeta", &qcdweightup_nodzeta, "qcdweightup_nodzeta/F");
  tree->Branch("qcdweightdown_nodzeta", &qcdweightdown_nodzeta, "qcdweightdown_nodzeta/F");

  tree->Branch("weight", &weight, "weight/F");

  tree->Branch("m_vis",        &m_vis,        "m_vis/F");
  tree->Branch("m_vis_muUp",   &m_vis_muUp,   "m_vis_muUp/F");
  tree->Branch("m_vis_muDown", &m_vis_muDown, "m_vis_muDown/F");
  tree->Branch("m_vis_eUp",    &m_vis_eUp,    "m_vis_eUp/F");
  tree->Branch("m_vis_eDown",  &m_vis_eDown,  "m_vis_eDown/F");
  tree->Branch("m_vis_scaleUp",   &m_vis_scaleUp,   "m_vis_scaleUp/F");
  tree->Branch("m_vis_scaleDown", &m_vis_scaleDown, "m_vis_scaleDown/F");
  tree->Branch("m_vis_resoUp",    &m_vis_resoUp,    "m_vis_resoUp/F");
  tree->Branch("m_vis_resoDown",  &m_vis_resoDown,  "m_vis_resoDown/F"); 

  tree->Branch("mTtot",        &mTtot,        "mTtot/F");
  tree->Branch("mTtot_muUp",   &mTtot_muUp,   "mTtot_muUp/F");
  tree->Branch("mTtot_muDown", &mTtot_muDown, "mTtot_muDown/F");
  tree->Branch("mTtot_eUp",    &mTtot_eUp,    "mTtot_eUp/F");
  tree->Branch("mTtot_eDown",  &mTtot_eDown,  "mTtot_eDown/F"); 
  tree->Branch("mTtot_scaleUp",   &mTtot_scaleUp,   "mTtot_scaleUp/F");
  tree->Branch("mTtot_scaleDown", &mTtot_scaleDown, "mTtot_scaleDown/F");
  tree->Branch("mTtot_resoUp",    &mTtot_resoUp,    "mTtot_resoUp/F");
  tree->Branch("mTtot_resoDown",  &mTtot_resoDown,  "mTtot_resoDown/F"); 

  tree->Branch("m_sv",    &m_sv,   "m_sv/F");
  tree->Branch("mt_sv",   &mt_sv,  "mt_sv/F");
  tree->Branch("pt_sv",   &pt_sv,  "pt_sv/F");
  tree->Branch("eta_sv",  &eta_sv, "eta_sv/F");
  tree->Branch("phi_sv",  &phi_sv, "phi_sv/F");

  tree->Branch("m_sv_scaleUp",   &m_sv_scaleUp,   "m_sv_scaleUp/F");
  tree->Branch("m_sv_scaleDown", &m_sv_scaleDown, "m_sv_scaleDown/F");
  tree->Branch("m_sv_resoUp",    &m_sv_resoUp,    "m_sv_resoUp/F");
  tree->Branch("m_sv_resoDown",  &m_sv_resoDown,  "m_sv_resoDown/F");
  tree->Branch("m_sv_eUp",       &m_sv_eUp,       "m_sv_eUp/F");
  tree->Branch("m_sv_eDown",     &m_sv_eDown,     "m_sv_eDown/F");
  tree->Branch("m_sv_muUp",      &m_sv_muUp,      "m_sv_muUp/F");
  tree->Branch("m_sv_muDown",    &m_sv_muDown,    "m_sv_muDown/F");

  tree->Branch("mt_sv_scaleUp",   &mt_sv_scaleUp,   "mt_sv_scaleUp/F");
  tree->Branch("mt_sv_scaleDown", &mt_sv_scaleDown, "mt_sv_scaleDown/F");
  tree->Branch("mt_sv_resoUp",    &mt_sv_resoUp,    "mt_sv_resoUp/F");
  tree->Branch("mt_sv_resoDown",  &mt_sv_resoDown,  "mt_sv_resoDown/F");
  tree->Branch("mt_sv_eUp",       &mt_sv_eUp,       "mt_sv_eUp/F");
  tree->Branch("mt_sv_eDown",     &mt_sv_eDown,     "mt_sv_eDown/F");
  tree->Branch("mt_sv_muUp",      &mt_sv_muUp,      "mt_sv_muUp/F");
  tree->Branch("mt_sv_muDown",    &mt_sv_muDown,    "mt_sv_muDown/F");

  tree->Branch("pt_1", &pt_1, "pt_1/F");
  tree->Branch("pt_Up_1", &pt_Up_1, "pt_Up_1/F");
  tree->Branch("pt_Down_1", &pt_Down_1, "pt_Down_1/F");
  tree->Branch("phi_1", &phi_1, "phi_1/F");
  tree->Branch("eta_1", &eta_1, "eta_1/F");
  tree->Branch("m_1", &m_1, "m_1/F");
  tree->Branch("q_1", &q_1, "q_1/I");
  tree->Branch("iso_1", &iso_1, "iso_1/F");
  tree->Branch("mva_1", &mva_1, "mva_1/F");
  tree->Branch("d0_1", &d0_1, "d0_1/F");
  tree->Branch("dZ_1", &dZ_1, "dZ_1/F");
  tree->Branch("mt_1", &mt_1, "mt_1/F");
  tree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");

  tree->Branch("pt_2", &pt_2, "pt_2/F");
  tree->Branch("pt_Up_2", &pt_Up_2, "pt_Up_2/F");
  tree->Branch("pt_Down_2", &pt_Down_2, "pt_Down_2/F");

  tree->Branch("phi_2", &phi_2, "phi_2/F");
  tree->Branch("eta_2", &eta_2, "eta_2/F");
  tree->Branch("m_2", &m_2, "m_2/F");
  tree->Branch("q_2", &q_2, "q_2/I");
  tree->Branch("iso_2", &iso_2, "iso_2/F");
  tree->Branch("d0_2", &d0_2, "d0_2/F");
  tree->Branch("dZ_2", &dZ_2, "dZ_2/F");
  tree->Branch("mva_2", &mva_2, "mva_2/F");
  tree->Branch("mt_2", &mt_2, "mt_2/F");
  tree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");

  tree->Branch("os", &os, "os/O");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/O");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/O");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/O");

  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("metcov00", &metcov00, "metcov00/F");
  tree->Branch("metcov01", &metcov01, "metcov01/F");
  tree->Branch("metcov10", &metcov10, "metcov10/F");
  tree->Branch("metcov11", &metcov11, "metcov11/F");

  tree->Branch("mvamet", &mvamet, "mvamet/F");
  tree->Branch("mvametphi", &mvametphi, "mvametphi/F");
  tree->Branch("mvacov00", &mvacov00, "mvacov00/F");
  tree->Branch("mvacov01", &mvacov01, "mvacov01/F");
  tree->Branch("mvacov10", &mvacov10, "mvacov10/F");
  tree->Branch("mvacov11", &mvacov11, "mvacov11/F");

  tree->Branch("mvamet_uncorr", &mvamet_uncorr, "mvamet_uncorr/F");
  tree->Branch("mvametphi_uncorr", &mvametphi_uncorr, "mvametphi_uncorr/F");
  
  tree->Branch("mvamet_resoUp", &mvamet_resoUp, "mvamet_resoUp/F");
  tree->Branch("mvametphi_resoUp", &mvametphi_resoUp, "mvametphi_resoUp/F");
  
  tree->Branch("mvamet_resoDown", &mvamet_resoDown, "mvamet_resoDown/F");
  tree->Branch("mvametphi_resoDown", &mvametphi_resoDown, "mvametphi_resoDown/F");
  
  tree->Branch("mvamet_scaleUp", &mvamet_scaleUp, "mvamet_scaleUp/F");
  tree->Branch("mvametphi_scaleUp", &mvametphi_scaleUp, "mvametphi_scaleUp/F");
  
  tree->Branch("mvamet_scaleDown", &mvamet_scaleDown, "mvamet_scaleDown/F");
  tree->Branch("mvametphi_scaleDown", &mvametphi_scaleDown, "mvametphi_scaleDown/F");
  
  tree->Branch("genmet", &genmet, "genmet/F");
  tree->Branch("genmetphi", &genmetphi, "genmetphi/F");

  tree->Branch("pt_tt", &pt_tt, "pt_tt/F");
  tree->Branch("dr_tt", &dr_tt, "dr_tt/F");
  tree->Branch("dphi_tt", &dphi_tt, "dphi_tt/F");

  tree->Branch("pzetavis", &pzetavis, "pzetavis/F");

  tree->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
  tree->Branch("dzeta",&dzeta,"dzeta/F");

  tree->Branch("pzetamiss_mvamet", &pzetamiss_mvamet, "pzetamiss_mvamet/F");
  tree->Branch("dzeta_mvamet",&dzeta_mvamet,"dzeta_mvamet/F");

  tree->Branch("pzetamiss_mvamet_uncorr", &pzetamiss_mvamet_uncorr, "pzetamiss_mvamet_uncorr/F");
  tree->Branch("dzeta_mvamet_uncorr",&dzeta_mvamet_uncorr,"dzeta_mvamet_uncorr/F");

  tree->Branch("pzetamiss_mvamet_resoUp", &pzetamiss_mvamet_resoUp, "pzetamiss_mvamet_resoUp/F");
  tree->Branch("dzeta_mvamet_resoUp",&dzeta_mvamet_resoUp,"dzeta_mvamet_resoUp/F");

  tree->Branch("pzetamiss_mvamet_resoDown", &pzetamiss_mvamet_resoDown, "pzetamiss_mvamet_resoDown/F");
  tree->Branch("dzeta_mvamet_resoDown",&dzeta_mvamet_resoDown,"dzeta_mvamet_resoDown/F");

  tree->Branch("pzetamiss_mvamet_scaleUp", &pzetamiss_mvamet_scaleUp, "pzetamiss_mvamet_scaleUp/F");
  tree->Branch("dzeta_mvamet_scaleUp",&dzeta_mvamet_scaleUp,"dzeta_mvamet_scaleUp/F");

  tree->Branch("pzetamiss_mvamet_scaleDown", &pzetamiss_mvamet_scaleDown, "pzetamiss_mvamet_scaleDown/F");
  tree->Branch("dzeta_mvamet_scaleDown",&dzeta_mvamet_scaleDown,"dzeta_mvamet_scaleDown/F");

  tree->Branch("pzetamiss_genmet", &pzetamiss_genmet, "pzetamiss_genmet/F");
  tree->Branch("dzeta_genmet",&dzeta_genmet,"dzeta_genmet/F");

  tree->Branch("mva_gf", &mva_gf, "mva_gf/F");

  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");

  tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
  tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
  tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
  tree->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
  tree->Branch("jptunc_1", &jptunc_1, "jptunc_1/F");
  tree->Branch("jmva_1", &jmva_1, "jmva_1/F");
  tree->Branch("jlrm_1", &jlrm_1, "jlrm_1/F");
  tree->Branch("jctm_1", &jctm_1, "jctm_1/I");

  tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
  tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
  tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
  tree->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
  tree->Branch("jptunc_2", &jptunc_2, "jptunc_2/F");
  tree->Branch("jmva_2", &jmva_2, "jlrm_2/F");
  tree->Branch("jlrm_2", &jlrm_2, "jlrm_2/F");
  tree->Branch("jctm_2", &jctm_2, "jctm_2/I");

  tree->Branch("mjj", &mjj, "mjj/F");
  tree->Branch("jdeta", &jdeta, "jdeta/F");
  tree->Branch("njetingap", &njetingap, "njetingap/I");

  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("bpt", &bpt, "bpt/F");
  tree->Branch("beta", &beta, "beta/F");
  tree->Branch("bphi", &bphi, "bphi/F");

  tree->Branch("nuPx",&nuPx,"nuPx/F");
  tree->Branch("nuPy",&nuPy,"nuPy/F");
  tree->Branch("nuPz",&nuPz,"nuPz/F");
  
  tree->Branch("lepPx",&lepPx,"lepPx/F");
  tree->Branch("lepPy",&lepPy,"lepPy/F");
  tree->Branch("lepPz",&lepPz,"lepPz/F");
  
  tree->Branch("bosonPx",&bosonPx,"bosonPx/F");
  tree->Branch("bosonPy",&bosonPy,"bosonPy/F");
  tree->Branch("bosonPz",&bosonPz,"bosonPz/F");
  tree->Branch("bosonMass",&bosonMass,"bosonMass/F");
  
  tree->Branch("npartons",&npartons,"npartons/i");

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  std::vector<Period> periods;
  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
  
  if (isData) {
  std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
  if (inputFileStream.fail()) {
       std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
       std::cout << "please check" << std::endl;
       std::cout << "quitting program" << std::endl;
       exit(-1);
     }
  for(std::string s; std::getline(inputFileStream, s); )
    {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
   }

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Feb02.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Fall15_PU25_V1.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Lepton Scale Factors
  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  ScaleFactor * SF_muon17 = new ScaleFactor();
  SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));
  ScaleFactor * SF_muon8 = new ScaleFactor();
  SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));
  ScaleFactor * SF_electronIdIso = new ScaleFactor();
  SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
  ScaleFactor * SF_electron17 = new ScaleFactor();
  SF_electron17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron17TriggerFile));
  ScaleFactor * SF_electron12 = new ScaleFactor();
  SF_electron12->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron12TriggerFile));

  RecoilCorrector recoilMvaMetCorrector(RecoilMvaFileName);
  MEtSys metSys(MetSysFileName);

  // SV fit mass
  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);
  TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());

  // qcd weight (dzeta cut)
  QCDModelForEMu qcdWeight("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu.root");
  // qcd weight DZeta cut
  QCDModelForEMu qcdWeightNoDzeta("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_nodzeta.root");

  // BTag scale factors
  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2.csv");
  BTagCalibrationReader reader_BC(&calib,BTagEntry::OP_MEDIUM,"mujets","central");           // systematics type
  BTagCalibrationReader reader_Light(&calib,BTagEntry::OP_MEDIUM,"incl","central");           // systematics type

  //  std::cout << "SF_light (eta=0.6,pt=20.1) : " << reader_Light.eval(BTagEntry::FLAV_UDSG, 0.5, 20.1) << std::endl;
  //  std::cout << "SF_light (eta=2.1,pt=20.1) : " << reader_Light.eval(BTagEntry::FLAV_UDSG, 2.1, 20.1) << std::endl;
  //  std::cout << "SF_bc    (eta=0.6,pt=30.1) : " << reader_BC.eval(BTagEntry::FLAV_B, 0.5, 30.1) << std::endl;
  //  std::cout << "SF_bc    (eta=2.1,pt=30.1) : " << reader_BC.eval(BTagEntry::FLAV_B, 2.1, 30.1) << std::endl;

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  float MaxBJetPt = 670.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 30.;

  //  exit(-1);

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
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
	if (isData)
	  histWeightsH->Fill(0.,1.);
	else
	  histWeightsH->Fill(0.,genweight);
      }
    }

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
      if (analysisTree.event_run>maxRun)
        maxRun = analysisTree.event_run;

      if (analysisTree.event_run<minRun)
        minRun = analysisTree.event_run;
     
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      isZLL = false;
      isZEE = false;
      isZMM = false;
      isZTT = false;

      //      bool isPrompMuPlus = false;
      //      bool isPrompMuMinus = false;
      //      bool isPrompElePlus = false;
      //      bool isPrompEleMinus = false;

      nuPx = 0;
      nuPy = 0;
      nuPz = 0;
      lepPx = 0;
      lepPy = 0;
      lepPz = 0;
      bosonPx = 0;
      bosonPy = 0;
      bosonPz = 0;
      bosonMass = -1;

      float topPt = -1;
      float antitopPt = -1;
      
      bool isZfound = false;
      bool isWfound = false;
      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0,0,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0,0,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);

      if (!isData) {
	

	for (unsigned int igentau=0; igentau < analysisTree.gentau_count; ++igentau) {
	  TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
					      analysisTree.gentau_py[igentau],
					      analysisTree.gentau_pz[igentau],
					      analysisTree.gentau_e[igentau]);
	  TLorentzVector tauVisLV; tauVisLV.SetXYZT(analysisTree.gentau_visible_px[igentau],
						    analysisTree.gentau_visible_py[igentau],
						    analysisTree.gentau_visible_pz[igentau],
						    analysisTree.gentau_visible_e[igentau]);
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
	    promptTausFirstCopy.push_back(tauLV);
	    promptTausLV += tauLV;
	    wDecayProductsLV += tauLV;
	  }
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {	
	    promptTausLastCopy.push_back(tauVisLV);
	    promptVisTausLV += tauVisLV;
	  }
	  
	}

	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

	  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);

	  if (analysisTree.genparticles_pdgid[igen]==6)
	    topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);

	  if (analysisTree.genparticles_pdgid[igen]==-6)
	    antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				    analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);


	  if (analysisTree.genparticles_pdgid[igen]==23) { 
	    isZfound = true;
	    zBosonLV = genLV;
	  }
	  if (abs(analysisTree.genparticles_pdgid[igen])==24) { 
	    isWfound = true;
	    wBosonLV = genLV;
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptElectrons.push_back(genLV);
	      promptElectronsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptMuons.push_back(genLV);
	      promptMuonsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==12||
	      fabs(analysisTree.genparticles_pdgid[igen])==14||
	      fabs(analysisTree.genparticles_pdgid[igen])==16)  {
	    if (analysisTree.genparticles_fromHardProcess[igen]&&
		!analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      promptNeutrinos.push_back(genLV);
	      promptNeutrinosLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  

	}

	if (isDY) {
	  if (promptTausFirstCopy.size()==2) {
	    isZTT = true; isZMM = false; isZEE = false;
	    bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz(); 
	    bosonMass = promptTausLV.M();
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	  }
	  else if (promptMuons.size()==2) {
	    isZTT = false; isZMM = true; isZEE = false;
	    bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz(); 
	    bosonMass = promptMuonsLV.M(); 
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	  }
	  else {
	    isZTT = false; isZMM = false; isZEE = true;
	    bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz(); 
	    bosonMass = promptElectronsLV.M();
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	  }
	}
	else if (isW) {
	  bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
	  bosonMass = wDecayProductsLV.M();
	  if (promptTausLastCopy.size()==1) { 
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	  }
	  else if (promptMuons.size()==1) { 
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	  }
	  else { 
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	  }
	}
	else {
	  TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
	  bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
	  TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
	  lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
	}

	/*
	std::cout << "Taus (first copy) : " << promptTausFirstCopy.size() << std::endl;
	std::cout << "Taus (last copy)  : " << promptTausLastCopy.size() << std::endl;
	std::cout << "Muons             : " << promptMuons.size() << std::endl;
	std::cout << "Electrons         : " << promptElectrons.size() << std::endl;
	std::cout << "Neutrinos         : " << promptNeutrinos.size() << std::endl;
	if (isZfound)
	  std::cout << "ZBoson  px = " << zBosonLV.Px() 
		    << "        py = " << zBosonLV.Py()
		    << "        pz = " << zBosonLV.Pz() << std::endl;
	if (isWfound) {
	  std::cout << "WBoson  px = " << wBosonLV.Px() 
		    << "        py = " << wBosonLV.Py()
		    << "        pz = " << wBosonLV.Pz() << std::endl;
	  std::cout << "W(Dec)  px = " << wDecayProductsLV.Px() 
		    << "        py = " << wDecayProductsLV.Py()
		    << "        pz = " << wDecayProductsLV.Pz() << std::endl;
	}      
	std::cout << "Taus    px = " << promptTausLV.Px() 
		  << "        py = " << promptTausLV.Py()
		  << "        pz = " << promptTausLV.Pz() << std::endl;
	std::cout << "VisTaus px = " << promptVisTausLV.Px() 
		  << "        py = " << promptVisTausLV.Py()
		  << "        pz = " << promptVisTausLV.Pz() << std::endl;
	std::cout << "Muons   px = " << promptMuonsLV.Px() 
		  << "        py = " << promptMuonsLV.Py()
		  << "        pz = " << promptMuonsLV.Pz() << std::endl;
	std::cout << "Elect.  px = " << promptElectronsLV.Px() 
		  << "        py = " << promptElectronsLV.Py()
		  << "        pz = " << promptElectronsLV.Pz() << std::endl;
	std::cout << "Neut.   pz = " << promptNeutrinosLV.Px()
		  << "        py = " << promptNeutrinosLV.Py()
		  << "        pz = " << promptNeutrinosLV.Pz() << std::endl;
	std::cout << "Full V  px = " << bosonPx 
		  << "        py = " << bosonPy
		  << "        pz = " << bosonPz << std::endl;
	std::cout << "Vis V   px = " << lepPx 
		  << "        py = " << lepPy
		  << "        pz = " << lepPz << std::endl;
		  std::cout << std::endl; */

	ALL->Fill(0.0);
	if (isZMM) ZMM->Fill(0.);
	if (isZEE) ZEE->Fill(0.);
	if (isZTT) ZTT->Fill(0.);
	isZLL = isZMM || isZEE;
      }
      
      //      if (isZEE&&isZMM) {
      //	cout << "Warning : isZEE && isZMM simultaneously" << endl;
      //	cout << "Mu+:" << isPrompMuPlus
      //	     << "  Mu-:"  << isPrompMuMinus
      //	     << "  E+:" << isPrompElePlus
      //	     << "  E-:" << isPrompEleMinus << endl;
      //      }
      //      if ((isPrompElePlus&&!isPrompEleMinus)||(!isPrompElePlus&&isPrompEleMinus))
      //	cout << "Warning : only one prompt electron!" << endl;
      //      if ((isPrompMuPlus&&!isPrompMuMinus)||(!isPrompMuPlus&&isPrompMuMinus))
      //	cout << "Warning : only one prompt muon!" << endl;

      run = int(analysisTree.event_run);
      lumi = int(analysisTree.event_luminosityblock);
      evt = int(analysisTree.event_nr);

      if (isData && applyGoodRunSelection) {
        bool lumi = false;
        int n=analysisTree.event_run;
        int lum = analysisTree.event_luminosityblock;

        std::string num = std::to_string(n);
        std::string lnum = std::to_string(lum);
        for(const auto& a : periods)
          {

            if ( num.c_str() ==  a.name ) {
              for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
                if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
              }
              auto last = std::prev(a.ranges.end());
              if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
            }

          }
        if (!lumi) continue;

      }

      // weights
      mcweight = analysisTree.genweight;
      puweight = 1;
      trigweight_1 = 1;
      trigweight_2 = 1;
      trigweight = 1;
      idweight_1 = 1;
      idweight_2 = 1;
      isoweight_1 = 1;
      isoweight_2 = 1;
      effweight = 1;
      fakeweight = 1;
      embeddedWeight = 1;
      signalWeight = 1;
      topptweight = 1;
      qcdweight = 1;
      qcdweightup = 1;
      qcdweightdown = 1;
      qcdweight_nodzeta = 1;
      qcdweightup_nodzeta = 1;
      qcdweightdown_nodzeta = 1;
      weight = 1;

      npv = analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
      rho = analysisTree.rho;

      npartons = analysisTree.genparticles_noutgoing;
      
      if (!isData) {
	puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      //      cout << "puweight = " << puweight << endl;
	if (topPt>0&&antitopPt>0) {
	  topptweight = topPtWeight(topPt,antitopPt);
	  //	  std::cout << "topPt = " << topPt 
	  //		    << "   antitopPt = " << antitopPt 
	  //		    << "   weight = " << topptweight << std::endl;
	}
	histWeightsSkimH->Fill(double(0),double(mcweight));
	histWeightsTTH->Fill(double(0),double(mcweight*topptweight));
      }
      else {
	histWeightsSkimH->Fill(double(0),double(1));
	histWeightsTTH->Fill(double(0),double(1));
      }

      unsigned int nLowPtLegElectron = 0;
      bool isLowPtLegElectron = false;
      
      unsigned int nHighPtLegElectron = 0;
      bool isHighPtLegElectron = false;

      unsigned int nLowPtLegMuon = 0;
      bool isLowPtLegMuon = false;
      
      unsigned int nHighPtLegMuon = 0;
      bool isHighPtLegMuon = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //      std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==LowPtLegElectron) {
	  nLowPtLegElectron = i;
	  isLowPtLegElectron = true;
	}
	if (HLTFilter==HighPtLegElectron) {
	  nHighPtLegElectron = i;
	  isHighPtLegElectron = true;
	}
	if (HLTFilter==LowPtLegMuon) {
	  nLowPtLegMuon = i;
	  isLowPtLegMuon = true;
	}
	if (HLTFilter==HighPtLegMuon) {
	  nHighPtLegMuon = i;
	  isHighPtLegMuon = true;
	}
      }
      if (!isLowPtLegElectron) {
	std::cout << "HLT filter " << LowPtLegElectron << " not found" << std::endl;
	exit(-1);
      }
      if (!isHighPtLegElectron) {
	std::cout << "HLT filter " << HighPtLegElectron << " not found" << std::endl;
	exit(-1);
      }
      if (!isLowPtLegMuon) {
	std::cout << "HLT filter " << LowPtLegMuon << " not found" << std::endl;
	exit(-1);
      }
      if (!isHighPtLegMuon) {
	std::cout << "HLT filter " << HighPtLegMuon << " not found" << std::endl;
	exit(-1);
      }
      unsigned int nBTagDiscriminant = 0;
      for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
	TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
	if (discr.Contains(BTagDiscriminator))
	  nBTagDiscriminant = iBTag;
      }


      //      std::cout << "LowPtE  : " << LowPtLegElectron << " : " << nLowPtLegElectron << std::endl;
      //      std::cout << "HighPtE : " << HighPtLegElectron << " : " << nHighPtLegElectron << std::endl;
      //      std::cout << "LowPtM  : " << LowPtLegMuon << " : " << nLowPtLegMuon << std::endl;
      //      std::cout << "HighPtM : " << HighPtLegMuon << " : " << nHighPtLegMuon << std::endl;
      //      std::cout << std::endl;
      //      continue;

      // vertex cuts no longer required
      //      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      //      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      //      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
      //		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      //      if (dVertex>dVertexCut) continue;


      //                         HLT_Mu17_Ele12                     ||      HLT_Mu8_Ele17       
      //      bool trigAccept = (analysisTree.hltriggerresults_second[5]==1)||(analysisTree.hltriggerresults_second[6]==1);
      //      if (!trigAccept) continue;

      //      if (nonOverlap&&checkOverlap)
      //      cout << "Muons = " << analysisTree.muon_count 
      //	   << "   Electrons = " << analysisTree.electron_count << std::endl;
      
      // electron selection
      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	//	bool electronMvaId = electronMvaIdWP80(analysisTree.electron_pt[ie],
	//					       analysisTree.electron_superclusterEta[ie],
	//					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
	bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie];
	if (!electronMvaId&&applyElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyElectronId) continue;
	electrons.push_back(ie);
      }

      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	muons.push_back(im);
      }

      //      cout << "  SelEle=" << electrons.size() 
      //	   << "  SelMu=" << muons.size() << std::endl;

      if (electrons.size()==0) continue;
      if (muons.size()==0) continue;
      
      // selecting muon and electron pair (OS or SS);
      int electronIndex = -1;
      int muonIndex = -1;

      float isoMuMin = 1e+10;
      float isoEleMin = 1e+10;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  electrons = " << electrons.size() << std::endl;
      bool isMuon17matched = false;
      bool isMuon8matched  = false;
      bool isElectron17matched = false;
      bool isElectron12matched = false;
      for (unsigned int im=0; im<muons.size(); ++im) {
	bool isMu17 = false;
	bool isMu8 = false;
	unsigned int mIndex  = muons.at(im);
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
	float photonIsoMu = analysisTree.muon_photonIso[mIndex];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
	float puIsoMu = analysisTree.muon_puIso[mIndex];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
	}
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig<deltaRTrigMatch) {
	    if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]&&
		analysisTree.muon_pt[mIndex]>ptMuonHighCut) { // Mu17 Leg
	      isMu17 = true;
	    }
	    if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]&&
		analysisTree.muon_pt[mIndex]>ptMuonLowCut) { // Mu8 Leg
	      isMu8 = true;
	    }
	  }
	}
	
	if (applyTriggerMatch && (!isMu17) && (!isMu8)) continue;

	for (unsigned int ie=0; ie<electrons.size(); ++ie) {

	  unsigned int eIndex = electrons.at(ie);

	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;

	  bool isEle17 = false;
	  bool isEle12 = false;

	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      if (analysisTree.trigobject_filters[iT][nHighPtLegElectron]&&
		  analysisTree.electron_pt[eIndex]>ptElectronHighCut) { // Ele17 Leg
		isEle17 = true;
	      }
	      if (analysisTree.trigobject_filters[iT][nLowPtLegElectron]&&
		  analysisTree.electron_pt[eIndex]>ptElectronLowCut) { // Ele12 Leg
		isEle12 = true;
	      }
	    }
	  }
	  
	  bool trigMatch = (isMu17&&isEle12) || (isMu8&&isEle17);
	  //	  std::cout << "Trigger match = " << trigMatch << std::endl;

	  if (applyTriggerMatch && !trigMatch) continue;
	  
	  //	  bool isKinematicMatch = false;
	  //	  if (isMu17&&isEle12) {
	  //	    if (analysisTree.muon_pt[mIndex]>ptMuonHighCut&&analysisTree.electron_pt[eIndex]>ptElectronLowCut)
	  //	      isKinematicMatch = true;
	  //	  }
	  //	  if (isMu8&&isEle17) {
	  //            if (analysisTree.muon_pt[mIndex]>ptMuonLowCut&&analysisTree.electron_pt[eIndex]>ptElectronHighCut)
	  //              isKinematicMatch = true;
	  //          }
	  //	  if (!isKinematicMatch) continue;

	  float neutralHadIsoEle = analysisTree.electron_neutralHadIso[eIndex];
	  float photonIsoEle = analysisTree.electron_photonIso[eIndex];
	  float chargedHadIsoEle = analysisTree.electron_chargedHadIso[eIndex];
	  float puIsoEle = analysisTree.electron_puIso[eIndex];
	  if (isIsoR03) {
	    neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
	    photonIsoEle = analysisTree.electron_r03_sumPhotonEt[eIndex];
	    chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
	    puIsoEle = analysisTree.electron_r03_sumPUPt[eIndex];
	  }
	  float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	  neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
	  float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	  float relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];

	  if (int(mIndex)!=muonIndex) {
	    if (relIsoMu==isoMuMin) {
	      if (analysisTree.muon_pt[mIndex]>analysisTree.muon_pt[muonIndex]) {
		isoMuMin  = relIsoMu;
		muonIndex = int(mIndex);
		isoEleMin = relIsoEle;
		electronIndex = int(eIndex);
		isMuon17matched = isMu17;
		isMuon8matched = isMu8;
		isElectron17matched = isEle17;
		isElectron12matched = isEle12;
	      }
	    }
	    else if (relIsoMu<isoMuMin) {
	      isoMuMin  = relIsoMu;
	      muonIndex = int(mIndex);
	      isoEleMin = relIsoEle;
	      electronIndex = int(eIndex);
	      isMuon17matched = isMu17;
	      isMuon8matched = isMu8;
	      isElectron17matched = isEle17;
	      isElectron12matched = isEle12;
	    }
	  }
	  else {
	    if (relIsoEle==isoEleMin) {
	      if (analysisTree.electron_pt[eIndex]>analysisTree.electron_pt[electronIndex]) {
		isoEleMin = relIsoEle;
		electronIndex = int(eIndex);
		isElectron17matched = isEle17;
		isElectron12matched = isEle12;
	      }
	    }
	    else if (relIsoEle<isoEleMin) {
	      isoEleMin = relIsoEle;
	      electronIndex = int(eIndex);
	      isElectron17matched = isEle17;
	      isElectron12matched = isEle12;
	    }
	  }
	  
	}
      }

      //      cout << "mIndex = " << muonIndex << "   eIndex = " << electronIndex << std::endl;

      if (electronIndex<0) continue;
      if (muonIndex<0) continue;
      //      std::cout << "Post synch selection " << std::endl;
      //      std::cout << std::endl;

      

      os = (analysisTree.muon_charge[muonIndex]*analysisTree.electron_charge[electronIndex]) < 0;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (int(ie)==electronIndex) continue;
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	//	bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[ie],
	//					       analysisTree.electron_superclusterEta[ie],
	//					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
	bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
	if (!electronMvaId&&applyVetoElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId&&applyVetoElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId&&applyVetoElectronId) continue;
	float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
	float photonIsoEle = analysisTree.electron_photonIso[ie];
	float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
	float puIsoEle = analysisTree.electron_puIso[ie];
	if (isIsoR03) {
	  neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
	  photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
	  chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
	  puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
	}
	float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
	float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];	
	if (relIsoEle>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      // looking for extra muon
      bool foundExtraMuon = false;
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (int(im)==muonIndex) continue;
	if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	float photonIsoMu = analysisTree.muon_photonIso[im];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	float puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[im];
	}
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      dilepton_veto = false;
      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

      if (applyInclusiveSelection) {
	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
	if (isoMuMin>isoMuonHighCut) continue;
	if (isoEleMin>isoElectronHighCut) continue;
      }

      //      cout << "dilepton_veto : " << dilepton_veto
      //      	   << "   extraelec_veto : " << extraelec_veto
      //      	   << "   extramuon_veto : " << extramuon_veto << endl;

      // filling muon variables
      pt_2 = analysisTree.muon_pt[muonIndex];
      pt_Up_2   = (1+muonScale)*pt_2;
      pt_Down_2 = (1-muonScale)*pt_2;
     

      eta_2 = analysisTree.muon_eta[muonIndex];
      phi_2 = analysisTree.muon_phi[muonIndex];
      q_2 = -1;
      if (analysisTree.muon_charge[muonIndex]>0)
	q_2 = 1;
      mva_2 = -9999;
      d0_2 = analysisTree.muon_dxy[muonIndex];
      dZ_2 = analysisTree.muon_dz[muonIndex];
      iso_2 = isoMuMin;
      m_2 = muonMass;
      
      float eleScale = eleScaleBarrel;
      if (fabs(analysisTree.electron_eta[electronIndex])>1.479) eleScale = eleScaleEndcap;
      // filling electron variables
      pt_1 = analysisTree.electron_pt[electronIndex];
      pt_Up_1 = (1+eleScale)*pt_1;
      pt_Down_1 = (1-eleScale)*pt_1;

      eta_1 = analysisTree.electron_eta[electronIndex];
      phi_1 = analysisTree.electron_phi[electronIndex];
      q_1 = -1;
      if (analysisTree.electron_charge[electronIndex]>0)
        q_1 = 1;
      mva_1 = analysisTree.electron_mva_id_nontrigPhys14[electronIndex];
      d0_1 = analysisTree.electron_dxy[electronIndex];
      dZ_1 = analysisTree.electron_dz[electronIndex];
      iso_1 = isoEleMin;
      m_1 = electronMass;


      // scale factors
      isoweight_1 = (float)SF_electronIdIso->get_ScaleFactor(double(pt_1),double(eta_1));
      isoweight_2 = (float)SF_muonIdIso->get_ScaleFactor(double(pt_2),double(eta_2));

      //      cout << "isoweight_1 = " << isoweight_1
      //	   << "isoweight_2 = " << isoweight_2 << endl;

      float Ele17EffData = (float)SF_electron17->get_EfficiencyData(double(pt_1),double(eta_1));
      float Ele17EffMC   = (float)SF_electron17->get_EfficiencyMC(double(pt_1),double(eta_1));

      float Ele12EffData = (float)SF_electron12->get_EfficiencyData(double(pt_1),double(eta_1));
      float Ele12EffMC   = (float)SF_electron12->get_EfficiencyMC(double(pt_1),double(eta_1));

      float Mu17EffData = (float)SF_muon17->get_EfficiencyData(double(pt_2),double(eta_2));
      float Mu17EffMC   = (float)SF_muon17->get_EfficiencyMC(double(pt_2),double(eta_2));

      float Mu8EffData = (float)SF_muon8->get_EfficiencyData(double(pt_2),double(eta_2));
      float Mu8EffMC   = (float)SF_muon8->get_EfficiencyMC(double(pt_2),double(eta_2));

      float trigWeightData = Mu17EffData*Ele12EffData + Mu8EffData*Ele17EffData - Mu17EffData*Ele17EffData;
      float trigWeightMC   = Mu17EffMC*Ele12EffMC     + Mu8EffMC*Ele17EffMC     - Mu17EffMC*Ele17EffMC;

      if (isMuon17matched && isElectron12matched) {
	trigweight_1 = (float)SF_electron12->get_ScaleFactor(double(pt_1),double(eta_1));
	trigweight_2 = (float)SF_muon17->get_ScaleFactor(double(pt_2),double(eta_2));
      }
      else if (isMuon8matched && isElectron17matched) {
	trigweight_1 = (float)SF_electron17->get_ScaleFactor(double(pt_1),double(eta_1));
	trigweight_2 = (float)SF_muon8->get_ScaleFactor(double(pt_2),double(eta_2));
      }
	
      if (trigWeightMC>1e-6)
	trigweight = trigWeightData / trigWeightMC;

      effweight = isoweight_1*isoweight_2*trigweight;

      //      cout << "effweight = " << effweight << endl;

      // dilepton system
      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector muonUpLV; muonUpLV.SetXYZM((1.0+muonScale)*analysisTree.muon_px[muonIndex],
						(1.0+muonScale)*analysisTree.muon_py[muonIndex],
						(1.0+muonScale)*analysisTree.muon_pz[muonIndex],
						muonMass);

      TLorentzVector muonDownLV; muonDownLV.SetXYZM((1.0-muonScale)*analysisTree.muon_px[muonIndex],
						    (1.0-muonScale)*analysisTree.muon_py[muonIndex],
						    (1.0-muonScale)*analysisTree.muon_pz[muonIndex],
						    muonMass);

      TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
						    analysisTree.electron_py[electronIndex],
						    analysisTree.electron_pz[electronIndex],
						    electronMass);

      TLorentzVector electronUpLV; electronUpLV.SetXYZM((1.0+eleScale)*analysisTree.electron_px[electronIndex],
							(1.0+eleScale)*analysisTree.electron_py[electronIndex],
							(1.0+eleScale)*analysisTree.electron_pz[electronIndex],
							electronMass);

      TLorentzVector electronDownLV; electronDownLV.SetXYZM((1.0-eleScale)*analysisTree.electron_px[electronIndex],
							    (1.0-eleScale)*analysisTree.electron_py[electronIndex],
							    (1.0-eleScale)*analysisTree.electron_pz[electronIndex],
							    electronMass);

      TLorentzVector dileptonLV = muonLV + electronLV;

      m_vis = dileptonLV.M();

      m_vis_muUp    = (muonUpLV+electronLV).M();
      m_vis_muDown  = (muonDownLV+electronLV).M();
      m_vis_eUp     = (muonLV+electronUpLV).M();
      m_vis_eDown   = (muonLV+electronDownLV).M();
      m_vis_scaleUp   = m_vis;
      m_vis_scaleDown = m_vis;
      m_vis_resoUp    = m_vis;
      m_vis_resoDown  = m_vis;

      // std::cout << "m_vis        = " << m_vis << std::endl;
      // std::cout << "m_vis_muUp   = " << m_vis_muUp << std::endl;
      // std::cout << "m_vis_muDown = " << m_vis_muDown << std::endl;
      // std::cout << "m_vis_eUp    = " << m_vis_eUp << std::endl;
      // std::cout << "m_vis_eDown  = " << m_vis_eDown << std::endl;
      // std::cout << std::endl;

      pt_tt = dileptonLV.Pt();

      dphi_tt = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
			   electronLV.Px(),electronLV.Py());

      dr_tt = deltaR(muonLV.Eta(),muonLV.Phi(),
		     electronLV.Eta(),electronLV.Phi());

      // qcd scale factor
      // no dzeta cut
      qcdweight     = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
      qcdweightup   = qcdWeight.getWeightUp(pt_1,pt_2,dr_tt);
      qcdweightdown = qcdWeight.getWeightDown(pt_1,pt_2,dr_tt);
      // dzeta cut
      qcdweight_nodzeta     = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt); 
      qcdweightup_nodzeta   = qcdWeightNoDzeta.getWeightUp(pt_1,pt_2,dr_tt);
      qcdweightdown_nodzeta = qcdWeightNoDzeta.getWeightDown(pt_1,pt_2,dr_tt);

      //      if (os<0.5) {
      //      	printf("QCD weights  : pt_1 = %6.1f ; pt_2 = %6.1f ; dr_tt = %4.2f\n",pt_1,pt_2,dr_tt);
      //      	printf("DZeta cut    : central = %4.2f ; up = %4.2f ; down = %4.2f\n",qcdweight,qcdweightup,qcdweightdown);
      //      	printf("No dZeta cut : central = %4.2f ; up = %4.2f ; down = %4.2f\n",qcdweight_nodzeta,qcdweightup_nodzeta,qcdweightdown_nodzeta);
      //      std::cout << std::endl;
      //     }
      
      // counting jets
      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;
      float ptLeadingBJet = -1;

      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	float jetEta = analysisTree.pfjet_eta[jet];
	if (absJetEta>jetEtaCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];
	if (jetPt<jetPtLowCut) continue;

	float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   eta_1,phi_1);
	if (dR1<dRJetLeptonCut) continue;

	float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           eta_2,phi_2);
        if (dR2<dRJetLeptonCut) continue;

	// jetId

	bool isPFJetId = looseJetiD(analysisTree,int(jet));
	if (!isPFJetId) continue;

	jetspt20.push_back(jet);

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

	  bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet

	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);

	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;

	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_BC.eval(BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_BC.eval(BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader_Light.eval(BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
	      tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
	    }
	    
	    if (tageff<1e-5)      tageff = 1e-5;
	    if (tageff>0.99999)   tageff = 0.99999;
	    rand.SetSeed((int)((jetEta+5)*100000));
	    double rannum = rand.Rndm();
	    
	    if (jet_scalefactor<1 && tagged) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		tagged = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !tagged) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) { 
		tagged = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }
	  }

	  if (tagged) {
	    bjets.push_back(jet);
	    if (jetPt>ptLeadingBJet) {
	      ptLeadingBJet = jetPt;
	      indexLeadingBJet = jet;
	    }
	  }
	}

	if (jetPt>jetPtHighCut)
	  jets.push_back(jet);

	if (indexLeadingJet>=0) {
	  if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
	    indexSubLeadingJet = jet;
	    ptSubLeadingJet = jetPt;
	  }
	}

	if (jetPt>ptLeadingJet) {
	  indexLeadingJet = jet;
	  ptLeadingJet = jetPt;
	}
      }
	
      njets = jets.size();
      njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      
      bpt = -9999;
      beta = -9999;
      bphi = -9999;
      
      if (indexLeadingBJet>=0) {
	bpt = analysisTree.pfjet_pt[indexLeadingBJet];
	beta = analysisTree.pfjet_eta[indexLeadingBJet];
	bphi = analysisTree.pfjet_phi[indexLeadingBJet];
      }
     
      jpt_1 = -9999;
      jeta_1 = -9999;
      jphi_1 = -9999;
      jptraw_1 = -9999;
      jptunc_1 = -9999;
      jmva_1 = -9999;
      jlrm_1 = -9999;
      jctm_1 = -9999;

      if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
	cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;

      if (indexLeadingJet>=0) {
	jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
	jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
	jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
	jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
	jmva_1 = analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
      }

      jpt_2 = -9999;
      jeta_2 = -9999;
      jphi_2 = -9999;
      jptraw_2 = -9999;
      jptunc_2 = -9999;
      jmva_2 = -9999;
      jlrm_2 = -9999;
      jctm_2 = -9999;

      if (indexSubLeadingJet>=0) {
	jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
	jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
	jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
	jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];
	jmva_2 = analysisTree.pfjet_pu_jet_full_mva[indexSubLeadingJet];
      }

      mjj =  -9999;
      jdeta =  -9999;
      njetingap = 0;

      if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {

	TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
					     analysisTree.pfjet_py[indexLeadingJet],
					     analysisTree.pfjet_pz[indexLeadingJet],
					     analysisTree.pfjet_e[indexLeadingJet]);

	TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
					     analysisTree.pfjet_py[indexSubLeadingJet],
					     analysisTree.pfjet_pz[indexSubLeadingJet],
					     analysisTree.pfjet_e[indexSubLeadingJet]);

	mjj = (jet1+jet2).M();
	jdeta = abs(analysisTree.pfjet_eta[indexLeadingJet]-
		    analysisTree.pfjet_eta[indexSubLeadingJet]);
 
	float etamax = analysisTree.pfjet_eta[indexLeadingJet];
	float etamin = analysisTree.pfjet_eta[indexSubLeadingJet];
	if (etamax<etamin) {
	  float tmp = etamax;
	  etamax = etamin;
	  etamin = tmp;
	}
	for (unsigned int jet=0; jet<jetspt20.size(); ++jet) {
	  int index = jetspt20.at(jet);
	  float etaX = analysisTree.pfjet_eta[index];
	  if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax) 
	    njetingap++;
	}


      }

      // METs
      met = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      metphi = TMath::ATan2(analysisTree.pfmet_ey,analysisTree.pfmet_ex);
      metcov00 = analysisTree.pfmet_sigxx;
      metcov01 = analysisTree.pfmet_sigxy;
      metcov10 = analysisTree.pfmet_sigyx;
      metcov11 = analysisTree.pfmet_sigyy;

      //      choosing mva met
      unsigned int metEMu = 0;
      bool mvaMetFound = false;
      //      cout << "MVA MET : " << analysisTree.mvamet_count << endl;
      for (unsigned int iMet=0; iMet<analysisTree.mvamet_count; ++iMet) {
	//      	cout << iMet << endl;
	if (analysisTree.mvamet_channel[iMet]==1) {
      	  if (int(analysisTree.mvamet_lep1[iMet])==muonIndex&&
      	      int(analysisTree.mvamet_lep2[iMet])==electronIndex) {
      	    metEMu = iMet;
      	    mvaMetFound = true;
      	  }
      	}
      }
      if (!mvaMetFound) {
	cout << "Warning : mva Met is not found..." << endl;
      }
      
      float mvamet_x = 0;
      float mvamet_y = 0;
      mvamet = 0;
      mvametphi = 0;
      mvacov00 = 10;
      mvacov01 = 10;
      mvacov10 = 10;
      mvacov11 = 10;
      if (analysisTree.mvamet_count>0) {
      	mvamet_x = analysisTree.mvamet_ex[metEMu];
      	mvamet_y = analysisTree.mvamet_ey[metEMu];
      	float mvamet_x2 = mvamet_x * mvamet_x;
      	float mvamet_y2 = mvamet_y * mvamet_y;

      	mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
      	mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
      	mvacov00 = analysisTree.mvamet_sigxx[metEMu];
      	mvacov01 = analysisTree.mvamet_sigxy[metEMu];
      	mvacov10 = analysisTree.mvamet_sigyx[metEMu];
      	mvacov11 = analysisTree.mvamet_sigyy[metEMu];
      }

      // Recoil corrections

      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;

      float mvamet_corr_x = mvamet_x;
      float mvamet_corr_y = mvamet_y;

      float mvamet_uncorr_x = mvamet_x;
      float mvamet_uncorr_y = mvamet_y;

      mvamet_uncorr = mvamet;
      mvametphi_uncorr = mvametphi;

      if ((isW||isDY)&&!isData)
	recoilMvaMetCorrector.CorrectByMeanResolution(mvamet_x,mvamet_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,mvamet_corr_x,mvamet_corr_y);

      //      printf("Before correction : mvamet_x = %7.2f   mvamet_y = %7.2f\n",mvamet_x,mvamet_y);
      //      printf("After correction  : mvamet_x = %7.2f   mvamet_y = %7.2f\n",mvamet_corr_x,mvamet_corr_y);
      //      std::cout << std::endl;

      mvamet_x = mvamet_corr_x;
      mvamet_y = mvamet_corr_y;

      mvamet = TMath::Sqrt(mvamet_x*mvamet_x+mvamet_y*mvamet_y); 
      mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
      
      // MEt related systematic uncertainties
      int bkgdType = 0;
      if (isDY||isW)
	bkgdType = MEtSys::ProcessType::BOSON;
      else if (isTOP)
	bkgdType = MEtSys::ProcessType::TOP;
      else 
	bkgdType = MEtSys::ProcessType::EWK; 

      float mvamet_scaleUp_x   = mvamet_x;
      float mvamet_scaleUp_y   = mvamet_y;
      float mvamet_scaleDown_x = mvamet_x;
      float mvamet_scaleDown_y = mvamet_y;
      float mvamet_resoUp_x    = mvamet_x;
      float mvamet_resoUp_y    = mvamet_y;
      float mvamet_resoDown_x  = mvamet_x;
      float mvamet_resoDown_y  = mvamet_y;

      if (!isData) {
	metSys.ApplyMEtSys(mvamet_x,mvamet_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Up,
			   mvamet_scaleUp_x,mvamet_scaleUp_y);
	metSys.ApplyMEtSys(mvamet_x,mvamet_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Down,
			   mvamet_scaleDown_x,mvamet_scaleDown_y);
	metSys.ApplyMEtSys(mvamet_x,mvamet_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Up,
			   mvamet_resoUp_x,mvamet_resoUp_y);
	metSys.ApplyMEtSys(mvamet_x,mvamet_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Down,
			   mvamet_resoDown_x,mvamet_resoDown_y);
      }

      
      mvamet_scaleUp = TMath::Sqrt(mvamet_scaleUp_x*mvamet_scaleUp_x+
				   mvamet_scaleUp_y*mvamet_scaleUp_y);
      mvametphi_scaleUp = TMath::ATan2(mvamet_scaleUp_y,mvamet_scaleUp_x);
      
      mvamet_scaleDown = TMath::Sqrt(mvamet_scaleDown_x*mvamet_scaleDown_x+
				     mvamet_scaleDown_y*mvamet_scaleDown_y);
      mvametphi_scaleDown = TMath::ATan2(mvamet_scaleDown_y,mvamet_scaleDown_x);
      
      mvamet_resoUp = TMath::Sqrt(mvamet_resoUp_x*mvamet_resoUp_x+
				  mvamet_resoUp_y*mvamet_resoUp_y);
      mvametphi_resoUp = TMath::ATan2(mvamet_resoUp_y,mvamet_resoUp_x);
      
      mvamet_resoDown = TMath::Sqrt(mvamet_resoDown_x*mvamet_resoDown_x+
				    mvamet_resoDown_y*mvamet_resoDown_y);
      mvametphi_resoDown = TMath::ATan2(mvamet_resoDown_y,mvamet_resoDown_x);


      //      printf("Uncorrected    :  %7.2f  %7.2f\n",mvamet_uncorr_x,mvamet_uncorr_y);
      //      printf("Central        :  %7.2f  %7.2f\n",mvamet_x,mvamet_y);
      //      printf("ScaleUp        :  %7.2f  %7.2f\n",mvamet_scaleUp_x,mvamet_scaleUp_y);
      //      printf("ScaleDown      :  %7.2f  %7.2f\n",mvamet_scaleDown_x,mvamet_scaleDown_y);
      //      printf("ResoUp         :  %7.2f  %7.2f\n",mvamet_resoUp_x,mvamet_resoUp_y);
      //      printf("ResoDown       :  %7.2f  %7.2f\n",mvamet_resoDown_x,mvamet_resoDown_y);
      //      std::cout << std::endl;


      genmet = TMath::Sqrt(analysisTree.genmet_ex*analysisTree.genmet_ex + analysisTree.genmet_ey*analysisTree.genmet_ey);
      genmetphi = TMath::ATan2(analysisTree.genmet_ey,analysisTree.genmet_ex);


      // bisector of electron and muon transverse momenta
      float electronUnitX = electronLV.Px()/electronLV.Pt();
      float electronUnitY = electronLV.Py()/electronLV.Pt();
	
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();

      float zetaX = electronUnitX + muonUnitX;
      float zetaY = electronUnitY + muonUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + electronLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + electronLV.Py();
      
      float vectorVisX = muonLV.Px() + electronLV.Px();
      float vectorVisY = muonLV.Py() + electronLV.Py();

      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;

      // computation of DZeta variable
      // pfmet
      computeDzeta(analysisTree.pfmet_ex,analysisTree.pfmet_ey,
		   zetaX,zetaY,pzetavis,pzetamiss,dzeta);

      // mvamet
      computeDzeta(mvamet_x,mvamet_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet,dzeta_mvamet); // corrected
      computeDzeta(mvamet_uncorr_x,mvamet_uncorr_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet_uncorr,dzeta_mvamet_uncorr); // uncorrected
      computeDzeta(mvamet_scaleUp_x,mvamet_scaleUp_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet_scaleUp,dzeta_mvamet_scaleUp); // scaleUp
      computeDzeta(mvamet_scaleDown_x,mvamet_scaleDown_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet_scaleDown,dzeta_mvamet_scaleDown); // scaleDown
      computeDzeta(mvamet_resoUp_x,mvamet_resoUp_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet_resoUp,dzeta_mvamet_resoUp); // resoUp
      computeDzeta(mvamet_resoDown_x,mvamet_resoDown_y,
		   zetaX,zetaY,pzetavis,pzetamiss_mvamet_resoDown,dzeta_mvamet_resoDown); // resoDown
       // genmet
      computeDzeta(analysisTree.genmet_ex,analysisTree.genmet_ey,
		   zetaX,zetaY,pzetavis,pzetamiss_genmet,dzeta_genmet);

      TLorentzVector mvametLV; mvametLV.SetXYZM(mvamet_x,mvamet_y,0.,0.);
      TLorentzVector mvametResoUpLV; mvametResoUpLV.SetXYZM(mvamet_resoUp_x,mvamet_resoUp_y,0.,0.);
      TLorentzVector mvametResoDownLV; mvametResoDownLV.SetXYZM(mvamet_resoDown_x,mvamet_resoDown_y,0.,0.);
      TLorentzVector mvametScaleUpLV; mvametScaleUpLV.SetXYZM(mvamet_scaleUp_x,mvamet_scaleUp_y,0.,0.);
      TLorentzVector mvametScaleDownLV; mvametScaleDownLV.SetXYZM(mvamet_scaleDown_x,mvamet_scaleDown_y,0.,0.);

      // computing total transverse mass
      mTtot = totalTransverseMass        ( muonLV ,     electronLV , mvametLV);

      mTtot_muUp   =  totalTransverseMass( muonUpLV ,   electronLV , mvametLV);
      mTtot_muDown =  totalTransverseMass( muonDownLV , electronLV , mvametLV);

      mTtot_eUp   =  totalTransverseMass ( muonLV , electronUpLV   , mvametLV);
      mTtot_eDown =  totalTransverseMass ( muonLV , electronDownLV , mvametLV);

      mTtot_scaleUp   =  totalTransverseMass ( muonLV , electronLV , mvametScaleUpLV);
      mTtot_scaleDown =  totalTransverseMass ( muonLV , electronLV , mvametScaleDownLV);

      mTtot_resoUp   =  totalTransverseMass ( muonLV , electronLV , mvametResoUpLV);
      mTtot_resoDown =  totalTransverseMass ( muonLV , electronLV , mvametResoDownLV);

      m_sv           = -9999;
      m_sv_muUp      = -9999;
      m_sv_muDown    = -9999;
      m_sv_eUp       = -9999;
      m_sv_eDown     = -9999;
      m_sv_scaleUp   = -9999;
      m_sv_scaleDown = -9999;
      m_sv_resoUp    = -9999;
      m_sv_resoDown  = -9999;

      mt_sv           = -9999;
      mt_sv_muUp      = -9999;
      mt_sv_muDown    = -9999;
      mt_sv_eUp       = -9999;
      mt_sv_eDown     = -9999;
      mt_sv_scaleUp   = -9999;
      mt_sv_scaleDown = -9999;
      mt_sv_resoUp    = -9999;
      mt_sv_resoDown  = -9999;

      if (computeSVFitMass && dzeta_mvamet>-30) {

	if (mvaMetFound) {
	  // covariance matrix MVAMET
	  TMatrixD covMVAMET(2, 2);
	  covMVAMET[0][0] =  mvacov00;
	  covMVAMET[1][0] =  mvacov01;
	  covMVAMET[0][1] =  mvacov10;
	  covMVAMET[1][1] =  mvacov11;

	  // mva met
	  // define electron 4-vector
	  svFitStandalone::MeasuredTauLepton svFitEle(svFitStandalone::kTauToElecDecay, 
						       pt_1, eta_1, phi_1, 0.51100e-3); 
	  svFitStandalone::MeasuredTauLepton svFitEleUp(svFitStandalone::kTauToElecDecay, 
							(1.0+eleScale)*pt_1, eta_1, phi_1, 0.51100e-3); 
	  svFitStandalone::MeasuredTauLepton svFitEleDown(svFitStandalone::kTauToElecDecay, 
							  (1.0-eleScale)*pt_1, eta_1, phi_1, 0.51100e-3); 
	  // define muon 4-vector
	  svFitStandalone::MeasuredTauLepton svFitMu(svFitStandalone::kTauToMuDecay,  
						      pt_2, eta_2, phi_2, 105.658e-3); 

  
	  // central value
	  SVfitStandaloneAlgorithm algo = SVFitMassComputation(svFitEle, svFitMu,
							       mvamet_x, mvamet_y, 
							       covMVAMET, inputFile_visPtResolution);

	  m_sv = algo.getMass(); // return value of svfit mass is in units of GeV
	  mt_sv = algo.transverseMass(); // return value of transverse svfit mass is in units of GeV
	  if ( !algo.isValidSolution() ) 
	    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
	  
	  pt_sv  = algo.pt(); 
	  eta_sv = algo.eta();
	  phi_sv = algo.phi();

	  m_sv_scaleUp   = m_sv;
	  m_sv_scaleDown = m_sv;
	  m_sv_resoUp    = m_sv;
	  m_sv_resoDown  = m_sv;
	  m_sv_eUp       = m_sv;
	  m_sv_eDown     = m_sv;
	  m_sv_muUp      = m_sv;
	  m_sv_muDown    = m_sv;

	  mt_sv_scaleUp   = mt_sv;
	  mt_sv_scaleDown = mt_sv;
	  mt_sv_resoUp    = mt_sv;
	  mt_sv_resoDown  = mt_sv;
	  mt_sv_eUp       = mt_sv;
	  mt_sv_eDown     = mt_sv;
	  mt_sv_muUp      = mt_sv;
	  mt_sv_muDown    = mt_sv;

	  if (!isData) { 

	    SVfitStandaloneAlgorithm algo_eUp = SVFitMassComputation(svFitEleUp, svFitMu,
								     mvamet_x, mvamet_y,
								     covMVAMET, inputFile_visPtResolution);
	    SVfitStandaloneAlgorithm algo_eDown = SVFitMassComputation(svFitEleDown, svFitMu,
								       mvamet_x, mvamet_y,
								       covMVAMET, inputFile_visPtResolution);
	    

	    m_sv_eUp    = algo_eUp.getMass();
	    mt_sv_eUp   = algo_eUp.transverseMass();
	    m_sv_eDown  = algo_eDown.getMass();
	    mt_sv_eDown = algo_eDown.transverseMass();

	    if (isDY) {
	      SVfitStandaloneAlgorithm algo_scaleUp = SVFitMassComputation(svFitEle, svFitMu,
									   mvamet_scaleUp_x, mvamet_scaleUp_y,
									   covMVAMET, inputFile_visPtResolution);

	      SVfitStandaloneAlgorithm algo_scaleDown = SVFitMassComputation(svFitEle, svFitMu,
									     mvamet_scaleDown_x, mvamet_scaleDown_y,
									     covMVAMET, inputFile_visPtResolution);

	      // SVfitStandaloneAlgorithm algo_resoUp = SVFitMassComputation(svFitEle, svFitMu,
	      // 								mvamet_resoUp_x, mvamet_resoUp_y,
	      // 								covMVAMET, inputFile_visPtResolution);
	      
	      // SVfitStandaloneAlgorithm algo_resoDown = SVFitMassComputation(svFitEle, svFitMu,
	      // 								  mvamet_resoDown_x, mvamet_resoDown_y,
	      // 								  covMVAMET, inputFile_visPtResolution);
	      
	      m_sv_scaleUp    = algo_scaleUp.getMass();
	      mt_sv_scaleUp   = algo_scaleUp.transverseMass();
	      m_sv_scaleDown  = algo_scaleDown.getMass();
	      mt_sv_scaleDown = algo_scaleDown.transverseMass();

	      // m_sv_resoUp     = algo_resoUp.getMass();
	      // mt_sv_resoUp    = algo_resoUp.transverseMass();
	      // m_sv_resoDown   = algo_resoDown.getMass();
	      // mt_sv_resoDown  = algo_resoDown.transverseMass();

	    }
	  }
	  //	  std::cout << "eta(e) = " << eta_1 << "   escale = " << eleScale << std::endl;
	  //	  std::cout << "msv = " << m_sv << std::endl;
	  //	  std::cout << "msv(eES)      : up = " << m_sv_eUp << "   down = " << m_sv_eDown << std::endl;
	  //	  std::cout << "msv(metScale) : up = " << m_sv_scaleUp << "   down = " << m_sv_scaleDown << std::endl;
	  //	  std::cout << "mtsv = " << mt_sv << std::endl;
	  //	  std::cout << "mtsv(eES)      : up = " << mt_sv_eUp << "   down = " << mt_sv_eDown << std::endl;
	  //	  std::cout << "mtsv(metScale) : up = " << mt_sv_scaleUp << "   down = " << mt_sv_scaleDown << std::endl;
	  //	  std::cout << std::endl;
	}
      }
      gen_match_1 = 6;
      gen_match_2 = 6;

      // isZTT = false;
      // isZL  = false;
      // isZJ  = false;

      float minDR_1 = 0.2;
      float minDR_2 = 0.2;
      if (!isData && isDY) {
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);
	  float ptGen = genLV.Pt();
	  bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool isAnyType = type1 || type2 || type3 || type4;
	  if (isAnyType) {
	    float etaGen = genLV.Eta();
	    float phiGen = genLV.Phi();
	    float deltaR_1 = deltaR(eta_1,phi_1,
				    etaGen,phiGen);
	    if (deltaR_1<minDR_1) {
	      minDR_1 = deltaR_1;
	      if (type1) gen_match_1 = 1;
	      else if (type2) gen_match_1 = 2;
	      else if (type3) gen_match_1 = 3;
	      else if (type4) gen_match_1 = 4;
	    }
	    
	    float deltaR_2 = deltaR(eta_2,phi_2,
				    etaGen,phiGen);
	    if (deltaR_2<minDR_2) {
	      minDR_2 = deltaR_2;
	      if (type1) gen_match_2 = 1;
	      else if (type2) gen_match_2 = 2;
	      else if (type3) gen_match_2 = 3;
	      else if (type4) gen_match_2 = 4;
	    }
	  }
	}

	if (gen_match_1>2&&gen_match_2>3) { 
	  isZTT = true;
	  isZLL = false;
	}
	else {
	  isZTT = false;
	  isZLL = true;
	}

      }

      tree->Fill();
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
  file->Close();
  delete file;
  
}



