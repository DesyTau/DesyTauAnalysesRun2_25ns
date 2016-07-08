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

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

float totalTransverseMass(TLorentzVector l1, 
			  TLorentzVector l2,
			  TLorentzVector l3) {

  TLorentzVector totalLV = l1 + l2 + l3;
  float    totalET = l1.Pt() +  l2.Pt() + l3.Pt();
  float     mTtot = TMath::Sqrt(totalET*totalET-totalLV.Pt()*totalLV.Pt());
  return    mTtot;

}

SVfitStandaloneAlgorithm SVFitMassComputation(svFitStandalone::MeasuredTauLepton svFitMuon,
					      svFitStandalone::MeasuredTauLepton svFitTau,
					      double measuredMVAMETx,
					      double measuredMVAMETy,
					      TMatrixD covMVAMET,
					      int tau_decayMode,
					      TFile * inputFile_visPtResolution
					      ) {

  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitMuon);
  measuredTauLeptons.push_back(svFitTau);
  bool validDecayMode = tau_decayMode<=4 || tau_decayMode>=10; // excluding 2prongs
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMVAMETx, measuredMVAMETy, covMVAMET, 0);
  algo.addLogM(false);  
  if (validDecayMode)
    algo.shiftVisPt(true, inputFile_visPtResolution);
  else
    algo.shiftVisPt(false,inputFile_visPtResolution);
  algo.integrateMarkovChain();
  
  return algo;

}

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool isDY   = cfg.get<bool>("IsDY");
  const bool isW    = cfg.get<bool>("IsW");
  const bool isTOP  = cfg.get<bool>("IsTOP");

  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");
  const bool applyInclusiveSelection = cfg.get<bool>("ApplyInclusiveSelection");
  const bool computeDitauMass        = cfg.get<bool>("ComputeDitauMass");

  // kinematic cuts on electrons
  const float ptTauCut   = cfg.get<float>("ptTauCut");
  const float etaTauCut     = cfg.get<float>("etaTauCut");
  const float dzTauCut      = cfg.get<float>("dzTauCut");
  const float isoTauCut     = cfg.get<float>("isoTauCut");

  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // kinematic cuts on muons
  const float ptMuonCut       = cfg.get<float>("ptMuonCut");
  const float ptMuonHighCut   = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut      = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut      = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut       = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut   = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut  = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId      = cfg.get<bool>("ApplyMuonId");

  // HLT filters
  const string isoMuonLeg   = cfg.get<string>("IsoMuonLeg");
  const string muonTauMuonLeg = cfg.get<string>("MuonTauMuonLeg");
  const string muonTauOverlap = cfg.get<string>("MuonTauOverlap");
  const string muonTauTauLeg  = cfg.get<string>("MuonTauTauLeg");
  const float ptMuonTriggerCut = cfg.get<float>("PtMuonTriggerCut");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // 
  const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
  const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
  const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
  const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
  const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");

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

  TString IsoMuonLeg(isoMuonLeg);
  TString MuonTauMuonLeg(muonTauMuonLeg);
  TString MuonTauOverlap(muonTauOverlap);
  TString MuonTauTauLeg(muonTauTauLeg);
  TString BTagDiscriminator(bTagDiscriminator);

  const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  const string MuonTriggerFile = cfg.get<string>("MuonTriggerEff");

  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  const string recoilPuppiFileName   = cfg.get<string>("RecoilPuppiFileName");
  TString RecoilPuppiFileName(recoilPuppiFileName);

  const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
  TString RecoilMvaFileName(recoilMvaFileName);

  const string metSysFileName   = cfg.get<string>("MetSysFileName");
  TString MetSysFileName(metSysFileName);

  // hard-coded for the moment
  float tauScale = 0.03;

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
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
  TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
  TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
  TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
  TH1D * GAM = new TH1D("GAM","",1,-0.5,0.5);

  TTree * tree = new TTree("TauCheck","TauCheck");
   // Declaration of leaf types

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
  Float_t         idweight_1;
  Float_t         idweight_2;
  Float_t         isoweight_1;
  Float_t         isoweight_2;
  Float_t         effweight;
  Float_t         fakeweight;
  Float_t         embeddedWeight;
  Float_t         signalWeight;
  Float_t         weight;

  Float_t         m_vis;
  Float_t         m_vis_tauUp;
  Float_t         m_vis_tauDown;
  Float_t         m_vis_scaleUp;
  Float_t         m_vis_scaleDown;
  Float_t         m_vis_resoUp;
  Float_t         m_vis_resoDown;

  Float_t         mTtot;
  Float_t         mTtot_tauUp;
  Float_t         mTtot_tauDown;
  Float_t         mTtot_resoUp;
  Float_t         mTtot_resoDown;
  Float_t         mTtot_scaleUp;
  Float_t         mTtot_scaleDown;

  //  Float_t         m_sv;
  //  Float_t         pt_sv;
  //  Float_t         eta_sv;
  //  Float_t         phi_sv;

  Float_t         pt_1;
  Float_t         phi_1;
  Float_t         eta_1;
  Float_t         m_1;
  Int_t           q_1;
  Float_t         iso_1;
  Float_t         mva_1;
  Float_t         d0_1;
  Float_t         dZ_1;

  Float_t         mt_1;
  Float_t         mt_puppi_1;
  Float_t         mt_mva_1;

  Float_t         mt_mva_uncorr_1;
  Float_t         mt_mva_scaleUp_1;
  Float_t         mt_mva_scaleDown_1;
  Float_t         mt_mva_resoUp_1;
  Float_t         mt_mva_resoDown_1;

  Int_t           gen_match_1; 

  Float_t         pt_2;
  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         m_2;
  Int_t           q_2;
  Float_t         iso_2;
  Float_t         d0_2;
  Float_t         dZ_2;
  Float_t         mva_2;
  Float_t         mt_2;
  Float_t         mt_puppi_2;
  Float_t         mt_mva_2;
  Int_t           decayMode_2;
  Int_t           gen_match_2;

  Float_t         genmet;
  Float_t         genmetphi;

  Bool_t          os;
  Bool_t          dilepton_veto;
  Bool_t          extraelec_veto;
  Bool_t          extramuon_veto;

  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t         againstElectronLooseMVA5_1;
  Float_t         againstElectronMediumMVA5_1;
  Float_t         againstElectronTightMVA5_1;
  Float_t         againstElectronVLooseMVA5_1;
  Float_t         againstElectronVTightMVA5_1;
  Float_t         againstMuonLoose3_1;
  Float_t         againstMuonTight3_1;

  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         againstElectronLooseMVA5_2;
  Float_t         againstElectronMediumMVA5_2;
  Float_t         againstElectronTightMVA5_2;
  Float_t         againstElectronVLooseMVA5_2;
  Float_t         againstElectronVTightMVA5_2;
  Float_t         againstMuonLoose3_2;
  Float_t         againstMuonTight3_2;
  Float_t         againstElectronVLooseMVA6_2;
  Float_t         byIsolationMVArun2v1DBoldDMwLTraw_2;
  Float_t         byTightIsolationMVArun2v1DBoldDMwLT_2;

  Float_t         met;
  Float_t         metphi;
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;

  Float_t         puppimet;
  Float_t         puppimetphi;
  //  Float_t         m_sv_puppi;
  //  Float_t         pt_sv_puppi;
  //  Float_t         eta_sv_puppi;
  //  Float_t         phi_sv_puppi;

  Float_t         mvamet;
  Float_t         mvametphi;

  Float_t         mvamet_uncorr;
  Float_t         mvametphi_uncorr;

  Float_t         mvamet_scaleUp;
  Float_t         mvametphi_scaleUp;

  Float_t         mvamet_scaleDown;
  Float_t         mvametphi_scaleDown;

  Float_t         mvamet_resoUp;
  Float_t         mvametphi_resoUp;

  Float_t         mvamet_resoDown;
  Float_t         mvametphi_resoDown;

  Float_t         mvacov00;
  Float_t         mvacov01;
  Float_t         mvacov10;
  Float_t         mvacov11;

  Float_t         m_sv_mvamet;
  Float_t         pt_sv_mvamet;
  Float_t         eta_sv_mvamet;
  Float_t         phi_sv_mvamet;

  Float_t         m_sv_mvamet_uncorr;
  Float_t         m_sv_mvamet_resoUp;
  Float_t         m_sv_mvamet_resoDown;
  Float_t         m_sv_mvamet_scaleUp;
  Float_t         m_sv_mvamet_scaleDown;
  Float_t         m_sv_mvamet_tauUp;
  Float_t         m_sv_mvamet_tauDown;

  Float_t         pt_tt;
  Float_t         pzetavis;
  Float_t         pzetamiss;
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
  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jptraw_2;
  Float_t         jptunc_2;
  Float_t         jmva_2;
  Float_t         jlrm_2;
  Int_t           jctm_2;
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
  Float_t         bosonPt;

  Float_t   Hresponse;
  Float_t   Hparal;
  Float_t   Hperp;

  Float_t   Hresponse_Fast;
  Float_t   Hparal_Fast;
  Float_t   Hperp_Fast;

  Float_t Hresponse_Reco;
  Float_t Hparal_Reco;
  Float_t Hperp_Reco;
  

  Bool_t isZLL;
  Bool_t isZMM;
  Bool_t isZEE;
  Bool_t isZTT;
  Bool_t isZL;
  Bool_t isZJ;

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
  tree->Branch("isZL",&isZL,"isZL/O");
  tree->Branch("isZJ",&isZJ,"isZL/O");

  tree->Branch("mcweight", &mcweight, "mcweight/F");
  tree->Branch("puweight", &puweight, "puweight/F");
  tree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
  tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  tree->Branch("idweight_1", &idweight_1, "idweight_1/F");
  tree->Branch("idweight_2", &idweight_2, "idweight_2/F");
  tree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
  tree->Branch("isoweight_2", &isoweight_2, "isoweight_2/F");
  tree->Branch("effweight", &effweight, "effweight/F");
  tree->Branch("fakeweight", &fakeweight, "fakeweight/F");
  tree->Branch("embeddedWeight", &embeddedWeight, "embeddedWeight/F");
  tree->Branch("signalWeight", &signalWeight, "signalWeight/F");
  tree->Branch("weight", &weight, "weight/F");

  tree->Branch("m_vis",         &m_vis,         "m_vis/F");
  tree->Branch("m_vis_tauUp",   &m_vis_tauUp,   "m_vis_tauUp/F");
  tree->Branch("m_vis_tauDown", &m_vis_tauDown, "m_vis_tauDown/F");
  tree->Branch("m_vis_resoUp",   &m_vis_resoUp,   "m_vis_resoUp/F");
  tree->Branch("m_vis_resoDown", &m_vis_resoDown, "m_vis_resoDown/F");
  tree->Branch("m_vis_scaleUp",   &m_vis_scaleUp,   "m_vis_scaleUp/F");
  tree->Branch("m_vis_scaleDown", &m_vis_scaleDown, "m_vis_scaleDown/F");

  tree->Branch("mTtot",         &mTtot,         "mTtot/F");
  tree->Branch("mTtot_tauUp",   &mTtot_tauUp,   "mTtot_tauUp/F");
  tree->Branch("mTtot_tauDown", &mTtot_tauDown, "mTtot_tauDown/F");
  tree->Branch("mTtot_resoUp",   &mTtot_resoUp,   "mTtot_resoUp/F");
  tree->Branch("mTtot_resoDown", &mTtot_resoDown, "mTtot_resoDown/F");
  tree->Branch("mTtot_scaleUp",   &mTtot_scaleUp,   "mTtot_scaleUp/F");
  tree->Branch("mTtot_scaleDown", &mTtot_scaleDown, "mTtot_scaleDown/F");

  tree->Branch("pt_1", &pt_1, "pt_1/F");
  tree->Branch("phi_1", &phi_1, "phi_1/F");
  tree->Branch("eta_1", &eta_1, "eta_1/F");
  tree->Branch("m_1", &m_1, "m_1/F");
  tree->Branch("q_1", &q_1, "q_1/I");
  tree->Branch("iso_1", &iso_1, "iso_1/F");
  tree->Branch("mva_1", &mva_1, "mva_1/F");
  tree->Branch("d0_1", &d0_1, "d0_1/F");
  tree->Branch("dZ_1", &dZ_1, "dZ_1/F");
  tree->Branch("mt_1", &mt_1, "mt_1/F");
  tree->Branch("mt_puppi_1",&mt_puppi_1,"mt_puppi_1/F");
  tree->Branch("mt_mva_1",&mt_mva_1,"mt_mva_1/F");

  tree->Branch("mt_mva_uncorr_1",&mt_mva_uncorr_1,"mt_mva_uncorr_1/F");
  tree->Branch("mt_mva_scaleUp_1",&mt_mva_scaleUp_1,"mt_mva_scaleUp_1/F");
  tree->Branch("mt_mva_scaleDown_1",&mt_mva_scaleDown_1,"mt_mva_scaleDown_1/F");
  tree->Branch("mt_mva_resoUp_1",&mt_mva_resoUp_1,"mt_mva_resoUp_1/F");
  tree->Branch("mt_mva_resoDown_1",&mt_mva_resoDown_1,"mt_mva_resoDown_1/F");

  tree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");

  tree->Branch("pt_2", &pt_2, "pt_2/F");
  tree->Branch("phi_2", &phi_2, "phi_2/F");
  tree->Branch("eta_2", &eta_2, "eta_2/F");
  tree->Branch("m_2", &m_2, "m_2/F");
  tree->Branch("q_2", &q_2, "q_2/I");
  tree->Branch("iso_2", &iso_2, "iso_2/F");
  tree->Branch("d0_2", &d0_2, "d0_2/F");
  tree->Branch("dZ_2", &dZ_2, "dZ_2/F");
  tree->Branch("mva_2", &mva_2, "mva_2/F");
  tree->Branch("mt_2", &mt_2, "mt_2/F");
  tree->Branch("mt_puppi_2",&mt_puppi_2,"mt_puppi_2/F");
  tree->Branch("mt_mva_2",&mt_mva_2,"mt_mva_2/F");
  tree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");
  tree->Branch("decayMode_2",&decayMode_2,"decayMode_2/I");

  tree->Branch("os", &os, "os/O");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/O");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/O");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/O");

  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
  tree->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2",&byLooseCombinedIsolationDeltaBetaCorr3Hits_2,"byLooseCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2",&byMediumCombinedIsolationDeltaBetaCorr3Hits_2,"byMediumCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byIsolationMVArun2v1DBoldDMwLTraw_2",&byIsolationMVArun2v1DBoldDMwLTraw_2,"byIsolationMVArun2v1DBoldDMwLTraw_2/F");
  tree->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2",&byTightCombinedIsolationDeltaBetaCorr3Hits_2,"byTightCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2",&byTightIsolationMVArun2v1DBoldDMwLT_2,"byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2, "againstElectronLooseMVA5_2/F");
  tree->Branch("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2, "againstElectronMediumMVA5_2/F");
  tree->Branch("againstElectronTightMVA5_2", &againstElectronTightMVA5_2, "againstElectronTightMVA5_2/F");
  tree->Branch("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2, "againstElectronVLooseMVA5_2/F");
  tree->Branch("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2, "againstElectronVTightMVA5_2/F");
  tree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/F");
  tree->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");
  tree->Branch("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2,"againstElectronVLooseMVA6_2/F");

  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("metcov00", &metcov00, "metcov00/F");
  tree->Branch("metcov01", &metcov01, "metcov01/F");
  tree->Branch("metcov10", &metcov10, "metcov10/F");
  tree->Branch("metcov11", &metcov11, "metcov11/F");

  tree->Branch("puppimet", &puppimet, "puppimet/F");
  tree->Branch("puppimetphi", &puppimetphi, "puppimetphi/F");

  tree->Branch("mvamet", &mvamet, "mvamet/F");
  tree->Branch("mvametphi", &mvametphi, "mvametphi/F");
  tree->Branch("mvamet_uncorr", &mvamet_uncorr, "mvamet_uncorr/F");
  tree->Branch("mvametphi_uncorr", &mvametphi_uncorr, "mvametphi_uncorr/F");

  tree->Branch("mvamet_scaleUp", &mvamet_scaleUp, "mvamet_scaleUp/F");
  tree->Branch("mvametphi_scaleUp", &mvametphi_scaleUp, "mvametphi_scaleUp/F");
  
  tree->Branch("mvamet_scaleDown", &mvamet_scaleDown, "mvamet_scaleDown/F");
  tree->Branch("mvametphi_scaleDown", &mvametphi_scaleDown, "mvametphi_scaleDown/F");

  tree->Branch("mvamet_resoUp", &mvamet_resoUp, "mvamet_resoUp/F");
  tree->Branch("mvametphi_resoUp", &mvametphi_resoUp, "mvametphi_resoUp/F");
  
  tree->Branch("mvamet_resoDown", &mvamet_resoDown, "mvamet_resoDown/F");
  tree->Branch("mvametphi_resoDown", &mvametphi_resoDown, "mvametphi_resoDown/F");

  tree->Branch("mvacov00", &mvacov00, "mvacov00/F");
  tree->Branch("mvacov01", &mvacov01, "mvacov01/F");
  tree->Branch("mvacov10", &mvacov10, "mvacov10/F");
  tree->Branch("mvacov11", &mvacov11, "mvacov11/F");

  tree->Branch("genmet",&genmet,"genmet/F");
  tree->Branch("genmetphi",&genmetphi,"genmetphi/F");

  //  tree->Branch("m_sv", &m_sv, "m_sv/F");
  //  tree->Branch("pt_sv", &pt_sv, "pt_sv/F");
  //  tree->Branch("eta_sv", &eta_sv, "eta_sv/F");
  //  tree->Branch("phi_sv", &phi_sv, "phi_sv/F");

  //  tree->Branch("m_sv_puppi", &m_sv_puppi, "m_sv_puppi/F");
  //  tree->Branch("pt_sv_puppi", &pt_sv_puppi, "pt_sv_puppi/F");
  //  tree->Branch("eta_sv_puppi", &eta_sv_puppi, "eta_sv_puppi/F");
  //  tree->Branch("phi_sv_puppi", &phi_sv_puppi, "phi_sv_puppi/F");

  tree->Branch("m_sv_mvamet",   &m_sv_mvamet,   "m_sv_mvamet/F");
  tree->Branch("pt_sv_mvamet",  &pt_sv_mvamet,  "pt_sv_mvamet/F");
  tree->Branch("eta_sv_mvamet", &eta_sv_mvamet, "eta_sv_mvamet/F");
  tree->Branch("phi_sv_mvamet", &phi_sv_mvamet, "phi_sv_mvamet/F");

  tree->Branch("m_sv_mvamet_uncorr",     &m_sv_mvamet_uncorr,    "m_sv_mvamet_uncorr/F");
  tree->Branch("m_sv_mvamet_tauUp",      &m_sv_mvamet_tauUp,     "m_sv_mvamet_tauUp/F");
  tree->Branch("m_sv_mvamet_tauDown",    &m_sv_mvamet_tauDown,   "m_sv_mvamet_tauDown/F");
  tree->Branch("m_sv_mvamet_scaleUp",    &m_sv_mvamet_scaleUp,   "m_sv_mvamet_scaleUp/F");
  tree->Branch("m_sv_mvamet_scaleDown",  &m_sv_mvamet_scaleDown, "m_sv_mvamet_scaleDown/F");
  tree->Branch("m_sv_mvamet_resoUp",     &m_sv_mvamet_resoUp,    "m_sv_mvamet_resoUp/F");
  tree->Branch("m_sv_mvamet_resoDown",   &m_sv_mvamet_resoDown,  "m_sv_mvamet_resoDown/F");

  tree->Branch("pt_tt", &pt_tt, "pt_tt/F");
  tree->Branch("pzetavis", &pzetavis, "pzetavis/F");
  tree->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
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
  tree->Branch("bosonPt",&bosonPt,"bosonPt/F");

  //  tree->Branch("Hresponse",&Hresponse,"Hresponse/F");
  //  tree->Branch("Hparal",&Hparal,"Hparal/F");
  //  tree->Branch("Hperp",&Hperp,"Hperp/F");

  //  tree->Branch("Hresponse_Fast",&Hresponse_Fast,"Hresponse_Fast/F");
  //  tree->Branch("Hparal_Fast",&Hparal_Fast,"Hparal_Fast/F");
  //  tree->Branch("Hperp_Fast",&Hperp_Fast,"Hperp_Fast/F");

  //  tree->Branch("Hresponse_Reco",&Hresponse_Reco,"Hresponse_Reco/F");
  //  tree->Branch("Hparal_Reco",&Hparal_Reco,"Hparal_Reco/F");
  //  tree->Branch("Hperp_Reco",&Hperp_Reco,"Hperp_Reco/F");

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

  // Official PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Feb02.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Fall15_PU25_V1.root","read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Muon scale factors
  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  ScaleFactor * SF_muonTrig = new ScaleFactor();
  SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerFile));

  // Recoil corrector and met systematics
  RecoilCorrector recoilMvaMetCorrector(RecoilMvaFileName);
  MEtSys metSys(MetSysFileName);

  // SV fit mass
  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());

    // BTag scale factors
  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2.csv");
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,  // operating point
			       "central");            // systematics type
  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_B,    // btag flavour
	      "mujets");            // measurement type
  
  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_C,    // btag flavour
	      "mujets");             // measurement type

  reader.load(calib,               // calibration instance
	      BTagEntry::FLAV_UDSG,    // btag flavour
	      "incl");             // measurement type
  
  //BTagCalibrationReader reader_BC(&calib,BTagEntry::OP_MEDIUM,"mujets","central");           // systematics type
  //BTagCalibrationReader reader_Light(&calib,BTagEntry::OP_MEDIUM,"incl","central");           // systematics type

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


  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  // Looping over RooT files 
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));

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
    
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
    int NE = int(histoInputEvents->GetEntries());
    std::cout << "      number of input events    = " << NE << std::endl;
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;
    AC1B analysisTree(_tree);
    Long64_t numberOfEntries = analysisTree.GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 

      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

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

      isZLL = false;
      isZEE = false;
      isZMM = false;
      isZTT = false;
      bool isZfound = false;
      bool isWfound = false;
      bool isHfound = false;

      nuPx = 0;
      nuPy = 0;
      nuPz = 0;

      lepPx = 0;
      lepPy = 0;
      lepPz = 0;

      bosonPx = 0;
      bosonPy = 0;
      bosonPz = 0;

      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      std::vector<TLorentzVector> nonpromptNeutrinos; nonpromptNeutrinos.clear();
      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0,0,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0,0,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0,0,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector nonpromptNeutrinosLV;  nonpromptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);
      if (!isData) {
	for (unsigned int igentau=0; igentau<analysisTree.gentau_count; ++igentau) {
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
	  // bool isLepton =
	  //   abs(analysisTree.genparticles_pdgid[igen])==11 ||
	  //   abs(analysisTree.genparticles_pdgid[igen])==13;
	  // bool isHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]>0 && analysisTree.genparticles_status[igen]>0;
	  // bool isNeutrino = 
	  //   abs(analysisTree.genparticles_pdgid[igen])==12 ||
	  //   abs(analysisTree.genparticles_pdgid[igen])==14 ||
	  //   abs(analysisTree.genparticles_pdgid[igen])==16;
	  // bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]>0;
	  // bool addFullV = ((isLepton||isNeutrino) &&isHardProcessFinalState) || isDirectHardProcessTauDecayProduct;
	  // bool addVisV  = (isLepton && isHardProcessFinalState) || (isDirectHardProcessTauDecayProduct &&!isNeutrino);

	  // if (addFullV) fullVLV += genLV;
	  // if (addVisV) visVLV += genLV;

	  if (analysisTree.genparticles_pdgid[igen]==23) { 
	    isZfound = true;
	    zBosonLV = genLV;
	  }
	  if (abs(analysisTree.genparticles_pdgid[igen])==24) { 
	    isWfound = true;
	    wBosonLV = genLV;
	  }
	  if (abs(analysisTree.genparticles_pdgid[igen])==25||
	      abs(analysisTree.genparticles_pdgid[igen])==35||
	      abs(analysisTree.genparticles_pdgid[igen])==36) {
            isHfound = true;
            hBosonLV = genLV;
          }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	    if (analysisTree.genparticles_isPrompt[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptElectrons.push_back(genLV);
	      promptElectronsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	    if (analysisTree.genparticles_isPrompt[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptMuons.push_back(genLV);
	      promptMuonsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==12||
	      fabs(analysisTree.genparticles_pdgid[igen])==14||
	      fabs(analysisTree.genparticles_pdgid[igen])==16)  {
	    if (analysisTree.genparticles_isPrompt[igen]&&
		!analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      promptNeutrinos.push_back(genLV);
	      promptNeutrinosLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	    if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      nonpromptNeutrinos.push_back(genLV);
	      nonpromptNeutrinosLV += genLV;
	    }
	  }
	}
      }
      /*      std::cout << "Taus (first copy) : " << promptTausFirstCopy.size() << std::endl;
      std::cout << "Taus (last copy)  : " << promptTausLastCopy.size() << std::endl;
      std::cout << "Muons             : " << promptMuons.size() << std::endl;
      std::cout << "Electrons         : " << promptElectrons.size() << std::endl;
      std::cout << "Neutrinos         : " << promptNeutrinos.size() << std::endl;
      std::cout << "Neutrinos (Taus)  : " << nonpromptNeutrinos.size() << std::endl;
      if (isZfound) 
	std::cout << "ZBoson    px = " << zBosonLV.Px() 
		  << "          py = " << zBosonLV.Py()
		  << "          pz = " << zBosonLV.Pz() 
		  << "          mass = " << zBosonLV.M() << std::endl;
      if (isWfound) {
	std::cout << "WBoson    px = " << wBosonLV.Px() 
		  << "          py = " << wBosonLV.Py()
		  << "          pz = " << wBosonLV.Pz() 
		  << "          mass = " << wBosonLV.M() << std::endl;
	std::cout << "W(Dec)    px = " << wDecayProductsLV.Px() 
		  << "          py = " << wDecayProductsLV.Py()
		  << "          pz = " << wDecayProductsLV.Pz() 
		  << "          mass = " << wDecayProductsLV.M() << std::endl;
      }
      if (isHfound) 
	std::cout << "HBoson    px = " << hBosonLV.Px()
		  << "          py = " << hBosonLV.Py()
		  << "          pz = " << hBosonLV.Pz() 
		  << "          mass = " << hBosonLV.M() << std::endl;
      std::cout << "Taus       px = " << promptTausLV.Px() 
      		<< "           py = " << promptTausLV.Py()
      		<< "           pz = " << promptTausLV.Pz() << std::endl;
      std::cout << "VisTaus    px = " << promptVisTausLV.Px() 
      		<< "           py = " << promptVisTausLV.Py()
      		<< "           pz = " << promptVisTausLV.Pz() << std::endl;
      std::cout << "Muons      px = " << promptMuonsLV.Px() 
      		<< "           py = " << promptMuonsLV.Py()
      		<< "           pz = " << promptMuonsLV.Pz() << std::endl;
      std::cout << "Elect.     px = " << promptElectronsLV.Px() 
      		<< "           py = " << promptElectronsLV.Py()
      		<< "           pz = " << promptElectronsLV.Pz() << std::endl;
      std::cout << "Neut.      pz = " << promptNeutrinosLV.Px()
                << "           py = " << promptNeutrinosLV.Py()
                << "           pz = " << promptNeutrinosLV.Pz() << std::endl;
      std::cout << "Nu(tau)    px = " << nonpromptNeutrinosLV.Px()
                << "           py = " << nonpromptNeutrinosLV.Py()
                << "           pz = " << nonpromptNeutrinosLV.Pz() << std::endl;
      TLorentzVector totTausLV = promptVisTausLV+nonpromptNeutrinosLV;
      std::cout << "Taus(test) px = " << totTausLV.Px()
		<< "           py = " << totTausLV.Py()
		<< "           pz = " << totTausLV.Pz() << std::endl;
      std::cout << std::endl; */

      if (isDY) {
	if (promptTausFirstCopy.size()==2) {
	  bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
	  lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	}
	else if (promptMuons.size()==2) {
	  bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
	  lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	}
	else {
	  bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
	  lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	}
      }
      else if (isW) {
	if (isWfound) {
	  bosonPx = wBosonLV.Px(); bosonPy = wBosonLV.Py(); bosonPz = wBosonLV.Pz();
	}
	else {
	  bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
	}
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
      bosonPt = TMath::Sqrt(bosonPx*bosonPx+bosonPy*bosonPy);
      isZLL = isZEE || isZMM;
      ALL->Fill(0.);
	
      if (isZEE) ZEE->Fill(0.);
      if (isZMM) ZMM->Fill(0.);
      if (isZTT) ZTT->Fill(0.);
      if (!isZfound) GAM->Fill(0.);
      
      //      cout << "Ok2" << endl;

      run = int(analysisTree.event_run);
      lumi = int(analysisTree.event_luminosityblock);
      evt = int(analysisTree.event_nr);

      // weights
      mcweight = analysisTree.genweight;
      puweight = 1;
      trigweight_1 = 1;
      trigweight_2 = 1;
      idweight_1 = 1;
      idweight_2 = 1;
      isoweight_1 = 1;
      isoweight_2 = 1;
      effweight = 1;
      fakeweight = 1;
      embeddedWeight = 1;
      signalWeight = 1;
      weight = 1;
      
      npv = analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
      rho = analysisTree.rho;
      if (!isData)
	puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      
      unsigned int nIsoMuonLeg = 0;
      bool isIsoMuonLeg = false;

      unsigned int nMuonTauMuonLeg = 0;
      bool isMuonTauMuonLeg = false;
      
      unsigned int nMuonTauOverlap = 0;
      bool isMuonTauOverlap = false;

      unsigned int nMuonTauTauLeg = 0;
      bool isMuonTauTauLeg = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //      std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==IsoMuonLeg) {
	  nIsoMuonLeg = i;
	  isIsoMuonLeg = true;
	}
	if (HLTFilter==MuonTauMuonLeg) {
	  nMuonTauMuonLeg = i;
	  isMuonTauMuonLeg = true;
	}
	if (HLTFilter==MuonTauOverlap) {
	  nMuonTauOverlap = i;
	  isMuonTauOverlap = true;
	}
	if (HLTFilter==MuonTauTauLeg) {
	  nMuonTauTauLeg = i;
	  isMuonTauTauLeg = true;
	}
      }
      if (!isIsoMuonLeg) {
	std::cout << "HLT filter " << IsoMuonLeg << " not found" << std::endl;
	exit(-1);
      }
      if (!isMuonTauMuonLeg) {
	std::cout << "HLT filter " << MuonTauMuonLeg << " not found" << std::endl;
	exit(-1);
      }
      if (!isMuonTauOverlap) {
	std::cout << "HLT filter " << MuonTauOverlap << " not found" << std::endl;
	exit(-1);
      }
      if (!isMuonTauTauLeg) {
	std::cout << "HLT filter " << MuonTauTauLeg << " not found" << std::endl;
	exit(-1);
      }
      unsigned int nBTagDiscriminant = 0;
      for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
	TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
	if (discr.Contains(BTagDiscriminator))
	  nBTagDiscriminant = iBTag;
      }

      // tau selection
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	if (analysisTree.tau_decayModeFinding[it]<=0.5) continue;
	if (analysisTree.tau_pt[it]<ptTauCut) continue;
	if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
	if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
	if (fabs(analysisTree.tau_charge[it])<0.5||
	    fabs(analysisTree.tau_charge[it])>1.5) continue;
	taus.push_back(it);
      }

      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	muons.push_back(im);
      }

      //      cout << "Ok3" << endl;

      if (taus.size()==0) continue;
      if (muons.size()==0) continue;
      
      // selecting muon and electron pair (OS or SS);
      int tauIndex = -1;
      int muonIndex = -1;

      float isoMuMin = 1e+10;
      float isoTauMin = -10;
      float ptMu = 0;
      float ptTau = 0;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;
      for (unsigned int im=0; im<muons.size(); ++im) {
	bool isIsoMuonLegMatch = false;
	//	bool isMuonTauMuonLegMatch = false;
	//	bool isMuonTauOverlapMuonMatch = false;
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
	  if (analysisTree.trigobject_filters[iT][nIsoMuonLeg]
	      &&analysisTree.muon_pt[mIndex]>ptMuonHighCut&&
	      analysisTree.trigobject_pt[iT]>ptMuonTriggerCut) { // IsoMu Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isIsoMuonLegMatch = true;
	    }
	  }
	  //	  if (analysisTree.trigobject_filters[iT][nMuonTauMuonLeg]) { // MuonTau Muon Leg
	  //	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
	  //				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  //	    if (dRtrig<deltaRTrigMatch) {
	  //	      isMuonTauMuonLegMatch = true;
	  //	    }
	  //	  }
	  //	  if (analysisTree.trigobject_filters[iT][nMuonTauOverlap]&&analysisTree.trigobject_isMuon[iT]) { // MuonTau Overlap Muon 
	  //	    float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
	  //				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  //	    if (dRtrig<deltaRTrigMatch) {
	  //	      isMuonTauOverlapMuonMatch = true;
	  //	    }
	  //	  }
	}

	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;

	  //	  bool isMuonTauOverlapTauMatch = false;
	  //	  bool isMuonTauTauLegMatch = false;

	  //	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  //	    if (analysisTree.trigobject_filters[iT][nMuonTauOverlap]&&analysisTree.trigobject_isTau[iT]) { // MuonTau Overlap Tau
	  //	      float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
	  //				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  //	      if (dRtrig<deltaRTrigMatch) {
	  //		isMuonTauOverlapTauMatch = true;
	  //	      }
	  //	    }
	  //	    if (analysisTree.trigobject_filters[iT][nMuonTauTauLeg]) { // MuonTau Tau Leg
	  //	      float dRtrig = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
	  //				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  //	      if (dRtrig<deltaRTrigMatch) {
	  //		isMuonTauTauLegMatch = true;
	  //	      }
	  //	    }
	  //	  }

	  bool trigMatch = 
	    (isIsoMuonLegMatch);
	  //	  || (isMuonTauOverlapMuonMatch&&isMuonTauMuonLegMatch&&isMuonTauOverlapTauMatch&&isMuonTauTauLegMatch);
	  //	  std::cout << "Trigger match = " << trigMatch << std::endl;

	  if (applyTriggerMatch && !trigMatch) continue;
	  
	  //	  bool isKinematicMatch = false;
	  //	  if (isMu23&&isEle12) {
	  //	    if (analysisTree.muon_pt[mIndex]>ptMuonHighCut&&analysisTree.electron_pt[eIndex]>ptElectronLowCut)
	  //	      isKinematicMatch = true;
	  //	  }
	  //	  if (isMu8&&isEle23) {
	  //            if (analysisTree.muon_pt[mIndex]>ptMuonLowCut&&analysisTree.electron_pt[eIndex]>ptElectronHighCut)
	  //              isKinematicMatch = true;
	  //          }
	  //	  if (!isKinematicMatch) continue;

	  //	  float isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex];
	  float isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];

	  if (int(mIndex)!=muonIndex) {
	    if (relIsoMu==isoMuMin) {
	      if (analysisTree.muon_pt[mIndex]>ptMu) {
		isoMuMin  = relIsoMu;
		ptMu = analysisTree.muon_pt[mIndex];
		muonIndex = int(mIndex);
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tauIndex = int(tIndex);
	      }
	    }
	    else if (relIsoMu<isoMuMin) {
	      isoMuMin  = relIsoMu;
	      ptMu = analysisTree.muon_pt[mIndex];
	      muonIndex = int(mIndex);
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tauIndex = int(tIndex);
	    }
	  }
	  else {
	    if (isoTau==isoTauMin) {
	      if (analysisTree.tau_pt[tIndex]>ptTau) {
		ptTau = analysisTree.tau_pt[tIndex];
		isoTauMin = isoTau;
		tauIndex = int(tIndex);
	      }
	    }
	    else if (isoTau>isoTauMin) {
	      ptTau = analysisTree.tau_pt[tIndex];
	      isoTauMin = isoTau;
	      tauIndex = int(tIndex);
	    }
	  }
	  
	}
      }

      //      std::cout << "mIndex = " << muonIndex << "   tauIndex = " << tauIndex << std::endl;

      if (tauIndex<0) continue;
      if (muonIndex<0) continue;
      //      std::cout << "Ok4 " << std::endl;

      os = (analysisTree.muon_charge[muonIndex]*analysisTree.tau_charge[tauIndex]) < 0;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	//	bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[ie],
	//					       analysisTree.electron_superclusterEta[ie],
	//					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
	bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
	if (!electronMvaId&&applyVetoElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
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

      // looking for extra muon's (dimuon veto)
      bool foundExtraMuon = false;
      vector<int> mu_dimuons; mu_dimuons.clear(); 
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

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

	if (analysisTree.muon_pt[im]>ptDilepMuonCut&&
	    fabs(analysisTree.muon_eta[im])<etaDilepMuonCut&&
	    analysisTree.muon_isGlobal[im]&&
	    analysisTree.muon_isTracker[im]&&
	    analysisTree.muon_isPF[im]&&
	    fabs(analysisTree.muon_dxy[im])<dxyDilepMuonCut&&
	    fabs(analysisTree.muon_dz[im])<dzDilepMuonCut&&
	    relIsoMu<isoDilepMuonCut)
	  mu_dimuons.push_back(im);

	if (int(im)==muonIndex) continue;
	if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

      dilepton_veto = false;
      if (mu_dimuons.size()>1) {
	for (unsigned int i1=0; i1<mu_dimuons.size()-1; ++i1) {
	  unsigned int indx1 = mu_dimuons[i1];
	  for (unsigned int i2=i1+1; i2<mu_dimuons.size(); ++i2 ) {
	    unsigned int indx2 = mu_dimuons[i2];
	    float dRmuons = deltaR(analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1],
				   analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	    if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[indx1]*analysisTree.muon_charge[indx2]<0)) dilepton_veto = true;
 	  }
	}
      }
      //      cout << analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex] << endl;
      //      cout << analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex] << endl;

      // applying inclusive selection
      if (applyInclusiveSelection) {

	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
	if (dilepton_veto)  continue;

	bool tauPass = 
	  analysisTree.tau_againstElectronVLooseMVA6[tauIndex]>0.5 &&
	  // analysisTree.tau_againstElectronVLooseMVA5[tauIndex]>0.5
	  analysisTree.tau_againstMuonTight3[tauIndex]>0.5 &&
	  //	  analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex]<isoTauCut;
	  //	  analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex] > 0.5;
	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex] > 0.5;
	if (!tauPass) continue;
	
	bool muonPass = isoMuMin<isoMuonHighCut;
	if (!muonPass) continue;
	//	if (!os) continue;
	
      }
      
      //      cout << "Ok5" << endl;

      // met
      float pfmet_ex = analysisTree.pfmet_ex;
      float pfmet_ey = analysisTree.pfmet_ey;
      met = TMath::Sqrt(pfmet_ex*pfmet_ex + pfmet_ey*pfmet_ey);
      metphi = TMath::ATan2(pfmet_ey,pfmet_ex);

      // puppimet
      float puppimet_ex = analysisTree.puppimet_ex;
      float puppimet_ey = analysisTree.puppimet_ey;
      puppimet = TMath::Sqrt(puppimet_ex*puppimet_ex + puppimet_ey*puppimet_ey);
      puppimetphi = TMath::ATan2(puppimet_ey,puppimet_ex);

      // genmet
      float genmet_x = analysisTree.genmet_ex;
      float genmet_y = analysisTree.genmet_ey;
      genmet = TMath::Sqrt(genmet_x*genmet_x+genmet_y*genmet_y);
      genmetphi = TMath::ATan2(genmet_y,genmet_x);

      metcov00 = analysisTree.pfmet_sigxx;
      metcov01 = analysisTree.pfmet_sigxy;
      metcov10 = analysisTree.pfmet_sigyx;
      metcov11 = analysisTree.pfmet_sigyy;

      float mvamet_ex = 0;
      float mvamet_ey = 0;
      mvamet = 0;
      mvametphi = 0;
      mvacov00 = 0;
      mvacov01 = 0;
      mvacov10 = 0;
      mvacov11 = 0;

      unsigned int metIndex = 0;
      bool mvaFound = false;
      if (analysisTree.mvamet_count>0) {
       	for (unsigned int imva=0; imva<analysisTree.mvamet_count; ++imva) {
       	  if (analysisTree.mvamet_channel[imva]==3) {
       	    if (int(analysisTree.mvamet_lep1[imva])==tauIndex&&int(analysisTree.mvamet_lep2[imva])==muonIndex) {
       	      metIndex = imva;
       	      mvaFound = true;
       	    }
       	  }
       	}
      }
      if (!mvaFound) {
	cout << "Warning... mva Met is not found" << endl;
      }

      if (analysisTree.mvamet_count>0) {
      	mvamet_ex = analysisTree.mvamet_ex[metIndex];
      	mvamet_ey = analysisTree.mvamet_ey[metIndex];
      	float mvamet_x2 = mvamet_ex * mvamet_ex;
      	float mvamet_y2 = mvamet_ey * mvamet_ey;
	
      	mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
      	mvametphi = TMath::ATan2(mvamet_ey,mvamet_ex);
      	mvacov00 = analysisTree.mvamet_sigxx[metIndex];
      	mvacov01 = analysisTree.mvamet_sigxy[metIndex];
      	mvacov10 = analysisTree.mvamet_sigyx[metIndex];
      	mvacov11 = analysisTree.mvamet_sigyy[metIndex];
      }

      //      std::cout << "Ok6" << std::endl;

      // filling muon variables
      pt_1 = analysisTree.muon_pt[muonIndex];
      eta_1 = analysisTree.muon_eta[muonIndex];
      phi_1 = analysisTree.muon_phi[muonIndex];
      q_1 = -1;
      if (analysisTree.muon_charge[muonIndex]>0)
	q_1 = 1;
      mva_1 = -9999;
      d0_1 = analysisTree.muon_dxy[muonIndex];
      dZ_1 = analysisTree.muon_dz[muonIndex];
      iso_1 = isoMuMin;
      m_1 = muonMass;

      // filling tau variables
      pt_2 = analysisTree.tau_pt[tauIndex];
      eta_2 = analysisTree.tau_eta[tauIndex];
      phi_2 = analysisTree.tau_phi[tauIndex];
      q_2 = -1;
      if (analysisTree.tau_charge[tauIndex]>0)
        q_2 = 1;
      mva_2 = -9999;
      d0_2 = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex];
      dZ_2 = analysisTree.tau_leadchargedhadrcand_dz[tauIndex];
      iso_2 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];
      m_2 = analysisTree.tau_mass[tauIndex];

      // counting jets
      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetshad; jetshad.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> hadjets; hadjets.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;
      float ptLeadingBJet = -1;

      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	float jetEta = analysisTree.pfjet_eta[jet];

	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta>jetEtaCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];
	if (jetPt<jetPtLowCut) continue;

	// jetId
	bool isPFJetId = looseJetiD(analysisTree,int(jet));
	if (!isPFJetId) continue;

	float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   eta_1,phi_1);
	if (dR1<dRJetLeptonCut) continue;

	float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           eta_2,phi_2);
        if (dR2<dRJetLeptonCut) continue;

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
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
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
      float mvamet_uncorr_ex =  mvamet_ex;
      float mvamet_uncorr_ey =  mvamet_ey;
      float mvamet_uncorr    =  mvamet;
      float mvametphi_uncorr =  mvametphi;

      float mvamet_corr_ex =  mvamet_ex;
      float mvamet_corr_ey =  mvamet_ey;
      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;
      if ((isDY||isW)&&!isData)
	recoilMvaMetCorrector.CorrectByMeanResolution(mvamet_ex,mvamet_ey,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,mvamet_corr_ex,mvamet_corr_ey);

      mvamet_ex = mvamet_corr_ex;
      mvamet_ex = mvamet_corr_ex;
      float mvamet_ex2 = mvamet_ex*mvamet_ex;
      float mvamet_ey2 = mvamet_ey*mvamet_ey;
      mvamet = TMath::Sqrt(mvamet_ex2+mvamet_ey2);
      mvametphi = TMath::ATan2(mvamet_ey,mvamet_ex);

      int bkgdType = 0;
      if (isDY||isW)
	bkgdType = MEtSys::ProcessType::BOSON;
      else if (isTOP)
	bkgdType = MEtSys::ProcessType::TOP;
      else 
	bkgdType = MEtSys::ProcessType::EWK; 


      float mvamet_scaleUp_ex = mvamet_ex;
      float mvamet_scaleUp_ey = mvamet_ey;
      float mvamet_scaleDown_ex = mvamet_ex;
      float mvamet_scaleDown_ey = mvamet_ey;
      float mvamet_resoUp_ex = mvamet_ex;
      float mvamet_resoUp_ey = mvamet_ey;
      float mvamet_resoDown_ex = mvamet_ex;
      float mvamet_resoDown_ey = mvamet_ey;

      if (!isData) {
	metSys.ApplyMEtSys(mvamet_ex,mvamet_ey,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Up,
			   mvamet_scaleUp_ex,mvamet_scaleUp_ey);
	metSys.ApplyMEtSys(mvamet_ex,mvamet_ey,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Down,
			   mvamet_scaleDown_ex,mvamet_scaleDown_ey);
	metSys.ApplyMEtSys(mvamet_ex,mvamet_ey,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Up,
			   mvamet_resoUp_ex,mvamet_resoUp_ey);
	metSys.ApplyMEtSys(mvamet_corr_ex,mvamet_corr_ey,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Down,
			   mvamet_resoDown_ex,mvamet_resoDown_ey);
      }
      mvamet_scaleUp = TMath::Sqrt(mvamet_scaleUp_ex*mvamet_scaleUp_ex+
				   mvamet_scaleUp_ey*mvamet_scaleUp_ey);
      mvametphi_scaleUp = TMath::ATan2(mvamet_scaleUp_ey,mvamet_scaleUp_ex);
      
      mvamet_scaleDown = TMath::Sqrt(mvamet_scaleDown_ex*mvamet_scaleDown_ex+
				     mvamet_scaleDown_ey*mvamet_scaleDown_ey);
      mvametphi_scaleDown = TMath::ATan2(mvamet_scaleDown_ey,mvamet_scaleDown_ex);
      
      mvamet_resoUp = TMath::Sqrt(mvamet_resoUp_ex*mvamet_resoUp_ex+
				  mvamet_resoUp_ey*mvamet_resoUp_ey);
      mvametphi_resoUp = TMath::ATan2(mvamet_resoUp_ey,mvamet_resoUp_ex);
      
      mvamet_resoDown = TMath::Sqrt(mvamet_resoDown_ex*mvamet_resoDown_ex+
				    mvamet_resoDown_ey*mvamet_resoDown_ey);
      mvametphi_resoDown = TMath::ATan2(mvamet_resoDown_ey,mvamet_resoDown_ex);


      // printf("Uncorrected    :  %7.2f  %7.2f\n",mvamet_ex,mvamet_ey);
      // printf("Central        :  %7.2f  %7.2f\n",mvamet_corr_ex,mvamet_corr_ey);
      // printf("ScaleUp        :  %7.2f  %7.2f\n",mvamet_scaleUp_ex,mvamet_scaleUp_ey);
      // printf("ScaleDown      :  %7.2f  %7.2f\n",mvamet_scaleDown_ex,mvamet_scaleDown_ey);
      // printf("ResoUp         :  %7.2f  %7.2f\n",mvamet_resoUp_ex,mvamet_resoUp_ey);
      // printf("ResoDown       :  %7.2f  %7.2f\n",mvamet_resoDown_ex,mvamet_resoDown_ey);
      // std::cout << std::endl;


      // met related variables

      float dPhiMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				     pfmet_ex,pfmet_ey);
      mt_1 = TMath::Sqrt(2*met*pt_1*(1-TMath::Cos(dPhiMETMuon)));

      float dPhiPuppiMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
					  puppimet_ex,puppimet_ey);
      mt_puppi_1 = TMath::Sqrt(2*puppimet*pt_1*(1-TMath::Cos(dPhiPuppiMETMuon)));

      float dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
       					mvamet_ex,mvamet_ey);
      mt_mva_1 = TMath::Sqrt(2*mvamet*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));

      dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				  mvamet_uncorr_ex,mvamet_uncorr_ey);
      mt_mva_uncorr_1 = TMath::Sqrt(2*mvamet_uncorr*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));

      
      dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				  mvamet_resoUp_ex,mvamet_resoUp_ey);
      mt_mva_resoUp_1 = TMath::Sqrt(2*mvamet_resoUp*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));      

      dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				  mvamet_resoDown_ex,mvamet_resoDown_ey);
      mt_mva_resoDown_1 = TMath::Sqrt(2*mvamet_resoDown*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));      

      dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				  mvamet_scaleUp_ex,mvamet_scaleUp_ey);
      mt_mva_scaleUp_1 = TMath::Sqrt(2*mvamet_scaleUp*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));      

      dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				  mvamet_scaleDown_ex,mvamet_scaleDown_ey);
      mt_mva_scaleDown_1 = TMath::Sqrt(2*mvamet_scaleDown*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));      



      float dPhiMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.tau_py[tauIndex],
				    pfmet_ex,pfmet_ey);
      mt_2 = TMath::Sqrt(2*met*pt_2*(1-TMath::Cos(dPhiMETTau)));
      float dPhiPuppiMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.tau_py[tauIndex],
					 puppimet_ex,puppimet_ey);
      mt_puppi_2 = TMath::Sqrt(2*puppimet*pt_2*(1-TMath::Cos(dPhiPuppiMETTau)));
      float dPhiMvaMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.muon_py[tauIndex],
       				       mvamet_ex,mvamet_ey);
      mt_mva_2 = TMath::Sqrt(2*mvamet*pt_2*(1-TMath::Cos(dPhiMvaMETTau)));
      mt_mva_2 = 0;
      decayMode_2 = analysisTree.tau_decayMode[tauIndex];

      //      std::cout << "Ok7" << std::endl;

      gen_match_1 = 6;
      gen_match_2 = 6;

      isZTT = false;
      isZL  = false;
      isZJ  = false;

      float minDR_1 = 0.2;
      float minDR_2 = 0.2;
      if (!isData) {
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	  float ptGen = PtoPt(analysisTree.genparticles_px[igen],
			      analysisTree.genparticles_py[igen]);
	  bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
	  bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
	  bool isAnyType = type1 || type2 || type3 || type4;
	  if (isAnyType) {
	    float etaGen = PtoEta(analysisTree.genparticles_px[igen],
				  analysisTree.genparticles_py[igen],
				  analysisTree.genparticles_pz[igen]);
	    float phiGen = PtoPhi(analysisTree.genparticles_px[igen],
				  analysisTree.genparticles_py[igen]);
	    float deltaR_1 = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
				    etaGen,phiGen);
	    if (deltaR_1<minDR_1) {
	      minDR_1 = deltaR_1;
	      if (type1) gen_match_1 = 1;
	      else if (type2) gen_match_1 = 2;
	      else if (type3) gen_match_1 = 3;
	      else if (type4) gen_match_1 = 4;
	    }
	    
	    float deltaR_2 = deltaR(analysisTree.tau_eta[tauIndex],analysisTree.tau_phi[tauIndex],
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
	for (unsigned int igen=0; igen < analysisTree.gentau_count; ++igen) {
	  if (analysisTree.gentau_visibleNoLep_pt[igen]>15.) {
	    float deltaR_1 = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
				    analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
	    if (deltaR_1<minDR_1) {
	      minDR_1 = deltaR_1;
	      gen_match_1 = 5;
	    }
	  float deltaR_2 = deltaR(analysisTree.tau_eta[tauIndex],analysisTree.tau_phi[tauIndex],
				  analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
	  if (deltaR_2<minDR_2) {
	    minDR_2 = deltaR_2;
	    gen_match_2 = 5;
	  }
	  }
	}
	if (gen_match_2==5) isZTT = true;
	else if (gen_match_2==6) isZJ = true;
	else isZL = true;
	isZLL = isZL || isZJ;
      }

      float mvamet_Reco_ex =  mvamet_ex;
      float mvamet_Reco_ey =  mvamet_ey;
      float lep_Reco_px = analysisTree.muon_px[muonIndex];
      float lep_Reco_py = analysisTree.muon_py[muonIndex];
      if (isDY&&isZTT&&!isData) {
	lep_Reco_px += analysisTree.tau_px[tauIndex];
	lep_Reco_py += analysisTree.tau_py[tauIndex];
      }
      
      // muon scale factors
      isoweight_1 = (float)SF_muonIdIso->get_ScaleFactor(double(pt_1),double(eta_1));
      trigweight_1 = (float)SF_muonTrig->get_ScaleFactor(double(pt_1),double(eta_1));
      effweight = isoweight_1*trigweight_1;

      int tau_decayMode = analysisTree.tau_decayMode[tauIndex];

      //      m_sv = -9999;
      //      pt_sv = -9999;
      //      eta_sv = -9999;
      //     phi_sv = -9999;

      //      m_sv_puppi = -9999;
      //      pt_sv_puppi = -9999;
      //      eta_sv = -9999;
      //      phi_sv = -9999;

      m_sv_mvamet = -9999;
      pt_sv_mvamet = -9999;
      eta_sv_mvamet = -9999;
      phi_sv_mvamet = -9999;

      m_sv_mvamet_uncorr    = -9999;
      m_sv_mvamet_scaleUp   = -9999;
      m_sv_mvamet_scaleDown = -9999;
      m_sv_mvamet_resoUp    = -9999;
      m_sv_mvamet_resoDown  = -9999;
      m_sv_mvamet_tauUp     = -9999;
      m_sv_mvamet_tauDown   = -9999;


      bool validDecayMode = tau_decayMode<=4 || tau_decayMode>=10; // excluding 2prongs
      float tauMass = m_2;
      //      if (tau_decayMode==0) tauMass = pionMass;
      //      std::cout << "MT = " << mt_mva_1 << std::endl;

      if (computeDitauMass && mt_mva_1<45.0 && validDecayMode) {

	if (mvaFound) {
	  // covariance matrix MVAMET
	  TMatrixD covMVAMET(2, 2);
	  covMVAMET[0][0] =  mvacov00;
	  covMVAMET[1][0] =  mvacov01;
	  covMVAMET[0][1] =  mvacov10;
	  covMVAMET[1][1] =  mvacov11;

	  // mva met
	  // define muon 4-vector
	  svFitStandalone::MeasuredTauLepton svFitMuon(svFitStandalone::kTauToMuDecay, 
						       pt_1, eta_1, phi_1, muonMass); 
	  // define tau 4-vector
	  svFitStandalone::MeasuredTauLepton svFitTau(svFitStandalone::kTauToHadDecay,  
						      pt_2, eta_2, phi_2, 
						      tauMass, tau_decayMode); 
	  svFitStandalone::MeasuredTauLepton svFitTauUp(svFitStandalone::kTauToHadDecay,
							(1.0+tauScale)*pt_2, eta_2, phi_2,
							(1.0+tauScale)*tauMass, tau_decayMode); 
	  svFitStandalone::MeasuredTauLepton svFitTauDown(svFitStandalone::kTauToHadDecay,
							  (1.0-tauScale)*pt_2, eta_2, phi_2,
							  (1.0-tauScale)*tauMass, tau_decayMode); 
  
	  // central value
	  SVfitStandaloneAlgorithm algo = SVFitMassComputation(svFitMuon, svFitTau,
							       mvamet_ex, mvamet_ey, 
							       covMVAMET, tau_decayMode, inputFile_visPtResolution);

	  m_sv_mvamet = algo.getMass(); // return value is in units of GeV
	  //	  double transverse_mass = algoY.transverseMass(); // Transverse SVFit mass
	  if ( !algo.isValidSolution() ) 
	    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
	  
	  pt_sv_mvamet  = algo.pt(); 
	  eta_sv_mvamet = algo.eta();
	  phi_sv_mvamet = algo.phi();

	  m_sv_mvamet_uncorr    = m_sv_mvamet;
	  m_sv_mvamet_scaleUp   = m_sv_mvamet;
	  m_sv_mvamet_scaleDown = m_sv_mvamet;
	  m_sv_mvamet_resoUp   = m_sv_mvamet;
	  m_sv_mvamet_resoDown = m_sv_mvamet;
	  m_sv_mvamet_tauUp   = m_sv_mvamet;
	  m_sv_mvamet_tauDown = m_sv_mvamet;

	  if (!isData) { // compute m_sv_mvamet
	    //	    if (isDY||isW) {
	    //	      SVfitStandaloneAlgorithm algo_uncorr = SVFitMassComputation(svFitMuon, svFitTau,
	    //									  mvamet_uncorr_ex, mvamet_uncorr_ey,
	    //									  covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    //	      m_sv_mvamet_uncorr = algo_uncorr.getMass();
	    //	    }
	    // SVfitStandaloneAlgorithm aglo_scaleUp = SVFitMassComputation(svFitMuon, svFitTau,
	    // 								 mvamet_scaleUp_ex, mvamet_scaleUp_ey,
	    // 								 covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    // SVfitStandaloneAlgorithm aglo_resoUp = SVFitMassComputation(svFitMuon, svFitTau,
	    // 								mvamet_resoUp_ex, mvamet_resoUp_ey,
	    // 								covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    SVfitStandaloneAlgorithm aglo_tauUp = SVFitMassComputation(svFitMuon, svFitTauUp,
								       mvamet_ex, mvamet_ey,
								       covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    // SVfitStandaloneAlgorithm aglo_scaleDown = SVFitMassComputation(svFitMuon, svFitTau,
	    // 								   mvamet_scaleDown_ex, mvamet_scaleDown_ey,
	    // 								   covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    // SVfitStandaloneAlgorithm aglo_resoDown = SVFitMassComputation(svFitMuon, svFitTau,
	    // 								  mvamet_resoDown_ex, mvamet_resoDown_ey,
	    // 								  covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    SVfitStandaloneAlgorithm aglo_tauDown = SVFitMassComputation(svFitMuon, svFitTauDown,
									 mvamet_ex, mvamet_ey,
									 covMVAMET, tau_decayMode, inputFile_visPtResolution);
	    

	    //	    m_sv_mvamet_scaleUp = aglo_scaleUp.getMass();
	    //	    m_sv_mvamet_resoUp  = aglo_resoUp.getMass();
	    m_sv_mvamet_tauUp   = aglo_tauUp.getMass();

	    //	    m_sv_mvamet_scaleDown = aglo_scaleDown.getMass();
	    //	    m_sv_mvamet_resoDown  = aglo_resoDown.getMass();
	    m_sv_mvamet_tauDown   = aglo_tauDown.getMass();

	    //	    std::cout << "msv = " << m_sv_mvamet << "   uncorrected : " << m_sv_mvamet_uncorr << std::endl;
	    //	    std::cout << "tauES    : up = " << m_sv_mvamet_tauUp << "   down = " << m_sv_mvamet_tauDown << std::endl;
	    //	    std::cout << "metScale : up = " << m_sv_mvamet_scaleUp << "   down = " << m_sv_mvamet_scaleDown << std::endl;
	    //	    std::cout << "metReso  : up = " << m_sv_mvamet_resoUp << "   down = " << m_sv_mvamet_resoDown << std::endl;


	  }
	}
	else {
	  std::cout << "MVA MET is not found : computation of svFit mass impossible" << std::endl;
	}
      }

      byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
      byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      againstElectronLooseMVA5_2 = analysisTree.tau_againstElectronLooseMVA5[tauIndex];
      againstElectronMediumMVA5_2 = analysisTree.tau_againstElectronMediumMVA5[tauIndex];
      againstElectronTightMVA5_2 = analysisTree.tau_againstElectronTightMVA5[tauIndex];
      againstElectronVLooseMVA5_2 = analysisTree.tau_againstElectronVLooseMVA5[tauIndex];
      againstElectronVTightMVA5_2 = analysisTree.tau_againstElectronVTightMVA5[tauIndex];
      againstMuonLoose3_2 = analysisTree.tau_againstMuonLoose3[tauIndex];
      againstMuonTight3_2 = analysisTree.tau_againstMuonTight3[tauIndex];
      againstElectronVLooseMVA6_2 = analysisTree.tau_againstElectronVLooseMVA6[tauIndex];
      byIsolationMVArun2v1DBoldDMwLTraw_2 = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tauIndex];
      byTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.tau_px[tauIndex],
					  analysisTree.tau_py[tauIndex],
					  analysisTree.tau_pz[tauIndex],
					  analysisTree.tau_e[tauIndex]);


      TLorentzVector tauUpLV; tauUpLV.SetXYZT((1.0+tauScale)*analysisTree.tau_px[tauIndex],
					      (1.0+tauScale)*analysisTree.tau_py[tauIndex],
					      (1.0+tauScale)*analysisTree.tau_pz[tauIndex],
					      (1.0+tauScale)*analysisTree.tau_e[tauIndex]);

      TLorentzVector tauDownLV; tauDownLV.SetXYZT((1.0-tauScale)*analysisTree.tau_px[tauIndex],
						  (1.0-tauScale)*analysisTree.tau_py[tauIndex],
						  (1.0-tauScale)*analysisTree.tau_pz[tauIndex],
						  (1.0-tauScale)*analysisTree.tau_e[tauIndex]);



      TLorentzVector dileptonLV = muonLV + tauLV;

      // visible mass
      m_vis         = dileptonLV.M();
      m_vis_tauUp   = (muonLV+tauUpLV).M();
      m_vis_tauDown = (muonLV+tauDownLV).M();
      m_vis_scaleUp   = dileptonLV.M();
      m_vis_scaleDown = dileptonLV.M();
      m_vis_resoUp    = dileptonLV.M();
      m_vis_resoDown  = dileptonLV.M();

      TLorentzVector mvametLV; mvametLV.SetXYZM(mvamet_ex,mvamet_ey,0.,0.);
      TLorentzVector mvametResoUpLV; mvametResoUpLV.SetXYZM(mvamet_resoUp_ex,mvamet_resoUp_ey,0.,0.);
      TLorentzVector mvametResoDownLV; mvametResoDownLV.SetXYZM(mvamet_resoDown_ex,mvamet_resoDown_ey,0.,0.);
      TLorentzVector mvametScaleUpLV; mvametScaleUpLV.SetXYZM(mvamet_scaleUp_ex,mvamet_scaleUp_ey,0.,0.);
      TLorentzVector mvametScaleDownLV; mvametScaleDownLV.SetXYZM(mvamet_scaleDown_ex,mvamet_scaleDown_ey,0.,0.);

      // computing total transverse mass
      mTtot = totalTransverseMass         ( muonLV ,     tauLV ,     mvametLV);

      mTtot_tauUp   =  totalTransverseMass( muonLV ,     tauUpLV ,   mvametLV);
      mTtot_tauDown =  totalTransverseMass( muonLV ,     tauDownLV , mvametLV);

      mTtot_scaleUp   =  totalTransverseMass ( muonLV ,  tauLV , mvametScaleUpLV);
      mTtot_scaleDown =  totalTransverseMass ( muonLV ,  tauLV , mvametScaleDownLV);

      mTtot_resoUp    =  totalTransverseMass ( muonLV ,  tauLV , mvametResoUpLV);
      mTtot_resoDown  =  totalTransverseMass ( muonLV ,  tauLV , mvametResoDownLV);

      // visible ditau pt 
      pt_tt = dileptonLV.Pt();

      // bisector of electron and muon transverse momenta
      float tauUnitX = tauLV.Px()/tauLV.Pt();
      float tauUnitY = tauLV.Py()/tauLV.Pt();
	
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + tauLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + tauLV.Py();
      
      float vectorVisX = muonLV.Px() + tauLV.Px();
      float vectorVisY = muonLV.Py() + tauLV.Py();

      // computation of DZeta variable
      pzetamiss = mvamet_ex*zetaX + mvamet_ey*zetaY;
      pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;

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



