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
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

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
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/CalibrationOfImpactParameters.h"
#include "DesyTauAnalyses/NTupleMaker/interface/helper_functions_em.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ggF_qcd_uncertainty_2017.h"

#include "TCut.h"

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "DesyTauAnalyses/NTupleMaker/interface/btagSF.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"


int main(int argc, char * argv[]) {
    
   // first argument - config file
    // second argument - filelist
    
    bool sync = false;
   
    using namespace std;
    
    // **** configuration
    Config cfg(argv[1]);
    
    const bool computeSVFitMass = cfg.get<bool>("ComputeSVFitMass");
    const bool removeGammaStar = cfg.get<bool>("RemoveGammaStar");
    
    const bool isData = cfg.get<bool>("IsData");
    const bool isDY   = cfg.get<bool>("IsDY");
    const bool isW    = cfg.get<bool>("IsW");
    const bool isTOP    = cfg.get<bool>("IsTOP");
    const bool isEmbedded = cfg.get<bool>("IsEmbedded");
    const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
    const string jsonFile = cfg.get<string>("jsonFile");
    const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
    const bool isSignal = cfg.get<bool>("IsSignal");

    // kinematic cuts on electrons
    const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
    const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
    const float etaElectronCut     = cfg.get<float>("etaElectronCut");
    const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
    const float dzElectronCut      = cfg.get<float>("dzElectronCut");
    const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
    const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
    const bool applySpring16ElectronId     = cfg.get<bool>("ApplySpring16ElectronId");
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
    const bool applyICHEPMuonId     = cfg.get<bool>("ApplyICHEPMuonId");
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
    const bool isMuonIsoR03 = cfg.get<bool>("IsMuonIsoR03");
    const bool isElectronIsoR03 = cfg.get<bool>("IsElectronIsoR03");
    const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
    const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
    const bool applyDzFilterMatch = cfg.get<bool>("ApplyDzFilterMatch");
    const string mu23ele12DzFilter = cfg.get<string>("Mu23Ele12DzFilter");
    const string mu8ele23DzFilter = cfg.get<string>("Mu8Ele23DzFilter");
    
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
    
    TString Mu23Ele12DzFilter(mu23ele12DzFilter);
    TString Mu8Ele23DzFilter(mu8ele23DzFilter);

    TString BTagDiscriminator(bTagDiscriminator);
    
    const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
    const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");
    
    const string Muon23TriggerFile = cfg.get<string>("Muon23TriggerEff");
    const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");
    
    const string Electron23TriggerFile = cfg.get<string>("Electron23TriggerEff");
    const string Electron12TriggerFile = cfg.get<string>("Electron12TriggerEff");
    
    const string trackingSFFile = cfg.get<string>("TrackingSFFile");
    
    const float muonScale = cfg.get<float>("MuonScale");
    const float eleScaleBarrel = cfg.get<float>("EleScaleBarrel");
    const float eleScaleEndcap = cfg.get<float>("EleScaleEndcap");
    
    //const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
    //TString RecoilMvaFileName(recoilMvaFileName);
    const string recoilFileName   = cfg.get<string>("RecoilFileName");
    TString RecoilFileName(recoilFileName);
    
    const string metSysFileName   = cfg.get<string>("MetSysFileName");
    TString MetSysFileName(metSysFileName);
    
    const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
    TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
    
    const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
    TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
    
    const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
    const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");
    
    TString PileUpDataFile(pileUpDataFile);
    TString PileUpMCFile(pileUpMCFile);
    
    const string correctionWSFile_embedded = cfg.get<string>("CorrectionWSFile_embedded");
    const string correctionWSFile_embedded_trigger = cfg.get<string>("CorrectionWSFile_trigger");

    const bool apply_ggh_reweighting = cfg.get<bool>("ApplygghReweighting");
    // **** end of configuration

    const float a_jetMu = 0.902;
    const float b_jetMu = 0.0025;

    const float a_jetMuUp = 0.992;
    const float b_jetMuUp = 0.005;

    const float a_jetMuDown = 0.812;
    const float b_jetMuDown = 0.0;

    const float a_jetEle = 0.794;
    const float b_jetEle = 0.0015;

    const float a_jetEleUp = 0.883;
    const float b_jetEleUp = 0.003;

    const float a_jetEleDown = 0.702;
    const float b_jetEleDown = 0.0;

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
    Float_t         effweight_jetMuUp;
    Float_t         effweight_jetMuDown;
    Float_t         effweight_jetEUp;
    Float_t         effweight_jetEDown;
    Float_t         fakeweight;
    Float_t         embeddedWeight;
    Float_t         signalWeight;
    Float_t         topptweight;
    Float_t         topptweightRun2;

    Float_t         qcdweight;
    Float_t         qcdweightup;
    Float_t         qcdweightdown;
    
    Float_t         qcdweight_nodzeta;
    Float_t         qcdweightup_nodzeta;
    Float_t         qcdweightdown_nodzeta;
    
    Float_t         qcdweight_0jet_rate_up;
    Float_t         qcdweight_0jet_rate_down;
    Float_t         qcdweight_1jet_rate_up;
    Float_t         qcdweight_1jet_rate_down;
    Float_t         qcdweight_0jet_shape_up;
    Float_t         qcdweight_0jet_shape_down;
    Float_t         qcdweight_1jet_shape_up;
    Float_t         qcdweight_1jet_shape_down;
    
    Float_t         qcdweight_iso_up;
    Float_t         qcdweight_iso_down;


    Float_t         zptmassweight;
    Float_t         zptmassweight_esup;
    Float_t         zptmassweight_esdown;
    Float_t         zptmassweight_ttup;
    Float_t         zptmassweight_ttdown;
    Float_t         zptmassweight_statpt0up;
    Float_t         zptmassweight_statpt0down;
    Float_t         zptmassweight_statpt40up;
    Float_t         zptmassweight_statpt40down;
    Float_t         zptmassweight_statpt80up;
    Float_t         zptmassweight_statpt80down;

    Float_t         weight;
    Float_t         weight_ggh_NNLOPS;
    
    Float_t         m_vis;
    Float_t         pt_vis;

    Float_t         mTtot;
    Float_t         mCDF;

    Float_t         mTdileptonMET;
    Float_t         mTemu;

    Float_t         m_sv;
    Float_t         mt_sv;

    Float_t         pt_sv;
    Float_t         eta_sv;
    Float_t         phi_sv;

    Float_t         pt_sv_gen;
    Float_t         eta_sv_gen;
    Float_t         phi_sv_gen;

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
    
    Float_t         mtmax;
    
    Bool_t          os;
    Bool_t          dilepton_veto;
    Bool_t          extraelec_veto;
    Bool_t          extramuon_veto;

    Bool_t          trg_muonelectron;
    
    Float_t         met;
    Float_t         metphi;
    Float_t         metcov00;
    Float_t         metcov01;
    Float_t         metcov10;
    Float_t         metcov11;

    Float_t         met_recoilscaleUp;
    Float_t         met_recoilscaleDown;
    Float_t         met_recoilresoUp;
    Float_t         met_recoilresoDown;
    
    Float_t         met_uncorr;
    Float_t         metphi_uncorr;
    
    Float_t         met_unclMetUp;
    Float_t         metphi_unclMetUp;
    
    Float_t         met_unclMetDown;
    Float_t         metphi_unclMetDown;

    Float_t         genmet;
    Float_t         genmetphi;
    
    Float_t         msvmet;
    Float_t         msvmetphi;
    
    Float_t         pt_tt;
    Float_t         dr_tt;
    Float_t         dphi_tt;
    
    Float_t         pt_ttjj;

    Float_t         pzetavis;
    Float_t         pzetamiss;
    Float_t         dzeta;
    
    Float_t         pzetamiss_genmet;
    Float_t         dzeta_genmet;
    
    Float_t         mva_gf;
    
    Int_t           njets;
    Int_t           njets_jecUncEta0To5Up;
    Int_t           njets_jecUncEta0To5Down;
    Int_t           njets_jecUncEta0To3Up;
    Int_t           njets_jecUncEta0To3Down;
    Int_t           njets_jecUncEta3To5Up;
    Int_t           njets_jecUncEta3To5Down;
    Int_t           njets_jecUncRelativeBalUp;
    Int_t           njets_jecUncRelativeBalDown;

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
    Float_t         dijetphi;
    Float_t         dijetpt;
    Float_t         jdeta;
    Int_t           njetingap;
    
    Int_t           nbtag;
    Int_t           nbtag_mistagUp;
    Int_t           nbtag_mistagDown;
    Int_t           nbtag_btagUp;
    Int_t           nbtag_btagDown;
    Int_t           nbtag_noSF;
    Double_t         bpt;
    Double_t         beta;
    Double_t         bphi;
    
    Float_t         nuPx;  // neutrinos from t -> W(lv) + b
    Float_t         nuPy;  // or from tau->l+v+v events
    Float_t         nuPz;  // (x,y,z) components
    Float_t         nuPt;  // pT
    Float_t         nuPhi; // phi
    
    Float_t         nuPx_msv; // neutrinos after msv fit
    Float_t         nuPy_msv; //
    Float_t         nuPz_msv; //
    Float_t         nuPt_msv; //
    Float_t         nuPhi_msv; //
    
    Float_t         msv_gen;
    Float_t         mtsv_gen;
    Float_t         mtBoson_gen;
    
    Float_t         mTtot_gen;
    Float_t         mTemu_gen;
    Float_t         mTemet_gen;
    Float_t         mTmumet_gen;
    
    Float_t         bdt;
    Float_t         bdt_ggh;
    Float_t         bdt_bbh;
    
    Float_t         lepPx;
    Float_t         lepPy;
    Float_t         lepPz;
    
    Float_t         bosonPx;
    Float_t         bosonPy;
    Float_t         bosonPz;
    Float_t         bosonPt;
    Float_t         bosonMass;
    
    Float_t         dphi_mumet;
    Float_t         dphi_emet;
    
    UInt_t          npartons;
    
    Bool_t isZLL;
    Bool_t isZMM;
    Bool_t isZEE;
    Bool_t isZTT;
    
    Bool_t metFilters_;

    Bool_t badChargedCandidateFilter_;
    Bool_t badPFMuonFilter_;
    Bool_t badGlobalMuonFilter_;
    Bool_t muonBadTrackFilter_;
    Bool_t chargedHadronTrackResolutionFilter_;

    Bool_t badMuonFilter_;
    Bool_t duplicateMuonFilter_;

    Float_t weightScale1;
    Float_t weightScale2;
    Float_t weightScale3;
    Float_t weightScale4;
    Float_t weightScale5;
    Float_t weightScale6;
    Float_t weightScale7;
    Float_t weightScale8;

    Float_t weightPDFup;
    Float_t weightPDFdown;
    
    Float_t d0_1_cal;
    Float_t d0_2_cal;
    Float_t dZ_1_cal;
    Float_t dZ_2_cal;

    Bool_t veto_embedded;

    Float_t higgspt_HTXS;
    Int_t njets_HTXS;
    Int_t htxs_stage0cat;
    Int_t htxs_stage1cat;

    Float_t THU_ggH_Mu;
    Float_t THU_ggH_Res;
    Float_t THU_ggH_Mig01;
    Float_t THU_ggH_Mig12;
    Float_t THU_ggH_VBF2j;
    Float_t THU_ggH_VBF3j;
    Float_t THU_ggH_PT60;
    Float_t THU_ggH_PT120;
    Float_t THU_ggH_qmtop;

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
    tree->Branch("veto_embedded",&veto_embedded,"veto_embedded/O");
     
    tree->Branch("weightScale1",&weightScale1,"weightScale1/F");
    tree->Branch("weightScale2",&weightScale2,"weightScale2/F");
    tree->Branch("weightScale3",&weightScale3,"weightScale3/F");
    tree->Branch("weightScale4",&weightScale4,"weightScale4/F");
    tree->Branch("weightScale5",&weightScale5,"weightScale5/F");
    tree->Branch("weightScale6",&weightScale6,"weightScale6/F");
    tree->Branch("weightScale7",&weightScale7,"weightScale7/F");
    tree->Branch("weightScale8",&weightScale8,"weightScale8/F");

    tree->Branch("weightPDFup",&weightPDFup,"weightPDFup/F");
    tree->Branch("weightPDFdown",&weightPDFdown,"weightPDFdown/F");

    tree->Branch("weight_ggh_NNLOPS", &weight_ggh_NNLOPS, "weight_ggh_NNLOPS/F");
    tree->Branch("THU_ggH_Mu", &THU_ggH_Mu, "THU_ggH_Mu/F");
    tree->Branch("THU_ggH_Res", &THU_ggH_Res, "THU_ggH_Res/F");
    tree->Branch("THU_ggH_Mig01", &THU_ggH_Mig01, "THU_ggH_Mig01/F");
    tree->Branch("THU_ggH_Mig12", &THU_ggH_Mig12, "THU_ggH_Mig12/F");
    tree->Branch("THU_ggH_VBF2j", &THU_ggH_VBF2j, "THU_ggH_VBF2j/F");
    tree->Branch("THU_ggH_VBF3j", &THU_ggH_VBF3j, "THU_ggH_VBF3j/F");
    tree->Branch("THU_ggH_PT60" , &THU_ggH_PT60, "THU_ggH_PT60/F");
    tree->Branch("THU_ggH_PT120", &THU_ggH_PT120, "THU_ggH_PT120/F");
    tree->Branch("THU_ggH_qmtop", &THU_ggH_qmtop, "THU_ggH_qmtop/F");

    tree->Branch("higgspt_HTXS",&higgspt_HTXS,"higgspt_HTXS/F");
    tree->Branch("njets_HTXS",&njets_HTXS,"njets_HTXS/I");
    tree->Branch("htxs_stage0cat",&htxs_stage0cat,"htxs_stage0cat/I");
    tree->Branch("htxs_stage1cat",&htxs_stage1cat,"htxs_stage1cat/I");

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
    tree->Branch("effweight_jetMuUp", &effweight_jetMuUp, "effweight_jetMuUp/F");
    tree->Branch("effweight_jetMuDown", &effweight_jetMuDown, "effweight_jetMuDown/F");
    tree->Branch("effweight_jetEUp", &effweight_jetEUp, "effweight_jetEUp/F");
    tree->Branch("effweight_jetEDown", &effweight_jetEDown, "effweight_jetEDown/F");
    tree->Branch("fakeweight", &fakeweight, "fakeweight/F");
    tree->Branch("embeddedWeight", &embeddedWeight, "embeddedWeight/F");
    tree->Branch("signalWeight", &signalWeight, "signalWeight/F");
    tree->Branch("topptweight", &topptweight, "topptweight/F");
    tree->Branch("topptweightRun2", &topptweightRun2, "topptweightRun2/F");

    tree->Branch("qcdweight", &qcdweight, "qcdweight/F");
    tree->Branch("qcdweightup", &qcdweightup, "qcdweightup/F");
    tree->Branch("qcdweightdown", &qcdweightdown, "qcdweightdown/F");
    
    tree->Branch("qcdweight_nodzeta", &qcdweight_nodzeta, "qcdweight_nodzeta/F");
    tree->Branch("qcdweightup_nodzeta", &qcdweightup_nodzeta, "qcdweightup_nodzeta/F");
    tree->Branch("qcdweightdown_nodzeta", &qcdweightdown_nodzeta, "qcdweightdown_nodzeta/F");
    
    tree->Branch("qcdweight_0jet_rate_up",&qcdweight_0jet_rate_up,"qcdweight_0jet_rate_up/F");
    tree->Branch("qcdweight_0jet_rate_down ",&qcdweight_0jet_rate_down,"qcdweight_0jet_rate_down/F");
    tree->Branch("qcdweight_1jet_rate_up ",&qcdweight_1jet_rate_up,"qcdweight_1jet_rate_up/F");
    tree->Branch("qcdweight_1jet_rate_down",&qcdweight_1jet_rate_down,"qcdweight_1jet_rate_down/F");
    tree->Branch("qcdweight_0jet_shape_up",&qcdweight_0jet_shape_up,"qcdweight_0jet_shape_up/F");
    tree->Branch("qcdweight_0jet_shape_down",&qcdweight_0jet_shape_down,"qcdweight_0jet_shape_down/F");  
    tree->Branch("qcdweight_1jet_shape_up ",&qcdweight_1jet_shape_up,"qcdweight_1jet_shape_up/F");
    tree->Branch("qcdweight_1jet_shape_down",&qcdweight_1jet_shape_down,"qcdweight_1jet_shape_down/F");

    tree->Branch("qcdweight_iso_up",&qcdweight_iso_up,"qcdweight_iso_up/F");
    tree->Branch("qcdweight_iso_down",&qcdweight_iso_down,"qcdweight_iso_down/F");

    tree->Branch("zptmassweight",&zptmassweight,"zptmassweight/F");

    tree->Branch("zptmassweight_esup",&zptmassweight_esup,"zptmassweight_esup/F");
    tree->Branch("zptmassweight_esdown",&zptmassweight_esdown,"zptmassweight_esdown/F");

    tree->Branch("zptmassweight_ttup",&zptmassweight_ttup,"zptmassweight_ttup/F");
    tree->Branch("zptmassweight_ttdown",&zptmassweight_ttdown,"zptmassweight_ttdown/F");

    tree->Branch("zptmassweight_statpt0up",&zptmassweight_statpt0up,"zptmassweight_statpt0up/F");
    tree->Branch("zptmassweight_statpt0down",&zptmassweight_statpt0down,"zptmassweight_statpt0down/F");

    tree->Branch("zptmassweight_statpt40up",&zptmassweight_statpt40up,"zptmassweight_statpt40up/F");
    tree->Branch("zptmassweight_statpt40down",&zptmassweight_statpt40down,"zptmassweight_statpt40down/F");

    tree->Branch("zptmassweight_statpt80up",&zptmassweight_statpt80up,"zptmassweight_statpt80up/F");
    tree->Branch("zptmassweight_statpt80down",&zptmassweight_statpt80down,"zptmassweight_statpt80down/F");

    tree->Branch("weight", &weight, "weight/F");
    
    tree->Branch("metFilters",&metFilters_,"metFilters/O");
    tree->Branch("trg_muonelectron",&trg_muonelectron,"trg_muonelectron/O");

    tree->Branch("badChargedCandidateFilter",&badChargedCandidateFilter_,"badChargedCandidateFilter/O");
    tree->Branch("badPFMuonFilter",&badPFMuonFilter_,"badPFMuonFilter/O");
    tree->Branch("badGlobalMuonFilter",&badGlobalMuonFilter_,"badGlobalMuonFilter/O");
    tree->Branch("muonBadTrackFilter",&muonBadTrackFilter_,"muonBadTrackFilter/O");
    tree->Branch("chargedHadronTrackResolutionFilter",&chargedHadronTrackResolutionFilter_,"chargedHadronTrackResolutionFilter/O");

    tree->Branch("badMuonFilter",&badMuonFilter_,"badMuonFilter/O");
    tree->Branch("duplicateMuonFilter",&duplicateMuonFilter_,"duplicateMuonFilter/O");
    
    tree->Branch("m_vis",        &m_vis,        "m_vis/F");
    tree->Branch("pt_vis",        &pt_vis,        "pt_vis/F");
    tree->Branch("mTtot",        &mTtot,        "mTtot/F");
    tree->Branch("mCDF",         &mCDF,         "mCDF/F");
    
    tree->Branch("mTdileptonMET", &mTdileptonMET, "mTdileptonMET/F");

    tree->Branch("mTemu",        &mTemu,        "mTemu/F");

    tree->Branch("m_sv",    &m_sv,   "m_sv/F");
    tree->Branch("mt_sv",   &mt_sv,  "mt_sv/F");

    tree->Branch("pt_sv",   &pt_sv,  "pt_sv/F");
    tree->Branch("eta_sv",  &eta_sv, "eta_sv/F");
    tree->Branch("phi_sv",  &phi_sv, "phi_sv/F");
    
    tree->Branch("pt_sv_gen",   &pt_sv_gen,  "pt_sv_gen/F");
    tree->Branch("eta_sv_gen",  &eta_sv_gen, "eta_sv_gen/F");
    tree->Branch("phi_sv_gen",  &phi_sv_gen, "phi_sv_gen/F");
    
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
    tree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");
    
    tree->Branch("mtmax", &mtmax, "mtmax/F");
    
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

    tree->Branch("met_recoilscaleUp", &met_recoilscaleUp, "met_recoilscaleUp/F");
    tree->Branch("met_recoilscaleDown", &met_recoilscaleDown, "met_recoilscaleDown/F");
    tree->Branch("met_recoilresoUp", &met_recoilresoUp, "met_recoilresoUp/F");
    tree->Branch("met_recoilresoDown", &met_recoilresoDown, "met_recoilresoDown/F");
    
    tree->Branch("met_uncorr", &met_uncorr, "met_uncorr/F");
    tree->Branch("metphi_uncorr", &metphi_uncorr, "metphi_uncorr/F");
    
    tree->Branch("met_unclMetUp", &met_unclMetUp, "met_unclMetUp/F");
    tree->Branch("metphi_unclMetUp", &metphi_unclMetUp, "metphi_unclMetUp/F");
    
    tree->Branch("met_unclMetDown", &met_unclMetDown, "met_unclMetDown/F");
    tree->Branch("metphi_unclMetDown", &metphi_unclMetDown, "metphi_unclMetDown/F");

    tree->Branch("genmet", &genmet, "genmet/F");
    tree->Branch("genmetphi", &genmetphi, "genmetphi/F");
    
    tree->Branch("mTtot_gen",  &mTtot_gen,   "mTtot_gen/F");
    tree->Branch("mTemu_gen",  &mTemu_gen,   "mTemu_gen/F");
    tree->Branch("mTemet_gen", &mTemet_gen,  "mTemet_gen/F");
    tree->Branch("mTmumet_gen",&mTmumet_gen, "mTmumet_gen/F");
    tree->Branch("msv_gen",&msv_gen,"msv_gen/F");
    tree->Branch("mtsv_gen",&mtsv_gen,"mtsv_gen/F");
    tree->Branch("mtBoson_gen",&mtBoson_gen,"mtBoson_gen/F");
    
    tree->Branch("dphi_mumet",&dphi_mumet,"dphi_mumet/F");
    tree->Branch("dphi_emet",&dphi_emet,"dphi_emet/F");
    
    tree->Branch("msvmet", &msvmet, "msvmet/F");
    tree->Branch("msvmetphi", &msvmetphi, "msvmetphi/F");
    
    tree->Branch("pt_tt", &pt_tt, "pt_tt/F");

    tree->Branch("dr_tt", &dr_tt, "dr_tt/F");
    tree->Branch("dphi_tt", &dphi_tt, "dphi_tt/F");
    tree->Branch("pt_ttjj", &pt_ttjj, "pt_ttjj/F");
    tree->Branch("pzetavis", &pzetavis, "pzetavis/F");

    tree->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");

    tree->Branch("dzeta",&dzeta,"dzeta/F");    

    tree->Branch("pzetamiss_genmet", &pzetamiss_genmet, "pzetamiss_genmet/F");
    tree->Branch("dzeta_genmet",&dzeta_genmet,"dzeta_genmet/F");
    
    tree->Branch("mva_gf", &mva_gf, "mva_gf/F");
    
    tree->Branch("njets", &njets, "njets/I");
    tree->Branch("njets_jecUncEta0To5Up", &njets_jecUncEta0To5Up, "njets_jecUncEta0To5Up/I");
    tree->Branch("njets_jecUncEta0To5Down", &njets_jecUncEta0To5Down, "njets_jecUncEta0To5Down/I");
    tree->Branch("njets_jecUncEta0To3Up", &njets_jecUncEta0To3Up, "njets_jecUncEta0To3Up/I");
    tree->Branch("njets_jecUncEta0To3Down", &njets_jecUncEta0To3Down, "njets_jecUncEta0To3Down/I");
    tree->Branch("njets_jecUncEta3To5Up", &njets_jecUncEta3To5Up, "njets_jecUncEta3To5Up/I");
    tree->Branch("njets_jecUncEta3To5Down", &njets_jecUncEta3To5Down, "njets_jecUncEta3To5Down/I");
    tree->Branch("njets_jecUncRelativeBalUp", &njets_jecUncRelativeBalUp, "njets_jecUncRelativeBalUp/I");
    tree->Branch("njets_jecUncRelativeBalDown", &njets_jecUncRelativeBalDown, "njets_jecUncRelativeBalDown/I");

    tree->Branch("njetspt20", &njetspt20, "njetspt20/I");

    tree->Branch("bdt",&bdt,"bdt/F");
    tree->Branch("bdt_bbh",&bdt_bbh,"bdt_bbh/F");
    tree->Branch("bdt_ggh",&bdt_ggh,"bdt_ggh/F");
    
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
    tree->Branch("dijetphi", &dijetphi, "dijetphi/F");
    tree->Branch("dijetpt", &dijetpt, "dijetpt/F");
    tree->Branch("jdeta", &jdeta, "jdeta/F");
    tree->Branch("njetingap", &njetingap, "njetingap/I");
    
    tree->Branch("nbtag", &nbtag, "nbtag/I");
    tree->Branch("nbtag_mistagUp", &nbtag_mistagUp, "nbtag_mistagUp/I");
    tree->Branch("nbtag_mistagDown", &nbtag_mistagDown, "nbtag_mistagDown/I");
    tree->Branch("nbtag_btagUp", &nbtag_btagUp, "nbtag_btagUp/I");
    tree->Branch("nbtag_btagDown", &nbtag_btagDown, "nbtag_btagDown/I");
    tree->Branch("nbtag_noSF", &nbtag_noSF, "nbtag_noSF/I");
    tree->Branch("bpt_1",   &bpt,   "bpt/D");
    tree->Branch("beta_1",  &beta,  "beta/D");
    tree->Branch("bphi",  &bphi,  "bphi/F");
    
    tree->Branch("nuPx",&nuPx,"nuPx/F");
    tree->Branch("nuPy",&nuPy,"nuPy/F");
    tree->Branch("nuPz",&nuPz,"nuPz/F");
    tree->Branch("nuPt",&nuPt,"nuPt/F");
    tree->Branch("nuPhi",&nuPhi,"nuPhi/F");
    
    tree->Branch("nuPx_msv",&nuPx_msv,"nuPx_msv/F");
    tree->Branch("nuPy_msv",&nuPy_msv,"nuPy_msv/F");
    tree->Branch("nuPz_msv",&nuPz_msv,"nuPz_msv/F");
    tree->Branch("nuPt_msv",&nuPt_msv,"nuPt_msv/F");
    tree->Branch("nuPhi_msv",&nuPhi_msv,"nuPhi_msv/F");
    
    tree->Branch("lepPx",&lepPx,"lepPx/F");
    tree->Branch("lepPy",&lepPy,"lepPy/F");
    tree->Branch("lepPz",&lepPz,"lepPz/F");
    
    tree->Branch("bosonPx",&bosonPx,"bosonPx/F");
    tree->Branch("bosonPy",&bosonPy,"bosonPy/F");
    tree->Branch("bosonPz",&bosonPz,"bosonPz/F");
    tree->Branch("bosonPt",&bosonPt,"bosonPt/F");
    tree->Branch("bosonMass",&bosonMass,"bosonMass/F");
    
    tree->Branch("npartons",&npartons,"npartons/i");

    tree->Branch("d0_1_cal",&d0_1_cal,"d0_1_cal/F");
    tree->Branch("d0_2_cal",&d0_2_cal,"d0_2_cal/F");
    tree->Branch("dZ_1_cal",&dZ_1_cal,"dZ_1_cal/F");
    tree->Branch("dZ_2_cal",&dZ_2_cal,"dZ_2_cal/F");

    string cmsswBase = (getenv ("CMSSW_BASE"));

    // START : Prepare uncertainties (NEW) ================================================================================================================

    vector<TString> unc_vars= {"met",  // 0
			       "metphi", // 1
			       "mTtot",  // 2
			       "mTdileptonMET",  // 3
			       "pt_tt",  // 4
			       "pt_ttjj", // 5
			       "pzetamiss", // 6
			       "dzeta", // 7
			       "mt_1", // 8
			       "mt_2", // 9
			       "mtmax", // 10
			       "dphi_emet", // 11
			       "dphi_mumet", // 12
			       "pzetavis", // 13
			       "m_vis", // 14
			       "pt_vis", // 15
			       "pt_1", // 16
			       "pt_2", // 17
			       "jpt_1", // 18
			       "jpt_2", // 19
			       "mjj", // 20
			       "dijetphi", // 21
			       "dijetpt", // 22
			       "m_sv", // 23
			       "mTemu"}; // 24

    struct inputs {
      float container[50]; // this array is used as a container for the shifted variables in the propagate_uncertainty function
      TLorentzVector electronLV;
      TLorentzVector muonLV;
      TLorentzVector metLV;
      TLorentzVector jet1LV;
      TLorentzVector jet2LV;
    };

    inputs unclMetUp;
    inputs unclMetDown;
    inputs escaleUp;
    inputs escaleDown;
    //inputs mscaleUp;
    //inputs mscaleDown;
    inputs recoilscaleUp;
    inputs recoilscaleDown;
    inputs recoilresoUp;
    inputs recoilresoDown;
    inputs jecUncEta0To5Up;
    inputs jecUncEta0To5Down;
    inputs jecUncEta0To3Up;
    inputs jecUncEta0To3Down;
    inputs jecUncEta3To5Up;
    inputs jecUncEta3To5Down;
    inputs jecUncRelativeBalUp;
    inputs jecUncRelativeBalDown;

    map<TString, inputs> uncertainty_map = { { "unclMetUp" , unclMetUp },
					     { "unclMetDown" , unclMetDown },
					     { "escaleUp" , escaleUp },
					     { "escaleDown" , escaleDown },
                  //{ "mscaleUp" , mscaleUp },
					   //{ "mscaleDown" , mscaleDown },
					     { "recoilscaleUp" , recoilscaleUp },
					     { "recoilscaleDown" , recoilscaleDown },
					     { "recoilresoUp" , recoilresoUp },
					     { "recoilresoDown" , recoilresoDown },
					     { "jecUncEta0To5Up" , jecUncEta0To5Up },
					     { "jecUncEta0To5Down" , jecUncEta0To5Down },
					     { "jecUncEta0To3Up" , jecUncEta0To3Up },
					     { "jecUncEta0To3Down" , jecUncEta0To3Down },
					     { "jecUncEta3To5Up" , jecUncEta3To5Up },
					     { "jecUncEta3To5Down" , jecUncEta3To5Down },
					     { "jecUncRelativeBalUp" , jecUncRelativeBalUp },
					     { "jecUncRelativeBalDown" , jecUncRelativeBalDown },
    };

    for(auto &uncert : uncertainty_map){
      int count_unc = 0;
      for(auto &unc_var : unc_vars){
	tree->Branch(unc_var+"_"+uncert.first,
		     &uncert.second.container[count_unc],
		     unc_var+"_"+uncert.first+"/F");
	count_unc += 1;
      }
    }
    // START : Prepare uncertainties (NEW) ================================================================================================================

    // START : JEC uncertainties (NEW) ================================================================================================================
    const int nsrc_Eta0To5 = 13;
    const char* srcnames_Eta0To5[nsrc_Eta0To5] = {"SinglePionECAL",
						  "SinglePionHCAL",
						  "AbsoluteFlavMap",
						  "AbsoluteMPFBias",
						  "AbsoluteScale",
						  "AbsoluteStat",
						  "Fragmentation",
						  "FlavorQCD",
						  "TimePtEta",
						  "PileUpDataMC",
						  "RelativeFSR",
						  "RelativeStatFSR",
						  "PileUpPtRef"};
    const int nsrc_Eta0To3 = 9;
    const char* srcnames_Eta0To3[nsrc_Eta0To3] = {"PileUpPtEC1",
						  "PileUpPtEC2",
						  "PileUpPtBB",
						  "RelativeJEREC1",
						  "RelativeJEREC2",
						  "RelativePtEC1",
						  "RelativePtEC2",
						  "RelativeStatEC",
						  "RelativePtBB"};
    const int nsrc_Eta3To5 = 4;
    const char* srcnames_Eta3To5[nsrc_Eta3To5] = {"RelativeStatHF",
						  "RelativePtHF",
						  "PileUpPtHF",
						  "RelativeJERHF"};
    const int nsrc_RelativeBal = 1;
    const char* srcnames_RelativeBal[nsrc_RelativeBal] = {"RelativeBal"};

    std::vector<JetCorrectionUncertainty*> vsrc_Eta0To5(nsrc_Eta0To5);
    std::vector<JetCorrectionUncertainty*> vsrc_Eta0To3(nsrc_Eta0To3);
    std::vector<JetCorrectionUncertainty*> vsrc_Eta3To5(nsrc_Eta3To5);
    std::vector<JetCorrectionUncertainty*> vsrc_RelativeBal(nsrc_RelativeBal);

    for (int isrc = 0; isrc < nsrc_Eta0To5; isrc++) {
      const char *name = srcnames_Eta0To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/Summer16_23Sep2016V4_DATA_UncertaintySources_AK4PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To5[isrc] = unc;
    }
    for (int isrc = 0; isrc < nsrc_Eta0To3; isrc++) {
      const char *name = srcnames_Eta0To3[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/Summer16_23Sep2016V4_DATA_UncertaintySources_AK4PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To3[isrc] = unc;
    }
    for (int isrc = 0; isrc < nsrc_Eta3To5; isrc++) {
      const char *name = srcnames_Eta3To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/Summer16_23Sep2016V4_DATA_UncertaintySources_AK4PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta3To5[isrc] = unc;
    }
    for (int isrc = 0; isrc < nsrc_RelativeBal; isrc++) {
      const char *name = srcnames_RelativeBal[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/Summer16_23Sep2016V4_DATA_UncertaintySources_AK4PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_RelativeBal[isrc] = unc;
    }

    map< TString , vector<JetCorrectionUncertainty*> > jec_unc_map = {
      { "jecUncEta0To5"     , vsrc_Eta0To5 },
      { "jecUncEta0To3"     , vsrc_Eta0To3 },
      { "jecUncEta3To5"     , vsrc_Eta3To5 },
      { "jecUncRelativeBal" , vsrc_RelativeBal }};

    // END : JEC uncertainties (NEW) ================================================================================================================

    int nTotalFiles = 0;
    std::string dummy;
    // count number of files --->
    while (fileList0 >> dummy) nTotalFiles++;
    std::vector<Period> periods;
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

    // Load CrystalBallEfficiency class
    TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
    int openSuccessful = gSystem->Load( pathToCrystalLib );
    if (openSuccessful !=0 ) {
      cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
      exit( -1 );
    }

    //*****************
    //****** BDT ******
    TH1F * histMva =  new TH1F("MVA_BDT", "MVA_BDT",100 , -1.0, 1.0);
    //TH1F * histMva_massCut =  new TH1F("MVA_BDT_massCut", "MVA_BDT",100 , -1.0, 1.0);
    //This loads the library
    TMVA::Tools::Instance();
    
    //Create TMVA Reader Object
    TMVA::Reader *reader = new TMVA::Reader("!V:!Color");
    //create set of variables as declared in weight file and declared them to reader
    reader->AddVariable( "met", &met );
    reader->AddVariable( "dzeta", &dzeta );
    reader->AddVariable( "dphi_mumet:=AngularDistance(phi_2,metphi)", &dphi_mumet );
    reader->AddVariable("pt_2", &pt_2);
    reader->AddVariable("pt_1", &pt_1);
    reader->AddVariable("dr_tt", &dr_tt);
    //BookMethod
    //  reader->BookMVA("BDT", "TMVA/weights/TMVA_Roch_22032016_BDT.weights.xml");
    reader->BookMVA("BDT", "/nfs/dust/cms/user/rasp/CMSSW/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/data/hMSSM_BDT_GG_BBH_M_500_Set.weights.xml");
    
    // BBH
    TMVA::Reader *readerBBH = new TMVA::Reader("!V:!Color");
    readerBBH->AddVariable( "met", &met );
    readerBBH->AddVariable( "dzeta", &dzeta );
    readerBBH->AddVariable( "dphi_mumet:=AngularDistance(phi_2,metphi)", &dphi_mumet );
    readerBBH->AddVariable("pt_2", &pt_2);
    readerBBH->AddVariable("pt_1", &pt_1);
    readerBBH->AddVariable("dr_tt", &dr_tt);
    readerBBH->BookMVA("BDT", "/nfs/dust/cms/user/rasp/CMSSW/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/data/hMSSM_NoDY_BDT_GG_BBH_M_500_Set.weights.xml");
    
    // GGH
    TMVA::Reader *readerGGH = new TMVA::Reader("!V:!Color");
    readerGGH->AddVariable( "met", &met );
    readerGGH->AddVariable( "dzeta", &dzeta );
    readerGGH->AddVariable( "dphi_mumet:=AngularDistance(phi_2,metphi)", &dphi_mumet );
    readerGGH->AddVariable("pt_2", &pt_2);
    readerGGH->AddVariable("pt_1", &pt_1);
    readerGGH->AddVariable("dr_tt", &dr_tt);
    readerGGH->BookMVA("BDT", "/nfs/dust/cms/user/rasp/CMSSW/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/data/hMSSM_NoDY_BDT_GG_H_M_500_Set.weights.xml");
    
    
    // PU reweighting
    PileUp * PUofficial = new PileUp();
    if (!isData && !isEmbedded){
       TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read");
       TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read");
       TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
       TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
       PUofficial->set_h_data(PU_data);
       PUofficial->set_h_MC(PU_mc);
    }
    // Lepton Scale Factors
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_muon23 = new ScaleFactor();
    SF_muon23->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon23TriggerFile));
    ScaleFactor * SF_muon8 = new ScaleFactor();
    SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));
    ScaleFactor * SF_electronIdIso = new ScaleFactor();
    SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
    ScaleFactor * SF_electron23 = new ScaleFactor();
    SF_electron23->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron23TriggerFile));
    ScaleFactor * SF_electron12 = new ScaleFactor();
    SF_electron12->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron12TriggerFile));
    
    // tracking efficiency SF
    TFile * fileTrackingSF = new TFile(TString(cmsswBase)+"/src/"+TString(trackingSFFile));
    TH1D * trackEffMuonH = (TH1D*)fileTrackingSF->Get("effTrackingMu");
    TH1D * trackEffEleH = (TH1D*)fileTrackingSF->Get("effTrackingE");
    
    // MEt filters          
    std::vector<TString> metFlags; metFlags.clear();
    metFlags.push_back("Flag_HBHENoiseFilter");  
    metFlags.push_back("Flag_HBHENoiseIsoFilter"); 
    metFlags.push_back("Flag_globalTightHalo2016Filter");
    metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
    metFlags.push_back("Flag_goodVertices"); 
    if (isData)
       metFlags.push_back("Flag_eeBadScFilter");
    metFlags.push_back("Flag_BadPFMuonFilter"); 
    metFlags.push_back("Flag_BadChargedCandidateFilter"); 

    std::vector<TString> badChargedCandidateFlag; badChargedCandidateFlag.clear();
    badChargedCandidateFlag.push_back("Flag_BadChargedCandidateFilter");
    std::vector<TString> badPFMuonFlag; badPFMuonFlag.clear();
    badPFMuonFlag.push_back("Flag_BadPFMuonFilter"); 
    std::vector<TString> badMuonFlag; badMuonFlag.clear();
    badMuonFlag.push_back("Flag_badMuons"); 
    std::vector<TString> badGlobalMuonFlag; badGlobalMuonFlag.clear();
    badGlobalMuonFlag.push_back("Flag_BadGlobalMuonFilter"); //not applied
    std::vector<TString> muonBadTrackFlag; muonBadTrackFlag.clear();
    muonBadTrackFlag.push_back("Flag_muonBadTrackFilter"); //not applied
    std::vector<TString> chargedHadronTrackResolutionFlag; chargedHadronTrackResolutionFlag.clear();
    chargedHadronTrackResolutionFlag.push_back("Flag_chargedHadronTrackResolutionFilter");//not applied

    MEtSys metSys(MetSysFileName);
    
    RecoilCorrector recoilMetCorrector(RecoilFileName);
    
    CalibrationOfImpactParameters calibrateIP;

    // SV fit mass   // not needed for ClassicSV FIT
    edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
    TH1::AddDirectory(false);
    TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
    
    // qcd weight (dzeta cut)
    //QCDModelForEMu qcdWeight("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_2016BtoH.root");
    // qcd weight DZeta cut
    //QCDModelForEMu qcdWeightNoDzeta("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_2016BtoH.root");
    
    TString correctionsWorkspaceFileName_qcd = cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/htt_scalefactors_v16_5_2.root";
    TFile * correctionWorkSpaceFile_qcd = new TFile(correctionsWorkspaceFileName_qcd);
    RooWorkspace *correctionWS_qcd = (RooWorkspace*)correctionWorkSpaceFile_qcd->Get("w");
    // BTag scale factors
    BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_Moriond17_B_H.csv");
    BTagCalibrationReader reader_BTAG(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    reader_BTAG.load(calib,BTagEntry::FLAV_B,"comb");
    reader_BTAG.load(calib,BTagEntry::FLAV_C,"comb");
    reader_BTAG.load(calib,BTagEntry::FLAV_UDSG,"incl");
    
    float etaBTAG[2] = {0.5,2.1};
    float ptBTAG[5] = {25.,35.,50.,100.,200.};
    
    std::cout << std::endl;
    for (int iEta=0; iEta<2; ++iEta) {
        for (int iPt=0; iPt<5; ++iPt) {
            float sfB = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
            float sfC = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
            float sfLight = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
            printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
        }
    }
    std::cout << std::endl;
    
    TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_Moriond2017.root"));
    TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
    TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
    TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
    TRandom3 rand;
    
    float MaxBJetPt = 1000.;
    float MaxLJetPt = 1000.;
    float MinLJetPt = 20.;
    float MinBJetPt = 20.; // !!!!!
    
    TFile *file_ggh_reweighting = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/NNLOPS_reweight.root"));
    TGraph * gr_NNLOPSratio_pt_mcatnlo_0jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_mcatnlo_0jet");
    TGraph * gr_NNLOPSratio_pt_mcatnlo_1jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_mcatnlo_1jet");
    TGraph * gr_NNLOPSratio_pt_mcatnlo_2jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_mcatnlo_2jet");
    TGraph * gr_NNLOPSratio_pt_mcatnlo_3jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_mcatnlo_3jet");

    // Z pt mass weights 
    TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName);
    if (fileZMassPtWeights->IsZombie()) {
        std::cout << "File " << TString(cmsswBase) << "/src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
        exit(-1);
    }
    TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName);
    if (histZMassPtWeights==NULL) {
        std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName
        << std::endl;
        exit(-1);
    }
    
    // Z-pt weights from correction WS
    TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/HTT-utilities/CorrectionsWorkspace/htt_scalefactors_v16_5.root";
    TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
    RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");
    //  exit(-1);

    // correction WS
    TString correctionsWorkspaceFileName_embedded = TString(cmsswBase)+"/src/"+correctionWSFile_embedded;
    TFile * correctionWorkSpaceFile_embedded  = new TFile(correctionsWorkspaceFileName_embedded);
    RooWorkspace *correctionWS_embedded  = (RooWorkspace*)correctionWorkSpaceFile_embedded->Get("w");

    TString correctionsWorkspaceFileName_embedded_trigger = TString(cmsswBase)+"/src/"+correctionWSFile_embedded_trigger;
    TFile * correctionWorkSpaceFile_embedded_trigger  = new TFile(correctionsWorkspaceFileName_embedded_trigger);
    RooWorkspace *correctionWS_embedded_trigger  = (RooWorkspace*)correctionWorkSpaceFile_embedded_trigger->Get("w");



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
	// numberOfEntries = 10000;

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
            veto_embedded = false;
            //      bool isPrompMuPlus = false;
            //      bool isPrompMuMinus = false;
            //      bool isPrompElePlus = false;
            //      bool isPrompEleMinus = false;
            
            // weights
            mcweight = analysisTree.genweight;
            if (sync && mcweight<0) continue; // TO DO: needed for sync?
	    weightScale1 = analysisTree.weightScale1;
	    weightScale2 = analysisTree.weightScale2;
	    weightScale3 = analysisTree.weightScale3;
	    weightScale4 = analysisTree.weightScale4;
	    weightScale5 = analysisTree.weightScale5;
	    weightScale6 = analysisTree.weightScale6;
	    weightScale7 = analysisTree.weightScale7;
	    weightScale8 = analysisTree.weightScale8;

	    weightPDFup   = analysisTree.weightPDFup;
	    weightPDFdown = analysisTree.weightPDFdown;

       weight_ggh_NNLOPS = 1.;
            puweight = 1;
            trigweight_1 = 1;
            trigweight_2 = 1;
            trigweight = 1;
            idweight_1 = 1;
            idweight_2 = 1;
            isoweight_1 = 1;
            isoweight_2 = 1;
            effweight = 1;
	    effweight_jetMuUp = 1;
	    effweight_jetMuDown = 1;
	    effweight_jetEUp = 1;
	    effweight_jetEDown = 1;
            fakeweight = 1;
            embeddedWeight = 1;
            signalWeight = 1;
            topptweight = 1;
	    topptweightRun2 = 1;

       THU_ggH_Mu = 1.0;
       THU_ggH_Res = 1.0;
       THU_ggH_Mig01 = 1.0;
       THU_ggH_Mig12 = 1.0;
       THU_ggH_VBF2j = 1.0;
       THU_ggH_VBF3j = 1.0;
       THU_ggH_PT60 = 1.0;
       THU_ggH_PT120 = 1.0;
       THU_ggH_qmtop = 1.0;
            qcdweight = 1;
            qcdweightup = 1;
            qcdweightdown = 1;
            qcdweight_nodzeta = 1;
            qcdweightup_nodzeta = 1;
            qcdweightdown_nodzeta = 1;

            qcdweight_0jet_rate_up =  1;
            qcdweight_0jet_rate_down =  1;
            qcdweight_1jet_rate_up =  1;
            qcdweight_1jet_rate_down =  1;
            qcdweight_0jet_shape_up =  1;
            qcdweight_0jet_shape_down =  1;
            qcdweight_1jet_shape_up =  1;
            qcdweight_1jet_shape_down =  1;
            
            qcdweight_iso_up =  1;
            qcdweight_iso_down =  1;
            zptmassweight = 1;
	    zptmassweight_esup = 1;
	    zptmassweight_esdown = 1;
	    zptmassweight_ttup = 1;
	    zptmassweight_ttdown = 1;
	    zptmassweight_statpt0up = 1;
	    zptmassweight_statpt0down = 1;
	    zptmassweight_statpt40up = 1;
	    zptmassweight_statpt40down = 1;
	    zptmassweight_statpt80up = 1;
	    zptmassweight_statpt80down = 1;

            weight = 1;
            
            nuPx = 0;
            nuPy = 0;
            nuPz = 0;
            nuPt = 0;
            nuPhi = 0;
            
            nuPx_msv = 0;
            nuPy_msv = 0;
            nuPz_msv = 0;
            nuPt_msv = 0;
            nuPhi_msv = 0;
            
            lepPx = 0;
            lepPy = 0;
            lepPz = 0;
            bosonPx = 0;
            bosonPy = 0;
            bosonPz = 0;
            bosonPt = 0;
            bosonMass = -1;
            
            njets_HTXS = -1.;
            higgspt_HTXS = -1.;
            htxs_stage0cat = -1.;
            htxs_stage1cat = -1.;


            metFilters_ = true;

	    badChargedCandidateFilter_ = true;
            badPFMuonFilter_ = true;
	    badGlobalMuonFilter_ = true;
	    muonBadTrackFilter_ = true;
	    chargedHadronTrackResolutionFilter_ = true;

	    badMuonFilter_ = true;
	    duplicateMuonFilter_ = true;

            float topPt = -1;
            float antitopPt = -1;
            
            bool isGSfound = false;
            std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
            std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
            std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
            std::vector<TLorentzVector> promptMuons; promptMuons.clear();
            std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
            std::vector<TLorentzVector> tauNeutrinos; tauNeutrinos.clear();
            
            TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
            TLorentzVector tauNeutrinosLV;  tauNeutrinosLV.SetXYZT(0,0,0,0);
            TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
            TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
            TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);

	    TLorentzVector genBosonLV; genBosonLV.SetXYZT(0,0,0,0);
	    TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);
    
            if (!isData) {
                
	      // computing boson 4-vector

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
                    
                    if (analysisTree.genparticles_pdgid[igen]==22 && analysisTree.genparticles_status[igen]==44)
                        isGSfound = true;
                    
                    if (analysisTree.genparticles_pdgid[igen]==23) {
                        zBosonLV = genLV;
                    }
                    if (analysisTree.genparticles_pdgid[igen]==25||
                        analysisTree.genparticles_pdgid[igen]==35||
                        analysisTree.genparticles_pdgid[igen]==36) {
                        hBosonLV = genLV;
                    }
                    if (abs(analysisTree.genparticles_pdgid[igen])==24) {
                        wBosonLV = genLV;
                    }

		    bool fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1;
		    bool isMuon = false;
		    bool isElectron = false;
		    bool isNeutrino = false;
		    bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==11) {
		      isElectron = true;
                        if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                            promptElectrons.push_back(genLV);
                            promptElectronsLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                    }
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==13) {
		      isMuon = true;
                        if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                            promptMuons.push_back(genLV);
                            promptMuonsLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                    }
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==12||
                        abs(analysisTree.genparticles_pdgid[igen])==14||
                        abs(analysisTree.genparticles_pdgid[igen])==16)  {
		      isNeutrino = true;
                        if ((analysisTree.genparticles_fromHardProcess[igen]||analysisTree.genparticles_isPrompt[igen])&&
                            !analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
                            analysisTree.genparticles_status[igen]==1) {
                            promptNeutrinos.push_back(genLV);
                            promptNeutrinosLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                        if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
                            analysisTree.genparticles_status[igen]==1) {
                            tauNeutrinos.push_back(genLV);
                            tauNeutrinosLV += genLV;
                        }
                    }
                    
		    bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
		    bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);
		    
		    if (isBoson)
		      genBosonLV += genLV;
		    if (isVisibleBoson)
		      genVisBosonLV += genLV;

                }
                
                if (isGSfound) {
                    //	  std::cout << "gamma* found : " << std::endl;
                    if (removeGammaStar) continue;
                }
                
                if (isDY||isEmbedded) {
                   
                    if (promptTausFirstCopy.size()==2) {
                        isZTT = true; isZMM = false; isZEE = false;
                        bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
                        bosonMass = promptTausLV.M();
                        lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
                        mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);
                        
                        double gt1_pt  = promptTausFirstCopy[0].Pt();
                        double gt1_eta = promptTausFirstCopy[0].Eta();
                        double gt2_pt  = promptTausFirstCopy[1].Pt();
                        double gt2_eta = promptTausFirstCopy[1].Eta();
                        correctionWS_embedded->var("gt_pt")->setVal(gt1_pt);
                        correctionWS_embedded->var("gt_eta")->setVal(gt1_eta);
                        double id1_embed = correctionWS_embedded->function("m_sel_idEmb_ratio")->getVal();
                        correctionWS_embedded->var("gt_pt")->setVal(gt2_pt);
                        correctionWS_embedded->var("gt_eta")->setVal(gt2_eta);
                        double id2_embed = correctionWS_embedded->function("m_sel_idEmb_ratio")->getVal();
                        correctionWS_embedded->var("gt1_pt")->setVal(gt1_pt);
                        correctionWS_embedded->var("gt1_eta")->setVal(gt1_eta);
                        correctionWS_embedded->var("gt2_pt")->setVal(gt2_pt);
                        correctionWS_embedded->var("gt2_eta")->setVal(gt2_eta);
                        double trg_embed = correctionWS_embedded->function("m_sel_trg_ratio")->getVal();
                        embeddedWeight = id1_embed * id2_embed * trg_embed;
                        
                    }
                    else if (promptMuons.size()==2) {
                        isZTT = false; isZMM = true; isZEE = false;
                        bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
                        bosonMass = promptMuonsLV.M();
                        lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
                        mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
                    }
                    else {
                        isZTT = false; isZMM = false; isZEE = true;
                        bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
                        bosonMass = promptElectronsLV.M();
                        lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
                        if (promptElectrons.size()==2)
                            mtBoson_gen = mT(promptElectrons[0],promptElectrons[1]);
                    }
                    nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
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
                    nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
                }
                else {
                    TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
                    bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
                    TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
                    lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
                    nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
                }
            
		bosonPx = genBosonLV.Px();
		bosonPy = genBosonLV.Py();
		bosonPz = genBosonLV.Pz();
		bosonPt = genBosonLV.Pt();
		bosonMass = genBosonLV.M();

		lepPx = genVisBosonLV.Px();
		lepPy = genVisBosonLV.Py();
		lepPz = genVisBosonLV.Pz();

                nuPt = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
                nuPhi = TMath::ATan2(nuPy,nuPx);
                
                
                
                if (isDY) { // applying Z pt mass weights
                    zptmassweight = 1;
                    if (bosonMass>50.0) {
                        float bosonMassX = bosonMass;
                        float bosonPtX = bosonPt;
                        if (bosonMassX>1000.) bosonMassX = 1000.;
                        if (bosonPtX<1.)      bosonPtX = 1.;
                        if (bosonPtX>1000.)   bosonPtX = 1000.;
			zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
			histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
		    }
                }
                
                ALL->Fill(0.0);
                if (isZMM) ZMM->Fill(0.);
                if (isZEE) ZEE->Fill(0.);
                if (isZTT) ZTT->Fill(0.);
                isZLL = isZMM || isZEE;
            }

            if (isSignal){
               njets_HTXS = analysisTree.htxs_njets30;
               higgspt_HTXS = analysisTree.htxs_higgsPt;
               htxs_stage0cat = analysisTree.htxs_stage0cat;
               htxs_stage1cat = analysisTree.htxs_stage1cat;
               if (apply_ggh_reweighting)
                  {
                     if      (njets_HTXS==0) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_mcatnlo_0jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)125.0));
                     else if (njets_HTXS==1) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_mcatnlo_1jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)625.0));
                     else if (njets_HTXS==2) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_mcatnlo_2jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)800.0));
                     else if (njets_HTXS>=3) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_mcatnlo_3jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)925.0));
                     else weight_ggh_NNLOPS = 1.0;
                     std::cout<<weight_ggh_NNLOPS<<std::endl;
                     std::vector<double> ggF_unc = qcd_ggF_uncertSF_2017(njets_HTXS, higgspt_HTXS, htxs_stage1cat, 1.0);
                     THU_ggH_Mu = ggF_unc[0];
                     THU_ggH_Res = ggF_unc[1];
                     THU_ggH_Mig01 = ggF_unc[2];
                     THU_ggH_Mig12 = ggF_unc[3];
                     THU_ggH_VBF2j = ggF_unc[4];
                     THU_ggH_VBF3j = ggF_unc[5];
                     THU_ggH_PT60 = ggF_unc[6];
                     THU_ggH_PT120 = ggF_unc[7];
                     THU_ggH_qmtop = ggF_unc[8];
                     std::cout<<"THU 1: "<<ggF_unc[0]<<std::endl;
                     std::cout<<"THU 2: "<<ggF_unc[1]<<std::endl;
                     std::cout<<"THU 3: "<<ggF_unc[2]<<std::endl;
                  }
            }

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
            
            npv = analysisTree.primvertex_count;
            npu = analysisTree.numtruepileupinteractions;
            rho = analysisTree.rho;
            
            npartons = analysisTree.genparticles_noutgoing;
            
            if (!isData && !isEmbedded) {
                puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));

                if (topPt>0&&antitopPt>0) {
		  topptweight = topPtWeight(topPt,antitopPt,true);
		  topptweightRun2 = topPtWeight(topPt,antitopPt,false);
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

	    unsigned int nMu23Ele12DzFilter = 0;
	    bool isMu23Ele12DzFilter = false;

	    unsigned int nMu8Ele23DzFilter = 0;
            bool isMu8Ele23DzFilter = false;
            
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
		if (HLTFilter==Mu23Ele12DzFilter) {
		  nMu23Ele12DzFilter = i;
		  isMu23Ele12DzFilter = true;
		}
		if (HLTFilter==Mu8Ele23DzFilter) {
		  nMu8Ele23DzFilter = i;
		  isMu8Ele23DzFilter = true;
		}
            }
            if (applyTriggerMatch) {
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
		if (applyDzFilterMatch) {
		  if (!isMu23Ele12DzFilter) {
		    std::cout << "HLT filter " << Mu23Ele12DzFilter << " not found" << std::endl;
                    exit(-1);
		  }
		  if (!isMu8Ele23DzFilter) {
		    std::cout << "HLT filter " << Mu8Ele23DzFilter << " not found" << std::endl;
                    exit(-1);
		  }
		}
            }
            
            unsigned int nBTagDiscriminant = 0;
            for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
                TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
                if (discr.Contains(BTagDiscriminator))
                    nBTagDiscriminant = iBTag;
            }
            
            
            // MET Filters                    
	    metFilters_ = metFiltersPasses(analysisTree,metFlags);
	    badChargedCandidateFilter_ = metFiltersPasses(analysisTree,badChargedCandidateFlag);
	    badPFMuonFilter_ = metFiltersPasses(analysisTree,badPFMuonFlag);
       badMuonFilter_ = metFiltersPasses(analysisTree,badMuonFlag);
	    badGlobalMuonFilter_ = metFiltersPasses(analysisTree,badGlobalMuonFlag);
	    muonBadTrackFilter_ = metFiltersPasses(analysisTree,muonBadTrackFlag);
	    chargedHadronTrackResolutionFilter_ = metFiltersPasses(analysisTree,chargedHadronTrackResolutionFlag);


            // electron selection
            vector<int> electrons; electrons.clear();
            for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
               if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;  
               if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
               if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
               if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
               bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie];
               if (applySpring16ElectronId) electronMvaId = analysisTree.electron_mva_wp80_general_Spring16_v1[ie]>0.5; 
               if (!electronMvaId) continue;
               if (!analysisTree.electron_pass_conversion[ie]) continue;         
               if (analysisTree.electron_nmissinginnerhits[ie]>1) continue;      
               electrons.push_back(ie);
            }
        
            // muon selection
            vector<int> muons; muons.clear();
            for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
                if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;          
                if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
                if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
                if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
                if (analysisTree.muon_isBad[im]) badMuonFilter_ = false;                     
                if (analysisTree.muon_isDuplicate[im]) duplicateMuonFilter_ = false;                   
		bool muonId = analysisTree.muon_isMedium[im]; 
		if (applyICHEPMuonId) muonId = analysisTree.muon_isICHEP[im];
                if (!muonId) continue;
                muons.push_back(im);
            }
       if (electrons.size()==0) continue; 
       if (muons.size()==0) continue;          
       
            
            // selecting muon and electron pair (OS or SS);
            int electronIndex = -1;
            int muonIndex = -1;
            
            float isoMuMin = 1e+10;
            float isoEleMin = 1e+10;
            bool isMuon23matched = false;
            bool isMuon8matched  = false;
            bool isElectron23matched = false;
            bool isElectron12matched = false;
            for (unsigned int im=0; im<muons.size(); ++im) {
                unsigned int mIndex  = muons.at(im);
                float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
                float photonIsoMu = analysisTree.muon_photonIso[mIndex];
                float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
                float puIsoMu = analysisTree.muon_puIso[mIndex];
                if (isMuonIsoR03) {
                    neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
                    photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
                    chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
                    puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
                }
                float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
                neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
                float absIsoMu = chargedHadIsoMu + neutralIsoMu;
                float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];

                for (unsigned int ie=0; ie<electrons.size(); ++ie) {
                    
                    unsigned int eIndex = electrons.at(ie);
                    
                    float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
                                      analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);
                    
                    if (dR<dRleptonsCut) continue;

                    bool trigMatch =
                    (analysisTree.muon_pt[mIndex]>ptMuonHighCut) ||
                    (analysisTree.electron_pt[eIndex]>ptElectronHighCut);
                    
                    if (!trigMatch) continue;

                    float neutralHadIsoEle = analysisTree.electron_neutralHadIso[eIndex];
                    float photonIsoEle = analysisTree.electron_photonIso[eIndex];
                    float chargedHadIsoEle = analysisTree.electron_chargedHadIso[eIndex];
                    float puIsoEle = analysisTree.electron_puIso[eIndex];
                    if (isElectronIsoR03) {
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
			if (analysisTree.electron_pt[eIndex]>analysisTree.electron_pt[electronIndex]) {
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

            if (electronIndex<0) continue;
            if (muonIndex<0) continue;
       
            metFilters_ = metFilters_ && badMuonFilter_ && duplicateMuonFilter_ &&  badChargedCandidateFilter_ && badPFMuonFilter_;              
            if (badMuonFilter_ < 0.5) continue;
            if (duplicateMuonFilter_ < 0.5) continue; 

	    bool isMu23 = false;
	    bool isMu8 = false;
	    bool isMu23dz = false;
	    bool isMu8dz  = false;

	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      float dRtrig = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]&&
		    analysisTree.muon_pt[muonIndex]>ptMuonHighCut) { // Mu23 Leg
		  isMu23 = true;
		}
		if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]&&
		    analysisTree.muon_pt[muonIndex]>ptMuonLowCut) { // Mu8 Leg
		  isMu8 = true;
		}
		if (analysisTree.trigobject_filters[iT][nMu23Ele12DzFilter])
		  isMu23dz = true;
		if (analysisTree.trigobject_filters[iT][nMu8Ele23DzFilter])
		  isMu8dz = true;
		
	      }
	    }
	    if (applyDzFilterMatch) {
	      isMu23 = isMu23 && isMu23dz;
	      isMu8 = isMu8 && isMu8dz;
	    }

	    bool isEle23 = false;
	    bool isEle12 = false;
	    bool isEle23dz = false;
	    bool isEle12dz = false;
	    
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      float dRtrig = deltaR(analysisTree.electron_eta[electronIndex],analysisTree.electron_phi[electronIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<deltaRTrigMatch) {
		if (analysisTree.trigobject_filters[iT][nHighPtLegElectron]&&
		    analysisTree.electron_pt[electronIndex]>ptElectronHighCut) { // Ele23 Leg
		  isEle23 = true;
		}
		if (analysisTree.trigobject_filters[iT][nLowPtLegElectron]&&
		    analysisTree.electron_pt[electronIndex]>ptElectronLowCut) { // Ele12 Leg
		  isEle12 = true;
		}
		if (analysisTree.trigobject_filters[iT][nMu23Ele12DzFilter])
		  isEle12dz = true;
		if (analysisTree.trigobject_filters[iT][nMu8Ele23DzFilter])
		  isEle23dz = true;
	      }
	    }
            
	    if (applyDzFilterMatch) {
	      isEle23 = isEle23 && isEle23dz;
	      isEle12 = isEle12 && isEle12dz;
	    }		      
            

	    trg_muonelectron = 
	      (isMu23&&isEle12&&analysisTree.muon_pt[muonIndex]>ptMuonHighCut) ||
	      (isMu8&&isEle23&&analysisTree.electron_pt[electronIndex]>ptElectronHighCut);

	    if (!sync && applyTriggerMatch&&trg_muonelectron<0.5) continue;
            
            os = (analysisTree.muon_charge[muonIndex]*analysisTree.electron_charge[electronIndex]) < 0;
            // looking for extra electron  
            bool foundExtraElectron = false;
            for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
               if (int(ie)==electronIndex) continue;
               if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
               if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
               if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
               if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue; 
               bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
               if (applySpring16ElectronId) electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie]>0.5;
               if (!electronMvaId&&applyVetoElectronId) continue;
               if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
               
               if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
               float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
               float photonIsoEle = analysisTree.electron_photonIso[ie];
               float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
               float puIsoEle = analysisTree.electron_puIso[ie];
               if (isElectronIsoR03) {                                               
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
                bool muonId = analysisTree.muon_isMedium[im];
                if (applyICHEPMuonId) muonId = analysisTree.muon_isICHEP[im];
                if (!muonId&&applyVetoMuonId) continue;
                float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
                float photonIsoMu = analysisTree.muon_photonIso[im];
                float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
                float puIsoMu = analysisTree.muon_puIso[im];
                if (isMuonIsoR03) {
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
            
            // filling muon variables
            pt_2 = analysisTree.muon_pt[muonIndex];
            eta_2 = analysisTree.muon_eta[muonIndex];
            phi_2 = analysisTree.muon_phi[muonIndex];
            q_2 = -1;
            if (analysisTree.muon_charge[muonIndex]>0)
                q_2 = 1;
            mva_2 = -10;
            d0_2 = analysisTree.muon_dxy[muonIndex];
            dZ_2 = analysisTree.muon_dz[muonIndex];
            iso_2 = isoMuMin;
            m_2 =  classic_svFit::muonMass;
            
            float eleScale = eleScaleBarrel;
            if (fabs(analysisTree.electron_eta[electronIndex])>1.479) eleScale = eleScaleEndcap;
            // filling electron variables
            pt_1 = analysisTree.electron_pt[electronIndex];
            
            eta_1 = analysisTree.electron_eta[electronIndex];
            phi_1 = analysisTree.electron_phi[electronIndex];
            q_1 = -1;
            if (analysisTree.electron_charge[electronIndex]>0)
                q_1 = 1;
            mva_1 = analysisTree.electron_mva_id_nontrigPhys14[electronIndex];
            d0_1 = analysisTree.electron_dxy[electronIndex];
            dZ_1 = analysisTree.electron_dz[electronIndex];
            iso_1 = isoEleMin;
            m_1 =  classic_svFit::electronMass;
            
            
            isoweight_1 = 1;
            isoweight_2 = 1;
            trigweight_1 = 1;
            trigweight_2 = 1;
            trigweight = 1;
            effweight = 1;
            
            
            if (!isData || isEmbedded) {
           
               correctionWS_embedded->var("e_pt")->setVal(pt_1);
               correctionWS_embedded->var("e_eta")->setVal(eta_1);
               correctionWS_embedded->var("e_iso")->setVal(iso_1);
               correctionWS_embedded->var("m_pt")->setVal(pt_2);
               correctionWS_embedded->var("m_eta")->setVal(eta_2);
               correctionWS_embedded->var("m_iso")->setVal(iso_2);
                // scale factors
               if (isEmbedded) {
                  isoweight_1 = correctionWS_embedded->function("e_looseiso_ratio")->getVal() * correctionWS_embedded->function("e_id_ratio")->getVal();
                  isoweight_2 = correctionWS_embedded->function("m_looseiso_ratio")->getVal() * correctionWS_embedded->function("m_id_ratio")->getVal();
               }
               else {
                  isoweight_1 = (float)SF_electronIdIso->get_ScaleFactor(double(pt_1),double(eta_1));
                  isoweight_2 = (float)SF_muonIdIso->get_ScaleFactor(double(pt_2),double(eta_2));
               }
		
                correctionWS->var("e_pt")->setVal(pt_1);
                correctionWS->var("e_eta")->setVal(eta_1);
                idweight_1 = correctionWS->function("e_trk_ratio")->getVal();
		          
                correctionWS->var("m_eta")->setVal(eta_2);
                correctionWS->var("m_pt")->setVal(pt_2);
                idweight_2 = correctionWS->function("m_trk_ratio")->getVal();
                isoweight_1 *= idweight_1;
                isoweight_2 *= idweight_2;
                
                //		cout << "isoweight_1 = " << isoweight_1
                //		     << "isoweight_2 = " << isoweight_2 << endl;
                
                float Ele23EffData = 1;
                float Ele12EffData = 1;
                float Mu23EffData = 1;
                float Mu8EffData = 1;
                
                float Ele23EffMC = 1;
                float Ele12EffMC = 1;
                float Mu23EffMC = 1;
                float Mu8EffMC = 1;
                
                Ele23EffData = (float)SF_electron23->get_EfficiencyData(double(pt_1),double(eta_1));
                Ele12EffData = (float)SF_electron12->get_EfficiencyData(double(pt_1),double(eta_1));
                Mu23EffData = (float)SF_muon23->get_EfficiencyData(double(pt_2),double(eta_2));
                Mu8EffData = (float)SF_muon8->get_EfficiencyData(double(pt_2),double(eta_2));
              
                correctionWS_embedded_trigger->var("e_pt")->setVal(pt_1);
                correctionWS_embedded_trigger->var("e_eta")->setVal(eta_1);
                correctionWS_embedded_trigger->var("m_pt")->setVal(pt_2);
                correctionWS_embedded_trigger->var("m_eta")->setVal(eta_2);
                
                if (isEmbedded){
                   Ele23EffData = correctionWS_embedded_trigger->function("e_trg23_binned_ic_data")->getVal();
                   Ele12EffData = correctionWS_embedded_trigger->function("e_trg12_binned_ic_data")->getVal(); 
                   Mu23EffData  = correctionWS_embedded_trigger->function("m_trg23_binned_ic_data")->getVal();
                   Mu8EffData   = correctionWS_embedded_trigger->function("m_trg8_binned_ic_data")->getVal();
                }
                float trigWeightData = Mu23EffData*Ele12EffData + Mu8EffData*Ele23EffData - Mu23EffData*Ele23EffData;
                
                if (applyTriggerMatch && !isData) {
                   if (!isEmbedded){
                      Ele23EffMC   = (float)SF_electron23->get_EfficiencyMC(double(pt_1),double(eta_1));
                      Ele12EffMC   = (float)SF_electron12->get_EfficiencyMC(double(pt_1),double(eta_1));
                      Mu23EffMC   = (float)SF_muon23->get_EfficiencyMC(double(pt_2),double(eta_2));
                      Mu8EffMC   = (float)SF_muon8->get_EfficiencyMC(double(pt_2),double(eta_2));
                   }
                   else {
                      Ele23EffMC = correctionWS_embedded_trigger->function("e_trg23_binned_ic_embed")->getVal();
                      Ele12EffMC = correctionWS_embedded_trigger->function("e_trg12_binned_ic_embed")->getVal(); 
                      Mu23EffMC  = correctionWS_embedded_trigger->function("m_trg23_binned_ic_embed")->getVal();
                      Mu8EffMC   = correctionWS_embedded_trigger->function("m_trg8_binned_ic_embed")->getVal();
                   }
                   float trigWeightMC   = Mu23EffMC*Ele12EffMC     + Mu8EffMC*Ele23EffMC     - Mu23EffMC*Ele23EffMC;
           
                   // TO DO: SET ALSO FOR EMBEDDED SAMPLES
                    // if (isMuon23matched && isElectron12matched) {
                    //    if (!isEmbedded){
                    //       trigweight_1 = (float)SF_electron12->get_ScaleFactor(double(pt_1),double(eta_1));
                    //       trigweight_2 = (float)SF_muon23->get_ScaleFactor(double(pt_2),double(eta_2));
                    //    }
                    // }
                    // else if (isMuon8matched && isElectron23matched) {
                    //    if (!isEmbedded){
                    //       trigweight_1 = (float)SF_electron23->get_ScaleFactor(double(pt_1),double(eta_1));
                    //       trigweight_2 = (float)SF_muon8->get_ScaleFactor(double(pt_2),double(eta_2));
                    //    }     
                    // }
                    if (trigWeightMC>1e-6)
                        trigweight = trigWeightData / trigWeightMC;
                    
                }
                else {
                    trigweight = trigWeightData;
                    // if (isMuon23matched && isElectron12matched) {
                    //     trigweight_1 = (float)SF_electron12->get_EfficiencyData(double(pt_1),double(eta_1));
                    //     trigweight_2 = (float)SF_muon23->get_EfficiencyData(double(pt_2),double(eta_2));
                    // }
                    // else if (isMuon8matched && isElectron23matched) {
                    //     trigweight_1 = (float)SF_electron23->get_EfficiencyData(double(pt_1),double(eta_1));
                    //     trigweight_2 = (float)SF_muon8->get_EfficiencyData(double(pt_2),double(eta_2));
                    // }
                }
                
                effweight = trigweight;
            }
            
            // dilepton system
            TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
                                                  analysisTree.muon_py[muonIndex],
                                                  analysisTree.muon_pz[muonIndex],
                                                   classic_svFit::muonMass);
            
            TLorentzVector muonUpLV; muonUpLV.SetXYZM((1.0+muonScale)*analysisTree.muon_px[muonIndex],
                                                      (1.0+muonScale)*analysisTree.muon_py[muonIndex],
                                                      (1.0+muonScale)*analysisTree.muon_pz[muonIndex],
                                                       classic_svFit::muonMass);
            
            TLorentzVector muonDownLV; muonDownLV.SetXYZM((1.0-muonScale)*analysisTree.muon_px[muonIndex],
                                                          (1.0-muonScale)*analysisTree.muon_py[muonIndex],
                                                          (1.0-muonScale)*analysisTree.muon_pz[muonIndex],
                                                           classic_svFit::muonMass);
            
            TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
                                                          analysisTree.electron_py[electronIndex],
                                                          analysisTree.electron_pz[electronIndex],
                                                           classic_svFit::electronMass);
            
            TLorentzVector electronUpLV; electronUpLV.SetXYZM((1.0+eleScale)*analysisTree.electron_px[electronIndex],
                                                              (1.0+eleScale)*analysisTree.electron_py[electronIndex],
                                                              (1.0+eleScale)*analysisTree.electron_pz[electronIndex],
                                                               classic_svFit::electronMass);
            
            TLorentzVector electronDownLV; electronDownLV.SetXYZM((1.0-eleScale)*analysisTree.electron_px[electronIndex],
                                                                  (1.0-eleScale)*analysisTree.electron_py[electronIndex],
                                                                  (1.0-eleScale)*analysisTree.electron_pz[electronIndex],
                                                                   classic_svFit::electronMass);
            
            TLorentzVector dileptonLV = muonLV + electronLV;
            
            m_vis = dileptonLV.M();
            pt_vis = dileptonLV.Pt();
            
            dphi_tt = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
                                 electronLV.Px(),electronLV.Py());
            
            dr_tt = deltaR(muonLV.Eta(),muonLV.Phi(),
                           electronLV.Eta(),electronLV.Phi());
            
            // qcd scale factor
            // no dzeta cut
            //qcdweight     = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            //std::cout<<"qcdweight "<<qcdweight<<std::endl;
            //qcdweightup   = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            //std::cout<<"qcdweight Up"<<qcdweight<<std::endl;
            //qcdweightdown = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            //std::cout<<"qcdweight down"<<qcdweight<<std::endl;
            // dzeta cut
            //qcdweight_nodzeta     = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
            //qcdweightup_nodzeta   = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
            //qcdweightdown_nodzeta = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
         
            
            // counting jets
            vector<unsigned int> jets; jets.clear();
            vector<unsigned int> jetspt20; jetspt20.clear();
            vector<unsigned int> bjets; bjets.clear();
            vector<unsigned int> bjets_mistagUp; bjets_mistagUp.clear();
            vector<unsigned int> bjets_mistagDown; bjets_mistagDown.clear();
            vector<unsigned int> bjets_btagUp; bjets_btagUp.clear();
            vector<unsigned int> bjets_btagDown; bjets_btagDown.clear();
            vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
            vector<unsigned int> bjetsRaw; bjetsRaw.clear();

	    njets = 0;
	    njets_jecUncEta0To5Up   = 0;
	    njets_jecUncEta0To5Down = 0;
	    njets_jecUncEta0To3Up   = 0;
	    njets_jecUncEta0To3Down = 0;
	    njets_jecUncEta3To5Up   = 0;
	    njets_jecUncEta3To5Down = 0;
	    njets_jecUncRelativeBalUp   = 0;
	    njets_jecUncRelativeBalDown = 0;

	    TLorentzVector jetLV;
            TLorentzVector jet1;
            TLorentzVector jet2;
	    map<TString,TLorentzVector> jet1LV_jecUnc;
	    map<TString,TLorentzVector> jet2LV_jecUnc;
	    map<TString,TLorentzVector> metLV_jecUnc;
	    TLorentzVector metLV;
            float met_x = analysisTree.pfmetcorr_ex;
            float met_y = analysisTree.pfmetcorr_ey;
            met = TMath::Sqrt(met_x*met_x + met_y*met_y);
	    metLV.SetXYZT(met_x,met_y,0.,met);

            int indexLeadingJet = -1;
            float ptLeadingJet = -1;
            
            int indexSubLeadingJet = -1;
            float ptSubLeadingJet = -1;
            
            int indexLeadingBJet = -1;
            float ptLeadingBJet = -1;
            
            for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
                
               double absJetEta = fabs(analysisTree.pfjet_eta[jet]);
               double jetEta = analysisTree.pfjet_eta[jet];
               if (absJetEta>jetEtaCut) continue;
                
               float jetPt = analysisTree.pfjet_pt[jet];

	       jetLV.SetPxPyPzE(analysisTree.pfjet_px[jet],
				analysisTree.pfjet_py[jet],
				analysisTree.pfjet_pz[jet],
				analysisTree.pfjet_e[jet]);

	       map<TString,TLorentzVector> jetLV_jecUnc;

	       // Include variations for jec uncertainties
	       for (auto uncer_split : jec_unc_map) {
             float sum_unc   = 0;
             for (auto single_jec_unc : uncer_split.second){
                JetCorrectionUncertainty *unc = single_jec_unc;
                unc->setJetPt(jetPt);
                unc->setJetEta(jetEta);
                double unc_ = unc->getUncertainty(true);
                sum_unc  += pow(unc_,2);
             }
             float unc_total  = TMath::Sqrt(sum_unc);
             jetLV_jecUnc[uncer_split.first + "Up"]   = jetLV * ( 1 + unc_total);
             jetLV_jecUnc[uncer_split.first + "Down"] = jetLV * ( 1 - unc_total);
             // Propagate jec uncertainties to met
             if( metLV_jecUnc.find(uncer_split.first+"Up") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Up"] = metLV;
             if( metLV_jecUnc.find(uncer_split.first+"Down") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Down"] = metLV;
             if ( !((isW||isDY||isSignal)&&!isData && applySimpleRecoilCorrections)) {
                metLV_jecUnc[uncer_split.first + "Up"]   -= jetLV* unc_total;
                metLV_jecUnc[uncer_split.first + "Down"] += jetLV* unc_total;
             }
	       }

               float jetPt_tocheck = jetPt;
               if (sync) jetPt_tocheck = jetPt;
	       else{
	       	 float jetPtmin = jetPt;
	       	 if(jetLV_jecUnc.at("jecUncEta0To5Down").Pt() < jetPtmin) jetPtmin = jetLV_jecUnc.at("jecUncEta0To5Down").Pt();
	       	 if(jetLV_jecUnc.at("jecUncEta0To3Down").Pt() < jetPtmin) jetPtmin = jetLV_jecUnc.at("jecUncEta0To3Down").Pt();
	       	 if(jetLV_jecUnc.at("jecUncEta3To5Down").Pt() < jetPtmin) jetPtmin = jetLV_jecUnc.at("jecUncEta3To5Down").Pt();
	       	 if(jetLV_jecUnc.at("jecUncRelativeBalDown").Pt() < jetPtmin) jetPtmin = jetLV_jecUnc.at("jecUncRelativeBalDown").Pt();
	       	 jetPt_tocheck = jetPtmin;
	       }
	       if (jetPt_tocheck<jetPtLowCut) continue;   
             
               bool isPFJetId = looseJetiD(analysisTree,int(jet));
               if (!isPFJetId) continue;

               bool cleanedJet = true;
               
               float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                  eta_1,phi_1);
               if (dR1<dRJetLeptonCut) cleanedJet = false;
               
               float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                  eta_2,phi_2);
               if (dR2<dRJetLeptonCut) cleanedJet = false;
               
               // jetId
               
               if (!cleanedJet) continue;
              
               if (jetPt>jetPtLowCut)
                  jetspt20.push_back(jet);

               if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance
                  
                  bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                  bool tagged_mistagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                  bool tagged_mistagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                  bool tagged_btagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                  bool tagged_btagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                  bool taggedRaw = tagged;
                  
                  if (!isData) {
                     int flavor = abs(analysisTree.pfjet_flavour[jet]);
                     
                     double jet_scalefactor      = 1;
                     double jet_scalefactor_up   = 1;
                     double jet_scalefactor_down = 1;
                     double JetPtForBTag = jetPt;
                     double tageff = 1;
                     
                     if (flavor==5) {
		       if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
		       if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
		       jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
		       jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
		       jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
		       tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                     }
                     else if (flavor==4) {
		       if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
		       if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
		       jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
		       jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
		       jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
		       tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                     }
                     else {
		       if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
		       if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
		       jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
		       jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
		       jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
		       tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                     }
                     
                     if (tageff<1e-5)      tageff = 1e-5;
                     if (tageff>0.99999)   tageff = 0.99999;
                     rand.SetSeed((int)((jetEta+5)*100000));
                     double rannum = rand.Rndm();

                     if (tagged) { // demote
		       if(jet_scalefactor<1){
			 double fraction = 1-jet_scalefactor;
			 if (rannum<fraction) tagged = false;
		       }
		       if(jet_scalefactor_up<1){
			 double fraction_up = 1-jet_scalefactor_up;
			 if (rannum<fraction_up) tagged_mistagUp = false;
		       }
		       if(jet_scalefactor_down<1){
			 double fraction_down = 1-jet_scalefactor_down;
			 if (rannum<fraction_down) tagged_mistagDown = false;
		       }
		       tagged_btagUp   = tagged;
		       tagged_btagDown = tagged;
                     }
                     else if (!tagged) { // promote
		       if(jet_scalefactor>1){
			 double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
			 if (rannum<fraction) tagged = true;
		       }
		       if(jet_scalefactor_up>1){
			 double fraction_up = (jet_scalefactor_up-1.0)/(1.0/tageff-1.0);
			 if (rannum<fraction_up) tagged_btagUp = true;
		       }
		       if(jet_scalefactor_down>1){
			 double fraction_down = (jet_scalefactor_down-1.0)/(1.0/tageff-1.0);
			 if (rannum<fraction_down) tagged_btagDown = true;
		       }
		       tagged_mistagUp   = tagged;
		       tagged_mistagDown = tagged;
                     }
                  }

                  if (taggedRaw) bjetsRaw.push_back(jet);
                  
                  if (tagged) {
                     bjets.push_back(jet);
                     if (jetPt>ptLeadingBJet) {
                        ptLeadingBJet = jetPt;
                        indexLeadingBJet = jet;
                     }
                  }

		  if(tagged_mistagUp)   bjets_mistagUp.push_back(jet);
		  if(tagged_mistagDown) bjets_mistagDown.push_back(jet);
		  if(tagged_btagUp)     bjets_btagUp.push_back(jet);
		  if(tagged_btagDown)   bjets_btagDown.push_back(jet);
               }

	       if (jetLV_jecUnc.at("jecUncEta0To5Up").Pt()>jetPtHighCut) njets_jecUncEta0To5Up += 1;
	       if (jetLV_jecUnc.at("jecUncEta0To5Down").Pt()>jetPtHighCut) njets_jecUncEta0To5Down += 1;
	       if (jetLV_jecUnc.at("jecUncEta0To3Up").Pt()>jetPtHighCut) njets_jecUncEta0To3Up += 1;
	       if (jetLV_jecUnc.at("jecUncEta0To3Down").Pt()>jetPtHighCut) njets_jecUncEta0To3Down += 1;
	       if (jetLV_jecUnc.at("jecUncEta3To5Up").Pt()>jetPtHighCut) njets_jecUncEta3To5Up += 1;
	       if (jetLV_jecUnc.at("jecUncEta3To5Down").Pt()>jetPtHighCut) njets_jecUncEta3To5Down += 1;
	       if (jetLV_jecUnc.at("jecUncRelativeBalUp").Pt()>jetPtHighCut) njets_jecUncRelativeBalUp += 1;
	       if (jetLV_jecUnc.at("jecUncRelativeBalDown").Pt()>jetPtHighCut) njets_jecUncRelativeBalDown += 1;

               if (jetPt>jetPtHighCut)
		 jets.push_back(jet); 
               if (indexLeadingJet>=0) {
                  if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
                     indexSubLeadingJet = jet;
                     ptSubLeadingJet = jetPt;
                  }
               }
               
               if (jetPt>ptLeadingJet) {
                  indexSubLeadingJet = indexLeadingJet;
                  ptSubLeadingJet = ptLeadingJet;
                  indexLeadingJet = jet;
                  ptLeadingJet = jetPt;
               }
            }
            
            njets = jets.size();
            int njetsMax = njets;
            
            njetspt20 = jetspt20.size();
            nbtag = bjets.size();
            nbtag_mistagUp   = bjets_mistagUp.size();
            nbtag_mistagDown = bjets_mistagDown.size();
            nbtag_btagUp   = bjets_btagUp.size();
            nbtag_btagDown = bjets_btagDown.size();
            nbtag_noSF = bjetsRaw.size();
            
            bpt = -10;
            beta = -10;
            bphi = -10;
            
            if (indexLeadingBJet>=0) {
                bpt = analysisTree.pfjet_pt[indexLeadingBJet];
                beta = analysisTree.pfjet_eta[indexLeadingBJet];
                bphi = analysisTree.pfjet_phi[indexLeadingBJet];
            }
            
            jpt_1 = -10;

            jeta_1 = -10;
            jphi_1 = -10;
            jptraw_1 = -10;
            jptunc_1 = -10;
            jmva_1 = -10;
            jlrm_1 = -10;
            jctm_1 = -10;
            
            if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
                cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
            
            if (indexLeadingJet>=0) {
                jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
                jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
                jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
                jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
                jmva_1 = analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
            }
            
            jpt_2 = -10;
            jeta_2 = -10;
            jphi_2 = -10;
            jptraw_2 = -10;
            jptunc_2 = -10;
            jmva_2 = -10;
            jlrm_2 = -10;
            jctm_2 = -10;
            
            if (indexSubLeadingJet>=0) {
                jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
                jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
                jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
                jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];
                jmva_2 = analysisTree.pfjet_pu_jet_full_mva[indexSubLeadingJet];
            }
            
            mjj =  -10;
            dijetpt = -10;
            dijetphi = -10;
            jdeta =  -10;
            njetingap = 0;
            if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {

	      jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
			      analysisTree.pfjet_py[indexLeadingJet],
			      analysisTree.pfjet_pz[indexLeadingJet],
			      analysisTree.pfjet_e[indexLeadingJet]);

	      for (auto uncer_split : jec_unc_map) {
		float sum_unc   = 0;
		for (auto single_jec_unc : uncer_split.second){
		  JetCorrectionUncertainty *unc = single_jec_unc;
		  unc->setJetPt(jet1.Pt());
		  unc->setJetEta(jet1.Eta());
		  double unc_ = unc->getUncertainty(true);
		  sum_unc  += pow(unc_,2);
		}
		float unc_total = TMath::Sqrt(sum_unc);
		jet1LV_jecUnc[uncer_split.first+"Up"]   = jet1*(1+unc_total);
		jet1LV_jecUnc[uncer_split.first+"Down"] = jet1*(1-unc_total);
	      }

	      jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
			      analysisTree.pfjet_py[indexSubLeadingJet],
			      analysisTree.pfjet_pz[indexSubLeadingJet],
			      analysisTree.pfjet_e[indexSubLeadingJet]);

	      for (auto uncer_split : jec_unc_map) {
		float sum_unc   = 0;
		for (auto single_jec_unc : uncer_split.second){
		  JetCorrectionUncertainty *unc = single_jec_unc;
		  unc->setJetPt(jet2.Pt());
		  unc->setJetEta(jet2.Eta());
		  double unc_ = unc->getUncertainty(true);
		  sum_unc  += pow(unc_,2);
		}
		float unc_total = TMath::Sqrt(sum_unc);
		jet2LV_jecUnc[uncer_split.first+"Up"]   = jet2*(1+unc_total);
		jet2LV_jecUnc[uncer_split.first+"Down"] = jet2*(1-unc_total);
	      }

	      mjj      = (jet1+jet2).M();
	      dijetpt  = (jet1+jet2).Pt();
	      dijetphi = (jet1+jet2).Phi();
	      jdeta = fabs(analysisTree.pfjet_eta[indexLeadingJet]-analysisTree.pfjet_eta[indexSubLeadingJet]);

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
            float met_x_recoilscaleUp = analysisTree.pfmetcorr_ex;
            float met_x_recoilscaleDown = analysisTree.pfmetcorr_ex;
            float met_y_recoilscaleUp = analysisTree.pfmetcorr_ey;
            float met_y_recoilscaleDown = analysisTree.pfmetcorr_ey;
            float met_x_recoilresoUp = analysisTree.pfmetcorr_ex;
            float met_x_recoilresoDown = analysisTree.pfmetcorr_ex;
            float met_y_recoilresoUp = analysisTree.pfmetcorr_ey;
            float met_y_recoilresoDown = analysisTree.pfmetcorr_ey;

            float met_unclMetUp_x    = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            float met_unclMetUp_y    = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            float met_unclMetDown_x  = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            float met_unclMetDown_y  = analysisTree.pfmetcorr_ey_UnclusteredEnDown;

            met = TMath::Sqrt(met_x*met_x + met_y*met_y);
            metphi = TMath::ATan2(met_y,met_x);
            metcov00 = analysisTree.pfmetcorr_sigxx;
            metcov01 = analysisTree.pfmetcorr_sigxy;
            metcov10 = analysisTree.pfmetcorr_sigyx;
            metcov11 = analysisTree.pfmetcorr_sigyy;

            int njetsforrecoil = njets;
            if (isW) njetsforrecoil = njets + 1;
            
            met_uncorr = met;
            metphi_uncorr = metphi;
            
            float pfmet_corr_x = met_x;
            float pfmet_corr_y = met_y;
            
            if ((isW||isDY||isSignal)&&!isData) {
                if (applySimpleRecoilCorrections) {
                    recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
		    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
		    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
		    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
		    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
                }
                else {
                    recoilMetCorrector.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                }
            }

            met_x = pfmet_corr_x;
            met_y = pfmet_corr_y;
            met = TMath::Sqrt(met_x*met_x+met_y*met_y);
            metphi = TMath::ATan2(met_y,met_x);
	    // update metLV (regognize changes if present)
	    metLV.SetXYZT(met_x,met_y,0.,met);
            met_recoilscaleUp   = TMath::Sqrt(met_x_recoilscaleUp*met_x_recoilscaleUp+met_y_recoilscaleUp*met_y_recoilscaleUp);
            met_recoilscaleDown = TMath::Sqrt(met_x_recoilscaleDown*met_x_recoilscaleDown+met_y_recoilscaleDown*met_y_recoilscaleDown);
            met_recoilresoUp    = TMath::Sqrt(met_x_recoilresoUp*met_x_recoilresoUp+met_y_recoilresoUp*met_y_recoilresoUp);
            met_recoilresoDown  = TMath::Sqrt(met_x_recoilresoDown*met_x_recoilresoDown+met_y_recoilresoDown*met_y_recoilresoDown);
	    TLorentzVector metLV_recoilscaleUp; metLV_recoilscaleUp.SetXYZT(met_x_recoilscaleUp, met_y_recoilscaleUp, 0., met_recoilscaleUp);
	    TLorentzVector metLV_recoilscaleDown; metLV_recoilscaleDown.SetXYZT(met_x_recoilscaleDown, met_y_recoilscaleDown, 0., met_recoilscaleDown);
	    TLorentzVector metLV_recoilresoUp; metLV_recoilresoUp.SetXYZT(met_x_recoilresoUp, met_y_recoilresoUp, 0., met_recoilresoUp);
	    TLorentzVector metLV_recoilresoDown; metLV_recoilresoDown.SetXYZT(met_x_recoilresoDown, met_y_recoilresoDown, 0., met_recoilresoDown);

            met_unclMetUp = TMath::Sqrt(met_unclMetUp_x*met_unclMetUp_x+met_unclMetUp_y*met_unclMetUp_y);
            metphi_unclMetUp = TMath::ATan2(met_unclMetUp_y,met_unclMetUp_x);
            met_unclMetDown = TMath::Sqrt(met_unclMetDown_x*met_unclMetDown_x+met_unclMetDown_y*met_unclMetDown_y);
            metphi_unclMetDown = TMath::ATan2(met_unclMetDown_y,met_unclMetDown_x);

            float genmet_ex = analysisTree.genmet_ex;
            float genmet_ey = analysisTree.genmet_ey;

            genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
            genmetphi = TMath::ATan2(genmet_ey,genmet_ex);
       
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
            
            float vectorX = met_x + muonLV.Px() + electronLV.Px();
            float vectorY = met_y + muonLV.Py() + electronLV.Py();
            float vectorVisX = muonLV.Px() + electronLV.Px();
            float vectorVisY = muonLV.Py() + electronLV.Py();
            pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;

	    double px_escaleUp = (1+eleScale) * pt_1 * TMath::Cos(phi_1);
	    double py_escaleUp = (1+eleScale) * pt_1 * TMath::Sin(phi_1);
	    
	    double px_e = pt_1 * TMath::Cos(phi_1);
	    double py_e = pt_1 * TMath::Sin(phi_1);
	    
	    double px_escaleDown = (1-eleScale) * pt_1 * TMath::Cos(phi_1);
	    double py_escaleDown = (1-eleScale) * pt_1 * TMath::Sin(phi_1);
	    
	    double metx_escaleUp = met_x + px_e - px_escaleUp;
	    double mety_escaleUp = met_y + py_e - py_escaleUp;
	    double metx_escaleDown = met_x + px_e - px_escaleDown;
	    double mety_escaleDown = met_y + py_e - py_escaleDown;

            // computation of DZeta variable
            // pfmet
            computeDzeta(met_x,met_y,
                         zetaX,zetaY,pzetavis,pzetamiss,dzeta);
            
            // genmet
            computeDzeta(analysisTree.genmet_ex,analysisTree.genmet_ey,
                         zetaX,zetaY,pzetavis,pzetamiss_genmet,dzeta_genmet);
            
	    float met_escaleUp = TMath::Sqrt(metx_escaleUp*metx_escaleUp+mety_escaleUp*mety_escaleUp);
	    float met_escaleDown = TMath::Sqrt(metx_escaleDown*metx_escaleDown+mety_escaleDown*mety_escaleDown);
            
            TLorentzVector metLV_escaleUp; metLV_escaleUp.SetXYZT(metx_escaleUp,mety_escaleUp,0.,met_escaleUp);
            TLorentzVector metLV_escaleDown; metLV_escaleDown.SetXYZT(metx_escaleDown,mety_escaleDown,0.,met_escaleDown);

            mt_1 = mT(electronLV,metLV);
            mt_2 = mT(muonLV,metLV);
            mtmax = TMath::Max(float(mt_1),float(mt_2));
            
            mCDF = (muonLV+electronLV+metLV).M();

            TLorentzVector metLV_unclMetUp; metLV_unclMetUp.SetXYZT(met_unclMetUp_x,met_unclMetUp_y,0.,met_unclMetUp);
            TLorentzVector metLV_unclMetDown; metLV_unclMetDown.SetXYZT(met_unclMetDown_x,met_unclMetDown_y,0.,met_unclMetDown);
            // computing total transverse mass
            mTtot = totalTransverseMass        ( muonLV ,     electronLV , metLV);
            mTemu   = mT(electronLV,muonLV);
            dphi_mumet = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
                                    metLV.Px(),metLV.Py());
            dphi_emet  = dPhiFrom2P(electronLV.Px(),electronLV.Py(),
                                    metLV.Px(),metLV.Py());
            
            bdt     = reader->EvaluateMVA("BDT");
            bdt_ggh = readerGGH->EvaluateMVA("BDT");
            bdt_bbh = readerBBH->EvaluateMVA("BDT");
            
	    mTdileptonMET = mT(dileptonLV,metLV);

	    pt_ttjj = -10;
	    if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
	      pt_ttjj = (muonLV+electronLV+metLV+jet1+jet2).Pt();
	    }
	    pt_tt = (muonLV+electronLV+metLV).Pt();

            m_sv   = -10;
            mt_sv  = -10;
	    pt_sv  = -10;
            eta_sv = -10;
            phi_sv = -10;
            
            bool checkSV = false;
            if (sync) checkSV = computeSVFitMass;
            else checkSV = computeSVFitMass && dzeta>-50 && iso_1<0.15 && iso_2<0.2 && trg_muonelectron > 0.5;
            
	    TMatrixD covMET(2, 2);
            if (checkSV) {
                    // covariance matrix MET
                    covMET[0][0] =  metcov00;
                    covMET[1][0] =  metcov10;
                    covMET[0][1] =  metcov01;
                    covMET[1][1] =  metcov11;

                    // define electron 4-vector
                    classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, 
                                                 pt_1, eta_1, phi_1, 0.51100e-3); 
                    // define muon 4-vector
                    classic_svFit::MeasuredTauLepton svFitMu(classic_svFit::MeasuredTauLepton::kTauToMuDecay,  
                                                                pt_2, eta_2, phi_2, 105.658e-3); 
                    
                    // central value
                    ClassicSVfit algo = SVFitMassComputation(svFitEle, svFitMu,
                                                                         met_x, met_y,
                                                                         covMET, inputFile_visPtResolution);
                    
                    m_sv =static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass();
                    mt_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass(); // return value of transverse svfit mass is in units of GeV
                    if ( !algo.isValidSolution() ) 
                        std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
                    pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt(); 
                    // std::cout<<"pt: "<<pt_sv<<std::endl;
                    eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
                    phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
                    
                    float px_sv = pt_sv*TMath::Cos(phi_sv);
                    float py_sv = pt_sv*TMath::Sin(phi_sv);
                    
                    float msvmet_ex = px_sv - dileptonLV.Px();
                    float msvmet_ey = py_sv - dileptonLV.Py();

                    msvmet = TMath::Sqrt(msvmet_ex*msvmet_ex+msvmet_ey*msvmet_ey);
                    msvmetphi = TMath::ATan2(msvmet_ey,msvmet_ex);
                    
            }
            gen_match_1 = 6;
            gen_match_2 = 6;
            
            float minDR_1 = 0.2;
            float minDR_2 = 0.2;
            unsigned int gen_1 = 0;
            bool bgen_1 = false;
            unsigned int gen_2 = 0;
            bool bgen_2 = false;
            
            if (!isData) {
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
                    if (isAnyType && analysisTree.genparticles_status[igen]==1) {
                        float etaGen = genLV.Eta();
                        float phiGen = genLV.Phi();
                        float deltaR_1 = deltaR(eta_1,phi_1,
                                                etaGen,phiGen);
                        if (deltaR_1<minDR_1) {
                            minDR_1 = deltaR_1;
                            gen_1 = igen;
                            bgen_1 = true;
                            if (type1) gen_match_1 = 1;
                            else if (type2) gen_match_1 = 2;
                            else if (type3) gen_match_1 = 3;
                            else if (type4) gen_match_1 = 4;
                        }
                        
                        float deltaR_2 = deltaR(eta_2,phi_2,
                                                etaGen,phiGen);
                        if (deltaR_2<minDR_2) {
                            minDR_2 = deltaR_2;
                            gen_2 = igen;
                            bgen_2 = true;
                            if (type1) gen_match_2 = 1;
                            else if (type2) gen_match_2 = 2;
                            else if (type3) gen_match_2 = 3;
                            else if (type4) gen_match_2 = 4;
                        }
                    }
                }
                TLorentzVector genElectronLV; genElectronLV.SetPtEtaPhiM(pt_1,eta_1,phi_1,0.51100e-3);
                TLorentzVector genMuonLV; genMuonLV.SetPtEtaPhiM(pt_2,eta_2,phi_2,105.658e-3);
                if (bgen_1) genElectronLV.SetXYZT(analysisTree.genparticles_px[gen_1],
                                                  analysisTree.genparticles_py[gen_1],
                                                  analysisTree.genparticles_pz[gen_1],
                                                  analysisTree.genparticles_e[gen_1]);
                if (bgen_2) genMuonLV.SetXYZT(analysisTree.genparticles_px[gen_2],
                                              analysisTree.genparticles_py[gen_2],
                                              analysisTree.genparticles_pz[gen_2],
                                              analysisTree.genparticles_e[gen_2]);
                
                //	TLorentzVector genDileptonLV = genElectronLV + genMuonLV;
                TLorentzVector genMetLV; genMetLV.SetXYZM(genmet_ex,genmet_ey,
                                                          0,0);
                
                mTemu_gen = mT(genElectronLV,genMuonLV);
                mTemet_gen = mT(genElectronLV,genMetLV);
                mTmumet_gen = mT(genMuonLV,genMetLV);
                mTtot_gen = totalTransverseMass(genElectronLV,genMuonLV,genMetLV);
                
                if (gen_match_1>2&&gen_match_2>3) { 
                    isZTT = true;
                    isZLL = false;
                }
                else {
                    isZTT = false;
                    isZLL = true;
                }
                if (gen_match_1 == 3 && gen_match_2 ==4) veto_embedded = true;
    
            correctionWS_qcd->var("e_pt")->setVal(pt_1);
            correctionWS_qcd->var("m_pt")->setVal(pt_2);
            correctionWS_qcd->var("njets")->setVal(njets);
            correctionWS_qcd->var("dR")->setVal(dr_tt);
            double_t em_qcd_osss_binned = correctionWS_qcd->function("em_qcd_osss_binned")->getVal();
         
            double_t em_qcd_osss_binned_rate_up = correctionWS_qcd->function("em_qcd_osss_rateup_binned")->getVal();
            double_t em_qcd_osss_binned_rate_down = correctionWS_qcd->function("em_qcd_osss_ratedown_binned")->getVal();
            double_t em_qcd_osss_binned_shape_up = correctionWS_qcd->function("em_qcd_osss_shapeup_binned")->getVal();
            double_t em_qcd_osss_binned_shape_down = correctionWS_qcd->function("em_qcd_osss_shapedown_binned")->getVal();
                        
            qcdweight = em_qcd_osss_binned;
            
            if (njets==0) {
               qcdweight_0jet_rate_up =  em_qcd_osss_binned_rate_up;
               qcdweight_0jet_rate_down = em_qcd_osss_binned_rate_down;   
               qcdweight_0jet_shape_up = em_qcd_osss_binned_shape_up;
               qcdweight_0jet_shape_down =  em_qcd_osss_binned_shape_down;
               qcdweight_1jet_rate_up =  em_qcd_osss_binned;
               qcdweight_1jet_rate_down = em_qcd_osss_binned; 
               qcdweight_1jet_shape_up = em_qcd_osss_binned;
               qcdweight_1jet_shape_down = em_qcd_osss_binned;
            }   
            else {
               qcdweight_0jet_rate_up =  em_qcd_osss_binned;
               qcdweight_0jet_rate_down = em_qcd_osss_binned;
               qcdweight_0jet_shape_up = em_qcd_osss_binned;
               qcdweight_0jet_shape_down =  em_qcd_osss_binned;
               qcdweight_1jet_rate_up =  em_qcd_osss_binned_rate_up;
               qcdweight_1jet_rate_down =  em_qcd_osss_binned_rate_down;
               qcdweight_1jet_shape_up = em_qcd_osss_binned_shape_up;
               qcdweight_1jet_shape_down =  em_qcd_osss_binned_shape_down;
            }
            
            qcdweight_iso_up = correctionWS_qcd->function("em_qcd_bothaiso_extrap_up")->getVal();
            qcdweight_iso_down = correctionWS_qcd->function("em_qcd_bothaiso_extrap_down")->getVal();
         




		double weightE = 1;
		double weightEUp = 1;
		double weightEDown = 1;
		if (gen_match_1==6) {
		  double dRmin = 0.5;
		  double ptJet = pt_1;
		  for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
		    double etaJet = analysisTree.pfjet_eta[ijet];
		    double phiJet = analysisTree.pfjet_phi[ijet];
		    double dRjetE = deltaR(etaJet,phiJet,eta_1,phi_1);
		    if (dRjetE<dRmin) {
		      dRmin = dRjetE;
		      ptJet = analysisTree.pfjet_pt[ijet];
		    }
		  }
		  // if (isTOP) {
		  //   weightE = 1;
		  //   weightEUp = 1;
		  //   weightEDown = 1;
		  // }
		  
        weightE     = a_jetEle + b_jetEle*ptJet;
        weightEUp   = a_jetEleUp + b_jetEleUp*ptJet;
        weightEDown = a_jetEleDown + b_jetEleDown*ptJet;
		  
		  fakeweight *= weightE;
		}
		else {
		  weightE = isoweight_1;
		  weightEUp = isoweight_1;
		  weightEDown = isoweight_1;
      }

		double weightMu = 1;
		double weightMuUp = 1;
		double weightMuDown = 1;
		if (gen_match_2==6) {
		  double dRmin = 0.5;
		  double ptJet = pt_2;
		  for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
		    double etaJet = analysisTree.pfjet_eta[ijet];
		    double phiJet = analysisTree.pfjet_phi[ijet];
		    double dRjetM = deltaR(etaJet,phiJet,eta_2,phi_2);
		    if (dRjetM<dRmin) {
		      dRmin = dRjetM;
		      ptJet = analysisTree.pfjet_pt[ijet];
		    }
		  }
        weightMu     = a_jetMu     + b_jetMu*ptJet;
        weightMuUp   = a_jetMuUp   + b_jetMuUp*ptJet;
        weightMuDown = a_jetMuDown + b_jetMuDown*ptJet;
        
		  fakeweight *= weightMu;
		}
		else {
		  weightMu = isoweight_2;
		  weightMuUp = isoweight_2;
		  weightMuDown = isoweight_2;
      }
		float effweight0 = effweight;
		effweight = effweight0 * weightE * weightMu;
		effweight_jetMuUp   = effweight0 * weightE     * weightMuUp;
		effweight_jetMuDown = effweight0 * weightE     * weightMuDown; 
		effweight_jetEUp    = effweight0 * weightEUp   * weightMu;
		effweight_jetEDown  = effweight0 * weightEDown * weightMu; 
            
      if (sync) weight = effweight * puweight * 0.979;
      // std::cout<<"effweight: "<<effweight/(idweight_1*idweight_2)<<std::endl;
      // std::cout<<"puweight: "<<puweight<<std::endl;
      // std::cout<<"e tracking: "<<idweight_1<<std::endl;
      // std::cout<<"mu tracking: "<<idweight_2<<std::endl;
            }

	    // First set the calibrated variables to the uncalibrated versions -> necessary for data
	    d0_1_cal=d0_1;
	    dZ_1_cal=dZ_1;
	    d0_2_cal=d0_2;
	    dZ_2_cal=dZ_2;

	    // if (!isData){               
	    //    //check if particle is prompt or non-prompt
       //         if (gen_match_1 == 3) {  
       //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_1,eta_1,"d0_1",d0_1, d0_1_cal);
       //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_1,eta_1,"dZ_1",dZ_1, dZ_1_cal);
       //            //d0_1_cal = doquantileshift_d01_np.shift(d0_1, 0.55);
       //            //dZ_1_cal = doquantileshift_dZ1_np.shift(dZ_1, 0.55);
                  
       //         }
       //         else if (gen_match_1 == 1){
       //            calibrateIP.DoCalibrationForPromptElectrons(pt_1,eta_1,"d0_1",d0_1, d0_1_cal);
       //            calibrateIP.DoCalibrationForPromptElectrons(pt_1,eta_1,"dZ_1",dZ_1, dZ_1_cal);
       //            //d0_1_cal = doquantileshift_d01_pele.shift(d0_1, 0.55);
       //            //dZ_1_cal = doquantileshift_dZ1_pele.shift(dZ_1, 0.55);
       //         }
       //         else {
       //            d0_1_cal = d0_1;
       //            dZ_1_cal = dZ_1;
       //         }
       //         if (gen_match_2 == 4) {  
       //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_2,eta_2,"d0_2",d0_2,d0_2_cal);
       //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_2,eta_2,"dZ_2",dZ_2,dZ_2_cal);
       //            //d0_2_cal = doquantileshift_d02_np.shift(d0_2, 0.55);
       //            //dZ_2_cal = doquantileshift_dZ2_np.shift(dZ_2, 0.55);
       //         }
       //         else if (gen_match_2 == 2){
       //            calibrateIP.DoCalibrationForPromptMuons(pt_2,eta_2,"d0_2",d0_2, d0_2_cal);
       //            calibrateIP.DoCalibrationForPromptMuons(pt_2,eta_2,"dZ_2",dZ_2, dZ_2_cal);
       //            //d0_2_cal = doquantileshift_d02_pmu.shift(d0_2, 0.55);
       //            //dZ_2_cal = doquantileshift_dZ2_pmu.shift(dZ_2, 0.55);
       //         }
       //         else{
       //            d0_2_cal =d0_2;
       //            dZ_2_cal =dZ_2;
       //         }
       //      }       

	    // Add for all relevant variables the met uncertainty
	    for(auto &uncert : uncertainty_map){
	      uncert.second.electronLV = electronLV;
	      uncert.second.muonLV     = muonLV;
	      uncert.second.metLV      = metLV;
	      uncert.second.jet1LV     = jet1;
	      uncert.second.jet2LV     = jet2;
	    }

	    uncertainty_map.at("unclMetUp").metLV = metLV_unclMetUp;
	    uncertainty_map.at("unclMetDown").metLV = metLV_unclMetDown;
	    uncertainty_map.at("escaleUp").metLV = metLV_escaleUp;
	    uncertainty_map.at("escaleUp").electronLV = electronUpLV;
	    uncertainty_map.at("escaleDown").metLV = metLV_escaleDown;
	    uncertainty_map.at("escaleDown").electronLV = electronDownLV;
	    //uncertainty_map.at("mscaleUp").muonLV = muonUpLV;
	    //uncertainty_map.at("mscaleDown").muonLV = muonDownLV;
	    uncertainty_map.at("recoilscaleUp").metLV = metLV_recoilscaleUp;
	    uncertainty_map.at("recoilscaleDown").metLV = metLV_recoilscaleDown;
	    uncertainty_map.at("recoilresoUp").metLV = metLV_recoilresoUp;
	    uncertainty_map.at("recoilresoDown").metLV = metLV_recoilresoDown;
	    uncertainty_map.at("jecUncEta0To5Up").metLV = metLV_jecUnc.at("jecUncEta0To5Up");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Up").jet1LV = jet1LV_jecUnc.at("jecUncEta0To5Up");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To5Up").jet2LV = jet2LV_jecUnc.at("jecUncEta0To5Up");
	    uncertainty_map.at("jecUncEta0To5Down").metLV = metLV_jecUnc.at("jecUncEta0To5Down");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Down").jet1LV = jet1LV_jecUnc.at("jecUncEta0To5Down");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To5Down").jet2LV = jet2LV_jecUnc.at("jecUncEta0To5Down");
	    uncertainty_map.at("jecUncEta0To3Up").metLV = metLV_jecUnc.at("jecUncEta0To3Up");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Up").jet1LV = jet1LV_jecUnc.at("jecUncEta0To3Up");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To3Up").jet2LV = jet2LV_jecUnc.at("jecUncEta0To3Up");
	    uncertainty_map.at("jecUncEta0To3Down").metLV = metLV_jecUnc.at("jecUncEta0To3Down");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Down").jet1LV = jet1LV_jecUnc.at("jecUncEta0To3Down");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To3Down").jet2LV = jet2LV_jecUnc.at("jecUncEta0To3Down");
	    uncertainty_map.at("jecUncEta3To5Up").metLV = metLV_jecUnc.at("jecUncEta3To5Up");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Up").jet1LV = jet1LV_jecUnc.at("jecUncEta3To5Up");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta3To5Up").jet2LV = jet2LV_jecUnc.at("jecUncEta3To5Up");
	    uncertainty_map.at("jecUncEta3To5Down").metLV = metLV_jecUnc.at("jecUncEta3To5Down");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Down").jet1LV = jet1LV_jecUnc.at("jecUncEta3To5Down");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncEta3To5Down").jet2LV = jet2LV_jecUnc.at("jecUncEta3To5Down");
	    uncertainty_map.at("jecUncRelativeBalUp").metLV = metLV_jecUnc.at("jecUncRelativeBalUp");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalUp").jet1LV = jet1LV_jecUnc.at("jecUncRelativeBalUp");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeBalUp").jet2LV = jet2LV_jecUnc.at("jecUncRelativeBalUp");
	    uncertainty_map.at("jecUncRelativeBalDown").metLV = metLV_jecUnc.at("jecUncRelativeBalDown");
	    if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalDown").jet1LV = jet1LV_jecUnc.at("jecUncRelativeBalDown");
	    if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeBalDown").jet2LV = jet2LV_jecUnc.at("jecUncRelativeBalDown");

	    for(auto &uncert : uncertainty_map){

	      bool is_data_or_embedded = isData || isEmbedded;

	      propagate_uncertainty( uncert.first,
				     uncert.second.metLV, covMET, inputFile_visPtResolution,
				     uncert.second.muonLV,
				     uncert.second.electronLV,
				     uncert.second.jet1LV,
				     uncert.second.jet2LV,
				     uncert.second.container,
				     is_data_or_embedded, checkSV);
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

