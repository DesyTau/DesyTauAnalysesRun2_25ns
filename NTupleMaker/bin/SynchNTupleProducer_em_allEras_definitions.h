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
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"

#include "TSystem.h"
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
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "HTT-utilities/QCDModelingEMu/interface/QCDModelForEMu.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/helper_functions_em.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ggF_qcd_uncertainty_2017.h"
#include "DesyTauAnalyses/NTupleMaker/interface/qq2Hqq_uncert_scheme.cpp"

#include "TCut.h"

#include "DesyTauAnalyses/NTupleMaker/interface/btagSF.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "ZZMatrixElement/MELA/interface/TUtil.hh"

using namespace std;

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

LorentzVector tau1P4;
LorentzVector tau2P4;

unsigned int minRun = 99999999;
unsigned int maxRun = 0;
Int_t           run;
Int_t           lumi;
ULong64_t       evt;
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
Float_t         qcdweight_2jet_rate_up;
Float_t         qcdweight_2jet_rate_down;
Float_t         qcdweight_0jet_shape_up;
Float_t         qcdweight_0jet_shape_down;
Float_t         qcdweight_1jet_shape_up;
Float_t         qcdweight_1jet_shape_down;
Float_t         qcdweight_2jet_shape_up;
Float_t         qcdweight_2jet_shape_down;

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

Float_t         prefiringweight;
Float_t         prefiringweightup;
Float_t         prefiringweightdown;

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

Float_t         pt_1;
Float_t         phi_1;
Float_t         eta_1;
Float_t         m_1;
Int_t           q_1;
Float_t         iso_1;
//Float_t         mva_1;
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
Int_t           njets_jecUncRelativeSampleYearUp;
Int_t           njets_jecUncRelativeSampleYearDown;
Int_t           njets_jecUncEC2Up;
Int_t           njets_jecUncEC2Down;
Int_t           njets_jecUncFlavorQCDUp;
Int_t           njets_jecUncFlavorQCDDown;
Int_t           njets_jecUncEC2YearUp;
Int_t           njets_jecUncEC2YearDown;
Int_t           njets_jecUncHFYearUp;
Int_t           njets_jecUncHFYearDown;
Int_t           njets_jecUncAbsoluteYearUp;
Int_t           njets_jecUncAbsoluteYearDown;
Int_t           njets_jecUncBBEC1YearUp;
Int_t           njets_jecUncBBEC1YearDown;

Int_t           njetspt20;

Float_t         jpt_1;

Float_t         jeta_1;
Float_t         jphi_1;
Float_t         jptraw_1;
Int_t           gen_match_1;

Float_t         jpt_2;

Float_t         jeta_2;
Float_t         jphi_2;
Float_t         jptraw_2;
Int_t           gen_match_2;

Float_t         mjj;
Float_t         dijetphi;
Float_t         dijetpt;
Float_t         jdeta;

Int_t           nbtag;
Int_t           nbtag_mistagUp;
Int_t           nbtag_mistagDown;
Int_t           nbtag_btagUp;
Int_t           nbtag_btagDown;
Int_t           nbtag_noSF;
Double_t         bpt;
Double_t         beta_1;
Double_t         bphi;

Float_t         nuPx;  // neutrinos from t -> W(lv) + b
Float_t         nuPy;  // or from tau->l+v+v events
Float_t         nuPz;  // (x,y,z) components
Float_t         nuPt;  // pT
Float_t         nuPhi; // phi

Float_t         mtBoson_gen;
    
Float_t         mTtot_gen;
Float_t         mTemu_gen;
Float_t         mTemet_gen;
Float_t         mTmumet_gen;

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

Bool_t isPromptZMM;
Bool_t isPromptZEE;
Bool_t isZTTMM;
Bool_t isZTTEE;
Bool_t isZTTEM;

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
Int_t htxs_stage1p1cat;

Float_t THU_ggH_Mu;
Float_t THU_ggH_Res;
Float_t THU_ggH_Mig01;
Float_t THU_ggH_Mig12;
Float_t THU_ggH_VBF2j;
Float_t THU_ggH_VBF3j;
Float_t THU_ggH_PT60;
Float_t THU_ggH_PT120;
Float_t THU_ggH_qmtop;

Float_t THU_qqH_TOT;
Float_t THU_qqH_PTH200;
Float_t THU_qqH_Mjj60;
Float_t THU_qqH_Mjj120;
Float_t THU_qqH_Mjj350;
Float_t THU_qqH_Mjj700;
Float_t THU_qqH_Mjj1000;
Float_t THU_qqH_Mjj1500;
Float_t THU_qqH_25;
Float_t THU_qqH_JET01;

float MaxBJetPt = 1000.;
float MaxLJetPt = 1000.;
float MinLJetPt = 20.;
float MinBJetPt = 20.; // !!!!!

TH1F *tagEff_B;
TH1F *tagEff_C;
TH1F *tagEff_Light;
TRandom3 r;

TH2D * histZMassPtWeights;

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
                           "pt_sv", // 24
                           "eta_sv", // 25
                           "phi_sv", // 26
                           "mt_sv", // 27
                           "mTemu"}; // 28


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
inputs eresoUp;
inputs eresoDown;
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
inputs jecUncRelativeSampleYearUp;
inputs jecUncRelativeSampleYearDown;
inputs jecUncEC2Up;
inputs jecUncEC2Down;
inputs jecUncFlavorQCDUp;
inputs jecUncFlavorQCDDown;
inputs jecUncAbsoluteYearUp;
inputs jecUncAbsoluteYearDown;
inputs jecUncBBEC1YearUp;
inputs jecUncBBEC1YearDown;
inputs jecUncEC2YearUp;
inputs jecUncEC2YearDown;
inputs jecUncHFYearUp;
inputs jecUncHFYearDown;

map<TString, inputs> uncertainty_map = { { "unclMetUp" , unclMetUp },
                                         { "unclMetDown" , unclMetDown },
                                         { "escaleUp" , escaleUp },
                                         { "escaleDown" , escaleDown },
                                         { "eresoUp" , eresoUp },
                                         { "eresoDown" , eresoDown },
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
                                         { "jecUncRelativeSampleYearUp" , jecUncRelativeSampleYearUp },
                                         { "jecUncRelativeSampleYearDown" , jecUncRelativeSampleYearDown },
                                         { "jecUncEC2Up" , jecUncEC2Up},
                                         { "jecUncEC2Down" , jecUncEC2Down},
                                         { "jecUncFlavorQCDUp" , jecUncFlavorQCDUp},
                                         { "jecUncFlavorQCDDown" , jecUncFlavorQCDDown},
                                         { "jecUncAbsoluteYearUp" , jecUncAbsoluteYearUp},
                                         { "jecUncAbsoluteYearDown" , jecUncAbsoluteYearDown},
                                         { "jecUncBBEC1YearUp" , jecUncBBEC1YearUp},
                                         { "jecUncBBEC1YearDown" , jecUncBBEC1YearDown},
                                         { "jecUncEC2YearUp" , jecUncEC2YearUp},
                                         { "jecUncEC2YearDown" , jecUncEC2YearDown},
                                         { "jecUncHFYearUp" , jecUncHFYearUp},
                                         { "jecUncHFYearDown" , jecUncHFYearDown}                                         
};

const int nsrc_Eta0To5 = 8; // == Absolute according to https://docs.google.com/spreadsheets/d/1Feuj1n0MdotcPq19Mht7SUIgvkXkA4hiB0BxEuBShLw/edit#gid=1345121349
const char* srcnames_Eta0To5[nsrc_Eta0To5] = {"SinglePionECAL",
                                              "SinglePionHCAL",
                                              "AbsoluteMPFBias",
                                              "AbsoluteScale",
                                              "Fragmentation",
                                              "PileUpDataMC",
                                              "RelativeFSR",
                                              "PileUpPtRef"};
const int nsrc_AbsoluteYear = 3;
const char* srcnames_AbsoluteYear[nsrc_AbsoluteYear] = {"AbsoluteStat",
                                                        "RelativeStatFSR", 
                                                        "TimePtEta"};


const int nsrc_Eta0To3 = 3; // == BBEC1 according to https://docs.google.com/spreadsheets/d/1Feuj1n0MdotcPq19Mht7SUIgvkXkA4hiB0BxEuBShLw/edit#gid=1345121349
const char* srcnames_Eta0To3[nsrc_Eta0To3] = {"PileUpPtEC1",
                                              "PileUpPtBB",
                                              "RelativePtBB"};
const int nsrc_BBEC1Year= 3;
const char* srcnames_BBEC1Year[nsrc_BBEC1Year] = {"RelativeJEREC1",
                                                  "RelativePtEC1",
                                                  "RelativeStatEC"};

const int nsrc_Eta3To5 = 3;// == HF according to https://docs.google.com/spreadsheets/d/1Feuj1n0MdotcPq19Mht7SUIgvkXkA4hiB0BxEuBShLw/edit#gid=1345121349
const char* srcnames_Eta3To5[nsrc_Eta3To5] = {"RelativePtHF",
                                              "PileUpPtHF",
                                              "RelativeJERHF"};

const int nsrc_HFYear = 1;
const char* srcnames_HFYear[nsrc_HFYear] = {"RelativeStatHF"};

const int nsrc_RelativeBal = 1;
const char* srcnames_RelativeBal[nsrc_RelativeBal] = {"RelativeBal"};

const int nsrc_RelativeSampleYear = 1;
const char* srcnames_RelativeSampleYear[nsrc_RelativeSampleYear] = {"RelativeSample"};

const int nsrc_EC2 = 1;
const char*srcnames_EC2[nsrc_EC2] = {"PileUpPtEC2"};


const int nsrc_EC2Year = 2;
const char*srcnames_EC2Year[nsrc_EC2Year] = {"RelativeJEREC2", 
                                             "RelativePtEC2"};

const int nsrc_FlavorQCD = 1;
const char* srcnames_FlavorQCD[nsrc_FlavorQCD] = {"FlavorQCD"};

std::vector<JetCorrectionUncertainty*> vsrc_Eta0To5(nsrc_Eta0To5);
std::vector<JetCorrectionUncertainty*> vsrc_Eta0To3(nsrc_Eta0To3);
std::vector<JetCorrectionUncertainty*> vsrc_Eta3To5(nsrc_Eta3To5);
std::vector<JetCorrectionUncertainty*> vsrc_RelativeBal(nsrc_RelativeBal);
std::vector<JetCorrectionUncertainty*> vsrc_RelativeSampleYear(nsrc_RelativeSampleYear);
std::vector<JetCorrectionUncertainty*> vsrc_EC2(nsrc_EC2);
std::vector<JetCorrectionUncertainty*> vsrc_FlavorQCD(nsrc_FlavorQCD);
std::vector<JetCorrectionUncertainty*> vsrc_HFYear(nsrc_HFYear);
std::vector<JetCorrectionUncertainty*> vsrc_EC2Year(nsrc_EC2Year);
std::vector<JetCorrectionUncertainty*> vsrc_BBEC1Year(nsrc_BBEC1Year);
std::vector<JetCorrectionUncertainty*> vsrc_AbsoluteYear(nsrc_AbsoluteYear);



TTree *tree = new TTree("TauCheck","TauCheck");

float topPt = -1;
float antitopPt = -1;

bool isGSfound = false;
std::vector<TLorentzVector> promptTausFirstCopy; 
std::vector<TLorentzVector> promptTausLastCopy;
std::vector<TLorentzVector> promptElectrons;
std::vector<TLorentzVector> promptMuons;
std::vector<TLorentzVector> promptNeutrinos;
std::vector<TLorentzVector> tauNeutrinos;

TLorentzVector promptTausLV;
TLorentzVector promptVisTausLV;
TLorentzVector zBosonLV;
TLorentzVector wBosonLV;
TLorentzVector hBosonLV;
TLorentzVector promptElectronsLV;
TLorentzVector promptMuonsLV;
TLorentzVector promptNeutrinosLV;
TLorentzVector tauNeutrinosLV;
TLorentzVector wDecayProductsLV;
TLorentzVector fullVLV;
TLorentzVector visVLV;

TLorentzVector genBosonLV;
TLorentzVector genVisBosonLV;

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

bool isMuon23matched = false;
bool isMuon8matched  = false;
bool isElectron23matched = false;
bool isElectron12matched = false;
bool isMu23 = false;
bool isMu8 = false;
bool isMu23dz = false;
bool isMu8dz  = false;
bool isEle23 = false;
bool isEle12 = false;
bool isEle23dz = false;
bool isEle12dz = false;
float Ele23EffData = 1;
float Ele12EffData = 1;
float Mu23EffData = 1;
float Mu8EffData = 1;

float Ele23EffMC = 1;
float Ele12EffMC = 1;
float Mu23EffMC = 1;
float Mu8EffMC = 1;

vector<unsigned int> jets;
vector<unsigned int> jetspt20;
vector<unsigned int> bjets;
vector<unsigned int> bjets_mistagUp;
vector<unsigned int> bjets_mistagDown;
vector<unsigned int> bjets_btagUp;
vector<unsigned int> bjets_btagDown;
vector<unsigned int> bjets_nocleaned;
vector<unsigned int> bjetsRaw;

TLorentzVector jetLV;
TLorentzVector jet1;
TLorentzVector jet2;
map<TString,TLorentzVector> jet1LV_jecUnc;
map<TString,TLorentzVector> jet2LV_jecUnc;
map<TString,TLorentzVector> metLV_jecUnc;

int indexLeadingJet = -1;
float ptLeadingJet = -1;

int indexSubLeadingJet = -1;
float ptSubLeadingJet = -1;

int indexLeadingBJet = -1;
float ptLeadingBJet = -1;

bool checkSV = false;
bool checkFastMTT = false;
bool passesPreSel = false; 
bool passesFastMTTPreSel = false;

bool isSVFitUsed = false;
bool isFastMTTUsed = false;

bool isPuppiMETUsed = false;

// MELA outputs
// 1. Matrix element variables for different hypotheses (VBF Higgs, ggH + 2 jets, Z + 2 jets)
float ME_vbf, ME_ggh, ME_z2j_1, ME_z2j_2;
// 2. Energy transfer (Q^2) variables
float ME_q2v1, ME_q2v2;
// 3. Angle variables
float ME_costheta1, ME_costheta2, ME_phi, ME_costhetastar, ME_phi1;
// 4. Main BG vs. Higgs discriminators
float ME_vbf_vs_Z, ME_ggh_vs_Z, ME_vbf_vs_ggh;

void SetupTree(){
 
   tree->Branch("run", &run, "run/I");
   tree->Branch("lumi", &lumi, "lumi/I");
   tree->Branch("evt", &evt, "evt/l");
   tree->Branch("npv", &npv, "npv/I");
   tree->Branch("npu", &npu, "npu/I");
   tree->Branch("rho", &rho, "rho/F");
   
   tree->Branch("isZLL",&isZLL,"isZLL/O");
   tree->Branch("isZEE",&isZEE,"isZEE/O");
   tree->Branch("isZMM",&isZMM,"isZMM/O");
   tree->Branch("isZTT",&isZTT,"isZTT/O");
   tree->Branch("veto_embedded",&veto_embedded,"veto_embedded/O");
   tree->Branch("isPromptZMM",&isPromptZMM,"isPromptZMM/O");
   tree->Branch("isPromptZEE",&isPromptZEE,"isPromptZEE/O");
   tree->Branch("isZTTEM",&isZTTEM,"isZTTEM/O");
   tree->Branch("isZTTMM",&isZTTMM,"isZTTMM/O");
   tree->Branch("isZTTEE",&isZTTEE,"isZTTEE/O"); 
   
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
 
   tree->Branch("THU_qqH_TOT", &THU_qqH_TOT, "THU_qqH_TOT/F");
   tree->Branch("THU_qqH_PTH200", &THU_qqH_PTH200, "THU_qqH_PTH200/F");
   tree->Branch("THU_qqH_Mjj60", &THU_qqH_Mjj60, "THU_qqH_Mjj60/F");
   tree->Branch("THU_qqH_Mjj120", &THU_qqH_Mjj120, "THU_qqH_Mjj120/F");
   tree->Branch("THU_qqH_Mjj350", &THU_qqH_Mjj350, "THU_qqH_Mjj350/F");
   tree->Branch("THU_qqH_Mjj700", &THU_qqH_Mjj700, "THU_qqH_Mjj700/F");
   tree->Branch("THU_qqH_Mjj1000" , &THU_qqH_Mjj1000, "THU_qqH_Mjj1000/F");
   tree->Branch("THU_qqH_Mjj1500", &THU_qqH_Mjj1500, "THU_qqH_Mjj1500/F");
   tree->Branch("THU_qqH_25", &THU_qqH_25, "THU_qqH_25/F");
   tree->Branch("THU_qqH_JET01", &THU_qqH_JET01, "THU_qqH_JET01/F");
   
   tree->Branch("higgspt_HTXS",&higgspt_HTXS,"higgspt_HTXS/F");
   tree->Branch("njets_HTXS",&njets_HTXS,"njets_HTXS/I");
   tree->Branch("htxs_stage0cat",&htxs_stage0cat,"htxs_stage0cat/I");
   tree->Branch("htxs_stage1cat",&htxs_stage1cat,"htxs_stage1cat/I");
   tree->Branch("htxs_stage1p1cat",&htxs_stage1p1cat,"htxs_stage1p1cat/I");
   
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
   tree->Branch("qcdweight_0jet_rate_down",&qcdweight_0jet_rate_down,"qcdweight_0jet_rate_down/F");
   tree->Branch("qcdweight_1jet_rate_up",&qcdweight_1jet_rate_up,"qcdweight_1jet_rate_up/F");
   tree->Branch("qcdweight_1jet_rate_down",&qcdweight_1jet_rate_down,"qcdweight_1jet_rate_down/F");
   tree->Branch("qcdweight_2jet_rate_up",&qcdweight_2jet_rate_up,"qcdweight_2jet_rate_up/F");
   tree->Branch("qcdweight_2jet_rate_down",&qcdweight_2jet_rate_down,"qcdweight_2jet_rate_down/F");
   tree->Branch("qcdweight_0jet_shape_up",&qcdweight_0jet_shape_up,"qcdweight_0jet_shape_up/F");
   tree->Branch("qcdweight_0jet_shape_down",&qcdweight_0jet_shape_down,"qcdweight_0jet_shape_down/F");  
   tree->Branch("qcdweight_1jet_shape_up",&qcdweight_1jet_shape_up,"qcdweight_1jet_shape_up/F");
   tree->Branch("qcdweight_1jet_shape_down",&qcdweight_1jet_shape_down,"qcdweight_1jet_shape_down/F");
   tree->Branch("qcdweight_2jet_shape_up",&qcdweight_2jet_shape_up,"qcdweight_2jet_shape_up/F");
   tree->Branch("qcdweight_2jet_shape_down",&qcdweight_2jet_shape_down,"qcdweight_2jet_shape_down/F");
   
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

   tree->Branch("prefiringweight", &prefiringweight, "prefiringweight/F");
   tree->Branch("prefiringweightup", &prefiringweightup, "prefiringweightup/F");
   tree->Branch("prefiringweightdown", &prefiringweightdown, "prefiringweightdown/F");
   
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
   
   tree->Branch("isSVFitUsed",    &isSVFitUsed,   "isSVFitUsed/O");
   tree->Branch("isFastMTTUsed",    &isFastMTTUsed,   "isFastMTTUsed/O");
   tree->Branch("m_sv",    &m_sv,   "m_sv/F");
   tree->Branch("mt_sv",   &mt_sv,  "mt_sv/F");
   
   tree->Branch("pt_sv",   &pt_sv,  "pt_sv/F");
   tree->Branch("eta_sv",  &eta_sv, "eta_sv/F");
   tree->Branch("phi_sv",  &phi_sv, "phi_sv/F");
    
   tree->Branch("pt_1", &pt_1, "pt_1/F");
   tree->Branch("phi_1", &phi_1, "phi_1/F");
   tree->Branch("eta_1", &eta_1, "eta_1/F");
   tree->Branch("m_1", &m_1, "m_1/F");
   tree->Branch("q_1", &q_1, "q_1/I");
   tree->Branch("iso_1", &iso_1, "iso_1/F");
   //tree->Branch("mva_1", &mva_1, "mva_1/F");
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
   tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/O");
   tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/O");
   
   tree->Branch("isPuppiMETUsed", &isPuppiMETUsed, "isPuppiMETUsed/O");
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
   tree->Branch("mtBoson_gen",&mtBoson_gen,"mtBoson_gen/F");
    
   tree->Branch("dphi_mumet",&dphi_mumet,"dphi_mumet/F");
   tree->Branch("dphi_emet",&dphi_emet,"dphi_emet/F");
      
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
   tree->Branch("njets_jecUncRelativeSampleYearUp", &njets_jecUncRelativeSampleYearUp, "njets_jecUncRelativeSampleYearUp/I");
   tree->Branch("njets_jecUncRelativeSampleYearDown", &njets_jecUncRelativeSampleYearDown, "njets_jecUncRelativeSampleYearDown/I");
   tree->Branch("njets_jecUncEC2Up", &njets_jecUncEC2Up, "njets_jecUncEC2Up/I");
   tree->Branch("njets_jecUncEC2Down", &njets_jecUncEC2Down, "njets_jecUncEC2Down/I");
   tree->Branch("njets_jecUncFlavorQCDUp", &njets_jecUncFlavorQCDUp, "njets_jecUncFlavorQCDUp/I");
   tree->Branch("njets_jecUncFlavorQCDDown", &njets_jecUncFlavorQCDDown, "njets_jecUncFlavorQCDDown/I");
   tree->Branch("njets_jecUncEC2YearUp", &njets_jecUncEC2YearUp, "njets_jecUncEC2YearUp/I");
   tree->Branch("njets_jecUncEC2YearDown", &njets_jecUncEC2YearDown, "njets_jecUncEC2YearDown/I");
   tree->Branch("njets_jecUncHFYearUp", &njets_jecUncHFYearUp, "njets_jecUncHFYearUp/I");
   tree->Branch("njets_jecUncHFYearDown", &njets_jecUncHFYearDown, "njets_jecUncHFYearDown/I");
   tree->Branch("njets_jecUncAbsoluteYearUp", &njets_jecUncAbsoluteYearUp, "njets_jecUncAbsoluteYearUp/I");
   tree->Branch("njets_jecUncAbsoluteYearDown", &njets_jecUncAbsoluteYearDown, "njets_jecUncAbsoluteYearDown/I");
   tree->Branch("njets_jecUncBBEC1YearUp", &njets_jecUncBBEC1YearUp, "njets_jecUncBBEC1YearUp/I");
   tree->Branch("njets_jecUncBBEC1YearDown", &njets_jecUncBBEC1YearDown, "njets_jecUncBBEC1YearDown/I");


   tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
   
   tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
   
   tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
   tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
   tree->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
    
   tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
   tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
   tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
   tree->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
   
   tree->Branch("mjj", &mjj, "mjj/F");
   tree->Branch("dijetphi", &dijetphi, "dijetphi/F");
   tree->Branch("dijetpt", &dijetpt, "dijetpt/F");
   tree->Branch("jdeta", &jdeta, "jdeta/F");
   
   tree->Branch("nbtag", &nbtag, "nbtag/I");
   tree->Branch("nbtag_mistagUp", &nbtag_mistagUp, "nbtag_mistagUp/I");
   tree->Branch("nbtag_mistagDown", &nbtag_mistagDown, "nbtag_mistagDown/I");
   tree->Branch("nbtag_btagUp", &nbtag_btagUp, "nbtag_btagUp/I");
   tree->Branch("nbtag_btagDown", &nbtag_btagDown, "nbtag_btagDown/I");
   tree->Branch("nbtag_noSF", &nbtag_noSF, "nbtag_noSF/I");
   tree->Branch("bpt_1",   &bpt,   "bpt/D");
   tree->Branch("beta_1",  &beta_1,  "beta_1/D");
   tree->Branch("bphi",  &bphi,  "bphi/F");
   
   tree->Branch("nuPx",&nuPx,"nuPx/F");
   tree->Branch("nuPy",&nuPy,"nuPy/F");
   tree->Branch("nuPz",&nuPz,"nuPz/F");
   tree->Branch("nuPt",&nuPt,"nuPt/F");
   tree->Branch("nuPhi",&nuPhi,"nuPhi/F");
    
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

   // MELA outputs
   // 1. Matrix element variables for different hypotheses (VBF Higgs, ggH + 2 jets, Z + 2 jets)
   tree->Branch("ME_ggh", &ME_ggh, "ME_ggh/F");
   tree->Branch("ME_vbf", &ME_vbf, "ME_vbf/F");
   tree->Branch("ME_z2j_1", &ME_z2j_1, "ME_z2j_1/F");
   tree->Branch("ME_z2j_2", &ME_z2j_2, "ME_z2j_2/F");

   // 2. Energy transfer (Q^2) variables
   tree->Branch("ME_q2v1", &ME_q2v1, "ME_q2v1/F");
   tree->Branch("ME_q2v2", &ME_q2v2, "ME_q2v2/F");
   
   // 3. Angle variables
   tree->Branch("ME_costheta1", &ME_costheta1, "ME_costheta1/F");
   tree->Branch("ME_costheta2", &ME_costheta2, "ME_costheta2/F");
   tree->Branch("ME_phi", &ME_phi, "ME_phi/F");
   tree->Branch("ME_costhetastar", &ME_costhetastar, "ME_costhetastar/F");
   tree->Branch("ME_phi1", &ME_phi1, "ME_phi1/F");
   
   // 4. Main BG vs. Higgs discriminators
   tree->Branch("ME_vbf_vs_Z", &ME_vbf_vs_Z, "ME_vbf_vs_Z/F");
   tree->Branch("ME_ggh_vs_Z", &ME_ggh_vs_Z, "ME_ggh_vs_Z/F");
   tree->Branch("ME_vbf_vs_ggh", &ME_vbf_vs_ggh, "ME_vbf_vs_ggh/F");
}

void SetDefaultValues(){

   isZLL = false;
   isZEE = false;
   isZMM = false;
   isZTT = false;
   veto_embedded = false;
   isPromptZMM = false;
   isPromptZEE = false;
   isZTTMM = false;
   isZTTEE = false;
   isZTTEM = false;

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
   
   THU_qqH_TOT = 1.0;
   THU_qqH_PTH200 = 1.0;
   THU_qqH_Mjj60 = 1.0;
   THU_qqH_Mjj120 = 1.0;
   THU_qqH_Mjj350 = 1.0;
   THU_qqH_Mjj700 = 1.0;
   THU_qqH_Mjj1000 = 1.0;
   THU_qqH_Mjj1500 = 1.0;
   THU_qqH_25 = 1.0;
   THU_qqH_JET01 = 1.0;

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
   qcdweight_2jet_rate_up =  1;
   qcdweight_2jet_rate_down =  1;
   qcdweight_0jet_shape_up =  1;
   qcdweight_0jet_shape_down =  1;
   qcdweight_1jet_shape_up =  1;
   qcdweight_1jet_shape_down =  1;
   qcdweight_2jet_shape_up =  1;
   qcdweight_2jet_shape_down =  1;

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

   prefiringweight     = 1.;
   prefiringweightup   = 1.;
   prefiringweightdown = 1.;

   nuPx = 0;
   nuPy = 0;
   nuPz = 0;
   nuPt = 0;
   nuPhi = 0;
   
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
   htxs_stage1p1cat = -1.;
   
   metFilters_ = true;
   badChargedCandidateFilter_ = true;
   badPFMuonFilter_ = true;
   badMuonFilter_ = true;
   duplicateMuonFilter_ = true;
   
   topPt = -1;
   antitopPt = -1;
   
   isGSfound = false;
   promptTausFirstCopy.clear();
   promptTausLastCopy.clear();
   promptElectrons.clear();
   promptMuons.clear();
   promptNeutrinos.clear();
   tauNeutrinos.clear();
   
   promptTausLV.SetXYZT(0.001,0.001,0,0);
   promptVisTausLV.SetXYZT(0.001,0.001,0,0);
   zBosonLV.SetXYZT(0,0,0,0);
   wBosonLV.SetXYZT(0,0,0,0);
   hBosonLV.SetXYZT(0,0,0,0);
   promptElectronsLV.SetXYZT(0.001,0.001,0,0);
   promptMuonsLV.SetXYZT(0.001,0.001,0,0);
   promptNeutrinosLV.SetXYZT(0,0,0,0);
   tauNeutrinosLV.SetXYZT(0,0,0,0);
   wDecayProductsLV.SetXYZT(0,0,0,0);
   fullVLV.SetXYZT(0,0,0,0);
   visVLV.SetXYZT(0,0,0,0);
   
   genBosonLV.SetXYZT(0,0,0,0);
   genVisBosonLV.SetXYZT(0,0,0,0);

   nLowPtLegElectron = 0;
   isLowPtLegElectron = false;
         
   nHighPtLegElectron = 0;
   isHighPtLegElectron = false;
   
   nLowPtLegMuon = 0;
   isLowPtLegMuon = false;
   
   nHighPtLegMuon = 0;
   isHighPtLegMuon = false;

   nMu23Ele12DzFilter = 0;
   isMu23Ele12DzFilter = false;
     
   nMu8Ele23DzFilter = 0;
   isMu8Ele23DzFilter = false;

   isMuon23matched = false;
   isMuon8matched  = false;
   isElectron23matched = false;
   isElectron12matched = false;
   isMu23 = false;
   isMu8 = false;
   isMu23dz = false;
   isMu8dz  = false;
   isEle23 = false;
   isEle12 = false;
   isEle23dz = false;
   isEle12dz = false;

   isoweight_1 = 1;
   isoweight_2 = 1;
   trigweight = 1;
   effweight = 1;

   Ele23EffData = 1;
   Ele12EffData = 1;
   Mu23EffData = 1;
   Mu8EffData = 1;
   
   Ele23EffMC = 1;
   Ele12EffMC = 1;
   Mu23EffMC = 1;
   Mu8EffMC = 1;
   
   jets.clear();
   jetspt20.clear();
   bjets.clear();
   bjets_mistagUp.clear();
   bjets_mistagDown.clear();
   bjets_btagUp.clear();
   bjets_btagDown.clear();
   bjets_nocleaned.clear();
   bjetsRaw.clear();
   
   njets = 0;
   njets_jecUncEta0To5Up   = 0;
   njets_jecUncEta0To5Down = 0;
   njets_jecUncEta0To3Up   = 0;
   njets_jecUncEta0To3Down = 0;
   njets_jecUncEta3To5Up   = 0;
   njets_jecUncEta3To5Down = 0;
   njets_jecUncRelativeBalUp   = 0;
   njets_jecUncRelativeBalDown = 0;
   njets_jecUncRelativeSampleYearUp   = 0;
   njets_jecUncRelativeSampleYearDown = 0;
   njets_jecUncEC2Up   = 0;
   njets_jecUncEC2Down = 0;
   njets_jecUncFlavorQCDUp   = 0;
   njets_jecUncFlavorQCDDown = 0;
   njets_jecUncEC2YearUp   = 0;
   njets_jecUncEC2YearDown = 0;
   njets_jecUncHFYearUp   = 0;
   njets_jecUncHFYearDown = 0;
   njets_jecUncAbsoluteYearUp   = 0;
   njets_jecUncAbsoluteYearDown = 0;
   njets_jecUncBBEC1YearUp   = 0;
   njets_jecUncBBEC1YearDown = 0;

   indexLeadingJet = -1;
   ptLeadingJet = -1;
   
   indexSubLeadingJet = -1;
   ptSubLeadingJet = -1;
   
   indexLeadingBJet = -1;
   ptLeadingBJet = -1;

   bpt = -10;
   beta_1 = -10;
   bphi = -10;
            
   jpt_1 = -10;
   jeta_1 = -10;
   jphi_1 = -10;
   jptraw_1 = -10;
   
   jpt_2 = -10;
   jeta_2 = -10;
   jphi_2 = -10;
   jptraw_2 = -10;
   mjj =  -10;
   dijetpt = -10;
   dijetphi = -10;
   jdeta =  -10;

   m_sv   = -10;
   mt_sv  = -10;
   pt_sv  = -10;
   eta_sv = -10;
   phi_sv = -10;
        
   ME_vbf= -10;
   ME_q2v1 = -10;
   ME_q2v2 = -10;
   ME_costheta1 = -10;
   ME_costheta2 = -10;
   ME_phi = -10;
   ME_costhetastar = -10;
   ME_phi1 = -10;
   ME_z2j_1 = -10;
   ME_z2j_2 = -10;
   ME_vbf_vs_Z = -10;
   ME_ggh_vs_Z = -10;
   ME_vbf_vs_ggh = -10;
   ME_ggh  = -10;

}
