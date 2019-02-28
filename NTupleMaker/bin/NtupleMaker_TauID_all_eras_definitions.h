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
#include "TF1.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"


 // first argument - config file 
  // second argument - filelist


using namespace std;


// histograms  
// number of input events
// and sum of eventweights

Float_t mueffweight =1.;
Float_t mutrigweight =1;

int nFiles = 0;
int nEvents = 0;
int WJetEvents = 0;
int WProdEvents = 0;
int WMuNuEvents = 0;
int WTauNuEvents = 0;
int DiJetEvents = 0;
int SingleJetEvents = 0;
int JetTauEvents = 0;
int TrigEvents = 0;

int nTotalFiles = 0;


// booleans 
bool isWJet = false;
bool isWTauNu = false;
bool isWMuNu  = false;
bool isDiJet = false;


UInt_t run_;
UInt_t lumi_;
UInt_t event_;

Float_t puWeight_;
Float_t genWeight_;
Float_t trigWeight_;
Float_t weight_;

Bool_t trigger_;
Bool_t isWTrig_;
Bool_t isZTrig_;
Float_t metNoSelMu_;
Float_t mhtNoSelMu_;
UInt_t nMuonTrig_;
UInt_t nSelMuonTrig_;
Float_t dPhiMetMuon_;

UInt_t nVert_;
Float_t genHt_;

Bool_t  trig_;

Bool_t metFilters_;

Float_t wMass_;
Float_t wPt_;
Float_t wEta_;
Float_t wPhi_;
Int_t   wDecay_;
Int_t   wTauDecay_;
Float_t wCharge_;

Float_t lepWPt_;
Float_t lepWEta_;
Float_t lepWPhi_;
Float_t lepWE_;

Float_t nuWPt_;
Float_t nuWEta_;
Float_t nuWPhi_;

Float_t tauMatchJetPt_;
Float_t tauMatchJetEta_;
Float_t tauMatchJetPhi_;

Float_t met_;
Float_t metphi_;
Float_t mttau_;
Float_t mtgen_;
Float_t mtmuon_;

Float_t muonPt_;
Float_t muonPz_;
Float_t muonEta_;
Float_t muonPhi_;
Int_t   muonQ_;

Float_t muon2Pt_;
Float_t muon2Eta_;
Float_t muon2Phi_;
Int_t   muon2Q_;

Float_t recoilM_;
Float_t recoilPt_;
Float_t recoilEta_;
Float_t recoilPhi_;

Float_t tauPt_;
Float_t tauEta_;
Float_t tauPhi_;
Float_t tauMass_;
Int_t   tauQ_;
Float_t tauPz_;

Float_t genTauWPt_;
Float_t genTauWEta_;
Float_t genTauWPhi_;
Float_t genTauWE_;

Float_t tauJetPt_;
Float_t tauJetEta_;
Float_t tauJetPhi_;
Bool_t  tauJetTightId_;
Int_t tauJetFlavor_;

Float_t recoilDPhi_;

Float_t recoilJetRatio_;
Float_t recoilJetDPhi_;

Int_t   tauDecay_;
Int_t   tauGenDecay_;
Int_t   tauGenMatchDecay_;
UInt_t  tauNtrk05_;
UInt_t  tauNtrk08_;
UInt_t  tauNtrk1_;

UInt_t  tauGenMatch_;

Bool_t  tauDM_;
Bool_t  tauNewDM_;

Float_t tauIso_;
Bool_t  tauLooseIso_;
Bool_t  tauMediumIso_;
Bool_t  tauTightIso_;

Bool_t  tauVLooseMvaIso_;
Bool_t  tauLooseMvaIso_;
Bool_t  tauMediumMvaIso_;
Bool_t  tauTightMvaIso_;
Bool_t  tauVTightMvaIso_;
Bool_t  tauVVTightMvaIso_;

Bool_t  tauVVLooseMva2017v2Iso_;
Bool_t  tauVLooseMva2017v2Iso_;
Bool_t  tauLooseMva2017v2Iso_;
Bool_t  tauMediumMva2017v2Iso_;
Bool_t  tauTightMva2017v2Iso_;
Bool_t  tauVTightMva2017v2Iso_;
Bool_t  tauVVTightMva2017v2Iso_;

Bool_t tauAntiMuonLoose3_;
Bool_t tauAntiMuonTight3_;

Bool_t tauAntiElectronVLooseMVA6_;
Bool_t tauAntiElectronLooseMVA6_;
Bool_t tauAntiElectronMediumMVA6_;
Bool_t tauAntiElectronTightMVA6_;
Bool_t tauAntiElectronVTightMVA6_;

Bool_t tauSinglePFTau180Trk50_;
Bool_t tauSinglePFTau180Trk50oneprong_;

Float_t tauLeadingTrackPt_;
Float_t tauLeadingTrackEta_;
Float_t tauLeadingTrackPhi_;
Float_t tauLeadingTrackDz_;
Float_t tauLeadingTrackDxy_;

UInt_t nMuon_;
UInt_t nSelMuon_;
UInt_t nElec_;

UInt_t nJetsCentral20_;
UInt_t nJetsCentral30_;

UInt_t nJetsForward20_;
UInt_t nJetsForward30_;

Float_t jetPt_;
Float_t jetEta_;
Float_t jetPhi_;
Float_t jetBtag_;

UInt_t jetChargedMult_;
UInt_t jetNeutralMult_;
UInt_t jetChargedHadMult_;
Float_t jetNeutralEMEnergyFraction_;
Float_t jetNeutralHadEnergyFraction_;

Float_t jet2Pt_;
Float_t jet2Eta_;
Float_t jet2Phi_;
Float_t jet2Btag_;

UInt_t jet2ChargedMult_;
UInt_t jet2NeutralMult_;
UInt_t jet2ChargedHadMult_;
Float_t jet2NeutralEMEnergyFraction_;
Float_t jet2NeutralHadEnergyFraction_;

UInt_t nSelTaus_;
UInt_t nTaus20_;
UInt_t nTaus30_;

Float_t JetHt_            ; // sumJetPtCentral30 + sumJetPtForward30
Float_t SoftJetHt_        ; // sumJetPtCentral20 + sumJetPtForward30
Float_t Ht_               ; // sumJetPtCentral30 + sumJetPtForward30 + sumLeptonPt
Float_t SoftHt_           ; // sumJetPtCentral20 + sumJetPtForward30 + sumLeptonPt
Float_t HtNoRecoil_       ; // sumJetPtCentral30 + sumJetPtForward30 + sumLeptonPt - sumPtRecoil
Float_t SoftHtNoRecoil_   ; // sumJetPtCentral20 + sumJetPtForward30 + sumLeptonPt - sumPtRecoil 
Float_t mht_;
Float_t mhtNoMu_;
Float_t metNoMu_;
Float_t dPhiMetTau_;

Int_t selection_; 
UInt_t npartons_; 
UInt_t npartonsNLO_;
Float_t lheWPt_;
//  0 : Z->mumu+Jet, 
//  1 : W->muv+Jet
//  2 : W*->muv 
//  3 : W*->tauv
//  4 : dijet
// 10 : W->mu+v 
// 11 : single jet + MEt
// 12 : dijet sample

Bool_t pfJet60_;
Bool_t pfJet80_;
Bool_t pfJet140_;
Bool_t pfJet200_;
Bool_t pfJet260_;
Bool_t pfJet320_;
Bool_t pfJet400_;
Bool_t pfJet450_;

Bool_t pf2Jet60_;
Bool_t pf2Jet80_;
Bool_t pf2Jet140_;
Bool_t pf2Jet200_;

UInt_t nJets20_;
Float_t jet20Pt_[10];
Float_t jet20Eta_[10];
Float_t jet20Phi_[10];

unsigned int nSingleMuonHLTFilter = 0;
bool isSingleMuonHLTFilter = false;

unsigned int nPFJet60HLTFilter = 0;
bool isPFJet60HLTFilter = false;

unsigned int nPFJet80HLTFilter = 0;
bool isPFJet80HLTFilter = false;

unsigned int nPFJet140HLTFilter = 0;
bool isPFJet140HLTFilter = false;

unsigned int nPFJet200HLTFilter = 0;
bool isPFJet200HLTFilter = false;

unsigned int nPFJet260HLTFilter = 0;
bool isPFJet260HLTFilter = false;

unsigned int nPFJet320HLTFilter = 0;
bool isPFJet320HLTFilter = false;

unsigned int nPFJet400HLTFilter = 0;
bool isPFJet400HLTFilter = false;

unsigned int nPFJet450HLTFilter = 0;
bool isPFJet450HLTFilter = false;

unsigned int nSinglePFTau180Trk50Filter = 0;
bool isSinglePFTau180Trk50Filter = false;

unsigned int nSinglePFTau180Trk50oneprongFilter = 0;
bool isSinglePFTau180Trk50oneprongFilter = false;

// unsigned int nSinglePFTau180Trk50Filter = 0;
// bool isSinglePFTau180Trk50Filter = false;

// unsigned int nSinglePFTau180Trk50oneprongFilter = 0;
// bool isSinglePFTau180Trk50oneprongFilter = false;

TTree * ntuple_ = new TTree("NTuple","NTuple");
TTree * trigNTuple_ = new TTree("TriggerNTuple","TriggerNTuple");


void SetupTrees()
{
     
   ntuple_->Branch("event",&event_,"event/i");
   ntuple_->Branch("run",  &run_,  "run/i");
   ntuple_->Branch("luminosityBlock", &lumi_,  "luminosityBlock/i");
   
   ntuple_->Branch("puWeight",  &puWeight_,  "puWeight/F");
   ntuple_->Branch("genWeight", &genWeight_, "genWeight/F");
   ntuple_->Branch("trigWeight",&trigWeight_,"trigWeight/F");
   ntuple_->Branch("weight",    &weight_,    "weight/F");
   ntuple_->Branch("mueffweight", &mueffweight, "mueffweight/F");
   ntuple_->Branch("mutrigweight",&mutrigweight,"mutrigweight/F");
   
   ntuple_->Branch("NVert",&nVert_,"NVert/i");
   ntuple_->Branch("trigger",&trig_,"trigger/O");
   ntuple_->Branch("metFilters",&metFilters_,"metFilters/O");
   
   ntuple_->Branch("WMass",&wMass_,   "WMass/F");
   ntuple_->Branch("WPt",  &wPt_,     "WPt/F");
   ntuple_->Branch("WEta", &wEta_,    "WEta/F");
   ntuple_->Branch("WPhi", &wPhi_,    "WPhi/F");
   ntuple_->Branch("WDecay", &wDecay_,"WDecay/I");
   ntuple_->Branch("WCharge", &wCharge_, "WCharge/F");
   ntuple_->Branch("WTauDecay",&wTauDecay_, "WTauDecay/I");
   
   ntuple_->Branch("lepWPt",&lepWPt_,"lepWPt/F");
   ntuple_->Branch("lepWEta",&lepWEta_,"lepWEta/F");
   ntuple_->Branch("lepWPhi",&lepWPhi_,"lepWPhi/F");
   ntuple_->Branch("lepWE",&lepWE_,"lepWE/F");
   
   ntuple_->Branch("nuWPt",&nuWPt_,"nuWPt/F");
   ntuple_->Branch("nuWEta",&nuWEta_,"nuWEta/F");
   ntuple_->Branch("nuWPhi",&nuWPhi_,"nuWPhi/F");

   ntuple_->Branch("met",    &met_,   "met/F");
   ntuple_->Branch("metphi", &metphi_,"metphi/F");
   ntuple_->Branch("mttau",  &mttau_, "mttau/F");
   ntuple_->Branch("mtgen",  &mtgen_, "mtgen/F");
   ntuple_->Branch("mtmuon", &mtmuon_,"mtmuon/F");
   ntuple_->Branch("genHt",&genHt_,"genHt/F");
   
   ntuple_->Branch("muonPt",  &muonPt_,  "muonPt/F");
   ntuple_->Branch("muonPz",  &muonPz_,  "muonPz/F");
   ntuple_->Branch("muonEta", &muonEta_, "muonEta/F");
   ntuple_->Branch("muonPhi", &muonPhi_, "muonPhi/F");
   ntuple_->Branch("muonQ",   &muonQ_,   "muonQ/I");
   
   ntuple_->Branch("muon2Pt",  &muon2Pt_,  "muon2Pt/F");
   ntuple_->Branch("muon2Eta", &muon2Eta_, "muon2Eta/F");
   ntuple_->Branch("muon2Phi", &muon2Phi_, "muon2Phi/F");
   ntuple_->Branch("muon2Q",   &muon2Q_,   "muon2Q/I");
   
   ntuple_->Branch("tauPt",  &tauPt_,  "tauPt/F");
   ntuple_->Branch("tauPz",  &tauPz_,  "tauPz/F");
   ntuple_->Branch("tauEta", &tauEta_, "tauEta/F");
   ntuple_->Branch("tauPhi", &tauPhi_, "tauPhi/F");
   ntuple_->Branch("tauMass",&tauMass_,"tauMass/F");
   ntuple_->Branch("tauQ",   &tauQ_,   "tauQ/I");
   
   ntuple_->Branch("genTauWPt",  &genTauWPt_,  "genTauWPt/F");
   ntuple_->Branch("genTauWEta", &genTauWEta_, "genTauWEta/F");
   ntuple_->Branch("genTauWPhi", &genTauWPhi_, "genTauWPhi/F");
   ntuple_->Branch("genTauWE", &genTauWE_, "genTauWE/F");
   
   ntuple_->Branch("tauJetPt",  &tauJetPt_,  "tauJetPt/F");
   ntuple_->Branch("tauJetEta", &tauJetEta_, "tauJetEta/F");
   ntuple_->Branch("tauJetPhi", &tauJetPhi_, "tauJetPhi/F");
   ntuple_->Branch("tauJetTightId", &tauJetTightId_, "tauJetTightId/O");
   ntuple_->Branch("tauJetFlavor",&tauJetFlavor_,"tauJetFlavor/I");
   
   ntuple_->Branch("tauLeadingTrackPt",&tauLeadingTrackPt_,"tauLeadingTrackPt/F");
   ntuple_->Branch("tauLeadingTrackEta",&tauLeadingTrackEta_,"tauLeadingTrackEta/F");
   ntuple_->Branch("tauLeadingTrackPhi",&tauLeadingTrackPhi_,"tauLeadingTrackPhi/F");
   ntuple_->Branch("tauLeadingTrackDz",&tauLeadingTrackDz_,"tauLeadingTrackDz/F");
   ntuple_->Branch("tauLeadingTrackDxy",&tauLeadingTrackDxy_,"tauLeadingTrackDxy/F");
   
   ntuple_->Branch("recoilDPhi",&recoilDPhi_,"recoilDPhi/F");
   
   ntuple_->Branch("recoilJetRatio",&recoilJetRatio_,"recoilJetRatio/F");
   ntuple_->Branch("recoilJetDPhi",&recoilJetDPhi_,"recoilJetDPhi/F");
   
   ntuple_->Branch("recoilM",&recoilM_,"recoilM/F");
   ntuple_->Branch("recoilPt",&recoilPt_,"recoilPt/F");
   ntuple_->Branch("recoilEta",&recoilEta_,"recoilEta/F");
   ntuple_->Branch("recoilPhi",&recoilPhi_,"recoilPhi/F");
   
   ntuple_->Branch("tauDecay",   &tauDecay_,   "tauDecay/I");
   ntuple_->Branch("tauGenDecay",&tauGenDecay_,"tauGenDecay/I");
   ntuple_->Branch("tauGenMatchDecay",&tauGenMatchDecay_,"tauGenMatchDecay/I");
   ntuple_->Branch("tauGenMatch",&tauGenMatch_,"tauGenMatch/i");
   
   ntuple_->Branch("tauNtrk1", &tauNtrk1_, "tauNtrk1/i");
   ntuple_->Branch("tauNtrk08",&tauNtrk08_,"tauNtrk08/i");
   ntuple_->Branch("tauNtrk05",&tauNtrk05_,"tauNtrk05/i");
   
   ntuple_->Branch("tauDM",&tauDM_,"tauDM/O");
   ntuple_->Branch("tauNewDM",&tauNewDM_,"tauNewDM/O");
   
   ntuple_->Branch("tauIso", &tauIso_, "tauIso/F");
   ntuple_->Branch("tauLooseIso", &tauLooseIso_, "tauLooseIso/O");
   ntuple_->Branch("tauMediumIso",&tauMediumIso_,"tauMediumIso/O");
   ntuple_->Branch("tauTightIso", &tauTightIso_, "tauTightIso/O");
   
   ntuple_->Branch("tauVLooseMvaIso", &tauVLooseMvaIso_, "tauVLooseMvaIso/O");
   ntuple_->Branch("tauLooseMvaIso", &tauLooseMvaIso_, "tauLooseMvaIso/O");
   ntuple_->Branch("tauMediumMvaIso",&tauMediumMvaIso_,"tauMediumMvaIso/O");
   ntuple_->Branch("tauTightMvaIso", &tauTightMvaIso_, "tauTightMvaIso/O");
   ntuple_->Branch("tauVTightMvaIso", &tauVTightMvaIso_, "tauVTightMvaIso/O");
   ntuple_->Branch("tauVVTightMvaIso", &tauVVTightMvaIso_, "tauVVTightMvaIso/O");
   
   ntuple_->Branch("tauVVLooseMva2017v2Iso", &tauVVLooseMva2017v2Iso_, "tauVVLooseMva2017v2Iso/O");
   ntuple_->Branch("tauVLooseMva2017v2Iso", &tauVLooseMva2017v2Iso_, "tauVLooseMva2017v2Iso/O");
   ntuple_->Branch("tauLooseMva2017v2Iso", &tauLooseMva2017v2Iso_, "tauLooseMva2017v2Iso/O");
   ntuple_->Branch("tauMediumMva2017v2Iso",&tauMediumMva2017v2Iso_,"tauMediumMva2017v2Iso/O");
   ntuple_->Branch("tauTightMva2017v2Iso", &tauTightMva2017v2Iso_, "tauTightMva2017v2Iso/O");
   ntuple_->Branch("tauVTightMva2017v2Iso", &tauVTightMva2017v2Iso_, "tauVTightMva2017v2Iso/O");
   ntuple_->Branch("tauVVTightMva2017v2Iso", &tauVVTightMva2017v2Iso_, "tauVVTightMva2017v2Iso/O");
   
   ntuple_->Branch("tauAntiMuonLoose3",&tauAntiMuonLoose3_,"tauAntiMuonLoose3/O");
   ntuple_->Branch("tauAntiMuonTight3",&tauAntiMuonTight3_,"tauAntiMuonTight3/O");
   
   ntuple_->Branch("tauAntiElectronVLooseMVA6",&tauAntiElectronVLooseMVA6_,"tauAntiElectronVLooseMVA6/O");
   ntuple_->Branch("tauAntiElectronLooseMVA6", &tauAntiElectronLooseMVA6_, "tauAntiElectronLooseMVA6/O");
   ntuple_->Branch("tauAntiElectronMediumMVA6", &tauAntiElectronMediumMVA6_, "tauAntiElectronMediumMVA6/O");
   ntuple_->Branch("tauAntiElectronTightMVA6",&tauAntiElectronTightMVA6_,"tauAntiElectronTightMVA6/O");
   ntuple_->Branch("tauAntiElectronVTightMVA6", &tauAntiElectronVTightMVA6_, "tauAntiElectronVTightMVA6/O");
   
   ntuple_->Branch("tauSinglePFTau180Trk50",&tauSinglePFTau180Trk50_,"tauSinglePFTau180Trk50/O");
   ntuple_->Branch("tauSinglePFTau180Trk50oneprong",&tauSinglePFTau180Trk50oneprong_,"tauSinglePFTau180Trk50oneprong/O");
   
   
   ntuple_->Branch("nMuon",&nMuon_,"nMuon/i");
   ntuple_->Branch("nSelMuon",&nSelMuon_,"nSelMuon/i");
   ntuple_->Branch("nElec",&nElec_,"nElec/i");
   
   ntuple_->Branch("nJetsCentral20",&nJetsCentral20_,"nJetsCentral20/i");
   ntuple_->Branch("nJetsCentral30",&nJetsCentral30_,"nJetsCentral30/i");
   
   ntuple_->Branch("nJetsForward20",&nJetsForward20_,"nJetsForward20/i");
   ntuple_->Branch("nJetsForward30",&nJetsForward30_,"nJetsForward30/i");
   
   ntuple_->Branch("jetPt", &jetPt_, "jetPt/F");
   ntuple_->Branch("jetEta",&jetEta_,"jetEta/F");
   ntuple_->Branch("jetPhi",&jetPhi_,"jetPhi/F");
   ntuple_->Branch("jetBtag",&jetBtag_,"jetBtag/F");
   
   ntuple_->Branch("jetChargedMult",   &jetChargedMult_,   "jetChargedMult/i");
   ntuple_->Branch("jetNeutralMult",   &jetNeutralMult_,   "jetNeutralMult/i");
   ntuple_->Branch("jetChargedHadMult",&jetChargedHadMult_,"jetChargedHadMult/i");
   
   ntuple_->Branch("jetNeutralEMEnergyFraction", &jetNeutralEMEnergyFraction_, "jetNeutralEMEnergyFraction/F");
   ntuple_->Branch("jetNeutralHadEnergyFraction",&jetNeutralHadEnergyFraction_,"jetNeutralHadEnergyFraction/F");
   
   ntuple_->Branch("jet2Pt", &jet2Pt_, "jet2Pt/F");
   ntuple_->Branch("jet2Eta",&jet2Eta_,"jet2Eta/F");
   ntuple_->Branch("jet2Phi",&jet2Phi_,"jet2Phi/F");
   ntuple_->Branch("jet2Btag",&jet2Btag_,"jet2Btag/F");
   
   ntuple_->Branch("jet2ChargedMult",   &jet2ChargedMult_,   "jet2ChargedMult/i");
   ntuple_->Branch("jet2NeutralMult",   &jet2NeutralMult_,   "jet2NeutralMult/i");
   ntuple_->Branch("jet2ChargedHadMult",&jet2ChargedHadMult_,"jet2ChargedHadMult/i");
   
   ntuple_->Branch("jet2NeutralEMEnergyFraction", &jet2NeutralEMEnergyFraction_, "jet2NeutralEMEnergyFraction/F");
   ntuple_->Branch("jet2NeutralHadEnergyFraction",&jet2NeutralHadEnergyFraction_,"jet2NeutralHadEnergyFraction/F");
   
   ntuple_->Branch("nSelTaus",&nSelTaus_,"nSelTaus/i");
   ntuple_->Branch("nTaus20",&nTaus20_,"nTaus20/i");
   ntuple_->Branch("nTaus30",&nTaus30_,"nTaus30/i");
   
   ntuple_->Branch("JetHt",&JetHt_,"JetHt/F");
   ntuple_->Branch("SoftJetHt",&SoftJetHt_,"SoftJetHt/F");
   ntuple_->Branch("Ht",&Ht_,"Ht/F");
   ntuple_->Branch("SoftHt",&SoftHt_,"SoftHt/F");
   ntuple_->Branch("HtNoRecoil",&HtNoRecoil_,"HtNoRecoil/F");
   ntuple_->Branch("SoftHtNoRecoil",&SoftHtNoRecoil_,"SoftHtNoRecoil/F");
   ntuple_->Branch("mht",&mht_,"mht/F");
   ntuple_->Branch("mhtNoMu",&mhtNoMu_,"mhtNoMu/F");
   ntuple_->Branch("metNoMu",&metNoMu_,"metNoMu/F");
   ntuple_->Branch("dPhiMetTau",&dPhiMetTau_,"dPhiMetTau/F");
   
   ntuple_->Branch("pfJet60",&pfJet60_,"pfJet60/O");
   ntuple_->Branch("pfJet80",&pfJet80_,"pfJet80/O");
   ntuple_->Branch("pfJet140",&pfJet140_,"pfJet140/O");
   ntuple_->Branch("pfJet200",&pfJet200_,"pfJet200/O");
   ntuple_->Branch("pfJet260",&pfJet260_,"pfJet260/O");
   ntuple_->Branch("pfJet320",&pfJet320_,"pfJet320/O");
   ntuple_->Branch("pfJet400",&pfJet400_,"pfJet400/O");
   ntuple_->Branch("pfJet450",&pfJet450_,"pfJet450/O");
   
   ntuple_->Branch("pf2Jet60",&pf2Jet60_,"pf2Jet60/O");
   ntuple_->Branch("pf2Jet80",&pf2Jet80_,"pf2Jet80/O");
   ntuple_->Branch("pf2Jet140",&pf2Jet140_,"p2fJet140/O");
   ntuple_->Branch("pf2Jet200",&pf2Jet200_,"p2fJet200/O");
   
   ntuple_->Branch("Selection",&selection_,"Selection/I");
   
   ntuple_->Branch("npartons",&npartons_,"npartons/i");
   ntuple_->Branch("npartonsNLO",&npartonsNLO_,"npartonsNLO/i");
   ntuple_->Branch("lheWPt",&lheWPt_,"lheWPt/F");
   
   
   trigNTuple_->Branch("event",&event_,"event/i");
   trigNTuple_->Branch("run",  &run_,  "run/i");
   trigNTuple_->Branch("luminosityBlock", &lumi_,  "luminosityBlock/i");
   trigNTuple_->Branch("trigger",&trigger_,"trigger/O");
   trigNTuple_->Branch("NVert",&nVert_,"NVert/i");
   trigNTuple_->Branch("metNoMu",&metNoMu_,"metNoMu/F");
   trigNTuple_->Branch("mhtNoMu",&mhtNoMu_,"mhtNoMu/F");
   trigNTuple_->Branch("metNoSelMu",&metNoSelMu_,"metNoSelMu/F");
   trigNTuple_->Branch("mhtNoSelMu",&mhtNoSelMu_,"mhtNoSelMu/F");
   trigNTuple_->Branch("IsW",&isWTrig_,"IsW/O");
   trigNTuple_->Branch("IsZ",&isZTrig_,"IsZ/O");
   trigNTuple_->Branch("nMuon",&nMuonTrig_,"nMuon/i");
   trigNTuple_->Branch("nSelMuon",&nSelMuonTrig_,"nSelMuon/i");
   trigNTuple_->Branch("npartons",&npartons_,"npartons/i");
   trigNTuple_->Branch("npartonsNLO",&npartonsNLO_,"npartonsNLO/i");
   trigNTuple_->Branch("lheWPt",&lheWPt_,"lheWPt/F");
   trigNTuple_->Branch("met",&met_,"met/F");
   trigNTuple_->Branch("mht",&mht_,"mht/F");
   trigNTuple_->Branch("mtmuon",&mtmuon_,"mtmuon/F");
   trigNTuple_->Branch("muonPt",&muonPt_,"muonPt/F");
   trigNTuple_->Branch("muonEta",&muonEta_,"muonEta/F");
   trigNTuple_->Branch("dPhiMetMuon",&dPhiMetMuon_,"dPhiMetMuon/F");
   trigNTuple_->Branch("WMass",&wMass_,   "WMass/F");
   trigNTuple_->Branch("puWeight",  &puWeight_,  "puWeight/F");
   trigNTuple_->Branch("genWeight", &genWeight_, "genWeight/F");
   trigNTuple_->Branch("metFilters",&metFilters_,"metFilters/O");
   trigNTuple_->Branch("nJetsCentral20",&nJetsCentral20_,"nJetsCentral20/i");
   trigNTuple_->Branch("nJetsCentral30",&nJetsCentral30_,"nJetsCentral30/i");
   trigNTuple_->Branch("nJetsForward20",&nJetsForward20_,"nJetsForward20/i");
   trigNTuple_->Branch("nJetsForward30",&nJetsForward30_,"nJetsForward30/i");
   trigNTuple_->Branch("nElec",&nElec_,"nElec/i");
   trigNTuple_->Branch("nSelTaus",&nSelTaus_,"nSelTaus/i");
    
}



void SetDefaultValues(){

   weight_ = 1;
   genWeight_ = 1;
   trigWeight_ = 1;
   puWeight_ = 1;
   
   trig_ = false;
   metFilters_ = true;

   wMass_ = -1;
   wPt_ = -1;
   wEta_ = 0;
   wPhi_ = 0;
   wDecay_ = -1;
   wTauDecay_ = -1;
   wCharge_ = 0.0;
 
   nuWPt_ = -1;
   nuWEta_ = 0;
   nuWPhi_ = 0;

   lepWPt_ = -1;
   lepWEta_ = 0;
   lepWPhi_ = 0;
   lepWE_   = 0;

   met_ =  -1;
   metphi_ =  0;
   mttau_ = 0;
   mtgen_ = 0;
   mtmuon_ = 0;
   
   muonPt_ =  -1;
   muonPz_ =  -1;
   muonEta_ = 0;
   muonPhi_ = 0;
   muonQ_ = 0;
   dPhiMetMuon_ = 0;

   muon2Pt_ =  -1;
   muon2Eta_ = 0;
   muon2Phi_ = 0;
   muon2Q_ = 0;
   
   tauPt_ = 0;
   tauPz_ = 0;
   tauEta_ = 0;
   tauPhi_ = 0;
   tauMass_ = 0;
   tauQ_ = 0;

   tauJetPt_ = 0;
   tauJetEta_ = 0;
   tauJetPhi_ = 0;
   tauJetTightId_ = false;
   tauJetFlavor_ = 0;
   
   recoilDPhi_ = 0;
   
   recoilJetRatio_ = -1;
   recoilJetDPhi_ = 0;

   recoilM_ = -1;
   recoilPt_ = -1;
   recoilEta_ = 0;
   recoilPhi_ = 0;
   
   tauDecay_ = -1;
   tauGenDecay_ = -1;
   tauGenMatchDecay_ = -1;
   tauGenMatch_ = 6;

   tauNtrk1_  = 0;
   tauNtrk05_ = 0;
   tauNtrk08_ = 0;
   
   tauLeadingTrackPt_  = 0;
   tauLeadingTrackEta_ = 0;
   tauLeadingTrackPhi_ = 0;
   tauLeadingTrackDz_  = -999;
   tauLeadingTrackDxy_ = -999;
   
   tauDM_ = false;
   tauNewDM_ = false;

   tauIso_ = 0;
   tauLooseIso_ = false;
   tauMediumIso_ = false;
   tauTightIso_ = false;
   
   tauVLooseMvaIso_ = false;
   tauLooseMvaIso_ = false;
   tauMediumMvaIso_ = false;
   tauTightMvaIso_ = false;
   tauVTightMvaIso_ = false;
   tauVVTightMvaIso_ = false;
   
   tauVVLooseMva2017v2Iso_ = false;
   tauVLooseMva2017v2Iso_ = false;
   tauLooseMva2017v2Iso_ = false;
   tauMediumMva2017v2Iso_ = false;
   tauTightMva2017v2Iso_ = false;
   tauVTightMva2017v2Iso_ = false;
   tauVVTightMva2017v2Iso_ = false;
   
   tauAntiMuonLoose3_ = false;
   tauAntiMuonTight3_ = false;
   
   tauAntiElectronVLooseMVA6_ = false;
   tauAntiElectronLooseMVA6_ = false;
   tauAntiElectronMediumMVA6_ = false;
   tauAntiElectronTightMVA6_ = false;
   tauAntiElectronVTightMVA6_ = false;
   
   tauSinglePFTau180Trk50_ = false;
   tauSinglePFTau180Trk50oneprong_ = false;

   nMuon_ = 0;
   nSelMuon_ = 0;
   nElec_ = 0;
   
   nJetsCentral20_ = 0;
   nJetsCentral30_ = 0;
   
   nJetsForward20_ = 0;
   nJetsForward30_ = 0;

   jetPt_ = 0;
   jetEta_ = 0;
   jetPhi_ = 0;
   jetBtag_ = 0;
   
   jetChargedMult_ = 0;
   jetNeutralMult_ = 0;
   jetChargedHadMult_ = 0;
   
   jetNeutralEMEnergyFraction_ = 0;
   jetNeutralHadEnergyFraction_ = 0;
   
   jet2Pt_ = 0;
   jet2Eta_ = 0;
   jet2Phi_ = 0;
   jet2Btag_ = 0;
   
   jet2ChargedMult_ = 0;
   jet2NeutralMult_ = 0;
   jet2ChargedHadMult_ = 0;
   
   jet2NeutralEMEnergyFraction_ = 0;
   jet2NeutralHadEnergyFraction_ = 0;
  
   nJets20_ = 0;
   nTaus20_ = 0;
   nTaus30_ = 0;
   nSelTaus_ = 0;
   
   JetHt_ = 0;
   SoftJetHt_ = 0;
   Ht_ = 0;
   SoftHt_ = 0;
   
   HtNoRecoil_ = 0;
   SoftHtNoRecoil_ = 0;
   
   selection_ = -1;

   pfJet60_ = false;
   pfJet80_ = false;
   pfJet140_ = false;
   pfJet200_ = false;
   pfJet260_ = false;
   pfJet320_ = false;
   pfJet400_ = false;
   pfJet450_ = false;
   
   pf2Jet60_ = false;
   pf2Jet80_ = false;
   pf2Jet140_ = false;
   pf2Jet200_ = false;

   mueffweight =1.;
   mutrigweight =1;

   isWJet = false;
   isWTauNu = false;
   isWMuNu  = false;
   isDiJet = false;

   trigger_ = false;
   isWTrig_ = false;
   isZTrig_ = false;
   metNoMu_ = 0;
   dPhiMetTau_ = 0;
   mht_ = 0;
   mhtNoMu_ = 0;
   metNoSelMu_ = 0;
   mhtNoSelMu_ = 0;
   nMuonTrig_ = 0;
   nSelMuonTrig_ = 0;

   npartons_ = 9999;
   npartonsNLO_ = 9999;
   lheWPt_ = -1;

   nSingleMuonHLTFilter = 0;
   isSingleMuonHLTFilter = false;
   
   nPFJet60HLTFilter = 0;
   isPFJet60HLTFilter = false;
   
   nPFJet80HLTFilter = 0;
   isPFJet80HLTFilter = false;
   
   nPFJet140HLTFilter = 0;
   isPFJet140HLTFilter = false;
   
   nPFJet200HLTFilter = 0;
   isPFJet200HLTFilter = false;
   
   nPFJet260HLTFilter = 0;
   isPFJet260HLTFilter = false;
   
   nPFJet320HLTFilter = 0;
   isPFJet320HLTFilter = false;
   
   nPFJet400HLTFilter = 0;
   isPFJet400HLTFilter = false;
   
   nPFJet450HLTFilter = 0;
   isPFJet450HLTFilter = false;
   
   nSinglePFTau180Trk50Filter = 0;
   isSinglePFTau180Trk50Filter = false;
   
   nSinglePFTau180Trk50oneprongFilter = 0;
   isSinglePFTau180Trk50oneprongFilter = false;
   
   // nSinglePFTau180Trk50Filter = 0;
   // isSinglePFTau180Trk50Filter = false;

   // nSinglePFTau180Trk50oneprongFilter = 0;
   // isSinglePFTau180Trk50oneprongFilter = false;


}
