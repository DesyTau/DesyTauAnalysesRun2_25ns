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

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"

bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags) {

  bool passed = true;
  unsigned int nFlags = metFlags.size();
  //  std::cout << "MEt filters : " << std::endl;
  for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
    TString flagName(it->first);
    //    std::cout << it->first << " : " << it->second << std::endl;
    for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
      if (flagName.Contains(metFlags[iFilter])) {
	if (it->second==0) {
	  passed = false;
	  break;
	}
      }
    }
  }
  //  std::cout << "Passed : " << passed << std::endl;
  return passed;

}

double dPhiFromLV(TLorentzVector v1, TLorentzVector v2) {

  return dPhiFrom2P(v1.Px(),v1.Py(),v2.Px(),v2.Py());

}

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // general settings
  const bool debug = cfg.get<bool>("Debug");
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");
  
  // tau cuts
  const float ptTauCut  = cfg.get<float>("PtTauCut");
  const float etaTauCut = cfg.get<float>("EtaTauCut");

  // trigger information
  const string metHTLName        = cfg.get<string>("MetHLTName");
  const string singleMuonHLTName = cfg.get<string>("SingleMuonHLTName");
  const string singleMuonHLTFilterName = cfg.get<string>("SingleMuonHLTFilterName");
  const string pfJet40HLTFilterName = cfg.get<string>("PFJet40HLTFilterName"); 
  const string pfJet60HLTFilterName = cfg.get<string>("PFJet60HLTFilterName"); 
  const string pfJet80HLTFilterName = cfg.get<string>("PFJet80HLTFilterName"); 
  const string pfJet140HLTFilterName = cfg.get<string>("PFJet140HLTFilterName"); 

  TString MetHLTName(metHTLName);
  TString SingleMuonHLTName(singleMuonHLTName);
  TString SingleMuonHLTFilterName(singleMuonHLTFilterName);
  TString PFJet40HLTFilterName(pfJet40HLTFilterName);
  TString PFJet60HLTFilterName(pfJet60HLTFilterName);
  TString PFJet80HLTFilterName(pfJet80HLTFilterName);
  TString PFJet140HLTFilterName(pfJet140HLTFilterName);

  // muon selection
  const float ptMuCut       = cfg.get<float>("PtMuCut");
  const float etaMuCut      = cfg.get<float>("EtaMuCut");
  const float isoMuCut      = cfg.get<float>("IsoMuCut");
  const float dxyMuCut      = cfg.get<float>("dxyMuCut");  
  const float dzMuCut       = cfg.get<float>("dzMuCut");
  const float isoSelMuCut   = cfg.get<float>("IsoSelMuCut");
  const float ptSelMuCut    = cfg.get<float>("PtSelMuCut");
  const float ptTrigMuCut   = cfg.get<float>("PtTrigMuCut");
  const bool  isDRIso03 = cfg.get<bool>("IsDRIso03");
  const bool  isMuonIdICHEP = cfg.get<bool>("IsMuonIdICHEP"); 

  // electron selection
  const float ptEleCut   = cfg.get<float>("PtEleCut");
  const float etaEleCut  = cfg.get<float>("EtaEleCut");
  const float isoEleCut  = cfg.get<float>("IsoEleCut");
  const float dxyEleCut  = cfg.get<float>("dxyEleCut");
  const float dzEleCut   = cfg.get<float>("dzEleCut");
  
  // topological cuts (Z->mumu+Jet)
  const float ptLeadingMuCut_ZJet      = cfg.get<float>("PtLeadingMuCut_ZJet");
  const float ptTrailingMuCut_ZJet     = cfg.get<float>("PtTrailingMuCut_ZJet");
  const float ZMassLowerCut_ZJet       = cfg.get<float>("ZMassLowerCut_ZJet");
  const float ZMassUpperCut_ZJet       = cfg.get<float>("ZMassUpperCut_ZJet");
  const float ptJetZRatioLowerCut_ZJet = cfg.get<float>("PtJetZRatioLowerCut_ZJet");
  const float ptJetZRatioUpperCut_ZJet = cfg.get<float>("PtJetZRatioUpperCut_ZJet");
  const float deltaPhiZJetCut_ZJet     = cfg.get<float>("DeltaPhiZJetCut_ZJet");

  // topological cuts (W->muv+Jet)
  const float ptMuCut_WJet              = cfg.get<float>("PtMuCut_WJet");
  const float mtCut_WJet                = cfg.get<float>("MtCut_WJet");
  const float ptJetWRatioLowerCut_WJet  = cfg.get<float>("PtJetWRatioLowerCut_WJet");
  const float ptJetWRatioUpperCut_WJet  = cfg.get<float>("PtJetWRatioUpperCut_WJet");
  const float deltaPhiWJetCut_WJet      = cfg.get<float>("DeltaPhiWJetCut_WJet");

  // topological cuts (dijet)
  const float ptJetCut_DiJet              = cfg.get<float>("PtJetCut_DiJet");
  const float etaJetCut_DiJet             = cfg.get<float>("EtaJetCut_DiJet");
  const float ptTauJetRatioLowerCut_DiJet = cfg.get<float>("PtTauJetRatioLowerCut_DiJet");
  const float ptTauJetRatioUpperCut_DiJet = cfg.get<float>("PtTauJetRatioUpperCut_DiJet");
  const float deltaPhiTauJetCut_DiJet     = cfg.get<float>("DeltaPhiTauJetCut_DiJet");

  // topological cuts (trigger study)
  const float ZMassCut_Trig = cfg.get<float>("ZMassCut_Trig");
  const float mtCut_Trig = cfg.get<float>("MtCut_Trig");

  // momentum scales
  const float tauMomScale = cfg.get<float>("TauMomScale");
  const float muonMomScale = cfg.get<float>("MuonMomScale");
  const float eleMomScale = cfg.get<float>("EleMomScale");
  const int unclusteredES = cfg.get<int>("UnclusteredES");
  const int jetES = cfg.get<int>("JetES");

  const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  const string MuonTrigFile  = cfg.get<string>("MuonTrigEff");

  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);
  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string eventHistoName("eventCount/EventCount");
  std::string eventHistoNameData("makeroottree/nEvents");
  std::string weightsHistoName("eventCount/EventWeights");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");

  file->cd("");

  // histograms  
  // number of input events
  // and sum of eventweights
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);

  TH1D * WMassH     = new TH1D("WMassH","",300,0,3000);
  TH1D * WPtH       = new TH1D("WPtH",  "",100,0,1000);
  TH1D * WDecayH    = new TH1D("WTauDecayH","",5,-1.5,3.5);
  TH1D * WTauDecayH = new TH1D("WDecayH","",11,-1.5,9.5);
  TH1D * nVerticesH = new TH1D("nVerticesH","",50,-0.5,49.5);

  // ntuple variables

  UInt_t run_;
  UInt_t event_;
  
  Float_t puWeight_;
  Float_t genWeight_;
  Float_t weight_;
  Float_t genHt_;

  UInt_t nVert_;

  Bool_t metFilters_;

  Float_t wMass_;
  Float_t wPt_;
  Float_t wEta_;
  Float_t wPhi_;
  Int_t   wDecay_;
  Int_t   wTauDecay_;

  Float_t lepWPt_;
  Float_t lepWEta_;
  Float_t lepWPhi_;

  Float_t nuWPt_;
  Float_t nuWEta_;
  Float_t nuWPhi_;

  Float_t tauMatchJetPt_;
  Float_t tauMatchJetEta_;
  Float_t tauMatchJetPhi_;

  Float_t met_;
  Float_t metphi_;
  Float_t mttau_;
  Float_t mtmuon_;

  Float_t muonPt_;
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

  Float_t recoilRatio_;
  Float_t recoilDPhi_;

  Int_t   tauDecay_;
  Int_t   tauGenDecay_;
  Int_t   tauGenMatchDecay_;
  UInt_t  tauNtrk05_;
  UInt_t  tauNtrk08_;
  UInt_t  tauNtrk1_;

  Bool_t  tauDM_;
  Bool_t  tauNewDM_;

  Bool_t  tauLooseIso_;
  Bool_t  tauMediumIso_;
  Bool_t  tauTightIso_;

  Bool_t  tauLooseMvaIso_;
  Bool_t  tauMediumMvaIso_;
  Bool_t  tauTightMvaIso_;
  Bool_t  tauVTightMvaIso_;

  Bool_t tauAntiMuonLoose3_;
  Bool_t tauAntiMuonTight3_;

  Bool_t tauAntiElectronVLooseMVA5_;
  Bool_t tauAntiElectronLooseMVA5_;

  Bool_t tauAntiElectronVLooseMVA6_;
  Bool_t tauAntiElectronLooseMVA6_;

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

  Float_t jettauPt_;
  Float_t jettauEta_;
  Float_t jettauPhi_;

  Bool_t  jettauDM_;

  Bool_t  jettauLooseIso_;
  Bool_t  jettauMediumIso_;
  Bool_t  jettauTightIso_;

  Bool_t  jettauLooseMvaIso_;
  Bool_t  jettauMediumMvaIso_;
  Bool_t  jettauTightMvaIso_;
  Bool_t  jettauVTightMvaIso_;

  Bool_t jettauAntiMuonLoose3_;
  Bool_t jettauAntiMuonTight3_;

  Int_t   jettauDecay_;
  Int_t   jettauGenDecay_;

  Bool_t jettauAntiElectronVLooseMVA6_;
  Bool_t jettauAntiElectronLooseMVA6_;

  UInt_t jetChargedMult_;
  UInt_t jetNeutralMult_;
  UInt_t jetChargedHadMult_;
  Float_t jetNeutralEMEnergyFraction_;
  Float_t jetNeutralHadEnergyFraction_;

  Float_t jet2Pt_;
  Float_t jet2Eta_;
  Float_t jet2Phi_;

  UInt_t jet2ChargedMult_;
  UInt_t jet2NeutralMult_;
  UInt_t jet2ChargedHadMult_;
  Float_t jet2NeutralEMEnergyFraction_;
  Float_t jet2NeutralHadEnergyFraction_;

  Float_t mueffweight;
  Float_t mutrigweight;

  UInt_t nSelTaus_;
  UInt_t nTaus20_;
  UInt_t nTaus30_;

  Float_t JetHt_            ; // sumJetPtCentral30 + sumJetPtForward30
  Float_t SoftJetHt_        ; // sumJetPtCentral20 + sumJetPtForward30
  Float_t Ht_               ; // sumJetPtCentral30 + sumJetPtForward30 + sumLeptonPt
  Float_t SoftHt_           ; // sumJetPtCentral20 + sumJetPtForward30 + sumLeptonPt
  Float_t HtNoRecoil_       ; // sumJetPtCentral30 + sumJetPtForward30 + sumLeptonPt - sumPtRecoil
  Float_t SoftHtNoRecoil_   ; // sumJetPtCentral20 + sumJetPtForward30 + sumLeptonPt - sumPtRecoil 

  Int_t selection_; 
  //  0 : Z->mumu+Jet, 
  //  1 : W->muv+Jet
  //  2 : W*->muv 
  //  3 : W*->tauv
  //  4 : dijet
  // 10 : W->mu+v 
  // 11 : single jet + MEt
  // 12 : dijet sample

  Bool_t pfJet40_;
  Bool_t pfJet60_;
  Bool_t pfJet80_;
  Bool_t pfJet140_;

  Bool_t pf2Jet40_;
  Bool_t pf2Jet60_;
  Bool_t pf2Jet80_;
  Bool_t pf2Jet140_;
  
  ////////////
  Float_t taujetPt_;
  ///////////

  UInt_t nJets20_;
  Float_t jet20Pt_[10];
  Float_t jet20Eta_[10];
  Float_t jet20Phi_[10];

  TTree * wntuple_ = new TTree("WNTuple","WNTuple"); // small ntuple

  wntuple_->Branch("WMass",&wMass_,   "WMass/F");
  wntuple_->Branch("WPt",  &wPt_,     "WPt/F");
  wntuple_->Branch("WEta", &wEta_,    "WEta/F");
  wntuple_->Branch("WPhi", &wPhi_,    "WPhi/F");
  wntuple_->Branch("WDecay", &wDecay_,"WDecay/I");
  wntuple_->Branch("WTauDecay",&wTauDecay_, "WTauDecay/I");

  wntuple_->Branch("lepWPt",&lepWPt_,"lepWPt/F");
  wntuple_->Branch("lepWEta",&lepWEta_,"lepWEta/F");
  wntuple_->Branch("lepWPhi",&lepWPhi_,"lepWPhi/F");

  wntuple_->Branch("nuWPt",&nuWPt_,"nuWPt/F");
  wntuple_->Branch("nuWEta",&nuWEta_,"nuWEta/F");
  wntuple_->Branch("nuWPhi",&nuWPhi_,"nuWPhi/F");

  wntuple_->Branch("met",&met_,"met/F");
  wntuple_->Branch("metphi",&metphi_,"metphi/F");

  wntuple_->Branch("nJets20",&nJets20_,"nJets20/i");
  wntuple_->Branch("jet20Pt",jet20Pt_,"jet20Pt[nJets20]/F");
  wntuple_->Branch("jet20Eta",jet20Eta_,"jet20Eta[nJets20]/F");
  wntuple_->Branch("jet20Phi",jet20Phi_,"jet20Phi[nJets20]/F");

  wntuple_->Branch("tauPt",  &tauPt_,  "tauPt/F");
  wntuple_->Branch("tauEta", &tauEta_, "tauEta/F");
  wntuple_->Branch("tauPhi", &tauPhi_, "tauPhi/F");

  TTree * ntuple_ = new TTree("NTuple","NTuple");

  ntuple_->Branch("event",&event_,"event/i"); 
  ntuple_->Branch("run",  &run_,  "run/i");

  ntuple_->Branch("puWeight",  &puWeight_,  "puWeight/F");
  ntuple_->Branch("genWeight", &genWeight_, "genWeight/F");
  ntuple_->Branch("weight",    &weight_,    "weight/F");
  ntuple_->Branch("mueffweight", &mueffweight, "mueffweight/F");
  ntuple_->Branch("mutrigweight",&mutrigweight,"mutrigweight/F");

  ntuple_->Branch("genHt",&genHt_,"genHt/F");

  ntuple_->Branch("NVert",&nVert_,"NVert/i");
  ntuple_->Branch("metFilters",&metFilters_,"metFilters/O");

  ntuple_->Branch("WMass",&wMass_,   "WMass/F");
  ntuple_->Branch("WPt",  &wPt_,     "WPt/F");
  ntuple_->Branch("WEta", &wEta_,    "WEta/F");
  ntuple_->Branch("WPhi", &wPhi_,    "WPhi/F");
  ntuple_->Branch("WDecay", &wDecay_,"WDecay/I");
  ntuple_->Branch("WTauDecay",&wTauDecay_, "WTauDecay/I");

  ntuple_->Branch("lepWPt",&lepWPt_,"lepWPt/F");
  ntuple_->Branch("lepWEta",&lepWEta_,"lepWEta/F");
  ntuple_->Branch("lepWPhi",&lepWPhi_,"lepWPhi/F");

  ntuple_->Branch("nuWPt",&nuWPt_,"nuWPt/F");
  ntuple_->Branch("nuWEta",&nuWEta_,"nuWEta/F");
  ntuple_->Branch("nuWPhi",&nuWPhi_,"nuWPhi/F");

  ntuple_->Branch("met",    &met_,   "met/F");
  ntuple_->Branch("metphi", &metphi_,"metphi/F");
  ntuple_->Branch("mttau",  &mttau_, "mttau/F");
  ntuple_->Branch("mtmuon", &mtmuon_,"mtmuon/F");

  ntuple_->Branch("muonPt",  &muonPt_,  "muonPt/F");
  ntuple_->Branch("muonEta", &muonEta_, "muonEta/F");
  ntuple_->Branch("muonPhi", &muonPhi_, "muonPhi/F");
  ntuple_->Branch("muonQ",   &muonQ_,   "muonQ/I");

  ntuple_->Branch("muon2Pt",  &muon2Pt_,  "muon2Pt/F");
  ntuple_->Branch("muon2Eta", &muon2Eta_, "muon2Eta/F");
  ntuple_->Branch("muon2Phi", &muon2Phi_, "muon2Phi/F");
  ntuple_->Branch("muon2Q",   &muon2Q_,   "muon2Q/I");

  ntuple_->Branch("recoilRatio",&recoilRatio_,"recoilRatio/F");
  ntuple_->Branch("recoilDPhi",&recoilDPhi_,"recoilDPhi/F");

  ntuple_->Branch("recoilM",&recoilM_,"recoilM/F");
  ntuple_->Branch("recoilPt",&recoilPt_,"recoilPt/F");
  ntuple_->Branch("recoilEta",&recoilEta_,"recoilEta/F");
  ntuple_->Branch("recoilPhi",&recoilPhi_,"recoilPhi/F");

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

  ntuple_->Branch("jettauPt",&jettauPt_,"jettauPt/F");
  ntuple_->Branch("jettauEta",&jettauEta_,"jettauEta/F");
  ntuple_->Branch("jettauPhi",&jettauPhi_,"jettauPhi/F");

  ntuple_->Branch("jettauDM",&jettauDM_,"jettauDM/O");

  ntuple_->Branch("jettauLooseIso", &jettauLooseIso_, "jettauLooseIso/O");
  ntuple_->Branch("jettauMediumIso",&jettauMediumIso_,"jettauMediumIso/O");
  ntuple_->Branch("jettauTightIso", &jettauTightIso_, "jettauTightIso/O");

  ntuple_->Branch("jettauLooseMvaIso", &jettauLooseMvaIso_, "jettauLooseMvaIso/O");
  ntuple_->Branch("jettauMediumMvaIso",&jettauMediumMvaIso_,"jettauMediumMvaIso/O");
  ntuple_->Branch("jettauTightMvaIso", &jettauTightMvaIso_, "jettauTightMvaIso/O");
  ntuple_->Branch("jettauVTightMvaIso", &jettauVTightMvaIso_, "jettauVTightMvaIso/O");

  ntuple_->Branch("jettauAntiMuonLoose3",&jettauAntiMuonLoose3_,"jettauAntiMuonLoose3/O");
  ntuple_->Branch("jettauAntiMuonTight3",&jettauAntiMuonTight3_,"jettauAntiMuonTight3/O");

  ntuple_->Branch("jettauAntiElectronVLooseMVA6",&jettauAntiElectronVLooseMVA6_,"jettauAntiElectronVLooseMVA6/O");
  ntuple_->Branch("jettauAntiElectronLooseMVA6", &jettauAntiElectronLooseMVA6_, "jettauAntiElectronLooseMVA6/O");

  ntuple_->Branch("jettauDecay",&jettauDecay_,"jettauDecay/F");
  ntuple_->Branch("jettauGenDecay",&jettauGenDecay_,"jettauGenDecay/F");

  ntuple_->Branch("jetChargedMult",   &jetChargedMult_,   "jetChargedMult/i");
  ntuple_->Branch("jetNeutralMult",   &jetNeutralMult_,   "jetNeutralMult/i");
  ntuple_->Branch("jetChargedHadMult",&jetChargedHadMult_,"jetChargedHadMult/i");

  ntuple_->Branch("jetNeutralEMEnergyFraction", &jetNeutralEMEnergyFraction_, "jetNeutralEMEnergyFraction/F");
  ntuple_->Branch("jetNeutralHadEnergyFraction",&jetNeutralHadEnergyFraction_,"jetNeutralHadEnergyFraction/F");

  ntuple_->Branch("jet2Pt", &jet2Pt_, "jet2Pt/F");
  ntuple_->Branch("jet2Eta",&jet2Eta_,"jet2Eta/F");
  ntuple_->Branch("jet2Phi",&jet2Phi_,"jet2Phi/F");

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

  ntuple_->Branch("pfJet40",&pfJet40_,"pfJet40/O");
  ntuple_->Branch("pfJet60",&pfJet60_,"pfJet60/O");
  ntuple_->Branch("pfJet80",&pfJet80_,"pfJet80/O");
  ntuple_->Branch("pfJet140",&pfJet140_,"pfJet140/O");

  ntuple_->Branch("pf2Jet40",&pf2Jet40_,"pf2Jet40/O");
  ntuple_->Branch("pf2Jet60",&pf2Jet60_,"pf2Jet60/O");
  ntuple_->Branch("pf2Jet80",&pf2Jet80_,"pf2Jet80/O");
  ntuple_->Branch("pf2Jet140",&pf2Jet140_,"p2fJet140/O");

  ntuple_->Branch("Selection",&selection_,"Selection/I");

  TH1D * dRtauCentralJetH = new TH1D("dRtauCentralJetH","",50,0.,5.0);
  TH1D * dRtauForwardJetH = new TH1D("dRtauForwardJetH","",50,0.,5.0);

  // project directory
  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Good run selection
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
  
  // Official PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUDataFile,"read");
  TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUMCFile, "read");
  TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
  TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get("pileup");
  PUofficial->set_h_data(PUOfficial_data);
  PUofficial->set_h_MC(PUOfficial_mc);

  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  ScaleFactor * SF_muonTrig = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));

  // MEt filters
  std::vector<TString> metFlags; metFlags.clear();
  metFlags.push_back("Flag_HBHENoiseFilter");
  metFlags.push_back("Flag_HBHENoiseIsoFilter");
  metFlags.push_back("Flag_globalTightHalo2016Filter");
  metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  metFlags.push_back("Flag_goodVertices");
  metFlags.push_back("Flag_eeBadScFilter");
  metFlags.push_back("Flag_muonBadTrackFilter");
  metFlags.push_back("Flag_chargedHadronTrackResolutionFilter");

  int nFiles = 0;
  int nEvents = 0;
  int ZJetEvents = 0;
  int WJetEvents = 0;
  int WProdEvents = 0;
  int WMuNuEvents = 0;
  int WTauNuEvents = 0;
  int DiJetEvents = 0;
  int SingleJetEvents = 0;
  int JetTauEvents = 0;
  int TrigEvents = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  for (int iF=0; iF<nTotalFiles; ++iF) { // loop over files

    std::string filen;
    fileList >> filen;

    // opening file
    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    if (file_->IsZombie()) {
      cout << "Problems opening file : quitting program" << endl;
      exit(-1);
    }

    // accessing tree
    TTree * tree_ = (TTree*)file_->Get(TString(ntupleName));
    if (tree_==NULL) { 
      cout << "NTuple " << ntupleName << " is not found in file : quitting program" << endl;
      exit(-1);
    }
    AC1B analysisTree(tree_);
    Long64_t numberOfEntries = analysisTree.GetEntries();
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (nEvents%10000==0) cout << "Processed " << nEvents << endl;
      
      // ***************************
      // initialize ntuple variables
      // ***************************
      run_ = analysisTree.event_run;
      event_ = analysisTree.event_nr;
      nVert_ = analysisTree.primvertex_count;

      if (event_<=0) {
	std::cout << "Event : " << event_ << std::endl;
	std::cout << "From NTuple = " << analysisTree.event_nr << std::endl;
      }

      weight_ = 1;
      genWeight_ = 1;
      puWeight_ = 1;

      metFilters_ = metFiltersPasses(analysisTree,metFlags);

      wMass_ = -1;
      wPt_ = -1;
      wEta_ = 0;
      wPhi_ = 0;
      wDecay_ = -1;
      wTauDecay_ = -1;

      lepWPt_ = -1;
      lepWEta_ = 0;
      lepWPhi_ = 0;

      nuWPt_ = -1;
      nuWEta_ = 0;
      nuWPhi_ = 0;

      met_ =  -1;
      metphi_ =  0;
      mttau_ = 0;
      mtmuon_ = 0;

      muonPt_ =  -1;
      muonEta_ = 0;
      muonPhi_ = 0;
      muonQ_ = 0;

      muon2Pt_ =  -1;
      muon2Eta_ = 0;
      muon2Phi_ = 0;
      muon2Q_ = 0;

      tauPt_ = 0;
      tauEta_ = 0;
      tauPhi_ = 0;
      tauMass_ = 0;
      tauQ_ = 0;


      ///////////////
      taujetPt_ = 0;
      ///////////////

      genHt_ = analysisTree.genparticles_lheHt;

      recoilRatio_ = -1;
      recoilDPhi_ = 0;

      recoilM_ = -1;
      recoilPt_ = -1;
      recoilEta_ = 0;
      recoilPhi_ = 0;

      tauDecay_ = -1;
      tauGenDecay_ = -1;
      tauGenMatchDecay_ = -1;

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
      
      tauLooseIso_ = false;
      tauMediumIso_ = false;
      tauTightIso_ = false;

      tauLooseMvaIso_ = false;
      tauMediumMvaIso_ = false;
      tauTightMvaIso_ = false;
      tauVTightMvaIso_ = false;

      tauAntiMuonLoose3_ = false;
      tauAntiMuonTight3_ = false;

      tauAntiElectronVLooseMVA5_ = false;
      tauAntiElectronLooseMVA5_ = false;
      
      tauAntiElectronVLooseMVA6_ = false;
      tauAntiElectronLooseMVA6_ = false;
      
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

      jettauPt_ = 0;
      jettauEta_ = 0;
      jettauPhi_ = 0;

      jettauDM_ = false;
      jettauLooseIso_ = false;
      jettauMediumIso_ = false;
      jettauTightIso_ = false;
      jettauLooseMvaIso_ = false;
      jettauMediumMvaIso_ = false;
      jettauTightMvaIso_ = false;
      jettauVTightMvaIso_ = false;
      jettauAntiMuonLoose3_ = false;
      jettauAntiMuonTight3_ = false;

      jettauDecay_ = -1;
      jettauGenDecay_ = -1;

      jetChargedMult_ = 0;
      jetNeutralMult_ = 0;
      jetChargedHadMult_ = 0;

      jetNeutralEMEnergyFraction_ = 0;
      jetNeutralHadEnergyFraction_ = 0;

      jet2Pt_ = 0;
      jet2Eta_ = 0;
      jet2Phi_ = 0;

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

      pfJet40_ = false;
      pfJet60_ = false;
      pfJet80_ = false;
      pfJet140_ = false;

      pf2Jet40_ = false;
      pf2Jet60_ = false;
      pf2Jet80_ = false;
      pf2Jet140_ = false;

      mueffweight = 1;
      mutrigweight = 1;

      // booleans 
      bool isZJet = false;
      bool isWJet = false;
      bool isWTauNu = false;
      bool isWMuNu  = false;
      bool isDiJet = false;



      if (debug) {
	std::cout << "Run = " << analysisTree.event_nr << "    Event = " << analysisTree.event_run << std::endl; 
	std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
	std::cout << "Number of taus          = " << analysisTree.tau_count << std::endl;
	std::cout << "Number of jets          = " << analysisTree.pfjet_count << std::endl;
	std::cout << "Number of muons         = " << analysisTree.muon_count << std::endl;
	std::cout << "Number of electrons     = " << analysisTree.electron_count << std::endl;
      }	

      // applying genweight 
      if (!isData) {
	if (analysisTree.genweight<0)
	  genWeight_ = -1;
	else
	  genWeight_ = 1;
	weight_ *= genWeight_;
      }
      histWeightsH->Fill(double(0.),double(genWeight_));
      
      // **********************************
      // *** Analysis of generator info ***
      // **********************************
      int indexW  = -1;
      int indexNu = -1; // nu from W
      int indexMu = -1; // muon from W
      int indexE  = -1; // elec from W
      int indexTau = -1; // tau from W
      int indexTauE = -1; // W->tau->e
      int indexTauMu = -1; // W->tau->mu
      vector<TLorentzVector> gentauLV; gentauLV.clear();
      vector<int> gentauDecay; gentauDecay.clear();
      vector<TLorentzVector> genmuonLV; genmuonLV.clear();
      vector<TLorentzVector> genelecLV; genelecLV.clear();
      vector<TLorentzVector> gentaumuonLV; gentaumuonLV.clear();
      vector<TLorentzVector> gentauelecLV; gentauelecLV.clear();
      TLorentzVector wmuonLV; wmuonLV.SetXYZT(0,0,0,0);
      TLorentzVector welecLV; welecLV.SetXYZT(0,0,0,0);
      TLorentzVector wtauLV;  wtauLV.SetXYZT(0,0,0,0);
      TLorentzVector wnuLV;   wnuLV.SetXYZT(0,0,0,0);
      TLorentzVector wallnuLV; wallnuLV.SetXYZT(0,0,0,0);
      if (!isData) {
	for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	  
	  float pxGen = analysisTree.genparticles_px[igen];
	  float pyGen = analysisTree.genparticles_py[igen];
	  float pzGen = analysisTree.genparticles_pz[igen];
	  float etaGen = PtoEta(pxGen,pyGen,pzGen);
	  float ptGen  = PtoPt(pxGen,pyGen);

	  TLorentzVector genPartLV; genPartLV.SetXYZT(analysisTree.genparticles_px[igen],
						      analysisTree.genparticles_py[igen],
						      analysisTree.genparticles_pz[igen],
						      analysisTree.genparticles_e[igen]);
	  
	  if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
	    indexW = igen;

	  if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==12 || 
	      TMath::Abs(analysisTree.genparticles_pdgid[igen])==14 ||
	      TMath::Abs(analysisTree.genparticles_pdgid[igen])==16) { 

	    if (analysisTree.genparticles_info[igen]==(1<<1)) {
	      indexNu = igen;
	      wnuLV = genPartLV;
	    }
	    if (analysisTree.genparticles_info[igen]==(1<<1) ||
		analysisTree.genparticles_info[igen]==((1<<1)|(1<<2))) {
              wallnuLV += genPartLV;
            }

	  }
	  
	  if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==13) {
	    if ( analysisTree.genparticles_info[igen]==(1<<1) ||
		 analysisTree.genparticles_info[igen]==(1<<0) ) // W/Z->mu
	      genmuonLV.push_back(genPartLV);
	    if ( analysisTree.genparticles_info[igen]==((1<<0)|(1<<2)) ||
		 analysisTree.genparticles_info[igen]==((1<<1)|(1<<2))) // W/Z -> tau -> mu
	      gentaumuonLV.push_back(genPartLV);
	    if ( analysisTree.genparticles_info[igen]==(1<<1) ) { // W -> muv  
	      indexMu = igen;
	      wmuonLV = genPartLV;
	    }
	    if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->mu
	      indexTauMu = igen;
	  }

	  if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==11) { // electron
	    if ( analysisTree.genparticles_info[igen]==(1<<1) || 
		 analysisTree.genparticles_info[igen]==(1<<2) ) // W/Z->e 
	      genelecLV.push_back(genPartLV);
	    if ( analysisTree.genparticles_info[igen]==((1<<0)|(1<<2)) || 
		 analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W/Z -> tau -> e 
	      gentauelecLV.push_back(genPartLV);
	    if ( analysisTree.genparticles_info[igen]==(1<<1) ) { // W->ev
              indexE = igen;
	      welecLV = genPartLV;
	    }
            if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->e
	      indexTauE = igen;
	  }
	}

	for (unsigned int igentau=0; igentau<analysisTree.gentau_count;++igentau) {
	  TLorentzVector GenVisTau; GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
						      analysisTree.gentau_visible_py[igentau],
						      analysisTree.gentau_visible_pz[igentau],
						      analysisTree.gentau_visible_e[igentau]);
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau] ) { // W/Z->tau  
	    gentauLV.push_back(GenVisTau);
	    gentauDecay.push_back(analysisTree.gentau_decayMode[igentau]);
	    indexTau = igentau;
	    wtauLV = GenVisTau;
	  }
	}
      }

      TLorentzVector lorentzVectorGenW; lorentzVectorGenW.SetXYZT(0,0,0,0);
      if (indexW>=0) 
	lorentzVectorGenW.SetXYZT(analysisTree.genparticles_px[indexW],
				  analysisTree.genparticles_py[indexW],
				  analysisTree.genparticles_pz[indexW],
				  analysisTree.genparticles_e[indexW]);

      if (indexW>=0) {
	wMass_ = lorentzVectorGenW.M();
	wPt_   = lorentzVectorGenW.Pt();
	wEta_  = lorentzVectorGenW.Eta();
	wPhi_  = lorentzVectorGenW.Phi();
	nuWPt_ = wnuLV.Pt();
	nuWEta_ = wnuLV.Eta();
	nuWPhi_ = wnuLV.Phi();
	wDecay_ = 0;
	wTauDecay_ = -1;
	if (indexMu>=0) {
	  wDecay_ = 2;
	  lepWPt_  = wmuonLV.Pt();
	  lepWEta_ = wmuonLV.Eta(); 
	  lepWPhi_ = wmuonLV.Phi();
	}
	else if (indexE>=0) {
	  wDecay_ = 1;
	  lepWPt_  = welecLV.Pt();
          lepWEta_ = welecLV.Eta();
          lepWPhi_ = welecLV.Phi();
	}
	else if (indexTau>=0) {
	  lepWPt_  = wtauLV.Pt();
          lepWEta_ = wtauLV.Eta();
          lepWPhi_ = wtauLV.Phi();
	  wDecay_ = 3;
	  wTauDecay_ = analysisTree.gentau_decayMode[indexTau];
	  if (wTauDecay_<0) wTauDecay_ = -1;
	}
	WMassH->Fill(wMass_,weight_);
	WPtH->Fill(wPt_,weight_);
	WDecayH->Fill(wDecay_,weight_);
	WTauDecayH->Fill(wTauDecay_,weight_);
      }

      if (analysisTree.primvertex_count==0) continue; // at least one good primary vertex

      if (!isData) { 
	puWeight_ =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      	weight_ *= puWeight_; 
	if (debug)
	  cout << "nPU = " << analysisTree.numtruepileupinteractions << " --> puweight = " << puWeight_ << endl;
      }

      // ***********************************
      // applying good run selection on data
      // ***********************************
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


      // *********************************
      // ***** accessing trigger info ****
      // *********************************

      unsigned int nSingleMuonHLTFilter = 0;
      bool isSingleMuonHLTFilter = false;

      unsigned int nPFJet40HLTFilter = 0;
      bool isPFJet40HLTFilter = false;

      unsigned int nPFJet60HLTFilter = 0;
      bool isPFJet60HLTFilter = false;
      
      unsigned int nPFJet80HLTFilter = 0;
      bool isPFJet80HLTFilter = false;

      unsigned int nPFJet140HLTFilter = 0;
      bool isPFJet140HLTFilter = false;
      
      for (unsigned int i=0; i<analysisTree.run_hltfilters->size(); ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==SingleMuonHLTFilterName) {
	  nSingleMuonHLTFilter = i;
	  isSingleMuonHLTFilter = true;
	}
	if (HLTFilter==PFJet40HLTFilterName) {
	  nPFJet40HLTFilter = i;
	  isPFJet40HLTFilter = true;
	}
	if (HLTFilter==PFJet60HLTFilterName) {
	  nPFJet60HLTFilter = i;
	  isPFJet60HLTFilter = true;
	}
	if (HLTFilter==PFJet80HLTFilterName) {
	  nPFJet80HLTFilter = i;
	  isPFJet80HLTFilter = true;
	}
	if (HLTFilter==PFJet140HLTFilterName) {
	  nPFJet140HLTFilter = i;
	  isPFJet140HLTFilter = true;
	}
      }
      if (!isSingleMuonHLTFilter) {
	std::cout << "HLT filter " << SingleMuonHLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet40HLTFilter) {
	std::cout << "HLT filter " << PFJet40HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet60HLTFilter) {
	std::cout << "HLT filter " << PFJet60HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet80HLTFilter) {
	std::cout << "HLT filter " << PFJet80HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet140HLTFilter) {
	std::cout << "HLT filter " << PFJet140HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      // ************************************
      // **** end accessing trigger info ****
      // ************************************


      // ***************************************************
      // accessing PF MET and changing momentum scale of met
      // ***************************************************
      float pfmet_ex = analysisTree.pfmet_ex;
      float pfmet_ey = analysisTree.pfmet_ey;
      if (jetES<0) {
	pfmet_ex = analysisTree.pfmet_ex_JetEnDown;
	pfmet_ey = analysisTree.pfmet_ey_JetEnDown;
      }
      else if (jetES>0) {
	pfmet_ex = analysisTree.pfmet_ex_JetEnUp;
        pfmet_ey = analysisTree.pfmet_ey_JetEnUp;
      }
      else if (unclusteredES<0) {
	pfmet_ex = analysisTree.pfmet_ex_UnclusteredEnDown;
	pfmet_ey = analysisTree.pfmet_ey_UnclusteredEnDown;
      }
      else if (unclusteredES>0) {
	pfmet_ex = analysisTree.pfmet_ex_UnclusteredEnUp;
        pfmet_ey = analysisTree.pfmet_ey_UnclusteredEnUp;
      }
      else {
	pfmet_ex = analysisTree.pfmet_ex;
	pfmet_ey = analysisTree.pfmet_ey;
      }
      //      cout << endl;
      //      cout << "metx            = " << pfmet_ex << endl;
      //      cout << "metx(JetEnUp)   = " << analysisTree.pfmet_ex_JetEnUp << endl; 
      //      cout << "metx(JetEnDown) = " << analysisTree.pfmet_ex_JetEnDown << endl; 
      //      cout << "metx(UncEnUp)   = " << analysisTree.pfmet_ex_UnclusteredEnUp << endl; 
      //      cout << "metx(UncEnDown) = " << analysisTree.pfmet_ex_UnclusteredEnDown << endl; 
      //      cout << endl;
      //      cout << "mety            = " << pfmet_ey << endl;
      //      cout << "mety(JetEnUp)   = " << analysisTree.pfmet_ey_JetEnUp << endl; 
      //      cout << "mety(JetEnDown) = " << analysisTree.pfmet_ey_JetEnDown << endl; 
      //      cout << "mety(UncEnUp)   = " << analysisTree.pfmet_ey_UnclusteredEnUp << endl; 
      //      cout << "mety(UncEnDown) = " << analysisTree.pfmet_ey_UnclusteredEnDown << endl; 
      //      cout << endl;

      met_ = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      if (met_<1e-4) met_ = 1e-4;
      metphi_ = TMath::ATan2(pfmet_ey,pfmet_ex);
      TLorentzVector lorentzVectorMet; lorentzVectorMet.SetXYZM(pfmet_ex,pfmet_ey,0,0);

      // *************************
      // **** accessing muons ****
      // *************************
      TLorentzVector lorentzVectorAllMuons; lorentzVectorAllMuons.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorAllSelMuons; lorentzVectorAllSelMuons.SetXYZT(0,0,0,0);
      std::vector<unsigned int> muonIndexes; muonIndexes.clear();
      std::vector<unsigned int> selMuonIndexes; selMuonIndexes.clear();
      std::vector<TLorentzVector>  lorentzVectorMuons; lorentzVectorMuons.clear();
      std::vector<TLorentzVector>  lorentzVectorSelMuons; lorentzVectorSelMuons.clear();
      int indexTriggerMu = -1;
      float ptTriggerMu  = -1;
      float etaTriggerMu = -1; 
      float muonHt = 0;
      for (unsigned int imuon=0; imuon<analysisTree.muon_count; ++imuon) {
	analysisTree.muon_px[imuon] *= muonMomScale;
	analysisTree.muon_py[imuon] *= muonMomScale;
	analysisTree.muon_pz[imuon] *= muonMomScale;
	analysisTree.muon_pt[imuon] *= muonMomScale;
	if (analysisTree.muon_pt[imuon]<ptMuCut) continue;
	if (fabs(analysisTree.muon_eta[imuon])>etaMuCut) continue;
	bool passedId = analysisTree.muon_isMedium[imuon];
	if (isMuonIdICHEP) {
	  bool goodGlob =
	    analysisTree.muon_isGlobal[imuon] &&
	    analysisTree.muon_normChi2[imuon] < 3 &&
	    analysisTree.muon_combQ_chi2LocalPosition[imuon] < 12 &&
	    analysisTree.muon_combQ_trkKink[imuon] < 20;
	  passedId =
	    analysisTree.muon_isLoose[imuon] &&
	    analysisTree.muon_validFraction[imuon] >0.49 &&
	    analysisTree.muon_segmentComp[imuon] > (goodGlob ? 0.303 : 0.451);
	}
	if (!passedId) continue;
	bool passedIpCuts = 
	  fabs(analysisTree.muon_dxy[imuon]) < dxyMuCut &&
	  fabs(analysisTree.muon_dz[imuon]) < dzMuCut;
	if (!passedIpCuts) continue;
	float absIso = 0; 
	if (isDRIso03) {
	  absIso = analysisTree.muon_r03_sumChargedHadronPt[imuon];
	  float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[imuon] + 
	    analysisTree.muon_r03_sumPhotonEt[imuon] - 0.5*analysisTree.muon_r03_sumPUPt[imuon];
	  neutralIso = TMath::Max(float(0),neutralIso);
	  absIso += neutralIso;
	}
	else {
	  absIso = analysisTree.muon_chargedHadIso[imuon];
          float neutralIso = analysisTree.muon_neutralHadIso[imuon] +
            analysisTree.muon_photonIso[imuon] -
            0.5*analysisTree.muon_puIso[imuon];
          neutralIso = TMath::Max(float(0),neutralIso);
          absIso += neutralIso;
	} 
	float relIso = absIso/analysisTree.muon_pt[imuon];
	bool passedIso = relIso < isoMuCut;
	if (!passedIso) continue;

	// all muons -->
	muonIndexes.push_back(imuon);
	TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[imuon],
					      analysisTree.muon_py[imuon],
					      analysisTree.muon_pz[imuon],
					      muonMass);
	lorentzVectorAllMuons += muonLV;
	lorentzVectorMuons.push_back(muonLV);
	muonHt += analysisTree.muon_pt[imuon];

	// selected muons -->
	bool passedSelIso = relIso < isoSelMuCut;
	if (analysisTree.muon_pt[imuon]>ptSelMuCut && 
	    passedSelIso) {
	  selMuonIndexes.push_back(imuon);
	  TLorentzVector muonSelLV; muonSelLV.SetXYZM(analysisTree.muon_px[imuon],
						      analysisTree.muon_py[imuon],
						      analysisTree.muon_pz[imuon],
						      muonMass);
	  lorentzVectorAllSelMuons += muonSelLV;
	  lorentzVectorSelMuons.push_back(muonSelLV);
	}

	// triggering muon -->
	if (analysisTree.muon_pt[imuon]>ptTrigMuCut && 
	    passedSelIso) {
	  bool trigMatch = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[imuon],analysisTree.muon_phi[imuon],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>0.5) continue;
	    if (analysisTree.trigobject_filters[iT][nSingleMuonHLTFilter]) trigMatch = true;

	  }
	  if (trigMatch&&analysisTree.muon_pt[imuon]>ptTriggerMu) {
	    ptTriggerMu = analysisTree.muon_pt[imuon];
	    etaTriggerMu = analysisTree.muon_eta[imuon];
	    indexTriggerMu = int(imuon);
	  }
	}
      }

      nMuon_ = muonIndexes.size();
      nSelMuon_ = selMuonIndexes.size();

      float ptSecondMu  = -1;
      float etaSecondMu = -1;
      int indexSecondMu = -1;
      TLorentzVector lorentzVectorTriggerMu; lorentzVectorTriggerMu.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorSecondMu;  lorentzVectorSecondMu.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorZ; lorentzVectorZ.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorW; lorentzVectorW.SetXYZT(0,0,0,0);
      if (indexTriggerMu>=0) {
	lorentzVectorTriggerMu.SetXYZM(analysisTree.muon_px[indexTriggerMu],
				       analysisTree.muon_py[indexTriggerMu],
				       analysisTree.muon_pz[indexTriggerMu],
				       muonMass);
	muonPt_  = lorentzVectorTriggerMu.Pt();
	muonEta_ = lorentzVectorTriggerMu.Eta();
	muonPhi_ = lorentzVectorTriggerMu.Phi();
	muonQ_   = int(analysisTree.muon_charge[indexTriggerMu]);
	mtmuon_  = mT(lorentzVectorTriggerMu,lorentzVectorMet);
	lorentzVectorW = lorentzVectorTriggerMu + lorentzVectorMet;
	for (unsigned int iMu = 0; iMu < selMuonIndexes.size(); ++iMu) {
	  int indexMu = int(selMuonIndexes.at(iMu));
	  if (indexMu==indexTriggerMu) continue;
	  float netcharge = analysisTree.muon_charge[indexTriggerMu]*analysisTree.muon_charge[indexMu];
	  if (netcharge>0) continue;
	  if (analysisTree.muon_pt[indexMu]>ptSecondMu) {
	    ptSecondMu = analysisTree.muon_pt[indexMu];
	    etaSecondMu = analysisTree.muon_eta[indexMu];
	    indexSecondMu = int(indexMu);
	  }
	}
	if (indexSecondMu>=0) {
	  lorentzVectorSecondMu.SetXYZM(analysisTree.muon_px[indexSecondMu],
					analysisTree.muon_py[indexSecondMu],
					analysisTree.muon_pz[indexSecondMu],
					muonMass);
	  lorentzVectorZ = lorentzVectorTriggerMu + lorentzVectorSecondMu;
	  muon2Pt_  = lorentzVectorSecondMu.Pt();
	  muon2Eta_ = lorentzVectorSecondMu.Eta();
	  muon2Phi_ = lorentzVectorSecondMu.Phi();
	  muon2Q_   = int(analysisTree.muon_charge[indexSecondMu]);
	  float zmass = lorentzVectorZ.M();
	  if (zmass>60&&zmass<120) nVerticesH->Fill(double(analysisTree.primvertex_count),weight_); // vertex distribution in Z->mumu
	}
      }
      float ptLeadingMu = ptTriggerMu;
      float ptTrailingMu = ptSecondMu;
      float etaLeadingMu = etaTriggerMu;
      float etaTrailingMu = etaSecondMu;

      if (ptTrailingMu>ptLeadingMu) {
	ptLeadingMu = ptSecondMu;
	ptTrailingMu = ptTriggerMu;
	etaLeadingMu = etaSecondMu;
	etaTrailingMu = etaTriggerMu;
      }
      // *****************************
      // **** end accessing muons ****
      // *****************************

      if (debug)
	std::cout << "end accessing muons " << std::endl;

      // *****************************
      // **** accessing electrons ****
      // *****************************
      TLorentzVector lorentzVectorAllElectrons; lorentzVectorAllElectrons.SetXYZT(0,0,0,0);
      std::vector<unsigned int> eleIndexes; eleIndexes.clear();
      std::vector<TLorentzVector>  lorentzVectorElectrons; lorentzVectorElectrons.clear();
      float elecHt = 0;
      for (unsigned int ielec=0; ielec<analysisTree.electron_count; ++ielec) {
	analysisTree.electron_px[ielec] *= eleMomScale;
        analysisTree.electron_py[ielec] *= eleMomScale;
        analysisTree.electron_pz[ielec] *= eleMomScale;
        analysisTree.electron_pt[ielec] *= eleMomScale;
	bool passedId = 
	  analysisTree.electron_cutId_veto_Spring15[ielec] &&
          analysisTree.electron_pass_conversion[ielec] &&
          analysisTree.electron_nmissinginnerhits[ielec] <= 1;
	bool passedIpCuts = 
	  fabs(analysisTree.electron_dxy[ielec]) < dxyEleCut &&
	  fabs(analysisTree.electron_dz[ielec]) < dzEleCut;
	float absIso = analysisTree.electron_r03_sumChargedHadronPt[ielec];
	float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[ielec] + 
	  analysisTree.electron_r03_sumPhotonEt[ielec] - 0.5*analysisTree.electron_r03_sumPUPt[ielec];
	neutralIso = TMath::Max(float(0),neutralIso);
	absIso += neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ielec];
	bool passedIso = relIso < isoEleCut;
	if (analysisTree.electron_pt[ielec]>ptEleCut && 
	    fabs(analysisTree.electron_eta[ielec])<etaEleCut &&
	    passedId && passedIpCuts && passedIso) {
	  TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[ielec],
							analysisTree.electron_py[ielec],
							analysisTree.electron_pz[ielec],
							electronMass);
	  lorentzVectorAllElectrons += electronLV; 
	  eleIndexes.push_back(ielec);
	  elecHt = analysisTree.electron_pt[ielec];
	}
      }
      nElec_ = eleIndexes.size();
      // *********************************
      // **** end accessing electrons ****
      // *********************************     

      if (debug)
	std::cout << "end accessing electrons " << std::endl;

      // ************************
      // **** accessing jets ****
      // ************************
      TLorentzVector lorentzVectorAllJetsForMht; lorentzVectorAllJetsForMht.SetXYZT(0,0,0,0);
      std::vector<unsigned int> centralJets20Indexes; centralJets20Indexes.clear();
      std::vector<unsigned int> forwardJets20Indexes; forwardJets20Indexes.clear();
      std::vector<unsigned int> centralJets30Indexes; centralJets30Indexes.clear();
      std::vector<unsigned int> forwardJets30Indexes; forwardJets30Indexes.clear();
      float htCentral20 = 0;
      float htCentral30 = 0;
      float htForward20 = 0;
      float htForward30 = 0;
      std::vector<unsigned int> triggerJetsIndexes; triggerJetsIndexes.clear();
      std::vector<bool> jets40trigger; jets40trigger.clear();
      std::vector<bool> jets60trigger; jets60trigger.clear();
      std::vector<bool> jets80trigger; jets80trigger.clear();
      std::vector<bool> jets140trigger; jets140trigger.clear();
      //      if (analysisTree.pfjet_count>100)
	//	std::cout << "pfjet_count = " << analysisTree.pfjet_count << endl;
      unsigned int indexLeadingJet = 0; 
      float ptLeadingJet = -1;
      for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {

	float scaleJ = 1;
	//	cout << "Jet " << ijet << " Pt = " << analysisTree.pfjet_pt[ijet] << "  Unc = " << analysisTree.pfjet_jecUncertainty[ijet] << endl;

	if (jetES<0)
	  scaleJ = 1.0 - analysisTree.pfjet_jecUncertainty[ijet];
	else if (jetES>0)
	  scaleJ = 1.0 + analysisTree.pfjet_jecUncertainty[ijet];
	else
	  scaleJ = 1.0;

	//	std::cout << ijet << "  :  " << scaleJ << std::endl;

	analysisTree.pfjet_px[ijet] *= scaleJ;
	analysisTree.pfjet_py[ijet] *= scaleJ;
	analysisTree.pfjet_pz[ijet] *= scaleJ;
	analysisTree.pfjet_pt[ijet] *= scaleJ;
	analysisTree.pfjet_e[ijet]  *= scaleJ;

	float absJetEta = fabs(analysisTree.pfjet_eta[ijet]);

	if (absJetEta>5.2) continue;
	if (analysisTree.pfjet_pt[ijet]<20.0) continue; 
	
	// jetId
	bool isPFLooseJetId = looseJetiD(analysisTree,int(ijet));
	bool isPFTightJetId = tightJetiD(analysisTree,int(ijet));
	//	bool isPULooseJetId = looseJetPUiD(analysisTree,int(ijet));

	// jet four-vector
	TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[ijet],
					    analysisTree.pfjet_py[ijet],
					    analysisTree.pfjet_pz[ijet],
					    analysisTree.pfjet_e[ijet]);
	// counting jets for Mht
	if (isPFTightJetId) {
	  lorentzVectorAllJetsForMht += jetLV;
	}

	// accept only jets with looseId and loosePUId
	if (!isPFLooseJetId) continue;
	//	if (!isPULooseJetId) continue;

	// checking overlap with muons
	bool overlapWithMuon = false;
	for (unsigned int iMu=0; iMu<muonIndexes.size(); ++iMu) {
	  unsigned int indexMu = muonIndexes.at(iMu);
	  float dRJetLep = deltaR(analysisTree.muon_eta[indexMu],analysisTree.muon_phi[indexMu],
				  analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet]);
	  if (dRJetLep<0.4) {
	    overlapWithMuon = true;
	    break;
	  }
	}
	if (overlapWithMuon) continue;
	  
	// checking overlap with electrons
	bool overlapWithEle = false;
	for (unsigned int iEle=0; iEle<eleIndexes.size(); ++iEle) {
	  unsigned int indexEle = eleIndexes.at(iEle);
	  float dRJetLep = deltaR(analysisTree.electron_eta[indexEle],analysisTree.electron_phi[indexEle],
				  analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet]);
	  if (dRJetLep<0.4) {
	    overlapWithEle = true;
	    break;
	  }
	}
	if (overlapWithEle) continue;

	// pt > 20 GeV
	if (analysisTree.pfjet_pt[ijet]>20.0) {
	  if (absJetEta<2.4) {
	    centralJets20Indexes.push_back(ijet);
	    htCentral20 += analysisTree.pfjet_pt[ijet];
	    if (analysisTree.pfjet_pt[ijet]>ptLeadingJet)
	      indexLeadingJet = ijet;
	  }
	  else if (absJetEta<4.7) {
	    forwardJets20Indexes.push_back(ijet);
	    htForward20 += analysisTree.pfjet_pt[ijet]; 
	  }
	  if (nJets20_<10) {
	    jet20Pt_[nJets20_]  = analysisTree.pfjet_pt[ijet];
	    jet20Eta_[nJets20_] = analysisTree.pfjet_eta[ijet];
	    jet20Phi_[nJets20_] = analysisTree.pfjet_phi[ijet];
	    nJets20_++;
	  }
	}

	// pt > 30 GeV
	if (analysisTree.pfjet_pt[ijet]>30) {
	  if (absJetEta<2.4) {
	    centralJets30Indexes.push_back(ijet);
	    htCentral30 += analysisTree.pfjet_pt[ijet];
	  }
	  else if (absJetEta<4.7) {
	    forwardJets30Indexes.push_back(ijet);
	    htForward30 += analysisTree.pfjet_pt[ijet];
	  }
	}

	// triggering jets
	if (analysisTree.pfjet_pt[ijet]<ptJetCut_DiJet) continue;
	if (absJetEta>etaJetCut_DiJet) continue;
	bool trigMatch40 = false;
	bool trigMatch60 = false;
	bool trigMatch80 = false;
	bool trigMatch140 = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;
	  if ((analysisTree.trigobject_filters[iT][nPFJet40HLTFilter]&&isPFJet40HLTFilter))   trigMatch40 = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   trigMatch60 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   trigMatch80 = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) trigMatch140 = true;
	  
	}
	if (!isData) {
	  trigMatch40 = true;
	  trigMatch60 = true;
	  trigMatch80 = true;
	  trigMatch140 = true;
	}
	triggerJetsIndexes.push_back(ijet);
	jets40trigger.push_back(trigMatch40);
	jets60trigger.push_back(trigMatch60);
	jets80trigger.push_back(trigMatch80);
	jets140trigger.push_back(trigMatch140);
      }
      nJetsCentral20_ = centralJets20Indexes.size();
      nJetsCentral30_ = centralJets30Indexes.size();
      nJetsForward20_ = forwardJets20Indexes.size();
      nJetsForward30_ = forwardJets30Indexes.size();
      JetHt_     = htCentral30 + htForward30;
      SoftJetHt_ = htCentral20 + htForward30;
      Ht_        = JetHt_     + muonHt + elecHt;
      SoftHt_    = SoftJetHt_ + muonHt + elecHt;
      TLorentzVector lorentzVectorJet; lorentzVectorJet.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorJet2; lorentzVectorJet2.SetXYZT(0,0,0,0);
      if (nJetsCentral20_>0) {
	unsigned int indexJet0 = indexLeadingJet;
	jetPt_ = analysisTree.pfjet_pt[indexJet0];
	jetEta_ = analysisTree.pfjet_eta[indexJet0];
	jetPhi_ = analysisTree.pfjet_phi[indexJet0];
	jetChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet0];
	jetNeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet0];
	jetChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet0];
	jetNeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];
	jetNeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];  
	pfJet40_ = false;
	pfJet60_ = false;
	pfJet80_ = false;
	pfJet140_ = false;
	lorentzVectorJet.SetXYZT(analysisTree.pfjet_px[indexJet0],
				 analysisTree.pfjet_py[indexJet0],
				 analysisTree.pfjet_pz[indexJet0],
				 analysisTree.pfjet_e[indexJet0]);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet0],analysisTree.pfjet_phi[indexJet0],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;
	  if ((analysisTree.trigobject_filters[iT][nPFJet40HLTFilter]&&isPFJet40HLTFilter))   pfJet40_  = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pfJet60_  = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pfJet80_  = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pfJet140_ = true;
	}
	float dRmin = 0.6;
	unsigned int taujetIndex = 0;
	for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) {
	  TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[itau],
					      analysisTree.tau_py[itau],
					      analysisTree.tau_pz[itau],
					      analysisTree.tau_mass[itau]);

	  float dZ = fabs(analysisTree.tau_vertexz[itau]-analysisTree.primvertex_z);
	  if (dZ>1e-4) continue;

	  float deltaRjettau = deltaR(lorentzVectorJet.Eta(),lorentzVectorJet.Phi(),
				      tauLV.Eta(),tauLV.Phi());

	  if (deltaRjettau<dRmin) {
	    dRmin = deltaRjettau;
	    taujetIndex = itau;
	  }
	}
	if (dRmin<0.5) {
	  jettauPt_ = analysisTree.tau_pt[taujetIndex];
	  jettauEta_ = analysisTree.tau_eta[taujetIndex];
	  jettauPhi_ = analysisTree.tau_phi[taujetIndex];
	  
	  jettauDM_ = analysisTree.tau_decayModeFinding[taujetIndex] > 0.5;
	  jettauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  jettauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  jettauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  
	  jettauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  
	  jettauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[taujetIndex] > 0.5;
	  jettauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[taujetIndex] > 0.5;
	  
	  jettauAntiElectronVLooseMVA6_ = analysisTree.tau_againstElectronVLooseMVA6[taujetIndex] > 0.5;
	  jettauAntiElectronLooseMVA6_ = analysisTree.tau_againstElectronLooseMVA6[taujetIndex] > 0.5;
	  
	  jettauDecay_ = analysisTree.tau_decayMode[taujetIndex];
	  jettauGenDecay_ = analysisTree.tau_genDecayMode[taujetIndex];
	  
	}

      }
      if (nJetsCentral20_>1) {
	unsigned int indexJet1 = centralJets20Indexes.at(1);
	if (indexJet1 == indexLeadingJet) indexJet1 = centralJets20Indexes.at(0);
	jet2Pt_ = analysisTree.pfjet_pt[indexJet1];
	jet2Eta_ = analysisTree.pfjet_eta[indexJet1];
	jet2Phi_ = analysisTree.pfjet_phi[indexJet1];
	jet2ChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet1];
	jet2NeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet1];
	jet2ChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet1];
	jet2NeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];
	jet2NeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];  
	pf2Jet40_ = false;
	pf2Jet60_ = false;
	pf2Jet80_ = false;
	pf2Jet140_ = false;
	lorentzVectorJet2.SetXYZT(analysisTree.pfjet_px[indexJet1],
				  analysisTree.pfjet_py[indexJet1],
				  analysisTree.pfjet_pz[indexJet1],
				  analysisTree.pfjet_e[indexJet1]);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet1],analysisTree.pfjet_phi[indexJet1],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;
	  if ((analysisTree.trigobject_filters[iT][nPFJet40HLTFilter]&&isPFJet40HLTFilter))   pf2Jet40_  = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pf2Jet60_  = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pf2Jet80_  = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pf2Jet140_ = true;
	}
      }
      

      // ****************************
      // **** end accessing jets ****
      // ****************************

      if (debug)
	std::cout << "end accessing jets" << std::endl;

      // ***************************
      // ******* WJet selection ****
      // ***************************
      if (lorentzVectorW.Pt()>1e-4) {

	recoilRatio_ = jetPt_ / lorentzVectorW.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorJet);

	isWJet = ptTriggerMu>ptMuCut_WJet && jetPt_>20; 
	isWJet = isWJet && mtmuon_ > mtCut_WJet;
	isWJet = isWJet && recoilRatio_>ptJetWRatioLowerCut_WJet && recoilRatio_<ptJetWRatioUpperCut_WJet;
	isWJet = isWJet && recoilDPhi_>deltaPhiWJetCut_WJet;
	if (isWJet) {
	  mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
          mutrigweight = SF_muonTrig->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
	  HtNoRecoil_     = Ht_     - ptTriggerMu;
	  SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu;
	  recoilM_   = lorentzVectorW.M();
          recoilPt_  = lorentzVectorW.Pt();
          recoilEta_ = lorentzVectorW.Eta();
          recoilPhi_ = lorentzVectorW.Phi();
	  selection_ = 1;
	  ntuple_->Fill();
	  WJetEvents++;
	}
      }
      

    } // end of file processing (loop over events in one file)
    nFiles++;
    delete tree_;
    file_->Close();
    delete file_;
  }
  int allEvents   = int(inputEventsH->GetSumOfWeights());
  double sumOfWeights = histWeightsH->GetSumOfWeights();
  std::cout << "Total number of input events      = " << allEvents << std::endl;
  std::cout << "Total weight sum                  = " << sumOfWeights << std::endl;
  std::cout << "Total number of events in Tree    = " << nEvents << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}



