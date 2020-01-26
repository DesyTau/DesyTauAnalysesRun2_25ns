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

  // electron selection
  const float ptEleCut   = cfg.get<float>("PtEleCut");
  const float etaEleCut  = cfg.get<float>("EtaEleCut");
  const float isoEleCut  = cfg.get<float>("IsoEleCut");
  const float dxyEleCut  = cfg.get<float>("dxyEleCut");
  const float dzEleCut   = cfg.get<float>("dzEleCut");
  const float isoSelEleCut   = cfg.get<float>("IsoSelEleCut");

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
  const string puMCFile   = cfg.get<string>("PileUpMCFile");
  const string puMCHist   = cfg.get<string>("PileUpMCHist");

  const bool JetLeptonFake = cfg.get<bool>("JetLeptonFake");

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);
  TString PUMCHist(puMCHist);
  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
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
  Int_t nPartons_;
  Int_t nPartonsNLO_;
  Float_t genWPt_;

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

  UInt_t nMuon_;
  UInt_t nSelMuon_;
  UInt_t nElec_;
  UInt_t nSelElec_;

  UInt_t nJetsCentral10_;
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
  Bool_t  jettauNewDM_;

  Bool_t  jettauLooseIso_;
  Bool_t  jettauMediumIso_;
  Bool_t  jettauTightIso_;

  Bool_t  jettauVLooseMvaIso_;
  Bool_t  jettauLooseMvaIso_;
  Bool_t  jettauMediumMvaIso_;
  Bool_t  jettauTightMvaIso_;
  Bool_t  jettauVTightMvaIso_;

  Bool_t  jettauVLooseMva2017Iso_;
  Bool_t  jettauLooseMva2017Iso_;
  Bool_t  jettauMediumMva2017Iso_;
  Bool_t  jettauTightMva2017Iso_;
  Bool_t  jettauVTightMva2017Iso_;

  Bool_t jettauAntiMuonLoose3_;
  Bool_t jettauAntiMuonTight3_;

  Int_t   jettauDecay_;
  Int_t   jettauGenDecay_;
  Int_t   jettauGenMatch_;

  Bool_t jettauAntiElectronVLooseMVA6_;
  Bool_t jettauAntiElectronLooseMVA6_;
  
  Bool_t jettaubyVVVLooseDeepTau2017v2p1VSe_;
  Bool_t jettaubyVVLooseDeepTau2017v2p1VSe_;
  Bool_t jettaubyVLooseDeepTau2017v2p1VSe_;
  Bool_t jettaubyLooseDeepTau2017v2p1VSe_;
  Bool_t jettaubyMediumDeepTau2017v2p1VSe_;
  Bool_t jettaubyTightDeepTau2017v2p1VSe_;
  Bool_t jettaubyVTightDeepTau2017v2p1VSe_;
  Bool_t jettaubyVVTightDeepTau2017v2p1VSe_;

  Bool_t jettaubyVLooseDeepTau2017v2p1VSmu_;
  Bool_t jettaubyLooseDeepTau2017v2p1VSmu_;
  Bool_t jettaubyMediumDeepTau2017v2p1VSmu_;
  Bool_t jettaubyTightDeepTau2017v2p1VSmu_;

  Bool_t jettaubyVVVLooseDeepTau2017v2p1VSjet_;
  Bool_t jettaubyVVLooseDeepTau2017v2p1VSjet_;
  Bool_t jettaubyVLooseDeepTau2017v2p1VSjet_;
  Bool_t jettaubyLooseDeepTau2017v2p1VSjet_;
  Bool_t jettaubyMediumDeepTau2017v2p1VSjet_;
  Bool_t jettaubyTightDeepTau2017v2p1VSjet_;
  Bool_t jettaubyVTightDeepTau2017v2p1VSjet_;
  Bool_t jettaubyVVTightDeepTau2017v2p1VSjet_;

  Float_t jettauLeadingTrackPt_;
  Float_t jettauLeadingTrackEta_;
  Float_t jettauLeadingTrackPhi_;
  Float_t jettauLeadingTrackDz_;
  Float_t jettauLeadingTrackDxy_;

  UInt_t jetChargedMult_;
  UInt_t jetNeutralMult_;
  UInt_t jetChargedHadMult_;
  Float_t jetNeutralEMEnergyFraction_;
  Float_t jetNeutralHadEnergyFraction_;

  Float_t muonJetTauMass_;
  Float_t muonJetTauTrkMass_;

  Float_t jet2Pt_;
  Float_t jet2Eta_;
  Float_t jet2Phi_;

  UInt_t jet2ChargedMult_;
  UInt_t jet2NeutralMult_;
  UInt_t jet2ChargedHadMult_;
  Float_t jet2NeutralEMEnergyFraction_;
  Float_t jet2NeutralHadEnergyFraction_;

  Float_t jetmuonPt_;
  Float_t jetmuonEta_;
  Float_t jetmuonPhi_;

  Float_t jetelecPt_;
  Float_t jetelecEta_;
  Float_t jetelecPhi_;

  Float_t mueffweight;
  Float_t mutrigweight;

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
  
  TTree * ntuple_ = new TTree("NTuple","NTuple");

  ntuple_->Branch("event",&event_,"event/i"); 
  ntuple_->Branch("run",  &run_,  "run/i");

  ntuple_->Branch("puWeight",  &puWeight_,  "puWeight/F");
  ntuple_->Branch("genWeight", &genWeight_, "genWeight/F");
  ntuple_->Branch("weight",    &weight_,    "weight/F");
  ntuple_->Branch("mueffweight", &mueffweight, "mueffweight/F");
  ntuple_->Branch("mutrigweight",&mutrigweight,"mutrigweight/F");

  ntuple_->Branch("genHt",&genHt_,"genHt/F");
  ntuple_->Branch("genWPt",&genWPt_,"genWPt/F");
  ntuple_->Branch("nPartons",&nPartons_,"nPartons/I");
  ntuple_->Branch("nPartonsNLO",&nPartonsNLO_,"nPartonsNLO/I");

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
  ntuple_->Branch("nSelElec",&nSelElec_,"nSelElec/i");

  ntuple_->Branch("nJetsCentral10",&nJetsCentral10_,"nJetsCentral10/i");
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
  ntuple_->Branch("jettauNewDM",&jettauNewDM_,"jettauNewDM/O");

  ntuple_->Branch("jettauLooseIso", &jettauLooseIso_, "jettauLooseIso/O");
  ntuple_->Branch("jettauMediumIso",&jettauMediumIso_,"jettauMediumIso/O");
  ntuple_->Branch("jettauTightIso", &jettauTightIso_, "jettauTightIso/O");

  ntuple_->Branch("jettauVLooseMvaIso", &jettauVLooseMvaIso_, "jettauVLooseMvaIso/O");
  ntuple_->Branch("jettauLooseMvaIso", &jettauLooseMvaIso_, "jettauLooseMvaIso/O");
  ntuple_->Branch("jettauMediumMvaIso",&jettauMediumMvaIso_,"jettauMediumMvaIso/O");
  ntuple_->Branch("jettauTightMvaIso", &jettauTightMvaIso_, "jettauTightMvaIso/O");
  ntuple_->Branch("jettauVTightMvaIso", &jettauVTightMvaIso_, "jettauVTightMvaIso/O");

  ntuple_->Branch("jettauVLooseMva2017Iso", &jettauVLooseMva2017Iso_, "jettauVLooseMva2017Iso/O");
  ntuple_->Branch("jettauLooseMva2017Iso", &jettauLooseMva2017Iso_, "jettauLooseMva2017Iso/O");
  ntuple_->Branch("jettauMediumMva2017Iso",&jettauMediumMva2017Iso_,"jettauMediumMva2017Iso/O");
  ntuple_->Branch("jettauTightMva2017Iso", &jettauTightMva2017Iso_, "jettauTightMva2017Iso/O");
  ntuple_->Branch("jettauVTightMva2017Iso", &jettauVTightMva2017Iso_, "jettauVTightMva2017Iso/O");

  ntuple_->Branch("jettauAntiMuonLoose3",&jettauAntiMuonLoose3_,"jettauAntiMuonLoose3/O");
  ntuple_->Branch("jettauAntiMuonTight3",&jettauAntiMuonTight3_,"jettauAntiMuonTight3/O");

  ntuple_->Branch("jettauAntiElectronVLooseMVA6",&jettauAntiElectronVLooseMVA6_,"jettauAntiElectronVLooseMVA6/O");
  ntuple_->Branch("jettauAntiElectronLooseMVA6", &jettauAntiElectronLooseMVA6_, "jettauAntiElectronLooseMVA6/O");

  // **************** DeepTau discriminators **********************
  // **************************************************************

  // Deep antielectron discriminators
  ntuple_->Branch("jettaubyVVVLooseDeepTau2017v2p1VSe",
		  &jettaubyVVVLooseDeepTau2017v2p1VSe_,
		  "jettaubyVVVLooseDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyVVLooseDeepTau2017v2p1VSe",
		  &jettaubyVVLooseDeepTau2017v2p1VSe_,
		  "jettaubyVVLooseDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyVLooseDeepTau2017v2p1VSe",
		  &jettaubyVLooseDeepTau2017v2p1VSe_,
		  "jettaubyVLooseDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyLooseDeepTau2017v2p1VSe",
		  &jettaubyLooseDeepTau2017v2p1VSe_,
		  "jettaubyLooseDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyMediumDeepTau2017v2p1VSe",
		  &jettaubyMediumDeepTau2017v2p1VSe_,
		  "jettaubyMediumDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyTightDeepTau2017v2p1VSe",
		  &jettaubyTightDeepTau2017v2p1VSe_,
		  "jettaubyTightDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyVTightDeepTau2017v2p1VSe",
		  &jettaubyVTightDeepTau2017v2p1VSe_,
		  "jettaubyVTightDeepTau2017v2p1VSe/O");
  ntuple_->Branch("jettaubyVVTightDeepTau2017v2p1VSe",
		  &jettaubyVVTightDeepTau2017v2p1VSe_,
		  "jettaubyVVTightDeepTau2017v2p1VSe/O");

  ntuple_->Branch("jettaubyVVVLooseDeepTau2017v2p1VSjet",
		  &jettaubyVVVLooseDeepTau2017v2p1VSjet_,
		  "jettaubyVVVLooseDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyVVLooseDeepTau2017v2p1VSjet",
		  &jettaubyVVLooseDeepTau2017v2p1VSjet_,
		  "jettaubyVVLooseDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyVLooseDeepTau2017v2p1VSjet",
		  &jettaubyVLooseDeepTau2017v2p1VSjet_,
		  "jettaubyVLooseDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyLooseDeepTau2017v2p1VSjet",
		  &jettaubyLooseDeepTau2017v2p1VSjet_,
		  "jettaubyLooseDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyMediumDeepTau2017v2p1VSjet",
		  &jettaubyMediumDeepTau2017v2p1VSjet_,
		  "jettaubyMediumDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyTightDeepTau2017v2p1VSjet",
		  &jettaubyTightDeepTau2017v2p1VSjet_,
		  "jettaubyTightDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyVTightDeepTau2017v2p1VSjet",
		  &jettaubyVTightDeepTau2017v2p1VSjet_,
		  "jettaubyVTightDeepTau2017v2p1VSjet/O");
  ntuple_->Branch("jettaubyVVTightDeepTau2017v2p1VSjet",
		  &jettaubyVVTightDeepTau2017v2p1VSjet_,
		  "jettaubyVVTightDeepTau2017v2p1VSjet/O");

  ntuple_->Branch("jettaubyVLooseDeepTau2017v2p1VSmu",
		  &jettaubyVLooseDeepTau2017v2p1VSmu_,
		  "jettaubyVLooseDeepTau2017v2p1VSmu/O");
  ntuple_->Branch("jettaubyLooseDeepTau2017v2p1VSmu",
		  &jettaubyLooseDeepTau2017v2p1VSmu_,
		  "jettaubyLooseDeepTau2017v2p1VSmu/O");
  ntuple_->Branch("jettaubyMediumDeepTau2017v2p1VSmu",
		  &jettaubyMediumDeepTau2017v2p1VSmu_,
		  "jettaubyMediumDeepTau2017v2p1VSmu/O");
  ntuple_->Branch("jettaubyTightDeepTau2017v2p1VSmu",
		  &jettaubyTightDeepTau2017v2p1VSmu_,
		  "jettaubyTightDeepTau2017v2p1VSmu/O");


  ntuple_->Branch("jettauDecay",&jettauDecay_,"jettauDecay/I");
  ntuple_->Branch("jettauGenDecay",&jettauGenDecay_,"jettauGenDecay/I");
  ntuple_->Branch("jettauGenMatch",&jettauGenMatch_,"jettauGenMatch/I");

  ntuple_->Branch("jettauLeadingTrackPt",&jettauLeadingTrackPt_,"jettauLeadingTrackPt/F");
  ntuple_->Branch("jettauLeadingTrackEta",&jettauLeadingTrackEta_,"jettauLeadingTrackEta/F");
  ntuple_->Branch("jettauLeadingTrackPhi",&jettauLeadingTrackPhi_,"jettauLeadingTrackPhi/F");
  ntuple_->Branch("jettauLeadingTrackDz",&jettauLeadingTrackDz_,"jettauLeadingTrackDz/F");
  ntuple_->Branch("jettauLeadingTrackDxy",&jettauLeadingTrackDxy_,"jettauLeadingTrackDxy/F");

  ntuple_->Branch("muonjettauMass",&muonJetTauMass_,"muonjettauMass/F");
  ntuple_->Branch("muonjettautrkMass",&muonJetTauTrkMass_,"muonjettautrkMass/F");

  ntuple_->Branch("jetChargedMult",   &jetChargedMult_,   "jetChargedMult/i");
  ntuple_->Branch("jetNeutralMult",   &jetNeutralMult_,   "jetNeutralMult/i");
  ntuple_->Branch("jetChargedHadMult",&jetChargedHadMult_,"jetChargedHadMult/i");

  ntuple_->Branch("jetNeutralEMEnergyFraction", &jetNeutralEMEnergyFraction_, "jetNeutralEMEnergyFraction/F");
  ntuple_->Branch("jetNeutralHadEnergyFraction",&jetNeutralHadEnergyFraction_,"jetNeutralHadEnergyFraction/F");

  ntuple_->Branch("jetmuonPt",&jetmuonPt_,"jetmuonPt/F");
  ntuple_->Branch("jetmuonEta",&jetmuonEta_,"jetmuonEta/F");
  ntuple_->Branch("jetmuonPhi",&jetmuonPhi_,"jetmuonPhi/F");

  ntuple_->Branch("jetelecPt",&jetelecPt_,"jetelecPt/F");
  ntuple_->Branch("jetelecEta",&jetelecEta_,"jetelecEta/F");
  ntuple_->Branch("jetelecPhi",&jetelecPhi_,"jetelecPhi/F");

  ntuple_->Branch("jet2Pt", &jet2Pt_, "jet2Pt/F");
  ntuple_->Branch("jet2Eta",&jet2Eta_,"jet2Eta/F");
  ntuple_->Branch("jet2Phi",&jet2Phi_,"jet2Phi/F");

  ntuple_->Branch("jet2ChargedMult",   &jet2ChargedMult_,   "jet2ChargedMult/i");
  ntuple_->Branch("jet2NeutralMult",   &jet2NeutralMult_,   "jet2NeutralMult/i");
  ntuple_->Branch("jet2ChargedHadMult",&jet2ChargedHadMult_,"jet2ChargedHadMult/i");

  ntuple_->Branch("jet2NeutralEMEnergyFraction", &jet2NeutralEMEnergyFraction_, "jet2NeutralEMEnergyFraction/F");
  ntuple_->Branch("jet2NeutralHadEnergyFraction",&jet2NeutralHadEnergyFraction_,"jet2NeutralHadEnergyFraction/F");

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
  
  std::cout << "OK " << std::endl;

  // Official PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUDataFile,"read");
  TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PUMCFile, "read");
  TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
  TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get(PUMCHist);
  if (PUOfficial_mc->IsZombie()) {
    std::cout << "Histogram " << PUMCHist << "_pileup is not found in file " << PUMCFile << std::endl;
    exit(-1);
  }
  PUofficial->set_h_data(PUOfficial_data);
  PUofficial->set_h_MC(PUOfficial_mc);

  std::cout << "OK1 " << std::endl;

  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  ScaleFactor * SF_muonTrig = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));

  std::cout << "OK2 " << std::endl;

  // MEt filters
  std::vector<TString> metFlags; metFlags.clear();
  metFlags.push_back("Flag_HBHENoiseFilter");
  metFlags.push_back("Flag_HBHENoiseIsoFilter");
  //  metFlags.push_back("Flag_globalTightHalo2016Filter");
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

    TTree * _inittree = (TTree*)file_->Get(TString(initNtupleName));
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
	else {
	  Float_t GenWeight = 1;
	  if (genweight<0) GenWeight = -1;
	  histWeightsH->Fill(0.,GenWeight);
	}
      }
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

      genHt_ = analysisTree.genparticles_lheHt;
      genWPt_ = analysisTree.genparticles_lheWPt;
      metFilters_ = metFiltersPasses(analysisTree,metFlags);

      nPartons_ = int(analysisTree.genparticles_noutgoing);
      nPartonsNLO_ = int(analysisTree.genparticles_noutgoing_NLO);

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


      recoilRatio_ = -1;
      recoilDPhi_ = 0;

      recoilM_ = -1;
      recoilPt_ = -1;
      recoilEta_ = 0;
      recoilPhi_ = 0;

      nMuon_ = 0;
      nSelMuon_ = 0;
      nElec_ = 0;
      nSelElec_ = 0;

      nJetsCentral10_ = 0;
      nJetsCentral20_ = 0;
      nJetsCentral30_ = 0;

      nJetsForward20_ = 0;
      nJetsForward30_ = 0;

      jetPt_  = -9999;
      jetEta_ = -9999;
      jetPhi_ = -9999;

      jettauPt_ = 0;
      jettauEta_ = 0;
      jettauPhi_ = 0;

      jettauDM_ = false;
      jettauNewDM_ = false;
      jettauLooseIso_ = false;
      jettauMediumIso_ = false;
      jettauTightIso_ = false;

      jettauVLooseMvaIso_ = false;
      jettauLooseMvaIso_ = false;
      jettauMediumMvaIso_ = false;
      jettauTightMvaIso_ = false;
      jettauVTightMvaIso_ = false;

      jettauVLooseMva2017Iso_ = false;
      jettauLooseMva2017Iso_ = false;
      jettauMediumMva2017Iso_ = false;
      jettauTightMva2017Iso_ = false;
      jettauVTightMva2017Iso_ = false;

      jettauAntiMuonLoose3_ = false;
      jettauAntiMuonTight3_ = false;

      jettaubyVLooseDeepTau2017v2p1VSmu_ = false;
      jettaubyLooseDeepTau2017v2p1VSmu_ = false;
      jettaubyMediumDeepTau2017v2p1VSmu_ = false;
      jettaubyTightDeepTau2017v2p1VSmu_ = false;

      jettaubyVVVLooseDeepTau2017v2p1VSe_ = false;
      jettaubyVVLooseDeepTau2017v2p1VSe_ = false;
      jettaubyVLooseDeepTau2017v2p1VSe_ = false;
      jettaubyLooseDeepTau2017v2p1VSe_ = false;
      jettaubyMediumDeepTau2017v2p1VSe_ = false;
      jettaubyTightDeepTau2017v2p1VSe_ = false;
      jettaubyVTightDeepTau2017v2p1VSe_ = false;
      jettaubyVVTightDeepTau2017v2p1VSe_ = false;

      jettaubyVVVLooseDeepTau2017v2p1VSjet_ = false;
      jettaubyVVVLooseDeepTau2017v2p1VSjet_ = false;
      jettaubyVVLooseDeepTau2017v2p1VSjet_ = false;
      jettaubyVLooseDeepTau2017v2p1VSjet_ = false;
      jettaubyMediumDeepTau2017v2p1VSjet_ = false;
      jettaubyTightDeepTau2017v2p1VSjet_ = false;
      jettaubyVTightDeepTau2017v2p1VSjet_ = false;
      jettaubyVVTightDeepTau2017v2p1VSjet_ = false;

      jettauLeadingTrackPt_  = 0;
      jettauLeadingTrackEta_ = 0;
      jettauLeadingTrackPhi_ = 0;
      jettauLeadingTrackDz_  = -999;
      jettauLeadingTrackDxy_ = -999;

      jettauDecay_ = -1;
      jettauGenDecay_ = -1;
      jettauGenMatch_ = -1;

      jetChargedMult_ = 0;
      jetNeutralMult_ = 0;
      jetChargedHadMult_ = 0;

      jetNeutralEMEnergyFraction_ = 0;
      jetNeutralHadEnergyFraction_ = 0;

      jet2Pt_ = -9999;
      jet2Eta_ = -9999;
      jet2Phi_ = -9999;

      jet2ChargedMult_ = 0;
      jet2NeutralMult_ = 0;
      jet2ChargedHadMult_ = 0;

      jet2NeutralEMEnergyFraction_ = 0;
      jet2NeutralHadEnergyFraction_ = 0;

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
      if (isData) {
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
      }
      // ************************************
      // **** end accessing trigger info ****
      // ************************************


      // ***************************************************
      // accessing PF MET and changing momentum scale of met
      // ***************************************************
      float pfmet_ex = analysisTree.pfmetcorr_ex;
      float pfmet_ey = analysisTree.pfmetcorr_ey;
      if (jetES<0) {
	pfmet_ex = analysisTree.pfmetcorr_ex_JetEnDown;
	pfmet_ey = analysisTree.pfmetcorr_ey_JetEnDown;
      }
      else if (jetES>0) {
	pfmet_ex = analysisTree.pfmetcorr_ex_JetEnUp;
        pfmet_ey = analysisTree.pfmetcorr_ey_JetEnUp;
      }
      else if (unclusteredES<0) {
	pfmet_ex = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
	pfmet_ey = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
      }
      else if (unclusteredES>0) {
	pfmet_ex = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
        pfmet_ey = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
      }
      else {
	pfmet_ex = analysisTree.pfmetcorr_ex;
	pfmet_ey = analysisTree.pfmetcorr_ey;
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
      float phiTriggerMu = -1;
      float muonHt = 0;
      for (unsigned int imuon=0; imuon<analysisTree.muon_count; ++imuon) {
	analysisTree.muon_px[imuon] *= muonMomScale;
	analysisTree.muon_py[imuon] *= muonMomScale;
	analysisTree.muon_pz[imuon] *= muonMomScale;
	analysisTree.muon_pt[imuon] *= muonMomScale;
	if (analysisTree.muon_pt[imuon]<ptMuCut) continue;
	if (fabs(analysisTree.muon_eta[imuon])>etaMuCut) continue;
	bool passedId = analysisTree.muon_isMedium[imuon];
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
	    phiTriggerMu = analysisTree.muon_phi[imuon];
	    indexTriggerMu = int(imuon);
	  }
	}
      }

      nMuon_ = muonIndexes.size();
      nSelMuon_ = selMuonIndexes.size();

      float ptSecondMu  = -1;
      float etaSecondMu = -1;
      float phiSecondMu = -1;
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
	    phiSecondMu = analysisTree.muon_phi[indexMu];
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
      float phiLeadingMu = phiTriggerMu;
      float phiTrailingMu = phiSecondMu;

      if (ptTrailingMu>ptLeadingMu) {
	ptLeadingMu   = ptSecondMu;
	ptTrailingMu  = ptTriggerMu;
	etaLeadingMu  = etaSecondMu;
	etaTrailingMu = etaTriggerMu;
	phiLeadingMu  = phiSecondMu;
	phiTrailingMu = phiTriggerMu;
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
      std::vector<unsigned int> selEleIndexes; selEleIndexes.clear();
      std::vector<TLorentzVector>  lorentzVectorElectrons; lorentzVectorElectrons.clear();
      float elecHt = 0;
      for (unsigned int ielec=0; ielec<analysisTree.electron_count; ++ielec) {
	analysisTree.electron_px[ielec] *= eleMomScale;
        analysisTree.electron_py[ielec] *= eleMomScale;
        analysisTree.electron_pz[ielec] *= eleMomScale;
        analysisTree.electron_pt[ielec] *= eleMomScale;
	bool passedId = 
	  analysisTree.electron_mva_Loose_Iso_Fall17_v1[ielec] &&
          analysisTree.electron_pass_conversion[ielec] &&
          analysisTree.electron_nmissinginnerhits[ielec] <= 1;
	bool passedIdSel = analysisTree.electron_mva_wp80_Iso_Fall17_v1[ielec]>0.5 &&
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
	bool passedSelIso = relIso < isoSelEleCut;
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
	if (analysisTree.electron_pt[ielec]>ptEleCut &&
            fabs(analysisTree.electron_eta[ielec])<etaEleCut &&
            passedIdSel && passedIpCuts && passedSelIso)
	  selEleIndexes.push_back(ielec);
      }
      nElec_ = eleIndexes.size();
      nSelElec_ = selEleIndexes.size();

      // *********************************
      // **** end accessing electrons ****
      // *********************************     

      if (debug)
	std::cout << "end accessing electrons " << std::endl;

      // ************************
      // **** accessing jets ****
      // ************************
      TLorentzVector lorentzVectorAllJetsForMht; lorentzVectorAllJetsForMht.SetXYZT(0,0,0,0);
      std::vector<unsigned int> centralJets10Indexes; centralJets10Indexes.clear();
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
	
	// jetId
	bool isPFLooseJetId = looseJetiD(analysisTree,int(ijet));
	bool isPFTightJetId = tightJetiD(analysisTree,int(ijet));
	//	bool isPULooseJetId = looseJetPUiD(analysisTree,int(ijet));

	// jet four-vector
	TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[ijet],
					    analysisTree.pfjet_py[ijet],
					    analysisTree.pfjet_pz[ijet],
					    analysisTree.pfjet_e[ijet]);

	
	if (analysisTree.pfjet_pt[ijet]>10.&&absJetEta<2.5) 
	  centralJets10Indexes.push_back(ijet);

	if (analysisTree.pfjet_pt[ijet]<20) continue;


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
      nJetsCentral10_ = centralJets10Indexes.size();
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
      TLorentzVector lorentzVectorTau; lorentzVectorTau.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorTauTrk; lorentzVectorTauTrk.SetXYZT(0,0,0,0);
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
	  jettauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[taujetIndex] > 0.5;
	  jettauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  jettauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  jettauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[taujetIndex] > 0.5;
	  
	  jettauVLooseMvaIso_ = analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  jettauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[taujetIndex] > 0.5;
	  
	  jettauVLooseMva2017Iso_ = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[taujetIndex] > 0.5;
	  jettauLooseMva2017Iso_ = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[taujetIndex] > 0.5;
	  jettauMediumMva2017Iso_ = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[taujetIndex] > 0.5;
	  jettauTightMva2017Iso_ = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[taujetIndex] > 0.5;
	  jettauVTightMva2017Iso_ = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[taujetIndex] > 0.5;
	  
	  jettauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[taujetIndex] > 0.5;
	  jettauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[taujetIndex] > 0.5;
	  
	  jettauAntiElectronVLooseMVA6_ = analysisTree.tau_againstElectronVLooseMVA6[taujetIndex] > 0.5;
	  jettauAntiElectronLooseMVA6_ = analysisTree.tau_againstElectronLooseMVA6[taujetIndex] > 0.5;
	 
	  jettaubyVVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyVVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyVLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyLooseDeepTau2017v2p1VSjet_ = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyMediumDeepTau2017v2p1VSjet_ = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byTightDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[taujetIndex] > 0.5;
	  jettaubyVVTightDeepTau2017v2p1VSjet_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[taujetIndex] > 0.5;

	  jettaubyVVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyVVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyVLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyLooseDeepTau2017v2p1VSe_ = analysisTree.tau_byLooseDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyMediumDeepTau2017v2p1VSe_ = analysisTree.tau_byMediumDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyTightDeepTau2017v2p1VSe_ = analysisTree.tau_byTightDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVTightDeepTau2017v2p1VSe[taujetIndex] > 0.5;
	  jettaubyVVTightDeepTau2017v2p1VSe_ = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[taujetIndex] > 0.5;

	  jettaubyVLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[taujetIndex] > 0.5;
	  jettaubyLooseDeepTau2017v2p1VSmu_ = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[taujetIndex] > 0.5;
	  jettaubyMediumDeepTau2017v2p1VSmu_ = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[taujetIndex] > 0.5;
	  jettaubyTightDeepTau2017v2p1VSmu_ = analysisTree.tau_byTightDeepTau2017v2p1VSmu[taujetIndex] > 0.5;

	  jettauDecay_ = analysisTree.tau_decayMode[taujetIndex];
	  jettauGenDecay_ = analysisTree.tau_genDecayMode[taujetIndex];
	  
	  jettauLeadingTrackPt_ = PtoPt(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				     analysisTree.tau_leadchargedhadrcand_py[indexTau]);

	  jettauLeadingTrackEta_ = PtoEta(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				       analysisTree.tau_leadchargedhadrcand_py[indexTau],
				       analysisTree.tau_leadchargedhadrcand_pz[indexTau]);

	  jettauLeadingTrackPhi_ = PtoPhi(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				       analysisTree.tau_leadchargedhadrcand_py[indexTau]);
	  
	  jettauLeadingTrackDz_  = analysisTree.tau_leadchargedhadrcand_dz[indexTau];
	  jettauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];

	  lorentzVectorTau.SetXYZT(analysisTree.tau_px[indexTau],
				   analysisTree.tau_py[indexTau],
				   analysisTree.tau_pz[indexTau],
				   analysisTree.tau_e[indexTau]);

	  lorentzVectorTauTrk.SetXYZM(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				      analysisTree.tau_leadchargedhadrcand_py[indexTau],
				      analysisTree.tau_leadchargedhadrcand_pz[indexTau],
				      pionMass);

	  
	  jettauGenMatch_ = 6;
	  if (jettauGenDecay_>=0) jettauGenMatch_ = 5;
	  float minDR = 0.2;
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
		float dR = deltaR(jettauEta_,jettauPhi_,
				  etaGen,phiGen);
		if (dR<minDR) {
		  minDR = dR;
		  if (type1) jettauGenMatch_ = 1;
		  else if (type2) jettauGenMatch_ = 2;
		  else if (type3) jettauGenMatch_ = 3;
		  else if (type4) jettauGenMatch_ = 4;
		}
	      }
	    }
	    //	    std::cout << "Fake tau =>  pT = " << jettauPt_ << "   eta = " << jettauEta_ << "   phi = " << jettauPhi_ << "  genMatch = " << jettauGenMatch_ << std::endl;
	  }
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

      if (JetLeptonFake) {

      // ******************************
      // ********* ZJet selection *****
      // ******************************
	if (lorentzVectorZ.Pt()>1e-4) {

	  bool JetFound = false;
	  unsigned int jetFake = 0;
	  float deltaPhiMax = 0;
	  TLorentzVector lorentzFakeJet;
	  for (unsigned int iJet=0; iJet<centralJets10Indexes.size(); iJet++) {
	    unsigned int jetIndex = centralJets10Indexes.at(iJet);
	    TLorentzVector lorentzJet; lorentzJet.SetXYZT(analysisTree.pfjet_px[jetIndex],
						 analysisTree.pfjet_py[jetIndex],
						 analysisTree.pfjet_pz[jetIndex],
						 analysisTree.pfjet_e[jetIndex]);
	    float dRLead = deltaR(etaLeadingMu,phiLeadingMu,
			       lorentzJet.Eta(),lorentzJet.Phi());
	    if (dRLead<0.5) continue;
	    float dRTrail = deltaR(etaTrailingMu,phiTrailingMu,
			       lorentzJet.Eta(),lorentzJet.Phi());
	    if (dRTrail<0.5) continue;
	    
	    float deltaPhi = dPhiFromLV(lorentzVectorZ,lorentzJet);
	    if (deltaPhi>deltaPhiMax) {
	      JetFound = true;
	      jetFake = jetIndex;
	      deltaPhiMax = deltaPhi;
	      lorentzFakeJet = lorentzJet;
	    }

	  }
	  jetmuonPt_  = -9999;
	  jetmuonEta_ = -9999;
	  jetmuonPhi_ = -9999;
	  jetelecPt_  = -9999;
	  jetelecEta_ = -9999;
	  jetelecPhi_ = -9999;
	
	  if (JetFound) {

	    int fakeFound = 0;
	    unsigned int MuFake = 0;
	    unsigned int EleFake = 0;
	    float mindR = 0.5;
	    TLorentzVector lorentzLepFake;

	    for (unsigned int iSelMu=0; iSelMu<selMuonIndexes.size(); ++iSelMu) {
	      unsigned int indexMu = selMuonIndexes.at(iSelMu);
	      TLorentzVector lorentzMu; lorentzMu.SetXYZM(analysisTree.muon_px[indexMu],
							  analysisTree.muon_py[indexMu],
							  analysisTree.muon_pz[indexMu],
							  MuMass);
	      float dRl = deltaR(lorentzFakeJet.Eta(),lorentzFakeJet.Phi(),
				 lorentzMu.Eta(),lorentzMu.Phi());
	      if (dRl<mindR) {
		mindR = dRl;
		fakeFound = 1;
		MuFake = indexMu;
		lorentzLepFake = lorentzMu;
	      }
	    }

	    for (unsigned int iSelEle=0; iSelEle<selEleIndexes.size(); ++iSelEle) {
              unsigned int indexEle = selEleIndexes.at(iSelEle);
              TLorentzVector lorentzEle; lorentzEle.SetXYZM(analysisTree.electron_px[indexEle],
							    analysisTree.electron_py[indexEle],
							    analysisTree.electron_pz[indexEle],
							    electronMass);
              float dRl = deltaR(lorentzFakeJet.Eta(),lorentzFakeJet.Phi(),
                                 lorentzEle.Eta(),lorentzEle.Phi());
              if (dRl<mindR) {
                mindR = dRl;
                fakeFound = 2;
                EleFake = indexMu;
		lorentzLepFake = lorentzEle;
              }

            }
	    

	    recoilRatio_ =  lorentzLepFake.Pt() / lorentzVectorZ.Pt();
	    recoilDPhi_  = dPhiFromLV(lorentzVectorZ,lorentzLepFake);
	    isZJet = ptLeadingMu>ptLeadingMuCut_ZJet;
	    isZJet = isZJet && ptTrailingMu>ptTrailingMuCut_ZJet;
	    isZJet = isZJet && lorentzVectorZ.M()>ZMassLowerCut_ZJet && lorentzVectorZ.M()<ZMassUpperCut_ZJet;
	    isZJet = isZJet && recoilRatio_>ptJetZRatioLowerCut_ZJet && recoilRatio_<ptJetZRatioUpperCut_ZJet;
	    isZJet = isZJet && recoilDPhi_>deltaPhiZJetCut_ZJet;
	    if (isZJet) { 
	      float leadMuIso = SF_muonIdIso->get_ScaleFactor(ptLeadingMu, etaLeadingMu);
	      float trailMuIso = SF_muonIdIso->get_ScaleFactor(ptTrailingMu, etaTrailingMu);
	      mueffweight = leadMuIso*trailMuIso;
	      float leadMuTrig = SF_muonTrig->get_EfficiencyData(ptLeadingMu, etaLeadingMu);
	      float trailMuTrig = SF_muonTrig->get_EfficiencyData(ptTrailingMu, etaTrailingMu);
	      float leadMuTrigMC = SF_muonTrig->get_EfficiencyMC(ptLeadingMu, etaLeadingMu);
	      float trailMuTrigMC = SF_muonTrig->get_EfficiencyMC(ptTrailingMu, etaTrailingMu);
	      float trigData = 1 - (1-leadMuTrig)*(1-trailMuTrig);
	      float trigMC = 1 - (1-leadMuTrigMC)*(1-trailMuTrigMC);
	      mutrigweight = 1;
	      if (trigMC>0) mutrigweight = trigData / trigMC;
	      if (fakeFound==1) {
		jetmuonPt_ = lorentzLepFake.Pt();
		jetmuonEta_ = lorentzLepFake.Eta();
		jetmuonPhi_ = lorentzLepFake.Phi();
	      }
	      if (fakeFound==2) {
		jetelecPt_ = lorentzLepFake.Pt();
		jetelecEta_ = lorentzLepFake.Eta();
		jetelecPhi_ = lorentzLepFake.Phi();
	      }
	      jetPt_ = lorentzFakeJet.Pt();
	      jetEta_= lorentzFakeJet.Eta();
	      jetPhi_= lorentzFakeJet.Phi();

	      HtNoRecoil_     = Ht_     - ptTriggerMu - ptSecondMu;
	      SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu - ptSecondMu;
	      recoilM_   = lorentzVectorZ.M();
	      recoilPt_  = lorentzVectorZ.Pt();
	      recoilEta_ = lorentzVectorZ.Eta();
	      recoilPhi_ = lorentzVectorZ.Phi();
	      selection_ = 0;
	      ntuple_->Fill();
	      ZJetEvents++;
	    }
	  }
	}
      }
      else {
      
	// ***************************
	// ******* WJet selection ****
	// ***************************
	if (lorentzVectorW.Pt()>1e-4) {

	  muonJetTauMass_ = (lorentzVectorTriggerMu+lorentzVectorTau).M();
	  muonJetTauTrkMass_ = (lorentzVectorTriggerMu+lorentzVectorTauTrk).M();
	  
	  recoilRatio_ = jetPt_ / lorentzVectorW.Pt();
	  recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorJet);
	  
	  isWJet = ptTriggerMu>ptMuCut_WJet && jetPt_>20 && fabs(jetEta_)<2.4; 
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




//  LocalWords:  taujetIndex
