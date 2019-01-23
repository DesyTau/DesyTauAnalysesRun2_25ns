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

bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags, bool isData) {

  bool passed = true;
  unsigned int nFlags = metFlags.size();
  //  std::cout << "MEt filters : " << std::endl;
  for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
    TString flagName(it->first);
    //    std::cout << it->first << " : " << it->second << std::endl;
    for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
      if (flagName.Contains(metFlags[iFilter])) {
	if(flagName.Contains("eeBadScFilter") && !isData) continue;
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
  
  // trigger information
  const string metHTLName        = cfg.get<string>("MetHLTName");
  const string singleMuonHLTName = cfg.get<string>("SingleMuonHLTName");
  const string singleMuonHLTFilterName = cfg.get<string>("SingleMuonHLTFilterName");
  const string pfJet60HLTFilterName = cfg.get<string>("PFJet60HLTFilterName"); 
  const string pfJet80HLTFilterName = cfg.get<string>("PFJet80HLTFilterName"); 
  const string pfJet140HLTFilterName = cfg.get<string>("PFJet140HLTFilterName"); 
  const string pfJet200HLTFilterName = cfg.get<string>("PFJet200HLTFilterName"); 
  const string pfJet260HLTFilterName = cfg.get<string>("PFJet260HLTFilterName"); 
  const string pfJet320HLTFilterName = cfg.get<string>("PFJet320HLTFilterName"); 
  const string pfJet400HLTFilterName = cfg.get<string>("PFJet400HLTFilterName"); 
  const string pfJet450HLTFilterName = cfg.get<string>("PFJet450HLTFilterName"); 

  TString MetHLTName(metHTLName);
  TString SingleMuonHLTName(singleMuonHLTName);
  TString SingleMuonHLTFilterName(singleMuonHLTFilterName);
  TString PFJet60HLTFilterName(pfJet60HLTFilterName);
  TString PFJet80HLTFilterName(pfJet80HLTFilterName);
  TString PFJet140HLTFilterName(pfJet140HLTFilterName);
  TString PFJet200HLTFilterName(pfJet200HLTFilterName);
  TString PFJet260HLTFilterName(pfJet260HLTFilterName);
  TString PFJet320HLTFilterName(pfJet320HLTFilterName);
  TString PFJet400HLTFilterName(pfJet400HLTFilterName);
  TString PFJet450HLTFilterName(pfJet450HLTFilterName);

  // tau cuts
  const float ptTauCut  = cfg.get<float>("PtTauCut");
  const float etaTauCut = cfg.get<float>("EtaTauCut");

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
  
  // topological cuts (W*->tau+v)
  const float metCut_WTauNu                 = cfg.get<float>("MetCut_WTauNu");
  const float ptTauCut_WTauNu               = cfg.get<float>("PtTauCut_WTauNu");

  // topological cuts (W*->mu+v)
  const float ptMuCut_WMuNu               = cfg.get<float>("PtMuCut_WMuNu");
  const float metCut_WMuNu                = cfg.get<float>("MetCut_WMuNu");

  // topological cuts (W->muv+Jet)
  const float ptMuCut_WJet              = cfg.get<float>("PtMuCut_WJet");
  const float ptTauCut_WJet             = cfg.get<float>("PtTauCut_WJet");
  const float metCut_WJet               = cfg.get<float>("MetCut_WJet");

  // topological cuts (dijet)
  const float ptJetCut_DiJet              = cfg.get<float>("PtJetCut_DiJet");
  const float ptTauCut_DiJet              = cfg.get<float>("PtTauCut_DiJet");
  const float etaJetCut_DiJet             = cfg.get<float>("EtaJetCut_DiJet");
  const float deltaPhiTauJetCut_DiJet     = cfg.get<float>("DeltaPhiTauJetCut_DiJet");

  // topological cuts (trigger study)
  const float mtCut_Trig = cfg.get<float>("MtCut_Trig");

  // trigger eff filename
  const string trigEffFileName = cfg.get<string>("TrigEffFileName");

  // trigger weight cuts
  const float mhtNoMu_Trig              = cfg.get<float>("mhtNoMu_Trig");
  const float metNoMu_Trig              = cfg.get<float>("metNoMu_Trig");

  // momentum scales
  const float tauMomScale = cfg.get<float>("TauMomScale");
  const string tauDecayMode  = cfg.get<string>("TauDecayMode");
  const float muonMomScale = cfg.get<float>("MuonMomScale");
  const float eleMomScale = cfg.get<float>("EleMomScale");
  const int unclusteredES = cfg.get<int>("UnclusteredES");
  const int jetES = cfg.get<int>("JetES");

  const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  const string MuonTrigFile  = cfg.get<string>("MuonTrigEff");

  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");
  const string samplenameForPUHist = cfg.get<string>("SampleNameForPUHist");

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

  // ntuple variables

  UInt_t run_;
  UInt_t lumi_;
  UInt_t event_;
  
  Float_t puWeight_;
  Float_t genWeight_;
  Float_t trigWeight_;
  Float_t weight_;

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
  wntuple_->Branch("lepWE",&lepWE_,"lepWE/F");

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
  wntuple_->Branch("tauPz",  &tauPz_,  "tauPz/F");
  wntuple_->Branch("tauEta", &tauEta_, "tauEta/F");
  wntuple_->Branch("tauPhi", &tauPhi_, "tauPhi/F");

  wntuple_->Branch("genTauWPt",  &genTauWPt_,  "genTauWPt/F");
  wntuple_->Branch("genTauWEta", &genTauWEta_, "genTauWEta/F");
  wntuple_->Branch("genTauWPhi", &genTauWPhi_, "genTauWPhi/F");
  wntuple_->Branch("genTauWE", &genTauWE_, "genTauWE/F");

  TTree * ntuple_ = new TTree("NTuple","NTuple");

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
  ntuple_->Branch("WTauDecay",&wTauDecay_, "WTauDecay/I");

  ntuple_->Branch("lepWPt",&lepWPt_,"lepWPt/F");
  ntuple_->Branch("lepWEta",&lepWEta_,"lepWEta/F");
  ntuple_->Branch("lepWPhi",&lepWPhi_,"lepWPhi/F");
  ntuple_->Branch("lepWE",&lepWE_,"lepWE/F");

  ntuple_->Branch("nuWPt",&nuWPt_,"nuWPt/F");
  ntuple_->Branch("nuWEta",&nuWEta_,"nuWEta/F");
  ntuple_->Branch("nuWPhi",&nuWPhi_,"nuWPhi/F");

  TString numbers[20];

  for (int number=0; number<20; ++number) {
    char Number[20];
    int NN = number+1;
    if (NN<10)
      sprintf(Number,"%1i",NN);
    else
      sprintf(Number,"%2i",NN);

    numbers[number] = TString(Number);

  }

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

  TH1D * dRtauCentralJetH = new TH1D("dRtauCentralJetH","",50,0.,5.0);
  TH1D * dRtauForwardJetH = new TH1D("dRtauForwardJetH","",50,0.,5.0);

  Bool_t trigger_;
  Bool_t isWTrig_;
  Bool_t isZTrig_;
  Float_t metNoSelMu_;
  Float_t mhtNoSelMu_;
  UInt_t nMuonTrig_;
  UInt_t nSelMuonTrig_;
  Float_t dPhiMetMuon_;

  TTree * trigNTuple_ = new TTree("TriggerNTuple","TriggerNTuple");
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



  // project directory
  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Good run selection
  std::vector<Period> periods;
  string fullPathToJsonFile = cmsswBase + "/src/" + jsonFile;
  
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

  PileUp * PUofficial = new PileUp(); 
  if (!isData){
    // Official PU reweighting
    TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/"+PUDataFile,"read");
    TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/"+PUMCFile, "read");
    TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
    TString NamePUHistMC = samplenameForPUHist + "_pileup";
    TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get(NamePUHistMC);
    if( !PUOfficial_mc) {
      cout<<endl<<"MC pileup histogram "<<NamePUHistMC<<" does not exist in root file "<<PUMCFile<<". Exiting."<<endl<<endl;
      exit(-1);
    }
    PUofficial->set_h_data(PUOfficial_data);
    PUofficial->set_h_MC(PUOfficial_mc);
  }
  
  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  ScaleFactor * SF_muonTrig = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
  
  // Trigger efficiencies
  map<int,TGraphAsymmErrors*>  map_trigEffData;
  map<int,TGraphAsymmErrors*>  map_trigEffMC;
  TFile * trigEffFile = new TFile(TString(cmsswBase)+"/src/"+trigEffFileName,"read");
  if(trigEffFile->IsZombie()){
    cout<<"Trigger file "<<trigEffFileName<<" does not exists. Exiting."<<endl;
    exit(-1);
  }
  TIter next(trigEffFile->GetListOfKeys());
  TKey *key = 0;

  while ((key = (TKey*)next()))
    {
      TClass *c = gROOT->GetClass(key->GetClassName());
      if (!c->InheritsFrom("TGraphAsymmErrors")) continue;
      TGraphAsymmErrors *g = (TGraphAsymmErrors*) key->ReadObj();
      TString gName = g->GetName();
      Ssiz_t pos    = gName.First("To");
      Int_t  len    = gName.Length();
      int upperBound   = atoi( (TString) gName(pos+2,len));
      if(gName.Contains("data"))    map_trigEffData.insert( std::make_pair( upperBound , (TGraphAsymmErrors*) g->Clone() ) );
      else if(gName.Contains("mc")) map_trigEffMC.insert(   std::make_pair( upperBound , (TGraphAsymmErrors*) g->Clone() ) );
      delete g;
    }
  delete key;

  // MEt filters
  std::vector<TString> metFlags; metFlags.clear();
  metFlags.push_back("Flag_HBHENoiseFilter");
  metFlags.push_back("Flag_HBHENoiseIsoFilter");
  metFlags.push_back("Flag_globalTightHalo2016Filter");
  metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  metFlags.push_back("Flag_goodVertices"); 
  if (isData)
     metFlags.push_back("Flag_eeBadScFilter");
  //metFlags.push_back("Flag_muonBadTrackFilter");
  //metFlags.push_back("Flag_chargedHadronTrackResolutionFilter");
  metFlags.push_back("Flag_BadChargedCandidateFilter");
  metFlags.push_back("Flag_BadPFMuonFilter");
  metFlags.push_back("Flag_ecalBadCalibFilter");
  
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

      if(debug){
	cout<<"muon_count = "<<analysisTree.muon_count<<endl;
	cout<<"electron_count = "<<analysisTree.electron_count<<endl;
	cout<<"track_count = "<<analysisTree.track_count<<endl;
	cout<<"tau_count = "<<analysisTree.tau_count<<endl;
	cout<<"gentau_count = "<<analysisTree.gentau_count<<endl;
	cout<<"trigobject_count = "<<analysisTree.trigobject_count<<endl;
	cout<<"l1tau_count = "<<analysisTree.l1tau_count<<endl;
	cout<<"l1egamma_count = "<<analysisTree.l1egamma_count<<endl;
	cout<<"l1muon_count = "<<analysisTree.l1muon_count<<endl;
	cout<<"pfjet_count = "<<analysisTree.pfjet_count<<endl;
	cout<<"genparticles_count = "<<analysisTree.genparticles_count<<endl;
	cout<<"dimuon_count = "<<analysisTree.dimuon_count<<endl;
      }

      // ***************************
      // initialize ntuple variables
      // ***************************
      run_ = analysisTree.event_run;
      lumi_ = analysisTree.event_luminosityblock;
      event_ = analysisTree.event_nr;
      nVert_ = analysisTree.primvertex_count;

      if (event_<=0) {
	std::cout << "Event : " << event_ << std::endl;
	std::cout << "From NTuple = " << analysisTree.event_nr << std::endl;
      }

      weight_ = 1;
      genWeight_ = 1;
      trigWeight_ = 1;
      puWeight_ = 1;

      genHt_ = analysisTree.genparticles_lheHt;

      trig_ = false;
      metFilters_ = true;

      wMass_ = -1;
      wPt_ = -1;
      wEta_ = 0;
      wPhi_ = 0;
      wDecay_ = -1;
      wTauDecay_ = -1;

      lepWPt_ = -1;
      lepWEta_ = 0;
      lepWPhi_ = 0;
      lepWE_   = 0;

      nuWPt_ = -1;
      nuWEta_ = 0;
      nuWPhi_ = 0;

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

      mueffweight = 1;
      mutrigweight = 1;

      // booleans 
      bool isWJet = false;
      bool isWTauNu = false;
      bool isWMuNu  = false;
      bool isDiJet = false;

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

      if (debug) {
	std::cout << "Run = " << analysisTree.event_nr << "    Event = " << analysisTree.event_run << std::endl; 
	std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
	std::cout << "Number of taus          = " << analysisTree.tau_count << std::endl;
	std::cout << "Number of jets          = " << analysisTree.pfjet_count << std::endl;
	std::cout << "Number of muons         = " << analysisTree.muon_count << std::endl;
	std::cout << "Number of electrons     = " << analysisTree.electron_count << std::endl;
      }	

      // Set MC relevant variables
      if (!isData) {
	if (analysisTree.genweight<0) genWeight_ = -1;
	else                          genWeight_ = 1;

	weight_ *= genWeight_;

	npartons_ = analysisTree.genparticles_noutgoing;
	npartonsNLO_ = analysisTree.genparticles_noutgoing_NLO;
	lheWPt_   = analysisTree.genparticles_lheWPt;
      }
      histWeightsH->Fill(double(0.),double(genWeight_));
      
      // **********************************
      // *** Analysis of generator info ***
      // **********************************
      int indexW  = -1;
      //int indexNu = -1; // nu from W
      int indexMu = -1; // muon from W
      int indexE  = -1; // elec from W
      int indexTau = -1; // tau from W
      //int indexTauE = -1; // W->tau->e
      //int indexTauMu = -1; // W->tau->mu
      vector<TLorentzVector> gentauLV; gentauLV.clear();
      vector<int> gentauDecay; gentauDecay.clear();
      vector<TLorentzVector> genmuonLV; genmuonLV.clear();
      vector<TLorentzVector> genelecLV; genelecLV.clear();
      vector<TLorentzVector> gentaumuonLV; gentaumuonLV.clear();
      vector<TLorentzVector> gentauelecLV; gentauelecLV.clear();
      TLorentzVector wmuonLV; wmuonLV.SetXYZT(0,0,0,0);
      TLorentzVector welecLV; welecLV.SetXYZT(0,0,0,0);
      TLorentzVector wgenvistauLV;  wgenvistauLV.SetXYZT(0,0,0,0);
      TLorentzVector wgentauLV;  wgentauLV.SetXYZT(0,0,0,0);
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
	      //indexNu = igen;
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
	    //if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->mu
	    //  indexTauMu = igen;
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
            //if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->e
	    //  indexTauE = igen;
	  }
	}

	for (unsigned int igentau=0; igentau<analysisTree.gentau_count;++igentau) {
	  TLorentzVector GenVisTau; GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
						      analysisTree.gentau_visible_py[igentau],
						      analysisTree.gentau_visible_pz[igentau],
						      analysisTree.gentau_visible_e[igentau]);
	  TLorentzVector GenTau; GenTau.SetXYZT(analysisTree.gentau_px[igentau],
						analysisTree.gentau_py[igentau],
						analysisTree.gentau_pz[igentau],
						analysisTree.gentau_e[igentau]);

	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau] ) { // W/Z->tau  
	    gentauLV.push_back(GenVisTau);
	    gentauDecay.push_back(analysisTree.gentau_decayMode[igentau]);
	    indexTau = igentau;
	    wgenvistauLV = GenVisTau;
	    wgentauLV = GenTau;
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
	  lepWE_   = wmuonLV.E();
	}
	else if (indexE>=0) {
	  wDecay_ = 1;
	  lepWPt_  = welecLV.Pt();
          lepWEta_ = welecLV.Eta();
          lepWPhi_ = welecLV.Phi();
          lepWE_   = welecLV.E();
	}
	else if (indexTau>=0) {
	  lepWPt_  = wgenvistauLV.Pt();
          lepWEta_ = wgenvistauLV.Eta();
          lepWPhi_ = wgenvistauLV.Phi();
          lepWE_   = wgenvistauLV.E();
	  
	  genTauWPt_  = wgentauLV.Pt();
          genTauWEta_ = wgentauLV.Eta();
          genTauWPhi_ = wgentauLV.Phi();
          genTauWE_   = wgentauLV.E();
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
      //      bool isSingleMuonHLT = false;
      bool isMetHLT = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(MetHLTName)) {
	  if (it->second==1) isMetHLT = true;
	}
      }
      trigger_ = isMetHLT;
      trig_ = isMetHLT;

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

      for (unsigned int i=0; i<analysisTree.run_hltfilters->size(); ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==SingleMuonHLTFilterName) {
	  nSingleMuonHLTFilter = i;
	  isSingleMuonHLTFilter = true;
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
	if (HLTFilter==PFJet200HLTFilterName) {
	  nPFJet200HLTFilter = i;
	  isPFJet200HLTFilter = true;
	}
	if (HLTFilter==PFJet260HLTFilterName) {
	  nPFJet260HLTFilter = i;
	  isPFJet260HLTFilter = true;
	}
	if (HLTFilter==PFJet320HLTFilterName) {
	  nPFJet320HLTFilter = i;
	  isPFJet320HLTFilter = true;
	}
	if (HLTFilter==PFJet400HLTFilterName) {
	  nPFJet400HLTFilter = i;
	  isPFJet400HLTFilter = true;
	}
	if (HLTFilter==PFJet450HLTFilterName) {
	  nPFJet450HLTFilter = i;
	  isPFJet450HLTFilter = true;
	}
      }
      //      if (isData) {
      if (!isSingleMuonHLTFilter) {
	std::cout << "HLT filter " << SingleMuonHLTFilterName << " not found" << std::endl;
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
      if (!isPFJet200HLTFilter) {
	std::cout << "HLT filter " << PFJet200HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet260HLTFilter) {
	std::cout << "HLT filter " << PFJet260HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet320HLTFilter) {
	std::cout << "HLT filter " << PFJet320HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet400HLTFilter) {
	std::cout << "HLT filter " << PFJet400HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      if (!isPFJet450HLTFilter) {
	std::cout << "HLT filter " << PFJet450HLTFilterName << " not found" << std::endl;
	exit(-1);
      }
      // ************************************
      // **** end accessing trigger info ****
      // ************************************


      // ***************************************************
      // accessing PF MET and changing momentum scale of met
      // ***************************************************


      float pfmetcorr_ex = analysisTree.pfmetcorr_ex;
      float pfmetcorr_ey = analysisTree.pfmetcorr_ey;

      if (!isData) {
	if (jetES<0) {
	  pfmetcorr_ex = analysisTree.pfmetcorr_ex_JetEnDown;
	  pfmetcorr_ey = analysisTree.pfmetcorr_ey_JetEnDown;
	}
	else if (jetES>0) {
	  pfmetcorr_ex = analysisTree.pfmetcorr_ex_JetEnUp;
	  pfmetcorr_ey = analysisTree.pfmetcorr_ey_JetEnUp;
	}
	else if (unclusteredES<0) {
	  pfmetcorr_ex = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
	  pfmetcorr_ey = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
	}
	else if (unclusteredES>0) {
	  pfmetcorr_ex = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
	  pfmetcorr_ey = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
	}
	else {
	  pfmetcorr_ex = analysisTree.pfmetcorr_ex;
	  pfmetcorr_ey = analysisTree.pfmetcorr_ey;
	}
      }

      met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
      if (met_<1e-4) met_ = 1e-4;
      metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
      TLorentzVector lorentzVectorMet; lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);

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
      metNoMu_    =  (lorentzVectorMet+lorentzVectorAllMuons).Pt();
      metNoSelMu_ =  (lorentzVectorMet+lorentzVectorAllSelMuons).Pt();

      nMuon_ = muonIndexes.size();
      nMuonTrig_ = nMuon_;
      nSelMuon_ = selMuonIndexes.size();
      nSelMuonTrig_ = nSelMuon_;

      float ptSecondMu  = -1;
      //float etaSecondMu = -1;
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
   muonPz_  = lorentzVectorTriggerMu.Pz();
	muonEta_ = lorentzVectorTriggerMu.Eta();
	muonPhi_ = lorentzVectorTriggerMu.Phi();
	muonQ_   = int(analysisTree.muon_charge[indexTriggerMu]);
	pfmetcorr_ex = pfmetcorr_ex + lorentzVectorTriggerMu.Px() - lorentzVectorTriggerMu.Px();
	pfmetcorr_ey = pfmetcorr_ey + lorentzVectorTriggerMu.Py() - lorentzVectorTriggerMu.Py();
	met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
	metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
	lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);
	mtmuon_  = mT(lorentzVectorTriggerMu,lorentzVectorMet);
	dPhiMetMuon_ = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet); 
	lorentzVectorW = lorentzVectorTriggerMu + lorentzVectorMet;
	for (unsigned int iMu = 0; iMu < selMuonIndexes.size(); ++iMu) {
	  int indexMu = int(selMuonIndexes.at(iMu));
	  if (indexMu==indexTriggerMu) continue;
	  float netcharge = analysisTree.muon_charge[indexTriggerMu]*analysisTree.muon_charge[indexMu];
	  if (netcharge>0) continue;
	  if (analysisTree.muon_pt[indexMu]>ptSecondMu) {
	    ptSecondMu = analysisTree.muon_pt[indexMu];
	    //etaSecondMu = analysisTree.muon_eta[indexMu];
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
	}
	isWTrig_ = mtmuon_ >  mtCut_Trig;
      }
      float ptLeadingMu = ptTriggerMu;
      float ptTrailingMu = ptSecondMu;

      if (ptTrailingMu>ptLeadingMu) {
	ptLeadingMu = ptSecondMu;
	ptTrailingMu = ptTriggerMu;
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
	  analysisTree.electron_cutId_veto_Summer16[ielec] &&
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
      std::vector<bool> jets60trigger; jets60trigger.clear();
      std::vector<bool> jets80trigger; jets80trigger.clear();
      std::vector<bool> jets140trigger; jets140trigger.clear();
      std::vector<bool> jets200trigger; jets200trigger.clear();
      std::vector<bool> jets260trigger; jets260trigger.clear();
      std::vector<bool> jets320trigger; jets320trigger.clear();
      std::vector<bool> jets400trigger; jets400trigger.clear();
      std::vector<bool> jets450trigger; jets450trigger.clear();

      for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {

	// Clean for EEnoise jets
	if( analysisTree.pfjet_pt[ijet] < 50 && fabs(analysisTree.pfjet_eta[ijet]) < 3.139 && fabs(analysisTree.pfjet_eta[ijet]) > 2.65) continue;

	// Scale for sys uncertainties
	float scaleJ = 1;

	if (jetES<0)      scaleJ = 1.0 - analysisTree.pfjet_jecUncertainty[ijet];
	else if (jetES>0) scaleJ = 1.0 + analysisTree.pfjet_jecUncertainty[ijet];
	else 	          scaleJ = 1.0;

	analysisTree.pfjet_px[ijet] *= scaleJ;
	analysisTree.pfjet_py[ijet] *= scaleJ;
	analysisTree.pfjet_pz[ijet] *= scaleJ;
	analysisTree.pfjet_pt[ijet] *= scaleJ;
	analysisTree.pfjet_e[ijet]  *= scaleJ;

	float absJetEta = fabs(analysisTree.pfjet_eta[ijet]);

	if (absJetEta>5.2) continue;
	if (analysisTree.pfjet_pt[ijet]<20.0) continue; 



	// jetID : accept only jets with tight id for 2017 analysis
	bool isPFTightJetId = tightJetiD_2017(analysisTree,int(ijet));
	if (!isPFTightJetId) continue;

	// jet four-vector
	TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[ijet],
					    analysisTree.pfjet_py[ijet],
					    analysisTree.pfjet_pz[ijet],
					    analysisTree.pfjet_e[ijet]);

	// counting jets for Mht
	lorentzVectorAllJetsForMht += jetLV;

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
	bool trigMatch60 = false;
	bool trigMatch80 = false;
	bool trigMatch140 = false;
	bool trigMatch200 = false;
	bool trigMatch260 = false;
	bool trigMatch320 = false;
	bool trigMatch400 = false;
	bool trigMatch450 = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;

	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   trigMatch60 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   trigMatch80 = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) trigMatch140 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) trigMatch200 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet260HLTFilter]&&isPFJet260HLTFilter)) trigMatch260 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet320HLTFilter]&&isPFJet320HLTFilter)) trigMatch320 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet400HLTFilter]&&isPFJet400HLTFilter)) trigMatch400 = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet450HLTFilter]&&isPFJet450HLTFilter)) trigMatch450 = true;
	}
	if (!isData) {
	  trigMatch60 = true;
	  trigMatch80 = true;
	  trigMatch140 = true;
	  trigMatch200 = true;
	  trigMatch260 = true;
	  trigMatch320 = true;
	  trigMatch400 = true;
	  trigMatch450 = true;
	}
	triggerJetsIndexes.push_back(ijet);
	jets60trigger.push_back(trigMatch60);
	jets80trigger.push_back(trigMatch80);
	jets140trigger.push_back(trigMatch140);
	jets200trigger.push_back(trigMatch200);
	jets260trigger.push_back(trigMatch260);
	jets320trigger.push_back(trigMatch320);
	jets400trigger.push_back(trigMatch400);
	jets450trigger.push_back(trigMatch450);
      }
      nJetsCentral20_ = centralJets20Indexes.size();
      nJetsCentral30_ = centralJets30Indexes.size();
      nJetsForward20_ = forwardJets20Indexes.size();
      nJetsForward30_ = forwardJets30Indexes.size();
      JetHt_     = htCentral30 + htForward30;
      SoftJetHt_ = htCentral20 + htForward30;
      Ht_        = JetHt_     + muonHt + elecHt;
      SoftHt_    = SoftJetHt_ + muonHt + elecHt;
      mht_        = (lorentzVectorAllJetsForMht).Pt();
      mhtNoMu_    = (lorentzVectorAllJetsForMht - lorentzVectorAllMuons).Pt();
      mhtNoSelMu_ = (lorentzVectorAllJetsForMht - lorentzVectorAllSelMuons).Pt();
      TLorentzVector lorentzVectorJet; lorentzVectorJet.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorJet2; lorentzVectorJet2.SetXYZT(0,0,0,0);
      if (nJetsCentral30_>0) {
	unsigned int indexJet0 = centralJets30Indexes.at(0);
	jetPt_ = analysisTree.pfjet_pt[indexJet0];
	jetEta_ = analysisTree.pfjet_eta[indexJet0];
	jetPhi_ = analysisTree.pfjet_phi[indexJet0];
	jetBtag_ = analysisTree.pfjet_btag[indexJet0][0];
	jetChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet0];
	jetNeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet0];
	jetChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet0];
	jetNeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];
	jetNeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];  
	pfJet60_ = false;
	pfJet80_ = false;
	pfJet140_ = false;
	pfJet200_ = false;
	pfJet260_ = false;
	pfJet320_ = false;
	pfJet400_ = false;
	pfJet450_ = false;
	lorentzVectorJet.SetXYZT(analysisTree.pfjet_px[indexJet0],
				 analysisTree.pfjet_py[indexJet0],
				 analysisTree.pfjet_pz[indexJet0],
				 analysisTree.pfjet_e[indexJet0]);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet0],analysisTree.pfjet_phi[indexJet0],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;
	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pfJet60_  = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pfJet80_  = true; 
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pfJet140_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) pfJet200_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet260HLTFilter]&&isPFJet260HLTFilter)) pfJet260_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet320HLTFilter]&&isPFJet320HLTFilter)) pfJet320_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet400HLTFilter]&&isPFJet400HLTFilter)) pfJet400_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet450HLTFilter]&&isPFJet450HLTFilter)) pfJet450_ = true;
	}
      }
      if (nJetsCentral30_>1) {
	unsigned int indexJet1 = centralJets30Indexes.at(1);
	jet2Pt_ = analysisTree.pfjet_pt[indexJet1];
	jet2Eta_ = analysisTree.pfjet_eta[indexJet1];
	jet2Phi_ = analysisTree.pfjet_phi[indexJet1];
	jet2Btag_ = analysisTree.pfjet_btag[indexJet1][0];
	jet2ChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet1];
	jet2NeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet1];
	jet2ChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet1];
	jet2NeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];
	jet2NeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];  
	pf2Jet60_ = false;
	pf2Jet80_ = false;
	pf2Jet140_ = false;
	pf2Jet200_ = false;
	lorentzVectorJet2.SetXYZT(analysisTree.pfjet_px[indexJet1],
				  analysisTree.pfjet_py[indexJet1],
				  analysisTree.pfjet_pz[indexJet1],
				  analysisTree.pfjet_e[indexJet1]);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet1],analysisTree.pfjet_phi[indexJet1],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>0.5) continue;
	  if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pf2Jet60_  = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pf2Jet80_  = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pf2Jet140_ = true;
	  if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) pf2Jet200_ = true;
	}
      }
      

      // ****************************
      // **** end accessing jets ****
      // ****************************

      if (debug)
	std::cout << "end accessing jets" << std::endl;

      // ************************
      // **** accessing taus ****
      // ************************
      std::vector<unsigned int> tauIndexes; tauIndexes.clear();
      std::vector<unsigned int> tau20Indexes; tau20Indexes.clear();
      std::vector<unsigned int> tau30Indexes; tau30Indexes.clear();
      std::vector<int> tauGenMatchDecay; tauGenMatchDecay.clear();
      for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) { // loop over taus

	if( (tauDecayMode == "1prong0pizeros"     && analysisTree.tau_decayMode[itau]==0) ||
	    (tauDecayMode == "1prongUpTo4pizeros" && analysisTree.tau_decayMode[itau]>=1 && analysisTree.tau_decayMode[itau]<=4) ||
	    (tauDecayMode == "3prong0pizeros"     && analysisTree.tau_decayMode[itau]==10) ){

	  analysisTree.tau_px[itau]   *= tauMomScale;
	  analysisTree.tau_py[itau]   *= tauMomScale;
	  analysisTree.tau_pz[itau]   *= tauMomScale;
	  analysisTree.tau_pt[itau]   *= tauMomScale;
	  analysisTree.tau_e[itau]    *= tauMomScale;
	  analysisTree.tau_mass[itau] *= tauMomScale;
	}
	
	if (fabs(analysisTree.tau_eta[itau])>2.4) continue; // loose eta cut
	if (analysisTree.tau_pt[itau]<20.) continue; // loose pt cut

	//	if (!(analysisTree.tau_decayModeFindingNewDMs[itau]>0.5||analysisTree.tau_decayModeFinding[itau]>0.5)) 
	//	  cout << "Zombie tau " << itau
	//	       << "   pt(tau) = " << analysisTree.tau_pt[itau]
	//	       << "   eta(tau) = " << analysisTree.tau_eta[itau] << endl;

	float dZ = fabs(analysisTree.tau_vertexz[itau]-analysisTree.primvertex_z);
	if (dZ>1e-4) continue; // dz criterion

	bool foundByDecayMode = analysisTree.tau_decayModeFindingNewDMs[itau]>0.5 || analysisTree.tau_decayModeFinding[itau]>0.5;
	if (!foundByDecayMode) continue; // DM finding

	// finding matching jet -->
	//	int matchedJetIndex = -1;
	//	float dRmin = 0.4;
	//	for (unsigned int ijet=0; ijet<centralJets20Indexes.size(); ijet++) {
	//	  unsigned int indexJet = centralJets20Indexes.at(ijet);
	//	  float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
	//			    analysisTree.pfjet_eta[indexJet],analysisTree.pfjet_phi[indexJet]);
	//	  if (dR<dRmin) {
	//	    dRmin = dR;
	//	    matchedJetIndex = int(indexJet);
	//	  }
	//	}
	//	if (matchedJetIndex>=0) {
	//	  cout << "jet matching tau " << itau 
	//	       << "  pt(tau) = " << analysisTree.tau_pt[itau]
	//	       << "  pt(jet) = " << analysisTree.pfjet_pt[matchedJetIndex] 
	//	       << "  dR(tau,jet) = " << dRmin << endl;
	//	}
	//	else {
	//	  cout << "no matching jet found for tau " << itau 
	//	       << "  pt(tau) = " << analysisTree.tau_pt[itau] << endl;
	//	}

	// finding matching mu -->
	int matchedMuIndex = -1;
	float dRmin = 0.4;
	for (unsigned int im=0; im<muonIndexes.size(); im++) {
	  unsigned int indexMu = muonIndexes.at(im);
	  float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
			    analysisTree.muon_eta[indexMu],analysisTree.muon_phi[indexMu]);
	  if (dR<dRmin) {
	    dRmin = dR;
	    matchedMuIndex = int(indexMu);
	  }
	}
	//	if (matchedMuIndex>=0) {
	//	  cout << "muon matching tau " << itau 
	//	       << "  pt(tau) = " << analysisTree.tau_pt[itau]
	//	       << "  pt(mu) = " << analysisTree.muon_pt[matchedMuIndex] 
	//	       << "  dR(tau,mu) = " << dRmin << endl;
	//	}
	

	// finding matching e -->
	int matchedEleIndex = -1;
	dRmin = 0.4;
	for (unsigned int ie=0; ie<eleIndexes.size(); ie++) {
	  unsigned int indexEle = eleIndexes.at(ie);
	  float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
			    analysisTree.electron_eta[indexEle],analysisTree.electron_phi[indexEle]);
	  if (dR<dRmin) {
	    dRmin = dR;
	    matchedEleIndex = int(indexEle);
	  }
	}
	//	if (matchedEleIndex>=0) {
	//	  cout << "electron matching tau " << itau 
	//	       << "  pt(tau) = " << analysisTree.tau_pt[itau]
	//	       << "  pt(e) = " << analysisTree.electron_pt[matchedEleIndex] 
	//	       << "  dR(tau,e) = " << dRmin << endl;
	//	}


	if (analysisTree.tau_pt[itau]>20.&&matchedMuIndex<0&&matchedEleIndex<0) {
	  tau20Indexes.push_back(itau);
	}

	if (analysisTree.tau_pt[itau]>30.&&matchedMuIndex<0&&matchedEleIndex<0) {
	  tau30Indexes.push_back(itau);
	}

	if (analysisTree.tau_pt[itau]>ptTauCut&&fabs(analysisTree.tau_eta[itau])<etaTauCut&&matchedMuIndex<0&&matchedEleIndex<0) { 
	  tauIndexes.push_back(itau);
	  int genMatchDecay = -1;
	  float dRmin = 0.2;
	  //	  std::cout << "size = " << gentauLV.size() << std::endl;
	  for (unsigned int igentau=0; igentau<gentauLV.size(); ++igentau) {
	    TLorentzVector genTauLV = gentauLV.at(igentau);
	    int decaymode = gentauDecay.at(igentau);
	    float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
			      genTauLV.Eta(),genTauLV.Phi());
	    if (dR<dRmin) {
	      dRmin = dR;
	      genMatchDecay = decaymode;
	    }
	  }
	  tauGenMatchDecay.push_back(genMatchDecay);
	}
      } // end loop over taus
      nTaus20_ = tau20Indexes.size();
      nTaus30_ = tau30Indexes.size();
      nSelTaus_ = tauIndexes.size();
      TLorentzVector lorentzVectorTau; lorentzVectorTau.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorTauJet; lorentzVectorTauJet.SetXYZT(0,0,0,0);
      if (nSelTaus_>0) {

	unsigned int indexTau = tauIndexes.at(0);
	lorentzVectorTau.SetXYZM(analysisTree.tau_px[indexTau],
				 analysisTree.tau_py[indexTau],
				 analysisTree.tau_pz[indexTau],
				 analysisTree.tau_mass[indexTau]);
	pfmetcorr_ex = pfmetcorr_ex + lorentzVectorTau.Px() - lorentzVectorTau.Px();
	pfmetcorr_ey = pfmetcorr_ey + lorentzVectorTau.Py() - lorentzVectorTau.Py();
	met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
	metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
	lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);
	mttau_ = mT(lorentzVectorTau,lorentzVectorMet);
	mtgen_ = mT(wgentauLV,wnuLV);
	tauPt_ = analysisTree.tau_pt[indexTau];
   tauPz_ = analysisTree.tau_pz[indexTau];
	tauEta_ = analysisTree.tau_eta[indexTau];
	tauPhi_ = analysisTree.tau_phi[indexTau];
	tauMass_ = analysisTree.tau_mass[indexTau];
	tauQ_ = int(analysisTree.tau_charge[indexTau]);
	tauNtrk1_ = analysisTree.tau_ntracks_pt1[indexTau];
	tauNtrk05_ = analysisTree.tau_ntracks_pt05[indexTau];
	tauNtrk08_ = analysisTree.tau_ntracks_pt08[indexTau];

	tauLeadingTrackPt_ = PtoPt(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				   analysisTree.tau_leadchargedhadrcand_py[indexTau]);

	tauLeadingTrackEta_ = PtoEta(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				    analysisTree.tau_leadchargedhadrcand_py[indexTau],
				    analysisTree.tau_leadchargedhadrcand_pz[indexTau]);

	tauLeadingTrackPhi_ = PtoPhi(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				     analysisTree.tau_leadchargedhadrcand_py[indexTau]);

	tauLeadingTrackDz_  = analysisTree.tau_leadchargedhadrcand_dz[indexTau];
	tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];

	tauDecay_ = analysisTree.tau_decayMode[indexTau];
	tauGenDecay_ = analysisTree.tau_genDecayMode[indexTau];
	tauGenMatchDecay_ = tauGenMatchDecay.at(0);

	if (tauDecay_<0) tauDecay_ = -1;
	if (tauGenDecay_<0) tauGenDecay_ = -1;
	if (tauGenMatchDecay_<0) tauGenMatchDecay_ = -1;

	tauGenMatch_ = 6;
	if (tauGenMatchDecay_>=0) tauGenMatch_ = 5;
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
	      float dR = deltaR(tauEta_,tauPhi_,
				etaGen,phiGen);
	      if (dR<minDR) {
		minDR = dR;
		if (type1) tauGenMatch_ = 1;
		else if (type2) tauGenMatch_ = 2;
		else if (type3) tauGenMatch_ = 3;
		else if (type4) tauGenMatch_ = 4;
	      }

	    }
	  }
	}

	tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
	tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;

   tauIso_ = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[indexTau];
	tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;

	tauVLooseMvaIso_ = analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauVVTightMvaIso_ = analysisTree.tau_byVVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;

	tauVVLooseMva2017v2Iso_ = analysisTree.tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauVLooseMva2017v2Iso_ = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauLooseMva2017v2Iso_ = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauMediumMva2017v2Iso_ = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauTightMva2017v2Iso_ = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauVTightMva2017v2Iso_ = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	tauVVTightMva2017v2Iso_ = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;

	tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
	tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;

	tauAntiElectronVLooseMVA6_ = analysisTree.tau_againstElectronVLooseMVA6[indexTau] > 0.5;
	tauAntiElectronLooseMVA6_  = analysisTree.tau_againstElectronLooseMVA6[indexTau] > 0.5;
	tauAntiElectronMediumMVA6_  = analysisTree.tau_againstElectronMediumMVA6[indexTau] > 0.5;
	tauAntiElectronTightMVA6_  = analysisTree.tau_againstElectronTightMVA6[indexTau] > 0.5;
	tauAntiElectronVTightMVA6_ = analysisTree.tau_againstElectronVTightMVA6[indexTau] > 0.5;

	// finding matching jet
	bool jetFound = false;
	float dRmin = 0.4;
	int indexMatchingJet = -1;
	for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
	  TLorentzVector lorentzVectorJ; lorentzVectorJ.SetXYZT(analysisTree.pfjet_px[ijet],
								analysisTree.pfjet_py[ijet],
								analysisTree.pfjet_pz[ijet],
								analysisTree.pfjet_e[ijet]);
	  float drJetTau = deltaR(lorentzVectorJ.Eta(),lorentzVectorJ.Phi(),
				  lorentzVectorTau.Eta(),lorentzVectorTau.Phi());

	  if (drJetTau<dRmin) {
	    dRmin = drJetTau;
	    jetFound = true;
	    indexMatchingJet = ijet;
	    lorentzVectorTauJet = lorentzVectorJ;
	  }

	}
	if (!jetFound) {
	  lorentzVectorTauJet = lorentzVectorTau;
	  continue;
	}

	tauJetPt_  = lorentzVectorTauJet.Pt();
	tauJetEta_ = lorentzVectorTauJet.Eta();
	tauJetPhi_ = lorentzVectorTauJet.Phi();
	tauJetTightId_ = tightJetiD_2017(analysisTree,indexMatchingJet);

      }
      // ****************************
      // **** end accessing taus ****
      // ****************************
      
      if (debug)
	std::cout << "end of accessing taus" << std::endl;

      // ****************************
      // ****** trigger weight ******
      // ****************************
      float trigEffData = 1.0;
      float trigEffMC   = 1.0;

      trigWeight_ = 1;

      for (auto const& it : map_trigEffMC)
	{
	  if( mhtNoMu_ < it.first)
	    {
	      trigEffMC    =  it.second->Eval(metNoMu_);
	      trigEffData  =  map_trigEffData[it.first]->Eval(metNoMu_);
	      break;
	    }
	}

      if(trigEffMC !=0 ) trigWeight_ = trigEffData / trigEffMC;
      if(trigWeight_ < 0) trigWeight_=0;
      if(mhtNoMu_<mhtNoMu_Trig || metNoMu_< metNoMu_Trig) trigWeight_=0;

      if (debug) {
	cout << "MetNoMu  = " << metNoMu_
	     << "  MhtNoMu  = " << mhtNoMu_
	     << "  trigWeight = " << trigWeight_ << endl;
      }
      weight_ *= trigWeight_;

      // setting met filters
      metFilters_ = metFiltersPasses(analysisTree,metFlags,isData);

      // ********************************
      // **** filling trigger ntuple ****
      // ********************************
      if (ptTriggerMu>ptTrigMuCut && metNoMu_>metNoMu_Trig && mhtNoMu_>mhtNoMu_Trig) {
	trigNTuple_->Fill();
	TrigEvents++;
      }
      
      // ***************************
      // ******* WJet selection ****
      // ***************************
      if (lorentzVectorW.Pt()>1e-4) {
	recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorTau);
	recoilJetRatio_ = lorentzVectorTauJet.Pt()/lorentzVectorW.Pt();
	recoilJetDPhi_ = dPhiFromLV(lorentzVectorW,lorentzVectorTauJet);
	isWJet = ptTriggerMu>ptMuCut_WJet; 
	isWJet = isWJet && met_>metCut_WJet;
	isWJet = isWJet && tauPt_>ptTauCut_WJet;
	isWJet = isWJet && nMuon_ == 1;
	isWJet = isWJet && nElec_ == 0;
	isWJet = isWJet && nSelTaus_ == 1;
	isWJet = isWJet && nJetsCentral30_ == 1;
	isWJet = isWJet && abs(muonEta_)<2.1;
	isWJet = isWJet && tauDM_;


	if (isWJet) {
	  if (!isData) {
	    mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
	    mutrigweight = SF_muonTrig->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
	  }
	  HtNoRecoil_     = Ht_     - ptTriggerMu;
	  SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu;
	  recoilM_   = lorentzVectorW.M();
          recoilPt_  = lorentzVectorW.Pt();
          recoilEta_ = lorentzVectorW.Eta();
          recoilPhi_ = lorentzVectorW.Phi();
	  dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
	  selection_ = 1;
	  ntuple_->Fill();
	  WJetEvents++;
	}
      }
      
      // ********************************
      // ******* W*->MuNu selection *****
      // ********************************
      if (lorentzVectorMet.Pt()>1e-4) {
	recoilDPhi_  = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet);
	recoilJetRatio_ = -1;
	recoilJetDPhi_  = 0;
	isWMuNu = ptTriggerMu>ptMuCut_WMuNu;
	isWMuNu = isWMuNu && met_>metCut_WMuNu;
	isWMuNu = isWMuNu && nMuon_ == 1;
	isWMuNu = isWMuNu && nElec_ == 0;
	isWMuNu = isWMuNu && nSelTaus_ == 0;
	isWMuNu = isWMuNu && nJetsCentral30_ == 0;
	isWMuNu = isWMuNu && abs(muonEta_)<2.1;

	if (isWMuNu) {
	  if (!isData) {
	    mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
	    mutrigweight = SF_muonTrig->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
	  }
	  HtNoRecoil_     = Ht_;
	  SoftHtNoRecoil_ = SoftHt_;
	  recoilM_   = lorentzVectorMet.M();
	  recoilPt_  = lorentzVectorMet.Pt();
	  recoilEta_ = lorentzVectorMet.Eta();
	  recoilPhi_ = lorentzVectorMet.Phi();
	  selection_ = 2;
	  ntuple_->Fill();
	  WMuNuEvents++;
	}
      }

      // ********************************
      // ****** W*->TauNu selection *****
      // ******************************** 
      if (lorentzVectorMet.Pt()>1e-4) {
	recoilDPhi_  = dPhiFromLV(lorentzVectorTau,lorentzVectorMet);
	recoilJetRatio_ = lorentzVectorTauJet.Pt()/lorentzVectorMet.Pt();
	recoilJetDPhi_ = dPhiFromLV(lorentzVectorMet,lorentzVectorTauJet);
	isWTauNu = met_>metCut_WTauNu;
	isWTauNu = isWTauNu && tauPt_>ptTauCut_WTauNu;
	isWTauNu = isWTauNu && nSelTaus_ >= 1;
	isWTauNu = isWTauNu && nJetsCentral30_==1;
	isWTauNu = isWTauNu && nMuon_ == 0;
	isWTauNu = isWTauNu && nElec_ == 0;
	isWTauNu = isWTauNu && tauDM_;
	isWTauNu = isWTauNu && trigger_;

	if (isWTauNu) {
	  HtNoRecoil_     = Ht_;
	  SoftHtNoRecoil_ = SoftHt_;
	  recoilM_   = lorentzVectorMet.M();
          recoilPt_  = lorentzVectorMet.Pt();
          recoilEta_ = lorentzVectorMet.Eta();
          recoilPhi_ = lorentzVectorMet.Phi();
	  dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
	  selection_ = 3;
	  ntuple_->Fill();
	  WTauNuEvents++;
	}
      }


      // *********************************
      // ****** Jet+Tau selection ********
      // *********************************
      isDiJet = nSelTaus_>0 && triggerJetsIndexes.size()>0;
      bool foundJetTauPair = false;
      if (isDiJet) {
	for (unsigned int iTau=0; iTau<tauIndexes.size(); ++iTau) { // loop over taus
	  
	  unsigned int indexTau = tauIndexes.at(iTau);

	  if(nMuon_!=0)          continue;
	  if(nElec_!=0)          continue;
	  if(nSelTaus_!=1)       continue;
	  if(nJetsCentral30_!=2) continue;
	  if(analysisTree.tau_pt[indexTau]<ptTauCut_DiJet) continue;
	  if(!pfJet60_ && !pfJet80_ && !pfJet140_ && !pfJet200_ && !pfJet260_ && !pfJet320_ && !pfJet400_ && !pfJet450_ ) continue;

	  TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[indexTau],
					      analysisTree.tau_py[indexTau],
					      analysisTree.tau_pz[indexTau],
					      analysisTree.tau_mass[indexTau]);
	  // finding recoiling jet 
	  int acceptedJetIndex = -1;
	  int acceptedJetDirIndex = -1;
	  TLorentzVector recoilJetLV; recoilJetLV.SetXYZT(0,0,0,0);
	  float dPhiTauJetMax = deltaPhiTauJetCut_DiJet;

	  for (unsigned int iJet=0; iJet<triggerJetsIndexes.size(); ++iJet) { // loop over jets
	    unsigned int indexJet = triggerJetsIndexes.at(iJet);
	    TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[indexJet],
						analysisTree.pfjet_py[indexJet],
						analysisTree.pfjet_pz[indexJet],
						analysisTree.pfjet_e[indexJet]);
	    if (jetLV.Pt()<ptJetCut_DiJet) continue;
	    if (fabs(jetLV.Eta())>etaJetCut_DiJet) continue;
	    float ptTauJetRatio = tauLV.Pt() / jetLV.Pt();
	    //if (ptTauJetRatio<ptTauJetRatioLowerCut_DiJet) continue;
	    //if (ptTauJetRatio>ptTauJetRatioUpperCut_DiJet) continue;
	    float dPhiTauJet = dPhiFromLV(tauLV,jetLV);
	    if (dPhiTauJet>dPhiTauJetMax) {
	      acceptedJetIndex = int(indexJet);
	      acceptedJetDirIndex = int(iJet);
	      dPhiTauJetMax = dPhiTauJet;
	      recoilJetLV = jetLV;
	    }
	  } // end loop over jet
	  if (acceptedJetIndex>=0) { // recoil jet found

	    foundJetTauPair =  true;

	    mttau_ = mT(tauLV,lorentzVectorMet);
	    mtgen_ = mT(wgentauLV,wnuLV);
	    tauPt_ = analysisTree.tau_pt[indexTau];
       tauPz_ = analysisTree.tau_pz[indexTau];
	    tauEta_ = analysisTree.tau_eta[indexTau];
	    tauPhi_ = analysisTree.tau_phi[indexTau];
	    tauMass_ = analysisTree.tau_mass[indexTau];
	    tauQ_ = int(analysisTree.tau_charge[indexTau]);
	    tauNtrk1_ = analysisTree.tau_ntracks_pt1[indexTau];
	    tauNtrk05_ = analysisTree.tau_ntracks_pt05[indexTau];
	    tauNtrk08_ = analysisTree.tau_ntracks_pt08[indexTau];

	    // finding matching jet
	    bool jetFound = false;
	    float dRmin = 0.4;
	    int indexMatchingJet = -1;
	    for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
	      TLorentzVector lorentzVectorJ; lorentzVectorJ.SetXYZT(analysisTree.pfjet_px[ijet],
								    analysisTree.pfjet_py[ijet],
								    analysisTree.pfjet_pz[ijet],
								    analysisTree.pfjet_e[ijet]);
	      float drJetTau = deltaR(lorentzVectorJ.Eta(),lorentzVectorJ.Phi(),
				      tauEta_,tauPhi_);
	      
	      if (drJetTau<dRmin) {
		dRmin = drJetTau;
		jetFound = true;
		indexMatchingJet = ijet;
		lorentzVectorTauJet = lorentzVectorJ;
	      }
	      
	    }
	    if (!jetFound) {
	      lorentzVectorTauJet = lorentzVectorTau;
	      continue;
	    }

	    tauJetPt_  = lorentzVectorTauJet.Pt();
	    tauJetEta_ = lorentzVectorTauJet.Eta();
	    tauJetPhi_ = lorentzVectorTauJet.Phi();
	    tauJetTightId_ = tightJetiD_2017(analysisTree,indexMatchingJet);

	    tauLeadingTrackPt_ = PtoPt(analysisTree.tau_leadchargedhadrcand_px[indexTau],
				       analysisTree.tau_leadchargedhadrcand_py[indexTau]);
	    
	    tauLeadingTrackEta_ = PtoEta(analysisTree.tau_leadchargedhadrcand_px[indexTau],
					 analysisTree.tau_leadchargedhadrcand_py[indexTau],
					 analysisTree.tau_leadchargedhadrcand_pz[indexTau]);
	    
	    tauLeadingTrackPhi_ = PtoPhi(analysisTree.tau_leadchargedhadrcand_px[indexTau],
					 analysisTree.tau_leadchargedhadrcand_py[indexTau]);
	    
	    tauLeadingTrackDz_  = analysisTree.tau_leadchargedhadrcand_dz[indexTau];
	    tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];
	    
	    tauDecay_ = analysisTree.tau_decayMode[indexTau];
	    tauGenDecay_ = analysisTree.tau_genDecayMode[indexTau];
	    tauGenMatchDecay_ = tauGenMatchDecay.at(iTau);

	    if (tauDecay_<0) tauDecay_ = -1;
	    if (tauGenDecay_<0) tauGenDecay_ = -1;
	    if (tauGenMatchDecay_<0) tauGenMatchDecay_ = -1;

	    pfJet60_ = jets60trigger.at(acceptedJetDirIndex);
	    pfJet80_ = jets80trigger.at(acceptedJetDirIndex);
	    pfJet140_ = jets140trigger.at(acceptedJetDirIndex);
	    pfJet200_ = jets200trigger.at(acceptedJetDirIndex);
	    pfJet260_ = jets260trigger.at(acceptedJetDirIndex);
	    pfJet320_ = jets320trigger.at(acceptedJetDirIndex);
	    pfJet400_ = jets400trigger.at(acceptedJetDirIndex);
	    pfJet450_ = jets450trigger.at(acceptedJetDirIndex);

	    tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
	    tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;

       tauIso_ =  analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[indexTau];
	    tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	    tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	    tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;

	    tauVLooseMvaIso_ = analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauVVTightMvaIso_ = analysisTree.tau_byVVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;

	    tauVVLooseMva2017v2Iso_ = analysisTree.tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauVLooseMva2017v2Iso_ = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauLooseMva2017v2Iso_ = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauMediumMva2017v2Iso_ = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauTightMva2017v2Iso_ = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauVTightMva2017v2Iso_ = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
	    tauVVTightMva2017v2Iso_ = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
 
	    tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
	    tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;

	    recoilDPhi_ = dPhiFromLV(tauLV,recoilJetLV);
	    recoilJetRatio_ = lorentzVectorTauJet.Pt()/recoilJetLV.Pt();
	    recoilJetDPhi_ = dPhiFromLV(lorentzVectorTauJet,recoilJetLV);
	    recoilM_ = recoilJetLV.M();
	    recoilPt_ = recoilJetLV.Pt();
	    recoilEta_ = recoilJetLV.Eta();
	    recoilPhi_ = recoilJetLV.Phi();
	    HtNoRecoil_     = Ht_     - recoilJetLV.Pt();
	    SoftHtNoRecoil_ = SoftHt_ - recoilJetLV.Pt();
	    dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
	    selection_ = 4;
	    if(!tauDM_) continue;
	    ntuple_->Fill();
	  }
	} // END: loop over tau
      }
      if (foundJetTauPair) JetTauEvents++;

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
  std::cout << "Total number of trigger events    = " << TrigEvents << std::endl;
  std::cout << "Total number of W+Jet events      = " << WJetEvents << std::endl;
  std::cout << "Total number of W->muv events     = " << WMuNuEvents << std::endl;
  std::cout << "Total number of W->tauv events    = " << WTauNuEvents << std::endl;
  std::cout << "Total number of single jet events = " << SingleJetEvents << std::endl;
  std::cout << "Total number of jet+tau events    = " << JetTauEvents << std::endl;
  std::cout << "Total number of dijet events      = " << DiJetEvents << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;

  //  delete fakerate;
  //  delete fakerateE;
  
}



