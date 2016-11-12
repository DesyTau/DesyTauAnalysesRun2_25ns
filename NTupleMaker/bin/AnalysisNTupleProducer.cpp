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
void FakeRateMva(float pt, float * mean, float * error, int &iso) {

  // 80X Run2016BCD
  if (pt<150) {
    mean[0]  = 0.0725; error[0] = 0.0033;
    mean[1]  = 0.0422; error[1] = 0.0024;
    mean[2]  = 0.0271; error[2] = 0.0020;
    mean[3]  = 0.019; error[3] = 0.004;
    iso = 0;
  }
  else if (pt>=150&&pt<200) {
    mean[0]  = 0.0818; error[0] = 0.0090;
    mean[1]  = 0.0492; error[1] = 0.0065;
    mean[2]  = 0.0290; error[2] = 0.0051;
    mean[3]  = 0.009; error[3] = 0.008;
    iso = 1;
  }
  else {
    mean[0]  = 0.0640; error[0] = 0.0101;
    mean[1]  = 0.0369; error[1] = 0.0079;
    mean[2]  = 0.0217; error[2] = 0.0065;
    mean[3]  = 0.029; error[3] = 0.027;
    iso = 2;
  }

  // 80X Run2016B
  /*
  if (pt<150) {
    mean[0]  = 0.072; error[0] = 0.004;
    mean[1]  = 0.043; error[1] = 0.003;
    mean[2]  = 0.026; error[2] = 0.003;
    mean[3]  = 0.019; error[3] = 0.004;
    iso = 0;
  }
  else if (pt>=150&&pt<200) {
    mean[0]  = 0.065; error[0] = 0.011;
    mean[1]  = 0.037; error[1] = 0.008;
    mean[2]  = 0.021; error[2] = 0.006;
    mean[3]  = 0.009; error[3] = 0.008;
    iso = 1;
  }
  else {
    mean[0]  = 0.069; error[0] = 0.014;
    mean[1]  = 0.042; error[1] = 0.012;
    mean[2]  = 0.026; error[2] = 0.011;
    mean[3]  = 0.029; error[3] = 0.027;
    iso = 2;
  }
  */
  // 76X MEt filters + pt(leadtrk) > 4 GeV
  //  if (pt<150) {
  //    mean[0]  = 0.080; error[0] = 0.009;
  //    mean[1]  = 0.048; error[1] = 0.007;
  //    mean[2]  = 0.035; error[2] = 0.006;
  //    iso = 0;
  //  }
  //  else if (pt>=150&&pt<200) {
  //    mean[0]  = 0.061; error[0] = 0.020;
  //    mean[1]  = 0.049; error[1] = 0.019;
  //    mean[2]  = 0.017; error[2] = 0.012;
  //    iso = 1;
  //  }
  //  else {
  //    mean[0]  = 0.078; error[0] = 0.035;
  //    mean[1]  = 0.072; error[1] = 0.035;
  //    mean[2]  = 0.032; error[2] = 0.027;
  //    iso = 2;
  //  }

}

void FakeRate(float pt, float * mean, float * error, int &iso) {

  // 80X Run2016BCD
  if (pt<150) {
    mean[0]  = 0.0949; error[0] = 0.0037;
    mean[1]  = 0.0596; error[1] = 0.0029;
    mean[2]  = 0.0376; error[2] = 0.0024;
    iso = 0;
  }
  else if (pt>=150&&pt<200) {
    mean[0]  = 0.1153; error[0] = 0.0096;
    mean[1]  = 0.0735; error[1] = 0.0075;
    mean[2]  = 0.0508; error[2] = 0.0071;
    iso = 1;
  }
  else {
    mean[0]  = 0.0924; error[0] = 0.0121;
    mean[1]  = 0.0629; error[1] = 0.0102;
    mean[2]  = 0.0453; error[2] = 0.0092;
    iso = 2;
  } 

  // 80X Run2016B
  /*
  if (pt<150) {
    mean[0]  = 0.092; error[0] = 0.004;
    mean[1]  = 0.056; error[1] = 0.004;
    mean[2]  = 0.035; error[2] = 0.003;
    iso = 0;
  }
  else if (pt>=150&&pt<200) {
    mean[0]  = 0.103; error[0] = 0.013;
    mean[1]  = 0.060; error[1] = 0.010;
    mean[2]  = 0.041; error[2] = 0.009;
    iso = 1;
  }
  else {
    mean[0]  = 0.109; error[0] = 0.019;
    mean[1]  = 0.076; error[1] = 0.017;
    mean[2]  = 0.058; error[2] = 0.015;
    iso = 2;
  } 
  */
   // 76X MEt filters + pt(leadtrk) > 4 GeV
   //   if (pt<150) {
   //    mean[0]  = 0.100; error[0] = 0.010;
   //    mean[1]  = 0.061; error[1] = 0.008;
   //    mean[2]  = 0.041; error[2] = 0.007;
   //    iso = 0;
   //  }
   //  else if (pt>=150&&pt<200) {
   //    mean[0]  = 0.097; error[0] = 0.025;
   //    mean[1]  = 0.081; error[1] = 0.023;
   //    mean[2]  = 0.045; error[2] = 0.018;
   //    iso = 1;
   //  }
   //  else {
   //    mean[0]  = 0.072; error[0] = 0.035;
   //    mean[1]  = 0.066; error[1] = 0.035;
   //    mean[2]  = 0.049; error[2] = 0.031;
   //    iso = 2;
   //  } 

   // 76X no MEt filters, 
   /*  if (pt<150) {
       mean[0]  = 0.110; error[0] = 0.010;
       mean[1]  = 0.067; error[1] = 0.008;
       mean[2]  = 0.046; error[2] = 0.007;
       iso = 0;
       }
       else if (pt>=150&&pt<200) {
       mean[0]  = 0.100; error[0] = 0.024;
       mean[1]  = 0.081; error[1] = 0.022;
       mean[2]  = 0.045; error[2] = 0.017;
       iso = 1;
       }
       else {
       mean[0]  = 0.113; error[0] = 0.047;
       mean[1]  = 0.087; error[1] = 0.043;
       mean[2]  = 0.068; error[2] = 0.040;
       iso = 2;
       } 
       
       // 74X
       if (pt<150) {
       mean[0]  = 0.100; error[0] = 0.007;
       mean[1]  = 0.067; error[1] = 0.006;
       mean[2]  = 0.043; error[2] = 0.005;
       iso = 0;
       }
       else if (pt>=150&&pt<200) {
       mean[0]  = 0.100; error[0] = 0.020;
       mean[1]  = 0.080; error[1] = 0.018;
       mean[2]  = 0.052; error[2] = 0.015;
       iso = 1;
       }
       else {
       mean[0]  = 0.086; error[0] = 0.044;
       mean[1]  = 0.028; error[1] = 0.020;
       mean[2]  = 0.022; error[2] = 0.022;
       iso = 2;
       }
   */

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
  const bool applyHT100Cut = cfg.get<bool>("ApplyHT100Cut");
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
  
  // topological cuts (W*->tau+v)
  const float ptTauMetRatioLowerCut_WTauNu  = cfg.get<float>("PtTauMetRatioLowerCut_WTauNu"); 
  const float ptTauMetRatioUpperCut_WTauNu  = cfg.get<float>("PtTauMetRatioUpperCut_WTauNu"); 
  const float deltaPhiTauMetCut_WTauNu      = cfg.get<float>("DeltaPhiTauMetCut_WTauNu");
  const float metCut_WTauNu                 = cfg.get<float>("MetCut_WTauNu");

  // topological cuts (W*->mu+v)
  const float ptMuCut_WMuNu               = cfg.get<float>("PtMuCut_WMuNu");
  const float ptMuMetRatioLowerCut_WMuNu  = cfg.get<float>("PtMuMetRatioLowerCut_WMuNu");
  const float ptMuMetRatioUpperCut_WMuNu  = cfg.get<float>("PtMuMetRatioUpperCut_WMuNu");
  const float deltaPhiMuMetCut_WMuNu      = cfg.get<float>("DeltaPhiMuMetCut_WMuNu");
  const float metCut_WMuNu                = cfg.get<float>("MetCut_WMuNu");

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

  // weighting (tau fake rate)
  const string tauFakeRateFileName = cfg.get<string>("TauFakeRateFileName");
  const string tauAntiLFakeRateFileName = cfg.get<string>("TauAntiLFakeRateFileName");

  // trigger eff filename
  const string trigEffFileName = cfg.get<string>("TrigEffFileName");
  const string trigEffFileName74X = cfg.get<string>("TrigEffFileName74X");

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
  Float_t trigWeight_;
  Float_t trigWeight74X_;
  Float_t weight_;

  UInt_t nVert_;

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

  Float_t nuWPt_;
  Float_t nuWEta_;
  Float_t nuWPhi_;

  Float_t tauMatchJetPt_;
  Float_t tauMatchJetEta_;
  Float_t tauMatchJetPhi_;

  //  Float_t effLoose_;
  //  Float_t effMedium_;
  //  Float_t effTight_;

  //  Float_t effAntiLLoose_;
  //  Float_t effAntiLMedium_;
  //  Float_t effAntiLTight_;

  //  Float_t fakeLoose_;
  //  Float_t fakeMedium_;
  //  Float_t fakeTight_;

  Float_t fakeAntiLLoose_;
  Float_t fakeAntiLMedium_;
  Float_t fakeAntiLTight_;

  Float_t fakeAntiLLooseUp1_;
  Float_t fakeAntiLMediumUp1_;
  Float_t fakeAntiLTightUp1_;

  Float_t fakeAntiLLooseUp2_;
  Float_t fakeAntiLMediumUp2_;
  Float_t fakeAntiLTightUp2_;

  Float_t fakeAntiLLooseUp3_;
  Float_t fakeAntiLMediumUp3_;
  Float_t fakeAntiLTightUp3_;

  Float_t fakeAntiLLooseMva_;
  Float_t fakeAntiLMediumMva_;
  Float_t fakeAntiLTightMva_;
  Float_t fakeAntiLVTightMva_;

  Float_t fakeAntiLLooseMvaUp1_;
  Float_t fakeAntiLMediumMvaUp1_;
  Float_t fakeAntiLTightMvaUp1_;
  Float_t fakeAntiLVTightMvaUp1_;

  Float_t fakeAntiLLooseMvaUp2_;
  Float_t fakeAntiLMediumMvaUp2_;
  Float_t fakeAntiLTightMvaUp2_;
  Float_t fakeAntiLVTightMvaUp2_;

  Float_t fakeAntiLLooseMvaUp3_;
  Float_t fakeAntiLMediumMvaUp3_;
  Float_t fakeAntiLTightMvaUp3_;
  Float_t fakeAntiLVTightMvaUp3_;
 

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
  ntuple_->Branch("trigWeight",&trigWeight_,"trigWeight/F");
  ntuple_->Branch("trigWeight74X",&trigWeight74X_,"trigWeight74X/F");
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

  ntuple_->Branch("nuWPt",&nuWPt_,"nuWPt/F");
  ntuple_->Branch("nuWEta",&nuWEta_,"nuWEta/F");
  ntuple_->Branch("nuWPhi",&nuWPhi_,"nuWPhi/F");

  ntuple_->Branch("fakeAntiLLoose", &fakeAntiLLoose_, "fakeAntiLLoose/F");
  ntuple_->Branch("fakeAntiLMedium",&fakeAntiLMedium_,"fakeAntiLMedium/F");
  ntuple_->Branch("fakeAntiLTight", &fakeAntiLTight_, "fakeAntiLTight/F");

  ntuple_->Branch("fakeAntiLLooseUp1", &fakeAntiLLooseUp1_, "fakeAntiLLooseUp1/F");
  ntuple_->Branch("fakeAntiLMediumUp1",&fakeAntiLMediumUp1_,"fakeAntiLMediumUp1/F");
  ntuple_->Branch("fakeAntiLTightUp1", &fakeAntiLTightUp1_, "fakeAntiLTightUp1/F");

  ntuple_->Branch("fakeAntiLLooseUp2", &fakeAntiLLooseUp2_, "fakeAntiLLooseUp2/F");
  ntuple_->Branch("fakeAntiLMediumUp2",&fakeAntiLMediumUp2_,"fakeAntiLMediumUp2/F");
  ntuple_->Branch("fakeAntiLTightUp2" ,&fakeAntiLTightUp2_, "fakeAntiLTightUp2/F");

  ntuple_->Branch("fakeAntiLLooseUp3", &fakeAntiLLooseUp3_, "fakeAntiLLooseUp3/F");
  ntuple_->Branch("fakeAntiLMediumUp3",&fakeAntiLMediumUp3_,"fakeAntiLMediumUp3/F");
  ntuple_->Branch("fakeAntiLTightUp3", &fakeAntiLTightUp3_, "fakeAntiLTightUp3/F");

  ntuple_->Branch("fakeAntiLLooseMva", &fakeAntiLLooseMva_, "fakeAntiLLooseMva/F");
  ntuple_->Branch("fakeAntiLMediumMva",&fakeAntiLMediumMva_,"fakeAntiLMediumMva/F");
  ntuple_->Branch("fakeAntiLTightMva", &fakeAntiLTightMva_, "fakeAntiLTightMva/F");
  ntuple_->Branch("fakeAntiLVTightMva", &fakeAntiLVTightMva_, "fakeAntiLVTightMva/F");

  ntuple_->Branch("fakeAntiLLooseMvaUp1", &fakeAntiLLooseMvaUp1_, "fakeAntiLLooseMvaUp1/F");
  ntuple_->Branch("fakeAntiLMediumMvaUp1",&fakeAntiLMediumMvaUp1_,"fakeAntiLMediumMvaUp1/F");
  ntuple_->Branch("fakeAntiLTightMvaUp1", &fakeAntiLTightMvaUp1_, "fakeAntiLTightMvaUp1/F");
  ntuple_->Branch("fakeAntiLVTightMvaUp1", &fakeAntiLVTightMvaUp1_, "fakeAntiLVTightMvaUp1/F");

  ntuple_->Branch("fakeAntiLLooseMvaUp2", &fakeAntiLLooseMvaUp2_, "fakeAntiLLooseMvaUp2/F");
  ntuple_->Branch("fakeAntiLMediumMvaUp2",&fakeAntiLMediumMvaUp2_,"fakeAntiLMediumMvaUp2/F");
  ntuple_->Branch("fakeAntiLTightMvaUp2" ,&fakeAntiLTightMvaUp2_, "fakeAntiLTightMvaUp2/F");
  ntuple_->Branch("fakeAntiLVTightMvaUp2" ,&fakeAntiLVTightMvaUp2_, "fakeAntiLVTightMvaUp2/F");

  ntuple_->Branch("fakeAntiLLooseMvaUp3", &fakeAntiLLooseMvaUp3_, "fakeAntiLLooseMvaUp3/F");
  ntuple_->Branch("fakeAntiLMediumMvaUp3",&fakeAntiLMediumMvaUp3_,"fakeAntiLMediumMvaUp3/F");
  ntuple_->Branch("fakeAntiLTightMvaUp3", &fakeAntiLTightMvaUp3_, "fakeAntiLTightMvaUp3/F");
  ntuple_->Branch("fakeAntiLVTightMvaUp3", &fakeAntiLVTightMvaUp3_, "fakeAntiLVTightMvaUp3/F");

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

  ntuple_->Branch("tauPt",  &tauPt_,  "tauPt/F");
  ntuple_->Branch("tauEta", &tauEta_, "tauEta/F");
  ntuple_->Branch("tauPhi", &tauPhi_, "tauPhi/F");
  ntuple_->Branch("tauMass",&tauMass_,"tauMass/F");
  ntuple_->Branch("tauQ",   &tauQ_,   "tauQ/I");

  ntuple_->Branch("tauLeadingTrackPt",&tauLeadingTrackPt_,"tauLeadingTrackPt/F");
  ntuple_->Branch("tauLeadingTrackEta",&tauLeadingTrackEta_,"tauLeadingTrackEta/F");
  ntuple_->Branch("tauLeadingTrackPhi",&tauLeadingTrackPhi_,"tauLeadingTrackPhi/F");
  ntuple_->Branch("tauLeadingTrackDz",&tauLeadingTrackDz_,"tauLeadingTrackDz/F");
  ntuple_->Branch("tauLeadingTrackDxy",&tauLeadingTrackDxy_,"tauLeadingTrackDxy/F");

  ntuple_->Branch("recoilRatio",&recoilRatio_,"recoilRatio/F");
  ntuple_->Branch("recoilDPhi",&recoilDPhi_,"recoilDPhi/F");

  ntuple_->Branch("recoilM",&recoilM_,"recoilM/F");
  ntuple_->Branch("recoilPt",&recoilPt_,"recoilPt/F");
  ntuple_->Branch("recoilEta",&recoilEta_,"recoilEta/F");
  ntuple_->Branch("recoilPhi",&recoilPhi_,"recoilPhi/F");

  ntuple_->Branch("tauDecay",   &tauDecay_,   "tauDecay/I");
  ntuple_->Branch("tauGenDecay",&tauGenDecay_,"tauGenDecay/I");
  ntuple_->Branch("tauGenMatchDecay",&tauGenMatchDecay_,"tauGenMatchDecay/I");

  ntuple_->Branch("tauNtrk1", &tauNtrk1_, "tauNtrk1/i");
  ntuple_->Branch("tauNtrk08",&tauNtrk08_,"tauNtrk08/i");
  ntuple_->Branch("tauNtrk05",&tauNtrk05_,"tauNtrk05/i");

  ntuple_->Branch("tauDM",&tauDM_,"tauDM/O");
  ntuple_->Branch("tauNewDM",&tauNewDM_,"tauNewDM/O");

  ntuple_->Branch("tauLooseIso", &tauLooseIso_, "tauLooseIso/O");
  ntuple_->Branch("tauMediumIso",&tauMediumIso_,"tauMediumIso/O");
  ntuple_->Branch("tauTightIso", &tauTightIso_, "tauTightIso/O");

  ntuple_->Branch("tauLooseMvaIso", &tauLooseMvaIso_, "tauLooseMvaIso/O");
  ntuple_->Branch("tauMediumMvaIso",&tauMediumMvaIso_,"tauMediumMvaIso/O");
  ntuple_->Branch("tauTightMvaIso", &tauTightMvaIso_, "tauTightMvaIso/O");
  ntuple_->Branch("tauVTightMvaIso", &tauVTightMvaIso_, "tauVTightMvaIso/O");

  ntuple_->Branch("tauAntiMuonLoose3",&tauAntiMuonLoose3_,"tauAntiMuonLoose3/O");
  ntuple_->Branch("tauAntiMuonTight3",&tauAntiMuonTight3_,"tauAntiMuonTight3/O");

  ntuple_->Branch("tauAntiElectronVLooseMVA5",&tauAntiElectronVLooseMVA5_,"tauAntiElectronVLooseMVA5/O");
  ntuple_->Branch("tauAntiElectronLooseMVA5", &tauAntiElectronLooseMVA5_, "tauAntiElectronLooseMVA5/O");

  ntuple_->Branch("tauAntiElectronVLooseMVA6",&tauAntiElectronVLooseMVA6_,"tauAntiElectronVLooseMVA6/O");
  ntuple_->Branch("tauAntiElectronLooseMVA6", &tauAntiElectronLooseMVA6_, "tauAntiElectronLooseMVA6/O");

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

  Bool_t trigger_;
  Bool_t isWTrig_;
  Bool_t isZTrig_;
  Float_t metNoMu_;
  Float_t mhtNoMu_;
  Float_t metNoSelMu_;
  Float_t mhtNoSelMu_;
  UInt_t nMuonTrig_;
  UInt_t nSelMuonTrig_;

  TTree * trigNTuple_ = new TTree("TriggerNTuple","TriggerNTuple");
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

 // Trigger efficiencies
  TFile * trigEffFile = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+trigEffFileName);
  TF1 * trigEffDataLowerMt = (TF1*)trigEffFile->Get("MhtLt100_data");
  TF1 * trigEffDataUpperMt = (TF1*)trigEffFile->Get("MhtGt100_data");
  TF1 * trigEffMCLowerMt   = (TF1*)trigEffFile->Get("MhtLt100_mc");
  TF1 * trigEffMCUpperMt   = (TF1*)trigEffFile->Get("MhtGt100_mc");

  TFile * trigEffFile_74X = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+trigEffFileName74X);
  TF1 * trigEffDataLowerMt_74X = (TF1*)trigEffFile_74X->Get("MhtLt100_data");
  TF1 * trigEffDataUpperMt_74X = (TF1*)trigEffFile_74X->Get("MhtGt100_data");
  TF1 * trigEffMCLowerMt_74X   = (TF1*)trigEffFile_74X->Get("MhtLt100_mc");
  TF1 * trigEffMCUpperMt_74X   = (TF1*)trigEffFile_74X->Get("MhtGt100_mc");

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
      trigWeight_ = 1;
      trigWeight74X_ = 1;
      puWeight_ = 1;

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

      nuWPt_ = -1;
      nuWEta_ = 0;
      nuWPhi_ = 0;

      //      fakeLoose_ = 1.;
      //      fakeMedium_ = 1.;
      //      fakeTight_ = 1.;

      fakeAntiLLooseUp1_ = 1.;
      fakeAntiLMediumUp1_ = 1.;
      fakeAntiLTightUp1_ = 1.;

      fakeAntiLLooseUp2_ = 1.;
      fakeAntiLMediumUp2_ = 1.;
      fakeAntiLTightUp2_ = 1.;

      fakeAntiLLooseUp3_ = 1.;
      fakeAntiLMediumUp3_ = 1.;
      fakeAntiLTightUp3_ = 1.;

      fakeAntiLLoose_ = 1.;
      fakeAntiLMedium_ = 1.;
      fakeAntiLTight_ = 1.;

      fakeAntiLLooseMvaUp1_ = 1.;
      fakeAntiLMediumMvaUp1_ = 1.;
      fakeAntiLTightMvaUp1_ = 1.;
      fakeAntiLVTightMvaUp1_ = 1.;

      fakeAntiLLooseMvaUp2_ = 1.;
      fakeAntiLMediumMvaUp2_ = 1.;
      fakeAntiLTightMvaUp2_ = 1.;
      fakeAntiLVTightMvaUp2_ = 1.;

      fakeAntiLLooseMvaUp3_ = 1.;
      fakeAntiLMediumMvaUp3_ = 1.;
      fakeAntiLTightMvaUp3_ = 1.;
      fakeAntiLVTightMvaUp3_ = 1.;

      fakeAntiLLooseMva_ = 1.;
      fakeAntiLMediumMva_ = 1.;
      fakeAntiLTightMva_ = 1.;
      fakeAntiLVTightMva_ = 1.;

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

      trigger_ = false;
      isWTrig_ = false;
      isZTrig_ = false;
      metNoMu_ = 0;
      mhtNoMu_ = 0;
      metNoSelMu_ = 0;
      mhtNoSelMu_ = 0;
      nMuonTrig_ = 0;
      nSelMuonTrig_ = 0;

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
      
      if (applyHT100Cut && !isData && analysisTree.genparticles_lheHt>100) continue; 

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
      //      bool isSingleMuonHLT = false;
      bool isMetHLT = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	//	if (trigName.Contains(SingleMuonHLTName)) {
	  //	  std::cout << it->first << " : " << it->second << std::endl;
	//	  if (it->second==1)
	//	    isSingleMuonHLT = true;
	//	}
	if (trigName.Contains(MetHLTName)) {
	  //	  std::cout << it->first << " : " << it->second << std::endl;
          if (it->second==1)
            isMetHLT = true;
        }
      }
      trigger_ = isMetHLT;
      trig_ = isMetHLT;

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
      if (isData) {
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
	  if (!isData) trigMatch = true;
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
	isZTrig_ = lorentzVectorZ.M()  > ZMassCut_Trig;
	isWTrig_ = mtmuon_ >  mtCut_Trig;
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
      mhtNoMu_    = (lorentzVectorAllJetsForMht - lorentzVectorAllMuons).Pt();
      mhtNoSelMu_ = (lorentzVectorAllJetsForMht - lorentzVectorAllSelMuons).Pt();
      TLorentzVector lorentzVectorJet; lorentzVectorJet.SetXYZT(0,0,0,0);
      TLorentzVector lorentzVectorJet2; lorentzVectorJet2.SetXYZT(0,0,0,0);
      if (nJetsCentral30_>0) {
	unsigned int indexJet0 = centralJets30Indexes.at(0);
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
      }
      if (nJetsCentral30_>1) {
	unsigned int indexJet1 = centralJets30Indexes.at(1);
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

      // ************************
      // **** accessing taus ****
      // ************************
      std::vector<unsigned int> tauIndexes; tauIndexes.clear();
      std::vector<unsigned int> tau20Indexes; tau20Indexes.clear();
      std::vector<unsigned int> tau30Indexes; tau30Indexes.clear();
      std::vector<int> tauGenMatchDecay; tauGenMatchDecay.clear();
      for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) { // loop over taus

	analysisTree.tau_px[itau]   *= tauMomScale;
	analysisTree.tau_py[itau]   *= tauMomScale;
	analysisTree.tau_pz[itau]   *= tauMomScale;
	analysisTree.tau_pt[itau]   *= tauMomScale;
	analysisTree.tau_e[itau]    *= tauMomScale;
	analysisTree.tau_mass[itau] *= tauMomScale;

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
      bool isSingleJet = nSelTaus_==1 && nJetsCentral30_<=1 && nJetsForward30_==0;
      bool isNoJets    = nSelTaus_==0 && nJetsCentral30_==0 && nJetsForward30_==0;
      TLorentzVector lorentzVectorTau; lorentzVectorTau.SetXYZT(0,0,0,0);
      if (nSelTaus_>0) {

	unsigned int indexTau = tauIndexes.at(0);
	lorentzVectorTau.SetXYZM(analysisTree.tau_px[indexTau],
				 analysisTree.tau_py[indexTau],
				 analysisTree.tau_pz[indexTau],
				 analysisTree.tau_mass[indexTau]);

	mttau_ = mT(lorentzVectorTau,lorentzVectorMet);
	tauPt_ = analysisTree.tau_pt[indexTau];
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

	tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
	tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;

	tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;

	tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;

	//	std::cout << "TauMva : Loose = " << tauLooseMvaIso_
	//		  << "  Medium = " << tauMediumMvaIso_
	//		  << "  Tight = " << tauTightMvaIso_ 
	//		  << "  VTight = " << tauVTightMvaIso_ << std::endl;

	tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
	tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;

	tauAntiElectronVLooseMVA5_ = analysisTree.tau_againstElectronVLooseMVA5[indexTau] > 0.5;
	tauAntiElectronLooseMVA5_ = analysisTree.tau_againstElectronLooseMVA5[indexTau] > 0.5;

	tauAntiElectronVLooseMVA6_ = analysisTree.tau_againstElectronVLooseMVA6[indexTau] > 0.5;
	tauAntiElectronLooseMVA6_ = analysisTree.tau_againstElectronLooseMVA6[indexTau] > 0.5;

	float fake[4]; 
	float efake[4];
	int iso;

	FakeRate(tauPt_,fake, efake,iso);

	fakeAntiLLoose_  = fake[0]; 
	fakeAntiLMedium_ = fake[1];
	fakeAntiLTight_  = fake[2];
	
	fakeAntiLLooseUp1_  = fake[0]; 
	fakeAntiLMediumUp1_ = fake[1];
	fakeAntiLTightUp1_  = fake[2];

	fakeAntiLLooseUp2_  = fake[0];
	fakeAntiLMediumUp2_ = fake[1];
	fakeAntiLTightUp2_  = fake[2];

	fakeAntiLLooseUp3_  = fake[0];
	fakeAntiLMediumUp3_ = fake[1];
	fakeAntiLTightUp3_  = fake[2];

	if (iso==0) {
	  fakeAntiLLooseUp1_  = fake[0] + efake[0]; 
	  fakeAntiLMediumUp1_ = fake[1] + efake[1];
	  fakeAntiLTightUp1_  = fake[2] + efake[2];
	}
	else if (iso==1) {
	  fakeAntiLLooseUp2_  = fake[0] + efake[0]; 
	  fakeAntiLMediumUp2_ = fake[1] + efake[1];
	  fakeAntiLTightUp2_  = fake[2] + efake[2];
	}
	else {
	  fakeAntiLLooseUp3_  = fake[0] + efake[0]; 
	  fakeAntiLMediumUp3_ = fake[1] + efake[1];
	  fakeAntiLTightUp3_  = fake[2] + efake[2];
	}

	// MVA
	FakeRateMva(tauPt_,fake, efake,iso);

	fakeAntiLLooseMva_   = fake[0]; 
	fakeAntiLMediumMva_  = fake[1];
	fakeAntiLTightMva_   = fake[2];
	fakeAntiLVTightMva_  = fake[3];
	
	fakeAntiLLooseMvaUp1_   = fake[0]; 
	fakeAntiLMediumMvaUp1_  = fake[1];
	fakeAntiLTightMvaUp1_   = fake[2];
	fakeAntiLVTightMvaUp1_  = fake[3];

	fakeAntiLLooseMvaUp2_   = fake[0];
	fakeAntiLMediumMvaUp2_  = fake[1];
	fakeAntiLTightMvaUp2_   = fake[2];
	fakeAntiLVTightMvaUp2_  = fake[3];

	fakeAntiLLooseMvaUp3_   = fake[0];
	fakeAntiLMediumMvaUp3_  = fake[1];
	fakeAntiLTightMvaUp3_   = fake[2];
	fakeAntiLVTightMvaUp3_  = fake[3];

	if (iso==0) {
	  fakeAntiLLooseMvaUp1_  = fake[0] + efake[0]; 
	  fakeAntiLMediumMvaUp1_ = fake[1] + efake[1];
	  fakeAntiLTightMvaUp1_  = fake[2] + efake[2];
	  fakeAntiLVTightMvaUp1_ = fake[3] + efake[3];
	}
	else if (iso==1) {
	  fakeAntiLLooseMvaUp2_  = fake[0] + efake[0]; 
	  fakeAntiLMediumMvaUp2_ = fake[1] + efake[1];
	  fakeAntiLTightMvaUp2_  = fake[2] + efake[2];
	  fakeAntiLVTightMvaUp2_ = fake[3] + efake[3];
	}
	else {
	  fakeAntiLLooseMvaUp3_  = fake[0] + efake[0]; 
	  fakeAntiLMediumMvaUp3_ = fake[1] + efake[1];
	  fakeAntiLTightMvaUp3_  = fake[2] + efake[2];
	  fakeAntiLVTightMvaUp3_ = fake[3] + efake[3];
	}
	

	//	cout << "fake  Loose = " << fakeAntiLLoose_
	//	     << "   Medium = " << fakeAntiLMedium_
	//	     << "   Tight  = " << fakeAntiLTight_ << endl;

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
      float trigEffData_74X = 1.0;
      float trigEffMC_74X   = 1.0;

      if (metNoMu_>100&&metNoMu_<300) {
	if (mhtNoMu_<100) {
	  trigEffData = trigEffDataLowerMt->Eval(metNoMu_);
	  trigEffMC   = trigEffMCLowerMt->Eval(metNoMu_);
	  trigEffData_74X = trigEffDataLowerMt_74X->Eval(metNoMu_);
	  trigEffMC_74X   = trigEffMCLowerMt_74X->Eval(metNoMu_);
	}
	else {
	  trigEffData = trigEffDataUpperMt->Eval(metNoMu_);
          trigEffMC   = trigEffMCUpperMt->Eval(metNoMu_);
	  trigEffData_74X = trigEffDataUpperMt_74X->Eval(metNoMu_);
          trigEffMC_74X   = trigEffMCUpperMt_74X->Eval(metNoMu_);
	}
	
      }
      if (trigEffData<0) trigEffData = 0;
      if (trigEffData_74X<0) trigEffData_74X = 0;
      trigWeight_ = trigEffData;
      trigWeight74X_ = trigEffData_74X;
      if (debug) {
	cout << "MetNoMu = " << metNoMu_ 
	     << "  MhtNoMu = " << mhtNoMu_ 
	     << "  trigWeight = " << trigWeight_ << endl;
      }
      weight_ *= trigWeight_;

      // ********************************
      // **** filling trigger ntuple ****
      // ********************************
      if (ptTriggerMu>ptTrigMuCut) {
	trigNTuple_->Fill();
	TrigEvents++;
      }
      

      // ******************************
      // ********* ZJet selection *****
      // ******************************
      /*
      if (lorentzVectorZ.Pt()>1e-4) {
	recoilRatio_ = tauPt_ / lorentzVectorZ.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorZ,lorentzVectorTau);
	isZJet = ptLeadingMu>ptLeadingMuCut_ZJet;
	isZJet = isZJet && ptTrailingMu>ptTrailingMuCut_ZJet;
	isZJet = isZJet && lorentzVectorZ.M()>ZMassLowerCut_ZJet && lorentzVectorZ.M()<ZMassUpperCut_ZJet;
	isZJet = isZJet && recoilRatio_>ptJetZRatioLowerCut_ZJet && recoilRatio_<ptJetZRatioUpperCut_ZJet;
	isZJet = isZJet && recoilDPhi_>deltaPhiZJetCut_ZJet;
	if (isZJet) { 
	  float leadMuIso = SF_muonIdIso->get_ScaleFactor(ptLeadingMu, etaLeadingMu);
	  float trailMuIso = SF_muonIdIso->get_ScaleFactor(ptTrailingMu, etaTrailingMu);
	  mueffweight = leadMuIso*trailMuIso;
	  //	  float leadMuTrig = SF_muonTrig->get_EfficiencyData(ptLeadingMu, etaLeadingMu);
	  //	  float trailMuTrig = SF_muonTrig->get_EfficiencyData(ptTrailingMu, etaTrailingMu);
	  //	  mutrigweight = 1 - (1-leadMuTrig)*(1-trailMuTrig);
	  mutrigweight = SF_muonTrig->get_EfficiencyData(ptTriggerMu, etaTriggerMu);
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
      */
      // ***************************
      // ******* WJet selection ****
      // ***************************
      if (lorentzVectorW.Pt()>1e-4) {
	recoilRatio_ = tauPt_ / lorentzVectorW.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorTau);
	isWJet = ptTriggerMu>ptMuCut_WJet; 
	isWJet = isWJet && mtmuon_ > mtCut_WJet;
	isWJet = isWJet && recoilRatio_>ptJetWRatioLowerCut_WJet && recoilRatio_<ptJetWRatioUpperCut_WJet;
	isWJet = isWJet && recoilDPhi_>deltaPhiWJetCut_WJet;
	if (isWJet) {
	  mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
          mutrigweight = SF_muonTrig->get_EfficiencyData(ptTriggerMu, etaTriggerMu);
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
      
      // ***************************
      // ******* W+Jet selection ****
      // ***************************
      /*
      if (lorentzVectorW.Pt()>1e-4&&nJetsCentral30_>=1) {
	recoilRatio_ = jetPt_ / lorentzVectorW.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorJet);
	isWJet = ptTriggerMu>ptMuCut_WJet; 
	isWJet = isWJet && mtmuon_ > mtCut_WJet;
	isWJet = isWJet && recoilRatio_>ptJetWRatioLowerCut_WJet && recoilRatio_<ptJetWRatioUpperCut_WJet;
	isWJet = isWJet && recoilDPhi_>deltaPhiWJetCut_WJet;
	if (isWJet) {
	  mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
          mutrigweight = SF_muonTrig->get_EfficiencyData(ptTriggerMu, etaTriggerMu);
	  HtNoRecoil_     = Ht_     - ptTriggerMu;
	  SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu;
	  recoilM_   = lorentzVectorW.M();
          recoilPt_  = lorentzVectorW.Pt();
          recoilEta_ = lorentzVectorW.Eta();
          recoilPhi_ = lorentzVectorW.Phi();
	  selection_ = 13;
	  ntuple_->Fill();
	  WJetEvents++;
	}
      }
      */
      // ********************************
      // ****** W->mu+v selection *******
      // ********************************
      /*
      if (lorentzVectorW.Pt()>1e-4) {
	bool isWprod = ptTriggerMu>ptMuCut_WJet;
        isWprod = isWprod && mtmuon_ > mtCut_WJet;
	if (isWprod) {
	  mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
          mutrigweight = SF_muonTrig->get_EfficiencyData(ptTriggerMu, etaTriggerMu);
	  recoilRatio_ = lorentzVectorTriggerMu.Pt() / lorentzVectorMet.Pt();
	  recoilDPhi_  = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet);
	  HtNoRecoil_     = Ht_     - ptTriggerMu;
	  SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu;
	  recoilM_   = lorentzVectorW.M();
          recoilPt_  = lorentzVectorW.Pt();
          recoilEta_ = lorentzVectorW.Eta();
          recoilPhi_ = lorentzVectorW.Phi();
	  selection_ = 10;
	  ntuple_->Fill();
	  WProdEvents++;
	}
      }
      */
      // ********************************
      // ******* W*->MuNu selection *****
      // ********************************
      if (lorentzVectorMet.Pt()>1e-4) {
	recoilRatio_ = ptTriggerMu/lorentzVectorMet.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet);
	isWMuNu = ptTriggerMu>ptMuCut_WMuNu;
	isWMuNu = isWMuNu && met_>metCut_WMuNu;
	isWMuNu = isWMuNu && recoilRatio_>ptMuMetRatioLowerCut_WMuNu && recoilRatio_<ptMuMetRatioUpperCut_WMuNu;
	isWMuNu = isWMuNu && recoilDPhi_>deltaPhiMuMetCut_WMuNu;
	if (isWMuNu) {
	  mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
          mutrigweight = SF_muonTrig->get_EfficiencyData(ptTriggerMu, etaTriggerMu);
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

      // *********************************
      // ** Single JET + MET selection ***
      // *********************************
      /*
      if (lorentzVectorMet.Pt()>1e-4&&nJetsCentral30_==1) {
	recoilRatio_ = jetPt_/lorentzVectorMet.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorJet,lorentzVectorMet);
	if (recoilRatio_>0.7&&recoilRatio_<1.3&&recoilDPhi_>2.4) {
	  HtNoRecoil_     = Ht_;
	  SoftHtNoRecoil_ = SoftHt_;
	  recoilM_   = lorentzVectorMet.M();
	  recoilPt_  = lorentzVectorMet.Pt();
	  recoilEta_ = lorentzVectorMet.Eta();
	  recoilPhi_ = lorentzVectorMet.Phi();
	  selection_ = 11;
	  ntuple_->Fill();
	  SingleJetEvents++;
	}

      }
      */
      // *************************
      // **** Dijet selection ****
      // *************************
      /*
      if (nJetsCentral30_>=2) {
	recoilRatio_ = jet2Pt_/jetPt_;
	recoilDPhi_  = dPhiFromLV(lorentzVectorJet,lorentzVectorJet2);	
	if (recoilRatio_>0.7&&recoilRatio_<1.3&&recoilDPhi_>2.4) {
	  HtNoRecoil_     = Ht_;
	  SoftHtNoRecoil_ = SoftHt_;
	  recoilM_   = lorentzVectorMet.M();
	  recoilPt_  = lorentzVectorMet.Pt();
	  recoilEta_ = lorentzVectorMet.Eta();
	  recoilPhi_ = lorentzVectorMet.Phi();
	  selection_ = 12;
	  ntuple_->Fill();
	  DiJetEvents++;
	}
      }
      */

      // 

      // ********************************
      // ****** W*->TauNu selection *****
      // ******************************** 
      if (lorentzVectorMet.Pt()>1e-4) {
	recoilRatio_ = tauPt_ / lorentzVectorMet.Pt();
	recoilDPhi_  = dPhiFromLV(lorentzVectorTau,lorentzVectorMet);
	isWTauNu = met_>metCut_WTauNu;
	isWTauNu = isWTauNu && recoilRatio_>ptTauMetRatioLowerCut_WTauNu && recoilRatio_<ptTauMetRatioUpperCut_WTauNu;
	isWTauNu = isWTauNu && recoilDPhi_>deltaPhiTauMetCut_WTauNu;
	if (isWTauNu) {
	  HtNoRecoil_     = Ht_;
	  SoftHtNoRecoil_ = SoftHt_;
	  recoilM_   = lorentzVectorMet.M();
          recoilPt_  = lorentzVectorMet.Pt();
          recoilEta_ = lorentzVectorMet.Eta();
          recoilPhi_ = lorentzVectorMet.Phi();
	  metFilters_ = metFiltersPasses(analysisTree,metFlags);
	  selection_ = 3;
	  ntuple_->Fill();
	  WTauNuEvents++;
	  // filling special ntuple ->
	  // bool fillWNTuple = nMuon_==0 && nElec_==0 && isSingleJet && wMass_>0;
	  // fillWNTuple = fillWNTuple && tauDM_>0.5 && tauLooseIso_>0.5 && tauAntiMuonLoose3_>0.5 && tauAntiElectronLooseMVA6_>0.5;
	  // fillWNTuple = fillWNTuple && tauGenMatchDecay_<0;
	  // if (fillWNTuple) wntuple_->Fill();
	  // end filling special ntuple 

	  // bool isInterestingEvent =  
	  //   nMuon_==0 && 
	  //   nElec_==0 && 
	  //   nSelTaus_==1 && 
	  //   tauDM_>0.5 && 
	  //   tauLooseIso_>0.5 && 
	  //   tauAntiMuonLoose3_>0.5 && 
	  //   tauAntiElectronLooseMVA5_>0.5 && 
	  //   tauPt_>100 && met_>110 && 
	  //   trigger_>0.5 && 
	  //   nJetsCentral30_ <= 1 && 
	  //   nJetsForward30_ == 0;
	  // if (isInterestingEvent) {
	  //   ntuple_->Fill();
	  // std::cout << "Selected taus = " << nSelTaus_
	  // 	      << "  tauPt = " << tauPt_
	  // 	      << "  tauEta = " << tauEta_
	  // 	      << "  Forward jets (30) = " << nJetsForward30_
	  // 	      << "  Central jets (30) = " << nJetsCentral30_ << std::endl;
	  // for (unsigned int iF=0; iF<centralJets30Indexes.size();++iF) {
	  //   unsigned int indexJet = centralJets30Indexes.at(iF);
	  //   float deltaRJetTau = deltaR(analysisTree.pfjet_eta[indexJet],
	  // 				  analysisTree.pfjet_phi[indexJet],
	  // 				  tauEta_,
	  // 				  tauPhi_);
	  
	  //   std::cout << "  Central jet " << iF 
	  // 		<< "   pT = " << analysisTree.pfjet_pt[indexJet] 
	  // 		<< "   eta = " << analysisTree.pfjet_eta[indexJet] 
	  // 		<< " dR(jet,tau) = " << deltaRJetTau << std::endl;
	  
	  
	  // }
	  // for (unsigned int iF=0; iF<forwardJets30Indexes.size();++iF) {
	  //   unsigned int indexJet = forwardJets30Indexes.at(iF);
	  //   float deltaRJetTau = deltaR(analysisTree.pfjet_eta[indexJet],
	  // 				  analysisTree.pfjet_phi[indexJet],
	  // 				  tauEta_,
	  // 				  tauPhi_);
	  
	  //   std::cout << "  Forward jet " << iF 
	  // 		<< "   pT = " << analysisTree.pfjet_pt[indexJet] 
	  // 		<< "   eta = " << analysisTree.pfjet_eta[indexJet] 
	  // 		<< " dR(jet,tau) = " << deltaRJetTau << std::endl;
	  
	  
	  // }
	  // std::cout << std::endl;
	  //	  }
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
	    if (ptTauJetRatio<ptTauJetRatioLowerCut_DiJet) continue;
	    if (ptTauJetRatio>ptTauJetRatioUpperCut_DiJet) continue;
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
	    tauPt_ = analysisTree.tau_pt[indexTau];
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
	    tauGenMatchDecay_ = tauGenMatchDecay.at(iTau);

	    if (tauDecay_<0) tauDecay_ = -1;
	    if (tauGenDecay_<0) tauGenDecay_ = -1;
	    if (tauGenMatchDecay_<0) tauGenMatchDecay_ = -1;

	    pfJet40_ = jets40trigger.at(acceptedJetDirIndex);
	    pfJet60_ = jets60trigger.at(acceptedJetDirIndex);
	    pfJet80_ = jets80trigger.at(acceptedJetDirIndex);
	    pfJet140_ = jets140trigger.at(acceptedJetDirIndex);

	    tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
	    tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;

	    tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	    tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
	    tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;

	    tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
	    tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;

	    tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
	    tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;
	    tauAntiElectronVLooseMVA5_ = analysisTree.tau_againstElectronVLooseMVA5[indexTau] > 0.5;
	    tauAntiElectronLooseMVA5_ = analysisTree.tau_againstElectronLooseMVA5[indexTau] > 0.5;

	    recoilRatio_ = tauPt_/recoilJetLV.Pt();
	    recoilDPhi_ = dPhiFromLV(tauLV,recoilJetLV);
	    recoilM_ = recoilJetLV.M();
	    recoilPt_ = recoilJetLV.Pt();
	    recoilEta_ = recoilJetLV.Eta();
	    recoilPhi_ = recoilJetLV.Phi();
	    HtNoRecoil_     = Ht_     - recoilJetLV.Pt();
	    SoftHtNoRecoil_ = SoftHt_ - recoilJetLV.Pt();
	    selection_ = 4;
	    ntuple_->Fill();
	  }
	}
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
  std::cout << "Total number of Z+Jet events      = " << ZJetEvents << std::endl;
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
  
}



