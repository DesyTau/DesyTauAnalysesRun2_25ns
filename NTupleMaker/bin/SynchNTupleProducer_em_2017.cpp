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
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

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

float totalTransverseMass(TLorentzVector l1,
                          TLorentzVector l2,
                          TLorentzVector l3) {
   
   TLorentzVector totalLV = l1 + l2 + l3;
   float    totalET = l1.Pt() +  l2.Pt() + l3.Pt();
   float     mTtot = TMath::Sqrt(totalET*totalET-totalLV.Pt()*totalLV.Pt());
   return    mTtot;
   
}

void computeDzeta(float metX,  float metY,           //for signal ~0, for bg large negativ values because MET mot aligned with tau decay products
                  float zetaX, float zetaY,
                  float pzetavis,
                  float & pzetamiss,
                  float & dzeta) {
   
   pzetamiss = metX*zetaX + metY*zetaY;
   dzeta = pzetamiss - 0.85*pzetavis;
}


float topPtWeight(float pt1,                         //top pT re-weighting either Run1 oder Run2 values
                  float pt2,
		  bool run1) {
    
  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter

  if (run1) {
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
    a = 0.156;    // Run1 a parameter
    b = -0.00137;  // Run1 b parameter
  }
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
  
  return TMath::Sqrt(w1*w2);
    
}

ClassicSVfit SVFitMassComputation(classic_svFit::MeasuredTauLepton svFitEle,
                                              classic_svFit::MeasuredTauLepton svFitMu,
                                              double measuredMVAMETx,
                                              double measuredMVAMETy,
                                              TMatrixD covMVAMET,
                                              TFile * inputFile_visPtResolution
                                              ) {
    
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitEle);
    measuredTauLeptons.push_back(svFitMu);
    
    int verbosity = 1;
    ClassicSVfit svFitAlgo(verbosity);
    double kappa = 3.; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
    svFitAlgo.addLogM_fixed(true, kappa);
    svFitAlgo.integrate(measuredTauLeptons, measuredMVAMETx, measuredMVAMETy, covMVAMET);
    
    return svFitAlgo;
    
}

bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags) {  // metFlags: vector with event filters to be checked, 
                                                                     // compared with filter results stored in n-tuple
   bool passed = true;                                               // returns true if all filters are passed
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

//struct myclass {
//  bool operator() (int i,int j) { return (i<j);}
//} myobject, myobjectX;

//const float electronMass = 0;
//const float muonMass = 0.10565837;
//const float pionMass = 0.1396;

// float getEffectiveArea(float eta) {  //effective areas for v2 Id
//     float effArea = 0.1440;
//     float absEta = fabs(eta);
//     if (absEta<1.0) effArea = 0.1440;
//     else if (absEta < 1.4790) effArea = 0.1562;
//     else if (absEta < 2.0) effArea = 0.1032;
//     else if (absEta < 2.2) effArea = 0.0859;
//     else if (absEta < 2.3) effArea = 0.1116;
//     else if (absEta < 2.4) effArea = 0.1321;
//     else if (absEta < 2.5) effArea = 0.1654;
//     return effArea;

//   } 



float getEffectiveArea(float eta) {
    float effArea = 0.1566;
    float absEta = fabs(eta);
    if (absEta<1.0) effArea = 0.1566;
    else if (absEta < 1.4790) effArea = 0.1626;
    else if (absEta < 2.0) effArea = 0.1073;
    else if (absEta < 2.2) effArea = 0.0854;
    else if (absEta < 2.3) effArea = 0.1051;
    else if (absEta < 2.4) effArea = 0.1204;
    else if (absEta < 5.0) effArea = 0.1524;
    return effArea;

  } 

int main(int argc, char * argv[]) {
   
   // first argument - config file
   // second argument - filelist
   
   using namespace std;
   
   // **** configuration                                                      //read from config file, all values have to be given
   bool debug = false;
   
   Config cfg(argv[1]);                                                       //why argv[1]?
    
   const bool computeSVFitMass = cfg.get<bool>("ComputeSVFitMass");
   const bool removeGammaStar = cfg.get<bool>("RemoveGammaStar");
    
   const bool isData = cfg.get<bool>("IsData");
   const bool isDY   = cfg.get<bool>("IsDY");
   const bool isSignal   = cfg.get<bool>("IsSignal");
   const bool isEmbedded = cfg.get<bool>("IsEmbedded");
   const bool isW    = cfg.get<bool>("IsW");
   const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
   const string jsonFile = cfg.get<string>("jsonFile");
   const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
   const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
   
   // kinematic cuts on electrons
   const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
   const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
   const float etaElectronCut     = cfg.get<float>("etaElectronCut");
   const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
   const float dzElectronCut      = cfg.get<float>("dzElectronCut");
   const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
   const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
   const bool applyIsoElectronId  = cfg.get<bool>("ApplyIsoElectronId");
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
   
   
   const float muonScale = cfg.get<float>("MuonScale");
   const float eleScaleBarrel = cfg.get<float>("EleScaleBarrel");
   const float eleScaleEndcap = cfg.get<float>("EleScaleEndcap");
   
   const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
   TString RecoilMvaFileName(recoilMvaFileName);
   const string recoilFileName   = cfg.get<string>("RecoilFileName");
   TString RecoilFileName(recoilFileName);
    
   const string metSysFileName   = cfg.get<string>("MetSysFileName");
   TString MetSysFileName(metSysFileName);
    
   const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
   TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
    
   const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
   TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
    
   const string pileUpDataFile = cfg.get<string>("PileUpDataFile");
   const string pileUpMCFile = cfg.get<string>("PileUpMCFile");
   const string samplenameForPUHist = cfg.get<string>("SampleNameForPUHist");
   TString PileUpDataFile(pileUpDataFile);
   TString PileUpMCFile(pileUpMCFile);

   const bool applyWSCorr = cfg.get<bool>("applyWorkspaceCorrection");
   const string correctionWSFile = cfg.get<string>("CorrectionWSFile");
 
    // **** end of configuration

   const float a_jetMu = 0.902;                                //values for jets faking electron and muon, dominant in high mass region
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
    
    TString TStrName(rootFileName);                           // name of file list 
    std::cout <<TStrName <<std::endl;
    
    // output fileName with histograms
    TFile * file = new TFile(TStrName+TString(".root"),"recreate");  //create outputfilename, name of file list +.root
    file->cd("");
    
    TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);   //define output histograms to be stored in output file
    TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
    TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
    TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
    TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
    TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
    TH1D * histWeightsSkimH = new TH1D("histWeightsSkimH","",1,-0.5,0.5);
    TH1D * histWeightsTTH = new TH1D("histWeightsTTH","",1,-0.5,0.5);
    /*    
    TTree * treeGen = new TTree("GenHiggs","GenHiggs");
    */    

    Float_t higgsMass;
    Float_t higgsPt;
    Float_t higgsEta;
    
    TTree * tree = new TTree("TauCheck","TauCheck");               // define output tree with leefs
    // Declaration of leaf types
    unsigned int minRun = 99999999;
    unsigned int maxRun = 0;
    Int_t           run;
    Int_t           lumi;
    Int_t           evt;
    Int_t           npv;
    Float_t         npu;
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

    Float_t         btag0weight;
    Float_t         btag0weight_Up;
    Float_t         btag0weight_Down;

    Float_t         btag1weight;
    Float_t         btag1weight_Up;
    Float_t         btag1weight_Down;
    
    Float_t         btag2weight;
    Float_t         btag2weight_Up;
    Float_t         btag2weight_Down;

    Float_t         qcdweight;
    Float_t         qcdweightup;
    Float_t         qcdweightdown;
    
    Float_t         qcdweight_nodzeta;
    Float_t         qcdweightup_nodzeta;
    Float_t         qcdweightdown_nodzeta;
    
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
    
    Float_t         m_vis;
    Float_t         m_vis_muUp;
    Float_t         m_vis_muDown;
    Float_t         m_vis_escaleUp;
    Float_t         m_vis_escaleDown;
    Float_t         m_vis_jesUp;
    Float_t         m_vis_jesDown;
    Float_t         m_vis_resoUp;
    Float_t         m_vis_resoDown;
    
    Float_t         mTtot;
    Float_t         mTtot_muUp;
    Float_t         mTtot_muDown;
    Float_t         mTtot_escaleUp;
    Float_t         mTtot_escaleDown;
    Float_t         mTtot_jesUp;
    Float_t         mTtot_jesDown;
    Float_t         mTtot_resoUp;
    Float_t         mTtot_resoDown;
    Float_t         mCDF;

    Float_t         mTdileptonMET;
    Float_t         mTdileptonMET_muUp;
    Float_t         mTdileptonMET_muDown;
    Float_t         mTdileptonMET_escaleUp;
    Float_t         mTdileptonMET_escaleDown;
    Float_t         mTdileptonMET_jesUp;
    Float_t         mTdileptonMET_jesDown;
    Float_t         mTdileptonMET_resoUp;
    Float_t         mTdileptonMET_resoDown;
    
    Float_t         mTemu;
    Float_t         mTemet;
    Float_t         mTmumet;
    
    Float_t         m_sv;
    Float_t         mt_sv;

    Float_t         pt_sv;
    Float_t         eta_sv;
    Float_t         phi_sv;

    Float_t         pt_sv_escaleUp;
    Float_t         eta_sv_escaleUp;
    Float_t         phi_sv_escaleUp;

    Float_t         pt_sv_escaleDown;
    Float_t         eta_sv_escaleDown;
    Float_t         phi_sv_escaleDown;
    
    Float_t         pt_sv_gen;
    Float_t         eta_sv_gen;
    Float_t         phi_sv_gen;
    
    Float_t         m_sv_escaleUp;
    Float_t         m_sv_escaleDown;
    Float_t         m_sv_muUp;
    Float_t         m_sv_muDown;
    Float_t         m_sv_jesUp;
    Float_t         m_sv_jesDown;
    Float_t         m_sv_resoUp;
    Float_t         m_sv_resoDown;
    
    Float_t         mt_sv_escaleUp;
    Float_t         mt_sv_escaleDown;
    Float_t         mt_sv_muUp;
    Float_t         mt_sv_muDown;
    Float_t         mt_sv_jesUp;
    Float_t         mt_sv_jesDown;
    Float_t         mt_sv_resoUp;
    Float_t         mt_sv_resoDown;
    
    
    Float_t         pt_1;
    Float_t         pt_1_escaleUp;
    Float_t         pt_1_escaleDown;
    
    Float_t         phi_1;
    Float_t         eta_1;
    Float_t         m_1;
    Int_t           q_1;
    Float_t         iso_1;
    Float_t         mva_1;
    Float_t         d0_1;
    Float_t         d0err_1;
    Float_t         d0sig_1;
    Float_t         d0sig_1_abs;
    Float_t         d0_diff;
    Float_t         dZ_diff;
    Float_t         dZ_1;
    Float_t         dZerr_1;
    Float_t         dZsig_1;
    Float_t         dZsig_1_abs;
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
    Float_t         d0err_2;
    Float_t         d0sig_2;
    Float_t         d0sig_2_abs;
    Float_t         dZ_2;
    Float_t         dZerr_2;
    Float_t         dZsig_2;
    Float_t         dZsig_2_abs;
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
    
    Float_t         met_uncorr;
    Float_t         metphi_uncorr;
    
    Float_t         met_resoUp;
    Float_t         metphi_resoUp;
    
    Float_t         met_resoDown;
    Float_t         metphi_resoDown;
    
    Float_t         met_jesUp;
    Float_t         metphi_jesUp;
    
    Float_t         met_jesDown;
    Float_t         metphi_jesDown;
    
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
    
    Float_t         mvamet_jesUp;
    Float_t         mvametphi_jesUp;
    
    Float_t         mvamet_jesDown;
    Float_t         mvametphi_jesDown;
    
    Float_t         genmet;
    Float_t         genmetphi;
    
    Float_t         msvmet;
    Float_t         msvmetphi;
    
    Float_t         pt_tt;
    Float_t         dr_tt;
    Float_t         dphi_tt;
    Float_t         m_leptons;
    
    Float_t         pzetavis;
    Float_t         pzetamiss;
    Float_t         dzeta;
    
    Float_t         pzetamiss_mvamet;
    Float_t         dzeta_mvamet;
    
    Float_t         pzetamiss_mvamet_uncorr;
    Float_t         dzeta_mvamet_uncorr;
    
    Float_t         pzetamiss_mvamet_jesUp;
    Float_t         dzeta_mvamet_jesUp;
    
    Float_t         pzetamiss_mvamet_jesDown;
    Float_t         dzeta_mvamet_jesDown;
    
    Float_t         pzetamiss_mvamet_resoUp;
    Float_t         dzeta_mvamet_resoUp;
    
    Float_t         pzetamiss_mvamet_resoDown;
    Float_t         dzeta_mvamet_resoDown;

    Float_t         pzetamiss_jesUp;
    Float_t         dzeta_jesUp;
    
    Float_t         pzetamiss_jesDown;
    Float_t         dzeta_jesDown;
    
    Float_t         pzetamiss_resoUp;
    Float_t         dzeta_resoUp;
    
    Float_t         pzetamiss_resoDown;
    Float_t         dzeta_resoDown;

    Float_t         dzeta_escaleUp;
    Float_t         dzeta_escaleDown;
    
    Float_t         pzetavis_escaleUp;
    Float_t         pzetavis_escaleDown;

    Float_t         pzetamiss_escaleUp;
    Float_t         pzetamiss_escaleDown;

    Float_t         pzetamiss_genmet;
    Float_t         dzeta_genmet;
    
    Float_t         mva_gf;
    
    Int_t           njets;
    Int_t           njets_jesUp;
    Int_t           njets_jesDown;

    Int_t           njetspt20;


    Float_t         jpt_1;
    Float_t         jpt_1_jesUp;
    Float_t         jpt_1_jesDown;

    Float_t         jeta_1;
    Float_t         jphi_1;
    Float_t         jptraw_1;
    Float_t         jptunc_1;
    Float_t         jmva_1;
    Float_t         jlrm_1;
    Int_t           jctm_1;
    Int_t           gen_match_1;
    
    Float_t         jpt_2;
    Float_t         jpt_2_jesUp;
    Float_t         jpt_2_jesDown;


    Float_t         jeta_2;
    Float_t         jphi_2;
    Float_t         jptraw_2;
    Float_t         jptunc_2;
    Float_t         jmva_2;
    Float_t         jlrm_2;
    Int_t           jctm_2;
    Int_t           gen_match_2;
    
    Float_t         mjj;
    Float_t         mjj_jesUp;
    Float_t         mjj_jesDown;

    Float_t         jdeta;
    Int_t           njetingap;
    
    Int_t           nbtag;
    Int_t           nbtag_noSF;
    Float_t         bpt;
    Double_t         beta;
    Float_t         bphi;
    
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
    Float_t         bosonEta;
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
    


    tree->Branch("run", &run, "run/I");                                  // I=int, F=Float, O=bool
    tree->Branch("lumi", &lumi, "lumi/I");
    tree->Branch("evt", &evt, "evt/I");
    tree->Branch("npv", &npv, "npv/I");
    tree->Branch("npu", &npu, "npu/F");
    tree->Branch("rho", &rho, "rho/F");
    
    tree->Branch("isZLL",&isZLL,"isZLL/O");
    tree->Branch("isZEE",&isZEE,"iEE/O");
    tree->Branch("isZMM",&isZMM,"isZMM/O");
    tree->Branch("isZTT",&isZTT,"isZTT/O");
    
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

    tree->Branch("btag0weight",&btag0weight,"btag0weight/F");
    tree->Branch("btag0weight_Up",&btag0weight_Up,"btag0weight_Up/F");
    tree->Branch("btag0weight_Down",&btag0weight_Down,"btag0weight_Down/F");
    tree->Branch("btag1weight",&btag1weight,"btag1weight/F");
    tree->Branch("btag1weight_Up",&btag1weight_Up,"btag1weight_Up/F");
    tree->Branch("btag1weight_Down",&btag1weight_Down,"btag1weight_Down/F");
    tree->Branch("btag2weight",&btag2weight,"btag2weight/F");
    tree->Branch("btag2weight_Up",&btag2weight_Up,"btag2weight_Up/F");
    tree->Branch("btag2weight_Down",&btag2weight_Down,"btag2weight_Down/F");

    tree->Branch("qcdweight", &qcdweight, "qcdweight/F");
    tree->Branch("qcdweightup", &qcdweightup, "qcdweightup/F");
    tree->Branch("qcdweightdown", &qcdweightdown, "qcdweightdown/F");
    
    tree->Branch("qcdweight_nodzeta", &qcdweight_nodzeta, "qcdweight_nodzeta/F");
    tree->Branch("qcdweightup_nodzeta", &qcdweightup_nodzeta, "qcdweightup_nodzeta/F");
    tree->Branch("qcdweightdown_nodzeta", &qcdweightdown_nodzeta, "qcdweightdown_nodzeta/F");
    
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
    tree->Branch("muonBadTrackFilter",&muonBadTrackFilter_,"muonBadTrackFilter/O");
    tree->Branch("chargedHadronTrackResolutionFilter",&chargedHadronTrackResolutionFilter_,"chargedHadronTrackResolutionFilter/O");

    tree->Branch("badMuonFilter",&badMuonFilter_,"badMuonFilter/O");
    tree->Branch("duplicateMuonFilter",&duplicateMuonFilter_,"duplicateMuonFilter/O");
    
    tree->Branch("m_vis",        &m_vis,        "m_vis/F");
    tree->Branch("m_vis_muUp",   &m_vis_muUp,   "m_vis_muUp/F");
    tree->Branch("m_vis_muDown", &m_vis_muDown, "m_vis_muDown/F");
    tree->Branch("m_vis_escaleUp",    &m_vis_escaleUp,    "m_vis_escaleUp/F");
    tree->Branch("m_vis_escaleDown",  &m_vis_escaleDown,  "m_vis_escaleDown/F");
    tree->Branch("m_vis_jesUp",   &m_vis_jesUp,   "m_vis_jesUp/F");
    tree->Branch("m_vis_jesDown", &m_vis_jesDown, "m_vis_jesDown/F");
    tree->Branch("m_vis_resoUp",    &m_vis_resoUp,    "m_vis_resoUp/F");
    tree->Branch("m_vis_resoDown",  &m_vis_resoDown,  "m_vis_resoDown/F");
    
    tree->Branch("mTtot",        &mTtot,        "mTtot/F");
    tree->Branch("mCDF",         &mCDF,         "mCDF/F");
    tree->Branch("mTtot_muUp",   &mTtot_muUp,   "mTtot_muUp/F");
    tree->Branch("mTtot_muDown", &mTtot_muDown, "mTtot_muDown/F");
    tree->Branch("mTtot_escaleUp",    &mTtot_escaleUp,    "mTtot_escaleUp/F");
    tree->Branch("mTtot_escaleDown",  &mTtot_escaleDown,  "mTtot_escaleDown/F");
    tree->Branch("mTtot_jesUp",   &mTtot_jesUp,   "mTtot_jesUp/F");
    tree->Branch("mTtot_jesDown", &mTtot_jesDown, "mTtot_jesDown/F");
    tree->Branch("mTtot_resoUp",    &mTtot_resoUp,    "mTtot_resoUp/F");
    tree->Branch("mTtot_resoDown",  &mTtot_resoDown,  "mTtot_resoDown/F");
    
    tree->Branch("mTdileptonMET", &mTdileptonMET, "mTdileptonMET/F");
    tree->Branch("mTdileptonMET_muUp", &mTdileptonMET_muUp, "mTdileptonMET_muUp/F");
    tree->Branch("mTdileptonMET_muDown", &mTdileptonMET_muDown, "mTdileptonMET_muDown/F");
    tree->Branch("mTdileptonMET_escaleUp", &mTdileptonMET_escaleUp, "mTdileptonMET_escaleUp/F");
    tree->Branch("mTdileptonMET_escaleDown", &mTdileptonMET_escaleDown, "mTdileptonMET_escaleDown/F");
    tree->Branch("mTdileptonMET_jesUp", &mTdileptonMET_jesUp, "mTdileptonMET_jesUp/F");
    tree->Branch("mTdileptonMET_jesDown", &mTdileptonMET_jesDown, "mTdileptonMET_jesDown/F");
    tree->Branch("mTdileptonMET_resoUp", &mTdileptonMET_resoUp, "mTdileptonMET_resoUp/F");
    tree->Branch("mTdileptonMET_resoDown", &mTdileptonMET_resoDown, "mTdileptonMET_resoDown/F");

    tree->Branch("mTemu",        &mTemu,        "mTemu/F");
    tree->Branch("mTemet",       &mTemet,       "mTemet/F");
    tree->Branch("mTmumet",      &mTmumet,      "mTmumet/F");
    
    tree->Branch("m_sv",    &m_sv,   "m_sv/F");
    tree->Branch("mt_sv",   &mt_sv,  "mt_sv/F");

    tree->Branch("pt_sv",   &pt_sv,  "pt_sv/F");
    tree->Branch("eta_sv",  &eta_sv, "eta_sv/F");
    tree->Branch("phi_sv",  &phi_sv, "phi_sv/F");
    
    tree->Branch("pt_sv_escaleUp",   &pt_sv_escaleUp,  "pt_sv_escaleUp/F");
    tree->Branch("eta_sv_escaleUp",  &eta_sv_escaleUp, "eta_sv_escaleUp/F");
    tree->Branch("phi_sv_escaleUp",  &phi_sv_escaleUp, "phi_sv_escaleUp/F");
    
    tree->Branch("pt_sv_escaleDown",   &pt_sv_escaleDown,  "pt_sv_escaleDown/F");
    tree->Branch("eta_sv_escaleDown",  &eta_sv_escaleDown, "eta_sv_escaleDown/F");
    tree->Branch("phi_sv_escaleDown",  &phi_sv_escaleDown, "phi_sv_escaleDown/F");
    
    tree->Branch("pt_sv_gen",   &pt_sv_gen,  "pt_sv_gen/F");
    tree->Branch("eta_sv_gen",  &eta_sv_gen, "eta_sv_gen/F");
    tree->Branch("phi_sv_gen",  &phi_sv_gen, "phi_sv_gen/F");
    
    tree->Branch("m_sv_jesUp",   &m_sv_jesUp,   "m_sv_jesUp/F");
    tree->Branch("m_sv_jesDown", &m_sv_jesDown, "m_sv_jesDown/F");
    tree->Branch("m_sv_resoUp",    &m_sv_resoUp,    "m_sv_resoUp/F");
    tree->Branch("m_sv_resoDown",  &m_sv_resoDown,  "m_sv_resoDown/F");
    tree->Branch("m_sv_escaleUp",       &m_sv_escaleUp,       "m_sv_escaleUp/F");
    tree->Branch("m_sv_escaleDown",     &m_sv_escaleDown,     "m_sv_escaleDown/F");
    tree->Branch("m_sv_muUp",      &m_sv_muUp,      "m_sv_muUp/F");
    tree->Branch("m_sv_muDown",    &m_sv_muDown,    "m_sv_muDown/F");
    
    tree->Branch("mt_sv_jesUp",   &mt_sv_jesUp,   "mt_sv_jesUp/F");
    tree->Branch("mt_sv_jesDown", &mt_sv_jesDown, "mt_sv_jesDown/F");
    tree->Branch("mt_sv_resoUp",    &mt_sv_resoUp,    "mt_sv_resoUp/F");
    tree->Branch("mt_sv_resoDown",  &mt_sv_resoDown,  "mt_sv_resoDown/F");
    tree->Branch("mt_sv_escaleUp",       &mt_sv_escaleUp,       "mt_sv_escaleUp/F");
    tree->Branch("mt_sv_escaleDown",     &mt_sv_escaleDown,     "mt_sv_escaleDown/F");
    tree->Branch("mt_sv_muUp",      &mt_sv_muUp,      "mt_sv_muUp/F");
    tree->Branch("mt_sv_muDown",    &mt_sv_muDown,    "mt_sv_muDown/F");
    
    tree->Branch("pt_1", &pt_1, "pt_1/F");
    tree->Branch("pt_1_escaleUp", &pt_1_escaleUp, "pt_1_escaleUp/F");
    tree->Branch("pt_1_escaleDown", &pt_1_escaleDown, "pt_1_escaleDown/F");
    tree->Branch("phi_1", &phi_1, "phi_1/F");
    tree->Branch("eta_1", &eta_1, "eta_1/F");
    tree->Branch("m_1", &m_1, "m_1/F");
    tree->Branch("q_1", &q_1, "q_1/I");
    tree->Branch("iso_1", &iso_1, "iso_1/F");
    tree->Branch("mva_1", &mva_1, "mva_1/F");
    tree->Branch("d0_1", &d0_1, "d0_1/F");
    tree->Branch("d0err_1", &d0err_1, "d0err_1/F");
    tree->Branch("d0sig_1", &d0sig_1, "d0sig_1/F");
    tree->Branch("d0sig_1_abs", &d0sig_1_abs, "d0sig_1_abs/F");
    tree->Branch("d0_diff", &d0_diff, "d0_diff/F");
    tree->Branch("dZ_1", &dZ_1, "dZ_1/F");
    tree->Branch("dZerr_1", &dZerr_1, "dZerr_1/F");
    tree->Branch("dZsig_1", &dZsig_1, "dZsig_1/F");
    tree->Branch("dZsig_1_abs", &dZsig_1_abs, "dZsig_1_abs/F");
    tree->Branch("dZ_diff", &dZ_diff, "dZ_diff/F");
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
    tree->Branch("d0err_2", &d0err_2, "d0err_2/F");
    tree->Branch("d0sig_2", &d0sig_2, "d0sig_2/F");
    tree->Branch("d0sig_2_abs", &d0sig_2_abs, "d0sig_2_abs/F");
    tree->Branch("dZ_2", &dZ_2, "dZ_2/F");
    tree->Branch("dZerr_2", &dZerr_2, "dZerr_2/F");
    tree->Branch("dZsig_2", &dZsig_2, "dZsig_2/F");
    tree->Branch("dZsig_2_abs", &dZsig_2_abs, "dZsig_2_abs/F");
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
    
    tree->Branch("met_uncorr", &met_uncorr, "met_uncorr/F");
    tree->Branch("metphi_uncorr", &metphi_uncorr, "metphi_uncorr/F");
    
    tree->Branch("met_resoUp", &met_resoUp, "met_resoUp/F");
    tree->Branch("metphi_resoUp", &metphi_resoUp, "metphi_resoUp/F");
    
    tree->Branch("met_resoDown", &met_resoDown, "met_resoDown/F");
    tree->Branch("metphi_resoDown", &metphi_resoDown, "metphi_resoDown/F");
    
    tree->Branch("met_jesUp", &met_jesUp, "met_jesUp/F");
    tree->Branch("metphi_jesUp", &metphi_jesUp, "metphi_jesUp/F");
    
    tree->Branch("met_jesDown", &met_jesDown, "met_jesDown/F");
    tree->Branch("metphi_jesDown", &metphi_jesDown, "metphi_jesDown/F");
    
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
    
    tree->Branch("mvamet_jesUp", &mvamet_jesUp, "mvamet_jesUp/F");
    tree->Branch("mvametphi_jesUp", &mvametphi_jesUp, "mvametphi_jesUp/F");
    
    tree->Branch("mvamet_jesDown", &mvamet_jesDown, "mvamet_jesDown/F");
    tree->Branch("mvametphi_jesDown", &mvametphi_jesDown, "mvametphi_jesDown/F");
    
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
    tree->Branch("m_leptons", &m_leptons, "m_leptons/F");
    
    tree->Branch("pzetavis", &pzetavis, "pzetavis/F");
    tree->Branch("pzetavis_escaleUp", &pzetavis_escaleUp, "pzetavis_escaleUp/F");
    tree->Branch("pzetavis_escaleDown", &pzetavis_escaleDown, "pzetavis_escaleDown/F");
    

    tree->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
    tree->Branch("pzetamiss_escaleUp", &pzetamiss_escaleUp, "pzetamiss_escaleUp/F");
    tree->Branch("pzetamiss_escaleDown", &pzetamiss_escaleDown, "pzetamiss_escaleDown/F");

    tree->Branch("dzeta",&dzeta,"dzeta/F");    
    tree->Branch("dzeta_escaleUp",&dzeta_escaleUp,"dzeta_escaleUp/F");
    tree->Branch("dzeta_escaleDown",&dzeta_escaleDown,"dzeta_escaleDown/F");


    tree->Branch("pzetamiss_resoUp", &pzetamiss_resoUp, "pzetamiss_resoUp/F");
    tree->Branch("dzeta_resoUp",&dzeta_resoUp,"dzeta_resoUp/F");
    
    tree->Branch("pzetamiss_resoDown", &pzetamiss_resoDown, "pzetamiss_resoDown/F");
    tree->Branch("dzeta_resoDown",&dzeta_resoDown,"dzeta_resoDown/F");
    
    tree->Branch("pzetamiss_jesUp", &pzetamiss_jesUp, "pzetamiss_jesUp/F");
    tree->Branch("dzeta_jesUp",&dzeta_jesUp,"dzeta_jesUp/F");
    
    tree->Branch("pzetamiss_jesDown", &pzetamiss_jesDown, "pzetamiss_jesDown/F");
    tree->Branch("dzeta_jesDown",&dzeta_jesDown,"dzeta_jesDown/F");
    


    tree->Branch("pzetamiss_mvamet", &pzetamiss_mvamet, "pzetamiss_mvamet/F");
    tree->Branch("dzeta_mvamet",&dzeta_mvamet,"dzeta_mvamet/F");
    
    tree->Branch("pzetamiss_mvamet_uncorr", &pzetamiss_mvamet_uncorr, "pzetamiss_mvamet_uncorr/F");
    tree->Branch("dzeta_mvamet_uncorr",&dzeta_mvamet_uncorr,"dzeta_mvamet_uncorr/F");
    
    tree->Branch("pzetamiss_mvamet_resoUp", &pzetamiss_mvamet_resoUp, "pzetamiss_mvamet_resoUp/F");
    tree->Branch("dzeta_mvamet_resoUp",&dzeta_mvamet_resoUp,"dzeta_mvamet_resoUp/F");
    
    tree->Branch("pzetamiss_mvamet_resoDown", &pzetamiss_mvamet_resoDown, "pzetamiss_mvamet_resoDown/F");
    tree->Branch("dzeta_mvamet_resoDown",&dzeta_mvamet_resoDown,"dzeta_mvamet_resoDown/F");
    
    tree->Branch("pzetamiss_mvamet_jesUp", &pzetamiss_mvamet_jesUp, "pzetamiss_mvamet_jesUp/F");
    tree->Branch("dzeta_mvamet_jesUp",&dzeta_mvamet_jesUp,"dzeta_mvamet_jesUp/F");
    
    tree->Branch("pzetamiss_mvamet_jesDown", &pzetamiss_mvamet_jesDown, "pzetamiss_mvamet_jesDown/F");
    tree->Branch("dzeta_mvamet_jesDown",&dzeta_mvamet_jesDown,"dzeta_mvamet_jesDown/F");
    

    tree->Branch("pzetamiss_genmet", &pzetamiss_genmet, "pzetamiss_genmet/F");
    tree->Branch("dzeta_genmet",&dzeta_genmet,"dzeta_genmet/F");
    
    tree->Branch("mva_gf", &mva_gf, "mva_gf/F");
    
    tree->Branch("njets", &njets, "njets/I");
    tree->Branch("njets_jesUp", &njets_jesUp, "njets_jesUp/I");
    tree->Branch("njets_jesDown", &njets_jesDown, "njets_jesDown/I");


    tree->Branch("njetspt20", &njetspt20, "njetspt20/I");

    tree->Branch("bdt",&bdt,"bdt/F");
    tree->Branch("bdt_bbh",&bdt_bbh,"bdt_bbh/F");
    tree->Branch("bdt_ggh",&bdt_ggh,"bdt_ggh/F");
    
    tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
    tree->Branch("jpt_1_jesUp", &jpt_1_jesUp, "jpt_1_jesUp/F");
    tree->Branch("jpt_1_jesDown",&jpt_1_jesDown, "jpt_1_jesDown/F");

    tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
    tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
    tree->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
    tree->Branch("jptunc_1", &jptunc_1, "jptunc_1/F");
    tree->Branch("jmva_1", &jmva_1, "jmva_1/F");
    tree->Branch("jlrm_1", &jlrm_1, "jlrm_1/F");
    tree->Branch("jctm_1", &jctm_1, "jctm_1/I");
    
    tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
    tree->Branch("jpt_2_jesUp", &jpt_2_jesUp, "jpt_2_jesUp/F");
    tree->Branch("jpt_2_jesDown", &jpt_2_jesDown, "jpt_2_jesDown/F");

    tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
    tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
    tree->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
    tree->Branch("jptunc_2", &jptunc_2, "jptunc_2/F");
    tree->Branch("jmva_2", &jmva_2, "jlrm_2/F");
    tree->Branch("jlrm_2", &jlrm_2, "jlrm_2/F");
    tree->Branch("jctm_2", &jctm_2, "jctm_2/I");
    
    tree->Branch("mjj", &mjj, "mjj/F");
    tree->Branch("mjj_jesUp", &mjj_jesUp, "mjj_jesUp/F");
    tree->Branch("mjj_jesDown", &mjj_jesDown, "mjj_jesDown/F");

    tree->Branch("jdeta", &jdeta, "jdeta/F");
    tree->Branch("njetingap", &njetingap, "njetingap/I");
    
    tree->Branch("nbtag", &nbtag, "nbtag/I");
    tree->Branch("nbtag_noSF", &nbtag_noSF, "nbtag_noSF/I");
    tree->Branch("bpt_1",   &bpt,   "bpt/F");
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

    //set JES uncertainties, TO DO: UPDATE FOR 2017
    JESUncertainties * jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Summer16_UncertaintySources_AK4PFchs.txt");
    std::vector<std::string> uncertNames = jecUncertainties->getUncertNames();

    std::cout << "Number of uncertainties = " << uncertNames.size() << std::endl;

    int njetsUncUp[30];
    int njetsUncDown[30];
    float mjjUncUp[30];
    float mjjUncDown[30];
    
    int iUncert = 0;     //jec uncertainties, leaf for njets and mjj distributions varied up&down for each unc 
    if (!isData) {
       for (auto const& Name : uncertNames) {
          TString name(Name);
          cout << name << endl;
          tree->Branch("njets_"+name+"Up",&njetsUncUp[iUncert],"njets_"+name+"Up/I");
          tree->Branch("njets_"+name+"Down",&njetsUncDown[iUncert],"njets_"+name+"Down/I");
          tree->Branch("mjj_"+name+"Up",&mjjUncUp[iUncert],"njets_"+name+"Up/F");
          tree->Branch("mjj_"+name+"Down",&mjjUncDown[iUncert],"njets_"+name+"Down/F");
          iUncert++;
       }
       cout << endl;
    }

    /*    
    treeGen->Branch("bosonPt",&bosonPt,"bosonPt/F");
    treeGen->Branch("bosonMass",&bosonMass,"bosonMass/F");
    treeGen->Branch("bosonEta",&bosonEta,"bosonEta/F");
    
    treeGen->Branch("higgsPt",&higgsPt,"higgsPt/F");
    treeGen->Branch("higgsMass",&higgsMass,"higgsMass/F");
    treeGen->Branch("higgsEta",&higgsEta,"higgsEta/F");
    */
    
    int nTotalFiles = 0;
    std::string dummy;
    // count number of files --->
    while (fileList0 >> dummy) nTotalFiles++;
    std::vector<Period> periods;
    string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
    
    if (isData) {                                                     //read in periods from JSON file
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
       TString NamePUHistMC = samplenameForPUHist;
       TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
       TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(NamePUHistMC);
       PUofficial->set_h_data(PU_data);               
       PUofficial->set_h_MC(PU_mc);   
    }                   //histograms are set, ratio built, no weight calculated so far
       
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
    
    // MEt filters                                 //set MET Filters
   std::vector<TString> metFlags; metFlags.clear();
   metFlags.push_back("Flag_HBHENoiseFilter");//ok
   metFlags.push_back("Flag_HBHENoiseIsoFilter");//ok
   metFlags.push_back("Flag_globalTightHalo2016Filter"); //ok
   metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");//ok
   metFlags.push_back("Flag_goodVertices");  //ok
   if (isData)
     metFlags.push_back("Flag_eeBadScFilter");//ok
   metFlags.push_back("Flag_BadPFMuonFilter");//ok
   metFlags.push_back("Flag_BadChargedCandidateFilter");//ok
   metFlags.push_back("Flag_ecalBadCalibFilter");//ok

   std::vector<TString> badGlobalMuonFlag; badGlobalMuonFlag.clear();
   badGlobalMuonFlag.push_back("Flag_BadGlobalMuonFilter");
   std::vector<TString> muonBadTrackFlag; muonBadTrackFlag.clear();
   muonBadTrackFlag.push_back("Flag_muonBadTrackFilter");
   std::vector<TString> chargedHadronTrackResolutionFlag; chargedHadronTrackResolutionFlag.clear();
   chargedHadronTrackResolutionFlag.push_back("Flag_chargedHadronTrackResolutionFilter");
   
   RecoilCorrector recoilMvaMetCorrector(RecoilMvaFileName); //improve description of MET in MC events without genuine MET
   MEtSys metSys(MetSysFileName);
    
   RecoilCorrector recoilMetCorrector(RecoilFileName);
    
   // SV fit mass
   edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root"); //TO DO: NOT NEEDED ANYMORE?
   TH1::AddDirectory(false);
   TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
    
    // qcd weight (dzeta cut)
   QCDModelForEMu qcdWeight("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_2016BtoH.root"); 
   // qcd weight DZeta cut
   QCDModelForEMu qcdWeightNoDzeta("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_2016BtoH.root");
   
   // BTag scale factors
   BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_Moriond17_B_H.csv");  //read out SF for btagging
   BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central",{"up","down"});
   BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central",{"up","down"});
   BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central",{"up","down"});
   reader_B.load(calib,BTagEntry::FLAV_B,"comb");
   reader_C.load(calib,BTagEntry::FLAV_C,"comb");
   reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");
    
   float etaBTAG[2] = {0.5,2.1};
   float ptBTAG[5] = {25.,35.,50.,100.,200.};
   
   std::cout << std::endl;
   for (int iEta=0; iEta<2; ++iEta) {
     for (int iPt=0; iPt<5; ++iPt) {
       float sfB = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
       float sfC = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
       float sfLight = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
       printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
     }
   }
   std::cout << std::endl;
   // read out efficiencies for btagging
   TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_Moriond2017.root"));
   TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
   TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
   TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
   TRandom3 rand;
   
   float MaxBJetPt = 1000.;
   float MaxLJetPt = 1000.;
   float MinLJetPt = 20.;
   float MinBJetPt = 20.; // !!!!!
   
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
   
   // correction WS
   TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/"+correctionWSFile;
   TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
   RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");
    
    int nEvents = 0;
    int selEvents = 0;
    int nFiles = 0;
        
    for (int iF=0; iF<nTotalFiles; ++iF) { // loop over input file names
        
       std::string filen;
       fileList >> filen;
       
       std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
       TFile * file_ = TFile::Open(TString(filen));
        
       TTree * _inittree = NULL;                       // fill genweights  for MC
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
       _tree = (TTree*)file_->Get(TString(ntupleName));  //check if makeroortree is available
       if (_tree==NULL) continue;
       
       TH1D * histoInputEvents = NULL;
       histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
       if (histoInputEvents==NULL) continue;
       int NE = int(histoInputEvents->GetEntries());
       std::cout << "      number of input events    = " << NE << std::endl;
       for (int iE=0;iE<NE;++iE)                        //fill histogram with number of events
          inputEventsH->Fill(0.);
       
       AC1B analysisTree(_tree);
       
       Long64_t numberOfEntries = analysisTree.GetEntries();
       
       std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
       
       for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {
          
          analysisTree.GetEntry(iEntry);                //store maxRun and minRun
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
          
          // weights
          mcweight = analysisTree.genweight;           //store genweights
          
          weightScale1 = analysisTree.weightScale1;    //weights for scale uncertainty
          weightScale2 = analysisTree.weightScale2;
          weightScale3 = analysisTree.weightScale3;
          weightScale4 = analysisTree.weightScale4;
          weightScale5 = analysisTree.weightScale5;
          weightScale6 = analysisTree.weightScale6;
          weightScale7 = analysisTree.weightScale7;
          weightScale8 = analysisTree.weightScale8;
          
          weightPDFup   = analysisTree.weightPDFup;   //store pdf weights
          weightPDFdown = analysisTree.weightPDFdown;
          
          puweight = 1;                               //set default values for variables
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
          
          btag0weight = 1;
          btag0weight_Up = 1;
          btag0weight_Down = 1;
          btag1weight = 1;
          btag1weight_Up = 1;
          btag1weight_Down = 1;
          btag2weight = 1;
          btag2weight_Up = 1;
          btag2weight_Down = 1;
          
          qcdweight = 1;
          qcdweightup = 1;
          qcdweightdown = 1;
          qcdweight_nodzeta = 1;
          qcdweightup_nodzeta = 1;
          qcdweightdown_nodzeta = 1;
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
          bosonEta = 0;
          bosonMass = -1;
          higgsPt = 0;
          higgsEta = 0;
          higgsMass = -1;

          
          metFilters_ = true;
          
          badGlobalMuonFilter_ = true;
          muonBadTrackFilter_ = true;
          chargedHadronTrackResolutionFilter_ = true;
          
          badMuonFilter_ = true;
          duplicateMuonFilter_ = true;
          
          float topPt = -1;
          float antitopPt = -1;
          
          bool isZfound = false;
          bool isWfound = false;
          bool isHfound = false;
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
             // store generator taus from LV, also for visible tau
	     //	std::cout << "gentaus : " << analysisTree.gentau_count << std::endl;
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
             //store genparticles
             for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
                
                TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                    analysisTree.genparticles_py[igen],
                                                    analysisTree.genparticles_pz[igen],
                                                    analysisTree.genparticles_e[igen]);
                //store pt of top and antitop 
                if (analysisTree.genparticles_pdgid[igen]==6)
                   topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                       analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
                
                if (analysisTree.genparticles_pdgid[igen]==-6)
                   antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                           analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
                //store GStar
                if (analysisTree.genparticles_pdgid[igen]==22 && analysisTree.genparticles_status[igen]==44)
                   isGSfound = true;
                //store Z boson
                if (analysisTree.genparticles_pdgid[igen]==23) {
                   isZfound = true;
                   zBosonLV = genLV;
                }
                //store higgs bosons
                if (analysisTree.genparticles_pdgid[igen]==25||
                    analysisTree.genparticles_pdgid[igen]==35||
                    analysisTree.genparticles_pdgid[igen]==36) {
                   isHfound = true;
                   hBosonLV = genLV;
                }
                //store W boson
                if (abs(analysisTree.genparticles_pdgid[igen])==24) {
                   isWfound = true;
                   wBosonLV = genLV;
                }
                
                bool fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1;
                bool isMuon = false;
                bool isElectron = false;
                bool isNeutrino = false;
                bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];
                //store prompt electrons, add electron LV to promptElectronsLV and wDecayProductsLV
                if (abs(analysisTree.genparticles_pdgid[igen])==11) {
                   isElectron = true;
                   if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                      promptElectrons.push_back(genLV);
                      promptElectronsLV += genLV;
                      wDecayProductsLV += genLV;
                   }
                }
                //store prompt muons, add electron LV to promptElectronsLV and wDecayProductsLV
                if (abs(analysisTree.genparticles_pdgid[igen])==13) {
                   isMuon = true;
                   if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                      promptMuons.push_back(genLV);
                      promptMuonsLV += genLV;
                      wDecayProductsLV += genLV;
                   }
                }
                //store prompt muons, add electron LV to promptElectronsLV and wDecayProductsLV, distiguish between neutrinos from taus and others
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
                
             } //end of loop over genparticles
             
             if (isGSfound) {      //remove events where a gamma star has been found if removeGammaStar=true
                //	  std::cout << "gamma* found : " << std::endl;
                if (removeGammaStar) continue;
             }
             
             if (isDY||isEmbedded) { //if DY than check which decay mode and store gen-Info of the Boson
                
                if (promptTausFirstCopy.size()==2) {
                   isZTT = true; isZMM = false; isZEE = false;
                   bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
                   bosonMass = promptTausLV.M();
                   bosonEta  = promptTausLV.Eta();
                   lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
                   mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);
		   double gt1_pt  = promptTausFirstCopy[0].Pt();
		   double gt1_eta = promptTausFirstCopy[0].Eta();
		   double gt2_pt  = promptTausFirstCopy[1].Pt();
		   double gt2_eta = promptTausFirstCopy[1].Eta();
		   
		   correctionWS->var("gt_pt")->setVal(gt1_pt);
		   correctionWS->var("gt_eta")->setVal(gt1_eta);
		   double id1_embed = correctionWS->function("m_sel_idEmb_ratio")->getVal();
		   correctionWS->var("gt_pt")->setVal(gt2_pt);
                   correctionWS->var("gt_eta")->setVal(gt2_eta);
		   double id2_embed = correctionWS->function("m_sel_idEmb_ratio")->getVal();
		   correctionWS->var("gt1_pt")->setVal(gt1_pt);
                   correctionWS->var("gt1_eta")->setVal(gt1_eta);
		   correctionWS->var("gt2_pt")->setVal(gt2_pt);
                   correctionWS->var("gt2_eta")->setVal(gt2_eta);
		   double trg_embed = correctionWS->function("m_sel_trg_ratio")->getVal();
		   embeddedWeight = id1_embed * id2_embed * trg_embed;
		   //              std::cout << "id1 embedded  = " << id1_embed << std::endl;
		   //              std::cout << "id2 embedded  = " << id2_embed << std::endl;
		   //              std::cout << "trig embedded = " << trg_embed << std::endl;
		   //              std::cout << "embedded weight = " << embeddedWeight << std::endl;
                }
                else if (promptMuons.size()==2) {
                   isZTT = false; isZMM = true; isZEE = false;
                   bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
                   bosonMass = promptMuonsLV.M();
                   bosonEta = promptMuonsLV.Eta();
                   lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
                   mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
                }
                else {
                   isZTT = false; isZMM = false; isZEE = true;
                   bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
                   bosonMass = promptElectronsLV.M();
                   bosonEta = promptElectronsLV.Eta();
                   lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
                   if (promptElectrons.size()==2)
                      mtBoson_gen = mT(promptElectrons[0],promptElectrons[1]);
                }
                //store momentum from neutrino of tau lepton decays
                nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
                }
             else if (isW) { //store geninfo for W+jets events
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
             else { //store Higgs boson kinematics
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
             
             // bosonPt = TMath::Sqrt(bosonPx*bosonPx+bosonPy*bosonPy);
             /*
               if (isHfound) {
               higgsMass = hBosonLV.M();
               higgsPt   = hBosonLV.Pt();
               higgsEta  = hBosonLV.Eta();
               }
               else {
               bosonPt = -9999;
               higgsMass = -9999;
               higgsPt = -9999;
               higgsEta = -9999;
               }
             */
             //                treeGen->Fill();
             
             /*
               std::cout << "Taus (first copy) : " << promptTausFirstCopy.size() << std::endl;
               std::cout << "Taus (last copy)  : " << promptTausLastCopy.size() << std::endl;
               std::cout << "Muons             : " << promptMuons.size() << std::endl;
               std::cout << "Electrons         : " << promptElectrons.size() << std::endl;
               std::cout << "Neutrinos         : " << promptNeutrinos.size() << std::endl;
               std::cout << "Neutrinos (taus)  : " << tauNeutrinos.size() << std::endl;
               if (isZfound)
               std::cout << "ZBoson  px = " << zBosonLV.Px()
               << "        py = " << zBosonLV.Py()
               << "        pz = " << zBosonLV.Pz() << std::endl;
               
               if (isHfound)
               std::cout << "HBoson  px = " << hBosonLV.Px()
               << "        py = " << hBosonLV.Py()
               << "        pz = " << hBosonLV.Pz() << std::endl;
               std::cout << "  Px = " << bosonPx 
               << "  Py = " << bosonPy
               << "  Pz = " << bosonPz << std::endl;
               std::cout << "  lepPx = " << lepPx 
               << "  lepPy = " << lepPy
               << "  lepPz = " << lepPz << std::endl;
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
               std::cout << "Nu(W)   pz = " << promptNeutrinosLV.Px()
               << "        py = " << promptNeutrinosLV.Py()
               << "        pz = " << promptNeutrinosLV.Pz() << std::endl;
               std::cout << "Nu(tau) px = " << tauNeutrinosLV.Px()
               << "        py = " << tauNeutrinosLV.Py()
               << "        pz = " << tauNeutrinosLV.Px() << std::endl;
               std::cout << "Nu(x)   px = " << nuPx
               << "        py = " << nuPy
               << "        pz = " << nuPz << std::endl;
             */
             
             
             
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
                   // correctionWS->var("z_gen_pt")->setVal(bosonPtX);
                   // correctionWS->var("z_gen_mass")->setVal(bosonMassX);
                   // zptmassweight = correctionWS->function("zpt_weight_nom")->getVal();
                   // zptmassweight_esup = correctionWS->function("zpt_weight_esup")->getVal();
                   // zptmassweight_esdown = correctionWS->function("zpt_weight_esdown")->getVal();
                   // zptmassweight_ttup = correctionWS->function("zpt_weight_ttup")->getVal();
                   // zptmassweight_ttdown = correctionWS->function("zpt_weight_ttdown")->getVal();
                   // zptmassweight_statpt0up = correctionWS->function("zpt_weight_statpt0up")->getVal();
                   // zptmassweight_statpt0down = correctionWS->function("zpt_weight_statpt0down")->getVal();
                   // zptmassweight_statpt40up = correctionWS->function("zpt_weight_statpt40up")->getVal();
                   // zptmassweight_statpt40down = correctionWS->function("zpt_weight_statpt40down")->getVal();
                   // zptmassweight_statpt80up = correctionWS->function("zpt_weight_statpt80up")->getVal();
                   // zptmassweight_statpt80down = correctionWS->function("zpt_weight_statpt80down")->getVal();
                    /*
                     if (bosonPtX>500||bosonMassX>500) {
                     std::cout << "boson pt = " << bosonPtX << "  mass = " << bosonMassX << std::endl;
                     std::cout << "  weight = " << zptmassweight << std::endl;
                     std::cout << "  weight_esup = " << zptmassweight_esup << "  weight_esdown = " << zptmassweight_esdown << std::endl;
                     std::cout << "  weight_ttup = " << zptmassweight_ttup << "  weight_ttdown = " << zptmassweight_ttdown << std::endl;
                     std::cout << "  weight_statpt0up = " << zptmassweight_statpt0up << "  weight_statpt0down = " << zptmassweight_statpt0down << std::endl;
                     std::cout << "  weight_statpt40up = " << zptmassweight_statpt40up << "  weight_statpt40down = " << zptmassweight_statpt40down << std::endl;
                     std::cout << "  weight_statpt80up = " << zptmassweight_statpt80up << "  weight_statpt80down = " << zptmassweight_statpt80down << std::endl;
                     } 
                   */
                }
                
                
             }
             //fill the Z boson decay mode into a histogram
             ALL->Fill(0.0);
             if (isZMM) ZMM->Fill(0.);
             if (isZEE) ZEE->Fill(0.);
             if (isZTT) ZTT->Fill(0.);
             isZLL = isZMM || isZEE;
             zptmassweight =1.0;
          }
            
            //  if (isZEE&&isZMM) {
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
            
            //  cout << "Ok" << endl;
            
	    


          run = int(analysisTree.event_run);
          lumi = int(analysisTree.event_luminosityblock);
          evt = int(analysisTree.event_nr);

	  // embedded weight

          if (debug) std::cout<<"check good run selection"<<std::endl;
          if (isData && applyGoodRunSelection) { //apply good run selection, periods: vector filled from JSON file
             if (debug) std::cout<<"in good run list"<<std::endl;
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
             if (debug) std::cout<<"Event passed good run selection"<<std::endl;
          }
          //std::cout << "passed lumi" << endl;
            
          npv = analysisTree.primvertex_count;
          npu = analysisTree.numtruepileupinteractions;
          rho = analysisTree.rho;
          
          npartons = analysisTree.genparticles_noutgoing;
          
          if (!isData && !isEmbedded) { //PU re-weighting, top pT re-reweighting, store weights
             puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
                     
             //std::cout<<"puweight applied"<<std::endl;
             //      cout << "puweight = " << puweight << endl;
             if (topPt>0&&antitopPt>0) {
                topptweight = topPtWeight(topPt,antitopPt,true);
                topptweightRun2 = topPtWeight(topPt,antitopPt,false);
                	  // std::cout << "topPt = " << topPt
                	  //      << "   antitopPt = " << antitopPt
                	  //      << "   weight = " << topptweight << std::endl;
             }
             histWeightsSkimH->Fill(double(0),double(mcweight)); 
             histWeightsTTH->Fill(double(0),double(mcweight*topptweight));
          }
          else {
             histWeightsSkimH->Fill(double(0),double(1));
             histWeightsTTH->Fill(double(0),double(1));
          }
          //std::cout<<"before seraching filters"<<std::endl;
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
          //std::cout << "nfiltres = " << nfilters << std::endl;
          for (unsigned int i=0; i<nfilters; ++i) {
             //	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
             //search for triggers in HLTFilter
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
                    continue;
                }
                if (!isHighPtLegElectron) {
                    std::cout << "HLT filter " << HighPtLegElectron << " not found" << std::endl;
                    continue;
                }
                if (!isLowPtLegMuon) {
                    std::cout << "HLT filter " << LowPtLegMuon << " not found" << std::endl;
                    continue;
                }
                if (!isHighPtLegMuon) {
                    std::cout << "HLT filter " << HighPtLegMuon << " not found" << std::endl;
                    continue;
                }
                if (applyDzFilterMatch) {
                   if (!isMu23Ele12DzFilter) {
                      std::cout << "HLT filter " << Mu23Ele12DzFilter << " not found" << std::endl;
                      continue;
                   }
                   if (!isMu8Ele23DzFilter) {
                      std::cout << "HLT filter " << Mu8Ele23DzFilter << " not found" << std::endl;
                      continue;
                   }
                }
            }
            //std::cout<<"before btag"<<std::endl;
            //find the right b tag discriminant  
            unsigned int nBTagDiscriminant = 0;
            for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
               TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
               if (discr.Contains(BTagDiscriminator))
                  nBTagDiscriminant = iBTag;
            }
            
            //std::cout<<"before met filters"<<std::endl; 
            // MET Filters //store if event passes the met filters
            metFilters_ = metFiltersPasses(analysisTree,metFlags);
            //badGlobalMuonFilter_ = metFiltersPasses(analysisTree,badGlobalMuonFlag);
            //muonBadTrackFilter_ = metFiltersPasses(analysisTree,muonBadTrackFlag);
            //chargedHadronTrackResolutionFilter_ = metFiltersPasses(analysisTree,chargedHadronTrackResolutionFlag);
   
            /*
              for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
              bool lowPtMuon = false;
              bool highPtMuon = false;
              if (analysisTree.trigobject_filters[iT][nHighPtLegMuon])
              highPtMuon = true;
              if (analysisTree.trigobject_filters[iT][nLowPtLegMuon])
              lowPtMuon = true;
              
              if (highPtMuon)
              std::cout << "High pT Muon leg ";
              if (highPtMuon)
              std::cout << "Low pT Muon leg ";
              
              if (highPtMuon||lowPtMuon)
              std::cout << "TO : pT = " << analysisTree.trigobject_pt[iT] 
              << "  eta = " << analysisTree.trigobject_eta[iT]
              << "  phi = " << analysisTree.trigobject_phi[iT] << std::endl;
              
              }
	    
              
              std::cout << "LowPtE  : " << LowPtLegElectron << " : " << nLowPtLegElectron << std::endl;
              std::cout << "HighPtE : " << HighPtLegElectron << " : " << nHighPtLegElectron << std::endl;
              std::cout << "LowPtM  : " << LowPtLegMuon << " : " << nLowPtLegMuon << std::endl;
              std::cout << "HighPtM : " << HighPtLegMuon << " : " << nHighPtLegMuon << std::endl;
              std::cout << BTagDiscriminator << " : "<< nBTagDiscriminant << std::endl;
              std::cout << std::endl;
            */

            //      continue;
            
            // vertex cuts no longer required
            //      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
            //      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
            //      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
            //		       analysisTree.primvertex_y*analysisTree.primvertex_y);
            //      if (dVertex>dVertexCut) continue;
            
            
            //                         HLT_Mu23_Ele12                     ||      HLT_Mu8_Ele23
            //      bool trigAccept = (analysisTree.hltriggerresults_second[5]==1)||(analysisTree.hltriggerresults_second[6]==1);
            //      if (!trigAccept) continue;
            
            //      if (nonOverlap&&checkOverlap)
            //      cout << "Muons = " << analysisTree.muon_count
            //	   << "   Electrons = " << analysisTree.electron_counte << std::endl;
            
            
            // electron selection
            // loop over electrons and check if event passes electron kinematic cuts and ID criteria
            vector<int> electrons; electrons.clear();
            //	    cout << "Number of electrons = " << analysisTree.electron_count << endl;
            for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
               if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
               if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
               if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
               if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
               //	bool electronMvaId = electronMvaIdWP80(analysisTree.electron_pt[ie],
               //					       analysisTree.electron_superclusterEta[ie],
               //					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
               bool electronMvaId = true;
               if (applyIsoElectronId) electronMvaId = analysisTree.electron_mva_wp80_Iso_Fall17_v1[ie]>0.5;
              else electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[ie]>0.5;
               if (!electronMvaId) continue;
               if (!analysisTree.electron_pass_conversion[ie]) continue;
               if (analysisTree.electron_nmissinginnerhits[ie]>1) continue;
               electrons.push_back(ie);
            }
            
            // muon selection
            // loop over muons and check if event passes electron kinematic cuts and ID criteria
            vector<int> muons; muons.clear();
            //	    cout << "Number of muons = " << analysisTree.muon_count << std::endl;
            for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
               if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
               if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
               if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
               if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
               /*
                bool goodGlobal =
                analysisTree.muon_isGlobal[im] &&
                analysisTree.muon_normChi2[im] < 3 &&
                analysisTree.muon_combQ_chi2LocalPosition[im] < 12 &&
                analysisTree.muon_combQ_trkKink[im] < 20;
                bool isICHEPmedium  =
                analysisTree.muon_isLoose[im] &&
                analysisTree.muon_validFraction[im] >0.49 &&
                analysisTree.muon_segmentComp[im] > (goodGlobal ? 0.303 : 0.451);
               */
               if (analysisTree.muon_isBad[im]) badMuonFilter_ = false;
               if (analysisTree.muon_isDuplicate[im]) duplicateMuonFilter_ = false;
               bool muonId = analysisTree.muon_isMedium[im];
               if (applyICHEPMuonId) muonId = analysisTree.muon_isICHEP[im];
               //	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
               if (!muonId) continue;
               muons.push_back(im);
            }
            //only store events that pass the met filters and that contain an electron or muon
           metFilters_ = metFilters_ ;////&& badMuonFilter_; //&& duplicateMuonFilter_;
          //if (metFilters_<0.5) continue;
          
          if (debug) std::cout<<"Event passed met filters"<<std::endl;
            //	    cout << "  SelEle=" << electrons.size()
            //		 << "  SelMu=" << muons.size() << std::endl;
            if (electrons.size()==0) continue;
            if (debug) std::cout<<"Event has at least one electron"<<std::endl;
            if (muons.size()==0) continue;
            if (debug) std::cout<<"Event has at least one muon"<<std::endl;
            // selecting muon and electron pair (OS or SS);
            int electronIndex = -1;
            int muonIndex = -1;
            
            float isoMuMin = 1e+10;
            float isoEleMin = 1e+10;
            //      if (muons.size()>1||electrons.size()>1)
            //      std::cout << "muons = " << muons.size() << "  electrons = " << electrons.size() << std::endl;
            bool isMuon23matched = false;
            bool isMuon8matched  = false;
            bool isElectron23matched = false;
            bool isElectron12matched = false;
            //calcul;ate isolation of the muon in default cone and in deltaR =0.3
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
               
               //		cout << "muon " << im
               //		     << "   pt = " << analysisTree.muon_pt[im] << "   eta = " << analysisTree.muon_eta[im]
               //		     << " : isMu23 = " << isMu23 << "  isMu8 = " << isMu8 << std::endl;
               
               for (unsigned int ie=0; ie<electrons.size(); ++ie) {
                  
                  unsigned int eIndex = electrons.at(ie);
                  
                  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
                                    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);
                  
                  
                  //		    std::cout << "deltaR = " << dR << std::endl; 
                  
                  if (dR<dRleptonsCut) continue;
                  //cut on Delta R between the electron an the muons
                  if (debug) std::cout<<"Event passes distance cut between electrons and muons"<<std::endl;
                  
                  //		    cout << "electron " << ie 
                  //			 << "   pt = " << analysisTree.electron_pt[ie] << "   eta = " << analysisTree.electron_eta[ie]
                  //			 << " : isEle23 = " << isEle23 << "  isEle12 = " << isEle12 << std::endl;
                  
                  //Is there a lepton passing the high lepton pT cuts
                  bool trigMatch =
                     (analysisTree.muon_pt[mIndex]>ptMuonHighCut) ||
                     (analysisTree.electron_pt[eIndex]>ptElectronHighCut);
                  
                  //if (!trigMatch) continue;
                  if (debug)  std::cout<<"Event has been triggered"<<std::endl;
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
                  
                  //isolation variable for electrons, switched to rho isolation
                  // float neutralHadIsoEle = analysisTree.electron_neutralHadIso[eIndex];
                  // float photonIsoEle = analysisTree.electron_photonIso[eIndex];
                  // float chargedHadIsoEle = analysisTree.electron_chargedHadIso[eIndex];
                  // float puIsoEle = analysisTree.electron_puIso[eIndex];
                  // if (isElectronIsoR03) {
                  //    neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
                  //    photonIsoEle = analysisTree.electron_r03_sumPhotonEt[eIndex];
                  //    chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
                  //    puIsoEle = analysisTree.electron_r03_sumPUPt[eIndex];
                  // }
                  // float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
                  // neutralIsoEle = TMath::Max(float(0),neutralIsoEle);
                  // float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
                  // float relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];
                               
                  float rhoNeutral = analysisTree.rho;
                  float  eA = getEffectiveArea( fabs(analysisTree.electron_superclusterEta[eIndex]) );
                  float absIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex] +
		    TMath::Max(0.0f,analysisTree.electron_r03_sumNeutralHadronEt[eIndex]+analysisTree.electron_r03_sumPhotonEt[eIndex]-eA*rhoNeutral);
                  float relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];
                  
                  //store most isolated muons and electrons
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
            
            //	    cout << "mIndex = " << muonIndex << "   eIndex = " << electronIndex << std::endl;
            
            
            if (electronIndex<0) continue;
            if (muonIndex<0) continue;
            if (debug) std::cout<<"Event has an isolated lepton"<<std::endl;
            
            bool isMu23 = false;
            bool isMu8 = false;
            bool isMu23dz = false;
            bool isMu8dz  = false;
          
            //loop over trigger objects
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
               float dRtrig = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
                                     analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
               //matched to trigger object? either low or high leg
               if (dRtrig<deltaRTrigMatch) {
                  if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]&&
                      analysisTree.muon_pt[muonIndex]>ptMuonHighCut) { // Mu23 Leg
                     isMu23 = true;
                  }
                  if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]&&
                      analysisTree.muon_pt[muonIndex]>ptMuonLowCut) { // Mu8 Leg
                     isMu8 = true;
                  }
                  if (analysisTree.trigobject_filters[iT][nMu23Ele12DzFilter]){
                     isMu23dz = true;
                  }
                  if (analysisTree.trigobject_filters[iT][nMu8Ele23DzFilter]){
                     isMu8dz = true;
                  }
               }
            }
            //if applyDzFilterMatch=true check if matched to either Mu or Mudz 
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
               //matched to trigger object? either low or high leg
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
                     {
                        isEle12dz = true;
                     }
                  if (analysisTree.trigobject_filters[iT][nMu8Ele23DzFilter]){
                     isEle23dz = true;
                  }
               }
            }
            //if applyDzFilterMatch=true check if matched to either Ele or Eledz 
            if (applyDzFilterMatch) {
               isEle23 = isEle23 && isEle23dz;
               isEle12 = isEle12 && isEle12dz;
            }		      
            trg_muonelectron = 
               (isMu23&&isEle12&&analysisTree.muon_pt[muonIndex]>ptMuonHighCut) ||
               (isMu8&&isEle23&&analysisTree.electron_pt[electronIndex]>ptElectronHighCut);
            //event stored if trigger matching is true and a match to one of the triggers has been found
            //if (applyTriggerMatch&&trg_muonelectron<0.5) continue;
            if (debug) std::cout<<"A trigger match has been found"<<std::endl;
            
            //OS Muon-Ele pair?
            os = (analysisTree.muon_charge[muonIndex]*analysisTree.electron_charge[electronIndex]) < 0;
    
            // looking for extra electron
            bool foundExtraElectron = false;
            for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
               if (int(ie)==electronIndex) continue;
               if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
               if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
               if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
               if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
               bool electronMvaId = true;
               if (applyIsoElectronId) electronMvaId = analysisTree.electron_mva_wp90_Iso_Fall17_v1[ie]>0.5;
	       else electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[ie]>0.5;
               if (!electronMvaId&&applyVetoElectronId) continue;
               if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
               if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
               // float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
               // float photonIsoEle = analysisTree.electron_photonIso[ie];
               // float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
               // float puIsoEle = analysisTree.electron_puIso[ie];
               // if (isElectronIsoR03) {
               //    neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
               //    photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
               //    chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
               //    puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
               // }
               // float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
               // neutralIsoEle = TMath::Max(float(0),neutralIsoEle);
               // float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
               // float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];
               
               float rhoNeutral = analysisTree.rho;
               float  eA = getEffectiveArea( fabs(analysisTree.electron_superclusterEta[ie]) );
               float absIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie] +
                  TMath::Max(0.0f,analysisTree.electron_r03_sumNeutralHadronEt[ie]+analysisTree.electron_r03_sumPhotonEt[ie]-eA*rhoNeutral);
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
            
            //	    if (extraelec_veto) continue;
            //	    if (extramuon_veto) continue;
            
            //      cout << "dilepton_veto : " << dilepton_veto
            //	   << "   extraelec_veto : " << extraelec_veto
            //	   << "   extramuon_veto : " << extramuon_veto << endl;
            
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
            d0err_2 = analysisTree.muon_dxyerr[muonIndex];
            d0sig_2 = analysisTree.muon_dxy[muonIndex]/analysisTree.muon_dxyerr[muonIndex];
            d0sig_2_abs = fabs(analysisTree.muon_dxy[muonIndex]/analysisTree.muon_dxyerr[muonIndex]);
            dZ_2 = analysisTree.muon_dz[muonIndex];
            dZerr_2 = analysisTree.muon_dzerr[muonIndex];
            dZsig_2 = analysisTree.muon_dz[muonIndex]/analysisTree.muon_dzerr[muonIndex];
            dZsig_2_abs = fabs(analysisTree.muon_dz[muonIndex]/analysisTree.muon_dzerr[muonIndex]);
            iso_2 = isoMuMin;
            m_2 =classic_svFit::muonMass;
            
            float eleScale = eleScaleBarrel;
            if (fabs(analysisTree.electron_eta[electronIndex])>1.479) eleScale = eleScaleEndcap;

            // filling electron variables
            pt_1 = analysisTree.electron_pt[electronIndex];
            pt_1_escaleUp = (1+eleScale)*pt_1;
            pt_1_escaleDown = (1-eleScale)*pt_1;
            
            eta_1 = analysisTree.electron_eta[electronIndex];
            phi_1 = analysisTree.electron_phi[electronIndex];
            q_1 = -1;
            if (analysisTree.electron_charge[electronIndex]>0)
               q_1 = 1;
            mva_1 = analysisTree.electron_mva_wp80_Iso_Fall17_v1[electronIndex];
            d0_1 = analysisTree.electron_dxy[electronIndex];
            d0err_1 = analysisTree.electron_dxyerr[electronIndex];
            d0sig_1 = analysisTree.electron_dxy[electronIndex]/analysisTree.electron_dxyerr[electronIndex];
            d0sig_1_abs = fabs(analysisTree.electron_dxy[electronIndex]/analysisTree.electron_dxyerr[electronIndex]);
            dZ_1 = analysisTree.electron_dz[electronIndex];
            dZerr_1 = analysisTree.electron_dzerr[electronIndex];
            dZsig_1 = analysisTree.electron_dz[electronIndex]/analysisTree.electron_dzerr[electronIndex];
            dZsig_1_abs = fabs(analysisTree.electron_dz[electronIndex]/analysisTree.electron_dzerr[electronIndex]);
            d0_diff = analysisTree.muon_dxy[muonIndex] - analysisTree.electron_dxy[electronIndex];
            dZ_diff = analysisTree.muon_dz[muonIndex] - analysisTree.electron_dz[electronIndex];
            
            iso_1 = isoEleMin;
            m_1 = classic_svFit::electronMass;
            
            
            isoweight_1 = 1;
            isoweight_2 = 1;
            trigweight_1 = 1;
            trigweight_2 = 1;
            trigweight = 1;
            effweight = 1;
            
            //isoweight also includes tracking SF
            if (!isData || isEmbedded) {

             correctionWS->var("e_pt")->setVal(pt_1);
             correctionWS->var("e_eta")->setVal(eta_1);
             correctionWS->var("e_iso")->setVal(iso_1);
             correctionWS->var("m_pt")->setVal(pt_2);
             correctionWS->var("m_eta")->setVal(eta_2);
             correctionWS->var("m_iso")->setVal(iso_2);
               
               // scale factors
               if (applyWSCorr) {
                  if (isEmbedded) {
                     isoweight_1 = correctionWS->function("e_idiso_binned_embed_ratio")->getVal();
                     isoweight_2 = correctionWS->function("m_looseiso_binned_embed_ratio")->getVal() * correctionWS->function("m_id_embed_ratio")->getVal();
                  }
                  else {
                     isoweight_1 = correctionWS->function("e_idiso_binned_ratio")->getVal();
                     isoweight_2 = correctionWS->function("m_idiso_binned_ratio")->getVal();
                  }
               }
               else {
                  isoweight_1 = (float)SF_electronIdIso->get_ScaleFactor(double(pt_1),double(eta_1));
                  isoweight_2 = (float)SF_muonIdIso->get_ScaleFactor(double(pt_2),double(eta_2));
               }

               idweight_1 = correctionWS->function("e_trk_ratio")->getVal();
               idweight_2 = correctionWS->function("m_trk_ratio")->getVal();

               isoweight_1 *= idweight_1;
               isoweight_2 *= idweight_2;

	       //	       cout << "isoweight_1 = " << isoweight_1
	       //		    << "   isoweight_2 = " << isoweight_2 << endl;
	       //
	       float Ele23EffData = 1;
	       float Ele12EffData = 1;
	       float Mu23EffData = 1;
	       float Mu8EffData = 1;

	       float Ele23EffMC = 1;
	       float Ele12EffMC = 1;
	       float Mu23EffMC = 1;
	       float Mu8EffMC = 1;

	       if (applyWSCorr) {
             Ele23EffData = correctionWS->function("e_trg_binned_23_data")->getVal();
             Ele12EffData = correctionWS->function("e_trg_binned_12_data")->getVal(); 
             Mu23EffData  = correctionWS->function("m_trg_binned_23_data")->getVal();
             Mu8EffData   = correctionWS->function("m_trg_binned_8_data")->getVal();
             if (isEmbedded) {
                Ele23EffMC = correctionWS->function("e_trg_binned_23_embed")->getVal();
                Ele12EffMC = correctionWS->function("e_trg_binned_12_embed")->getVal(); 
                Mu23EffMC  = correctionWS->function("m_trg_binned_23_embed")->getVal();
                Mu8EffMC   = correctionWS->function("m_trg_binned_8_embed" )->getVal();
             }
             else {
                Ele23EffMC = correctionWS->function("e_trg_binned_23_mc")->getVal();
                Ele12EffMC = correctionWS->function("e_trg_binned_12_mc")->getVal();
                Mu23EffMC  = correctionWS->function("m_trg_binned_23_mc")->getVal();
                Mu8EffMC   = correctionWS->function("m_trg_binned_8_mc" )->getVal();
             }
	       }
	       else {
             Ele23EffData = (float)SF_electron23->get_EfficiencyData(double(pt_1),double(eta_1));
             Ele12EffData = (float)SF_electron12->get_EfficiencyData(double(pt_1),double(eta_1));
             Mu23EffData = (float)SF_muon23->get_EfficiencyData(double(pt_2),double(eta_2));
             Mu8EffData = (float)SF_muon8->get_EfficiencyData(double(pt_2),double(eta_2));
             Ele23EffMC   = (float)SF_electron23->get_EfficiencyMC(double(pt_1),double(eta_1));
             Ele12EffMC   = (float)SF_electron12->get_EfficiencyMC(double(pt_1),double(eta_1));
             Mu23EffMC   = (float)SF_muon23->get_EfficiencyMC(double(pt_2),double(eta_2));
             Mu8EffMC   = (float)SF_muon8->get_EfficiencyMC(double(pt_2),double(eta_2));
	       }

	       float trigWeightData = Mu23EffData*Ele12EffData + Mu8EffData*Ele23EffData - Mu23EffData*Ele23EffData;
               
	       if (applyTriggerMatch && !isData) {
		 float trigWeightMC   = Mu23EffMC*Ele12EffMC     + Mu8EffMC*Ele23EffMC     - Mu23EffMC*Ele23EffMC;
               
		 /*
		 if (isMuon23matched && isElectron12matched) {
		 trigweight_1 = (float)SF_electron12->get_ScaleFactor(double(pt_1),double(eta_1));
		   trigweight_2 = (float)SF_muon23->get_ScaleFactor(double(pt_2),double(eta_2));
		 }
		 else if (isMuon8matched && isElectron23matched) {
		   trigweight_1 = (float)SF_electron23->get_ScaleFactor(double(pt_1),double(eta_1));
		   trigweight_2 = (float)SF_muon8->get_ScaleFactor(double(pt_2),double(eta_2));
		 }
		 */
		 if (trigWeightMC>1e-6)
		   trigweight = trigWeightData / trigWeightMC;
               
	       }
	       else {
		 trigweight = trigWeightData;
		 //    if (isMuon23matched && isElectron12matched) {
		 //       trigweight_1 = (float)SF_electron12->get_EfficiencyData(double(pt_1),double(eta_1));
		 //       trigweight_2 = (float)SF_muon23->get_EfficiencyData(double(pt_2),double(eta_2));
		 //    }
		 //    else if (isMuon8matched && isElectron23matched) {
		 //       trigweight_1 = (float)SF_electron23->get_EfficiencyData(double(pt_1),double(eta_1));
		 //       trigweight_2 = (float)SF_muon8->get_EfficiencyData(double(pt_2),double(eta_2));
	       }
	       
	       //	       std::cout << "trigweight = " << trigweight << std::endl;
            }
	    effweight = trigweight*isoweight_1*isoweight_2;

            
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
            
            //store variables for visible mass
            m_vis = dileptonLV.M();
            
            m_vis_muUp    = (muonUpLV+electronLV).M();
            m_vis_muDown  = (muonDownLV+electronLV).M();
            m_vis_escaleUp     = (muonLV+electronUpLV).M();
            m_vis_escaleDown   = (muonLV+electronDownLV).M();
            m_vis_jesUp   = m_vis;
            m_vis_jesDown = m_vis;
            m_vis_resoUp    = m_vis;
            m_vis_resoDown  = m_vis;
            
            // std::cout << "m_vis        = " << m_vis << std::endl;
            // std::cout << "m_vis_muUp   = " << m_vis_muUp << std::endl;
            // std::cout << "m_vis_muDown = " << m_vis_muDown << std::endl;
            // std::cout << "m_vis_escaleUp    = " << m_vis_escaleUp << std::endl;
            // std::cout << "m_vis_escaleDown  = " << m_vis_escaleDown << std::endl;
            // std::cout << std::endl;
            
            //            pt_tt = dileptonLV.Pt();
            
            //delta R and delta Phi, muon and electron
            dphi_tt = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
                                 electronLV.Px(),electronLV.Py());
            
            dr_tt = deltaR(muonLV.Eta(),muonLV.Phi(),
                           electronLV.Eta(),electronLV.Phi());
            
            // qcd scale factor
            // no dzeta cut
            qcdweight     = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            qcdweightup   = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            qcdweightdown = qcdWeight.getWeight(pt_1,pt_2,dr_tt);
            // dzeta cut
            qcdweight_nodzeta     = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
            qcdweightup_nodzeta   = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
            qcdweightdown_nodzeta = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
            
            qcdweight     = 1.0;
            qcdweightup   = 1.0;
            qcdweightdown = 1.0;
            // dzeta cut
            qcdweight_nodzeta     = 1.0;
            qcdweightup_nodzeta   = 1.0;
            qcdweightdown_nodzeta = 1.0;


            //      if (os<0.5) {
            //	printf("QCD weights  : pt_1 = %6.1f ; pt_2 = %6.1f ; dr_tt = %4.2f\n",pt_1,pt_2,dr_tt);
            //	printf("2016         : central = %4.2f ; up = %4.2f ; down = %4.2f\n",qcdweight,qcdweightup,qcdweightdown);
            //	printf("2015         : central = %4.2f ; up = %4.2f ; down = %4.2f\n",qcdweight_nodzeta,qcdweightup_nodzeta,qcdweightdown_nodzeta);
            //	std::cout << std::endl;
            //      }
            
            //      std::cout << "Jets" << std::endl;
            
            // counting jets
            vector<unsigned int> jets; jets.clear();
            vector<unsigned int> jetsUp; jetsUp.clear();
            vector<unsigned int> jetsDown; jetsDown.clear();
            vector<unsigned int> jetspt20; jetspt20.clear();
            vector<unsigned int> bjets; bjets.clear();
            vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
            vector<unsigned int> bjetsRaw; bjetsRaw.clear();
            
            int indexLeadingJet = -1;
            float ptLeadingJet = -1;
            
            int indexSubLeadingJet = -1;
            float ptSubLeadingJet = -1;
            
            int indexLeadingBJet = -1;
            float ptLeadingBJet = -1;
          
            //impose cuts on jets
            for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
               
               double absJetEta = fabs(analysisTree.pfjet_eta[jet]);
               float jetEta = analysisTree.pfjet_eta[jet];
               if (absJetEta>jetEtaCut) continue;
               
               float jetPt = analysisTree.pfjet_pt[jet];
               float jetPtDown = analysisTree.pfjet_pt[jet]*(1.0-analysisTree.pfjet_jecUncertainty[jet]);
               float jetPtUp   = analysisTree.pfjet_pt[jet]*(1.0+analysisTree.pfjet_jecUncertainty[jet]);
               //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
               //if (jetPtDown<jetPtLowCut) continue;
               if (jetPtDown<jetPtLowCut) continue;

               bool isPFJetId = tightJetiD_2017(analysisTree,int(jet));
               if (!isPFJetId) continue;
         
               //check distance of jet to leptons
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
               //b-jet tagged?
                if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance
                   
                   bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                   bool taggedRaw = tagged;
                   
                   if (!isData) {
                      int flavor = abs(analysisTree.pfjet_flavour[jet]);
                      
                      double jet_scalefactor = 1;
                      double JetPtForBTag = jetPt;
                      double tageff = 1;
                      //store b-tagging eff
                      if (flavor==5) {
                         if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                         if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                         jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                         tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                      }
                        else if (flavor==4) {
                           if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                           if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                           jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                           tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                        }
                        else {
                           if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                           if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                           jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
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
                   
                   if (taggedRaw)
                      bjetsRaw.push_back(jet);
                   //store leading b jet
                   if (tagged) {
                      bjets.push_back(jet);
                      if (jetPt>ptLeadingBJet) {
                         ptLeadingBJet = jetPt;
                         indexLeadingBJet = jet;
                      }
                   }
                   
                }
                //jets with down and up variations of JEC passing pT cuts?
                if (jetPtUp>jetPtHighCut)
                   jetsUp.push_back(jet);
                
                if (jetPtDown>jetPtHighCut)
                   jetsDown.push_back(jet);

                if (jetPt>jetPtHighCut)
                   jets.push_back(jet);
                //also store subleading jets
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
            njets_jesUp = jetsUp.size();
            njets_jesDown = jetsDown.size();
            
            int njetsMax = njets;
            //check if number of jets changes when uncertainties are applied, if yes store this as njetsmax
            if (!isData) {
               jecUncertainties->runOnEvent(analysisTree,eta_1,phi_1,eta_2,phi_2);
               
               iUncert = 0;
               for (auto const& Name : uncertNames) {
                  njetsUncUp[iUncert] = jecUncertainties->getNJets(Name,true);
                  njetsUncDown[iUncert] = jecUncertainties->getNJets(Name,false);
                  mjjUncUp[iUncert] = jecUncertainties->getMjj(Name,true);
                  mjjUncDown[iUncert] = jecUncertainties->getMjj(Name,false);
                  if (njetsUncUp[iUncert]>njetsMax) njetsMax = njetsUncUp[iUncert];
                  if (njetsUncDown[iUncert]>njetsMax) njetsMax = njetsUncDown[iUncert];
                  iUncert++;
               }
            }


            int njetsCheckup = jecUncertainties->getNJets();
            
            njetspt20 = jetspt20.size();
            nbtag = bjets.size();
            nbtag_noSF = bjetsRaw.size();
            
            if (!isData) {
               int nnbtag = nbtag_noSF;
               btag0weight = 0;
               btag0weight_Up = 0;
               btag0weight_Down = 0;
            
               btag1weight = 0;
               btag1weight_Up = 0;
               btag1weight_Down = 0;

               btag2weight = 0;
               btag2weight_Up = 0;
               btag2weight_Down = 0;


               //if less than two btags sttore their pt and flavor, btag weight
               if (nnbtag<=2) {
                  
                  double b1Pt = 1;
                  double b2Pt = 1;
                  int b1Flav = 0;
                  int b2Flav = 0;
                  if (nnbtag>=1) {
                     int b1index = bjetsRaw.at(0);
                     b1Pt = analysisTree.pfjet_pt[b1index];
                     b1Flav = analysisTree.pfjet_flavour[b1index];
                  }
                  if (nnbtag==2) {
                     int b2index = bjetsRaw.at(1);
                     b2Pt = analysisTree.pfjet_pt[b2index];
                     b2Flav = analysisTree.pfjet_flavour[b2index];
                  }
                  
                  btag0weight = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,0,0));
                  btag0weight_Up = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,1,0));
                  btag0weight_Down = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,-1,0));

                  btag1weight = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,0,1));
                  btag1weight_Up = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,1,1));
                  btag1weight_Down = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,-1,1));

                  btag2weight = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,0,2));
                  btag2weight_Up = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,1,2));
                  btag2weight_Down = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,-1,2));
                  
               }
	      
               //	      cout << "nbtag(raw) = " << nnbtag << "  weight(central) = " << btag0weight 
               //		   << "   weight(up) = " << btag0weight_Up  
               //		   << "   weight(down) = " << btag0weight_Down << endl;
            }
            
            bpt = -10;
            beta = -10;
            bphi = -10;
            
            //store pt and eta of the leading bjet
            if (indexLeadingBJet>=0) {
               bpt = analysisTree.pfjet_pt[indexLeadingBJet];
               beta = analysisTree.pfjet_eta[indexLeadingBJet];
               bphi = analysisTree.pfjet_phi[indexLeadingBJet];
            }
            
            jpt_1 = -10;
            jpt_1_jesUp = -10;
            jpt_1_jesDown = -10;
            
            jeta_1 = -10;
            jphi_1 = -10;
            jptraw_1 = -10;
            jptunc_1 = -10;
            jmva_1 = -10;
            jlrm_1 = -10;
            jctm_1 = -10;
            
            if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
                cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
            
            //store pt and eta of the leading jet
            if (indexLeadingJet>=0) {
               jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
               jpt_1_jesUp = analysisTree.pfjet_pt[indexLeadingJet]*(1+analysisTree.pfjet_jecUncertainty[indexLeadingJet]);
               jpt_1_jesDown = analysisTree.pfjet_pt[indexLeadingJet]*(1-analysisTree.pfjet_jecUncertainty[indexLeadingJet]);
               jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
               jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
               jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
               jmva_1 = -1 ;//analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
               //		cout << "Leading jet pt = " << jpt_1 << "   eta = " << jeta_1 << endl;
               //		for (auto const& Name : uncertNames) {
               //		  cout << "    " << Name << " : " << jecUncertainties->getUncertainty(Name,jpt_1,jeta_1) << endl;
               //		}
            }
            //store pt and eta of the subleading jet
            jpt_2 = -10;
            jpt_2_jesUp = -10;
            jpt_2_jesDown = -10;
            
            jeta_2 = -10;
            jphi_2 = -10;
            jptraw_2 = -10;
            jptunc_2 = -10;
            jmva_2 = -10;
            jlrm_2 = -10;
            jctm_2 = -10;
            
            if (indexSubLeadingJet>=0) {
               jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
               jpt_2_jesUp = analysisTree.pfjet_pt[indexSubLeadingJet]*(1+analysisTree.pfjet_jecUncertainty[indexSubLeadingJet]);
               jpt_2_jesDown = analysisTree.pfjet_pt[indexSubLeadingJet]*(1-analysisTree.pfjet_jecUncertainty[indexSubLeadingJet]);
               jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
               jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
               jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];
               jmva_2 = -1;//analysisTree.pfjet_pu_jet_full_mva[indexSubLeadingJet];
            }
            
            mjj =  -10;
            mjj_jesUp = -10;
            mjj_jesDown = -10;
            jdeta =  -10;
            njetingap = 0;
            if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
               
               float unc1Up   = 1 + analysisTree.pfjet_jecUncertainty[indexLeadingJet]; 
               float unc1Down = 1 - analysisTree.pfjet_jecUncertainty[indexLeadingJet];
               
               float unc2Up   = 1 + analysisTree.pfjet_jecUncertainty[indexSubLeadingJet];
               float unc2Down = 1 - analysisTree.pfjet_jecUncertainty[indexSubLeadingJet];
               
               TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
                                                    analysisTree.pfjet_py[indexLeadingJet],
                                                    analysisTree.pfjet_pz[indexLeadingJet],
                                                    analysisTree.pfjet_e[indexLeadingJet]);
               
               TLorentzVector jet1Up; jet1Up.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet]*unc1Up,
                                                        analysisTree.pfjet_py[indexLeadingJet]*unc1Up,
                                                        analysisTree.pfjet_pz[indexLeadingJet]*unc1Up,
                                                        analysisTree.pfjet_e[indexLeadingJet]*unc1Up);
               
               TLorentzVector jet1Down; jet1Down.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet]*unc1Down,
                                                            analysisTree.pfjet_py[indexLeadingJet]*unc1Down,
                                                            analysisTree.pfjet_pz[indexLeadingJet]*unc1Down,
                                                            analysisTree.pfjet_e[indexLeadingJet]*unc1Down);
               
               TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
                                                    analysisTree.pfjet_py[indexSubLeadingJet],
                                                    analysisTree.pfjet_pz[indexSubLeadingJet],
                                                    analysisTree.pfjet_e[indexSubLeadingJet]);
               
               
               TLorentzVector jet2Up; jet2Up.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet]*unc2Up,
                                                        analysisTree.pfjet_py[indexSubLeadingJet]*unc2Up,
                                                        analysisTree.pfjet_pz[indexSubLeadingJet]*unc2Up,
                                                        analysisTree.pfjet_e[indexSubLeadingJet]*unc2Up);
               
               TLorentzVector jet2Down; jet2Down.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet]*unc2Down,
                                                            analysisTree.pfjet_py[indexSubLeadingJet]*unc2Down,
                                                            analysisTree.pfjet_pz[indexSubLeadingJet]*unc2Down,
                                                            analysisTree.pfjet_e[indexSubLeadingJet]*unc2Down);
               //store mjj, jdeta, etamax etamin
               mjj = (jet1+jet2).M();
               mjj_jesUp = (jet1Up+jet2Up).M();
               mjj_jesDown = (jet1Down+jet2Down).M();
               
               jdeta = fabs(analysisTree.pfjet_eta[indexLeadingJet]-
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
            float met_x = analysisTree.pfmetcorr_ex;
            float met_y = analysisTree.pfmetcorr_ey;
            //            if (!isData) {
            //                met_x = analysisTree.pfmet_ex;
            //                met_y = analysisTree.pfmet_ey;
            //            }
            met_x = analysisTree.pfmetcorr_ex;
            met_y = analysisTree.pfmetcorr_ey;
            float met_jesUp_x   = analysisTree.pfmetcorr_ex_JetEnUp;
            float met_jesUp_y   = analysisTree.pfmetcorr_ey_JetEnUp;
            float met_jesDown_x = analysisTree.pfmetcorr_ex_JetEnDown;
            float met_jesDown_y = analysisTree.pfmetcorr_ey_JetEnDown;
            float met_resoUp_x    = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            float met_resoUp_y    = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            float met_resoDown_x  = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            float met_resoDown_y  = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
            
            float metcorr_jesUp_x   = analysisTree.pfmetcorr_ex_JetEnUp;
            float metcorr_jesUp_y   = analysisTree.pfmetcorr_ey_JetEnUp;
            float metcorr_jesDown_x = analysisTree.pfmetcorr_ex_JetEnDown;
            float metcorr_jesDown_y = analysisTree.pfmetcorr_ey_JetEnDown;
            float metcorr_resoUp_x    = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            float metcorr_resoUp_y    = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            float metcorr_resoDown_x  = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            float metcorr_resoDown_y  = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
            
            met = TMath::Sqrt(met_x*met_x + met_y*met_y);
            metphi = TMath::ATan2(met_y,met_x);
            metcov00 = analysisTree.pfmetcorr_sigxx;
            metcov01 = analysisTree.pfmetcorr_sigxy;
            metcov10 = analysisTree.pfmetcorr_sigyx;
            metcov11 = analysisTree.pfmetcorr_sigyy;
        
            /*
              if(!isData)
              {
              metcov00 = analysisTree.pfmet_sigxx;
              metcov01 = analysisTree.pfmet_sigxy;
              metcov10 = analysisTree.pfmet_sigyx;
              metcov11 = analysisTree.pfmet_sigyy;
              }
            */
            /*
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
            */            
            //MVAMET =MET
            float mvamet_x = met_x;
            float mvamet_y = met_y;
            mvamet = met;
            mvametphi = metphi;
            /*
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
            */
            // Recoil corrections
            
            int njetsforrecoil = njets;
            if (isW) njetsforrecoil = njets + 1;  
            
            float mvamet_corr_x = mvamet_x;
            float mvamet_corr_y = mvamet_y;
            
            float mvamet_uncorr_x = mvamet_x;
            float mvamet_uncorr_y = mvamet_y;
            
            mvamet_uncorr = mvamet;
            mvametphi_uncorr = mvametphi;
            
            met_uncorr = met;
            metphi_uncorr = metphi;
            
            float pfmet_corr_x = met_x;
            float pfmet_corr_y = met_y;
            
            if ( (isW||isDY||isSignal) &&!isData && applyRecoilCorrections ) {
               //store corrected MET
               if (applySimpleRecoilCorrections) {
                  //                    recoilMvaMetCorrector.CorrectByMeanResolution(mvamet_x,mvamet_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,mvamet_corr_x,mvamet_corr_y);
                  recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                  recoilMetCorrector.CorrectByMeanResolution(metcorr_jesUp_x,  metcorr_jesUp_y,  bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_jesUp_x,  met_jesUp_y);
                  recoilMetCorrector.CorrectByMeanResolution(metcorr_jesDown_x,metcorr_jesDown_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_jesDown_x,met_jesDown_y);
                  recoilMetCorrector.CorrectByMeanResolution(metcorr_resoUp_x,   metcorr_resoUp_y,   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_resoUp_x,   met_resoUp_y);
                  recoilMetCorrector.CorrectByMeanResolution(metcorr_resoDown_x, metcorr_resoDown_y, bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_resoDown_x, met_resoDown_y);
                  
               }
               else {
                  //                    recoilMvaMetCorrector.Correct(mvamet_x,mvamet_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,mvamet_corr_x,mvamet_corr_y);
                  recoilMetCorrector.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
               }
            }

            //for DY and W+jets corrected MET is stored
            met_x = pfmet_corr_x;
            met_y = pfmet_corr_y;
            met = TMath::Sqrt(met_x*met_x+met_y*met_y);
            metphi = TMath::ATan2(met_y,met_x);
            
            //      printf("Before correction : mvamet_x = %7.2f   mvamet_y = %7.2f\n",mvamet_x,mvamet_y);
            //      printf("After correction  : mvamet_x = %7.2f   mvamet_y = %7.2f\n",mvamet_corr_x,mvamet_corr_y);
            //      std::cout << std::endl;
            
            mvamet_x = mvamet_corr_x;
            mvamet_y = mvamet_corr_y;
            
            mvamet = TMath::Sqrt(mvamet_x*mvamet_x+mvamet_y*mvamet_y); 
            mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
        
            float mvamet_jesUp_x   = mvamet_x;
            float mvamet_jesUp_y   = mvamet_y;
            float mvamet_jesDown_x = mvamet_x;
            float mvamet_jesDown_y = mvamet_y;
            float mvamet_resoUp_x    = mvamet_x;
            float mvamet_resoUp_y    = mvamet_y;
            float mvamet_resoDown_x  = mvamet_x;
            float mvamet_resoDown_y  = mvamet_y;

            mvamet_jesUp = TMath::Sqrt(mvamet_jesUp_x*mvamet_jesUp_x+
                                         mvamet_jesUp_y*mvamet_jesUp_y);
            mvametphi_jesUp = TMath::ATan2(mvamet_jesUp_y,mvamet_jesUp_x);
            
            mvamet_jesDown = TMath::Sqrt(mvamet_jesDown_x*mvamet_jesDown_x+
                                           mvamet_jesDown_y*mvamet_jesDown_y);
            mvametphi_jesDown = TMath::ATan2(mvamet_jesDown_y,mvamet_jesDown_x);
            
            mvamet_resoUp = TMath::Sqrt(mvamet_resoUp_x*mvamet_resoUp_x+
                                        mvamet_resoUp_y*mvamet_resoUp_y);
            mvametphi_resoUp = TMath::ATan2(mvamet_resoUp_y,mvamet_resoUp_x);
            
            mvamet_resoDown = TMath::Sqrt(mvamet_resoDown_x*mvamet_resoDown_x+
                                          mvamet_resoDown_y*mvamet_resoDown_y);
            mvametphi_resoDown = TMath::ATan2(mvamet_resoDown_y,mvamet_resoDown_x);
            
 
            met_jesUp = TMath::Sqrt(met_jesUp_x*met_jesUp_x+
				      met_jesUp_y*met_jesUp_y);
            metphi_jesUp = TMath::ATan2(met_jesUp_y,met_jesUp_x);
            
            met_jesDown = TMath::Sqrt(met_jesDown_x*met_jesDown_x+
					met_jesDown_y*met_jesDown_y);
            metphi_jesDown = TMath::ATan2(met_jesDown_y,met_jesDown_x);
            
            met_resoUp = TMath::Sqrt(met_resoUp_x*met_resoUp_x+
				     met_resoUp_y*met_resoUp_y);
            metphi_resoUp = TMath::ATan2(met_resoUp_y,met_resoUp_x);
            
            met_resoDown = TMath::Sqrt(met_resoDown_x*met_resoDown_x+
				       met_resoDown_y*met_resoDown_y);
            metphi_resoDown = TMath::ATan2(met_resoDown_y,met_resoDown_x);
            
            
            //      printf("Uncorrected    :  %7.2f  %7.2f\n",mvamet_uncorr_x,mvamet_uncorr_y);
            //      printf("Central        :  %7.2f  %7.2f\n",mvamet_x,mvamet_y);
            //      printf("ScaleUp        :  %7.2f  %7.2f\n",mvamet_jesUp_x,mvamet_jesUp_y);
            //      printf("ScaleDown      :  %7.2f  %7.2f\n",mvamet_jesDown_x,mvamet_jesDown_y);
            //      printf("ResoUp         :  %7.2f  %7.2f\n",mvamet_resoUp_x,mvamet_resoUp_y);
            //      printf("ResoDown       :  %7.2f  %7.2f\n",mvamet_resoDown_x,mvamet_resoDown_y);
            //      std::cout << std::endl;
           
            //store genmet
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

            float vectorVisX_escaleUp = muonLV.Px() + electronUpLV.Px();
            float vectorVisY_escaleUp = muonLV.Py() + electronUpLV.Py();

            float vectorVisX_escaleDown = muonLV.Px() + electronDownLV.Px();
            float vectorVisY_escaleDown = muonLV.Py() + electronDownLV.Py();
            
            pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
            pzetavis_escaleUp = vectorVisX_escaleUp*zetaX + vectorVisY_escaleUp*zetaY;
            pzetavis_escaleDown = vectorVisX_escaleDown*zetaX + vectorVisY_escaleDown*zetaY;
            
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

            computeDzeta(met_jesUp_x,met_jesUp_y,
                         zetaX,zetaY,pzetavis,pzetamiss_jesUp,dzeta_jesUp); // scaleUp
            computeDzeta(met_jesDown_x,met_jesDown_y,
                         zetaX,zetaY,pzetavis,pzetamiss_jesDown,dzeta_jesDown); // scaleDown
            computeDzeta(met_resoUp_x,met_resoUp_y,
                         zetaX,zetaY,pzetavis,pzetamiss_resoUp,dzeta_resoUp); // resoUp
            computeDzeta(met_resoDown_x,met_resoDown_y,
                         zetaX,zetaY,pzetavis,pzetamiss_resoDown,dzeta_resoDown); // resoDown
            computeDzeta(metx_escaleUp,mety_escaleUp,
                         zetaX,zetaY,pzetavis_escaleUp,pzetamiss_escaleUp,dzeta_escaleUp); // eUp
            computeDzeta(metx_escaleDown,mety_escaleDown,
                         zetaX,zetaY,pzetavis_escaleDown,pzetamiss_escaleDown,dzeta_escaleDown); // eDown
            
            // mvamet
            computeDzeta(mvamet_x,mvamet_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet,dzeta_mvamet); // corrected
            computeDzeta(mvamet_uncorr_x,mvamet_uncorr_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet_uncorr,dzeta_mvamet_uncorr); // uncorrected
            computeDzeta(mvamet_jesUp_x,mvamet_jesUp_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet_jesUp,dzeta_mvamet_jesUp); // scaleUp
            computeDzeta(mvamet_jesDown_x,mvamet_jesDown_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet_jesDown,dzeta_mvamet_jesDown); // scaleDown
            computeDzeta(mvamet_resoUp_x,mvamet_resoUp_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet_resoUp,dzeta_mvamet_resoUp); // resoUp
            computeDzeta(mvamet_resoDown_x,mvamet_resoDown_y,
                         zetaX,zetaY,pzetavis,pzetamiss_mvamet_resoDown,dzeta_mvamet_resoDown); // resoDown
            
            // genmet
            computeDzeta(analysisTree.genmet_ex,analysisTree.genmet_ey,
                         zetaX,zetaY,pzetavis,pzetamiss_genmet,dzeta_genmet);
            
	    //            float ETmis = TMath::Sqrt(mvamet_x*mvamet_x+mvamet_y*mvamet_y);
            TLorentzVector mvametLV; mvametLV.SetXYZT(mvamet_x,mvamet_y,0.,mvamet);
	    //            float PFETmis = TMath::Sqrt(met_x*met_x+met_y*met_y);
            TLorentzVector metLV; metLV.SetXYZT(met_x,met_y,0.,met);

            float met_escaleUp = TMath::Sqrt(metx_escaleUp*metx_escaleUp+mety_escaleUp*mety_escaleUp);
            float met_escaleDown = TMath::Sqrt(metx_escaleDown*metx_escaleDown+mety_escaleDown*mety_escaleDown);
            
            TLorentzVector metEleUpLV; metEleUpLV.SetXYZT(metx_escaleUp,mety_escaleUp,0.,met_escaleUp);
            TLorentzVector metEleDownLV; metEleDownLV.SetXYZT(metx_escaleDown,mety_escaleDown,0.,met_escaleDown);
            //calculate tarnsverse masses
            mt_1 = mT(electronLV,metLV);
            mt_2 = mT(muonLV,metLV);
            mtmax = TMath::Max(float(mt_1),float(mt_2));
            //transverse mass of muon + electron +met
            mCDF = (muonLV+electronLV+metLV).M();
            //transverse momentum of muon + electron +met
            pt_tt = (muonLV+electronLV+metLV).Pt();
            m_leptons = (muonLV+electronLV).M();
            
            
            TLorentzVector metResoUpLV; metResoUpLV.SetXYZT(met_resoUp_x,met_resoUp_y,0.,met_resoUp);
            TLorentzVector metResoDownLV; metResoDownLV.SetXYZT(met_resoDown_x,met_resoDown_y,0.,met_resoDown);
            TLorentzVector metScaleUpLV; metScaleUpLV.SetXYZT(met_jesUp_x,met_jesUp_y,0.,met_jesUp);
            TLorentzVector metScaleDownLV; metScaleDownLV.SetXYZT(met_jesDown_x,met_jesDown_y,0.,met_jesDown);
            // computing total transverse mass
            mTtot = totalTransverseMass        ( muonLV ,     electronLV , metLV);
            
            mTtot_muUp   =  totalTransverseMass( muonUpLV ,   electronLV , metLV);
            mTtot_muDown =  totalTransverseMass( muonDownLV , electronLV , metLV);
            
            mTtot_escaleUp   =  totalTransverseMass ( muonLV , electronUpLV   , metLV);
            mTtot_escaleDown =  totalTransverseMass ( muonLV , electronDownLV , metLV);
            
            mTtot_jesUp   =  totalTransverseMass ( muonLV , electronLV , metScaleUpLV);
            mTtot_jesDown =  totalTransverseMass ( muonLV , electronLV , metScaleDownLV);
           
            mTtot_resoUp   =  totalTransverseMass ( muonLV , electronLV , metResoUpLV);
            mTtot_resoDown =  totalTransverseMass ( muonLV , electronLV , metResoDownLV);
            
            mTemu   = mT(electronLV,muonLV);
            mTemet  = mT(electronLV,metLV);
            mTmumet = mT(muonLV,metLV);
            
            dphi_mumet = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
                                    metLV.Px(),metLV.Py());
            dphi_emet  = dPhiFrom2P(electronLV.Px(),electronLV.Py(),
                                    metLV.Px(),metLV.Py());
            
            bdt     = reader->EvaluateMVA("BDT");   
            bdt_ggh = readerGGH->EvaluateMVA("BDT");
            bdt_bbh = readerBBH->EvaluateMVA("BDT");
            
            mTdileptonMET = mT(dileptonLV,metLV);

            mTdileptonMET_jesUp = mT(dileptonLV,metScaleUpLV);
            mTdileptonMET_jesDown = mT(dileptonLV,metScaleDownLV);
            mTdileptonMET_resoUp = mT(dileptonLV,metResoUpLV);
            mTdileptonMET_resoDown = mT(dileptonLV,metResoDownLV);
            
            mTdileptonMET_muUp = mT(muonUpLV+electronLV,metLV);
            mTdileptonMET_muDown = mT(muonDownLV+electronLV,metLV);
            mTdileptonMET_escaleUp = mT(muonLV+electronUpLV,metEleUpLV);
            mTdileptonMET_escaleDown = mT(muonLV+electronDownLV,metEleDownLV);
            
            //	    std::cout << "mT(ll,MET) = " << mTdileptonMET << "  mTtot = " << mTtot << std::endl;
            //      std::cout << "BDT       = " << bdt << std::endl;
            //      std::cout << "BDT (bbH) = " << bdt_bbh << std::endl;
            //      std::cout << "BDT (ggH) = " << bdt_ggh << std::endl;
            
            m_sv           = -10;
            m_sv_muUp      = -10;
            m_sv_muDown    = -10;
            m_sv_escaleUp       = -10;
            m_sv_escaleDown     = -10;
            m_sv_jesUp   = -10;
            m_sv_jesDown = -10;
            m_sv_resoUp    = -10;
            m_sv_resoDown  = -10;
            
            mt_sv           = -10;
            mt_sv_muUp      = -10;
            mt_sv_muDown    = -10;
            mt_sv_escaleUp       = -10;
            mt_sv_escaleDown     = -10;
            mt_sv_jesUp   = -10;
            mt_sv_jesDown = -10;
            mt_sv_resoUp    = -10;
            mt_sv_resoDown  = -10;
            
            pt_sv   = -10;
            eta_sv  = -10;
            phi_sv  = -10;
            
            pt_sv_escaleUp  = -10;
            eta_sv_escaleUp = -10;
            phi_sv_escaleUp = -10;
            
            pt_sv_escaleDown  = -10;
            eta_sv_escaleDown = -10;
            phi_sv_escaleDown = -10;

            //calculate SV mass only for certain events
            //if (computeSVFitMass && dzeta>-40 && iso_1<0.5 && iso_2<0.5 && njetsMax>0) {
                if (computeSVFitMass) {
               //                if (mvaMetFound) {
               // covariance matrix MET
                    TMatrixD covMET(2, 2);
                    covMET[0][0] =  metcov00;
                    covMET[1][0] =  metcov10;
                    covMET[0][1] =  metcov01;
                    covMET[1][1] =  metcov11;

		    //		    std::cout << "Eta = " << eta_1 << "   eleScale = " << eleScale << std::endl;

                    
                    // mva met
                    // define electron 4-vector
                    classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, 
                                                                pt_1, eta_1, phi_1, 0.51100e-3); 
                    classic_svFit::MeasuredTauLepton svFitEleUp(classic_svFit::MeasuredTauLepton::kTauToElecDecay, 
                                                                  (1.0+eleScale)*pt_1, eta_1, phi_1, 0.51100e-3); 
                    classic_svFit::MeasuredTauLepton svFitEleDown(classic_svFit::MeasuredTauLepton::kTauToElecDecay, 
                                                                    (1.0-eleScale)*pt_1, eta_1, phi_1, 0.51100e-3); 
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
                    eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
                    phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
                                        
                    float px_sv = pt_sv*TMath::Cos(phi_sv);
                    float py_sv = pt_sv*TMath::Sin(phi_sv);
                    
                    float msvmet_ex = px_sv - dileptonLV.Px();
                    float msvmet_ey = py_sv - dileptonLV.Py();
                    
                    msvmet = TMath::Sqrt(msvmet_ex*msvmet_ex+msvmet_ey*msvmet_ey);
                    msvmetphi = TMath::ATan2(msvmet_ey,msvmet_ex);
                    
                    //	  SVfitStandaloneAlgorithm algoGen = SVFitMassComputation(svFitEle, svFitMu,
                    //								  genmet_ex, genmet_ey, 
                    //								  covMVAMET, inputFile_visPtResolution);
                    
                    //	  msv_gen  = algoGen.getMass(); // return value of svfit mass is in units of GeV
                    //	  mtsv_gen = algoGen.transverseMass(); // return value of transverse svfit mass is in units of GeV
                    //	  if ( !algoGen.isValidSolution() ) 
                    //	    std::cout << "sorry -- status of NLL is not valid [" << algoGen.isValidSolution() << "]" << std::endl;
                    
                    //		    std::cout << "msv = " << m_sv << "   msv_gen = " << msv_gen << std::endl;
                    //	  std::cout << std::endl;
                    
                    //	  pt_sv_gen = algoGen.pt();
                    //	  eta_sv_gen = algoGen.eta();
                    //	  phi_sv_gen = algo.phi();
                    
                    //	  float px_sv_gen = pt_sv_gen*TMath::Cos(phi_sv_gen);
                    //	  float py_sv_gen = pt_sv_gen*TMath::Sin(phi_sv_gen);
                    
                    //	  nuPx_msv = px_sv_gen - dileptonLV.Px();
                    //	  nuPy_msv = py_sv_gen - dileptonLV.Px();
                    //	  nuPz_msv = nuPz;
                    //	  nuPt_msv = TMath::Sqrt(nuPx_msv*nuPx_msv+nuPy_msv*nuPy_msv);
                    //	  nuPhi_msv = TMath::ATan2(nuPy_msv,nuPx_msv);
                    
                    m_sv_jesUp   = m_sv;
                    m_sv_jesDown = m_sv;
                    m_sv_resoUp    = m_sv;
                    m_sv_resoDown  = m_sv;
                    m_sv_escaleUp       = m_sv;
                    m_sv_escaleDown     = m_sv;
                    m_sv_muUp      = m_sv;
                    m_sv_muDown    = m_sv;
                    
                    mt_sv_jesUp   = mt_sv;
                    mt_sv_jesDown = mt_sv;
                    mt_sv_resoUp    = mt_sv;
                    mt_sv_resoDown  = mt_sv;
                    mt_sv_escaleUp       = mt_sv;
                    mt_sv_escaleDown     = mt_sv;
                    mt_sv_muUp      = mt_sv;
                    mt_sv_muDown    = mt_sv;
                
                    bool applyMSVvariations = true;
                   
                    if (!isData && applyMSVvariations) { 
                        

		      /*		      		      
                                    std::cout << "Event = " << analysisTree.event_nr << std::endl;
                                    std::cout << "MetX = " << met_x << "   MetY = " << met_y << std::endl;
                                    std::cout << "MetX_escaleUp = " << metx_escaleUp << "   MetY_escaleUp = " << mety_escaleUp << std::endl;
                                    std::cout << "MetX_ = " << metx_escaleDown << "   MetY_escaleDown = " << mety_escaleDown << std::endl;
                                    
                                    std::cout << "eta1 = " << eta_1 << std::endl;
                                    std::cout << "pt1 = " << pt_1 << "  pt1_escaleUp = " << (1.0+eleScale)*pt_1 << "  pt1_escaleDown = " << (1.0-eleScale)*pt_1 << std::endl;
                                    std::cout << "cov00 = " << metcov00 << "  cov01 = " << metcov01 << "   cov11 = " << metcov11 << std::endl; 
		      */

                       ClassicSVfit algo_escaleUp = SVFitMassComputation(svFitEleUp, svFitMu,
                                                        metx_escaleUp, mety_escaleUp,
                                                        covMET, inputFile_visPtResolution);
                       ClassicSVfit algo_escaleDown = SVFitMassComputation(svFitEleDown, svFitMu,
                                                                                  metx_escaleDown, mety_escaleDown,
                                                                                   covMET, inputFile_visPtResolution);
                        
                        
                       m_sv_escaleUp    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleUp.getHistogramAdapter())->getMass();
                       mt_sv_escaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleUp.getHistogramAdapter())->getTransverseMass();
                       
                       m_sv_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleDown.getHistogramAdapter())->getMass();
                       mt_sv_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleDown.getHistogramAdapter())->getTransverseMass();

                       pt_sv_escaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleUp.getHistogramAdapter())->getPt(); 
                       eta_sv_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleUp.getHistogramAdapter())->getEta(); 
                       phi_sv_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleUp.getHistogramAdapter())->getPhi(); 

                       pt_sv_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleDown.getHistogramAdapter())->getPt(); 
                       eta_sv_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleDown.getHistogramAdapter())->getEta(); 
                       phi_sv_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo_escaleDown.getHistogramAdapter())->getPhi();  

                        /*
                          std::cout << "msv = " << m_sv 
                          << "   msv_escaleUp = " << m_sv_escaleUp 
                          << "   msv_escaleDown = " << m_sv_escaleDown << std::endl;
                          std::cout << std::endl;
                        */
                        
                        /*                        
                                                  if (isDY || isW) {
                                                  SVfitStandaloneAlgorithm algo_jesUp = SVFitMassComputation(svFitEle, svFitMu,
                                                  mvamet_jesUp_x, mvamet_jesUp_y,
                                                  covMET, inputFile_visPtResolution);
                                                  
                                                  SVfitStandaloneAlgorithm algo_jesDown = SVFitMassComputation(svFitEle, svFitMu,
                                                  mvamet_jesDown_x, mvamet_jesDown_y,
                                                  covMET, inputFile_visPtResolution);
                                                  
                                                  // SVfitStandaloneAlgorithm algo_resoUp = SVFitMassComputation(svFitEle, svFitMu,
                                                  // 								mvamet_resoUp_x, mvamet_resoUp_y,
                                                  // 								covMVAMET, inputFile_visPtResolution);
                                                  
                                                  // SVfitStandaloneAlgorithm algo_resoDown = SVFitMassComputation(svFitEle, svFitMu,
                                                  // 								  mvamet_resoDown_x, mvamet_resoDown_y,
                                                  // 								  covMVAMET, inputFile_visPtResolution);
                                                  
                                                  m_sv_jesUp    = algo_jesUp.getMass();
                                                  mt_sv_jesUp   = algo_jesUp.transverseMass();
                                                  m_sv_jesDown  = algo_jesDown.getMass();
                                                  mt_sv_jesDown = algo_jesDown.transverseMass();
                                                  
                                                  // m_sv_resoUp     = algo_resoUp.getMass();
                                                  // mt_sv_resoUp    = algo_resoUp.transverseMass();
                                                  // m_sv_resoDown   = algo_resoDown.getMass();
                                                  // mt_sv_resoDown  = algo_resoDown.transverseMass();
                                                  
                                                  }
                        */
                        //                    }
                        //	  std::cout << "eta(e) = " << eta_1 << "   escale = " << eleScale << std::endl;
                        //	  std::cout << "msv = " << m_sv << std::endl;
                        //	  std::cout << "msv(eES)      : up = " << m_sv_escaleUp << "   down = " << m_sv_escaleDown << std::endl;
                        //	  std::cout << "msv(metScale) : up = " << m_sv_jesUp << "   down = " << m_sv_jesDown << std::endl;
                        //	  std::cout << "mtsv = " << mt_sv << std::endl;
                        //	  std::cout << "mtsv(eES)      : up = " << mt_sv_escaleUp << "   down = " << mt_sv_escaleDown << std::endl;
                        //	  std::cout << "mtsv(metScale) : up = " << mt_sv_jesUp << "   down = " << mt_sv_jesDown << std::endl;
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
            unsigned int gen_1 = 0;
            bool bgen_1 = false;
            unsigned int gen_2 = 0;
            bool bgen_2 = false;
            
            //store genElectron, genMuon, genmet
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
               
	       /*
               double weightE = 1;   //weights for jet electron fake rate
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
                  //    weightE = 1;
                  //    weightEUp = 1;
                  //    weightEDown = 1;
                  // }
                  //else {
                     weightE     = a_jetEle + b_jetEle*ptJet;
                     weightEUp   = a_jetEleUp + b_jetEleUp*ptJet;
                     weightEDown = a_jetEleDown + b_jetEleDown*ptJet;
                     //}
                  fakeweight *= weightE;
               }
               else {
                  weightE = isoweight_1;
                  weightEUp = isoweight_1;
                  weightEDown = isoweight_1;
               }
               
               double weightMu = 1;   //weights for jet muon fake rate
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
                  // if (isTOP) {
                  //    weightMu = 1;
                  //    weightMuUp = 1;
                  //    weightMuDown = 1;
                  // }
                  //else {
                     weightMu     = a_jetMu     + b_jetMu*ptJet;
                     weightMuUp   = a_jetMuUp   + b_jetMuUp*ptJet;
                     weightMuDown = a_jetMuDown + b_jetMuDown*ptJet;
                     // }
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
	       */
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

