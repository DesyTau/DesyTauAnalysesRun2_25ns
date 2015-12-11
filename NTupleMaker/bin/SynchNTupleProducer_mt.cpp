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

int binNumber(float x, int nbins, float * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double PtoPhi(double Px, double Py) {
  return TMath::ATan2(Py,Px);
}

double PtoPt(double Px, double Py) {
  return TMath::Sqrt(Px*Px+Py*Py);
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;

}
bool electronMvaIdTight(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.73) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.57) passed = true;
  }
  else {
    if (mva>0.05) passed = true;
  }

  return passed;

}
bool electronMvaIdLoose(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.35) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.20) passed = true;
  }
  else {
    if (mva>-0.52) passed = true;
  }

  return passed;

}
bool electronMvaIdWP80(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.253;
    else 
      passed = mva > 0.965;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > 0.081;
    else
      passed = mva > 0.917;
  }
  else {
    if (pt<10)
      passed = mva > -0.081;
    else
      passed = mva > 0.683;
  }

  return passed;

}
bool electronMvaIdWP90(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.483;
    else 
      passed = mva > 0.933;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > -0.267;
    else
      passed = mva > 0.825;
  }
  else {
    if (pt<10)
      passed = mva > -0.323;
    else
      passed = mva > 0.337;
  }

  return passed;

}


const float electronMass = 0.51100e-3;
const float muonMass = 0.10565837;
const float tauMass  = 1.7768;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");
  const bool isDY = cfg.get<bool>("IsDY");
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

  // check overlap
  const bool checkOverlap = cfg.get<bool>("CheckOverlap");
  
  TString IsoMuonLeg(isoMuonLeg);
  TString MuonTauMuonLeg(muonTauMuonLeg);
  TString MuonTauOverlap(muonTauOverlap);
  TString MuonTauTauLeg(muonTauTauLeg);
  TString BTagDiscriminator(bTagDiscriminator);

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
  Float_t         m_sv;
  Float_t         pt_sv;
  Float_t         eta_sv;
  Float_t         phi_sv;

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

  Float_t         met;
  Float_t         metphi;
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;

  Float_t         puppimet;
  Float_t         puppimetphi;
  Float_t         m_sv_puppi;
  Float_t         pt_sv_puppi;
  Float_t         eta_sv_puppi;
  Float_t         phi_sv_puppi;

  Float_t         mvamet;
  Float_t         mvametphi;
  Float_t         mvacov00;
  Float_t         mvacov01;
  Float_t         mvacov10;
  Float_t         mvacov11;

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

  tree->Branch("m_vis", &m_vis, "m_vis/F");

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

  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, "byCombinedIsolationDeltaBetaCorrRaw3Hits_1/F");
  tree->Branch("againstElectronLooseMVA5_1", &againstElectronLooseMVA5_1, "againstElectronLooseMVA5_1/F");
  tree->Branch("againstElectronMediumMVA5_1", &againstElectronMediumMVA5_1, "againstElectronMediumMVA5_1/F");
  tree->Branch("againstElectronTightMVA5_1", &againstElectronTightMVA5_1, "againstElectronTightMVA5_1/F");
  tree->Branch("againstElectronVLooseMVA5_1", &againstElectronVLooseMVA5_1, "againstElectronVLooseMVA5_1/F");
  tree->Branch("againstElectronVTightMVA5_1", &againstElectronVTightMVA5_1, "againstElectronVTightMVA5_1/F");
  tree->Branch("againstMuonLoose3_1", &againstMuonLoose3_1, "againstMuonLoose3_1/F");
  tree->Branch("againstMuonTight3_1", &againstMuonTight3_1, "againstMuonTight3_1/F");

  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
  tree->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2",&byLooseCombinedIsolationDeltaBetaCorr3Hits_2,"byLooseCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2",&byMediumCombinedIsolationDeltaBetaCorr3Hits_2,"byMediumCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2",&byTightCombinedIsolationDeltaBetaCorr3Hits_2,"byTightCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2, "againstElectronLooseMVA5_2/F");
  tree->Branch("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2, "againstElectronMediumMVA5_2/F");
  tree->Branch("againstElectronTightMVA5_2", &againstElectronTightMVA5_2, "againstElectronTightMVA5_2/F");
  tree->Branch("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2, "againstElectronVLooseMVA5_2/F");
  tree->Branch("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2, "againstElectronVTightMVA5_2/F");
  tree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/F");
  tree->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");

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
  tree->Branch("mvacov00", &mvacov00, "mvacov00/F");
  tree->Branch("mvacov01", &mvacov01, "mvacov01/F");
  tree->Branch("mvacov10", &mvacov10, "mvacov10/F");
  tree->Branch("mvacov11", &mvacov11, "mvacov11/F");

  tree->Branch("genmet",&genmet,"genmet/F");
  tree->Branch("genmetphi",&genmetphi,"genmetphi/F");

  tree->Branch("m_sv", &m_sv, "m_sv/F");
  tree->Branch("pt_sv", &pt_sv, "pt_sv/F");
  tree->Branch("eta_sv", &eta_sv, "eta_sv/F");
  tree->Branch("phi_sv", &phi_sv, "phi_sv/F");

  tree->Branch("m_sv_puppi", &m_sv_puppi, "m_sv_puppi/F");
  tree->Branch("pt_sv_puppi", &pt_sv_puppi, "pt_sv_puppi/F");
  tree->Branch("eta_sv_puppi", &eta_sv_puppi, "eta_sv_puppi/F");
  tree->Branch("phi_sv_puppi", &phi_sv_puppi, "phi_sv_puppi/F");

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

  PileUp * PUofficial = new PileUp();

  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Nov17.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring15_PU25_Startup.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  int nonOverlap = 0;

  vector<int> runList; runList.clear();
  vector<int> eventList; eventList.clear();

  if (checkOverlap) {
    std::ifstream fileEvents("overlap.txt");
    int Run, Event, Lumi;
    std::cout << "Non-overlapping events ->" << std::endl;
    while (fileEvents >> Run >> Event >> Lumi) {
      runList.push_back(Run);
      eventList.push_back(Event);
      std::cout << Run << ":" << Event << std::endl;
    }
    std::cout << std::endl;
  }
  std::ofstream fileOutput("overlap.out");

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

      //      cout << "Ok1" << endl;

      if (isDY) {
	for (unsigned int igentau=0; igentau<analysisTree.gentau_count; ++igentau) {
	  if (analysisTree.gentau_isPrompt[igentau]) isZTT = true;
	}
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
	  if (analysisTree.genparticles_pdgid[igen]==23) isZfound = true;
	  //	if (analysisTree.genparticles_pdgid[igen]==23 && analysisTree.genparticles_status[igen]==62) isZwithCode = true;
	  //	if (analysisTree.genparticles_pdgid[igen]==23 && analysisTree.genparticles_info[igen]==11) isZEE = true;
	  //	if (analysisTree.genparticles_pdgid[igen]==23 && analysisTree.genparticles_info[igen]==13) isZMM = true;
	  //	if (analysisTree.genparticles_pdgid[igen]==23 && analysisTree.genparticles_info[igen]==15) isZTT = true;
	  if (analysisTree.genparticles_isPrompt[igen]&&fabs(analysisTree.genparticles_pdgid[igen])==11) isZEE = true;
	  if (analysisTree.genparticles_isPrompt[igen]&&fabs(analysisTree.genparticles_pdgid[igen])==13) isZMM = true;
	}
	
	isZLL = isZEE || isZMM;
	
	ALL->Fill(0.);
	
	if (isZEE) ZEE->Fill(0.);
	if (isZMM) ZMM->Fill(0.);
	if (isZTT) ZTT->Fill(0.);
	if (!isZfound) GAM->Fill(0.);
      }

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
	if (discr=BTagDiscriminator)
	  nBTagDiscriminant = iBTag;
      }

      //      std::cout << "Isomu      : " << IsoMuonLeg << " : " << nIsoMuonLeg << std::endl;
      //      std::cout << "MuTauMuLeg   : " << MuonTauMuonLeg << " : " << nMuonTauMuonLeg << std::endl;
      //      std::cout << "MuTauTauLeg  : " << MuonTauTauLeg << " : " << nMuonTauTauLeg << std::endl;
      //      std::cout << "MuTauOverlap : " << MuonTauOverlap << " : " << nMuonTauOverlap << std::endl;
      //      std::cout << std::endl;
      //      continue;

      // vertex cuts no needed
      //      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      //      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      //      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
      //		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      //      if (dVertex>dVertexCut) continue;

      // tau selection
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	if (analysisTree.tau_decayModeFindingNewDMs[it]<=0.5) continue;
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
      float isoTauMin = 1e+10;
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

	  float isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex];

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
	    else if (isoTau<isoTauMin) {
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


      // applying inclusive selection
      if (applyInclusiveSelection) {

	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
	if (dilepton_veto)  continue;

	bool tauPass = analysisTree.tau_againstElectronVLooseMVA5[tauIndex]>0.5 &&
	  analysisTree.tau_againstMuonTight3[tauIndex]>0.5 &&
	  //	  analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex]<isoTauCut;
	  analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex] > 0.5;
	if (!tauPass) continue;
	
	bool muonPass = isoMuMin<isoMuonHighCut;
	if (!muonPass) continue;
	//	if (!os) continue;
	
      }
      
      //      cout << "Ok5" << endl;

      // met
      met = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      metphi = TMath::ATan2(analysisTree.pfmet_ey,analysisTree.pfmet_ex);

      // puppimet
      puppimet = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex + analysisTree.puppimet_ey*analysisTree.puppimet_ey);
      puppimetphi = TMath::ATan2(analysisTree.puppimet_ey,analysisTree.puppimet_ex);

      metcov00 = analysisTree.pfmet_sigxx;
      metcov01 = analysisTree.pfmet_sigxy;
      metcov10 = analysisTree.pfmet_sigyx;
      metcov11 = analysisTree.pfmet_sigyy;

      float mvamet_x = -999999;
      float mvamet_y = -999999;
      mvamet = -999999;
      mvametphi = -999999;
      mvacov00 = -999999;
      mvacov01 = -999999;
      mvacov10 = -999999;
      mvacov11 = -999999;

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
      
      float genmet_x = analysisTree.genmet_ex;
      float genmet_y = analysisTree.genmet_ey;
      genmet = TMath::Sqrt(genmet_x*genmet_x+genmet_y*genmet_y);
      genmetphi = TMath::ATan2(genmet_y,genmet_x);

      //      if (analysisTree.mvamet_count>0) {
      //	mvamet_x = analysisTree.mvamet_ex[metIndex];
      //	mvamet_y = analysisTree.mvamet_ey[metIndex];
      //	float mvamet_x2 = mvamet_x * mvamet_x;
      //	float mvamet_y2 = mvamet_y * mvamet_y;
      //
      //	mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
      //	mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
      //	mvacov00 = analysisTree.mvamet_sigxx[metIndex];
      //	mvacov01 = analysisTree.mvamet_sigxy[metIndex];
      //	mvacov10 = analysisTree.mvamet_sigyx[metIndex];
      //	mvacov11 = analysisTree.mvamet_sigyy[metIndex];
      //	std::cout << "MVA MET found" << std::endl;
      //      }
      //      else {
      //	std::cout << "No MVA MET found" << std::endl;
      //      }

      //      cout << "Ok6" << endl;

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
      float dPhiMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
				     analysisTree.pfmet_ex,analysisTree.pfmet_ey);
      mt_1 = TMath::Sqrt(2*met*pt_1*(1-TMath::Cos(dPhiMETMuon)));
      float dPhiPuppiMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
					  analysisTree.puppimet_ex,analysisTree.puppimet_ey);
      mt_puppi_1 = TMath::Sqrt(2*puppimet*pt_1*(1-TMath::Cos(dPhiPuppiMETMuon)));
      //      float dPhiMvaMETMuon = dPhiFrom2P(analysisTree.muon_px[muonIndex],analysisTree.muon_py[muonIndex],
      //					mvamet_x,mvamet_y);
      //      mt_mva_1 = TMath::Sqrt(2*mvamet*pt_1*(1-TMath::Cos(dPhiMvaMETMuon)));



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
      iso_2 = isoTauMin;
      m_2 = analysisTree.tau_mass[tauIndex];
      float dPhiMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.tau_py[tauIndex],
				    analysisTree.pfmet_ex,analysisTree.pfmet_ey);
      mt_2 = TMath::Sqrt(2*met*pt_2*(1-TMath::Cos(dPhiMETTau)));
      float dPhiPuppiMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.tau_py[tauIndex],
					 analysisTree.puppimet_ex,analysisTree.puppimet_ey);
      mt_puppi_2 = TMath::Sqrt(2*puppimet*pt_2*(1-TMath::Cos(dPhiPuppiMETTau)));
      //      float dPhiMvaMETTau = dPhiFrom2P(analysisTree.tau_px[tauIndex],analysisTree.muon_py[tauIndex],
      //				       mvamet_x,mvamet_y);
      //      mt_mva_2 = TMath::Sqrt(2*mvamet*pt_2*(1-TMath::Cos(dPhiMvaMETTau)));
      decayMode_2 = analysisTree.tau_decayMode[tauIndex];

      gen_match_1 = 6;
      gen_match_2 = 6;

      //      cout << "Ok7" << endl;

      if (isDY) {

	isZTT = false;
	isZL  = false;
	isZJ  = false;

	float minDR_1 = 0.2;
	float minDR_2 = 0.2;
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

      int tau_decayMode = analysisTree.tau_decayMode[tauIndex];

      m_sv = -9999;
      pt_sv = -9999;
      eta_sv = -9999;
      phi_sv = -9999;

      m_sv_puppi = -9999;
      pt_sv_puppi = -9999;
      eta_sv = -9999;
      phi_sv = -9999;


      bool validDecayMode = tau_decayMode<=4 || tau_decayMode>=10; // excluding 2prongs

      //      if (computeDitauMass && validDecayMode) {
      if (computeDitauMass) {
	// computing SV fit mass

	// define MET
	double measuredMETx = analysisTree.pfmet_ex;
	double measuredMETy = analysisTree.pfmet_ey;
	// define MET covariance
	TMatrixD covMET(2, 2);
	covMET[0][0] =  analysisTree.pfmet_sigxx;
	covMET[1][0] =  analysisTree.pfmet_sigxy;
	covMET[0][1] =  analysisTree.pfmet_sigyx;
	covMET[1][1] =  analysisTree.pfmet_sigyy;
	// define lepton four vectors
	std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, 
									pt_1, eta_1, phi_1, muonMass)); // tau -> muon decay (Pt, eta, phi, mass)
	//	if (tau_decayMode>11) tau_decayMode = 11;
	//	if (tau_decayMode>2&&tau_decayMode<10) tau_decayMode = 2;
	float tauMass = m_2;
	if (tau_decayMode==0) tauMass = pionMass;
	//	std::cout << "Tau decay mode = " << tau_decayMode << std::endl;
	measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay,  
									pt_2, eta_2, phi_2, tauMass, tau_decayMode)); // tau -> hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())

	// svfit mass pfmet --->
	SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
	algo.addLogM(false);  
	edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
	TH1::AddDirectory(false);  
	TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
	if (validDecayMode)
	  algo.shiftVisPt(true, inputFile_visPtResolution);
	else
	  algo.shiftVisPt(false,inputFile_visPtResolution);
	algo.integrateMarkovChain();
	
	m_sv = algo.getMass(); // return value is in units of GeV
	if ( algo.isValidSolution() ) {
	  //	  std::cout << "SVfit mass (pfmet)    = " << m_sv << std::endl;
	} else {
	  std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
	}
	
	pt_sv = algo.pt(); 
	eta_sv = algo.eta();
	phi_sv = algo.phi();

	delete inputFile_visPtResolution;

	// svfit mass puppimet --->
	double measuredPuppiMETx = analysisTree.puppimet_ex;
	double measuredPuppiMETy = analysisTree.puppimet_ey;
	SVfitStandaloneAlgorithm algoX(measuredTauLeptons, measuredPuppiMETx, measuredPuppiMETy, covMET, 0);
	//	SVfitStandaloneAlgorithm algoX(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
	algoX.addLogM(false);  
	//	edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
	TH1::AddDirectory(false);  
	inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
	if (validDecayMode)
	  algoX.shiftVisPt(true, inputFile_visPtResolution);
	else
	  algoX.shiftVisPt(false,inputFile_visPtResolution);
	algoX.integrateMarkovChain();
	
	m_sv_puppi = algoX.getMass(); // return value is in units of GeV
	if ( algo.isValidSolution() ) {
	  //	  std::cout << "SVfit mass (puppimet) = " << m_sv << std::endl;
	} else {
	  std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
	}
	
	pt_sv_puppi = algoX.pt(); 
	eta_sv_puppi = algoX.eta();
	phi_sv_puppi = algoX.phi();

	delete inputFile_visPtResolution;

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

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.tau_px[tauIndex],
					  analysisTree.tau_py[tauIndex],
					  analysisTree.tau_pz[tauIndex],
					  analysisTree.tau_e[tauIndex]);

      TLorentzVector dileptonLV = muonLV + tauLV;

      // visible mass
      m_vis = dileptonLV.M();
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

      // choosing mva met
      //      unsigned int metEMu = 0;
      //      for (unsigned int iMet=0; iMet<analysisTree.mvamet_count; ++iMet) {
      //	if (analysisTree.mvamet_channel[iMet]==1) metEMu = iMet;
      //      }


      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + tauLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + tauLV.Py();
      
      float vectorVisX = muonLV.Px() + tauLV.Px();
      float vectorVisY = muonLV.Py() + tauLV.Py();

      // computation of DZeta variable
      pzetamiss = analysisTree.pfmet_ex*zetaX + analysisTree.pfmet_ey*zetaY;
      pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;

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
	float energy = analysisTree.pfjet_e[jet];
        float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
        float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
        float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
        float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
        unsigned int chm = analysisTree.pfjet_chargedmulti[jet];
	unsigned int nem = analysisTree.pfjet_neutralmulti[jet];
        unsigned int npr = chm + nem;
	bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>2.4 || (elf<0.99 && chf>0 && chm>0 && absJetEta<=2.4)) && absJetEta<=3.0;
	if (absJetEta>3.0)
	  isPFJetId = phf<0.90 && nem>10 && absJetEta>3.0;
	if (!isPFJetId) continue;

	jetspt20.push_back(jet);

	if (absJetEta<bJetEtaCut && analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut) { // b-jet
	  bjets.push_back(jet);
	  if (jetPt>ptLeadingBJet) {
	    ptLeadingBJet = jetPt;
	    indexLeadingBJet = jet;
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



