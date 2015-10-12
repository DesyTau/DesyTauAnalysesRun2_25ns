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
#include "TH1D.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "TGraphAsymmErrors.h"

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

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

struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject, myobjectX;


int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");


  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonLowCut = cfg.get<float>("etaMuonLowCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");
  const bool  applyTauTauSelection = cfg.get<bool>("ApplyTauTauSelection");
  const bool  selectZToTauTauMuMu = cfg.get<bool>("SelectZToTauTauMuMu");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");

  // trigger
  const string muonTriggerName  = cfg.get<string>("MuonTriggerName");
  const string muonFilterName   = cfg.get<string>("MuonFilterName");
  const string muon17FilterName = cfg.get<string>("Muon17FilterName"); 
  const string muon8FilterName = cfg.get<string>("Muon8FilterName"); 

  TString MuonTriggerName(muonTriggerName);
  TString MuonFilterName(muonFilterName);

  TString Muon17FilterName(muon17FilterName);
  TString Muon8FilterName(muon8FilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // jet related cuts
  const float jetEtaCut      = cfg.get<float>("JetEtaCut");
  const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");
  const float jetPtHighCut   = cfg.get<float>("JetPtHighCut");
  const float jetPtLowCut    = cfg.get<float>("JetPtLowCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  // Run range
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  // lepton scale factors
  const string muonSfDataBarrel = cfg.get<string>("MuonSfDataBarrel");
  const string muonSfDataEndcap = cfg.get<string>("MuonSfDataEndcap");
  const string muonSfMcBarrel = cfg.get<string>("MuonSfMcBarrel");
  const string muonSfMcEndcap = cfg.get<string>("MuonSfMcEndcap");

  // **** end of configuration

  // Run-lumi selector
  std::vector<Period> periods;
    
  std::fstream inputFileStream("temp", std::ios::in);
  for(std::string s; std::getline(inputFileStream, s); )
    {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }


  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * histWeightsH = new TH1F("histWeightsH","",1,-0.5,0.5);

  // Histograms after selecting unique dimuon pair
  TH1F * ptLeadingMuSelH = new TH1F("ptLeadingMuSelH","",100,0,200);
  TH1F * ptTrailingMuSelH = new TH1F("ptTrailingMuSelH","",100,0,200);
  TH1F * etaLeadingMuSelH = new TH1F("etaLeadingMuSelH","",50,-2.5,2.5);
  TH1F * etaTrailingMuSelH = new TH1F("etaTrailingMuSelH","",50,-2.5,2.5);
  TH1F * massSelH = new TH1F("massSelH","",200,0,200);
  TH1F * metSelH  = new TH1F("metSelH","",200,0,400);

  TString scales[21] = {"M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","0",
			"P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};
  
  TH1F * massSelScaleH[21];
  TH1F * metSelScaleH[21];
  for (int iScale=0; iScale<21; ++iScale) {    
    massSelScaleH[iScale] = new TH1F("massSel"+scales[iScale]+"H","",200,0,200);
    metSelScaleH[iScale] = new TH1F("metSel"+scales[iScale]+"H","",200,0,400);
  }


  TH1F * nJets30SelH    = new TH1F("nJets30SelH","",11,-0.5,10.5);
  TH1F * nJets30etaCutSelH = new TH1F("nJets30etaCutSelH","",11,-0.5,10.5);
  TH1F * nJets20SelH    = new TH1F("nJets20SelH","",11,-0.5,10.5);
  TH1F * nJets20etaCutSelH = new TH1F("nJets20etaCutSelH","",11,-0.5,10.5);

  TH1F * HT30SelH       = new TH1F("HT30SelH","",50,0,500);
  TH1F * HT30etaCutSelH = new TH1F("HT30etaCutSelH","",50,0,500);
  TH1F * HT20SelH       = new TH1F("HT20SelH","",50,0,500);
  TH1F * HT20etaCutSelH = new TH1F("HT20etaCutSelH","",50,0,500);
  
  //Discriminiant histos
  TH1F * h_dimuonEta = new TH1F("dimuonEta","",50,-6,+6);//21Aug
  TH1F * h_dimuonEta_genMuMatch = new TH1F("dimuonEta_genMuMatch","",50,-6,+6);//21Aug
  TH1F * h_ptRatio = new TH1F ("ptRatio","",50,0,1);
  TH1F * h_dxy_muon1 =new TH1F ("dxy_muon1","",50,-0.02,0.02);
  TH1F * h_dxy_muon2 =new TH1F ("dxy_muon2","",50,-0.02,0.02);
  TH1F * h_dz_muon1 = new TH1F ("dz_muon1","",50,-0.1,0.1);
  TH1F * h_dz_muon2 = new TH1F ("dz_muon2","",50,-0.1,0.1);
  TH1F * h_dcaSigdxy_mu1 = new TH1F ("dcaSigdxy_mu1","",50,-4,4);
  TH1F * h_dcaSigdxy_mu2 = new TH1F ("dcaSigdxy_mu2","",50,-4,4);
  TH1F * h_dcaSigdz_mu1 = new TH1F ("dcaSigdz_mu1","",50,-4,4);
  TH1F * h_dcaSigdz_mu2 = new TH1F ("dcaSigdz_mu2","",50,-4,4);
  TH1F * h_dcaSigdxy_muon1 = new TH1F ("dcaSigdxy_muon1","",50,-4,4);
  TH1F * h_dcaSigdxy_muon2 = new TH1F ("dcaSigdxy_muon2","",50,-4,4);
  TH1F * h_dcaSigdz_muon1 = new TH1F ("dcaSigdz_muon1","",50,-4,4);
  TH1F * h_dcaSigdz_muon2 = new TH1F ("dcaSigdz_muon2","",50,-4,4);
  TH1F * h_dcaSigdxy_mu1_genMuMatch = new TH1F ("dcaSigdxy_mu1_genMuMatch","",50,-4,4);
  TH1F * h_dcaSigdxy_mu2_genMuMatch = new TH1F ("dcaSigdxy_mu2_genMuMatch","",50,-4,4);
  TH1F * h_dcaSigdz_mu1_genMuMatch = new TH1F ("dcaSigdz_mu1_genMuMatch","",50,-4,4);
  TH1F * h_dcaSigdz_mu2_genMuMatch = new TH1F ("dcaSigdz_mu2_genMuMatch","",50,-4,4);
  TH1F * h_phi_leadingMu_MET =new TH1F ("phi_leadingMu_MET","",50,0,3);
  TH1F * h_phi_trailingMu_MET =new TH1F ("phi_trailingMu_MET","",50,0,3);
  TH1F * h_dxy_muon1_dcaCut =new TH1F ("dxy_muon1_dcaCut","",50,-0.02,0.02);
  TH1F * h_dxy_muon2_dcaCut =new TH1F ("dxy_muon2_dcaCut","",50,-0.02,0.02);
  TH1F * h_dimuonEta_dcaCut = new TH1F("dimuonEta_dcaCut","",50,-6,+6);
  TH1F * h_ptRatio_dcaCut = new TH1F ("ptRatio_dcaCut","",50,0,1);

  TH1F * NumberOfVerticesH = new TH1F("NumberOfVerticesH","",51,-0.5,50.5);
  
  int nPtBins = 6;
  float ptBins[7] = {10,15,20,30,40,60,1000};

  int nPtBinsTrig = 7;
  float ptBinsTrig[8] = {10,15,19,24,30,40,60,1000};  
  
  int nEtaBins = 2;
  float etaBins[3] = {0,1.48,2.5}; 
  
  TString PtBins[6] = {"Pt10to15",
		       "Pt15to20",
		       "Pt20to30",
		       "Pt30to40",
		       "Pt40to60",
		       "PtGt60"};
  
  TString PtBinsTrig[7] = {"Pt10to15",
			   "Pt15to19",
			   "Pt19to24",
			   "Pt24to30",
			   "Pt30to40",
			   "Pt40to60",
			   "PtGt60"};

  TString EtaBins[2] = {"Barrel",
			"Endcap"};

  TH1F * ZMassEtaPtPass[2][6];
  TH1F * ZMassEtaPtFail[2][6];

  TH1F * ZMassMu23EtaPtPass[2][7];
  TH1F * ZMassMu23EtaPtFail[2][7];

  TH1F * ZMassMu8EtaPtPass[2][7];
  TH1F * ZMassMu8EtaPtFail[2][7];

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassEtaPtPass[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",60,60,120);
      ZMassEtaPtFail[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",60,60,120);
    }
    for (int iPt=0; iPt<nPtBinsTrig; ++iPt) {
      ZMassMu23EtaPtPass[iEta][iPt] = new TH1F("ZMassMu23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",60,60,120);
      ZMassMu23EtaPtFail[iEta][iPt] = new TH1F("ZMassMu23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",60,60,120);
      ZMassMu8EtaPtPass[iEta][iPt] = new TH1F("ZMassMu8"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",60,60,120);
      ZMassMu8EtaPtFail[iEta][iPt] = new TH1F("ZMassMu8"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",60,60,120);
    }
  }

  TH1F * ZMassMu23BarrelPass = new TH1F("ZMassMu23BarrelPass","",60,60,120);
  TH1F * ZMassMu23BarrelFail = new TH1F("ZMassMu23BarrelFail","",60,60,120);
  TH1F * ZMassMu23EndcapPass = new TH1F("ZMassMu23EndcapPass","",60,60,120);
  TH1F * ZMassMu23EndcapFail = new TH1F("ZMassMu23EndcapFail","",60,60,120);

  TH1F * ZMassMu8BarrelPass  = new TH1F("ZMassMu8BarrelPass", "",60,60,120);
  TH1F * ZMassMu8BarrelFail  = new TH1F("ZMassMu8BarrelFail", "",60,60,120);
  TH1F * ZMassMu8EndcapPass  = new TH1F("ZMassMu8EndcapPass", "",60,60,120);
  TH1F * ZMassMu8EndcapFail  = new TH1F("ZMassMu8EndcapFail", "",60,60,120);

  
  int nJetBins = 3;
  int nZPtBins = 5;
  float zPtBins[6] = {0,10,20,30,50,1000};

  TString NJetBins[3] = {"NJet0","NJet1","NJetGe2"};
  TString ZPtBins[5] = {"Pt0to10",
			"Pt10to20",
			"Pt20to30",
			"Pt30to50",
			"PtGt50"};

  TH1F * recoilZParalH[3];
  TH1F * recoilZPerpH[3];

  for (int iBin=0; iBin<nJetBins; ++iBin) {
    recoilZParalH[iBin] = new TH1F("recoilZParal"+NJetBins[iBin]+"H","",100,-200,200);
    recoilZPerpH[iBin] = new TH1F("recoilZPerp"+NJetBins[iBin]+"H","",100,-200,200);
  }
  

  int nFiles = 0;
  int nEvents = 0;
  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // reading vertex weights
  TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertDataFileName);
  TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertMcFileName);

  TH1F * vertexDataH = (TH1F*)fileDataNVert->Get(TString(vertHistName));
  TH1F * vertexMcH   = (TH1F*)fileMcNVert->Get(TString(vertHistName));

  float normVertexData = vertexDataH->GetSumOfWeights();
  float normVertexMc   = vertexMcH->GetSumOfWeights();

  vertexDataH->Scale(1/normVertexData);
  vertexMcH->Scale(1/normVertexMc);

  TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataBarrel);  // mu SF barrel data
  TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataEndcap); // mu SF endcap data
  TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcBarrel);  // mu SF barrel MC
  TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcEndcap); // mu SF endcap MC 
  
  TGraphAsymmErrors *hEffBarrelData = (TGraphAsymmErrors*)f10->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapData = (TGraphAsymmErrors*)f11->Get("ZMassEndcap");
  TGraphAsymmErrors *hEffBarrelMC = (TGraphAsymmErrors*)f12->Get("ZMassBarrel");
  TGraphAsymmErrors *hEffEndcapMC = (TGraphAsymmErrors*)f13->Get("ZMassEndcap");
  
  double * dataEffBarrel = new double[10];
  double * dataEffEndcap = new double[10];
  double * mcEffBarrel = new double[10];
  double * mcEffEndcap = new double[10];
  
  dataEffBarrel = hEffBarrelData->GetY();
  dataEffEndcap = hEffEndcapData->GetY();
  mcEffBarrel = hEffBarrelMC->GetY();
  mcEffEndcap = hEffEndcapMC->GetY();

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
 
  //----Attention----//
  //if(XSec!=1) nTotalFiles=20;
  //nTotalFiles=5;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
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
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      //------------------------------------------------

      if (!isData) 
	weight *=analysisTree.genweight;

      histWeightsH->Fill(float(0),weight);


      if (!isData) {
	if (applyPUreweighting) {
	  int binNvert = vertexDataH->FindBin(analysisTree.primvertex_count);
	  float_t dataNvert = vertexDataH->GetBinContent(binNvert);
	  float_t mcNvert = vertexMcH->GetBinContent(binNvert);
	  if (mcNvert < 1e-10){mcNvert=1e-10;}
	  float_t vertWeight = dataNvert/mcNvert;
	  weight *= vertWeight;
	  //	  cout << "NVert = " << analysisTree.primvertex_count << "   weight = " << vertWeight << endl;
	}
	if (applyTauTauSelection) {
	  unsigned int nTaus = 0;
	  if (analysisTree.gentau_count>0) {
	    //	  cout << "Generated taus present" << endl;
	    for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++itau) {
	      //	    cout << itau << "  : pt = " 
	      //		 << analysisTree.gentau_visible_pt[itau] 
	      //		 << "   eta = " <<  analysisTree.gentau_visible_eta[itau]
	      //		 << "   mother = " << int(analysisTree.gentau_mother[itau]) << endl;
	      if (int(analysisTree.gentau_mother[itau])==3) nTaus++;
	      
	    }
	  }
	  bool notTauTau = nTaus < 2;
	  //	  std::cout << "nTaus = " << nTaus << std::endl;
	  
	  if (selectZToTauTauMuMu&&notTauTau) { 
	    //	    std::cout << "Skipping event..." << std::endl;
	    //	    cout << endl;
	    continue;
	  }
	  if (!selectZToTauTauMuMu&&!notTauTau) { 
	    //	    std::cout << "Skipping event..." << std::endl;
	    //	    cout << endl;
	    continue;
	  }
	  //	  cout << endl;
	}
      }

      if (isData){
	
	
	bool lumi = false;
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;
	
	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {
	    
	    if ( num.c_str() ==  a.name ) {
	      //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      
	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		
		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	      //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	      
	      
	    }
	    
	  }
    
	if (!lumi) continue;
	//if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
	//std::remove("myinputfile");
      }     

      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;
      
      //std::cout << " Run : " << analysisTree.event_run << std::endl;
      
      bool isNewRun = true;
      if (allRuns.size()>0) {
	for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	  if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }
	
      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);
      
      bool isTriggerMuon = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(MuonTriggerName)) {
	  if (it->second==1)
	    isTriggerMuon = true;
	}
      }
      
      if (!isTriggerMuon) continue;
      
      unsigned int nMuonFilter = 0;
      bool isMuonFilter = false;

      unsigned int nMuon17Filter = 0;
      bool isMuon17Filter = false;

      unsigned int nMuon8Filter = 0;
      bool isMuon8Filter = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MuonFilterName) {
	  nMuonFilter = i;
	  isMuonFilter = true;
	}
	if (HLTFilter==Muon17FilterName) {
	  nMuon17Filter = i;
	  isMuon17Filter = true;
	}
	if (HLTFilter==Muon8FilterName) {
	  nMuon8Filter = i;
	  isMuon8Filter = true;
	}
      }
      if (!isMuonFilter) {
	cout << "Filter " << MuonFilterName << " not found " << endl;
	exit(-1);
      }
      
      // vertex cuts
      
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      
      // muon selection
      
      vector<unsigned int> allMuons; allMuons.clear();
      vector<unsigned int> idMuons; idMuons.clear();
      vector<unsigned int> isoMuons; isoMuons.clear();
      vector<float> isoMuonsValue; isoMuonsValue.clear();
      vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
      vector<bool> isMuonMatched23Filter; isMuonMatched23Filter.clear();
      vector<bool> isMuonMatched8Filter; isMuonMatched8Filter.clear();

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	allMuons.push_back(im);
	bool muPassed    = true;
	bool mu23Matched = false;
	bool mu8Matched  = false;
	if (analysisTree.muon_pt[im]<ptMuonLowCut) muPassed = false;
	if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) muPassed = false;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) muPassed = false;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) muPassed = false;
	if (!analysisTree.muon_isMedium[im]) muPassed = false;
	if (muPassed) idMuons.push_back(im);
	float absIso = analysisTree.muon_r03_sumChargedHadronPt[im];
	float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] + 
	  analysisTree.muon_r03_sumPhotonEt[im] - 
	  0.5*analysisTree.muon_r03_sumPUPt[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	absIso += neutralIso;
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut) muPassed = false;
	if (muPassed) { 
	  isoMuons.push_back(im);
	  isoMuonsValue.push_back(relIso);
	}
	//	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
	isMuonPassedIdIso.push_back(muPassed);
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nMuon17Filter] && analysisTree.trigobject_pt[iT]>23) 
	    mu23Matched = true;
	  if (analysisTree.trigobject_filters[iT][nMuon8Filter] && analysisTree.trigobject_pt[iT]>8) 
	    mu8Matched = true;
	  
	}
	isMuonMatched23Filter.push_back(mu23Matched);
	isMuonMatched8Filter.push_back(mu8Matched);
      }

      unsigned int indx1 = 0;
      unsigned int indx2 = 0;
      bool isIsoMuonsPair = false;
      float isoMin = 9999;
      if (isoMuons.size()>0) {
	for (unsigned int im1=0; im1<isoMuons.size(); ++im1) {
	  unsigned int index1 = isoMuons[im1];
	  bool isMu1matched = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMuonFilter] && 
		analysisTree.muon_pt[index1] > ptMuonHighCut &&
		fabs(analysisTree.muon_eta[index1]) < etaMuonHighCut) 
	      isMu1matched = true;
	  }
	  if (isMu1matched) {
	    for (unsigned int iMu=0; iMu<allMuons.size(); ++iMu) {
	      unsigned int indexProbe = allMuons[iMu];
	      if (index1==indexProbe) continue;
	      float q1 = analysisTree.muon_charge[index1];
	      float q2 = analysisTree.muon_charge[indexProbe];
	      if (q1*q2>0) continue;
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[indexProbe],analysisTree.muon_phi[indexProbe]);
	      if (dR<dRleptonsCut) continue; 
	      float ptProbe = TMath::Min(float(analysisTree.muon_pt[indexProbe]),float(ptBins[nPtBins]-0.1));
	      float absEtaProbe = fabs(analysisTree.muon_eta[indexProbe]);
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      int ptBinTrig = binNumber(ptProbe,nPtBinsTrig,ptBinsTrig);
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[indexProbe],
						  analysisTree.muon_py[indexProbe],
						  analysisTree.muon_pz[indexProbe],
						  muonMass);
	      float mass = (muon1+muon2).M();
	      if (isMuonPassedIdIso[iMu]) { 
		ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight); 
		// muon23 filter
		if (isMuonMatched23Filter[iMu]) 
		  ZMassMu23EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight); 
		else
		  ZMassMu23EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight); 
		if (ptProbe>24&&absEtaProbe<2.4) {
		  if (absEtaProbe<1.48) {
		    if (isMuonMatched23Filter[iMu]) 
		      ZMassMu23BarrelPass->Fill(mass,weight); 
		    else
		      ZMassMu23BarrelFail->Fill(mass,weight);
		  }
		  else {
		    if (isMuonMatched23Filter[iMu])
                      ZMassMu23EndcapPass->Fill(mass,weight);
                    else
                      ZMassMu23EndcapFail->Fill(mass,weight);
		  }
		}
		// muon8 filter
		if (isMuonMatched8Filter[iMu]) 
		  ZMassMu8EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight); 
		else
		  ZMassMu8EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight); 
		if (ptProbe>10&&absEtaProbe<2.4) {
		  if (absEtaProbe<1.48) {
		    if (isMuonMatched8Filter[iMu]) 
		      ZMassMu8BarrelPass->Fill(mass,weight); 
		    else
		      ZMassMu8BarrelFail->Fill(mass,weight);
		  }
		  else {
		    if (isMuonMatched8Filter[iMu])
                      ZMassMu8EndcapPass->Fill(mass,weight);
                    else
                      ZMassMu8EndcapFail->Fill(mass,weight);
		  }
		}
	      }
	      else
		ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
	    }
	  }
	  for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
	    unsigned int index2 = isoMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    bool isMu2matched = false;
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig>DRTrigMatch) continue;
	      if (analysisTree.trigobject_filters[iT][nMuonFilter] && 
		  analysisTree.muon_pt[index2] > ptMuonHighCut &&
		  fabs(analysisTree.muon_eta[index2]) < etaMuonHighCut) 
		isMu2matched = true;
	    }
	    bool isPairSelected = q1*q2 > 0;
	    if (oppositeSign) isPairSelected = q1*q2 < 0;
	    bool isTriggerMatch = (isMu1matched || isMu2matched);
	    float dRmumu = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				  analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	    if (isTriggerMatch && isPairSelected && dRmumu>dRleptonsCut) {
	      bool sumIso = isoMuonsValue[im1]+isoMuonsValue[im2];
	      if (sumIso<isoMin) {
		isIsoMuonsPair = true;
		isoMin = sumIso;
		if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]) {
		  indx1 = index1;
		  indx2 = index2;
		}
		else {
		  indx2 = index1;
		  indx1 = index2;
		}
	      }
	    }
	  }
	}
      }

      if (isIsoMuonsPair) {      
	//match to genparticles
	bool genmatch_m1 = false, genmatch_m2 = false;
	if (!isData) {
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	    if (fabs(analysisTree.genparticles_pdgid[igen])==13 &&
		analysisTree.genparticles_status[igen]==1) {

	       TLorentzVector gen_mu1; gen_mu1.SetPxPyPzE(analysisTree.genparticles_px[igen],
							 analysisTree.genparticles_py[igen],
							 analysisTree.genparticles_pz[igen],
							  analysisTree.genparticles_e[igen]);
	       TLorentzVector gen_mu2; gen_mu2.SetPxPyPzE(analysisTree.genparticles_px[igen],
							 analysisTree.genparticles_py[igen],
							 analysisTree.genparticles_pz[igen],
							  analysisTree.genparticles_e[igen]);
	      double deltaR_mu1 = deltaR(gen_mu1.Eta(), gen_mu1.Phi(),
					 analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1]);
	      if(deltaR_mu1 < 0.3)genmatch_m1 = true;

	      double deltaR_mu2 = deltaR(gen_mu2.Eta(), gen_mu2.Phi(),
					 analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	      if(deltaR_mu2 < 0.3)genmatch_m2 = true;
   
	    }
	  }
	}

	// selecting good jets --->

	// HT variables
	float HT30 = 0;
	float HT20 = 0;
	float HT30etaCut = 0;
	float HT20etaCut = 0;

	// number of jets
	int nJets30 = 0;
	int nJets20 = 0;
	int nJets30etaCut = 0;
	int nJets20etaCut = 0;
	
	for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	  float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	  if (absJetEta>jetEtaCut) continue;
	  
	  float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			     analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1]);
	  if (dR1<dRJetLeptonCut) continue;
	  
	  float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			     analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	  
	  if (dR2<dRJetLeptonCut) continue;
	  
	  // pfJetId

	  float energy = analysisTree.pfjet_e[jet];
	  float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
	  float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
	  float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
	  float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
	  float chm = analysisTree.pfjet_chargedmulti[jet];
	  float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
	  bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>2.4 || (elf<0.99 && chf>0 && chm>0));
	  if (!isPFJetId) continue;
	  
	  if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
	    nJets30++;
	    HT30 += analysisTree.pfjet_pt[jet];
	    if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
	      HT30etaCut += analysisTree.pfjet_pt[jet];
	      nJets30etaCut++;
	    }
	  }

	  if (analysisTree.pfjet_pt[jet]>jetPtLowCut) {
	    nJets20++;
	    HT20 += analysisTree.pfjet_pt[jet]; 
	    if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
              HT20etaCut += analysisTree.pfjet_pt[jet];
              nJets20etaCut++;
            }
	  }	  


	}

	if (!isData && applyLeptonSF) {
	  float ptMu1 = TMath::Min(float(59),analysisTree.muon_pt[indx1]);
	  float ptMu2 = TMath::Min(float(59),analysisTree.muon_pt[indx2]);
	  int ptBinMu1 = binNumber(ptMu1,nPtBins,ptBins);
	  int ptBinMu2 = binNumber(ptMu2,nPtBins,ptBins);
	  float absEtaMu1 = fabs(analysisTree.muon_eta[indx1]);
	  float absEtaMu2 = fabs(analysisTree.muon_eta[indx2]);
	  float dataMu1 = 1;
	  float mcMu1 = 1;
	  float dataMu2 = 1;
	  float mcMu2 = 1;
	  if(absEtaMu1 < 1.48){
	    dataMu1 = dataEffBarrel[ptBinMu1];
	    mcMu1 = mcEffBarrel[ptBinMu1];
	  }
	  else {
	    dataMu1 = dataEffEndcap[ptBinMu1];
	    mcMu1 = mcEffEndcap[ptBinMu1];
	  }
	  if(absEtaMu2 < 1.48){
	    dataMu2 = dataEffBarrel[ptBinMu2];
	    mcMu2 = mcEffBarrel[ptBinMu2];
	  }
	  else {
	    dataMu2 = dataEffEndcap[ptBinMu2];
	    mcMu2 = mcEffEndcap[ptBinMu2];
	  }
	  float wMu1 = dataMu1/mcMu1;
	  float wMu2 = dataMu2/mcMu2;
	  //	  cout << "Muons SF : Mu1 = " << wMu1 << "   Mu2 = " << wMu2 << endl;
	  weight = weight*wMu1*wMu2;
	}

	TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[indx1],
					analysisTree.muon_py[indx1],
					analysisTree.muon_pz[indx1],
					muonMass);

	TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[indx2],
					analysisTree.muon_py[indx2],
					analysisTree.muon_pz[indx2],
					muonMass);

	TLorentzVector dimuon = mu1 + mu2;
	float massSel = dimuon.M();

	massSelH->Fill(massSel,weight);
	for (int iScale=0; iScale<21; ++ iScale) {
	  float scaleFactor = 0.98 + 0.002*float(iScale);
	  massSelScaleH[iScale]->Fill(massSel*scaleFactor,weight);
	}


	if (massSel>20) {
	  ptLeadingMuSelH->Fill(analysisTree.muon_pt[indx1],weight);
	  ptTrailingMuSelH->Fill(analysisTree.muon_pt[indx2],weight);
	  etaLeadingMuSelH->Fill(analysisTree.muon_eta[indx1],weight);
	  etaTrailingMuSelH->Fill(analysisTree.muon_eta[indx2],weight);
	  
	  nJets30SelH->Fill(double(nJets30),weight);
	  nJets20SelH->Fill(double(nJets20),weight);
	  nJets30etaCutSelH->Fill(double(nJets30etaCut),weight);
	  nJets20etaCutSelH->Fill(double(nJets20etaCut),weight);
	  
	  HT30SelH->Fill(double(HT30),weight);
	  HT20SelH->Fill(double(HT20),weight);
	  HT30etaCutSelH->Fill(double(HT30etaCut),weight);
	  HT20etaCutSelH->Fill(double(HT20etaCut),weight);

	  // cout << "dxy (mu1) = " << analysisTree.muon_dxy[indx1] << "   error = " << analysisTree.muon_dxyerr[indx1] << std::endl;
	  //cout << "dxy (mu2) = " << analysisTree.muon_dxy[indx2] << "   error = " << analysisTree.muon_dxyerr[indx2] << std::endl;
	  //	 vector<int>indexMuTau; indexMuTau.clear();  
	  //	 if (selectZToTauTauMuMu && indexMuTau.size()!=2) continue;
	  float dimuonEta = dimuon.Eta();
	  float dimuonPt = dimuon.Pt();
	  float sumMuonPt = (mu1.Pt()+mu2.Pt());
	  float ptRatio = 0.0;
	  if (sumMuonPt != 0)
	    ptRatio = (dimuonPt/sumMuonPt);
	  float dcaSigdxy_muon1 = 0.0;
	  float dcaSigdxy_muon2 = 0.0;
	  float dcaSigdz_muon1  = 0.0;
	  float dcaSigdz_muon2  =0.0;
	  float dcaSigdxy_mu1 = 0.0;
	  float dcaSigdxy_mu2 = 0.0;
          float dcaSigdz_mu1  = 0.0;
          float dcaSigdz_mu2  =0.0;
	  float dcaSigdxy_mu1_genMuMatch = 0.0;
          float dcaSigdxy_mu2_genMuMatch = 0.0;
          float dcaSigdz_mu1_genMuMatch  = 0.0;
          float dcaSigdz_mu2_genMuMatch  =0.0;
	  if (analysisTree.muon_dxyerr[indx1] != 0){
	    dcaSigdxy_mu1=(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]);
	    dcaSigdxy_muon1= log10(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]);
	  }
		      
	  if (analysisTree.muon_dxyerr[indx2] != 0){
	    dcaSigdxy_mu2=(analysisTree.muon_dxy[indx2]/analysisTree.muon_dxyerr[indx2]);
	    dcaSigdxy_muon2= log10(analysisTree.muon_dxy[indx2]/analysisTree.muon_dxyerr[indx2]);
	  }

	  if (analysisTree.muon_dzerr[indx1] != 0){
	    dcaSigdz_mu1 =(analysisTree.muon_dz[indx1]/analysisTree.muon_dzerr[indx1]);
	    dcaSigdz_muon1 = log10(analysisTree.muon_dz[indx1]/analysisTree.muon_dzerr[indx1]);
	  }
	  
	  if (analysisTree.muon_dzerr[indx2] != 0){
	    dcaSigdz_mu2 = (analysisTree.muon_dz[indx2]/analysisTree.muon_dzerr[indx2]);
	    dcaSigdz_muon2 = log10(analysisTree.muon_dz[indx2]/analysisTree.muon_dzerr[indx2]);
	  }
	  	
	  //filling the histograms for discriminators
	  h_dimuonEta->Fill(dimuonEta,weight);
          if (genmatch_m1 && genmatch_m2) h_dimuonEta_genMuMatch->Fill(dimuonEta,weight);
	  h_ptRatio->Fill(ptRatio,weight);
	  h_dxy_muon1->Fill(analysisTree.muon_dxy[indx1],weight);
	  h_dxy_muon2->Fill(analysisTree.muon_dxy[indx2],weight);
	  h_dz_muon1->Fill(analysisTree.muon_dz[indx1],weight);
	  h_dz_muon2->Fill(analysisTree.muon_dz[indx2],weight);
	  h_dcaSigdxy_muon1->Fill(dcaSigdxy_muon1,weight);
	  h_dcaSigdxy_muon2->Fill(dcaSigdxy_muon2,weight);
	  h_dcaSigdz_muon1->Fill(dcaSigdz_muon1,weight);
	  h_dcaSigdz_muon2->Fill(dcaSigdz_muon2,weight);
	  h_dcaSigdxy_mu1->Fill(dcaSigdxy_mu1,weight);
          h_dcaSigdxy_mu2->Fill(dcaSigdxy_mu2,weight);
          h_dcaSigdz_mu1->Fill(dcaSigdz_mu1,weight);
          h_dcaSigdz_mu2->Fill(dcaSigdz_mu2,weight);
	  // if(genmatch_m1 && genmatch_mu2)
	  //h_dcaSigdxy_mu1_genMuMatch->Fill(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]);
	  if (genmatch_m1 && genmatch_m2) h_dcaSigdxy_mu1_genMuMatch->Fill(dcaSigdxy_mu1,weight);
          if (genmatch_m1 && genmatch_m2) h_dcaSigdxy_mu2_genMuMatch->Fill(dcaSigdxy_mu2,weight);
          if (genmatch_m1 && genmatch_m2) h_dcaSigdz_mu1_genMuMatch->Fill(dcaSigdz_mu1,weight);
          if (genmatch_m1 && genmatch_m2) h_dcaSigdz_mu2_genMuMatch->Fill(dcaSigdz_mu2,weight);
	  
	  h_phi_leadingMu_MET->Fill((analysisTree.muon_phi[indx1]-analysisTree.pfmet_phi),weight);
	  h_phi_trailingMu_MET->Fill((analysisTree.muon_phi[indx2]-analysisTree.pfmet_phi),weight);
	  
	  if(!isData){
	    if (dcaSigdxy_muon1<1.8 && dcaSigdxy_muon1>0.2) h_dxy_muon1_dcaCut->Fill(analysisTree.muon_dxy[indx1],weight);
	    if (dcaSigdxy_muon1<1.8 && dcaSigdxy_muon1>0.2) h_dxy_muon2_dcaCut->Fill(analysisTree.muon_dxy[indx2],weight);
	    if (dcaSigdxy_muon1<1.8 && dcaSigdxy_muon1>0.2) h_dimuonEta_dcaCut->Fill(dimuonEta,weight);
	    if (dcaSigdxy_muon1<1.8 && dcaSigdxy_muon1>0.2) h_ptRatio_dcaCut->Fill(ptRatio,weight);
	  }
	  
	  NumberOfVerticesH->Fill(float(analysisTree.primvertex_count),weight);
	  float metSel = sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex+analysisTree.pfmet_ey*analysisTree.pfmet_ey);
	  metSelH->Fill(metSel,weight);
	  for (int iScale=0; iScale<21; ++ iScale) {
	    float scaleFactor = 0.9 + 0.01*float(iScale);
	    metSelScaleH[iScale]->Fill(metSel*scaleFactor,weight);
	  }


	  if (massSel>70&&massSel<110) {
	    float unitX = dimuon.Px()/dimuon.Pt();
	    float unitY = dimuon.Py()/dimuon.Pt();
	    float phiUnit = TMath::ATan2(unitY,unitX);
	    float recoilParal = analysisTree.pfmet_ex*unitX + analysisTree.pfmet_ey*unitY;
	    float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
	    float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
	    float recoilPerp = analysisTree.pfmet_ex*perpUnitX + analysisTree.pfmet_ey*perpUnitY;
	    int jetBin = 0;
	    if (nJets30==1)
	      jetBin = 1;
	    else if (nJets30>1)
	      jetBin = 2;
	    recoilZParalH[jetBin]->Fill(recoilParal,weight);
	    recoilZPerpH[jetBin]->Fill(recoilPerp,weight);

	  }

	}

      }

      if (isIsoMuonsPair) selEventsIsoMuons++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                     = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  //cout << "weight used:" << weight << std::endl;

  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  
  file->Write();
  file->Close();
  delete file;
  
  }



