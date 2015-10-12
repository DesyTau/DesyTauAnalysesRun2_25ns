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


const float MuMass = 0.105658367;

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
} myobject;

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;
   float ht=0;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

  // kinematic cuts on electron
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float etaElectronHighCut = cfg.get<float>("etaElectronHighCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronCut     = cfg.get<float>("isoElectronCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool  oppositeSign   = cfg.get<bool>("OppositeSign");

  //trigger
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
  const string electronHLTName = cfg.get<string>("ElectronHLTName");
  const string electronHLTFilterName = cfg.get<string>("ElectronHLTFilterName");
  const string electronEle23FilterName = cfg.get<string>("ElectronEle23FilterName");
  const string electronEle12FilterName = cfg.get<string>("ElectronEle12FilterName");

  TString ElectronHLTName(electronHLTName);
  TString ElectronHLTFilterName(electronHLTFilterName);
  TString ElectronEle23FilterName(electronEle23FilterName);
  TString ElectronEle12FilterName(electronEle12FilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // Run range
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  // lepton scale factors
  const string eleSfDataBarrel = cfg.get<string>("EleSfDataBarrel");
  const string eleSfDataEndcap = cfg.get<string>("EleSfDataEndcap");
  const string eleSfMcBarrel = cfg.get<string>("EleSfMcBarrel");
  const string eleSfMcEndcap = cfg.get<string>("EleSfMcEndcap");

  std::vector<Period> periods;
    
  std::fstream inputFileStream("temp", std::ios::in);
  for(std::string s; std::getline(inputFileStream, s); )
    {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  
  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1F * histWeightsH = new TH1F("histWeightsH","",1,-0.5,0.5);
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

   TString scales[21] = {"M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","0",
			 "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};

  TH1F * MetH = new TH1F("MetH","",200,0.,400.);
  TH1F * MetScaleH[21];
  for (int iScale=0; iScale<21; ++iScale) {    
    MetScaleH[iScale] = new TH1F("Met"+scales[iScale]+"H","",200,0,400);
  }


  TH1F * JPsiMassAllElectronsH = new TH1F("JPsiMassAllElectronsH","",200,2,4);
  TH1F * JPsiMassAllElectronsDRCutH = new TH1F("JPsiMassAllElectronsDRCutH","",200,2,4);
  TH1F * JPsiMassIdLooseElectronsH = new TH1F("JPsiMassIdLooseElectronsH","",200,2,4);
  TH1F * JPsiMassIdLooseElectronsDRCutH = new TH1F("JPsiMassIdLooseElectronsDRCutH","",200,2,4);
  TH1F * JPsiMassIdTightElectronsH = new TH1F("JPsiMassIdTightElectronsH","",200,2,4);
  TH1F * JPsiMassIdTightElectronsDRCutH = new TH1F("JPsiMassIdTightElectronsDRCutH","",200,2,4);
  
  TH1F * YpsilonMassAllElectronsH = new TH1F("YpsilonMassAllElectronsH","",400,8,12);
  TH1F * YpsilonMassAllElectronsDRCutH = new TH1F("YpsilonMassAllElectronsDRCutH","",400,8,12);
  TH1F * YpsilonMassIdLooseElectronsH = new TH1F("YpsilonMassIdLooseElectronsH","",400,8,12);
  TH1F * YpsilonMassIdLooseElectronsDRCutH = new TH1F("YpsilonMassIdLooseElectronsDRCutH","",400,8,12);
  TH1F * YpsilonMassIdTightElectronsH = new TH1F("YpsilonMassIdTightElectronsH","",400,8,12);
  TH1F * YpsilonMassIdTightElectronsDRCutH = new TH1F("YpsilonMassIdTightElectronsDRCutH","",400,8,12);
  
  TH1F * ZMassAllElectronsH = new TH1F("ZMassAllElectronsH","",200,0,200);
  TH1F * ZMassAllElectronsDRCutH = new TH1F("ZMassAllElectronsDRCutH","",400,60,120);
  TH1F * ZMassIdLooseElectronsH = new TH1F("ZMassIdLooseElectronsH","",200,0,200);
  TH1F * ZMassIdLooseElectronsDRCutH = new TH1F("ZMassIdLooseElectronsDRCutH","",200,0,200);
  TH1F * ZMassIdTightElectronsH = new TH1F("ZMassIdTightElectronsH","",200,0,200);
  TH1F * ZMassIdTightElectronsDRCutH = new TH1F("ZMassIdTightElectronsDRCutH","",200,0,200);
  TH1F * ZMassIsoElectronsH = new TH1F("ZMassIsoElectronsH","",200,0,200);
  TH1F * ZMassIsoElectronsDRCutH = new TH1F("ZMassIsoElectronsDRCutH","",200,0,200);
  TH1F * ZMassIsoTightElectronsH = new TH1F("ZMassIsoTightElectronsH","",200,0,200);
  TH1F * ZMassIsoTightElectronsDRCutH = new TH1F("ZMassIsoTightElectronsDRCutH","",200,0,200);
  TH1F * ZMassIsoTightElectronsUniqueH = new TH1F("ZMassIsoTightElectronsUniqueH","",200,0,200);

  TH1F * ZMassIsoTightElectronsUniqueScaleH[21];
  for (int iScale=0; iScale<21; ++iScale) {    
    ZMassIsoTightElectronsUniqueScaleH[iScale]  = new TH1F("ZMassIsoTightElectronsUnique"+scales[iScale]+"H","",200,0,200);
  }

  //new histograms
  TH1F * ptLeadingAllElectronsH = new TH1F("ptLeadingAllElectronsH","",200,0,200);
  TH1F * etaLeadingAllElectronsH = new TH1F("etaLeadigAllElectrons","",60,-3,3);
  TH1F * ptLeadingIdLooseElectronH = new TH1F("ptLeadingIdLooseElectronH","",200,0,200);
  TH1F * etaLeadingIdLooseElectronH = new TH1F("etaLeadingIdLooseElectronH","",60,-3,3);
  TH1F * ptLeadingIdLooseElectronDRCutH = new TH1F("ptLeadingIdLooseElectronDRCutH","",200,0,200);
  TH1F * etaLeadingIdLooseElectronDRCutH = new TH1F("etaLeadingIdLooseElectronDRCutH","",60,-3,3);
  TH1F * ptLeadingIdTightElectrons = new TH1F("ptLeadingIdTightElectrons","",200,0,200);
  TH1F * etaLeadingIdTightElectrons = new TH1F("etaLeadingIdTightElectrons","",60,-3,3);
  TH1F * ptLeadingIdTightElectronsDRCutH = new TH1F("ptLeadingIdTightElectronsDRCutH","",200,0,200);
  TH1F * etaLeadingIdTightElectronsDRCutH = new TH1F("etaLeadingIdTightElectronsDRCutH","",60,-3,3);
  TH1F * ptLeadingIsoTightElectronsH = new TH1F("ptLeadingIsoTightElectronsH","",200,0,200);
  TH1F * etaLeadingIsoTightElectronsH = new TH1F("etaLeadingIsoTightElectronsH","",60,-3,3);
  TH1F * ptLeadingIsoTightElectronsDRCutH = new TH1F("ptLeadingIsoTightElectronsDRCutH","",200,0,200);
  TH1F * etaLeadingIsoTightElectronsDRCutH = new TH1F("etaLeadingIsoTightElectronsDRCutH","",60,-3,3);
  TH1F * ptLeadingIsoTightElectronsUniqueH = new TH1F("ptLeadingIsoTightElectronsUniqueH","",200,0,200);
  TH1F * etaLeadingIsoTightElectronsUniqueH = new TH1F("etaLeadingIsoTightElectronsUniqueH","",60,-3,3);
  
  TH1F * ptTrailingAllElectronsH = new TH1F("ptAllElectronsH","",200,0,200);
  TH1F * etaTrailingAllElectronsH = new TH1F("etaAllElectrons","",60,-3,3);
  TH1F * ptTrailingIdLooseElectronH = new TH1F("ptTrailingIdLooseElectronH","",200,0,200);
  TH1F * etaTrailingIdLooseElectronH = new TH1F("etaTrailingIdLooseElectronH","",60,-3,3);
  TH1F * ptTrailingIdLooseElectronDRCutH = new TH1F("ptTrailingIdLooseElectronDRCutH","",200,0,200);
  TH1F * etaTrailingIdLooseElectronDRCutH = new TH1F("etaTrailingIdLooseElectronDRCutH","",60,-3,3);
  TH1F * ptTrailingIdTightElectrons = new TH1F("ptTrailingIdTightElectrons","",200,0,200);
  TH1F * etaTrailingIdTightElectrons = new TH1F("etaTrailingIdTightElectrons","",60,-3,3);
  TH1F * ptTrailingIdTightElectronsDRCutH = new TH1F("ptTrailingIdTightElectronsDRCutH","",200,0,200);
  TH1F * etaTrailingIdTightElectronsDRCutH = new TH1F("etaTrailingIdTightElectronsDRCutH","",60,-3,3);
  TH1F * ptTrailingIsoTightElectronsH = new TH1F("ptTrailingIsoTightElectronsH","",200,0,200);
  TH1F * etaTrailingIsoTightElectronsH = new TH1F("etaTrailingIsoTightElectronsH","",60,-3,3);
  TH1F * ptTrailingIsoTightElectronsDRCutH = new TH1F("ptTrailingIsoTightElectronsDRCutH","",200,0,200);
  TH1F * etaTrailingIsoTightElectronsDRCutH = new TH1F("etaTrailingIsoTightElectronsDRCutH","",60,-3,3);
  TH1F * ptTrailingIsoTightElectronsUniqueH = new TH1F("ptTrailingIsoTightElectronsUniqueH","",200,0,200);
  TH1F * etaTrailingIsoTightElectronsUniqueH = new TH1F("etaTrailingIsoTightElectronsUniqueH","",60,-3,3);

  TH1F * NumberOfVerticesH = new TH1F("NumberOfVerticesH","",51,-0.5,50.5);
  
  int nPtBins = 7;
  float ptBins[8] = {10,15,20,25,30,40,60,1000};

  int nPtBinsTrig = 7;
  float ptBinsTrig[8] = {10,13,18,24,30,40,60,1000};

  int nEtaBins = 2;
  float etaBins[3] = {0,1.48,2.5}; 

  TString PtBins[7] = {"Pt10to15",
		       "Pt15to20",
		       "Pt20to25",
		       "Pt25to30",
		       "Pt30to40",
		       "Pt40to60",
		       "PtGt60"};

  TString PtBinsTrig[7] = {"Pt10to13",
			   "Pt13to18",
			   "Pt18to24",
			   "Pt24to30",
			   "Pt30to40",
			   "Pt40to60",
			   "PtGt60"};

  TString EtaBins[2] = {"Barrel",
			"Endcap"};

  TH1F * ZMassEtaPtPass[2][7];
  TH1F * ZMassEtaPtFail[2][7];

  TH1F * ZMassEle23EtaPtPass[2][7];
  TH1F * ZMassEle23EtaPtFail[2][7];

  TH1F * ZMassEle12EtaPtPass[2][7];
  TH1F * ZMassEle12EtaPtFail[2][7];


  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassEtaPtPass[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",60,60,120);
      ZMassEtaPtFail[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",60,60,120);
    }
    for (int iPt=0; iPt<nPtBinsTrig; ++iPt) {
      ZMassEle23EtaPtPass[iEta][iPt] = new TH1F("ZMassEle23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",60,60,120);
      ZMassEle23EtaPtFail[iEta][iPt] = new TH1F("ZMassEle23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",60,60,120);
      ZMassEle12EtaPtPass[iEta][iPt] = new TH1F("ZMassEle12"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",60,60,120);
      ZMassEle12EtaPtFail[iEta][iPt] = new TH1F("ZMassEle12"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",60,60,120);
    }
  }

  TH1F * ZMassEle23BarrelPass = new TH1F("ZMassEle23BarrelPass","",60,60,120);
  TH1F * ZMassEle23BarrelFail = new TH1F("ZMassEle23BarrelFail","",60,60,120);

  TH1F * ZMassEle23EndcapPass = new TH1F("ZMassEle23EndcapPass","",60,60,120);
  TH1F * ZMassEle23EndcapFail = new TH1F("ZMassEle23EndcapFail","",60,60,120);

  TH1F * ZMassEle12BarrelPass = new TH1F("ZMassEle12BarrelPass","",60,60,120);
  TH1F * ZMassEle12BarrelFail = new TH1F("ZMassEle12BarrelFail","",60,60,120);

  TH1F * ZMassEle12EndcapPass = new TH1F("ZMassEle12EndcapPass","",60,60,120);
  TH1F * ZMassEle12EndcapFail = new TH1F("ZMassEle12EndcapFail","",60,60,120);


  // reweighting for vertices
  string cmsswBase = (getenv ("CMSSW_BASE"));

  // reading vertex weights
  TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertDataFileName);
  TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+vertMcFileName);

  TH1F * hvertWeight = new TH1F("hvertWeight","",40,0,1);
  TH1F * hNvert = new TH1F("hNvert","",51,-0.5,50.5);
  
  TH1D *hNvertData = (TH1D*)fileDataNVert->Get(TString(vertHistName));
  TH1D *hNvertMC =(TH1D*)fileMcNVert->Get(TString(vertHistName));
  Float_t normData = hNvertData->GetSumOfWeights();   
  Float_t normMC =  hNvertMC->GetSumOfWeights();
  hNvertData->Scale(1/normData);
  hNvertMC->Scale(1/normMC);
  
  TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+eleSfDataBarrel);  // ele SF barrel data
  TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+eleSfDataEndcap); // ele SF endcap data
  TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+eleSfMcBarrel);  // ele SF barrel MC
  TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+eleSfMcEndcap); // ele SF endcap MC 
  
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
  
  
  // Float_t normDataBarrel =hEffBarrelData->GetSumOfWeights(); 
  // Float_t normDataEndcap =hEffEndcapData->GetSumOfWeights(); 
  // Float_t normMCBarrel =hEffBarrelMC->GetSumOfWeights(); 
  // Float_t normMCEndcap =hEffEndcapMC->GetSumOfWeights();
  
  // hEffBarrelData->Scale(1/normDataBarrel);
  // hEffEndcapData->Scale(1/normDataEndcap);
  // hEffBarrelMC->Scale(1/normMCBarrel);
  // hEffEndcapMC->Scale(1/normMCEndcap);
  
  
  int nFiles = 0;
  int nEvents = 0;
  int selEventsAllElectrons = 0;
  int selEventsIdLooseElectrons = 0;
  int selEventsIdTightElectrons = 0;
  int selEventsIsoElectrons = 0;
  int selEventsIsoTightElectrons = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;
		
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

      if (!isData) {
	weight *= analysisTree.genweight;
	hNvert->Fill(analysisTree.primvertex_count,weight);
      }
      histWeightsH->Fill(float(0),weight);

      if (!isData && applyPUreweighting) {
	//reweighting
	int binNvert = hNvert->FindBin(analysisTree.primvertex_count);
	float_t dataNvert = hNvertData->GetBinContent(binNvert);
	float_t mcNvert = hNvertMC->GetBinContent(binNvert);
	if (mcNvert < 1e-10){mcNvert=1e-10;}
	float_t vertWeight = dataNvert/mcNvert;
	hvertWeight->Fill(vertWeight);
	weight *= vertWeight;
	//	cout << "NVert = " << analysisTree.primvertex_count << "  : " << vertWeight << endl;
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
	      //	      std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      
	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		
		//   cout<<b->lower<<"  "<<b->bigger<<endl;
                if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	      // std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	      

	    }
	    
	  }
	if (!lumi) continue;
	

	if (analysisTree.event_run<RunMin)
	  RunMin = analysisTree.event_run;
	
	if (analysisTree.event_run>RunMax)
	  RunMax = analysisTree.event_run;
	
      }

      bool isTriggerElectron = false;
      
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(ElectronHLTName)) {
	  //	  std::cout << ElectronHLTName << " : " << it->second << std::endl;
	  if (it->second==1)
            isTriggerElectron = true;
	}
      }
      
      unsigned int nElectronFilter = 0;
      bool isElectronFilter = false;
      
      unsigned int nEle23Filter = 0;
      bool isEle23Filter = false;
      
      unsigned int nEle12Filter = 0;
      bool isEle12Filter = false;
      
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==ElectronHLTFilterName) {
	  nElectronFilter = i;
	  isElectronFilter = true;
	}
	if (HLTFilter==ElectronEle23FilterName) {
	  nEle23Filter = i;
	  isEle23Filter = true;
	}
	if (HLTFilter==ElectronEle12FilterName) {
	  nEle12Filter = i;
	  isEle12Filter = true;
	}
      }

      if (!isElectronFilter) {
	cout << "Filter " << ElectronHLTFilterName << " not found" << endl;
	exit(-1);
      }

      if (!isEle23Filter) {
	//	cout << "Filter " << ElectronEle23FilterName << " not found" << endl;
      }

      if (!isEle12Filter) {
	//	cout << "Filter " << ElectronEle12FilterName << " not found" << endl;
      }

      if (applyTrigger && !isTriggerElectron) continue;

      // vertex cuts
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
     
      // electron selection
      vector<unsigned int> allElectrons; allElectrons.clear();
      vector<unsigned int> idTightElectrons; idTightElectrons.clear();
      vector<unsigned int> idLooseElectrons; idLooseElectrons.clear();
      vector<unsigned int> isoElectrons; isoElectrons.clear();
      vector<unsigned int> isoTightElectrons; isoTightElectrons.clear();
      vector<bool> allElectronsIsLeg23Matched; allElectronsIsLeg23Matched.clear();
      vector<bool> allElectronsIsLeg12Matched; allElectronsIsLeg12Matched.clear();

      vector<bool> allElectronsIsTriggerMatched; allElectronsIsTriggerMatched.clear();
      vector<bool> idTightElectronsIsTriggerMatched; idTightElectronsIsTriggerMatched.clear();
      vector<bool> idLooseElectronsIsTriggerMatched; idLooseElectronsIsTriggerMatched.clear();
      vector<bool> isoElectronsIsTriggerMatched; isoElectronsIsTriggerMatched.clear();
      vector<bool> isoTightElectronsIsTriggerMatched; isoTightElectronsIsTriggerMatched.clear();

      vector<float> allElectronsIso; allElectronsIso.clear();
      vector<float> idTightElectronsIso; idTightElectronsIso.clear();
      vector<float> idLooseElectronsIso; idLooseElectronsIso.clear();
      vector<float> isoElectronsIso; isoElectronsIso.clear();
      vector<float> isoTightElectronsIso; isoTightElectronsIso.clear();

      vector<bool> isElectronPassedIdIso; isElectronPassedIdIso.clear();

      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {

	bool electronTriggerMatch = false;
	bool electronEle23Match = false;
	bool electronEle12Match = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.electron_eta[im],analysisTree.electron_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nElectronFilter] && 
	      analysisTree.electron_pt[im]>ptElectronHighCut && 
	      fabs(analysisTree.electron_eta[im])<etaElectronHighCut) { // Electron Leg of single electron trigger
	    electronTriggerMatch = true;
	  }
	  if (analysisTree.trigobject_filters[iT][nEle23Filter])
	    electronEle23Match = true;
	  if (analysisTree.trigobject_filters[iT][nEle12Filter])
	    electronEle12Match = true;
	    
	}

        allElectrons.push_back(im);
	allElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	allElectronsIsLeg23Matched.push_back(electronEle23Match);
	allElectronsIsLeg12Matched.push_back(electronEle12Match);
	bool isPassed = true;

	if (analysisTree.electron_pt[im]<ptElectronLowCut) isPassed = false;
	if (fabs(analysisTree.electron_eta[im])>etaElectronCut) isPassed = false;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) isPassed = false ;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut)  isPassed = false;
	bool electronMvaIdTightX = electronMvaIdWP80(analysisTree.electron_pt[im],
						     analysisTree.electron_superclusterEta[im],
						     analysisTree.electron_mva_id_nontrigPhys14[im]) && 
	  analysisTree.electron_pass_conversion[im] &&
	  analysisTree.electron_nmissinginnerhits[im] <= 1;
	bool electronMvaIdLooseX = electronMvaIdWP90(analysisTree.electron_pt[im],
						     analysisTree.electron_superclusterEta[im],
						     analysisTree.electron_mva_id_nontrigPhys14[im]) &&
	  analysisTree.electron_pass_conversion[im] &&
	  analysisTree.electron_nmissinginnerhits[im] <= 1;

	// isolation
	float absIso = analysisTree.electron_r03_sumChargedHadronPt[im];
	float neutralIso = 
	  analysisTree.electron_r03_sumNeutralHadronEt[im] + 
	  analysisTree.electron_r03_sumPhotonEt[im] - 
	  0.5*analysisTree.electron_r03_sumPUPt[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	absIso += neutralIso;
	float relIso = absIso/analysisTree.electron_pt[im];

	if (electronMvaIdLooseX && isPassed) {
	  idLooseElectrons.push_back(im);
	  idLooseElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	  idLooseElectronsIso.push_back(relIso);
	}
	if (electronMvaIdTightX && isPassed) {
	  idTightElectrons.push_back(im);
	  idTightElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	  idTightElectronsIso.push_back(relIso);
	}
	if (relIso<isoElectronCut && electronMvaIdLooseX && isPassed) {
	  isoElectrons.push_back(im);
	  isoElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	  isoElectronsIso.push_back(relIso);
	}
	if (relIso<isoElectronCut && electronMvaIdTightX && isPassed) {
	  isoTightElectrons.push_back(im);
	  isoTightElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	  isoTightElectronsIso.push_back(relIso);
	}
	isPassed = isPassed && electronMvaIdTightX && relIso<isoElectronCut;
	isElectronPassedIdIso.push_back(isPassed);
      }

      // end of electron selection

      // std::cout << "allElectrons : " << allElectrons.size() << std::endl;
      // std::cout << "idElectrons  : " << idElectrons.size() << std::endl;
      // std::cout << "isoElectrons : " << isoElectrons.size() << std::endl;

      bool isAllElectronsPair = false;
      if (allElectrons.size()>1) {
	//	std::cout << "allElectrons : " << allElectrons.size() << std::endl;
	for (unsigned int im1=0; im1<allElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<allElectrons.size(); ++im2) {
	    unsigned int index1 = allElectrons[im1];
	    unsigned int index2 = allElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    bool isTriggerMatched = allElectronsIsTriggerMatched[im1] || allElectronsIsTriggerMatched[im2]; 
	    if (q1*q2<0 && isTriggerMatched) {
	      isAllElectronsPair = true;
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[index2],
							  analysisTree.electron_py[index2],
							  analysisTree.electron_pz[index2],
							  electronMass);
	      TLorentzVector dielectron = electron1 + electron2;
	      float mass = dielectron.M();
	      JPsiMassAllElectronsH->Fill(mass,weight);
	      YpsilonMassAllElectronsH->Fill(mass,weight);
	      ZMassAllElectronsH->Fill(mass,weight);
              float ptLeadingE = analysisTree.electron_pt[index1];
	      float etaLeadingE = analysisTree.electron_eta[index1];
	      float ptTrailingE = analysisTree.electron_pt[index2];
	      float etaTrailingE = analysisTree.electron_eta[index2];
	      if (ptTrailingE>ptLeadingE) {
		float temp = ptLeadingE;
		ptLeadingE = ptTrailingE;
		ptTrailingE = temp;
		temp = etaLeadingE;
		etaLeadingE = etaTrailingE;
		etaTrailingE = temp;
	      }
	      if (mass>60&&mass<120) {
		ptLeadingAllElectronsH->Fill(ptLeadingE,weight);
		etaLeadingAllElectronsH->Fill(etaLeadingE,weight);
		ptTrailingAllElectronsH->Fill(ptTrailingE,weight);
		etaTrailingAllElectronsH->Fill(etaTrailingE,weight);
	      }
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassAllElectronsDRCutH->Fill(mass,weight);
		YpsilonMassAllElectronsDRCutH->Fill(mass,weight);
		ZMassAllElectronsDRCutH->Fill(mass,weight);
	      }
	    }
	  }
	}
      }

      bool isIdLooseElectronsPair = false;
      if (idLooseElectrons.size()>1) {
	for (unsigned int im1=0; im1<idLooseElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<idLooseElectrons.size(); ++im2) {
	    unsigned int index1 = idLooseElectrons[im1];
	    unsigned int index2 = idLooseElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
            bool isTriggerMatched = allElectronsIsTriggerMatched[im1] || allElectronsIsTriggerMatched[im2]; 
	    if (q1*q2<0 && isTriggerMatched) {
	      isIdLooseElectronsPair = true;
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[index2],
							  analysisTree.electron_py[index2],
							  analysisTree.electron_pz[index2],
							  electronMass);
	      TLorentzVector dielectron = electron1 + electron2;
	      float mass = dielectron.M();
	      //	      std::cout << "Mass = " << mass << std::endl;
	      JPsiMassIdLooseElectronsH->Fill(mass,weight);
	      YpsilonMassIdLooseElectronsH->Fill(mass,weight);
	      ZMassIdLooseElectronsH->Fill(mass,weight);
	      float ptLeadingE = analysisTree.electron_pt[index1];
	      float etaLeadingE = analysisTree.electron_eta[index1];
	      float ptTrailingE = analysisTree.electron_pt[index2];
	      float etaTrailingE = analysisTree.electron_eta[index2];
	      if (ptTrailingE>ptLeadingE) {
		float temp = ptLeadingE;
		ptLeadingE = ptTrailingE;
		ptTrailingE = temp;
		temp = etaLeadingE;
		etaLeadingE = etaTrailingE;
		etaTrailingE = temp;
	      }
	      if (mass>60&&mass<120) {
		ptLeadingIdLooseElectronH->Fill(ptLeadingE,weight);
		etaLeadingIdLooseElectronH->Fill(etaLeadingE,weight);
		ptTrailingIdLooseElectronH->Fill(ptTrailingE,weight);
		etaTrailingIdLooseElectronH->Fill(etaTrailingE,weight);
	      }
	      
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassIdLooseElectronsDRCutH->Fill(mass,weight);
		YpsilonMassIdLooseElectronsDRCutH->Fill(mass,weight);
		ZMassIdLooseElectronsDRCutH->Fill(mass,weight);
		float ptLeadingE = analysisTree.electron_pt[index1];
		float etaLeadingE = analysisTree.electron_eta[index1];
		float ptTrailingE = analysisTree.electron_pt[index2];
		float etaTrailingE = analysisTree.electron_eta[index2];
		if (ptTrailingE>ptLeadingE) {
		  float temp = ptLeadingE;
		  ptLeadingE = ptTrailingE;
		  ptTrailingE = temp;
		  temp = etaLeadingE;
		  etaLeadingE = etaTrailingE;
		  etaTrailingE = temp;
		}
          	if (mass>60&&mass<120) {
		  ptLeadingIdLooseElectronDRCutH->Fill(ptLeadingE,weight);
		  etaLeadingIdLooseElectronDRCutH->Fill(etaLeadingE,weight);
		  ptTrailingIdLooseElectronDRCutH ->Fill(ptTrailingE,weight);
		  etaTrailingIdLooseElectronDRCutH->Fill(etaTrailingE,weight);
		}
	      }
	    }
	  }
	}
      }
      
      bool isIdTightElectronsPair = false;
      if (idTightElectrons.size()>1) {
	for (unsigned int im1=0; im1<idTightElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<idTightElectrons.size(); ++im2) {
	    unsigned int index1 = idTightElectrons[im1];
	    unsigned int index2 = idTightElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
            bool isTriggerMatched = allElectronsIsTriggerMatched[im1] || allElectronsIsTriggerMatched[im2]; 
	    if (q1*q2<0 && isTriggerMatched) {
	      isIdTightElectronsPair = true;
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[index2],
							  analysisTree.electron_py[index2],
							  analysisTree.electron_pz[index2],
							  electronMass);
	      TLorentzVector dielectron = electron1 + electron2;
	      float mass = dielectron.M();
	      //	      std::cout << "Mass = " << mass << std::endl;
	      JPsiMassIdTightElectronsH->Fill(mass,weight);
	      YpsilonMassIdTightElectronsH->Fill(mass,weight);
	      ZMassIdTightElectronsH->Fill(mass,weight);
              float ptLeadingE = analysisTree.electron_pt[index1];
	      float etaLeadingE = analysisTree.electron_eta[index1];
	      float ptTrailingE = analysisTree.electron_pt[index2];
	      float etaTrailingE = analysisTree.electron_eta[index2];
	      if (ptTrailingE>ptLeadingE) {
		float temp = ptLeadingE;
		ptLeadingE = ptTrailingE;
		ptTrailingE = temp;
		temp = etaLeadingE;
		etaLeadingE = etaTrailingE;
		etaTrailingE = temp;
	      }
	      if (mass>60&&mass<120) {
		ptLeadingIdTightElectrons->Fill(ptLeadingE,weight);
		etaLeadingIdTightElectrons->Fill(etaLeadingE,weight);
		ptTrailingIdTightElectrons ->Fill(ptTrailingE,weight);
		etaTrailingIdTightElectrons->Fill(etaTrailingE,weight);
	      }
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassIdTightElectronsDRCutH->Fill(mass,weight);
		YpsilonMassIdTightElectronsDRCutH->Fill(mass,weight);
		ZMassIdTightElectronsDRCutH->Fill(mass,weight);
		float ptLeadingE = analysisTree.electron_pt[index1];
		float etaLeadingE = analysisTree.electron_eta[index1];
		float ptTrailingE = analysisTree.electron_pt[index2];
		float etaTrailingE = analysisTree.electron_eta[index2];
		if (ptTrailingE>ptLeadingE) {
		  float temp = ptLeadingE;
		  ptLeadingE = ptTrailingE;
		  ptTrailingE = temp;
		  temp = etaLeadingE;
		  etaLeadingE = etaTrailingE;
		  etaTrailingE = temp;
		}
         	if (mass>60&&mass<120) {
		  ptLeadingIdTightElectronsDRCutH->Fill(ptLeadingE,weight);
		  etaLeadingIdTightElectronsDRCutH->Fill(etaLeadingE,weight);
		  ptTrailingIdTightElectronsDRCutH ->Fill(ptTrailingE,weight);
		  etaTrailingIdTightElectronsDRCutH->Fill(etaTrailingE,weight);
		}
	      }
	    }
	  }
	}
      }     

      bool isIsoElectronsPair = false;
      if (isoElectrons.size()>1) {
	for (unsigned int im1=0; im1<isoElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<isoElectrons.size(); ++im2) {
	    unsigned int index1 = isoElectrons[im1];
	    unsigned int index2 = isoElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    bool isTriggerMatched = isoElectronsIsTriggerMatched[im1] || isoElectronsIsTriggerMatched[im2]; 
	    if (q1*q2<0  && isTriggerMatched){
	      isIsoElectronsPair = true;
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[index2],
							  analysisTree.electron_py[index2],
							  analysisTree.electron_pz[index2],
							  electronMass);
	      TLorentzVector dielectron = electron1 + electron2;
	      float mass = dielectron.M();
	      ZMassIsoElectronsH->Fill(mass,weight);
	      
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		ZMassIsoElectronsDRCutH->Fill(mass,weight);
	      }
	    }
	  }
	}
	  }
      
      bool isIsoTightElectronsPair = false;
      if (isoTightElectrons.size()>0) {
	unsigned int iE1 = 0;
	unsigned int iE2 = 0;
	bool isPairFound = false;
	float minIso = 999999;
	for (unsigned int im1=0; im1<isoTightElectrons.size(); ++im1) {
	  unsigned int index1 = isoTightElectrons[im1];
	  bool isTriggerMatched = isoTightElectronsIsTriggerMatched[im1];
	  if (isTriggerMatched &&
	      analysisTree.electron_pt[index1]>ptElectronHighCut&&
	      fabs(analysisTree.electron_eta[index1])<etaElectronHighCut) {
	    for (unsigned int iE=0; iE<allElectrons.size(); ++iE) {
	      unsigned int indexProbe = allElectrons[iE];
	      if (index1==indexProbe) continue;
	      float q1 = analysisTree.electron_charge[index1];
	      float q2 = analysisTree.electron_charge[indexProbe];
	      if (q1*q2>0) continue;
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[indexProbe],analysisTree.electron_phi[indexProbe]);
	      if (dR<dRleptonsCut) continue;
	      float ptProbe = TMath::Min(float(analysisTree.electron_pt[indexProbe]),float(999.9));
	      float absEtaProbe = fabs(analysisTree.electron_eta[indexProbe]);
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      int ptBinTrig = binNumber(ptProbe,nPtBinsTrig,ptBinsTrig);
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[indexProbe],
							  analysisTree.electron_py[indexProbe],
							  analysisTree.electron_pz[indexProbe],
							  electronMass);
	      float mass = (electron1+electron2).M();
	      if (isElectronPassedIdIso[iE]) {
		ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight);
		// ele23 filter
		if (allElectronsIsLeg23Matched[iE])
		  ZMassEle23EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		else 
		  ZMassEle23EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		if (ptProbe>24&&absEtaProbe<2.5) {
		  if (absEtaProbe<1.48) {
		    if (allElectronsIsLeg23Matched[iE])
		      ZMassEle23BarrelPass->Fill(mass,weight);
		    else 
		      ZMassEle23BarrelFail->Fill(mass,weight);
		  }
		  else {
		    if (allElectronsIsLeg23Matched[iE])
                      ZMassEle23EndcapPass->Fill(mass,weight);
                    else 
                      ZMassEle23EndcapFail->Fill(mass,weight);
		  }
		}
		// ele12 filter
		if (allElectronsIsLeg12Matched[iE])
		  ZMassEle12EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		else 
		  ZMassEle12EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		if (ptProbe>13&&absEtaProbe<2.5) {
		  if (absEtaProbe<1.48) {
		    if (allElectronsIsLeg12Matched[iE])
		      ZMassEle12BarrelPass->Fill(mass,weight);
		    else 
		      ZMassEle12BarrelFail->Fill(mass,weight);
		  }
		  else {
		    if (allElectronsIsLeg12Matched[iE])
                      ZMassEle12EndcapPass->Fill(mass,weight);
                    else 
                      ZMassEle12EndcapFail->Fill(mass,weight);
		  }
		}
	      }
	      else
		ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
	    }
	  }
	  for (unsigned int im2=im1+1; im2<isoTightElectrons.size(); ++im2) {
	    unsigned int index2 = isoTightElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    bool isTriggerMatched1 = isoTightElectronsIsTriggerMatched[im1] && 
	      analysisTree.electron_pt[index1]>ptElectronHighCut &&
	      fabs(analysisTree.electron_eta[index1])<etaElectronHighCut;
	    bool isTriggerMatched2 = isoTightElectronsIsTriggerMatched[im2] && 
	      analysisTree.electron_pt[index2]>ptElectronHighCut &&
	      fabs(analysisTree.electron_eta[index2])<etaElectronHighCut; 
	    bool isTriggerMatched = (isTriggerMatched1 && analysisTree.electron_pt[index2]>ptElectronLowCut) ||
	      (isTriggerMatched2 && analysisTree.electron_pt[index1]>ptElectronLowCut);
	    bool sign = q1*q2>0;
	    float deltaREE = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				    analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	    if (oppositeSign)
	      sign = q1*q2 < 0;
	    if (sign && isTriggerMatched && deltaREE>dRleptonsCut) {
	      float isoSum = isoTightElectronsIso[im1] + isoTightElectronsIso[im2];
	      if (isoSum<minIso) {
		minIso = isoSum;
		isPairFound = true;
		iE1 = index1;
		iE2 = index2;
	      }

	      isIsoTightElectronsPair = true;
	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[index2],
							  analysisTree.electron_py[index2],
							  analysisTree.electron_pz[index2],
							  electronMass);
	      TLorentzVector dielectron = electron1 + electron2;
	      float mass = dielectron.M();
	      ZMassIsoTightElectronsH->Fill(mass,weight);
	      float ptLeadingE = analysisTree.electron_pt[index1];
	      float etaLeadingE = analysisTree.electron_eta[index1];
	      float ptTrailingE = analysisTree.electron_pt[index2];
	      float etaTrailingE = analysisTree.electron_eta[index2];
	      if (ptTrailingE>ptLeadingE) {
		float temp = ptLeadingE;
		ptLeadingE = ptTrailingE;
		ptTrailingE = temp;
		temp = etaLeadingE;
		etaLeadingE = etaTrailingE;
		etaTrailingE = temp;
	      }
	      if (mass>60&&mass<120) {
		ptLeadingIsoTightElectronsH->Fill(ptLeadingE,weight);
		etaLeadingIsoTightElectronsH->Fill(etaLeadingE,weight);
		ptTrailingIsoTightElectronsH ->Fill(ptTrailingE,weight);
		etaTrailingIsoTightElectronsH ->Fill(etaTrailingE,weight);
	      }

	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		ZMassIsoTightElectronsDRCutH->Fill(mass,weight);
		float ptLeadingE = analysisTree.electron_pt[index1];
		float etaLeadingE = analysisTree.electron_eta[index1];
		float ptTrailingE = analysisTree.electron_pt[index2];
		float etaTrailingE = analysisTree.electron_eta[index2];
		if (ptTrailingE>ptLeadingE) {
		  float temp = ptLeadingE;
		  ptLeadingE = ptTrailingE;
		  ptTrailingE = temp;
		  temp = etaLeadingE;
		  etaLeadingE = etaTrailingE;
		  etaTrailingE = temp;
		}
		ptLeadingIsoTightElectronsDRCutH ->Fill(ptLeadingE,weight);
		etaLeadingIsoTightElectronsDRCutH->Fill(etaLeadingE,weight);
		ptTrailingIsoTightElectronsDRCutH ->Fill(ptTrailingE,weight);
		etaTrailingIsoTightElectronsDRCutH->Fill(etaTrailingE,weight);
	      }
	      
	    }
	  }
	}
	if (isPairFound) {
	  TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[iE1],
						      analysisTree.electron_py[iE1],
						      analysisTree.electron_pz[iE1],
						      electronMass);
	  TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[iE2],
						      analysisTree.electron_py[iE2],
						      analysisTree.electron_pz[iE2],
						      electronMass);
	  float mass = (electron1+electron2).M();

          float ptLeadingE = analysisTree.electron_pt[iE1];
	  float etaLeadingE = analysisTree.electron_eta[iE1];
	  float ptTrailingE = analysisTree.electron_pt[iE2];
	  float etaTrailingE = analysisTree.electron_eta[iE2];
	  if (ptTrailingE>ptLeadingE) {
	    float temp = ptLeadingE;
	    ptLeadingE = ptTrailingE;
	    ptTrailingE = temp;
	    temp = etaLeadingE;
	    etaLeadingE = etaTrailingE;
	    etaTrailingE = temp;
	  } 
	  //reweighting for efficiency
	  
	  float ptLeadingEX = TMath::Min(float(59),ptLeadingE);
	  float ptTrailingEX = TMath::Min(float(59),ptTrailingE);
	  int ptBinLeading = binNumber(ptLeadingEX,nPtBins,ptBins);
	  int ptBinTrailing = binNumber(ptTrailingEX,nPtBins,ptBins); 


	  // int binDataEff = hEffBarrelData->FindBin(ptLeadingE);
	  // int binMCEff = hEffBarrelMC->FindBin(ptLeadingE);
	  // int binDataEff_1 = hEffEndcapData->FindBin(ptLeadingE);
	  // int binMCEff_1 = hEffEndcapMC->FindBin(ptLeadingE);

	  if (!isData && applyLeptonSF) {
	    float dataLeading = 1;
	    float mcLeading = 1;
	    float dataTrailing = 1;
	    float mcTrailing = 1;
	    if(abs(etaLeadingE) < 1.48){
	      dataLeading = double(dataEffBarrel[ptBinLeading]);
	      mcLeading = double(mcEffBarrel[ptBinLeading]);
	    }
	    else {
	      dataLeading = double(dataEffEndcap[ptBinLeading]);
	      mcLeading = double(mcEffEndcap[ptBinLeading]);
	    }
	    if(abs(etaTrailingE) < 1.48){
	      dataTrailing = double(dataEffBarrel[ptBinTrailing]);
	      mcTrailing = double(mcEffBarrel[ptBinTrailing]);
	    }
	    else {
	      dataTrailing = double(dataEffEndcap[ptBinTrailing]);
	      mcTrailing= double(mcEffEndcap[ptBinTrailing]);
	    }
	    float_t wLeading  = dataLeading/mcLeading;
	    float_t wTrailing = dataTrailing/mcTrailing ;
	    
	    //	    cout << "El SF : lead = " << wLeading << "   trail = " << wTrailing << endl;
	    weight = weight*wLeading*wTrailing;
	  }

	  ZMassIsoTightElectronsUniqueH->Fill(mass,weight);
	  for (int iScale=0; iScale<21; ++ iScale) {
	    float scaleFactor = 0.98 + 0.002*float(iScale);
	    ZMassIsoTightElectronsUniqueScaleH[iScale]->Fill(mass*scaleFactor,weight);
	  }


	  if (mass>30) {
	    ptLeadingIsoTightElectronsUniqueH ->Fill(ptLeadingE,weight);
	    etaLeadingIsoTightElectronsUniqueH->Fill(etaLeadingE,weight);
	    ptTrailingIsoTightElectronsUniqueH ->Fill(ptTrailingE,weight);
	    etaTrailingIsoTightElectronsUniqueH->Fill(etaTrailingE,weight);
	    float pfMetX = analysisTree.pfmet_ex;
	    float pfMetY = analysisTree.pfmet_ey;
	    float pfMet = sqrt(pfMetX*pfMetX+pfMetY*pfMetY);
	    NumberOfVerticesH->Fill(float(analysisTree.primvertex_count),weight);
	    MetH->Fill(pfMet,weight);
	    for (int iScale=0; iScale<21; ++ iScale) {
	      float scaleFactor = 0.9 + 0.01*float(iScale);
	      MetScaleH[iScale]->Fill(pfMet*scaleFactor,weight);
	    }
	  }
	}
      }
      
      if (isAllElectronsPair) selEventsAllElectrons++;
      if (isIdLooseElectronsPair)  selEventsIdLooseElectrons++;
      if (isIdTightElectronsPair)  selEventsIdTightElectrons++;
      if (isIsoElectronsPair) selEventsIsoElectrons++;
      if (isIsoTightElectronsPair) selEventsIsoTightElectrons++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                              = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                            = " << nEvents << std::endl;
  std::cout << "Total number of selected events (electron pairs)          = " << selEventsAllElectrons << std::endl;
  std::cout << "Total number of selected events (idLoose electron pairs)  = " << selEventsIdLooseElectrons << std::endl;
  std::cout << "Total number of selected events (idTight electron pairs)  = " << selEventsIdTightElectrons << std::endl;
  std::cout << "Total number of selected events (iso electron pairs)      = " << selEventsIsoElectrons << std::endl;
  std::cout << "Total number of selected events (iso tight electron pairs)= " << selEventsIsoTightElectrons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  file->Write();
  file->Close();
  delete file;
  
  
}
