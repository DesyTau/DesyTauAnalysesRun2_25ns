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

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

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

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on probed electrons
  const float ptElectronCut   = cfg.get<float>("ptElectronCut");
  const float ptElectronHighCut   = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");
  const bool applyElectronIso    = cfg.get<bool>("ApplyElectronIso");
  const bool applyElectronIPcuts = cfg.get<bool>("ApplyElectronIPcuts");
  
  // kinematic cuts on tag electrons
  const float ptElectronTagCut  = cfg.get<float>("ptElectronTagCut");
  const float etaElectronTagCut = cfg.get<float>("etaElectronTagCut");
  const float dxyElectronTagCut      = cfg.get<float>("dxyElectronTagCut");
  const float dzElectronTagCut      = cfg.get<float>("dzElectronTagCut");
  const float isoElectronLowTagCut  = cfg.get<float>("isoElectronLowTagCut");
  const float isoElectronHighTagCut = cfg.get<float>("isoElectronHighTagCut");
  const bool applyElectronIdTag     = cfg.get<bool>("ApplyElectronIdTag");
  const bool applyElectronIsoTag    = cfg.get<bool>("ApplyElectronIsoTag");
  const bool applyElectronIPcutsTag = cfg.get<bool>("ApplyElectronIPcutsTag");
  const bool applyElectronTriggerMatchTag = cfg.get<bool>("ApplyElectronTriggerMatchTag");

  // kinematic cuts on probed muons
  const float ptMuonCut       = cfg.get<float>("ptMuonCut");
  const float ptMuonHighCut   = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const int applyMuonId      = cfg.get<int>("ApplyMuonId");
  const bool applyMuonIso    = cfg.get<bool>("ApplyMuonIso");
  const bool applyMuonIPcuts = cfg.get<bool>("ApplyMuonIPcuts");

  // kinematic cuts on tag muons
  const float ptMuonTagCut  = cfg.get<float>("ptMuonTagCut");
  const float etaMuonTagCut  = cfg.get<float>("etaMuonTagCut");
  const float dxyMuonTagCut     = cfg.get<float>("dxyMuonTagCut");
  const float dzMuonTagCut      = cfg.get<float>("dzMuonTagCut");
  const float isoMuonLowTagCut  = cfg.get<float>("isoMuonLowTagCut");
  const float isoMuonHighTagCut = cfg.get<float>("isoMuonHighTagCut");
  const int applyMuonIdTag      = cfg.get<int>("ApplyMuonIdTag");
  const bool applyMuonIsoTag    = cfg.get<bool>("ApplyMuonIsoTag");
  const bool applyMuonIPcutsTag = cfg.get<bool>("ApplyMuonIPcutsTag");
  const bool applyMuonTriggerMatchTag = cfg.get<bool>("ApplyMuonTriggerMatchTag");

  // topological cuts
  const float dRLeptonsCut   = cfg.get<float>("dRLeptonsCut");
  const float dPhiLeptonsLowCut  = cfg.get<float>("dPhiLeptonsLowCut");
  const float dPhiLeptonsHighCut = cfg.get<float>("dPhiLeptonsHighCut");
  const float dZetaCut       = cfg.get<float>("dZetaCut");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");

  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const float DRTrigMatch      = cfg.get<float>("DRTrigMatch");
  const bool singleMuTrigger   = cfg.get<bool>("SingleMuTrigger");

  const bool useGeneratorAsTag = cfg.get<bool>("UseGeneratorAsTag"); 
  const bool isTauAcceptance = cfg.get<bool>("IsTauAcceptance"); 

  // vertex cuts
  const float vertexNdofCut = cfg.get<float>("VertexNdofCut");
  const float vertexZCut    = cfg.get<float>("VertexZCut");
  const float vertexDCut    = cfg.get<float>("VertexDCut");

  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;
  bool isTTJets = false;
  if (TStrName.Contains("TTJets"))
    isTTJets = true;

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

  // vertex information
  TH1F * vertexNdofH = new TH1F("vertexNdofH","",100,-0.5,99.5);
  TH1F * vertexXH = new TH1F("vertexXH","",250,0.,0.5);
  TH1F * vertexYH = new TH1F("vertexYH","",250,0.0,0.5);
  TH1F * vertexZH = new TH1F("vertexZH","",200,-50.,50.);
  TH1F * vertexDH = new TH1F("vertexDH","",100,0.,5.);
  TH1F * vertexChi2H = new TH1F("vertexChi2H","",200,0,20);
  // ******************

  int nEtaBins = 3;
  float etaBins[4] = {-0.1, 0.8, 1.479, 5};

  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};

  TString EtaBinName[3] = {"0To0p8",
			   "0p8To1p5",
			   "1p5To2p3"};

  TString PtBinName[7] = {"pt10to15",
			  "pt15to20",
			  "pt20to30",
			  "pt30to40",
			  "pt40to50",
			  "pt50to70",
			  "pt70to100"};

  TString PassLevel[15] = {"True",//0
			   "TrueAll",//1
			   "TruePassed",//2
			   "RecoAll",//3    
			   "RecoPassed",//4
			   "TrigDenHigh",//5
			   "TrigNumHigh",//6
			   "TrigDenLow",//7
			   "TrigNumLow",//8
			   "TPAll",//9
			   "TPPassed",//10
			   "RecoAllTag",//11
			   "RecoPassedTag",//12 
			   "RecoPassedId",//13
			   "RecoPassedIdTag"//14
  };

  TString Type[2] = {"Prompt","NonPrompt"};
  TString LepQ[2] = {"Pos","Neg"};

  TString ptCoarseBins[2] = {"PtLess20","PtGreater20"};

  // histograms for lepton Id
  TH1F * MuonPtH[2][2][15][3];
  TH1F * ElectronPtH[2][2][15][3];

  TH1F * MuonNvertEtaH[2][2][15][3];
  TH1F * ElectronNvertEtaH[2][2][15][3];

  TH1F * MuonNvertPtH[2][2][15][2];
  TH1F * ElectronNvertPtH[2][2][15][2];

  TH1F * MuonNvertH[2][2][15];
  TH1F * ElectronNvertH[2][2][15];

  TH1F * MuonEtaH[2][2][15][2];
  TH1F * ElectronEtaH[2][2][15][2];

  TH1F * MuonChargedIsoH[2][2][15][3];
  TH1F * MuonNeutralIsoH[2][2][15][3];
  TH1F * MuonNeutralCorrIsoH[2][2][15][3];
  TH1F * MuonIsoH[2][2][15][3];
  TH1F * MuonRelIsoH[2][2][15][3];

  TH1F * ElectronChargedIsoH[2][2][15][3];
  TH1F * ElectronNeutralIsoH[2][2][15][3];
  TH1F * ElectronNeutralCorrIsoH[2][2][15][3];
  TH1F * ElectronIsoH[2][2][15][3];
  TH1F * ElectronRelIsoH[2][2][15][3];

  TH1F * MuonVertexDxH[2][2][15][3];
  TH1F * MuonVertexDyH[2][2][15][3];
  TH1F * MuonVertexDzH[2][2][15][3];
  
  TH1F * ElectronVertexDxH[2][2][15][3];
  TH1F * ElectronVertexDyH[2][2][15][3];
  TH1F * ElectronVertexDzH[2][2][15][3];

  TH2F * Mu23Ele12H = new TH2F("Mu23Ele12H","",2,-0.5,1.5,2,-0.5,1.5);
  TH2F * Mu8Ele23H = new TH2F("Mu8Ele23H","",2,-0.5,1.5,2,-0.5,1.5);

  TH1F * MuonPt23H       = new TH1F("MuonPt23H","",100,0,100);
  TH1F * MuonPt23PassedH = new TH1F("MuonPt23PassedH","",100,0,100);

  TH1F * MuonPt8H       = new TH1F("MuonPt8H","",100,0,100);
  TH1F * MuonPt8PassedH = new TH1F("MuonPt8PassedH","",100,0,100);

  TH1F * ElectronPt23H       = new TH1F("ElectronPt23H","",100,0,100);
  TH1F * ElectronPt23PassedH = new TH1F("ElectronPt23PassedH","",100,0,100);

  TH1F * ElectronPt12H       = new TH1F("ElectronPt12H","",100,0,100);
  TH1F * ElectronPt12PassedH = new TH1F("ElectronPt12PassedH","",100,0,100);


  TH1F * MuonPt23RecoH       = new TH1F("MuonPt23RecoH","",100,0,100);
  TH1F * MuonPt23RecoPassedH = new TH1F("MuonPt23RecoPassedH","",100,0,100);

  TH1F * MuonPt8RecoH       = new TH1F("MuonPt8RecoH","",100,0,100);
  TH1F * MuonPt8RecoPassedH = new TH1F("MuonPt8RecoPassedH","",100,0,100);

  TH1F * ElectronPt23RecoH       = new TH1F("ElectronPt23RecoH","",100,0,100);
  TH1F * ElectronPt23RecoPassedH = new TH1F("ElectronPt23RecoPassedH","",100,0,100);

  TH1F * ElectronPt12RecoH       = new TH1F("ElectronPt12RecoH","",100,0,100);
  TH1F * ElectronPt12RecoPassedH  = new TH1F("ElectronPt12RecoPassedH","",100,0,100);


  TH1F * MuonPt23RecoXH       = new TH1F("MuonPt23RecoXH","",100,0,100);
  TH1F * MuonPt23RecoPassedXH = new TH1F("MuonPt23RecoPassedXH","",100,0,100);

  TH1F * MuonPt8RecoXH       = new TH1F("MuonPt8RecoXH","",100,0,100);
  TH1F * MuonPt8RecoPassedXH = new TH1F("MuonPt8RecoPassedXH","",100,0,100);

  TH1F * ElectronPt23RecoXH       = new TH1F("ElectronPt23RecoXH","",100,0,100);
  TH1F * ElectronPt23RecoPassedXH = new TH1F("ElectronPt23RecoPassedXH","",100,0,100);

  TH1F * ElectronPt12RecoXH       = new TH1F("ElectronPt12RecoXH","",100,0,100);
  TH1F * ElectronPt12RecoPassedXH  = new TH1F("ElectronPt12RecoPassedXH","",100,0,100);


  TH1F * MuonPt23RecoYH       = new TH1F("MuonPt23RecoYH","",100,0,100);
  TH1F * MuonPt23RecoPassedYH = new TH1F("MuonPt23RecoPassedYH","",100,0,100);

  TH1F * MuonPt8RecoYH       = new TH1F("MuonPt8RecoYH","",100,0,100);
  TH1F * MuonPt8RecoPassedYH = new TH1F("MuonPt8RecoPassedYH","",100,0,100);

  TH1F * ElectronPt23RecoYH       = new TH1F("ElectronPt23RecoYH","",100,0,100);
  TH1F * ElectronPt23RecoPassedYH = new TH1F("ElectronPt23RecoPassedYH","",100,0,100);

  TH1F * ElectronPt12RecoYH       = new TH1F("ElectronPt12RecoYH","",100,0,100);
  TH1F * ElectronPt12RecoPassedYH  = new TH1F("ElectronPt12RecoPassedYH","",100,0,100);

  TH1F * ZToTauTauH = new TH1F("ZToTauTauH","",3,-0.5,2.5);
  TH1F * numPartonsH = new TH1F("numPartonsH","",10,-0.5,9.5);

  TH1F * dRMuonLeadingPartonH[2][2][2][3];
  TH1F * dRElectronLeadingPartonH[2][2][2][3];
  
  TH1F * dRMuonLeadingPartonTagH[2][2][2][3];
  TH1F * dRElectronLeadingPartonTagH[2][2][2][3];

  TH2F * dRMuonPartonDPhiTagH[2][2][2];
  TH2F * dRElectronPartonDPhiTagH[2][2][2];

  for (int iLepQ=0; iLepQ<2; ++iLepQ) {
    for (int iType=0; iType<2; ++iType) {
      for (int iPt=0; iPt<2; ++iPt) {
	dRMuonPartonDPhiTagH[iLepQ][iType][iPt] = new TH2F("dRMuonPartonDPhiTag"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt],"",50,0.,5,50,0,TMath::Pi());
	dRElectronPartonDPhiTagH[iLepQ][iType][iPt] = new TH2F("dRElectronPartonDPhiTag"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt],"",50,0.,5,50,0,TMath::Pi());


	for (int iEta=0; iEta<3; ++iEta) {
	  dRMuonLeadingPartonH[iLepQ][iType][iPt][iEta] = new TH1F("dRMuonLeadingParton"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt]+EtaBinName[iEta],"",50,0,5);
	  dRElectronLeadingPartonH[iLepQ][iType][iPt][iEta] = new TH1F("dRElectronLeadingParton"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt]+EtaBinName[iEta],"",50,0,5);
	  dRMuonLeadingPartonTagH[iLepQ][iType][iPt][iEta] = new TH1F("dRMuonLeadingPartonTag"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt]+EtaBinName[iEta],"",50,0,5);
	  dRElectronLeadingPartonTagH[iLepQ][iType][iPt][iEta] = new TH1F("dRElectronLeadingPartonTag"+LepQ[iLepQ]+Type[iType]+ptCoarseBins[iPt]+EtaBinName[iEta],"",50,0,5);

	}	
      }

      for (int iPassLevel=0; iPassLevel<15; ++iPassLevel) {

	MuonNvertH[iLepQ][iType][iPassLevel] = new TH1F("MuonNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel],"",50,-0.5,49.5);
	ElectronNvertH[iLepQ][iType][iPassLevel] = new TH1F("ElectronNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel],"",50,-0.5,49.5);

	for (int iPt=0; iPt<2; ++iPt) {

	  MuonEtaH[iLepQ][iType][iPassLevel][iPt] = new TH1F("MuonEta_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+ptCoarseBins[iPt],"",48,-2.4,2.4);
	  ElectronEtaH[iLepQ][iType][iPassLevel][iPt] = new TH1F("ElectronEta_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+ptCoarseBins[iPt],"",48,-2.4,2.4);

	  MuonNvertPtH[iLepQ][iType][iPassLevel][iPt] = new TH1F("MuonNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+ptCoarseBins[iPt],"",50,-0.5,49.5);
	  ElectronNvertPtH[iLepQ][iType][iPassLevel][iPt] = new TH1F("ElectronNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+ptCoarseBins[iPt],"",50,-0.5,49.5);

	}
	for (int iEta=0; iEta<3; ++iEta) {

	  MuonPtH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonPt_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0.,100);
	  ElectronPtH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronPt_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0.,100);

	  MuonNvertEtaH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",50,-0.5,49.5);
	  ElectronNvertEtaH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronNvert_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",50,-0.5,49.5);

	  MuonChargedIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonChargedIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  MuonNeutralIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonNeutralIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  MuonNeutralCorrIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonNeutralCorrIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  MuonIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  MuonRelIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonRelIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",40,0,2);

	  ElectronChargedIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronChargedIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  ElectronNeutralIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronNeutralIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  ElectronNeutralCorrIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronNeutralCorrIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  ElectronIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0,20);
	  ElectronRelIsoH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronRelIso_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",40,0,2);

	  MuonVertexDxH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonVertexDx_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-1,1);
	  MuonVertexDyH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonVertexDy_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-1,1);
	  MuonVertexDzH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonVertexDz_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-20,20);
	  
	  ElectronVertexDxH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronVertexDx_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-1,1);
	  ElectronVertexDyH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronVertexDy_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-1,1);
	  ElectronVertexDzH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronVertexDz_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,-20,20);

	}
      }
    }
  }

  TString passName[2] = {"Pass","Fail"};

  TH1F * massMuMuId[2][3][7][2];
  TH1F * massEEId[2][3][7][2];
  
  for (int iQ=0; iQ<2; ++iQ) {
    for (int iEta=0; iEta<3; ++iEta) {
      for (int iPt=0; iPt<7; ++iPt) {
	for (int iPass=0; iPass<2; ++iPass) {
	  massMuMuId[iQ][iEta][iPt][iPass] = new TH1F("massMuMuId_"+LepQ[iQ]+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass],"",60,60,120);
	  massEEId[iQ][iEta][iPt][iPass] = new TH1F("massEEId_"+LepQ[iQ]+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass],"",60,60,120);
	}
      }
    }
  }

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

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

      // **** vertex quality cuts
      vertexNdofH->Fill(analysisTree.primvertex_ndof,weight);
      if (analysisTree.primvertex_ndof<vertexNdofCut) continue;

      vertexZH->Fill(analysisTree.primvertex_z,weight);
      if (fabs(analysisTree.primvertex_z)>vertexZCut) continue;

      vertexXH->Fill(analysisTree.primvertex_x,weight);
      vertexYH->Fill(analysisTree.primvertex_y,weight);
      float vertexD = TMath::Sqrt(analysisTree.primvertex_x*analysisTree.primvertex_x+
				  analysisTree.primvertex_y*analysisTree.primvertex_y);
      vertexDH->Fill(vertexD,weight);
      if (fabs(vertexD)>vertexDCut) continue;

      float normChi2 = analysisTree.primvertex_chi2/analysisTree.primvertex_ndof;
      if (normChi2>20) normChi2 = 19.99;
      vertexChi2H->Fill(normChi2,weight);

      // **** end vertex quality cuts 


      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      vector<unsigned int> indexZ; indexZ.clear();
      vector<unsigned int> indexW; indexW.clear();
      vector<unsigned int> indexMu; indexMu.clear();
      vector<unsigned int> indexE; indexE.clear();
      vector<unsigned int> indexTau; indexTau.clear();
      vector<int> qTau; qTau.clear();
      vector<unsigned int> indexPartons; indexPartons.clear();

      float genVertX = 0;
      float genVertY = 0;
      float genVertZ = 0;

      //      bool count_partons = false;

      for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

       	float pxGen = analysisTree.genparticles_px[igen];
       	float pyGen = analysisTree.genparticles_py[igen];
       	float pzGen = analysisTree.genparticles_pz[igen];
       	float etaGen = PtoEta(pxGen,pyGen,pzGen);
       	float ptGen  = PtoPt(pxGen,pyGen);
	float phiGen = TMath::ATan2(pyGen,pxGen);

	if (analysisTree.genparticles_status[igen]==23) {
	  int pdgP = analysisTree.genparticles_pdgid[igen];
	  if (TMath::Abs(pdgP)<6||pdgP==21)
	    indexPartons.push_back(igen);

	}

       	if (fabs(analysisTree.genparticles_pdgid[igen])==23 && 
	    (analysisTree.genparticles_status[igen]==62||analysisTree.genparticles_status[igen]==52)) { 
       	  indexZ.push_back(igen);
	  //	  std::cout << "Z boson :  pt = " << ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //	     	    << "   status = " << analysisTree.genparticles_status[igen] <<  std::endl; 
	  genVertX = analysisTree.genparticles_vx[igen];
	  genVertY = analysisTree.genparticles_vy[igen];
	  genVertZ = analysisTree.genparticles_vz[igen];
	  //	  count_partons = true;
	}
       	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && 
	    (analysisTree.genparticles_status[igen]==62||analysisTree.genparticles_status[igen]==52)) { 
       	  indexW.push_back(igen);
	  // std::cout << "W boson :  pt = " << ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  // 	    << "   status = " << analysisTree.genparticles_status[igen] <<  std::endl; 
	}
       	if (fabs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_status[igen]==1) {
	  // std::cout << "muon Id = " <<  analysisTree.genparticles_pdgid[igen]
	  //  	    << "   pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //  	    << "   decay = " << analysisTree.genparticles_info[igen] << std::endl; 
	  if (analysisTree.genparticles_info[igen]==1||
	      analysisTree.genparticles_info[igen]==5||
	      analysisTree.genparticles_info[igen]==2||
	      analysisTree.genparticles_info[igen]==6) 
	    indexMu.push_back(igen);
       	}
       	if (fabs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_status[igen]==1) {
	  // std::cout << "electron Id = " <<  analysisTree.genparticles_pdgid[igen]
	  //      	    << "   pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //      	    << "   decay = " << analysisTree.genparticles_info[igen] << std::endl; 
	  if (analysisTree.genparticles_info[igen]==1||
	      analysisTree.genparticles_info[igen]==5||
	      analysisTree.genparticles_info[igen]==2||
	      analysisTree.genparticles_info[igen]==6)
	    indexE.push_back(igen);
       	}
      }

      float ptPartonMax = -1;
      int nPartonLeading = -1;
      // std::cout << "Number of partons = " << indexPartons.size() << std::endl;
      for (unsigned int iP=0; iP<indexPartons.size(); ++iP) {
	unsigned int index = indexPartons[iP];
	int pdgP = analysisTree.genparticles_pdgid[iP];
	float pxGen = analysisTree.genparticles_px[iP];
        float pyGen = analysisTree.genparticles_py[iP];
        float pzGen = analysisTree.genparticles_pz[iP];
        float etaGen = PtoEta(pxGen,pyGen,pzGen);
        float phiGen = TMath::ATan2(pyGen,pxGen);
        float ptGen  = PtoPt(pxGen,pyGen);
	// std::cout << iP << " " << pdgP 
	// 	  << "   pt = " << ptGen
	// 	  << "   eta = " << etaGen
	// 	  << "   phi = " << phiGen << std::endl;
	if (ptGen>ptPartonMax) {
	  ptPartonMax = ptGen;
	  nPartonLeading = iP;
	}
      }
      numPartonsH->Fill(indexPartons.size(),weight);


      // std::cout << std::endl;
      // continue;

      float dVertX = analysisTree.primvertex_x - genVertX;
      float dVertY = analysisTree.primvertex_y - genVertY;
      float dVertZ = analysisTree.primvertex_z - genVertZ;

      bool isZToEE = false;
      bool isZToEEopen = false;
      bool isZToMM = false;
      bool isZToMMopen = false;

      bool posMuonInAcceptance = false;
      bool negMuonInAcceptance = false;

      bool posElectronInAcceptance = false;
      bool negElectronInAcceptance = false;

      if (indexMu.size()==2) {
	unsigned int index1 = indexMu.at(0);
	unsigned int index2 = indexMu.at(1);
	int q1 = analysisTree.genparticles_pdgid[index1];
	int q2 = analysisTree.genparticles_pdgid[index2];
	if ( (q1*q2<0) && 
	     (analysisTree.genparticles_info[index1]==1) &&
	     (analysisTree.genparticles_info[index2]==1) ) { 
	  isZToMM = true;
	  
	  float pxGen1 = analysisTree.genparticles_px[index1];
	  float pyGen1 = analysisTree.genparticles_py[index1];
	  float pzGen1 = analysisTree.genparticles_pz[index1];
	  float ptGen1 = PtoPt(pxGen1,pyGen1);
	  float etaGen1 = PtoEta(pxGen1,pyGen1,pzGen1);
	  float phiGen1 = TMath::ATan2(pyGen1,pxGen1);
	  
	  float pxGen2 = analysisTree.genparticles_px[index2];
	  float pyGen2 = analysisTree.genparticles_py[index2];
	  float pzGen2 = analysisTree.genparticles_pz[index2];
	  float ptGen2 = PtoPt(pxGen2,pyGen2);
	  float etaGen2 = PtoEta(pxGen2,pyGen2,pzGen2);
	  float phiGen2 = TMath::ATan2(pyGen2,pxGen2);

	  float DR = deltaR(etaGen1,phiGen1,
			    etaGen2,phiGen2);
	  float DPhi = dPhiFrom2P(pxGen1, pyGen1,
				  pxGen2, pyGen2);

	  if (ptGen1>ptMuonTagCut && fabs(etaGen1)<etaMuonTagCut) {
	    if (q1<0) posMuonInAcceptance = true;
	    if (q1>0) negMuonInAcceptance = true;
	  }
	  if (ptGen2>ptMuonTagCut && fabs(etaGen2)<etaMuonTagCut) {
	    if (q2<0) posMuonInAcceptance = true;
	    if (q2>0) negMuonInAcceptance = true;
	  }

	  // for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	  //   TLorentzVector muonReco4P; muonReco4P.SetXYZM(analysisTree.muon_px[im],
	  // 						  analysisTree.muon_py[im],
	  // 						  analysisTree.muon_pz[im],
	  // 						  muonMass);
	  //   float dP = (muonGen4P-muonReco4P).P();
	  //   float dPoverP = dP/muonGen4P.P();
	  //   if (dPoverP<0.1) { 
	  // 	fabs(analysisTree.muon_eta[im])<etaMuonCut &&
	  // 	analysisTree.muon_pt[im]>ptMuonCut &&
	  // 	analysisTree.muon_pt[im]<99.99) {
	  //     dPMin = dPoverP;
	  //     matchedMuon = im;
	  //   }
	  // }

	  if (DR>dRLeptonsCut && DPhi>dPhiLeptonsLowCut && DPhi<dPhiLeptonsHighCut)
	    isZToMMopen = true;

	}
      }

      if (indexE.size()==2) {
	unsigned int index1 = indexE.at(0);
	unsigned int index2 = indexE.at(1);
	int q1 = analysisTree.genparticles_pdgid[index1];
	int q2 = analysisTree.genparticles_pdgid[index2];
	if ( (q1*q2<0) && 
	     (analysisTree.genparticles_info[index1]==1) &&
	     (analysisTree.genparticles_info[index2]==1) ) { 
	  isZToEE = true;
	  
	  float pxGen1 = analysisTree.genparticles_px[index1];
	  float pyGen1 = analysisTree.genparticles_py[index1];
	  float pzGen1 = analysisTree.genparticles_pz[index1];
	  float ptGen1 = PtoPt(pxGen1,pyGen1);
	  float etaGen1 = PtoEta(pxGen1,pyGen1,pzGen1);
	  float phiGen1 = TMath::ATan2(pyGen1,pxGen1);
	  
	  float pxGen2 = analysisTree.genparticles_px[index2];
	  float pyGen2 = analysisTree.genparticles_py[index2];
	  float pzGen2 = analysisTree.genparticles_pz[index2];
	  float ptGen2 = PtoPt(pxGen2,pyGen2);
	  float etaGen2 = PtoEta(pxGen2,pyGen2,pzGen2);
	  float phiGen2 = TMath::ATan2(pyGen2,pxGen2);

	  float DR = deltaR(etaGen1,phiGen1,
			    etaGen2,phiGen2);
	  float DPhi = dPhiFrom2P(pxGen1, pyGen1,
				  pxGen2, pyGen2);

	  if (ptGen1>ptMuonTagCut && fabs(etaGen1)<etaMuonTagCut) {
	    if (q1<0) posElectronInAcceptance = true;
	    if (q1>0) negElectronInAcceptance = true;
	  }
	  if (ptGen2>ptMuonTagCut && fabs(etaGen2)<etaMuonTagCut) {
	    if (q2<0) posElectronInAcceptance = true;
	    if (q2>0) negElectronInAcceptance = true;
	  }

	  if (DR>dRLeptonsCut && DPhi>dPhiLeptonsLowCut && DPhi<dPhiLeptonsHighCut)
	    isZToEEopen = true;

	}
      }
      bool isZToTauTau = (!isZToEE) && (!isZToMM) && (!isTTJets);

      bool posElectronFromTauInAcceptance = false;
      bool negElectronFromTauInAcceptance = false;

      bool posMuonFromTauInAcceptance = false;
      bool negMuonFromTauInAcceptance = false;

      bool posTauInAcceptance = false;
      bool negTauInAcceptance = false;
      
      if (isZToTauTau) {
	if (indexMu.size()==1 && indexE.size()==1) {

	  for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++ itau) {
	    float pxGen = analysisTree.gentau_px[itau];
	    float pyGen = analysisTree.gentau_py[itau];
	    float pzGen = analysisTree.gentau_pz[itau];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float ptGen  = PtoPt(pxGen,pyGen);
	    indexTau.push_back(itau);
	    float dRMax = 1;
	    int chargeTau = 0;
	    bool isTauToMu = true;

	    for (unsigned int iM=0; iM<indexMu.size(); ++iM) { // over muons
	      unsigned int index = indexMu[iM];
	      float pxGenMu = analysisTree.genparticles_px[index];
	      float pyGenMu = analysisTree.genparticles_py[index];
	      float pzGenMu = analysisTree.genparticles_pz[index];
	      float etaGenMu = PtoEta(pxGen,pyGenMu,pzGenMu);
	      float phiGenMu = TMath::ATan2(pyGenMu,pxGenMu);
	      float dRTauMu = deltaR(etaGen,phiGen,
				     etaGenMu,phiGenMu);
	      //	      std::cout << "dRTauMu = " << dRTauMu << std::endl;
	      if (dRTauMu<dRMax) {
		dRMax = dRTauMu;
		if (analysisTree.genparticles_pdgid[index]>0)
		  chargeTau = -1;
		else 
		  chargeTau = 1;
	      }
	      
	    }
	    
	    for (unsigned int iM=0; iM<indexE.size(); ++iM) { // over electrons
	      unsigned int index = indexE[iM];
	      float pxGenMu = analysisTree.genparticles_px[index];
	      float pyGenMu = analysisTree.genparticles_py[index];
	      float pzGenMu = analysisTree.genparticles_pz[index];
	      float etaGenMu = PtoEta(pxGen,pyGenMu,pzGenMu);
	      float phiGenMu = TMath::ATan2(pyGenMu,pxGenMu);
	      float dRTauMu = deltaR(etaGen,phiGen,
				     etaGenMu,phiGenMu);
	      //	      std::cout << "dRTauE = " << dRTauMu << std::endl;
	      if (dRTauMu<dRMax) {
		dRMax = dRTauMu;
		isTauToMu = false;
		if (analysisTree.genparticles_pdgid[index]>0)
		  chargeTau = -1;
		else 
		  chargeTau = 1;
	      }
	      
	    }
	    float minTauPt = ptElectronCut;
	    float maxTauEta = etaElectronCut;
	    if (isTauToMu) {
	      minTauPt = ptMuonCut;
	      maxTauEta = etaMuonCut;
	    }

	    qTau.push_back(chargeTau);
	    if (ptGen>minTauPt && fabs(etaGen)<maxTauEta) {
	      if (chargeTau>0)
		posTauInAcceptance = true;
	      else 
		negTauInAcceptance = true;
	    }

	    //	    std::cout << "tau    pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen
	    //		      << "   charge = " << chargeTau  << std::endl; 
	  }
	  //	  std::cout << "++++++++++++++++++++++++++++++" << std::endl;
	  //	  std::cout << std::endl;
	}

	if (indexE.size()==1) {
	  unsigned int index = indexE.at(0);
	  int info = analysisTree.genparticles_info[index];
	  if (info==5) {
	    int q = analysisTree.genparticles_pdgid[index];
	    float pxGen = analysisTree.genparticles_px[index];
	    float pyGen = analysisTree.genparticles_py[index];
	    float pzGen = analysisTree.genparticles_pz[index];
	    float ptGen = PtoPt(pxGen,pyGen);
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    if (ptGen>ptElectronTagCut && fabs(etaGen)<etaElectronTagCut) {
	      if (q<0)
		posElectronFromTauInAcceptance = true;
	      else
		negElectronFromTauInAcceptance = true;
	    }
	  }	  
	}
	if (indexMu.size()==1) {
	  unsigned int index = indexMu.at(0);
	  int info = analysisTree.genparticles_info[index];
	  if (info==5) {
	    int q = analysisTree.genparticles_pdgid[index];
	    float pxGen = analysisTree.genparticles_px[index];
	    float pyGen = analysisTree.genparticles_py[index];
	    float pzGen = analysisTree.genparticles_pz[index];
	    float ptGen = PtoPt(pxGen,pyGen);
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    if (ptGen>ptMuonTagCut && fabs(etaGen)<etaMuonTagCut) {
	      if (q<0)
		posMuonFromTauInAcceptance = true;
	      else
		negMuonFromTauInAcceptance = true;
	    }
	  }	  
	}	

      }

      // negMuonInAcceptance = true;
      // posMuonInAcceptance = true;
      
      // negElectronInAcceptance = true;
      // posElectronInAcceptance = true;

      // ****************
      // trigger matching
      // ****************

      bool isIsoMu24fired = false;
      bool isMu23fired    = false;
      bool isMu8fired     = false;
      bool isEle27fired   = false;
      bool isEle23fired   = false;
      bool isEle12fired   = false;

      
      bool isAny = false;
      if (isAny) {
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  isAny = analysisTree.trigobject_filters[iT][6]||analysisTree.trigobject_filters[iT][7];
	  if (isAny) break;
	}
	if (isAny) {
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    
	    bool trigFilter = false;
	    TString filterName;
	    float etaTrig = analysisTree.trigobject_eta[iT];
	    float phiTrig = analysisTree.trigobject_phi[iT];
	    float ptTrig  = analysisTree.trigobject_pt[iT];
	    if (analysisTree.trigobject_filters[iT][2]||analysisTree.trigobject_filters[iT][6]||
		analysisTree.trigobject_filters[iT][3]||analysisTree.trigobject_filters[iT][7]) {
	      printf("%2i  %5.1f  %5.1f  %6.1f  : ",iT,
		     analysisTree.trigobject_pt[iT],
		     analysisTree.trigobject_eta[iT],
		     analysisTree.trigobject_phi[iT]);
	      cout << " Mu23=" << analysisTree.trigobject_filters[iT][2]
		   << " Mu8=" << analysisTree.trigobject_filters[iT][3]
		   << " El23=" << analysisTree.trigobject_filters[iT][7]
		   << " El12=" << analysisTree.trigobject_filters[iT][6] << endl;
	      
	    }
	  }
	}
	if (isAny) { 
	  cout << "Mu23Ele12 : " << analysisTree.hltriggerresults_second[5] << endl;
	  cout << "Mu8Ele23  : " << analysisTree.hltriggerresults_second[6] << endl;
	  cout << endl;
	}      
      }

      //      continue;
      for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	bool trigFilter = false;
	TString filterName;
	float etaTrig = analysisTree.trigobject_eta[iT];
	float phiTrig = analysisTree.trigobject_phi[iT];
	float ptTrig  = analysisTree.trigobject_pt[iT];

	if (analysisTree.trigobject_filters[iT][0]) { // IsoMu24 leg
	  trigFilter = true;
	  filterName = "IsoMu24";
	  for (unsigned int J=0; J<indexMu.size(); ++J) {
	    int indexJ = indexMu.at(J);
	    //	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isIsoMu24fired = true;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][2]) { // Mu23 Leg
	  trigFilter = true;
	  filterName = "Mu23";
	  for (unsigned int J=0; J<indexMu.size(); ++J) {
	    int indexJ = indexMu.at(J);
	    //	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isMu23fired = true;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][3]) { // Mu8 Leg
	  trigFilter = true;
	  filterName = "Mu8";
	  for (unsigned int J=0; J<indexMu.size(); ++J) {
	    int indexJ = indexMu.at(J);
	    //	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isMu8fired = true;
	    }
	  }
	}
	  
	if (analysisTree.trigobject_filters[iT][4]) { // Ele27 leg
	  trigFilter = true;
	  filterName = "Ele27";
	  for (unsigned int J=0; J<indexE.size(); ++J) {
	    int indexJ = indexE.at(J);
	    //	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isEle27fired = true;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][7]) { // Ele23 leg
	  trigFilter = true;
	  filterName = "Ele23";
	  for (unsigned int J=0; J<indexE.size(); ++J) {
	    int indexJ = indexE.at(J);
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isEle23fired = true;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][6]) { // Ele12 leg
	  trigFilter = true;
	  filterName = "Ele12";
	  for (unsigned int J=0; J<indexE.size(); ++J) {
	    int indexJ = indexE.at(J);
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isEle12fired = true;
	    }
	  }

	}

	if (trigFilter) {
	//   std::cout << filterName 
	// 	    << "  pt = " << ptTrig
	// 	    << "  eta = " << etaTrig
	// 	    << "  phi = " << phiTrig << std::endl;
	}

      }

      float ele23fired = 0;
      float ele12fired = 0;
      float mu23fired = 0;
      float mu8fired = 0;
      

      if (isEle23fired) 
	ele23fired = 1;
      if (isEle12fired)
	ele12fired = 1;
      if (isMu23fired)
	mu23fired = 1;
      if (isMu8fired)
	mu8fired = 1;

      
      if (isZToTauTau) {
	ZToTauTauH->Fill(float(0),weight);
	if (indexMu.size()!=1&&indexE.size()>0) {
	  for (unsigned int i=0; i<indexMu.size(); ++i) {
	    unsigned int index = indexMu.at(i);
	  }
	}
	if (indexE.size()!=1&&indexMu.size()>0) {
	  for (unsigned int i=0; i<indexE.size(); ++i) {
	    unsigned int index = indexE.at(i);
	  }
	}
      }

      if (indexMu.size()==1&&indexE.size()==1) {
	ZToTauTauH->Fill(float(1),weight);
	if (isZToTauTau) ZToTauTauH->Fill(float(2),weight);
	int indexJ = indexMu.at(0);
	float pxGen = analysisTree.genparticles_px[indexJ];
	float pyGen = analysisTree.genparticles_py[indexJ];
	float pzGen = analysisTree.genparticles_pz[indexJ];
	float etaGen = PtoEta(pxGen,pyGen,pzGen);
	float ptGen =  PtoPt(pxGen,pyGen);
	if (fabs(etaGen)<etaMuonCut) {
	  if (isEle27fired) {
	    MuonPt23H->Fill(ptGen,weight);
	    MuonPt8H->Fill(ptGen,weight);
	      if (isMu23fired)
		MuonPt23PassedH->Fill(ptGen,weight);
	    if (isMu8fired)
	      MuonPt8PassedH->Fill(ptGen,weight);
	  }
	}

	indexJ = indexE.at(0);
        pxGen = analysisTree.genparticles_px[indexJ];
        pyGen = analysisTree.genparticles_py[indexJ];
        pzGen = analysisTree.genparticles_pz[indexJ];
        etaGen = PtoEta(pxGen,pyGen,pzGen);
        ptGen =  PtoPt(pxGen,pyGen);
	if (fabs(etaGen)<etaElectronCut) {
	  if (isMu23fired) {
	    ElectronPt12H->Fill(ptGen,weight);
	    if (isEle12fired)
              ElectronPt12PassedH->Fill(ptGen,weight);
	  }
	  if (isMu8fired) {
	    ElectronPt23H->Fill(ptGen,weight);
	    if (isEle23fired)
	      ElectronPt23PassedH->Fill(ptGen,weight);
	  }
	}

      }
      
      if (posMuonFromTauInAcceptance && negElectronFromTauInAcceptance) {      
	Mu23Ele12H->Fill(mu23fired,ele12fired,weight);
	Mu8Ele23H->Fill(mu8fired,ele23fired,weight);
      }
       if (negMuonFromTauInAcceptance && posElectronFromTauInAcceptance) {      
	Mu23Ele12H->Fill(mu23fired,ele12fired,weight);
	Mu8Ele23H->Fill(mu8fired,ele23fired,weight);
      }

      
      // std::cout << "Z Decay : " << isZToEE << isZToMM << isZToTauTau << std::endl; 
      // std::cout << std::endl;

      // if (isTTJets) {
      // 	isZToEE = true;
      // 	isZToMM = true;
      // 	isZToEEopen = true;
      // 	isZToMMopen = true;
      // }

      // *********************************************
      // selecting tag and probe electron (highest pt)
      // *********************************************

      int indexTagPosE = -1;
      int indexTagNegE = -1;

      float ptMinTagPosE = 0;
      float ptMinTagNegE = 0;

      int indexPassPosE = -1;
      int indexPassNegE = -1;

      float ptMinPassPosE = 0;
      float ptMinPassNegE = 0;

      int indexProbePosE = -1;
      int indexProbeNegE = -1;

      float ptMinProbePosE = 0;
      float ptMinProbeNegE = 0;
     
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]>ptElectronCut&&
	    analysisTree.electron_pt[ie]<99.99&&
	    fabs(analysisTree.electron_eta[ie])<etaElectronCut
	    //	    &&
	    //	    fabs(analysisTree.electron_dz[ie])<0.5&&
	    //	    fabs(analysisTree.electron_dxy[ie])<0.2
	    ) {
	  if (analysisTree.electron_pt[ie]>ptMinProbePosE&&analysisTree.electron_charge[ie]>0.5) {
	    indexProbePosE = ie;
	    ptMinProbePosE = analysisTree.electron_pt[ie];
	  }
	  if (analysisTree.electron_pt[ie]>ptMinProbeNegE&&analysisTree.electron_charge[ie]<-0.5) {
	    indexProbeNegE = ie;
	    ptMinProbeNegE = analysisTree.electron_pt[ie];
	  }
	}
	if (applyElectronTriggerMatchTag && analysisTree.hltriggerresults_second[0]==0) continue; // Ele27
	if (analysisTree.electron_pt[ie]<ptElectronTagCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronTagCut) continue;
	if (applyElectronIPcutsTag && fabs(analysisTree.electron_dxy[ie])>dxyElectronTagCut) continue;
	if (applyElectronIPcutsTag && fabs(analysisTree.electron_dz[ie])>dzElectronTagCut) continue;
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ie];
	if (applyElectronIsoTag && relIso>isoElectronHighTagCut) continue;
	if (applyElectronIsoTag && relIso<isoElectronLowTagCut) continue;
	bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[ie],
						analysisTree.electron_mva_id_nontrigPhys14[ie]);
	if (applyElectronIdTag && !electronMvaId) continue;
	if (applyElectronIdTag && !analysisTree.electron_pass_conversion[ie]) continue;
	if (applyElectronIdTag && analysisTree.electron_nmissinginnerhits[ie]!=0) continue;
	if (analysisTree.electron_pt[ie]>ptMinPassPosE&&analysisTree.electron_charge[ie]>0.5) {
	  indexPassPosE = ie;
	  ptMinPassPosE = analysisTree.electron_pt[ie];
	}
	if (analysisTree.electron_pt[ie]>ptMinPassNegE&&analysisTree.electron_charge[ie]<-0.5) {
	  indexPassNegE = ie;
	  ptMinPassNegE = analysisTree.electron_pt[ie];
	}
	bool trigMatch = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][4]) { // Single Ele27 Leg
	    float dRtrig = deltaR(analysisTree.electron_eta[ie],analysisTree.electron_phi[ie],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<DRTrigMatch) {
	      trigMatch = true;
	      break;
	    }
	  }
	}
	if (applyElectronTriggerMatchTag && !trigMatch) continue;
      	if (analysisTree.electron_pt[ie]>ptMinTagPosE&&analysisTree.electron_charge[ie]>0.5) {
      	  indexTagPosE = ie;
      	  ptMinTagPosE = analysisTree.electron_pt[ie];
      	}
      	if (analysisTree.electron_pt[ie]>ptMinTagNegE&&analysisTree.electron_charge[ie]<-0.5) {
      	  indexTagNegE = ie;
      	  ptMinTagNegE = analysisTree.electron_pt[ie];
      	}
      }

      if ( !useGeneratorAsTag ) {
	negElectronInAcceptance = indexTagNegE>=0;
	posElectronInAcceptance = indexTagPosE>=0;
      }

      // ******************************************
      // selecting tag and probe muons (highest pt)
      // ******************************************

      int indexTagPosMu = -1;
      int indexTagNegMu = -1;

      float ptMinTagPosMu = 0;
      float ptMinTagNegMu = 0;

      int indexPassPosMu = -1;
      int indexPassNegMu = -1;

      float ptMinPassPosMu = 0;
      float ptMinPassNegMu = 0;

      int indexProbePosMu = -1;
      int indexProbeNegMu = -1;

      float ptMinProbePosMu = 0;
      float ptMinProbeNegMu = 0;

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]>ptMuonCut&&
	    analysisTree.muon_pt[im]<99.99&&
	    fabs(analysisTree.muon_eta[im])<etaMuonCut
	    //	    &&
	    //	    fabs(analysisTree.muon_dz[im])<0.5&&
	    //	    fabs(analysisTree.muon_dxy[im])<0.2
	    ) {
	  if (analysisTree.muon_pt[im]>ptMinProbePosMu&&analysisTree.muon_charge[im]>0.5 ) {
	    indexProbePosMu = im;
	    ptMinProbePosMu = analysisTree.muon_pt[im];
	  }
	  if (analysisTree.muon_pt[im]>ptMinProbeNegMu&&analysisTree.muon_charge[im]<-0.5) {
	    indexProbeNegMu = im;
	    ptMinProbeNegMu = analysisTree.muon_pt[im];
	  }
	}
	if (applyMuonTriggerMatchTag && analysisTree.hltriggerresults_second[2]==0) continue; // IsoMu24
	if (analysisTree.muon_pt[im]<ptMuonTagCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonTagCut) continue;
	if (applyMuonIPcutsTag && fabs(analysisTree.muon_dxy[im])>dxyMuonTagCut) continue;
	if (applyMuonIPcutsTag && fabs(analysisTree.muon_dz[im])>dzMuonTagCut) continue;
	float neutralIso = 
	  analysisTree.muon_neutralHadIso[im] + 
	  analysisTree.muon_photonIso[im] - 
	  0.5*analysisTree.muon_puIso[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.muon_chargedHadIso[im] + neutralIso;
	float relIso = absIso/analysisTree.muon_pt[im];
	if (applyMuonIsoTag && relIso>isoMuonHighTagCut) continue;
	if (applyMuonIsoTag && relIso<isoMuonLowTagCut) continue;
	if (applyMuonIdTag==1 && !analysisTree.muon_isLoose[im]) continue;
	if (applyMuonIdTag==2 && !analysisTree.muon_isMedium[im]) continue;
	if (applyMuonIdTag==3 && !analysisTree.muon_isTight[im]) continue;
	if (analysisTree.muon_pt[im]>ptMinPassPosMu&&analysisTree.muon_charge[im]>0.5) {
          indexPassPosMu = im;
          ptMinPassPosMu = analysisTree.muon_pt[im];
        }
        if (analysisTree.muon_pt[im]>ptMinPassNegMu&&analysisTree.muon_charge[im]<-0.5) {
          indexPassNegMu = im;
          ptMinPassNegMu = analysisTree.muon_pt[im];
        }
	bool trigMatch = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][0]) { // Single IsoMu24 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<DRTrigMatch) {
	      trigMatch = true;
	      break;
	    }
	  }
	}
	if (applyMuonTriggerMatchTag && !trigMatch) continue;
      	if (analysisTree.muon_pt[im]>ptMinTagPosMu&&analysisTree.muon_charge[im]>0.5) {
      	  indexTagPosMu = im;
      	  ptMinTagPosMu = analysisTree.muon_pt[im];
      	}
      	if (analysisTree.muon_pt[im]>ptMinTagNegMu&&analysisTree.muon_charge[im]<-0.5) {
      	  indexTagNegMu = im;
      	  ptMinTagNegMu = analysisTree.muon_pt[im];
      	}
      }

      if ( !useGeneratorAsTag ) {
	negMuonInAcceptance = indexTagNegMu>=0;
	posMuonInAcceptance = indexTagPosMu>=0;
      }
      
      // MEt computation
      float MEt = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex+analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // *****************************
      // Loop over generator electrons 
      // *****************************
      for (unsigned int iE=0; iE<indexE.size(); ++iE) {
	unsigned int index = indexE.at(iE);
	
	int electronQ = 0; 
	if (analysisTree.genparticles_pdgid[index]>0) 
	  electronQ = 1;
	TLorentzVector electronGen4P; electronGen4P.SetXYZM(analysisTree.genparticles_px[index],
							    analysisTree.genparticles_py[index],
							    analysisTree.genparticles_pz[index],
							    electronMass);
	bool isPrompt = true;
	if (analysisTree.genparticles_info[index]==5||
	    analysisTree.genparticles_info[index]==6) 
	  isPrompt = false;

	float electronGenEta = electronGen4P.Eta();
	float electronGenPt  = electronGen4P.Pt();
	float electronGenPhi = electronGen4P.Phi();

	bool isEventSelected = isZToEE && ((electronQ==0&&negElectronInAcceptance) || (electronQ==1&&posElectronInAcceptance)) && isZToEEopen;

	if (isZToTauTau || isTTJets) {
	  isEventSelected = false;
	  if (indexMu.size()==1) {
	    int indexMuon = indexMu[0];
	    TLorentzVector muonGen4P; muonGen4P.SetXYZM(analysisTree.genparticles_px[indexMuon],
							analysisTree.genparticles_py[indexMuon],
							analysisTree.genparticles_pz[indexMuon],
							muonMass);
	    float muonGenPt  = muonGen4P.Pt();
	    float muonGenEta = muonGen4P.Eta();
	    float muonGenPhi = muonGen4P.Phi();
	    
	    float deltaRMuEle = deltaR(muonGenEta,muonGenPhi,
				       electronGenEta,electronGenPhi);

	    if (fabs(muonGenEta)<etaMuonCut && deltaRMuEle>dRLeptonsCut) {
	      isEventSelected = (muonGenPt>ptMuonCut && electronGenPt>ptElectronHighCut) ||
		(muonGenPt>ptMuonHighCut && electronGenPt>ptElectronCut);
	    }

	  }

	}	  

	int etaBinGen = binNumber(fabs(electronGenEta),nEtaBins,etaBins);

	float dPMin = 0.05;
	int matchedElectron = -1;
	if (isPrompt) {
	  ElectronPtH[electronQ][0][0][etaBinGen]->Fill(electronGenPt,weight);
	}
	else {
	  ElectronPtH[electronQ][1][0][etaBinGen]->Fill(electronGenPt,weight);
	}

	for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	  TLorentzVector electronReco4P; electronReco4P.SetXYZM(analysisTree.electron_px[ie],
								analysisTree.electron_py[ie],
								analysisTree.electron_pz[ie],
								electronMass);
	  float dP = (electronGen4P-electronReco4P).P();
	  float dPoverP = dP/electronGen4P.P();
	  if (dPoverP<dPMin && 
	      fabs(analysisTree.electron_eta[ie])<etaElectronCut &&
	      analysisTree.electron_pt[ie]>ptElectronCut &&
	      analysisTree.electron_pt[ie]<99.99 
	      //	      && 
	      //	      fabs(analysisTree.electron_dz[ie])<0.5 && 
	      //	      fabs(analysisTree.electron_dxy[ie])<0.2
	      ) {
	    dPMin = dPoverP;
	    matchedElectron = ie;
	  }
	}
	if (matchedElectron>=0) {

	  //	  std::cout << " Electron -> IsPrompt = " << isPrompt << "  IsTTJets = " << isTTJets << std::endl;

	  float electronRecoPt = analysisTree.electron_pt[matchedElectron];
	  float electronRecoEta = analysisTree.electron_eta[matchedElectron];

	  int etaBinReco = binNumber(fabs(electronRecoEta),nEtaBins,etaBins);
	  //	  std::cout << "electron eta = " << electronRecoEta << " : Bin = " << etaBinReco << std::endl;

	  int iPt = 0;
	  if (electronRecoPt>20) 
	    iPt = 1;

	  float dRLepParton = 4.99;
	  if (nPartonLeading>=0) {
	    TLorentzVector parton4P; parton4P.SetXYZM(analysisTree.genparticles_px[nPartonLeading],
						      analysisTree.genparticles_py[nPartonLeading],
						      analysisTree.genparticles_pz[nPartonLeading],0.);
	    float partonGenEta = parton4P.Eta();
	    float partonGenPhi = parton4P.Phi();
	    dRLepParton = deltaR(electronGenEta,electronGenPhi,
				 partonGenEta,partonGenPhi);
	    
	  }
	  
	  if (isPrompt) { 

	    ElectronPtH[electronQ][0][1][etaBinGen]->Fill(electronGenPt,weight);

	    ElectronPtH[electronQ][0][3][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][0][3][iPt]->Fill(electronRecoEta,weight);
	      
	    ElectronVertexDxH[electronQ][0][3][etaBinReco]->Fill(dVertX,weight);
	    ElectronVertexDyH[electronQ][0][3][etaBinReco]->Fill(dVertY,weight);
	    ElectronVertexDzH[electronQ][0][3][etaBinReco]->Fill(dVertZ,weight);
	    
	    ElectronNvertH[electronQ][0][3]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertEtaH[electronQ][0][3][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertPtH[electronQ][0][3][iPt]->Fill(float(analysisTree.primvertex_count),weight);

	    dRElectronLeadingPartonH[electronQ][0][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    if (isEventSelected) { 

	      ElectronPtH[electronQ][0][11][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][0][11][iPt]->Fill(electronRecoEta,weight);
	      
	      ElectronVertexDxH[electronQ][0][11][etaBinReco]->Fill(dVertX,weight);
	      ElectronVertexDyH[electronQ][0][11][etaBinReco]->Fill(dVertY,weight);
	      ElectronVertexDzH[electronQ][0][11][etaBinReco]->Fill(dVertZ,weight);
	      
	      ElectronNvertH[electronQ][0][11]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertEtaH[electronQ][0][11][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertPtH[electronQ][0][11][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	      dRElectronLeadingPartonTagH[electronQ][0][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    }
	    

	  }
	  else {

	    ElectronPtH[electronQ][1][1][etaBinGen]->Fill(electronGenPt,weight);

	    ElectronPtH[electronQ][1][3][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][1][3][iPt]->Fill(electronRecoEta,weight);

	    dRElectronLeadingPartonH[electronQ][1][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    if (isEventSelected) { 

	      ElectronPtH[electronQ][1][11][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][1][11][iPt]->Fill(electronRecoEta,weight);

	      dRElectronLeadingPartonTagH[electronQ][1][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    }

	  }

	  if (fabs(analysisTree.electron_dxy[matchedElectron])>dxyElectronCut&&applyElectronIPcuts) continue;
	  if (fabs(analysisTree.electron_dz[matchedElectron])>dzElectronCut&&applyElectronIPcuts) continue;
	  bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[matchedElectron],
						  analysisTree.electron_mva_id_nontrigPhys14[matchedElectron]);
	  if (!electronMvaId&&applyElectronId) continue;
	  if (!analysisTree.electron_pass_conversion[matchedElectron]&&applyElectronId) continue;
	  if (analysisTree.electron_nmissinginnerhits[matchedElectron]!=0&&applyElectronId) continue;

	  float neutralIso = 
	    analysisTree.electron_neutralHadIso[matchedElectron] + 
	    analysisTree.electron_photonIso[matchedElectron];
	  float neutralIsoCorr = neutralIso - 0.5*analysisTree.electron_puIso[matchedElectron];
	  neutralIsoCorr = TMath::Max(float(0),neutralIsoCorr); 
	  float absIso = analysisTree.electron_chargedHadIso[matchedElectron] + neutralIsoCorr;
	  float relIso = absIso/analysisTree.electron_pt[matchedElectron];

	  if (isPrompt) {

	    ElectronChargedIsoH[electronQ][0][13][etaBinReco]->Fill(analysisTree.electron_chargedHadIso[matchedElectron],weight);
	    ElectronNeutralIsoH[electronQ][0][13][etaBinReco]->Fill(neutralIso,weight);
	    ElectronNeutralCorrIsoH[electronQ][0][13][etaBinReco]->Fill(neutralIsoCorr,weight);
	    ElectronIsoH[electronQ][0][13][etaBinReco]->Fill(absIso,weight);
	    ElectronRelIsoH[electronQ][0][13][etaBinReco]->Fill(relIso,weight);

	    ElectronPtH[electronQ][0][13][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][0][13][iPt]->Fill(electronRecoEta,weight);
	    
	    ElectronVertexDxH[electronQ][0][13][etaBinReco]->Fill(dVertX,weight);
	    ElectronVertexDyH[electronQ][0][13][etaBinReco]->Fill(dVertY,weight);
	    ElectronVertexDzH[electronQ][0][13][etaBinReco]->Fill(dVertZ,weight);
	    
	    ElectronNvertH[electronQ][0][13]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertEtaH[electronQ][0][13][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertPtH[electronQ][0][13][iPt]->Fill(float(analysisTree.primvertex_count),weight);

	    if (isEventSelected) {

	      ElectronChargedIsoH[electronQ][0][14][etaBinReco]->Fill(analysisTree.electron_chargedHadIso[matchedElectron],weight);
	      ElectronNeutralIsoH[electronQ][0][14][etaBinReco]->Fill(neutralIso,weight);
	      ElectronNeutralCorrIsoH[electronQ][0][14][etaBinReco]->Fill(neutralIsoCorr,weight);
	      ElectronIsoH[electronQ][0][14][etaBinReco]->Fill(absIso,weight);
	      ElectronRelIsoH[electronQ][0][14][etaBinReco]->Fill(relIso,weight);
	      
	      ElectronPtH[electronQ][0][14][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][0][14][iPt]->Fill(electronRecoEta,weight);
	      
	      ElectronVertexDxH[electronQ][0][14][etaBinReco]->Fill(dVertX,weight);
	      ElectronVertexDyH[electronQ][0][14][etaBinReco]->Fill(dVertY,weight);
	      ElectronVertexDzH[electronQ][0][14][etaBinReco]->Fill(dVertZ,weight);
	      
	      ElectronNvertH[electronQ][0][14]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertEtaH[electronQ][0][14][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertPtH[electronQ][0][14][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	    }

	  }
	  else {
	    ElectronPtH[electronQ][1][13][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][1][13][iPt]->Fill(electronRecoEta,weight);

	    if (isEventSelected) { 
	      ElectronPtH[electronQ][1][14][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][1][14][iPt]->Fill(electronRecoEta,weight);
	    }

	  }

	  if (relIso>isoElectronHighCut&&applyElectronIso) continue;
	  if (relIso<isoElectronLowCut&&applyElectronIso) continue;

	  if (isPrompt) {

	    ElectronPtH[electronQ][0][2][etaBinGen]->Fill(electronGenPt,weight);

	    ElectronPtH[electronQ][0][4][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][0][4][iPt]->Fill(electronRecoEta,weight);
	    
	    ElectronNvertH[electronQ][0][4]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertEtaH[electronQ][0][4][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    ElectronNvertPtH[electronQ][0][4][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	    
	    if (singleMuTrigger) {
	      if (isIsoMu24fired) {
		ElectronPtH[electronQ][0][5][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPtH[electronQ][0][7][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt23RecoXH->Fill(electronRecoPt,weight);
		ElectronPt12RecoXH->Fill(electronRecoPt,weight);
		if (isEle23fired) {
		  ElectronPtH[electronQ][0][6][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt23RecoPassedXH->Fill(electronRecoPt,weight);
		}
		if (isEle12fired) {
		  ElectronPtH[electronQ][0][8][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt12RecoPassedXH->Fill(electronRecoPt,weight);
		}
	      }
	      
	    }
	    else {
	      if (isMu23fired) {
		ElectronPtH[electronQ][0][7][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt12RecoXH->Fill(electronRecoPt,weight);
		if (isEle12fired) {
		  ElectronPtH[electronQ][0][8][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt12RecoPassedXH->Fill(electronRecoPt,weight);
		}
	      }
	      if (isMu8fired) {
		ElectronPtH[electronQ][0][5][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt23RecoXH->Fill(electronRecoPt,weight);
		if (isEle23fired) {
		  ElectronPtH[electronQ][0][6][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt23RecoPassedXH->Fill(electronRecoPt,weight);
		}
	      }
	    }

 	    if (isEventSelected) { 
	      
	      ElectronPtH[electronQ][0][12][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][0][12][iPt]->Fill(electronRecoEta,weight);
	      
	      ElectronNvertH[electronQ][0][12]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertEtaH[electronQ][0][12][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      ElectronNvertPtH[electronQ][0][12][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	      if (MEt>80) {
		if (singleMuTrigger) {
		  if (isIsoMu24fired) {
		    ElectronPt23RecoYH->Fill(electronRecoPt,weight);
		    ElectronPt12RecoYH->Fill(electronRecoPt,weight);
		    if (isEle23fired) {
		      ElectronPt23RecoPassedYH->Fill(electronRecoPt,weight);
		    }
		    if (isEle12fired) {
		      ElectronPt12RecoPassedYH->Fill(electronRecoPt,weight);
		    }
		  }
		}
		else {
		  if (isMu23fired) {
		    ElectronPt12RecoYH->Fill(electronRecoPt,weight);
		    if (isEle12fired)
		      ElectronPt12RecoPassedYH->Fill(electronRecoPt,weight);
		  }
		  if (isMu8fired) {
		    ElectronPt23RecoYH->Fill(electronRecoPt,weight);
		    if (isEle23fired)
		      ElectronPt23RecoPassedYH->Fill(electronRecoPt,weight);
		  }
		}
	      }
	    }
	  }
	  else {
	    
	    ElectronPtH[electronQ][1][2][etaBinGen]->Fill(electronGenPt,weight);
	    
	    ElectronPtH[electronQ][1][4][etaBinReco]->Fill(electronRecoPt,weight);
	    ElectronEtaH[electronQ][1][4][iPt]->Fill(electronRecoEta,weight);
	    
	    if (isEventSelected) { 
	      ElectronPtH[electronQ][1][12][etaBinReco]->Fill(electronRecoPt,weight);
	      ElectronEtaH[electronQ][1][12][iPt]->Fill(electronRecoEta,weight);
	    }

	    if (singleMuTrigger) {
	      if (isIsoMu24fired) {
		ElectronPtH[electronQ][1][5][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPtH[electronQ][1][7][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt23RecoH->Fill(electronRecoPt,weight);
		ElectronPt12RecoH->Fill(electronRecoPt,weight);
		if (isEle23fired) {
		  ElectronPtH[electronQ][1][6][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt23RecoPassedH->Fill(electronRecoPt,weight);
		}
		if (isEle12fired) {
		  ElectronPtH[electronQ][1][8][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt12RecoPassedH->Fill(electronRecoPt,weight);
		}
	      }
	    }
	    else {
	      if (isMu23fired) {
		ElectronPtH[electronQ][1][7][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt12RecoH->Fill(electronRecoPt,weight);
		if (isEle12fired) {
		  ElectronPtH[electronQ][1][8][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt12RecoPassedH->Fill(electronRecoPt,weight);
		}
	      }
	      if (isMu8fired) {
		ElectronPtH[electronQ][1][5][etaBinReco]->Fill(electronRecoPt,weight);
		ElectronPt23RecoH->Fill(electronRecoPt,weight);
		if (isEle12fired) {
		  ElectronPtH[electronQ][1][6][etaBinReco]->Fill(electronRecoPt,weight);
		  ElectronPt23RecoPassedH->Fill(electronRecoPt,weight);
		}
	      }
	    }

	  }
	  
	}
	
      }

      // *************************
      // Loop over generator muons
      // *************************
      for (unsigned int iM=0; iM<indexMu.size(); ++iM) {
	unsigned int index = indexMu.at(iM);
	int muonQ = 0; // positive muon
	if (analysisTree.genparticles_pdgid[index]>0) // negative muon 
	  muonQ = 1;
	TLorentzVector muonGen4P; muonGen4P.SetXYZM(analysisTree.genparticles_px[index],
						    analysisTree.genparticles_py[index],
						    analysisTree.genparticles_pz[index],
						    muonMass);
	bool isPrompt = true;
	if (analysisTree.genparticles_info[index]==5||
	    analysisTree.genparticles_info[index]==6)
	  isPrompt = false;

	float muonGenEta = muonGen4P.Eta();
	float muonGenPt  = muonGen4P.Pt();
	float muonGenPhi = muonGen4P.Phi();

	bool isEventSelected = isZToMM && ((muonQ==0&&negMuonInAcceptance) || (muonQ==1&&posMuonInAcceptance)) && isZToMMopen;

	float deltaPhiMuons = 0;
	if (isZToMM) {
	  unsigned int index1 = indexMu.at(0);
	  unsigned int index2 = indexMu.at(1);
	  deltaPhiMuons = dPhiFrom2P(analysisTree.genparticles_px[index1],analysisTree.genparticles_py[index1],
				     analysisTree.genparticles_px[index2],analysisTree.genparticles_py[index2]);
	}
	
	if (isZToTauTau || isTTJets) {
	  isEventSelected = false;
	  if (indexE.size()==1) {
	    int indexEle = indexE[0];
	    TLorentzVector electronGen4P; electronGen4P.SetXYZM(analysisTree.genparticles_px[indexEle],
								analysisTree.genparticles_py[indexEle],
								analysisTree.genparticles_pz[indexEle],
								muonMass);
	    float electronGenPt  = electronGen4P.Pt();
	    float electronGenEta = electronGen4P.Eta();
	    float electronGenPhi = electronGen4P.Phi();

	    float deltaREleMu = deltaR(electronGenEta,electronGenPhi,
				       muonGenEta,muonGenPhi);
	    
	    if (fabs(electronGenEta)<etaElectronCut && deltaREleMu>dRLeptonsCut) {
	      isEventSelected = (muonGenPt>ptMuonCut && electronGenPt>ptElectronHighCut) ||
		(muonGenPt>ptMuonHighCut && electronGenPt>ptElectronCut);

	    }

	  }

	}	  


	int etaBinGen = binNumber(fabs(muonGenEta),nEtaBins,etaBins);

	float dPMin = 0.05;
	int matchedMuon = -1;
	if (isPrompt) {
	  MuonPtH[muonQ][0][0][etaBinGen]->Fill(muonGenPt,weight);
	}
	else {
	  MuonPtH[muonQ][1][0][etaBinGen]->Fill(muonGenPt,weight);
	}

	for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	  TLorentzVector muonReco4P; muonReco4P.SetXYZM(analysisTree.muon_px[im],
							analysisTree.muon_py[im],
							analysisTree.muon_pz[im],
							muonMass);
	  float dP = (muonGen4P-muonReco4P).P();
	  float dPoverP = dP/muonGen4P.P();

	  // std::cout << "dPoverP = " << dPoverP 
	  // 	    << "  id = " << analysisTree.genparticles_pdgid[index]
	  // 	    << "   q = " << analysisTree.muon_charge[im] << std::endl;

	  if (dPoverP<dPMin && 
	      fabs(analysisTree.muon_eta[im])<etaMuonCut &&
	      analysisTree.muon_pt[im]>ptMuonCut &&
	      analysisTree.muon_pt[im]<99.99 
	      //	      && 
	      //	      fabs(analysisTree.muon_dz[im])<0.5 &&
	      //	      fabs(analysisTree.muon_dxy[im])<0.2
	      ) {
	    dPMin = dPoverP;
	    matchedMuon = im;
	  }
	}
	if (matchedMuon>=0) {

	  //	  std::cout << " Muon -> IsPrompt = " << isPrompt << "  IsTTJets = " << isTTJets << std::endl;

	  float muonRecoPt = analysisTree.muon_pt[matchedMuon];
	  float muonRecoEta = analysisTree.muon_eta[matchedMuon];

	  int etaBinReco = binNumber(fabs(muonRecoEta),nEtaBins,etaBins);
	  //	  std::cout << "muon eta = " << muonRecoEta << " : Bin = " << etaBinReco << std::endl;
	  
	  int iPt = 0;
	  if (muonRecoPt>20)
	    iPt = 1;
	  
	  float dRLepParton = 4.99;
	  if (nPartonLeading>=0) {
	    TLorentzVector parton4P; parton4P.SetXYZM(analysisTree.genparticles_px[nPartonLeading],
						      analysisTree.genparticles_py[nPartonLeading],
						      analysisTree.genparticles_pz[nPartonLeading],0.);
	    float partonGenEta = parton4P.Eta();
	    float partonGenPhi = parton4P.Phi();
	    dRLepParton = deltaR(muonGenEta,muonGenPhi,
				 partonGenEta,partonGenPhi);
	  }

	  if (isPrompt) {

	    MuonPtH[muonQ][0][1][etaBinGen]->Fill(muonGenPt,weight);

	    MuonPtH[muonQ][0][3][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][0][3][iPt]->Fill(muonRecoEta,weight);

	    MuonVertexDxH[muonQ][0][3][etaBinReco]->Fill(dVertX,weight);
	    MuonVertexDyH[muonQ][0][3][etaBinReco]->Fill(dVertY,weight);
	    MuonVertexDzH[muonQ][0][3][etaBinReco]->Fill(dVertZ,weight);

	    MuonNvertH[muonQ][0][3]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertEtaH[muonQ][0][3][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertPtH[muonQ][0][3][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	    
	    dRMuonLeadingPartonH[muonQ][0][iPt][etaBinGen]->Fill(dRLepParton,weight);
	    

	    if (isEventSelected) { 
	      
	      MuonPtH[muonQ][0][11][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][0][11][iPt]->Fill(muonRecoEta,weight);
	      
	      MuonVertexDxH[muonQ][0][11][etaBinReco]->Fill(dVertX,weight);
	      MuonVertexDyH[muonQ][0][11][etaBinReco]->Fill(dVertY,weight);
	      MuonVertexDzH[muonQ][0][11][etaBinReco]->Fill(dVertZ,weight);
	      
	      MuonNvertH[muonQ][0][11]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertEtaH[muonQ][0][11][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertPtH[muonQ][0][11][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	      dRMuonLeadingPartonTagH[muonQ][0][iPt][etaBinGen]->Fill(dRLepParton,weight);
	      if (nPartonLeading>=0)
		dRMuonPartonDPhiTagH[muonQ][0][iPt]->Fill(dRLepParton,deltaPhiMuons,weight);

	    }
	  }
	  else {

	    MuonPtH[muonQ][1][1][etaBinGen]->Fill(muonGenPt,weight);

	    MuonPtH[muonQ][1][3][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][1][3][iPt]->Fill(muonRecoEta,weight);

	    dRMuonLeadingPartonH[muonQ][0][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    if (isEventSelected) {

	      MuonPtH[muonQ][1][11][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][1][11][iPt]->Fill(muonRecoEta,weight);

	      dRMuonLeadingPartonTagH[muonQ][1][iPt][etaBinGen]->Fill(dRLepParton,weight);

	    }

	  }

	  if (fabs(analysisTree.muon_dxy[matchedMuon])>dxyMuonCut&&applyMuonIPcuts) continue;
	  if (fabs(analysisTree.muon_dz[matchedMuon])>dzMuonCut&&applyMuonIPcuts) continue;
	  if (applyMuonId==1 && !analysisTree.muon_isLoose[matchedMuon]) continue;
	  if (applyMuonId==2 && !analysisTree.muon_isMedium[matchedMuon]) continue;
	  if (applyMuonId==3 && !analysisTree.muon_isTight[matchedMuon]) continue;

	  float neutralIso = 
	    analysisTree.muon_neutralHadIso[matchedMuon] + 
	    analysisTree.muon_photonIso[matchedMuon]; 
	  // std::cout << "Muon CH = " << analysisTree.muon_chargedHadIso[matchedMuon]
	  // 	    << "  NH = " << analysisTree.muon_neutralHadIso[matchedMuon]
	  // 	    << "  PH = " << analysisTree.muon_photonIso[matchedMuon]
	  // 	    << "  PU = " << analysisTree.muon_puIso[matchedMuon] 
	  // 	    << "  NPV = " << analysisTree.primvertex_count << std::endl;
	  float neutralIsoCorr = neutralIso - 0.5*analysisTree.muon_puIso[matchedMuon];
	  neutralIsoCorr = TMath::Max(float(0),neutralIsoCorr); 
	  float absIso = analysisTree.muon_chargedHadIso[matchedMuon] + neutralIsoCorr;
	  float relIso = absIso/analysisTree.muon_pt[matchedMuon];

	  if (isPrompt) {
	    
	    MuonPtH[muonQ][0][13][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][0][13][iPt]->Fill(muonRecoEta,weight);

	    MuonChargedIsoH[muonQ][0][13][etaBinReco]->Fill(analysisTree.muon_chargedHadIso[matchedMuon],weight);
	    MuonNeutralIsoH[muonQ][0][13][etaBinReco]->Fill(neutralIso,weight);
	    MuonNeutralCorrIsoH[muonQ][0][13][etaBinReco]->Fill(neutralIsoCorr,weight);
	    MuonIsoH[muonQ][0][13][etaBinReco]->Fill(absIso,weight);
	    MuonRelIsoH[muonQ][0][13][etaBinReco]->Fill(relIso,weight);
	    
	    MuonVertexDxH[muonQ][0][13][etaBinReco]->Fill(dVertX,weight);
	    MuonVertexDyH[muonQ][0][13][etaBinReco]->Fill(dVertY,weight);
	    MuonVertexDzH[muonQ][0][13][etaBinReco]->Fill(dVertZ,weight);

	    MuonNvertH[muonQ][0][13]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertEtaH[muonQ][0][13][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertPtH[muonQ][0][13][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	    
	    if (isEventSelected) {

	      MuonPtH[muonQ][0][14][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][0][14][iPt]->Fill(muonRecoEta,weight);
	      
	      MuonChargedIsoH[muonQ][0][14][etaBinReco]->Fill(analysisTree.muon_chargedHadIso[matchedMuon],weight);
	      MuonNeutralIsoH[muonQ][0][14][etaBinReco]->Fill(neutralIso,weight);
	      MuonNeutralCorrIsoH[muonQ][0][14][etaBinReco]->Fill(neutralIsoCorr,weight);
	      MuonIsoH[muonQ][0][14][etaBinReco]->Fill(absIso,weight);
	      MuonRelIsoH[muonQ][0][14][etaBinReco]->Fill(relIso,weight);
	      
	      MuonVertexDxH[muonQ][0][14][etaBinReco]->Fill(dVertX,weight);
	      MuonVertexDyH[muonQ][0][14][etaBinReco]->Fill(dVertY,weight);
	      MuonVertexDzH[muonQ][0][14][etaBinReco]->Fill(dVertZ,weight);
	      
	      MuonNvertH[muonQ][0][14]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertEtaH[muonQ][0][14][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertPtH[muonQ][0][14][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	    }
          }
	  else {
	    MuonPtH[muonQ][1][13][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][1][13][iPt]->Fill(muonRecoEta,weight);

	    if (isEventSelected) {
	      MuonPtH[muonQ][1][14][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][1][14][iPt]->Fill(muonRecoEta,weight);
	    }

	  }

	  if (relIso>isoMuonHighCut&&applyMuonIso) continue;
	  if (relIso<isoMuonLowCut&&applyMuonIso) continue;

	  if (isPrompt) {

	    MuonPtH[muonQ][0][2][etaBinGen]->Fill(muonGenPt,weight);

	    MuonPtH[muonQ][0][4][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][0][4][iPt]->Fill(muonRecoEta,weight);
	    
	    MuonVertexDxH[muonQ][0][4][etaBinReco]->Fill(dVertX,weight);
	    MuonVertexDyH[muonQ][0][4][etaBinReco]->Fill(dVertY,weight);
	    MuonVertexDzH[muonQ][0][4][etaBinReco]->Fill(dVertZ,weight);
	    
	    MuonNvertH[muonQ][0][4]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertEtaH[muonQ][0][4][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	    MuonNvertPtH[muonQ][0][4][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	    
	    if (isEle27fired) {
	      MuonPtH[muonQ][0][5][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonPtH[muonQ][0][7][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonPt23RecoXH->Fill(muonRecoPt,weight);
	      MuonPt8RecoXH->Fill(muonRecoPt,weight);
	      if (isMu23fired) {
		MuonPtH[muonQ][0][6][etaBinReco]->Fill(muonRecoPt,weight);
		MuonPt23RecoPassedXH->Fill(muonRecoPt,weight);
	      }
	      if (isMu8fired) {
		MuonPtH[muonQ][0][8][etaBinReco]->Fill(muonRecoPt,weight);
		MuonPt8RecoPassedXH->Fill(muonRecoPt,weight);
	      }
	    }

	    if (isEventSelected) {

	      MuonPtH[muonQ][0][12][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][0][12][iPt]->Fill(muonRecoEta,weight);

	      MuonVertexDxH[muonQ][0][12][etaBinReco]->Fill(dVertX,weight);
	      MuonVertexDyH[muonQ][0][12][etaBinReco]->Fill(dVertY,weight);
	      MuonVertexDzH[muonQ][0][12][etaBinReco]->Fill(dVertZ,weight);
	      
	      MuonNvertH[muonQ][0][12]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertEtaH[muonQ][0][12][etaBinReco]->Fill(float(analysisTree.primvertex_count),weight);
	      MuonNvertPtH[muonQ][0][12][iPt]->Fill(float(analysisTree.primvertex_count),weight);
	      
	      if (MEt>80) {
		if (isEle27fired) {
		  MuonPt23RecoYH->Fill(muonRecoPt,weight);
		  MuonPt8RecoYH->Fill(muonRecoPt,weight);
		  if (isMu23fired) {
		    MuonPt23RecoPassedYH->Fill(muonRecoPt,weight);
		  }
		  if (isMu8fired) {
		    MuonPt8RecoPassedYH->Fill(muonRecoPt,weight);
		  }
		}
	      }
	      
	    }
	  }
	  else {

	    MuonPtH[muonQ][1][2][etaBinGen]->Fill(muonGenPt,weight);

	    MuonPtH[muonQ][1][4][etaBinReco]->Fill(muonRecoPt,weight);
	    MuonEtaH[muonQ][1][4][iPt]->Fill(muonRecoEta,weight);

	    if (isEventSelected) {
	      MuonPtH[muonQ][1][12][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonEtaH[muonQ][1][12][iPt]->Fill(muonRecoEta,weight);
	    }

	    if (isEle27fired) {
	      MuonPtH[muonQ][1][5][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonPtH[muonQ][1][7][etaBinReco]->Fill(muonRecoPt,weight);
	      MuonPt23RecoH->Fill(muonRecoPt,weight);
	      MuonPt8RecoH->Fill(muonRecoPt,weight);
	      if (isMu23fired) {
		MuonPtH[muonQ][1][6][etaBinReco]->Fill(muonRecoPt,weight);
		MuonPt23RecoPassedH->Fill(muonRecoPt,weight);
	      }
	      if (isMu8fired) {
		MuonPtH[muonQ][1][8][etaBinReco]->Fill(muonRecoPt,weight);
		MuonPt8RecoPassedH->Fill(muonRecoPt,weight);
	      }
	    }

	  }

	}

      }

      // ******************************
      // electron Id with tag-and-probe
      // ******************************
      if (analysisTree.hltriggerresults_second[0]==1) { // Single Ele27 trigger
	for (int iQ=0; iQ<2; ++iQ) {

	  int indexTag = indexTagPosE;
	  int indexProbe = indexProbeNegE;
	  if (iQ==1) {
	    indexTag = indexTagNegE;
	    indexProbe = indexProbePosE;
	  }
	  

	  if (indexTag>=0&&indexProbe>=0) { 

	    TLorentzVector tagP4; tagP4.SetXYZM(analysisTree.electron_px[indexTag],
						analysisTree.electron_py[indexTag],
						analysisTree.electron_pz[indexTag],
						electronMass);
	    TLorentzVector probeP4; probeP4.SetXYZM(analysisTree.electron_px[indexProbe],
						    analysisTree.electron_py[indexProbe],
						    analysisTree.electron_pz[indexProbe],
						    electronMass);
	    
	    // artificial matching to Z->EE events 
	    //	    if (!isZToEE) continue;
	    //	    if (!isZToEEopen) continue;
	    //	    float dPMin = 1;
	    //	    for (unsigned int genE=0; genE<indexE.size(); ++genE) {
	    //	      int index = indexE.at(genE);
	    //	      TLorentzVector genP4; genP4.SetXYZM(analysisTree.genparticles_px[index],
	    //						  analysisTree.genparticles_py[index],
	    //						  analysisTree.genparticles_pz[index],
	    //					  electronMass);

	    //	      float dP = (genP4-probeP4).P();
	    //	      float dPoverP = dP/genP4.P();
	    //	      if (dPoverP<dPMin) {
	    //		dPMin = dPoverP;
	    //	      }
	    //	    }
	    //	    if (dPMin>0.05) continue;

	    float deltaREE = deltaR(analysisTree.electron_eta[indexTag],
				    analysisTree.electron_phi[indexTag],
				    analysisTree.electron_eta[indexProbe],
				    analysisTree.electron_phi[indexProbe]);
	    float deltaPhiEE = dPhiFrom2P(analysisTree.electron_px[indexTag],
					  analysisTree.electron_py[indexTag],
					  analysisTree.electron_px[indexProbe],
					  analysisTree.electron_py[indexProbe]);

	    if (deltaREE<dRLeptonsCut) continue;
	    if (deltaPhiEE<dPhiLeptonsLowCut) continue;
	    if (deltaPhiEE>dPhiLeptonsHighCut) continue;

	    float mass = (tagP4+probeP4).M();
	    bool probePassed = true;
	    if (fabs(analysisTree.electron_dxy[indexProbe])>dxyElectronCut&&applyElectronIPcuts) probePassed = false;
	    if (fabs(analysisTree.electron_dz[indexProbe])>dzElectronCut&&applyElectronIPcuts) probePassed = false;
	    float neutralIso = 
	      analysisTree.electron_neutralHadIso[indexProbe] + 
	      analysisTree.electron_photonIso[indexProbe] - 
	      0.5*analysisTree.electron_puIso[indexProbe];
	    neutralIso = TMath::Max(float(0),neutralIso); 
	    float absIso = analysisTree.electron_chargedHadIso[indexProbe] + neutralIso;
	    float relIso = absIso/analysisTree.electron_pt[indexProbe];
	    if (relIso>isoElectronHighCut&&applyElectronIso) probePassed = false;
	    if (relIso<isoElectronLowCut&&applyElectronIso) probePassed = false;
	    bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[indexProbe],
						    analysisTree.electron_mva_id_nontrigPhys14[indexProbe]);
	    if (!electronMvaId&&applyElectronId) probePassed = false;
	    if (!analysisTree.electron_pass_conversion[indexProbe]&&applyElectronId) probePassed = false;
	    if (analysisTree.electron_nmissinginnerhits[indexProbe]!=0&&applyElectronId) probePassed = false; 
	    int etaBin = binNumber(fabs(analysisTree.electron_eta[indexProbe]),nEtaBins,etaBins);
	    int ptBin = binNumber(analysisTree.electron_pt[indexProbe],nPtBins,ptBins);
	    ElectronPtH[iQ][0][9][etaBin]->Fill(analysisTree.electron_pt[indexProbe],weight);
	    if (probePassed) {
	      massEEId[iQ][etaBin][ptBin][0]->Fill(mass,weight);
	      ElectronPtH[iQ][0][10][etaBin]->Fill(analysisTree.electron_pt[indexProbe],weight);
	    }
	    else
	      massEEId[iQ][etaBin][ptBin][1]->Fill(mass,weight);
	  }
	}
      }

      // **************************
      // muon Id with tag-and-probe
      // **************************
      if (analysisTree.hltriggerresults_second[2]==1||analysisTree.hltriggerresults_second[3]==1) { // Single IsoMu24 || IsoTkMu24 trigger
	for (int iQ=0; iQ<2; ++iQ) {

	  int indexTag = indexTagPosMu;
	  int indexProbe = indexProbeNegMu;
	  if (iQ==1) {
	    indexTag = indexTagNegMu;
	    indexProbe = indexProbePosMu;
	  }
	  

	  if (indexTag>=0&&indexProbe>=0) { 
	    
	    TLorentzVector tagP4; tagP4.SetXYZM(analysisTree.muon_px[indexTag],
						analysisTree.muon_py[indexTag],
						analysisTree.muon_pz[indexTag],
						muonMass);
	    TLorentzVector probeP4; probeP4.SetXYZM(analysisTree.muon_px[indexProbe],
						    analysisTree.muon_py[indexProbe],
						    analysisTree.muon_pz[indexProbe],
						    muonMass);
	    // artificial matching to Z->MM events
	    //	    if (!isZToMM) continue;
	    //	    if (!isZToMMopen) continue;
	    //	    float dPMin = 1;
	    //	    for (unsigned int genMu=0; genMu<indexMu.size(); ++genMu) {
	    //	      int index = indexMu.at(genMu);
	    //	      TLorentzVector genP4; genP4.SetXYZM(analysisTree.genparticles_px[index],
	    //						  analysisTree.genparticles_py[index],
	    //						  analysisTree.genparticles_pz[index],
	    //						  muonMass);

	    //	      float dP = (genP4-probeP4).P();
	    //	      float dPoverP = dP/genP4.P();
	    //	      if (dPoverP<dPMin) {
	    //		dPMin = dPoverP;
	    //	      }
	    //	    }
	    //	    if (dPMin>0.05) continue;

	    float deltaRMuMu = deltaR(analysisTree.muon_eta[indexTag],
				      analysisTree.muon_phi[indexTag],
				      analysisTree.muon_eta[indexProbe],
				      analysisTree.muon_phi[indexProbe]);
	    float deltaPhiMuMu = dPhiFrom2P(analysisTree.muon_px[indexTag],
					    analysisTree.muon_py[indexTag],
					    analysisTree.muon_px[indexProbe],
					    analysisTree.muon_py[indexProbe]);

	    if (deltaRMuMu<dRLeptonsCut) continue;
	    if (deltaPhiMuMu<dPhiLeptonsLowCut) continue;
	    if (deltaPhiMuMu>dPhiLeptonsHighCut) continue;

	    float mass = (tagP4+probeP4).M();
	    bool probePassed = true;
	    if (fabs(analysisTree.muon_dxy[indexProbe])>dxyMuonCut&&applyMuonIPcuts) probePassed = false;
	    if (fabs(analysisTree.muon_dz[indexProbe])>dzMuonCut&&applyMuonIPcuts) probePassed = false;
	    float neutralIso = 
	      analysisTree.muon_neutralHadIso[indexProbe] + 
	      analysisTree.muon_photonIso[indexProbe] - 
	      0.5*analysisTree.muon_puIso[indexProbe];
	    neutralIso = TMath::Max(float(0),neutralIso); 
	    float absIso = analysisTree.muon_chargedHadIso[indexProbe] + neutralIso;
	    float relIso = absIso/analysisTree.muon_pt[indexProbe];
	    if (relIso>isoMuonHighCut&&applyMuonIso) probePassed = false;
	    if (relIso<isoMuonLowCut&&applyMuonIso) probePassed = false;
	    if (applyMuonId==1 && !analysisTree.muon_isLoose[indexProbe]) probePassed = false;
	    if (applyMuonId==2 && !analysisTree.muon_isMedium[indexProbe]) probePassed = false;
	    if (applyMuonId==3 && !analysisTree.muon_isTight[indexProbe]) probePassed = false;
	    int etaBin = binNumber(fabs(analysisTree.muon_eta[indexProbe]),nEtaBins,etaBins);
	    int ptBin = binNumber(analysisTree.muon_pt[indexProbe],nPtBins,ptBins);
	    MuonPtH[iQ][0][9][etaBin]->Fill(analysisTree.muon_pt[indexProbe],weight);
	    if (probePassed) {
	      massMuMuId[iQ][etaBin][ptBin][0]->Fill(mass,weight);
	      MuonPtH[iQ][0][10][etaBin]->Fill(analysisTree.muon_pt[indexProbe],weight);
	    }
	    else
	      massMuMuId[iQ][etaBin][ptBin][1]->Fill(mass,weight);
	  }
	}
      }

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



