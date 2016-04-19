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
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include "TRandom.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#include "DesyTauAnalyses/NTupleMaker/interface/rochcor2015.h"
#include "DesyTauAnalyses/NTupleMaker/interface/RoccoR.h"

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

float topPtWeight(float pt1,
		  float pt2) {

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter
  
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  if (pt1>400) w1 = 1;
  if (pt2>400) w2 = 1;

  return TMath::Sqrt(w1*w2);

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
 	//const bool applyInclusiveSelection = cfg.get<bool>("ApplyInclusiveSelection");
	const bool isData = cfg.get<bool>("IsData");
	
	//Run-lumi selector
	const string jsonFile = cfg.get<string>("jsonFile");
	
	//Apply corrections****
 	//pile up reweighting
 	const bool applyPUreweighting_vertices = cfg.get<bool>("ApplyPUreweighting_vertices");
	const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");
 
	const bool applyMEtRecoilCorrections = cfg.get<bool>("ApplyMEtRecoilCorrections");
	const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

	const bool applyTopPtReweighting = cfg.get<bool>("ApplyTopPtReweighting");
	
	const bool applyRochCorr = cfg.get<bool>("ApplyRochCorr");
 	
	
	//ztotautautomumu selection
	const bool  applyTauTauSelection = cfg.get<bool>("ApplyTauTauSelection");
	const bool  selectZToTauTauMuMu = cfg.get<bool>("SelectZToTauTauMuMu");
 
	// kinematic cuts on muons
	const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
 	const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
 	const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
 	const float etaMuonLowCut = cfg.get<float>("etaMuonLowCut");
 	const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
 	const float dzMuonCut      = cfg.get<float>("dzMuonCut");
	const float isoMuonCut     = cfg.get<float>("isoMuonCut");
	
	// vertex cuts
 	const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
 	const float zVertexCut     = cfg.get<float>("ZVertexCut");
	const float dVertexCut     = cfg.get<float>("DVertexCut");
	
	// topological cuts
 	const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
 	const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
	const bool oppositeSign    = cfg.get<bool>("OppositeSign");
	
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
  	const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");

	TString BTagDiscriminator(bTagDiscriminator);
  
	const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  
	// trigger
  	const bool doubleMuonTrigger = cfg.get<bool>("DoubleMuonTrigger");
	
  	const string muonTriggerName = cfg.get<string>("MuonTriggerName");
  	const string muonFilterName = cfg.get<string>("MuonFilterName");
  	const string doubleMuonHighPtFilterName = cfg.get<string>("DoubleMuonHighPtFilterName");
  	const string doubleMuonLowPtFilterName = cfg.get<string>("DoubleMuonLowPtFilterName");
	
	TString MuonTriggerName(muonTriggerName);
  	TString MuonFilterName(muonFilterName);
  	TString DoubleMuonHighPtFilterName(doubleMuonHighPtFilterName);
  	TString DoubleMuonLowPtFilterName(doubleMuonLowPtFilterName);
	
  	const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
  	const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");
	const string MuonTrigFile = cfg.get<string>("MuonTrigEff");

  	const string recoilFileName   = cfg.get<string>("RecoilFileName");
	TString RecoilFileName(recoilFileName);
	
    const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
    TString RecoilMvaFileName(recoilMvaFileName);		
	
	//const string recoilPuppiFileName   = cfg.get<string>("RecoilPuppiFileName");	
	//TString RecoilPuppiFileName(recoilPuppiFileName);
	
	//Vertex distribution filenames and histname for pileup calculation
	const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  	const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  	const string vertHistName     = cfg.get<string>("VertexHistName");
	
	// Run range
  	const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  	const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
	
	const bool fillBDTNTuple = cfg.get<bool>("FillBDTNTuple");
	// **** end of configuration
	
  	string cmsswBase = (getenv ("CMSSW_BASE"));
	
 	TString fullDir = TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/");
	
  	// int nDPtBins = 9;
  	// float DPtBins[10] = {0,10,15,20,25,30,40,50,60,1000};
  	// TString dPtBins[9] = {"Pt0to10",
	// "Pt10to15",
	// "Pt15to20",
	// "Pt20to25",
	// "Pt25to30",
	// "Pt30to40",
	// "Pt40to50",
	// "Pt50to60",
	// "PtGt60"};
	
	// int nDEtaBins = 4;
  	// float DEtaBins[5]= {0,0.9,1.2,2.1,2.4};
  	// TString dEtaBins[4] = {"Eta0to0p9",
	// "Eta0p9to1p2",
	// "Eta1p2to2p1",
	// "Eta2p1to2p4"};
	
  	// //TH1D * dxyMu1[4][9];
  	// //TH1D * dxyMu2[4][9];
  	// //TH1D * dzMu1[4][9];
  	// //TH1D * dzMu2[4][9];
  
  	// int nJetBins = 3;
  	// int nZPtBins = 5;
  	// float zPtBins[6] = {0,10,20,30,50,1000};
  	// float jetBins[4] = {-0.5,0.5,1.5,2.5};

  	// TString NJetBins[3] = {"NJet0","NJet1","NJetGe2"};
  	// TString ZPtBins[5] = {"Pt0to10H",
	// "Pt10to20H",
	// "Pt20to30H",
	// "Pt30to50H",
	// "PtGt50H"};
	
	// // Saving Z pt bins
	// TH1D * ZPtBinsH = new TH1D("ZPtBinsH","ZPtBinsH",nZPtBins,zPtBins);
	// for (int iB=0; iB<nZPtBins; ++iB) 
	// 	ZPtBinsH->GetXaxis()->SetBinLabel(iB+1,ZPtBins[iB]);
	
	// // Saving jet bins
	// TH1D * JetBinsH = new TH1D("JetBinsH","JetBinsH",nJetBins,jetBins);
	// for (int iB=0; iB<nJetBins; ++iB)
	// 	JetBinsH->GetXaxis()->SetBinLabel(iB+1,NJetBins[iB]);
	
  	// TH1D * recoilZParalH[3];
  	// TH1D * recoilZPerpH[3];
  	// TH1D * recoilZParal_Ptbins_nJetsH[3][5];
  	// TH1D * recoilZPerp_Ptbins_nJetsH[3][5];
	  
  	// TH1D * ptRatio_nJetsH[3];
	
	// TString RecoilZParal("recoilZParal_");
  	// TString RecoilZPerp("recoilZPerp_");
	
	int nPtBins = 7;
        float ptBins[8] = {10,15,20,25,30,40,60,1000};
	
	 int nEtaBins = 3;
	 float etaBins[4] = {-0.1,0.9,1.2,2.4}; 
	TString PtBins[7] = {"Pt10to15",
	"Pt15to20",
	"Pt20to25",
	"Pt25to30",
	"Pt30to40",
	"Pt40to60",
	"PtGt60"};
	
	TString EtaBins[3] = {"EtaLt0p9",
	"Eta0p9to1p2",
	"EtaGt1p2"};
	
	// file name and tree name
  	std::string rootFileName(argv[2]);
  	std::ifstream fileList(argv[2]);
  	std::ifstream fileList0(argv[2]);
  	std::string ntupleName("makeroottree/AC1B");
	
	TString TStrName(rootFileName);
	std::cout <<TStrName <<std::endl;
	
	TFile * file = new TFile(TStrName+TString(".root"),"recreate");
	file->cd("");
	
	TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
   	TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
	
	// Histograms after selecting unique dimuon pair
   	TH1D * ptLeadingMuSelH = new TH1D("ptLeadingMuSelH","",100,0,200);
   	TH1D * ptTrailingMuSelH = new TH1D("ptTrailingMuSelH","",100,0,200);
   	TH2F * ptScatter =new TH2F("ptScatter","",100,0,200,100,0,200);
   	TProfile2D *hprof2D_pt = new TProfile2D("hprof2D_pt","",100,0,200,100,0,200);
   	TH1D * etaLeadingMuSelH = new TH1D("etaLeadingMuSelH","",50,-2.5,2.5);
   	TH1D * etaTrailingMuSelH = new TH1D("etaTrailingMuSelH","",50,-2.5,2.5);
	TH1D * h_dimuonPt = new TH1D ("dimuonPt","",100,0,200);
   	TH1D * massSelH = new TH1D("massSelH","",200,0,200);
	TH1D * mass_sv = new TH1D("mass_sv","",200,0,200);
	//TH1D * massGenSelH = new TH1D("massGenSelH","",200,0,200);
	TH1D * massSelGen1H = new TH1D("massSelGen1H","",200,0,200);
	TH1D * dimuonMass_dca = new TH1D ("dimuonMass_dca","",200,0,200);
   	TH1D * metSelH  = new TH1D("metSelH","",200,0,400);
	TH1D * mvametSelH = new TH1D("mvametSelH","",200,0,400);

   	TH1D * nJets30SelH    = new TH1D("nJets30SelH","",11,-0.5,10.5);
   	TH1D * nJets30etaCutSelH = new TH1D("nJets30etaCutSelH","",11,-0.5,10.5);
   	TH1D * nJets20SelH    = new TH1D("nJets20SelH","",11,-0.5,10.5);
   	TH1D * nJets20etaCutSelH = new TH1D("nJets20etaCutSelH","",11,-0.5,10.5);

   	TH1D * HT30SelH       = new TH1D("HT30SelH","",50,0,500);
   	TH1D * HT30etaCutSelH = new TH1D("HT30etaCutSelH","",50,0,500);
   	TH1D * HT20SelH       = new TH1D("HT20SelH","",50,0,500);
   	TH1D * HT20etaCutSelH = new TH1D("HT20etaCutSelH","",50,0,500);
	
	//Discriminiant histos
   	TH1D * h_dimuonEta = new TH1D("dimuonEta","",50,-6,+6);//21Aug
   	TH1D * h_ptRatio = new TH1D ("ptRatio","",50,0,1);
   	TH1D * h_ptRatio_test = new TH1D ("ptRatio_test","",50,0,1);
	
   	TH1D * h_dxy_muon1 =new TH1D ("dxy_muon1","",50,-0.02,0.02);
   	TH1D * h_dxy_muon2 =new TH1D ("dxy_muon2","",50,-0.02,0.02);
   	TH1D * h_dz_muon1 = new TH1D ("dz_muon1","",50,-0.1,0.1);
   	TH1D * h_dz_muon2 = new TH1D ("dz_muon2","",50,-0.1,0.1);
  
   	TH1D * h_dcaSigdxy_muon1 = new TH1D ("dcaSigdxy_muon1","",50,-4,4);
   	TH1D * h_dcaSigdxy_muon2 = new TH1D ("dcaSigdxy_muon2","",50,-4,4);
   	TH1D * h_dcaSigdz_muon1 = new TH1D ("dcaSigdz_muon1","",50,-4,4);
   	TH1D * h_dcaSigdz_muon2 = new TH1D ("dcaSigdz_muon2","",50,-4,4);
	
	TH1D * h_dcaSig2Mu2D = new TH1D ("dcaSig2Mu2D","",50,-4,4);
	TH1D * h_dcaSig2Mu2D_Jet0 = new TH1D ("dcaSig2Mu2D_Jet0","",50, -4,4);
	TH1D * h_dcaSig2Mu3D = new TH1D ("dcaSig2Mu3D","",50,-4,4);
	TH1D * h_sig2Mu2D = new TH1D ("sig2Mu2D","",50,0,100);
	TH1D * h_sig2Mu3D = new TH1D ("sig2Mu3D","",50,0,3000);
	
   	TH1D * h_phi_leadingMu_MET =new TH1D ("phi_leadingMu_MET","",50,0,3);
   	TH1D * h_phi_trailingMu_MET =new TH1D ("phi_trailingMu_MET","",50,0,3);
	TH1D * h_phi_PosMu_MET =new TH1D ("phi_PosMu_MET","",50,0,3);
   	TH1D * h_phi_TwoMu = new TH1D ("phi_TwoMu","",50,0,3); 
	
   	TH1D * h_DZeta = new TH1D("DZeta","",100,-400,200);
   	TH1D * h_Njets = new TH1D("Njets","",10,0,10);
   	TH1D * NumberOfVerticesH = new TH1D("NumberOfVerticesH","",50,-0.5,50.5);
	
   	// TString scales[21] = {"M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","0",
	// "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};
	
	// TH1D * massSelScaleH[21];
  	// TH1D * metSelScaleH[21];
   	// for (int iScale=0; iScale<21; ++iScale) {
	// 	massSelScaleH[iScale] = new TH1D("massSel"+scales[iScale]+"H","",200,0,200);
     	// metSelScaleH[iScale] = new TH1D("metSel"+scales[iScale]+"H","",200,0,400);
	// }
   	TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);
   	TH1D * MuSF_IdIso_Mu2H = new TH1D("MuIdIsoSF_Mu2H", "MuIdIsoSF_Mu2", 100, 0.5,1.5);
	
   	// TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
   	// TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);
	
	 //BDT histos
	TH1F * histBdt =  new TH1F("MVA_BDT", "MVA_BDT",800 , -0.8, 0.8);
	TH1F * histMva_massCut =  new TH1F("MVA_BDT_massCut", "MVA_BDT",1000 , -1.0, 1.0);
	TH1F * dimuonMass_BDT0p5 =new TH1F("dimuonMass_BDT0p5", "", 200,0,200);
	
	//Booking Invariant mass for different BDT cutt
	// const int iCut =15;
	// const char Nname = 100;
	// TH1F* InvMass[iCut];
	// char name[Nname];
	// for (int i=0;i<iCut;i++){
	//   sprintf(name,"InvMass_%i",i);
	//   InvMass[i]=new TH1F(name,"",200,20,200);
	// }

	 //Booking variables in Tree for BDT training  
	Float_t n_genWeight;
	Float_t n_dimuonEta;
   	Float_t n_ptRatio;
   	Float_t n_dcaSigdxy1;
   	Float_t n_dcaSigdxy2;
	Float_t n_dcaSigdz1;
	Float_t n_dcaSigdz2;
	Float_t n_phiangle;
   	Float_t n_twomuPhi;
   	Float_t n_dxy_muon1;
   	Float_t n_dxy_muon2;
   	Float_t n_dz_muon1;
   	Float_t n_dz_muon2;
   	Float_t n_dcaSigdxy_mu1;
   	Float_t n_dcaSigdxy_mu2;
   	Float_t n_dcaSigdz_mu1;
   	Float_t n_dcaSigdz_mu2;
	Float_t n_dcaSig2Mu2D;
	Float_t n_dcaSig2Mu3D;
	Float_t n_sig2Mu2D;
	Float_t n_sig2Mu3D;
   	Float_t n_MissingEt;
	Float_t n_DZeta;
   	Float_t n_dimuonMass;
   	Float_t n_met;
	Float_t n_mvamet;
	Float_t n_mvamet_ex;
	Float_t n_mvamet_ey;
	Float_t n_covmet_xx;
	Float_t n_covmet_xy;
	Float_t n_covmet_yy;
	Float_t n_leadingPt;
	Float_t n_trailingPt;
	Float_t n_leadingEta;
	Float_t n_trailingEta;
	Float_t n_leadingPhi;
	Float_t n_trailingPhi;
	Float_t n_jets;
	Float_t n_noOfvertices;
	Bool_t n_genAccept;
	Float_t n_genZ;
	Float_t n_mvaBDT;
	Float_t n_m_sv;
	Float_t n_pt_sv; 
	Float_t	n_eta_sv;
	Float_t	n_phi_sv; 

	// TTree * TW = new TTree("TW","Weights");
   	// TW->Branch("genWeight",&n_genWeight,"n_genWeight/F");
	
	TTree * T = new TTree("T","Discriminant variables for BDT");
	T->Branch("dimuonEta",&n_dimuonEta,"n_dimuonEta/F");
   	T->Branch("ptRatio",&n_ptRatio,"n_ptRatio/F");
   	T->Branch("dxy_muon1",&n_dxy_muon1,"n_dxy_muon1/F");
   	T->Branch("dxy_muon2",&n_dxy_muon2,"n_dxy_muon2/F");
   	T->Branch("dz_muon1",&n_dz_muon1,"n_dz_muon1/F");
   	T->Branch("dz_muon2",&n_dz_muon2,"n_dz_muon2/F");
   	T->Branch("dcaSigdxy_muon1",&n_dcaSigdxy1,"n_dcaSigdxy1/F");
   	T->Branch("dcaSigdxy_muon2",&n_dcaSigdxy2,"n_dcaSigdxy2/F");
   	T->Branch("dcaSigdz_muon1",&n_dcaSigdz1,"n_dcaSigdz1/F");
   	T->Branch("dcaSigdz_muon2",&n_dcaSigdz2,"n_dcaSigdz2/F");
	T->Branch("dcaSig2Mu2D", &n_dcaSig2Mu2D, "n_dcaSig2Mu2D/F");
	T->Branch("dcaSig2Mu3D", &n_dcaSig2Mu3D, "n_dcaSig2Mu3D/F");
	T->Branch("sig2Mu2D", &n_sig2Mu2D, "n_sig2Mu2D/F");
	T->Branch("sig2Mu3D", &n_sig2Mu3D, "n_sig2Mu3D/F");
   	T->Branch("MissingET",&n_MissingEt,"n_MissingEt/F");
	T->Branch("phi_PosMu_MET",&n_phiangle, "n_phiangle/F");
   	T->Branch("phi_TwoMu",&n_twomuPhi,"n_twomuPhi/F");
	T->Branch("DZeta",&n_DZeta,"n_DZeta/F");
   	T->Branch("genWeight",&n_genWeight,"n_genWeight/F");
   	T->Branch("dimuonMass",&n_dimuonMass,"n_dimuonMass/F");
	T->Branch("met",&n_met,"n_met/F");
	T->Branch("mvamet_ex",&n_mvamet_ex,"n_mvamet_ex/F");
	T->Branch("mvamet_ey",&n_mvamet_ey,"n_mvamet_ey/F");
	T->Branch("mvamet",&n_mvamet,"n_mvamet/F");
	T->Branch("covmetxx",&n_covmet_xx,"n_covmet_xx/F");
	T->Branch("covmetxy",&n_covmet_xy,"n_covmet_xy/F");
	T->Branch("covmetyy",&n_covmet_yy,"n_covmet_yy/F");
	T->Branch("leadingPt", &n_leadingPt, "n_leadingPt/F");
	T->Branch("trailingPt", &n_trailingPt, "n_trailingPt/F");
	T->Branch("leadingEta", &n_leadingEta, "n_leadingEta/F");
	T->Branch("trailingEta", &n_trailingEta, "n_trailingEta/F");
	T->Branch("leadingPhi", &n_leadingPhi, "n_leadingPhi/F");
	T->Branch("trailingPhi", &n_trailingPhi, "n_trailingPhi/F");
	T->Branch("jets", &n_jets,"n_jets/F");
	T->Branch("noOfvertices", &n_noOfvertices, " n_noOfvertices/F");
	T->Branch("genAccept", &n_genAccept, "genAccept/O");
	T->Branch("genZ", &n_genZ, "n_genZ/F");
	T->Branch("mvaBDT", &n_mvaBDT, "&n_mvaBDT/F");
	T->Branch("m_sv", &n_m_sv, "n_m_sv/F");
	T->Branch("pt_sv", &n_pt_sv, "n_pt_sv/F");
	T->Branch("eta_sv", &n_eta_sv, "n_eta_sv/F");
	T->Branch("phi_sv", &n_phi_sv, "n_phi_sv/F");


	//This loads the library
	TMVA::Tools::Instance();
	
	//Create TMVA Reader Object
	TMVA::Reader *reader = new TMVA::Reader("!V:!Color");
	
	//create set of variables as declared in weight file and declared them to reader
	reader->AddVariable( "dimuonEta",&n_dimuonEta);
	reader->AddVariable( "phi_PosMu_MET",&n_phiangle);
	reader->AddVariable("DZeta",&n_DZeta);
	reader->AddVariable( "mvamet",&n_mvamet);
	
	//BookMethod
	reader->BookMVA("BDT", "rootFiles/TMVA/weights/TMVA_Roch_24032016_BDT.weights.xml");
	
	
	// TH1D * ZMassEtaPtPass[3][7];
   	// TH1D * ZMassEtaPtFail[3][7];
	
	// for (int iEta=0; iEta<nEtaBins; ++iEta) {
	// 	for (int iPt=0; iPt<nPtBins; ++iPt) {
	// 		ZMassEtaPtPass[iEta][iPt] = new TH1D("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",60,60,120);
	// 		ZMassEtaPtFail[iEta][iPt] = new TH1D("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",60,60,120);
	// 	}
	// }
	
	// for (int iBin=0; iBin<nJetBins; ++iBin) {
	// 	recoilZParalH[iBin] = new TH1D(RecoilZParal+NJetBins[iBin]+"H","",100,-200,200);
	// 	recoilZPerpH[iBin] = new TH1D(RecoilZPerp+NJetBins[iBin]+"H","",100,-200,200);
	// }
	
	// for (int iJets=0; iJets<nJetBins; ++iJets) {
	// 	for (int iPtBins=0; iPtBins<nZPtBins; ++iPtBins){
	// 		recoilZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);
	// 		recoilZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);
	// 	}
	// }
	
	// for (int iJetBin=0; iJetBin<nJetBins; ++iJetBin) {
	// 	ptRatio_nJetsH[iJetBin]= new TH1D("ptRatio_nJets"+NJetBins[iJetBin]+"H","",100,0,1);
	// }
	
	//for (int iEta=0; iEta<nDEtaBins; ++iEta){
	//	for (int iPt=0; iPt<nDPtBins; ++iPt){
			//std::cout<< "iEta =" <<iEta<< " and iPt =" << iPt<<std::endl;
	//		dxyMu1[iEta][iPt]= new TH1D("dxyMu1_"+dEtaBins[iEta]+"_"+dPtBins[iPt],"",50,-0.02,0.02);
	//		dxyMu2[iEta][iPt]= new TH1D("dxyMu2_"+dEtaBins[iEta]+"_"+dPtBins[iPt],"",50,-0.02,0.02);
	//		dzMu1[iEta][iPt]= new TH1D("dzMu1_"+dEtaBins[iEta]+"_"+dPtBins[iPt],"",50,-0.2,0.2);
	//		dzMu2[iEta][iPt]= new TH1D("dzMu2_"+dEtaBins[iEta]+"_"+dPtBins[iPt],"",50,-0.2,0.2);
	//	}
	//}
	
	// Pile UP reweighting Options (vertices or ofiicial)
	
	if (applyPUreweighting_vertices and applyPUreweighting_official)
	  {std::cout<<"ERROR: Choose only ONE PU reweighting method (vertices or official, not both!!)"<<std::endl; exit(-1);}
	
	//reweighting with vertices
	
	// reading vertex weights
	TFile * fileDataNVert = new TFile(fullDir+vertDataFileName);
	TFile * fileMcNVert   = new TFile(fullDir+vertMcFileName);
	
	bool vertexFileFound = true;
	if (fileDataNVert->IsZombie()) {
		//    std::cout << "File " << fullDir << vertDataFileName << " is not found" << std::endl;
		vertexFileFound = false;
	}
	
	if (fileMcNVert->IsZombie()) {
		//    std::cout << "File " << fullDir << vertMcFileName << " is not found" << std::endl;
		vertexFileFound = false;
	}
   	if (!vertexFileFound) 
		exit(-1);
	
	TH1D * vertexDataH = (TH1D*)fileDataNVert->Get(TString(vertHistName));
   	TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get(TString(vertHistName));
   	if (vertexDataH==NULL||vertexMcH==NULL) {
		std::cout << "Vertex distribution histogram " << vertHistName << " is not found" << std::endl;
		exit(-1);
	}
	
	float normVertexData = vertexDataH->GetSumOfWeights();
   	float normVertexMc   = vertexMcH->GetSumOfWeights();
	
   	vertexDataH->Scale(1/normVertexData);
   	vertexMcH->Scale(1/normVertexMc);
	
	// PU reweighting official receipe
	//Initialize Pile up object   
   	PileUp * PUofficial = new PileUp();
	
	if (applyPUreweighting_official) {
		TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Feb02.root","read");
	   	TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Fall15_PU25_V1.root", "read");
	   	TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
	   	TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
	   	PUofficial->set_h_data(PU_data);
	   	PUofficial->set_h_MC(PU_mc);
	}
	
	// Lepton Scale Factors 
	
   	ScaleFactor * SF_muonIdIso = new ScaleFactor(); 
   	SF_muonIdIso = new ScaleFactor();
   	//std::cout<<"test1"<<std::endl;
   	SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
	ScaleFactor * SF_muonTrig;
	SF_muonTrig = new ScaleFactor();
	SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
   	ScaleFactor *SF_muon17 = new ScaleFactor();
   	//std::cout<<"test12"<<std::endl;
   	SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));
   	ScaleFactor *SF_muon8 = new ScaleFactor();
   	//std::cout<<"test13"<<std::endl;
   	SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));
	
	//loading recoil resolution
	RecoilCorrector recoilPFMetCorrector("HTT-utilities/RecoilCorrections/data/"+RecoilFileName);
	RecoilCorrector recoilMvaMetCorrector("HTT-utilities/RecoilCorrections/data/"+RecoilMvaFileName);
	//   	RecoilCorrector recoilPuppiMetCorrector("HTT-utilities/RecoilCorrections/data/recoilPuppiMEt_amcatnlo.root");
	
	//Rochester Correction Z mass
	rochcor2015 *rmcor = new rochcor2015();
	//	rmcor->init(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/RoccoR_13TeV.txt");
	//**SVFit**
	edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
	TH1::AddDirectory(false);  
	TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
	
	
   	//std::cout<<"test14"<<std::endl;
   	int nFiles = 0;
   	int nEvents = 0;
   	int selEventsAllMuons = 0;
   	int selEventsIdMuons = 0;
   	int selEventsIsoMuons = 0;
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
			//std::cout << "number of events   = "<< nEvents<<std::endl; 
			nEvents++;
			if (nEvents%10000==0) 
				cout << "      processed " << nEvents << " events" << endl;
			
			float weight = 1;
			
			float pfmet_ex = analysisTree.pfmet_ex;
			float pfmet_ey = analysisTree.pfmet_ey;
			float pfmet_phi = analysisTree.pfmet_phi;
			float pfmet = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
			float  pfmetcorr_ex =0.0;
			float  pfmetcorr_ey = 0.0;
			
			//std::cout << "before Roch correction "<<std::endl;

			if (applyRochCorr){
			  for (unsigned int iM = 0; iM<analysisTree.muon_count; ++iM){
			    TLorentzVector mu; //TLorentzVeccor object of the reconstructed muon.
			    mu.SetPtEtaPhiM(analysisTree.muon_pt[iM],
					 analysisTree.muon_eta[iM],
					 analysisTree.muon_phi[iM],
					 muonMass);

			    float charge = analysisTree.muon_charge[iM];
			    float qter = 1.0;
			    
			    if (charge > 0.5) {
			      qter = 1.0;
			      // std::cout<< "charge = "<<charge<< "/t qter =" << qter<<std::endl;
			    }
			    else{ 
			      qter = -1.0;
			      // std::cout<< "charge = "<<charge<< "/t qter =" << qter<<std::endl;
			    }
			    float ntrk =0.0;
			    float runopt = 0.0;
			    //			    std::cout << "Uncorrected : " << mu.Px() << "  "  
			    //				      << mu.Py() << "  "
			    //				      << mu.Pz() << "  " << std::endl;
			    if (!isData){
			      rmcor->momcor_mc(mu,charge,ntrk=0, qter);
			    }
			    if (isData){
			      rmcor->momcor_data(mu, charge, runopt=0, qter); 
			    }
			    //			    std::cout << "Corrected : " << mu.Px() << "  "  
			    //				      << mu.Py() << "  "
			    //				      << mu.Pz() << "  " << std::endl;
			    analysisTree.muon_px[iM] = mu.Px();
			    analysisTree.muon_py[iM] = mu.Py();
			    analysisTree.muon_pz[iM] = mu.Pz();
			    analysisTree.muon_pt[iM] = mu.Pt();
			    analysisTree.muon_eta[iM] = mu.Eta();
			    analysisTree.muon_phi[iM] = mu.Phi();
			    //			    analysisTree.muon_e[iM] = mu.E();
			    
			  }
			}
			//std::cout << "After roch correction "<<std::endl;

			//------------------------------------------------
			// if (!isData) {
			// 	n_genWeight = 1;
				
			// 	if (analysisTree.genweight<0)
			// 		n_genWeight = -1;
				
			// 	//	std::cout << "GenWeight = " << analysisTree.genweight << std::endl;
			// 	weight *= n_genWeight;
			// 	TW->Fill();
			// }
			
			histWeightsH->Fill(float(0),weight);
			
			//weights
			
			//Float_t mcweight = analysisTree.genweight;
			Float_t PUweight = 1;
			Float_t trigweight_1 = 1;
			Float_t trigweight_2 = 1;
			Float_t trigweight = 1;
			Float_t idweight_1 = 1;
			Float_t idweight_2 = 1;
			Float_t isoweight_1 = 1;
			Float_t isoweight_2 = 1;
			Float_t effweight = 1;
			Float_t fakeweight = 1;
			Float_t embeddedWeight = 1;
			Float_t signalWeight = 1;
			//Float_t weight = 1;
			
			TLorentzVector genZ; genZ.SetXYZT(0,0,0,91.2); 
			TLorentzVector genV; genV.SetXYZT(0,0,0,0);
			TLorentzVector genL; genL.SetXYZT(0,0,0,0);
			if (!isData) {
			  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
			    //	  cout << igen << "   pdgId = " << analysisTree.genparticles_pdgid[igen] << endl;
			    TLorentzVector genPart; genPart.SetXYZT(analysisTree.genparticles_px[igen],
								    analysisTree.genparticles_py[igen],
								    analysisTree.genparticles_pz[igen],
								    analysisTree.genparticles_e[igen]);
			    if (analysisTree.genparticles_pdgid[igen]==23) {
			      genZ.SetXYZT(analysisTree.genparticles_px[igen],
					   analysisTree.genparticles_py[igen],
					   analysisTree.genparticles_pz[igen],
					   analysisTree.genparticles_e[igen]);
			    }
			    bool isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
			    bool isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
			    bool isLepton = isMuon || isElectron;
			    bool isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
			      fabs(analysisTree.genparticles_pdgid[igen])==14||
			      fabs(analysisTree.genparticles_pdgid[igen])==16;
			    bool isPrompt = analysisTree.genparticles_isPrompt[igen]||
			      analysisTree.genparticles_isPromptTauDecayProduct[igen];
			    
			    if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
			      if (isLepton&&
			     	  fabs(genPart.Eta())<2.4&&
			     	  genPart.Pt()>10) {
			     	genV += genPart;
			     	genL += genPart;
			      }
			      //			      if(isMuon) {
			      //				genV += genPart;
			      //				genL += genPart;
			      //			      }
			      
			      if (isNeutrino) 
				genV += genPart;
			    }
			  }
			  if (genV.Pt()<0.1) genV.SetXYZM(0.1,0.1,0.,0.);
			  
			  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
			    if (analysisTree.genparticles_pdgid[igen]==23 &&  analysisTree.genparticles_status[igen]==62)
			      TLorentzVector genZ; genZ.SetXYZT(analysisTree.genparticles_px[igen],
					   analysisTree.genparticles_py[igen],
					   analysisTree.genparticles_pz[igen],
					   analysisTree.genparticles_e[igen]);
			    float  visZPx= genZ.Px();
			    float visZPy= genZ.Py();
						
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
			      // TLorentzVector TwoMu = gen_mu1 +gen_mu2;
			      if ((gen_mu1.Pt()> 10 && gen_mu2.Pt()> 20)||(gen_mu1.Pt()> 20 && gen_mu2.Pt()> 10))
				massSelGen1H->Fill(genZ.M(),weight);
			    }
			  }
			  
			 
						  
			  
				
				if (applyPUreweighting_vertices) {
					int binNvert = vertexDataH->FindBin(analysisTree.primvertex_count);
					float_t dataNvert = vertexDataH->GetBinContent(binNvert);
	   				float_t mcNvert = vertexMcH->GetBinContent(binNvert);
					if (mcNvert < 1e-10){mcNvert=1e-10;}
					float_t vertWeight = dataNvert/mcNvert;
					weight *= vertWeight;
				}
				
				if (applyPUreweighting_official) {
				  //nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
					PUweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
					weight *= float(PUweight);
					//PUweightsOfficialH->Fill(PUweight);
				}
				
				if (applyTauTauSelection) {
					unsigned int nTaus = 0;
					if (analysisTree.gentau_count>0) {
						//	  cout << "Generated taus present" << endl;
						for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++itau) {
							// cout << itau << endl; "  : pt = " 
							//		 << analysisTree.gentau_visible_pt[itau] 
							//		 << "   eta = " <<  analysisTree.gentau_visible_eta[itau]
							//		 << "   mother = " << int(analysisTree.gentau_mother[itau]) << endl;
							if (int(analysisTree.gentau_mother[itau])==3) nTaus++;
						}
						// std::cout << "nTaus = " << nTaus << std::endl;//check
					}
					bool notTauTau = nTaus < 2;
					
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
					
					// Number muons from tau decays
					unsigned int nMuons = 0;
					for (unsigned int iMu=0; iMu<analysisTree.genparticles_count; ++iMu) {
						
						if (fabs(analysisTree.genparticles_pdgid[iMu])==13 && 
							analysisTree.genparticles_status[iMu]==1 && fabs(analysisTree.genparticles_mother[iMu]==5)) nMuons++;
						//selEventsAllMuons+=nMuons;
					} 
				}
				
				if (applyTopPtReweighting) {
				  float topPt = -1;
				  float antitopPt = -1;
				  for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
				    
				    if (analysisTree.genparticles_pdgid[igen]==6)
				      topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    
				    if (analysisTree.genparticles_pdgid[igen]==-6)
	      antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    
	    
				  }
				  if (topPt>0&&antitopPt>0) {
				    float topptweight = topPtWeight(topPt,antitopPt);
				    //	    cout << "toppt = " << topPt << "   antitoppt = " << antitopPt << "   weight = " << topptweight << endl;
				    weight *= topptweight;
				  }
				}
			}
			
			std:: vector<Period> periods;
			
			string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
			
			if (isData){
				
				std:: fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
				if (inputFileStream.fail()) {
					std:: cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
					std:: cout << "please check" << std::endl;
					std:: cout << "quitting program" << std::endl;
					exit(-1);
				}
				
				for(std::string s; std::getline(inputFileStream, s); ){
					periods.push_back(Period());
					std:: stringstream ss(s);
					ss >>  periods.back();
				}
				
				bool lumi = false;
				int n=analysisTree.event_run;
				int lum = analysisTree.event_luminosityblock;
				
				std::string num = std::to_string(n);
				std::string lnum = std::to_string(lum);
				for(const auto& a : periods)
				{
					
					if ( num.c_str() ==  a.name ) {
						//  std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
						//   std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
						for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
							
							//	cout<<b->lower<<"  "<<b->bigger<<endl;
							if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
						}
						auto last = std::prev(a.ranges.end());
						//std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
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
			
			unsigned int nBTagDiscriminant = 0;
			for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
				TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
				if (discr=BTagDiscriminator)
					nBTagDiscriminant = iBTag;
			}
			
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
			
			unsigned int nDoubleMuonHighPtFilter = 0;
			bool isDoubleMuonHighPtFilter = false;
			
			unsigned int nDoubleMuonLowPtFilter = 0;
			bool isDoubleMuonLowPtFilter = false;
			
			unsigned int nfilters = analysisTree.run_hltfilters->size();
			for (unsigned int i=0; i<nfilters; ++i) {
				TString HLTFilter(analysisTree.run_hltfilters->at(i));
				if (HLTFilter==MuonFilterName) {
					nMuonFilter = i;
					isMuonFilter = true;
				}
				if (HLTFilter==DoubleMuonHighPtFilterName) {
					nDoubleMuonHighPtFilter = i;
					isDoubleMuonHighPtFilter = true;
				}
				if (HLTFilter==DoubleMuonLowPtFilterName) {
					nDoubleMuonLowPtFilter = i;
					isDoubleMuonLowPtFilter = true;
				}
			}
			if (!isMuonFilter) {
				cout << "Filter " << MuonFilterName << " not found " << endl;
				exit(-1);
			}
			if (!isDoubleMuonHighPtFilter) {
				cout << "Filter " << DoubleMuonHighPtFilterName << " not found " << endl;
				exit(-1);
			}
			if (!isDoubleMuonLowPtFilter) {
				cout << "Filter " << DoubleMuonLowPtFilterName << " not found " << endl;
				exit(-1);
			}
			
			// muon selection
			vector<unsigned int> allMuons; allMuons.clear();
			vector<unsigned int> idMuons; idMuons.clear();
			vector<unsigned int> isoMuons; isoMuons.clear();
			vector<float> isoMuonsValue; isoMuonsValue.clear();
			vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
			for (unsigned int im = 0; im<analysisTree.muon_count; ++im){
				allMuons.push_back(im);
				//totalMuons++;
				//cout << "All muons number = "<< allMuons.size()<<endl; //check4
				bool muPassed = true;
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
			}
			
			unsigned int indx1 = 0;
			unsigned int indx2 = 0;
			bool isIsoMuonsPair = false;
			bool isMu1matched = false;
			bool isMu1HighPtmatched = false;
			bool isMu1LowPtmatched = false;
			bool isMu2matched = false;
			bool isMu2HighPtmatched = false;
			bool isMu2LowPtmatched = false;
			float isoMin = 9999;
			if (isoMuons.size()>0) {
				for (unsigned int im1=0; im1<isoMuons.size(); ++im1) {
					unsigned int index1 = isoMuons[im1];
					//bool isMu1matched = false;
					//bool isMu1HighPtmatched = false;
					//bool isMu1LowPtmatched = false;
					for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
						float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
						analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
						if (dRtrig>DRTrigMatch) continue;
						//  if (dRtrig<DRTrigMatch) totalMuonsTrigmatch++;//check
						if (analysisTree.trigobject_filters[iT][nMuonFilter] && 
							analysisTree.muon_pt[index1] > ptMuonHighCut &&
								fabs(analysisTree.muon_eta[index1]) < etaMuonHighCut) 
									isMu1matched = true;// nMu1match++;//check just a counter
						if (analysisTree.trigobject_filters[iT][nDoubleMuonHighPtFilter] &&
							analysisTree.muon_pt[index1] > ptMuonHighCut &&
								fabs(analysisTree.muon_eta[index1]) < etaMuonHighCut)
									isMu1HighPtmatched = true;
						if (analysisTree.trigobject_filters[iT][nDoubleMuonLowPtFilter] &&
							analysisTree.muon_pt[index1] > ptMuonLowCut &&
								fabs(analysisTree.muon_eta[index1]) < etaMuonLowCut)
									isMu1LowPtmatched = true;
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
							// if (dR>dRleptonsCut) totalMuonsDRmatch++;//check
							float ptProbe = TMath::Min(float(analysisTree.muon_pt[indexProbe]),float(ptBins[nPtBins]-0.1));
							float absEtaProbe = fabs(analysisTree.muon_eta[indexProbe]);
							int ptBin = binNumber(ptProbe,nPtBins,ptBins);
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
							// if (isMuonPassedIdIso[iMu]) ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight);
							// else ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
						}
					}
					//bool isMu2matched = false;
					//bool isMu2HighPtmatched = false;
					//bool isMu2LowPtmatched = false;
					for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
						unsigned int index2 = isoMuons[im2];
						float q1 = analysisTree.muon_charge[index1];
						float q2 = analysisTree.muon_charge[index2];
						//bool isMu2matched = false;
						//bool isMu2HighPtmatched = false;
						//bool isMu2LowPtmatched = false;
						for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
							float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
							analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
							if (dRtrig>DRTrigMatch) continue;
							if (analysisTree.trigobject_filters[iT][nMuonFilter] && 
								analysisTree.muon_pt[index2] > ptMuonHighCut && 
									fabs(analysisTree.muon_eta[index2]) < etaMuonHighCut) 
										isMu2matched = true;
							if (analysisTree.trigobject_filters[iT][nDoubleMuonHighPtFilter] && 
								analysisTree.muon_pt[index2] > ptMuonHighCut && 
									fabs(analysisTree.muon_eta[index2]) < etaMuonHighCut) 
										isMu2HighPtmatched = true;
							if (analysisTree.trigobject_filters[iT][nDoubleMuonLowPtFilter] && 
								analysisTree.muon_pt[index2] > ptMuonLowCut && 
									fabs(analysisTree.muon_eta[index2]) < etaMuonLowCut) 
										isMu2LowPtmatched = true;
						}
						bool isPairSelected = q1*q2 > 0;
						if (oppositeSign) isPairSelected = q1*q2 < 0;
						bool isTriggerMatch = (isMu1matched || isMu2matched);
						if (doubleMuonTrigger) 
							isTriggerMatch = (isMu1LowPtmatched && isMu2HighPtmatched) || (isMu2LowPtmatched && isMu1HighPtmatched);
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
				TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[indx1],
								analysisTree.muon_py[indx1],
								analysisTree.muon_pz[indx1],
								muonMass);
				
				TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[indx2],
								analysisTree.muon_py[indx2],
								analysisTree.muon_pz[indx2],
								muonMass);

				TLorentzVector genZ;// genZ.SetXYZM(0,0,0,91.2);
				TLorentzVector gen_mu1;
				TLorentzVector gen_mu2;
				if (!isData) {
					for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
						if (analysisTree.genparticles_pdgid[igen]==23 &&  analysisTree.genparticles_status[igen]==62)
							genZ.SetXYZT(analysisTree.genparticles_px[igen],
								     analysisTree.genparticles_py[igen],
								     analysisTree.genparticles_pz[igen],
								     analysisTree.genparticles_e[igen]);
						float visZPx= genZ.Px();
						float visZPy= genZ.Py();
						
						if (fabs(analysisTree.genparticles_pdgid[igen])==13 && 
						    analysisTree.genparticles_status[igen]==1) {
							
						  TLorentzVector genPart; genPart.SetPxPyPzE(analysisTree.genparticles_px[igen],
											     analysisTree.genparticles_py[igen],
											     analysisTree.genparticles_pz[igen],
											     analysisTree.genparticles_e[igen]);

						  float dRmu1 = deltaR(mu1.Eta(),mu1.Phi(),
								       genPart.Eta(),genPart.Phi());
						  if (dRmu1<0.5)
						    gen_mu1 = genPart;
							

						  float dRmu2 = deltaR(mu2.Eta(),mu2.Phi(),
								       genPart.Eta(),genPart.Phi()); 

						  if (dRmu2<0.5)
						    gen_mu2 = genPart;

						}
					}
				}
				

				n_genAccept = genV.M()>60 && genV.M()<120;
				float ptLeadGen = TMath::Max(gen_mu1.Pt(),gen_mu2.Pt());
				float ptTrailGen = TMath::Min(gen_mu1.Pt(),gen_mu2.Pt());
				n_genAccept = n_genAccept && ptLeadGen>20;
				n_genAccept = n_genAccept && ptTrailGen>10;
				n_genAccept = n_genAccept && fabs(gen_mu1.Eta())<2.4;
				n_genAccept = n_genAccept && fabs(gen_mu2.Eta())<2.4;



			// accessing Mva Met
				bool mvaMetFound = false;
				unsigned int metMuMu = 0; 
				for (unsigned int iMet=0; iMet<analysisTree.mvamet_count; ++iMet) {
				  if (analysisTree.mvamet_channel[iMet]==5) {
				    if (analysisTree.mvamet_lep1[iMet]==indx1&&
					analysisTree.mvamet_lep2[iMet]==indx2) {
				      metMuMu = iMet;
				      mvaMetFound = true;
				    }
				  }
				}
				float mvamet = 0;
				float mvamet_phi = 0;
				float mvamet_ex = 0;
				float mvamet_ey = 0;
				float n_covmet_xx =0;
				float n_covmet_xy =0;
				float n_covmet_yy =0;
				if (analysisTree.mvamet_count>0) {
				  mvamet_ex = analysisTree.mvamet_ex[metMuMu];
				  mvamet_ey = analysisTree.mvamet_ey[metMuMu];
				  float mvamet_ex2 = mvamet_ex * mvamet_ex;
				  float mvamet_ey2 = mvamet_ey * mvamet_ey;
				  n_covmet_xx = analysisTree.mvamet_sigxx[metMuMu];
				  n_covmet_xy = analysisTree.mvamet_sigxy[metMuMu];
				  n_covmet_yy = analysisTree.mvamet_sigyy[metMuMu];
				  
				   // std::cout << "xx = " << n_covmet_xx 
				   // 	    << "   xy = " << n_covmet_xy
				   // 	    << "   yy = " << n_covmet_yy << std::endl;

				  mvamet = TMath::Sqrt(mvamet_ex2+mvamet_ey2);
				  mvamet_phi = TMath::ATan2(mvamet_ey,mvamet_ex);
				}
    
				// selecting good jets ---
				
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
					
					// float energy = analysisTree.pfjet_e[jet];
					// float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
					// float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
					// float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
					// float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
					// float chm = analysisTree.pfjet_chargedmulti[jet];
					// float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
					// bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>2.4 || (elf<0.99 && chf>0 && chm>0));
					bool isPFJetId = looseJetiD(analysisTree,int(jet));
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
				
				
				TLorentzVector dimuon = mu1 + mu2;
				float visZPx=dimuon.Px();
				float visZPy=dimuon.Py();
				
				//	std::cout << "Before corrections : MetX = " << pfmet_ex << "   MetY = " << pfmet_ey << "   Met_Phi = " << pfmet_phi << std::endl;
				
				
				if (!isData) {
					if (applyLeptonSF) {
						//float pt1 = mu1.Pt(); 
						//if (pt1>1000) pt1 = 999;
						//float eta1 = TMath::Abs(mu1.Eta());
						//if (eta1>2.4) eta1 = 2.39;
						//int ptBin1 = binNumber(pt1,nPtBinsSF,ptBinsSF);
						//int etaBin1 = binNumber(eta1,nEtaBins,etaBins);
						//float sf1 = mueffSF[etaBin1][ptBin1]; // old sf method
						
						//Official Scale factor method
						double ptMu1 = (double)analysisTree.muon_pt[indx1];
						double ptMu2 = (double)analysisTree.muon_pt[indx2];
						double etaMu1 = (double)analysisTree.muon_eta[indx1];
						double etaMu2 = (double)analysisTree.muon_eta[indx2];
						
						//double IdIsoSF_mu1 = SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
						//double IdIsoSF_mu2 = SF_muonIdIso->get_ScaleFactor(ptMu2, etaMu2);
						
						isoweight_1=(float)SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
						isoweight_2=(float)SF_muonIdIso->get_ScaleFactor(ptMu2, etaMu2);
						
						MuSF_IdIso_Mu1H->Fill(isoweight_1);
						MuSF_IdIso_Mu2H->Fill(isoweight_2);
						
						weight = weight*isoweight_1*isoweight_2;
						
						double effDataTrig1 = SF_muonTrig->get_EfficiencyData(ptMu1, etaMu1);  
						double effDataTrig2 = SF_muonTrig->get_EfficiencyData(ptMu2, etaMu2);  

						double effMcTrig1 = SF_muonTrig->get_EfficiencyMC(ptMu1, etaMu1);  
						double effMcTrig2 = SF_muonTrig->get_EfficiencyMC(ptMu2, etaMu2);  

						double effTrigData = 1 - (1-effDataTrig1)*(1-effDataTrig2);
						double effMcTrig = 1 - (1-effMcTrig1)*(1-effMcTrig2);

						
						float Mu17EffData1 = (float)SF_muon17->get_EfficiencyData(ptMu1,etaMu1);
						float Mu17EffMC1   = (float)SF_muon17->get_EfficiencyMC(ptMu1,etaMu1);
						
						float Mu8EffData1 = (float)SF_muon8->get_EfficiencyData(ptMu1,etaMu1);
						float Mu8EffMC1   = (float)SF_muon8->get_EfficiencyMC(ptMu1,etaMu1);
						
						float Mu17EffData2 = (float)SF_muon17->get_EfficiencyData(ptMu2,etaMu2);
						float Mu17EffMC2   = (float)SF_muon17->get_EfficiencyMC(ptMu2,etaMu2);
						
						float Mu8EffData2 = (float)SF_muon8->get_EfficiencyData(ptMu2,etaMu2);
						float Mu8EffMC2   = (float)SF_muon8->get_EfficiencyMC(ptMu2,etaMu2);
						
						float trigWeightData = Mu17EffData1*Mu8EffData2 + Mu8EffData1*Mu17EffData2 - Mu17EffData1*Mu17EffData2;
						float trigWeightMC = Mu17EffMC1*Mu8EffMC2 + Mu8EffMC1*Mu17EffMC2 - Mu17EffMC1*Mu17EffMC2;
						
						if (doubleMuonTrigger){
							if (isMu1LowPtmatched && isMu2HighPtmatched) {
								trigweight_1 = (float)SF_muon17->get_ScaleFactor(ptMu1,etaMu1);
								trigweight_2 = (float)SF_muon8->get_ScaleFactor(ptMu2,etaMu2);
							}
							
							else if (isMu2LowPtmatched && isMu1HighPtmatched){
								trigweight_1 = (float)SF_muon17->get_ScaleFactor(ptMu1,etaMu1);
								trigweight_2 = (float)SF_muon8->get_ScaleFactor(ptMu2,etaMu2);
							}
							
							if (trigWeightMC>1e-6)
								trigweight = trigWeightData / trigWeightMC;
							
							effweight = isoweight_1*isoweight_2*trigweight;
							weight = weight*effweight;
						}
						
						else
						  if (effTrigData>0 && effMcTrig>0) {
						  double weightTrig = effTrigData/effMcTrig;
						  // std::cout << "mu 1 ->  pt = " << ptMu1 << "   eta = " << etaMu1 << std::endl;
						  // std::cout << "mu 2 ->  pt = " << ptMu2 << "   eta = " << etaMu2 << std::endl;
						  // std::cout << "WeightTrig = " << weightTrig << std::endl;
						  weight = weight*weightTrig;
						}
						
						
						
						//	    std::cout << "mu1 : pt=" << mu1.Pt() << " (" << ptBin1 << ")"
						//		      << "  eta=" << mu1.Eta() << " (" << etaBin1 << ")"
						//		      << "  eff(data)=" << mueffData[etaBin1][ptBin1] 
						//		      << "  eff(MC)=" << mueffMC[etaBin1][ptBin1]
						//		      << "  SF=" << mueffSF[etaBin1][ptBin1] << std::endl;
						
						// float pt2 = mu2.Pt(); 
						//if (pt2>1000) pt2 = 999;
						//float eta2 = TMath::Abs(mu2.Eta());
						//if (eta2>2.4) eta2 = 2.39;
						//int ptBin2 = binNumber(pt2,nPtBinsSF,ptBinsSF);
						//int etaBin2 = binNumber(eta2,nEtaBins,etaBins);
						//float sf2 = mueffSF[etaBin2][ptBin2];
						
 						//	    std::cout << "mu2 : pt=" << mu2.Pt() << " (" << ptBin2 << ")"
						//		      << "  eta=" << mu2.Eta() << " (" << etaBin2 << ")"
						//		      << "  eff(data)=" << mueffData[etaBin2][ptBin2] 
						//		      << "  eff(MC)=" << mueffMC[etaBin2][ptBin2]
						//		      << "  SF=" << mueffSF[etaBin2][ptBin2] << std::endl;
						
						// float sf = sf1 * sf2;
						//weight *= sf;
						// std::cout << std::endl;
						
					}
					
					if (applyMEtRecoilCorrections) {
					  float pfmetcorr_ex = pfmet_ex;
					  // std::cout << "V : px = " << genV.Px()
					  // 	    << "    py = " << genV.Py()
					  // 	    << "    pz = " << genV.Pz() 
					  // 	    << "    mass = " << genV.M() << std::endl;
					  // std::cout << "Z : px = " << genZ.Px()
					  // 	    << "    py = " << genZ.Py()
					  // 	    << "    pz = " << genZ.Pz() 
					  // 	    << "    mass = " << genZ.M() << std::endl;
					  // std::cout << "L : px = " << genL.Px()
					  // 	    << "    py = " << genL.Py()
					  // 	    << "    pz = " << genL.Pz() 
					  // 	    << "    mass = " << genL.M() << std::endl;
					  // std::cout << std::endl;
					  float pfmetcorr_ey = pfmet_ey;
					  recoilPFMetCorrector.CorrectByMeanResolution(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),nJets30,pfmetcorr_ex,pfmetcorr_ey);
					  //					  std::cout << "PFMet : (" << pfmet_ex << "," << pfmet_ey << ")  "
					  //						    << "  (" << pfmetcorr_ex << "," << pfmetcorr_ey << ")" << std::endl; 
					  pfmet_phi = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
					  pfmet = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
					  pfmet_ex = pfmetcorr_ex;
					  pfmet_ey = pfmetcorr_ey;

	  
	  
					  //					  float puppimetcorr_ex = puppimet_ex;
					  //					  float puppimetcorr_ey = puppimet_ey;
					  //					  recoilPuppiMetCorrector.Correct(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,puppimetcorr_ex,puppimetcorr_ey);
					  //	  std::cout << "PuppiMet : (" << puppimet_ex << "," << puppimet_ey << ")  "
					  //		    << "  (" << puppimetcorr_ex << "," << puppimetcorr_ey << ")" << std::endl; 
					  //					  puppimet_phi = TMath::ATan2(puppimetcorr_ey,puppimetcorr_ex);
					  //					  puppimet = TMath::Sqrt(puppimetcorr_ex*puppimetcorr_ex+puppimetcorr_ey*puppimetcorr_ey);
					  
					  //					  puppimet_ex = puppimetcorr_ex;
					  //					  puppimet_ey = puppimetcorr_ey;
					  
					  float mvametcorr_ex = mvamet_ex;
					  float mvametcorr_ey = mvamet_ey;
					  recoilMvaMetCorrector.CorrectByMeanResolution(mvamet_ex,mvamet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),nJets30,mvametcorr_ex,mvametcorr_ey);
					  // 	  std::cout << "MvaMet : (" << mvamet_ex << "," << mvamet_ey << ")  "
					  //	                      << "  (" << mvametcorr_ex << "," << mvametcorr_ey << ")" << std::endl;
					  mvamet_phi = TMath::ATan2(mvametcorr_ey,mvametcorr_ex);
					  mvamet = TMath::Sqrt(mvametcorr_ex*mvametcorr_ex+mvametcorr_ey*mvametcorr_ey);
					  mvamet_ex = mvametcorr_ex;
					  mvamet_ey = mvametcorr_ey; 

					}
					//	std::cout << "After  corrections : MetX = " << pfmet_ex << "   MetY = " << pfmet_ey << "   Met_Phi = " << pfmet_phi << std::endl;
				}
				
				
				
				float massSel = dimuon.M();
				
				//bisector of dimuon transerve momenta for defining the dzeta variables
				float mu1UnitX = mu1.Px()/mu1.Pt();
				float mu1UnitY = mu1.Py()/mu1.Pt();

				float mu2UnitX = mu2.Px()/mu2.Pt();
				float mu2UnitY = mu2.Py()/mu2.Pt();
				
				float zetaX = mu1UnitX + mu2UnitX;
				float zetaY = mu1UnitY + mu2UnitY;
		  
				float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
				zetaX = zetaX/normZeta;
				zetaY = zetaY/normZeta;
				
				//float vectorX = pfmet_ex + mu2.Px() + mu1.Px();
				//float vectorY = pfmet_ey + mu2.Py() + mu1.Py();
				float vectorX = mvamet_ex + mu2.Px() + mu1.Px();
				float vectorY = mvamet_ey + mu2.Py() + mu1.Py();
				
				float vectorVisX = mu2.Px() + mu1.Px();
				float vectorVisY = mu2.Py() + mu1.Py();
				
				// computation of DZeta variable
				float PZeta = vectorX*zetaX + vectorY*zetaY;
				float PVisZeta = vectorVisX*zetaX + vectorVisY*zetaY;
				float DZeta = PZeta - 1.85*PVisZeta;
				
				if (massSel>20) {
					massSelH->Fill(massSel,weight);
					/*for (int iScale=0; iScale<21; ++ iScale) {
						float scaleFactor = 0.98 + 0.002*float(iScale);
						massSelScaleH[iScale]->Fill(massSel*scaleFactor,weight);
						}*/
					massSelGen1H->Fill(genZ.M(),weight);
					
					 n_leadingPt = analysisTree.muon_pt[indx1];
					 n_trailingPt = analysisTree.muon_pt[indx2];
					 n_leadingEta = analysisTree.muon_eta[indx1];
					 n_trailingEta = analysisTree.muon_eta[indx2];
					 n_leadingPhi = analysisTree.muon_phi[indx1];
					 n_trailingPhi = analysisTree.muon_phi[indx2];
					 if (analysisTree.muon_pt[indx1]<analysisTree.muon_pt[indx2]) {
					   n_leadingPt = analysisTree.muon_pt[indx2];
					   n_trailingPt = analysisTree.muon_pt[indx1];
					   n_leadingEta = analysisTree.muon_eta[indx2];
					   n_trailingEta = analysisTree.muon_eta[indx1];
					   n_leadingPhi = analysisTree.muon_phi[indx2];
					   n_trailingPhi = analysisTree.muon_phi[indx1];
					 }
					
					ptLeadingMuSelH->Fill(n_leadingPt,weight);
					ptTrailingMuSelH->Fill(n_trailingPt,weight);
					etaLeadingMuSelH->Fill(n_leadingEta,weight);
					etaTrailingMuSelH->Fill(n_trailingEta,weight);
					ptScatter->Fill(n_leadingPt,n_trailingPt,weight);
					hprof2D_pt->Fill(n_leadingPt,n_trailingPt,weight);
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
					float ptRatio_test = 0.0; 
					if (sumMuonPt != 0){
						ptRatio = (dimuonPt/sumMuonPt);
						ptRatio_test = (dimuonPt/sumMuonPt);
					}
					float dcaSigdxy_muon1 = 0.0;
					float dcaSigdxy_muon2 = 0.0;
					float dcaSigdz_muon1  = 0.0;
					float dcaSigdz_muon2  =0.0;
					float dcaSig2Mu2D = 0.0;
					float dcaSig2Mu3D = 0.0;
					float sig2Mu2D = 0.0;
					float sig2Mu3D = 0.0;
					float phi_LeadingMu_MET= 0.0;
					float phi_TrailingMu_MET =0.0;
					float phi_PosMu_MET =0.0;
					float phi_TwoMu =0.0;
					float leadingPt = 0.0;
					float trailingPt = 0.0;
					float leadingEta = 0.0;
					float trailing = 0.0;
					float jets = 0.0;
					float noOfvertices = 0.0;
					n_m_sv = -9999;
					n_pt_sv = -9999;
					n_eta_sv = -9999;
					n_phi_sv = -9999;
					 
					 
					if (analysisTree.muon_dxyerr[indx1] != 0){
						//dcaSigdxy_mu1=(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]);
						dcaSigdxy_muon1= log10(fabs(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]));
						//std::cout << "dcaSigdxy_muon1 is "<< dcaSigdxy_muon1 <<" and before log is "<< dcaSigdxy_mu1<< " nan is "  << isnan(dcaSigdxy_muon1)<< endl;
					}
					
					if (analysisTree.muon_dxyerr[indx2] != 0){
						//dcaSigdxy_mu2=(analysisTree.muon_dxy[indx2]/analysisTree.muon_dxyerr[indx2]);
						dcaSigdxy_muon2= log10(fabs(analysisTree.muon_dxy[indx2]/analysisTree.muon_dxyerr[indx2]));
					}
					
					if (analysisTree.muon_dzerr[indx1] != 0){
						//dcaSigdz_mu1 =(analysisTree.muon_dz[indx1]/analysisTree.muon_dzerr[indx1]);
						dcaSigdz_muon1 = log10(fabs(analysisTree.muon_dz[indx1]/analysisTree.muon_dzerr[indx1]));
					}
			  
					if (analysisTree.muon_dzerr[indx2] != 0){
						//dcaSigdz_mu2 = (analysisTree.muon_dz[indx2]/analysisTree.muon_dzerr[indx2]);
						dcaSigdz_muon2 = log10(fabs(analysisTree.muon_dz[indx2]/analysisTree.muon_dzerr[indx2]));
					}
					for(unsigned int dimu=0; dimu<analysisTree.dimuon_count; ++dimu){
					  if (analysisTree.dimuon_dist2DE != 0){
					    sig2Mu2D = (analysisTree.dimuon_dist2D[dimu]/analysisTree.dimuon_dist2DE[dimu]);
					    dcaSig2Mu2D = log10(fabs(analysisTree.dimuon_dist2D[dimu]/analysisTree.dimuon_dist2DE[dimu]));
					  }
						
					  if (analysisTree.dimuon_dist3DE != 0){
					    sig2Mu3D = (analysisTree.dimuon_dist3D[dimu]/analysisTree.dimuon_dist3DE[dimu]);
					    dcaSig2Mu3D = log10(fabs(analysisTree.dimuon_dist3D[dimu]/analysisTree.dimuon_dist3DE[dimu]));
					  }
					}	
					
					//filling the histograms for discriminators
					h_dimuonEta->Fill(dimuonEta,weight);
					h_dimuonPt->Fill(dimuonPt,weight);
					// if (genmatch_m1 && genmatch_m2) h_dimuonEta_genMuMatch->Fill(dimuonEta,weight);
					h_ptRatio->Fill(ptRatio,weight);
					h_dxy_muon1->Fill(analysisTree.muon_dxy[indx1],weight);
					h_dxy_muon2->Fill(analysisTree.muon_dxy[indx2],weight);
					h_dz_muon1->Fill(analysisTree.muon_dz[indx1],weight);
					h_dz_muon2->Fill(analysisTree.muon_dz[indx2],weight);
					
					//int iEta=0, iPt=0;
					//if (dimuonEta < 0.9) iEta = 0;
					//if (dimuonEta>0.9 && dimuonEta<1.2) iEta = 1;
					//if (dimuonEta>1.2 && dimuonEta<2.1) iEta = 2;
					//if (dimuonEta>2.1 && dimuonEta<2.4) iEta = 3;
					////if (dimuonEta>2.4) iEta = 4;
					//if (dimuonPt < 10) iPt = 0;
					//if (dimuonPt>10 && dimuonPt<15) iPt = 1;
					//if (dimuonPt>15 && dimuonPt<20) iPt = 2;
					//if (dimuonPt>20 && dimuonPt<25) iPt = 3;
					//if (dimuonPt>25 && dimuonPt<30) iPt = 4;
					//if (dimuonPt>30 && dimuonPt<40) iPt = 5;
					//if (dimuonPt>40 && dimuonPt<50) iPt = 6;
					//if (dimuonPt>50 && dimuonPt<60) iPt = 7;
					//if (dimuonPt>60) iPt = 8;
					
					//dxyMu1[iEta][iPt]->Fill(analysisTree.muon_dxy[indx1],weight);
					//dxyMu2[iEta][iPt]->Fill(analysisTree.muon_dxy[indx2],weight);
					//dzMu1[iEta][iPt]->Fill(analysisTree.muon_dz[indx1],weight);
					//dzMu2[iEta][iPt]->Fill(analysisTree.muon_dz[indx2],weight);
					
					h_dcaSigdxy_muon1->Fill(dcaSigdxy_muon1,weight);
					h_dcaSigdxy_muon2->Fill(dcaSigdxy_muon2,weight);
					h_dcaSigdz_muon1->Fill(dcaSigdz_muon1,weight);
					h_dcaSigdz_muon2->Fill(dcaSigdz_muon2,weight);
					
					h_dcaSig2Mu2D->Fill(dcaSig2Mu2D,weight);
					if (nJets30==0) h_dcaSig2Mu2D_Jet0->Fill(dcaSig2Mu2D,weight);
					h_dcaSig2Mu3D->Fill(dcaSig2Mu3D,weight);
					h_sig2Mu2D->Fill(sig2Mu2D,weight);
					h_sig2Mu3D->Fill(sig2Mu3D,weight);
					
					float q1 = analysisTree.muon_charge[indx1];
					float q2 = analysisTree.muon_charge[indx2];
					
					//	phi_LeadingMu_MET = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],pfmet_ex,pfmet_ey);
					phi_LeadingMu_MET = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],mvamet_ex,mvamet_ey);
					h_phi_leadingMu_MET->Fill(fabs(phi_LeadingMu_MET),weight);
					
					//phi_TrailingMu_MET = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],pfmet_ex,pfmet_ey);
					phi_TrailingMu_MET = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],mvamet_ex,mvamet_ey);
					h_phi_trailingMu_MET->Fill(fabs(phi_TrailingMu_MET),weight);
					
					if (q1>0)
						phi_PosMu_MET = phi_LeadingMu_MET;
					else
						phi_PosMu_MET = phi_TrailingMu_MET;
					
					h_phi_PosMu_MET->Fill(fabs(phi_PosMu_MET),weight);
					
					phi_TwoMu = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],analysisTree.muon_px[indx2],analysisTree.muon_py[indx2]);
					h_phi_TwoMu->Fill(phi_TwoMu,weight);
					if (fabs(phi_TwoMu)>2) 
						h_ptRatio_test-> Fill (ptRatio_test, weight);
					
					h_DZeta->Fill(DZeta,weight);
					h_Njets->Fill(analysisTree.pfjet_count,weight);
					
					n_phiangle=fabs(phi_PosMu_MET);
					n_dimuonEta=dimuonEta;
					n_mvamet= mvamet;
					n_DZeta = DZeta;
				
					n_mvamet_ex = mvamet_ex;
					n_mvamet_ey = mvamet_ey;
					
					// int iJetBin=0;
					// if (nJets30==1)
					//   iJetBin = 1;
					// else if (nJets30>1)
					//   iJetBin = 2;
					// ptRatio_nJetsH[iJetBin]->Fill(ptRatio,weight);
					
					if (dcaSig2Mu2D < 0.5) dimuonMass_dca-> Fill(massSel,weight);
					NumberOfVerticesH->Fill(float(analysisTree.primvertex_count),weight);
					float metSel = sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
					metSelH->Fill(pfmet,weight);
					mvametSelH->Fill(mvamet,weight);
					/*for (int iScale=0; iScale<21; ++ iScale) {
						float scaleFactor = 0.7 + 0.03*float(iScale);
						metSelScaleH[iScale]->Fill(pfmet*scaleFactor,weight);
						}*/
					
					 float mvaValue_ =  reader->EvaluateMVA("BDT");
					 //std::cout<< "mvaBDT   "<<mvaValue_<<std::endl;
					 histBdt->Fill(mvaValue_, weight);
					 if(massSel > 20 && massSel < 70)histMva_massCut->Fill(mvaValue_, weight);
					 // for (int i=0;i<iCut;i++){
					 //   double cut = i*0.05;
					 //   if(mvaValue_ > cut) InvMass[i]->Fill(n_dimuonMass, weight);
					 // }
					 if (mvaValue_ >0.5) dimuonMass_BDT0p5->Fill(massSel,weight);
					 
					

					 /////*****svFit*******
					 if (mvaValue_ > 0.5 && massSel < 80){
					   //define MET
					   double measuredMETx = n_mvamet_ex;
					   double measuredMETy = n_mvamet_ey;
					   
					   // define MET covariance
					   TMatrixD covMET(2, 2);
					 
					   // std::cout << "covmetxx " << n_covmet_xx << "\t covmetxy " << n_covmet_xy  << "\t covmetyy " << n_covmet_yy <<std::endl;
					   
					   covMET[0][0] = n_covmet_xx;
					   covMET[1][0] = n_covmet_xy;
					   covMET[0][1] = n_covmet_xy;
					   covMET[1][1] = n_covmet_yy;
					 
					   // define lepton four vectors
					   std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
					   
					   measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, n_leadingPt, n_leadingEta, n_leadingPhi, 105.658e-3)); 
					   measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, n_trailingPt ,n_trailingEta, n_trailingPhi, 105.658e-3));
					 
					   SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
					   algo.addLogM(false);  
					  
					   
					   algo.shiftVisPt(true, inputFile_visPtResolution);
					   
					   algo.integrateMarkovChain();
					   n_m_sv = algo.getMass(); // return value is in units of GeV
					   
					   mass_sv->Fill(n_m_sv,weight);
					   
					   bool algoVerify = false;
					   algoVerify = algo.isValidSolution();
					   if ( algoVerify ) {
					     //	  std::cout << "SVfit mass (pfmet)    = " << m_sv << std::endl;
					   } else {
					     std::cout << "sorry -- status of NLL is not valid [" << algoVerify << "]" << std::endl;
					   }
					   
					   n_pt_sv = algo.pt(); 
					   n_eta_sv = algo.eta();
					   n_phi_sv = algo.phi();
					   
					   
					 }
				       //fill ntuples for BDT training
					 //n_dimuonEta=dimuonEta;
					n_ptRatio=ptRatio;
					n_dxy_muon1= analysisTree.muon_dxy[indx1];
					n_dxy_muon2= analysisTree.muon_dxy[indx2]; 
					n_dz_muon1= analysisTree.muon_dz[indx1]; 
					n_dz_muon2= analysisTree.muon_dz[indx2]; 
					n_dcaSigdxy1=dcaSigdxy_muon1;
					n_dcaSigdxy2=dcaSigdxy_muon2;
					//n_dcaSigdxy_mu1=dcaSigdxy_mu1;
					//n_dcaSigdxy_mu2=dcaSigdxy_mu2;
					n_dcaSigdz1=dcaSigdz_muon1;
					n_dcaSigdz2=dcaSigdz_muon2;
					//n_dcaSigdz_mu1=dcaSigdz_mu1;
					//n_dcaSigdz_mu2=dcaSigdz_mu2;
					n_dcaSig2Mu2D = dcaSig2Mu2D;
					n_dcaSig2Mu3D = dcaSig2Mu3D;
					n_sig2Mu2D = sig2Mu2D;
					n_sig2Mu3D = sig2Mu3D;
					n_MissingEt=metSel;
					//n_phiangle=fabs(phi_PosMu_MET);
					n_twomuPhi=fabs(phi_TwoMu);
					//n_DZeta = DZeta;
					n_genWeight = weight;
					n_dimuonMass = massSel;
					n_met= pfmet;
					//n_mvamet= mvamet;
					// n_leadingPt = analysisTree.muon_pt[indx1];
					// n_trailingPt = analysisTree.muon_pt[indx2];
					// n_leadingEta = analysisTree.muon_eta[indx1];
					// n_trailingEta = analysisTree.muon_eta[indx2];
					// n_leadingPhi = analysisTree.muon_phi[indx1];
					// n_trailingPhi = analysisTree.muon_phi[indx2];
					// if (analysisTree.muon_pt[indx1]<analysisTree.muon_pt[indx2]) {
					//   n_leadingPt = analysisTree.muon_pt[indx2];
					//   n_trailingPt = analysisTree.muon_pt[indx1];
					//   n_leadingEta = analysisTree.muon_eta[indx2];
					//   n_trailingEta = analysisTree.muon_eta[indx1];
					//   n_leadingPhi = analysisTree.muon_phi[indx2];
					//   n_trailingPhi = analysisTree.muon_phi[indx1];
					// }
					n_jets = double(nJets30);
					n_noOfvertices = analysisTree.primvertex_count;
					n_genZ = genZ.M();
					n_mvaBDT = mvaValue_;
					n_mvamet_ex = mvamet_ex;
					n_mvamet_ey = mvamet_ey;
					
					 if (fillBDTNTuple) {
					  
					  
					   T->Fill();
					 }
					
					// if (massSel>70&&massSel<110) {
					// 	float unitX = dimuon.Px()/dimuon.Pt();
					// 	float unitY = dimuon.Py()/dimuon.Pt();
					// 	float phiUnit = TMath::ATan2(unitY,unitX);
					// 	float recoilParal = pfmet_ex*unitX + pfmet_ey*unitY;
					// 	float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
					// 	float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
					// 	float recoilPerp = pfmet_ex*perpUnitX + pfmet_ey*perpUnitY;
					// 	int jetBin = 0;
					// 	if (nJets30==1)
					// 		jetBin = 1;
					// 	else if (nJets30>1)
					// 		jetBin = 2;
					// 	recoilZParalH[jetBin]->Fill(recoilParal,weight);
					// 	recoilZPerpH[jetBin]->Fill(recoilPerp,weight);
						
					// 	int iJets=0, iPtBins=0;
					// 	if (nJets30==1) iJets = 1;
					// 	if (nJets30>=2) iJets = 2;
					// 	if (dimuonPt < 10) iPtBins = 0;
					// 	if (dimuonPt>10 && dimuonPt<20) iPtBins = 1;
					// 	if (dimuonPt>20 && dimuonPt<30) iPtBins = 2;
					// 	if (dimuonPt>30 && dimuonPt<50) iPtBins = 3;
					// 	if (dimuonPt>50) iPtBins = 4;
					// 	recoilZParal_Ptbins_nJetsH[iJets][iPtBins]->Fill(recoilParal,weight);
					// 	recoilZPerp_Ptbins_nJetsH[iJets][iPtBins]->Fill(recoilPerp,weight);
						
					// 	int iJetBin=0;
					// 	if (nJets30==1)
					// 		iJetBin = 1;
					// 	else if (nJets30>1)
					// 		iJetBin = 2;
					// 	ptRatio_nJetsH[iJetBin]->Fill(ptRatio,weight);
					// }
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
	delete inputFile_visPtResolution;
}



   
