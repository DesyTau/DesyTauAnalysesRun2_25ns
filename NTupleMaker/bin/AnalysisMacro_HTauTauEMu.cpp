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

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");
  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float MEtCut         = cfg.get<float>("MEtCut");
  const float dZetaCut       = cfg.get<float>("dZetaCut");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");

  const float MEtCutTTJets   = cfg.get<float>("MEtCutTTJets");
  // cuts on jets
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagLooseCut = cfg.get<float>("btagLooseCut");
  const float btagMediumCut = cfg.get<float>("btagMediumCut");

  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");


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
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

  TH1F * muonPtAllH = new TH1F("muonPtAllH","",40,0,200);
  TH1F * electronPtAllH = new TH1F("electronPtAllH","",40,0,200);

  // histograms (dilepton selection)
  TH1F * electronPtH  = new TH1F("electronPtH","",40,0,200);
  TH1F * electronEtaH = new TH1F("electronEtaH","",50,-2.5,2.5); 
  TH1F * muonPtH  = new TH1F("muonPtH","",40,0,200);
  TH1F * muonEtaH = new TH1F("muonEtaH","",50,-2.5,2.5); 
  TH1F * dileptonMassH = new TH1F("dileptonMassH","",40,0,200);
  TH1F * dileptonPtH = new TH1F("dileptonPtH","",40,0,200);
  TH1F * dileptonEtaH = new TH1F("dileptonEtaH","",100,-5,5);
  TH1F * dileptondRH = new TH1F("dileptondRH","",60,0,6);
  TH1F * ETmissH = new TH1F("ETmissH","",40,0,200);
  TH1F * MtH = new TH1F("MtH","",40,0,200);
  TH1F * DZetaH = new TH1F("DZetaH","",60,-400,200);
  TH1F * nJets30H = new TH1F("nJets30H","",10,-0.5,9.5);
  TH1F * nJets20H = new TH1F("nJets20H","",10,-0.5,9.5);
  TH1F * nBLooseJetsH = new TH1F("nBLooseJetsH","",10,-0.5,9.5);
  TH1F * nBMediumJetsH = new TH1F("nBMediumJetsH","",10,-0.5,9.5);

  // histograms (dilepton selection + DZeta)
  TH1F * electronPtDZetaCutH  = new TH1F("electronPtDZetaCutH","",40,0,200);
  TH1F * electronEtaDZetaCutH = new TH1F("electronEtaDZetaCutH","",50,-2.5,2.5); 
  TH1F * muonPtDZetaCutH  = new TH1F("muonPtDZetaCutH","",40,0,200);
  TH1F * muonEtaDZetaCutH = new TH1F("muonEtaDZetaCutH","",50,-2.5,2.5); 
  TH1F * dileptonMassDZetaCutH = new TH1F("dileptonMassDZetaCutH","",40,0,200);
  TH1F * dileptonPtDZetaCutH = new TH1F("dileptonPtDZetaCutH","",40,0,200);
  TH1F * dileptonEtaDZetaCutH = new TH1F("dileptonEtaDZetaCutH","",100,-5,5);
  TH1F * dileptondRDZetaCutH = new TH1F("dileptondRDZetaCutH","",60,0,6);
  TH1F * ETmissDZetaCutH = new TH1F("ETmissDZetaCutH","",40,0,200);
  TH1F * MtDZetaCutH = new TH1F("MtDZetaCutH","",40,0,200);
  TH1F * DZetaDZetaCutH = new TH1F("DZetaDZetaCutH","",60,-400,200);
  TH1F * nJets30DZetaCutH = new TH1F("nJets30DZetaCutH","",10,-0.5,9.5);
  TH1F * nJets20DZetaCutH = new TH1F("nJets20DZetaCutH","",10,-0.5,9.5);
  TH1F * nBLooseJetsDZetaCutH = new TH1F("nBLooseJetsDZetaCutH","",10,-0.5,9.5);
  TH1F * nBMediumJetsDZetaCutH = new TH1F("nBMediumJetsDZetaCutH","",10,-0.5,9.5);


  // histograms (dilepton selection + MEtCut)
  TH1F * electronPtMEtCutH  = new TH1F("electronPtMEtCutH","",40,0,200);
  TH1F * electronEtaMEtCutH = new TH1F("electronEtaMEtCutH","",50,-2.5,2.5); 
  TH1F * muonPtMEtCutH  = new TH1F("muonPtMEtCutH","",40,0,200);
  TH1F * muonEtaMEtCutH = new TH1F("muonEtaMEtCutH","",50,-2.5,2.5); 
  TH1F * dileptonMassMEtCutH = new TH1F("dileptonMassMEtCutH","",40,0,200);
  TH1F * dileptonPtMEtCutH = new TH1F("dileptonPtMEtCutH","",40,0,200);
  TH1F * dileptonEtaMEtCutH = new TH1F("dileptonEtaMEtCutH","",100,-5,5);
  TH1F * dileptondRMEtCutH = new TH1F("dileptondRMEtCutH","",60,0,6);
  TH1F * ETmissMEtCutH = new TH1F("ETmissMEtCutH","",40,0,200);
  TH1F * MtMEtCutH = new TH1F("MtMEtCutH","",40,0,200);
  TH1F * DZetaMEtCutH = new TH1F("DZetaMEtCutH","",60,-400,200);
  TH1F * nJets30MEtCutH = new TH1F("nJets30MEtCutH","",10,-0.5,9.5);
  TH1F * nJets20MEtCutH = new TH1F("nJets20MEtCutH","",10,-0.5,9.5);
  TH1F * nBLooseJetsMEtCutH = new TH1F("nBLooseJetsMEtCutH","",10,-0.5,9.5);
  TH1F * nBMediumJetsMEtCutH = new TH1F("nBMediumJetsMEtCutH","",10,-0.5,9.5);

  // histograms (dilepton selection + DZeta + MEtCut)
  TH1F * electronPtSelH  = new TH1F("electronPtSelH","",40,0,200);
  TH1F * electronEtaSelH = new TH1F("electronEtaSelH","",50,-2.5,2.5); 
  TH1F * muonPtSelH  = new TH1F("muonPtSelH","",40,0,200);
  TH1F * muonEtaSelH = new TH1F("muonEtaSelH","",50,-2.5,2.5); 
  TH1F * dileptonMassSelH = new TH1F("dileptonMassSelH","",40,0,200);
  TH1F * dileptonPtSelH = new TH1F("dileptonPtSelH","",40,0,200);
  TH1F * dileptonEtaSelH = new TH1F("dileptonEtaSelH","",100,-5,5);
  TH1F * dileptondRSelH = new TH1F("dileptondRSelH","",60,0,6);
  TH1F * ETmissSelH = new TH1F("ETmissSelH","",40,0,200);
  TH1F * MtSelH = new TH1F("MtSelH","",40,0,200);
  TH1F * DZetaSelH = new TH1F("DZetaSelH","",60,-400,200);
  TH1F * nJets30SelH = new TH1F("nJets30SelH","",10,-0.5,9.5);
  TH1F * nJets20SelH = new TH1F("nJets20SelH","",10,-0.5,9.5);
  TH1F * nBLooseJetsSelH = new TH1F("nBLooseJetsSelH","",10,-0.5,9.5);
  TH1F * nBMediumJetsSelH = new TH1F("nBMediumJetsSelH","",10,-0.5,9.5);

  // histograms (dilepton selection + TTJets selection)
  TH1F * electronPtTTJetsSelH  = new TH1F("electronPtTTJetsSelH","",40,0,200);
  TH1F * electronEtaTTJetsSelH = new TH1F("electronEtaTTJetsSelH","",50,-2.5,2.5); 
  TH1F * muonPtTTJetsSelH  = new TH1F("muonPtTTJetsSelH","",40,0,200);
  TH1F * muonEtaTTJetsSelH = new TH1F("muonEtaTTJetsSelH","",50,-2.5,2.5); 
  TH1F * dileptonMassTTJetsSelH = new TH1F("dileptonMassTTJetsSelH","",40,0,200);
  TH1F * dileptonPtTTJetsSelH = new TH1F("dileptonPtTTJetsSelH","",40,0,200);
  TH1F * dileptonEtaTTJetsSelH = new TH1F("dileptonEtaTTJetsSelH","",100,-5,5);
  TH1F * dileptondRTTJetsSelH = new TH1F("dileptondRTTJetsSelH","",60,0,6);
  TH1F * ETmissTTJetsSelH = new TH1F("ETmissTTJetsSelH","",40,0,200);
  TH1F * MtTTJetsSelH = new TH1F("MtTTJetsSelH","",40,0,200);
  TH1F * DZetaTTJetsSelH = new TH1F("DZetaTTJetsSelH","",60,-400,200);
  TH1F * nJets30TTJetsSelH = new TH1F("nJets30TTJetsSelH","",10,-0.5,9.5);
  TH1F * nJets20TTJetsSelH = new TH1F("nJets20TTJetsSelH","",10,-0.5,9.5);
  TH1F * nBLooseJetsTTJetsSelH = new TH1F("nBLooseJetsTTJetsSelH","",10,-0.5,9.5);
  TH1F * nBMediumJetsTTJetsSelH = new TH1F("nBMediumJetsTTJetsSelH","",10,-0.5,9.5);

  TH1F * TTJetsDiLeptonH = new TH1F("TTJetsDiLeptonH","",1,-0.5,0.5);
  TH1F * TTJetsDiLeptonElePtH = new TH1F("TTJetsDiLeptonElePtH","",1000,0,1000);
  TH1F * TTJetsDiLeptonMuPtH = new TH1F("TTJetsDiLeptonMuPtH","",1000,0,1000);
  TH1F * TTJetsDiLeptonEleEtaH = new TH1F("TTJetsDiLeptonElePtH","",100,-5,5);
  TH1F * TTJetsDiLeptonMuEtaH = new TH1F("TTJetsDiLeptonMuPtH","",100,-5,5);
  TH1F * TTJetsDiLeptonAccH = new TH1F("TTJetsDiLeptonAccH","",1,-0.5,0.5);

  TH1F * DYJetsDiLeptonH = new TH1F("DYJetsDiLeptonH","",1,-0.5,0.5);
  TH1F * DYJetsDiLeptonElePtH = new TH1F("DYJetsDiLeptonElePtH","",1000,0,1000);
  TH1F * DYJetsDiLeptonMuPtH = new TH1F("DYJetsDiLeptonMuPtH","",1000,0,1000);
  TH1F * DYJetsDiLeptonEleEtaH = new TH1F("DYJetsDiLeptonElePtH","",100,-5,5);
  TH1F * DYJetsDiLeptonMuEtaH = new TH1F("DYJetsDiLeptonMuPtH","",100,-5,5);
  TH1F * DYJetsDiLeptonAccH = new TH1F("DYJetsDiLeptonAccH","",1,-0.5,0.5);


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

      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      vector<int> indexW; indexW.clear();
      vector<int> indexNu; indexNu.clear(); 
      vector<int> indexMu; indexMu.clear();
      vector<int> indexE; indexE.clear();
      vector<int> indexMuTau; indexMuTau.clear();
      vector<int> indexETau; indexETau.clear();
      vector<TLorentzVector> lvW; lvW.clear();
      vector<TLorentzVector> lvNu; lvNu.clear(); 
      vector<TLorentzVector> lvMu; lvMu.clear();
      vector<TLorentzVector> lvE; lvE.clear();
      vector<TLorentzVector> lvMuTau; lvMuTau.clear();
      vector<TLorentzVector> lvETau; lvETau.clear();
      
      // int nGenMuons = 0;
      // int nGenElectrons = 0;
      for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

       	float pxGen = analysisTree.genparticles_px[igen];
       	float pyGen = analysisTree.genparticles_py[igen];
       	float pzGen = analysisTree.genparticles_pz[igen];
	float enGen = analysisTree.genparticles_e[igen];

	TLorentzVector lorentzVector; lorentzVector.SetPxPyPzE(pxGen,pyGen,pzGen,enGen);

       	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && 
	    (analysisTree.genparticles_status[igen]==62||analysisTree.genparticles_status[igen]==52)) { 
       	  indexW.push_back(igen);
	  lvW.push_back(lorentzVector);
	}

       	if ((fabs(analysisTree.genparticles_pdgid[igen])==12 
       	     ||fabs(analysisTree.genparticles_pdgid[igen])==14)
       	    && analysisTree.genparticles_info[igen]==2 && analysisTree.genparticles_status[igen]==1) {
       	  indexNu.push_back(igen);
	  lvNu.push_back(lorentzVector);
	}

       	if (fabs(analysisTree.genparticles_pdgid[igen])==13 &&
	    analysisTree.genparticles_info[igen]==2 && analysisTree.genparticles_status[igen]==1) {
	  indexMu.push_back(igen);
	  lvMu.push_back(lorentzVector);
	}

       	if (fabs(analysisTree.genparticles_pdgid[igen])==13 &&
	    analysisTree.genparticles_info[igen]==5 && analysisTree.genparticles_status[igen]==1) {
	  indexMuTau.push_back(igen);
	  lvMuTau.push_back(lorentzVector);
	}

     	if (fabs(analysisTree.genparticles_pdgid[igen])==11 &&
	    analysisTree.genparticles_info[igen]==2 && analysisTree.genparticles_status[igen]==1 ) {
	  indexE.push_back(igen);
	  lvE.push_back(lorentzVector);
	}

     	if (fabs(analysisTree.genparticles_pdgid[igen])==11 &&
	    analysisTree.genparticles_info[igen]==5 && analysisTree.genparticles_status[igen]==1 ) {
	  indexETau.push_back(igen);
	  lvETau.push_back(lorentzVector);
	}

      }
      // w bosons;
      //      std::cout << "Number of W " << indexW.size() << std::endl;

      vector<int> indexEfromW; indexEfromW.clear();
      vector<int> indexMufromW; indexMufromW.clear();

      //      for (unsigned int iW=0; iW<indexW.size(); ++iW) {
      //	int index = indexW.at(iW);
      //	int pdgId = analysisTree.genparticles_pdgid[index];
      //	if (pdgId==24)
      //	  printf("W +    : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
      //		 lvW[iW].Px(),lvW[iW].Py(),lvW[iW].Pz(),lvW[iW].M());
      //	else
      //	  printf("W -    : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
      //                 lvW[iW].Px(),lvW[iW].Py(),lvW[iW].Pz(),lvW[iW].M());
      //      }

      for (unsigned int iNu=0; iNu<indexNu.size(); ++iNu) {

	int IndexNu = indexNu.at(iNu);
	int PdgIdNu = analysisTree.genparticles_pdgid[IndexNu];

	for (unsigned int iE=0; iE<indexE.size(); ++iE) {
	  int IndexE = indexE.at(iE);
	  int PdgIdE = analysisTree.genparticles_pdgid[IndexE];
	  bool isWplus = (PdgIdE==-11) && (PdgIdNu==12);
	  bool isWminus = (PdgIdE==11) && (PdgIdNu==-12);
	  if (isWplus||isWminus) {
	    TLorentzVector sumLV = lvNu[iNu]+lvE[iE];
	    if (sumLV.M()>60) indexEfromW.push_back(IndexE);
	    //	    if (isWplus)
	    //	      printf("(e,v)+ : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
	    //		     sumLV.Px(),sumLV.Py(),sumLV.Pz(),sumLV.M());
	    //	    if (isWminus)
	    //              printf("(e,v)- : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
	    //                     sumLV.Px(),sumLV.Py(),sumLV.Pz(),sumLV.M());
	  }
	}

	for (unsigned int iM=0; iM<indexMu.size(); ++iM) {
	  int IndexM = indexMu.at(iM);
	  int PdgIdM = analysisTree.genparticles_pdgid[IndexM];
	  bool isWplus = (PdgIdM==-13) && (PdgIdNu==14);
	  bool isWminus = (PdgIdM==13) && (PdgIdNu==-14);
	  if (isWplus||isWminus) {
	    TLorentzVector sumLV = lvNu[iNu]+lvMu[iM];
	    if (sumLV.M()>60) indexMufromW.push_back(IndexM);
	    //	    if (isWplus)
	    //	      printf("(m,v)+ : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
	    //		     sumLV.Px(),sumLV.Py(),sumLV.Pz(),sumLV.M());
	    //	    if (isWminus)
	    //              printf("(m,v)- : p = %7.1f, %7.1f, %7.1f,  mass = %4.1f\n",
	    //                     sumLV.Px(),sumLV.Py(),sumLV.Pz(),sumLV.M());
	  }
	}

      }	

      if (indexETau.size()>0&&indexMuTau.size()>0) {
	DYJetsDiLeptonH->Fill(0.);

	DYJetsDiLeptonElePtH->Fill(lvETau[0].Pt());
        DYJetsDiLeptonEleEtaH->Fill(lvETau[0].Eta());
        DYJetsDiLeptonMuPtH->Fill(lvMuTau[0].Pt());
        DYJetsDiLeptonMuEtaH->Fill(lvMuTau[0].Eta());

	if (fabs(lvMuTau[0].Eta())<etaMuonCut && fabs(lvETau[0].Eta())<etaElectronCut) {
          bool isAccepted = ( lvMuTau[0].Pt()>ptMuonLowCut && lvETau[0].Pt()>ptElectronHighCut) ||
            ( lvMuTau[0].Pt()>ptMuonHighCut && lvETau[0].Pt()>ptElectronLowCut ) ;
          if (isAccepted)
            DYJetsDiLeptonAccH->Fill(0.);
        }



      }


      if (indexEfromW.size()>0&&indexMufromW.size()>0) {
	TTJetsDiLeptonH->Fill(0.);
	int IndexM = indexMufromW.at(0);
	int IndexE = indexEfromW.at(0);

	float pxGenM = analysisTree.genparticles_px[IndexM];
        float pyGenM = analysisTree.genparticles_py[IndexM];
        float pzGenM = analysisTree.genparticles_pz[IndexM];
        float enGenM = analysisTree.genparticles_e[IndexM];

        TLorentzVector fourVecM; fourVecM.SetPxPyPzE(pxGenM,pyGenM,pzGenM,enGenM);

	float pxGenE = analysisTree.genparticles_px[IndexE];
        float pyGenE = analysisTree.genparticles_py[IndexE];
        float pzGenE = analysisTree.genparticles_pz[IndexE];
        float enGenE = analysisTree.genparticles_e[IndexE];

        TLorentzVector fourVecE; fourVecE.SetPxPyPzE(pxGenE,pyGenE,pzGenE,enGenE);

	TTJetsDiLeptonElePtH->Fill(fourVecE.Pt());
	TTJetsDiLeptonEleEtaH->Fill(fourVecE.Eta());
	TTJetsDiLeptonMuPtH->Fill(fourVecM.Pt());
	TTJetsDiLeptonMuEtaH->Fill(fourVecM.Eta());

	if (fabs(fourVecM.Eta())<etaMuonCut && fabs(fourVecE.Eta())<etaElectronCut) {
	  bool isAccepted = ( fourVecM.Pt()>ptMuonLowCut && fourVecE.Pt()>ptElectronHighCut) ||
	    ( fourVecM.Pt()>ptMuonHighCut && fourVecE.Pt()>ptElectronLowCut ) ;
	  if (isAccepted)
	    TTJetsDiLeptonAccH->Fill(0.);
	}


      }

      continue;
      // vertex cuts

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;

      
      // electron selection

      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	electronPtAllH->Fill(analysisTree.electron_pt[ie],weight);
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ie];
	if (relIso>isoElectronHighCut) continue;
	if (relIso<isoElectronLowCut) continue;
	bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[ie],
						analysisTree.electron_mva_id_nontrigPhys14[ie]);
	if (!electronMvaId&&applyElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]!=0&&applyElectronId) continue;
	electrons.push_back(ie);
      }

      // muon selection

      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	muonPtAllH->Fill(analysisTree.muon_pt[im],weight);
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	float neutralIso = 
	  analysisTree.muon_neutralHadIso[im] + 
	  analysisTree.muon_photonIso[im] - 
	  0.5*analysisTree.muon_puIso[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.muon_chargedHadIso[im] + neutralIso;
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonHighCut) continue;
	if (relIso<isoMuonLowCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	muons.push_back(im);
      }

      if (electrons.size()==0) continue;
      if (muons.size()==0) continue;

      // selecting muon and electron pair (OS or SS);
      float ptScalarSum = -1;
      float dRleptons = -1;
      int electronIndex = -1;
      int muonIndex = -1;
      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	int eIndex = electrons[ie];
	for (unsigned int im=0; im<muons.size(); ++im) {
	  int mIndex = muons[im];
	  float qProd = analysisTree.electron_charge[eIndex]*analysisTree.muon_charge[mIndex];
	  if (oppositeSign && qProd>0) continue;
	  if (!oppositeSign && qProd<0) continue;
	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;

	  bool kinematics = 
	    (analysisTree.electron_pt[eIndex]>ptElectronLowCut && analysisTree.muon_pt[mIndex]>ptMuonHighCut) || 
	    (analysisTree.electron_pt[eIndex]>ptElectronHighCut && analysisTree.muon_pt[mIndex]>ptMuonLowCut);

	  if (!kinematics) continue;

	  float sumPt = analysisTree.electron_pt[eIndex] + analysisTree.muon_pt[mIndex];
	  if (sumPt>ptScalarSum) {
	    ptScalarSum = sumPt;
	    dRleptons = dR;
	    electronIndex = ie;
	    muonIndex = im;
	  }
	  
	}
      }

      if (ptScalarSum<0) continue;

      // Trigger matching

      bool triggerMatch = false;

      if (analysisTree.hltriggerresults_second[5]==1) { // HLT_Mu23_Ele12
	if (analysisTree.muon_pt[muonIndex]>ptMuonHighCut&&analysisTree.electron_pt[electronIndex]>ptElectronLowCut) {
	  bool muonMatch = false;
	  bool electronMatch = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][2]) { // Mu23 Leg
	      float dRtrig = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<0.3) {
		muonMatch = true;
	      }
	    }
	    if (analysisTree.trigobject_filters[iT][6]) { // Ele12 Leg
	      float dRtrig = deltaR(analysisTree.electron_eta[electronIndex],analysisTree.electron_phi[electronIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<0.3) {
		electronMatch = true;
	      }
	    }
	  }
	  if (muonMatch&&electronMatch)
	    triggerMatch = true;
	}
      }

      if (analysisTree.hltriggerresults_second[6]==1) { // HLT_Mu8_Ele23
	if (analysisTree.muon_pt[muonIndex]>ptMuonLowCut&&analysisTree.electron_pt[electronIndex]>ptElectronHighCut) {
	  bool muonMatch = false;
	  bool electronMatch = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    if (analysisTree.trigobject_filters[iT][3]) { // Mu8 Leg
	      float dRtrig = deltaR(analysisTree.muon_eta[muonIndex],analysisTree.muon_phi[muonIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<0.3) {
		muonMatch = true;
	      }
	    }
	    if (analysisTree.trigobject_filters[iT][7]) { // Ele23 Leg
	      float dRtrig = deltaR(analysisTree.electron_eta[electronIndex],analysisTree.electron_phi[electronIndex],
				    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	      if (dRtrig<0.3) {
		electronMatch = true;
	      }
	    }
	  }
	  if (muonMatch&&electronMatch)
	    triggerMatch = true;
	}
      }

      if (!triggerMatch&&applyTriggerMatch) continue;

      // computation of kinematic variables

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
						    analysisTree.electron_py[electronIndex],
						    analysisTree.electron_pz[electronIndex],
						    electronMass);

      TLorentzVector dileptonLV = muonLV + electronLV;
      float dileptonMass = dileptonLV.M();
      float dileptonPt = dileptonLV.Pt();
      float dileptonEta = dileptonLV.Eta();

      float ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

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

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + electronLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + electronLV.Py();
      
      float vectorVisX = muonLV.Px() + electronLV.Px();
      float vectorVisY = muonLV.Py() + electronLV.Py();

      // computation of DZeta variable
      float PZeta = vectorX*zetaX + vectorY*zetaY;
      float PVisZeta = vectorVisX*zetaX + vectorVisY*zetaY;
      float DZeta = PZeta - 1.85*PVisZeta;

      // computation of MT variable
      float dPhi = dPhiFrom2P( dileptonLV.Px(), dileptonLV.Py(),
			       analysisTree.pfmet_ex,  analysisTree.pfmet_ey );

      float MT = TMath::Sqrt(2*dileptonPt*ETmiss*(1-TMath::Cos(dPhi)));

      // number of jets
      int nJets30 = 0;
      int nJets20 = 0;
      int nBLooseJets = 0;
      int nBMediumJets = 0;

      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta>jetEtaCut) continue;

	float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   muonLV.Eta(),muonLV.Phi());
	if (dR1<dRJetLeptonCut) continue;

	float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   electronLV.Eta(),electronLV.Phi());

	if (dR2<dRJetLeptonCut) continue;

	// puJetId
	float energy = analysisTree.pfjet_e[jet];
        float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
        float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
        float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
        float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
        float chm = analysisTree.pfjet_chargedmulti[jet];
        float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
	bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>2.4 || (elf<0.99 && chf>0 && chm>0));
	if (!isPFJetId) continue;

	if (analysisTree.pfjet_pt[jet]>jetPtHighCut)
	  nJets30++;

	if (analysisTree.pfjet_pt[jet]>jetPtLowCut)
	  nJets20++;

	if (analysisTree.pfjet_pt[jet]>jetPtLowCut && absJetEta<bJetEtaCut) {
	  if (analysisTree.pfjet_btag[jet][6]>btagLooseCut)
	    nBLooseJets++;
	  if (analysisTree.pfjet_btag[jet][6]>btagMediumCut)
	    nBMediumJets++;
	}

      }

      // filling histograms after dilepton selection

      electronPtH->Fill(electronLV.Pt(),weight);
      electronEtaH->Fill(electronLV.Eta(),weight);
      
      muonPtH->Fill(muonLV.Pt(),weight);
      muonEtaH->Fill(muonLV.Eta(),weight);
      
      dileptonMassH->Fill(dileptonMass,weight);
      dileptonPtH->Fill(dileptonPt,weight);
      dileptonEtaH->Fill(dileptonEta,weight);
      dileptondRH->Fill(dRleptons,weight);
      
      ETmissH->Fill(ETmiss,weight);
      MtH->Fill(MT,weight);
      DZetaH->Fill(DZeta,weight);

      nJets30H->Fill(float(nJets30),weight);
      nJets20H->Fill(float(nJets20),weight);
      nBLooseJetsH->Fill(float(nBLooseJets),weight);
      nBMediumJetsH->Fill(float(nBMediumJets),weight);

      // topological cut
      if (DZeta>dZetaCut) {
      
	electronPtDZetaCutH->Fill(electronLV.Pt(),weight);
	electronEtaDZetaCutH->Fill(electronLV.Eta(),weight);
      
	muonPtDZetaCutH->Fill(muonLV.Pt(),weight);
	muonEtaDZetaCutH->Fill(muonLV.Eta(),weight);
	
	dileptonMassDZetaCutH->Fill(dileptonMass,weight);
	dileptonPtDZetaCutH->Fill(dileptonPt,weight);
	dileptonEtaDZetaCutH->Fill(dileptonEta,weight);
	dileptondRDZetaCutH->Fill(dRleptons,weight);
	
	ETmissDZetaCutH->Fill(ETmiss,weight);
	MtDZetaCutH->Fill(MT,weight);
	DZetaDZetaCutH->Fill(DZeta,weight);

	nJets30DZetaCutH->Fill(float(nJets30),weight);
	nJets20DZetaCutH->Fill(float(nJets20),weight);
	nBLooseJetsDZetaCutH->Fill(float(nBLooseJets),weight);
	nBMediumJetsDZetaCutH->Fill(float(nBMediumJets),weight);

      }

      if (ETmiss<MEtCut) {

	electronPtMEtCutH->Fill(electronLV.Pt(),weight);
	electronEtaMEtCutH->Fill(electronLV.Eta(),weight);
      
	muonPtMEtCutH->Fill(muonLV.Pt(),weight);
	muonEtaMEtCutH->Fill(muonLV.Eta(),weight);
	
	dileptonMassMEtCutH->Fill(dileptonMass,weight);
	dileptonPtMEtCutH->Fill(dileptonPt,weight);
	dileptonEtaMEtCutH->Fill(dileptonEta,weight);
	dileptondRMEtCutH->Fill(dRleptons,weight);
	
	ETmissMEtCutH->Fill(ETmiss,weight);
	MtMEtCutH->Fill(MT,weight);
	DZetaMEtCutH->Fill(DZeta,weight);

	nJets30MEtCutH->Fill(float(nJets30),weight);
	nJets20MEtCutH->Fill(float(nJets20),weight);
	nBLooseJetsMEtCutH->Fill(float(nBLooseJets),weight);
	nBMediumJetsMEtCutH->Fill(float(nBMediumJets),weight);

      }

      if (ETmiss<MEtCut&&DZeta>dZetaCut) {

	electronPtSelH->Fill(electronLV.Pt(),weight);
	electronEtaSelH->Fill(electronLV.Eta(),weight);
      
	muonPtSelH->Fill(muonLV.Pt(),weight);
	muonEtaSelH->Fill(muonLV.Eta(),weight);
	
	dileptonMassSelH->Fill(dileptonMass,weight);
	dileptonPtSelH->Fill(dileptonPt,weight);
	dileptonEtaSelH->Fill(dileptonEta,weight);
	dileptondRSelH->Fill(dRleptons,weight);
	
	ETmissSelH->Fill(ETmiss,weight);
	MtSelH->Fill(MT,weight);
	DZetaSelH->Fill(DZeta,weight);

	nJets30SelH->Fill(float(nJets30),weight);
	nJets20SelH->Fill(float(nJets20),weight);
	nBLooseJetsSelH->Fill(float(nBLooseJets),weight);
	nBMediumJetsSelH->Fill(float(nBMediumJets),weight);


      }
      
      if (ETmiss>MEtCutTTJets&&nBLooseJets>0) {
	
	electronPtTTJetsSelH->Fill(electronLV.Pt(),weight);
	electronEtaTTJetsSelH->Fill(electronLV.Eta(),weight);
      
	muonPtTTJetsSelH->Fill(muonLV.Pt(),weight);
	muonEtaTTJetsSelH->Fill(muonLV.Eta(),weight);
	
	dileptonMassTTJetsSelH->Fill(dileptonMass,weight);
	dileptonPtTTJetsSelH->Fill(dileptonPt,weight);
	dileptonEtaTTJetsSelH->Fill(dileptonEta,weight);
	dileptondRTTJetsSelH->Fill(dRleptons,weight);
	
	ETmissTTJetsSelH->Fill(ETmiss,weight);
	MtTTJetsSelH->Fill(MT,weight);
	DZetaTTJetsSelH->Fill(DZeta,weight);

	nJets30TTJetsSelH->Fill(float(nJets30),weight);
	nJets20TTJetsSelH->Fill(float(nJets20),weight);
	nBLooseJetsTTJetsSelH->Fill(float(nBLooseJets),weight);
	nBMediumJetsTTJetsSelH->Fill(float(nBMediumJets),weight);

      }


      //      std::cout << std::endl;
      
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



