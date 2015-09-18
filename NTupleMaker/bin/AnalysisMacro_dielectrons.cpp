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

struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject;

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on muons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronCut     = cfg.get<float>("isoElectronCut");
  const float isoElectronTightCut = cfg.get<float>("isoElectronTightCut");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // Run range
  const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // **** end of configuration


  // file name and tree name
	char ff[100];

	sprintf(ff,"%s/%s",argv[3],argv[2]);

  // file name and tree name
  std::string rootFileName(argv[2]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  string SelectionSign="elel";

  TString era=argv[3];
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  


  // output fileName with histograms
  TFile * file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());
  TH1D * hxsec = new TH1D("xsec","",1,0,0);
  TH1D * histWeights = new TH1D("histWeights","",1,0,0);
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

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

  TH1F * ZMassAllElectronsH = new TH1F("ZMassAllElectronsH","",60,60,120);
  TH1F * ZMassAllElectronsDRCutH = new TH1F("ZMassAllElectronsDRCutH","",400,60,120);
  TH1F * ZMassIdLooseElectronsH = new TH1F("ZMassIdLooseElectronsH","",60,60,120);
  TH1F * ZMassIdLooseElectronsDRCutH = new TH1F("ZMassIdLooseElectronsDRCutH","",60,60,120);
  TH1F * ZMassIdTightElectronsH = new TH1F("ZMassIdTightElectronsH","",60,60,120);
  TH1F * ZMassIdTightElectronsDRCutH = new TH1F("ZMassIdTightElectronsDRCutH","",60,60,120);
  TH1F * ZMassIsoElectronsH = new TH1F("ZMassIsoElectronsH","",60,60,120);
  TH1F * ZMassIsoElectronsDRCutH = new TH1F("ZMassIsoElectronsDRCutH","",60,60,120);
  TH1F * ZMassIsoTightElectronsH = new TH1F("ZMassIsoTightElectronsH","",60,60,120);
  TH1F * ZMassIsoTightElectronsDRCutH = new TH1F("ZMassIsoTightElectronsDRCutH","",60,60,120);

  Float_t XSec=-1;
  Float_t xs,fact,fact2;
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

 
  ifstream ifs("xsecs");
  string line;

  while(std::getline(ifs, line)) // read one line from ifs
    {
		
      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt;
      iss >> dt >> xs >> fact >> fact2;
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      //if (dt==argv[2]) {
      //if (std::string::npos != dt.find(argv[2])) {
      if (  dt == argv[2]) {
	XSec= xs*fact*fact2;
	cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
      }
        
    }

  if (XSec<0) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}

  std::vector<unsigned int> allRuns; allRuns.clear();

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
	bool isData = false;
     
	if (XSec == 1) isData = true;
	
	if (isData){
      if (analysisTree.event_run<RunRangeMin) continue;
      if (analysisTree.event_run>RunRangeMax) continue;

      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;

      bool isNewRun = true;
      for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	if (analysisTree.event_run==allRuns.at(iR)) {
	  isNewRun = false;
	  break;
	}
      }

      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);
	}
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
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {
	allElectrons.push_back(im);
	if (analysisTree.electron_pt[im]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut) continue;
	if (!analysisTree.electron_pass_conversion[im]) continue;
	if (analysisTree.electron_nmissinginnerhits[im]!=0) continue;
	bool electronMvaIdTightX = electronMvaIdTight(analysisTree.electron_superclusterEta[im],
						      analysisTree.electron_mva_id_nontrigPhys14[im]);
	bool electronMvaIdLooseX = electronMvaIdLoose(analysisTree.electron_superclusterEta[im],
						      analysisTree.electron_mva_id_nontrigPhys14[im]);
	if (!electronMvaIdLooseX) continue;
	idLooseElectrons.push_back(im);
	if (!electronMvaIdTightX) continue;
	idTightElectrons.push_back(im);
	float absIso = analysisTree.electron_chargedHadIso[im];
	float relIso = absIso/analysisTree.electron_pt[im];
	if (relIso>isoElectronCut) continue;
	isoElectrons.push_back(im);
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[im] + 
	  analysisTree.electron_photonIso[im] - 
	  0.5*analysisTree.electron_puIso[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	absIso += neutralIso;
	relIso = absIso/analysisTree.electron_pt[im];
	if (relIso>isoElectronTightCut) continue;
	isoTightElectrons.push_back(im);
      }

      // std::cout << "allElectrons : " << allElectrons.size() << std::endl;
      // std::cout << "idElectrons  : " << idElectrons.size() << std::endl;
      // std::cout << "isoElectrons : " << isoElectrons.size() << std::endl;

      //      continue;

      bool isAllElectronsPair = false;
      if (allElectrons.size()>1) {
	//	std::cout << "allElectrons : " << allElectrons.size() << std::endl;
	for (unsigned int im1=0; im1<allElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<allElectrons.size(); ++im2) {
	    unsigned int index1 = allElectrons[im1];
	    unsigned int index2 = allElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    if (q1*q2<0) {
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
	//	std::cout << "allElectrons : " << allElectrons.size() << "  :  OK" << std::endl;
      }

      //      continue;

      bool isIdLooseElectronsPair = false;
      if (idLooseElectrons.size()>1) {
	for (unsigned int im1=0; im1<idLooseElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<idLooseElectrons.size(); ++im2) {
	    unsigned int index1 = idLooseElectrons[im1];
	    unsigned int index2 = idLooseElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    if (q1*q2<0) {
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
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassIdLooseElectronsDRCutH->Fill(mass,weight);
		YpsilonMassIdLooseElectronsDRCutH->Fill(mass,weight);
		ZMassIdLooseElectronsDRCutH->Fill(mass,weight);
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
	    if (q1*q2<0) {
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
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		JPsiMassIdTightElectronsDRCutH->Fill(mass,weight);
		YpsilonMassIdTightElectronsDRCutH->Fill(mass,weight);
		ZMassIdTightElectronsDRCutH->Fill(mass,weight);
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
	    if (q1*q2<0) {
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
      if (isoTightElectrons.size()>1) {
	for (unsigned int im1=0; im1<isoTightElectrons.size()-1; ++im1) {
	  for (unsigned int im2=im1+1; im2<isoTightElectrons.size(); ++im2) {
	    unsigned int index1 = isoTightElectrons[im1];
	    unsigned int index2 = isoTightElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    if (q1*q2<0) {
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
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	      if (dR>dRleptonsCut) {
		ZMassIsoTightElectronsDRCutH->Fill(mass,weight);
	      }

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
    // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
 
  file->cd(SelectionSign.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  histWeights->Write();
  file->Write();
  file->Close();
  delete file;
  
  
}



