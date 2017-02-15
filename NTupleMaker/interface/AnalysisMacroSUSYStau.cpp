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
#include "Riostream.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

using namespace std;
const float MuMass = 0.105658367;




bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


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


double DeltaPhi(TLorentzVector METV, TLorentzVector LeptonV){
  TLorentzVector Ws = METV + LeptonV;

        //Delta phi between W and Lep
        //standard root defintion (+ fabs)takes care of getting values between 0 and pi
        double DelPhiWLep = fabs(Ws.DeltaPhi(LeptonV));
        //alternative definiton with the same result, if you want to cross check
      	 Double_t DelPhiWLepAlt = (Ws.Phi() - LeptonV.Phi());
        if (DelPhiWLepAlt > TMath::Pi()) DelPhiWLepAlt -= 2*TMath::Pi();
        if (DelPhiWLepAlt <= -TMath::Pi()) DelPhiWLepAlt += 2*TMath::Pi();
        DelPhiWLepAlt = fabs(DelPhiWLepAlt);
		
	return DelPhiWLep;

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








const int CutNumb = 5;
TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutNumb,0.5,CutNumb+0.5);
TH1D *hHT[CutNumb];
TH1D *hST[CutNumb];
TH1D *h0JetpT[CutNumb];
TH1D *hnJet[CutNumb];
TH1D *hnBJet[CutNumb];
TH1D *hnEl[CutNumb];
TH1D *hnMu[CutNumb];
TH1D *hLeppt[CutNumb];
TH1D *hElpt[CutNumb];
TH1D *hMupt[CutNumb];
TH1D *hLepeta[CutNumb];
TH1D *hEleta[CutNumb];
TH1D *hMueta[CutNumb];
TH1D *hMET[CutNumb];
TH1D *hnOver[CutNumb];
TH1D *hdPhiMETLep[CutNumb];
TH1D *hdPhiJMET[CutNumb];
TH1D *hToppT[CutNumb];
TH1D *hnTop[CutNumb];
TH1D *hnW[CutNumb];
TH1D *hWTagpT[CutNumb];
TH1D *hWmassTagpT[CutNumb];
TH1D *hWTagMass[CutNumb];
TH1D *hWmassTagMass[CutNumb];
TH1D *hnWmass[CutNumb];
TH1D *hMT[CutNumb];
TH1D *hDZeta[CutNumb];
TH1D *hel_miniISO[CutNumb];
TH1D *hmu_miniISO[CutNumb];
  


string CutList[CutNumb];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNumber){
    for(int cj = 0; cj < CutNumber; cj++)
    {
        CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
        TString cutName=CutList[cj];
        TString nCut;
        nCut.Form("%d",cj);

        hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,400,0.0,4000.0);
        hHT[cj]->Sumw2();
        hST[cj] = new TH1D ("ST_"+nCut,"ST "+cutName,400,0.0,4000.0);
        hST[cj]->Sumw2();
        h0JetpT[cj] = new TH1D ("0JetpT_"+nCut,"0JetpT "+cutName,200,0.0,2000.0);
        h0JetpT[cj]->Sumw2();
        hnJet[cj] = new TH1D ("nJet_"+nCut,"nJet "+cutName,20,0,20);
        hnJet[cj]->Sumw2();
        hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,20,0,20);
        hnBJet[cj]->Sumw2();
        hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,0,10);
        hnEl[cj]->Sumw2();
        hLeppt[cj] = new TH1D ("LeppT_"+nCut,"Lep pT "+cutName,100,0,1000);
        hLeppt[cj]->Sumw2();
        hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,100,0,1000);
        hElpt[cj]->Sumw2();
        hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,0,10);
        hnMu[cj]->Sumw2();
        hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,100,0,1000);
        hMupt[cj]->Sumw2();
        hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
        hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,100,-4,4);
        hEleta[cj]->Sumw2();
        hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,100,-4,4);
        hMueta[cj]->Sumw2();
        hLepeta[cj] = new TH1D ("Lepeta_"+nCut,"Lep eta "+cutName,100,-4,4);
        hLepeta[cj]->Sumw2();
        hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,200.0,0.0,4000.0);
        hMET[cj]->Sumw2();
        hdPhiMETLep[cj] = new TH1D("dPhiMETLep_"+nCut,"dPhiMETLep "+cutName,64,0.0,3.2);
        hdPhiMETLep[cj]->Sumw2();
        hdPhiJMET[cj] = new TH1D("dPhiJMET_"+nCut,"dPhiJMET "+cutName,64,0.0,3.2);
        hdPhiJMET[cj]->Sumw2();
        hToppT[cj] = new TH1D ("ToppT_"+nCut,"ToppT "+cutName,200,0.0,2000.0);
        hToppT[cj]->Sumw2();
        hnTop[cj] = new TH1D ("nTop_"+nCut,"nTop "+cutName,20,0,20);
        hnTop[cj]->Sumw2();
        hnW[cj] = new TH1D ("nW_"+nCut,"nW "+cutName,20,0,20);
        hnW[cj]->Sumw2();
        hWTagpT[cj] = new TH1D ("WpT_"+nCut,"WpT "+cutName,200,0.0,2000.0);
        hWTagpT[cj]->Sumw2();
        hWmassTagpT[cj] = new TH1D ("WmasspT_"+nCut,"WmasspT "+cutName,200,0.0,2000.0);
        hWmassTagpT[cj]->Sumw2();
        hWTagMass[cj] = new TH1D ("WMass_"+nCut,"WMass "+cutName,200,0.0,2000.0);
        hWTagMass[cj]->Sumw2();
        hWmassTagMass[cj] = new TH1D ("WmassMass_"+nCut,"WmassMass "+cutName,200,0.0,2000.0);
        hWmassTagMass[cj]->Sumw2();
        hnWmass[cj] = new TH1D ("nWmass_"+nCut,"nWmass "+cutName,20,0,20);
        hnWmass[cj]->Sumw2();
        hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,40,0,200);
        hMT[cj]->Sumw2();
 /*       hDZeta[CutNumb]= new TH1D ("DZeta_"+nCut,"DZeta "+cutName,300,-400,200);
        hDZeta[cj]->Sumw2();
*/
        hel_miniISO[cj]= new TH1D ("elminiISO_"+nCut,"elminiISO "+cutName,50,0,5);;
        hel_miniISO[cj]->Sumw2();
        hmu_miniISO[cj]= new TH1D ("muminiISO_"+nCut,"muminiISO "+cutName,50,0,5);;
        hmu_miniISO[cj]->Sumw2();
    }
}


//void //FillMainHists(int CutIndex, Double_t EvWeight, bool FillBJets = true){
//void FillMainHists(int CutIndex, Double_t EvWeight, int nEl, int nMu, int nJets){
//void FillMainHists(int CutIndex, Double_t EvWeight, TLorentz ElV, TLorentz MuV, JetsV){
//void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  *ElV, vector<TLorentzVector>  *MuV,vector<TLorentzVector>  *JetsV,Float_t MET){
void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV,vector<TLorentzVector>  JetsV, TLorentzVector  MetV){
	//void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  *JetsV){

	hnJet[CutIndex]->Fill(JetsV.size(),EvWeight);
        hnEl[CutIndex]->Fill(ElV.size(),EvWeight);
        hnMu[CutIndex]->Fill(MuV.size(),EvWeight);
   //     if (JetsV.size() > 0) h0JetpT[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);
  //  if(FillBJets){
  //      hnBJet[CutIndex]->Fill(Obj.nBJetGood,EvWeight);
  //  }
    if (ElV.size() > 0)
    {
     hElpt[CutIndex]->Fill(ElV.at(0).Pt(),EvWeight);
     hLeppt[CutIndex]->Fill(ElV.at(0).Pt(),EvWeight);
     
 //    for (unsigned int nel = 0; nel<ElV.size();nel++){
 //    hel_miniISO[CutIndex]->Fill(analysisTree.electron_miniISO[nel],EvWeight);
 //    }
     Float_t dPhi=dPhiFrom2P( ElV.at(0).Px(), ElV.at(0).Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      		
      Float_t MT = TMath::Sqrt(2*ElV.at(0).Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
    
    }


    if (MuV.size() > 0)
    {
     TLorentzVector WBos = MetV + MuV.at(0);
     hMupt[CutIndex]->Fill(MuV.at(0).Pt(),EvWeight);
     hLeppt[CutIndex]->Fill(MuV.at(0).Pt(),EvWeight);
  //   for (unsigned int nmu = 0; nmu<MuV.size();nmu++){
  //   hmu_miniISO[CutIndex]->Fill(analysisTree.muon_miniISO[nmu],EvWeight);
  //  }
     Float_t dPhi=dPhiFrom2P( MuV.at(0).Px(), MuV.at(0).Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      Float_t MT = TMath::Sqrt(2*MuV.at(0).Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
     //cout<<"  "<<dPhi<<"  "<<fabs(WBos.DeltaPhi(MuV.at(0)))<<endl;
    }


    float sumpT=0;
   if (JetsV.size()>0){
    for (unsigned int ij=0;ij<JetsV.size();ij++){
         sumpT+=JetsV.at(ij).Pt();
     Float_t dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), MetV.Px(),  MetV.Py() );
     hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);
	 }
     hHT[CutIndex]->Fill(sumpT,EvWeight);
    }

      hMET[CutIndex]->Fill(MetV.Pt(),EvWeight);
     // hMT[CutIndex]->Fill(MT,EvWeight);
    //  hDZeta->Fill(DZeta,EvWeight);

   /*
    hST[CutIndex]->Fill(Obj.ST,EvWeight);
    hnTop[CutIndex]->Fill(Obj.nTopTagJetGood,EvWeight);
    hdPhiJMET[CutIndex]->Fill(Obj.minDelPhiJMet,EvWeight);
    if(Obj.nTopTagJetGood>0)hToppT[CutIndex]->Fill(Obj.goodTopTagJet[0].Pt(),EvWeight);
    hnW[CutIndex]->Fill(Obj.nWTagJetGood,EvWeight);
    if(Obj.nWTagJetGood>0){
    hWTagpT[CutIndex]->Fill(Obj.goodWTagJet[0].Pt(),EvWeight);
    hWTagMass[CutIndex]->Fill(Obj.goodWTagJet[0].M(),EvWeight);

    }
    hnWmass[CutIndex]->Fill(Obj.nWmassTagJetGood,EvWeight);
    if(Obj.nWmassTagJetGood>0){
     hWmassTagpT[CutIndex]->Fill(Obj.goodWmassTagJet[0].Pt(),EvWeight);
    hWmassTagMass[CutIndex]->Fill(Obj.goodWmassTagJet[0].M(),EvWeight);
   }*/
	
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
  string SelectionSign=argv[3];

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
  
  // kinematic cuts on Jets
  const float etaJetCut   = cfg.get<float>("etaJetCut");
  const float ptJetCut   = cfg.get<float>("ptJetCut");
  
  
  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZetaCut       = cfg.get<float>("dZetaCut");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  
  const float Lumi   = cfg.get<float>("Lumi");

  const float bTag 	     = cfg.get<float>("bTag");
  // **** end of configuration
 
   //TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   //dir.ReplaceAll("basic.C","");
   //dir.ReplaceAll("/./","/");
   //ifstream in;
   Float_t XSec=-1;
	ifstream ifs("xsecs");

	string line;

	while(std::getline(ifs, line)) // read one line from ifs
	{
    	istringstream iss(line); // access line as a stream

    // we only need the first two columns
    string dt;
    Float_t xs;

    ifs >> dt >> xs; // no need to read further

	
    cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<endl;
     //if (dt==argv[2]) {
     if (std::string::npos != dt.find(argv[2])) {
	     XSec=xs; 
	     cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<endl;
 		}
        
	}

	if (XSec<0) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}





	

string	CutList[CutNumb]={"No cut","Trigger",SelectionSign+" dR > "+to_string(dRleptonsCut),"1 >=b-Jets ("+to_string(ptJetCut)+" GeV)","lep SumpT> 0","dPhi > 1"};

//string	CutListt[CutNumb]={"No cut","Trigger",SelectionSign+" dR > "+to_string(dRleptonsCut),"3 >=Jets ("+to_string(ptJetCut)+" GeV)","lep SumpT> 0","dZeta >"+to_string(dZetaCut)};
	
if (SelectionSign == "1l" || SelectionSign == "1L"){
CutList[2] = "==1l";
CutList[4]= "dPhi > 1";
}


//CutList[CutNumb]=CutListt[CutNumb];

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"update");
  file->mkdir(SelectionSign.c_str());
  file->cd(SelectionSign.c_str());
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);

  TH1D * muonPtAllH = new TH1D("muonPtAllH","",40,0,200);
  TH1D * electronPtAllH = new TH1D("electronPtAllH","",40,0,200);

  // histograms (dilepton selection)
  TH1D * electronPtH  = new TH1D("electronPtH","",40,0,200);
  TH1D * electronEtaH = new TH1D("electronEtaH","",50,-2.5,2.5); 
  TH1D * muonPtH  = new TH1D("muonPtH","",40,0,200);
  TH1D * muonEtaH = new TH1D("muonEtaH","",50,-2.5,2.5); 

  TH1D * dileptonMassH = new TH1D("dileptonMassH","",40,0,200);
  TH1D * dileptonPtH = new TH1D("dileptonPtH","",40,0,200);
  TH1D * dileptonEtaH = new TH1D("dileptonEtaH","",100,-5,5);
  TH1D * dileptondRH = new TH1D("dileptondRH","",60,0,6);
  TH1D * ETmissH = new TH1D("ETmissH","",40,0,200);
  TH1D * MtH = new TH1D("MtH","",40,0,200);
  TH1D * DZetaH = new TH1D("DZetaH","",60,-400,200);

  // histograms (dilepton selection + DZeta cut DZeta)
  TH1D * electronPtSelH  = new TH1D("electronPtSelH","",40,0,200);
  TH1D * electronEtaSelH = new TH1D("electronEtaSelH","",50,-2.5,2.5); 
  TH1D * muonPtSelH  = new TH1D("muonPtSelH","",40,0,200);
  TH1D * muonEtaSelH = new TH1D("muonEtaSelH","",50,-2.5,2.5); 

  TH1D * dileptonMassSelH = new TH1D("dileptonMassSelH","",40,0,200);
  TH1D * dileptonPtSelH = new TH1D("dileptonPtSelH","",40,0,200);
  TH1D * dileptonEtaSelH = new TH1D("dileptonEtaSelH","",100,-5,5);
  TH1D * dileptondRSelH = new TH1D("dileptondRSelH","",60,0,6);
  TH1D * ETmissSelH = new TH1D("ETmissSelH","",40,0,200);
  TH1D * MtSelH = new TH1D("MtSelH","",40,0,200);
  TH1D * DZetaSelH = new TH1D("DZetaSelH","",60,-400,200);

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
  int iCut=0;
    double CFCounter[CutNumb];
    double statUnc[CutNumb];
    int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; i++){
          CFCounter[i] = 0;
         iCFCounter[i] = 0;
         statUnc[i] =0;
        }

  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
 
  SetupHists(CutNumb); 
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
    //NE=1000;

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
     //numberOfEntries = 1000;
    
     std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
     TLorentzVector ElV, MuV, JetsV, METV;

     vector<TLorentzVector> JetsMV;
     vector<TLorentzVector>  ElMV;
     vector<TLorentzVector>  MuMV;
     vector<TLorentzVector>  LeptMV;

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
     
      analysisTree.GetEntry(iEntry);
      nEvents++;
      //JetsMV = new vector<TLorentzVector>;
      //ElMV = new vector<TLorentzVector>;
      //MuMV = new vector<TLorentzVector>;
      //LeptMV = new vector<TLorentzVector>;
      
       

      JetsMV.clear();
      ElMV.clear();
      MuMV.clear();


      Float_t MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      
      METV.SetPx(analysisTree.pfmet_ex);	      
      METV.SetPy(analysisTree.pfmet_ey);
 
      //Float_t BTAGTree = analysisTree.pfjet_btag[6];
      //MetV.SetPxPyPzE(analysisTree.pfmet_ex, analysisTree.pfmet_ey[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
 

     // ElV.SetPxPyPzE(analysisTree.electron_px[ij], analysisTree.ele_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
     // MuV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      for (unsigned int ij = 0; ij<analysisTree.pfjet_count; ++ij) {
//     if (analysisTree.pfjet_pt[ij]>30) { 
	     JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      JetsMV.push_back(JetsV);
  //     }
      } 




      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
    //  if (analysisTree.muon_pt[im]>10) { 
	      MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
       MuMV.push_back(MuV);

      //}
      }

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
  //    if ( analysisTree.electron_pt[ie]>10){ 
	      ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
       ElMV.push_back(ElV);
    //  }
      }


      float weight = 1;
       iCut = 0;
      
      Double_t EvWeight = 1.0;
      EvWeight *= weight ;
      
      FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
     // hnJets[->Fill(pfjet_count);

      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      // int indexW  = -1;
      // int indexNu = -1; 
      // int indexMu = -1;
      // int indexE  = -1;
      // int nGenMuons = 0;
      // int nGenElectrons = 0;
      // for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      // 	float pxGen = analysisTree.genparticles_px[igen];
      // 	float pyGen = analysisTree.genparticles_py[igen];
      // 	float pzGen = analysisTree.genparticles_pz[igen];
      // 	float etaGen = PtoEta(pxGen,pyGen,pzGen);
      // 	float ptGen  = PtoPt(pxGen,pyGen);

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
      // 	  indexW = igen;
      // 	if ((fabs(analysisTree.genparticles_pdgid[igen])==12 
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==14
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==16) 
      // 	    && analysisTree.genparticles_info[igen]== (1<<1) )
      // 	  indexNu = igen;

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==13) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexMu = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenMuons++;
      // 	  }
      // 	}
      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==11) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexE = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenElectrons++;
      // 	  }
      // 	}
      // }

      // trigger selection
 
      //selecTable.Fill(1,0, weight );      
      bool trigAccept = false;

      for (int i=0; i<kMaxhltriggerresults; ++i) {
	if ((i==5||i==6)&&analysisTree.hltriggerresults_second[i]==1) {
	  //	  std::cout << analysisTree.run_hltnames->at(i) << " : " << analysisTree.hltriggerresults_second[i] << std::endl;
	  trigAccept = true;
	}
      }

      if (!trigAccept) continue;
      //Trigger
      FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      CFCounter[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
      /////now clear the Mu.El.Jets again to fill them again after cleaning
      MuMV.clear();
      ElMV.clear();
     // JetsMV.clear();
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
	electrons.push_back(ie);
        ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
        ElMV.push_back(ElV);
	LeptMV.push_back(ElV);
        hel_miniISO[1]->Fill(analysisTree.electron_miniISO[ie],weight);
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
	MuV.SetPtEtaPhiM(analysisTree.muon_pt[im], analysisTree.muon_eta[im], analysisTree.muon_phi[im], muonMass);
        MuMV.push_back(MuV);
	LeptMV.push_back(MuV);
        hmu_miniISO[1]->Fill(analysisTree.muon_miniISO[im],weight);
      }


      ///sort leptons vector
      sort(LeptMV.begin(), LeptMV.end(),ComparePt); 
      
	//for (unsigned int j=0;j<LeptMV.size();++j) cout<<" j "<<j<<"  "<<LeptMV.at(j).Pt()<<endl;
	//cout<<""<<endl;
      ////////jets cleaning 
      Float_t DRmax=0.4;
      vector<int> jets; jets.clear();
      TLorentzVector leptonsV, muonJ, jetsLV;
      
      
      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      for (unsigned int il = 0; il<LeptMV.size(); ++il) {
      
	 for (unsigned int ij = 0; ij<JetsMV.size(); ++ij) {
        
		 if(fabs(JetsMV.at(ij).Eta())>etaJetCut) continue;
                 if(fabs(JetsMV.at(ij).Pt())<ptJetCut) continue;
      
       Float_t Dr= deltaR(LeptMV.at(il).Eta(), LeptMV.at(il).Phi(),JetsMV.at(ij).Eta(),JetsMV.at(ij).Phi());

     if (  Dr  < DRmax) {
	     
	     JetsMV.erase (JetsMV.begin()+ij);
    		 }	
		       
	 }
      }
      

      if (electrons.size()==0) continue;
      if (muons.size()==0) continue;
      //if (jets.size()<4) continue;

      // selecting muon and electron pair (OS or SS);
      float ptScalarSum = -1;
      float dRleptons = -1;
      int electronIndex = -1;
      int muonIndex = -1;
      

      if (SelectionSign == "1l" || SelectionSign == "1L"){
      	//cout<<" LeptMV.at(0).Pt() "<<LeptMV.at(0).Pt()<<endl;
	if (LeptMV.at(0).Pt()< ptMuonLowCut) continue;  /// need to change that actually
	//if (LeptMV.size()!=1) continue;
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
      }



      if (SelectionSign == "SS" || SelectionSign == "OS"){
      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	int eIndex = electrons[ie];
	    
	if ( analysisTree.electron_pt[eIndex]< ptElectronLowCut ) ElMV.erase(ElMV.begin()+ie);
	
			for (unsigned int im=0; im<muons.size(); ++im) {
	  int mIndex = muons[im];
	  float qProd = analysisTree.electron_charge[eIndex]*analysisTree.muon_charge[mIndex];
	
	  if ( analysisTree.muon_pt[mIndex]<ptMuonLowCut ) MuMV.erase(MuMV.begin()+im);

	  //// This is for the SS to be true
	  if (SelectionSign == "SS") {
	  if (oppositeSign || qProd<0) continue;
	  //if (!oppositeSign && qProd<0) continue;
	  
	  }
	  if (SelectionSign == "OS") {
	  if (!oppositeSign || qProd>0) continue;
	  //if (oppositeSign &&  qProd<0) continue;
	  }
           
	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;
          //CFCounter[iCut]+= weight;
          //iCFCounter[iCut]++;
          //iCut++;

		
           Float_t leadLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] > analysisTree.muon_pt[mIndex]	: analysisTree.muon_pt[mIndex];
           Float_t trailLep = analysisTree.electron_pt[eIndex] ? analysisTree.electron_pt[eIndex] < analysisTree.muon_pt[mIndex] : analysisTree.muon_pt[mIndex];

 
	   //cout<<"  Leading "<<leadLep<<" trailing "<<trailLep<<endl;
	   //cout<<"analysisTree.electron_pt[eIndex] "<<analysisTree.electron_pt[eIndex]<<"  nalysisTree.muon_pt[mIndex] "<<analysisTree.muon_pt[mIndex]<<endl;
           //cout<<""<<endl;

	  bool kinematics = 
            (analysisTree.electron_pt[eIndex]> ptElectronLowCut && analysisTree.muon_pt[mIndex]>ptMuonHighCut) || 
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
      	}/// dR cut

          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	
      }//OS or SS selection
	  //if (JetsMV.size()<3) continue;

	  /*for (unsigned int ij = 0; ij < JetsMV.size();ij++){
              if (JetsMV.at(ij).Pt()<30) continue;
	  }*/

 	  bool btagged= false;
	  for (unsigned int ib = 0; ib <analysisTree.pfjet_count;ib++){
            if (analysisTree.pfjet_btag[ib][6]  > bTag) btagged = true;
  		  //cout<<" pfjet_b "<<ib<<"  "<<analysisTree.pfjet_btag[ib][6]<<endl;
	  }
	  if (btagged) continue;

          // Jets
	  FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
          // pt Scalar
	  if (SelectionSign == "OS" || SelectionSign == "SS") {
    	  if (ptScalarSum<0  ) continue;
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
	  }
      // computations of kinematic variables

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
      float dPhi=-999; 
      float MT = -999;

	if (SelectionSign !="1l" && SelectionSign !="1L") {
	   	dPhi=dPhiFrom2P( dileptonLV.Px(), dileptonLV.Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      		MT=TMath::Sqrt(2*dileptonPt*ETmiss*(1-TMath::Cos(dPhi)));
	}

		if (SelectionSign =="1l" || SelectionSign =="1L")  {
		    	dPhi=dPhiFrom2P( LeptMV.at(0).Px(), LeptMV.at(0).Py(), analysisTree.pfmet_ex,  analysisTree.pfmet_ey );
      		MT = TMath::Sqrt(2*LeptMV.at(0).Pt()*ETmiss*(1-TMath::Cos(dPhi)));
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

      // topological cut
      //if (DZeta<dZetaCut) continue;
       if (dPhi<1) continue; 
       
          FillMainHists(iCut, EvWeight, ElMV, MuMV, JetsMV,METV);
      	  CFCounter[iCut]+= weight;
          iCFCounter[iCut]++;
          iCut++;
      
	  
     if (analysisTree.electron_count>0)  hel_miniISO[CutNumb-1]->Fill(analysisTree.electron_miniISO[0],EvWeight);
     if (analysisTree.muon_count>0)  hmu_miniISO[CutNumb-1]->Fill(analysisTree.muon_miniISO[0],EvWeight);

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


      //      std::cout << std::endl;
      
      selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    
    file_->Close();
    delete file_;
}
 
    cout << endl << "Finished event loop" << endl;
    for (int i=0;i<CutNumb;++i){
	         CFCounter[i] *= float(XSec*Lumi/inputEventsH->GetSum());
                 if (iCFCounter[i] <0.2) statUnc[i] =0;
                else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
        }

    //write out cutflow
    ofstream tfile;
   // TString outname = argv[argc-1];
    TString outname=argv[2];
    TString textfilename = "cutflow_"+outname+"_"+SelectionSign+".txt";
    tfile.open(textfilename);
    tfile << "########################################" << endl;
    tfile << "Cut efficiency numbers:" << endl;

        tfile << " Cut "<<"\t & \t"<<"#Evnts for "<<Lumi/1000<<" fb-1 & \t"<<" Uncertainty \t"<<" cnt\t"<<endl;
       for(int ci = 0; ci < CutNumb; ci++)
        {
                tfile << CutList[ci]<<"\t & \t"
                      << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
                CutFlow->SetBinContent(1+ci,CFCounter[ci]);
        }

    tfile.close();
    //ofstream tfile1;
    //TString textfile_Con = "CMG_cutflow_Con_Mu_"+outname+".txt";
    //tfile1.open(textfile_Con);
    //tfile1 << "########################################" << endl;
    //tfile1 << "RCS:" << endl;






  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");

/*
 for(int cj = 0; cj < CutNumb; cj++)
    {
        file->cd("");
        //outf->mkdir(CutList[cj]);
        //outf->cd(CutList[cj]);
        h0JetpT[cj]->Write();
        hnJet[cj]->Write();
        hnOver[cj]->Write();
        hnBJet[cj]->Write();
        hnEl[cj]->Write();
        hElpt[cj]->Write();
        hnMu[cj]->Write();
        hMupt[cj]->Write();
        hLepeta[cj]->Write();
        hMET[cj]->Write();
        hHT[cj]->Write();
        hST[cj]->Write();
        hToppT[cj]->Write();
        hnTop[cj]->Write();
        hWTagpT[cj]->Write();
        hWTagMass[cj]->Write();
        hnW[cj]->Write();
        hWmassTagpT[cj]->Write();
        hWmassTagMass[cj]->Write();
        hnWmass[cj]->Write();
        hdPhiMETLep[cj]->Write();
        hdPhiJMET[cj]->Write();

    }
*/

  file->Write();
  file->Close();
  
  delete file;
  
}



