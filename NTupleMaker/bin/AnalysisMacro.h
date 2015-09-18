
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

using namespace std;
const Float_t MuMass = 0.105658367;
const Float_t tauMass = 1.776;


const Float_t electronMass = 0;
const Float_t muonMass = 0.10565837;
const Float_t pionMass = 0.1396;
const  int CutN=15;


unsigned int RunMin = 9999999;
unsigned int RunMax = 0;
     
unsigned int tau_index=-1;
unsigned int mu_index=-1;
unsigned int el_index=-1;

bool isData = false;

Float_t XSec=-1;
Float_t xs,fact,fact2;
   //string CutList[10];
vector<string> CutList;

//TH1D * histWeights = new TH1D("histWeights","",1,-0.5,0.5);
TH1D * histWeights = new TH1D("histWeights","",1,0,0);
TH1D * histWeights2 = new TH1D("histWeights2","",1,0,0);

TH1D * hWeights [CutN];

TH1D *hHT[CutN];
TH1D *hHT2[CutN];
TH1D *hHT3[CutN];
TH1D *hHT4[CutN];
//TH1D *hST[CutN];
//TH1D *h0JetpT[CutN];
TH1D *hnJet[CutN];
TH1D *hnBJet[CutN];

TH1D *hInvMassMuTau[CutN];
TH1D *hInvMassMuEl[CutN];
TH1D *hInvMassTauTau[CutN];
TH1D *hInvMassElEl[CutN];
TH1D *hInvMassElTau[CutN];
TH1D *hInvMassMuMu[CutN];


TH1D *hnEl[CutN];
TH1D *hElpt[CutN];
TH1D *hEleta[CutN];
TH1D *hel_relISO[CutN];
TH1D *hel_relISOL[CutN];
TH1D *hel_miniISO[CutN];
TH1D *hel_miniISOL[CutN];

TH1D *hnLep[CutN];
TH1D *hLeppt[CutN];
TH1D *hLepeta[CutN];

TH1D *hnMu[CutN];
TH1D *hMupt[CutN];
TH1D *hMueta[CutN];

TH1D *hmu_relISO[CutN];
TH1D *hmu_relISOL[CutN];
TH1D *hmu_miniISO[CutN];
TH1D *hmu_miniISOL[CutN];




TH1D *hnTau[CutN];
TH1D *hTaupt[CutN];
TH1D *hTaueta[CutN];


TH1D *hMET[CutN];
//TH1D *hnOver[CutN];
TH1D *hdPhiMETLep[CutN];
TH1D *hdPhiJMET[CutN];

TH1D *hMT[CutN];
TH1D *hMTel[CutN];
TH1D *hMTmu[CutN];
TH1D *hMTtau[CutN];
TH1D *hDZeta[CutN];

TH1D *hMTmutau[CutN];
TH1D *hMTmuel[CutN];
TH1D *hMTeltau[CutN];
TH1D *hMTtautau[CutN];


TH1D *hdR_mutau[CutN];
TH1D *hdR_eltau[CutN];
TH1D *hdR_tautau[CutN];
TH1D *hdR_muel[CutN];


TH1D *hnpv[CutN];
TH1D *hnpu[CutN];
TH1D *hnrho[CutN];


TH2D *hmet_MT[CutN];
TH2D *hmet_MTel[CutN];
TH2D *hmet_MTmu[CutN];
TH2D *hmet_MTtau[CutN];

TH2D *hmet_dPhi[CutN];
TH2D *hmet_dPhiel[CutN];
TH2D *hmet_dPhimu[CutN];
TH2D *hmet_dPhitau[CutN];

TH2D *hMT_dPhi[CutN];
TH2D *hMT_dPhiel[CutN];
TH2D *hMT_dPhimu[CutN];
TH2D *hMT_dPhitau[CutN];
  
  TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutN,0.5,CutN+0.5);

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);

  TH1D * muonPtAllH = new TH1D("muonPtAllH","",40,0,200);
  TH1D * electronPtAllH = new TH1D("electronPtAllH","",40,0,200);
  TH1D * tauPtAllH = new TH1D("tauPtAllH","",40,0,200);

  // histograms (dilepton selection)
  TH1D * electronPtH  = new TH1D("electronPtH","",40,0,200);
  TH1D * electronEtaH = new TH1D("electronEtaH","",50,-2.5,2.5); 
  TH1D * muonPtH  = new TH1D("muonPtH","",40,0,200);
  TH1D * muonEtaH = new TH1D("muonEtaH","",50,-2.5,2.5); 
  TH1D * tauEtaAllH = new TH1D("tauEtaAllH","",50,-2.5,2.5); 
 
  TH1D * dileptonMassH = new TH1D("dileptonMassH","",40,0,200);
  TH1D * dileptonPtH = new TH1D("dileptonPtH","",40,0,200);
  TH1D * dileptonEtaH = new TH1D("dileptonEtaH","",100,-5,5);
  TH1D * dileptondRH = new TH1D("dileptondRH","",60,0,6);
  TH1D * ETmissH = new TH1D("ETmissH","",40,0,200);
  TH1D * MtH = new TH1D("MtH_2l","",40,0,200);
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
  TH1D * MtSelH = new TH1D("MtSelH_2l","",40,0,200);
  TH1D * DZetaSelH = new TH1D("DZetaSelH","",60,-400,200);

     TLorentzVector ElV, MuV, TauV, JetsV, METV;

     vector<TLorentzVector> JetsMV;
     vector<TLorentzVector>  ElMV;
     vector<TLorentzVector>  MuMV;
     vector<TLorentzVector>  TauMV;
     vector<TLorentzVector>  LeptMV;

bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


int binNumber(Float_t x, int nbins, Float_t * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

Float_t effBin(Float_t x, int nbins, Float_t * bins, Float_t * eff) {

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



bool electronMvaIdTight(Float_t eta, Float_t mva) {

  Float_t absEta = fabs(eta);

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


bool electronVetoTight(Float_t SuperClusterEta, Float_t eta, Float_t phi, Float_t full5, Float_t hOverE, Float_t d0, Float_t dZ, Float_t ooE, Float_t pfISO, Float_t nMissing, bool convVeto) {


  bool passed = false;


if (fabs(SuperClusterEta) <= 1.479 ){
	
	if (  fabs(eta)<0.013625 && fabs(phi)< 0.230374 && fabs(full5) < 0.011586 && hOverE < 0.181130 && fabs(d0) < 0.094095 && fabs(dZ) < 0.713070 && ooE < 0.295751 && pfISO < 0.158721 && nMissing <= 2 && convVeto) 
	       passed = true;
}

else if (      1.479  < fabs(SuperClusterEta)  && fabs(SuperClusterEta) < 2.5 ){
 
  	if ( fabs(eta)<0.011932 && fabs(phi)< 0.255450 && fabs(full5) < 0.031849 && hOverE < 0.223870 && fabs(d0) < 0.342293 && fabs(dZ) < 0.953461 && ooE < 0.155501 && pfISO < 0.177032 && nMissing <= 3 && convVeto) 

	       passed = true;
}
       else  
	       passed = false;

  return passed;
}









//string CutList[CutN];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNer){
    for(int cj = 0; cj < CutNer; cj++)
    {
        CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
        TString cutName=CutList[cj];
        TString nCut;
        nCut.Form("%d",cj);
	///generic variables
	//
	//Jets
        hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,400,0.0,4000.0);
        hHT[cj]->Sumw2();
        hHT2[cj] = new TH1D ("HT2_"+nCut,"HT2 "+cutName,400,0.0,4000.0);
        hHT2[cj]->Sumw2();
        hHT3[cj] = new TH1D ("HT3_"+nCut,"HT3 "+cutName,400,0.0,4000.0);
        hHT3[cj]->Sumw2();
        hHT4[cj] = new TH1D ("HT4_"+nCut,"HT4 "+cutName,400,0.0,4000.0);
        hHT4[cj]->Sumw2();


        //h0JetpT[cj] = new TH1D ("0JetpT_"+nCut,"0JetpT "+cutName,200,0.0,2000.0);
        //h0JetpT[cj]->Sumw2();
        hnJet[cj] = new TH1D ("nJet_"+nCut,"nJet "+cutName,20,0,20);
        hnJet[cj]->Sumw2();
        hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,20,0,20);
        hnBJet[cj]->Sumw2();

        hWeights[cj] = new TH1D ("hWeights_"+nCut,"hWeights "+cutName,20,0,0);
        hWeights[cj]->Sumw2();
	
        hInvMassMuTau[cj] = new TH1D ("hInvMassMuTau_"+nCut,"hInvMassMuTau "+cutName,80,0,160);
        hInvMassMuMu[cj] = new TH1D ("hInvMassMuMu_"+nCut,"hInvMassMuMu "+cutName,80,0,160);
        hInvMassElTau[cj] = new TH1D ("hInvMassElTau_"+nCut,"hInvMassElTau "+cutName,80,0,160);
        hInvMassElEl[cj] = new TH1D ("hInvMassElEl_"+nCut,"hInvMassElEl "+cutName,80,0,160);
        hInvMassTauTau[cj] = new TH1D ("hInvMassTauTau_"+nCut,"hInvMassTauTau "+cutName,80,0,160);
        
	//Leptons
	//
	//
	hnLep[cj] = new TH1D ("nLep_"+nCut,"nLep "+cutName,10,0,10);
        hnLep[cj]->Sumw2();
	hLeppt[cj] = new TH1D ("LeppT_"+nCut,"Lep pT "+cutName,100,0,1000);
        hLeppt[cj]->Sumw2();
        hLepeta[cj] = new TH1D ("Lepeta_"+nCut,"Lep eta "+cutName,100,-4,4);
        hLepeta[cj]->Sumw2();
	//hST[cj] = new TH1D ("ST_"+nCut,"ST "+cutName,400,0.0,4000.0);
        //hST[cj]->Sumw2();
        
	//Muons
	//
	//
	hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,0,10);
        hnMu[cj]->Sumw2();
        hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,100,0,1000);
        hMupt[cj]->Sumw2();
        hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,100,-4,4);
        hMueta[cj]->Sumw2();
        
	//Taus
	//
	//
        hnTau[cj] = new TH1D ("nTau_"+nCut,"nTau "+cutName,10,0,10);
        hnTau[cj]->Sumw2();
        hTaupt[cj] = new TH1D ("TaupT_"+nCut,"Tau pT "+cutName,100,0,1000);
        hTaupt[cj]->Sumw2();
        hTaueta[cj] = new TH1D ("Taueta_"+nCut,"Tau eta "+cutName,100,-4,4);
        hTaueta[cj]->Sumw2();
	
	//hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
        //Electrons
	//
	//
	hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,0,10);
        hnEl[cj]->Sumw2();
        hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,100,0,1000);
        hElpt[cj]->Sumw2();
	hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,100,-4,4);
        hEleta[cj]->Sumw2();
       
       
	hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,200.0,0.0,2000.0);
        hMET[cj]->Sumw2();
        //dPhi
        //
        //
       	hdPhiMETLep[cj] = new TH1D("dPhiMETLep_"+nCut,"dPhiMETLep "+cutName,64,0.0,3.2);
        hdPhiMETLep[cj]->Sumw2();
        hdPhiJMET[cj] = new TH1D("dPhiJMET_"+nCut,"dPhiJMET "+cutName,64,0.0,3.2);
        hdPhiJMET[cj]->Sumw2();

	//MT
	//
	//
        hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,40,0,200);
        hMT[cj]->Sumw2();
        hMTel[cj] = new TH1D ("MTel_"+nCut,"MTel "+cutName,40,0,200);
        hMTel[cj]->Sumw2();
        hMTmu[cj] = new TH1D ("MTmu_"+nCut,"MTmu "+cutName,40,0,200);
        hMTmu[cj]->Sumw2();

        hMTtau[cj] = new TH1D ("MTtau_"+nCut,"MTtau "+cutName,40,0,200);
        hMTtau[cj]->Sumw2();

        hMTmutau[cj] = new TH1D ("MTmutau_"+nCut,"MTmutau "+cutName,40,0,200);
        hMTmutau[cj]->Sumw2();
     
     	hMTmuel[cj] = new TH1D ("MTmuel_"+nCut,"MTmuel "+cutName,40,0,200);
        hMTmuel[cj]->Sumw2();
       
        hMTeltau[cj] = new TH1D ("MTeltau_"+nCut,"MTmeltau "+cutName,40,0,200);
        hMTeltau[cj]->Sumw2();
      
      	hMTtautau[cj] = new TH1D ("MTtautau_"+nCut,"MTtautau "+cutName,40,0,200);
        hMTtautau[cj]->Sumw2();

	/*       hDZeta[CutN]= new TH1D ("DZeta_"+nCut,"DZeta "+cutName,300,-400,200);
        hDZeta[cj]->Sumw2();
*/
        hel_miniISO[cj]= new TH1D ("elminiISO_"+nCut,"elminiISO "+cutName,50,0,5);;
        hel_miniISO[cj]->Sumw2();
        hel_miniISOL[cj]= new TH1D ("elminiISOL_"+nCut,"elminiISOL "+cutName,50,0,5);;
        hel_miniISOL[cj]->Sumw2();

        hel_relISO[cj]= new TH1D ("elrelISO_"+nCut,"elrelISO "+cutName,50,0,5);;
        hel_relISO[cj]->Sumw2();
	hel_relISOL[cj]= new TH1D ("elrelISOL_"+nCut,"murelISOL "+cutName,50,0,5);;
        hel_relISOL[cj]->Sumw2();
        
        
        hmu_miniISO[cj]= new TH1D ("muminiISO_"+nCut,"muminiISO "+cutName,50,0,5);;
        hmu_miniISO[cj]->Sumw2();
    	hmu_miniISOL[cj]= new TH1D ("muminiISOL_"+nCut,"muminiISOL "+cutName,50,0,5);;
        hmu_miniISOL[cj]->Sumw2();

	hmu_relISO[cj]= new TH1D ("murelISO_"+nCut,"murelISO "+cutName,50,0,5);;
        hmu_relISO[cj]->Sumw2();
	hmu_relISOL[cj]= new TH1D ("murelISOL_"+nCut,"murelISOL "+cutName,50,0,5);;
        hmu_relISOL[cj]->Sumw2();
 
 	hdR_eltau[cj]= new TH1D ("dR_eltau_"+nCut,"dR_eltau "+cutName,60,0,6);;
        hdR_eltau[cj]->Sumw2();
        
	hdR_mutau[cj]= new TH1D ("dR_mutau_"+nCut,"dR_mutau "+cutName,60,0,6);;
        hdR_mutau[cj]->Sumw2();

	hdR_tautau[cj]= new TH1D ("dR_tautau_"+nCut,"dR_tautau "+cutName,60,0,6);;
        hdR_tautau[cj]->Sumw2();
	
	hdR_muel[cj]= new TH1D ("dR_muel_"+nCut,"dR_muel "+cutName,60,0,6);;
        hdR_muel[cj]->Sumw2();

	hnpv[cj]= new TH1D ("npv_"+nCut,"npv "+cutName,50,0,50);;
        hnpv[cj]->Sumw2();
	hnpu[cj]= new TH1D ("npu_"+nCut,"npu "+cutName,50,0,50);;
        hnpu[cj]->Sumw2();
	hnrho[cj]= new TH1D ("nrho_"+nCut,"nrho "+cutName,50,0,50);;
        hnrho[cj]->Sumw2();
 
 
	hmet_dPhi[cj] = new TH2D ("met_dPhi_"+nCut,"met_dPhi "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
	hmet_dPhi[cj]->Sumw2();
        hmet_MT[cj] = new TH2D ("met_MT_"+nCut,"met_MT "+cutName,200.0,0.0,2000.0,40,0,200);
	hmet_MT[cj]->Sumw2();

	hmet_dPhiel[cj] = new TH2D ("met_dPhiel_"+nCut,"met_dPhiel "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
	hmet_dPhiel[cj]->Sumw2();
        hmet_MTel[cj] = new TH2D ("met_MTel_"+nCut,"met_MTel "+cutName,200.0,0.0,2000.0,40,0,200);
	hmet_MTel[cj]->Sumw2();
 

	hmet_dPhimu[cj] = new TH2D ("met_dPhimu_"+nCut,"met_dPhimu "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
	hmet_dPhimu[cj]->Sumw2();
	hmet_MTmu[cj] = new TH2D ("met_MTmu_"+nCut,"met_MTmu "+cutName,200.0,0.0,2000.0,40,0,200);
	hmet_MTmu[cj]->Sumw2();

	hmet_dPhitau[cj] = new TH2D ("met_dPhitau_"+nCut,"met_dPhitau "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
	hmet_dPhitau[cj]->Sumw2();
        hmet_MTtau[cj] = new TH2D ("met_MTtau_"+nCut,"met_MTtau "+cutName,200.0,0.0,2000.0,40,0,200);
	hmet_MTtau[cj]->Sumw2();


 
	hMT_dPhi[cj]= new TH2D ("MT_dPhi_"+nCut,"MT_dPhi "+cutName,40,0,200,64,0.0,3.2);
	hMT_dPhi[cj]->Sumw2();
        hMT_dPhiel[cj]= new TH2D ("MTel_dPhi_"+nCut,"MTel_dPhi "+cutName,40,0,200,64,0.0,3.2);
	hMT_dPhiel[cj]->Sumw2();
        hMT_dPhimu[cj]= new TH2D ("MTmu_dPhi_"+nCut,"MTmu_dPhi "+cutName,40,0,200,64,0.0,3.2);
	hMT_dPhimu[cj]->Sumw2();
        hMT_dPhitau[cj]= new TH2D ("MTtau_dPhi_"+nCut,"MTtau_dPhi "+cutName,40,0,200,64,0.0,3.2);
	hMT_dPhitau[cj]->Sumw2();

    }
}

void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV,vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_){}


void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel){}


void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV, vector<TLorentzVector>  TauV, vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_, string & Sel, int  mIndex, int eIndex, int  tIndex){	
	//void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  *JetsV){
if (JetsV.size()>0 )	hnJet[CutIndex]->Fill(JetsV.size(),EvWeight);
if (MuV.size() >0  )      hnMu[CutIndex]->Fill(MuV.size(),EvWeight);
if (TauV.size() >0 )      hnTau[CutIndex]->Fill(TauV.size(),EvWeight);
if (ElV.size() >0  )     hnEl[CutIndex]->Fill(ElV.size(),EvWeight);
        hnLep[CutIndex]->Fill(ElV.size()+MuV.size()+TauV.size(),EvWeight);
	hnpu[CutIndex]->Fill(tree_.numpileupinteractions,EvWeight);
	hnpv[CutIndex]->Fill(tree_.primvertex_count,EvWeight);
	hnrho[CutIndex]->Fill(tree_.rho,EvWeight);
	hWeights[CutIndex]->Fill(EvWeight);

	TLorentzVector muV ; muV.SetPtEtaPhiM(tree_.muon_pt[mIndex], tree_.muon_eta[mIndex], tree_.muon_phi[mIndex], muonMass);
	TLorentzVector tauV; tauV.SetPtEtaPhiM(tree_.tau_pt[tIndex], tree_.tau_eta[tIndex], tree_.tau_phi[tIndex], tauMass);
	TLorentzVector elV;  elV.SetPtEtaPhiM(tree_.electron_pt[eIndex], tree_.electron_eta[eIndex], tree_.electron_phi[eIndex], tauMass);

     	//     if (JetsV.size() > 0) h0JetpT[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);
  //  if(FillBJets){
  //      hnBJet[CutIndex]->Fill(Obj.nBJetGood,EvWeight);
  //  }
  //
  //cout << " Passing arguments befor" << JetsV.size()<<" muV "<<MuV.size()<<" tauV "<<TauV.size()<<" elV "<<ElV.size()<<" elI "<<eIndex<<" mI "<<mIndex<<"  tI "<<tIndex<<endl;
	if (Sel=="mutau" && MuV.size()>0 && TauV.size()>0 && mIndex >-1 && tIndex >-1 ){
	
//	     	Float_t Dr= deltaR(MuV.at(mIndex).Eta(), MuV.at(mIndex).Phi(),TauV.at(mIndex).Eta(),TauV.at(mIndex).Phi());
	     	//Float_t Dr= deltaR(MuV.at(0).Eta(), MuV.at(0).Phi(),TauV.at(0).Eta(),TauV.at(0).Phi());
		//Float_t Dr = MuV.at(mIndex).DeltaR(TauV.at(tIndex));
		Float_t Dr = muV.DeltaR(tauV);
//

        TLorentzVector DiL = muV  + tauV;
        Float_t dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
     	Float_t MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     	hMTmutau[CutIndex]->Fill(MT,EvWeight);
	hInvMassMuTau[CutIndex]->Fill(DiL.M(),EvWeight);
	hdR_mutau[CutIndex]->Fill(Dr,EvWeight);
	}

	if (Sel=="eltau" && ElV.size()>0 && TauV.size()>0  && mIndex>-1 && eIndex > -1){
        //Float_t Dr= deltaR(ElV.at(eIndex).Eta(), ElV.at(eIndex).Phi(),TauV.at(tIndex).Eta(),TauV.at(tIndex).Phi());
		Float_t Dr = elV.DeltaR(tauV);
	hdR_eltau[CutIndex]->Fill(Dr,EvWeight);
     
     	TLorentzVector DiL = elV + tauV;
        Float_t dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
     	Float_t MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     	hMTeltau[CutIndex]->Fill(MT,EvWeight);
	hInvMassElTau[CutIndex]->Fill(DiL.M(),EvWeight);
	
		}

	if (Sel=="muel" && MuV.size()>0 && ElV.size()>0  && mIndex>-1 && eIndex > -1 ){
        //Flaat_t Dr= deltaR(MuV.at(mIndex).Eta(), MuV.at(mIndex).Phi(),ElV.at(eIndex).Eta(),ElV.at(eIndex).Phi());
		Float_t Dr = ElV.at(0).DeltaR(MuV.at(0));
	hdR_muel[CutIndex]->Fill(Dr,EvWeight);
       TLorentzVector DiL = ElV.at(0) + MuV.at(0);
        Float_t dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
     	Float_t MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     	hMTmuel[CutIndex]->Fill(MT,EvWeight);
	hInvMassMuEl[CutIndex]->Fill(DiL.M(),EvWeight);
	}
//////////////// FIX THE INDEXES FOR TWO TAU CANDIDATES!!!!
	if (Sel=="tautau" &&  TauV.size()>1 && tIndex > -1){
		Float_t Dr = TauV.at(0).DeltaR(TauV.at(0));
        //Float_t Dr= deltaR(TauV.at(mIndex).Eta(), TauV.at(mIndex).Phi(),TauV.at(tIndex).Eta(),TauV.at(tIndex).Phi());
	hdR_mutau[CutIndex]->Fill(Dr,EvWeight);
        TLorentzVector DiL = TauV.at(0) + TauV.at(0);
        Float_t dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), MetV.Px(),  MetV.Py() );
     	Float_t MT = TMath::Sqrt(2*DiL.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     	hMTtautau[CutIndex]->Fill(MT,EvWeight);
	if (mIndex >-1 && tIndex >-1){
	TLorentzVector dimass = TauV.at(mIndex) + TauV.at(tIndex);
	hInvMassTauTau[CutIndex]->Fill(dimass.M(),EvWeight);
		}
	}

	if (ElV.size() > 0 && eIndex>-1)
    {
     hElpt[CutIndex]->Fill(elV.Pt(),EvWeight);
     hEleta[CutIndex]->Fill(elV.Eta(),EvWeight);

     Float_t dPhi=dPhiFrom2P( elV.Px(), elV.Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      		
     Float_t MT = TMath::Sqrt(2*elV.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     hMT[CutIndex]->Fill(MT,EvWeight);
     hMTel[CutIndex]->Fill(MT,EvWeight);

     hmet_dPhi[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);
     hmet_dPhiel[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);

     hmet_MT[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);
     hmet_MTel[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);

     hMT_dPhi[CutIndex]->Fill(MT,dPhi,EvWeight);
     hMT_dPhiel[CutIndex]->Fill(MT,dPhi,EvWeight);

	   
	hel_miniISOL[CutIndex]->Fill(tree_.electron_miniISO[eIndex],EvWeight);
	Float_t neutralIso = 
	  tree_.electron_neutralHadIso[eIndex] + 
	  tree_.electron_photonIso[eIndex] - 
	  0.5*tree_.electron_puIso[eIndex];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = tree_.electron_chargedHadIso[eIndex] + neutralIso;
	Float_t relIso = absIso/tree_.electron_pt[eIndex];
        hel_relISOL[CutIndex]->Fill(relIso,EvWeight);

    for (unsigned int ie=0;ie<ElV.size();ie++){

     hLeppt[CutIndex]->Fill(ElV.at(ie).Pt(),EvWeight);
     hLepeta[CutIndex]->Fill(ElV.at(ie).Eta(),EvWeight);

	hel_miniISO[CutIndex]->Fill(tree_.electron_miniISO[ie],EvWeight);
	Float_t neutralIso = 
	  tree_.electron_neutralHadIso[ie] + 
	  tree_.electron_photonIso[ie] - 
	  0.5*tree_.electron_puIso[ie];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = tree_.electron_chargedHadIso[ie] + neutralIso;
	Float_t relIso = absIso/tree_.electron_pt[ie];
        hel_relISO[CutIndex]->Fill(relIso,EvWeight);
	
   	}
    }




    if (MuV.size() > 0 && mIndex>-1)
    {
     TLorentzVector WBos = MetV + muV;
     hMupt[CutIndex]->Fill(muV.Pt(),EvWeight);
     hMueta[CutIndex]->Fill(muV.Eta(),EvWeight);
  

     Float_t dPhi=dPhiFrom2P( muV.Px(), muV.Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      
     Float_t MT = TMath::Sqrt(2*muV.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
     hMT[CutIndex]->Fill(MT,EvWeight);
     hMTmu[CutIndex]->Fill(MT,EvWeight);

     hmet_dPhi[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);
     hmet_dPhimu[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);

     hmet_MT[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);
     hmet_MTmu[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);

     hMT_dPhi[CutIndex]->Fill(MT,dPhi,EvWeight);
     hMT_dPhimu[CutIndex]->Fill(MT,dPhi,EvWeight);
 
	   
	hmu_miniISOL[CutIndex]->Fill(tree_.muon_miniISO[mIndex],EvWeight);
	Float_t neutralIso = 
	  tree_.muon_neutralHadIso[mIndex] + 
	  tree_.muon_photonIso[mIndex] - 
	  0.5*tree_.muon_puIso[mIndex];
	neutralIso = TMath::Max(Float_t(mIndex),neutralIso); 
	Float_t absIso = tree_.muon_chargedHadIso[mIndex] + neutralIso;
	Float_t relIso = absIso/tree_.muon_pt[mIndex];
        hmu_relISOL[CutIndex]->Fill(relIso,EvWeight);
     
     for (unsigned int im=0;im<MuV.size();im++){

     hLeppt[CutIndex]->Fill(MuV.at(im).Pt(),EvWeight);
     hLepeta[CutIndex]->Fill(MuV.at(im).Eta(),EvWeight);

	hmu_miniISO[CutIndex]->Fill(tree_.muon_miniISO[im],EvWeight);
   
	Float_t neutralIso = 
	  tree_.muon_neutralHadIso[im] + 
	  tree_.muon_photonIso[im] - 
	  0.5*tree_.muon_puIso[im];
	neutralIso = TMath::Max(Float_t(im),neutralIso); 
	Float_t absIso = tree_.muon_chargedHadIso[im] + neutralIso;
	Float_t relIso = absIso/tree_.muon_pt[im];
        hmu_relISO[CutIndex]->Fill(relIso,EvWeight);
   		}
    }
  
    if (TauV.size() > 0 && tIndex>-1)
    {
     hTaupt[CutIndex]->Fill(tauV.Pt(),EvWeight);
     hTaueta[CutIndex]->Fill(tauV.Eta(),EvWeight);

     for (unsigned int it=0;it<TauV.size();it++){
     hLeppt[CutIndex]->Fill(TauV.at(it).Pt(),EvWeight);
     hLepeta[CutIndex]->Fill(TauV.at(it).Eta(),EvWeight);

    }
    }

    Float_t sumpT=0;
    int bjets=0;
    
    if (JetsV.size()>0){
  
	    if (JetsV.size()>1)
	hHT2[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt(),EvWeight);
	    if (JetsV.size()>2)
	hHT3[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt(),EvWeight);
	    if (JetsV.size()>3)
	hHT4[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt()+JetsV.at(3).Pt(),EvWeight);

	    for (unsigned int ij=0;ij<JetsV.size();ij++){
         sumpT+=JetsV.at(ij).Pt();
         Float_t dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), MetV.Px(),  MetV.Py() );
    
     hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);
	 
     for (unsigned int ib = 0; ib <tree_.pfjet_count;ib++){
        if (float(tree_.pfjet_pt[ib]) == float(JetsV.at(ij).Pt())  &&  tree_.pfjet_btag[ib][8]  > 0.814) {
		bjets++;
	// cout<<tree_.pfjet_pt[ib] <<"  "<<JetsV.at(ij).Pt()<<"  "<< tree_.pfjet_btag[ib][8]  <<endl;
     	}
     }
      hnBJet[CutIndex]->Fill(bjets,EvWeight);
	    
		    }
     hHT[CutIndex]->Fill(sumpT,EvWeight);
    }
      hMET[CutIndex]->Fill(MetV.Pt(),EvWeight);
      
		                      //cout<<" pfjet_b "<<ib<<"  "<<analysisTree.pfjet_btag[ib][6]<<endl;
				      //          }
				      //                    if (btagged) continue;


     // hMT[CutIndex]->Fill(MT,EvWeight);
    //  hDZeta->Fill(DZeta,EvWeight);
}









