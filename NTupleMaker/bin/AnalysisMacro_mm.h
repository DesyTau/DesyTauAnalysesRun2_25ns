
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"


const float MuMass = 0.105658367;
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


double Module(double x, double y, double z){

  double module2=x*x+y*y+z*z;
  double module=TMath::Sqrt(module2);
  return module;

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

// histograms
 
//xsection
TH1D * hxsec = new TH1D("xsec","",1,0,0);
TH1D * histWeights = new TH1D("histWeights","",1,0,0);

TH1D * JPsiMassAllMuonsH = new TH1D("JPsiMassAllMuonsH","",200,2,4);
TH1D * JPsiMassAllMuonsDRCutH = new TH1D("JPsiMassAllMuonsDRCutH","",200,2,4);
TH1D * JPsiMassIdMuonsH = new TH1D("JPsiMassIdMuonsH","",200,2,4);
TH1D * JPsiMassIdMuonsDRCutH = new TH1D("JPsiMassIdMuonsDRCutH","",200,2,4);

TH1D * YpsilonMassAllMuonsH = new TH1D("YpsilonMassAllMuonsH","",400,8,12);
TH1D * YpsilonMassAllMuonsDRCutH = new TH1D("YpsilonMassAllMuonsDRCutH","",400,8,12);
TH1D * YpsilonMassIdMuonsH = new TH1D("YpsilonMassIdMuonsH","",400,8,12);
TH1D * YpsilonMassIdMuonsDRCutH = new TH1D("YpsilonMassIdMuonsDRCutH","",400,8,12);

TH1D * ZMassAllMuonsH = new TH1D("ZMassAllMuonsH","",60,60,120);
TH1D * ZMassAllMuonsDRCutH = new TH1D("ZMassAllMuonsDRCutH","",400,60,120);
TH1D * ZMassIdMuonsH = new TH1D("ZMassIdMuonsH","",60,60,120);
TH1D * ZMassIdMuonsDRCutH = new TH1D("ZMassIdMuonsDRCutH","",60,60,120);
TH1D * ZMassIsoMuonsH = new TH1D("ZMassIsoMuonsH","",60,60,120);
TH1D * ZMassIsoMuonsDRCutH = new TH1D("ZMassIsoMuonsDRCutH","",60,60,120);
TH2D * MassvsPt = new TH2D("MassvsPt","",100,0,200,200,0,400);


//cutsflow histrograms

const int Ncuts=20;
const int Nname=100;

TH1D* h_PrimaryVertexN[Ncuts];
TH1D* h_PrimaryVertexZ[Ncuts];
TH1D* h_muon_pt[Ncuts];
TH1D* h_muon_eta[Ncuts];
TH1D* h_muon_dxy[Ncuts];
TH1D* h_muon_dz[Ncuts];
TH1D* h_pt_muon1[Ncuts];
TH1D* h_pt_muon2[Ncuts];
TH1D* h_eta_muon1[Ncuts];
TH1D* h_eta_muon2[Ncuts];
TH1D* h_phi_muon1[Ncuts];
TH1D* h_phi_muon2[Ncuts];
TH1I* h_Njets[Ncuts];
TH1D* h_pfmet[Ncuts];
//TH1D* h_Ht[Ncuts];


TH1D * h_tot_m1m2 = new TH1D("h_tot_m1m2","",100,0.,200.);
TH1D * h_tot_pt1pt2 = new TH1D("h_tot_pt1pt2","",100,0.,500.);
TH1D * h_Ht = new TH1D("h_Ht","",100,0,1500);
//TH1I * h_primvert = new TH1I("h_primvert","",100,0,30);
TH1D * h_deltaR = new TH1D("h_deltaR","",100,-10,10);
TH1D * h_tot_mass = new TH1D("h_tot_mass","",100,0.,200.);
//TH1D * h_massTo12 = new TH1D("h_massTo12","",100,0.,200.);
//TH1D * h_massUpto12 = new TH1D("h_massUpto12","",100,0.,200.);
TH1D * h_coupleHPt = new TH1D("h_coupleHPt","",100,0.,200.);

//in the z mass window
TH1D * h_tot_m1m2_Z = new TH1D("h_tot_m1m2_Z","",100,50.,150.);
TH1D * h_tot_pt1pt2_Z = new TH1D("h_tot_pt1pt2_Z","",100,0.,500.);
//TH1D * h_Ht_Z = new TH1D("h_Ht_Z","",100,0,1500);
//TH1D * h_deltaR_Z = new TH1D("h_deltaR_Z","",100,-10,10);


void InitAllHisto(){

  char name[Nname];

  for(int i=0;i<Nname;i++) name[i]=0;//initalization

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_Njets_%i",i);
    h_Njets[i] = new TH1I(name,"",20,0.,20.);
  }
  /*
  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_Ht_%i",i);
    h_Ht[i] = new TH1D(name,"",100,0,1500);
  }
  */
  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_pfmet_%i",i);
    h_pfmet[i] = new TH1D(name,"",100,-100.,300.);
  }
  
  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_PrimaryVertexN_%i",i);
    h_PrimaryVertexN[i] = new TH1D(name,"",100,0,50);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_PrimaryVertexZ_%i",i);
    h_PrimaryVertexZ[i] = new TH1D(name,"",100,-20,20);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_muon_pt_%i",i);
    h_muon_pt[i] = new TH1D(name,"",750,0,1500);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_muon_eta_%i",i);
    h_muon_eta[i] = new TH1D(name,"",100,-10,10);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_muon_dxy_%i",i);
    h_muon_dxy[i] = new TH1D(name,"",100,-20,20);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_muon_dz_%i",i);
    h_muon_dz[i] = new TH1D(name,"",100,-20,20);
  }

  //histograms of 2 muons

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_pt_muon1_%i",i);
    h_pt_muon1[i] = new TH1D(name,"",750,0,1500);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_pt_muon2_%i",i);
    h_pt_muon2[i] = new TH1D(name,"",750,0,1500);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_eta_muon1_%i",i);
    h_eta_muon1[i] = new TH1D(name,"",100,-10,10);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_eta_muon2_%i",i);
    h_eta_muon2[i] = new TH1D(name,"",100,-10,10);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_phi_muon1_%i",i);
    h_phi_muon1[i] = new TH1D(name,"",100,-10,10);
  }

  for(int i=0;i<Ncuts;i++){
    sprintf(name,"h_phi_muon2_%i",i);
    h_phi_muon2[i] = new TH1D(name,"",100,-10,10);
  }

}



void FillAllGeneralHisto(Int_t CutsIndex, float Weight, Int_t Nmuon, float* Muon_pt, float* Muon_eta, float* Muon_dxy, float* Muon_dz){
  //void FillAllHisto(Int_t CutsIndex, float Weight, UInt_t PrimVertex){
  //h_Ht[CutsIndex]->Fill(Ht,Weight);
  for(Int_t i=0;i<Nmuon;i++){ 
    h_muon_pt[CutsIndex]->Fill(*(Muon_pt+i),Weight); 
    h_muon_eta[CutsIndex]->Fill(*(Muon_eta+i),Weight); 
    h_muon_dxy[CutsIndex]->Fill(*(Muon_dxy+i),Weight); 
    h_muon_dz[CutsIndex]->Fill(*(Muon_dz+i),Weight); 
  } 
}


void FillAllHisto(Int_t CutsIndex, float Weight, Float_t PrimVertexN, Float_t PrimVertexZ, Float_t Njets, Float_t pfmet){
  h_PrimaryVertexN[CutsIndex]->Fill(PrimVertexN,Weight);
  h_PrimaryVertexZ[CutsIndex]->Fill(PrimVertexZ,Weight);
  h_Njets[CutsIndex]->Fill(Njets,Weight);
  h_pfmet[CutsIndex]->Fill(pfmet,Weight);
}


void FillHistosMuon1(Int_t CutsIndex, float Weight, Float_t pt, Float_t eta, Float_t phi){
  h_pt_muon1[CutsIndex]->Fill(pt,Weight);
  h_eta_muon1[CutsIndex]->Fill(eta,Weight);
  h_phi_muon1[CutsIndex]->Fill(phi,Weight);
}


void FillHistosMuon2(Int_t CutsIndex, float Weight, Float_t pt, Float_t eta, Float_t phi){
  h_pt_muon2[CutsIndex]->Fill(pt,Weight);
  h_eta_muon2[CutsIndex]->Fill(eta,Weight);
  h_phi_muon2[CutsIndex]->Fill(phi,Weight);
}


void WriteAllHisto(){

  char name[Nname];

  for(int i=0;i<Nname;i++) name[i]=0;//initalization

  for(int i=0;i<Ncuts;i++){
    if(h_Njets[i]->GetEntries()>0) h_Njets[i]->Write();
  }

    for(int i=0;i<Ncuts;i++){
    if(h_pfmet[i]->GetEntries()>0) h_pfmet[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_PrimaryVertexN[i]->GetEntries()>0) h_PrimaryVertexN[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_PrimaryVertexZ[i]->GetEntries()>0) h_PrimaryVertexZ[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_muon_pt[i]->GetEntries()>0) h_muon_pt[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_muon_eta[i]->GetEntries()>0) h_muon_eta[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_muon_dxy[i]->GetEntries()>0) h_muon_dxy[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_muon_dz[i]->GetEntries()>0) h_muon_dz[i]->Write();
  }

  //histograms of 2 muons

  for(int i=0;i<Ncuts;i++){
    if(h_pt_muon1[i]->GetEntries()>0) h_pt_muon1[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_pt_muon2[i]->GetEntries()>0) h_pt_muon2[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_eta_muon1[i]->GetEntries()>0) h_eta_muon1[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_eta_muon2[i]->GetEntries()>0) h_eta_muon2[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_phi_muon1[i]->GetEntries()>0) h_phi_muon1[i]->Write();
  }

  for(int i=0;i<Ncuts;i++){
    if(h_phi_muon2[i]->GetEntries()>0) h_phi_muon2[i]->Write();
  }

}



struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject;


