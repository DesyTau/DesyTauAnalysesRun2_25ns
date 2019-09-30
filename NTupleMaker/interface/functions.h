#ifndef NTupleMakerFunctions_h
#define NTupleMakerFunctions_h

#include "TMath.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

const double MuMass = 0.105658367;
const double tauMass = 1.776;
const double electronMass = 0;
const double muonMass = 0.10565837;
const double pionMass = 0.1396;
using namespace std;



double mctcorr(const double v1[4],const double v2[4]
                         ,const double vds[4],const double ptm[2]
                         ,const double ecm=14000.0,const double mxlo=0.0);
double mct(const double v1[4],const double v2[4]);
double mt2(const double v1[4],const double v2[4]
                     ,const double vds[4],const double ptm[2]
                     ,const double ecm=14000.0,const double mxlo=0.0);
double mctminmt2(const double mctsqr,const double m1sqr,
                           const double m2sqr,const double chisqr=0.0);
double mt2neg(const double v1[4],const double v2[4]
                        ,const double ptm[2],const double mxlo=0.0);
double mcy(const double v1[4],const double v2[4]
                     ,const double vds[4],const double ptm[2]);
double mcx(const double v1[4],const double v2[4]
                     ,const double vds[4],const double ptm[2]);

// private:

  double m_mctecm, m_mctehat, m_pb;








vector <int > runvect_;
vector <int > lumivect_;
vector <int > eventvect_;
struct Event {
    std::string name;
    int run;
    int lumi;
    int eventrn;
};

std::vector<std::string>
split(const std::string &s, char delim = ':')
{

    std::vector<std::string> elems;

    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
            elems.push_back(item);
    }
    return elems;
}
  
bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }

int binNumber(float x, int nbins, float * bins) {

  int binN = nbins-1;

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

double mT (TLorentzVector v1, TLorentzVector Metv){

  double dPhimT=dPhiFrom2P( v1.Px(), v1.Py(), Metv.Px(),  Metv.Py() );

  double MT = TMath::Sqrt(2*v1.Pt()*Metv.Pt()*(1-TMath::Cos(dPhimT)));
  return MT;
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


bool electronVetoTight(double SuperClusterEta, double eta, double phi, double full5, double hOverE, double d0, double dZ, double ooE, double pfISO, double nMissing, bool convVeto) {


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


void ComputeMetFromHadRecoil(float Hparal,
			     float Hperp,
			     float genVPx, 
			     float genVPy,
			     float visVPx,
			     float visVPy,
			     float & metX,
			     float & metY) {
  
  float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
  float unitX = genVPx/genVPt;
  float unitY = genVPy/genVPt;

  float unitPhi = TMath::ATan2(unitY,unitX);
  float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
  float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

  float det = unitX*unitPerpY - unitY*unitPerpX;
  float Hx = (Hparal*unitPerpY - Hperp*unitY)/det;
  float Hy = (Hperp*unitX - Hparal*unitPerpX)/det;

  metX = -Hx - visVPx;
  metY = -Hy - visVPy;

}

void ComputeHadRecoilFromMet(float metX,
			     float metY,
			     float genVPx, 
			     float genVPy,
			     float visVPx,
			     float visVPy,
			     float & Hparal,
			     float & Hperp) {

  float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
  float unitX = genVPx/genVPt;
  float unitY = genVPy/genVPt;

  float unitPhi = TMath::ATan2(unitY,unitX);
  float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
  float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

  float Hx = -metX - visVPx;
  float Hy = -metY - visVPy;

  Hparal = Hx*unitX + Hy*unitY;
  Hperp = Hx*unitPerpX + Hy*unitPerpY;

}


struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject, myobjectX;


namespace genTools{
  float topPtWeight(float pt1, float pt2) {
  
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
  
    float a = 0.0615;    // Run2 a parameter
    float b = -0.0005;   // Run2 b parameter
  
    float w1 = TMath::Exp(a+b*pt1);
    float w2 = TMath::Exp(a+b*pt2);
  
    return TMath::Sqrt(w1*w2);  
  }
  
  float topPtWeight_Run1(float pt1, float pt2) {
  
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
  
    float a = 0.156;    // Run1 a parameter
    float b = -0.00137;  // Run1 b parameter
  
    float w1 = TMath::Exp(a+b*pt1);
    float w2 = TMath::Exp(a+b*pt2);
  
    return TMath::Sqrt(w1*w2);  
  }
  
  
  float topPtWeight(const AC1B& analysisTree, unsigned int LHC_runEra=2){
	// argument LHC_runEra allows to switch between Run1 and Run2 weights. Default is Run2.
	// for Run1 weights use "LHC_runEra =1", for Run2 weights use "LHC_runEra =2". 
  
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
  
    if(topPt > 0. && antitopPt > 0.){
	  if (LHC_runEra ==2)
      	return topPtWeight(topPt, antitopPt);
	  else if (LHC_runEra ==1)
		return topPtWeight_Run1(topPt, antitopPt);
	}
  
    return 1.;
  };
  
  // these are functionally the same as above, logic rewritten
  float calculate_topPtWeight(float pt1, float pt2, float a, float b) {
    if (pt1 > 400) pt1 = 400; // approximation works only up to 400 GeV
    if (pt2 > 400) pt2 = 400;
    float w1 = TMath::Exp(a + b*pt1);
    float w2 = TMath::Exp(a + b*pt2);
    
    return TMath::Sqrt(w1*w2);  
  }

  float return_topPtWeight(const AC1B &analysisTree, float a, float b){
    float topPt = -1;
    float antitopPt = -1;
    
    for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) { 
      if (analysisTree.genparticles_pdgid[igen] == 6)
        topPt = TMath::Sqrt(analysisTree.genparticles_px[igen] * analysisTree.genparticles_px[igen] + analysisTree.genparticles_py[igen] * analysisTree.genparticles_py[igen]);
      if (analysisTree.genparticles_pdgid[igen] == -6)
        antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen] * analysisTree.genparticles_px[igen] + analysisTree.genparticles_py[igen] * analysisTree.genparticles_py[igen]);
    }
    
    if(topPt > 0. && antitopPt > 0.)
      return calculate_topPtWeight(topPt, antitopPt, a, b);
    else 
      return 1.;
  };
  
  
  TLorentzVector genZ(const AC1B& analysisTree){
    TLorentzVector genZ; genZ.SetXYZM(0,0,0,91.2);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      if (analysisTree.genparticles_pdgid[igen]==23||analysisTree.genparticles_pdgid[igen]==22) {
	if (analysisTree.genparticles_fromHardProcess[igen])
	  genZ.SetXYZT(analysisTree.genparticles_px[igen],
		       analysisTree.genparticles_py[igen],
		       analysisTree.genparticles_pz[igen],
		       analysisTree.genparticles_e[igen]);
      }
    }
    return genZ;
  }

  TLorentzVector genL(const AC1B& analysisTree){
    TLorentzVector genL; genL.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
     
    bool isMuon = 0;
    bool isElectron = 0;
    bool isChargedLepton = 0;
    //bool isNeutrino = 0;
    bool fromHardProcessFinalState = 0;
    //bool isDirectHardProcessTauDecayProduct = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
      isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
      isChargedLepton = isMuon || isElectron;
      //isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||	fabs(analysisTree.genparticles_pdgid[igen])==14||	fabs(analysisTree.genparticles_pdgid[igen])==16;
      fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen] && analysisTree.genparticles_status[igen]==1;
      //isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

      /*if((fromHardProcessFinalState && isChargedLepton) || (isDirectHardProcessTauDecayProduct && (!isNeutrino)))
	genL += genPart;
    }

    return genL;*/

      if(fromHardProcessFinalState && isChargedLepton)
	genL += genPart;
    }
    
    bool fromHardProcess = 0;
    bool isLastCopy = 0;
    for (unsigned int itau=0; itau<analysisTree.gentau_count; ++itau) {
      genPart.SetXYZT(analysisTree.gentau_visible_px[itau],
		      analysisTree.gentau_visible_py[itau],
		      analysisTree.gentau_visible_pz[itau],
		      analysisTree.gentau_visible_e[itau]);
      
      fromHardProcess = analysisTree.gentau_fromHardProcess[itau];
      isLastCopy = analysisTree.gentau_isLastCopy[itau];

      if (fromHardProcess && isLastCopy)
	genL += genPart;
    }
      
    return genL;
      
  }

  TLorentzVector genNu(const AC1B& analysisTree){
    TLorentzVector genNu; genNu.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    bool isNeutrino = 0;
    bool isPrompt = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	fabs(analysisTree.genparticles_pdgid[igen])==14||
	fabs(analysisTree.genparticles_pdgid[igen])==16;
      isPrompt = analysisTree.genparticles_isPrompt[igen]||
	analysisTree.genparticles_isPromptTauDecayProduct[igen];
      
      if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
	if (isNeutrino) 
	  genNu += genPart;
      }
    }

    return genNu;
  }

  TLorentzVector genV(const AC1B& analysisTree){
    TLorentzVector genV; genV.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    bool isMuon = 0;
    bool isElectron = 0;
    bool isChargedLepton = 0;
    bool isNeutrino = 0;
    bool fromHardProcessFinalState = 0;
    bool isDirectHardProcessTauDecayProduct = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
      isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
      isChargedLepton = isMuon || isElectron;
      isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	           fabs(analysisTree.genparticles_pdgid[igen])==14||
	           fabs(analysisTree.genparticles_pdgid[igen])==16;
      fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen] && analysisTree.genparticles_status[igen]==1;
      isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

      if((fromHardProcessFinalState && (isChargedLepton || isNeutrino)) || isDirectHardProcessTauDecayProduct){
	genV += genPart;}
    }
    
    if (genV.Pt()<0.1) genV.SetXYZM(0.1,0.1,0.,0.);
    
    return genV;
  }

  int nJetsHad(const AC1B& analysisTree){
    int njetshad = 0;
    bool isChargedLepton = 0;
    bool fromHardProcess = 0;
    
    float genEta = 0.;
    float genPhi = 0.;
    
    for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

      if (analysisTree.pfjet_pt[jet]<=30.) continue;
      float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
      if (absJetEta >= 4.7) continue;
      
      // jetId
      float energy = analysisTree.pfjet_e[jet];
      energy *= analysisTree.pfjet_energycorr[jet];
      float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
      float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
      float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
      float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
      float muf = analysisTree.pfjet_muonenergy[jet]/energy;
      float chm = analysisTree.pfjet_chargedmulti[jet];
      float nm = analysisTree.pfjet_neutralmulti[jet];
      float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
      //bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      bool isPFJetId = false;
      if (absJetEta<=3.0)
	isPFJetId = (nhf < 0.99 && phf < 0.99 && npr > 1) && (absJetEta>2.4 || (chf>0 && chm > 0 && elf < 0.99));
      else
	isPFJetId = phf < 0.9 && nm > 10;
      //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99 && muf < 0.8) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      
      if (!isPFJetId) continue;

      int overlap = 0;
      
      for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {	
	isChargedLepton = ( fabs(analysisTree.genparticles_pdgid[igen])==11 ||
			    fabs(analysisTree.genparticles_pdgid[igen])==13);
	fromHardProcess = analysisTree.genparticles_fromHardProcess[igen];

	if (!isChargedLepton) continue;
	if (!fromHardProcess) continue;

	genEta = PtoEta( analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen]); 
	genPhi = PtoPhi( analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen]);
	
	float dR = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			  genEta, genPhi);
	if (dR>0.5) continue;
	
	overlap = 1;
	break;
      }

      if (overlap == 1) continue;

      for (unsigned int itau=0; itau<analysisTree.gentau_count; ++itau) {
	fromHardProcess = analysisTree.gentau_fromHardProcess[itau];
	
      	if (!fromHardProcess) continue;

	genEta = PtoEta( analysisTree.gentau_visible_px[itau], analysisTree.gentau_visible_py[itau], analysisTree.gentau_visible_pz[itau]);
	genPhi = PtoPhi( analysisTree.gentau_visible_px[itau], analysisTree.gentau_visible_py[itau]);
			 	
	float dR = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			  genEta, genPhi);
	if (dR>0.5) continue;
	
	overlap = 1;
	break;
      }

      if (overlap == 1) continue;
      
      njetshad++;
    }

    return njetshad;
  } 
  
  enum RecoilCorrectionsMethod{QuantileRemap=1, MeanResolution};

  int RecoilCorrections( RecoilCorrector& corr, int method,
			 float met, float metphi,
			 float vx, float vy,
			 float lx, float ly,
			 int njets,
			 float& metcorr, float& metphicorr ){
    float metx = met*TMath::Cos(metphi);
    float mety = met*TMath::Sin(metphi);
    float metcorrx = metx;
    float metcorry = mety;
    if (method == 1)
      corr.Correct(metx, mety, vx, vy, lx, ly, njets, metcorrx, metcorry);
    else if(method == 2)
      corr.CorrectByMeanResolution(metx, mety, vx, vy, lx, ly, njets, metcorrx, metcorry);

    metcorr = TMath::Sqrt(metcorrx * metcorrx + metcorry * metcorry);
    metphicorr = TMath::ATan2(metcorry, metcorrx);
            
    return method;
  }
};

namespace utils{
  enum channel{UNKNOWN = 0, ETAU, MUTAU, TAUTAU, EE, MUMU, EMU};
}

namespace calc{
  float mt( float lpt, float lphi, float met, float metphi){
    return sqrt(2*lpt*met*(1.-TMath::Cos(lphi-metphi)));
  }

  float pzetavis( float zx, float zy, float visx, float visy){
    return zx*visx+zy*visy;
  }

  float pzetavis(const TLorentzVector& v1, const TLorentzVector& v2){
    float v1ux = v1.Px()/v1.Pt();
    float v1uy = v1.Py()/v1.Pt();

    float v2ux = v2.Px()/v2.Pt();
    float v2uy = v2.Py()/v2.Pt();    
    
    float zx = v1ux + v2ux;
    float zy = v1uy + v2uy;
      
    float modz = TMath::Sqrt(zx*zx+zy*zy);
    zx = zx/modz;
    zy = zy/modz;
    
    return pzetavis( zx, zy, v1.Px()+v2.Px(), v1.Py()+v2.Py());
  }

  float pzetamiss( float zx, float zy, float met, float metphi){
    return zx*met*TMath::Cos(metphi)+zy*met*TMath::Sin(metphi);
  }

  float pzetamiss(const TLorentzVector& v1, const TLorentzVector& v2, const TLorentzVector& met){
    float v1ux = v1.Px()/v1.Pt();
    float v1uy = v1.Py()/v1.Pt();

    float v2ux = v2.Px()/v2.Pt();
    float v2uy = v2.Py()/v2.Pt();    
    
    float zx = v1ux + v2ux;
    float zy = v1uy + v2uy;
      
    float modz = TMath::Sqrt(zx*zx+zy*zy);
    zx = zx/modz;
    zy = zy/modz;
    
    return pzetamiss( zx, zy, met.Pt(), met.Phi()); 
  }

  std::shared_ptr<SVfitStandaloneAlgorithm> svFit(const TLorentzVector& v1, int dm_1,
						  const TLorentzVector& v2, int dm_2,
						  const TLorentzVector& met, const TMatrixD& covMET,
						  utils::channel ch, TFile* visPtRes){
    
    using namespace utils;
    
    if(ch == UNKNOWN){
      return 0;
    }

    svFitStandalone::kDecayType t1_decay = svFitStandalone::kTauToHadDecay;
    svFitStandalone::kDecayType t2_decay = svFitStandalone::kTauToHadDecay;

    if (ch == ETAU || ch == EE || ch == EMU)
      t1_decay = svFitStandalone::kTauToElecDecay;
    if (ch == EE)
      t2_decay = svFitStandalone::kTauToElecDecay;

    if (ch == MUTAU || ch == MUMU)
      t1_decay = svFitStandalone::kTauToMuDecay;
    if (ch == EMU || ch == MUMU)
      t2_decay = svFitStandalone::kTauToMuDecay;
    
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton( t1_decay, v1.Pt(), v1.Eta(), v1.Phi(), v1.M(), dm_1));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton( t2_decay, v2.Pt(), v2.Eta(), v2.Phi(), v2.M(), dm_2));    

    std::shared_ptr<SVfitStandaloneAlgorithm> algo = std::make_shared<SVfitStandaloneAlgorithm>(measuredTauLeptons, met.Px(), met.Py(), covMET, 0);
    algo->addLogM(false);
    algo->shiftVisPt(true,visPtRes);
    algo->integrateMarkovChain();

    return algo;
  }
}


double MCT(double Pt1, double Phi1, double Pt2, double Phi2){

	double mct = 0;
	mct = sqrt( 2*Pt1*Pt2*(1+cos(Phi1-Phi2)) );
	return mct;
}


double mct(const double v1[4],const double v2[4])
{
  double et1 = sqrt(fmax(v1[0]*v1[0]-v1[3]*v1[3],0.0));
  double et2 = sqrt(fmax(v2[0]*v2[0]-v2[3]*v2[3],0.0));
  return sqrt(fmax(pow(et1+et2,2)-pow(v1[1]-v2[1],2)
                      -pow(v1[2]-v2[2],2),0.0));
}

double mctcorr(const double v1[4],const double v2[4]
		       ,const double vds[4],const double ptm[2]
		       ,const double ecm,const double mxlo)
{
// v1, v2 are the 4-vectors of the two aggregated decay products
// vds is the 4-vector of the downstream objects (excluding v1 and v2)
// ptm is the pTmiss 2-vector {pxmiss,pymiss}
// ecm is the centre-of-mass energy (D=14TeV)
// mxlo is a lower bound on the mass of each invisible decay product (D=0)

  double ptus[2] = {-v1[1]-v2[1]-vds[1]-ptm[0],-v1[2]-v2[2]-vds[2]-ptm[1]};
  m_pb = sqrt(pow(ptus[0],2)+pow(ptus[1],2));

  if (m_pb==0) {

    return mct(v1,v2);    

  } else {

// Transform to new basis in transverse plane
// ISR goes in +ve x direction

    double vb1[4] = {v1[0],(v1[1]*ptus[0]+v1[2]*ptus[1])/m_pb
	      ,(v1[1]*ptus[1]-v1[2]*ptus[0])/m_pb,v1[3]};
    double vb2[4] = {v2[0],(v2[1]*ptus[0]+v2[2]*ptus[1])/m_pb
	      ,(v2[1]*ptus[1]-v2[2]*ptus[0])/m_pb,v2[3]};
    double vey1 = sqrt(fmax(pow(vb1[0],2)-pow(vb1[1],2)-pow(vb1[3],2),0.0));
    double vey2 = sqrt(fmax(pow(vb2[0],2)-pow(vb2[1],2)-pow(vb2[3],2),0.0));
    double ax = vb1[1]*vey2+vb2[1]*vey1;

// Boost v1 and v2 with Ecm
// v1 and v2 were boosted in -ve x direction,
// so to correct we boost v1 and v2 in +ve x direction 

    double beta = m_pb/ecm;
    double gamma = 1.0/sqrt(1.0-beta*beta);
    double vb1p[4] = {gamma*(vb1[0]+vb1[1]*beta),gamma*(vb1[1]+vb1[0]*beta)
		      ,vb1[2],vb1[3]};
    double vb2p[4] = {gamma*(vb2[0]+vb2[1]*beta),gamma*(vb2[1]+vb2[0]*beta)
		      ,vb2[2],vb2[3]};
    double vey1p = sqrt(fmax(pow(vb1p[0],2)-pow(vb1p[1],2)
			     -pow(vb1p[3],2),0.0));
    double vey2p = sqrt(fmax(pow(vb2p[0],2)-pow(vb2p[1],2)
			     -pow(vb2p[3],2),0.0));
    double axecm = vb1p[1]*vey2p+vb2p[1]*vey1p;
    m_mctecm = mct(vb1p,vb2p);

// Boost v1 and v2 with e0
// v1 and v2 were boosted in -ve x direction,
// so to correct we boost v1 and v2 in +ve x direction  

    beta = m_pb/(v1[0]+v2[0]+vds[0]
	       +sqrt(pow(ptm[0],2)+pow(ptm[1],2)+4.0*pow(mxlo,2)));
    gamma = 1.0/sqrt(1.0-beta*beta);
    double vb1pp[4] = {gamma*(vb1[0]+vb1[1]*beta),gamma*(vb1[1]+vb1[0]*beta)
		      ,vb1[2],vb1[3]};
    double vb2pp[4] = {gamma*(vb2[0]+vb2[1]*beta),gamma*(vb2[1]+vb2[0]*beta)
		      ,vb2[2],vb2[3]};
    double vey1pp = sqrt(fmax(pow(vb1pp[0],2)-pow(vb1pp[1],2)
			      -pow(vb1pp[3],2),0.0));
    double vey2pp = sqrt(fmax(pow(vb2pp[0],2)-pow(vb2pp[1],2)
			      -pow(vb2pp[3],2),0.0));
    double axehat = vb1pp[1]*vey2pp+vb2pp[1]*vey1pp;
    m_mctehat = mct(vb1pp,vb2pp);

    if ( (ax>=0) || (axecm>=0) ) {
      return m_mctecm;
    } else {
      if (axehat<0) {
	return m_mctehat;
      } else {    
	return sqrt(fmax(pow(vey1+vey2,2)-pow(vb1[2]-vb2[2],2),0.0));
      }
    }
  }
}


double mt2(const double v1[4],const double v2[4]
		       ,const double vds[4],const double ptm[2]
		       ,const double ecm,const double mxlo)
{
  // Boost corrected analytical mt2
  double m1sqr = pow(v1[0],2)-pow(v1[1],2)-pow(v1[2],2)-pow(v1[3],2);
  double m2sqr = pow(v2[0],2)-pow(v2[1],2)-pow(v2[2],2)-pow(v2[3],2);
  double chisqr = pow(mxlo,2);
  double m1 = sqrt(fmax(m1sqr,0.0));
  double m2 = sqrt(fmax(m2sqr,0.0));
  double chi = mxlo;

  double qtest[3] = {0.0,ptm[0]-(chi/m1)*v1[1],ptm[1]-(chi/m1)*v1[2]};
  qtest[0] = sqrt(chisqr+pow(qtest[1],2)+pow(qtest[2],2));
  double ptest[3] = {0.0,ptm[0]-(chi/m2)*v2[1],ptm[1]-(chi/m2)*v2[2]};
  ptest[0] = sqrt(chisqr+pow(ptest[1],2)+pow(ptest[2],2));

  double et1 = sqrt(fmax(v1[0]*v1[0]-v1[3]*v1[3],0.0));
  double et2 = sqrt(fmax(v2[0]*v2[0]-v2[3]*v2[3],0.0));

  double bq[3] = {et2+qtest[0],v2[1]+qtest[1],v2[2]+qtest[2]};
  double ap[3] = {et1+ptest[0],v1[1]+ptest[1],v1[2]+ptest[2]};

  // Use unbalanced solutions
  // This is an approximation as we haven't boost-corrected bq and ap

  if (pow(m1+chi,2)>=pow(bq[0],2)-pow(bq[1],2)-pow(bq[2],2)) {
    return m1+chi;
  } else if (pow(m2+chi,2)>=pow(ap[0],2)-pow(ap[1],2)-pow(ap[2],2)) {
    return m2+chi;
  }

  // Else use balanced solution
  // Note that call to mctcorr also sets m_mctehat, m_mctecm and m_pb

  double mctminsqr = pow(mctcorr(v1,v2,vds,ptm,ecm,mxlo),2);
  double mdmin = mctminmt2(mctminsqr,m1sqr,m2sqr,chisqr);
  if (m_pb == 0) return mdmin;

  double mctecmsqr = pow(m_mctecm,2);
  double mctehatsqr = pow(m_mctehat,2);
  double mctmaxsqr = fmax(mctecmsqr,mctehatsqr);

  double mctdsqr[4] = {(chi*(3.0*m1sqr+m2sqr)+2.0*m1*(m1sqr+m2sqr))/(chi+m1),
		       (chi*(3.0*m1sqr+m2sqr)-2.0*m1*(m1sqr+m2sqr))/(chi-m1),
		       (chi*(3.0*m2sqr+m1sqr)+2.0*m2*(m2sqr+m1sqr))/(chi+m2),
		       (chi*(3.0*m2sqr+m1sqr)-2.0*m2*(m2sqr+m1sqr))/(chi-m2)};
  mdmin = fmin(mdmin,mctminmt2(mctmaxsqr,m1sqr,m2sqr,chisqr));

  for (int i=0;i<4;i++)
    if ( (mctdsqr[i]>mctminsqr) && (mctdsqr[i]<mctmaxsqr) )
      mdmin = fmin(mdmin,mctminmt2(mctdsqr[i],m1sqr,m2sqr,chisqr));
  
  return fmax(fmax(mdmin,m1+chi),m2+chi);
}

double mctminmt2(const double mctsqr,const double m1sqr,
			  const double m2sqr, const double chisqr)
{
  double at = 0.5*(mctsqr - m1sqr - m2sqr);
  return sqrt(fmax(chisqr+at+sqrt(fmax((1.0+(4.0*chisqr)/(2.0*at-m1sqr-m2sqr))
				       *(pow(at,2)-m1sqr*m2sqr),0.0)),0.0));
}

double mt2neg(const double v1[4],const double v2[4]
		       ,const double ptm[2],const double mxlo)
{
  // Non boost-corrected analytical mt2
  double m1sqr = pow(v1[0],2)-pow(v1[1],2)-pow(v1[2],2)-pow(v1[3],2);
  double m2sqr = pow(v2[0],2)-pow(v2[1],2)-pow(v2[2],2)-pow(v2[3],2);
  double chisqr = pow(mxlo,2);
  double m1 = sqrt(fmax(m1sqr,0.0));
  double m2 = sqrt(fmax(m2sqr,0.0));
  double chi = mxlo;

  double qtest[3] = {0.0,ptm[0]-(chi/m1)*v1[1],ptm[1]-(chi/m1)*v1[2]};
  qtest[0] = sqrt(chisqr+pow(qtest[1],2)+pow(qtest[2],2));
  double ptest[3] = {0.0,ptm[0]-(chi/m2)*v2[1],ptm[1]-(chi/m2)*v2[2]};
  ptest[0] = sqrt(chisqr+pow(ptest[1],2)+pow(ptest[2],2));

  double et1 = sqrt(fmax(v1[0]*v1[0]-v1[3]*v1[3],0.0));
  double et2 = sqrt(fmax(v2[0]*v2[0]-v2[3]*v2[3],0.0));

  double bq[3] = {et2+qtest[0],v2[1]+qtest[1],v2[2]+qtest[2]};
  double ap[3] = {et1+ptest[0],v1[1]+ptest[1],v1[2]+ptest[2]};

  // Use unbalanced solutions
  if (pow(m1+chi,2)>=pow(bq[0],2)-pow(bq[1],2)-pow(bq[2],2)) {
    return m1+chi;
  } else if (pow(m2+chi,2)>=pow(ap[0],2)-pow(ap[1],2)-pow(ap[2],2)) {
    return m2+chi;
  }

  // Else use balanced solution
  double mctminsqr = pow(mct(v1,v2),2);
  double mdmin = mctminmt2(mctminsqr,m1sqr,m2sqr,chisqr);
  
  return fmax(fmax(mdmin,m1+chi),m2+chi);
}


double mcy(const double v1[4],const double v2[4]
		       ,const double vds[4],const double ptm[2])
{
// v1, v2 are the 4-vectors of the two aggregated decay products
// vds is the 4-vector of the downstream objects (excluding v1 and v2)
// ptm is the pTmiss 2-vector {pxmiss,pymiss}

  double ptus[2] = {-v1[1]-v2[1]-vds[1]-ptm[0],-v1[2]-v2[2]-vds[2]-ptm[1]};
  double pb = sqrt(pow(ptus[0],2)+pow(ptus[1],2));

  if (pb==0) {

    return mct(v1,v2);    

  } else {

// Transform to new basis in transverse plane
// ISR goes in +ve x direction

    double vb1[4] = {v1[0],(v1[1]*ptus[0]+v1[2]*ptus[1])/pb
	      ,(v1[1]*ptus[1]-v1[2]*ptus[0])/pb,v1[3]};
    double vb2[4] = {v2[0],(v2[1]*ptus[0]+v2[2]*ptus[1])/pb
	      ,(v2[1]*ptus[1]-v2[2]*ptus[0])/pb,v2[3]};
    double vey1 = sqrt(fmax(pow(vb1[0],2)-pow(vb1[1],2)-pow(vb1[3],2),0.0));
    double vey2 = sqrt(fmax(pow(vb2[0],2)-pow(vb2[1],2)-pow(vb2[3],2),0.0));
    return sqrt(fmax(pow(vey1+vey2,2)-pow(vb1[2]-vb2[2],2),0.0));
  }
}

double mcx(const double v1[4],const double v2[4]
		       ,const double vds[4],const double ptm[2])
{
// v1, v2 are the 4-vectors of the two aggregated decay products
// vds is the 4-vector of the downstream objects (excluding v1 and v2)
// ptm is the pTmiss 2-vector {pxmiss,pymiss}

  double ptus[2] = {-v1[1]-v2[1]-vds[1]-ptm[0],-v1[2]-v2[2]-vds[2]-ptm[1]};
  double pb = sqrt(pow(ptus[0],2)+pow(ptus[1],2));

  if (pb==0) {

    return mct(v1,v2);    

  } else {

// Transform to new basis in transverse plane
// ISR goes in +ve x direction

    double vb1[4] = {v1[0],(v1[1]*ptus[0]+v1[2]*ptus[1])/pb
	      ,(v1[1]*ptus[1]-v1[2]*ptus[0])/pb,v1[3]};
    double vb2[4] = {v2[0],(v2[1]*ptus[0]+v2[2]*ptus[1])/pb
	      ,(v2[1]*ptus[1]-v2[2]*ptus[0])/pb,v2[3]};
    double vex1 = sqrt(fmax(pow(vb1[0],2)-pow(vb1[2],2)-pow(vb1[3],2),0.0));
    double vex2 = sqrt(fmax(pow(vb2[0],2)-pow(vb2[2],2)-pow(vb2[3],2),0.0));
    return sqrt(fmax(pow(vex1+vex2,2)-pow(vb1[1]-vb2[1],2),0.0));
  }
}


bool metFiltersPasses2(AC1B &tree_, std::vector<TString> metFlags) {

  bool passed = true;
  unsigned int nFlags = metFlags.size();
  //  std::cout << "MEt filters : " << std::endl;
  for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
    TString flagName(it->first);
    //    std::cout << it->first << " : " << it->second << std::endl;
    for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
      if (flagName.Contains(metFlags[iFilter])) {
	if (it->second==0) {
	  passed = false;
	  break;
	}
      }
    }
  }
  //  std::cout << "Passed : " << passed << std::endl;
  return passed;

}
#endif
