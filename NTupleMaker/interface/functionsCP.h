#include "TMath.h"
#include "TLorentzVector.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);
TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector neutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector ipVec(const AC1B * analysisTree, int tauIndex);
int chargedPiIndex(const AC1B * analysisTree, int tauIndex);
void gen_acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);
TLorentzVector gen_chargedPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector gen_neutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex);
int gen_chargedPiIndex(const AC1B * analysisTree, int tauIndex);

double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2) {

  double y1 = 1;
  double y2 = 1;
  y1 = Pi1.E() - ref1.E();
  y2 = Pi2.E() - ref2.E();

  double y = y1*y2; 



  TLorentzVector Prongsum = Pi1 + Pi2;
  TVector3 boost = -Prongsum.BoostVector();
  Pi1.Boost(boost);
  Pi2.Boost(boost);
  ref1.Boost(boost);
  ref2.Boost(boost);

  // get 3-vectors
  TVector3 vecPi1 = Pi1.Vect();
  TVector3 vecPi2 = Pi2.Vect();
  TVector3 vecRef1 = ref1.Vect();
  TVector3 vecRef2 = ref2.Vect();

  // normalize them 
  vecPi1 *= 1/vecPi1.Mag();
  vecPi2 *= 1/vecPi2.Mag();
  vecRef1 *= 1/vecRef1.Mag();
  vecRef2 *= 1/vecRef2.Mag();

  // transverse components  
  TVector3 vecRef1transv = vecRef1 - vecPi1*(vecPi1*vecRef1);
  TVector3 vecRef2transv = vecRef2 - vecPi2*(vecPi2*vecRef2);

  double acop = TMath::ACos(vecRef1transv*vecRef1transv);

  if (y<0) {
    acop = acop + TMath::Pi();
  }
  if (acop>2*TMath::Pi()) {
    acop = acop - 2*TMath::Pi();
  }

  double sign = vecPi2 * vecRef1transv.Cross(vecRef2transv);
  if (sign<0) acop = 2.0*TMath::Pi() - acop;


  return acop;

}

void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2){

  bool correctDecay1 = analysisTree->tau_decayMode[tauIndex1]<10;
  bool correctDecay2 = analysisTree->tau_decayMode[tauIndex2]<10;
  bool correctDecay = correctDecay1 && correctDecay2;

  if(!correctDecay){
    otree->acotautau_00 = -9999;
    otree->acotautau_10 = -9999;
    otree->acotautau_01 = -9999;
    otree->acotautau_11 = -9999;
    return;
  }
  //4-momenta of charged and neutral Pi
  TLorentzVector tau1Prong;

  tau1Prong=chargedPivec(analysisTree,tauIndex1);
  TLorentzVector tau2Prong;
  tau2Prong=chargedPivec(analysisTree,tauIndex2);

  TLorentzVector tau1IP;
  tau1IP = ipVec(analysisTree,tauIndex1);
  TLorentzVector tau1Pi0;
  tau1Pi0.SetXYZT(0.,0.,0.,0.);
  TLorentzVector tau2IP;
  tau2IP = ipVec(analysisTree,tauIndex2);
  TLorentzVector tau2Pi0;
  tau2Pi0.SetXYZT(0.,0.,0.,0.);
 
  if (analysisTree->tau_decayMode[tauIndex1]==1) tau1Pi0 = neutralPivec(analysisTree,tauIndex1);
  
  if (analysisTree->tau_decayMode[tauIndex2]==1) tau2Pi0 = neutralPivec(analysisTree,tauIndex2);

  double acop = 0.; 
  
  otree->acotautau_00=acoCP(tau1Prong,tau2Prong,tau1IP,tau2IP);

  if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==0)
    otree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP);

  if (analysisTree->tau_decayMode[tauIndex1]==0&&analysisTree->tau_decayMode[tauIndex2]==1)
    otree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0);

  if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==1){
    otree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP);
    otree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0);
    otree->acotautau_11=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2Pi0);
  }
};

TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex){
  int piIndex=-1;
  piIndex=chargedPiIndex(analysisTree,tauIndex);
  TLorentzVector chargedPi;
  chargedPi.SetXYZT(0,0,0,0);
  if (piIndex>-1) {
    chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			 analysisTree->tau_constituents_py[tauIndex][piIndex],
			 analysisTree->tau_constituents_pz[tauIndex][piIndex],
			 analysisTree->tau_constituents_e[tauIndex][piIndex]);
  }  

  return chargedPi;
};

int chargedPiIndex(const AC1B * analysisTree, int tauIndex){

  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex=-1;
  float maxPt=-1;
  int sign = -1;
  if(analysisTree->tau_charge[tauIndex]>0) sign = 1; 
  for(int i=0;i<ncomponents;i++){//selects the highest energy Pi with the same sign of the tau 
    if((analysisTree->tau_constituents_pdgId[tauIndex][i]*sign)==211){
      TLorentzVector lvector; lvector.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
					      analysisTree->tau_constituents_py[tauIndex][i],
					      analysisTree->tau_constituents_pz[tauIndex][i],
					      analysisTree->tau_constituents_e[tauIndex][i]);
      double Pt = lvector.Pt();
      if(Pt>maxPt){
	piIndex=i;
	maxPt = Pt;
      }
    }
  }
  /*
  TLorentzVector chargedPi;
  chargedPi.SetXYZT(0,0,0,0);
  if (piIndex>-1) {
    chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			 analysisTree->tau_constituents_py[tauIndex][piIndex],
			 analysisTree->tau_constituents_pz[tauIndex][piIndex],
			 analysisTree->tau_constituents_e[tauIndex][piIndex]);
  }  
  */
  return piIndex;
};

TLorentzVector neutralPivec(const AC1B * analysisTree, int tauIndex){
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex=-1;
  TLorentzVector neutralPi; neutralPi.SetXYZT(0.,0.,0.,0.);
  
  
  for(int i=0;i<ncomponents;i++){
    if(analysisTree->tau_constituents_pdgId[tauIndex][i]==22||abs(analysisTree->tau_constituents_pdgId[tauIndex][i])==11){

      TLorentzVector neutralpart;      //momenta for photons, electrons and positrons

      neutralpart.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][i],
			     analysisTree->tau_constituents_py[tauIndex][i],
			     analysisTree->tau_constituents_pz[tauIndex][i],
			     analysisTree->tau_constituents_e[tauIndex][i]);
      neutralPi+=neutralpart;
      }
  }

  return neutralPi;
};

TLorentzVector ipVec(const AC1B * analysisTree, int tauIndex) {

  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);
  int piIndex=-1;
  piIndex=chargedPiIndex(analysisTree,tauIndex);



  if (piIndex>-1) {
    /*//Both methods are equivalent and correct, I simply prefer working with vectors compared to arrays. I can switch back to the previous version any time.
    double vert[3]  = {analysisTree->primvertex_x, analysisTree->primvertex_y, analysisTree->primvertex_z} ; 
    double vpart[3] = {analysisTree->tau_constituents_vx[tauIndex][piIndex],
		       analysisTree->tau_constituents_vy[tauIndex][piIndex],
		       analysisTree->tau_constituents_vz[tauIndex][piIndex]} ;
    double ppart[3] = {analysisTree->tau_constituents_px[tauIndex][piIndex],
                       analysisTree->tau_constituents_py[tauIndex][piIndex],
                       analysisTree->tau_constituents_pz[tauIndex][piIndex]} ;

    double magP = 0;
    for (int ic=0; ic<3; ++ic)
      magP += ppart[ic]*ppart[ic];    
    magP = TMath::Sqrt(magP);
    for (int ic=0; ic<3; ++ic)
      ppart[ic] /= magP;

    double time = 0;
    for (int ic=0; ic<3; ++ic)
      time += (vert[ic]-vpart[ic])*ppart[ic];

    double ip[3] = {0,0,0};
    for (int ic=0; ic<3; ++ic)
      ip[ic] = vpart[ic] + ppart[ic]*time - vert[ic]; 

    vec.SetXYZT(ip[0],ip[1],ip[2],0.);
    */

    TVector3 vertex(analysisTree->primvertex_x,
		    analysisTree->primvertex_y,
		    analysisTree->primvertex_z);
    
    TVector3 secvertex(analysisTree->tau_constituents_vx[tauIndex][piIndex],
		       analysisTree->tau_constituents_vy[tauIndex][piIndex],
		       analysisTree->tau_constituents_vz[tauIndex][piIndex]);
    
    TVector3 momenta(analysisTree->tau_constituents_px[tauIndex][piIndex],
		     analysisTree->tau_constituents_py[tauIndex][piIndex],
		     analysisTree->tau_constituents_pz[tauIndex][piIndex]);

    TVector3 r(0.,0.,0.);
    r=secvertex-vertex;
    double projection=r*momenta/momenta.Mag2();
    TVector3 IP;    
    IP=r-momenta*projection;
    vec.SetVect(IP);
    vec.SetT(0.);
  }

  return vec;

};

 
void gen_acott(const AC1B * analysisTree, Synch17GenTree *gentree, int tauIndex1, int tauIndex2){

  bool correctDecay1 = analysisTree->gentau_decayMode[tauIndex1]<3;
  bool correctDecay2 = analysisTree->gentau_decayMode[tauIndex2]<3;
  bool correctDecay = correctDecay1 && correctDecay2;

  if(!correctDecay){
    gentree->acotautau_00 = -9999;
    gentree->acotautau_10 = -9999;
    gentree->acotautau_01 = -9999;
    gentree->acotautau_11 = -9999;
    return;
  }
  //4-momenta of charged and neutral Pi
  TLorentzVector tau1Prong;

  tau1Prong=gen_chargedPivec(analysisTree,tauIndex1);
  TLorentzVector tau2Prong;
  tau2Prong=gen_chargedPivec(analysisTree,tauIndex2);

  TLorentzVector tau1IP;
  tau1IP = gen_ipVec(analysisTree,tauIndex1);
  TLorentzVector tau1Pi0;
  tau1Pi0.SetXYZT(0.,0.,0.,0.);
  TLorentzVector tau2IP;
  tau2IP = gen_ipVec(analysisTree,tauIndex2);
  TLorentzVector tau2Pi0;
  tau2Pi0.SetXYZT(0.,0.,0.,0.);
 
  if (analysisTree->gentau_decayMode[tauIndex1]==1||analysisTree->gentau_decayMode[tauIndex1]==2) tau1Pi0 = gen_neutralPivec(analysisTree,tauIndex1);
  
  if (analysisTree->gentau_decayMode[tauIndex2]==1||analysisTree->gentau_decayMode[tauIndex2]==2) tau2Pi0 = gen_neutralPivec(analysisTree,tauIndex2);

  double acop = 0.; 
  
  gentree->acotautau_00=acoCP(tau1Prong,tau2Prong,tau1IP,tau2IP);

  if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==0)
    gentree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP);

  if (analysisTree->tau_decayMode[tauIndex1]==0&&analysisTree->tau_decayMode[tauIndex2]==1)
    gentree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0);

  if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==1)
    gentree->acotautau_11=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2Pi0);

};


TLorentzVector gen_chargedPivec(const AC1B * analysisTree, int tauIndex){
  int piIndex=-1;
  piIndex=gen_chargedPiIndex(analysisTree,tauIndex);
  TLorentzVector chargedPi;
  chargedPi.SetXYZT(0,0,0,0);
  if (piIndex>-1) {
    chargedPi.SetPxPyPzE(analysisTree->genparticles_px[piIndex],
			 analysisTree->genparticles_py[piIndex],
			 analysisTree->genparticles_pz[piIndex],
			 analysisTree->genparticles_e[piIndex]);
  }  

  return chargedPi;
};

int gen_chargedPiIndex(const AC1B * analysisTree, int tauIndex){

  int npart = analysisTree->genparticles_count;
  int piIndex=-1;
  float maxPt=-1;
  TLorentzVector Tau;
  double dR;
  const double dRcut=0.3;
  Tau.SetPxPyPzE(analysisTree->gentau_visible_px[tauIndex],
		 analysisTree->gentau_visible_py[tauIndex],
		 analysisTree->gentau_visible_pz[tauIndex],
		 analysisTree->gentau_visible_e[tauIndex]);
  
  for(int i=0;i<npart;i++){//selects the highest energy Pi with the same sign of the tau 
    if(analysisTree->genparticles_pdgid[i]==211&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      double Pt = lvector.Pt();
      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(Pt>maxPt&&dR<dRcut){
	piIndex=i;
	maxPt = Pt;
      }
    }
  }
  /*
  TLorentzVector chargedPi;
  chargedPi.SetXYZT(0,0,0,0);
  if (piIndex>-1) {
    chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			 analysisTree->tau_constituents_py[tauIndex][piIndex],
			 analysisTree->tau_constituents_pz[tauIndex][piIndex],
			 analysisTree->tau_constituents_e[tauIndex][piIndex]);
  }  
  */
  return piIndex;
};

TLorentzVector gen_neutralPivec(const AC1B * analysisTree, int tauIndex){

  int npart = analysisTree->genparticles_count;
  int piIndex=-1;
  TLorentzVector Tau;
  TLorentzVector neutralPi; neutralPi.SetXYZT(0.,0.,0.,0.);
  double dR;
  const double dRcut=0.3;
  Tau.SetPxPyPzE(analysisTree->gentau_visible_px[tauIndex],
		 analysisTree->gentau_visible_py[tauIndex],
		 analysisTree->gentau_visible_pz[tauIndex],
		 analysisTree->gentau_visible_e[tauIndex]);
  
  for(int i=0;i<npart;i++){
    if(analysisTree->genparticles_pdgid[i]==111&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      double Pt = lvector.Pt();
      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut)neutralPi+=lvector;   
    }
  }
  return neutralPi;
};

TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex) {

  int ncomponents = analysisTree->genparticles_count;
  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);
  int piIndex=-1;
  piIndex=gen_chargedPiIndex(analysisTree,tauIndex);



  if (piIndex>-1) {
    /*//Both methods are equivalent and correct, I simply prefer working with vectors compared to arrays. I can switch back to the previous version any time.
    double vert[3]  = {analysisTree->primvertex_x, analysisTree->primvertex_y, analysisTree->primvertex_z} ; 
    double vpart[3] = {analysisTree->tau_constituents_vx[tauIndex][piIndex],
		       analysisTree->tau_constituents_vy[tauIndex][piIndex],
		       analysisTree->tau_constituents_vz[tauIndex][piIndex]} ;
    double ppart[3] = {analysisTree->tau_constituents_px[tauIndex][piIndex],
                       analysisTree->tau_constituents_py[tauIndex][piIndex],
                       analysisTree->tau_constituents_pz[tauIndex][piIndex]} ;

    double magP = 0;
    for (int ic=0; ic<3; ++ic)
      magP += ppart[ic]*ppart[ic];    
    magP = TMath::Sqrt(magP);
    for (int ic=0; ic<3; ++ic)
      ppart[ic] /= magP;

    double time = 0;
    for (int ic=0; ic<3; ++ic)
      time += (vert[ic]-vpart[ic])*ppart[ic];

    double ip[3] = {0,0,0};
    for (int ic=0; ic<3; ++ic)
      ip[ic] = vpart[ic] + ppart[ic]*time - vert[ic]; 

    vec.SetXYZT(ip[0],ip[1],ip[2],0.);
    */

    TVector3 vertex(analysisTree->primvertex_x,
		    analysisTree->primvertex_y,
		    analysisTree->primvertex_z);
    
    TVector3 secvertex(analysisTree->genparticles_vx[piIndex],
		       analysisTree->genparticles_vy[piIndex],
		       analysisTree->genparticles_vz[piIndex]);
    
    TVector3 momenta(analysisTree->genparticles_px[piIndex],
		     analysisTree->genparticles_py[piIndex],
		     analysisTree->genparticles_pz[piIndex]);

    TVector3 r(0.,0.,0.);
    r=secvertex-vertex;
    double projection=r*momenta/momenta.Mag2();
    TVector3 IP;    
    IP=r-momenta*projection;
    vec.SetVect(IP);
    vec.SetT(0.);
  }

  return vec;

};
