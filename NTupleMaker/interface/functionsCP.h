#include "TMath.h"
#include "TLorentzVector.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);
TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector neutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector ipVec(const AC1B * analysisTree, int tauIndex);

double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2) {


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

  double sign = vecPi2 * vecRef1transv.Cross(vecRef2transv);
  if (sign<0)
    acop = 2.0*TMath::Pi() - acop;

  return acop;

}

void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2){

  bool correctDecay1 = analysisTree->tau_decayMode[tauIndex1]<10;
  bool correctDecay2 = analysisTree->tau_decayMode[tauIndex2]<10;
  bool correctDecay = correctDecay1 && correctDecay2;

  if(!correctDecay){
    otree->acotautau = -9999;
    return;
  }
  //4-momenta of charged and neutral Pi
  TLorentzVector tau1Prong=chargedPivec(analysisTree,tauIndex1);
  TLorentzVector tau2Prong=chargedPivec(analysisTree,tauIndex2);
  TLorentzVector tau1ref;
  TLorentzVector tau2ref;
  double y1 = 1;
  double y2 = 1;
  if (analysisTree->tau_decayMode[tauIndex1]==0) 
    tau1ref = ipVec(analysisTree,tauIndex1);
  else {
    tau1ref = neutralPivec(analysisTree,tauIndex1);
    y1 = tau1Prong.E() - tau1ref.E();
  }
  if (analysisTree->tau_decayMode[tauIndex2]==0)
    tau2ref = ipVec(analysisTree,tauIndex1);
  else {
    tau2ref = neutralPivec(analysisTree,tauIndex2);
    y2 = tau2Prong.E() - tau2ref.E();
  }

  double y = y1*y2; 

  //filling the variables for the CP measurement
  double acop = acoCP(tau1Prong,tau2Prong,
		      tau1ref,tau2ref);
  if (y<0) {
    acop = acop + TMath::Pi();
  }
  if (acop>2*TMath::Pi()) {
    acop = acop - 2*TMath::Pi();
  }

  otree->acotautau = acop;

};

TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex){

  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex=-1;
  float maxPt=-1;
  for(int i=0;i<ncomponents;i++){//selects the highest energy Pi with the same sign of the tau
    TLorentzVector lvector; lvector.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
					    analysisTree->tau_constituents_py[tauIndex][i],
					    analysisTree->tau_constituents_pz[tauIndex][i],
					    analysisTree->tau_constituents_e[tauIndex][i]);
    double Pt = lvector.Pt();
    if(fabs(double(analysisTree->tau_constituents_charge[tauIndex][i]))>0.5&&Pt>maxPt){
      piIndex=i;
      maxPt = Pt;
    }
  }
  TLorentzVector chargedPi; chargedPi.SetXYZT(0,0,0,0);
  if (piIndex>-1) {
    chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			 analysisTree->tau_constituents_py[tauIndex][piIndex],
			 analysisTree->tau_constituents_pz[tauIndex][piIndex],
			 analysisTree->tau_constituents_e[tauIndex][piIndex]);
  }  

  return chargedPi;
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
  TLorentzVector vec; vec.SetXYZT(0.,0.,0.,0.);

  int piIndex = -1;
  float maxPt = -1;
  for(int i=0;i<ncomponents;i++){//selects the highest energy Pi with the same sign of the tau

    TLorentzVector lvector; lvector.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
                                            analysisTree->tau_constituents_py[tauIndex][i],
                                            analysisTree->tau_constituents_pz[tauIndex][i],
                                            analysisTree->tau_constituents_e[tauIndex][i]);
    double Pt = lvector.Pt();
    if(fabs(double(analysisTree->tau_constituents_charge[tauIndex][i]))>0.5&&Pt>maxPt){
      piIndex=i;
      maxPt = Pt;
    }
  }

  if (piIndex>-1) {
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

  }

  return vec;

} 
