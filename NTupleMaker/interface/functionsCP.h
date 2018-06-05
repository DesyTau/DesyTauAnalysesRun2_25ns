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




void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2){
  if(analysisTree->tau_decayMode[tauIndex1]>1||analysisTree->tau_decayMode[tauIndex2]>1){
    std::cout << "WRONG DECAY MODE" << std::endl;
    return;
  }
  //4-momenta of charged and neutral Pi
  TLorentzVector tau1Prong,tau1Pi0orIP;
  TLorentzVector tau2Prong,tau2Pi0orIP;
  tau1Prong=chargedPivec(analysisTree,tauIndex1);
  tau2Prong=chargedPivec(analysisTree,tauIndex2);

  if(analysisTree->tau_decayMode[tauIndex1]==1||analysisTree->tau_decayMode[tauIndex2]==1){
    tau1Pi0orIP=neutralPivec(analysisTree,tauIndex1);
    tau2Pi0orIP=neutralPivec(analysisTree,tauIndex2);
  }
  TLorentzVector Prongsum;
  //finding the boosting vector to the charged Pi ZMF
  Prongsum=tau1Prong+tau2Prong;
  TVector3 boost;
  boost=Prongsum.BoostVector();
  boost*=-1;
  //boosting every 4-momenta in the charged Pi ZMF
  tau1Prong.Boost(boost);
  tau1Pi0orIP.Boost(boost);
  tau2Prong.Boost(boost);
  tau2Pi0orIP.Boost(boost);

  //getting the spatial components of the pi momenta
  TVector3 norm1,norm2;
  TVector3 vec1chPi,vec1Pi0;
  TVector3 vec2chPi,vec2Pi0;
  vec1chPi=tau1Prong.Vect();
  vec1Pi0=tau1Pi0orIP.Vect();
  vec2chPi=tau2Prong.Vect();
  vec2Pi0=tau2Pi0orIP.Vect();

  norm1=vec1chPi.Cross(vec1Pi0);
  norm1*=(1/norm1.Mag());
  norm2=vec2chPi.Cross(vec2Pi0);
  norm2*=(1/norm2.Mag());

  //filling the variables for the CP measurement
  otree->acotautau=acos(norm1*norm2);

  otree->tau1DecayPlaneX=norm1.X();
  otree->tau1DecayPlaneY=norm1.Y();
  otree->tau1DecayPlaneZ=norm1.Z();
  otree->tau2DecayPlaneX=norm2.X();
  otree->tau2DecayPlaneY=norm2.Y();
  otree->tau2DecayPlaneZ=norm2.Z();
};

TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex){
  int sign = -1;
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  if(analysisTree->tau_charge[tauIndex]>0) sign = 1;
  int piIndex=-1;
  float maxE=-1;
  for(int i=0;i<ncomponents;i++){//selects the highest energy Pi with the same sign of the tau
    if((analysisTree->tau_constituents_pdgId[tauIndex][i]*sign)==211&&analysisTree->tau_constituents_e[tauIndex][i]>maxE){
      piIndex=i;
      maxE=analysisTree->tau_constituents_e[tauIndex][piIndex];
    }
  }
  TLorentzVector chargedPi;
  chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
		       analysisTree->tau_constituents_py[tauIndex][piIndex],
		       analysisTree->tau_constituents_pz[tauIndex][piIndex],
		       analysisTree->tau_constituents_e[tauIndex][piIndex]);
  
  return chargedPi;
};

TLorentzVector neutralPivec(const AC1B * analysisTree, int tauIndex){
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex=-1;
  TLorentzVector neutralPi;
  
  
  for(int i=0;i<ncomponents;i++){
    if(analysisTree->tau_constituents_pdgId[tauIndex][i]==22||abs(analysisTree->tau_constituents_pdgId[tauIndex][i])==11){

      TLorentzVector neutralpart;      //momenta for photons, electrons and positrons

      neutralpart.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			     analysisTree->tau_constituents_py[tauIndex][piIndex],
			     analysisTree->tau_constituents_pz[tauIndex][piIndex],
			     analysisTree->tau_constituents_e[tauIndex][piIndex]);
      neutralPi+=neutralpart;
      }
  }

  return neutralPi;
};
