#ifndef Jets_h
#define Jets_h
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
using namespace std;


bool tightJetiD(AC1B &tree_ ,int jet){
	bool tightJetID;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
                float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
                float nem = tree_.pfjet_neutralemenergy[jet]/energy;
                float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
                float chm = tree_.pfjet_chargedmulti[jet] ;
                float muf = tree_.pfjet_muonenergy[jet]/energy;
                float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
                float elf = tree_.pfjet_chargedemenergy[jet]/energy;
		float nnpart = tree_.pfjet_neutralmulti[jet];
		if (fabs(eta)<=3.)
                
			tightJetID = (nhf<0.90 && nem<0.90 && npr>1) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && elf<0.99) || fabs(eta)>2.4) && fabs(eta)<=3.0; ///
		
		else    
			tightJetID = (nem<0.90 && nnpart>10 && fabs(eta)>3.0 ) ;
	return tightJetID;

}

bool looseJetiD(AC1B &tree_, int jet){

	bool looseJetID = false;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
                float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
                float nem = tree_.pfjet_neutralemenergy[jet]/energy;
                float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
                float chm = tree_.pfjet_chargedmulti[jet] ;
                float muf = tree_.pfjet_muonenergy[jet]/energy;
                float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
                float elf = tree_.pfjet_chargedemenergy[jet]/energy;
		float nnpart = tree_.pfjet_neutralmulti[jet];
		if (fabs(eta)<=3.)
                
			looseJetID = (nhf<0.99 && nem<0.99 && npr>1) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && elf<0.99) || fabs(eta)>2.4) && fabs(eta)<=3.0; ///
		
		else  
			looseJetID = (nem<0.90 && nnpart>10 && fabs(eta)>3.0 ) ;

	return looseJetID;

}


bool tightLepVetoJetiD(AC1B &tree_, int jet){
	bool tightLepVetoJetID = false;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
                float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
                float nem = tree_.pfjet_neutralemenergy[jet]/energy;
                float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
                float chm = tree_.pfjet_chargedmulti[jet] ;
                float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
                float muf = tree_.pfjet_muonenergy[jet]/energy;
                float elf = tree_.pfjet_chargedemenergy[jet]/energy;

 	        tightLepVetoJetID = (nhf<0.90 && nem<0.90 && npr>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && elf<0.90) || fabs(eta)>2.4) && fabs(eta)<=3.0 ;
	return tightLepVetoJetID;

}


double JERSF(double eta, double weight){

  double absEta = fabs(eta);

  if ( absEta<0.8) { weight *=1.061;}
  if ( absEta>0.8 && absEta<1.3) { weight *=1.088;}
  if ( absEta>1.3 && absEta<1.9) { weight *=1.106;}
  if ( absEta>1.9 && absEta<2.5) { weight *=1.126;}
  if ( absEta>2.5 && absEta<3.0) { weight *=1.343;}
  if ( absEta>3.0 && absEta<3.2) { weight *=1.303;}
  if ( absEta>3.2 && absEta<3.5) { weight *=1.320;}

  return weight;

}	




#endif
