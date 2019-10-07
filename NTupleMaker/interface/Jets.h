// source : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016

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

bool tightJetiD_2016(AC1B &tree_ ,int jet){
        // source : https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	bool tightJetID = false;
	float energy = tree_.pfjet_e[jet];
	float eta = tree_.pfjet_eta[jet];
	float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
	float nem = tree_.pfjet_neutralemenergy[jet]/energy;
	float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
	float chm = tree_.pfjet_chargedmulti[jet] ;
	float muf = tree_.pfjet_muonenergy[jet]/energy;
	float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
	float elf = tree_.pfjet_chargedemenergy[jet]/energy;
	float nm  = tree_.pfjet_neutralmulti[jet];
	float nnpart = tree_.pfjet_neutralmulti[jet];
	if (fabs(eta)<=2.7)
	  {
            tightJetID = (nhf<0.90 && nem<0.90 && npr>1) && ((abs(eta)<=2.4 && chf>0 && chm>0 && elf<0.99) || abs(eta)>2.4);
	  }
	else if (fabs(eta)<=3.0)
	  {
            tightJetID = nhf<0.98 && nem>0.01 && nm>2;
	  }
	else
	  {
            tightJetID = nem<0.90 && nm>10;
	  }
	return tightJetID;

}

bool looseJetiD_2016(AC1B &tree_ ,int jet){
        // source : https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	bool tightJetID = false;
	float energy = tree_.pfjet_e[jet];
	float eta = tree_.pfjet_eta[jet];
	float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
	float nem = tree_.pfjet_neutralemenergy[jet]/energy;
	float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
	float chm = tree_.pfjet_chargedmulti[jet] ;
	float muf = tree_.pfjet_muonenergy[jet]/energy;
	float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
	float elf = tree_.pfjet_chargedemenergy[jet]/energy;
	float nm  = tree_.pfjet_neutralmulti[jet];
	float nnpart = tree_.pfjet_neutralmulti[jet];
	if (fabs(eta)<=2.7)
	  {
            tightJetID = (nhf<0.99 && nem<0.99 && npr>1) && ((abs(eta)<=2.4 && chf>0 && chm>0 && elf<0.99) || abs(eta)>2.4);
	  }
	else if (fabs(eta)<=3.0)
	  {
            tightJetID = nhf<0.98 && nem>0.01 && nm>2;
	  }
	else
	  {
            tightJetID = nem<0.90 && nm>10;
	  }
	return tightJetID;

}

bool tightJetiD_2017(AC1B &tree_ ,int jet){
	bool tightJetID_2017 = false;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
                float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
                float nem = tree_.pfjet_neutralemenergy[jet]/energy;
                float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
                float chm = tree_.pfjet_chargedmulti[jet] ;
                float muf = tree_.pfjet_muonenergy[jet]/energy;
                float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
                float elf = tree_.pfjet_chargedemenergy[jet]/energy;
                float nm  = tree_.pfjet_neutralmulti[jet];
		float nnpart = tree_.pfjet_neutralmulti[jet];
      if (fabs(eta)<=2.7)
         {
            tightJetID_2017 = (nhf < 0.9 && nem < 0.9 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0));
         }
      else if (fabs(eta)<=3.0)
         {
            tightJetID_2017 = (nem < 0.99 && nem > 0.02 && nm > 2);
         }
      else
         {
            tightJetID_2017 = nem < 0.9 && nhf > 0.02 && nm > 10;
         }   
      return tightJetID_2017;

}

//Merijn: overload with const 2019 8 2
bool tightJetiD_2017(const AC1B &tree_ ,int jet){
	bool tightJetID_2017 = false;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
                float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
                float nem = tree_.pfjet_neutralemenergy[jet]/energy;
                float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
                float chm = tree_.pfjet_chargedmulti[jet] ;
                float muf = tree_.pfjet_muonenergy[jet]/energy;
                float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
                float elf = tree_.pfjet_chargedemenergy[jet]/energy;
                float nm  = tree_.pfjet_neutralmulti[jet];
		float nnpart = tree_.pfjet_neutralmulti[jet];
      if (fabs(eta)<=2.7)
         {
            tightJetID_2017 = (nhf < 0.9 && nem < 0.9 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0));
         }
      else if (fabs(eta)<=3.0)
         {
            tightJetID_2017 = (nem < 0.99 && nem > 0.02 && nm > 2);
         }
      else
         {
            tightJetID_2017 = nem < 0.9 && nhf > 0.02 && nm > 10;
         }   
      return tightJetID_2017;

}



bool tightJetiD_2018(AC1B &tree_ ,int jet){
	bool tightJetID_2018 = false;
		float energy = tree_.pfjet_e[jet];
		float eta = tree_.pfjet_eta[jet];
      float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
      float nem = tree_.pfjet_neutralemenergy[jet]/energy;
      float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
      float chm = tree_.pfjet_chargedmulti[jet] ;
      float muf = tree_.pfjet_muonenergy[jet]/energy;
      float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
      float elf = tree_.pfjet_chargedemenergy[jet]/energy;
      float nm  = tree_.pfjet_neutralmulti[jet];
		float nnpart = tree_.pfjet_neutralmulti[jet];

      if (fabs(eta)<=2.6)
         {
            tightJetID_2018 = (nhf < 0.9 && nem < 0.9 && npr > 1 && chf>0 && chm > 0);
         }
      else if (fabs(eta)<=2.7)
         {
            tightJetID_2018 = (nhf < 0.9 && nem < 0.99 && chm>0);
         }
      else if (fabs(eta)<=3.0)
         {
            tightJetID_2018 = (nem < 0.99 && nem > 0.02 && nm > 2);
         }
      else
         {
            tightJetID_2018 = nem < 0.9 && nhf > 0.02 && nm > 10;
         }
      return tightJetID_2018;

}


bool tightJetID(const AC1B &tree_, int jet, int era){
	bool tightJetID = false;
	float energy = tree_.pfjet_e[jet];
	float eta = tree_.pfjet_eta[jet];
	float nhf = tree_.pfjet_neutralhadronicenergy[jet] / energy;
	float nem = tree_.pfjet_neutralemenergy[jet] / energy;
	float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];
	float chm = tree_.pfjet_chargedmulti[jet] ;
	float muf = tree_.pfjet_muonenergy[jet] / energy;
	float chf = tree_.pfjet_chargedhadronicenergy[jet] / energy;
	float elf = tree_.pfjet_chargedemenergy[jet] / energy;
	float nm  = tree_.pfjet_neutralmulti[jet];
	float nnpart = tree_.pfjet_neutralmulti[jet];
	
	if (era == 2016){
		if (fabs(eta) <= 2.7)
			tightJetID = (nhf < 0.90 && nem < 0.90 && npr > 1) && ((abs(eta) <= 2.4 && chf > 0 && chm > 0 && elf < 0.99) || abs(eta) > 2.4);
		else if (fabs(eta) <= 3.0)
			tightJetID = nhf < 0.98 && nem > 0.01 && nm > 2;
		else
			tightJetID = nem < 0.90 && nm > 10;
	}
	else if (era == 2017){
		if (fabs(eta) <= 2.7)
			tightJetID = (nhf < 0.90 && nem < 0.90 && npr > 1) && ((abs(eta) <= 2.4 && chf > 0 && chm > 0) || abs(eta) > 2.4);
		else if (fabs(eta) <= 3.0)
			tightJetID = nem < 0.99 && nem > 0.02 && nm > 2;
		else
			tightJetID = nem < 0.90 && nhf > 0.02 && nm > 10;
	}
	return tightJetID;

}


bool looseJetiD(AC1B &tree_, int jet){  // updated recipe for 74x,76x,80x

        bool looseJetID = false;
	float energy = tree_.pfjet_e[jet];
	energy *= tree_.pfjet_energycorr[jet]; // uncorrected energy must be used
	float eta = tree_.pfjet_eta[jet];
	float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
	float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
	float nem = tree_.pfjet_neutralemenergy[jet]/energy;
	float elf = tree_.pfjet_chargedemenergy[jet]/energy;
	float muf = tree_.pfjet_muonenergy[jet]/energy;
	float chm = tree_.pfjet_chargedmulti[jet] ;
	float nm  = tree_.pfjet_neutralmulti[jet];
	float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];

	if (fabs(eta)<=2.7)
	  {
	    looseJetID = (nhf < 0.99 && nem < 0.99 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0 && elf < 0.99));
	  }
	else if (fabs(eta)<=3.0)
	  {
	    looseJetID = (nem < 0.9 && nm > 2);
	  }
	else
	  {
	    looseJetID = nem < 0.9 && nm > 10;
	  }

	return looseJetID;
}

//Merijn 2019 8 2: overload the function with const tree. Note that this function needs to be maintained and kept in synch!
bool looseJetiD(const AC1B &tree_, int jet){  // updated recipe for 74x,76x,80x

        bool looseJetID = false;
	float energy = tree_.pfjet_e[jet];
	energy *= tree_.pfjet_energycorr[jet]; // uncorrected energy must be used
	float eta = tree_.pfjet_eta[jet];
	float chf = tree_.pfjet_chargedhadronicenergy[jet]/energy;
	float nhf = tree_.pfjet_neutralhadronicenergy[jet]/energy;
	float nem = tree_.pfjet_neutralemenergy[jet]/energy;
	float elf = tree_.pfjet_chargedemenergy[jet]/energy;
	float muf = tree_.pfjet_muonenergy[jet]/energy;
	float chm = tree_.pfjet_chargedmulti[jet] ;
	float nm  = tree_.pfjet_neutralmulti[jet];
	float npr = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];

	if (fabs(eta)<=2.7)
	  {
	    looseJetID = (nhf < 0.99 && nem < 0.99 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0 && elf < 0.99));
	  }
	else if (fabs(eta)<=3.0)
	  {
	    looseJetID = (nem < 0.9 && nm > 2);
	  }
	else
	  {
	    looseJetID = nem < 0.9 && nm > 10;
	  }

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
