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


// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
bool tightJetID(const AC1B &tree_, int jet, int era){
	bool tightJetID = false;
	float energy = tree_.pfjet_e[jet];
	float eta = tree_.pfjet_eta[jet];
	float NHF = tree_.pfjet_neutralhadronicenergy[jet] / energy;									 
	float NEMF = tree_.pfjet_neutralemenergy[jet] / energy;  											 
	float NumConst = tree_.pfjet_chargedmulti[jet] + tree_.pfjet_neutralmulti[jet];     
	float CHM = tree_.pfjet_chargedmulti[jet];																		 											
	float MUF = tree_.pfjet_muonenergy[jet] / energy;															 
	float CHF = tree_.pfjet_chargedhadronicenergy[jet] / energy;									 
	float CEMF = tree_.pfjet_chargedemenergy[jet] / energy;												 
	float NumNeutralParticle  = tree_.pfjet_neutralmulti[jet]; 																		 													
	
	if (era == 2016){
		if (fabs(eta) <= 2.7)
			tightJetID = (NHF < 0.90 && NEMF < 0.90 && NumConst > 1) && ((abs(eta) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.99) || abs(eta) > 2.4);
		else if (fabs(eta) <= 3.0)
			tightJetID = NHF < 0.98 && NEMF > 0.01 && NumNeutralParticle > 2;
		else
			tightJetID = NEMF < 0.90 && NumNeutralParticle > 10;
	}
	else if (era == 2017){
		if (fabs(eta) <= 2.7)
			tightJetID = (NHF < 0.90 && NEMF < 0.90 && NumConst > 1) && ((abs(eta) <= 2.4 && CHF > 0 && CHM > 0) || abs(eta) > 2.4);
		else if (fabs(eta) <= 3.0)
			tightJetID = NEMF < 0.99 && NEMF > 0.02 && NumNeutralParticle > 2;
		else
			tightJetID = NEMF < 0.90 && NHF > 0.02 && NumNeutralParticle > 10;
	}
	else if (era == 2018){
		if (fabs(eta) <= 2.6)
			tightJetID = CEMF < 0.8 && CHM > 0 && CHF > 0 && NumConst > 1 && NEMF < 0.9 && MUF < 0.8 && NHF < 0.9;
		else if (fabs(eta) <= 2.7)
			tightJetID = CEMF < 0.8 && CHM > 0 && NEMF < 0.99 && MUF < 0.8 && NHF < 0.9;
		else if (fabs(eta) <= 3.0)
			tightJetID = NEMF > 0.02 && NEMF < 0.99 && NumNeutralParticle > 2;
		else
			tightJetID = NEMF < 0.90 && NHF > 0.2 && NumNeutralParticle > 10;
	}
	else
	{
		std::cout << "era is not 2016, 2017, 2018, exiting" << '\n';
		exit(-1);
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

bool jetPUID(AC1B &tree_, int jet, TString wp="Tight"){
  // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#Working_points
  //4 Eta Categories 0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
  //Tight Id
  size_t etaIndex=0;
  vector<float> Pt010 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt1020 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt2030 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt3050 = { 0.86, -0.10, -0.05, -0.01};
  if(wp=="Medium"){
    //Medium Id
    Pt010 = { 0.18, -0.55, -0.42, -0.36};
    Pt1020 = { 0.18, -0.55, -0.42, -0.36};
    Pt2030 = { 0.18, -0.55, -0.42, -0.36};
    Pt3050 = { 0.61, -0.35, -0.23, -0.17};
  }else if(wp=="Loose"){
    //Loose Id
    Pt010 = {-0.97, -0.68, -0.53, -0.47};
    Pt1020 = {-0.97, -0.68, -0.53, -0.47};
    Pt2030 = {-0.97, -0.68, -0.53, -0.47};
    Pt3050 = {-0.89, -0.52, -0.38, -0.30};
  }


  float pt = tree_.pfjet_pt[jet];
  float eta = tree_.pfjet_eta[jet];
  float mva =  tree_.pfjet_pu_jet_fullDisc_mva[jet];

  if(fabs(eta) < 2.5) etaIndex=0;
  else if(fabs(eta) < 2.75) etaIndex=1;
  else if(fabs(eta) < 3.0)  etaIndex=2;
  else if(fabs(eta) < 5.0)  etaIndex=3;

  float cut = -1.;

  if(pt<10) cut = Pt010[etaIndex];
  else if(pt<20) cut = Pt1020[etaIndex];
  else if(pt<30) cut = Pt2030[etaIndex];
  else if(pt<50) cut = Pt3050[etaIndex];


  bool passedPUID = (bool) mva > cut;
  if (pt >= 50) passedPUID = true;
  if (eta >=5.0)passedPUID = false;

  return passedPUID;

}




#endif
