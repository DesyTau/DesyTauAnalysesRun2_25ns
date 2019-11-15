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
#include "HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h"

#define PI_MASS     0.13957 //All energies and momenta are expressed in GeV
#define MUON_MASS   0.105658

/*Updates Merijn 2018 11 13
-added acott_Impr, which is improved version of acott. Note it takes the decay channel as input also. 
If the channel is e-t or mu-t, it will calculate the lepton vx etc. 
-added ipVec_Lepton. it uses the reco vx, vy ,and vz to calculate the leptonic tau impact vector.
-looking into potential issue with acoCPCOUT
*/

void acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);

void acott_Impr(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2,TString ch);
TLorentzVector chargedPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector neutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector charged_constituents_P4(const AC1B * analysisTree, int tauIndex);
//TLorentzVector ipVec(const AC1B * analysisTree, int tauIndex, Synch17Tree *otree);
TLorentzVector ipVec(const AC1B * analysisTree, int tauIndex); //merijn 2019 4 2: the otree argument is superfluous after all..
int chargedPiIndex(const AC1B * analysisTree, int tauIndex);

void gen_acott(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);
TLorentzVector gen_chargedPivec(const AC1B * analysisTree, int tauIndex, int partId);
TLorentzVector gen_neutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector gen_ThreeProngVec(const AC1B * analysisTree, int tauIndex);
TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex, int piIndex, TVector3 vertex);
TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex, TVector3 vertex);
TVector3 gen_ThreeProngSVertex(const AC1B * analysisTree, int tauIndex);
int gen_chargedPiIndex(const AC1B * analysisTree, int tauIndex, int partId);
vector<int> gen_ThreeProngIndices(const AC1B * analysisTree, int tauIndex);

TLorentzVector ipVec_Lepton(const AC1B * analysisTree, int tauIndex, TString ch);

double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2,
	     bool firstNegative, bool pi01, bool pi02, Synch17Tree *otree);

double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2,
	     bool firstNegative, bool pi01, bool pi02, Synch17GenTree* gentree);

TVector3 calculate_IP_helix_mu(const AC1B * analysisTree, int muIndex, TVector3 vertex_coord);
TVector3 calculate_IP_helix_tauh(const AC1B * analysisTree, int tauIndex, TVector3 vertex_coord);
TVector3 get_refitted_PV_with_BS(const AC1B * analysisTree, int leptonIndex, int tauIndex, int &is_refitted_PV_with_BS);

//Merijn: updated function to do CP calculations
//called with: acott_Impr(&analysisTree,otree,tauIndex,leptonIndex, ch);
//Update call instead with: acott_Impr(&analysisTree,otree,leptonIndex, tauIndex, ch). Throughout code we'll switch between mt, et and tt case!
//in later case, the lepton index must be treated as a tau instead of leptonic index
//currently we'll ignore e-mu

void acott_Impr(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2, TString channel){
  // cout<<"Start acott_Impr"<<endl;
 
  //Merijn 2019 1 10: there may be situations where we pass correctdecay, but ultimately the acotau does NOT get calculated. For these situations, currently acotau seems not initialised.
  otree->acotautau_00 = -9999;
  otree->acotautau_10 = -9999;
  otree->acotautau_01 = -9999;
  otree->acotautau_11 = -9999;

  otree->acotautauPsi_00=-9999;
  otree->acotautauPsi_01=-9999;
  otree->acotautauPsi_10=-9999;
  otree->acotautauPsi_11=-9999;

  //currently needed since cp initialiser only called once..
  otree->VxConstitTau1=-9999;
  otree->VyConstitTau1=-9999;
  otree->VzConstitTau1=-9999;
  
  otree->VxConstitTau2=-9999;
  otree->VyConstitTau2=-9999;
  otree->VzConstitTau2=-9999;

  bool decay1haspion=false;
  bool decay2haspion=false;

  //make here a scan if the first decay is hadronic. Only do for tt case, to avoid possible segfaults if tau not present

  if(channel=="tt"){
  for(unsigned int pindex=0;pindex<analysisTree->tau_constituents_count[tauIndex1];pindex++){
    if(abs(analysisTree->tau_constituents_pdgId[tauIndex1][pindex])==211){
      decay1haspion=true;
      break;}}

  //Vinay
  for(unsigned int pindex=0;pindex<analysisTree->tau_constituents_count[tauIndex2];pindex++){
    if(abs(analysisTree->tau_constituents_pdgId[tauIndex2][pindex])==211){
      decay2haspion=true;
      break;}}
  }
  
  //make here a scan if the second decay is hadronic. Only for e-mu would not require this..
  if(channel=="mt"||channel=="et"){
  for(unsigned int pindex=0;pindex<analysisTree->tau_constituents_count[tauIndex2];pindex++){
    if(abs(analysisTree->tau_constituents_pdgId[tauIndex2][pindex])==211){
    //if(abs(analysisTree->tau_constituents_pdgId[tauIndex2][pindex])==13){
      decay2haspion=true;
      break;}}}
  
  
//Merijn: these definitions can be discussed and adjusted to include 3 prong decays, the allow for mode 10 as well
  bool correctDecay1=false;
  bool correctDecay2=false;
  //now fix for the tree cases. Note: currently we EXCLUDE 3 prong, since NOT implemented yet. 

  if(channel=="mt"||channel=="et"){
    correctDecay1=true;
    correctDecay2 = analysisTree->tau_decayMode[tauIndex2]<10&&decay2haspion;}

  if(channel=="tt"){//just require both to be hadronic
    correctDecay1 = analysisTree->tau_decayMode[tauIndex1]<10&&decay1haspion; 
    correctDecay2 = analysisTree->tau_decayMode[tauIndex2]<10&&decay2haspion;}
  
  bool correctDecay = correctDecay1 && correctDecay2;

//properly inintialise to -9999 if decay not correct
  if(!correctDecay){
    return;
  }
  
  //Merijn: here set 4-momentum of tau1. For tt case, currently take the highest pt charged  pion
  TLorentzVector tau1Prong; //the mectric convention used elsehwere for reference. Its in form px, py, pz, E.

  //Merijn: tauindex1 may point to electron or muon. So use the switch and index to obtain correct kinematics, assuming electron and muon masses
  if(channel=="mt"){    
    double muenergy=TMath::Sqrt(TMath::Power(analysisTree->muon_px[tauIndex1],2)+TMath::Power(analysisTree->muon_py[tauIndex1],2)+TMath::Power(analysisTree->muon_pz[tauIndex1],2)+TMath::Power(0.105658,2));//reconstruct energy from 3 momenta + mass
    tau1Prong.SetPxPyPzE(analysisTree->muon_px[tauIndex1],analysisTree->muon_py[tauIndex1],analysisTree->muon_pz[tauIndex1],muenergy);
    //save the dx, dy, dz for tau1:
    otree->VxConstitTau1=analysisTree->muon_vx[tauIndex1];
    otree->VyConstitTau1=analysisTree->muon_vy[tauIndex1];
    otree->VzConstitTau1=analysisTree->muon_vz[tauIndex1];

    otree->chconst_1_pt=analysisTree->muon_pt[tauIndex1];
    otree->chconst_1_eta=analysisTree->muon_eta[tauIndex1];
    otree->chconst_1_phi=analysisTree->muon_phi[tauIndex1];
  }

  
  if(channel=="et"){
    double elecenergy=TMath::Sqrt(TMath::Power(analysisTree->electron_px[tauIndex1],2)+TMath::Power(analysisTree->electron_py[tauIndex1],2)+TMath::Power(analysisTree->electron_pz[tauIndex1],2)+TMath::Power(0.000511,2));
    tau1Prong.SetPxPyPzE(analysisTree->electron_px[tauIndex1],analysisTree->electron_py[tauIndex1],analysisTree->electron_pz[tauIndex1],elecenergy);
    //save the dx, dy, dz for tau1:
    otree->VxConstitTau1=analysisTree->electron_vx[tauIndex1];
    otree->VyConstitTau1=analysisTree->electron_vy[tauIndex1];
    otree->VzConstitTau1=analysisTree->electron_vz[tauIndex1];

    otree->chconst_1_pt=analysisTree->electron_pt[tauIndex1];
    otree->chconst_1_eta=analysisTree->electron_eta[tauIndex1];
    otree->chconst_1_phi=analysisTree->electron_phi[tauIndex1];
  }

  if(channel=="tt"){
    tau1Prong=chargedPivec(analysisTree,tauIndex1);//Merijn: changed to index1. Only works if call for tt!
    int piIndex_=chargedPiIndex(analysisTree,tauIndex1);
    if(piIndex_>-1){
      
      /*
      otree->VxConstitTau1=analysisTree->tau_constituents_vx[tauIndex1][piIndex_];
      otree->VyConstitTau1=analysisTree->tau_constituents_vy[tauIndex1][piIndex_];   
      otree->VzConstitTau1=analysisTree->tau_constituents_vz[tauIndex1][piIndex_];*/

      //Merijn 2019 2 9: updated to new definition
      otree->VxConstitTau1=analysisTree->tau_pca3D_x[tauIndex1];
      otree->VyConstitTau1=analysisTree->tau_pca3D_y[tauIndex1];
      otree->VzConstitTau1=analysisTree->tau_pca3D_z[tauIndex1];

      //Merijn 2019 4 2: comment out and replace with vector from above:    
      otree->chconst_1_pt=tau1Prong.Pt();
      otree->chconst_1_phi=tau1Prong.Phi();
      otree->chconst_1_eta=tau1Prong.Eta();
    }
  }
  
//merijn: we vetoed already if the second tau didn't contain a pion, so can safely calculate the momentum of 2nd prong from hadrons
  TLorentzVector tau2Prong;
  tau2Prong=chargedPivec(analysisTree,tauIndex2);
  int piIndexfortau2=chargedPiIndex(analysisTree,tauIndex2);
  if(piIndexfortau2>-1){
    /*
    otree->VxConstitTau2=analysisTree->tau_constituents_vx[tauIndex2][piIndexfortau2];
    otree->VyConstitTau2=analysisTree->tau_constituents_vy[tauIndex2][piIndexfortau2];   
    otree->VzConstitTau2=analysisTree->tau_constituents_vz[tauIndex2][piIndexfortau2]; */

    //Merijn 2019 2 9: updated to new definition
    otree->VxConstitTau2=analysisTree->tau_pca3D_x[tauIndex2];
    otree->VyConstitTau2=analysisTree->tau_pca3D_y[tauIndex2];
    otree->VzConstitTau2=analysisTree->tau_pca3D_z[tauIndex2];
    
    //Merijn 2019 4 3: replace with the tau2Prong vector. Keep old code in case something broke:
    otree->chconst_2_pt=tau2Prong.Pt();
    otree->chconst_2_phi=tau2Prong.Phi();
    otree->chconst_2_eta=tau2Prong.Eta();
    
  }
    
  
  bool firstNegative = false;//Merijn 2019 1 10: instead I ask if the second is positive. WON'T WORK FOR E-MU CASE!
//comment 2019 2 5: I need to re-evaluate why I implemented this. I think it is due to fact that initially lepton and muon index were innertwined
  if (analysisTree->tau_charge[tauIndex2]>0.0)
    firstNegative = true;

  TLorentzVector tau1IP;
  //Merijn: for leptonic decay calculate in different way than for hadronic
  if(channel=="et"||channel=="mt") tau1IP = ipVec_Lepton(analysisTree,tauIndex1,channel);
  else{tau1IP = ipVec(analysisTree,tauIndex1);}//tt: treat as tau index..
  TLorentzVector tau1Pi0;
  tau1Pi0.SetXYZT(0.,0.,0.,0.);
  
  //Merijn: have to assume here we WON'T look into e-mu case! 
  TLorentzVector tau2IP;
  tau2IP = ipVec(analysisTree,tauIndex2);
  TLorentzVector tau2Pi0;
  tau2Pi0.SetXYZT(0.,0.,0.,0.); 
  
  //Merijn: here calculate neutral pion vectors. Only look for pions in first tau if channel is tt.
  if(channel=="tt"){ if(analysisTree->tau_decayMode[tauIndex1]>=1&&analysisTree->tau_decayMode[tauIndex1]<=3) tau1Pi0 = neutralPivec(analysisTree,tauIndex1);}
  
  if(analysisTree->tau_decayMode[tauIndex2]>=1&&analysisTree->tau_decayMode[tauIndex2]<=3){ 
	  tau2Pi0 = neutralPivec(analysisTree,tauIndex2); }
 
  otree->acotautau_00=acoCP(tau1Prong,tau2Prong,tau1IP,tau2IP,firstNegative,false,false,otree);
  //I think it should work since everything assigned for 3 cases. Note: aco_00 will be filled with mu x 1-prong and mu x 1.1-prong
  
  
  if(channel=="et"||channel=="mt"){//we only fill aco_01 for et and mt
    if (analysisTree->tau_decayMode[tauIndex2]==1){
      otree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true,otree);
    }    
  }
  //NOTE: Commented till DPG
  if(channel=="tt"){//merijn 2019 1 10: only for tt we can use the index to assess a tau, otherwise its a lepton!
    cout<<"accident"<<endl;
    if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==0)
      otree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP,firstNegative,true,false,otree);
    
    if (analysisTree->tau_decayMode[tauIndex1]==0&&analysisTree->tau_decayMode[tauIndex2]==1){
      otree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true,otree);
    }
    
    if (analysisTree->tau_decayMode[tauIndex1]==1&&analysisTree->tau_decayMode[tauIndex2]==1){
      otree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP,firstNegative,true,false,otree);
      otree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true,otree);
      otree->acotautau_11=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2Pi0,firstNegative,true,true,otree);
    }
  }
  
  //  cout<<"End acott_Impr"<<endl;


  //Merijn 2 25: here calculate the alpha-minus observable
  //potentially some calculations could crash when divide by 0..
  //note: for anything else than mu-tau it is not very well defined.. we just pick the negative charged prong its track! this may not make sense in the end..
  //an alternative is to just only fill if the pion of tau2 is negative

  TVector3 IPVEC, PVEC;
  TVector3 ZVEC(0.,0.,1.);

  if(firstNegative){
     IPVEC=tau1IP.Vect();
     PVEC=tau1Prong.Vect();
  }
  else{
    IPVEC=tau2IP.Vect();
    PVEC=tau2Prong.Vect();
  }
  
  //normalise
  IPVEC *= 1/IPVEC.Mag();
  PVEC *= 1/PVEC.Mag();

  TVector3 v1=ZVEC.Cross(PVEC);
  v1*= 1/v1.Mag();

  TVector3 v2=IPVEC.Cross(PVEC);
  v2*= 1/v2.Mag();
  double v3=v1.Dot(v2);
  otree->alphaminus=TMath::ACos(abs(v3));
  
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
	// selects the highest energy Pi with the same sign of the tau
	
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex = -1;
  float maxPt = -1;
  int sign = -1;
  if(analysisTree->tau_charge[tauIndex] > 0) sign = 1; 
  for(int i = 0; i < ncomponents; i++){ 
		if((analysisTree->tau_constituents_pdgId[tauIndex][i]*sign) == 211){
      TLorentzVector lvector; 
			lvector.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
								      analysisTree->tau_constituents_py[tauIndex][i],
								      analysisTree->tau_constituents_pz[tauIndex][i],
								      analysisTree->tau_constituents_e[tauIndex][i]);
      double Pt = lvector.Pt();
      if(Pt > maxPt){
				piIndex = i;
				maxPt = Pt;
      }
    }
  }
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

       
    //Merijn: temporarily add gen vertex instead.. please leave this code for future reference
/*
    TVector3 vertex;
    for (unsigned int igen=0; igen<analysisTree->genparticles_count; ++igen) {
      if (analysisTree->genparticles_pdgid[igen]==23||analysisTree->genparticles_pdgid[igen]==24||
	  analysisTree->genparticles_pdgid[igen]==25||analysisTree->genparticles_pdgid[igen]==35||analysisTree->genparticles_pdgid[igen]==36) {
	vertex.SetX(analysisTree->genparticles_vx[igen]);
	vertex.SetY(analysisTree->genparticles_vy[igen]);
	vertex.SetZ(analysisTree->genparticles_vz[igen]);
	break;
      }
    }
*/
    
    TVector3 secvertex(analysisTree->tau_pca3D_x[tauIndex],
		       analysisTree->tau_pca3D_y[tauIndex],
		       analysisTree->tau_pca3D_z[tauIndex]);
    
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
  else{cout<<"GENUINELY BIZAR: THERE WAS NO CHARGED PION FOUND.. "<<endl;}

  return vec;

};


//Merijn: added function to calculate the ipvec for leptons
TLorentzVector ipVec_Lepton(const AC1B * analysisTree, int tauIndex, TString ch) {
//now the tauindex is for mu or electron!

  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);

  

TVector3 vertex(analysisTree->primvertex_x,
		    analysisTree->primvertex_y,
		    analysisTree->primvertex_z);
  

  /*
  //Merijn: temporarily replace vertex with gen level info
  TVector3 vertex;
  for (unsigned int igen=0; igen<analysisTree->genparticles_count; ++igen) {
    if (analysisTree->genparticles_pdgid[igen]==23||analysisTree->genparticles_pdgid[igen]==24||
	analysisTree->genparticles_pdgid[igen]==25||analysisTree->genparticles_pdgid[igen]==35||analysisTree->genparticles_pdgid[igen]==36) {
      vertex.SetX(analysisTree->genparticles_vx[igen]);
      vertex.SetY(analysisTree->genparticles_vy[igen]);
      vertex.SetZ(analysisTree->genparticles_vz[igen]);
      break;
    }
  }
*/  
    
TVector3 secvertex(0.,0.,0.);
TVector3 momenta(0.,0.,0.);    

if(ch=="et"){
  secvertex.SetXYZ(analysisTree->electron_vx[tauIndex],
		   analysisTree->electron_vy[tauIndex],
		   analysisTree->electron_vz[tauIndex]);
  
  
  momenta.SetXYZ(analysisTree->electron_px[tauIndex],
		 analysisTree->electron_py[tauIndex],
		 analysisTree->electron_pz[tauIndex]);}
 
 if(ch=="mt"){
   secvertex.SetXYZ(analysisTree->muon_vx[tauIndex], //Merijn try somethin 2019 2 9
		    analysisTree->muon_vy[tauIndex],
		    analysisTree->muon_vz[tauIndex]);
   
   
   momenta.SetXYZ(analysisTree->muon_px[tauIndex],
		  analysisTree->muon_py[tauIndex],
		  analysisTree->muon_pz[tauIndex]);
   
 }


    TVector3 r(0.,0.,0.);
    r=secvertex-vertex;
    
    double projection=r*momenta/momenta.Mag2();
    TVector3 IP;    
    IP=r-momenta*projection;
    vec.SetVect(IP);   
    vec.SetT(0.);

  return vec;
};

void gen_acott(const AC1B * analysisTree, Synch17GenTree *gentree, int tauIndex1, int tauIndex2){

  bool correctDecay1 = analysisTree->gentau_decayMode[tauIndex1]<=6||analysisTree->gentau_decayMode[tauIndex1]==8||analysisTree->gentau_decayMode[tauIndex1]==9;
  bool correctDecay2 = analysisTree->gentau_decayMode[tauIndex2]<=6||analysisTree->gentau_decayMode[tauIndex2]==8||analysisTree->gentau_decayMode[tauIndex2]==9;
  bool correctDecay = correctDecay1 && correctDecay2;
  
  bool oneProngPi01 = analysisTree->gentau_decayMode[tauIndex1]<=2 && analysisTree->gentau_decayMode[tauIndex1]>=1;
  bool oneProngPi02 = analysisTree->gentau_decayMode[tauIndex2]<=2 && analysisTree->gentau_decayMode[tauIndex2]>=1;

  bool threeProng1 = analysisTree->gentau_decayMode[tauIndex1]==4;
  bool threeProng2 = analysisTree->gentau_decayMode[tauIndex2]==4;
  
  gentree->acotautau_00 = -9999;
  gentree->acotautau_10 = -9999;
  gentree->acotautau_01 = -9999;
  gentree->acotautau_11 = -9999;

  gentree->acotautau_20 = -9999;
  gentree->acotautau_02 = -9999;
  gentree->acotautau_21 = -9999;
  gentree->acotautau_12 = -9999;
  gentree->acotautau_22 = -9999;

  gentree->VxConstitTau1=-9999;
  gentree->VyConstitTau1=-9999;
  gentree->VzConstitTau1=-9999;

  gentree->VxConstitTau2=-9999;
  gentree->VyConstitTau2=-9999;
  gentree->VzConstitTau2=-9999;

  gentree->decaymode_1 = analysisTree->gentau_decayMode[tauIndex1];
  gentree->decaymode_2 = analysisTree->gentau_decayMode[tauIndex2];

  if (!correctDecay)
    return;
  


  TVector3 vertex;
  for (unsigned int igen=0; igen<analysisTree->genparticles_count; ++igen) {
    if ((analysisTree->genparticles_pdgid[igen]==23||analysisTree->genparticles_pdgid[igen]==24||
	 analysisTree->genparticles_pdgid[igen]==25||analysisTree->genparticles_pdgid[igen]==35||analysisTree->genparticles_pdgid[igen]==36)
	 &&analysisTree->genparticles_isLastCopy[igen]&&analysisTree->genparticles_isPrompt[igen]) {
      vertex.SetX(analysisTree->genparticles_vx[igen]);
      vertex.SetY(analysisTree->genparticles_vy[igen]);
      vertex.SetZ(analysisTree->genparticles_vz[igen]);
      break;
    }
  }
  int partId1 = 211;
  int partId2 = 211;
  if (analysisTree->gentau_decayMode[tauIndex1]==8) partId1 = -13;
  if (analysisTree->gentau_decayMode[tauIndex1]==9) partId1 = -11;
  if (analysisTree->gentau_decayMode[tauIndex2]==8) partId2 = -13;
  if (analysisTree->gentau_decayMode[tauIndex2]==9) partId2 = -11;


  //4-momenta of charged and neutral Pi
  TLorentzVector tau1Prong=gen_chargedPivec(analysisTree,tauIndex1,partId1);
  int piIndex1 = gen_chargedPiIndex(analysisTree,tauIndex1,partId1);

  //Merijn 2019 2 25: add the kinematic information of the constituents
  /*
  TLorentzVector lvector1;
  lvector1.SetXYZT(analysisTree->genparticles_px[piIndex1],
		  analysisTree->genparticles_py[piIndex1],
		  analysisTree->genparticles_pz[piIndex1],
		  analysisTree->genparticles_e[piIndex1]);

  gentree->chconst_1_pt=lvector1.Pt();
  gentree->chconst_1_phi=lvector1.Phi();
  gentree->chconst_1_eta=lvector1.Eta();*/

  //Merijn 2019 4 2: replace with vector already defined..
  gentree->chconst_1_pt=tau1Prong.Pt();
  gentree->chconst_1_phi=tau1Prong.Phi();
  gentree->chconst_1_eta=tau1Prong.Eta();
  
  
  TLorentzVector tau2Prong=gen_chargedPivec(analysisTree,tauIndex2,partId2);
  int piIndex2 = gen_chargedPiIndex(analysisTree,tauIndex2,partId2);

  //Merijn 2019 2 25: add the kinematic information of the constituents
  /*
  TLorentzVector lvector2;
  lvector2.SetXYZT(analysisTree->genparticles_px[piIndex2],
		  analysisTree->genparticles_py[piIndex2],
		  analysisTree->genparticles_pz[piIndex2],
		  analysisTree->genparticles_e[piIndex2]);
  gentree->chconst_2_pt=lvector2.Pt();
  gentree->chconst_2_eta=lvector2.Eta();
  gentree->chconst_2_phi=lvector2.Phi();
  */
  
  //Merijn 2019 4 2: replace with vector already defined..
  gentree->chconst_2_pt=tau2Prong.Pt();
  gentree->chconst_2_eta=tau2Prong.Eta();
  gentree->chconst_2_phi=tau2Prong.Phi();

  //  std::cout << "Mode1 = " << analysisTree->gentau_decayMode[tauIndex1] << "  Mode2 = " << analysisTree->gentau_decayMode[tauIndex2] << std::endl;
  //  std::cout << "pion1 = " << analysisTree->genparticles_pdgid[piIndex1] << "    pion2 = " << analysisTree->genparticles_pdgid[piIndex2] << std::endl;

  TLorentzVector tau1IP;
  tau1IP = gen_ipVec(analysisTree,tauIndex1,piIndex1,vertex);
  TLorentzVector tau1Pi0;
  tau1Pi0.SetXYZT(0.,0.,0.,0.);
  TLorentzVector tau2IP;
  tau2IP = gen_ipVec(analysisTree,tauIndex2,piIndex2,vertex);
  TLorentzVector tau2Pi0;
  tau2Pi0.SetXYZT(0.,0.,0.,0.);

  if (oneProngPi01) tau1Pi0 = gen_neutralPivec(analysisTree,tauIndex1);  
  if (oneProngPi02) tau2Pi0 = gen_neutralPivec(analysisTree,tauIndex2);

  TLorentzVector tau1_3ProngVec;
  tau1_3ProngVec.SetXYZT(0.,0.,0.,0.);
  if (threeProng1){
    tau1_3ProngVec=gen_ThreeProngVec(analysisTree,tauIndex1);
    tau1IP = gen_ipVec(analysisTree,tauIndex1,vertex);
  }
  TLorentzVector tau2_3ProngVec;
  tau2_3ProngVec.SetXYZT(0.,0.,0.,0.);
  if (threeProng2){
    tau2_3ProngVec=gen_ThreeProngVec(analysisTree,tauIndex2);
    tau2IP = gen_ipVec(analysisTree,tauIndex2,vertex);
  }

  gentree->VxConstitTau1=analysisTree->genparticles_vx[piIndex1];
  gentree->VyConstitTau1=analysisTree->genparticles_vy[piIndex1];
  gentree->VzConstitTau1=analysisTree->genparticles_vz[piIndex1];

  gentree->VxConstitTau2=analysisTree->genparticles_vx[piIndex2];
  gentree->VyConstitTau2=analysisTree->genparticles_vy[piIndex2];
  gentree->VzConstitTau2=analysisTree->genparticles_vz[piIndex2];
 
  bool firstNegative = false;
  if (analysisTree->genparticles_pdgid[piIndex1]==-211) firstNegative = true;
  if (analysisTree->genparticles_pdgid[piIndex1]==11) firstNegative = true;
  if (analysisTree->genparticles_pdgid[piIndex1]==13) firstNegative = true;
  
  gentree->acotautau_00=acoCP(tau1Prong,tau2Prong,tau1IP,tau2IP,firstNegative,false,false, gentree);

  if (oneProngPi01)
    gentree->acotautau_10=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP,firstNegative,true,false, gentree);

  //Merijn: I think we need the other configuration also..
    if (oneProngPi02)
    gentree->acotautau_01=acoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true, gentree);
  

  //Merijn 2019 2 9: presume that this is code supposed to keep by Andrea..
  //Merijn 2019 2 9: the gentree arugment was kept. I remove it now..
  if (oneProngPi01&&oneProngPi02)
    gentree->acotautau_11=acoCP(tau1Prong,tau2Prong,tau1Pi0,tau2Pi0,firstNegative,true,true, gentree);

  if(threeProng1)
    gentree->acotautau_20=acoCP(tau1_3ProngVec,tau2Prong,tau1IP,tau2IP,firstNegative,false,false, gentree);

  if(threeProng2)
    gentree->acotautau_02=acoCP(tau1Prong,tau2_3ProngVec,tau1IP,tau2IP,firstNegative,false,false, gentree);

  if (oneProngPi01&&threeProng2)
    gentree->acotautau_12=acoCP(tau1Prong,tau2_3ProngVec,tau1Pi0,tau2IP,firstNegative,true,false, gentree);

  if (threeProng1&&oneProngPi02)
    gentree->acotautau_21=acoCP(tau1_3ProngVec,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true, gentree);

  if(threeProng1&&threeProng2)

    gentree->acotautau_22=acoCP(tau1_3ProngVec,tau2_3ProngVec,tau1IP,tau2IP,firstNegative,false,false, gentree);

  //Merijn 2 25: here calculate the alpha-minus observable
  //potentially some calculations could crash when divide by 0..
  TVector3 IPVEC, PVEC;
  TVector3 ZVEC(0.,0.,1.);

  if(firstNegative){
     IPVEC=tau1IP.Vect();
     PVEC=tau1Prong.Vect();}
  else{
    IPVEC=tau2IP.Vect();
    PVEC=tau2Prong.Vect();}
  
  //normalise
  IPVEC *= 1/IPVEC.Mag();
  PVEC *= 1/PVEC.Mag();

  TVector3 v1=ZVEC.Cross(PVEC);
  v1*= 1/v1.Mag();

  TVector3 v2=IPVEC.Cross(PVEC);
  v2*= 1/v2.Mag();
  double v3=v1.Dot(v2);
  gentree->alphaminus=TMath::ACos(abs(v3));
 
};
 

TLorentzVector gen_chargedPivec(const AC1B * analysisTree, int tauIndex, int partId){
  int piIndex=-1;
  piIndex=gen_chargedPiIndex(analysisTree,tauIndex,partId);
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

int gen_chargedPiIndex(const AC1B * analysisTree, int tauIndex, int partId){
  //TO FIX: once the tau charge at gen level has been added, select the Prong with the correct charge
  int npart = analysisTree->genparticles_count;
  int piIndex=-1;
  float maxPt=-1;
  TLorentzVector Tau;
  double dR;
  double dRcut=0.3;
  Tau.SetPxPyPzE(analysisTree->gentau_visible_px[tauIndex],
		 analysisTree->gentau_visible_py[tauIndex],
		 analysisTree->gentau_visible_pz[tauIndex],
		 analysisTree->gentau_visible_e[tauIndex]);
  
  for(int i=0;i<npart;i++){//selects the highest energy Pi
    if((analysisTree->genparticles_pdgid[i]*analysisTree->gentau_charge[tauIndex])==partId&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)&&analysisTree->genparticles_isLastCopy[i]){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      double Pt = lvector.Pt();
      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut){
	piIndex=i;
	dRcut = dR;
	maxPt = Pt;
      }
    }
  }
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
  int npi0=0;
  for(int i=0;i<npart;i++){
    if(analysisTree->genparticles_pdgid[i]==111&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)&&analysisTree->genparticles_isLastCopy[i]){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      double Pt = lvector.Pt();
      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut) { 
	neutralPi+=lvector;   
	npi0++;
      }
    }
  }
  return neutralPi;
};

TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex, int piIndex, TVector3 vertex) {

  int ncomponents = analysisTree->genparticles_count;
  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);

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
    
    //  vec.SetVect(r);
    vec.SetT(0.);
  }

  return vec;

};

TLorentzVector gen_ipVec(const AC1B * analysisTree, int tauIndex, TVector3 vertex) {

  int ncomponents = analysisTree->genparticles_count;
  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);

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

    TVector3 secvertex = gen_ThreeProngSVertex(analysisTree,tauIndex);
    TLorentzVector LVec = gen_ThreeProngVec(analysisTree,tauIndex);
    TVector3 momenta= LVec.Vect();

    TVector3 r(0.,0.,0.);
    r=secvertex-vertex;
    double projection=r*momenta/momenta.Mag2();
    TVector3 IP;    
    IP=r-momenta*projection;
    vec.SetVect(IP);
    vec.SetT(0.);

  return vec;

};


vector<int> gen_ThreeProngIndices(const AC1B * analysisTree, int tauIndex){
  int npart = analysisTree->genparticles_count;
  TLorentzVector Tau;
  double dR;
  double dRcut=0.5;
  vector<int> ThreeProngIndices;
  Tau.SetPxPyPzE(analysisTree->gentau_visible_px[tauIndex],
		 analysisTree->gentau_visible_py[tauIndex],
		 analysisTree->gentau_visible_pz[tauIndex],
		 analysisTree->gentau_visible_e[tauIndex]);
  ThreeProngIndices.push_back(gen_chargedPiIndex(analysisTree,tauIndex,211));
  for(int i=0;i<npart;i++){
    if((analysisTree->genparticles_pdgid[i]*analysisTree->gentau_charge[tauIndex]==211)&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)&&analysisTree->genparticles_isLastCopy[i]){
      if(i==ThreeProngIndices.at(0))continue;
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut){
	ThreeProngIndices.push_back(i);
      }
    }
  }
  for(int i=0;i<npart;i++){
    if((analysisTree->genparticles_pdgid[i]*analysisTree->gentau_charge[tauIndex]==-211)&&(analysisTree->genparticles_info[i]==12||analysisTree->genparticles_info[i]==5)&&analysisTree->genparticles_isLastCopy[i]){
      TLorentzVector lvector;
      lvector.SetXYZT(analysisTree->genparticles_px[i],
		      analysisTree->genparticles_py[i],
		      analysisTree->genparticles_pz[i],
		      analysisTree->genparticles_e[i]);

      dR=deltaR(Tau.Eta(),Tau.Phi(),lvector.Eta(),lvector.Phi());
      if(dR<dRcut){
	ThreeProngIndices.push_back(i);
      }
    }
  }

  return ThreeProngIndices;
};

TLorentzVector gen_ThreeProngVec(const AC1B * analysisTree, int tauIndex){
  vector<int> ThreeProngIndices = gen_ThreeProngIndices(analysisTree,tauIndex);
  int npart = ThreeProngIndices.size();
  TLorentzVector ThreeProngVec;
  ThreeProngVec.SetXYZT(0.,0.,0.,0.);
  
  if(npart!=3){
    //Meirjn 2010 2 9: Andrea and Merijn commented out a warning message here..
    
    return ThreeProngVec;
  }
  
  for(int piIndex : ThreeProngIndices){//selects the highest energy Pi
    TLorentzVector lvector;    
    lvector.SetXYZT(analysisTree->genparticles_px[piIndex],
		    analysisTree->genparticles_py[piIndex],
		    analysisTree->genparticles_pz[piIndex],
		    analysisTree->genparticles_e[piIndex]);
    
    ThreeProngVec+=lvector;
  }

  return ThreeProngVec;
};

TVector3 gen_ThreeProngSVertex(const AC1B * analysisTree, int tauIndex){
  vector<int> ThreeProngIndices = gen_ThreeProngIndices(analysisTree,tauIndex);

  //TO FIX: once the tau charge at gen level has been added, select the Prong with the correct charge
  int npart = ThreeProngIndices.size();
  TVector3 ThreeProngSVertex={0.,0.,0.};
  
  if(npart!=3){ 
    //cout << "ERROR: found more than 3 prongs!" << endl;
    return ThreeProngSVertex;
  }  
  for(int piIndex : ThreeProngIndices){//selects the highest energy Pi
    TVector3 secvertex(analysisTree->genparticles_vx[piIndex]/(float)npart,
		       analysisTree->genparticles_vy[piIndex]/(float)npart,
		       analysisTree->genparticles_vz[piIndex]/(float)npart);
    ThreeProngSVertex+=secvertex;
  }

  return ThreeProngSVertex;
};


//WORK IN PROGRESS
float Bfunction(TLorentzVector pi, TLorentzVector a1){
  float B = 0.;
  B=(pow(pi.E()*a1.E()-(a1.Vect()).Dot(pi.Vect()),2)-a1.Mag2()*PI_MASS*PI_MASS)/a1.Mag2();
  return B;
};

Float_t gen_A1Polarization(const AC1B * analysisTree, int tauIndex){
  vector<int> ThreeProngIndices = gen_ThreeProngIndices(analysisTree,tauIndex);
  float pol=-9999.;
  int npart = ThreeProngIndices.size();
  if(analysisTree->gentau_decayMode[tauIndex]!=4)return pol;
  else if(npart!=3){
    cout << "ERROR: found " << npart << " prongs!" << endl;
    return npart*1000.;
  }
  TLorentzVector LVec[3];
  TVector3 Vec[3];
  TLorentzVector LVecSum;
  LVecSum.SetXYZT(0.,0.,0.,0.);
  int i=0;
  for(int piIndex : ThreeProngIndices){//selects the highest energy Pi
    LVec[i].SetPxPyPzE(analysisTree->genparticles_px[piIndex],
			 analysisTree->genparticles_py[piIndex],
			 analysisTree->genparticles_pz[piIndex],
			 analysisTree->genparticles_e[piIndex]);

    //Merijn 2019 2 9: commented originally this out due to compilation issue.
    Vec[i]=LVec[i].Vect();
    LVecSum+=LVec[i];
    i++;
  }
  TVector3 VecSum=LVecSum.Vect();
  //for(int j=0;j<3;j++) cout << "Pi" << j+1 << ": " << LVec[j].Px() << ": " << LVec[j].Py() << ": " << LVec[j].Pz() << ": " << LVec[j].E() << endl;
  //cout << "a1: " << LVecSum.Px() << ": " << LVecSum.Py() << ": " << LVecSum.Pz() << ": " << LVecSum.E() << endl;

  pol=Vec[2].Dot(Vec[0].Cross(Vec[1]))/VecSum.Mag();

  //Lambda=1/2*sqrt(2 B1 B2 + 2 B2 B3 + 2 B3 B1 - B1^2 - B2^2 - B3^2)
  float Lambda=0.;
  bool problem=false;
  for(int j=0;j<3;j++)Lambda+=Bfunction(LVec[j],LVecSum); 
  Lambda=pow(Lambda,2);
  //cout << endl;// << "B1+B2+B3: " << Lambda << endl;
  //cout << "(B1+B2+B3)^2: " << Lambda << endl;
  for(int j=0;j<3;j++)Lambda-=2*pow(Bfunction(LVec[j],LVecSum),2);
  if(Lambda<0.) problem=true;

  //Merijn 2019 2 9: comment out cout, since annoying when run over larger samples..
  /*
  if(problem){
    cout << "PROBLEM" <<endl;
    for(int j=0;j<3;j++)cout << " B" << j+1 << ": " << Bfunction(LVec[j],LVecSum);
    cout << endl << "Lambda: " << Lambda <<endl;
    }*/
  Lambda=0.5*pow(Lambda,0.5);

  /*
  if(problem){
    cout << "pol: " << pol <<endl;
    cout << "1/2 * sqrt(Lambda): " << Lambda <<endl;
    cout << "cos(beta)= "<< pol/Lambda <<endl<<endl;
    return 9999.;
  }
  */
  return pol/Lambda;
};


//Merijn: adjust to take otree as well, conventient for debguggin..
double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2,
	     bool firstNegative, bool pi01, bool pi02, Synch17Tree *otree) {

  double y1 = 1;
  double y2 = 1;
  if (pi01)
    y1 = Pi1.E() - ref1.E();
  if (pi02)
    y2 = Pi2.E() - ref2.E();


  double y = y1*y2; 
  //double y = -1; 

  TLorentzVector Prongsum = Pi1 + Pi2;
  TVector3 boost = -Prongsum.BoostVector();
  Pi1.Boost(boost);
  Pi2.Boost(boost);
  ref1.Boost(boost);
  ref2.Boost(boost);

  //  std::cout << "First negative = " << firstNegative << "  pi01 = " << pi01 << "   pi02 = " << pi02 << std::endl;
  //  std::cout << "Px(1) = " << Pi1.Px() << "  Py(1) = " << Pi1.Py() << "  Pz(1) = " << Pi1.Pz() << std::endl;
  //  std::cout << "Px(2) = " << Pi2.Px() << "  Py(2) = " << Pi2.Py() << "  Pz(2) = " << Pi2.Pz() << std::endl;
  //  std::cout << "Ux(1) = " << ref1.Px() << "  Uy(1) = " << ref1.Py() << "  Uz(1) = " << ref1.Pz() << std::endl;
  //  std::cout << "Ux(2) = " << ref2.Px() << "  Uy(2) = " << ref2.Py() << "  Uz(2) = " << ref2.Pz() << std::endl;
  //  std::cout << std::endl;
  
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

  vecRef1transv *= 1/vecRef1transv.Mag();
  vecRef2transv *= 1/vecRef2transv.Mag();

  double acop = TMath::ACos(vecRef1transv*vecRef2transv);
  double acoporiginal=acop;
  double sign = vecPi2 * vecRef1transv.Cross(vecRef2transv);
  double psioriginal =sign;//

  if (firstNegative)
    sign = vecPi1 * vecRef2transv.Cross(vecRef1transv);
  
  if (sign<0) acop = 2.0*TMath::Pi() - acop;

  if (y<0) {
    acop = acop + TMath::Pi();
    if (acop>2*TMath::Pi()) {
      acop = acop - 2*TMath::Pi();
    }
  }

  //Merijn added a flag for events where sign is bad. in generator level generally ~5%. 
  if(isinf(sign)) sign=-3;
  if(sign!=sign) sign=-3;

  //Merijn: I think these assignments are correct..
  if(!pi01&&!pi02)  otree->acotautauPsi_00=sign;

  if (pi01){ otree->acotautauPsi_10=sign;
  }
  if (pi02){
    otree->acotautauPsi_01=sign;
  }
  if (pi01&&pi02){ otree->acotautauPsi_11=sign;
  }

  return acop;
}


double acoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2,
	     bool firstNegative, bool pi01, bool pi02, Synch17GenTree* otree) {

  double y1 = 1;
  double y2 = 1;
  if (pi01)
    y1 = Pi1.E() - ref1.E();
  if (pi02)
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

  vecRef1transv *= 1/vecRef1transv.Mag();
  vecRef2transv *= 1/vecRef2transv.Mag();

  double acop = TMath::ACos(vecRef1transv*vecRef2transv);

  double sign = vecPi2 * vecRef1transv.Cross(vecRef2transv);

  if (firstNegative)
    sign = vecPi1 * vecRef2transv.Cross(vecRef1transv);
  
  if (sign<0) acop = 2.0*TMath::Pi() - acop;

  if (y<0) {
    acop = acop + TMath::Pi();
    if (acop>2*TMath::Pi()) {
      acop = acop - 2*TMath::Pi();
    }
  }

  //Merijn 2019 2 9: I'd like to save the observables for inspection and thus oncommented this again. I added gentree to the argument everywhere.
  
  if(isinf(sign)) sign=-3;
  if(sign!=sign) sign=-3;
      
  //Merijn: I think this is correct assignment
  //Merijn: I think this is correct..
  if(!pi01&&!pi02)  otree->acotautauPsi_00=sign;
  if (pi01){ otree->acotautauPsi_10=sign; }
  if (pi02){ otree->acotautauPsi_01=sign;}
  if (pi01&&pi02){ otree->acotautauPsi_11=sign;}
  
  return acop;
}

TVector3 calculate_IP_helix_mu(const AC1B * analysisTree, int muIndex, TVector3 vertex_coord){
	// helical IP for tau decaying into muon

	double B = analysisTree->muon_Bfield[muIndex];
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> p4_mu;  
	TLorentzVector p4_mu_auxil; 
	std::vector<float> h_param_mu = {};
	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> ref_mu;
	
	ref_mu.SetX(analysisTree->muon_referencePoint[muIndex][0]);
	ref_mu.SetY(analysisTree->muon_referencePoint[muIndex][1]);
	ref_mu.SetZ(analysisTree->muon_referencePoint[muIndex][2]);
	p4_mu_auxil.SetXYZM(analysisTree->muon_px[muIndex], analysisTree->muon_py[muIndex], analysisTree->muon_pz[muIndex], MUON_MASS);
	p4_mu.SetPxPyPzE(p4_mu_auxil.Px(),p4_mu_auxil.Py(),p4_mu_auxil.Pz(),p4_mu_auxil.E());
	for(auto i:  analysisTree->muon_helixparameters[muIndex]) h_param_mu.push_back(i);	
	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> pv(vertex_coord.X(), vertex_coord.Y(), vertex_coord.Z());

	ImpactParameter IP;
	TVector3 IP_helix_mu = IP.CalculatePCA(B, h_param_mu, ref_mu, pv, p4_mu);        

	return IP_helix_mu;
}

TVector3 calculate_IP_helix_tauh(const AC1B * analysisTree, int tauIndex, TVector3 vertex_coord){
	// helical IP for tau_h
	// NB: for calculation will take the 4-momentum of the leading charged hadron, same sign as tau

	double B = analysisTree->tau_Bfield[tauIndex];
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> p4_tau;  
	std::vector<float> h_param_tau = {};
	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> ref_tau;
	
	ref_tau.SetX(analysisTree->tau_referencePoint[tauIndex][0]);
	ref_tau.SetY(analysisTree->tau_referencePoint[tauIndex][1]);
	ref_tau.SetZ(analysisTree->tau_referencePoint[tauIndex][2]);
	
	int leading_pi_index = chargedPiIndex(analysisTree, tauIndex);
	p4_tau.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][leading_pi_index],
										analysisTree->tau_constituents_py[tauIndex][leading_pi_index],
										analysisTree->tau_constituents_pz[tauIndex][leading_pi_index],
										analysisTree->tau_constituents_e[tauIndex][leading_pi_index]);
	
	for(auto i:  analysisTree->tau_helixparameters[tauIndex]) h_param_tau.push_back(i);	
	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> pv(vertex_coord.X(), vertex_coord.Y(), vertex_coord.Z());
	
	ImpactParameter IP;
	TVector3 IP_helix_tau = IP.CalculatePCA(B, h_param_tau, ref_tau, pv, p4_tau);        
	return IP_helix_tau;
}


TVector3 get_refitted_PV_with_BS(const AC1B * analysisTree, int leptonIndex, int tauIndex, bool &is_refitted_PV_with_BS){
	float vtx_x = analysisTree->primvertexwithbs_x; // by default store non-refitted PV with BS constraint if refitted one is not found
	float vtx_y = analysisTree->primvertexwithbs_y;
	float vtx_z = analysisTree->primvertexwithbs_z;
	is_refitted_PV_with_BS = false;

	for(unsigned int i = 0; i < analysisTree->refitvertexwithbs_count; i++)
	{
    if( (leptonIndex == analysisTree->refitvertexwithbs_muIndex[i][0] || leptonIndex == analysisTree->refitvertexwithbs_muIndex[i][1]) &&
        (tauIndex == analysisTree->refitvertexwithbs_tauIndex[i][0] || tauIndex == analysisTree->refitvertexwithbs_tauIndex[i][1]))
    {
      vtx_x = analysisTree->refitvertexwithbs_x[i];
      vtx_y = analysisTree->refitvertexwithbs_y[i];
      vtx_z = analysisTree->refitvertexwithbs_z[i];
			is_refitted_PV_with_BS = true;
    }
	}
	TVector3 vertex_coord(vtx_x, vtx_y, vtx_z);
	return vertex_coord;
}

TLorentzVector charged_constituents_P4(const AC1B * analysisTree, int tauIndex){
	int ncomponents = analysisTree->tau_constituents_count[tauIndex];
	TLorentzVector sum_charged_P4; sum_charged_P4.SetXYZT(0.,0.,0.,0.);
	TLorentzVector constituent_P4;   
	
	for(int i = 0; i < ncomponents; i++){
		if(analysisTree->tau_constituents_charge[tauIndex][i] != 0 && analysisTree->tau_constituents_pdgId[tauIndex][i] > 0){
			constituent_P4.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][i],
																analysisTree->tau_constituents_py[tauIndex][i],
																analysisTree->tau_constituents_pz[tauIndex][i],
																analysisTree->tau_constituents_e[tauIndex][i]);
			sum_charged_P4 += constituent_P4;
			}
	}
	return sum_charged_P4;
}
