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
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#define PI_MASS         0.13957 //All energies and momenta are expressed in GeV
#define PI0_MASS        0.13498
#define MUON_MASS       0.105658
#define ELECTRON_MASS   0.000511
#define RHO_MASS        0.775260

void acott_Impr_tt(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2);
TLorentzVector IpVec(const AC1B * analysisTree, int Index, int tauIndex, int vertexType);
int ChargedPiIndex(const AC1B * analysisTree, int tauIndex);
TLorentzVector ChargedPivec(const AC1B * analysisTree, int tauIndex);
double AcoCP(TLorentzVector Pi1, TLorentzVector Pi2,TLorentzVector ref1, TLorentzVector ref2,bool firstNegative, bool pi01, bool pi02, Synch17Tree *otree);
TVector3 Get_refitted_PV_with_BS(const AC1B * analysisTree, int firstIndex, int secondIndex,TString ch, bool &is_refitted_PV_with_BS);
TLorentzVector NeutralPivec(const AC1B * analysisTree, int tauIndex);
TLorentzVector IP_helix_Tauh(const AC1B * analysisTree, int tauIndex, TVector3 PV_coord);
double IP_significance_helix_Tauh(const AC1B * analysisTree, int tauIndex, TVector3 PV_coord, const std::vector<float> &PV_cov_components,ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> & ipCovariance, TVector3 & ip);
std::vector<TLorentzVector> A1_rho_pi(const AC1B * analysisTree, int tauIndex);
std::vector<float> Get_refitted_PVBS_cov(const AC1B * analysisTree, int firstIndex, int secondIndex, 
					 TString ch, bool &is_refitted_PV_with_BS);

void acott_Impr_tt(const AC1B * analysisTree, Synch17Tree *otree, int tauIndex1, int tauIndex2){
  otree->acotautau_00 = -9999;
  otree->acotautau_10 = -9999;
  otree->acotautau_01 = -9999;
  otree->acotautau_11 = -9999;

  otree->acotautau_00 = -9999;
  otree->acotautau_10 = -9999;
  otree->acotautau_01 = -9999;
  otree->acotautau_11 = -9999;

  otree->acotautau_helix_00= -9999;
  otree->acotautau_helix_01= -9999;
  otree->acotautau_helix_10= -9999;
  otree->acotautau_helix_11= -9999;

  otree->acotautau_helix_IPIP_pi40_pionA1= -9999;
  otree->acotautau_helix_IPDP_pionA1= -9999;

  otree->acotautau_helix_IPIP_pi40_rhoA1= -9999;
  otree->acotautau_helix_DPIP_pi40_rhoA1= -9999;
  otree->acotautau_helix_DPDP_rhoA1= -9999;

  otree->acotautau_helix_IPIP_pi40_A1A1 = -9999;
  otree->acotautau_helix_IPDP_pi40_A1A1 = -9999;
  otree->acotautau_helix_DPDP_A1A1 = -9999;

  otree->acotautau_helix_IPIP_pi40_pionRho = -9999;
  otree->acotautau_helix_IPDP_pionRho = -9999;
  otree->acotautau_helix_IPIP_pi40_RhoRho  = -9999;
  otree->acotautau_helix_IPDP_pi40_RhoRho = -9999;
  otree->acotautau_helix_DPDP_pionRho = -9999;
  
  otree->ip0x_1 = -9999;
  otree->ip0y_1 = -9999;
  otree->ip0z_1 = -9999;

  otree->ip0x_2 = -9999;
  otree->ip0y_2 = -9999;
  otree->ip0z_2 = -9999;
  otree->ipnx_uncorr_1 = -9999;
  otree->ipny_uncorr_1 = -9999;
  otree->ipnz_uncorr_1 = -9999;

  otree->ipnx_uncorr_2 = -9999;
  otree->ipny_uncorr_2 = -9999;
  otree->ipnz_uncorr_2 = -9999;

  TLorentzVector tau1IP = IpVec(analysisTree,tauIndex2,tauIndex1,2);
  TLorentzVector tau2IP = IpVec(analysisTree,tauIndex1,tauIndex2,2);
  bool is_refitted_PV_with_BS = true;
  TVector3 PV_coord = Get_refitted_PV_with_BS(analysisTree, tauIndex1, tauIndex2, "tt", is_refitted_PV_with_BS);

  std::vector<float> PV_covariance; PV_covariance.clear();
  PV_covariance =  Get_refitted_PVBS_cov(analysisTree, tauIndex1, tauIndex2, "tt", is_refitted_PV_with_BS);


  TLorentzVector tau1IP_helix = IP_helix_Tauh(analysisTree,tauIndex1,PV_coord);
  TLorentzVector tau2IP_helix = IP_helix_Tauh(analysisTree,tauIndex2,PV_coord);
 
  TLorentzVector tau2Prong;
  TLorentzVector tau1Prong;
  TLorentzVector tau1Pi0;
  TLorentzVector tau2Pi0;
  tau1Pi0.SetXYZT(0.,0.,0.,0.); 
  tau2Pi0.SetXYZT(0.,0.,0.,0.); 
  
  tau1Prong=ChargedPivec(analysisTree,tauIndex1);
  tau2Prong=ChargedPivec(analysisTree,tauIndex2);
  int piIndexfortau1=ChargedPiIndex(analysisTree,tauIndex1);
  int piIndexfortau2=ChargedPiIndex(analysisTree,tauIndex2);

  otree->ipnx_uncorr_1 = tau1IP.X()/tau1IP.Mag();
  otree->ipny_uncorr_1 = tau1IP.Y()/tau1IP.Mag();
  otree->ipnz_uncorr_1 = tau1IP.Z()/tau1IP.Mag();

  otree->ipnx_uncorr_2 = tau2IP.X()/tau2IP.Mag();
  otree->ipny_uncorr_2 = tau2IP.Y()/tau2IP.Mag();
  otree->ipnz_uncorr_2 = tau2IP.Z()/tau2IP.Mag();

  otree->ip0x_1 = tau1IP_helix.X();
  otree->ip0y_1 = tau1IP_helix.Y();
  otree->ip0z_1 = tau1IP_helix.Z();
 
  otree->ip0x_2 = tau2IP_helix.X();
  otree->ip0y_2 = tau2IP_helix.Y();
  otree->ip0z_2 = tau2IP_helix.Z();
 

  bool firstNegative = false;
  
  std::vector<TLorentzVector> partLV; partLV.clear();
  if (analysisTree->tau_charge[tauIndex1]>0.0 || analysisTree->tau_charge[tauIndex2]>0.0)
    firstNegative = true;
  
  if((analysisTree->tau_decayMode[tauIndex1]>=1&&analysisTree->tau_decayMode[tauIndex1]<=3)||
     (analysisTree->tau_MVADM2017v1[tauIndex1]>=0.5&&analysisTree->tau_MVADM2017v1[tauIndex1]<=3.5)){

    tau1Prong=ChargedPivec(analysisTree,tauIndex1);   
    tau1Pi0 = NeutralPivec(analysisTree,tauIndex1); 
    
  }
  if (analysisTree->tau_decayMode[tauIndex1]==10||analysisTree->tau_decayMode[tauIndex1]==11||
      (analysisTree->tau_MVADM2017v1[tauIndex1]>=9.5&&analysisTree->tau_MVADM2017v1[tauIndex1]<=11.5)){

    partLV = A1_rho_pi(analysisTree,tauIndex1);
    tau1Prong=partLV.at(0);
    tau1Pi0 = partLV.at(1);
  }
  if((analysisTree->tau_decayMode[tauIndex2]>=1&&analysisTree->tau_decayMode[tauIndex2]<=3)||
     (analysisTree->tau_MVADM2017v1[tauIndex2]>=0.5&&analysisTree->tau_MVADM2017v1[tauIndex2]<=3.5)){
    
    tau2Pi0 = NeutralPivec(analysisTree,tauIndex2); 

  }
  if (analysisTree->tau_decayMode[tauIndex2]==10||analysisTree->tau_decayMode[tauIndex2]==11||
      (analysisTree->tau_MVADM2017v1[tauIndex2]>=9.5&&analysisTree->tau_MVADM2017v1[tauIndex2]<=11.5)){

    partLV = A1_rho_pi(analysisTree,tauIndex2);
    tau2Prong=partLV.at(0);
    tau2Pi0 = partLV.at(1);
  }

 otree->acotautau_00 = AcoCP(tau1Prong,tau2Prong,tau1IP,tau2IP,firstNegative,false,false,otree);
 otree->acotautau_01 = AcoCP(tau1Prong,tau2Prong,tau1IP,tau2Pi0,firstNegative,false,true,otree);
 otree->acotautau_10 = AcoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP,firstNegative,true,false,otree);

 otree->acotautau_helix_00 = AcoCP(tau1Prong,tau2Prong,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);
 otree->acotautau_helix_01 = AcoCP(tau1Prong,tau2Prong,tau1IP_helix,tau2Pi0,firstNegative,false,true,otree);
 otree->acotautau_helix_10 = AcoCP(tau1Prong,tau2Prong,tau1Pi0,tau2IP_helix,firstNegative,true,false,otree);

  ////////////////////////////////////////////////
  //Begin on pion pT studies on A1 variables
  ////////////////////////////////////////////////
  TLorentzVector pion1LV;
  TLorentzVector RhoPion1LV;
  TLorentzVector RhoPion2LV;
  TLorentzVector A1Pion1LV;
  TLorentzVector A1Pion2LV;
  pion1LV.SetXYZT(0.,0.,0.,0.); 
  RhoPion1LV.SetXYZT(0.,0.,0.,0.); 
  RhoPion2LV.SetXYZT(0.,0.,0.,0.); 
  A1Pion1LV.SetXYZT(0.,0.,0.,0.); 
  A1Pion2LV.SetXYZT(0.,0.,0.,0.); 
  ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> ipCov2;
  TVector3 IP2;
  double ip_sig_tau1 = IP_significance_helix_Tauh(analysisTree,tauIndex1,PV_coord,PV_covariance,ipCov2,IP2);
  double ip_sig_tau2 = IP_significance_helix_Tauh(analysisTree,tauIndex2,PV_coord,PV_covariance,ipCov2,IP2);
  if(analysisTree->tau_decayMode[tauIndex1]==0 && analysisTree->tau_decayMode[tauIndex2]==2){

    pion1LV   = ChargedPivec(analysisTree,tauIndex1);
    A1Pion1LV = ChargedPivec(analysisTree,tauIndex2);
    double IPIP_pionA1 = AcoCP(pion1LV,A1Pion1LV,tau1IP_helix,tau2Pi0,firstNegative,false,true,otree);
    if(A1Pion1LV.Pt() >= 40 && ip_sig_tau1>1.5 && ip_sig_tau2>1.5){
      otree->acotautau_helix_IPIP_pi40_pionA1 = AcoCP(pion1LV,A1Pion1LV,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);
    }
    else{
      otree->acotautau_helix_IPDP_pionA1 = IPIP_pionA1;
    }
  }
 if(analysisTree->tau_decayMode[tauIndex1]==2 && analysisTree->tau_decayMode[tauIndex2]==0){
    pion1LV   = ChargedPivec(analysisTree,tauIndex2);
    A1Pion1LV = ChargedPivec(analysisTree,tauIndex1);
    double IPIP_pionA1 = AcoCP(pion1LV,A1Pion1LV,tau2Pi0,tau1IP_helix,firstNegative,true,false,otree);
    if(A1Pion1LV.Pt() >= 40 && ip_sig_tau2>1.5 && ip_sig_tau1>1.5){
      otree->acotautau_helix_IPIP_pi40_pionA1 = AcoCP(pion1LV,A1Pion1LV,tau2IP_helix,tau1IP_helix,firstNegative,false,false,otree);}
    else{
      otree->acotautau_helix_IPDP_pionA1 = IPIP_pionA1;
    }
 }    
 if(analysisTree->tau_decayMode[tauIndex1]==1 && analysisTree->tau_decayMode[tauIndex2]==2){
   RhoPion1LV   = ChargedPivec(analysisTree,tauIndex1);
   A1Pion1LV    = ChargedPivec(analysisTree,tauIndex2);
   double IPIP_rhoA1 = AcoCP(RhoPion1LV,A1Pion1LV,tau1Pi0,tau2Pi0,firstNegative,true,true,otree);
   if(A1Pion1LV.Pt() >= 40 && ip_sig_tau2>1.5){
     otree->acotautau_helix_DPIP_pi40_rhoA1 = AcoCP(RhoPion1LV,A1Pion1LV,tau1Pi0,tau2IP_helix,firstNegative,true,false,otree);}
   else{
     otree->acotautau_helix_DPDP_rhoA1 = IPIP_rhoA1;
   }
 }
 if(analysisTree->tau_decayMode[tauIndex1]==2 && analysisTree->tau_decayMode[tauIndex2]==1){
   RhoPion1LV   = ChargedPivec(analysisTree,tauIndex2);
   A1Pion1LV    = ChargedPivec(analysisTree,tauIndex1);
   double IPIP_rhoA1 = AcoCP(RhoPion1LV,A1Pion1LV,tau1Pi0,tau2Pi0,firstNegative,true,true,otree);
   if(A1Pion1LV.Pt() >= 40 && ip_sig_tau1>1.5){
     otree->acotautau_helix_IPIP_pi40_rhoA1 = AcoCP(RhoPion1LV,A1Pion1LV,tau2IP_helix,tau1IP_helix,firstNegative,false,false,otree);}
   else{
     otree->acotautau_helix_DPDP_rhoA1 = IPIP_rhoA1;
   }
 }
if(analysisTree->tau_decayMode[tauIndex1]==2 && analysisTree->tau_decayMode[tauIndex2]==2){
   A1Pion1LV   = ChargedPivec(analysisTree,tauIndex1);
   A1Pion2LV    = ChargedPivec(analysisTree,tauIndex2);
   double  IPIP_A1A1 = AcoCP(A1Pion1LV,A1Pion2LV,tau1Pi0,tau2Pi0,firstNegative,true,true,otree);
   if(A1Pion1LV.Pt() >= 40 && A1Pion2LV.Pt()>=40 && ip_sig_tau1>1.5 && ip_sig_tau2>1.5){
     otree->acotautau_helix_IPIP_pi40_A1A1 = AcoCP(A1Pion1LV,A1Pion2LV,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);}
   else if(A1Pion1LV.Pt() >= 40 && A1Pion2LV.Pt()<40 && ip_sig_tau1>1.5){
     otree->acotautau_helix_IPDP_pi40_A1A1 = AcoCP(A1Pion1LV,A1Pion2LV,tau1IP_helix,tau2Pi0,firstNegative,false,true,otree);}
   else if(A1Pion1LV.Pt() < 40 && A1Pion2LV.Pt()>=40 && ip_sig_tau2>1.5){
     otree->acotautau_helix_IPDP_pi40_A1A1 = AcoCP(A1Pion1LV,A1Pion2LV,tau1Pi0,tau2IP_helix,firstNegative,true,false,otree);}
   else{
     otree->acotautau_helix_DPDP_A1A1 = IPIP_A1A1;
   }
 }
if(analysisTree->tau_decayMode[tauIndex1]==0 && analysisTree->tau_decayMode[tauIndex2]==1){
  pion1LV = ChargedPivec(analysisTree,tauIndex1);
  RhoPion1LV = ChargedPivec(analysisTree,tauIndex2);
  double IPIP_pionA1 = AcoCP(pion1LV,RhoPion1LV,tau1IP_helix,tau2Pi0,firstNegative,false,true,otree);
  if(RhoPion1LV.Pt() >=40 && ip_sig_tau1>1.5&& ip_sig_tau2>1.5){
    otree->acotautau_helix_IPIP_pi40_pionRho = AcoCP(pion1LV,RhoPion1LV,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);}
  else{
    otree->acotautau_helix_IPDP_pionRho = IPIP_pionA1;}
 }
if(analysisTree->tau_decayMode[tauIndex1]==1 && analysisTree->tau_decayMode[tauIndex2]==0){
  pion1LV = ChargedPivec(analysisTree,tauIndex2);
  RhoPion1LV = ChargedPivec(analysisTree,tauIndex1);
  double IPIP_pionA1 = AcoCP(pion1LV,RhoPion1LV,tau1Pi0,tau2IP_helix,firstNegative,true,false,otree);
  if(RhoPion1LV.Pt() >=40 && ip_sig_tau1>1.5&& ip_sig_tau2>1.5){
    otree->acotautau_helix_IPIP_pi40_pionRho = AcoCP(pion1LV,RhoPion1LV,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);}
  else{
    otree->acotautau_helix_IPDP_pionRho = IPIP_pionA1;}
 }
if(analysisTree->tau_decayMode[tauIndex1]==1 && analysisTree->tau_decayMode[tauIndex2]==1){
  RhoPion1LV = ChargedPivec(analysisTree,tauIndex1);
  RhoPion2LV = ChargedPivec(analysisTree,tauIndex2);
  double IPIP_pionA1 = AcoCP(RhoPion1LV,RhoPion2LV,tau1Pi0,tau2Pi0,firstNegative,true,true,otree);
  if(RhoPion1LV.Pt() >=40 && RhoPion2LV.Pt() >=40 && ip_sig_tau1>1.5&& ip_sig_tau2>1.5){
    otree->acotautau_helix_IPIP_pi40_RhoRho = AcoCP(RhoPion1LV,RhoPion2LV,tau1IP_helix,tau2IP_helix,firstNegative,false,false,otree);}
  else if(RhoPion1LV.Pt() >=40 && RhoPion2LV.Pt() <40 && ip_sig_tau1>1.5){
    otree->acotautau_helix_IPDP_pi40_RhoRho = AcoCP(RhoPion1LV,RhoPion2LV,tau1IP_helix,tau2Pi0,firstNegative,false,true,otree);}
  else if(RhoPion1LV.Pt() <40 && RhoPion2LV.Pt() >=40 && ip_sig_tau2>1.5){
    otree->acotautau_helix_IPDP_pi40_RhoRho = AcoCP(RhoPion1LV,RhoPion2LV,tau1Pi0,tau2IP_helix,firstNegative,true,false,otree);}
  else{
    otree->acotautau_helix_DPDP_pionRho = IPIP_pionA1;}
 }
///////////////////////////////////////////
 

}
TLorentzVector IP_helix_Tauh(const AC1B * analysisTree, int tauIndex, TVector3 PV_coord){
  // helical IP for tau_h
  // NB: for calculation will take the 4-momentum of the leading charged hadron, same sign as tau
  TLorentzVector LVIP={0.,0.,0.,0.};

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
  ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> PV(PV_coord.X(), PV_coord.Y(), PV_coord.Z());
  
  ImpactParameter IP;
  //TVector3 IP_helix_tau = IP.CalculatePCA(B, h_param_tau, ref_tau, pv, p4_tau);     //kept for retrocompatibility           
  TVector3 IP_helix_tau = IP.CalculatePCA(B, h_param_tau, ref_tau, PV);              
  LVIP.SetVect(IP_helix_tau); 
  return LVIP;
};

TLorentzVector IpVec(const AC1B * analysisTree, int Index, int tauIndex, int vertexType) {
  TLorentzVector vec;
  vec.SetXYZT(0.,0.,0.,0.);
  int piIndex=-1;
  piIndex=ChargedPiIndex(analysisTree,tauIndex);

  if (piIndex>-1) {

    TVector3 vertex(analysisTree->primvertex_x,
		    analysisTree->primvertex_y,
		    analysisTree->primvertex_z);
    
    if (vertexType==1) {
      vertex.SetX(analysisTree->primvertexwithbs_x);
      vertex.SetY(analysisTree->primvertexwithbs_y);
      vertex.SetZ(analysisTree->primvertexwithbs_z);
    }
    else if (vertexType==2) {
      bool refitted_PV_with_BS = false;
      vertex = Get_refitted_PV_with_BS(analysisTree, Index, tauIndex, "tt", refitted_PV_with_BS);
    }
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
    
}
TLorentzVector ChargedPivec(const AC1B * analysisTree, int tauIndex){
  int piIndex = -1;
  piIndex = ChargedPiIndex(analysisTree, tauIndex);
  TLorentzVector chargedPi;
  chargedPi.SetXYZT(0,0,0,0);

  if (piIndex > -1) {
    chargedPi.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][piIndex],
			 analysisTree->tau_constituents_py[tauIndex][piIndex],
			 analysisTree->tau_constituents_pz[tauIndex][piIndex],
			 analysisTree->tau_constituents_e[tauIndex][piIndex]);
  }  
  return chargedPi;
};

int ChargedPiIndex(const AC1B * analysisTree, int tauIndex){
  // selects the highest energy Pi with the same sign of the tau
  
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex = -1;
  float maxPt = -1;
  int sign = -1;
  if(analysisTree->tau_charge[tauIndex] > 0) sign = 1; 
  for(int i = 0; i < ncomponents; i++){ 
    if((analysisTree->tau_constituents_pdgId[tauIndex][i]*sign) == 321 ||
       (analysisTree->tau_constituents_pdgId[tauIndex][i]*sign) == 211) {
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
  if (piIndex>=0) return piIndex;

  for(int i = 0; i < ncomponents; i++){ 
    if((analysisTree->tau_constituents_pdgId[tauIndex][i]*sign) == -13 || 
       (analysisTree->tau_constituents_pdgId[tauIndex][i]*sign) == -11) {
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
  if (piIndex>=0) return piIndex;

  for(int i = 0; i < ncomponents; i++){ 
    if(TMath::Abs(analysisTree->tau_constituents_pdgId[tauIndex][i]) == 321 ||
       TMath::Abs(analysisTree->tau_constituents_pdgId[tauIndex][i]) == 211) { 
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
  if (piIndex>=0) return piIndex;

  for(int i = 0; i < ncomponents; i++){ 
    if(TMath::Abs(analysisTree->tau_constituents_pdgId[tauIndex][i]) == 13 || 
       TMath::Abs(analysisTree->tau_constituents_pdgId[tauIndex][i]) == 11){
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

  if (piIndex<0) {
    std::cout << "EVT = " << analysisTree->event_nr << std::endl;
    std::cout << "tau charge = " << analysisTree->tau_charge[tauIndex] 
            << "    dmHPS = " << analysisTree->tau_decayMode[tauIndex] 
	      << "    dmMVA = " << analysisTree->tau_MVADM2017v1[tauIndex] << std::endl;
    for(int i = 0; i < ncomponents; i++){ 
      std::cout << "constituent " << i << "  pdgId = " << analysisTree->tau_constituents_pdgId[tauIndex][i] << "   q = " << analysisTree->tau_constituents_charge[tauIndex][i] << std::endl;
    }
    std::cout << std::endl;
  }

  return piIndex;
};
TLorentzVector NeutralPivec(const AC1B * analysisTree, int tauIndex){
  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  int piIndex=-1;
  TLorentzVector neutralPi; neutralPi.SetXYZT(0.,0.,0.,0.);
  double etot = 0.;
  double emax = 0.;
  float maxPT=-100;int leadingIndex;float sumPT =0,sumE=0;
  TLorentzVector neutralpart; neutralpart.SetXYZT(0.,0.,0.,0.);     //momentum of the leading particle
  for(int i=0;i<ncomponents;i++){
    if(analysisTree->tau_constituents_pdgId[tauIndex][i]==22||abs(analysisTree->tau_constituents_pdgId[tauIndex][i])==11){
                                               
      etot += analysisTree->tau_constituents_e[tauIndex][i];
      if(analysisTree->tau_constituents_e[tauIndex][i] > emax){
	neutralpart.SetPxPyPzE(analysisTree->tau_constituents_px[tauIndex][i],
			       analysisTree->tau_constituents_py[tauIndex][i],
			       analysisTree->tau_constituents_pz[tauIndex][i],
			       analysisTree->tau_constituents_e[tauIndex][i]);
	
        maxPT = neutralpart.Pt();
	emax = analysisTree->tau_constituents_e[tauIndex][i];
        leadingIndex = i;
      }
    }
  }
  if (etot>PI0_MASS) {
    double mom = neutralpart.P();
    double momtot = TMath::Sqrt(etot*etot-PI0_MASS*PI0_MASS);
    double scale = momtot/mom;
    double px = scale*neutralpart.Px();
    double py = scale*neutralpart.Py();
    double pz = scale*neutralpart.Pz();
    neutralpart.SetPxPyPzE(px,py,pz,etot);
  }
  else {
    double mom = neutralpart.P();
    double scale = etot/mom;
    double px = scale*neutralpart.Px();
    double py = scale*neutralpart.Py();
    double pz = scale*neutralpart.Pz();
    neutralpart.SetPxPyPzE(px,py,pz,etot);    
  }

  return neutralpart;
  
};

TVector3 Get_refitted_PV_with_BS(const AC1B * analysisTree, int firstIndex, int secondIndex, 
				 TString ch, bool &is_refitted_PV_with_BS) {
   float vtx_x = analysisTree->primvertex_x; // by default store non-refitted PV with BS constraint if refitted one is not found
   float vtx_y = analysisTree->primvertex_y;
   float vtx_z = analysisTree->primvertex_z;
   is_refitted_PV_with_BS = false;
   std::vector<int> first_indices(2, -787);

   for(unsigned int i = 0; i < analysisTree->refitvertexwithbs_count; i++)
     {    
       if (ch == "mt")
	 {
	   first_indices[0] = analysisTree->refitvertexwithbs_muIndex[i][0];
	   first_indices[1] = analysisTree->refitvertexwithbs_muIndex[i][1];
	 }
       else if (ch == "et")
	 {
	   first_indices[0] = analysisTree->refitvertexwithbs_eleIndex[i][0];
	   first_indices[1] = analysisTree->refitvertexwithbs_eleIndex[i][1];
	 }
       else if (ch == "tt")
	 {
	   first_indices[0] = analysisTree->refitvertexwithbs_tauIndex[i][0];
	   first_indices[1] = analysisTree->refitvertexwithbs_tauIndex[i][1];
	 }
        
       // secondIndex is assumed to be tau_h
       if( (firstIndex == first_indices[0] || firstIndex == first_indices[1]) &&
	   (secondIndex == analysisTree->refitvertexwithbs_tauIndex[i][0] || secondIndex == analysisTree->refitvertexwithbs_tauIndex[i][1]))
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
std::vector<float> Get_refitted_PVBS_cov(const AC1B * analysisTree, int firstIndex, int secondIndex, 
					 TString ch, bool &is_refitted_PV_with_BS) {
  
  std::vector<float> PV_covariance; PV_covariance.clear();
  float vtx_x = analysisTree->primvertex_x; // by default store non-refitted PV with BS constraint if refitted one is not found
  float vtx_y = analysisTree->primvertex_y;
  float vtx_z = analysisTree->primvertex_z;
  is_refitted_PV_with_BS = false;
  std::vector<int> first_indices(2, -787);

  for(unsigned int i = 0; i < analysisTree->refitvertexwithbs_count; i++)
    {    
      if (ch == "mt")
	{
	  first_indices[0] = analysisTree->refitvertexwithbs_muIndex[i][0];
	  first_indices[1] = analysisTree->refitvertexwithbs_muIndex[i][1];
	}
      else if (ch == "et")
	{
	  first_indices[0] = analysisTree->refitvertexwithbs_eleIndex[i][0];
	  first_indices[1] = analysisTree->refitvertexwithbs_eleIndex[i][1];
	}
      else if (ch == "tt")
	{
	  first_indices[0] = analysisTree->refitvertexwithbs_tauIndex[i][0];
	  first_indices[1] = analysisTree->refitvertexwithbs_tauIndex[i][1];
	}
        
      // secondIndex is assumed to be tau_h
      if( (firstIndex == first_indices[0] || firstIndex == first_indices[1]) &&
	  (secondIndex == analysisTree->refitvertexwithbs_tauIndex[i][0] || secondIndex == analysisTree->refitvertexwithbs_tauIndex[i][1]))
	{
	  vtx_x = analysisTree->refitvertexwithbs_x[i];
	  vtx_y = analysisTree->refitvertexwithbs_y[i];
	  vtx_z = analysisTree->refitvertexwithbs_z[i];
	  is_refitted_PV_with_BS = true;
	  for (int j=0; j<6; ++j) 
	    PV_covariance[j] = analysisTree->refitvertexwithbs_cov[i][j];
	
	}
    }
 
  return PV_covariance;
}
double AcoCP(TLorentzVector Pi1, TLorentzVector Pi2, 
	     TLorentzVector ref1, TLorentzVector ref2,
	     bool firstNegative, bool pi01, bool pi02, Synch17Tree *otree) {

  double y1 = 1;
  double y2 = 1;
  if (pi01)
    y1 = Pi1.E() - ref1.E();
  if (pi02)
    y2 = Pi2.E() - ref2.E();

  otree->y1_LF=(Pi1.E() - ref1.E())/(Pi1.E() + ref1.E());
  otree->y2_LF=(Pi2.E() - ref2.E())/(Pi2.E() + ref2.E());

  double y = y1*y2; 
  //double y = -1; 

  //save copies before boost in charged zmf:
  TLorentzVector Pi1Lab=Pi1;
  TLorentzVector Pi2Lab=Pi2;
  TLorentzVector ref1Lab=ref1;
  TLorentzVector ref2Lab=ref2;

  TLorentzVector Prongsum = Pi1 + Pi2;
  TVector3 boost = -Prongsum.BoostVector();
  Pi1.Boost(boost);
  Pi2.Boost(boost);
  ref1.Boost(boost);
  ref2.Boost(boost);

  otree->y1_ZMF=(Pi1.E() - ref1.E())/(Pi1.E() + ref1.E());
  otree->y2_ZMF=(Pi2.E() - ref2.E())/(Pi2.E() + ref2.E());

  // 
  // if (pi01&&pi02){

  TLorentzVector Rho1, Rho2, Rho1Boost, Rho2Boost;
  Rho1=Pi1+ref1;
  Rho2=Pi2+ref2;

  TLorentzVector ProngsumRho = Rho1 + Rho2;
  TVector3 boostrho = -ProngsumRho.BoostVector();
  Rho1.Boost(boostrho);
  Rho2.Boost(boostrho);

  TVector3 vecRho1 = Rho1.Vect();
  TVector3 vecRho2 = Rho2.Vect();
  double vecRho1Mag=vecRho1.Mag();
  double vecRho2Mag=vecRho2.Mag();  
  double CRho1, CRho2;
  double mH=125.18;
  double mtau=1.776;

  if(vecRho1Mag!=0&&vecRho2Mag!=0){
    CRho1=TMath::Sqrt( (TMath::Power(0.5*mH,2)-TMath::Power(mtau,2))/TMath::Power(vecRho1Mag,2) );
    CRho2=TMath::Sqrt( (TMath::Power(0.5*mH,2)-TMath::Power(mtau,2))/TMath::Power(vecRho2Mag,2) );
    vecRho1*=CRho1;
    vecRho2*=CRho2;
    //define new 4 vectors
    // TLorentzVector P1(vecRho1[0],vecRho1[1],vecRho1[2],(0.5*mH)); 
    TLorentzVector P1(vecRho1,(0.5*mH));
    TLorentzVector P2(vecRho2,(0.5*mH));

    //  cout<<" cros check P1.Mag() "<<P1.Mag() <<endl; 
    //  cout<<" cros check P2.Mag() "<<P2.Mag() <<endl;
    TVector3 boostvec1=-P1.BoostVector();
    TVector3 boostvec2=-P2.BoostVector();

    Pi1Lab.Boost(boostvec1);
    ref1Lab.Boost(boostvec1);
    Pi2Lab.Boost(boostvec2);
    ref2Lab.Boost(boostvec2);

    otree->y1_TMF=(Pi1Lab.E() - ref1Lab.E())/(Pi1Lab.E() + ref1Lab.E());
    otree->y2_TMF=(Pi2Lab.E() - ref2Lab.E())/(Pi2Lab.E() + ref2Lab.E());

  }
  //  cout<<"y2_T1F "<<y2_T1F<<endl;
  //  cout<<"y2_ZMF "<<otree->y2_ZMF<<endl;

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
double IP_significance_helix_Tauh(const AC1B * analysisTree, int tauIndex, TVector3 PV_coord, const std::vector<float> &PV_cov_components, ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >> & ipCovariance, TVector3 & ip)
{
  ImpactParameter IP;
  std::pair <TVector3, ROOT::Math::SMatrix<float,3,3, ROOT::Math::MatRepStd< float, 3, 3 >>> ipAndCov;
  std::vector<float> h_param_tau = {};
  RMPoint ref_tau;
  SMatrixSym3D PV_covariance;
  SMatrixSym5F helix_params_covariance;
  
  int k = 0;
  double B = analysisTree->tau_Bfield[tauIndex];
  ref_tau.SetX(analysisTree->tau_referencePoint[tauIndex][0]);
  ref_tau.SetY(analysisTree->tau_referencePoint[tauIndex][1]);
  ref_tau.SetZ(analysisTree->tau_referencePoint[tauIndex][2]);
  RMPoint PV(PV_coord.X(), PV_coord.Y(), PV_coord.Z());
  for(auto i:  analysisTree->tau_helixparameters[tauIndex]) h_param_tau.push_back(i);
  
  // !! check that in NTupleMaker the logic of filling PV_cov_components is the same 
  // for more on how to fill SMatrices see: https://root.cern/doc/master/SMatrixDoc.html 
  for (size_t i = 0; i < 5; i++)
    for (size_t j = i; j < 5; j++) // should be symmetrically completed automatically
      helix_params_covariance[i][j] = analysisTree->tau_helixparameters_covar[tauIndex][i][j];
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = i; j < 3; j++) // should be symmetrically completed automatically
	{
	  PV_covariance[i][j] = PV_cov_components[k];
	  k++;
	}
    }
  
  ipAndCov = IP.CalculateIPandCovariance(
					 B, // (double)
					 h_param_tau, // (std::vector<float>)
					 ref_tau, // (RMPoint)
					 PV, // (RMPoint)
					 helix_params_covariance, // (ROOT::Math::SMatrix<float,5,5, ROOT::Math::MatRepSym<float,5>>)
					 PV_covariance // (SMatrixSym3D)
					 );
  
  ip = ipAndCov.first;
  ipCovariance = ipAndCov.second;
  double ipSignificance = IP.CalculateIPSignificanceHelical(ip, ipCovariance);
  
  return ipSignificance;
}
std::vector<TLorentzVector> A1_rho_pi(const AC1B * analysisTree, int tauIndex) {

  int ncomponents = analysisTree->tau_constituents_count[tauIndex];
  std::vector<TLorentzVector> partLV; partLV.clear();
  TLorentzVector lv_rho; lv_rho.SetXYZT(0.,0.,0.,0.);

  int index_pi_opposite = -1; 
  TLorentzVector lv_pi_opposite; 
  lv_pi_opposite.SetXYZT(0.,0.,0.,0.);

  int index_pi_from_rho = -1;
  TLorentzVector lv_pi_from_rho;
  lv_pi_from_rho.SetXYZT(0.,0.,0.,0.);

  int index_pi_lead = -1;
  TLorentzVector lv_pi_lead;
  lv_pi_lead.SetXYZT(0.,0.,0.,0.);

  int tau_charge = -1; 
  if (analysisTree->tau_charge[tauIndex]>0) tau_charge = 1;

  for(int i = 0; i < ncomponents; i++){
    if(analysisTree->tau_constituents_charge[tauIndex][i] != 0 && abs(analysisTree->tau_constituents_pdgId[tauIndex][i])==211){
      if (tau_charge*analysisTree->tau_constituents_charge[tauIndex][i]<0) { 
	lv_pi_opposite.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
			       analysisTree->tau_constituents_py[tauIndex][i],
			       analysisTree->tau_constituents_pz[tauIndex][i],
			       analysisTree->tau_constituents_e[tauIndex][i]);
	index_pi_opposite = i;
	//break;
      }
    }
  }


  double massdiff = 1e+10;
  for (int i = 0; i < ncomponents; i++) {
    if(analysisTree->tau_constituents_charge[tauIndex][i] != 0 && abs(analysisTree->tau_constituents_pdgId[tauIndex][i])==211){
      if (tau_charge*analysisTree->tau_constituents_charge[tauIndex][i]>0) {
        TLorentzVector lv; lv.SetXYZT(analysisTree->tau_constituents_px[tauIndex][i],
				      analysisTree->tau_constituents_py[tauIndex][i],
				      analysisTree->tau_constituents_pz[tauIndex][i],
				      analysisTree->tau_constituents_e[tauIndex][i]);
	double mass = (lv_pi_opposite+lv).M();
	double diff = fabs(mass-RHO_MASS);
	if (diff<massdiff) {
	  index_pi_from_rho = i;
	  lv_pi_from_rho = lv;
	  massdiff = diff;
	}  
      }
    }
  }

  partLV.push_back(lv_pi_from_rho);
  partLV.push_back(lv_pi_opposite);

  return partLV;

}
