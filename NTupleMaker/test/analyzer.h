//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 23 09:55:01 2016 by ROOT version 6.02/05
// from TChain eltau/T/
//////////////////////////////////////////////////////////

#ifndef analyzer_h
#define analyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "plots.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <map>
#include <cmath>
#include <string>


//using namespace Mt2;
// Header file for the classes stored in the TTree if any.

class analyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  vector<TLorentzVector> AllJets_Lepton_noMet;
  vector<TLorentzVector> JetsMV;
  vector<TLorentzVector>  ElMV;
  vector<TLorentzVector>  MuMV;
  vector<TLorentzVector>  TauMV;
  vector<TLorentzVector>  LeptMV;
  TLorentzVector  ElV;
  TLorentzVector  MuV;
  TLorentzVector  TauV;
  TLorentzVector  JetsV;
  TLorentzVector  LeptV1;
  TLorentzVector  LeptV2;
  TLorentzVector  METV;

  vector <int> btag_index;


  const double MuMass = 0.105658367;
  double tauMass = 1.776;
  const double electronMass = 0.51100e-3;
  const double muonMass = 0.10565837;
  const double pionMass = 0.1396;
  const double bTagCut  = 0.8000;

    double dzeta=-9999;
    double pzetavis=-9999;
    double pzetamiss=-9999;
    double Mt2as=-1;
    double met=0;
    double metx=0;
    double mety=0;
    double metphi=0;
   Float_t         wScale0;
   Float_t         wScale1;
   Float_t         wScale2;
   Float_t         wScale3;
   Float_t         wScale4;
   Float_t         wScale5;
   Float_t         wScale6;
   Float_t         wScale7;
   Float_t         wScale8;
   Float_t         wPDFUp;
   Float_t         wPDFDown;

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types

  Int_t	   nbtag;
  Int_t	   njets;
   int 	   	   event_lumi;
   int 	   	   event_run;
   Float_t         genmet_Ex;
   Float_t         genmet_Ey;
   Float_t         NuPx;
   Float_t         NuPy;
   Float_t         NuPz;
   Float_t         NuPt;
   Float_t         NuPhi;
  Float_t         met_ex;
  Float_t         met_ey;
  Float_t         met_ez;
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_ex_JetEnUp;
  Float_t         met_ey_JetEnUp;
  Float_t         met_ex_JetEnDown;
  Float_t         met_ey_JetEnDown;
  Float_t		met_ex_UnclusteredEnUp;
  Float_t		met_ey_UnclusteredEnUp;
  Float_t		met_ex_UnclusteredEnDown;
  Float_t		met_ey_UnclusteredEnDown;

  Float_t         met_ex_JetEnUp_recoil;
  Float_t         met_ey_JetEnUp_recoil;
  Float_t         met_ex_JetEnDown_recoil;
  Float_t         met_ey_JetEnDown_recoil;
  Float_t         met_ex_UnclusteredEnUp_recoil;
  Float_t         met_ey_UnclusteredEnUp_recoil;
  Float_t         met_ex_UnclusteredEnDown_recoil;
  Float_t         met_ey_UnclusteredEnDown_recoil;
  Float_t         met_ex_recoil;
  Float_t         met_ey_recoil;
  Float_t         genmet;
  Float_t         genmetphi;
  Float_t         met_scaleUp;
  Float_t         metphi_scaleUp;
  Float_t         met_scaleDown;
  Float_t         metphi_scaleDown;
  Float_t         met_resoUp;
  Float_t         met_resoDown;
  Float_t         metphi_resoUp;
  Float_t         metphi_resoDown;

  Float_t         gen_weight;
  Float_t         pu_weight;
  Float_t         LSF_weight;
  Float_t         LSF_weight_el;
  Float_t         LSF_weight_mu;
  Float_t         LSF_weight_2;
  Float_t         LSF_weight_1;
  Float_t         TFR_weight;
  Float_t         top_weight;
  Float_t         all_weight;
  Float_t         trig_weight;
  Float_t         trig_weight_1;
  Float_t         trig_weight_2;
  Float_t         zptmassweight;
  Float_t         xsecs;
  Int_t           muon_index;
  Int_t           muon_index_1;
  Int_t           muon_index_2;
  Int_t           electron_index;
  Int_t           taus_index;
  Int_t           primvert_count;
  Float_t         primvert_x;
  Float_t         primvert_y;
  Float_t         primvert_z;
  Int_t           mu_count;
  Float_t         mu_px[20];
  Float_t         mu_py[20];
  Float_t         mu_pz[20];
  Float_t         mu_pt[20];
  Float_t         mu_eta[20];
  Float_t         mu_phi[20];
  Float_t         mu_charge[20];
  Float_t         mu_miniISO[20];
  Float_t         mu_dxy[20];
  Float_t         mu_dz[20];
  Float_t         mu_dzerr[20];
  Float_t         mu_dxyerr[20];
  Float_t         mu_relIso[20];
		
  Float_t     mu_neutralHadIso[20]; 
  Float_t     mu_photonIso[20]; 
  Float_t     mu_chargedHadIso[20]; 
  Float_t     mu_puIso[20]; 
  Float_t     mu_neutralIso[20];
  Float_t     mu_absIsoMu[20]; 
  Float_t     mu_relIsoMu[20]; 

  Float_t     el_neutralHadIso[20]; 
  Float_t     el_photonIso[20]; 
  Float_t     el_chargedHadIso[20]; 
  Float_t     el_puIso[20]; 
  Float_t     el_neutralIso[20];
  Float_t     el_absIsoEl[20]; 
  Float_t     el_relIsoEl[20]; 
  Float_t     el_isMVA[20]; 
               
  Int_t           jet_count;
  Int_t           npv;
  int           genTauDecayMode1;
  int           genTauDecayMode2;
  int           genTauDecayModeAll;
  Int_t           npu;
  Int_t           jets_cleaned[30];
  Float_t         jet_e[30];
  Float_t         jet_px[30];
  Float_t         jet_py[30];
  Float_t         jet_pz[30];
  Float_t         jet_pt[30];
  Float_t         jet_eta[30];
  Float_t         jet_phi[30];
  Float_t         jet_flavour[30];
  Float_t         jet_btag[30];
  Float_t         CFCounter_[30];
  Int_t           jet_isLoose[30];
  Int_t           el_count;
  Float_t         el_px[20];
  Float_t         el_py[20];
  Float_t         el_pz[20];
  Float_t         el_pt[20];
  Float_t         el_eta[20];
  Float_t         el_phi[20];
  Float_t         el_miniISO[20];
  Float_t         el_dxy[20];
  Float_t         el_dxyerr[20];
  Float_t         el_dz[20];
  Float_t         el_dzerr[20];
  Float_t         el_charge[20];
  Float_t         el_relIso[20];
  Int_t           ta_count;
  Float_t         ta_px[30];
  Float_t         ta_py[30];
  Float_t         ta_pz[30];
  Float_t         ta_mass[30];
  Float_t         ta_eta[30];
  Float_t         ta_phi[30];
  Float_t         ta_pt[30];
  Float_t         ta_dxy[30];
  Float_t         ta_dz[30];
  Float_t         ta_charge[30];
  Float_t         ta_relIso[30];
  Float_t         ta_IsoFlag;
  Float_t         event_sign;
  Float_t         met_flag;
  Float_t         event_secondLeptonVeto;
  Float_t         eleMVA;
  Float_t         event_thirdLeptonVeto;
  Float_t         event_leptonDrTrigger;
  //   string 	   datasetName;
  string * datasetName = new std::string(); 
  string          *regionName;

  Float_t 	genTauMatched;
  Float_t 	genLeptonMatched;
  Float_t 	genLeptonMatchedEl;
  Float_t 	genLeptonMatchedMu;
  Float_t 	genTauDecayedMuMatched;
  Float_t 	genTauDecayedElMatched;
  Float_t 	genLeptonPromptElMatched;
  Float_t 	genLeptonPromptMuMatched;
  Float_t         genLeptonMatchedPrompEl;
  Float_t         genLeptonMatchedPrompMu;
  Float_t         genElMatchedToTauDecay;
  Float_t         genMuMatchedToTauDecay;
  Float_t         genTauMatchedToTauDecay;
  Float_t	        genLeptonMatchedPromptEl;
  Float_t	        genLeptonMatchedPromptMu;
  Float_t	        genLeptonMatchedPromptTau;
  Float_t	   genElMatchedHadrDecay;
  Float_t	   genMuMatchedHadrDecay;
  Float_t	   genTauMatchedHadrDecay;
  Float_t	   genLeptonMatchedGluon;
  Float_t	   genLeptonMatchedHFQ;
  Float_t	   genLeptonMatchedLFQ;

  Float_t         matchedTauToPromptEl;
  Float_t         matchedTauToPromptMu;
  Float_t         matchedTauToPromptTau;
  Float_t         matchedTauToTauDecEl;
  Float_t         matchedTauToTauDecMu;
  Float_t         matchedTauToTauDecTau;

  Float_t	   matchedTauToElHadronDec;
  Float_t	   matchedTauToMuHadronDec;
  Float_t	   matchedTauToTauHadronDec;
   Float_t	   matchedTauToGluon;
   Float_t	   matchedTauToHFQ;
   Float_t	   matchedTauToLFQ;
  Float_t 	isDYTT;
  Float_t 	isDYEE;
  Float_t 	isDYMM;
  Float_t 	isDYNuNu;

  Float_t	   qcdweight;
  Float_t	   qcdweightup;
  Float_t	   qcdweightdown;
  Int_t 	   npartons;

  // List of branches
  TBranch        *b_met_ex;   //!
  TBranch        *b_met_ey;   //!
  TBranch        *b_met_ez;   //!
  TBranch        *b_met_pt;   //!
  TBranch        *b_met_phi;   //!
  TBranch        *b_met_ex_JetEnUp;
  TBranch        *b_met_ey_JetEnUp;
  TBranch        *b_met_ex_JetEnDown;
  TBranch        *b_met_ey_JetEnDown;
  TBranch	       *b_met_ex_UnclusteredEnUp;
  TBranch	       *b_met_ey_UnclusteredEnUp;
  TBranch	       *b_met_ex_UnclusteredEnDown;
  TBranch	       *b_met_ey_UnclusteredEnDown;
  TBranch        *b_met_ex_JetEnUp_recoil;
  TBranch        *b_met_ey_JetEnUp_recoil;
  TBranch        *b_met_ex_JetEnDown_recoil;
  TBranch        *b_met_ey_JetEnDown_recoil;
  TBranch	       *b_met_ex_UnclusteredEnUp_recoil;
  TBranch	       *b_met_ey_UnclusteredEnUp_recoil;
  TBranch	       *b_met_ex_UnclusteredEnDown_recoil;
  TBranch	       *b_met_ey_UnclusteredEnDown_recoil;
  TBranch        *b_met_ex_recoil;
  TBranch        *b_met_ey_recoil;
  TBranch        *b_genmet;
  TBranch        *b_event_lumi;
  TBranch        *b_event_run;
  TBranch        *b_genmet_Ex;
  TBranch        *b_genmet_Ey;
  TBranch        *b_NuPx;
  TBranch        *b_NuPy;
  TBranch        *b_NuPz;
  TBranch        *b_NuPt;
  TBranch        *b_NuPhi;
  TBranch        *b_genmetphi;
  TBranch        *b_met_scaleUp;
  TBranch        *b_metphi_scaleUp;
  TBranch        *b_met_scaleDown;
  TBranch        *b_metphi_scaleDown;
  TBranch        *b_met_resoUp;
  TBranch        *b_met_resoDown;
  TBranch        *b_metphi_resoUp;
  TBranch        *b_metphi_resoDown;

  TBranch        *b_gen_weight;   //!
  TBranch	 *b_wScale0;
  TBranch	 *b_wScale1;
  TBranch	 *b_wScale2;
  TBranch	 *b_wScale3;
  TBranch	 *b_wScale4;
  TBranch	 *b_wScale5;
  TBranch	 *b_wScale6;
  TBranch	 *b_wScale7;
  TBranch	 *b_wScale8;
  TBranch	 *b_wPDFUp;
  TBranch	 *b_wPDFDown;
  TBranch        *b_pu_weight;   //!
  TBranch        *b_LSF_weight;   //!
  TBranch        *b_LSF_weight_mu;   //!
  TBranch        *b_LSF_weight_el;   //!
  TBranch        *b_LSF_weight_1;   //!
  TBranch        *b_LSF_weight_2;   //!
  TBranch        *b_TFR_weight;   //!
  TBranch        *b_top_weight;   //!
  TBranch        *b_all_weight;   //!
  TBranch        *b_trig_weight;   //!
  TBranch        *b_trig_weight_1;   //!
  TBranch        *b_trig_weight_2;   //!
  TBranch        *b_zptmassweight;   //!
  TBranch        *b_xsecs;   //!
  TBranch        *b_muon_index;   //!
  TBranch        *b_muon_index_1;   //!
  TBranch        *b_muon_index_2;   //!
  TBranch        *b_electron_index;   //!
  TBranch        *b_taus_index;   //!
  TBranch        *b_primvert_count;   //!
  TBranch        *b_primvert_x;   //!
  TBranch        *b_primvert_y;   //!
  TBranch        *b_primvert_z;   //!
  TBranch        *b_mu_count;   //!
  TBranch        *b_mu_px;   //!
  TBranch        *b_mu_py;   //!
  TBranch        *b_mu_pz;   //!
  TBranch        *b_mu_pt;   //!
  TBranch        *b_mu_eta;   //!
  TBranch        *b_mu_phi;   //!
  TBranch        *b_mu_charge;   //!
  TBranch        *b_mu_miniISO;   //!
  TBranch        *b_mu_dxy;   //!
  TBranch        *b_mu_dz;   //!
  TBranch        *b_mu_dxyerr;   //!
  TBranch        *b_mu_dzerr;   //!
  TBranch     *b_mu_neutralHadIso; 
  TBranch     *b_mu_photonIso; 
  TBranch     *b_mu_chargedHadIso; 
  TBranch     *b_mu_puIso; 
  TBranch     *b_mu_neutralIso;
  TBranch     *b_mu_absIsoMu; 
  TBranch     *b_mu_relIsoMu; 

  TBranch     *b_el_neutralHadIso; 
  TBranch     *b_el_photonIso; 
  TBranch     *b_el_chargedHadIso; 
  TBranch     *b_el_puIso; 
  TBranch     *b_el_neutralIso;
  TBranch     *b_el_absIsoEl; 
  TBranch     *b_el_relIsoEl; 
  TBranch     *b_el_isMVA; 
  TBranch        *b_mu_relIso;   //!
  TBranch        *b_jet_count;   //!
  TBranch        *b_npv;   //!
  TBranch        *b_genTauDecayMode1;   //!
  TBranch        *b_genTauDecayMode2;   //!
  TBranch        *b_genTauDecayModeAll;   //!
  TBranch        *b_npu;   //!
  TBranch        *b_jets_cleaned;   //!
  TBranch        *b_jet_e;   //!
  TBranch        *b_jet_px;   //!
  TBranch        *b_jet_py;   //!
  TBranch        *b_jet_pz;   //!
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_eta;   //!
  TBranch        *b_jet_phi;   //!
  TBranch        *b_jet_flavour;   //!
  TBranch        *b_jet_btag;   //!
  TBranch        *b_jet_isLoose;   //!
  TBranch        *b_el_count;   //!
  TBranch        *b_el_px;   //!
  TBranch        *b_el_py;   //!
  TBranch        *b_el_pz;   //!
  TBranch        *b_el_pt;   //!
  TBranch        *b_el_eta;   //!
  TBranch        *b_el_phi;   //!
  TBranch        *b_el_miniISO;   //!
  TBranch        *b_el_dxy;   //!
  TBranch        *b_el_dz;   //!
  TBranch        *b_el_dxyerr;   //!
  TBranch        *b_el_dzerr;   //!
  TBranch        *b_el_charge;   //!
  TBranch        *b_el_relIso;   //!
  TBranch        *b_ta_count;   //!
  TBranch        *b_ta_px;   //!
  TBranch        *b_ta_py;   //!
  TBranch        *b_ta_pz;   //!
  TBranch        *b_ta_mass;   //!
  TBranch        *b_ta_eta;   //!
  TBranch        *b_ta_phi;   //!
  TBranch        *b_ta_pt;   //!
  TBranch        *b_ta_dxy;   //!
  TBranch        *b_ta_dz;   //!
  TBranch        *b_ta_charge;   //!
  TBranch        *b_ta_IsoFlag;   //!
  TBranch        *b_ta_Iso;   //!
  TBranch        *b_ta_relIso;   //!
  TBranch        *b_datasetName;   //!
  TBranch        *b_CFCounter_;  //!
  TBranch        *b_regionName;   //!
  TBranch        *b_event_sign;   //!
  TBranch        *b_met_flag;   //!
  TBranch        *b_event_secondLeptonVeto;   //!
  TBranch        *b_eleMVA;   //!
  TBranch        *b_event_thirdLeptonVeto;   //!
  TBranch        *b_event_leptonDrTrigger;   //!
  TBranch        *b_genTauMatched;   //!
  TBranch        *b_genLeptonMatched;   //!
  TBranch        *b_genLeptonMatchedEl;   //!
  TBranch        *b_genLeptonMatchedMu;   //!
  TBranch        *b_genTauDecayedMuMatched;   //!
  TBranch        *b_genTauDecayedElMatched;   //!
  TBranch        *b_genLeptonPromptElMatched;   //!
  TBranch        *b_genLeptonPromptMuMatched;   //!
  TBranch        *b_genLeptonMatchedPrompEl;   //!
  TBranch        *b_genLeptonMatchedPrompMu;   //!
  TBranch        *b_genTauMatchedToTauDecay;   //!
  TBranch        *b_matchedTauToPromptEl;   //!
  TBranch        *b_matchedTauToPromptMu;   //!
  TBranch        *b_matchedTauToPromptTau;   //!
  TBranch        *b_matchedTauToTauDecEl;   //!
  TBranch        *b_matchedTauToTauDecMu;   //!
  TBranch        *b_matchedTauToTauDecTau;   //!
  TBranch        *b_matchedTauToElHadronDec;   //!
  TBranch        *b_matchedTauToMuHadronDec;   //!
  TBranch        *b_matchedTauToTauHadronDec;   //!
  TBranch        *b_matchedTauToGluon;   //!
  TBranch        *b_matchedTauToHFQ;   //!
  TBranch        *b_matchedTauToLFQ;   //!
  TBranch        *b_genLeptonMatchedPromptEl;   //!
  TBranch        *b_genLeptonMatchedPromptMu;   //!
  TBranch        *b_genLeptonMatchedPromptTau;   //!
  TBranch        *b_genElMatchedHadrDecay;   //!
  TBranch        *b_genMuMatchedHadrDecay;   //!
  TBranch        *b_genTauMatchedHadrDecay;   //!
  TBranch        *b_genLeptonMatchedGluon;   //!
  TBranch        *b_genLeptonMatchedLFQ;   //!
  TBranch        *b_genLeptonMatchedHFQ;   //!
  TBranch        *b_genElMatchedToTauDecay;   //!
  TBranch        *b_genMuMatchedToTauDecay;   //!
  
  TBranch        *b_isDYTT;   //!
  TBranch        *b_isDYEE;   //!
  TBranch        *b_isDYMM;   //!
  TBranch        *b_isDYNuNu;   //!





  TBranch        *b_qcdweight;   //!
  TBranch        *b_qcdweightup;   //!
  TBranch        *b_qcdweightdown;   //!
  TBranch        *b_npartons;   //!
  TBranch        *b_nbtag;   //!
  TBranch        *b_njets;   //!

  analyzer(TTree *tree=0);
  virtual ~analyzer();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  //virtual void main(int argc, char * argv[]) ;
  // virtual void Loop(int argc, char * argv[]) ;
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  /*		void TopHist(Double_t EvWeight){

		histTopPt->Fill(0.,EvWeight);
		}
  */



////////////////////fake study
  void FillHists(int CutIndex, Double_t EvWeight, TLorentzVector  elV, TLorentzVector  muV, TLorentzVector tauV, vector<TLorentzVector>  &JetsV, TLorentzVector  &MetV, double Chimass, double mintermediate, string & Sel, int  mIndex, int eIndex, int  tIndex, int &category, int &categoryJet, bool isdata, bool dogenMET){


if (CutIndex==17){
   int nSR = GetSRIndex(muV, tauV, JetsV, MetV, Chimass, mIndex, tIndex, CutIndex,dogenMET);


//cout<<" cat "<<category<<" catJet "<<categoryJet<<endl;

    CutFlowUnWFakeRate[CutIndex][0]->Fill(category+1,EvWeight);
    CutFlowUnWFakeRateJet[CutIndex][0]->Fill(categoryJet+1,EvWeight);

if (nSR!=-1){
    CutFlowUnWFakeRate[CutIndex][nSR]->Fill(category+1,EvWeight);
    CutFlowUnWFakeRateJet[CutIndex][nSR]->Fill(categoryJet+1,EvWeight);
}
    
}


    hnJetFake[CutIndex][0]->Fill(njets,EvWeight);
    if (!isdata)  hnJetFake[CutIndex][category]->Fill(njets,EvWeight);
			
    Mt2as = -1.;
    met = MetV.Pt();
    metx=MetV.Px();
    mety=MetV.Py();
    metphi=MetV.Phi();

if (dogenMET) { 
	met = sqrt(genmet_Ex*genmet_Ex + genmet_Ey*genmet_Ey);
	metx = genmet_Ex;
	mety = genmet_Ey;
	metphi = genmetphi;

		}



//cout<<""<<endl;

    ///////////////////mutau

    if (Sel=="mutau" && mIndex >-1 && tIndex>-1 ){

      tauMass = tauV.M();

      TLorentzVector t1 ; t1 = MetV+muV;

      double v1[4] = {muV.Px(),muV.Py(),muV.Pz(),muonMass};
      double v2[4] = {tauV.Px(),tauV.Py(),tauV.Pz(),tauMass};
      double ecm = 13000.0;
      double mxlo = Chimass;
      double ptm[2] = {metx,mety};
      double vds[4] = {0,0,0,0};

      Mt2as = asymm_mt2_lester_bisect::get_mT2(muonMass, muV.Px(), muV.Py(),tauMass,tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);

      double mcta = sqrt( 2*muV.Pt()*tauV.Pt()*(1+cos(muV.Phi()-tauV.Phi())) );

      double pa[3] ={muonMass, muV.Px(),muV.Py()};
      double pb[3] ={tauMass, tauV.Px(),tauV.Py()};
      double pmiss[3] ={0., metx,mety};
				 
      hMt2lesterFake[CutIndex][0]->Fill(Mt2as,EvWeight);
      if (!isdata)	hMt2lesterFake[CutIndex][category]->Fill(Mt2as,EvWeight);

      float tauUnitX = tauV.Px()/tauV.Pt();
      float tauUnitY = tauV.Py()/tauV.Pt();

      float muonUnitX = muV.Px()/muV.Pt();
      float muonUnitY = muV.Py()/muV.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + muV.Px() + tauV.Px();
      float vectorY = mety + muV.Py() + tauV.Py();

      float vectorVisX = muV.Px() + tauV.Px();
      float vectorVisY = muV.Py() + tauV.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;


//   int nSRr = GetSRIndex(muV, tauV, JetsV, MetV, Chimass, mIndex, tIndex,CutIndex);
//if (CutIndex==17)
//    cout<<" inside founder  again ..............MT2 "<<Mt2as<<" met "<<met<<" DZ "<<dzeta<<"  m "<<mIndex<<" t "<<tIndex<<"  "<<JetsV.size()<<"  Index  "<<CutIndex<<endl;

      hDZetaFake[CutIndex][0]->Fill(dzeta,EvWeight);
      if (!isdata)		hDZetaFake[CutIndex][category]->Fill(dzeta,EvWeight);


    }


    if (mIndex !=-1 )
      {


	hMuptFake[CutIndex][0]->Fill(mu_pt[mIndex],EvWeight);
	hMuetaFake[CutIndex][0]->Fill(mu_eta[mIndex],EvWeight);
	hIsoMuFake[CutIndex][0]->Fill(mu_relIso[0],EvWeight);

	if (!isdata){	hMuptFake[CutIndex][category]->Fill(mu_pt[mIndex],EvWeight);
	  hMuetaFake[CutIndex][category]->Fill(mu_eta[mIndex],EvWeight);
	  hIsoMuFake[CutIndex][category]->Fill(mu_relIso[0],EvWeight);
	}

      }


    if (tIndex !=-1)
      {
	hTauptFake[CutIndex][0]->Fill(ta_pt[tIndex],EvWeight);
	hTauetaFake[CutIndex][0]->Fill(ta_eta[tIndex],EvWeight);
	hIsoTauFake[CutIndex][0]->Fill(ta_relIso[0],EvWeight);
        if (genTauDecayMode1>-1)  hTauDecayMode1[CutIndex]->Fill(genTauDecayMode1,EvWeight);
        if (genTauDecayMode2>-1)  hTauDecayMode2[CutIndex]->Fill(genTauDecayMode2,EvWeight);
	hTauDecayModeAll[CutIndex]->Fill(genTauDecayMode1,EvWeight);
        hTauDecayModeAll[CutIndex]->Fill(genTauDecayMode2,EvWeight);
	if (!isdata){
	  hTauptFake[CutIndex][category]->Fill(ta_pt[tIndex],EvWeight);
	  hTauetaFake[CutIndex][category]->Fill(ta_eta[tIndex],EvWeight);
	  hIsoTauFake[CutIndex][category]->Fill(ta_relIso[0],EvWeight);
	}
      }




    /////////////global variables

    if (met>0. ) {

      hMETFake[CutIndex][0]->Fill(met,EvWeight);
      if (!isdata)		hMETFake[CutIndex][category]->Fill(met,EvWeight);
    }


      hGenMETFB[CutIndex]->Fill(genmet,EvWeight);



  }






  void FillHists(int CutIndex, Double_t EvWeight, TLorentzVector  elV, TLorentzVector  muV, TLorentzVector tauV, vector<TLorentzVector>  &JetsV, TLorentzVector  &MetV, double Chimass, double mintermediate, string & Sel, int  mIndex, int eIndex, int  tIndex, bool dogenMET){

    met = MetV.Pt();
    metx=MetV.Px();
    mety=MetV.Py();
    metphi = MetV.Phi();

//cout<<" before "<<met<<"  "<<metx<<"  "<<mety<<"  "<<metphi<<endl;

if (dogenMET) { met = sqrt(genmet_Ex*genmet_Ex + genmet_Ey*genmet_Ey);
	metx = genmet_Ex;
	mety = genmet_Ey;
	metphi = genmetphi;
	}
//cout<<" after "<<met<<"  "<<metx<<"  "<<mety<<"  "<<metphi<<endl;

    AllJets_Lepton_noMet.clear();
    for (unsigned int i = 0; i <   JetsV.size(); i++) AllJets_Lepton_noMet.push_back (JetsV.at(i));


    if (Sel=="mutau" )    AllJets_Lepton_noMet.push_back (muV);
    if (Sel=="WJetsmu" )    AllJets_Lepton_noMet.push_back (muV);
    if (Sel=="eltau" )    AllJets_Lepton_noMet.push_back (elV);
    if (Sel=="muel" )   { AllJets_Lepton_noMet.push_back (muV);AllJets_Lepton_noMet.push_back (elV);};

    //cout<<"  "<<muon_index<<"  "<<electron_index<<"  "<<taus_index<<"  "<<jet_count<<"  "<<met_ex<<"  "<<met_ey<<endl;

    //cout<<" Init  "<<mIndex<<"  "<<eIndex<<"  "<<tIndex<<endl;

    sumpT=0;

    double sumMuonpT=0;
    double sumElpT=0;
    double sumTaupT=0;  
    hnJet[CutIndex]->Fill(njets,EvWeight);
    hnpartons[CutIndex]->Fill(npartons,EvWeight);
    hnMu[CutIndex]->Fill(mu_count,EvWeight);
    hnTau[CutIndex]->Fill(ta_count,EvWeight);
    hnEl[CutIndex]->Fill(el_count,EvWeight);

    hWeights[CutIndex]->Fill(EvWeight);
    hnpv[CutIndex]->Fill(npv,EvWeight);
    hnpu[CutIndex]->Fill(npu,EvWeight);
    hEventSign[CutIndex]->Fill(event_sign,EvWeight);
			
    Mt2as = 1.;
    dzeta=-9999;

    ///////////////////mutau
    ////////////////////

    if (Sel=="mutau" && mIndex >-1 && tIndex>-1 ){


      tauMass = tauV.M();

      double deta=deltaEta(muV.Px(),  muV.Py(),muV.Pz(), tauV.Px(),  tauV.Py(),tauV.Pz() );
      hdEtaDil[CutIndex]->Fill(deta,EvWeight);

      TLorentzVector t1 ; t1 = MetV+muV;
      double detaM= muV.Eta()-t1.Eta();
      hdEtaMuMET[CutIndex]->Fill(detaM,EvWeight);

      double dPhiD=dPhiFrom2P( muV.Px(),  muV.Py(),tauV.Px(),  tauV.Py());
      hdPhiDil[CutIndex]->Fill(dPhiD,EvWeight);
				
      //   Mt2::Basic_Mt2_332_Calculator mt2Calculator;
      //   Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(muV.Px(), muV.Py()), muonMass);
      //   Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(tauV.Px(), tauV.Py()), tauMass);
      //   Mt2::TwoVector                pT_Miss(metx, mety);    
      //double mt2  = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, Chimass);

      double v1[4] = {muV.Px(),muV.Py(),muV.Pz(),muonMass};
      double v2[4] = {tauV.Px(),tauV.Py(),tauV.Pz(),tauMass};
      double ecm = 13000.0;
      double mxlo = Chimass;
      double ptm[2] = {metx,mety};
      double vds[4] = {0,0,0,0};
      // TMctLib t;

      Mt2as = asymm_mt2_lester_bisect::get_mT2(muonMass, muV.Px(), muV.Py(),tauMass,tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);
      //double Mt2as = 0; Mt2as = asymm_mt2_lester_bisect::get_mT2(muV.M(), muV.Px(), muV.Py(),tauV.M(),tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);
      //double Mt2ass = 0; Mt2ass = asymm_mt2_lester_bisect::get_mT2(muV.M(), muV.Px(), muV.Py(),tauMass,tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);
      double mcta = sqrt( 2*muV.Pt()*tauV.Pt()*(1+cos(muV.Phi()-tauV.Phi())) );

      //cout<<"  old "<<muonMass<<"  "<<muV.M()<<"  "<<tauMass<<"  "<<tauV.M()<<"  "<<Mt2as<<"  "<<endl;

				
      //   cout<<"MCTa  = "<<mcta<<" MCT = "<<mct(v1,v2)
      //   <<", MCTcorr = "<<mctcorr(v1,v2,vds,ptm,ecm,mxlo)
      //   <<", MT2 = "<<mt2(v1,v2,vds,ptm,ecm,mxlo)<<"  Mt2 lester "<<Mt2as
      //   <<", MCy = "<<mcy(v1,v2,vds,ptm)
      //  <<", MCx = "<<mcx(v1,v2,vds,ptm)
      //  <<endl;
      double pa[3] ={muonMass, muV.Px(),muV.Py()};
      double pb[3] ={tauMass, tauV.Px(),tauV.Py()};
      double pmiss[3] ={0., metx,mety};
      //mt2_bisect::mt2 mt2_event;
      //mt2_event.set_momenta( pa, pb, pmiss );
      //mt2_event.set_mn(Chimass);
      //double mt2_value = 0;//= mt2_event.get_mt2();
				 
      float MT2v =0 ; MT2v = mt2(v1,v2,vds,ptm,ecm,mxlo);

      //cout<<" new mt2 "<<mt2_value<<"  "<<MT2v<<endl;

      hMt2mutau[CutIndex]->Fill(MT2v,EvWeight);
      //hMt2[CutIndex]->Fill(mt2_value,EvWeight);
      hMt2lestermutau[CutIndex]->Fill(Mt2as,EvWeight);
      hMt2lestermutauFB[CutIndex]->Fill(Mt2as,EvWeight);
      hMCTmutau[CutIndex]->Fill(mct(v1,v2),EvWeight);
      hMCTxmutau[CutIndex]->Fill(mcx(v1,v2,vds,ptm),EvWeight);
      hMCTymutau[CutIndex]->Fill(mcy(v1,v2,vds,ptm),EvWeight);
      hMCTbmutau[CutIndex]->Fill(mcta,EvWeight);

      double MCTcorr = mctcorr(v1,v2,vds,ptm,ecm,mxlo);
      hMCTcor[CutIndex]->Fill(MCTcorr,EvWeight);

      //cout<<" asym "<<Mt2as<<"  "<<mt2<<endl;

      double Centr = Centrality(AllJets_Lepton_noMet);
      hCentrality[CutIndex]->Fill(Centr,EvWeight);
      //cout<<"Centr "<<Centr<<"  "<< EvWeight <<endl;
      //if(CutIndex>18){
      double mTB =0;//= Lester::mTBound(muV.M(), muV.Px(), muV.Py(), muV.Pz(), tauMass,tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);

      double mTtrue = Lester::mTTrue(muV.M(), muV.Px(), muV.Py(), muV.Pz(), tauMass,tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);


      hMTBoundmutau[CutIndex]->Fill(mTB,EvWeight);
      hMTTruemutau[CutIndex]->Fill(mTtrue,EvWeight);
      //	}
      //double mTB = Lester::mTTrue(muonMass, muV.Px(), muV.Py(), muV.Pz(), tauMass,tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);
      //			cout<<"  "<<mTB<<endl;
      //	     	double Dr= deltaR(MuV.at(mIndex).Eta(), MuV.at(mIndex).Phi(),tauV.at(mIndex).Eta(),tauV.at(mIndex).Phi());
      //double Dr= deltaR(muV.Eta(), muV.Phi(),tauV.Eta(),tauV.Phi());
      double Dr = muV.DeltaR(tauV);
      //
      TLorentzVector DiL = muV  + tauV;
      double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), metx,  mety );
      float MTDil = TMath::Sqrt(2*DiL.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTmutau[CutIndex]->Fill(MTDil,EvWeight);
      hInvMassMuTau[CutIndex]->Fill(DiL.M(),EvWeight);
      hdR_mutau[CutIndex]->Fill(Dr,EvWeight);




      dPhi=dPhiFrom2P( mu_px[mIndex], mu_py[mIndex], metx,  mety );
      //dPhi=dPhiFrom2P( muV.Px(), muV.Py(), metx,  mety );
      hdPhiMuMET[CutIndex]->Fill(dPhi,EvWeight);

      float MT = TMath::Sqrt(2*mu_pt[mIndex]*met*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hMTmu[CutIndex]->Fill(MT,EvWeight);
      hMTmuFineBin[CutIndex]->Fill(MT,EvWeight);
      dPhi=dPhiFrom2P( tauV.Px(), tauV.Py(), metx,  mety );
      float MTt = TMath::Sqrt(2*tauV.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTsum[CutIndex]->Fill(MT+MTt,EvWeight);
      hMTdif[CutIndex]->Fill(fabs(MT-MTt),EvWeight);
      hMTmax[CutIndex]->Fill(max(MT,MTt),EvWeight);
      hMTmin[CutIndex]->Fill(min(MT,MTt),EvWeight);
      float mttotal = sqrt(MT*MT + MTt*MTt);
      hMTtot[CutIndex]->Fill(mttotal,EvWeight);


      float tauUnitX = tauV.Px()/tauV.Pt();
      float tauUnitY = tauV.Py()/tauV.Pt();

      float muonUnitX = muV.Px()/muV.Pt();
      float muonUnitY = muV.Py()/muV.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + muV.Px() + tauV.Px();
      float vectorY = mety + muV.Py() + tauV.Py();

      float vectorVisX = muV.Px() + tauV.Px();
      float vectorVisY = muV.Py() + tauV.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;

			
if (CutIndex==17 &&  (int)event_lumi==101925 && (int) event_run ==1) 
	cout<< " More Info  vectorVisX "<<  vectorVisX<<" zetaX "<<zetaX<<" vectorVisY "<<vectorVisY<<" zetaY "<<zetaY<<endl;

      hDZeta[CutIndex]->Fill(dzeta,EvWeight);
      hDZetaFB[CutIndex]->Fill(dzeta,EvWeight);



      /*
	hmet_DZeta_MT2lester[CutIndex]->Fill(met,dzeta,Mt2as,EvWeight);
	hmet_MT[CutIndex]->Fill(met,MT,EvWeight);
	hmet_MTsum[CutIndex]->Fill(met,MT+MTt,EvWeight);
	hmet_MTtot[CutIndex]->Fill(met,mttotal,EvWeight);
	hmet_MCTb[CutIndex]->Fill(met,mcta,EvWeight);
	hmet_MT2lester[CutIndex]->Fill(met,Mt2as,EvWeight);
	hmet_TauPt[CutIndex]->Fill(met,tauV.Pt(),EvWeight);
	hmet_DZeta[CutIndex]->Fill(met,dzeta,EvWeight);
	hmet_dR[CutIndex]->Fill(met,Dr,EvWeight);
	hmet_MTDil[CutIndex]->Fill(met,MTDil,EvWeight);

	hMT_MTsum[CutIndex]->Fill(MT,MT+MTt,EvWeight);
	hMT_MTtot[CutIndex]->Fill(MT,mttotal,EvWeight);
	hMT_MCTb[CutIndex]->Fill(MT,mcta,EvWeight);
	hMT_MT2lester[CutIndex]->Fill(MT,Mt2as,EvWeight);
	hMT_TauPt[CutIndex]->Fill(MT,tauV.Pt(),EvWeight);
	hMT_DZeta[CutIndex]->Fill(MT,dzeta,EvWeight);
	hMT_dR[CutIndex]->Fill(MT,Dr,EvWeight);
	hMT_MTDil[CutIndex]->Fill(MT,MTDil,EvWeight);

	hMTsum_MTtot[CutIndex]->Fill(MT+MTt,mttotal,EvWeight);
	hMTsum_MCTb[CutIndex]->Fill(MT+MTt,mcta,EvWeight);
	hMTsum_MT2lester[CutIndex]->Fill(MT+MTt,Mt2as,EvWeight);
	hMTsum_TauPt[CutIndex]->Fill(MT+MTt,tauV.Pt(),EvWeight);
	hMTsum_DZeta[CutIndex]->Fill(MT+MTt,dzeta,EvWeight);
	hMTsum_dR[CutIndex]->Fill(MT+MTt,Dr,EvWeight);
	hMTsum_MTDil[CutIndex]->Fill(MT+MTt,MTDil,EvWeight);

	hMTtot_MCTb[CutIndex]->Fill(mttotal,mcta,EvWeight);
	hMTtot_MT2lester[CutIndex]->Fill(mttotal,Mt2as,EvWeight);
	hMTtot_TauPt[CutIndex]->Fill(mttotal,tauV.Pt(),EvWeight);
	hMTtot_DZeta[CutIndex]->Fill(mttotal,dzeta,EvWeight);
	hMTtot_dR[CutIndex]->Fill(mttotal,Dr,EvWeight);
	hMTtot_MTDil[CutIndex]->Fill(mttotal,MTDil,EvWeight);

	hMCTb_MT2lester[CutIndex]->Fill(mcta,Mt2as,EvWeight);
	hMCTb_TauPt[CutIndex]->Fill(mcta,tauV.Pt(),EvWeight);
	hMCTb_DZeta[CutIndex]->Fill(mcta,dzeta,EvWeight);
	hMCTb_dR[CutIndex]->Fill(mcta,Dr,EvWeight);
	hMCTb_MTDil[CutIndex]->Fill(mcta,MTDil,EvWeight);

	hMT2lester_TauPt[CutIndex]->Fill(Mt2as,tauV.Pt(),EvWeight);
	hMT2lester_DZeta[CutIndex]->Fill(Mt2as,dzeta,EvWeight);
	hMT2lester_dR[CutIndex]->Fill(Mt2as,Dr,EvWeight);
	hMT2lester_MTDil[CutIndex]->Fill(Mt2as,MTDil,EvWeight);

	hTauPt_DZeta[CutIndex]->Fill(tauV.Pt(),dzeta,EvWeight);
	hTauPt_dR[CutIndex]->Fill(tauV.Pt(),Dr,EvWeight);
	hTauPt_MTDil[CutIndex]->Fill(tauV.Pt(),MTDil,EvWeight);

	hdR_MTDil[CutIndex]->Fill(Dr,MTDil,EvWeight);
      */

      if (JetsV.size()>0){
	double dPhiJT=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), tauV.Px(),  tauV.Py() );
	hdPhiJ0Tau[CutIndex]->Fill(dPhiJT,EvWeight);
	double dPhiJL=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), muV.Px(),  muV.Py() );
	hdPhiJ0Mu[CutIndex]->Fill(dPhiJL,EvWeight);
				
	double deta = muV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Mu[CutIndex]->Fill(deta,EvWeight);

	double deta2 = tauV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Tau[CutIndex]->Fill(deta2,EvWeight);

	double Drtj = tauV.DeltaR(JetsV.at(0));
	hdR_taujet[CutIndex]->Fill(Drtj,EvWeight);
      }

    }
    if (Sel=="eltau" && eIndex > -1 && tIndex>-1){

      tauMass = tauV.M();

      double detael=deltaEta(elV.Px(),  elV.Py(),elV.Pz(), tauV.Px(),  tauV.Py(),tauV.Pz() );
      hdEtaDil[CutIndex]->Fill(detael,EvWeight);

      TLorentzVector t1 ; t1 = MetV+elV;
      double detaM= elV.Eta()-t1.Eta();
      hdEtaElMET[CutIndex]->Fill(detaM,EvWeight);

      double dPhiD=dPhiFrom2P( elV.Px(),  elV.Py(),tauV.Px(),  tauV.Py());
      hdPhiDil[CutIndex]->Fill(dPhiD,EvWeight);

      //double Dr= deltaR(elV.Eta(), elV.Phi(),tauV.Eta(),tauV.Phi());
      double v1[4] = {elV.Px(),elV.Py(),elV.Pz(),electronMass};
      double v2[4] = {tauV.Px(),tauV.Py(),tauV.Pz(),tauMass};
      double ecm = 13000.0;
      double mxlo = Chimass;
      double ptm[2] = {metx,mety};
      double vds[4] = {0,0,0,0};
			
      Mt2as = asymm_mt2_lester_bisect::get_mT2(electronMass, elV.Px(), elV.Py(),tauMass,tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);

      //   cout<<" Will call with Mt2 "<<Mt2<<"   "<<Chimass<<endl;
      hMt2eltau[CutIndex]->Fill(Mt2as,EvWeight);
      double mcta = sqrt( 2*elV.Pt()*tauV.Pt()*(1+cos(elV.Phi()-tauV.Phi())) );
      float MT2v =0 ; MT2v = mt2(v1,v2,vds,ptm,ecm,mxlo);
      hMt2eltau[CutIndex]->Fill(MT2v,EvWeight);
      hMt2lestereltau[CutIndex]->Fill(Mt2as,EvWeight);
      hMt2lestereltauFB[CutIndex]->Fill(Mt2as,EvWeight);
			
      hMCTeltau[CutIndex]->Fill(mct(v1,v2),EvWeight);
      hMCTxeltau[CutIndex]->Fill(mcx(v1,v2,vds,ptm),EvWeight);
      hMCTyeltau[CutIndex]->Fill(mcy(v1,v2,vds,ptm),EvWeight);
      hMCTbeltau[CutIndex]->Fill(mcta,EvWeight);
      double MCTcorr = mctcorr(v1,v2,vds,ptm,ecm,mxlo);
      hMCTcor[CutIndex]->Fill(MCTcorr,EvWeight);


      double Centr = Centrality(AllJets_Lepton_noMet);
      hCentrality[CutIndex]->Fill(Centr,EvWeight);
      //cout<<"Centr "<<Centr<<"  "<< EvWeight <<endl;
      //double mTB = Lester::mTBound(elV.E(), elV.Px(), elV.Py(), elV.Pz(), tauV.E(),tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);
      //hTBoundeltau[CutIndex]->Fill(mTB);
      //double mTB = Lester::mTTrue(elV.E(), muV.Px(), muV.Py(), muV.Pz(), tauMass,tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);
      //hTBoundeltau[CutIndex]->Fill(mTB,EvWeight);
      //cout<<"  "<<mTB<<endl;
      //	double mTtrue = Lester::mTTrue(elV.M(), elV.Px(), elV.Py(), elV.Pz(), tauMass,tauV.Px(),tauV.Py(),tauV.Pz(),metx,mety, mintermediate);
      double Dr = elV.DeltaR(tauV);


      TLorentzVector DiL = elV + tauV;
      double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), metx,  mety );
      float MTDil = TMath::Sqrt(2*DiL.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTeltau[CutIndex]->Fill(MTDil,EvWeight);
      hInvMassElTau[CutIndex]->Fill(DiL.M(),EvWeight);
      hdR_eltau[CutIndex]->Fill(Dr,EvWeight);



      dPhi=dPhiFrom2P( elV.Px(), elV.Py(), metx,  mety );
      hdPhiElMET[CutIndex]->Fill(dPhi,EvWeight);

      float MT = TMath::Sqrt(2*elV.Pt()*met*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hMTel[CutIndex]->Fill(MT,EvWeight);
      dPhi=dPhiFrom2P( tauV.Px(), tauV.Py(), metx,  mety );
      float MTt = TMath::Sqrt(2*tauV.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTsum[CutIndex]->Fill(MT+MTt,EvWeight);
      hMTdif[CutIndex]->Fill(fabs(MT-MTt),EvWeight);
      hMTmax[CutIndex]->Fill(max(MT,MTt),EvWeight);
      hMTmin[CutIndex]->Fill(min(MT,MTt),EvWeight);
      float mttotal = sqrt(MT*MT + MTt*MTt);
      hMTtot[CutIndex]->Fill(mttotal,EvWeight);


      float tauUnitX = tauV.Px()/tauV.Pt();
      float tauUnitY = tauV.Py()/tauV.Pt();

      float elUnitX = elV.Px()/elV.Pt();
      float elUnitY = elV.Py()/elV.Pt();

      float zetaX = tauUnitX + elUnitX;
      float zetaY = tauUnitY + elUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + elV.Px() + tauV.Px();
      float vectorY = mety + elV.Py() + tauV.Py();

      float vectorVisX = elV.Px() + tauV.Px();
      float vectorVisY = elV.Py() + tauV.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;


      hDZeta[CutIndex]->Fill(dzeta,EvWeight);
      hDZetaFB[CutIndex]->Fill(dzeta,EvWeight);
      /*
	hmet_DZeta_MT2lester[CutIndex]->Fill(met,dzeta,Mt2as,EvWeight);

	hmet_MT[CutIndex]->Fill(met,MT,EvWeight);
	hmet_MTsum[CutIndex]->Fill(met,MT+MTt,EvWeight);
	hmet_MTtot[CutIndex]->Fill(met,mttotal,EvWeight);
	hmet_MCTb[CutIndex]->Fill(met,mcta,EvWeight);
	hmet_MT2lester[CutIndex]->Fill(met,Mt2as,EvWeight);
	hmet_TauPt[CutIndex]->Fill(met,tauV.Pt(),EvWeight);
	hmet_DZeta[CutIndex]->Fill(met,dzeta,EvWeight);
	hmet_dR[CutIndex]->Fill(met,Dr,EvWeight);
	hmet_MTDil[CutIndex]->Fill(met,MTDil,EvWeight);

	hMT_MTsum[CutIndex]->Fill(MT,MT+MTt,EvWeight);
	hMT_MTtot[CutIndex]->Fill(MT,mttotal,EvWeight);
	hMT_MCTb[CutIndex]->Fill(MT,mcta,EvWeight);
	hMT_MT2lester[CutIndex]->Fill(MT,Mt2as,EvWeight);
	hMT_TauPt[CutIndex]->Fill(MT,tauV.Pt(),EvWeight);
	hMT_DZeta[CutIndex]->Fill(MT,dzeta,EvWeight);
	hMT_dR[CutIndex]->Fill(MT,Dr,EvWeight);
	hMT_MTDil[CutIndex]->Fill(MT,MTDil,EvWeight);

	hMTsum_MTtot[CutIndex]->Fill(MT+MTt,mttotal,EvWeight);
	hMTsum_MCTb[CutIndex]->Fill(MT+MTt,mcta,EvWeight);
	hMTsum_MT2lester[CutIndex]->Fill(MT+MTt,Mt2as,EvWeight);
	hMTsum_TauPt[CutIndex]->Fill(MT+MTt,tauV.Pt(),EvWeight);
	hMTsum_DZeta[CutIndex]->Fill(MT+MTt,dzeta,EvWeight);
	hMTsum_dR[CutIndex]->Fill(MT+MTt,Dr,EvWeight);
	hMTsum_MTDil[CutIndex]->Fill(MT+MTt,MTDil,EvWeight);

	hMTtot_MCTb[CutIndex]->Fill(mttotal,mcta,EvWeight);
	hMTtot_MT2lester[CutIndex]->Fill(mttotal,Mt2as,EvWeight);
	hMTtot_TauPt[CutIndex]->Fill(mttotal,tauV.Pt(),EvWeight);
	hMTtot_DZeta[CutIndex]->Fill(mttotal,dzeta,EvWeight);
	hMTtot_dR[CutIndex]->Fill(mttotal,Dr,EvWeight);
	hMTtot_MTDil[CutIndex]->Fill(mttotal,MTDil,EvWeight);

	hMCTb_MT2lester[CutIndex]->Fill(mcta,Mt2as,EvWeight);
	hMCTb_TauPt[CutIndex]->Fill(mcta,tauV.Pt(),EvWeight);
	hMCTb_DZeta[CutIndex]->Fill(mcta,dzeta,EvWeight);
	hMCTb_dR[CutIndex]->Fill(mcta,Dr,EvWeight);
	hMCTb_MTDil[CutIndex]->Fill(mcta,MTDil,EvWeight);

	hMT2lester_TauPt[CutIndex]->Fill(Mt2as,tauV.Pt(),EvWeight);
	hMT2lester_DZeta[CutIndex]->Fill(Mt2as,dzeta,EvWeight);
	hMT2lester_dR[CutIndex]->Fill(Mt2as,Dr,EvWeight);
	hMT2lester_MTDil[CutIndex]->Fill(Mt2as,MTDil,EvWeight);

	hTauPt_DZeta[CutIndex]->Fill(tauV.Pt(),dzeta,EvWeight);
	hTauPt_dR[CutIndex]->Fill(tauV.Pt(),Dr,EvWeight);
	hTauPt_MTDil[CutIndex]->Fill(tauV.Pt(),MTDil,EvWeight);

	hdR_MTDil[CutIndex]->Fill(Dr,MTDil,EvWeight);
      */



      if (JetsV.size()>0){
	double dPhiJT=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), tauV.Px(),  tauV.Py() );
	hdPhiJ0Tau[CutIndex]->Fill(dPhiJT,EvWeight);
	double dPhiJL=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), elV.Px(),  elV.Py() );
	hdPhiJ0El[CutIndex]->Fill(dPhiJL,EvWeight);

	double deta = elV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0El[CutIndex]->Fill(deta,EvWeight);

	double deta2 = tauV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Tau[CutIndex]->Fill(deta2,EvWeight);

	double Drtj = tauV.DeltaR(JetsV.at(0));
	hdR_taujet[CutIndex]->Fill(Drtj,EvWeight);
      }

    }////////////end of eltau channe


    if (Sel=="muel" &&  eIndex >-1 &&  mIndex>-1){

      double deta=deltaEta(muV.Px(),  muV.Py(),muV.Pz(), elV.Px(),  elV.Py(),elV.Pz() );
      hdEtaDil[CutIndex]->Fill(deta,EvWeight);

      // of course the following variables do not make any sense....
      TLorentzVector t1 ; t1 = MetV+muV;
      double detaM= muV.Eta()-t1.Eta();
      //			hdEtaMuMET[CutIndex]->Fill(detaM,EvWeight);

      TLorentzVector t2 ; t2 = MetV+elV;
      double detaM2= muV.Eta()-t2.Eta();
      //			hdEtaElMET[CutIndex]->Fill(detaM2,EvWeight);


      double dPhiD=dPhiFrom2P( muV.Px(),  muV.Py(), elV.Px(),  elV.Py());
      hdPhiDil[CutIndex]->Fill(dPhiD,EvWeight);

      TLorentzVector DiL = muV  + elV;

      double detaM3= muV.Eta()-DiL.Eta();
      //			hdEtaDilMET[CutIndex]->Fill(detaM2,EvWeight);

      double v1[4] = {muV.Px(),muV.Py(),muV.Pz(),muonMass};
      double v2[4] = {elV.Px(),elV.Py(),elV.Pz(),electronMass};
      double ecm = 13000.0;
      double mxlo = Chimass;
      double ptm[2] = {metx,mety};
      double vds[4] = {0,0,0,0};
      // TMctLib t;

      Mt2as = asymm_mt2_lester_bisect::get_mT2(muonMass, muV.Px(), muV.Py(),electronMass,elV.Px(),elV.Py(),metx,mety,Chimass,Chimass,0);
      double mcta = sqrt( 2*muV.Pt()*elV.Pt()*(1+cos(muV.Phi()-elV.Phi())) );

      float MT2v =0 ; MT2v = mt2(v1,v2,vds,ptm,ecm,mxlo);

      hMt2muel[CutIndex]->Fill(MT2v,EvWeight);
      hMt2lestermuel[CutIndex]->Fill(Mt2as,EvWeight);
      hMt2lestermuelFB[CutIndex]->Fill(Mt2as,EvWeight);
      hMCTmuel[CutIndex]->Fill(mct(v1,v2),EvWeight);
      hMCTxmuel[CutIndex]->Fill(mcx(v1,v2,vds,ptm),EvWeight);
      hMCTymuel[CutIndex]->Fill(mcy(v1,v2,vds,ptm),EvWeight);
      hMCTbmuel[CutIndex]->Fill(mcta,EvWeight);
      double MCTcorr = mctcorr(v1,v2,vds,ptm,ecm,mxlo);
      hMCTcor[CutIndex]->Fill(MCTcorr,EvWeight);

      //cout<<" asym "<<Mt2as<<"  "<<mt2<<endl;

      double Centr = Centrality(AllJets_Lepton_noMet);
      hCentrality[CutIndex]->Fill(Centr,EvWeight);
      //cout<<"Centr "<<Centr<<"  "<< EvWeight <<endl;
      //double mTB = Lester::mTBound(muV.E(), muV.Px(), muV.Py(), muV.Pz(), elV.E(),elV.Px(),elV.Py(),elV.Pz(),metx,mety, mintermediate);
      //hTBoundmuel[CutIndex]->Fill(mTB,EvWeight);
      //cout<<"  "<<mTB<<endl;

      double Dr = muV.DeltaR(elV);

      double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), metx,  mety );
      float MTDil = TMath::Sqrt(2*DiL.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTmuel[CutIndex]->Fill(MTDil,EvWeight);
      hInvMassMuEl[CutIndex]->Fill(DiL.M(),EvWeight);
      hdR_muel[CutIndex]->Fill(Dr,EvWeight);


      dPhi=dPhiFrom2P( mu_px[mIndex], mu_py[mIndex], metx,  mety );
      hdPhiMuMET[CutIndex]->Fill(dPhi,EvWeight);

      float MT = TMath::Sqrt(2*mu_pt[mIndex]*met*(1-TMath::Cos(dPhi)));
      hMTmu[CutIndex]->Fill(MT,EvWeight);
      hMTmuFineBin[CutIndex]->Fill(MT,EvWeight);
      dPhi=dPhiFrom2P( elV.Px(), elV.Py(), metx,  mety );
      float MTel = TMath::Sqrt(2*elV.Pt()*met*(1-TMath::Cos(dPhi)));
					
      dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), metx,  mety );
      float MTelmu = TMath::Sqrt(2*DiL.Pt()*met*(1-TMath::Cos(dPhi)));

      hMTsum[CutIndex]->Fill(MT+MTel+MTelmu,EvWeight);
      hMTdif[CutIndex]->Fill(fabs(MT-MTel-MTelmu),EvWeight);
      hMT[CutIndex]->Fill(MTelmu,EvWeight);

      float mttotal = sqrt(MT*MT + MTel*MTel+ MTelmu*MTelmu);
      hMTtot[CutIndex]->Fill(mttotal,EvWeight);
				
      hMTel[CutIndex]->Fill(MTel,EvWeight);



      float tauUnitX = elV.Px()/elV.Pt();
      float tauUnitY = elV.Py()/elV.Pt();

      float muonUnitX = muV.Px()/muV.Pt();
      float muonUnitY = muV.Py()/muV.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + muV.Px() + elV.Px();
      float vectorY = mety + muV.Py() + elV.Py();

      float vectorVisX = muV.Px() + elV.Px();
      float vectorVisY = muV.Py() + elV.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;


      hDZeta[CutIndex]->Fill(dzeta,EvWeight);
      hDZetaFB[CutIndex]->Fill(dzeta,EvWeight);




      if (JetsV.size()>0){
	double dPhiJL=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), muV.Px(),  muV.Py() );
	hdPhiJ0Mu[CutIndex]->Fill(dPhiJL,EvWeight);
	dPhiJL=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), elV.Px(),  elV.Py() );
	hdPhiJ0El[CutIndex]->Fill(dPhiJL,EvWeight);

	double deta = muV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Mu[CutIndex]->Fill(deta,EvWeight);
	double deta1 = elV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0El[CutIndex]->Fill(deta1,EvWeight);
      }

      /*
	hmet_DZeta_MT2lester[CutIndex]->Fill(met,dzeta,Mt2as,EvWeight);

	hmet_MT[CutIndex]->Fill(met,MT,EvWeight);
	hmet_MTsum[CutIndex]->Fill(met,MT+MTel,EvWeight);
	hmet_MTtot[CutIndex]->Fill(met,mttotal,EvWeight);
	hmet_MCTb[CutIndex]->Fill(met,mcta,EvWeight);
	hmet_MT2lester[CutIndex]->Fill(met,Mt2as,EvWeight);
	hmet_DZeta[CutIndex]->Fill(met,dzeta,EvWeight);
	hmet_dR[CutIndex]->Fill(met,Dr,EvWeight);
	hmet_MTDil[CutIndex]->Fill(met,MTDil,EvWeight);

	hMT_MTsum[CutIndex]->Fill(MT,MT+MTel,EvWeight);
	hMT_MTtot[CutIndex]->Fill(MT,mttotal,EvWeight);
	hMT_MCTb[CutIndex]->Fill(MT,mcta,EvWeight);
	hMT_MT2lester[CutIndex]->Fill(MT,Mt2as,EvWeight);
	hMT_DZeta[CutIndex]->Fill(MT,dzeta,EvWeight);
	hMT_dR[CutIndex]->Fill(MT,Dr,EvWeight);
	hMT_MTDil[CutIndex]->Fill(MT,MTDil,EvWeight);

	hMTsum_MTtot[CutIndex]->Fill(MT+MTel,mttotal,EvWeight);
	hMTsum_MCTb[CutIndex]->Fill(MT+MTel,mcta,EvWeight);
	hMTsum_MT2lester[CutIndex]->Fill(MT+MTel,Mt2as,EvWeight);
	hMTsum_DZeta[CutIndex]->Fill(MT+MTel,dzeta,EvWeight);
	hMTsum_dR[CutIndex]->Fill(MT+MTel,Dr,EvWeight);
	hMTsum_MTDil[CutIndex]->Fill(MT+MTel,MTDil,EvWeight);

	hMTtot_MCTb[CutIndex]->Fill(mttotal,mcta,EvWeight);
	hMTtot_MT2lester[CutIndex]->Fill(mttotal,Mt2as,EvWeight);
	hMTtot_DZeta[CutIndex]->Fill(mttotal,dzeta,EvWeight);
	hMTtot_dR[CutIndex]->Fill(mttotal,Dr,EvWeight);
	hMTtot_MTDil[CutIndex]->Fill(mttotal,MTDil,EvWeight);

	hMCTb_MT2lester[CutIndex]->Fill(mcta,Mt2as,EvWeight);
	hMCTb_TauPt[CutIndex]->Fill(mcta,tauV.Pt(),EvWeight);
	hMCTb_DZeta[CutIndex]->Fill(mcta,dzeta,EvWeight);
	hMCTb_dR[CutIndex]->Fill(mcta,Dr,EvWeight);
	hMCTb_MTDil[CutIndex]->Fill(mcta,MTDil,EvWeight);

	hMT2lester_TauPt[CutIndex]->Fill(Mt2as,tauV.Pt(),EvWeight);
	hMT2lester_DZeta[CutIndex]->Fill(Mt2as,dzeta,EvWeight);
	hMT2lester_dR[CutIndex]->Fill(Mt2as,Dr,EvWeight);
	hMT2lester_MTDil[CutIndex]->Fill(Mt2as,MTDil,EvWeight);


	hdR_MTDil[CutIndex]->Fill(Dr,MTDil,EvWeight);
      */


    }/////////////////end of muel channel




    /////////////////////////////////////////////1D distributions
    /////////////////////////////////////////////
    ////////////////////////////////////////////////

      	int nB = -1;
	int histoIndex = 0;
      if (JetsV.size()==0) nB = 0;
      if (JetsV.size()==1) nB = 24;



    if (dzeta > -500 && met>0 && Mt2as>0 && JetsV.size()==0){
	
	if (met<40){
				
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+1,EvWeight);histoIndex=nB+1; }
	    if (dzeta>-100 && dzeta<0) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+2,EvWeight);histoIndex=nB+2;}
	    if (dzeta>0) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+3,EvWeight);histoIndex=nB+3;}
	  }
	  if (Mt2as>40) 	{
	    { hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+4,EvWeight);histoIndex=nB+4;}
	  }


	}//////met < 40 bin
	if (met>40 && met<80){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+5,EvWeight);histoIndex=nB+5; 
	   if (CutIndex==17) {
		   cout<<" FOUND THIS event here.............EvntWeight "<<EvWeight<<" run no "<<event_run<<" lumi "<<event_lumi<<" MT2 "<<Mt2as<<" DZeta "<<dzeta<<" GenMET "<<genmet<<" pfMET "<<met<<" met "<<endl;
		   cout<<" Pzetamiss  "<<pzetamiss<<" Pzetavisible "<<pzetavis<<"  "<<endl;
		   cout<<" SumPx(neu) "<<NuPx<<" genMet_ex "<<genmet_Ex<<" pfMet_Ex "<<met_ex_recoil<< endl; 
		   cout<<" SumPy(neu) "<<NuPy<<" genMet_ey "<<genmet_Ey<<" pfMet_Ey "<<met_ey_recoil<< endl;
		   cout<<" SumPt(neu) "<<NuPt<<" genMet "<<genmet<<" pfMet "<<met<<" met "<<met<<" muPt "<<muV.Pt()<<" TaupT " <<tauV.Pt()<<endl;
		   
		   }
	    
	    }
	    if (dzeta>-100 && dzeta<50) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+6,EvWeight);histoIndex=nB+6;}
	    if (dzeta>50) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+7,EvWeight);histoIndex=nB+7;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+8,EvWeight);histoIndex=nB+8;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+9,EvWeight);histoIndex=nB+9;}
	  }

	  if (Mt2as>80 ) 	{
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+10,EvWeight);histoIndex=nB+10;}
	  }


	}///met 40-80 bin

	if (met>80 && met<120){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+11,EvWeight);histoIndex=nB+11;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+12,EvWeight);histoIndex=nB+12;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+13,EvWeight);histoIndex=nB+13;}
	    if (dzeta>-150 ) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+14,EvWeight);histoIndex=nB+14;}
	  }

	  if (Mt2as>80) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+15,EvWeight);histoIndex=nB+15;}
	  
	  


	}///met 80-120 bin

	if (met>120 && met < 250){
							
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+16,EvWeight);histoIndex=nB+16;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+17,EvWeight);histoIndex=nB+17;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+18,EvWeight);histoIndex=nB+18;}
	    if (dzeta>-150 && dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+19,EvWeight);histoIndex=nB+19;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+20,EvWeight);histoIndex=nB+20;}
	  }

	  if (Mt2as>80 && Mt2as<100) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+21,EvWeight);histoIndex=nB+21;}
	  
	  if (Mt2as>100 && Mt2as<120) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+22,EvWeight);histoIndex=nB+22;}
	  
	  if (Mt2as>120) 	{
	    hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+23,EvWeight);histoIndex=nB+23; }


	}///met 120 bin & < 250
	if (met>250){
							
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+24,EvWeight);histoIndex=nB+24;}
	  

	}///met 250 bin /// Jets == 0
	
//	cout<<" MET here  "<<genmet<<"  "<<MetDif<<"  "<<MetDif/genmet<<endl;



if (histoIndex>0){
	float MetDif = genmet -met;
	hDeltaMET[CutIndex][histoIndex]->Fill(MetDif,EvWeight);
	hDeltaMETRel[CutIndex][histoIndex]->Fill(MetDif/genmet,EvWeight);
	hmet1D[CutIndex]->Fill(MetDif,EvWeight);
	}

    }/////////dzeta>-9999




    if (dzeta > -500 && met>0 && Mt2as>0 && JetsV.size()==1){
	
	if (met<40){
				
	  if (Mt2as<40) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+1,EvWeight);histoIndex=nB+1;}
	    if (dzeta>-150 && dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+2,EvWeight);histoIndex=nB+2; }
	    if (dzeta>-100 && dzeta<0) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+3,EvWeight);histoIndex=nB+3;}
	    if (dzeta>0) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+4,EvWeight);histoIndex=nB+4;}
	  }
	  if (Mt2as>40) 	{
	    { hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+5,EvWeight);histoIndex=nB+5;}
	  }


	}//////met < 40 bin
	if (met>40 && met<80){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+6,EvWeight);histoIndex=nB+6;}
	    if (dzeta>-100 && dzeta<50) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+7,EvWeight);histoIndex=nB+7;}
	    if (dzeta>50) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+8,EvWeight);histoIndex=nB+8;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+9,EvWeight);histoIndex=nB+9;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+10,EvWeight);histoIndex=nB+10;}
	  }

	  if (Mt2as>80 ) 	{
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+11,EvWeight);histoIndex=nB+11;}
	  }


	}///met 40-80 bin

	if (met>80 && met<120){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+12,EvWeight);histoIndex=nB+12;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+13,EvWeight);histoIndex=nB+13;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+14,EvWeight);histoIndex=nB+14;;}
	    if (dzeta>-150 ) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+15,EvWeight);histoIndex=nB+15;}
	  }

	  if (Mt2as>80 && Mt2as<120) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+16,EvWeight);histoIndex=nB+16;}
	  
	  if (Mt2as>120) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+17,EvWeight);histoIndex=nB+17;}
	  


	}///met 80-120 bin

	if (met>120 && met < 250){
							
	  if (Mt2as<40) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+18,EvWeight);histoIndex=nB+18;}
	    if (dzeta>-150 && dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+19,EvWeight);histoIndex=nB+19;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+20,EvWeight);histoIndex=nB+20;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+21,EvWeight);histoIndex=nB+21;}
	    if (dzeta>-150 && dzeta<-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+22,EvWeight);histoIndex=nB+22;}
	    if (dzeta>-100) {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+23,EvWeight);histoIndex=nB+23;}
	  }

	  if (Mt2as>80 && Mt2as<100) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+24,EvWeight);histoIndex=nB+24;}
	  
	  if (Mt2as>100 && Mt2as<120) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+25,EvWeight);histoIndex=nB+25;}
	  
	  if (Mt2as>120) 	{
	    hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+26,EvWeight);histoIndex=nB+26; }


	}///met 120 bin & < 250
	if (met>250){
							
	  if (Mt2as>80 && Mt2as<100) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+27,EvWeight);histoIndex=nB+27;}
	  
	  if (Mt2as>100 && Mt2as<120) 	
	    {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+28,EvWeight);histoIndex=nB+28;}
	  
	  if (Mt2as>120) 	
	  {hmet_MT2lester_DZeta_01J1D[CutIndex]->Fill(nB+29,EvWeight);histoIndex=nB+29; }

	}///met 250 bin
	
//	cout<<" MET here  "<<genmet<<"  "<<MetDif<<"  "<<MetDif/genmet<<endl;



if (histoIndex>0){
	float MetDif = genmet -met;
	hDeltaMET[CutIndex][histoIndex]->Fill(MetDif,EvWeight);
	hDeltaMETRel[CutIndex][histoIndex]->Fill(MetDif/genmet,EvWeight);
	hmet1D[CutIndex]->Fill(MetDif,EvWeight);
}

    }/////////dzeta>-9999

      ///////////////////////
      ///////////////////////
      ///////////////////////

    if ( (Sel=="WJetsmu" || Sel=="Channel") && mIndex >-1 ){


      TLorentzVector t1 ; t1 = MetV+muV;
      double detaM= muV.Eta()-t1.Eta();
      hdEtaMuMET[CutIndex]->Fill(detaM,EvWeight);
			

      double Centr = Centrality(AllJets_Lepton_noMet);
      hCentrality[CutIndex]->Fill(Centr,EvWeight);


      double dPhi=dPhiFrom2P( mu_px[mIndex], mu_py[mIndex], metx,  mety );
      //dPhi=dPhiFrom2P( muV.Px(), muV.Py(), metx,  mety );
      hdPhiMuMET[CutIndex]->Fill(dPhi,EvWeight);

      float MT = TMath::Sqrt(2*mu_pt[mIndex]*met*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hMTmu[CutIndex]->Fill(MT,EvWeight);
      hMTmuFineBin[CutIndex]->Fill(MT,EvWeight);



      if (JetsV.size()>0){
	double dPhiJL=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), muV.Px(),  muV.Py() );
	hdPhiJ0Mu[CutIndex]->Fill(dPhiJL,EvWeight);
				
	double deta = muV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Mu[CutIndex]->Fill(deta,EvWeight);

	double deta2 = tauV.Eta()-JetsV.at(0).Eta();
	hdEtaJ0Tau[CutIndex]->Fill(deta2,EvWeight);
      }

    }////////////WJetsmu
    ///////////End of 
    ////////////////////////// fill generic variables
    if (eIndex != -1)
      {

	hElneutralHadIso[CutIndex]->Fill(el_neutralHadIso[eIndex],EvWeight);
	hElphotonIso[CutIndex]->Fill(el_photonIso[eIndex],EvWeight);
	hElchargedHadIso[CutIndex]->Fill(el_chargedHadIso[eIndex],EvWeight);
	hElpuIso[CutIndex]->Fill(el_puIso[eIndex],EvWeight);
	hElneutralIso[CutIndex]->Fill(el_neutralIso[eIndex],EvWeight);
	hElabsIsoEl[CutIndex]->Fill(el_absIsoEl[eIndex],EvWeight);
	hElrelIsoEl[CutIndex]->Fill(el_relIsoEl[eIndex],EvWeight);

	hElpt[CutIndex]->Fill(el_pt[eIndex],EvWeight);
	hEleta[CutIndex]->Fill(el_eta[eIndex],EvWeight);
	hEldxy[CutIndex]->Fill(el_dxy[eIndex],EvWeight);
	hEldz[CutIndex]->Fill(el_dz[eIndex],EvWeight);
	hEldxyerr[CutIndex]->Fill(el_dxyerr[eIndex],EvWeight);
	hEldzerr[CutIndex]->Fill(el_dzerr[eIndex],EvWeight);

	hElIPsigxy[CutIndex]->Fill(fabs(el_dxy[eIndex]/el_dxyerr[eIndex]),EvWeight);
	hElIPsigz[CutIndex]->Fill(fabs(el_dz[eIndex]/el_dzerr[eIndex]),EvWeight);

	sumElpT +=elV.Pt();
	hIsoEl[CutIndex]->Fill(el_relIso[0],EvWeight);

      }

    if (mIndex !=-1 )
      {




	hMuneutralHadIso[CutIndex]->Fill(mu_neutralHadIso[mIndex],EvWeight);
	hMuphotonIso[CutIndex]->Fill(mu_photonIso[mIndex],EvWeight);
	hMuchargedHadIso[CutIndex]->Fill(mu_chargedHadIso[mIndex],EvWeight);
	hMupuIso[CutIndex]->Fill(mu_puIso[mIndex],EvWeight);
	hMuneutralIso[CutIndex]->Fill(mu_neutralIso[mIndex],EvWeight);
	hMuabsIsoMu[CutIndex]->Fill(mu_absIsoMu[mIndex],EvWeight);
	hMurelIsoMu[CutIndex]->Fill(mu_relIsoMu[mIndex],EvWeight);

	hMupt[CutIndex]->Fill(mu_pt[mIndex],EvWeight);
	hMueta[CutIndex]->Fill(mu_eta[mIndex],EvWeight);
	hMudxy[CutIndex]->Fill(mu_dxy[mIndex],EvWeight);
	hMudz[CutIndex]->Fill(mu_dz[mIndex],EvWeight);
	hMudxyerr[CutIndex]->Fill(mu_dxyerr[mIndex],EvWeight);
	hMudzerr[CutIndex]->Fill(mu_dzerr[mIndex],EvWeight);

	hMuIPsigxy[CutIndex]->Fill(fabs(mu_dxy[mIndex]/mu_dxyerr[mIndex]),EvWeight);
	hMuIPsigz[CutIndex]->Fill(fabs(mu_dz[mIndex]/mu_dzerr[mIndex]),EvWeight);

	sumMuonpT +=muV.Pt();
	hIsoMu[CutIndex]->Fill(mu_relIso[0],EvWeight);


      }


    if (tIndex >-1 && Sel !="muel")
      {
	//cout<< "  taus "<<tIndex<<endl;
	hTaupt[CutIndex]->Fill(ta_pt[tIndex],EvWeight);
	hTaueta[CutIndex]->Fill(ta_eta[tIndex],EvWeight);
	htau_dxy[CutIndex]->Fill(ta_dxy[tIndex],EvWeight);
	htau_dz[CutIndex]->Fill(ta_dz[tIndex],EvWeight);

	double dPhiJ=dPhiFrom2P( ta_px[mIndex], ta_py[mIndex], metx,  mety );
	hdPhiTauMET[CutIndex]->Fill(dPhiJ,EvWeight);

	sumTaupT +=muV.Pt();
	sumTaupT +=elV.Pt();
	hIsoTau[CutIndex]->Fill(ta_relIso[0],EvWeight);
        if (genTauDecayMode1>-1)  hTauDecayMode1[CutIndex]->Fill(genTauDecayMode1,EvWeight);
        if (genTauDecayMode2>-1)  hTauDecayMode2[CutIndex]->Fill(genTauDecayMode2,EvWeight);
	hTauDecayModeAll[CutIndex]->Fill(genTauDecayMode1,EvWeight);
	hTauDecayModeAll[CutIndex]->Fill(genTauDecayMode2,EvWeight);
      }


    if (JetsV.size()>0){

      double dPhiJ=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), metx,  mety );
      hdPhiJ0MET[CutIndex]->Fill(dPhiJ,EvWeight);
      hCosdPhiJ0MET[CutIndex]->Fill(TMath::Cos(dPhiJ),EvWeight);

      hPtJ0[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);


    }

    if (JetsV.size()>1){
      hPtJ1[CutIndex]->Fill(JetsV.at(1).Pt(),EvWeight);
      hHT2[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(1).Px(), JetsV.at(1).Py(), metx,  mety );
      hdPhiJ1MET[CutIndex]->Fill(dPhiJ,EvWeight);

      double dPhiJ0J1=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), JetsV.at(1).Px(), JetsV.at(1).Py() );
      hdPhiJ0J1[CutIndex]->Fill(dPhiJ0J1,EvWeight);


    }
    if (JetsV.size()>2){
      hPtJ2[CutIndex]->Fill(JetsV.at(2).Pt(),EvWeight);
      hHT3[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(2).Px(), JetsV.at(2).Py(), metx,  mety );
      hdPhiJ2MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }
    if (JetsV.size()>3){
      hPtJ3[CutIndex]->Fill(JetsV.at(3).Pt(),EvWeight);
      hHT4[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt()+JetsV.at(3).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(3).Px(), JetsV.at(3).Py(), metx,  mety );
      hdPhiJ3MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }

    for (unsigned int ij=0;ij<JetsV.size();ij++){
      //sumpT+=jet_pt[ij];
      sumpT+=JetsV.at(ij).Pt();
      //double dPhiJ=dPhiFrom2P( jet_px[ij], jet_py[ij], metx,  mety );
      double dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), metx,  mety );
      hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);
    }

    hnBJet[CutIndex]->Fill(nbtag,EvWeight);

    hHT[CutIndex]->Fill(sumpT,EvWeight);

    double HText =sumpT;
    if (Sel=="mutau" && mIndex !=-1)    HText += mu_pt[mIndex];
    if (Sel=="eltau" && eIndex !=-1)    HText += elV.Pt();
    if (Sel=="muel" && eIndex !=-1 &&  mIndex !=-1)     {HText += elV.Pt();HText += mu_pt[mIndex];}

    hHText[CutIndex]->Fill(HText,EvWeight);

    if (Sel=="mutau" && mIndex !=-1 && JetsV.size()>0) {
      hRht[CutIndex]->Fill(mu_pt[mIndex]/HText,EvWeight);
      hPtOHT[CutIndex]->Fill(  mu_pt[mIndex]/sumpT,EvWeight);
    }
    if (Sel=="eltau" && eIndex !=-1 && JetsV.size()>0) {
      hRht[CutIndex]->Fill(elV.Pt()/HText,EvWeight);
      hPtOHT[CutIndex]->Fill(  elV.Pt()/sumpT,EvWeight);

    }
    if (Sel=="muel" && eIndex !=-1 && mIndex !=-1 && JetsV.size()>0) {
      hRht[CutIndex]->Fill( (elV.Pt()+mu_pt[mIndex])/HText,EvWeight);
      hPtOHT[CutIndex]->Fill(  (elV.Pt()+mu_pt[mIndex])/sumpT,EvWeight);

    }

    /////////////global variables

    hMeffMuon[CutIndex]->Fill(sumMuonpT+sumpT+ met,EvWeight);
    hMeffEl[CutIndex]->Fill(sumElpT+sumpT+ met,EvWeight);
    hMeff[CutIndex]->Fill(sumTaupT+sumpT+sumMuonpT+ met,EvWeight);
    if (met>0. ) {

      hMET[CutIndex]->Fill(met,EvWeight);
      hMETFB[CutIndex]->Fill(met,EvWeight);
      hMETphi[CutIndex]->Fill(metphi,EvWeight);
      hHTOsqrMET[CutIndex]->Fill(  sumpT/sqrt(met),EvWeight);
      hMeffMuonOsqrMET[CutIndex]->Fill( (sumMuonpT+sumpT+ met)/sqrt(met),EvWeight);
      hMeffElOsqrMET[CutIndex]->Fill( (sumElpT+sumpT+ met)/sqrt(met),EvWeight);
      hMeffOsqrMET[CutIndex]->Fill( (sumTaupT+sumpT+ met)/sqrt(met),EvWeight);
    }

      hGenMETFB[CutIndex]->Fill(genmet,EvWeight);

  }


/////////////////////////////////////////////////////////////////////////////
  int GetSRIndex(TLorentzVector  muV, TLorentzVector tauV, vector<TLorentzVector>  &JetsV, TLorentzVector  &MetV, double Chimass, int  mIndex, int  tIndex, int CutIndex,bool dogenMET){

    int binIndex=-1;
			
    if (mIndex >-1 && tIndex>-1 ){

    Mt2as = -1.;
    met = MetV.Pt();
    metx = MetV.Px();
    mety = MetV.Py();
    metphi = MetV.Phi();
    dzeta=-9999;

if (dogenMET) { met = sqrt(genmet_Ex*genmet_Ex + genmet_Ey*genmet_Ey);
	metx = genmet_Ex;
	mety = genmet_Ey;
	metphi = genmetphi;
	}

      tauMass = tauV.M();

      Mt2as = asymm_mt2_lester_bisect::get_mT2(muonMass, muV.Px(), muV.Py(),tauMass,tauV.Px(),tauV.Py(),metx,mety,Chimass,Chimass,0);

      float tauUnitX = tauV.Px()/tauV.Pt();
      float tauUnitY = tauV.Py()/tauV.Pt();

      float muonUnitX = muV.Px()/muV.Pt();
      float muonUnitY = muV.Py()/muV.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + muV.Px() + tauV.Px();
      float vectorY = mety + muV.Py() + tauV.Py();

      float vectorVisX = muV.Px() + tauV.Px();
      float vectorVisY = muV.Py() + tauV.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;



    /////////////////////////////////////////////1D distributions

    int nB = 0;
      if (JetsV.size()==0) nB=0; 
      if (JetsV.size()==1) nB=24; 
    if (dzeta > -500 && met>0 && Mt2as>0 && JetsV.size()==0){
	

	if (met<40){
				
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+1;}
	    if (dzeta>-100 && dzeta<0) {binIndex=nB+2;}
	    if (dzeta>0) {binIndex=nB+3;}
	  }
	  if (Mt2as>40) 	{
	    {binIndex=nB+4;}
	  }


	}//////met < 40 bin
	if (met>40 && met<80){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+5;}
	    if (dzeta>-100 && dzeta<50) {binIndex=nB+6;}
	    if (dzeta>50) {binIndex=nB+7;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-100) {binIndex=nB+8;}
	    if (dzeta>-100) {binIndex=nB+9;}
	  }

	  if (Mt2as>80 ) 	{
	    {binIndex=nB+10;}
	  }


	}///met 40-80 bin

	if (met>80 && met<120){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+11;};
	    if (dzeta>-100) {binIndex=nB+12;};
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {binIndex=nB+13;};
	    if (dzeta>-150 ) {binIndex=nB+14;};
	  }

	  if (Mt2as>80 ) 	{
	    {binIndex=nB+15;};
	  }


	}///met 80-120 bin

	if (met>120 && met < 250){
							
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+16;};
	    if (dzeta>-100) {binIndex=nB+17;};
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {binIndex=nB+18;};
	    if (dzeta>-150 && dzeta<-100) {binIndex=nB+19;};
	    if (dzeta>-100) {binIndex=nB+20;};
	  }

	  if (Mt2as>80 && Mt2as<100) 	{
	    binIndex=nB+21;
	  }
	  if (Mt2as>100 && Mt2as<120) 	{
	    binIndex=nB+22;
	  }
	  if (Mt2as>120) 	{
	    binIndex=nB+23;
	  }


	}///met 120 bin & < 250
	if (met>250){
							
	  if (Mt2as>80 && Mt2as<100) 	{
	    binIndex=nB+24;
	  }

	}///met 250 bin
	


    }/////////dzeta>-9999



    if (dzeta > -500 && met>0 && Mt2as>0 && JetsV.size()==1){
	

	if (met<40){
				
	  if (Mt2as<40) 	{
	    if (dzeta<-150) {binIndex=nB+1;}
	    if (dzeta>-150 && dzeta<-100) {binIndex=nB+2;}
	    if (dzeta>-100 && dzeta<0) {binIndex=nB+3;}
	    if (dzeta>0) {binIndex=nB+4;}
	  }
	  if (Mt2as>40) 	{
	    {binIndex=nB+5;}
	  }


	}//////met < 40 bin
	if (met>40 && met<80){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+6;}
	    if (dzeta>-100 && dzeta<50) {binIndex=nB+7;}
	    if (dzeta>50) {binIndex=nB+8;}
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-100) {binIndex=nB+9;}
	    if (dzeta>-100) {binIndex=nB+10;}
	  }

	  if (Mt2as>80 ) 	{
	    {binIndex=nB+11;}
	  }


	}///met 40-80 bin

	if (met>80 && met<120){
	  if (Mt2as<40) 	{
	    if (dzeta<-100) {binIndex=nB+12;};
	    if (dzeta>-100) {binIndex=nB+13;};
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {binIndex=nB+14;};
	    if (dzeta>-150 ) {binIndex=nB+15;};
	  }

	  if (Mt2as>80 && Mt2as<120) 	{
	    binIndex=nB+16;
	  }
	  if (Mt2as>120) 	{
	    binIndex=nB+17;
	  }


	}///met 80-120 bin

	if (met>120 && met < 250){
							
	  if (Mt2as<40) 	{
	    if (dzeta<-150) {binIndex=nB+18;};
	    if (dzeta>-150 && dzeta<-100) {binIndex=nB+19;};
	    if (dzeta>-100) {binIndex=nB+20;};
	  }
	  if (Mt2as>40 && Mt2as<80) 	{
	    if (dzeta<-150) {binIndex=nB+21;};
	    if (dzeta>-150 && dzeta<-100) {binIndex=nB+22;};
	    if (dzeta>-100) {binIndex=nB+23;};
	  }

	  if (Mt2as>80 && Mt2as<100) 	{
	    binIndex=nB+24;
	  }
	  if (Mt2as>100 && Mt2as<120) 	{
	    binIndex=nB+25;
	  }
	  if (Mt2as>120) 	{
	    binIndex=nB+26;
	  }


	}///met 120 bin & < 250
	if (met>250){
							
	  if (Mt2as>80 && Mt2as<100)  binIndex=nB+27;
	  if (Mt2as>100 && Mt2as<120) binIndex=nB+28;
	  if (Mt2as>120) binIndex=nB+29;

	}///met 250 bin
	


    }/////////dzeta>-9999









//    cout<<" inside founder  MT2 "<<Mt2as<<" met "<<met<<" DZ "<<dzeta<<"  m "<<mIndex<<" t "<<tIndex<<" Binindex "<<binIndex<<" Jets "<<JetsV.size()<<"  CutIndex  "<<CutIndex<<endl;

    }

    return binIndex;


  }

///////////////////////////////////////////////////////////////



  void FillHistsDiL(int CutIndex, Double_t EvWeight, TLorentzVector  LeptV1, TLorentzVector  LeptV2,  vector<TLorentzVector>  &JetsV, TLorentzVector  &MetV, double Chimass, double mintermediate, string & Sel, int  mIndex_1, int mIndex_2){
    AllJets_Lepton_noMet.clear();
    for (unsigned int i = 0; i <   JetsV.size(); i++) AllJets_Lepton_noMet.push_back (JetsV.at(i));


    if (Sel=="mumu" )  {AllJets_Lepton_noMet.push_back (LeptV1);AllJets_Lepton_noMet.push_back (LeptV2);};

    if (mIndex_1 >-1 && mIndex_2 >-1){
      double dPhiDil=dPhiFrom2P(LeptV1.Px(),  LeptV1.Py(),LeptV2.Px(),  LeptV2.Py() );
      hdPhiDil[CutIndex]->Fill(dPhiDil,EvWeight);
      double detal=deltaEta(LeptV1.Px(),  LeptV1.Py(),LeptV1.Pz(), LeptV2.Px(),  LeptV2.Py(),LeptV2.Pz() );
      hdEtaDil[CutIndex]->Fill(detal,EvWeight);
    }

    sumpT=0;

    double sumMuonpT=0;
    double sumElpT=0;
    double sumTaupT=0;  

    hnJet[CutIndex]->Fill(njets,EvWeight);
    hnpartons[CutIndex]->Fill(npartons,EvWeight);
    hnMu[CutIndex]->Fill(mu_count,EvWeight);
    hnTau[CutIndex]->Fill(ta_count,EvWeight);
    hnEl[CutIndex]->Fill(el_count,EvWeight);

    hWeights[CutIndex]->Fill(EvWeight);
    hnpv[CutIndex]->Fill(npv,EvWeight);
    hnpu[CutIndex]->Fill(npu,EvWeight);
    hEventSign[CutIndex]->Fill(event_sign,EvWeight);
	


    if (Sel=="mumu"){
           		
      Mt2as = -1.;
      met = MetV.Pt();
      dzeta=-9999;
				

      double v1[4] = {LeptV1.Px(),LeptV1.Py(),LeptV1.Pz(),muonMass};
      double v2[4] = {LeptV2.Px(),LeptV2.Py(),LeptV2.Pz(),muonMass};
      double ecm = 13000.0;
      double mxlo = Chimass;
      double ptm[2] = {metx,mety};
      double vds[4] = {0,0,0,0};
      // TMctLib t;
      Mt2as = asymm_mt2_lester_bisect::get_mT2(muonMass, LeptV1.Px(), LeptV1.Py(), muonMass, LeptV2.Px(), LeptV2.Py(),metx,mety,Chimass,Chimass,0);
      double mcta = sqrt( 2* LeptV1.Pt()* LeptV2.Pt()*(1+cos( LeptV1.Phi()- LeptV2.Phi())) );
				
				 
      float MT2v =0 ; MT2v = mt2(v1,v2,vds,ptm,ecm,mxlo);

      hMt2Dil[CutIndex]->Fill(MT2v,EvWeight);
      hMt2lesterDil[CutIndex]->Fill(Mt2as,EvWeight);
      hMt2lesterDilFB[CutIndex]->Fill(Mt2as,EvWeight);
      hMCTDil[CutIndex]->Fill(mct(v1,v2),EvWeight);
      hMCTbDil[CutIndex]->Fill(mcta,EvWeight);

      double Centr = Centrality(AllJets_Lepton_noMet);
      hCentrality[CutIndex]->Fill(Centr,EvWeight);
      double Dr =  LeptV1.DeltaR( LeptV2);
      //
      TLorentzVector DiL =  LeptV1  + LeptV2;
      double dPhi=dPhiFrom2P( DiL.Px(), DiL.Py(), metx,  mety );
      float MT = TMath::Sqrt(2*DiL.Pt()*met*(1-TMath::Cos(dPhi)));
      hMTDil[CutIndex]->Fill(MT,EvWeight);
      hInvMassDil[CutIndex]->Fill(DiL.M(),EvWeight);
      hInvMassDilFineBin[CutIndex]->Fill(DiL.M(),EvWeight);
      hdR_Dil[CutIndex]->Fill(Dr,EvWeight);
      hPtDil[CutIndex]->Fill(DiL.Pt(),EvWeight);


      dPhi=dPhiFrom2P( mu_px[mIndex_1], mu_py[mIndex_1], metx,  mety );
      hdPhiLept1MET[CutIndex]->Fill(dPhi,EvWeight);
      float MTlept1 = TMath::Sqrt(2*mu_pt[mIndex_1]*met*(1-TMath::Cos(dPhi)));
      hMTlept1[CutIndex]->Fill(MTlept1,EvWeight);


      dPhi=dPhiFrom2P( mu_px[mIndex_2], mu_py[mIndex_2], metx,  mety );
      hdPhiLept2MET[CutIndex]->Fill(dPhi,EvWeight);
      float MTlept2 = TMath::Sqrt(2*mu_pt[mIndex_2]*met*(1-TMath::Cos(dPhi)));
      hMTlept2[CutIndex]->Fill(MTlept2,EvWeight);
					
      hMTsum[CutIndex]->Fill(MTlept1+MTlept2,EvWeight);
      float mttotal = sqrt(MTlept1*MTlept1 + MTlept2*MTlept2);
      hMTtot[CutIndex]->Fill(mttotal,EvWeight);


      float tauUnitX = LeptV2.Px()/LeptV2.Pt();
      float tauUnitY = LeptV2.Py()/LeptV2.Pt();

      //	cout<<" CHECK =========== "<<tauV.Pt()<<"  "<<ta_pt[tIndex]<<endl;	
      float muonUnitX = LeptV1.Px()/LeptV1.Pt();
      float muonUnitY = LeptV1.Py()/LeptV1.Pt();

      float zetaX = tauUnitX + muonUnitX;
      float zetaY = tauUnitY + muonUnitY;

      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = metx + LeptV1.Px() + LeptV2.Px();
      float vectorY = mety + LeptV1.Py() + LeptV2.Py();

      float vectorVisX = LeptV1.Px() + LeptV2.Px();
      float vectorVisY = LeptV1.Py() + LeptV2.Py();

      // computation of DZeta variable
      // pfmet
      pzetamiss = metx*zetaX + mety*zetaY;
      pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
      dzeta = pzetamiss - 0.85*pzetavis;
			

      hDZeta[CutIndex]->Fill(dzeta,EvWeight);
      hDZetaFB[CutIndex]->Fill(dzeta,EvWeight);

      if (JetsV.size()>0){
	double dPhiJT=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), LeptV1.Px(),  LeptV1.Py() );
	hdPhiJ0Lept1[CutIndex]->Fill(dPhiJT,EvWeight);
	dPhiJT=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), LeptV2.Px(),  LeptV2.Py() );
	hdPhiJ0Lept2[CutIndex]->Fill(dPhiJT,EvWeight);
      }

    }///mumu sel


    if (mIndex_1 !=-1 )
      {


	int mIndex = mIndex_1;
	hLept1neutralHadIso[CutIndex]->Fill(mu_neutralHadIso[mIndex],EvWeight);
	hLept1photonIso[CutIndex]->Fill(mu_photonIso[mIndex],EvWeight);
	hLept1chargedHadIso[CutIndex]->Fill(mu_chargedHadIso[mIndex],EvWeight);
	hLept1puIso[CutIndex]->Fill(mu_puIso[mIndex],EvWeight);
	hLept1neutralIso[CutIndex]->Fill(mu_neutralIso[mIndex],EvWeight);
	hLept1absIsoMu[CutIndex]->Fill(mu_absIsoMu[mIndex],EvWeight);
	hLept1relIsoMu[CutIndex]->Fill(mu_relIsoMu[mIndex],EvWeight);

	hLept1pt[CutIndex]->Fill(mu_pt[mIndex],EvWeight);
	hLept1eta[CutIndex]->Fill(mu_eta[mIndex],EvWeight);
	hLept1dxy[CutIndex]->Fill(mu_dxy[mIndex],EvWeight);
	hLept1dz[CutIndex]->Fill(mu_dz[mIndex],EvWeight);
	hLept1dxyerr[CutIndex]->Fill(mu_dxyerr[mIndex],EvWeight);
	hLept1dzerr[CutIndex]->Fill(mu_dzerr[mIndex],EvWeight);

	hLept1IPsigxy[CutIndex]->Fill(fabs(mu_dxy[mIndex]/mu_dxyerr[mIndex]),EvWeight);
	hLept1IPsigz[CutIndex]->Fill(fabs(mu_dz[mIndex]/mu_dzerr[mIndex]),EvWeight);

	double dPhiJ=dPhiFrom2P( mu_px[mIndex], mu_py[mIndex], metx,  mety );
	hdPhiLept1MET[CutIndex]->Fill(dPhiJ,EvWeight);


	sumMuonpT +=LeptV1.Pt();
	hIsoLept1[CutIndex]->Fill(mu_relIso[0],EvWeight);


	//cout<<" Inside  "<<relIsoL<<"  Cut"<<"  "<<CutIndex<<endl;


      }

    if (mIndex_2 !=-1 )
      {

	int mIndex = mIndex_2;
	hLept2neutralHadIso[CutIndex]->Fill(mu_neutralHadIso[mIndex],EvWeight);
	hLept2photonIso[CutIndex]->Fill(mu_photonIso[mIndex],EvWeight);
	hLept2chargedHadIso[CutIndex]->Fill(mu_chargedHadIso[mIndex],EvWeight);
	hLept2puIso[CutIndex]->Fill(mu_puIso[mIndex],EvWeight);
	hLept2neutralIso[CutIndex]->Fill(mu_neutralIso[mIndex],EvWeight);
	hLept2absIsoMu[CutIndex]->Fill(mu_absIsoMu[mIndex],EvWeight);
	hLept2relIsoMu[CutIndex]->Fill(mu_relIsoMu[mIndex],EvWeight);

	hLept2pt[CutIndex]->Fill(mu_pt[mIndex],EvWeight);
	hLept2eta[CutIndex]->Fill(mu_eta[mIndex],EvWeight);
	hLept2dxy[CutIndex]->Fill(mu_dxy[mIndex],EvWeight);
	hLept2dz[CutIndex]->Fill(mu_dz[mIndex],EvWeight);
	hLept2dxyerr[CutIndex]->Fill(mu_dxyerr[mIndex],EvWeight);
	hLept2dzerr[CutIndex]->Fill(mu_dzerr[mIndex],EvWeight);

	hLept2IPsigxy[CutIndex]->Fill(fabs(mu_dxy[mIndex]/mu_dxyerr[mIndex]),EvWeight);
	hLept2IPsigz[CutIndex]->Fill(fabs(mu_dz[mIndex]/mu_dzerr[mIndex]),EvWeight);

	double dPhiJ=dPhiFrom2P( mu_px[mIndex], mu_py[mIndex], metx,  mety );
	hdPhiLept2MET[CutIndex]->Fill(dPhiJ,EvWeight);


	sumMuonpT +=LeptV2.Pt();
	hIsoLept2[CutIndex]->Fill(mu_relIso[1],EvWeight);


	//cout<<" Inside  "<<relIsoL<<"  Cut"<<"  "<<CutIndex<<endl;


      }


    /////////////global variables

    if (JetsV.size()>0){

      double dPhiJ=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), metx,  mety );
      hdPhiJ0MET[CutIndex]->Fill(dPhiJ,EvWeight);
      hCosdPhiJ0MET[CutIndex]->Fill(TMath::Cos(dPhiJ),EvWeight);

      hPtJ0[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);


    }

    if (JetsV.size()>1){
      hPtJ1[CutIndex]->Fill(JetsV.at(1).Pt(),EvWeight);
      hHT2[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(1).Px(), JetsV.at(1).Py(), metx,  mety );
      hdPhiJ1MET[CutIndex]->Fill(dPhiJ,EvWeight);

      double dPhiJ0J1=dPhiFrom2P( JetsV.at(0).Px(), JetsV.at(0).Py(), JetsV.at(1).Px(), JetsV.at(1).Py() );
      hdPhiJ0J1[CutIndex]->Fill(dPhiJ0J1,EvWeight);

      TLorentzVector JetMass = JetsV.at(0) + JetsV.at(1) ;
      hDiJetMass_J0J1[CutIndex]->Fill(JetMass.M(),EvWeight);
      double detaj=deltaEta(JetsV.at(0).Px(),  JetsV.at(0).Py(), JetsV.at(0).Pz(), JetsV.at(1).Px(),  JetsV.at(1).Py(), JetsV.at(1).Pz() );
      //double detaj = JetsV.at(0).Eta() - JetsV.at(1).Eta();
      hdEtaJ0J1[CutIndex]->Fill(detaj,EvWeight);

      //cout<<" Jet1 properties Eta  "<<JetsV.at(0).Eta()<<"  "<<JetsV.at(1).Eta()<<" Phi "<<JetsV.at(0).Phi()<<"   "<<JetsV.at(1).Phi()<<endl;

    }
    if (JetsV.size()>2){
      hPtJ2[CutIndex]->Fill(JetsV.at(2).Pt(),EvWeight);
      hHT3[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(2).Px(), JetsV.at(2).Py(), metx,  mety );
      hdPhiJ2MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }
    if (JetsV.size()>3){
      hPtJ3[CutIndex]->Fill(JetsV.at(3).Pt(),EvWeight);
      hHT4[CutIndex]->Fill(JetsV.at(0).Pt()+JetsV.at(1).Pt()+JetsV.at(2).Pt()+JetsV.at(3).Pt(),EvWeight);
      double dPhiJ=dPhiFrom2P( JetsV.at(3).Px(), JetsV.at(3).Py(), metx,  mety );
      hdPhiJ3MET[CutIndex]->Fill(dPhiJ,EvWeight);
    }

    for (unsigned int ij=0;ij<JetsV.size();ij++){
      //sumpT+=jet_pt[ij];
      sumpT+=JetsV.at(ij).Pt();
      //double dPhiJ=dPhiFrom2P( jet_px[ij], jet_py[ij], metx,  mety );
      double dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), metx,  mety );
      hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);
    }

    hnBJet[CutIndex]->Fill(nbtag,EvWeight);

    hHT[CutIndex]->Fill(sumpT,EvWeight);

    double HText =sumpT;
    if (Sel=="mumu" && mIndex_1 !=-1 &&  mIndex_2 !=-1)     {HText += LeptV1.Pt();HText += LeptV2.Pt();}

    hHText[CutIndex]->Fill(HText,EvWeight);



    hMeffMuon[CutIndex]->Fill(sumMuonpT+sumpT+ met,EvWeight);
    if (met>0. ) {

      hMET[CutIndex]->Fill(met,EvWeight);
      hMETFB[CutIndex]->Fill(met,EvWeight);
      hMETphi[CutIndex]->Fill(metphi,EvWeight);
      hHTOsqrMET[CutIndex]->Fill(  sumpT/sqrt(met),EvWeight);
      hMeffMuonOsqrMET[CutIndex]->Fill( (sumMuonpT+sumpT+ met)/sqrt(met),EvWeight);
    }

      hGenMETFB[CutIndex]->Fill(genmet,EvWeight);

  }





};

#endif

#ifdef analyzer_cxx
analyzer::analyzer(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {

#ifdef SINGLE_TREE
    // The following code should be used if you want this class to access
    // a single tree instead of a chain
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f || !f->IsOpen()) {
      f = new TFile("Memory Directory");
    }
    f->GetObject("muel/T",tree);

#else // SINGLE_TREE

    // The following code should be used if you want this class to access a chain
    // of trees.
    TChain * chain = new TChain("muel/T","");

    TString ChainName = "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/muel/FILEIN/muel/T";
    TString Tsystematic="SYSTEMATICHERE";
    if (Tsystematic == "JetEnUp" || Tsystematic == "JetEnDown" || Tsystematic == "UnclEnDown" || Tsystematic == "UnclEnDown" || Tsystematic == "TauEnUp" || Tsystematic == "TauEnDown" || Tsystematic == "ElEnUp" || Tsystematic == "ElEnDown" || Tsystematic == "MuEnUp" || Tsystematic == "MuEnDown" || Tsystematic == "BTagUp" || Tsystematic == "BTagDown") 
      ChainName ="/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/muel_"+Tsystematic+"/FILEIN/muel/T";				

    cout<<"  Will now chain "<<ChainName<<endl;
    chain->Add(ChainName);

    cout<<"  chained alright..  "<<ChainName<<endl;
    tree = chain;
 //   tree->Print();
 //   chain->Print();
#endif // SINGLE_TREE

  }
  Init(tree);
}

analyzer::~analyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t analyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t analyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}







void analyzer::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  //		tree->SetMaxVirtualSize(1000000);
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("met_ex", &met_ex, &b_met_ex);
  fChain->SetBranchAddress("met_ey", &met_ey, &b_met_ey);
  fChain->SetBranchAddress("met_ez", &met_ez, &b_met_ez);
  fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);

  fChain->SetBranchAddress("met_ex_recoil", &met_ex_recoil, &b_met_ex_recoil);
  fChain->SetBranchAddress("met_ey_recoil", &met_ey_recoil, &b_met_ey_recoil);
  fChain->SetBranchAddress("genmet", &genmet, &b_genmet);
  fChain->SetBranchAddress("genmetphi", &genmetphi, &b_genmetphi);
  fChain->SetBranchAddress("met_scaleUp", &met_scaleUp, &b_met_scaleUp);
  fChain->SetBranchAddress("met_scaleDown", &met_scaleDown, &b_met_scaleDown);
  fChain->SetBranchAddress("metphi_scaleUp", &metphi_scaleUp, &b_metphi_scaleUp);
  fChain->SetBranchAddress("metphi_scaleDown", &metphi_scaleDown, &b_metphi_scaleDown);

  fChain->SetBranchAddress("genmet_Ex", &genmet_Ex, &b_genmet_Ex);
  fChain->SetBranchAddress("genmet_Ey", &genmet_Ey, &b_genmet_Ey);
  fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
  fChain->SetBranchAddress("event_run", &event_run, &b_event_run);
  fChain->SetBranchAddress("NuPx", &NuPx, &b_NuPx);
  fChain->SetBranchAddress("NuPy", &NuPy, &b_NuPy);
  fChain->SetBranchAddress("NuPz", &NuPz, &b_NuPz);
  fChain->SetBranchAddress("NuPt", &NuPt, &b_NuPt);
  fChain->SetBranchAddress("NuPhi", &NuPhi, &b_NuPhi);

  fChain->SetBranchAddress("met_resoUp", &met_resoUp, &b_met_resoUp);
  fChain->SetBranchAddress("met_resoDown", &met_resoDown, &b_met_resoDown);
  fChain->SetBranchAddress("metphi_resoUp", &metphi_resoUp, &b_metphi_resoUp);
  fChain->SetBranchAddress("metphi_resoDown", &metphi_resoDown, &b_metphi_resoDown);

  fChain->SetBranchAddress("met_ex_JetEnUp_recoil", &met_ex_JetEnUp_recoil, &b_met_ex_JetEnUp_recoil);
  fChain->SetBranchAddress("met_ey_JetEnUp_recoil", &met_ey_JetEnUp_recoil, &b_met_ey_JetEnUp_recoil);

  fChain->SetBranchAddress("met_ex_JetEnDown_recoil", &met_ex_JetEnDown_recoil, &b_met_ex_JetEnDown_recoil);
  fChain->SetBranchAddress("met_ey_JetEnDown_recoil", &met_ey_JetEnDown_recoil, &b_met_ey_JetEnDown_recoil);

  fChain->SetBranchAddress("met_ex_UnclusteredEnDown_recoil", &met_ex_UnclusteredEnDown_recoil, &b_met_ex_UnclusteredEnDown_recoil);
  fChain->SetBranchAddress("met_ey_UnclusteredEnDown_recoil", &met_ey_UnclusteredEnDown_recoil, &b_met_ey_UnclusteredEnDown_recoil);

  fChain->SetBranchAddress("met_ex_UnclusteredEnUp_recoil", &met_ex_UnclusteredEnUp_recoil, &b_met_ex_UnclusteredEnUp_recoil);
  fChain->SetBranchAddress("met_ey_UnclusteredEnUp_recoil", &met_ey_UnclusteredEnUp_recoil, &b_met_ey_UnclusteredEnUp_recoil);

  fChain->SetBranchAddress("met_ex_JetEnUp", &met_ex_JetEnUp, &b_met_ex_JetEnUp);
  fChain->SetBranchAddress("met_ey_JetEnUp", &met_ey_JetEnUp, &b_met_ey_JetEnUp);

  fChain->SetBranchAddress("met_ex_JetEnDown", &met_ex_JetEnDown, &b_met_ex_JetEnDown);
  fChain->SetBranchAddress("met_ey_JetEnDown", &met_ey_JetEnDown, &b_met_ey_JetEnDown);

  fChain->SetBranchAddress("met_ex_UnclusteredEnDown", &met_ex_UnclusteredEnDown, &b_met_ex_UnclusteredEnDown);
  fChain->SetBranchAddress("met_ey_UnclusteredEnDown", &met_ey_UnclusteredEnDown, &b_met_ey_UnclusteredEnDown);

  fChain->SetBranchAddress("met_ex_UnclusteredEnUp", &met_ex_UnclusteredEnUp, &b_met_ex_UnclusteredEnUp);
  fChain->SetBranchAddress("met_ey_UnclusteredEnUp", &met_ey_UnclusteredEnUp, &b_met_ey_UnclusteredEnUp);


  fChain->SetBranchAddress("wScale0", &wScale0, &b_wScale0);
  fChain->SetBranchAddress("wScale1", &wScale1, &b_wScale1);
  fChain->SetBranchAddress("wScale2", &wScale2, &b_wScale2);
  fChain->SetBranchAddress("wScale3", &wScale3, &b_wScale3);
  fChain->SetBranchAddress("wScale4", &wScale4, &b_wScale4);
  fChain->SetBranchAddress("wScale5", &wScale5, &b_wScale5);
  fChain->SetBranchAddress("wScale6", &wScale6, &b_wScale6);
  fChain->SetBranchAddress("wScale7", &wScale7, &b_wScale7);
  fChain->SetBranchAddress("wScale8", &wScale8, &b_wScale8);
  fChain->SetBranchAddress("wPDFUp", &wPDFUp, &b_wPDFUp);
  fChain->SetBranchAddress("wPDFDown", &wPDFDown, &b_wPDFDown);
  fChain->SetBranchAddress("gen_weight", &gen_weight, &b_gen_weight);
  fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  fChain->SetBranchAddress("LSF_weight", &LSF_weight, &b_LSF_weight);
  fChain->SetBranchAddress("LSF_weight_mu", &LSF_weight_mu, &b_LSF_weight_mu);
  fChain->SetBranchAddress("LSF_weight_el", &LSF_weight_el, &b_LSF_weight_el);
  fChain->SetBranchAddress("LSF_weight_1", &LSF_weight_1, &b_LSF_weight_1);
  fChain->SetBranchAddress("LSF_weight_2", &LSF_weight_2, &b_LSF_weight_2);
  fChain->SetBranchAddress("TFR_weight", &TFR_weight, &b_TFR_weight);
  fChain->SetBranchAddress("top_weight", &top_weight, &b_top_weight);
  fChain->SetBranchAddress("all_weight", &all_weight, &b_all_weight);
  fChain->SetBranchAddress("trig_weight", &trig_weight, &b_trig_weight);
  fChain->SetBranchAddress("trig_weight_1", &trig_weight_1, &b_trig_weight_1);
  fChain->SetBranchAddress("trig_weight_2", &trig_weight_2, &b_trig_weight_2);
  fChain->SetBranchAddress("zptmassweight", &zptmassweight, &b_zptmassweight);
  fChain->SetBranchAddress("xsecs", &xsecs, &b_xsecs);
  fChain->SetBranchAddress("muon_index", &muon_index, &b_muon_index);
  fChain->SetBranchAddress("muon_index_1", &muon_index_1, &b_muon_index_1);
  fChain->SetBranchAddress("muon_index_2", &muon_index_2, &b_muon_index_2);
  fChain->SetBranchAddress("electron_index", &electron_index, &b_electron_index);
  fChain->SetBranchAddress("taus_index", &taus_index, &b_taus_index);
  fChain->SetBranchAddress("primvert_count", &primvert_count, &b_primvert_count);
  fChain->SetBranchAddress("primvert_x", &primvert_x, &b_primvert_x);
  fChain->SetBranchAddress("primvert_y", &primvert_y, &b_primvert_y);
  fChain->SetBranchAddress("primvert_z", &primvert_z, &b_primvert_z);
  fChain->SetBranchAddress("mu_count", &mu_count, &b_mu_count);
  fChain->SetBranchAddress("mu_px", mu_px, &b_mu_px);
  fChain->SetBranchAddress("mu_py", mu_py, &b_mu_py);
  fChain->SetBranchAddress("mu_pz", mu_pz, &b_mu_pz);
  fChain->SetBranchAddress("mu_pt", mu_pt, &b_mu_pt);
  fChain->SetBranchAddress("mu_eta", mu_eta, &b_mu_eta);
  fChain->SetBranchAddress("mu_phi", mu_phi, &b_mu_phi);
  fChain->SetBranchAddress("mu_charge", mu_charge, &b_mu_charge);
  fChain->SetBranchAddress("mu_miniISO", mu_miniISO, &b_mu_miniISO);
  fChain->SetBranchAddress("mu_dxy", mu_dxy, &b_mu_dxy);
  fChain->SetBranchAddress("mu_dz", mu_dz, &b_mu_dz);
  fChain->SetBranchAddress("mu_dxyerr", mu_dxyerr, &b_mu_dxyerr);
  fChain->SetBranchAddress("mu_dzerr", mu_dzerr, &b_mu_dzerr);

  fChain->SetBranchAddress("mu_neutralHadIso", mu_neutralHadIso, &b_mu_neutralHadIso);
  fChain->SetBranchAddress("mu_photonIso", mu_photonIso, &b_mu_photonIso);
  fChain->SetBranchAddress("mu_chargedHadIso", mu_chargedHadIso, &b_mu_chargedHadIso);
  fChain->SetBranchAddress("mu_puIso", mu_puIso, &b_mu_puIso);
  fChain->SetBranchAddress("mu_neutralIso", mu_neutralIso, &b_mu_neutralIso);
  fChain->SetBranchAddress("mu_absIsoMu", mu_absIsoMu, &b_mu_absIsoMu);
  fChain->SetBranchAddress("mu_relIsoMu", mu_relIsoMu, &b_mu_relIsoMu);
  fChain->SetBranchAddress("el_neutralHadIso", el_neutralHadIso, &b_el_neutralHadIso);
  fChain->SetBranchAddress("el_photonIso", el_photonIso, &b_el_photonIso);
  fChain->SetBranchAddress("el_chargedHadIso", el_chargedHadIso, &b_el_chargedHadIso);
  fChain->SetBranchAddress("el_puIso", el_puIso, &b_el_puIso);
  fChain->SetBranchAddress("el_neutralIso", el_neutralIso, &b_el_neutralIso);
  fChain->SetBranchAddress("el_absIsoEl", el_absIsoEl, &b_el_absIsoEl);
  fChain->SetBranchAddress("el_relIsoEl", el_relIsoEl, &b_el_relIsoEl);
  fChain->SetBranchAddress("el_isMVA", el_isMVA, &b_el_isMVA);

  fChain->SetBranchAddress("mu_relIso", mu_relIso, &b_mu_relIso);
  fChain->SetBranchAddress("jet_count", &jet_count, &b_jet_count);
  fChain->SetBranchAddress("npv", &npv, &b_npv);
  fChain->SetBranchAddress("genTauDecayMode1", &genTauDecayMode1, &b_genTauDecayMode1);
  fChain->SetBranchAddress("genTauDecayMode2", &genTauDecayMode2, &b_genTauDecayMode2);
  fChain->SetBranchAddress("genTauDecayModeAll", &genTauDecayModeAll, &b_genTauDecayModeAll);
  fChain->SetBranchAddress("npu", &npu, &b_npu);
  fChain->SetBranchAddress("jets_cleaned", &jets_cleaned, &b_jets_cleaned);
  fChain->SetBranchAddress("njets", &njets, &b_njets);
  fChain->SetBranchAddress("jet_e", jet_e, &b_jet_e);
  fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
  fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
  fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
  fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  fChain->SetBranchAddress("jet_flavour", jet_flavour, &b_jet_flavour);
  fChain->SetBranchAddress("jet_btag", jet_btag, &b_jet_btag);
  fChain->SetBranchAddress("jet_isLoose", jet_isLoose, &b_jet_isLoose);
  fChain->SetBranchAddress("el_count", &el_count, &b_el_count);
  fChain->SetBranchAddress("el_px", el_px, &b_el_px);
  fChain->SetBranchAddress("el_py", el_py, &b_el_py);
  fChain->SetBranchAddress("el_pz", el_pz, &b_el_pz);
  fChain->SetBranchAddress("el_pt", el_pt, &b_el_pt);
  fChain->SetBranchAddress("el_eta", el_eta, &b_el_eta);
  fChain->SetBranchAddress("el_phi", el_phi, &b_el_phi);
  fChain->SetBranchAddress("el_miniISO", el_miniISO, &b_el_miniISO);
  fChain->SetBranchAddress("el_dxy", el_dxy, &b_el_dxy);
  fChain->SetBranchAddress("el_dz", el_dz, &b_el_dz);
  fChain->SetBranchAddress("el_dxyerr", el_dxyerr, &b_el_dxyerr);
  fChain->SetBranchAddress("el_dzerr", el_dzerr, &b_el_dzerr);
  fChain->SetBranchAddress("el_charge", el_charge, &b_el_charge);
  fChain->SetBranchAddress("el_relIso", el_relIso, &b_el_relIso);
  fChain->SetBranchAddress("ta_count", &ta_count, &b_ta_count);
  fChain->SetBranchAddress("ta_px", ta_px, &b_ta_px);
  fChain->SetBranchAddress("ta_py", ta_py, &b_ta_py);
  fChain->SetBranchAddress("ta_pz", ta_pz, &b_ta_pz);
  fChain->SetBranchAddress("ta_mass", ta_mass, &b_ta_mass);
  fChain->SetBranchAddress("ta_eta", ta_eta, &b_ta_eta);
  fChain->SetBranchAddress("ta_phi", ta_phi, &b_ta_phi);
  fChain->SetBranchAddress("ta_pt", ta_pt, &b_ta_pt);
  fChain->SetBranchAddress("ta_dxy", ta_dxy, &b_ta_dxy);
  fChain->SetBranchAddress("ta_dz", ta_dz, &b_ta_dz);
  fChain->SetBranchAddress("ta_charge", ta_charge, &b_ta_charge);
  fChain->SetBranchAddress("ta_IsoFlag", &ta_IsoFlag, &b_ta_IsoFlag);
  fChain->SetBranchAddress("ta_relIso", ta_relIso, &b_ta_relIso);
  fChain->SetBranchAddress("datasetName", &datasetName);
  fChain->SetBranchAddress("CFCounter_", CFCounter_,&b_CFCounter_);
  fChain->SetBranchAddress("regionName", &regionName, &b_regionName);
  fChain->SetBranchAddress("event_sign", &event_sign, &b_event_sign);
  fChain->SetBranchAddress("met_flag", &met_flag, &b_met_flag);
  fChain->SetBranchAddress("eleMVA", &eleMVA, &b_eleMVA);
  fChain->SetBranchAddress("event_secondLeptonVeto", &event_secondLeptonVeto, &b_event_secondLeptonVeto);
  fChain->SetBranchAddress("event_thirdLeptonVeto", &event_thirdLeptonVeto, &b_event_thirdLeptonVeto);
  fChain->SetBranchAddress("event_leptonDrTrigger", &event_leptonDrTrigger, &b_event_leptonDrTrigger);
  fChain->SetBranchAddress("genTauMatched", &genTauMatched, &b_genTauMatched);
  fChain->SetBranchAddress("genLeptonMatched", &genLeptonMatched, &b_genLeptonMatched);
  fChain->SetBranchAddress("genLeptonMatchedEl", &genLeptonMatchedEl, &b_genLeptonMatchedEl);
  fChain->SetBranchAddress("genLeptonMatchedMu", &genLeptonMatchedMu, &b_genLeptonMatchedMu);
  fChain->SetBranchAddress("genTauDecayedMuMatched", &genTauDecayedMuMatched, &b_genTauDecayedMuMatched);
  fChain->SetBranchAddress("genTauDecayedElMatched", &genTauDecayedElMatched, &b_genTauDecayedElMatched);
  fChain->SetBranchAddress("genLeptonPromptElMatched", &genLeptonPromptElMatched, &b_genLeptonPromptElMatched);
  fChain->SetBranchAddress("genLeptonPromptMuMatched", &genLeptonPromptMuMatched, &b_genLeptonPromptMuMatched);

  fChain->SetBranchAddress("genLeptonMatchedPromptEl", &genLeptonMatchedPromptMu, &b_genLeptonMatchedPromptEl);
  fChain->SetBranchAddress("genLeptonMatchedPromptMu", &genLeptonMatchedPromptMu, &b_genLeptonMatchedPromptMu);
  //			fChain->SetBranchAddress("genLeptonMatchedPromptTau", &genLeptonMatchedPromptMu, &b_genLeptonMatchedPromptTau);

  fChain->SetBranchAddress("genLeptonMatchedPrompEl", &genLeptonMatchedPrompEl, &b_genLeptonMatchedPrompEl);
  fChain->SetBranchAddress("genLeptonMatchedPrompMu", &genLeptonMatchedPrompMu, &b_genLeptonMatchedPrompMu);
  fChain->SetBranchAddress("genElMatchedToTauDecay", &genElMatchedToTauDecay, &b_genElMatchedToTauDecay);
  fChain->SetBranchAddress("genMuMatchedToTauDecay", &genMuMatchedToTauDecay, &b_genMuMatchedToTauDecay);
  fChain->SetBranchAddress("genTauMatchedToTauDecay", &genTauMatchedToTauDecay, &b_genTauMatchedToTauDecay);
  fChain->SetBranchAddress("matchedTauToPromptEl", &matchedTauToPromptEl, &b_matchedTauToPromptEl);
  fChain->SetBranchAddress("matchedTauToPromptMu", &matchedTauToPromptMu, &b_matchedTauToPromptMu);
  fChain->SetBranchAddress("matchedTauToPromptTau", &matchedTauToPromptTau, &b_matchedTauToPromptTau);
  fChain->SetBranchAddress("matchedTauToTauDecEl", &matchedTauToTauDecEl, &b_matchedTauToTauDecEl);
  fChain->SetBranchAddress("matchedTauToTauDecMu", &matchedTauToTauDecMu, &b_matchedTauToTauDecMu);
  fChain->SetBranchAddress("matchedTauToTauDecTau", &matchedTauToTauDecTau, &b_matchedTauToTauDecTau);
  fChain->SetBranchAddress("matchedTauToElHadronDec", &matchedTauToElHadronDec, &b_matchedTauToElHadronDec);
  fChain->SetBranchAddress("matchedTauToMuHadronDec", &matchedTauToMuHadronDec, &b_matchedTauToMuHadronDec);
  fChain->SetBranchAddress("matchedTauToTauHadronDec", &matchedTauToTauHadronDec, &b_matchedTauToTauHadronDec);
  fChain->SetBranchAddress("matchedTauToGluon", &matchedTauToGluon, &b_matchedTauToGluon);
  fChain->SetBranchAddress("matchedTauToLFQ", &matchedTauToLFQ, &b_matchedTauToLFQ);
  fChain->SetBranchAddress("matchedTauToHFQ", &matchedTauToHFQ, &b_matchedTauToHFQ);
  fChain->SetBranchAddress("genLeptonMatchedGluon", &genLeptonMatchedGluon, &b_genLeptonMatchedGluon);
  fChain->SetBranchAddress("genLeptonMatchedLFQ", &genLeptonMatchedLFQ, &b_genLeptonMatchedLFQ);
  fChain->SetBranchAddress("genLeptonMatchedHFQ", &genLeptonMatchedHFQ, &b_genLeptonMatchedHFQ);
  
  fChain->SetBranchAddress("isDYTT", &isDYTT, &b_isDYTT);
  fChain->SetBranchAddress("isDYEE", &isDYEE, &b_isDYEE);
  fChain->SetBranchAddress("isDYMM", &isDYMM, &b_isDYMM);
  fChain->SetBranchAddress("isDYNuNu", &isDYNuNu, &b_isDYNuNu);



			
  fChain->SetBranchAddress("qcdweight", &qcdweight, &b_qcdweight);
  fChain->SetBranchAddress("qcdweightup", &qcdweightup, &b_qcdweightup);
  fChain->SetBranchAddress("qcdweightdown", &qcdweightdown, &b_qcdweightdown);
  fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);

  fChain->SetBranchAddress("npartons",&npartons,&b_npartons);

  Notify();
}

Bool_t analyzer::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}






void analyzer::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t analyzer::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}




#endif // #ifdef analyzer_cxx
