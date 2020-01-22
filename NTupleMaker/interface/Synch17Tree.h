/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza
//////////////////////////////////////////////////////////

#ifndef Synch17Tree_h
#define Synch17Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//We´ll store for now 5 CP mixing scenarios: "sm_htt125", "ps_htt125", "mm_htt125" "minusmm_htt125", "mix0p375_htt125". Names chose to stay somewhat consistent with choices IC in the CH branch

#include <vector>
#include <string>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Synch17Tree {
public :
  
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain
  
  // Declaration of leaf types
  //Event ID
  UInt_t          run; 
  UInt_t          lumi;
  ULong64_t       evt;
  //Pile up
  Int_t           npv;
  Float_t         npu;
  Float_t         rho;
  
  // MET filters
  Bool_t          passedAllMetFilters;
  
  //Leptons
  Float_t         pt_1;
  Float_t         phi_1;
  Float_t         eta_1;
  Float_t         chconst_1_pt;
  Float_t         chconst_1_eta;
  Float_t         chconst_1_phi;   
  Float_t         m_1;
  Int_t           gen_match_1; 
  Int_t           q_1;
  Float_t         iso_1;
  Float_t         mva_1;
  Float_t         mva17_1;
  Float_t         d0_1;
  Float_t         dZ_1;
  Float_t         d0err_1;
  Float_t         dZerr_1;
  Float_t         ip0x_1;
  Float_t         ip0y_1;
  Float_t         ip0z_1;
  Float_t         ipx_1;
  Float_t         ipy_1;
  Float_t         ipz_1;
  Float_t         ipx_uncorr_1;
  Float_t         ipy_uncorr_1;
  Float_t         ipz_uncorr_1;
  Float_t         ipxy_1;
  Double_t         IP_signif_PV_with_BS_1;
  Float_t         ipn_1;
  Float_t         drip_1;
  Float_t         detaip_1;
  Float_t         dphiip_1;
  Float_t         ipxy_uncorr_1;
  Float_t         ipn_uncorr_1;
  Float_t         drip_uncorr_1;
  Float_t         detaip_uncorr_1;
  Float_t         dphiip_uncorr_1;
  Float_t         mt_1;
  Float_t         puppimt_1;
  Int_t 	  tau_decay_mode_1;
  Float_t         dm_1;
  Float_t         dmMVA_1;
  Float_t 	      chpt_1;
  Float_t 	      cheta_1;
  Float_t 	      chphi_1;
  Float_t 	      chm_1;
  Float_t 	      npt_1;
  Float_t 	      neta_1;
  Float_t 	      nphi_1;
  Float_t 	      nm_1;


  Float_t         pt_2;
  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         chconst_2_pt;
  Float_t         chconst_2_eta;
  Float_t         chconst_2_phi;     
  Float_t         m_2;
  Int_t           gen_match_2; 
  Int_t           q_2;
  Float_t         iso_2;
  Float_t         mva_2;
  Float_t         mva17_2;
  Float_t         d0_2;
  Float_t         dZ_2;
  Float_t         d0err_2;
  Float_t         dZerr_2;
  Float_t         ip0x_2;
  Float_t         ip0y_2;
  Float_t         ip0z_2;
  Float_t         ipx_2;
  Float_t         ipy_2;
  Float_t         ipz_2;
  Float_t         ipx_uncorr_2;
  Float_t         ipy_uncorr_2;
  Float_t         ipz_uncorr_2;
  Float_t         ipxy_2;
  Double_t        IP_signif_PV_with_BS_2;
  Float_t         ipn_2;
  Float_t         drip_2;
  Float_t         detaip_2;
  Float_t         dphiip_2;
  Float_t         ipxy_uncorr_2;
  Float_t         ipn_uncorr_2;
  Float_t         drip_uncorr_2;
  Float_t         detaip_uncorr_2;
  Float_t         dphiip_uncorr_2;
  Float_t         mt_2;
  Float_t         puppimt_2;
  Int_t 	  tau_decay_mode_2;
  Float_t         dm_2;
  Float_t         dmMVA_2;
  Float_t 	      chpt_2;
  Float_t 	      cheta_2;
  Float_t 	      chphi_2;
  Float_t 	      chm_2;
  Float_t 	      npt_2;
  Float_t 	      neta_2;
  Float_t 	      nphi_2;
  Float_t 	      nm_2;

  //TO FIX
  Float_t         againstElectronLooseMVA6_1;
  Float_t         againstElectronMediumMVA6_1;
  Float_t         againstElectronTightMVA6_1;
  Float_t         againstElectronVLooseMVA6_1;
  Float_t         againstElectronVTightMVA6_1;
  Float_t         againstMuonLoose3_1;
  Float_t         againstMuonTight3_1;
  Float_t         chargedIsoPtSum_1;
  Float_t         neutralIsoPtSum_1;
  Float_t         decayModeFindingOldDMs_1;
  Float_t         puCorrPtSum_1;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t         byIsolationMVA3newDMwoLTraw_1;
  Float_t         byIsolationMVA3oldDMwoLTraw_1;
  Float_t         byIsolationMVA3newDMwLTraw_1;
  Float_t         byIsolationMVA3oldDMwLTraw_1;
  Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t 	  byIsolationMVArun2017v2DBoldDMwLTraw2017_1;
  Float_t         byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t         byLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t         byMediumIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t 	  byTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t 	  byVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t 	  byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t         idisoweight_1;
  Float_t         idisoweight_antiiso_1;
  Float_t         trigweight_1;
  Float_t         trigweight_antiiso_1;

  Float_t         againstElectronLooseMVA6_2;
  Float_t         againstElectronMediumMVA6_2;
  Float_t         againstElectronTightMVA6_2;
  Float_t         againstElectronVLooseMVA6_2;
  Float_t         againstElectronVTightMVA6_2;
  Float_t         againstMuonLoose3_2;
  Float_t         againstMuonTight3_2;
  Float_t         chargedIsoPtSum_2;
  Float_t         decayModeFindingOldDMs_2;
  Float_t         neutralIsoPtSum_2;
  Float_t         puCorrPtSum_2;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t         byIsolationMVA3newDMwoLTraw_2;
  Float_t         byIsolationMVA3oldDMwoLTraw_2;
  Float_t         byIsolationMVA3newDMwLTraw_2;
  Float_t         byIsolationMVA3oldDMwLTraw_2;
  Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t 	  byIsolationMVArun2017v2DBoldDMwLTraw2017_2;
  Float_t         byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t         byLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t         byMediumIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t 	  byTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t 	  byVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t 	  byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
    
  Float_t	   deepTauVsEleRaw_1;	
  Float_t	   deepTauVsJetRaw_1;	
  Float_t	   deepTauVsMuRaw_1;	
  Float_t	   deepTauVsEleRaw_2;	
  Float_t	   deepTauVsJetRaw_2;	
  Float_t	   deepTauVsMuRaw_2;	
  Float_t	   byLooseDeepTau2017v2p1VSe_2;	
  Float_t	   byLooseDeepTau2017v2p1VSjet_2;	
  Float_t	   byLooseDeepTau2017v2p1VSmu_2;	
  Float_t	   byMediumDeepTau2017v2p1VSe_2;	
  Float_t	   byMediumDeepTau2017v2p1VSjet_2;	
  Float_t	   byMediumDeepTau2017v2p1VSmu_2;	
  Float_t	   byTightDeepTau2017v2p1VSe_2;	
  Float_t	   byTightDeepTau2017v2p1VSjet_2;	
  Float_t	   byTightDeepTau2017v2p1VSmu_2;	
  Float_t	   byVLooseDeepTau2017v2p1VSe_2;	
  Float_t	   byVLooseDeepTau2017v2p1VSjet_2;	
  Float_t	   byVLooseDeepTau2017v2p1VSmu_2;	
  Float_t	   byVTightDeepTau2017v2p1VSe_2;	
  Float_t	   byVTightDeepTau2017v2p1VSjet_2;	
  Float_t	   byVVLooseDeepTau2017v2p1VSe_2;	
  Float_t	   byVVLooseDeepTau2017v2p1VSjet_2;	
  Float_t	   byVVTightDeepTau2017v2p1VSe_2;	
  Float_t	   byVVTightDeepTau2017v2p1VSjet_2;	
  Float_t	   byVVVLooseDeepTau2017v2p1VSe_2;	
  Float_t	   byVVVLooseDeepTau2017v2p1VSjet_2;	

  Float_t    MVADM2017v1DM0raw_2;
  Float_t    MVADM2017v1DM10raw_2;
  Float_t    MVADM2017v1DM11raw_2;
  Float_t    MVADM2017v1DM1raw_2;
  Float_t    MVADM2017v1DM2raw_2;
  Float_t    MVADM2017v1DMotherraw_2;

  // new 
  Float_t         idisoweight_2;
  Float_t         idisoweight_antiiso_2;
  Float_t         trigweight_2;
  Float_t         trigweight_antiiso_2;
  Float_t         tauvsjetweightMedium_2;

  /////////////////////////////////////////////////////////////// NEW NEW
  Float_t	 efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
//  Float_t	 efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
//  Float_t	 efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1;
  Float_t	 efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1;
//  Float_t	 efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
//  Float_t	 efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  ////////////////////////////////////////////////////////////
  Float_t	 efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
//  Float_t	 efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
//  Float_t	 efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2;
  Float_t	 efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2;
//  Float_t	 efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
//  Float_t	 efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  ////////////////////////////////////////////////////////////
  Float_t	 correction_againstElectronVLooseMVA6_1;
//  Float_t	 correction_againstElectronLooseMVA6_1;
//  Float_t	 correction_againstElectronMediumMVA6_1;
  Float_t	 correction_againstElectronTightMVA6_1;
//  Float_t	 correction_againstElectronVTightMVA6_1;
  ////////////////////////////////////////////////////////////
  Float_t	 correction_againstElectronVLooseMVA6_2;
//  Float_t	 correction_againstElectronLooseMVA6_2;
//  Float_t	 correction_againstElectronMediumMVA6_2;
  Float_t	 correction_againstElectronTightMVA6_2;
//  Float_t	 correction_againstElectronVTightMVA6_2;
  ////////////////////////////////////////////////////////////
  Float_t	 correction_againstMuonLoose3_1;
  Float_t	 correction_againstMuonTight3_1;
  ////////////////////////////////////////////////////////////
  Float_t	 correction_againstMuonLoose3_2;
  Float_t	 correction_againstMuonTight3_2;
  ////////////////////////////////////////////////////////////
  //Trig and weights
  Float_t         weight;
  Float_t         mcweight;
  Float_t         puweight;
  Float_t         effweight;
  Float_t         trigweight;
  Float_t         embweight;

  Float_t         topptweight;
  Double_t 	  zptweight;
  Double_t        trkeffweight;
  Float_t         etaufakeweight;
  Float_t         mutaufakeweight;
  
  Bool_t          trg_singlemuon;
  Bool_t	  trg_singleelectron;
  Bool_t          singleLepTrigger;
  Bool_t          trg_mutaucross;
  Bool_t          trg_mutaucross_mu;
  Bool_t          trg_mutaucross_tau;
  Bool_t          trg_doubletau;
  Bool_t          ditauTrigger;
  Bool_t          xTriggerLep;
  Bool_t          xTriggerTau;
  Bool_t          xTrigger;
  //MET
  Float_t         met;
  Float_t         metphi;  
  Float_t         met_uncorr;
  Float_t         metphi_uncorr;  
  Float_t         met_rcmr;
  Float_t         metphi_rcmr;  
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;
  Float_t         pzetavis;
  Float_t         pzetamiss;

  //PUPPI MET
  Float_t         puppimet;
  Float_t         puppimetphi;  
  Float_t         puppimet_rcmr;
  Float_t         puppimetphi_rcmr;  
  Float_t         puppimetcov00;
  Float_t         puppimetcov01;
  Float_t         puppimetcov10;
  Float_t         puppimetcov11;
  Float_t         puppipzetamiss;

  //di tau system
  Float_t         pt_tt;
  Float_t         m_vis;
  Float_t         mt_tot;
  Float_t         m_sv;
  Float_t         pt_sv;
  Float_t         eta_sv;
  Float_t         phi_sv;
  Float_t         met_sv;
  Float_t         mt_sv;
  
  Float_t         m_fast;
  Float_t         pt_fast;
  Float_t         eta_fast;
  Float_t         phi_fast;
  Float_t         mt_fast;

  //VBF
  Float_t         mjj;
  Float_t         jdeta;
  Float_t         jdphi;
  Float_t         dijetpt; 
  Float_t         dijeteta; 
  Float_t         dijetphi; 
  Int_t           njetingap;
  Int_t           njetingap20;

  //Jets
  Int_t           njets;
  Int_t           njetshad;
  Int_t           njetspt20;
  Float_t         jpt_1;
  Float_t         jeta_1;
  Float_t         jphi_1;
  Float_t         jcsv_1;
  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jcsv_2;
  
  //b-jets
  Int_t           nbtag;
  Float_t         bpt_1;
  Float_t         beta_1;
  Float_t         bphi_1;
  Float_t         bcsv_1;
  Float_t         bpt_2;
  Float_t         beta_2;
  Float_t         bphi_2;
  Float_t         bcsv_2;
  
  //Trigger 

  Int_t           gen_noutgoing;
  Int_t           os;
  Int_t           dilepton_veto;
  Int_t           extraelec_veto;
  Int_t           extramuon_veto;

  //CP measurement

  Float_t         tau1DecayPlaneX;
  Float_t         tau1DecayPlaneY;
  Float_t         tau1DecayPlaneZ;
  Float_t         tau2DecayPlaneX;
  Float_t         tau2DecayPlaneY;
  Float_t         tau2DecayPlaneZ;

  Float_t         acotautau_00;
  Float_t         acotautau_10;
  Float_t         acotautau_01;
  Float_t         acotautau_11;

  Float_t         acotautau_bs_00;
  Float_t         acotautau_bs_10;
  Float_t         acotautau_bs_01;
  Float_t         acotautau_bs_11;

  Float_t         acotautau_refitbs_00;
  Float_t         acotautau_refitbs_10;
  Float_t         acotautau_refitbs_01;
  Float_t         acotautau_refitbs_11;
  Float_t         acotautau_refitbs_02;

  Float_t         acotautau_helix_00;
  Float_t         acotautau_helix_10;
  Float_t         acotautau_helix_01;
  Float_t         acotautau_helix_11;

  Float_t         acotautau_refitbs_uncorr_00;
  Float_t         acotautau_refitbs_uncorr_10;
  Float_t         acotautau_refitbs_uncorr_01;
  Float_t         acotautau_refitbs_uncorr_11;
 
  Float_t         acotautau_helix_uncorr_00;
  Float_t         acotautau_helix_uncorr_10;
  Float_t         acotautau_helix_uncorr_01;
  Float_t         acotautau_helix_uncorr_11;

  //Merijn add acotau for psi:
  Float_t         acotautauPsi_00;
  Float_t         acotautauPsi_10;
  Float_t         acotautauPsi_01;
  Float_t         acotautauPsi_11;

  Int_t		  pdgcodetau2;

  //Points of closest approach
  Float_t         tau_pca2D_x_1;
  Float_t         tau_pca2D_y_1;
  Float_t         tau_pca2D_z_1;
  Float_t         tau_pca3D_x_1;
  Float_t         tau_pca3D_y_1;
  Float_t         tau_pca3D_z_1;
  Float_t         tau_pca2D_x_2;
  Float_t         tau_pca2D_y_2;
  Float_t         tau_pca2D_z_2;
  Float_t         tau_pca3D_x_2;
  Float_t         tau_pca3D_y_2;
  Float_t         tau_pca3D_z_2;

  //Secondary vertices
  Float_t         tau_SV_x_1;
  Float_t         tau_SV_y_1;
  Float_t         tau_SV_z_1;
  Float_t         tau_SV_covxx_1;
  Float_t         tau_SV_covyx_1;
  Float_t         tau_SV_covzx_1;
  Float_t         tau_SV_covyy_1;
  Float_t         tau_SV_covzy_1;
  Float_t         tau_SV_covzz_1;
  Float_t         tau_SV_x_2;
  Float_t         tau_SV_y_2;
  Float_t         tau_SV_z_2;
  Float_t         tau_SV_covxx_2;
  Float_t         tau_SV_covyx_2;
  Float_t         tau_SV_covzx_2;
  Float_t         tau_SV_covyy_2;
  Float_t         tau_SV_covzy_2;
  Float_t         tau_SV_covzz_2;

  Float_t alpha_IP_1;
  Float_t alpha_IP_uncorr_1;
  Float_t alpha_plane_1;

  Float_t alpha_IP_2;
  Float_t alpha_IP_uncorr_2;
  Float_t alpha_plane_2;

  //reco vertices
  Float_t RecoVertexX;
  Float_t RecoVertexY;
  Float_t RecoVertexZ;

  Float_t pvx;
  Float_t pvy;
  Float_t pvz;
  Bool_t  is_refitted_PV_with_BS;
  Int_t   v_tracks;
  
  Float_t GenVertexX;
  Float_t GenVertexY;
  Float_t GenVertexZ;

  //Merijn: add the vx of the tau decay products
  Float_t VxConstitTau1;
  Float_t VyConstitTau1;
  Float_t VzConstitTau1;
  
  Float_t VxConstitTau2;
  Float_t VyConstitTau2;
  Float_t VzConstitTau2;
  Float_t alphaminus;
  Float_t alphaminus_uncorr;

  Double_t TauSpinnerWeightsEven;
  Double_t TauSpinnerWeightsOdd;
  Double_t TauSpinnerWeightsMaxMix;
  Double_t TauSpinnerWeightsMinusMaxMix;
  Double_t TauSpinnerWeightsMix0p375;

  //Vinay: ditau_vis_pT + MET
  Float_t Prompt_pT;

  Bool_t isrefitBS;

  Bool_t apply_recoil;

  //////////////////////////////////////////////
  //            List of branches              //
  //////////////////////////////////////////////
  TBranch        *b_run;   //!
  TBranch	 *b_lumi;
  TBranch	 *b_evt;
  //Pile up
  TBranch	 *b_npv;
  TBranch	 *b_npu;
  TBranch	 *b_rho;
  
  // MET filters
  TBranch	 *b_passedAllMetFilters;
  
  //Leptons
  TBranch	 *b_pt_1;
  TBranch	 *b_phi_1;
  TBranch	 *b_eta_1;
  TBranch        *b_chconst_1_pt;
  TBranch        *b_chconst_1_eta;
  TBranch        *b_chconst_1_phi;     
  TBranch	 *b_m_1;
  TBranch	 *b_gen_match_1; 
  TBranch	 *b_q_1;
  TBranch	 *b_iso_1;
  TBranch	 *b_mva_1;
  TBranch	 *b_mva17_1;
  TBranch	 *b_d0_1;
  TBranch	 *b_dZ_1;
  TBranch	 *b_d0err_1;
  TBranch	 *b_dZerr_1;
  TBranch	 *b_ip0x_1;
  TBranch	 *b_ip0y_1;
  TBranch	 *b_ip0z_1;
  TBranch	 *b_ipx_1;
  TBranch	 *b_ipy_1;
  TBranch	 *b_ipz_1;
  TBranch	 *b_ipx_uncorr_1;
  TBranch	 *b_ipy_uncorr_1;
  TBranch	 *b_ipz_uncorr_1;
  TBranch	 *b_ipxy_1;
  TBranch	 *b_IP_signif_PV_with_BS_1;
  TBranch	 *b_ipn_1;
  TBranch	 *b_drip_1;
  TBranch	 *b_detaip_1;
  TBranch	 *b_dphiip_1;
  TBranch	 *b_ipxy_uncorr_1;
  TBranch	 *b_ipn_uncorr_1;
  TBranch	 *b_drip_uncorr_1;
  TBranch	 *b_detaip_uncorr_1;
  TBranch	 *b_dphiip_uncorr_1;
  TBranch	 *b_mt_1;
  TBranch	 *b_puppimt_1;
  TBranch  *b_tau_decay_mode_1;
  TBranch  *b_dm_1;
  TBranch  *b_dmMVA_1;
  TBranch  *b_chpt_1;
  TBranch  *b_cheta_1;
  TBranch  *b_chphi_1;
  TBranch  *b_chm_1;
  TBranch  *b_npt_1;
  TBranch  *b_neta_1;
  TBranch  *b_nphi_1;
  TBranch  *b_nm_1;

  TBranch	 *b_pt_2;
  TBranch	 *b_phi_2;
  TBranch	 *b_eta_2;
  TBranch  *b_chconst_2_pt;
  TBranch  *b_chconst_2_eta;
  TBranch  *b_chconst_2_phi; 
  TBranch	 *b_m_2;
  TBranch	 *b_gen_match_2; 
  TBranch	 *b_q_2;
  TBranch	 *b_iso_2;
  TBranch	 *b_mva_2;
  TBranch	 *b_mva17_2;
  TBranch	 *b_d0_2;
  TBranch	 *b_dZ_2;
  TBranch	 *b_d0err_2;
  TBranch	 *b_dZerr_2;
  TBranch	 *b_ip0x_2;
  TBranch	 *b_ip0y_2;
  TBranch	 *b_ip0z_2;
  TBranch	 *b_ipx_2;
  TBranch	 *b_ipy_2;
  TBranch	 *b_ipz_2;
  TBranch	 *b_ipx_uncorr_2;
  TBranch	 *b_ipy_uncorr_2;
  TBranch	 *b_ipz_uncorr_2;
  TBranch	 *b_ipxy_2;
  TBranch	 *b_IP_signif_PV_with_BS_2;
  TBranch	 *b_ipn_2;
  TBranch	 *b_drip_2;
  TBranch	 *b_detaip_2;
  TBranch	 *b_dphiip_2;
  TBranch	 *b_ipxy_uncorr_2;
  TBranch	 *b_ipn_uncorr_2;
  TBranch	 *b_drip_uncorr_2;
  TBranch	 *b_detaip_uncorr_2;
  TBranch	 *b_dphiip_uncorr_2;
  TBranch	 *b_mt_2;
  TBranch	 *b_puppimt_2;
  TBranch  *b_tau_decay_mode_2;
  TBranch  *b_dm_2;
  TBranch  *b_dmMVA_2;
  TBranch  *b_chpt_2;
  TBranch  *b_cheta_2;
  TBranch  *b_chphi_2;
  TBranch  *b_chm_2;
  TBranch  *b_npt_2;
  TBranch  *b_neta_2;
  TBranch  *b_nphi_2;
  TBranch  *b_nm_2;

  //TO FIX
  TBranch	 *b_againstElectronLooseMVA6_1;
  TBranch	 *b_againstElectronMediumMVA6_1;
  TBranch	 *b_againstElectronTightMVA6_1;
  TBranch	 *b_againstElectronVLooseMVA6_1;
  TBranch	 *b_againstElectronVTightMVA6_1;
  TBranch	 *b_againstMuonLoose3_1;
  TBranch	 *b_againstMuonTight3_1;
  TBranch	 *b_chargedIsoPtSum_1;
  TBranch	 *b_decayModeFindingOldDMs_1;
  TBranch	 *b_neutralIsoPtSum_1;
  TBranch	 *b_puCorrPtSum_1;
  TBranch	 *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  TBranch	 *b_byIsolationMVA3newDMwoLTraw_1;
  TBranch	 *b_byIsolationMVA3oldDMwoLTraw_1;
  TBranch	 *b_byIsolationMVA3newDMwLTraw_1;
  TBranch	 *b_byIsolationMVA3oldDMwLTraw_1;
  TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch 	 *b_byIsolationMVArun2017v2DBoldDMwLTraw2017_1;
  TBranch	 *b_byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_byTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1;

  TBranch	   *b_deepTauVsEleRaw_1;	
  TBranch	   *b_deepTauVsJetRaw_1;	
  TBranch	   *b_deepTauVsMuRaw_1;	
  TBranch	   *b_deepTauVsEleRaw_2;	
  TBranch	   *b_deepTauVsJetRaw_2;	
  TBranch	   *b_deepTauVsMuRaw_2;	
  TBranch	   *b_byLooseDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byLooseDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byLooseDeepTau2017v2p1VSmu_2;	
  TBranch	   *b_byMediumDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byMediumDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byMediumDeepTau2017v2p1VSmu_2;	
  TBranch	   *b_byTightDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byTightDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byTightDeepTau2017v2p1VSmu_2;	
  TBranch	   *b_byVLooseDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byVLooseDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byVLooseDeepTau2017v2p1VSmu_2;	
  TBranch	   *b_byVTightDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byVTightDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byVVLooseDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byVVLooseDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byVVTightDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byVVTightDeepTau2017v2p1VSjet_2;	
  TBranch	   *b_byVVVLooseDeepTau2017v2p1VSe_2;	
  TBranch	   *b_byVVVLooseDeepTau2017v2p1VSjet_2;	

  TBranch    *b_MVADM2017v1DM0raw_2;
  TBranch    *b_MVADM2017v1DM10raw_2;
  TBranch    *b_MVADM2017v1DM11raw_2;
  TBranch    *b_MVADM2017v1DM1raw_2;
  TBranch    *b_MVADM2017v1DM2raw_2;
  TBranch    *b_MVADM2017v1DMotherraw_2;

  ///////////////////////////////////////////////////////////NEW
  TBranch	 *b_efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
//  TBranch	 *b_efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
//  TBranch	 *b_efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1;
  TBranch	 *b_efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1;
//  TBranch	 *b_efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
//  TBranch	 *b_efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
//  TBranch	 *b_efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
//  TBranch	 *b_efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2;
//  TBranch	 *b_efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
//  TBranch	 *b_efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_correction_againstElectronVLooseMVA6_1;
//  TBranch	 *b_correction_againstElectronLooseMVA6_1;
//  TBranch	 *b_correction_againstElectronMediumMVA6_1;
  TBranch	 *b_correction_againstElectronTightMVA6_1;
//  TBranch	 *b_correction_againstElectronVTightMVA6_1;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_correction_againstElectronVLooseMVA6_2;
//  TBranch	 *b_correction_againstElectronLooseMVA6_2;
//  TBranch	 *b_correction_againstElectronMediumMVA6_2;
  TBranch	 *b_correction_againstElectronTightMVA6_2;
//  TBranch	 *b_correction_againstElectronVTightMVA6_2;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_correction_againstMuonLoose3_1;
  TBranch	 *b_correction_againstMuonTight3_1;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_correction_againstMuonLoose3_2;
  TBranch	 *b_correction_againstMuonTight3_2;
  ////////////////////////////////////////////////////////////
  TBranch	 *b_idisoweight_1;
  TBranch	 *b_idisoweight_antiiso_1;
  TBranch	 *b_trigweight_1;
  TBranch	 *b_trigweight_antiiso_1;
  TBranch	 *b_againstElectronLooseMVA6_2;
  TBranch	 *b_againstElectronMediumMVA6_2;
  TBranch	 *b_againstElectronTightMVA6_2;
  TBranch	 *b_againstElectronVLooseMVA6_2;
  TBranch	 *b_againstElectronVTightMVA6_2;
  TBranch	 *b_againstMuonLoose3_2;
  TBranch	 *b_againstMuonTight3_2;
  TBranch	 *b_chargedIsoPtSum_2;
  TBranch	 *b_decayModeFindingOldDMs_2;
  TBranch	 *b_neutralIsoPtSum_2;
  TBranch	 *b_puCorrPtSum_2;
  TBranch	 *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  TBranch	 *b_byIsolationMVA3newDMwoLTraw_2;
  TBranch	 *b_byIsolationMVA3oldDMwoLTraw_2;
  TBranch	 *b_byIsolationMVA3newDMwLTraw_2;
  TBranch	 *b_byIsolationMVA3oldDMwLTraw_2;
  TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch 	 *b_byIsolationMVArun2017v2DBoldDMwLTraw2017_2;
  TBranch	 *b_byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_byTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  TBranch	 *b_idisoweight_2;
  TBranch	 *b_idisoweight_antiiso_2;
  TBranch	 *b_trigweight_2;
  TBranch	 *b_trigweight_antiiso_2;
  TBranch	 *b_tauvsjetweightMedium_2;

  //Trig and weights
  TBranch	 *b_weight;
  TBranch	 *b_mcweight;
  TBranch	 *b_puweight;
  TBranch	 *b_effweight;
  TBranch        *b_trigweight;
  TBranch        *b_embweight;

  TBranch	 *b_topptweight;
  TBranch	 *b_zptweight;
  TBranch	 *b_trkeffweight;
  TBranch	 *b_etaufakeweight;
  TBranch	 *b_mutaufakeweight;

  TBranch	 *b_trg_singlemuon;
  TBranch	 *b_trg_singleelectron;
  TBranch	 *b_singleLepTrigger;
  TBranch	 *b_trg_mutaucross;
  TBranch  *b_trg_mutaucross_mu;
  TBranch  *b_trg_mutaucross_tau;
  TBranch  *b_trg_doubletau;
  TBranch  *b_ditauTrigger;
  TBranch	 *b_xTrigger;
  TBranch  *b_xTriggerLep;
  TBranch  *b_xTriggerTau;
  //MET
  TBranch	 *b_met;
  TBranch	 *b_metphi;
  TBranch	 *b_met_uncorr;
  TBranch	 *b_metphi_uncorr;
  TBranch	 *b_met_rcmr;
  TBranch	 *b_metphi_rcmr;  
  TBranch	 *b_metcov00;
  TBranch	 *b_metcov01;
  TBranch	 *b_metcov10;
  TBranch	 *b_metcov11;
  TBranch	 *b_pzetavis;
  TBranch	 *b_pzetamiss;

  //PUPPI MET
  TBranch	 *b_puppimet;
  TBranch	 *b_puppimetphi;
  TBranch	 *b_puppimet_rcmr;
  TBranch	 *b_puppimetphi_rcmr;  
  TBranch	 *b_puppimetcov00;
  TBranch	 *b_puppimetcov01;
  TBranch	 *b_puppimetcov10;
  TBranch	 *b_puppimetcov11;
  TBranch	 *b_puppipzetamiss;

  //di tau system
  TBranch	 *b_pt_tt;
  TBranch	 *b_m_vis;
  TBranch	 *b_mt_tot;
  TBranch	 *b_m_sv;
  TBranch	 *b_pt_sv;
  TBranch	 *b_eta_sv;
  TBranch	 *b_phi_sv;
  TBranch	 *b_met_sv;
  TBranch	 *b_mt_sv;

  TBranch	 *b_m_fast;
  TBranch	 *b_pt_fast;
  TBranch	 *b_eta_fast;
  TBranch	 *b_phi_fast;
  TBranch	 *b_mt_fast;

  //VBF
  TBranch	 *b_mjj;
  TBranch	 *b_jdeta;
  TBranch	 *b_dijetpt;
  TBranch	 *b_dijeteta;
  TBranch	 *b_dijetphi;
  TBranch	 *b_jdphi;
  TBranch	 *b_njetingap;
  TBranch	 *b_njetingap20;

  //Jets
  TBranch	 *b_njets;
  TBranch	 *b_njetshad;
  TBranch	 *b_njetspt20;
  TBranch	 *b_jpt_1;
  TBranch	 *b_jeta_1;
  TBranch	 *b_jphi_1;
  TBranch	 *b_jcsv_1;
  TBranch	 *b_jpt_2;
  TBranch	 *b_jeta_2;
  TBranch	 *b_jphi_2;
  TBranch	 *b_jcsv_2;
  
  //b-jets
  TBranch	 *b_nbtag;
  TBranch	 *b_bpt_1;
  TBranch	 *b_beta_1;
  TBranch	 *b_bphi_1;
  TBranch	 *b_bcsv_1;
  TBranch	 *b_bpt_2;
  TBranch	 *b_beta_2;
  TBranch	 *b_bphi_2;
  TBranch	 *b_bcsv_2;
  
  //Misc
  

  TBranch	 *b_gen_noutgoing;
  TBranch	 *b_os;
  TBranch	 *b_dilepton_veto;
  TBranch	 *b_extraelec_veto;
  TBranch	 *b_extramuon_veto;

  //CP measurement

  TBranch        *b_tau1DecayPlaneX;
  TBranch        *b_tau1DecayPlaneY;
  TBranch        *b_tau1DecayPlaneZ;
  TBranch        *b_tau2DecayPlaneX;
  TBranch        *b_tau2DecayPlaneY;
  TBranch        *b_tau2DecayPlaneZ;

  TBranch        *b_acotautau_00;
  TBranch        *b_acotautau_10;
  TBranch        *b_acotautau_01;
  TBranch        *b_acotautau_11;

  TBranch        *b_acotautau_bs_00;
  TBranch        *b_acotautau_bs_10;
  TBranch        *b_acotautau_bs_01;
  TBranch        *b_acotautau_bs_11;

  TBranch        *b_acotautau_refitbs_00;
  TBranch        *b_acotautau_refitbs_10;
  TBranch        *b_acotautau_refitbs_01;
  TBranch        *b_acotautau_refitbs_11;
  TBranch        *b_acotautau_refitbs_02;

  TBranch        *b_acotautau_helix_00;
  TBranch        *b_acotautau_helix_10;
  TBranch        *b_acotautau_helix_01;
  TBranch        *b_acotautau_helix_11;

  TBranch        *b_acotautau_refitbs_uncorr_00;
  TBranch        *b_acotautau_refitbs_uncorr_10;
  TBranch        *b_acotautau_refitbs_uncorr_01;
  TBranch        *b_acotautau_refitbs_uncorr_11;

  TBranch        *b_acotautau_helix_uncorr_00;
  TBranch        *b_acotautau_helix_uncorr_10;
  TBranch        *b_acotautau_helix_uncorr_01;
  TBranch        *b_acotautau_helix_uncorr_11;

  TBranch        *b_acotautauPsi_00;
  TBranch        *b_acotautauPsi_10;
  TBranch        *b_acotautauPsi_01;
  TBranch        *b_acotautauPsi_11;
  TBranch	 *b_pdgcodetau2;

  //Points of closest approach
  TBranch        *b_tau_pca2D_x_1;
  TBranch        *b_tau_pca2D_y_1;
  TBranch        *b_tau_pca2D_z_1;
  TBranch        *b_tau_pca3D_x_1;
  TBranch        *b_tau_pca3D_y_1;
  TBranch        *b_tau_pca3D_z_1;
  TBranch        *b_tau_pca2D_x_2;
  TBranch        *b_tau_pca2D_y_2;
  TBranch        *b_tau_pca2D_z_2;
  TBranch        *b_tau_pca3D_x_2;
  TBranch        *b_tau_pca3D_y_2;
  TBranch        *b_tau_pca3D_z_2;

  //Secondary vertices
  TBranch        *b_tau_SV_x_1;
  TBranch        *b_tau_SV_y_1;
  TBranch        *b_tau_SV_z_1;
  TBranch        *b_tau_SV_covxx_1;
  TBranch        *b_tau_SV_covyx_1;
  TBranch        *b_tau_SV_covzx_1;
  TBranch        *b_tau_SV_covyy_1;
  TBranch        *b_tau_SV_covzy_1;
  TBranch        *b_tau_SV_covzz_1;

  TBranch        *b_tau_SV_x_2;
  TBranch        *b_tau_SV_y_2;
  TBranch        *b_tau_SV_z_2;
  TBranch        *b_tau_SV_covxx_2;
  TBranch        *b_tau_SV_covyx_2;
  TBranch        *b_tau_SV_covzx_2;
  TBranch        *b_tau_SV_covyy_2;
  TBranch        *b_tau_SV_covzy_2;
  TBranch        *b_tau_SV_covzz_2;

  //reco vertices
  //RECO vertex info is practical to have

  TBranch        *b_RecoVertexX;
  TBranch        *b_RecoVertexY;
  TBranch        *b_RecoVertexZ;  

  TBranch        *b_pvx;
  TBranch        *b_pvy;
  TBranch        *b_pvz;  
  TBranch        *b_is_refitted_PV_with_BS;  
  TBranch        *b_v_tracks;

//gen vertex info is practical to have
  TBranch        *b_GenVertexX;
  TBranch        *b_GenVertexY;
  TBranch        *b_GenVertexZ;
  
  TBranch        *b_VxConstitTau1;
  TBranch        *b_VyConstitTau1;
  TBranch        *b_VzConstitTau1;

  TBranch        *b_VxConstitTau2;
  TBranch        *b_VyConstitTau2;
  TBranch        *b_VzConstitTau2;
  TBranch        *b_alphaminus;
  TBranch        *b_alphaminus_uncorr;

  TBranch        *b_alpha_IP_1;
  TBranch        *b_alpha_IP_uncorr_1;
  TBranch        *b_alpha_plane_1;

  TBranch        *b_alpha_IP_2;
  TBranch        *b_alpha_IP_uncorr_2;
  TBranch        *b_alpha_plane_2;

  TBranch        *b_TauSpinnerWeightsEven;
  TBranch        *b_TauSpinnerWeightsOdd;
  TBranch        *b_TauSpinnerWeightsMaxMix;
  TBranch        *b_TauSpinnerWeightsMinusMaxMix;
  TBranch        *b_TauSpinnerWeightsMix0p375;

  TBranch       *b_Prompt_pT;

  TBranch       *b_isrefitBS;

  Synch17Tree(TTree *tree=0);
  virtual ~Synch17Tree();

  virtual void Init(TTree *tree);
  
  //Read methods
  virtual void     ReadInit(TTree *tree);
  virtual void     ReadReset();
  virtual Long64_t GetEntries();
  virtual Long64_t LoadTree(Long64_t entry);
  virtual Long64_t GetEntry(Long64_t entry);
  virtual Long64_t LoadedEntryId();
  virtual void     Show(Long64_t entry = -1);
  virtual Int_t    Cut(Long64_t entry);

  //Write methods
  virtual void WriteInit(TTree *tree);
  virtual void Fill();

protected:
  bool lock;
  Long64_t ientry;
};

#endif //!endif Synch17Tree_h
