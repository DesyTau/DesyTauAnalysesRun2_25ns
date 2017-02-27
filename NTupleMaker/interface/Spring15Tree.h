/////////////////////////////////////////////////////////////
// Read and Write Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
// 
// Based on a class automatically generated on
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
// from TTree TauCheck/TauCheck
// found on file: VBFHToTauTau_M-125_Synch.root
//////////////////////////////////////////////////////////

#ifndef Spring15Tree_h
#define Spring15Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Spring15Tree {
public :
  
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain
  
  // Declaration of leaf types
  UInt_t          run;
  UInt_t          lumi;
  ULong64_t       evt;
  Int_t           npv;
  Float_t         npu;
  Float_t         rho;
  Float_t         xs;
  Float_t         mcweight;
  Float_t         pu_weight;
  Float_t         trigweight_1;
  Float_t         trigweight_antiiso_1;
  Float_t         trigweight_2;
  Float_t         idisoweight_1;
  Float_t         idisoweight_antiiso_1;
  Float_t         idisoweight_2;
  Float_t         topptweight;
  Double_t 	  zptweight;
  Double_t        trkeffweight_1;
  Float_t         effweight;
  Float_t         etaufakeweight;
  Float_t         mutaufakeweight;
  Float_t         fakeweight;
  Float_t         embeddedWeight;
  Float_t         signalWeight;
  Float_t         weight;
  Float_t         lheHt;
  Int_t           gen_noutgoing;
  Int_t           njetshad;
  Float_t         genV_px;
  Float_t         genV_py; 
  Float_t         genV_pz;
  Float_t         genV_e;
  Float_t         genL_px;
  Float_t         genL_py; 
  Float_t         genL_pz;
  Float_t         genL_e;
  Float_t         m_vis;
  Float_t         mt_tot;
  Float_t         m_sv;
  Float_t         pt_sv;
  Float_t         eta_sv;
  Float_t         phi_sv;
  Float_t         met_sv;
  Float_t         mt_sv;
  Float_t         pt_1;
  Float_t         phi_1;
  Float_t         eta_1;
  Float_t         m_1;
  Int_t           gen_match_1; 
  Int_t           q_1;
  Float_t         iso_1;
  Float_t         mva_1;
  Float_t         d0_1;
  Float_t         dZ_1;
  Float_t         mt_1;
  Float_t         pfmt_1;
  Float_t         puppimt_1;
  Float_t         mt_rcqr_1;
  Float_t         pfmt_rcqr_1;
  Float_t         puppimt_rcqr_1;
  Float_t         mt_rcmr_1;
  Float_t         pfmt_rcmr_1;
  Float_t         puppimt_rcmr_1;
  Float_t         mt_rc_njetsreco_1;
  Float_t         pfmt_rc_njetsreco_1;
  Float_t         mt_rc_visreco_1;
  Float_t         pfmt_rc_visreco_1;
  Int_t 	  tau_decay_mode_1;
  Float_t         pt_2;
  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         m_2;
  Int_t           gen_match_2; 
  Int_t           q_2;
  Float_t         iso_2;
  Float_t         d0_2;
  Float_t         dZ_2;
  Float_t         mva_2;
  Float_t         mt_2;
  Float_t         pfmt_2;
  Float_t         puppimt_2;
  Float_t         mt_rcqr_2;
  Float_t         pfmt_rcqr_2;
  Float_t         puppimt_rcqr_2;
  Float_t         mt_rcmr_2;
  Float_t         pfmt_rcmr_2;
  Float_t         mt_rc_njetsreco_2;
  Float_t         pfmt_rc_njetsreco_2;
  Float_t         mt_rc_visreco_2;
  Float_t         pfmt_rc_visreco_2;  
  Float_t         puppimt_rcmr_2; 
  Int_t 	      tau_decay_mode_2;
  Int_t           os;
  Int_t           dilepton_veto;
  Int_t           extraelec_veto;
  Int_t           extramuon_veto;
  Float_t 		  chargedIsoPtSum_2;
  Float_t         neutralIsoPtSum_2;
  Float_t         puCorrPtSum_2;
  UInt_t          isolationGammaCands_size_2;
  UInt_t          signalGammaCands_size_2;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_1;
  Float_t         againstElectronLooseMVA5_1;
  Float_t         againstElectronMediumMVA5_1;
  Float_t         againstElectronTightMVA5_1;
  Float_t         againstElectronVLooseMVA5_1;
  Float_t         againstElectronVTightMVA5_1;
  Float_t 	  againstElectronVLooseMVA6_2;
  Float_t 	  againstElectronTightMVA6_2;
  Float_t         againstMuonLoose3_1;
  Float_t         againstMuonTight3_1;
  UInt_t          n_badmuons;
  UInt_t          n_duplicatemuons;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t 		  byIsolationMVArun2v1DBoldDMwLTraw_2;
  Float_t         byVLooseIsolationMVArun2v1DBoldDMwLT_2;
  Float_t         byLooseIsolationMVArun2v1DBoldDMwLT_2;
  Float_t         byMediumIsolationMVArun2v1DBoldDMwLT_2;
  Float_t 	      byTightIsolationMVArun2v1DBoldDMwLT_2;
  Float_t 	      byVTightIsolationMVArun2v1DBoldDMwLT_2;
  Float_t 	      byVVTightIsolationMVArun2v1DBoldDMwLT_2;
  Float_t         againstElectronLooseMVA5_2;
  Float_t         againstElectronMediumMVA5_2;
  Float_t         againstElectronTightMVA5_2;
  Float_t         againstElectronVLooseMVA5_2;
  Float_t         againstElectronVTightMVA5_2;
  Float_t         againstMuonLoose3_2;
  Float_t         againstMuonTight3_2;
  Float_t         met;
  Float_t         metphi;  
  Float_t         metcov00;
  Float_t         metcov01;
  Float_t         metcov10;
  Float_t         metcov11;
  Float_t         met_rcqr;
  Float_t         metphi_rcqr;
  Float_t         met_rcmr;
  Float_t         metphi_rcmr;  
  Float_t         met_rc_njetsreco;
  Float_t         metphi_rc_njetsreco;
  Float_t         met_rc_visreco;
  Float_t         metphi_rc_visreco;  
  Float_t         mvamet;
  Float_t         mvametphi;
  Float_t         mvacov00;
  Float_t         mvacov01;
  Float_t         mvacov10;
  Float_t         mvacov11;
  Float_t         mvamet_rcqr;
  Float_t         mvametphi_rcqr;
  Float_t         mvamet_rcmr;
  Float_t         mvametphi_rcmr;
  Float_t         mvamet_rc_njetsreco;
  Float_t         mvametphi_rc_njetsreco;
  Float_t         mvamet_rc_visreco;
  Float_t         mvametphi_rc_visreco;  
  Float_t         puppimet;
  Float_t         puppimetphi;
  Float_t         puppimet_rcqr;
  Float_t         puppimetphi_rcqr;
  Float_t         puppimet_rcmr;
  Float_t         puppimetphi_rcmr;
  Float_t         pt_tt;
  Float_t         pzetavis;
  Float_t         pzetamiss;
  Float_t         pfpzetamiss;
  Float_t         puppipzetamiss;
  Float_t         pzetamiss_rcqr;
  Float_t         pfpzetamiss_rcqr;
  Float_t         puppipzetamiss_rcqr;
  Float_t         pzetamiss_rcmr;
  Float_t         pfpzetamiss_rcmr;
  Float_t         puppipzetamiss_rcmr;
  Float_t         pzetamiss_rc_njetsreco;
  Float_t         pfpzetamiss_rc_njetsreco;
  Float_t         pzetamiss_rc_visreco;
  Float_t         pfpzetamiss_rc_visreco;  
  Float_t         mva_gf;
  Int_t           njets;
  Int_t           njetspt20;
  Float_t         jpt_1;
  Float_t         jeta_1;
  Float_t         jphi_1;
  Float_t         jptraw_1;
  Float_t         jrawf_1;
  Float_t         jptunc_1;
  Float_t         jmva_1;
  Float_t         jlrm_1;
  Int_t           jctm_1;
  Bool_t          jpuid_loose_1;
  Bool_t          jpuid_medium_1;
  Bool_t          jpuid_tight_1;
  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jptraw_2;
  Float_t         jrawf_2;
  Float_t         jptunc_2;
  Float_t         jmva_2;
  Float_t         jlrm_2;
  Int_t           jctm_2;
  Bool_t          jpuid_loose_2;
  Bool_t          jpuid_medium_2;
  Bool_t          jpuid_tight_2;
  Float_t         mjj;
  Float_t         jdeta;
  Int_t           njetingap;
  Int_t           nbtag;
  Float_t         bpt_1;
  Float_t         beta_1;
  Float_t         bphi_1;
  Float_t         brawf_1;
  Float_t         bmva_1;
  Float_t         bcsv_1;
  Bool_t          bpuid_loose_1;
  Bool_t          bpuid_medium_1;
  Bool_t          bpuid_tight_1;
  Float_t         bpt_2;
  Float_t         beta_2;
  Float_t         bphi_2;
  Float_t         brawf_2;
  Float_t         bmva_2;
  Float_t         bcsv_2;
  Bool_t          bpuid_loose_2;
  Bool_t          bpuid_medium_2;
  Bool_t          bpuid_tight_2;
  Bool_t          singleLepTrigger;
  Bool_t          xTrigger;
  
  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_evt;   //!
  TBranch        *b_npv;   //!
  TBranch        *b_npu;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_xs;   //!
  TBranch        *b_mcweight;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_trigweight_1;   //!
  TBranch        *b_trigweight_antiiso_1;   //!
  TBranch        *b_trigweight_2;   //!
  TBranch        *b_idisoweight_1;   //!
  TBranch        *b_idisoweight_antiiso_1;   //!
  TBranch        *b_idisoweight_2;   //!
  TBranch        *b_topptweight;   //! 
  TBranch        *b_zptweight;
  TBranch        *b_trkeffweight_1;
  TBranch        *b_effweight;   //! 
  TBranch        *b_etaufakeweight;   //!
  TBranch        *b_mutaufakeweight;   //!
  TBranch        *b_fakeweight;   //!
  TBranch        *b_embeddedWeight;   //!
  TBranch        *b_signalWeight;   //!
  TBranch        *b_weight;   //!
  TBranch        *b_lheHt;   //!
  TBranch        *b_gen_noutgoing;   //!  
  TBranch        *b_njetshad;   //!  
  TBranch        *b_genV_px;    //!
  TBranch        *b_genV_py;    //!  
  TBranch        *b_genV_pz;    //!
  TBranch        *b_genV_e;    //!
  TBranch        *b_genL_px;    //!
  TBranch        *b_genL_py;    //!  
  TBranch        *b_genL_pz;    //!
  TBranch        *b_genL_e;    //!
  TBranch        *b_m_vis;   //!
  TBranch 		 *b_mt_tot;
  TBranch        *b_m_sv;   //!
  TBranch        *b_pt_sv;   //!
  TBranch        *b_eta_sv;   //!
  TBranch        *b_phi_sv;   //!
  TBranch        *b_met_sv;   //!
  TBranch        *b_mt_sv;   //!  
  TBranch        *b_pt_1;   //!
  TBranch        *b_phi_1;   //!
  TBranch        *b_eta_1;   //!
  TBranch        *b_m_1;   //!
  TBranch        *b_gen_match_1;   //!
  TBranch        *b_q_1;   //!
  TBranch        *b_iso_1;   //!
  TBranch        *b_mva_1;   //!
  TBranch        *b_d0_1;   //!
  TBranch        *b_dZ_1;   //!
  TBranch        *b_mt_1;   //!
  TBranch        *b_pfmt_1;   //!
  TBranch        *b_puppimt_1;   //!
  TBranch        *b_mt_rcqr_1;   //!
  TBranch        *b_pfmt_rcqr_1;   //!
  TBranch        *b_puppimt_rcqr_1;   //!
  TBranch        *b_mt_rcmr_1;   //!
  TBranch        *b_pfmt_rcmr_1;   //!
  TBranch        *b_puppimt_rcmr_1;   //!
  TBranch        *b_mt_rc_njetsreco_1;   //!
  TBranch        *b_pfmt_rc_njetsreco_1;   //!
  TBranch        *b_mt_rc_visreco_1;   //!
  TBranch        *b_pfmt_rc_visreco_1;   //!
  TBranch        *b_tau_decay_mode_1;   //!
  TBranch        *b_pt_2;   //!
  TBranch        *b_phi_2;   //!
  TBranch        *b_eta_2;   //!
  TBranch        *b_m_2;   //!
  TBranch        *b_gen_match_2;   //!
  TBranch        *b_q_2;   //!
  TBranch        *b_iso_2;   //!
  TBranch        *b_d0_2;   //!
  TBranch        *b_dZ_2;   //!
  TBranch        *b_mva_2;   //!
  TBranch        *b_mt_2;   //!
  TBranch        *b_pfmt_2;   //!
  TBranch        *b_puppimt_2;   //!
  TBranch        *b_mt_rcqr_2;   //!
  TBranch        *b_pfmt_rcqr_2;   //!
  TBranch        *b_puppimt_rcqr_2;   //!
  TBranch        *b_mt_rcmr_2;   //!
  TBranch        *b_pfmt_rcmr_2;   //!
  TBranch        *b_puppimt_rcmr_2;   //!
  TBranch        *b_mt_rc_njetsreco_2;   //!
  TBranch        *b_pfmt_rc_njetsreco_2;   //!
  TBranch        *b_mt_rc_visreco_2;   //!
  TBranch        *b_pfmt_rc_visreco_2;   //!
  TBranch        *b_tau_decay_mode_2;   //!
  TBranch        *b_os;   //!
  TBranch        *b_dilepton_veto;   //!
  TBranch        *b_extraelec_veto;   //!
  TBranch        *b_extramuon_veto;   //!
  TBranch 		 *b_chargedIsoPtSum_2; 
  TBranch        *b_neutralIsoPtSum_2;
  TBranch        *b_puCorrPtSum_2;
  TBranch        *b_isolationGammaCands_size_2;
  TBranch        *b_signalGammaCands_size_2;
  TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;   //!
  TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_1;   //!
  TBranch        *b_againstElectronLooseMVA5_1;   //!
  TBranch        *b_againstElectronMediumMVA5_1;   //!
  TBranch        *b_againstElectronTightMVA5_1;   //!
  TBranch        *b_againstElectronVLooseMVA5_1;   //!
  TBranch        *b_againstElectronVTightMVA5_1;   //!
  TBranch        *b_againstMuonLoose3_1;   //!
  TBranch        *b_againstMuonTight3_1;   //!
  TBranch        *b_n_badmuons;
  TBranch        *b_n_duplicatemuons;
  TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
  TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_2;   //!
  TBranch 		 *b_byIsolationMVArun2v1DBoldDMwLTraw_2;
  TBranch	     *b_byVLooseIsolationMVArun2v1DBoldDMwLT_2;
  TBranch	     *b_byLooseIsolationMVArun2v1DBoldDMwLT_2;
  TBranch	     *b_byMediumIsolationMVArun2v1DBoldDMwLT_2;
  TBranch	     *b_byTightIsolationMVArun2v1DBoldDMwLT_2;
  TBranch	     *b_byVTightIsolationMVArun2v1DBoldDMwLT_2;
  TBranch	     *b_byVVTightIsolationMVArun2v1DBoldDMwLT_2;
  TBranch        *b_againstElectronLooseMVA5_2;   //!
  TBranch        *b_againstElectronMediumMVA5_2;   //!
  TBranch        *b_againstElectronTightMVA5_2;   //!
  TBranch        *b_againstElectronVLooseMVA5_2;   //!
  TBranch        *b_againstElectronVTightMVA5_2;   //!
  TBranch	     *b_againstElectronVLooseMVA6_2;
  TBranch        *b_againstElectronTightMVA6_2;
  TBranch        *b_againstMuonLoose3_2;   //!
  TBranch        *b_againstMuonTight3_2;   //!
  TBranch        *b_met;   //!
  TBranch        *b_metphi;   //!
  TBranch        *b_metcov00;   //!
  TBranch        *b_metcov01;   //!
  TBranch        *b_metcov10;   //!
  TBranch        *b_metcov11;   //!
  TBranch        *b_met_rcqr;   //!
  TBranch        *b_metphi_rcqr;   //!
  TBranch        *b_met_rcmr;   //!
  TBranch        *b_metphi_rcmr;   //!
  TBranch        *b_met_rc_njetsreco;   //!
  TBranch        *b_metphi_rc_njetsreco;   //!
  TBranch        *b_met_rc_visreco;   //!
  TBranch        *b_metphi_rc_visreco;   //!
  TBranch        *b_mvamet;   //!
  TBranch        *b_mvametphi;   //!
  TBranch        *b_mvacov00;   //!
  TBranch        *b_mvacov01;   //!
  TBranch        *b_mvacov10;   //!
  TBranch        *b_mvacov11;   //!
  TBranch        *b_mvamet_rcqr;   //!
  TBranch        *b_mvametphi_rcqr;   //!
  TBranch        *b_mvamet_rcmr;   //!
  TBranch        *b_mvametphi_rcmr;   //!
  TBranch        *b_mvamet_rc_njetsreco;   //!
  TBranch        *b_mvametphi_rc_njetsreco;   //!
  TBranch        *b_mvamet_rc_visreco;   //!
  TBranch        *b_mvametphi_rc_visreco;   //!
  TBranch        *b_puppimet;   //!
  TBranch        *b_puppimetphi;   //!
  TBranch        *b_puppimet_rcqr;   //!
  TBranch        *b_puppimetphi_rcqr;   //!
  TBranch        *b_puppimet_rcmr;   //!
  TBranch        *b_puppimetphi_rcmr;   //!  
  TBranch        *b_pt_tt;   //!
  TBranch        *b_pzetavis;   //!
  TBranch        *b_pzetamiss;   //!
  TBranch        *b_pfpzetamiss;   //!
  TBranch        *b_puppipzetamiss;   //!
  TBranch        *b_pzetamiss_rcqr;   //!
  TBranch        *b_pfpzetamiss_rcqr;   //!
  TBranch        *b_puppipzetamiss_rcqr;   //!
  TBranch        *b_pzetamiss_rcmr;   //!
  TBranch        *b_pfpzetamiss_rcmr;   //!
  TBranch        *b_puppipzetamiss_rcmr;   //!
  TBranch        *b_pzetamiss_rc_njetsreco;   //!
  TBranch        *b_pfpzetamiss_rc_njetsreco;   //!
  TBranch        *b_pzetamiss_rc_visreco;   //!
  TBranch        *b_pfpzetamiss_rc_visreco;   //!
  TBranch        *b_mva_gf;   //!
  TBranch        *b_njets;   //!
  TBranch        *b_njetspt20;   //!
  TBranch        *b_jpt_1;   //!
  TBranch        *b_jeta_1;   //!
  TBranch        *b_jphi_1;   //!
  TBranch        *b_jptraw_1;   //!
  TBranch        *b_jrawf_1;   //!
  TBranch        *b_jptunc_1;   //!
  TBranch        *b_jmva_1;   //!
  TBranch        *b_jlrm_1;   //!
  TBranch        *b_jctm_1;   //!
  TBranch        *b_jpuid_loose_1;
  TBranch        *b_jpuid_medium_1;
  TBranch        *b_jpuid_tight_1;
  TBranch        *b_jpt_2;   //!
  TBranch        *b_jeta_2;   //!
  TBranch        *b_jphi_2;   //!
  TBranch        *b_jptraw_2;   //!
  TBranch        *b_jrawf_2;   //!
  TBranch        *b_jptunc_2;   //!
  TBranch        *b_jlrm_2;   //!
  TBranch        *b_jctm_2;   //!
  TBranch        *b_jpuid_loose_2;
  TBranch        *b_jpuid_medium_2;
  TBranch        *b_jpuid_tight_2;
  TBranch        *b_mjj;   //!
  TBranch        *b_jdeta;   //!
  TBranch        *b_njetingap;   //!
  TBranch        *b_nbtag;   //!
  TBranch        *b_bpt_1;   //!
  TBranch        *b_beta_1;   //!
  TBranch        *b_bphi_1;   //!
  TBranch        *b_brawf_1;   //!
  TBranch        *b_bmva_1;   //!
  TBranch        *b_bcsv_1;   //!
  TBranch        *b_bpuid_loose_1;
  TBranch        *b_bpuid_medium_1;
  TBranch        *b_bpuid_tight_1;
  TBranch        *b_bpt_2;   //!
  TBranch        *b_beta_2;   //!
  TBranch        *b_bphi_2;   //!
  TBranch        *b_brawf_2;   //!
  TBranch        *b_bmva_2;   //!
  TBranch        *b_bcsv_2;   //!
  TBranch        *b_bpuid_loose_2;
  TBranch        *b_bpuid_medium_2;
  TBranch        *b_bpuid_tight_2;
  TBranch        *b_singleLepTrigger;
  TBranch        *b_xTrigger;
  
  Spring15Tree(TTree *tree=0);
  virtual ~Spring15Tree();

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

#endif //!endif Spring15Tree_h
