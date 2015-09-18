/////////////////////////////////////////////////////////////
// Read and Write Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
// 
// Based on a class automatically generated on
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
// from TTree TauCheck/TauCheck
// found on file: VBFHToTauTau_M-125_Synch.root
//////////////////////////////////////////////////////////

#ifndef Phys14Tree_h
#define Phys14Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Phys14Tree {
public :
  
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain
  
  // Declaration of leaf types
  Int_t           run;
  Int_t           lumi;
  Int_t           evt;
  Int_t           npv;
  Int_t           npu;
  Float_t         rho;
  Float_t         mcweight;
  Float_t         puweight;
  Float_t         trigweight_1;
  Float_t         trigweight_2;
  Float_t         idweight_1;
  Float_t         idweight_2;
  Float_t         isoweight_1;
  Float_t         isoweight_2;
  Float_t         effweight;
  Float_t         fakeweight;
  Float_t         embeddedWeight;
  Float_t         signalWeight;
  Float_t         weight;
  Float_t         m_vis;
  Float_t         m_sv;
  Float_t         pt_sv;
  Float_t         eta_sv;
  Float_t         phi_sv;
  Float_t         pt_1;
  Float_t         phi_1;
  Float_t         eta_1;
  Float_t         m_1;
  Int_t           q_1;
  Float_t         iso_1;
  Float_t         mva_1;
  Float_t         d0_1;
  Float_t         dZ_1;
  Float_t         mt_1;
  Float_t         pt_2;
  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         m_2;
  Int_t           q_2;
  Float_t         iso_2;
  Float_t         d0_2;
  Float_t         dZ_2;
  Float_t         mva_2;
  Float_t         mt_2;
  Char_t          os;
  Char_t          dilepton_veto;
  Char_t          extraelec_veto;
  Char_t          extramuon_veto;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t         againstElectronLooseMVA5_1;
  Float_t         againstElectronMediumMVA5_1;
  Float_t         againstElectronTightMVA5_1;
  Float_t         againstElectronVLooseMVA5_1;
  Float_t         againstElectronVTightMVA5_1;
  Float_t         againstMuonLoose3_1;
  Float_t         againstMuonTight3_1;
  Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
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
  Float_t         mvamet;
  Float_t         mvametphi;
  Float_t         mvacov00;
  Float_t         mvacov01;
  Float_t         mvacov10;
  Float_t         mvacov11;
  Float_t         pt_tt;
  Float_t         pzetavis;
  Float_t         pzetamiss;
  Float_t         mva_gf;
  Int_t           njets;
  Int_t           njetspt20;
  Float_t         jpt_1;
  Float_t         jeta_1;
  Float_t         jphi_1;
  Float_t         jptraw_1;
  Float_t         jptunc_1;
  Float_t         jmva_1;
  Float_t         jlrm_1;
  Int_t           jctm_1;
  Float_t         jpt_2;
  Float_t         jeta_2;
  Float_t         jphi_2;
  Float_t         jptraw_2;
  Float_t         jptunc_2;
  Float_t         jmva_2;
  Float_t         jlrm_2;
  Int_t           jctm_2;
  Float_t         mjj;
  Float_t         jdeta;
  Int_t           njetingap;
  Int_t           nbtag;
  Float_t         bpt;
  Float_t         beta;
  Float_t         bphi;
  
  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_evt;   //!
  TBranch        *b_npv;   //!
  TBranch        *b_npu;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_mcweight;   //!
  TBranch        *b_puweight;   //!
  TBranch        *b_trigweight_1;   //!
  TBranch        *b_trigweight_2;   //!
  TBranch        *b_idweight_1;   //!
  TBranch        *b_idweight_2;   //!
  TBranch        *b_isoweight_1;   //!
  TBranch        *b_isoweight_2;   //!
  TBranch        *b_effweight;   //!
  TBranch        *b_fakeweight;   //!
  TBranch        *b_embeddedWeight;   //!
  TBranch        *b_signalWeight;   //!
  TBranch        *b_weight;   //!
  TBranch        *b_m_vis;   //!
  TBranch        *b_m_sv;   //!
  TBranch        *b_pt_sv;   //!
  TBranch        *b_eta_sv;   //!
  TBranch        *b_phi_sv;   //!
  TBranch        *b_pt_1;   //!
  TBranch        *b_phi_1;   //!
  TBranch        *b_eta_1;   //!
  TBranch        *b_m_1;   //!
  TBranch        *b_q_1;   //!
  TBranch        *b_iso_1;   //!
  TBranch        *b_mva_1;   //!
  TBranch        *b_d0_1;   //!
  TBranch        *b_dZ_1;   //!
  TBranch        *b_mt_1;   //!
  TBranch        *b_pt_2;   //!
  TBranch        *b_phi_2;   //!
  TBranch        *b_eta_2;   //!
  TBranch        *b_m_2;   //!
  TBranch        *b_q_2;   //!
  TBranch        *b_iso_2;   //!
  TBranch        *b_d0_2;   //!
  TBranch        *b_dZ_2;   //!
  TBranch        *b_mva_2;   //!
  TBranch        *b_mt_2;   //!
  TBranch        *b_os;   //!
  TBranch        *b_dilepton_veto;   //!
  TBranch        *b_extraelec_veto;   //!
  TBranch        *b_extramuon_veto;   //!
  TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;   //!
  TBranch        *b_againstElectronLooseMVA5_1;   //!
  TBranch        *b_againstElectronMediumMVA5_1;   //!
  TBranch        *b_againstElectronTightMVA5_1;   //!
  TBranch        *b_againstElectronVLooseMVA5_1;   //!
  TBranch        *b_againstElectronVTightMVA5_1;   //!
  TBranch        *b_againstMuonLoose3_1;   //!
  TBranch        *b_againstMuonTight3_1;   //!
  TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
  TBranch        *b_againstElectronLooseMVA5_2;   //!
  TBranch        *b_againstElectronMediumMVA5_2;   //!
  TBranch        *b_againstElectronTightMVA5_2;   //!
  TBranch        *b_againstElectronVLooseMVA5_2;   //!
  TBranch        *b_againstElectronVTightMVA5_2;   //!
  TBranch        *b_againstMuonLoose3_2;   //!
  TBranch        *b_againstMuonTight3_2;   //!
  TBranch        *b_met;   //!
  TBranch        *b_metphi;   //!
  TBranch        *b_metcov00;   //!
  TBranch        *b_metcov01;   //!
  TBranch        *b_metcov10;   //!
  TBranch        *b_metcov11;   //!
  TBranch        *b_mvamet;   //!
  TBranch        *b_mvametphi;   //!
  TBranch        *b_mvacov00;   //!
  TBranch        *b_mvacov01;   //!
  TBranch        *b_mvacov10;   //!
  TBranch        *b_mvacov11;   //!
  TBranch        *b_pt_tt;   //!
  TBranch        *b_pzetavis;   //!
  TBranch        *b_pzetamiss;   //!
  TBranch        *b_mva_gf;   //!
  TBranch        *b_njets;   //!
  TBranch        *b_njetspt20;   //!
  TBranch        *b_jpt_1;   //!
  TBranch        *b_jeta_1;   //!
  TBranch        *b_jphi_1;   //!
  TBranch        *b_jptraw_1;   //!
  TBranch        *b_jptunc_1;   //!
  TBranch        *b_jmva_1;   //!
  TBranch        *b_jlrm_1;   //!
  TBranch        *b_jctm_1;   //!
  TBranch        *b_jpt_2;   //!
  TBranch        *b_jeta_2;   //!
  TBranch        *b_jphi_2;   //!
  TBranch        *b_jptraw_2;   //!
  TBranch        *b_jptunc_2;   //!
  TBranch        *b_jlrm_2;   //!
  TBranch        *b_jctm_2;   //!
  TBranch        *b_mjj;   //!
  TBranch        *b_jdeta;   //!
  TBranch        *b_njetingap;   //!
  TBranch        *b_nbtag;   //!
  TBranch        *b_bpt;   //!
  TBranch        *b_beta;   //!
  TBranch        *b_bphi;   //!
  
  Phys14Tree(TTree *tree=0);
  virtual ~Phys14Tree();

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

#endif //!endif Phys14Tree_h
