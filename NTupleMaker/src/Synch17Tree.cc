/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza
//////////////////////////////////////////////////////////


#define Synch17Tree_cxx
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"

// Initialization
Synch17Tree::Synch17Tree(TTree *tree) : fChain(0) {
  lock = false;
  ientry = -2;
  
  if (!tree) return;
  
  Init(tree);
}

void Synch17Tree::Init(TTree *tree){
  if(!tree) return;
	      
  if(tree->GetEntries())
    ReadInit(tree);
  else
    WriteInit(tree);
}

// Destructors
Synch17Tree::~Synch17Tree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// Read methods
void Synch17Tree::ReadInit(TTree *tree)
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
   if (lock) return;
   ReadReset();
   
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("rho", &rho, &b_rho);

   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);

   fChain->SetBranchAddress("chconst_1_pt", &chconst_1_pt, &b_chconst_1_pt);
   fChain->SetBranchAddress("chconst_1_eta", &chconst_1_eta, &b_chconst_1_eta);
   fChain->SetBranchAddress("chconst_1_phi", &chconst_1_phi, &b_chconst_1_phi); 
   
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("mva_1", &mva_1, &b_mva_1);
   fChain->SetBranchAddress("mva17_1", &mva17_1, &b_mva17_1);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("d0err_1", &d0err_1, &b_d0err_1);
   fChain->SetBranchAddress("dZerr_1", &dZerr_1, &b_dZerr_1);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("tau_decay_mode_1", &tau_decay_mode_1, &b_tau_decay_mode_1); 
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);

   fChain->SetBranchAddress("chconst_2_pt", &chconst_2_pt, &b_chconst_2_pt);
   fChain->SetBranchAddress("chconst_2_eta", &chconst_2_eta, &b_chconst_2_eta);
   fChain->SetBranchAddress("chconst_2_phi", &chconst_2_phi, &b_chconst_2_phi); 
   
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("mva_2", &mva_2, &b_mva_2);
   fChain->SetBranchAddress("mva17_2", &mva17_2, &b_mva17_2);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("d0err_2", &d0err_2, &b_d0err_2);
   fChain->SetBranchAddress("dZerr_2", &dZerr_2, &b_dZerr_2);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("tau_decay_mode_2", &tau_decay_mode_2, &b_tau_decay_mode_2); 
   
   fChain->SetBranchAddress("trigweight_1", &trigweight_1, &b_trigweight_1);
   fChain->SetBranchAddress("trigweight_antiiso_1", &trigweight_antiiso_1, &b_trigweight_antiiso_1);
   fChain->SetBranchAddress("idisoweight_1", &idisoweight_1, &b_idisoweight_1);
   fChain->SetBranchAddress("idisoweight_antiiso_1", &idisoweight_antiiso_1, &b_idisoweight_antiiso_1);

   fChain->SetBranchAddress("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1, &b_againstElectronVLooseMVA6_1);
   fChain->SetBranchAddress("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1, &b_againstElectronVTightMVA6_1);
   fChain->SetBranchAddress("againstElectronTightMVA6_1", &againstElectronTightMVA6_1, &b_againstElectronTightMVA6_1);
   fChain->SetBranchAddress("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1, &b_againstElectronMediumMVA6_1);
   fChain->SetBranchAddress("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1, &b_againstElectronLooseMVA6_1);
   fChain->SetBranchAddress("againstMuonLoose3_1", &againstMuonLoose3_1, &b_againstMuonLoose3_1);
   fChain->SetBranchAddress("againstMuonTight3_1", &againstMuonTight3_1, &b_againstMuonTight3_1);
   fChain->SetBranchAddress("chargedIsoPtSum_1", &chargedIsoPtSum_1, &b_chargedIsoPtSum_1);
   fChain->SetBranchAddress("neutralIsoPtSum_1", &neutralIsoPtSum_1, &b_neutralIsoPtSum_1);
   fChain->SetBranchAddress("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1, &b_decayModeFindingOldDMs_1);
   fChain->SetBranchAddress("puCorrPtSum_1", &puCorrPtSum_1, &b_puCorrPtSum_1);

   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1, &b_byIsolationMVA3newDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1, &b_byIsolationMVA3oldDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1, &b_byIsolationMVA3newDMwLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1, &b_byIsolationMVA3oldDMwLTraw_1);

   fChain->SetBranchAddress("byIsolationMVArun2017v2DBoldDMwLTraw2017_1", &byIsolationMVArun2017v2DBoldDMwLTraw2017_1, &b_byIsolationMVArun2017v2DBoldDMwLTraw2017_1);
   fChain->SetBranchAddress("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byTightIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
//////////////////////////////////////////////////////////////NEW///////////////////////
   fChain->SetBranchAddress("efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
//   fChain->SetBranchAddress("efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
//   fChain->SetBranchAddress("efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1);
   fChain->SetBranchAddress("efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1);
//   fChain->SetBranchAddress("efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
//   fChain->SetBranchAddress("efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1, &b_efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
//////////////////////////////////////////////////////////////
   fChain->SetBranchAddress("efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
//   fChain->SetBranchAddress("efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
//   fChain->SetBranchAddress("efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2);
//   fChain->SetBranchAddress("efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
//   fChain->SetBranchAddress("efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
//////////////////////////////////////////////////////////////////////////////////////////
  fChain->SetBranchAddress("correction_againstElectronVLooseMVA6_1", &correction_againstElectronVLooseMVA6_1, &b_correction_againstElectronVLooseMVA6_1);
//  fChain->SetBranchAddress("correction_againstElectronLooseMVA6_1", &correction_againstElectronLooseMVA6_1, &b_correction_againstElectronLooseMVA6_1);
//  fChain->SetBranchAddress("correction_againstElectronMediumMVA6_1", &correction_againstElectronMediumMVA6_1, &b_correction_againstElectronMediumMVA6_1);
  fChain->SetBranchAddress("correction_againstElectronTightMVA6_1", &correction_againstElectronTightMVA6_1, &b_correction_againstElectronTightMVA6_1);
//  fChain->SetBranchAddress("correction_againstElectronVTightMVA6_1", &correction_againstElectronVTightMVA6_1, &b_correction_againstElectronVTightMVA6_1);
//////////////////////////////////////////////////////////////////////////////////////////
  fChain->SetBranchAddress("correction_againstElectronVLooseMVA6_2", &correction_againstElectronVLooseMVA6_2, &b_correction_againstElectronVLooseMVA6_2);
//  fChain->SetBranchAddress("correction_againstElectronLooseMVA6_2", &correction_againstElectronLooseMVA6_2, &b_correction_againstElectronLooseMVA6_2);
//  fChain->SetBranchAddress("correction_againstElectronMediumMVA6_2", &correction_againstElectronMediumMVA6_2, &b_correction_againstElectronMediumMVA6_2);
  fChain->SetBranchAddress("correction_againstElectronTightMVA6_2", &correction_againstElectronTightMVA6_2, &b_correction_againstElectronTightMVA6_2);
//  fChain->SetBranchAddress("correction_againstElectronVTightMVA6_2", &correction_againstElectronVTightMVA6_2, &b_correction_againstElectronVTightMVA6_2);
//////////////////////////////////////////////////////////////////////////////////////////
  fChain->SetBranchAddress("correction_againstMuonLoose3_1", &correction_againstMuonLoose3_1, &b_correction_againstMuonLoose3_1);
  fChain->SetBranchAddress("correction_againstMuonTight3_1", &correction_againstMuonTight3_1, &b_correction_againstMuonTight3_1);
//////////////////////////////////////////////////////////////////////////////////////////
  fChain->SetBranchAddress("correction_againstMuonLoose3_2", &correction_againstMuonLoose3_2, &b_correction_againstMuonLoose3_2);
  fChain->SetBranchAddress("correction_againstMuonTight3_2", &correction_againstMuonTight3_2, &b_correction_againstMuonTight3_2);
//////////////////////////////////////////////////////////////////////////////////////////
   fChain->SetBranchAddress("trigweight_2", &trigweight_2, &b_trigweight_2);
   fChain->SetBranchAddress("trigweight_antiiso_2", &trigweight_antiiso_2, &b_trigweight_antiiso_2);
   fChain->SetBranchAddress("idisoweight_2", &idisoweight_2, &b_idisoweight_2);
   fChain->SetBranchAddress("idisoweight_antiiso_2", &idisoweight_antiiso_2, &b_idisoweight_antiiso_2);

   fChain->SetBranchAddress("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, &b_againstElectronVLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, &b_againstElectronVTightMVA6_2);
   fChain->SetBranchAddress("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, &b_againstElectronTightMVA6_2);
   fChain->SetBranchAddress("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, &b_againstElectronMediumMVA6_2);
   fChain->SetBranchAddress("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, &b_againstElectronLooseMVA6_2);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("chargedIsoPtSum_2", &chargedIsoPtSum_2, &b_chargedIsoPtSum_2);
   fChain->SetBranchAddress("neutralIsoPtSum_2", &neutralIsoPtSum_2, &b_neutralIsoPtSum_2);
   fChain->SetBranchAddress("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, &b_decayModeFindingOldDMs_2);
   fChain->SetBranchAddress("puCorrPtSum_2", &puCorrPtSum_2, &b_puCorrPtSum_2);
   
   fChain->SetBranchAddress("byDeepTau2017v2p1VSeraw_2", &byDeepTau2017v2p1VSeraw_2, &b_byDeepTau2017v2p1VSeraw_2);
   fChain->SetBranchAddress("byDeepTau2017v2p1VSjetraw_2", &byDeepTau2017v2p1VSjetraw_2, &b_byDeepTau2017v2p1VSjetraw_2);
   fChain->SetBranchAddress("byDeepTau2017v2p1VSmuraw_2", &byDeepTau2017v2p1VSmuraw_2, &b_byDeepTau2017v2p1VSmuraw_2);
   fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSe_2", &byLooseDeepTau2017v2p1VSe_2, &b_byLooseDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSjet_2", &byLooseDeepTau2017v2p1VSjet_2, &b_byLooseDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSmu_2", &byLooseDeepTau2017v2p1VSmu_2, &b_byLooseDeepTau2017v2p1VSmu_2);
   fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSe_2", &byMediumDeepTau2017v2p1VSe_2, &b_byMediumDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSjet_2", &byMediumDeepTau2017v2p1VSjet_2, &b_byMediumDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSmu_2", &byMediumDeepTau2017v2p1VSmu_2, &b_byMediumDeepTau2017v2p1VSmu_2);
   fChain->SetBranchAddress("byTightDeepTau2017v2p1VSe_2", &byTightDeepTau2017v2p1VSe_2, &b_byTightDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byTightDeepTau2017v2p1VSjet_2", &byTightDeepTau2017v2p1VSjet_2, &b_byTightDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byTightDeepTau2017v2p1VSmu_2", &byTightDeepTau2017v2p1VSmu_2, &b_byTightDeepTau2017v2p1VSmu_2);
   fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSe_2", &byVLooseDeepTau2017v2p1VSe_2, &b_byVLooseDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSjet_2", &byVLooseDeepTau2017v2p1VSjet_2, &b_byVLooseDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSmu_2", &byVLooseDeepTau2017v2p1VSmu_2, &b_byVLooseDeepTau2017v2p1VSmu_2);
   fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSe_2", &byVTightDeepTau2017v2p1VSe_2, &b_byVTightDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSjet_2", &byVTightDeepTau2017v2p1VSjet_2, &b_byVTightDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSe_2", &byVVLooseDeepTau2017v2p1VSe_2, &b_byVVLooseDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSjet_2", &byVVLooseDeepTau2017v2p1VSjet_2, &b_byVVLooseDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSe_2", &byVVTightDeepTau2017v2p1VSe_2, &b_byVVTightDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSjet_2", &byVVTightDeepTau2017v2p1VSjet_2, &b_byVVTightDeepTau2017v2p1VSjet_2);
   fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSe_2", &byVVVLooseDeepTau2017v2p1VSe_2, &b_byVVVLooseDeepTau2017v2p1VSe_2);
   fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSjet_2", &byVVVLooseDeepTau2017v2p1VSjet_2, &b_byVVVLooseDeepTau2017v2p1VSjet_2);

   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, &b_byIsolationMVA3newDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, &b_byIsolationMVA3oldDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, &b_byIsolationMVA3newDMwLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, &b_byIsolationMVA3oldDMwLTraw_2);

   fChain->SetBranchAddress("byIsolationMVArun2017v2DBoldDMwLTraw2017_2", &byIsolationMVArun2017v2DBoldDMwLTraw2017_2, &b_byIsolationMVArun2017v2DBoldDMwLTraw2017_2);
   fChain->SetBranchAddress("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byTightIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
   fChain->SetBranchAddress("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2, &b_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("trigweight", &trigweight, &b_trigweight);

   fChain->SetBranchAddress("topptweight", &topptweight, &b_topptweight);
   fChain->SetBranchAddress("zptweight", &zptweight, &b_zptweight);
   fChain->SetBranchAddress("trkeffweight", &trkeffweight, &b_trkeffweight);
   fChain->SetBranchAddress("effweight", &effweight, &b_effweight); 
   fChain->SetBranchAddress("etaufakeweight", &etaufakeweight, &b_etaufakeweight);
   fChain->SetBranchAddress("mutaufakeweight", &mutaufakeweight, &b_mutaufakeweight);

   fChain->SetBranchAddress("xTrigger",  &xTrigger, &b_xTrigger);
   fChain->SetBranchAddress("xTriggerLep",  &xTriggerLep, &b_xTriggerLep);
   fChain->SetBranchAddress("xTriggerTau",  &xTriggerTau, &b_xTriggerTau);
   fChain->SetBranchAddress("trg_singlemuon",  &trg_singlemuon, &b_trg_singlemuon);
   fChain->SetBranchAddress("trg_singleelectron",  &trg_singleelectron, &b_trg_singleelectron);
   fChain->SetBranchAddress("singleLepTrigger",  &singleLepTrigger, &b_singleLepTrigger);
   fChain->SetBranchAddress("ditauTrigger", &ditauTrigger, &b_ditauTrigger);

   //MET
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("met_rcmr", &met_rcmr, &b_met_rcmr);
   fChain->SetBranchAddress("metphi_rcmr", &metphi_rcmr, &b_metphi_rcmr);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("pzetavis", &pzetavis, &b_pzetavis);
   fChain->SetBranchAddress("pzetamiss", &pzetamiss, &b_pzetamiss);
   
   //di tau system
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("mt_tot", &mt_tot, &b_mt_tot);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
   fChain->SetBranchAddress("eta_sv", &eta_sv, &b_eta_sv);
   fChain->SetBranchAddress("phi_sv", &phi_sv, &b_phi_sv);
   fChain->SetBranchAddress("met_sv", &met_sv, &b_met_sv);
   fChain->SetBranchAddress("mt_sv", &mt_sv, &b_mt_sv);   

   //VBF
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("jdeta", &jdeta, &b_jdeta);
   fChain->SetBranchAddress("dijetpt", &dijetpt, &b_dijetpt);
   fChain->SetBranchAddress("jdphi", &jdphi, &b_jdphi);
   fChain->SetBranchAddress("njetingap", &njetingap, &b_njetingap);
   fChain->SetBranchAddress("njetingap20", &njetingap20, &b_njetingap20);

   //jets
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njetshad", &njetshad, &b_njetshad);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   
   //b-jets
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("bpt_1", &bpt_1, &b_bpt_1);
   fChain->SetBranchAddress("beta_1", &beta_1, &b_beta_1);
   fChain->SetBranchAddress("bphi_1", &bphi_1, &b_bphi_1);
   fChain->SetBranchAddress("bcsv_1", &bcsv_1, &b_bcsv_1);
   fChain->SetBranchAddress("bpt_2", &bpt_2, &b_bpt_2);
   fChain->SetBranchAddress("beta_2", &beta_2, &b_beta_2);
   fChain->SetBranchAddress("bphi_2", &bphi_2, &b_bphi_2);
   fChain->SetBranchAddress("bcsv_2", &bcsv_2, &b_bcsv_2);

   //Misc   
   fChain->SetBranchAddress("gen_noutgoing", &gen_noutgoing, &b_gen_noutgoing);
   fChain->SetBranchAddress("os", &os, &b_os);
   fChain->SetBranchAddress("dilepton_veto", &dilepton_veto, &b_dilepton_veto);
   fChain->SetBranchAddress("extraelec_veto", &extraelec_veto, &b_extraelec_veto);
   fChain->SetBranchAddress("extramuon_veto", &extramuon_veto, &b_extramuon_veto);

   //CP measurement   
   fChain->SetBranchAddress("tau1DecayPlaneX", &tau1DecayPlaneX, &b_tau1DecayPlaneX);
   fChain->SetBranchAddress("tau1DecayPlaneY", &tau1DecayPlaneY, &b_tau1DecayPlaneY);
   fChain->SetBranchAddress("tau1DecayPlaneZ", &tau1DecayPlaneZ, &b_tau1DecayPlaneZ);
   fChain->SetBranchAddress("tau2DecayPlaneX", &tau2DecayPlaneX, &b_tau2DecayPlaneX);
   fChain->SetBranchAddress("tau2DecayPlaneY", &tau2DecayPlaneY, &b_tau2DecayPlaneY);
   fChain->SetBranchAddress("tau2DecayPlaneZ", &tau2DecayPlaneZ, &b_tau2DecayPlaneZ);
   fChain->SetBranchAddress("acotautau_00", &acotautau_00, &b_acotautau_00);
   fChain->SetBranchAddress("acotautau_10", &acotautau_10, &b_acotautau_10);
   fChain->SetBranchAddress("acotautau_01", &acotautau_01, &b_acotautau_01);
   fChain->SetBranchAddress("acotautau_11", &acotautau_11, &b_acotautau_11);

   fChain->SetBranchAddress("acotautauPsi_00", &acotautauPsi_00, &b_acotautauPsi_00);
   fChain->SetBranchAddress("acotautauPsi_10", &acotautauPsi_10, &b_acotautauPsi_10);
   fChain->SetBranchAddress("acotautauPsi_01", &acotautauPsi_01, &b_acotautauPsi_01);
   fChain->SetBranchAddress("acotautauPsi_11", &acotautauPsi_11, &b_acotautauPsi_11);


   fChain->SetBranchAddress("pdgcodetau2", &pdgcodetau2, &b_pdgcodetau2);

  //Points of closest approach

   fChain->SetBranchAddress("tau_pca2D_x_1", &tau_pca2D_x_1, &b_tau_pca2D_x_1);
   fChain->SetBranchAddress("tau_pca2D_y_1", &tau_pca2D_y_1, &b_tau_pca2D_y_1);
   fChain->SetBranchAddress("tau_pca2D_z_1", &tau_pca2D_z_1, &b_tau_pca2D_z_1);
   fChain->SetBranchAddress("tau_pca3D_x_1", &tau_pca3D_x_1, &b_tau_pca3D_x_1);
   fChain->SetBranchAddress("tau_pca3D_y_1", &tau_pca3D_y_1, &b_tau_pca3D_y_1);
   fChain->SetBranchAddress("tau_pca3D_z_1", &tau_pca3D_z_1, &b_tau_pca3D_z_1);

   fChain->SetBranchAddress("tau_pca2D_x_2", &tau_pca2D_x_2, &b_tau_pca2D_x_2);
   fChain->SetBranchAddress("tau_pca2D_y_2", &tau_pca2D_y_2, &b_tau_pca2D_y_2);
   fChain->SetBranchAddress("tau_pca2D_z_2", &tau_pca2D_z_2, &b_tau_pca2D_z_2);
   fChain->SetBranchAddress("tau_pca3D_x_2", &tau_pca3D_x_2, &b_tau_pca3D_x_2);
   fChain->SetBranchAddress("tau_pca3D_y_2", &tau_pca3D_y_2, &b_tau_pca3D_y_2);
   fChain->SetBranchAddress("tau_pca3D_z_2", &tau_pca3D_z_2, &b_tau_pca3D_z_2);

   //Secondary vertex
   fChain->SetBranchAddress("tau_SV_x_1", &tau_SV_x_1, &b_tau_SV_x_1);
   fChain->SetBranchAddress("tau_SV_y_1", &tau_SV_y_1, &b_tau_SV_y_1);
   fChain->SetBranchAddress("tau_SV_z_1", &tau_SV_z_1, &b_tau_SV_z_1);
   fChain->SetBranchAddress("tau_SV_covxx_1", &tau_SV_covxx_1, &b_tau_SV_covxx_1);
   fChain->SetBranchAddress("tau_SV_covyx_1", &tau_SV_covyx_1, &b_tau_SV_covyx_1);
   fChain->SetBranchAddress("tau_SV_covzx_1", &tau_SV_covzx_1, &b_tau_SV_covzx_1);
   fChain->SetBranchAddress("tau_SV_covyy_1", &tau_SV_covyy_1, &b_tau_SV_covyy_1);
   fChain->SetBranchAddress("tau_SV_covzy_1", &tau_SV_covzy_1, &b_tau_SV_covzy_1);
   fChain->SetBranchAddress("tau_SV_covzz_1", &tau_SV_covzz_1, &b_tau_SV_covzz_1);

   fChain->SetBranchAddress("tau_SV_x_2", &tau_SV_x_2, &b_tau_SV_x_2);
   fChain->SetBranchAddress("tau_SV_y_2", &tau_SV_y_2, &b_tau_SV_y_2);
   fChain->SetBranchAddress("tau_SV_z_2", &tau_SV_z_2, &b_tau_SV_z_2);
   fChain->SetBranchAddress("tau_SV_covxx_2", &tau_SV_covxx_2, &b_tau_SV_covxx_2);
   fChain->SetBranchAddress("tau_SV_covyx_2", &tau_SV_covyx_2, &b_tau_SV_covyx_2);
   fChain->SetBranchAddress("tau_SV_covzx_2", &tau_SV_covzx_2, &b_tau_SV_covzx_2);
   fChain->SetBranchAddress("tau_SV_covyy_2", &tau_SV_covyy_2, &b_tau_SV_covyy_2);
   fChain->SetBranchAddress("tau_SV_covzy_2", &tau_SV_covzy_2, &b_tau_SV_covzy_2);
   fChain->SetBranchAddress("tau_SV_covzz_2", &tau_SV_covzz_2, &b_tau_SV_covzz_2);


  //RECO vertex info useful to have
   fChain->SetBranchAddress("RecoVertexX", &RecoVertexX, &b_RecoVertexX);
   fChain->SetBranchAddress("RecoVertexY", &RecoVertexY, &b_RecoVertexY);
   fChain->SetBranchAddress("RecoVertexZ", &RecoVertexZ, &b_RecoVertexZ);

  //gen vertex info useful to have
   fChain->SetBranchAddress("GenVertexX", &GenVertexX, &b_GenVertexX);
   fChain->SetBranchAddress("GenVertexY", &GenVertexY, &b_GenVertexY);
   fChain->SetBranchAddress("GenVertexZ", &GenVertexZ, &b_GenVertexZ);
   
   fChain->SetBranchAddress("VxConstitTau1", &VxConstitTau1, &b_VxConstitTau1);
   fChain->SetBranchAddress("VyConstitTau1", &VyConstitTau1, &b_VyConstitTau1);
   fChain->SetBranchAddress("VzConstitTau1", &VzConstitTau1, &b_VzConstitTau1);

   fChain->SetBranchAddress("VxConstitTau2", &VxConstitTau2, &b_VxConstitTau2);
   fChain->SetBranchAddress("VyConstitTau2", &VyConstitTau2, &b_VyConstitTau2);
   fChain->SetBranchAddress("VzConstitTau2", &VzConstitTau2, &b_VzConstitTau2);
   fChain->SetBranchAddress("alphaminus", &alphaminus, &b_alphaminus);


   lock=true;
}

void Synch17Tree::ReadReset(){
  ientry = -2;
}

Long64_t Synch17Tree::GetEntries(){
  if (!fChain) return -1;
  
  return fChain->GetEntriesFast();
}

Long64_t Synch17Tree::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent)
    fCurrent = fChain->GetTreeNumber();

  ientry = entry;
  return centry;
}

Long64_t Synch17Tree::GetEntry(Long64_t entry){
  // Read contents of entry.
  // Sequential reading if entry < 0
  
  if (!fChain) return -1;

  if (entry < 0){
    if(ientry < 0) ientry = -1;

    ++ientry;
    entry = ientry;
  }

  if (LoadTree(entry) < 0) return -1;

  return fChain->GetEntry(entry);
}

Long64_t Synch17Tree::LoadedEntryId(){
  return ientry;
}
    

void Synch17Tree::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t Synch17Tree::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


//Write Methods
void Synch17Tree::WriteInit(TTree *tree) {
  if (!tree) return;
   
  fChain = tree;
  fCurrent = -1;

  fChain->Branch("run", &run, "run/i");
  fChain->Branch("lumi", &lumi, "lumi/i");
  fChain->Branch("evt", &evt, "evt/l");
  fChain->Branch("npv", &npv, "npv/I");
  fChain->Branch("npu", &npu, "npu/F");
  fChain->Branch("rho", &rho, "rho/F");

   fChain->Branch("pt_1", &pt_1, "pt_1/F");
   fChain->Branch("phi_1", &phi_1, "phi_1/F");
   fChain->Branch("eta_1", &eta_1, "eta_1/F");

   fChain->Branch("chconst_1_pt", &chconst_1_pt, "chconst_1_pt/F");
   fChain->Branch("chconst_1_eta", &chconst_1_eta, "chconst_1_eta/F");
   fChain->Branch("chconst_1_phi", &chconst_1_phi, "chconst_1_phi/F");
   
   fChain->Branch("m_1", &m_1, "m_1/F");
   fChain->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
   fChain->Branch("q_1", &q_1, "q_1/I");
   fChain->Branch("iso_1", &iso_1, "iso_1/F");
   fChain->Branch("mva_1", &mva_1, "mva_1/F");
   fChain->Branch("mva17_1", &mva17_1, "mva17_1/F");
   fChain->Branch("d0_1", &d0_1, "d0_1/F");
   fChain->Branch("dZ_1", &dZ_1, "dZ_1/F");
   fChain->Branch("d0err_1", &d0err_1, "d0err_1/F");
   fChain->Branch("dZerr_1", &dZerr_1, "dZerr_1/F");
   fChain->Branch("mt_1", &mt_1, "mt_1/F");
   fChain->Branch("tau_decay_mode_1", &tau_decay_mode_1, "tau_decay_mode_1/I"); 
   fChain->Branch("pt_2", &pt_2, "pt_2/F");
   fChain->Branch("phi_2", &phi_2, "phi_2/F");
   fChain->Branch("eta_2", &eta_2, "eta_2/F");

   fChain->Branch("chconst_2_pt", &chconst_2_pt, "chconst_2_pt/F");
   fChain->Branch("chconst_2_eta", &chconst_2_eta, "chconst_2_eta/F");
   fChain->Branch("chconst_2_phi", &chconst_2_phi, "chconst_2_phi/F");
   
   fChain->Branch("m_2", &m_2, "m_2/F");
   fChain->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
   fChain->Branch("q_2", &q_2, "q_2/I");
   fChain->Branch("iso_2", &iso_2, "iso_2/F");
   fChain->Branch("mva_2", &mva_2, "mva_2/F");
   fChain->Branch("mva17_2", &mva17_2, "mva17_2/F");
   fChain->Branch("d0_2", &d0_2, "d0_2/F");
   fChain->Branch("dZ_2", &dZ_2, "dZ_2/F");
   fChain->Branch("d0err_2", &d0err_2, "d0err_2/F");
   fChain->Branch("dZerr_2", &dZerr_2, "dZerr_2/F");
   fChain->Branch("mt_2", &mt_2, "mt_2/F");
   fChain->Branch("tau_decay_mode_2", &tau_decay_mode_2, "tau_decay_mode_2/I"); 
   
   fChain->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1, "againstElectronVLooseMVA6_1/F");
   fChain->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1, "againstElectronVTightMVA6_1/F");
   fChain->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1, "againstElectronTightMVA6_1/F");
   fChain->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1, "againstElectronMediumMVA6_1/F");
   fChain->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1, "againstElectronLooseMVA6_1/F");
   fChain->Branch("againstMuonLoose3_1", &againstMuonLoose3_1, "againstMuonLoose3_1/F");
   fChain->Branch("againstMuonTight3_1", &againstMuonTight3_1, "againstMuonTight3_1/F");
   fChain->Branch("chargedIsoPtSum_1", &chargedIsoPtSum_1, "chargedIsoPtSum_1/F");
   fChain->Branch("neutralIsoPtSum_1", &neutralIsoPtSum_1, "neutralIsoPtSum_1/F");
   fChain->Branch("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1, "decayModeFindingOldDMs_1/F");
   fChain->Branch("puCorrPtSum_1", &puCorrPtSum_1, "puCorrPtSum_1/F");

   fChain->Branch("byDeepTau2017v2p1VSeraw_2", &byDeepTau2017v2p1VSeraw_2, "byDeepTau2017v2p1VSeraw_2/F");
   fChain->Branch("byDeepTau2017v2p1VSjetraw_2", &byDeepTau2017v2p1VSjetraw_2, "byDeepTau2017v2p1VSjetraw_2/F");
   fChain->Branch("byDeepTau2017v2p1VSmuraw_2", &byDeepTau2017v2p1VSmuraw_2, "byDeepTau2017v2p1VSmuraw_2/F");
   fChain->Branch("byLooseDeepTau2017v2p1VSe_2", &byLooseDeepTau2017v2p1VSe_2, "byLooseDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byLooseDeepTau2017v2p1VSjet_2", &byLooseDeepTau2017v2p1VSjet_2, "byLooseDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byLooseDeepTau2017v2p1VSmu_2", &byLooseDeepTau2017v2p1VSmu_2, "byLooseDeepTau2017v2p1VSmu_2/F");
   fChain->Branch("byMediumDeepTau2017v2p1VSe_2", &byMediumDeepTau2017v2p1VSe_2, "byMediumDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byMediumDeepTau2017v2p1VSjet_2", &byMediumDeepTau2017v2p1VSjet_2, "byMediumDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byMediumDeepTau2017v2p1VSmu_2", &byMediumDeepTau2017v2p1VSmu_2, "byMediumDeepTau2017v2p1VSmu_2/F");
   fChain->Branch("byTightDeepTau2017v2p1VSe_2", &byTightDeepTau2017v2p1VSe_2, "byTightDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byTightDeepTau2017v2p1VSjet_2", &byTightDeepTau2017v2p1VSjet_2, "byTightDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byTightDeepTau2017v2p1VSmu_2", &byTightDeepTau2017v2p1VSmu_2, "byTightDeepTau2017v2p1VSmu_2/F");
   fChain->Branch("byVLooseDeepTau2017v2p1VSe_2", &byVLooseDeepTau2017v2p1VSe_2, "byVLooseDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byVLooseDeepTau2017v2p1VSjet_2", &byVLooseDeepTau2017v2p1VSjet_2, "byVLooseDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byVLooseDeepTau2017v2p1VSmu_2", &byVLooseDeepTau2017v2p1VSmu_2, "byVLooseDeepTau2017v2p1VSmu_2/F");
   fChain->Branch("byVTightDeepTau2017v2p1VSe_2", &byVTightDeepTau2017v2p1VSe_2, "byVTightDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byVTightDeepTau2017v2p1VSjet_2", &byVTightDeepTau2017v2p1VSjet_2, "byVTightDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byVVLooseDeepTau2017v2p1VSe_2", &byVVLooseDeepTau2017v2p1VSe_2, "byVVLooseDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byVVLooseDeepTau2017v2p1VSjet_2", &byVVLooseDeepTau2017v2p1VSjet_2, "byVVLooseDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byVVTightDeepTau2017v2p1VSe_2", &byVVTightDeepTau2017v2p1VSe_2, "byVVTightDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byVVTightDeepTau2017v2p1VSjet_2", &byVVTightDeepTau2017v2p1VSjet_2, "byVVTightDeepTau2017v2p1VSjet_2/F");
   fChain->Branch("byVVVLooseDeepTau2017v2p1VSe_2", &byVVVLooseDeepTau2017v2p1VSe_2, "byVVVLooseDeepTau2017v2p1VSe_2/F");
   fChain->Branch("byVVVLooseDeepTau2017v2p1VSjet_2", &byVVVLooseDeepTau2017v2p1VSjet_2, "byVVVLooseDeepTau2017v2p1VSjet_2/F");

   fChain->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, "byCombinedIsolationDeltaBetaCorrRaw3Hits_1/F");
   fChain->Branch("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1, "byIsolationMVA3newDMwoLTraw_1/F");
   fChain->Branch("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1, "byIsolationMVA3oldDMwoLTraw_1/F");
   fChain->Branch("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1, "byIsolationMVA3newDMwLTraw_1/F");
   fChain->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1, "byIsolationMVA3oldDMwLTraw_1/F");

  fChain->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_1", &byIsolationMVArun2017v2DBoldDMwLTraw2017_1, "byIsolationMVArun2017v2DBoldDMwLTraw2017_1/F");
  fChain->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1, "byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1/F");
  fChain->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_1, "byLooseIsolationMVArun2017v2DBoldDMwLT2017_1/F");
  fChain->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_1, "byMediumIsolationMVArun2017v2DBoldDMwLT2017_1/F");  
  fChain->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byTightIsolationMVArun2017v2DBoldDMwLT2017_1, "byTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
  fChain->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_1, "byVTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
  fChain->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1, "byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//////////////////////////////////////////////////////////////////////////////NEW
  fChain->Branch("efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//  fChain->Branch("efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//  fChain->Branch("efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1/F");
  fChain->Branch("efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//  fChain->Branch("efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//  fChain->Branch("efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1, "efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1/F");
//////////////////////////////////////////////////////////////////////////////
  fChain->Branch("efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2/F");
//  fChain->Branch("efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2/F");
//  fChain->Branch("efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2/F");
  fChain->Branch("efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");
//  fChain->Branch("efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");
//  fChain->Branch("efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2, "efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");
//////////////////////////////////////////////////////////////////////////////
  fChain->Branch("correction_againstElectronVLooseMVA6_1", &correction_againstElectronVLooseMVA6_1, "correction_againstElectronVLooseMVA6_1/F");
//  fChain->Branch("correction_againstElectronLooseMVA6_1", &correction_againstElectronLooseMVA6_1, "correction_againstElectronLooseMVA6_1/F");
//  fChain->Branch("correction_againstElectronMediumMVA6_1", &correction_againstElectronMediumMVA6_1, "correction_againstElectronMediumMVA6_1/F");
  fChain->Branch("correction_againstElectronTightMVA6_1", &correction_againstElectronTightMVA6_1, "correction_againstElectronTightMVA6_1/F");
//  fChain->Branch("correction_againstElectronVTightMVA6_1", &correction_againstElectronVTightMVA6_1, "correction_againstElectronVTightMVA6_1/F");
//////////////////////////////////////////////////////////////////////////////
  fChain->Branch("correction_againstElectronVLooseMVA6_2", &correction_againstElectronVLooseMVA6_2, "correction_againstElectronVLooseMVA6_2/F");
//  fChain->Branch("correction_againstElectronLooseMVA6_2", &correction_againstElectronLooseMVA6_2, "correction_againstElectronLooseMVA6_2/F");
//  fChain->Branch("correction_againstElectronMediumMVA6_2", &correction_againstElectronMediumMVA6_2, "correction_againstElectronMediumMVA6_2/F");
  fChain->Branch("correction_againstElectronTightMVA6_2", &correction_againstElectronTightMVA6_2, "correction_againstElectronTightMVA6_2/F");
//  fChain->Branch("correction_againstElectronVTightMVA6_2", &correction_againstElectronVTightMVA6_2, "correction_againstElectronVTightMVA6_2/F");
//////////////////////////////////////////////////////////////////////////////
  fChain->Branch("correction_againstMuonLoose3_1", &correction_againstMuonLoose3_1, "correction_againstMuonLoose3_1/F");
  fChain->Branch("correction_againstMuonTight3_1", &correction_againstMuonTight3_1, "correction_againstMuonTight3_1/F");
//////////////////////////////////////////////////////////////////////////////
  fChain->Branch("correction_againstMuonLoose3_2", &correction_againstMuonLoose3_2, "correction_againstMuonLoose3_2/F");
  fChain->Branch("correction_againstMuonTight3_2", &correction_againstMuonTight3_2, "correction_againstMuonTight3_2/F");
//////////////////////////////////////////////////////////////////////////////

   fChain->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
   fChain->Branch("trigweight_antiiso_1", &trigweight_antiiso_1, "trigweight_antiiso_1/F");
   fChain->Branch("idisoweight_1", &idisoweight_1, "idisoweight_1/F");
   fChain->Branch("idisoweight_antiiso_1", &idisoweight_antiiso_1, "idisoweight_antiiso_1/F");
   fChain->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
   fChain->Branch("trigweight_antiiso_2", &trigweight_antiiso_2, "trigweight_antiiso_2/F");
   fChain->Branch("idisoweight_2", &idisoweight_2, "idisoweight_2/F");
   fChain->Branch("idisoweight_antiiso_2", &idisoweight_antiiso_2, "idisoweight_antiiso_2/F");
   fChain->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, "againstElectronVLooseMVA6_2/F");
   fChain->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, "againstElectronVTightMVA6_2/F");
   fChain->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, "againstElectronTightMVA6_2/F");
   fChain->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, "againstElectronMediumMVA6_2/F");
   fChain->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, "againstElectronLooseMVA6_2/F");
   fChain->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/F");
   fChain->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");
   fChain->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2, "chargedIsoPtSum_2/F");
   fChain->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2, "neutralIsoPtSum_2/F");
   fChain->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, "decayModeFindingOldDMs_2/F");
   fChain->Branch("puCorrPtSum_2", &puCorrPtSum_2, "puCorrPtSum_2/F");

   fChain->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
   fChain->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, "byIsolationMVA3newDMwoLTraw_2/F");
   fChain->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, "byIsolationMVA3oldDMwoLTraw_2/F");
   fChain->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, "byIsolationMVA3newDMwLTraw_2/F");
   fChain->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, "byIsolationMVA3oldDMwLTraw_2/F");

   fChain->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_2", &byIsolationMVArun2017v2DBoldDMwLTraw2017_2, "byIsolationMVArun2017v2DBoldDMwLTraw2017_2/F");
   fChain->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2, "byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2/F");
   fChain->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_2, "byLooseIsolationMVArun2017v2DBoldDMwLT2017_2/F");
   fChain->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_2, "byMediumIsolationMVArun2017v2DBoldDMwLT2017_2/F");  
   fChain->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byTightIsolationMVArun2017v2DBoldDMwLT2017_2, "byTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");
   fChain->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_2, "byVTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");
   fChain->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2, "byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2/F");

   fChain->Branch("weight", &weight, "weight/F");
   fChain->Branch("mcweight", &mcweight, "mcweight/F");
   fChain->Branch("puweight", &puweight, "puweight/F");
   fChain->Branch("trigweight", &trigweight, "trigweight/F");

   fChain->Branch("topptweight", &topptweight, "topptweight/F");
   fChain->Branch("zptweight", &zptweight, "zptweight/D");
   fChain->Branch("trkeffweight", &trkeffweight, "trkeffweight/D");
   fChain->Branch("effweight", &effweight, "effweight/F"); 
   fChain->Branch("etaufakeweight", &etaufakeweight, "etaufakeweight/F");
   fChain->Branch("mutaufakeweight", &mutaufakeweight, "mutaufakeweight/F");

   fChain->Branch("xTrigger",  &xTrigger, "xTrigger/O");
   fChain->Branch("xTriggerLep",  &xTriggerLep, "xTriggerLep/O");
   fChain->Branch("xTriggerTau",  &xTriggerTau, "xTriggerTau/O");
   fChain->Branch("trg_singlemuon",  &trg_singlemuon, "trg_singlemuon/O");
   fChain->Branch("trg_singleelectron",  &trg_singleelectron, "trg_singleelectron/O");
   fChain->Branch("singleLepTrigger",  &singleLepTrigger, "singleLepTrigger/O");
   fChain->Branch("ditauTrigger", &ditauTrigger,"ditauTrigger/O");

   //MET
   fChain->Branch("met", &met, "met/F");
   fChain->Branch("metphi", &metphi, "metphi/F");
   fChain->Branch("met_rcmr", &met_rcmr, "met_rcmr/F");
   fChain->Branch("metphi_rcmr", &metphi_rcmr, "metphi_rcmr/F");
   fChain->Branch("metcov00", &metcov00, "metcov00/F");
   fChain->Branch("metcov01", &metcov01, "metcov01/F");
   fChain->Branch("metcov10", &metcov10, "metcov10/F");
   fChain->Branch("metcov11", &metcov11, "metcov11/F");
   fChain->Branch("pzetavis", &pzetavis, "pzetavis/F");
   fChain->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
   
   //di tau system
   fChain->Branch("pt_tt", &pt_tt, "pt_tt/F");
   fChain->Branch("m_vis", &m_vis, "m_vis/F");
   fChain->Branch("mt_tot", &mt_tot, "mt_tot/F");
   fChain->Branch("m_sv", &m_sv, "m_sv/F");
   fChain->Branch("pt_sv", &pt_sv, "pt_sv/F");
   fChain->Branch("eta_sv", &eta_sv, "eta_sv/F");
   fChain->Branch("phi_sv", &phi_sv, "phi_sv/F");
   fChain->Branch("met_sv", &met_sv, "met_sv/F");
   fChain->Branch("mt_sv", &mt_sv, "mt_sv/F");   

   //VBF
   fChain->Branch("mjj", &mjj, "mjj/F");
   fChain->Branch("jdeta", &jdeta, "jdeta/F");
   fChain->Branch("dijetpt", &dijetpt, "dijetpt/F");
   fChain->Branch("jdphi", &jdphi, "jdphi/F");
   fChain->Branch("njetingap", &njetingap, "njetingap/I");
   fChain->Branch("njetingap20", &njetingap20, "njetingap20/I");

   //jets
   fChain->Branch("njets", &njets, "njets/I");
   fChain->Branch("njetshad", &njetshad, "njetshad/I");
   fChain->Branch("njetspt20", &njetspt20, "njetspt20/I");
   fChain->Branch("jpt_1", &jpt_1, "jpt_1/F");
   fChain->Branch("jeta_1", &jeta_1, "jeta_1/F");
   fChain->Branch("jphi_1", &jphi_1, "jphi_1/F");
   fChain->Branch("jpt_2", &jpt_2, "jpt_2/F");
   fChain->Branch("jeta_2", &jeta_2, "jeta_2/F");
   fChain->Branch("jphi_2", &jphi_2, "jphi_2/F");
   
   //b-jets
   fChain->Branch("nbtag", &nbtag, "nbtag/I");
   fChain->Branch("bpt_1", &bpt_1, "bpt_1/F");
   fChain->Branch("beta_1", &beta_1, "beta_1/F");
   fChain->Branch("bphi_1", &bphi_1, "bphi_1/F");
   fChain->Branch("bcsv_1", &bcsv_1, "bcsv_1/F");
   fChain->Branch("bpt_2", &bpt_2, "bpt_2/F");
   fChain->Branch("beta_2", &beta_2, "beta_2/F");
   fChain->Branch("bphi_2", &bphi_2, "bphi_2/F");
   fChain->Branch("bcsv_2", &bcsv_2, "bcsv_2/F");

   //Misc   
   fChain->Branch("gen_noutgoing", &gen_noutgoing, "gen_noutgoing/I");
   fChain->Branch("os", &os, "os/I");
   fChain->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/I");
   fChain->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/I");
   fChain->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/I");

   //CP measurement

   fChain->Branch("tau1DecayPlaneX", &tau1DecayPlaneX, "tau1DecayPlaneX/F");
   fChain->Branch("tau1DecayPlaneY", &tau1DecayPlaneY, "tau1DecayPlaneY/F");
   fChain->Branch("tau1DecayPlaneZ", &tau1DecayPlaneZ, "tau1DecayPlaneZ/F");
   fChain->Branch("tau2DecayPlaneX", &tau2DecayPlaneX, "tau2DecayPlaneX/F");
   fChain->Branch("tau2DecayPlaneY", &tau2DecayPlaneY, "tau2DecayPlaneY/F");
   fChain->Branch("tau2DecayPlaneZ", &tau2DecayPlaneZ, "tau2DecayPlaneZ/F");
   fChain->Branch("acotautau_00", &acotautau_00, "acotautau_00/F");
   fChain->Branch("acotautau_10", &acotautau_10, "acotautau_10/F");
   fChain->Branch("acotautau_01", &acotautau_01, "acotautau_01/F");
   fChain->Branch("acotautau_11", &acotautau_11, "acotautau_11/F");

   fChain->Branch("acotautauPsi_00", &acotautauPsi_00, "acotautauPsi_00/F");
   fChain->Branch("acotautauPsi_10", &acotautauPsi_10, "acotautauPsi_10/F");
   fChain->Branch("acotautauPsi_01", &acotautauPsi_01, "acotautauPsi_01/F");
   fChain->Branch("acotautauPsi_11", &acotautauPsi_11, "acotautauPsi_11/F");

   fChain->Branch("pdgcodetau2", &pdgcodetau2, "pdgcodetau2/F");

  //Points of closest approach

   fChain->Branch("tau_pca2D_x_1", &tau_pca2D_x_1, "tau_pca2D_x_1/F");
   fChain->Branch("tau_pca2D_y_1", &tau_pca2D_y_1, "tau_pca2D_y_1/F");
   fChain->Branch("tau_pca2D_z_1", &tau_pca2D_z_1, "tau_pca2D_z_1/F");
   fChain->Branch("tau_pca3D_x_1", &tau_pca3D_x_1, "tau_pca3D_x_1/F");
   fChain->Branch("tau_pca3D_y_1", &tau_pca3D_y_1, "tau_pca3D_y_1/F");
   fChain->Branch("tau_pca3D_z_1", &tau_pca3D_z_1, "tau_pca3D_z_1/F");
   fChain->Branch("tau_pca2D_x_2", &tau_pca2D_x_2, "tau_pca2D_x_2/F");
   fChain->Branch("tau_pca2D_y_2", &tau_pca2D_y_2, "tau_pca2D_y_2/F");
   fChain->Branch("tau_pca2D_z_2", &tau_pca2D_z_2, "tau_pca2D_z_2/F");
   fChain->Branch("tau_pca3D_x_2", &tau_pca3D_x_2, "tau_pca3D_x_2/F");
   fChain->Branch("tau_pca3D_y_2", &tau_pca3D_y_2, "tau_pca3D_y_2/F");
   fChain->Branch("tau_pca3D_z_2", &tau_pca3D_z_2, "tau_pca3D_z_2/F");

   fChain->Branch("tau_SV_x_1", &tau_SV_x_1, "tau_SV_x_1/F");
   fChain->Branch("tau_SV_y_1", &tau_SV_y_1, "tau_SV_y_1/F");
   fChain->Branch("tau_SV_z_1", &tau_SV_z_1, "tau_SV_z_1/F");
   fChain->Branch("tau_SV_covxx_1", &tau_SV_covxx_1, "tau_SV_covxx_1/F");
   fChain->Branch("tau_SV_covyx_1", &tau_SV_covyx_1, "tau_SV_covyx_1/F");
   fChain->Branch("tau_SV_covzx_1", &tau_SV_covzx_1, "tau_SV_covzx_1/F");
   fChain->Branch("tau_SV_covyy_1", &tau_SV_covyy_1, "tau_SV_covyy_1/F");
   fChain->Branch("tau_SV_covzy_1", &tau_SV_covzy_1, "tau_SV_covzy_1/F");
   fChain->Branch("tau_SV_covzz_1", &tau_SV_covzz_1, "tau_SV_covzz_1/F");

   fChain->Branch("tau_SV_x_2", &tau_SV_x_2, "tau_SV_x_2/F");
   fChain->Branch("tau_SV_y_2", &tau_SV_y_2, "tau_SV_y_2/F");
   fChain->Branch("tau_SV_z_2", &tau_SV_z_2, "tau_SV_z_2/F");
   fChain->Branch("tau_SV_covxx_2", &tau_SV_covxx_2, "tau_SV_covxx_2/F");
   fChain->Branch("tau_SV_covyx_2", &tau_SV_covyx_2, "tau_SV_covyx_2/F");
   fChain->Branch("tau_SV_covzx_2", &tau_SV_covzx_2, "tau_SV_covzx_2/F");
   fChain->Branch("tau_SV_covyy_2", &tau_SV_covyy_2, "tau_SV_covyy_2/F");
   fChain->Branch("tau_SV_covzy_2", &tau_SV_covzy_2, "tau_SV_covzy_2/F");
   fChain->Branch("tau_SV_covzz_2", &tau_SV_covzz_2, "tau_SV_covzz_2/F");

   fChain->Branch("RecoVertexX", &RecoVertexX, "RecoVertexX/F");
   fChain->Branch("RecoVertexY", &RecoVertexY, "RecoVertexY/F");
   fChain->Branch("RecoVertexZ", &RecoVertexZ, "RecoVertexZ/F"); 

   fChain->Branch("GenVertexX", &GenVertexX, "GenVertexX/F");
   fChain->Branch("GenVertexY", &GenVertexY, "GenVertexY/F");
   fChain->Branch("GenVertexZ", &GenVertexZ, "GenVertexZ/F");

   fChain->Branch("VxConstitTau1", &VxConstitTau1, "VxConstitTau1/F");
   fChain->Branch("VyConstitTau1", &VyConstitTau1, "VyConstitTau1/F");
   fChain->Branch("VzConstitTau1", &VzConstitTau1, "VzConstitTau1/F");

   fChain->Branch("VxConstitTau2", &VxConstitTau2, "VxConstitTau2/F");
   fChain->Branch("VyConstitTau2", &VyConstitTau2, "VyConstitTau2/F");
   fChain->Branch("VzConstitTau2", &VzConstitTau2, "VzConstitTau2/F");
   fChain->Branch("alphaminus", &alphaminus, "alphaminus/F");

   
}

void Synch17Tree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}
