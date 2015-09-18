/////////////////////////////////////////////////////////////
// Read and Write Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
// 
// Based on a class automatically generated on
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
// from TTree TauCheck/TauCheck
// found on file: VBFHToTauTau_M-125_Synch.root
/////////////////////////////////////////////////////////////

#define Phys14Tree_cxx
#include "DesyTauAnalyses/NTupleMaker/interface/Phys14Tree.h"

// Initialization
Phys14Tree::Phys14Tree(TTree *tree) : fChain(0) {
  lock = false;
  ientry = -2;
  
  if (!tree) return;
  
  Init(tree);
}

void Phys14Tree::Init(TTree *tree){
  if(!tree) return;
	      
  if(tree->GetEntries())
    ReadInit(tree);
  else
    WriteInit(tree);
}

// Destructors
Phys14Tree::~Phys14Tree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// Read methods
void Phys14Tree::ReadInit(TTree *tree)
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
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("trigweight_1", &trigweight_1, &b_trigweight_1);
   fChain->SetBranchAddress("trigweight_2", &trigweight_2, &b_trigweight_2);
   fChain->SetBranchAddress("idweight_1", &idweight_1, &b_idweight_1);
   fChain->SetBranchAddress("idweight_2", &idweight_2, &b_idweight_2);
   fChain->SetBranchAddress("isoweight_1", &isoweight_1, &b_isoweight_1);
   fChain->SetBranchAddress("isoweight_2", &isoweight_2, &b_isoweight_2);
   fChain->SetBranchAddress("effweight", &effweight, &b_effweight);
   fChain->SetBranchAddress("fakeweight", &fakeweight, &b_fakeweight);
   fChain->SetBranchAddress("embeddedWeight", &embeddedWeight, &b_embeddedWeight);
   fChain->SetBranchAddress("signalWeight", &signalWeight, &b_signalWeight);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
   fChain->SetBranchAddress("eta_sv", &eta_sv, &b_eta_sv);
   fChain->SetBranchAddress("phi_sv", &phi_sv, &b_phi_sv);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("mva_1", &mva_1, &b_mva_1);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("mva_2", &mva_2, &b_mva_2);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("os", &os, &b_os);
   fChain->SetBranchAddress("dilepton_veto", &dilepton_veto, &b_dilepton_veto);
   fChain->SetBranchAddress("extraelec_veto", &extraelec_veto, &b_extraelec_veto);
   fChain->SetBranchAddress("extramuon_veto", &extramuon_veto, &b_extramuon_veto);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
   fChain->SetBranchAddress("againstElectronLooseMVA5_1", &againstElectronLooseMVA5_1, &b_againstElectronLooseMVA5_1);
   fChain->SetBranchAddress("againstElectronMediumMVA5_1", &againstElectronMediumMVA5_1, &b_againstElectronMediumMVA5_1);
   fChain->SetBranchAddress("againstElectronTightMVA5_1", &againstElectronTightMVA5_1, &b_againstElectronTightMVA5_1);
   fChain->SetBranchAddress("againstElectronVLooseMVA5_1", &againstElectronVLooseMVA5_1, &b_againstElectronVLooseMVA5_1);
   fChain->SetBranchAddress("againstElectronVTightMVA5_1", &againstElectronVTightMVA5_1, &b_againstElectronVTightMVA5_1);
   fChain->SetBranchAddress("againstMuonLoose3_1", &againstMuonLoose3_1, &b_againstMuonLoose3_1);
   fChain->SetBranchAddress("againstMuonTight3_1", &againstMuonTight3_1, &b_againstMuonTight3_1);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2, &b_againstElectronLooseMVA5_2);
   fChain->SetBranchAddress("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2, &b_againstElectronMediumMVA5_2);
   fChain->SetBranchAddress("againstElectronTightMVA5_2", &againstElectronTightMVA5_2, &b_againstElectronTightMVA5_2);
   fChain->SetBranchAddress("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2, &b_againstElectronVLooseMVA5_2);
   fChain->SetBranchAddress("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2, &b_againstElectronVTightMVA5_2);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("mvamet", &mvamet, &b_mvamet);
   fChain->SetBranchAddress("mvametphi", &mvametphi, &b_mvametphi);
   fChain->SetBranchAddress("mvacov00", &mvacov00, &b_mvacov00);
   fChain->SetBranchAddress("mvacov01", &mvacov01, &b_mvacov01);
   fChain->SetBranchAddress("mvacov10", &mvacov10, &b_mvacov10);
   fChain->SetBranchAddress("mvacov11", &mvacov11, &b_mvacov11);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("pzetavis", &pzetavis, &b_pzetavis);
   fChain->SetBranchAddress("pzetamiss", &pzetamiss, &b_pzetamiss);
   fChain->SetBranchAddress("mva_gf", &mva_gf, &b_mva_gf);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jptraw_1", &jptraw_1, &b_jptraw_1);
   fChain->SetBranchAddress("jptunc_1", &jptunc_1, &b_jptunc_1);
   fChain->SetBranchAddress("jmva_1", &jmva_1, &b_jmva_1);
   fChain->SetBranchAddress("jlrm_1", &jlrm_1, &b_jlrm_1);
   fChain->SetBranchAddress("jctm_1", &jctm_1, &b_jctm_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("jptraw_2", &jptraw_2, &b_jptraw_2);
   fChain->SetBranchAddress("jptunc_2", &jptunc_2, &b_jptunc_2);
   fChain->SetBranchAddress("jmva_2", &jmva_2, &b_jlrm_2);
   fChain->SetBranchAddress("jlrm_2", &jlrm_2, &b_jlrm_2);
   fChain->SetBranchAddress("jctm_2", &jctm_2, &b_jctm_2);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("jdeta", &jdeta, &b_jdeta);
   fChain->SetBranchAddress("njetingap", &njetingap, &b_njetingap);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("bpt", &bpt, &b_bpt);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("bphi", &bphi, &b_bphi);

   lock=true;
}

void Phys14Tree::ReadReset(){
  ientry = -2;
}

Long64_t Phys14Tree::GetEntries(){
  if (!fChain) return -1;
  
  return fChain->GetEntriesFast();
}

Long64_t Phys14Tree::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent)
    fCurrent = fChain->GetTreeNumber();

  ientry = entry;
  return centry;
}

Long64_t Phys14Tree::GetEntry(Long64_t entry){
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

Long64_t Phys14Tree::LoadedEntryId(){
  return ientry;
}
    

void Phys14Tree::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t Phys14Tree::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


//Write Methods
void Phys14Tree::WriteInit(TTree *tree) {
  if (!tree) return;
   
  fChain = tree;
  fCurrent = -1;

  fChain->Branch("run", &run, "run/I");
  fChain->Branch("lumi", &lumi, "lumi/I");
  fChain->Branch("evt", &evt, "evt/I");
  fChain->Branch("npv", &npv, "npv/I");
  fChain->Branch("npu", &npu, "npu/I");
  fChain->Branch("rho", &rho, "rho/F");
  
  fChain->Branch("mcweight", &mcweight, "mcweight/F");
  fChain->Branch("puweight", &puweight, "puweight/F");
  fChain->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
  fChain->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  fChain->Branch("idweight_1", &idweight_1, "idweight_1/F");
  fChain->Branch("idweight_2", &idweight_2, "idweight_2/F");
  fChain->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
  fChain->Branch("isoweight_2", &isoweight_2, "isoweight_2/F");
  fChain->Branch("effweight", &effweight, "effweight/F");
  fChain->Branch("fakeweight", &fakeweight, "fakeweight/F");
  fChain->Branch("embeddedWeight", &embeddedWeight, "embeddedWeight/F");
  fChain->Branch("signalWeight", &signalWeight, "signalWeight/F");
  fChain->Branch("weight", &weight, "weight/F");

  fChain->Branch("m_vis", &m_vis, "m_vis/F");
  fChain->Branch("m_sv", &m_sv, "m_sv/F");
  fChain->Branch("pt_sv", &pt_sv, "pt_sv/F");
  fChain->Branch("eta_sv", &eta_sv, "eta_sv/F");
  fChain->Branch("phi_sv", &phi_sv, "phi_sv/F");

  fChain->Branch("pt_1", &pt_1, "pt_1/F");
  fChain->Branch("phi_1", &phi_1, "phi_1/F");
  fChain->Branch("eta_1", &eta_1, "eta_1/F");
  fChain->Branch("m_1", &m_1, "m_1/F");
  fChain->Branch("q_1", &q_1, "q_1/I");
  fChain->Branch("iso_1", &iso_1, "iso_1/F");
  fChain->Branch("mva_1", &mva_1, "mva_1/F");
  fChain->Branch("d0_1", &d0_1, "d0_1/F");
  fChain->Branch("dZ_1", &dZ_1, "dZ_1/F");
  fChain->Branch("mt_1", &mt_1, "mt_1/F");

  fChain->Branch("pt_2", &pt_2, "pt_2/F");
  fChain->Branch("phi_2", &phi_2, "phi_2/F");
  fChain->Branch("eta_2", &eta_2, "eta_2/F");
  fChain->Branch("m_2", &m_2, "m_2/F");
  fChain->Branch("q_2", &q_2, "q_2/I");
  fChain->Branch("iso_2", &iso_2, "iso_2/F");
  fChain->Branch("d0_2", &d0_2, "d0_2/F");
  fChain->Branch("dZ_2", &dZ_2, "dZ_2/F");
  fChain->Branch("mva_2", &mva_2, "mva_2/F");
  fChain->Branch("mt_2", &mt_2, "mt_2/F");
  
  fChain->Branch("os", &os, "os/C");
  fChain->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/C");
  fChain->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/C");
  fChain->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/C");
  
  fChain->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, "byCombinedIsolationDeltaBetaCorrRaw3Hits_1/F");
  fChain->Branch("againstElectronLooseMVA5_1", &againstElectronLooseMVA5_1, "againstElectronLooseMVA5_1/F");
  fChain->Branch("againstElectronMediumMVA5_1", &againstElectronMediumMVA5_1, "againstElectronMediumMVA5_1/F");
  fChain->Branch("againstElectronTightMVA5_1", &againstElectronTightMVA5_1, "againstElectronTightMVA5_1/F");
  fChain->Branch("againstElectronVLooseMVA5_1", &againstElectronVLooseMVA5_1, "againstElectronVLooseMVA5_1/F");
  fChain->Branch("againstElectronVTightMVA5_1", &againstElectronVTightMVA5_1, "againstElectronVTightMVA5_1/F");
  fChain->Branch("againstMuonLoose3_1", &againstMuonLoose3_1, "againstMuonLoose3_1/F");
  fChain->Branch("againstMuonTight3_1", &againstMuonTight3_1, "againstMuonTight3_1/F");
  
  fChain->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
  fChain->Branch("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2, "againstElectronLooseMVA5_2/F");
  fChain->Branch("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2, "againstElectronMediumMVA5_2/F");
  fChain->Branch("againstElectronTightMVA5_2", &againstElectronTightMVA5_2, "againstElectronTightMVA5_2/F");
  fChain->Branch("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2, "againstElectronVLooseMVA5_2/F");
  fChain->Branch("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2, "againstElectronVTightMVA5_2/F");
  fChain->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/F");
  fChain->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");
    
  fChain->Branch("met", &met, "met/F");
  fChain->Branch("metphi", &metphi, "metphi/F");
  fChain->Branch("metcov00", &metcov00, "metcov00/F");
  fChain->Branch("metcov01", &metcov01, "metcov01/F");
  fChain->Branch("metcov10", &metcov10, "metcov10/F");
  fChain->Branch("metcov11", &metcov11, "metcov11/F");
  
  fChain->Branch("mvamet", &mvamet, "mvamet/F");
  fChain->Branch("mvametphi", &mvametphi, "mvametphi/F");
  fChain->Branch("mvacov00", &mvacov00, "mvacov00/F");
  fChain->Branch("mvacov01", &mvacov01, "mvacov01/F");
  fChain->Branch("mvacov10", &mvacov10, "mvacov10/F");
  fChain->Branch("mvacov11", &mvacov11, "mvacov11/F");
  
  fChain->Branch("pt_tt", &pt_tt, "pt_tt/F");
  fChain->Branch("pzetavis", &pzetavis, "pzetavis/F");
  fChain->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
  fChain->Branch("mva_gf", &mva_gf, "mva_gf/F");
  
  fChain->Branch("njets", &njets, "njets/I");
  fChain->Branch("njetspt20", &njetspt20, "njetspt20/I");
  
  fChain->Branch("jpt_1", &jpt_1, "jpt_1/F");
  fChain->Branch("jeta_1", &jeta_1, "jeta_1/F");
  fChain->Branch("jphi_1", &jphi_1, "jphi_1/F");
  fChain->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
  fChain->Branch("jptunc_1", &jptunc_1, "jptunc_1/F");
  fChain->Branch("jmva_1", &jmva_1, "jmva_1/F");
  fChain->Branch("jlrm_1", &jlrm_1, "jlrm_1/F");
  fChain->Branch("jctm_1", &jctm_1, "jctm_1/I");
  
  fChain->Branch("jpt_2", &jpt_2, "jpt_2/F");
  fChain->Branch("jeta_2", &jeta_2, "jeta_2/F");
  fChain->Branch("jphi_2", &jphi_2, "jphi_2/F");
  fChain->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
  fChain->Branch("jptunc_2", &jptunc_2, "jptunc_2/F");
  fChain->Branch("jmva_2", &jmva_2, "jlrm_2/F");
  fChain->Branch("jlrm_2", &jlrm_2, "jlrm_2/F");
  fChain->Branch("jctm_2", &jctm_2, "jctm_2/I");
  
  fChain->Branch("mjj", &mjj, "mjj/F");
  fChain->Branch("jdeta", &jdeta, "jdeta/F");
  fChain->Branch("njetingap", &njetingap, "njetingap/I");
    
  fChain->Branch("nbtag", &nbtag, "nbtag/I");
  fChain->Branch("bpt", &bpt, "bpt/F");
  fChain->Branch("beta", &beta, "beta/F");
  fChain->Branch("bphi", &bphi, "bphi/F");    
}

void Phys14Tree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}
