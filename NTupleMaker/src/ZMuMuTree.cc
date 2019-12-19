#define ZMuMuTree_cxx
#include "DesyTauAnalyses/NTupleMaker/interface/ZMuMuTree.h"

//Inizialization
ZMuMuTree::ZMuMuTree(TTree *tree) : fChain(0) {
  lock = false;
  ientry = -2;
  
  if (!tree) return;
  
  Init(tree);
}

void ZMuMuTree::Init(TTree *tree){
  if(!tree) return;
	      
  if(tree->GetEntries())
    ReadInit(tree);
  else
    WriteInit(tree);
}

// Destructors
ZMuMuTree::~ZMuMuTree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// Read methods
void ZMuMuTree::ReadInit(TTree *tree)
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

   // leading and trailing muons indicated with _1 and _2, respectively
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1); 
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1); 
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1); 
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);

   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2); 
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2); 
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2); 
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);

   fChain->SetBranchAddress("os", &os, &b_os); 
   fChain->SetBranchAddress("pt_ll", &pt_ll, &b_pt_ll); 
   fChain->SetBranchAddress("px_ll", &px_ll, &b_px_ll); 
   fChain->SetBranchAddress("py_ll", &py_ll, &b_py_ll); 
   fChain->SetBranchAddress("m_ll", &m_ll, &b_m_ll); 
   fChain->SetBranchAddress("eta_ll", &eta_ll, &b_eta_ll); 
   fChain->SetBranchAddress("phi_ll", &phi_ll, &b_phi_ll); 
   fChain->SetBranchAddress("delta_phi_ll", &delta_phi_ll, &b_delta_phi_ll); 

   fChain->SetBranchAddress("genV_px", &genV_px, &b_genV_px); 
   fChain->SetBranchAddress("genV_py", &genV_py, &b_genV_py); 
   fChain->SetBranchAddress("genV_pt", &genV_pt, &b_genV_pt); 
   fChain->SetBranchAddress("genV_mass", &genV_mass, &b_genV_mass); 
   fChain->SetBranchAddress("genV_eta", &genV_eta, &b_genV_eta); 
   fChain->SetBranchAddress("genV_phi", &genV_phi, &b_genV_phi); 
   fChain->SetBranchAddress("genZ_pt", &genZ_pt, &b_genZ_pt); 
   fChain->SetBranchAddress("genZ_mass", &genZ_mass, &b_genZ_mass); 
   fChain->SetBranchAddress("is_gen_Zee", &is_gen_Zee, &b_is_gen_Zee); 
   fChain->SetBranchAddress("is_gen_Zmm", &is_gen_Zmm, &b_is_gen_Zmm);  
   fChain->SetBranchAddress("is_gen_Ztt", &is_gen_Ztt, &b_is_gen_Ztt); 

   // leading and trailing jets indicated with _1 and _2, respectively
   fChain->SetBranchAddress("pt_jet_1", &pt_jet_1, &b_pt_jet_1); 
   fChain->SetBranchAddress("eta_jet_1", &eta_jet_1, &b_eta_jet_1); 
   fChain->SetBranchAddress("phi_jet_1", &phi_jet_1, &b_phi_jet_1); 

   fChain->SetBranchAddress("pt_jet_2", &pt_jet_2, &b_pt_jet_2); 
   fChain->SetBranchAddress("eta_jet_2", &eta_jet_2, &b_eta_jet_2); 
   fChain->SetBranchAddress("phi_jet_2", &phi_jet_2, &b_phi_jet_2); 

   fChain->SetBranchAddress("m_jj", &m_jj, &b_m_jj); 
   fChain->SetBranchAddress("delta_phi_jj", &delta_phi_jj, &b_delta_phi_jj); 
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30); 
   fChain->SetBranchAddress("njets20", &njets20, &b_njets20); 
   fChain->SetBranchAddress("HT30", &HT30, &b_HT30); 
   fChain->SetBranchAddress("HT20", &HT20, &b_HT20); 

   fChain->SetBranchAddress("met", &met, &b_met); 
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi); 
   fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet); 
   fChain->SetBranchAddress("puppimet_phi", &puppimet_phi, &b_puppimet_phi); 

   fChain->SetBranchAddress("met_rcorr", &met_rcorr, &b_met_rcorr);
   fChain->SetBranchAddress("met_phi_rcorr", &met_phi_rcorr, &b_met_phi_rcorr);
   fChain->SetBranchAddress("puppimet_rcorr", &puppimet_rcorr, &b_puppimet_rcorr);
   fChain->SetBranchAddress("puppimet_phi_rcorr", &puppimet_phi_rcorr, &b_puppimet_phi_rcorr);

   fChain->SetBranchAddress("pt_ratio_ll", &pt_ratio_ll, &b_pt_ratio_ll);

   fChain->SetBranchAddress("met_response_gen", &met_response_gen, &b_met_response_gen);
   fChain->SetBranchAddress("met_recoil_paral", &met_recoil_paral, &b_met_recoil_paral);
   fChain->SetBranchAddress("met_recoil_perp", &met_recoil_perp, &b_met_recoil_perp);
   fChain->SetBranchAddress("met_response_had", &met_response_had, &b_met_response_had);
   fChain->SetBranchAddress("puppimet_response_gen", &puppimet_response_gen, &b_puppimet_response_gen);
   fChain->SetBranchAddress("puppimet_recoil_paral", &puppimet_recoil_paral, &b_puppimet_recoil_paral);
   fChain->SetBranchAddress("puppimet_recoil_perp", &puppimet_recoil_perp, &b_puppimet_recoil_perp);
   fChain->SetBranchAddress("puppimet_response_had", &puppimet_response_had, &b_puppimet_response_had);

   fChain->SetBranchAddress("met_response_gen_rcorr", &met_response_gen_rcorr, &b_met_response_gen_rcorr);
   fChain->SetBranchAddress("met_recoil_paral_rcorr", &met_recoil_paral_rcorr, &b_met_recoil_paral_rcorr);
   fChain->SetBranchAddress("met_recoil_perp_rcorr", &met_recoil_perp_rcorr, &b_met_recoil_perp_rcorr);
   fChain->SetBranchAddress("met_response_had_rcorr", &met_response_had_rcorr, &b_met_response_had_rcorr);
   fChain->SetBranchAddress("puppimet_response_gen_rcorr", &puppimet_response_gen_rcorr, &b_puppimet_response_gen_rcorr);
   fChain->SetBranchAddress("puppimet_recoil_paral_rcorr", &puppimet_recoil_paral_rcorr, &b_puppimet_recoil_paral_rcorr);
   fChain->SetBranchAddress("puppimet_recoil_perp_rcorr", &puppimet_recoil_perp_rcorr, &b_puppimet_recoil_perp_rcorr);
   fChain->SetBranchAddress("puppimet_response_had_rcorr", &puppimet_response_had_rcorr, &b_puppimet_response_had_rcorr);

   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight); 
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight); 
   fChain->SetBranchAddress("toppt_weight", &toppt_weight, &b_toppt_weight); 
   fChain->SetBranchAddress("zpt_weight", &zpt_weight, &b_zpt_weight); 
   fChain->SetBranchAddress("track_weight_1", &track_weight_1, &b_track_weight_1); 
   fChain->SetBranchAddress("track_weight_2", &track_weight_2, &b_track_weight_2); 
   fChain->SetBranchAddress("idiso_weight_1", &idiso_weight_1, &b_idiso_weight_1); 
   fChain->SetBranchAddress("idiso_weight_2", &idiso_weight_2, &b_idiso_weight_2); 
   fChain->SetBranchAddress("trig_weight", &trig_weight, &b_trig_weight); 

   lock=true;

}

void ZMuMuTree::ReadReset(){
  ientry = -2;
}

Long64_t ZMuMuTree::GetEntries(){
  if (!fChain) return -1;
  
  return fChain->GetEntriesFast();
}

Long64_t ZMuMuTree::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent)
    fCurrent = fChain->GetTreeNumber();

  ientry = entry;
  return centry;
}

Long64_t ZMuMuTree::GetEntry(Long64_t entry){
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

Long64_t ZMuMuTree::LoadedEntryId(){
  return ientry;
}

void ZMuMuTree::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t ZMuMuTree::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void ZMuMuTree::WriteInit(TTree *tree) {
  if (!tree) return;
   
  fChain = tree;
  fCurrent = -1;

  fChain->Branch("run", &run, "run/l");
  fChain->Branch("lumi", &lumi, "lumi/I");
  fChain->Branch("evt", &evt, "evt/I");
  fChain->Branch("npv", &npv, "npv/I");

  fChain->Branch("pt_1", &pt_1, "pt_1/F");
  fChain->Branch("eta_1", &eta_1, "eta_1/F");
  fChain->Branch("phi_1", &phi_1, "phi_1/F");
  fChain->Branch("iso_1", &iso_1, "iso_1/F");
  fChain->Branch("q_1", &q_1, "q_1/F");

  fChain->Branch("pt_2", &pt_2, "pt_2/F");
  fChain->Branch("eta_2", &eta_2, "eta_2/F");
  fChain->Branch("phi_2", &phi_2, "phi_2/F");
  fChain->Branch("iso_2", &iso_2, "iso_2/F");
  fChain->Branch("q_2", &q_2, "q_2/F");

  fChain->Branch("delta_phi_ll", &delta_phi_ll, "delta_phi_ll/F");
  fChain->Branch("m_ll", &m_ll, "m_ll/F");
  fChain->Branch("pt_ll", &pt_ll, "pt_ll/F");
  fChain->Branch("px_ll", &px_ll, "px_ll/F");
  fChain->Branch("py_ll", &py_ll, "py_ll/F");
  fChain->Branch("eta_ll", &eta_ll, "eta_ll/F");
  fChain->Branch("phi_ll", &phi_ll, "phi_ll/F");
  fChain->Branch("os", &os, "os/F");

  fChain->Branch("genV_px", &genV_px, "genV_px/F");
  fChain->Branch("genV_py", &genV_py, "genV_py/F");
  fChain->Branch("genV_pt", &genV_pt, "genV_pt/F");
  fChain->Branch("genV_mass", &genV_mass, "genV_mass/F");
  fChain->Branch("genV_eta", &genV_eta, "genV_eta/F");
  fChain->Branch("genV_phi", &genV_phi, "genV_phi/F");
  fChain->Branch("genZ_pt", &genZ_pt, "genZ_pt/F");
  fChain->Branch("genZ_mass", &genZ_mass, "genZ_mass/F");
  fChain->Branch("is_gen_Zee", &is_gen_Zee, "is_gen_Zee/F");
  fChain->Branch("is_gen_Zmm", &is_gen_Zmm, "is_gen_Zmm/F");
  fChain->Branch("is_gen_Ztt", &is_gen_Ztt, "is_gen_Ztt/F");

  fChain->Branch("pt_jet_1", &pt_jet_1, "pt_jet_1/F");
  fChain->Branch("eta_jet_1", &eta_jet_1, "eta_jet_1/F");
  fChain->Branch("phi_jet_1", &phi_jet_1, "phi_jet_1/F");

  fChain->Branch("pt_jet_2", &pt_jet_2, "pt_jet_2/F");
  fChain->Branch("eta_jet_2", &eta_jet_2, "eta_jet_2/F");
  fChain->Branch("phi_jet_2", &phi_jet_2, "phi_jet_2/F");

  fChain->Branch("m_jj", &m_jj, "m_jj/F");
  fChain->Branch("delta_phi_jj", &delta_phi_jj, "delta_phi_jj/F");
  fChain->Branch("njets30", &njets30, "njets30/I");
  fChain->Branch("njets20", &njets20, "njets20/I");
  fChain->Branch("HT30", &HT30, "HT30/F");
  fChain->Branch("HT20", &HT20, "HT20/F");

  fChain->Branch("met", &met, "met/F");
  fChain->Branch("met_phi", &met_phi, "met_phi/F");
  fChain->Branch("puppimet", &puppimet, "puppimet/F");
  fChain->Branch("puppimet_phi", &puppimet_phi, "puppimet_phi/F");
  fChain->Branch("met_rcorr", &met_rcorr, "met_rcorr/F");
  fChain->Branch("met_phi_rcorr", &met_phi_rcorr, "met_phi_rcorr/F");
  fChain->Branch("puppimet_rcorr", &puppimet_rcorr, "puppimet_rcorr/F");
  fChain->Branch("puppimet_phi_rcorr", &puppimet_phi_rcorr, "puppimet_phi_rcorr/F");
  fChain->Branch("pt_ratio_ll", &pt_ratio_ll, "pt_ratio_ll/F");
  fChain->Branch("met_response_gen", &met_response_gen, "met_response_gen/F");
  fChain->Branch("puppimet_response_gen", &puppimet_response_gen, "puppimet_response_gen/F");
  fChain->Branch("met_recoil_paral", &met_recoil_paral, "met_recoil_paral/F");
  fChain->Branch("met_recoil_perp", &met_recoil_perp, "met_recoil_perp/F");
  fChain->Branch("met_response_had", &met_response_had, "met_response_had/F");
  fChain->Branch("puppimet_recoil_paral", &puppimet_recoil_paral, "puppimet_recoil_paral/F");
  fChain->Branch("puppimet_recoil_perp", &puppimet_recoil_perp, "puppimet_recoil_perp/F");
  fChain->Branch("puppimet_response_had", &puppimet_response_had, "puppimet_response_had/F");

  fChain->Branch("met_response_gen_rcorr", &met_response_gen_rcorr, "met_response_gen_rcorr/F");
  fChain->Branch("puppimet_response_gen_rcorr", &puppimet_response_gen_rcorr, "puppimet_response_gen_rcorr/F");
  fChain->Branch("met_recoil_paral_rcorr", &met_recoil_paral_rcorr, "met_recoil_paral_rcorr/F");
  fChain->Branch("met_recoil_perp_rcorr", &met_recoil_perp_rcorr, "met_recoil_perp_rcorr/F");
  fChain->Branch("met_response_had_rcorr", &met_response_had_rcorr, "met_response_had_rcorr/F");
  fChain->Branch("puppimet_recoil_paral_rcorr", &puppimet_recoil_paral_rcorr, "puppimet_recoil_paral_rcorr/F");
  fChain->Branch("puppimet_recoil_perp_rcorr", &puppimet_recoil_perp_rcorr, "puppimet_recoil_perp_rcorr/F");
  fChain->Branch("puppimet_response_had_rcorr", &puppimet_response_had_rcorr, "puppimet_response_had_rcorr/F");

  fChain->Branch("mc_weight", &mc_weight, "mc_weight/F");
  fChain->Branch("pu_weight", &pu_weight, "pu_weight/F");
  fChain->Branch("toppt_weight", &toppt_weight, "toppt_weight/D");
  fChain->Branch("zpt_weight", &zpt_weight, "zpt_weight/D");
  fChain->Branch("track_weight_1", &track_weight_1, "track_weight_1/D");
  fChain->Branch("track_weight_2", &track_weight_2, "track_weight_2/D");
  fChain->Branch("idiso_weight_1", &idiso_weight_1, "idiso_weight_1/D");
  fChain->Branch("idiso_weight_2", &idiso_weight_2, "idiso_weight_2/D");
  fChain->Branch("trig_weight", &trig_weight, "trig_weight/D");
}

void ZMuMuTree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}
