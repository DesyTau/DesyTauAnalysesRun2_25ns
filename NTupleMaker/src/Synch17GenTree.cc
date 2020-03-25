/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza

//Merijn van de Klundert <merijn.van.de.klundert@desy.de>

//////////////////////////////////////////////////////////


#define Synch17Tree_cxx
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"

// Initialization
Synch17GenTree::Synch17GenTree(TTree *tree) : fChain(0) {
  lock = false;
  ientry = -2;
  
  if (!tree) return;
  
  Init(tree);
}

void Synch17GenTree::Init(TTree *tree){
  if(!tree) return;
	      
  if(tree->GetEntries())
    ReadInit(tree);
  else
    WriteInit(tree);
}

// Destructors
Synch17GenTree::~Synch17GenTree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// Read methods
void Synch17GenTree::ReadInit(TTree *tree)
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

   fChain->SetBranchAddress("Higgs_pt", &Higgs_pt, &b_Higgs_pt);
   fChain->SetBranchAddress("Higgs_eta", &Higgs_eta, &b_Higgs_eta);
   fChain->SetBranchAddress("Higgs_phi", &Higgs_phi, &b_Higgs_phi);
   fChain->SetBranchAddress("Higgs_mass", &Higgs_mass, &b_Higgs_mass);

   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);

   fChain->SetBranchAddress("chconst_1_pt", &chconst_1_pt, &b_chconst_1_pt);
   fChain->SetBranchAddress("chconst_1_eta", &chconst_1_eta, &b_chconst_1_eta);
   fChain->SetBranchAddress("chconst_1_phi", &chconst_1_phi, &b_chconst_1_phi); 
   
   fChain->SetBranchAddress("decaymode_1", &decaymode_1, &b_decaymode_1);

   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);

   fChain->SetBranchAddress("chconst_2_pt", &chconst_2_pt, &b_chconst_2_pt);
   fChain->SetBranchAddress("chconst_2_eta", &chconst_2_eta, &b_chconst_2_eta);
   fChain->SetBranchAddress("chconst_2_phi", &chconst_2_phi, &b_chconst_2_phi); 
   
   fChain->SetBranchAddress("decaymode_2", &decaymode_2, &b_decaymode_2);

   fChain->SetBranchAddress("acotautau_00", &acotautau_00, &b_acotautau_00);
   fChain->SetBranchAddress("acotautau_10", &acotautau_10, &b_acotautau_10);
   fChain->SetBranchAddress("acotautau_01", &acotautau_01, &b_acotautau_01);
   fChain->SetBranchAddress("acotautau_11", &acotautau_11, &b_acotautau_11);   
   fChain->SetBranchAddress("acotautau_02", &acotautau_02, &b_acotautau_02);
   fChain->SetBranchAddress("acotautau_12", &acotautau_12, &b_acotautau_12);
   fChain->SetBranchAddress("acotautau_20", &acotautau_20, &b_acotautau_20);
   fChain->SetBranchAddress("acotautau_21", &acotautau_21, &b_acotautau_21);
   fChain->SetBranchAddress("acotautau_22", &acotautau_22, &b_acotautau_22);

   fChain->SetBranchAddress("a1polarization_1", &a1polarization_1, &b_a1polarization_1);
   fChain->SetBranchAddress("a1polarization_2", &a1polarization_2, &b_a1polarization_2);

   fChain->SetBranchAddress("acotautauPsi_00", &acotautauPsi_00, &b_acotautauPsi_00);
   fChain->SetBranchAddress("acotautauPsi_10", &acotautauPsi_10, &b_acotautauPsi_10);
   fChain->SetBranchAddress("acotautauPsi_01", &acotautauPsi_01, &b_acotautauPsi_01);
   fChain->SetBranchAddress("acotautauPsi_11", &acotautauPsi_11, &b_acotautauPsi_11);   
   fChain->SetBranchAddress("acotautauPsi_02", &acotautauPsi_02, &b_acotautauPsi_02);
   fChain->SetBranchAddress("acotautauPsi_12", &acotautauPsi_12, &b_acotautauPsi_12);
   fChain->SetBranchAddress("acotautauPsi_20", &acotautauPsi_20, &b_acotautauPsi_20);
   fChain->SetBranchAddress("acotautauPsi_21", &acotautauPsi_21, &b_acotautauPsi_21);
   fChain->SetBranchAddress("acotautauPsi_22", &acotautauPsi_22, &b_acotautauPsi_22);

  //gen vertex info useful to have
   fChain->SetBranchAddress("VertexX", &VertexX, &b_VertexX);
   fChain->SetBranchAddress("VertexY", &VertexY, &b_VertexY);
   fChain->SetBranchAddress("VertexZ", &VertexZ, &b_VertexZ);

   //need to have vx tau constituents..
   fChain->SetBranchAddress("VxConstitTau1", &VxConstitTau1, &b_VxConstitTau1);
   fChain->SetBranchAddress("VyConstitTau1", &VyConstitTau1, &b_VyConstitTau1);
   fChain->SetBranchAddress("VzConstitTau1", &VzConstitTau1, &b_VzConstitTau1);
   
   fChain->SetBranchAddress("VxConstitTau2", &VxConstitTau2, &b_VxConstitTau2);
   fChain->SetBranchAddress("VyConstitTau2", &VyConstitTau2, &b_VyConstitTau2);
   fChain->SetBranchAddress("VzConstitTau2", &VzConstitTau2, &b_VzConstitTau2);

   fChain->SetBranchAddress("alphaminus", &alphaminus, &b_alphaminus);
   
   fChain->SetBranchAddress("sm_htt125", &sm_htt125, &b_sm_htt125);
   fChain->SetBranchAddress("ps_htt125", &ps_htt125, &b_ps_htt125);
   fChain->SetBranchAddress("mm_htt125", &mm_htt125, &b_mm_htt125);
   fChain->SetBranchAddress("minusmm_htt125", &minusmm_htt125, &b_minusmm_htt125);
   fChain->SetBranchAddress("mix0p375_htt125", &mix0p375_htt125, &b_mix0p375_htt125);

   fChain->SetBranchAddress("y1_LF", &y1_LF, &b_y1_LF);
   fChain->SetBranchAddress("y2_LF", &y2_LF, &b_y2_LF);
   fChain->SetBranchAddress("y1_ZMF", &y1_ZMF, &b_y1_ZMF);
   fChain->SetBranchAddress("y2_ZMF", &y1_ZMF, &b_y1_ZMF);

   lock=true;
}

void Synch17GenTree::ReadReset(){
  ientry = -2;
}

Long64_t Synch17GenTree::GetEntries(){
  if (!fChain) return -1;
  
  return fChain->GetEntriesFast();
}

Long64_t Synch17GenTree::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent)
    fCurrent = fChain->GetTreeNumber();

  ientry = entry;
  return centry;
}

Long64_t Synch17GenTree::GetEntry(Long64_t entry){
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

Long64_t Synch17GenTree::LoadedEntryId(){
  return ientry;
}
    

void Synch17GenTree::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t Synch17GenTree::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


//Write Methods
void Synch17GenTree::WriteInit(TTree *tree) {
  if (!tree) return;
   
  fChain = tree;
  fCurrent = -1;

   fChain->Branch("gen_Higgs_pt", &Higgs_pt, "gen_Higgs_pt/F");
   fChain->Branch("gen_Higgs_eta", &Higgs_eta, "gen_Higgs_eta/F");
   fChain->Branch("gen_Higgs_phi", &Higgs_phi, "gen_Higgs_phi/F");
   fChain->Branch("gen_Higgs_mass", &Higgs_mass, "gen_Higgs_mass/F");

   fChain->Branch("gen_pt_1", &pt_1, "gen_pt_1/F");
   fChain->Branch("gen_phi_1", &phi_1, "gen_phi_1/F");
   fChain->Branch("gen_eta_1", &eta_1, "gen_eta_1/F");
   fChain->Branch("gen_chconst_1_pt", &chconst_1_pt, "gen_chconst_1_pt/F");
   fChain->Branch("gen_chconst_1_phi", &chconst_1_phi, "gen_chconst_1_phi/F");
   fChain->Branch("gen_chconst_1_eta", &chconst_1_eta, "gen_chconst_1_eta/F");
   fChain->Branch("gen_decaymode_1", &decaymode_1, "gen_decaymode_1/I");

   fChain->Branch("gen_pt_2", &pt_2, "gen_pt_2/F");
   fChain->Branch("gen_phi_2", &phi_2, "gen_phi_2/F");
   fChain->Branch("gen_eta_2", &eta_2, "gen_eta_2/F");

   fChain->Branch("gen_chconst_2_pt", &chconst_2_pt, "gen_chconst_2_pt/F");
   fChain->Branch("gen_chconst_2_phi", &chconst_2_phi, "gen_chconst_2_phi/F");
   fChain->Branch("gen_chconst_2_eta", &chconst_2_eta, "gen_chconst_2_eta/F");
   
   fChain->Branch("gen_decaymode_2", &decaymode_2, "gen_decaymode_2/I");
   fChain->Branch("gen_a1polarization_1", &a1polarization_1, "a1polarization_1/F");
   fChain->Branch("gen_a1polarization_2", &a1polarization_2, "a1polarization_2/F");

   fChain->Branch("gen_acotautau_00", &acotautau_00, "gen_acotautau_00/F");
   fChain->Branch("gen_acotautau_10", &acotautau_10, "gen_acotautau_10/F");
   fChain->Branch("gen_acotautau_01", &acotautau_01, "gen_acotautau_01/F");
   fChain->Branch("gen_acotautau_11", &acotautau_11, "gen_acotautau_11/F");
   fChain->Branch("gen_acotautau_02", &acotautau_02, "gen_acotautau_02/F");
   fChain->Branch("gen_acotautau_12", &acotautau_12, "gen_acotautau_12/F");
   fChain->Branch("gen_acotautau_20", &acotautau_20, "gen_acotautau_20/F");
   fChain->Branch("gen_acotautau_21", &acotautau_21, "gen_acotautau_21/F");
   fChain->Branch("gen_acotautau_22", &acotautau_22, "gen_acotautau_22/F");

   fChain->Branch("gen_acotautauPsi_00", &acotautauPsi_00, "gen_acotautauPsi_00/F");
   fChain->Branch("gen_acotautauPsi_10", &acotautauPsi_10, "gen_acotautauPsi_10/F");
   fChain->Branch("gen_acotautauPsi_01", &acotautauPsi_01, "gen_acotautauPsi_01/F");
   fChain->Branch("gen_acotautauPsi_11", &acotautauPsi_11, "gen_acotautauPsi_11/F");
   fChain->Branch("gen_acotautauPsi_02", &acotautauPsi_02, "gen_acotautauPsi_02/F");
   fChain->Branch("gen_acotautauPsi_12", &acotautauPsi_12, "gen_acotautauPsi_12/F");
   fChain->Branch("gen_acotautauPsi_20", &acotautauPsi_20, "gen_acotautauPsi_20/F");
   fChain->Branch("gen_acotautauPsi_21", &acotautauPsi_21, "gen_acotautauPsi_21/F");
   fChain->Branch("gen_acotautauPsi_22", &acotautauPsi_22, "gen_acotautauPsi_22/F");

   fChain->Branch("gen_VertexX", &VertexX, "gen_VertexX/F");
   fChain->Branch("gen_VertexY", &VertexY, "gen_VertexY/F");
   fChain->Branch("gen_VertexZ", &VertexZ, "gen_VertexZ/F");

   fChain->Branch("gen_VxConstitTau1", &VxConstitTau1, "gen_VxConstitTau1/F");
   fChain->Branch("gen_VyConstitTau1", &VyConstitTau1, "gen_VyConstitTau1/F");
   fChain->Branch("gen_VzConstitTau1", &VzConstitTau1, "gen_VzConstitTau1/F");
   
   fChain->Branch("gen_VxConstitTau2", &VxConstitTau2, "gen_VxConstitTau2/F");
   fChain->Branch("gen_VyConstitTau2", &VyConstitTau2, "gen_VyConstitTau2/F");
   fChain->Branch("gen_VzConstitTau2", &VzConstitTau2, "gen_VzConstitTau2/F");

   fChain->Branch("gen_alphaminus", &alphaminus, "gen_alphaminus/F");

  fChain->Branch("gen_sm_htt125", &sm_htt125, "gen_sm_htt125/D");
  fChain->Branch("gen_ps_htt125", &ps_htt125, "gen_ps_htt125/D");
  fChain->Branch("gen_mm_htt125", &mm_htt125, "gen_mm_htt125/D");
  fChain->Branch("gen_minusmm_htt125", &minusmm_htt125, "gen_minusmm_htt125/D");
  fChain->Branch("gen_mix0p375_htt125", &mix0p375_htt125, "gen_mix0p375_htt125/D");

   fChain->Branch("gen_y1_LF", &y1_LF, "gen_y1_LF/D");
   fChain->Branch("gen_y2_LF", &y2_LF, "gen_y2_LF/D");
   fChain->Branch("gen_y1_ZMF", &y1_ZMF, "gen_y1_ZMF/D");
   fChain->Branch("gen_y2_ZMF", &y2_ZMF, "gen_y2_ZMF/D");

}

void Synch17GenTree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}
