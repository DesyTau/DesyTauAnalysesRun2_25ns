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
   fChain->SetBranchAddress("decaymode_1", &decaymode_1, &b_decaymode_1);

   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
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

   fChain->Branch("Higgs_pt", &Higgs_pt, "Higgs_pt/F");
   fChain->Branch("Higgs_eta", &Higgs_eta, "Higgs_eta/F");
   fChain->Branch("Higgs_phi", &Higgs_phi, "Higgs_phi/F");
   fChain->Branch("Higgs_mass", &Higgs_mass, "Higgs_mass/F");

   fChain->Branch("pt_1", &pt_1, "pt_1/F");
   fChain->Branch("phi_1", &phi_1, "phi_1/F");
   fChain->Branch("eta_1", &eta_1, "eta_1/F");
   fChain->Branch("decaymode_1", &decaymode_1, "decaymode_1/I");

   fChain->Branch("pt_2", &pt_2, "pt_2/F");
   fChain->Branch("phi_2", &phi_2, "phi_2/F");
   fChain->Branch("eta_2", &eta_2, "eta_2/F");
   fChain->Branch("decaymode_2", &decaymode_2, "decaymode_2/I");

   fChain->Branch("acotautau_00", &acotautau_00, "acotautau_00/F");
   fChain->Branch("acotautau_10", &acotautau_10, "acotautau_10/F");
   fChain->Branch("acotautau_01", &acotautau_01, "acotautau_01/F");
   fChain->Branch("acotautau_11", &acotautau_11, "acotautau_11/F");
   fChain->Branch("acotautau_02", &acotautau_02, "acotautau_02/F");
   fChain->Branch("acotautau_12", &acotautau_12, "acotautau_12/F");
   fChain->Branch("acotautau_20", &acotautau_20, "acotautau_20/F");
   fChain->Branch("acotautau_21", &acotautau_21, "acotautau_21/F");
   fChain->Branch("acotautau_22", &acotautau_22, "acotautau_22/F");

   fChain->Branch("acotautauPsi_00", &acotautauPsi_00, "acotautauPsi_00/F");
   fChain->Branch("acotautauPsi_10", &acotautauPsi_10, "acotautauPsi_10/F");
   fChain->Branch("acotautauPsi_01", &acotautauPsi_01, "acotautauPsi_01/F");
   fChain->Branch("acotautauPsi_11", &acotautauPsi_11, "acotautauPsi_11/F");
   fChain->Branch("acotautauPsi_02", &acotautauPsi_02, "acotautauPsi_02/F");
   fChain->Branch("acotautauPsi_12", &acotautauPsi_12, "acotautauPsi_12/F");
   fChain->Branch("acotautauPsi_20", &acotautauPsi_20, "acotautauPsi_20/F");
   fChain->Branch("acotautauPsi_21", &acotautauPsi_21, "acotautauPsi_21/F");
   fChain->Branch("acotautauPsi_22", &acotautauPsi_22, "acotautauPsi_22/F");

   fChain->Branch("VertexX", &VertexX, "VertexX/F");
   fChain->Branch("VertexY", &VertexY, "VertexY/F");
   fChain->Branch("VertexZ", &VertexZ, "VertexZ/F");

   fChain->Branch("VxConstitTau1", &VxConstitTau1, "VxConstitTau1/F");
   fChain->Branch("VyConstitTau1", &VyConstitTau1, "VyConstitTau1/F");
   fChain->Branch("VzConstitTau1", &VzConstitTau1, "VzConstitTau1/F");
   
   fChain->Branch("VxConstitTau2", &VxConstitTau2, "VxConstitTau2/F");
   fChain->Branch("VyConstitTau2", &VyConstitTau2, "VyConstitTau2/F");
   fChain->Branch("VzConstitTau2", &VzConstitTau2, "VzConstitTau2/F");

}

void Synch17GenTree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}

