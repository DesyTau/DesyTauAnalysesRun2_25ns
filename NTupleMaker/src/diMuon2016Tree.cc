#define diMuon2016Tree_cxx
#include "DesyTauAnalyses/NTupleMaker/interface/diMuon2016Tree.h"

//Inizialization
diMuon2016Tree::diMuon2016Tree(TTree *tree) : fChain(0) {
  lock = false;
  ientry = -2;
  
  if (!tree) return;
  
  Init(tree);
}

void diMuon2016Tree::Init(TTree *tree){
  if(!tree) return;
	      
  if(tree->GetEntries())
    ReadInit(tree);
  else
    WriteInit(tree);
}

// Destructors
diMuon2016Tree::~diMuon2016Tree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// Read methods
void diMuon2016Tree::ReadInit(TTree *tree)
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
   
   fChain->SetBranchAddress("pt_0", &pt_0, &b_pt_0); 
   fChain->SetBranchAddress("eta_0", &eta_0, &b_eta_0); 
   fChain->SetBranchAddress("phi_0", &phi_0, &b_phi_0); 

   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1); 
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1); 
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1); 

    fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis); 

   lock=true;

}

void diMuon2016Tree::ReadReset(){
  ientry = -2;
}

Long64_t diMuon2016Tree::GetEntries(){
  if (!fChain) return -1;
  
  return fChain->GetEntriesFast();
}

Long64_t diMuon2016Tree::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  
  if (fChain->GetTreeNumber() != fCurrent)
    fCurrent = fChain->GetTreeNumber();

  ientry = entry;
  return centry;
}

Long64_t diMuon2016Tree::GetEntry(Long64_t entry){
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

Long64_t diMuon2016Tree::LoadedEntryId(){
  return ientry;
}

void diMuon2016Tree::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t diMuon2016Tree::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void diMuon2016Tree::WriteInit(TTree *tree) {
  if (!tree) return;
   
  fChain = tree;
  fCurrent = -1;

  fChain->Branch("run", &run, "run/I");
  fChain->Branch("lumi", &lumi, "lumi/I");
  fChain->Branch("evt", &evt, "evt/I");

  fChain->Branch("npv", &npv, "npv/I");

  fChain->Branch("pt_0", &pt_0, "pt_0/F");
  fChain->Branch("eta_0", &eta_0, "eta_0/F");
  fChain->Branch("phi_0", &phi_0, "phi_0/F");

  fChain->Branch("pt_1", &pt_1, "pt_1/F");
  fChain->Branch("eta_1", &eta_1, "eta_1/F");
  fChain->Branch("phi_1", &phi_1, "phi_1/F");

  fChain->Branch("m_vis", &m_vis, "m_vis/F");

}

void diMuon2016Tree::Fill(){
  if(!fChain) return;
  if(lock) return;

  fChain->Fill();
}
