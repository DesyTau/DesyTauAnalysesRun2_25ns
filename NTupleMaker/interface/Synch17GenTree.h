/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza
//////////////////////////////////////////////////////////

#ifndef Synch17GenTree_h
#define Synch17GenTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Synch17GenTree {
public :
  
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain
  
  // Declaration of leaf types
  //Event ID
  //Leptons
  Float_t         Higgs_pt;
  Float_t         Higgs_eta;
  Float_t         Higgs_phi;
  Float_t         Higgs_mass;

  Float_t         genpt_1;
  Float_t         genphi_1;
  Float_t         geneta_1;
  Int_t           genmode_1;
  Float_t         genpt_2;
  Float_t         genphi_2;
  Float_t         geneta_2;
  Int_t           genmode_2;
  Float_t         acotautau_00;
  Float_t         acotautau_10;
  Float_t         acotautau_01;
  Float_t         acotautau_11;
  
  

  //////////////////////////////////////////////
  //            List of branches              //
  //////////////////////////////////////////////
  TBranch        *b_Higgs_pt;
  TBranch        *b_Higgs_eta;
  TBranch        *b_Higgs_phi;
  TBranch        *b_Higgs_mass;
  TBranch	 *b_genpt_1;
  TBranch	 *b_genphi_1;
  TBranch	 *b_geneta_1;
  TBranch	 *b_genmode_1;
  TBranch	 *b_genpt_2;
  TBranch	 *b_genphi_2;
  TBranch	 *b_geneta_2;
  TBranch	 *b_genmode_2;
  TBranch        *b_acotautau_00;
  TBranch        *b_acotautau_10;
  TBranch        *b_acotautau_01;
  TBranch        *b_acotautau_11;
  
  Synch17GenTree(TTree *tree=0);
  virtual ~Synch17GenTree();

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

#endif //!endif Synch17GenTree_h
