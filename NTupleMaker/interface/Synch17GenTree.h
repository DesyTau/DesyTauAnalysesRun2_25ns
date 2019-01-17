/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza

//Merijn added some changes <merijn.van.de.klundert@desy.de>
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
  Float_t         acotautau_02;
  Float_t         acotautau_20;
  Float_t         acotautau_12;
  Float_t         acotautau_21;
  Float_t         acotautau_22;

  Float_t         acotautauPsi_00;
  Float_t         acotautauPsi_10;
  Float_t         acotautauPsi_01;
  Float_t         acotautauPsi_11;
  Float_t         acotautauPsi_02;
  Float_t         acotautauPsi_20;
  Float_t         acotautauPsi_12;
  Float_t         acotautauPsi_21;
  Float_t         acotautauPsi_22;


//gen vertex info is practical to have
  Float_t GenVertexX;
  Float_t GenVertexY;
  Float_t GenVertexZ;

    Float_t VxConstitTau1;
    Float_t VyConstitTau1;
    Float_t VzConstitTau1;

    Float_t VxConstitTau2;
    Float_t VyConstitTau2;
    Float_t VzConstitTau2;
  

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
  TBranch        *b_acotautau_02;
  TBranch        *b_acotautau_20;
  TBranch        *b_acotautau_12;
  TBranch        *b_acotautau_21;
  TBranch        *b_acotautau_22;

  TBranch        *b_acotautauPsi_00;
  TBranch        *b_acotautauPsi_10;
  TBranch        *b_acotautauPsi_01;
  TBranch        *b_acotautauPsi_11;
  TBranch        *b_acotautauPsi_02;
  TBranch        *b_acotautauPsi_20;
  TBranch        *b_acotautauPsi_12;
  TBranch        *b_acotautauPsi_21;
  TBranch        *b_acotautauPsi_22;

//gen vertex info is practical to have
  TBranch        *b_GenVertexX;
  TBranch        *b_GenVertexY;
  TBranch        *b_GenVertexZ;

  //Merijn: branches for the tau decay product constituent vx etc.
  TBranch        *b_VxConstitTau1;
  TBranch        *b_VyConstitTau1;
  TBranch        *b_VzConstitTau1;
  TBranch        *b_VxConstitTau2;
  TBranch        *b_VyConstitTau2;
  TBranch        *b_VzConstitTau2;

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
