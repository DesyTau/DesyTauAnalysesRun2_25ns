/////////////////////////////////////////////////////////////
// Read and Write Synch Ntuple for CP measurement in h->tau tau
// Author: Andrea Cardini <andrea.cardini@desy.de>
// 
// Based on Spring15Tree by Francesco Costanza

//Merijn added changes and tauspinner weights.. <merijn.van.de.klundert@desy.de>
//////////////////////////////////////////////////////////

#ifndef Synch17GenTree_h
#define Synch17GenTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//We´ll store for now 5 CP mixing scenarios: "sm_htt125", "ps_htt125", "mm_htt125" "minusmm_htt125", "mix0p375_htt125". Names chose to stay somewhat consistent with choices IC in the CH branch
#include <vector>
#include <string>

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

  Float_t         pt_1;
  Float_t         phi_1;
  Float_t         eta_1;
  Int_t           decaymode_1;
  Float_t         chconst_1_pt;
  Float_t         chconst_1_eta;
  Float_t         chconst_1_phi; 
  Float_t         pt_2;
  Float_t         phi_2;
  Float_t         eta_2;
  Float_t         chconst_2_pt;
  Float_t         chconst_2_eta;
  Float_t         chconst_2_phi;   
  Int_t           decaymode_2;
  Float_t         acotautau_00;
  Float_t         acotautau_10;
  Float_t         acotautau_01;
  Float_t         acotautau_11;
  Float_t         acotautau_02;
  Float_t         acotautau_20;
  Float_t         acotautau_12;
  Float_t         acotautau_21;
  Float_t         acotautau_22;

  Float_t         a1polarization_1;
  Float_t         a1polarization_2;

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
  Float_t VertexX;
  Float_t VertexY;
  Float_t VertexZ;
  
  Float_t VxConstitTau1;
  Float_t VyConstitTau1;
  Float_t VzConstitTau1;
  
  Float_t VxConstitTau2;
  Float_t VyConstitTau2;
  Float_t VzConstitTau2;
  
  Float_t alphaminus;
  //tauspinner weights
  Double_t sm_htt125;
  Double_t ps_htt125;
  Double_t mm_htt125;
  Double_t minusmm_htt125;
  Double_t mix0p375_htt125;

  //////////////////////////////////////////////
  //            List of branches              //
  //////////////////////////////////////////////
  TBranch        *b_Higgs_pt;
  TBranch        *b_Higgs_eta;
  TBranch        *b_Higgs_phi;
  TBranch        *b_Higgs_mass;
  TBranch	 *b_pt_1;
  TBranch	 *b_phi_1;
  TBranch	 *b_eta_1;
  TBranch        *b_chconst_1_pt;
  TBranch        *b_chconst_1_eta;
  TBranch        *b_chconst_1_phi;   
  TBranch	 *b_decaymode_1;
  TBranch	 *b_pt_2;
  TBranch	 *b_phi_2;
  TBranch	 *b_eta_2;
  TBranch        *b_chconst_2_pt;
  TBranch        *b_chconst_2_eta;
  TBranch        *b_chconst_2_phi;    
  TBranch	 *b_decaymode_2;
  TBranch        *b_acotautau_00;
  TBranch        *b_acotautau_10;
  TBranch        *b_acotautau_01;
  TBranch        *b_acotautau_11;
  TBranch        *b_acotautau_02;
  TBranch        *b_acotautau_20;
  TBranch        *b_acotautau_12;
  TBranch        *b_acotautau_21;
  TBranch        *b_acotautau_22;

  TBranch        *b_a1polarization_1;
  TBranch        *b_a1polarization_2;

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
  TBranch        *b_VertexX;
  TBranch        *b_VertexY;
  TBranch        *b_VertexZ;

  //Merijn: branches for the tau decay product constituent vx etc.
  TBranch        *b_VxConstitTau1;
  TBranch        *b_VyConstitTau1;
  TBranch        *b_VzConstitTau1;
  TBranch        *b_VxConstitTau2;
  TBranch        *b_VyConstitTau2;
  TBranch        *b_VzConstitTau2;
  TBranch        *b_alphaminus;

  TBranch        *b_sm_htt125;
  TBranch        *b_ps_htt125;
  TBranch        *b_mm_htt125;
  TBranch        *b_minusmm_htt125;
  TBranch        *b_mix0p375_htt125;

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
