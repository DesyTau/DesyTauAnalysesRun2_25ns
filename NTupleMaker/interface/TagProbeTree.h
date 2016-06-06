#ifndef TagProbeTree_h
#define TagProbeTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TagProbeTree {
public :
	
	TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain


  // Declaration of leaf types
  Int_t           run;
  Int_t           lumi;
  Int_t           evt;

  Float_t         pt_tag;
  Float_t         eta_tag;
  Float_t         phi_tag;

  Float_t         pt_probe;
  Float_t         eta_probe;
  Float_t         phi_probe;

  Float_t         m_vis;

  Bool_t          id_probe;
  Bool_t          iso_probe;

  Int_t          hlt_1_probe;
  Int_t          hlt_2_probe;
  Int_t          hlt_3_probe;
  Int_t          hlt_4_probe;
  Int_t          hlt_5_probe;
  Int_t          hlt_6_probe;

  Float_t        mcweight;
  Float_t        PUweight;

    // List of branches
  TBranch         *b_run;
  TBranch         *b_lumi;
  TBranch         *b_evt;

  TBranch         *b_pt_tag;
  TBranch         *b_eta_tag;
  TBranch         *b_phi_tag;

  TBranch         *b_pt_probe;
  TBranch         *b_eta_probe;
  TBranch         *b_phi_probe;

  TBranch         *b_m_vis;

  TBranch         *b_id_probe;
  TBranch         *b_iso_probe;
  TBranch         *b_hlt_probe;

  TBranch          *b_hlt_1_probe;
  TBranch          *b_hlt_2_probe;
  TBranch          *b_hlt_3_probe;
  TBranch          *b_hlt_4_probe;
  TBranch          *b_hlt_5_probe;
  TBranch          *b_hlt_6_probe;

  TBranch          *b_mcweight;
  TBranch          *b_PUweight;


  TagProbeTree(TTree *tree=0);
  virtual ~TagProbeTree();

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

#endif //!endif TagProbeTree_h
