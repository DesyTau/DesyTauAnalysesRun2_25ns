#ifndef diMuon2016Tree_h
#define diMuon2016Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class diMuon2016Tree {
public :
	
	TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain


  // Declaration of leaf types
  Int_t           run;
  Int_t           lumi;
  Int_t           evt;

  Int_t           npv;

  Float_t         pt_0;
  Float_t         eta_0;
  Float_t         phi_0;

  Float_t         pt_1;
  Float_t         eta_1;
  Float_t         phi_1;

  Float_t         m_vis;

    // List of branches
  TBranch         *b_run;
  TBranch         *b_lumi;
  TBranch         *b_evt;

  TBranch         *b_npv;

  TBranch         *b_pt_0;
  TBranch         *b_eta_0;
  TBranch         *b_phi_0;

  TBranch         *b_pt_1;
  TBranch         *b_eta_1;
  TBranch         *b_phi_1;

  TBranch         *b_m_vis;


  diMuon2016Tree(TTree *tree=0);
  virtual ~diMuon2016Tree();

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

#endif //!endif diMuon2016Tree_h
