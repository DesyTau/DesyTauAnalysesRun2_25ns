#ifndef ZMuMuTree_h
#define ZMuMuTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ZMuMuTree {
public :
	
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //!current Tree number in a TChain


  // Declaration of leaf types
  ULong64_t       run; //
  Int_t           lumi; //
  Int_t           evt; //
  Int_t 		  npv; //

   // leading and trailing muons indicated with _1 and _2, respectively
  Float_t         pt_1; //
  Float_t         eta_1; //
  Float_t         phi_1; //
  Float_t         iso_1; //
  Float_t         q_1; //

  Float_t         pt_2; //
  Float_t         eta_2; //
  Float_t         phi_2; // 
  Float_t         iso_2; //
  Float_t         q_2; //

  Float_t         os; //
  Float_t         pt_ll; //
  Float_t         px_ll; //
  Float_t         py_ll; //
  Float_t         m_ll; // 
  Float_t         eta_ll; //
  Float_t         phi_ll; //
  Float_t         pt_ratio_ll; // 
  Float_t         delta_phi_ll; //

  Float_t         genV_px; //
  Float_t         genV_py; //
  Float_t         genV_pt; //
  Float_t         genV_mass; //
  Float_t         genV_eta; //
  Float_t         genV_phi; //
  Float_t         genZ_pt; //
  Float_t         genZ_mass; //
  Float_t         is_gen_Zee; //
  Float_t         is_gen_Zmm; //
  Float_t         is_gen_Ztt; //

   // leading and trailing jets indicated with _1 and _2, respectively
  Float_t         pt_jet_1; //
  Float_t         eta_jet_1; //
  Float_t         phi_jet_1; //

  Float_t         pt_jet_2; //
  Float_t         eta_jet_2; //
  Float_t         phi_jet_2; //

  Float_t         m_jj; //
  Float_t         delta_phi_jj; //

  Int_t           njets30; //
  Int_t           njets20; //
  Float_t         HT30; //
  Float_t         HT20; // 

  // _rcorr for values after recoil corrections are applied. 
  Float_t         met; //
  Float_t         met_phi; //
  Float_t         puppimet; //
  Float_t         puppimet_phi; //
  Float_t         met_rcorr; //
  Float_t         met_phi_rcorr; //
  Float_t         puppimet_rcorr; //
  Float_t         puppimet_phi_rcorr; // 
  Float_t         met_response_gen; //
  Float_t         puppimet_response_gen; //

  Float_t         met_recoil_paral; //
  Float_t         met_recoil_perp; //
  Float_t         met_response_had; //
  Float_t         puppimet_recoil_paral; //
  Float_t         puppimet_recoil_perp; //
  Float_t         puppimet_response_had;  //

  Float_t         met_response_gen_rcorr;
  Float_t         puppimet_response_gen_rcorr;
  Float_t         met_recoil_paral_rcorr;
  Float_t         met_recoil_perp_rcorr;
  Float_t         met_response_had_rcorr; 
  Float_t         puppimet_recoil_paral_rcorr;
  Float_t         puppimet_recoil_perp_rcorr;
  Float_t         puppimet_response_had_rcorr; 

  Float_t         mc_weight; //
  Float_t         pu_weight; // 
  Double_t        toppt_weight; //
  Double_t        zpt_weight; //
  Double_t        track_weight_1; //
  Double_t        track_weight_2; //
  Double_t        idiso_weight_1; //
  Double_t        idiso_weight_2; //
  Double_t        trig_weight; //

    // List of branches
  TBranch         *b_run;
  TBranch         *b_lumi;
  TBranch         *b_evt;
  TBranch 		  *b_npv;

  TBranch         *b_pt_1;
  TBranch         *b_eta_1;
  TBranch         *b_phi_1;
  TBranch         *b_iso_1;
  TBranch         *b_q_1;

  TBranch         *b_pt_2;
  TBranch         *b_eta_2;
  TBranch         *b_phi_2;
  TBranch         *b_iso_2;
  TBranch         *b_q_2;

  TBranch         *b_os; 
  TBranch         *b_pt_ll;
  TBranch         *b_px_ll;
  TBranch         *b_py_ll;
  TBranch         *b_m_ll;
  TBranch         *b_eta_ll;
  TBranch         *b_phi_ll;
  TBranch         *b_delta_phi_ll;

  TBranch         *b_genV_px;
  TBranch         *b_genV_py;
  TBranch         *b_genV_pt;
  TBranch         *b_genV_mass;
  TBranch         *b_genV_eta;
  TBranch         *b_genV_phi;
  TBranch         *b_genZ_pt;
  TBranch         *b_genZ_mass;
  TBranch         *b_is_gen_Zee;
  TBranch         *b_is_gen_Zmm;
  TBranch         *b_is_gen_Ztt;

  TBranch         *b_pt_jet_1;
  TBranch         *b_eta_jet_1;
  TBranch         *b_phi_jet_1;
  TBranch         *b_pt_jet_2;
  TBranch         *b_eta_jet_2;
  TBranch         *b_phi_jet_2;

  TBranch         *b_m_jj;
  TBranch         *b_delta_phi_jj;
  TBranch         *b_njets30;
  TBranch         *b_njets20;
  TBranch         *b_HT30;
  TBranch         *b_HT20;

  TBranch         *b_met;
  TBranch         *b_met_phi;
  TBranch         *b_puppimet;
  TBranch         *b_puppimet_phi;
  TBranch         *b_met_rcorr;
  TBranch         *b_met_phi_rcorr;
  TBranch         *b_puppimet_rcorr;
  TBranch         *b_puppimet_phi_rcorr;
  TBranch         *b_pt_ratio_ll;
  TBranch         *b_met_response_gen;
  TBranch         *b_puppimet_response_gen;
  TBranch         *b_met_recoil_paral;
  TBranch         *b_met_recoil_perp;
  TBranch         *b_met_response_had;
  TBranch         *b_puppimet_recoil_paral;
  TBranch         *b_puppimet_recoil_perp;
  TBranch         *b_puppimet_response_had;


  TBranch         *b_met_response_gen_rcorr;
  TBranch         *b_puppimet_response_gen_rcorr;
  TBranch         *b_met_recoil_paral_rcorr;
  TBranch         *b_met_recoil_perp_rcorr;
  TBranch         *b_met_response_had_rcorr; 
  TBranch         *b_puppimet_recoil_paral_rcorr;
  TBranch         *b_puppimet_recoil_perp_rcorr;
  TBranch         *b_puppimet_response_had_rcorr; 

  TBranch         *b_mc_weight;
  TBranch         *b_pu_weight;
  TBranch         *b_toppt_weight;
  TBranch         *b_zpt_weight;
  TBranch         *b_track_weight_1;
  TBranch         *b_track_weight_2;
  TBranch         *b_idiso_weight_1;
  TBranch         *b_idiso_weight_2;
  TBranch         *b_trig_weight;

  ZMuMuTree(TTree *tree=0);
  virtual ~ZMuMuTree();

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

#endif //!endif ZMuMuTree_h
