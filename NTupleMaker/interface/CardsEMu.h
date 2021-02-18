#ifndef DataCardsEMu_h
#define DataCardsEMu_h

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

// Samples = 
// Data - data_obs
// DYJetsToLL - ZL
// WJets - W
// TTbar - TTL
// EWK - VV
// EMB - EMB 
// ggHWW125, 
// qqHWW125, 
// WHWW125,
// ZHWW125,
// ggH125, 
// qqH125 
// WH125,
// ZH125,  
// bbH_$mass 
// ggH_t_$mass
// ggH_b_$mass
// ggH_i_$mass
// ggA_t_$mass
// ggA_b_$mass
// ggA_i_$mass

// sampleToProcess = Data (Data and QCD)
//                 = DYJetsToLL
//                 = EMB 
//                 = MC (EWK WJets)
//                 = SMggH (gg->H H->tautau and H->WW)
//                 = SMothers (qqH, WH, ZH)
//                 = bbH (MSSM)
//                 = ggH_t
//                 = ggH_b
//                 = ggH_i
//                 = ggA_t
//                 = ggA_b
//                 = ggA_i

using namespace std;

template<class T1, class T2>
  vector<T1> extract_first(const map<T1, T2>& v) {
  vector<T1> vFirst;
  for (auto element : v) {
    vFirst.push_back(element.first);
  }
  return vFirst;
};

template<class T1, class T2>
  vector<T2> extract_second(const vector<pair<T1, T2> >& v) {
  vector<T2> vSecond;
  for (auto element : v) {
    vSecond.push_back(element.second);
  }
  return vSecond;
};

class CardsEMu {

 public:

  CardsEMu(TString Sample,
	   TString Era,
	   TString Category,
	   TString InputDir,
	   TString OutputDir,
	   TString Variable,
	   int nbins,
	   double xmin,
	   double xmax,
	   int TriggerOption,
	   bool RunWithSystematics,
	   bool RunOnEmbedded); 

  void SetVariableToPlot(TString var, int nbins, double xmin, double xmax);
  void SetVariableToPlot(TString var, int nbins, double * bins);
  void SetVariableToPlot(TString var, int nbins, vector<double> bins);

  bool Run();
  bool CloseFile();
  void PrintSamples();
  void PrintShapeSystematics();
  void PrintWeightSystematics();
  ~CardsEMu();

 private:

  const map<TString, double> xsecs = {
    {"WJetsToLNu"  , 61526.7}, // NNLO (1)
    {"W1JetsToLNu" , 9370.5}, // NNLO (2)
    {"W2JetsToLNu" , 3170.9}, // NNLO (3)
    {"W3JetsToLNu" , 1132.5}, // NNLO (4)
    {"W4JetsToLNu" , 631.5 }, // NNLO (5)
    {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
    {"DY1JetsToLL_M-50"      , 977.1}, // NNLO (20a)
    {"DY2JetsToLL_M-50"      , 347.3}, // NNLO (20b)
    {"DY3JetsToLL_M-50"      , 126.1}, // NNLO (20c)
    {"DY4JetsToLL_M-50"      , 71.67}, // NNLO (20d)
    {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
    {"TTToHadronic"     , 377.96}, // NNLO (22)
    {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
    {"ST_t-channel_top_4f"      , 136.02}, // ? (24) -> could be improved
    {"ST_t-channel_antitop_4f"  , 80.95}, // ? (25) -> could be improved
    {"ST_tW_top_5f"             , 35.85}, // ? (26) -> could be improved
    {"ST_tW_antitop_5f"         , 35.85}, // ? (27) -> could be improved
    {"VVTo2L2Nu"                , 13.84},
    {"WWToLNuQQ"                , 49.997},
    {"WZTo2L2Q"                 , 5.52},
    {"WZTo1L1Nu2Q"              , 10.71},
    {"WZTo1L3Nu"                , 3.05},
    {"WZJToLLLNu"               , 4.708},
    {"WZTo3LNu"                 , 4.43},
    {"ZZTo4L"                   , 1.26},
    {"ZZTo2L2Q"                 , 3.38},
    {"GluGluHToTauTau_M125"     , 3.00},
    {"VBFHToTauTau_M125"        , 0.237},
    {"WplusHToTauTau_M125"      , 0.0527},
    {"WminusHToTauTau_M125"     , 0.0334},
    {"ZHToTauTau_M125_13TeV"    , 0.0477},
    {"GluGluHToWWTo2L2Nu_M125"  , 1.09},
    {"VBFHToWWTo2L2Nu_M125"     , 0.0850},
    {"HWminusJ_HToWW_M125"      , 0.114},
    {"HWplusJ_HToWW_M125"       , 0.18},
    {"ZHJ_HToWW_M125"           , 0.163}
  };

  const vector<TString> SingleMuon_2018 = {
    "SingleMuon_Run2018A",
    "SingleMuon_Run2018B",
    "SingleMuon_Run2018C",
    "SingleMuon_Run2018D"
  };
  
  const vector<TString> SingleElectron_2018 = {
    "EGamma_Run2018A",
    "EGamma_Run2018B",
    "EGamma_Run2018C",
    "EGamma_Run2018D"
  };
  
  const vector<TString> MuonEG_2018 = {
    "MuonEG_Run2018A",
    "MuonEG_Run2018B",
    "MuonEG_Run2018C",
    "MuonEG_Run2018D"
  };
  
  const vector<TString> EmbeddedElMu_2018 = {
    "EmbeddedElMu_Run2018A",
    "EmbeddedElMu_Run2018B",
    "EmbeddedElMu_Run2018C",
    "EmbeddedElMu_Run2018D"
  };
  
  // ****** 2017 ********

  const vector<TString> SingleMuon_2017 = {
    "SingleMuon_Run2017B",
    "SingleMuon_Run2017C",
    "SingleMuon_Run2017D",
    "SingleMuon_Run2017E",
    "SingleMuon_Run2017F"
  };
  
  const vector<TString> SingleElectron_2017 = {
    "SingleElectron_Run2017B",
    "SingleElectron_Run2017C",
    "SingleElectron_Run2017D",
    "SingleElectron_Run2017E",
    "SingleElectron_Run2017F",
  };

  const vector<TString> MuonEG_2017 = {
    "MuonEG_Run2017B",
    "MuonEG_Run2017C",
    "MuonEG_Run2017D",
    "MuonEG_Run2017E",
    "MuonEG_Run2017F"
  };

  const vector<TString> EmbeddedElMu_2017 = {
    "EmbeddedElMu_Run2017B",
    "EmbeddedElMu_Run2017C",
    "EmbeddedElMu_Run2017D",
    "EmbeddedElMu_Run2017E",
    "EmbeddedElMu_Run2017F",
  };
  
  // ******* 2016 ******

  const vector<TString> SingleMuon_2016 = {
    "SingleMuon_Run2016B",
    "SingleMuon_Run2016C",
    "SingleMuon_Run2016D",
    "SingleMuon_Run2016E",
    "SingleMuon_Run2016F",
    "SingleMuon_Run2016G",
    "SingleMuon_Run2016H"
  };
  
  const vector<TString> SingleElectron_2016 = {
    "SingleElectron_Run2016B",
    "SingleElectron_Run2016C",
    "SingleElectron_Run2016D",
    "SingleElectron_Run2016E",
    "SingleElectron_Run2016F",
    "SingleElectron_Run2016G",
    "SingleElectron_Run2016H"
  };

  const vector<TString> MuonEG_2016 = {
    "MuonEG_Run2016B",
    "MuonEG_Run2016C",
    "MuonEG_Run2016D",
    "MuonEG_Run2016E",
    "MuonEG_Run2016F",
    "MuonEG_Run2016G",
    "MuonEG_Run2016H"
  };

  const vector<TString> EmbeddedElMu_2016 = {
    "EmbeddedElMu_Run2016B",
    "EmbeddedElMu_Run2016C",
    "EmbeddedElMu_Run2016D",
    "EmbeddedElMu_Run2016E",
    "EmbeddedElMu_Run2016F",
    "EmbeddedElMu_Run2016G",
    "EmbeddedElMu_Run2016H"
  };
  
  const vector<TString> WJets = {
    "WJetsToLNu_0",
    "WJetsToLNu_1",
    "WJetsToLNu_2",
    "WJetsToLNu_3",
    "WJetsToLNu_4",
    "W1JetsToLNu",
    "W2JetsToLNu",
    "W3JetsToLNu",
    "W4JetsToLNu"
  };
  
  vector<TString> DYJetsToLL = {
    "DYJetsToLL_M-50_0",
    "DYJetsToLL_M-50_1",
    "DYJetsToLL_M-50_2",
    "DYJetsToLL_M-50_3",
    "DYJetsToLL_M-50_4",
    "DY1JetsToLL_M-50",
    "DY2JetsToLL_M-50",
    "DY3JetsToLL_M-50",
    "DY4JetsToLL_M-50"  
  };

  vector<TString> DYJetsToTT = {
    "DYJetsToTT_M-50_0",
    "DYJetsToTT_M-50_1",
    "DYJetsToTT_M-50_2",
    "DYJetsToTT_M-50_3",
    "DYJetsToTT_M-50_4",
    "DY1JetsToTT_M-50",
    "DY2JetsToTT_M-50",
    "DY3JetsToTT_M-50",
    "DY4JetsToTT_M-50"  
  };
  
  /*
    vector<TString> DYJets = {
    "DYJetsToLL_M-50",
    "DY1JetsToLL_M-50",
    "DY2JetsToLL_M-50",
    "DY3JetsToLL_M-50",
    "DY4JetsToLL_M-50"  
    };
  */
  
  // needed for stitching 
  map<TString,TString> WJetsFiles = {
    {"WJetsToLNu_0","WJetsToLNu"},
    {"WJetsToLNu_1","WJetsToLNu"},
    {"WJetsToLNu_2","WJetsToLNu"},
    {"WJetsToLNu_3","WJetsToLNu"},
    {"WJetsToLNu_4","WJetsToLNu"},
    {"W1JetsToLNu","W1JetsToLNu"},
    {"W2JetsToLNu","W2JetsToLNu"},
    {"W3JetsToLNu","W3JetsToLNu"},
    {"W4JetsToLNu","W4JetsToLNu"}
  };

  map<TString,TString> DYJetsLLFiles = {
    {"DYJetsToLL_M-50_0","DYJetsToLL_M-50"},
    {"DYJetsToLL_M-50_1","DYJetsToLL_M-50"},
    {"DYJetsToLL_M-50_2","DYJetsToLL_M-50"},
    {"DYJetsToLL_M-50_3","DYJetsToLL_M-50"},
    {"DYJetsToLL_M-50_4","DYJetsToLL_M-50"},
    {"DY1JetsToLL_M-50","DY1JetsToLL_M-50"},
    {"DY2JetsToLL_M-50","DY2JetsToLL_M-50"},
    {"DY3JetsToLL_M-50","DY3JetsToLL_M-50"},
    {"DY4JetsToLL_M-50","DY4JetsToLL_M-50"},
  };

  map<TString,TString> DYJetsTTFiles = {
    {"DYJetsToTT_M-50_0","DYJetsToLL_M-50"},
    {"DYJetsToTT_M-50_1","DYJetsToLL_M-50"},
    {"DYJetsToTT_M-50_2","DYJetsToLL_M-50"},
    {"DYJetsToTT_M-50_3","DYJetsToLL_M-50"},
    {"DYJetsToTT_M-50_4","DYJetsToLL_M-50"},
    {"DY1JetsToTT_M-50","DY1JetsToLL_M-50"},
    {"DY2JetsToTT_M-50","DY2JetsToLL_M-50"},
    {"DY3JetsToTT_M-50","DY3JetsToLL_M-50"},
    {"DY4JetsToTT_M-50","DY4JetsToLL_M-50"},
  };
  
  vector<TString> EWK = {
    "ST_t-channel_top_4f",
    "ST_t-channel_antitop_4f",
    "ST_tW_top_5f",
    "ST_tW_antitop_5f",
    "VVTo2L2Nu",
    "WZTo2L2Q",
    "WZTo3LNu",
    "ZZTo2L2Q",
    "ZZTo4L",
  };

  vector<TString> TT = {
    "TTTo2L2Nu",
    "TTToHadronic",
    "TTToSemiLeptonic"
  };
  
  vector<TString> GluGluHToTauTau = {
    "GluGluHToTauTau_M125"
  };
  
  vector<TString> VBFHToTauTau = {
    "VBFHToTauTau_M125"
  };

  vector<TString> WHToTauTau = {
    "WplusHToTauTau_M125",
    "WminusHToTauTau_M125"
  };

  vector<TString> ZHToTauTau = {
    "ZHToTauTau_M125_13TeV"
  };

  vector<TString> GluGluHToWW = {
    "GluGluHToWWTo2L2Nu_M125"
  };
  
  vector<TString> VBFHToWW = {
    "VBFHToWWTo2L2Nu_M125"
  };

  vector<TString> WHToWW = {
    "HWminusJ_HToWW_M125",
    "HWplusJ_HToWW_M125"
  };

  vector<TString> ZHToWW = {
    "ZHJ_HToWW_M125"
  };

  const vector<TString> masses = {
    "110",
    "120",
    "125",
    "130",
    "140",
    "160",
    "180",
    "200",
    "250",
    "300",
    "350",
    "400",
    "450",
    "500",
    "600",
    "700",
    "800",
    "900",
    "1000",
    "1200",
    "1400",
    "1600",
    "1800",
    "2000",
    "2300",
    "2600",
    "2900",
    "3200"
  };

  const TString BaseTreeName = "TauCheck"; 
  const TString baseNameBBH  = "SUSYGluGluToBBHToTauTau_M-";
  const TString baseNameGGH  = "SUSYGluGluToHToTauTau_M-"; 
  const TString notTauTau = "&&!(gen_match_1==3&&gen_match_2==4)";
  const TString TauTau = "&&(gen_match_1==3&&gen_match_2==4)";

  map<TString,TString> generatorName = 
    {
      {"SUSYGluGluToBBHToTauTau_M-","amcatnlo"},
      {"SUSYGluGluToHToTauTau_M-","pythia"}
    };

  // *****************************
  // ***** Shape systematics *****
  // *****************************

  const map<TString,TString> EmbeddedShapeSystematics = {
    {"CMS_scale_e_emb","CMS_scale_e_13TeV"},
    {"CMS_scale_m","CMS_scale_mu_13TeV"},
  };

  const map<TString,TString> EmbeddedMetShapeSystematics = {
    {"scale_embed_metDown","CMS_scale_met_embedded_13TeVUp"},
    {"scale_embed_metUp","CMS_scale_met_embedded_13TeVDown"}
  };

  const map<TString,TString> ShapeSystematics_Common = {
    {"CMS_scale_e","CMS_scale_e_13TeV"},
    {"CMS_scale_m","CMS_scale_mu_13TeV"},
    {"CMS_scale_j_FlavorQCD","CMS_scale_j_FlavorQCD_13TeV"},
    {"CMS_scale_j_RelativeBal","CMS_scale_j_RelativeBal_13TeV"},
    {"CMS_scale_j_HF","CMS_scale_j_HF_13TeV"},
    {"CMS_scale_j_BBEC1","CMS_scale_j_BBEC1_13TeV"},
    {"CMS_scale_j_EC2","CMS_scale_j_EC2_13TeV"},
    {"CMS_scale_j_Absolute","CMS_scale_j_Absolute_13TeV"},    
  };

  const map<TString,TString> ShapeSystematics_2016 = {
    {"CMS_scale_met_unclustered_2016","CMS_scale_met_unclustered_13TeV"},
    {"CMS_htt_boson_res_met_2016","CMS_htt_boson_reso_met_13TeV"},
    {"CMS_htt_boson_scale_met_2016","CMS_htt_boson_scale_met_13TeV"},
    {"CMS_htt_eff_b_2016","CMS_eff_b_13TeV"},
    {"CMS_htt_mistag_b_2016","CMS_mistag_b_13TeV"},
    {"CMS_scale_j_Absolute_2016","CMS_scale_j_Absolute_2016_13TeV"},
    {"CMS_scale_j_HF_2016","CMS_scale_j_HF_2016_13TeV"},
    {"CMS_scale_j_EC2_2016","CMS_scale_j_EC2_2016_13TeV"},
    {"CMS_scale_j_RelativeSample_2016","CMS_scale_j_RelativeSample_2016_13TeV"},
    {"CMS_scale_j_BBEC1_2016","CMS_scale_j_BBEC1_2016_13TeV"},
    {"CMS_res_j_2016","CMS_res_j_13TeV"},
  };  

  const map<TString,TString> ShapeSystematics_2017 = {
    {"CMS_scale_met_unclustered_2017","CMS_scale_met_unclustered_13TeV"},
    {"CMS_htt_boson_res_met_2017","CMS_htt_boson_reso_met_13TeV"},
    {"CMS_htt_boson_scale_met_2017","CMS_htt_boson_scale_met_13TeV"},
    {"CMS_htt_eff_b_2017","CMS_eff_b_13TeV"},
    {"CMS_htt_mistag_b_2017","CMS_mistag_b_13TeV"},
    {"CMS_scale_j_Absolute_2017","CMS_scale_j_Absolute_2017_13TeV"},
    {"CMS_scale_j_HF_2017","CMS_scale_j_HF_2017_13TeV"},
    {"CMS_scale_j_EC2_2017","CMS_scale_j_EC2_2017_13TeV"},
    {"CMS_scale_j_RelativeSample_2017","CMS_scale_j_RelativeSample_2017_13TeV"},
    {"CMS_scale_j_BBEC1_2017","CMS_scale_j_BBEC1_2017_13TeV"},
    {"CMS_res_j_2017","CMS_res_j_13TeV"},
  };  

  const map<TString,TString> ShapeSystematics_2018 = {
    {"CMS_scale_met_unclustered_2018","CMS_scale_met_unclustered_13TeV"},
    {"CMS_htt_boson_res_met_2018","CMS_htt_boson_reso_met_13TeV"},
    {"CMS_htt_boson_scale_met_2018","CMS_htt_boson_scale_met_13TeV"},
    {"CMS_htt_eff_b_2018","CMS_eff_b_13TeV"},
    {"CMS_htt_mistag_b_2018","CMS_mistag_b_13TeV"},
    {"CMS_scale_j_Absolute_2018","CMS_scale_j_Absolute_2018_13TeV"},
    {"CMS_scale_j_HF_2018","CMS_scale_j_HF_2018_13TeV"},
    {"CMS_scale_j_EC2_2018","CMS_scale_j_EC2_2018_13TeV"},
    {"CMS_scale_j_RelativeSample_2018","CMS_scale_j_RelativeSample_2018_13TeV"},
    {"CMS_scale_j_BBEC1_2018","CMS_scale_j_BBEC1_2018_13TeV"},
    {"CMS_res_j_2018","CMS_res_j_13TeV"},
  };  

  // ******************************
  // *** Weight systematics *******
  // ******************************

  const map<TString,TString> PrefiringSystematics = {
    {"CMS_prefiringUp","(prefiringweightUp/prefiringweight)"},
    {"CMS_prefiringDown","(prefiringweightDown/prefiringweight)"},    
  };

  const map<TString,TString> DYShapeSystematics_2016 = {
    {"CMS_htt_dyShape_2016Up","zptweight"},
    {"CMS_htt_dyShape_2016Down","(1/zptweight)"},
  };

  const map<TString,TString> DYShapeSystematics = {
    {"CMS_htt_dyShapeUp","zptweight"},
    {"CMS_htt_dyShapeDown","(1/zptweight)"},
  };

  const map<TString,TString> TopShapeSystematics = {
    {"CMS_htt_ttbarShapeUp","topptweight"},
    {"CMS_htt_ttbarShapeDown","(1/topptweight)"},
  };

  const map<TString,TString> SignalSystematics = {
    {"CMS_scale_ggH","weight_CMS_scale_gg_13TeV"},
    {"CMS_PS_ISR_ggH","weight_CMS_PS_ISR_ggH_13TeV"},
    {"CMS_PS_FSR_ggH","weight_CMS_PS_FSR_ggH_13TeV"},
  };

  const map<TString,TString> QCDWeightSystematics = {
    {"CMS_htt_qcd_0jet_rate","qcdweight_deltaR_0jet_Par0"},
    {"CMS_htt_qcd_0jet_shape","qcdweight_deltaR_0jet_Par1"},
    {"CMS_htt_qcd_0jet_shape2","qcdweight_deltaR_0jet_Par2"},
    {"CMS_htt_qcd_1jet_rate","qcdweight_deltaR_1jet_Par0"},
    {"CMS_htt_qcd_1jet_shape","qcdweight_deltaR_1jet_Par1"},
    {"CMS_htt_qcd_1jet_shape2","qcdweight_deltaR_1jet_Par2"},
    {"CMS_htt_qcd_2jet_rate","qcdweight_deltaR_2jet_Par0"},
    {"CMS_htt_qcd_2jet_shape","qcdweight_deltaR_2jet_Par1"},
    {"CMS_htt_qcd_2jet_shape2","qcdweight_deltaR_0jet_Par2"},
  };

  const map<TString,TString> QCDIsoSystematics = {
    {"CMS_htt_qcd_isoUp","qcdweight_isolationcorrection"},
    {"CMS_htt_qcd_isoDown","(1.0/qcdweight_isolationcorrection)"},    
  };

  map<TString,TString> weight_ggH = {
    {"ggH_t","H_t_ratio"},
    {"ggH_b","H_b_ratio"},
    {"ggH_i","H_i_ratio"},
    {"ggh_t","H_t_ratio"},
    {"ggh_b","H_b_ratio"},
    {"ggh_i","H_i_ratio"},
    {"ggA_t","A_t_ratio"},
    {"ggA_b","A_b_ratio"},
    {"ggA_i","A_i_ratio"}
  };

  map<TString, double> scaleQCD = {
    {"2016",0.71},
    {"2017",0.69},
    {"2018",0.67}
  };

  const map<TString,TString> systematics_ggH = {
    {"QCDscale_ggH_REWEIGHT","scale"},
    {"Hdamp_ggH_REWEIGHT","hdamp"}
  };

  // ******************************************************* 
  // ************* CLASS VARIABLES *************************
  // *******************************************************

  TString TriggerCut; // trigger cut for MC
  TString TriggerCutEMu; // trigger cut for MuonEG
  TString TriggerCutSingleE; // trigger cut for SingleE
  TString TriggerCutSingleMu; // trigger cut for SingleMuon

  map<TString, TString> weightSystematicsMap;
  map<TString, TString> shapeSystematicsMap;
  map<TString, double> sampleNormMap; // sample-norm map
  map<TString, TFile*> sampleFileMap; // sample-file map
  map<TString, vector<TString> > nameSampleMap; // sample map
  map<TString, TString> sampleSpecificCutMap; // sample-specific-cut map
  map<TString, double> sampleXSecMap;
  map<TString, double> sampleNeventsMap;
  map<TString, TString> nameHistoMap;
  map<TString, TH1D*> nameTH1DMap;
  vector<TString> samplesContainer;

  int numberOfShapeSystematics;
  int numberOfWeightSystematics;

  // for variables which are different from mt_tot; 
  int nBins;
  double Bins[100];
  double xMin;
  double xMax;

  TString sampleToProcess;
  TString input_dir; // includes era
  TString output_filename; // outputdir + era / sampleToProcess + Category;
  TFile * outputFile; 

  TString globalWeight;
  TString commonCuts;
  TString category;
  TString era;

  bool block;
  bool isBTag;
  bool isEquidistantBinning;
  bool runWithSystematics;
  bool runOnEmbedded;

  TString variable;
  double lumi;
  int triggerOption; // 0 - emu, 1 - singleLep, 2 - comb

  // **************************************************
  // ************* Internal member functions **********
  // **************************************************

  TString SampleSpecificCut(TString name, 
			    TString sampleName);

  void InitializeSample(TString name);

  void CreateMSSMList(TString baseName, TString templName);
  void InitializeSamples();

  TH1D * ProcessSample(TString name,
		       TString sysName,
		       bool OS,
		       bool weightSys);

  TH1D * createHisto(TString histName);

  int CreateShapeSystematicsMap();
  int CreateWeightSystematicsMap();
  bool RunData();
  bool RunModel();
  void zeroBins(TH1D * hist);

};

#endif
