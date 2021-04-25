#ifndef Settings_h
#define Settings_h

#include "TString.h"
#include <map>
#include <vector>

using namespace std;

map<TString, double> xsecs = {
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


std::vector<TString> SingleMuon_2018 = {
  "SingleMuon_Run2018A",
  "SingleMuon_Run2018B",
  "SingleMuon_Run2018C",
  "SingleMuon_Run2018D"
};

std::vector<TString> SingleElectron_2018 = {
  "EGamma_Run2018A",
  "EGamma_Run2018B",
  "EGamma_Run2018C",
  "EGamma_Run2018D"
};

std::vector<TString> MuonEG_2018 = {
  "MuonEG_Run2018A",
  "MuonEG_Run2018B",
  "MuonEG_Run2018C",
  "MuonEG_Run2018D"
};

std::vector<TString> EmbeddedElMu_2018 = {
  "EmbeddedElMu_Run2018A_new",
  "EmbeddedElMu_Run2018B_new",
  "EmbeddedElMu_Run2018C_new",
  "EmbeddedElMu_Run2018D_new"
};

// ****** 2017 ********

std::vector<TString> SingleMuon_2017 = {
  "SingleMuon_Run2017B",
  "SingleMuon_Run2017C",
  "SingleMuon_Run2017D",
  "SingleMuon_Run2017E",
  "SingleMuon_Run2017F"
};

std::vector<TString> SingleElectron_2017 = {
  "SingleElectron_Run2017B",
  "SingleElectron_Run2017C",
  "SingleElectron_Run2017D",
  "SingleElectron_Run2017E",
  "SingleElectron_Run2017F",
};

std::vector<TString> MuonEG_2017 = {
  "MuonEG_Run2017B",
  "MuonEG_Run2017C",
  "MuonEG_Run2017D",
  "MuonEG_Run2017E",
  "MuonEG_Run2017F"
};

std::vector<TString> EmbeddedElMu_2017 = {
  "EmbeddedElMu_Run2017B_new",
  "EmbeddedElMu_Run2017C_new",
  "EmbeddedElMu_Run2017D_new",
  "EmbeddedElMu_Run2017E_new",
  "EmbeddedElMu_Run2017F_new",
};

// ******* 2016 ******

std::vector<TString> SingleMuon_2016 = {
  "SingleMuon_Run2016B",
  "SingleMuon_Run2016C",
  "SingleMuon_Run2016D",
  "SingleMuon_Run2016E",
  "SingleMuon_Run2016F",
  "SingleMuon_Run2016G",
  "SingleMuon_Run2016H"
};

std::vector<TString> SingleElectron_2016 = {
  "SingleElectron_Run2016B",
  "SingleElectron_Run2016C",
  "SingleElectron_Run2016D",
  "SingleElectron_Run2016E",
  "SingleElectron_Run2016F",
  "SingleElectron_Run2016G",
  "SingleElectron_Run2016H"
};

std::vector<TString> MuonEG_2016 = {
  "MuonEG_Run2016B",
  "MuonEG_Run2016C",
  "MuonEG_Run2016D",
  "MuonEG_Run2016E",
  "MuonEG_Run2016F",
  "MuonEG_Run2016G",
  "MuonEG_Run2016H"
};

std::vector<TString> EmbeddedElMu_2016 = {
  "EmbeddedElMu_Run2016B_new",
  "EmbeddedElMu_Run2016C_new",
  "EmbeddedElMu_Run2016D_new",
  "EmbeddedElMu_Run2016E_new",
  "EmbeddedElMu_Run2016F_new",
  "EmbeddedElMu_Run2016G_new",
  "EmbeddedElMu_Run2016H_new"
};

std::vector<TString> WJets = {
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
 
std::vector<TString> DYJets = {
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

std::map<TString,TString> DYJetsFiles = {
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

std::map<TString,TString> WJetsFiles = {
  {"WJetsToLNu_0","WJetsToLNu"},
  {"WJetsToLNu_1","WJetsToLNu"},
  {"WJetsToLNu_2","WJetsToLNu"},
  {"WJetsToLNu_3","WJetsToLNu"},
  {"WJetsToLNu_4","WJetsToLNu"},
  {"W1JetsToLNu","W1JetsToLNu"},
  {"W2JetsToLNu","W2JetsToLNu"},
  {"W3JetsToLNu","W3JetsToLNu"},
  {"W4JetsToLNu","W4JetsToLNu"},
};

std::vector<TString> EWK = {
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

std::vector<TString> TT_EXCL = {
  "TTTo2L2Nu",
  "TTToHadronic",
  "TTToSemiLeptonic"
};

std::vector<TString> TT_INCL = {
  "TT_INCL"
};

TString baseNameBBH = "SUSYGluGluToBBHToTauTau_M-";
TString baseNameGGH = "SUSYGluGluToHToTauTau_M-"; 

std::vector<TString> masses = {
  "60",
  "80",
  "100",
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
  "3200",
  "3500"
};

const TString BaseTreeName = "TauCheck"; 
vector<TString> SystematicsNames = {"",  
				    "CMS_scale_met_unclustered_13TeVUp",
				    "CMS_scale_met_unclustered_13TeVDown",
				    "CMS_scale_met_boson_resolution_13TeVUp",
				    "CMS_scale_met_boson_resolution_13TeVDown",
				    "CMS_scale_met_boson_response_13TeVUp",
				    "CMS_scale_met_boson_response_13TeVDown",
				    "CMS_htt_boson_reso_met_13TeVUp",
				    "CMS_htt_boson_reso_met_13TeVDown",
				    "CMS_htt_boson_scale_met_13TeVUp",
				    "CMS_htt_boson_scale_met_13TeVDown",
				    "CMS_eff_b_13TeVUp",
				    "CMS_eff_b_13TeVDown",
				    "CMS_eff_mistag_13TeVUp",
				    "CMS_eff_mistag_13TeVDown",
				    "CMS_scale_j_FlavorQCD_13TeVUp",
				    "CMS_scale_j_FlavorQCD_13TeVDown",
				    "CMS_scale_j_RelativeBal_13TeVUp",
				    "CMS_scale_j_RelativeBal_13TeVDown",
				    "CMS_scale_j_HF_13TeVUp",
				    "CMS_scale_j_HF_13TeVDown",
				    "CMS_scale_j_BBEC1_13TeVUp",
				    "CMS_scale_j_BBEC1_13TeVDown",
				    "CMS_scale_j_EC2_13TeVUp",
				    "CMS_scale_j_EC2_13TeVDown",
				    "CMS_scale_j_Absolute_13TeVUp",
				    "CMS_scale_j_Absolute_13TeVDown",
				    "CMS_scale_j_Absolute_2016_13TeVUp",
				    "CMS_scale_j_Absolute_2016_13TeVDown",
				    "CMS_scale_j_HF_2016_13TeVUp",
				    "CMS_scale_j_HF_2016_13TeVDown",
				    "CMS_scale_j_EC2_2016_13TeVUp",
				    "CMS_scale_j_EC2_2016_13TeVDown",
				    "CMS_scale_j_RelativeSample_2016_13TeVUp",
				    "CMS_scale_j_RelativeSample_2016_13TeVDown",
				    "CMS_scale_j_BBEC1_2016_13TeVUp",
				    "CMS_scale_j_BBEC1_2016_13TeVDown",
				    "CMS_scale_j_Absolute_2017_13TeVUp",
				    "CMS_scale_j_Absolute_2017_13TeVDown",
				    "CMS_scale_j_HF_2017_13TeVUp",
				    "CMS_scale_j_HF_2017_13TeVDown",
				    "CMS_scale_j_EC2_2017_13TeVUp",
				    "CMS_scale_j_EC2_2017_13TeVDown",
				    "CMS_scale_j_RelativeSample_2017_13TeVUp",
				    "CMS_scale_j_RelativeSample_2017_13TeVDown",
				    "CMS_scale_j_BBEC1_2017_13TeVUp",
				    "CMS_scale_j_BBEC1_2017_13TeVDown",
				    "CMS_scale_j_Absolute_2018_13TeVUp",
				    "CMS_scale_j_Absolute_2018_13TeVDown",
				    "CMS_scale_j_HF_2018_13TeVUp",
				    "CMS_scale_j_HF_2018_13TeVDown",
				    "CMS_scale_j_EC2_2018_13TeVUp",
				    "CMS_scale_j_EC2_2018_13TeVDown",
				    "CMS_scale_j_RelativeSample_2018_13TeVUp",
				    "CMS_scale_j_RelativeSample_2018_13TeVDown",
				    "CMS_scale_j_BBEC1_2018_13TeVUp",
				    "CMS_scale_j_BBEC1_2018_13TeVDown",
				    "CMS_scale_mu_13TeVUp",
				    "CMS_scale_mu_13TeVDown",
				    "CMS_scale_e_13TeVUp",
				    "CMS_scale_e_13TeVDown",
				    "CMS_res_j_13TeVUp",
				    "CMS_res_j_13TeVDown"
};

int nBinsNoBTag = 31;
float binsNoBTag[32] = 
  {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};

int nBinsBTag = 17;
float binsBTag[18] =
  {20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};

#endif
