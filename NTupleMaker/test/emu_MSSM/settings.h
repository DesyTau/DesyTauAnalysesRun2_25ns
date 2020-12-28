#ifndef Settings_h
#define Settings_h

#include "TString.h"
#include <map>
#include <vector>

using namespace std;

map<TString, double> xsecs = {
  {"WJetsToLNu"  , 52760*1.166}, // NNLO (1)
  {"W1JetsToLNu" , 1.166*8104.}, // NNLO (2)
  {"W2JetsToLNu" , 1.166*2796.}, // NNLO (3)
  {"W3JetsToLNu" , 1.166*993.5}, // NNLO (4)
  {"W4JetsToLNu" , 1.166*544.4}, // NNLO (5)
  {"ZZ" , 12.19},  // LO (17) -> could be improved
  {"WW" , 118.7},  // NNLO QCD (18)
  {"WZ" , 27.68},  // LO (19) -> could be improved
  {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
  {"DY1JetsToLL_M-50"      , 878.7*1.079}, // NNLO (20a)
  {"DY2JetsToLL_M-50"      , 304.4*1.079}, // NNLO (20b)
  {"DY3JetsToLL_M-50"      , 111.5*1.079}, // NNLO (20c)
  {"DY4JetsToLL_M-50"      , 44.03*1.079}, // NNLO (20d)
  {"TT"               , 831.76}, // NNLO (21 inclusive)
  {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
  {"TTToHadronic"     , 377.96}, // NNLO (22)
  {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
  {"ST_t-channel_top_4f"     , 136.02}, // ? (24) -> could be improved
  {"ST_t-channel_antitop_4f" , 80.95}, // ? (25) -> could be improved
  {"ST_tW_top_5f"            , 35.85}, // ? (26) -> could be improved
  {"ST_tW_antitop_5f"        , 35.85}, // ? (27) -> could be improved
  {"VVTo2L2Nu"               , 11.95 },
  {"WWToLNuQQ"               , 49.997 },
  {"WZTo2L2Q"                , 5.595 },
  {"WZTo1L1Nu2Q"             , 10.71 },
  {"WZTo1L3Nu"               , 3.05 },
  {"WZJToLLLNu"               , 4.708 },
  {"WZTo3LNu"                 , 4.43 },
  {"ZZTo4L"                   , 1.212 },
  {"ZZTo2L2Q"                 , 3.22 },  
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
  "EmbeddedElMu_Run2018A",
  "EmbeddedElMu_Run2018B",
  "EmbeddedElMu_Run2018C",
  "EmbeddedElMu_Run2018D"
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
  "EmbeddedElMu_Run2017B",
  "EmbeddedElMu_Run2017C",
  "EmbeddedElMu_Run2017D",
  "EmbeddedElMu_Run2017E",
  "EmbeddedElMu_Run2017F",
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
  "EmbeddedElMu_Run2016B",
  "EmbeddedElMu_Run2016C",
  "EmbeddedElMu_Run2016D",
  "EmbeddedElMu_Run2016E",
  "EmbeddedElMu_Run2016F",
  "EmbeddedElMu_Run2016G",
  "EmbeddedElMu_Run2016H"
};

std::vector<TString> WJets = {
  "WJetsToLNu",
  "W1JetsToLNu",
  "W2JetsToLNu",
  "W3JetsToLNu",
  "W4JetsToLNu"
}; 

std::vector<TString> DYJets = {
  "DYJetsToLL_M-50",
  "DY1JetsToLL_M-50",
  "DY2JetsToLL_M-50",
  "DY3JetsToLL_M-50",
  "DY4JetsToLL_M-50"  
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
  "TT"
};

int nBinsNoBTag = 31;
float binsNoBTag[32] = 
  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};

int nBinsBTag = 17;
float binsBTag[18] =
  {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};

std::vector<TString> masses = {
  "140",
  "180",
  "200",
  "250",
  "350",
  "400",
  "450",
  "600",
  "700",
  "800",
  "900",
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

std::map<TString,double> Masses = {
  {"140",140},
  {"180",180},
  {"200",200},
  {"250",250},
  {"350",350},
  {"400",400},
  {"450",450},
  {"600",600},
  {"700",700},
  {"800",800},
  {"900",900},
  {"1200",1200},
  {"1400",1400},
  {"1600",1600},
  {"1800",1800},
  {"2000",2000},
  {"2300",2300},
  {"2600",2600},
  {"2900",2900},
  {"3200",3200}
};


#endif
