#include <iostream>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

double getNEventsProcessed(TFile * file)
{
  TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
  double nevents = histWeightsH->GetSumOfWeights();
  return nevents;
}

void createDNNinput(TString inputDir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu2016/NTuples"){

  vector<TString> SingleTop;
  SingleTop.push_back("ST_t-channel_antitop");
  SingleTop.push_back("ST_t-channel_top");
  SingleTop.push_back("ST_tW_antitop");
  SingleTop.push_back("ST_tW_antitop");
  vector<TString> TTBar;
  TTBar.push_back("TTBar");
  vector<TString> Diboson;
  Diboson.push_back("VVTo2L2Nu");
  Diboson.push_back("WZJToLLLNu");
  Diboson.push_back("WZTo1L1Nu2Q");
  Diboson.push_back("WZTo1L3Nu");
  Diboson.push_back("WZTo2L2Q");
  Diboson.push_back("ZZTo2L2Q");
  Diboson.push_back("ZZTo4L");
  Diboson.push_back("WWToLNuQQ");
  vector<TString> WJets;
  WJets.push_back("W1JetsToLNu");
  WJets.push_back("W2JetsToLNu");
  WJets.push_back("W3JetsToLNu");
  WJets.push_back("W4JetsToLNu");
  WJets.push_back("WJetsToLNu");
  vector<TString> DYJets;
  DYJets.push_back("DY1JetsToLL_M-50");
  DYJets.push_back("DY2JetsToLL_M-50");
  DYJets.push_back("DY3JetsToLL_M-50");
  DYJets.push_back("DY4JetsToLL_M-50");
  DYJets.push_back("DYJetsToLL_M-50");
  DYJets.push_back("DYJetsToLL_M-50");
  vector<TString> GluGluHToTauTau_M125;
  GluGluHToTauTau_M125.push_back("GluGluHToTauTau_M125");
  vector<TString> VBFHToTauTau_M125;
  VBFHToTauTau_M125.push_back("VBFHToTauTau_M125");
  
  map< TString , vector<TString> > samples_map = {
    { "SingleTop"            , SingleTop },
    // { "TTBar"                , TTBar },
    // { "Diboson"              , Diboson },
    // { "WJets"                , WJets  },
    // { "DYJets"               , DYJets  },
    // { "GluGluHToTauTau_M125" , GluGluHToTauTau_M125  },
    // { "VBFHToTauTau_M125"    , VBFHToTauTau_M125  }
  };

  map<TString, double> xsec_map = {
    { "DYJetsToLL_M-50"          , 5765  },
    { "DY1JetsToLL_M-50"         , -1  },
    { "DY2JetsToLL_M-50"         , -1  },
    { "DY3JetsToLL_M-50"         , -1  },
    { "DY4JetsToLL_M-50"         , -1  },
    { "DYJetsToLL_M-10to50"      , 18610 },
    { "WJetsToLNu"               , 61526.7},
    { "W1JetsToLNu"              , -1},
    { "W2JetsToLNu"              , -1},
    { "W3JetsToLNu"              , -1},
    { "W4JetsToLNu"              , -1},
    { "TTBar"                    , 831.76},
    { "ST_t-channel_antitop"     , -1 },
    { "ST_t-channel_top"         , -1 },
    { "ST_tW_antitop"            , -1 },
    { "ST_tW_top"                , -1 },
    { "VVTo2L2Nu"                , -1 },
    { "WWToLNuQQ"                , -1 },
    { "WZTo2L2Q"                 , -1 },
    { "WZTo1L1Nu2Q"              , -1 },
    { "WZTo1L3Nu"                , -1 },
    { "WZJToLLLNu"               , -1 },
    { "ZZTo4L"                   , -1 },
    { "ZZTo2L2Q"                 , -1 },
    { "GluGluHToTauTau_M125"     , -1 },
    { "VBFHToTauTau_M125"        , -1 }
  };

  double luminosity = 10000;

  for (auto const& sample : samples_map){
    
    cout << endl << sample.first << "  :  " << endl ;
    
    TFile *outFile = new TFile(sample.first + ".root","RECREATE");
    TTree *outTree = new TTree("TauCheck", "tree created as DNN input");
    
    bool firstTree = true;
    
    for(auto const& subsample: sample.second) {
      
      cout << "  - " << subsample << endl;
      
      TFile *inFile = new TFile(inputDir+"/"+subsample+"_files/"+subsample+"_1.root","READ");
      TTree *inTree = (TTree*) inFile->Get("TauCheck");
      
      double nevents = getNEventsProcessed(inFile);

      outFile->cd();

      float lumi_xsec_weight;
      if(firstTree){
	outTree = inTree->CloneTree(0);
	firstTree = false;
	TBranch *w = outTree->Branch("lumi_xsec_weight", &lumi_xsec_weight, "lumi_xsec_weight/F");
      }
      outTree->SetBranchAddress("lumi_xsec_weight",&lumi_xsec_weight);
      
      for (int i=0; i<inTree->GetEntries(); i++) {
	inTree->GetEntry(i);

	if( xsec_map.find(subsample) == xsec_map.end() ){
	  cout << endl << endl << "Sample " << subsample << " is missing in xsec_map. Exit code." << endl << endl ;
	  exit(-1);
	}

	lumi_xsec_weight = xsec_map[subsample]*luminosity/nevents;
	outTree->Fill();
      }

      outTree->AutoSave();
    
    }
  }

  cout << endl; 
}
