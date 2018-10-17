#include <iostream>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

double luminosity = 10000;

double getNEventsProcessed(TString filename)
{
  TFile * file = new TFile(filename);
  TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
  double nevents = histWeightsH->GetSumOfWeights();
  file->Close();
  delete file;
  return nevents;
}

void createDNNinput(TString inputDir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu2016/NTuples"){

  // Define the subsamples that belong to a certain proccess
  vector<TString> SingleTop       = { "ST_t-channel_antitop" , "ST_t-channel_top" , "ST_tW_antitop" , "ST_tW_antitop"};
  vector<TString> TTbar           = {"TTbar"};
  vector<TString> Diboson         = { "VVTo2L2Nu" , "WZJToLLLNu" , "WZTo1L1Nu2Q" , "WZTo1L3Nu" , "WZTo2L2Q" , "ZZTo2L2Q" , "ZZTo4L" , "WWToLNuQQ"};
  vector<TString> WJets           = { "W1JetsToLNu" , "W2JetsToLNu" , "W3JetsToLNu" , "W4JetsToLNu" , "WJetsToLNu" };
  vector<TString> DYJets          = { "DY1JetsToLL_M-50" , "DY2JetsToLL_M-50" , "DY3JetsToLL_M-50" , "DY4JetsToLL_M-50" , "DYJetsToLL_M-50" , "DYJetsToLL_M-10to50"};
  vector<TString> GluGluHToTauTau = { "GluGluHToTauTau_M125" };
  vector<TString> VBFHToTauTau    = { "VBFHToTauTau_M125" };

  // Mapping of subsamples to output root-file
  map< TString , vector<TString> > samples_map = {
    { "SingleTop"       , SingleTop },
    { "TTbar"           , TTbar },
    { "Diboson"         , Diboson },
    { "WJets"           , WJets  },
    { "DYJets"          , DYJets  },
    { "GluGluHToTauTau" , GluGluHToTauTau  },
    { "VBFHToTauTau"    , VBFHToTauTau  }
  };

  // Cross-section map
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

  // Needed for stitching
  double xsecWIncl      = xsec_map["WJetsToLNu"];
  double xsecW1Jets     = xsec_map["W1JetsToLNu"];
  double xsecW2Jets     = xsec_map["W2JetsToLNu"];
  double xsecW3Jets     = xsec_map["W3JetsToLNu"];
  double xsecW4Jets     = xsec_map["W4JetsToLNu"];
  double xsecDYIncl     = xsec_map["DYJetsToLL_M-50"];
  double xsecDY1Jets    = xsec_map["DY1JetsToLL_M-50"];
  double xsecDY2Jets    = xsec_map["DY2JetsToLL_M-50"];
  double xsecDY3Jets    = xsec_map["DY3JetsToLL_M-50"];
  double xsecDY4Jets    = xsec_map["DY4JetsToLL_M-50"];
  double neventsWIncl   = getNEventsProcessed(inputDir+"/WJetsToLNu.root");
  double neventsW1Jets  = getNEventsProcessed(inputDir+"/W1JetsToLNu.root");
  double neventsW2Jets  = getNEventsProcessed(inputDir+"/W2JetsToLNu.root");
  double neventsW3Jets  = getNEventsProcessed(inputDir+"/W3JetsToLNu.root");
  double neventsW4Jets  = getNEventsProcessed(inputDir+"/W4JetsToLNu.root");
  double neventsDYIncl  = getNEventsProcessed(inputDir+"/DYJetsToLL_M-50.root");
  double neventsDY1Jets = getNEventsProcessed(inputDir+"/DY1JetsToLL_M-50.root");
  double neventsDY2Jets = getNEventsProcessed(inputDir+"/DY2JetsToLL_M-50.root");
  double neventsDY3Jets = getNEventsProcessed(inputDir+"/DY3JetsToLL_M-50.root");
  double neventsDY4Jets = getNEventsProcessed(inputDir+"/DY4JetsToLL_M-50.root");
  
  // Loop over all samples
  for (auto const& sample : samples_map){
    
    cout << endl << sample.first << "  :  " << endl ;
    
    TFile *outFile = new TFile(sample.first + ".root","RECREATE");
    TTree *outTree = new TTree("TauCheck", "tree created as DNN input");    
    bool firstTree = true;

    for(auto const& subsample: sample.second) {
      
      cout << "  - " << subsample << endl;
      
      TFile *inFile  = new TFile( inputDir + "/" + subsample + ".root" ,"READ");
      TTree *inTree  = (TTree*) inFile -> Get("TauCheck");
      double nevents = getNEventsProcessed( inputDir + "/" + subsample + ".root" );
      outFile->cd();

      // SetBranchAddress for variables that need to be added and that are needed for preselection
      float lumi_xsec_weight;
      float iso_1;
      float iso_2;
      bool extraelec_veto;
      bool extramuon_veto;
      float pt_1;
      float pt_2;
      bool metFilters;
      bool trg_muonelectron;
      unsigned int npartons;
      if(firstTree){
	outTree    = inTree->CloneTree(0);
	TBranch *w = outTree->Branch("lumi_xsec_weight", &lumi_xsec_weight, "lumi_xsec_weight/F");
	firstTree  = false;
      }
      outTree->SetBranchAddress("lumi_xsec_weight",&lumi_xsec_weight);
      outTree->SetBranchAddress("npartons",&npartons);
      outTree->SetBranchAddress("iso_1",&iso_1);
      outTree->SetBranchAddress("iso_2",&iso_2);
      outTree->SetBranchAddress("extraelec_veto",&extraelec_veto);
      outTree->SetBranchAddress("extramuon_veto",&extramuon_veto);
      outTree->SetBranchAddress("pt_1",&pt_1);
      outTree->SetBranchAddress("pt_2",&pt_2);
      outTree->SetBranchAddress("metFilters",&metFilters);
      outTree->SetBranchAddress("trg_muonelectron",&trg_muonelectron);

      for (int i=0; i<inTree->GetEntries(); i++) {

	inTree->GetEntry(i);

	// Add here preselection if necessary

	if( xsec_map.find(subsample) == xsec_map.end() ){
	  cout << endl << endl << "Sample " << subsample << " is missing in xsec_map. Exit code." << endl << endl ;
	  exit(-1);
	}

	// lumi-xsec-weight added
	lumi_xsec_weight = xsec_map[subsample]*luminosity/nevents;
	
	// Stitching only for wjets MC in n-jet binned samples in npartons
	if( sample.first.Contains("WJets") ){
	  if(npartons == 1)      lumi_xsec_weight = luminosity / ( neventsW1Jets/xsecW1Jets + neventsWIncl/xsecWIncl );
	  else if(npartons == 2) lumi_xsec_weight = luminosity / ( neventsW2Jets/xsecW2Jets + neventsWIncl/xsecWIncl );
	  else if(npartons == 3) lumi_xsec_weight = luminosity / ( neventsW3Jets/xsecW3Jets + neventsWIncl/xsecWIncl );
	  else if(npartons == 4) lumi_xsec_weight = luminosity / ( neventsW4Jets/xsecW4Jets + neventsWIncl/xsecWIncl );
	  else                   lumi_xsec_weight = luminosity / ( neventsWIncl/xsecWIncl );
	}
	else if( sample.first.Contains("DYJets") ){
	  if(npartons == 1)      lumi_xsec_weight = luminosity / ( neventsDY1Jets/xsecDY1Jets + neventsDYIncl/xsecDYIncl );
	  else if(npartons == 2) lumi_xsec_weight = luminosity / ( neventsDY2Jets/xsecDY2Jets + neventsDYIncl/xsecDYIncl );
	  else if(npartons == 3) lumi_xsec_weight = luminosity / ( neventsDY3Jets/xsecDY3Jets + neventsDYIncl/xsecDYIncl );
	  else if(npartons == 4) lumi_xsec_weight = luminosity / ( neventsDY4Jets/xsecDY4Jets + neventsDYIncl/xsecDYIncl );
	  else                   lumi_xsec_weight = luminosity / ( neventsDYIncl/xsecDYIncl );
	}

	outTree->Fill();
      }

      outTree->AutoSave();

    }
  }

  cout << endl; 
}
