#include "HttStylesNew.cc"
void TestTrackCollection() {

  SetStyle();
  
  //  TString fileName = "output_DATA.root";
  //  TString fileName = "/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/mc/SUSYGluGluToBBHToTauTau_M-900/SUSYGluGluToBBHToTauTau_M-900_1857.root";
  TString fileName = "/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2018/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_8077.root";

  TFile * file = new TFile(fileName);
  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");
  tree->Draw("track_count");


}
