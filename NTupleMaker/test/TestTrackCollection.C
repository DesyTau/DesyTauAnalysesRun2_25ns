#include "HttStylesNew.cc"
void TestTrackCollection() {

  SetStyle();
  
  //  TString fileName = "/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/data/DoubleMuon_Run2018D/DoubleMuon_Run2018D_614.root";
  TString fileName = "output_MC.root";

  TFile * file = new TFile(fileName);
  TTree * tree = (TTree*)file->Get("makeroottree/AC1B");
  TCanvas * canv = new TCanvas("canv","");
  tree->Draw("track_pt","track_pt<30.");
  canv->Update();

}
