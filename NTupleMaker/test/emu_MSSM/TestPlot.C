#include "HttStylesNew.cc"
void TestPlot(
	      //	      TString fileName = "2018/EmbeddedElMu_Run2018A.root",
	      TString fileName = "2018/SUSYGluGluToHToTauTau_M-1200.root",
	      //	      TString var = "trg_singlemuon||trg_singleelectron",
	      TString var = "trg_muhigh_elow||trg_ehigh_mulow",
	      int Nbins = 2,
	      float xmin = -0.5,
	      float xmax = 1.5) {

  SetStyle();

  TFile * file = new TFile(fileName);
  TTree * tree = (TTree*)file->Get("TauCheck");
  TH1D * hist = new TH1D("hist","",Nbins,xmin,xmax);
  TCanvas * canv = new TCanvas("canv","",600,600);
  tree->Draw(var+">>hist","iso_1<0.15&&iso_2<0.2");
  canv->Update();

}
