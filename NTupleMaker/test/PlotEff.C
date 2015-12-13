#include "HttStylesNew.cc"
void PlotEff(TString fileName1 = "SingleMuon_2015D_iso0p1",
	     TString fileName2 = "SingleMuon_2015D_iso0p15",
	     bool isData1 = true,
	     bool isData2 = true,
	     TString leg1 = "relIso<0.1",
	     TString leg2 = "relIso<0.15",
	     TString histBaseName = "ZMassEta0p9to1p2",
	     TString xtit = "muon p_{T} [GeV]",
	     TString header = "  Data, 0.9<|#eta|<1.2") {
 
  SetStyle();
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);

  TFile * file1 = new TFile(fileName1+"_"+histBaseName+".root");
  TFile * file2 = new TFile(fileName2+"_"+histBaseName+".root");

  TString postfix1("_MC");
  if (isData1) postfix1 = "_Data";
  TString postfix2("_MC");
  if (isData2) postfix2 = "_Data";

  TGraphAsymmErrors * eff1 = (TGraphAsymmErrors*)file1->Get(histBaseName+postfix1);
  TGraphAsymmErrors * eff2 = (TGraphAsymmErrors*)file2->Get(histBaseName+postfix2);

  eff1->SetLineColor(2);
  eff1->SetMarkerColor(2);

  eff2->SetLineColor(4);
  eff2->SetMarkerColor(4);
  eff2->SetMarkerStyle(21);
  eff2->GetYaxis()->SetRangeUser(0.0,1.02);

  eff2->GetXaxis()->SetTitle(xtit);

  TCanvas * canv = new TCanvas("canv","",700,600);
  eff2->Draw("APE");
  eff1->Draw("PESame");
  TLegend * leg = new TLegend(0.55,0.17,0.94,0.40);
  leg->SetFillColor(0);
  leg->SetHeader(header);
  leg->SetTextSize(0.06);
  leg->AddEntry(eff1,leg1,"lp");
  leg->AddEntry(eff2,leg2,"lp");
  leg->Draw();
  canv->SetGridx();
  canv->SetGridy();
  canv->Print("figures/"+histBaseName+"_eff.png");
  canv->Update();


}
