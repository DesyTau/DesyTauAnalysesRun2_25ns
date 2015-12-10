#include "HttStylesNew.cc"
void PlotLeptonEff(TString effFileName = "Muon_IdIso0p10_eff",
		   TString histBaseName = "ZMassEtaGt1p2",
		   TString xtit = "muon p_{T} [GeV]",
		   TString header = "|#eta|>1.2") {
 
  SetStyle();
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);

  TFile * file = new TFile(effFileName+".root");

  TGraphAsymmErrors * effData = (TGraphAsymmErrors*)file->Get(histBaseName+"_Data");
  TGraphAsymmErrors * effMC = (TGraphAsymmErrors*)file->Get(histBaseName+"_MC");

  effData->SetLineColor(2);
  effData->SetMarkerColor(2);

  effMC->SetLineColor(4);
  effMC->SetMarkerColor(4);
  effMC->SetMarkerStyle(21);
  effMC->GetYaxis()->SetRangeUser(0.0,1.02);

  effMC->GetXaxis()->SetTitle(xtit);

  TCanvas * canv = new TCanvas("canv","",700,600);
  effMC->Draw("APE");
  effData->Draw("PESame");
  TLegend * leg = new TLegend(0.6,0.17,0.94,0.40);
  leg->SetFillColor(0);
  leg->SetHeader(header);
  leg->SetTextSize(0.06);
  leg->AddEntry(effData,"Data","lp");
  leg->AddEntry(effMC,"Z MC","lp");
  leg->Draw();
  canv->SetGridx();
  canv->SetGridy();
  canv->Print("figures/"+histBaseName+"_eff.png");
  canv->Update();


}
