#include "HttStylesNew.cc"
void PlotEff(TString dataFileName = "DataEE_eleEndcap",
	     TString mcFileName = "DYJetsToEE_M-50_eleEndcap",
	     TString histBaseName = "ZMassEndcap") {
 
  SetStyle();

  TFile * fileData = new TFile(dataFileName+".root");
  TFile * fileMC = new TFile(mcFileName+".root");

  TGraphAsymmErrors * effData = (TGraphAsymmErrors*)fileData->Get(histBaseName);
  TGraphAsymmErrors * effMC = (TGraphAsymmErrors*)fileMC->Get(histBaseName);

  effData->SetLineColor(2);
  effData->SetMarkerColor(2);
  effData->SetMarkerSize(1.2);
  effData->SetMarkerStyle(20);

  effMC->SetLineColor(4);
  effMC->SetMarkerColor(4);
  effMC->SetMarkerStyle(21);
  effMC->GetXaxis()->SetTitleOffset(1.1);
  effMC->GetXaxis()->SetTitleSize(0.05);
  effMC->GetYaxis()->SetTitleOffset(1.1);
  effMC->GetYaxis()->SetTitleSize(0.05);
  effMC->GetYaxis()->SetRangeUser(0.0,1.001);
  TString EtaRegion("Endcap");
  if (histBaseName.Contains("Barrel"))
    EtaRegion = "Barrel";

  TCanvas * canv = new TCanvas("canv","",700,600);
  effMC->Draw("APE");
  effData->Draw("PESame");
  TLegend * leg = new TLegend(0.69,0.17,0.94,0.40);
  leg->SetFillColor(0);
  leg->SetHeader(EtaRegion);
  leg->SetTextSize(0.06);
  leg->AddEntry(effData,"Data","lp");
  leg->AddEntry(effMC,"MC","lp");
  leg->Draw();
  canv->SetGridx();
  canv->SetGridy();
  canv->Print(histBaseName+"_eff.png");
  canv->Update();


}
