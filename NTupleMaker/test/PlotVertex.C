#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotVertex(TString dataFileName = "SingleMuon_Run2015D",
		TString mcFileName = "DYJetsToMuMu_M-50_13TeV-amcatnloFXFX-pythia8",
		TString histName = "NumberOfVerticesH",
		TString legMC = "Z#rightarrow #mu#mu") {

  SetStyle();

  TFile * dataFile = new TFile(dataFileName+".root");
  TFile * mcFile = new TFile(mcFileName+".root");

  TH1F * dataHist = (TH1F*)dataFile->Get(histName);
  TH1F * mcHist = (TH1F*)mcFile->Get(histName);
 
  float normData = 1/dataHist->GetSumOfWeights();
  float normMC = 1/mcHist->GetSumOfWeights();

  int nBins = dataHist->GetNbinsX();
  for (int iB=1; iB<=nBins; ++iB) {
    float x = dataHist->GetBinContent(iB);
    float e = dataHist->GetBinError(iB);
    dataHist->SetBinContent(iB,x*normData);
    dataHist->SetBinError(iB,e*normData);
    x = mcHist->GetBinContent(iB);
    mcHist->SetBinContent(iB,x*normMC);
    mcHist->SetBinError(iB,0);
  }


  dataHist->SetLineColor(2);
  dataHist->SetLineWidth(2);
  dataHist->SetMarkerColor(2);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(1.2);
  mcHist->SetLineColor(4);
  mcHist->SetLineWidth(2);

  dataHist->GetXaxis()->SetRangeUser(0,25);
  dataHist->GetXaxis()->SetTitle("N of primary vertices");
  dataHist->GetYaxis()->SetTitle("a.u.");
  dataHist->GetXaxis()->SetTitleOffset(1.);
  dataHist->GetYaxis()->SetTitleOffset(1.3);
  dataHist->GetXaxis()->SetTitleSize(0.055);
  dataHist->GetYaxis()->SetTitleSize(0.055);
  

  TCanvas * canv = new TCanvas("canv","",600,600);

  dataHist->Draw("e1");
  mcHist->Draw("same");
  dataHist->Draw("e1same");

  TLegend * leg = new TLegend(0.6,0.7,0.92,0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(dataHist,"Data","lp");
  leg->AddEntry(mcHist,legMC+" (MC)","lp");
  leg->Draw();
 

  canv->Update();
  canv->Print(dataFileName+"_"+histName+".png");

}
