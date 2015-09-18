#include "tdrstyle.C"
#include "HtoH.h"

Double_t TurnOn(Double_t * x, Double_t * par) {

  return 0.5*par[0]*(TMath::Erf((x[0]-par[1])/par[2])+1);

}

// .x PlotTrigEffGenPt.C("DYJets","TTJets","ElectronPt23",0.4,1);
void PlotTrigEffPt(TString fileName1 = "DYJets_TP",
		   TString fileName2 = "TTJets_TP",
		   TString TrigName = "ElectronPt23",
		   float yMin = 0.,
		   float yMax = 1.0) {

  // inputs
  // fileName1  : filename of DYJets root file produced by AnalysisMacroEff (w/o .root extension)
  // fileName2  : filename of TTJets root file produced by AnalysisMacroEff (w/o .root extension)
  // TrigName   : name of leg, possible options : MuonPt8, MuonPt23, ElectronPt12, ElectronPt23
  // yMin, yMax : ranges of y axis 

  TString yTitle("Efficiency  (Ele12 Leg)");
  if (TrigName.Contains("ElectronPt23"))
    yTitle = "Efficiency  (Ele23 Leg)";
  if (TrigName.Contains("MuonPt23"))
    yTitle = "Efficiency  (Mu23 Leg)";
  if (TrigName.Contains("MuonPt8"))
    yTitle = "Efficiency  (Mu8 Leg)";

  setTDRStyle();
  //  gStyle->SetOptStat(0000);

  TFile * file1 = new TFile(fileName1+".root");
  TFile * file2 = new TFile(fileName2+".root");

  int nPtBins = 12;
  float ptBins[13] = {6,8,10,12,14,16,18,20,25,30,40,50,100};
  if (TrigName.Contains("Pt23")) {
    nPtBins = 13;
    ptBins[0] = 10;
    ptBins[1] = 15;
    ptBins[2] = 20;
    ptBins[3] = 21;
    ptBins[4] = 22;
    ptBins[5] = 23;
    ptBins[6] = 24;
    ptBins[7] = 25;
    ptBins[8] = 26;
    ptBins[9] = 28;
    ptBins[10] = 30;
    ptBins[11] = 40;
    ptBins[12] = 50;
    ptBins[13] = 100;
  }

  TString numName1 = TrigName + "RecoPassedH";
  TString denName1 = TrigName + "RecoH";

  TString numName2 = TrigName + "RecoPassedXH";
  TString denName2 = TrigName + "RecoXH";

  TH1F * histNum1  = (TH1F*)file1->Get(numName1);
  TH1F * histNumX1 = TH1toTH1(histNum1,nPtBins,ptBins,false,"_new");
  TH1F * histDen1  = (TH1F*)file1->Get(denName1);
  TH1F * histDenX1 = TH1toTH1(histDen1,nPtBins,ptBins,false,"_new");

  TH1F * histNum2  = (TH1F*)file2->Get(numName2);
  TH1F * histNumX2 = TH1toTH1(histNum2,nPtBins,ptBins,false,"_new");
  TH1F * histDen2  = (TH1F*)file2->Get(denName2);
  TH1F * histDenX2 = TH1toTH1(histDen2,nPtBins,ptBins,false,"_new");
  TH1F * inputEventsH = (TH1F*)file1->Get("inputEventsH");

  TGraphAsymmErrors * eff1 = new TGraphAsymmErrors();
  eff1->Divide(histNumX1,histDenX1);
  eff1->SetMarkerSize(1.5);
  eff1->SetLineWidth(2);
 
  TGraphAsymmErrors * eff2 = new TGraphAsymmErrors();
  eff2->Divide(histNumX2,histDenX2);
  eff2->SetMarkerSize(2);
  eff2->SetLineWidth(2);
  eff2->SetMarkerStyle(21);
  eff2->SetMarkerColor(2);
  eff2->SetLineColor(2);

  TH2F * frame = new TH2F("frame","",2,10,50,2,yMin,yMax);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
  if (TrigName.Contains("Electron"))
    frame->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
  TCanvas * canv = new TCanvas("canv","",900,700);
  frame->Draw();
  TF1 * fitFunc = new TF1("fitFunc",TurnOn,9,50,3);
  fitFunc->SetParameter(0,0.9);
  fitFunc->SetParameter(1,23);
  fitFunc->SetParameter(2,2);
  fitFunc->SetLineWidth(2);
  //  eff2->Fit("fitFunc","LR");
  //  fitFunc->Draw("SAME");
  eff2->Draw("PESAME");
  eff1->Draw("PESAME");

  TLegend * leg = new TLegend(0.6,0.2,0.92,0.4);
  leg->SetFillColor(0);
  leg->AddEntry(eff1,"Z#rightarrow#tau#tau (e+#mu)","lp");
  leg->AddEntry(eff2,"t#bar{t} (e+#mu)","lp");
  leg->Draw();
  canv->SetGridx();
  canv->SetGridy();
  
  canv->Update();
  canv->Print(fileName1+"_"+TrigName+".png");

}
