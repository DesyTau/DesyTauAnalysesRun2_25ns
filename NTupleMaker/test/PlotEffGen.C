#include "tdrstyle.C"
#include "HtoH.h"
void PlotEffGen(TString fileName = "DYJetsEff",
		// TString numName = "TPPassed",
		// TString denName = "TPAll",
		//		TString numName = "RecoPassedTag",
		//		TString denName = "RecoAllTag",
		TString numName = "RecoPassed",
		TString denName = "RecoAll",
		bool posTag = true,
		bool isMuId = true,
		bool isPrompt = true) {


  setTDRStyle();

  TString name("ElectronPt_");
  if (isMuId)
    name = "MuonPt_";

  TString LepQ("Neg");
  if (posTag)
    LepQ = "Pos";

  TString Type("NonPrompt");
  if (isPrompt)
    Type = "Prompt";

  int nEtaBins = 3;
  TString EtaBinName[3] = {"0To0p8",
			   "0p8To1p5",
			   "1p5To2p3"};
  TFile * file = new TFile(fileName+".root");

  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};

  TH1F * histos[3][2];
  for (int iEta=0; iEta<3; ++iEta) {
    TH1F * histNum = (TH1F*)file->Get(name+LepQ+Type+numName+EtaBinName[iEta]);
    histos[iEta][0] = TH1toTH1(histNum,nPtBins,ptBins,true,"_new");
    TH1F * histDen = (TH1F*)file->Get(name+LepQ+Type+denName+EtaBinName[iEta]);
    histos[iEta][1] = TH1toTH1(histDen,nPtBins,ptBins,true,"_new");
    
  }

  int color[3] = {1,2,4};
  int symbol[3] = {20,21,22};

  TGraphAsymmErrors * eff[3];
  for (int iEta=0; iEta<3; ++iEta) {
    eff[iEta] = new TGraphAsymmErrors();
    eff[iEta]->SetLineColor(color[iEta]);
    eff[iEta]->SetMarkerColor(color[iEta]);
    eff[iEta]->SetLineWidth(2);
    eff[iEta]->SetMarkerSize(1.6);
    eff[iEta]->SetMarkerStyle(symbol[iEta]);
    eff[iEta]->Divide(histos[iEta][0],histos[iEta][1],"Pois");
  }  

 
  TH2F * frame = new TH2F("frame","",2,0,100,2,0,1);
  TCanvas * canv = new TCanvas("canv","",700,700);
  frame->Draw();
  eff[0]->Draw("PS");
  eff[1]->Draw("PS");
  eff[2]->Draw("PS");
  canv->Update();

}
