#include "tdrstyle.C"
#include "HtoH.h"
// .x PlotEffPt.C("DYJets","RecoPassedTag","RecoAllTag","RecoPassedTag","RecoAllTag",1,"NonPrompt","Prompt",2,0.4,1.0,"Id+Iso efficiency")
void PlotEffPt(TString fileName = "DYJets_TP",
	       TString numName0 = "RecoPassedTag",
	       TString denName0 = "RecoAllTag",
	       TString numName1 = "RecoPassedTag",
	       TString denName1 = "RecoAllTag",
	       bool isMuId = true,
	       TString namePrompt0 = "NonPrompt",
	       TString namePrompt1 = "Prompt",
	       int EtaBin = 2,
	       float yMin = 0.4,
	       float yMax = 1.0,
	       TString yTitle="Id+Iso efficiency") {
  // inputs
  // fileName - the name of the RooT file produced by AnalysisMacroEff (without .root extension)
  // numName0 - name of the numerator histogram for non-prompt leptons
  // denName0 - name of the denominator histogram for non-prompt leptons
  // numName1 - name of the numerator histogram for prompt leptons
  // denName1 - name of the denominator histogram for prompt leptons
  // isMuId - false : electron , true : muon
  // namePrompt0 - name of the histogram for prompt leptons
  // namePrompt1 - name of the histogram for non-prompt leptons
  // EtaBin - eta bin  0 : 0-0.8,  1 : 0.8-1.5, 2 : 1.5-2.4
  // yMin, yMax - ranges of the y axis (efficiency)
  // yTitle : title of y axis

  setTDRStyle();

  int lepSign = 0;

  TString name("ElectronPt_");
  if (isMuId)
    name = "MuonPt_";

  TString LepQ("Neg");
  if (lepSign>0)
    LepQ = "Pos";

  TString Type[2];
  Type[0] = namePrompt0;
  Type[1] = namePrompt1;


  int nEtaBins = 3;
  TString EtaBinName[3] = {"0To0p8",
			   "0p8To1p5",
			   "1p5To2p3"};

  TString EtaLeg[3] = {"#eta = [0,0.8]","#eta = [0.8,1.5]","#eta = [1.5,2.5]"};
  if (isMuId)
    EtaLeg[2] = "#eta = [1.5,2.4]";

  TFile * file = new TFile(fileName+".root");

  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};
  TString numName[2];
  TString denName[2];
  numName[0] = numName0;
  numName[1] = numName1;
  denName[0] = denName0;
  denName[1] = denName1;  

  TH1F * histos[3][2][2];
  for (int iEta=0; iEta<3; ++iEta) {
    for (int iType=0; iType<2; ++iType) {
      TH1F * histNum = (TH1F*)file->Get(name+LepQ+Type[iType]+numName[iType]+EtaBinName[iEta]);
      if (lepSign==0) {
	TH1F * histNumX = (TH1F*)file->Get(name+"Pos"+Type[iType]+numName[iType]+EtaBinName[iEta]);
	histNum->Add(histNum,histNumX);
      }
      histos[iEta][iType][0] = TH1toTH1(histNum,nPtBins,ptBins,false,"_new");
      TH1F * histDen = (TH1F*)file->Get(name+LepQ+Type[iType]+denName[iType]+EtaBinName[iEta]);
      if (lepSign==0) {
	TH1F * histDenX = (TH1F*)file->Get(name+"Pos"+Type[iType]+denName[iType]+EtaBinName[iEta]);
	histDen->Add(histDen,histDenX);
      }
      histos[iEta][iType][1] = TH1toTH1(histDen,nPtBins,ptBins,false,"_new");
    }
  }

  int color[3] = {1,2,4};
  int symbol[3] = {20,21,22};

  TGraphAsymmErrors * eff[3][2];
  for (int iEta=0; iEta<3; ++iEta) {
    for (int iType=0; iType<2; ++iType) {
      eff[iEta][iType] = new TGraphAsymmErrors();
      eff[iEta][iType]->SetLineColor(color[iEta]);
      eff[iEta][iType]->SetMarkerColor(color[iEta]);
      eff[iEta][iType]->SetLineWidth(2);
      eff[iEta][iType]->SetMarkerSize(2.5);
      eff[iEta][iType]->SetMarkerStyle(symbol[iEta]);
      if (iType==1)
	eff[iEta][iType]->SetMarkerStyle(symbol[iEta]+4);
      eff[iEta][iType]->Divide(histos[iEta][iType][0],histos[iEta][iType][1]);
    }
  }  
 
  // **************
  // Tag-and-Probe
  // **************

  TString mass("massEEId_");
  if (isMuId) 
    mass = "massMuMuId_";

  TString TagQ("Pos");

  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};
  TString PtBinName[7] = {"pt10to15",
			  "pt15to20",
			  "pt20to30",
			  "pt30to40",
			  "pt40to50",
			  "pt50to70",
			  "pt70to100"};
  
  TString passName[2] = {"Pass","Fail"};

  TH1F * histosTP[3][2];
  for (int iEta=0; iEta<3; ++iEta) {
    for (int iPass=0; iPass<2; ++iPass) {
      histosTP[iEta][iPass] = new TH1F(passName[iPass]+EtaBinName[iEta],"",nPtBins,ptBins);
      for (int iPt=0; iPt<nPtBins; ++iPt) {
	TH1F * hist = (TH1F*)file->Get(mass+TagQ+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass]);
	float yield = hist->GetSum();
	histosTP[iEta][iPass]->SetBinContent(iPt+1,yield);
      }
    }
  }
  
  TGraphAsymmErrors * effTP[3];
  for (int iEta=0; iEta<3; ++iEta) {
    effTP[iEta] = new TGraphAsymmErrors();
    effTP[iEta]->SetLineColor(color[iEta]);
    effTP[iEta]->SetMarkerColor(color[iEta]);
    effTP[iEta]->SetLineWidth(2);
    effTP[iEta]->SetMarkerSize(3);
    effTP[iEta]->SetMarkerStyle(28);
    histosTP[iEta][1]->Add(histosTP[iEta][0],histosTP[iEta][1]);
    effTP[iEta]->Divide(histosTP[iEta][0],histosTP[iEta][1]);
  }  

  TH2F * frame = new TH2F("frame","",2,10,50,2,yMin,yMax);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
  //  frame->GetYaxis()->SetNdivisions(505);
  TString LegName("Z#rightarrowee");
  if (isMuId) {
    frame->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
    LegName = "Z#rightarrow#mu#mu";
  }
  TCanvas * canv = new TCanvas("canv","",900,700);
  frame->Draw();
  eff[EtaBin][0]->Draw("EPS");
  eff[EtaBin][1]->Draw("EPS");
  //  effTP[EtaBin]->Draw("EPS");

  TLegend * leg = new TLegend(0.6,0.2,0.9,0.5);
  leg->SetFillColor(0);
  leg->SetHeader(EtaLeg[EtaBin]);
  leg->AddEntry(eff[EtaBin][0],"Z#rightarrow#tau#tau#rightarrowe+#mu","lp");
  leg->AddEntry(eff[EtaBin][1],LegName,"lp");
  //  leg->AddEntry(effTP[EtaBin] ,LegName+" (T&P)","lp");
  leg->Draw();
  canv->Update();
  canv->Print(fileName+"_"+name+numName0+"_"+denName0+"_"+EtaBinName[EtaBin]+".gif");

}
