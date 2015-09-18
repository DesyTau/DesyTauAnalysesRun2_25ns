#include "HttStylesNew.cc"
#include "HtoH.h"

void Plot(TString histName = "dileptonMassSelH",
	  TString xtitle = "M_{e#mu} [GeV]",
	  TString ytitle = "Entries / 5 GeV",
	  float xLower = 0,
	  float xUpper = 200,
	  float yLower = 0,
	  float yUpper = -1,
	  bool logY = false) {

  SetStyle();
  
  TString samples[4] = {"QCD_MuEnriched","WJets","TTJets","DYJets"};
  float xsec[4] = {866600000*0.00044,50100,424,4746};

  float lumi = 5000; // 5 fb-1

  TFile * file = new TFile("DYJets.root");

  TH1F * histOld = (TH1F*)file->Get(histName);

  int nBins = histOld->GetNbinsX();
  float xMin = histOld->GetBinLowEdge(1);
  float xMax = histOld->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  float bins[100];
  int nBinsNew = nBins;
  std::cout << "New number of bins : ";
  std::cin >> nBinsNew;

  if (nBins % nBinsNew >0) { 
    std::cout << "new number of bins = " << nBinsNew 
	      << "  not multiple of " << nBins << std::endl;
    return;
  }
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB)
    bins[iB] = xMin + float(iB)*binWidth;

  TH1F * qcdHist = new TH1F("qcdHist","",nBinsNew,bins);
  TH1F * wHist  = new TH1F("wHist","",nBinsNew,bins);
  TH1F * ttHist = new TH1F("ttHist","",nBinsNew,bins);
  TH1F * zHist = new TH1F("zHist","",nBinsNew,bins);

  int nSamples = 4;

  //  return;

  for (int iS=0; iS<nSamples; ++iS) {
    //    std::cout << "Sample = " << iS << std::endl;
    TFile * fileMC = new TFile(samples[iS]+".root");
    TH1F * histOld = (TH1F*)fileMC->Get(histName);
    TH1F * hist = TH1toTH1(histOld,nBinsNew,bins,true,"_new_"+samples[iS]);
    TH1F * eventCount = (TH1F*)fileMC->Get("inputEventsH");
    float nGen = eventCount->GetSum();
    float norm = xsec[iS]*lumi/nGen;
    TH1F * tempHist = qcdHist;
    if (iS==1)
      tempHist = wHist;
    if (iS==2)
      tempHist = ttHist;
    if (iS==3)
      tempHist = zHist;
    tempHist->Add(tempHist,hist,1.,norm);
  }

  //  float dataEvents = 0;
  //  float ttEvents = 0;

  wHist->Add(wHist,qcdHist);
  ttHist->Add(ttHist,wHist);
  zHist->Add(zHist,ttHist);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    wHist->SetBinError(iB,0);
    ttHist->SetBinError(iB,0);
    qcdHist->SetBinError(iB,0);
    zHist->SetBinError(iB,0);
  }

  InitHist(qcdHist,"","",TColor::GetColor("#ffccff"),1001);
  InitHist(ttHist ,"","",TColor::GetColor("#9999cc"),1001);
  InitHist(wHist  ,"","",TColor::GetColor("#de5a6a"),1001);
  InitHist(zHist  ,"","",TColor::GetColor("#ffcc66"),1001);

  TCanvas * canv1 = MakeCanvas("canv1", "", 700, 700);

  zHist->GetXaxis()->SetTitle(xtitle);
  zHist->GetYaxis()->SetTitle(ytitle);

  zHist->Draw("h");
  ttHist->Draw("sameh");
  wHist->Draw("sameh");
  qcdHist->Draw("sameh");

  TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
  SetLegendStyle(leg);
  leg->AddEntry(zHist,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(ttHist,"t#bar{t}","f");
  leg->AddEntry(wHist,"W+Jets","f");
  leg->AddEntry(qcdHist,"QCD multijets","f");
  leg->Draw();

  TPad * pad = canv1->GetPad(0);
  pad->RedrawAxis();
  if (logY)
    canv1->SetLogy();

  canv1->Update();

}
