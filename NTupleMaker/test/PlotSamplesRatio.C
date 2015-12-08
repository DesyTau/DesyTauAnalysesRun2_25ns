#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotSamplesRatio(TString histName = "puppiMetSelH",
		      TString xtitle = "puppi E_{T}^{mis} [GeV]",
		      TString ytitle = "Events",
		      float xLower = 0,
		      float xUpper = 400,
		      bool logY = true) {

  SetStyle();

  TFile * file = new TFile("SingleMuon_2015D.root");
  
  TString samples[16] = {"VVTo2L2Nu", // (0)
			 "WWToLNuQQ", // (1) 
			 "WZTo1L1Nu2Q", // (2)
			 "WZTo1L3Nu", // (3)
			 "WZTo2L2Q", // (5) 
			 "WZJets", // (5)
			 "ZZTo2L2Q", // (6)
			 "ZZTo4L", // (7)
			 "WJetsToLNu_MG", // (8)
			 "ST_t-channel_antitop_4f_leptonDecays", // (9)
			 "ST_t-channel_top_4f_leptonDecays", // (10)
			 "ST_tW_antitop_5f_inclusiveDecays", // (11)
			 "ST_tW_top_5f_inclusiveDecays", // (12)
			 "TTPowHeg", // (13)
			 "DYJetsToLL_M-10to50", // (14)
			 "DYJetsToLL_M-50" // (15)
  };

  float xsec[16] = {11.95,   // VVTo2L2Nu (0)
		   49.997,   // WWToLNuQQ (1)
		    10.71,   // WZTo1L1Nu2Q (2)
		    3.05,    // WZTo1L3Nu (3)
		    5.595,   // WZTo2L2Q (4)
		    5.26,    // WZJets (3L11Nu) (5) 
		    3.22,    // ZZTo2L2Q (6)
		    1.212,   // ZZTo4L (7)
		    61526,   // WJets (8)
		    80.95,   // ST_t-channel_antitop_4f_leptonDecays (9)
		    136.95,  // ST_t-channel_top_4f_leptonDecays (10)
		    35.6,    // ST_tW_antitop_5f_inclusiveDecays (11)
		    35.6,    // ST_tW_top_5f_inclusiveDecays (12)
		    831.8,   // TTbarPowHeg (13)
		    //		    71310,   // DYJetsToLL_M-5to50 (14) 
		    18610, // DYJetsToLL_M-10to50 (14)
		    6025   // DYJetsToLL_M-50_MG (15)
  };

  float lumi = 2090;

  TH1F * histDataOld = (TH1F*)file->Get(histName);

  int nBins = histDataOld->GetNbinsX();
  float xMin = histDataOld->GetBinLowEdge(1);
  float xMax = histDataOld->GetBinLowEdge(nBins+1);

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

  TH1F * histData = TH1toTH1(histDataOld,nBinsNew,bins,true,"_Data_new");

  TH1F * ewkHist = new TH1F("ewkHist","",nBinsNew,bins);
  TH1F * ttHist  = new TH1F("ttHist","",nBinsNew,bins);
  TH1F * wHist = new TH1F("wHist","",nBinsNew,bins);
  TH1F * zHist = new TH1F("zHist","",nBinsNew,bins);

  int nSamples = 16;

  //  return;

  for (int iS=0; iS<nSamples; ++iS) {
    //    std::cout << "Sample = " << iS << std::endl;
    TFile * fileMC = new TFile(samples[iS]+".root");
    TH1F * histOld = (TH1F*)fileMC->Get(histName);
    TH1F * hist = TH1toTH1(histOld,nBinsNew,bins,true,"_new_"+samples[iS]);
    TH1F * eventCount = (TH1F*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[iS]*lumi/nGen;
    TH1F * tempHist = ewkHist;
    if (iS==8)
      tempHist = wHist;
    if (iS>8&&iS<14)
      tempHist = ttHist;
    if (iS==14||iS==15)
      tempHist = zHist;
    tempHist->Add(tempHist,hist,1.,norm);

  }

  //  float dataEvents = 0;
  //  float ttEvents = 0;

  wHist->Add(wHist,ewkHist);
  ttHist->Add(ttHist,wHist);
  zHist->Add(zHist,ttHist);

  TH1F * bkgdErr = (TH1F*)zHist->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    ewkHist->SetBinError(iB,0);
    ttHist->SetBinError(iB,0);
    wHist->SetBinError(iB,0);
    zHist->SetBinError(iB,0);
  }

  InitData(histData);
  InitHist(wHist,"","",kMagenta,1001);
  InitHist(ttHist,"","",kCyan,1001);
  InitHist(ewkHist,"","",kBlue-4,1001);
  InitHist(zHist,"","",kYellow,1001);
  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  histData->GetXaxis()->SetRangeUser(xLower,xUpper);
  float yUpper = histData->GetMaximum();
  if (logY)
    histData->GetYaxis()->SetRangeUser(0.5,2*yUpper);
  else
    histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);

  histData->SetMarkerSize(1.5);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  //  nData = histData->GetSum();
  //  float nMC   = ttHist->GetSum();
  //  float eData = TMath::Sqrt(nData);


  TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);

  histData->Draw("e1");
  zHist->Draw("sameh");
  ttHist->Draw("sameh");
  wHist->Draw("sameh");
  ewkHist->Draw("sameh");
  histData->Draw("e1same");

  float chi2 = 0;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = zHist->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.65,0.6,0.9,0.88);
  SetLegendStyle(leg);
  leg->SetTextSize(0.05);
  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(ewkHist,"dibosons","f");
  leg->AddEntry(wHist,"W+Jets","f");
  leg->AddEntry(ttHist,"t#bar{t}+single top","f");
  leg->AddEntry(zHist,"Z#rightarrow#mu#mu","f");
  leg->Draw();

  TLatex * cms = new TLatex(0.25,0.94,"CMS Preliminary   L = 41 pb^{-1} at #sqrt{s} = 13 TeV");

  cms->SetNDC();
  cms->SetTextSize(0.05);
  cms->Draw();

  if (logY) upper->SetLogy(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1F * ratioH = (TH1F*)histData->Clone("ratioH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.5);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0,2.2);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.5);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = zHist->GetBinContent(iB);
    if (x1>0&&x2>0) {
      float e1 = histData->GetBinError(iB);
      float ratio = x1/x2;
      float eratio = e1/x2;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }


  // ------------>Primitives in pad: lower
  lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);

  ratioH->Draw("e1");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);

  canv1->Print(histName+".png");
  canv1->Print(histName+".pdf","Portrait pdf");


}
