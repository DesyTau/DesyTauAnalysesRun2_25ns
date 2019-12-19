#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"
// massSelH
//  TH1D * ptLeadingMuSelH = new TH1D("ptLeadingMuSelH","",100,0,200);
//  TH1D * ptTrailingMuSelH = new TH1D("ptTrailingMuSelH","",100,0,200);
//  TH1D * etaLeadingMuSelH = new TH1D("etaLeadingMuSelH","",50,-2.5,2.5);
//  TH1D * etaTrailingMuSelH = new TH1D("etaTrailingMuSelH","",50,-2.5,2.5);
//  TH1D * massSelH = new TH1D("massSelH","",200,0,200);
//  TH1D * massExtendedSelH = new TH1D("massExtendedSelH","",200,0,200);
//  TH1D * dimuonPtSelH = new TH1D("dimuonPtSelH","",200,0,200);
//  TH1D * metSelH  = new TH1D("metSelH","",200,0,400);

void PlotSamplesRatio(TString histName ="metSelH",
		      TString xtitle = "E_{T}^{mis} [GeV]",
		      TString ytitle = "Events / 10 GeV",
		      float xLower =      0,
		      float xUpper =    400,
		      float yLower =      1,
		      bool logY = true,
		      bool drawLeg = true) {

  float qcdScale = 2;
  float ttScale = 1;
  float DYscaleLow = 1.0;
  float DYscaleHigh = 1;
  TString suffix = "_recoil";

  TString dir("./final/");

  SetStyle();

  TFile * file = new TFile(dir+"SingleMuon_Run2016.root");
  TFile * fileSS = new TFile(dir+"SingleMuon_Run2016_ss.root");  
  
  TString samples[10] = {"WW_13TeV-pythia8",                      // (0)
			 "WZ_13TeV-pythia8",                      // (1) 
			 "ZZ_13TeV-pythia8",                      // (2)
			 "WJetsToLNu_13TeV-madgraphMLM-pythia8",  // (3)
			 "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8",           // (4)
			 "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8",       // (5)
			 "ST_t-channel_top_4f_InclusiveDecays_13TeV-powheg-pythia8",    // (6)
			 "ST_t-channel_antitop_4f_InclusiveDecays_13TeV-powheg-pythia8",// (7)
			 "TT_13TeV-powheg-pythia8",                           // (8)
			 "DYJetsToLL_M-50_13TeV_madgraphMLM_pythia8"+suffix   // (9)
  };

  float xsec[10] = {118.7,       // WW   (0)
		    27.68,       // WZ   (1)
		    12.19,       // ZZ   (2)
		    52760*1.166, // W+Jets (3)
		    35.85,       // tW_top (4)
		    35.85,       // tW_antitop (5) 
		    44.33,       // top-tchannel (6)
		    26.38,       // antitop-tchannel (7)
		    831.8,       // TT      (8)
		    5345*1.079   // DYJets  (9)
  };

  float lumi = 35900;

  TH1D * histDataOld = (TH1D*)file->Get(histName);
  TH1D * histDataOldSS = (TH1D*)fileSS->Get(histName);

  int nBins = histDataOld->GetNbinsX();
  float xMin = histDataOld->GetBinLowEdge(1);
  float xMax = histDataOld->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  float bins[300];
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

  TH1D * histData = TH1DtoTH1D(histDataOld,nBinsNew,bins,true,"_Data_new");
  TH1D * qcdHist  = TH1DtoTH1D(histDataOldSS,nBinsNew,bins,true,"_qcd");

  TH1D * ewkHist = new TH1D("ewkHist","",nBinsNew,bins);
  TH1D * ttHist  = new TH1D("ttHist","",nBinsNew,bins);
  TH1D * zHist = new TH1D("zHist","",nBinsNew,bins);

  int nSamples = 10;

  //  return;
  for (int iS=0; iS<nSamples; ++iS) {
    //    std::cout << "Sample = " << iS << std::endl;
    TFile * fileMC = new TFile(dir+samples[iS]+".root");
    TH1D * histOld = (TH1D*)fileMC->Get(histName);
    TH1D * hist = TH1DtoTH1D(histOld,nBinsNew,bins,true,"_new_"+samples[iS]);
    TH1D * eventCount = (TH1D*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[iS]*lumi/nGen;
    TH1D * tempHist = ewkHist;
    if (iS==8)
      tempHist = ttHist;
    if (iS==9)
      tempHist = zHist;
    std::cout << samples[iS] << " nGen = " << nGen << " " << hist->GetSumOfWeights() << std::endl;
    tempHist->Add(tempHist,hist,1.,norm);

  }

  //  float dataEvents = 0;
  //  float ttEvents = 0;
  std::cout << "QCD : " << qcdScale*qcdHist->GetSumOfWeights() << std::endl;
  std::cout << "EWK : " << ewkHist->GetSumOfWeights() << std::endl;
  std::cout << "TTJ : " << ttHist->GetSumOfWeights() << std::endl;
  std::cout << "ZMM : " << zHist->GetSumOfWeights() << std::endl;

  float lumiSys = 0.02;
  float lepSys = 0.06;
  float ttSys = 0.1;
  float ewkSys = 0.15;
  float qcdSys = 0.2;

  for (int iB=1; iB<=nBinsNew; ++iB) {

    float qcdX = qcdScale*qcdHist->GetBinContent(iB);
    float qcdE = qcdScale*qcdHist->GetBinError(iB);
    
    float ewkX = ewkHist->GetBinContent(iB);
    float ewkE = ewkHist->GetBinError(iB);
    
    float ttX  = ttHist->GetBinContent(iB);
    float ttE  = ttHist->GetBinError(iB);
    
    float zX  = zHist->GetBinContent(iB);
    float zE  = zHist->GetBinError(iB);
    
    if (zX<0) zX = 0;
    
    float ttErr   = ttX*ttSys;
    float ewkErr  = ewkX*ewkSys;
    float qcdErr  = qcdX*qcdSys;
    
    ewkX += qcdX;
    ttX  += ewkX;
    zX   += ttX;

    float lumiErr = zX*lumiSys;
    float lepErr  = zX*lepSys;

    float totErr = TMath::Sqrt(lumiErr*lumiErr+
			       lepErr*lepErr+
			       ttErr*ttErr+
			       ewkErr*ewkErr+
			       qcdErr*qcdErr+
			       ewkE*ewkE+
			       ttE*ttE+
			       qcdE*qcdE+
			       zE*zE);
    if (totErr>zX) totErr=0.8*zX;
    
    qcdHist->SetBinContent(iB,qcdX);
    ewkHist->SetBinContent(iB,ewkX);
    ttHist->SetBinContent(iB,ttX);
    zHist->SetBinContent(iB,zX);
    zHist->SetBinError(iB,totErr);

  }

  std::cout << "BKG : " << zHist->GetSumOfWeights() << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << std::endl;

  TH1D * bkgdErr = (TH1D*)zHist->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3344);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    qcdHist->SetBinError(iB,0);
    ewkHist->SetBinError(iB,0);
    ttHist->SetBinError(iB,0);
    zHist->SetBinError(iB,0);
  }

  InitData(histData);
  InitHist(qcdHist,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(ttHist,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(ewkHist,"","",TColor::GetColor("#6F2D35"),1001);
  InitHist(zHist,"","",TColor::GetColor("#4496C8"),1001);
  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.6);
  //  histData->GetYaxis()->SetTitleSize(0.06);
  histData->GetXaxis()->SetRangeUser(xLower,xUpper);
  float yUpper = histData->GetMaximum();
  if (logY)
    histData->GetYaxis()->SetRangeUser(yLower,10*yUpper);
  else
    histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);

  histData->SetMarkerSize(1.2);
  histData->GetXaxis()->SetLabelSize(0);
  //  histData->GetYaxis()->SetLabelSize(0.06);

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
  ewkHist->Draw("sameh");
  qcdHist->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
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
  char Chi2[100];
  sprintf(Chi2,"%6.1f",chi2/float(nBinsNew));

  TLegend * leg = new TLegend(0.6,0.5,0.9,0.85);
  SetLegendStyle(leg);
  leg->SetTextSize(0.046);
  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(zHist,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(ttHist,"t#bar{t}","f");
  leg->AddEntry(ewkHist,"electroweak","f");
  leg->AddEntry(qcdHist,"QCD multijets","f");
  if (drawLeg) leg->Draw();
  writeExtraText = true;
  extraText   = "Preliminary";
  CMS_lumi(upper,4,33); 
  //  plotchannel("#mu#mu");
  TLatex chi2Latex;
  chi2Latex.SetNDC();
  TString chi2Label = "#chi^{2}/ndof = " + TString(Chi2);
  chi2Latex.SetTextSize(0.045);
  //  chi2Latex.DrawLatex(0.64,0.4,chi2Label);

  if (logY) upper->SetLogy(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  //ratioH->GetYaxis()->SetRangeUser(0.01,1.99);
  ratioH->GetYaxis()->SetRangeUser(0.5,1.5);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleOffset(4.0);
  ratioH->GetXaxis()->SetTickLength(0.07);

  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  //  ratioH->GetYaxis()->SetTitleSize(0.8);
  ratioH->GetYaxis()->SetTitleOffset(1.5);
  ratioH->GetYaxis()->SetTickLength(0.04);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = zHist->GetBinContent(iB);
    ratioErrH->SetBinContent(iB,1.0);
    ratioErrH->SetBinError(iB,0.0);
    float xBkg = bkgdErr->GetBinContent(iB);
    float errBkg = bkgdErr->GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      float err = TMath::Sqrt(relErr*relErr);
      ratioErrH->SetBinError(iB,err);

    }
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
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
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
  ratioErrH->Draw("e2same");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Print(histName+suffix+".png");
  //  canv1->Print(histName+suffix+".pdf","Portrait pdf");

  //  TFile * fileO = new TFile("met.root","recreate");
  //  fileO->cd("");
  //  canv1->Print("figures/"+histName+suffix+".png");
  //  canv1->Print("figures/"+histName+suffix+".pdf","Portrait pdf");
  //  canv1->Write("canv");
  //  fileO->Close();


}
