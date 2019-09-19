Double_t DoubleAsymGauss(Double_t * x, Double_t * par) {

  Double_t aLeft  = (x[0]-par[1])/par[2];
  Double_t bLeft  = (x[0]-par[1])/par[3];
  Double_t aRight = (x[0]-par[1])/par[4];
  Double_t bRight = (x[0]-par[1])/par[5];

  Double_t result = 1.0;

  if (x[0]<par[1]) 
    result = par[6]*TMath::Exp(-0.5*aLeft*aLeft)+(1-par[6])*TMath::Exp(-0.5*bLeft*bLeft);
  else 
    result = par[7]*TMath::Exp(-0.5*aRight*aRight)+(1-par[7])*TMath::Exp(-0.5*bRight*bRight);

  
  return par[0]*result;

}

Double_t DoubleSymGauss(Double_t * x, Double_t * par) {

  Double_t a = x[0]/par[1];
  Double_t b = x[0]/par[2];

  Double_t result = 1.0;

  result = par[3]*TMath::Exp(-0.5*a*a)+(1-par[3])*TMath::Exp(-0.5*b*b);
  
  return par[0]*result;

}

// recoilZPerp_${NJets}${Pt}H
// recoilZParal_${NJets}${Pt}H
// $NJets = NJet0, NJet1, NJetGe2
// $Pt =  Pt0to10, Pt10to20, Pt20to30, Pt50to50, PtGt50
void FitRecoil(TH1D * histOld1,
	       TH1D * histOld2,
	       TFile * outputFile,
	       TString histName,
	       float xminF = -120,
	       float xmaxF =  120,
	       TString xtit = "U_{2} [GeV]",
	       TString ytit = "1/N #times dN/dU_{2} [ GeV^{-1} ]",
	       float xminA = -160,
	       float xmaxA =  160,
	       bool asymGauss = false,
	       bool rebin = false,
	       bool logY = true) {
  
  //  setTDRStyle();
  //  gStyle->SetOptStat(0000);
  //  gStyle->SetOptFit(0000);
  
  //  TFile * file1 = new TFile("SingleMuon_2015D.root");
  //  TFile * file2 = new TFile("DYJetsToLL_M-50_ZMuMu_10232015.root");
  //  TH1D * histOld1 = (TH1D*)file1->Get(histName);
  //  TH1D * histOld2 = (TH1D*)file2->Get(histName);
  
  int nBins = histOld1->GetNbinsX();
  float xmin = histOld1->GetBinLowEdge(1);
  float xmax = histOld1->GetBinLowEdge(nBins+1);
  std::cout << histName << "  nBins = " << nBins
	    << "   xmin = " << xmin 
	    << "   xmax = " << xmax << std::endl;
  float bins[200];
  int nBinsNew = nBins;
  //  std::cout << "Enter new number of bins : ";
  //  std::cin >> nBinsNew;
  
  if (rebin)
   nBinsNew = nBins / 2;

  float width = (xmax-xmin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB)
    bins[iB] = xmin + width*float(iB);

  TH1D * hist1 = (TH1D*)TH1DtoTH1D(histOld1, nBinsNew, bins, true, "_new");
  TH1D * hist2 = (TH1D*)TH1DtoTH1D(histOld2, nBinsNew, bins, true, "_new");
  float norm1 = hist1->GetSum();
  float norm2 = hist2->GetSum();
  //  for (int iB=1; iB<=nBinsNew; ++iB) {
  //    hist1->SetBinContent(iB,hist1->GetBinContent(iB)/norm1);
  //    hist1->SetBinError(iB,hist1->GetBinError(iB)/norm1);
  //    hist2->SetBinContent(iB,hist2->GetBinContent(iB)/norm2);
  //    hist2->SetBinError(iB,hist2->GetBinError(iB)/norm2);
  //  }

  TF1 * fitFunc1 = NULL; 
  TF1 * fitFunc2 = NULL; 
  if (asymGauss) {
    fitFunc1 = new TF1("fitFunc1",DoubleAsymGauss,xminF,xmaxF,8);
    fitFunc2 = new TF1("fitFunc2",DoubleAsymGauss,xminF,xmaxF,8);
  }
  else {
    fitFunc1 = new TF1("fitFunc1",DoubleSymGauss,xminF,xmaxF,4);
    fitFunc2 = new TF1("fitFunc2",DoubleSymGauss,xminF,xmaxF,4);
  }

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  float xMax = hist1->GetMaximum();
  float rms = hist1->GetRMS();
  float fract = 0.5;
  float x0 = hist1->GetMean();
  if (asymGauss) {
    fitFunc1->SetParameter(0,xMax);
    fitFunc1->SetParameter(1,x0);
    fitFunc1->SetParameter(2,0.6*rms);
    fitFunc1->SetParameter(3,2*rms);
    fitFunc1->SetParameter(4,0.6*rms);
    fitFunc1->SetParameter(5,2*rms);
    fitFunc1->SetParameter(6,fract);
    fitFunc1->SetParLimits(6,0.05,0.95);
    fitFunc1->SetParameter(7,fract);
    fitFunc1->SetParLimits(7,0.05,0.95);
  }
  else {
    fitFunc1->SetParameter(0,xMax);
    fitFunc1->SetParameter(1,0.6*rms);
    fitFunc1->SetParameter(2,2*rms);
    fitFunc1->SetParameter(3,fract);
    fitFunc1->SetParLimits(3,0.05,0.95);
  }
  fitFunc1->SetLineWidth(2);
  fitFunc1->SetLineColor(2);
  hist1->Fit("fitFunc1","LR");


  xMax = hist2->GetMaximum();
  rms = hist2->GetRMS();
  x0 = hist2->GetMean();
  if (asymGauss) {
    fitFunc2->SetParameter(0,xMax);
    fitFunc2->SetParameter(1,x0);
    fitFunc2->SetParameter(2,0.6*rms);
    fitFunc2->SetParameter(3,2*rms);
    fitFunc2->SetParameter(4,0.6*rms);
    fitFunc2->SetParameter(5,2*rms);
    fitFunc2->SetParameter(6,fract);
    fitFunc2->SetParLimits(6,0.05,0.95);
    fitFunc2->SetParameter(7,fract);
    fitFunc2->SetParLimits(7,0.05,0.95);
  }
  else {
    fitFunc2->SetParameter(0,xMax);
    fitFunc2->SetParameter(1,0.6*rms);
    fitFunc2->SetParameter(2,2*rms);
    fitFunc2->SetParameter(3,fract);
    fitFunc2->SetParLimits(3,0.05,0.95);
  }
  fitFunc2->SetLineWidth(2);
  fitFunc2->SetLineColor(4);
  hist2->Fit("fitFunc2","LR");
  
  hist1->SetLineColor(2);
  hist1->SetMarkerColor(2);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(1.2);
  hist1->GetXaxis()->SetRangeUser(xminA,xmaxA);

  hist2->SetLineColor(4);
  hist2->SetMarkerColor(4);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerSize(1.2);
  
  float yUpper = hist1->GetMaximum();
  if (hist2->GetMaximum()>yUpper) yUpper = hist2->GetMaximum();
  hist1->GetYaxis()->SetRangeUser(0.,1.1*yUpper);
  if (logY) hist1->GetYaxis()->SetRangeUser(1.,2*yUpper);

  delete dummyCanv;

  TCanvas * canv = new TCanvas("canv","",700,700);
   
  hist1->GetXaxis()->SetTitle(xtit);
  hist1->GetXaxis()->SetTitleOffset(1.1);
  hist1->GetYaxis()->SetTitle(ytit);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->Draw("e1");
  hist2->Draw("e1same");

  canv->SetGridx();
  canv->SetGridy();
  canv->SetLogy(logY);
  canv->Update();
  TLegend * leg = new TLegend(0.67,0.72,0.95,0.92);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  TString jetLeg("    N_{jets} = 0 ");
  if (histName.Contains("NJet1"))
    jetLeg = "    N_{jets} = 1 ";
  if (histName.Contains("NJetGe2"))
    jetLeg = "    N_{jets} > 1 ";
  TString ptLeg("[0,10] GeV");
  if (histName.Contains("Pt10to20"))
    ptLeg = "[10,20] GeV";
  if (histName.Contains("Pt20to30"))
    ptLeg = "[20,30] GeV";
  if (histName.Contains("Pt30to50"))
    ptLeg = "[30,50] GeV";
  if (histName.Contains("PtGt50"))
    ptLeg = "> 50 GeV";
  TString header = jetLeg + ptLeg;
  TH1D * dummy = new TH1D("dummy","",1,0,1);
  dummy->SetLineColor(0);
  leg->SetHeader(jetLeg);
  leg->AddEntry(dummy,ptLeg,"l");
  leg->AddEntry(hist1,"Data","lp");
  leg->AddEntry(hist2,"MC","lp");
  leg->Draw();
  canv->Print(histName+"_fit.png");


  std::cout << "Integral (data) = " << fitFunc1->Integral(xminF,xmaxF) << std::endl;
  std::cout << "Integral (MC)   = " << fitFunc2->Integral(xminF,xmaxF) << std::endl;

  fitFunc1->SetParameter(0,fitFunc1->GetParameter(0)/fitFunc1->Integral(xminF,xmaxF));
  fitFunc2->SetParameter(0,fitFunc2->GetParameter(0)/fitFunc2->Integral(xminF,xmaxF));

  std::cout << "Integral (data) = " << fitFunc1->Integral(xminF,xmaxF) << std::endl;
  std::cout << "Integral (MC)   = " << fitFunc2->Integral(xminF,xmaxF) << std::endl;

  outputFile->cd("");
  fitFunc1->Write(histName+"_data");
  fitFunc2->Write(histName+"_mc");
  histOld1->Write(histName+"_hist_data");
  histOld2->Write(histName+"_hist_mc");


}
