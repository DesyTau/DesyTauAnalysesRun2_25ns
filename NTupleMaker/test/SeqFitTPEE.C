Double_t PolinBackground(Double_t * x, Double_t * par) {
  
  //  Double_t result = par[0] + par[1]*x[0];
  Double_t result = par[0]*exp(-par[1]*x[0]);
  return result;

}

Double_t DoubleGauss(Double_t * x, Double_t * par) {

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

Double_t SignalPlusPolinBackground(Double_t * x,
                                   Double_t * par) {

  Double_t parBkg[2];
  parBkg[0] = par[8];
  parBkg[1] = par[9];

  Double_t signal = DoubleGauss(x,par);
  Double_t bkgd   = PolinBackground(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}

#include "/nfs/dust/cms/user/bottav/CMSSW_7_4_14/src/DesyTauAnalyses/NTupleMaker/test/HtoH.h"

void FitTP(TString SampleName,
	   TString Name,
	   TH1F * histPassOld,
	   TH1F * histFailOld,
	   bool fitPass,
	   bool fitFail,
	   bool reduceBinsPass,
	   bool reduceBinsFail,
	   TCanvas * c1,
	   TCanvas * c2,
	   float scale,
	   float * output, 
	   bool savePngPlot) {

  int nBins = histPassOld->GetNbinsX();
  float xmin = histPassOld->GetBinLowEdge(1);
  float xmax = histPassOld->GetBinLowEdge(nBins+1);
  float bins1[200];
  float bins2[200];
  int nBinsNew = nBins;
  //  std::cout << "Enter new number of bins : ";
  //  std::cin >> nBinsNew;
  
  int nBinsNew1 = nBinsNew;
  int nBinsNew2 = nBinsNew;
  if (reduceBinsPass)
    nBinsNew1 = nBinsNew / 2;
  if (reduceBinsFail)
    nBinsNew2 = nBinsNew / 2;

  float width1 = (xmax-xmin)/float(nBinsNew1);
  for (int iB=0; iB<=nBinsNew1; ++iB)
    bins1[iB] = xmin + width1*float(iB);

  float width2 = (xmax-xmin)/float(nBinsNew2);
  for (int iB=0; iB<=nBinsNew2; ++iB)
    bins2[iB] = xmin + width2*float(iB);

  float xminF = 60;
  float xmaxF = 120;

  TH1F * histPass = (TH1F*)TH1toTH1(histPassOld, nBinsNew1, bins1, true, "_new");
  TH1F * histFail = (TH1F*)TH1toTH1(histFailOld, nBinsNew2, bins2, true, "_new");
  TF1 * fitFuncPass = new TF1("fitFuncPass",SignalPlusPolinBackground,xminF,xmaxF,10);
  TF1 * fitFuncFail = new TF1("fitFuncFail",SignalPlusPolinBackground,xminF,xmaxF,10);

  float xMax = histPass->GetMaximum();
  float rms = histPass->GetRMS();
  float fract = 0.5;
  float x0 = 91.2;
  float Width1 = 2.4;
  float Width2 = 4;
  //  std::cout << "Alpha = " << Alpha << std::endl;
  //  std::cout << "Beta  = " << Beta << std::endl;

  fitFuncPass->SetParameter(0,xMax);
  fitFuncPass->SetParLimits(0,0,9999999);
  fitFuncPass->SetParameter(1,x0);
  fitFuncPass->SetParameter(2,Width1);
  fitFuncPass->SetParLimits(2,0.5*Width1,2*Width1);
  fitFuncPass->SetParameter(3,Width2);
  fitFuncPass->SetParLimits(3,0.5*Width2,2*Width2);
  fitFuncPass->SetParameter(4,Width1);
  fitFuncPass->SetParLimits(4,0.5*Width1,2*Width1);
  fitFuncPass->SetParameter(5,Width2);
  fitFuncPass->SetParLimits(5,0.5*Width2,2*Width2);
  fitFuncPass->SetParameter(6,fract);
  fitFuncPass->SetParLimits(6,0.05,0.95);
  fitFuncPass->SetParameter(7,fract);
  fitFuncPass->SetParLimits(7,0.05,0.95);
  float y1 = histPass->GetBinContent(1);
  float y2 = histPass->GetBinContent(nBinsNew1);
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  float x1 = histPass->GetBinCenter(1);
  float x2 = histPass->GetBinCenter(nBinsNew1);
  float bPar = TMath::Log(y1/y2)/(x2-x1);
  float aPar = y1/TMath::Exp(-bPar*x1);
  fitFuncPass->SetParameter(8,aPar);
  fitFuncPass->SetParameter(9,bPar);
  fitFuncPass->SetLineWidth(2);

  c1->cd();
  if (histPass->GetSum()<10)
    fitPass = false;
  histPass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  histPass->GetYaxis()->SetTitle("Events / 2 GeV");
  histPass->GetYaxis()->SetTitleOffset(1.1);
  histPass->GetYaxis()->SetRangeUser(0,1.1*histPass->GetMaximum());
  if (fitPass) histPass->Fit("fitFuncPass","LR");
  histPass->SetName("Pass");
  histPass->Draw("e1");
  TF1 * bkgFuncPass = new TF1("bkgFuncPass",PolinBackground,xminF,xmaxF,2);
  bkgFuncPass->SetParameter(0,fitFuncPass->GetParameter(8));
  bkgFuncPass->SetParameter(1,fitFuncPass->GetParameter(9));
  bkgFuncPass->SetLineColor(4);
  bkgFuncPass->SetLineWidth(2);
  bkgFuncPass->SetLineStyle(2);
  if (fitPass) bkgFuncPass->Draw("lsame");
  c1->Update();
  if (savePngPlot) c1->Print(SampleName+"_"+Name+"_pass.png");
  TF1 * sigFuncPass = new TF1("sigFuncPass",DoubleGauss,xminF,xmaxF,8);
  for (int iP=0; iP<8; ++iP)
    sigFuncPass->SetParameter(iP,fitFuncPass->GetParameter(iP));

  float nPass = 0.1;
  if (fitPass)
    nPass = sigFuncPass->Integral(xminF,xmaxF)/width1;
  else 
    nPass = histPass->Integral(10,50);

  std::cout << "nPass = " << nPass << std::endl;

  //
  histFail->SetName("Fail");
  xMax = histFail->GetMaximum();
  fitFuncFail->SetParameter(0,xMax);
  fitFuncFail->SetParLimits(0,0,99999999);
  fitFuncFail->SetParameter(1,fitFuncPass->GetParameter(1));
  fitFuncFail->SetParameter(2,fitFuncPass->GetParameter(2));
  fitFuncFail->SetParLimits(2,0.5*fitFuncPass->GetParameter(2),2*fitFuncPass->GetParameter(2));
  fitFuncFail->SetParameter(3,fitFuncPass->GetParameter(3));
  fitFuncFail->SetParLimits(3,0.5*fitFuncPass->GetParameter(3),2*fitFuncPass->GetParameter(3));
  fitFuncFail->SetParameter(4,fitFuncPass->GetParameter(4));
  fitFuncFail->SetParLimits(4,0.5*fitFuncPass->GetParameter(4),2*fitFuncPass->GetParameter(4));
  fitFuncFail->SetParameter(5,fitFuncPass->GetParameter(5));
  fitFuncFail->SetParLimits(5,0.5*fitFuncPass->GetParameter(5),2*fitFuncPass->GetParameter(5));
  fitFuncFail->SetParameter(6,fitFuncPass->GetParameter(6));
  fitFuncFail->SetParLimits(6,0.05,0.95);
  fitFuncFail->SetParameter(7,fitFuncPass->GetParameter(7));
  fitFuncFail->SetParLimits(7,0.05,0.95);
  y1 = histFail->GetBinContent(1);
  y2 = histFail->GetBinContent(nBinsNew2);
  if (y1==0) y1 = 1;
  if (y2==0) y2 = 1;
  x1 = histFail->GetBinCenter(1);
  x2 = histFail->GetBinCenter(nBinsNew2);
  bPar = TMath::Log(y1/y2)/(x2-x1);
  aPar = y1/TMath::Exp(-bPar*x1);

  fitFuncFail->SetParameter(8,aPar);
  fitFuncFail->SetParameter(9,bPar);
  fitFuncFail->SetLineWidth(2);
  c2->cd();
  if (histFail->GetSum()<10)
    fitFail = false;
  histFail->GetXaxis()->SetTitle("m_{ee} [GeV]");
  histFail->GetYaxis()->SetTitle("Events");
  histFail->GetYaxis()->SetTitleOffset(1.1);
  histFail->GetYaxis()->SetRangeUser(0,1.1*histFail->GetMaximum());
  if (fitFail) histFail->Fit("fitFuncFail","LR");
  else  histFail->Draw("e1");
  TF1 * bkgFuncFail = new TF1("bkgFuncFail",PolinBackground,xminF,xmaxF,2);
  bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(8));
  bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(9));
  bkgFuncFail->SetLineColor(4);
  bkgFuncFail->SetLineWidth(2);
  bkgFuncFail->SetLineStyle(2);
  if (fitFail)
    bkgFuncFail->Draw("lsame");
  c2->Update();
  if (savePngPlot) c2->Print(SampleName+"_"+Name+"_fail.png");

  TF1 * sigFuncFail = new TF1("sigFuncFail",DoubleGauss,xminF,xmaxF,8);
  for (int iP=0; iP<8; ++iP)
    sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  
  float nFail = 0.1;
  if (fitFail)
    nFail = sigFuncFail->Integral(xminF,xmaxF)/width2;
  else
    nFail = histFail->Integral(10,50);

  std::cout << "Eff = " << int(nPass) << "/" << int(nPass+nFail) << " = " << nPass/(nPass+nFail) << std::endl;

  TH1F * numH = new TH1F("numH","",1,-0.5,0.5);
  TH1F * denH = new TH1F("denH","",1,-0.5,0.5);

  float EFF = 1;
  numH->SetBinContent(1,nPass);
  denH->SetBinContent(1,nFail+nPass);
  EFF = nPass / (nFail+nPass);
  output[0] = scale*nPass;
  output[1] = scale*nFail;
  EFF *= 100;

  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
  eff->Divide(numH,denH);
  
  float errUp   = float(100*eff->GetErrorYhigh(0));
  float errDown = float(100*eff->GetErrorYlow(0));
  std::cout << std::endl;
  printf("Eff = %5.2f + %4.2f - %4.2f\n",EFF,errUp,errDown);

}

void SeqFitTPEE(//TString fileName = "SingleElectron_2015D_05Oct",
		TString fileName = "DYJetsToLLM-50_MG",
		TString histBaseName = "ZMassEndcap",
		//TString histBaseName = "ZMassEle17Barrel",
		TString xTitle = "electron p_{T} [GeV]",
		TString yTitle = "Efficiency",
		float scale = 1, 
		bool savePng = false // set to true if want to save all fail pass plots
		) {
  
  TString SampleName("_MC");
  if (fileName.Contains("SingleElectron")) SampleName = "_Data";



  // binning for ID+ISO efficiency
  const int nPtBins = 7; float ptBins[nPtBins+1] = {10, 15, 20, 25, 30, 40, 60, 100};
  TString PtBins[nPtBins+1]  = {"Pt10to15","Pt15to20","Pt20to25","Pt25to30","Pt30to40","Pt40to60", "PtGt60"};
  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins);

  //binning for trigger efficiency
  /*
  const int nPtBinsTrig = 16; float ptBinsTrig[nPtBinsTrig+1]= {10,13,16,19,22,25,28,31,34,37,40,45,50,60,70,100,1000};
  TString PtBinsTrig[nPtBinsTrig+1] = {"Pt10to13","Pt13to16","Pt16to19","Pt19to22","Pt22to25","Pt25to28","Pt28to31","Pt31to34","Pt34to37",
			    "Pt37to40","Pt40to45","Pt45to50","Pt50to60","Pt60to70","Pt70to100","PtGt100"};
  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBinsTrig,ptBinsTrig);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBinsTrig,ptBinsTrig);
  */

  

  TFile * file = new TFile(fileName+".root");
  TH1F * weightsH = (TH1F*)file->Get("histWeightsH");
  float nGen = weightsH->GetSumOfWeights();
  float norm = 6000*41/nGen; 
  if (fileName.Contains("SingleElectron")) norm = 1;

  for (int iPt=0; iPt<nPtBins; ++iPt) {
    TH1F * histPassOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Pass");
    TH1F * histFailOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Fail");
    int nBinsX = histPassOld->GetNbinsX();
    for (int iB=1;iB<=nBinsX;++iB) {
      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
    }
    float output[2];
    TCanvas * c1 = new TCanvas("c1","",700,600);
    TCanvas * c2 = new TCanvas("c2","",700,600);
    bool fitPass = true;
    bool fitFail = true;
    bool rebinPass = false;
    bool rebinFail = true;
    FitTP(fileName,histBaseName+PtBins[iPt],histPassOld,histFailOld,fitPass,fitFail,rebinPass,rebinFail,c1,c2,scale,output,savePng);
    c1->cd();
    c1->Update();
    c2->cd();
    c2->Update();
    int pass = 0;
    numeratorH->SetBinContent(iPt+1,output[0]);
    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);
  }
  // prepare output file
  TFile * outputFile = new TFile(fileName+"_"+histBaseName+".root","recreate");
  outputFile->cd();
  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
  eff->Divide(numeratorH,denominatorH);
  eff->GetXaxis()->SetTitle(xTitle);
  eff->GetXaxis()->SetRangeUser(10.01,59.99);
  eff->GetYaxis()->SetRangeUser(0,1.0);
  eff->GetYaxis()->SetTitle(yTitle);
  eff->SetMarkerSize(1.7);
  eff->SetLineWidth(2);
  TCanvas * canv = new TCanvas("canv","",700,600);
  eff->Draw("APE");
  canv->Update();
  eff->Write(histBaseName+SampleName);

  // copy the eta bins histogram in the output file
  file->cd();
  TH1F * etaBinsH = (TH1F*)file->Get("etaBinsH");
  outputFile->cd(); 
  etaBinsH->Write();
  outputFile->Close();
}
