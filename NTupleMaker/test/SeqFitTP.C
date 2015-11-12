Double_t Background(Double_t * x, Double_t * par) {
  
  //  Double_t result = par[0] + par[1]*x[0];
  Double_t result = par[0]*exp(-par[1]*x[0]);
  return result;

}

Double_t BackgroundFSR(Double_t * x, Double_t * par) {
  
  //  Double_t result = par[0] + par[1]*x[0];
  Double_t result = par[0]*exp(-par[1]*x[0]);
  Double_t bC = (x[0]-par[3])/par[4];
  Double_t gauss = par[2]*TMath::Exp(-0.5*bC*bC);
  return result+gauss;

}

Double_t Signal(Double_t * x, Double_t * par) {
  
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

Double_t SignalFSR(Double_t * x, Double_t * par) {

  Double_t leftCore  = (x[0]-par[1])/par[2];
  Double_t aRight = (x[0]-par[1])/par[3];
  Double_t bRight = (x[0]-par[1])/par[4];

  Double_t result = 1.0;
  Double_t shift = 0;
  if (x[0]<par[1]) {
    result = TMath::Exp(-0.5*leftCore*leftCore);
  }
  else 
    result = par[5]*TMath::Exp(-0.5*aRight*aRight)+(1-par[5])*TMath::Exp(-0.5*bRight*bRight);
  
  return par[0]*result + shift;


}

Double_t SignalPlusBackground(Double_t * x,
			      Double_t * par) {

  Double_t parBkg[2];
  parBkg[0] = par[8];
  parBkg[1] = par[9];

  Double_t signal = Signal(x,par);
  Double_t bkgd   = Background(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}

Double_t SignalPlusBackgroundFSR(Double_t * x,
				 Double_t * par) {

  Double_t parBkg[5];
  parBkg[0] = par[6];
  parBkg[1] = par[7];
  parBkg[2] = par[8];
  parBkg[3] = par[9];
  parBkg[4] = par[10];

  Double_t signal = SignalFSR(x,par);
  Double_t bkgd   = BackgroundFSR(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}


#include "HttStylesNew.cc"
#include "HtoH.h"

const Float_t minMass = 85;
const Float_t maxMass = 97;
const Float_t xminF = 60;
const Float_t xmaxF = 120;
const Float_t xminInt = 85;
const Float_t xmaxInt = 97;

void FitTP(TString SampleName,
	   TString Name,
	   TString xTitle,
	   TH1F * histPassOld,
	   TH1F * histFailOld,
	   bool fitPass,
	   bool fitFail,
	   bool fitFailCore,
	   bool reduceBinsPass,
	   bool reduceBinsFail,
	   TCanvas * c1,
	   TCanvas * c2,
	   float * output) {

  
  SetStyle();

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


  TH1F * histPass = (TH1F*)TH1toTH1(histPassOld, nBinsNew1, bins1, true, "_new");
  TH1F * histFail = (TH1F*)TH1toTH1(histFailOld, nBinsNew2, bins2, true, "_new");
  TF1 * fitFuncPass = new TF1("fitFuncPass",SignalPlusBackground,xminF,xmaxF,10);
  TF1 * fitFuncFail;
  if (fitFailCore)
    fitFuncFail = new TF1("fitFuncFail",SignalPlusBackground,xminF,xmaxF,10);
  else
    fitFuncFail = new TF1("fitFuncFail",SignalPlusBackgroundFSR,xminF,xmaxF,11);

  float xMax = histPass->GetMaximum();
  float rms = histPass->GetRMS();
  float fract = 0.5;
  float massSF = 0.9;
  float x0 = 91.2;
  float Width1 = 3;
  float Width2 = 5;
  //  std::cout << "Alpha = " << Alpha << std::endl;
  //  std::cout << "Beta  = " << Beta << std::endl;

  float y1 = histPass->GetBinContent(1);
  float y2 = histPass->GetBinContent(nBinsNew1);
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  float x1 = histPass->GetBinCenter(1);
  float x2 = histPass->GetBinCenter(nBinsNew1);
  float bPar = TMath::Log(y1/y2)/(x2-x1);
  float aPar = y1/TMath::Exp(-bPar*x1);

  float xBkg = aPar*exp(-bPar*x0);
  float xTot = histPass->GetBinContent(histPass->FindBin(x0));
  float xSig = xTot - xBkg;

  fitFuncPass->SetParameter(0,xSig);
  fitFuncPass->SetParLimits(0,-1e+10,1e+10);
  fitFuncPass->SetParameter(1,x0);
  fitFuncPass->SetParameter(2,Width1);
  fitFuncPass->SetParLimits(2,0.3*Width1,5*Width1);
  fitFuncPass->SetParameter(3,Width2);
  fitFuncPass->SetParLimits(3,0.3*Width2,5*Width2);
  fitFuncPass->SetParameter(4,Width1);
  fitFuncPass->SetParLimits(4,0.3*Width1,5*Width1);
  fitFuncPass->SetParameter(5,Width2);
  fitFuncPass->SetParLimits(5,0.3*Width2,5*Width2);
  fitFuncPass->SetParameter(6,fract);
  fitFuncPass->SetParLimits(6,0.05,0.95);
  fitFuncPass->SetParameter(7,fract);
  fitFuncPass->SetParLimits(7,0.05,0.95);
  fitFuncPass->SetParameter(8,aPar);
  fitFuncPass->SetParameter(9,bPar);
  fitFuncPass->SetParLimits(8,0.8*aPar,1.2*aPar);
  if (bPar>0)
    fitFuncPass->SetParLimits(9,0.8*bPar,1.2*bPar);
  else
    fitFuncPass->SetParLimits(9,1.2*bPar,0.8*bPar);

  fitFuncPass->SetLineWidth(2);
  fitFuncPass->SetLineColor(2);

  float xBkg = aPar*exp(-bPar*91.2);
  float xTot = histPass->GetBinContent(histPass->FindBin(91.2));
  float xSig = xTot - xBkg;
  fitFuncPass->SetParameter(0,xSig);

  c1->cd();
  if (histPass->GetSum()<10)
    fitPass = false;
  histPass->GetXaxis()->SetTitle(xTitle);
  //  histPass->GetYaxis()->SetTitle("Events");
  histPass->GetXaxis()->SetNdivisions(210);
  histPass->GetXaxis()->SetTitleOffset(1.1);
  histPass->GetYaxis()->SetRangeUser(0,1.1*histPass->GetMaximum());
  if (fitPass) histPass->Fit("fitFuncPass","LR");
  histPass->SetName("Pass");
  histPass->Draw("e1");
  TF1 * bkgFuncPass = new TF1("bkgFuncPass",Background,xminF,xmaxF,2);
  bkgFuncPass->SetParameter(0,fitFuncPass->GetParameter(8));
  bkgFuncPass->SetParameter(1,fitFuncPass->GetParameter(9));
  bkgFuncPass->SetLineColor(4);
  bkgFuncPass->SetLineWidth(2);
  bkgFuncPass->SetLineStyle(2);
  if (fitPass) bkgFuncPass->Draw("lsame");
  c1->Update();
  c1->Print(SampleName+"_"+Name+"_pass.png");
  TF1 * sigFuncPass = new TF1("sigFuncPass",Signal,xminF,xmaxF,8);
  for (int iP=0; iP<8; ++iP)
    sigFuncPass->SetParameter(iP,fitFuncPass->GetParameter(iP));

  float nPass = 0.1;
  if (fitPass)
    nPass = sigFuncPass->Integral(xminInt,xmaxInt)/width1;
  else 
    nPass = histPass->Integral(histPass->FindBin(minMass),histPass->FindBin(maxMass));

  // ***************************************
  histFail->SetName("Fail");
  y1 = histFail->GetBinContent(1);
  y2 = histFail->GetBinContent(nBinsNew2);
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  x1 = histFail->GetBinCenter(1);
  x2 = histFail->GetBinCenter(nBinsNew2);
  bPar = TMath::Log(y1/y2)/(x2-x1);
  aPar = y1/TMath::Exp(-bPar*x1);

  xBkg = aPar*exp(-bPar*x0);
  xTot = histFail->GetBinContent(histFail->FindBin(x0));
  xSig = xTot - xBkg;

  float sigmaShift = 10;
  float xShift = histFail->GetBinContent(histFail->FindBin(massSF*x0)) - aPar*exp(-bPar*massSF*x0);

  fitFuncFail->SetParameter(0,xSig);
  fitFuncFail->SetParLimits(0,-4e+10,4e+10);
  fitFuncFail->SetParameter(1,x0);
  if (fitFailCore) {
    fitFuncFail->SetParameter(2,Width1);
    fitFuncFail->SetParameter(3,Width2);
    fitFuncFail->SetParameter(4,Width1);
    fitFuncFail->SetParameter(5,Width2);
    fitFuncFail->SetParLimits(2,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(3,0.5*Width2,3*Width2);
    fitFuncFail->SetParLimits(4,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(5,0.5*Width2,3*Width2);
    fitFuncFail->SetParameter(6,fract);
    fitFuncFail->SetParLimits(6,0.05,0.95);
    fitFuncFail->SetParameter(7,fract);
    fitFuncFail->SetParLimits(7,0.05,0.95);
    fitFuncFail->SetParameter(8,aPar);
    fitFuncFail->SetParameter(9,bPar);
    fitFuncFail->SetParLimits(8,0.8*aPar,1.2*aPar);
    if (bPar>0)
      fitFuncFail->SetParLimits(9,0.8*bPar,1.2*bPar);
    else
      fitFuncFail->SetParLimits(9,1.2*bPar,0.8*bPar);
  }
  else {
    fitFuncFail->SetParameter(2,Width1);
    fitFuncFail->SetParameter(3,Width1);
    fitFuncFail->SetParameter(4,Width2);
    fitFuncFail->SetParLimits(2,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(3,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(4,0.5*Width2,3*Width2);
    fitFuncFail->SetParameter(5,fract);
    fitFuncFail->SetParLimits(5,0.05,0.95);
    fitFuncFail->SetParameter(6,aPar);
    fitFuncFail->SetParameter(7,bPar);
    fitFuncFail->SetParLimits(6,0.5*aPar,2*aPar);
    if (bPar>0)
      fitFuncFail->SetParLimits(7,0.8*bPar,1.2*bPar);
    else
      fitFuncFail->SetParLimits(7,1.2*bPar,0.8*bPar); 
    fitFuncFail->SetParameter(8,xShift);
    fitFuncFail->SetParameter(9,massSF*x0);
    fitFuncFail->SetParameter(10,sigmaShift);
  }


  fitFuncFail->SetLineWidth(2);
  fitFuncFail->SetLineColor(2);
  c2->cd();
  if (histFail->GetSum()<10)
    fitFail = false;
  histFail->GetXaxis()->SetTitle(xTitle);
  //  histFail->GetYaxis()->SetTitle("Events");
  histPass->GetXaxis()->SetNdivisions(210);
  histFail->GetXaxis()->SetTitleOffset(1.1);
  histFail->GetYaxis()->SetRangeUser(0,1.1*histFail->GetMaximum());
  if (fitFail) histFail->Fit("fitFuncFail","LR");
  else  histFail->Draw("e1");
  TF1 * bkgFuncFail = NULL;
  if (fitFailCore) {
    bkgFuncFail = new TF1("bkgFuncFail",Background,xminF,xmaxF,2);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(8));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(9));
  }
  else {
    bkgFuncFail = new TF1("bkgFuncFail",BackgroundFSR,xminF,xmaxF,5);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(6));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(7));
    bkgFuncFail->SetParameter(2,fitFuncFail->GetParameter(8));
    bkgFuncFail->SetParameter(3,fitFuncFail->GetParameter(9));
    bkgFuncFail->SetParameter(4,fitFuncFail->GetParameter(10));
  }
  bkgFuncFail->SetLineColor(4);
  bkgFuncFail->SetLineWidth(2);
  bkgFuncFail->SetLineStyle(2);
  if (fitFail)
    bkgFuncFail->Draw("lsame");
  c2->Update();
  c2->Print(SampleName+"_"+Name+"_fail.png");
  TF1 * sigFuncFail;
  if (fitFailCore) {
    sigFuncFail = new TF1("sigFuncFail",Signal,xminF,xmaxF,8);
    for (int iP=0; iP<8; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }  
  else {
    sigFuncFail = new TF1("sigFuncFail",SignalFSR,xminF,xmaxF,6);
    for (int iP=0; iP<6; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }

  float nFail = 0.1;
  if (fitFail)
    nFail = sigFuncFail->Integral(xminInt,xmaxInt)/width2;
  else
    nFail = histFail->Integral(histFail->FindBin(minMass),histFail->FindBin(maxMass));

  std::cout << "Eff = " << int(nPass) << "/" << int(nPass+nFail) << " = " << nPass/(nPass+nFail) << std::endl;

  TH1F * numH = new TH1F("numH","",1,-0.5,0.5);
  TH1F * denH = new TH1F("denH","",1,-0.5,0.5);

  float EFF = 1;
  numH->SetBinContent(1,nPass);
  denH->SetBinContent(1,nFail+nPass);
  EFF = nPass / (nFail+nPass);
  output[0] = nPass;
  output[1] = nFail;
  EFF *= 100;

  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
  eff->Divide(numH,denH);
  
  float errUp   = float(100*eff->GetErrorYhigh(0));
  float errDown = float(100*eff->GetErrorYlow(0));
  std::cout << std::endl;
  printf("Eff = %5.2f + %4.2f - %4.2f\n",EFF,errUp,errDown);

}

// DYJetsToEE_M-50_13TeV-amcatnloFXFX-pythia8
// DYJetsToMuMu_M-50_13TeV-amcatnloFXFX-pythia8
// SingleMuon_Run2015D
// SingleElectron_Run2015D
void SeqFitTP(TString fileName = "DYJetsToMuMu_M-50_13TeV-amcatnloFXFX-pythia8",
	      TString histBaseName = "ZMassIsoMuEtaLt0p9",
	      TString xTitle = "muon p_{T} [GeV]",
	      TString yTitle = "efficiency",
	      TString xtit = "m_{#mu#mu} [GeV]") {
  
  TString SampleName("MC");
  if (fileName.Contains("Run2015"))
    SampleName = "Data";

  int nPtBins = 6;
  float ptBins[30]; 
  ptBins[0] = 10;
  ptBins[1] = 15;
  ptBins[2] = 20;
  ptBins[3] = 25;
  ptBins[4] = 30;
  ptBins[5] = 40;
  ptBins[6] = 60;
  ptBins[7] = 1000;

  bool fitFailCore[30];
  for (int i=0; i<30; ++i)
    fitFailCore[i] = true;

  //  fitFailCore[2] = false;
  //  fitFailCore[3] = false;
  //  fitFailCore[4] = false;
  //  fitFailCore[5] = false;


  TString PtBins[30];
  PtBins[0] = "Pt10to15";
  PtBins[1] = "Pt15to20";
  PtBins[2] = "Pt20to25";
  PtBins[3] = "Pt25to30";
  PtBins[4] = "Pt30to40";
  PtBins[5] = "Pt40to60";
  PtBins[6] = "PtGt60";

  if (histBaseName.Contains("ZMassMu")||histBaseName.Contains("ZMassIsoMu")||
      histBaseName.Contains("ZMassEle")||histBaseName.Contains("ZMassIsoEle")) {
    nPtBins = 14;
    ptBins[0] = 10;
    ptBins[1] = 13;
    ptBins[2] = 16;
    ptBins[3] = 19;
    ptBins[4] = 22;
    ptBins[5] = 25;
    ptBins[6] = 28;
    ptBins[7] = 31;
    ptBins[8] = 34;
    ptBins[9] = 37;
    ptBins[10] = 40;
    ptBins[11] = 45;
    ptBins[12] = 50;
    ptBins[13] = 60;
    ptBins[14] = 70;
    ptBins[15] = 100;
    ptBins[16] = 1000;
    PtBins[0] = "Pt10to13";
    PtBins[1] = "Pt13to16";
    PtBins[2] = "Pt16to19";
    PtBins[3] = "Pt19to22";
    PtBins[4] = "Pt22to25";
    PtBins[5] = "Pt25to28";
    PtBins[6] = "Pt28to31";
    PtBins[7] = "Pt31to34";
    PtBins[8] = "Pt34to37";
    PtBins[9] = "Pt37to40";
    PtBins[10] = "Pt40to45";
    PtBins[11] = "Pt45to50";
    PtBins[12] = "Pt50to60";
    PtBins[13] = "Pt60to70";
    PtBins[14] = "Pt70to100";
    PtBins[15] = "PtGt100";

  }

  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins);

  TFile * file = new TFile(fileName+".root");
  TH1F * weightsH = (TH1F*)file->Get("histWeightsH");
  float nGen = weightsH->GetSumOfWeights();
  float norm = 6025*1200/nGen; 
  if (fileName.Contains("SingleMu")||fileName.Contains("SingleEle")) norm = 1;

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
    bool rebinFail = false;
    FitTP(SampleName,histBaseName+PtBins[iPt],xtit,histPassOld,histFailOld,fitPass,fitFail,fitFailCore[iPt],rebinPass,rebinFail,c1,c2,output);
    c1->cd();
    c1->Update();
    c2->cd();
    c2->Update();
    numeratorH->SetBinContent(iPt+1,output[0]);
    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);
  }
  TFile * outputFile = new TFile(fileName+"_"+histBaseName+".root","recreate");
  outputFile->cd();
  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
  eff->Divide(numeratorH,denominatorH);
  eff->GetXaxis()->SetTitle(xTitle);
  eff->GetXaxis()->SetRangeUser(10.01,59.99);
  eff->GetYaxis()->SetRangeUser(0,1.0);
  eff->GetYaxis()->SetTitle(yTitle);
  eff->GetXaxis()->SetTitleOffset(1.1);
  eff->GetYaxis()->SetTitleOffset(1.1);
  eff->SetMarkerSize(1.7);
  eff->SetLineWidth(2);
  TCanvas * canv = new TCanvas("canv","",700,600);
  eff->Draw("APE");
  canv->Update();
  eff->Write(histBaseName);
  outputFile->Close();


}
