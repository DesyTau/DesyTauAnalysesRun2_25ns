// Exponential background
Double_t Background(Double_t * x, Double_t * par) {
  
  Double_t result = par[0]*exp(-par[1]*x[0]);
  return result;

}

// Exponential background with FSR component
// for the Z->ll+gamma events
// FSR : additional gaussian with peak shifted
// towards lower values of m(ll)
Double_t BackgroundFSR(Double_t * x, Double_t * par) {
  
  //  Double_t result = par[0] + par[1]*x[0];
  Double_t result = par[3]*exp(-par[4]*x[0]);
  Double_t bC = (x[0]-par[1])/par[2];
  Double_t gaussFSR = par[0]*TMath::Exp(-0.5*bC*bC);
  return result+gaussFSR;

}

// Signal : double asymmetric gaussian
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

// asymmetric gaussian
// Left side  : single gaussian
// Right side : double gaussian
Double_t SignalX(Double_t * x, Double_t * par) {

  Double_t leftCore  = (x[0]-par[1])/par[2];
  Double_t aRight = (x[0]-par[1])/par[3];
  Double_t bRight = (x[0]-par[1])/par[4];

  Double_t result = 1.0;
  if (x[0]<par[1]) {
    result = TMath::Exp(-0.5*leftCore*leftCore);
  }
  else
    result = par[5]*TMath::Exp(-0.5*aRight*aRight)+(1-par[5])*TMath::Exp(-0.5*bRight*bRight);

  return par[0]*result;

}

// asymmetric gaussian with FSR component
// Left side  : single gaussian
// Right side : double gaussian
// + 1 gaussian with shifted peak value to account for Z->ll+gamma events
Double_t SignalFSR(Double_t * x, Double_t * par) {

  Double_t leftCore  = (x[0]-par[1])/par[2];
  Double_t aRight = (x[0]-par[1])/par[3];
  Double_t bRight = (x[0]-par[1])/par[4];

  Double_t result = 1.0;  
  Double_t bC = (x[0]-par[7])/par[8];
  Double_t gaussFSR = par[6]*TMath::Exp(-0.5*bC*bC);

  if (x[0]<par[1]) {
    result = TMath::Exp(-0.5*leftCore*leftCore);
  }
  else 
    result = par[5]*TMath::Exp(-0.5*aRight*aRight)+(1-par[5])*TMath::Exp(-0.5*bRight*bRight);
  
  return par[0]*result + gaussFSR;

}

// Signal + background 
// (double asymmetric gaussian + exponent)
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

// Signal + exponential background with FSR component
Double_t SignalPlusBackgroundFSR(Double_t * x,
				 Double_t * par) {

  Double_t parBkg[5];
  parBkg[0] = par[6];
  parBkg[1] = par[7];
  parBkg[2] = par[8];
  parBkg[3] = par[9];
  parBkg[4] = par[10];

  Double_t signal = SignalX(x,par);
  Double_t bkgd   = BackgroundFSR(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}

// Signal with FSR component + exponential background 
Double_t SignalFSRPlusBackground(Double_t * x,
				 Double_t * par) {

  Double_t parBkg[2];
  parBkg[0] = par[9];
  parBkg[1] = par[10];

  Double_t signal = SignalFSR(x,par);
  Double_t bkgd   = Background(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}



#include "HttStylesNew.cc"
#include "HtoH.h"

// ranges of m(ll) distribution
// to compute yields in samples
// of passing and failing probes
const Float_t minMass = 85; // lower range to count events (if fit not applied) 
const Float_t maxMass = 97; // upper range to count events (if fit not applied)
const Float_t xminF = 60; // lower range of fitting function
const Float_t xmaxF = 120; // upper range of fitting function
const Float_t xminInt = 85; // lower range to integrate signal function
const Float_t xmaxInt = 97; // upper range to integrate signal function
// FitTagAndProbe : subroutine performing fits of dilepton mass
//                  distributions in the sample of passing and
//                  failing probes
void FitPassAndFail(TString SampleName, // Data or MC
		    TString Name, // generic name (Eta and Pt bin) 
		                  // used for graphic files (png)
		    TString xTitle, // X-axis title of histograms
		    TH1F * histPassOld, // m(ll) histogram with passing probes
		    TH1F * histFailOld, // m(ll) histogram with failing probes
		    bool fitPass, // fit histogram of passing probes
		    bool fitFail, // fit histogram of failing probes
		    bool fitFailWithFSR, // fit histogram of failing
		                         // probes with FSR component
		                         // in the signal function 
		                         // (Gaussian accounting for Z->ll+gamma events)       
		    bool reduceBinsPass, // reduce number of bins by 2 in histogram of passing probes
		    bool reduceBinsFail, // reduce number of bins by 2 in histogram of failing probes
		    TCanvas * c1, // canvas with distribution of passing probes
		    TCanvas * c2, // canvas with distribution of failing probes
		    float * output // output[0] - yield of passing probes
		                   // output[1] - yield of failing probes
		                   // the efficiency is computed as output[0]/(output[1]+output[0]) 
		    ) {
  
  SetStyle();

  int nBins = histPassOld->GetNbinsX();
  float xmin = histPassOld->GetBinLowEdge(1);
  float xmax = histPassOld->GetBinLowEdge(nBins+1);
  float bins1[200];
  float bins2[200];
  int nBinsNew = nBins;
  
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

  // rebinning of histograms 
  TH1F * histPass = (TH1F*)TH1toTH1(histPassOld, nBinsNew1, bins1, true, "_new");
  TH1F * histFail = (TH1F*)TH1toTH1(histFailOld, nBinsNew2, bins2, true, "_new");

  histPass->GetXaxis()->SetRangeUser(xminF+0.01,xmaxF-0.01);
  histFail->GetXaxis()->SetRangeUser(xminF+0.01,xmaxF-0.01);
  // fitting function for passing probes
  TF1 * fitFuncPass = new TF1("fitFuncPass",SignalPlusBackground,xminF,xmaxF,10);

  // fitting function for failing probes
  TF1 * fitFuncFail;
  if (fitFailWithFSR) // signal function with FSR component for failing probes 
    fitFuncFail = new TF1("fitFuncFail",SignalFSRPlusBackground,xminF,xmaxF,11);
  else // standard dou
    fitFuncFail = new TF1("fitFuncFail",SignalPlusBackground,xminF,xmaxF,10);

  //****************************************
  //****** Passing probes ******************
  //**************************************** 

  float xMax = histPass->GetMaximum();
  float rms = histPass->GetRMS();

  // Starting values of fitted parameters
  float x0 = 91.2;  // mZ (peak of distribution)
  float Width1 = 3; // width of narrower gaussian 
  float Width2 = 5; // width of wider gaussian
  float fract = 0.5;  // fractions of Gaussians
                      // with narrower width

  // computing initial parameters for
  // exponential function (background)
  // y = aPar*exp(-bPar*x)
  // from histogram endpoints
  //
  float y1 = histPass->GetBinContent(histPass->FindBin(xminF+0.01));
  float y2 = histPass->GetBinContent(histPass->FindBin(xmaxF-0.01));
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  float x1 = histPass->GetBinCenter(histPass->FindBin(xminF+0.01));
  float x2 = histPass->GetBinCenter(histPass->FindBin(xmaxF-0.01));
  float bPar = TMath::Log(y1/y2)/(x2-x1);
  float aPar = y1/TMath::Exp(-bPar*x1);

  // computing starting value of parameter
  // defining height of double asymmetric
  // Gaussian for signal (xSig)
  float xBkg = aPar*exp(-bPar*x0);
  float xTot = histPass->GetBinContent(histPass->FindBin(x0));
  float xSig = xTot - xBkg;

  fitFuncPass->SetParameter(0,xSig);
  fitFuncPass->SetParLimits(0,0.,1e+10);
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

  c1->cd();
  if (histPass->GetSum()<10) // too few entries to fit histogram with passing probes
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

  // signal component (passing probes)
  // used to compute yield
  TF1 * sigFuncPass = new TF1("sigFuncPass",Signal,xminF,xmaxF,8);
  for (int iP=0; iP<8; ++iP)
    sigFuncPass->SetParameter(iP,fitFuncPass->GetParameter(iP));

  float nPass = 0.1;
  if (fitPass)
    nPass = sigFuncPass->Integral(xminInt,xmaxInt)/width1;
  else 
    nPass = histPass->Integral(histPass->FindBin(minMass),histPass->FindBin(maxMass));

  // *******************************************
  // *********** Failing probes ****************
  // *******************************************

  histFail->SetName("Fail");
  // setting starting points of 
  // the background exponent parameters 
  y1 = histFail->GetBinContent(histFail->FindBin(xminF+0.01));
  y2 = histFail->GetBinContent(histFail->FindBin(xmaxF-0.01));
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  x1 = histFail->GetBinCenter(histFail->FindBin(xminF+0.01));
  x2 = histFail->GetBinCenter(histFail->FindBin(xmaxF-0.01));
  bPar = TMath::Log(y1/y2)/(x2-x1);
  aPar = y1/TMath::Exp(-bPar*x1);

  // height of double gaussian
  xBkg = aPar*exp(-bPar*x0);
  xTot = histFail->GetBinContent(histFail->FindBin(x0));
  xSig = xTot - xBkg;

  // parameters of the FSR component 
  // in signal (simple gaussian)
  float sigmaFSR = 5; // width of FSR component
  float massFSR = 0.9; // central value of FSR component (mZ*massFSR)
  // height of FSR component
  float heightFSR = histFail->GetBinContent(histFail->FindBin(massFSR*x0)) - aPar*exp(-bPar*massFSR*x0);
  if (heightFSR<0) heightFSR = 1e-2;

  fitFuncFail->SetParameter(0,xSig);
  fitFuncFail->SetParLimits(0,-4e+10,4e+10);
  fitFuncFail->SetParameter(1,x0);
  if (fitFailWithFSR) {
    fitFuncFail->SetParameter(2,Width1);
    fitFuncFail->SetParameter(3,Width1);
    fitFuncFail->SetParameter(4,Width2);
    fitFuncFail->SetParLimits(2,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(3,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(4,0.5*Width2,3*Width2);
    fitFuncFail->SetParameter(5,fract);
    fitFuncFail->SetParLimits(5,0.05,0.95);
    fitFuncFail->SetParameter(6,heightFSR);
    fitFuncFail->SetParLimits(6,0,xSig);
    fitFuncFail->SetParameter(7,massFSR*x0);
    fitFuncFail->SetParameter(8,sigmaFSR);
    fitFuncFail->SetParLimits(8,0.2*sigmaFSR,1.5*sigmaFSR);
    fitFuncFail->SetParameter(9,aPar);
    fitFuncFail->SetParameter(10,bPar);
    fitFuncFail->SetParLimits(9,0.5*aPar,2*aPar);
    if (bPar>0)
      fitFuncFail->SetParLimits(10,0.8*bPar,1.2*bPar);
    else
      fitFuncFail->SetParLimits(10,1.2*bPar,0.8*bPar); 
  }
  else {
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
  if (fitFailWithFSR) {
  //    bkgFuncFail = new TF1("bkgFuncFail",BackgroundFSR,xminF,xmaxF,5);
  //    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(6));
  //    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(7));
  //    bkgFuncFail->SetParameter(2,fitFuncFail->GetParameter(8));
  //    bkgFuncFail->SetParameter(3,fitFuncFail->GetParameter(9));
  //    bkgFuncFail->SetParameter(4,fitFuncFail->GetParameter(10));
    bkgFuncFail = new TF1("bkgFuncFail",Background,xminF,xmaxF,2);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(9));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(10));
  }
  else {
    bkgFuncFail = new TF1("bkgFuncFail",Background,xminF,xmaxF,2);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(8));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(9));
  }
  bkgFuncFail->SetLineColor(4);
  bkgFuncFail->SetLineWidth(2);
  bkgFuncFail->SetLineStyle(2);
  if (fitFail)
    bkgFuncFail->Draw("lsame");
  c2->Update();
  c2->Print(SampleName+"_"+Name+"_fail.png");

  // signal function to compute
  // yield in the sample of failing probes
  TF1 * sigFuncFail;
  if (fitFailWithFSR) {
    sigFuncFail = new TF1("sigFuncFail",SignalFSR,xminF,xmaxF,9);
    for (int iP=0; iP<9; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }
  else {
    sigFuncFail = new TF1("sigFuncFail",Signal,xminF,xmaxF,8);
    for (int iP=0; iP<8; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }

  float nFail = 0.0;
  if (fitFail)
    nFail = sigFuncFail->Integral(xminInt,xmaxInt)/width2;
  else
    nFail = histFail->Integral(histFail->FindBin(minMass),histFail->FindBin(maxMass));

  std::cout << "Eff = " << int(nPass) << "/" << int(nPass+nFail) << " = " << nPass/(nPass+nFail) << std::endl;

  // compute efficiency as the ratio
  // of Npass/(Npass+Nfail)
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

// FitTP performs fitting of tag-and-probe histograms
// sequentially for different lepton pt pins and given 
// eta bin (defined by TString histBaseName)
// ****************************************************
// Tag-and-probe histogram naming convention 
// ZMass${etaBin}${ptBin}Pass - m(ll) distribution
//                              in sample of passing probes 
//                              in given (pt,eta) bin
//                              of probed lepton
// ZMass${etaBin}${ptBin}Fail - m(ll) distribution 
//                              in sample of failing probes
//                              in given (pt,eta) bin
//                              of probed lepton 
// **************************************************
// Files :
// SingleMuon_Run201D.root (used for muon Id/Iso)
// SingleElectron_Run2015D (used for electron Id/Iso)
// 
// Names of eta bins :
// ${etaBin} = EtaLt0p9, Eta0p9to1p2, EtaGt1p2 (muons)
// ${etaBin} = Barrel, Endcap (electrons) 
// Names pt bins (both muons and electrons) :
// ${ptBin} = Pt10to15, Pt15to20, Pt20to25,
//            Pt25to30, Pt30to40, Pt40to60
// Names of pt bins are specified in the macro


void FitTP(
	   TString fileName = "SingleMuon_Run2015D", // RooT file with tag-&-probe histograms (w/o *.root extension)
 	   TString histBaseName = "ZMassEta0p9to1p2", // Basename of the histograms to fit (eta bin) 
	   TString xTitle = "muon p_{T} [GeV]", // x axis title (efficiency plot)
	   TString yTitle = "efficiency", // y axis title (efficiency plot)
	   TString xtit = "m_{#mu#mu} [GeV]", // title for dilepton mass distributions
	   float norm = 1 // luminosity normalization factor (1 for data) 
	   ) {
  
  int nPtBins = 6;
  float ptBins[7]; 
  ptBins[0] = 10;
  ptBins[1] = 15;
  ptBins[2] = 20;
  ptBins[3] = 25;
  ptBins[4] = 30;
  ptBins[5] = 40;
  ptBins[6] = 60;

  // define if in the fit of failing probes
  // the FSR component will be used in the 
  // signal function
  bool fitWithFSR[6];
  for (int i=0; i<6; ++i)
    fitWithFSR[i] = false;

  TString PtBins[6];
  PtBins[0] = "Pt10to15";
  PtBins[1] = "Pt15to20";
  PtBins[2] = "Pt20to25";
  PtBins[3] = "Pt25to30";
  PtBins[4] = "Pt30to40";
  PtBins[5] = "Pt40to60";

  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins);

  TFile * file = new TFile(fileName+".root");

  for (int iPt=0; iPt<nPtBins; ++iPt) { // loop ove pt bins
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
    FitPassAndFail(fileName,histBaseName+PtBins[iPt],xtit,histPassOld,histFailOld,fitPass,fitFail,fitWithFSR[iPt],rebinPass,rebinFail,c1,c2,output);
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
  //  eff->GetXaxis()->SetRangeUser(10.01,59.99);
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
