#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH1F.h>

Double_t Background(Double_t * x, Double_t * par) {
  
  //  Double_t result = par[0] + par[1]*x[0];
  Double_t result = par[0]*exp(-par[1]*x[0]);
  return result;

}

Double_t FSR(Double_t * x, Double_t * par) {

  Double_t bC = (x[0]-par[1])/par[2];
  Double_t gauss = par[0]*TMath::Exp(-0.5*bC*bC);
  return gauss;

}


Double_t BackgroundFSR(Double_t * x, Double_t * par) {
  
  Double_t fsr = FSR(x,par);
  Double_t result = par[3]*exp(-par[4]*x[0]);
  return fsr+result;

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

Double_t SignalFSR(Double_t * x, Double_t * par) {

  Double_t signal = SignalX(x,par);
  Double_t parFSR[3];
  parFSR[0] = par[6];
  parFSR[1] = par[7];
  parFSR[2] = par[8];
  Double_t fsr = FSR(x,parFSR);
  return signal + fsr;


}

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

  Double_t signal = SignalX(x,par);
  Double_t bkgd   = BackgroundFSR(x,parBkg);

  Double_t result = signal + bkgd;

  return result;

}



void FitPassAndFail(TString SampleName,
		    TString Name,
		    TString xTitle,
		    TH1F * histPassOld,
		    TH1F * histFailOld,
		    bool fitPass,
		    bool fitFail,
		    bool fitWithFSR,
		    bool reduceBinsPass,
		    bool reduceBinsFail,
		    TCanvas * c1,
		    TCanvas * c2,
		    float * output,
        TString dir_name) {

  gErrorIgnoreLevel = kFatal;

  const Float_t minMass = 75;
  const Float_t maxMass = 105;
  const Float_t xminF = 55;
  const Float_t xmaxF = 120;
  const Float_t xminInt = 75;
  const Float_t xmaxInt = 105;
  
  SetStyle();

  int nBins = histPassOld->GetNbinsX();
  float xmin = histPassOld->GetBinLowEdge(1);
  float xmax = histPassOld->GetBinLowEdge(nBins+1);

  int nBinsNew = nBins;
  //  std::cout << "Enter new number of bins : ";
  //  std::cin >> nBinsNew;
  
  int nBinsNew1 = nBinsNew;
  int nBinsNew2 = nBinsNew;
  if (reduceBinsPass)
    nBinsNew1 = nBinsNew / 2;
  if (reduceBinsFail)
    nBinsNew2 = nBinsNew / 2;

  float *bins1 = new float[nBinsNew1+1];
  float *bins2 = new float[nBinsNew2+1];

  float width1 = (xmax-xmin)/float(nBinsNew1);
  for (int iB=0; iB<=nBinsNew1; ++iB){
    bins1[iB] = xmin + width1*float(iB);
  }

  float width2 = (xmax-xmin)/float(nBinsNew2);
  for (int iB=0; iB<=nBinsNew2; ++iB)
    bins2[iB] = xmin + width2*float(iB);

  TH1F * histPass = new TH1F("histPass","",nBinsNew1,bins1);
  TH1F * histFail = new TH1F("histFail","",nBinsNew2,bins2);


  for (int iB=0;iB<nBins;++iB) {

    float xB = 0.5*(histPassOld->GetBinLowEdge(iB+1)+histPassOld->GetBinLowEdge(iB+2));
    float xC = histPassOld->GetBinContent(iB+1);
    float xE = histPassOld->GetBinError(iB+1);
    int binX = histPass->FindBin(xB);
    float yC = histPass->GetBinContent(binX);
    float yE = histPass->GetBinError(binX);
    float content = xC + yC;
    float error = TMath::Sqrt(xE*xE + yE*yE);
    histPass->SetBinContent(binX,content);
    if (true)
      histPass->SetBinError(binX,error);
    else 
    histPass->SetBinError(binX,0);
    
  }

  for (int iB=0;iB<nBins;++iB) {

    float xB = 0.5*(histFailOld->GetBinLowEdge(iB+1)+histFailOld->GetBinLowEdge(iB+2));
    float xC = histFailOld->GetBinContent(iB+1);
    float xE = histFailOld->GetBinError(iB+1);
    int binX = histFail->FindBin(xB);
    float yC = histFail->GetBinContent(binX);
    float yE = histFail->GetBinError(binX);
    float content = xC + yC;
    float error = TMath::Sqrt(xE*xE + yE*yE);
    histFail->SetBinContent(binX,content);
    if (true)
      histFail->SetBinError(binX,error);
    else 
    histFail->SetBinError(binX,0);
    
  }

  histPass->SetName("Pass");
  histFail->SetName("Fail");


  histPass->GetXaxis()->SetRangeUser(xminF+0.01,xmaxF-0.01);
  histFail->GetXaxis()->SetRangeUser(xminF+0.01,xmaxF-0.01);
  histPass->SetMarkerSize(1.2);
  histFail->SetMarkerSize(1.2);

  TF1 * fitFuncPass = new TF1("fitFuncPass",SignalPlusBackground,xminF,xmaxF,10);
  TF1 * fitFuncFail;
  if (fitWithFSR)
    fitFuncFail = new TF1("fitFuncFail",SignalFSRPlusBackground,xminF,xmaxF,11);
  else
    fitFuncFail = new TF1("fitFuncFail",SignalPlusBackground,xminF,xmaxF,10);

  float xMax = histPass->GetMaximum();
  float rms = histPass->GetRMS();
  float fract = 0.5;
  float massSF = 0.9;
  float x0 = 91.2;
  float Width1 = 3;
  float Width2 = 5;
  //  std::cout << "Alpha = " << Alpha << std::endl;
  //  std::cout << "Beta  = " << Beta << std::endl;

  float y1 = histPass->GetBinContent(histPass->FindBin(xminF+0.01));
  float y2 = histPass->GetBinContent(histPass->FindBin(xmaxF-0.01));
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  float x1 = histPass->GetBinCenter(histPass->FindBin(xminF+0.01));
  float x2 = histPass->GetBinCenter(histPass->FindBin(xmaxF-0.01));
  float bPar = TMath::Log(y1/y2)/(x2-x1);
  float aPar = y1/TMath::Exp(-bPar*x1);

  float xBkg = aPar*exp(-bPar*x0);
  float xTot = histPass->GetBinContent(histPass->FindBin(x0));
  float xSig = xTot - xBkg;


  fitFuncPass->SetParameter(0,xSig);
  fitFuncPass->SetParLimits(0,-1,1e+10);
  fitFuncPass->SetParameter(1,x0);
  fitFuncPass->SetParameter(2,Width1);
  fitFuncPass->SetParLimits(2,0.2*Width1,5*Width1);
  fitFuncPass->SetParameter(3,Width2);
  fitFuncPass->SetParLimits(3,0.2*Width2,5*Width2);
  fitFuncPass->SetParameter(4,Width1);
  fitFuncPass->SetParLimits(4,0.2*Width1,5*Width1);
  fitFuncPass->SetParameter(5,Width2);
  fitFuncPass->SetParLimits(5,0.2*Width2,5*Width2);
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

  //float xBkg = aPar*exp(-bPar*91.2);
  //float xTot = histPass->GetBinContent(histPass->FindBin(91.2));
  //float xSig = xTot - xBkg;
  fitFuncPass->SetParameter(0,xSig);


  c1->cd();
  if (histPass->GetSum()<10)
    fitPass = false;
  histPass->GetXaxis()->SetTitle(xTitle);
  //  histPass->GetYaxis()->SetTitle("Events");
  histPass->GetXaxis()->SetNdivisions(210);
  histPass->GetXaxis()->SetTitleOffset(1.1);
  histPass->GetYaxis()->SetRangeUser(0,1.1*histPass->GetMaximum());
  if (fitPass) histPass->Fit("fitFuncPass","LRQ");
  histPass->Draw("e1");
  TF1 * bkgFuncPass = new TF1("bkgFuncPass",Background,xminF,xmaxF,2);
  bkgFuncPass->SetParameter(0,fitFuncPass->GetParameter(8));
  bkgFuncPass->SetParameter(1,fitFuncPass->GetParameter(9));
  bkgFuncPass->SetLineColor(4);
  bkgFuncPass->SetLineWidth(2);
  bkgFuncPass->SetLineStyle(2);
  if (fitPass) bkgFuncPass->Draw("lsame");
  c1->Update();
  c1->Print(dir_name + "/" + SampleName+"_"+Name+"_pass.png");
  TF1 * sigFuncPass = new TF1("sigFuncPass",Signal,xminF,xmaxF,8);
  for (int iP=0; iP<8; ++iP)
    sigFuncPass->SetParameter(iP,fitFuncPass->GetParameter(iP));

  float nPass = 0.1;
  if (fitPass)
    nPass = sigFuncPass->Integral(xminInt,xmaxInt)/width1;
  else 
    nPass = histPass->Integral(histPass->FindBin(minMass),histPass->FindBin(maxMass)) 
                               -histPass->Integral(histPass->FindBin(60),histPass->FindBin(75))
                               -histPass->Integral(histPass->FindBin(105),histPass->FindBin(120));




  // ***************************************
  y1 = histFail->GetBinContent(histFail->FindBin(xminF+0.01));
  y2 = histFail->GetBinContent(histFail->FindBin(xmaxF-0.01));
  if (y1<0.1) y1 = 0.1;
  if (y2<0.1) y2 = 0.1;
  x1 = histFail->GetBinCenter(histFail->FindBin(xminF+0.01));
  x2 = histFail->GetBinCenter(histFail->FindBin(xmaxF-0.01));
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
  if (fitWithFSR) {
    fitFuncFail->SetParameter(2,Width1);
    fitFuncFail->SetParameter(3,Width1);
    fitFuncFail->SetParameter(4,Width2);
    fitFuncFail->SetParLimits(2,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(3,0.5*Width1,3*Width1);
    fitFuncFail->SetParLimits(4,0.5*Width2,3*Width2);
    fitFuncFail->SetParameter(5,fract);
    fitFuncFail->SetParLimits(5,0.05,0.95);
    fitFuncFail->SetParameter(6,xShift);
    fitFuncFail->SetParameter(7,massSF*x0);
    fitFuncFail->SetParameter(8,sigmaShift);
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
  if (fitFail) histFail->Fit("fitFuncFail","LRQ");
  else  histFail->Draw("e1");
  TF1 * bkgFuncFail = NULL;
  if (fitWithFSR) {
    bkgFuncFail = new TF1("bkgFuncFail",Background,xminF,xmaxF,2);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(9));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(10));
/*    bkgFuncFail = new TF1("bkgFuncFail",BackgroundFSR,xminF,xmaxF,5);
    bkgFuncFail->SetParameter(0,fitFuncFail->GetParameter(6));
    bkgFuncFail->SetParameter(1,fitFuncFail->GetParameter(7));
    bkgFuncFail->SetParameter(2,fitFuncFail->GetParameter(8));
    bkgFuncFail->SetParameter(3,fitFuncFail->GetParameter(9));
    bkgFuncFail->SetParameter(4,fitFuncFail->GetParameter(10));*/
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
//  c2->Print(dir_name + "/" + SampleName+"_"+Name+"_fail.png");
  TF1 * sigFuncFail;
  if (fitWithFSR) {
    sigFuncFail = new TF1("sigFuncFail",SignalFSR,xminF,xmaxF,9);
    for (int iP=0; iP<9; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }  
  else {
    sigFuncFail = new TF1("sigFuncFail",Signal,xminF,xmaxF,8);
    for (int iP=0; iP<8; ++iP)
      sigFuncFail->SetParameter(iP,fitFuncFail->GetParameter(iP));
  }

  float nFail = 0.1;
  if (fitFail)
    nFail = sigFuncFail->Integral(xminInt,xmaxInt)/width2;
  else
    nFail = histFail->Integral(histFail->FindBin(minMass),histFail->FindBin(maxMass))
                    -histFail->Integral(histFail->FindBin(60),histFail->FindBin(75))
                    -histFail->Integral(histFail->FindBin(105),histFail->FindBin(120));

//sanity check for low stat bins
  if (nPass<=0 ) nPass = 0.;
  if (nFail<=0 ) nFail = 0.;                  
//  if (nPass<0.1 && nFail<0.1) nFail = 0.1; 
//  sigFuncFail->SetLineWidth(5);
//  sigFuncFail->SetLineColor(kGreen);
//  sigFuncFail->Draw("lsame");
  c2->Update();
  c2->Print(dir_name + "/" + SampleName+"_"+Name+"_fail.png");

  cout<<endl;
  if(fitWithFSR) cout<<"FSR";
  std::cout << "Eff = " << int(nPass) << "/(" << int(nPass)<<" + "<<int(nFail) << ") = " << nPass/(nPass+nFail) << std::endl;

  TH1F * numH = new TH1F("numH","",1,-0.5,0.5);
  TH1F * denH = new TH1F("denH","",1,-0.5,0.5);

  float EFF = 1;
  numH->SetBinContent(1,nPass);
  denH->SetBinContent(1,nFail+nPass);
  EFF = nPass / (nFail+nPass);
  output[0] = nPass;
  output[1] = nFail;
  EFF *= 100;

  TGraphAsymmErrors * __eff = new TGraphAsymmErrors();
  __eff->Divide(numH,denH);
  
  float errUp   = float(100*__eff->GetErrorYhigh(0));
  float errDown = float(100*__eff->GetErrorYlow(0));
  printf("Eff = %5.2f + %4.2f - %4.2f\n",EFF,errUp,errDown);

}
