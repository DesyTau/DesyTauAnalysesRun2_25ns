#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>

TH1F * TH1toTH1(TH1F * histo, int nBinsX, float * xBins, bool setErr,TString suffix) {

  TString name = TString(histo->GetName()) + suffix;

  TH1F * newHisto = new TH1F(name,"",nBinsX,xBins);

  int nBins = histo->GetNbinsX();
  
  for (int iB=0;iB<nBins;++iB) {
    float xB = 0.5*(histo->GetBinLowEdge(iB+1)+histo->GetBinLowEdge(iB+2));
    float xC = histo->GetBinContent(iB+1);
    float xE = histo->GetBinError(iB+1);
    int binX = newHisto->FindBin(xB);
    float yC = newHisto->GetBinContent(binX);
    float yE = newHisto->GetBinError(binX);
    float content = xC + yC;
    float error = TMath::Sqrt(xE*xE + yE*yE);
    newHisto->SetBinContent(binX,content);
    if (setErr)
      newHisto->SetBinError(binX,error); 
    //    else 
    //      newHisto->SetBinError(binX,0);
  }

  return newHisto;

}

TH2F * TH2toTH2(TH2F * histo, int nBinsX, float * xBins, int nBinsY, float * yBins, TString suffix) {

  TString name = TString(histo->GetName()) + "_" + suffix;

  TH2F * newHisto = new TH2F(name,"",nBinsX,xBins,nBinsY,yBins);

  int iBins = histo->GetNbinsX();
  int jBins = histo->GetNbinsY();
  
  for (int iB=0;iB<iBins;++iB) {
    for (int jB=0;jB<jBins;++jB) {
      float xB = 0.5*(histo->GetXaxis()->GetBinLowEdge(iB+1)+histo->GetXaxis()->GetBinLowEdge(iB+2));
      float yB = 0.5*(histo->GetYaxis()->GetBinLowEdge(jB+1)+histo->GetYaxis()->GetBinLowEdge(jB+2));
      float xC = histo->GetBinContent(iB+1,jB+1);
      float xE = histo->GetBinError(iB+1,jB+1);
      int binX = newHisto->GetXaxis()->FindBin(xB);
      int binY = newHisto->GetYaxis()->FindBin(yB);
      float yC = newHisto->GetBinContent(binX,binY);
      float yE = newHisto->GetBinError(binX,binY);
      float content = xC + yC;
      float error = TMath::Sqrt(xE*xE + yE*yE);
      newHisto->SetBinContent(binX,binY,content);
      newHisto->SetBinError(binX,binY,error); 
    }
  }
  return newHisto;

}
