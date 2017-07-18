#include "DesyTauAnalyses/NTupleMaker/interface/QCDModel.h"

QCDModel::QCDModel(TString fileName) {

  file = new TFile(fileName);

  deltaRMuParton = (TH1D*)file->Get("deltaRPartonMuH");

  TString partonFlavor[4] = {"uds","g","c","b"};
  TString muonPartonNetCharge[2] = {"opposite","same"};
  TString partonMomRange[3] = {"Lt50","50to100","Gt100"};
  TString muonType[2] = {"HighMu","LowMu"}; 
  TString region[3] = {"Iso","LooseIso","Sb"};

    for (int imom=0; imom<3; ++imom) {
    for (int imu=0; imu<2; ++imu) {
      ProbPartonMu[imom][imu] = (TH2D*)file->Get("ProbPartonMu_"+partonMomRange[imom]+"_"+muonType[imu]);
      for (int ireg=0; ireg<3; ++ireg) {
	ProbIsoMu[imom][imu][ireg] = (TH2D*)file->Get("ProbIsoMu_"+partonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	pdfMass[imom][imu][ireg] = (TH3D*)file->Get("pdfMass_"+partonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
      }
    }
  }



}

QCDModel::~QCDModel() {

}

double QCDModel::getProbPartonMu(int iMom, int muType, int iflav, int inet) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;

  double output = ProbPartonMu[iMom][muType]->GetBinContent(iflav+1,inet+1);

  return output;

}

double QCDModel::getProbIsoMu(int iMom, int muType, int ireg, int iflav, int inet) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  double output = ProbIsoMu[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1);

  return output;

}

double QCDModel::getProbPartonIsoMu(int iMom, int muType, int ireg, int iflav, int inet) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  double output = 
    ProbPartonMu[iMom][muType]->GetBinContent(iflav+1,inet+1) * 
    ProbIsoMu[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1);

  return output;

}

double QCDModel::getMassPdf(int iMom, int muType, int ireg, int iflav, int inet, double mass) {
  
  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  int imass = int(mass);

  double output = pdfMass[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1,imass+1);
  return output;

}

double QCDModel::getPartonMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  int imass = int(mass);
  double output = getProbPartonIsoMu(iMom,muType,ireg,iflav,inet) *
    pdfMass[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1,imass+1);

  return output;

}


double QCDModel::getMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  int imass = int(mass);
  double output = getProbIsoMu(iMom,muType,ireg,iflav,inet) *
    pdfMass[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1,imass+1);

  return output;

}


TH1D * QCDModel::getDRHist() {

  return deltaRMuParton;

}
