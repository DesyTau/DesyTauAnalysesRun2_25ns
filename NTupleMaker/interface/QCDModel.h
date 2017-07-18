#ifndef QCDModel_h
#define QCDModel_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

using namespace std;

class QCDModel {

 public: 

 QCDModel(TString fileName);
 ~QCDModel();

 double getProbPartonMu(int iMom, int muType, int iflav, int inet); // jet -> mu
 double getProbIsoMu(int iMom, int muType, int ireg, int iflav,int inet); // mu -> iso mu
 double getProbPartonIsoMu(int iMom, int muType, int ireg, int iflav,int inet); // jet -> iso mu
 double getMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass); // pdf(mass) 
 double getPartonMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass); // prob(jet->isoMu)*pdf(mass) 
 double getMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass); // prob(Mu->isoMu)*pdf(mass)
 TH1D * getDRHist(); // dR histo


 private:

 TFile * file;
 // 3 momentum bins x 2 muon types                                                                                                                                    
 TH2D * ProbPartonMu[3][2]; // 2D : [flavor,net-charge]                                                                                                               
 // 3 momentum bins x 2 muon types x 3 regions                                                                                                                        
 TH2D * ProbIsoMu[3][2][3]; // 2D : [flavor,net-charge]                                                                                                               
 // 3 momentum bins x 2 muon types x 3 regions                                                                                                                        
 TH3D * pdfMass[3][2][3]; // 3D : [flavor,net-charge,mass]                                                                                                            
 // deltaR(parton,muon)                                                                                                                                               
 TH1D * deltaRMuParton;


};

#endif
