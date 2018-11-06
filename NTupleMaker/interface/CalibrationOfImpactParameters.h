#ifndef CalibrationOfImpactParameters_h
#define CalibrationOfImpactParameters_h

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>

class CalibrationOfImpactParameters{

public:
   
   CalibrationOfImpactParameters();
   ~CalibrationOfImpactParameters();

   void DoCalibrationForNonPromptLeptons(Float_t pt,Float_t eta,TString IP_Name,Float_t IP, Float_t &IP_cal);
   void DoCalibrationForPromptElectrons(Float_t pt,Float_t eta,TString IP_Name,Float_t IP, Float_t &IP_cal);
   void DoCalibrationForPromptMuons(Float_t pt,Float_t eta,TString IP_Name,Float_t IP, Float_t &IP_cal);
   
private:

   Float_t Calibrate(TH1D * Data, TH1D* MC, Float_t IP);
   TFile *f_np;
   TFile *f_p_ele;
   TFile *f_p_muon;
};

#endif
