/*
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys_WIP.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"
#include "HTT-utilities/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"

#define pi 	3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 		   0.105658
#define tauMass 		   1.77682
#define pionMass 		   0.1396


#include <iostream>
#include <vector>
#include <map>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
//#include "Plotting.h"
//#include "Plotting_Style.h"
//#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"

*/

//int PlotMacro(int argc, char * argv[]){
int PlotMacro(){

  gStyle->SetOptStat(0);
  
  int sample=0;
  int GenReco=1; //1 is RECO
  int Prong=1;//specifies calculation of hadronic vertex acotauta_0*prong*. 0 is impact param, 1 is rho or a particle.
  int DecayMode=1; //Note: genmode1=8. For first reco tau no need to specify anything per def.. 
  int PhiorPsi=0;	  

 TString CutReco="&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&mva17_2>0.5&&mt_1<60&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5||xTrigger>0.5)&&(os>0.5)";
//  TString CutReco="";

  TString CutGen="gen_decaymode_1==8"; 
  
  TString samplename="";   
 // if(sample==0)samplename="ggH_125"; //ggH_125_0_mt_Sync
 // if(sample==0)samplename="ggH_125_0_mt_Sync"; //for checking use line below..
 if(sample==0)samplename="ggh_Update_IP_0_mt_Sync";

 //  if(sample==0)samplename="ggH_125_SingleFile_PubLoc_0_mt_Sync";	
 //  if(sample==0)samplename="ggH_125_SingleFile_0_mt_Sync";	
 //  if(sample==1)samplename="SUSYGluGluHTauTau_120";
 if(sample==1)samplename="SUSYGluGluHTauTau_120_0_mt_Sync";

 //  if(sample==2)samplename="DYJetsToLL"; //old DY sample..
 // if(sample==2)samplename="DYJetsToLL_2019_1_14";
 // if(sample==2)samplename="DYJetsToLL_2019_1_14_SingleFile_0_mt_Sync";

  if(sample==2)samplename="DYJetsToLL_2019_1_14_10File_0_mt_Sync";

  TString GenRecoString="";
  if(GenReco==0) GenRecoString="GEN";
  if(GenReco==1) GenRecoString="RECO";  

  TString RECOObs;
  if(GenReco==0){
   if(PhiorPsi==0) RECOObs="gen_acotautau_0";
   if(PhiorPsi==1) RECOObs="gen_acotautauPsi_0";}

  if(GenReco==1){
   if(PhiorPsi==0) RECOObs="acotautau_0";
   if(PhiorPsi==1) RECOObs="acotautauPsi_0";}

  RECOObs+=Prong;
  RECOObs+=">>CPhist";

  TString RECOCUTString="tau_decay_mode_2==";
  RECOCUTString+=DecayMode;
  RECOCUTString+=CutReco;

  TString CutGenString=CutGen;
  CutGenString+="&&gen_decaymode_2==";
  CutGenString+=DecayMode;
  cout<<"CutGenString "<<CutGenString<<endl;
 
  TString inputfile="/nfs/dust/cms/user/klundert/HiggsCPTauProject/DTSoft_2018_9_18/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/mutau/";  
  inputfile+=samplename;
  inputfile+=".root";

  TString outputfile="/nfs/dust/cms/user/klundert/HiggsCPTauProject/DTSoft_2018_9_18/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/Plots/2019_2_09/";  

cout<<"inputfile "<<inputfile<<endl;

  TFile * f=new TFile(inputfile,"open");
  TTree * tree;

  if(GenReco==0) tree= (TTree*)f->Get("GenTauCheck");
  if(GenReco==1) tree= (TTree*)f->Get("TauCheck");

  
  TCanvas* CPCanvas=new TCanvas("CPCanvas","CPCanvas",700,500);

  int NBins=10;
  TH1D * CPhist;
if(PhiorPsi==0) CPhist= new TH1D("CPhist",samplename,NBins,0,2*TMath::Pi());
if(PhiorPsi==1) CPhist= new TH1D("CPhist",samplename,NBins,-1.2,1.2);

  CPhist->GetXaxis()->SetTitle("#Phi_{CP}");
  CPhist->GetYaxis()->SetTitle("N/Integral");

  if(GenReco==0)  tree->Draw(RECOObs,CutGenString);
  if(GenReco==1)  tree->Draw(RECOObs,RECOCUTString);
//  if(GenReco==1)  tree->Draw(RECOObs);

cout<<"Obsservable "<<RECOObs<<endl;
cout<<"RECOCUTString "<<RECOCUTString<<endl;
cout<<"GEN CUTString "<<CutGenString<<endl;

//  if(GenReco==0)  tree->Draw("acotautau_01>>CPhist","tau_decay_mode_2==1&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&mva17_2>0.5&&mt_1<60&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5||xTrigger>0.5)&&(os>0.5)");


  double integral=CPhist->Integral(1,NBins+1);
  CPhist->Scale(1./integral);
  CPhist->GetYaxis()->SetRangeUser(0,CPhist->GetMaximum()*1.5);
  cout<<"CPhist->GetMaximum() "<<CPhist->GetBinContent(CPhist->GetMaximumBin())<<" CPhist->GetMinimum() "<<CPhist->GetBinContent(CPhist->GetMinimumBin()) <<endl;
  double AmpBaseRatio=( CPhist->GetBinContent(CPhist->GetMaximumBin())-CPhist->GetBinContent(CPhist->GetMinimumBin()) )/(2*CPhist->GetBinContent(CPhist->GetMinimumBin()));
  TString LegendString=GenRecoString;
  LegendString+=" ";
  LegendString+=AmpBaseRatio;

TString LegendString2="Integral ";
double integralhist=CPhist->Integral();
LegendString2+=integral;

  
  cout<<"LegendString 1 "<<LegendString<<endl;
  Int_t dot = LegendString.First('.'); Int_t len = LegendString.Length(); LegendString.Remove(dot+3,len-dot);
  cout<<"LegendString after cut: "<<LegendString<<endl;



  TLegend *Legend=new TLegend(0.7,0.7,0.9,0.9);
  Legend->SetHeader("Ampl/baseline");
  Legend->AddEntry(CPhist,LegendString,"l");
  Legend->AddEntry(CPhist,LegendString2,"l");

  Legend->Draw();
  CPCanvas->Update();



  //save canvas:
  TString pdfstring=outputfile;
  pdfstring+=samplename;
  pdfstring+=GenRecoString;
  pdfstring+="_aco0";
  pdfstring+=Prong;
  pdfstring+="_DecayMode_";
  pdfstring+=DecayMode;
  pdfstring+=".pdf";
  CPCanvas->SaveAs(pdfstring);
  
  
    return 0; 
}


