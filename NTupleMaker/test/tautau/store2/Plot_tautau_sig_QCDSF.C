
#include <iostream>
#include <vector>
#include <map>
#include <iomanip>


#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"

using namespace std;
TH1D* histdata(TString variable,int nBins,float xmin,float xmax,bool isQCD=true){
  TFile* f_data = new TFile("DATA_Tau_QCD_tt_Sync.root");
  TTree* tree_data = (TTree*)f_data->Get("TauCheck");
  TH1D* h_data = new TH1D("h_data","",nBins,xmin,xmax);h_data->Sumw2();
  TString Cuts = "&&pt_1>20&&pt_2>20&&Prompt_pT>50&&extraelec_veto<0.5&&extramuon_veto<0.5";
  if(!isQCD){
    tree_data->Draw(variable+">>h_data","os>0.5"+Cuts);
  }
  else{
    tree_data->Draw(variable+">>h_data","os<0.5"+Cuts);
  }
  return h_data;
}
//"DYJetsToLL_M-50_TuneCP5.root","TTTo2L2Nu_TuneCP5.root",
void Plot_tautau_sig(TString variable,int nBins,float xmin,float xmax,TString Weight = "weight*"){
  double lumi = 41465;double Wnorm=1.0;
  TString zptmassweight="zptweight*"; 
  TString file_names[6] = {"DYJetsToLL_tt_Sync.root","DYJetsToLL_tt_Sync.root","WJetsToLNu_tt_Sync.root","TTTo2L2Nu_tt_Sync.root","DYJetsToLL_tt_Sync.root"};
  TString Cuts = "&&pt_1>20&&pt_2>20&&Prompt_pT>50&&extraelec_veto<0.5&&extramuon_veto<0.5";
  TString isZTT = "&&((gen_match_1==5&&gen_match_2==5)||(gen_match_2==5&&gen_match_1==5))";
  TString isZLL = "&&((gen_match_1==1&&gen_match_2==1)||(gen_match_2==2&&gen_match_1==2))";
  TH1D* h_Bkg[6];TH1D* h_Bkg_SS[6];double x_sec[5]={5765,5765,61526,83.31,5765};
  int fill_color[5] = {42,8,2,3,4};
  TString cut[6];TString cutSS[6];
  cut[0]=Weight+zptmassweight+"(os>0.5"+isZTT+Cuts+")";cutSS[0]=Weight+"(os<0.5"+isZTT+Cuts+")";
  cut[1]=Weight+"(os>0.5&&((gen_match_1==5&&gen_match_2==6)||(gen_match_2==5&&gen_match_1==6))"+Cuts+")";
  cutSS[1]=Weight+"(os>0.5&&((gen_match_1==5&&gen_match_2==6)||(gen_match_2==5&&gen_match_1==6))"+Cuts+")";
  cut[2]=Weight+"(os>0.5"+Cuts+")";cutSS[2]=Weight+"(os<0.5"+Cuts+")";
  cut[3]=Weight+"(os>0.5"+Cuts+")";cutSS[3]=Weight+"(os<0.5"+Cuts+")";
  cut[4]=Weight+"(os>0.5"+isZLL+Cuts+")";cutSS[4]=Weight+"(os<0.5"+isZLL+Cuts+")";
  for(int ifile =0;ifile<5;ifile++){
    TFile* f_Bkg = new TFile(file_names[ifile]);
    TString histname = file_names[ifile];
    TH1D * nWeightedEvents = (TH1D*)f_Bkg->Get("nWeightedEvents"); //inputEventsH
    double norm = x_sec[ifile]*lumi/nWeightedEvents->GetSumOfWeights();
    h_Bkg[ifile]= new TH1D("h_Bkg","",nBins,xmin,xmax);
    h_Bkg_SS[ifile]= new TH1D("h_Bkg_SS","",nBins,xmin,xmax); 
    TTree* tree_ = (TTree*)f_Bkg ->Get("TauCheck");
    tree_->Draw(variable+">>h_Bkg",cut[ifile]);
    tree_->Draw(variable+">>h_Bkg_SS",cutSS[ifile]);
    h_Bkg[ifile]->Scale(norm);
    h_Bkg_SS[ifile]->Scale(norm);
  }
  
  // ********************************************
  // ***** Stitching of Drell-Yan samples *******
  // ********************************************
  TString refSamples[6];
  double refXSec[6];
  double refEvents[6] = {0,0,0,0,0,0};
  int nSamplesDY = 9;
  
  TH1D * histZtt[9];
  TH1D * histZll[9];
  TH1D * histZttSS[9];
  refSamples[0] = "DYJetsToLL_tt_Sync";
  refSamples[1] = "DY1JetsToLL_tt_Sync";
  refSamples[2] = "DY2JetsToLL_tt_Sync";
  refSamples[3] = "DY3JetsToLL_tt_Sync";
  refSamples[4] = "DY4JetsToLL_tt_Sync";

  refXSec[0] = 5765;
  refXSec[1] = 1.164*1012.5;
  refXSec[2] = 1.164*332.8;
  refXSec[3] = 1.164*101.8;
  refXSec[4] = 1.164*54.8;

  for (int iDY=0; iDY<5; ++iDY) {
    TFile * file = new TFile(refSamples[iDY]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iDY] = nWeightedEvents->GetSumOfWeights();
  }
  TString npartonCutsDY[9] = {"&&(gen_noutgoing==0||gen_noutgoing>4)", //cut on inclusive sample
  			      "&&gen_noutgoing==1",//cut on inclusive sample
  			      "&&gen_noutgoing==2",//cut on inclusive sample
  			      "&&gen_noutgoing==3",//cut on inclusive sample
  			      "&&gen_noutgoing==4",//cut on inclusive sample
  			      "",
  			      "",
  			      "",
  			      ""
  };
  TString dySampleNames[9] = {"DYJetsToLL_tt_Sync",
  			      "DYJetsToLL_tt_Sync",
  			      "DYJetsToLL_tt_Sync",
  			      "DYJetsToLL_tt_Sync",
  			      "DYJetsToLL_tt_Sync",
  			      "DY1JetsToLL_tt_Sync",
  			      "DY2JetsToLL_tt_Sync",
  			      "DY3JetsToLL_tt_Sync",
  			      "DY4JetsToLL_tt_Sync"
  };
  double dyNorm[9];
  dyNorm[0] = lumi*refXSec[0]/refEvents[0];
  dyNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  dyNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);

  TString cutsZtt[9];TString cutsZll[9];
  TString cutsZttSS[9];
  for (int iDY=0; iDY<nSamplesDY; ++iDY) {
    cutsZtt[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
    cutsZll[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZLL+")";
    cutsZttSS[iDY] = Weight+zptmassweight+"(os<0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
  }
  // filling histograms for DY samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples
    //cout<<dySampleNames[i]<<endl;
    TFile * file = new TFile(dySampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = dyNorm[i];
    TString histNameZtt   = dySampleNames[i] + variable + "_ztt_os";
    TString histNameZll   = dySampleNames[i] + variable + "_zll_os";
    TString histNameZttSS = dySampleNames[i] + variable + "_ztt_ss";
    histZtt[i]   = new TH1D(histNameZtt,"",nBins,xmin,xmax);
    histZll[i]   = new TH1D(histNameZll,"",nBins,xmin,xmax);
    histZttSS[i] = new TH1D(histNameZttSS,"",nBins,xmin,xmax);
    tree -> Draw(variable+">>"+histNameZtt,  cutsZtt[i]);
    tree -> Draw(variable+">>"+histNameZll,  cutsZll[i]);
    tree -> Draw(variable+">>"+histNameZttSS,cutsZttSS[i]);
    histZtt[i]   -> Scale(norm);
    histZll[i]   -> Scale(norm);
    histZttSS[i] -> Scale(norm);
  }
  h_Bkg[0]   = histZtt[0];
  h_Bkg[4]   = histZll[0];
  //h_Bkg_SS[0] = histZttSS[0];
  for (int iDY=1; iDY<nSamplesDY; ++iDY) {
    h_Bkg[0]-> Add(h_Bkg[0],histZtt[iDY]);
    h_Bkg[4]-> Add(h_Bkg[4],histZll[iDY]);
    //h_Bkg_SS[0]-> Add(h_Bkg_SS[0],histZtt[iDY]);
  }
  
  // *****************************************
  // ***** Stitching of W+Jets samples *******
  // *****************************************
  TH1D * histW[10];
  TH1D * histWSS[10];
  int nSamplesW = 9;
    
  TString npartonCuts[9] = {"&&(gen_noutgoing==0||gen_noutgoing>4)", //cut on inclusive sample
  			    "&&gen_noutgoing==1",//cut on inclusive sample
  			    "&&gen_noutgoing==2",//cut on inclusive sample
  			    "&&gen_noutgoing==3",//cut on inclusive sample
  			    "&&gen_noutgoing==4",//cut on inclusive sample
  			    "",
  			    "",
  			    "",
  			    ""
  };
  // redefine reference cross sections and reference samples
  refSamples[0] = "WJetsToLNu_tt_Sync";
  refSamples[1] = "W1JetsToLNu_tt_Sync";
  refSamples[2] = "W2JetsToLNu_tt_Sync";
  refSamples[3] = "W3JetsToLNu_tt_Sync";
  refSamples[4] = "W4JetsToLNu_tt_Sync";
  refXSec[0] = Wnorm*61527;
  refXSec[1] = Wnorm*1.221*9644.5;
  refXSec[2] = Wnorm*1.221*3144.5;
  refXSec[3] = Wnorm*1.221*954.8;
  refXSec[4] = Wnorm*1.221*485.8;
  refEvents[0] = 0;
  refEvents[1] = 0;
  refEvents[2] = 0;
  refEvents[3] = 0;
  refEvents[4] = 0;
  for (int iW=0; iW<5; ++iW) {
    TFile * file = new TFile(refSamples[iW]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iW] = nWeightedEvents->GetSumOfWeights();   //number of events with amc@NLO weight
  }

  TString wSampleNames[9] = {"WJetsToLNu_tt_Sync",
  			     "WJetsToLNu_tt_Sync",
  			     "WJetsToLNu_tt_Sync",
  			     "WJetsToLNu_tt_Sync",
  			     "WJetsToLNu_tt_Sync",
  			     "W1JetsToLNu_tt_Sync",
  			     "W2JetsToLNu_tt_Sync",
  			     "W3JetsToLNu_tt_Sync",
  			     "W4JetsToLNu_tt_Sync"
  };

  double wNorm[9];
  wNorm[0] = lumi*refXSec[0]/refEvents[0];  //norm for inlcusive with all events
  wNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]); //add inclusive sample (all events with one parton) with W1Jet sample
  wNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  wNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  wNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);

  TString cutsW[9];
  TString cutsWSS[9];

  for (int iW=0; iW<9; ++iW) {
    cutsW[iW]   = Weight+"(os>0.5"+Cuts+npartonCuts[iW]+")"; //apply parton cuts
    cutsWSS[iW] = Weight+"(os<0.5"+Cuts+npartonCuts[iW]+")";
  }
  // filling histograms for WJets samples
  for (int i=0; i<nSamplesW; ++i) { // run over W+Jets samples
    //cout<< wSampleNames[i] <<endl;
    TFile * file = new TFile(wSampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = wNorm[i];

    TString histNameW   = wSampleNames[i] + variable + "_w_os";
    TString histNameWSS = wSampleNames[i] + variable + "_w_ss";
    histW[i]   = new TH1D(histNameW,"",nBins,xmin,xmax);
    histWSS[i] = new TH1D(histNameWSS,"",nBins,xmin,xmax);

    tree->Draw(variable+">>"+histNameW,  cutsW[i]); //fill histogram with cuts applied
    tree->Draw(variable+">>"+histNameWSS,cutsWSS[i]);

    histW[i]   -> Scale(norm);
    histWSS[i] -> Scale(norm);
  }
  h_Bkg[2]   = histW[0];
  h_Bkg_SS[2] = histWSS[0];
  for (int iDY=1; iDY<nSamplesW; ++iDY) {
    h_Bkg[2]-> Add(h_Bkg[2],histW[iDY]);
    h_Bkg_SS[2]-> Add(h_Bkg_SS[2],histWSS[iDY]);
  }
  /**************************************************************************************/
  // QCD Estimation
  /*--------------------------*/
  TH1F* hist_QCD=(TH1F*)histdata(variable,nBins,xmin,xmax,true)->Clone();
  TH1F* hist_DYSS=(TH1F*)h_Bkg_SS[0]->Clone();
  TH1F* hist_WSS=(TH1F*)h_Bkg_SS[2]->Clone();
  TH1F* hist_ttSS=(TH1F*)h_Bkg_SS[3]->Clone();
  // hist_QCD->Add(hist_DYSS,-1);
  // hist_QCD->Add(hist_WSS,-1);
  // hist_QCD->Add(hist_ttSS,-1);
  TH1F* hist_Data=(TH1F*)histdata(variable,nBins,xmin,xmax,false)->Clone();
  //hist_QCD->Scale(1.6);
  /**************************************************************************************/
  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************
  
  float N_OS = hist_QCD->GetEntries();
  float N_SS = hist_Data->GetEntries();
  cout<<"qcd sf "<<N_OS/N_SS<<endl;
  
}

void Plot_all(){
  Plot_tautau_sig("m_vis",30,0,250);
  // Plot_tautau_sig("jdeta",8,0,6);
  // Plot_tautau_sig("njets",7,0,7);
  // Plot_tautau_sig("nbtag",6,0,6);
  // Plot_tautau_sig("pt_tt",30,0,250);
  // Plot_tautau_sig("met",30,0,150);
  // Plot_tautau_sig("metphi",30,-3,3);
  // Plot_tautau_sig("bpt_1",30,0,150);
  // Plot_tautau_sig("bpt_2",30,0,150);
  // Plot_tautau_sig("jpt_2",30,0,150);
  // Plot_tautau_sig("jpt_1",30,0,150);
  // Plot_tautau_sig("pt_1",30,0,150);
  // Plot_tautau_sig("pt_2",30,0,150);
  // Plot_tautau_sig("eta_1",30,-3,3);
  // Plot_tautau_sig("eta_2",30,-3,3);
  // Plot_tautau_sig("dijetpt",35,0,150);
  // Plot_tautau_sig("mjj",35,0,200);
  // Plot_tautau_sig("pt_tt",35,20,300);
  // Plot_tautau_sig("m_sv",35,20,300);
}
