#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"


void FixObservable(int Observable, TString& Variable ,TString& xtitle,int& nBins,float& xmin,float& xmax);



void Plot_lept_mutau_Updated_Seq(int Observable,
			     TString Weight = "puweight*effweight*mcweight*",
			     TString Cuts = "&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>30&&mva17_2>0.5&&mt_1<50&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5||xTrigger>0.5)",//&&mva17_2>0.5&&mt_1<60&&(m_vis>60&&m_vis<90)
			     TString ytitle = "Events",
			     TString DataFile = "DATA_SingleMuon",//"SingleMuon_Run2017",
			     TString directory = "./mutau_2019_5_9_SVFit_DijetpT/",//"./mutau_2019_4_8/",
			     TString outputDir = "./Plots/",
			     TString Suffix = "MuTau_",        // for name of pdf
			     TString suffix = "_nocut",                      // for name of pdf
			     bool logY = false, 
			     //double lumi = 14350  //RunsBC
			     //double lumi = 13463  //MuF
			     //double lumi = 27835 //RunsBCDE
			     double lumi = 41465 //SingleMuon
			     ){

  TString Variable = "m_vis";
  TString xtitle = "m_{vis} [GeV]";
  int nBins  = 30;
  float xmin =    0;
  float xmax =  300;
  cout<<"Observable"<<Observable<<endl;
  FixObservable(Observable, Variable, xtitle, nBins, xmin,xmax);

  cout<<"Variable "<<Variable<<endl;     
  cout<<"xtitle "<<xtitle<<endl;

  

  TH1::SetDefaultSumw2();
  SetStyle();
  const int nSamples = 16;

  // some settings
  bool blindData = false;  
  int nbMin = 4;    //bins used for blinding
  int nbMax = 11;
  bool plotLeg = true;
  int position = 0; // 0 - right, 1 - left, 2 - central
  bool showSignal = true;

  // scale factors
  double TTnorm = 1.0;   //scale factors for normalization of ttbar and wjets
  double Wnorm  = 1.0;  // TO DO: this need to be determined !!!

  // weights
  TString topweight("1*");   
  TString qcdweight("1.06*");
  TString zptmassweight="1.0*";                  //TO DO: CHANGE WEIGHTs
   
  // These samples names should match the root file names (-> sampleNames[i].root is read in later)
  TString sampleNames[16] = {
        DataFile, // data (0)
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",// (1)Drell-Yan Z->TT
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (2)Drell-Yan Z->LL
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (3)WJets
        "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8",//(4)TTbar leptonic
        "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8",//(5) hadronic
        "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8",//(6) semileptonic
        "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (7) SingleTop tW tbar
        "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (8) SingleTop tW t
        "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",// (9) SingleTop t antitop
        "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",// (10) SingleTop t top
        "WW_TuneCP5_13TeV-pythia8",// (11)WW
        "WZ_TuneCP5_13TeV-pythia8",// (12)WZ
        "ZZ_TuneCP5_13TeV-pythia8",// (13)ZZ
	"GluGluHToTauTau_M125_13TeV_powheg_pythia8", // (14) Scalar
	"GluGluHToTauTau_M125_13TeV_powheg_pythia8" // (14) Scalar
	//	"HPseudoScalar" // (15) Pseudoscalar 
  };

  // Corresponding cross sections
  double xsec[16] = {1, // data (0)
		     5765.4,  // DY(50) (1)
		     5765.4,  // DY(50) (2)
		     Wnorm*61526.7,// WJets (3)
		     TTnorm*87.31,  // TT  (4)
		     TTnorm*380.1,  // TT Hadronic  (5)
		     TTnorm*364.4,  // TT Semilept  (6)
		     38.09,           // ST_tW_antitop (7)
		     38.09,           // ST_tW_top_5f (8)
		     80.95,           // ST_t-channel_antitop (9)
		     136.95,           // ST_t-channel_top (10)
		     63.21, // WW   (11)
		     22.82,  // WZ    (12)
		     10.32,  // ZZ      (13)
		     3.05*50,//43.92*0.0632*100,  // signal gg->Higgs times 10 ! (14)
		     3.05*50//43.92*0.0632*100  // signal gg->Higgs times 10 ! (15)
  };     


  // *******************************
  // ***** Selection Cuts    *******
  // *******************************

  TString cuts[16];
  TString cutsSS[16];

  // MC specific cuts to select certain type of particle
  TString isZTT="&&(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";
  TString isZLL="&&!(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";

  // Selection cuts applied to all samples
  for (int i=0; i<nSamples; ++i) {
    cuts[i]   = Weight+"(os>0.5"+Cuts+")";
    cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
  }

  // Some special cuts/weights applied to only a few samples
  cuts[0] = "(os>0.5"+Cuts+"&&metFilters>0.5)"; //DATA
  cuts[0] = "(os>0.5"+Cuts+")"; //DATA
  cuts[1] = Weight+zptmassweight+"(os>0.5"+Cuts+isZTT+")";
  cuts[2] = Weight+zptmassweight+"(os>0.5"+Cuts+isZLL+")";
  cuts[4] = Weight+topweight+"(os>0.5"+Cuts+")";
  cuts[5] = Weight+topweight+"(os>0.5"+Cuts+")";
  cuts[6] = Weight+topweight+"(os>0.5"+Cuts+")";
   
  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+"&&metFilters>0.5)";
  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZTT+")";
  cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZLL+")";
  cutsSS[4] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[5] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[6] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";

  // *******************************
  // ***** Filling Histograms ******
  // *******************************

  TH1D * hist[16];
  TH1D * histSS[16];


  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  // Draw main selection for all histograms in sampleNames
  for (int i=0; i<nSamples; ++i) {

   // cout << endl << sampleNames[i] << ":" << endl;
    // Reading input file
    TFile * file = new TFile( directory + sampleNames[i] + ".root");

    // Get tree and one further important histogram from input file
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents"); //inputEventsH
    TTree * tree = (TTree*)file->Get("TauCheck");

    // Calculate normalization of this sample
    double norm = xsec[i]*lumi/nWeightedEvents->GetSumOfWeights(); 

/*
    cout << "xsec: " << xsec[i] << endl;
    cout << "lumi: " << lumi << endl;
    cout << "norm: " << norm << endl;*/

    // Name and initialize histograms
    TString histName   = sampleNames[i] + Variable + "_ss";
    TString histNameSS = sampleNames[i] + Variable + "_os";
    hist[i]   = new TH1D(histName,"",nBins,xmin,xmax);
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);

   // cout << "Drawing ..." << endl;
    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);

    if (i>0) // if sample is MC sample -> Scale to xsec and luminosity
      {
	hist[i]   -> Scale(norm);
	histSS[i] -> Scale(norm);
	cout << "sample " << sampleNames[i] << "  norm: "<<norm << endl;
      }
	//    cout << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }
  cout << endl;
  delete dummyCanv;


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

  TString refSamples[6];
  double refXSec[6];
  double refEvents[6] = {0,0,0,0,0,0};

  // redefine reference cross sections and reference samples
  refSamples[0] = "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[1] = "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[2] = "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[3] = "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[4] = "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8";
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
    TFile * file = new TFile(directory+refSamples[iW]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iW] = nWeightedEvents->GetSumOfWeights();   //number of events with amc@NLO weight
  }

  TString wSampleNames[9] = {"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
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
    cutsWSS[iW] = Weight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iW]+")";
  }

  // filling histograms for WJets samples
  for (int i=0; i<nSamplesW; ++i) { // run over W+Jets samples
    cout<< wSampleNames[i] <<endl;
    TFile * file = new TFile(directory+wSampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = wNorm[i];

    TString histNameW   = wSampleNames[i] + Variable + "_w_os";
    TString histNameWSS = wSampleNames[i] + Variable + "_w_ss";
    histW[i]   = new TH1D(histNameW,"",nBins,xmin,xmax);
    histWSS[i] = new TH1D(histNameWSS,"",nBins,xmin,xmax);

    tree->Draw(Variable+">>"+histNameW,  cutsW[i]); //fill histogram with cuts applied
    tree->Draw(Variable+">>"+histNameWSS,cutsWSS[i]);

    histW[i]   -> Scale(norm);
    histWSS[i] -> Scale(norm);
  }

  hist[3]   = histW[0];
  histSS[3] = histWSS[0];

  for (int iW=1; iW<9; ++iW)
    {
      hist[3]   -> Add(hist[3],histW[iW]);
      histSS[3] -> Add(histSS[3],histWSS[iW]);
    }


  // ********************************************
  // ***** Stitching of Drell-Yan samples *******
  // ********************************************
  int nSamplesDY = 9;

  TH1D * histZtt[9];
  TH1D * histZttSS[9];
  TH1D * histZll[9];
  TH1D * histZllSS[9];

  refSamples[0] = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[1] = "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[2] = "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[3] = "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8";
  refSamples[4] = "DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8";
  //refSamples[5] = "DYJetsToLL_M-10to50";

  refXSec[0] = 5765;
  refXSec[1] = 1.164*1012.5;
  refXSec[2] = 1.164*332.8;
  refXSec[3] = 1.164*101.8;
  refXSec[4] = 1.164*54.8;
  //refXSec[5] = 15820;

  for (int iDY=0; iDY<5; ++iDY) {
    TFile * file = new TFile(directory+refSamples[iDY]+".root");
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

  TString dySampleNames[9] = {"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
			       "DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"//,
			      //"DYJetsToLL_M-10to50"
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
  //dyNorm[9] = lumi*refXSec[5]/refEvents[5];

  TString cutsZtt[9];
  TString cutsZttSS[9];
  TString cutsZll[9];
  TString cutsZllSS[9];

  for (int iDY=0; iDY<nSamplesDY; ++iDY) {
    cutsZtt[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
    cutsZttSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
    cutsZll[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZLL+")";
    cutsZllSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCutsDY[iDY]+isZLL+")";
  }


  // filling histograms for DY samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples
    cout<<dySampleNames[i]<<endl;
    TFile * file = new TFile(directory+dySampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = dyNorm[i];

    TString histNameZtt   = dySampleNames[i] + Variable + "_ztt_os";
    TString histNameZttSS = dySampleNames[i] + Variable + "_ztt_ss";
    TString histNameZll   = dySampleNames[i] + Variable + "_zll_os";
    TString histNameZllSS = dySampleNames[i] + Variable + "_zll_ss";
    histZtt[i]   = new TH1D(histNameZtt,"",nBins,xmin,xmax);
    histZttSS[i] = new TH1D(histNameZttSS,"",nBins,xmin,xmax);
    histZll[i]   = new TH1D(histNameZll,"",nBins,xmin,xmax);
    histZllSS[i] = new TH1D(histNameZllSS,"",nBins,xmin,xmax);

    tree -> Draw(Variable+">>"+histNameZtt,  cutsZtt[i]);
    tree -> Draw(Variable+">>"+histNameZttSS,cutsZttSS[i]);
    tree -> Draw(Variable+">>"+histNameZll,  cutsZll[i]);
    tree -> Draw(Variable+">>"+histNameZllSS,cutsZllSS[i]);

    histZtt[i]   -> Scale(norm);
    histZttSS[i] -> Scale(norm);
    histZll[i]   -> Scale(norm);
    histZllSS[i] -> Scale(norm);

    // cout << dySampleNames[i] << " -> ZTT : Entries = " << histZtt[i]->GetEntries()
    // 	 << " : Sum of weights = " << histZtt[i]->GetSumOfWeights()
    // 	 << "    ZLL : Entries = " << histZll[i]->GetEntries() << " : Sum of weights = " << histZll[i]->GetSumOfWeights() << endl;
  }

  hist[1]   = histZtt[0];
  histSS[1] = histZttSS[0];
  hist[2]   = histZll[0];
  histSS[2] = histZllSS[0];

  for (int iDY=1; iDY<9; ++iDY) {
    hist[1]  -> Add(hist[1],histZtt[iDY]);
    hist[2]  -> Add(hist[2],histZll[iDY]);
    histSS[1]-> Add(histSS[1],histZttSS[iDY]);
    histSS[2]-> Add(histSS[2],histZllSS[iDY]);
  }

  // ********************************************
  // ***** Adding similar backgrounds     *******
  // ********************************************

  // Adding up single top and VV backgrounds
  for (int iH=8; iH<14; ++iH) {
    hist[7]->Add(hist[7],hist[iH]);
    histSS[7]->Add(histSS[7],histSS[iH]);
  }
  // Adding SS backgrounds
  for (int iH=2; iH<8; ++iH) histSS[1]->Add(histSS[1],histSS[iH]);
  // Adding top
  hist[4]->Add(hist[4],hist[5]);
  hist[4]->Add(hist[4],hist[6]);


  // ********************************************
  // ***** QCD background estimation      *******
  // ********************************************

  float dataSS     = histSS[0]->GetSumOfWeights();
  float dataSSfull = histSS[0]->Integral(0,nBins+1);
  
  // Subtracting background from SS
  histSS[0]->Add(histSS[0],histSS[1],1,-1);

  float nonQCD             = histSS[1]->GetSumOfWeights();
  float nonQCDfull         = histSS[1]->Integral(0,nBins+1);
  float nonQCDfraction     = nonQCD/dataSS;
  float nonQCDfractionFull = nonQCDfull/dataSSfull;

  cout << endl;
  cout << "SS region :    " << endl;
  cout << "W  (MC)      : " << histSS[4]->GetSumOfWeights() << " : "<< histSS[4]->Integral(0,nBins+1) << endl;
  cout << "non-QCD (MC) : " << nonQCD << " : " << nonQCDfull << endl;
  cout << "data         : " << dataSS << " : " << dataSSfull << endl;
  cout << "non-QCD fraction : " << nonQCDfraction << " : " << nonQCDfractionFull << endl; 
  cout << endl;

  // ************************************
  // ***** Summarize backgrounds  *******
  // ************************************

  TH1D * histData = (TH1D*)hist[0]   -> Clone("data_obs");
  TH1D * QCD      = (TH1D*)histSS[0] -> Clone("QCD");
  TH1D * ZTT      = (TH1D*)hist[1]   -> Clone("ZTT");
  TH1D * ZLL      = (TH1D*)hist[2]   -> Clone("ZLL");
  TH1D * W        = (TH1D*)hist[3]   -> Clone("W");
  TH1D * TT       = (TH1D*)hist[4]   -> Clone("TT");
  TH1D * VV       = (TH1D*)hist[7]   -> Clone("VV");
  TH1D * SMH      = (TH1D*)hist[14]  -> Clone("SMH");
  TH1D * BSMH     = (TH1D*)hist[15]  -> Clone("BSMH");
  for(int i=0;i<nSamples;i++){
    cout << setw(nSamples) << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }

  cout << endl;
  cout << "QCD : Sum of weights = " << QCD->GetSumOfWeights() << " : Integral = " << QCD->Integral(1,nBins+1) << endl;
  cout << "VV  : Sum of weights = " << VV->GetSumOfWeights()  << " : Integral = " << VV->Integral(1,nBins+1)  << endl;
  cout << "W   : Sum of weights = " << W->GetSumOfWeights()   << " : Integral = " << W->Integral(1,nBins+1)   << endl;
  cout << "TT  : Sum of weights = " << TT->GetSumOfWeights()  << " : Integral = " << TT->Integral(1,nBins+1)  << endl;
  cout << "ZLL : Sum of weights = " << ZLL->GetSumOfWeights() << " : Integral = " << ZLL->Integral(1,nBins+1) << endl;
  cout << "ZTT : Sum of weights = " << ZTT->GetSumOfWeights() << " : Integral = " << ZTT->Integral(1,nBins+1) << endl;

  float nData    = histData->GetSumOfWeights();
  float nTT      = TT->GetSumOfWeights();
  float nW       = W->GetSumOfWeights();
  float eData    = TMath::Sqrt(nData);
  float nNonTT   = VV->GetSumOfWeights() + ZTT->GetSumOfWeights() + ZLL->GetSumOfWeights() + QCD->GetSumOfWeights() + W->GetSumOfWeights();
  float nNonW    = VV->GetSumOfWeights() + ZTT->GetSumOfWeights() + ZLL->GetSumOfWeights() + QCD->GetSumOfWeights() + TT->GetSumOfWeights();
  float ttScale  = (nData-nNonTT)/nTT;
  float ttScaleE = eData/nTT;
  float bkgE     = 0.3*nNonTT/nTT;
  float WScale   = (nData-nNonW)/nW;
  float WScaleE  = eData/nW;
  float WbkgE    = 0.3*nNonW/nW;

  cout << endl;
  cout << "************************" << endl;
  cout << "TT scale factor = " << ttScale << " +/- " << ttScaleE << " +/- " << bkgE << endl;
  cout << "W scale factor = " << WScale << " +/- " << WScaleE << " +/- " << WbkgE << endl;
  cout << "************************" << endl;
  cout << endl;

  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************

  TH1D * dummy = (TH1D*)ZTT->Clone("dummy");
  float errQCD = 0.15; // ad-hoc sys uncertainty of QCD background
  float errVV  = 0.15; // ad-hoc sys uncertainty of VV background
  float errW   = 0.10; // ad-hoc sys uncertainty of W+Jets background
  float errTT  = 0.07; // ad-hoc sys uncertainty of TT background

  for (int iB=1; iB<=nBins; ++iB) {   // Add general systematic uncertainties to each bin as error
    float eQCD   = errQCD*QCD->GetBinContent(iB);
    float eVV    = errVV*VV->GetBinContent(iB);
    float eW     = errW*W->GetBinContent(iB);
    float eTT    = errTT*TT->GetBinContent(iB);
    float err2   = eQCD*eQCD+eVV*eVV + eW*eW + eTT*eTT;     //TO DO: WAS IST MIT DEM FEHLER AUF DY?
    float errTot = TMath::Sqrt(err2);
    cout << "eQCD: " << eQCD << "  eVV: " << eVV << "  eW: " << eW << "  eTT: " << eTT << "  eTotal: " << errTot << endl;
    dummy -> SetBinError(iB,errTot);
    SMH   -> SetBinError(iB,0);
    BSMH  -> SetBinError(iB,0);
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING 
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();
  bool LargeScale = true;  
  TCanvas* canv1;
  if(LargeScale)canv1 = new TCanvas("c1", "c1", 2000,800);
  else canv1 = new TCanvas("c1", "c1", 1200,1000);
  canv1->cd();
  vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
  pads[0]->SetLogy(logY);

  vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
  h[0] -> Draw();

  string units="";
  string xtitle_ = (string) xtitle;
  size_t pos = xtitle_.find("[");
  if(pos!=string::npos) {
    units = xtitle_.substr(pos+1, xtitle_.find("]") - pos -1 );
    xtitle_ = xtitle_.substr(0, pos);
  }

  pads[1] -> cd();
  h[1]    -> Draw();
  SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.65, 1.35);
  StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
  h[1] -> GetYaxis()->SetNdivisions(4);
  h[1] -> GetXaxis()->SetTitleOffset(0.95);
  h[1] -> GetXaxis()->SetNdivisions(505);
  h[1] -> GetYaxis()->SetTitleOffset(1.1);
  if(LargeScale)  h[1] -> GetYaxis()->SetTitleOffset(0.9);
  pads[0] -> cd();
  h[0] -> GetYaxis()->SetTitleOffset(1.6);
  if(LargeScale)h[0] -> GetYaxis()->SetTitleOffset(1);
  pads[1] -> SetGrid(0,1);
  //it complains if the minimum is set to 0 and you try to set log y scale
  if(logY) h[0] -> SetMinimum(1);
  pads[0] -> cd();

  // Setup legend
  TLegend *legend = PositionedLegend(0.25, 0.30, 3, 0.03);
  legend -> SetTextFont(42);

  histData -> SetMarkerColor(1);
  histData -> SetLineColor(1);
  histData -> SetFillColor(1);
  histData -> SetFillStyle(0);
  histData -> SetLineWidth(2);
  histData -> SetMarkerStyle(20);
  histData -> SetMarkerSize(1.1);

  InitHist(QCD,TColor::GetColor("#FFCCFF"));
  InitHist(ZLL,TColor::GetColor("#DE5A6A"));
  InitHist(TT,TColor::GetColor("#9999CC"));
  InitHist(VV,TColor::GetColor("#6F2D35"));
  InitHist(ZTT,TColor::GetColor("#FFCC66"));
  InitHist(W,TColor::GetColor("#4496C8"));

  legend -> AddEntry(histData, "Observed", "ple");
  legend -> AddEntry(ZTT,"Z#rightarrow #tau#tau","f");
  legend -> AddEntry(TT,"t#bar{t}","f");
  legend -> AddEntry(ZLL,"Z#rightarrow #mu#mu/ee","f");
  legend -> AddEntry(W,"W+jets","f");
  legend -> AddEntry(VV,"single top + diboson","f");
  legend -> AddEntry(QCD,"QCD","f");

  // Add all bkg contributions to one stack plot
  THStack *stack = new THStack("Background","");
  stack -> Add(QCD);
  stack -> Add(VV);
  stack -> Add(W);
  stack -> Add(ZLL);
  stack -> Add(TT);
  stack -> Add(ZTT);
  stack -> Draw("hsame");

  canv1->Update();

  InitSignal(SMH ,2);
  InitSignal(BSMH,4);
  if (showSignal)
    {
      legend->AddEntry(SMH,"SM Higgs(125) #times 50","f");
      SMH->Draw("hsame");
      legend->AddEntry(BSMH,"BSM Higgs(125) #times 50","f");
      BSMH->Draw("hsame");
    }

  canv1->Update();

  if (blindData)
    {
      for (int iB=nbMin; iB<=nbMax; ++iB)
        {
	  histData->SetBinContent(iB,-1);
	  histData->SetBinError(iB,0);
        }
    }

  // Initialize a histogram which adds all error up
  TH1D * bkgdErr = (TH1D*)stack->GetStack()->Last()->Clone("bkgdErr");
  float errLumi = 0.03;
  float errMuon = 0.03;
  float errElectron = 0.04;
  for (int iB=1; iB<=nBins; ++iB) {
    QCD->SetBinError(iB,0);
    VV->SetBinError(iB,0);
    TT->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    float eStat =  bkgdErr->GetBinError(iB);
    float X = bkgdErr->GetBinContent(iB);
    float eLumi = errLumi * X;
    float eMuon = errMuon * X;
    float eElectron = errElectron * X;
    float eBkg = dummy->GetBinError(iB);
    float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg+eMuon*eMuon+eElectron*eElectron);
    bkgdErr->SetBinError(iB,Err);
    cout << "eStat = " << eStat << " : eLumi = "<< eLumi <<" : eBkg = " << eBkg << endl;
  }

  bkgdErr -> SetMarkerSize(0);
  int new_idx = 923;//CreateTransparentColor(13,1.0);
  bkgdErr -> SetFillColor(new_idx);
  bkgdErr -> SetFillStyle(3004);
  bkgdErr -> SetLineWidth(1);
  bkgdErr -> Draw("e2same");
  legend  -> AddEntry(bkgdErr, "Bkg. uncertainty" , "F" );
  canv1   -> Update();

  TH1D * ratioH    = (TH1D*)histData -> Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr  -> Clone("ratioErrH");
  ratioH -> Divide((TH1D*)stack->GetStack()->Last()); // Divide by the sum of the THStack

  // Set error of MC bkg correctly in ratio
  for (int iB=1; iB<=nBins; ++iB) {
    ratioErrH -> SetBinContent(iB,1.0);
    ratioErrH -> SetBinError(iB,0.0);
    float xBkg   = bkgdErr -> GetBinContent(iB);
    float errBkg = bkgdErr -> GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);
    }
  }

  pads[1]->cd();
  ratioErrH->Draw("e2same");
  ratioH->Draw("pe0same");

  pads[0]->cd();
  histData->Draw("pesame");

  FixTopRange(pads[0], GetPadYMax(pads[0]), 0.115);
  //DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
  if (lumi== 28686) DrawTitle(pads[0], "28.7 fb^{-1} (13 TeV, 2017)", 3);
  else if (lumi== 13960) DrawTitle(pads[0], "14.0 fb^{-1} (13 TeV, 2017)", 3);
  else if (lumi== 41465) DrawTitle(pads[0], "41.5 fb^{-1} (13 TeV, 2017)", 3);
  else DrawTitle(pads[0], "42.8 fb^{-1} (13 TeV, 2017)", 3);
  DrawTitle(pads[0], "#scale[1.2]{#bf{CMS} Work in progress}", 1);
  FixBoxPadding(pads[0], legend, 0.05);
  //legend->SetNColumns(2);
  legend->Draw();
  FixOverlay();
  canv1->Update();
  pads[0]->GetFrame()->Draw();

  canv1 -> Print( outputDir + Suffix + DataFile + "_" + Variable + suffix + ".pdf" );
  canv1 -> Print( outputDir + Suffix + DataFile + "_" + Variable + suffix + ".eps" );
  canv1 -> Print( outputDir + Suffix + DataFile + "_" + Variable + suffix + ".png" );
}


void FixObservable(int Observable, TString& Variable ,TString& xtitle,int& nBins,float& xmin,float& xmax){
  cout<<"Observable from func"<<Observable<<endl;

  //ordered according to appearance in the config file
  if(Observable==0){Variable = "pt_2"; xtitle = "pt_2 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==1){Variable = "jpt_1"; xtitle = " jpt_1 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==2){Variable = "jpt_2"; xtitle = " jpt_2 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==3){Variable = "bpt_1"; xtitle = "bpt_1 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==4){Variable = "bpt_2"; xtitle = "bpt_2 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==5){Variable = "njets"; xtitle = "njets"; nBins  = 10; xmin=0;xmax=10;}
  if(Observable==6){Variable = "nbtag"; xtitle = "nbtag"; nBins  = 10; xmin=0;xmax=10;}
  if(Observable==7){Variable = "m_sv"; xtitle = "m_sv [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==8){Variable = "mt_1"; xtitle = "mt_1 [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==9){Variable = "pt_tt"; xtitle = "pt_tt [GeV]"; nBins  = 20; xmin=0;xmax=200;}
  if(Observable==10){Variable = "mjj"; xtitle = "mjj [GeV]"; nBins  = 30; xmin=0;xmax=300;}
  if(Observable==11){Variable = "jdeta"; xtitle = "j#Delta#eta"; nBins  = 12; xmin=0;xmax=6;}
  if(Observable==12){Variable = "dijetpt"; xtitle = "dijetpt [GeV]"; nBins  = 30; xmin=0;xmax=300;}

  //things from the AN:
 if(Observable==13){Variable = "m_vis"; xtitle = "m_{vis} [GeV]"; nBins  = 30; xmin=0;xmax=300;}
 if(Observable==14){Variable = "pt_1"; xtitle = "muon pt [GeV]"; nBins  = 30; xmin=0;xmax=300;}
 if(Observable==15){Variable = "eta_1"; xtitle = "muon eta [GeV]"; nBins  = 12; xmin=-6;xmax=6;}
 if(Observable==16){Variable = "eta_2"; xtitle = "tau eta [GeV]"; nBins  = 12; xmin=-6;xmax=6;}

 if(Observable==17){Variable = "met_sv"; xtitle = "met_sv [GeV]"; nBins  = 30; xmin=0;xmax=300;}
 if(Observable==18){Variable = "pt_sv"; xtitle = "pt_sv [GeV]"; nBins  = 30; xmin=0;xmax=300;}
 //?? pt_ttjj?
 
  // if(Observable==14){Variable = "mt_tot"; xtitle = " mt_tot [GeV]"; nBins  = 30; xmin=0;xmax=300;}

}
