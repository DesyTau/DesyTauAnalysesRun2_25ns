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


/*2018 8 31: Merijn adjusted Plot_lept_mutau_Updated.C to produce the qcdweight factor. Adjustments w.r.t. Plot_lept_mutau_Updated:
-Currently the qcdweight factor is set to 1
-Put iso_1>0.15
-At end of macro we produce the data-estimate of QCD background (i.e. data-simulated EW background). The factor is produced as the ratio OS/SS
-

*/

void Plot_lept_mutau_Updated_QCDSF(TString Variable = "m_vis",
			     TString xtitle = "p_{T} (muon) [GeV]",
			     int nBins  =   30,
			     float xmin =    0,
			     float xmax =  300,
			     TString Weight = "puweight*effweight*mcweight*",
			     TString Cuts = "&&iso_1>0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>20&&mva17_2>0.5&&mt_1<60&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5||xTrigger>0.5)",//&&mva17_2>0.5&&mt_1<60
			     TString ytitle = "Events",
			     TString DataFile = "DATA_SingleMuon",
			     TString directory = "./mutau/",
			     TString outputDir = "./mutau/output/",
			     TString Suffix = "MuTau_",        // for name of pdf
			     TString suffix = "",                      // for name of pdf
			     bool logY = false, 
			     //double lumi = 14350  //RunsBC
			     //double lumi = 13463  //MuF
			     //double lumi = 27835 //RunsBCDE
			     double lumi = 41465 //SingleMuon
			     ){

  TH1::SetDefaultSumw2();
  SetStyle();

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
//  TString qcdweight("1.06*"); //2018 8 30: Merijn corrected to 1.16. New, but preliminary value.s
//  TString qcdweight("1.156104*");
  TString qcdweight("1.164685*");
  TString zptmassweight="1.0*";                  //TO DO: CHANGE WEIGHTs
   
  // These samples names should match the root file names (-> sampleNames[i].root is read in later)
  TString sampleNames[15] = {
    DataFile, // data (0)
    "DYJetsToLL",     // isZTT  (1)
    "DYJetsToLL",     // !isZTT (2)
    "WJetsToLNu",      // (3)
    "TTTo2L2Nu",               // (4)
    "TTToHadronic",               // (5)
    "TTToSemiLeptonic", //(6)
    "ST_tW_antitop", // (7)
    "ST_tW_top",     // (8)
    "ST_t_antitop",//(9)
    "ST_t_top",//(10)
    "WW",  // (11)
    "WZ",     // (12)
    "ZZ",             // (13)
    "ggH_125" // (14) 
  };

  // Corresponding cross sections
  double xsec[15] = {1, // data (0)
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
		     43.92*0.0632*10  // signal gg->Higgs times 10 !!! (14)
  };     


  // *******************************
  // ***** Selection Cuts    *******
  // *******************************

  TString cuts[15];
  TString cutsSS[15];
	
  // MC specific cuts to select certain type of particle
  TString isZTT="&&(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";
  TString isZLL="&&!(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";

  // Selection cuts applied to all samples
  for (int i=0; i<15; ++i) {
    cuts[i]   = Weight+"(os>0.5"+Cuts+")";
    cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
  }

  // Some special cuts/weights applied to only a few samples
  cuts[0] = "(os>0.5"+Cuts+")"; //DATA
  cuts[1] = Weight+zptmassweight+"(os>0.5"+Cuts+isZTT+")";
  cuts[2] = Weight+zptmassweight+"(os>0.5"+Cuts+isZLL+")";
  cuts[4] = Weight+topweight+"(os>0.5"+Cuts+")";
  cuts[5] = Weight+topweight+"(os>0.5"+Cuts+")";
  cuts[6] = Weight+topweight+"(os>0.5"+Cuts+")";
   
  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZTT+")";
  cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZLL+")";
  cutsSS[4] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[5] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
  cutsSS[6] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";

  // *******************************
  // ***** Filling Histograms ******
  // *******************************

  TH1D * hist[15];
  TH1D * histSS[15];

  int nSamples = 15;

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  // Draw main selection for all histograms in sampleNames
  for (int i=0; i<nSamples; ++i) {

    cout << endl << sampleNames[i] << ":" << endl;
    // Reading input file
    TFile * file = new TFile( directory + sampleNames[i] + ".root");

    // Get tree and one further important histogram from input file
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents"); //inputEventsH
    TTree * tree = (TTree*)file->Get("TauCheck");

    // Calculate normalization of this sample
    double norm = xsec[i]*lumi/nWeightedEvents->GetSumOfWeights(); 

    cout << "xsec: " << xsec[i] << endl;
    cout << "lumi: " << lumi << endl;
    cout << "norm: " << norm << endl;

    // Name and initialize histograms
    TString histName   = sampleNames[i] + Variable + "_ss";
    TString histNameSS = sampleNames[i] + Variable + "_os";
    hist[i]   = new TH1D(histName,"",nBins,xmin,xmax);
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);

    cout << "Drawing ..." << endl;
    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);

    if (i>0) // if sample is MC sample -> Scale to xsec and luminosity
      {
	hist[i]   -> Scale(norm);
	histSS[i] -> Scale(norm);
      }

    cout << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
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
  refSamples[0] = "WJetsToLNu";
  refSamples[1] = "W1JetsToLNu";
  refSamples[2] = "W2JetsToLNu";
  refSamples[3] = "W3JetsToLNu";
  refSamples[4] = "W4JetsToLNu";
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

  TString wSampleNames[9] = {"WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "W1JetsToLNu",
			     "W2JetsToLNu",
			     "W3JetsToLNu",
			     "W4JetsToLNu"
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

  TH1D * histZtt[10];
  TH1D * histZttSS[10];
  TH1D * histZll[10];
  TH1D * histZllSS[10];

  refSamples[0] = "DYJetsToLL";
  refSamples[1] = "DY1JetsToLL";
  refSamples[2] = "DY2JetsToLL";
  refSamples[3] = "DY3JetsToLL";
  refSamples[4] = "DY4JetsToLL";
  refSamples[5] = "DYJetsToLL_M-10to50";

  refXSec[0] = 5765;
  refXSec[1] = 1.164*1012.5;
  refXSec[2] = 1.164*332.8;
  refXSec[3] = 1.164*101.8;
  refXSec[4] = 1.164*54.8;
  refXSec[5] = 15820;

  for (int iDY=0; iDY<6; ++iDY) {
    TFile * file = new TFile(directory+refSamples[iDY]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iDY] = nWeightedEvents->GetSumOfWeights();
  }

  TString npartonCutsDY[10] = {"&&(gen_noutgoing==0||gen_noutgoing>4)", //cut on inclusive sample
			       "&&gen_noutgoing==1",//cut on inclusive sample
			       "&&gen_noutgoing==2",//cut on inclusive sample
			       "&&gen_noutgoing==3",//cut on inclusive sample
			       "&&gen_noutgoing==4",//cut on inclusive sample
			       "",
			       "",
			       "",
			       "",
			       ""
  };

  TString dySampleNames[10] = {"DYJetsToLL",
			       "DYJetsToLL",
			       "DYJetsToLL",
			       "DYJetsToLL",
			       "DYJetsToLL",
			       "DY1JetsToLL",
			       "DY2JetsToLL",
			       "DY3JetsToLL",
			       "DY4JetsToLL",
			       "DYJetsToLL_M-10to50"
  };

  double dyNorm[10];
  dyNorm[0] = lumi*refXSec[0]/refEvents[0];
  dyNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  dyNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  dyNorm[9] = lumi*refXSec[5]/refEvents[5];

  TString cutsZtt[10];
  TString cutsZttSS[10];
  TString cutsZll[10];
  TString cutsZllSS[10];

  for (int iDY=0; iDY<10; ++iDY) {
    cutsZtt[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
    cutsZttSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCutsDY[iDY]+isZTT+")";
    cutsZll[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCutsDY[iDY]+isZLL+")";
    cutsZllSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCutsDY[iDY]+isZLL+")";
  }

  int nSamplesDY = 10;

  // filling histograms for DY samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples

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


	cout << dySampleNames[i] << " -> ZTT : Entries = " << histZtt[i]->GetEntries()
    	 << " : Sum of weights = " << histZtt[i]->GetSumOfWeights()
     	 << "    ZLL : Entries = " << histZll[i]->GetEntries() << " : Sum of weights = " << histZll[i]->GetSumOfWeights() << endl;
  }

  hist[1]   = histZtt[0];
  histSS[1] = histZttSS[0];
  hist[2]   = histZll[0];
  histSS[2] = histZllSS[0];

  for (int iDY=1; iDY<10; ++iDY) {
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

  cout<<"data SS "<<dataSS<<endl;  

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

  //2018 8 28: Merijn: add all OS background processes:
  TH1D * AllOSBG=(TH1D*)ZTT-> Clone("AllOSBG");
	//cout<<"1 AllOSBG->GetSumOfWeights() "<< AllOSBG->GetSumOfWeights()<<endl;
   AllOSBG->Add(ZLL);
	//cout<<"2 AllOSBG->GetSumOfWeights() "<< AllOSBG->GetSumOfWeights()<<endl;
   AllOSBG->Add(W);
	//cout<<"3 AllOSBG->GetSumOfWeights() "<< AllOSBG->GetSumOfWeights()<<endl;
   AllOSBG->Add(TT);
	//cout<<"4 AllOSBG->GetSumOfWeights() "<< AllOSBG->GetSumOfWeights()<<endl;
   AllOSBG->Add(VV);
cout<<"All Opposite sign BG GetSumOfWeights() "<< AllOSBG->GetSumOfWeights()<<endl;

  //make copy of OS data:
  TH1D * DATAQCDOS = (TH1D*)hist[0]   -> Clone("DATAQCDOS");
//cout<<"1 DATAQCDOS->GetSumOfWeights() "<< DATAQCDOS->GetSumOfWeights()<<endl;
  
  //subtract the OS background from data:
  DATAQCDOS->Add(AllOSBG,-1);  
  //cout<<"2 DATAQCDOS->GetSumOfWeights() "<< DATAQCDOS->GetSumOfWeights()<<endl;

  cout<<"We calculate the updatedremo ratio DATAQCDOS->GetSumOfWeights()/QCD->GetSumOfWeights() "<<DATAQCDOS->GetSumOfWeights()/QCD->GetSumOfWeights()<<endl;
cout<<endl;
  //cout<<"The QCD weight factor is given as DATAQCDOS->GetSumOfWeights()/QCD->GetSumOfWeights() "<<DATAQCDOS->GetSumOfWeights()/QCD->GetSumOfWeights() <<endl;

  for(int i=0;i<15;i++){
    cout << setw(15) << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }

  cout << endl;
  cout << "QCD : Sum of weights = " << QCD->GetSumOfWeights() << " : Integral = " << QCD->Integral(1,nBins+1) << endl;
  cout << "VV  : Sum of weights = " << VV->GetSumOfWeights()  << " : Integral = " << VV->Integral(1,nBins+1)  << endl;
  cout << "W   : Sum of weights = " << W->GetSumOfWeights()   << " : Integral = " << W->Integral(1,nBins+1)   << endl;
  cout << "TT  : Sum of weights = " << TT->GetSumOfWeights()  << " : Integral = " << TT->Integral(1,nBins+1)  << endl;
  cout << "ZLL : Sum of weights = " << ZLL->GetSumOfWeights() << " : Integral = " << ZLL->Integral(1,nBins+1) << endl;
  cout << "ZTT : Sum of weights = " << ZTT->GetSumOfWeights() << " : Integral = " << ZTT->Integral(1,nBins+1) << endl;
  cout<<"QCD OS : Sum of weights = " << DATAQCDOS->GetSumOfWeights() << " : Integral = " << DATAQCDOS->Integral(1,nBins+1) << endl;

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
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING 
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();

  TCanvas* canv1 = new TCanvas("c1", "c1");
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
  h[1] -> GetXaxis()->SetTitleOffset(1.1);
  h[1] -> GetXaxis()->SetNdivisions(505);
  h[1] -> GetYaxis()->SetTitleOffset(1.7);
  pads[0] -> cd();
  h[0] -> GetYaxis()->SetTitleOffset(2.1);
  pads[1] -> SetGrid(0,1);
  //it complains if the minimum is set to 0 and you try to set log y scale
  if(logY) h[0] -> SetMinimum(1);
  pads[0] -> cd();

  // Setup legend
  TLegend *legend = PositionedLegend(0.40, 0.30, 3, 0.03);
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

  InitSignal(SMH,2);
  if (showSignal)
    {
      legend->AddEntry(SMH,"SM Higgs(125) #times 10","f");
      SMH->Draw("hsame");
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
  DrawTitle(pads[0], "#mu#tau", 1);
  FixBoxPadding(pads[0], legend, 0.05);
  legend->Draw();
  FixOverlay();
  canv1->Update();
  pads[0]->GetFrame()->Draw();

  canv1 -> Print( outputDir + Suffix + DataFile + "_" + Variable + suffix + ".pdf" );
  canv1 -> Print( outputDir + Suffix + DataFile + "_" + Variable + suffix + ".eps" );
}
