#include <iostream>
#include <vector>
#include <map>
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
void Plot_lept_mutau(TString Variable = "m_vis",
		     TString xtitle = "m_{vis} [GeV]",
		     int nBins  =   30,
		     float xmin =    0,
		     float xmax =  300,
		     TString Weight = "puweight*effweight*mcweight*",
		     TString Cuts = "&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>20&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5||xTrigger>0.5)",//&&mva17_2>0.5&&mt_1<60
		     TString ytitle = "Events",
		     TString DataFile = "DATA_SingleMuon",
		     TString directory = "/nfs/dust/cms/user/cardinia/HtoTauTau/CMSSW_9_4_0_patch1/src/DesyTauAnalyses/NTupleMaker/test/mutau/",
		     TString Suffix = "MuTau_",        // for name of pdf
		     TString suffix = "",                      // for name of pdf
		     TString category = "inclusive",
		     bool logY = false, 
		     //double lumi = 14350  //RunsBC
		     //double lumi = 13463  //MuF
		     //double lumi = 27835 //RunsBCDE
		     double lumi = 41465 //SingleMuon
		     ){
  //ModTDRStyle();
  
   bool blindData = false;  
   int nbMin = 4;    //bins used for blinding
   int nbMax = 11;
   bool plotLeg = true;
   int position = 0; // 0 - right, 1 - left, 2 - central
   bool showSignal = false;

   TH1::SetDefaultSumw2();
   //TH2::SetDefaultSumw2();
   
   //double lumi =42810;
   double TTnorm = 1.0;   //scale factors for normalization of ttbar and wjets
   //double Wnorm  = 1.28; //from 80 GeV < mt_1 < 150 GeV
   double Wnorm  = 1.033;

   //weights
   //TString topweight("topptweightRun2*");
   TString topweight("1*");
   //TString topweight("topptweight*");
  
   
   TString qcdweight("1.06*");//TO FIX
   //TString qcdweight("qcdweight*");
   // TString zptmassweight("zptmassweight*");
   //  TString qcdweight("(qcdweight*qcdweight/qcdweightup)*");
   //  TString qcdweight("qcdweight_nodzeta*");
   
  
   TString sampleNames[14] = {
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
      "ZZ"             // (13)
   };
   
   std::cout << "Sample Names" << std::endl;

   double xsec[14] = {1, // data (0)
                      5765.4,  // DY(50) (1)
                      5765.4,  // DY(50) (2)
                      Wnorm*61526.7,// WJets (3)
                      TTnorm*87.31,  // TT  (4)
                      TTnorm*380.1,  // TT Hadronic  (5)
                      //TTnorm*365.35,  // TT Semilept  (6)
                      TTnorm*364.4,  // TT Semilept  (6)
                      38.09,           // ST_tW_antitop (7)
                      38.09,           // ST_tW_top_5f (8)
                      80.95,           // ST_t-channel_antitop (9)
                      136.95,           // ST_t-channel_top (10)
                      63.21, // WW   (11)
                      22.82,  // WZ    (12)
                      10.32  // ZZ      (13)
                      //          0   // dummy 
   };     
  
   TString cuts[14];
   TString cutsSS[14];
   
   //Work in Progress
   TString isZTT="&&(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";
   //TString isZLL="&&((gen_match_1==1&&gen_match_2==1)||(gen_match_1==2&&gen_match_2==2))";
   TString isZLL="&&!(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))";
  
   if(category=="0jet_low")
      {
         Cuts += "&&pt_2>15 && pt_2<35 && njets==0 && dzeta>-35";
      }
   
   if(category=="0jet_high")
      {
         Cuts += "&&pt_2>35 && njets==0 && dzeta>-35";
      }
   
   if(category=="1jet_low")
      {
         Cuts += "&&pt_2>15 && pt_2<35 && (njets==1 || (njets==2 && mjj<500)) && dzeta>-35";
      }
   
   if(category=="1jet_high")
      {
         Cuts += "&&pt_2>35 && (njets==1 || (njets==2 && mjj<500)) && dzeta>-35";
      }
   
   if(category=="vbf_low")
      {
         Cuts += "&&pt_2>15 && njets==2 && mjj>500 && mjj<800 && dzeta>-10";
      }
   
   if(category=="vbf_high")
      {
         Cuts += "&&pt_2>15 && njets==2 && mjj>800 && dzeta>-10";
      }
   
   for (int i=0; i<14; ++i) {
      cuts[i] = Weight+"(os>0.5"+Cuts+")";
      cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
   }
   TString zptmassweight="1.0*";                  //TO DO: CHANGE WEIGHTs
   cuts[0] = "(os>0.5"+Cuts+"&&metFilters>0.5)"; //DATA
   //cuts[0] = "(os>0.5&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>20&&mva_2>0.5&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5))"; //DATA
   cuts[0] = "(os>0.5"+Cuts+")"; //DATA
   cuts[1] = Weight+zptmassweight+"(os>0.5"+Cuts+isZTT+")";
   cuts[2] = Weight+zptmassweight+"(os>0.5"+Cuts+isZLL+")";
   cuts[4]  = Weight+topweight+"(os>0.5"+Cuts+")";
   cuts[5]  = Weight+topweight+"(os>0.5"+Cuts+")";
   cuts[6]  = Weight+topweight+"(os>0.5"+Cuts+")";
   
   cutsSS[0] = qcdweight+"(os<0.5"+Cuts+"&&metFilters>0.5)";
   //cutsSS[0] = qcdweight+"(os<0.5&&iso_1<0.15&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>20&&mva_2>0.5&&againstMuonTight3_2>0.5&&againstElectronVLooseMVA6_2>0.5&&(singleLepTrigger>0.5))";
   cutsSS[0] = qcdweight+"(os<0.5"+Cuts+")";
   cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZTT+")";
   cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+isZLL+")";
   cutsSS[4]  = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
   cutsSS[5]  = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
   cutsSS[6]  = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";


   TH1D * hist[14];
   TH1D * histSS[14];

   int nSamples = 14;
   
   TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
   
   // filling histograms
   for (int i=0; i<nSamples; ++i) {
     //std::cout << i << std::endl;

      std::cout << sampleNames[i] << std::endl;
      TFile * file = new TFile(directory+sampleNames[i]+".root");
      //std::cout << "file" << std::endl;

      TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");//inputEventsH
      TTree * tree = (TTree*)file->Get("TauCheck");
      std::cout << "tree" << std::endl;
      std::cout << "xsec: " << xsec[i] << std::endl;
      std::cout << "lumi: " << lumi << std::endl;

      double norm = xsec[i]*lumi/nWeightedEvents->GetSumOfWeights(); 
      std::cout << "norm: " << norm << std::endl;

      TString histName = sampleNames[i] + Variable + "_ss";
      TString histNameSS = sampleNames[i] + Variable + "_os";
      hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
      histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
      hist[i]->Sumw2();
      histSS[i]->Sumw2();
      std::cout << "Draw" << std::endl;
      tree->Draw(Variable+">>"+histName,cuts[i]);
      tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
      std::cout << sampleNames[i] << " " << hist[i]->GetEntries() << " " << hist[i]->Integral(0,nBins+1) << std::endl;
      if (i>0) {
         for (int iB=1; iB<=nBins; ++iB) {
            double x = hist[i]->GetBinContent(iB);
            double e = hist[i]->GetBinError(iB);
            hist[i]->SetBinContent(iB,norm*x);
            hist[i]->SetBinError(iB,norm*e);
            double xSS = histSS[i]->GetBinContent(iB);
            double eSS = histSS[i]->GetBinError(iB);
            histSS[i]->SetBinContent(iB,norm*xSS);
            histSS[i]->SetBinError(iB,norm*eSS);
         }
      }
   }
   
   delete dummyCanv;
   
   std::cout << "W+Jets samples" << std::endl;
    
    
    // *******************************
    // ***** W+Jets samples *******
    // *******************************
    
    TH1D * histW[10];
    TH1D * histWSS[10];
    int nSamplesDY = 9;
    
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
    // redefine reference cross sections
    // and reference samples
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
    
    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file = new TFile(directory+refSamples[iDY]+".root");
        TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
        refEvents[iDY] = nWeightedEvents->GetSumOfWeights();   //number of events with amc@NLO weight
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

    for (int iDY=0; iDY<9; ++iDY) {
       cutsW[iDY]   = Weight+"(os>0.5"+Cuts+npartonCuts[iDY]+")"; //apply parton cuts
       cutsWSS[iDY] = Weight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+")";
    }
    
    // filling histograms for WJets samples
    for (int i=0; i<nSamplesDY; ++i) { // run over samples
        TFile * file = new TFile(directory+wSampleNames[i]+".root");
        TTree * tree = (TTree*)file->Get("TauCheck");
        double norm = wNorm[i];
        TString histNameW   = wSampleNames[i] + Variable + "_w_os";
        TString histNameWSS = wSampleNames[i] + Variable + "_w_ss";
        histW[i]   = new TH1D(histNameW,"",nBins,xmin,xmax);
        histWSS[i] = new TH1D(histNameWSS,"",nBins,xmin,xmax);
        histW[i]->Sumw2();
        histWSS[i]->Sumw2();
        tree->Draw(Variable+">>"+histNameW,  cutsW[i]); //fill histogram with cuts applied
        tree->Draw(Variable+">>"+histNameWSS,cutsWSS[i]);
        for (int iB=1; iB<=nBins; ++iB)
        {
            double x = histW[i]->GetBinContent(iB);
            double e = histW[i]->GetBinError(iB);
            histW[i]->SetBinContent(iB,norm*x);   //set proper norm
            histW[i]->SetBinError(iB,norm*e);
	    //std::cout << "( "<< iB << " : " << norm*x << " ) ";
            x = histWSS[i]->GetBinContent(iB);
            e = histWSS[i]->GetBinError(iB);
            histWSS[i]->SetBinContent(iB,norm*x);
            histWSS[i]->SetBinError(iB,norm*e);
        }
        //std::cout << wSampleNames[i] << " -> W = " << histW[i]->GetEntries() << " : " << histW[i]->GetSumOfWeights() << " : " <<  norm << "    bin content 9: " << histW[i]->GetBinContent(9) << "   bin content 10: " << histW[i]->GetBinContent(10)
	//                  << std::endl;
        //    delete file;
    }
    hist[3]   = histW[0];
    histSS[3] = histWSS[0];
    
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[3]->Add(hist[3],histW[iDY]);
        histSS[3]->Add(histSS[3],histWSS[iDY]);
    }
    
    
    
    // *******************************
    // ***** Drell-Yan samples *******
    // *******************************
    
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
	//std::cout << iDY << ": " << cutsZtt[iDY] << "      "  << cutsZttSS[iDY] << "      "  << cutsZll[iDY] << "      "  << cutsZllSS[iDY] << "      " << std::endl;
    }
    
    nSamplesDY = 10;
    
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
        histZtt[i]->Sumw2();
        histZttSS[i]->Sumw2();
        histZll[i]->Sumw2();
        histZllSS[i]->Sumw2();
        tree->Draw(Variable+">>"+histNameZtt,  cutsZtt[i]);
        tree->Draw(Variable+">>"+histNameZttSS,cutsZttSS[i]);
        tree->Draw(Variable+">>"+histNameZll,  cutsZll[i]);
        tree->Draw(Variable+">>"+histNameZllSS,cutsZllSS[i]);
                
        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZtt[i]->GetBinContent(iB);
            double e = histZtt[i]->GetBinError(iB);
            histZtt[i]->SetBinContent(iB,norm*x);
            histZtt[i]->SetBinError(iB,norm*e);
            x = histZttSS[i]->GetBinContent(iB);
            e = histZttSS[i]->GetBinError(iB);
            histZttSS[i]->SetBinContent(iB,norm*x);
            histZttSS[i]->SetBinError(iB,norm*e);           
            
            x = histZll[i]->GetBinContent(iB);
            e = histZll[i]->GetBinError(iB);
            histZll[i]->SetBinContent(iB,norm*x);
            histZll[i]->SetBinError(iB,norm*e);
            x = histZllSS[i]->GetBinContent(iB);
            e = histZllSS[i]->GetBinError(iB);
            histZllSS[i]->SetBinContent(iB,norm*x);
            histZllSS[i]->SetBinError(iB,norm*e);
        }
        std::cout << dySampleNames[i] << " -> ZTT = " << histZtt[i]->GetEntries() << " : " << histZtt[i]->GetSumOfWeights()
        << "    ZLL = " << histZll[i]->GetEntries() << " : " << histZll[i]->GetSumOfWeights() << std::endl;
        //    delete file;
    }
    
    hist[1]   = histZtt[0];
    histSS[1] = histZttSS[0];
    hist[2]   = histZll[0];
    histSS[2] = histZllSS[0];
    
    for (int iDY=1; iDY<10; ++iDY) {
        hist[1]->Add(hist[1],histZtt[iDY]);
        hist[2]->Add(hist[2],histZll[iDY]);
        histSS[1]->Add(histSS[1],histZttSS[iDY]);
        histSS[2]->Add(histSS[2],histZllSS[iDY]);
    }

  // hist[1]->Add(hist[1],hist[3]);
  // hist[2]->Add(hist[2],hist[4]);

  //  adding up single top and VV backgrounds
  for (int iH=8; iH<14; ++iH) {
    hist[7]->Add(hist[7],hist[iH]);
    histSS[7]->Add(histSS[7],histSS[iH]);
  }
  //adding SS backgrounds
  for (int iH=2; iH<8; ++iH)
     histSS[1]->Add(histSS[1],histSS[iH]);
  

  //adding top
  hist[4]->Add(hist[4],hist[5]);
  hist[4]->Add(hist[4],hist[6]);

  float dataSS = histSS[0]->GetSumOfWeights();
  float dataSSfull = histSS[0]->Integral(0,nBins+1);
  
  // subtracting background from SS
  histSS[0]->Add(histSS[0],histSS[1],1,-1);
    
  float nonQCD = 
     histSS[1]->GetSumOfWeights();
  
  float nonQCDfull =
     histSS[1]->Integral(0,nBins+1);
  
  // float Wfraction     = histSS[5]->GetSumOfWeights()/dataSS;
  // float WfractionFull = histSS[5]->Integral(0,nBins+1)/dataSSfull;
  
  float nonQCDfraction = nonQCD/dataSS;
  float nonQCDfractionFull = nonQCDfull/dataSSfull;
  
  std::cout << "SS region" << std::endl;
  //std::cout << "VV (MC)      : " << histSS[7]->GetSumOfWeights() << " : "<< histSS[7]->Integral(0,nBins+1) << std::endl;
  std::cout << "W  (MC)      : " << histSS[4]->GetSumOfWeights() << " : "<< histSS[4]->Integral(0,nBins+1) << std::endl;
  // std::cout << "TT (MC)      : " << histSS[6]->GetSumOfWeights() << " : "<< histSS[6]->Integral(0,nBins+1) <<  std::endl;
  // std::cout << "DY (MC)      : " << histSS[1]->GetSumOfWeights() << " : "<< histSS[1]->Integral(0,nBins+1) << std::endl;
  std::cout << "non-QCD (MC) : " << nonQCD << " : " << nonQCDfull << std::endl;
  std::cout << "data         : " << dataSS << " : " << dataSSfull << std::endl;
  // std::cout << "W+Jets  fraction : " << Wfraction << " : " << WfractionFull << std::endl;
  std::cout << "non-QCD fraction : " << nonQCDfraction << " : " << nonQCDfractionFull << std::endl; 
  std::cout << std::endl;
 
 
  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
  TH1D * ZTT = (TH1D*)hist[1]->Clone("ZTT");
  TH1D * ZLL = (TH1D*)hist[2]->Clone("ZLL");
  TH1D * W   = (TH1D*)hist[3]->Clone("W");
  TH1D * TT  = (TH1D*)hist[4]->Clone("TT");
  TH1D * VV  = (TH1D*)hist[7]->Clone("VV");
  //TH1D * SMH = (TH1D*)hist[22]->Clone("SMH");
  for(int i=0;i<14;i++)std::cout << sampleNames[i] << " " << hist[i]->GetEntries() << " " << hist[i]->Integral(0,nBins+1) << std::endl;
 
 std::cout << "QCD : " << QCD->GetSumOfWeights() << " : " << QCD->Integral(1,nBins+1) << std::endl;
  std::cout << "VV  : " << VV->GetSumOfWeights() << " : " << VV->Integral(1,nBins+1) << std::endl;
  std::cout << "W   : " << W->GetSumOfWeights() << " : " << W->Integral(1,nBins+1) << std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : " << TT->Integral(1,nBins+1) << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << " : " << ZLL->Integral(1,nBins+1) << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << " : " << ZTT->Integral(1,nBins+1) << std::endl;

  float nData = histData->GetSumOfWeights();
  float nTT   = TT->GetSumOfWeights();
  float nW   = W->GetSumOfWeights();
  float eData = TMath::Sqrt(nData);
  float nNonTT = 
    VV->GetSumOfWeights() + 
    ZTT->GetSumOfWeights() +
    ZLL->GetSumOfWeights() + 
     QCD->GetSumOfWeights() +  
     W->GetSumOfWeights(); 
  float nNonW = 
    VV->GetSumOfWeights() + 
    ZTT->GetSumOfWeights() +
    ZLL->GetSumOfWeights() + 
     QCD->GetSumOfWeights() +  
     TT->GetSumOfWeights(); 

  float ttScale = (nData-nNonTT)/nTT;
  float ttScaleE = eData/nTT;
  float bkgE = 0.3*nNonTT/nTT;
  std::cout << "************************" << std::endl;
  std::cout << "TT scale factor = " << ttScale << " +/- " << ttScaleE << " +/- " << bkgE << std::endl;
  
  float WScale = (nData-nNonW)/nW;
  float WScaleE = eData/nW;
  float WbkgE = 0.3*nNonW/nW;
  
  std::cout << "W scale factor = " << WScale << " +/- " << WScaleE << " +/- " << WbkgE << std::endl;
  std::cout << "************************" << std::endl;
 

  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************
  TH1D * dummy = (TH1D*)ZTT->Clone("dummy");
  float errQCD = 0.15; // normalization of QCD background
  float errVV = 0.15; // normalization of VV background
  float errW = 0.1; // normalization of W+Jets background
  float errTT = 0.07; // normalization of TT background
  for (int iB=1; iB<=nBins; ++iB) {
     float eQCD = errQCD*QCD->GetBinContent(iB);
     float eVV = errVV*VV->GetBinContent(iB);
     float eW = errW*W->GetBinContent(iB);
     float eTT = errTT*TT->GetBinContent(iB);
     //float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eTT*eTT;
     float err2 = eQCD*eQCD+eVV*eVV + eW*eW + eTT*eTT;     //TO DO: WAS IST MIT DEM FEHLER AUF DY?
     std::cout<<"eQCD: "<<eQCD<<"eVV: "<<eVV<<"eW: "<<eW<<"eTT: "<<eTT<<std::endl;
     float errTot = TMath::Sqrt(err2);  //w
     std::cout<<errTot<<std::endl;
     dummy->SetBinError(iB,errTot); 
     //SMH->SetBinError(iB,0);
  }
 
  VV->Add(VV,QCD);
  W->Add(W,VV);
  ZLL->Add(ZLL,W);
  TT->Add(TT,ZLL);
  ZTT->Add(ZTT,TT);
  std::cout << "BKG : " << ZTT->GetSumOfWeights() << " : " << ZTT->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/ZTT->GetSumOfWeights() << "+/-" 
	    << TMath::Sqrt(histData->GetSumOfWeights())/ZTT->GetSumOfWeights() << std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////
    //ModTDRStyle();
    TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
    ModTDRStyle();

  
    TCanvas* canv1 = new TCanvas("c1", "c1");
    canv1->cd();
    std::vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
    pads[0]->SetLogy(logY);
    
    std::vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
    h[0]->Draw();
    
    std::string units="";
    std::string xtitle_ = (std::string) xtitle;
    size_t pos = xtitle_.find("[");
    if(pos!=std::string::npos) {
        units = xtitle_.substr(pos+1, xtitle_.find("]") - pos -1 );
        xtitle_ = xtitle_.substr(0, pos);
    }
    
    pads[1]->cd();
    h[1]->Draw();
    SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.65, 1.35);
    StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
    h[1]->GetYaxis()->SetNdivisions(4);
    h[1]->GetXaxis()->SetTitleOffset(1.1);
    h[1]->GetXaxis()->SetNdivisions(505);
    h[1]->GetYaxis()->SetTitleOffset(1.7);
    pads[0]->cd();
    h[0]->GetYaxis()->SetTitleOffset(2.1);
    pads[1]->SetGrid(0,1);
    //it complains if the minimum is set to 0 and you try to set log y scale
    if(logY) h[0]->SetMinimum(1);
    pads[0]->cd();
    
    // Setup legend
    TLegend *legend = PositionedLegend(0.40, 0.30, 3, 0.03);
    legend->SetTextFont(42);
    
    histData->SetMarkerColor(1);
    histData->SetLineColor(1);
    histData->SetFillColor(1);
    histData->SetFillStyle(0);
    histData->SetLineWidth(2);
    histData->SetMarkerStyle(20);
    histData->SetMarkerSize(1.1);
    
    legend->AddEntry(histData, "Observed", "ple");
    
    InitHist(QCD,TColor::GetColor("#FFCCFF"));
    InitHist(ZLL,TColor::GetColor("#DE5A6A"));
    InitHist(TT,TColor::GetColor("#9999CC"));
    InitHist(VV,TColor::GetColor("#6F2D35"));
    //InitHist(ZLL,TColor::GetColor("#4496C8"));
    InitHist(ZTT,TColor::GetColor("#FFCC66"));
    InitHist(W,TColor::GetColor("#4496C8"));

    legend->AddEntry(ZTT,"Z#rightarrow #tau#tau","f");
    //legend->AddEntry(ZLL,"Z#rightarrow ll","f");
    legend->AddEntry(TT,"t#bar{t}","f");
    legend->AddEntry(ZLL,"Z#rightarrow #mu#mu/ee","f");
    legend->AddEntry(VV,"single top + diboson","f");
    legend->AddEntry(W,"W+jets","f");
    legend->AddEntry(QCD,"QCD","f");
    
    ZTT->Draw("sameh");
    //ZLL->Draw("sameh");
    TT->Draw("sameh");
    ZLL->Draw("sameh");
    W->Draw("sameh");
    VV->Draw("sameh");
    QCD->Draw("sameh");
    
    canv1->Update();
    
   // InitSignal(vbfH,2);
    //InitSignal(SMH,2);
    
    // if (showSignal)
    // {
    //     legend->AddEntry(SMH,"SM Higgs(125) #times 10","f");
    // }
    
    // if (showSignal)
    // {
    //     SMH->Draw("hsame");
    // }
    
    canv1->Update();
    
    if (blindData)
    {
        for (int iB=nbMin; iB<=nbMax; ++iB)
        {
            histData->SetBinContent(iB,-1);
            histData->SetBinError(iB,0);
        }
    }
    
    
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
        std::cout<<"eStat: "<<eStat<<" eLumi:"<<eLumi<<"eBkg: "<<eBkg<<std::endl;
    }

    
    
    bkgdErr->SetMarkerSize(0);
    int new_idx = 923;//CreateTransparentColor(13,1.0);
    bkgdErr->SetFillColor(new_idx);
    bkgdErr->SetFillStyle(3004);
    bkgdErr->SetLineWidth(1);
    bkgdErr->Draw("e2same");
    legend->AddEntry(bkgdErr, "Bkg. uncertainty" , "F" );
    canv1->Update();

    TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
    TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    
    for (int iB=1; iB<=nBins; ++iB) {
        float x1 = histData->GetBinContent(iB);
        float x2 = ZTT->GetBinContent(iB);
        ratioErrH->SetBinContent(iB,1.0);
        ratioErrH->SetBinError(iB,0.0);
        float xBkg = bkgdErr->GetBinContent(iB);
        float errBkg = bkgdErr->GetBinError(iB);
        if (xBkg>0) {
            float relErr = errBkg/xBkg;
            ratioErrH->SetBinError(iB,relErr);
            
        }
        if (x1>0&&x2>0) {
            float e1 = histData->GetBinError(iB);
            float ratio = x1/x2;
            float eratio = e1/x2;
            ratioH->SetBinContent(iB,ratio);
            ratioH->SetBinError(iB,eratio);
        }
        else {
            ratioH->SetBinContent(iB,1000);
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
    
    canv1->Print("/nfs/dust/cms/user/cardinia/HtoTauTau/CMSSW_9_4_0_patch1/src/DesyTauAnalyses/NTupleMaker/test/Plots/"+Suffix+DataFile+"_"+Variable+suffix+".pdf");
    canv1->Print("/nfs/dust/cms/user/cardinia/HtoTauTau/CMSSW_9_4_0_patch1/src/DesyTauAnalyses/NTupleMaker/test/Plots/"+Suffix+DataFile+"_"+Variable+suffix+".eps");
}
