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
void Plot(TString directory = "/nfs/dust/cms/user/rasp/Run/Run2016/EMu_2016BCD/",
	  TString suffix = "",
	  TString DataFile = "MuonEG_Run2016",
	  TString Variable = "mTtot",
	  TString MA =  "300",
	  TString Suffix = "_mttot_0jetlow",
	  int nBins  =   30,
	  float xmin =    0,
	  float xmax =  300,
	  TString Weight = "mcweight*puweight*effweight*",
	  //TString Cuts = "&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",
      TString Cuts = "&&dzeta_mvamet>-20&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",
	  //	  TString Cuts = "&&dzeta_mvamet<-60&&mvamet>80&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",
	  TString xtitle = "m_{T}^{tot} [GeV]",
	  TString ytitle = "Events",
      TString category = "inclusive",
      bool logY = false
          )
{
  //ModTDRStyle();

  bool blindData = true;
  int nbMin = 4;
  int nbMax = 11;
  bool plotLeg = true;
  int position = 0; // 0 - right, 1 - left, 2 - central
  bool showSignal = true;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  double lumi = 12892;
  double TTnorm = 1.0;
  double Wnorm  = 1.0;

  //  TString topweight("");
  TString topweight2("");
  TString topweight("topptweight*");
  //  TString topweight2("topptweight*topptweight*");

  //  TString qcdweight("2.02*");
  TString qcdweight("qcdweight*");
  TString zptmassweight("zptmassweight*");
  //  TString qcdweight("(qcdweight*qcdweight/qcdweightup)*");
  //  TString qcdweight("qcdweight_nodzeta*");

  TString sampleNames[30] = {
    DataFile, // data (0)
    "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",     // isZTT  (1)
    "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",     // !isZTT (2)
    "DYJetsToLL_M-10to50_13TeV-madgraphMLM", // isZTT  (3)
    "DYJetsToLL_M-10to50_13TeV-madgraphMLM", // !isZTT (4)
    "WJetsToLNu_13TeV-madgraphMLM",      // (5)
    "TTJets_13TeV-powheg",               // (6)
    "ST_t-channel_top_4f_leptonDecays_13TeV-powheg", // (7)
    "ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg",    // (8)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg", // (9)
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg",     // (10)
    "VVTo2L2Nu_13TeV_amcatnloFXFX",    // (11)
    "WWTo1L1Nu2Q_13TeV_amcatnloFXFX",  // (12)
    "WZTo2L2Q_13TeV_amcatnloFXFX",     // (13)
    "WZTo1L1Nu2Q_13TeV_amcatnloFXFX",  // (14)
    "WZTo1L3Nu_13TeV_amcatnloFXFX",    // (15)
    "WZTo3LNu_13TeV-powheg",           // (16)
    "ZZTo4L_13TeV_powheg",             // (17)
    "ZZTo2L2Q_13TeV_amcatnloFXFX",     // (18)
    "WGToLNuG_13TeV-madgraphMLM",             // (19)
    "WGstarToLNuMuMu_012Jets_13TeV-madgraph", // (20)
    "WGstarToLNuEE_012Jets_13TeV-madgraph",   // (21)
    "GluGluHToTauTau_M125_13TeV_powheg",    // (22)
    "VBFHToTauTau_M125_13TeV_powheg", // (23)
    ""    // (24);
  };


  double xsec[30] = {1, // data (0)
		     5765,  // DY(50) (1)
		     5765,  // DY(50) (2)
		     18610,   // DY(10to50)     (3)
		     18610,   // DY(10to50)     (4)
		     Wnorm*61526.7,// WJets (5)
		     TTnorm*831.76,  // TT  (6)
		     136.95*3*0.108, // ST_t-channel_top (7)
		     80.95*3*0.108,  // ST_t-channel_antitop (8)
		     35.6,           // ST_tW_antitop (9)
		     35.6,           // ST_tW_top_5f (10)
		     11.95,  // VV          (11)
		     49.997, // WWToLNuQQ   (12)
		     5.595,  // WZTo2L2Q    (13)
		     10.71,  // WZTo1L1Nu2Q (14)
		     3.05,   // WZTo1L3Nu   (15)
		     5.26,   // WZTo3L1Nu   (16)
		     1.212,  // ZZTo4L      (17)
		     3.22,   // ZZTo2L2Q    (18)
             489.0,  // WGToLNuG        (19)
             2.793,  // WGstarToLNuMuMu (20)
             3.526,  // WGstarToLNuEE   (21)
		     43.92*0.0632,  // gg->H (22)
		     3.748*0.0632,  // VBF H   (23)
		     0   // dummy 
  };     
  
  TString cuts[30];
  TString cutsSS[30];
  
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

  for (int i=0; i<30; ++i) {
    cuts[i] = Weight+"(os>0.5"+Cuts+")";
    cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
  }
  cuts[0] = "(os>0.5"+Cuts+"&&metFilters>0.5)";
  cuts[1] = Weight+zptmassweight+"(os>0.5"+Cuts+"&&isZTT)";
  cuts[2] = Weight+zptmassweight+"(os>0.5"+Cuts+"&&!isZTT)";
  cuts[3] = Weight+"(os>0.5"+Cuts+"&&isZTT)";
  cuts[4] = Weight+"(os>0.5"+Cuts+"&&!isZTT)";

  cuts[6]  = Weight+topweight+"(os>0.5"+Cuts+")";

  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+"&&metFilters>0.5)";
  cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";
  cutsSS[3] = Weight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[4] = Weight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";

  cutsSS[6]  = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";


  TH1D * hist[30];
  TH1D * histSS[30];

  int nSamples = 24;

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
  
  // filling histograms
  for (int i=0; i<nSamples; ++i) {
    //    std::cout << sampleNames[i] << std::endl;
    TFile * file = new TFile(directory+sampleNames[i]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
    TString histName = sampleNames[i] + Variable + "_ss";
    TString histNameSS = sampleNames[i] + Variable + "_os";
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
    hist[i]->Sumw2();
    histSS[i]->Sumw2();
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
    
    
    
    // *******************************
    // ***** W+Jets samples *******
    // *******************************
    
    TH1D * histW[10];
    TH1D * histWSS[10];
    int nSamplesDY = 9;
    
    TString npartonCuts[9] = {"&&(npartons==0||npartons>4)",
        "&&npartons==1",
        "&&npartons==2",
        "&&npartons==3",
        "&&npartons==4",
        "",
        "",
        "",
        ""
    };
    
    TString refSamples[5];
    double refXSec[5];
    double refEvents[5] = {0,0,0,0,0};
    // redefine reference cross sections
    // and reference samples
    refSamples[0] = "WJetsToLNu_13TeV-madgraphMLM";
    refSamples[1] = "W1JetsToLNu_13TeV-madgraphMLM";
    refSamples[2] = "W2JetsToLNu_13TeV-madgraphMLM";
    refSamples[3] = "W3JetsToLNu_13TeV-madgraphMLM";
    refSamples[4] = "W4JetsToLNu_13TeV-madgraphMLM";
    
    refXSec[0] = 61527;
    refXSec[1] = 1.221*9644.5;
    refXSec[2] = 1.221*3144.5;
    refXSec[3] = 1.221*954.8;
    refXSec[4] = 1.221*485.6;
    
    refEvents[0] = 0;
    refEvents[1] = 0;
    refEvents[2] = 0;
    refEvents[3] = 0;
    refEvents[4] = 0;
    
    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file = new TFile(directory+refSamples[iDY]+".root");
        TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
        refEvents[iDY] = histWeightsH->GetSumOfWeights();
    }
    TString wSampleNames[9] = {"WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "W1JetsToLNu_13TeV-madgraphMLM",
			     "W2JetsToLNu_13TeV-madgraphMLM",
			     "W3JetsToLNu_13TeV-madgraphMLM",
			     "W4JetsToLNu_13TeV-madgraphMLM"
    };

    double wNorm[9];
    wNorm[0] = lumi*refXSec[0]/refEvents[0];
    wNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
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
        cutsW[iDY]   = Weight+"(os>0.5"+Cuts+npartonCuts[iDY]+")";
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
        tree->Draw(Variable+">>"+histNameW,  cutsW[i]);
        tree->Draw(Variable+">>"+histNameWSS,cutsWSS[i]);
        
        for (int iB=1; iB<=nBins; ++iB)
        {
            double x = histW[i]->GetBinContent(iB);
            double e = histW[i]->GetBinError(iB);
            histW[i]->SetBinContent(iB,norm*x);
            histW[i]->SetBinError(iB,norm*e);
            x = histWSS[i]->GetBinContent(iB);
            e = histWSS[i]->GetBinError(iB);
            histWSS[i]->SetBinContent(iB,norm*x);
            histWSS[i]->SetBinError(iB,norm*e);
        }
        std::cout << wSampleNames[i] << " -> W = " << histW[i]->GetEntries() << " : " << histW[i]->GetSumOfWeights() 
        << std::endl;
        //    delete file;
    }
    hist[5]   = histW[0];
    histSS[5] = histWSS[0];
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[5]->Add(hist[5],histW[iDY]);
        histSS[5]->Add(histSS[5],histWSS[iDY]);
    }
    
    
    
    // *******************************
    // ***** Drell-Yan samples *******
    // *******************************
    
    TH1D * histZtt[10];
    TH1D * histZttSS[10];
    
    TH1D * histZll[10];
    TH1D * histZllSS[10];
    
    refSamples[0] = "DYJetsToLL_M-50_13TeV-madgraphMLM_ext";
    refSamples[1] = "DY1JetsToLL_M-50_13TeV-madgraphMLM";
    refSamples[2] = "DY2JetsToLL_M-50_13TeV-madgraphMLM";
    refSamples[3] = "DY3JetsToLL_M-50_13TeV-madgraphMLM";
    refSamples[4] = "DY4JetsToLL_M-50_13TeV-madgraphMLM";

    refXSec[0] = 5765;
    refXSec[1] = 1.164*1012.5;
    refXSec[2] = 1.164*332.8;
    refXSec[3] = 1.164*101.8;
    refXSec[4] = 1.164*54.8;
    
    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file = new TFile(directory+refSamples[iDY]+".root");
        TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
        refEvents[iDY] = histWeightsH->GetSumOfWeights();
    }
    TString dySampleNames[9] = {"DYJetsToLL_M-50_13TeV-madgraphMLM_ext",
        "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",
        "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",
        "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",
        "DYJetsToLL_M-50_13TeV-madgraphMLM_ext",
        "DY1JetsToLL_M-50_13TeV-madgraphMLM",
        "DY2JetsToLL_M-50_13TeV-madgraphMLM",
        "DY3JetsToLL_M-50_13TeV-madgraphMLM",
        "DY4JetsToLL_M-50_13TeV-madgraphMLM"
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
    
    
    TString cutsZtt[9];
    TString cutsZttSS[9];
    
    TString cutsZll[9];
    TString cutsZllSS[9];
    
    for (int iDY=0; iDY<9; ++iDY) {
        cutsZtt[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCuts[iDY]+"&&isZTT>0.5)";
        cutsZttSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+"&&isZTT>0.5)";
        cutsZll[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCuts[iDY]+"&&isZTT<0.5)";
        cutsZllSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+"&&isZTT<0.5)";
    }
    
    nSamplesDY = 9;
    
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
    
    
    for (int iDY=1; iDY<9; ++iDY) {
        hist[1]->Add(hist[1],histZtt[iDY]);
        hist[2]->Add(hist[2],histZll[iDY]);
        histSS[1]->Add(histSS[1],histZttSS[iDY]);
        histSS[2]->Add(histSS[2],histZllSS[iDY]);
    }

    

  hist[1]->Add(hist[1],hist[3]);
  hist[2]->Add(hist[2],hist[4]);

  //  adding up single top and VV backgrounds
  for (int iH=8; iH<19; ++iH) {
    hist[7]->Add(hist[7],hist[iH]);
    histSS[7]->Add(histSS[7],histSS[iH]);
  }

  // adding up WJets and WGamma samples
    for (int iH = 19; iH < 22; ++iH)
    {
        hist[5]->Add(hist[5],hist[iH]);
        histSS[5]->Add(histSS[5],histSS[iH]);
    }
  // adding up SM Higgs!
 
    hist[22]->Add(hist[22],hist[23]);
    
    for (int iB=1; iB<=nBins; ++iB)
    {
        double x = hist[22]->GetBinContent(iB);
        hist[22]->SetBinContent(iB,10*x);
    }
    
  for (int iH=2; iH<5; ++iH)
    histSS[1]->Add(histSS[1],histSS[iH]);

  float dataSS = histSS[0]->GetSumOfWeights();
  float dataSSfull = histSS[0]->Integral(0,nBins+1);

  // subtracting background from SS
  histSS[0]->Add(histSS[0],histSS[1],1,-1);
  histSS[0]->Add(histSS[0],histSS[5],1,-1);
  histSS[0]->Add(histSS[0],histSS[6],1,-1);
  histSS[0]->Add(histSS[0],histSS[7],1,-1);

  float nonQCD = 
    histSS[7]->GetSumOfWeights() + 
    histSS[5]->GetSumOfWeights() + 
    histSS[6]->GetSumOfWeights() + 
    histSS[1]->GetSumOfWeights();

  float nonQCDfull =
    histSS[7]->Integral(0,nBins+1) +
    histSS[5]->Integral(0,nBins+1) +
    histSS[6]->Integral(0,nBins+1) +
    histSS[1]->Integral(0,nBins+1);

  float Wfraction     = histSS[5]->GetSumOfWeights()/dataSS;
  float WfractionFull = histSS[5]->Integral(0,nBins+1)/dataSSfull;

  float nonQCDfraction = nonQCD/dataSS;
  float nonQCDfractionFull = nonQCDfull/dataSSfull;

  std::cout << "SS region" << std::endl;
  std::cout << "VV (MC)      : " << histSS[7]->GetSumOfWeights() << " : "<< histSS[7]->Integral(0,nBins+1) << std::endl;
  std::cout << "W  (MC)      : " << histSS[5]->GetSumOfWeights() << " : "<< histSS[5]->Integral(0,nBins+1) << std::endl;
  std::cout << "TT (MC)      : " << histSS[6]->GetSumOfWeights() << " : "<< histSS[6]->Integral(0,nBins+1) <<  std::endl;
  std::cout << "DY (MC)      : " << histSS[1]->GetSumOfWeights() << " : "<< histSS[1]->Integral(0,nBins+1) << std::endl;
  std::cout << "non-QCD (MC) : " << nonQCD << " : " << nonQCDfull << std::endl;
  std::cout << "data         : " << dataSS << " : " << dataSSfull << std::endl;
  std::cout << "W+Jets  fraction : " << Wfraction << " : " << WfractionFull << std::endl;
  std::cout << "non-QCD fraction : " << nonQCDfraction << " : " << nonQCDfractionFull << std::endl; 
  std::cout << std::endl;

  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
  TH1D * ZTT = (TH1D*)hist[1]->Clone("ZTT");
  TH1D * ZLL = (TH1D*)hist[2]->Clone("ZLL");
  TH1D * W   = (TH1D*)hist[5]->Clone("W");
  TH1D * TT  = (TH1D*)hist[6]->Clone("TT");
  TH1D * VV  = (TH1D*)hist[7]->Clone("VV");
  TH1D * SMH = (TH1D*)hist[22]->Clone("SMH");

  std::cout << "QCD : " << QCD->GetSumOfWeights() << " : " << QCD->Integral(1,nBins+1) << std::endl;
  std::cout << "VV  : " << VV->GetSumOfWeights() << " : " << VV->Integral(1,nBins+1) << std::endl;
  std::cout << "W   : " << W->GetSumOfWeights() << " : " << W->Integral(1,nBins+1) << std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : " << TT->Integral(1,nBins+1) << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << " : " << ZLL->Integral(1,nBins+1) << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << " : " << ZTT->Integral(1,nBins+1) << std::endl;

  float nData = histData->GetSumOfWeights();
  float nTT   = TT->GetSumOfWeights();
  float eData = TMath::Sqrt(nData);
  float nNonTT = 
    VV->GetSumOfWeights() + 
    ZTT->GetSumOfWeights() +
    ZLL->GetSumOfWeights() + 
    QCD->GetSumOfWeights() +
    W->GetSumOfWeights(); 

  float ttScale = (nData-nNonTT)/nTT;
  float ttScaleE = eData/nTT;
  float bkgE = 0.3*nNonTT/nTT;

  //  std::cout << "TT scale factor = " << ttScale << " +/- " << ttScaleE << " +/- " << bkgE << std::endl;


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
    float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eTT*eTT;
    float errTot = TMath::Sqrt(err2);
    dummy->SetBinError(iB,errTot);
    SMH->SetBinError(iB,0);
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
    h[1]->GetXaxis()->SetTitleOffset(1.2);
    h[1]->GetYaxis()->SetTitleOffset(2.0);
    pads[0]->cd();
    h[0]->GetYaxis()->SetTitleOffset(2.0);
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
    
    legend->AddEntry(ZTT,"Z#rightarrow #tau#tau","f");
    //legend->AddEntry(ZLL,"Z#rightarrow ll","f");
    legend->AddEntry(TT,"t#bar{t}","f");
    legend->AddEntry(ZLL,"Electroweak","f");
    //  leg->AddEntry(VV,"VV+VVV","f");
    legend->AddEntry(QCD,"Misidentified e/#mu","f");
    
    ZTT->Draw("sameh");
    //ZLL->Draw("sameh");
    TT->Draw("sameh");
    ZLL->Draw("sameh");
    //W->Draw("sameh");
    //  VV->Draw("sameh");
    QCD->Draw("sameh");
    
    canv1->Update();
    
   // InitSignal(vbfH,2);
    InitSignal(SMH,2);
    
    if (showSignal)
    {
        legend->AddEntry(SMH,"SM Higgs(125) #times 10","f");
    }
    
    if (showSignal)
    {
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
    }

    
    bkgdErr->SetMarkerSize(0);
    int new_idx = CreateTransparentColor(13,0.4);
    bkgdErr->SetFillColor(new_idx);
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
    
    FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
    DrawTitle(pads[0], "12.9 fb^{-1} (13 TeV)", 3);
    DrawTitle(pads[0], "e#mu", 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();
    
    canv1->Print("/nfs/dust/cms/user/ywen/Taus2EMu_80X/workSpacev2/figures/"+Variable+Suffix+suffix+".pdf");
}
