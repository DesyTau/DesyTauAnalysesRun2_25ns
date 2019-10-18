//----------------------Version 2.0-------------------------//
//Plot variables for E -> Tau Fake Rates study
//Author: Yiwen Wen & Andrea Cardini
//DESY
//----------------------------------------------------------//
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
#include "DesyTauAnalyses/NTupleMaker/test/Plotting.h"
#include "DesyTauAnalyses/NTupleMaker/test/Plotting_Style.h"
#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
void PlotNtupleVariables_ETauFRwithDeepTau(
                                 TString Variable= "m_vis",
                                 TString xtitle = "visible mass [GeV] ",
                                 int nbins = 30,
                                 float xmin = 0,
                                 float xmax = 250,
				 TString wp = "Medium",
				 TString wpIso = "Tight",
				 TString Year="2018",
				 bool DeepTau = true,
				 TString MCweight="effweight*mcweight*",
				 TString Cut="(iso_1<0.1&&mt_1<30)", 
                                 TString ytitle = "Events",
                                 bool passProbe = true,
                                 bool logY = false,
                                 bool legLeft = false)
{  
  const int nSamples = 20;

  bool applyPU = true;
  TString tauIso = "tauby"+wpIso+"IsolationMVArun2v1DBoldDMwLT";
  TString againstMu = "tauagainstMuonLoose3";
  TString DeepTauString ="";

if(DeepTau>0.5){
  DeepTauString ="DeepTau";
  tauIso = "tauby"+wpIso+"DeepTau2017v2VSjet";
  againstMu = "tauagainstMuLooseDeepTau";
  }
  

  TString directory="./";
  TH1::SetDefaultSumw2();
  SetStyle();
    TString samples[nSamples] =
    {
        "EGamma_Run2018",//(0)data
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",// (1)Drell-Yan Z->EE
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (2)Drell-Yan ZJ
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (3)Drell-Yan ZTT(tau -> lepton)
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", // (4)Drell-Yan ZTT(hadronic tau)
        "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",//(5)TTbar leptonic
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",//(6) semileptonic
        "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",//(7) hadronic
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",// (8)WJets
        "WW_TuneCP5_13TeV-pythia8",// (9)WW
        "WZ_TuneCP5_13TeV-pythia8",// (10)WZ
        "ZZ_TuneCP5_13TeV-pythia8",// (11)ZZ
        "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (12) SingleTop tW tbar
        "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (13) SingleTop tW t
        "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",// (14) SingleTop t antitop
        "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"// (15) SingleTop t top
    };
    double WJetsweight=0.8369;

    //double WJetsweight=1.;
    double DYJetsweight=1.;

    double xsec[nSamples] = {
        1, // (0)data
        DYJetsweight,// (1)Drell-Yan Z->EE
        DYJetsweight,// (2)Drell-Yan ZJ
        DYJetsweight,// (3)Drell-Yan ZTT(tau->lepton)
        DYJetsweight,// (4)Drell-Yan ZTT(hadronic tau)
        88.29, // (5)TTPowHeg leptonic
        377.96, // (6)TTBar semileptonic
        365.35, // (7)TTBar hadronic
        WJetsweight, // (6)WJetsToLNu_MG
        WJetsweight, // (7)WJetsToLNu_MG
        WJetsweight, // (8)WJetsToLNu_MG
        WJetsweight, // (9)WJetsToLNu_MG
        WJetsweight, // (10)WJetsToLNu_MG
        63.21,  // WW (11)
        22.82,   // WZ (12)
        10.32,  // ZZ (13)
        38.06,   // ST_tW_antitop_5f_inclusiveDecays (14)
        38.09,   // ST_tW_top_5f_inclusiveDecays (15)
        80.95, // (16) SingleTop t antitop
        136.02// (17) SingleTop t top
    };
	float lumi;
	if(Year=="2018") lumi = 59970;
	else if(Year=="2017") lumi = 41860;
	else if(Year=="2016") lumi = 36773;
	else exit(EXIT_FAILURE);

	//Deal with bins and binning
	float xMin = xmin;
	float xMax = xmax;
	int nBins = nbins;
	float bins[100];
	float binWidth = (xMax-xMin)/float(nBins);
        for (int iB=0; iB<=nBins; ++iB)
        	bins[iB] = xMin + float(iB)*binWidth;

	TH1D * hist[nSamples];
	TH1D * histSS[nSamples];

	//inistiating cuts
	TString cuts[nSamples];
	TString cutsSS[nSamples];
	
	TString qcdweight ="1.06*";
    	
	Cut+="*("+againstMu+">0.5&&"+tauIso+">0.5)";

	if (DeepTau){
		Cut="(DecayModeProbe!=5&&DecayModeProbe!=6)*"+Cut;
			
	}else{
		Cut="(tau_decayModeFinding>0.5)*"+Cut;
		
	}	
				

	for (int i=0; i<nSamples; ++i)
	{
	  cuts[i] = MCweight+Cut+"*(os>0.5) ";
        cutsSS[i] = qcdweight+MCweight+Cut+"*(os<0.5)";
  	}
	TString ZEEcut    = MCweight+Cut+"* (gen_match_1 == 1 && gen_match_2 == 1) ";
	TString ZJcut     = MCweight+Cut+"* (gen_match_2 == 6) ";
	TString ZTT_elcut = MCweight+Cut+"* (gen_match_1 == 3 && (gen_match_2 == 3 || gen_match_2 == 4)) ";
	TString ZTT_etcut = MCweight+Cut+"* (gen_match_1 == 3 && gen_match_2 == 5)"; 
	
	//1= prompt ele at gen level 2=prompt mu 3= non prompt ele (from tau) 4= non promt mu (from tau) 5= tau at gen level 6= jet

	if (DeepTau){
		ZTT_etcut += "*0.83";
	}else{
		ZTT_etcut += "*0.9";
	}

  	cuts[0] = Cut       +"*(os>0.5)";
	cuts[1] = ZEEcut    +"*(os>0.5)";
	cuts[2] = ZJcut     +"*(os>0.5)";
	cuts[3] = ZTT_elcut +"*(os>0.5)";
	cuts[4] = ZTT_etcut +"* (os>0.5)";

    
	cutsSS[0] = qcdweight+Cut       +"*(os<0.5)";
	cutsSS[1] = qcdweight+ZEEcut    +"*(os<0.5)";
	cutsSS[2] = qcdweight+ZJcut     +"*(os<0.5)";
	cutsSS[3] = qcdweight+ZTT_elcut +"*(os<0.5)";
	cutsSS[4] = qcdweight+ZTT_etcut +"*(os<0.5)";
	
									//actual pass-fail selection
	if(passProbe)
	  {
	    for (int i=0; i<nSamples; ++i)
	      {
		cuts[i] ="(tauagainstEle"+wp+DeepTauString+">0.5)*"+cuts[i];
		cutsSS[i] ="(tauagainstEle"+wp+DeepTauString+">0.5)*"+cutsSS[i];
		}
	    
	  }
	else
	  {
	    for (int i=0; i<nSamples; ++i)
	      {
		cuts[i] ="(tauagainstEle"+wp+DeepTauString+"<0.5)*"+cuts[i];
		cutsSS[i] ="(tauagainstEle"+wp+DeepTauString+"<0.5)*"+cutsSS[i];
		}
	  }
	
	if(applyPU)
	  {
	    for (int i=1; i<nSamples; ++i)
	      {
		cuts[i] ="puweight*"+cuts[i];
		cutsSS[i] ="puweight*"+cutsSS[i];
	      }
	  }
	
	for (int i=0; i<nSamples; ++i){
	  	cout << directory+samples[i]+".root" << endl;
	        TFile * file = new TFile(directory+samples[i]+".root");
		TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
		TTree * tree = (TTree*)file->Get("ETauFR");
		
		double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
		if(samples[i].Contains("Jets")){
		  norm=xsec[i]*lumi;
		  cuts[i]+="*stitchweight";
		  cutsSS[i]+="*stitchweight";
		}
		cout<< norm <<endl;
		TString histName = samples[i] + "_"+Variable;
		TString histNameSS = samples[i] + "_"+Variable+"_ss";
		hist[i] = new TH1D(histName,"",nBins,xMin,xMax);
		hist[i]->Sumw2();
		histSS[i] = new TH1D(histNameSS,"",nBins,xMin,xMax);
		histSS[i]->Sumw2();
		tree->Draw(Variable+">>"+histName,cuts[i]);
		tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
		if (i>0) // if sample is MC sample -> Scale to xsec and luminosity
		  {
		    hist[i]   -> Scale(norm);
		    histSS[i] -> Scale(norm);
		  }
		cout << samples[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
		
	}

    
	

  // ********************************************
  // ***** Stitching of Drell-Yan samples *******
  // ********************************************
    
	TH1D * histZEE[4];
	TH1D * histZJ[4];
	TH1D * histZTT_EL[4];
	TH1D * histZTT_ET[4];
	TH1D * histSSZEE[4];
	TH1D * histSSZJ[4];
	TH1D * histSSZTT_EL[4];
	TH1D * histSSZTT_ET[4];
	TString refSamplesDY=samples[1];

	
	
	// filling histograms for DY samples
	for (int i=0; i<4; i++) { // run over samples
	  TString name=refSamplesDY;
	  TString njets=TString::Itoa(i+1,10);
	  njets+="Jets";
	  TString filename=name.ReplaceAll("Jets",njets);
	  cout<<filename<<endl;
	
	  TFile * file = new TFile(directory+filename+".root");
	  TTree * tree = (TTree*)file->Get("ETauFR");
	  double norm = lumi;
	  
	  TString histNameZEE   = filename + Variable + "_zee_os";
	  TString histNameSSZEE = filename + Variable + "_zee_ss";
	  TString histNameZJ   = filename + Variable + "_zj_os";
	  TString histNameSSZJ = filename + Variable + "_zj_ss";
	  TString histNameZTT_EL   = filename + Variable + "_ztt_el_os";
	  TString histNameSSZTT_EL = filename + Variable + "_ztt_el_ss";
	  TString histNameZTT_ET   = filename + Variable + "_ztt_et_os";
	  TString histNameSSZTT_ET = filename + Variable + "_ztt_et_ss";
	  histZEE[i]   = new TH1D(histNameZEE,"",nBins,xmin,xmax);
	  histSSZEE[i] = new TH1D(histNameSSZEE,"",nBins,xmin,xmax);
	  histZJ[i]   = new TH1D(histNameZJ,"",nBins,xmin,xmax);
	  histSSZJ[i] = new TH1D(histNameSSZJ,"",nBins,xmin,xmax);
	  histZTT_EL[i]   = new TH1D(histNameZTT_EL,"",nBins,xmin,xmax);
	  histSSZTT_EL[i] = new TH1D(histNameSSZTT_EL,"",nBins,xmin,xmax);
	  histZTT_ET[i]   = new TH1D(histNameZTT_ET,"",nBins,xmin,xmax);
	  histSSZTT_ET[i] = new TH1D(histNameSSZTT_ET,"",nBins,xmin,xmax);
	  

	  tree -> Draw(Variable+">>"+histNameZEE,     cuts[1]);
	  tree -> Draw(Variable+">>"+histNameSSZEE,   cutsSS[1]);  
	  tree -> Draw(Variable+">>"+histNameZJ,      cuts[2]);
	  tree -> Draw(Variable+">>"+histNameSSZJ,    cutsSS[2]);
	  tree -> Draw(Variable+">>"+histNameZTT_EL,  cuts[3]);
	  tree -> Draw(Variable+">>"+histNameSSZTT_EL,cutsSS[3]);
	  tree -> Draw(Variable+">>"+histNameZTT_ET,  cuts[4]);
	  tree -> Draw(Variable+">>"+histNameSSZTT_ET,cutsSS[4]);
	
	  histZEE[i]   ->Scale(norm); 
	  histSSZEE[i] ->Scale(norm); 
	  histZJ[i]   ->Scale(norm); 
	  histSSZJ[i] ->Scale(norm); 
	  histZTT_EL[i]   ->Scale(norm);
	  histSSZTT_EL[i] ->Scale(norm);
	  histZTT_ET[i]   ->Scale(norm);
	  histSSZTT_ET[i] ->Scale(norm);
	  
  }
	
	for (int iDY=0; iDY<4; ++iDY) {
	  hist[1]  -> Add(hist[1],histZEE[iDY]);
	  hist[2]  -> Add(hist[2],histZJ[iDY]);
	  hist[3]  -> Add(hist[3],histZTT_EL[iDY]);
	  hist[4]  -> Add(hist[4],histZTT_ET[iDY]);
	  histSS[1]-> Add(histSS[1],histSSZEE[iDY]);
	  histSS[2]-> Add(histSS[2],histSSZJ[iDY]);
	  histSS[3]  -> Add(histSS[3],histSSZTT_EL[iDY]);
	  histSS[4]  -> Add(histSS[4],histSSZTT_ET[iDY]);
	}
	
    //adding up ttbar backgrounds
    for (int iH=6; iH<8; ++iH)
    {
        hist[5]->Add(hist[5],hist[iH]);
    }
	
	for (int iH=9; iH<13; ++iH)
	  {
	    hist[8]->Add(hist[8],hist[iH]);
	  }

	//  adding up VV backgrounds
	for (int iH=14; iH<nSamples; ++iH)
	{
        hist[13]->Add(hist[13],hist[iH]);
    }
    for (int iH=1; iH<nSamples; ++iH)
    {
        histSS[0]->Add(histSS[0],histSS[iH],1,-1);
    }
    
    //Getting rid of negative bin in QCD
    for (int iB=1; iB<=nBins; ++iB)
    {
        float ySS = histSS[0]->GetBinContent(iB);
        if (ySS<0)
        {
            histSS[0]->SetBinContent(iB,0.);
            histSS[0]->SetBinError(iB,0.);
        }
        
    }
    
        TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
	TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
	TH1D * ZEE = (TH1D*)hist[1]->Clone("ZEE");
	TH1D * ZJ = (TH1D*)hist[2]->Clone("ZJ");
	TH1D * ZTT_el = (TH1D*)hist[3]->Clone("ZTT_el");
	TH1D * ZTT_et = (TH1D*)hist[4]->Clone("ZTT_et");
	TH1D * W = (TH1D*)hist[8]->Clone("W");
	TH1D * TT = (TH1D*)hist[5]->Clone("TT");
	TH1D * VV = (TH1D*)hist[13]->Clone("VV");
    
    float totData = data_obs->GetSumOfWeights();
    float totZEE = ZEE->GetSumOfWeights();
    float totZJ = ZJ->GetSumOfWeights();
    float totZTT_el = ZTT_el->GetSumOfWeights();
    float totZTT_et = ZTT_et->GetSumOfWeights();
    float totQCD = QCD->GetSumOfWeights();
    float totTT = TT->GetSumOfWeights();
    float totVV = VV->GetSumOfWeights();
    float totW = W->GetSumOfWeights();
    float totMC = totZEE+totZJ+totZTT_el+totZTT_et+totTT+totW+totVV+totQCD;

    float sfW = (totData-totZEE-totZJ-totZTT_el-totZTT_et-totTT-totQCD-totVV)/totW;
    float sfWerr = TMath::Sqrt(totData)/totW;
    std::cout << "Data: " << data_obs->GetSumOfWeights() << std::endl;
    std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
    std::cout << "VV  : " << VV->GetSumOfWeights() << std::endl;
    std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
    std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
    std::cout << "ZEE : " << ZEE->GetSumOfWeights() << std::endl;
    std::cout << "ZJ : " << ZJ->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_el : " << ZTT_el->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_et : " << ZTT_et->GetSumOfWeights() << std::endl;
    std::cout << "MC : " << totMC << std::endl;

    std::cout << "sfW = " << sfW << " +/- " << sfWerr << std::endl;

	TH1D * dummy = (TH1D*)ZEE->Clone("dummy");
    float errLumi = 0.03;
    float errQCD = 0.1;
    float errDY=0.1;
    float errVV = 0.15;
    float errW = 0.15;
    float errTT = 0.15;
    
    for (int iB=1; iB<=nBins; ++iB)
    {
        float eQCD = errQCD*QCD->GetBinContent(iB);
        float eVV = errVV*VV->GetBinContent(iB);
        float eDYMM = errDY*ZEE->GetBinContent(iB);
        float eDYTT_el = errDY*ZTT_el->GetBinContent(iB);
        float eDYTT_et = errDY*ZTT_et->GetBinContent(iB);
        float eDYJ = errDY*ZJ->GetBinContent(iB);
        float eW = errW*W->GetBinContent(iB);
        float eTT = errTT*TT->GetBinContent(iB);
        float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eTT*eTT+eDYMM*eDYMM+eDYTT_el*eDYTT_el+eDYTT_et*eDYTT_et+eDYJ*eDYJ;
        float errTot = TMath::Sqrt(err2);
        dummy->SetBinError(iB,errTot);
    }
    
	W->Add(W,QCD);
	TT->Add(TT,W);
	VV->Add(VV,TT);
	ZTT_et->Add(ZTT_et,VV);
	ZTT_el->Add(ZTT_el,ZTT_et);
	ZJ->Add(ZJ,ZTT_el);
	ZEE->Add(ZEE,ZJ);
	
	TH1D * bkgdErr = (TH1D*)ZEE->Clone("bkgdErr");
  	bkgdErr->SetFillStyle(3013);
  	bkgdErr->SetFillColor(1);
  	bkgdErr->SetMarkerStyle(21);
  	bkgdErr->SetMarkerSize(0);

  	for (int iB=1; iB<=nBins; ++iB)
    {  
        W->SetBinError(iB,0);
        QCD->SetBinError(iB,0);
        VV->SetBinError(iB,0);
        TT->SetBinError(iB,0);
        ZTT_el->SetBinError(iB,0);
        ZTT_et->SetBinError(iB,0);
        ZJ->SetBinError(iB,0);
        ZEE->SetBinError(iB,0);
        float eStat =  bkgdErr->GetBinError(iB);
        float X = bkgdErr->GetBinContent(iB);
        float eLumi = errLumi * X;
        float eBkg = dummy->GetBinError(iB);
        float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg);
        bkgdErr->SetBinError(iB,Err);
    }  
   	//Colors 
 	Int_t colorZEE = TColor::GetColor("#ffcc66");
	Int_t colorZOther = TColor::GetColor("#4496C8");
	Int_t colorTT = TColor::GetColor("#9999CC");
	Int_t colorVV = TColor::GetColor("#6F2D35");
	Int_t colorW = TColor::GetColor("#DE5A6A");
	Int_t colorQCD = TColor::GetColor("#FFCCFF");
	InitHist(QCD,TColor::GetColor("#FFCCFF"));
	InitHist(ZEE,TColor::GetColor("#DE5A6A"));
	InitHist(TT,TColor::GetColor("#9999CC"));
	InitHist(VV,TColor::GetColor("#6F2D35"));
	InitHist(ZJ,TColor::GetColor("#FFCC66"));
	InitHist(W,TColor::GetColor("#4496C8"));

	InitData(data_obs);
	data_obs->GetXaxis()->SetTitle(xtitle);
 	data_obs->GetYaxis()->SetTitle(ytitle);
	data_obs->GetYaxis()->SetTitleOffset(1.5);
	data_obs->GetYaxis()->SetTitleSize(0.06);
	data_obs->GetXaxis()->SetRangeUser(xmin,xmax);
	float yUpper = data_obs->GetMaximum();
	if (logY)
	data_obs->GetYaxis()->SetRangeUser(0.5,2*yUpper);
	else
	data_obs->GetYaxis()->SetRangeUser(0,1.2*yUpper);
	data_obs->SetMarkerSize(1.5);
	data_obs->GetXaxis()->SetLabelSize(0);
	data_obs->GetYaxis()->SetLabelSize(0.06);
	
	TCanvas * canv1 = MakeCanvas("canv1", "", 1000, 800);
	
	TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
	upper->Draw();
	upper->cd();
	upper->SetFillColor(0);
	upper->SetBorderMode(0);
	upper->SetBorderSize(10);
	upper->SetTickx(1);
	upper->SetTicky(1);
	upper->SetLeftMargin(0.17);
	upper->SetRightMargin(0.05);
	upper->SetBottomMargin(0.02);
 	upper->SetFrameFillStyle(0);
	upper->SetFrameLineStyle(0);
  	upper->SetFrameLineWidth(2);
  	upper->SetFrameBorderMode(0);
  	upper->SetFrameBorderSize(10);
  	upper->SetFrameFillStyle(0);
  	upper->SetFrameLineStyle(0);
  	upper->SetFrameLineWidth(2);
  	upper->SetFrameBorderMode(0);
  	upper->SetFrameBorderSize(10);

	//Drawing histogram
  	data_obs->Draw("e1");
	ZEE->Draw("sameh");
  	ZJ->Draw("sameh");//ZJ is ZOther, including ZJ+ZTT
  	VV->Draw("sameh");
  	TT->Draw("sameh");
	W->Draw("sameh");
	QCD->Draw("sameh");
  	data_obs->Draw("e1same");
	bkgdErr->Draw("e2same");
	//Calculating chi2
	float chi2 = 0;
	for (int iB=1; iB<=nBins; ++iB) 
	{
		float xData = data_obs->GetBinContent(iB);
		float xMC = ZEE->GetBinContent(iB);
		if (xMC>1e-1) 
		{
      			float diff2 = (xData-xMC)*(xData-xMC);
      			chi2 += diff2/xMC;
    		}
  	}
  	std::cout << std::endl;
  	std::cout << "Chi2 = " << chi2 << std::endl;
  	std::cout << std::endl;

  	float x1Leg = 0.70;
  	float x2Leg = 0.90;
  	if (legLeft) 
	{
    		x1Leg = 0.20;
    		x2Leg = 0.45;
  	}
	TLegend * leg = new TLegend(x1Leg,0.6,x2Leg,0.88);
  	SetLegendStyle(leg);
  	leg->SetTextSize(0.05);
  	leg->AddEntry(data_obs,"Data","lp");
  	leg->AddEntry(VV,"Dibosons","f");
  	leg->AddEntry(W,"WJets","f");
 	leg->AddEntry(QCD,"QCD","f");
  	leg->AddEntry(TT,"t#bar{t}","f");
	leg->AddEntry(ZJ,"DY others","f");
  	leg->AddEntry(ZEE,"Z#rightarrow e e","f");
  	leg->Draw();
  	//plotchannel("e#mu");
    //	if (!applyPU) suffix = "_noPU";
	TString title;
	if (Year== "2016") DrawTitle(pads[0], "35.9 fb^{-1} (13 TeV, 2016)", 3);
	if (Year== "2017") DrawTitle(pads[0], "41.9 fb^{-1} (13 TeV, 2017)", 3);
	if (Year== "2018") DrawTitle(pads[0], "59.9 fb^{-1} (13 TeV, 2018)", 3);
  	//TLatex * cms = new TLatex(0.20,0.94,"CMS 2018, L = 59.9 fb^{-1} at #sqrt{s} = 13 TeV");

  	//cms->SetNDC();
  	//cms->SetTextSize(0.05);
  	//cms->Draw();
	DrawTitle(pads[0], "#scale[1.2]{         #bf{CMS} Work in progress}", 1);
	//TLatex * workinprogress = new TLatex(x1Leg,0.55,"Work in progress");
	//workinprogress->SetNDC();
	//workinprogress->SetTextSize(0.05);
	//workinprogress->Draw();
    
    
  	if (logY) upper->SetLogy(true);

 	upper->Draw("SAME");
  	upper->RedrawAxis();
  	upper->Modified();
  	upper->Update();
 	canv1->cd();
	
	TH1D * ratioH = (TH1D*)data_obs->Clone("ratioH");
	TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  	ratioH->SetMarkerColor(1);
 	ratioH->SetMarkerStyle(20);
  	ratioH->SetMarkerSize(1.5);
  	ratioH->SetLineColor(1);
  	ratioH->GetYaxis()->SetRangeUser(0.3,1.8);
 	ratioH->GetYaxis()->SetNdivisions(505);
  	ratioH->GetXaxis()->SetLabelFont(42);
  	ratioH->GetXaxis()->SetLabelOffset(0.04);
  	ratioH->GetXaxis()->SetLabelSize(0.14);
  	ratioH->GetXaxis()->SetTitleSize(0.13);
  	ratioH->GetXaxis()->SetTitleOffset(1.2);
  	ratioH->GetYaxis()->SetTitle("obs/exp");
 	ratioH->GetYaxis()->SetLabelFont(42);
  	ratioH->GetYaxis()->SetLabelOffset(0.015);
  	ratioH->GetYaxis()->SetLabelSize(0.13);
  	ratioH->GetYaxis()->SetTitleSize(0.14);
  	ratioH->GetYaxis()->SetTitleOffset(0.5);
  	ratioH->GetXaxis()->SetTickLength(0.07);
  	ratioH->GetYaxis()->SetTickLength(0.04);
  	ratioH->GetYaxis()->SetLabelOffset(0.01);
	
	for (int iB=1; iB<=nBins; ++iB) 
	{
    		float x1 = data_obs->GetBinContent(iB);
    		float x2 = ZEE->GetBinContent(iB);
    		ratioErrH->SetBinContent(iB,1.0);
    		ratioErrH->SetBinError(iB,0.0);
		float xBkg = bkgdErr->GetBinContent(iB);
   		float errBkg = bkgdErr->GetBinError(iB);
		if (xBkg>0) 
		{
     			float relErr = errBkg/xBkg;
      			ratioErrH->SetBinError(iB,relErr);
    		}
		if (x1>0&&x2>0) 
		{
      			float e1 = data_obs->GetBinError(iB);
      			float ratio = x1/x2;
      			float eratio = e1/x2;
      			ratioH->SetBinContent(iB,ratio);
     			ratioH->SetBinError(iB,eratio);
    		}
   		else 
		{
      			ratioH->SetBinContent(iB,1000);
    		}
  	}


	TPad *lower = new TPad("lower", "pad",0,0,1,0.30);
  	lower->Draw();
  	lower->cd();
  	lower->SetFillColor(0);
  	lower->SetBorderMode(0);
  	lower->SetBorderSize(10);
 	lower->SetGridy();
  	lower->SetTickx(1);
  	lower->SetTicky(1);
  	lower->SetLeftMargin(0.17);
 	lower->SetRightMargin(0.05);
  	lower->SetTopMargin(0.026);
  	lower->SetBottomMargin(0.35);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);

  	ratioH->Draw("e1");
	ratioErrH->Draw("e2same");
  	lower->Modified();
  	lower->RedrawAxis();
  	canv1->cd();
  	canv1->Modified();
  	canv1->cd();
  	canv1->SetSelected(canv1);
	canv1->Print("./Plots/"+Variable+"_"+wp+DeepTauString+"Iso"+wpIso+".png");
	canv1->Print("./Plots/"+Variable+"_"+wp+DeepTauString+"Iso"+wpIso+".pdf","Portrait pdf");


}

