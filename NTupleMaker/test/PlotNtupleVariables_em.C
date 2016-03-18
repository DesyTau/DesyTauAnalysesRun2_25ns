//----------------------Version 0.01-------------------------//
//Plotting Macro for Z ->tautau-> e+mu study 
//Author: Yiwen Wen
//DESY
//-----------------------------------------------------------//
#include "HttStylesNew.cc"
#include "HtoH.h"
void PlotNtupleVariables(TString varName= "m_vis",
			 TString xtitle = "m_{vis} [GeV]",
			 TString ytitle = "Events",
			 float xLower = 0,
			 float xUpper = 400,
			 int numberofbins = 40,
			 bool logY = true,
			 bool legLeft = false)
{
	SetStyle();
	TString samples[17] = {"SynNTuples_MuonEG",//(0) data
			       "SynNTuples_DYJetsToLL_M-50_MG", // (1) Drell-Yan 
                               "SynNTuples_DYJetsToLL_M-50_MG", // (2) Drell-Yan
			       "SynNTuples_VVTo2L2Nu", // (3) 
                               "SynNTuples_ZZTo2L2Q", // (4)
                               "SynNTuples_ZZTo4L", // (5)
                               "SynNTuples_WWTo1L1Nu2Q", // (6) 
                               "SynNTuples_WZTo2L2Q", // (7)
                               "SynNTuples_WZJets", // (8)
                               "SynNTuples_WZTo1L3Nu", // (9)
                               "SynNTuples_WZTo1L1Nu2Q", // (10)
                               "SynNTuples_ST_tW_antitop_5f_inclusiveDecays", // (11)
                               "SynNTuples_ST_tW_top_5f_inclusiveDecays", // (12)
                               "SynNTuples_ST_t-channel_top_4f_leptonDecays", // (13)
                               "SynNTuples_ST_t-channel_antitop_4f_leptonDecays", // (14)
                               "SynNTuples_WJetsToLNu_MG",//(15)
                               "SynNTuples_TTPowHeg"// (16)
                              };


	float xsec[17] = {1, //data(0)
			  6025,   // DYJetsToLL_M-50_MG (1)
                          6025,   // DYJetsToLL_M-50_MG (2)(not Z to tautau)
			  11.95,  // VVTo2L2Nu (3)
                          3.22,   // ZZTo2L2Q (4)
                          1.212,    // ZZTo4L (5)
                          49.997,   // WWTo1L1Nu2Q (6)
                          5.595,    // WZTo2L2Q (7) 
                          5.26,    // WZJets (8)
                          3.05,   // WZTo1L3Nu (9)
                          10.71,   // WZTo1L1Nu2Q (10)
                          35.6,   // ST_tW_antitop_5f_inclusiveDecays (11)
                          35.6,   // ST_tW_top_5f_inclusiveDecays (12)
                          136.95,   // ST_t-channel_top_4f_leptonDecays (13)
                          80.95,    // ST_t-channel_antitop_4f_leptonDecays (14)
                          61526,   // WJetsToLNu_MG (15)
                          831.8 // TTPowHeg (16)	
			  };

	float lumi = 2090;

	//Deal with bins and binning
	float xMin = xLower;
	float xMax = xUpper;
	int nBins = numberofbins;
	float bins[100];
	float binWidth = (xMax-xMin)/float(nBins);
        for (int iB=0; iB<=nBins; ++iB)
        	bins[iB] = xMin + float(iB)*binWidth;

	int nSamples = 17;
	TH1D * hist[17];
	TH1D * histSS[17];
	//inistiating cuts
	TString cuts[17];
	TString cutsSS[17];

	for (int i=0; i<17; ++i) 
	{
    		cuts[i] = "mcweight*(os>0.5)";
    		cutsSS[i] = "mcweight*(os<0.5)";
  	}

  	cuts[0] = "os>0.5";
	cuts[1] = "mcweight*(os>0.5&&isZTT)";
	cuts[2] = "mcweight*(os>0.5&&!isZTT)";

  	cutsSS[0] = "os<0.5";
	cutsSS[1] = "mcweight*(os<0.5&&isZTT)";
  	cutsSS[2] = "mcweight*(os<0.5&&!isZTT)";



std::cout <<"Here1"<<std::endl;
	for (int i=0; i<nSamples; ++i)
	{
		TFile * file = new TFile(samples[i]+".root");
		TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
		TTree * tree = (TTree*)file->Get("TauCheck");
		double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
		TString histName = samples[i] + "_"+varName;
		TString histNameSS = samples[i] + "_"+varName+"_ss";
		hist[i] = new TH1D(histName,"",nBins,xMin,xMax);
		hist[i]->Sumw2(); 
		histSS[i] = new TH1D(histNameSS,"",nBins,xMin,xMax);
		histSS[i]->Sumw2();

		tree->Draw("m_vis>>"+histName,cuts[i]);
		tree->Draw("m_vis>>"+histNameSS,cutsSS[i]);
		
		if(i>0)
		{	
			for (int iB=1; iB<=nBins; ++iB) 
			{
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
std::cout <<"Here2"<<std::endl;
	//  adding up single top and VV backgrounds
	for (int iH=4; iH<15; ++iH) 
	{
    		hist[3]->Add(hist[3],hist[iH]);
      	}	
	// subtracting background from SS
	for (int iH=1; iH<4; ++iH) 
	{
		histSS[0]->Add(histSS[0],histSS[iH],1,-1);
	}
	histSS[0]->Add(histSS[0],histSS[15],1,-1);
	histSS[0]->Add(histSS[0],histSS[16],1,-1);

	TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
	TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
	TH1D * ZTT = (TH1D*)hist[1]->Clone("ZTT");
	TH1D * ZLL = (TH1D*)hist[2]->Clone("ZLL");
	TH1D * W = (TH1D*)hist[15]->Clone("W");
	TH1D * TT = (TH1D*)hist[16]->Clone("TT");
	TH1D * VV = (TH1D*)hist[3]->Clone("VV");

	W->Add(W,QCD);
	TT->Add(TT,W);
	VV->Add(VV,TT);
	ZLL->Add(ZLL,VV);
	ZTT->Add(ZTT,ZLL);

	InitData(data_obs);
	InitHist(QCD,"","",kMagenta,1001);
	InitHist(W,"","",kGreen,1001);
	InitHist(TT,"","",kBlue-4,1001);
	InitHist(VV,"","",kRed,1001);
	InitHist(ZLL,"","",kCyan,1001);
	InitHist(ZTT,"","",kYellow,1001);	
	data_obs->GetXaxis()->SetTitle(xtitle);
 	data_obs->GetYaxis()->SetTitle(ytitle);
	data_obs->GetYaxis()->SetTitleOffset(1.3);
	data_obs->GetYaxis()->SetTitleSize(0.06);
	data_obs->GetXaxis()->SetRangeUser(xLower,xUpper);
	float yUpper = data_obs->GetMaximum();
	if (logY)
	data_obs->GetYaxis()->SetRangeUser(0.5,2*yUpper);
	else
	data_obs->GetYaxis()->SetRangeUser(0,1.2*yUpper);
	data_obs->SetMarkerSize(1.5);
	data_obs->GetXaxis()->SetLabelSize(0);
	data_obs->GetYaxis()->SetLabelSize(0.06);
	
	TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
	
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
	ZTT->Draw("sameh");
  	ZLL->Draw("sameh");
  	VV->Draw("sameh");
  	TT->Draw("sameh");
	W->Draw("sameh");
	QCD->Draw("sameh");
  	data_obs->Draw("e1same");

	//Calculating chi2
	float chi2 = 0;
	for (int iB=1; iB<=nBins; ++iB) 
	{
		float xData = data_obs->GetBinContent(iB);
		float xMC = ZTT->GetBinContent(iB);
		if (xMC>1e-1) 
		{
      			float diff2 = (xData-xMC)*(xData-xMC);
      			chi2 += diff2/xMC;
    		}
  	}
  	std::cout << std::endl;
  	std::cout << "Chi2 = " << chi2 << std::endl;
  	std::cout << std::endl;

  	float x1Leg = 0.65;
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
  	leg->AddEntry(VV,"SingleTop+VV","f");
  	leg->AddEntry(W,"WJets","f");
 	leg->AddEntry(QCD,"QCD","f");
  	leg->AddEntry(TT,"t#bar{t}","f");
	leg->AddEntry(ZLL,"Z#rightarrowll","f");
  	leg->AddEntry(ZTT,"Z#rightarrow#tau#tau#rightarrow e+#mu","f");
  	leg->Draw();

  	TLatex * cms = new TLatex(0.25,0.94,"CMS Preliminary L = 2.1 fb^{-1} at #sqrt{s} = 13 TeV");

  	cms->SetNDC();
  	cms->SetTextSize(0.05);
  	cms->Draw();

  	if (logY) upper->SetLogy(true);

 	upper->Draw("SAME");
  	upper->RedrawAxis();
  	upper->Modified();
  	upper->Update();
 	canv1->cd();
	
	TH1D * ratioH = (TH1D*)data_obs->Clone("ratioH");
  	ratioH->SetMarkerColor(1);
 	ratioH->SetMarkerStyle(20);
  	ratioH->SetMarkerSize(1.5);
  	ratioH->SetLineColor(1);
  	ratioH->GetYaxis()->SetRangeUser(0.5,1.5);
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
    		float x2 = ZTT->GetBinContent(iB);
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


	lower = new TPad("lower", "pad",0,0,1,0.30);
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

  	lower->Modified();
  	lower->RedrawAxis();
  	canv1->cd();
  	canv1->Modified();
  	canv1->cd();
  	canv1->SetSelected(canv1);
	canv1->Print(varName+".png");




}
