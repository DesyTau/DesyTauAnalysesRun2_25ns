#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotDataCard_em(TString filename = "htt_em.inputs-sm-13TeV-mvis_noPU.root",
		  TString xtitle = "m_{vis} [GeV]",
		  TString ytitle = "Events",
		  float xLower = 0,
		  float xUpper = 150,
		  int numberofbins = 30,
		  bool logY = false,
		  bool legLeft = false) 

{
	SetStyle();
	TFile * file = new TFile(filename);
	TString histName[7] = {"data_obs",
	"ZTT",
	"ZLL", 
	"TT",
	"W",
	"VV",
	"QCD"
	};
	TH1D * histDataOld = (TH1D*)file->Get("em_inclusive/data_obs");
	int nBins = histDataOld->GetNbinsX();
	float xMin = histDataOld->GetBinLowEdge(1);
	float xMax = histDataOld->GetBinLowEdge(nBins+1);
	float bins[100];
	int nBinsNew = numberofbins;
        float binWidth = (xMax-xMin)/float(nBinsNew);
	for (int iB=0; iB<=nBinsNew; ++iB)
		bins[iB] = xMin + float(iB)*binWidth;
	TH1D * histData = new TH1D("histData","",nBinsNew,bins);

	TH1D * qcdHist = new TH1D("qcdHist","",nBinsNew,bins);
	TH1D * vvHist = new TH1D("vvHist","",nBinsNew,bins);
	TH1D * wHist = new TH1D("wHist","",nBinsNew,bins);
	TH1D * ttHist  = new TH1D("ttHist","",nBinsNew,bins);
        TH1D * zllHist = new TH1D("zllHist","",nBinsNew,bins);
	TH1D * ztautauHist = new TH1D("ztautauHist","",nBinsNew,bins);

	int nSamples = 7;

  //  return;

	for (int iS=0; iS<nSamples; ++iS) 
	{
		//std::cout << "Sample = " << iS << std::endl;
		TH1D * histOld = (TH1D*)file->Get("em_inclusive/"+histName[iS]);
		TH1D * hist = TH1DtoTH1D(histOld,nBinsNew,bins,true,"_new_"+histName[iS]);
    		if(iS==0) histData = hist;
		if(iS==1) qcdHist = hist;
		if(iS==2) vvHist = hist;
		if(iS==3) wHist = hist;
		if(iS==4) ttHist = hist;
		if(iS==5) zllHist = hist;
		if(iS==6) ztautauHist = hist;  
	}

	//  float dataEvents = 0;
	//  float ttEvents = 0;

	float totData = histData->GetSumOfWeights();
	float totZTauTau = ztautauHist->GetSumOfWeights();
	float totZLL = zllHist->GetSumOfWeights();
	float totTT = ttHist->GetSumOfWeights();
	float totW = wHist->GetSumOfWeights();
	float totVV = vvHist->GetSumOfWeights();
	float totQCD = qcdHist->GetSumOfWeights();

	float sfTT = (totData-totZTauTau-totZLL-totW-totVV-totQCD)/totTT;
	float sfTTerr = TMath::Sqrt(totData)/totTT;
	std::cout << "TTJets SF = " << sfTT << " +/- " << sfTTerr << std::endl; 


	vvHist->Add(vvHist,qcdHist);
	wHist->Add(wHist,vvHist);
	ttHist->Add(ttHist,wHist);
	zllHist->Add(zllHist,ttHist);
	ztautauHist->Add(ztautauHist,zllHist);

	InitData(histData);
	InitHist(qcdHist,"","",kMagenta,1001);
	InitHist(vvHist,"","",kRed-5,1001);
	InitHist(wHist,"","",kRed,1001);
	InitHist(ttHist,"","",kBlue-6,1001);
	InitHist(zllHist,"","",kGreen+2,1001);
	InitHist(ztautauHist,"","",kOrange,1001);
	histData->GetXaxis()->SetTitle(xtitle);
	histData->GetYaxis()->SetTitle(ytitle);
	histData->GetYaxis()->SetTitleOffset(1.3);
	histData->GetYaxis()->SetTitleSize(0.06);
	histData->GetXaxis()->SetRangeUser(xLower,xUpper);
	float yUpper = histData->GetMaximum();
	if (logY)
	histData->GetYaxis()->SetRangeUser(0.5,2*yUpper);
	else
	histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);
	histData->SetMarkerSize(1.5);
	histData->GetXaxis()->SetLabelSize(0);
	histData->GetYaxis()->SetLabelSize(0.06);

	//  nData = histData->GetSum();
	//  float nMC   = ttHist->GetSum();
	//  float eData = TMath::Sqrt(nData);

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

	histData->Draw("e1");
	ztautauHist->Draw("sameh");
	zllHist->Draw("sameh");
	ttHist->Draw("sameh");
	wHist->Draw("sameh");
	vvHist->Draw("sameh");
	qcdHist->Draw("sameh");
	histData->Draw("e1same");

	float chi2 = 0;
	for (int iB=1; iB<=nBinsNew; ++iB) 
	{
		float xData = histData->GetBinContent(iB);
		float xMC = ztautauHist->GetBinContent(iB);
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
	leg->AddEntry(histData,"Data","lp");
	leg->AddEntry(qcdHist,"QCD","f");
	leg->AddEntry(vvHist,"VV","f");
	leg->AddEntry(wHist,"W","f");
	leg->AddEntry(ttHist,"t#bar{t}","f");
        leg->AddEntry(zllHist,"Z#rightarrow ll","f");
	leg->AddEntry(ztautauHist,"Z#rightarrow#tau#tau","f");
	leg->Draw();

	TLatex * cms = new TLatex(0.25,0.94,"CMS Preliminary L = 1280 pb^{-1} at #sqrt{s} = 13 TeV");

	cms->SetNDC();
	cms->SetTextSize(0.05);
	cms->Draw();

	if (logY) upper->SetLogy(true);
	upper->Draw("SAME");
	upper->RedrawAxis();
	upper->Modified();
	upper->Update();
	canv1->cd();

	TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
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

	for (int iB=1; iB<=nBinsNew; ++iB) 
	{
		float x1 = histData->GetBinContent(iB);
		float x2 = ztautauHist->GetBinContent(iB);
		if (x1>0&&x2>0) 
		{
			float e1 = histData->GetBinError(iB);
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


  // ------------>Primitives in pad: lower
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

	canv1->Print("Ztautau_em.png");
	canv1->Print("Ztautau_em.pdf","Portrait pdf");


}
