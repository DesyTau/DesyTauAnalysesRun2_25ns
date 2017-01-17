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
//WithoutNPV/
void PlotsForUnrolledDistr(TString directory = "/nfs/dust/cms/user/bobovnii/new/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/Processing/01Jets/MuTau/",
	  TString suffix = "",
	  TString File = "Templates3D_met_DZeta_MT2lester16_35invfb_mt_C1N2.root",
	  TString variable = "met_DZeta_MT2lester16",
	  TString Suffix = "",
	  TString xtitle = "",
	  TString ytitle = "Events",
      bool logY = true,
      bool BlindData = false

          )


{
	  vector<TString> vars_;vars_.clear();
	  vars_.push_back("met_DZeta_MT2lester16");
	  vars_.push_back("met_DZeta_MT2lester16");
  //ModTDRStyle();


int mycolor=TColor::GetColor("#ffcc66");
int mycolorvv=TColor::GetColor("#6F2D35");
//int mycolorvv=TColor::GetColor("#FF6633");
int mycolorqcd=TColor::GetColor("#ffccff");
int mycolortt=TColor::GetColor("#9999cc");
//int mycolorttx=TColor::GetColor("#bbccdd");
int mycolorttx=TColor::GetColor("#33CCFF");
int mycolorwjet=TColor::GetColor("#de5a6a");
//int mycolorwjet=TColor::GetColor("#66CC66");
int mycolordyj=TColor::GetColor("#ffcc66");

  bool plotLeg = true;
  int position = 2; // 0 - right, 1 - left, 2 - central
  int npos = 1;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  //  TString topweight("");

TString BinLabels[100] = {
"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31",
"32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
"61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93"
};


	ifstream BinLabelsFile;
	BinLabelsFile.open("LabelListTESR");
	int i=0;
	string name;
   	 while(getline(BinLabelsFile, name))
		{
//		BinLabels[i]=name;
		i++;
		}


	directory="";
    TFile * file = new TFile(File);
for (int svar=0;svar<vars_.size();++svar)
{
	TString Variable=vars_.at(svar);

    TH1D * TT = (TH1D*)file->Get("1D_tt_"+Variable);
    TH1D * WJ = (TH1D*)file->Get("1D_wj_"+Variable);
    TH1D * DY = (TH1D*)file->Get("1D_dyj_"+Variable);
    TH1D * ST = (TH1D*)file->Get("1D_sT_"+Variable);
    TH1D * VV = (TH1D*)file->Get("1D_dib_"+Variable);
    TH1D * TTX = (TH1D*)file->Get("1D_ttx_"+Variable);
    TH1D * QCD = (TH1D*)file->Get("1D_qcd_"+Variable);
    TH1D * allbkg = (TH1D*)file->Get("1D_tt_"+Variable);
   // allbkg->Reset();
    TH1D * histData = (TH1D*)file->Get("data_obs");
    TH1D * histSignal = (TH1D*)file->Get("1D_C1N2_100_LSP50_B_"+Variable);
    //TH1D * hist = (TH1D*)file->Get("");
int nBins  =   histData->GetNbinsX();
if (BlindData) {for (int iB=1; iB<=nBins; ++iB) {histData->SetBinContent(iB,0);}}


  std::cout << "ST   : " << ST->GetSumOfWeights() << " : "  <<" GetNbinsX= " << ST->GetNbinsX()<< std::endl;
  std::cout << "VV   : " << VV->GetSumOfWeights() << " : "  <<" GetNbinsX= " << VV->GetNbinsX() <<std::endl;
  std::cout << "QCD   : " << QCD->GetSumOfWeights() << " : " <<" GetNbinsX= " << QCD->GetNbinsX()<< std::endl;
  std::cout << "WJ   : " << WJ->GetSumOfWeights() << " : "  <<" GetNbinsX= " << WJ->GetNbinsX()<< std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : "  <<" GetNbinsX= " << TT->GetNbinsX()<< std::endl;
  std::cout << "TTX  : " << TTX->GetSumOfWeights() << " : "  <<" GetNbinsX= " << TTX->GetNbinsX()<< std::endl;
  std::cout << "DY : " << DY->GetSumOfWeights() << " : "  <<" GetNbinsX= " << DY->GetNbinsX()<< std::endl;
  std::cout << "histData : " << histData->GetSumOfWeights() << " : "  <<" GetNbinsX= " << histData->GetNbinsX()<< std::endl;
//  std::cout << "REST : " << REST->GetSumOfWeights() << " : "  <<" GetNbinsX= " << REST->GetNbinsX()<< std::endl;
  std::cout << "histSignal : " << histSignal->GetSumOfWeights() << " : "  <<" GetNbinsX= " << histSignal->GetNbinsX()<< std::endl;

/*
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

*/
  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************
 /* TH1D * dummy = (TH1D*)ZTT->Clone("dummy");
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
   */ 
  
//   allbkg->Add(TT);
   allbkg->Add(ST);
   allbkg->Add(TTX);
   allbkg->Add(WJ);
   allbkg->Add(DY);
   allbkg->Add(QCD);
   allbkg->Add(VV);
 // VV->Add(REST,VV);
 // QCD->Add(QCD,VV);
 // DY->Add(DY,QCD);
//  TT->Add(TT,DY);
   TT->Add(ST);
//  WJ->Add(WJ,TT);
  std::cout << "BKG : " << allbkg->GetSumOfWeights() << " : " << allbkg->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/allbkg->GetSumOfWeights() << "+/-" 
	    << TMath::Sqrt(histData->GetSumOfWeights())/allbkg->GetSumOfWeights() << std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////
    //ModTDRStyle();
    TH1D * bkgdErr = (TH1D*)allbkg->Clone("bkgdErr");
    ModTDRStyle();

  
    TCanvas* canv1 = new TCanvas("c1", "c1");
    canv1->cd();
    std::vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
    pads[0]->SetLogy(logY);
    histData->GetXaxis()->SetLabelSize(0.025);
    histData->GetXaxis()->SetTimeFormat("#splitline{%s}{%s}");
     for (int iB=1; iB<=nBins; ++iB) {histData->GetXaxis()->SetBinLabel(iB,BinLabels[iB-1]);}
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
    SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.4, 1.6);
    StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
    h[1]->GetYaxis()->SetNdivisions(4);
    h[1]->GetXaxis()->SetTitleOffset(1.2);
    h[1]->GetYaxis()->SetTitleOffset(1.0);
    pads[0]->cd();
    h[0]->GetYaxis()->SetTitleOffset(1.0);
    pads[1]->SetGrid(0,1);
    //it complains if the minimum is set to 0 and you try to set log y scale
    if(logY) h[0]->SetMinimum(0.09);
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
    histData->SetMarkerSize(1);


    int color = CreateTransparentColor(13,0.);
    histSignal->SetMarkerColor(2);
    histSignal->SetLineColor(color);
    histSignal->SetFillColor(color);
    histSignal->SetFillStyle(0);
    //histSignal->SetLineWidth(2);
    histSignal->SetMarkerStyle(31);
    histSignal->SetMarkerSize(1.5);
    
    if (!BlindData)legend->AddEntry(histData, "Observed", "ple");
    



    InitHist(QCD,mycolorqcd);
    InitHist(DY,mycolordyj);
    InitHist(TT,mycolortt);
    InitHist(ST,mycolortt);
    InitHist(TTX,mycolorttx);
    InitHist(WJ,mycolorwjet);
    InitHist(VV,mycolorvv);
    
    legend->AddEntry(TT,"TTJets/singleT","f");
    legend->AddEntry(TTX,"TTX/TG","f");
    legend->AddEntry(WJ,"WJets","f");
    legend->AddEntry(QCD,"QCD","f");
    legend->AddEntry(DY,"DYJets","f");
    legend->AddEntry(VV,"VV/VG/VVV","f");
    legend->AddEntry(histSignal,"C1N2_100_LSP_50","ple");



for (int iB=1; iB<=nBins; ++iB) {WJ->SetBinError(iB,0);TT->SetBinError(iB,0);DY->SetBinError(iB,0);QCD->SetBinError(iB,0);VV->SetBinError(iB,0);TTX->SetBinError(iB,0);}
    WJ->SetMarkerStyle(0);
    TT->SetMarkerStyle(0);
    DY->SetMarkerStyle(0);
    QCD->SetMarkerStyle(0);
    VV->SetMarkerStyle(0);
    TTX->SetMarkerStyle(0);

    THStack *hs = new THStack("","");
    hs->Add(VV);
    hs->Add(DY);
    hs->Add(QCD);
    hs->Add(WJ);
    hs->Add(TTX);
    hs->Add(TT);
    TH1D *hsum = ((TH1D*)(hs->GetStack()->Last()));
    for (int nb=0;nb<hsum->GetNbinsX()+1;++nb){if (hsum->GetBinContent(nb)==0) hsum->SetBinContent(nb,0.1);}
   // hsum->Draw("sameh");
    TT->Draw("sameh");
    TTX->Draw("sameh");
    WJ->Draw("sameh");
    QCD->Draw("sameh");   
    DY->Draw("sameh");
    VV->Draw("sameh");
    histSignal->Draw("sameh");

    canv1->Update();

    canv1->Update();
 /*   
    if (blindData)
    {
        for (int iB=nbMin; iB<=nbMax; ++iB)
        {
            histData->SetBinContent(iB,-1);
            histData->SetBinError(iB,0);
        }
    }
   */ 
    
    float errLumi = 0.062;
    //float errMuon = 0.03;
    //float errElectron = 0.04;
    for (int iB=1; iB<=nBins; ++iB) {
        //QCD->SetBinError(iB,0);
        //VV->SetBinError(iB,0);
        //TT->SetBinError(iB,0);
        WJ->SetBinError(iB,0);
        //DY->SetBinError(iB,0);
        //ZTT->SetBinError(iB,0);
        float eStat =  bkgdErr->GetBinError(iB);
        float X = bkgdErr->GetBinContent(iB);
        float eLumi = errLumi * X;
        //float eMuon = errMuon * X;
        //float eElectron = errElectron * X;
        //float eBkg = dummy->GetBinError(iB);
        //float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg+eMuon*eMuon+eElectron*eElectron);
	float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi);
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
        float x2 = WJ->GetBinContent(iB);
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
    ratioErrH->GetYaxis()->SetRangeUser(0.4,1.6);

    ratioErrH->Draw("e2same");
    ratioH->Draw("pe0same");

    pads[0]->cd();
    histData->Draw("pesame");
    
    FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
    DrawTitle(pads[0], "36.5 fb^{-1} (13 TeV)", 3);
    DrawTitle(pads[0], Variable, 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();
    canv1->Print(Variable+Suffix+suffix+".pdf");
}
}
