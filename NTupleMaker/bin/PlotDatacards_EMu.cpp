using namespace std;

#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
#include "HttStylesNew.cc"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <map>
#include <vector>

void Plot(
	  TFile * file,
	  TString  era,
	  TString category
	  ) {

  gStyle->SetOptStat(0000);

  bool plotLegend = true;
  bool legRight = true;
  bool logY = false;
  bool logX = true;
  float yLower = 20;
  float scaleYUpper = 50;
  TString xtitle = "m_{T}^{tot}";
  TString ytitle = "";


  TH1D * histData = (TH1D*)file->Get(category+"/data_obs");
  if (histData==NULL) {
    std::cout << "  era/category " << era << "/" << category << " does not exist" << std::endl;
    std::cout << "  nothing is done... " << std::endl;
    std::cout << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;
    return;  
  }
  else {
    std::cout << std::endl;
    std::cout << "processing category era/category : " << era << "/" << category << std::endl;
    std::cout << std::endl;
  }
  TH1D * TT  = (TH1D*)file->Get(category+"/TTL");
  TH1D * ZTT = (TH1D*)file->Get(category+"/EMB");
  TH1D * ZLL = (TH1D*)file->Get(category+"/ZL");
  TH1D * W   = (TH1D*)file->Get(category+"/W");
  TH1D * EWK = (TH1D*)file->Get(category+"/VVL"); 
  TH1D * QCD = (TH1D*)file->Get(category+"/QCD");

  int nBins = histData->GetNbinsX();

  std::cout << "Top : " << TT->GetSumOfWeights() << std::endl;
  std::cout << "EWK : " << EWK->GetSumOfWeights() << std::endl;
  std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
  std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << std::endl;

  //  adding normalization systematics
  double ZTT_norm = 0.04; //  normalization ZTT :  4% (EMBEDDED)
  double EWK_norm = 0.05; //  normalization EWK :  5%
  double QCD_norm = 0.10; //  normalization Fakes : 10%
  double ZLL_mtau = 0.02; //  mu->tau fake rate ZLL : 2%
  double TT_norm  = 0.06; //  normalization TT  :  6%
  double W_norm   = 0.06; //  normalization W   :  6%

  double eff_Emb = 0.04;
  double eff_MC  = 0.04;

  for (int iB=1; iB<=nBins; ++iB) {

    float ztt  = ZTT->GetBinContent(iB);
    float ztte = ZTT->GetBinError(iB);
    ztte = TMath::Sqrt(ztte*ztte+ztt*ztt*(ZTT_norm*ZTT_norm+eff_Emb*eff_Emb));
    ZTT->SetBinError(iB,ztte);

    float ewk  = EWK->GetBinContent(iB);
    float ewke = EWK->GetBinError(iB);
    ewke = TMath::Sqrt(ewke*ewke+ewk*ewk*(EWK_norm*EWK_norm+eff_MC*eff_MC));
    EWK->SetBinError(iB,ewke);

    float qcd  = QCD->GetBinContent(iB);
    float qcde = QCD->GetBinError(iB);
    qcde = TMath::Sqrt(qcde*qcde+qcd*qcd*QCD_norm*QCD_norm);
    QCD->SetBinError(iB,qcde);
    if (qcd<0) {
      QCD->SetBinContent(iB,0.);
      QCD->SetBinError(iB,0.);
    }
    //    std::cout << "bin : " << iB << " : " << QCD->GetBinContent(iB) << std::endl;

    float w = W->GetBinContent(iB);
    float we = W->GetBinError(iB);
    we = TMath::Sqrt(we*we+w*w*(W_norm*W_norm+eff_MC*eff_MC));
    W->SetBinError(iB,we);

    float tt  = TT->GetBinContent(iB);
    float tte = TT->GetBinError(iB);
    tte = TMath::Sqrt(tte*tte+tt*tt*TT_norm*TT_norm);
    TT->SetBinError(iB,tte);

    float zll  = ZLL->GetBinContent(iB);
    float zlle = ZLL->GetBinError(iB);
    zlle = TMath::Sqrt(zlle*zlle+zll*zll*eff_MC*eff_MC);
    ZLL->SetBinError(iB,zlle);

    /*
    std::cout << iB << " : " 
	      << "qcd = " << qcd << " +/- " << qcde 
	      << "  w = " << w << " +/- " << we 
	      << "  ztt = " << ztt << " +/- " << ztte
	      << "  zll = " << zll << " +/- " << zlle 
	      << "  ewk = " << ewk << " +/- " << ewke << std::endl;
    */
  }

  EWK->Add(EWK,TT);
  W->Add(W,EWK);
  QCD->Add(QCD,W);
  ZLL->Add(ZLL,QCD);
  ZTT->Add(ZTT,ZLL);

  std::cout << std::endl;
  std::cout << "Model : " << ZTT->GetSumOfWeights() << std::endl;
  std::cout << "Data  : " << histData->GetSumOfWeights() << std::endl;
  std::cout << "ratio = " << histData->GetSumOfWeights()/ZTT->GetSumOfWeights() << std::endl;
  std::cout << std::endl;

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  //  TH1
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  InitData(histData);

  InitHist(QCD,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(EWK,"","",TColor::GetColor("#6F2D35"),1001);
  InitHist(W,"","",TColor::GetColor("#DE5A6A"),1001);
  InitHist(ZLL,"","",TColor::GetColor("#4496C8"),1001);
  InitHist(ZTT,"","",TColor::GetColor("#FFCC66"),1001);

  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetXaxis()->SetRangeUser(61,3900);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  float yUpper = histData->GetMaximum();
  if (ZTT->GetMaximum()>yUpper)
    yUpper = ZTT->GetMaximum();
  histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);
  ZTT->GetYaxis()->SetRangeUser(0,1.2*ZTT->GetMaximum());
  if (logY) {
    histData->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
    ZTT->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
  }

  histData->SetMarkerSize(1.4);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
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
  ZTT->Draw("sameh");
  ZLL->Draw("sameh");
  QCD->Draw("sameh");
  W->Draw("sameh");
  //  EWK->Draw("sameh");
  TT->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  histData->Draw("e1same");
  float chi2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = W->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2/ndof = " << chi2 << "/" << nBins << std::endl;
  std::cout << "Prob = " << TMath::Prob(chi2,double(nBins)) << std::endl;
  std::cout << std::endl;

  TLegend * leg;
  if (legRight) leg = new TLegend(0.7,0.40,0.90,0.70);
  else leg = new TLegend(0.20,0.40,0.50,0.70);
  SetLegendStyle(leg);
  leg->SetTextSize(0.04);
  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(ZTT,"embedded #tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow ll","f");
  leg->AddEntry(QCD,"QCD","f");
  leg->AddEntry(W,"electroweak","f");
  //  leg->AddEntry(EWK,"electroweak","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  //  leg->AddEntry(histBBH,"bbH(1200)","l");
  //  leg->AddEntry(histGGH,"ggH(1200)","l");
  if (plotLegend) leg->Draw();

  if (logY) upper->SetLogy(true);
  if (logX) upper->SetLogx(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.301,1.699);
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

  // ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
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
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  if (logX) lower->SetLogx(true);
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Update();
  canv1->Print(category+"_"+era+".png");

  std::cout << std::endl;
  std::cout << "Signal yields for 1 pb... " << std::endl;
  for (auto mass : masses) {
    TH1D * bbH = (TH1D*)file->Get(category+"/bbH_"+mass);
    TH1D * ggHt = (TH1D*)file->Get(category+"/ggH_t_"+mass);
    TH1D * ggHb = (TH1D*)file->Get(category+"/ggH_b_"+mass);
    TH1D * ggHi = (TH1D*)file->Get(category+"/ggH_i_"+mass);
    TH1D * ggAt = (TH1D*)file->Get(category+"/ggA_t_"+mass);
    TH1D * ggAb = (TH1D*)file->Get(category+"/ggA_b_"+mass);
    TH1D * ggAi = (TH1D*)file->Get(category+"/ggA_i_"+mass);
    std::cout << "mass=" << mass << " : ";
    if (bbH!=NULL) 
      std::cout << "  bbH = " << bbH->GetSumOfWeights();
    else 
      std::cout << "  bbH --- ";
    if (ggHt!=NULL&&ggHb!=NULL&&ggHi!=NULL) 
      std::cout << "  ggH_{t,b,i} = " 
		<< "{" << ggHt->GetSumOfWeights() 
		<< "," << ggHb->GetSumOfWeights()
		<< "," << ggHi->GetSumOfWeights()
		<< "}";
    else 
      std::cout << "  ggH --- ";

    if (ggAt!=NULL&&ggAb!=NULL&&ggAi!=NULL) 
      std::cout << "  ggA_{t,b,i} = " 
		<< "{" << ggAt->GetSumOfWeights() 
		<< "," << ggAb->GetSumOfWeights()
		<< "," << ggAi->GetSumOfWeights()
		<< "}";
    else 
      std::cout << "  ggA --- ";


    std::cout << std::endl;
  }

  std::cout << std::endl; 
  std::cout << "+++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;

}

int main(int argc, char * argv[]) {

  TString input_dir = "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/datacards";

  vector<TString> eras = {"2016","2017","2018"};
  vector<TString> categories = 
    {
      //      "em_inclusive",
      "em_DZetaLtm35",
      "em_Nbtag0_DZetam35Tom10",
      "em_Nbtag0_DZetam10To30",
      "em_Nbtag0_DZetaGt30",
      "em_NbtagGt1_DZetam35Tom10",
      "em_NbtagGt1_DZetam10To30",
      "em_NbtagGt1_DZetaGt30"
    };

  for (auto era : eras) {
    TString rootFileName = input_dir+"/"+era+"/htt_em_mssm.root";
    TFile * file = new TFile(rootFileName);
    if (file->IsZombie()) {
      std::cout << "cannot open RooT file " << rootFileName << std::endl;
      std::cout << "skipping..." << std::endl;
      continue;
    }
    else {
      std::cout << "processing era " << era << std::endl;
    }
    for (auto category : categories) {
      Plot(file,era,category);	
    }
    delete file;
  }

}
