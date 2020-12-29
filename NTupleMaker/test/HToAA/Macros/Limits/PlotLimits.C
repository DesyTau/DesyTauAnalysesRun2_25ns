#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"

#include "Plotting.h"
#include "Plotting_Style.h"

void PlotLimits(char* fileList = "limits",
		bool blindData = true) {

  SetStyle();
  gStyle->SetOptFit(0000);
  gStyle->SetErrorX(0.5);
  gROOT->SetBatch();
  const int nPoints = 15;

  // signal strength limits sigma*BR / sigma*BR (at tanb=30)
  double mA[nPoints];      
  double minus2R[nPoints]; 
  double minus1R[nPoints]; 
  double medianR[nPoints]; 
  double plus1R[nPoints];  
  double plus2R[nPoints];  
  double obsR[nPoints];    

  double obs[nPoints];
  double minus2[nPoints];
  double minus1[nPoints];
  double median[nPoints];
  double plus1[nPoints];
  double plus2[nPoints];

  std::ifstream inputList(fileList);

  TString FileList(fileList);

  TString fileName;

  double MH;
  double LIMIT;

  int counter = 0;

  while (inputList >> fileName) {

    //    std::cout << fileName << std::endl;

    TFile * file = new TFile(fileName);

    TTree * tree = (TTree*)file->Get("limit");

    //    std::cout << "file : " << file << std::endl;
    //    std::cout << "tree : " << tree << std::endl;

    tree->SetBranchAddress("limit",&LIMIT);
    tree->SetBranchAddress("mh",&MH);

    tree->GetEntry(0);
    mA[counter] = float(MH);
    minus2R[counter] = float(LIMIT);

    //    std::cout << mA[counter] << std::endl;
    
    tree->GetEntry(1);
    minus1R[counter] = float(LIMIT);

    tree->GetEntry(2);
    medianR[counter] = float(LIMIT);

    tree->GetEntry(3);
    plus1R[counter] = float(LIMIT);

    tree->GetEntry(4);
    plus2R[counter] = float(LIMIT);

    tree->GetEntry(5);
    obsR[counter] = float(LIMIT);
    if (blindData)
      obsR[counter] = medianR[counter];

    counter++; 
      
  }


  std::cout << " m(Phi1)  -2s   -1s   exp   +1s   +2s   obs " << std::endl; 
  //           "100  24.1  28.2  33.8  40.8  48.2  23.0

  auto f_out = TFile::Open("Limits_tree.root","RECREATE");
  TNtuple Limits("Limits","Limits","m_a:minus2:minus1:median:plus1:plus2:obs");

  for (int i=0; i<counter; ++i) {

    obs[i]    = obsR[i];
    minus2[i] = minus2R[i];
    minus1[i] = minus1R[i];
    median[i] = medianR[i];
    plus1[i]  = plus1R[i];
    plus2[i]  = plus2R[i];

    Limits.Fill(int(mA[i]),minus2[i],minus1[i],median[i],plus1[i],plus2[i],obs[i]);

    char strOut[200];
    sprintf(strOut,"%3i  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f",
	    int(mA[i]),minus2[i],minus1[i],median[i],plus1[i],plus2[i],obs[i]);
    std::cout << strOut << std::endl;

  }
  f_out->Write();
  
  double zeros[nPoints];
  double upper[nPoints];
  double lower[nPoints];
  double central[nPoints];
  for (int i=0; i<counter; ++i) {
    zeros[i] = 0;
    central[i] = 15; 
    minus2[i] = median[i] - minus2[i];
    minus1[i] = median[i] - minus1[i];
    plus1[i]  = plus1[i]  - median[i];
    plus2[i]  = plus2[i]  - median[i];
    upper[i] = 15 - central[i];
    lower[i] = central[i] - obs[i];
  }
  
  
  int nPointsX = counter;

  TGraph * obsG = new TGraph(nPointsX, mA, obs);
  obsG->SetLineWidth(3);
  obsG->SetLineColor(1);
  obsG->SetLineWidth(3);
  obsG->SetMarkerColor(1);
  obsG->SetMarkerStyle(0);
  obsG->SetMarkerSize(0);

  TGraph * expG = new TGraph(nPointsX, mA, median);
  expG->SetLineWidth(3);
  expG->SetLineColor(2);
  expG->SetLineStyle(2);
  
  TGraphAsymmErrors * observed = new TGraphAsymmErrors(nPointsX, mA, central, zeros, zeros, lower, upper);
  observed->SetFillColor(kCyan-4);
  observed->SetLineWidth(3);

  TGraphAsymmErrors * innerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus1, plus1);
  innerBand->SetFillColor(kGreen+1);
  innerBand->SetLineColor(0);

  TGraphAsymmErrors * outerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus2, plus2);
  outerBand->SetFillColor(kOrange);
  outerBand->SetLineColor(0);

  double mAExcl[2] = {4,21};
  double x0Excl[2] = {0.6,0.6};
  double zeroExcl[2] = {0,0};
  double lowerExcl[2] = {0.26,0.26};
  double upperExcl[2] = {0.3,0.3};

  TGraphAsymmErrors * excluded = new TGraphAsymmErrors(2,mAExcl,x0Excl,zeroExcl,zeroExcl,lowerExcl,upperExcl);
  excluded->SetMarkerStyle(0);
  excluded->SetMarkerSize(0);
  int new_idx = CreateTransparentColor(kCyan,0.4);
  excluded->SetFillColor(new_idx);
  excluded->SetFillStyle(1001);
  excluded->SetLineColor(0);

  TH2F * frame = NULL;

  //  frame = new TH2F("frame","",2,100,500,2,0,70);
  frame = new TH2F("frame","",2,4,21,2,0,0.7);
  frame->GetXaxis()->SetTitle("m_{a_{1}} [GeV]");
  frame->GetYaxis()->SetTitle(" #sigma B / #sigma_{SM}");
  frame->SetMarkerStyle(0);
  frame->SetMarkerColor(0);
  frame->SetMarkerSize(0);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNdivisions(206);
  frame->GetYaxis()->SetTitleOffset(1.55);  
  frame->GetYaxis()->SetTitleSize(0.048);  
  

  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  frame->Draw();

  outerBand->Draw("3same");
  innerBand->Draw("3same");
  expG->Draw("lsame");
  if (!blindData)
    obsG->Draw("lpsame");
  excluded->Draw("3same");

  float xLeg = 0.18;
  float yLeg = 0.83;
  float xLegend = 0.57;
  float yLegend = 0.41;
  float sizeLeg = 0.27;

  TLegend * leg = new TLegend(0.25,0.55,0.6,0.85);
  leg->SetHeader("95% CL upper limits"); //,"C");
  leg->SetFillColor(0);
  leg->SetTextSize(0.02);
  leg->SetBorderSize(0);
  if (!blindData) 
    leg->AddEntry(obsG,"Observed","lp");
  leg->AddEntry(expG,"Expected","l");
  leg->AddEntry(innerBand,"68% expected","f");
  leg->AddEntry(outerBand,"95% expected","f");
  leg->AddEntry(excluded,"Excluded by ATLAS-CMS","f");
  leg->AddEntry(frame,"combined coupling analysis","p");
  leg->Draw();

  extraText = "Work in Progress";
  writeExtraText = true;
  CMS_lumi(canv,4,33); 
  canv->RedrawAxis();

  leg->Draw();
  canv->Update();
  //  TString suffix(fileList);
  canv->Print("BR_limits.pdf","Portrait pdf");
  canv->Print("BR_limits.png");

}
