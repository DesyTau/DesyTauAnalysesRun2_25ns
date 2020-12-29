#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotDistribution(TString VariableName = "MuLTrk_DR", // variable name
                      TString filename = "test", // name of file with output (w/o ".root")
                      TString xtitle = "DR", // title of x axis
                      TString ytitle = "Normalized to unity", // title of y axis
                      float xLower = 0.001, // lower boundary of x axis
                      float xUpper = 0.999 // upper boundary of x axis
                      ) {
  
  xLower=xLower+0.001;
  xUpper=xUpper-0.001;
   
  TFile * file = new TFile(filename+".root");
  
  
  TH1D * SignalH_Old = (TH1D*)file->Get("InputVariables_Id/"+VariableName+"__Signal_Id");
  TH1D * BackgroundH_Old = (TH1D*)file->Get("InputVariables_Id/"+VariableName+"__Background_Id");
  
  ////// Rebinning /////
  int nBins = BackgroundH_Old->GetNbinsX();
  float xMin = BackgroundH_Old->GetBinLowEdge(1);
  float xMax = BackgroundH_Old->GetBinLowEdge(nBins+1);

  //std::cout << std::endl;
  //std::cout << "Histogram for Variable " << VariableName << " : " << "nbins = " << nBins
  //          << " , min = " << xMin
  //          << " , max = " << xMax << std::endl;
  //std::cout << std::endl;

  float bins[300];
  int nBinsNew = nBins;
  //std::cout << "New number of bins : ";
  //std::cin >> nBinsNew;

  //if (nBins%nBinsNew!=0) {
  //  std::cout << "new number of bins = " << nBinsNew
  //            << "  not multiple of " << nBins << std::endl;
  //  return;
  //}
  
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB) bins[iB] = xMin + float(iB)*binWidth;

  
  TH1D * SignalH = (TH1D*)TH1DtoTH1D(SignalH_Old,nBinsNew,bins,true,"_new");
  TH1D * BackgroundH = (TH1D*)TH1DtoTH1D(BackgroundH_Old,nBinsNew,bins,true,"_new");
  
  SignalH->Scale(1./SignalH->GetSumOfWeights());
  BackgroundH->Scale(1./BackgroundH->GetSumOfWeights());
  
  ////// Histos Colors ////////
  int tcolor_dark_2 = kAzure;
  int tcolor_light_2 = kAzure-2;
  
  int tcolor_dark_1 = kGreen+3;
  int tcolor_light_1 = kGreen-6;


  for (int iB=1; iB<=nBinsNew; iB++) 
  {

    SignalH->SetBinError(iB,0.);
    BackgroundH->SetBinError(iB,0.);
 
  }
  
  
  SignalH->Draw("h");
  BackgroundH->Draw("hsame");
  SignalH->Draw("hsame");
  
  
  SignalH->SetStats(0);
  SignalH->GetXaxis()->SetRangeUser(xLower,xUpper); 
  SignalH->GetXaxis()->SetLabelSize(0.04);
  SignalH->GetXaxis()->SetTitleOffset(1.2);
  SignalH->GetXaxis()->SetTitleSize(0.055);
  SignalH->GetXaxis()->SetTitle(xtitle);
  SignalH->GetXaxis()->SetTickLength(0.055);
  SignalH->GetXaxis()->SetTickSize(0.013);
  SignalH->GetYaxis()->SetTitleOffset(1.5);
  SignalH->GetYaxis()->SetTitle(ytitle);
  SignalH->GetYaxis()->SetTitleSize(0.045);
  SignalH->GetYaxis()->SetLabelSize(0.04);
  SignalH->GetYaxis()->SetTickLength(0.055);
  SignalH->GetYaxis()->SetTickSize(0.013);
  if(SignalH->GetMaximum() > BackgroundH->GetMaximum()) SignalH->GetYaxis()->SetRangeUser(0,1.1*SignalH->GetMaximum());
  else SignalH->GetYaxis()->SetRangeUser(0,1.1*BackgroundH->GetMaximum());
  SignalH->SetFillColorAlpha(tcolor_light_1,0.35);
  SignalH->SetLineColor(tcolor_dark_1);
  SignalH->SetLineWidth(1);
  SignalH->SetLineStyle(1);
  SignalH->SetFillStyle(3344);
  
  BackgroundH->SetStats(0);
  BackgroundH->SetFillColorAlpha(tcolor_light_2,0.35);
  BackgroundH->SetLineColor(tcolor_dark_2);
  BackgroundH->SetLineWidth(1);
  BackgroundH->SetLineStyle(1);
  BackgroundH->SetFillStyle(1001);
  
  


}

void PlotAll(TString filename = "test")
{

  int nvariables = 4;

  TString VarName[]= {"MuLMuT_DR", "MuLTrk_Mass", "MuLTrk_DR", "MuLTrkMuTTrk_Mass"}; 

  TString xTitle[] = {"#DeltaR_{ #mu_{L},#mu_{S}}", "m_{ #mu_{L},trk_{L}} [GeV]", 
                      "#DeltaR_{ #mu_{S},trk_{S}}", "m_{ #mu_{L},trk_{L},#mu_{S},trk_{S}} [GeV]"
                     };

  float	xLow[]     = {1.5, 0., 0.,0.};

  float xUp[]      = {5., 15., 1.2, 120.};

  gROOT->SetBatch();

  // Drawing in Canvas 
  TCanvas * c = new TCanvas("c","c",3000,3000) ;
  c->Divide(1,2);

  TPad *upper = (TPad*)c->cd(1);
  upper->SetPad("upper","upper",0,0.81,1,1,0,0,0);
  upper->Draw();
  upper->cd();

  
  TBox * PadBorder = new TBox(0.01,0.01,0.99,0.95);
  PadBorder->SetLineColor(kBlack);
  PadBorder->SetFillColor(kWhite);
  PadBorder->SetLineWidth(1);
  PadBorder->Draw("l");
  
  TBox * Signal = new TBox(0.1,0.1,0.9,0.9);
  Signal->SetLineColor(kGreen+3);
  Signal->SetFillColorAlpha(kGreen-6,0.35);
  Signal->SetFillStyle(3344);
  
  TBox * Background = new TBox(0.1,0.85,0.3,0.9);
  Background->SetLineColor(kAzure);
  Background->SetFillColorAlpha(kAzure-2,0.35);
  Background->SetFillStyle(1001);
  
  TText * CMS = new TText(.8,.5,"CMS");
  CMS->SetTextFont(61);
  CMS->SetTextSize(0.4);
  CMS->Draw("same");
  
  TText * Label = new TText(.8,.2,"Private work");
  Label->SetTextFont(52);
  Label->SetTextSize(0.15);
  Label->Draw("same");
 
  TLegend * leg = new TLegend(0.05,0.1,0.6,0.9);
  leg->SetTextSize(0.15);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->SetHeader("           #color[436]{#scale[1.2]{m_{a_{1}} = 10 GeV}}");
  leg->AddEntry(Signal,"#color[635]{Signal}","f");
  leg->AddEntry(Background,"#color[635]{Background}","f");
  leg->Draw();
 
  
  //upper->Modified();
  //upper->SetTicks();
  //upper->SetLeftMargin(0.13);
  //upper->SetRightMargin(0.05);
  //upper->SetBottomMargin(0.02);
  upper->Update();


  TPad * lower = (TPad*)c->cd(2);
  lower->SetPad("lower","lower",0,0,1,0.8,kWhite,0,0);
  lower->SetLeftMargin(0.35);
  lower->SetRightMargin(0.35);
  lower->Draw();

  
  TPad * lowPads = (TPad*)lower->cd();
  lowPads->Divide(2,2,0.001,0.001);

  for(int i=0;i<nvariables;i++) {

      lowPads->cd(i+1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetTopMargin(0.01);
      if(i==0 || i==2) gPad->SetRightMargin(0.01);
      if(i==1 || i==3) gPad->SetRightMargin(0.03);


      PlotDistribution(VarName[i], // variable name
                                                 filename, // name of file with output (w/o ".root")
                                                 xTitle[i], // title of x axis
                                                 "Normalized to unity", // title of y axis
                                                 xLow[i], // lower boundary of x axis
                                                 xUp[i] // upper boundary of x axis
                                                 );
  }
  

  //lower->Modified();
  //lower->SetTicks();
  //lower->SetLeftMargin(0.15);
  //lower->SetRightMargin(0.15);
  //lower->SetBottomMargin(0.15);
  

  c->Update();
  c->Print("Signal_vs_Background.pdf","Portrait pdf");


}
