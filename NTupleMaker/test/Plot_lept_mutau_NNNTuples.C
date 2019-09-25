///////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Plotting macro for NN NTuples
///    Features: FF method, NN output categorization, unrolled distribution for acoplanarity angle
///    Author: Andrea Cardini ( andrea.cardini@desy.de , andrea.cardini@cern.ch )
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"



void Plot_lept_mutau_NNNTuples(TString Variable = "mt_1",
			       TString xtitle = "m_{T} [GeV]",
			       int nBins  =   30,
			       float xmin =    0,
			       float xmax =  150,
			       TString Weight = "xsec_lumi_weight*puweight*effweight*mcweight*",
			       TString Cut="(mt_1<50)*",
			       int categoryIndex=-1,
			       TString ytitle = "Events",
			       TString directory = "HtautauCP_mutau/Inputs/NTuples_mt_",
			       TString outputDir = "./Plots/",
			       int year=2017,
			       bool FFmethod = false,  
			       bool LargeScale = false,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = true,
			       int scaleSignal = 1.,
			       bool blindData = false,
			       bool FORCE = false
			       )
{
  using namespace std;
  if(categoryIndex>=0){
    Cut+="(predicted_class=="+TString::Itoa(categoryIndex,10)+")*";
  }
  TString yVar,xVar;
  bool var2D=false;
  if(Variable.Contains(":")){
    TObjArray *array = Variable.Tokenize(":");
    yVar=array->At(0)->GetName();
    xVar=array->At(1)->GetName();
    TString range=to_string(xmax-xmin).c_str();
    nBins=nBins*3;
    xmax+=2*(xmax-xmin);
    Variable=xVar+"*("+yVar+"<0.4)+("+range+"+"+xVar+"*("+yVar+">=0.4&&"+yVar+"<0.7))+(2*"+range+"+"+xVar+"*("+yVar+">=0.7))";
    var2D=true;
  }
  TString IsoCut=Cut+"(mva17_2>0.5)*";
  TString AntiIsoCut=Cut+"(mva17_2<0.5)*";
  directory="$CMSSW_BASE/src/"+directory;
  directory+=year;
  directory+="/";
  TH1::SetDefaultSumw2();
  SetStyle();
  const int nSamples = 10; //DY is used twice, for Zll and Ztt

  // Simple check to prevent accidental unblinding in SR, to unblind set FORCE to True
  if(!FORCE&&(categoryIndex==0||categoryIndex==1)){
    blindData=true;
  }

  // weights
  TString topweight("*1");
  TString Wjets_weight("*0.96");
  TString qcdweight("1.06*");
  TString zptmassweight="*1.0";                  //TO DO: CHANGE WEIGHTs
 
  TString sampleNames[nSamples] = {
        "mt-NOMINAL_ntuple_Data", // data (0)
        "mt-NOMINAL_ntuple_DY",// (1)Drell-Yan Z->TT
        "mt-NOMINAL_ntuple_DY", // (2)Drell-Yan Z->LL
        "mt-NOMINAL_ntuple_WJets",// (3)WJets
        "mt-NOMINAL_ntuple_TT",//(4)TTbar leptonic, hadronic, + semileptonic
        "mt-NOMINAL_ntuple_SingleTop", // (5) SingleTop tW tbar, SingleTop tW t, SingleTop t antitop, SingleTop t top
        "mt-NOMINAL_ntuple_VV",// (6) WW, WZ, ZZ
	"mt-NOMINAL_ntuple_ggH", // (7) Scalar ggH
	"mt-NOMINAL_ntuple_qqH", // (8) Scalar VBF H
	"mt-NOMINAL_ntuple_CPodd" // (9) Pseudoscalar 
  };

  cout<<"this are the samples"<<endl;
  for (int i=0; i<nSamples; ++i) {
    cout << endl << sampleNames[i] << ":" << endl;}

  // *******************************
  // ***** Selection Cuts    *******
  // *******************************

  TString cuts[nSamples];
  TString cutsSS[nSamples];
  TString cutsaIso[nSamples];

  // MC specific cuts to select certain type of particle
  TString isZTT="*(gen_match_2==5)";
  TString isZLL="*!(gen_match_2==5)";

  // Selection cuts applied to all samples
  for (int i=0; i<nSamples; ++i){
    cuts[i]   = IsoCut+Weight+"(os>0.5)";
    cutsSS[i] = IsoCut+Weight+qcdweight+"(os<0.5)";
    cutsaIso[i] = AntiIsoCut+Weight+"(os>0.5)*ff_nom*((gen_match_2==5)*0.88+(gen_match_2==2||gen_match_2==4)*correction_againstMuonTight3_2+(gen_match_2==1||gen_match_2==3))";
  }

  //specific selection weights for data, DY and top
  cuts[0] = IsoCut+"(os>0.5)"; //DATA
  cuts[1] += zptmassweight+isZTT;
  cuts[2] += zptmassweight+isZLL;
  cuts[3] += Wjets_weight;
  cuts[4] += topweight; 
  cuts[7] += "*"+TString::Itoa(scaleSignal,10);
  cuts[8] += "*"+TString::Itoa(scaleSignal,10);
  cuts[9] += "*"+TString::Itoa(scaleSignal,10);
  if(FFmethod) for(int i=2; i<nSamples; ++i) cuts[i] += "*(gen_match_2!=6)";
   
  cutsSS[0] = IsoCut+qcdweight+"(os<0.5)";
  cutsSS[1] += zptmassweight+isZTT;
  cutsSS[2] += zptmassweight+isZLL;
  cutsSS[3] += Wjets_weight; 
  cutsSS[4] += topweight; 
  cutsSS[7] += "*"+TString::Itoa(scaleSignal,10);
  cutsSS[8] += "*"+TString::Itoa(scaleSignal,10);
  cutsSS[9] += "*"+TString::Itoa(scaleSignal,10);
 
  cutsaIso[0] = AntiIsoCut+"(os>0.5)*ff_nom";
  cutsaIso[1] += zptmassweight+isZTT;
  cutsaIso[2] += zptmassweight+isZLL;
  cutsaIso[3] += Wjets_weight; 
  cutsaIso[4] += topweight; 
  cutsaIso[7] += "*"+TString::Itoa(scaleSignal,10);
  cutsaIso[8] += "*"+TString::Itoa(scaleSignal,10);
  cutsaIso[9] += "*"+TString::Itoa(scaleSignal,10);
 
  // *******************************
  // ***** Filling Histograms ******
  // *******************************

  TH1D * hist[nSamples];
  TH1D * histSS[nSamples];
  TH1D * hist_AntiIso[nSamples];
  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  // Draw main selection for all histograms in sampleNames
  for (int i=0; i<nSamples; ++i) {

    // Reading input file
    TFile * file = new TFile( directory + sampleNames[i] + ".root");
    TTree * tree = (TTree*)file->Get("TauCheck"); 
    
    // Name and initialize histograms
    TString histName   = sampleNames[i] + Variable + "_os";
    TString histNameSS = sampleNames[i] + Variable + "_ss";
    TString histNameaIso   = sampleNames[i] + Variable + "_aIso";

    hist[i]   = new TH1D(histName,"",nBins,xmin,xmax);
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
    hist_AntiIso[i] = new TH1D(histNameaIso,"",nBins,xmin,xmax);
    cout << "Drawing ..." << endl;
    tree->Draw(Variable+">>"+histName,cuts[i]);
    cout << cuts[i] <<endl;
    tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
    tree->Draw(Variable+">>"+histNameaIso,cutsaIso[i]);
    cout << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }


  cout << endl;
  delete dummyCanv;



  // ********************************************
  // ***** Adding similar backgrounds     *******
  // ********************************************

  // Adding up single top and VV backgrounds
  hist[5]->Add(hist[5],hist[6]);
  
  // Adding SS backgrounds
  for (int iH=2; iH<6; iH++) histSS[1]->Add(histSS[iH]);

  // Adding antiIso backgrounds
  for (int iH=2; iH<6; iH++) hist_AntiIso[1]->Add(hist_AntiIso[iH]);
  
  // ********************************************
  // ***** QCD background estimation      *******
  // ********************************************

  float dataSS     = histSS[0]->GetSumOfWeights();
  float dataSSfull = histSS[0]->Integral(0,nBins+1);
  
  // Subtracting background from SS
  histSS[0]->Add(histSS[0],histSS[1],1,-1);

  // Subtracting true taus from anti-Iso data
  hist_AntiIso[0]->Add(hist_AntiIso[0],hist_AntiIso[1],1,-1);

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

  TH1D * histData   = (TH1D*)hist[0]         -> Clone("data_obs");
  TH1D * QCD        = (TH1D*)histSS[0]       -> Clone("QCD");
  TH1D * ZTT        = (TH1D*)hist[1]         -> Clone("ZTT");
  TH1D * ZLL        = (TH1D*)hist[2]         -> Clone("ZLL");
  TH1D * W          = (TH1D*)hist[3]         -> Clone("W");
  TH1D * TT         = (TH1D*)hist[4]         -> Clone("TT"); 
  TH1D * VV         = (TH1D*)hist[5]         -> Clone("VV"); 
  TH1D * ggH        = (TH1D*)hist[7]         -> Clone("ggH");
  TH1D * qqH        = (TH1D*)hist[8]         -> Clone("qqH");
  TH1D * CPoddH     = (TH1D*)hist[9]         -> Clone("CPoddH");
  TH1D * Fakes      = (TH1D*)hist_AntiIso[0] -> Clone("Fakes");
  
  for(int i=0;i<nSamples;i++){
    cout << setw(nSamples) << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }

  cout << endl;
  cout << "QCD : Sum of weights = " << QCD->GetSumOfWeights() << " : Integral = " << QCD->Integral(1,nBins+1) << endl;
  cout << "VV  : Sum of weights = " << VV->GetSumOfWeights()  << " : Integral = " << VV->Integral(1,nBins+1)  << endl;
  cout << "W   : Sum of weights = " << W->GetSumOfWeights()   << " : Integral = " << W->Integral(1,nBins+1)   << endl;
  cout << "TT  : Sum of weights = " << TT->GetSumOfWeights()  << " : Integral = " << TT->Integral(1,nBins+1)  << endl;
  cout << "ZLL : Sum of weights = " << ZLL->GetSumOfWeights() << " : Integral = " << ZLL->Integral(1,nBins+1) << endl;
  cout << "ZTT : Sum of weights = " << ZTT->GetSumOfWeights() << " : Integral = " << ZTT->Integral(1,nBins+1) << endl;
  cout << "Fakes : Sum of weights = " << Fakes->GetSumOfWeights() << " : Integral = " << Fakes->Integral(1,nBins+1) << endl;


  if(FFmethod){
    cout << endl;
  cout << "************************" << endl;
  cout << "True Taus = " << VV->Integral(1,nBins+1)+W->Integral(1,nBins+1)+TT->Integral(1,nBins+1)+ZLL->Integral(1,nBins+1)+ZTT->Integral(1,nBins+1) << endl;
  cout << "True Taus + Fakes = " << VV->Integral(1,nBins+1)+W->Integral(1,nBins+1)+TT->Integral(1,nBins+1)+ZLL->Integral(1,nBins+1)+ZTT->Integral(1,nBins+1)+Fakes->Integral(1,nBins+1) << endl;
  cout << "************************" << endl;
  cout << endl;

  }else {
    cout << endl;
  cout << "************************" << endl;
  cout << "MC = " << VV->Integral(1,nBins+1)+W->Integral(1,nBins+1)+TT->Integral(1,nBins+1)+ZLL->Integral(1,nBins+1)+ZTT->Integral(1,nBins+1) << endl;
  cout << "MC + QCD = " << VV->Integral(1,nBins+1)+W->Integral(1,nBins+1)+TT->Integral(1,nBins+1)+ZLL->Integral(1,nBins+1)+ZTT->Integral(1,nBins+1)+QCD->Integral(1,nBins+1) << endl;
  cout << "************************" << endl;
  cout << endl;


  }


  cout << endl;
  cout << "************************" << endl;
  cout << "Fake taus = " << VV->Integral(1,nBins+1)+W->Integral(1,nBins+1)+QCD->Integral(1,nBins+1) << endl;
  cout << "************************" << endl;
  cout << endl;

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
    ggH   -> SetBinError(iB,0);
    CPoddH  -> SetBinError(iB,0);
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING 
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();
  TCanvas* canv1;
  if(LargeScale)canv1 = new TCanvas("c1", "c1", 2000,800);
  else canv1 = new TCanvas("c1", "c1", 1200,1000);
  canv1->cd();
  vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
  pads[0]->SetLogy(logY);

  vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
  //h[0]->GetYaxis()->SetMaxDigits(3);
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
  h[1] -> GetXaxis()->SetTitleOffset(0.95);
  h[1] -> GetXaxis()->SetNdivisions(505);
  h[1] -> GetYaxis()->SetTitleOffset(1.1);
  if(LargeScale)h[1] -> GetYaxis()->SetTitleOffset(0.9);
  pads[0] -> cd();
  h[0] -> GetYaxis()->SetTitleOffset(1.6);
  if(LargeScale)h[0] -> GetYaxis()->SetTitleOffset(1);
  pads[1] -> SetGrid(0,1);
  //it complains if the minimum is set to 0 and you try to set log y scale
  if(logY) h[0] -> SetMinimum(1);
  pads[0] -> cd();

  // Setup legend
  TLegend *legend = PositionedLegend(0.25, 0.30, 3, 0.03);
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
  InitHist(Fakes,TColor::GetColor("#c6f74a"));

  legend -> AddEntry(histData, "Observed", "ple");
  legend -> AddEntry(ZTT,"Z#rightarrow #tau#tau","f");
  legend -> AddEntry(TT,"t#bar{t}","f");
  legend -> AddEntry(ZLL,"Z#rightarrow #mu#mu/ee","f");
  if(FFmethod)legend -> AddEntry(Fakes,"j#rightarrow #tau_{h}","f");
  else{
  legend -> AddEntry(QCD,"QCD","f");
  
  legend -> AddEntry(W,"W+jets","f");
  }
  legend -> AddEntry(VV,"single top + diboson","f");
  
  // Add all bkg contributions to one stack plot
  THStack *stack = new THStack("Background","");
  if(FFmethod)stack->Add(Fakes);
  else{
    stack -> Add(QCD);
    stack -> Add(W);
  }
  stack -> Add(VV);
  stack -> Add(ZLL);
  stack -> Add(TT);
  stack -> Add(ZTT);
  stack -> Draw("hist same");

  canv1->Update();

  InitSignal(ggH ,2);
  InitSignal(qqH ,3);
  InitSignal(CPoddH,4);
  if (showSignal){
    if(compareCP){
      legend->AddEntry(ggH,"CP-even Higgs(125) #times "+TString::Itoa(scaleSignal,10),"f");
      ggH->Draw("hist same");
      legend->AddEntry(CPoddH,"CP-odd Higgs(125) #times "+TString::Itoa(scaleSignal,10),"f");
      CPoddH->Draw("hist same");
    }else{
      legend->AddEntry(ggH,"ggH Higgs(125) #times 1","f");
      ggH->Draw("hist same");
      legend->AddEntry(qqH,"VBF Higgs(125) #times 1","f");
      qqH->Draw("hist same");
    }
  }

  canv1->Update();

  if (blindData)
    {
      for (int iB=0; iB<=nBins; ++iB)
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
    Fakes->SetBinError(iB,0);
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

  if (year== 2016) DrawTitle(pads[0], "35.9 fb^{-1} (13 TeV, 2016)", 3);
  if (year== 2017) DrawTitle(pads[0], "41.9 fb^{-1} (13 TeV, 2017)", 3);
  if (year== 2018) DrawTitle(pads[0], "59.7 fb^{-1} (13 TeV, 2018)", 3);
  
  DrawTitle(pads[0], "#scale[1.2]{         #bf{CMS} Work in progress}", 1);
  FixBoxPadding(pads[0], legend, 0.05);
  //legend->SetNColumns(2);
  legend->Draw();
  FixOverlay();
  canv1->Update();
  pads[0]->GetFrame()->Draw();
  if(var2D)Variable="Unrolled_"+xVar;
  canv1 -> Print( outputDir + "MuTau_" + Variable + (categoryIndex>=0 ? ((TString)"_Cat"+TString::Itoa(categoryIndex,10)+"_") : (TString)"_") +(FFmethod ? (TString)"fakes" : (TString)"MC") + ".pdf" );
  canv1 -> Print( outputDir + "MuTau_" + Variable + (categoryIndex>=0 ? ((TString)"_Cat"+TString::Itoa(categoryIndex,10)+"_") : (TString)"_") +(FFmethod ? (TString)"fakes" : (TString)"MC") + ".png" );
}
