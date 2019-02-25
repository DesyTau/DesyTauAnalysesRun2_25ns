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
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TPaveLabel.h"

double roundto3digits(double value){
    if (value == 0.0) return 0.0;

    double factor = pow(10.0, 3 - ceil(log10(fabs(value))));
    return round(value * factor ) / factor;   
}
void Plot_acott(
		 int firstDecay = 0,
		 int secondDecay = 0,
		 TString Cuts = "&&iso_1<0.15&&extraelec_veto<0.5&&byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5&&dilepton_veto<0.5&&extramuon_veto<0.5&&pt_1>20&&pt_2>20&&(singleLepTrigger>0.5||xTrigger>0.5)",//&&mva17_2>0.5&&mt_1<60
		 int nBins  =   20,
		 float xmin =    0,
		 float xmax =  2*TMath::Pi(),
		 TString xtitle = "#phi_{CP}",
		 TString ytitle = "Events",
		 TString directory = "/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/mutau/",
		 TString outputDir = "/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/Plots/",
		 TString Suffix = "GenCP_",        // for name of pdf
		 TString suffix = "",                      // for name of pdf
		 double lumi = 41465 //SingleMuon
		 ){

  TH1::SetDefaultSumw2();
  TString Variable = "acotautau";

  // weights
  // These samples names should match the root file names (-> sampleNames[i].root is read in later)
  TString sampleNames[3] = {
    "ggH_125", // signal (0)
    "ggA_130", // signal (1)
    "DYJetsToLL_M-50" // DY (2)
  };

  // Corresponding cross sections
  double xsec[3] = {1.,1., 5765.4}; //TO FIX, for now we are looking at normalized plots    


  // *******************************
  // ***** Selection Cuts    *******
  // *******************************


  TString Cut="";
  if(firstDecay==0)Cut+= "(genmode_1==0)";
  if(firstDecay==-1)Cut+= "(genmode_1<=6||genmode_1==8||genmode_1==9)";
  if(firstDecay==8)Cut+= "(genmode_1==8)";
  if(firstDecay==9)Cut+= "(genmode_1==9)";
  if(firstDecay==1)Cut+= "(genmode_1>=1&&genmode_1<=2)";
  if(firstDecay==2)Cut+= "(genmode_1>=4&&genmode_1<=5)";

  Cut+= "&&";
  if(secondDecay==0)Cut+= "(genmode_2==0)";
  if(secondDecay==-1)Cut+= "(genmode_2<=6||genmode_2==8||genmode_2==9)";
  if(secondDecay==8)Cut+= "(genmode_2==8)";
  if(secondDecay==9)Cut+= "(genmode_2==9)";
  if(secondDecay==1)Cut+= "(genmode_2>=1&&genmode_2<=2)";
  if(secondDecay==2)Cut+= "(genmode_2>=4&&genmode_2<=5)";

  Variable+="_";
  if(firstDecay==1)Variable+="1";
  else if(firstDecay==2)Variable+="2";
  else Variable+="0";
  if(secondDecay==1)Variable+="1";
  else if(secondDecay==2)Variable+="2";
  else Variable+="0";


  // *******************************
  // ***** Filling Histograms ******
  // *******************************

  TH1D * hist[3];

  int nSamples = 3;

  TCanvas * dummyCanv = new TCanvas("dummy","",600,600);

  TString ampl[3];


  // Draw main selection for all histograms in sampleNames
  for (int i=0; i<nSamples; ++i) {

    cout << endl << sampleNames[i] << ":" << endl;
    // Reading input file
    TFile * file = new TFile( directory + sampleNames[i] + ".root");

    // Get tree and one further important histogram from input file
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents"); //inputEventsH
    TTree * tree = (TTree*)file->Get("GenTauCheck");

    // Calculate normalization of this sample

    // Name and initialize histograms
    TString histName       = sampleNames[i] + Variable;
    hist[i]   = new TH1D(histName,"",nBins,xmin,xmax);

    cout << "Drawing ..." << endl;
    //    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histName,Cut);
    double norm = hist[i]->Integral(); 
    hist[i]->Scale((float)1./norm);
    cout << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1);
    double amplitude= 2* ( hist[i]->GetMaximum()-hist[i]->GetMinimum() ) / ( hist[i]->GetMaximum()+hist[i]->GetMinimum() );
    cout << " : Amplitude/Baseline = " << amplitude << endl;
    ampl[i]="";
    ampl[i]+=roundto3digits(amplitude);
  
  }
  cout << endl;
  delete dummyCanv;


  // ************************************
  // ***** Summarize backgrounds  *******
  // ************************************

  TH1D * ggH = (TH1D*)hist[0]   -> Clone("Scalar");
  TH1D * ggA = (TH1D*)hist[1]   -> Clone("Pseudoscalar");
  TH1D * DY  = (TH1D*)hist[2]   -> Clone("DY");
  /* This can be used in case of rescaling being applied to the sample, at the moment we work with normalized histograms so it has been removed
  for(int i=0;i<3;i++)
  cout << setw(2) << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) <<endl; 
  */

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING 
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();

  TCanvas* canv1 = new TCanvas("c1", "CP measurement in #tau_h#tau_h");
  canv1->cd();
  //vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
  

  

  // Setup legend
  TLegend *legend = PositionedLegend(0.40, 0.20, 3, 0.0);
  legend -> SetTextFont(42);

  ggH -> SetLineColor(2);
  ggH -> SetFillStyle(0);
  ggH -> SetLineWidth(3);
  ggA -> SetLineColor(4);
  ggA -> SetFillStyle(0);
  ggA -> SetLineWidth(3);
  DY  -> SetLineColor(1);
  DY  -> SetFillStyle(0);
  DY  -> SetLineWidth(3);

  legend->SetHeader("Sample: ampl/baseline");
  legend -> AddEntry(ggH, "Scalar :"       +ampl[0] , "l");
  legend -> AddEntry(ggA, "Pseudoscalar : "+ampl[1] , "l");
  legend -> AddEntry(DY , "DY : "          +ampl[2] , "l");

  // Add all bkg contributions to one stack plot

  canv1->Update();
  ggH->SetTitle("CP measurement in #tau#tau");
  ggA->SetTitle("CP measurement in #tau#tau");
  xtitle+= " channel ";
  xtitle+= TString::Itoa(firstDecay,10);
  xtitle+= TString::Itoa(secondDecay,10);
  ggH->GetXaxis()->SetTitle(xtitle);
  ggH->GetYaxis()->SetTitle("arbitrary units");
  ggH->GetYaxis()->SetRangeUser(0.,ggH->GetMaximum()*1.5);

  ggH->Draw();
  ggA->Draw("same");
  DY->Draw("same");
  legend->Draw();

  


  canv1->Update();




  //DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);

  //legend->Draw();
  FixOverlay();
  canv1->Update();

  canv1 -> Print( outputDir + Suffix + "acotautau_" + TString::Itoa(firstDecay,10) + TString::Itoa(secondDecay,10) + suffix + ".pdf" );
  canv1 -> Print( outputDir + Suffix + "acotautau_" + TString::Itoa(firstDecay,10) + TString::Itoa(secondDecay,10) + suffix + ".eps" );
  canv1 -> Print( outputDir + Suffix + "acotautau_" + TString::Itoa(firstDecay,10) + TString::Itoa(secondDecay,10) + suffix + ".png" );
}
