#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
// #include "boost/lexical_cast.hpp"
// #include "boost/algorithm/string.hpp"
// #include "boost/format.hpp"
// #include "boost/program_options.hpp"
// #include "boost/range/algorithm.hpp"
// #include "boost/range/algorithm_ext.hpp"
// #include "Plotting.h"
// #include "Plotting_Style.h"
// #include "HttStylesNew.cc"
#include "TPad.h"
#include "TH2.h"
#include "TH1.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"



void ComputeFakeFractions_2dbin(TString directory = "./store2/",
				TString outputDir = "./control_plots/",
				TString Variable = "m_vis",
				TString suffix = "_antiIso",
				double lumi = 41465 ,
				TString weight = "weight*",//Ditau,
				int nbins = 8,
				float bins[10],
				int hbin = 150
 				) {
  int nSamples = 4;
  bins = {0,50,80,110,150,200,250,300,1000};
  TString samples[4]={"DATA_Tau_tt_Sync.root","DYJetsToLL_tt_Sync.root","WJetsToLNu_tt_Sync.root","TTTo2L2Nu_tt_Sync.root"};
  TString AntiIso_1 = "os>0.5 &&extraelec_veto<0.5 && extramuon_veto<0.5 && pt_1>20&&pt_2>20 && mva17_1<0.5 && mva17_2>0.5 && ";
  TString AntiIso_2 = "os>0.5 &&extraelec_veto<0.5 && extramuon_veto<0.5 && pt_1>20&&pt_2>20 && mva17_2<0.5 && mva17_1>0.5 && ";
  TString AntiIso_1_d = "(os>0.5 &&extraelec_veto<0.5 && extramuon_veto<0.5 && pt_1>20&&pt_2>20 && mva17_1<0.5 && mva17_2>0.5 )";
  TString AntiIso_2_d = "(os>0.5 &&extraelec_veto<0.5 && extramuon_veto<0.5 && pt_1>20&&pt_2>20 && mva17_2<0.5 && mva17_1>0.5 )";
  double xsec[4] = {1,6225.42,61526,88.29};
  
  TString cutsF[4],cutsT[4];
  cutsF[0] ="("+AntiIso_1_d+"||"+AntiIso_2_d+")";
  cutsF[1] =weight+"0.5*(("+AntiIso_1+"gen_match_1==6"+")||("+AntiIso_2+"gen_match_2==6"+"))";
  cutsF[2] =weight+"0.5*(("+AntiIso_1+"gen_match_1==6"+")||("+AntiIso_2+"gen_match_2==6"+"))";
  cutsF[3] =weight+"0.5*(("+AntiIso_1+"gen_match_1==6"+")||("+AntiIso_2+"gen_match_2==6"+"))";
  cutsT[1] =weight+"0.5*0.89*(("+AntiIso_1+"gen_match_1!=6"+")||("+AntiIso_2+"gen_match_2!=6"+"))";
  cutsT[2] =weight+"0.5*0.89*(("+AntiIso_1+"gen_match_1!=6"+")||("+AntiIso_2+"gen_match_2!=6"+"))";
  cutsT[3] =weight+"0.5*0.89*(("+AntiIso_1+"gen_match_1!=6"+")||("+AntiIso_2+"gen_match_2!=6"+"))";
  TH1D * histFakes[nSamples];
  TH1D * histTrues[nSamples];
  TH1D * histData[nSamples];

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
  for (int i=0; i<nSamples; ++i) {
        TFile * file = new TFile( directory + samples[i]);
	TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents"); //inputEventsH
	TTree * tree = (TTree*)file->Get("TauCheck");
	double norm = xsec[i]*lumi/nWeightedEvents->GetSumOfWeights();
	cout<<"norm = "<<norm<<endl;
	TString histName   = samples[i] + Variable;
	TString hdata = histName+"_data";
	TString histfake = histName+"_fake";
	TString histtrue = histName+"_true";
	histData[i]    = new TH1D(hdata,"",nbins,bins);
	histFakes[i]   = new TH1D(histfake,"",nbins,bins);
	histTrues[i]   = new TH1D(histtrue,"",nbins,bins);
	if(i==0){
	  tree->Draw(Variable+">>"+hdata,cutsF[i]);
	}
	else{
	  tree->Draw(Variable+">>"+histfake,cutsF[i]);
	  tree->Draw(Variable+">>"+histtrue,cutsT[i]);
	  histFakes[i]->Scale(norm);
	  histTrues[i]->Scale(norm);
	}

  }

  TFile* check = new TFile("check.root","recreate");
  histData[0]->Write();
  histFakes[1]->Write();
  histFakes[2]->Write();
  histFakes[3]->Write();
  check->Close();
  //Space for Stiching

	cout << "data" << " : Entries = " << histData[0]->GetEntries() << " : Integral = " << histData[0]->Integral() << endl;
	cout << samples[1] << " : Entries = " << histFakes[1]->GetEntries() << " : Integral = " << histFakes[1]->Integral() << endl;
	cout << samples[2] << " : Entries = " << histFakes[2]->GetEntries() << " : Integral = " << histFakes[2]->Integral() << endl;
	cout << samples[3] << " : Entries = " << histFakes[3]->GetEntries() << " : Integral = " << histFakes[3]->Integral() << endl;
	cout << samples[1] << " : Entries = " << histTrues[1]->GetEntries() << " : Integral = " << histTrues[1]->Integral() << endl;
	cout << samples[2] << " : Entries = " << histTrues[2]->GetEntries() << " : Integral = " << histTrues[2]->Integral() << endl;
	cout << samples[3] << " : Entries = " << histTrues[3]->GetEntries() << " : Integral = " << histTrues[3]->Integral() << endl;
  
  //W+jet category
  TH1D* WCat_Fake = (TH1D*)histFakes[2]->Clone("WCat_Fake");
  TH1D* WCat_True = (TH1D*)histTrues[2]->Clone("WCat_True");
  WCat_Fake->Add(histFakes[1]);
  WCat_True->Add(histTrues[1]);
  TH1D* WCat_Full = (TH1D*)WCat_Fake->Clone("WCat_Full");
  WCat_Full->Add(WCat_True);
 
  //tt+ll category
  TH1D* TTCat_Fake = (TH1D*)histFakes[3]->Clone("TTCat_Fake");
  TH1D* TTCat_True = (TH1D*)histTrues[3]->Clone("TTCat_True");
  // TTCat_Fake->Add(histFakes[1]);
  // TTCat_True->Add(histTrues[1]);
  TH1D* TTCat_Full = (TH1D*)TTCat_Fake->Clone("TTCat_Full");
  TTCat_Full->Add(TTCat_True);

  //Calculating QCD
  TH1D* QCD = (TH1D*)histData[0]->Clone("QCD");
  QCD->Add(WCat_Full,-1);
  QCD->Add(TTCat_Full,-1);

  cout << "QCD" << " : Entries = " << QCD->GetEntries() << " : Integral = " << QCD->Integral() << endl;

  //Estimation of Background fraction
  TH1D* ff_W = new TH1D("ff_W","",nbins,0,1);
  TH1D* ff_TT = new TH1D("ff_TT","",nbins,0,1);
  TH1D* ff_QCD = new TH1D("ff_QCD","",nbins,0,1);
  for(int ib=0;ib<nbins;ib++){
    double nume_w = -1;
    nume_w = WCat_Fake->GetBinContent(ib);
    double nume_tt = -1;
    nume_tt= TTCat_Fake->GetBinContent(ib);
    double nume_qcd = -1;
    nume_qcd = QCD->GetBinContent(ib);
    double deno = -1;
    deno = (WCat_Fake->GetBinContent(ib)+TTCat_Fake->GetBinContent(ib)+QCD->GetBinContent(ib));
    if(nume_qcd<0)continue;
    double ff_w = nume_w/deno;
    double ff_tt = nume_tt/deno;
    double ff_qcd = nume_qcd/deno;
    ff_W->Fill(ff_w);
    ff_TT->Fill(ff_tt);
    ff_QCD->Fill(ff_qcd);
    cout<<ff_w+ff_tt+ff_qcd<<endl;
  }
  ff_QCD->SetMarkerStyle(21);
  ff_W->SetMarkerStyle(22);
  ff_TT->SetMarkerStyle(23);
  ff_QCD->SetMarkerColor(kRed);
  ff_W->SetMarkerColor(kGreen);
  ff_TT->SetMarkerColor(kBlue);
  //ff_QCD->Draw("p");
  //ff_TT->Draw("psame");
  //ff_W->Draw("psame");
}
