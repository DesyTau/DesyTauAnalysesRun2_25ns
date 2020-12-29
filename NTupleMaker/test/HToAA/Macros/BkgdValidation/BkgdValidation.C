#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"


//////////////////////**** SideBands Names ****////////////////////////////////////////////////////
////// SemiIso ---- LooseIso ---- LooseSemiIso ---- LeadingSemiIso ---- LeadingLooseIso ---- Sel //
///////////////////////////////////////////////////////////////////////////////////////////////////


void BkgdValidation(TString mass="15", TString sideband_1 = "LooseIso_N23", TString sideband_2 = "LooseSemiIso", bool logY = true) {
  
  double massD = atof( mass );

  std::cout << std::endl;
  std::cout << "Mass = " << mass << std::endl;
  std::cout << std::endl;
  
  TString legnedstr_1 = "Sel";
  TString legnedstr_2 = "Semi-Iso";
  if(sideband_1 == "LooseSemiIso") legnedstr_1 = "Loose-Semi-Iso";
  if(sideband_1 == "LooseIso") legnedstr_1 = "Loose-Iso";
  if(sideband_1 == "LeadingSemiIso") legnedstr_1 = "Leading-Semi-Iso";
  if(sideband_1 == "LeadingLooseIso") legnedstr_1 = "Leading-Loose-Iso";
  if(sideband_1 == "MixingMCandDATA") legnedstr_1 = "MC&DATA";
  	
   
  TFile * file = new TFile("../DoubleMuon_BDTOutput_M-"+mass+".root");
  
  
  TH1D * Sideband1H_Old = (TH1D*)file->Get("MVABDTOutput"+sideband_1+"H");
  TH1D * Sideband2H_Old = (TH1D*)file->Get("MVABDTOutput"+sideband_2+"H");
  
  
  Sideband1H_Old->Rebin(100);
  Sideband2H_Old->Rebin(100);
 

  int nBins = Sideband1H_Old->GetNbinsX();
  float iniBin =-1.;
  float endBin = 1.;
  float bins[nBins] ; bins[0]=iniBin;
  for(int i=1;i<=nBins;i++) bins[i] = bins[0]+i*(endBin-iniBin)/nBins;
  
  for(int iBins=1;iBins<Sideband1H_Old->GetNbinsX();iBins++)
  {
   if(Sideband1H_Old->GetBinContent(iBins)==0. || Sideband2H_Old->GetBinContent(iBins)==0.)
   {
    float temp=bins[iBins-(Sideband1H_Old->GetNbinsX()-nBins)];
    for(int jArr=iBins-(Sideband1H_Old->GetNbinsX()-nBins);jArr<nBins;jArr++)
    {
     bins[jArr]=bins[jArr+1];
    }
    bins[nBins]=temp;
    nBins=nBins-1;
   }
  }
   if(Sideband1H_Old->GetBinContent(Sideband1H_Old->GetNbinsX())==0. || Sideband2H_Old->GetBinContent(Sideband2H_Old->GetNbinsX())==0.) 
   {
	   bins[nBins-1]=bins[nBins];
	   nBins=nBins-1;
   }
  
  
  TH1D * Sideband1H = (TH1D*)TH1DtoTH1D(Sideband1H_Old,nBins,bins,true,"_new");
  TH1D * Sideband2H = (TH1D*)TH1DtoTH1D(Sideband2H_Old,nBins,bins,true,"_new");
  
  Sideband2H->Scale(Sideband1H->GetSumOfWeights()/Sideband2H->GetSumOfWeights());
  
  int tcolor = kGreen+4;
  
  ////// Histos Colors ////////
  int tcolor_2 = 4;
  if(sideband_1 == "LooseSemiIso") tcolor_2 = kRed+3;
  if(sideband_1 == "LooseIso") tcolor_2 = kViolet+2;
  if(sideband_1 == "LeadingSemiIso") tcolor_2 = kAzure+7;
  if(sideband_1 == "LeadingLooseIso") tcolor_2 = kRed-9;
  if(sideband_1 == "MixingMCandDATA") tcolor_2 = kPink;
  
  TGraphAsymmErrors * data = new TGraphAsymmErrors();
  data->SetMarkerStyle(kCircle);
  data->SetLineColor(tcolor_2);
  data->SetLineWidth(2);
  data->SetMarkerColor(tcolor_2);
  data->SetMarkerSize(3.5);
  
  TGraphAsymmErrors * ratio = new TGraphAsymmErrors();
  ratio->SetMarkerStyle(kCircle);
  ratio->SetLineColor(tcolor_2);
  ratio->SetLineWidth(2);
  ratio->SetMarkerColor(tcolor_2);
  ratio->SetMarkerSize(3.5);
  
  TGraphAsymmErrors * ratioERR = new TGraphAsymmErrors();
  ratioERR->SetFillColorAlpha(tcolor, 0.55);


  for (int iB=1; iB<=nBins; iB++) 
  {
    
    double y2 = Sideband2H->GetBinContent(iB);
    double y2_err = Sideband2H->GetBinError(iB);
    
    double x1 = Sideband1H->GetBinCenter(iB);
    double x1_err =Sideband1H->GetBinWidth(iB)/2.;
    double y1 = Sideband1H->GetBinContent(iB);
    double y1_errlow = -0.5 + TMath::Sqrt(y1+0.25);
    double y1_errhigh = 0.5 + TMath::Sqrt(y1+0.25);
    
    double r;
    double r_errlow;
    double r_errhigh;
    
    double u;
    double u_err;
    
    if(y2<0.01) 
    {
		r=0.;
		r_errlow=0.;
		r_errhigh=0.;
		
		u=1.;
		u_err=0.;
    }
	else
	{
		r = y1/y2;
        r_errlow = y1_errlow/y2;
        r_errhigh = y1_errlow/y2;
        
        u=1.;
        u_err = y2_err/y2;
	}
    
    data->SetPoint(iB-1,x1,y1);
    data->SetPointError(iB-1,x1_err,x1_err,y1_errlow,y1_errhigh);
    
    ratio->SetPoint(iB-1,x1,r);
    ratio->SetPointError(iB-1,x1_err,x1_err,r_errlow,r_errhigh);
    
    ratioERR->SetPoint(iB-1,x1,u);
    ratioERR->SetPointError(iB-1,x1_err,x1_err,u_err,u_err);
    
    Sideband2H->SetBinError(iB,0);
 
  }
  
  
    
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  
  upper->cd();

  Sideband2H->Draw();
  data->Draw("pZ");
  
  upper->Modified();
  upper->SetTicks();
  upper->SetLeftMargin(0.13);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  
  Sideband2H->SetStats(0);
  Sideband2H->GetXaxis()->SetRangeUser(-1.,1.); 
  Sideband2H->GetXaxis()->SetLabelSize(0.);
  Sideband2H->GetXaxis()->SetLabelOffset(0.);
  Sideband2H->GetYaxis()->SetTitleOffset(1.2);
  Sideband2H->GetYaxis()->SetTitle("Events / Bin");
  Sideband2H->GetYaxis()->SetTitleSize(0.045);
  Sideband2H->GetYaxis()->SetLabelSize(0.04);
  Sideband2H->GetYaxis()->SetTickLength(0.055);
  Sideband2H->GetYaxis()->SetTickSize(0.013);
  Sideband2H->GetYaxis()->SetRangeUser(0,3*Sideband2H->GetMaximum());
  if (logY) Sideband2H->GetYaxis()->SetRangeUser(0.01,10000*Sideband2H->GetMaximum());
  Sideband2H->SetLineColor(tcolor);
  Sideband2H->SetLineWidth(2);
  Sideband2H->SetLineStyle(1);
  
  TLegend * leg = new TLegend(0.15,0.65,0.65,0.85);
  leg->SetTextSize(0.035);
  leg->SetHeader("                  #color[436]{  m_{a_{1}} = "+mass+" GeV }");
  leg->SetLineColor(kOrange+3);
  leg->SetLineWidth(2);
  leg->AddEntry(data,"Sideband Region (#color[635]{"+legnedstr_1+"})","lp");
  leg->AddEntry(Sideband2H,"Sideband Region (#color[635]{"+legnedstr_2+"})","lf");
  leg->Draw();
 
  
  writeExtraText = true;
  extraText = "Work in Progress";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);

  upper->RedrawAxis();
  upper->Update();


  lower->cd();
  

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineColorAlpha(tcolor,0.45);
  
  TMultiGraph *mg = new TMultiGraph();

  mg->Add(ratio,"PZ");
  mg->Add(ratioERR,"2");
  
  mg->Draw("APL");
  line->Draw("same");
 
  lower->Modified();
  lower->SetTicks();
  lower->SetLeftMargin(0.13);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  
  mg->GetXaxis()->SetRangeUser(-1.,1.); 
  mg->GetXaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetLabelOffset(0.03);
  mg->GetXaxis()->SetLabelSize(0.14);
  mg->GetXaxis()->SetTitleSize(0.13);
  mg->GetXaxis()->SetTitleOffset(1.35);
  mg->GetXaxis()->SetTitle("BDT Output");
  mg->GetXaxis()->SetTickLength(0.025);
  mg->GetXaxis()->SetTickSize(0.08);
  mg->GetYaxis()->SetRangeUser(-1.,2.99);
  mg->GetYaxis()->SetLabelOffset(0.008);
  mg->GetYaxis()->SetLabelSize(0.08);
  mg->GetYaxis()->SetTitleSize(0.1);
  mg->GetYaxis()->SetNdivisions(6);
  mg->GetYaxis()->SetTitleOffset(0.45);
  mg->GetYaxis()->SetTitle("Ratio      ");
  mg->GetYaxis()->SetTickLength(0.025);
  mg->GetYaxis()->SetTickSize(0.02);
  
  
  lower->RedrawAxis();
  lower->Update();

}


void PlotValidation(TString mass="15", TString sideband_1 = "MixingMCandDATA", TString sideband_2 = "Sel", bool logY = true){

   gROOT->SetBatch();
   
  
   // Drawing in Canvas 
   TCanvas *c = new TCanvas("c","c",3000,3000) ;
   c->cd();
   
   
   BkgdValidation(mass,sideband_1,sideband_2,logY);

   
   
   c->Update();
   c->Print("Validation_"+sideband_1+"_vs_"+sideband_2+"_ma"+mass+".pdf","Portrait pdf");
   delete c;

}

void PlotAllValidations(){
   
  
  int ngenmass = 15;
  TString genmassString[] = {"4","5","6","7","8","9","10","11","12","13","14","15","17","19","21"};
  
  int nRegions = 7;
  
  TString regions[7] = {"SemiIso", "LooseIso", "LooseSemiIso", "LeadingSemiIso", "LeadingLooseIso", "MixingMCandDATA", "Sel"};

  //********************************************************//
  //**********Loop over the generated mass points***********//
  //********************************************************//

  for (int igenMass=0;igenMass<ngenmass;igenMass++) // Loop over the generated mass points
  {
	 for (int iregion=0;iregion<nRegions;iregion++) PlotValidation(genmassString[igenMass], regions[iregion], "SemiIso", true);
	 
  }
   
   
}

