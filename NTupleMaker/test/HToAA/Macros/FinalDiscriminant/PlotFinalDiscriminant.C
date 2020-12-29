#include "CMS_lumi.C"
#include "HttStylesNew.cc"


void PlotFinalDiscriminant(TString mass="5", bool blinddata = false, bool drawLeg = true, bool logY = true) {
  

  TFile * file  = new TFile("../DataCards/haa-13TeV_2016_ma"+mass+".root");
  TFile * fileFit = new TFile("../DataCards/fitDiagnostics_ma"+mass+".root");
  
  TH1D * histggh  = (TH1D*)file->Get("ggh");
  TH1D * histmmtt  = (TH1D*)file->Get("mmtt");
  TH1D * histvbf  = (TH1D*)file->Get("vbf");
  TH1D * histvh  = (TH1D*)file->Get("vh");
  TH1D * histtth  = (TH1D*)file->Get("tth");
  TH1D * histData = (TH1D*)file->Get("data_obs");
  TH1D * histBkgd = (TH1D*)file->Get("bkgd");
  TH1D * histBkgdX = (TH1D*)fileFit->Get("shapes_fit_b/haa/total_background")->Clone("histBkgdX");
  
  //////// Stacking Histograms on top of eachother ///////////
  histvh->Add(histtth);
  histvbf->Add(histvh);
  histggh->Add(histvbf);
  
  histggh->Scale(0.2);
  histmmtt->Scale(0.2);

  ////////////// Printing Yileds ////////////
  std::cout<< std::endl;
  std::cout<< "Yield H->aa->4tau: " << histggh->GetSumOfWeights() << std::endl;
  std::cout<< "Yield H->aa->2tau2mu: " << histmmtt->GetSumOfWeights() << std::endl;

  int nBins=histBkgd->GetNbinsX();

  /////// Observed Data Graph //////////
  TGraphAsymmErrors * data = new TGraphAsymmErrors();
  data->SetMarkerStyle(kCircle);
  data->SetLineColor(4);
  data->SetLineWidth(2);
  data->SetMarkerColor(4);
  data->SetMarkerSize(3.5);
  
  /////// Ratio observed data/bkgd Graph //////////
  TGraphAsymmErrors * ratio = new TGraphAsymmErrors();
  ratio->SetMarkerStyle(kCircle);
  ratio->SetLineColor(4);
  ratio->SetLineWidth(2);
  ratio->SetMarkerColor(4);
  ratio->SetMarkerSize(3.5);
 
  /////// Error ratio observed data/bkgd Graph //////////
  TGraphAsymmErrors * ratioERR = new TGraphAsymmErrors();
  ratioERR->SetFillColorAlpha(kOrange-3, 0.35);


  for (int iB=1; iB<=nBins; iB++) 
  {
    
    double xB = histBkgdX->GetBinContent(iB);
    
    double x = histData->GetBinCenter(iB);
    double ex =histData->GetBinWidth(iB)/2.;
    double y = histData->GetBinContent(iB);
    double eylow = -0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    double eyhigh = 0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    
    double r;
    double erlow;
    double erhigh;
    
    double erERR;
    
    if(xB<0.01) 
    {
		r=0.;
		erlow=0.;
		erhigh=0.;
		
		erERR=0.;
    }
	else
	{
		r = histData->GetBinContent(iB)/xB;
        erlow = eylow/xB;
        erhigh = eyhigh/xB;
        
        erERR = histBkgdX->GetBinError(iB)/xB;
	}
	
    if(blinddata && x>-0.2) {y=1e10;r=1000;}
    
    data->SetPoint(iB-1,x,y);
    data->SetPointError(iB-1,ex,ex,eylow,eyhigh);
    
    ratio->SetPoint(iB-1,x,r);
    ratio->SetPointError(iB-1,ex,ex,erlow,erhigh);
    
    ratioERR->SetPoint(iB-1,x,1.);
    ratioERR->SetPointError(iB-1,ex,ex,erERR,erERR);
    
    histBkgd->SetBinContent(iB,xB);
    histBkgd->SetBinError(iB,0);
    histggh->SetBinError(iB,0);
    histvbf->SetBinError(iB,0);
    histvh->SetBinError(iB,0);
    histtth->SetBinError(iB,0);
    histmmtt->SetBinError(iB,0);
    
  }
  
  TH1D * hist4t = (TH1D*)histggh->Clone("hist4t");
  

  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();

  upper->cd();

  histBkgd->Draw("h");
  hist4t->Draw("hsame"); 
  histmmtt->Draw("hsame");  
  data->Draw("pZ");

  upper->Modified();
  upper->SetTicks();
  upper->SetLeftMargin(0.13);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  
  histBkgd->GetXaxis()->SetRangeUser(-1.,1.); 
  histBkgd->GetXaxis()->SetLabelSize(0.);
  histBkgd->GetXaxis()->SetLabelOffset(0.);
  histBkgd->GetYaxis()->SetTitleOffset(1.1);
  histBkgd->GetYaxis()->SetTitle("Events / Bin");
  histBkgd->GetYaxis()->SetTitleSize(0.045);
  histBkgd->GetYaxis()->SetLabelSize(0.04);
  histBkgd->GetYaxis()->SetTickLength(0.055);
  histBkgd->GetYaxis()->SetTickSize(0.013);
  histBkgd->GetYaxis()->SetRangeUser(0,3*histBkgd->GetMaximum());
  if (logY) histBkgd->GetYaxis()->SetRangeUser(0.01,1000*histBkgd->GetMaximum());
  histBkgd->SetLineColor(kBlack);
  histBkgd->SetLineWidth(1);
  histBkgd->SetLineStyle(1);
  histBkgd->SetFillColorAlpha(kOrange-3,0.35);
  histBkgd->SetFillStyle(1001);
  
  TH1D * histBkgdLegend = (TH1D*)histBkgd->Clone("histBkgd");
  InitHist(histBkgdLegend,"","",TColor::GetColor(kBlue-7),1001);
  histBkgdLegend->SetLineColor(kBlack);
  histBkgdLegend->SetLineWidth(1);
  histBkgdLegend->SetLineStyle(1);
  histBkgdLegend->SetFillColorAlpha(kOrange-3, 0.35);

  
  hist4t->SetLineColor(kTeal+4);
  hist4t->SetLineWidth(3);
  hist4t->SetLineStyle(2);
  
  histmmtt->SetLineColor(kRed);
  histmmtt->SetLineWidth(3);
  histmmtt->SetLineStyle(3);
  
  TLegend * leg = new TLegend(0.2,0.60,0.55,0.85);
  leg->SetHeader("#color[436]{  m_{a_{1}} = "+mass+" GeV }"); //, "C");
  leg->SetTextSize(0.035);
  leg->SetLineColor(kOrange+3);
  leg->SetNColumns(2);
  leg->AddEntry(data,"Observed","lep");
  leg->AddEntry(histBkgdLegend,"Bkg","lf");
  leg->AddEntry(hist4t," 4#tau","l");
  leg->AddEntry(histmmtt," 2#mu2#tau","l");
  if (drawLeg) {
    leg->Draw();
  }
  
  
  writeExtraText = true;
  extraText = "Private work";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);

  upper->RedrawAxis();
  upper->Update();


  lower->cd();
  

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineColorAlpha(kBlue,0.65);
  
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
  mg->GetYaxis()->SetRangeUser(0.,2.49);
  mg->GetYaxis()->SetLabelOffset(0.008);
  mg->GetYaxis()->SetLabelSize(0.08);
  mg->GetYaxis()->SetTitleSize(0.1);
  mg->GetYaxis()->SetNdivisions(6);
  mg->GetYaxis()->SetTitleOffset(0.45);
  mg->GetYaxis()->SetTitle("Obs/Bkg   ");
  mg->GetYaxis()->SetTickLength(0.025);
  mg->GetYaxis()->SetTickSize(0.02);
  
  lower->RedrawAxis();
  lower->Update();
  
}

void PlotDiscriminant(TString mass="10", bool blinddata = false, bool drawLeg = true, bool logY = true){

   gROOT->SetBatch();
   
   
   // Drawing in Canvas 
   TCanvas *c = new TCanvas("c","c",3000,3000);
   c->cd();
	  
   PlotFinalDiscriminant(mass, blinddata,drawLeg,logY);
   
   c->Update();
   c->Print("FinalDiscriminant_ma"+mass+".pdf","Portrait pdf");

}

void PlotAllDiscriminants(bool blinddata = false){
   
   
   int ngenmass = 15;
   TString genmassString[] = {"4","5","6","7","8","9","10","11","12","13","14","15","17","19","21"};
  
   for (int igenMass=0;igenMass<ngenmass;igenMass++) PlotDiscriminant(genmassString[igenMass], blinddata, true, true);

}


  

