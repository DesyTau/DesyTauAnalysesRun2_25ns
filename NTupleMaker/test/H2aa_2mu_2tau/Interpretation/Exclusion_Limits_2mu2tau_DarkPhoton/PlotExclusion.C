#include <string>
#include <cstdlib>
#include <iostream>
#include "CMS_lumi.C"

void PlotExclusion()
{

	gROOT->SetBatch();	//to keep canva from popping up

	/////////////////////////////////////////////////////////////////////////
	/////******Loading H->aa->2mu2tau limits as funtion of mass ******//////
	/////////////////////////////////////////////////////////////////////////

	TFile *file_Limits = new TFile("../../Run2Combination/Limits_tree.root");
	TNtuple *Limits = (TNtuple*) file_Limits->Get("Limits");

	int nmass = Limits->GetEntries();

	Float_t massf, minus2f, minus1f, medianf, plus1f, plus2f, obsf;

	Limits->SetBranchAddress("m_a", &massf);
	Limits->SetBranchAddress("minus2", &minus2f);
	Limits->SetBranchAddress("minus1", &minus1f);
	Limits->SetBranchAddress("median", &medianf);
	Limits->SetBranchAddress("plus1", &plus1f);
	Limits->SetBranchAddress("plus2", &plus2f);
	Limits->SetBranchAddress("obs", &obsf);

	TGraph *gr_minus1 = new TGraph(nmass);
	TGraph *gr_median = new TGraph(nmass);
	TGraph *gr_plus1 = new TGraph(nmass);
	TGraph *gr_minus2 = new TGraph(nmass);
	TGraph *gr_obs = new TGraph(nmass);
	TGraph *gr_plus2 = new TGraph(nmass);

	float lowest_ma = 100.;
	float highest_ma = 0.;

	for (Int_t i = 0; i < nmass; i++)
	{
		Limits->GetEntry(i);

		gr_minus2->SetPoint(i, massf, minus2f / 1);
		gr_minus1->SetPoint(i, massf, minus1f / 1);
		gr_median->SetPoint(i, massf, medianf / 1);
		gr_plus1->SetPoint(i, massf, plus1f / 1);
		gr_plus2->SetPoint(i, massf, plus2f / 1);
		gr_obs->SetPoint(i, massf, obsf / 1);

		if (massf < lowest_ma) lowest_ma = massf;
		if (massf > highest_ma) highest_ma = massf;
	}

    
	cout << "Lowest Mass: " << lowest_ma << " GeV" << endl;
	cout << "Highest Mass: " << highest_ma << " GeV" << endl;

	/////******Defining TGraph to be used ******//////
	TGraph *Median_G = new TGraph();
	Median_G->SetLineColor(1);
	Median_G->SetLineWidth(1);
	Median_G->SetMarkerColor(2);
	Median_G->SetFillColor(kWhite);
	Median_G->SetMarkerSize(1);
	Median_G->SetMarkerStyle(9);

	TGraph *Obs_G = new TGraph();
	Obs_G->SetLineColor(kBlue);
	Obs_G->SetLineWidth(1);
	Obs_G->SetFillColorAlpha(kCyan, 0.35);

	TGraphAsymmErrors *Inner_Band = new TGraphAsymmErrors();
	Inner_Band->SetFillColorAlpha(kGreen + 1, 0.75);
	Inner_Band->SetLineColor(kWhite);

	TGraphAsymmErrors *Outer_Band = new TGraphAsymmErrors();
	Outer_Band->SetFillColorAlpha(kOrange, 0.85);
	Outer_Band->SetLineColor(kWhite);

	TGraphAsymmErrors *Shaded_Zone = new TGraphAsymmErrors();
	Shaded_Zone->SetFillColorAlpha(kCyan, 0.35);

	int nbins = 4;
	float bins[5] = { 3.6, 4.0, 9., 11., 21 };
	TH1F *Uncertain_Regions = new TH1F("Uncertain_Regions", "Uncertain_Regions", nbins, bins);
	Uncertain_Regions->SetBinContent(1, 2);
	Uncertain_Regions->SetBinError(1, 100);
	Uncertain_Regions->SetBinContent(3, 2);
	Uncertain_Regions->SetBinError(3, 100);
	Uncertain_Regions->SetFillColorAlpha(kGray, 0.5);
	Uncertain_Regions->SetLineColorAlpha(kGray, 0.5);

	/////////////////////////////////////////////////////////////////////////////////
	/////******Loading NTuple with 1D Distribution of BR (m_a) ******//////
	/////////////////////////////////////////////////////////////////////////////////

	float m_a, BR_a_ll;

	TFile *file = new TFile("BR.root");
	TNtuple *ntuple = (TNtuple*) file->Get("Ntuple_DarkPhoton");
	ntuple->SetBranchAddress("m_a", &m_a);
	ntuple->SetBranchAddress("BR_a_ll", &BR_a_ll);

    TGraph *BR_a_ll_Graph = new TGraph();
    int ipoint = -1;
	for (Int_t i = 0; i < ntuple->GetEntries(); i++)
	{
		ntuple->GetEntry(i);
		
		ipoint++;

		BR_a_ll_Graph->SetPoint(i, m_a, BR_a_ll);
	}
	
	int ipoints = -1;
	for (Int_t i = 0; i <= 200; i++)
	{
		float masspoint = 3.6+(21-3.6)/200*i;
		
		Float_t X, Y_Minus1, Y_Median, Y_Minus2, Y_Plus1, Y_Plus2, Y_Obs;

		float BR_H_aa_2mu2tau = 2 *BR_a_ll_Graph->Eval(masspoint, 0, "S") * BR_a_ll_Graph->Eval(masspoint, 0, "S");	//BR_aa floating 

		ipoints++;

        //cout << "ipoint:" << ipoint << " m_a:" << m_a << " BR_a_ll:"<< BR_a_ll<< endl;
		Y_Minus2 = gr_minus2->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);
		Y_Minus1 = gr_minus1->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);
		Y_Median = gr_median->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);
		Y_Plus1 = gr_plus1->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);
		Y_Plus2 = gr_plus2->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);
		Y_Obs = gr_obs->Eval(masspoint, 0, "S") / (BR_H_aa_2mu2tau *1000);

		Median_G->SetPoint(ipoints, masspoint, Y_Median);
		Obs_G->SetPoint(ipoints, masspoint, Y_Obs);

		Inner_Band->SetPoint(ipoints, masspoint, Y_Median);
		Inner_Band->SetPointError(ipoints, 0., 0., fabs(Y_Minus1 - Y_Median), fabs(Y_Plus1 - Y_Median));

		Outer_Band->SetPoint(ipoints, masspoint, Y_Median);
		Outer_Band->SetPointError(ipoints, 0., 0., fabs(Y_Minus2 - Y_Median), fabs(Y_Plus2 - Y_Median));

		Shaded_Zone->SetPoint(ipoints, masspoint, Y_Obs);
		Shaded_Zone->SetPointError(ipoints, 0., 0., 0., 100);
	}

	////////////////////////////////////////////////////////
	/////******Ploting Constraints in 1D plot ******//////
	////////////////////////////////////////////////////////

	TCanvas *c = new TCanvas("c", "c", 900, 900);
	c->SetLogy();

	TLine *line = new TLine(lowest_ma, 1., highest_ma - 1, 1.);
	line->SetLineStyle(9);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(Outer_Band, "3");
	mg->Add(Inner_Band, "3");
    mg->Add(Shaded_Zone, "3");
	mg->Add(Obs_G, "L");
	mg->Add(Median_G, "L");

	mg->Draw("APL");
	line->Draw("same");
	Uncertain_Regions->Draw("sameh2");

	TLegend *legend = new TLegend(0.58, 0.15, 0.87, 0.4);
	legend->SetTextSize(0.03);
	legend->AddEntry(Outer_Band, "#pm 2#sigma expected", "f");
	legend->AddEntry(Inner_Band, "#pm 1#sigma expected", "f");
	legend->AddEntry(Median_G, "Expected", "l");
	legend->AddEntry(Obs_G, "Observed", "lf");
	legend->SetHeader("Dark photon", "C");	// option "C" allows to center the header
	TLegendEntry *header = (TLegendEntry*) legend->GetListOfPrimitives()->First();
	header->SetTextSize(.04);
	header->SetTextFont(22);
	legend->Draw();

	gPad->Modified();
	gPad->SetTicks();
	mg->SetMinimum(0.000001);
	mg->SetMaximum(1.5);
	mg->GetXaxis()->SetTitle("m_{Z_{D}}[GeV]");
	mg->GetXaxis()->SetLimits(lowest_ma, highest_ma - 1);
	mg->GetXaxis()->SetLabelSize(0.03);
	mg->GetYaxis()->SetTitle("95% CL on  #font[42]{#frac{#sigma(h)}{#sigma_{SM}} BR(#font[42]{h#rightarrow Z_{D}Z_{D}})}  ");
	mg->GetYaxis()->SetTitleFont(132);
	mg->GetXaxis()->SetTitleFont(132);
	mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetXaxis()->SetTitleSize(0.045);
	mg->GetYaxis()->SetTitleOffset(1.6);
	mg->GetYaxis()->SetLabelSize(0.03);

	gPad->RedrawAxis();
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.1);

	extraText = "Private work";
	writeExtraText = false;
	CMS_lumi(c, 12, 33);
	c->Update();
	c->Print("Exclusion_Limits.pdf");
}
