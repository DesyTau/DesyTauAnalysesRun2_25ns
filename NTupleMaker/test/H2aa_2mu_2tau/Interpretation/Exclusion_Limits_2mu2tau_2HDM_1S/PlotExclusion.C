#include <string>
#include <cstdlib>
#include <iostream>
#include "CMS_lumi.C"

void Exclusion(TString Model_Type = "I", TString Tan_Beta = "--")
{

	const int MassEntries = 499,
	TanBetaEntries = 24;

	float BR_H_aa = 1.0;
	double Tan_Betaf = atof(Tan_Beta);
	int itan_beta = 0;

	if (Tan_Betaf <= 1.0) itan_beta = int(Tan_Betaf / 0.1) - 5;
	if (Tan_Betaf > 1.0) itan_beta = int(Tan_Betaf / 0.5) - 2 + 5;
	if (itan_beta < 0) itan_beta = 0;

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

		gr_minus2->SetPoint(i, massf, minus2f / 1000);
		gr_minus1->SetPoint(i, massf, minus1f / 1000);
		gr_median->SetPoint(i, massf, medianf / 1000);
		gr_plus1->SetPoint(i, massf, plus1f / 1000);
		gr_plus2->SetPoint(i, massf, plus2f / 1000);
		gr_obs->SetPoint(i, massf, obsf / 1000);

		if (massf < lowest_ma) lowest_ma = massf;
		if (massf > highest_ma) highest_ma = massf;
	}

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

	/////////////////////////////////////////////////////////////////////////////////
	/////******Loading NTuple with 1D Distribution of BR (m_a,tan_beta) ******//////
	/////////////////////////////////////////////////////////////////////////////////

	float m_a, tan_beta, BR_a_tautau, BR_a_mumu;

	TFile *file = new TFile("BR.root");
	TNtuple *ntuple = (TNtuple*) file->Get("Ntuple_" + Model_Type);
	ntuple->SetBranchAddress("m_a", &m_a);
	ntuple->SetBranchAddress("tan_beta", &tan_beta);
	ntuple->SetBranchAddress("BR_a_tautau", &BR_a_tautau);
	ntuple->SetBranchAddress("BR_a_mumu", &BR_a_mumu);

	int ipoint = -1;
	for (Int_t i = MassEntries * itan_beta; i < MassEntries *(itan_beta + 1); i++)
	{
		ntuple->GetEntry(i);
		Float_t X, Y_Minus1, Y_Median, Y_Minus2, Y_Plus1, Y_Plus2, Y_Obs;

		float BR_H_aa_2mu2tau = BR_H_aa *2 *BR_a_mumu * BR_a_tautau;

		if (m_a < lowest_ma || m_a > highest_ma) continue;
		ipoint++;

		Y_Minus2 = gr_minus2->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;
		Y_Minus1 = gr_minus1->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;
		Y_Median = gr_median->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;
		Y_Plus1 = gr_plus1->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;
		Y_Plus2 = gr_plus2->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;
		Y_Obs = gr_obs->Eval(m_a, 0, "S") / BR_H_aa_2mu2tau;

		Median_G->SetPoint(ipoint, m_a, Y_Median);
		Obs_G->SetPoint(ipoint, m_a, Y_Obs);

		Inner_Band->SetPoint(ipoint, m_a, Y_Median);
		Inner_Band->SetPointError(ipoint, 0., 0., fabs(Y_Minus1 - Y_Median), fabs(Y_Plus1 - Y_Median));

		Outer_Band->SetPoint(ipoint, m_a, Y_Median);
		Outer_Band->SetPointError(ipoint, 0., 0., fabs(Y_Minus2 - Y_Median), fabs(Y_Plus2 - Y_Median));

		Shaded_Zone->SetPoint(ipoint, m_a, Y_Obs);
		Shaded_Zone->SetPointError(ipoint, 0., 0., 0., 100);
	}

	////////////////////////////////////////////////////////
	/////******Ploting Constraints in 1D plot ******//////
	////////////////////////////////////////////////////////
	TLine *line = new TLine(lowest_ma, 1., highest_ma - 0.2, 1.);
	line->SetLineStyle(9);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(Outer_Band, "3");
	mg->Add(Inner_Band, "3");
	mg->Add(Shaded_Zone, "3");
	mg->Add(Obs_G, "L");
	mg->Add(Median_G, "C");

	mg->Draw("APL");
	line->Draw("same");

	TLegend *legend = new TLegend(0.38, 0.15, 0.87, 0.4);
	legend->SetTextSize(0.03);
	legend->AddEntry(Outer_Band, "#pm 2#sigma expected", "f");
	legend->AddEntry(Inner_Band, "#pm 1#sigma expected", "f");
	legend->AddEntry(Median_G, "Expected", "l");
	legend->AddEntry(Obs_G, "Observed", "lf");
	if (Model_Type != "I") legend->SetHeader("2HDM+S Type " + Model_Type + ", tan#beta = " + Tan_Beta, "C");	// option "C" allows to center the header
	else legend->SetHeader("2HDM+S Type " + Model_Type, "C");
	TLegendEntry *header = (TLegendEntry*) legend->GetListOfPrimitives()->First();
	header->SetTextSize(.04);
	header->SetTextFont(22);
	legend->Draw();

	gPad->Modified();
	gPad->SetTicks();
	mg->SetMinimum(0.0001);
	mg->SetMaximum(1.5);
	mg->GetXaxis()->SetTitle("m_{a_{1}}[GeV]");
	mg->GetXaxis()->SetLimits(lowest_ma, highest_ma - 0.2);
	mg->GetXaxis()->SetLabelSize(0.03);
	mg->GetYaxis()->SetTitle("95% CL on #font[42]{#frac{#sigma(h)}{#sigma_{SM}} BR(#font[42]{h#rightarrowa_{1}a_{1}})}  ");
	mg->GetYaxis()->SetTitleFont(132);
	mg->GetXaxis()->SetTitleFont(132);
	mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetXaxis()->SetTitleSize(0.045);
	mg->GetYaxis()->SetTitleOffset(1.6);
	mg->GetYaxis()->SetLabelSize(0.03);

	gPad->RedrawAxis();
}

void PlotExclusion()
{
	gROOT->SetBatch();	//to keep canva from popping up

	int nType = 4;
	TString Model_Type[] = { "I", "II", "III", "IV" };
	TString Tan_Beta[] = { "", "5", "2", "0.5" };

	// Drawing in Canvas 
	TCanvas *c = new TCanvas("c", "", 900, 900);

	for (int i = 0; i < nType; i++)
	{
		gPad->SetLogy();
		gPad->SetLeftMargin(0.15);
		Exclusion(Model_Type[i], Tan_Beta[i]);
		extraText = "Private work";
		writeExtraText = false;
		CMS_lumi(c, 12, 33);
		c->Print("Exclusion_Limits_" + Model_Type[i] + ".pdf");
	}

	c->Update();

}
