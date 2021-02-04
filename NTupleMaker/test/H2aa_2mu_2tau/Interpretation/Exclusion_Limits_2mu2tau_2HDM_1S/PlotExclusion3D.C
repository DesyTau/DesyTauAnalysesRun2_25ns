#include <string>
#include "CMS_lumi.C"

void PlotExclusion3D_type(TString Model_Type = "IV", TString BR_h_aa_Upper = "0.34")
{

	const int MassEntries = 499,
	TanBetaEntries = 24;

	float BR_H_aa = atof(BR_h_aa_Upper);

	//////////////////////////////////////////////////////////////////////////
	//////******Loading H->aa->2mu2tau limits as funtion of mass ******//////
	//////////////////////////////////////////////////////////////////////////

	TFile *file_Limits = new TFile("../../Run2Combination/Limits_tree.root");
	TNtuple *Limits = (TNtuple*) file_Limits->Get("Limits");

	int nmass = Limits->GetEntries();

	Float_t massf, obsf;

	Limits->SetBranchAddress("m_a", &massf);
	Limits->SetBranchAddress("obs", &obsf);
	//Limits->SetBranchAddress("median", &obsf);

	TGraph *gr_obs = new TGraph(nmass);

	float lowest_ma = 100.;
	float highest_ma = 0.;

	for (Int_t i = 0; i < nmass; i++)
	{
		Limits->GetEntry(i);
		gr_obs->SetPoint(i, massf, obsf / 1000);

		if (massf < lowest_ma) lowest_ma = massf;
		if (massf > highest_ma) highest_ma = massf;
	}

	/////////////////////////////////////////////////////////////////////////////////
	/////******Loading NTuple with 2D Distribution of BR (m_a,tan_beta) ******//////
	/////////////////////////////////////////////////////////////////////////////////

	float m_a, tan_beta, BR_a_tautau, BR_a_mumu;

	TFile *file = new TFile("BR.root");
	TNtuple *ntuple = (TNtuple*) file->Get("Ntuple_" + Model_Type);
	ntuple->SetBranchAddress("m_a", &m_a);
	ntuple->SetBranchAddress("tan_beta", &tan_beta);
	ntuple->SetBranchAddress("BR_a_tautau", &BR_a_tautau);
	ntuple->SetBranchAddress("BR_a_mumu", &BR_a_mumu);

	////// values of m_a and tan_beta	///////
	std::vector<Double_t> *tmp1;
	file->GetObject("mass_values", tmp1);
	std::vector<Double_t> m_a_values = *tmp1;

	std::vector<Double_t> *tmp2;
	file->GetObject("tanbeta_values", tmp2);
	std::vector<Double_t> tan_beta_values = *tmp2;

	/////// Creating TH2D	////////////
	TH2D * BRvsMaTanbeta_2DH;
	float tan_beta_bins[TanBetaEntries + 1];
	float mass_bins[MassEntries + 1];

	tan_beta_bins[0] = 0.45;
	tan_beta_bins[TanBetaEntries] = 10.25;
	mass_bins[0] = 0.97;
	mass_bins[MassEntries] = 70.15;

	for (Int_t i = 0; i < MassEntries - 1; i++) mass_bins[i + 1] = float(0.5 *(m_a_values[i] + m_a_values[i + 1]));
	for (Int_t i = 0; i < TanBetaEntries - 1; i++) tan_beta_bins[i + 1] = float(0.5 *(tan_beta_values[i] + tan_beta_values[i + 1]));

	BRvsMaTanbeta_2DH = new TH2D("BRvsMaTanbeta_2DH", "", MassEntries, mass_bins, TanBetaEntries, tan_beta_bins);

	/////// Setting limits on BR(h->aa)	//////
	for (Int_t i = 0; i < ntuple->GetEntries(); i++)
	{

		if (Model_Type == "I")
		{
			for (Int_t j = 0; j < TanBetaEntries; j++)
			{
				ntuple->GetEntry(i);

				Int_t mass_bin = BRvsMaTanbeta_2DH->GetXaxis()->FindBin(m_a);
				Int_t tanbeta_bin = BRvsMaTanbeta_2DH->GetYaxis()->FindBin(tan_beta_values[j]);

				float BR_aa_2mu2tau = 2 *BR_a_mumu * BR_a_tautau;

				if (m_a > lowest_ma && m_a < highest_ma) BRvsMaTanbeta_2DH->SetBinContent(mass_bin, tanbeta_bin, gr_obs->Eval(m_a, 0, "S") / BR_aa_2mu2tau);
			}
		}
		else
		{
			ntuple->GetEntry(i);

			Int_t mass_bin = BRvsMaTanbeta_2DH->GetXaxis()->FindBin(m_a);
			Int_t tanbeta_bin = BRvsMaTanbeta_2DH->GetYaxis()->FindBin(tan_beta);

			float BR_aa_2mu2tau = 2 *BR_a_mumu * BR_a_tautau;

			if (m_a > lowest_ma && m_a < highest_ma) BRvsMaTanbeta_2DH->SetBinContent(mass_bin, tanbeta_bin, gr_obs->Eval(m_a, 0, "S") / BR_aa_2mu2tau);
		}
	}

	//////////////////////////////////////////////////////////
	/////******Ploting Constraints in 2D+1 plane ******//////
	//////////////////////////////////////////////////////////

	BRvsMaTanbeta_2DH->SetStats(0);

	BRvsMaTanbeta_2DH->GetXaxis()->SetTitle("m_{a_{1}}[GeV]");
	BRvsMaTanbeta_2DH->GetXaxis()->SetRangeUser(lowest_ma + 0.05, highest_ma - 0.2);
	BRvsMaTanbeta_2DH->GetXaxis()->SetLabelSize(0.03);
	BRvsMaTanbeta_2DH->GetXaxis()->SetTitleFont(132);
	BRvsMaTanbeta_2DH->GetXaxis()->SetTitleSize(0.05);
	BRvsMaTanbeta_2DH->GetXaxis()->SetTitleOffset(0.9);

	BRvsMaTanbeta_2DH->GetYaxis()->SetTitle("tan#beta");
	BRvsMaTanbeta_2DH->GetYaxis()->SetRangeUser(0.5 + 0.05, 10 - 0.05);
	BRvsMaTanbeta_2DH->GetYaxis()->SetLabelSize(0.03);
	BRvsMaTanbeta_2DH->GetYaxis()->SetTitleFont(132);
	BRvsMaTanbeta_2DH->GetYaxis()->SetTitleSize(0.055);
	BRvsMaTanbeta_2DH->GetYaxis()->SetTitleOffset(0.8);

	BRvsMaTanbeta_2DH->GetZaxis()->SetTitle("#scale[0.6]{95% CL on #font[42]{#frac{#sigma(h)}{#sigma_{SM}} BR(#font[42]{h#rightarrowa_{1}a_{1}})}}");
	BRvsMaTanbeta_2DH->GetZaxis()->SetRangeUser(0.001, 10);
	BRvsMaTanbeta_2DH->GetZaxis()->SetLabelSize(0.03);
	BRvsMaTanbeta_2DH->GetZaxis()->SetTitleFont(132);
	BRvsMaTanbeta_2DH->GetZaxis()->SetTitleSize(0.065);
	BRvsMaTanbeta_2DH->GetZaxis()->SetTitleOffset(1.);

	BRvsMaTanbeta_2DH->DrawCopy("colz");

	TH2D *Contour1_2DH = (TH2D*) BRvsMaTanbeta_2DH->Clone("Contour1_2DH");
	TH2D *Contour2_2DH = (TH2D*) BRvsMaTanbeta_2DH->Clone("Contour2_2DH");

	double contour_1[1];
	contour_1[0] = BR_H_aa;

	double contour_2[1];
	contour_2[0] = 0.05;

	Contour1_2DH->SetContour(1, contour_1);
	Contour1_2DH->SetLineColor(kRed);
	Contour1_2DH->SetLineStyle(kDashed);
	Contour1_2DH->SetLineWidth(2);
	Contour1_2DH->Draw("cont3 same");

	Contour2_2DH->SetContour(1, contour_2);
	Contour2_2DH->SetLineColor(kBlue + 3);
	Contour2_2DH->SetLineStyle(kDashed);
	Contour2_2DH->SetLineWidth(2);
	Contour2_2DH->Draw("cont3 same");

	TLegend *legend = new TLegend(0.15, 0.68, 0.53, 0.86);
	//legend->SetFillStyle(0);
	//legend->SetBorderSize(0);
	legend->SetLineColor(kWhite);
	legend->SetTextSize(0.04);
	legend->AddEntry(Contour1_2DH, "#scale[0.55]{#font[42]{95% CL on #frac{#sigma}{#sigma_{SM}} BR(#font[42]{h#rightarrowa_{1}a_{1}}) = " + BR_h_aa_Upper + "}}", "L");
	legend->AddEntry(Contour2_2DH, "#scale[0.55]{#font[42]{95% CL on #frac{#sigma}{#sigma_{SM}} BR(#font[42]{h#rightarrowa_{1}a_{1}}) = 0.05}}", "L");
	legend->SetHeader("#color[1]{2HDM+S Type " + Model_Type + "}", "C");	// option "C" allows to center the header
	legend->SetMargin(0.15);
	legend->SetEntrySeparation(0.5);
	TLegendEntry *header = (TLegendEntry*) legend->GetListOfPrimitives()->First();
	header->SetTextSize(.04);
	header->SetTextFont(22);
	legend->Draw();

	gStyle->SetPalette(kViridis);
	gStyle->SetNumberContours(500);

	gPad->Modified();
	gPad->SetLogz();
	gPad->SetTicks();
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.2);

	gPad->RedrawAxis();

}

void PlotExclusion3D()
{

	gROOT->SetBatch();	//to keep canva from popping up
	int nType = 4;
	TString Model_Type[] = { "I", "II", "III", "IV" };

	// Drawing in Canvas 
	TCanvas *c = new TCanvas("c", "c", 600, 600);

	for (int i = 0; i < nType; i++)
	{
		PlotExclusion3D_type(Model_Type[i], "0.34");
		extraText = "Private work";
		writeExtraText = false;
		CMS_lumi(c, 12, 33);
		c->Print("Exclusion_Limits3D_" + Model_Type[i] + ".pdf");
	}

	c->Update();
}
