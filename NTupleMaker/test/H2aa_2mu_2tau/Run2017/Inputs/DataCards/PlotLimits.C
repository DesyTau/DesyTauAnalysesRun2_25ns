#include "HttStylesNew.cc"
#include "CMS_lumi.C"

void PlotLimits(bool blindData = true, char *fileList = "limits")
{

	gROOT->SetBatch();
	
	const int nPoints = 88;

	//// signal strength limits sigma*BR / sigma*BR (at tanb=30)	////
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

	while (inputList >> fileName)
	{

		//    std::cout << fileName << std::endl;

		TFile *file = new TFile(fileName);

		TTree *tree = (TTree*) file->Get("limit");

		//    std::cout << "file : " << file << std::endl;
		//    std::cout << "tree : " << tree << std::endl;

		tree->SetBranchAddress("limit", &LIMIT);
		tree->SetBranchAddress("mh", &MH);

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

		counter++;
	}

	std::cout << " m(Phi1)  -2s   -1s   exp   +1s   +2s   obs " << std::endl;
	//           "100  24.1  28.2  33.8  40.8  48.2  23.0

	auto f_out = TFile::Open("Limits_tree.root", "RECREATE");
	TNtuple Limits("Limits", "Limits", "m_a:minus2:minus1:median:plus1:plus2:obs");

	for (int i = 0; i < counter; ++i)
	{

		obs[i] = obsR[i];
		minus2[i] = minus2R[i];
		minus1[i] = minus1R[i];
		median[i] = medianR[i];
		plus1[i] = plus1R[i];
		plus2[i] = plus2R[i];

		Limits.Fill(mA[i], minus2[i], minus1[i], median[i], plus1[i], plus2[i], obs[i]);

		char strOut[200];
		sprintf(strOut, "%5.10f  %5.10f  %5.10f  %5.10f  %5.10f  %5.10f  %5.10f",
			mA[i], minus2[i], minus1[i], median[i], plus1[i], plus2[i], obs[i]);
		std::cout << strOut << std::endl;
		cout << "endl" << endl;
	}

	f_out->Write();

    double zeros[counter];
    
	for (int i = 0; i < counter; ++i)
	{
		zeros[i]=0;
		minus2[i] = median[i] - minus2[i];
		minus1[i] = median[i] - minus1[i];
		plus1[i] = plus1[i] - median[i];
		plus2[i] = plus2[i] - median[i];
	}      

	int nPointsX = counter;

	TGraph *obsG = new TGraph(nPointsX, mA, obs);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(1);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(0);
	obsG->SetMarkerSize(0);

	TGraph *expG = new TGraph(nPointsX, mA, median);
	expG->SetLineWidth(2);
	expG->SetLineColor(2);
	expG->SetLineStyle(2);

	TGraphAsymmErrors *innerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus1, plus1);
	innerBand->SetFillColor(kGreen + 1);
	innerBand->SetLineColor(0);

	TGraphAsymmErrors *outerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus2, plus2);
	outerBand->SetFillColor(kOrange);
	outerBand->SetLineColor(0);

    ////// Plotting /////
    
    TCanvas *canv = MakeCanvas("canv", "histograms", 900, 900);
    canv->SetLeftMargin(0.15);
    canv->SetRightMargin(0.1);
    canv->SetBottomMargin(0.15);
	canv->SetTicks();
	
	TH1F *frame = NULL;

	frame = new TH1F("frame", "", 2, 3.59, 21.01);
	frame->SetStats(0);
	frame->GetYaxis()->SetRangeUser(0.,2.);
	frame->GetXaxis()->SetTitle("m_{a_{1}} [GeV]");
	frame->GetYaxis()->SetTitle("#sigma B / #sigma_{SM}");
	frame->GetXaxis()->SetNdivisions(510);
	frame->GetYaxis()->SetNdivisions(510);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleOffset(1.25);
	frame->GetXaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTickLength(0.2);
	frame->GetYaxis()->SetTickLength(0.2);
    frame->GetXaxis()->SetTickSize(0.02);
    frame->GetYaxis()->SetTickSize(0.02);
    
	frame->Draw();
	outerBand->Draw("3same");
	innerBand->Draw("3same");
	expG->Draw("lsame");
	if (!blindData) obsG->Draw("lpsame");


	TLegend *leg = new TLegend(0.25, 0.6, 0.60, 0.84);
	leg->SetFillColor(0);
	leg->SetHeader("95% CL upper limits");	//,"C");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	if (!blindData) leg->AddEntry(obsG, "Observed", "lp");
	leg->AddEntry(expG, "Expected", "l");
	leg->AddEntry(innerBand, "68% expected", "f");
	leg->AddEntry(outerBand, "95% expected", "f");
	leg->Draw();

	extraText = "Private work";
	writeExtraText = false;
	CMS_lumi(canv, 5, 33);
	
    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(3.6,2.01,"#times 10^{-3}");
	
      
	canv->RedrawAxis();

	leg->Draw();
	canv->Update();
	canv->Print("BR_limits.pdf");

}
