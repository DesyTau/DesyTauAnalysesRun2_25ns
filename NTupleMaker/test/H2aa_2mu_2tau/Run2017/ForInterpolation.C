#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "TLegend.h"
using namespace RooFit;

void GetFittingPar()
{

	gROOT->SetBatch();
	// Silence INFO messages
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	// Silence additional MINUIT output 
	RooMsgService::instance().setSilentMode(true);
	// Silence Info in TCanvas::Print: messages
	gErrorIgnoreLevel = kWarning;

	//// Mass points	////
	int nSamples = 19;

	float points[19] = { 3.6, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };

	TString mpoints[19] = { "3p6", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21" };

	float TrkTrk_DR_Up[19] = { 0.3, 0.3, 0.4, 0.5, 0.6, 0.8, 0.8, 0.8, 0.8, 1, 1, 1.2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 };

	///////////////////////////////// Category Strings	//////////////////////////////////////////
	int nCategories = 3;
	TString categories[3] = { "lep_lep", "lep_had", "had_had" };

	///////////////////////////////// Fit parameters	////////////////////////////////////////////
	float Par_1[3] = { 0.0355363, 0.0353708, 0.0355238 };
	float Par_2[3] = { 0.00562332, 0.00581018, 0.00704332 };
	float Par_3[3] = { 0.0453863, 0.0477197, 0.0520917 };
	float Par_4[3] = { 0.0125874, 0.0123017, 0.0116557 };
	float Par_5[3] = { 0.496051, 0.566405, 0.716669 };

	float Par_6[3] = { 0.03798388, 0.04020731, 0.04233213 };
	float Par_7[3] = { 0.00327519, 0.0056714, 0.007209588 };
	float Par_8[3] = { 0, 0.0296614, 0.227886 };
	float Par_9[3] = { 0.004320131, 0.003069381, 0.002486288 };
	float Par_10[3] = { 0.03722824, 0.03760699, 0.03452946 };
	float Par_11[3] = { -0.07015738, -0.08169356, -0.09413662 };
	float Par_12[3] = { -0.01137329, -0.0177734, -0.02389334 };

	float Par_13[3] = { 128.693, 126.383, 120.025 };
	float Par_14[3] = { 6.10861, 4.28535, 5.72352 };
	float Par_15[3] = { 19.2103, 23.3101, 22.9058 };

	float df[3] = { 2, 3, 3 };

	///////////////////////////////// Parameter Strings	//////////////////////////////////////////
	int nParameters = 19;
	TString Parameters[19] = { "mean_V", "width_V", "sigma_V", "mean_G1", "sigma_G1", "mean_G2", "sigma_G2", "frac", "mean_L1", "sigma_L1", "mean_G3", "sigma_G3", "mean_L2", "sigma_L2", "frac1", "mean_L3", "sigma_L3", "mean_G4", "sigma_G4" };

	TGraphErrors *gr[nParameters];

	for (int i = 0; i < nParameters; i++)
	{
		gr[i] = new TGraphErrors(nSamples);
	}

	///////////////////////////////// Directory String	//////////////////////////////////////////             
	TString dir = "/nfs/dust/cms/user/consuegs/Analyses/H2aa_2mu2tau/Run2017/";

	ostringstream str;
	str << "Chi2_Summary.txt";
	string nn = str.str();
	const char *p = nn.c_str();
	std::ofstream textFile(p);

	for (int icategory = 0; icategory < nCategories; icategory++)	// Loop over categories
	{

		TFile *fileInput_10 = new TFile(dir + "MVA_BDT/MergeTrees/tree_s_ma10_file.root");

		std::cout << std::endl << std::endl;
		std::cout << "--- --------------------------------------------------------- --- " << std::endl;
		std::cout << "       *** Saving fit parameters for category: " << categories[icategory] << " ***" << std::endl;
		std::cout << "--- --------------------------------------------------------- --- " << std::endl;
		std::cout << std::endl << std::endl;

		TTree *tree_s_10 = (TTree*) fileInput_10->Get("tree_s_" + categories[icategory]);

		//// Declare observable for MuMuTrkTrkMET_Mass	//// 
		RooRealVar MuMuTrkTrkMET_Mass("MuMuTrkTrkMET_Mass", "MuMuTrkTrkMET_Mass", 50, 250);

		//// For Convolution of Landau and Gaussian (MuMuTrkTrkMET_Mass)	////      
		RooRealVar mean_L3("mean_L3", "mean_L3", Par_13[icategory], Par_13[icategory] - 2, Par_13[icategory] + 2);
		RooRealVar sigma_L3("sigma_L3", "sigma_L3", Par_14[icategory], Par_14[icategory] - 2, Par_14[icategory] + 2);
		RooRealVar mean_G4("mean_G4", "mean_G4", 0.1);
		RooRealVar sigma_G4("sigma_G4", "sigma_G4", Par_15[icategory], Par_15[icategory] - 2, Par_15[icategory] + 2);

		//// Declare Functions for MuMuTrkTrkMET_Mass	//// 
		RooLandau Landau3("Landau3", "Landau3", MuMuTrkTrkMET_Mass, mean_L3, sigma_L3);
		RooGaussian Gaussian4("Gaussian4", "Gaussian4", MuMuTrkTrkMET_Mass, mean_G4, sigma_G4);

		//// Add pdfs to be used in model for MuMuTrkTrkMET_Mass	//// 
		RooAbsPdf *Signal_Model_MuMuTrkTrkMET_Mass = new RooFFTConvPdf("Signal_Model_MuMuTrkTrkMET_Mass", "Signal_Model_MuMuTrkTrkMET_Mass", MuMuTrkTrkMET_Mass, Landau3, Gaussian4);

		//// Create a binned dataset that imports content of TH1 and associates its contents to observable	//// 
		RooDataSet dataSR_10("dataSR", "dataSR", RooArgSet(MuMuTrkTrkMET_Mass), Import(*tree_s_10));

		//// Fitting Signal Model to Histogram	//// 
		RooFitResult * r_MuMuTrkTrkMET_Mass;

		r_MuMuTrkTrkMET_Mass = Signal_Model_MuMuTrkTrkMET_Mass->fitTo(dataSR_10, Save());

		//// For computation of the chi square/ndof appearing in the plots	////  
		RooPlot *frame_MuMuTrkTrkMET_Mass = MuMuTrkTrkMET_Mass.frame(Bins(50));
		dataSR_10.plotOn(frame_MuMuTrkTrkMET_Mass);
		Signal_Model_MuMuTrkTrkMET_Mass->plotOn(frame_MuMuTrkTrkMET_Mass);

		textFile << "Chi2/NoF: " << categories[icategory] << "--> MuMuTrkTrkMET_Mass" << "-->" << mpoints[8] << " = " << frame_MuMuTrkTrkMET_Mass->chiSquare(3) << std::endl;

		for (int isample = 0; isample < nSamples; isample++)	// Loop over samples
		{

			TFile *fileInput = new TFile(dir + "MVA_BDT/MergeTrees/tree_s_ma" + mpoints[isample] + "_file.root");

			TTree *tree_s = (TTree*) fileInput->Get("tree_s_" + categories[icategory]);

			//// Declare observables	//// 
			RooRealVar MuMu_Mass("MuMu_Mass", "MuMu_Mass", 3.5, 22);
			RooRealVar TrkTrk_DR("TrkTrk_DR", "TrkTrk_DR", 0, 1.5);
			RooRealVar MuMu_DR("MuMu_DR", "MuMu_DR", 0., 1.5);

            //// For Voigtian (MuMu_Mass)	//// 
			RooRealVar mean_V("mean_V", "mean_V", points[isample], points[isample] - 0.02, points[isample] + 0.02);
			RooRealVar width_V("width_V", "width_V", 0.00635518 *points[isample], 0.00635518 *points[isample] - 0.0025, 0.00635518 *points[isample] + 0.0025);
			RooRealVar sigma_V("sigma_V", "sigma_V", 0.00802073 *points[isample], 0.00802073 *points[isample] - 0.0025, 0.00802073 *points[isample] + 0.0025);

			//// For 2 Gaussian (MuMu_DR)	////        
			RooRealVar mean_G1("mean_G1", "mean_G1", Par_1[icategory] *points[isample], Par_1[icategory] *points[isample] - 0.05, Par_1[icategory] *points[isample] + 0.05);
			RooRealVar sigma_G1("sigma_G1", "sigma_G1", Par_2[icategory] *points[isample], Par_2[icategory] *points[isample] - 0.005, Par_2[icategory] *points[isample] + 0.005);
			RooRealVar mean_G2("mean_G2", "mean_G2", Par_3[icategory] *points[isample], Par_3[icategory] *points[isample] - 0.05, Par_3[icategory] *points[isample] + 0.05);
			RooRealVar sigma_G2("sigma_G2", "sigma_G2", Par_4[icategory] *points[isample], Par_4[icategory] *points[isample] - 0.005, Par_4[icategory] *points[isample] + 0.005);
			RooRealVar frac("frac", "frac", Par_5[icategory]);

			//// For 2 Landau (TrkTrk_DR)	////
			RooRealVar mean_L1("mean_L1", "mean_L1", Par_6[icategory] *points[isample] + Par_11[icategory], Par_6[icategory] *points[isample] + Par_11[icategory] - 0.02, Par_6[icategory] *points[isample] + Par_11[icategory] + 0.02);
			RooRealVar sigma_L1("sigma_L1", "sigma_L1", Par_7[icategory] *points[isample] + Par_12[icategory], Par_7[icategory] *points[isample] + Par_12[icategory] - 0.002, Par_7[icategory] *points[isample] + Par_12[icategory] + 0.002);
			RooRealVar mean_G3("mean_G3", "mean_G3", 0);
			RooRealVar sigma_G3("sigma_G3", "sigma_G3", Par_9[icategory] *points[isample] + Par_10[icategory], Par_9[icategory] *points[isample] + Par_10[icategory] - 0.002, Par_9[icategory] *points[isample] + Par_10[icategory] + 0.002);
			RooRealVar mean_L2("mean_L2", "mean_L2", 0.0610851);
			RooRealVar sigma_L2("sigma_L2", "sigma_L2", 0.0264235);
			RooRealVar frac1("frac1", "frac1", Par_8[icategory]);

			//// Declare Functions	////
			RooVoigtian Voigtian("Voigtian", "Voigtian", MuMu_Mass, mean_V, width_V, sigma_V);
			RooGaussian Gaussian1("Gaussian1", "Gaussian1", MuMu_DR, mean_G1, sigma_G1);
			RooGaussian Gaussian2("Gaussian2", "Gaussian2", MuMu_DR, mean_G2, sigma_G2);
			RooLandau Landau1("Landau1", "Landau1", TrkTrk_DR, mean_L1, sigma_L1);
			RooGaussian Gaussian3("Gaussian3", "Gaussian3", TrkTrk_DR, mean_G3, sigma_G3);
			RooFFTConvPdf Convoluted("Convoluted", "Convoluted", TrkTrk_DR, Landau1, Gaussian3);
			RooLandau Landau2("Landau2", "Landau2", TrkTrk_DR, mean_L2, sigma_L2);

			//// Add pdfs to be used in model	//// 
			RooAbsPdf *Signal_Model_MuMu_Mass = new RooVoigtian("Signal_Model_MuMu_Mass", "Signal_Model_MuMu_Mass", MuMu_Mass, mean_V, width_V, sigma_V);
			RooAbsPdf *Signal_Model_MuMu_DR = new RooAddPdf("Signal_Model_MuMu_DR", "Signal_Model_MuMu_DR", RooArgList(Gaussian1, Gaussian2), RooArgList(frac));
			RooAbsPdf *Signal_Model_TrkTrk_DR = new RooAddPdf("Signal_Model_TrkTrk_DR", "Signal_Model_TrkTrk_DR", RooArgList(Landau2, Convoluted), RooArgList(frac1));

			//// Create a binned datasets that imports contents of TH1 and associates its contents to observable	////
			RooDataSet dataSR("dataSR", "dataSR", RooArgSet(MuMu_Mass, MuMu_DR, TrkTrk_DR, MuMuTrkTrkMET_Mass), Import(*tree_s));

			//// Fitting Signal Models to Histograms	////
			RooFitResult * r_MuMu_Mass;
			RooFitResult * r_TrkTrk_DR;
			RooFitResult * r_MuMu_DR;

			r_MuMu_Mass = Signal_Model_MuMu_Mass->fitTo(dataSR, Save());
			r_TrkTrk_DR = Signal_Model_TrkTrk_DR->fitTo(dataSR, Save());
			r_MuMu_DR = Signal_Model_MuMu_DR->fitTo(dataSR, Save());

			//// For computation of the chi square/ndof appearing in the plots	//// 
			RooPlot *frame_MuMu_Mass = MuMu_Mass.frame(points[isample] - 1., points[isample] + 1, 500);
			dataSR.plotOn(frame_MuMu_Mass);
			Signal_Model_MuMu_Mass->plotOn(frame_MuMu_Mass);
			RooPlot *frame_TrkTrk_DR = TrkTrk_DR.frame(0, TrkTrk_DR_Up[isample], 100);
			dataSR.plotOn(frame_TrkTrk_DR);
			Signal_Model_TrkTrk_DR->plotOn(frame_TrkTrk_DR);
			RooPlot *frame_MuMu_DR = MuMu_DR.frame(0, 1.5, 50);
			dataSR.plotOn(frame_MuMu_DR);
			Signal_Model_MuMu_DR->plotOn(frame_MuMu_DR);

			textFile << "Chi2/NoF: " << categories[icategory] << "--> MuMu_Mass" << "-->" << mpoints[isample] << " = " << frame_MuMu_Mass->chiSquare(3) << std::endl;
			textFile << "Chi2/NoF: " << categories[icategory] << "--> MuMu_DR" << "-->" << mpoints[isample] << " = " << frame_MuMu_DR->chiSquare(4) << std::endl;
			textFile << "Chi2/NoF: " << categories[icategory] << "--> TrkTrk_DR" << "-->" << mpoints[isample] << " = " << frame_TrkTrk_DR->chiSquare(df[icategory]) << std::endl;

			gr[0]->SetPoint(isample, points[isample], mean_V.getVal());
			gr[0]->SetPointError(isample, 0., mean_V.getError());
			gr[0]->GetYaxis()->SetRangeUser(0, 35);
			gr[1]->SetPoint(isample, points[isample], width_V.getVal());
			gr[1]->SetPointError(isample, 0., width_V.getError());
			gr[1]->GetYaxis()->SetRangeUser(0, width_V.getVal() + 0.1);
			gr[2]->SetPoint(isample, points[isample], sigma_V.getVal());
			gr[2]->SetPointError(isample, 0., sigma_V.getError());
			gr[2]->GetYaxis()->SetRangeUser(0, sigma_V.getVal() + 0.2);
			gr[3]->SetPoint(isample, points[isample], mean_G1.getVal());
			gr[3]->SetPointError(isample, 0., mean_G1.getError());
			gr[3]->GetYaxis()->SetRangeUser(0, mean_G1.getVal() + 0.5);
			gr[4]->SetPoint(isample, points[isample], sigma_G1.getVal());
			gr[4]->SetPointError(isample, 0., sigma_G1.getError());
			gr[4]->GetYaxis()->SetRangeUser(0, sigma_G1.getVal() + 0.1);
			gr[5]->SetPoint(isample, points[isample], mean_G2.getVal());
			gr[5]->SetPointError(isample, 0., mean_G2.getError());
			gr[5]->GetYaxis()->SetRangeUser(0, mean_G2.getVal() + 1);
			gr[6]->SetPoint(isample, points[isample], sigma_G2.getVal());
			gr[6]->SetPointError(isample, 0., sigma_G2.getError());
			gr[6]->GetYaxis()->SetRangeUser(0, sigma_G2.getVal() + 0.1);
			gr[7]->SetPoint(isample, points[isample], frac.getVal());
			gr[7]->SetPointError(isample, 0., frac.getError());
			gr[7]->GetYaxis()->SetRangeUser(frac.getVal() - 0.1, frac.getVal() + 0.1);
			gr[8]->SetPoint(isample, points[isample], mean_L1.getVal());
			gr[8]->SetPointError(isample, 0., mean_L1.getError());
			gr[8]->GetYaxis()->SetRangeUser(0, mean_L1.getVal() + 0.5);
			gr[9]->SetPoint(isample, points[isample], sigma_L1.getVal());
			gr[9]->SetPointError(isample, 0., sigma_L1.getError());
			gr[9]->GetYaxis()->SetRangeUser(0, sigma_L1.getVal() + 0.1);
			gr[10]->SetPoint(isample, points[isample], mean_G3.getVal());
			gr[10]->SetPointError(isample, 0., mean_G3.getError());
			gr[10]->GetYaxis()->SetRangeUser(mean_G3.getVal() - 0.1, mean_G3.getVal() + 0.1);
			gr[11]->SetPoint(isample, points[isample], sigma_G3.getVal());
			gr[11]->SetPointError(isample, 0., sigma_G3.getError());
			gr[11]->GetYaxis()->SetRangeUser(sigma_G3.getVal() - 0.1, sigma_G3.getVal() + 0.1);
			gr[12]->SetPoint(isample, points[isample], mean_L2.getVal());
			gr[12]->SetPointError(isample, 0., mean_L2.getError());
			gr[12]->GetYaxis()->SetRangeUser(mean_L2.getVal() - 0.1, mean_L2.getVal() + 0.1);
			gr[13]->SetPoint(isample, points[isample], sigma_L2.getVal());
			gr[13]->SetPointError(isample, 0., sigma_L2.getError());
			gr[13]->GetYaxis()->SetRangeUser(sigma_L2.getVal() - 0.1, sigma_L2.getVal() + 0.1);
			gr[14]->SetPoint(isample, points[isample], frac1.getVal());
			gr[14]->SetPointError(isample, 0., frac1.getError());
			gr[14]->GetYaxis()->SetRangeUser(Par_8[icategory] - 0.05, Par_8[icategory] + 0.05);
			gr[15]->SetPoint(isample, points[isample], mean_L3.getVal());
			gr[15]->SetPointError(isample, 0., mean_L3.getError());
			gr[15]->GetYaxis()->SetRangeUser(mean_L3.getVal() - 10, mean_L3.getVal() + 10);
			gr[16]->SetPoint(isample, points[isample], sigma_L3.getVal());
			gr[16]->SetPointError(isample, 0., sigma_L3.getError());
			gr[16]->GetYaxis()->SetRangeUser(sigma_L3.getVal() - 10, sigma_L3.getVal() + 10);
			gr[17]->SetPoint(isample, points[isample], mean_G4.getVal());
			gr[17]->SetPointError(isample, 0., mean_G4.getError());
			gr[17]->GetYaxis()->SetRangeUser(mean_G4.getVal() - 0.1, mean_G4.getVal() + 0.1);
			gr[18]->SetPoint(isample, points[isample], sigma_G4.getVal());
			gr[18]->SetPointError(isample, 0., sigma_G4.getError());
			gr[18]->GetYaxis()->SetRangeUser(sigma_G4.getVal() - 10, sigma_G4.getVal() + 10);
		}

		//// Save TGraphErrors which will serve as input for interpolation procedure for mass points between the 1 GeV step	////
		TFile f("Parameters_ForInterpolation.root", "UPDATE");
		TCanvas *canv = MakeCanvas("canv", "canv", 900, 900);
		canv->SetLeftMargin(0.12);
		canv->SetBottomMargin(0.1);

		for (int i = 0; i < nParameters; i++)
		{
			gr[i]->Write(categories[icategory] + "_" + Parameters[i]);
			gr[i]->SetTitle("");
			gr[i]->SetMarkerColor(kBlue);
			gr[i]->SetMarkerStyle(21);
			gr[i]->Draw("ALP");
			gr[i]->GetYaxis()->SetTitle("Parameter value");
			gr[i]->GetYaxis()->SetTitleOffset(1.5);
			gr[i]->GetXaxis()->SetTitle("m_{a_{1}}[GeV]");

			TLegend *leg = new TLegend(0.66, 0.72, 0.92, 0.74);
			leg->SetTextFont(42);
			leg->SetTextSize(0.04);
			leg->SetLineWidth(0);
			leg->AddEntry(gr[i], Parameters[i], "lp");
			leg->Draw();

			writeExtraText = true;
			extraText = "Work in progress";
			CMS_lumi(canv, 5, 33);
			plotchannel(" _" + categories[icategory] + " Channel");

			canv->SaveAs(dir + "GetFittingPar/" + Parameters[i] + "_" + categories[icategory] + ".pdf");
			canv->Update();
		}
		delete canv;
	}
}

void Wspacewrite()
{

	// Silence INFO messages
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	// Silence additional MINUIT output 
	RooMsgService::instance().setSilentMode(true);
	// Silence Info in TCanvas::Print: messages
	gErrorIgnoreLevel = kWarning;

	///////////////////////////////// Parameter Strings	//////////////////////////////////////////
	int nParameters = 19;
	TString Parameters[19] = { "mean_V", "width_V", "sigma_V", "mean_G1", "sigma_G1", "mean_G2", "sigma_G2", "frac", "mean_L1", "sigma_L1", "mean_G3", "sigma_G3", "mean_L2", "sigma_L2", "frac1", "mean_L3", "sigma_L3", "mean_G4", "sigma_G4" };

	//// Categories	////
	int nCategories = 3;
	TString categories[3] = { "lep_lep", "lep_had", "had_had" };

	//// Open root file containing Tgraphs for each of the 19 parameters	////
	TFile *inputfile = TFile::Open("Parameters_ForInterpolation.root");

	//// Workspace output file	////
	TFile *outputfile = new TFile("Workspace_Interpolation.root", "RECREATE");

	std::cout << std::endl << std::endl;
	std::cout << "--- ------------------------------------------------------------------------------- --- " << std::endl;
	std::cout << "       *** Saving workspace with fit parameters from interpolation procedure" << " ***" << std::endl;
	std::cout << "--- ------------------------------------------------------------------------------- --- " << std::endl;
	std::cout << std::endl << std::endl;

	for (int icategory = 0; icategory < nCategories; icategory++)	// Loop over categories
	{

		for (double mass = 3.60; mass <= 21.0; mass += 0.20)	// Loop over masses
		{

			//// To retrieve Tgraphs from root file	////
			TGraphErrors *gr[nParameters];
			TGraphErrors *gr_Up[nParameters];
			TGraphErrors *gr_Down[nParameters];

			for (int iparameter = 0; iparameter < nParameters; iparameter++)	// Loop over parameters
			{

				inputfile->GetObject(categories[icategory] + "_" + Parameters[iparameter], gr[iparameter]);

				gr_Up[iparameter] = (TGraphErrors*) gr[iparameter]->Clone();
				gr_Down[iparameter] = (TGraphErrors*) gr[iparameter]->Clone();

				for (int ipoint = 0; ipoint < gr[iparameter]->GetN(); ipoint++)	// Loop over points
				{
					double x, y;
					gr[iparameter]->GetPoint(ipoint, x, y);
					gr_Up[iparameter]->SetPoint(ipoint, x, y + gr[iparameter]->GetErrorY(ipoint));
					gr_Up[iparameter]->SetPointError(ipoint, 0., 0.);
					gr_Down[iparameter]->SetPoint(ipoint, x, y - gr[iparameter]->GetErrorY(ipoint));
					gr_Down[iparameter]->SetPointError(ipoint, 0., 0.);
				}	// End loop over points

			}	// End loop over parameters

			TString ws_sufix = std::to_string(int(mass + 0.1)) + "p" + std::to_string(int(int(mass *10 + 0.1) % 10));

			//// Declare observables	////
			RooRealVar MuMu_Mass("MuMu_Mass", "MuMu_Mass", 3.5, 22);
			RooRealVar TrkTrk_DR("TrkTrk_DR", "TrkTrk_DR", 0, 1.5);
			RooRealVar MuMu_DR("MuMu_DR", "MuMu_DR", 0., 1.5);
			RooRealVar MuMuTrkTrkMET_Mass("MuMuTrkTrkMET_Mass", "MuMuTrkTrkMET_Mass", 50, 250);

			//// For Voigtian (MuMu_Mass)	////
			RooRealVar mean_V("mean_V", "mean_V", gr[0]->Eval(mass, 0, "S"));
			mean_V.setError(fabs(gr_Up[0]->Eval(mass, 0, "S") - gr_Down[0]->Eval(mass, 0, "S")) / 2.);
			RooRealVar width_V("width_V", "width_V", gr[1]->Eval(mass, 0, "S"));
			width_V.setError(fabs(gr_Up[1]->Eval(mass, 0, "S") - gr_Down[1]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_V("sigma_V", "sigma_V", gr[2]->Eval(mass, 0, "S"));
			sigma_V.setError(fabs(gr_Up[2]->Eval(mass, 0, "S") - gr_Down[2]->Eval(mass, 0, "S")) / 2.);

			//// For Convolution of Landau and Gaussian (MuMu_DR)	////
			RooRealVar mean_G1("mean_G1", "mean_G1", gr[3]->Eval(mass, 0, "S"));
			mean_G1.setError(fabs(gr_Up[3]->Eval(mass, 0, "S") - gr_Down[3]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_G1("sigma_G1", "sigma_G1", gr[4]->Eval(mass, 0, "S"));
			sigma_G1.setError(fabs(gr_Up[4]->Eval(mass, 0, "S") - gr_Down[4]->Eval(mass, 0, "S")) / 2.);
			RooRealVar mean_G2("mean_G2", "mean_G2", gr[5]->Eval(mass, 0, "S"));
			mean_G2.setError(fabs(gr_Up[5]->Eval(mass, 0, "S") - gr_Down[5]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_G2("sigma_G2", "sigma_G2", gr[6]->Eval(mass, 0, "S"));
			sigma_G2.setError(fabs(gr_Up[6]->Eval(mass, 0, "S") - gr_Down[6]->Eval(mass, 0, "S")) / 2.);
			RooRealVar frac("frac", "frac", gr[7]->Eval(mass, 0, "S"));
			frac.setError(fabs(gr_Up[7]->Eval(mass, 0, "S") - gr_Down[7]->Eval(mass, 0, "S")) / 2.);

			//// For 2 Landau (TrkTrk_DR)	////
			RooRealVar mean_L1("mean_L1", "mean_L1", gr[8]->Eval(mass, 0, "S"));
			mean_L1.setError(fabs(gr_Up[8]->Eval(mass, 0, "S") - gr_Down[8]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_L1("sigma_L1", "sigma_L1", gr[9]->Eval(mass, 0, "S"));
			sigma_L1.setError(fabs(gr_Up[9]->Eval(mass, 0, "S") - gr_Down[9]->Eval(mass, 0, "S")) / 2.);
			RooRealVar mean_G3("mean_G3", "mean_G3", gr[10]->Eval(mass, 0, "S"));
			mean_G3.setError(fabs(gr_Up[10]->Eval(mass, 0, "S") - gr_Down[10]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_G3("sigma_G3", "sigma_G3", gr[11]->Eval(mass, 0, "S"));
			sigma_G3.setError(fabs(gr_Up[11]->Eval(mass, 0, "S") - gr_Down[11]->Eval(mass, 0, "S")) / 2.);
			RooRealVar mean_L2("mean_L2", "mean_L2", gr[12]->Eval(mass, 0, "S"));
			mean_L2.setError(fabs(gr_Up[12]->Eval(mass, 0, "S") - gr_Down[12]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_L2("sigma_L2", "sigma_L2", gr[13]->Eval(mass, 0, "S"));
			sigma_L2.setError(fabs(gr_Up[13]->Eval(mass, 0, "S") - gr_Down[13]->Eval(mass, 0, "S")) / 2.);
			RooRealVar frac1("frac1", "frac1", gr[14]->Eval(mass, 0, "S"));
			frac1.setError(fabs(gr_Up[14]->Eval(mass, 0, "S") - gr_Down[14]->Eval(mass, 0, "S")) / 2.);

			//// For Convolution of Landau and Gaussian (MuMuTrkTrkMET_Mass)	////
			RooRealVar mean_L3("mean_L3", "mean_L3", gr[15]->Eval(mass, 0, "S"));
			mean_L3.setError(fabs(gr_Up[15]->Eval(mass, 0, "S") - gr_Down[15]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_L3("sigma_L3", "sigma_L3", gr[16]->Eval(mass, 0, "S"));
			sigma_L3.setError(fabs(gr_Up[16]->Eval(mass, 0, "S") - gr_Down[16]->Eval(mass, 0, "S")) / 2.);
			RooRealVar mean_G4("mean_G4", "mean_G4", gr[17]->Eval(mass, 0, "S"));
			mean_G4.setError(fabs(gr_Up[17]->Eval(mass, 0, "S") - gr_Down[17]->Eval(mass, 0, "S")) / 2.);
			RooRealVar sigma_G4("sigma_G4", "sigma_G4", gr[18]->Eval(mass, 0, "S"));
			sigma_G4.setError(fabs(gr_Up[18]->Eval(mass, 0, "S") - gr_Down[18]->Eval(mass, 0, "S")) / 2.);

			//// Declare Functions	////
			RooVoigtian Voigtian("Voigtian", "Voigtian", MuMu_Mass, mean_V, width_V, sigma_V);
			RooGaussian Gaussian1("Gaussian1", "Gaussian1", MuMu_DR, mean_G1, sigma_G1);
			RooGaussian Gaussian2("Gaussian2", "Gaussian2", MuMu_DR, mean_G2, sigma_G2);
			RooLandau Landau1("Landau1", "Landau1", TrkTrk_DR, mean_L1, sigma_L1);
			RooGaussian Gaussian3("Gaussian3", "Gaussian3", TrkTrk_DR, mean_G3, sigma_G3);
			RooFFTConvPdf Convoluted("Convoluted", "Convoluted", TrkTrk_DR, Landau1, Gaussian3);
			RooLandau Landau2("Landau2", "Landau2", TrkTrk_DR, mean_L2, sigma_L2);
			RooLandau Landau3("Landau3", "Landau3", MuMuTrkTrkMET_Mass, mean_L3, sigma_L3);
			RooGaussian Gaussian4("Gaussian4", "Gaussian4", MuMuTrkTrkMET_Mass, mean_G4, sigma_G4);

			//// Add pdfs to be used in model	//// 
			RooAbsPdf *Signal_Model_MuMu_Mass = new RooVoigtian("Signal_Model_MuMu_Mass", "Signal_Model_MuMu_Mass", MuMu_Mass, mean_V, width_V, sigma_V);
			RooAbsPdf *Signal_Model_MuMu_DR = new RooAddPdf("Signal_Model_MuMu_DR", "Signal_Model_MuMu_DR", RooArgList(Gaussian1, Gaussian2), RooArgList(frac));
			RooAbsPdf *Signal_Model_TrkTrk_DR = new RooAddPdf("Signal_Model_TrkTrk_DR", "Signal_Model_TrkTrk_DR", RooArgList(Landau2, Convoluted), RooArgList(frac1));
			RooAbsPdf *Signal_Model_MuMuTrkTrkMET_Mass = new RooFFTConvPdf("Signal_Model_MuMuTrkTrkMET_Mass", "Signal_Model_MuMuTrkTrkMET_Mass", MuMuTrkTrkMET_Mass, Landau3, Gaussian4);

			//// Create multi-dimensional p.d.f., multiply Signal_Model_MuMu_Mass, Signal_Model_MuMu_DR, Signal_Model_TrkTrk_DR and Signal_Model_MuMuTrkTrkMET_Mass	////
			RooProdPdf model("model", "model", RooArgList(*Signal_Model_MuMu_Mass, *Signal_Model_MuMu_DR, *Signal_Model_TrkTrk_DR, *Signal_Model_MuMuTrkTrkMET_Mass));

			////Import model and all its components into the workspace	////
			RooWorkspace *ws = new RooWorkspace("ws_" + categories[icategory] + "_" + ws_sufix, "ws_" + categories[icategory] + "_" + ws_sufix);

			ws->import(model);
			ws->defineSet("parameters", "mean_V,width_V,sigma_V,mean_G1,sigma_G1,mean_G2,sigma_G2,frac,mean_L1,sigma_L1,mean_G3,sigma_G3,mean_L2,sigma_L2,frac1,mean_L3,sigma_L3,mean_G4,sigma_G4");
			ws->defineSet("observables", "MuMu_Mass,TrkTrk_DR,MuMu_DR,MuMuTrkTrkMET_Mass");

			//// Print workspace contents	////
			//ws->Print();

			//// Save the workspace into a ROOT file	////
			outputfile->cd();
			ws->Write();
		}	// End loop over masses

	}	// End loop over categories      

}
