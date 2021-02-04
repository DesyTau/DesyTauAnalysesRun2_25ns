#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////           It contains 5 void functions:                                                                                     	////
////                                                                                                                             	////                               
////            1-) BDTClassificationSignal -> Classification Output of BDT for all Signal Samples                               	////
////            2-) BDTClassificationData   -> Classification Output of BDT for Background Model, Background Tests and Data in SR	////
////            3-) InterpolateAcceptance   -> Create function of Acceptance vs Mass for Interpolation                           	//// 
////            4-) InputsDataCards         -> Creation of DataCards and corresponding input files for limit computation         	////
////            5-) CreateInputs            -> It runs for all samples the classifications and creates the datacards             	////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//********************************************************//
//**************Switch of input variables ***************//
//********************************************************//
bool MuMu_Mass_isIN = true;
bool MuMu_Pt_isIN = false;
bool MuMu_DR_isIN = true;
bool MuMuMET_DPhi_isIN = false;
bool TrkTrk_Mass_isIN = false;
bool TrkTrk_Pt_isIN = false;
bool TrkTrk_DR_isIN = true;
bool TrkTrkMET_DPhi_isIN = false;
bool MuMuTrkTrk_Mass_isIN = false;
bool MuMuTrkTrk_Pt_isIN = false;
bool MuMuTauTau_Mass_isIN = false;
bool MuMuTauTau_Pt_isIN = false;
bool MuMuTrkTrkMET_Mass_isIN = true;
bool MET_Pt_isIN = false;

void BDTClassificationSignal(TString mass = "15p0")
{

	float MuMu_Mass, MuMu_Pt, MuMu_DR, MuMuMET_DPhi, TrkTrk_Mass, TrkTrk_Pt, TrkTrk_DR, TrkTrkMET_DPhi,
	MuMuTrkTrk_Mass, MuMuTrkTrk_Pt, MuMuTauTau_Mass, MuMuTauTau_Pt, MET_Pt, MuMuTrkTrkMET_Mass;
	float mvavalue;

	///////////////////////////////// Category Strings	////////////////////////////////////////////
	int ntrackCat = 3;

	TString trackCat[] = { "_lep_lep", "_lep_had", "_had_had" };

	///////////////////////////////// Region Strings	////////////////////////////////////////////
	int nRegions = 1;

	TString regions[1] = { "Sel" };

	std::cout << std::endl << std::endl;
	std::cout << "                    Start BDT Classification" << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << "                    ****     ***     *******                   " << std::endl;
	std::cout << "                    *   *    *   *      *                      " << std::endl;
	std::cout << "                    * **     *   *      *                      " << std::endl;
	std::cout << "                    *   *    *   *      *                      " << std::endl;
	std::cout << "                    ****     ***        *                      " << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "--- ------------------------------------------------------------" << std::endl;
	std::cout << "***         Starting loop over the Signal Samples           ***" << std::endl;
	std::cout << "--- ------------------------------------------------------------" << std::endl;
	std::cout << std::endl << std::endl;

	//**********Input &Output Files ***********//

	TFile *fileIn = new TFile("../Workspace_Interpolation.root");
	TFile *fileOut = new TFile("SUSYGGH_BDTOutput_M-" + mass + ".root", "recreate");

	std::cout << std::endl;
	std::cout << "--- BDT Classification    : Using input file: " << fileIn->GetName() << std::endl;

	///////////////////////////////////////////////////////////////////////
	//**********BDT Classification in Categories and Regions ***********//
	///////////////////////////////////////////////////////////////////////

	for (int icat = 0; icat < ntrackCat; icat++)	// Loop over regions
	{

		fileOut->cd();

		//**********Loading TMVA Reader ***********//

		TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

		if (MuMu_Mass_isIN)
			reader->AddVariable("MuMu_Mass", &MuMu_Mass);
		if (MuMu_Pt_isIN)
			reader->AddVariable("MuMu_Pt", &MuMu_Pt);
		if (MuMu_DR_isIN)
			reader->AddVariable("MuMu_DR", &MuMu_DR);
		if (MuMuMET_DPhi_isIN)
			reader->AddVariable("MuMuMET_DPhi", &MuMuMET_DPhi);
		if (TrkTrk_Mass_isIN)
			reader->AddVariable("TrkTrk_Mass", &TrkTrk_Mass);
		if (TrkTrk_Pt_isIN)
			reader->AddVariable("TrkTrk_Pt", &TrkTrk_Pt);
		if (TrkTrk_DR_isIN)
			reader->AddVariable("TrkTrk_DR", &TrkTrk_DR);
		if (TrkTrkMET_DPhi_isIN)
			reader->AddVariable("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
		if (MuMuTrkTrk_Mass_isIN)
			reader->AddVariable("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
		if (MuMuTrkTrk_Pt_isIN)
			reader->AddVariable("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
		if (MuMuTauTau_Mass_isIN)
			reader->AddVariable("MuMuTauTau_Mass", &MuMuTauTau_Mass);
		if (MuMuTauTau_Pt_isIN)
			reader->AddVariable("MuMuTauTau_Pt", &MuMuTauTau_Pt);
		if (MuMuTrkTrkMET_Mass_isIN)
			reader->AddVariable("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);
		if (MET_Pt_isIN)
			reader->AddVariable("MET_Pt", &MET_Pt);

		reader->BookMVA("BDT", "/nfs/dust/cms/user/consuegs/Analyses/H2aa_2mu2tau/Run2016/MVA_BDT/ma" + mass + "/weights/TMVAClassification_BDT" + trackCat[icat] + ".weights.xml");

		//**********Generate and Classify ***********//

		for (int iregion = 0; iregion < nRegions; iregion++)	// Loop over regions
		{

			double MuMu_Mass_D, MuMu_Pt_D, MuMu_DR_D, MuMuMET_DPhi_D, TrkTrk_Mass_D, TrkTrk_Pt_D, TrkTrk_DR_D, TrkTrkMET_DPhi_D,
			MuMuTrkTrk_Mass_D, MuMuTrkTrk_Pt_D, MuMuTauTau_Mass_D, MuMuTauTau_Pt_D, MET_Pt_D, MuMuTrkTrkMET_Mass_D;

			////////// Retrieving Workspace for Event Generation	//////////

			RooWorkspace *ws = (RooWorkspace*) fileIn->Get("ws" + trackCat[icat] + "_" + mass);

			RooAbsPdf *SignalModel = ws->pdf("model");

			RooDataSet *Dataset_Old = SignalModel->generate(*ws->set("observables"), 10000);

			RooAbsData::setDefaultStorageType(RooAbsData::Tree);

			RooDataSet *Dataset = new RooDataSet("Dataset", "Dataset", Dataset_Old, *Dataset_Old->get());

			TTree *tree = (TTree*) Dataset->tree()->Clone();

			if (MuMu_Mass_isIN) tree->SetBranchAddress("MuMu_Mass", &MuMu_Mass_D);
			if (MuMu_Pt_isIN) tree->SetBranchAddress("MuMu_Pt", &MuMu_Pt_D);
			if (MuMu_DR_isIN) tree->SetBranchAddress("MuMu_DR", &MuMu_DR_D);
			if (MuMuMET_DPhi_isIN) tree->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi_D);
			if (TrkTrk_Mass_isIN) tree->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass_D);
			if (TrkTrk_Pt_isIN) tree->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt_D);
			if (TrkTrk_DR_isIN) tree->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR_D);
			if (TrkTrkMET_DPhi_isIN) tree->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi_D);
			if (MuMuTrkTrk_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass_D);
			if (MuMuTrkTrk_Pt_isIN) tree->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt_D);
			if (MuMuTauTau_Mass_isIN) tree->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass_D);
			if (MuMuTauTau_Pt_isIN) tree->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt_D);
			if (MET_Pt_isIN) tree->SetBranchAddress("MET_Pt", &MET_Pt_D);
			if (MuMuTrkTrkMET_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass_D);

			////////// Classifiying	//////////   

			TH1D *Histo = new TH1D("MVABDTOutput" + regions[iregion] + trackCat[icat] + "H", "", 2000, -1., 1.);

			for (int ievt = 0; ievt < tree->GetEntries(); ievt++)	// Loop over tree entries
			{
				tree->GetEntry(ievt);

				MuMu_Mass = float_t(MuMu_Mass_D);
				MuMu_Pt = float_t(MuMu_Pt_D);
				MuMu_DR = float_t(MuMu_DR_D);
				MuMuMET_DPhi = float_t(MuMuMET_DPhi_D);
				TrkTrk_Mass = float_t(TrkTrk_Mass_D);
				TrkTrk_Pt = float_t(TrkTrk_Pt_D);
				TrkTrk_DR = float_t(TrkTrk_DR_D);
				TrkTrkMET_DPhi = float_t(TrkTrkMET_DPhi_D);
				MuMuTrkTrk_Mass = float_t(MuMuTrkTrk_Mass_D);
				MuMuTrkTrk_Pt = float_t(MuMuTrkTrk_Pt_D);
				MuMuTauTau_Mass = float_t(MuMuTauTau_Mass_D);
				MuMuTauTau_Pt = float_t(MuMuTauTau_Pt_D);
				MET_Pt = float_t(MET_Pt_D);
				MuMuTrkTrkMET_Mass = float_t(MuMuTrkTrkMET_Mass_D);

				mvavalue = reader->EvaluateMVA("BDT");

				Histo->Fill(mvavalue);
			}	// End loop over tree entries

			fileOut->cd();
			Histo->Write("MVABDTOutput" + regions[iregion] + trackCat[icat] + "H");

			//**********Propagation of parameter uncertainties ***********//

			int nParameters = 13;

			TString Parameters[13] = { "mean_V", "width_V", "sigma_V", "mean_G1", "sigma_G1", "mean_G2", "sigma_G2", "mean_L1", "sigma_L1", "sigma_G3", "mean_L3", "sigma_L3", "sigma_G4" };

			for (int iparameter = 0; iparameter < nParameters; iparameter++)	// Loop over parameters
			{

				///////// Shift Up	///////////// 
				ws->var(Parameters[iparameter])->setVal(ws->var(Parameters[iparameter])->getValV() + ws->var(Parameters[iparameter])->getError());

				RooDataSet *DatasetUncUp_Old = SignalModel->generate(*ws->set("observables"), 10000);

				RooDataSet *DatasetUncUp = new RooDataSet("DatasetUncUp", "DatasetUncUp", DatasetUncUp_Old, *DatasetUncUp_Old->get());

				TTree *treeUncUp = (TTree*) DatasetUncUp->tree()->Clone();

				if (MuMu_Mass_isIN) treeUncUp->SetBranchAddress("MuMu_Mass", &MuMu_Mass_D);
				if (MuMu_Pt_isIN) treeUncUp->SetBranchAddress("MuMu_Pt", &MuMu_Pt_D);
				if (MuMu_DR_isIN) treeUncUp->SetBranchAddress("MuMu_DR", &MuMu_DR_D);
				if (MuMuMET_DPhi_isIN) treeUncUp->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi_D);
				if (TrkTrk_Mass_isIN) treeUncUp->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass_D);
				if (TrkTrk_Pt_isIN) treeUncUp->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt_D);
				if (TrkTrk_DR_isIN) treeUncUp->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR_D);
				if (TrkTrkMET_DPhi_isIN) treeUncUp->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi_D);
				if (MuMuTrkTrk_Mass_isIN) treeUncUp->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass_D);
				if (MuMuTrkTrk_Pt_isIN) treeUncUp->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt_D);
				if (MuMuTauTau_Mass_isIN) treeUncUp->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass_D);
				if (MuMuTauTau_Pt_isIN) treeUncUp->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt_D);
				if (MET_Pt_isIN) treeUncUp->SetBranchAddress("MET_Pt", &MET_Pt_D);
				if (MuMuTrkTrkMET_Mass_isIN) treeUncUp->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass_D);

				//------- Classifiying -------//   

				TH1D *HistoUp = new TH1D("MVABDTOutput" + regions[iregion] + trackCat[icat] + "_" + Parameters[iparameter] + "_UncUpH", "", 2000, -1., 1.);

				for (int ievt = 0; ievt < treeUncUp->GetEntries(); ievt++)	// Loop over tree entries
				{
					treeUncUp->GetEntry(ievt);

					MuMu_Mass = float_t(MuMu_Mass_D);
					MuMu_Pt = float_t(MuMu_Pt_D);
					MuMu_DR = float_t(MuMu_DR_D);
					MuMuMET_DPhi = float_t(MuMuMET_DPhi_D);
					TrkTrk_Mass = float_t(TrkTrk_Mass_D);
					TrkTrk_Pt = float_t(TrkTrk_Pt_D);
					TrkTrk_DR = float_t(TrkTrk_DR_D);
					TrkTrkMET_DPhi = float_t(TrkTrkMET_DPhi_D);
					MuMuTrkTrk_Mass = float_t(MuMuTrkTrk_Mass_D);
					MuMuTrkTrk_Pt = float_t(MuMuTrkTrk_Pt_D);
					MuMuTauTau_Mass = float_t(MuMuTauTau_Mass_D);
					MuMuTauTau_Pt = float_t(MuMuTauTau_Pt_D);
					MET_Pt = float_t(MET_Pt_D);
					MuMuTrkTrkMET_Mass = float_t(MuMuTrkTrkMET_Mass_D);

					mvavalue = reader->EvaluateMVA("BDT");

					HistoUp->Fill(mvavalue);
				}	// End loop over tree entries

				fileOut->cd();
				HistoUp->Write("MVABDTOutput" + regions[iregion] + trackCat[icat] + "_" + Parameters[iparameter] + "_UncUpH");

				///////// Shift Down	///////////// 
				ws->var(Parameters[iparameter])->setVal(ws->var(Parameters[iparameter])->getValV() - ws->var(Parameters[iparameter])->getError());

				RooDataSet *DatasetUncDown_Old = SignalModel->generate(*ws->set("observables"), 10000);

				RooDataSet *DatasetUncDown = new RooDataSet("DatasetUncDown", "DatasetUncDown", DatasetUncDown_Old, *DatasetUncDown_Old->get());

				TTree *treeUncDown = (TTree*) DatasetUncDown->tree()->Clone();

				if (MuMu_Mass_isIN) treeUncDown->SetBranchAddress("MuMu_Mass", &MuMu_Mass_D);
				if (MuMu_Pt_isIN) treeUncDown->SetBranchAddress("MuMu_Pt", &MuMu_Pt_D);
				if (MuMu_DR_isIN) treeUncDown->SetBranchAddress("MuMu_DR", &MuMu_DR_D);
				if (MuMuMET_DPhi_isIN) treeUncDown->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi_D);
				if (TrkTrk_Mass_isIN) treeUncDown->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass_D);
				if (TrkTrk_Pt_isIN) treeUncDown->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt_D);
				if (TrkTrk_DR_isIN) treeUncDown->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR_D);
				if (TrkTrkMET_DPhi_isIN) treeUncDown->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi_D);
				if (MuMuTrkTrk_Mass_isIN) treeUncDown->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass_D);
				if (MuMuTrkTrk_Pt_isIN) treeUncDown->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt_D);
				if (MuMuTauTau_Mass_isIN) treeUncDown->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass_D);
				if (MuMuTauTau_Pt_isIN) treeUncDown->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt_D);
				if (MET_Pt_isIN) treeUncDown->SetBranchAddress("MET_Pt", &MET_Pt_D);
				if (MuMuTrkTrkMET_Mass_isIN) treeUncDown->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass_D);

				//------- Classifiying -------//   

				TH1D *HistoDown = new TH1D("MVABDTOutput" + regions[iregion] + trackCat[icat] + "_" + Parameters[iparameter] + "_UncDownH", "", 2000, -1., 1.);

				for (int ievt = 0; ievt < treeUncDown->GetEntries(); ievt++)	// Loop over tree entries
				{
					treeUncDown->GetEntry(ievt);

					MuMu_Mass = float_t(MuMu_Mass_D);
					MuMu_Pt = float_t(MuMu_Pt_D);
					MuMu_DR = float_t(MuMu_DR_D);
					MuMuMET_DPhi = float_t(MuMuMET_DPhi_D);
					TrkTrk_Mass = float_t(TrkTrk_Mass_D);
					TrkTrk_Pt = float_t(TrkTrk_Pt_D);
					TrkTrk_DR = float_t(TrkTrk_DR_D);
					TrkTrkMET_DPhi = float_t(TrkTrkMET_DPhi_D);
					MuMuTrkTrk_Mass = float_t(MuMuTrkTrk_Mass_D);
					MuMuTrkTrk_Pt = float_t(MuMuTrkTrk_Pt_D);
					MuMuTauTau_Mass = float_t(MuMuTauTau_Mass_D);
					MuMuTauTau_Pt = float_t(MuMuTauTau_Pt_D);
					MET_Pt = float_t(MET_Pt_D);
					MuMuTrkTrkMET_Mass = float_t(MuMuTrkTrkMET_Mass_D);

					mvavalue = reader->EvaluateMVA("BDT");

					HistoDown->Fill(mvavalue);
				}	// End loop over tree entries

				fileOut->cd();
				HistoDown->Write("MVABDTOutput" + regions[iregion] + trackCat[icat] + "_" + Parameters[iparameter] + "_UncDownH");
			}	// End loop over parameters

		}	// End loop over regions

	}	// End loop over categories

	//**********End of Classification ***********//

	std::cout << "--- Created root file: \"" << fileOut->GetName() << "\" containing the BDT output histograms" << std::endl;
	std::cout << "==> BDT Classification is done!" << std::endl << std::endl;
	fileOut->Close();

	std::cout << std::endl;
	std::cout << "End loop over the Signal" << std::endl;
	std::cout << std::endl << std::endl;
	fileIn->Close();
	std::cout << "                     End BDT Classification" << std::endl;

}

void BDTClassificationData(TString mass = "15p0")
{

	float MuMu_Mass, MuMu_Pt, MuMu_DR, MuMuMET_DPhi, TrkTrk_Mass, TrkTrk_Pt, TrkTrk_DR, TrkTrkMET_DPhi,
	MuMuTrkTrk_Mass, MuMuTrkTrk_Pt, MuMuTauTau_Mass, MuMuTauTau_Pt, MET_Pt, MuMuTrkTrkMET_Mass;
	float mvavalue;

	///////////////////////////////// Category Strings	////////////////////////////////////////////
	int ntrackCat = 3;

	TString trackCat[] = { "_lep_lep", "_lep_had", "_had_had" };

	///////////////////////////////// Region Strings	////////////////////////////////////////////
    int nRegions = 7;

    TString regions[7] = {"NN00", "NNNN", "00NN", "00SemiIso", "SoftIso", "00SoftIso", "Sel"};

	std::cout << std::endl << std::endl;
	std::cout << "                    Start BDT Classification" << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << "                    ****     ***     *******                   " << std::endl;
	std::cout << "                    *   *    *   *      *                      " << std::endl;
	std::cout << "                    * **     *   *      *                      " << std::endl;
	std::cout << "                    *   *    *   *      *                      " << std::endl;
	std::cout << "                    ****     ***        *                      " << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "--- ------------------------------------------------------------" << std::endl;
	std::cout << "***         Starting loop over the Data Samples              ***" << std::endl;
	std::cout << "--- ------------------------------------------------------------" << std::endl;
	std::cout << std::endl << std::endl;

	//**********Input & Output Files ***********//

	TFile *fileIn = new TFile("../SingleMuon_Run2016.root");
	TFile *fileInSS = new TFile("../SingleMuon_Run2016_SS.root");
	TFile *fileOut = new TFile("SingleMuon_BDTOutput_M-" + mass + ".root", "recreate");

	std::cout << std::endl;
	std::cout << "--- BDT Classification    : Using input file: " << fileIn->GetName() << std::endl;

	///////////////////////////////////////////////////////////////////////
	//**********BDT Classification in Categories and Regions ***********//
	///////////////////////////////////////////////////////////////////////

	for (int icat = 0; icat < ntrackCat; icat++)	// Loop over regions
	{

		fileOut->cd();

		///////////// Channel Strings	/////////////////////
		TString channels[9] = { "_mu_mu", "_mu_ele", "_ele_mu", "_ele_ele", "_mu_had", "_had_mu", "_ele_had", "_had_ele", "_had_had" };
		int nchannels = 0;
		int startat = 0;

		if (trackCat[icat] == "_lep_lep")
		{
			nchannels = 4;
			startat = 0;
		}
		if (trackCat[icat] == "_lep_had")
		{
			nchannels = 4;
			startat = 4;
		}
		if (trackCat[icat] == "_had_had")
		{
			nchannels = 1;
			startat = 8;
		}

		//**********Loading TMVA Reader ***********//

		TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

		if (MuMu_Mass_isIN)
			reader->AddVariable("MuMu_Mass", &MuMu_Mass);
		if (MuMu_Pt_isIN)
			reader->AddVariable("MuMu_Pt", &MuMu_Pt);
		if (MuMu_DR_isIN)
			reader->AddVariable("MuMu_DR", &MuMu_DR);
		if (MuMuMET_DPhi_isIN)
			reader->AddVariable("MuMuMET_DPhi", &MuMuMET_DPhi);
		if (TrkTrk_Mass_isIN)
			reader->AddVariable("TrkTrk_Mass", &TrkTrk_Mass);
		if (TrkTrk_Pt_isIN)
			reader->AddVariable("TrkTrk_Pt", &TrkTrk_Pt);
		if (TrkTrk_DR_isIN)
			reader->AddVariable("TrkTrk_DR", &TrkTrk_DR);
		if (TrkTrkMET_DPhi_isIN)
			reader->AddVariable("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
		if (MuMuTrkTrk_Mass_isIN)
			reader->AddVariable("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
		if (MuMuTrkTrk_Pt_isIN)
			reader->AddVariable("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
		if (MuMuTauTau_Mass_isIN)
			reader->AddVariable("MuMuTauTau_Mass", &MuMuTauTau_Mass);
		if (MuMuTauTau_Pt_isIN)
			reader->AddVariable("MuMuTauTau_Pt", &MuMuTauTau_Pt);
		if (MuMuTrkTrkMET_Mass_isIN)
			reader->AddVariable("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);
		if (MET_Pt_isIN)
			reader->AddVariable("MET_Pt", &MET_Pt);

		reader->BookMVA("BDT", "/nfs/dust/cms/user/consuegs/Analyses/H2aa_2mu2tau/Run2016/MVA_BDT/ma" + mass + "/weights/TMVAClassification_BDT" + trackCat[icat] + ".weights.xml");

		//**********Classify ***********//

		for (int iregion = 0; iregion < nRegions; iregion++)	// Loop over regions
		{

			TH1D *Histo = new TH1D("MVABDTOutput" + regions[iregion] + trackCat[icat] + "H", "", 2000, -1., 1.);

			for (int ichann = startat; ichann < (startat + nchannels); ichann++)
			{

				TTree *tree = (TTree*) fileIn->Get("tree_" + regions[iregion] + "_0p2" + channels[ichann]);
				if (regions[iregion] == "Sel") tree = (TTree*) fileIn->Get("tree" + regions[iregion] + "_0p2" + channels[ichann]);

				if (MuMu_Mass_isIN) tree->SetBranchAddress("MuMu_Mass", &MuMu_Mass);
				if (MuMu_Pt_isIN) tree->SetBranchAddress("MuMu_Pt", &MuMu_Pt);
				if (MuMu_DR_isIN) tree->SetBranchAddress("MuMu_DR", &MuMu_DR);
				if (MuMuMET_DPhi_isIN) tree->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi);
				if (TrkTrk_Mass_isIN) tree->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass);
				if (TrkTrk_Pt_isIN) tree->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt);
				if (TrkTrk_DR_isIN) tree->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR);
				if (TrkTrkMET_DPhi_isIN) tree->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
				if (MuMuTrkTrk_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
				if (MuMuTrkTrk_Pt_isIN) tree->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
				if (MuMuTauTau_Mass_isIN) tree->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass);
				if (MuMuTauTau_Pt_isIN) tree->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt);
				if (MET_Pt_isIN) tree->SetBranchAddress("MET_Pt", &MET_Pt);
				if (MuMuTrkTrkMET_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);

				////////// Classifiying	//////////   

				for (int ievt = 0; ievt < tree->GetEntries(); ievt++)	// Loop over tree entries
				{
					tree->GetEntry(ievt);

					mvavalue = reader->EvaluateMVA("BDT");

					Histo->Fill(mvavalue);
				}	// End loop over tree entries

			}	// End loop over channels

			fileOut->cd();
			Histo->Write("MVABDTOutput" + regions[iregion] + trackCat[icat] + "H");
		}	// End loop over regions
		
		
		////////////////////SS Signal Region Data/////////////////////////////////
		TH1D *Histo = new TH1D("MVABDTOutputSel_SS"+ trackCat[icat] + "H", "", 2000, -1., 1.);

		for (int ichann = startat; ichann < (startat + nchannels); ichann++)
		{

			TTree *tree = (TTree*) fileInSS->Get("treeSel_0p2" + channels[ichann]);
				
			if (MuMu_Mass_isIN) tree->SetBranchAddress("MuMu_Mass", &MuMu_Mass);
			if (MuMu_Pt_isIN) tree->SetBranchAddress("MuMu_Pt", &MuMu_Pt);
			if (MuMu_DR_isIN) tree->SetBranchAddress("MuMu_DR", &MuMu_DR);
			if (MuMuMET_DPhi_isIN) tree->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi);
			if (TrkTrk_Mass_isIN) tree->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass);
			if (TrkTrk_Pt_isIN) tree->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt);
			if (TrkTrk_DR_isIN) tree->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR);
			if (TrkTrkMET_DPhi_isIN) tree->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
			if (MuMuTrkTrk_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
			if (MuMuTrkTrk_Pt_isIN) tree->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
			if (MuMuTauTau_Mass_isIN) tree->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass);
			if (MuMuTauTau_Pt_isIN) tree->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt);
			if (MET_Pt_isIN) tree->SetBranchAddress("MET_Pt", &MET_Pt);
			if (MuMuTrkTrkMET_Mass_isIN) tree->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);

			////////// Classifiying	//////////   

			for (int ievt = 0; ievt < tree->GetEntries(); ievt++)	// Loop over tree entries
			{
				tree->GetEntry(ievt);

				mvavalue = reader->EvaluateMVA("BDT");

				Histo->Fill(mvavalue);
			}	// End loop over tree entries

		}	// End loop over channels

		fileOut->cd();
		Histo->Write("MVABDTOutputSel_SS" + trackCat[icat] + "H");

	}	// End loop over categories

	//**********End of Classification ***********//

	std::cout << "--- Created root file: \"" << fileOut->GetName() << "\" containing the BDT output histograms" << std::endl;
	std::cout << "==> BDT Classification is done!" << std::endl << std::endl;
	fileOut->Close();

	std::cout << std::endl;
	std::cout << "End loop over the Data" << std::endl;
	fileIn->Close();
	fileInSS->Close();
	std::cout << std::endl << std::endl;
	std::cout << "                     End BDT Classification" << std::endl;

}

void InterpolateAcceptance()
{

	///// Mass points	/////
	int nSamples = 19;

	float points[19] = { 3.6, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };

	TString mpoints[19] = { "3p6", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21" };

	///// Categories	/////
	int nCategories = 3;

	TString categories[3] = { "_lep_lep", "_lep_had", "_had_had" };

	TFile *fileOut = new TFile("Acceptance_ForInterpolation.root", "recreate");

	for (int icategory = 0; icategory < nCategories; icategory++)	// Loop over categories
	{

		TGraphErrors *gr = new TGraphErrors(nSamples);

		for (int isample = 0; isample < nSamples; isample++)	// Loop over samples
		{

			TFile *fileIn = new TFile("../SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-" + mpoints[isample] + ".root");

			TH1D *Initial_EventsH = (TH1D*) fileIn->Get("histWeightsH");
			TH1D *Final_EventsH = new TH1D("counter_FinalEvents" + categories[icategory] + "H", "", 1, 0., 2.);

			///////////// Channel Strings	/////////////////////
			TString channels[9] = { "_mu_mu", "_mu_ele", "_ele_mu", "_ele_ele", "_mu_had", "_had_mu", "_ele_had", "_had_ele", "_had_had" };
			int nchannels = 0;
			int startat = 0;

			if (categories[icategory] == "_lep_lep")
			{
				nchannels = 4;
				startat = 0;
			}
			if (categories[icategory] == "_lep_had")
			{
				nchannels = 4;
				startat = 4;
			}
			if (categories[icategory] == "_had_had")
			{
				nchannels = 1;
				startat = 8;
			}

			for (int ichann = startat; ichann < (startat + nchannels); ichann++) Final_EventsH->Add(Final_EventsH, (TH1D*) fileIn->Get("counter_FinalEventsH_0p2" + channels[ichann]), 1., 1.);

			gr->SetPoint(isample, points[isample], Final_EventsH->GetSumOfWeights() / Initial_EventsH->GetSumOfWeights());
			gr->SetPointError(isample, 0., TMath::Sqrt(Final_EventsH->GetSumOfWeights()) / Initial_EventsH->GetSumOfWeights());
		}	// End loop over samples

		fileOut->cd();
		gr->Write("Acceptance_vs_Mass" + categories[icategory]);
	}	// End loop over categories

	fileOut->Close();
}

void InputsDataCards(double mass = 15.0, TString trackCat = "_lep_had")
{

	TString mass_string = std::to_string(int(mass + 0.1)) + "p" + std::to_string(int(int(mass *10 + 0.1) % 10));

	std::cout << std::endl;
	std::cout << "Mass = " << mass << std::endl;
	std::cout << std::endl;

	double lumi = 35922;
	double xsecGGH = 48.52;

	float muMass = 0.105658;
	float tauMass = 1.777;
	float muMassToAMass = 2 *muMass / mass;
	float tauMassToAMass = 2 *tauMass / mass;

	float GammaAmumu = muMass *muMass *TMath::Sqrt(1 - muMassToAMass *muMassToAMass);
	float GammaAtautau = tauMass *tauMass *TMath::Sqrt(1 - tauMassToAMass *tauMassToAMass);

	float ratio = GammaAmumu / GammaAtautau;

	xsecGGH = xsecGGH *0.001;	// -> 2*B(H->aa)*Br(a->tautau)*B(a->mumu)

	TFile *file = new TFile("SingleMuon_BDTOutput_M-" + mass_string + ".root");
	TFile *fileGGH = new TFile("SUSYGGH_BDTOutput_M-" + mass_string + ".root");
	TFile *fileAcc = new TFile("Acceptance_ForInterpolation.root");

	///// Acceptance	/////
	TGraphErrors * gr;
	fileAcc->GetObject("Acceptance_vs_Mass" + trackCat, gr);

    TString BkgdModel, BkgdModel_unc;
    if (trackCat == "_lep_lep") {BkgdModel = "NNNN"; BkgdModel_unc = "Sel_SS";}
    else if (trackCat == "_lep_had") {BkgdModel = "NNNN"; BkgdModel_unc = "Sel_SS";}
    else  {BkgdModel = "NNNN"; BkgdModel_unc = "Sel_SS";}
        
	///// Shape (Nominal)	/////
	TH1D *hist_Old = (TH1D*) file->Get("MVABDTOutput" + BkgdModel + trackCat + "H");

	TH1D *histGGH_Old = (TH1D*) fileGGH->Get("MVABDTOutputSel" + trackCat + "H");

	///// Shape (Uncertainty)	/////
	TH1D *histUnc_Old = (TH1D*) file->Get("MVABDTOutput" + BkgdModel_unc + trackCat + "H");

	///// Data Observed	/////
	TH1D *hist_obs_Old = (TH1D*) file->Get("MVABDTOutputSel" + trackCat + "H");

	///////////// Finding Optimal Binning	////////////////
	hist_Old->Rebin(100);
	histGGH_Old->Rebin(100);
	histUnc_Old->Rebin(100);
	hist_obs_Old->Rebin(100);

	int nBins = hist_Old->GetNbinsX();
	float iniBin = -1.;
	float endBin = 1.;
	float bins[nBins];
	bins[0] = iniBin;
	for (int i = 1; i <= nBins; i++) bins[i] = bins[0] + i *(endBin - iniBin) / nBins;

	for (int iBins = 1; iBins < hist_Old->GetNbinsX(); iBins++)
	{
		if (hist_Old->GetBinContent(iBins) == 0. || histUnc_Old->GetBinContent(iBins) == 0.)
		{
			float temp = bins[iBins - (hist_Old->GetNbinsX() - nBins)];
			for (int jArr = iBins - (hist_Old->GetNbinsX() - nBins); jArr < nBins; jArr++)
			{
				bins[jArr] = bins[jArr + 1];
			}
			bins[nBins] = temp;
			nBins = nBins - 1;
		}
	}
	if (hist_Old->GetBinContent(hist_Old->GetNbinsX()) == 0. || histUnc_Old->GetBinContent(hist_Old->GetNbinsX()) == 0.)
	{
		bins[nBins - 1] = bins[nBins];
		nBins = nBins - 1;
	}

	///////////// Re-binning Histos	////////////////
	TH1D *hist = (TH1D*) TH1DtoTH1D(hist_Old, nBins, bins, true, "_new_obs");
	TH1D *histGGH = (TH1D*) TH1DtoTH1D(histGGH_Old, nBins, bins, true, "_new_Sig");
	TH1D *hist_obs = (TH1D*) TH1DtoTH1D(hist_obs_Old, nBins, bins, true, "_new_Data");
	TH1D *histUnc = (TH1D*) TH1DtoTH1D(histUnc_Old, nBins, bins, true, "_new_Unc");

	///////////// Normalization for Histos	////////////////
	double gghNorm = xsecGGH *lumi *gr->Eval(mass, 0, "S");

	double bkgNorm = hist_obs->GetSumOfWeights();

	hist->Scale(1 / hist->GetSumOfWeights() *bkgNorm);
	histGGH->Scale(1 / histGGH->GetSumOfWeights() *gghNorm);

	histUnc->Scale(1 / histUnc->GetSumOfWeights() *bkgNorm);

	///////////// Backgound Shape Uncertainty	////////////////
	TH1D *histUncUp = (TH1D*) histUnc->Clone("histUncUp");
	TH1D *histUncDown = (TH1D*) histUnc->Clone("histUncDown");

	for (int iB = 1; iB <= nBins; ++iB)
	{

		if (histUnc->GetBinContent(iB) != 0) histUncUp->SetBinContent(iB, hist->GetBinContent(iB) *hist->GetBinContent(iB) / histUnc->GetBinContent(iB));
		else
		{
			histUncUp->SetBinContent(iB, 0.);
		}	
		if (trackCat == "_had_had") hist->SetBinError(iB, 3 * hist->GetBinError(iB));	
	}

	histUncUp->Scale(1 / histUncUp->GetSumOfWeights() *bkgNorm);
	histUncDown->Scale(1 / histUncDown->GetSumOfWeights() *bkgNorm);

	///////////// Printing Yields	////////////////
	std::cout << "Bkg  Norm : " << bkgNorm << std::endl;
	std::cout << "ggH  Norm : " << histGGH->GetSumOfWeights() << std::endl;
	std::cout << "Data      : " << hist_obs->GetSumOfWeights() << std::endl;

	///////////// Creating Datacards & Root Files for limit extraction	////////////////
	TString DirName = "DataCards/";
	TString BaseName = "haa-13TeV_2016_ma" + mass_string;

	TFile *fileInputs = new TFile(DirName + BaseName + trackCat + ".root", "recreate");
	hist_obs->Write("data_obs");
	hist->Write("bkgd");
	histGGH->Write("ggh");
	histUncUp->Write("bkgd_CMS_bdtbkg" + trackCat + "_2016Up");
	histUncDown->Write("bkgd_CMS_bdtbkg" + trackCat + "_2016Down");

	ostringstream str;
	str << DirName + BaseName + trackCat << ".txt";
	string nn = str.str();
	const char *p = nn.c_str();

	std::ofstream textFile(p);
    textFile << "imax 1   number of channels" << std::endl;
    textFile << "jmax *   number of backgrounds" << std::endl;
    textFile << "kmax *   number of nuisance parameters" << std::endl;
    textFile << "-----------------" << std::endl;
    textFile << "observation " << hist_obs->GetSumOfWeights() << std::endl;
    textFile << "-----------------" << std::endl;
    textFile << "shapes * * " << BaseName + trackCat + ".root" << "  $PROCESS    $PROCESS_$SYSTEMATIC " << std::endl;
    textFile << "-----------------" << std::endl;
    textFile << "bin";
    textFile << "             haa       haa  "   << std::endl;
    textFile << "process      ggh       bkgd " << std::endl;
    textFile << "process       0         1   " << std::endl;
    textFile << "rate      " <<
		histGGH->GetSumOfWeights() << "  " <<
		bkgNorm << std::endl;
    textFile << "-----------------------------" << std::endl;
    textFile << "CMS_lumi_2016                           lnN         1.025          -" << std::endl;
    textFile << "CMS_eff_m_2016                          lnN         1.04           -" << std::endl;
    textFile << "CMS_trkiso_2016                         lnN         1.19           -" << std::endl;
    textFile << "CMS_bdtbkg" + trackCat + "_2016        shape         -           1.00" << std::endl; 
    textFile << "QCDScale_ggH                lnN         1.046/0.933    -" << std::endl;
    textFile << "PDF_ggh                     lnN         1.032          -" << std::endl;
    textFile << "bkgNorm" + trackCat + "_2016  rateParam  haa  bkgd  1  [0.5,1.5]" << std::endl;
    textFile << "* autoMCStats 0" << std::endl;
    textFile << std::endl;

	///// Parameter Signal Uncertainties	/////

	int nParameters = 13;

	TString Parameters[13] = { "mean_V", "width_V", "sigma_V", "mean_G1", "sigma_G1", "mean_G2", "sigma_G2", "mean_L1", "sigma_L1", "sigma_G3", "mean_L3", "sigma_L3", "sigma_G4" };

	for (int iparameter = 0; iparameter < nParameters; iparameter++)	// Loop over parameters
	{

		TH1D *histGGHUncUp_Old = (TH1D*) fileGGH->Get("MVABDTOutputSel" + trackCat + "_" + Parameters[iparameter] + "_UncUpH");
		TH1D *histGGHUncDown_Old = (TH1D*) fileGGH->Get("MVABDTOutputSel" + trackCat + "_" + Parameters[iparameter] + "_UncDownH");

		histGGHUncUp_Old->Rebin(100);
		histGGHUncDown_Old->Rebin(100);

		TH1D *histGGHUncUp = (TH1D*) TH1DtoTH1D(histGGHUncUp_Old, nBins, bins, true, "_new_SigUncUp");
		TH1D *histGGHUncDown = (TH1D*) TH1DtoTH1D(histGGHUncDown_Old, nBins, bins, true, "_new_SigUncDown");

		histGGHUncUp->Scale(1 / histGGHUncUp->GetSumOfWeights() *gghNorm);
		histGGHUncDown->Scale(1 / histGGHUncDown->GetSumOfWeights() *gghNorm);

		fileInputs->cd();

		histGGHUncUp->Write("ggh_CMS_" + Parameters[iparameter] + trackCat + "_2016Up");
		histGGHUncDown->Write("ggh_CMS_" + Parameters[iparameter] + trackCat + "_2016Down");

		textFile << "CMS_" + Parameters[iparameter] + trackCat + "_2016            shape   1.00    -  " << std::endl;
	}	// End loop over parameters

        for (int iBin = 1; iBin <=hist->GetNbinsX(); iBin++)        // Loop over Bins
        {

                TH1D *histUp = (TH1D*) hist->Clone("histUp");
                TH1D *histDown = (TH1D*) hist->Clone("histDown");

                histUp->SetBinContent(iBin, histUp->GetBinContent(iBin) + 0.5 *histUp->GetBinError(iBin));
                histDown->SetBinContent(iBin, histDown->GetBinContent(iBin) - 0.5 *histDown->GetBinError(iBin));
                fileInputs->cd();

                TString binname;
                binname.Form("%d",iBin);

                //histUp->Write("bkgd_CMS_Bin_" + binname + trackCat + "_2016Up");
                //histDown->Write("bkgd_CMS_Bin_" + binname + trackCat + "_2016Down");

                //textFile << "CMS_Bin_" + binname + trackCat + "_2016            shape    -    1.00 " << std::endl;
        }	// End loop over bins

	textFile << std::endl;

	fileInputs->Close();
	file->Close();
    fileGGH->Close();
    fileAcc->Close();

	////////////////////////////////////////////////////////////////////////////////////////
	////////////////*********************END **************************////////////////////
	////////////////////////////////////////////////////////////////////////////////////////

}

void CreateInputs(double mass_ini=3.60, double mass_end=21.01)
{

	int ntrackCat = 3;
	TString trackCat[] = { "_lep_lep", "_lep_had", "_had_had" };

	InterpolateAcceptance();

	for (double mass = mass_ini; mass < mass_end; mass += 0.20)	// Loop over masses
	{
		TString mass_string = std::to_string(int(mass + 0.1)) + "p" + std::to_string(int(int(mass *10 + 0.1) % 10));

		BDTClassificationSignal(mass_string);
		BDTClassificationData(mass_string);

		for (int icat = 0; icat < ntrackCat; icat++) InputsDataCards(mass, trackCat[icat]);
	}	// End loop over masses

}
