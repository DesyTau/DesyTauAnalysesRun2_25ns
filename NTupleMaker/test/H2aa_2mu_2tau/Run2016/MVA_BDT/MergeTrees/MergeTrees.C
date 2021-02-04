#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"

void MergeTreesSignal(TString mass = "5")
{

	float MuMu_Mass, MuMu_Pt, MuMu_DR, MuMuMET_DPhi, TrkTrk_Mass, TrkTrk_Pt, TrkTrk_DR, TrkTrkMET_DPhi, MuMuTrkTrk_Mass, MuMuTrkTrk_Pt, MuMuTauTau_Mass, MuMuTauTau_Pt, MET_Pt, MuMuTrkTrkMET_Mass;

	int ntrackCat = 3;
	TString trackCat[] = { "_lep_lep", "_lep_had", "_had_had" };

	TString _Channels[9] = { "_mu_mu", "_mu_ele", "_ele_mu", "_ele_ele", "_mu_had", "_had_mu", "_ele_had", "_had_ele", "_had_had" };
	int nchannels = 0;
	int startat = 0;

	TFile *fileGGH = new TFile("../../SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-" + mass + ".root");

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Signal File	/////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl << std::endl;
	std::cout << "--- Merging Trees of Categories --> Using input file: " << fileGGH->GetName() << std::endl;

	TFile *tree_s_file = new TFile("tree_s_ma" + mass + "_file.root", "RECREATE");
	tree_s_file->cd("");

	for (int iTrk = 0; iTrk < ntrackCat; iTrk++)	////// loop over channels
	{

		if (trackCat[iTrk] == "_lep_lep")
		{
			nchannels = 4;
			startat = 0;
		}
		if (trackCat[iTrk] == "_lep_had")
		{
			nchannels = 4;
			startat = 4;
		}
		if (trackCat[iTrk] == "_had_had")
		{
			nchannels = 1;
			startat = 8;
		}

		TTree *tree_s = new TTree("tree_s" + trackCat[iTrk], "tree_s" + trackCat[iTrk]);
		tree_s->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
		tree_s->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
		tree_s->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
		tree_s->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
		tree_s->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
		tree_s->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
		tree_s->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
		tree_s->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
		tree_s->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
		tree_s->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
		tree_s->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
		tree_s->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
		tree_s->Branch("MET_Pt", &MET_Pt, "MET_Pt");
		tree_s->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");

		for (int ichann = startat; ichann < (startat + nchannels); ichann++)	////// loop over cat's inside a channel
		{

			TTree *treeS = (TTree*) fileGGH->Get("treeSel_0p2" + _Channels[ichann]);

			treeS->SetBranchAddress("MuMu_Mass", &MuMu_Mass);
			treeS->SetBranchAddress("MuMu_Pt", &MuMu_Pt);
			treeS->SetBranchAddress("MuMu_DR", &MuMu_DR);
			treeS->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi);
			treeS->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass);
			treeS->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt);
			treeS->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR);
			treeS->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
			treeS->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
			treeS->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
			treeS->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass);
			treeS->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt);
			treeS->SetBranchAddress("MET_Pt", &MET_Pt);
			treeS->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);

			int entriesS = treeS->GetEntries();

			for (int i = 0; i < entriesS; i++)	////// loop over tree entries
			{
				treeS->GetEntry(i);
				tree_s->Fill();
			}	////// end loop over tree entries

		}	////// end loop over cat's inside a channel

		tree_s->Write();
	}	////// end loop over channels

	tree_s_file->Close();

}

void MergeTreesBkgd()
{

	float MuMu_Mass, MuMu_Pt, MuMu_DR, MuMuMET_DPhi, TrkTrk_Mass, TrkTrk_Pt, TrkTrk_DR, TrkTrkMET_DPhi, MuMuTrkTrk_Mass, MuMuTrkTrk_Pt, MuMuTauTau_Mass, MuMuTauTau_Pt, MET_Pt, MuMuTrkTrkMET_Mass;

	int ntrackCat = 3;
	TString trackCat[] = { "_lep_lep", "_lep_had", "_had_had" };

	TString _Channels[9] = { "_mu_mu", "_mu_ele", "_ele_mu", "_ele_ele", "_mu_had", "_had_mu", "_ele_had", "_had_ele", "_had_had" };
	int nchannels = 0;
	int startat = 0;

	TFile *file = new TFile("../../SingleMuon_Run2016.root");

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Background File	/////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl << std::endl;
	std::cout << "--- Merging Trees of Categories --> Using input file: " << file->GetName() << std::endl;
	std::cout << std::endl << std::endl;

	TFile *tree_b_file = new TFile("tree_b_file.root", "RECREATE");
	tree_b_file->cd("");

	for (int iTrk = 0; iTrk < ntrackCat; iTrk++)	////// loop over channels
	{

		if (trackCat[iTrk] == "_lep_lep")
		{
			nchannels = 4;
			startat = 0;
		}
		if (trackCat[iTrk] == "_lep_had")
		{
			nchannels = 4;
			startat = 4;
		}
		if (trackCat[iTrk] == "_had_had")
		{
			nchannels = 1;
			startat = 8;
		}

		TTree *tree_b = new TTree("tree_b" + trackCat[iTrk], "tree_b" + trackCat[iTrk]);
		tree_b->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
		tree_b->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
		tree_b->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
		tree_b->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
		tree_b->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
		tree_b->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
		tree_b->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
		tree_b->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
		tree_b->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
		tree_b->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
		tree_b->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
		tree_b->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
		tree_b->Branch("MET_Pt", &MET_Pt, "MET_Pt");
		tree_b->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");

		for (int ichann = startat; ichann < (startat + nchannels); ichann++)	////// loop over cat's inside a channel
		{

			TTree *treeB = (TTree*) file->Get("tree_NNNN_0p2" + _Channels[ichann]);

			treeB->SetBranchAddress("MuMu_Mass", &MuMu_Mass);
			treeB->SetBranchAddress("MuMu_Pt", &MuMu_Pt);
			treeB->SetBranchAddress("MuMu_DR", &MuMu_DR);
			treeB->SetBranchAddress("MuMuMET_DPhi", &MuMuMET_DPhi);
			treeB->SetBranchAddress("TrkTrk_Mass", &TrkTrk_Mass);
			treeB->SetBranchAddress("TrkTrk_Pt", &TrkTrk_Pt);
			treeB->SetBranchAddress("TrkTrk_DR", &TrkTrk_DR);
			treeB->SetBranchAddress("TrkTrkMET_DPhi", &TrkTrkMET_DPhi);
			treeB->SetBranchAddress("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass);
			treeB->SetBranchAddress("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt);
			treeB->SetBranchAddress("MuMuTauTau_Mass", &MuMuTauTau_Mass);
			treeB->SetBranchAddress("MuMuTauTau_Pt", &MuMuTauTau_Pt);
			treeB->SetBranchAddress("MET_Pt", &MET_Pt);
			treeB->SetBranchAddress("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass);

			int entriesB = treeB->GetEntries();	////// loop over tree entries

			for (int i = 0; i < entriesB; i++)
			{
				treeB->GetEntry(i);
				tree_b->Fill();
			}	////// end loop over tree entries

		}	////// end loop over cat's inside a channel

		tree_b->Write();
	}	////// end loop over channels

	tree_b_file->Close();

}

void MergeAll()
{

	//Mass points
	int nSamples = 19;

	TString mpoints[19] = { "3p6", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21" };

	for (int isample = 0; isample < nSamples; isample++) MergeTreesSignal(mpoints[isample]);

	MergeTreesBkgd();

}
