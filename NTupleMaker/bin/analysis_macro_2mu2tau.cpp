#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <map>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/src/Config.cc"
#include "TRandom.h"
#include "TRandom3.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TSystem.h"
using namespace std;
double PtoEta(double Px, double Py, double Pz)
{

	double P = TMath::Sqrt(Px *Px + Py *Py + Pz *Pz);
	double cosQ = Pz / P;
	double Q = TMath::ACos(cosQ);
	double Eta = -TMath::Log(TMath::Tan(0.5 *Q));
	return Eta;

}

double dPhiFrom2P(double Px1, double Py1,
	double Px2, double Py2)
{

	double prod = Px1 *Px2 + Py1 * Py2;
	double mod1 = TMath::Sqrt(Px1 *Px1 + Py1 *Py1);
	double mod2 = TMath::Sqrt(Px2 *Px2 + Py2 *Py2);

	double cosDPhi = prod / (mod1 *mod2);

	return TMath::ACos(cosDPhi);

}

double deltaR(double Eta1, double Phi1,
	double Eta2, double Phi2)
{

	double Px1 = TMath::Cos(Phi1);
	double Py1 = TMath::Sin(Phi1);

	double Px2 = TMath::Cos(Phi2);
	double Py2 = TMath::Sin(Phi2);

	double dPhi = dPhiFrom2P(Px1, Py1, Px2, Py2);
	double dEta = Eta1 - Eta2;

	double dR = TMath::Sqrt(dPhi *dPhi + dEta *dEta);

	return dR;

}

double Function_1(double _alfa, double _beta,
	TVector3 _P_tau1, TVector3 _P_tau2, double _mumu_mass)
{

	double f1 = 2 * 1.77682 * 1.77682 - _mumu_mass *_mumu_mass + 2 *(TMath::Sqrt(1.77682 * 1.77682 + _alfa *_alfa *_P_tau1.Mag2()) *TMath::Sqrt(1.77682 * 1.77682 + _beta *_beta *_P_tau2.Mag2()) - _alfa *_beta *_P_tau1.Dot(_P_tau2));

	return f1;

}

double Function_2(double _alfa, double _beta,
	TVector3 _P_tau1, TVector3 _P_tau2, TVector3 _P_MET)
{

	double f2 = (_alfa - 1) *(_alfa - 1) *_P_tau1.Mag2() + (_beta - 1) *(_beta - 1) *_P_tau2.Mag2() + 2 *(_alfa - 1) *(_beta - 1) *_P_tau1.Dot(_P_tau2) - _P_MET.Mag2();

	return f2;

}

double MuonMomentumScale(double eta) {
  
  double absEta = abs(eta);
  double scale = 0.0;
  if (absEta<1.2)
    scale = 0.004;
  else if (absEta>=1.2&&absEta<2.1)
    scale = 0.009;
  else 
    scale = 0.027;

  return scale;
}

// obsolete parameterization
double MuLegEfficiency(double pt, double eta, double ptThres, double etaThres)
{

	double absEta = fabs(eta);

	if (absEta > etaThres) return 0;

	double effEtaLt0p9 = 0.932;
	double effEta0p9to1p2 = 0.922;
	double effEtaGt1p2 = 0.950;

	double eff = 1.0;
	double ptThresLow = ptThres - 1.0;
	double ptThresHigh = ptThres + 1.0;
	if (ptThres > 30)
	{
		ptThresLow = ptThres - 3.0;
		ptThresHigh = ptThres + 3.0;
		effEtaLt0p9 = 0.919;
		effEta0p9to1p2 = 0.841;
		effEtaGt1p2 = 0.866;
	}
	else if (ptThres > 16)
	{
		effEtaLt0p9 = 0.931;
		effEta0p9to1p2 = 0.919;
		effEtaGt1p2 = 0.926;
	}

	if (pt < ptThresLow)
	{
		eff = 0;
	}
	else if (pt >= ptThresLow && pt < ptThresHigh)
	{
		if (absEta < 0.9)
			eff = 0.5 * effEtaLt0p9;
		else if (absEta >= 0.9 && absEta < 1.2)
			eff = 0.5 * effEta0p9to1p2;
		else
			eff = 0.5 * effEtaGt1p2;
	}
	else
	{
		if (absEta < 0.9)
			eff = effEtaLt0p9;
		else if (absEta >= 0.9 && absEta < 1.2)
			eff = effEta0p9to1p2;
		else
			eff = effEtaGt1p2;
	}

	return eff;

}

const float MuMass = 0.105658367;
const float PionMass = 0.13957;

int main(int argc, char *argv[])
{

	if (argc < 2)
	{
		//  std::cout << "Usage of the program : Hto4TausAnalysis[file_list]" << std::endl;
		//  std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
		exit(1);
	}

	// ****configuration
	Config cfg(argv[1]);

	const bool isData = cfg.get<bool> ("IsData");
	const int year = cfg.get<int> ("year");
	const bool applyHiggsPtWeight = cfg.get<bool> ("ApplyHiggsPtWeight");

	// kinematic cuts on muons
	const float ptGoodMuonCut = cfg.get<float> ("ptGoodMuonCut");
	const float ptIsoMuonCut = cfg.get<float> ("ptIsoMuonCut");
	const float etaMuonCut = cfg.get<float> ("etaMuonCut");
	const float dxyMuonCut = cfg.get<float> ("dxyMuonCut");
	const float dzMuonCut = cfg.get<float> ("dzMuonCut");

	// kinematic cuts on Dimuon
	const float ptSumDiMuCut = cfg.get<float> ("ptSumDiMuCut");
	const float ptSumDiTrkCut = cfg.get<float> ("ptSumDiTrkCut");
	const float massDiMuCut = cfg.get<float> ("massDiMuCut");
	const float dxyDiMuCut = cfg.get<float> ("dxyDiMuCut");
	const float dzDiMuCut = cfg.get<float> ("dzDiMuCut");

	// jets
	const float btagCut = cfg.get<float> ("btagCut");

	// topological cuts
	const float dRMuonsCut = cfg.get<float> ("dRMuonsCut");

	// track selection
	const float ptTrkLooseCut = cfg.get<float> ("ptTrkLooseCut");
	const float ptTrkCut = cfg.get<float> ("ptTrkCut");
	const float ptTrkSoftCut = cfg.get<float>("ptTrkSoftCut");
	const float etaTrkCut = cfg.get<float> ("etaTrkCut");
	const float dxyTrkLooseCut = cfg.get<float> ("dxyTrkLooseCut");
	const float dxyTrkCut = cfg.get<float> ("dxyTrkCut");
	const float dzTrkLooseCut = cfg.get<float> ("dzTrkLooseCut");
	const float dzTrkCut = cfg.get<float> ("dzTrkCut");

	// Visible Mass Cut
	const float visibleMassCut = cfg.get<float> ("visibleMassCut");
	const float TotalMassWindow = cfg.get<float> ("TotalMassWindow");

	const bool applyGoodRunSelection = cfg.get<bool> ("ApplyGoodRunSelection");
	const string jsonFile = cfg.get<string> ("jsonFile");

	// trigger
	const bool applyTriggerMatch = cfg.get<bool> ("ApplyTriggerMatch");

	// IsoTkMu24 (2016), IsoMu27 (2017), IsoMu24 (2018) 
	const string isomuTriggerName = cfg.get<string> ("IsoMuTriggerName");
	const string isomuFilterName = cfg.get<string> ("IsoMuFilterName");

	// trigger matching
	const float DRTrigMatch = cfg.get<float> ("DRTrigMatch");
	const unsigned int numberOfMuons = cfg.get < unsigned int > ("NumberOfMuons");

	TString IsoMuTriggerName(isomuTriggerName);
	TString IsoMuFilterName(isomuFilterName);

	const string pileUpDataFile = cfg.get<string> ("PileUpDataFileName");
	const string pileUpMCFile = cfg.get<string> ("PileUpMCFileName");

	TString PileUpDataFile(pileUpDataFile);
	TString PileUpMCFile(pileUpMCFile);

	const string MuIsoTriggerFile = cfg.get<string> ("MuonTriggerEffFile");
	const string CorrectionWSFileName = cfg.get<string> ("CorrectionWorkspaceFileName");

	// Higgs pt reweighting
	const string higgsPtFileName = cfg.get<string> ("HiggsPtFileName");
	TString HiggsPtFileName(higgsPtFileName);
	const bool isVH = cfg.get<bool> ("IsVH");
	// **********end of configuration *******************

	float a0_MuMu_DR = 0.936;
	float a1_MuMu_DR = 0.135;
	float a0_TrkTrk_DR = 0.918;
	float a1_TrkTrk_DR = 0.359;

	if (year==2017) {
	  a0_MuMu_DR = 0.919;
	  a1_MuMu_DR = 0.164;
	  a0_TrkTrk_DR = 0.920;
	  a1_TrkTrk_DR = 0.331;
	}

	if (year==2016) {
	  a0_MuMu_DR = 0.932;
	  a1_MuMu_DR = 0.137;
	  a0_TrkTrk_DR = 0.929;
	  a1_TrkTrk_DR = 0.287;
	}

	std::ifstream fileList(argv[2]);

	// event info
	ULong64_t event_nr;
	unsigned int event_run;
	unsigned int event_luminosityblock;

	// tracks
	UInt_t track_count;
	int track_ID[1000];
	float track_px[1000];
	float track_py[1000];
	float track_pz[1000];
	float track_pt[1000];
	float track_eta[1000];
	float track_phi[1000];
	float track_charge[1000];
	float track_mass[1000];
	float track_dxy[1000];
	float track_dxyerr[1000];
	float track_dz[1000];
	float track_dzerr[1000];
	bool track_highPurity[1000];
	// muons
	UInt_t muon_count;
	UInt_t muon_nMuonStations[1000];
	UInt_t muon_nMuonHits[1000];
	UInt_t muon_nPixelHits[1000];
	UInt_t muon_nTrackerHits[1000];
	float muon_px[1000];
	float muon_py[1000];
	float muon_pz[1000];
	float muon_pt[1000];
	float muon_eta[1000];
	float muon_phi[1000];
	float muon_pterror[1000];
	float muon_chi2[1000];
	float muon_ndof[1000];
	float muon_charge[1000];
	float muon_dxy[1000];
	float muon_dxyerr[1000];
	float muon_dz[1000];
	float muon_dzerr[1000];
	float muon_chargedHadIso[1000];
	float muon_neutralHadIso[1000];
	float muon_photonIso[1000];
	float muon_puIso[1000];
	bool muon_isPF[1000];
	bool muon_isGlobal[1000];
	bool muon_isTracker[1000];
	bool muon_isTight[1000];
	bool muon_isLoose[1000];
	bool muon_isMedium[1000];
	bool muon_isICHEP[1000];

	UInt_t genparticles_count;
	Float_t genparticles_e[1000];
	Float_t genparticles_px[1000];
	Float_t genparticles_py[1000];
	Float_t genparticles_pz[1000];
	Int_t genparticles_pdgid[1000];
	Int_t genparticles_status[1000];
	UInt_t genparticles_info[1000];

	float genweight;

	float metx;
	float mety;
	float met;
	float metphi;

	float metx_JetEnUp;
	float mety_JetEnUp;

	float metx_JetEnDown;
	float mety_JetEnDown;

	float metx_UnclEnUp;
	float mety_UnclEnUp;

	float metx_UnclEnDown;
	float mety_UnclEnDown;


	// Trigger
	unsigned int trigobject_count;
	float trigobject_px[1000];
	float trigobject_py[1000];
	float trigobject_pz[1000];
	float trigobject_pt[1000];
	float trigobject_eta[1000];
	float trigobject_phi[1000];
	bool trigobject_filters[1000][200];

	// Jets
	unsigned int pfjet_count;
	float pfjet_e[200];
	float pfjet_px[200];
	float pfjet_py[200];
	float pfjet_pz[200];
	float pfjet_pt[200];
	float pfjet_eta[200];
	float pfjet_phi[200];
	int pfjet_flavour[200];
	float pfjet_btag[200][10];
	Bool_t pfjet_pu_jet_fullId_loose[100];
	Bool_t pfjet_pu_jet_fullId_medium[100];
	Bool_t pfjet_pu_jet_fullId_tight[100];

	float numtruepileupinteractions;

	//unsigned int iLeadingPosTrig = 0;
	//vector<bool> trigobject_filter; trigobject_filter.clear();

	std::map<std::string, int> *hltriggerresults = new std::map<std::string, int> ();
	std::map<std::string, int> *hltriggerprescales = new std::map<std::string, int> ();
	std::vector<std::string > *hltfilters = new std::vector<std::string > ();

	std::string rootFileName(argv[2]);

	std::string chainName("makeroottree/AC1B");
	std::string initNtupleName("initroottree/AC1B");
	TString TStrName(rootFileName);
	//std::cout <<TStrName <<std::endl;
	if (TStrName.Contains("Signal"))
	{
		//  std::cout << "=============================" << std::endl;
		//  std::cout << "=== Running on Signal MC ====" << std::endl;
		//  std::cout << "=============================" << std::endl;
		//  std::cout << std::endl;
	}

	TString FullName = TStrName;

	TFile *file = new TFile(FullName + TString(".root"), "recreate");

	file->cd("");

	string cmsswBase = (getenv("CMSSW_BASE"));
	///////////////////////////////////////////////////
	/////////////////////HISTOGRAMS////////////////////
	//////////////// Muons-PU-Trigger//////////////////
	///////////////////////////////////////////////////

	TH1D *muonCountH = new TH1D("muonCountH", "", 11, -0.5, 10.5);
	TH1D *nGoodMuonsH = new TH1D("nGoodMuonsH", "", 11, -0.5, 10.5);
	TH1D *puWeightH = new TH1D("puWeightH", "", 250, 0, 5);

	///////////////////////////////////////////////////
	/////// histograms after triple-muon selection/////
	///////////////////////////////////////////////////

	//////Categorie Denominations//////////////////////
	TString DRIntervalString[1] = { "_0p2" };
	int nDRInterval = 1;
	double DRInterval[1] = { 0.2 };

	TString TrackCatString[3] = { "_mu", "_ele", "_had" };
	int nTrackCat = 3;

	TH1D *triggerWeightH = new TH1D("triggerWeightH", "", 100, 0, 2);
	TH1D *idIsoWeightH = new TH1D("idIsoWeightH", "", 100, 0, 2);
	TH1D *HiggsPtWeightH = new TH1D("HiggsPtWeightH", "", 100, 0, 2);

	//////////////////Tracks///////////////////////////
	TH1D *nTracksAmumu1candidateH[nDRInterval];
	TH1D *nTracksAmumu2candidateH[nDRInterval];
	TH1D *nTracksAtautau1candidateH[nDRInterval];
	TH1D *nTracksAtautau2candidateH[nDRInterval];

	TH1D *nSoftTracksAmumu1candidateH[nDRInterval];
	TH1D *nSoftTracksAmumu2candidateH[nDRInterval];
	TH1D *nSoftTracksAtautau1candidateH[nDRInterval];
	TH1D *nSoftTracksAtautau2candidateH[nDRInterval];

	TH1D *nSignalTracksAmumu1candidateH[nDRInterval];
	TH1D *nSignalTracksAmumu2candidateH[nDRInterval];
	TH1D *nSignalTracksAtautau1candidateH[nDRInterval];
	TH1D *nSignalTracksAtautau2candidateH[nDRInterval];

	///////////////////////////////////////////////////
	/// Histograms For Diferent Regions To Be Used	////
	//~hist. w/o sufix are a combination of all Cat.~//
	///////////////////////////////////////////////////
	///Varables for Dataset
	float MuMu_Mass, MuMu_Mass_muScaleUp, MuMu_Mass_muScaleDown, MuMu_Pt, MuMu_Pt_muScaleUp, MuMu_Pt_muScaleDown, MuMu_DR, MuMuMET_DPhi, TrkTrk_Mass, TrkTrk_Pt, TrkTrk_DR, TrkTrkMET_DPhi, MuMuTrkTrk_Mass, MuMuTrkTrk_Pt, MuMuTauTau_Mass, MuMuTauTau_Pt, MET_Pt, MuMuTrkTrkMET_Mass, MuMuTrkTrkMET_Mass_muScaleUp, MuMuTrkTrkMET_Mass_muScaleDown, MuMuTrkTrkMET_Mass_JetEnUp, MuMuTrkTrkMET_Mass_JetEnDown, MuMuTrkTrkMET_Mass_UnclEnUp, MuMuTrkTrkMET_Mass_UnclEnDown;
	float Eventweight;
	float Correctionweight;
	float TrkTrk_DR_weight_Up;
	float TrkTrk_DR_weight_Down;
	float MuMu_DR_weight_Up;
	float MuMu_DR_weight_Down;


	////// Control Regions for Bkgd Estimation	////////
	TTree *tree_NNNN[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00NN[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00SemiIso[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_SoftIso[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00SoftIso[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_VerySoftIso[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00VerySoftIso[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_NN00[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00_0SemiIsoorSemiIso0[nDRInterval][nTrackCat][nTrackCat];

	TTree *tree_00_0NorN0[nDRInterval][nTrackCat][nTrackCat];

	////////////////// Signal Region	//////////////////
	TTree *treeSel[nDRInterval][nTrackCat][nTrackCat];

	/////////////////// Counters	/////////////////////
	TH1D *counter_MuonKinematicsH = new TH1D("counter_MuonKinematicsH", "", 1, 0., 2.);

	TH1D *counter_FinalEventsCatH[nDRInterval][nTrackCat][nTrackCat];

	TH1D *counter_InputEventsH = new TH1D("counter_InputEventsH", "", 1, 0., 2.);
	TH1D *counter_FinalEventsH = new TH1D("counter_FinalEventsH", "", 1, 0., 2.);
	TH1D *counter_MuonSizeGTE2H = new TH1D("counter_MuonSizeGTE2H", "", 1, 0., 2.);

	TH1D *histWeightsH = new TH1D("histWeightsH", "", 1, 0., 2.);
	TH1D *histWeightsSingleMuH = new TH1D("histWeightsSingleMuH", "", 1, 0., 2.);
	TH1D *histWeightsTripleMuH = new TH1D("histWeightsTripleMuH", "", 1, 0., 2.);
	TH1D *histWeightsDoubleMuSSH = new TH1D("histWeightsDoubleMuSSH", "", 1, 0., 2.);
	TH1D *histWeightsAllTriggersH = new TH1D("histWeightsAllTriggersH", "", 1, 0., 2.);

	TH1D *BjetsMultipH = new TH1D("BjetsMultipH", "", 21, -0.5, 20.5);
	TH1D *BTagEventWeightH = new TH1D("BTagEventWeightH", "", 400, -2., 2.);
	TH1D *histWeightsBTagH = new TH1D("histWeightsBTagH", "", 1., 0., 2.);
	TH1D *histWeightsZeroBTagH = new TH1D("histWeightsZeroBTagH", "", 1., 0., 2.);

	//////////// Monte Carlo information	/////////////
	TH1D *deltaRMuonPionH = new TH1D("deltaRMuonPionH", "", 200, 0, 2);
	TH1D *pionPtH = new TH1D("pionPtH", "", 1000, 0, 100);
	TH1D *deltaRMuonMuonH = new TH1D("deltaRMuonMuonH", "", 200, 0, 2);
	TH1D *deltaRSSMuonsH = new TH1D("deltaRSSMuonsH", "", 200, 0, 4);
	TH2D *deltaRMuonMuonvsdeltaRMuonPionH = new TH2D("deltaRMuonMuonvsdeltaRMuonPionH", "", 200, 0, 2, 200, 0, 2.);

	//////Kinematics of 3 Muons + 1 Prong////////////////////////
	TH1D *ptAmumu1candidateH = new TH1D("ptAmumu1candidateH", "", 400, 0, 400);
	TH1D *ptAmumu2candidateH = new TH1D("ptAmumu2candidateH", "", 400, 0, 400);
	TH1D *ptAtautau1candidateH = new TH1D("ptAtautau1candidateH", "", 400, 0, 400);
	TH1D *ptAtautau2candidateH = new TH1D("ptAtautau2candidateH", "", 400, 0, 400);
	TH1D *etaAmumu1candidateH = new TH1D("etaAmumu1candidateH", "", 48, -2.4, 2.4);
	TH1D *etaAmumu2candidateH = new TH1D("etaAmumu2candidateH", "", 48, -2.4, 2.4);
	TH1D *etaAtautau1candidateH = new TH1D("etaAtautau1candidateH", "", 48, -2.4, 2.4);
	TH1D *etaAtautau2candidateH = new TH1D("etaAtautau2candidateH", "", 48, -2.4, 2.4);

	TH1D *dxyAmumu1candidateH = new TH1D("dxyAmumu1candidateH", "", 200, -0.5, 0.5);
	TH1D *dxyAmumu2candidateH = new TH1D("dxyAmumu2candidateH", "", 200, -0.5, 0.5);
	TH1D *dxyAtautau1candidateH = new TH1D("dxyAtautau1candidateH", "", 200, -0.5, 0.5);
	TH1D *dxyAtautau2candidateH = new TH1D("dxyAtautau2candidateH", "", 200, -0.5, 0.5);
	TH1D *dzAmumu1candidateH = new TH1D("dzAmumu1candidateH", "", 200, -0.5, 0.5);
	TH1D *dzAmumu2candidateH = new TH1D("dzAmumu2candidateH", "", 200, -0.5, 0.5);
	TH1D *dzAtautau1candidateH = new TH1D("dzAtautau1candidateH", "", 200, -0.5, 0.5);
	TH1D *dzAtautau2candidateH = new TH1D("dzAtautau2candidateH", "", 200, -0.5, 0.5);

	TH1D *dimuonMassH = new TH1D("dimuonMassH", "", 2000, 0., 20.);
	TH1D *ditrackMassH = new TH1D("ditrackMassH", "", 2000, 0., 20.);
	TH1D *MetH = new TH1D("MetH", "", 400, 0., 400.);

	for (int iInt = 0; iInt < nDRInterval; ++iInt)	/////starting loop over all Intervals (Categories) used for Histograms
	{

		//============================================================================================================//

		nTracksAmumu1candidateH[iInt] = new TH1D("nTracksAmumu1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nTracksAmumu2candidateH[iInt] = new TH1D("nTracksAmumu2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nTracksAtautau1candidateH[iInt] = new TH1D("nTracksAtautau1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nTracksAtautau2candidateH[iInt] = new TH1D("nTracksAtautau2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);

		nSoftTracksAmumu1candidateH[iInt] = new TH1D("nSoftTracksAmumu1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSoftTracksAmumu2candidateH[iInt] = new TH1D("nSoftTracksAmumu2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSoftTracksAtautau1candidateH[iInt] = new TH1D("nSoftTracksAtautau1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSoftTracksAtautau2candidateH[iInt] = new TH1D("nSoftTracksAtautau2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);

		nSignalTracksAmumu1candidateH[iInt] = new TH1D("nSignalTracksAmumu1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSignalTracksAmumu2candidateH[iInt] = new TH1D("nSignalTracksAmumu2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSignalTracksAtautau1candidateH[iInt] = new TH1D("nSignalTracksAtautau1candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);
		nSignalTracksAtautau2candidateH[iInt] = new TH1D("nSignalTracksAtautau2candidateH" + DRIntervalString[iInt], "", 21, -0.5, 20.5);

		//============================================================================================================//
		//============================================================================================================//

		for (int jInt = 0; jInt < nTrackCat; ++jInt)	/////starting loop over track categories used for Histograms
		{
			//============================================================================================================//
			for (int kInt = 0; kInt < nTrackCat; ++kInt)	/////starting loop over track categories used for Histograms
			{
			 	//============================================================================================================//

				////////////////// NNNN Region	//////////////////
				tree_NNNN[iInt][jInt][kInt] = new TTree("tree_NNNN" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_NNNN" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_NNNN[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_NNNN[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_NNNN[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_NNNN[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_NNNN[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_NNNN[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_NNNN[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// 00NN Region	//////////////////
				tree_00NN[iInt][jInt][kInt] = new TTree("tree_00NN" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00NN" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00NN[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00NN[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00NN[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00NN[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00NN[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00NN[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00NN[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00NN[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// 00SemiIso Region	//////////////////
				tree_00SemiIso[iInt][jInt][kInt] = new TTree("tree_00SemiIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00SemiIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00SemiIso[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// SoftIso Region	//////////////////
				tree_SoftIso[iInt][jInt][kInt] = new TTree("tree_SoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_SoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_SoftIso[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_SoftIso[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// VerySoftIso Region	//////////////////
				tree_VerySoftIso[iInt][jInt][kInt] = new TTree("tree_VerySoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_VerySoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("Correctionweight",&Correctionweight,"Correctionweight");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR_weight_Up",&TrkTrk_DR_weight_Up,"TrkTrk_DR_weight_Up");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR_weight_Down",&TrkTrk_DR_weight_Down,"TrkTrk_DR_weight_Down");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_DR_weight_Up",&MuMu_DR_weight_Up,"MuMu_DR_weight_Up");
				tree_VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_DR_weight_Down",&MuMu_DR_weight_Down,"MuMu_DR_weight_Down");

				////////////////// 00SoftIso Region	//////////////////
				tree_00SoftIso[iInt][jInt][kInt] = new TTree("tree_00SoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00SoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00SoftIso[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// 00VerySoftIso Region	//////////////////
				tree_00VerySoftIso[iInt][jInt][kInt] = new TTree("tree_00VerySoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00VerySoftIso" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00VerySoftIso[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// NN00 Region	//////////////////
				tree_NN00[iInt][jInt][kInt] = new TTree("tree_NN00" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_NN00" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_NN00[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_NN00[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_NN00[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_NN00[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_NN00[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_NN00[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_NN00[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_NN00[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// 00_0SemiIsoorSemiIso0 Region	//////////////////
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt] = new TTree("tree_00_0SemiIsoorSemiIso0" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00_0SemiIsoorSemiIso0" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00_0SemiIsoorSemiIso0[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// 00_0NorN0 Region	//////////////////
				tree_00_0NorN0[iInt][jInt][kInt] = new TTree("tree_00_0NorN0" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "tree_00_0NorN0" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				tree_00_0NorN0[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				////////////////// Signal Region	//////////////////
				treeSel[iInt][jInt][kInt] = new TTree("treeSel" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "treeSel" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt]);
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Mass_muScaleUp", &MuMu_Mass_muScaleUp, "MuMu_Mass_muScaleUp");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Mass_muScaleDown", &MuMu_Mass_muScaleDown, "MuMu_Mass_muScaleDown");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Pt_muScaleUp", &MuMu_Pt_muScaleUp, "MuMu_Pt_muScaleUp");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_Pt_muScaleDown", &MuMu_Pt_muScaleDown, "MuMu_Pt_muScaleDown");
				treeSel[iInt][jInt][kInt]->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
				treeSel[iInt][jInt][kInt]->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
				treeSel[iInt][jInt][kInt]->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
				treeSel[iInt][jInt][kInt]->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
				treeSel[iInt][jInt][kInt]->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
				treeSel[iInt][jInt][kInt]->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Mass", &MuMuTrkTrk_Mass, "MuMuTrkTrk_Mass");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrk_Pt", &MuMuTrkTrk_Pt, "MuMuTrkTrk_Pt");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTauTau_Mass", &MuMuTauTau_Mass, "MuMuTauTau_Mass");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTauTau_Pt", &MuMuTauTau_Pt, "MuMuTauTau_Pt");
				treeSel[iInt][jInt][kInt]->Branch("MET_Pt", &MET_Pt, "MET_Pt");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass", &MuMuTrkTrkMET_Mass, "MuMuTrkTrkMET_Mass");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_muScaleUp", &MuMuTrkTrkMET_Mass_muScaleUp, "MuMuTrkTrkMET_Mass_muScaleUp");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_muScaleDown", &MuMuTrkTrkMET_Mass_muScaleDown, "MuMuTrkTrkMET_Mass_muScaleDown");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_JetEnUp", &MuMuTrkTrkMET_Mass_JetEnUp, "MuMuTrkTrkMET_Mass_JetEnUp");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_JetEnDown", &MuMuTrkTrkMET_Mass_JetEnDown, "MuMuTrkTrkMET_Mass_JetEnDown");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_UnclEnUp", &MuMuTrkTrkMET_Mass_UnclEnUp, "MuMuTrkTrkMET_Mass_UnclEnUp");
				treeSel[iInt][jInt][kInt]->Branch("MuMuTrkTrkMET_Mass_UnclEnDown", &MuMuTrkTrkMET_Mass_UnclEnDown, "MuMuTrkTrkMET_Mass_UnclEnDown");
				treeSel[iInt][jInt][kInt]->Branch("Eventweight", &Eventweight, "Eventweight");

				//============================================================================================================//

				counter_FinalEventsCatH[iInt][jInt][kInt] = new TH1D("counter_FinalEventsH" + DRIntervalString[iInt] + TrackCatString[jInt] + TrackCatString[kInt], "", 1, 0., 2.);
			}	/////end loop over all track categories used for Histograms
		}	/////end loop over all track categories used for Histograms

	}	/////end loop over all Intervals (Categories) used for Histograms

	float higgsPt;
	float higgsSMPt;
	bool isWplus;
	bool isWminus;
	bool isZ;
	float genWeight;

	////// Generator Level Distributions	////////
	TTree *tree_GenLev = new TTree("tree_GenLev", "tree_GenLev");
	tree_GenLev->Branch("MuMu_Mass", &MuMu_Mass, "MuMu_Mass");
	tree_GenLev->Branch("MuMu_Pt", &MuMu_Pt, "MuMu_Pt");
	tree_GenLev->Branch("MuMu_DR", &MuMu_DR, "MuMu_DR");
	tree_GenLev->Branch("MuMuMET_DPhi", &MuMuMET_DPhi, "MuMuMET_DPhi");
	tree_GenLev->Branch("TrkTrk_Mass", &TrkTrk_Mass, "TrkTrk_Mass");
	tree_GenLev->Branch("TrkTrk_Pt", &TrkTrk_Pt, "TrkTrk_Pt");
	tree_GenLev->Branch("TrkTrk_DR", &TrkTrk_DR, "TrkTrk_DR");
	tree_GenLev->Branch("TrkTrkMET_DPhi", &TrkTrkMET_DPhi, "TrkTrkMET_DPhi");
	tree_GenLev->Branch("MET_Pt", &MET_Pt, "MET_Pt");

	TTree *higgsTree = new TTree("higgsTree", "");
	higgsTree->Branch("HiggsPt", &higgsPt, "HiggsPt/F");
	higgsTree->Branch("isWplus", &isWplus, "isWplus/O");
	higgsTree->Branch("isWminus", &isWminus, "isWminus/O");
	higgsTree->Branch("isZ", &isZ, "isZ/O");
	higgsTree->Branch("genweight", &genWeight, "genweight/F");

	TTree *higgsSMTree = new TTree("higgsSMTree", "");
	higgsSMTree->Branch("HiggsSMPt", &higgsSMPt, "HiggsSMPt/F");
	higgsSMTree->Branch("genweight", &genWeight, "genweight/F");

	TH1D *triggerH = new TH1D("triggerH", "", 2, -0.5, 1.5);
	TH1D *triggerSigH = new TH1D("triggerSigH", "", 2, -0.5, 1.5);
	TH1D *nmuonsH = new TH1D("nmuonsH", "", 10, -0.5, 9.5);

	// Run-lumi selector
	string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
	std::vector<Period> periods;
	if (isData)
	{
		// read the good runs
		std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios:: in);
		if (inputFileStream.fail())
		{
			// std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
			// std::cout << "please check" << std::endl;
			// std::cout << "quitting program" << std::endl;
			exit(-1);
		}

		for (std::string s; std::getline(inputFileStream, s);)
		{
			periods.push_back(Period());
			std::stringstream ss(s);
			ss >> periods.back();
		}
	}

	// PU reweighting
	PileUp *PUofficial = new PileUp();
	TFile *filePUdistribution_data = new TFile(TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/" + PileUpDataFile, "read");
	TFile *filePUdistribution_MC = new TFile(TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/" + PileUpMCFile, "read");
	TH1D *PU_data = (TH1D*) filePUdistribution_data->Get("pileup");
	TH1D *PU_mc = (TH1D*) filePUdistribution_MC->Get("pileup");
	PUofficial->set_h_data(PU_data);
	PUofficial->set_h_MC(PU_mc);

	// Higgs reweighting
	TString fullpath_HiggsPtFile = TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/" + HiggsPtFileName;
	TFile *higgsPtFile = NULL;
	TH1D *higgsPtH = NULL;
	TH1D *higgsPt_WPlusH = NULL;
	TH1D *higgsPt_WMinusH = NULL;
	TH1D *higgsPt_ZH = NULL;
	if (applyHiggsPtWeight)
	{
		std::cout << "ApplyHiggsPtWeight = " << applyHiggsPtWeight << std::endl;
		higgsPtFile = new TFile(fullpath_HiggsPtFile);
		if (higgsPtFile->IsZombie())
		{
			std::cout << fullpath_HiggsPtFile << "  not found" << std::endl;
			exit(-1);
		}
		if (isVH)
		{
			std::cout << "IsVH = " << isVH << std::endl;
			higgsPt_WPlusH = (TH1D*) higgsPtFile->Get("kfactor_WplusH");
			higgsPt_WMinusH = (TH1D*) higgsPtFile->Get("kfactor_WminusH");
			higgsPt_ZH = (TH1D*) higgsPtFile->Get("kfactor_ZH");
			if (higgsPt_WPlusH == NULL)
			{
				std::cout << "histogram kfactor_WplusH is not found in file " << fullpath_HiggsPtFile << std::endl;
				exit(-1);
			}
			if (higgsPt_WMinusH == NULL)
			{
				std::cout << "histogram kfactor_WminusH is not found in file " << fullpath_HiggsPtFile << std::endl;
				exit(-1);
			}
			if (higgsPt_ZH == NULL)
			{
				std::cout << "histogram kfactor_ZH is not found in file " << fullpath_HiggsPtFile << std::endl;
				exit(-1);
			}
		}
		else
		{
			higgsPtH = (TH1D*) higgsPtFile->Get("kfactor");
			if (higgsPtH == NULL)
			{
				std::cout << "histogram kfactor is not found in file " << fullpath_HiggsPtFile << std::endl;
				exit(-1);
			}
		}
	}
	//  std::cout << "Higgs Pt histogram : " << higgsPtH << std::endl;

	std::string BTagCalibFileName;
	if (year == 2016) BTagCalibFileName = "DeepCSV_2016LegacySF_V1";
	else if (year == 2017) BTagCalibFileName = "DeepCSV_94XSF_V4_B_F";
	else if (year == 2018) BTagCalibFileName = "DeepCSV_102XSF_V1";
	else
	{
		std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n';
		exit(-1);
	}

	// B-Tagging reweighting
	BTagCalibration calib("DeepCSV", cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/data/" + BTagCalibFileName + ".csv");
	BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM, "central",
	{
		"up", "down" });
	BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM, "central",
	{
		"up", "down" });
	BTagCalibrationReader reader_L(BTagEntry::OP_MEDIUM, "central",
	{
		"up", "down" });
	reader_B.load(calib, BTagEntry::FLAV_B, "comb");
	reader_C.load(calib, BTagEntry::FLAV_C, "comb");
	reader_L.load(calib, BTagEntry::FLAV_UDSG, "comb");

	// Trigger efficiency
	ScaleFactor *SF_muon = new ScaleFactor();
	SF_muon->init_ScaleFactor(TString(cmsswBase) + "/src/" + TString(MuIsoTriggerFile));

	// Load CrystalBallEfficiency class
	TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
	int openSuccessful = gSystem->Load(pathToCrystalLib);
	if (openSuccessful != 0)
	{
		cout << pathToCrystalLib << " not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. " << endl;
		exit(-1);
	}

	// Correction workspace
	TString correctionsWorkspaceFileName = TString(cmsswBase) + "/src/" + TString(CorrectionWSFileName);
	TFile *correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
	if (correctionWorkSpaceFile->IsZombie())
	{
		std::cout << correctionsWorkspaceFileName << " does not exist " << std::endl;
	}
	RooWorkspace *correctionWS = (RooWorkspace*) correctionWorkSpaceFile->Get("w");
	if (correctionWS == NULL)
	{
		std::cout << "correction workspace not found " << std::endl;
	}

	TString filen;
	int iFiles = 0;
	int events = 0;
	while (fileList >> filen)
	{
		iFiles++;
		//cout << "file " << iFiles << " : " << filen << endl;

		TFile *file_ = TFile::Open(TString(filen));
		if (file_ == NULL) continue;

		TTree *tree_ = (TTree*) file_->Get(TString(chainName));
		TTree *inittree_ = (TTree*) file_->Get(TString(initNtupleName));

		if (tree_ == NULL) continue;
		if (inittree_ == NULL) continue;

		tree_->SetMaxVirtualSize(3000000);
		// event info
		tree_->SetBranchAddress("event_nr", &event_nr);
		tree_->SetBranchAddress("event_run", &event_run);
		tree_->SetBranchAddress("event_luminosityblock", &event_luminosityblock);

		// Muons
		tree_->SetBranchAddress("muon_count", &muon_count);
		tree_->SetBranchAddress("muon_nMuonStations", muon_nMuonStations);
		tree_->SetBranchAddress("muon_nMuonHits", muon_nMuonHits);
		tree_->SetBranchAddress("muon_nPixelHits", muon_nPixelHits);
		tree_->SetBranchAddress("muon_nTrackerHits", muon_nTrackerHits);
		tree_->SetBranchAddress("muon_px", muon_px);
		tree_->SetBranchAddress("muon_py", muon_py);
		tree_->SetBranchAddress("muon_pz", muon_pz);
		tree_->SetBranchAddress("muon_pt", muon_pt);
		tree_->SetBranchAddress("muon_eta", muon_eta);
		tree_->SetBranchAddress("muon_phi", muon_phi);
		tree_->SetBranchAddress("muon_pterror", muon_pterror);
		tree_->SetBranchAddress("muon_chi2", muon_chi2);
		tree_->SetBranchAddress("muon_ndof", muon_ndof);
		tree_->SetBranchAddress("muon_charge", muon_charge);
		tree_->SetBranchAddress("muon_dxy", muon_dxy);
		tree_->SetBranchAddress("muon_dxyerr", muon_dxyerr);
		tree_->SetBranchAddress("muon_dz", muon_dz);
		tree_->SetBranchAddress("muon_dzerr", muon_dzerr);
		tree_->SetBranchAddress("muon_chargedHadIso", muon_chargedHadIso);
		tree_->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso);
		tree_->SetBranchAddress("muon_photonIso", muon_photonIso);
		tree_->SetBranchAddress("muon_puIso", muon_puIso);
		tree_->SetBranchAddress("muon_isMedium", muon_isMedium);
		tree_->SetBranchAddress("muon_isICHEP", muon_isICHEP);

		// MET
		tree_->SetBranchAddress("pfmetcorr_ex", &metx);
		tree_->SetBranchAddress("pfmetcorr_ey", &mety);
		tree_->SetBranchAddress("pfmetcorr_pt", &met);
		tree_->SetBranchAddress("pfmetcorr_phi", &metphi);

		tree_->SetBranchAddress("pfmetcorr_ex_JetEnUp",&metx_JetEnUp);
		tree_->SetBranchAddress("pfmetcorr_ey_JetEnUp",&mety_JetEnUp);

		tree_->SetBranchAddress("pfmetcorr_ex_JetEnDown",&metx_JetEnDown);
		tree_->SetBranchAddress("pfmetcorr_ey_JetEnDown",&mety_JetEnDown);

		tree_->SetBranchAddress("pfmetcorr_ex_UnclusteredEnUp",&metx_UnclEnUp);
		tree_->SetBranchAddress("pfmetcorr_ey_UnclusteredEnUp",&mety_UnclEnUp);

		tree_->SetBranchAddress("pfmetcorr_ex_UnclusteredEnDown",&metx_UnclEnDown);
		tree_->SetBranchAddress("pfmetcorr_ey_UnclusteredEnDown",&mety_UnclEnDown);



		// Tracks
		tree_->SetBranchAddress("track_count", &track_count);
		tree_->SetBranchAddress("track_ID", track_ID);
		tree_->SetBranchAddress("track_px", track_px);
		tree_->SetBranchAddress("track_py", track_py);
		tree_->SetBranchAddress("track_pz", track_pz);
		tree_->SetBranchAddress("track_pt", track_pt);
		tree_->SetBranchAddress("track_eta", track_eta);
		tree_->SetBranchAddress("track_phi", track_phi);
		tree_->SetBranchAddress("track_mass", track_mass);
		tree_->SetBranchAddress("track_charge", track_charge);
		tree_->SetBranchAddress("track_dxy", track_dxy);
		tree_->SetBranchAddress("track_dxyerr", track_dxyerr);
		tree_->SetBranchAddress("track_dz", track_dz);
		tree_->SetBranchAddress("track_dzerr", track_dzerr);
		tree_->SetBranchAddress("track_highPurity", track_highPurity);

		// trigger objects
		tree_->SetBranchAddress("trigobject_count", &trigobject_count);
		tree_->SetBranchAddress("trigobject_px", trigobject_px);
		tree_->SetBranchAddress("trigobject_py", trigobject_py);
		tree_->SetBranchAddress("trigobject_pz", trigobject_pz);
		tree_->SetBranchAddress("trigobject_pt", trigobject_pt);
		tree_->SetBranchAddress("trigobject_eta", trigobject_eta);
		tree_->SetBranchAddress("trigobject_phi", trigobject_phi);
		tree_->SetBranchAddress("trigobject_filters", trigobject_filters);

		// jets
		tree_->SetBranchAddress("pfjet_count", &pfjet_count);
		tree_->SetBranchAddress("pfjet_e", pfjet_e);
		tree_->SetBranchAddress("pfjet_px", pfjet_px);
		tree_->SetBranchAddress("pfjet_py", pfjet_py);
		tree_->SetBranchAddress("pfjet_pz", pfjet_pz);
		tree_->SetBranchAddress("pfjet_pt", pfjet_pt);
		tree_->SetBranchAddress("pfjet_eta", pfjet_eta);
		tree_->SetBranchAddress("pfjet_phi", pfjet_phi);
		tree_->SetBranchAddress("pfjet_flavour", pfjet_flavour);
		tree_->SetBranchAddress("pfjet_btag", pfjet_btag);
		tree_->SetBranchAddress("pfjet_pu_jet_fullId_loose", pfjet_pu_jet_fullId_loose);
		tree_->SetBranchAddress("pfjet_pu_jet_fullId_medium", pfjet_pu_jet_fullId_medium);
		tree_->SetBranchAddress("pfjet_pu_jet_fullId_tight", pfjet_pu_jet_fullId_tight);

		// Additional trigger objects
		tree_->SetBranchAddress("run_hltfilters", &hltfilters);
		//   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
		tree_->SetBranchAddress("hltriggerresults", &hltriggerresults);
		tree_->SetBranchAddress("hltriggerprescales", &hltriggerprescales);

		if (!isData) tree_->SetBranchAddress("numtruepileupinteractions", &numtruepileupinteractions);

		if (!isData)
		{
			inittree_->SetBranchAddress("genweight", &genweight);
			tree_->SetBranchAddress("genweight", &genweight);
			tree_->SetBranchAddress("genparticles_count", &genparticles_count);
			tree_->SetBranchAddress("genparticles_e", genparticles_e);
			tree_->SetBranchAddress("genparticles_px", genparticles_px);
			tree_->SetBranchAddress("genparticles_py", genparticles_py);
			tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
			tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
			tree_->SetBranchAddress("genparticles_status", genparticles_status);
			tree_->SetBranchAddress("genparticles_info", genparticles_info);
		}

		TString sample(argv[2]);
		bool is4tau = sample.Contains("ToAA_AToTauTau");

		int numberOfCandidates = tree_->GetEntries();

		TRandom3 rand;

		for (int iCand = 0; iCand < inittree_->GetEntries(); iCand++)
		{

			inittree_->GetEntry(iCand);

			float weight = 1;
			if (!isData)
			{
				weight *= genweight;
			}
			histWeightsH->Fill(1.0, weight);
		}

		std::cout << "number of events (initroottree) = " << inittree_->GetEntries() << std::endl;
		std::cout << "number of events (makeroottree) = " << numberOfCandidates << std::endl;

		delete inittree_;

		for (int iCand = 0; iCand < numberOfCandidates; iCand++)
		{

			tree_->GetEntry(iCand);

			events++;
			if (events % 10000 == 0) cout << "   processed events : " << events << endl;

			float weight = 1;
			if (!isData)
			{
				weight *= genweight;
			}

			std::vector < unsigned int > posPion;
			posPion.clear();
			std::vector < unsigned int > negPion;
			negPion.clear();
			std::vector < unsigned int > posMuon;
			posMuon.clear();
			std::vector < unsigned int > negMuon;
			negMuon.clear();

			bool HiggsFound = false;
			unsigned int higgsIndex = 0;
			bool HiggsSMFound = false;
			unsigned int higgsSMIndex = 0;
			unsigned int nmuons = 0;

			isWplus = false;
			isWminus = false;
			isZ = false;

			if (!isData)
			{
			 	//       std::cout << "Generated particles = " << genparticles_count << std::endl;
				for (unsigned int iP = 0; iP < genparticles_count; ++iP)
				{
					if (genparticles_status[iP] == 1)
					{
					 			//&&genparticles_info[iP]==12
						if (genparticles_pdgid[iP] == 13) negMuon.push_back(iP);
						if (genparticles_pdgid[iP] == -13) posMuon.push_back(iP);
						if (genparticles_pdgid[iP] == 211) posPion.push_back(iP);
						if (genparticles_pdgid[iP] == -211) negPion.push_back(iP);
					}
					if (genparticles_pdgid[iP] == 35)
					{
						higgsIndex = iP;
						HiggsFound = true;
					}
					if (genparticles_pdgid[iP] == 25)
					{
						higgsSMIndex = iP;
						HiggsSMFound = true;
					}
					if (genparticles_pdgid[iP] == 24) isWplus = true;
					if (genparticles_pdgid[iP] == -24) isWminus = true;
					if (genparticles_pdgid[iP] == 23) isZ = true;
				}
				// if (posMuon.size()==1&&negMuon.size()==1&&posPion.size()==1&&negPion.size()==1) {         	//std::cout << "H->aa->2mu2tau : " << std::endl;
				//std::cout << "Number of mu-   : " << negMuon.size() << std::endl;
				//std::cout << "Number of mu+   : " << posMuon.size() << std::endl;
				//std::cout << "Number of pion- : " << negPion.size() << std::endl;
				//std::cout << "Number of pion+ : " << posPion.size() << std::endl;
				//std::cout << std::endl;
				// }
				if (posMuon.size() == 1 && negMuon.size() == 1 && negPion.size() == 1 && posPion.size() == 1)
				{

					unsigned int indexMu1 = posMuon.at(0);
					unsigned int indexMu2 = negMuon.at(0);

					unsigned int indexPi1 = posPion.at(0);
					unsigned int indexPi2 = negPion.at(0);

					TLorentzVector muon1LV;
					muon1LV.SetXYZT(genparticles_px[indexMu1],
						genparticles_py[indexMu1],
						genparticles_pz[indexMu1],
						genparticles_e[indexMu1]);

					TLorentzVector muon2LV;
					muon2LV.SetXYZT(genparticles_px[indexMu2],
						genparticles_py[indexMu2],
						genparticles_pz[indexMu2],
						genparticles_e[indexMu2]);

					TLorentzVector pion1LV;
					pion1LV.SetXYZT(genparticles_px[indexPi1],
						genparticles_py[indexPi1],
						genparticles_pz[indexPi1],
						genparticles_e[indexPi1]);

					TLorentzVector pion2LV;
					pion2LV.SetXYZT(genparticles_px[indexPi2],
						genparticles_py[indexPi2],
						genparticles_pz[indexPi2],
						genparticles_e[indexPi2]);

					TLorentzVector DiMuonLV = muon1LV + muon2LV;
					TLorentzVector DiTrackLV = pion1LV + pion2LV;
					TLorentzVector Met4LV;
					Met4LV.SetXYZM(metx, mety, 0, 0);

					MuMu_Mass = DiMuonLV.M();
					MuMu_Pt = DiMuonLV.Pt();
					MuMu_DR = deltaR(muon1LV.Eta(), muon1LV.Phi(), muon2LV.Eta(), muon2LV.Phi());
					MuMuMET_DPhi = dPhiFrom2P(DiMuonLV.Px(), DiMuonLV.Py(), Met4LV.Px(), Met4LV.Py());
					TrkTrk_Mass = DiTrackLV.M();
					TrkTrk_Pt = DiTrackLV.Pt();
					TrkTrk_DR = deltaR(pion1LV.Eta(), pion1LV.Phi(), pion2LV.Eta(), pion2LV.Phi());
					TrkTrkMET_DPhi = dPhiFrom2P(DiTrackLV.Px(), DiTrackLV.Py(), Met4LV.Px(), Met4LV.Py());
					MET_Pt = Met4LV.Pt();

					if (MuMu_DR > 1.5 || TrkTrk_DR > 1.5) continue;
					if (muon1LV.Pt() < 24 && muon2LV.Pt() < 24) continue;

					tree_GenLev->Fill();
				}
			}	//end "if" for MC studies

			if (HiggsFound)
			{
				TLorentzVector higgsLV;
				higgsLV.SetXYZT(genparticles_px[higgsIndex],
					genparticles_py[higgsIndex],
					genparticles_pz[higgsIndex],
					genparticles_e[higgsIndex]);
				higgsPt = higgsLV.Pt();
				higgsTree->Fill();

				if (applyHiggsPtWeight && is4tau)
				{
					double HiggsPtForWeighting = higgsPt;
					if (higgsPt > 500) HiggsPtForWeighting = 499;
					double higgsPtWeight = 1;
					if (isVH)
					{
						if (isWplus)
							higgsPtWeight = higgsPt_WPlusH->GetBinContent(higgsPt_WPlusH->FindBin(HiggsPtForWeighting));
						if (isWminus)
							higgsPtWeight = higgsPt_WMinusH->GetBinContent(higgsPt_WMinusH->FindBin(HiggsPtForWeighting));
						if (isZ)
							higgsPtWeight = higgsPt_ZH->GetBinContent(higgsPt_ZH->FindBin(HiggsPtForWeighting));
					}
					else
					{
						higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
					}
					weight *= higgsPtWeight;
					HiggsPtWeightH->Fill(higgsPtWeight);
					//	   std::cout << "HiggsPt weight (Pythia) = " << higgsPtWeight << std::endl;
				}
			}

			if (HiggsSMFound)
			{
				TLorentzVector higgsLV;
				higgsLV.SetXYZT(genparticles_px[higgsSMIndex],
					genparticles_py[higgsSMIndex],
					genparticles_pz[higgsSMIndex],
					genparticles_e[higgsSMIndex]);
				higgsSMPt = higgsLV.Pt();

				if (applyHiggsPtWeight && !is4tau)
				{
					double HiggsPtForWeighting = higgsSMPt;
					if (higgsSMPt > 500) HiggsPtForWeighting = 499;
					double higgsPtWeight = 1;
					higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
					weight *= higgsPtWeight;
					HiggsPtWeightH->Fill(higgsPtWeight);
					//	   std::cout << "HiggsPt weight (Madgraph) = " << higgsPtWeight << std::endl;
				}

				higgsSMTree->Fill();
			}

			counter_InputEventsH->Fill(1.0, weight);

			if (isData)
			{
				if (applyGoodRunSelection)
				{
					bool lumi = false;
					int n = event_run;
					int lum = event_luminosityblock;

					std::string num = std::to_string(n);
					std::string lnum = std::to_string(lum);
					for (const auto &a: periods)
					{
						if (num.c_str() == a.name)
						{
						 				//	      std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
							//std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

							for (auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b)
							{

								//   cout<<b->lower<<"  "<<b->bigger<<endl;
								if (lum >= b->lower && lum <= b->bigger) lumi = true;
							}
							auto last = std::prev(a.ranges.end());
							// std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
							if ((lum >= last->lower && lum <= last->bigger)) lumi = true;
						}
					}
					if (!lumi) continue;
				}
			}

			float puweight = 1;
			if (!isData)
			{
				puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
				//       std::cout << "n(true interactions) = " << numtruepileupinteractions << "   :  PU weight = " << puweight << std::endl;
			}
			puWeightH->Fill(puweight, 1.0);
			weight *= puweight;

			unsigned int nIsoMu = 0;

			bool isIsoMuFilter = false;

			unsigned int nfilters = hltfilters->size();
			for (unsigned int i = 0; i < nfilters; ++i)
			{
			 	//       std::cout << hltfilters->at(i) << std::endl;
				TString HLTFilter(hltfilters->at(i));
				if (HLTFilter == IsoMuFilterName)
				{
					nIsoMu = i;
					isIsoMuFilter = true;
				}
			}

			if (!isIsoMuFilter)
			{
				cout << "Filter " << IsoMuFilterName << " not found " << endl;
				exit(-1);
			}

			// ********************
			// selecting good muons
			// ********************
			vector < unsigned int > muons;
			muons.clear();
			vector<bool> passSingleMu;
			passSingleMu.clear();
			vector<bool> passMu17;
			passMu17.clear();
			vector<bool> passMu8;
			passMu8.clear();
			vector<bool> passMu12;
			passMu12.clear();
			vector<bool> passMu10;
			passMu10.clear();
			vector<bool> passMu5;
			passMu5.clear();
			vector<bool> passMu45;
			passMu45.clear();
			for (UInt_t i = 0; i < muon_count; i++)
			{
				bool muonID = muon_isMedium[i];	// MC

				if (!muonID) continue;
				if (fabs(muon_dxy[i]) > dxyMuonCut) continue;
				if (fabs(muon_dz[i]) > dzMuonCut) continue;
				if (muon_pt[i] < ptGoodMuonCut) continue;
				if (fabs(muon_eta[i]) > etaMuonCut) continue;
				//      cout << "muon pt = " << muon_pt[i] << endl;
				rand.SetSeed((int)((muon_eta[i] + 2.41) *100000));
				double rannum = rand.Rndm();
				// old staff for trigger efficiency studies
				double effMu17 = MuLegEfficiency(muon_pt[i], muon_eta[i], 17.0, 2.4);
				double effMu8 = MuLegEfficiency(muon_pt[i], muon_eta[i], 8.0, 2.4);
				double effMu12 = MuLegEfficiency(muon_pt[i], muon_eta[i], 12.0, 2.4);
				double effMu10 = MuLegEfficiency(muon_pt[i], muon_eta[i], 10.0, 2.4);
				double effMu5 = MuLegEfficiency(muon_pt[i], muon_eta[i], 5.0, 2.4);
				double effMu45 = MuLegEfficiency(muon_pt[i], muon_eta[i], 45.0, 2.1);
				bool isMu17fired = effMu17 > rannum;
				bool isMu8fired = effMu8 > rannum;
				bool isMu12fired = effMu12 > rannum;
				bool isMu10fired = effMu10 > rannum;
				bool isMu5fired = effMu5 > rannum;
				bool isMu45fired = effMu45 > rannum;
				passMu17.push_back(isMu17fired);
				passMu12.push_back(isMu12fired);
				passMu10.push_back(isMu10fired);
				passMu8.push_back(isMu8fired);
				passMu5.push_back(isMu5fired);
				passMu45.push_back(isMu45fired);
				//
				muons.push_back(i);
			}

			nGoodMuonsH->Fill(float(muons.size()), weight);

			if (muons.size() < 2) continue;	// quit event if number of good muons < 3
			counter_MuonSizeGTE2H->Fill(1., weight);

			// **********************
			// selecting good tracks
			// **********************
			vector < unsigned int > tracks;
			tracks.clear();

			for (unsigned int iTrk = 0; iTrk < track_count; ++iTrk)
			{
				if (fabs(track_charge[iTrk]) < 0.1) continue;	// make sure we are not taking neutral stuff
				if (fabs(track_dxy[iTrk]) > dxyTrkLooseCut) continue;
				if (fabs(track_dz[iTrk]) > dzTrkLooseCut) continue;
				if (fabs(track_eta[iTrk]) > etaTrkCut) continue;
				if (fabs(track_pt[iTrk]) < ptTrkLooseCut) continue;
				if (abs(track_ID[iTrk]) != 11 && abs(track_ID[iTrk]) != 13 && abs(track_ID[iTrk]) != 211) continue;

				tracks.push_back(iTrk);
			}	// end loop over good tracks

			// ************************************************
			// **********Some trigger studies ****************
			// ***TripleMu triggers vs. DoubleMuSS trigger ***
			// ************************************************

			// checking 4-muon final states
			if (muons.size() >= numberOfMuons)
			{
				bool isSingleMuon = false;
				for (unsigned int im = 0; im < muons.size(); ++im)
				{
					if (passMu45.at(im))
					{
						isSingleMuon = true;
						break;
					}
				}
				bool isDoubleMuSS = false;
				for (unsigned int im1 = 0; im1 < muons.size() - 1; ++im1)
				{
					for (unsigned int im2 = im1 + 1; im2 < muons.size(); ++im2)
					{
						int mu1 = muons.at(im1);
						int mu2 = muons.at(im2);
						if (muon_charge[mu1] *muon_charge[mu2] > 0.0)
						{
							if ((passMu17.at(im1) && passMu8.at(im2)) || (passMu17.at(im2) && passMu8.at(im1)))
							{
								isDoubleMuSS = true;
								break;
							}
						}
					}
					if (isDoubleMuSS) break;
				}
				// 012
				bool triple012 = (passMu12.at(0) && passMu10.at(1) && passMu5.at(2));
				bool triple021 = (passMu12.at(0) && passMu10.at(2) && passMu5.at(1));
				bool triple102 = (passMu12.at(1) && passMu10.at(0) && passMu5.at(2));
				bool triple120 = (passMu12.at(1) && passMu10.at(2) && passMu5.at(0));
				bool triple210 = (passMu12.at(2) && passMu10.at(1) && passMu5.at(0));
				bool triple201 = (passMu12.at(2) && passMu10.at(0) && passMu5.at(1));
				bool Triple012 = triple012 || triple021 || triple102 || triple120 || triple210 || triple201;
				// 013
				bool Triple013 = false;
				if (muons.size() > 3)
				{
					bool triple013 = (passMu12.at(0) && passMu10.at(1) && passMu5.at(3));
					bool triple031 = (passMu12.at(0) && passMu10.at(3) && passMu5.at(1));
					bool triple103 = (passMu12.at(1) && passMu10.at(0) && passMu5.at(3));
					bool triple130 = (passMu12.at(1) && passMu10.at(3) && passMu5.at(0));
					bool triple310 = (passMu12.at(3) && passMu10.at(1) && passMu5.at(0));
					bool triple301 = (passMu12.at(3) && passMu10.at(0) && passMu5.at(1));
					Triple013 = triple013 || triple031 || triple103 || triple130 || triple310 || triple301;
				}
				// 123
				bool Triple123 = false;
				if (muons.size() > 3)
				{
					bool triple123 = (passMu12.at(1) && passMu10.at(2) && passMu5.at(3));
					bool triple132 = (passMu12.at(1) && passMu10.at(3) && passMu5.at(2));
					bool triple213 = (passMu12.at(2) && passMu10.at(1) && passMu5.at(3));
					bool triple231 = (passMu12.at(2) && passMu10.at(3) && passMu5.at(1));
					bool triple312 = (passMu12.at(3) && passMu10.at(1) && passMu5.at(2));
					bool triple321 = (passMu12.at(3) && passMu10.at(2) && passMu5.at(1));
					Triple123 = triple123 || triple132 || triple213 || triple231 || triple312 || triple321;
				}
				// TripleMu
				bool isTripleMuon = Triple012 || Triple013 || Triple123;
				if (isSingleMuon) histWeightsSingleMuH->Fill(1.0, weight);
				if (isDoubleMuSS) histWeightsDoubleMuSSH->Fill(1.0, weight);
				if (isTripleMuon) histWeightsTripleMuH->Fill(1.0, weight);
				bool isAllTriggers = isTripleMuon || isDoubleMuSS || isSingleMuon;
				if (isAllTriggers) histWeightsAllTriggersH->Fill(1.0, weight);
			}

			// *************************b Jet Rejection*************************//
			std::vector < unsigned int > EveBJetMult;
			EveBJetMult.clear();	// all B Jets

			for (unsigned int ipfjet = 0; ipfjet < pfjet_count; ipfjet++)
			{
				TLorentzVector Jet4Momentum;
				Jet4Momentum.SetPtEtaPhiE(pfjet_pt[ipfjet], pfjet_eta[ipfjet], pfjet_phi[ipfjet], pfjet_e[ipfjet]);

				// b Jet counting
				if (Jet4Momentum.Pt() > 20 && Jet4Momentum.Pt() < 1000 && fabs(Jet4Momentum.Eta()) < 2.4 &&
					pfjet_btag[ipfjet][0] > btagCut) EveBJetMult.push_back(ipfjet);
			}

			BjetsMultipH->Fill(EveBJetMult.size(), weight);

			if (!isData)
			{
				double bTagWeight_0b = 1.;
				//// Here B-Tagging Scale Factors for 0 b-tagged jets	////////
				for (unsigned int ibjet = 0; ibjet < EveBJetMult.size(); ibjet++)
				{
					int indexb = EveBJetMult.at(ibjet);

					TLorentzVector Jet4Momentum;
					Jet4Momentum.SetPtEtaPhiE(pfjet_pt[indexb], pfjet_eta[indexb], pfjet_phi[indexb], pfjet_e[indexb]);

					float jetPt = Jet4Momentum.Pt();
					float jetAbsEta = fabs(Jet4Momentum.Eta());
					int flavor = abs(pfjet_flavour[indexb]);

					double SF_CSV;

					if (flavor == 5) SF_CSV = reader_B.eval_auto_bounds("central", BTagEntry::FLAV_B, jetAbsEta, jetPt);	//B

					else if (flavor == 4) SF_CSV = reader_C.eval_auto_bounds("central", BTagEntry::FLAV_C, jetAbsEta, jetPt);	//C

					else SF_CSV = reader_L.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, jetAbsEta, jetPt);	//UDSG

					bTagWeight_0b = bTagWeight_0b *(1. - SF_CSV);
				}

				histWeightsBTagH->Fill(1., bTagWeight_0b);
				BTagEventWeightH->Fill(bTagWeight_0b, 1.);
			}

			if (EveBJetMult.size() > 0) continue;

			if (!isData) histWeightsZeroBTagH->Fill(1., 1.);

			///////////////////////////////////////////////////////
			////////******A->mumu Candidates*********////////////
			///////////////////////////////////////////////////////

			float maxPtSumMuMu = ptSumDiMuCut;
			int iAmumu1candidate = -1;
			int iAmumu2candidate = -1;

			for (unsigned int i1 = 0; i1 < muons.size() - 1; ++i1)
			{
				int index1 = muons.at(i1);
				for (unsigned int i2 = i1 + 1; i2 < muons.size(); ++i2)
				{
					int index2 = muons.at(i2);

					float ptSum = TMath::Sqrt((muon_px[index1] + muon_px[index2]) *(muon_px[index1] + muon_px[index2]) + (muon_py[index1] + muon_py[index2]) *(muon_py[index1] + muon_py[index2]));
					float charge = muon_charge[index1] *muon_charge[index2];
					bool chargeSelection = charge < 0;
					float dRmuons = deltaR(muon_eta[index1], muon_phi[index1],
						muon_eta[index2], muon_phi[index2]);

					if (chargeSelection && dRmuons < dRMuonsCut && ptSum > maxPtSumMuMu && fabs(muon_dxy[index1]) < dxyDiMuCut &&
						fabs(muon_dxy[index2]) < dxyDiMuCut && fabs(muon_dz[index1]) < dzDiMuCut && fabs(muon_dz[index2]) < dzDiMuCut)
					{
						iAmumu1candidate = index1;
						iAmumu2candidate = index2;
						maxPtSumMuMu = ptSum;
					}
				}
			}

			if (iAmumu1candidate < 0) continue;
			if (iAmumu2candidate < 0) continue;

			////////////////////////////////////////////////////
			//////************Here SM Trigger ************//////
			////////////////////////////////////////////////////

			bool mu1MatchIsoMu = false;

			for (unsigned int iT = 0; iT < trigobject_count; ++iT)
			{
				float dRtrig = deltaR(muon_eta[iAmumu1candidate], muon_phi[iAmumu1candidate],
					trigobject_eta[iT], trigobject_phi[iT]);
				if (dRtrig > DRTrigMatch) continue;
				if (trigobject_filters[iT][nIsoMu])
					mu1MatchIsoMu = true;
			}

			bool mu2MatchIsoMu = false;

			for (unsigned int iT = 0; iT < trigobject_count; ++iT)
			{
				float dRtrig = deltaR(muon_eta[iAmumu2candidate], muon_phi[iAmumu2candidate],
					trigobject_eta[iT], trigobject_phi[iT]);
				if (dRtrig > DRTrigMatch) continue;
				if (trigobject_filters[iT][nIsoMu])
					mu2MatchIsoMu = true;
			}

			bool muLPtCut = muon_pt[iAmumu1candidate] > ptIsoMuonCut && fabs(muon_eta[iAmumu1candidate]) < etaMuonCut;
			bool muTPtCut = muon_pt[iAmumu2candidate] > ptIsoMuonCut && fabs(muon_eta[iAmumu2candidate]) < etaMuonCut;

			// trigger condition
			bool isTriggerMatched = true;
			if (applyTriggerMatch)
			{
				isTriggerMatched = (mu1MatchIsoMu && muLPtCut) || (mu2MatchIsoMu && muTPtCut);
			}
			else
			{
				isTriggerMatched = muLPtCut;
			}
			if (!isTriggerMatched) continue;

			double triggerWeight = 1;
			double idIsoWeightL = 1;
			double idIsoWeightT = 1;

			// *******************************************************
			//   Trigger efficiency, MuonID, and tracking efficiencies
			// *******************************************************

			if (!isData)
			{

				double effDataMu1 = SF_muon->get_EfficiencyData(muon_pt[iAmumu1candidate], muon_eta[iAmumu1candidate]);

				double effMCMu1 = SF_muon->get_EfficiencyMC(muon_pt[iAmumu1candidate], muon_eta[iAmumu1candidate]);

				double effDataMu2 = SF_muon->get_EfficiencyData(muon_pt[iAmumu2candidate], muon_eta[iAmumu2candidate]);

				double effMCMu2 = SF_muon->get_EfficiencyMC(muon_pt[iAmumu2candidate], muon_eta[iAmumu2candidate]);

				double trigWeightData = effDataMu1 + effDataMu2 - effDataMu1 * effDataMu2;

				double trigWeightMC = effMCMu1 + effMCMu2 - effMCMu1 * effMCMu2;

				if (applyTriggerMatch)
				{
					if (trigWeightMC > 0) triggerWeight = trigWeightData / trigWeightMC;
				}
				else triggerWeight = trigWeightData;

				correctionWS->var("m_pt")->setVal(muon_pt[iAmumu1candidate]);
				correctionWS->var("m_eta")->setVal(muon_eta[iAmumu1candidate]);
				correctionWS->var("m_iso")->setVal(0.0);
				idIsoWeightL = correctionWS->function("m_idiso_ic_ratio")->getVal();
				correctionWS->var("m_pt")->setVal(muon_pt[iAmumu2candidate]);
				correctionWS->var("m_eta")->setVal(muon_eta[iAmumu2candidate]);
				correctionWS->var("m_iso")->setVal(0.0);
				idIsoWeightT = correctionWS->function("m_idiso_ic_ratio")->getVal();
			}

			triggerWeightH->Fill(triggerWeight, 1.0);
			weight *= triggerWeight;
			idIsoWeightH->Fill(idIsoWeightL * idIsoWeightT, 1.0);
			weight = weight *idIsoWeightL * idIsoWeightT;

			TLorentzVector Amumu1candidate4;
			Amumu1candidate4.SetXYZM(muon_px[iAmumu1candidate],
						 muon_py[iAmumu1candidate],
						 muon_pz[iAmumu1candidate],
						 MuMass);

			TLorentzVector Amumu2candidate4;
			Amumu2candidate4.SetXYZM(muon_px[iAmumu2candidate],
						 muon_py[iAmumu2candidate],
						 muon_pz[iAmumu2candidate],
						 MuMass);

			double scale1 = MuonMomentumScale(Amumu1candidate4.Eta());
			double scale2 = MuonMomentumScale(Amumu2candidate4.Eta());

			//			std::cout << "scale1 = " << scale1 
			//				  << "   scale2 = " << scale2 << std::endl;

			TLorentzVector Amumu1candidate4_muScaleUp;
			Amumu1candidate4_muScaleUp.SetXYZM((1.0+scale1)*muon_px[iAmumu1candidate],
							   (1.0+scale1)*muon_py[iAmumu1candidate],
							   (1.0+scale1)*muon_pz[iAmumu1candidate],
							   MuMass);
			TLorentzVector Amumu1candidate4_muScaleDown;
			Amumu1candidate4_muScaleDown.SetXYZM((1.0-scale1)*muon_px[iAmumu1candidate],
							     (1.0-scale1)*muon_py[iAmumu1candidate],
							     (1.0-scale1)*muon_pz[iAmumu1candidate],
							   MuMass);

			TLorentzVector Amumu2candidate4_muScaleUp;
			Amumu2candidate4_muScaleUp.SetXYZM((1.0+scale2)*muon_px[iAmumu2candidate],
							   (1.0+scale2)*muon_py[iAmumu2candidate],
							   (1.0+scale2)*muon_pz[iAmumu2candidate],
							   MuMass);
			TLorentzVector Amumu2candidate4_muScaleDown;
			Amumu2candidate4_muScaleDown.SetXYZM((1.0-scale2)*muon_px[iAmumu2candidate],
							     (1.0-scale2)*muon_py[iAmumu2candidate],
							     (1.0-scale2)*muon_pz[iAmumu2candidate],
							     MuMass);
			/*
			std::cout << "pt1 = " << Amumu1candidate4.Pt() 
				  << "   up = " << Amumu1candidate4_muScaleUp.Pt()
				  << "   down"  << Amumu1candidate4_muScaleDown.Pt()
				  << std::endl;
			
			std::cout << "pt2 = " << Amumu2candidate4.Pt() 
				  << "   up = " << Amumu2candidate4_muScaleUp.Pt()
				  << "   down"  << Amumu2candidate4_muScaleDown.Pt()
				  << std::endl;
			*/

			TLorentzVector diMuon4 = Amumu1candidate4 + Amumu2candidate4;
			TLorentzVector diMuon4_muScaleUp = Amumu1candidate4_muScaleUp + Amumu2candidate4_muScaleUp;
			TLorentzVector diMuon4_muScaleDown = Amumu1candidate4_muScaleDown + Amumu2candidate4_muScaleDown;

			TLorentzVector Met4;
			Met4.SetXYZM(metx, mety, 0, 0);
			double metx_muScaleUp = metx + Amumu1candidate4.Px() + Amumu2candidate4.Px() - Amumu1candidate4_muScaleUp.Px() - Amumu2candidate4_muScaleUp.Px();
			double mety_muScaleUp = metx + Amumu1candidate4.Py() + Amumu2candidate4.Py() - Amumu1candidate4_muScaleUp.Py() - Amumu2candidate4_muScaleUp.Py();
			double metx_muScaleDown = metx + Amumu1candidate4.Px() + Amumu2candidate4.Px() - Amumu1candidate4_muScaleDown.Px() - Amumu2candidate4_muScaleDown.Px();
			double mety_muScaleDown = metx + Amumu1candidate4.Py() + Amumu2candidate4.Py() - Amumu1candidate4_muScaleDown.Py() - Amumu2candidate4_muScaleDown.Py();

			TLorentzVector Met4_muScaleUp;
			Met4_muScaleUp.SetXYZM(metx,mety,0.,0.);
			TLorentzVector Met4_muScaleDown;
			Met4_muScaleDown.SetXYZM(metx,mety,0.,0.);

			TLorentzVector Met4_JetEnUp;
			Met4_JetEnUp.SetXYZM(metx_JetEnUp,mety_JetEnUp,0.,0.);
			TLorentzVector Met4_JetEnDown;
			Met4_JetEnDown.SetXYZM(metx_JetEnDown,mety_JetEnDown,0.,0.);

			TLorentzVector Met4_UnclEnUp;
			Met4_UnclEnUp.SetXYZM(metx_UnclEnUp,mety_UnclEnUp,0.,0.);
			TLorentzVector Met4_UnclEnDown;
			Met4_UnclEnDown.SetXYZM(metx_UnclEnDown,mety_UnclEnDown,0.,0.);

			float dimuonMass = diMuon4.M();
			float dimuonMass_muScaleUp = diMuon4_muScaleUp.M();
			float dimuonMass_muScaleDown = diMuon4_muScaleDown.M();	
			float dimuonPt = diMuon4.Pt();
			float dimuonPt_muScaleUp = diMuon4_muScaleUp.Pt();
			float dimuonPt_muScaleDown = diMuon4_muScaleDown.Pt();
			float drMuMu = deltaR(Amumu1candidate4.Eta(), Amumu1candidate4.Phi(), Amumu2candidate4.Eta(), Amumu2candidate4.Phi());
			float metPt = Met4.Pt();
			float dPhiMuMuMet = dPhiFrom2P(diMuon4.Px(), diMuon4.Py(), Met4.Px(), Met4.Py());

			if (dimuonMass < massDiMuCut || dimuonMass > 22.) continue;	// requiring InvMass ("a" boson mass) to be greater than 2*mass_tau and less than 22

			// filling histograms (muon kinematics)
			dimuonMassH->Fill(dimuonMass, weight);
			ptAmumu1candidateH->Fill(Amumu1candidate4.Pt(), weight);
			ptAmumu2candidateH->Fill(Amumu2candidate4.Pt(), weight);
			etaAmumu1candidateH->Fill(Amumu1candidate4.Eta(), weight);
			etaAmumu2candidateH->Fill(Amumu2candidate4.Eta(), weight);

			// filling histograms (muon impact parameter)
			dxyAmumu1candidateH->Fill(muon_dxy[iAmumu1candidate], weight);
			dxyAmumu2candidateH->Fill(muon_dxy[iAmumu2candidate], weight);
			dzAmumu1candidateH->Fill(muon_dz[iAmumu1candidate], weight);
			dzAmumu2candidateH->Fill(muon_dz[iAmumu2candidate], weight);

			// filling histograms (MET and counter)
			MetH->Fill(Met4.Pt(), weight);
			counter_MuonKinematicsH->Fill(1.0, weight);

			///////////////////////////////////////////////////////
			////////******A->tautau Candidates*********////////////
			///////////////////////////////////////////////////////

			float maxPtSumTrkTrk = ptSumDiTrkCut;
			int iAtautau1candidate = -1;
			int iAtautau2candidate = -1;

			for (unsigned int i1 = 0; i1 < tracks.size() - 1; ++i1)
			{
				int index1 = tracks.at(i1);

				TLorentzVector trk1st4;
				trk1st4.SetXYZM(track_px[index1],
					track_py[index1],
					track_pz[index1],
					track_mass[index1]);

				TLorentzVector Amumu1candidateDiff1 = Amumu1candidate4 - trk1st4;
				TLorentzVector Amumu2candidateDiff1 = Amumu2candidate4 - trk1st4;

				if (Amumu1candidateDiff1.P() < 0.1 || Amumu2candidateDiff1.P() < 0.1) continue;

				for (unsigned int i2 = i1 + 1; i2 < tracks.size(); ++i2)
				{
					int index2 = tracks.at(i2);

					TLorentzVector trk2nd4;
					trk2nd4.SetXYZM(track_px[index2],
						track_py[index2],
						track_pz[index2],
						track_mass[index2]);

					TLorentzVector Amumu1candidateDiff2 = Amumu1candidate4 - trk2nd4;
					TLorentzVector Amumu2candidateDiff2 = Amumu2candidate4 - trk2nd4;

					TLorentzVector diTrack4 = trk1st4 + trk2nd4;
					TLorentzVector Visible4 = diTrack4 + diMuon4;

					float dRAbosonCandidates = deltaR(diMuon4.Eta(), diMuon4.Phi(),
						diTrack4.Eta(), diTrack4.Phi());

					if (Amumu1candidateDiff2.P() < 0.1 || Amumu2candidateDiff2.P() < 0.1 || dRAbosonCandidates < dRMuonsCut) continue;

					float ptSum = diTrack4.Pt();
					float charge = track_charge[index1] *track_charge[index2];
					bool chargeSelection = charge < 0;
					float drTrkTrk = deltaR(trk1st4.Eta(), trk1st4.Phi(),
						trk2nd4.Eta(), trk2nd4.Phi());

					if (drTrkTrk < dRMuonsCut && chargeSelection && track_pt[index1] > ptTrkCut && track_pt[index2] > ptTrkCut && ptSum > maxPtSumTrkTrk &&
						fabs(track_dxy[index1]) < dxyTrkCut && fabs(track_dxy[index2]) < dxyTrkCut && fabs(track_dz[index1]) < dzTrkCut && fabs(track_dz[index2]) < dzTrkCut)
					{
						iAtautau1candidate = index1;
						iAtautau2candidate = index2;
						maxPtSumTrkTrk = ptSum;
					}
				}
			}

			if (iAtautau1candidate < 0) continue;
			if (iAtautau2candidate < 0) continue;

			////Identifying tracks///
			int IDAtautau1candidate = -1;
			int IDAtautau2candidate = -1;

			if (abs(track_ID[iAtautau1candidate]) == 13) IDAtautau1candidate = 0;	//// is muon
			if (abs(track_ID[iAtautau1candidate]) == 11) IDAtautau1candidate = 1;	//// is electron
			if (abs(track_ID[iAtautau1candidate]) == 211) IDAtautau1candidate = 2;	//// is hadron

			if (abs(track_ID[iAtautau2candidate]) == 13) IDAtautau2candidate = 0;	//// is muon
			if (abs(track_ID[iAtautau2candidate]) == 11) IDAtautau2candidate = 1;	//// is electron
			if (abs(track_ID[iAtautau2candidate]) == 211) IDAtautau2candidate = 2;	//// is hadron

			TLorentzVector Atautau1candidate4;
			Atautau1candidate4.SetXYZM(track_px[iAtautau1candidate],
				track_py[iAtautau1candidate],
				track_pz[iAtautau1candidate],
				track_mass[iAtautau1candidate]);

			TLorentzVector Atautau2candidate4;
			Atautau2candidate4.SetXYZM(track_px[iAtautau2candidate],
				track_py[iAtautau2candidate],
				track_pz[iAtautau2candidate],
				track_mass[iAtautau2candidate]);

			TLorentzVector diTrack4 = Atautau1candidate4 + Atautau2candidate4;

			TLorentzVector Visible4 = diTrack4 + diMuon4;
			TLorentzVector Visible4_muScaleUp = diTrack4 + diMuon4_muScaleUp;
			TLorentzVector Visible4_muScaleDown = diTrack4 + diMuon4_muScaleDown;

			TLorentzVector VisibleandMet4 = Visible4 + Met4;
			TLorentzVector VisibleandMet4_muScaleUp = Visible4_muScaleUp + Met4_muScaleUp;
			TLorentzVector VisibleandMet4_muScaleDown = Visible4_muScaleDown + Met4_muScaleDown;
			TLorentzVector VisibleandMet4_JetEnUp = Visible4 + Met4_JetEnUp;
			TLorentzVector VisibleandMet4_JetEnDown = Visible4 + Met4_JetEnDown;
			TLorentzVector VisibleandMet4_UnclEnUp = Visible4 + Met4_UnclEnUp;
			TLorentzVector VisibleandMet4_UnclEnDown = Visible4 + Met4_UnclEnDown;

			float ditrackMass = diTrack4.M();
			float ditrackPt = diTrack4.Pt();
			float drTrkTrk = deltaR(Atautau1candidate4.Eta(), Atautau1candidate4.Phi(), Atautau2candidate4.Eta(), Atautau2candidate4.Phi());
			float dPhiTrkTrkMet = dPhiFrom2P(diTrack4.Px(), diTrack4.Py(), Met4.Px(), Met4.Py());
			float visibleMass = Visible4.M();
			float visiblePt = Visible4.Pt();
			float visibleandMetMass = VisibleandMet4.M();
			float visibleandMetMass_muScaleUp = VisibleandMet4_muScaleUp.M();
			float visibleandMetMass_muScaleDown = VisibleandMet4_muScaleDown.M();
			float visibleandMetMass_JetEnUp = VisibleandMet4_JetEnUp.M();
			float visibleandMetMass_JetEnDown = VisibleandMet4_JetEnDown.M();
			float visibleandMetMass_UnclEnUp = VisibleandMet4_UnclEnUp.M();
			float visibleandMetMass_UnclEnDown = VisibleandMet4_UnclEnDown.M();

			if (ditrackMass > 22.) continue;	// requiring InvMass ("a" boson mass) to be less than 22;

			// filling histograms (tracks kinematics)
			ditrackMassH->Fill(diTrack4.M(), weight);
			ptAtautau1candidateH->Fill(Atautau1candidate4.Pt(), weight);
			ptAtautau2candidateH->Fill(Atautau2candidate4.Pt(), weight);
			etaAtautau1candidateH->Fill(Atautau1candidate4.Eta(), weight);
			etaAtautau2candidateH->Fill(Atautau2candidate4.Eta(), weight);

			// filling histograms (tracks impact parameter)
			dxyAtautau1candidateH->Fill(track_dxy[iAtautau1candidate], weight);
			dxyAtautau2candidateH->Fill(track_dxy[iAtautau2candidate], weight);
			dzAtautau1candidateH->Fill(track_dz[iAtautau1candidate], weight);
			dzAtautau2candidateH->Fill(track_dz[iAtautau2candidate], weight);

			//////////////////////////////////////// ***********************************////////////////////////////////////////////
			////////////////////////////////////////  Beginning of analysis by Category 	////////////////////////////////////////////
			//////////////////////////////////////// ***********************************////////////////////////////////////////////

			const double restartweight = weight;

			for (int iInt = 0; iInt < nDRInterval; ++iInt)	/////starting loop over all Intervals (Categories) for Dimuon Pair
			{
			 	//============================================================================================================//
				weight = restartweight;

				/////////////////////////////////////////////////////////////////////
				//   counting tracks around Amumucandidates &Atautaucandidates  	///
				/////////////////////////////////////////////////////////////////////

				std::vector < unsigned int > trkAmumu1candidate;
				trkAmumu1candidate.clear();	// all tracks
				std::vector < unsigned int > trkAmumu2candidate;
				trkAmumu2candidate.clear();	// all tracks
				std::vector < unsigned int > trkAtautau1candidate;
				trkAtautau1candidate.clear();	// all tracks
				std::vector < unsigned int > trkAtautau2candidate;
				trkAtautau2candidate.clear();	// all tracks

				std::vector < unsigned int > soft_trkAmumu1candidate;
				soft_trkAmumu1candidate.clear();	// all soft tracks
				std::vector < unsigned int > soft_trkAmumu2candidate;
				soft_trkAmumu2candidate.clear();	// all soft tracks
				std::vector < unsigned int > soft_trkAtautau1candidate;
				soft_trkAtautau1candidate.clear();	// all soft tracks
				std::vector < unsigned int > soft_trkAtautau2candidate;
				soft_trkAtautau2candidate.clear();	// all soft tracks

				std::vector < unsigned int > sig_trkAmumu1candidate;
				sig_trkAmumu1candidate.clear();	// all signal tracks
				std::vector < unsigned int > sig_trkAmumu2candidate;
				sig_trkAmumu2candidate.clear();	// all signal tracks
				std::vector < unsigned int > sig_trkAtautau1candidate;
				sig_trkAtautau1candidate.clear();	// all signal tracks
				std::vector < unsigned int > sig_trkAtautau2candidate;
				sig_trkAtautau2candidate.clear();	// all signal tracks

				for (unsigned int iTrk = 0; iTrk < track_count; ++iTrk)
				{
					if (fabs(track_charge[iTrk]) < 0.1) continue;	// make sure we are not taking neutral stuff
					if (fabs(track_dxy[iTrk]) > dxyTrkLooseCut) continue;
					if (fabs(track_dz[iTrk]) > dzTrkLooseCut) continue;
					if (fabs(track_eta[iTrk]) > etaTrkCut) continue;
					if (fabs(track_pt[iTrk]) < ptTrkLooseCut) continue;

					TLorentzVector trk4;
					trk4.SetXYZM(track_px[iTrk],
						track_py[iTrk],
						track_pz[iTrk],
						track_mass[iTrk]);

					TLorentzVector Atautau1candidateDiff = Atautau1candidate4 - trk4;
					TLorentzVector Atautau2candidateDiff = Atautau2candidate4 - trk4;
					TLorentzVector Amumu1candidateDiff = Amumu1candidate4 - trk4;
					TLorentzVector Amumu2candidateDiff = Amumu2candidate4 - trk4;

					if (Atautau1candidateDiff.P() > 0.1 && Atautau2candidateDiff.P() > 0.1 && Amumu1candidateDiff.P() > 0.1 && Amumu2candidateDiff.P() > 0.1)	// track is not any of the selected objects
					{

						float drTrkAmumu1 = deltaR(Amumu1candidate4.Eta(), Amumu1candidate4.Phi(),
							track_eta[iTrk], track_phi[iTrk]);
						float drTrkAmumu2 = deltaR(Amumu2candidate4.Eta(), Amumu2candidate4.Phi(),
							track_eta[iTrk], track_phi[iTrk]);
						float drTrkATauTau1 = deltaR(Atautau1candidate4.Eta(), Atautau1candidate4.Phi(),
							track_eta[iTrk], track_phi[iTrk]);
						float drTrkATauTau2 = deltaR(Atautau2candidate4.Eta(), Atautau2candidate4.Phi(),
							track_eta[iTrk], track_phi[iTrk]);

						if (drTrkAmumu1 < DRInterval[iInt])
						{
							trkAmumu1candidate.push_back(iTrk);
							if (fabs(track_pt[iTrk]) < ptTrkSoftCut) soft_trkAmumu1candidate.push_back(iTrk);
							if (track_pt[iTrk] > ptTrkCut && fabs(track_dxy[iTrk]) < dxyTrkCut && fabs(track_dz[iTrk]) < dzTrkCut) sig_trkAmumu1candidate.push_back(iTrk);
						}
						if (drTrkAmumu2 < DRInterval[iInt])
						{
							trkAmumu2candidate.push_back(iTrk);
							if (fabs(track_pt[iTrk]) < ptTrkSoftCut) soft_trkAmumu2candidate.push_back(iTrk);
							if (track_pt[iTrk] > ptTrkCut && fabs(track_dxy[iTrk]) < dxyTrkCut && fabs(track_dz[iTrk]) < dzTrkCut) sig_trkAmumu2candidate.push_back(iTrk);
						}
						if (drTrkATauTau1 < DRInterval[iInt])
						{
							trkAtautau1candidate.push_back(iTrk);
							if (fabs(track_pt[iTrk]) < ptTrkSoftCut) soft_trkAtautau1candidate.push_back(iTrk);
							if (track_pt[iTrk] > ptTrkCut && fabs(track_dxy[iTrk]) < dxyTrkCut && fabs(track_dz[iTrk]) < dzTrkCut) sig_trkAtautau1candidate.push_back(iTrk);
						}
						if (drTrkATauTau2 < DRInterval[iInt])
						{
							trkAtautau2candidate.push_back(iTrk);
							if (fabs(track_pt[iTrk]) < ptTrkSoftCut) soft_trkAtautau2candidate.push_back(iTrk);
							if (track_pt[iTrk] > ptTrkCut && fabs(track_dxy[iTrk]) < dxyTrkCut && fabs(track_dz[iTrk]) < dzTrkCut) sig_trkAtautau2candidate.push_back(iTrk);
						}
					}
				}	// end loop over tracks

				nTracksAmumu1candidateH[iInt]->Fill(float(trkAmumu1candidate.size()), weight);
				nTracksAmumu2candidateH[iInt]->Fill(float(trkAmumu2candidate.size()), weight);
				nTracksAtautau1candidateH[iInt]->Fill(float(trkAtautau1candidate.size()), weight);
				nTracksAtautau2candidateH[iInt]->Fill(float(trkAtautau2candidate.size()), weight);

				nSoftTracksAmumu1candidateH[iInt]->Fill(float(soft_trkAmumu1candidate.size()), weight);
				nSoftTracksAmumu2candidateH[iInt]->Fill(float(soft_trkAmumu2candidate.size()), weight);
				nSoftTracksAtautau1candidateH[iInt]->Fill(float(soft_trkAtautau1candidate.size()), weight);
				nSoftTracksAtautau2candidateH[iInt]->Fill(float(soft_trkAtautau2candidate.size()), weight);

				nSignalTracksAmumu1candidateH[iInt]->Fill(float(sig_trkAmumu1candidate.size()), weight);
				nSignalTracksAmumu2candidateH[iInt]->Fill(float(sig_trkAmumu2candidate.size()), weight);
				nSignalTracksAtautau1candidateH[iInt]->Fill(float(sig_trkAtautau1candidate.size()), weight);
				nSignalTracksAtautau2candidateH[iInt]->Fill(float(sig_trkAtautau2candidate.size()), weight);

				/////////////////////////////////////////////////////////////////////
				//     reconstructing Higgs 4-Vector from kinematic fitting      	///
				/////////////////////////////////////////////////////////////////////

				double alfa = 1;
				double beta = 1;

				double estimator = 1000000;
				double mass_mumu = diMuon4.M();

				TVector3 P_tau1;
				P_tau1.SetXYZ(Atautau1candidate4.Px(), Atautau1candidate4.Py(), Atautau1candidate4.Pz());
				TVector3 P_tau2;
				P_tau2.SetXYZ(Atautau2candidate4.Px(), Atautau2candidate4.Py(), Atautau2candidate4.Pz());

				TVector3 PT_tau1;
				PT_tau1.SetXYZ(Atautau1candidate4.Px(), Atautau1candidate4.Py(), 0.);
				TVector3 PT_tau2;
				PT_tau2.SetXYZ(Atautau2candidate4.Px(), Atautau2candidate4.Py(), 0.);
				TVector3 PT_MET;
				PT_MET.SetXYZ(Met4.Px(), Met4.Py(), 0.);

				for (unsigned int ialfa = 0; ialfa < 100; ++ialfa)
				{

					for (unsigned int ibeta = 0; ibeta < 100; ++ibeta)
					{

						double Estimate = fabs(Function_1(1 + ialfa *0.2, 1 + ibeta *0.2, P_tau1, P_tau2, mass_mumu)) +
							fabs(Function_2(1 + ialfa *0.2, 1 + ibeta *0.2, PT_tau1, PT_tau2, PT_MET));

						if (Estimate <= estimator)
						{
							alfa = 1 + ialfa *0.2;
							beta = 1 + ibeta *0.2;
							estimator = Estimate;
						}
					}
				}

				TLorentzVector FirstTau4;
				FirstTau4.SetXYZM(alfa *Atautau1candidate4.Px(),
					alfa *Atautau1candidate4.Py(),
					alfa *Atautau1candidate4.Pz(),
					1.77682);

				TLorentzVector SecondTau4;
				SecondTau4.SetXYZM(beta *Atautau2candidate4.Px(),
					beta *Atautau2candidate4.Py(),
					beta *Atautau2candidate4.Pz(),
					1.77682);

				TLorentzVector diTau4 = FirstTau4 + SecondTau4;

				TLorentzVector Total4 = diTau4 + diMuon4;

				float totalMass = Total4.M();
				float totalPt = Total4.Pt();

				// Visible Mass Cut
				if (Visible4.M() > visibleMassCut) continue;

				// Diagonal Cut Cut
				if (diTrack4.M() > diMuon4.M()) continue;

				// Total Mass Cut
				if (fabs(Total4.M() - 125.) > TotalMassWindow) continue;

				//////////////////////////////////////// **************************************////////////////////////////////////////////
				////////////////////////////////////////  Defining Sidebands and Signal Region 	////////////////////////////////////////////
				//////////////////////////////////////// **************************************////////////////////////////////////////////

				// Signal Region
				bool signalRegion = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 && trkAtautau1candidate.size() == 0 && trkAtautau2candidate.size() == 0;

				// Control Regions for Background Modeling 

				bool controlNNNN = !signalRegion;

				bool control00NN = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					!signalRegion;

				bool control00SemiIso = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					(soft_trkAtautau1candidate.size() == trkAtautau1candidate.size()) &&
					(soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) &&
					!signalRegion;

				bool controlSoftIso = (trkAmumu1candidate.size() <= 2 && soft_trkAmumu1candidate.size() == trkAmumu1candidate.size()) &&
					(trkAmumu2candidate.size() <= 2 && soft_trkAmumu2candidate.size() == trkAmumu2candidate.size()) &&
					(trkAtautau1candidate.size() <= 2 && soft_trkAtautau1candidate.size() == trkAtautau1candidate.size()) &&
					(trkAtautau2candidate.size() <= 2 && soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) &&
					!signalRegion;

				bool controlVerySoftIso = (trkAmumu1candidate.size() <= 1 && soft_trkAmumu1candidate.size() == trkAmumu1candidate.size()) &&
					(trkAmumu2candidate.size() <= 1 && soft_trkAmumu2candidate.size() == trkAmumu2candidate.size()) &&
					(trkAtautau1candidate.size() <= 1 && soft_trkAtautau1candidate.size() == trkAtautau1candidate.size()) &&
					(trkAtautau2candidate.size() <= 1 && soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) &&
					!signalRegion;

				bool control00SoftIso = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					(trkAtautau1candidate.size() <= 2 && soft_trkAtautau1candidate.size() == trkAtautau1candidate.size()) &&
					(trkAtautau2candidate.size() <= 2 && soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) &&
					!signalRegion;

				bool control00VerySoftIso = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					(trkAtautau1candidate.size() <= 1 && soft_trkAtautau1candidate.size() == trkAtautau1candidate.size()) &&
					(trkAtautau2candidate.size() <= 1 && soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) &&
					!signalRegion;

				bool controlNN00 = trkAtautau1candidate.size() == 0 && trkAtautau2candidate.size() == 0 &&
					!signalRegion;

				bool control00_0SemiIsoorSemiIso0 = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					((trkAtautau1candidate.size() == 0 &&
							soft_trkAtautau2candidate.size() == trkAtautau2candidate.size()) ||
						(trkAtautau2candidate.size() == 0 &&
							soft_trkAtautau1candidate.size() == trkAtautau1candidate.size())
				) &&
					!signalRegion;

				bool control00_0NorN0 = trkAmumu1candidate.size() == 0 && trkAmumu2candidate.size() == 0 &&
					(trkAtautau1candidate.size() == 0 || trkAtautau2candidate.size() == 0) &&
					!signalRegion;

				/////////////////////////****************///////////////////////////////
				//////********filling histograms by Regions &Categories ********//////
				/////////////////////////***************///////////////////////////////

				////****Setting the values of the variables *****/////
				MuMu_Mass = dimuonMass;
				MuMu_Mass_muScaleUp = dimuonMass_muScaleUp;
				MuMu_Mass_muScaleDown = dimuonMass_muScaleDown;
				//				std::cout << "MuMu_Mass = " << MuMu_Mass
				//					  << "    up = " << MuMu_Mass_muScaleUp
				//					  << "    down = " << MuMu_Mass_muScaleDown 
				//					  << std::endl;
				MuMu_Pt = dimuonPt;
				MuMu_Pt_muScaleUp = dimuonPt_muScaleUp;
				MuMu_Pt_muScaleDown = dimuonPt_muScaleDown;
				//				std::cout << "MuMu_Pt = " << MuMu_Pt
				//					  << "    up = " << MuMu_Pt_muScaleUp
				//					  << "    down = " << MuMu_Pt_muScaleDown 
				//					  << std::endl;
				
				MuMu_DR = drMuMu;
				MuMuMET_DPhi = dPhiMuMuMet;
				TrkTrk_Mass = ditrackMass;
				TrkTrk_Pt = ditrackPt;
				TrkTrk_DR = drTrkTrk;
				TrkTrkMET_DPhi = dPhiTrkTrkMet;
				MuMuTrkTrk_Mass = visibleMass;
				MuMuTrkTrk_Pt = visiblePt;
				MuMuTauTau_Mass = totalMass;
				MuMuTauTau_Pt = totalPt;
				MET_Pt = metPt;
				MuMuTrkTrkMET_Mass = visibleandMetMass;
				MuMuTrkTrkMET_Mass_muScaleUp = visibleandMetMass_muScaleUp;
				MuMuTrkTrkMET_Mass_muScaleDown = visibleandMetMass_muScaleDown;
				MuMuTrkTrkMET_Mass_JetEnUp = visibleandMetMass_JetEnUp;
				MuMuTrkTrkMET_Mass_JetEnDown = visibleandMetMass_JetEnDown;
				MuMuTrkTrkMET_Mass_UnclEnUp = visibleandMetMass_UnclEnUp;
				MuMuTrkTrkMET_Mass_UnclEnDown = visibleandMetMass_UnclEnDown;

				//				std::cout << "MuMuTrkTrkMET_Mass = " << MuMuTrkTrkMET_Mass
				//					  << "  muScaleUp = " << MuMuTrkTrkMET_Mass_muScaleUp
				//					  << "  muScaleDown = " << MuMuTrkTrkMET_Mass_muScaleDown 
				//					  << "  JetEnUp = " << MuMuTrkTrkMET_Mass_JetEnUp
				//					  << "  JetEnDown = " << MuMuTrkTrkMET_Mass_JetEnDown 
				//					  << "  UnclEnUp = " << MuMuTrkTrkMET_Mass_UnclEnUp
				//					  << "  UnclEnDown = " << MuMuTrkTrkMET_Mass_UnclEnDown 
				//					  << std::endl;
				Eventweight = weight;

				if (controlNNNN)
				{
					tree_NNNN[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00NN)
				{
					tree_00NN[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00SemiIso)
				{
					tree_00SemiIso[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (controlSoftIso)
				{
					tree_SoftIso[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (controlVerySoftIso)
				{
				  float mumu_DR = TMath::Min(float(1.5),MuMu_DR);
				  float tautau_DR = TMath::Min(float(1.5),TrkTrk_DR);
				  float weight_MuMu_DR = a0_MuMu_DR + a1_MuMu_DR*mumu_DR;
				  float weight_TrkTrk_DR = a0_TrkTrk_DR + a1_TrkTrk_DR*tautau_DR;
				  Correctionweight = weight_MuMu_DR * weight_TrkTrk_DR;
				  TrkTrk_DR_weight_Down = Correctionweight/weight_TrkTrk_DR;
				  MuMu_DR_weight_Down = Correctionweight/weight_MuMu_DR;
				  TrkTrk_DR_weight_Up = Correctionweight*weight_TrkTrk_DR;
				  MuMu_DR_weight_Up = Correctionweight*weight_MuMu_DR;
				  /*
				  std::cout << "Correction weight = " << Correctionweight << std::endl;
				  std::cout << "TrkTrk_DR_weight_Down = " << TrkTrk_DR_weight_Down << std::endl; 
				  std::cout << "TrkTrk_DR_weight_Up = " << TrkTrk_DR_weight_Up << std::endl; 
				  std::cout << "MuMu_DR_weight_Down = " << MuMu_DR_weight_Down << std::endl; 
				  std::cout << "MuMu_DR_weight_Up = " << MuMu_DR_weight_Up << std::endl; 
				  */
				  tree_VerySoftIso[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00SoftIso)
				{
					tree_00SoftIso[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00VerySoftIso)
				{
					tree_00VerySoftIso[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (controlNN00)
				{
					tree_NN00[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00_0SemiIsoorSemiIso0)
				{
					tree_00_0SemiIsoorSemiIso0[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (control00_0NorN0)
				{
					tree_00_0NorN0[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();
				}

				if (signalRegion)
				{

					//event info
					//std::cout << event_run <<":"<< event_luminosityblock <<":"<< event_nr << std::endl;

					////////////////////Reweighting Track and Muon Iso/////////////////////////////
					double scaleFIso;
					if (year == 2016) scaleFIso = 0.933245;
					else if (year == 2017) scaleFIso = 0.963655;
					else if (year == 2018) scaleFIso = 0.909517;
					else
					{
						std::cout << "year is not 2016, 2017, 2018 - exiting" << '\n';
						exit(-1);
					}

					if (!isData) weight *= scaleFIso * scaleFIso;
					////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////

					//////// now filling histograms	///////
					Eventweight = weight;
					treeSel[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill();

					counter_FinalEventsCatH[iInt][IDAtautau1candidate][IDAtautau2candidate]->Fill(1., weight);
					counter_FinalEventsH->Fill(1., weight);
				}

				//============================================================================================================//
			}	/////end loop over all Intervals (Categories) used for Mu-Trk Pair

		}	// icand loop

		delete tree_;
		file_->Close();
		delete file_;
	}	// filelist loop

	file->cd("");
	file->Write();
	file->Close();

	//delete file;
}	// int main loop
