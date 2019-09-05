#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"



int main(int argc, char * argv[]) {

	// first argument - config file 
	// second argument - filelist

	using namespace std;

	// **** configuration
	Config cfg(argv[1]);
	string Channel="muel";


	bool fillplots= true;
	//  const bool applyInclusiveSelection = cfg.get<bool>("ApplyInclusiveSelection");
	const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
	const bool applyPUreweighting_vertices = cfg.get<bool>("ApplyPUreweighting_vertices");
	const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");
	const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
	const bool isData = cfg.get<bool>("IsData");
	const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
	// kinematic cuts on electrons
	const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
	const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
	const float etaElectronCut     = cfg.get<float>("etaElectronCut");
	const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
	const float dzElectronCut      = cfg.get<float>("dzElectronCut");
	const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
	const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
	const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");
	const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
	const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");

	// veto electrons
	const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
	const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
	const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
	const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
	const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
	const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

	// kinematic cuts on muons
	const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
	const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
	const float etaMuonCut     = cfg.get<float>("etaMuonCut");
	const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
	const float dzMuonCut      = cfg.get<float>("dzMuonCut");
	const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
	const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
	const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
	const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
	const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");

	// veto muons
	const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
	const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
	const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
	const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
	const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
	const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

	// vertex cuts
	const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
	const float zVertexCut     = cfg.get<float>("ZVertexCut");
	const float dVertexCut     = cfg.get<float>("DVertexCut");


	// topological cuts
	const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
	const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
	const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
	const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");

	// jets
	const double etaJetCut   = cfg.get<double>("etaJetCut");
	const double ptJetCut   = cfg.get<double>("ptJetCut");


	//triggers
	TString LowPtLegElectron(lowPtLegElectron);
	TString HighPtLegElectron(highPtLegElectron);

	TString LowPtLegMuon(lowPtLegMuon);
	TString HighPtLegMuon(highPtLegMuon);


	const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
	const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");

	const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
	const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");

	const string Electron17TriggerFile = cfg.get<string>("Electron17TriggerEff");
	const string Electron12TriggerFile = cfg.get<string>("Electron12TriggerEff");

	const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");

	const string dataBaseDir = cfg.get<string>("DataBaseDir");
	// **** end of configuration


	const string jsonFile = cfg.get<string>("jsonFile");
	// vertex distributions filenames and histname
	const string vertDataFileName = cfg.get<string>("VertexDataFileName");
	const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
	const string vertHistName     = cfg.get<string>("VertexHistName");

	// lepton scale factors
	const string muonSfDataBarrel = cfg.get<string>("MuonSfDataBarrel");
	const string muonSfDataEndcap = cfg.get<string>("MuonSfDataEndcap");
	const string muonSfMcBarrel = cfg.get<string>("MuonSfMcBarrel");
	const string muonSfMcEndcap = cfg.get<string>("MuonSfMcEndcap");

	const double Lumi   = cfg.get<double>("Lumi");
	const double bTag   = cfg.get<double>("bTag");
	const double metcut = cfg.get<double>("metcut");

	const string Region  = cfg.get<string>("Region");
	const string Sign  = cfg.get<string>("Sign");
	string cmsswBase = (getenv ("CMSSW_BASE"));
	string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;


	// Run-lumi selector
	std::vector<Period> periods;  
	if (isData) { // read the good runs 
		std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
		if (inputFileStream.fail() ) {
			std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
			std::cout << "please check" << std::endl;
			std::cout << "quitting program" << std::endl;
			exit(-1);
		}

		for(std::string s; std::getline(inputFileStream, s); ) {
			//std::fstream inputFileStream("temp", std::ios::in);
			periods.push_back(Period());
			std::stringstream ss(s);
			ss >> periods.back();
		}
	}

	// CutList

	CutList.clear();
	CutList.push_back("No cut");
	CutList.push_back("No cut after PU");
	CutList.push_back("Trigger");
	CutList.push_back("OS$\\mu-e$");
	CutList.push_back("3rd $\\ell$-Veto");
	CutList.push_back("Lepton SF");
	CutList.push_back("topPtRwgt");
	CutList.push_back("Cleaned Jets");
	CutList.push_back("$ E_T^{\\rm miss}>$ 100");
	CutList.push_back("Jets $>$2");
	CutList.push_back("$>0>b-tag");
	CutList.push_back("$1.5<\\Delta R<4$");


	int CutNumb = int(CutList.size());
	xs=1;fact=1;fact2=1;

	unsigned int RunMin = 9999999;
	unsigned int RunMax = 0;

	ifstream ifs("xsecs");
	string line;

	while(std::getline(ifs, line)) // read one line from ifs
	{

		fact=fact2=1;
		istringstream iss(line); // access line as a stream

		// we only need the first two columns
		string dt,st1,st2;st1="stau2_1";st2="stau5_2";
		iss >> dt >> xs >> fact >> fact2;
		//ifs >> dt >> xs; // no need to read further
		//cout<< " "<<dt<<"  "<<endl;
		//cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
		//if (dt==argv[2]) {
		//if (std::string::npos != dt.find(argv[2])) {
		if (  dt == argv[2]) {
			XSec= xs*fact*fact2;
			cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
		}
		/*
		   if ( argv[2] == st1) {ChiMass=100;mIntermediate=200;}
		   else if (argv[2] == st2) {ChiMass=200;mIntermediate=500;}
		   */
		if (isData) XSec=1.;
		ChiMass=0.0;
	}

	if (XSec<0&& !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}

	std::vector<unsigned int> allRuns; allRuns.clear();

	cout<<" ChiMass is "<<ChiMass<<"  "<<mIntermediate<<endl;
	bool extraelec_veto = true;
	bool extramuon_veto = true;
	char ff[100];


	sprintf(ff,"%s/%s",argv[3],argv[2]);

	if (applyPUreweighting_vertices and applyPUreweighting_official) 

	{std::cout<<"ERROR: Choose only ONE PU reweighting method (vertices or official, not both!) " <<std::endl; exit(-1);}

	// reweighting with vertices

	// reading vertex weights
	TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertDataFileName);
	TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertMcFileName);

	TH1D * vertexDataH = (TH1D*)fileDataNVert->Get(TString(vertHistName));
	TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get(TString(vertHistName));

	float normVertexData = vertexDataH->GetSumOfWeights();
	float normVertexMc   = vertexMcH->GetSumOfWeights();


	vertexDataH->Scale(1/normVertexData);
	vertexMcH->Scale(1/normVertexMc);

	PileUp * PUofficial = new PileUp();

	TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Nov17.root","read");
	TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring15_PU25_Startup.root", "read");
	TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
	TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
	PUofficial->set_h_data(PU_data);
	PUofficial->set_h_MC(PU_mc);


	TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataBarrel);  // mu SF barrel data
	TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataEndcap); // mu SF endcap data
	TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcBarrel);  // mu SF barrel MC
	TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcEndcap); // mu SF endcap MC 

	TGraphAsymmErrors *hEffBarrelData = (TGraphAsymmErrors*)f10->Get("ZMassBarrel");
	TGraphAsymmErrors *hEffEndcapData = (TGraphAsymmErrors*)f11->Get("ZMassEndcap");
	TGraphAsymmErrors *hEffBarrelMC = (TGraphAsymmErrors*)f12->Get("ZMassBarrel");
	TGraphAsymmErrors *hEffEndcapMC = (TGraphAsymmErrors*)f13->Get("ZMassEndcap");

	double * dataEffBarrel = new double[10];
	double * dataEffEndcap = new double[10];
	double * mcEffBarrel = new double[10];
	double * mcEffEndcap = new double[10];

	dataEffBarrel = hEffBarrelData->GetY();
	dataEffEndcap = hEffEndcapData->GetY();
	mcEffBarrel = hEffBarrelMC->GetY();
	mcEffEndcap = hEffEndcapMC->GetY();


	cout<<" Will initialize lepton SF...."<<endl;
	// Lepton Scale Factors 

	TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);

	cout<<TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile)<<endl;

	ScaleFactor * SF_muonIdIso = new ScaleFactor();
	SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));

	cout<<" Will initialize lepton SF 1...."<<endl;
	ScaleFactor * SF_muon17 = new ScaleFactor();
	SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon17TriggerFile));


	cout<<" Will initialize lepton SF 2...."<<endl;
	ScaleFactor * SF_muon8 = new ScaleFactor();
	SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));

	cout<<" Will initialize lepton SF 3...."<<endl;
	ScaleFactor * SF_electronIdIso = new ScaleFactor();
	SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));

	cout<<" Will initialize lepton SF 4...."<<endl;
	ScaleFactor * SF_electron17 = new ScaleFactor();
	SF_electron17->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron17TriggerFile));

	ScaleFactor * SF_electron12 = new ScaleFactor();
	SF_electron12->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron12TriggerFile));

	cout<<" Initialized lepton SF...."<<endl;

	ScaleFactor * SF_TFR; 
	bool applyTFR = false;
	if (applyTFR) {
		SF_TFR = new ScaleFactor();
		SF_TFR->init_ScaleFactorb(TString(cmsswBase)+"/src/"+TString(TauFakeRateFile),applyTFR);
	}

	////////


	double Weight=0;
	int nTotalFiles = 0;
	int iCut=0;
	double CFCounter[CutNumb];
	double statUnc[CutNumb];
	int iCFCounter[CutNumb];
	for (int i=0;i < CutNumb; i++){
		CFCounter[i] = 0;
		iCFCounter[i] = 0;
		statUnc[i] =0;
	}
	// file name and tree name
	std::string rootFileName(argv[2]);
	//std::ifstream fileList(argv[2]);
	std::ifstream fileList(ff);
	//std::ifstream fileList0(argv[2]);
	std::ifstream fileList0(ff);
	std::string ntupleName("makeroottree/AC1B");
	std::string initNtupleName("initroottree/AC1B");

	TString era=argv[3];

	TString TStrName(rootFileName+"_"+Region+"_"+Sign);
	std::cout <<" The filename will be "<<TStrName <<std::endl;  

	// output fileName with histograms
	TFile * file;
	if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
	if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
	file->mkdir(Channel.c_str());
	file->cd(Channel.c_str());



	/*  


	// qcd weight 
	TFile * fileQCDdRLt2 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_dRLt2.root","read");
	TH2D * qcdHistWeightsdRLt2H = (TH2D*)fileQCDdRLt2->Get("QCDraio");
	TFile * fileQCDdR2to4 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_dR2to4.root","read");
	TH2D * qcdHistWeightsdR2to4H = (TH2D*)fileQCDdR2to4->Get("QCDraio");
	TFile * fileQCDdRGt4 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_dRGt4.root","read");
	TH2D * qcdHistWeightsdRGt4H = (TH2D*)fileQCDdRGt4->Get("QCDraio");
	// qcd weight DZeta cut
	TFile * fileQCDDZetadRLt2 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_DZeta_dRLt2.root","read");
	TH2D * qcdHistWeightsDZetadRLt2H = (TH2D*)fileQCDDZetadRLt2->Get("QCDraio");
	TFile * fileQCDDZetadR2to4 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_DZeta_dR2to4.root","read");
	TH2D * qcdHistWeightsDZetadR2to4H = (TH2D*)fileQCDDZetadR2to4->Get("QCDraio");
	TFile * fileQCDDZetadRGt4 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/QCD_IsoE_AntiIsoMu_DZeta_dRGt4.root","read");
	TH2D * qcdHistWeightsDZetadRGt4H = (TH2D*)fileQCDDZetadRGt4->Get("QCDraio");

	int nBinsQCD = qcdHistWeightsdRLt2H->GetNbinsX();
	float qcdMin = qcdHistWeightsdRLt2H->GetXaxis()->GetBinLowEdge(1);
	float qcdMax = qcdHistWeightsdRLt2H->GetXaxis()->GetBinLowEdge(1+nBinsQCD);

	std::cout << "nBins(QCD) = " << nBinsQCD 
	<< "    min(QCD) = " << qcdMin
	<< "    max(QCD) = " << qcdMax << std::endl; 
	//  exit(-1);
	*/

	int nFiles = 0;
	int nEvents = 0;
	int selEvents = 0;

	int selEventsAllMuons = 0;
	int selEventsIdMuons = 0;
	int selEventsIsoMuons = 0;
	bool lumi=false;
	bool isLowIsoMu=false;
	bool isHighIsoMu = false;
	bool isLowIsoTau=false;
	bool isHighIsoTau = false;


	std::string dummy;
	// count number of files --->
	while (fileList0 >> dummy) nTotalFiles++;
	SetupTree();
	SetupHists(CutNumb); 
	if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
	//if (nTotalFiles>50) nTotalFiles=50;
	//nTotalFiles = 10;
	for (int iF=0; iF<nTotalFiles; ++iF) {

		std::string filen;
		fileList >> filen;

		std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
		TFile * file_ = TFile::Open(TString(filen));


		TH1D * histoInputEvents = NULL;
		histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
		if (histoInputEvents==NULL) continue;
		int NE = int(histoInputEvents->GetEntries());
		for (int iE=0;iE<NE;++iE)
			inputEventsH->Fill(0.);
		std::cout << "      number of input events         = " << NE << std::endl;

		TTree * _inittree = NULL;
		_inittree = (TTree*)file_->Get(TString(initNtupleName));
		if (_inittree==NULL) continue;
		Float_t genweight;
		if (!isData)
			_inittree->SetBranchAddress("genweight",&genweight);
		Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
		std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
		for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
			_inittree->GetEntry(iEntry);
			if (isData)
				histWeightsH->Fill(0.,1.);
			else
				histWeightsH->Fill(0.,genweight);
		}

		TTree * _tree = NULL;
		_tree = (TTree*)file_->Get(TString(ntupleName));
		if (_tree==NULL) continue;
		Long64_t numberOfEntries = _tree->GetEntries();
		std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
		AC1B analysisTree(_tree);



		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

			Float_t weight = 1;
			Float_t puweight = 1;
			//float topptweight = 1;
			analysisTree.GetEntry(iEntry);
			nEvents++;

			iCut = 0;

			//std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

			if (nEvents%50000==0) 
				cout << "      processed " << nEvents << " events" << endl; 

			if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
			if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
			double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
					analysisTree.primvertex_y*analysisTree.primvertex_y);
			if (dVertex>dVertexCut) continue;
			if (analysisTree.primvertex_count<2) continue;  

			//isData= false;
			bool lumi=false;
			isLowIsoMu=false;
			isHighIsoMu = false;
			isLowIsoTau=false;
			isHighIsoTau = false;

			Float_t genweights;
			float topPt = 0;
			float antitopPt = 0;
			LSF_weight = 1.;
			TFR_weight = 1.;
			top_weight = 1.;
			all_weight = 1.;
			pu_weight = 1.;
			gen_weight = 1.;
			trig_weight = 1.;
			bool isZTT = false;

			if(!isData) 
			{
				/*   TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");

				     genweightsTree->SetBranchAddress("genweight",&genweights);
				     Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
				     for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) { 
				     genweightsTree->GetEntry(iEntryInit);
				     histWeightsH->Fill(0.,genweights);
				     }
				     */ /*
					   for (unsigned int igent=0; igent < analysisTree.gentau_count; ++igent) {
					   if (analysisTree.gentau_isPrompt[igent]) isZTT = true; 
					   }
					   */
			if (!isData && ( string::npos != filen.find("TTJets")  || string::npos != filen.find("TTPowHeg")) ) {
				for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
					// 		cout<< "  info = " <<  int(analysisTree.genparticles_count) <<"  "<<int(analysisTree.genparticles_pdgid[igen])<<endl;

					if (analysisTree.genparticles_pdgid[igen]==6)
						topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
								analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
					if (analysisTree.genparticles_pdgid[igen]==-6)
						antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
								analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);

				}    
			}
				weight *= analysisTree.genweight;
				gen_weight *=analysisTree.genweight;
				lumi=true;
				//cout<<"  weight from init "<<genweights<< "  "<<analysisTree.genweight<<"  "<<weight<<endl;
			}



			if (isData)  {
				XSec = 1.;
				histRuns->Fill(analysisTree.event_run);
				///////////////according to dimuons
				int n=analysisTree.event_run;
				int lum = analysisTree.event_luminosityblock;

				std::string num = std::to_string(n);
				std::string lnum = std::to_string(lum);
				for(const auto& a : periods)
				{

					if ( num.c_str() ==  a.name ) {
						//std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
						//     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

						for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

							//	cout<<b->lower<<"  "<<b->bigger<<endl;
							if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
						}
						auto last = std::prev(a.ranges.end());
						//    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
						if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


					}

				}

				if (!lumi) continue;
				//if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
			}


			if (analysisTree.event_run<RunMin)
				RunMin = analysisTree.event_run;

			if (analysisTree.event_run>RunMax)
				RunMax = analysisTree.event_run;

			//std::cout << " Run : " << analysisTree.event_run << std::endl;

			bool isNewRun = true;
			if (allRuns.size()>0) {
				for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
					if (analysisTree.event_run==allRuns.at(iR)) {
						isNewRun = false;
						break;
					}
				}
			}

			if (isNewRun) 
				allRuns.push_back(analysisTree.event_run);




			if (!lumi) continue;


			JetsMV.clear();
			ElMV.clear();
			TauMV.clear();
			MuMV.clear();
			LeptMV.clear();
			mu_index=-1;
			tau_index=-1;
			el_index=-1;


			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;



			double MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

			METV.SetPx(analysisTree.pfmet_ex);	      
			METV.SetPy(analysisTree.pfmet_ey);


			for (unsigned int ijj = 0; ijj<analysisTree.pfjet_count; ++ijj) {
				JetsV.SetPxPyPzE(analysisTree.pfjet_px[ijj], analysisTree.pfjet_py[ijj], analysisTree.pfjet_pz[ijj], analysisTree.pfjet_e[ijj]);
				JetsMV.push_back(JetsV);
			} 


			for (unsigned int imm = 0; imm<analysisTree.muon_count; ++imm) {
				MuV.SetPtEtaPhiM(analysisTree.muon_pt[imm], analysisTree.muon_eta[imm], analysisTree.muon_phi[imm], muonMass);
				MuMV.push_back(MuV);
				//	mu_index=0;
			}

			for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
				ElV.SetPtEtaPhiM(analysisTree.electron_pt[ie], analysisTree.electron_eta[ie], analysisTree.electron_phi[ie], electronMass);
				ElMV.push_back(ElV);
				//	el_index=0;
			}

			for (unsigned int itt = 0; itt<analysisTree.tau_count; ++itt) {
				TauV.SetPtEtaPhiM(analysisTree.tau_pt[itt], analysisTree.tau_eta[itt], analysisTree.tau_phi[itt], tauMass);
				TauMV.push_back(TauV);
				//	tau_index=0;
			}

			if (!isData) 
			{
				if (applyPUreweighting)	 {
					puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
					weight *=puweight;      
					pu_weight = puweight;
				}
			}



			// vector <string> ss; ss.push_back(.c_str());
			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;





			unsigned int nLowPtLegElectron = 0;
			bool isLowPtLegElectron = false;

			unsigned int nHighPtLegElectron = 0;
			bool isHighPtLegElectron = false;

			unsigned int nLowPtLegMuon = 0;
			bool isLowPtLegMuon = false;

			unsigned int nHighPtLegMuon = 0;
			bool isHighPtLegMuon = false;

			unsigned int nfilters = analysisTree.run_hltfilters->size();
			//      std::cout << "nfiltres = " << nfilters << std::endl;
			for (unsigned int i=0; i<nfilters; ++i) {
				//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
				TString HLTFilter(analysisTree.run_hltfilters->at(i));
				if (HLTFilter==LowPtLegElectron) {
					nLowPtLegElectron = i;
					isLowPtLegElectron = true;
				}
				if (HLTFilter==HighPtLegElectron) {
					nHighPtLegElectron = i;
					isHighPtLegElectron = true;
				}
				if (HLTFilter==LowPtLegMuon) {
					nLowPtLegMuon = i;
					isLowPtLegMuon = true;
				}
				if (HLTFilter==HighPtLegMuon) {
					nHighPtLegMuon = i;
					isHighPtLegMuon = true;
				}
			}
			if (!isLowPtLegElectron) {
				std::cout << "HLT filter " << LowPtLegElectron << " not found" << std::endl;
				exit(-1);
			}
			if (!isHighPtLegElectron) {
				std::cout << "HLT filter " << HighPtLegElectron << " not found" << std::endl;
				exit(-1);
			}
			if (!isLowPtLegMuon) {
				std::cout << "HLT filter " << LowPtLegMuon << " not found" << std::endl;
				exit(-1);
			}
			if (!isHighPtLegMuon) {
				std::cout << "HLT filter " << HighPtLegMuon << " not found" << std::endl;
				exit(-1);
			}/*
			    unsigned int nBTagDiscriminant = 0;
			    for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
			    TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
			    if (discr=BTagDiscriminator)
			    nBTagDiscriminant = iBTag;
			    }
			    */


			//cout<<" The trigger is ok  "<<isLowPtLegElectron<<" "<<isHighPtLegElectron<<"  "<<isHighPtLegMuon<<"  "<<isLowPtLegMuon<<endl;


			// electron selection
			vector<int> electrons; electrons.clear();
			for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
				if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
				if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
				if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
				if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
				//	bool electronMvaId = electronMvaIdWP80(analysisTree.electron_pt[ie],
				//					       analysisTree.electron_superclusterEta[ie],
				//					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
				bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie];
				if (!electronMvaId&&applyElectronId) continue;
				if (!analysisTree.electron_pass_conversion[ie]&&applyElectronId) continue;
				if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyElectronId) continue;
				electrons.push_back(ie);
				el_index=int(ie);
			}

			// muon selection
			vector<int> muons; muons.clear();
			for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
				if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
				if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
				if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
				if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
				if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
				muons.push_back(im);
				mu_index=int(im);
			}


			//          cout << "  SelEle=" << electrons.size() 
			//    	   << "  SelMu=" << muons.size() << std::endl;

			if (electrons.size()==0) continue;
			if (muons.size()==0) continue;



			// selecting muon and electron pair (OS or SS);

			float isoMuMin = 1e+10;
			float isoEleMin = 1e+10;
			//      if (muons.size()>1||electrons.size()>1)
			//      std::cout << "muons = " << muons.size() << "  electrons = " << electrons.size() << std::endl;
			bool isMuon17matched = false;
			bool isMuon8matched  = false;
			bool isElectron17matched = false;
			bool isElectron12matched = false;
			for (unsigned int im=0; im<muons.size(); ++im) {
				bool isMu17 = false;
				bool isMu8 = false;
				int mIndex  = muons.at(im);
				float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
				float photonIsoMu = analysisTree.muon_photonIso[mIndex];
				float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
				float puIsoMu = analysisTree.muon_puIso[mIndex];
				if (isIsoR03) {
					neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[mIndex];
					photonIsoMu = analysisTree.muon_r03_sumPhotonEt[mIndex];
					chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[mIndex];
					puIsoMu = analysisTree.muon_r03_sumPUPt[mIndex];
				}
				float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
				neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
				float absIsoMu = chargedHadIsoMu + neutralIsoMu;
				float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];
				for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
					float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
							analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
					if (dRtrig<deltaRTrigMatch) {
						if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]&&
								analysisTree.muon_pt[mIndex]>ptMuonHighCut) { // Mu17 Leg
							isMu17 = true;
						}
						if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]&&
								analysisTree.muon_pt[mIndex]>ptMuonLowCut) { // Mu8 Leg
							isMu8 = true;
						}
					}
				}

				if (applyTriggerMatch && (!isMu17) && (!isMu8)) continue;

				for (unsigned int ie=0; ie<electrons.size(); ++ie) {

					unsigned int eIndex = electrons.at(ie);

					float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
							analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

					if (dR<dRleptonsCut) continue;

					bool isEle17 = false;
					bool isEle12 = false;

					for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
						float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
								analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
						if (dRtrig<deltaRTrigMatch) {
							if (analysisTree.trigobject_filters[iT][nHighPtLegElectron]&&
									analysisTree.electron_pt[eIndex]>ptElectronHighCut) { // Ele17 Leg
								isEle17 = true;
							}
							if (analysisTree.trigobject_filters[iT][nLowPtLegElectron]&&
									analysisTree.electron_pt[eIndex]>ptElectronLowCut) { // Ele12 Leg
								isEle12 = true;
							}
						}
					}

					bool trigMatch = (isMu17&&isEle12) || (isMu8&&isEle17);
					//	  std::cout << "Trigger match = " << trigMatch << std::endl;

					if (applyTriggerMatch && !trigMatch) continue;

					if(fillplots)
						FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
					CFCounter[iCut]+= weight;
					iCFCounter[iCut]++;
					iCut++;


					//	  bool isKinematicMatch = false;
					//	  if (isMu17&&isEle12) {
					//	    if (analysisTree.muon_pt[mIndex]>ptMuonHighCut&&analysisTree.electron_pt[eIndex]>ptElectronLowCut)
					//	      isKinematicMatch = true;
					//	  }
					//	  if (isMu8&&isEle17) {
					//            if (analysisTree.muon_pt[mIndex]>ptMuonLowCut&&analysisTree.electron_pt[eIndex]>ptElectronHighCut)
					//              isKinematicMatch = true;
					//          }
					//	  if (!isKinematicMatch) continue;

					float neutralHadIsoEle = analysisTree.electron_neutralHadIso[eIndex];
					float photonIsoEle = analysisTree.electron_photonIso[eIndex];
					float chargedHadIsoEle = analysisTree.electron_chargedHadIso[eIndex];
					float puIsoEle = analysisTree.electron_puIso[eIndex];
					if (isIsoR03) {
						neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
						photonIsoEle = analysisTree.electron_r03_sumPhotonEt[eIndex];
						chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
						puIsoEle = analysisTree.electron_r03_sumPUPt[eIndex];
					}
					float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
					neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
					float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
					float relIsoEle = absIsoEle/analysisTree.electron_pt[eIndex];

					if (int(mIndex)!=int(mu_index)) {
						if (relIsoMu==isoMuMin) {
							if (analysisTree.muon_pt[mIndex]>analysisTree.muon_pt[mu_index]) {
								isoMuMin  = relIsoMu;
								mu_index = int(mIndex);
								isoEleMin = relIsoEle;
								el_index = int(eIndex);
								isMuon17matched = isMu17;
								isMuon8matched = isMu8;
								isElectron17matched = isEle17;
								isElectron12matched = isEle12;
							}
						}
						else if (relIsoMu<isoMuMin) {
							isoMuMin  = relIsoMu;
							mu_index = int(mIndex);
							isoEleMin = relIsoEle;
							el_index = int(eIndex);
							isMuon17matched = isMu17;
							isMuon8matched = isMu8;
							isElectron17matched = isEle17;
							isElectron12matched = isEle12;
						}
					}
					else {
						if (relIsoEle==isoEleMin) {
							if (analysisTree.electron_pt[eIndex]>analysisTree.electron_pt[el_index]) {
								isoEleMin = relIsoEle;
								el_index = int(eIndex);
								isElectron17matched = isEle17;
								isElectron12matched = isEle12;
							}
						}
						else if (relIsoEle<isoEleMin) {
							isoEleMin = relIsoEle;
							el_index = int(eIndex);
							isElectron17matched = isEle17;
							isElectron12matched = isEle12;
						}
					}

				}
			}


			if ((int)el_index<0) continue;
			if ((int)mu_index<0) continue;
			//    cout << "im = " << mu_index << "   eIndex = " << el_index << std::endl;
			//      std::cout << "Post synch selection " << std::endl;
			//      std::cout << std::endl;



			bool os = (analysisTree.muon_charge[mu_index]*analysisTree.electron_charge[el_index]) < 0;

			if (!os) continue;
			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;



			//////triger weight
			double ptMu1 = (double)analysisTree.muon_pt[mu_index];
			double etaMu1 = (double)analysisTree.muon_eta[mu_index];

			double ptEl1 = (double)analysisTree.electron_pt[el_index];
			double etaEl1 = (double)analysisTree.electron_eta[el_index];


			float trigweight=1.;
			float trigweight_1=1.;
			float trigweight_2=1.;

			float Ele17EffData = (float)SF_electron17->get_EfficiencyData(double(ptEl1),double(etaEl1));
			float Ele17EffMC   = (float)SF_electron17->get_EfficiencyMC(double(ptEl1),double(etaEl1));

			float Ele12EffData = (float)SF_electron12->get_EfficiencyData(double(ptEl1),double(etaEl1));
			float Ele12EffMC   = (float)SF_electron12->get_EfficiencyMC(double(ptEl1),double(etaEl1));

			float Mu17EffData = (float)SF_muon17->get_EfficiencyData(double(ptMu1),double(etaMu1));
			float Mu17EffMC   = (float)SF_muon17->get_EfficiencyMC(double(ptMu1),double(etaMu1));

			float Mu8EffData = (float)SF_muon8->get_EfficiencyData(double(ptMu1),double(etaMu1));
			float Mu8EffMC   = (float)SF_muon8->get_EfficiencyMC(double(ptMu1),double(etaMu1));

			float trigWeightData = Mu17EffData*Ele12EffData + Mu8EffData*Ele17EffData - Mu17EffData*Ele17EffData;
			float trigWeightMC   = Mu17EffMC*Ele12EffMC     + Mu8EffMC*Ele17EffMC     - Mu17EffMC*Ele17EffMC;

			if (isMuon17matched && isElectron12matched) {
				trigweight_1 = (float)SF_electron12->get_ScaleFactor(double(ptEl1),double(etaEl1));
				trigweight_2 = (float)SF_muon17->get_ScaleFactor(double(ptMu1),double(etaMu1));
			}
			else if (isMuon8matched && isElectron17matched) {
				trigweight_1 = (float)SF_electron17->get_ScaleFactor(double(ptEl1),double(etaEl1));
				trigweight_2 = (float)SF_muon8->get_ScaleFactor(double(ptMu1),double(etaMu1));
			}


			if (!isData) {
				if (trigWeightMC>1e-6)
					trigweight = trigWeightData / trigWeightMC;
				weight *= trigweight;
				trig_weight = trigweight;
			}

			/*     if(fillplots)
			       FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			       CFCounter[iCut]+= weight;
			       iCFCounter[iCut]++;
			       iCut++;
			       */




			// looking for extra electron
			bool foundExtraElectron = false;
			for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
				if (int(ie)==int(el_index)) continue;
				if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
				if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
				if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
				if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
				//	bool electronMvaId = electronMvaIdWP90(analysisTree.electron_pt[ie],
				//					       analysisTree.electron_superclusterEta[ie],
				//					       analysisTree.electron_mva_id_nontrigPhys14[ie]);
				bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
				if (!electronMvaId&&applyVetoElectronId) continue;
				if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId&&applyVetoElectronId) continue;
				if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId&&applyVetoElectronId) continue;
				float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
				float photonIsoEle = analysisTree.electron_photonIso[ie];
				float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
				float puIsoEle = analysisTree.electron_puIso[ie];
				if (isIsoR03) {
					neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
					photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
					chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
					puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
				}
				float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
				neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
				float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
				float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];	
				if (relIsoEle>isoVetoElectronCut) continue;
				foundExtraElectron = true;
			}

			// looking for extra muon
			bool foundExtraMuon = false;
			for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
				if (int(im)==int(mu_index)) continue;
				if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
				if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
				if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
				if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
				if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
				float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
				float photonIsoMu = analysisTree.muon_photonIso[im];
				float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
				float puIsoMu = analysisTree.muon_puIso[im];
				if (isIsoR03) {
					neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
					photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
					chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
					puIsoMu = analysisTree.muon_r03_sumPUPt[im];
				}
				float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
				neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
				float absIsoMu = chargedHadIsoMu + neutralIsoMu;
				float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
				if (relIsoMu>isoVetoMuonCut) continue;
				foundExtraMuon = true;
			}

			//dilepton_veto = false;
			extraelec_veto = foundExtraElectron;
			extramuon_veto = foundExtraMuon;

			bool applyInclusiveSelection = true;
			if (applyInclusiveSelection) {
				if (extraelec_veto) continue;
				if (extramuon_veto) continue;
				if (isoMuMin>isoMuonHighCut) continue;
				if (isoEleMin>isoElectronHighCut) continue;
			}

			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;


			double IdIsoSF_mu1 = SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
			double IdIsoSF_El1 = SF_electronIdIso->get_ScaleFactor(ptEl1, etaEl1);




			if (!isData && applyLeptonSF) {

				MuSF_IdIso_Mu1H->Fill(IdIsoSF_mu1);
				//	  if (ptMu1<20||ptMu2<20) {
				//	    std::cout << "mu 1 ->  pt = " << ptMu1 << "   eta = " << etaMu1 << std::endl;
				//	    std::cout << "eff data mu 1 = " << SF_muonIdIso->get_EfficiencyData(ptMu1, etaMu1)<< " |  eff mc mu 1 = " << SF_muonIdIso->get_EfficiencyMC(ptMu1, etaMu1)<<std::endl;
				//	    std::cout << "mu 2 ->  pt = " << ptMu2 << "   eta = " << etaMu2 << std::endl;
				//	    std::cout << "eff data mu 2 = " << SF_muonIdIso->get_EfficiencyData(ptMu2, etaMu2)<< " |  eff mc mu 2 = " << SF_muonIdIso->get_EfficiencyMC(ptMu2, etaMu2)<<std::endl;
				//	    std::cout << "SF mu1 = " << IdIsoSF_mu1 << std::endl;
				//	    std::cout << "SF mu2 = " << IdIsoSF_mu2 << std::endl;
				//	    
				//	    std::cout << " mass = " << massSel << std::endl;
				//	    std::cout << std::endl;
				//	  }
				weight  *=IdIsoSF_mu1*IdIsoSF_El1;
				LSF_weight = IdIsoSF_mu1*IdIsoSF_El1;
			}

			//cout<<" ok with SF "<<weight<<endl;
			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;


				if (!isData && applyTFR) {

					//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	

					double ptTau1 = (double)analysisTree.tau_pt[tau_index];
					double etaTau1 = (double)analysisTree.tau_eta[tau_index];
					double TFRSF_mu1 = SF_TFR->get_ScaleFactor(ptTau1, etaTau1);

					MuSF_IdIso_Mu1H->Fill(TFRSF_mu1);
					weight *= TFRSF_mu1;
					TFR_weight  = TFRSF_mu1;
					//cout<<"  "<<TFRSF_mu1<<"  for  eta  "<<etaTau1<<  " pT  "<< ptTau1<<endl;
				}




			if (!isData && ( string::npos != filen.find("TTJets")  || string::npos != filen.find("TTPowHeg")) ) 
			{
				if (topPt>0.&&antitopPt>0.) {
					float topptweight = topPtWeight(topPt,antitopPt);
					//cout<<"  "<<topPt<<"  "<<antitopPt<<"  "<<topptweight<<endl;
					weight = weight*topptweight;
					top_weight = topptweight;
				}
			}

			//cout<<" ok with PU "<<weight<<endl;

			if(fillplots)
				FillMainHists(iCut, weight, ElMV, MuMV, TauMV,JetsMV,METV, ChiMass,mIntermediate,analysisTree, Channel, mu_index,el_index,tau_index);
			CFCounter[iCut]+= weight;
			iCFCounter[iCut]++;
			iCut++;

				muon_index = (int)mu_index;
				electron_index = (int)el_index;
				taus_index = (int)tau_index;

				mu_count= (int)analysisTree.muon_count;
				//cout<<" here ============================> "<<iEntry<<"  "<<mu_count<<"  "<<(int)analysisTree.muon_count<<"  "<<analysisTree.muon_count<<endl;
				for (unsigned int im=0;im<analysisTree.muon_count;im++){
				mu_px[im]=analysisTree.muon_px[im];
				mu_py[im]=analysisTree.muon_py[im];
				mu_pz[im]=analysisTree.muon_pz[im];
				mu_eta[im]=analysisTree.muon_eta[im];
				mu_pt[im]=analysisTree.muon_pt[im];
				mu_phi[im]=analysisTree.muon_phi[im];
				mu_charge[im]=analysisTree.muon_charge[im];
				mu_dxy[im]=analysisTree.muon_dxy[im];
				mu_dz[im]=analysisTree.muon_dz[im];
				}
				

				el_count=(int)analysisTree.electron_count;
				for (unsigned int ie=0;ie<analysisTree.electron_count;ie++){
				el_px[ie]=analysisTree.electron_px[ie];
				el_py[ie]=analysisTree.electron_py[ie];
				el_pz[ie]=analysisTree.electron_pz[ie];
				el_eta[ie]=analysisTree.electron_eta[ie];
				el_pt[ie]=analysisTree.electron_pt[ie];
				el_phi[ie]=analysisTree.electron_phi[ie];
				el_charge[ie]=analysisTree.electron_charge[ie];
				el_dxy[ie]=analysisTree.electron_dxy[ie];
				el_dz[ie]=analysisTree.electron_dz[ie];
				}

				
				ta_count=(int)analysisTree.tau_count;
				for (unsigned int it=0;it<analysisTree.tau_count;it++){
				ta_px[it]=analysisTree.tau_px[it];
				ta_py[it]=analysisTree.tau_py[it];
				ta_pz[it]=analysisTree.tau_pz[it];
				ta_eta[it]=analysisTree.tau_eta[it];
				ta_pt[it]=analysisTree.tau_pt[it];
				ta_phi[it]=analysisTree.tau_phi[it];
				ta_charge[it]=analysisTree.tau_charge[it];
				ta_dxy[it]=analysisTree.tau_dxy[it];
				ta_dz[it]=analysisTree.tau_dz[it];
				}
				jet_count=(int)analysisTree.pfjet_count;
				for (unsigned int jj=0;jj<analysisTree.pfjet_count;jj++){

				jet_e[jj] = analysisTree.pfjet_e[jj];
				jet_px[jj] = analysisTree.pfjet_px[jj];
				jet_py[jj] = analysisTree.pfjet_py[jj];
				jet_pz[jj] = analysisTree.pfjet_pz[jj];
				jet_pt[jj] = analysisTree.pfjet_pt[jj];
				jet_eta[jj] = analysisTree.pfjet_eta[jj];
				jet_phi[jj] = analysisTree.pfjet_phi[jj];
				jet_flavour[jj] = analysisTree.pfjet_flavour[jj];
				jet_btag[jj] = analysisTree.pfjet_btag[jj][0];
				}
				met_ex = analysisTree.pfmet_ex;
				met_ey = analysisTree.pfmet_ey;
				met_ez = analysisTree.pfmet_ez;
				met_pt = analysisTree.pfmet_pt;
				met_phi = analysisTree.pfmet_phi;

				all_weight = weight;





			double ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

			double ptScalarSum = -1;


			bool btagged= false;


			JetsMV.clear();

			float jetEtaCut = 2.4;
			float DRmax = 0.5;
			int countjets = 0;
			for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
				float absJetEta = fabs(analysisTree.pfjet_eta[jet]);

				if (absJetEta > etaJetCut) continue;
				if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;


				//double Dr= deltaR(LeptMV.at(il).Eta(), LeptMV.at(il).Phi(),


				bool isPFJetId = false ; 
				isPFJetId =looseJetiD(analysisTree,jet);


				if (!isPFJetId) continue;
				jet_isLoose[jet] = isPFJetId;

				double Drmu=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
						analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
				double Drel=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
						analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
				if (  Drmu  < DRmax)  continue;
				if (  Drel  < DRmax)  continue;


				if (analysisTree.pfjet_btag[jet][0]  > bTag) btagged = true;
				JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
				JetsMV.push_back(JetsV);
				countjets++;
			}	 

				T->Fill();

				continue;

			selEvents++;

			} // end of file processing (loop over events in one file)
			nFiles++;
			delete _tree;
			file_->Close();
			delete file_;
		}

		for(int ci = 0; ci < CutNumb; ci++)
		{
			// tfile << CutList[ci]<<"\t & \t"
			//	    << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;

			CutFlowUnW->SetBinContent(1+ci,0);
			CutFlow->SetBinContent(1+ci,0);
			CutFlowUnW->SetBinContent(1+ci,float(CFCounter[ci]) );


			CFCounter[ci] *= double(XSec*Lumi/( histWeightsH->GetSumOfWeights()));

			CutFlow->SetBinContent(1+ci,float(CFCounter[ci]));

			cout << " i "<<ci<<" "<<iCFCounter[ci]<<"  "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<"  "<<CutFlowUnW->GetBinContent(1+ci)<<"  "<<CutFlow->GetBinContent(1+ci)<<endl;   
			if (iCFCounter[ci] <0.2) statUnc[ci] =0;
			//else statUnc[i] = CFCounter[i]/sqrt(iCFCounter[i]);
			else statUnc[ci] = sqrt(CFCounter[ci]);
		}
		std::cout << std::endl;
		int allEvents = int(inputEventsH->GetEntries());
		std::cout << "Total number of input events    = " << allEvents << std::endl;
		std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
		std::cout << "Total number of selected events = " << selEvents << std::endl;
		std::cout << std::endl;
		file->cd(Channel.c_str());
		hxsec->Fill(XSec);
		hxsec->Write();
		inputEventsH->Write();
		histWeightsH->Write();
		histRuns->Write();
		CutFlowUnW->Write();
		CutFlow->Write();
		file->Write();
		file->Close();
		delete file;

	}



