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
#include "TRFIOFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
int main(int argc, char * argv[]) {



  // **** configuration
  Config cfg(argv[1]);
  string Channel="mutau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");

  //////////systematics
  /*
    const bool ApplyTauEnergyScaleUnc     = cfg.get<bool>("ApplyTauEnergyScaleUnc");
    const double TauEnergyScaleUnc   = cfg.get<double>("TauEnergyScaleUnc");
    const bool ApplyTauCorrectionUncSignPositive     = cfg.get<bool>("ApplyTauCorrectionUncSignPositive");

    const bool ApplyMuEnergyScaleUnc     = cfg.get<bool>("ApplyMuEnergyScaleUnc");
    const double MuEnergyScaleUnc   = cfg.get<double>("MuEnergyScaleUnc");
    const bool ApplyMuonCorrectionUncSignPositive     = cfg.get<bool>("ApplyMuonCorrectionUncSignPositive");

    const bool ApplyElEnergyScaleUnc     = cfg.get<bool>("ApplyElEnergyScaleUnc");
    const double ElEnergyScaleUncBarrel   = cfg.get<double>("ElEnergyScaleUncBarrel");
    const double ElEnergyScaleUncEndcaps   = cfg.get<double>("ElEnergyScaleUncEndcaps");
    const bool ApplyElectronCorrectionUncSignPositive     = cfg.get<bool>("ApplyElectronCorrectionUncSignPositive");

    const bool ApplyJetEnergyCorrectionUnc     = cfg.get<bool>("ApplyJetEnergyCorrectionUnc");
    const bool ApplyJetEnergyCorrectionUncSignPositive     = cfg.get<bool>("ApplyJetEnergyCorrectionUncSignPositive");
  */

  bool ApplyTauEnergyScaleUnc = false;
  bool ApplyTauCorrectionUncSignPositive = false;
  bool ApplyElEnergyScaleUnc = false;
  bool ApplyElectronCorrectionUncSignPositive = false;
  bool ApplyMuEnergyScaleUnc = false;
  bool ApplyMuonCorrectionUncSignPositive = false;
  bool ApplyJetEnergyCorrectionUnc = false;
  bool ApplyJetEnergyCorrectionUncSignPositive = false;

  const double TauEnergyScaleUnc   = cfg.get<double>("TauEnergyScaleUnc");
  const double MuEnergyScaleUnc   = cfg.get<double>("MuEnergyScaleUnc");
  const double ElEnergyScaleUncBarrel   = cfg.get<double>("ElEnergyScaleUncBarrel");
  const double ElEnergyScaleUncEndcaps   = cfg.get<double>("ElEnergyScaleUncEndcaps");
  //_Nominal _JetEnUp _JetEnDown  _ElEnUp _ElEnDown _MuEnUp _MuEnDown

  string BTag_ = "central";

  string Systematic=argv[5];
  if (Systematic=="1" || Systematic=="" || isData) Systematic = "Nominal";

  if (string::npos != Systematic.find("TauEnUp")){ ApplyTauEnergyScaleUnc = true; ApplyTauCorrectionUncSignPositive = true;}
  if (string::npos != Systematic.find("TauEnDown")){ ApplyTauEnergyScaleUnc = true; ApplyTauCorrectionUncSignPositive = false;}

  if (string::npos != Systematic.find("ElEnUp")){ ApplyElEnergyScaleUnc = true; ApplyElectronCorrectionUncSignPositive = true;}
  if (string::npos != Systematic.find("ElEnDown")){ ApplyElEnergyScaleUnc = true; ApplyElectronCorrectionUncSignPositive = false;}

  if (string::npos != Systematic.find("MuEnUp")){ ApplyMuEnergyScaleUnc = true; ApplyMuonCorrectionUncSignPositive = true;}
  if (string::npos != Systematic.find("MuEnDown")){ ApplyMuEnergyScaleUnc = true; ApplyMuonCorrectionUncSignPositive = false;}

  if (string::npos != Systematic.find("JetEnUp")){ ApplyJetEnergyCorrectionUnc = true; ApplyJetEnergyCorrectionUncSignPositive = true;}
  if (string::npos != Systematic.find("JetEnDown")){ ApplyJetEnergyCorrectionUnc = true; ApplyJetEnergyCorrectionUncSignPositive = false;}

  if (string::npos != Systematic.find("BTagUp")){ BTag_ = "up";}
  if (string::npos != Systematic.find("BTagDown")){ BTag_ = "down";}

  ////////////muons

  const double ptMuonCut   = cfg.get<double>("ptMuonCut");
  const double etaMuonCut     = cfg.get<double>("etaMuonCut");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonLowCut  = cfg.get<double>("isoMuonLowCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

  ////// tau
  const double ptTauCut = cfg.get<double>("ptTauCut"); 
  const double etaTauCut = cfg.get<double>("etaTauCut"); 
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double  againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
  const double  againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double  vertexz =  cfg.get<double>("vertexz");
  const double  byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // dilemuon veto 
  const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
  const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
  const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
  const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
  const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");


  //veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");



  const string dataBaseDir = cfg.get<string>("DataBaseDir");


  const string MuonidIsoEffFile = cfg.get<string>("MuonidIsoEffFile");
  const string MuontrigEffFile = cfg.get<string>("MuontrigEffFile");


  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");

  // topSingleMuonTriggerFile
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");

  const string TrigLeg  = cfg.get<string>("SingleMuonFilterName") ;
  const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");


  // vertex distributions filenames and histname


  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

  RecoilCorrector recoilMetCorrector("DesyTauAnalyses/NTupleMaker/data/TypeI-PFMet_Run2016BtoH.root");

  MEtSys metSys("HTT-utilities/RecoilCorrections/data/MEtSys.root");

  //  const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");


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
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }
  TString MainTrigger(TrigLeg);


  const double bTag   = cfg.get<double>("bTag");



  xs=1;fact=1;fact2=1;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  ifstream ifs("xsecs");
  string line;

  XSec=1.;
  xsecs=XSec;
  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;
  bool SUSY = false;
  float SusyMotherMassF;
  float SusyLSPMassF;
  char ff[200];
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  std::ifstream fileList(ff);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  int nTotalFiles = 0;

  string SaveDir=argv[3];


  if (string::npos == Systematic.find("Nominal")) {SaveDir.append("_");SaveDir.append(argv[5]);}

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<"  The systematic will be "<<Systematic<<"  and save dir will be  "<<SaveDir<<endl;


  if (string::npos != datasetName.find("SMS-") || string::npos != datasetName.find("stau") || string::npos != datasetName.find("C1") || string::npos != datasetName.find("Chi")) SUSY = true;
  if (string::npos != datasetName.find("Stop")) SUSY = true;

  cout<<" Is it a NewPhysics dataset ?  "<<SUSY<<endl;

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2016_271036-284044_80bins.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Moriond17_PU25ns_V1.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);


  string BtagCVS = "CSVv2Moriond17_2017_1_26_BtoH.csv" ;  
  if (SUSY) BtagCVS = "fastsim_csvv2_ttbar_26_1_2017.csv";

  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/"+BtagCVS);


  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,BTag_);
  if (!SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"comb");
    reader_C.load(calib,BTagEntry::FLAV_C,"comb");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");}
  if (SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"fastsim");
    reader_C.load(calib,BTagEntry::FLAV_C,"fastsim");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"fastsim");}



  float etaBTAG[2] = {0.5,2.1};
  float ptBTAG[5] = {25.,35.,50.,100.,200.};

  std::cout << std::endl;
  for (int iEta=0; iEta<2; ++iEta) {
    for (int iPt=0; iPt<5; ++iPt) {
      float sfB = reader_B.eval_auto_bounds(BTag_,BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC = reader_C.eval_auto_bounds(BTag_,BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight = reader_Light.eval_auto_bounds(BTag_,BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
      printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
    }
  }
  std::cout << std::endl;

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_ichep2016.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  float MaxBJetPt = 1000.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;

  // Z pt mass weights
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016_BtoH.root"); 
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016.root" << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get("zptmass_histo"); 
  if (histZMassPtWeights==NULL) {
    std::cout << " ZMassPT Weights histogram cannot found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << std::endl;
    exit(-1);
  }


  // Lepton Scale Factors 


  cout<<"  Initializing iD SF files....."<<endl;
  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));

  cout<<"  Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuontrigEffFile));


  cout<<" ended initialization here "<<endl;


  TFile * file;
  if (isData) file = new TFile(SaveDir+"/"+TStrName+TString("_DataDriven.root"),"update");
  else file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"update");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  SetupTree(); 

  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    //////////// for SUSY!!!
    TFile * file_ = TFile::Open(TString(filen));


    bool WithInit = true;
    if (SUSY) WithInit=false;

    if (WithInit) cout << "With initroottree"<<endl;
    if (!WithInit) cout << "Without initroottree ================="<<endl;


    TTree * _inittree = NULL;
    if (!WithInit)  _inittree = (TTree*)file_->Get(TString(ntupleName));
    if (WithInit)  _inittree = (TTree*)file_->Get(TString(initNtupleName));

    if (_inittree==NULL) continue;
    Float_t genweight;
    if (!isData)
      _inittree->SetBranchAddress("genweight",&genweight);
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
    std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; ++iEntry) {
      _inittree->GetEntry(iEntry);
      if (isData)
	histWeightsH->Fill(0.,1.);
    }


    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);

    TLorentzVector BlobA; BlobA.SetXYZT(0,0,0,0);
    TLorentzVector BlobB; BlobB.SetXYZT(0,0,0,0);
    TLorentzVector PairLV; PairLV.SetXYZT(0,0,0,0);

    if (!isData && !WithInit)
      //if (!isData)
      {    
	for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
	  {
	    analysisTree.GetEntry(iEntry);
	    histWeightsH->Fill(0.,analysisTree.genweight);

	    if (SUSY){
	      BlobA.SetXYZT(0,0,0,0);
	      BlobB.SetXYZT(0,0,0,0);
	      PairLV.SetXYZT(0,0,0,0);
	      for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

		TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
						    analysisTree.genparticles_py[igen],
						    analysisTree.genparticles_pz[igen],
						    analysisTree.genparticles_e[igen]);

		if (string::npos != datasetName.find("C1N2") || string::npos != datasetName.find("Chi")) {
		  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
		  if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;

		  //			  cout<<"   BlobA "<<BlobA.Pt()<<"  BlobB  "<<BlobB.Pt()<<endl;
		}

		if (string::npos != datasetName.find("C1C1") ) {
		  if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
		  if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
		}

		if (string::npos != datasetName.find("left") || string::npos != datasetName.find("max") ) {
		  if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
		  if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
		}
		if (string::npos != datasetName.find("right")) {
		  if (analysisTree.genparticles_pdgid[igen]==2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
		  if (analysisTree.genparticles_pdgid[igen]==-2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
		}
		if (string::npos != datasetName.find("Stop")) {
		  if (abs(analysisTree.genparticles_pdgid[igen])==1000006 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
		  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
		}

	      }

	      if (BlobA.M()>0 && BlobB.M()>0) PairLV = BlobA+BlobB;


	      if (PairLV.Pt()>0) histPt->Fill(PairLV.Pt());


	    }

	  }
      }



    float genweights=1.;

    if(!isData && WithInit) 
      {
	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	genweightsTree->SetBranchAddress("genweight",&genweights);



	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) { 
	  genweightsTree->GetEntry(iEntryInit);
	  histWeightsH->Fill(0.,genweights);


	}

      }


    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  


      bool lumi=false;
      bool CutBasedTauId = false;

      ///////////////////////////////////////////////////////////////////////////// systematic study
      if (ApplyTauEnergyScaleUnc && !isData)
	{
	  double ApplyTauCorrectionUncSign=1;
	  if (!ApplyTauCorrectionUncSignPositive) ApplyTauCorrectionUncSign = -1;

	  for (unsigned int it = 0; it<analysisTree.tau_count; ++it) 
	    {

	      analysisTree.tau_pt[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
	      analysisTree.tau_px[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
	      analysisTree.tau_py[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
	      analysisTree.tau_pz[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
	      analysisTree.tau_e[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);

	      //analysisTree.pfmet_ex = analysisTree.pfmet_ex+((analysisTree.tau_px[it]/TauEnergyScaleUnc)-analysisTree.tau_px[it]);
	      //analysisTree.pfmet_ey = analysisTree.pfmet_ey+((analysisTree.tau_py[it]/TauEnergyScaleUnc)-analysisTree.tau_py[it]);
	    }

	} 

      if (ApplyJetEnergyCorrectionUnc && !isData)
	{
	  double ApplyJetEnergyCorrectionUncSign=1;
	  for (unsigned int it = 0; it<analysisTree.pfjet_count; ++it) 
	    {
	      if (!ApplyJetEnergyCorrectionUncSignPositive) ApplyJetEnergyCorrectionUncSign = -1;
	      analysisTree.pfjet_pt[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
	      analysisTree.pfjet_px[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
	      analysisTree.pfjet_py[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
	      analysisTree.pfjet_pz[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
	      analysisTree.pfjet_e[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
	    }

	} 

      if (ApplyElEnergyScaleUnc && !isData)
	{
	  double ElEnergyScaleUnc=1;
	  double ApplyElectronCorrectionUncSign=1;
	  if (!ApplyElectronCorrectionUncSignPositive) ApplyElectronCorrectionUncSign = -1;

	  for (unsigned int it = 0; it<analysisTree.electron_count; ++it) 
	    {

	      if (analysisTree.electron_eta[it] < 1.48)  ElEnergyScaleUnc = ElEnergyScaleUncBarrel;
	      if (analysisTree.electron_eta[it] > 1.48)  ElEnergyScaleUnc = ElEnergyScaleUncEndcaps;

	      analysisTree.electron_pt[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
	      analysisTree.electron_px[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
	      analysisTree.electron_py[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
	      analysisTree.electron_pz[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);

	    }

	} 

      if (ApplyMuEnergyScaleUnc && !isData)
	{
	  double ApplyMuonCorrectionUncSign=1;
	  if (!ApplyMuonCorrectionUncSignPositive) ApplyMuonCorrectionUncSign = -1;

	  for (unsigned int it = 0; it<analysisTree.muon_count; ++it) 
	    {

	      analysisTree.muon_pt[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
	      analysisTree.muon_px[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
	      analysisTree.muon_py[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
	      analysisTree.muon_pz[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
	    }

	} 

      ///////////////////////////////////////////////////////////////////////////// systematic study end

      float topPt = 0;
      float antitopPt = 0;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;
      /////////needed for Recoil
      bool isW = false;
      bool isDY = false;
      bool isZTT = false;
      bool isZMM = false;
      bool isZEE = false;
      bool isTOP = false;
      if (!isData &&  string::npos != filen.find("JetsToLNu") ) isW=true;
      if (!isData &&  string::npos != filen.find("TTWJetsToLNu") ) isW=false;
      if (!isData &&  string::npos != filen.find("JetsToLL_M") )  isDY=true;
      if (!isData &&  string::npos != filen.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") ) isTOP=true;

      float nuPx = 0;
      float nuPy = 0;
      float nuPz = 0;
      float nuPt = 0;
      float nuPhi = 0;

      float nuPx_msv = 0;
      float nuPy_msv = 0;
      float nuPz_msv = 0;
      float nuPt_msv = 0;
      float nuPhi_msv = 0;

      float lepPx = 0;
      float lepPy = 0;
      float lepPz = 0;
      float bosonPx = 0;
      float bosonPy = 0;
      float bosonPz = 0;
      float bosonPt = 0;
      float bosonEta = 0;
      float bosonMass = -1;

      bool isZfound = false;
      bool isWfound = false;
      bool isHfound = false;
      bool isGSfound = false;
      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      std::vector<TLorentzVector> tauNeutrinos; tauNeutrinos.clear();

      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector tauNeutrinosLV;  tauNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);

      if (!isData) {

	for (unsigned int igentau=0; igentau < analysisTree.gentau_count; ++igentau) {
	  TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
					      analysisTree.gentau_py[igentau],
					      analysisTree.gentau_pz[igentau],
					      analysisTree.gentau_e[igentau]);
	  TLorentzVector tauVisLV; tauVisLV.SetXYZT(analysisTree.gentau_visible_px[igentau],
						    analysisTree.gentau_visible_py[igentau],
						    analysisTree.gentau_visible_pz[igentau],
						    analysisTree.gentau_visible_e[igentau]);
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
	    promptTausFirstCopy.push_back(tauLV);
	    promptTausLV += tauLV;
	    wDecayProductsLV += tauLV;
	  }
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {	
	    promptTausLastCopy.push_back(tauVisLV);
	    promptVisTausLV += tauVisLV;
	  }

	}
	BlobA.SetXYZT(0,0,0,0);
	BlobB.SetXYZT(0,0,0,0);
	PairLV.SetXYZT(0,0,0,0);

	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {


	  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);

	  if (SUSY){

	    if (string::npos != datasetName.find("C1N2") || string::npos != datasetName.find("Chi")) {
	      if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	      if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
	    }

	    if (string::npos != datasetName.find("C1C1")) {
	      if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	      if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
	    }


	    if (string::npos != datasetName.find("left") || string::npos != datasetName.find("max") ) {
	      if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	      if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
	    }
	    if (string::npos != datasetName.find("right")) {
	      if (analysisTree.genparticles_pdgid[igen]==2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	      if (analysisTree.genparticles_pdgid[igen]==-2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
	    }
	    if (string::npos != datasetName.find("Stop")) {
	      if (abs(analysisTree.genparticles_pdgid[igen])==1000006 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	      if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
	    }

	    if (BlobA.M()>0 && BlobB.M()>0) PairLV = BlobA+BlobB;


	    PtSystem = PairLV.Pt();
	  }

	  if (analysisTree.genparticles_pdgid[igen]==22 && analysisTree.genparticles_status[igen]==44)
	    isGSfound = true;

	  if (analysisTree.genparticles_pdgid[igen]==23) { 
	    isZfound = true;
	    zBosonLV = genLV;
	  }
	  if (analysisTree.genparticles_pdgid[igen]==25||
	      analysisTree.genparticles_pdgid[igen]==35||
	      analysisTree.genparticles_pdgid[igen]==36) { 
	    isHfound = true;
	    hBosonLV = genLV;
	  }
	  if (abs(analysisTree.genparticles_pdgid[igen])==24) { 
	    isWfound = true;
	    wBosonLV = genLV;
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptElectrons.push_back(genLV);
	      promptElectronsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptMuons.push_back(genLV);
	      promptMuonsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==12||
	      fabs(analysisTree.genparticles_pdgid[igen])==14||
	      fabs(analysisTree.genparticles_pdgid[igen])==16)  {
	    if ((analysisTree.genparticles_fromHardProcess[igen]||analysisTree.genparticles_isPrompt[igen])&&
		!analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      promptNeutrinos.push_back(genLV);
	      promptNeutrinosLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	    if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      tauNeutrinos.push_back(genLV);
	      tauNeutrinosLV += genLV;
	    }
	  }

	  /////////Matching ISR Jets

	}

	/*	if (isGSfound) {
	//	  std::cout << "gamma* found : " << std::endl;
	if (removeGammaStar) continue;
	}
	*/
	isDYTT=false;
	isDYLL=false;
	isDYLL=false;
	isDYEE=false;
	isDYMM=false;
	isFromTauTau=false;
	isFromTauMuon=false;
	isFromTauEl=false;
	isFromMuonMuon=false;
	isFromElEl=false;
	isFromOther=false;
//study the origin of mu-tau
	if (!isData){

	  if (promptTausFirstCopy.size()==2) isFromTauTau = true;
	  else if (promptTausFirstCopy.size()==1 && promptMuons.size()==1) isFromTauMuon = true;
	  else if (promptTausFirstCopy.size()==1 && promptElectrons.size()==1) isFromTauEl = true;
	  else if (promptMuons.size()==2) isFromMuonMuon = true;
	  else if (promptElectrons.size()==2) isFromElEl = true;
	  else isFromOther = true;


	}

	if (isDY) {

	  if (promptTausFirstCopy.size()==2) {
	    isZTT = true; isZMM = false; isZEE = false;
	    isDYTT=true;
	    bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz(); 
	    bosonMass = promptTausLV.M();
	    bosonEta  = promptTausLV.Eta();
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	    //mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);
	  }
	  else if (promptMuons.size()==2) {
	    isZTT = false; isZMM = true; isZEE = false;
	    isDYMM=true;
	    bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz(); 
	    bosonMass = promptMuonsLV.M(); 
	    bosonEta = promptMuonsLV.Eta();
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	    //mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
	  }
	  else {
	    isZTT = false; isZMM = false; isZEE = true;
	    isDYEE=true;
	    bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz(); 
	    bosonMass = promptElectronsLV.M();
	    bosonEta = promptElectronsLV.Eta();
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	    //if (promptElectrons.size()==2)
	    //  mtBoson_gen = mT(promptElectrons[0],promptElectrons[1]);
	  }
	  nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
	}

	else if (isW) {
	  bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
	  bosonMass = wDecayProductsLV.M();
	  if (promptTausLastCopy.size()==1) { 
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	  }
	  else if (promptMuons.size()==1) { 
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	  }
	  else { 
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	  }
	  nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
	}
	else {
	  TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
	  bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
	  TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
	  lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
	  nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
	}

	nuPt = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
	nuPhi = TMath::ATan2(nuPy,nuPx);

	bosonPt = TMath::Sqrt(bosonPx*bosonPx+bosonPy*bosonPy);

      }

      if (isDY) { // applying Z pt mass weights
	zptmassweight = 1;
	if (bosonMass>50.0) {
	  float bosonMassX = bosonMass;
	  float bosonPtX = bosonPt;
	  if (bosonMassX>1000.) bosonMassX = 1000.;
	  if (bosonPtX<1.)      bosonPtX = 1.;
	  if (bosonPtX>1000.)   bosonPtX = 1000.;
	  zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
							    histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
	}
      }



      if (!isData && ( string::npos != filen.find("TTJets")  || string::npos != filen.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") || string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8")) ) 
	{
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	    // 		cout<< "  info = " <<  int(analysisTree.genparticles_count) <<"  "<<int(analysisTree.genparticles_pdgid[igen])<<endl;

	    if (analysisTree.genparticles_pdgid[igen]==6)
	      topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    if (analysisTree.genparticles_pdgid[igen]==-6)
	      antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);


	  }

	  if (topPt>0.&&antitopPt>0.) {
	    topptweight = topPtWeight(topPt,antitopPt);
	    weight *= topptweight;
	    top_weight = topptweight;
	    // cout<<"  "<<topPt<<"  "<<antitopPt<<"  "<<topptweight<<endl;
	  }


	  histTopPt->Fill(0.,topptweight*analysisTree.genweight);
	  histTopPtSq->Fill(0.,topptweight*topptweight*analysisTree.genweight);

	}
      if (!isData ) {
	weight *= analysisTree.genweight;
	gen_weight *=analysisTree.genweight;
	// std::cout <<"analysisTree.genweight "<< float(analysisTree.genweight) << std::endl;

	double cLower, cUpper;
	vector <double> ScalesV; ScalesV.clear();
	//	ScalesV.push_back(wScale0);
	ScalesV.push_back(analysisTree.weightScale1);
	ScalesV.push_back(analysisTree.weightScale2);
	ScalesV.push_back(analysisTree.weightScale3);
	ScalesV.push_back(analysisTree.weightScale4);
	ScalesV.push_back(analysisTree.weightScale5);
	ScalesV.push_back(analysisTree.weightScale6);
	ScalesV.push_back(analysisTree.weightScale7);
	ScalesV.push_back(analysisTree.weightScale8);

	cLower = *min_element(ScalesV.begin(), ScalesV.end());
	cUpper = *max_element(ScalesV.begin(), ScalesV.end());
	histWeightsScalesUp->Fill(0.,analysisTree.genweight*cUpper);
	histWeightsScalesDown->Fill(0.,analysisTree.genweight*cLower);

	histWeightsPDFUp->Fill(0.,analysisTree.genweight*analysisTree.weightPDFup);
	histWeightsPDFDown->Fill(0.,analysisTree.genweight*analysisTree.weightPDFdown);

	lumi=true;
	wScale0 = analysisTree.weightScale0;
	wScale1 = analysisTree.weightScale1;
	wScale2 = analysisTree.weightScale2;
	wScale3 = analysisTree.weightScale3;
	wScale4 = analysisTree.weightScale4;
	wScale5 = analysisTree.weightScale5;
	wScale6 = analysisTree.weightScale6;
	wScale7 = analysisTree.weightScale7;
	wScale8 = analysisTree.weightScale8;
	wPDFUp = analysisTree.weightPDFup;
	wPDFDown = analysisTree.weightPDFdown;

      }


      if (isData)  {
	XSec = 1.;
	histRuns->Fill(analysisTree.event_run);
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

      bool Run2016A, Run2016B, Run2016C, Run2016D, Run2016E, Run2016F, Run2016G,Run2016H;
      bool RunBCDEF = false;
      bool RunGH = false;
      Run2016A=false;
      Run2016B=false;
      Run2016C=false;
      Run2016D=false;
      Run2016E=false;
      Run2016F=false;
      Run2016G=false;
      Run2016H=false;


      int RunNo = analysisTree.event_run;
      if (isData){
	if (RunNo >  271036-1 &&  RunNo < 271658+1 ) Run2016A = true;
	if (RunNo >  272007-1 &&  RunNo < 275376+1 ) Run2016B = true;
	if (RunNo >  275657-1 &&  RunNo < 276283+1 ) Run2016C = true;
	if (RunNo >  276315-1 &&  RunNo < 276811+1 ) Run2016D = true;
	if (RunNo >  276831-1 &&  RunNo < 277420+1 ) Run2016E = true;
	if (RunNo >  277772-1 &&  RunNo < 278808+1 ) Run2016F = true;
	if (RunNo >  278820-1 &&  RunNo < 280385+1 ) Run2016G = true;
	if (RunNo >  280919-1 &&  RunNo < 284044+1 ) Run2016H = true;
	//cout<<Run2016A<<"  "<<Run2016B<<"  "<<Run2016E<<endl;

	if (Run2016B || Run2016C || Run2016D || Run2016E || Run2016F) RunBCDEF = true;
	if (Run2016G || Run2016H) RunGH = true;
      }

      std::vector<TString> metFlags; metFlags.clear();
      //////////////MET filters flag

      metFlags.push_back("Flag_HBHENoiseFilter");
      metFlags.push_back("Flag_HBHENoiseIsoFilter");
      metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
      metFlags.push_back("Flag_goodVertices");
      metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
      // metFlags.push_back("Flag_METFilters");
      metFlags.push_back("Flag_eeBadScFilter");
      metFlags.push_back("Flag_BadChargedCandidateFilter");
      metFlags.push_back("Flag_BadPFMuonFilter");
      metFlags.push_back("Flag_muonBadTrackFilter");
      metFlags.push_back("Flag_chargedHadronTrackResolutionFilter");


      bool METflag = metFiltersPasses2(analysisTree, metFlags);
      met_flag = METflag;
      //if (!METflag && isData) continue;



      if (!isData) 
	{
	  puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	  //	puweight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
	  weight *=puweight; 
	  pu_weight = puweight;
	}


      bool trigAccept = false;

      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;


      if (!SUSY){
	unsigned int nfilters = analysisTree.run_hltfilters->size();
	//  std::cout << "nfiltres = " << nfilters << std::endl;
	for (unsigned int i=0; i<nfilters; ++i) {
	  //	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	  TString HLTFilter(analysisTree.run_hltfilters->at(i));
	  if (HLTFilter==MainTrigger) {
	    nMainTrigger = i;
	    isMainTrigger = true;
	  }

	}
      }//if isData check for filters

      if (SUSY) isMainTrigger = true;

      if (!isMainTrigger) {
	std::cout << "HLT filter for Mu Trigger " << MainTrigger << " not found" << std::endl;
	return(-1);
      }



      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

	if (isData && analysisTree.muon_isDuplicate[im]) continue;
	if (isData && analysisTree.muon_isBad[im]) continue;
	if (analysisTree.muon_pt[im]<SingleMuonTriggerPtCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	if ( fabs(analysisTree.muon_charge[im]) != 1) continue;
	muons.push_back((int)im);


      }
      if (muons.size()==0) continue;

      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (analysisTree.tau_decayModeFinding[it]<decayModeFinding) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
	if ( fabs(analysisTree.tau_charge[it]) != 1 ) continue;
	taus.push_back((int)it);

      }

      if (taus.size()==0)  continue;


      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float isoMuMin  = 1e+10;
      float isoTauMin = 1.; 
      float isoTau = 1.; 
      if (CutBasedTauId) isoTauMin = 1e+10;
      if (!CutBasedTauId) isoTauMin = -10;
      float ptMu = 0;
      float ptTau = 0;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;

      bool isLegMatch = false;
      for (unsigned int im=0; im<muons.size(); ++im) {
	isLegMatch = false;
	unsigned int mIndex  = muons.at(im);
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
	float photonIsoMu = analysisTree.muon_photonIso[mIndex];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
	float puIsoMu = analysisTree.muon_puIso[mIndex];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[mIndex];
	  photonIsoMu = analysisTree.muon_r04_sumPhotonEt[mIndex];
	  chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[mIndex];
	  puIsoMu = analysisTree.muon_r04_sumPUPt[mIndex];
	}
	double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	double neutralIsoMu = max(double(0),neutralIsoMuN); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];

	if (!SUSY)
	  {	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	      if (analysisTree.trigobject_filters[iT][nMainTrigger]
		  && analysisTree.muon_pt[mIndex]>SingleMuonTriggerPtCut &&
		  analysisTree.trigobject_pt[iT]>SingleMuonTriggerPtCut) { // IsoMu Leg
		float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				      analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
		if (dRtrig<deltaRTrigMatch) 
		  isLegMatch = true;

	      }

	    }
	  }
	if (SUSY && analysisTree.muon_pt[mIndex]>SingleMuonTriggerPtCut) isLegMatch = true;

	if (!isLegMatch) continue;


	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCutmutau) continue;



	  if (CutBasedTauId){
	    isoTau= analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex];

	    if ((int)mIndex!=(int)mu_index) {
	      if (relIsoMu==isoMuMin) {
		if (analysisTree.muon_pt[mIndex]>ptMu) {
		  isoMuMin  = relIsoMu;
		  ptMu = analysisTree.muon_pt[mIndex];
		  mu_index =(int)mIndex;
		  isoTauMin = isoTau;
		  ptTau = analysisTree.tau_pt[tIndex];
		  tau_index =(int)tIndex;
		}
	      }
	      else if (relIsoMu<isoMuMin) {
		isoMuMin  = relIsoMu;
		ptMu = analysisTree.muon_pt[mIndex];
		mu_index =(int)mIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index =(int)tIndex;
	      }
	    }
	    else {
	      if (isoTau==isoTauMin) {
		if (analysisTree.tau_pt[tIndex]>ptTau) {
		  ptTau = analysisTree.tau_pt[tIndex];
		  isoTauMin = isoTau;
		  tau_index =(int)tIndex;
		}
	      }
	      else if (isoTau<isoTauMin) {
		ptTau = analysisTree.tau_pt[tIndex];
		isoTauMin = isoTau;
		tau_index =(int)tIndex;
	      }
	    }

	  }	 

	  if (!CutBasedTauId){
	    isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
	    //isoTau = analysisTree.tau_chargedIsoPtSum[tIndex];

	    if ((int)mIndex!=(int)mu_index) {
	      if (relIsoMu==isoMuMin) {
		if (analysisTree.muon_pt[mIndex]>ptMu) {
		  isoMuMin  = relIsoMu;
		  ptMu = analysisTree.muon_pt[mIndex];
		  mu_index =(int)mIndex;
		  isoTauMin = isoTau;
		  ptTau = analysisTree.tau_pt[tIndex];
		  tau_index =(int)tIndex;
		}
	      }
	      else if (relIsoMu<isoMuMin) {
		isoMuMin  = relIsoMu;
		ptMu = analysisTree.muon_pt[mIndex];
		mu_index =(int)mIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index =(int)tIndex;
	      }
	    }
	    else {
	      if (isoTau==isoTauMin) {
		if (analysisTree.tau_pt[tIndex]>ptTau) {
		  ptTau = analysisTree.tau_pt[tIndex];
		  isoTauMin = isoTau;
		  tau_index =(int)tIndex;
		}
	      }
	      else if (isoTau>isoTauMin) {
		ptTau = analysisTree.tau_pt[tIndex];
		isoTauMin = isoTau;
		tau_index =(int)tIndex;
	      }
	    }
	  }

	}
      }


      bool TauId = false;

      if ( analysisTree.tau_againstElectronVLooseMVA6[tau_index]>0.5 &&   analysisTree.tau_againstMuonTight3[tau_index]>0.5) TauId = true;

      if (!TauId) continue;

      if ((int)tau_index<0) continue;
      if ((int)mu_index<0) continue;

      bool tauPass =false;

      if (!CutBasedTauId){
	tauPass=
	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index] > 0.5;
	//analysisTree.tau_chargedIsoPtSum[tau_index] < 0.8;


	isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
	ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];
	ta_IsoFlagVTight[0]=analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[tau_index];
	ta_IsoFlagTight[0]=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];
	ta_IsoFlagMedium[0]=analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[tau_index];
	ta_IsoFlagLoose[0]=analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[tau_index];
	ta_IsoFlagVLoose[0]=analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[tau_index];
	//isoTau = analysisTree.tau_chargedIsoPtSum[tau_index];
	//ta_IsoFlag=analysisTree.tau_chargedIsoPtSum[tau_index];

      }

      if (CutBasedTauId){
	tauPass=
	  analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] > 0.5;

	isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
	ta_IsoFlag=analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index];
      }

      //	if (!tauPass) continue;

      ta_relIso[0]=isoTauMin;
      mu_relIso[0]=isoMuMin;

      double q = analysisTree.tau_charge[tau_index] * analysisTree.muon_charge[mu_index];
      event_sign  = q;

      double dRmutau = deltaR(analysisTree.tau_eta[(int)tau_index],analysisTree.tau_phi[(int)tau_index],
			      analysisTree.muon_eta[(int)mu_index],analysisTree.muon_phi[(int)mu_index]);
      if (dRmutau < 0.5) continue;



      bool          dilepton_veto=false;
      bool          extraelec_veto=false;
      bool          extramuon_veto=false;

      event_secondLeptonVeto = false;
      event_thirdLeptonVeto = false;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie];
	if (!electronMvaId&&applyVetoElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
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
	double neutralIsoEleN = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	double neutralIsoEle = max(double(0),neutralIsoEleN); 
	float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];	
	if (relIsoEle>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      // looking for extra muon's (dimuon veto)
      bool foundExtraMuon = false;
      vector<int> mu_dimuons; mu_dimuons.clear(); 
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if ((int)im==(int)mu_index) continue;

	if (isData && analysisTree.muon_isDuplicate[im]) continue;
	if (isData && analysisTree.muon_isBad[im]) continue;
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	float photonIsoMu = analysisTree.muon_photonIso[im];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	float puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r04_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r04_sumPUPt[im];
	}
	double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	double neutralIsoMu = max(double(0),neutralIsoMuN); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];

	if (analysisTree.muon_pt[im]>ptDilepMuonCut&&
	    fabs(analysisTree.muon_eta[im])<etaDilepMuonCut&&
	    analysisTree.muon_isGlobal[im]&&
	    analysisTree.muon_isTracker[im]&&
	    analysisTree.muon_isPF[im]&&
	    fabs(analysisTree.muon_dxy[im])<dxyDilepMuonCut&&
	    fabs(analysisTree.muon_dz[im])<dzDilepMuonCut&&
	    relIsoMu<isoDilepMuonCut &&
	    fabs(analysisTree.muon_charge[im]) ==1)  
	  {
	    float dRmuons = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				   analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);
	    if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[im]*analysisTree.muon_charge[mu_index]<0)) dilepton_veto = true;
	  }
	// mu_dimuons.push_back(im);

	if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

      /*
	if (mu_dimuons.size()>1) {
	for (unsigned int i1=0; i1<mu_dimuons.size()-1; ++i1) {
	1unsigned int indx1 = mu_dimuons[i1];
	for (unsigned int i2=i1+1; i2<mu_dimuons.size(); ++i2 ) {
	unsigned int indx2 = mu_dimuons[i2];
	float dRmuons = deltaR(analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1],
	analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[indx1]*analysisTree.muon_charge[indx2]<0)) dilepton_veto = true;
	}
	}
	}
      */
      //      cout << analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] << endl;
      //      cout << analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index] << endl;



      // applying inclusive selection

      event_secondLeptonVeto = dilepton_veto;
      //	if (dilepton_veto)  continue;


      //	if (extraelec_veto) continue;
      //	if (extramuon_veto) continue;

      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;




      ///////////////Trigger weight 
      double ptMu1 = (double)analysisTree.muon_pt[mu_index];
      double etaMu1 = (double)analysisTree.muon_eta[mu_index];
      float trigweight=1.;

      float EffFromData = 1.;
      float EffFromMC = 1.;


      if (!isData) {
	EffFromData = (float)SF_muonTrigger->get_EfficiencyData(double(ptMu1),double(etaMu1));
	EffFromMC   = (float)SF_muonTrigger->get_EfficiencyMC(double(ptMu1),double(etaMu1));

	if (EffFromMC>1e-6)
	  trigweight = EffFromData / EffFromMC;
	if (SUSY)  trigweight = EffFromData ;

	weight *= trigweight;
	trig_weight = trigweight;


	double IdIsoSF_mu = 1;

	IdIsoSF_mu =  SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);

	LSF_weight = IdIsoSF_mu;
	weight *= LSF_weight;
      }

      ///////////////Check if the selected tau has a gen matched - if not, apply Tau Fake Rate


      if (!isData){
	genTauMatched = false;
	genLeptonMatched = false;
	genLeptonMatchedPromptEl = false;
	genLeptonMatchedPromptMu = false;
	genLeptonMatchedPromptTau = false;
	genElMatchedToTauDecay = false;
	genMuMatchedToTauDecay = false;
	genTauMatchedToTauDecay = false;
	genElMatchedHadrDecay = false;
	genMuMatchedHadrDecay = false;
	genTauMatchedHadrDecay = false;
	genLeptonMatchedHFQ = false;
	genLeptonMatchedLFQ = false;
	genLeptonMatchedGluon =false;
	matchedTauToPromptEl = false;
	matchedTauToPromptMu = false;
	matchedTauToTauDecEl =false;
	matchedTauToTauDecMu =false;
	matchedTauToElHadronDec = false;
	matchedTauToMuHadronDec = false;
	matchedTauToTauHadronDec = false;
	matchedTauToGluon = false;
	matchedTauToHFQ = false;
	matchedTauToLFQ = false;
	genTauDecayMode1=-1;
	genTauDecayMode2=-1;
	TLorentzVector genTauV;  
	TLorentzVector genPartV;  

	bool FoundFirstMatchedTau = false;

	for (unsigned int gt = 0 ; gt < analysisTree.gentau_count; ++gt){

	  genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);

	  double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    genTauV.Eta(), genTauV.Phi());

	  if (Drr < 0.5 && ( analysisTree.gentau_isPrompt[gt] > 0.5  ) && genTauV.Pt() > 15. ) genTauMatched = true;
	  if (genTauMatched && !FoundFirstMatchedTau) {genTauDecayMode1 = analysisTree.gentau_decayMode[gt]; FoundFirstMatchedTau = true;}
	  if (genTauMatched && FoundFirstMatchedTau)   genTauDecayMode2 = analysisTree.gentau_decayMode[gt];

	}

	for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

	  if ( (abs(analysisTree.genparticles_pdgid[igen])==11 || abs(analysisTree.genparticles_pdgid[igen])==13 || abs(analysisTree.genparticles_pdgid[igen])==15 || abs(analysisTree.genparticles_pdgid[igen])<6 || abs(analysisTree.genparticles_pdgid[igen])==21)){

	    genPartV.SetXYZT(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen], analysisTree.genparticles_e[igen]);

	    double Drl=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			      genPartV.Eta(),genPartV.Phi());

	    double DrTauLepton=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
				      genPartV.Eta(),genPartV.Phi());

	    if (Drl < 0.5 && genPartV.Pt() > 8){
	      //if ( abs(analysisTree.genparticles_pdgid[igen])==13)
	      //cout<<analysisTree.genparticles_pdgid[igen]<<" isDirectHadronDecayProduct "<<analysisTree.genparticles_isDirectHadronDecayProduct[igen]<<" isTauDecay "<<analysisTree.genparticles_isTauDecayProduct[igen]<<" isDirectTau  "<<analysisTree.genparticles_isDirectTauDecayProduct[igen]<<"  isPrompt  "<<analysisTree.genparticles_isPrompt[igen]<<"  isDecayedLeptonHadron "<<analysisTree.genparticles_isDecayedLeptonHadron[igen]<<" is there gentau matched ?  "<<genTauMatched<<endl;
	      genLeptonMatched = true;

	      if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptEl = true;
	      if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptMu = true;
	      if ( abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptTau = true;

	      if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genElMatchedToTauDecay = true;
	      if (abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genMuMatchedToTauDecay = true;
	      if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genTauMatchedToTauDecay = true;

	      if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genElMatchedHadrDecay = true;
	      if (abs(analysisTree.genparticles_pdgid[igen])==13 && (analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ||  analysisTree.genparticles_isDecayedLeptonHadron[igen] >0.5 )) genMuMatchedHadrDecay = true;
	      if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genTauMatchedHadrDecay = true;

	      if ( (abs(analysisTree.genparticles_pdgid[igen])==4 || abs(analysisTree.genparticles_pdgid[igen])==5) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedHFQ = true;

	      if ( (abs(analysisTree.genparticles_pdgid[igen])==1 || abs(analysisTree.genparticles_pdgid[igen])==2 || abs(analysisTree.genparticles_pdgid[igen])==3) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedLFQ = true;
	      if ( (abs(analysisTree.genparticles_pdgid[igen])==21 ) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedGluon = true;

	    }

	    if (DrTauLepton < 0.5 && genPartV.Pt() > 8. ) 

	      {

		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToElHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToMuHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToTauHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==21 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToGluon= true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])==4  || abs(analysisTree.genparticles_pdgid[igen])==5 ) && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToHFQ= true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])>0  && abs(analysisTree.genparticles_pdgid[igen])<4 ) && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToLFQ= true;

	      }

	  }//loop for gen 11 13 15

	}


      }//!isData



      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)muons.size();

      for (unsigned int im=0;im<muons.size(); ++im){
	unsigned int mIndex = muons[im];
	mu_px[im]=analysisTree.muon_px[mIndex];
	mu_py[im]=analysisTree.muon_py[mIndex];
	mu_pz[im]=analysisTree.muon_pz[mIndex];
	mu_eta[im]=analysisTree.muon_eta[mIndex];
	mu_pt[im]=analysisTree.muon_pt[mIndex];
	mu_phi[im]=analysisTree.muon_phi[mIndex];
	mu_charge[im]=analysisTree.muon_charge[mIndex];
	mu_dxy[im]=analysisTree.muon_dxy[mIndex];
	mu_dz[im]=analysisTree.muon_dz[mIndex];
	mu_dxyerr[im]=analysisTree.muon_dxyerr[mIndex];
	mu_dzerr[im]=analysisTree.muon_dzerr[mIndex];

	mu_neutralHadIso[im] = analysisTree.muon_r04_sumNeutralHadronEt[mIndex];
	mu_photonIso[im] = analysisTree.muon_r04_sumPhotonEt[mIndex];
	mu_chargedHadIso[im] = analysisTree.muon_r04_sumChargedHadronPt[mIndex];
	mu_puIso[im] = analysisTree.muon_r04_sumPUPt[mIndex];

	double neutralIso = mu_neutralHadIso[mIndex] + mu_photonIso[mIndex] - 0.5*mu_puIso[mIndex];
	neutralIso = max(double(0),neutralIso);
	mu_neutralIso[mIndex] = neutralIso;
	mu_absIsoMu[im] = mu_chargedHadIso[mIndex] + neutralIso;
	mu_relIsoMu[im]  = mu_absIsoMu[mIndex]/mu_pt[mIndex] ;

      }



      ta_count=(int)taus.size();

      for (unsigned int it=0;it<taus.size(); ++it){
	unsigned int itt = taus[it];
	ta_px[it]=analysisTree.tau_px[itt];
	ta_py[it]=analysisTree.tau_py[itt];
	ta_pz[it]=analysisTree.tau_pz[itt];
	ta_eta[it]=analysisTree.tau_eta[itt];
	ta_pt[it]=analysisTree.tau_pt[itt];
	ta_phi[it]=analysisTree.tau_phi[itt];
	ta_charge[it]=analysisTree.tau_charge[itt];
	ta_dxy[it]=analysisTree.tau_dxy[itt];
	ta_dz[it]=analysisTree.tau_dz[itt];
	//
	ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[itt];
	ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[itt];
	ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[itt];

      }
      jet_count=(int)analysisTree.pfjet_count;
      for (unsigned int jj=0;jj<analysisTree.pfjet_count; ++jj){

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




      float jetEta = 2.4;
      float DRmax = 0.5;

      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjetsTight; bjetsTight.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;

      int indexLeadingBJet = -1;

      int counter_cleaned_jets = 0;

      for (int n=0;n<30;n++){


	jets_cleaned[n]=-1;
	bjets_cleaned[n]=-1;
	bjets_cleanedTight[n]=-1;
      }

      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	double drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);



	if (drr < 0.5 && !isData) 
	  {

	    if (analysisTree.pfjet_flavour[jet] == 21) matchedTauToGluon = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 4 || abs(analysisTree.pfjet_flavour[jet]) == 5) matchedTauToHFQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) < 4) matchedTauToLFQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 1 ) matchedTauToDownQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 2 ) matchedTauToUpQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 3 ) matchedTauToStrangeQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 4 ) matchedTauToCharmQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 5 ) matchedTauToBottomQ = true;

	  }

	double drl=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);


	if (drl < 0.5 && !isData) 
	  {

	    if (analysisTree.pfjet_flavour[jet] == 21) genLeptonMatchedGluon = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) == 4 || abs(analysisTree.pfjet_flavour[jet]) == 5) genLeptonMatchedHFQ = true;
	    if (abs(analysisTree.pfjet_flavour[jet]) < 4) genLeptonMatchedLFQ = true;

	  }

	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


	bool isPFJetId = false ; 
	isPFJetId =looseJetiD(analysisTree,jet);
	//isPFJetId =tightLepVetoJetiD(analysisTree,jet);

	if (!isPFJetId) continue;
	bool cleanedJet = true;
	bool cleanedJetTau = true;
	bool cleanedJetMu = true;

	double Dr=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Dr  < DRmax) {
	  cleanedJet=false;
	  cleanedJetMu=false;

	}

	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if ( Drr < DRmax) {

	  cleanedJet=false;
	  cleanedJetTau=false;
	}

	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

	  bool btaggedTight  (analysisTree.pfjet_btag[jet][0]  > 0.9535) ;
	  bool btagged  (analysisTree.pfjet_btag[jet][0]  > bTag) ;

	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);

	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;

	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_B.eval_auto_bounds(BTag_,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_C.eval_auto_bounds(BTag_,BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader_Light.eval_auto_bounds(BTag_,BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
	      tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
	    }

	    if (tageff<1e-5)      tageff = 1e-5;
	    if (tageff>0.99999)   tageff = 0.99999;
	    rand.SetSeed((int)((jetEta+5)*100000));
	    double rannum = rand.Rndm();

	    if (jet_scalefactor<1 && btagged) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		btagged = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !btagged) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) { 
		btagged = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }

	    if (jet_scalefactor<1 && btaggedTight) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		btaggedTight = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !btaggedTight) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) { 
		btaggedTight = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }


	  } //is Data

	  if (btagged && cleanedJet) 	   bjets.push_back(jet);
	  if (btaggedTight && cleanedJet) 	   bjetsTight.push_back(jet);

	  if (btaggedTight && cleanedJet) bjets_cleanedTight[counter_cleaned_jets]=(int)jet;
	  if (btagged && cleanedJet) bjets_cleaned[counter_cleaned_jets]=(int)jet;

	}


	if (cleanedJet){

	  //	cout<<"  will push to save now cleaned jet  "<<(int)jet<<"  for counter_cleaned_jet "<<(int)counter_cleaned_jets<<" event "<<iEntry<<endl;

	  jets.push_back((int)jet);
	  jets_cleaned[counter_cleaned_jets]=(int)jet;
	  jet_jecUn[counter_cleaned_jets] = analysisTree.pfjet_jecUncertainty[jet];
	  counter_cleaned_jets++;
	}


      }///loop in all jets

      njets = jets.size();
      jet_count = jets.size();
      nbtag = bjets.size();
      nbtagTight = bjetsTight.size();

      npv =  analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;

      // SusyMother = SusyMotherMassF;
      // SusyLSP = SusyLSPMassF;


      /////////////////// Recoil corrections

      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;

      ////while using old MC ntuples, need to use proper MET collection
      float pfmet_corr_x = 1.;
      float pfmet_corr_y = 1.;
      float met_x = 1.;
      float met_y = 1.;

      pfmet_corr_x = analysisTree.pfmetcorr_ex;
      pfmet_corr_y = analysisTree.pfmetcorr_ey;
      met_x = analysisTree.pfmetcorr_ex;
      met_y = analysisTree.pfmetcorr_ey;


      if ((isW || isDY) && !isData) {

	recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);

	float met_corr_x=1.;
	float met_corr_y=1.;		
	met_x= analysisTree.pfmetcorr_ex_JetEnUp;
	met_y= analysisTree.pfmetcorr_ey_JetEnUp;
	met_corr_x= analysisTree.pfmetcorr_ex_JetEnUp;
	met_corr_y= analysisTree.pfmetcorr_ey_JetEnUp;

	recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
	met_ex_JetEnUp_recoil = met_corr_x;
	met_ey_JetEnUp_recoil = met_corr_y;


	met_x= analysisTree.pfmetcorr_ex_JetEnDown;
	met_y= analysisTree.pfmetcorr_ey_JetEnDown;
	met_corr_x= analysisTree.pfmetcorr_ex_JetEnDown;
	met_corr_y= analysisTree.pfmetcorr_ey_JetEnDown;
	recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
	met_ex_JetEnDown_recoil = met_corr_x;
	met_ey_JetEnDown_recoil = met_corr_y;


	met_x=analysisTree.pfmetcorr_ex_UnclusteredEnUp;
	met_y=analysisTree.pfmetcorr_ey_UnclusteredEnUp;
	met_corr_x=analysisTree.pfmetcorr_ex_UnclusteredEnUp;
	met_corr_y=analysisTree.pfmetcorr_ey_UnclusteredEnUp;
	recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
	met_ex_UnclusteredEnUp_recoil = met_corr_x;
	met_ey_UnclusteredEnUp_recoil = met_corr_y;


	met_x=analysisTree.pfmetcorr_ex_UnclusteredEnDown;
	met_y=analysisTree.pfmetcorr_ey_UnclusteredEnDown;
	met_corr_x=analysisTree.pfmetcorr_ex_UnclusteredEnDown;
	met_corr_y=analysisTree.pfmetcorr_ey_UnclusteredEnDown;
	recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
	met_ex_UnclusteredEnDown_recoil = met_corr_x;
	met_ey_UnclusteredEnDown_recoil = met_corr_y;


	met_x = pfmet_corr_x;
	met_y = pfmet_corr_y;

	// MEt related systematic uncertainties
	int bkgdType = 0;
	if (isDY||isW)
	  bkgdType = MEtSys::ProcessType::BOSON;
	else if (isTOP)
	  bkgdType = MEtSys::ProcessType::TOP;
	else 
	  bkgdType = MEtSys::ProcessType::EWK; 

	float met_scaleUp_x   = met_x;
	float met_scaleUp_y   = met_y;
	float met_scaleDown_x = met_x;
	float met_scaleDown_y = met_y;
	float met_resoUp_x    = met_x;
	float met_resoUp_y    = met_y;
	float met_resoDown_x  = met_x;
	float met_resoDown_y  = met_y;

	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Up,
			   met_scaleUp_x,met_scaleUp_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Down,
			   met_scaleDown_x,met_scaleDown_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Up,
			   met_resoUp_x,met_resoUp_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Down,
			   met_resoDown_x,met_resoDown_y);


	met_scaleUp = TMath::Sqrt(met_scaleUp_x*met_scaleUp_x+
				  met_scaleUp_y*met_scaleUp_y);
	metphi_scaleUp = TMath::ATan2(met_scaleUp_y,met_scaleUp_x);

	met_scaleDown = TMath::Sqrt(met_scaleDown_x*met_scaleDown_x+
				    met_scaleDown_y*met_scaleDown_y);
	metphi_scaleDown = TMath::ATan2(met_scaleDown_y,met_scaleDown_x);

	met_resoUp = TMath::Sqrt(met_resoUp_x*met_resoUp_x+
				 met_resoUp_y*met_resoUp_y);
	metphi_resoUp = TMath::ATan2(met_resoUp_y,met_resoUp_x);

	met_resoDown = TMath::Sqrt(met_resoDown_x*met_resoDown_x+
				   met_resoDown_y*met_resoDown_y);
	metphi_resoDown = TMath::ATan2(met_resoDown_y,met_resoDown_x);

	met_ex_recoil = pfmet_corr_x;
	met_ey_recoil = pfmet_corr_y;
	//revert back to uncorrected met
	if(!isData)
	  {      met_x = analysisTree.pfmetcorr_ex;
	    met_y = analysisTree.pfmetcorr_ey;
	  }


      }//if isW, isDY !isData

      met_ex = met_x;
      met_ey = met_y;
      met_ez = analysisTree.pfmet_ez;
      met_pt = TMath::Sqrt(met_ex*met_ex + met_ey*met_ey);
      met_phi = TMath::ATan2(met_y,met_x);

      met_ex_JetEnUp = analysisTree.pfmetcorr_ex_JetEnUp;
      met_ey_JetEnUp = analysisTree.pfmetcorr_ey_JetEnUp;

      met_ex_JetEnDown = analysisTree.pfmetcorr_ex_JetEnDown;
      met_ey_JetEnDown = analysisTree.pfmetcorr_ey_JetEnDown;

      met_ex_UnclusteredEnUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
      met_ey_UnclusteredEnUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;

      met_ex_UnclusteredEnDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
      met_ey_UnclusteredEnDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;


      float genmet_ex = analysisTree.genmet_ex;
      float genmet_ey = analysisTree.genmet_ey;

      genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
      genmetphi = TMath::ATan2(genmet_ey,genmet_ex);

      if (!isData) npartons = analysisTree.genparticles_noutgoing;

      all_weight = weight;
      event_run = analysisTree.event_run;
      event_lumi = analysisTree.event_luminosityblock;
      NuPx = nuPx;
      NuPy = nuPy;
      NuPz = nuPz;
      NuPt = nuPt;
      NuPhi = nuPhi;
      genmet_Ex = analysisTree.genmet_ex;
      genmet_Ey = analysisTree.genmet_ey;

      genHT = analysisTree.genparticles_lheHt;



      T->Fill();

      selEvents++;
      /////////////////////////////////////////////////



    } // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }


  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;


  std::cout << std::endl;
  int allEvents = (int)inputEventsH->GetEntries();
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;


  file->cd(Channel.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histWeightsScalesUp->Write();
  histWeightsScalesDown->Write();
  histWeightsPDFUp->Write();
  histWeightsPDFDown->Write();
  histTopPt->Write();
  histTopPtSq->Write();
  histRuns->Write();
  CutFlowUnW->Write();
  histPt->Write();
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}
