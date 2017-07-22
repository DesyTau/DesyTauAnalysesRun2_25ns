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
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>
#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "HTT-utilities/QCDModelingEMu/interface/QCDModelForEMu.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"

int main(int argc, char * argv[]) {

  // **** configuration
  Config cfg(argv[1]);
  string Channel="muel";

  // kinematic cuts on electrons
  bool fillplots= false;
  const bool isData = cfg.get<bool>("IsData");

//////////systematics
  /*
  const bool ApplyTauEnergyScaleUnc     = cfg.get<bool>("ApplyTauEnergyScaleUnc");
  const double TauEnergyScaleUnc   = cfg.get<double>("TauEnergyScaleUnc");
  const bool ApplyMuEnergyScaleUnc     = cfg.get<bool>("ApplyMuEnergyScaleUnc");
  const double MuEnergyScaleUnc   = cfg.get<double>("MuEnergyScaleUnc");
  const bool ApplyElEnergyScaleUnc     = cfg.get<bool>("ApplyElEnergyScaleUnc");
  const double ElEnergyScaleUncBarrel   = cfg.get<double>("ElEnergyScaleUncBarrel");
  const double ElEnergyScaleUncEndcaps   = cfg.get<double>("ElEnergyScaleUncEndcaps");
  const bool ApplyJetEnergyCorrectionUnc     = cfg.get<bool>("ApplyJetEnergyCorrectionUnc");
  const bool ApplyJetEnergyCorrectionUncSignPositive     = cfg.get<bool>("ApplyJetEnergyCorrectionUncSignPositive");
  const bool ApplyElectronCorrectionUncSignPositive     = cfg.get<bool>("ApplyElectronCorrectionUncSignPositive");
  const bool ApplyMuonCorrectionUncSignPositive     = cfg.get<bool>("ApplyMuonCorrectionUncSignPositive");
  const bool ApplyTauCorrectionUncSignPositive     = cfg.get<bool>("ApplyTauCorrectionUncSignPositive");
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

  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCutmuel");
  const double ptMuonHighCut  = cfg.get<double>("ptMuonHighCutmuel");
  const double etaMuonCut     = cfg.get<double>("etaMuonCutmuel");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonHighCut = cfg.get<double>("isoMuonHighCutmuel");

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCutmuel");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCutmuel");
  const float etaElectronCut     = cfg.get<float>("etaElectronCutmuel");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCutmuel");
  const float dzElectronCut      = cfg.get<float>("dzElectronCutmuel");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCutmuel");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCutmuel");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronIdmuel");
  const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
  const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");
  const string mu8Ele23DZ = cfg.get<string>("Mu8Ele23DZ");
  const string mu23Ele12DZ = cfg.get<string>("Mu23Ele12DZ");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");
  const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
  const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");


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



  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");
  // topSingleMuonTriggerFile
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCutmuel");
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");



  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
 

  RecoilCorrector recoilMetCorrector("DesyTauAnalyses/NTupleMaker/data/TypeI-PFMet_Run2016BtoH.root");

  MEtSys metSys("HTT-utilities/RecoilCorrections/data/MEtSys.root");


  //const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");

////////// Triggers 
  TString LowPtLegElectron(lowPtLegElectron);
  TString HighPtLegElectron(highPtLegElectron);
  TString Mu23Ele12DZ(mu23Ele12DZ);
  TString Mu8Ele23DZ(mu8Ele23DZ);
  
  TString LowPtLegMuon(lowPtLegMuon);
  TString HighPtLegMuon(highPtLegMuon);

  
  const string Muon8TriggerFile0p15 = cfg.get<string>("Muon8TriggerEff0p15");
  const string Muon23TriggerFile0p15 = cfg.get<string>("Muon23TriggerEff0p15");
  const string Electron12TriggerFile0p1 = cfg.get<string>("Electron12TriggerEff0p1");
  const string Electron23TriggerFile0p1 = cfg.get<string>("Electron23TriggerEff0p1");


  const string MuonidIsoEffFile = cfg.get<string>("MuonidIsoEffFile");
  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEffFile");

  /////////////////////////////////////////////






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
  char ff[100];
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



if (string::npos != datasetName.find("SMS-") || string::npos != datasetName.find("stau") || string::npos != datasetName.find("C1") || string::npos != rootFileName.find("Chi")) SUSY = true;


  // qcd weight (dzeta cut)
  QCDModelForEMu qcdWeight("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu.root");
  // qcd weight DZeta cut
QCDModelForEMu qcdWeightNoDzeta("HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_nodzeta.root");



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



	cout<<" Initializing iD SF files....."<<endl;

  ScaleFactor * SF_muon80p15 = new ScaleFactor();
  ScaleFactor * SF_muon230p15 = new ScaleFactor();

  SF_muon80p15->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile0p15));
  SF_muon230p15->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon23TriggerFile0p15));

  ScaleFactor * SF_electron120p1 = new ScaleFactor();
  ScaleFactor * SF_electron230p1 = new ScaleFactor();

  SF_electron120p1->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron12TriggerFile0p1));
  SF_electron230p1->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron23TriggerFile0p1));



  cout<<" Initializing iD SF files....."<<endl;
    ScaleFactor * SF_muonIdIso0p15 = new ScaleFactor();
    SF_muonIdIso0p15->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));
  
    ScaleFactor * SF_electronIdIso0p1 = new ScaleFactor();
SF_electronIdIso0p1->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));



  cout<<" ended initialization here "<<endl;



  TFile * file;
  if (isData) file = new TFile(SaveDir+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"update");
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
if (!WithInit) cout << "Without initroottree"<<endl;


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

	  if (string::npos != datasetName.find("C1N2")) {
			  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("C1C1") || string::npos != datasetName.find("Chi")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("stau")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
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


      Float_t weight = 1;
      Float_t puweight = 1;
      Float_t topptweight = 1;
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
      LSF_weight_mu = 1.;
      LSF_weight_el = 1.;
      LSF_weight_1 = 1.;
      LSF_weight_2 = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight_1 = 1.;
      trig_weight_2 = 1.;
      trig_weight = 1.;

//needed for Recoil
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

	  if (string::npos != datasetName.find("C1N2")) {
			  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("C1C1") || string::npos != datasetName.find("Chi")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("stau")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
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


	//cout<<" BCDEF  "<<RunBCDEF<<" GH "<<RunGH<<endl;

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


      unsigned int nLowPtLegElectron = 0;
      bool isLowPtLegElectron = false;
      
      unsigned int nHighPtLegElectron = 0;
      bool isHighPtLegElectron = false;

      unsigned int nLowPtLegMuon = 0;
      bool isLowPtLegMuon = false;
      
      unsigned int nHighPtLegMuon = 0;
      bool isHighPtLegMuon = false;

      unsigned int nMu8Ele23DZ = 0;
      bool isMu8Ele23DZ = false;

      unsigned int nMu23Ele12DZ = 0;
      bool isMu23Ele12DZ = false;

      if (!SUSY){
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

	if (HLTFilter==Mu23Ele12DZ) {
	  nMu23Ele12DZ = i;
	  isMu23Ele12DZ = true;
	}
	if (HLTFilter==Mu8Ele23DZ) {
	  nMu8Ele23DZ = i;
	  isMu8Ele23DZ = true;
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
      }

      if (RunNo > 278820-1 && !isMu23Ele12DZ) {
	std::cout << "HLT filter " << Mu23Ele12DZ << " not found" << std::endl;
	exit(-1);
      }
      }


	

      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp80_general_Spring16_v1[ie];
	if (!electronMvaId&&applyElectronId) continue;
	if (applyElectronId && !analysisTree.electron_pass_conversion[ie]) continue;
	if (applyElectronId && analysisTree.electron_nmissinginnerhits[ie]>1) continue;
	if (fabs(analysisTree.electron_charge[ie]) !=1) continue;
	electrons.push_back(ie);
      }
      if (electrons.size()==0) continue;

      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (isData && analysisTree.muon_isDuplicate[im]) continue;
	if (isData && analysisTree.muon_isBad[im]) continue;
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (!analysisTree.muon_isMedium[im]) continue;
	if (fabs(analysisTree.muon_charge[im]) != 1) continue;
	muons.push_back(im);
      }


      if (muons.size()==0) continue;

      int el_index=-1;
      int mu_index=-1;
      int tau_index=-1;

      float isoMuMin = 1e+10;
      float isoEleMin = 1e+10;

      bool isMuon23matched = false;
      bool isMuon8matched  = false;
      bool isElectron23matched = false;
      bool isElectron12matched = false;
      for (unsigned int im=0; im<muons.size(); ++im) {
	bool isMu23 = false;
	bool isMu8 = false;
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
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];

	if (!SUSY){
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig<deltaRTrigMatch) {
	    if (analysisTree.trigobject_filters[iT][nHighPtLegMuon]&&
		analysisTree.muon_pt[mIndex]>ptMuonHighCut) { // Mu23 Leg
	      isMu23 = true;
	    }
	    if (analysisTree.trigobject_filters[iT][nLowPtLegMuon]&&
		analysisTree.muon_pt[mIndex]>ptMuonLowCut) { // Mu8 Leg
	      isMu8 = true;
	    }
	  }
	}
      }//isData

	  if (SUSY && analysisTree.muon_pt[mIndex]>ptMuonHighCut) isMu23 = true;
	  if (SUSY && analysisTree.muon_pt[mIndex]>ptMuonLowCut) isMu8 = true;

	  bool trigMuMatch = (isMu23) || (isMu8);

	if (applyTriggerMatch && !isMu23 && !isMu8) continue;
	  bool isEle23 = false;
	  bool isEle12 = false;


	for (unsigned int ie=0; ie<electrons.size(); ++ie) {

	  unsigned int eIndex = electrons.at(ie);

	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;


	if (!SUSY){
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      if (analysisTree.trigobject_filters[iT][nHighPtLegElectron]&&
		  analysisTree.electron_pt[eIndex]>ptElectronHighCut) { // Ele23 Leg
		isEle23 = true;
	      }
	      if (analysisTree.trigobject_filters[iT][nLowPtLegElectron]&&
		  analysisTree.electron_pt[eIndex]>ptElectronLowCut) { // Ele12 Leg
		isEle12 = true;
	      }
	    }
	  }
	}//SUSY
	  
	  if (SUSY && analysisTree.electron_pt[eIndex]>ptElectronHighCut) isEle23 = true;
	  if (SUSY && analysisTree.electron_pt[eIndex]>ptElectronLowCut) isEle12 = true;
	  bool trigMatchA = (isMu23 && isEle12) ;
	  bool trigMatchB = (isMu8 && isEle23);

	  if (applyTriggerMatch && !isEle23 && !isEle12) continue;

	  if (!trigMatchA && !trigMatchB) continue;

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

	  if ((int)mIndex!=(int)mu_index) {
	    if (relIsoMu==isoMuMin) {
	      if (analysisTree.muon_pt[mIndex]>analysisTree.muon_pt[mu_index]) {
		isoMuMin  = relIsoMu;
		mu_index = (int)mIndex;
		isoEleMin = relIsoEle;
		el_index = (int)eIndex;
		isMuon23matched = isMu23;
		isMuon8matched = isMu8;
		isElectron23matched = isEle23;
		isElectron12matched = isEle12;
	      }
	    }
	    else if (relIsoMu<isoMuMin) {
	      isoMuMin  = relIsoMu;
	      mu_index = (int)mIndex;
	      isoEleMin = relIsoEle;
	      el_index = (int)eIndex;
	      isMuon23matched = isMu23;
	      isMuon8matched = isMu8;
	      isElectron23matched = isEle23;
	      isElectron12matched = isEle12;
	    }
	  }
	  else {
	    if (relIsoEle==isoEleMin) {
	      if (analysisTree.electron_pt[eIndex]>analysisTree.electron_pt[el_index]) {
		isoEleMin = relIsoEle;
		el_index = (int)eIndex;
		isElectron23matched = isEle23;
		isElectron12matched = isEle12;
	      }
	    }
	    else if (relIsoEle<isoEleMin) {
	      isoEleMin = relIsoEle;
	      el_index = (int)eIndex;
	      isElectron23matched = isEle23;
	      isElectron12matched = isEle12;
	    }
	  }
	  
	}
      }

       //     cout << "mIndex = " << mu_index << "   eIndex = " << el_index << std::endl;

      if ((int)el_index<0) continue;
      if ((int)mu_index<0) continue;



      el_relIso[0]=isoEleMin;
      mu_relIso[0]=isoMuMin;

      if (isoEleMin>0.3 || isoMuMin>0.3) continue;
      
      double q = analysisTree.electron_charge[el_index] * analysisTree.muon_charge[mu_index];
      event_sign  = q;


  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;

  event_secondLeptonVeto = false;
  event_thirdLeptonVeto = false;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (int(ie)==(int)el_index) continue;
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie];
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
	if (int(im)==(int)mu_index) continue;
	if (isData && analysisTree.muon_isDuplicate[im]) continue;
	if (isData && analysisTree.muon_isBad[im]) continue;
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
	  neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r04_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r04_sumPUPt[im];
	}
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }


      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;


      // applying inclusive selection

   	event_secondLeptonVeto = false;

	
	if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;

//	if (event_thirdLeptonVeto) continue;


///////////////Trigger weight 
//
//
     float  pt_2 = analysisTree.muon_pt[mu_index];
     float  eta_2 = analysisTree.muon_eta[mu_index];
     float  phi_2 = analysisTree.muon_phi[mu_index];
      
     float pt_1 = analysisTree.electron_pt[el_index];
     float eta_1 = analysisTree.electron_eta[el_index];
     float phi_1 = analysisTree.electron_phi[el_index];

     float Ele23EffData0p1 = 1.;
     float Ele12EffData0p1 = 1.;
     float Mu8EffData0p15 = 1.;
     float Mu23EffData0p15 = 1.;
     float Ele23EffMC0p1 = 1.;
     float Ele12EffMC0p1 = 1.;
     float Mu8EffMC0p15 = 1.;
     float Mu23EffMC0p15 = 1.;

	trig_weight_1 = 0;
	trig_weight_2 = 0;
     float LumiA, LumiB, Lumi;
      Lumi = 36590.; LumiA = 20233.; LumiB = 16357.;

      if (!isData) {
     //// Ele_0p1 Muon_0p15

      Ele12EffData0p1 = (float)SF_electron120p1->get_EfficiencyData(double(pt_1),double(eta_1));
      Ele23EffData0p1 = (float)SF_electron230p1->get_EfficiencyData(double(pt_1),double(eta_1));
      Ele12EffMC0p1 = (float)SF_electron120p1->get_EfficiencyMC(double(pt_1),double(eta_1));
      Ele23EffMC0p1 = (float)SF_electron230p1->get_EfficiencyMC(double(pt_1),double(eta_1));

      Mu8EffData0p15 = (float)SF_muon80p15->get_EfficiencyData(double(pt_2),double(eta_2));
      Mu23EffData0p15 = (float)SF_muon230p15->get_EfficiencyData(double(pt_2),double(eta_2));
      Mu8EffMC0p15 = (float)SF_muon80p15->get_EfficiencyMC(double(pt_2),double(eta_2));
      Mu23EffMC0p15 = (float)SF_muon230p15->get_EfficiencyMC(double(pt_2),double(eta_2));

	float termA =0;   termA = (Mu23EffData0p15*Ele12EffData0p1) ; 
	float termB =0;   termB = (Mu8EffData0p15*Ele23EffData0p1); 
	float termC =0;   termC = (Mu23EffData0p15*Ele23EffData0p1);

	float termAMC =0;   termAMC =  (Mu23EffMC0p15*Ele12EffMC0p1) ; 
	float termBMC =0;   termBMC = (Mu8EffMC0p15*Ele23EffMC0p1) ; 
	float termCMC =0;   termCMC = (Mu23EffMC0p15*Ele23EffMC0p1);

	
	 float trigWeight = 0;//
	 if (termAMC + termBMC - termCMC >0.001) trigWeight = (termA + termB - termC)/ (termAMC + termBMC - termCMC); 

	 //if ( (termA==0 && termB==0) || (termAMC==0 && termBMC==0)) trigWeight =0;

	// cout<<"  "<<trigWeight<<"  "<<termA<<"  "<<termAMC<<"  "<<termB<<"  "<<termBMC<<"  "<<termC<<"  "<<termCMC<<endl;

         if (SUSY) trigWeight = Mu23EffData0p15*Ele12EffData0p1 + Mu8EffData0p15*Ele23EffData0p1 - Mu23EffData0p15*Ele23EffData0p1;

	trig_weight = (trigWeight)/Lumi + 0.953 * (trigWeight)/Lumi;
	trig_weight_1 = (trigWeight) ;
	trig_weight_2 =  0.953*(trigWeight) ;
	weight *= trig_weight;



	///LSF 

 	float  isoweight_1 = 1;
      	float  isoweight_2 = 1;

        isoweight_1 = SF_muonIdIso0p15->get_ScaleFactor(pt_2,eta_2);
        isoweight_2 = SF_electronIdIso0p1->get_ScaleFactor(pt_1,eta_1);

	LSF_weight_mu = isoweight_1 ;
	LSF_weight_el = isoweight_2 ;
	LSF_weight = LSF_weight_mu * LSF_weight_el;
      }

      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)analysisTree.muon_count;
      //cout<<" here ============================> "<<iEntry<<"  "<<mu_count<<"  "<<(int)analysisTree.muon_count<<"  "<<analysisTree.muon_count<<endl;
      for (unsigned int im=0;im<analysisTree.muon_count; ++im){
	mu_px[im]=analysisTree.muon_px[im];
	mu_py[im]=analysisTree.muon_py[im];
	mu_pz[im]=analysisTree.muon_pz[im];
	mu_eta[im]=analysisTree.muon_eta[im];
	mu_pt[im]=analysisTree.muon_pt[im];
	mu_phi[im]=analysisTree.muon_phi[im];
	mu_charge[im]=analysisTree.muon_charge[im];
	mu_dxy[im]=analysisTree.muon_dxy[im];
	mu_dz[im]=analysisTree.muon_dz[im];
	mu_dxyerr[im]=analysisTree.muon_dxyerr[im];
	mu_dzerr[im]=analysisTree.muon_dzerr[im];

        mu_neutralHadIso[im] = analysisTree.muon_r04_sumNeutralHadronEt[im];
        mu_photonIso[im] = analysisTree.muon_r04_sumPhotonEt[im];
        mu_chargedHadIso[im] = analysisTree.muon_r04_sumChargedHadronPt[im];
        mu_puIso[im] = analysisTree.muon_r04_sumPUPt[im];
 
        double neutralIso = mu_neutralHadIso[im] + mu_photonIso[im] - 0.5*mu_puIso[im];
        neutralIso = max(double(0),neutralIso);
	mu_neutralIso[im] = neutralIso;
        mu_absIsoMu[im] = mu_chargedHadIso[im] + neutralIso;
	mu_relIsoMu[im]  = mu_absIsoMu[im]/mu_pt[im] ;
   
     	}



      el_count=(int)analysisTree.electron_count;
      for (unsigned int ie=0;ie<analysisTree.electron_count; ++ie){
	el_px[ie]=analysisTree.electron_px[ie];
	el_py[ie]=analysisTree.electron_py[ie];
	el_pz[ie]=analysisTree.electron_pz[ie];
	el_eta[ie]=analysisTree.electron_eta[ie];
	el_pt[ie]=analysisTree.electron_pt[ie];
	el_phi[ie]=analysisTree.electron_phi[ie];
	el_charge[ie]=analysisTree.electron_charge[ie];
	el_dxy[ie]=analysisTree.electron_dxy[ie];
	el_dz[ie]=analysisTree.electron_dz[ie];
	el_dxyerr[ie]=analysisTree.electron_dxyerr[ie];
	el_dzerr[ie]=analysisTree.electron_dzerr[ie];

        el_neutralHadIso[ie] = analysisTree.electron_r03_sumNeutralHadronEt[ie];
        el_photonIso[ie] = analysisTree.electron_r03_sumPhotonEt[ie];
        el_chargedHadIso[ie] = analysisTree.electron_r03_sumChargedHadronPt[ie];
        el_puIso[ie] = analysisTree.electron_r03_sumPUPt[ie];
 
        double neutralIso = el_neutralHadIso[ie] + el_photonIso[ie] - 0.5*el_puIso[ie];
        neutralIso = max(double(0),neutralIso);
        el_neutralIso[ie] = neutralIso ;
        el_absIsoEl[ie] = el_chargedHadIso[ie] + el_neutralIso[ie];
	el_relIsoEl[ie]  = el_absIsoEl[ie]/el_pt[ie] ;

      }

				
      ta_count=(int)analysisTree.tau_count;
      for (unsigned int it=0;it<analysisTree.tau_count; ++it){
	ta_mass[it]=analysisTree.tau_mass[it];
	ta_px[it]=analysisTree.tau_px[it];
	ta_py[it]=analysisTree.tau_py[it];
	ta_pz[it]=analysisTree.tau_pz[it];
	ta_eta[it]=analysisTree.tau_eta[it];
	ta_pt[it]=analysisTree.tau_pt[it];
	ta_phi[it]=analysisTree.tau_phi[it];
	ta_charge[it]=analysisTree.tau_charge[it];
	ta_dxy[it]=analysisTree.tau_dxy[it];
	ta_dz[it]=analysisTree.tau_dz[it];
	//
     	ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[it];
     	ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[it];
     	ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[it];



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

      float dr_tt = deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
		    analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index]);
      qcdweight     = qcdWeightNoDzeta.getWeight(pt_1,pt_2,dr_tt);
      qcdweightup   = qcdWeightNoDzeta.getWeightUp(pt_1,pt_2,dr_tt);
      qcdweightdown = qcdWeightNoDzeta.getWeightDown(pt_1,pt_2,dr_tt);


      TLorentzVector leptonsV, muonJ, jetsLV;


      float jetEta = 2.4;
      float DRmax = 0.5;
      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;

	int counter_cleaned_jets = 0;


      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


	bool isPFJetId = false ; 
	isPFJetId =looseJetiD(analysisTree,jet);
	//isPFJetId =tightLepVetoJetiD(analysisTree,jet);

	if (!isPFJetId) continue;
	bool cleanedJet = true;


	double Dr=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Dr  < DRmax)  cleanedJet=false;

	double Dre=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if (  Dre  < DRmax )  cleanedJet=false;

	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

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
	  } //is Data

	  if (btagged && cleanedJet) bjets.push_back(jet);
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
      //njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      //nbtag_nocleaned = bjets_nocleaned.size();

      npv =  analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
	
      SusyMother = SusyMotherMassF;
      SusyLSP = SusyLSPMassF;


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


/*
     float m = sqrt(met_ex*met_ex+met_ey*met_ey);
     float m1 =sqrt(met_ex_JetEnUp*met_ex_JetEnUp+met_ey_JetEnUp*met_ey_JetEnUp);
     float m2 =sqrt(met_ex_JetEnDown*met_ex_JetEnDown+met_ey_JetEnDown*met_ey_JetEnDown);


cout<<"  ex  "<<met_ex<<"  "<<met_ey<<" exJU " <<met_ex_JetEnUp<<"  "<<met_ey_JetEnUp<<" exUncU "<<met_ex_UnclusteredEnUp<<"  "<<met_ey_UnclusteredEnUp<<endl;

cout<<" "<<m<<"  "<<m1<<"  "<<m2<<endl;

cout<<""<<endl;
*/
      float genmet_ex = analysisTree.genmet_ex;
      float genmet_ey = analysisTree.genmet_ey;

      genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
      genmetphi = TMath::ATan2(genmet_ey,genmet_ex);

      if (!isData) npartons = analysisTree.genparticles_noutgoing;

      all_weight = weight;

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
      T->Fill();
	
      selEvents++;
      /////////////////////////////////////////////////


    } // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;

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
  histPt->Write();
  CutFlowUnW->Write();
  CutFlow->Write();
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}


