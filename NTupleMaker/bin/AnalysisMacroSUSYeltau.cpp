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
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/ScaleFactor.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  // **** configuration
  Config cfg(argv[1]);
  string Channel="eltau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");


  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

  // kinematics electrons
  const float  ptElectronCut       = cfg.get<float>("ptElectronCuteltau");
  const double etaElectronCut     = cfg.get<double>("etaElectronCuteltau");
  const double dxyElectronCut     = cfg.get<double>("dxyElectronCuteltau");
  const double dzElectronCut      = cfg.get<double>("dzElectronCuteltau");
  const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCuteltau");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

 // tau kinematics
  const float ptTauCut   = cfg.get<float>("ptTauCut");
  const float etaTauCut     = cfg.get<float>("etaTauCut");
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double  againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
  const double   againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double   againstMuonLoose3  = cfg.get<double>("againstMuonLoose3");
  const double   vertexz =  cfg.get<double>("vertexz");
  const double   byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");
 
 

  // veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // dilepton veto
  const float ptDilepElectronCut = cfg.get<float>("ptDilepElectronCuteltau");
  const float etaDilepElectronCut = cfg.get<float>("etaDilepElectronCuteltau");
  const float dxyDilepElectronCut = cfg.get<float>("dxyDilepElectronCuteltau");
  const float dzDilepElectronCut = cfg.get<float>("dzDilepElectronCuteltau");
  const float isoDilepElectronCut = cfg.get<float>("isoDilepElectronCuteltau");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCuteltau");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");


  const string dataBaseDir = cfg.get<string>("DataBaseDir");



  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEffFile");
  const string SingleElectronTriggerFile = cfg.get<string>("ElectrontrigEffFile");


  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");

  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");




  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");


  // topological cuts
  const double dRleptonsCuteltau   = cfg.get<double>("dRleptonsCuteltau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  string TrigLeg = cfg.get<string>("SingleElectronFilterName");
  const float SingleElectronTriggerPtCut = cfg.get<float>("SingleElectronTriggerPtCut");


  //if (isData) TrigLeg  = cfg.get<string>("El23LegData");
 // const float singleElectronTriggerEtaCut = cfg.get<float>("SingleElectronTriggerEtaCuteltau");


 // const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
  //const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

  // vertex distributions filenames and histname

  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

  //const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");

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

  TString MainTrigger(TrigLeg);



  const double Lumi   = cfg.get<double>("Lumi");
  const double bTag   = cfg.get<double>("bTag");


  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("No cut after PU");
  CutList.push_back("el-tau");
  CutList.push_back("2nd lepV");
  CutList.push_back("3rd lepV");
  CutList.push_back("Trigger");
  CutList.push_back("Lepton SF");
  CutList.push_back("TauFakeRate");
  CutList.push_back("topPtRwgt");
  CutList.push_back("Cleaned jets");
  CutList.push_back("b-Veto");


  int CutNumb = int(CutList.size());
  xs=1;fact=1;fact2=1;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  ifstream ifs("xsecs");
  string line;
 
 
/*
  while(std::getline(ifs, line)) // read one line from ifs
    {

      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt;
      iss >> dt >> xs >> fact >> fact2;
      //datasetName = dt.c_str();
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      if (  dt == argv[2]) {
	XSec= xs*fact*fact2;
	cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<" number of expected events for Lumi "<<Lumi <<" /pb  = " <<XSec*Lumi<<endl;
      }
      
	if ( argv[2] == st1) {ChiMass=100;mIntermediate=200;}
	else if (argv[2] == st2) {ChiMass=200;mIntermediate=500;}
      
      if (isData) XSec=1.;
      ChiMass=0.0;
    }

  if (XSec<0&& !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;XSec = 1;}*/


	XSec=1.;
	ChiMass=0.0;	
  xsecs=XSec;
  std::vector<unsigned int> allRuns; allRuns.clear();

  cout<<" ChiMass is "<<ChiMass<<"  "<<mIntermediate<<endl;
  bool doThirdLeptVeto=true;
  bool doDiElVeto=true;


// PU reweighting
  PileUp * PUofficial = new PileUp();
   //TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_2016_Cert_Cert_271036-276811_NoL1T_xsec63mb.root","read");

  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_Cert_271036-277148_13TeV_PromptReco_Collisions16_xsec69p2mb.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring16_PU25ns_V1.root", "read");


  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);


  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv");
  BTagCalibrationReader reader_BC(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central");
  reader_BC.load(calib,BTagEntry::FLAV_B,"mujets");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");

  //  std::cout << "SF_light (eta=0.6,pt=20.1) : " << reader_Light.eval(BTagEntry::FLAV_UDSG, 0.5, 20.1) << std::endl;
  //  std::cout << "SF_light (eta=2.1,pt=20.1) : " << reader_Light.eval(BTagEntry::FLAV_UDSG, 2.1, 20.1) << std::endl;
  //  std::cout << "SF_bc    (eta=0.6,pt=30.1) : " << reader_BC.eval(BTagEntry::FLAV_B, 0.5, 30.1) << std::endl;
  //  std::cout << "SF_bc    (eta=2.1,pt=30.1) : " << reader_BC.eval(BTagEntry::FLAV_B, 2.1, 30.1) << std::endl;

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;
	
  float MaxBJetPt = 670.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;




  // Lepton Scale Factors 

  TH1D * ElSF_IdIso_El1H = new TH1D("ElIdIsoSF_El1H", "ElIdIsoSF_El1", 100, 0.5,1.5);

	cout<<" Initializing iD SF files....."<<endl;

  ScaleFactor * SF_elIdIso; 
  if (applyLeptonSF) {
    SF_elIdIso = new ScaleFactor();
    SF_elIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
  }

	cout<<" Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_electronTrigger = new ScaleFactor();
  SF_electronTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(SingleElectronTriggerFile));

/*  cout<<" Will try to initialize the TFR now.... "<<endl;
  ScaleFactor * SF_TFR; 
  bool applyTFR = true;
  if (applyTFR) {
    SF_TFR = new ScaleFactor();
    SF_TFR->init_ScaleFactorb(TString(cmsswBase)+"/src/"+TString(TauFakeRateFile),applyTFR);
  }
*/

 char ff[100];
 sprintf(ff,"%s/%s",argv[3],argv[2]);

  double Weight=0;
  int nTotalFiles = 0;
  int iCut=0;
  double CFCounter[CutNumb];
  double statUnc[CutNumb];
  int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; i++){
    CFCounter[i] = 0;
    CFCounter_[i] = 0;
    iCFCounter[i] = 0;
    statUnc[i] =0;
  }
  // file name and tree name
  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString era=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  //regionName = Region.c_str()+"_"+Sign.c_str();
  std::cout <<" The filename will be "<<TStrName <<std::endl;  
  // output fileName with histograms
std::string st1,st2;
bool SUSY = false;
float SusyMotherMassF;
float SusyLSPMassF;
if (string::npos != rootFileName.find("SMS-") || string::npos != rootFileName.find("stau") || string::npos != rootFileName.find("C1"))
	{
	//st1 =  rootFileName.substr(4,3);
	SusyMotherMassF = stof(argv[5]);
	//st1=string(argv[5]);
	//st2 =  rootFileName.substr(11);
	//st2=string(argv[6]);
	SusyLSPMassF = stof(argv[6]);
	SUSY = true;
	  std::cout <<" SUSY "<< " SusyMotherMassF= "<<SusyMotherMassF <<" SusyLSPMassF= "<<SusyLSPMassF <<std::endl;  
	}
/*
if (string::npos != rootFileName.find("SMS-TChiStauStau"))
	{
	st1 =  rootFileName.substr(5,3);
	SusyMotherMassF = stof(st1);
	st2 =  rootFileName.substr(12);
	SusyLSPMassF = stof(st2);
	SUSY = true;
	  std::cout <<" SUSY "<< " SusyMotherMassF= "<<SusyMotherMassF <<" SusyLSPMassF= "<<SusyLSPMassF <<std::endl;  
	}
*/



  TFile * file;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
/*  if (SUSY)
	{  
	TString TStrNameS(rootFileName+invMuStr+invTauStr+invMETStr+"_"+Region+"_"+Sign+"_"+st1+"_"+st2);
	file = new TFile(era+"/"+TStrNameS+TString(".root"),"update");

	}*/
  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->  
  while (fileList0 >> dummy) nTotalFiles++;
  //string treename = rootFileName+"_tree.root";
  
  SetupTree(); 
  SetupHists(CutNumb); 
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;


 for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
//////////// for SUSY!!!
//if (SUSY){
//if (iF+1 != nTotalFiles) continue;}
    TFile * file_ = TFile::Open(TString(filen));
/*
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
    int NE = (int)histoInputEvents->GetEntries();
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
    std::cout << "      number of input events         = " << NE << std::endl;
*/

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
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
      _inittree->GetEntry(iEntry);
      if (isData)
	histWeightsH->Fill(0.,1.);
      //else
      //histWeightsH->Fill(0.,genweight);
    }




    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);

	if (!isData && !WithInit)
	//if (!isData)
		{    
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
			{
			analysisTree.GetEntry(iEntry);
		/*	if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
			histWeightsH->Fill(0.,analysisTree.genweight);
			}
		}
  	float genweights=1.;

	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	genweightsTree->SetBranchAddress("genweight",&genweights);

    if(!isData && WithInit) 
      {

	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) { 
	  genweightsTree->GetEntry(iEntryInit);
/*			if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
	  histWeightsH->Fill(0.,genweights);
	}
    
      }

    float topPt = -1;
    float antitopPt = -1;






    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
/*	if (SUSY)////FOR SUSY!!!!!!!!
		{
		if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
		&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
		}*/
      nEvents++;
      iCut = 0;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 
      /*
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  
*/

      //isData= false;
      bool lumi=false;
      bool CutBasedTauId = false;

      float topPt = 0;
      float antitopPt = 0;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;

	  

      if (!isData && ( string::npos != filen.find("TTJets")  || string::npos != filen.find("TTPowHeg") || string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8")) ) 
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


      histTopPt->Fill(0.,topptweight);

	}
	  if (!isData ) {
	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;
	   // std::cout <<"analysisTree.genweight "<< float(analysisTree.genweight) << std::endl;
	  lumi=true;
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

	std::vector<TString> metFlags; metFlags.clear();
     //////////////MET filters flag
//      if (isData){

	 
	 metFlags.push_back("Flag_HBHENoiseFilter");
	 metFlags.push_back("Flag_HBHENoiseIsoFilter");
	 metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	 metFlags.push_back("Flag_goodVertices");
	 metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
	// metFlags.push_back("Flag_METFilters");
	 metFlags.push_back("Flag_eeBadScFilter");

  //    }


	bool METflag = metFiltersPasses2(analysisTree, metFlags);
	met_flag = METflag;
//	if (!METflag && isData) continue;


      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();

     
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;




      double MET = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      METV.SetPx(analysisTree.pfmet_ex);	      
      METV.SetPy(analysisTree.pfmet_ey);


      if (!isData) 
	{
	  if (applyPUreweighting)	 {
	    puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	//	puweight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
	    weight *=puweight;      
	    pu_weight = puweight;
	  }

	}

  
      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      bool trigAccept = false;
	

      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;



      if (isData){
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//		std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}

      }

      }

      if (!isData ) isMainTrigger = true;


      if (!isMainTrigger) {
	std::cout << "HLT filter for SingleEle " << MainTrigger << " not found" << std::endl;
	return(-1);
      }

      /////now clear the Mu.El.Jets again to fill them again after cleaning
      MuMV.clear();
      ElMV.clear();
      TauMV.clear();
      LeptMV.clear();


    vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
//	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
//	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	if (applyElectronId && !analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie]) continue;
	if (applyElectronId && !analysisTree.electron_pass_conversion[ie]) continue;
	if (applyElectronId && analysisTree.electron_nmissinginnerhits[ie]>1) continue;
	if (fabs(analysisTree.electron_charge[ie]) !=1) continue;
	 electrons.push_back((int)ie);

      }
      if (electrons.size()==0 ) continue;



      vector<int> taus; taus.clear();

      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (analysisTree.tau_decayModeFinding[it]<decayModeFinding) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
        if ( fabs(analysisTree.tau_charge[it]) != 1 ) continue;
	  taus.push_back((int)it);

	}


      if (taus.size()==0) continue;


      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;


      float isoElecMin  = 1e+10;
      float isoTauMin = 1.; 
      if (CutBasedTauId) isoTauMin = 10e+10;
      if (!CutBasedTauId) isoTauMin = -10;
      float ptEl = 0;
      float ptTau = 0;
      //      if (electrons.size()>1||electrons.size()>1)
      //      std::cout << "electrons = " << electrons.size() << "  taus = " << taus.size() << std::endl;
      //
      for (unsigned int im=0; im<electrons.size(); ++im) {
	bool isEleLegMatch = false;
	//	bool isElectronTauSingleElectronFilterNameMatch = false;
	//	bool isElectronTauOverlapElectronMatch = false;
	unsigned int eIndex  = electrons.at(im);
	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[eIndex];
	float photonIsoElec = analysisTree.electron_photonIso[eIndex];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[eIndex];
	float puIsoElec = analysisTree.electron_puIso[eIndex];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[eIndex];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[eIndex];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/float(analysisTree.electron_pt[eIndex]);


	if (isData)
	{ for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      &&analysisTree.electron_pt[eIndex]>ptElectronCut&&
	      analysisTree.trigobject_pt[iT]>SingleElectronTriggerPtCut) { // IsoElec Leg
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) 
	      isEleLegMatch = true;
	    
	  }
	 }
	}
      if (!isData) isEleLegMatch = true;

      if (!isEleLegMatch) continue;

	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex]);

	  if (dR<dRleptonsCuteltau) continue;

	
	  float isoTau =1.;

	if (CutBasedTauId){
	isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex];

	  if ( (int) eIndex!= (int)el_index) {
	    if (relIsoElec==isoElecMin) {
	      if (analysisTree.electron_pt[eIndex]>ptEl) {
		isoElecMin  = relIsoElec;
		ptEl = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }
	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptEl = analysisTree.electron_pt[eIndex];
	      el_index = (int)eIndex;
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tau_index = (int)tIndex;
	    }
	  }
	  else {
	    if (isoTau==isoTauMin) {
	      if (analysisTree.tau_pt[tIndex]>ptTau) {
		ptTau = analysisTree.tau_pt[tIndex];
		isoTauMin = isoTau;
		tau_index = (int)tIndex;
	      }
	    }
	    else if (isoTau<isoTauMin) {
	      ptTau = analysisTree.tau_pt[tIndex];
	      isoTauMin = isoTau;
	      tau_index = (int)tIndex;
	    }
	  }
	  
	}



	if (!CutBasedTauId){
   isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
   //isoTau = analysisTree.tau_chargedIsoPtSum[tIndex];

	  if ( (int)eIndex!= (int)el_index) {
	    if (relIsoElec==isoElecMin) {
	      if (analysisTree.electron_pt[eIndex]>ptEl) {
		isoElecMin  = relIsoElec;
		ptEl = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }

	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptEl = analysisTree.electron_pt[eIndex];
	      el_index = (int)eIndex;
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tau_index = (int)tIndex;
	    }
          }
          else {
            if (isoTau==isoTauMin) {
              if (analysisTree.tau_pt[tIndex]>ptTau) {
                ptTau = analysisTree.tau_pt[tIndex];
                isoTauMin = isoTau;
                tau_index = (int)tIndex;
              }
            }
            else if (isoTau>isoTauMin) {
              ptTau = analysisTree.tau_pt[tIndex];
              isoTauMin = isoTau;
              tau_index = (int)tIndex;
            }
          }

        }

      }
      }

      //      std::cout << "eIndex = " << el_index << "   tau_index = " << tau_index << std::endl;
      //

      bool TauId = false;

            if ( analysisTree.tau_againstElectronTightMVA6[tau_index]>0.5 &&   analysisTree.tau_againstMuonLoose3[tau_index]>0.5) TauId = true;

	if (!TauId) continue;

      if ((int)tau_index<0) continue;
      if ((int)el_index<0) continue;

	bool tauPass =false;float isoTau = -999;

	if (!CutBasedTauId)
	{
		tauPass=
	  	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index] > 0.5;
	  	  //analysisTree.tau_chargedIsoPtSum[tau_index] < 0.8;

 
       isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
       ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];
       //isoTau = analysisTree.tau_chargedIsoPtSum[tau_index];
       //ta_IsoFlag=analysisTree.tau_chargedIsoPtSum[tau_index];
	}

	if (CutBasedTauId)
	{
		tauPass=
	         analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] > 0.5;

          isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
          ta_IsoFlag=analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index];
	}
	
     if (!tauPass) continue;

      ta_relIso[0]= isoTauMin;
      el_relIso[0] = isoElecMin;

      double q = analysisTree.tau_charge[tau_index] * analysisTree.electron_charge[el_index];
      event_sign  = q;


	double dReltau = deltaR(analysisTree.tau_eta[(int)tau_index],analysisTree.tau_phi[(int)tau_index],
				analysisTree.electron_eta[(int)el_index],analysisTree.electron_phi[(int)el_index]);
	if (dReltau < 0.5) continue;

      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;
      
  event_secondLeptonVeto = false;
  event_thirdLeptonVeto = false;



      // looking for extra muon
      bool foundExtraMuon = false;
      for (unsigned int ie = 0; ie<analysisTree.muon_count; ++ie) {
	if (analysisTree.muon_pt[ie]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[ie])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[ie])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[ie])>dzVetoMuonCut) continue;
	//if (applyVetoMuonId && !analysisTree.muon_isMedium[ie]) continue;
	if (applyVetoMuonId && !analysisTree.muon_isICHEP[ie]) continue;
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[ie];
        float photonIsoMu = analysisTree.muon_photonIso[ie];
        float chargedHadIsoMu = analysisTree.muon_chargedHadIso[ie];
        float puIsoMu = analysisTree.muon_puIso[ie];
        if (isIsoR03) {
          neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[ie];
          photonIsoMu = analysisTree.muon_r04_sumPhotonEt[ie];
          chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[ie];
          puIsoMu = analysisTree.muon_r04_sumPUPt[ie];
        }
        double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
        double neutralIsoMu = max(double(0),neutralIsoMuN);
        float absIsoMu = chargedHadIsoMu + neutralIsoMu;
        float relIsoMu = absIsoMu/analysisTree.muon_pt[ie];
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      // looking for extra electron's (dielectron veto)
      bool foundExtraElectron = false;
      vector<unsigned int> e_dielectrons; e_dielectrons.clear(); 
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {

      	      if ((int)im==(int)el_index) continue;

	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[im];
	float photonIsoElec = analysisTree.electron_photonIso[im];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[im];
	float puIsoElec = analysisTree.electron_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[im];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[im];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[im];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[im];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/analysisTree.electron_pt[im];

	if (analysisTree.electron_pt[im]>ptDilepElectronCut&&
	    fabs(analysisTree.electron_eta[im])<etaDilepElectronCut&&
	    fabs(analysisTree.electron_dxy[im])<dxyDilepElectronCut&&
	    fabs(analysisTree.electron_dz[im])<dzDilepElectronCut&&
	    analysisTree.electron_cutId_veto_Spring15[im]&&
	    relIsoElec<isoDilepElectronCut && 
	    fabs(analysisTree.electron_charge[im]) ==1)
	{
	    
	float dRelectrons = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			   analysisTree.electron_eta[im],analysisTree.electron_phi[im]);

	    if (dRelectrons>dRDilepVetoCut && (analysisTree.electron_charge[el_index]*analysisTree.electron_charge[im]<0.)) 
	      dilepton_veto = true;

	}
	 // e_dielectrons.push_back(im);

	if (analysisTree.electron_pt[im]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[im];
	if (applyVetoElectronId && !electronMvaId) continue;
	if (applyVetoElectronId && !analysisTree.electron_pass_conversion[im]) continue;
	if (applyVetoElectronId && analysisTree.electron_nmissinginnerhits[im]>1) continue;
	if (relIsoElec>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

//      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;
/*
      dilepton_veto = false;
      if (e_dielectrons.size()>1) {
	for (unsigned int i1=0; i1<e_dielectrons.size()-1; ++i1) {
	  unsigned int indx1 = e_dielectrons[i1];
	  for (unsigned int i2=i1+1; i2<e_dielectrons.size(); ++i2 ) {
	    unsigned int indx2 = e_dielectrons[i2];
	    float dRelectrons = deltaR(analysisTree.electron_eta[indx1],analysisTree.electron_phi[indx1],
				   analysisTree.electron_eta[indx2],analysisTree.electron_phi[indx2]);
	    if (dRelectrons>dRDilepVetoCut && (analysisTree.electron_charge[indx1]*analysisTree.electron_charge[indx2]<0)) 
	      dilepton_veto = true;
 	  }
	}
      }
      //      cout << analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] << endl;
      //      cout << analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index] << endl;
*/
   	event_secondLeptonVeto = dilepton_veto;
//	if (dilepton_veto)  continue;


      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;

      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

      /////////////////Make sure that the selected electron fired the trigger and it is within a Dr < 0.5


      /////////////////////////////////
      // Apply trigger SF
      // ////////////////////////////
      double ptEl1 = analysisTree.electron_pt[el_index];
      double etaEl1 = analysisTree.electron_eta[el_index];
      float trigweight=1.;


      float EffFromData = (float)SF_electronTrigger->get_EfficiencyData(double(ptEl1),double(etaEl1));
      /*float El22EffMC   = (float)SF_electronTrigger->get_EfficiencyMC(double(ptEl1),double(etaEl1));*/

     /* if (!isData) {
	if (El22EffMC>1e-6)
	  trigweight = EffFromData / El22EffMC;
	if (!isData &&  ( string::npos != filen.find("stau") || string::npos != filen.find("C1")) )  trigweight = EffFromData;
	weight *= trigweight;
	trig_weight = trigweight;
      }
	*/

	if (!isData) trigweight = EffFromData;
	weight *= trigweight;
	trig_weight = trigweight;


      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;

	///////////Lepton SF
      if (!isData && applyLeptonSF) {

	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
      double ptEl1 = analysisTree.electron_pt[el_index];
      double etaEl1 = analysisTree.electron_eta[el_index];
	double IdIsoSF_el1 = SF_elIdIso->get_ScaleFactor(ptEl1, etaEl1);

	ElSF_IdIso_El1H->Fill(IdIsoSF_el1);
	weight *= IdIsoSF_el1;
	LSF_weight = IdIsoSF_el1;
      }

      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;
 
      ///////////////Check if the selected tau has a gen matched - if not, apply Tau Fake Rate
      bool isTauMatched = false;
      bool isGenLeptonMatched = false;
	if (!isData){
      TLorentzVector genTauV;  
	TLorentzVector genLepV;  

	for (unsigned int gt = 0 ; gt < analysisTree.gentau_count;++gt){

	 // genTauV.SetXYZT(0.,0.,0.,0.);
	  genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);


	  double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    genTauV.Eta(), genTauV.Phi());


	  if (Drr < 0.2 && analysisTree.gentau_isPrompt[gt] > 0.5  && genTauV.Pt() > 15. ) isTauMatched = true;

	}
      
      
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      		  if ( (abs(analysisTree.genparticles_pdgid[igen])==11 || abs(analysisTree.genparticles_pdgid[igen])==13) && analysisTree.genparticles_isPrompt[igen] > 0.5){

	  genLepV.SetXYZT(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen], analysisTree.genparticles_e[igen]);

	  double Drm=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  genLepV.Eta(),genLepV.Phi());

		if (Drm < 0.2 && genLepV.Pt() > 8. ) isGenLeptonMatched = true;

		}
      
	  }
      
      
      }//!isData

	genTauMatched = isTauMatched;
	genLeptonMatched = isGenLeptonMatched;

	/////////TFR
 /*     if (!isData && applyTFR && !isTauMatched) {


	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	

	double ptTau1 = (double)analysisTree.tau_pt[tau_index];
	double etaTau1 = (double)analysisTree.tau_eta[tau_index];
	double TFRSF_el1 = SF_TFR->get_ScaleFactor(ptTau1, etaTau1);

	//MuSF_IdIso_Mu1H->Fill(TFRSF_mu1);
	//weight *= TFRSF_el1;
	TFR_weight  = TFRSF_el1;
	//cout<<"  "<<TFRSF_mu1<<"  for  eta  "<<etaTau1<<  " pT  "<< ptTau1<<endl;
      	}*/

      CFCounter[iCut]+= weight;
      CFCounter_[iCut]+= weight;
      iCFCounter[iCut]++;
      iCut++;


      ///////////////////////////////////////////////////////////
      //////////////////////////////////////////////
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

        el_neutralHadIso[ie] = analysisTree.electron_r03_sumNeutralHadronEt[ie];
        el_photonIso[ie] = analysisTree.electron_r03_sumPhotonEt[ie];
        el_chargedHadIso[ie] = analysisTree.electron_r03_sumChargedHadronPt[ie];
        el_puIso[ie] = analysisTree.electron_r03_sumPUPt[ie];
 
        double neutralIso = el_neutralHadIso[ie] + el_photonIso[ie] - 0.5*el_puIso[ie];
        neutralIso = max(double(0),neutralIso);
        el_neutralIso[ie] = neutralIso ;
        el_absIsoEl[ie] = el_chargedHadIso[ie] + mu_neutralIso[ie];
	el_relIsoEl[ie]  = el_absIsoEl[ie]/el_pt[ie] ;

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
     	ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[it];
     	ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[it];
     	ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[it];
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

      if (!isData) npartons = analysisTree.genparticles_noutgoing;

      all_weight = weight;


//systematics and other 

     met_ex_JetEnUp = analysisTree.pfmet_ex_JetEnUp;
     met_ey_JetEnUp = analysisTree.pfmet_ey_JetEnUp;

     met_ex_JetEnDown = analysisTree.pfmet_ex_JetEnDown;
     met_ey_JetEnDown = analysisTree.pfmet_ey_JetEnDown;

     met_ex_UnclusteredEnUp = analysisTree.pfmet_ex_UnclusteredEnUp;
     met_ey_UnclusteredEnUp = analysisTree.pfmet_ey_UnclusteredEnUp;
   
     met_ex_UnclusteredEnDown = analysisTree.pfmet_ex_UnclusteredEnDown;
     met_ey_UnclusteredEnDown = analysisTree.pfmet_ey_UnclusteredEnDown;

      TLorentzVector leptonsV, muonJ, jetsLV;

      //      continue;

      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      bool btagged= false;


      JetsMV.clear();

      float jetEtaCut = 2.4;
      float jetEta = 2.4;
      float DRmax = 0.5;
      int countjets = 0;
      bool dRmuJet = false;
      bool dRtauJet = false;
      float bJetEtaCut = ptJetCut;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;
      float ptLeadingBJet = -1;

	int counter_cleaned_jets = 0;



      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);

	if (absJetEta > etaJetCut) continue;
	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


	bool isPFJetId = false ; 
	isPFJetId =looseJetiD(analysisTree,jet);

	//				cout<<" 1- jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  iEntry "<<iEntry<<endl;
	if (!isPFJetId) continue;
	jet_isLoose[jet] = isPFJetId;
	bool cleanedJet = true;
	//				cout<<"  jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  "<<iEntry<<endl;

	double Dre=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if (  Dre  < DRmax )  cleanedJet=false;

	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if ( Drr < DRmax) cleanedJet=false;


	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

	if (analysisTree.pfjet_btag[jet][0]  > bTag) btagged = true;
	


	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);

	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;

	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_BC.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_BC.eval_auto_bounds("central",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
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
	  }
//if (btagged)
//	cout<<"  what here "<<btagged<<"  "<<cleanedJet<<endl;
	  if (btagged) {
	    if (cleanedJet) { 
	      bjets.push_back(jet);
	      if (ptJetCut>ptLeadingBJet) {
		ptLeadingBJet = ptJetCut;
		indexLeadingBJet = jet;
	      }
	    }
	  }
	}


	if (!cleanedJet) continue;

	   counter_cleaned_jets++;

	  jets.push_back(jet);
	  jets_cleaned[counter_cleaned_jets]=jet;

	if (indexLeadingJet>=0) {
	  if (ptJetCut<ptLeadingJet&&ptJetCut>ptSubLeadingJet) {
	    indexSubLeadingJet = jet;
	    ptSubLeadingJet = ptJetCut;
	  }
	}

	if (ptJetCut>ptLeadingJet) {
	  indexLeadingJet = jet;
	  ptLeadingJet = ptJetCut;
	}
		JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
	JetsMV.push_back(JetsV);
	//jet jetc
	jet_jecUn[counter_cleaned_jets++] = analysisTree.pfjet_jecUncertainty[jet];

      }
      njets = jets.size();
      countjets = jets.size();
      jet_count = jets.size();
      //njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      //nbtag_nocleaned = bjets_nocleaned.size();

      npv =  analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
	
      SusyMother = SusyMotherMassF;
      SusyLSP = SusyLSPMassF;

      T->Fill();


      continue;
      /////////////////////////////////////////////////



      selEvents++;
    } // end of file processing (loop over events in one file)



    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;


  cout<<" Total events  "<<nEvents<<"  Will use weight  "<<histWeightsH->GetSumOfWeights()<<" Norm Factor for a Lumi of "<<Lumi<<"/pb is "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<endl;
  cout<<" First content "<<CFCounter[0]<<endl;
  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;

  //write out cutflow
  ofstream tfile;
  TString outname=argv[2];
  TString textfilename = "cutflow_"+outname+"_"+Channel+"_"+argv[3]+".txt";
  //  tfile.open(textfilename);
  cout<<" iCFCounter is the pure event count, CFCounter is normalized to "<<Lumi<<" /pb and CFCounter_ is the sum of weights for a given bin"<<endl;
  //  tfile << "########################################" << endl;
  for(int ci = 0; ci < CutNumb; ci++)
    {
      // tfile << CutList[ci]<<"\t & \t"
      //	    << CFCounter[ci]  <<"\t & \t"<< statUnc[ci] <<"\t & \t"<< iCFCounter[ci] << endl;
      CFCounter_[iCut] = iCFCounter[iCut];
      CutFlowUnW->SetBinContent(1+ci,0);
      CutFlow->SetBinContent(1+ci,0);
      CutFlowUnW->SetBinContent(1+ci,float(CFCounter[ci]) );


      CFCounter[ci] *= float(XSec*Lumi/( histWeightsH->GetSumOfWeights()));

      CutFlow->SetBinContent(1+ci,float(CFCounter[ci]));

      cout << " i "<<ci<<" iCFCounter "<<iCFCounter[ci]<<" CFCounter  "<<CFCounter[ci]<<" CFCounter_ " <<CFCounter_[ci]<<"  "<<XSec*Lumi/( histWeightsH->GetSumOfWeights())<<" UnWBinContent "<<CutFlowUnW->GetBinContent(1+ci)<<" Norm to Lumi Content "<<CutFlow->GetBinContent(1+ci)<<" Cut "<<CutList[ci]<<endl;   
      if (iCFCounter[ci] <0.2) statUnc[ci] =0;
      else statUnc[ci] = sqrt(CFCounter[ci]);
    }



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
  histTopPt->Write();
  histRuns->Write();
  CutFlowUnW->Write();
  CutFlow->Write();
  ElSF_IdIso_El1H->Write();
  file->Write();
  file->Close();


  delete file;

}




