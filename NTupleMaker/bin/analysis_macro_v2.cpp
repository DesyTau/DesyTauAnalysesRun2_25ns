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
  string Channel="mumu";



  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  

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


////////////muons

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");

  // topological cuts
  const float dRMuonsCut   = cfg.get<float>("dRMuonsCut");
  const bool sameSign      = cfg.get<bool>("SameSignMuons");

  // track selection
  const float dRIsoMuon       = cfg.get<float>("dRIsoMuon");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  // trigger
  const string dimuonTriggerName = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName1 = cfg.get<string>("MuonLowPtFilterName1");
  const string muonLowPtFilterName2 = cfg.get<string>("MuonLowPtFilterName2");
  const string dimuonDzFilterName = cfg.get<string>("DimuonDzFilterName");
  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 

  TString DiMuonTriggerName(dimuonTriggerName);
  TString MuonHighPtFilterName(muonHighPtFilterName);
  TString MuonLowPtFilterName1(muonLowPtFilterName1);
  TString MuonLowPtFilterName2(muonLowPtFilterName2);
  TString DiMuonDzFilterName(dimuonDzFilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");




  const string dataBaseDir = cfg.get<string>("DataBaseDir");

/*
  const string MuonidIsoEffFile = cfg.get<string>("MuonidIsoEffFile");
  const string MuontrigEffFile = cfg.get<string>("MuontrigEffFile");



  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");

  // topSingleMuonTriggerFile

  const string TrigLeg  = cfg.get<string>("SingleMuonFilterName") ;
  const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");
*/

  // vertex distributions filenames and histname


  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
 

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
//  TString MainTrigger(TrigLeg);



  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();

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

  TString TStrName(rootFileName);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<"  The systematic will be "<<Systematic<<"  and save dir will be  "<<SaveDir<<endl;



// PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2016_271036-284044_80bins.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Moriond17_PU25ns_V1.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);


  // Lepton Scale Factors 

/*
  cout<<"  Initializing iD SF files....."<<endl;
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));

  cout<<"  Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuontrigEffFile));


  cout<<" ended initialization here "<<endl;
*/

  TFile * file;
  	file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"recreate");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());


  TH1F * muonCountH = new TH1F("muonCountH","",11,-0.5,10.5);
  TH1F * nGoodMuonsH = new TH1F("nGoodMuonsH","",11,-0.5,10.5);

  // histograms after dimuon selection
  TH1F * ptLeadingMuH = new TH1F("ptLeadingMuH","",50,0,100);
  TH1F * ptTrailingMuH = new TH1F("ptTrailingMuH","",50,0,100);
  TH1F * etaLeadingMuH = new TH1F("etaLeadingMuH","",50,-2.5,2.5);
  TH1F * etaTrailingMuH = new TH1F("etaTrailingMuH","",50,-2.5,2.5);
  TH1F * dimuonMassH = new TH1F("dimuonMassH","",500,0,500);
  TH1F * nTracksLeadingMuH = new TH1F("nTracksLeadingMuH","",21,-0.5,20.5);
  TH1F * nTracksTrailingMuH = new TH1F("nTracksTrailingMuH","",21,-0.5,20.5);
  TH1F * nSigTracksLeadingMuH = new TH1F("nSigTracksLeadingMuH","",21,-0.5,20.5);
  TH1F * nSigTracksTrailingMuH = new TH1F("nSigTracksTrailingMuH","",21,-0.5,20.5);

  TH1F * ptTrackLeadingMuH = new TH1F("ptTrackLeadingMuH","",100,0,100);
  TH1F * etaTrackLeadingMuH = new TH1F("etaTrackLeadingMuH","",50,-2.5,2.5);
  TH1F * dxyTrackLeadingMuH = new TH1F("dxyTrackLeadingMuH","",200,-0.5,0.5);
  TH1F * dzTrackLeadingMuH = new TH1F("dzTrackLeadingMuH","",200,-1,1);

  TH1F * ptTrackTrailingMuH = new TH1F("ptTrackTrailingMuH","",100,0,100);
  TH1F * etaTrackTrailingMuH = new TH1F("etaTrackTrailingMuH","",50,-2.5,2.5);
  TH1F * dxyTrackTrailingMuH = new TH1F("dxyTrackTrailingMuH","",200,-0.5,0.5);
  TH1F * dzTrackTrailingMuH = new TH1F("dzTrackTrailingMuH","",200,-1,1);

  TH1F * ptLeadingMuSelH = new TH1F("ptLeadingMuSelH","",50,0,100);
  TH1F * ptTrailingMuSelH = new TH1F("ptTrailingMuSelH","",50,0,100);
  TH1F * etaLeadingMuSelH = new TH1F("etaLeadingMuSelH","",50,-2.5,2.5);
  TH1F * etaTrailingMuSelH = new TH1F("etaTrailingMuSelH","",50,-2.5,2.5);
  TH1F * dimuonMassSelH = new TH1F("dimuonMassSelH","",500,0,500);
  TH1F * invMass2Mu2TrkSelH = new TH1F("invMass2Mu2TrkSelH","",500,0,500);

  TH1F * massLeadingMuTrkH = new TH1F("massLeadingMuTrkH","",500,0,50);
  TH1F * massTrailingMuTrkH = new TH1F("massTrailingMuTrkH","",500,0,50);
  TH1F * massMuTrkH = new TH1F("massMuTrkH","",500,0,50);
  TH2F * m1m2SelH = new TH2F("m1m2SelH","",20,0,20,20,0,20);

  // Counters
  TH1F * counter_InputEventsH=new TH1F("counter_InputEventsH","",1,0.,2.);
  TH1F * counterMuonIDH=new TH1F("counterMuonIDH","",1,0.,2.);
  TH1F * counter_MuonIPH=new TH1F("counter_MuonIPH","",1,0.,2.);
  TH1F * counter_MuonEtaH=new TH1F("counter_MuonEtaH","",1,0.,2.);
  TH1F * counter_MuonSizeGTE2H=new TH1F("counter_MuonSizeGTE2H","",1,0.,2.);
  TH1F * counter_MuonSameSignH=new TH1F("counter_MuonSameSignH","",1,0.,2.);
  TH1F * counter_dRMuonH=new TH1F("counter_dRMuonH","",1,0.,2.);
  TH1F * counter_Mu1_Mu2TrigH=new TH1F("counter_Mu1_Mu2TrigH","",1,0.,2.);
  TH1F * counter_MaxMuon_ptSumH=new TH1F("counter_Muon_ptSumH","",1,0.,2.);
  TH1F * counter_MuonKinematicsH=new TH1F("counter_MuonKinematicsH","",1,0.,2.);         
  TH1F * counter_nTracksH=new TH1F("counter_nTracksH","",1,0.,2.);
  TH1F * counter_OneTrackOnlyH=new TH1F("counter_OneTrackOnlyH","",1,0.,2.);                  
  TH1F * counter_ChargeReqLeadMuonTrkH=new TH1F("counter_ChargeReqLeadMuonTrkH","",1,0.,2.);         
  TH1F * counter_ChargeReqTrailMuonTrkH=new TH1F("counter_ChargeReqTrailMuonTrkH","",1,0.,2.);         
  TH1F * counter_TightTrkIPReqH=new TH1F("counter_TightTrkIPReqH","",1,0.,2.);         
  TH1F * counter_TightTrkptReqH=new TH1F("counter_TightTrkptReqH","",1,0.,2.);         
  TH1F * counter_FinalEventsH=new TH1F("counter_FinalEventsH","",1,0.,2.);         

   // Background studies
  TH1F * InvMassN23leadingH = new TH1F("InvMassN23leadingH","",10,0.,10.);
  TH1F * InvMassN23trailingH = new TH1F("InvMassN23trailingH","",10,0.,10.);
  TH1F * InvMassN23H = new TH1F("InvMassN23H","",10,0.,10.);

  TH1F * InvMassHardestNtrk23leadingH = new TH1F("InvMassHardestNtrk23leadingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk23trailingH = new TH1F("InvMassHardestNtrk23trailingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk23H = new TH1F("InvMassHardestNtrk23H","",10,0.,10.);

  TH1F * InvMassSoftestNtrk23leadingH = new TH1F("InvMassSoftestNtrk23leadingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk23trailingH = new TH1F("InvMassSoftestNtrk23trailingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk23H = new TH1F("InvMassSoftestNtrk23H","",10,0.,10.);

  TH1F * InvMassHardestNtrk1leadingH = new TH1F("InvMassHardestNtrk1leadingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk1trailingH = new TH1F("InvMassHardestNtrk1trailingH","",10,0.,10.);
  TH1F * InvMassHardestNtrk1H = new TH1F("InvMassHardestNtrk1H","",10,0.,10.);

  TH1F * InvMassSoftestNtrk1leadingH = new TH1F("InvMassSoftestNtrk1leadingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk1trailingH = new TH1F("InvMassSoftestNtrk1trailingH","",10,0.,10.);
  TH1F * InvMassSoftestNtrk1H = new TH1F("InvMassSoftestNtrk1H","",10,0.,10.);

  TH1F * InvMassLeadingH = new TH1F("InvMassLeadingH","",10,0.,10.);
  TH1F * InvMassTrailingH = new TH1F("InvMassTrailingH","",10,0.,10.);
  TH1F * InvMassH = new TH1F("InvMassH","",10,0.,10.);
  TH2F * InvMass2DH = new TH2F("InvMass2DH","",10,0.,10.,10,0.,10.);

  // Correlation Plots
  TH1F * InvMassTrackPlusMuon1D_ControlH = new TH1F("InvMassTrackPlusMuon1D_ControlH","",100,0.,10.); 
  TH2F * InvMassTrackPlusMuon2D_ControlLeadH = new TH2F("InvMassTrackPlusMuon2D_ControlLeadH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_ControlTrailH = new TH2F("InvMassTrackPlusMuon2D_ControlTrailH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_ControlBothH = new TH2F("InvMassTrackPlusMuon2D_ControlBothH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_ControlH = new TH2F("InvMassTrackPlusMuon2D_ControlH","",100,0.,10.,100,0.,10.);












  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  
  //SetupTree(); 

  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
 
for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));


bool WithInit = true;

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
    for (Long64_t iCand=0; iCand<numberOfEntriesInitTree; ++iCand) {
      _inittree->GetEntry(iCand);
      if (isData)
	histWeightsH->Fill(0.,1.);
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
		for (Long64_t iCand=0; iCand<numberOfEntries; ++iCand) 
			{
			analysisTree.GetEntry(iCand);
			histWeightsH->Fill(0.,analysisTree.genweight);
			}
		}
  	float genweights=1.;

    if(!isData && WithInit) 
      {
	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	genweightsTree->SetBranchAddress("genweight",&genweights);


	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iCandInit=0; iCandInit<numberOfEntriesInit; ++iCandInit) { 
	  genweightsTree->GetEntry(iCandInit);
	  histWeightsH->Fill(0.,genweights);
	}
    
      }




    for (Long64_t iCand=0; iCand<numberOfEntries; ++iCand) { 


      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iCand);
      nEvents++;

     counter_InputEventsH->Fill(1.);

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


      histTopPt->Fill(0.,topptweight);
      histTopPtSq->Fill(0.,topptweight*topptweight);

	}
	  if (!isData ) {
	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;
	   // std::cout <<"analysisTree.genweight "<< float(analysisTree.genweight) << std::endl;
	  lumi=true;
	  }
					

      if (isData)  {
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
	if (!METflag && isData) continue;



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

     unsigned int nMu8Leg   = 0;
     unsigned int nMu17Leg  = 0;
     unsigned int nDZFilter = 0;
     bool isMu8Leg = false;
     bool isMu17Leg = false;
     bool isDZFilter = false;


      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));

       if (HLTFilter==MuonHighPtFilterName) {
   nMu17Leg = i;
   isMu17Leg = true;
  }
       if (HLTFilter==MuonLowPtFilterName1 || HLTFilter==MuonLowPtFilterName2) {
   nMu8Leg = i;
   isMu8Leg = true;
  }
       if (HLTFilter==DiMuonDzFilterName) {
   nDZFilter = i;
   isDZFilter = true;
  }

      }

      if (!isMu17Leg) {
	cout << "Filter " << MuonHighPtFilterName << " not found " << endl;
	exit(-1);
      }
      if (!isMu8Leg) {
	cout << "Filter " << MuonLowPtFilterName1 << " and  "<<MuonLowPtFilterName2<< "  not found " << endl;
	exit(-1);
      }
      if (!isDZFilter) {
	cout << "Filter " << DiMuonDzFilterName << " not found " << endl;
	exit(-1);
      }
      //      std::cout << std::endl;

      muonCountH->Fill(float(analysisTree.muon_count),weight);


      vector<unsigned int> muons; muons.clear();
      for(UInt_t i=0;i<analysisTree.muon_count;i++){
	if(!analysisTree.muon_isLoose[i]) continue;
	if(!analysisTree.muon_isPF[i]) continue;
	counterMuonIDH->Fill(1.);
	if(fabs(analysisTree.muon_dxy[i])>dxyMuonCut) continue;
	if(fabs(analysisTree.muon_dz[i])>dzMuonCut) continue;
	counter_MuonIPH->Fill(1.);
	if(analysisTree.muon_pt[i]<ptMuonLowCut) continue;
	if(fabs(analysisTree.muon_eta[i])>etaMuonLowCut) continue;
	counter_MuonEtaH->Fill(1.);
	//      cout << "muon pt = " << muon_pt[i] << endl;
	muons.push_back(i);
      }

      nGoodMuonsH->Fill(float(muons.size()),weight);
      
      if (muons.size()<2) continue; // quit event if number of good muons < 2
           counter_MuonSizeGTE2H->Fill(1.);
      float maxPtSum = -1;
      int iLeading = -1;
      int iTrailing = -1;
      for (unsigned int i1=0; i1<muons.size()-1; ++i1) {
	int index1 = muons.at(i1);
	for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {
	  int index2 = muons.at(i2);
	  float ptSum = analysisTree.muon_pt[index1] + analysisTree.muon_pt[index2];
	  float charge = analysisTree.muon_charge[index1] *  analysisTree.muon_charge[index2];
	  bool chargeSelection = charge<0;
	  if (sameSign)
	    chargeSelection = charge>0;
	  if (!chargeSelection) continue;
	  counter_MuonSameSignH->Fill(1.);
	  float dRmuons = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				 analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);

	  if (dRmuons<dRMuonsCut) continue;
	  counter_dRMuonH->Fill(1.);    
	  bool mu1MatchMu17 = false;
	  bool mu1MatchMu8  = false;
	  bool mu1MatchDz   = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMu17Leg] && analysisTree.muon_pt[index1]>ptMuonHighCut && fabs(analysisTree.muon_eta[index1])<etaMuonHighCut)
	      mu1MatchMu17 = true;
	    if (analysisTree.trigobject_filters[iT][nMu8Leg] && analysisTree.muon_pt[index1]>ptMuonLowCut && fabs(analysisTree.muon_eta[index1])<etaMuonLowCut)
	      mu1MatchMu8 = true;
	    if (analysisTree.trigobject_filters[iT][nDZFilter])
	      mu1MatchDz = true;
	  }
	  bool mu2MatchMu17 = false;
	  bool mu2MatchMu8  = false;
	  bool mu2MatchDz   = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMu17Leg] && analysisTree.muon_pt[index2]>ptMuonHighCut && fabs(analysisTree.muon_eta[index2])<etaMuonHighCut)
	      mu2MatchMu17 = true;
	    if (analysisTree.trigobject_filters[iT][nMu8Leg] && analysisTree.muon_pt[index2]>ptMuonLowCut && fabs(analysisTree.muon_eta[index2])<etaMuonLowCut)
	      mu2MatchMu8 = true;
	    if (analysisTree.trigobject_filters[iT][nDZFilter])
	      mu2MatchDz = true;
	  }
	  
	  // trigger condition
	  bool isTriggerMatched = ((mu1MatchMu17&&mu2MatchMu8)||(mu1MatchMu8&&mu2MatchMu17)) && mu1MatchDz && mu2MatchDz;
	  if (!isTriggerMatched) continue;
	  counter_Mu1_Mu2TrigH->Fill(1.);                 
	  if (ptSum>maxPtSum) {
	    maxPtSum = ptSum;
	    if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]) {
	      iLeading = index1;
	      iTrailing = index2;
	    }
	    else {
	      iLeading = index2;
	      iTrailing = index1;
	    }
	  }
	}
      }

      if (iLeading<0) continue;
      if (iTrailing<0) continue;
          counter_MaxMuon_ptSumH->Fill(1.);                 

      // dimuon selection passed 
      TLorentzVector LeadingMuon4; LeadingMuon4.SetXYZM(analysisTree.muon_px[iLeading],
              analysisTree.muon_py[iLeading],
              analysisTree.muon_pz[iLeading],
              MuMass);

      TLorentzVector TrailingMuon4; TrailingMuon4.SetXYZM(analysisTree.muon_px[iTrailing],
                analysisTree.muon_py[iTrailing],
                analysisTree.muon_pz[iTrailing],
                MuMass);

      TLorentzVector diMuon4 = LeadingMuon4 + TrailingMuon4;

      float dimuonMass = diMuon4.M();

      // filling histograms (muon kinematics)
      dimuonMassH->Fill(dimuonMass,weight);
      ptLeadingMuH->Fill(analysisTree.muon_pt[iLeading],weight);
      ptTrailingMuH->Fill(analysisTree.muon_pt[iTrailing],weight);
      etaLeadingMuH->Fill(analysisTree.muon_eta[iLeading],weight);
      etaTrailingMuH->Fill(analysisTree.muon_eta[iTrailing],weight);
      counter_MuonKinematicsH->Fill(1.);                  

      // counting tracks around each muon
      std::vector<unsigned int> trkLeadingMu; trkLeadingMu.clear(); // all tracks
      std::vector<unsigned int> trkTrailingMu; trkTrailingMu.clear(); // all tracks
      std::vector<unsigned int> trkSigLeadingMu; trkSigLeadingMu.clear(); // signal tracks
      std::vector<unsigned int> trkSigTrailingMu; trkSigTrailingMu.clear(); // signal tracks
      unsigned int hardestTrkLeading = 0; // index of hardest track around leading mu
      unsigned int hardestTrkTrailing = 0; // index of hardest track around trailing mu
      unsigned int softestTrkLeading = 0; // index of softest track around leading mu
      unsigned int softestTrkTrailing = 0; // index of softest track around trailing mu

      float ptHardestLeading = 0;
      float ptHardestTrailing = 0;
      float ptSoftestLeading = 1e+10;
      float ptSoftestTrailing = 1e+10;
      std::vector<unsigned int> Soft_trkLeadingMu; Soft_trkLeadingMu.clear(); // Lead Muon tracks control region
      std::vector<unsigned int> Soft_trkTrailingMu; Soft_trkTrailingMu.clear(); // Trail Muon tracks control region


      for (unsigned int iTrk=0; iTrk<analysisTree.track_count; ++iTrk) {
	if (fabs(analysisTree.track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
	if (fabs(analysisTree.track_dxy[iTrk])>dxyTrkLooseCut) continue;
	if (fabs(analysisTree.track_dz[iTrk])>dzTrkLooseCut) continue;
	if (fabs(analysisTree.track_eta[iTrk])>etaTrkCut) continue;
	if (fabs(analysisTree.track_pt[iTrk])<ptTrkLooseCut) continue;
  
	TLorentzVector trk4; trk4.SetXYZM(analysisTree.track_px[iTrk],
					  analysisTree.track_py[iTrk],
					  analysisTree.track_pz[iTrk],
					  analysisTree.track_mass[iTrk]);
	
	TLorentzVector leadingMuDiff = LeadingMuon4 - trk4;
	if (leadingMuDiff.P()>0.1) { // track is not leading muon
	  float drTrkMu = deltaR(analysisTree.muon_eta[iLeading],analysisTree.muon_phi[iLeading],
				 analysisTree.track_eta[iTrk],   analysisTree.track_phi[iTrk]);
          float qTrkLeadingMu = analysisTree.track_charge[iTrk]*analysisTree.muon_charge[iLeading];
	  if (drTrkMu<dRIsoMuon){
	    trkLeadingMu.push_back(iTrk);
            if (analysisTree.track_pt[iTrk]>ptTrkLooseCut && analysisTree.track_pt[iTrk]< ptTrkCut)
              Soft_trkLeadingMu.push_back(iTrk);
	  }
	  if (drTrkMu<dRIsoMuon && qTrkLeadingMu<0 && fabs(analysisTree.track_dxy[iTrk])<dxyTrkCut && fabs(analysisTree.track_dz[iTrk])<dzTrkCut && analysisTree.track_pt[iTrk]>ptTrkCut) {
	    trkSigLeadingMu.push_back(iTrk);
	    if (analysisTree.track_pt[iTrk]>ptHardestLeading) {
	      ptHardestLeading = analysisTree.track_pt[iTrk];
	      hardestTrkLeading = iTrk;
	    }
	    if (analysisTree.track_pt[iTrk]<ptSoftestLeading) {
	      ptSoftestLeading = analysisTree.track_pt[iTrk];
	      softestTrkLeading = iTrk;
	    }
	  }
	}

	TLorentzVector trailingMuDiff = TrailingMuon4 - trk4;
	if (trailingMuDiff.P()>0.1) { // track is not trailing muon
	  float drTrkMu = deltaR(analysisTree.muon_eta[iTrailing],analysisTree.muon_phi[iTrailing],
				 analysisTree.track_eta[iTrk],analysisTree.track_phi[iTrk]);
          float qTrkTrailingMu = analysisTree.track_charge[iTrk]*analysisTree.muon_charge[iTrailing];         
	  if (drTrkMu<dRIsoMuon){
            trkTrailingMu.push_back(iTrk);
            if (analysisTree.track_pt[iTrk] > ptTrkLooseCut && analysisTree.track_pt[iTrk]< ptTrkCut)
              Soft_trkTrailingMu.push_back(iTrk);
	  }
	  if (drTrkMu<dRIsoMuon && qTrkTrailingMu<0 && fabs(analysisTree.track_dxy[iTrk])<dxyTrkCut && fabs(analysisTree.track_dz[iTrk])<dzTrkCut && analysisTree.track_pt[iTrk]>ptTrkCut) {
	    trkSigTrailingMu.push_back(iTrk);
	    if (analysisTree.track_pt[iTrk]>ptHardestTrailing) {
	      ptHardestTrailing = analysisTree.track_pt[iTrk];
	      hardestTrkTrailing = iTrk;
	    }
	    if (analysisTree.track_pt[iTrk]<ptSoftestTrailing) {
	      ptSoftestTrailing = analysisTree.track_pt[iTrk];
	      softestTrkTrailing = iTrk;
	    }
	    
	  }
	}

      }
      
      nTracksLeadingMuH->Fill(float(trkLeadingMu.size()),weight);
      nTracksTrailingMuH->Fill(float(trkTrailingMu.size()),weight);
      nSigTracksLeadingMuH->Fill(float(trkSigLeadingMu.size()),weight);
      nSigTracksTrailingMuH->Fill(float(trkSigTrailingMu.size()),weight);

      // defining sidebands and signal region

      // sideband N23
      bool isN23leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
      bool isN23trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==2||trkLeadingMu.size()==3); 

      // sidebands Ntrk23
      bool isNtrk23leading  = trkSigLeadingMu.size()>0 && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
      bool isNtrk23trailing = trkSigTrailingMu.size()>0 && (trkLeadingMu.size()==2||trkLeadingMu.size()==3);

      // sidebands Ntrk1
      bool isNtrk1leading  = trkSigLeadingMu.size()>0 && trkTrailingMu.size()==1;
      bool isNtrk1trailing = trkSigTrailingMu.size()>0 && trkLeadingMu.size()==1;

      // signal region
      bool signalRegion = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkSigTrailingMu.size()==1&&trkTrailingMu.size()==1);
      
      counter_nTracksH->Fill(1.);                 

      // sidebands N23
      if (isN23leading) {
	int iTrk = trkSigLeadingMu[0];
	TLorentzVector Track4; Track4.SetXYZM(analysisTree.track_px[iTrk],
					      analysisTree.track_py[iTrk],
					      analysisTree.track_pz[iTrk],
					      analysisTree.track_mass[iTrk]);
	TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
	float mass = TrackPlusMuon4.M();
	InvMassN23leadingH->Fill(mass,weight);
	InvMassN23H->Fill(mass,weight);
      }      
      if (isN23trailing) {
	int iTrk = trkSigTrailingMu[0];
	TLorentzVector Track4; Track4.SetXYZM(analysisTree.track_px[iTrk],
					      analysisTree.track_py[iTrk],
					      analysisTree.track_pz[iTrk],
					      analysisTree.track_mass[iTrk]);
	TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
	float mass = TrackPlusMuon4.M();
	InvMassN23trailingH->Fill(mass,weight);
	InvMassN23H->Fill(mass,weight);
      }      

      // sidebands Ntrk23
      if (isNtrk23leading) {
	TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(analysisTree.track_px[hardestTrkLeading],
							    analysisTree.track_py[hardestTrkLeading],
							    analysisTree.track_pz[hardestTrkLeading],
							    analysisTree.track_mass[hardestTrkLeading]);
	TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
	float mass = HardestTrackPlusMuon4.M();
	InvMassHardestNtrk23leadingH->Fill(mass,weight);
	InvMassHardestNtrk23H->Fill(mass,weight);
	TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(analysisTree.track_px[softestTrkLeading],
							    analysisTree.track_py[softestTrkLeading],
							    analysisTree.track_pz[softestTrkLeading],
							    analysisTree.track_mass[softestTrkLeading]);
	TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
	mass = SoftestTrackPlusMuon4.M();
	InvMassSoftestNtrk23leadingH->Fill(mass,weight);
	InvMassSoftestNtrk23H->Fill(mass,weight);
      }      
      if (isNtrk23trailing) {
	TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(analysisTree.track_px[hardestTrkTrailing],
							    analysisTree.track_py[hardestTrkTrailing],
							    analysisTree.track_pz[hardestTrkTrailing],
							    analysisTree.track_mass[hardestTrkTrailing]);
	TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
	float mass = HardestTrackPlusMuon4.M();
	InvMassHardestNtrk23trailingH->Fill(mass,weight);
	InvMassHardestNtrk23H->Fill(mass,weight);
	TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(analysisTree.track_px[softestTrkTrailing],
							    analysisTree.track_py[softestTrkTrailing],
							    analysisTree.track_pz[softestTrkTrailing],
							    analysisTree.track_mass[softestTrkTrailing]);
	TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
	mass = SoftestTrackPlusMuon4.M();
	InvMassSoftestNtrk23trailingH->Fill(mass,weight);
	InvMassSoftestNtrk23H->Fill(mass,weight);
      }      
      
      // sidebands Ntrk1
      if (isNtrk1leading) {
	TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(analysisTree.track_px[hardestTrkLeading],
							    analysisTree.track_py[hardestTrkLeading],
							    analysisTree.track_pz[hardestTrkLeading],
							    analysisTree.track_mass[hardestTrkLeading]);
	TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
	float mass = HardestTrackPlusMuon4.M();
	InvMassHardestNtrk1leadingH->Fill(mass,weight);
	InvMassHardestNtrk1H->Fill(mass,weight);
	TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(analysisTree.track_px[softestTrkLeading],
							    analysisTree.track_py[softestTrkLeading],
							    analysisTree.track_pz[softestTrkLeading],
							    analysisTree.track_mass[softestTrkLeading]);
	TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
	mass = SoftestTrackPlusMuon4.M();
	InvMassSoftestNtrk1leadingH->Fill(mass,weight);
	InvMassSoftestNtrk1H->Fill(mass,weight);
      }      
      if (isNtrk1trailing) {
	TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(analysisTree.track_px[hardestTrkTrailing],
							    analysisTree.track_py[hardestTrkTrailing],
							    analysisTree.track_pz[hardestTrkTrailing],
							    analysisTree.track_mass[hardestTrkTrailing]);
	TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
	float mass = HardestTrackPlusMuon4.M();
	InvMassHardestNtrk1trailingH->Fill(mass,weight);
	InvMassHardestNtrk1H->Fill(mass,weight);
	TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(analysisTree.track_px[softestTrkTrailing],
							    analysisTree.track_py[softestTrkTrailing],
							    analysisTree.track_pz[softestTrkTrailing],
							    analysisTree.track_mass[softestTrkTrailing]);
	TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
	mass = SoftestTrackPlusMuon4.M();
	InvMassSoftestNtrk1trailingH->Fill(mass,weight);
	InvMassSoftestNtrk1H->Fill(mass,weight);
      }      
      
      // signal region
      if (signalRegion) {
	counter_FinalEventsH->Fill(1.,weight);
	int iTrkLeading = trkSigTrailingMu[0];
	TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(analysisTree.track_px[iTrkLeading],
                                                            analysisTree.track_py[iTrkLeading],
                                                            analysisTree.track_pz[iTrkLeading],
                                                            analysisTree.track_mass[iTrkLeading]);
	TLorentzVector MuonTrackLeading4 = LeadingMuon4 + TrackLeading4;
	float massTrkMuLeading = MuonTrackLeading4.M();
	
	int iTrkTrailing = trkSigTrailingMu[0];
        TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(analysisTree.track_px[iTrkTrailing],
							      analysisTree.track_py[iTrkTrailing],
							      analysisTree.track_pz[iTrkTrailing],
							      analysisTree.track_mass[iTrkTrailing]);
        TLorentzVector MuonTrackTrailing4 = TrailingMuon4 + TrackTrailing4;
	float massTrkMuTrailing = MuonTrackTrailing4.M();
	
	InvMassLeadingH->Fill(massTrkMuLeading,weight);
	InvMassTrailingH->Fill(massTrkMuTrailing,weight);
	InvMassH->Fill(massTrkMuLeading,weight);
	InvMassH->Fill(massTrkMuTrailing,weight);
	InvMass2DH->Fill(massTrkMuLeading,massTrkMuTrailing,weight);
      }
      
      // Correlation distributions
      // ****
      // *** Amithabh, you should add conditions that there are
      // *** no additional tracks apart from signal track and soft tracks
      // *** I modified code below
      // ****
      // Old code
      //      bool LeadMuonControl_2Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1) && trkSigTrailingMu.size()==1;
      //      bool TrailMuonControl_2Tracks = (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1) && trkSigLeadingMu.size()==1;
      //      bool LeadMuonControl_3Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2) && trkSigTrailingMu.size()==1;
      //      bool TrailMuonControl_3Tracks = (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2) && trkSigLeadingMu.size()==1;
      //      bool BothMuonControl_2Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1) && (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1) ;
      //      bool BothMuonControl_3Tracks = (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2) && (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2) ;
      //      bool ControlAll = LeadMuonControl_2Tracks || TrailMuonControl_2Tracks || LeadMuonControl_3Tracks || TrailMuonControl_3Tracks || BothMuonControl_2Tracks || BothMuonControl_3Tracks;
      // New code
      bool signalLeadingMu = trkSigLeadingMu.size()==1 && trkLeadingMu.size()==1;

      bool signalTrailingMu = trkSigTrailingMu.size()==1 && trkTrailingMu.size()==1;

      bool bkgdLeadingMu = 
	(trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
	(trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3);

      bool bkgdTrailingMu = 
	(trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
	(trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3);

      bool ControlAll = (signalLeadingMu&&bkgdTrailingMu) || 
	(signalTrailingMu&&bkgdLeadingMu) || (bkgdLeadingMu&&bkgdTrailingMu);

      // * Now we can use this boolean to select bkg sideband
      // * where correlation coefficients are computed
      // * It is sufficient to use only one boolean - ControlAll
    if(ControlAll){
      // leading muon and associated track
      int iTrkLeading = trkSigLeadingMu[0];
      TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(analysisTree.track_px[iTrkLeading],
							  analysisTree.track_py[iTrkLeading],
							  analysisTree.track_pz[iTrkLeading],
							  analysisTree.track_mass[iTrkLeading]);
      TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
      float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();

      // trailing muon and associated track
      int iTrkTrailing = trkSigTrailingMu[0];
      TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(analysisTree.track_px[iTrkTrailing],
							    analysisTree.track_py[iTrkTrailing],
							    analysisTree.track_pz[iTrkTrailing],
							    analysisTree.track_mass[iTrkTrailing]);
      TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
      float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();

      float masshigh = massLeadingMuonTrk;
      float masslow = massTrailingMuonTrk;
      
      if (masshigh<masslow) {
	masshigh = massTrailingMuonTrk;
	masslow = massLeadingMuonTrk;
      }

      // filling histograms
      InvMassTrackPlusMuon1D_ControlH->Fill(massLeadingMuonTrk,weight);
      InvMassTrackPlusMuon1D_ControlH->Fill(massTrailingMuonTrk,weight);
      InvMassTrackPlusMuon2D_ControlH->Fill(masslow, masshigh, weight);
    }

    selEvents++;
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
  histRuns->Write();
  histWeightsH->Write();
  histTopPt->Write(); 
  histTopPtSq->Write(); 
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}
