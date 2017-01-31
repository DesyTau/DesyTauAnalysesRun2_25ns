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
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
using namespace std;

double TauiD(string sel,string working_point){

float SF = 1;

if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Loose" )
        SF = 0.81;

if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Medim" )
        SF = 0.79;
if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased_Tight" )
        SF = 0.79;


if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vloose" )
        SF = 0.83;
if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vloose" )
        SF= 0.84;
if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Medium" )
        SF = 0.84;

if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Tight" )
        SF = 0.83;

if ((sel =="mutau" || sel == "eltau")  && working_point == "MVA_Vtight" )
        SF = 0.80;



return SF;

}



double TauFakeRate(float pt,float eta, string sel,string working_point){

float SF = 1;

//80x MVAid

if ( working_point == "MVA"){



if ( sel =="mutau"){

if (  fabs(eta) < 0.9 )
        {
                if (pt>20 && pt<30) SF = 1.07968;
                if (pt>30 && pt<40) SF = 0.847355;
                if (pt>40 ) SF = 0.833833;
        }
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 )
        {

                if (pt>20 && pt<30) SF = 1.03273;
                if (pt>30 && pt<40) SF = 0.943387;
                if (pt>40 ) SF = 1.07113;
        }

if (  fabs(eta) > 1.2 && fabs(eta) < 2.1 )
        {

                if (pt>20 && pt<30) SF = 1.08635;
                if (pt>30 && pt<40) SF = 1.13754;
                if (pt>40) SF = 1.0387;
        }
if (  fabs(eta) > 2.1 && fabs(eta) < 2.4 )
        {

                if (pt>20 && pt<30) SF = 0.977868;
                if (pt>30 && pt<40) SF = 0.974665;
                if (pt>40) SF = 0.902185;
        }

}//mutau MVA

if ( sel =="eltau"){

if (  fabs(eta) < 0.8 )
        {
                if (pt>20 && pt<30) SF = 1.05966;
                if (pt>30 && pt<40) SF = 0.91867;
                if (pt>40 ) SF = 0.800276;
        }
if (  fabs(eta) > 0.8 && fabs(eta) < 1.44 )
        {

                if (pt>20 && pt<30) SF = 1.02955;
                if (pt>30 && pt<40) SF = 1.00145;
                if (pt>40 ) SF = 0.874349;
        }

if (  fabs(eta) > 1.44 && fabs(eta) < 1.566 )
	{


		if (pt>20 && pt<30) SF = 1.20397;
		if (pt>30 && pt<40) SF = 0.514959;
		if (pt>40) SF = 0.690121;
	}
if (  fabs(eta) > 1.566 && fabs(eta) < 2.1 )
	{


		if (pt>20 && pt<30) SF = 0.960575;
		if (pt>30 && pt<40) SF = 1.09704;
		if (pt>40) SF = 0.839154;
	}
}//eltau MVA

}








// Charged
if ( (sel =="mutau" || sel == "eltau") && working_point == "ChargedIso" ){

if (  fabs(eta) < 0.9 ) 
	{
		if (pt>20 && pt<30) SF = 1.26544;
		if (pt>30 && pt<50) SF = 1.25239;
		if (pt>50 ) SF = 1.38857;
	}
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 ) 
	{

		if (pt>20 && pt<30) SF = 1.21749;
		if (pt>30 && pt<50) SF = 1.10979;
		if (pt>50 ) SF = 1.60393;
	}

if (  fabs(eta) > 1.2 && fabs(eta) < 2.4 ) 
	{

		if (pt>20 && pt<30) SF = 1.27961;
		if (pt>30 && pt<50) SF = 1.14411;
		if (pt>50) SF = 1.2188;
	}
}


//CutBased
if ((sel =="mutau" || sel == "eltau")  && working_point == "CutBased" ){

if (  fabs(eta) < 0.9 ) 
	{
		if (pt>20 && pt<30) SF = 0.898437;
		if (pt>30 && pt<50) SF = 0.946704;
		if (pt>50 ) SF = 0.96842;
	}
if (  fabs(eta) > 0.9 && fabs(eta) < 1.2 ) 
	{

		if (pt>20 && pt<30) SF = 1.27757;
		if (pt>30 && pt<50) SF = 1.26811;
		if (pt>50 ) SF = 0.75345;
	}

if (  fabs(eta) > 1.2 && fabs(eta) < 2.4 ) 
	{

		if (pt>20 && pt<30) SF = 1.0773;
		if (pt>30 && pt<50) SF = 1.00049;
		if (pt>50) SF = 0.820108;
	}
}


return SF;

}


int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist



  // **** configuration
  Config cfg(argv[1]);
  string Channel="mutau";

  // kinematic cuts on electrons
  
  // kinematic cuts on electrons
  bool fillplots= false;
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

  RecoilCorrector recoilMetCorrector("HTT-utilities/RecoilCorrections/data/PFMET_MG_2016BCD_RooT_5.2.root");

  MEtSys metSys("HTT-utilities/RecoilCorrections/data/MEtSys.root");

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
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }
  TString MainTrigger(TrigLeg);

  const double bTag   = cfg.get<double>("bTag");

  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("No cut after PU");
  CutList.push_back("topPtReweighting");
  CutList.push_back("METFilters");
  CutList.push_back("Trigger");
  CutList.push_back("el");
  CutList.push_back("tau");
  CutList.push_back("el-tau OS");
  CutList.push_back("diEl Veto");
  CutList.push_back("3rd Lep Veto");
  CutList.push_back("Trigger eff");
  CutList.push_back("LSF");
  CutList.push_back("gt 0 jets");
  CutList.push_back("Loose/Tight");
  CutList.push_back("MET");
  CutList.push_back("MT");
  CutList.push_back("DeltaPhi");
  CutList.push_back("Ratio");

  int CutNumb = int(CutList.size());

  TH1D *CutFlowUnWLoose= new TH1D("CutFlowUnWLoose","Cut Flow",CutN,1,CutN+1);
  TH1D *CutFlowUnWTight= new TH1D("CutFlowUnWTight","Cut Flow",CutN,1,CutN+1);

  for(int cj = 0; cj < CutNumb; ++cj)    {
    CutFlowUnWLoose->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
    CutFlowUnWTight->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
  }

  xs=1;fact=1;fact2=1;
  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  ifstream ifs("xsecs");
  string line;

  XSec = 1;
  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;

  char ff[100];


  sprintf(ff,"%s/%s",argv[3],argv[2]);


  cout<<"  Initializing PU files....."<<endl;

// PU reweighting
  PileUp * PUofficial = new PileUp();
   //TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_2016_Cert_Cert_271036-276811_NoL1T_xsec63mb.root","read");

//  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_Cert_271036-277148_13TeV_PromptReco_Collisions16_xsec69p2mb.root","read");
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_RunBCDE_ReReco.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring16_PU25ns_V1.root", "read");


  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);



  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv");
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central");
  reader_B.load(calib,BTagEntry::FLAV_B,"comb");
  reader_C.load(calib,BTagEntry::FLAV_C,"comb");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");


  float etaBTAG[2] = {0.5,2.1};
  float ptBTAG[5] = {25.,35.,50.,100.,200.};

  std::cout << std::endl;
  for (int iEta=0; iEta<2; ++iEta) {
    for (int iPt=0; iPt<5; ++iPt) {
      float sfB = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
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
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016.root"); 
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016.root" << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get("zptmass_histo"); 
  if (histZMassPtWeights==NULL) {
    std::cout << " ZMassPT Weights histogram cannot found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << std::endl;
    exit(-1);
  }





  TH1D * ElSF_IdIso_El1H = new TH1D("ElIdIsoSF_El1H", "ElIdIsoSF_El1", 100, 0.5,1.5);

	cout<<" Initializing iD SF files....."<<endl;

  ScaleFactor * SF_electronIdIso; 
  if (applyLeptonSF) {
    SF_electronIdIso = new ScaleFactor();
    SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
  }

	cout<<" Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_electronTrigger = new ScaleFactor();
  SF_electronTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(SingleElectronTriggerFile));

  int nTotalFiles = 0;
  int nominator=-1;int denominator = -1;

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
  TFile * file ;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
  TH1::SetDefaultSumw2(true);
  //file->mkdir(Channel.c_str());
  //file->cd();
  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
  TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);
  TH1D * NumberOfVerticesH = new TH1D("NumberOfVerticesH","",51,-0.5,50.5);
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * histWeightsSkimmedH = new TH1D("histWeightsSkimmedH","",1,-0.5,0.5);
  TH1D *hnpv = new TH1D ("hnpv","npv ",100,-0.5,99.5);; 
  TH1D *hnpu = new TH1D ("hnpu","npu ",100,-0.5,99.5);;
  TH1D *hnpvS = new TH1D ("hnpvS","npv ",100,-0.5,99.5);; 
  TH1D *hnpuS = new TH1D ("hnpuS","npu ",100,-0.5,99.5);;

  TH1D * hMT = new TH1D("hMT","",20,0,200);
  TH1D * hRatioSum = new TH1D("hRatioSum","",10,0,1);
  TH1D * hDPhi = new TH1D("hDPhi","",70,0,3.5);
  TH1D * hMET = new TH1D("hMET","",10,0,200);
  TH1D * hnJets = new TH1D("hnJets","",15,-0.5,14.5);
  TH1D * hIsoEl = new TH1D("hIsoEl","",100,0,0.5);
  TH1D * hIsoElSel = new TH1D("hIsoElSel","",100,0,0.5);

  TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);
  TH1D * hMTCut1L = new TH1D("hMTCut1L","",20,0,200);
  TH1D * hMTCut2L = new TH1D("hMTCut2L","",20,0,200);
  TH1D * hMTCut3L = new TH1D("hMTCut3L","",20,0,200);
  TH1D * hMTCut4L = new TH1D("hMTCut4L","",20,0,200);
  TH1D * hMTCutTFRL = new TH1D("hMTCutTFRL","",20,0,200);
  TH1D * hMTCut1T = new TH1D("hMTCut1T","",20,0,200);
  TH1D * hMTCut2T = new TH1D("hMTCut2T","",20,0,200);
  TH1D * hMTCut3T = new TH1D("hMTCut3T","",20,0,200);
  TH1D * hMTCut4T = new TH1D("hMTCut4T","",20,0,200);
  TH1D * hMTCutTFRT = new TH1D("hMTCutTFRT","",20,0,200);

  TH1D * hRatioSum1L = new TH1D("hRatioSum1L","",10,0,1);
  TH1D * hRatioSum2L = new TH1D("hRatioSum2L","",10,0,1);
  TH1D * hRatioSum3L = new TH1D("hRatioSum3L","",10,0,1);
  TH1D * hRatioSum4L = new TH1D("hRatioSum4L","",10,0,1);
  TH1D * hRatioSumTFRL = new TH1D("hRatioSumTFRL","",10,0,1);
  TH1D * hRatioSum1T = new TH1D("hRatioSum1T","",10,0,1);
  TH1D * hRatioSum2T = new TH1D("hRatioSum2T","",10,0,1);
  TH1D * hRatioSum3T = new TH1D("hRatioSum3T","",10,0,1);
  TH1D * hRatioSum4T = new TH1D("hRatioSum4T","",10,0,1);
  TH1D * hRatioSumTFRT = new TH1D("hRatioSumTFRT","",10,0,1);

  TH1D * hDPhiCut1L = new TH1D("hDPhiCut1L","",70,0,3.5);
  TH1D * hDPhiCut2L = new TH1D("hDPhiCut2L","",70,0,3.5);
  TH1D * hDPhiCut3L = new TH1D("hDPhiCut3L","",70,0,3.5);
  TH1D * hDPhiCut4L = new TH1D("hDPhiCut4L","",70,0,3.5);
  TH1D * hDPhiCutTFRL = new TH1D("hDPhiCutTFRL","",70,0,3.5);
  TH1D * hDPhiCut1T = new TH1D("hDPhiCut1T","",70,0,3.5);
  TH1D * hDPhiCut2T = new TH1D("hDPhiCut2T","",70,0,3.5);
  TH1D * hDPhiCut3T = new TH1D("hDPhiCut3T","",70,0,3.5);
  TH1D * hDPhiCut4T = new TH1D("hDPhiCut4T","",70,0,3.5);
  TH1D * hDPhiCutTFRT = new TH1D("hDPhiCutTFRT","",70,0,3.5);
  TH1D * hMETCut1L = new TH1D("hMETCut1L","",10,0,200);
  TH1D * hMETCut2L = new TH1D("hMETCut2L","",10,0,200);
  TH1D * hMETCut3L = new TH1D("hMETCut3L","",10,0,200);
  TH1D * hMETCut4L = new TH1D("hMETCut4L","",10,0,200);
  TH1D * hMETCutTFRL = new TH1D("hMETCutTFRL","",10,0,200);
  TH1D * hMETCut1T = new TH1D("hMETCut1T","",10,0,200);
  TH1D * hMETCut2T = new TH1D("hMETCut2T","",10,0,200);
  TH1D * hMETCut3T = new TH1D("hMETCut3T","",10,0,200);
  TH1D * hMETCut4T = new TH1D("hMETCut4T","",10,0,200);
  TH1D * hMETCutTFRT = new TH1D("hMETCutTFRT","",10,0,200);

  TH1D * hLooseIndex = new TH1D("hLooseIndex","",5,0,5);
  TH1D * hTightIndex = new TH1D("hTightIndex","",5,0,5);

 const int nPtBins = 3;
  float ptBins[4] = {20,30,40,1000};
	

  TString PtBins[3] = {
    "Pt20to30",
    "Pt30to40",
    "PtGt40"};//,
    //"PtGt60"};//,		       "Pt100to150",		       "Pt150to200",		        "PtGt200"};

 const int nPtBins2 = 2;
  float ptBins2[3] = {20,40,1000};
	

  TString PtBins2[2] = {
    "Pt20to40",
    "PtGt40"};//,


const    int nEtaBins = 4;
    float etaBins[5] = {0, 0.8, 1.444, 1.566, 2.3}; 

    TString EtaBins[4] = {"EtaLt0p8",
    "Eta0p8to1p44",
    "Eta1p44to1p566",
    "EtaGt1p566"};
		


  const int nCuts = 4;

  TString Cuts[4] = {"MET","mT","DPhi","All"};
  /////first stands for the Eta bin, second array for the cut 
  TH1D * FakeRatePtIncLoose[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTight[nEtaBins][nCuts];



  TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
 // etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; ++i){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  
  TH1D * PtBinsH = new TH1D("PtBinsH", "PtBinsH", nPtBins, ptBins);
 // PtBinsH->Draw();
  PtBinsH->GetXaxis()->Set(nPtBins, ptBins);
  for (int i=0; i<nPtBins; ++i){ PtBinsH->GetXaxis()->SetBinLabel(i+1, PtBins[i]);}

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iCut=0; iCut<nCuts; ++iCut) {
    	if (iEta!=2){  
      		FakeRatePtIncLoose[iEta][iCut] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTight[iEta][iCut] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
	}
  
   	if (iEta==2) {
      		FakeRatePtIncLoose[iEta][iCut] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTight[iEta][iCut] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
    		}
    }


    //for (int iPt=0; iPt<nPtBins; ++iPt) {
    //  FakeRatePt[iEta][iPt] = new TH1D("FakeRatePt"+EtaBins[iEta]+PtBins[iPt],"",100,0,1000);
    //  FakeRateNV[iEta][iPt] = new TH1D("FakeRateNV"+EtaBins[iEta]+PtBins[iPt],"",50,0,50);
    //  FakeRateEta[iEta][iPt] = new TH1D("FakeRateEta"+EtaBins[iEta]+PtBins[iPt],"",80,-4,4);
    // }

  }

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int selEventsAllElectrons = 0;
  int iCutL=0;
  int iCutT=0;
  double LooseCFCounter[CutNumb];
  double TightCFCounter[CutNumb];
  double statUnc[CutNumb];
  int iCFCounter[CutNumb];
  for (int i=0;i < CutNumb; ++i){
  LooseCFCounter[i] = 0;
  TightCFCounter[i] = 0;
    iCFCounter[i] = 0;
    //statUnc[i] =0;
  }
  bool lumi=false;
  bool CutBasedTauId = true;

  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));

bool SUSY = false;
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
	/*		if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
	  histWeightsH->Fill(0.,genweights);
	}
    
      }



    bool isWJ = false;
    bool isTT = false;
    bool isDY = false;
    bool isDYhigh = false;
    bool isDYlow = false;
    bool isDYNJ = false;
    bool isWNJ = false;
    bool isNJ = false;
	string tt = "TT_TuneCUETP8M1_13TeV-powheg-pythia8";
	string wj = "WJetsToLNu";
	string wj1 = "W1JetsToLNu";
	string wj2 = "W2JetsToLNu";
	string wj3 = "W3JetsToLNu";
	string wj4 = "W4JetsToLNu";
	string dyj = "DYJetsToLL";
	string dyjhigh = "DYJetsToLL_M-50";
	string dyjlow = "DYJetsToLL_M-5to";
	string dyjlow2 = "DYJetsToLL_M-10to";
	string dyj1 = "DY1JetsToLL";
	string dyj2 = "DY2JetsToLL";
	string dyj3 = "DY3JetsToLL";
	string dyj4 = "DY4JetsToLL";
	if (string::npos != filen.find(tt)) isTT= true;
	if (string::npos != filen.find(wj)) isWJ= true;
	if (string::npos != filen.find(dyj)) isDY= true;
	if (string::npos != filen.find(dyjhigh)) isDYhigh= true;
	if (string::npos != filen.find(dyjlow) || string::npos != filen.find(dyjlow2)) isDYlow= true;
	if (string::npos != filen.find(dyj1) || string::npos != filen.find(dyj2) || string::npos != filen.find(dyj3) || string::npos != filen.find(dyj4)) isDYNJ = true;
	if (string::npos != filen.find(wj1) || string::npos != filen.find(wj2) || string::npos != filen.find(wj3) || string::npos != filen.find(wj4)) isWNJ= true;
	isNJ = isDYNJ && isWNJ;
	




    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      float topPt = 0;
      float antitopPt = 0;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;

      analysisTree.GetEntry(iEntry);
      nEvents++;

      iCutT = 0;
      iCutL = 0;

      int nparton = analysisTree.genparticles_noutgoing;
		
      if (isWJ && (nparton>0 && nparton<5)) continue;
      if ( isDYhigh  && (nparton>0 && nparton<5)) continue;

/////////needed for Recoil
      bool isW = false;
      bool isDY = false;
      bool isZTT = false;
      bool isZMM = false;
      bool isZEE = false;
      bool isTOP = false;
      if (!isData &&  string::npos != filen.find("JetsToLNu") ) isW=true;
      if (!isData &&  string::npos != filen.find("JetsToLL_M") ) isDY=true;
      if (!isData &&  string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8") ) isTOP=true;

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

	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

	  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);

	  if (analysisTree.genparticles_pdgid[igen]==6)
	    topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);

	  if (analysisTree.genparticles_pdgid[igen]==-6)
	    antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				    analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);

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
	if (isDY) {
	  
	  if (promptTausFirstCopy.size()==2) {
	    isZTT = true; isZMM = false; isZEE = false;
	    bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz(); 
	    bosonMass = promptTausLV.M();
	    bosonEta  = promptTausLV.Eta();
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	    //mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);
	  }
	  else if (promptMuons.size()==2) {
	    isZTT = false; isZMM = true; isZEE = false;
	    bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz(); 
	    bosonMass = promptMuonsLV.M(); 
	    bosonEta = promptMuonsLV.Eta();
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	    //mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
	  }
	  else {
	    isZTT = false; isZMM = false; isZEE = true;
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



      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  


      bool lumi=false;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

      if (!isData )  { 

	weight *= analysisTree.genweight;
	gen_weight *=analysisTree.genweight;
	lumi=true;
	//cout<<"  weight from init "<<genweights<< "  "<<analysisTree.genweight<<"  "<<weight<<endl;


	if (applyPUreweighting)	 {
	  puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	  weight *=puweight; 
	  //pu_weight = puweight;
	}

      }

      LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

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

	  if (topPt>0.&& antitopPt>0. && !isData) {
	    topptweight = topPtWeight(topPt,antitopPt);
	    weight *= topptweight;
	    top_weight = topptweight;
	     // cout<<"  "<<topPt<<"  "<<antitopPt<<"  "<<topptweight<<endl;
	  }


      histTopPt->Fill(0.,topptweight);

	}

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

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
	//    	if (!lumi) cout<< " Failed to find good lumi "<<endl;
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


      //lumi=true;
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



	bool METflag = metFiltersPasses2(analysisTree, metFlags);
	if (!METflag && isData) continue;


      LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;



      JetsMV.clear();
      ElMV.clear();
      TauMV.clear();
      MuMV.clear();
      LeptMV.clear();
      mu_index=-1;
      tau_index=-1;
      el_index=-1;


      bool trigAccept = false;

      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;


      if (isData){
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
      }

	if (!isData) isMainTrigger = true;

      if (!isMainTrigger) {
	std::cout << "Fail on main HLT Filter " << MainTrigger << " not found" << std::endl;
	return(-1);
      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
      /////now clear the Mu.El.Jets again to fill them again after cleaning


      vector<int> electrons; electrons.clear();
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {
	if (analysisTree.electron_pt[im]<ptElectronCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut) continue;
	if (applyElectronId && !analysisTree.electron_mva_wp80_nontrig_Spring15_v1[im]) continue;
	if (applyElectronId && !analysisTree.electron_pass_conversion[im]) continue;
	if (applyElectronId && analysisTree.electron_nmissinginnerhits[im]>1) continue;
	if (fabs(analysisTree.electron_charge[im]) !=1) continue;
	 electrons.push_back(im);

      }
      if (electrons.size()==0 ) continue;


      LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (analysisTree.tau_decayModeFinding[it]<decayModeFinding) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
        if ( fabs(analysisTree.tau_charge[it]) != 1 ) continue;
	  taus.push_back((int)it);

	}

      if (taus.size()!=1)  continue;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

      int tau_index = -1;
      int mu_index = -1;

      float isoElecMin  = 1e+10;
      float isoTauMin = 1; 
      if (CutBasedTauId) isoTauMin = 1e+10;
      if (!CutBasedTauId) isoTauMin = -10;
      float ptMu = 0;
      float ptTau = 0;
      //      if (electrons.size()>1||electrons.size()>1)
      //      std::cout << "electrons = " << electrons.size() << "  taus = " << taus.size() << std::endl;
      for (unsigned int im=0; im<electrons.size(); ++im) {
	bool isLegMatch = false;
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
	float neutralIsoElec = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	neutralIsoElec = TMath::Max(float(0),neutralIsoElec); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/analysisTree.electron_pt[eIndex];


  	hIsoEl->Fill(relIsoElec,weight);

	if (relIsoElec > 0.1) continue;

	if (isData)
	{
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      &&analysisTree.electron_pt[eIndex]>ptElectronCut&&
	      analysisTree.trigobject_pt[iT]>SingleElectronTriggerPtCut) { // IsoElec Leg
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) {
	      isLegMatch = true;
	    }
	  }
	 }
	}
      if (!isData) isLegMatch = true;

      if (!isLegMatch) continue;

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
	      if (analysisTree.electron_pt[eIndex]>ptMu) {
		isoElecMin  = relIsoElec;
		ptMu = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }
	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptMu = analysisTree.electron_pt[eIndex];
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

	  if ( (int) eIndex!= (int)el_index) {
	    if (relIsoElec==isoElecMin) {
	      if (analysisTree.electron_pt[eIndex]>ptMu) {
		isoElecMin  = relIsoElec;
		ptMu = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }

	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptMu = analysisTree.electron_pt[eIndex];
	      el_index = (int)eIndex;
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tau_index = (int)tIndex;
	    }
          }/*
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
          }*/

        }

      }
      }
      //      std::cout << "eIndex = " << el_index << "   tau_index = " << tau_index << std::endl;
      //
      bool TauId = false;

            //if ( analysisTree.tau_againstElectronVLooseMVA6[tau_index]>0.5 &&   analysisTree.tau_againstMuonTight3[tau_index]>0.5) TauId = true;
            if ( analysisTree.tau_againstElectronTightMVA6[tau_index]>0.5 &&   analysisTree.tau_againstMuonLoose3[tau_index]>0.5) TauId = true;

	if (!TauId) continue;

      if ((int)tau_index<0) continue;
      if ((int)el_index<0) continue;



      //      std::cout << "Ok4 " << std::endl;
      bool isLoose = false;
      unsigned int tau_loose=-1;
      unsigned int tau_tight=-1;
      vector<int> tau; tau.clear();

 	if ((int)tau_index>-1) isLoose  = true;
	tau_loose = (int)tau_index;
	
	if (!isLoose) continue;
	 	 

      double q = analysisTree.tau_charge[tau_loose] * analysisTree.electron_charge[(int)el_index];
      if (q > 0) continue;

	double dReltau = deltaR(analysisTree.tau_eta[tau_loose],analysisTree.tau_phi[tau_loose],
				analysisTree.electron_eta[(int)el_index],analysisTree.electron_phi[(int)el_index]);
	if (dReltau < 0.5) continue;

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;


  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;
      

      // looking for extra muon
      bool foundExtraMuon = false;
      for (unsigned int ie = 0; ie<analysisTree.muon_count; ++ie) {
	if (analysisTree.muon_pt[ie]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[ie])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[ie])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[ie])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[ie]) continue;
	//if (applyVetoMuonId && !analysisTree.muon_isICHEP[ie]) continue;
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
        double neutralIsoMu = TMath::Max(double(0),neutralIsoMuN);
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
	double neutralIsoElec = TMath::Max(double(0),neutralIsoElecN); 
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
	if (dilepton_veto)  continue;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

      if(extraelec_veto || extramuon_veto)   continue;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;




///////////////Trigger weight 
      double ptEl1 = (double)analysisTree.electron_pt[el_index];
      double etaEl1 = (double)analysisTree.electron_eta[el_index];
      float trigweight = 1.;
	//cout<<" this is what goes for trigger  "<<ptEl1<<"  "<<etaEl1<<endl;
      float EffFromData = (float)SF_electronTrigger->get_EfficiencyData(double(ptEl1),double(etaEl1));
      /*float Mu22EffMC   = (float)SF_muonTrigger->get_EfficiencyMC(double(ptEl1),double(etaEl1));*/
	

	//if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) Signal=true;
     /* if (!isData) {
	if (Mu22EffMC>1e-6)
	  trigweight = EffFromData / Mu22EffMC;
	if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) trigweight = EffFromData;
	weight *= trigweight;
	trig_weight = trigweight;
	//	cout<<" Trigger weight "<<trigweight<<endl;
      }*/
	if (!isData) trigweight = EffFromData;
	weight *= trigweight;

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

	///LSF 
      if (!isData && applyLeptonSF) {

	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
	double ptEl1 = (double)analysisTree.electron_pt[el_index];
	double etaEl1 = (double)analysisTree.electron_eta[el_index];
	double IdIsoSF_el1 = SF_electronIdIso->get_ScaleFactor(ptEl1, etaEl1);
	weight *= IdIsoSF_el1;
      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
 ////////////////////////////// now check if it is also Tight
 	if (isLoose) denominator++;

	bool isTight =false;float isoTau = 0;

       if (!CutBasedTauId){
		isTight=
	  	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index] > 0.5;
	  	  //analysisTree.tau_chargedIsoPtSum[tau_index] < 0.8;

 
       isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
       //ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];
       //isoTau = analysisTree.tau_chargedIsoPtSum[tau_index];

	 }

	if (CutBasedTauId){
		isTight=
	         analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] > 0.5;

          isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
	}


	tau_tight = (int)tau_index;

      JetsMV.clear();
      int countjets=0;
      int countbjets=0;
      float DRmax = 0.5;


      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);

	if (absJetEta > etaJetCut) continue;
	if (fabs(analysisTree.pfjet_pt[jet])<20.) continue;
	float jetPt = analysisTree.pfjet_pt[jet];


	bool isPFJetId = false ; 
	isPFJetId =looseJetiD(analysisTree,jet);
	//isPFJetId =tightLepVetoJetiD(analysisTree,jet);

	//				cout<<" 1- jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  iEntry "<<iEntry<<endl;
	if (!isPFJetId) continue;
	jet_isLoose[jet] = isPFJetId;
	bool cleanedJet = true;
	//				cout<<"  jet is Loose "<<isPFJetId<<"  "<<jet_isLoose[jet]<<"  "<<iEntry<<endl;
	double Dr=deltaR(analysisTree.electron_eta[(int)el_index],analysisTree.electron_phi[(int)el_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if (  Dr  < DRmax)  cleanedJet=false;


	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if ( Drr > 0.2) cleanedJet=false;

	if (!cleanedJet) continue;

	JetsV.SetPxPyPzE(0.,0.,0.,0.);
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
	JetsMV.push_back(JetsV);	

	countjets++;

      }

      if (countjets==0) continue;

      LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
      hnJets->Fill(countjets,weight);
      hIsoElSel->Fill(isoElecMin,weight);

      if (fabs(analysisTree.tau_charge[tau_loose]) !=1) continue;
      if (fabs(analysisTree.electron_charge[el_index]) !=1) continue;



      if (isLoose) { hLooseIndex->Fill((int)tau_loose,weight);
	}
      if (isTight) { hTightIndex->Fill((int)tau_tight,weight);
	}

      //	if ((int)tau_tight != (int)tau_loose) continue; ////require that the tau that is identified a loose is the same that is the tight one


/////////////////// Recoil corrections

      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;

      float pfmet_corr_x = analysisTree.pfmetcorr_ex;
      float pfmet_corr_y = analysisTree.pfmetcorr_ey;
      float met_x = analysisTree.pfmetcorr_ex;
      float met_y = analysisTree.pfmetcorr_ey;

      if ((isW||isDY) && !isData) {

	  recoilMetCorrector.CorrectByMeanResolution(analysisTree.pfmetcorr_ex,analysisTree.pfmetcorr_ey,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
 
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

      }//if isW, isDY !isData

      met_ex = pfmet_corr_x;
      met_ey = pfmet_corr_y;
      met_ez = 0;//analysisTree.pfmet_ez;
      //met_pt = analysisTree.pfmet_pt;
      met_pt = TMath::Sqrt(pfmet_corr_x*pfmet_corr_x+pfmet_corr_y*pfmet_corr_y);
      //met_phi = analysisTree.pfmet_phi;
      met_phi = TMath::ATan2(pfmet_corr_y,pfmet_corr_x);







      double dPhi=-1;double MT=-1 ; double RatioSums=-1;


      double met = sqrt ( met_ex*met_ex + met_ey*met_ey);
      // w = mu+MET
      // ptW - ptJ/ptW+ptJ      
      double SumPtMuMET = sqrt ( ( met_ex +analysisTree.electron_px[el_index])*(met_ex+analysisTree.electron_px[el_index]) + 
				 (met_ey+analysisTree.electron_py[el_index])*(met_ey+analysisTree.electron_py[el_index]));

      RatioSums = (SumPtMuMET - JetsMV.at(0).Pt() )/ (SumPtMuMET +JetsMV.at(0).Pt() );


      TLorentzVector MetV; 
      MetV.SetPx(met_ex); 
      MetV.SetPy(met_ey);


      TLorentzVector elV ;  elV.SetPtEtaPhiM(analysisTree.electron_pt[el_index], analysisTree.electron_eta[el_index], analysisTree.electron_phi[el_index], electronMass);
      TLorentzVector tauV;  tauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_loose], analysisTree.tau_eta[tau_loose], analysisTree.tau_phi[tau_loose], tauMass);

      TLorentzVector Wb = elV  + MetV;

      dPhi = dPhiFrom2P( elV.Px(), elV.Py(), MetV.Px(),  MetV.Py() );
      MT = TMath::Sqrt(2*elV.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));

			
      float dPhiW = dPhiFrom2P( Wb.Px(), Wb.Py(), JetsMV.at(0).Px(),  JetsMV.at(0).Py() );

      hRatioSum->Fill(RatioSums,weight);
      hMT->Fill(MT,weight);
      hDPhi->Fill(dPhiW, weight);
      hMET->Fill(met, weight);
      hnpv->Fill(analysisTree.primvertex_count,weight);
      hnpu->Fill(analysisTree.numtruepileupinteractions,weight);



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

	double tauId = TauiD("mutau","MVA_Tight");
		
	if (isTauMatched && !isGenLeptonMatched) weight *= tauId;


      //cout<<"  "<<endl;
      if (isTight) nominator++;

      float ptProbe = TMath::Min(float(analysisTree.tau_pt[tau_loose]),float(ptBins[nPtBins]-0.1));
      float absEtaProbe = fabs(analysisTree.tau_eta[tau_loose]);
      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
      if (ptBin<0) continue;
      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
      if (etaBin<0) continue;

      //cout<< "filling here  "<<analysisTree.tau_pt[tau_loose]<<"  "<<ptBin<<"  "<<etaBin<<"  "<<weight<<endl;
      //FakeRatePt[etaBin][ptBin]->Fill(analysisTree.tau_pt[tau_loose],weight);

      bool MTb = (MT>60 && MT<120);
      double ptTau1 = (double)analysisTree.tau_pt[tau_loose];
      double etaTau1 = (double)analysisTree.tau_eta[tau_loose];
      string WP = "MVA";
      double tfr  = TauFakeRate(ptTau1,etaTau1,Channel,WP);
      if (isTauMatched || isData) tfr=1.;
      if (isLoose)

      {
      LooseCFCounter[iCutL]+= weight;
      iCutL++;
	      if (met>40){

	    FakeRatePtIncLoose[etaBin][0]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum1L->Fill(RatioSums,weight);
	    hMTCut1L->Fill(MT,weight);
	    hDPhiCut1L->Fill(dPhiW, weight);
	    hMETCut1L->Fill(met, weight);
      LooseCFCounter[iCutL]+= weight;
      iCutL++;

	  if (MTb ){ 
	    FakeRatePtIncLoose[etaBin][1]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum2L->Fill(RatioSums,weight);
	    hMTCut2L->Fill(MT,weight);
	    hDPhiCut2L->Fill(dPhiW, weight);
	    hMETCut2L->Fill(met, weight);

      LooseCFCounter[iCutL]+= weight;
      iCutL++;

	    if (dPhiW > 2.5){
	    FakeRatePtIncLoose[etaBin][2]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum3L->Fill(RatioSums,weight);
	    hMTCut3L->Fill(MT,weight);
	    hDPhiCut3L->Fill(dPhiW, weight);
	    hMETCut3L->Fill(met, weight);
      LooseCFCounter[iCutL]+= weight;
      iCutL++;

	if (RatioSums < 0.3 ) {
	    FakeRatePtIncLoose[etaBin][3]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum4L->Fill(RatioSums,weight);
	    hMTCut4L->Fill(MT,weight);
	    hDPhiCut4L->Fill(dPhiW, weight);
	    hMETCut4L->Fill(met, weight);

      LooseCFCounter[iCutL]+= weight;
      iCutL++;
/////////// corrected for TFR
//
	    hRatioSumTFRL->Fill(RatioSums,tfr*weight);
	    hMTCutTFRL->Fill(MT,tfr*weight);
	    hDPhiCutTFRL->Fill(dPhiW, tfr*weight);
	    hMETCutTFRL->Fill(met, tfr*weight);

	      }//met<80

	    }//dPhiW	
				
	  }//MTb
	}//ratio

      }//Loose 

      if (isTight) 

	{	
      TightCFCounter[iCutT]+= weight;
      iCutT++;

		if (met>40){
	    FakeRatePtIncTight[etaBin][0]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum1T->Fill(RatioSums,weight);
	    hMTCut1T->Fill(MT,weight);
	    hDPhiCut1T->Fill(dPhiW, weight);
	    hMETCut1T->Fill(met, weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;

	    if (MTb ){ 
	    FakeRatePtIncTight[etaBin][1]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum2T->Fill(RatioSums,weight);
	    hMTCut2T->Fill(MT,weight);
	    hDPhiCut2T->Fill(dPhiW, weight);
	    hMETCut2T->Fill(met, weight);
    
      TightCFCounter[iCutT]+= weight;
      iCutT++;

	      if (dPhiW > 2.5){
	    FakeRatePtIncTight[etaBin][2]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum3T->Fill(RatioSums,weight);
	    hMTCut3T->Fill(MT,weight);
	    hDPhiCut3T->Fill(dPhiW, weight);
	    hMETCut3T->Fill(met, weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;

	  if (RatioSums < 0.3 ) {

	    FakeRatePtIncTight[etaBin][3]->Fill(analysisTree.tau_pt[tau_loose],weight);
	    hRatioSum4T->Fill(RatioSums,weight);
	    hMTCut4T->Fill(MT,weight);
	    hDPhiCut4T->Fill(dPhiW, weight);
	    hMETCut4T->Fill(met, weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;
      
	    hRatioSumTFRT->Fill(RatioSums,tfr*weight);
	    hMTCutTFRT->Fill(MT,tfr*weight);
	    hDPhiCutTFRT->Fill(dPhiW, tfr*weight);
	    hMETCutTFRT->Fill(met, tfr*weight);


		}//met<80

	      }//dPhiW	
				
	    }//MTb
	  }//ratio

	}//Tight



      selEvents++;
    } // end of file processing (loop over events in one file)

    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<"done"<<endl;
  double fr =  double(nominator/denominator) ;
  if (denominator>0) cout <<" the tau fake rate "<<fr<<endl;
  cout <<" the tau fake rate "<<nominator<<"  "<<denominator<<endl;

  for(int ci = 0; ci < CutNumb; ++ci)
    {
      CutFlowUnWLoose->SetBinContent(1+ci,0);
      CutFlowUnWTight->SetBinContent(1+ci,0);
      CutFlowUnWLoose->SetBinContent(1+ci,float(LooseCFCounter[ci]) );
      CutFlowUnWTight->SetBinContent(1+ci,float(TightCFCounter[ci]) );

    }

  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  file->cd();
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histTopPt->Write();
  hnJets->Write();
  CutFlowUnWTight->Write();
  CutFlowUnWLoose->Write();
  file->Write();
  file->Close();

  delete file;

}


