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
#include "TMath.h"
#include "Riostream.h"
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


if (  fabs(eta) < 0.8 )
        {
                if (pt>20 && pt<30) SF = 1.12667;
                if (pt>30 && pt<40) SF = 1.0178;
                if (pt>40 ) SF = 0.961808;
        }
if (  fabs(eta) > 0.8 && fabs(eta) < 1.44 )
        {

                if (pt>20 && pt<30) SF = 1.12651;
                if (pt>30 && pt<40) SF = 1.08668;
                if (pt>40 ) SF = 0.832438;
        }

if (  fabs(eta) > 1.44 && fabs(eta) < 1.566 )
        {

                if (pt>20 && pt<40) SF = 0.907089;
                if (pt>40) SF = 1.70447;
        }
if (  fabs(eta) > 1.566 && fabs(eta) < 2.3 )
        {

                if (pt>20 && pt<40) SF = 0.875675;
                if (pt>40) SF = 0.91729;
        }

	
}
return SF;

}


int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist



  // **** configuration
  Config cfg(argv[1]);
  const bool isData = cfg.get<bool>("IsData");

  string Channel="mutau";

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


  // kinematic cuts on electrons
  

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

  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("topPt");
  CutList.push_back("ZPt");
  CutList.push_back("No cut after PU");
  CutList.push_back("METFilters");
  CutList.push_back("Trigger");
  CutList.push_back("mu");
  CutList.push_back("tau");
  CutList.push_back("$mu-tau OS");
  CutList.push_back("Trigger eff");
  CutList.push_back("LSF");
  CutList.push_back("gt 0 jets");
  CutList.push_back("Loose/Tight");
  CutList.push_back("MET");
  CutList.push_back("MT");
  CutList.push_back("DeltaPhi");
  CutList.push_back("nJets");

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

  PileUp * PUofficial = new PileUp();
  

  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2016_271036-284044_80bins.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Moriond17_PU25ns_V1.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  cout<<" done with PU...."<<endl;
  
  string BtagCVS = "CSVv2Moriond17_2017_1_26_BtoH.csv" ;  

  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/"+BtagCVS);

  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,BTag_);
  reader_B.load(calib,BTagEntry::FLAV_B,"comb");
  reader_C.load(calib,BTagEntry::FLAV_C,"comb");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");


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


  cout<<"  Initializing iD SF files....."<<endl;
  ScaleFactor * SF_muonIdIso = new ScaleFactor(); 
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));


  cout<<"  Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuontrigEffFile));


  int nTotalFiles = 0;
  int nominator=-1;int denominator = -1;

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  std::ifstream fileList(ff);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");


  string SaveDir=argv[3];


  if (string::npos == Systematic.find("Nominal")) {SaveDir.append("_");SaveDir.append(argv[5]);}

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<"  The systematic will be "<<Systematic<<"  and save dir will be  "<<SaveDir<<endl;

  // output fileName with histograms
  
  TFile * file;
  if (isData) file = new TFile(SaveDir+"/"+TStrName+TString("_DataDriven.root"),"update");
  else file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"update");

  TH1::SetDefaultSumw2(true);

  //file->mkdir(Channel.c_str());
  //file->cd();
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D *hnpv = new TH1D ("hnpv","npv ",100,-0.5,99.5);; 
  TH1D *hnpu = new TH1D ("hnpu","npu ",100,-0.5,99.5);;
  TH1D *hnpvS = new TH1D ("hnpvS","npv ",100,-0.5,99.5);; 
  TH1D *hnpuS = new TH1D ("hnpuS","npu ",100,-0.5,99.5);;

  TH1D * hMT = new TH1D("hMT","",20,0,200);
  TH1D * hRatioSum = new TH1D("hRatioSum","",10,0,1);
  TH1D * hDPhi = new TH1D("hDPhi","",70,0,3.5);
  TH1D * hMET = new TH1D("hMET","",10,0,200);
  TH1D * hnJets = new TH1D("hnJets","",15,-0.5,14.5);
  TH1D * hnMatchedJets = new TH1D("hnMatchedJets","",15,-0.5,14.5);
  TH1D * hnMatchedJetsT = new TH1D("hnMatchedJetsT","",15,-0.5,14.5);
  TH1D * hnbJets = new TH1D("hnbJets","",15,-0.5,14.5);
  TH1D * hIsoMu = new TH1D("hIsoMu","",100,0,0.5);
  TH1D * hIsoMuSel = new TH1D("hIsoMuSel","",100,0,0.5);

  TH1D * hxsec = new TH1D("xsec","",1,0,10e+20);

  TH1D * hnJets1L = new TH1D("hnJets1L","",15,-0.5,14.5);
  TH1D * hnbJets1L = new TH1D("hnbJets1L","",15,-0.5,14.5);
  TH1D * hnJets2L = new TH1D("hnJets2L","",15,-0.5,14.5);
  TH1D * hnbJets2L = new TH1D("hnbJets2L","",15,-0.5,14.5);
  TH1D * hnJets3L = new TH1D("hnJets3L","",15,-0.5,14.5);
  TH1D * hnbJets3L = new TH1D("hnbJets3L","",15,-0.5,14.5);
  TH1D * hnJets4L = new TH1D("hnJets4L","",15,-0.5,14.5);
  TH1D * hnbJets4L = new TH1D("hnbJets4L","",15,-0.5,14.5);
  TH1D * hnJets1T = new TH1D("hnJets1T","",15,-0.5,14.5);
  TH1D * hnbJets1T = new TH1D("hnbJets1T","",15,-0.5,14.5);
  TH1D * hnJets2T = new TH1D("hnJets2T","",15,-0.5,14.5);
  TH1D * hnbJets2T = new TH1D("hnbJets2T","",15,-0.5,14.5);
  TH1D * hnJets3T = new TH1D("hnJets3T","",15,-0.5,14.5);
  TH1D * hnbJets3T = new TH1D("hnbJets3T","",15,-0.5,14.5);
  TH1D * hnJets4T = new TH1D("hnJets4T","",15,-0.5,14.5);
  TH1D * hnbJets4T = new TH1D("hnbJets4T","",15,-0.5,14.5);
  TH1D * hnJetsTFRL = new TH1D("hnJetsTFRL","",15,-0.5,14.5);
  TH1D * hnbJetsTFRL = new TH1D("hnbJetsTFRL","",15,-0.5,14.5);
  TH1D * hnJetsTFRT = new TH1D("hnJetsTFRT","",15,-0.5,14.5);
  TH1D * hnbJetsTFRT = new TH1D("hnbJetsTFRT","",15,-0.5,14.5);

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

  TH1D * FakeRatePtIncLooseUpQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightUpQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseDownQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightDownQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseCharmQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightCharmQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseStrangeQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightStrangeQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseGluon[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightGluon[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseBottomQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightBottomQ[nEtaBins][nCuts];
  TH1D * FakeRatePtIncLooseNothing[nEtaBins][nCuts];
  TH1D * FakeRatePtIncTightNothing[nEtaBins][nCuts];



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
    	if (iEta<2){  
      		FakeRatePtIncLoose[iEta][iCut] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTight[iEta][iCut] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);

      		FakeRatePtIncLooseUpQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseUpQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightUpQ[iEta][iCut] = new TH1D("FakeRatePtIncTightUpQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseDownQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseDownQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightDownQ[iEta][iCut] = new TH1D("FakeRatePtIncTightDownQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseCharmQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseCharmQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightCharmQ[iEta][iCut] = new TH1D("FakeRatePtIncTightCharmQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseStrangeQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseStrangeQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightStrangeQ[iEta][iCut] = new TH1D("FakeRatePtIncTightStrangeQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseBottomQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseBottomQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightBottomQ[iEta][iCut] = new TH1D("FakeRatePtIncTightBottomQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseGluon[iEta][iCut] = new TH1D("FakeRatePtIncLooseGluon"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightGluon[iEta][iCut] = new TH1D("FakeRatePtIncTightGluon"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncLooseNothing[iEta][iCut] = new TH1D("FakeRatePtIncLooseNothing"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
      		FakeRatePtIncTightNothing[iEta][iCut] = new TH1D("FakeRatePtIncTightNothing"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);

	}
  
   	if (iEta>1) {
      		FakeRatePtIncLoose[iEta][iCut] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTight[iEta][iCut] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);

      		FakeRatePtIncLooseUpQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseUpQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightUpQ[iEta][iCut] = new TH1D("FakeRatePtIncTightUpQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseDownQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseDownQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightDownQ[iEta][iCut] = new TH1D("FakeRatePtIncTightDownQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseCharmQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseCharmQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightCharmQ[iEta][iCut] = new TH1D("FakeRatePtIncTightCharmQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseStrangeQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseStrangeQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightStrangeQ[iEta][iCut] = new TH1D("FakeRatePtIncTightStrangeQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseBottomQ[iEta][iCut] = new TH1D("FakeRatePtIncLooseBottomQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightBottomQ[iEta][iCut] = new TH1D("FakeRatePtIncTightBottomQ"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseGluon[iEta][iCut] = new TH1D("FakeRatePtIncLooseGluon"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightGluon[iEta][iCut] = new TH1D("FakeRatePtIncTightGluon"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncLooseNothing[iEta][iCut] = new TH1D("FakeRatePtIncLooseNothing"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
      		FakeRatePtIncTightNothing[iEta][iCut] = new TH1D("FakeRatePtIncTightNothing"+EtaBins[iEta]+Cuts[iCut],"",nPtBins2,ptBins2);
    		}
    }

    //for (int iPt=0; iPt<nPtBins; ++iPt) {
    //  FakeRatePt[iEta][iPt] = new TH1D("FakeRatePt"+EtaBins[iEta]+PtBins[iPt],"",100,0,1000);
    //  FakeRateNV[iEta][iPt] = new TH1D("FakeRateNV"+EtaBins[iEta]+PtBins[iPt],"",50,0,50);
    //  FakeRateEta[iEta][iPt] = new TH1D("FakeRateEta"+EtaBins[iEta]+PtBins[iPt],"",80,-4,4);
    // }

  }

std::string st1,st2;
bool SUSY = false;
float SusyMotherMassF;
float SusyLSPMassF;

//SMS-TChiSlepSnu_x0p5_TuneCUETP8M1_13TeV-madgraphMLM-pythia8   SMS-TChiStauStau_x0p5_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  SMS-TStauStau_TuneCUETP8M1_13TeV-madgraphMLM-pythia8

if (string::npos != rootFileName.find("SMS-") || string::npos != rootFileName.find("stau") || string::npos != rootFileName.find("C1"))
	{
	//st1 =  rootFileName.substr(4,3);
	////SusyMotherMassF = stof(argv[5]);
	//st1=string(argv[5]);
	//st2 =  rootFileName.substr(11);
	//st2=string(argv[6]);
	//SusyLSPMassF = stof(argv[6]);
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



  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

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
  bool CutBasedTauId = false;

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
	genTauMatched = false;
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
	matchedTauToNothing = false;
	matchedTauToUpQ = false;
	matchedTauToDownQ = false;
	matchedTauToStrangeQ = false;
	matchedTauToCharmQ = false;
	matchedTauToBottom = false;
	matchedTauToHFQ = false;
	matchedTauToLFQ = false;
	
	string tt = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8";
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
	string ttw = "TTWJetsToLNu";
	string signal = "stau";
	string signalC1 = "C1";
	string signalChi = "Chi";
	if (string::npos != filen.find(tt)) isTT= true;
	if (string::npos != filen.find(wj)) isWJ= true;
	if (string::npos != filen.find(ttw)) isWJ= false;
	if (string::npos != filen.find(dyj)) isDY= true;
	if (string::npos != filen.find(dyjhigh)) isDYhigh= true;
	if (string::npos != filen.find(dyjlow) || string::npos != filen.find(dyjlow2)) isDYlow= true;
	if (string::npos != filen.find(dyj1) || string::npos != filen.find(dyj2) || string::npos != filen.find(dyj3) || string::npos != filen.find(dyj4)) isDYNJ = true;
	if ((string::npos != filen.find(wj1) || string::npos != filen.find(wj2) || string::npos != filen.find(wj3) || string::npos != filen.find(wj4)) && string::npos == filen.find(ttw)) isWNJ= true;
	isNJ = isDYNJ && isWNJ;


    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);

      iCutT = 0;
      iCutL = 0;


      int nparton = analysisTree.genparticles_noutgoing;
		
      if ( (isWJ  || isDYhigh) && (nparton>0 && nparton<5)) continue;

      if (iEntry%50000==0) 
	cout << "      processed " << iEntry << " events" << endl; 

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


      if (!isData && ( string::npos != filen.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8")  || string::npos != filen.find("TTPowHeg") || string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8")) ) 
	{
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

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
	  weight *=zptmassweight;
	}

      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

      if (!isData )  { 

	weight *= analysisTree.genweight;
	gen_weight *=analysisTree.genweight;
	lumi=true;

	  puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
	  weight *=puweight; 

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

      if (!lumi) continue;

	int RunNo = analysisTree.event_run;


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

  //    }


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
	std::cout << "Fail on main HLT Filter " << MainTrigger << " not found" << std::endl;
	return(-1);
      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

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
      	if ( analysisTree.tau_againstElectronVLooseMVA6[it]<0.5) continue; 
      	if ( analysisTree.tau_againstMuonTight3[it]<0.5) continue;

	  taus.push_back((int)it);

	}

      //if (taus.size()!=1)  continue;
      if (taus.size()==0)  continue;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

      int tau_index = -1;
      int mu_index = -1;

      float isoMuMin  = 1e+10;
      float ptMu = 0;
      float ptTau = 0;
      float isoTauMin = 1; 
      float isoTau = 1.; 
      if (CutBasedTauId) isoTauMin = 1e+10;
      if (!CutBasedTauId) isoTauMin = -10;
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
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];

  	hIsoMu->Fill(relIsoMu,weight);

	if (relIsoMu > 0.15) continue;

	if (!SUSY)
	{	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  	if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      	&& analysisTree.muon_pt[mIndex]>SingleMuonTriggerPtCut &&
	      	analysisTree.trigobject_pt[iT]>SingleMuonTriggerPtCut) { // IsoMu Leg
	    	float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    	if (dRtrig<deltaRTrigMatch){
	      	isLegMatch = true;
		}

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


	  //cout<<"made it to pass the trigger "<<endl;
if (!CutBasedTauId){
   isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
   //isoTau = analysisTree.tau_chargedIsoPtSum[tIndex];

          if (int(mIndex)!=mu_index) {
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
	  //cout<<"made it to pass the TauId "<<endl;

      if ((int)tau_index<0) continue;
      if ((int)mu_index<0) continue;

      bool isLoose = false;
      unsigned int tau_loose=-1;
      unsigned int tau_tight=-1;
      vector<int> tau; tau.clear();

 	if ((int)tau_index>-1 && analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[tau_index] > 0.5) isLoose  = true;

	tau_loose = (int)tau_index;
	
	if (!isLoose) continue;

       ta_IsoFlagVTight[0]=analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[tau_index];
       ta_IsoFlagTight[0]=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];
       ta_IsoFlagMedium[0]=analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[tau_index];
       ta_IsoFlagLoose[0]=analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[tau_index];

      double q = analysisTree.tau_charge[tau_loose] * analysisTree.muon_charge[mu_index];
      if (q > 0) continue;

	double dRmutau = deltaR(analysisTree.tau_eta[tau_loose],analysisTree.tau_phi[tau_loose],
				analysisTree.muon_eta[(int)mu_index],analysisTree.muon_phi[(int)mu_index]);
	if (dRmutau < 0.5) continue;
	  //cout<<"made it to pass the dR "<<endl;

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;


      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie];
	if (!electronMvaId && applyVetoElectronId) continue;
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

	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;

	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

	if (foundExtraMuon || foundExtraElectron || dilepton_veto)  continue;



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
      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;


    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;

	///LSF 
      if (!isData) {
	double IdIsoSF_mu = 1;
	double IdIsoSF_mu1 = 1;
	double IdIsoSF_mu2 = 1;
		
	IdIsoSF_mu=  SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);

	LSF_weight = IdIsoSF_mu;
	weight *= LSF_weight;

      }

    LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
 ////////////////////////////// now check if it is also Tight
 	if (isLoose) denominator++;

	bool isTight =false;

       if (!CutBasedTauId){
		isTight=
	  	  analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index] > 0.5;
	  	  //analysisTree.tau_chargedIsoPtSum[tau_index] < 0.8;

 
       isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
       ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tau_index];

       //isoTau = analysisTree.tau_chargedIsoPtSum[tau_index];


	 }

	if (CutBasedTauId){
		isTight=
	         analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] > 0.5;

          isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
	}

if (isTight)
	tau_tight = (int)tau_index;
      int matchedJets=0;

      float jetEta = 2.4;
      float DRmax = 0.5;
      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> bjets; bjets.clear();

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


	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if ( Drr < DRmax) cleanedJet=false;

	
	if ( Drr < 0.2) matchedJets++;

		if (Drr < 0.5 && !genTauMatched && !isData) 
		{
			
	     if (analysisTree.pfjet_flavour[jet] == 21) matchedTauToGluon = true;

	     if (abs(analysisTree.pfjet_flavour[jet]) == 1 ) matchedTauToDownQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 2 ) matchedTauToUpQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 3 ) matchedTauToStrangeQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 4 ) matchedTauToCharmQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 5 ) matchedTauToBottomQ = true;

	     if (abs(analysisTree.pfjet_flavour[jet]) == 1 || abs(analysisTree.pfjet_flavour[jet]) == 5) matchedTauToHFQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) > 1 && abs(analysisTree.pfjet_flavour[jet]) < 5) matchedTauToLFQ = true;

			}

		if (!matchedTauToGluon && !matchedTauToDownQ && !matchedTauToUpQ && !matchedTauToStrangeQ && !matchedTauToCharmQ && !matchedTauToBottomQ) matchedTauToNothing = true;

	if (!cleanedJet) continue;
	JetsV.SetPxPyPzE(0.,0.,0.,0.);
	JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
	JetsMV.push_back(JetsV);	

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

	jets.push_back((int)jet);
	counter_cleaned_jets++;
	}


      }///loop in all jets

      njets = jets.size();
      jet_count = jets.size();
      nbtag = bjets.size();

	int countjets = njets;
	int countbjets = nbtag;

      if (countjets==0) continue;
      if (countbjets>0) continue;

      LooseCFCounter[iCutL]+= weight;
      TightCFCounter[iCutT]+= weight;
      iCutL++;
      iCutT++;
      hnJets->Fill(countjets,weight);
      hnMatchedJets->Fill(matchedJets,weight);
  if( isTight)    hnMatchedJetsT->Fill(matchedJets,weight);
      hnbJets->Fill(countbjets,weight);
      hIsoMuSel->Fill(isoMuMin,weight);

      if (fabs(analysisTree.tau_charge[tau_loose]) !=1) continue;
      if (fabs(analysisTree.muon_charge[mu_index]) !=1) continue;



      if (isLoose)  hLooseIndex->Fill((int)tau_loose,weight);
	
      if (isTight)  hTightIndex->Fill((int)tau_tight,weight);
	

/////////////////// Recoil corrections
	njets = countjets;

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

	     
      if (isWJ || isDY || isWNJ || isDYNJ) {

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
 
      met_ex_recoil = pfmet_corr_x;
      met_ey_recoil = pfmet_corr_y;

      }//if isW, isDY !isData

if ((isW || isDY) && !isData) 
	{
	

	if (string::npos != Systematic.find("JetEnUp")) {
      	met_x = met_ex_JetEnUp_recoil;
      	met_y = met_ey_JetEnUp_recoil;
	}

	if (string::npos != Systematic.find("JetEnDown")) {
      	met_x = met_ex_JetEnDown_recoil;
      	met_y = met_ey_JetEnDown_recoil;
	}

	if (string::npos != Systematic.find("UnclEnUp")) {
      	met_x = met_ex_UnclusteredEnUp_recoil;
      	met_y = met_ey_UnclusteredEnUp_recoil;
	}

	if (string::npos != Systematic.find("UnclEnDown")) {
      	met_x = met_ex_UnclusteredEnDown_recoil;
      	met_y = met_ey_UnclusteredEnDown_recoil;
	}

	}

if ((!isW && !isDY) && !isData) 
	{
	

	if (string::npos != Systematic.find("JetEnUp")) {
      	met_x = met_ex_JetEnUp;
      	met_y = met_ey_JetEnUp;
	}

	if (string::npos != Systematic.find("JetEnDown")) {
      	met_x = met_ex_JetEnDown;
      	met_y = met_ey_JetEnDown;
	}

	if (string::npos != Systematic.find("UnclEnUp")) {
      	met_x = met_ex_UnclusteredEnUp;
      	met_y = met_ey_UnclusteredEnUp;
	}

	if (string::npos != Systematic.find("UnclEnDown")) {
      	met_x = met_ex_UnclusteredEnDown;
      	met_y = met_ey_UnclusteredEnDown;
	}

	}

      met_ex = met_x;
      met_ey = met_y;
    //  met_ez = analysisTree.pfmet_ez;


      met_pt = TMath::Sqrt(met_ex*met_ex + met_ey*met_ey);
      met_phi = TMath::ATan2(met_ey,met_ex);



      double dPhi=-1;double MT=-1 ; double RatioSums=-1;

      double met = sqrt ( met_ex*met_ex + met_ey*met_ey);
      

      TLorentzVector MetV; 
      MetV.SetPx(met_ex); 
      MetV.SetPy(met_ey);


      TLorentzVector muV ;  muV.SetPtEtaPhiM(analysisTree.muon_pt[mu_index], analysisTree.muon_eta[mu_index], analysisTree.muon_phi[mu_index], muonMass);
      TLorentzVector tauV;  tauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_loose], analysisTree.tau_eta[tau_loose], analysisTree.tau_phi[tau_loose], tauMass);

      TLorentzVector Wb = muV  + MetV;

      double SumPtMuMET = sqrt ( ( met_ex +analysisTree.muon_px[mu_index])*(met_ex+analysisTree.muon_px[mu_index]) + 
				 (met_ey+analysisTree.muon_py[mu_index])*(met_ey+analysisTree.muon_py[mu_index]));

      RatioSums = (Wb.Pt() - JetsMV.at(0).Pt() )/ (Wb.Pt() +JetsMV.at(0).Pt() );


      dPhi = dPhiFrom2P( muV.Px(), muV.Py(), MetV.Px(),  MetV.Py() );
      MT = TMath::Sqrt(2*muV.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));

			
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

	for (unsigned int gt = 0 ; gt < analysisTree.gentau_count; ++gt){

	  genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);

	  double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    genTauV.Eta(), genTauV.Phi());

	  if (Drr < 0.2 && analysisTree.gentau_isPrompt[gt] > 0.5  && genTauV.Pt() > 15. ) isTauMatched = true;

	}
      
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      		  if ( (abs(analysisTree.genparticles_pdgid[igen])==11 || abs(analysisTree.genparticles_pdgid[igen])==13 || abs(analysisTree.genparticles_pdgid[igen])==15 || abs(analysisTree.genparticles_pdgid[igen])<6 || abs(analysisTree.genparticles_pdgid[igen])==21)){

	  genLepV.SetXYZT(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen], analysisTree.genparticles_e[igen]);

	  double Drl=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			  genLepV.Eta(),genLepV.Phi());

      
	double DrTauLepton=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  genLepV.Eta(),genLepV.Phi());

		if (Drl < 0.2 && genLepV.Pt() > 8){

		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptTau = true;

		if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genElMatchedToTauDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genMuMatchedToTauDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genTauMatchedToTauDecay = true;
		
		if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genElMatchedHadrDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genMuMatchedHadrDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genTauMatchedHadrDecay = true;
		
		if ( (abs(analysisTree.genparticles_pdgid[igen])==1 || abs(analysisTree.genparticles_pdgid[igen])==5) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedHFQ = true;
	
		if ( (abs(analysisTree.genparticles_pdgid[igen])==2 || abs(analysisTree.genparticles_pdgid[igen])==3 || abs(analysisTree.genparticles_pdgid[igen])==4) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedLFQ = true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])==21 ) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedGluon = true;



		}

		if (DrTauLepton < 0.2 && genLepV.Pt() > 8. ) 
			
		{
			genLeptonMatched = true;

		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToElHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToMuHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToTauHadronDec= true;
				}
      			}

		}
      
      }//!isData
	
	double tauId = TauiD("mutau","MVA_Tight");
		
//	if (!isData && isTauMatched) continue;
//	if (isTauMatched && !isGenLeptonMatched) weight *= 0.9;

      //cout<<"  "<<endl;
      if (isTight) nominator++;

 //     float ptProbe = min(float(analysisTree.tau_pt[tau_loose]),float(ptBins[nPtBins]-0.1));
      float absEtaProbe = fabs(analysisTree.tau_eta[tau_loose]);
   //   int ptBin = binNumber(ptProbe,nPtBins,ptBins);
     // if (ptBin<0) continue;
     
      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
      if (etaBin<0) continue;

 //     cout<<   " found the pt "<<ptProbe<<"  "<<analysisTree.tau_pt[tau_loose]<<" "<<ptBin<<"  "<<absEtaProbe<<"  "<<etaBin<<endl;

      //cout<< "filling here  "<<analysisTree.tau_pt[tau_loose]<<"  "<<ptBin<<"  "<<etaBin<<"  "<<weight<<endl;
      //FakeRatePt[etaBin][ptBin]->Fill(analysisTree.tau_pt[tau_loose],weight);

      bool MTb = (MT>60 && MT<120);
      double ptTau1 = (double)analysisTree.tau_pt[(int)tau_loose];
      double etaTau1 = (double)analysisTree.tau_eta[(int)tau_loose];
      string WP = "MVA";
      double tfr  = TauFakeRate(ptTau1,etaTau1,Channel,WP);
      if (isTauMatched || isData) tfr = 1.;
      if (isLoose)

      {
      LooseCFCounter[iCutL]+= weight;
      iCutL++;
	      if (met>40){

if (matchedTauToDownQ) FakeRatePtIncLooseDownQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncLooseUpQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncLooseStrangeQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncLooseCharmQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncLooseBottomQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncLooseGluon[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncLooseNothing[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);


	    FakeRatePtIncLoose[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	    hRatioSum1L->Fill(RatioSums,weight);
	    hMTCut1L->Fill(MT,weight);
	    hDPhiCut1L->Fill(dPhiW, weight);
	    hMETCut1L->Fill(met, weight);
      	    hnJets1L->Fill(countjets,weight);
	    hnbJets1L->Fill(countbjets,weight);
      LooseCFCounter[iCutL]+= weight;
      iCutL++;

if (isData) { 
	FakeRatePtIncLooseDownQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseUpQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseStrangeQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseCharmQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseBottomQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseGluon[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	  if (MTb ){ 

if (matchedTauToDownQ) FakeRatePtIncLooseDownQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncLooseUpQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncLooseStrangeQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncLooseCharmQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncLooseBottomQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncLooseGluon[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncLooseNothing[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncLoose[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	    hRatioSum2L->Fill(RatioSums,weight);
	    hMTCut2L->Fill(MT,weight);
	    hDPhiCut2L->Fill(dPhiW, weight);
	    hMETCut2L->Fill(met, weight);
      	    hnJets2L->Fill(countjets,weight);
	    hnbJets2L->Fill(countbjets,weight);

      LooseCFCounter[iCutL]+= weight;
      iCutL++;
if (isData) { 
	FakeRatePtIncLooseDownQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseUpQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseStrangeQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseCharmQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseBottomQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseGluon[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	    if (dPhiW > 2.5){

if (matchedTauToDownQ) FakeRatePtIncLooseDownQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncLooseUpQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncLooseStrangeQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncLooseCharmQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncLooseBottomQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncLooseGluon[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncLooseNothing[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncLoose[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	    hRatioSum3L->Fill(RatioSums,weight);
	    hMTCut3L->Fill(MT,weight);
	    hDPhiCut3L->Fill(dPhiW, weight);
	    hMETCut3L->Fill(met, weight);
      	    hnJets3L->Fill(countjets,weight);
	    hnbJets3L->Fill(countbjets,weight);
      LooseCFCounter[iCutL]+= weight;
      iCutL++;
if (isData) { 
	FakeRatePtIncLooseDownQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseUpQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseStrangeQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseCharmQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseBottomQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseGluon[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	//if (RatioSums < 0.3 ) {
	if (countjets<3 ) {

if (matchedTauToDownQ) FakeRatePtIncLooseDownQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncLooseUpQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncLooseStrangeQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncLooseCharmQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncLooseBottomQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncLooseGluon[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncLooseNothing[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncLoose[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	    hRatioSum4L->Fill(RatioSums,weight);
	    hMTCut4L->Fill(MT,weight);
	    hDPhiCut4L->Fill(dPhiW, weight);
	    hMETCut4L->Fill(met, weight);
      	    hnJets4L->Fill(countjets,weight);
	    hnbJets4L->Fill(countbjets,weight);

      LooseCFCounter[iCutL]+= weight;
      iCutL++;
if (isData) { 
	FakeRatePtIncLooseDownQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseUpQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseStrangeQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseCharmQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseBottomQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncLooseGluon[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}
/////////// corrected for TFR
//
	    hRatioSumTFRL->Fill(RatioSums,tfr*weight);
	    hMTCutTFRL->Fill(MT,tfr*weight);
	    hDPhiCutTFRL->Fill(dPhiW, tfr*weight);
	    hMETCutTFRL->Fill(met, tfr*weight);
      	    hnJetsTFRL->Fill(countjets,weight);
	    hnbJetsTFRL->Fill(countbjets,weight);

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
if (matchedTauToDownQ) FakeRatePtIncTightDownQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncTightUpQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncTightStrangeQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncTightCharmQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncTightBottomQ[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncTightGluon[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncTightNothing[etaBin][0]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncTight[etaBin][0]->Fill(analysisTree.tau_pt[(int)tau_loose],weight);
	    hRatioSum1T->Fill(RatioSums,weight);
	    hMTCut1T->Fill(MT,weight);
	    hDPhiCut1T->Fill(dPhiW, weight);
	    hMETCut1T->Fill(met, weight);
      	    hnJets1T->Fill(countjets,weight);
	    hnbJets1T->Fill(countbjets,weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;
if (isData) { 
	FakeRatePtIncTightDownQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightUpQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightStrangeQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightCharmQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightBottomQ[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightGluon[etaBin][0] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	    if (MTb ){ 
if (matchedTauToDownQ) FakeRatePtIncTightDownQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncTightUpQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncTightStrangeQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncTightCharmQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncTightBottomQ[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncTightGluon[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncTightNothing[etaBin][1]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	    FakeRatePtIncTight[etaBin][1]->Fill(analysisTree.tau_pt[(int)tau_loose],weight);
	    hRatioSum2T->Fill(RatioSums,weight);
	    hMTCut2T->Fill(MT,weight);
	    hDPhiCut2T->Fill(dPhiW, weight);
	    hMETCut2T->Fill(met, weight);
      	    hnJets2T->Fill(countjets,weight);
	    hnbJets2T->Fill(countbjets,weight);
    
      TightCFCounter[iCutT]+= weight;
      iCutT++;
if (isData) { 
	FakeRatePtIncTightDownQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightUpQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightStrangeQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightCharmQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightBottomQ[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightGluon[etaBin][1] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	      if (dPhiW > 2.5){
if (matchedTauToDownQ) FakeRatePtIncTightDownQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncTightUpQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncTightStrangeQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncTightCharmQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncTightBottomQ[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncTightGluon[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncTightNothing[etaBin][2]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncTight[etaBin][2]->Fill(analysisTree.tau_pt[(int)tau_loose],weight);
	    hRatioSum3T->Fill(RatioSums,weight);
	    hMTCut3T->Fill(MT,weight);
	    hDPhiCut3T->Fill(dPhiW, weight);
	    hMETCut3T->Fill(met, weight);
      	    hnJets3T->Fill(countjets,weight);
	    hnbJets3T->Fill(countbjets,weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;
if (isData) { 
	FakeRatePtIncTightDownQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightUpQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightStrangeQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightCharmQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightBottomQ[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightGluon[etaBin][2] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}

	//  if (RatioSums < 0.3 ) {
	if (countjets<3 ) {
if (matchedTauToDownQ) FakeRatePtIncTightDownQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToUpQ) FakeRatePtIncTightUpQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToStrangeQ) FakeRatePtIncTightStrangeQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToCharmQ) FakeRatePtIncTightCharmQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToBottomQ) FakeRatePtIncTightBottomQ[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToGluon) FakeRatePtIncTightGluon[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
if (matchedTauToNothing) FakeRatePtIncTightNothing[etaBin][3]->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

	    FakeRatePtIncTight[etaBin][3]->Fill(analysisTree.tau_pt[(int)tau_loose],weight);
	    hRatioSum4T->Fill(RatioSums,weight);
	    hMTCut4T->Fill(MT,weight);
	    hDPhiCut4T->Fill(dPhiW, weight);
	    hMETCut4T->Fill(met, weight);
      	    hnJets4T->Fill(countjets,weight);
	    hnbJets4T->Fill(countbjets,weight);

      TightCFCounter[iCutT]+= weight;
      iCutT++;
if (isData) { 
	FakeRatePtIncTightDownQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightUpQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightStrangeQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightCharmQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightBottomQ[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);
	FakeRatePtIncTightGluon[etaBin][3] ->Fill(double(analysisTree.tau_pt[(int)tau_loose]),weight);

}
      
	    hRatioSumTFRT->Fill(RatioSums,tfr*weight);
	    hMTCutTFRT->Fill(MT,tfr*weight);
	    hDPhiCutTFRT->Fill(dPhiW, tfr*weight);
	    hMETCutTFRT->Fill(met, tfr*weight);
      	    hnJetsTFRT->Fill(countjets,weight);
	    hnbJetsTFRT->Fill(countbjets,weight);


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
  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;

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
  hnMatchedJets->Write();
  hnMatchedJetsT->Write();
  hnJets->Write();
  hnbJets->Write();
  hnJets1L->Write();
  hnbJets1L->Write();
  hnJets2L->Write();
  hnbJets2L->Write();
  hnJets3L->Write();
  hnbJets3L->Write();
  hnJets4L->Write();
  hnbJets4L->Write();
  hnJets1T->Write();
  hnbJets1T->Write();
  hnJets2T->Write();
  hnbJets2T->Write();
  hnJets3T->Write();
  hnbJets3T->Write();
  hnJets4T->Write();
  hnbJetsTFRL->Write();
  hnJetsTFRL->Write();
  hnbJetsTFRT->Write();
  hnJetsTFRT->Write();
  hnbJets4T->Write();
  CutFlowUnWTight->Write();
  CutFlowUnWLoose->Write();
  file->Write();
  file->Close();

  delete file;

}


