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
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"

float topPtWeight(float pt1,
		  float pt2) {

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter
  
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  if (pt1>400) w1 = 1;
  if (pt2>400) w2 = 1;

  return TMath::Sqrt(w1*w2);

}

float nJetsWeight(int nJets) {

  float weight = 1;
  if (nJets==0)
    weight = 1.02;
  else if (nJets==1)
    weight = 0.95;
  else 
    weight = 0.93;

  return weight;

}

void computeRecoil(float metx, float mety,
		   float unitX,float unitY,
		   float perpUnitX, float perpUnitY,
		   float dimuonPt,
		   float & recoilParal,
		   float & recoilPerp,
		   float & responseHad) {

  recoilParal = metx*unitX + mety*unitY;
  recoilPerp = metx*perpUnitX + mety*perpUnitY;
  responseHad = 1 + recoilParal/dimuonPt;
}

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

  // pile up reweighting
  const bool applyPUreweighting_vertices = cfg.get<bool>("ApplyPUreweighting_vertices");
  const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");

  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
  const bool applyRecoilCorrections =  cfg.get<bool>("ApplyRecoilCorrections");
  const bool applyRecoilOnGenerator = cfg.get<bool>("ApplyRecoilOnGenerator");
  const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
  const bool applyTopPtReweighting = cfg.get<bool>("ApplyTopPtReweighting");
  const bool applyNJetReweighting = cfg.get<bool>("ApplyNJetReweighting");
  const bool applyZMassPtReweighting = cfg.get<bool>("ApplyZMassPtReweighting");
  const bool interpolateZMassPtWeight = cfg.get<bool>("InterpolateZMassPtWeight");

  // kinematic cuts on electron
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float etaElectronHighCut = cfg.get<float>("etaElectronHighCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float dxyElectronLooseCut     = cfg.get<float>("dxyElectronLooseCut");
  const float dzElectronLooseCut      = cfg.get<float>("dzElectronLooseCut");
  const float isoElectronCut     = cfg.get<float>("isoElectronCut");
  const float isoElectronProbeCut = cfg.get<float>("isoElectronProbeCut");
  const unsigned int electronIdType  = cfg.get<unsigned int>("ElectronIdType");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool  oppositeSign   = cfg.get<bool>("OppositeSign");
  const float dielectronMassCut = cfg.get<float>("DielectronMassCut");
  const bool isoDR03         = cfg.get<bool>("IsoDR03");

  // jet related cuts
  const float jetEtaCut      = cfg.get<float>("JetEtaCut");
  const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");
  const float jetPtHighCut   = cfg.get<float>("JetPtHighCut");
  const float jetPtLowCut    = cfg.get<float>("JetPtLowCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  //trigger
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
  const string electronHLTName = cfg.get<string>("ElectronHLTName");
  const string electronHLTFilterName = cfg.get<string>("ElectronHLTFilterName");
  const string electronEle23FilterName = cfg.get<string>("ElectronEle23FilterName");
  const string electronEle17FilterName = cfg.get<string>("ElectronEle17FilterName");
  const string electronEle12FilterName = cfg.get<string>("ElectronEle12FilterName");
  const string electronSingleEleFilterName = cfg.get<string>("ElectronSingleEleFilterName");

  const float singleEleTriggerPtCut = cfg.get<float>("SingleEleTriggerPtCut");
  const float singleEleTriggerEtaCut = cfg.get<float>("SingleEleTriggerEtaCut");
  const float eleTriggerPtCut = cfg.get<float>("eleTriggerPtCut");
  const float eleTriggerEtaCut = cfg.get<float>("eleTriggerEtaCut");

  TString ElectronHLTName(electronHLTName);
  TString ElectronHLTFilterName(electronHLTFilterName);
  TString ElectronEle23FilterName(electronEle23FilterName);
  TString ElectronEle17FilterName(electronEle17FilterName);
  TString ElectronEle12FilterName(electronEle12FilterName);
  TString ElectronSingleEleFilterName(electronSingleEleFilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");
  const string ElectronTrigFile  = cfg.get<string>("ElectronTrigEff"); 
  const string trackingSFFile = cfg.get<string>("TrackingSFFile");
  

  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  const string recoilPuppiFileName   = cfg.get<string>("RecoilPuppiFileName");
  TString RecoilPuppiFileName(recoilPuppiFileName);

  const string recoilMvaFileName   = cfg.get<string>("RecoilMvaFileName");
  TString RecoilMvaFileName(recoilMvaFileName);

  // systematics
  const string metsysFileName   = cfg.get<string>("MetSysFileName");
  TString MetSysFileName(metsysFileName);

  const string metsysPuppiFileName   = cfg.get<string>("MetSysPuppiFileName");
  TString MetSysPuppiFileName(metsysPuppiFileName);

  const string metsysMvaFileName   = cfg.get<string>("MetSysMvaFileName");
  TString MetSysMvaFileName(metsysMvaFileName);

  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);

  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

  const int bkgdType =  cfg.get<float>("BkgdType");

  const string jsonFile = cfg.get<string>("jsonFile");

  const int applyJES = cfg.get<int>("applyJES");

  const float eleMomScaleBarrel = cfg.get<float>("EleMomScaleBarrel");
  const float eleMomScaleEndcap = cfg.get<float>("EleMomScaleEndcap");

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  // **** end of configuration
  

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
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl; 

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsSkimmedH = new TH1D("histWeightsSkimmedH","",1,-0.5,0.5);

  TH1D * massZH = new TH1D("massZH","",1000,0,1000);
  TH1D * ptZH = new TH1D("ptZH","",1000,0,1000);
  TH2D * massPtZH = new TH2D("massPtZH","",100,0,1000,100,0,1000);


  TH1D * massSelH = new TH1D("massSelH","",200,0,200);
  TH1D * massExtendedSelH = new TH1D("massExtendedSelH","",500,0,5000);
  // Histograms after final selection
  TH1D * ptLeadingEleSelH = new TH1D("ptLeadingEleSelH","",100,0,200);
  TH1D * ptTrailingEleSelH = new TH1D("ptTrailingEleSelH","",100,0,200);
  TH1D * etaLeadingEleSelH = new TH1D("etaLeadingEleSelH","",50,-2.5,2.5);
  TH1D * etaTrailingEleSelH = new TH1D("etaTrailingEleSelH","",50,-2.5,2.5);
  TH1D * dielectronPtSelH = new TH1D("dielectronPtSelH","",100,0,1000);
  TH1D * dielectronEtaSelH = new TH1D("dielectronEtaSelH","",120,-6,6);

  TH1D * metSelH = new TH1D("metSelH","",200,0.,400.);
  TH1D * puppimetSelH = new TH1D("puppimetSelH","",200,0.,400.);
  TH1D * mvametSelH = new TH1D("mvametSelH","",200,0.,400.);

  TH1D * metZSelH  = new TH1D("metZSelH","",200,0,400);
  TH1D * puppimetZSelH = new TH1D("puppimetZSelH","",200,0,400);
  TH1D * mvametZSelH = new TH1D("mvametZSelH","",200,0,400);

  TString scales[21] = {"M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","0",
			"P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};

  TH1D * massSelScaleH[21];
  TH1D * metSelScaleH[21];
  TH1D * puppimetSelScaleH[21];
  TH1D * mvametSelScaleH[21];
  for (int iScale=0; iScale<21; ++iScale) {    
    massSelScaleH[iScale] = new TH1D("massSel"+scales[iScale]+"H","",200,0,200);
    metSelScaleH[iScale] = new TH1D("metSel"+scales[iScale]+"H","",200,0,400);
    puppimetSelScaleH[iScale] = new TH1D("puppimetSel"+scales[iScale]+"H","",200,0,400);
    mvametSelScaleH[iScale] = new TH1D("mvametSel"+scales[iScale]+"H","",200,0,400);
  }

  int nJetBins = 3;
  int nZPtBins = 5;
  float zPtBins[6] = {0,10,20,30,50,1000};
  float jetBins[4] = {-0.5,0.5,1.5,2.5};

  TString NJetBins[3] = {"NJet0","NJet1","NJetGe2"};
  TString ZPtBins[5] = {"Pt0to10",
			"Pt10to20",
			"Pt20to30",
			"Pt30to50",
			"PtGt50"};

  TString RecoilZParal("recoilZParal_");
  TString RecoilZPerp("recoilZPerp_");
  TString RecoilPuppiZParal("recoilPuppiZParal_");
  TString RecoilPuppiZPerp("recoilPuppiZPerp_");
  TString RecoilMvaZParal("recoilMvaZParal_");
  TString RecoilMvaZPerp("recoilMvaZPerp_");

  TString RecoilTopParal("recoilTopParal_");
  TString RecoilTopPerp("recoilTopPerp_");
  TString RecoilPuppiTopParal("recoilPuppiTopParal_");
  TString RecoilPuppiTopPerp("recoilPuppiTopPerp_");
  TString RecoilMvaTopParal("recoilMvaTopParal_");
  TString RecoilMvaTopPerp("recoilMvaTopPerp_");

  // Saving Z pt bins
  TH1D * ZPtBinsH = new TH1D("ZPtBinsH","ZPtBinsH",nZPtBins,zPtBins);
  for (int iB=0; iB<nZPtBins; ++iB) 
    ZPtBinsH->GetXaxis()->SetBinLabel(iB+1,ZPtBins[iB]);
  
  // Saving jet bins
  TH1D * JetBinsH = new TH1D("JetBinsH","JetBinsH",nJetBins,jetBins);
  for (int iB=0; iB<nJetBins; ++iB)
    JetBinsH->GetXaxis()->SetBinLabel(iB+1,NJetBins[iB]);

  TH1D * metSelNJets[3];
  TH1D * puppimetSelNJets[3];
  TH1D * mvametSelNJets[3];
  
  TH1D * metZSelNJets[3];
  TH1D * puppimetZSelNJets[3];
  TH1D * mvametZSelNJets[3];

  TH1D * metTopSelNJets[3];
  TH1D * puppimetTopSelNJets[3];
  TH1D * mvametTopSelNJets[3];

  TH1D * recoilZParalH[3];
  TH1D * recoilZPerpH[3];
  TH1D * recoilPuppiZParalH[3];
  TH1D * recoilPuppiZPerpH[3];
  TH1D * recoilMvaZParalH[3];
  TH1D * recoilMvaZPerpH[3];

  TH1D * recoilTopParalH[3];
  TH1D * recoilTopPerpH[3];
  TH1D * recoilPuppiTopParalH[3];
  TH1D * recoilPuppiTopPerpH[3];
  TH1D * recoilMvaTopParalH[3];
  TH1D * recoilMvaTopPerpH[3];

  TH1D * recoilZParal_Ptbins_nJetsH[3][5];
  TH1D * recoilZPerp_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiZParal_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiZPerp_Ptbins_nJetsH[3][5];
  TH1D * recoilMvaZParal_Ptbins_nJetsH[3][5];
  TH1D * recoilMvaZPerp_Ptbins_nJetsH[3][5];

  TH1D * recoilResponse[3];
  TH1D * recoilResponse_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiResponse[3];
  TH1D * recoilPuppiResponse_Ptbins_nJetsH[3][5];
  TH1D * recoilMvaResponse[3];
  TH1D * recoilMvaResponse_Ptbins_nJetsH[3][5];
  
    // systematics ===>

  TH1D * recoilMvaZParalSysH[3][4][2];
  TH1D * recoilMvaZPerpSysH[3][4][2];

  TH1D * recoilMvaTopParalSysH[3][4][2];
  TH1D * recoilMvaTopPerpSysH[3][4][2];

  TH1D * mvametZSelNJetsSysH[3][4][2];
  TH1D * mvametTopSelNJetsSysH[3][4][2];

  TString sysNames[4] = {"JES","MEtScale","MEtReso","MuonMom"};
  TString UpDown[2] = {"Up","Down"};


  for (int iBin=0; iBin<nJetBins; ++iBin) {
    recoilZParalH[iBin] = new TH1D(RecoilZParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilZPerpH[iBin] = new TH1D(RecoilZPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilPuppiZParalH[iBin] = new TH1D(RecoilPuppiZParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilPuppiZPerpH[iBin] = new TH1D(RecoilPuppiZPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilMvaZParalH[iBin] = new TH1D(RecoilMvaZParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilMvaZPerpH[iBin] = new TH1D(RecoilMvaZPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilTopParalH[iBin] = new TH1D(RecoilTopParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilTopPerpH[iBin] = new TH1D(RecoilTopPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilPuppiTopParalH[iBin] = new TH1D(RecoilPuppiTopParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilPuppiTopPerpH[iBin] = new TH1D(RecoilPuppiTopPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilMvaTopParalH[iBin] = new TH1D(RecoilMvaTopParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilMvaTopPerpH[iBin] = new TH1D(RecoilMvaTopPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilResponse[iBin] = new TH1D("recoilResponse"+NJetBins[iBin]+"H","",200,-10,10);
    recoilPuppiResponse[iBin] = new TH1D("recoilPuppiResponse"+NJetBins[iBin]+"H","",200,-10,10);
    recoilMvaResponse[iBin] = new TH1D("recoilMvaResponse"+NJetBins[iBin]+"H","",200,-10,10);

    metSelNJets[iBin] = new TH1D("metSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetSelNJets[iBin] = new TH1D("puppimetSel"+NJetBins[iBin]+"H","",200,0.,400.);
    mvametSelNJets[iBin] = new TH1D("mvametSel"+NJetBins[iBin]+"H","",200,0.,400.);

    metZSelNJets[iBin] = new TH1D("metZSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetZSelNJets[iBin] = new TH1D("puppimetZSel"+NJetBins[iBin]+"H","",200,0.,400.);
    mvametZSelNJets[iBin] = new TH1D("mvametZSel"+NJetBins[iBin]+"H","",200,0.,400.);

    metTopSelNJets[iBin] = new TH1D("metTopSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetTopSelNJets[iBin] = new TH1D("puppimetTopSel"+NJetBins[iBin]+"H","",200,0.,400.);
    mvametTopSelNJets[iBin] = new TH1D("mvametTopSel"+NJetBins[iBin]+"H","",200,0.,400.);

    for (int iSys=0; iSys<4; ++iSys) {
      for (int iUD=0; iUD<2; ++iUD) {

	mvametZSelNJetsSysH[iBin][iSys][iUD] = new TH1D("mvametZSel"+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,0.,400.);
	mvametTopSelNJetsSysH[iBin][iSys][iUD] = new TH1D("mvametTopSel"+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,0.,400.);
	
	recoilMvaZParalSysH[iBin][iSys][iUD] = new TH1D(RecoilMvaZParal+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,-400,400);
	recoilMvaZPerpSysH[iBin][iSys][iUD] = new TH1D(RecoilMvaZPerp+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,-400,400);

	recoilMvaTopParalSysH[iBin][iSys][iUD] = new TH1D(RecoilMvaTopParal+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,-400,400);
	recoilMvaTopPerpSysH[iBin][iSys][iUD] = new TH1D(RecoilMvaTopPerp+NJetBins[iBin]+sysNames[iSys]+UpDown[iUD]+"H","",200,-400,400);

      }
    }

  }
  
  for (int iJets=0; iJets<nJetBins; ++iJets) {
    for (int iPtBins=0; iPtBins<nZPtBins; ++iPtBins){

      recoilZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);
      recoilZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);

      recoilPuppiZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);
      recoilPuppiZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);

      recoilMvaZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilMvaZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);
      recoilMvaZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilMvaZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);

      recoilResponse_Ptbins_nJetsH[iJets][iPtBins] = new TH1D("recoilResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-10,10);
      recoilPuppiResponse_Ptbins_nJetsH[iJets][iPtBins] = new TH1D("recoilPuppiResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-10,10);
      recoilMvaResponse_Ptbins_nJetsH[iJets][iPtBins] = new TH1D("recoilMvaResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-10,10);

    }
  }

  
  TH1D * nJets30SelH    = new TH1D("nJets30SelH","",11,-0.5,10.5);
  TH1D * nJets30etaCutSelH = new TH1D("nJets30etaCutSelH","",11,-0.5,10.5);
  TH1D * nJets20SelH    = new TH1D("nJets20SelH","",11,-0.5,10.5);
  TH1D * nJets20etaCutSelH = new TH1D("nJets20etaCutSelH","",11,-0.5,10.5);

  TH1D * HT30SelH       = new TH1D("HT30SelH","",50,0,500);
  TH1D * HT30etaCutSelH = new TH1D("HT30etaCutSelH","",50,0,500);
  TH1D * HT20SelH       = new TH1D("HT20SelH","",50,0,500);
  TH1D * HT20etaCutSelH = new TH1D("HT20etaCutSelH","",50,0,500);
  
  TH1F * NumberOfVerticesH = new TH1F("NumberOfVerticesH","",51,-0.5,50.5);
  
  TH1D * EleSF_IdIso_Ele1H = new TH1D("EleIdIsoSF_Ele1H", "EleIdIsoSF_Ele1", 100, 0.5,1.5);
  TH1D * EleSF_IdIso_Ele2H = new TH1D("EleIdIsoSF_Ele2H", "EleIdIsoSF_Ele2", 100, 0.5,1.5);

  int nPtBins = 6;
  float ptBins[7] = {13,20,25,30,40,60,1000};

  int nPtBinsTrig = 16;
  float ptBinsTrig[17] = {10,
			  13,
			  16,
			  19,
			  22,
			  25,
			  28,
			  31,
			  34,
			  37,
			  40,
			  45,
			  50,
			  60,
			  70,
			  100,
			  1000};  
  
  int nEtaBins = 3;
  float etaBins[4] = {0,1.48,2.1,2.5}; 
  
  TString PtBins[6] = {"Pt13to20",
		       "Pt20to25",
		       "Pt25to30",
		       "Pt30to40",
		       "Pt40to60",
		       "PtGt60"};
  

  TString EtaBins[3] = {"EtaLt1p48",
			"Eta1p48to2p1",
			"EtaGt2p1"};

  TString PtBinsTrig[16] = {"Pt10to13",
			    "Pt13to16",
			    "Pt16to19",
			    "Pt19to22",
			    "Pt22to25",
			    "Pt25to28",
			    "Pt28to31",
			    "Pt31to34",
			    "Pt34to37",
			    "Pt37to40",
			    "Pt40to45",
			    "Pt45to50",
			    "Pt50to60",
			    "Pt60to70",
			    "Pt70to100",
			    "PtGt100"};

  TString JetBins[3] = {"Jet0","Jet1","JetGe2"};

  //*****  create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //



  TH1F * etaBinsH = new TH1F("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  etaBinsH->Draw();
  file->cd();
  etaBinsH->Write("etaBinsH");

  //*****  create pt histogram_s with pt ranges associated to their names (eg. Pt10to13, ..)   ***** //
  //*****  two different pT binning, one for IdIso and one for trigger   ***** //

  TH1D * ptBinsH =  new TH1D("ptBinsH", "ptBinsH", nPtBins, ptBins);
  ptBinsH->Draw();
  ptBinsH->GetXaxis()->Set(nPtBins, ptBins);
  for (int i=0; i<nPtBins; i++){ ptBinsH->GetXaxis()->SetBinLabel(i+1, PtBins[i]);}
  ptBinsH->Draw();
  file->cd();
  ptBinsH->Write("ptBinsH");

  TH1D * ptBinsTrigH =  new TH1D("ptBinsTrigH", "ptBinsTrigH", nPtBinsTrig, ptBinsTrig);
  ptBinsTrigH->Draw();
  ptBinsTrigH->GetXaxis()->Set(nPtBinsTrig, ptBinsTrig);
  for (int i=0; i<nPtBinsTrig; i++){ ptBinsTrigH->GetXaxis()->SetBinLabel(i+1, PtBinsTrig[i]);}
  ptBinsTrigH->Draw();
  file->cd();
  ptBinsTrigH->Write("ptBinsTrigH");

  TH1D * ZMassPass = new TH1D("ZMassPass","",80,50,130);
  TH1D * ZMassFail = new TH1D("ZMassFail","",80,50,130);

  TH1F * ZMassJetEtaPtPass[3][3][6];
  TH1F * ZMassJetEtaPtFail[3][3][6];

  TH1F * ZMassEtaPtPass[3][6];
  TH1F * ZMassEtaPtFail[3][6];

  TH1F * PromptPtPass[3];
  TH1F * PromptPtFail[3];

  TH1F * NonPromptPtPass[3];
  TH1F * NonPromptPtFail[3];

  TH1F * PromptSelPtPass[3];
  TH1F * PromptSelPtFail[3];

  TH1F * NonPromptSelPtPass[3];
  TH1F * NonPromptSelPtFail[3];

  TH1F * PromptSelJetPtPass[3][3];
  TH1F * PromptSelJetPtFail[3][3];

  TH1F * NonPromptSelJetPtPass[3][3];
  TH1F * NonPromptSelJetPtFail[3][3];

  TH1F * ZMassEle23EtaPtPass[3][16];
  TH1F * ZMassEle23EtaPtFail[3][16];

  TH1F * ZMassEle17EtaPtPass[3][16];
  TH1F * ZMassEle17EtaPtFail[3][16];

  TH1F * ZMassEle12EtaPtPass[3][16];
  TH1F * ZMassEle12EtaPtFail[3][16];

  TH1F * ZMassIsoEleEtaPtPass[3][16];
  TH1F * ZMassIsoEleEtaPtFail[3][16];

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    PromptPtPass[iEta] = new TH1F("PromptPt"+EtaBins[iEta]+"Pass","",nPtBins,ptBins);
    PromptPtFail[iEta] = new TH1F("PromptPt"+EtaBins[iEta]+"Fail","",nPtBins,ptBins);
    NonPromptPtPass[iEta] = new TH1F("NonPromptPt"+EtaBins[iEta]+"Pass","",nPtBins,ptBins);
    NonPromptPtFail[iEta] = new TH1F("NonPromptPt"+EtaBins[iEta]+"Fail","",nPtBins,ptBins);
    PromptSelPtPass[iEta] = new TH1F("PromptSelPt"+EtaBins[iEta]+"Pass","",nPtBins,ptBins);
    PromptSelPtFail[iEta] = new TH1F("PromptSelPt"+EtaBins[iEta]+"Fail","",nPtBins,ptBins);
    NonPromptSelPtPass[iEta] = new TH1F("NonPromptSelPt"+EtaBins[iEta]+"Pass","",nPtBins,ptBins);
    NonPromptSelPtFail[iEta] = new TH1F("NonPromptSelPt"+EtaBins[iEta]+"Fail","",nPtBins,ptBins);
    for (int iJet=0; iJet<3; ++iJet) {
      PromptSelJetPtPass[iEta][iJet] = new TH1F("PromptSelPt"+EtaBins[iEta]+JetBins[iJet]+"Pass","",nPtBins,ptBins);
      PromptSelJetPtFail[iEta][iJet] = new TH1F("PromptSelPt"+EtaBins[iEta]+JetBins[iJet]+"Fail","",nPtBins,ptBins);
      NonPromptSelJetPtPass[iEta][iJet] = new TH1F("NonPromptSelPt"+EtaBins[iEta]+JetBins[iJet]+"Pass","",nPtBins,ptBins);
      NonPromptSelJetPtFail[iEta][iJet] = new TH1F("NonPromptSelPt"+EtaBins[iEta]+JetBins[iJet]+"Fail","",nPtBins,ptBins);
    }
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassEtaPtPass[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
      ZMassEtaPtFail[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      for (int iJet=0; iJet<3; ++iJet) {
	ZMassJetEtaPtPass[iEta][iJet][iPt] = new TH1F("ZMass"+EtaBins[iEta]+JetBins[iJet]+PtBins[iPt]+"Pass","",80,50,130);
	ZMassJetEtaPtFail[iEta][iJet][iPt] = new TH1F("ZMass"+EtaBins[iEta]+JetBins[iJet]+PtBins[iPt]+"Fail","",80,50,130);
      }
    }
    for (int iPt=0; iPt<nPtBinsTrig; ++iPt) {
      ZMassEle23EtaPtPass[iEta][iPt] = new TH1F("ZMassEle23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",80,50,130);
      ZMassEle23EtaPtFail[iEta][iPt] = new TH1F("ZMassEle23"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",80,50,130);
      ZMassEle17EtaPtPass[iEta][iPt] = new TH1F("ZMassEle17"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",80,50,130);
      ZMassEle17EtaPtFail[iEta][iPt] = new TH1F("ZMassEle17"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",80,50,130);
      ZMassEle12EtaPtPass[iEta][iPt] = new TH1F("ZMassEle12"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",80,50,130);
      ZMassEle12EtaPtFail[iEta][iPt] = new TH1F("ZMassEle12"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",80,50,130);
      ZMassIsoEleEtaPtPass[iEta][iPt] = new TH1F("ZMassIsoEle"+EtaBins[iEta]+PtBinsTrig[iPt]+"Pass","",80,50,130);
      ZMassIsoEleEtaPtFail[iEta][iPt] = new TH1F("ZMassIsoEle"+EtaBins[iEta]+PtBinsTrig[iPt]+"Fail","",80,50,130);
    }
  }



  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);

  ULong64_t evt;
  float pt_1;
  float pt_2;
  float eta_1;
  float eta_2;
  float phi_1;
  float phi_2;

  float dilepton_pt;
  float dilepton_eta;
  float dilepton_phi;

  float mvametX;
  float metX;
  float metuncorrX;
  float puppimetX;

  float u1_mvametX;
  float u1_metX;
  float u1_metuncorrX;

  float u2_mvametX;
  float u2_metX;
  float u2_metuncorrX;

  int njets;
  bool match_1;
  bool match_2;

  TTree * synch_tree = new TTree("TauCheck","TauChek");

  synch_tree->Branch("evt", &evt, "evt/I");

  synch_tree->Branch("dilepton_pt",&dilepton_pt,"dilepton_pt/F");
  synch_tree->Branch("dilepton_eta",&dilepton_eta,"dilepton_eta/F");
  synch_tree->Branch("dilepton_phi",&dilepton_phi,"dilepton_phi/F");

  synch_tree->Branch("pt_1",&pt_1,"pt_1/F");
  synch_tree->Branch("pt_2",&pt_2,"pt_2/F");

  synch_tree->Branch("eta_1",&eta_1,"eta_1/F");
  synch_tree->Branch("eta_2",&eta_2,"eta_2/F");

  synch_tree->Branch("phi_1",&phi_1,"phi_1/F");
  synch_tree->Branch("phi_2",&phi_2,"phi_2/F");
  
  synch_tree->Branch("match_1",&match_1,"match_1/O");
  synch_tree->Branch("match_2",&match_2,"match_2/O");

  synch_tree->Branch("mvamet",&mvametX,"mvamet/F");
  synch_tree->Branch("met",&metX,"met/F");
  synch_tree->Branch("metuncorr",&metuncorrX,"metuncorr/F");
  synch_tree->Branch("puppimet",&puppimetX,"puppimet/F");

  synch_tree->Branch("u1_mvamet",&u1_mvametX,"u1_mvamet/F");
  synch_tree->Branch("u1_met",&u1_metX,"u1_met/F");
  synch_tree->Branch("u1_metuncorr",&u1_metuncorrX,"u1_metuncorr/F");

  synch_tree->Branch("u2_mvamet",&u2_mvametX,"u2_mvamet/F");
  synch_tree->Branch("u2_met",&u2_metX,"u2_met/F");
  synch_tree->Branch("u2_metuncorr",&u2_metuncorrX,"u2_metuncorr/F");

  synch_tree->Branch("njets", &njets, "njets/I");

  // PILE UP REWEIGHTING - OPTIONS

  if (applyPUreweighting_vertices and applyPUreweighting_official) 
	{std::cout<<"ERROR: Choose only ONE PU reweighting method (vertices or official, not both!) " <<std::endl; exit(-1);}

  // reweighting with vertices

  // reading vertex weights
  TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/"+vertDataFileName);
  TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/"+vertMcFileName);

  TH1F * hvertWeight = new TH1F("hvertWeight","",40,0,1);
  TH1F * hNvert = new TH1F("hNvert","",51,-0.5,50.5);
  
  TH1D *hNvertData = (TH1D*)fileDataNVert->Get(TString(vertHistName));
  TH1D *hNvertMC =(TH1D*)fileMcNVert->Get(TString(vertHistName));
  Float_t normData = hNvertData->GetSumOfWeights();   
  Float_t normMC =  hNvertMC->GetSumOfWeights();
  hNvertData->Scale(1/normData);
  hNvertMC->Scale(1/normMC);


  // reweighting official recipe 
  // initialize pile up object
  PileUp * PUofficial = new PileUp();
  
  if (applyPUreweighting_official) {
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read"); 
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read"); 
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUofficial->set_h_data(PU_data); 
    PUofficial->set_h_MC(PU_mc);
  }

  // Met recoil corrections
  RecoilCorrector recoilPFMetCorrector(RecoilFileName);
  RecoilCorrector recoilMvaMetCorrector(RecoilMvaFileName);
  RecoilCorrector recoilPuppiMetCorrector(RecoilPuppiFileName);

  // MetResponse
  MEtSys metSys(MetSysFileName);
  MEtSys metSysPuppi(MetSysPuppiFileName);
  MEtSys metSysMva(MetSysMvaFileName);


  // Lepton Scale Factors 
  ScaleFactor * SF_electronIdIso;
  ScaleFactor * SF_electronTrig;
  if (applyLeptonSF) {
    SF_electronIdIso = new ScaleFactor();
    SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
    SF_electronTrig = new ScaleFactor();
    SF_electronTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronTrigFile));
    
  }
  // tracking efficiency SF                                                                                                                                                           
  TFile * fileTrackingSF = new TFile(TString(cmsswBase)+"/src/"+TString(trackingSFFile));
  TH1D * trackEffEleH = (TH1D*)fileTrackingSF->Get("effTrackingE");



  // Z mass pt weights
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName); 
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "/src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName); 
  if (histZMassPtWeights==NULL) {
    std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName << std::endl;
    exit(-1);
  }
  
  int nFiles = 0;
  int nEvents = 0;
  int selEventsDielectrons = 0;
  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 99999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
		
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;
    
    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    if (file_->IsZombie()) {std::cout << "Failed to open file "<< filen << std::endl; exit(-1);} 

    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) {
      cout << "No histogram makeroottree/nEvents is found, skipping file" << endl;
      continue;
    }
    int NE = int(histoInputEvents->GetEntries());
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
    std::cout << "      number of input events         = " << NE << std::endl;

    TTree * _inittree = NULL;
    _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree==NULL) { 
      cout << "No " << initNtupleName << " is found, skipping file" << endl;
      continue;
    }
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
    if (_tree==NULL) { 
      cout << "No " << ntupleName << " is found, skipping file" << endl;
      continue;
    }
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);
  
    // EVENT LOOP //
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 

      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      
      float weight = 1;

      if (!isData) {
	weight *= analysisTree.genweight;
	hNvert->Fill(analysisTree.primvertex_count,weight);
      }
      histWeightsSkimmedH->Fill(float(0),weight);

      TLorentzVector genZ; genZ.SetXYZM(0,0,0,91.2); 
      TLorentzVector genV; genV.SetXYZM(0,0,0,0);
      TLorentzVector genL; genL.SetXYZM(0,0,0,0);
      TLorentzVector genMuonsLV; genMuonsLV.SetXYZT(0,0,0,0);
      TLorentzVector genElectronsLV; genElectronsLV.SetXYZT(0,0,0,0);
      TLorentzVector genTausLV; genTausLV.SetXYZT(0,0,0,0);
      TLorentzVector genVisTausLV; genVisTausLV.SetXYZT(0,0,0,0);
      std::vector<unsigned int> genMuons; genMuons.clear();
      std::vector<unsigned int> genElectrons; genElectrons.clear();
      std::vector<unsigned int> genTaus; genTaus.clear();
      std::vector<unsigned int> genVisTaus; genVisTaus.clear();
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
	    genTaus.push_back(igentau);
	    genTausLV += tauLV;
	  }
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {	
	    genVisTaus.push_back(igentau);
	    genVisTausLV += tauVisLV;
	  }
	  
	}

	for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	  //	  cout << igen << "   pdgId = " << analysisTree.genparticles_pdgid[igen] << endl;
	  TLorentzVector genPart; genPart.SetXYZT(analysisTree.genparticles_px[igen],
						  analysisTree.genparticles_py[igen],
						  analysisTree.genparticles_pz[igen],
						  analysisTree.genparticles_e[igen]);
	  if (analysisTree.genparticles_pdgid[igen]==23) {
	    if (analysisTree.genparticles_fromHardProcess[igen])
	      genZ.SetXYZT(analysisTree.genparticles_px[igen],
			   analysisTree.genparticles_py[igen],
			   analysisTree.genparticles_pz[igen],
			   analysisTree.genparticles_e[igen]);
	  }
	  bool isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
	  bool isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
	  bool isLepton = isMuon || isElectron;
	  bool isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	    fabs(analysisTree.genparticles_pdgid[igen])==14||
	    fabs(analysisTree.genparticles_pdgid[igen])==16;
	  bool isPrompt = analysisTree.genparticles_isPrompt[igen]||
	    analysisTree.genparticles_isPromptTauDecayProduct[igen];

	  if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      genElectrons.push_back(igen);
	      genElectronsLV += genPart;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      genMuons.push_back(igen);
	      genMuonsLV += genPart;
	    }
	  }
	  
	  if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
	    if (isLepton&&
		fabs(genPart.Eta())<2.4&&
		genPart.Pt()>10) {
	      genV += genPart;
	      genL += genPart;
	    }
	    if (isNeutrino) 
	      genV += genPart;
	  }
	}


	if (genV.Pt()<0.1) genV.SetXYZM(0.1,0.1,0.,0.);

	/*
	std::cout << "Prompt electrons = " << genElectrons.size() << std::endl;
	std::cout << "Prompt muons     = " << genMuons.size() << std::endl;
	std::cout << "Prompt taus      = " << genTaus.size() << std::endl;
	std::cout << "Prompt vis taus  = " << genVisTaus.size() << std::endl;
	printf("gen Z Boson   -> px = %7.1f  py = %7.1f  pz = %7.1f\n",genZ.Px(),genZ.Py(),genZ.Pz());
	printf("gen muons     -> px = %7.1f  py = %7.1f  pz = %7.1f\n",genMuonsLV.Px(),genMuonsLV.Py(),genMuonsLV.Pz());
	printf("gen electrons -> px = %7.1f  py = %7.1f  pz = %7.1f\n",genElectronsLV.Px(),genElectronsLV.Py(),genElectronsLV.Pz());
	printf("gen taus      -> px = %7.1f  py = %7.1f  pz = %7.1f\n",genTausLV.Px(),genTausLV.Py(),genTausLV.Pz());
	printf("gen vis taus  -> px = %7.1f  py = %7.1f  pz = %7.1f\n",genVisTausLV.Px(),genVisTausLV.Py(),genVisTausLV.Pz());
	std::cout << std::endl;
	*/	

	float genZPt = -1.0;
	float genZMass = -1.0;

	if (genElectrons.size()==2) {
	  genZPt   = genElectronsLV.Pt();
	  genZMass = genElectronsLV.M();
	}
	else if (genMuons.size()==2) {
	  genZPt   = genMuonsLV.Pt();
          genZMass = genMuonsLV.M();
	}
	else if (genTaus.size()==2) {
          genZPt   = genTausLV.Pt();
          genZMass = genTausLV.M();
        }

	massZH->Fill(genZMass,weight);
	ptZH->Fill(genZPt,weight);
	massPtZH->Fill(genZMass,genZPt,weight);

	if (applyZMassPtReweighting) {
	  if (genZMass>1000) genZMass=999;
	  if (genZPt>1000) genZPt=999;
	  if (genZMass>50.0&&genZPt>0.0) {
	    float dyWeight = 1;
	    if (interpolateZMassPtWeight) 
	      dyWeight = histZMassPtWeights->Interpolate(genZMass,genZPt);
	    else 
	      dyWeight *= histZMassPtWeights->GetBinContent(histZMassPtWeights->FindBin(genZMass,genZPt));
	    
	    //	    std::cout << "Z mass = " << genZMass << "   Z Pt = " << genZPt << "   weight = " << dyWeight << std::endl;
	    weight *= dyWeight;
	  }
	  
	}

	if (applyPUreweighting_vertices) {
	
	  int binNvert = hNvert->FindBin(analysisTree.primvertex_count);
	  float_t dataNvert = hNvertData->GetBinContent(binNvert);
	  float_t mcNvert = hNvertMC->GetBinContent(binNvert);
	  if (mcNvert < 1e-10){mcNvert=1e-10;}
	  float_t vertWeight = dataNvert/mcNvert;
	  hvertWeight->Fill(vertWeight);
	  weight *= vertWeight;
	  //	cout << "NVert = " << analysisTree.primvertex_count << "  : " << vertWeight << endl;
	}

	// PU reweighting with Ninteractions (official recipe) 
	if (applyPUreweighting_official) {
	  double Ninteractions = analysisTree.numtruepileupinteractions;
	  double PUweight = PUofficial->get_PUweight(Ninteractions);
	  weight *= PUweight;
	  PUweightsOfficialH->Fill(PUweight);
	  
	}
	if (applyTopPtReweighting) {
	  float topPt = -1;
	  float antitopPt = -1;
	  for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

	    if (analysisTree.genparticles_pdgid[igen]==6)
	      topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    
	    if (analysisTree.genparticles_pdgid[igen]==-6)
	      antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    
	    
	  }
	  if (topPt>0&&antitopPt>0) {
	    float topptweight = topPtWeight(topPt,antitopPt);
	    weight *= topptweight;
	  }
	}
      }

      if (isData && applyGoodRunSelection){

     	bool lumi = false;
        int n=analysisTree.event_run;
        int lum = analysisTree.event_luminosityblock;
	
	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {
	    
	    if ( num.c_str() ==  a.name ) {
	      //	      std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      
	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		
		//   cout<<b->lower<<"  "<<b->bigger<<endl;
                if (lum  >= b->lower && lum <= b->bigger ) lumi = true;

	      }
	      auto last = std::prev(a.ranges.end());
	      // std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	      

	    }
	    
	  }
	if (!lumi) continue;
	
      }

      
      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;
      
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

      bool isTriggerElectron = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(ElectronHLTName)) {
	  //	  std::cout << it->first << " : " << it->second << std::endl;
	  if (it->second==1)
            isTriggerElectron = true;
	}
      }

      if (applyTrigger && !isTriggerElectron) continue;

      unsigned int nElectronFilter = 0;
      bool isElectronFilter = false;
      
      unsigned int nEle23Filter = 0;
      bool isEle23Filter = false;
      
      unsigned int nEle17Filter = 0;
      bool isEle17Filter = false;
      
      unsigned int nEle12Filter = 0;
      bool isEle12Filter = false;
      
      unsigned int nSingleEleFilter = 0;
      bool isSingleEleFilter = false;
      
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==ElectronHLTFilterName) {
	  nElectronFilter = i;
	  isElectronFilter = true;
	}
	if (HLTFilter==ElectronEle23FilterName) {
	  nEle23Filter = i;
	  isEle23Filter = true;
	}
	if (HLTFilter==ElectronEle17FilterName) {
	  nEle17Filter = i;
	  isEle17Filter = true;
	}
	if (HLTFilter==ElectronEle12FilterName) {
	  nEle12Filter = i;
	  isEle12Filter = true;
	}
	if (HLTFilter==ElectronSingleEleFilterName) {
	  nSingleEleFilter = i;
	  isSingleEleFilter = true;
	}
      }

      if (applyTrigger) {
	if (!isElectronFilter) {
	  cout << "Filter " << ElectronHLTFilterName << " not found" << endl;
	  exit(-1);
	}
	if (!isEle23Filter) {
	  cout << "Filter " << ElectronEle23FilterName << " not found" << endl;
	  exit(-1);
	}
	if (!isEle17Filter) {
	  cout << "Filter " << ElectronEle17FilterName << " not found" << endl;
	  exit(-1);
	}
	if (!isEle12Filter) {
	  cout << "Filter " << ElectronEle12FilterName << " not found" << endl;
	  exit(-1);
	}
	if (!isSingleEleFilter) {
	  cout << "Filter " << ElectronSingleEleFilterName << " not found" << endl;
	  exit(-1);
	}
      }
      //      std::cout << "Filters detected..." << std::endl;

      float pfmet_ex = analysisTree.pfmetcorr_ex;
      float pfmet_ey = analysisTree.pfmetcorr_ey;
      float pfmet_phi = analysisTree.pfmetcorr_phi;
      if (!isData) {
	pfmet_ex = analysisTree.pfmet_ex;
	pfmet_ey = analysisTree.pfmet_ey;
	pfmet_phi = analysisTree.pfmet_phi;
      }
      float pfmet = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      float puppimet_ex = analysisTree.puppimet_ex;
      float puppimet_ey = analysisTree.puppimet_ey;
      float puppimet_phi = analysisTree.puppimet_phi;
      float puppimet = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);
      //      std::cout << "Mets..." << std::endl;

      // vertex cuts
      //      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      //      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      //      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
      //		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      //      if (dVertex>dVertexCut) continue;
     
      // electron selection
      vector<unsigned int> allElectrons; allElectrons.clear();
      vector<unsigned int> isoIdElectrons; isoIdElectrons.clear();
      vector<bool> allElectronsIsLeg23Matched; allElectronsIsLeg23Matched.clear();
      vector<bool> allElectronsIsLeg17Matched; allElectronsIsLeg17Matched.clear();
      vector<bool> allElectronsIsLeg12Matched; allElectronsIsLeg12Matched.clear();
      vector<bool> allElectronsIsLegSingleEleMatched; allElectronsIsLegSingleEleMatched.clear();


      vector<bool> allElectronsIsTriggerMatched; allElectronsIsTriggerMatched.clear();
      vector<bool> isoIdElectronsIsTriggerMatched; isoIdElectronsIsTriggerMatched.clear();

      vector<float> allElectronsIso; allElectronsIso.clear();
      vector<float> isoIdElectronsIso; isoIdElectronsIso.clear();

      vector<bool> isElectronPassedIdIso; isElectronPassedIdIso.clear();

      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {

	float absEta = fabs(analysisTree.electron_eta[im]);
	float momScale = eleMomScaleBarrel;
	if (absEta>1.48) momScale = eleMomScaleEndcap;

	if (!isData) {
	  analysisTree.electron_px[im] *= momScale;
	  analysisTree.electron_py[im] *= momScale;
	  analysisTree.electron_pz[im] *= momScale;
	  analysisTree.electron_pt[im] *= momScale;
	}

	// selecting sample of probes
	if (analysisTree.electron_pt[im]<ptBins[0]) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronLooseCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElectronLooseCut) continue;

	bool electronTriggerMatch = false;
	bool electronEle23Match = false;
	bool electronEle17Match = false;
	bool electronEle12Match = false;
	bool electronSingleEleMatch = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.electron_eta[im],analysisTree.electron_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nElectronFilter] &&
	      analysisTree.trigobject_pt[iT]>eleTriggerPtCut&&
	      fabs(analysisTree.trigobject_eta[iT])<eleTriggerEtaCut)  // Electron Leg of electron trigger
	    electronTriggerMatch = true;
	  if (analysisTree.trigobject_filters[iT][nEle23Filter])
	    electronEle23Match = true;
	  if (analysisTree.trigobject_filters[iT][nEle17Filter])
	    electronEle17Match = true;
	  if (analysisTree.trigobject_filters[iT][nEle12Filter])
	    electronEle12Match = true;
	  if (analysisTree.trigobject_filters[iT][nSingleEleFilter]&&
	      analysisTree.trigobject_pt[iT]>singleEleTriggerPtCut&&
	      fabs(analysisTree.trigobject_eta[iT])<singleEleTriggerEtaCut)
	    electronSingleEleMatch = true;
	    
	}
	if (!applyTrigger) {
	  electronTriggerMatch = true;
	  electronEle23Match = true;
	  electronEle17Match = true;
	  electronEle12Match = true;
	  electronSingleEleMatch = true;
	}
	
        allElectrons.push_back(im);
	allElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	allElectronsIsLeg23Matched.push_back(electronEle23Match);
	allElectronsIsLeg17Matched.push_back(electronEle17Match);
	allElectronsIsLeg12Matched.push_back(electronEle12Match);
	allElectronsIsLegSingleEleMatched.push_back(electronSingleEleMatch);
	
	bool isPassed = true;

	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) isPassed = false;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut)  isPassed = false;
	// isolation
	float absIso = 0; 
	if(isoDR03) {
	  absIso = analysisTree.electron_r03_sumChargedHadronPt[im];
	  float neutralIso = 
	  analysisTree.electron_r03_sumNeutralHadronEt[im] + 
	  analysisTree.electron_r03_sumPhotonEt[im] - 
	  0.5*analysisTree.electron_r03_sumPUPt[im];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  absIso += neutralIso;
	}
	else {
	  absIso = analysisTree.electron_chargedHadIso[im];
          float neutralIso = analysisTree.electron_neutralHadIso[im] +
            analysisTree.electron_photonIso[im] -
            0.5*analysisTree.electron_puIso[im];
          neutralIso = TMath::Max(float(0),neutralIso);
          absIso += neutralIso;
	}
	float relIso = absIso/analysisTree.electron_pt[im];
	allElectronsIso.push_back(relIso);
	bool electronId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[im];
	if (electronIdType==1)
	  electronId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[im];
	else if (electronIdType==2)
	  electronId = analysisTree.electron_mva_wp80_trig_Spring15_v1[im];
	else if (electronIdType==3)
	  electronId = analysisTree.electron_mva_wp90_trig_Spring15_v1[im];
	else if (electronIdType==4)
	  electronId = analysisTree.electron_cutId_loose_Spring15[im];
	else if (electronIdType==5)
	  electronId = analysisTree.electron_cutId_medium_Spring15[im];
	else if (electronIdType==6)
	  electronId = analysisTree.electron_cutId_tight_Spring15[im];
	else if (electronIdType==7)
	  electronId = analysisTree.electron_cutId_veto_Spring15[im];
	if (!analysisTree.electron_pass_conversion[im]) electronId = false;
	if (analysisTree.electron_nmissinginnerhits[im]>1) electronId = false;
	bool isPassedProbe = isPassed && electronId && relIso<isoElectronProbeCut;
	isElectronPassedIdIso.push_back(isPassedProbe);
	isPassed = isPassed && electronId && relIso<isoElectronCut;
	if (isPassed) { 
	  isoIdElectrons.push_back(im);
	  isoIdElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	  isoIdElectronsIso.push_back(relIso);
	}
      }

      // end of electron selection

      //      std::cout << "Electrons..." << std::endl;

      // Monte Carlo analysis
      vector<unsigned int> promptElectrons; promptElectrons.clear();
      vector<unsigned int> nonpromptElectrons; nonpromptElectrons.clear();
      if (!isData) {
	for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	  if (fabs(analysisTree.genparticles_pdgid[igen])==11&&analysisTree.genparticles_status[igen]==1) {
	    if (analysisTree.genparticles_info[igen]==1)
	      promptElectrons.push_back(igen);
	    if (analysisTree.genparticles_info[igen]==5)
	      nonpromptElectrons.push_back(igen);
	  }
	}
      }

      vector<bool> isElectronPrompt; isElectronPrompt.clear();
      for (unsigned int iRecoElectrons=0; iRecoElectrons<allElectrons.size(); iRecoElectrons++) {
	unsigned int irec = allElectrons[iRecoElectrons];
	bool isPrompt = false;
	TLorentzVector recElectron; recElectron.SetXYZM(analysisTree.electron_px[irec],
							analysisTree.electron_py[irec],
							analysisTree.electron_pz[irec],
							electronMass);

	for (unsigned int iElectrons=0; iElectrons<promptElectrons.size(); ++iElectrons) {
	  unsigned int igen = promptElectrons[iElectrons];
	  TLorentzVector genElectron; genElectron.SetXYZM(analysisTree.genparticles_px[igen],
							  analysisTree.genparticles_py[igen],
							  analysisTree.genparticles_pz[igen],
							  electronMass);
	  

	  float relativeDifference = (genElectron-recElectron).P()/genElectron.P();
	  if (relativeDifference<0.05) {
	    unsigned int iEta = 0;
	    isPrompt = true;
	    if (TMath::Abs(genElectron.Eta())<1.48) 
	      iEta = 1;
	    if (isElectronPassedIdIso[iRecoElectrons]) 
	      PromptPtPass[iEta]->Fill(genElectron.Pt(),weight);
	    else
	      PromptPtFail[iEta]->Fill(genElectron.Pt(),weight);
	  }
	}
	isElectronPrompt.push_back(isPrompt);
      }

      vector<bool> isElectronNonPrompt; isElectronNonPrompt.clear();
      for (unsigned int iRecoElectrons=0; iRecoElectrons<allElectrons.size(); iRecoElectrons++) {
	unsigned int irec = allElectrons[iRecoElectrons];
	bool isNonPrompt = false;
	TLorentzVector recElectron; recElectron.SetXYZM(analysisTree.electron_px[irec],
							analysisTree.electron_py[irec],
							analysisTree.electron_pz[irec],
							electronMass);

	for (unsigned int iElectrons=0; iElectrons<nonpromptElectrons.size(); ++iElectrons) {
	  unsigned int igen = nonpromptElectrons[iElectrons];
	  TLorentzVector genElectron; genElectron.SetXYZM(analysisTree.genparticles_px[igen],
							  analysisTree.genparticles_py[igen],
							  analysisTree.genparticles_pz[igen],
							  electronMass);
	  
	  float relativeDifference = (genElectron-recElectron).P()/genElectron.P();
	  if (relativeDifference<0.05) {
	    unsigned int iEta = 0;
	    isNonPrompt = true;
	    if (TMath::Abs(genElectron.Eta())<1.48) 
	      iEta = 1;
	    if (isElectronPassedIdIso[iRecoElectrons]) 
	      NonPromptPtPass[iEta]->Fill(genElectron.Pt(),weight);
	    else
	      NonPromptPtFail[iEta]->Fill(genElectron.Pt(),weight);
	  }
	}
	isElectronNonPrompt.push_back(isNonPrompt);
      }

      bool isElectronsPair = false;
      bool firstTrigger = true;
      if (isoIdElectrons.size()>0) {
	unsigned int iE1 = 0;
	unsigned int iE2 = 0;
	bool isPairFound = false;
	float minIso = 999999;
	for (unsigned int im1=0; im1<isoIdElectrons.size(); ++im1) {
	  unsigned int index1 = isoIdElectrons[im1];
	  bool isTriggerMatched = isoIdElectronsIsTriggerMatched[im1];
	  if (isTriggerMatched && analysisTree.electron_pt[index1]>ptElectronHighCut && fabs(analysisTree.electron_eta[index1])<etaElectronHighCut) {
	    for (unsigned int iE=0; iE<allElectrons.size(); ++iE) {
	      unsigned int indexProbe = allElectrons[iE];
	      if (index1==indexProbe) continue;
	      float q1 = analysisTree.electron_charge[index1];
	      float q2 = analysisTree.electron_charge[indexProbe];
	      if (q1*q2>0) continue;
	      float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				analysisTree.electron_eta[indexProbe],analysisTree.electron_phi[indexProbe]);
	      if (dR<dRleptonsCut) continue;
	      float dPhi = dPhiFrom2P(analysisTree.electron_px[index1],analysisTree.electron_py[index1],
				analysisTree.electron_px[indexProbe],analysisTree.electron_py[indexProbe]);
	      if (dPhi>dPhileptonsCut) continue;
	      float ptProbe = TMath::Min(float(analysisTree.electron_pt[indexProbe]),float(ptBins[nPtBins]-0.01));
	      float absEtaProbe = TMath::Min(fabs(analysisTree.electron_eta[indexProbe]),float(etaBins[nEtaBins]-0.01));
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      int ptBinTrig = binNumber(ptProbe,nPtBinsTrig,ptBinsTrig);

	      if (ptBin<0) continue;
	      if (etaBin<0) continue;
	      if (ptBinTrig<0) continue;

	      TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[index1],
							  analysisTree.electron_py[index1],
							  analysisTree.electron_pz[index1],
							  electronMass);
	      TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[indexProbe],
							  analysisTree.electron_py[indexProbe],
							  analysisTree.electron_pz[indexProbe],
							  electronMass);

	      // number of jets
	      int nJets30 = 0;
	      int nJets30etaCut = 0;
	
	      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
		float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	  
		float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.electron_eta[index1],analysisTree.electron_phi[index1]);
		if (dR1<dRJetLeptonCut) continue;
		
		float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.electron_eta[indexProbe],analysisTree.electron_phi[indexProbe]);		
		if (dR2<dRJetLeptonCut) continue;
	  
		// pfJetId
		bool isPFJetId = looseJetiD(analysisTree,int(jet));
		if (!isPFJetId) continue;

		if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
		  nJets30++;
		  if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		    nJets30etaCut++;
		  }
		}
	      }	

	      int JetBin = nJets30etaCut;
	      if (JetBin>2) JetBin = 2;

	      float mass = (electron1+electron2).M();
	      //	      cout << "probe electron : eta = " << analysisTree.electron_eta[indexProbe]
	      //		   << "   pt = " << analysisTree.electron_eta[indexProbe] << endl;
	      if (isElectronPassedIdIso[iE]) {
		ZMassPass->Fill(mass,weight);
		ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtPass[etaBin][JetBin][ptBin]->Fill(mass,weight);
		if (isElectronPrompt[iE]) {
		  PromptSelPtPass[etaBin]->Fill(ptProbe,weight); 
		  PromptSelJetPtPass[etaBin][JetBin]->Fill(ptProbe,weight);
		}
		if (isElectronNonPrompt[iE]) {
		  NonPromptSelPtPass[etaBin]->Fill(ptProbe,weight);
		  NonPromptSelJetPtPass[etaBin][JetBin]->Fill(ptProbe,weight);
		}

		// ele23 filter
		if (allElectronsIsLeg23Matched[iE])
		  ZMassEle23EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		else 
		  ZMassEle23EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		// ele17 filter
		if (allElectronsIsLeg17Matched[iE])
		  ZMassEle17EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		else 
		  ZMassEle17EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		// ele12 filter
		if (allElectronsIsLeg12Matched[iE])
		  ZMassEle12EtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		else 
		  ZMassEle12EtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		// single ele filter
		if (absEtaProbe<singleEleTriggerEtaCut) {
		  if (allElectronsIsLegSingleEleMatched[iE])
		    ZMassIsoEleEtaPtPass[etaBin][ptBinTrig]->Fill(mass,weight);
		  else 
		    ZMassIsoEleEtaPtFail[etaBin][ptBinTrig]->Fill(mass,weight);
		}
	      }
	      else {
		ZMassFail->Fill(mass,weight);
		ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtFail[etaBin][JetBin][ptBin]->Fill(mass,weight);
		if (isElectronPrompt[iE]) {
		  PromptSelPtFail[etaBin]->Fill(ptProbe,weight); 
		  PromptSelJetPtFail[etaBin][JetBin]->Fill(ptProbe,weight);
		}
		if (isElectronNonPrompt[iE]) {
		  NonPromptSelPtFail[etaBin]->Fill(ptProbe,weight);
		  NonPromptSelJetPtFail[etaBin][JetBin]->Fill(ptProbe,weight);
		}
	      }
	    }
	  }
	  for (unsigned int im2=im1+1; im2<isoIdElectrons.size(); ++im2) {
	    unsigned int index2 = isoIdElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    bool isTriggerMatched1 = isoIdElectronsIsTriggerMatched[im1] && 
	      analysisTree.electron_pt[index1]>ptElectronHighCut&&
	      fabs(analysisTree.electron_eta[index1])<etaElectronHighCut;
	    bool isTriggerMatched2 = isoIdElectronsIsTriggerMatched[im2]&& 
	      analysisTree.electron_pt[index2]>ptElectronHighCut &&
	      fabs(analysisTree.electron_eta[index2])<etaElectronHighCut; 
	    /*	    bool isTriggerMatched = 
	      (isTriggerMatched1 
	       && analysisTree.electron_pt[index2]>ptElectronLowCut 
	       && fabs(analysisTree.electron_eta[index2])<etaElectronCut) ||
	      (isTriggerMatched2 
	       && analysisTree.electron_pt[index1]>ptElectronLowCut
	       && fabs(analysisTree.electron_eta[index1])<etaElectronCut);
	    */
	    bool isTriggerMatched = false;
	    if (analysisTree.electron_pt[index1]>analysisTree.electron_pt[index2]) {
	      isTriggerMatched = isTriggerMatched1;
	      firstTrigger = true;
	    }
	    else {
	      isTriggerMatched = isTriggerMatched2;
              firstTrigger = true;
	    }
	    bool sign = q1*q2>0;
	    float deltaREE = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],
				    analysisTree.electron_eta[index2],analysisTree.electron_phi[index2]);
	    if (oppositeSign)
	      sign = q1*q2 < 0;
	    if (sign && isTriggerMatched && deltaREE>dRleptonsCut) {
	      float isoSum = isoIdElectronsIso[im1] + isoIdElectronsIso[im2];
	      if (isoSum<minIso) {
		minIso = isoSum;
		isPairFound = true;
		iE1 = index1;
		iE2 = index2;
		isElectronsPair = true;
	      }
	    }
	  }
	}

	//	std::cout << "Pair found : " << isPairFound << std::endl;
	
	if (isPairFound) {
	  
	  TLorentzVector electron1; electron1.SetXYZM(analysisTree.electron_px[iE1],
						      analysisTree.electron_py[iE1],
						      analysisTree.electron_pz[iE1],
						      electronMass);
	  TLorentzVector electron2; electron2.SetXYZM(analysisTree.electron_px[iE2],
						      analysisTree.electron_py[iE2],
						      analysisTree.electron_pz[iE2],
						      electronMass);
	  TLorentzVector dielectron = electron1+electron2;
	  float mass = dielectron.M();

	  float mvamet = analysisTree.pfmetcorr_pt;
	  float mvamet_phi = analysisTree.pfmetcorr_phi;
	  float mvamet_ex = analysisTree.pfmetcorr_ex;
	  float mvamet_ey = analysisTree.pfmetcorr_ey;
	  
	  // applying electron scale factors
	  if (!isData&&applyLeptonSF) {
	    // insert code for leptons SF here
	    // SF1 - scale factor for first electron index = iE1
	    // pt1 = analysisTree.electron_pt[iE1];
	    // eta2 = analysisTree.electron_eta[iE1];
	    // SF2 - scale factor for the second electron index = iE2
	    // pt2 = analysisTree.electron_pt[iE2];
	    // eta2 = analysisTree.electron_eta[iE2];
	    // weight = weight*SF1*SF2
	    double ptEle1 = (double)analysisTree.electron_pt[iE1];
	    double ptEle2 = (double)analysisTree.electron_pt[iE2];
	    double etaEle1 = (double)analysisTree.electron_eta[iE1];
	    double etaEle2 = (double)analysisTree.electron_eta[iE2];
	    double IdIsoSF_ele1 = SF_electronIdIso->get_ScaleFactor(ptEle1, etaEle1);
	    double IdIsoSF_ele2 = SF_electronIdIso->get_ScaleFactor(ptEle2, etaEle2);
	
	    EleSF_IdIso_Ele1H->Fill(IdIsoSF_ele1);
	    EleSF_IdIso_Ele2H->Fill(IdIsoSF_ele2);
	    
	    double trkSF1 = trackEffEleH->GetBinContent(trackEffEleH->FindBin(etaEle1));
	    double trkSF2 = trackEffEleH->GetBinContent(trackEffEleH->FindBin(etaEle2));
	    

	    //	    if (ptEle1<20||ptEle2<20) {
	    //	      std::cout << "ele 1 ->  pt = " << ptEle1 << "   eta = " << etaEle1 << std::endl;
	    //	      std::cout << "eff data ele 1 = " << SF_electronIdIso->get_EfficiencyData(ptEle1, etaEle1)<< " |  eff mc ele 1 = " << SF_electronIdIso->get_EfficiencyMC(ptEle1, etaEle1)<<std::endl;
	    //	      std::cout << "ele 2 ->  pt = " << ptEle2 << "   eta = " << etaEle2 << std::endl;
	    //	      std::cout << "eff data ele 2 = " << SF_electronIdIso->get_EfficiencyData(ptEle2, etaEle2)<< " |  eff mc ele 2 = " << SF_electronIdIso->get_EfficiencyMC(ptEle2, etaEle2)<<std::endl;
	    //	      std::cout << "SF ele1 = " << IdIsoSF_ele1 << std::endl;
	    //	      std::cout << "SF ele2 = " << IdIsoSF_ele2 << std::endl;
	    //	      
	    //	      std::cout << " mass = " << mass << std::endl;
	    //	      std::cout << std::endl;
	    //	    }
	    weight = weight*IdIsoSF_ele1*IdIsoSF_ele2*trkSF1*trkSF2;

	    double effDataTrig1 = SF_electronTrig->get_EfficiencyData(ptEle1, etaEle1);  
	    double effDataTrig2 = SF_electronTrig->get_EfficiencyData(ptEle2, etaEle2);  
	    //	    double effTrigData = 1 - (1-effDataTrig1)*(1-effDataTrig2);
	    double effTrigData = 1;
	    if (firstTrigger) 
	      effTrigData = effDataTrig1;
	    else
	      effTrigData = effDataTrig2;

	    if (applyTrigger) {

	      double effMcTrig1 = SF_electronTrig->get_EfficiencyMC(ptEle1, etaEle1);  
	      double effMcTrig2 = SF_electronTrig->get_EfficiencyMC(ptEle2, etaEle2);  
	      double effMcTrig = 1 - (1-effMcTrig1)*(1-effMcTrig2);
	    
	      if (effTrigData>0&&effMcTrig>0) {
		double weightTrig = effTrigData/effMcTrig;
		// std::cout << "ele 1 ->  pt = " << ptEle1 << "   eta = " << etaEle1 << std::endl;
		// std::cout << "ele 2 ->  pt = " << ptEle2 << "   eta = " << etaEle2 << std::endl;
		// std::cout << "WeightTrig = " << weightTrig << std::endl;
		weight = weight*weightTrig;
	      }
	    }
	    else {
	      weight = weight*effTrigData;
	    }
	    
	  }
	   
	  // selecting good jets --->
	  
	  // HT variables
	  float HT30 = 0;
	  float HT20 = 0;
	  float HT30etaCut = 0;
	  float HT20etaCut = 0;
	  
	  // number of jets
	  int nJets30 = 0;
	  int nJets20 = 0;
	  int nJets30etaCut = 0;
	  int nJets20etaCut = 0;

	  int nJets30Up = 0;
	  int nJets30Down = 0;


	  for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	    float scale = 1;
	    if (applyJES>0) scale = 1 + analysisTree.pfjet_jecUncertainty[jet];
	    if (applyJES<0) scale = 1 - analysisTree.pfjet_jecUncertainty[jet];
	    if (applyJES!=0&&!isData) {
	      analysisTree.pfjet_pt[jet] *= scale;
	      analysisTree.pfjet_px[jet] *= scale;
	      analysisTree.pfjet_py[jet] *= scale;
	      analysisTree.pfjet_pz[jet] *= scale;
	      analysisTree.pfjet_e[jet] *= scale;
	    }
	    float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	    if (absJetEta>jetEtaCut) continue;
	    
	    float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			       analysisTree.electron_eta[iE1],analysisTree.electron_phi[iE1]);
	    if (dR1<dRJetLeptonCut) continue;
	    
	    float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			       analysisTree.electron_eta[iE2],analysisTree.electron_phi[iE2]);
	    if (dR2<dRJetLeptonCut) continue;
	  
	    // pfJetId
	    bool isPFJetId = looseJetiD(analysisTree,int(jet));
	    if (!isPFJetId) continue;
	    
	    if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
	      nJets30++;
	      HT30 += analysisTree.pfjet_pt[jet];
	      if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		HT30etaCut += analysisTree.pfjet_pt[jet];
		nJets30etaCut++;
	      }
	    }
	    float ptJetUp = analysisTree.pfjet_pt[jet]*(1 + analysisTree.pfjet_jecUncertainty[jet]);
	    float ptJetDown = analysisTree.pfjet_pt[jet]*(1 - analysisTree.pfjet_jecUncertainty[jet]);
	    if (ptJetUp>jetPtHighCut)
	      nJets30Up++;
	    if (ptJetDown>jetPtHighCut)
	      nJets30Down++;
	    
	    if (analysisTree.pfjet_pt[jet]>jetPtLowCut) {
	      nJets20++;
	      HT20 += analysisTree.pfjet_pt[jet]; 
	      if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		HT20etaCut += analysisTree.pfjet_pt[jet];
		nJets20etaCut++;
	      }
	    }	  
	    

	  }

	  if (!isData&&applyNJetReweighting) {
	    float njetWeight = nJetsWeight(nJets30);
	    weight *= njetWeight;
	  }

	  float visiblePx = dielectron.Px();
	  float visiblePy = dielectron.Py();
	  if (applyRecoilOnGenerator) {
	    visiblePx = genL.Px();
	    visiblePy = genL.Py();
	  }

	  if (!isData && applyRecoilCorrections) {
	    
	    //	    std::cout << "applying recoil corrections " << std::endl;

	    float pfmetcorr_ex = pfmet_ex;
	    float pfmetcorr_ey = pfmet_ey;
	    if (applySimpleRecoilCorrections)
	      recoilPFMetCorrector.CorrectByMeanResolution(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);
	    else 
	      recoilPFMetCorrector.Correct(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);
	    pfmet_phi = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
	    pfmet = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
	    pfmet_ex = pfmetcorr_ex;
	    pfmet_ey = pfmetcorr_ey;
	    
	    float puppimetcorr_ex = puppimet_ex;
	    float puppimetcorr_ey = puppimet_ey;
	    if (applySimpleRecoilCorrections)
	      recoilPuppiMetCorrector.CorrectByMeanResolution(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,puppimetcorr_ex,puppimetcorr_ey);
	    else
	      recoilPuppiMetCorrector.Correct(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,puppimetcorr_ex,puppimetcorr_ey);
	    puppimet_phi = TMath::ATan2(puppimetcorr_ey,puppimetcorr_ex);
	    puppimet = TMath::Sqrt(puppimetcorr_ex*puppimetcorr_ex+puppimetcorr_ey*puppimetcorr_ey);
	    puppimet_ex = puppimetcorr_ex;
	    puppimet_ey = puppimetcorr_ey;

	    float mvametcorr_ex = mvamet_ex;
	    float mvametcorr_ey = mvamet_ey;
	    if (applySimpleRecoilCorrections)
	      recoilMvaMetCorrector.CorrectByMeanResolution(mvamet_ex,mvamet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,mvametcorr_ex,mvametcorr_ey);
	    else
	      recoilMvaMetCorrector.Correct(mvamet_ex,mvamet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,mvametcorr_ex,mvametcorr_ey);
	    mvamet_phi = TMath::ATan2(mvametcorr_ey,mvametcorr_ex);
	    mvamet = TMath::Sqrt(mvametcorr_ex*mvametcorr_ex+mvametcorr_ey*mvametcorr_ey);
	    mvamet_ex = mvametcorr_ex;
	    mvamet_ey = mvametcorr_ey;

	  }

	    
          float ptLeadingE = analysisTree.electron_pt[iE1];
	  float etaLeadingE = analysisTree.electron_eta[iE1];
	  float phiLeadingE = analysisTree.electron_phi[iE1];

	  float ptTrailingE = analysisTree.electron_pt[iE2];
	  float etaTrailingE = analysisTree.electron_eta[iE2];
	  float phiTrailingE = analysisTree.electron_phi[iE2];

	  if (ptTrailingE>ptLeadingE) {
	    float temp = ptLeadingE;
	    ptLeadingE = ptTrailingE;
	    ptTrailingE = temp;
	    temp = etaLeadingE;
	    etaLeadingE = etaTrailingE;
	    etaTrailingE = temp;
	    temp = phiLeadingE;
	    phiLeadingE = phiTrailingE;
	    phiTrailingE = temp;
	    float itemp = iE1;
	    iE1 = iE2;
	    iE2 = itemp;
	  } 

	  if (mass>dielectronMassCut) {

	    massSelH->Fill(mass,weight);
	    massExtendedSelH->Fill(mass,weight);
	    for (int iScale=0; iScale<21; ++ iScale) {
	      float scaleFactor = 0.98 + 0.002*float(iScale);
	      massSelScaleH[iScale]->Fill(mass*scaleFactor,weight);
	    }

	    nJets30SelH->Fill(nJets30,weight);
	    nJets20SelH->Fill(nJets20,weight);
	    nJets30etaCutSelH->Fill(nJets30etaCut,weight);
	    nJets20etaCutSelH->Fill(nJets20etaCut,weight);
	    ptLeadingEleSelH->Fill(ptLeadingE,weight);
	    ptTrailingEleSelH->Fill(ptTrailingE,weight);
	    etaLeadingEleSelH->Fill(etaLeadingE,weight);
	    etaTrailingEleSelH->Fill(etaTrailingE,weight);
	    dielectronPtSelH->Fill(dielectron.Pt(),weight);
	    dielectronEtaSelH->Fill(dielectron.Eta(),weight);
	    NumberOfVerticesH->Fill(float(analysisTree.primvertex_count),weight);
	    metSelH->Fill(pfmet,weight);
	    puppimetSelH->Fill(puppimet,weight);
	    mvametSelH->Fill(mvamet,weight);

	    // ************************
	    // ** Filling MET NTuple **
	    // ************************

	    evt = analysisTree.event_nr;
	    pt_1 = ptLeadingE;
	    eta_1 = etaLeadingE;
	    phi_1 = phiLeadingE;
	    match_1 = false;
	    pt_2 = ptTrailingE;
	    eta_2 = etaTrailingE;
	    phi_2 = phiTrailingE;
	    match_2 = false;
	    njets = nJets30;
	    mvametX = mvamet;
	    metX = pfmet;
	    metuncorrX = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex+analysisTree.pfmet_ey*analysisTree.pfmet_ey);
	    puppimetX = puppimet;
	    if (genElectrons.size()==2) {
	      for (unsigned int iGen=0; iGen<genElectrons.size(); ++iGen) {
		unsigned int igen = genElectrons[iGen];
		TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
						    analysisTree.genparticles_py[igen],
						    analysisTree.genparticles_pz[igen],
						    analysisTree.genparticles_e[igen]);
		float deltaRLeading = deltaR(etaLeadingE,phiLeadingE,
					     genLV.Eta(),genLV.Phi());
		float deltaRTrailing = deltaR(etaTrailingE,phiTrailingE,
					     genLV.Eta(),genLV.Phi());
		if (deltaRLeading<0.2)
		  match_1 = true;
		if (deltaRTrailing<0.2)
		  match_2 = true;

	      }
	    }

	    float dielectronPt = dielectron.Pt();
	    float dielectronEta = dielectron.Eta();
	    float dielectronPhi = dielectron.Phi();

	    dilepton_pt = dielectronPt;
	    dilepton_eta = dielectronEta;
	    dilepton_phi = dielectronPhi;

	    float unitX = dielectron.Px()/dielectron.Pt();
	    float unitY = dielectron.Py()/dielectron.Pt();
	    float phiUnit = TMath::ATan2(unitY,unitX);
	    float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
	    float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
	    
	    float pfmetcorr_ex = analysisTree.pfmetcorr_ex;
	    float pfmetcorr_ey = analysisTree.pfmetcorr_ey;

	    float u1 = 0;
	    float u2  = 0;
	    float rH = 0;

	    computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,u1,u2,rH);
	    u1_metuncorrX = u1;
	    u2_metuncorrX = u2;

	    computeRecoil(pfmetcorr_ex,pfmetcorr_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,u1,u2,rH);
            u1_metX = u1;
            u2_metX = u2;

	    computeRecoil(mvamet_ex,mvamet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,u1,u2,rH);
            u1_mvametX = u1;
            u2_mvametX = u2;

	    synch_tree->Fill();

	    // jet pt bin
	    int jetBin = 0;
	    if (nJets30==1)
	      jetBin = 1;
	    else if (nJets30>1)
	      jetBin = 2;

	    int jetBinUp = 0;
	    if (nJets30Up==1)
	      jetBinUp = 1;
	    if (nJets30Up>1)
	      jetBinUp = 2;
	    
	    int jetBinDown = 0;
	    if (nJets30Down==1)
	      jetBinDown = 1;
	    if (nJets30Down>1)
	      jetBinDown = 2;

	    metSelNJets[jetBin]->Fill(pfmet,weight);
	    puppimetSelNJets[jetBin]->Fill(puppimet,weight);
	    mvametSelNJets[jetBin]->Fill(mvamet,weight);

	    for (int iScale=0; iScale<21; ++ iScale) {
	      float scaleFactor = 0.9 + 0.01*float(iScale);
	      metSelScaleH[iScale]->Fill(pfmet*scaleFactor,weight);
	      puppimetSelScaleH[iScale]->Fill(puppimet*scaleFactor,weight);
	      mvametSelScaleH[iScale]->Fill(mvamet*scaleFactor,weight);
	    }


	    
	    // systematic shifts in mva met --->

	    float mvamet_ex_ScaleUp   = mvamet_ex;
	    float mvamet_ex_ScaleDown = mvamet_ex;
	    float mvamet_ex_ResoUp    = mvamet_ex;
	    float mvamet_ex_ResoDown  = mvamet_ex;
	    float mvamet_ex_ElectronUp    = mvamet_ex - 0.01*dielectron.Px();
	    float mvamet_ex_ElectronDown  = mvamet_ex + 0.01*dielectron.Px();
	    
	    float mvamet_ey_ScaleUp   = mvamet_ey;
	    float mvamet_ey_ScaleDown = mvamet_ey;
	    float mvamet_ey_ResoUp    = mvamet_ey;
	    float mvamet_ey_ResoDown  = mvamet_ey;
	    float mvamet_ey_ElectronUp    = mvamet_ey - 0.01*dielectron.Py();
	    float mvamet_ey_ElectronDown  = mvamet_ey + 0.01*dielectron.Py();

	    float dielectronPtUp   = 1.01*dielectronPt;
	    float dielectronPtDown = 0.99*dielectronPt;
	

	    metSysMva.ShiftMEt(mvamet_ex,mvamet_ey,
			       genV.Px(),genV.Py(),
			       visiblePx,visiblePy,
			       nJets30,
			       bkgdType,
			       MEtSys::SysType::Response,
			       1.10,
			       mvamet_ex_ScaleUp,mvamet_ey_ScaleUp);

	    metSysMva.ShiftMEt(mvamet_ex,mvamet_ey,
			       genV.Px(),genV.Py(),
			       visiblePx,visiblePy,
			       nJets30,
			       bkgdType,
			       MEtSys::SysType::Response,
			       0.90,
			       mvamet_ex_ScaleDown,mvamet_ey_ScaleDown);
	
	    metSysMva.ShiftMEt(mvamet_ex,mvamet_ey,
			       genV.Px(),genV.Py(),
			       visiblePx,visiblePy,
			       nJets30,
			       bkgdType,
			       MEtSys::SysType::Resolution,
			       1.10,
			       mvamet_ex_ResoUp,mvamet_ey_ResoUp);
	    
	    metSysMva.ShiftMEt(mvamet_ex,mvamet_ey,
			       genV.Px(),genV.Py(),
			       visiblePx,visiblePy,
			       nJets30,
			       bkgdType,
			       MEtSys::SysType::Resolution,
			       0.90,
			       mvamet_ex_ResoDown,mvamet_ey_ResoDown);
	    
	    float mvamet_ScaleUp = TMath::Sqrt(mvamet_ex_ScaleUp*mvamet_ex_ScaleUp+
					       mvamet_ey_ScaleUp*mvamet_ey_ScaleUp);
	    
	    float mvamet_ScaleDown = TMath::Sqrt(mvamet_ex_ScaleDown*mvamet_ex_ScaleDown+
						 mvamet_ey_ScaleDown*mvamet_ey_ScaleDown);
	    
	    float mvamet_ResoUp = TMath::Sqrt(mvamet_ex_ResoUp*mvamet_ex_ResoUp+
					      mvamet_ey_ResoUp*mvamet_ey_ResoUp);
	    
	    float mvamet_ResoDown = TMath::Sqrt(mvamet_ex_ResoDown*mvamet_ex_ResoDown+
						mvamet_ey_ResoDown*mvamet_ey_ResoDown);
	    
	    float mvamet_ElectronUp = TMath::Sqrt(mvamet_ex_ElectronUp*mvamet_ex_ElectronUp+
						  mvamet_ey_ElectronUp*mvamet_ey_ElectronUp);
	    
	    float mvamet_ElectronDown = TMath::Sqrt(mvamet_ex_ElectronDown*mvamet_ex_ElectronDown+
						    mvamet_ey_ElectronDown*mvamet_ey_ElectronDown);
	    int ptBin = binNumber(TMath::Min(float(dielectronPt),float(999)),nZPtBins,zPtBins);

						
	  
	    float recoilParal = 0;
	    float recoilPerp  = 0;
	    float responseHad = 0;
	    computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilParal,recoilPerp,responseHad);

	    float recoilPuppiParal = 0;
	    float recoilPuppiPerp  = 0;
	    float responsePuppiHad = 0;
	    computeRecoil(puppimet_ex,puppimet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilPuppiParal,recoilPuppiPerp,responsePuppiHad);

	    float recoilMvaParal = 0;
	    float recoilMvaPerp  = 0;
	    float responseMvaHad = 0;
	    computeRecoil(mvamet_ex,mvamet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilMvaParal,recoilMvaPerp,responseMvaHad);

	    float recoilMvaScaleUpParal = 0;
	    float recoilMvaScaleUpPerp  = 0;
	    float responseMvaScaleUpHad = 0;
	    computeRecoil(mvamet_ex_ScaleUp,mvamet_ey_ScaleUp,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilMvaScaleUpParal,recoilMvaScaleUpPerp,responseMvaScaleUpHad);

	    float recoilMvaScaleDownParal = 0;
	    float recoilMvaScaleDownPerp  = 0;
	    float responseMvaScaleDownHad = 0;
	    computeRecoil(mvamet_ex_ScaleDown,mvamet_ey_ScaleDown,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilMvaScaleDownParal,recoilMvaScaleDownPerp,responseMvaScaleDownHad);

	    float recoilMvaResoUpParal = 0;
	    float recoilMvaResoUpPerp  = 0;
	    float responseMvaResoUpHad = 0;
	    computeRecoil(mvamet_ex_ResoUp,mvamet_ey_ResoUp,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilMvaResoUpParal,recoilMvaResoUpPerp,responseMvaResoUpHad);

	    float recoilMvaResoDownParal = 0;
	    float recoilMvaResoDownPerp  = 0;
	    float responseMvaResoDownHad = 0;
	    computeRecoil(mvamet_ex_ResoDown,mvamet_ey_ResoDown,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilMvaResoDownParal,recoilMvaResoDownPerp,responseMvaResoDownHad);

	    float recoilMvaElectronUpParal = 0;
	    float recoilMvaElectronUpPerp  = 0;
	    float responseMvaElectronUpHad = 0;
	    computeRecoil(mvamet_ex_ElectronUp,mvamet_ey_ElectronUp,unitX,unitY,perpUnitX,perpUnitY,dielectronPtUp,recoilMvaElectronUpParal,recoilMvaElectronUpPerp,responseMvaElectronUpHad);
	    
	    float recoilMvaElectronDownParal = 0;
	    float recoilMvaElectronDownPerp  = 0;
	    float responseMvaElectronDownHad = 0;
	    computeRecoil(mvamet_ex_ElectronDown,mvamet_ey_ElectronDown,unitX,unitY,perpUnitX,perpUnitY,dielectronPtDown,recoilMvaElectronDownParal,recoilMvaElectronDownPerp,responseMvaElectronDownHad);

	    if (mass>70&&mass<110) {
	      metZSelH->Fill(pfmet,weight);
	      puppimetZSelH->Fill(puppimet,weight);
	      mvametZSelH->Fill(mvamet,weight);
	      
	      metZSelNJets[jetBin]->Fill(pfmet,weight);
	      puppimetZSelNJets[jetBin]->Fill(puppimet,weight);
	      mvametZSelNJets[jetBin]->Fill(mvamet,weight);


	      // pfmet

	      recoilZParalH[jetBin]->Fill(recoilParal,weight);
	      recoilZPerpH[jetBin]->Fill(recoilPerp,weight);
	      recoilZParal_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilParal,weight);
	      recoilZPerp_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPerp,weight);
	      recoilResponse[jetBin]->Fill(responseHad,weight);
	      recoilResponse_Ptbins_nJetsH[jetBin][ptBin]->Fill(responseHad,weight);
	      
	      // puppimet
	      recoilPuppiZParalH[jetBin]->Fill(recoilPuppiParal,weight);
	      recoilPuppiZPerpH[jetBin]->Fill(recoilPuppiPerp,weight);
	      recoilPuppiZParal_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPuppiParal,weight);
	      recoilPuppiZPerp_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPuppiPerp,weight);
	      recoilPuppiResponse[jetBin]->Fill(responsePuppiHad,weight);
	      recoilPuppiResponse_Ptbins_nJetsH[jetBin][ptBin]->Fill(responsePuppiHad,weight);
	      
	      // mvamet
	      recoilMvaZParalH[jetBin]->Fill(recoilMvaParal,weight);
	      recoilMvaZPerpH[jetBin]->Fill(recoilMvaPerp,weight);
	      recoilMvaZParal_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilMvaParal,weight);
	      recoilMvaZPerp_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilMvaPerp,weight);
	      recoilMvaResponse[jetBin]->Fill(responseMvaHad,weight);
	      recoilMvaResponse_Ptbins_nJetsH[jetBin][ptBin]->Fill(responseMvaHad,weight);

	      // systematic changes

	      // JES
	      mvametZSelNJetsSysH[jetBinUp][0][0]->Fill(mvamet,weight);
	      mvametZSelNJetsSysH[jetBinDown][0][1]->Fill(mvamet,weight);
	      
	      recoilMvaZParalSysH[jetBinUp][0][0]->Fill(recoilMvaParal,weight);
	      recoilMvaZParalSysH[jetBinDown][0][1]->Fill(recoilMvaParal,weight);
	      
	      recoilMvaZPerpSysH[jetBinUp][0][0]->Fill(recoilMvaPerp,weight);
	      recoilMvaZPerpSysH[jetBinDown][0][1]->Fill(recoilMvaPerp,weight);
	      
	      // MEt scale 
	      mvametZSelNJetsSysH[jetBin][1][0]->Fill(mvamet_ScaleUp,weight);
	      mvametZSelNJetsSysH[jetBin][1][1]->Fill(mvamet_ScaleDown,weight);
	      
	      recoilMvaZParalSysH[jetBin][1][0]->Fill(recoilMvaScaleUpParal,weight);
	      recoilMvaZParalSysH[jetBin][1][1]->Fill(recoilMvaScaleDownParal,weight);

	      recoilMvaZPerpSysH[jetBin][1][0]->Fill(recoilMvaScaleUpPerp,weight);
	      recoilMvaZPerpSysH[jetBin][1][1]->Fill(recoilMvaScaleDownPerp,weight);

	      // MEt resolution 
	      mvametZSelNJetsSysH[jetBin][2][0]->Fill(mvamet_ResoUp,weight);
	      mvametZSelNJetsSysH[jetBin][2][1]->Fill(mvamet_ResoDown,weight);

	      recoilMvaZParalSysH[jetBin][2][0]->Fill(recoilMvaResoUpParal,weight);
	      recoilMvaZParalSysH[jetBin][2][1]->Fill(recoilMvaResoDownParal,weight);
	      
	      recoilMvaZPerpSysH[jetBin][2][0]->Fill(recoilMvaResoUpPerp,weight);
	      recoilMvaZPerpSysH[jetBin][2][1]->Fill(recoilMvaResoDownPerp,weight);
	      
	      // Electron momentum
	      mvametZSelNJetsSysH[jetBin][3][0]->Fill(mvamet_ElectronUp,weight);
	      mvametZSelNJetsSysH[jetBin][3][1]->Fill(mvamet_ElectronDown,weight);

	      recoilMvaZParalSysH[jetBin][3][0]->Fill(recoilMvaElectronUpParal,weight);
	      recoilMvaZParalSysH[jetBin][3][1]->Fill(recoilMvaElectronDownParal,weight);
	      
	      recoilMvaZPerpSysH[jetBin][3][0]->Fill(recoilMvaElectronUpPerp,weight);
	      recoilMvaZPerpSysH[jetBin][3][1]->Fill(recoilMvaElectronDownPerp,weight);
	      
	    }

	    if (mass<70||mass>110) { // Z-depleted region
	      
	      metTopSelNJets[jetBin]->Fill(pfmet,weight);
	      puppimetTopSelNJets[jetBin]->Fill(puppimet,weight);
	      mvametTopSelNJets[jetBin]->Fill(mvamet,weight);
	      
	      // pfmet
	      recoilTopParalH[jetBin]->Fill(recoilParal,weight);
	      recoilTopPerpH[jetBin]->Fill(recoilPerp,weight);
	      
	      // puppimet
	      recoilPuppiTopParalH[jetBin]->Fill(recoilPuppiParal,weight);
	      recoilPuppiTopPerpH[jetBin]->Fill(recoilPuppiPerp,weight);
	      
	      // mvamet
	      recoilMvaTopParalH[jetBin]->Fill(recoilMvaParal,weight);
	      recoilMvaTopPerpH[jetBin]->Fill(recoilMvaPerp,weight);
	      

	      // systematic changes
	      
	      // JES
	      mvametTopSelNJetsSysH[jetBinUp][0][0]->Fill(mvamet,weight);
	      mvametTopSelNJetsSysH[jetBinDown][0][1]->Fill(mvamet,weight);
	    
	      recoilMvaTopParalSysH[jetBinUp][0][0]->Fill(recoilMvaParal,weight);
	      recoilMvaTopParalSysH[jetBinDown][0][1]->Fill(recoilMvaParal,weight);

	      recoilMvaTopPerpSysH[jetBinUp][0][0]->Fill(recoilMvaPerp,weight);
	      recoilMvaTopPerpSysH[jetBinDown][0][1]->Fill(recoilMvaPerp,weight);

	      // MEt scale 
	      mvametTopSelNJetsSysH[jetBin][1][0]->Fill(mvamet_ScaleUp,weight);
	      mvametTopSelNJetsSysH[jetBin][1][1]->Fill(mvamet_ScaleDown,weight);

	      recoilMvaTopParalSysH[jetBin][1][0]->Fill(recoilMvaScaleUpParal,weight);
	      recoilMvaTopParalSysH[jetBin][1][1]->Fill(recoilMvaScaleDownParal,weight);

	      recoilMvaTopPerpSysH[jetBin][1][0]->Fill(recoilMvaScaleUpPerp,weight);
	      recoilMvaTopPerpSysH[jetBin][1][1]->Fill(recoilMvaScaleDownPerp,weight);

	      // MEt resolution 
	      mvametTopSelNJetsSysH[jetBin][2][0]->Fill(mvamet_ResoUp,weight);
	      mvametTopSelNJetsSysH[jetBin][2][1]->Fill(mvamet_ResoDown,weight);
	      
	      recoilMvaTopParalSysH[jetBin][2][0]->Fill(recoilMvaResoUpParal,weight);
	      recoilMvaTopParalSysH[jetBin][2][1]->Fill(recoilMvaResoDownParal,weight);

	      recoilMvaTopPerpSysH[jetBin][2][0]->Fill(recoilMvaResoUpPerp,weight);
	      recoilMvaTopPerpSysH[jetBin][2][1]->Fill(recoilMvaResoDownPerp,weight);

	      // Electron momentum
	      mvametTopSelNJetsSysH[jetBin][3][0]->Fill(mvamet_ElectronUp,weight);
	      mvametTopSelNJetsSysH[jetBin][3][1]->Fill(mvamet_ElectronDown,weight);

	      recoilMvaTopParalSysH[jetBin][3][0]->Fill(recoilMvaElectronUpParal,weight);
	      recoilMvaTopParalSysH[jetBin][3][1]->Fill(recoilMvaElectronDownParal,weight);
	      
	      recoilMvaTopPerpSysH[jetBin][3][0]->Fill(recoilMvaElectronUpPerp,weight);
	      recoilMvaTopPerpSysH[jetBin][3][1]->Fill(recoilMvaElectronDownPerp,weight);
	      
	    }
	  }
	}
      }
      
      if (isElectronsPair) selEventsDielectrons++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                              = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                            = " << nEvents << std::endl;
  std::cout << "Total number of selected events (electron pairs)          = " << selEventsDielectrons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs   :";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  
  file->Write();
  file->Close();
  delete file;
  
  
}
