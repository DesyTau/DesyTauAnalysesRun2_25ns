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

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"

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
  const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");

  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
  const bool applyRecoilCorrections =  cfg.get<bool>("ApplyRecoilCorrections");
  const bool applyRecoilOnGenerator = cfg.get<bool>("ApplyRecoilOnGenerator");
  const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
  const bool applyTopPtReweighting = cfg.get<bool>("ApplyTopPtReweighting");
  const bool applyZMassPtReweighting = cfg.get<bool>("ApplyZMassPtReweighting");
  const bool interpolateZMassPtWeight = cfg.get<bool>("InterpolateZMassPtWeight");

  // kinematic cuts on electron
  const float ptElectronCut       = cfg.get<float>("ptElectronCut");
  const float ptElectronProbeCut  = cfg.get<float>("ptElectronProbeCut");
  const float etaElectronCut      = cfg.get<float>("etaElectronCut");
  const float etaElectronProbeCut = cfg.get<float>("etaElectronProbeCut");
  const float dxyElectronCut      = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut       = cfg.get<float>("dzElectronCut");
  const float dxyElectronProbeCut = cfg.get<float>("dxyElectronProbeCut");
  const float dzElectronProbeCut  = cfg.get<float>("dzElectronProbeCut");
  const float isoElectronCut      = cfg.get<float>("isoElectronCut");
  const float isoElectronProbeCut = cfg.get<float>("isoElectronProbeCut");
  const unsigned int electronIdType = cfg.get<unsigned int>("ElectronIdType");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool  oppositeSign   = cfg.get<bool>("OppositeSign");
  const float dielectronMassCut = cfg.get<float>("DielectronMassCut");
  const bool isoDR03         = cfg.get<bool>("IsoDR03");
  const bool applyRhoCorrectedIso = cfg.get<bool>("ApplyRhoCorrectedIso");

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

  TString ElectronHLTName(electronHLTName);
  TString ElectronHLTFilterName(electronHLTFilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // electron SFs
  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");
  const string ElectronTrigFile  = cfg.get<string>("ElectronTrigEff"); 
  const string correctionWorkspaceFile = cfg.get<string>("CorrectionWorkspaceFile");
  TString CorrectionWorkspaceFile(correctionWorkspaceFile);

  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  // Z Mass pT reweighting
  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

  // json file
  const string jsonFile = cfg.get<string>("jsonFile");

  // corrections
  const int applyJES = cfg.get<int>("applyJES");
  const float eleMomScaleBarrel = cfg.get<float>("EleMomScaleBarrel");
  const float eleMomScaleEndcap = cfg.get<float>("EleMomScaleEndcap");

  // pileup 
  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");
  const string pileUpMCHist = cfg.get<string>("PileUpMCHistName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);
  TString PileUpMCHist(pileUpMCHist);

  float isoCone = 0.4;
  if (isoDR03) isoCone = 0.3;

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

  TH1D * nvertH = new TH1D("nvertH","",100,-0.5,99.5);
  TH1D * nvertSelH = new TH1D("nvertSelH","",100,-0.5,99.5);

  // Histograms after final selection
  TH1D * massSelH = new TH1D("massSelH","",200,0,200);
  TH1D * massExtendedSelH = new TH1D("massExtendedSelH","",500,0,5000);
  TH1D * ptLeadingEleSelH = new TH1D("ptLeadingEleSelH","",100,0,200);
  TH1D * ptTrailingEleSelH = new TH1D("ptTrailingEleSelH","",100,0,200);
  TH1D * etaLeadingEleSelH = new TH1D("etaLeadingEleSelH","",50,-2.5,2.5);
  TH1D * etaTrailingEleSelH = new TH1D("etaTrailingEleSelH","",50,-2.5,2.5);
  TH1D * dielectronPtSelH = new TH1D("dielectronPtSelH","",100,0,1000);
  TH1D * dielectronEtaSelH = new TH1D("dielectronEtaSelH","",120,-6,6);

  TH1D * metSelH = new TH1D("metSelH","",200,0.,400.);
  TH1D * puppimetSelH = new TH1D("puppimetSelH","",200,0.,400.);

  TH1D * metZSelH  = new TH1D("metZSelH","",200,0,400);
  TH1D * puppimetZSelH = new TH1D("puppimetZSelH","",200,0,400);

  TString scales[21] = {"M10","M9","M8","M7","M6","M5","M4","M3","M2","M1","0",
			"P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};

  TH1D * massSelScaleH[21];
  TH1D * metSelScaleH[21];
  TH1D * puppimetSelScaleH[21];
  for (int iScale=0; iScale<21; ++iScale) {    
    massSelScaleH[iScale] = new TH1D("massSel"+scales[iScale]+"H","",200,0,200);
    metSelScaleH[iScale] = new TH1D("metSel"+scales[iScale]+"H","",200,0,400);
    puppimetSelScaleH[iScale] = new TH1D("puppimetSel"+scales[iScale]+"H","",200,0,400);
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

  TString RecoilTopParal("recoilTopParal_");
  TString RecoilTopPerp("recoilTopPerp_");
  TString RecoilPuppiTopParal("recoilPuppiTopParal_");
  TString RecoilPuppiTopPerp("recoilPuppiTopPerp_");

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
  
  TH1D * metZSelNJets[3];
  TH1D * puppimetZSelNJets[3];

  TH1D * metTopSelNJets[3];
  TH1D * puppimetTopSelNJets[3];

  TH1D * recoilZParalH[3];
  TH1D * recoilZPerpH[3];
  TH1D * recoilPuppiZParalH[3];
  TH1D * recoilPuppiZPerpH[3];

  TH1D * recoilTopParalH[3];
  TH1D * recoilTopPerpH[3];
  TH1D * recoilPuppiTopParalH[3];
  TH1D * recoilPuppiTopPerpH[3];

  TH1D * recoilZParal_Ptbins_nJetsH[3][5];
  TH1D * recoilZPerp_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiZParal_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiZPerp_Ptbins_nJetsH[3][5];

  TH1D * recoilResponse[3];
  TH1D * recoilResponse_Ptbins_nJetsH[3][5];
  TH1D * recoilPuppiResponse[3];
  TH1D * recoilPuppiResponse_Ptbins_nJetsH[3][5];
  
  for (int iJets=0; iJets<nJetBins; ++iJets) {

    metZSelNJets[iJets] = new TH1D("metZSel_"+NJetBins[iJets],"",400,0,400);
    metTopSelNJets[iJets] = new TH1D("metTopSel_"+NJetBins[iJets],"",400,0,400);
    metSelNJets[iJets] = new TH1D("metSel_"+NJetBins[iJets],"",400,0,400);

    puppimetZSelNJets[iJets] = new TH1D("puppimetZSel_"+NJetBins[iJets],"",400,0,400);
    puppimetTopSelNJets[iJets] = new TH1D("puppimetTopSel_"+NJetBins[iJets],"",400,0,400);
    puppimetSelNJets[iJets] = new TH1D("puppimetSel_"+NJetBins[iJets],"",400,0,400);

    recoilZParalH[iJets] = new TH1D(RecoilZParal+NJetBins[iJets],"",200,-400,400);
    recoilZPerpH[iJets] = new TH1D(RecoilZPerp+NJetBins[iJets],"",200,-400,400);
    recoilPuppiZParalH[iJets] = new TH1D(RecoilPuppiZParal+NJetBins[iJets],"",200,-400,400);
    recoilPuppiZPerpH[iJets] = new TH1D(RecoilPuppiZPerp+NJetBins[iJets],"",200,-400,400);

    recoilTopParalH[iJets] = new TH1D(RecoilTopParal+NJetBins[iJets],"",200,-400,400);
    recoilTopPerpH[iJets] = new TH1D(RecoilTopPerp+NJetBins[iJets],"",200,-400,400);
    recoilPuppiTopParalH[iJets] = new TH1D(RecoilPuppiTopParal+NJetBins[iJets],"",200,-400,400);
    recoilPuppiTopPerpH[iJets] = new TH1D(RecoilPuppiTopPerp+NJetBins[iJets],"",200,-400,400);

    recoilResponse[iJets] = new TH1D("recoilResponse"+NJetBins[iJets],"",200,-10,10);
    recoilPuppiResponse[iJets] = new TH1D("recoilPuppiResponse"+NJetBins[iJets],"",200,-10,10);

    for (int iPtBins=0; iPtBins<nZPtBins; ++iPtBins){

      recoilZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);
      recoilZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);

      recoilPuppiZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);
      recoilPuppiZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-400,400);

      recoilResponse_Ptbins_nJetsH[iJets][iPtBins] = new TH1D("recoilResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-10,10);
      recoilPuppiResponse_Ptbins_nJetsH[iJets][iPtBins] = new TH1D("recoilPuppiResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",200,-10,10);

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

  int nPtBins = 9;
  float ptBins[10] = {10,15,20,25,30,40,50,100,200,1000};

  int nEtaBins = 5;
  float etaBins[6] = {0, 1.0, 1.479, 1.653, 2.1, 2.5}; 

  TString PtBins[9] = {"Pt10to15",
		       "Pt15to20",
		       "Pt20to25",
		       "Pt25to30",		       
		       "Pt30to40",
		       "Pt40to50",
		       "Pt50to100",
		       "Pt100to200",
		       "PtGt200"};

  TString EtaBins[5] = {"EtaLt1p0",
			"Eta1p0to1p48",
			"Eta1p48to1p65",
			"Eta1p65to2p1",
			"EtaGt2p1"};

  TString JetBins[3] = {"Jet0","Jet1","JetGe2"};

  int nIso1Bins = 4;
  float iso1Bins[5] = {0.,0.1,0.2,0.5,1.0}; 
  TString Iso1Bins[4] = {"Iso1Lt0p1",
			 "Iso10p1to0p2",
			 "Iso10p2to0p5",
			 "Iso1Gt0p5"};

  int nIso2Bins = 4;
  float iso2Bins[5] = {0.,0.15,0.3,0.5,1.0}; 
  TString Iso2Bins[4] = {"Iso2Lt0p15",
			 "Iso20p15to0p3",
			 "Iso20p3to0p5",
			 "Iso2Gt0p5"};

  //*****  create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //
  TH1F * etaBinsH = new TH1F("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  etaBinsH->Draw();
  file->cd();
  etaBinsH->Write("etaBinsH");

  //*****  create pt histograms with pt ranges associated to their names ***** //
  TH1D * ptBinsH =  new TH1D("ptBinsH", "ptBinsH", nPtBins, ptBins);
  ptBinsH->Draw();
  ptBinsH->GetXaxis()->Set(nPtBins, ptBins);
  for (int i=0; i<nPtBins; i++){ ptBinsH->GetXaxis()->SetBinLabel(i+1, PtBins[i]);}
  ptBinsH->Draw();
  file->cd();
  ptBinsH->Write("ptBinsH");

  // create histograms with isolation binning 1 (single e channel)
  TH1D * iso1BinsH =  new TH1D("iso1BinsH", "iso1BinsH", nIso1Bins, iso1Bins);
  iso1BinsH->Draw();
  iso1BinsH->GetXaxis()->Set(nIso1Bins, iso1Bins);
  for (int i=0; i<nIso1Bins; i++){ iso1BinsH->GetXaxis()->SetBinLabel(i+1, Iso1Bins[i]);}
  iso1BinsH->Draw();
  file->cd();
  iso1BinsH->Write("iso1BinsH");

  // create histograms with isolation binning 2 (e+mu channel)
  TH1D * iso2BinsH =  new TH1D("iso2BinsH", "iso2BinsH", nIso2Bins, iso2Bins);
  iso2BinsH->Draw();
  iso2BinsH->GetXaxis()->Set(nIso2Bins, iso2Bins);
  for (int i=0; i<nIso2Bins; i++){ iso2BinsH->GetXaxis()->SetBinLabel(i+1, Iso2Bins[i]);}
  iso2BinsH->Draw();
  file->cd();
  iso2BinsH->Write("iso2BinsH");

  TH1D * ZMassPass = new TH1D("ZMassPass","",80,50,130);
  TH1D * ZMassFail = new TH1D("ZMassFail","",80,50,130);

  TH1F * ZMassJetEtaPtPass[3][5][9];
  TH1F * ZMassJetEtaPtFail[3][5][9];

  TH1F * ZMassEtaPtPass[5][9];
  TH1F * ZMassEtaPtFail[5][9];

  TH1F * ZMassIso1EtaPtPass[4][5][9];
  TH1F * ZMassIso1EtaPtFail[4][5][9];

  TH1F * ZMassIso2EtaPtPass[4][5][9];
  TH1F * ZMassIso2EtaPtFail[4][5][9];

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassEtaPtPass[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
      ZMassEtaPtFail[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      for (int iJet=0; iJet<3; ++iJet) {
	ZMassJetEtaPtPass[iJet][iEta][iPt] = new TH1F("ZMass"+JetBins[iJet]+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
	ZMassJetEtaPtFail[iJet][iEta][iPt] = new TH1F("ZMass"+JetBins[iJet]+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      }
      for (int iIso=0; iIso<nIso1Bins; iIso++) {
	ZMassIso1EtaPtPass[iIso][iEta][iPt] = new TH1F("ZMass"+Iso1Bins[iIso]+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
        ZMassIso1EtaPtFail[iIso][iEta][iPt] = new TH1F("ZMass"+Iso1Bins[iIso]+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      }
      for (int iIso=0; iIso<nIso2Bins; iIso++) {
	ZMassIso2EtaPtPass[iIso][iEta][iPt] = new TH1F("ZMass"+Iso2Bins[iIso]+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
        ZMassIso2EtaPtFail[iIso][iEta][iPt] = new TH1F("ZMass"+Iso2Bins[iIso]+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      }
    }
  }

  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);

  // PILE UP REWEIGHTING - OPTIONS

  // reweighting official recipe 
  // initialize pile up object
  PileUp * PUofficial = new PileUp();
  
  if (applyPUreweighting_official) {
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read"); 
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read"); 
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(PileUpMCHist);
    PUofficial->set_h_data(PU_data); 
    PUofficial->set_h_MC(PU_mc);
  }

  // Met recoil corrections
  RecoilCorrector recoilPFMetCorrector(RecoilFileName);

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
  TFile * workspaceFile = new TFile(TString(cmsswBase)+"/src/"+CorrectionWorkspaceFile);
  RooWorkspace *correctionWS = (RooWorkspace*)workspaceFile->Get("w");

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
	nvertH->Fill(analysisTree.primvertex_count,weight);
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
      
      unsigned int nSingleEleFilter = 0;
      bool isSingleEleFilter = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==ElectronHLTFilterName) {
	  nElectronFilter = i;
	  isElectronFilter = true;
	}
      }

      if (applyTrigger) {
	if (!isElectronFilter) {
	  cout << "Filter " << ElectronHLTFilterName << " not found" << endl;
	  exit(-1);
	}
      }

      //      std::cout << "HLT Filters " << std::endl; 

      float pfmet_ex = analysisTree.pfmetcorr_ex;
      float pfmet_ey = analysisTree.pfmetcorr_ey;
      float pfmet_phi = analysisTree.pfmetcorr_phi;
      float pfmet = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);

      float puppimet_ex = analysisTree.puppimet_ex;
      float puppimet_ey = analysisTree.puppimet_ey;
      float puppimet_phi = analysisTree.puppimet_phi;
      float puppimet = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);

      // electron selection

      // probe electrons
      vector<unsigned int> allElectrons; allElectrons.clear();
      vector<float> allElectronsIso; allElectronsIso.clear();
      vector<bool> allElectronsPassedId; allElectronsPassedId.clear();
      vector<bool> allElectronsPassedIdIso; allElectronsPassedIdIso.clear();

      // tag electrons
      vector<unsigned int> isoIdElectrons; isoIdElectrons.clear();
      vector<bool> isoIdElectronsIsTriggerMatched; isoIdElectronsIsTriggerMatched.clear();
      vector<float> isoIdElectronsIso; isoIdElectronsIso.clear();

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
	if (analysisTree.electron_pt[im]<ptElectronProbeCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaElectronProbeCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronProbeCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzElectronProbeCut) continue;

	bool electronTriggerMatch = false;

	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.electron_eta[im],analysisTree.electron_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nElectronFilter])
	    electronTriggerMatch = true;
	}

	if (!applyTrigger) 
	  electronTriggerMatch = true;
	
	// isolation
	float absIso = 0; 
	if (applyRhoCorrectedIso) {
	  absIso = abs_Iso_et(im, &analysisTree, isoCone);
	}
	else {
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
	}
	float relIso = absIso/analysisTree.electron_pt[im];

	bool isPassed = true;
	if (fabs(analysisTree.electron_dxy[im])>dxyElectronCut) isPassed = false;
	if (fabs(analysisTree.electron_dz[im])>dzElectronCut)  isPassed = false;
	bool electronId = true;
	if (electronIdType==1)
	  electronId = analysisTree.electron_mva_wp80_Iso_Fall17_v1[im]>0.5;
	else if (electronIdType==2)
	  electronId = analysisTree.electron_mva_wp90_Iso_Fall17_v1[im]>0.5;
	else if (electronIdType==3)
	  electronId = analysisTree.electron_mva_wp80_general_Spring16_v1[im]>0.5;
	else if (electronIdType==4)
	  electronId = analysisTree.electron_mva_wp90_general_Spring16_v1[im]>0.5;
	if (!analysisTree.electron_pass_conversion[im]) electronId = false;
	if (analysisTree.electron_nmissinginnerhits[im]>1) electronId = false;
	bool isPassedId = isPassed && electronId;
	bool isPassedIdIso = isPassed && electronId && relIso<isoElectronCut;

	allElectrons.push_back(im);
	allElectronsPassedId.push_back(isPassedId);
	allElectronsPassedIdIso.push_back(isPassedIdIso);
	allElectronsIso.push_back(relIso);

	if (isPassedIdIso) { 
	  isoIdElectrons.push_back(im);
	  isoIdElectronsIso.push_back(relIso);
	  isoIdElectronsIsTriggerMatched.push_back(electronTriggerMatch);
	}
      }
      // end of electron selection

      //      std::cout << "Electron selection -> " << std::endl;

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
	  if (isTriggerMatched && analysisTree.electron_pt[index1]>ptElectronCut && fabs(analysisTree.electron_eta[index1])<etaElectronCut) {
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
	      float isoProbe = TMath::Min(allElectronsIso[iE],float(iso1Bins[nIso1Bins]-0.001));
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      int iso1Bin = binNumber(isoProbe,nIso1Bins,iso1Bins);
	      int iso2Bin = binNumber(isoProbe,nIso2Bins,iso2Bins);

	      bool passId = allElectronsPassedId[iE];

	      if (ptBin<0) continue;
	      if (etaBin<0) continue;


	      //	      std::cout << "ptBin = " << ptBin
	      //			<< "   etaBin = " << etaBin
	      //			<< "   iso1Bin = " << iso1Bin
	      //			<< "   iso2Bin = " << iso2Bin << std::endl;

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
	  
		// tight pfJetId
		bool isPFJetId = tightJetiD_2017(analysisTree,int(jet));
		if (!isPFJetId) continue;

		if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
		  nJets30++;
		  if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		    nJets30etaCut++;
		  }
		}
	      }	

	      int jetBin = nJets30etaCut;
	      if (jetBin>2) jetBin = 2;

	      float mass = (electron1+electron2).M();
	      //	      cout << "probe electron : eta = " << analysisTree.electron_eta[indexProbe]
	      //		   << "   pt = " << analysisTree.electron_eta[indexProbe] << endl;
	      if (allElectronsPassedIdIso[iE]) {
		ZMassPass->Fill(mass,weight);
		ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtPass[jetBin][etaBin][ptBin]->Fill(mass,weight);
	      }
	      else {
		ZMassFail->Fill(mass,weight);
		ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtFail[jetBin][etaBin][ptBin]->Fill(mass,weight);
	      }
	      for (int iIso=0; iIso<nIso1Bins; ++iIso) {
		if (iIso==iso1Bin && passId) 
		  ZMassIso1EtaPtPass[iso1Bin][etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassIso1EtaPtFail[iso1Bin][etaBin][ptBin]->Fill(mass,weight);
	      }
	      for (int iIso=0; iIso<nIso2Bins; ++iIso) {
		if (iIso==iso2Bin && passId) 
		  ZMassIso2EtaPtPass[iso2Bin][etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassIso2EtaPtFail[iso2Bin][etaBin][ptBin]->Fill(mass,weight);
	      }
	    }
	  }

	  for (unsigned int im2=im1+1; im2<isoIdElectrons.size(); ++im2) {
	    unsigned int index2 = isoIdElectrons[im2];
	    float q1 = analysisTree.electron_charge[index1];
	    float q2 = analysisTree.electron_charge[index2];
	    bool isTriggerMatched1 = isoIdElectronsIsTriggerMatched[im1] && 
	      analysisTree.electron_pt[index1]>ptElectronCut&&
	      fabs(analysisTree.electron_eta[index1])<etaElectronCut;
	    bool isTriggerMatched2 = isoIdElectronsIsTriggerMatched[im2]&& 
	      analysisTree.electron_pt[index2]>ptElectronCut &&
	      fabs(analysisTree.electron_eta[index2])<etaElectronCut; 
	    bool isTriggerMatched = isTriggerMatched1 || isTriggerMatched2;
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
	    
	    correctionWS->var("e_pt")->setVal(ptEle1);
	    correctionWS->var("e_eta")->setVal(etaEle1);
	    double trkSF1 = correctionWS->function("e_reco_ratio")->getVal();

	    correctionWS->var("e_pt")->setVal(ptEle2);
	    correctionWS->var("e_eta")->setVal(etaEle2);
	    double trkSF2 = correctionWS->function("e_reco_ratio")->getVal();
	    

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
	    double effTrigData = 1 - (1-effDataTrig1)*(1-effDataTrig2);

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
	   
	  //	  std::cout << "Ok1" << std::endl;

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
	    bool isPFJetId = tightJetiD_2017(analysisTree,int(jet));
	    if (!isPFJetId) continue;

	    bool noisyJet = analysisTree.pfjet_pt[jet]<50 && absJetEta > 2.65 && absJetEta < 3.139;

	    if (noisyJet) continue;


	    
	    if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
	      nJets30++;
	      HT30 += analysisTree.pfjet_pt[jet];
	      if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		HT30etaCut += analysisTree.pfjet_pt[jet];
		nJets30etaCut++;
	      }
	    }
	    if (analysisTree.pfjet_pt[jet]>jetPtLowCut) {
	      nJets20++;
	      HT20 += analysisTree.pfjet_pt[jet]; 
	      if (fabs(analysisTree.pfjet_eta[jet])<jetEtaTrkCut) {
		HT20etaCut += analysisTree.pfjet_pt[jet];
		nJets20etaCut++;
	      }
	    }	  
	  }

	  float visiblePx = dielectron.Px();
	  float visiblePy = dielectron.Py();
	  if (applyRecoilOnGenerator) {
	    visiblePx = genL.Px();
	    visiblePy = genL.Py();
	  }

	  //	  std::cout << "Ok2" << std::endl;

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

	    //	    std::cout << "Ok3" << std::endl;

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

	    float dielectronPt = dielectron.Pt();
	    float dielectronEta = dielectron.Eta();
	    float dielectronPhi = dielectron.Phi();

	    float unitX = dielectron.Px()/dielectron.Pt();
	    float unitY = dielectron.Py()/dielectron.Pt();
	    float phiUnit = TMath::ATan2(unitY,unitX);
	    float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
	    float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
	    
	    // jet pt bin
	    int jetBin = 0;
	    if (nJets30==1)
	      jetBin = 1;
	    else if (nJets30>1)
	      jetBin = 2;

	    // Z pt bin
	    int ptBin = binNumber(TMath::Min(float(dielectronPt),float(999)),nZPtBins,zPtBins);

	    metSelNJets[jetBin]->Fill(pfmet,weight);
	    puppimetSelNJets[jetBin]->Fill(puppimet,weight);

	    float recoilParal = 0;
	    float recoilPerp  = 0;
	    float responseHad = 0;
	    computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilParal,recoilPerp,responseHad);

	    float recoilPuppiParal = 0;
	    float recoilPuppiPerp  = 0;
	    float responsePuppiHad = 0;
	    computeRecoil(puppimet_ex,puppimet_ey,unitX,unitY,perpUnitX,perpUnitY,dielectronPt,recoilPuppiParal,recoilPuppiPerp,responsePuppiHad);

	    for (int iScale=0; iScale<21; ++ iScale) {
	      float scaleFactor = 0.9 + 0.01*float(iScale);
	      metSelScaleH[iScale]->Fill(pfmet*scaleFactor,weight);
	      puppimetSelScaleH[iScale]->Fill(puppimet*scaleFactor,weight);
	    }

	    if (mass>70&&mass<110) {
	      metZSelH->Fill(pfmet,weight);
	      puppimetZSelH->Fill(puppimet,weight);
	      
	      metZSelNJets[jetBin]->Fill(pfmet,weight);
	      puppimetZSelNJets[jetBin]->Fill(puppimet,weight);

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
	      
	    }

	    if (mass<70||mass>110) { // Z-depleted region
	      
	      metTopSelNJets[jetBin]->Fill(pfmet,weight);
	      puppimetTopSelNJets[jetBin]->Fill(puppimet,weight);
	      
	      // pfmet
	      recoilTopParalH[jetBin]->Fill(recoilParal,weight);
	      recoilTopPerpH[jetBin]->Fill(recoilPerp,weight);
	      
	      // puppimet
	      recoilPuppiTopParalH[jetBin]->Fill(recoilPuppiParal,weight);
	      recoilPuppiTopPerpH[jetBin]->Fill(recoilPuppiPerp,weight);
	      
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
