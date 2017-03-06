#include "DesyTauAnalyses/NTupleMaker/interface/Unfold.C"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>

// categories
// em_inclusive
// em_0jet
// em_boosted
// em_vbf
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"

int main(int argc, char * argv[]) {
  
  using namespace std;
    
  // **** configuration
  Config cfg(argv[1]);
  TString directory(argv[2]);
  TString category(argv[3]);
  TString Suffix(argv[4]);

  bool applyBTagVeto = cfg.get<bool>("ApplyBTagVeto");
  bool applyBTagWeight = cfg.get<bool>("ApplyBTagWeight");
  bool useRun2TopPt = cfg.get<bool>("UseRun2TopPt");
  bool applyMTCut = cfg.get<bool>("ApplyMTCut");

  double lumi = 35900;
  double DYNorm = 1.02;

  TString DataFile("MuonEG_Run2016");
  TString TauCheck("TauCheck");
  
  TString btagVeto("");
  TString btagVetoUp("");
  TString btagVetoDown("");
  TString mistagVetoUp("");
  TString mistagVetoDown("");
  TString btagVetoForData("");
  
  if (applyBTagVeto) {
    btagVeto = "&&nbtag==0";
    btagVetoUp = "&&nbtag==0";
    btagVetoDown = "&&nbtag==0";
    mistagVetoUp = "&&nbtag==0";
    mistagVetoDown = "&&nbtag==0";
    if (applyBTagWeight) {
      btagVetoForData = "&&nbtag==0";
      btagVeto = "";
      btagVetoUp = "";
      btagVetoDown = "";
      mistagVetoUp = "";
      mistagVetoDown = "";
    }
  }

  TString mTCut("");
  TString mTCut_eUp("");
  TString mTCut_eDown("");
  if (applyMTCut) {
    mTCut = "&&mTdileptonMET<60";
    mTCut_eUp   = "&&mTdileptonMET_eUp<60";
    mTCut_eDown = "&&mTdileptonMET_eDown<60";

  }  

  TString CutsKine         = "&&pt_1>13&&pt_2>15&&TMath::Max(pt_1,pt_2)>24";
  TString CutsKine_eUp     = "&&pt_Up_1>13&&pt_2>15&&TMath::Max(pt_Up_1,pt_2)>24";
  TString CutsKine_eDown   = "&&pt_Down_1>13&&pt_2>15&&TMath::Max(pt_Down_1,pt_2)>24";

  CutsKine += mTCut;
  CutsKine_eUp += mTCut_eUp;
  CutsKine_eDown += mTCut_eDown;

  TString CutsIso          = "&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5";
  TString CutsIsoSS        = "&&iso_1<0.50&&iso_2>0.2&&iso_2<0.5&&extraelec_veto<0.5&&extramuon_veto<0.5";
  TString CutsCategory             = "&&dzeta>-35";
  TString CutsCategoryJesUp        = "&&dzeta>-35";
  TString CutsCategoryJesDown      = "&&dzeta>-35";
  TString CutsCategoryEUp          = "&&dzeta_eUp>-35";
  TString CutsCategoryEDown        = "&&dzeta_eDown>-35";

  // definition of 2D variables
  // and 2D histograms

  int nBinsX = 100;
  int nBinsY = 100;
  float binsX[50];
  float binsY[50];

  TString VariableX       = "m_vis";
  TString VariableX_eUp   = "m_vis_eUp";
  TString VariableX_eDown = "m_vis_eDown";

  TString VariableY         = "pt_2";
  TString VariableY_eUp     = "pt_2";
  TString VariableY_eDown   = "pt_2";
  TString VariableY_jesUp   = "pt_2";
  TString VariableY_jesDown = "pt_2";

  // 0jet
  int nBinsX_0jet = 12;
  float binsX_0jet[13] = {0,50,55,60,65,70,75,80,85,90,95,100,400};
  int nBinsY_0jet = 6;
  float binsY_0jet[7] = {15,20,25,30,35,40,300};

  // boosted 
  int nBinsX_boosted = 10;
  float binsX_boosted[11] = {0,80,90,100,110,120,130,140,150,160,300};
  int nBinsY_boosted = 6;
  float binsY_boosted[7] = {0,100,150,200,250,300,5000};

  // vbf
  int nBinsX_vbf = 5;
  float binsX_vbf[6] = {0,95,115,135,155,400}; 
  int nBinsY_vbf = 4;
  float binsY_vbf[5] = {300,700,1100,1500,10000};

  nBinsX = nBinsX_0jet;
  nBinsY = nBinsY_0jet;
  for (int iB=0; iB<=nBinsX; ++iB) binsX[iB] = binsX_0jet[iB];
  for (int iB=0; iB<=nBinsY; ++iB) binsY[iB] = binsY_0jet[iB];

  // Weights
  TString Weight = "mcweight*puweight*effweight*0.978824*";
  if (applyBTagWeight)
    Weight += "btag0weight*";

  TString qcdweight("2.30*");

  TString topweight("topptweight*");
  if (useRun2TopPt)
    topweight = "topptweightRun2*";

  TString zptmassweight("zptmassweight*");

  TString ggScaleWeightUp = "(0.9421 - 0.00001699*pt_2)*";
  TString ggScaleWeightDown = "(1.0579 + 0.00001699*pt_2)*";

  // Drell-Yan corrections (in slices of mjj/pt_sv)
  float dyCorr2D[10];

  for (int iDY=0; iDY<10; ++iDY)
    dyCorr2D[iDY] = 1;

  TString JesUncName[29] = {"AbsoluteFlavMap",
			    "AbsoluteMPFBias",
			    "AbsoluteScale",
			    "AbsoluteStat",
			    "FlavorQCD",
			    "Fragmentation",
			    "PileUpDataMC",
			    "PileUpEnvelope",
			    "PileUpPtBB",
			    "PileUpPtEC1",
			    "PileUpPtEC2",
			    "PileUpPtHF",
			    "PileUpPtRef",
			    "RelativeFSR",
			    "RelativeJEREC1",
			    "RelativeJEREC2",
			    "RelativeJERHF",
			    "RelativePtBB",
			    "RelativePtEC1",
			    "RelativePtEC2",
			    "RelativePtHF",
			    "RelativeStatEC",
			    "RelativeStatFSR",
			    "RelativeStatHF",
			    "RelativeBal",
			    "SinglePionECAL",
			    "SinglePionHCAL",
			    "TimePtEta",
			    "Total"};


  TString CutsCategoryJesSys[29][2];
  TString VarYJesSys[29][2];

  if (category=="em_0jet") {

    CutsCategory   = "&&njets==0&&dzeta>-35";
    CutsIsoSS = "&&iso_1<0.3&&iso_2>0.1&&iso_2<0.3&&extraelec_veto<0.5&&extramuon_veto<0.5";
    CutsCategoryJesUp   = "&&njets_Up==0&&dzeta>-35";
    CutsCategoryJesDown = "&&njets_Down==0&&dzeta>-35";
    CutsCategoryEUp     = "&&njets==0&&dzeta_eUp>-35";
    CutsCategoryEDown   = "&&njets==0&&dzeta_eDown>-35";

    for (int iJes=0; iJes<29; ++iJes) {      
      CutsCategoryJesSys[iJes][0] = "&&njets_"+JesUncName[iJes]+"Up==0&&dzeta>-35";
      CutsCategoryJesSys[iJes][1] = "&&njets_"+JesUncName[iJes]+"Down==0&&dzeta>-35";
      VarYJesSys[iJes][0] = "pt_2";
      VarYJesSys[iJes][1] = "pt_2";
    }

    VariableX       = "m_vis";
    VariableX_eUp   = "m_vis_eUp";
    VariableX_eDown = "m_vis_eDown";

    VariableY         = "pt_2";
    VariableY_eUp     = "pt_2";
    VariableY_eDown   = "pt_2";
    VariableY_jesUp   = "pt_2";
    VariableY_jesDown = "pt_2";

    nBinsX = nBinsX_0jet;
    nBinsY = nBinsY_0jet;
    for (int iB=0; iB<=nBinsX; ++iB) binsX[iB] = binsX_0jet[iB];
    for (int iB=0; iB<=nBinsY; ++iB) binsY[iB] = binsY_0jet[iB];

    qcdweight = "2.26*";

    ggScaleWeightUp = "(0.9421 - 0.00001699*pt_2)*";
    ggScaleWeightDown = "(1.0579 + 0.00001699*pt_2)*";

  }
  else if (category=="em_boosted") {

    CutsCategory   = "&&(njets==1 || (njets==2 && mjj<300) || njets>2)&&dzeta>-35";
    CutsIsoSS = "&&iso_1<0.3&&iso_2>0.1&&iso_2<0.3&&extraelec_veto<0.5&&extramuon_veto<0.5";
    CutsCategoryJesUp   = "&&(njets_Up==1 || (njets_Up==2 && mjj_Up<300) || njets_Up>2)&&dzeta>-35";
    CutsCategoryJesDown = "&&(njets_Down==1 || (njets_Down==2 && mjj_Down<300) || njets_Down>2)&&dzeta>-35";
    CutsCategoryEUp     = "&&(njets==1 || (njets==2 && mjj<300) || njets>2)&&dzeta_eUp>-35";
    CutsCategoryEDown   = "&&(njets==1 || (njets==2 && mjj<300) || njets>2)&&dzeta_eDown>-35";

    for (int iJes=0; iJes<29; ++iJes) {      
      CutsCategoryJesSys[iJes][0] =  "&&(njets_"+JesUncName[iJes]+"Up==1 || (njets_"+JesUncName[iJes]+"Up==2 && mjj_"+JesUncName[iJes]+"Up<300) || njets_"+JesUncName[iJes]+"Up>2)&&dzeta>-35";
      CutsCategoryJesSys[iJes][1] =  "&&(njets_"+JesUncName[iJes]+"Down==1 || (njets_"+JesUncName[iJes]+"Down==2 && mjj_"+JesUncName[iJes]+"Down<300) || njets_"+JesUncName[iJes]+"Down>2)&&dzeta>-35";
      VarYJesSys[iJes][0] = "pt_sv";
      VarYJesSys[iJes][1] = "pt_sv";
    }

    VariableX       = "m_sv";
    VariableX_eUp   = "m_sv_eUp";
    VariableX_eDown = "m_sv_eDown";

    VariableY         = "pt_sv";
    VariableY_eUp     = "pt_sv";
    VariableY_eDown   = "pt_sv";
    VariableY_jesUp   = "pt_sv";
    VariableY_jesDown = "pt_sv";

    nBinsX = nBinsX_boosted;
    nBinsY = nBinsY_boosted;
    for (int iB=0; iB<=nBinsX; ++iB) binsX[iB] = binsX_boosted[iB];
    for (int iB=0; iB<=nBinsY; ++iB) binsY[iB] = binsY_boosted[iB];
    
    ggScaleWeightUp = "(0.9358 + 0.00088712 * pt_sv)*";
    ggScaleWeightDown = "(1.0642 - 0.00088712 * pt_sv)*";
    
    qcdweight = "2.25*";

  }
  else if (category=="em_vbf") {

    CutsCategory = "&& njets==2 && mjj>=300 && dzeta>-10";
    CutsIsoSS = "&&iso_1<0.5&&iso_2>0.2&&iso_2<0.5&&extraelec_veto<0.5&&extramuon_veto<0.5";
    CutsCategoryJesUp   = "&& njets_Up==2 && mjj_Up>=300 && dzeta>-10";
    CutsCategoryJesDown = "&& njets_Down==2 && mjj_Down>=300 && dzeta>-10";
    CutsCategoryEUp   = "&& njets==2 && mjj>=300 && dzeta_eUp>-10";
    CutsCategoryEDown = "&& njets==2 && mjj>=300 && dzeta_eDown>-10";

    for (int iJes=0; iJes<29; ++iJes) {      
      CutsCategoryJesSys[iJes][0] = "&& njets_"+JesUncName[iJes]+"Up==2 && mjj_"+JesUncName[iJes]+"Up>=300 && dzeta>-10"; 
      CutsCategoryJesSys[iJes][1] = "&& njets_"+JesUncName[iJes]+"Down==2 && mjj_"+JesUncName[iJes]+"Down>=300 && dzeta>-10"; 
      VarYJesSys[iJes][0] = "mjj_"+JesUncName[iJes]+"Up";
      VarYJesSys[iJes][1] = "mjj_"+JesUncName[iJes]+"Down";
    }

    VariableX       = "m_sv";
    VariableX_eUp   = "m_sv_eUp";
    VariableX_eDown = "m_sv_eDown";

    VariableY         = "mjj";
    VariableY_eUp     = "mjj";
    VariableY_eDown   = "mjj";
    VariableY_jesUp   = "mjj_Up";
    VariableY_jesDown = "mjj_Down";

    nBinsX = nBinsX_vbf;
    nBinsY = nBinsY_vbf;
    for (int iB=0; iB<=nBinsX; ++iB) binsX[iB] = binsX_vbf[iB];
    for (int iB=0; iB<=nBinsY; ++iB) binsY[iB] = binsY_vbf[iB];

    qcdweight = "2.84*";

    ggScaleWeightUp = "(1.032 + 0.00010196 * mjj)*";
    ggScaleWeightDown = "(0.968 - 0.00010196 * mjj)*";
    
    dyCorr2D[0] = 1.06;
    dyCorr2D[1] = 0.98;
    dyCorr2D[2] = 0.95;
    dyCorr2D[3] = 0.95;

  }

  // applying btag veto
  CutsCategory += btagVeto;
  CutsCategoryEUp += btagVeto;
  CutsCategoryEDown += btagVeto;
  CutsCategoryJesUp += btagVeto;
  CutsCategoryJesDown += btagVeto;
  for (int iJes=0; iJes<29; ++iJes ) {
    for (int iUp=0; iUp<2; ++iUp) {
      CutsCategoryJesSys[iJes][iUp] += btagVeto;
    }
  }

  // Binning
  int nBins = nBinsX * nBinsY;
  float xmin = 0.01;
  float xmax = float(nBins) - 0.01;

  // Central cuts and variable
  TString Cuts   = CutsKine + CutsIso   + CutsCategory;
  TString CutsSS = CutsKine + CutsIsoSS + CutsCategory;
  TString Variable = VariableY+":"+VariableX;

  // **********************************************
  // systematics
  // 0 - topUp
  // 1 - topDown
  // 2 - eUp
  // 3 - eDown
  // 4 - jesUp
  // 5 - jesDown
  // 6 - ZPtMassUp
  // 7 - ZPtMassDown
  // 8 - QCD scale gg->H
  // 9 - QCD scale gg->H
  // 10 - 2D Zmumu Up (pt_sv,mjj)  
  // 11 - 2D Zmumu Down (pt_sv,mjj)
  // ***********************************************

  int nSys = 8; // basic shape systematic uncertainties

  TString CutsKineSys[20];
  TString CutsCategorySys[20];
  TString VariableSys[20];
  TString CutsSys[20];
  TString CutsSSSys[20];

  for (int iSys=0; iSys<20; ++iSys) {
    CutsKineSys[iSys] = CutsKine;
    CutsCategorySys[iSys] = CutsCategory;
    VariableSys[iSys] = Variable;
  }

  CutsKineSys[2] = CutsKine_eUp;
  CutsKineSys[3] = CutsKine_eDown; 
  
  CutsCategorySys[2] = CutsCategoryEUp;
  CutsCategorySys[3] = CutsCategoryEDown;
  CutsCategorySys[4] = CutsCategoryJesUp;
  CutsCategorySys[5] = CutsCategoryJesDown;

  VariableSys[2] = VariableY_eUp     + ":" + VariableX_eUp;
  VariableSys[3] = VariableY_eDown   + ":" + VariableX_eDown;
  VariableSys[4] = VariableY_jesUp   + ":" + VariableX;
  VariableSys[5] = VariableY_jesDown + ":" + VariableX;

  for (int iSys=0; iSys<nSys; ++iSys) { 
    CutsSys[iSys] = CutsKineSys[iSys] + CutsIso + CutsCategorySys[iSys];
    CutsSSSys[iSys] = CutsKineSys[iSys] + CutsIsoSS + CutsCategorySys[iSys];
  }

  // JES uncertainties
  TString VariableJesSys[29][2];
  TString CutsJesSys[29][2]; 
  TString CutsSSJesSys[29][2];

  for (int iJes=0; iJes<29; ++iJes ) {
    for (int iUp=0; iUp<2; ++iUp) {
      CutsJesSys[iJes][iUp] = CutsKine + CutsIso +  CutsCategoryJesSys[iJes][iUp];
      CutsSSJesSys[iJes][iUp] = CutsKine + CutsIsoSS +  CutsCategoryJesSys[iJes][iUp];
      VariableJesSys[iJes][iUp] = VarYJesSys[iJes][iUp] + ":" + VariableX;
    }
  }

  // Top weight systematics
  TString topweightSys[20]; 
  for (int iTop=0; iTop<20; ++iTop) {
    if (useRun2TopPt) 
      topweightSys[iTop] = "topptweightRun2*";
    else
      topweightSys[iTop] = "topptweight*";
  }
  topweightSys[0] = "topptweight*topptweight*";
  topweightSys[1] = "";
  if (useRun2TopPt)
    topweightSys[0] = "topptweightRun2*topptweightRun2*";

  // Z Pt,Mass systematics
  TString zptmassweightSys[20];
  for (int iDY=0; iDY<20; ++iDY)
    zptmassweightSys[iDY] = "zptmassweight*";
  zptmassweightSys[6] = "(1.0+1.1*(zptmassweight-1))*";
  zptmassweightSys[7] = "(1.0+0.9*(zptmassweight-1))*";

  TString sysName[20] = {"_CMS_htt_ttbarShape_13TeVUp", // 0
			 "_CMS_htt_ttbarShape_13TeVDown", // 1
			 "_CMS_scale_e_em_13TeVUp", // 2
			 "_CMS_scale_e_em_13TeVDown", // 3
			 "_CMS_scale_j_13TeVUp", // 4
			 "_CMS_scale_j_13TeVDown", // 5
			 "_CMS_htt_dyShape_13TeVUp", // 6
			 "_CMS_htt_dyShape_13TeVDown", // 7
			 "_CMS_scale_gg_13TeVUp", // 8
			 "_CMS_scale_gg_13TeVDown", // 9
			 "_CMS_htt_zmumuShape_0jet_13TeVUp", // 10
			 "_CMS_htt_zmumuShape_0jet_13TeVDown", // 11
			 "",
			 "",
			 "",
			 "",
			 "",
			 "",
			 "",
			 ""
  };

  if (category=="em_boosted") {
    sysName[10] = "_CMS_htt_zmumuShape_boosted_13TeVUp";
    sysName[11] = "_CMS_htt_zmumuShape_boosted_13TeVDown";
  }

  if (category=="em_vbf") {
    sysName[10] = "_CMS_htt_zmumuShape_VBF_13TeVUp";
    sysName[11] = "_CMS_htt_zmumuShape_VBF_13TeVDown";
  }

  TString jesSysName[29][2];
  for (int iJes=0; iJes<29; ++iJes) {
    jesSysName[iJes][0] = "_CMS_scale_j_"+JesUncName[iJes]+"_13TeVUp";
    jesSysName[iJes][1] = "_CMS_scale_j_"+JesUncName[iJes]+"_13TeVDown";
  }

  // defining samples

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TH1D * hist[30];
  TH1D * histSS[30];
  TH1D * histSSrelaxed[30];

  TH1D * histSys[30][20];
  TH1D * histSSSys[30][20];
  TH1D * histSSrelaxedSys[30][20];

  TH1D * histJesSys[30][29][2];
  TH1D * histSSJesSys[30][29][2];

  TH1D * histSSrelaxedJesSys[30][29][2];
  TH1D * histSignal[5][5];
  TH1D * histSignalSys[5][5][20];
  TH1D * histSignalJesSys[5][5][29][2];


  TString sampleNames[30] = {
    DataFile, // data (0)
    "DYJetsToLL_M-50_13TeV-madgraphMLM", // isZTT (1)
    "DYJetsToLL_M-50_13TeV-madgraphMLM", // !isZTT (2)
    "DYJetsToLL_M-10to50_13TeV-madgraphMLM", // isZTT (3)
    "DYJetsToLL_M-10to50_13TeV-madgraphMLM", // !isZTT (4)
    "WJetsToLNu_13TeV-madgraphMLM",    // (5)
    "TTJets_13TeV-powheg",             // (6)
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg",     // (7)
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg", // (8)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg",     // (9)
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg",         // (10)
    "VVTo2L2Nu_13TeV_amcatnloFXFX",    // (11)
    "WWToLNuQQ_13TeV_powheg",          // (12)
    "WZTo2L2Q_13TeV_amcatnloFXFX",     // (13)
    "WZTo1L1Nu2Q_13TeV_amcatnloFXFX",  // (14)
    "WZTo1L3Nu_13TeV_amcatnloFXFX",    // (15)
    "WZJToLLLNu_13TeV_amcatnloFXFX",   // (16)
    "ZZTo4L_13TeV_powheg",             // (17)
    "ZZTo2L2Q_13TeV_amcatnloFXFX",     // (18)
    "WGToLNuG_13TeV-madgraphMLM-pythia8",     // (19)
    "WGstarToLNuMuMu_012Jets_13TeV-madgraph", // (20)
    "WGstarToLNuEE_012Jets_13TeV-madgraph",   // (21)
    "EWKZ2Jets_ZToLL_M-50_13TeV-madgraph",    // (22) 
    "EWKWPlus2Jets_WToLNu_13TeV-madgraph",    // (23)
    "EWKWMinus2Jets_WToLNu_13TeV-madgraph",   // (24)
    "GluGluHToWWTo2L2Nu_M125_13TeV_powheg",      // (25)
    "VBFHToWWTo2L2Nu_M125_13TeV_powheg",         // (26)
    "", // (27)
    "", // (28)
    ""  // (29)
  };

  TString signalSamples[5] = {"GluGluHToTauTau_M",
			      "VBFHToTauTau_M",
			      "ZHToTauTau_M",
			      "WplusHToTauTau_M",
			      "WminusHToTauTau_M"};
  
  TString massH[5] = {"110","120","125","130","140"};

  double xsecSignal[5][5];
  // gg->H
  xsecSignal[0][0] = 57.90*0.0791;
  xsecSignal[0][1] = 52.22*0.0698;
  xsecSignal[0][2] = 48.58*0.0627;
  xsecSignal[0][3] = 45.31*0.0541;
  xsecSignal[0][4] = 36.00*0.0360;

  // VBF
  xsecSignal[1][0] = 4.434*0.0791;
  xsecSignal[1][1] = 3.935*0.0698;
  xsecSignal[1][2] = 3.782*0.0627;
  xsecSignal[1][3] = 3.637*0.0541;
  xsecSignal[1][4] = 3.492*0.0360;

  // ZH
  xsecSignal[2][0] = 1.309*0.0791;
  xsecSignal[2][1] = 0.994*0.0698;
  xsecSignal[2][2] = 0.884*0.0627;
  xsecSignal[2][3] = 0.790*0.0541;
  xsecSignal[2][4] = 0.6514*0.0360;

  // WplusH
  xsecSignal[3][0] = 1.3350*0.0791;
  xsecSignal[3][1] = 0.9558*0.0698;
  xsecSignal[3][2] = 0.8400*0.0627;
  xsecSignal[3][3] = 0.7414*0.0541;
  xsecSignal[3][4] = 0.6308*0.0360;
  
  // WminusH
  xsecSignal[4][0] = 0.8587*0.0791;
  xsecSignal[4][1] = 0.6092*0.0698;
  xsecSignal[4][2] = 0.5328*0.0627;
  xsecSignal[4][3] = 0.4676*0.0541;
  xsecSignal[4][4] = 0.3940*0.0360;

  double xsec[30] = {1, // data (0)
		     DYNorm*5765,  // DY (1)
		     DYNorm*5765,  // DY (2)
		     18610, // DY low mass (3)
		     18610, // DY low mass (4)
		     61526.7, // WJets (5)
		     831.76,  // TT (6)
		     136.95, // ST_t-channel_top (7)
		     80.95,  // ST_t-channel_antitop (8)
		     35.6,   // ST_tW_antitop (9)
		     35.6,   // ST_tW_top_5f (10)
		     11.95,  // VV (11)
		     49.997, // WWToLNuQQ (12)
		     5.595,  // WZTo2L2Q (13)
		     10.71,  // WZTo1L1Nu2Q (14)
		     3.05,   // WZTo1L3Nu (15)
		     4.708,  // WZJets (3L1Nu) (16)
		     1.212,  // ZZTo4L (17)
		     3.22,   // ZZTo2L2Q (18)
		     489.0,  // WGToLNuG (19)
		     2.793,  // WGToLNuMuMu (20)
		     3.526,  // WGToLNuEE   (21)
		     3.987,  // EWK Z    (22)
		     25.62,  // EWK W+   (23)
		     20.20,  // EWK W-   (24)
		     48.58*0.215*0.324*0.324,   // ggHWW  (25)
		     3.782*0.215*0.324*0.324,   // qqHWW  (26)
		     0, // (27)
		     0, // (28)
		     0  // (29)
  };      
  
  TString cuts[30];
  TString cutsSS[30];
  TString cutsSSrelaxed[30];

  TString cutsSys[30][20];
  TString cutsSSSys[30][20];
  TString cutsSSrelaxedSys[30][20];

  TString cutsJesSys[30][29][2];
  TString cutsSSJesSys[30][29][2];
  TString cutsSSrelaxedJesSys[30][29][2];

  TString cutsSignal;
  TString cutsSignalSys[20];
  TString cutsSignalJesSys[20][2];

  // Central cuts ---->
  for (int i=0; i<30; ++i) {
    cuts[i]   = Weight+"(os>0.5"+Cuts+")";
    cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
    cutsSSrelaxed[i] = Weight+qcdweight+"(os<0.5"+CutsSS+")";
  }
  cuts[0] = "(os>0.5"+Cuts+btagVetoForData+")";
  cuts[1] = Weight+zptmassweight+"(os>0.5"+Cuts+"&&isZTT)";
  cuts[2] = Weight+zptmassweight+"(os>0.5"+Cuts+"&&!isZTT)";
  cuts[3] = Weight+"(os>0.5"+Cuts+"&&isZTT)";
  cuts[4] = Weight+"(os>0.5"+Cuts+"&&!isZTT)";
  cuts[6] = Weight+topweight+"(os>0.5"+Cuts+")";

  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+btagVetoForData+")";
  cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";
  cutsSS[3] = Weight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[4] = Weight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";
  cutsSS[6] = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";

  cutsSSrelaxed[0] = "(os<0.5"+CutsSS+btagVetoForData+")";
  cutsSSrelaxed[1] = Weight+zptmassweight+"(os<0.5"+CutsSS+"&&isZTT)";
  cutsSSrelaxed[2] = Weight+zptmassweight+"(os<0.5"+CutsSS+"&&!isZTT)";
  cutsSSrelaxed[3] = Weight+"(os<0.5"+CutsSS+"&&isZTT)";
  cutsSSrelaxed[4] = Weight+"(os<0.5"+CutsSS+"&&!isZTT)";
  cutsSSrelaxed[6] = Weight+topweight+"(os<0.5"+CutsSS+")";

  cutsSignal = Weight+"(os>0.5"+Cuts+")"; 


  // Systematics ->
  for (int iSys=0; iSys<nSys; ++iSys) {

    for (int i=0; i<30; ++i) {
      cutsSys[i][iSys] = Weight+"(os>0.5"+CutsSys[iSys]+")";
      cutsSSSys[i][iSys] = Weight+"(os<0.5"+CutsSys[iSys]+")";
      cutsSSrelaxedSys[i][iSys] = Weight+"(os<0.5"+CutsSSSys[iSys]+")";
    }

    cutsSys[0][iSys] = "(os>0.5"+CutsSys[iSys]+btagVetoForData+")";
    cutsSys[1][iSys] = Weight+zptmassweightSys[iSys]+"(os>0.5"+CutsSys[iSys]+"&&isZTT)";
    cutsSys[2][iSys] = Weight+zptmassweightSys[iSys]+"(os>0.5"+CutsSys[iSys]+"&&!isZTT)";
    cutsSys[3][iSys] = Weight+"(os>0.5"+CutsSys[iSys]+"&&isZTT)";
    cutsSys[4][iSys] = Weight+"(os>0.5"+CutsSys[iSys]+"&&!isZTT)";
    cutsSys[6][iSys] = Weight+topweightSys[iSys]+"(os>0.5"+CutsSys[iSys]+")";

    cutsSSSys[0][iSys] = qcdweight+"(os<0.5"+CutsSys[iSys]+btagVetoForData+")";
    cutsSSSys[1][iSys] = Weight+zptmassweightSys[iSys]+qcdweight+"(os<0.5"+CutsSys[iSys]+"&&isZTT)";
    cutsSSSys[2][iSys] = Weight+zptmassweightSys[iSys]+qcdweight+"(os<0.5"+CutsSys[iSys]+"&&!isZTT)";
    cutsSSSys[3][iSys] = Weight+qcdweight+"(os<0.5"+CutsSys[iSys]+"&&isZTT)";
    cutsSSSys[4][iSys] = Weight+qcdweight+"(os<0.5"+CutsSys[iSys]+"&&!isZTT)";
    cutsSSSys[6][iSys] = Weight+topweightSys[iSys]+qcdweight+"(os<0.5"+CutsSys[iSys]+")";

    cutsSSrelaxedSys[0][iSys] = "(os<0.5"+CutsSSSys[iSys]+btagVetoForData+")";
    cutsSSrelaxedSys[1][iSys] = Weight+zptmassweightSys[iSys]+"(os<0.5"+CutsSSSys[iSys]+"&&isZTT)";
    cutsSSrelaxedSys[2][iSys] = Weight+zptmassweightSys[iSys]+"(os<0.5"+CutsSSSys[iSys]+"&&!isZTT)";
    cutsSSrelaxedSys[3][iSys] = Weight+"(os<0.5"+CutsSSSys[iSys]+"&&isZTT)";
    cutsSSrelaxedSys[4][iSys] = Weight+"(os<0.5"+CutsSSSys[iSys]+"&&!isZTT)";
    cutsSSrelaxedSys[6][iSys] = Weight+topweightSys[iSys]+"(os<0.5"+CutsSSSys[iSys]+")";

    cutsSignalSys[iSys] = Weight+"(os>0.5"+CutsSys[iSys]+")";

  }
  cutsSignalSys[8] = Weight+ggScaleWeightUp+"(os>0.5"+Cuts+")";
  cutsSignalSys[9] = Weight+ggScaleWeightDown+"(os>0.5"+Cuts+")";

  // JES Systematics ->
  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {

      cutsSignalJesSys[iSys][iUp] = Weight+"(os>0.5"+CutsJesSys[iSys][iUp]+")";

      for (int i=0; i<30; ++i) {
	cutsJesSys[i][iSys][iUp] = Weight+"(os>0.5"+CutsJesSys[iSys][iUp]+")";
	cutsSSJesSys[i][iSys][iUp] = Weight+"(os<0.5"+CutsJesSys[iSys][iUp]+")";
	cutsSSrelaxedJesSys[i][iSys][iUp] = Weight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+")";
      }
      
      cutsJesSys[0][iSys][iUp] = "(os>0.5"+CutsJesSys[iSys][iUp]+")";
      cutsJesSys[1][iSys][iUp] = Weight+zptmassweight+"(os>0.5"+CutsJesSys[iSys][iUp]+"&&isZTT)";
      cutsJesSys[2][iSys][iUp] = Weight+zptmassweight+"(os>0.5"+CutsJesSys[iSys][iUp]+"&&!isZTT)";
      cutsJesSys[3][iSys][iUp] = Weight+"(os>0.5"+CutsJesSys[iSys][iUp]+"&&isZTT)";
      cutsJesSys[4][iSys][iUp] = Weight+"(os>0.5"+CutsJesSys[iSys][iUp]+"&&!isZTT)";
      cutsJesSys[6][iSys][iUp] = Weight+topweight+"(os>0.5"+CutsJesSys[iSys][iUp]+")";
      
      cutsSSJesSys[0][iSys][iUp] = qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+")";
      cutsSSJesSys[1][iSys][iUp] = Weight+zptmassweight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+"&&isZTT)";
      cutsSSJesSys[2][iSys][iUp] = Weight+zptmassweight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+"&&!isZTT)";
      cutsSSJesSys[3][iSys][iUp] = Weight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+"&&isZTT)";
      cutsSSJesSys[4][iSys][iUp] = Weight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+"&&!isZTT)";
      cutsSSJesSys[6][iSys][iUp] = Weight+topweight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+")";
      
      cutsSSrelaxedJesSys[0][iSys][iUp] = "(os<0.5"+CutsSSJesSys[iSys][iUp]+")";
      cutsSSrelaxedJesSys[1][iSys][iUp] = Weight+zptmassweight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+"&&isZTT)";
      cutsSSrelaxedJesSys[2][iSys][iUp] = Weight+zptmassweight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+"&&!isZTT)";
      cutsSSrelaxedJesSys[3][iSys][iUp] = Weight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+"&&isZTT)";
      cutsSSrelaxedJesSys[4][iSys][iUp] = Weight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+"&&!isZTT)";
      cutsSSrelaxedJesSys[6][iSys][iUp] = Weight+topweight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+")";
      
    }
  }


  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  // ********************************************
  // ******** Data and MC background samples ****
  // ********************************************
  
  int nSamples = 27;

  for (int i=0; i<nSamples; ++i) { // run over samples
    TFile * file = new TFile(directory+sampleNames[i]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
    TTree * tree = (TTree*)file->Get(TauCheck);
    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();

    TString histName = sampleNames[i] + "_os";
    TString histNameSS = sampleNames[i] + "_ss";
    TString histNameSSrelaxed = sampleNames[i] + "_ss_relaxed";

    TH2D * hist2D = new TH2D(histName,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * hist2DSS = new TH2D(histNameSS,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * hist2DSSrelaxed = new TH2D(histNameSSrelaxed,"",nBinsX,binsX,nBinsY,binsY);

    hist2D->Sumw2();
    hist2DSS->Sumw2();
    hist2DSSrelaxed->Sumw2();

    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
    tree->Draw(Variable+">>"+histNameSSrelaxed,cutsSSrelaxed[i]);

    hist[i] = (TH1D*)Unfold(hist2D);
    histSS[i] = (TH1D*)Unfold(hist2DSS);
    histSSrelaxed[i] = (TH1D*)Unfold(hist2DSSrelaxed);

    std::cout << sampleNames[i] << std::endl;

    if (i==0) {
      for (int iSys=0; iSys<nSys; ++iSys ) {
	histSys[i][iSys] = (TH1D*)hist[i]->Clone(histName+sysName[iSys]);
	histSSSys[i][iSys] = (TH1D*)histSS[i]->Clone(histNameSS+sysName[iSys]);
	histSSrelaxedSys[i][iSys] = (TH1D*)histSSrelaxed[i]->Clone(histNameSSrelaxed+sysName[iSys]);
      }
      for (int iSys=0; iSys<29; ++iSys ) {
	for (int iUp=0; iUp<2; ++iUp) {
	  histJesSys[i][iSys][iUp] = (TH1D*)hist[i]->Clone(histName+jesSysName[iSys][iUp]);
	  histSSJesSys[i][iSys][iUp] = (TH1D*)histSS[i]->Clone(histNameSS+jesSysName[iSys][iUp]);
	  histSSrelaxedJesSys[i][iSys][iUp] = (TH1D*)histSSrelaxed[i]->Clone(histNameSSrelaxed+jesSysName[iSys][iUp]);
	}
      }
    }

    if (i>0) { // MC samples (systematcs and normalization)

      // systematics 
      for (int iSys=0; iSys<nSys; ++iSys ) {


	TH2D * hist2DSys = new TH2D(histName+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
	hist2DSys->Sumw2();
	tree->Draw(VariableSys[iSys]+">>"+histName+sysName[iSys],cutsSys[i][iSys]);
	histSys[i][iSys] = (TH1D*)Unfold(hist2DSys);

	TH2D * hist2DSSSys = new TH2D(histNameSS+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
	hist2DSSSys->Sumw2();
	tree->Draw(VariableSys[iSys]+">>"+histNameSS+sysName[iSys],cutsSSSys[i][iSys]);
	histSSSys[i][iSys] = (TH1D*)Unfold(hist2DSSSys);

	TH2D * hist2DSSrelaxedSys = new TH2D(histNameSSrelaxed+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
	hist2DSSrelaxedSys->Sumw2();
	tree->Draw(VariableSys[iSys]+">>"+histNameSSrelaxed+sysName[iSys],cutsSSrelaxedSys[i][iSys]);
	histSSrelaxedSys[i][iSys] = (TH1D*)Unfold(hist2DSSrelaxedSys);

	std::cout << sampleNames[i] << sysName[iSys] << std::endl;

      }

      // JES systematics
      for (int iSys=0; iSys<29; ++iSys ) {
	for (int iUp = 0; iUp<2; ++iUp) {

	  TH2D * hist2DJesSys = new TH2D(histName+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	  hist2DJesSys->Sumw2();
	  tree->Draw(VariableJesSys[iSys][iUp]+">>"+histName+jesSysName[iSys][iUp],cutsJesSys[i][iSys][iUp]);
	  histJesSys[i][iSys][iUp] = (TH1D*)Unfold(hist2DJesSys);

	  TH2D * hist2DSSJesSys = new TH2D(histNameSS+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	  hist2DSSJesSys->Sumw2();
	  tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameSS+jesSysName[iSys][iUp],cutsSSJesSys[i][iSys][iUp]);
	  histSSJesSys[i][iSys][iUp] = (TH1D*)Unfold(hist2DSSJesSys);
	  
	  TH2D * hist2DSSrelaxedJesSys = new TH2D(histNameSSrelaxed+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	  hist2DSSrelaxedJesSys->Sumw2();
	  tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameSSrelaxed+jesSysName[iSys][iUp],cutsSSrelaxedJesSys[i][iSys][iUp]);
	  histSSrelaxedJesSys[i][iSys][iUp] = (TH1D*)Unfold(hist2DSSrelaxedJesSys);
	  
	  std::cout << sampleNames[i] << jesSysName[iSys][iUp] << std::endl;

	}
      }


      for (int iB=1; iB<=nBins; ++iB) {
	double x = hist[i]->GetBinContent(iB);
	double e = hist[i]->GetBinError(iB);
    	hist[i]->SetBinContent(iB,norm*x);
    	hist[i]->SetBinError(iB,norm*e);
	x = histSS[i]->GetBinContent(iB);
	e = histSS[i]->GetBinError(iB);
    	histSS[i]->SetBinContent(iB,norm*x);
    	histSS[i]->SetBinError(iB,norm*e);
	x = histSSrelaxed[i]->GetBinContent(iB);
	e = histSSrelaxed[i]->GetBinError(iB);
    	histSSrelaxed[i]->SetBinContent(iB,norm*x);
    	histSSrelaxed[i]->SetBinError(iB,norm*e);
	for (int iSys=0; iSys<nSys; ++iSys) {
	  x = histSys[i][iSys]->GetBinContent(iB);
	  e = histSys[i][iSys]->GetBinError(iB);
	  histSys[i][iSys]->SetBinContent(iB,norm*x);
	  histSys[i][iSys]->SetBinError(iB,norm*e);
	  x = histSSSys[i][iSys]->GetBinContent(iB);
	  e = histSSSys[i][iSys]->GetBinError(iB);
	  histSSSys[i][iSys]->SetBinContent(iB,norm*x);
	  histSSSys[i][iSys]->SetBinError(iB,norm*e);
	  x = histSSrelaxedSys[i][iSys]->GetBinContent(iB);
	  e = histSSrelaxedSys[i][iSys]->GetBinError(iB);
	  histSSrelaxedSys[i][iSys]->SetBinContent(iB,norm*x);
	  histSSrelaxedSys[i][iSys]->SetBinError(iB,norm*e);
	}
	for (int iSys=0; iSys<29; ++iSys) {
	  for (int iUp=0; iUp<2; ++iUp) {
	    x = histJesSys[i][iSys][iUp]->GetBinContent(iB);
	    e = histJesSys[i][iSys][iUp]->GetBinError(iB);
	    histJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	    histJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	    x = histSSJesSys[i][iSys][iUp]->GetBinContent(iB);
	    e = histSSJesSys[i][iSys][iUp]->GetBinError(iB);
	    histSSJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	    histSSJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	    x = histSSrelaxedJesSys[i][iSys][iUp]->GetBinContent(iB);
	    e = histSSrelaxedJesSys[i][iSys][iUp]->GetBinError(iB);
	    histSSrelaxedJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	    histSSrelaxedJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	  }
	}
      }
    }

    cout << sampleNames[i] << " -> OS = " << hist[i]->GetSumOfWeights() << " : SS = " << histSS[i]->GetSumOfWeights() << endl;

  }

  // *******************************
  // **** Signal samples ***********
  // *******************************

  for (int i=0; i<5; ++i) { // run over samples
    for (int iM=0; iM<5; ++iM) { // run over masses
      TString sampleName = signalSamples[i] + massH[iM] + "_13TeV_powheg"; 
      TFile * file = new TFile(directory+sampleName+".root");
      TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
      TTree * tree = (TTree*)file->Get(TauCheck);
      double norm = xsecSignal[i][iM]*lumi/histWeightsH->GetSumOfWeights();
      
      TString histName = sampleName + "_os";
      TH2D * hist2D = new TH2D(histName,"",nBinsX,binsX,nBinsY,binsY);
      
      hist2D->Sumw2();
      
      tree->Draw(Variable+">>"+histName,cuts[i]);
      histSignal[i][iM] = (TH1D*)Unfold(hist2D);
      
      std::cout << sampleName << std::endl;

      // systematics 
      int nSysSignal = nSys;
      if (i == 0) nSysSignal = nSys+2; // ggH
      
      for (int iSys=0; iSys<nSysSignal; ++iSys ) {
	TH2D * hist2DSys = new TH2D(histName+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
	hist2DSys->Sumw2();
	tree->Draw(VariableSys[iSys]+">>"+histName+sysName[iSys],cutsSignalSys[iSys]);
	histSignalSys[i][iM][iSys] = (TH1D*)Unfold(hist2DSys);
	std::cout << sampleName << sysName[iSys] << std::endl;
      }

      
      // JES systematics
      for (int iSys=0; iSys<29; ++iSys ) {
	for (int iUp = 0; iUp<2; ++iUp) {
	  TH2D * hist2DJesSys = new TH2D(histName+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	  hist2DJesSys->Sumw2();
	  tree->Draw(VariableJesSys[iSys][iUp]+">>"+histName+jesSysName[iSys][iUp],cutsSignalJesSys[iSys][iUp]);
	  histSignalJesSys[i][iM][iSys][iUp] = (TH1D*)Unfold(hist2DJesSys);
	  std::cout << sampleName << jesSysName[iSys][iUp] << std::endl;
	}
      }
      
      
      for (int iB=1; iB<=nBins; ++iB) {
	double x = histSignal[i][iM]->GetBinContent(iB);
	double e = histSignal[i][iM]->GetBinError(iB);
    	histSignal[i][iM]->SetBinContent(iB,norm*x);
    	histSignal[i][iM]->SetBinError(iB,norm*e);
	for (int iSys=0; iSys<nSysSignal; ++iSys) {
	  x = histSignalSys[i][iM][iSys]->GetBinContent(iB);
	  e = histSignalSys[i][iM][iSys]->GetBinError(iB);
	  histSignalSys[i][iM][iSys]->SetBinContent(iB,norm*x);
	  histSignalSys[i][iM][iSys]->SetBinError(iB,norm*e);
	}
	for (int iSys=0; iSys<29; ++iSys) {
	  for (int iUp=0; iUp<2; ++iUp) {
	    x = histSignalJesSys[i][iM][iSys][iUp]->GetBinContent(iB);
	    e = histSignalJesSys[i][iM][iSys][iUp]->GetBinError(iB);
	    histSignalJesSys[i][iM][iSys][iUp]->SetBinContent(iB,norm*x);
	    histSignalJesSys[i][iM][iSys][iUp]->SetBinError(iB,norm*e);
	  }
	}
      } 

      cout << sampleName << " -> " << histSignal[i][iM]->GetSumOfWeights() << endl;

    }
  }

  // *******************************
  // ***** Drell-Yan samples *******
  // *******************************

  TH1D * histZtt[10];
  TH1D * histZttSys[10][20];
  TH1D * histZttJesSys[10][29][2];

  TH1D * histZttSS[10];
  TH1D * histZttSSSys[10][20];
  TH1D * histZttSSJesSys[10][29][2];

  TH1D * histZttSSrelaxed[10];
  TH1D * histZttSSrelaxedSys[10][20];
  TH1D * histZttSSrelaxedJesSys[10][29][2];

  TH1D * histZll[10];
  TH1D * histZllSys[10][20];
  TH1D * histZllJesSys[10][29][2];

  TH1D * histZllSS[10];
  TH1D * histZllSSSys[10][20];
  TH1D * histZllSSJesSys[10][29][2];

  TH1D * histZllSSrelaxed[10];
  TH1D * histZllSSrelaxedSys[10][20];
  TH1D * histZllSSrelaxedJesSys[10][29][2];

  TString refSamples[5] = {"DYJetsToLL_M-50_13TeV-madgraphMLM",
			   "DY1JetsToLL_M-50_13TeV-madgraphMLM",
			   "DY2JetsToLL_M-50_13TeV-madgraphMLM",
			   "DY3JetsToLL_M-50_13TeV-madgraphMLM",
			   "DY4JetsToLL_M-50_13TeV-madgraphMLM"
  };
  double refXSec[5] = {DYNorm*5765,
		       DYNorm*1.164*1012.5,
		       DYNorm*1.164*332.8,
		       DYNorm*1.164*101.8,
		       DYNorm*1.164*54.8};

  double refEvents[5] = {0,0,0,0,0};

  for (int iDY=0; iDY<5; ++iDY) {
    TFile * file = new TFile(directory+refSamples[iDY]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
    refEvents[iDY] = histWeightsH->GetSumOfWeights();
  }
  TString dySampleNames[9] = {"DYJetsToLL_M-50_13TeV-madgraphMLM",
			      "DYJetsToLL_M-50_13TeV-madgraphMLM",
			      "DYJetsToLL_M-50_13TeV-madgraphMLM",
			      "DYJetsToLL_M-50_13TeV-madgraphMLM",
			      "DYJetsToLL_M-50_13TeV-madgraphMLM",
			      "DY1JetsToLL_M-50_13TeV-madgraphMLM",
			      "DY2JetsToLL_M-50_13TeV-madgraphMLM",
			      "DY3JetsToLL_M-50_13TeV-madgraphMLM",
			      "DY4JetsToLL_M-50_13TeV-madgraphMLM"
  };

  double dyNorm[9];
  dyNorm[0] = lumi*refXSec[0]/refEvents[0];
  dyNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  dyNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);

  TString npartonCuts[9] = {"&&(npartons==0||npartons>4)",
			    "&&npartons==1",
			    "&&npartons==2",
			    "&&npartons==3",
			    "&&npartons==4",
			    "",
			    "",
			    "",
			    ""
  };

  TString cutsZtt[9];
  TString cutsZttSys[9][20];
  TString cutsZttJesSys[9][29][2];

  TString cutsZttSS[9];
  TString cutsZttSSSys[9][20];
  TString cutsZttSSJesSys[9][29][2];

  TString cutsZttSSrelaxed[9];
  TString cutsZttSSrelaxedSys[9][20];
  TString cutsZttSSrelaxedJesSys[9][29][2];
  
  TString cutsZll[9];
  TString cutsZllSys[9][20];
  TString cutsZllJesSys[9][29][2];

  TString cutsZllSS[9];
  TString cutsZllSSSys[9][20];
  TString cutsZllSSJesSys[9][29][2];

  TString cutsZllSSrelaxed[9];
  TString cutsZllSSrelaxedSys[9][20];
  TString cutsZllSSrelaxedJesSys[9][29][2];
  

  for (int iDY=0; iDY<9; ++iDY) {

    cutsZtt[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCuts[iDY]+"&&isZTT>0.5)";
    cutsZttSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+"&&isZTT>0.5)";
    cutsZttSSrelaxed[iDY] = Weight+zptmassweight+"(os<0.5"+CutsSS+npartonCuts[iDY]+"&&isZTT>0.5)";

    cutsZll[iDY]   = Weight+zptmassweight+"(os>0.5"+Cuts+npartonCuts[iDY]+"&&isZTT<0.5)";
    cutsZllSS[iDY] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+"&&isZTT<0.5)";
    cutsZllSSrelaxed[iDY] = Weight+zptmassweight+"(os<0.5"+CutsSS+npartonCuts[iDY]+"&&isZTT<0.5)";

    for (int iSys=0; iSys<nSys; ++iSys) {

      cutsZttSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+"(os>0.5"+CutsSys[iSys]+npartonCuts[iDY]+"&&isZTT>0.5)";
      cutsZllSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+"(os>0.5"+CutsSys[iSys]+npartonCuts[iDY]+"&&isZTT<0.5)";

      cutsZttSSSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+qcdweight+"(os<0.5"+CutsSys[iSys]+npartonCuts[iDY]+"&&isZTT>0.5)";
      cutsZllSSSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+qcdweight+"(os<0.5"+CutsSys[iSys]+npartonCuts[iDY]+"&&isZTT<0.5)";

      cutsZttSSrelaxedSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+"(os<0.5"+CutsSSSys[iSys]+npartonCuts[iDY]+"&&isZTT>0.5)";
      cutsZllSSrelaxedSys[iDY][iSys] = Weight+zptmassweightSys[iSys]+"(os<0.5"+CutsSSSys[iSys]+npartonCuts[iDY]+"&&isZTT<0.5)";

    }

    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {

	cutsZttJesSys[iDY][iSys][iUp] = Weight+zptmassweight+"(os>0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT>0.5)";
	cutsZllJesSys[iDY][iSys][iUp] = Weight+zptmassweight+"(os>0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT<0.5)";

	cutsZttSSJesSys[iDY][iSys][iUp] = Weight+zptmassweight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT>0.5)";
	cutsZllSSJesSys[iDY][iSys][iUp] = Weight+zptmassweight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT<0.5)";

	cutsZttSSrelaxedJesSys[iDY][iSys][iUp] = Weight+zptmassweight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT>0.5)";
	cutsZllSSrelaxedJesSys[iDY][iSys][iUp] = Weight+zptmassweight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+npartonCuts[iDY]+"&&isZTT<0.5)";

      }
    }

  }

  int nSamplesDY = 9;

  // filling histograms for DY samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples

    TFile * file = new TFile(directory+dySampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get(TauCheck);
    double norm = dyNorm[i];


    TString histNameZtt          = dySampleNames[i] + "_ztt_os";
    TString histNameZttSS        = dySampleNames[i] + "_ztt_ss";
    TString histNameZttSSrelaxed = dySampleNames[i] + "_ztt_ss_relaxed";
    TString histNameZll          = dySampleNames[i] + "_zll_os";
    TString histNameZllSS        = dySampleNames[i] + "_zll_ss";
    TString histNameZllSSrelaxed = dySampleNames[i] + "_zll_ss_relaxed";
    
    TH2D * histZtt2D          = new TH2D(histNameZtt,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histZttSS2D        = new TH2D(histNameZttSS,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histZttSSrelaxed2D = new TH2D(histNameZttSSrelaxed,"",nBinsX,binsX,nBinsY,binsY);

    TH2D * histZll2D          = new TH2D(histNameZll,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histZllSS2D        = new TH2D(histNameZllSS,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histZllSSrelaxed2D = new TH2D(histNameZllSSrelaxed,"",nBinsX,binsX,nBinsY,binsY);

    histZtt2D->Sumw2();
    histZttSS2D->Sumw2();
    histZttSSrelaxed2D->Sumw2();
    histZll2D->Sumw2();
    histZllSS2D->Sumw2();
    histZllSSrelaxed2D->Sumw2();

    tree->Draw(Variable+">>"+histNameZtt,  cutsZtt[i]);
    tree->Draw(Variable+">>"+histNameZttSS,cutsZttSS[i]);
    tree->Draw(Variable+">>"+histNameZttSSrelaxed,cutsZttSSrelaxed[i]);
    tree->Draw(Variable+">>"+histNameZll,  cutsZll[i]);
    tree->Draw(Variable+">>"+histNameZllSS,cutsZllSS[i]);
    tree->Draw(Variable+">>"+histNameZllSSrelaxed,cutsZllSSrelaxed[i]);

    histZtt[i]          = (TH1D*)Unfold(histZtt2D);
    histZttSS[i]        = (TH1D*)Unfold(histZttSS2D);
    histZttSSrelaxed[i] = (TH1D*)Unfold(histZttSSrelaxed2D);
    histZll[i]          = (TH1D*)Unfold(histZll2D);
    histZllSS[i]        = (TH1D*)Unfold(histZllSS2D);
    histZllSSrelaxed[i] = (TH1D*)Unfold(histZllSSrelaxed2D);

    std::cout << dySampleNames[i] << std::endl;

    // systematics
    for (int iSys=0; iSys<nSys; ++iSys ) {


      TH2D * histZttSys2D = new TH2D(histNameZtt+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      TH2D * histZllSys2D = new TH2D(histNameZll+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      histZttSys2D->Sumw2();
      histZllSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameZtt+sysName[iSys],cutsZttSys[i][iSys]);
      tree->Draw(VariableSys[iSys]+">>"+histNameZll+sysName[iSys],cutsZllSys[i][iSys]);
      histZttSys[i][iSys] = (TH1D*)Unfold(histZttSys2D);
      histZllSys[i][iSys] = (TH1D*)Unfold(histZllSys2D);


      TH2D * histZttSSSys2D = new TH2D(histNameZttSS+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      TH2D * histZllSSSys2D = new TH2D(histNameZllSS+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      histZttSSSys2D->Sumw2();
      histZllSSSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameZttSS+sysName[iSys],cutsZttSSSys[i][iSys]);
      tree->Draw(VariableSys[iSys]+">>"+histNameZllSS+sysName[iSys],cutsZllSSSys[i][iSys]);
      histZttSSSys[i][iSys] = (TH1D*)Unfold(histZttSSSys2D);
      histZllSSSys[i][iSys] = (TH1D*)Unfold(histZllSSSys2D);


      TH2D * histZttSSrelaxedSys2D = new TH2D(histNameZttSSrelaxed+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      TH2D * histZllSSrelaxedSys2D = new TH2D(histNameZllSSrelaxed+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      histZttSSrelaxedSys2D->Sumw2();
      histZllSSrelaxedSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameZttSSrelaxed+sysName[iSys],cutsZttSSrelaxedSys[i][iSys]);
      tree->Draw(VariableSys[iSys]+">>"+histNameZllSSrelaxed+sysName[iSys],cutsZllSSrelaxedSys[i][iSys]);
      histZttSSrelaxedSys[i][iSys] = (TH1D*)Unfold(histZttSSrelaxedSys2D);
      histZllSSrelaxedSys[i][iSys] = (TH1D*)Unfold(histZllSSrelaxedSys2D);

      std::cout << dySampleNames[i] << sysName[iSys] << std::endl;

    }

    // JES Systematics 
    for (int iSys=0; iSys<29; ++iSys ) {
      for (int iUp=0; iUp<2; ++iUp) {

	TH2D * histZttSys2D = new TH2D(histNameZtt+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	TH2D * histZllSys2D = new TH2D(histNameZll+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histZttSys2D->Sumw2();
	histZllSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZtt+jesSysName[iSys][iUp],cutsZttJesSys[i][iSys][iUp]);
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZll+jesSysName[iSys][iUp],cutsZllJesSys[i][iSys][iUp]);
	histZttJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZttSys2D);
	histZllJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZllSys2D);


	TH2D * histZttSSSys2D = new TH2D(histNameZttSS+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	TH2D * histZllSSSys2D = new TH2D(histNameZllSS+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histZttSSSys2D->Sumw2();
	histZllSSSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZttSS+jesSysName[iSys][iUp],cutsZttSSJesSys[i][iSys][iUp]);
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZllSS+jesSysName[iSys][iUp],cutsZllSSJesSys[i][iSys][iUp]);
	histZttSSJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZttSSSys2D);
	histZllSSJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZllSSSys2D);


	TH2D * histZttSSrelaxedSys2D = new TH2D(histNameZttSSrelaxed+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	TH2D * histZllSSrelaxedSys2D = new TH2D(histNameZllSSrelaxed+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histZttSSrelaxedSys2D->Sumw2();
	histZllSSrelaxedSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZttSSrelaxed+jesSysName[iSys][iUp],cutsZttSSrelaxedJesSys[i][iSys][iUp]);
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameZllSSrelaxed+jesSysName[iSys][iUp],cutsZllSSrelaxedJesSys[i][iSys][iUp]);
	histZttSSrelaxedJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZttSSrelaxedSys2D);
	histZllSSrelaxedJesSys[i][iSys][iUp] = (TH1D*)Unfold(histZllSSrelaxedSys2D);
	
	std::cout << dySampleNames[i] << jesSysName[iSys][iUp] << std::endl;

      }
    }

    for (int iB=1; iB<=nBins; ++iB) {

      double x = histZtt[i]->GetBinContent(iB);
      double e = histZtt[i]->GetBinError(iB);
      histZtt[i]->SetBinContent(iB,norm*x);
      histZtt[i]->SetBinError(iB,norm*e);

      x = histZttSS[i]->GetBinContent(iB);
      e = histZttSS[i]->GetBinError(iB);
      histZttSS[i]->SetBinContent(iB,norm*x);
      histZttSS[i]->SetBinError(iB,norm*e);

      x = histZttSSrelaxed[i]->GetBinContent(iB);
      e = histZttSSrelaxed[i]->GetBinError(iB);
      histZttSSrelaxed[i]->SetBinContent(iB,norm*x);
      histZttSSrelaxed[i]->SetBinError(iB,norm*e);

      x = histZll[i]->GetBinContent(iB);
      e = histZll[i]->GetBinError(iB);
      histZll[i]->SetBinContent(iB,norm*x);
      histZll[i]->SetBinError(iB,norm*e);

      x = histZllSS[i]->GetBinContent(iB);
      e = histZllSS[i]->GetBinError(iB);
      histZllSS[i]->SetBinContent(iB,norm*x);
      histZllSS[i]->SetBinError(iB,norm*e);

      x = histZllSSrelaxed[i]->GetBinContent(iB);
      e = histZllSSrelaxed[i]->GetBinError(iB);
      histZllSSrelaxed[i]->SetBinContent(iB,norm*x);
      histZllSSrelaxed[i]->SetBinError(iB,norm*e);
      
      for (int iSys=0; iSys<nSys; ++iSys) {

	x = histZttSys[i][iSys]->GetBinContent(iB);
	e = histZttSys[i][iSys]->GetBinError(iB);
	histZttSys[i][iSys]->SetBinContent(iB,norm*x);
	histZttSys[i][iSys]->SetBinError(iB,norm*e);

	x = histZllSys[i][iSys]->GetBinContent(iB);
	e = histZllSys[i][iSys]->GetBinError(iB);
	histZllSys[i][iSys]->SetBinContent(iB,norm*x);
	histZllSys[i][iSys]->SetBinError(iB,norm*e);


	x = histZttSSSys[i][iSys]->GetBinContent(iB);
	e = histZttSSSys[i][iSys]->GetBinError(iB);
	histZttSSSys[i][iSys]->SetBinContent(iB,norm*x);
	histZttSSSys[i][iSys]->SetBinError(iB,norm*e);

	x = histZllSSSys[i][iSys]->GetBinContent(iB);
	e = histZllSSSys[i][iSys]->GetBinError(iB);
	histZllSSSys[i][iSys]->SetBinContent(iB,norm*x);
	histZllSSSys[i][iSys]->SetBinError(iB,norm*e);


	x = histZttSSrelaxedSys[i][iSys]->GetBinContent(iB);
	e = histZttSSrelaxedSys[i][iSys]->GetBinError(iB);
	histZttSSrelaxedSys[i][iSys]->SetBinContent(iB,norm*x);
	histZttSSrelaxedSys[i][iSys]->SetBinError(iB,norm*e);

	x = histZllSSrelaxedSys[i][iSys]->GetBinContent(iB);
	e = histZllSSrelaxedSys[i][iSys]->GetBinError(iB);
	histZllSSrelaxedSys[i][iSys]->SetBinContent(iB,norm*x);
	histZllSSrelaxedSys[i][iSys]->SetBinError(iB,norm*e);

      }

      for (int iSys=0; iSys<29; ++iSys) {
	for (int iUp=0; iUp<2; ++iUp) {
	  
	  x = histZttJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZttJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZttJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZttJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	  
	  x = histZllJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZllJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZllJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZllJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	  

	  x = histZttSSJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZttSSJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZttSSJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZttSSJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);

	  x = histZllSSJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZllSSJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZllSSJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZllSSJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);


	  x = histZttSSrelaxedJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZttSSrelaxedJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZttSSrelaxedJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZttSSrelaxedJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);

	  x = histZllSSrelaxedJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histZllSSrelaxedJesSys[i][iSys][iUp]->GetBinError(iB);
	  histZllSSrelaxedJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histZllSSrelaxedJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	  
	}
      }

    }
    cout << dySampleNames[i] << " ZTT = " << histZtt[i]->GetSumOfWeights() 
	 << "    ZLL = " <<  histZll[i]->GetSumOfWeights() << endl;
  }

  hist[1]   = histZtt[0];
  histSS[1] = histZttSS[0];
  histSSrelaxed[1] = histZttSSrelaxed[0];

  hist[2]   = histZll[0];
  histSS[2] = histZllSS[0];
  histSSrelaxed[2] = histZllSSrelaxed[0];
  
  for (int iSys=0; iSys<nSys; ++iSys) {

    histSys[1][iSys] = histZttSys[0][iSys];
    histSys[2][iSys] = histZllSys[0][iSys];

    histSSSys[1][iSys] = histZttSSSys[0][iSys];
    histSSSys[2][iSys] = histZllSSSys[0][iSys];

    histSSrelaxedSys[1][iSys] = histZttSSrelaxedSys[0][iSys];
    histSSrelaxedSys[2][iSys] = histZllSSrelaxedSys[0][iSys];


  }

  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {

      histJesSys[1][iSys][iUp] = histZttJesSys[0][iSys][iUp];
      histJesSys[2][iSys][iUp] = histZllJesSys[0][iSys][iUp];
      
      histSSJesSys[1][iSys][iUp] = histZttSSJesSys[0][iSys][iUp];
      histSSJesSys[2][iSys][iUp] = histZllSSJesSys[0][iSys][iUp];
      
      histSSrelaxedJesSys[1][iSys][iUp] = histZttSSrelaxedJesSys[0][iSys][iUp];
      histSSrelaxedJesSys[2][iSys][iUp] = histZllSSrelaxedJesSys[0][iSys][iUp];
      
    }
  }

  for (int iDY=1; iDY<9; ++iDY) {
    hist[1]->Add(hist[1],histZtt[iDY]);
    hist[2]->Add(hist[2],histZll[iDY]);
    histSS[1]->Add(histSS[1],histZttSS[iDY]);
    histSS[2]->Add(histSS[2],histZllSS[iDY]);
    histSSrelaxed[1]->Add(histSSrelaxed[1],histZttSSrelaxed[iDY]);
    histSSrelaxed[2]->Add(histSSrelaxed[2],histZllSSrelaxed[iDY]);

    for (int iSys=0; iSys<nSys; ++iSys) {
      histSys[1][iSys]->Add(histSys[1][iSys],histZttSys[iDY][iSys]);
      histSys[2][iSys]->Add(histSys[2][iSys],histZllSys[iDY][iSys]);
      histSSSys[1][iSys]->Add(histSSSys[1][iSys],histZttSSSys[iDY][iSys]);
      histSSSys[2][iSys]->Add(histSSSys[2][iSys],histZllSSSys[iDY][iSys]);
      histSSrelaxedSys[1][iSys]->Add(histSSrelaxedSys[1][iSys],histZttSSrelaxedSys[iDY][iSys]);
      histSSrelaxedSys[2][iSys]->Add(histSSrelaxedSys[2][iSys],histZllSSrelaxedSys[iDY][iSys]);
    }

    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histJesSys[1][iSys][iUp]->Add(histJesSys[1][iSys][iUp],histZttJesSys[iDY][iSys][iUp]);
	histJesSys[2][iSys][iUp]->Add(histJesSys[2][iSys][iUp],histZllJesSys[iDY][iSys][iUp]);
	histSSJesSys[1][iSys][iUp]->Add(histSSJesSys[1][iSys][iUp],histZttSSJesSys[iDY][iSys][iUp]);
	histSSJesSys[2][iSys][iUp]->Add(histSSJesSys[2][iSys][iUp],histZllSSJesSys[iDY][iSys][iUp]);
	histSSrelaxedJesSys[1][iSys][iUp]->Add(histSSrelaxedJesSys[1][iSys][iUp],histZttSSrelaxedJesSys[iDY][iSys][iUp]);
	histSSrelaxedJesSys[2][iSys][iUp]->Add(histSSrelaxedJesSys[2][iSys][iUp],histZllSSrelaxedJesSys[iDY][iSys][iUp]);
      }
    }

  }


  // *******************************
  // ***** W+Jets samples *******
  // *******************************

  TH1D * histW[10];
  TH1D * histWSys[10][20];
  TH1D * histWJesSys[10][29][2];

  TH1D * histWSS[10];
  TH1D * histWSSSys[10][20];
  TH1D * histWSSJesSys[10][29][2]; 

  TH1D * histWSSrelaxed[10];
  TH1D * histWSSrelaxedSys[10][20];
  TH1D * histWSSrelaxedJesSys[10][29][2];

  // redefine reference cross sections
  // and reference samples

  refSamples[0] = "WJetsToLNu_13TeV-madgraphMLM";
  refSamples[1] = "W1JetsToLNu_13TeV-madgraphMLM";
  refSamples[2] = "W2JetsToLNu_13TeV-madgraphMLM";
  refSamples[3] = "W3JetsToLNu_13TeV-madgraphMLM";
  refSamples[4] = "W4JetsToLNu_13TeV-madgraphMLM";

  refXSec[0] = 61527;
  refXSec[1] = 1.221*9644.5;
  refXSec[2] = 1.221*3144.5;
  refXSec[3] = 1.221*954.8;
  refXSec[4] = 1.221*485.6;

  refEvents[0] = 0;
  refEvents[1] = 0;
  refEvents[2] = 0;
  refEvents[3] = 0;
  refEvents[4] = 0;

  for (int iDY=0; iDY<5; ++iDY) {
    TFile * file = new TFile(directory+refSamples[iDY]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
    refEvents[iDY] = histWeightsH->GetSumOfWeights();
  }
  TString wSampleNames[9] = {"WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "WJetsToLNu_13TeV-madgraphMLM",
			     "W1JetsToLNu_13TeV-madgraphMLM",
			     "W2JetsToLNu_13TeV-madgraphMLM",
			     "W3JetsToLNu_13TeV-madgraphMLM",
			     "W4JetsToLNu_13TeV-madgraphMLM"
  };

  double wNorm[9];
  wNorm[0] = lumi*refXSec[0]/refEvents[0];
  wNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  wNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  wNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  wNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);

  TString cutsW[9];
  TString cutsWSys[9][20];
  TString cutsWJesSys[9][29][2];

  TString cutsWSS[9];
  TString cutsWSSSys[9][20];
  TString cutsWSSJesSys[9][29][2];

  TString cutsWSSrelaxed[9];
  TString cutsWSSrelaxedSys[9][20];
  TString cutsWSSrelaxedJesSys[9][29][2];

  for (int iDY=0; iDY<9; ++iDY) {
    cutsW[iDY]   = Weight+"(os>0.5"+Cuts+npartonCuts[iDY]+")";
    cutsWSS[iDY] = Weight+qcdweight+"(os<0.5"+Cuts+npartonCuts[iDY]+")";
    cutsWSSrelaxed[iDY] = Weight+"(os<0.5"+CutsSS+npartonCuts[iDY]+")";
    for (int iSys=0; iSys<nSys; ++iSys) {
      cutsWSys[iDY][iSys] = Weight+"(os>0.5"+CutsSys[iSys]+npartonCuts[iDY]+")";
      cutsWSSSys[iDY][iSys] = Weight+qcdweight+"(os<0.5"+CutsSys[iSys]+npartonCuts[iDY]+")";
      cutsWSSrelaxedSys[iDY][iSys] = Weight+"(os<0.5"+CutsSSSys[iSys]+npartonCuts[iDY]+")";
    }
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	cutsWJesSys[iDY][iSys][iUp] = Weight+"(os>0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+")";
	cutsWSSJesSys[iDY][iSys][iUp] = Weight+qcdweight+"(os<0.5"+CutsJesSys[iSys][iUp]+npartonCuts[iDY]+")";
	cutsWSSrelaxedJesSys[iDY][iSys][iUp] = Weight+"(os<0.5"+CutsSSJesSys[iSys][iUp]+npartonCuts[iDY]+")";
      }
    }
  }

  // filling histograms for WJets samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples

    TFile * file = new TFile(directory+wSampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get(TauCheck);
    double norm = wNorm[i];

    TString histNameW   = wSampleNames[i] + Variable + "_w_os";
    TString histNameWSS = wSampleNames[i] + Variable + "_w_ss";
    TString histNameWSSrelaxed = wSampleNames[i] + Variable + "_w_ss_relaxed";

    TH2D * histW2D = new TH2D(histNameW,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histWSS2D = new TH2D(histNameWSS,"",nBinsX,binsX,nBinsY,binsY);
    TH2D * histWSSrelaxed2D = new TH2D(histNameWSSrelaxed,"",nBinsX,binsX,nBinsY,binsY);

    histW2D->Sumw2();
    histWSS2D->Sumw2();
    histWSSrelaxed2D->Sumw2();

    tree->Draw(Variable+">>"+histNameW,  cutsW[i]);
    tree->Draw(Variable+">>"+histNameWSS,cutsWSS[i]);
    tree->Draw(Variable+">>"+histNameWSSrelaxed,cutsWSSrelaxed[i]);

    histW[i]   = (TH1D*)Unfold(histW2D);
    histWSS[i] = (TH1D*)Unfold(histWSS2D);
    histWSSrelaxed[i] = (TH1D*)Unfold(histWSSrelaxed2D);

    std::cout << wSampleNames[i] << std::endl;

    // systematics
    for (int iSys=0; iSys<nSys; ++iSys ) {
      

      TH2D * histWSys2D = new TH2D(histNameW+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);

      histWSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameW+sysName[iSys],cutsWSys[i][iSys]);
      histWSys[i][iSys] = (TH1D*)Unfold(histWSys2D);

      TH2D * histWSSSys2D = new TH2D(histNameWSS+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      histWSSSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameWSS+sysName[iSys],cutsWSSSys[i][iSys]);
      histWSSSys[i][iSys] = (TH1D*)Unfold(histWSSSys2D);

      TH2D * histWSSrelaxedSys2D = new TH2D(histNameWSSrelaxed+sysName[iSys],"",nBinsX,binsX,nBinsY,binsY);
      histWSSrelaxedSys2D->Sumw2();
      tree->Draw(VariableSys[iSys]+">>"+histNameWSSrelaxed+sysName[iSys],cutsWSSrelaxedSys[i][iSys]);
      histWSSrelaxedSys[i][iSys] = (TH1D*)Unfold(histWSSrelaxedSys2D);

      std::cout << wSampleNames[i] << sysName[iSys] << std::endl;

    }

    // JES systematics
    for (int iSys=0; iSys<29; ++iSys ) {
      for (int iUp=0; iUp<2; ++iUp) {

	TH2D * histWSys2D = new TH2D(histNameW+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histWSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameW+jesSysName[iSys][iUp],cutsWJesSys[i][iSys][iUp]);
	histWJesSys[i][iSys][iUp] = (TH1D*)Unfold(histWSys2D);
	
	TH2D * histWSSSys2D = new TH2D(histNameWSS+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histWSSSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameWSS+jesSysName[iSys][iUp],cutsWSSJesSys[i][iSys][iUp]);
	histWSSJesSys[i][iSys][iUp] = (TH1D*)Unfold(histWSSSys2D);
	
	TH2D * histWSSrelaxedSys2D = new TH2D(histNameWSSrelaxed+jesSysName[iSys][iUp],"",nBinsX,binsX,nBinsY,binsY);
	histWSSrelaxedSys2D->Sumw2();
	tree->Draw(VariableJesSys[iSys][iUp]+">>"+histNameWSSrelaxed+jesSysName[iSys][iUp],cutsWSSrelaxedJesSys[i][iSys][iUp]);
	histWSSrelaxedJesSys[i][iSys][iUp] = (TH1D*)Unfold(histWSSrelaxedSys2D);

	std::cout << wSampleNames[i] << jesSysName[iSys][iUp] << std::endl;

      }
    }

    for (int iB=1; iB<=nBins; ++iB) {

      double x = histW[i]->GetBinContent(iB);
      double e = histW[i]->GetBinError(iB);
      histW[i]->SetBinContent(iB,norm*x);
      histW[i]->SetBinError(iB,norm*e);

      x = histWSS[i]->GetBinContent(iB);
      e = histWSS[i]->GetBinError(iB);
      histWSS[i]->SetBinContent(iB,norm*x);
      histWSS[i]->SetBinError(iB,norm*e);

      x = histWSSrelaxed[i]->GetBinContent(iB);
      e = histWSSrelaxed[i]->GetBinError(iB);
      histWSSrelaxed[i]->SetBinContent(iB,norm*x);
      histWSSrelaxed[i]->SetBinError(iB,norm*e);
      
      for (int iSys=0; iSys<nSys; ++iSys) {

	x = histWSys[i][iSys]->GetBinContent(iB);
	e = histWSys[i][iSys]->GetBinError(iB);
	histWSys[i][iSys]->SetBinContent(iB,norm*x);
	histWSys[i][iSys]->SetBinError(iB,norm*e);

	x = histWSSSys[i][iSys]->GetBinContent(iB);
	e = histWSSSys[i][iSys]->GetBinError(iB);
	histWSSSys[i][iSys]->SetBinContent(iB,norm*x);
	histWSSSys[i][iSys]->SetBinError(iB,norm*e);

	x = histWSSrelaxedSys[i][iSys]->GetBinContent(iB);
	e = histWSSrelaxedSys[i][iSys]->GetBinError(iB);
	histWSSrelaxedSys[i][iSys]->SetBinContent(iB,norm*x);
	histWSSrelaxedSys[i][iSys]->SetBinError(iB,norm*e);

      }

      for (int iSys=0; iSys<29; ++iSys) {
	for (int iUp=0; iUp<2; ++iUp) {

	  x = histWJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histWJesSys[i][iSys][iUp]->GetBinError(iB);
	  histWJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histWJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);
	  
	  x = histWSSJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histWSSJesSys[i][iSys][iUp]->GetBinError(iB);
	  histWSSJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histWSSJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);

	  x = histWSSrelaxedJesSys[i][iSys][iUp]->GetBinContent(iB);
	  e = histWSSrelaxedJesSys[i][iSys][iUp]->GetBinError(iB);
	  histWSSrelaxedJesSys[i][iSys][iUp]->SetBinContent(iB,norm*x);
	  histWSSrelaxedJesSys[i][iSys][iUp]->SetBinError(iB,norm*e);

	}
      }


    }
    std::cout << wSampleNames[i] << " = " << histW[i]->GetSumOfWeights() << std::endl;
    //    delete file;
  }

  hist[5]   = histW[0];
  histSS[5] = histWSS[0];
  histSSrelaxed[5] = histWSSrelaxed[0];
  
  for (int iSys=0; iSys<nSys; ++iSys) {
    histSys[5][iSys] = histWSys[0][iSys];
    histSSSys[5][iSys] = histWSSSys[0][iSys];
    histSSrelaxedSys[5][iSys] = histWSSrelaxedSys[0][iSys];
  }

  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {
      histJesSys[5][iSys][iUp] = histWJesSys[0][iSys][iUp];
      histSSJesSys[5][iSys][iUp] = histWSSJesSys[0][iSys][iUp];
      histSSrelaxedJesSys[5][iSys][iUp] = histWSSrelaxedJesSys[0][iSys][iUp];
    }
  }


  for (int iDY=1; iDY<9; ++iDY) {
    hist[5]->Add(hist[5],histW[iDY]);
    histSS[5]->Add(histSS[5],histWSS[iDY]);
    histSSrelaxed[5]->Add(histSSrelaxed[5],histWSSrelaxed[iDY]);
    for (int iSys=0; iSys<nSys; ++iSys) {
      histSys[5][iSys]->Add(histSys[5][iSys],histWSys[iDY][iSys]);
      histSSSys[5][iSys]->Add(histSSSys[5][iSys],histWSSSys[iDY][iSys]);
      histSSrelaxedSys[5][iSys]->Add(histSSrelaxedSys[5][iSys],histWSSrelaxedSys[iDY][iSys]);
    }
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histJesSys[5][iSys][iUp]->Add(histJesSys[5][iSys][iUp],histWJesSys[iDY][iSys][iUp]);
	histSSJesSys[5][iSys][iUp]->Add(histSSJesSys[5][iSys][iUp],histWSSJesSys[iDY][iSys][iUp]);
	histSSrelaxedJesSys[5][iSys][iUp]->Add(histSSrelaxedJesSys[5][iSys][iUp],histWSSrelaxedJesSys[iDY][iSys][iUp]);
      }
    }
  }

  // ********************************
  // ********** FINAL PART **********
  // ********************************


  // Dreall-Yan low mass + high mass
  hist[1]->Add(hist[1],hist[3]); 
  hist[2]->Add(hist[2],hist[4]); 
  for (int iSys=0; iSys<nSys; ++iSys) {
    histSys[1][iSys]->Add(histSys[1][iSys],histSys[3][iSys]);
    histSys[2][iSys]->Add(histSys[2][iSys],histSys[4][iSys]);
  }
  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {
      histJesSys[1][iSys][iUp]->Add(histJesSys[1][iSys][iUp],histJesSys[3][iSys][iUp]);
      histJesSys[2][iSys][iUp]->Add(histJesSys[2][iSys][iUp],histJesSys[4][iSys][iUp]);
    }
  }

  //  adding up single top and VV backgrounds
  for (int iH=8; iH<19; ++iH) {
    hist[7]->Add(hist[7],hist[iH]);
    for (int iSys=0; iSys<nSys; ++iSys)
      histSys[7][iSys]->Add(histSys[7][iSys],histSys[iH][iSys]);
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histJesSys[7][iSys][iUp]->Add(histJesSys[7][iSys][iUp],histJesSys[iH][iSys][iUp]);
      }
    }
  }

  // adding up W+Jets and W+gamma samples
  for (int iH=19; iH<22; ++iH) {
    hist[5]->Add(hist[5],hist[iH]);
    for (int iSys=0; iSys<nSys; ++iSys)
      histSys[5][iSys]->Add(histSys[5][iSys],histSys[iH][iSys]);
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histJesSys[5][iSys][iUp]->Add(histJesSys[5][iSys][iUp],histJesSys[iH][iSys][iUp]);
      }
    }
  }
  for (int iH=23; iH<25; ++iH) {
    hist[5]->Add(hist[5],hist[iH]);
    for (int iSys=0; iSys<nSys; ++iSys)
      histSys[5][iSys]->Add(histSys[5][iSys],histSys[iH][iSys]);
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
        histJesSys[5][iSys][iUp]->Add(histJesSys[5][iSys][iUp],histJesSys[iH][iSys][iUp]);
      } 
    }
  }

  std::cout << "Data in SS region = " << histSS[0]->GetSumOfWeights() << std::endl;
 
  // subtracting background from SS
  for (int iH=1; iH<27; ++iH) {
    histSS[0]->Add(histSS[0],histSS[iH],1,-1);
    histSSrelaxed[0]->Add(histSSrelaxed[0],histSSrelaxed[iH],1,-1);
    for (int iSys=0; iSys<nSys; ++iSys) {
      histSSSys[0][iSys]->Add(histSSSys[0][iSys],histSSSys[iH][iSys],1,-1);
      histSSrelaxedSys[0][iSys]->Add(histSSrelaxedSys[0][iSys],histSSrelaxedSys[iH][iSys],1,-1);
    }
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histSSJesSys[0][iSys][iUp]->Add(histSSJesSys[0][iSys][iUp],histSSJesSys[iH][iSys][iUp],1,-1);
	histSSrelaxedJesSys[0][iSys][iUp]->Add(histSSrelaxedJesSys[0][iSys][iUp],histSSrelaxedJesSys[iH][iSys][iUp],1,-1);
      }
    }
  }

  // ensuring that there are no negative bins in SSrelaxed
  for (int iB=1; iB<=nBins; ++iB) {
    float xSSrelaxed = histSSrelaxed[0]->GetBinError(iB);
    if (xSSrelaxed<0)
      histSSrelaxed[0]->SetBinContent(iB,0.);
    for (int iSys=0; iSys<nSys; ++iSys) {
      xSSrelaxed = histSSrelaxedSys[0][iSys]->GetBinContent(iB);
      if (xSSrelaxed<0)
	histSSrelaxedSys[0][iSys]->SetBinContent(iB,0);
    }
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	xSSrelaxed = histSSrelaxedJesSys[0][iSys][iUp]->GetBinContent(iB);
	if (xSSrelaxed<0)
	  histSSrelaxedJesSys[0][iSys][iUp]->SetBinContent(iB,0);
      }
    }
  }

  // normalizing QCD templates
  float qcdNorm = histSS[0]->GetSumOfWeights();
  float qcdNormRelaxed = histSSrelaxed[0]->GetSumOfWeights(); 
  float ratioQCD = qcdNorm / qcdNormRelaxed;

  std::cout << "qcdNorm = " << qcdNorm << "   QCDnormRelaxed = " << qcdNormRelaxed << "   ratio = " << ratioQCD << std::endl;

  float qcdNormE2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xSSrelaxed = histSSrelaxed[0]->GetBinContent(iB);
    float xSSErelaxed = histSSrelaxed[0]->GetBinError(iB);
    float error = histSS[0]->GetBinError(iB);
    qcdNormE2 += error*error;
    float xHist = ratioQCD*histSSrelaxed[0]->GetBinContent(iB);
    float eHist = ratioQCD*histSSrelaxed[0]->GetBinError(iB);
    histSSrelaxed[0]->SetBinContent(iB,xHist);
    histSSrelaxed[0]->SetBinError(iB,eHist);
  }

  std::cout << "Cross-check : qcdNorm = " << histSS[0]->GetSumOfWeights() 
	    << "   QCDnormRelaxed = " << histSSrelaxed[0]->GetSumOfWeights() 
	    << "   ratio = " << histSS[0]->GetSumOfWeights()/histSSrelaxed[0]->GetSumOfWeights() << std::endl;

  for (int iSys=0; iSys<nSys; ++iSys) {
    qcdNorm = histSSSys[0][iSys]->GetSumOfWeights();
    qcdNormRelaxed = histSSrelaxedSys[0][iSys]->GetSumOfWeights();
    ratioQCD = qcdNorm / qcdNormRelaxed;
    for (int iB=1; iB<=nBins; ++iB) {
      float xHist = ratioQCD*histSSrelaxedSys[0][iSys]->GetBinContent(iB);
      float eHist = ratioQCD*histSSrelaxedSys[0][iSys]->GetBinError(iB);
      histSSrelaxedSys[0][iSys]->SetBinContent(iB,xHist);
      histSSrelaxedSys[0][iSys]->SetBinError(iB,eHist);
    }
  }

  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {
      qcdNorm = histSSJesSys[0][iSys][iUp]->GetSumOfWeights();
      qcdNormRelaxed = histSSrelaxedJesSys[0][iSys][iUp]->GetSumOfWeights();
      ratioQCD = qcdNorm / qcdNormRelaxed;
      for (int iB=1; iB<=nBins; ++iB) {
	float xHist = ratioQCD*histSSrelaxedJesSys[0][iSys][iUp]->GetBinContent(iB);
	float eHist = ratioQCD*histSSrelaxedJesSys[0][iSys][iUp]->GetBinError(iB);
	histSSrelaxedJesSys[0][iSys][iUp]->SetBinContent(iB,xHist);
	histSSrelaxedJesSys[0][iSys][iUp]->SetBinError(iB,eHist);
      }
    }
  }

  // ******************************************************************
  // * applying corrections to DY in slices of the second 2D variable *
  // ******************************************************************
  histSys[1][10] = (TH1D*)hist[1]->Clone("ZTT_2DUp");
  histSys[1][11] = (TH1D*)hist[1]->Clone("ZTT_2DDown");
  histSys[2][10] = (TH1D*)hist[2]->Clone("ZL_2DUp");
  histSys[2][11] = (TH1D*)hist[2]->Clone("ZL_2DDown");
  if (category=="em_vbf") {
    for (int iB=1; iB<=nBins; ++iB) {
      int binN = int((iB-1)/nBinsX);
      float dySF = dyCorr2D[binN];
      for (int iDY=1; iDY<=2; ++iDY) {
	float xDY = hist[iDY]->GetBinContent(iB);
	float eDY = hist[iDY]->GetBinError(iB);
	hist[iDY]->SetBinContent(iB,dySF*xDY);
	hist[iDY]->SetBinError(iB,dySF*eDY);
	histSys[iDY][10]->SetBinContent(iB,dySF*dySF*xDY);
	histSys[iDY][10]->SetBinError(iB,dySF*dySF*eDY);
	histSys[iDY][11]->SetBinContent(iB,xDY);
	histSys[iDY][11]->SetBinError(iB,eDY);
	for (int iSys=0; iSys<nSys; ++iSys) {
	  float xDYsys = histSys[iDY][iSys]->GetBinContent(iB);
	  float eDYsys = histSys[iDY][iSys]->GetBinError(iB);
	  histSys[iDY][iSys]->SetBinContent(iB,dySF*xDYsys);
	  histSys[iDY][iSys]->SetBinError(iB,dySF*eDYsys);
	}
	for (int iSys=0; iSys<29; ++iSys) {
	  for (int iUp=0; iUp<2; ++iUp) {
	    float xDYsys = histJesSys[iDY][iSys][iUp]->GetBinContent(iB);
	    float eDYsys = histJesSys[iDY][iSys][iUp]->GetBinError(iB);
	    histJesSys[iDY][iSys][iUp]->SetBinContent(iB,dySF*xDYsys);
	    histJesSys[iDY][iSys][iUp]->SetBinError(iB,dySF*eDYsys);
	  }
	}
      }
    }
  }

  // ********************************************
  // ********** adding ap WplusH and WminusH ****
  // ********************************************
  for (int iM=0; iM<5; ++iM) {
    histSignal[3][iM]->Add(histSignal[3][iM],histSignal[4][iM]);
    for (int iSys=0; iSys<nSys; ++iSys)
      histSignalSys[3][iM][iSys]->Add(histSignalSys[3][iM][iSys],histSignalSys[4][iM][iSys]);
    for (int iSys=0; iSys<29; ++iSys) {
      for (int iUp=0; iUp<2; ++iUp) {
	histSignalJesSys[3][iM][iSys][iUp]->Add(histSignalJesSys[3][iM][iSys][iUp],histSignalJesSys[4][iM][iSys][iUp]);
      }
    }
  }

  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * QCD = (TH1D*)histSSrelaxed[0]->Clone("QCD");
  TH1D * ZTT = (TH1D*)hist[1]->Clone("ZTT");
  TH1D * ZLL = (TH1D*)hist[2]->Clone("ZLL");
  TH1D * W   = (TH1D*)hist[5]->Clone("W");
  TH1D * TT  = (TH1D*)hist[6]->Clone("TT");
  TH1D * VV  = (TH1D*)hist[7]->Clone("VV");
  TH1D * EWKZ = (TH1D*)hist[22]->Clone("EWKZ");
  TH1D * ggHWW125 = (TH1D*)hist[25]->Clone("HWW_gg125");
  TH1D * qqHWW125 = (TH1D*)hist[26]->Clone("HWW_qq125");

  TString BaseName = "htt_em.inputs-sm-13TeV_" +Suffix;
  TString rootFileName = BaseName+".root";
  TFile * fileInputs = new TFile(rootFileName,"recreate"); 
  fileInputs->mkdir(category);
  fileInputs->cd(category);
  histData->Write("data_obs");
  TT->Write("TT");
  ZTT->Write("ZTT");
  ZLL->Write("ZL");
  W->Write("W");
  QCD->Write("QCD");
  VV->Write("VV");
  EWKZ->Write("EWKZ");
  ggHWW125->Write("HWW_gg125");
  qqHWW125->Write("HWW_qq125");

  for (int iSys=2; iSys<6; ++iSys) {
    histSSrelaxedSys[0][iSys]->Write("QCD"+sysName[iSys]);
    histSys[1][iSys]->Write("ZTT"+sysName[iSys]);
    histSys[2][iSys]->Write("ZL"+sysName[iSys]);
    histSys[5][iSys]->Write("W"+sysName[iSys]);
    histSys[6][iSys]->Write("TT"+sysName[iSys]);
    histSys[7][iSys]->Write("VV"+sysName[iSys]);
    histSys[22][iSys]->Write("EWKZ"+sysName[iSys]);
    histSys[25][iSys]->Write("HWW_gg125"+sysName[iSys]);
    histSys[26][iSys]->Write("HWW_qq125"+sysName[iSys]);
  }
  for (int iSys=0; iSys<2; ++iSys) {
    histSys[6][iSys]->Write("TT"+sysName[iSys]);
  }
  for (int iSys=6; iSys<7; ++iSys) {
    histSys[1][iSys]->Write("ZTT"+sysName[iSys]);
    histSys[2][iSys]->Write("ZL"+sysName[iSys]);
  }
  for (int iSys=10; iSys<11; ++iSys) {
    histSys[1][iSys]->Write("ZTT"+sysName[iSys]);
    histSys[2][iSys]->Write("ZL"+sysName[iSys]);
  }
  for (int iSys=0; iSys<29; ++iSys) {
    for (int iUp=0; iUp<2; ++iUp) {
      histSSrelaxedJesSys[0][iSys][iUp]->Write("QCD"+jesSysName[iSys][iUp]);
      histJesSys[1][iSys][iUp]->Write("ZTT"+jesSysName[iSys][iUp]);
      histJesSys[2][iSys][iUp]->Write("ZL"+jesSysName[iSys][iUp]);
      histJesSys[5][iSys][iUp]->Write("W"+jesSysName[iSys][iUp]);
      histJesSys[6][iSys][iUp]->Write("TT"+jesSysName[iSys][iUp]);
      histJesSys[7][iSys][iUp]->Write("VV"+jesSysName[iSys][iUp]);
      histJesSys[22][iSys][iUp]->Write("EWKZ"+jesSysName[iSys][iUp]);
      histJesSys[25][iSys][iUp]->Write("HWW_gg125"+jesSysName[iSys][iUp]);
      histJesSys[26][iSys][iUp]->Write("HWW_qq125"+jesSysName[iSys][iUp]);
    }
  }

  TString signalTemplates[4] = {"ggH","qqH","WH","ZH"};
  for (int iSig=0; iSig<4; ++iSig) {
    for (int iM=0; iM<5; ++iM) {
      TString sigTemplate = signalTemplates[iSig] + massH[iM];
      histSignal[iSig][iM]->Write(sigTemplate);
      for (int iSys=2; iSys<6; ++iSys) 
	histSignalSys[iSig][iM][iSys]->Write(sigTemplate+sysName[iSys]);
      if (iSig==0) {
	histSignalSys[iSig][iM][8]->Write(sigTemplate+sysName[8]);
	histSignalSys[iSig][iM][9]->Write(sigTemplate+sysName[9]);
      }
      for (int iSys=0; iSys<29; ++iSys) {
	for (int iUp=0; iUp<2; ++iUp) {
	  histSignalJesSys[iSig][iM][iSys][iUp]->Write(sigTemplate+jesSysName[iSys][iUp]);
	}
      }

    }
  }

  fileInputs->Close();

}
