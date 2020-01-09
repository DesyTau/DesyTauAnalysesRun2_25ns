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
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

float topPtWeight(float pt1,
		  float pt2) {

  float a = 0.156;    // Run1 a parameter
  float b = -0.00137;  // Run1 b parameter

  if (pt1>400) pt1 = 400;
  if (pt2>400) pt2 = 400;

  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  return TMath::Sqrt(w1*w2);

}

TVector3 get_refitted_PV_with_BS(const AC1B * analysisTree, int leptonIndex1, int leptonIndex2, bool &is_refitted_PV_with_BS){
	float vtx_x = analysisTree->primvertexwithbs_x; // by default store non-refitted PV with BS constraint if refitted one is not found
	float vtx_y = analysisTree->primvertexwithbs_y;
	float vtx_z = analysisTree->primvertexwithbs_z;
	is_refitted_PV_with_BS = false;

	for(unsigned int i = 0; i < analysisTree->refitvertexwithbs_count; i++)
	{
	  if( (leptonIndex1 == analysisTree->refitvertexwithbs_muIndex[i][0] || leptonIndex1 == analysisTree->refitvertexwithbs_muIndex[i][1]) &&
	      (leptonIndex2 == analysisTree->refitvertexwithbs_muIndex[i][0] || leptonIndex2 == analysisTree->refitvertexwithbs_muIndex[i][1]))
	    {
	      vtx_x = analysisTree->refitvertexwithbs_x[i];
	      vtx_y = analysisTree->refitvertexwithbs_y[i];

	      vtx_z = analysisTree->refitvertexwithbs_z[i];
	      is_refitted_PV_with_BS = true;
	    }
	}
	TVector3 vertex_coord(vtx_x, vtx_y, vtx_z);
	return vertex_coord;
}

TVector3 ipVec_Lepton(const AC1B * analysisTree, int muonIndex, TVector3 vertex) {

  TVector3 secvertex(0.,0.,0.);
  TVector3 momenta(0.,0.,0.);    

  secvertex.SetXYZ(analysisTree->muon_vx[muonIndex], 
		   analysisTree->muon_vy[muonIndex],
		   analysisTree->muon_vz[muonIndex]);
  
   
  momenta.SetXYZ(analysisTree->muon_px[muonIndex],
		 analysisTree->muon_py[muonIndex],
		 analysisTree->muon_pz[muonIndex]);
   
  TVector3 r(0.,0.,0.);
  r=secvertex-vertex;
  
  double projection=r*momenta/momenta.Mag2();
  TVector3 IP;    
  IP=r-momenta*projection;
  
  return IP;
};

/*

float topPtWeight(float pt1,
		  float pt2) {
    
  if (pt1>400) pt1 = 400;
  if (pt2>400) pt2 = 400;
    
  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter
    
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
    
  return TMath::Sqrt(w1*w2);  
}

*/

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

bool isICHEPmed(unsigned int Index, AC1B analysisTree) {
        bool goodGlob = 
	  analysisTree.muon_isGlobal[Index] && 
	  analysisTree.muon_normChi2[Index] < 3 && 
	  analysisTree.muon_combQ_chi2LocalPosition[Index] < 12 && 
	  analysisTree.muon_combQ_trkKink[Index] < 20;
        bool isICHEPmedium  = 
	  analysisTree.muon_isLoose[Index] &&
	  analysisTree.muon_validFraction[Index] >0.49 &&
	  analysisTree.muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
        return isICHEPmedium;
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
  const bool applyKITCorrection = cfg.get<bool>("ApplyKITCorrection");
  const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  const bool applyRecoilOnGenerator = cfg.get<bool>("ApplyRecoilOnGenerator");
  const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
  const bool applyTopPtReweighting = cfg.get<bool>("ApplyTopPtReweighting");
  const bool applyZMassPtReweighting = cfg.get<bool>("ApplyZMassPtReweighting");
  const bool interpolateZMassPtWeight = cfg.get<bool>("InterpolateZMassPtWeight");

  // kinematic cuts on muons
  const float ptMuonCut        = cfg.get<float>("ptMuonCut");
  const float ptMuonProbeCut   = cfg.get<float>("ptMuonProbeCut");
  const float etaMuonCut       = cfg.get<float>("etaMuonCut");
  const float etaMuonProbeCut  = cfg.get<float>("etaMuonProbeCut");
  const float dxyMuonCut       = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut        = cfg.get<float>("dzMuonCut");
  const float dxyMuonProbeCut  = cfg.get<float>("dxyMuonProbeCut");
  const float dzMuonProbeCut   = cfg.get<float>("dzMuonProbeCut");
  const float isoMuonCut       = cfg.get<float>("isoMuonCut");
  const bool  applyTauTauSelection = cfg.get<bool>("ApplyTauTauSelection");
  const bool  selectZToTauTauMuMu = cfg.get<bool>("SelectZToTauTauMuMu");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");
  const float dimuonMassCut = cfg.get<float>("DimuonMassCut");
  const bool isoDR03         = cfg.get<bool>("IsoDR03");

  // trigger
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
  const string muonTriggerName  = cfg.get<string>("MuonTriggerName");
  const string muonFilterName   = cfg.get<string>("MuonFilterName");

  TString MuonTriggerName(muonTriggerName);
  TString MuonFilterName(muonFilterName);

  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // jet related cuts
  const float jetEtaCut      = cfg.get<float>("JetEtaCut");
  const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");
  const float jetPtHighCut   = cfg.get<float>("JetPtHighCut");
  const float jetPtLowCut    = cfg.get<float>("JetPtLowCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  // vertex distributions filenames and histname
  const string vertDataFileName = cfg.get<string>("VertexDataFileName");
  const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
  const string vertHistName     = cfg.get<string>("VertexHistName");

  const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
  const string MuonTrigFile  = cfg.get<string>("MuonTrigEff"); 
  const string correctionWorkspaceFile = cfg.get<string>("CorrectionWorkspaceFile");
  TString CorrectionWorkspaceFile(correctionWorkspaceFile);

  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  // Z (mass,pt) reweighting
  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

  const bool applyIpCorrection = cfg.get<bool>("ApplyIpCorrection"); 
  const string ipCorrFileName  = cfg.get<string>("IpCorrectionFileName");
  TString IpCorrFileName(ipCorrFileName);

  const string jsonFile = cfg.get<string>("jsonFile");
  const int applyJES = cfg.get<int>("applyJES");
  const float muonMomScale = cfg.get<float>("MuonMomentumScale");

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");
  const string pileUpMCHist = cfg.get<string>("PileUpMCHistName");


  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);
  TString PileUpMCHist(pileUpMCHist);

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

  TH1::SetDefaultSumw2(true);

  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
  TH1D * histWeightsSkimmedH = new TH1D("histWeightsSkimmedH","",1,-0.5,0.5);

  TH1D * massZH = new TH1D("massZH","",1000,0,1000);
  TH1D * ptZH = new TH1D("ptZH","",1000,0,1000);
  TH2D * massPtZH = new TH2D("massPtZH","",100,0,1000,100,0,1000);

  // Histograms after final selection
  TH1D * ptLeadingMuSelH = new TH1D("ptLeadingMuSelH","",100,0,200);
  TH1D * ptTrailingMuSelH = new TH1D("ptTrailingMuSelH","",100,0,200);
  TH1D * etaLeadingMuSelH = new TH1D("etaLeadingMuSelH","",50,-2.5,2.5);
  TH1D * etaTrailingMuSelH = new TH1D("etaTrailingMuSelH","",50,-2.5,2.5);
  TH1D * massSelH = new TH1D("massSelH","",200,0,200);
  TH1D * massExtendedSelH =  new TH1D("massExtendedSelH","",500,0,5000);
  

  TH2D * dimuonMassPtH = new TH2D("dimuonMassPtH","",100,0,1000,100,0,1000);

  TH1D * dimuonPtSelH = new TH1D("dimuonPtSelH","",100,0,1000);
  TH1D * dimuonEtaSelH = new TH1D("dimuonEtaSelH","",120,-6,6);
  TH1D * metSelH  = new TH1D("metSelH","",200,0,400);
  TH1D * puppimetSelH = new TH1D("puppimetSelH","",200,0,400);

  TH1D * metZSelH  = new TH1D("metZSelH","",200,0,400);
  TH1D * puppimetZSelH = new TH1D("puppimetZSelH","",200,0,400);

  TH1D * metTopSelH  = new TH1D("metTopSelH","",200,0,400);
  TH1D * puppimetTopSelH = new TH1D("puppimetTopSelH","",200,0,400);

  TH1D * numberOfVerticesH = new TH1D("numberOfVerticesH","",100,-0.5,99.5);

  // Event categories
  TH1D * mass0jetH = new TH1D("mass0jetH","",200,0,200);
  TH1D * mass0jetBvetoH = new TH1D("mass0jetBvetoH","",200,0,200);

  TH1D * massBoostedH = new TH1D("massBoostedH","",200,0,200);
  TH1D * massBoostedBvetoH = new TH1D("massBoostedBvetoH","",200,0,200);

  TH1D * massVBFH = new TH1D("massVBFH","",200,0,200);
  TH1D * massVBFBvetoH = new TH1D("massVBFBvetoH","",200,0,200);

  // Applying selection and mass cut

  TH1D * nJets20BTagMedium0jetH = new TH1D("nJets20BTagMedium0jetH","",11,-0.5,10.5);
  TH1D * nJets20BTagMediumBoostedH = new TH1D("nJets20BTagMediumBoostedH","",11,-0.5,10.5);
  TH1D * nJets20BTagMediumVBFH = new TH1D("nJets20BTagMediumVBFH","",11,-0.5,10.5);

  TH1D * mjjVBFH = new TH1D("mjjVBFH","",30,0,3000);
  TH1D * mjjVBFBvetoH = new TH1D("mjjVBFBvetoH","",30,0,3000);

  TH1D * dimuonPtBoostedH = new TH1D("dimuonPtBoostedH","",100,0,1000);
  TH1D * dimuonPtBoostedBvetoH = new TH1D("dimuonPtBoostedBvetoH","",100,0,1000);

  TH1D * leadingJetPtH  = new TH1D("leadingJetPtH","",50,0,500);
  TH1D * leadingJetEtaH = new TH1D("leadingJetEtaH","",100,-5,5);
  TH1D * leadingJetPhiH = new TH1D("leadingJetPhiH","",100,-TMath::Pi(),TMath::Pi());
  
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


  TH1D * nJets30SelH    = new TH1D("nJets30SelH","",11,-0.5,10.5);
  TH1D * nJets30etaCutSelH = new TH1D("nJets30etaCutSelH","",11,-0.5,10.5);
  TH1D * nJets20SelH    = new TH1D("nJets20SelH","",11,-0.5,10.5);
  TH1D * nJets20etaCutSelH = new TH1D("nJets20etaCutSelH","",11,-0.5,10.5);

  TH1D * nJets20BTagLooseSelH    = new TH1D("nJets20BTagLooseSelH","",11,-0.5,10.5);
  TH1D * nJets20BTagMediumSelH   = new TH1D("nJets20BTagMediumSelH","",11,-0.5,10.5); 

  TH1D * nJets30ZSelH    = new TH1D("nJets30ZSelH","",11,-0.5,10.5);
  TH1D * nJets30etaCutZSelH = new TH1D("nJets30etaCutZSelH","",11,-0.5,10.5);
  TH1D * nJets20ZSelH    = new TH1D("nJets20ZSelH","",11,-0.5,10.5);
  TH1D * nJets20etaCutZSelH = new TH1D("nJets20etaCutZSelH","",11,-0.5,10.5);

  TH1D * nJets20BTagLooseZSelH    = new TH1D("nJets20BTagLooseZSelH","",11,-0.5,10.5);
  TH1D * nJets20BTagMediumZSelH   = new TH1D("nJets20BTagMediumZSelH","",11,-0.5,10.5); 

  TH1D * HT30SelH       = new TH1D("HT30SelH","",50,0,500);
  TH1D * HT30etaCutSelH = new TH1D("HT30etaCutSelH","",50,0,500);
  TH1D * HT20SelH       = new TH1D("HT20SelH","",50,0,500);
  TH1D * HT20etaCutSelH = new TH1D("HT20etaCutSelH","",50,0,500);
  
  TH1D * NumberOfVerticesH = new TH1D("NumberOfVerticesH","",51,-0.5,50.5);

  TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);
  TH1D * MuSF_IdIso_Mu2H = new TH1D("MuIdIsoSF_Mu2H", "MuIdIsoSF_Mu2", 100, 0.5,1.5);

  TH1D * ipxH = new TH1D("ipxH","ipxH",100,-0.01,0.01);
  TH1D * ipyH = new TH1D("ipyH","ipyH",100,-0.01,0.01);
  TH1D * ipzH = new TH1D("ipzH","ipzH",100,-0.01,0.01);

  TH1D * ipdxH = new TH1D("ipdxH","ipdxH",100,-0.02,0.02);
  TH1D * ipdyH = new TH1D("ipdyH","ipdyH",100,-0.02,0.02);
  TH1D * ipdzH = new TH1D("ipdzH","ipdzH",100,-0.02,0.02);

  TH1D * nxyipH = new TH1D("nxyipH","nxyipH",200,0.0,0.1);
  TH1D * dphiipH = new TH1D("dphiipH","dphiipH",200,-1,1);
  TH1D * detaipH = new TH1D("detaipH","detaipH",200,-1,1);

  TH1D * nxyip_BSH = new TH1D("nxyip_BSH","nxyipH",200,0.0,0.1);
  TH1D * dphiip_BSH = new TH1D("dphiip_BSH","dphiipH",200,-1,1);
  TH1D * detaip_BSH = new TH1D("detaip_BSH","detaipH",200,-1,1);

  TH1D * ipxEtaH[4];
  TH1D * ipyEtaH[4];
  TH1D * ipzEtaH[4];
  
  // tag and probe 
  int nPtBins = 7;
  float ptBins[8] = {10., 20., 30., 40., 50., 60., 100., 1000.};

  int nEtaBins = 4;
  float etaBins[5] = {0,0.9,1.2,2.1,2.4}; 
  
  TString PtBins[7] = {"Pt10to20",
		       "Pt20to30",
		       "Pt30to40",
		       "Pt40to50",
		       "Pt50to60",
		       "Pt60to100",
		       "PtGt100"};
  
  TString EtaBins[4] = {"EtaLt0p9",
			"Eta0p9to1p2",
			"Eta1p2to2p1",
			"EtaGt2p1"};

  for (int iEta=0; iEta<4; ++iEta) {
    ipxEtaH[iEta] = new TH1D("ipx"+EtaBins[iEta],"",200,-0.02,0.02);
    ipyEtaH[iEta] = new TH1D("ipy"+EtaBins[iEta],"",200,-0.02,0.02);
    ipzEtaH[iEta] = new TH1D("ipz"+EtaBins[iEta],"",200,-0.02,0.02);
  }

  TString JetBins[3] = {"Jet0","Jet1","JetGe2"};

  int nIso1Bins = 4;
  float iso1Bins[5] = {-0.01,0.15,0.2,0.5,1.0}; 
  TString Iso1Bins[4] = {"Iso1Lt0p15",
			 "Iso10p15to0p2",
			 "Iso10p2to0p5",
			 "Iso1Gt0p5"};

  int nIso2Bins = 4;
  float iso2Bins[5] = {-0.01,0.2,0.3,0.5,1.0}; 
  TString Iso2Bins[4] = {"Iso2Lt0p2",
			 "Iso20p2to0p3",
			 "Iso20p3to0p5",
			 "Iso2Gt0p5"};

  //*****  create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //
  TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  etaBinsH->Draw();
  file->cd();
  etaBinsH->Write("etaBinsH");


  //*****  create pt histogram_s with pt ranges associated to their names ***** //
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

  TH1D * ZMassEtaPtPass[4][11];
  TH1D * ZMassEtaPtFail[4][11];

  TH1D * ZMassJetEtaPtPass[3][4][11];
  TH1D * ZMassJetEtaPtFail[3][4][11];

  TH1F * ZMassIso1EtaPtPass[4][4][11];
  TH1F * ZMassIso1EtaPtFail[4][4][11];

  TH1F * ZMassIso2EtaPtPass[4][4][11];
  TH1F * ZMassIso2EtaPtFail[4][4][11];

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassEtaPtPass[iEta][iPt] = new TH1D("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
      ZMassEtaPtFail[iEta][iPt] = new TH1D("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
      for (int iJet=0; iJet<3; ++iJet) {
	ZMassJetEtaPtPass[iJet][iEta][iPt] = new TH1D("ZMass"+JetBins[iJet]+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
	ZMassJetEtaPtFail[iJet][iEta][iPt] = new TH1D("ZMass"+JetBins[iJet]+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
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
  TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

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

  TH1D * metSelNJetsH[3];
  TH1D * puppimetSelNJetsH[3];

  TH1D * metZSelNJetsH[3];
  TH1D * puppimetZSelNJetsH[3];

  TH1D * metTopSelNJetsH[3];
  TH1D * puppimetTopSelNJetsH[3];

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

  TH1D * recoilResponse_nJets[3];
  TH1D * recoilResponse_Ptbins_nJets[3][5];

  TH1D * recoilPuppiResponse_nJets[3];
  TH1D * recoilPuppiResponse_Ptbins_nJets[3][5];

  TH1D * recoilResponseMC_Ptbins[5];
  TH1D * recoilResponseMC_nJets[3];
  TH1D * recoilResponseMC_Ptbins_nJets[3][5];

  TH1D * recoilPuppiResponseMC_Ptbins[5];
  TH1D * recoilPuppiResponseMC_nJets[3];
  TH1D * recoilPuppiResponseMC_Ptbins_nJets[3][5];

  for (int iBin=0; iBin<nJetBins; ++iBin) {

    recoilZParalH[iBin] = new TH1D(RecoilZParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilZPerpH[iBin] = new TH1D(RecoilZPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilPuppiZParalH[iBin] = new TH1D(RecoilPuppiZParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilPuppiZPerpH[iBin] = new TH1D(RecoilPuppiZPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilTopParalH[iBin] = new TH1D(RecoilTopParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilTopPerpH[iBin] = new TH1D(RecoilTopPerp+NJetBins[iBin]+"H","",200,-400,400);

    recoilPuppiTopParalH[iBin] = new TH1D(RecoilPuppiTopParal+NJetBins[iBin]+"H","",200,-400,400);
    recoilPuppiTopPerpH[iBin] = new TH1D(RecoilPuppiTopPerp+NJetBins[iBin]+"H","",200,-400,400);

    metSelNJetsH[iBin] = new TH1D("metSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetSelNJetsH[iBin] = new TH1D("puppimetSel"+NJetBins[iBin]+"H","",200,0.,400.);

    metZSelNJetsH[iBin] = new TH1D("metZSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetZSelNJetsH[iBin] = new TH1D("puppimetZSel"+NJetBins[iBin]+"H","",200,0.,400.);

    metTopSelNJetsH[iBin] = new TH1D("metTopSel"+NJetBins[iBin]+"H","",200,0.,400.);
    puppimetTopSelNJetsH[iBin] = new TH1D("puppimetTopSel"+NJetBins[iBin]+"H","",200,0.,400.);

    recoilResponse_nJets[iBin] = new TH1D("recoilResponse"+NJetBins[iBin],"",400,-20,20);
    recoilPuppiResponse_nJets[iBin] = new TH1D("recoilPuppiResponse"+NJetBins[iBin],"",400,-20,20);

    recoilResponseMC_nJets[iBin] = new TH1D("recoilResponseMC"+NJetBins[iBin],"",400,-20,20);
    recoilPuppiResponseMC_nJets[iBin] = new TH1D("recoilPuppiResponseMC"+NJetBins[iBin],"",400,-20,20);

  }
  for (int iPtBins=0; iPtBins<nZPtBins; ++iPtBins){
    recoilResponseMC_Ptbins[iPtBins] = new TH1D("recoilResponseMC"+ZPtBins[iPtBins],"",400,-20,20);
    recoilPuppiResponseMC_Ptbins[iPtBins] = new TH1D("recoilPuppiResponseMC"+ZPtBins[iPtBins],"",400,-20,20);
  }

  for (int iJets=0; iJets<nJetBins; ++iJets) {
    for (int iPtBins=0; iPtBins<nZPtBins; ++iPtBins){

      recoilZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);
      recoilZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);

      recoilPuppiZParal_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZParal+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);
      recoilPuppiZPerp_Ptbins_nJetsH[iJets][iPtBins] = new TH1D(RecoilPuppiZPerp+NJetBins[iJets]+ZPtBins[iPtBins],"",100,-200,200);

      recoilResponse_Ptbins_nJets[iJets][iPtBins] = new TH1D("recoilResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",400,-20,20);
      recoilPuppiResponse_Ptbins_nJets[iJets][iPtBins] = new TH1D("recoilPuppiResponse"+NJetBins[iJets]+ZPtBins[iPtBins],"",400,-20,20);

      recoilResponseMC_Ptbins_nJets[iJets][iPtBins] = new TH1D("recoilResponseMC"+NJetBins[iJets]+ZPtBins[iPtBins],"",400,-20,20);
      recoilPuppiResponseMC_Ptbins_nJets[iJets][iPtBins] = new TH1D("recoilPuppiResponseMC"+NJetBins[iJets]+ZPtBins[iPtBins],"",400,-20,20);

    }
  }


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

  // HTT Met recoil corrections
  RecoilCorrector recoilPFMetCorrector(RecoilFileName);

  // Lepton Scale Factors 

  ScaleFactor * SF_muonIdIso; 
  ScaleFactor * SF_muonTrig;
  if (applyLeptonSF) {
    SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
  }

  // tracking efficiency SF
  TFile * workspaceFile = new TFile(TString(cmsswBase)+"/src/"+CorrectionWorkspaceFile);
  RooWorkspace *correctionWS = (RooWorkspace*)workspaceFile->Get("w");

  // IP correction
  IpCorrection * ipCorrection = new IpCorrection(TString(cmsswBase)+"/src/"+IpCorrFileName);

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

  // BTag scale factors                                                                                                                                                    
  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv");
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  reader_B.load(calib,BTagEntry::FLAV_B,"comb");
  reader_C.load(calib,BTagEntry::FLAV_C,"comb");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");

  float etaBTAG[2] = {0.5,2.1};
  float ptBTAG[5] = {25.,35.,50.,100.,200.};

  for (int iEta=0; iEta<2; ++iEta) {
    for (int iPt=0; iPt<5; ++iPt) {
      float sfB = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
      
      float sfB_Up = reader_B.eval_auto_bounds("up",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC_Up = reader_C.eval_auto_bounds("up",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight_Up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);

      float sfB_Down = reader_B.eval_auto_bounds("down",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC_Down = reader_C.eval_auto_bounds("down",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight_Down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
      
      printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f + %5.3f - %5.3f ; SFc = %5.3f + %5.3f - %5.3f ; SFl = %5.3f + %5.3f - %5.3f\n",
	     ptBTAG[iPt],etaBTAG[iEta],
	     sfB,sfB_Up,sfB_Down,
	     sfC,sfC_Up,sfC_Down,
	     sfLight,sfLight_Up,sfLight_Down);
      
        }
  }

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_ichep2016.root"));
  TH2F * tagEff_B = (TH2F*)fileTagging->Get("btag_eff_b");
  TH2F * tagEff_C = (TH2F*)fileTagging->Get("btag_eff_c");
  TH2F * tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");

  float MaxBJetPt = 1000;
  float MinBJetPt = 20;

  TRandom3 rand;

  int nFiles = 0;
  int nEvents = 0;
  int selEventsAllMuons = 0;
  int selEventsIdMuons = 0;
  int selEventsIsoMuons = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();
 
  for (int iF=0; iF<nTotalFiles; ++iF) {
  
    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));

    TTree * _inittree = NULL; 
    _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree != NULL) {
      Float_t genweight;
      if (!isData)
	_inittree->SetBranchAddress("genweight",&genweight);
      Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
      std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
      for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
	_inittree->GetEntry(iEntry);
	if (isData)
	  histWeightsH->Fill(0.,1.);
	else {
	  Float_t GenWeight = 1;
	  if (genweight<0) GenWeight = -1;
	  histWeightsH->Fill(0.,GenWeight);
	}
      }
    }

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
    if (_tree==NULL) continue;
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

      //------------------------------------------------

      if (!isData) { 
	float genweight = 1;
	if (analysisTree.genweight<0)
	  genweight = -1;
	weight *= genweight;
      }
      
      TLorentzVector genZ; genZ.SetXYZT(0,0,0,0); 
      TLorentzVector genV; genV.SetXYZT(0,0,0,0);
      TLorentzVector genL; genL.SetXYZT(0,0,0,0);
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
	      dyWeight = histZMassPtWeights->GetBinContent(histZMassPtWeights->FindBin(genZMass,genZPt));

	    //	    std::cout << "Z mass = " << genZMass << "   Z Pt = " << genZPt << "   weight = " << dyWeight << std::endl;
	    weight *= dyWeight;
	  }

	}
	
	//	cout << "PUinteractions = " << analysisTree.numtruepileupinteractions << endl;

        if (applyPUreweighting_official) {
	  nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
	  double Ninteractions = analysisTree.numtruepileupinteractions;
	  double PUweight = PUofficial->get_PUweight(Ninteractions);
	  weight *= float(PUweight);
	  PUweightsOfficialH->Fill(PUweight);
	  //	  cout << PUweight << endl;
        }

	if (applyTauTauSelection) {
	  unsigned int nTaus = 0;
	  if (analysisTree.gentau_count>0) {
	    //	  cout << "Generated taus present" << endl;
	    for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++itau) {
	      //	    cout << itau << "  : pt = " 
	      //		 << analysisTree.gentau_visible_pt[itau] 
	      //		 << "   eta = " <<  analysisTree.gentau_visible_eta[itau]
	      //		 << "   mother = " << int(analysisTree.gentau_mother[itau]) << endl;
	      if (int(analysisTree.gentau_mother[itau])==3) nTaus++;
	      
	    }
	  }
	  bool notTauTau = nTaus < 2;
	  //	  std::cout << "nTaus = " << nTaus << std::endl;
	  
	  if (selectZToTauTauMuMu&&notTauTau) { 
	    //	    std::cout << "Skipping event..." << std::endl;
	    //	    cout << endl;
	    continue;
	  }
	  if (!selectZToTauTauMuMu&&!notTauTau) { 
	    //	    std::cout << "Skipping event..." << std::endl;
	    //	    cout << endl;
	    continue;
	  }
	  //	  cout << endl;
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
	    //	    cout << "toppt = " << topPt << "   antitoppt = " << antitopPt << "   weight = " << topptweight << endl;
	    weight *= topptweight;
	  }
	}
      }

      
      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;
      
      //      continue;
      //      std::cout << "Event " << analysisTree.event_nr << std::endl;

      if (isData && applyGoodRunSelection){
	
	
	bool lumi = false;
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
	//std::remove("myinputfile");
      }     

      //      std::cout << "Event passed good run selection " << std::endl;

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
      
      bool isTriggerMuon = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(MuonTriggerName)) {
	  //	  std::cout << it->first << " : " << it->second << std::endl;
	  if (it->second==1)
	    isTriggerMuon = true;
	}
      }

      if (applyTrigger && !isTriggerMuon) continue;
      //      std::cout << "Event passed trigger" << std::endl;


      unsigned int nMuonFilter = 0;
      bool isMuonFilter = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	//	std::cout << HLTFilter << std::endl;
	if (HLTFilter==MuonFilterName) {
	  nMuonFilter = i;
	  isMuonFilter = true;
	}
      }

      if (applyTrigger) {
	if (!isMuonFilter) {
	  cout << "Filter " << MuonFilterName << " not found " << endl;
	  continue;
	}
      }

      float pfmet_ex = analysisTree.pfmetcorr_ex;
      float pfmet_ey = analysisTree.pfmetcorr_ey;
      float pfmet_phi = analysisTree.pfmetcorr_phi;
      float pfmet = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);

      float puppimet_ex = analysisTree.puppimet_ex;
      float puppimet_ey = analysisTree.puppimet_ey;
      float puppimet_phi = analysisTree.puppimet_phi;
      float puppimet = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);
      
      // muon selection

      // probe muons
      vector<unsigned int> allMuons; allMuons.clear();
      vector<float> allMuonsIso; allMuonsIso.clear();
      vector<bool> allMuonsPassedId; allMuonsPassedId.clear();
      vector<bool> allMuonsPassedIdIso; allMuonsPassedIdIso.clear();

      // tag muons
      vector<unsigned int> isoIdMuons; isoIdMuons.clear();
      vector<bool> isoIdMuonsIsTriggerMatched; isoIdMuonsIsTriggerMatched.clear();
      vector<float> isoIdMuonsIso; isoIdMuonsIso.clear();

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (!isData&&fabs(muonMomScale-1.0)>0.001) {
	  analysisTree.muon_px[im] *= muonMomScale;
	  analysisTree.muon_py[im] *= muonMomScale;
	  analysisTree.muon_pz[im] *= muonMomScale;
	  analysisTree.muon_pt[im] *= muonMomScale;
	} 

	if (analysisTree.muon_pt[im]<ptMuonProbeCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonProbeCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonProbeCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonProbeCut) continue;

	bool isMedium = analysisTree.muon_isMedium[im];
	bool isPassed = true;
        if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) isPassed = false;
        if (fabs(analysisTree.muon_dz[im])>dzMuonCut)  isPassed = false;
	float absIso = 0;
	if (isoDR03) { 
	  absIso = analysisTree.muon_r03_sumChargedHadronPt[im];
	  float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] + 
	    analysisTree.muon_r03_sumPhotonEt[im] - 
	    0.5*analysisTree.muon_r03_sumPUPt[im];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  absIso += neutralIso;
	}
	else {
	  absIso = analysisTree.muon_chargedHadIso[im];
          float neutralIso = analysisTree.muon_neutralHadIso[im] +
            analysisTree.muon_photonIso[im] -
            0.5*analysisTree.muon_puIso[im];
          neutralIso = TMath::Max(float(0),neutralIso);
          absIso += neutralIso;
	}
	float relIso = absIso/analysisTree.muon_pt[im];
	bool isPassedId = isMedium && isPassed;
	bool isPassedIdIso = isMedium && isPassed && relIso<isoMuonCut; 
	allMuons.push_back(im);
        allMuonsPassedId.push_back(isPassedId);
        allMuonsPassedIdIso.push_back(isPassedIdIso);
        allMuonsIso.push_back(relIso);

	bool isTriggerMatched = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nMuonFilter])
	    isTriggerMatched = true;
	}
	
	if (isPassedIdIso) {
	  isoIdMuons.push_back(im);
	  isoIdMuonsIso.push_back(relIso);
	  isoIdMuonsIsTriggerMatched.push_back(isTriggerMatched);
	}

      }

      unsigned int indx1 = 0;
      unsigned int indx2 = 0;
      bool isIsoMuonsPair = false;
      bool firstTrigger = true;
      float isoMin = 9999;
      if (isoIdMuons.size()>0) {
	for (unsigned int im1=0; im1<isoIdMuons.size(); ++im1) {
	  unsigned int index1 = isoIdMuons[im1];
	  bool isTriggerMatched = isoIdMuonsIsTriggerMatched[im1];
	  if (isTriggerMatched &&  analysisTree.muon_pt[index1] > ptMuonCut && fabs(analysisTree.muon_eta[index1]) < etaMuonCut) {
	    for (unsigned int iMu=0; iMu<allMuons.size(); ++iMu) {
	      unsigned int indexProbe = allMuons[iMu];
	      if (index1==indexProbe) continue;
	      float q1 = analysisTree.muon_charge[index1];
	      float q2 = analysisTree.muon_charge[indexProbe];
	      if (q1*q2>0) continue;
	      float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				analysisTree.muon_eta[indexProbe],analysisTree.muon_phi[indexProbe]);
	      if (dR<dRleptonsCut) continue; 
	      float dPhi = dPhiFrom2P(analysisTree.muon_px[index1],analysisTree.muon_py[index1],
				      analysisTree.muon_px[indexProbe],analysisTree.muon_py[indexProbe]);
	      if (dPhi>dPhileptonsCut) continue;
	      float ptProbe = TMath::Min(float(analysisTree.muon_pt[indexProbe]),float(ptBins[nPtBins]-0.1));
	      float absEtaProbe = fabs(analysisTree.muon_eta[indexProbe]);
	      float isoProbe = TMath::Min(allMuonsIso[iMu],float(iso1Bins[nIso1Bins]-0.001));
	      int ptBin = binNumber(ptProbe,nPtBins,ptBins);
	      int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
	      int iso1Bin = binNumber(isoProbe,nIso1Bins,iso1Bins);
              int iso2Bin = binNumber(isoProbe,nIso2Bins,iso2Bins);

	      bool passId = allMuonsPassedId[iMu];

	      if (ptBin<0) continue;
	      if (etaBin<0) continue;
	      if (iso1Bin<0) continue;
	      if (iso2Bin<0) continue;
	      
	      TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
						  analysisTree.muon_py[index1],
						  analysisTree.muon_pz[index1],
						  muonMass);
	      TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[indexProbe],
						  analysisTree.muon_py[indexProbe],
						  analysisTree.muon_pz[indexProbe],
						  muonMass);

	      // number of jets
	      int nJets30 = 0;
	      int nJets30etaCut = 0;

	      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
		float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
		if (absJetEta>jetEtaCut) continue;
	  
		float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.muon_eta[index1],analysisTree.muon_phi[index1]);
		if (dR1<dRJetLeptonCut) continue;
		
		float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
				   analysisTree.muon_eta[indexProbe],analysisTree.muon_phi[indexProbe]);
		
		if (dR2<dRJetLeptonCut) continue;
	  
		// pfJetId
		bool isPFJetId = tightJetiD_2017(analysisTree,int(jet));
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
 
	      float mass = (muon1+muon2).M();
	      if (allMuonsPassedIdIso[iMu]) { 
		ZMassPass->Fill(mass,weight);
		ZMassEtaPtPass[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtPass[JetBin][etaBin][ptBin]->Fill(mass,weight);
	      }
	      else {
		ZMassFail->Fill(mass,weight);
		ZMassEtaPtFail[etaBin][ptBin]->Fill(mass,weight);
		ZMassJetEtaPtFail[JetBin][etaBin][ptBin]->Fill(mass,weight);
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
	  for (unsigned int im2=im1+1; im2<isoIdMuons.size(); ++im2) {
	    unsigned int index2 = isoIdMuons[im2];
	    float q1 = analysisTree.muon_charge[index1];
	    float q2 = analysisTree.muon_charge[index2];
	    bool isTriggerMatched1 = isoIdMuonsIsTriggerMatched[im1] &&
              analysisTree.muon_pt[index1]>ptMuonCut&&
	      fabs(analysisTree.electron_eta[index1])<etaMuonCut;
            bool isTriggerMatched2 = isoIdMuonsIsTriggerMatched[im2]&&
              analysisTree.muon_pt[index2]>ptMuonCut &&
              fabs(analysisTree.muon_eta[index2])<etaMuonCut;
            bool isTriggerMatched = isTriggerMatched1 || isTriggerMatched2;
	    bool isPairSelected = q1*q2 > 0;
	    if (oppositeSign) isPairSelected = q1*q2 < 0;
	    float dRmumu = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
				  analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
	    if (isTriggerMatched && isPairSelected && dRmumu>dRleptonsCut) {
	      bool sumIso = isoIdMuonsIso[im1]+isoIdMuonsIso[im2];
	      if (sumIso<isoMin) {
		isIsoMuonsPair = true;
		isoMin = sumIso;
		if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]) {
		  indx1 = index1;
		  indx2 = index2;
		}
		else {
		  indx2 = index1;
		  indx1 = index2;
		}
	      }
	    }
	  }
	}
      }
      //      std::cout << "isolated pair " << isIsoMuonsPair << std::endl;


      if (isIsoMuonsPair) {      
	//match to genparticles
	bool genmatch_m1 = false, genmatch_m2 = false;

	TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[indx1],
					analysisTree.muon_py[indx1],
					analysisTree.muon_pz[indx1],
					muonMass);

	TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[indx2],
					analysisTree.muon_py[indx2],
					analysisTree.muon_pz[indx2],
					muonMass);

	TLorentzVector dimuon = mu1 + mu2;

	float massSel = dimuon.M();

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
	
	// number of btag jets
	int nJets20BTagLoose = 0;
	int nJets20BTagMedium = 0;

	int indexLeadingJet = -1;
	float ptLeadingJet = -1;

	int indexSubLeadingJet = -1;
	float ptSubLeadingJet = -1;

	for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
	  //	  std::cout << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
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

	  float jetEta = analysisTree.pfjet_eta[jet];
	  float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	  float jetPt = analysisTree.pfjet_pt[jet];
	  if (absJetEta>jetEtaCut) continue;
	  
	  float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			     analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1]);
	  if (dR1<dRJetLeptonCut) continue;
	  
	  float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			     analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	  if (dR2<dRJetLeptonCut) continue;
	  
	  // pfJetId
	  bool isPFJetId = tightJetiD_2017(analysisTree,int(jet));
	  if (!isPFJetId) continue;

	  bool noisyJet = analysisTree.pfjet_pt[jet]<50 && absJetEta > 2.65 && absJetEta < 3.139;

	  if (noisyJet) continue;
	  
	  if (analysisTree.pfjet_pt[jet]>jetPtHighCut) {
	    nJets30++;
	    HT30 += analysisTree.pfjet_pt[jet];

	    if (indexLeadingJet>=0) {
	      if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
		indexSubLeadingJet = jet;
		ptSubLeadingJet = jetPt;
	      }
	    }

	    if (jetPt>ptLeadingJet) {
	      indexSubLeadingJet = indexLeadingJet;
	      ptSubLeadingJet = ptLeadingJet;
	      indexLeadingJet = jet;
	      ptLeadingJet = jetPt; 
	    }

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

	      bool tagged = analysisTree.pfjet_btag[jet][0]>0.80;

	      float JetPtForBTag = analysisTree.pfjet_pt[jet];
	      float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      float jet_scalefactor = 1;
	      float tageff = 1;
	      if (!isData) {
		int flavor = abs(analysisTree.pfjet_flavour[jet]);
		if (flavor==5) { // b-quark
		  jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
		  tageff = tagEff_B->GetBinContent(tagEff_B->FindBin(JetPtForBTag,absJetEta));
		}
		else if (flavor==4) { // c-quark
		  jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                  tageff = tagEff_C->GetBinContent(tagEff_C->FindBin(JetPtForBTag,absJetEta));
		}
		else {
		  jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
		  tageff = tagEff_Light->GetBinContent(tagEff_Light->FindBin(JetPtForBTag,absJetEta));
		}

		//		std::cout << "Flavor = " << flavor << " pt = " <<  JetPtForBTag << "  eta = " << jetEta 
		//			  << "  SF = " << jet_scalefactor << " tageff = " << tageff << "  tagged = " << tagged << std::endl;
		//
		if (tageff<1e-5)      tageff = 1e-5;
		if (tageff>0.99999)   tageff = 0.99999;
		rand.SetSeed((int)((jetEta+5)*100000));
		double rannum = rand.Rndm();
		if (jet_scalefactor<1 && tagged) { // downgrading
		  double fraction = 1-jet_scalefactor;
		  if (rannum<fraction) {
		    tagged = false;
		  }
		}
		if (jet_scalefactor>1 && !tagged) { // upgrading
		  double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
		  if (rannum<fraction) {
		    tagged = true;
		  }
		}

		//		std::cout << "Updated taggged = " << tagged << std::endl;
	      }


	      if (analysisTree.pfjet_btag[jet][0]>0.46) 
		nJets20BTagLoose++;
	      if (tagged) 
		nJets20BTagMedium++;

            }
	  }	  


	}

	float mjj = -1;
	float etaLeadingJet = -999;
	float phiLeadingJet = -999;

	if (indexLeadingJet>=0) {
	  ptLeadingJet = analysisTree.pfjet_pt[indexLeadingJet];
	  etaLeadingJet = analysisTree.pfjet_eta[indexLeadingJet];
	  phiLeadingJet = analysisTree.pfjet_phi[indexLeadingJet];
	}

	if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {


	  TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
					       analysisTree.pfjet_py[indexLeadingJet],
					       analysisTree.pfjet_pz[indexLeadingJet],
					       analysisTree.pfjet_e[indexLeadingJet]);
	  
	  TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
                                               analysisTree.pfjet_py[indexSubLeadingJet],
                                               analysisTree.pfjet_pz[indexSubLeadingJet],
                                               analysisTree.pfjet_e[indexSubLeadingJet]);
	  mjj = (jet1+jet2).M();

	}

	bool is0Jet    = nJets30==0;
	bool isBoosted = nJets30==1 || (nJets30==2 && mjj<300) || nJets30>2;
	bool isVBF     = nJets30==2 && mjj>300;

	//	std::cout << "Jet processed" << std::endl;

	bool isRefittedVtx = false;
	TVector3 vertex = get_refitted_PV_with_BS(&analysisTree,indx1,indx2,isRefittedVtx);
	TVector3 ip1    = ipVec_Lepton(&analysisTree,indx1,vertex);
	TVector3 ip2    = ipVec_Lepton(&analysisTree,indx2,vertex);
	double ptMu1 = (double)analysisTree.muon_pt[indx1];
	double ptMu2 = (double)analysisTree.muon_pt[indx2];
	double etaMu1 = (double)analysisTree.muon_eta[indx1];
	double etaMu2 = (double)analysisTree.muon_eta[indx2];

	if (!isData && applyLeptonSF) {

	  //leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
	  double IdIsoSF_mu1 = SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
	  double IdIsoSF_mu2 = SF_muonIdIso->get_ScaleFactor(ptMu2, etaMu2);


	  /*
	  std::cout << "isRefittedVtx = " << isRefittedVtx << std::endl;
	  std::cout << "IP1 (x,y,z) = (" << ip1.X() << "," << ip1.Y() << "," << ip1.Z() << ")" << std::endl; 
	  std::cout << "IP2 (x,y,z) = (" << ip2.X() << "," << ip2.Y() << "," << ip2.Z() << ")" << std::endl; 
	  std::cout << std::endl;
	  */
	  correctionWS->var("m_eta")->setVal(etaMu1);
	  correctionWS->var("m_pt")->setVal(ptMu1);
	  double trackSF_mu1 = correctionWS->function("m_trk_ratio")->getVal();
	  if (applyKITCorrection)
	    IdIsoSF_mu1 = correctionWS->function("m_id_ratio")->getVal()*correctionWS->function("m_iso_ratio")->getVal();


	  correctionWS->var("m_eta")->setVal(etaMu2);
          correctionWS->var("m_pt")->setVal(ptMu2);
	  double trackSF_mu2 = correctionWS->function("m_trk_ratio")->getVal();
	  if (applyKITCorrection)
	    IdIsoSF_mu2 = correctionWS->function("m_id_ratio")->getVal()*correctionWS->function("m_iso_ratio")->getVal();

	  MuSF_IdIso_Mu1H->Fill(IdIsoSF_mu1);
	  MuSF_IdIso_Mu2H->Fill(IdIsoSF_mu2);
	  /*	  if (ptMu1<20||ptMu2<20) {
	  	  std::cout << "mu 1 ->  pt = " << ptMu1 << "   eta = " << etaMu1 << std::endl;
	  	  std::cout << "eff data mu 1 = " << SF_muonIdIso->get_EfficiencyData(ptMu1, etaMu1)<< " |  eff mc mu 1 = " << SF_muonIdIso->get_EfficiencyMC(ptMu1, etaMu1)<<std::endl;
	  	  std::cout << "mu 2 ->  pt = " << ptMu2 << "   eta = " << etaMu2 << std::endl;
		  std::cout << "eff data mu 2 = " << SF_muonIdIso->get_EfficiencyData(ptMu2, etaMu2)<< " |  eff mc mu 2 = " << SF_muonIdIso->get_EfficiencyMC(ptMu2, etaMu2)<<std::endl;
		  std::cout << "SF mu1 = " << IdIsoSF_mu1 << std::endl;
		  std::cout << "SF mu2 = " << IdIsoSF_mu2 << std::endl;
		  
		  std::cout << " mass = " << massSel << std::endl;
		  std::cout << std::endl;
		  } 
	  */
	  weight = weight*IdIsoSF_mu1*IdIsoSF_mu2*trackSF_mu1*trackSF_mu2;

	  double effDataTrig1 = SF_muonTrig->get_EfficiencyData(ptMu1, etaMu1);  
	  double effDataTrig2 = SF_muonTrig->get_EfficiencyData(ptMu2, etaMu2);  
	  double effTrigData = 1.0 - (1-effDataTrig1)*(1-effDataTrig2);

	  double weightTrig = 1;
	  if (applyTrigger) {
	    double effMcTrig1 = SF_muonTrig->get_EfficiencyMC(ptMu1, etaMu1);
	    double effMcTrig2 = SF_muonTrig->get_EfficiencyMC(ptMu2, etaMu2);
	    double effMcTrig = 1 - (1-effMcTrig1)*(1-effMcTrig2);
	    if (effTrigData>0&&effMcTrig>0) {
	      weightTrig = effTrigData/effMcTrig;
	      weight = weight*weightTrig;
	    }
	  }
	  else {
	    weightTrig = effTrigData;
	    weight = weight*weightTrig;
	  }
	}

	float visiblePx = dimuon.Px();
	float visiblePy = dimuon.Py();
	if (applyRecoilOnGenerator) {
	  visiblePx = genL.Px();
	  visiblePy = genL.Py();
	}

	if (!isData && applyRecoilCorrections) {

	  float pfmetcorr_ex = pfmet_ex;
	  float pfmetcorr_ey = pfmet_ey;

	  if (applySimpleRecoilCorrections)
	    recoilPFMetCorrector.CorrectByMeanResolution(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);
	  else 
	    recoilPFMetCorrector.Correct(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),visiblePx,visiblePy,nJets30,pfmetcorr_ex,pfmetcorr_ey);
	  //	  std::cout << "PFMet : (" << pfmet_ex << "," << pfmet_ey << ")  "
	  //		    << "  Quantile : (" << pfmetcorr_ex << "," << pfmetcorr_ey << ")" 
	  //		    << "  MeanReso : (" << pfmetcorr_ex_1 << "," << pfmetcorr_ey_1 << ")" 
	  //		    << std::endl; 
 	  pfmet_phi = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
          pfmet = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
	  pfmet_ex = pfmetcorr_ex;
	  pfmet_ey = pfmetcorr_ey;
	  
	}

	// selection on mass
	if (massSel>dimuonMassCut) {

	  numberOfVerticesH->Fill(double(analysisTree.primvertex_count),weight);

	  massSelH->Fill(massSel,weight);
	  massExtendedSelH->Fill(massSel,weight);
	  for (int iScale=0; iScale<21; ++ iScale) {
	    float scaleFactor = 0.98 + 0.002*float(iScale);
	    massSelScaleH[iScale]->Fill(massSel*scaleFactor,weight);
	  }

	  if (is0Jet) {
	    mass0jetH->Fill(double(massSel),weight);
	    if (nJets20BTagMedium==0) 
	      mass0jetBvetoH->Fill(double(massSel),weight);
	  }

	  if (isBoosted) {
	    massBoostedH->Fill(double(massSel),weight);
	    if (nJets20BTagMedium==0)
	      massBoostedBvetoH->Fill(double(massSel),weight);
	  }

	  if (isVBF) {
	    massVBFH->Fill(double(massSel),weight);
	    if (nJets20BTagMedium==0)
	      massVBFBvetoH->Fill(double(massSel),weight);
	  }

	  ptLeadingMuSelH->Fill(analysisTree.muon_pt[indx1],weight);
	  ptTrailingMuSelH->Fill(analysisTree.muon_pt[indx2],weight);
	  etaLeadingMuSelH->Fill(analysisTree.muon_eta[indx1],weight);
	  etaTrailingMuSelH->Fill(analysisTree.muon_eta[indx2],weight);
	  
	  nJets30SelH->Fill(double(nJets30),weight);
	  nJets20SelH->Fill(double(nJets20),weight);
	  nJets30etaCutSelH->Fill(double(nJets30etaCut),weight);
	  nJets20etaCutSelH->Fill(double(nJets20etaCut),weight);

	  HT30SelH->Fill(double(HT30),weight);
	  HT20SelH->Fill(double(HT20),weight);
	  HT30etaCutSelH->Fill(double(HT30etaCut),weight);
	  HT20etaCutSelH->Fill(double(HT20etaCut),weight);

	  nJets20BTagLooseSelH->Fill(double(nJets20BTagLoose),weight);
	  nJets20BTagMediumSelH->Fill(double(nJets20BTagMedium),weight);

	  // cout << "dxy (mu1) = " << analysisTree.muon_dxy[indx1] << "   error = " << analysisTree.muon_dxyerr[indx1] << std::endl;
	  //cout << "dxy (mu2) = " << analysisTree.muon_dxy[indx2] << "   error = " << analysisTree.muon_dxyerr[indx2] << std::endl;
	  //	 vector<int>indexMuTau; indexMuTau.clear();  
	  //	 if (selectZToTauTauMuMu && indexMuTau.size()!=2) continue;
	  float dimuonEta = dimuon.Eta();
	  float dimuonPt  = dimuon.Pt();
	  float dimuonPhi = dimuon.Phi();
	  float sumMuonPt = (mu1.Pt()+mu2.Pt());
	  float ptRatio = 0.0;
	  if (sumMuonPt != 0)
	    ptRatio = (dimuonPt/sumMuonPt);

	  dimuonMassPtH->Fill(massSel,dimuonPt,weight);
	  dimuonPtSelH->Fill(dimuon.Pt(),weight);
	  dimuonEtaSelH->Fill(dimuon.Eta(),weight);

	  NumberOfVerticesH->Fill(float(analysisTree.primvertex_count),weight);
	  metSelH->Fill(pfmet,weight);
	  puppimetSelH->Fill(puppimet,weight);

	  float unitX = dimuon.Px()/dimuon.Pt();
	  float unitY = dimuon.Py()/dimuon.Pt();
	  float phiUnit = TMath::ATan2(unitY,unitX);
	  float perpUnitX = TMath::Cos(phiUnit+0.5*TMath::Pi());
	  float perpUnitY = TMath::Sin(phiUnit+0.5*TMath::Pi());
	  
	  for (int iScale=0; iScale<21; ++ iScale) {
	    float scaleFactor = 0.7 + 0.03*float(iScale);
	    metSelScaleH[iScale]->Fill(pfmet*scaleFactor,weight);
	    puppimetSelScaleH[iScale]->Fill(puppimet*scaleFactor,weight);
	  }

	  // jet pt bin
	  int jetBin = 0;
	  if (nJets30==1)
	    jetBin = 1;
	  else if (nJets30>1)
	    jetBin = 2;

	  int ptBin = binNumber(TMath::Min(float(dimuonPt),float(999)),nZPtBins,zPtBins);

	  metSelNJetsH[jetBin]->Fill(pfmet,weight);
	  puppimetSelNJetsH[jetBin]->Fill(puppimet,weight);

	  if (!isData) {

	    int ptBinV = binNumber(TMath::Min(float(genV.Pt()),float(999)),nZPtBins,zPtBins);
	    float Hparal = 0;
	    float Hperp  = 0;

	    // pfmet ->
	    ComputeHadRecoilFromMet(pfmet_ex,pfmet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	    float response = - Hparal/genV.Pt();
	    recoilResponseMC_Ptbins[ptBinV]->Fill(response,weight);
	    recoilResponseMC_nJets[jetBin]->Fill(response,weight);
	    recoilResponseMC_Ptbins_nJets[jetBin][ptBinV]->Fill(response,weight);
	    //	    std::cout << "PFMet : Had(paral) = " << Hparal << "  Had(perp) = " << Hperp << std::endl;

	    // puppi met ->
	    ComputeHadRecoilFromMet(puppimet_ex,puppimet_ey,genV.Px(),genV.Py(),genL.Px(),genL.Py(),Hparal,Hperp);
	    response = - Hparal/genV.Pt();
	    recoilPuppiResponseMC_Ptbins[ptBinV]->Fill(response,weight);
	    recoilPuppiResponseMC_nJets[jetBin]->Fill(response,weight);
	    recoilPuppiResponseMC_Ptbins_nJets[jetBin][ptBinV]->Fill(response,weight);
	    //	    std::cout << "PuppiMet : Had(paral) = " << Hparal << "  Had(perp) = " << Hperp << std::endl;

	  }

	  float recoilParal = 0;
	  float recoilPerp  = 0;
	  float responseHad = 0;
	  computeRecoil(pfmet_ex,pfmet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilParal,recoilPerp,responseHad);

	  float recoilPuppiParal = 0;
	  float recoilPuppiPerp  = 0;
	  float responsePuppiHad = 0;
	  computeRecoil(puppimet_ex,puppimet_ey,unitX,unitY,perpUnitX,perpUnitY,dimuonPt,recoilPuppiParal,recoilPuppiPerp,responsePuppiHad);

	  if (massSel>70&&massSel<110) { // Z Region

	    int etaBin1 = binNumber(TMath::Min(float(fabs(etaMu1)),float(2.4)),nEtaBins,etaBins);
	    int etaBin2 = binNumber(TMath::Min(float(fabs(etaMu2)),float(2.4)),nEtaBins,etaBins);
	    
	    //	    cout << "eta(Mu1) = " << etaMu1 << " : " << etaBin1 << std::endl;
	    //	    cout << "eta(Mu2) = " << etaMu2 << " : " << etaBin2 << std::endl;

	    double ipx1 = ip1.X();
	    double ipx2 = ip2.X();

	    double ipy1 = ip1.Y();
	    double ipy2 = ip2.Y();

	    double ipz1 = ip1.Z();
	    double ipz2 = ip2.Z();
	    
	    if (applyIpCorrection) {

	      ipx1 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipx,ip1.X(),mu1.Eta());
	      ipy1 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipy,ip1.Y(),mu1.Eta());
	      ipz1 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipz,ip1.Z(),mu1.Eta());

	      ipx2 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipx,ip2.X(),mu2.Eta());
	      ipy2 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipy,ip2.Y(),mu2.Eta());
	      ipz2 = ipCorrection->correctIp(IpCorrection::Coordinate::Ipz,ip2.Z(),mu2.Eta());

	      /*
	      std::cout << "IPx (1) = " << ipx1 << " : " << ip1.X() << std::endl;
	      std::cout << "IPy (1) = " << ipy1 << " : " << ip1.Y() << std::endl;
	      std::cout << "IPz (1) = " << ipz1 << " : " << ip1.Z() << std::endl;

	      std::cout << "IPx (2) = " << ipx2 << " : " << ip2.X() << std::endl;
	      std::cout << "IPy (2) = " << ipy2 << " : " << ip2.Y() << std::endl;
	      std::cout << "IPz (2) = " << ipz2 << " : " << ip2.Z() << std::endl;
	      */
	    }

	    ipxH->Fill(ipx1,weight);
	    ipxH->Fill(ipx2,weight);

	    ipyH->Fill(ipy1,weight);
	    ipyH->Fill(ipy2,weight);

	    ipzH->Fill(ipz1,weight);
	    ipzH->Fill(ipz2,weight);

	    TLorentzVector n1; n1.SetXYZM(ipx1,ipy1,ipz1,0);
	    TLorentzVector n2; n2.SetXYZM(ipx2,ipy2,ipz2,0);

	    nxyipH->Fill(n1.Pt(),weight);
	    nxyipH->Fill(n2.Pt(),weight);
	    if (isRefittedVtx) {
	      nxyip_BSH->Fill(n1.Pt(),weight);
	      nxyip_BSH->Fill(n2.Pt(),weight);
	    }

	    double cosdphi1 = (mu1.Px()*n1.Px()+mu1.Py()*n1.Py())/(mu1.Pt()*n1.Pt());
	    double dphi1 = TMath::ACos(cosdphi1);
	    double sign = mu1.Px()*n1.Py()-mu1.Py()*n1.Px();
	    if (sign<0) 
	      dphi1 = -dphi1;
	    dphiipH->Fill(dphi1,weight);
	    double deta1 = n1.Eta()-mu1.Eta();
	    detaipH->Fill(deta1,weight);

	    double cosdphi2 = (mu2.Px()*n2.Px()+mu2.Py()*n2.Py())/(mu2.Pt()*n2.Pt());
	    double dphi2 = TMath::ACos(cosdphi2);
	    sign = mu2.Px()*n2.Py()-mu2.Py()*n2.Px();
	    if (sign<0) 
	      dphi2 = -dphi2;
	    dphiipH->Fill(dphi2,weight);
	    double deta2 = n2.Eta()-mu2.Eta();
	    detaipH->Fill(deta2,weight);

	    if (isRefittedVtx) {	      
	      dphiip_BSH->Fill(dphi1,weight);
	      detaip_BSH->Fill(deta1,weight);
	      dphiip_BSH->Fill(dphi2,weight);
	      detaip_BSH->Fill(deta2,weight);
	    }

	    ipxEtaH[etaBin1]->Fill(ipx1,weight);
	    ipxEtaH[etaBin2]->Fill(ipx2,weight);

	    ipyEtaH[etaBin1]->Fill(ipy1,weight);
	    ipyEtaH[etaBin2]->Fill(ipy2,weight);

	    ipzEtaH[etaBin1]->Fill(ipz1,weight);
	    ipzEtaH[etaBin2]->Fill(ipz2,weight);

	    ipdxH->Fill(ipx1-ipx2,weight);
	    ipdyH->Fill(ipy1-ipy2,weight);
	    ipdzH->Fill(ipz1-ipz2,weight);

	    metZSelH->Fill(pfmet,weight);
	    puppimetZSelH->Fill(puppimet,weight);

	    metZSelNJetsH[jetBin]->Fill(pfmet,weight);
	    puppimetZSelNJetsH[jetBin]->Fill(puppimet,weight);

	    nJets20BTagLooseZSelH->Fill(double(nJets20BTagLoose),weight);
	    nJets20BTagMediumZSelH->Fill(double(nJets20BTagMedium),weight);

	    nJets30ZSelH->Fill(double(nJets30),weight);
	    nJets20ZSelH->Fill(double(nJets20),weight);
	    nJets30etaCutZSelH->Fill(double(nJets30etaCut),weight);
	    nJets20etaCutZSelH->Fill(double(nJets20etaCut),weight);

	    if (nJets30>0) {
	      leadingJetPtH->Fill(ptLeadingJet,weight);
	      leadingJetEtaH->Fill(etaLeadingJet,weight);
	      leadingJetPhiH->Fill(phiLeadingJet,weight);
	    }
	    
	    if (is0Jet) {
	      nJets20BTagMedium0jetH->Fill(double(nJets20BTagMedium),weight);
	    }
	    if (isBoosted) {
	      nJets20BTagMediumBoostedH->Fill(double(nJets20BTagMedium),weight);
	      dimuonPtBoostedH->Fill(dimuonPt,weight);
	      if (nJets20BTagMedium==0) 
		dimuonPtBoostedBvetoH->Fill(dimuonPt,weight);
	    }
	    if (isVBF) {
	      nJets20BTagMediumVBFH->Fill(double(nJets20BTagMedium),weight);
	      mjjVBFH->Fill(mjj,weight);
	      if ( nJets20BTagMedium==0 )
		mjjVBFBvetoH->Fill(mjj,weight);
	    }

	    // pfmet
	    recoilZParalH[jetBin]->Fill(recoilParal,weight);
	    recoilZPerpH[jetBin]->Fill(recoilPerp,weight);
	    recoilZParal_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilParal,weight);
	    recoilZPerp_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPerp,weight);
	    recoilResponse_nJets[jetBin]->Fill(responseHad,weight);
	    recoilResponse_Ptbins_nJets[jetBin][ptBin]->Fill(responseHad,weight);

	    // puppimet
	    recoilPuppiZParalH[jetBin]->Fill(recoilPuppiParal,weight);
            recoilPuppiZPerpH[jetBin]->Fill(recoilPuppiPerp,weight);
            recoilPuppiZParal_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPuppiParal,weight);
            recoilPuppiZPerp_Ptbins_nJetsH[jetBin][ptBin]->Fill(recoilPuppiPerp,weight);
	    recoilPuppiResponse_nJets[jetBin]->Fill(responsePuppiHad,weight);
	    recoilPuppiResponse_Ptbins_nJets[jetBin][ptBin]->Fill(responsePuppiHad,weight);

	  }

	  if (massSel<70||massSel>110) { // Z-depleted region

	    metTopSelH->Fill(pfmet,weight);
	    puppimetTopSelH->Fill(puppimet,weight);

	    metTopSelNJetsH[jetBin]->Fill(pfmet,weight);
	    puppimetTopSelNJetsH[jetBin]->Fill(puppimet,weight);

	    // pfmet
	    recoilTopParalH[jetBin]->Fill(recoilParal,weight);
	    recoilTopPerpH[jetBin]->Fill(recoilPerp,weight);

	    // puppimet
	    recoilPuppiTopParalH[jetBin]->Fill(recoilPuppiParal,weight);
	    recoilPuppiTopPerpH[jetBin]->Fill(recoilPuppiPerp,weight);

	  }

	}

      }

      if (isIsoMuonsPair) selEventsIsoMuons++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Sum of weights                                   = " << histWeightsSkimmedH->GetSumOfWeights() << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  //cout << "weight used:" << weight << std::endl;

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



