#include "HttStylesNew.cc"
#include "HtoH.h"
#include "FitRecoil.C"

void FitRecoilSeq() {

  SetStyle();
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);

  TString dir("./final/");

  TString fileNameData(dir+"SingleMuon_Run2016.root");
  TString fileNameMC(dir+"DYJetsToLL_M-50_13TeV_madgraphMLM_pythia8_zpt.root");

  TString baseString("recoilZ");

  TString Proj[2] = {baseString+"Perp",
		     baseString+"Paral"};
  
  TString xtit[2] = {"U_{2} [GeV]","U_{1} [GeV]"};
  TString ytit[2] = {"Events",
		     "Events"};

  TString jetBins[3] = {"NJet0",
			"NJet1",
			"NJetGe2"};

  TString ptBins[5] = {"Pt0to10",
		       "Pt10to20",
		       "Pt20to30",
		       "Pt30to50",
		       "PtGt50"};

  int nZPtBins = 5;
  float zPtBinsX[6] = {0,10,20,30,50,1000};
  int nJetBins = 3;
  float jetBinsX[4] = {-0.5,0.5,1.5,2.5};

  TFile * fileOutput = new TFile(baseString+".root","recreate");
  fileOutput->cd("");
  TH1D * projH = new  TH1D("projH","",2,-1,1);
  TH1D * zPtBinsH = new TH1D("ZPtBinsH","",nZPtBins,zPtBinsX);
  TH1D * nJetBinsH = new TH1D("nJetBinsH","",nJetBins,jetBinsX);
  for (int i=0; i<2; ++i)
    projH->GetXaxis()->SetBinLabel(i+1,Proj[i]);
  
  for (int i=0; i<nJetBins; ++i)
    nJetBinsH->GetXaxis()->SetBinLabel(i+1,jetBins[i]);

  for (int i=0; i<nZPtBins; ++i)
    zPtBinsH->GetXaxis()->SetBinLabel(i+1,ptBins[i]);
  
  projH->Write("projH");
  zPtBinsH->Write("ZPtBinsH");
  nJetBinsH->Write("nJetBinsH");


  TString samples[10] = {
    "WW_13TeV-pythia8",                
    "WZ_13TeV-pythia8",                
    "ZZ_13TeV-pythia8",                
    "WJetsToLNu_13TeV-madgraphMLM-pythia8", 
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8", 
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8",
    "ST_t-channel_top_4f_InclusiveDecays_13TeV-powheg-pythia8", 
    "ST_t-channel_antitop_4f_InclusiveDecays_13TeV-powheg-pythia8",
    "TT_13TeV-powheg-pythia8",        
    "SingleMuon_Run2016_ss"
  };

  float xsec[10] = {118.7,
                    27.68,
                    12.19,
                    52760*1.166,
                    35.85,
                    35.85,
                    44.33,
                    26.38,
                    831.8,
		    2};    

  TFile * fileData   = new TFile(fileNameData);
  TFile * fileMC     = new TFile(fileNameMC);
  
  float xsecDY = 5345*1.079; // Drell Yan cross section
  float lumi = 35900;  // integrated luminosity

  TH1D * weightMC = (TH1D*)fileMC->Get("histWeightsH");
  float normMC = xsecDY*lumi/weightMC->GetSumOfWeights();

  for (int iProj=0; iProj<2; ++iProj) {
    for (int iJet=0; iJet<3; ++iJet) {
      for (int iPt=0; iPt<5; ++iPt) {
	TString histName = Proj[iProj] + "_" + jetBins[iJet] + ptBins[iPt];
	cout << histName << endl;
	TH1D * histData = (TH1D*)fileData->Get(histName);
	TH1D * histMC = (TH1D*)fileMC->Get(histName);
	histMC->Scale(normMC);

	for (int iS=0; iS<10; ++iS) {
	  
	  TFile * fileSample = new TFile(dir+samples[iS]+".root");
	  TH1D * weightsH = (TH1D*)fileSample->Get("histWeightsH");
	  double normSample = lumi*xsec[iS]/weightsH->GetSumOfWeights();
	  TH1D * histSample = (TH1D*)fileSample->Get(histName);
	  if (iS==11) normSample = 2;
	  histData->Add(histData,histSample,1,-normSample);

	  delete fileSample;

	}
	bool rebin = false;
	float xrange = 120; 
	bool logY = false;
	bool asymGauss = true;
	//	if (iJet==2) rebin = true;
	if (iProj==0) asymGauss = false;
	if (iJet==0) xrange = 120;
	FitRecoil(histData,histMC,fileOutput,histName,-xrange,xrange,xtit[iProj],ytit[iProj],-180,180,asymGauss,rebin,logY);
      }
    }
  }

  fileOutput->Close();

}
