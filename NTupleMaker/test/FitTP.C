#include "HttStylesNew.cc"
#include "HtoH.h"
#include "FitPassAndFail.C"
// FitTP performs fitting of tag-and-probe histograms
// sequentially for different lepton pt pins and given 
// eta bin (defined by TString histBaseName)
// ****************************************************
// Tag-and-probe histogram naming convention 
// ZMass${etaBin}${ptBin}Pass - m(ll) distribution
//                              in sample of passing probes 
//                              in given (pt,eta) bin
//                              of probed lepton
// ZMass${etaBin}${ptBin}Fail - m(ll) distribution 
//                              in sample of failing probes
//                              in given (pt,eta) bin
//                              of probed lepton 
// **************************************************
// Files :
// SingleMuon_Run201D.root (used for muon Id/Iso)
// SingleElectron_Run2015D (used for electron Id/Iso)
// 
// Names of eta bins :
// ${etaBin} = EtaLt0p9, Eta0p9to1p2, EtaGt1p2 (muons)
// ${etaBin} = EtaLt1p48, EtaGt1p48 (electrons) 
// Names pt bins (both muons and electrons) :
// Id/Iso
// ${ptBin} = Pt10to15, Pt15to20, Pt20to25,
//            Pt25to30, Pt30to40, Pt40to60, PtGt60
// Names of pt bins are specified in the macro


void FitTP(TString fileName = "DYJetsToLL_M-50_MG", // RooT file with tag-&-probe histograms 
 	   TString histBaseName = "ZMassEtaGt1p2", // Basename of the histograms to fit (eta bin) 
	   TString xTitle = "muon p_{T}[GeV]", // x axis title (efficiency plot)
	   TString yTitle = "efficiency", // y axis title (efficiency plot)
	   TString xtit = "m_{#mu#mu}[GeV]", // title for dilepton mass distributions
	   float norm = 1 // luminosity normalization factor (1 for data) 
	   ) {
  
  TString SampleName("_MC");
  if (fileName.Contains("SingleMuon")) SampleName = "_Data";
  if (fileName.Contains("SingleElectron")) SampleName = "_Data";

   //std::cout << "ok 1" << std::endl;

  // define if in the fit of failing probes
  // the FSR component will be used in the 
  // signal function
  bool fitWithFSR[17];
  for (int i=0; i<17; ++i)
    fitWithFSR[i] = true;

  // binning for ID+ISO efficiency
  int nPtBins = 8; 
  float ptBins[17];
  ptBins[0] = 10;	 
  ptBins[1] = 13;
  ptBins[2] = 16;
  ptBins[3] = 20;
  ptBins[4] = 25; 
  ptBins[5] = 30; 
  ptBins[6] = 40;
  ptBins[7] = 60;
  ptBins[8] = 80;

  TString PtBins[17];
  PtBins[0]= "Pt10to13"; PtBins[1]= "Pt13to16"; PtBins[2]= "Pt16to20"; PtBins[3]= "Pt20to25"; PtBins[4]="Pt25to30"; 
  PtBins[5]= "Pt30to40"; PtBins[6]="Pt40to60"; PtBins[7]="PtGt60";

  //binning for trigger efficiency

	//std::cout << histBaseName << std::endl;
	//std::cout << histBaseName.Contains("ZMassEle") << std::endl;

  if ( histBaseName.Contains("ZMassEle") || histBaseName.Contains("ZMassIsoEle") || histBaseName.Contains("ZMassMu") || histBaseName.Contains("ZMassIsoMu") ) {
    nPtBins = 16;
    ptBins[0] = 10; ptBins[1] = 13; ptBins[2] = 16; ptBins[3] = 19; ptBins[4] = 22; ptBins[5] = 25; ptBins[6] = 28; ptBins[7] = 31; ptBins[8] = 34; ptBins[9] = 37;
    ptBins[10] = 40; ptBins[11] = 45; ptBins[12] = 50; ptBins[13] = 60; ptBins[14] = 70; ptBins[15] = 100; ptBins[16] = 150;     
    
    PtBins[0] = "Pt10to13"; PtBins[1] = "Pt13to16"; PtBins[2] = "Pt16to19"; PtBins[3] ="Pt19to22"; PtBins[4] ="Pt22to25"; PtBins[5] ="Pt25to28"; PtBins[6] ="Pt28to31"; 
    PtBins[7] ="Pt31to34"; PtBins[8] ="Pt34to37"; PtBins[9] ="Pt37to40"; PtBins[10] ="Pt40to45"; PtBins[11] ="Pt45to50"; PtBins[12] ="Pt50to60"; PtBins[13] ="Pt60to70";
    PtBins[14] ="Pt70to100"; PtBins[15] ="PtGt100";
  };

  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins);

  TFile * file = new TFile(fileName+".root");

   //std::cout << "ok 2" << std::endl;

  for (int iPt=0; iPt<nPtBins; ++iPt) { // loop over pt bins

    TH1F * histPassOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Pass");
 	//std::cout << histBaseName+PtBins[iPt]+"Pass" <<std::endl;
	int nBinsX = histPassOld->GetNbinsX();
    TH1F * histFailOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Fail");
 	//std::cout << histBaseName+PtBins[iPt]+"Fail" <<std::endl;
	//std::cout << histPassOld << std::endl;
	//std::cout << histFailOld << std::endl;
    //int nBinsX = histPassOld->GetNbinsX();
	//std::cout << "ok 3" << std::endl;

    for (int iB=1;iB<=nBinsX;++iB) {
      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
    }
	//std::cout << "ok 4" << std::endl;
    float output[2];
    TCanvas * c1 = new TCanvas("c1","",700,600);
    TCanvas * c2 = new TCanvas("c2","",700,600);
    bool fitPass = true;
    bool fitFail = true;
    bool rebinPass = false;
    bool rebinFail = false;
    //std::cout << "ok before fitting" << std::endl;
    FitPassAndFail(fileName,histBaseName+PtBins[iPt],xtit,histPassOld,histFailOld,fitPass,fitFail,fitWithFSR[iPt],rebinPass,rebinFail,c1,c2,output);
    c1->cd();
    c1->Update();
    c2->cd();
    c2->Update();
    numeratorH->SetBinContent(iPt+1,output[0]);
    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);
  }
  TFile * outputFile = new TFile(fileName+"_"+histBaseName+".root","recreate");
  outputFile->cd();
  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
  eff->Divide(numeratorH,denominatorH);
  eff->GetXaxis()->SetTitle(xTitle);
  //  eff->GetXaxis()->SetRangeUser(10.01,59.99);
  eff->GetYaxis()->SetRangeUser(0,1.0);
  eff->GetYaxis()->SetTitle(yTitle);
  eff->GetXaxis()->SetTitleOffset(1.1);
  eff->GetYaxis()->SetTitleOffset(1.1);
  eff->SetMarkerSize(1.7);
  eff->SetLineWidth(2);
  TCanvas * canv = new TCanvas("canv","",700,600);
  eff->Draw("APE");
  canv->Update();
  canv->SaveAs(fileName+"_"+histBaseName+".png");
  eff->Write(histBaseName+SampleName);

 // copy the eta bins histogram in the output file
  file->cd();
  TH1F * etaBinsH = (TH1F*)file->Get("etaBinsH");
  outputFile->cd(); 
  etaBinsH->Write();
  outputFile->Close();


}
