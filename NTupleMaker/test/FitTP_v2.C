#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "DesyTauAnalyses/NTupleMaker/test/HtoH.h"
#include "DesyTauAnalyses/NTupleMaker/test/FitPassAndFail.C"
// FitTP performs fitting of tag-and-probe histograms
// sequentially for different lepton pt bins and eta bins. 
// ****************************************************
// Tag-and-probe histogram naming convention 
// For IdIso efficiencies
// ZMass${etaBin}${ptBin}Pass - m(ll) distribution
//                              in sample of passing probes 
//                              in given (pt,eta) bin
//                              of probed lepton
// ZMass${etaBin}${ptBin}Fail - m(ll) distribution 
//                              in sample of failing probes
//                              in given (pt,eta) bin
//                              of probed lepton 
// For trigger efficiencies, 
// ZMass${triggerName}${etaBin}${ptBin}
// **************************************************
// Names of pt and eta bins are taken from histograms in 
// the input file. 
// A maximum of 30 bins is allowed, to increase it modify  
// const int NMAX below. 
// **************************************************


// **************************************************
// Input Arguments
// fileName : input file name with TP histos, ptbins and eta bins
// what     : name of the efficiency you want to calculate. Please use 
//	      IdIso, Ele17, Ele12, Mu8, Mu17, IsoMu
// lepton   : Muon or Electron (case sensitive!)
// **************************************************


void FitTP_v2(TString fileName = "SingleMuon_2015D", // RooT file with tag-&-probe histograms 
 	   TString what = "IdIso", // or Mu8, IsoMu, ...
	   TString lepton = "Muon", //or Electron 
	   float norm = 1 // luminosity normalization factor (1 for data) 
	   ) 
{

	// name of output file 
	TString OutFileName = fileName+"_"+lepton+"_"+what+"_eff";


	// Title of axis in plots  
	TString yTitle = "Efficiency";
	TString xTitle = lepton+"  p_{T}[GeV]";
	TString xtit; 
	if (lepton == "Muon") xtit = "m_{#mu#mu}[GeV]"; 
	if (lepton == "Electron") xtit = "m_{ee}[GeV]"; 

  //names of final graphs - suffix
  TString SampleName("_MC");
  if (fileName.Contains("SingleMuon")) SampleName = "_Data";
  if (fileName.Contains("SingleElectron")) SampleName = "_Data";

  //open input file
  TFile * file = new TFile(fileName+".root");
  file->cd();

  //take eta and pt bins
  TH1F * etaBinsH = (TH1F*)file->Get("etaBinsH");
  int nEtaBins = etaBinsH->GetNbinsX();
  const int NMAX = 30;
  TString EtaBins[NMAX]; for (int i=0; i<nEtaBins; i++) {EtaBins[i]=etaBinsH->GetXaxis()->GetBinLabel(i+1);}

  TH1F * ptBinsTrigH = (TH1F*)file->Get("ptBinsTrigH");
  TH1F * ptBinsH = (TH1F*)file->Get("ptBinsH");
  int nPtBins;

  if (what == "IdIso") nPtBins = ptBinsH->GetNbinsX();
  else nPtBins = ptBinsTrigH->GetNbinsX();
  TString PtBins[NMAX]; 
  float ptBins_edges[nPtBins+1];

  if (what == "IdIso") {
			for (int i=0; i<nPtBins; i++) 
			{ PtBins[i]=ptBinsH->GetXaxis()->GetBinLabel(i+1); ptBins_edges[i]=ptBinsH->GetBinLowEdge(i+1); }
			ptBins_edges[nPtBins]= ptBinsH->GetBinLowEdge(nPtBins+1); 
  }
  else {
	for (int i=0; i<nPtBins; i++) 
	{ PtBins[i]=ptBinsTrigH->GetXaxis()->GetBinLabel(i+1); ptBins_edges[i]=ptBinsTrigH->GetBinLowEdge(i+1);  }
	ptBins_edges[nPtBins]= ptBinsTrigH->GetBinLowEdge(nPtBins+1); 
  }     

  // define if in the fit of failing probes
  // the FSR component will be used in the 
  // signal function 
  bool fitWithFSR[NMAX];	
  if (lepton == "Muon") { for (int i=0; i<nPtBins; i++)  fitWithFSR[i] = true;}

  if (lepton == "Electron") { for (int i=0; i<nPtBins; i++)  fitWithFSR[i] = false;  fitWithFSR[2] = true;  fitWithFSR[3] = true;}

	  

// building the histogram base name
  TString prefix = "ZMass";
  TString which = what; if (what == "IdIso") which = "";
  TString histBaseName; 

// define output file

 TFile * outputFile = new TFile(OutFileName+".root","recreate");
  

for (int iEta = 0; iEta < nEtaBins; iEta++) {  // loop over eta regions

  histBaseName = prefix+which+EtaBins[iEta];	
  //std::cout << "histBaseName" << histBaseName << std::endl;

  TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins_edges);
  TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins_edges);

  for (int iPt=0; iPt<nPtBins; ++iPt) { // loop over pt bins
    //std::cout <<" iPt = " << iPt << " | " << histBaseName+PtBins[iPt] <<std::endl;
    
    TH1F * histPassOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Pass");   // prefix+what+etaBin+ptBin+pass, fail
    int nBinsX = histPassOld->GetNbinsX();
    TH1F * histFailOld = (TH1F*)file->Get(histBaseName+PtBins[iPt]+"Fail");

    for (int iB=1;iB<=nBinsX;++iB) {
      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
    }

    float output[2];
    TCanvas * c1 = new TCanvas("c1","",700,600);
    TCanvas * c2 = new TCanvas("c2","",700,600);
    bool fitPass = true;
    bool fitFail = true;
    bool rebinPass = false;
    bool rebinFail = false;

    FitPassAndFail(fileName,histBaseName+PtBins[iPt],xtit,histPassOld,histFailOld,fitPass,fitFail,fitWithFSR[iPt],rebinPass,rebinFail,c1,c2,output);
    c1->cd();
    c1->Update();
    c2->cd();
    c2->Update();
    numeratorH->SetBinContent(iPt+1,output[0]);
    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);

  } // close loop over pt bins

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

} // close loop over eta regions

  outputFile->cd(); 
  etaBinsH->Write();
  outputFile->Close();


}
