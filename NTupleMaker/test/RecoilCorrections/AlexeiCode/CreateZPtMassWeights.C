#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"
// 

void CreateZPtMassWeights(TString histName = "dimuonMassPtH",
			  TString xtitle = "dimuon mass [GeV]",
			  TString ytitle = "dimuon p_{T} [GeV]",
			  bool logY = true,
			  bool drawLeg = false) {
  float qcdScale = 2;
  float ttScale = 1;
  float DYscaleLow = 1.0;
  float DYscaleHigh = 1;
  TString suffixHigh = "-madgraphMLM";
  TString suffixLow  = "-madgraphMLM";
  TString suffix = "";
  TString dir("./lepSF/");


  int nMassBins = 4;
  double massBins[5] = {50,100,200,500,1000};
  int nPtBins = 9;
  double ptBins[10] = {0,10,20,50,100,150,200,300,400,1000};

  //  SetStyle();
  gStyle->SetOptStat(0000);

  TFile * file = new TFile(dir+"SingleMuon_Run2018.root");
  TFile * fileSS = new TFile(dir+"SingleMuon_Run2018_ss.root");  
  
  TString samples[12] = {"WW_TuneCP5_13TeV-pythia8",                 // (0)
			 "WZ_TuneCP5_13TeV-pythia8",                 // (1) 
			 "ZZ_TuneCP5_13TeV-pythia8",                 // (2)
			 "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",  // (3)
			 "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",          // (4)
			 "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",      // (5)
			 "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",           // (6)
			 "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",       // (7)
			 "TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8",          // (8)
			 "TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8",   // (9)
			 "TTToHadronic_TuneCP5_13TeV_powheg_pythia8",       // (10)
			 "DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"+suffix   // (11)
  };

  float xsec[12] = {118.7,       // WW   (0)
		    27.68,       // WZ   (1)
		    12.19,       // ZZ   (2)
		    52760*1.166, // W+Jets (3)
		    35.85,       // tW_top (4)
		    35.85,       // tW_antitop (5) 
		    44.33,       // top-tchannel (6)
		    26.38,       // antitop-tchannel (7)
		    88.29,       // TTTo2L2Nu      (8)
		    365.34,      // TTToSemiLeptonic (9)
	            377.96,      // TTToHadronic (10)
		    5345*1.079   // DYJets (11)
  };

  float lumi = 60000;

  TH2D * histDataOld = (TH2D*)file->Get(histName);
  TH2D * histDataOldSS = (TH2D*)fileSS->Get(histName);

  int nBinsX = histDataOld->GetNbinsX();
  int nBinsY = histDataOld->GetNbinsY();

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << " nBinsX = " << nBinsX << "  nBinsY = " << nBinsY << std::endl;
  std::cout << std::endl;

  TH2D * dataHist = (TH2D*)TH2DtoTH2D(histDataOld,  nMassBins,massBins,nPtBins,ptBins,"_Data_new");
  TH2D * qcdHist  = (TH2D*)TH2DtoTH2D(histDataOldSS,nMassBins,massBins,nPtBins,ptBins,"_qcd");

  TH2D * mcHist = new TH2D("mcHist","",nMassBins,massBins,nPtBins,ptBins);
  TH2D * zHist  = new TH2D("zHist","" ,nMassBins,massBins,nPtBins,ptBins);

  int nSamples = 12;

  //  return;
  for (int iS=0; iS<nSamples; ++iS) {
    TFile * fileMC = new TFile(dir+samples[iS]+".root");
    TH2D * histOld = (TH2D*)fileMC->Get(histName);
    TH2D * hist = (TH2D*)TH2DtoTH2D(histOld,nMassBins,massBins,nPtBins,ptBins,"_new_"+samples[iS]);
    TH1D * eventCount = (TH1D*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[iS]*lumi/nGen;
    TH2D * tempHist = mcHist;
    if (iS==11)
      tempHist = zHist;
    tempHist->Add(tempHist,hist,1.,norm);
  }

  //  float dataEvents = 0;
  //  float ttEvents = 0;
  std::cout << "QCD        = " << qcdScale*qcdHist->GetSumOfWeights() << std::endl;
  std::cout << "non-QCD    = " << mcHist->GetSumOfWeights() << std::endl;
  std::cout << "Z->ll      = " << zHist->GetSumOfWeights() << std::endl;
  float allMC = qcdScale*qcdHist->GetSumOfWeights() + mcHist->GetSumOfWeights() + zHist->GetSumOfWeights();
  float allData = dataHist->GetSumOfWeights();
  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "Total bkgd = " << allMC << std::endl;
  std::cout << "Total data = " << allData << std::endl;
  std::cout << std::endl;

  mcHist->Add(mcHist,qcdHist,1,qcdScale);
  dataHist->Add(dataHist,mcHist,1,-1);
  double dataHistNorm = dataHist->GetSumOfWeights();
  double zHistNorm = zHist->GetSumOfWeights();

  TH2D * ratioHist = (TH2D*)dataHist->Clone("ratioHist");
  ratioHist->GetXaxis()->SetTitle(xtitle);
  ratioHist->GetYaxis()->SetTitle(ytitle);

  for (int iMass=1; iMass<=nMassBins; ++iMass) {
    for (int iPt=1; iPt<=nPtBins; ++iPt) {
      float xData = dataHist->GetBinContent(iMass,iPt)/dataHistNorm;
      float xZ = zHist->GetBinContent(iMass,iPt)/zHistNorm;
      float ratio = xData/xZ;
      ratioHist->SetBinContent(iMass,iPt,ratio);
      float xL = ratioHist->GetXaxis()->GetBinLowEdge(iMass);
      float xU = ratioHist->GetXaxis()->GetBinLowEdge(iMass+1);
      float yL = ratioHist->GetYaxis()->GetBinLowEdge(iPt);
      float yU = ratioHist->GetYaxis()->GetBinLowEdge(iPt+1);
      printf("[%3i-%3i,%3i-%3i] = %5.3f\n",int(xL+0.1),int(xU+0.1),int(yL+0.1),int(yU+0.1),ratio);
    }
  }

  

  TCanvas * canv1 = MakeCanvas("canv1", "", 900, 700);
  canv1->SetRightMargin(0.2);
  ratioHist->Draw("colztext");
  


  canv1->SetLogy(true);
  canv1->SetLogx(true);
  canv1->Print(histName+"_ratio.png");

  TFile * fileO = new TFile("zpt_weights_2018.root","recreate");
  fileO->cd("");
  ratioHist->Write("zptmass_histo");
  fileO->Close();

}
