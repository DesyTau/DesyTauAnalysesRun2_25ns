void ProduceDatacards_em(int nBins = 35,
			 float xmin = 0,
			 float xmax = 350,
			 bool pileUp = false) {

  double lumi = 2092;
  double normSS = 1.06;
    
  // sample names
  TString sampleNames[20] = {
    "MuEG_2015D", // data
    "DYJetsToLL_M-50_MG", // isZTT (ZTT)
    "DYJetsToLL_M-50_MG", // !isZTT (ZLL)
    "WJetsToLNu_MG", 
    "TTPowHeg",  
    "ST_t-channel_top_4f_leptonDecays",  
    "ST_t-channel_antitop_4f_leptonDecays",
    "ST_tW_antitop_5f_inclusiveDecays",
    "ST_tW_top_5f_inclusiveDecays",
    "VVTo2L2Nu",
    "WWToLNuQQ",
    "WZTo2L2Q",
    "WZTo1L1Nu2Q",
    "WZTo1L3Nu",
    "WZJets",
    "ZZTo4L",
    "ZZTo2L2Q",
    "",
    "",
    ""};

  double xsec[20] = {1, // data (0)
		     6025.2, // DY (1)
		     6025.2, // DY (2)
		     61526.7, // WJets (3)
		     831.76, // TT (4)
		     136.95, // ST_t-channel_top (5)
		     80.95,  // ST_t-channel_antitop (6)
		     35.6, // ST_tW_antitop (7)
		     35.6, // ST_tW_top_5f (8)
		     11.95,  // VV (9)
		     49.997, // WWToLNuQQ (10)
		     5.595, // WZTo2L2Q (11)
		     10.71, // WZTo1L1Nu2Q (12)
		     3.05, // WZTo1L3Nu (13)
		     5.26, // WZJets (3L1Nu) (14)
		     1.212, // ZZTo4L (15)
		     3.22, // ZZTo2L2Q (16)
		     0, //
		     0, //
		     0}; // 
  
  TString cuts[20];
  TString cutsSS[20];
  for (int i=0; i<20; ++i) {
    cuts[i] = "mcweight*(os>0.5)";
    cutsSS[i] = "mcweight*(os<0.5)";
  }

  cuts[0] = "os>0.5";
  cuts[1] = "mcweight*(os>0.5&&isZTT)";
  cuts[2] = "mcweight*(os>0.5&&!isZTT)";

  cutsSS[0] = "os<0.5";
  cutsSS[1] = "mcweight*(os<0.5&&isZTT)";
  cutsSS[2] = "mcweight*(os<0.5&&!isZTT)";
  
  if (pileUp) {
    for (int i=1; i<20; ++i) {
      cuts[i] = "puweight*" + cuts[i];
      cutsSS[i] = "puweight*" + cutsSS[i];
    }
  }

  TH1D * hist[20];
  TH1D * histSS[20];

  for (int i=0; i<17; ++i) {
    std::cout << sampleNames[i] << std::endl;
    TFile * file = new TFile(sampleNames[i]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
    TString histName = sampleNames[i] + "_mvis";
    TString histNameSS = sampleNames[i] + "_mvis_os";
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
    hist[i]->Sumw2();
    histSS[i]->Sumw2();
    tree->Draw("m_vis>>"+histName,cuts[i]);
    tree->Draw("m_vis>>"+histNameSS,cutsSS[i]);
    if (i>0) {
      for (int iB=1; iB<=nBins; ++iB) {
	double x = hist[i]->GetBinContent(iB);
	double e = hist[i]->GetBinError(iB);
    	hist[i]->SetBinContent(iB,norm*x);
    	hist[i]->SetBinError(iB,norm*e);
	double xSS = histSS[i]->GetBinContent(iB);
	double eSS = histSS[i]->GetBinError(iB);
    	histSS[i]->SetBinContent(iB,norm*xSS);
    	histSS[i]->SetBinError(iB,norm*eSS);
      }
    }
  }

  //  adding up single top and VV backgrounds
  for (int iH=6; iH<17; ++iH) {
    hist[5]->Add(hist[5],hist[iH]);
    histSS[5]->Add(histSS[5],histSS[iH]);
  }

  // subtracting background from SS
  for (int iH=1; iH<6; ++iH) {
    histSS[0]->Add(histSS[0],histSS[iH],1,-1);
  }

  for (int iB=1; iB<=nBins; ++iB) {
    double x = histSS[0]->GetBinContent(iB);
    double e = histSS[0]->GetBinError(iB);
    histSS[0]->SetBinContent(iB,normSS*x);
    histSS[0]->SetBinError(iB,normSS*e);
  }
  TString suffix("");
  if (!pileUp)
    suffix = "_noPU";

  TFile * file = new TFile("htt_em.inputs-sm-13TeV-mvis"+suffix+".root","recreate"); 
  file->mkdir("em_inclusive");
  file->cd("em_inclusive");
  TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
  TH1D * ZTT = (TH1D*)hist[1]->Clone("ZTT");
  TH1D * ZLL = (TH1D*)hist[2]->Clone("ZLL");
  TH1D * W = (TH1D*)hist[3]->Clone("W");
  TH1D * TT = (TH1D*)hist[4]->Clone("TT");
  TH1D * VV = (TH1D*)hist[5]->Clone("VV");
  file->Write();
  file->Close();
  

}
