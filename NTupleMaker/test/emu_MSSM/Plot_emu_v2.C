#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include "settings.h"
TString SpecificCut(TString sample) {
  TString cut("");
  
  if (sample.Contains("MuonEG")) {
    if (sample.Contains("2016"))
      cut = "&&!((trg_singleelectron<0.5&&pt_1>26.)||(trg_singlemuon&&pt_2>23.))";
    if (sample.Contains("2017"))
      cut = "&&!((trg_singleelectron<0.5&&pt_1>28.)||(trg_singlemuon&&pt_2>25.))";
    if (sample.Contains("2018"))
      cut = "&&!((trg_singleelectron<0.5&&pt_1>33.)||(trg_singlemuon&&pt_2>25.))";
  }
  if (sample.Contains("EGamma")||sample.Contains("SingleElectron")) {
    if (sample.Contains("2018"))
      cut = "&&!(trg_singlemuon>0.5&&pt_2>25.)";
    if (sample.Contains("2017"))
      cut = "&&!(trg_singlemuon>0.5&&pt_2>25.)";
    if (sample.Contains("2016"))
      cut = "&&!(trg_singlemuon>0.5&&pt_2>23.)";
  }
  
  if (sample.Contains("WJetsToLNu")||sample.Contains("DYJetsToLL_M-50"))
    cut = "&&gen_noutgoing==0";

  return cut;

}

void Plot_emu_v2( bool embedded = true,
		  TString era = "2018") {

  SetStyle();

  // ****************************************
  // ****** Variable to plot ****************
  // ****************************************
  TString Variable = "mt_tot_puppi";
  //  TString Variable = "dmMVA_2";
  //  TString xtitle = "IP signif. (#mu) (pv+bs)";
  TString xtitle = "mT(tot) [GeV]";
  TString ytitle = "Events / 10 GeV";
  int nBins  =                   1;
  float xmin =                   0;
  float xmax =                4000;
  float yLower =               10;
  float scaleYUpper =          10;

  TString suffix = "";
  if (embedded) suffix = "_embedded";

  bool plotLegend = true;
  bool legRight = true;
  bool logY = false;
  bool logX = false;

  // ******** end of settings *********
  TString dirPostfix("etau");
  TString outputGraphics("figures");

  //  TString dir = "/nfs/dust/cms/user/rasp/HiggsCP/"+dirPostfix+"/"+era;
  //  TString dirData = "/nfs/dust/cms/user/rasp/HiggsCP/etau/"+era;

  TString dir = "/nfs/dust/cms/user/rasp/grid-jobs/emu_MSSM/"+era;

  std::cout << dir << std::endl;

  lumi_13TeV = "2018, 59.7 fb^{-1}";
  if (era=="2017")
    lumi_13TeV = "2017, 41.0 fb^{-1}";
  if (era=="2016")
    lumi_13TeV = "2016, 35.9 fb^{-1}";

  TString Weight("puweight*mcweight*prefiringweight*");
  TString WeightEmb("mcweight*");
  TString WeightEff("(effweightExcl/trigweightExcl)*");
  Weight += WeightEff;
  WeightEmb += WeightEff;

  TString WeightDY = Weight + "zptweight*";
  TString WeightTop = Weight + "topptweight*";

  //  TString WeightQCD("TMath::Min(3.0,qcdweight)*");
  TString WeightQCD("2.3*");

  TString WeightSS = Weight + WeightQCD;
  TString WeightEmbSS = WeightEmb + WeightQCD;
  TString WeightDYSS = WeightDY + WeightQCD;
  TString WeightTopSS = WeightTop + WeightQCD;

  TString CutsEMu("((trg_muhigh_elow>0.5&&pt_2>24.0&&pt_1>13.0)||(trg_ehigh_mulow>0.5&&pt_1>24.0&&pt_2>10.0))");

  TString CutsSingleLep = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>33.))";
  if (era=="2017")
    CutsSingleLep = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>28.))";
  if (era=="2016")
    CutsSingleLep = "((trg_singlemuon>0.5&&pt_2>23.)||(trg_singleelectron>0.5&&pt_1>26.))";

  TString Cuts = "(" + CutsEMu + "||" + CutsSingleLep + ")";
  //  TString Cuts = CutsEMu;
  //  TString Cuts = CutsSingleLep;

  Cuts += TString("&&iso_1<0.15&&iso_2<0.20&&extraelec_veto<0.5&&extramuon_veto<0.5&&dr_tt>0.3&&pt_1>15.0&&pt_2>15.&&puppipzeta<-40.0");

  TString CutsOS = Cuts + TString("&&os>0.5");
  TString CutsSS = Cuts + TString("&&os<0.5");
  TString CutsZTT_OS  = CutsOS + TString("&&gen_match_1==3&&gen_match_2==4");
  TString CutsZLL_OS  = CutsOS + TString("&&!(gen_match_1==3&&gen_match_2==4)");
  TString CutsZTT_SS  = CutsSS + TString("&&gen_match_1==3&&gen_match_2==4");
  TString CutsZLL_SS  = CutsSS + TString("&&!(gen_match_1==3&&gen_match_2==4)");


  double lumi = 59740;
  if (era=="2017")
    lumi = 41900;
  if (era=="2016")
    lumi = 35890;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  std::vector<TString> SingleElectron = SingleElectron_2018;
  std::vector<TString> SingleMuon = SingleMuon_2018;
  std::vector<TString> MuonEG = MuonEG_2018;

  std::vector<TString> EmbedSamples = EmbeddedElMu_2018;
  std::vector<TString> DYSamples = DYJets;
  std::vector<TString> WJetsSamples = WJets;
  std::vector<TString> EWKSamples = EWK;
  std::vector<TString> TTSamples = TT_EXCL;

  if (era=="2017") {
    SingleElectron = SingleElectron_2017;
    SingleMuon = SingleMuon_2017;
    MuonEG = MuonEG_2017;
    EmbedSamples = EmbeddedElMu_2017;
  }
  if (era=="2016") {
    SingleElectron = SingleElectron_2016;
    SingleMuon = SingleMuon_2016;
    MuonEG = MuonEG_2016;
    EmbedSamples = EmbeddedElMu_2016;
    TTSamples = TT_INCL;
  }

  std::vector<TString> DataSamples; DataSamples.clear();
				      
  for (unsigned int i=0; i<SingleElectron.size(); ++i)
    DataSamples.push_back(SingleElectron.at(i));
  
  for (unsigned int i=0; i<SingleMuon.size(); ++i)
    DataSamples.push_back(SingleMuon.at(i));

  for (unsigned int i=0; i<MuonEG.size(); ++i)
    DataSamples.push_back(MuonEG.at(i));

  std::cout << "Data samples : " << std::endl;
  for (unsigned int i=0; i<DataSamples.size(); ++i)
    std::cout << DataSamples.at(i) << std::endl;

  struct SampleAttributes {
    TString name;
    std::vector<TString> sampleNames;
    TString weight;
    TString weightSS;
    TString cuts;
    TString cutsSS;
    TH1D * hist;
    TH1D * histSS;
  };

  TH1D * dataH   = new TH1D("data_os_Hist","",nBins,xmin,xmax);
  TH1D * dataSSH = new TH1D("data_ss_Hist","",nBins,xmin,xmax);

  TH1D * zttH   = new TH1D("ztt_os_Hist","",nBins,xmin,xmax);
  TH1D * zttSSH = new TH1D("ztt_ss_Hist","",nBins,xmin,xmax);

  TH1D * zllH   = new TH1D("zll_os_Hist","",nBins,xmin,xmax);
  TH1D * zllSSH = new TH1D("zll_ss_Hist","",nBins,xmin,xmax);

  TH1D * wH   = new TH1D("w_os_Hist","",nBins,xmin,xmax);
  TH1D * wSSH = new TH1D("w_ss_Hist","",nBins,xmin,xmax);

  TH1D * ewkH   = new TH1D("ewk_os_Hist","",nBins,xmin,xmax);
  TH1D * ewkSSH = new TH1D("ewk_ss_Hist","",nBins,xmin,xmax);

  TH1D * ttH   = new TH1D("tt_os_Hist","",nBins,xmin,xmax);
  TH1D * ttSSH = new TH1D("tt_ss_Hist","",nBins,xmin,xmax);

  // data
  SampleAttributes DataAttr;
  DataAttr.name = "Data";
  DataAttr.sampleNames = DataSamples;
  DataAttr.weight = TString("1.0*");
  DataAttr.weightSS = WeightQCD;
  DataAttr.cuts = CutsOS;
  DataAttr.cutsSS = CutsSS;
  DataAttr.hist = dataH;
  DataAttr.histSS = dataSSH;

  // ZTT
  SampleAttributes ZttAttr;
  ZttAttr.name = "ZTT";
  if (embedded) {
    ZttAttr.sampleNames = EmbedSamples;
    ZttAttr.weight = WeightEmb;
    ZttAttr.weightSS = WeightEmbSS;
    ZttAttr.cuts = CutsOS;
    ZttAttr.cutsSS = CutsSS;
  }
  else {
    ZttAttr.sampleNames = DYSamples;
    ZttAttr.weight = WeightDY;
    ZttAttr.weightSS = WeightDYSS;
    ZttAttr.cuts = CutsZTT_OS;
    ZttAttr.cutsSS = CutsZTT_SS;
  }
  ZttAttr.hist = zttH;
  ZttAttr.histSS = zttSSH;

  // ZLL
  SampleAttributes ZllAttr;
  ZllAttr.name = "ZLL";
  ZllAttr.sampleNames = DYSamples;
  ZllAttr.weight = WeightDY;
  ZllAttr.weightSS = WeightDYSS;
  ZllAttr.cuts = CutsZLL_OS;
  ZllAttr.cutsSS = CutsZLL_SS;
  ZllAttr.hist = zllH;
  ZllAttr.histSS = zllSSH;

  // WJets
  SampleAttributes WJetsAttr;
  WJetsAttr.name = "WJets";
  WJetsAttr.sampleNames = WJetsSamples;
  WJetsAttr.weight = Weight;
  WJetsAttr.weightSS = WeightSS;
  if (embedded) {
    WJetsAttr.cuts = CutsZLL_OS;
    WJetsAttr.cutsSS = CutsZLL_SS;
  }
  else {
    WJetsAttr.cuts = CutsOS;
    WJetsAttr.cutsSS = CutsSS;
  }
  WJetsAttr.hist = wH;
  WJetsAttr.histSS = wSSH;

  // EWK
  SampleAttributes EWKAttr = WJetsAttr;
  EWKAttr.name = "EWK";
  EWKAttr.sampleNames = EWKSamples;
  EWKAttr.hist = ewkH;
  EWKAttr.histSS = ewkSSH;

  // TTBar
  SampleAttributes TTAttr = WJetsAttr;
  TTAttr.name = "TTBar";
  TTAttr.sampleNames = TTSamples;
  TTAttr.weight = WeightTop;
  TTAttr.weightSS = WeightTopSS;
  TTAttr.hist = ttH;
  TTAttr.histSS = ttSSH;

  std::vector<SampleAttributes> AllSamples;
  AllSamples.push_back(DataAttr);
  AllSamples.push_back(ZttAttr);
  AllSamples.push_back(ZllAttr);
  AllSamples.push_back(WJetsAttr);
  AllSamples.push_back(EWKAttr);
  AllSamples.push_back(TTAttr);

  TCanvas * dummy = new TCanvas("dummy","",400,400);

  // filling histograms
  for (unsigned int i=0; i<AllSamples.size(); ++i) {
    SampleAttributes sampleAttr = AllSamples.at(i);
    TString name = sampleAttr.name; 
    std::vector<TString> Samples = sampleAttr.sampleNames;
    for (unsigned int j=0; j<Samples.size(); ++j) {
      TString sampleName = Samples.at(j);
      std::cout << "Processing sample : " << sampleName << std::endl;
      TFile * file = new TFile(dir+"/"+sampleName+".root");
      TTree * tree = (TTree*)file->Get("TauCheck");
      TH1D * histWeightsH = (TH1D*)file->Get("nWeightedEvents");
      TString histName = sampleName + "_os";
      TString histNameSS = sampleName + "_ss";
      TH1D * histSample = new TH1D(histName,"",nBins,xmin,xmax);
      TH1D * histSampleSS = new TH1D(histNameSS,"",nBins,xmin,xmax);
      histSample->Sumw2();
      histSampleSS->Sumw2();
      TString CutsSample = sampleAttr.cuts + SpecificCut(sampleName);
      TString CutsSampleSS = sampleAttr.cutsSS + SpecificCut(sampleName);
      tree->Draw(Variable+">>"+histName,sampleAttr.weight+"("+CutsSample+")");
      tree->Draw(Variable+">>"+histNameSS,sampleAttr.weightSS+"("+CutsSampleSS+")");
      double norm = 1.0;
      if (name=="Data"||sampleName.Contains("Embed")) {
	norm = 1.;
	//	std::cout << sampleName //<< " : " << CutsSample 
	//		  << " : "  << histSample->GetSumOfWeights() << std::endl;
	  //	<< " : " << histSampleSS->GetSumOfWeights() << std::endl;
      }
      else { 
	double xsec = xsecs[sampleName];
	double nevents = histWeightsH->GetSumOfWeights();
	norm = xsec*lumi/nevents;
      }
      sampleAttr.hist->Add(sampleAttr.hist,histSample,1.,norm);
      sampleAttr.histSS->Add(sampleAttr.histSS,histSampleSS,1.,norm);
      //      delete file;
      //      delete histSample;
      //      delete histSampleSS;
      
    }
  }  
  
  //  return;

  // creating QCD
  DataAttr.histSS->Add(DataAttr.histSS,ZttAttr.histSS,1,-1);
  DataAttr.histSS->Add(DataAttr.histSS,ZllAttr.histSS,1,-1);
  DataAttr.histSS->Add(DataAttr.histSS,WJetsAttr.histSS,1,-1);
  DataAttr.histSS->Add(DataAttr.histSS,EWKAttr.histSS,1,-1);
  DataAttr.histSS->Add(DataAttr.histSS,TTAttr.histSS,1,-1);
  
  std::cout << "HERE" << std::endl;

  TFile * fileBBH = new TFile(dir+"/SUSYGluGluToBBHToTauTau_M-1200.root");
  TTree * treeBBH = (TTree*)fileBBH->Get("TauCheck");
  TH1D * histWeightsBBH = (TH1D*)fileBBH->Get("nWeightedEvents");
  TH1D * histBBH = new TH1D("bbH1200","",nBins,xmin,xmax);
  treeBBH->Draw(Variable+">>bbH1200",Weight+"*("+CutsOS+")");

  TFile * fileGGH = new TFile(dir+"/SUSYGluGluToHToTauTau_M-1200.root");
  TTree * treeGGH = (TTree*)fileGGH->Get("TauCheck");
  TH1D * histWeightsGGH = (TH1D*)fileGGH->Get("nWeightedEvents");
  TH1D * histGGH = new TH1D("ggH1200","",nBins,xmin,xmax);
  treeGGH->Draw(Variable+">>ggH1200",Weight+"*("+CutsOS+")");

  histBBH->Scale(5000000.0/histWeightsBBH->GetSumOfWeights());
  histBBH->SetLineColor(2);
  histGGH->Scale(5000000.0/histWeightsGGH->GetSumOfWeights());
  histGGH->SetLineColor(4);
  std::cout << "histBBH : " << histBBH->GetSumOfWeights() << std::endl;
  std::cout << "histGGH : " << histGGH->GetSumOfWeights() << std::endl;  

  for (int iB=1; iB<=nBins; ++iB) {
    histBBH->SetBinError(iB,0);
    histGGH->SetBinError(iB,0);
  }

  delete dummy;
  
  TH1D * histData = DataAttr.hist;
  TH1D * W        = WJetsAttr.hist;
  TH1D * TT       = TTAttr.hist;
  TH1D * EWK      = EWKAttr.hist;
  TH1D * ZLL      = ZllAttr.hist;
  TH1D * ZTT      = ZttAttr.hist;
  TH1D * QCD      = DataAttr.histSS;

  std::cout << "Top : " << TT->GetSumOfWeights() << std::endl;
  std::cout << "EWK : " << EWK->GetSumOfWeights() << std::endl;
  std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
  std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << std::endl;

  //  adding normalization systematics
  double ZTT_norm = 0.04; //  normalization ZTT :  5% (EMBEDDED)
  double EWK_norm = 0.10; //  normalization EWK : 10%
  double QCD_norm = 0.15; //  normalization Fakes : 15%
  double ZLL_mtau = 0.10; //  mu->tau fake rate ZLL : 10%
  double TT_norm  = 0.07; //  normalization TT  :  7%
  double mu_ID    = 0.02; //  muon trigger ID   :  2%
  double tau_ID   = 0.02; //  tau ID            :  5%
  double W_norm   = 0.15; //  normalization W   : 15%

  for (int iB=1; iB<=nBins; ++iB) {

    float ztt  = ZTT->GetBinContent(iB);
    float ztte = ZTT->GetBinError(iB);
    ztte = TMath::Sqrt(ztte*ztte+ztt*ztt*(ZTT_norm*ZTT_norm+mu_ID*mu_ID+tau_ID*tau_ID));
    ZTT->SetBinError(iB,ztte);

    float ewk  = EWK->GetBinContent(iB);
    float ewke = EWK->GetBinError(iB);
    ewke = TMath::Sqrt(ewke*ewke+ewk*ewk*EWK_norm*EWK_norm);
    EWK->SetBinError(iB,ewke);

    float qcd  = QCD->GetBinContent(iB);
    float qcde = QCD->GetBinError(iB);
    qcde = TMath::Sqrt(qcde*qcde+qcd*qcd*QCD_norm*QCD_norm);
    QCD->SetBinError(iB,qcde);
    if (qcd<0) {
      QCD->SetBinContent(iB,0.);
      QCD->SetBinError(iB,0.);
    }

    float w = W->GetBinContent(iB);
    float we = W->GetBinError(iB);
    we = TMath::Sqrt(we*we+w*w*(W_norm*W_norm+mu_ID*mu_ID));
    W->SetBinError(iB,we);

    float tt  = TT->GetBinContent(iB);
    float tte = TT->GetBinError(iB);
    tte = TMath::Sqrt(tte*tte+tt*tt*TT_norm*TT_norm);
    TT->SetBinError(iB,tte);

    float zll  = ZLL->GetBinContent(iB);
    float zlle = ZLL->GetBinError(iB);
    zlle = TMath::Sqrt(zlle*zlle+zll*zll*(mu_ID*mu_ID+ZLL_mtau*ZLL_mtau));
    ZLL->SetBinError(iB,zlle);

    /*
    std::cout << iB << " : " 
	      << "qcd = " << qcd << " +/- " << qcde 
	      << "  w = " << w << " +/- " << we 
	      << "  ztt = " << ztt << " +/- " << ztte
	      << "  zll = " << zll << " +/- " << zlle 
	      << "  ewk = " << ewk << " +/- " << ewke << std::endl;
    */
  }

  EWK->Add(EWK,TT);
  W->Add(W,EWK);
  QCD->Add(QCD,W);
  ZLL->Add(ZLL,QCD);
  ZTT->Add(ZTT,ZLL);

  std::cout << std::endl;
  std::cout << "MC : " << ZTT->GetSumOfWeights() << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << std::endl;

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  InitData(histData);

  InitHist(QCD,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(EWK,"","",TColor::GetColor("#6F2D35"),1001);
  InitHist(W,"","",TColor::GetColor("#DE5A6A"),1001);
  InitHist(ZLL,"","",TColor::GetColor("#4496C8"),1001);
  InitHist(ZTT,"","",TColor::GetColor("#FFCC66"),1001);

  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  float yUpper = histData->GetMaximum();
  if (ZTT->GetMaximum()>yUpper)
    yUpper = ZTT->GetMaximum();
  histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);
  ZTT->GetYaxis()->SetRangeUser(0,1.2*ZTT->GetMaximum());
  if (logY) {
    histData->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
    ZTT->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
  }

  histData->SetMarkerSize(1.4);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);

  histData->Draw("e1");
  ZTT->Draw("sameh");
  ZLL->Draw("sameh");
  QCD->Draw("sameh");
  W->Draw("sameh");
  EWK->Draw("sameh");
  TT->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  histData->Draw("e1same");
  histBBH->Draw("same");
  histGGH->Draw("same");
  float chi2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = W->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg;
  if (legRight) leg = new TLegend(0.61,0.37,0.90,0.77);
  else leg = new TLegend(0.20,0.47,0.50,0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  leg->AddEntry(histData,"Data","lp");
  if (embedded) 
    leg->AddEntry(ZTT,"embedded Z#rightarrow#tau#tau","f");
  else 
    leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow ll","f");
  leg->AddEntry(QCD,"QCD","f");
  leg->AddEntry(W,"W#rightarrow e#nu","f");
  leg->AddEntry(EWK,"electroweak","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  leg->AddEntry(histBBH,"bbH(1200)","l");
  leg->AddEntry(histGGH,"ggH(1200)","l");
  if (plotLegend) leg->Draw();
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);
  if (logX) upper->SetLogx(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.601,1.399);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.5);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  for (int iB=1; iB<=nBins; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = ZTT->GetBinContent(iB);
    ratioErrH->SetBinContent(iB,1.0);
    ratioErrH->SetBinError(iB,0.0);
    float xBkg = bkgdErr->GetBinContent(iB);
    float errBkg = bkgdErr->GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);

    }
    if (x1>0&&x2>0) {
      float e1 = histData->GetBinError(iB);
      float ratio = x1/x2;
      float eratio = e1/x2;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }

  // ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);

  ratioH->Draw("e1");
  ratioErrH->Draw("e2same");
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  if (logX) lower->SetLogx(true);
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Update();
  canv1->Print("plot_"+Variable+suffix+"_"+era+"_emu.png");

}
