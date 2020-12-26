TH1D * Make1DHistoFromTree(TTree * tree, 
			   TString variable,
			   TString histName,
			   int nbins,
			   float * bins,
			   TString Cuts) {

  TH1D * hist = new TH1D(histName,"",nbins,bins);
  tree->Draw(variable+">>"+histName,Cuts);
  return hist;

}

TH1D * make1DEfficiency(TH1D * passH,
			TH1D * failH,
			TString name) {

  TH1D * effH = (TH1D*)passH->Clone(name);

  int nBins = effH->GetNbinsX();

  for (int iB=1; iB<=nBins; ++iB) {
    double xPass = passH->GetBinContent(iB);
    double ePass = passH->GetBinError(iB);
    double xFail = failH->GetBinContent(iB);
    double eFail = failH->GetBinError(iB);
    double xTotal = xPass + xFail;
    double err = 0;
    double eff = 0;
    if (xTotal>0) {
      double xTotal2 = xTotal*xTotal;
      double e1 = xFail*ePass/xTotal2;
      double e2 = eFail*xPass/xTotal2;
      err = TMath::Sqrt(e1*e1+e2*e2);
      eff = xPass / xTotal;
    }
    effH->SetBinContent(iB,eff);
    effH->SetBinError(iB,err);
  }
  return effH;

}

TString etaCut[3] = {
  "TMath::Abs(tauEta)<2.1",
  "TMath::Abs(tauEta)<1.48",
  "TMath::Abs(tauEta)>1.48"
};

TString etaLabel[3] = {
  "_etaLt2p1",
  "_barrel",
  "_endcap"
};

TString WPCut[3] = {
  "taubyVVLooseDeepTau2017v2p1VSe&&taubyLooseDeepTau2017v2p1VSmu", // tau+tau
  "taubyVVLooseDeepTau2017v2p1VSe&&taubyTightDeepTau2017v2p1VSmu", // mu+tau
  "taubyTightDeepTau2017v2p1VSe&&taubyLooseDeepTau2017v2p1VSmu" // e+tau
};

TString WPLabel[3] = {
  "_tt",
  "_mt",
  "_et"
};
