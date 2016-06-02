void ComputeAcceptance(TString fileName = "DYJetsToLL_M-50_amca_forced_ztautau.root") {

  TFile * file = new TFile(fileName);

  TH1D * histZTTGenWeightsH = (TH1D*)file->Get("histZTTGenWeightsH");
  TH1D * histGenCutsWeightsH = (TH1D*)file->Get("histGenCutsWeightsH");
  TH1D * histGenCutsGenWeightsH = (TH1D*)file->Get("histGenCutsGenWeightsH");

  TH1D * histRecCutsGenWeightsH = (TH1D*)file->Get("histRecCutsGenWeightsH");
  TH1D * histRecCutsWeightsH = (TH1D*)file->Get("histRecCutsWeightsH");
  TH1D * histBDTCutGenWeightsH = (TH1D*)file->Get("histBDTCutGenWeightsH");
  TH1D * histBDTCutWeightsH = (TH1D*)file->Get("histBDTCutWeightsH");


  double ninput = histZTTGenWeightsH->GetBinContent(1);
  double naccept = histGenCutsGenWeightsH->GetBinContent(1);
  double nacceptE = histGenCutsGenWeightsH->GetBinError(1);
  double nacceptCorr = histGenCutsWeightsH->GetBinContent(1);
  double nacceptCorrE = histGenCutsWeightsH->GetBinError(1);

  double nreccuts = histRecCutsGenWeightsH->GetBinContent(1);
  double nreccutsE = histRecCutsGenWeightsH->GetBinError(1);
  double nreccutsCorr = histRecCutsWeightsH->GetBinContent(1);
  double nreccutsCorrE = histRecCutsWeightsH->GetBinError(1);

  double nbdtcut = histBDTCutGenWeightsH->GetBinContent(1);
  double nbdtcutE = histBDTCutGenWeightsH->GetBinError(1);
  double nbdtcutCorr = histBDTCutWeightsH->GetBinContent(1);
  double nbdtcutCorrE = histBDTCutWeightsH->GetBinError(1);


  double acceptance = naccept/ninput;
  double acceptanceE = nacceptE/ninput;

  double acceptanceCorr  = nacceptCorr/ninput;
  double acceptanceCorrE = nacceptCorrE/ninput;

  double effreccuts = nreccuts/naccept;
  double effreccutsE = nreccutsE/naccept;

  double effreccutsCorr = nreccutsCorr/naccept;
  double effreccutsCorrE = nreccutsCorrE/naccept;

  double effbdtcut = nbdtcut/naccept;
  double effbdtcutE = nbdtcutE/naccept;

  double effbdtcutCorr  = nbdtcutCorr/naccept;
  double effbdtcutCorrE = nbdtcutCorrE/naccept;




  std::cout << std::endl;
  std::cout << "Z->tautau->mumu ---> " << std::endl;
  std::cout << std::endl;
  std::cout << "acceptance     (w/o corr) = " << acceptance << " +/- " << acceptanceE << std::endl;
  std::cout << "acceptance     (w corr)   = " << acceptanceCorr << " +/- " << acceptanceCorrE << std::endl;
  std::cout << std::endl;
  std::cout << "eff(kin. cuts) (w/o corr) = " << effreccuts << " +/- " << effreccutsE << std::endl;
  std::cout << "eff(kin. cuts) (w corr)   = " << effreccutsCorr << " +/- " << effreccutsCorrE << std::endl;
  std::cout << std::endl;
  std::cout << "eff(BDT > 0.5) (w/o corr) = " << effbdtcut << " +/- " << effbdtcutE << std::endl;
  std::cout << "eff(BDT > 0.5) (w corr)   = " << effbdtcutCorr << " +/- " << effbdtcutCorrE << std::endl;
  std::cout << std::endl;

}
