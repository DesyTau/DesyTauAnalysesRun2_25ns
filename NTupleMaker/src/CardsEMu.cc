#ifndef CardsEMu_ccx
#define CardsEMu_ccx

#include "DesyTauAnalyses/NTupleMaker/interface/CardsEMu.h"

CardsEMu::CardsEMu(TString Sample,
		   TString Era,
		   TString Category,
		   TString InputDir,
		   TString OutputDir,
		   TString Variable,
		   int nbins,
		   double xmin,
		   double xmax,
		   int TriggerOption,
		   bool RunWithSystematics,
		   bool RunOnEmbedded) {

  era = Era;
  input_dir = InputDir+"/"+era;
  output_filename = OutputDir+"/"+era+"/"+Sample+"_"+Category+".root";
  sampleToProcess = Sample;
  category = Category;
  triggerOption = TriggerOption;
  variable = Variable;
  runWithSystematics = RunWithSystematics;
  runOnEmbedded = RunOnEmbedded;

  block = false;

  if (xmin>xmax) {
    std::cout << "warning -> xmin for histogram  : " << xmin << " is greater than xmax : " << xmax << std::endl;
    std::cout << "nothing will be done..." << std::endl;
    block = true;
    return;
  }

  nBins = nbins;
  double binWidth = (xmax-xmin)/double(nBins);
  for (int iB=0; iB<=nBins; ++iB)
    Bins[iB] = xmin + double(iB)*binWidth;


  if (sampleToProcess!="Data"&&
      sampleToProcess!="EWK"&&
      sampleToProcess!="TTbar"&&
      sampleToProcess!="DYJets"&&
      sampleToProcess!="EMB"&&
      sampleToProcess!="SMHiggs"&&
      sampleToProcess!="MSSMHiggs") {
    std::cout << "Unknown sample specified : " << sampleToProcess << std::endl;
    std::cout << "Nothing will be done" << std::endl;
    block = true;
    return;
  }  

  if (category!="em_inclusive"&&
      category!="em_btag"&&
      category!="em_nobtag"&&
      category!="em_nobtag_highmsv"&&
      category!="em_ttbar_control"&&
      category!="em_ttbar_btag"&&
      category!="em_ttbar_nobtag"&&
      category!="em_nobtag_highdzeta"&&
      category!="em_nobtag_mediumdzeta"&&
      category!="em_nobtag_lowdzeta"&&
      category!="em_nobtag_highmsv_highdzeta"&&
      category!="em_nobtag_highmsv_mediumdzeta"&&
      category!="em_nobtag_highmsv_lowdzeta"&&
      category!="em_btag_highdzeta"&&
      category!="em_btag_mediumdzeta"&&
      category!="em_btag_lowdzeta") {
    std::cout << "Unknown category specified : " << category << std::endl;
    std::cout << "nothing will be done" << std::endl;
    block = true;
    return;
  }

  std::map<TString,TString> categoryCuts = {
    {"em_inclusive","&&pzeta>-35.0"},
    // ---
    {"em_btag","&&pzeta>-35.0&&nbtag>0"},
    {"em_nobtag","&&pzeta>-35.0&&nbtag==0"},
    {"em_nobtag_highmsv","&&pzeta>-35.0&&nbtag==0&&m_sv>250"},
    // ---
    {"em_ttbar_control","&&pzeta<-35.0"},
    {"em_ttbar_btag","&&pzeta<-35.0&&nbtag>0"},
    {"em_ttbar_nobtag","&&pzeta<-35.0&&nbtag==0"},
    // ---
    {"em_nobtag_highdzeta","&&pzeta>30.0&&nbtag==0"},
    {"em_nobtag_mediumdzeta","&&pzeta>-10.0&&pzeta<30.0&&nbtag==0"},
    {"em_nobtag_lowdzeta","&&pzeta>-35.0&&pzeta<-10.0&&nbtag==0"},
    // ---
    {"em_nobtag_highmsv_highdzeta","&&pzeta>30.0&&nbtag==0&&m_sv>250.0"},
    {"em_nobtag_highmsv_mediumdzeta","&&pzeta>-10.0&&pzeta<30.0&&nbtag==0&&m_sv>250.0"},
    {"em_nobtag_highmsv_lowdzeta","&&pzeta>-35.0&&pzeta<-10.0&&nbtag==0&&m_sv>250.0"},
    // ---
    {"em_btag_highdzeta","&&pzeta>30.0&&nbtag>0"},
    {"em_btag_mediumdzeta","&&pzeta>-10.0&&pzeta<30.0&&nbtag>0"},
    {"em_btag_lowdzeta","&&pzeta>-35.0&&pzeta<-10.0&&nbtag>0"},
  };

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  isEquidistantBinning = true;
  if (variable=="mt_tot")
    isEquidistantBinning = false;

  isBTag = false;
  if (Category.Contains("_btag"))
    isBTag = true;

  outputFile = new TFile(output_filename,"recreate");
  if (outputFile->IsZombie()) {
    std::cout << "Cannot open file" << output_filename << std::endl;
    std::cout << "Nothing will be done" << std::endl;
    block = true;
    return;
  }
  outputFile->mkdir(category);

  std::vector<TString> SingleMuon = SingleMuon_2018;
  std::vector<TString> SingleElectron = SingleElectron_2018;
  std::vector<TString> MuonEG = MuonEG_2018;
  std::vector<TString> DataSample;
  std::vector<TString> EmbeddedSample = EmbeddedElMu_2018;
  lumi = 59740;



  commonCuts = "pt_1>15.&&pt_2>15.&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&dr_tt>0.3";
  commonCuts += categoryCuts[category];

  TString TriggerEMu = "(trg_muhigh_elow>0.5&&pt_2>24.0)||(trg_ehigh_mulow>0.5&&pt_1>24.0)";
  TString TriggerSingleE = "trg_singleelectron>0.5&&pt_1>33.";
  TString TriggerSingleMu = "trg_singlemuon>0.5&&pt_2>25.";
  globalWeight = "1.0*";

  if (era=="2017") {
    SingleMuon = SingleMuon_2017;
    SingleElectron = SingleElectron_2017;
    MuonEG = MuonEG_2017;
    EmbeddedSample = EmbeddedElMu_2017;
    TriggerSingleE = "trg_singleelectron>0.5&&pt_1>28.";
    TriggerSingleMu = "trg_singlemuon>0.5&&pt_2>25.";    
    lumi = 41900;
  }
  if (era=="2016") {
    SingleMuon = SingleMuon_2016;
    SingleElectron = SingleElectron_2016;
    MuonEG = MuonEG_2016;
    EmbeddedSample = EmbeddedElMu_2016;
    lumi = 35890;
    TriggerSingleE = "trg_singleelectron>0.5&&pt_1>26.";
    TriggerSingleMu = "trg_singlemuon>0.5&&pt_2>25.";
  }
  if (triggerOption==0) {
    DataSample = MuonEG;
    TriggerCut = "&&(" + TriggerEMu + ")";
    TriggerCutEMu = "&&(" + TriggerEMu + ")";
    globalWeight = "weightEMu*";
  }
  else if (triggerOption==1) {

    for (unsigned int i=0; i<SingleMuon.size(); ++i) 
      DataSample.push_back(SingleMuon.at(i));
    for (unsigned int i=0; i<SingleElectron.size(); ++i) 
      DataSample.push_back(SingleElectron.at(i));

    TriggerCut = "&&(("+TriggerSingleE+")||("+TriggerSingleMu+"))";
    TriggerCutSingleMu = "&&(" + TriggerSingleMu + ")";
    TriggerCutSingleE = "&&(" + TriggerSingleE + ")&&!(" + TriggerSingleMu + ")";
    
    globalWeight = "weightSingle*";

  }
  else  {
    for (unsigned int i=0; i<SingleMuon.size(); ++i)
      DataSample.push_back(SingleMuon.at(i));
    for (unsigned int i=0; i<SingleElectron.size(); ++i)
      DataSample.push_back(SingleElectron.at(i));
    for (unsigned int i=0; i<MuonEG.size(); ++i)
      DataSample.push_back(MuonEG.at(i));

    TriggerCut = "&&(("+TriggerSingleE+")||("+TriggerSingleMu+")||("+TriggerEMu+"))";
    TriggerCutSingleMu = "&&(" + TriggerSingleMu + ")";
    TriggerCutSingleE = "&&(" + TriggerSingleE + ")&&!(" + TriggerSingleMu + ")";
    TriggerCutEMu = "&&("+TriggerEMu+")&&!("+TriggerSingleE+")&&!("+TriggerSingleMu+")";

    globalWeight = "weight*";

  }

  sampleXSecMap = xsecs;

  nameSampleMap["Data"] = DataSample; nameHistoMap["Data"] = "data_obs";
  nameSampleMap["DYJetsToLL"] = DYJetsToLL; nameHistoMap["DYJetsToLL"] = "ZLL"; 
  nameSampleMap["DYJetsToTT"] = DYJetsToTT; nameHistoMap["DYJetsToTT"] = "ZTT";
  nameSampleMap["WJets"] = WJets; nameHistoMap["WJets"] = "W";
  nameSampleMap["TTbar"] = TT; nameHistoMap["TTbar"] = "TTL";
  nameSampleMap["EWK"] = EWK; nameHistoMap["EWK"] = "VVL";
  nameSampleMap["EMB"] = EmbeddedSample; nameHistoMap["EMB"] = "EMB";
  nameSampleMap["ggHWW125"] = GluGluHToWW; nameHistoMap["ggHWW125"] = "ggHWW125";
  nameSampleMap["qqHWW125"] = VBFHToWW; nameHistoMap["qqHWW125"] = "qqHWW125";
  nameSampleMap["WHWW125"] = WHToWW; nameHistoMap["WHWW125"] = "WHWW125";
  nameSampleMap["ZHWW125"] = ZHToWW; nameHistoMap["ZHWW125"] = "ZHWW125";
  nameSampleMap["ggHTT125"] = GluGluHToTauTau; nameHistoMap["ggHTT125"] = "ggH125";
  nameSampleMap["qqHTT125"] = VBFHToTauTau; nameHistoMap["qqHTT125"] = "qqH125";
  nameSampleMap["WHTT125"] = WHToTauTau; nameHistoMap["WHTT125"] = "WH125";
  nameSampleMap["ZHTT125"] = ZHToTauTau; nameHistoMap["ZHTT125"] = "ZH125";
  histNamesMSSM[baseNameBBH] = "bbH";
  histNamesMSSM[baseNameGGH] = "ggH";

  samplesContainer.clear();
  if (sampleToProcess=="Data") {
    InitializeSample("Data");
    if (runOnEmbedded)
      samplesContainer.push_back("EMB");
    else 
      samplesContainer.push_back("DYJetsToTT");
    samplesContainer.push_back("DYJetsToLL");
    samplesContainer.push_back("WJets");
    samplesContainer.push_back("EWK");
    samplesContainer.push_back("TTbar");
  }
  if (sampleToProcess=="TTbar") {
    samplesContainer.push_back("TTbar");
  }
  if (sampleToProcess=="EMB") {
    samplesContainer.push_back("EMB");
  }
  if (sampleToProcess=="DYJets") {
    samplesContainer.push_back("DYJetsToLL");
    if (!runOnEmbedded) 
      samplesContainer.push_back("DYJetsToTT");
  }
  if (sampleToProcess=="EWK") {
    samplesContainer.push_back("WJets");
    samplesContainer.push_back("EWK");
  }
  if (sampleToProcess=="SMHiggs") {
    samplesContainer.push_back("ggHWW125");
    samplesContainer.push_back("qqHWW125");
    samplesContainer.push_back("WHWW125");
    samplesContainer.push_back("ZHWW125");
    samplesContainer.push_back("ggHTT125");
    samplesContainer.push_back("qqHTT125");
    samplesContainer.push_back("WHTT125");
    samplesContainer.push_back("ZHTT125");
  }
  if (sampleToProcess=="MSSMHiggs") {
    CreateMSSMList(baseNameBBH);
    CreateMSSMList(baseNameGGH);
  }

  InitializeSamples();

  numberOfShapeSystematics = CreateShapeSystematicsMap();
  numberOfWeightSystematics = CreateWeightSystematicsMap();
  

}

TString CardsEMu::SampleSpecificCut(TString name, TString sampleName) {
  TString triggerCut = TriggerCut;
  TString mcDiTauVeto("");
  TString ngenPartonsCut("");
  if (sampleName.Contains("MuonEG"))
    triggerCut = TriggerCutEMu;
  else if (sampleName.Contains("SingleElectron")||sampleName.Contains("EGamma"))
    triggerCut = TriggerCutSingleE;
  else if (sampleName.Contains("SingleMuon"))
    triggerCut = TriggerCutSingleMu;
  else
    triggerCut = TriggerCut;

  if (sampleName.Contains("WJetsToLNu_0")||
      sampleName.Contains("DYJetsToLL_M-50_0")||
      sampleName.Contains("DYJetsToTT_M-50_0"))
    ngenPartonsCut = "&&gen_noutgoing==0";

  if (sampleName.Contains("WJetsToLNu_1")||
      sampleName.Contains("DYJetsToLL_M-50_1")||
      sampleName.Contains("DYJetsToTT_M-50_1"))
    ngenPartonsCut = "&&gen_noutgoing==1";

  if (sampleName.Contains("WJetsToLNu_2")||
      sampleName.Contains("DYJetsToLL_M-50_2")||
      sampleName.Contains("DYJetsToTT_M-50_2"))
    ngenPartonsCut = "&&gen_noutgoing==2";

  if (sampleName.Contains("WJetsToLNu_3")||
      sampleName.Contains("DYJetsToLL_M-50_3")||
      sampleName.Contains("DYJetsToTT_M-50_3"))
    ngenPartonsCut = "&&gen_noutgoing==3";

  if (sampleName.Contains("WJetsToLNu_4")||
      sampleName.Contains("DYJetsToLL_M-50_4")||
      sampleName.Contains("DYJetsToTT_M-50_4"))
    ngenPartonsCut = "&&gen_noutgoing==4";

  if ((name=="TTbar"||name=="EWK"||name=="WJets")&&runOnEmbedded)
    mcDiTauVeto = "&&!(gen_match_1==3&&gen_match_2==4)";

  if (name=="DYJetsToLL")
    mcDiTauVeto = "&&!(gen_match_1==3&&gen_match_2==4)";
  if (name=="DYJetsToTT")
    mcDiTauVeto = "&&(gen_match_1==3&&gen_match_2==4)";

  TString Cut = triggerCut+ngenPartonsCut+mcDiTauVeto;

  return Cut;
}

void CardsEMu::InitializeSample(TString name) {

  std::vector<TString> sampleNames = nameSampleMap[name];
  for (auto sampleName : sampleNames) {
    TString baseFileName = sampleName;
    if (name=="DYJetsToLL") 
      baseFileName = DYJetsLLFiles[sampleName];
    if (name=="DYJetsToTT") 
      baseFileName = DYJetsTTFiles[sampleName];
    if (name=="WJets") 
      baseFileName = WJetsFiles[sampleName];
    TString fullPathName = input_dir + "/" + baseFileName + ".root";
    TFile * file = new TFile(fullPathName);
    if (file->IsZombie()) {
      std::cout << "cannot open file " << fullPathName << std::endl;
      std::cout << "nothing will be done " << std::endl;
      block = true;
      return;
    }
    sampleFileMap[sampleName] = file;
    TH1D * histEvents = (TH1D*)file->Get("nWeightedEvents");
    sampleNeventsMap[sampleName] = histEvents->GetSumOfWeights();
    sampleSpecificCutMap[sampleName] = SampleSpecificCut(name,sampleName);
  }
  
  if (name.Contains("DYJetsToLL")) {

    float nIncl = sampleNeventsMap["DYJetsToLL_M-50_0"];
    float xsecIncl = sampleXSecMap["DYJetsToLL_M-50"];

    float n1Jet = sampleNeventsMap["DY1JetsToLL_M-50"];
    float xsec1Jet = sampleXSecMap["DY1JetsToLL_M-50"];

    float n2Jet = sampleNeventsMap["DY2JetsToLL_M-50"];
    float xsec2Jet = sampleXSecMap["DY2JetsToLL_M-50"];

    float n3Jet = sampleNeventsMap["DY3JetsToLL_M-50"];
    float xsec3Jet = sampleXSecMap["DY3JetsToLL_M-50"];

    float n4Jet = sampleNeventsMap["DY4JetsToLL_M-50"];
    float xsec4Jet = sampleXSecMap["DY4JetsToLL_M-50"];

    sampleNormMap["DYJetsToLL_M-50_0"] = lumi*xsecIncl/nIncl;
    sampleNormMap["DYJetsToLL_M-50_1"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["DYJetsToLL_M-50_2"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["DYJetsToLL_M-50_3"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["DYJetsToLL_M-50_4"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

    sampleNormMap["DY1JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["DY2JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["DY3JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["DY4JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

  }
  else if (name.Contains("DYJetsToTT")) {

    float nIncl = sampleNeventsMap["DYJetsToTT_M-50_0"];
    float xsecIncl = sampleXSecMap["DYJetsToLL_M-50"];

    float n1Jet = sampleNeventsMap["DY1JetsToTT_M-50"];
    float xsec1Jet = sampleXSecMap["DY1JetsToLL_M-50"];

    float n2Jet = sampleNeventsMap["DY2JetsToTT_M-50"];
    float xsec2Jet = sampleXSecMap["DY2JetsToLL_M-50"];

    float n3Jet = sampleNeventsMap["DY3JetsToTT_M-50"];
    float xsec3Jet = sampleXSecMap["DY3JetsToLL_M-50"];

    float n4Jet = sampleNeventsMap["DY4JetsToLL_M-50"];
    float xsec4Jet = sampleXSecMap["DY4JetsToLL_M-50"];

    sampleNormMap["DYJetsToTT_M-50_0"] = lumi*xsecIncl/nIncl;
    sampleNormMap["DYJetsToTT_M-50_1"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["DYJetsToTT_M-50_2"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["DYJetsToTT_M-50_3"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["DYJetsToTT_M-50_4"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

    sampleNormMap["DY1JetsToTT_M-50"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["DY2JetsToTT_M-50"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["DY3JetsToTT_M-50"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["DY4JetsToTT_M-50"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

  }
  else if (name.Contains("WJets")) {

    float nIncl = sampleNeventsMap["WJetsToLNu_0"];
    float xsecIncl = sampleXSecMap["WJetsToLNu"];

    float n1Jet = sampleNeventsMap["W1JetsToLNu"];
    float xsec1Jet = sampleXSecMap["W1JetsToLNu"];

    float n2Jet = sampleNeventsMap["W2JetsToLNu"];
    float xsec2Jet = sampleXSecMap["W2JetsToLNu"];

    float n3Jet = sampleNeventsMap["W3JetsToLNu"];
    float xsec3Jet = sampleXSecMap["W3JetsToLNu"];

    float n4Jet = sampleNeventsMap["W4JetsToLNu"];
    float xsec4Jet = sampleXSecMap["W4JetsToLNu"];

    sampleNormMap["WJetsToLNu_0"] = lumi*xsecIncl/nIncl;
    sampleNormMap["WJetsToLNu_1"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["WJetsToLNu_2"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["WJetsToLNu_3"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["WJetsToLNu_4"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

    sampleNormMap["W1JetsToLNu"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    sampleNormMap["W2JetsToLNu"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    sampleNormMap["W3JetsToLNu"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    sampleNormMap["W4JetsToLNu"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

  }
  else {
    for (auto sampleName : sampleNames) {
      if (name=="Data"||name=="EMB")
	sampleNormMap[sampleName] = 1.0;
      else if (name.Contains("SUSY")) 
	sampleNormMap[sampleName] = lumi/sampleNeventsMap[sampleName];
      else
	sampleNormMap[sampleName] = lumi*sampleXSecMap[sampleName]/sampleNeventsMap[sampleName];
    }
  }

}
void CardsEMu::InitializeSamples() {
  for (auto sample : samplesContainer) {
    InitializeSample(sample);
  }
}

void CardsEMu::CreateMSSMList(TString baseName) {

  for (auto mass : masses) {
    TString susyInputFile = input_dir + "/" + baseName + mass + ".root";
    TFile * file = new TFile(susyInputFile);
    if (file->IsZombie()) {
      std::cout << "CardsEMu::CreateMSSMList() -> " << std::endl;
      std::cout << "Cannot open file " << susyInputFile << std::endl;
      std::cout << "skipping mass = " << mass << " for " << baseName << std::endl;
      continue;
    }
    TString sampleName = baseName + mass;
    vector<TString> sampleNames; 
    sampleNames.clear();
    sampleNames.push_back(sampleName);
    nameSampleMap[sampleName] = sampleNames;
    sampleXSecMap[sampleName] = 1.0; // 1 pb
    nameHistoMap[sampleName] = histNamesMSSM[baseName]+mass;
    samplesContainer.push_back(sampleName);
  }
  
}

TH1D * CardsEMu::ProcessSample(TString name,
			       TString sysName,
			       bool OS,
			       bool weightSys) {

  std::cout << std::endl;
  std::cout << "processing " << name;
  if (sysName!="") 
    std::cout << " with systematics " << sysName;
  std::cout << std::endl;
  if (!OS) std::cout << "  same-sign control region" << std::endl; 
  TString suffix = "_SS";
  TString weight = globalWeight;
  TString osCut = "&&os>0.5";
  if (OS) {
    suffix = "_OS";
    osCut= "&&os>0.5";
  }
  else {
    weight += "qcdweight*";
    osCut = "&&os<0.5";
  }
  TString histName = name + suffix + "_" + sysName;
  if (sysName=="") histName = name + suffix;
  std::cout << "   creating histo " << histName << std::endl;
  TH1D * hist = createHisto(histName);
  
  TString treeName = BaseTreeName;
  if (sysName=="") {
    treeName = BaseTreeName;
  }
  else {
    if (weightSys) {
      treeName = BaseTreeName;
      weight += weightSystematicsMap[sysName];
    }
    else {
      treeName = BaseTreeName + "_" + shapeSystematicsMap[sysName];
    }
  }

  vector<TString> sampleNames = nameSampleMap[name];
  for (auto sampleName : sampleNames) {
    TFile * file = sampleFileMap[sampleName];

    //    std::cout << "  weight = " << weight << "  " << treeName << std::endl;
    TTree * tree = (TTree*)file->Get(treeName);
    if (tree==NULL) {
      std::cout << "Tree named " << treeName << " does not exist for sample " << sampleName << std::endl;
      std::cout << "returning null pointer" << std::endl;
      delete hist;
      return NULL;
    }
    std::cout << "      " << sampleName << "  :  entries in tree = " << tree->GetEntries() << std::endl;
    TString cut = commonCuts+osCut+sampleSpecificCutMap[sampleName];
    double norm = sampleNormMap[sampleName];
    TH1D * histSample = createHisto("hist");
    TString Cuts = weight + "(" + cut + ")";
    //    std::cout << "       cuts = " << Cuts << std::endl;
    tree->Draw(variable+">>hist",Cuts);
    hist->Add(hist,histSample,1.,norm);
    delete histSample;
  }

  zeroBins(hist);
  return hist;

}

void CardsEMu::SetVariableToPlot(TString var, int nbins, double xmin, double xmax) {
  if (xmin>xmax) {
    std::cout << "CardsEMu::SetVariableToPlot() -> " << std::endl; 
    std::cout << "warning ->  xmin : " << xmin << " is greater than xmax : " << xmax << std::endl;
    std::cout << "nothing will be done..." << std::endl;
    block = true;
    return;
  }

  variable = var;
  nBins = nbins;
  double binWidth = (xmax-xmin)/double(nBins);
  for (int iB=0; iB<=nBins; ++iB)
    Bins[iB] = xmin + double(iB)*binWidth;

}

void CardsEMu::SetVariableToPlot(TString var, int nbins, double * bins) {

  bool notOrdered = false;
  for (int iB=0; iB<nbins; ++iB) {
    if (bins[iB]>bins[iB+1]) {
      std::cout << iB << " : " << bins[iB] << " vs. " << bins[iB+1] << std::endl;
      notOrdered = true;
      break;
    }
  }

  if (notOrdered) {
    std::cout << "CardsEMu::SetVariableToPlot() -> " << std::endl; 
    std::cout << "warning -> bins are not in increasing order  : " << std::endl;
    std::cout << "nothing will be done..." << std::endl;
    block = true;
    return;
  }

  variable = var;
  nBins = nbins;
  for (int iB=0; iB<=nBins; ++iB)
    Bins[iB] = bins[iB];

}

void CardsEMu::SetVariableToPlot(TString var, int nbins, vector<double> bins) {

  bool notOrdered = false;
  for (int iB=0; iB<nbins; ++iB) {
    if (bins[iB]>bins[iB+1]) {
      notOrdered = true;
      break;
    }
  }

  if (notOrdered) {
    std::cout << "CardsEMu::SetVariableToPlot() -> " << std::endl; 
    std::cout << "warning -> bins are not in increasing order  : " << std::endl;
    std::cout << "nothing will be done..." << std::endl;
    block = true;
    return;
  }

  variable = var;
  nBins = nbins;
  for (int iB=0; iB<=nBins; ++iB)
    Bins[iB] = bins[iB];

}

TH1D * CardsEMu::createHisto(TString histName) {
  TH1D * hist = new TH1D(histName,"",nBins,Bins);
  return hist;
}

int CardsEMu::CreateShapeSystematicsMap() {

  int numberOfShapeSystematics = 0;

  if (sampleToProcess=="Data")
    return numberOfShapeSystematics;

  if (sampleToProcess=="EMB") {
    for (auto mapIter : EmbeddedShapeSystematics) {
      TString sysName = mapIter.first;
      TString treeName = mapIter.second;
      shapeSystematicsMap[sysName+"Up"] = treeName+"Up";
      shapeSystematicsMap[sysName+"Down"] = treeName+"Down";
      numberOfShapeSystematics += 2;
    }       
    return numberOfShapeSystematics;
  }
  
  std::map<TString, TString> yearSysMap = ShapeSystematics_2016;
  if (era=="2016")
    yearSysMap = ShapeSystematics_2016;
  if (era=="2017")
    yearSysMap = ShapeSystematics_2017;
  if (era=="2018")
    yearSysMap = ShapeSystematics_2018;
  
  for (auto mapIter : yearSysMap) {
    TString sysName = mapIter.first;
    TString treeName = mapIter.second;
    shapeSystematicsMap[sysName+"Up"] = treeName+"Up";
    shapeSystematicsMap[sysName+"Down"] = treeName+"Down";
    numberOfShapeSystematics += 2;
  }

  for (auto mapIter : ShapeSystematics_Common) {
    TString sysName = mapIter.first;
    TString treeName = mapIter.second;
    shapeSystematicsMap[sysName+"Up"] = treeName+"Up";
    shapeSystematicsMap[sysName+"Down"] = treeName+"Down";
    numberOfShapeSystematics += 2;
  }

  return numberOfShapeSystematics;

}

int CardsEMu::CreateWeightSystematicsMap() {
   
  int numberOfWeightSystematics = 0; 
  if (sampleToProcess=="EMB") {
    return numberOfWeightSystematics;
  }

  // Data ---->
  if (sampleToProcess=="Data") {
    
    for (auto mapIter : QCDWeightSystematics) {
      TString sysName = mapIter.first;
      TString weightName = mapIter.second;
      weightSystematicsMap[sysName+"_"+era+"Up"] = weightName + "_up*";
      weightSystematicsMap[sysName+"_"+era+"Down"] = weightName + "_down*";
      numberOfWeightSystematics += 2;
    }
  
    std::map<TString,TString> QCDIsoSyst = QCDIsoSystematics_2016;
    if (era=="2017") 
      QCDIsoSyst = QCDIsoSystematics_2017;
    if (era=="2018")
      QCDIsoSyst = QCDIsoSystematics_2018;

    for (auto mapIter : QCDIsoSyst) {
      TString sysName = mapIter.first;
      TString weightName = mapIter.second;
      weightSystematicsMap[sysName] = weightName + "*";
      numberOfWeightSystematics++;
    }

    return numberOfWeightSystematics;

  }
  
  // Prefiring weight;
  if (era!="2018") {
    for (auto mapIter : PrefiringSystematics) {
      TString sysName = mapIter.first;
      TString weightName = mapIter.second;
      weightSystematicsMap[sysName] = weightName+"*";
      numberOfWeightSystematics++;
    }
  }

  // MC ---->
  if (sampleToProcess=="EWK")
    return numberOfWeightSystematics;
  
  // TTbar ---->
  if (sampleToProcess=="TTbar") {
    for (auto mapIter : TopShapeSystematics) {
      TString sysName = mapIter.first;
      TString weightName = mapIter.second;
      weightSystematicsMap[sysName] = weightName+"*";
      numberOfWeightSystematics++;
    }      
    return numberOfWeightSystematics;
  }

  if (sampleToProcess=="DYJets") {
    map<TString,TString> DYSys = DYShapeSystematics;
    if (era=="2016") DYSys = DYShapeSystematics_2016;
    for (auto mapIter : DYSys) {
      TString sysName = mapIter.first;
      TString weightName = mapIter.second;
      weightSystematicsMap[sysName] = weightName+"*";
      numberOfWeightSystematics++;
    }
    return numberOfWeightSystematics;
  }

  return numberOfWeightSystematics;

}

bool CardsEMu::RunData() {

  bool status = true;

  TH1D * histData = ProcessSample("Data","",true,false);
  nameTH1DMap["data_obs"] = histData;
  
  TH1D * QCD = ProcessSample("Data","",false,false);
  for (auto name : samplesContainer) {
    TH1D * hist = ProcessSample(name,"",false,false);
    QCD->Add(QCD,hist,1.,-1.);
    delete hist;
  }
  nameTH1DMap["QCD"] = QCD;

  if (!runWithSystematics) { 
    status = false;
    return status;
  }

  for (auto weightSys : weightSystematicsMap) {
    TString sysName = weightSys.first;
    TH1D * QCD_sys = ProcessSample("Data",sysName,false,true);
    for (auto name : samplesContainer) {
      TH1D * hist = ProcessSample(name,sysName,false,true);
      QCD_sys->Add(QCD_sys,hist,1.,-1.);
      delete hist;
    }
    nameTH1DMap["QCD_"+sysName] = QCD_sys;
  }

  return status;

}

bool CardsEMu::RunModel() {

  bool status = true;

  // central templates ---->
  for (auto name : samplesContainer) {
    TH1D * hist = ProcessSample(name,"",true,false);
    TString sampleHistName = nameHistoMap[name];
    nameTH1DMap[sampleHistName] = hist;
  }

  if (!runWithSystematics) {
    status = false;
    return status;
  }

  // systematics ----->
  for (auto name : samplesContainer) {
    for (auto weightSys : weightSystematicsMap) {
      TString sysName = weightSys.first;
      TString sampleHistName = nameHistoMap[name] + "_" + sysName;
      TH1D * hist = ProcessSample(name,sysName,true,true);
      nameTH1DMap[sampleHistName] = hist;
    }
    for (auto shapeSys : shapeSystematicsMap) {
      TString sysName = shapeSys.first;
      TH1D * hist = ProcessSample(name,sysName,true,false);
      TString sampleHistName = nameHistoMap[name] + "_" + sysName;
      // Systematics may not exist, be careful
      if (hist!=NULL) nameTH1DMap[sampleHistName] = hist;
    }
  }

  return status;

}



void CardsEMu::PrintWeightSystematics() {
  if (block) return;

  std::cout << std::endl;
  std::cout << "CardsEMu::PrintWeightSysMap() " << std::endl;
  std::cout << "Number of weight systematics : " << numberOfWeightSystematics << std::endl;
  for (auto weightSys : weightSystematicsMap) {
    TString sysName = weightSys.first;
    TString treeName = weightSys.second;
    std::cout << sysName << " : " << treeName << std::endl;
  }
  std::cout << std::endl;

}

void CardsEMu::PrintShapeSystematics() {
  if (block) return;

  std::cout << std::endl;
  std::cout << "CardsEMu::PrintShapeSysMap() " << std::endl;
  std::cout << "Number of shape systematics : " << numberOfShapeSystematics << std::endl;
  for (auto shapeSys : shapeSystematicsMap) {
    TString sysName = shapeSys.first;
    TString treeName = shapeSys.second;
    std::cout << sysName << " : " << treeName << std::endl;
  }
  std::cout << std::endl;

}

void CardsEMu::PrintSamples() {
  if (block) return;

  std::cout << std::endl;
  std::cout << "CardsEMu::PrintNameSampleMap() " << std::endl;
  std::cout << "printing sampleName : sampleNorm : sampleSpecificCut" << std::endl;

  if (sampleToProcess=="Data") {
    vector<TString> sampleNames = nameSampleMap["Data"];
    std::cout << "Data -> " << nameHistoMap["Data"] << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    for (auto sampleName : sampleNames) {
      std::cout << "    " << sampleName << " : " << sampleNormMap[sampleName] << " : "<< sampleSpecificCutMap[sampleName] << std::endl;
    }
  }

  for (auto name : samplesContainer) {
    std::cout << name << " -> " << nameHistoMap[name] << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    vector<TString> sampleNames = nameSampleMap[name];
    for (auto sampleName : sampleNames) {
      std::cout << "    " << sampleName << " : " << sampleNormMap[sampleName] << "---" << lumi/sampleNeventsMap[sampleName] << " : " << sampleSpecificCutMap[sampleName]<< std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Cuts to be applied for category " << category << std::endl;
  std::cout << "  -> " << commonCuts << std::endl;

  std::cout << std::endl;

}

bool CardsEMu::Run() {
  if (block) return false;

  bool status = true;

  if (sampleToProcess=="Data")
    status = RunData();
  else 
    status = RunModel();

  // save histograms ---->
  outputFile->cd(category);
  for (auto histogram : nameTH1DMap) {
    TString histoName = histogram.first;
    TH1D * histo = histogram.second;
    histo->Write(histoName);
  }
  return status;
};

bool CardsEMu::CloseFile() {

  bool result = false;
  if (block) return result;

  if (outputFile!=NULL) {
    outputFile->Close();
    delete outputFile;
    result = true;
  }
  return result;

}
void CardsEMu::zeroBins(TH1D * hist) {

  int nbins = hist->GetNbinsX();
  for (int iB=1; iB<=nbins; ++iB) {
    double x = hist->GetBinContent(iB);
    if (x<0.0) {
      hist->SetBinContent(iB,0.);
      hist->SetBinError(iB,0.);
    }
  }

}

CardsEMu::~CardsEMu() {

}

#endif
