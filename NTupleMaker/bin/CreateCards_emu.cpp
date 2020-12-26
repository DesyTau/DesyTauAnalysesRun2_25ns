#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"

bool CombineTriggers;

using namespace std;

TString SpecificCut(TString sample) {
  TString cut("");
  if (CombineTriggers) {
    TString EMu("(trg_muhigh_elow>0.5&&pt_2>24.)||(trg_ehigh_mulow>0.5&&pt_1>24.)");
    if (sample.Contains("SingleMuon")) {
      cut = "&&!(" + EMu + ")";
    }    
    if (sample.Contains("EGamma")) {
      if (sample.Contains("2016"))
	cut = "&&!((trg_singlemuon>0.5&&pt_2>23.)||"+EMu+")";
      if (sample.Contains("2017"))
	cut = "&&!((trg_singlemuon>0.5&&pt_2>25.)||"+EMu+")";
      if (sample.Contains("2018"))
	cut = "&&!((trg_singlemuon>0.5&&pt_2>25.)||"+EMu+")";
    }
  }
  else {
    if (sample.Contains("EGamma")||sample.Contains("SingleElectron")) {
      if (sample.Contains("2016"))
	cut = "&&!(trg_singlemuon>0.5&&pt_2>23.)";
      if (sample.Contains("2017"))
	cut = "&&!(trg_singlemuon>0.5&&pt_2>25.)";
      if (sample.Contains("2018"))
	cut = "&&!(trg_singlemuon>0.5&&pt_2>25.)";
    }
  }
  if (sample.Contains("WJetsToLNu")||sample.Contains("DYJetsToLL_M-50"))
    cut = "&&gen_noutgoing==0";

  return cut;

}

struct SampleAttributes {
  TString name;
  TString histname;
  std::vector<TString> sampleNames;
  std::map<TString, TFile* > fileMap;
  TString weight;
  TString cuts;
};

bool IsSystematicPresent(SampleAttributes sampleAttr, TString SysName) {
  TFile * file = sampleAttr.fileMap[sampleAttr.sampleNames[0]];
  bool isPresent = true;
  TTree * tree = (TTree*)file->Get(BaseTreeName+"_"+SysName);
  if (tree==NULL)
    isPresent = false;
  return isPresent;
}

SampleAttributes SetSampleAttributes(TString name, 
				     TString histname,
				     std::vector<TString> sampleNames,
				     TString weight,
				     TString cuts,
				     TString dir) {
  SampleAttributes sampleAttr;
  sampleAttr.name = name;
  sampleAttr.histname = histname;
  sampleAttr.sampleNames = sampleNames;
  sampleAttr.weight = weight;
  sampleAttr.cuts = cuts;
  std::map<TString,TFile*> fileMap;
  for (auto sample : sampleNames) {
    TFile * file = new TFile(dir+"/"+sample+".root");
    if (file->IsZombie()) {
      std::cout << " cannot open file " << dir << "/" << sample << ".root" << std::endl;
      exit(-1);
    }
    fileMap[sample] = file;
  }
  sampleAttr.fileMap = fileMap;

  return sampleAttr;

}

TH1D * ProcessSample(SampleAttributes sample,
		     TString suffix,
		     TString sysName,
		     TString Variable,		     
		     bool isBTag,
		     double lumi) {

  std::map<TString,TFile*> sampleMap = sample.fileMap;  
  TString name = sample.name;

  int nBins = nBinsNoBTag;
  if (isBTag) nBins = nBinsBTag ;
  float bins[100];
  for (int iB=0; iB<=nBins; ++iB) {
    if (isBTag)
      bins[iB] = binsBTag[iB];
    else 
      bins[iB] = binsNoBTag[iB];
  }

  TString SysName("");
  if (sysName!="")
      SysName = "_" + sysName;

  TH1D * hist = new TH1D(sample.histname+suffix,"",nBins,bins);
  hist->Sumw2();

  std::cout << "Processing sample : " << sample.name << std::endl;

  for (auto mapIter : sampleMap) {
    TString sampleName = mapIter.first; 
    TFile * file = mapIter.second;
    TTree * tree = (TTree*)file->Get(BaseTreeName+SysName);
    TH1D * histWeightsH = (TH1D*)file->Get("nWeightedEvents");
    TString histName = sampleName + "_hist";
    TH1D * histSample = new TH1D(histName,"",nBins,bins);
    histSample->Sumw2();
    TString CutsSample = sample.cuts + SpecificCut(sampleName);
    tree->Draw(Variable+">>"+histName,sample.weight+"("+CutsSample+")");
    double norm = 1.0;
    if (name=="Data"||name=="EmbedZTT") {
      norm = 1.;
    }
    else {
      if (name.Contains("SUSY")) {
	double nevents = histWeightsH->GetSumOfWeights();
	norm = lumi/nevents;
      }
      else {
	double xsec = xsecs[sampleName];
	double nevents = histWeightsH->GetSumOfWeights();
	norm = xsec*lumi/nevents;
      }
    }
    hist->Add(hist,histSample,1.,norm);
    delete histSample;
  }
  return hist;
}

int main(int argc, char * argv[]) {
  
  // argument - config file
  if (argc!=5) {
    std::cout << "usage of the scripts : CreateCards_emu [config] [era=2016,2017,2018] [btag] [category]" << std::endl;
    exit(-1);
  }

  // **** configuration
  Config cfg(argv[1]);
  
  const bool embedded = cfg.get<bool>("IsEmbedded");
  const int triggerOption = cfg.get<int>("TriggerOption");
  const string variable = cfg.get<string>("Variable");
  const vector<string> pzetaRanges = cfg.get<vector<string> >("PZetaRanges");
  
  TString era(argv[2]);
  TString BTAG(argv[3]);
  TString CAT(argv[4]);

  if (era!="2016"&&era!="2017"&&era!="2018") {
    std::cout << "Unknown era : " << era << std::endl;
    exit(-1);
  }

  if (BTAG!="btag"&&BTAG!="nobtag") {
    std::cout << "Unknown btag category : " << BTAG << std::endl;
    exit(-1);
  }
  if (CAT!="ttbar"&&CAT!="low_pzeta"&&CAT!="medium_pzeta"&&CAT!="high_pzeta"&&CAT!="inclusive") {
    std::cout << "Unknown category : " << CAT << std::endl;
    exit(-1);
  }

  bool isBTag = false;
  if (BTAG=="btag")
    isBTag = true;

  //  std::cout << "pzeta ranges : " << std::endl;
  //  for (auto pzetaRange : pzetaRanges) 
  //    std::cout << pzetaRange << std::endl;

  std::map<TString,TString> pzetaCut;
  pzetaCut["ttbar"] = "&&puppipzeta<" + TString(pzetaRanges[0]);
  pzetaCut["low_pzeta"] = "&&puppipzeta>=" + TString(pzetaRanges[0]) + "&&puppipzeta<" + TString(pzetaRanges[1]);
  pzetaCut["medium_pzeta"] = "&&puppipzeta>=" + TString(pzetaRanges[1]) + "&&puppipzeta<" + TString(pzetaRanges[2]);
  pzetaCut["high_pzeta"] = "&&puppipzeta>=" + TString(pzetaRanges[2]);
  pzetaCut["inclusive"] = "&&puppipzeta>=" + TString(pzetaRanges[0]);

  for (auto pzeta_cut : pzetaCut) 
    std::cout << pzeta_cut.first << " : " << pzeta_cut.second << std::endl;

  TString OutputFileName = CAT;
  if (CAT!="ttbar") {
    if (isBTag) 
      OutputFileName += "_btag";
    else 
      OutputFileName += "_nobtag";
  }
  OutputFileName += "_" + era + ".root";

  TString Variable(variable);

  TString dir = "/nfs/dust/cms/user/rasp/grid-jobs/emu_MSSM/"+era;

  std::cout << "Base directory : " << dir << std::endl;

  TString Weight("puweight*mcweight*prefiringweight*");
  TString WeightEmb("mcweight*");
  TString WeightEff("effweight*");

  if (triggerOption==1)
    WeightEff = "effweightExcl*";
  if (triggerOption==2)
    WeightEff = "(effweight/trigweight)*";
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

  TString CutsSingle = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>33.))";
  if (era=="2017")
    CutsSingle = "((trg_singlemuon>0.5&&pt_2>25.)||(trg_singleelectron>0.5&&pt_1>28.))";
  if (era=="2016")
    CutsSingle = "((trg_singlemuon>0.5&&pt_2>23.)||(trg_singleelectron>0.5&&pt_1>26.))";
  
  TString Cuts = CutsSingle;
  CombineTriggers = false;
  if (triggerOption==1) 
    Cuts = CutsEMu;
  if (triggerOption==2) {
    Cuts = "(" + CutsSingle + "||" + CutsEMu + ")";
    CombineTriggers = true;
  }

  Cuts += TString("&&iso_1<0.15&&iso_2<0.20&&extraelec_veto<0.5&&extramuon_veto<0.5&&dr_tt>0.3&&pt_1>15.0&&pt_2>15.");

  if (CAT!="ttbar") {
    if (isBTag) 
      Cuts += "&&nbtag>=1&&njets<=1";
    else 
      Cuts += "&&!(nbtag>=1&&njets<=1)";
  }

  TString CutsCategory = pzetaCut[CAT];

  Cuts += CutsCategory;

  TString CutsOS = Cuts + TString("&&os>0.5");
  TString CutsSS = Cuts + TString("&&os<0.5");
  TString CutsZTT_OS  = CutsOS + TString("&&gen_match_1==3&&gen_match_2==4");
  TString CutsZLL_OS  = CutsOS + TString("&&!(gen_match_1==3&&gen_match_2==4)");
  TString CutsZTT_SS  = CutsSS + TString("&&gen_match_1==3&&gen_match_2==4");
  TString CutsZLL_SS  = CutsSS + TString("&&!(gen_match_1==3&&gen_match_2==4)");

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

  double lumi = 59740;
  if (era=="2017")
    lumi = 41900;
  if (era=="2016")
    lumi = 35890;

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
				      
  if (triggerOption==0) {

    for (unsigned int i=0; i<SingleElectron.size(); ++i)
      DataSamples.push_back(SingleElectron.at(i));

    for (unsigned int i=0; i<SingleMuon.size(); ++i)
      DataSamples.push_back(SingleMuon.at(i));

  }
  if (triggerOption==1) {
    DataSamples = MuonEG;
  }
  if (triggerOption==2) {

    for (unsigned int i=0; i<SingleElectron.size(); ++i)
      DataSamples.push_back(SingleElectron.at(i));

    for (unsigned int i=0; i<SingleMuon.size(); ++i)
      DataSamples.push_back(SingleMuon.at(i));

    for (unsigned int i=0; i<MuonEG.size(); ++i)
      DataSamples.push_back(MuonEG.at(i));

  }

  std::cout << "Data samples : " << std::endl;
  for (unsigned int i=0; i<DataSamples.size(); ++i)
    std::cout << DataSamples.at(i) << std::endl;

  // data
  SampleAttributes DataAttr = SetSampleAttributes("Data","data_obs",DataSamples,"1.0*",CutsOS,dir);

  // ZTT
  SampleAttributes ZttAttr;
  if (embedded) 
    ZttAttr = SetSampleAttributes("EmbedZTT","EmbedZTT",EmbedSamples,WeightEmb,CutsOS,dir);
  else 
    ZttAttr = SetSampleAttributes("ZTT","ZTT",DYSamples,WeightDY,CutsZTT_OS,dir);

  // ZLL
  SampleAttributes ZllAttr = SetSampleAttributes("ZLL","ZLL",DYSamples,WeightDY,CutsZLL_OS,dir);

  // WJets
  SampleAttributes WJetsAttr;
  if (embedded) 
    WJetsAttr = SetSampleAttributes("WJets","W",WJetsSamples,Weight,CutsZLL_OS,dir);
  else
    WJetsAttr = SetSampleAttributes("WJets","W",WJetsSamples,Weight,CutsOS,dir);

  // VV
  SampleAttributes VVAttr;
  if (embedded) 
    VVAttr = SetSampleAttributes("VV","VV",EWKSamples,Weight,CutsZLL_OS,dir);
  else 
    VVAttr = SetSampleAttributes("VV","VV",EWKSamples,Weight,CutsOS,dir);

  // TT
  SampleAttributes TTAttr;
  if (embedded) 
    TTAttr = SetSampleAttributes("TT","TT",TTSamples,WeightTop,CutsZLL_OS,dir);
  else 
    VVAttr = SetSampleAttributes("TT","TT",TTSamples,WeightTop,CutsOS,dir);

  std::vector<SampleAttributes> BkgSamples;
  BkgSamples.push_back(ZttAttr);
  BkgSamples.push_back(ZllAttr);
  BkgSamples.push_back(WJetsAttr);
  BkgSamples.push_back(VVAttr);
  BkgSamples.push_back(TTAttr);

  std::vector<SampleAttributes> GGH_Samples;
  std::vector<SampleAttributes> BBH_Samples;
  for (auto mass: masses) {
    // ggH
    TString ggh_histname = "ggH"+mass;
    TString ggh_name = "SUSYGluGluToHToTauTau_M-"+mass;
    std::vector<TString> ggh_sample_names;
    ggh_sample_names.push_back("SUSYGluGluToHToTauTau_M-"+mass);
    SampleAttributes ggh_sample = SetSampleAttributes(ggh_name,ggh_histname,ggh_sample_names,Weight,CutsOS,dir);
    GGH_Samples.push_back(ggh_sample);

    // bbH
    TString bbh_histname = "bbH"+mass;
    TString bbh_name = "SUSYGluGluToBBHToTauTau_M-"+mass;
    std::vector<TString> bbh_sample_names;
    bbh_sample_names.push_back("SUSYGluGluToBBHToTauTau_M-"+mass);
    SampleAttributes bbh_sample = SetSampleAttributes(bbh_name,bbh_histname,bbh_sample_names,Weight,CutsOS,dir);
    BBH_Samples.push_back(bbh_sample);

  }

  std::cout << "Processing central templates ----> " << std::endl;

  TFile * outputFile = new TFile(OutputFileName,"recreate");
  outputFile->cd("");

  TH1D * dataHist = ProcessSample(DataAttr,"","",Variable,isBTag,lumi);
  outputFile->cd("");
  dataHist->Write("data_obs");
  double dataYield = dataHist->GetSumOfWeights();

  // creating bkgd and data
  double bkgdYield = 0;
  std::map<TString,double> bkg_yields;
  for (auto attr : BkgSamples) {
    TH1D * hist = ProcessSample(attr,"","",Variable,isBTag,lumi);
    outputFile->cd("");
    hist->Write(attr.histname);
    bkgdYield += hist->GetSumOfWeights();
    bkg_yields[attr.histname] = hist->GetSumOfWeights();
  }

  // creating GGH 
  std::map<TString,double> ggh_yields;
  for (auto attr : GGH_Samples) {
    TH1D * hist = ProcessSample(attr,"","",Variable,isBTag,lumi);
    outputFile->cd("");
    hist->Write(attr.histname);
    ggh_yields[attr.histname] = hist->GetSumOfWeights();
  }

  // creating BBH
  std::map<TString,double> bbh_yields;
  for (auto attr : BBH_Samples) {
    TH1D * hist = ProcessSample(attr,"","",Variable,isBTag,lumi);
    outputFile->cd("");
    hist->Write(attr.histname);    
    bbh_yields[attr.histname] = hist->GetSumOfWeights();
  }

  std::cout << std::endl;
  std::cout << "Creating QCD template -----> " << std::endl;
  std::cout << std::endl;

  // *******************************
  // creating QCD
  // *******************************
  DataAttr.cuts = CutsSS;
  DataAttr.weight = WeightQCD;

  // Ztt
  if (embedded) {
    ZttAttr.cuts = CutsSS;
    ZttAttr.weight = WeightEmbSS;
    WJetsAttr.cuts = CutsZLL_SS;
    VVAttr.cuts = CutsZLL_SS;
    TTAttr.cuts = CutsZLL_SS;
  }
  else {
    ZttAttr.cuts = CutsZTT_SS;
    ZttAttr.weight = WeightDYSS;
    WJetsAttr.cuts = CutsSS;
    VVAttr.cuts = CutsSS;
    TTAttr.cuts = CutsSS;    
  }

  ZllAttr.cuts   = CutsZLL_SS;
  ZllAttr.weight = WeightDYSS;
  WJetsAttr.weight = WeightSS;
  VVAttr.weight = WeightSS;
  TTAttr.weight = WeightTopSS;

  BkgSamples.clear();
  BkgSamples.push_back(ZttAttr);
  BkgSamples.push_back(ZllAttr);
  BkgSamples.push_back(WJetsAttr);
  BkgSamples.push_back(VVAttr);
  BkgSamples.push_back(TTAttr);
  
  TH1D * dataSS = ProcessSample(DataAttr,"_ss","",Variable,isBTag,lumi);

  for (auto attr : BkgSamples) {
    //    std::cout << attr.name << "  :  " << attr.weight << " * " << attr.cuts << std::endl;
    TH1D * histSS = ProcessSample(attr,"_ss","",Variable,isBTag,lumi);
    dataSS->Add(dataSS,histSS,1.,-1.);
  }
  bkg_yields["QCD"] = dataSS->GetSumOfWeights();
  bkgdYield += dataSS->GetSumOfWeights();

  std::cout << std::endl;
  for (auto bkg_yield : bkg_yields) {
    std::cout << bkg_yield.first << " : " << bkg_yield.second << std::endl;
  }
  std::cout << "----------------------" << std::endl;
  std::cout << "Total bkg  = " << bkgdYield << std::endl;
  std::cout << "Total data = " << dataYield << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  for (auto ggh_yield : ggh_yields) {
    std::cout << ggh_yield.first << " : " << ggh_yield.second << std::endl;
  }

  std::cout << std::endl;
  for (auto bbh_yield : bbh_yields) {
    std::cout << bbh_yield.first << " : " << bbh_yield.second << std::endl;
  }

  outputFile->cd("");
  dataSS->Write("QCD");
  outputFile->Close();
  delete outputFile;

}
