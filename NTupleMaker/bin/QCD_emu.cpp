#include "DesyTauAnalyses/NTupleMaker/interface/settings.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"

RooWorkspace * w;

struct SampleAttributes {
  TString name;
  TString histname;
  std::vector<TString> sampleNames;
  std::map<TString, TString> baseFileNamesMap;
  std::map<TString, TFile* > fileMap;
  std::map<TString, float> normMap;
  int selectDiTau;
};

bool is2016;
bool is2017;
bool is2018;

bool applySF;

float ptTriggerEmbed2017 = 40;
float etaTriggerEmbed2017 = 1.479;

float dzFilterEff_mc;
float dzFilterEff_data;

float pt_SingleE;
float pt_SingleMu;
float lumi;

int nBinsDR = 6;
float xminDR = 0.3;
float xmaxDR = 6.3;

int nBinsPt = 8;
float binsPt[9] = {15,20,25,30,40,50,70,100,200}; 

int nBinsMTtot = 13;
float binsMTtot[14] = {20,40,60,80,100,125,150,200,250,300,500,700,1000,4000};

int nBinsPt_2D = 4;
float binsPt_2D[5] = {15,25,40,70,200};

int nBinsMvis = 60;
float xminMvis = 0;
float xmaxMvis = 300;

int nBinsPZeta = 30;
float xminPZeta = -150;
float xmaxPZeta = 150;

int triggerOption;
TString input_dir;
TString output_dir;


SampleAttributes SetSampleAttributes(TString name, 
				     TString histname,
				     std::vector<TString> sampleNames,
				     int selectDiTau) {
  SampleAttributes sampleAttr;
  sampleAttr.name = name;
  sampleAttr.histname = histname;
  sampleAttr.sampleNames = sampleNames;
  std::map<TString,TFile*> fileMap;
  std::map<TString,TString> baseFileNamesMap;
  std::map<TString,float> normMap;
  std::map<TString,float> neventsMap;
  for (auto sample : sampleNames) {
    TString baseFileName = sample;
    if (name.Contains("WJets"))
      baseFileName = WJetsFiles[sample];
    if (name.Contains("DYJets"))
      baseFileName = DYJetsFiles[sample]; 
    baseFileNamesMap[sample] = baseFileName;
    TFile * file = new TFile(input_dir+"/"+baseFileName+".root");
    if (file->IsZombie()) {
      std::cout << " cannot open file " << input_dir << "/" << sample << ".root" << std::endl;
      exit(-1);
    }
    TH1D * histEvents = (TH1D*)file->Get("nWeightedEvents");
    neventsMap[sample] = histEvents->GetSumOfWeights();
    fileMap[sample] = file;
  }
  sampleAttr.fileMap = fileMap;
  sampleAttr.baseFileNamesMap = baseFileNamesMap;
  sampleAttr.selectDiTau = selectDiTau;

  if (name.Contains("DYJets")) {

    float nIncl = neventsMap["DYJetsToLL_M-50_0"];
    float xsecIncl = xsecs["DYJetsToLL_M-50"];

    float n1Jet = neventsMap["DY1JetsToLL_M-50"];
    float xsec1Jet = xsecs["DY1JetsToLL_M-50"];

    float n2Jet = neventsMap["DY2JetsToLL_M-50"];
    float xsec2Jet = xsecs["DY2JetsToLL_M-50"];

    float n3Jet = neventsMap["DY3JetsToLL_M-50"];
    float xsec3Jet = xsecs["DY3JetsToLL_M-50"];

    float n4Jet = neventsMap["DY4JetsToLL_M-50"];
    float xsec4Jet = xsecs["DY4JetsToLL_M-50"];

    normMap["DYJetsToLL_M-50_0"] = lumi*xsecIncl/nIncl;
    normMap["DYJetsToLL_M-50_1"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    normMap["DYJetsToLL_M-50_2"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    normMap["DYJetsToLL_M-50_3"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    normMap["DYJetsToLL_M-50_4"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

    normMap["DY1JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    normMap["DY2JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    normMap["DY3JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    normMap["DY4JetsToLL_M-50"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

  }
  else if (name.Contains("WJets")) {

    float nIncl = neventsMap["WJetsToLNu_0"];
    float xsecIncl = xsecs["WJetsToLNu"];

    float n1Jet = neventsMap["W1JetsToLNu"];
    float xsec1Jet = xsecs["W1JetsToLNu"];

    float n2Jet = neventsMap["W2JetsToLNu"];
    float xsec2Jet = xsecs["W2JetsToLNu"];

    float n3Jet = neventsMap["W3JetsToLNu"];
    float xsec3Jet = xsecs["W3JetsToLNu"];

    float n4Jet = neventsMap["W4JetsToLNu"];
    float xsec4Jet = xsecs["W4JetsToLNu"];

    normMap["WJetsToLNu_0"] = lumi*xsecIncl/nIncl;
    normMap["WJetsToLNu_1"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    normMap["WJetsToLNu_2"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    normMap["WJetsToLNu_3"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    normMap["WJetsToLNu_4"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);
    normMap["W1JetsToLNu"] = lumi/(nIncl/xsecIncl+n1Jet/xsec1Jet);
    normMap["W2JetsToLNu"] = lumi/(nIncl/xsecIncl+n2Jet/xsec2Jet);
    normMap["W3JetsToLNu"] = lumi/(nIncl/xsecIncl+n3Jet/xsec3Jet);
    normMap["W4JetsToLNu"] = lumi/(nIncl/xsecIncl+n4Jet/xsec4Jet);

  }
  else {
    for (auto sample : sampleNames) {
      if (name=="Data"||name=="Embedded")
	normMap[sample] = 1.0;
      else {
	normMap[sample] = lumi*xsecs[sample]/neventsMap[sample];
      }
    }
  }
  sampleAttr.normMap = normMap;
  return sampleAttr;

}

std::vector<TString> th1_basename = {

  "deltaR_0jet_OS",
  "deltaR_0jet_SS",

  "deltaR_1jet_OS",
  "deltaR_1jet_SS",

  "deltaR_2jet_OS",
  "deltaR_2jet_SS",  

  "deltaR_inclusive_OS",
  "deltaR_btag_OS",
  "deltaR_nobtag_OS",

  "deltaR_inclusive_SS",
  "deltaR_btag_SS",
  "deltaR_nobtag_SS",

  "pt1_inclusive_OS",
  "pt1_btag_OS",
  "pt1_nobtag_OS",

  "pt1_inclusive_SS",
  "pt1_btag_SS",
  "pt1_nobtag_SS",

  "pt2_inclusive_OS",
  "pt2_btag_OS",
  "pt2_nobtag_OS",

  "pt2_inclusive_SS",
  "pt2_btag_SS",
  "pt2_nobtag_SS",

  "deltaR_inclusive_isoMu_antiisoE_OS",
  "pt1_inclusive_isoMu_antiisoE_OS",
  "pt2_inclusive_isoMu_antiisoE_OS",

  "deltaR_inclusive_isoMu_antiisoE_SS",
  "pt1_inclusive_isoMu_antiisoE_SS",
  "pt2_inclusive_isoMu_antiisoE_SS",

  "deltaR_inclusive_antiisoMu_antiisoE_OS",
  "pt1_inclusive_antiisoMu_antiisoE_OS",
  "pt2_inclusive_antiisoMu_antiisoE_OS",

  "deltaR_inclusive_antiisoMu_antiisoE_SS",
  "pt1_inclusive_antiisoMu_antiisoE_SS",
  "pt2_inclusive_antiisoMu_antiisoE_SS",

  "mTtot_inclusive_OS",
  "mTtot_btag_OS",
  "mTtot_nobtag_OS",
  "mTtot_ttbar_OS",

  "mTtot_inclusive_SS",
  "mTtot_btag_SS",
  "mTtot_nobtag_SS",
  "mTtot_ttbar_SS",

  "mTtot_inclusive_looseisoMu_SS",
  "mTtot_btag_looseisoMu_SS",
  "mTtot_nobtag_looseisoMu_SS",
  "mTtot_ttbar_looseisoMu_SS",

  "mTtot_inclusive_isoMu_SS",
  "mTtot_btag_isoMu_SS",
  "mTtot_nobtag_isoMu_SS",
  "mTtot_ttbar_isoMu_SS",

  "mvis_sel_signal_SS",
  "mvis_sel_signal_OS",

  "mvis_sel_ttbar_SS",
  "mvis_sel_ttbar_OS",

  "mTtot_sel_signal_SS",
  "mTtot_sel_signal_OS",

  "mTtot_sel_ttbar_SS",
  "mTtot_sel_ttbar_OS",

  "pt1_sel_signal_SS",
  "pt1_sel_signal_OS",

  "pt1_sel_ttbar_SS",
  "pt1_sel_ttbar_OS",

  "pt2_sel_signal_SS",
  "pt2_sel_signal_OS",

  "pt2_sel_ttbar_SS",
  "pt2_sel_ttbar_OS",

  "pzeta_sel_SS",
  "pzeta_sel_OS"

};

std::vector<TString> th2_basename = {
  "pt1_pt2_inclusive_OS",
  "pt1_pt2_btag_OS",
  "pt1_pt2_nobtag_OS",

  "pt1_pt2_inclusive_SS",
  "pt1_pt2_btag_SS",
  "pt1_pt2_nobtag_SS",
   
  "pt1_pt2_inclusive_isoMu_antiisoE_OS",
  "pt1_pt2_inclusive_isoMu_antiisoE_SS",
 
  "pt1_pt2_inclusive_antiisoMu_antiisoE_OS",
  "pt1_pt2_inclusive_antiisoMu_antiisoE_SS",

};

void CreateMaps(std::map<TString,TH1D*> & th1, 
		std::map<TString,TH2D*> & th2,
		TString sampleName) {

  for (auto name : th1_basename) {    
    TH1D * hist;
    if (name.Contains("deltaR"))
      hist = new TH1D(name+"_"+sampleName,"",nBinsDR,xminDR,xmaxDR);
    else if (name.Contains("mTtot"))
      hist = new TH1D(name+"_"+sampleName,"",nBinsMTtot,binsMTtot);
    else if (name.Contains("mvis"))
      hist = new TH1D(name+"_"+sampleName,"",nBinsMvis,xminMvis,xmaxMvis);
    else if (name.Contains("pzeta"))
      hist = new TH1D(name+"_"+sampleName,"",nBinsPZeta,xminPZeta,xmaxPZeta);
    else
      hist = new TH1D(name+"_"+sampleName,"",nBinsPt,binsPt);
    hist->Sumw2();
    th1[name] = hist;
  }
  for (auto name : th2_basename) {
    TH2D * hist = new TH2D(name+"_"+sampleName,"",nBinsPt_2D,binsPt_2D,nBinsPt_2D,binsPt_2D);
    hist->Sumw2();
    th2[name] = hist;
  }

}

void RunOnTree(TTree * tree,
	       TString name,
	       TString sampleName,
	       std::map<TString,TH1D*> & th1, 
	       std::map<TString,TH2D*> & th2,
	       int selectDiTau,
	       float norm) {

  Synch17Tree * tree17 =  new Synch17Tree(tree, true);
  Long64_t nentries = tree17->GetEntries();
  /*
  std::cout << "Number of events in tree = " << nentries << std::endl;

  for (auto mapIter : th1) 
    std::cout << mapIter.first << ":" << mapIter.second << ":" << mapIter.second->GetName() << std::endl;
  for (auto mapIter : th2) 
    std::cout << mapIter.first << ":" << mapIter.second << ":" << mapIter.second->GetName() << std::endl;
  */

  bool isData = name=="Data";
  bool isEmbedded = name=="Embedded";
  bool isTTbar = name=="TT";
  bool isDY = name=="DYJetsLL";

  int count = 0;


  for (Long64_t ie = 0; ie<nentries; ++ie) {
    tree17->GetEntry(ie);
    if (ie!=0&&ie%100000==0)
      std::cout << "    processed " << ie << " events" << std::endl;

    /*
    std::cout << "pt_1 = " << pt_1 << std::endl;
    std::cout << "pt_2 = " << pt_2 << std::endl;
    std::cout << "iso_1 = " << iso_1 << std::endl;
    std::cout << "iso_2 = " << iso_2 << std::endl;
    std::cout << "extraelec_veto = " << extraelec_veto << std::endl;
    std::cout << "extramuon_veto = " << extramuon_veto << std::endl;
    std::cout << "pt_SingleE = " << pt_SingleE << std::endl;
    std::cout << "pt_SingleMu = " << pt_SingleMu << std::endl;
    std::cout << std::endl;
    */

    if (tree17->pt_1<15.0) continue;
    if (tree17->pt_2<15.0) continue;
    if (tree17->iso_1>0.5) continue;
    if (tree17->iso_2>0.5) continue;
    if (tree17->extraelec_veto>0.5) continue;
    if (tree17->extramuon_veto>0.5) continue;

    bool isSingleE  = (tree17->pt_1 > pt_SingleE  && tree17->trg_singleelectron>0.5);
    bool isSingleMu = (tree17->pt_2 > pt_SingleMu && tree17->trg_singlemuon>0.5);
    bool isSingleLep = isSingleE || isSingleMu;
    bool isEMu = (tree17->pt_1>24.&&tree17->trg_ehigh_mulow) || (tree17->pt_2>24.&&tree17->trg_muhigh_elow);
    if (!isData) {
      if (triggerOption==0&&!isEMu) continue;
      if (triggerOption==1&&!isSingleLep) continue;
      if (triggerOption==2&&!(isEMu||isSingleLep)) continue;
    }
    else {
      bool trigger = true;
      if (sampleName.Contains("EGamma")||sampleName.Contains("SingleElectron")) {
	trigger = isSingleE && !isSingleMu;
	if (!trigger) continue;
      }
      if (sampleName.Contains("SingleMuon")) {
	trigger = isSingleMu;
	if (!trigger) continue;
      }
      if (sampleName.Contains("MuonEG"))  {
	trigger = isEMu;
	if (triggerOption==2)
	  trigger = isEMu && !isSingleE && !isSingleMu;
	if (!trigger) continue;
      }

    }
    if (sampleName.Contains("WJetsToLNu_0")||sampleName.Contains("DYJetsToLL_M-50_0"))
      if (tree17->gen_noutgoing!=0) continue;
    if (sampleName.Contains("WJetsToLNu_1")||sampleName.Contains("DYJetsToLL_M-50_1"))
      if (tree17->gen_noutgoing!=1) continue;
    if (sampleName.Contains("WJetsToLNu_2")||sampleName.Contains("DYJetsToLL_M-50_2"))
      if (tree17->gen_noutgoing!=2) continue;
    if (sampleName.Contains("WJetsToLNu_3")||sampleName.Contains("DYJetsToLL_M-50_3"))
      if (tree17->gen_noutgoing!=3) continue;
    if (sampleName.Contains("WJetsToLNu_4")||sampleName.Contains("DYJetsToLL_M-50_4"))
      if (tree17->gen_noutgoing!=4) continue;

    //    std::cout << "selectDiTau = " << selectDiTau << std::endl;
    bool isDiTau = tree17->gen_match_1==3&&tree17->gen_match_2==4;
    if (selectDiTau==1)
      if (!isDiTau) continue;
    if (selectDiTau==-1)
      if (isDiTau) continue;
    //    std::cout << "Cuts passed " << std::endl;

    // weight to fill histo
    
    float effweight = tree17->effweight;    
    float weight = 1.0;    
    if (triggerOption==0)
      effweight = tree17->effweightEMu;
    if (triggerOption==1)
      effweight = tree17->effweightSingle;
    
    float effweightEMu = 1.0;//tree17->effweightEMu;
    float effweightSingle = 1.0;//tree17->effweightSingle;
    float effweightComb = 1.0;//tree17->effweight;
    float trigweightSingle = 1.0;//tree17->trigweightSingle;
    float trigweightEMu = 1.0;//tree17->trigweightEMu;
    float trigweight_comb = 1.0;//tree17->trigweight;

    if (!isData&&applySF) {
      TString suffix = "mc";
      TString suffixRatio = "ratio";
      if (isEmbedded) {suffix = "embed"; suffixRatio = "embed_ratio";}

      w->var("m_pt")->setVal(tree17->pt_2);
      w->var("m_eta")->setVal(tree17->eta_2);
      w->var("m_iso")->setVal(tree17->iso_2);

      float eff_data_trig_m = w->function("m_trg_ic_data")->getVal();
      float eff_mc_trig_m = w->function("m_trg_ic_" + suffix)->getVal();

      float eff_data_trig_mhigh = w->function("m_trg_23_ic_data")->getVal();
      float eff_data_trig_mlow = w->function("m_trg_8_ic_data")->getVal();
      float eff_mc_trig_mhigh = w->function("m_trg_23_ic_"+suffix)->getVal();
      float eff_mc_trig_mlow = w->function("m_trg_8_ic_"+suffix)->getVal();

      float idisoweight_2 = w->function("m_idiso_ic_" + suffixRatio)->getVal();
      float trkeffweight_2 = w->function("m_trk_ratio")->getVal(); //  may be wrong

      // electron weights
      w->var("e_pt")->setVal(tree17->pt_1);
      w->var("e_eta")->setVal(tree17->eta_1);
      w->var("e_iso")->setVal(tree17->iso_1);
      
      float eff_data_trig_e = w->function("e_trg_ic_data")->getVal();
      float eff_mc_trig_e = w->function("e_trg_ic_" + suffix)->getVal();

      float eff_data_trig_ehigh = w->function("e_trg_23_ic_data")->getVal();
      float eff_data_trig_elow = w->function("e_trg_12_ic_data")->getVal();
      float eff_mc_trig_ehigh = w->function("e_trg_23_ic_"+suffix)->getVal();
      float eff_mc_trig_elow = w->function("e_trg_12_ic_"+suffix)->getVal();

      float idisoweight_1 = w->function("e_idiso_ic_" + suffixRatio)->getVal();
      float trkeffweight_1 = w->function("e_trk_" + suffixRatio)->getVal();

      if (is2017&&isEmbedded) {
	if (tree17->pt_1<ptTriggerEmbed2017&&fabs(tree17->eta_1)>etaTriggerEmbed2017) {
	  eff_mc_trig_e = 1.0;
	}
      }

      if (tree17->pt_1<24.) {
	eff_data_trig_ehigh = 0;
	eff_mc_trig_ehigh = 0;
      }
	
      if (tree17->pt_2<24.) {
	eff_data_trig_mhigh = 0;
	eff_mc_trig_mhigh = 0;
      }

      if (tree17->pt_1<pt_SingleE) {
	eff_data_trig_e = 0;
	eff_mc_trig_e = 0;
      }

      if (tree17->pt_2<pt_SingleMu) {
	eff_data_trig_m = 0;
	eff_mc_trig_m = 0;
      }

      float eff_single_data = eff_data_trig_e + eff_data_trig_m - eff_data_trig_e*eff_data_trig_m;
      float eff_single_mc   = eff_mc_trig_e + eff_mc_trig_m - eff_mc_trig_e*eff_mc_trig_m;


      if (eff_single_mc<1e-3||eff_single_data<1e-3) 
	trigweightSingle = 0.0;
      else
	trigweightSingle = eff_single_data/eff_single_mc;
      
      
      float eff_emu_data = 
			eff_data_trig_mhigh*eff_data_trig_elow + 
			eff_data_trig_mlow*eff_data_trig_ehigh -
			eff_data_trig_mhigh*eff_data_trig_ehigh;
      float eff_emu_mc = 
			eff_mc_trig_mhigh*eff_mc_trig_elow + 
			eff_mc_trig_mlow*eff_mc_trig_ehigh -
			eff_mc_trig_mhigh*eff_mc_trig_ehigh;	
      
      if (eff_emu_mc<1e-3||eff_emu_data<1e-3)
	trigweightEMu = 0.0;
      else
	trigweightEMu = eff_emu_data/eff_emu_mc;
      
      float eff_comb_data = 
		     eff_data_trig_e + 
		     eff_data_trig_m + 
		     eff_emu_data*dzFilterEff_data - 
		     eff_data_trig_e*eff_data_trig_m - 
		     eff_data_trig_e*eff_data_trig_mlow*dzFilterEff_data -
		     eff_data_trig_m*eff_data_trig_elow*dzFilterEff_data +
		     eff_data_trig_e*eff_data_trig_m*dzFilterEff_data;
      float eff_comb_mc =
		     eff_mc_trig_e +
		     eff_mc_trig_m +
		     eff_emu_mc*dzFilterEff_mc -
		     eff_mc_trig_e*eff_mc_trig_m -
		     eff_mc_trig_e*eff_mc_trig_mlow*dzFilterEff_mc -
		     eff_mc_trig_m*eff_mc_trig_elow*dzFilterEff_mc +
		     eff_mc_trig_e*eff_mc_trig_m*dzFilterEff_mc;

      if (eff_comb_mc<1e-3||eff_comb_data<1e-4)
	trigweight_comb = 0.0;
      else
	trigweight_comb = eff_comb_data/eff_comb_mc;
      
      float idiso_eff = idisoweight_1 * trkeffweight_1 * idisoweight_2 * trkeffweight_2;
      effweightEMu = idiso_eff * trigweightEMu;
      effweightSingle = idiso_eff * trigweightSingle;
      effweightComb = idiso_eff * trigweight_comb;
      if (triggerOption==0)
	effweight = effweightEMu;
      else if (triggerOption==1)
	effweight = effweightSingle;
      else
	effweight = trigweight_comb;
    }

    /*
    if (tree17->iso_1<0.15&&tree17->iso_2<0.2&&(tree17->pt_1>100.0||tree17->pt_2>100.0)) {

      std::cout << "pt1 = " << tree17->pt_1 << "    pt2 = " << tree17->pt_2 << std::endl;

      std::cout << "Effweight (e-mu)      : KIT = " << tree17->effweightEMu << "   IC = " << effweightEMu << std::endl;

      std::cout << "Effweight (singleLep) : KIT = " << tree17->effweightSingle << "   IC = " << effweightSingle << std::endl;

      std::cout << "Effweight (comb)      : KIT = " << tree17->effweight << "   IC = " << effweightComb << std::endl;
      // ******
      // ******
      std::cout << "Trgweight (e-mu)      : KIT = " << tree17->trigweightEMu << "   IC = " << trigweightEMu << std::endl;

      std::cout << "Trgweight (singleLep) : KIT = " << tree17->trigweightSingle << "   IC = " << trigweightSingle << std::endl;

      std::cout << "Trgweight (comb)      : KIT = " << tree17->trigweight << "   IC = " << trigweight_comb << std::endl;

    }
    */
    if (isEmbedded) {
      weight = tree17->mcweight*tree17->embweight*effweight;
    }
    else if (isData) {
      weight = 1;
    }
    else {
      weight = tree17->mcweight*tree17->puweight*effweight;
      if (is2016||is2017)
	weight *= tree17->prefiringweight;
    }
    if (isTTbar) 
      weight *= tree17->topptweight;
    if (isDY)
      weight *= tree17->zptweight;

    if (tree17->os<0.5)
      weight *= tree17->qcdweight;

    weight *= norm;

    TString iSign = "SS";
    if (tree17->os>0.5)
      iSign = "OS";

    TString ijet = "0jet";
    if (tree17->njets==1)
      ijet = "1jet";
    if (tree17->njets>1)
      ijet = "2jet";

    TString ibtag = "nobtag";
    if (tree17->nbtag>=1)
      ibtag = "btag";

    // determination region
    if (tree17->iso_1<0.15&&tree17->iso_2>0.2&&tree17->pzeta>-35.0) {

      th1["deltaR_"+ibtag+"_"+iSign]->Fill(tree17->dr_tt,weight);
      th1["deltaR_"+ijet+"_"+iSign]->Fill(tree17->dr_tt,weight);
      th1["deltaR_inclusive_"+iSign]->Fill(tree17->dr_tt,weight);

      th1["pt1_"+ibtag+"_"+iSign]->Fill(tree17->pt_1,weight);
      th1["pt1_inclusive_"+iSign]->Fill(tree17->pt_1,weight);

      th1["pt2_"+ibtag+"_"+iSign]->Fill(tree17->pt_2,weight);
      th1["pt2_inclusive_"+iSign]->Fill(tree17->pt_2,weight);

      th1["mTtot_"+ibtag+"_"+iSign]->Fill(tree17->mt_tot,weight);
      th1["mTtot_inclusive_"+iSign]->Fill(tree17->mt_tot,weight);

      th2["pt1_pt2_inclusive_"+iSign]->Fill(tree17->pt_1,tree17->pt_2,weight);
      th2["pt1_pt2_"+ibtag+"_"+iSign]->Fill(tree17->pt_1,tree17->pt_2,weight);

    }

    // ttbar region
    if (tree17->iso_1<0.15&&tree17->iso_2>0.2&&tree17->pzeta<=-35.0) {
      th1["mTtot_ttbar_"+iSign]->Fill(tree17->mt_tot,weight);
    }

    //    std::cout << "m_vis = " << m_vis << std::endl;
    //    std::cout << "pzeta = " << pzeta << std::endl;
    //    std::cout << "mTtot = " << mt_tot << std::endl;

    // isolated muon, antiisolated electron
    if (tree17->iso_2<0.2&&tree17->iso_1>0.15&&tree17->pzeta>-35.0) {

      th1["deltaR_inclusive_isoMu_antiisoE_"+iSign]->Fill(tree17->dr_tt,weight);     
      th1["pt1_inclusive_isoMu_antiisoE_"+iSign]->Fill(tree17->pt_1,weight);
      th1["pt2_inclusive_isoMu_antiisoE_"+iSign]->Fill(tree17->pt_2,weight);
      th2["pt1_pt2_inclusive_isoMu_antiisoE_"+iSign]->Fill(tree17->pt_1,tree17->pt_2,weight);

    }

    // antiisolated muon, antiisolated electron
    if (tree17->iso_2>0.2&&tree17->iso_1>0.15&&tree17->pzeta>-35.0) {

      th1["deltaR_inclusive_antiisoMu_antiisoE_"+iSign]->Fill(tree17->dr_tt,weight);
      th1["pt1_inclusive_antiisoMu_antiisoE_"+iSign]->Fill(tree17->pt_1,weight);
      th1["pt2_inclusive_antiisoMu_antiisoE_"+iSign]->Fill(tree17->pt_2,weight);
      th2["pt1_pt2_inclusive_antiisoMu_antiisoE_"+iSign]->Fill(tree17->pt_1,tree17->pt_2,weight);

    }
    // shapes comparison 
    if (tree17->os<0.5) {
      // direct muon iso
      if (tree17->iso_1<0.15&&tree17->iso_2<0.2&&tree17->pzeta>-35.0) {
	th1["mTtot_inclusive_isoMu_SS"]->Fill(tree17->mt_tot,weight);
	th1["mTtot_"+ibtag+"_isoMu_SS"]->Fill(tree17->mt_tot,weight);
      }
      if  (tree17->iso_1<0.15&&tree17->iso_2<0.2&&tree17->pzeta<-35.0) {
	th1["mTtot_ttbar_isoMu_SS"]->Fill(tree17->mt_tot,weight);
      }
      // loose muon iso
      if (tree17->iso_1<0.15&&tree17->iso_2<0.4&&tree17->pzeta>-35.0) {
	th1["mTtot_inclusive_looseisoMu_SS"]->Fill(tree17->mt_tot,weight);
	th1["mTtot_"+ibtag+"_looseisoMu_SS"]->Fill(tree17->mt_tot,weight);
      }
      if  (tree17->iso_1<0.15&&tree17->iso_2<0.4&&tree17->pzeta<-35.0) {
	th1["mTtot_ttbar_looseisoMu_SS"]->Fill(tree17->mt_tot,weight);
      }      
    }

    if (tree17->iso_1<0.15&&tree17->iso_2<0.2&&tree17->pzeta>-35.0) { 
      th1["mvis_sel_signal_"+iSign]->Fill(tree17->m_vis,weight);
      th1["pt1_sel_signal_"+iSign]->Fill(tree17->m_vis,weight);
      th1["pt2_sel_signal_"+iSign]->Fill(tree17->m_vis,weight);
      th1["mTtot_sel_signal_"+iSign]->Fill(tree17->m_vis,weight);
      if (tree17->os>0.5) count += 1;
    }

    if (tree17->iso_1<0.15&&tree17->iso_2<0.2&&tree17->pzeta<=-35.0) { 
      th1["mvis_sel_ttbar_"+iSign]->Fill(tree17->m_vis,weight);
      th1["pt1_sel_ttbar_"+iSign]->Fill(tree17->m_vis,weight);
      th1["pt2_sel_ttbar_"+iSign]->Fill(tree17->m_vis,weight);
      th1["mTtot_sel_ttbar_"+iSign]->Fill(tree17->m_vis,weight);
    }

    if (tree17->iso_1<0.15&&tree17->iso_2<0.2) {
      th1["pzeta_sel_"+iSign]->Fill(tree17->pzeta,weight);
    }

  }
  //  std::cout << "   " << sampleName << "  :  " << count << std::endl;
  delete tree17;

} 

void ProcessSample(SampleAttributes sample,
		   TFile * output) {

  std::map<TString,TFile*> sampleMap = sample.fileMap;  
  std::map<TString,float> normMap = sample.normMap;
  TString name = sample.name;

  std::cout << "Processing sample : " << sample.name << std::endl;

  std::map<TString,TH1D*> th1;
  std::map<TString,TH2D*> th2;
  
  CreateMaps(th1,th2,name);

  for (auto mapIter : sampleMap) {
    TString sampleName = mapIter.first; 
    TFile * file = mapIter.second;
    TTree * tree = (TTree*)file->Get(BaseTreeName);
    float norm = normMap[sampleName];
    int selectDiTau = sample.selectDiTau;
    std::cout << "       " << sampleName << "      events in tree : " << tree->GetEntries() << std::endl;
    //    std::cout << "selectDiTau = " << selectDiTau << std::endl;
    //    std::cout << "norm = " << norm << std::endl;
    RunOnTree(tree,name,sampleName,th1,th2,selectDiTau,norm);
  }

  output->cd("");
  for (auto th1_iter : th1) {
    TString histname = th1_iter.first;
    TH1D * hist = th1_iter.second;
    hist->Write(name+"_"+histname);
  }
  for (auto th2_iter : th2) {
    TString histname = th2_iter.first;
    TH2D * hist = th2_iter.second;
    hist->Write(name+"_"+histname);
  }

}


int main(int argc, char * argv[]) {
  
  // argument - config file
  if (argc!=5) {
    std::cout << "usage of the scripts : QCD_emu [era=2016,2017,2018] [trigger Option = SingleLep, EMu, Comb] [Sample = Data, Embedded, MC, TT, All] [SF = SFnew, SFold]" << std::endl;
    exit(-1);
  }

  // **** input parameters
  TString era(argv[1]);
  TString Trigger(argv[2]);
  TString Sample(argv[3]);
  TString SF(argv[4]);

  if (era!="2016"&&era!="2017"&&era!="2018") {
    std::cout << "Unknown era : " << era << std::endl;
    exit(-1);
  }
  if (Trigger!="SingleLep"&&Trigger!="EMu"&&Trigger!="Comb") {
    std::cout << "Unknown trigger option : " << Trigger << std::endl;
    exit(-1);
  }
  if (Sample!="Data"&&Sample!="Embedded"&&Sample!="MC"&&Sample!="All"&&Sample!="TT") {
    std::cout << "Unknown sample option : " << Sample << std::endl;
  }

  if (SF=="SFold") 
    applySF = false;
  if (SF=="SFnew")
    applySF = true;

  if (era=="2018") {
    lumi = 59740;
    pt_SingleE = 33.0;
    pt_SingleMu = 25.0;
  }
  if (era=="2017") {
    lumi = 41900;
    pt_SingleE = 28.0;
    pt_SingleMu = 25.0;
  }
  if (era=="2016") {
    lumi = 35866;
    pt_SingleE = 26.0;
    pt_SingleMu = 23.0;
  }

  if (Trigger=="SingleLep")
    triggerOption = 1;
  else if (Trigger=="EMu") 
    triggerOption = 0;
  else if (Trigger=="Comb")
    triggerOption = 2;

  is2017 = false;
  is2018 = false;
  is2016 = false;
  dzFilterEff_mc = 0.98;
  dzFilterEff_data = 1.0;
  if (era=="2017") {
    is2017 = true;
    dzFilterEff_mc = 0.95;
    dzFilterEff_data = 0.95;
  }
  if (era=="2018") {
    dzFilterEff_mc = 0.95;
    dzFilterEff_data = 0.95;
  }

  
  string cmsswBase = (getenv("CMSSW_BASE"));

  TString workspace_filename = TString(cmsswBase) + "/src/DesyTauAnalyses/NTupleMaker/data/CorrectionWS_IC/htt_scalefactors_legacy_"+era+".root";

  TFile * f_workspace = new TFile(workspace_filename);
  w = (RooWorkspace*)f_workspace->Get("w");

  input_dir = "/nfs/dust/cms/user/rasp/grid-jobs/emu_MSSM_Dec20/"+era;
  output_dir = "./";

  std::cout << "Base directory : " << input_dir << std::endl;

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
				      
  if (triggerOption==1) {

    for (unsigned int i=0; i<SingleElectron.size(); ++i)
      DataSamples.push_back(SingleElectron.at(i));

    for (unsigned int i=0; i<SingleMuon.size(); ++i)
      DataSamples.push_back(SingleMuon.at(i));

  }
  if (triggerOption==0) {
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

  SampleAttributes DataAttr = SetSampleAttributes("Data","data",DataSamples,0);
  SampleAttributes EmbedAttr = SetSampleAttributes("Embedded","EmbedZTT",EmbedSamples,0);
  SampleAttributes ZllAttr = SetSampleAttributes("DYJetsLL","ZLL",DYSamples,-1);
  SampleAttributes WJetsAttr = SetSampleAttributes("WJets","W",WJetsSamples,-1);
  SampleAttributes VVAttr = SetSampleAttributes("VV","VV",EWKSamples,-1);
  SampleAttributes TTAttr = SetSampleAttributes("TT","TT",TTSamples,-1);

  TString OutputFileName = Sample + "_" + Trigger + "_" + SF + "_" + era + ".root";
  TFile * outputFile = new TFile(OutputFileName,"recreate");
  outputFile->cd("");
  
  if (Sample=="Data")
    ProcessSample(DataAttr,outputFile);
  if (Sample=="Embedded")
    ProcessSample(EmbedAttr,outputFile);
  if (Sample=="MC") {
    ProcessSample(ZllAttr,outputFile);
    ProcessSample(WJetsAttr,outputFile);
    ProcessSample(VVAttr,outputFile);
  }
  if (Sample=="TT") 
    ProcessSample(TTAttr,outputFile);
  if (Sample=="All") {
    ProcessSample(DataAttr,outputFile);
    ProcessSample(EmbedAttr,outputFile);
    ProcessSample(ZllAttr,outputFile);
    ProcessSample(WJetsAttr,outputFile);
    ProcessSample(VVAttr,outputFile);
    ProcessSample(TTAttr,outputFile);
  }


  outputFile->Close();
  delete outputFile;

}
