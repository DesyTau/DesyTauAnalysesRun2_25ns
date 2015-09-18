#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TEntryList.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"

//#include "HiggsAnalysis/CombinedLimit/interface/TH1Keys.h"
using namespace std;

#define MSSM             false
#define VERBOSE          false
#define DEBUG            true
//#define LOOP             true
#define scaleByBinWidth  false
#define DOSPLIT          true
#define useZDataMC       false

typedef map<TString, TChain* >  mapchain;

///////////////////////////////////////////////////////////////////////////////////////////////
TH1F* blindHistogram(TH1F* h, float xmin, float xmax, TString name) {

  TH1F* hOut = (TH1F*)h->Clone(name);
  int nBins  = hOut->GetNbinsX();
  float x    = 0;

  for(int i=1; i<nBins+1 ; i++) {
    x = hOut->GetBinCenter(i);
    if( x>xmin && x<xmax ) hOut->SetBinContent(i,-10);
  }

  return hOut;

}
///////////////////////////////////////////////////////////////////////////////////////////////

void makeHistoFromDensity(TH1* hDensity, TH1* hHistogram){

  if(hDensity->GetNbinsX() != hHistogram->GetNbinsX()){
    cout << "makeHistoFromDensity: different binning" << endl;
    return;
  }

  for(int k = 1 ; k <= hDensity->GetNbinsX(); k++){
    float bink   = hDensity->GetBinContent(k);
    float widthk = hHistogram->GetBinWidth(k);
    hDensity->SetBinContent(k, bink*widthk );
  }
  hDensity->Scale(hHistogram->Integral()/hDensity->Integral());
}

///////////////////////////////////////////////////////////////////////////////////////////////

void chooseSelection(TString variable_,
		     TString version_,
		     string selection_, 
		     TCut& tiso,
		     TCut& tdecaymode,
		     TCut& ltiso,
		     TCut& atiso,
		     TCut& antimu,
		     TCut& antiele,
		     TCut& selection,
		     TCut& tpt
		     )
{
  cout<<"VERSION in chooseSelection :"<< version_<<endl;

  if(version_.Contains("SVfitMassCut")) selection = selection && "diTauNSVfitMass>50.";

  else if(version_.Contains("taupt45")){
    if(version_.Contains("taupt4560"))
      tpt="ptL1>45 && ptL1<60 && ptL2>45 && ptL2<60";
    else
      tpt="ptL2>45";
  }
  else                            
    tpt="ptL1>45 && ptL2>45";      
  
  if(version_.Contains("2bTagged")) 
    selection =" nJets20BTagged>1";

  // Anti-Mu discriminator //
  if(version_.Contains("AntiMu3Loose"))    antimu = "tightestAntiMu3WPL1>0 && tightestAntiMu3WPL2>0";
  else if(version_.Contains("AntiMu3Tight"))    antimu = "tightestAntiMu3WPL1>1 && tightestAntiMu3WPL2>1";
  else if(version_.Contains("AntiMuMVALoose"))  antimu = "tightestAntiMuMVAWPL1>0 && tightestAntiMuMVAWPL2>0";
  else if(version_.Contains("AntiMuMVAMedium")) antimu = "tightestAntiMuMVAWPL1>1 && tightestAntiMuMVAWPL2>1";
  else if(version_.Contains("AntiMuMVATight"))  antimu = "tightestAntiMuMVAWPL1>2 && tightestAntiMuMVAWPL2>2";

  // Anti-Ele discriminator //
  if(version_.Contains("AntiEleLoose"))        antiele = "tightestAntiECutWPL1 > 0 && tightestAntiECutWPL2 > 0";
  else if(version_.Contains("AntiEleMedium"))  antiele = "tightestAntiECutWPL1 > 1 && tightestAntiECutWPL2 > 1";
  else if(version_.Contains("AntiEleTight"))   antiele = "tightestAntiECutWPL1 > 2 && tightestAntiECutWPL2 > 2";
  else if(version_.Contains("AntiEle5VLoose")) antiele = "tightestAntiEMVA5WPL1 > 0 && tightestAntiEMVA5WPL2 > 0";
  else if(version_.Contains("AntiEle5Loose"))  antiele = "tightestAntiEMVA5WPL1 > 1 && tightestAntiEMVA5WPL2 > 1";
  else if(version_.Contains("AntiEle5Medium")) antiele = "tightestAntiEMVA5WPL1 > 2 && tightestAntiEMVA5WPL2 > 2";
  else if(version_.Contains("AntiEle5Tight"))  antiele = "tightestAntiEMVA5WPL1 > 3 && tightestAntiEMVA5WPL2 > 3";
  else if(version_.Contains("AntiEle5VTight")) antiele = "tightestAntiEMVA5WPL1 > 4 && tightestAntiEMVA5WPL2 > 4";

  //TauID
  if(version_.Contains("TauOldDM")) tdecaymode = "decayModeFindingOldDM1>0.5 && decayModeFindingOldDM2>0.5";
  else if(version_.Contains("TauNewDM")) tdecaymode = "decayModeFindingNewDM1>0.5 && decayModeFindingNewDM2>0.5";

  // TauIso
  // TauIso DB3Hits cut-based //
  if(version_.Contains("HPSDB3H")) {
    tiso   = "hpsDB3HL1<1.0 && hpsDB3HL2<1.0" ;
    ltiso  = "hpsDB3HL1<10.0 && hpsDB3HL2<10.0" ;
    atiso  = "hpsDB3HL1<10.0 && hpsDB3HL2<10.0 && hpsDB3HL1>1.0 && hpsDB3HL2>1.0"; 
  }
  else if(version_.Contains("HPSMVA3oldDMwLTLoose")) {
    tiso   = "tightestHPSMVA3oldDMwLTWPL1>0 && tightestHPSMVA3oldDMwLTWPL2>0" ;//Tight 3
    ltiso   = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1" ;//Loose 1
    atiso  = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1 && tightestHPSMVA3oldDMwLTWPL1<=0 && tightestHPSMVA3oldDMwLTWPL2<=0";
  }
  else if(version_.Contains("HPSMVA3oldDMwLTMedium")) {
    tiso   = "tightestHPSMVA3oldDMwLTWPL1>1 && tightestHPSMVA3oldDMwLTWPL2>1" ;//Tight 3
    ltiso   = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1" ;//Loose 1
    atiso  = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1 && tightestHPSMVA3oldDMwLTWPL1<=1 && tightestHPSMVA3oldDMwLTWPL2<=1";
  }
  else if(version_.Contains("HPSMVA3oldDMwLTTight")) {
    tiso   = "tightestHPSMVA3oldDMwLTWPL1>2 && tightestHPSMVA3oldDMwLTWPL2>2" ;//Tight 3
    ltiso   = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1" ;//Loose 1
    atiso  = "tightestHPSMVA3oldDMwLTWPL1>-1 && tightestHPSMVA3oldDMwLTWPL2>-1 && tightestHPSMVA3oldDMwLTWPL1<=2 && tightestHPSMVA3oldDMwLTWPL2<=2";
  }
  else if(variable_.Contains("hpsMVA3oldDMwLT")) {
    tiso   = "etaL1<999 && etaL2<999" ;//Tight 3
    ltiso   = "etaL1<999 && etaL2<999" ;//Loose 1
    atiso  = "etaL1<999 && etaL2<999" ;//Loose 1
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////

void drawHistogram(TCut sbinCat,
		   TString type,
		   TString version_,
		   TString analysis_,
		   TChain* tree, 
		   TString variable, 
		   float& normalization, 
		   float& normalizationError, 
		   float scaleFactor,
		   TH1F* h,
		   TCut cut,
		   int verbose = 0){

  if(DEBUG) 
    cout << "Start drawHistogram : " << type << " " << version_ << endl;

  // Start processing
  if(tree!=0 && h!=0){

    // HLT matching //
    //TCut hltMatch("run>0");
    // GenMass cut //
    TCut genMass("run>0");

    // Reweighting
    TCut weight         = "run>0";
    //TCut sampleWeight   = "run>0";
    //TCut weightDY = "run>0";
    //TCut weightW  = "run>0";
    //TCut weightEmb= "run>0";

    if(type.Contains("MC")) {

      weight = "(weight*mcweight)";

      /* //to be used when weight is available
	 if(!type.Contains("SUSY"))
	 {
	 if(type.Contains("GGFHUp")) 
	 weight *= "HqTWeightUp"; 
	 else if(type.Contains("GGFHDown"))  
	 weight *= "HqTWeightDown";
	 else
	 weight *= "HqTWeight";
	 }
	 else
	 {
	 if(type.Contains("GGHUp")) 
	 weight *= "mssmHiggsPtReweightGluGlu_mhmax_tanBetaUp"; 
	 else if(type.Contains("GGHDown"))  
	 weight *= "mssmHiggsPtReweightGluGlu_mhmax_tanBetaDown";
	 else
	 weight *= "mssmHiggsPtReweightGluGlu_mhmax";
	 }
      */	 
    }
    /* //To be used when Embedded sample is available
      else if(type.Contains("Embed")) {
      genMass     = "genDiTauMass>50 && HLTxMu17Mu8>0.5"; // HLTxMu17Mu8
      weightEmb   = "embeddingWeight";
      
      weight = "(HLTTau*HLTMu*SFTau*SFMuID*ZmmWeight*weightDecayMode)";
      }
      
      if(version_.Contains("DecayModeNoCorr"))
      {
      weight = "(HLTTau*HLTMu*SFTau*SFMuID*ZmmWeight)";
      }
      }

    
    if(type.Contains("TTJetsUp")) 
      weight *= "topPtWeightUp"; 
    else if(type.Contains("TTJetsDown")) 
      weight *= "topPtWeightDown"; 
    else if(type.Contains("TTJets")) 
      weight *= "topPtWeightNom"; 
    */
    //     if(type.Contains("SUSY")) cout<<"weight : "<<weight<<endl;
    cout<<"weight : "<<weight<<endl;
    
    TCut pairIndex="pairIndex<1";

    tree->Draw(variable+">>"+TString(h->GetName()),cut*weight*sbinCat*genMass*pairIndex);

    // Scale the histogram, compute norm and err
    h->Scale(scaleFactor);
    normalization      = h->Integral();
    normalizationError = TMath::Sqrt(h->GetEntries()) * (normalization/h->GetEntries());

    //cout << h->GetEntries() << " entries => integral=" << normalization << endl;
    if(DEBUG) cout << h->GetEntries() << " entries => integral=" << normalization << endl;
    if(verbose==0) h->Reset();

  }
  else{
    cout << "ERROR : drawHistogram => null pointers to tree or histogram. Exit." << endl;
    return;
  }
}

TArrayF createBins(int nBins_ = 80 ,
		   float xMin_ = 0.,
		   float xMax_ = 400.,
		   int& nBins = *(new int()),
		   string selection_   = "inclusive",
		   TString variable_   = "diTauVisMass",
		   TString location    = "/nfs/dust/cms/user/anayak/CMS/OnSLC6/CMSSW_746p6_htt/src/DesyTauAnalyses/NTupleMaker/bin/binning/"
		   ){

  // input txt file with bins
  ifstream is;

  TArrayF dummy(2);
  dummy[0] = xMin_; dummy[1] = xMax_;
  
  char* c = new char[10];
  string filename = Form(location+"/bins_diTau_%s_%s.txt",variable_.Data(), selection_.c_str());
  if(MSSM)
    filename = Form(location+"/bins_diTau_%s_%s_mssm.txt",variable_.Data(), selection_.c_str());
  
  is.open(filename);
  if(nBins_<0 &&  !is.good()){
    cout << "Bins file not found" << endl;
    return dummy;
  }

  int nBinsFromFile = 0;
  while (is.good())     
    {
      is.getline(c,999,',');     
      if (is.good()){
	nBinsFromFile++;
	//cout << c << endl;
      }
    }

  // choose the number of bins
  nBins =  nBins_>0 ? nBins_ : nBinsFromFile-1 ;
  TArrayF bins(nBins+1);
  cout << "Making histograms with " << nBins << " bins:" << endl;

  is.close();
  //is.open(Form(location+"/bins_tauTau_%s_%s.txt",variable_.Data(), selection_.c_str())); 
  is.open(filename);
  
  nBinsFromFile = 0;

  if(nBins_>0){
    for( ; nBinsFromFile <= nBins ; nBinsFromFile++){
      bins[nBinsFromFile] =  xMin_ + nBinsFromFile*(xMax_-xMin_)/nBins_;
    }
  }
  else{
    while (is.good())  
      {
	is.getline(c,999,',');     
	if (is.good() && nBinsFromFile<=nBins) {
	  bins[nBinsFromFile] = atof(c);
	  cout << bins[nBinsFromFile] << ", " ;
	}
	nBinsFromFile++;
      }
    cout << endl;
  }

  return bins;

}

void evaluateQCD(mapchain mapAllTrees, TString version_, TString analysis_,
		 TH1F* qcdHisto, TH1F* ssHisto, string sign,
		 float& SSQCDinSignalRegionDATAIncl_, float& scaleFactorTTSSIncl,
		 TH1F* hExtrap, TString variable = "",
		 float scaleFactor=0., float TTxsectionRatio=0., float lumiCorrFactor = 0.,
		 float ExtrapolationFactorZDataMC = 0.,
		 float JtoTauCorrectionFactor=0., float LtoTauCorrectionFactor=0.,
		 float OStoSSRatioQCD = 0.,
		 TCut sbin = "", TCut sbinCat="", 
		 TCut sbinSS="", TCut sbinSSCat="",
		 TCut sbinAiso = "",
		 bool subtractTT=true, bool subtractVV=true){

  

  //Data SS in signal region
  float Error = 0.;
  float SSQCDinSignalRegionDATAIncl = 0.;
  drawHistogram(sbinSSCat,"Data", version_,analysis_, mapAllTrees["Data"], variable,              SSQCDinSignalRegionDATAIncl,        Error, 1.0,         hExtrap, sbinSS, 1);
  cout<<"********* HERE IS STANDARD MEASUREMENT OF QCD IN SS *********"<<endl;
  cout<<"      sbinSSCat = "<<sbinSSCat<<endl;
  cout<<"      sbin = "<<sbin<<endl;
  cout<<"      sbinSS = "<<sbinSS<<endl;
  cout<<"********* END *********"<<endl;
  if(ssHisto !=0) ssHisto->Add( hExtrap,  1.0);

  //W+jets SS
  float SSWJetsinSidebandRegionMCIncl    = 0.;
  drawHistogram(sbinSSCat, "MC",version_,analysis_, mapAllTrees["WJets"],     variable, SSWJetsinSidebandRegionMCIncl,      Error, scaleFactor, hExtrap, sbinSS,1);
  if(ssHisto !=0) ssHisto->Add( hExtrap,  1.0);

  //TTbar SS MC
  float SSTTbarinSidebandRegionMCIncl    = 0.;
  if(subtractTT) {
    drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["TTbar"],     variable, SSTTbarinSidebandRegionMCIncl,      Error, scaleFactor*TTxsectionRatio*scaleFactorTTSSIncl,       hExtrap, sbinSS,1);
    if(ssHisto !=0) ssHisto->Add( hExtrap, -1.0);
  }

  //Diboson SS MC
  float SSOthersinSidebandRegionMCIncl    = 0.;
  if(subtractVV) {
    drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["Others"],     variable, SSOthersinSidebandRegionMCIncl,     Error, scaleFactor,       hExtrap, sbinSS,1);
    if(ssHisto !=0) ssHisto->Add( hExtrap, -1.0);
  }

  //DY SS MC L To Tau
  float SSDYLtoTauinSidebandRegionMCIncl = 0.;
  drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["DYLtoTau"], variable, SSDYLtoTauinSidebandRegionMCIncl,  Error, lumiCorrFactor*scaleFactor*LtoTauCorrectionFactor,    hExtrap, sbinSS,1);
  if(ssHisto !=0) ssHisto->Add( hExtrap, -1.0);

  //DY SS MC Tau Tau
  float SSDYtoTauinSidebandRegionMCIncl = 0.;
  if(useZDataMC)
    drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["DYToTauTau"],  variable, SSDYtoTauinSidebandRegionMCIncl,    Error, lumiCorrFactor*scaleFactor*ExtrapolationFactorZDataMC, hExtrap, sbinSS,1);
  else
    drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["DYToTauTau"],  variable, SSDYtoTauinSidebandRegionMCIncl,    Error, lumiCorrFactor*scaleFactor, hExtrap, sbinSS,1);
  if(ssHisto !=0) ssHisto->Add( hExtrap, -1.0);

  //DY SS MC Jet To Tau
  float SSDYJtoTauinSidebandRegionMCIncl = 0.;
  drawHistogram(sbinSSCat,"MC", version_,analysis_, mapAllTrees["DYJtoTau"],  variable, SSDYJtoTauinSidebandRegionMCIncl,   Error, lumiCorrFactor*scaleFactor*JtoTauCorrectionFactor,     hExtrap, sbinSS,1);
  if(ssHisto !=0) ssHisto->Add( hExtrap, -1.0);

  cout << "Selected events in inclusive " << sign << " data " << SSQCDinSignalRegionDATAIncl << endl;

  SSQCDinSignalRegionDATAIncl  -= SSWJetsinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl  -= SSTTbarinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl  -= SSOthersinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl  -= SSDYLtoTauinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl  -= SSDYJtoTauinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl  -= SSDYtoTauinSidebandRegionMCIncl;
  SSQCDinSignalRegionDATAIncl *= OStoSSRatioQCD;
  //if(qcdHisto!=0) qcdHisto->Scale(OStoSSRatioQCD);

  //All non-QCD contributions subtracted from data, taken in the SR (OS, iso) --> dump
  cout << "- expected from WJets          " << SSWJetsinSidebandRegionMCIncl << endl;
  cout << "- expected from TTbar          " << SSTTbarinSidebandRegionMCIncl << endl;
  cout << "- expected from Others         " << SSOthersinSidebandRegionMCIncl << endl;
  cout << "- expected from DY->tautau     " << SSDYtoTauinSidebandRegionMCIncl << endl;
  cout << "- expected from DY->ll, l->tau " << SSDYLtoTauinSidebandRegionMCIncl << endl;
  cout << "- expected from DY->ll, j->tau " << SSDYJtoTauinSidebandRegionMCIncl  << endl;
  cout << "QCD in inclusive SS region is estimated to be " << SSQCDinSignalRegionDATAIncl/OStoSSRatioQCD  << "*" << OStoSSRatioQCD
       << " = " <<  SSQCDinSignalRegionDATAIncl << endl;
  SSQCDinSignalRegionDATAIncl_ = SSQCDinSignalRegionDATAIncl;

  //Data OS anti-iso
  float OSQCDinSignalRegionDATAInclaIso = 0.;
  drawHistogram(sbinCat,"Data", version_,analysis_, mapAllTrees["Data"], variable, OSQCDinSignalRegionDATAInclaIso,    Error, 1.0, hExtrap, sbinAiso);
  if(qcdHisto!=0) qcdHisto->Add(hExtrap,  1.0);

  //W+jets OS anti-iso
  float OSWJetsinSignalRegionMCInclaIso  = 0.;
  drawHistogram(sbinCat, "MC",version_,analysis_, mapAllTrees["WJets"],     variable, OSWJetsinSignalRegionMCInclaIso,      Error, 1.0, hExtrap, sbinAiso,1);
  if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);

  //TTbar OS MC anti-iso
  float OSTTbarinSignalRegionMCInclaIso  = 0.;
  if(subtractTT) {
    drawHistogram(sbinCat,"MCTTJets", version_,analysis_, mapAllTrees["TTbar"],     variable, OSTTbarinSignalRegionMCInclaIso,    Error, scaleFactor*TTxsectionRatio,       hExtrap, sbinAiso,1);
    if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);
  }

  //Diboson OS MC anti-iso
  float OSOthersinSignalRegionMCInclaIso    = 0.;
  if(subtractVV) {
    drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["Others"],     variable, OSOthersinSignalRegionMCInclaIso,     Error, scaleFactor,       hExtrap, sbinAiso,1);
    if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);
  }

  //DY OS MC Mu To Tau
  float OSDYLtoTauinSignalRegionMCInclaIso = 0.;
  drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["DYLtoTau"], variable, OSDYLtoTauinSignalRegionMCInclaIso,  Error, lumiCorrFactor*scaleFactor*LtoTauCorrectionFactor,    hExtrap, sbinAiso,1);
  if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);

  //DY OS MC Tau Tau
  float OSDYtoTauinSignalRegionMCInclaIso = 0.;
  if(useZDataMC)
    drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["DYToTauTau"],  variable, OSDYtoTauinSignalRegionMCInclaIso,    Error, lumiCorrFactor*scaleFactor*ExtrapolationFactorZDataMC, hExtrap, sbinAiso,1);
  else
    drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["DYToTauTau"],  variable, OSDYtoTauinSignalRegionMCInclaIso,    Error, lumiCorrFactor*scaleFactor, hExtrap, sbinAiso,1);
  if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);

  //DY OS MC Jet To Tau
  float OSDYJtoTauinSidebandRegionMCInclaIso = 0.;
  drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["DYJtoTau"],  variable, OSDYJtoTauinSidebandRegionMCInclaIso,   Error, lumiCorrFactor*scaleFactor*JtoTauCorrectionFactor,     hExtrap, sbinAiso,1);
  if(qcdHisto!=0) qcdHisto->Add(hExtrap, -1.0);

  if(qcdHisto!=0)qcdHisto->Scale(SSQCDinSignalRegionDATAIncl/qcdHisto->Integral());

  //Do separately for VBF
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////7

void plotTauTau( Int_t mH_           = 120,
		Int_t useEmbedding_ = 0,
		string selection_   = "inclusive",
		TString analysis_   = "",		  
		TString variable_   = "diTauVisMass",
		TString XTitle_     = "Visible mass",
		TString Unities_    = "GeV",
		TString outputDir   = "test",
		Int_t nBins_ = 80, Float_t xMin_=0, Float_t xMax_=400,
		Float_t magnifySgn_ = 1.0,
		Int_t logy_         = 0,
		Float_t maxY_       = 1.2,
		TString version_    = "v1",
		TString location    = ""
		) 
{   

  location = gSystem->pwd();
  location += "/results/";

  cout << endl;
  cout << "@@@@@@@@@@@@@@@@@@ Category  = " << selection_     <<  endl;
  cout << "@@@@@@@@@@@@@@@@@@ Variable  = " << string(variable_.Data()) <<  endl;
  cout << endl;

  ostringstream ossiTmH("");
  ossiTmH << mH_ ;
  TString TmH_ = ossiTmH.str();

  //const int nProd=3;
  //const int nMasses=0;
  //TString nameProd[nProd]={"GGFH","VBFH","VH"};
  ////int hMasses[nMasses]={90,95,100,105,110,115,120,125,130,135,140,145,150,155,160};
  //int hMasses[nMasses]={};
  //TString nameMasses[nMasses];
  //
  //if(DEBUG) cout << "build masses string" << endl;
  //for(int iM=0 ; iM<nMasses ; iM++) nameMasses[iM]=TString(Form("%d",hMasses[iM]));

  //const int nProdWW=2;
  //const int nMassesWW=11;
  //TString nameProdWW[nProdWW]={"GGFHWW","VBFHWW"};
  //int hMassesWW[nMassesWW]={110,115,120,125,130,135,140,145,150,155,160};
  //TString nameMassesWW[nMassesWW];
  //
  //if(DEBUG) cout << "build masses string" << endl;
  //for(int iM=0 ; iM<nMassesWW ; iM++) nameMassesWW[iM]=TString(Form("%d",hMassesWW[iM]));

  const int nProdS=2;
  const int nMassesS=1;
  TString nameProdS[nProdS]={"GGH","BBH"};
  //int hMassesS[nMassesS]={80,90,100,110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000};
  int hMassesS[nMassesS]={160};
  TString nameMassesS[nMassesS];
  if(DEBUG) cout << "build masses string" << endl;
  for(int iM=0 ; iM<nMassesS ; iM++) nameMassesS[iM]=TString(Form("%d",hMassesS[iM]));

  if(DEBUG) cout << "prepare yields file" << endl;
  ofstream out(Form(location+"/%s/yields/yieldsTauTau_mH%d_%s_%s_%s.txt",
		    outputDir.Data(),mH_,selection_.c_str(), analysis_.Data(),variable_.Data() ),ios_base::out); 
  out.precision(5);
  int nBins = nBins_;

  if(DEBUG) cout << "createBins" << endl;
  TArrayF bins = createBins(nBins_, xMin_, xMax_, nBins, selection_, variable_);
  cout<<"Bins : "<<endl;
  for(int i=0 ; i<bins.GetSize() ; i++)cout<<"bin "<<i<<"   "<<bins[i]<<endl;
  
  // LUMINOSITY //
  float Lumi;
  //
  Lumi = 40.24;

  /////////////////

  float lumiCorrFactor                     = 1 ;    
  float TTxsectionRatio                    = 1.0; 
  float OStoSSRatioQCD                     = 1.0; //1.06 at 8 TeV
  float SSIsoToSSAIsoRatioQCD              = 1.0;
  float LtoTauCorrectionFactor            = 1.0;
  float JtoTauCorrectionFactor             = 1.0;
  float ExtrapolationFactorZ               = 1.0;
  float ErrorExtrapolationFactorZ          = 1.0;
  float ExtrapolationFactorZDataMC         = 1.0;
  float ExtrapolationFactorSidebandZDataMC = 1.0;
  float ExtrapolationFactorZFromSideband   = 1.0;
  float scaleFactorTTOS                    = 1.0;
  float scaleFactorTTSS                    = 1.0;
  float scaleFactorTTSSIncl                = 1.0;

  cout << endl;
  cout << "Input: " << endl;
  cout << " > Lumi           = " << Lumi/1000. << " fb-1" << endl;
  cout << " > DY xsection SF = " << lumiCorrFactor << endl;
  cout << " > TTbar SF       = " << TTxsectionRatio << endl;
  cout << " > QCD OS/SS SF   = " << OStoSSRatioQCD << endl;
  cout << " > J->tau SF      = " << JtoTauCorrectionFactor << endl;
  cout << " > l->tau SF     = " << LtoTauCorrectionFactor << endl;
  cout << endl;

  /////////////////  change SVfit mass here ///////////////////

  string variableStr = "";
  TString variable(variableStr.c_str());
  variable = variable_;

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);
  //   if(variable=="dxyErrTau")c1->SetLogx(1);

  TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);
 
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(logy_);
  //   if(variable=="dxyErrTau")pad1->SetLogx(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TLegend* leg = new TLegend(0.53,0.48,0.85,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=13 TeV}{%.2f fb^{-1} #tau_{had}#tau_{had}}", Lumi/1000. ));

  THStack* aStack = new THStack("aStack","");

  TH1F* hSiml     = new TH1F( "hSiml"   ,"all"               , nBins , bins.GetArray());
  TH1F* hSgn      = new TH1F( "hSgn "   ,"vbf+ggf"           , nBins , bins.GetArray());         hSgn->SetFillColor(0); hSgn->SetLineColor(kBlue);hSgn->SetLineWidth(2);hSgn->SetLineStyle(kDashed);
  TH1F* hSgn1     = new TH1F( "hSgn1"   ,"vbf"               , nBins , bins.GetArray());         hSgn1->SetLineWidth(2);
  TH1F* hSgn2     = new TH1F( "hSgn2"   ,"ggf"               , nBins , bins.GetArray());         hSgn2->SetLineWidth(2);
  TH1F* hSgn3     = new TH1F( "hSgn3"   ,"vh"                , nBins , bins.GetArray());         hSgn3->SetLineWidth(2);
  TH1F* hData     = new TH1F( "hData"   ,"        "          , nBins , bins.GetArray());         hData->SetMarkerStyle(20);hData->SetMarkerSize(1.2);hData->SetMarkerColor(kBlack);hData->SetLineColor(kBlack);hData->SetXTitle(XTitle_+" ("+Unities_+")");hData->SetYTitle(Form(" Events/(%.1f %s)", hData->GetBinWidth(1), Unities_.Data() ) );hData->SetTitleSize(0.04,"X");hData->SetTitleSize(0.05,"Y");hData->SetTitleOffset(0.95,"Y");
  TH1F* hZttEmb   = new TH1F( "hZttEmb","Embedded"          , nBins , bins.GetArray());         hZttEmb->SetFillColor(kOrange-4);
  TH1F* hW        = new TH1F( "hW"      ,"W+jets"            , nBins , bins.GetArray());         hW->SetFillColor(kRed+2);
  TH1F* hEWK      = new TH1F( "hEWK"    ,"EWK"               , nBins , bins.GetArray());         hEWK->SetFillColor(kRed+2);
  TH1F* hZtt      = new TH1F( "hZtt"    ,"Ztautau"           , nBins , bins.GetArray());         hZtt->SetFillColor(kOrange-4);
  TH1F* hZll      = new TH1F( "hZll"    ,"Zll"               , nBins , bins.GetArray());         hZll->SetFillColor(kBlue-2);
  TH1F* hZj       = new TH1F( "hZj"    ,"Z+jets, jet to tau", nBins , bins.GetArray());         hZj->SetFillColor(kBlue-2);
  TH1F* hZfakes   = new TH1F( "hZfakes" ,"Z+jets, jet to tau", nBins , bins.GetArray());         hZfakes->SetFillColor(kBlue-2);
  TH1F* hTTb      = new TH1F( "hTTb"    ,"ttbar"             , nBins , bins.GetArray());         hTTb->SetFillColor(kBlue-8); 
  TH1F* hTTbUp    = new TH1F( "hTTbUp"    ,"ttbarUp"         , nBins , bins.GetArray());         hTTbUp->SetFillColor(kBlue-8); 
  TH1F* hTTbDown  = new TH1F( "hTTbDown"    ,"ttbarDown"     , nBins , bins.GetArray());         hTTbDown->SetFillColor(kBlue-8); 
  TH1F* hVV       = new TH1F( "hVV"     ,"Diboson"           , nBins , bins.GetArray());         hVV->SetFillColor(kRed+2);

  TH1F* hSS       = new TH1F( "hSS"     ,"same-sign"         , nBins , bins.GetArray());         hSS->SetFillColor(kMagenta-10);
  TH1F* hQCD      = new TH1F( "hQCD"    ,"QCD"      , nBins , bins.GetArray());                  hQCD->SetFillColor(kMagenta-10);
  TH1F* hDataAntiTauIsoQCD = new TH1F( "hDataAntiTauIsoQCD", "AntiTauIsoQCD", nBins , bins.GetArray());  hDataAntiTauIsoQCD->SetFillColor(kMagenta-10);
  TH1F* hDataAntiTauIso = new TH1F( "hDataAntiTauIso", "AntiTauIso", nBins , bins.GetArray());  hDataAntiTauIso->SetFillColor(kMagenta-10);

  /* // not needed for now
  //W Up/Down
  TH1F* hW_TauFakeUp        = new TH1F( "hW_TauFakeUp"      ,"W+jets"            , nBins , bins.GetArray());         hW_TauFakeUp->SetFillColor(kRed+2);
  TH1F* hW_TauFakeDown        = new TH1F( "hW_TauFakeDown"      ,"W+jets"            , nBins , bins.GetArray());    hW_TauFakeDown->SetFillColor(kRed+2);
  */
  //histograms with fine binning for MSSM
  TH1F* hZttEmb_fb  = new TH1F( "hZttEmb_fb","Embedded"          , 400, 0., 2000.); hZttEmb_fb->SetFillColor(kOrange-4); 
  TH1F* hW_fb        = new TH1F( "hW_fb"      ,"W+jets"            , 400, 0., 2000.); hW_fb->SetFillColor(kRed+2);
  TH1F* hZtt_fb      = new TH1F( "hZtt_fb"    ,"Ztautau"           , 400, 0., 2000.); hZtt_fb->SetFillColor(kOrange-4); 
  TH1F* hTTb_fb      = new TH1F( "hTTb_fb"    ,"ttbar"             , 400, 0., 2000.); hTTb_fb->SetFillColor(kBlue-8);
  TH1F* hVV_fb       = new TH1F( "hVV_fb"     ,"Diboson"           , 400, 0., 2000.); hVV_fb->SetFillColor(kRed+2);
  TH1F* hQCD_fb      = new TH1F( "hQCD_fb"    ,"QCD full vbf"      , 400, 0., 2000.); hQCD_fb->SetFillColor(kMagenta-10);
  TH1F* hZll_fb      = new TH1F( "hZll_fb"    ,"Zll"               , 400, 0., 2000.); hZll_fb->SetFillColor(kBlue-2);
  TH1F* hZj_fb       = new TH1F( "hZj_fb"    ,"Z+jets, jet to tau" , 400, 0., 2000.); hZj_fb->SetFillColor(kBlue-2);
  /*
  TH1F* hSignal[nProd][nMasses];
  for(int iP=0 ; iP<nProd ; iP++) {
    for(int iM=0 ; iM<nMasses ; iM++) {
      hSignal[iP][iM] = new TH1F("h"+nameProd[iP]+nameMasses[iM], nameProd[iP]+nameMasses[iM], nBins , bins.GetArray());
      hSignal[iP][iM]->SetLineWidth(2);
    }
  }
  //GGH Higgs pT weights up/down
  TH1F* hGGFHUp[nMasses]; TH1F* hGGFHDown[nMasses];
  for(int iM=0 ; iM<nMasses ; iM++) {
    hGGFHUp[iM] = new TH1F("hGGFH"+nameMasses[iM]+"Up", "GGFH"+nameMasses[iM]+"Up", nBins , bins.GetArray()); 
    hGGFHUp[iM]->SetLineWidth(2);
    hGGFHDown[iM] = new TH1F("hGGFH"+nameMasses[iM]+"Down", "GGFH"+nameMasses[iM]+"Down", nBins , bins.GetArray());  
    hGGFHDown[iM]->SetLineWidth(2); 
  }
  */
  //GGH Higgs pT weights up/down
  TH1F* hSUSYGGHUp[nMassesS]; TH1F* hSUSYGGHDown[nMassesS];
  for(int iM=0 ; iM<nMassesS ; iM++) {
    hSUSYGGHUp[iM] = new TH1F("hSUSYGGH"+nameMassesS[iM]+"Up", "SUSYGGH"+nameMassesS[iM]+"Up", nBins , bins.GetArray()); 
    hSUSYGGHUp[iM]->SetLineWidth(2);
    hSUSYGGHDown[iM] = new TH1F("hSUSYGGH"+nameMassesS[iM]+"Down", "SUSYGGH"+nameMassesS[iM]+"Down", nBins , bins.GetArray());  
    hSUSYGGHDown[iM]->SetLineWidth(2); 
  }
  /*
  TH1F* hSignalWW[nProdWW][nMassesWW];
  for(int iP=0 ; iP<nProdWW ; iP++) {
    for(int iM=0 ; iM<nMassesWW ; iM++) {
      hSignalWW[iP][iM] = new TH1F("h"+nameProdWW[iP]+nameMassesWW[iM], nameProdWW[iP]+nameMassesWW[iM], nBins , bins.GetArray());
      hSignalWW[iP][iM]->SetLineWidth(2);
    }
  }  
  */
  
  TH1F* hSusy[nProdS][nMassesS];
  
  if(MSSM) {
    for(int iP=0 ; iP<nProdS ; iP++) {
      for(int iM=0 ; iM<nMassesS ; iM++) {
	hSusy[iP][iM] = new TH1F("hSUSY"+nameProdS[iP]+nameMassesS[iM], nameProdS[iP]+nameMassesS[iM], nBins , bins.GetArray());
	hSusy[iP][iM]->SetLineWidth(2);
      }
    }
  }
  
  TH1F* hParameters   = new TH1F( "hParameters", "" ,30, 0, 30);
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

  TString pathToFile    = "/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74_50ns/";
  TString pathToFileSUSY  = "/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74/";
  cout<<"**********************************"<<endl;
  cout<<"**********************************"<<endl;
  cout<<"**********************************"<<endl;
  cout<<"**********************************"<<endl;


  // DATA
  TChain *data = new TChain("outTree");
  data->Add(pathToFile+"/nTupleRun2015B-Data_TauTau.root");
  if(!data) cout << "### DATA NTUPLE NOT FOUND ###" << endl;

  /* // EMBEDDED //
  TString treeEmbedded,fileAnalysisEmbedded;
  if(analysis_.Contains("TauUp") || analysis_.Contains("TauDown") ){
    treeEmbedded = "outTree"+analysis_;
    fileAnalysisEmbedded = analysis_; 
  }
  else {
    treeEmbedded = "outTree";
    fileAnalysisEmbedded = "nominal"; 
  }
  TChain *dataEmbedded = new TChain(treeEmbedded);
  //
  dataEmbedded->Add(pathToFile+"/nTupleRun2012A*Embedded_TauTau_"+fileAnalysisEmbedded+".root");
  */

  // BACKGROUNDS //
  TString treeMC,fileAnalysis;
  if(analysis_.Contains("TauUp") || analysis_.Contains("TauDown") ) {
    treeMC = "outTree"+analysis_;
    fileAnalysis=analysis_;
  }
  else {
    treeMC = "outTree";
    fileAnalysis="nominal";
  }
  //
  TChain *backgroundDY         = new TChain(treeMC);
  TChain *backgroundDYTauTau   = new TChain(treeMC);
  TChain *backgroundDYLtoTau   = new TChain(treeMC);
  TChain *backgroundDYJtoTau   = new TChain(treeMC);
  TChain *backgroundTTbar      = new TChain(treeMC);
  TChain *backgroundOthers     = new TChain(treeMC);
  TChain *backgroundWJets      = new TChain(treeMC);
  //
  
  backgroundDY      ->Add(pathToFile+"nTupleDYJets_TauTau_"+fileAnalysis+".root");

  backgroundTTbar   ->Add(pathToFile+"nTupleTTJets_TauTau_"+fileAnalysis+".root");
  //
  //backgroundOthers  ->Add(pathToFile+"nTupleSAntiTopT_TauTau_"+fileAnalysis+".root");
  backgroundOthers  ->Add(pathToFile+"nTupleSTopT_TauTau_"+fileAnalysis+".root");
  backgroundOthers  ->Add(pathToFile+"nTupleSAntiTopTW_TauTau_"+fileAnalysis+".root");
  backgroundOthers  ->Add(pathToFile+"nTupleSTopTW_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleWWTo2L2Nu_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleWWToLNuQQ_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleWWTo4Q_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleWZTo1L1Nu2Q_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleWZTo3LNu_TauTau_"+fileAnalysis+".root");
  //backgroundOthers  ->Add(pathToFile+"nTupleZZ_TauTau_"+fileAnalysis+".root");
  
  backgroundWJets   ->Add(pathToFile+"nTupleWJets_TauTau_"+fileAnalysis+".root");

  if(!backgroundDY)    cout << "###  NTUPLE DY NOT FOUND ###" << endl;  
  if(!backgroundTTbar) cout << "###  NTUPLE TT NOT FOUND ###" << endl;  
  if(!backgroundOthers)cout << "###  NTUPLE VVt NOT FOUND ###" << endl;  
  if(!backgroundWJets) cout << "###  NTUPLE W NOT FOUND ###" << endl;  
  /*
  TChain *signal[nProd][nMasses];
  for(int iP=0 ; iP<nProd ; iP++) {
    for(int iM=0 ; iM<nMasses ; iM++) {
      signal[iP][iM] = new TChain(treeMC);
      signal[iP][iM]->Add(pathToFile+"/nTuple"+nameProd[iP]+nameMasses[iM]+"_TauTau_"+fileAnalysis+".root");
      if(!signal[iP][iM])cout << "###  NTUPLE Signal " << nameProd[iP]+nameMasses[iM] << " NOT FOUND ###" << endl;  
    }
  }
  TChain *signalWW[nProdWW][nMassesWW];
  for(int iP=0 ; iP<nProdWW ; iP++) {
    for(int iM=0 ; iM<nMassesWW ; iM++) {
      signalWW[iP][iM] = new TChain(treeMC);
      signalWW[iP][iM]->Add(pathToFileHWW+"/nTuple"+nameProdWW[iP]+nameMassesWW[iM]+"_TauTau_"+fileAnalysis+".root");
      if(!signalWW[iP][iM])cout << "###  NTUPLE Signal " << nameProdWW[iP]+nameMassesWW[iM] << " NOT FOUND ###" << endl;
    }
    }*/
  TChain *signalSusy[nProdS][nMassesS];
  if(MSSM) {
    for(int iP=0 ; iP<nProdS ; iP++) {
      for(int iM=0 ; iM<nMassesS ; iM++) {
	signalSusy[iP][iM] = new TChain(treeMC);
	signalSusy[iP][iM]->Add(pathToFileSUSY+"/nTupleSUSY"+nameProdS[iP]+nameMassesS[iM]+"_TauTau_"+fileAnalysis+".root");
      }
    }
  }

  // Split DY into 3 sub-samples (Z->TauTau, Z->ll, ZJ)
  TFile *dummy1;
  if(DOSPLIT) {
    cout << "SPLIT DY SAMPLE ON THE FLY" << endl;
    //
    dummy1 = new TFile("dummy2.root","RECREATE");
    //
    cout << "Now copying g/Z -> tau+ tau- " << endl;
    backgroundDYTauTau = (TChain*)backgroundDY->CopyTree("isZtt==1");                 // g/Z -> tau+ tau-
    //
    cout << "Now copying g/Z -> l+l- l->tau" << endl;
    backgroundDYLtoTau     = (TChain*)backgroundDY->CopyTree("isZll==1"); // g/Z -> ll, l->tau
    //
    cout << "Now copying g/Z, j->tau" << endl;
    backgroundDYJtoTau  = (TChain*)backgroundDY->CopyTree("isZj==1"); // g/Z -> ll, jet->tau
  }
  else {
    cout << "USE DY SEPARATE SUB-SAMPLES" << endl;
    //
    backgroundDYTauTau  ->Add(pathToFile+"/nTupleDYJetsTauTau_TauTau_"  +fileAnalysis+".root");
    
    backgroundDYLtoTau->Add(pathToFile+"/nTupleDYJets*ZLtoTau_TauTau_"    +fileAnalysis+".root");
    
    backgroundDYJtoTau->Add(pathToFile+"/nTupleDYJets*ZJ_TauTau_"    +fileAnalysis+".root");
    
    
    cout << backgroundDYTauTau->GetEntries()  << " come from DY->tautau"         << endl;
    cout << backgroundDYJtoTau->GetEntries()  << " come from DY, jet->tau" << endl;
    cout << backgroundDYLtoTau->GetEntries()  << " come from DY-> LL"         << endl;
  }
  
  cout << backgroundDYTauTau->GetEntries()  << " come from DY->tautau"         << endl;
  cout << backgroundDYJtoTau->GetEntries()  << " come from DY, jet->tau" << endl;
  cout << backgroundDYLtoTau->GetEntries()  << " come from DY-> LL"         << endl;

  // MAPPING //
  if(VERBOSE) cout << "-- gather the trees" << endl;

  const int nVarious = 8;
  //const int nChainsSM = nVarious  + nProd*nMasses;
  //const int nChainsSM = nChainsSM1 + nProdWW*nMassesWW;
  const int nChains  = nVarious + nProdS*nMassesS; //nChainsSM + nProdS*nMassesS;
  TString treeNamesVarious[nVarious]={"SS", "Data", "WJets","TTbar","Others","DYToTauTau", "DYLtoTau", "DYJtoTau"};
  TChain* chainsVarious[nVarious]   ={data, data, backgroundWJets,backgroundTTbar,backgroundOthers, backgroundDYTauTau, backgroundDYLtoTau, backgroundDYJtoTau};
  TString  treeNames[nChains];
  TChain*  chains[nChains];
  mapchain mapAllTrees;

  for(int iCh=0 ; iCh<nChains ; iCh++) {
    
    if(iCh<nVarious) { // fill various tree names and trees
      treeNames[iCh] = treeNamesVarious[iCh];
      chains[iCh]    = chainsVarious[iCh];
    }    
    //else if(iCh<nChainsSM){ // fill signal names and trees
    //  treeNames[iCh] = nameProd[ int((iCh-nVarious)/nMasses) ] + nameMasses[ int((iCh-nVarious)%nMasses) ];
    //  chains[iCh]    = signal[ int((iCh-nVarious)/nMasses) ][ int((iCh-nVarious)%nMasses) ];
    //}
    //else if(iCh<nChainsSM){ // fill signal names and trees
    //  treeNames[iCh] = nameProdWW[ int((iCh-nChainsSM1)/nMassesWW) ] + nameMassesWW[ int((iCh-nChainsSM1)%nMassesWW) ];
    //  chains[iCh]    = signalWW[ int((iCh-nChainsSM1)/nMassesWW) ][ int((iCh-nChainsSM1)%nMassesWW) ];
    //}
    else { // fill signal names and trees
      treeNames[iCh] = "SUSY"+nameProdS[ int((iCh-nVarious)/nMassesS) ] + nameMassesS[ int((iCh-nVarious)%nMassesS) ];
      chains[iCh]    = signalSusy[ int((iCh-nVarious)/nMassesS) ][ int((iCh-nVarious)%nMassesS) ];
    }
    mapAllTrees[ treeNames[iCh] ] = chains[iCh]; // create an entry in the map
  }

  //   if(DEBUG) {
  cout << "######################" << endl
       << "### LIST OF CHAINS ###" << endl;
  for(int iCh=0 ; iCh<nChains ; iCh++) cout << treeNames[iCh] << endl;
  cout << "######################" << endl;
  //   }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ////// TAU PT+ID+ISO //////
  TCut tpt("ptL1>45 && ptL2>45");
  TCut antimu("tightestAntiMu3WPL1>0 && tightestAntiMu3WPL2>0");
  TCut antiele("tightestAntiEMVA5WPL1 > 0 && tightestAntiEMVA5WPL2 > 0");
  TCut tiso("hpsDB3HL1<1.0 && hpsDB3HL2<1.0");
  TCut tdecaymode("decayModeFindingNewDML1>0.5 && decayModeFindingNewDML2>0.5");
  TCut ltiso("hpsDB3HL1<10.0 && hpsDB3HL2<10.0");
  TCut atiso("hpsDB3HL1<10.0 && hpsDB3HL2<10.0 && hpsDB3HL1>1.0 && hpsDB3HL2>1.0 ");
  TCut selection("etaL1<999"); 

  // Choose selection wrt version_
  //chooseSelection(variable_,version_, selection_, tiso, tdecaymode, ltiso, atiso, antimu, antiele, selection,tpt);

  ////// EVENT WISE //////
  TCut lveto="nVetoMuon<=0 && nVetoElectron<=0";
  TCut SS("diTauCharge>0");
  TCut OS("diTauCharge<0");
  
  bool invertDiTauSign = bool(selection_.find("SS")!=string::npos);
  TCut diTauCharge = invertDiTauSign ? SS : OS; 
  string sign = invertDiTauSign ? "SS" : "OS";
  // HLT matching //
  TCut hltevent("HLTx==1 && HLTmatchL1==1 && HLTmatchL2==1");

  ////// CATEGORIES ///
  TCut oneJet("nJets30>=1");
  TCut twoJets("nJets30>=2");
  TCut vbf = TCut("nJets30>=2 && ptj1>30 && ptj2>30 && Mjj>500 && Detajj>3.5 && nVetoJets<=0 && nJets20BTagged<1 && diTauRecoPt>100");
  TCut vbfLoose = TCut("nJets30>=2 && pt1>30 && pt2>30 && Mjj>200 && Deta>2.0 && nVetoJets<=0");

  TCut oneJetHigh("nJets30>0 && ptj1>30 && nJets20BTagged<1 && diTauRecoPt>170");
  TCut oneJetMedium("nJets30>0 && ptj1>30 && nJets20BTagged<1 && diTauRecoPt>100 && diTauRecoPt<170");

  oneJetHigh = oneJetHigh && !vbf;
  oneJetMedium = oneJetMedium && !vbf;

  TCut bTag("nJets30<2 && nJets20BTagged>0");
  TCut bTagLoose("nJets30<2 && nJets20BTaggedLoose>0"); //for W shape in b-Category
  TCut nobTag("nJets20BTagged==0");
  TCut inclusive("");

  TCut novbf("nJets30<1 && nJets20BTagged==0");

  TCut sbinCatIncl("etaL1<999");
  TCut sbinCat("");
  if(     selection_.find("inclusive")!=string::npos) sbinCat = inclusive&&TCut("etaL1<999");
  else if(selection_.find("oneJet")!=string::npos)    sbinCat = oneJet;
  else if(selection_.find("twoJets")!=string::npos)   sbinCat = twoJets;
  else if(selection_.find("oneJetHigh")!=string::npos)sbinCat = oneJetHigh;
  else if(selection_.find("oneJetMedium")!=string::npos)sbinCat = oneJetMedium;
  else if(selection_.find("nobTag")!=string::npos)    sbinCat = nobTag;
  else if(selection_.find("vbf")!=string::npos)       sbinCat = vbf;
  else if(selection_.find("bTag")!=string::npos && selection_.find("nobTag")==string::npos) sbinCat = bTag;


  //TCut sbinPresel           = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && hltevent   ;
  //TCut sbinEmbeddingPresel  = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection               ;
  //TCut sbinLtisoPresel      = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection && hltevent   ;
  //TCut sbinAtisoPresel      = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection && hltevent   ;

  TCut sbinInclusive                     = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;
  TCut sbinChargeRelInclusive            = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection                 && hltevent   ;
  TCut sbinChargeRelLtisoInclusive       = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection                && hltevent   ;
  TCut sbinChargeRelAtisoInclusive       = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection                && hltevent   ;
  TCut sbinEmbeddingInclusive            = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && diTauCharge                ;
  TCut sbinSSInclusive                   = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;

  TCut sbinSSLtisoInclusive              = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;
  TCut sbinSSAtisoInclusive              = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;
  TCut sbinLtisoInclusive                = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;
  TCut sbinAtisoInclusive                = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;

  TCut sbin                   = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;
  TCut sbinEmbedding          = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && diTauCharge                ;
  TCut sbinSS                 = tpt && tiso  && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;
  TCut sbinSSltiso            = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;
  TCut sbinSSatiso            = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection && SS           && hltevent   ;
  TCut sbinLtiso              = tpt && ltiso && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;
  TCut sbinAtiso              = tpt && atiso && tdecaymode && antimu && antiele && lveto && selection && diTauCharge  && hltevent   ;

  cout << sbin << endl;
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  cout << endl;
  cout << "#############################################################" << endl;
  cout << ">>>>>>>>>>> BEGIN Compute inclusive informations <<<<<<<<<<<<" << endl;
  cout << "#############################################################" << endl;
  cout << endl;

  cout << "******** Extrapolation factors for Z->tautau normalization ********" << endl;
  // inclusive DY->tautau:
  TH1F* hExtrap = new TH1F("hExtrap","",nBins , bins.GetArray());
  float Error = 0.;

  //Monte Carlo ZTT Inclusive
  float ExtrapDYInclusive = 0.;
  drawHistogram(sbinCatIncl,"MC", version_,analysis_, mapAllTrees["DYToTauTau"], variable, ExtrapDYInclusive,   Error,   Lumi*lumiCorrFactor/1000., hExtrap, sbinInclusive);
  cout << "All Z->tautau             = " << ExtrapDYInclusive << " +/- " <<  Error << endl; 
  
  //Monte Carlo ZTT in the category
  float ExtrapDYNum = 0.;
  drawHistogram(sbinCat,"MC", version_,analysis_, mapAllTrees["DYToTauTau"], variable, ExtrapDYNum, Error,   Lumi*lumiCorrFactor/1000., hExtrap, sbin);

  //Monte Carlo ZTT: compute the probability for an event to be in a category
  float ExtrapolationFactorMC      = ExtrapDYNum/ExtrapDYInclusive;//probability to enter a category taken from MC
  float ErrorExtrapolationFactorMC = TMath::Sqrt(ExtrapolationFactorMC*(1-ExtrapolationFactorMC)/ExtrapDYInclusive);
  cout << "Extrap. factor using MadGraph            = " << ExtrapolationFactorMC << " +/- " << ErrorExtrapolationFactorMC << endl;

  //Embedded ZTT in the category (Embedding sample not available now)
  float ExtrapEmbedNum = 0.;
  //drawHistogram(sbinEmbeddingPresel,sbinCat,"Embed", version_,analysis_, mapAllTrees["Embedded"], variable, ExtrapEmbedNum, Error, 1.0, hExtrap, sbinEmbedding);

  if(DEBUG) cout << "drawHistogram Embed sbinEmbeddingInclusive" << endl;
  //Embedded ZTT Inclusive
  float ExtrapEmbedDen = 0.;
  //drawHistogram(sbinEmbeddingPresel,sbinCatIncl,"Embed", version_,analysis_, mapAllTrees["Embedded"], variable, ExtrapEmbedDen, Error, 1.0, hExtrap, sbinEmbeddingInclusive);

  //Probability to enter category, from Embedded
  ExtrapolationFactorZ             = ExtrapEmbedDen!=0 ? ExtrapEmbedNum/ExtrapEmbedDen : 0; 
  ErrorExtrapolationFactorZ        = ExtrapEmbedDen!=0 ? TMath::Sqrt(ExtrapolationFactorZ*(1-ExtrapolationFactorZ)/ExtrapEmbedDen) : 0;
  //Ratio of the probabilities from Embedded to MC
  ExtrapolationFactorZDataMC       = ExtrapolationFactorMC!=0 ? ExtrapolationFactorZ/ExtrapolationFactorMC : 0;
  if(!useEmbedding_) {
    ExtrapolationFactorZ = ExtrapolationFactorZDataMC = 1;
  }

  cout << "Extrap. factor using embedded sample     = " << ExtrapolationFactorZ << " +/- " << ErrorExtrapolationFactorZ << endl;
  cout << " ==> data/MC (signal region)             = " << ExtrapolationFactorZDataMC << " +/- " 
       << ExtrapolationFactorZDataMC*(ErrorExtrapolationFactorMC/ExtrapolationFactorMC + 
				      ErrorExtrapolationFactorZ/ExtrapolationFactorZ) 
       << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  cout << endl;
  cout << "******** Extrapolation factors for QCD normalization ********" << endl;
  
  float SSQCDinSignalRegionDATAIncl = 0.; 
  if(invertDiTauSign) OStoSSRatioQCD = 1.0;

  evaluateQCD(mapAllTrees, version_,analysis_, 0, 0, "SS", 
	      SSQCDinSignalRegionDATAIncl, scaleFactorTTSSIncl,
 	      hExtrap, variable,
 	      Lumi/1000,  TTxsectionRatio, lumiCorrFactor,
	      ExtrapolationFactorZDataMC,
 	      JtoTauCorrectionFactor, LtoTauCorrectionFactor,
	      OStoSSRatioQCD,
	      sbinInclusive, sbinCat,
	      sbinSSInclusive, sbinCat,
	      sbinAtisoInclusive,
	      true, true);

  delete hExtrap;

  cout << endl;
  cout << "#############################################################" << endl;
  cout << ">>>>>>>>>>> END Compute inclusive informations <<<<<<<<<<<<<<" << endl;
  cout << "#############################################################" << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
 
  float SSQCDinSignalRegionDATA = 0.; 

  mapchain::iterator it;
  TString currentName, h1Name, h2Name ;
  TChain* currentTree;

  for(int iCh=0 ; iCh<nChains ; iCh++) {

    cout        << endl        << ">>>> Dealing with sample ## " ;
    currentName = treeNames[iCh];
    cout        << currentName << " ##" << endl;

    it = mapAllTrees.find(currentName);
    if(it!=mapAllTrees.end()) currentTree = it->second;
    else {
      cout << "ERROR : sample not found in the map ! continue;" << endl;
      continue;
    }

    if(!currentTree) {
      cout << "ERROR : no such tree" << endl;
      continue;
    }

    if(!MSSM && currentName.Contains("SUSY")) continue;

    h1Name         = "h1_"+currentName;
    h2Name         = "h2_"+currentName;
    TH1F* h1       = new TH1F( h1Name ,"" , nBins , bins.GetArray());
    TH1F* h2       = new TH1F( h2Name ,"" , nBins , bins.GetArray());

    TH1F* hCleaner = new TH1F("hCleaner","",nBins , bins.GetArray());
    TH1F* hCleanerfb = new TH1F("hCleanerfb","",400, 0., 2000.); //fine binning hostogram for MSSM
    if ( !h1->GetSumw2N() )       h1->Sumw2(); 
    if ( !hCleaner->GetSumw2N() ) hCleaner->Sumw2();
    if ( !hCleanerfb->GetSumw2N() ) hCleanerfb->Sumw2();

    if(currentName.Contains("SS")){
      
      cout << "************** BEGIN QCD evaluation using SS events *******************" << endl;

      TH1F* hExtrapSS = new TH1F("hExtrapSS","",nBins , bins.GetArray());
      float dummy1 = 0.;      

      evaluateQCD(mapAllTrees, version_,analysis_, h1, hCleaner, "SS",
		  SSQCDinSignalRegionDATA, scaleFactorTTSS,
		  hExtrapSS, variable,
		  Lumi/1000,  TTxsectionRatio, lumiCorrFactor,
		  ExtrapolationFactorZDataMC,
		  JtoTauCorrectionFactor, LtoTauCorrectionFactor,
		  OStoSSRatioQCD,
		  sbin, sbinCat,
                  sbinSS, sbinCat,
		  sbinAtiso,
		  true, true);

      hExtrapSS->Reset();
      hQCD->Add(h1, 1.0);
      hSS->Add(hCleaner, OStoSSRatioQCD);
      hCleaner->Reset();

      //fine binning for MSSM
      if(selection_.find("nobTag")!=string::npos){
	TH1F* h1fb = new TH1F("h1fb","",400, 0., 2000.); //fine binning hostogram for MSSM 
	if ( !h1fb->GetSumw2N() ) h1fb->Sumw2();
	TH1F* hExtrapSSfb = new TH1F("hExtrapSSfb","",400, 0., 2000.);
	if ( !hExtrapSSfb->GetSumw2N() ) hExtrapSSfb->Sumw2();
	
	evaluateQCD(mapAllTrees, version_,analysis_, h1fb, hCleanerfb, "SS",
		    SSQCDinSignalRegionDATA, scaleFactorTTSS,
		    hExtrapSSfb, variable,
		    Lumi/1000,  TTxsectionRatio, lumiCorrFactor,
		    ExtrapolationFactorZDataMC,
		    JtoTauCorrectionFactor, LtoTauCorrectionFactor,
		    OStoSSRatioQCD,
		    sbin, sbinCat,
		    sbinSS, sbinCat,
		    sbinAtiso,
		    true, true);

	cout << "************** END QCD evaluation using SS events *******************" << endl;
	hQCD_fb->Add(h1fb, h1->Integral()/h1fb->Integral()); 
	hCleanerfb->Reset();
	delete hExtrapSSfb;
      }
      
    }
    else{

      if(!currentName.Contains("Embed")){

	float Error = 0.;

	if(currentName.Contains("DYToTauTau")){
	  float NormDYToTauTau = 0.;
	  drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYToTauTau, Error,   Lumi*lumiCorrFactor/1000., h1, sbin, 1);
	  hZtt->Add(h1, 1.0);

	  //fine binning for MSSM
	  if(selection_.find("bTag")!=string::npos){
	    hCleanerfb->Reset(); float NormDYToTauTau_fb = 0.;
	    drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYToTauTau_fb, Error,  Lumi*lumiCorrFactor/1000., hCleanerfb, sbin, 1);
	    hZtt_fb->Add(hCleanerfb, 1.0);
	    hCleanerfb->Reset();
	  }
	}
	if(currentName.Contains("TTbar")){
	  if(selection_.find("vbf")!=string::npos){ 
	    float NormTTjets = 0.; 
	    drawHistogram(vbfLoose,"MC", version_,analysis_, currentTree, variable, NormTTjets,     Error,   Lumi*TTxsectionRatio/1000., hCleaner, sbinInclusive, 1);

	    hTTb->Add(hCleaner, 1.0);
	    cout << "--- TTbar shape histogram has Integral=" << hTTb->Integral() << endl;
	    NormTTjets = 0.;

	    drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormTTjets,     Error,   Lumi*TTxsectionRatio/1000., h1, sbin, 1);
	    cout << "--- TTbar Norm histogram has Integral=" << h1->Integral() << endl;
	    hTTb->Scale( hTTb->Integral()!=0 ? h1->Integral()/hTTb->Integral() : 1.0 ); 
	  } 
	  else{
	    float NormTTjets = 0.;
	    drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormTTjets,     Error,   Lumi*TTxsectionRatio/1000., h1, sbin, 1);
	    hTTb->Add(h1, 1.0);

	    //fine binning for MSSM 
	    if(selection_.find("bTag")!=string::npos){
	      hCleanerfb->Reset(); float NormTTjets_fb = 0.;
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormTTjets_fb,     Error,  Lumi*TTxsectionRatio/1000., hCleanerfb, sbin, 1);
	      hTTb_fb->Add(hCleanerfb, hTTb->Integral()/hCleanerfb->Integral());
	      hCleanerfb->Reset();
	    }
	  }
	}
	else if(currentName.Contains("WJets")){

	  //For the shape of W in opposite sign region --> taken from MC
	  float NormWJets = 0.;
	  drawHistogram(sbinCat,"MC",version_,analysis_, currentTree, variable, NormWJets, Error,   Lumi/1000., h1, sbin, 1);
	  hW->Add(h1, 1.0);
	  	  
	  //fine binning for MSSM
	  if(selection_.find("bTag")!=string::npos){
	    hCleanerfb->Reset(); float NormWJets_fb = 0.;
	    hW_fb->Reset();
	    drawHistogram(sbinCat,"MC",version_,analysis_, currentTree, variable, NormWJets_fb, Error,   Lumi/1000., hCleanerfb, sbin, 1);
	    hW_fb->Add(hCleanerfb, h1->Integral()/hCleanerfb->Integral());
	  }
	  if(selection_.find("oneJet")!=string::npos){
	    hCleaner->Reset();
	    float NormWJets_ltiso = 0.;
	    hW->Reset();
	    drawHistogram(sbinCat,"MC",version_,analysis_, currentTree, variable, NormWJets_ltiso, Error,   Lumi/1000., hCleaner, sbinLtiso, 1);
	    hW->Add(hCleaner, h1->Integral()/hCleaner->Integral());
	  }
	  else if(selection_.find("vbf")!=string::npos){
	    hCleaner->Reset();
            float NormWJets_ltiso = 0.;
            hW->Reset();
            drawHistogram(vbfLoose,"MC",version_,analysis_, currentTree, variable, NormWJets_ltiso, Error,   Lumi/1000., hCleaner, sbinLtiso, 1);
            hW->Add(hCleaner, h1->Integral()/hCleaner->Integral());
          }
	}
	else if(currentName.Contains("DYJtoTau")){
	  if(selection_.find("vbf")!=string::npos){   
            float NormDYJtoTau = 0.; 
	    hCleaner->Reset();
	    drawHistogram(vbfLoose,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor/1000., hCleaner, sbinInclusive, 1);
	    NormDYJtoTau = 0.; 
	    if(useZDataMC)
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor*ExtrapolationFactorZDataMC/1000., h1, sbin, 1); 
	    else
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor/1000., h1, sbin, 1);
	    hCleaner->Scale( hCleaner->Integral()!=0 ? h1->Integral()/hCleaner->Integral() : 1.0 );
	    if(hCleaner) {
	      hZj->Add(hCleaner, 1.0);
	      hZfakes->Add(hCleaner,1.0);
	      hEWK->Add(hCleaner,1.0); 
	    }
	  }
	  else{
	    float NormDYJtoTau = 0.;
	    if(useZDataMC)
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor*ExtrapolationFactorZDataMC/1000., h1, sbin, 1);
	    else
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor/1000., h1, sbin, 1);
	    hZj->Add(h1, 1.0);
	    hZfakes->Add(h1,1.0);
	    hEWK->Add(h1,1.0);

	    //fine binning for MSSM  
            if(selection_.find("bTag")!=string::npos){ 
              hCleanerfb->Reset();hZj_fb->Reset(); float NormDYJtoTau_fb = 0.; 
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYJtoTau_fb, Error, Lumi*lumiCorrFactor*JtoTauCorrectionFactor/1000., hCleanerfb, sbin, 1);
	      if(hCleanerfb)hZj_fb->Add(hCleanerfb, h1->Integral()/hCleanerfb->Integral());
	    }
	  }
	}
	else if(currentName.Contains("DYLtoTau")){
          if(selection_.find("vbf")!=string::npos){
            float NormDYLtoTau = 0.;
            hCleaner->Reset();
            drawHistogram(vbfLoose,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor/1000., hCleaner, sbinInclusive, 1);
            NormDYLtoTau = 0.;
            if(useZDataMC)
              drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor*ExtrapolationFactorZDataMC/1000., h1, sbin, 1);
            else
              drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor/1000., h1, sbin, 1);
            hCleaner->Scale( hCleaner->Integral()!=0 ? h1->Integral()/hCleaner->Integral() : 1.0 );
            if(hCleaner) {
              hZll->Add(hCleaner, 1.0);
              hZfakes->Add(hCleaner,1.0);
              hEWK->Add(hCleaner,1.0);
            }
          }
          else{
            float NormDYLtoTau = 0.;
            if(useZDataMC)
              drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor*ExtrapolationFactorZDataMC/1000., h1, sbin, 1);
            else
              drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor/1000., h1, sbin, 1);
            hZll->Add(h1, 1.0);
            hZfakes->Add(h1,1.0);
            hEWK->Add(h1,1.0);

            //fine binning for MSSM                                                                                                                                                
            if(selection_.find("bTag")!=string::npos){
              hCleanerfb->Reset();hZll_fb->Reset(); float NormDYLtoTau_fb = 0.;
              drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormDYLtoTau_fb, Error, Lumi*lumiCorrFactor*LtoTauCorrectionFactor/1000., hCleanerfb, sbin, 1);
              if(hCleanerfb)hZll_fb->Add(hCleanerfb, h1->Integral()/hCleanerfb->Integral());
            }
	  }
	}
	else if(currentName.Contains("Others")){
	  if(selection_.find("vbf")!=string::npos){
            float NormOthers = 0.;  
	    drawHistogram(vbfLoose,"MC", version_,analysis_, currentTree, variable, NormOthers, Error,     Lumi/1000., hCleaner, sbinInclusive, 1);
	    NormOthers = 0.; 
            drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormOthers , Error,     Lumi/1000., h1, sbin, 1);
            hCleaner->Scale( hCleaner->Integral()!=0 ? h1->Integral()/hCleaner->Integral() : 1.0 ); 
	    if(hCleaner) {
	      hVV->Add(hCleaner, 1.0);  
	      hEWK->Add(hCleaner,1.0);  
	    }
          } 
	  else{
	    float NormOthers = 0.;
	    drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormOthers , Error,     Lumi/1000., h1, sbin, 1);
	    hVV->Add(h1, 1.0);
	    hEWK->Add(h1,1.0);

	    //fine binning for MSSM    
            if(selection_.find("bTag")!=string::npos){ 
              hCleanerfb->Reset();hVV_fb->Reset(); float NormOthers_fb = 0.;
	      drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormOthers_fb, Error,  Lumi/1000., hCleanerfb, sbin, 1);
              hVV_fb->Add(hCleanerfb, hVV->Integral()/hCleanerfb->Integral());
	    }
	  }
	}
	else if(currentName.Contains("Data")){

	  float NormData = 0.;
	  drawHistogram(sbinCat,"Data", version_,analysis_, currentTree, variable, NormData,  Error, 1.0 , h1, sbin, 1);
	  if(VERBOSE) cout << "DATA YIELD : " << NormData << "  +/- " << Error << endl;
	  hData->Add(h1, 1.0);
	  if ( !hData->GetSumw2N() )hData->Sumw2();
	  
	  if(selection_.find("vbf")!=string::npos){
	    drawHistogram(vbfLoose,"Data", version_,analysis_, currentTree, variable, NormData,  Error, 1.0 , hCleaner, sbinSSatiso ,1);
	    hDataAntiTauIso->Add(hCleaner);
	    //need to subtract other background contributions
	    float NormDYJtoTau = 0.; 
            drawHistogram(vbfLoose,"MC", version_,analysis_, mapAllTrees["DYJtoTau"], variable, NormDYJtoTau, Error,    Lumi*lumiCorrFactor*JtoTauCorrectionFactor/1000., hCleaner, sbinSSatiso, 1);
	    hDataAntiTauIso->Add(hCleaner, -1.0);
	    float NormDYLtoTau = 0.;
            drawHistogram(vbfLoose,"MC", version_,analysis_, mapAllTrees["DYLtoTau"], variable, NormDYLtoTau, Error,    Lumi*lumiCorrFactor*LtoTauCorrectionFactor/1000., hCleaner, sbinSSatiso, 1);
            hDataAntiTauIso->Add(hCleaner, -1.0);
	    float NormTTjets = 0.; 
            drawHistogram(vbfLoose,"MC", version_,analysis_, mapAllTrees["TTbar"], variable, NormTTjets,     Error,   Lumi*TTxsectionRatio*scaleFactorTTSS/1000., hCleaner, sbinSSatiso, 1);
	    hDataAntiTauIso->Add(hCleaner, -1.0);
	    float NormDYtoTauTau = 0.;
            drawHistogram(vbfLoose,"MC", version_,analysis_, mapAllTrees["DYtoTauTau"], variable, NormDYtoTauTau, Error,    Lumi*lumiCorrFactor/1000., hCleaner, sbinSSatiso, 1);
            hDataAntiTauIso->Add(hCleaner, -1.0);
            float NormWjets = 0.;
            drawHistogram(vbfLoose,"MC", version_,analysis_, mapAllTrees["WJet"], variable, NormWjets,     Error,   Lumi/1000., hCleaner, sbinSSatiso, 1);
            hDataAntiTauIso->Add(hCleaner, -1.0);

	    hDataAntiTauIsoQCD->Add(hDataAntiTauIso, hQCD->Integral()/hDataAntiTauIso->Integral());
	  }

	} 
	else if(currentName.Contains("VBFH")  || 
		currentName.Contains("GGFH")  ||
		currentName.Contains("VH")    ||
		currentName.Contains("SUSY")){

	  float NormSign = 0.;
	  if(VERBOSE) cout << "SIGNAL " << currentName << " ENTRIES " << currentTree->GetEntries() << endl;
	  drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormSign, Error,    Lumi/1000., h1, sbin, 1);

	  if(currentName.Contains("VBFH"+TmH_)){
	    hSgn1->Add(h1,1.0);
	    hSgn1->Scale(magnifySgn_);
	    hSgn->Add(hSgn1,1.0);
	  }
	  else if(currentName.Contains("GGFH"+TmH_)){
	    hSgn2->Add(h1,1.0);
	    hSgn2->Scale(magnifySgn_);
	    hSgn->Add(hSgn2,1.0);
	  }
	  else  if(currentName.Contains("VH"+TmH_)){
	    hSgn3->Add(h1,1.0);
	    hSgn3->Scale(magnifySgn_);
	    hSgn->Add(hSgn3,1.0);
	  }	   
	  /*
	  for(int iP=0 ; iP<nProd ; iP++)
	    for(int iM=0 ; iM<nMasses ; iM++)
	      if(currentName.Contains(nameProd[iP]+nameMasses[iM]))
		hSignal[iP][iM]->Add(h1,1.0);
	  //to be used when weight is available
	  for(int iM=0 ; iM<nMasses ; iM++){
	    if(currentName.Contains("GGFH"+nameMasses[iM])){
	      hCleaner->Reset(); float NormSignUp = 0.;
	      drawHistogram(sbinCat,"MCGGFHUp", version_, analysis_,currentTree, variable, NormSignUp, Error,   Lumi/1000., hCleaner, sbin, 1);
	      hGGFHUp[iM]->Add(hCleaner,1.0);
	      hCleaner->Reset(); float NormSignDown = 0.;
	      drawHistogram(sbinCat,"MCGGFHDown", version_, analysis_,currentTree, variable, NormSignDown, Error,   Lumi/1000., hCleaner, sbin, 1);
	      hGGFHDown[iM]->Add(hCleaner,1.0);
	    }
	  }
	 
	  for(int iP=0 ; iP<nProdWW ; iP++)
            for(int iM=0 ; iM<nMassesWW ; iM++)
              if(currentName.Contains(nameProdWW[iP]+nameMassesWW[iM]))
                hSignalWW[iP][iM]->Add(h1,1.0);
	  */
	  if(MSSM) {
	    if(currentName.Contains("SUSY"))
	      {
		/* //select events within 30% of Higgs mass
		TString sampleName = currentName;
		if(sampleName.Contains("SUSYGGH"))sampleName.ReplaceAll("SUSYGGH", "");
		else if(sampleName.Contains("SUSYBBH"))sampleName.ReplaceAll("SUSYBBH", "");
		float mA = atof(sampleName.Data());
		TCut HWidth(Form("genVMass > 0.7*%f && genVMass < 1.3*%f", mA, mA));  
		*/
		TCut HWidth("");
		float NormSign = 0.; 
		drawHistogram(sbinCat,"MC", version_,analysis_, currentTree, variable, NormSign, Error,    Lumi/1000., h1, (sbin&&HWidth), 1);

		for(int iP=0 ; iP<nProdS ; iP++)
		  {
		    for(int iM=0 ; iM<nMassesS ; iM++)
		      {			
			TString ProcessName("SUSY"+nameProdS[iP]+nameMassesS[iM]);
			if(currentName==ProcessName)
			  {
			    hSusy[iP][iM]->Add(h1,1.0);
			  }
		      }
		  }
		/*
		for(int iM=0 ; iM<nMassesS ; iM++)
		  {
		    TString ProcessName("SUSYGGH"+nameMassesS[iM]);
		    if(currentName==ProcessName)
		      {
			h1->Reset(); float NormSignUp = 0.;
			drawHistogram(sbinCat,"MCSUSYGGHUp", version_, analysis_,currentTree, variable, NormSignUp, Error,   Lumi/1000., h1, (sbin&&HWidth), 1);
			hSUSYGGHUp[iM]->Add(h1,1.0);

			h1->Reset(); float NormSignDown = 0.;
			drawHistogram(sbinCat,"MCSUSYGGHDown", version_, analysis_,currentTree, variable, NormSignDown, Error,   Lumi/1000., h1, (sbin&&HWidth), 1);
			hSUSYGGHDown[iM]->Add(h1,1.0);

		      }//End SUSY GGH
		  }//End Masses
		*/
	      }//End SUSY
	  }//End MSSM
	}//End MC Signal
      }
      else{
	if(selection_.find("vbf")!=string::npos){
	  //Moriond method
	  float NormEmbed = 0.; 
	  drawHistogram(sbinCat,"Embed", version_,analysis_, currentTree, variable, NormEmbed,  Error, 1.0 , h1,  sbinEmbedding  ,1); 
	  h1->Scale( (ExtrapolationFactorZ*ExtrapolationFactorZFromSideband)/h1->Integral()); 
	  hZttEmb->Add(h1, 1.0); 
	}
	else{
	  float NormEmbed = 0.;
	  drawHistogram(sbinCat,"Embed", version_,analysis_, currentTree, variable, NormEmbed,  Error, 1.0 , h1,  sbinEmbedding  ,1);
	  h1->Scale( (ExtrapolationFactorZ*ExtrapolationFactorZFromSideband)/h1->Integral());
	  hZttEmb->Add(h1, 1.0);
	  
	  //fine binning for MSSM
	  if(selection_.find("bTag")!=string::npos){ 
	    hCleanerfb->Reset(); float NormEmbed_fb = 0.;
	    drawHistogram(sbinCat,"Embed", version_,analysis_, currentTree, variable, NormEmbed_fb, Error, 1.0 , hCleanerfb,  sbinEmbedding  ,1);
	    hZttEmb_fb->Add(hCleanerfb, (ExtrapolationFactorZ*ExtrapolationFactorZFromSideband));
	  }
	}
      }
    }
  
    /////////////////////////////////////////////////////////////////////////////////////

    if(VERBOSE) cout<<(it->first) << " ==> " 
		    << h1->Integral() << " +/- " 
		    << TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries())
		    << endl;

    //out << (it->first) << "  " << h1->Integral() << " $\\pm$ " <<  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()) << endl;
    char* c = new char[50];
    if(h1->Integral()>=10) 
      sprintf(c,"$%.0f\\pm%.0f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=1)
      sprintf(c,"$%.1f\\pm%.1f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.1)
      sprintf(c,"$%.2f\\pm%.2f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.01)
      sprintf(c,"$%.3f\\pm%.3f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else
      sprintf(c,"$%.5f\\pm%.5f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    out << string(c) << "  //" << (it->first) << endl;
    delete c;

    delete hCleaner; delete hCleanerfb;
  }

  cout << endl;
  cout << "All samples done. Filling hParameters..." << endl;
  hParameters->SetBinContent(1, ExtrapolationFactorZ);               hParameters->GetXaxis()->SetBinLabel(1,"ExtrapolationFactorZ");
  hParameters->SetBinContent(2, ErrorExtrapolationFactorZ);          hParameters->GetXaxis()->SetBinLabel(2,"ErrorExtrapolationFactorZ");
  hParameters->SetBinContent(3, ExtrapolationFactorZDataMC);         hParameters->GetXaxis()->SetBinLabel(3,"ExtrapolationFactorZDataMC");
  hParameters->SetBinContent(4, ExtrapolationFactorZFromSideband);   hParameters->GetXaxis()->SetBinLabel(4,"ExtrapolationFactorZFromSideband");
  hParameters->SetBinContent(5, ExtrapolationFactorSidebandZDataMC); hParameters->GetXaxis()->SetBinLabel(5,"ExtrapolationFactorSidebandZDataMC");
  hParameters->SetBinContent(6, SSQCDinSignalRegionDATAIncl);        hParameters->GetXaxis()->SetBinLabel(9,"SSQCDinSignalRegionDATAIncl");
  hParameters->SetBinContent(7,SSQCDinSignalRegionDATA);            hParameters->GetXaxis()->SetBinLabel(13,"SSQCDinSignalRegionDATA");
  hParameters->SetBinContent(8,scaleFactorTTOS);                    hParameters->GetXaxis()->SetBinLabel(20,"scaleFactorTTOS");
  hParameters->SetBinContent(9,scaleFactorTTSS);                    hParameters->GetXaxis()->SetBinLabel(21,"scaleFactorTTSS");
  hParameters->SetBinContent(10,scaleFactorTTSSIncl);                hParameters->GetXaxis()->SetBinLabel(22,"scaleFactorTTSSIncl");
  hParameters->SetBinContent(11,SSIsoToSSAIsoRatioQCD);              hParameters->GetXaxis()->SetBinLabel(23,"SSIsoToSSAIsoRatioQCD");
  hParameters->SetBinContent(12,ExtrapDYInclusive);                  hParameters->GetXaxis()->SetBinLabel(24,"ExtrapDYInclusive");
  hParameters->SetBinContent(13,ExtrapolationFactorMC);             hParameters->GetXaxis()->SetBinLabel(28,"ExtrapolationFactorMC");
  hParameters->SetBinContent(14,ErrorExtrapolationFactorMC);        hParameters->GetXaxis()->SetBinLabel(29,"ErrorExtrapolationFactorMC");

  hParameters->GetXaxis()->LabelsOption("v");

  //YIELDS
  if(selection_.find("vbf")!=string::npos){
    out<<"Yields for VBF :"<<endl;
    out<<"VBF data : hData -> "<<hData->Integral()<<endl;
    //out<<"VBF Ztt : hZttEmb -> "<<hZttEmb->Integral()<<endl;
    out<<"VBF Ztt : hZtt -> "<<hZtt->Integral()<<endl;
    out<<"VBF QCD : hDataAntiTauIsoQCD -> "<<hDataAntiTauIsoQCD->Integral()<<endl;
    out<<"VBF Z, j->t : hZj -> "<<hZj->Integral()<<endl;
    out<<"VBF Z, l->t : hZll -> "<<hZll->Integral()<<endl;
    out<<"VBF TTb : hTTb -> "<<hTTb->Integral()<<endl;
    out<<"VBF VV : hVV -> "<<hVV->Integral()<<endl;
  }
  else if(selection_.find("oneJet")!=string::npos){
    out<<"Yields for oneJet :"<<endl;
    out<<"oneJet data : hData -> "<<hData->Integral()<<endl;
    //out<<"oneJet Ztt : hZttEmb -> "<<hZttEmb->Integral()<<endl;
    out<<"VBF Ztt : hZtt -> "<<hZtt->Integral()<<endl;
    out<<"oneJet QCD : hQCD -> "<<hQCD->Integral()<<endl;
    out<<"oneJet W : hW -> "<<hW->Integral()<<endl;
    out<<"oneJet Z, j->t : hZj -> "<<hZj->Integral()<<endl;
    out<<"VBF Z, l->t : hZll -> "<<hZll->Integral()<<endl;
    out<<"oneJet TTb : hTTb -> "<<hTTb->Integral()<<endl;
    out<<"oneJet VV : hVV -> "<<hVV->Integral()<<endl;
  }
  else if(selection_.find("bTag")!=string::npos){
    out<<"Yields for bTag :"<<endl;
    out<<"bTag data : hData -> "<<hData->Integral()<<endl;
    //out<<"bTag Ztt : hZttEmb -> "<<hZttEmb->Integral()<<endl;
    out<<"VBF Ztt : hZtt -> "<<hZtt->Integral()<<endl;
    out<<"bTag QCD : hQCD -> "<<hQCD->Integral()<<endl;
    out<<"bTag W : hW -> "<<hW->Integral()<<endl;
    out<<"bTag Z, j->t : hZj -> "<<hZj->Integral()<<endl;
    out<<"VBF Z, l->t : hZll -> "<<hZll->Integral()<<endl;
    out<<"bTag TTb : hTTb -> "<<hTTb->Integral()<<endl;
    out<<"bTag VV : hVV -> "<<hVV->Integral()<<endl;
  }
  else{
    out<<"Yields for "<<selection_<<" : "<<endl;
    out<<selection_<<" data : hData -> "<<hData->Integral()<<endl;
    //out<<selection_<<" Ztt : hZttEmb -> "<<hZttEmb->Integral()<<endl;
    out<<"VBF Ztt : hZtt -> "<<hZtt->Integral()<<endl;
    out<<selection_<<" QCD : hQCD -> "<<hQCD->Integral()<<endl;
    out<<selection_<<" W : hW -> "<<hW->Integral()<<endl;
    out<<selection_<<" Z, j->t : hZj -> "<<hZj->Integral()<<endl;
    out<<"VBF Z, l->t : hZll -> "<<hZll->Integral()<<endl;
    out<<selection_<<" TTb : hTTb -> "<<hTTb->Integral()<<endl;
    out<<selection_<<" VV : hVV -> "<<hVV->Integral()<<endl;
  }
  out.close();

  if(scaleByBinWidth && variable_.Contains("diTauNSVfitMass")){ 
    hData->Scale(1.0, "width");
    hTTb->Scale(1.0, "width");
    hZttEmb->Scale(1.0, "width");
    hZtt->Scale(1.0, "width");
    hDataAntiTauIsoQCD->Scale(1.0, "width");
    hQCD->Scale(1.0, "width");
    hSS->Scale(1.0, "width");
    hEWK->Scale(1.0, "width");
    hSgn->Scale(1.0, "width");
  }

  hSiml->Add(hTTb,1.0);
  if(useEmbedding_)
    hSiml->Add(hZttEmb,1.0);
  else 
    hSiml->Add(hZtt);

  if(selection_.find("vbf")!=string::npos)
    hSiml->Add(hDataAntiTauIsoQCD,1.0);
  else
    hSiml->Add(hSS,1.0); //hSiml->Add(hQCD,1.0);

  //VV + W + ZJ
  hSiml->Add(hEWK,1.0);

  if(selection_.find("vbf")!=string::npos)
    aStack->Add(hDataAntiTauIsoQCD);
  else
    aStack->Add(hSS); //aStack->Add(hQCD);
  
  // TT
  aStack->Add(hTTb);

  // VV + W + ZJ
  aStack->Add(hEWK);

  // ZTT
  if(useEmbedding_)
    aStack->Add(hZttEmb);
  else
    aStack->Add(hZtt);

  if(!logy_ && !MSSM)
    aStack->Add(hSgn);

  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hSgn,Form("(%.0fx) H#rightarrow#tau#tau m_{H}=%d",magnifySgn_,mH_),"F");
  
  if(useEmbedding_)
    leg->AddEntry(hZttEmb,"Z#rightarrow#tau#tau (embedded)","F");
  else
    leg->AddEntry(hZtt,"Z#rightarrow#tau#tau","F"); 
  leg->AddEntry(hEWK,"Electroweak","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTTb,"t#bar{t}","F");
  
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  if(scaleByBinWidth && variable_.Contains("diTauNSVfitMass")){
    hStack->SetYTitle(Form(" Events/(%.0f %s)", 1.0, Unities_.Data() ) );
    hData->SetYTitle(Form(" Events/(%.0f %s)", 1.0, Unities_.Data() ) );
  }
  else 
    hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  if(!logy_)
    hData->SetAxisRange(0.0, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  else
    hData->SetAxisRange(0.1, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  if(logy_ && !MSSM)
    hSgn->Draw("HISTSAME");

  leg->Draw();

  pad2->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

   TH1F* hRatio = new TH1F( "hRatio" ," ; ; #frac{(DATA-MC)}{#sqrt{DATA}}" , nBins , bins.GetArray());
   hRatio->Reset();
   hRatio->SetXTitle("");
   hRatio->SetYTitle("#frac{(DATA-MC)}{MC}");

   hRatio->SetMarkerStyle(kFullCircle);
   hRatio->SetMarkerSize(0.8);
   hRatio->SetLabelSize(0.12,"X");
   hRatio->SetLabelSize(0.10,"Y");
   hRatio->SetTitleSize(0.12,"Y");
   hRatio->SetTitleOffset(0.36,"Y");

   float maxPull = 0.;
   for(int k = 0 ; k < hRatio->GetNbinsX(); k++){
     float pull = hData->GetBinContent(k) - hSiml->GetBinContent(k);
     if(hSiml->GetBinContent(k)>0)
       pull /= hSiml->GetBinContent(k);
     hRatio->SetBinContent(k, pull);
     if(TMath::Abs(pull) > maxPull)
       maxPull = TMath::Abs(pull);
   }
   hRatio->SetAxisRange(-1.2*maxPull,1.2*maxPull,"Y");
   hRatio->Draw("P");

   TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hStack->GetXaxis()->GetXmax());
   line->SetLineStyle(3);
   line->SetLineWidth(1.5);
   line->SetLineColor(kBlack);
   line->Draw("SAME");

  //return;

  if(logy_){
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_log.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_log.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_log.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
  }
  else{
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
  }

  TH1F* hDataBlind  = 0; 
  TH1F* hRatioBlind = 0; 

  // Plot blinded histogram
  if(variable_.Contains("Mass")) {

    if(MSSM) {
      hDataBlind  = blindHistogram(hData,  100, 2000, "hDataBlind");
      hRatioBlind = blindHistogram(hRatio, 100, 2000, "hRatioBlind");
    }
    else{
      hDataBlind  = blindHistogram(hData,  100, 160, "hDataBlind");
      hRatioBlind = blindHistogram(hRatio, 100, 160, "hRatioBlind");
    }

    c1 = new TCanvas("c2","",5,30,650,600);
    c1->SetGrid(0,0);
    c1->SetFillStyle(4000);
    c1->SetFillColor(10);
    c1->SetTicky();
    c1->SetObjectStat(0);
    c1->SetLogy(logy_);

    pad1 = new TPad("pad2_1DEta","",0.05,0.22,0.96,0.97);
    pad2 = new TPad("pad2_2DEta","",0.05,0.02,0.96,0.20);
    
    pad1->SetFillColor(0);
    pad2->SetFillColor(0);
    pad1->Draw();
    pad2->Draw();
    
    pad1->cd();
    pad1->SetLogy(logy_);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleH(0.07);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleOffset(1.3,"y");

    hDataBlind->Draw("P");
    aStack->Draw("HISTSAME");
    hDataBlind->Draw("PSAME");
    if(logy_ && !MSSM) hSgn->Draw("HISTSAME");
    //if(logy_ && MSSM) hSgnSUSY->Draw("HISTSAME");
    leg->Draw();

    pad2->cd();
    gStyle->SetOptStat(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleH(0.07);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleOffset(1.3,"y");

    hRatioBlind->Draw("P");  

    if(logy_){
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind_log.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind_log.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind_log.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    }
    else{
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
      c1->SaveAs(Form(location+"/%s/plots/plot_tauTau_mH%d_%s_%s_%s_blind.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()));
    }
  }

  // templates for fitting
  TFile* fout = new TFile(Form(location+"/%s/histograms/tauTau_mH%d_%s_%s_%s.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.Data(),variable_.Data()),"RECREATE");
  fout->cd();

  hSiml->Write();
  hQCD->Write();
  hSS->Write();
  hZj->Write();
  hZfakes->Write();
  hTTb->Write();
  hTTbUp->Write();
  hTTbDown->Write();
  hZtt->Write();
  hZttEmb->Write();
  hW->Write();
  hVV->Write();
  hEWK->Write();
  hSgn1->Write();
  hSgn2->Write();
  hSgn3->Write();
  hDataAntiTauIsoQCD->Write();
  hData->Write();
  hParameters->Write();

  //fine binning histograms for MSSM
  hZttEmb_fb->Write(); hW_fb->Write(); hZtt_fb->Write(); 
  hZj_fb->Write(); hTTb_fb->Write(); 
  hQCD_fb->Write(); hVV_fb->Write(); 
  /*
  for(int iP=0 ; iP<nProd ; iP++)
    for(int iM=0 ; iM<nMasses ; iM++)
      if(hSignal[iP][iM]) hSignal[iP][iM]->Write();
  for(int iM=0 ; iM<nMasses ; iM++){
    hGGFHUp[iM]->Write();
    hGGFHDown[iM]->Write();
  }
  for(int iP=0 ; iP<nProdWW ; iP++)
    for(int iM=0 ; iM<nMassesWW ; iM++)
      if(hSignalWW[iP][iM]) hSignalWW[iP][iM]->Write();
  */
  if(MSSM) {
    for(int iP=0 ; iP<nProdS ; iP++)
      for(int iM=0 ; iM<nMassesS ; iM++)
	{
	  if(hSusy[iP][iM]) hSusy[iP][iM]->Write(); 
	}

    for(int iM=0 ; iM<nMassesS ; iM++){
      hSUSYGGHUp[iM]->Write();
      hSUSYGGHDown[iM]->Write();
    }
  }

  if(variable_.Contains("Mass")) hDataBlind->Write();

  fout->Write();
  fout->Close();
  std::cout<<" before deleting histograms"<<std::endl;
  delete hQCD; delete hSS; delete hZj; delete hZll; delete hZfakes; delete hTTb; delete hTTbUp; delete hTTbDown; delete hZtt; 
  delete hW; delete hVV; delete hSgn; delete hSgn1; delete hSgn2; delete hSgn3; delete hData; delete hParameters;
  delete hDataAntiTauIsoQCD; delete hDataAntiTauIso;

  delete hZttEmb_fb; delete hZtt_fb; delete hW_fb; delete hZj_fb; delete hZll_fb; delete hTTb_fb;
  delete hQCD_fb; delete hVV_fb; 
  
  /*
  for(int iP=0 ; iP<nProd ; iP++)
    for(int iM=0 ; iM<nMasses ; iM++)
      if(hSignal[iP][iM]) delete hSignal[iP][iM];
  for(int iM=0 ; iM<nMasses ; iM++){
    delete hGGFHUp[iM]; delete hGGFHDown[iM];
    }*/
  for(int iM=0 ; iM<nMassesS ; iM++){
    delete hSUSYGGHUp[iM]; delete hSUSYGGHDown[iM];
  }
  /*for(int iP=0 ; iP<nProdWW ; iP++)
    for(int iM=0 ; iM<nMassesWW ; iM++)
      if(hSignalWW[iP][iM]) delete hSignalWW[iP][iM];
  */
  if(MSSM) {
    for(int iP=0 ; iP<nProdS ; iP++)
      for(int iM=0 ; iM<nMassesS ; iM++)
	{
	  if(hSusy[iP][iM]) delete hSusy[iP][iM];
	}
  }

  delete aStack;  delete hEWK; delete hSiml; delete hZttEmb;  delete hRatio; delete line;
  delete fout;

  /*
  for(int iP=0 ; iP<nProd ; iP++) {
    for(int iM=0 ; iM<nMasses ; iM++) {
      //signal[iP][iM]->Close();
      delete signal[iP][iM];
    }
  }
  for(int iP=0 ; iP<nProdWW ; iP++) {
    for(int iM=0 ; iM<nMassesWW ; iM++) {
      delete signalWW[iP][iM];
    }
    }*/
  if(MSSM){
    for(int iP=0 ; iP<nProdS ; iP++) {
      for(int iM=0 ; iM<nMassesS ; iM++) {
	delete signalSusy[iP][iM];
      }
    }
  }

  delete data; 
  delete backgroundDY;      
  delete backgroundDYTauTau;
  delete backgroundDYLtoTau;
  delete backgroundDYJtoTau;
  delete backgroundTTbar;   
  delete backgroundOthers;  
  delete backgroundWJets;   
  std::cout<<" after deleting histograms"<<std::endl;  
}


///////////////////////////////////////////////////////////////////////////////////////////////



void plotTauTauAll( Int_t useEmbedded = 1, TString outputDir = "DiTauV1"){
      
  vector<string> variables;
  vector<int> mH;

  variables.push_back("diTauVisMass");
  variables.push_back("diTauNSVfitMass");
 
  mH.push_back(125);
  
  //plotTauTau(125,0,"inclusive",""   ,"decayModeL1",     "#tau_{1} decay mode","units"   ,outputDir,10,0,10, 5.0,0,1.4);
  //plotTauTau(125,0,"inclusive",""   ,"visibleTauMassL1","visible #tau_{1} mass","GeV"   ,outputDir,40,0,2,5.0,0,1.2);  
  //plotTauTau(125,0,"inclusive",""   ,"decayModeL2",     "#tau_{2} decay mode","units"   ,outputDir,10,0,10, 5.0,0,1.4);
  //plotTauTau(125,0,"inclusive",""   ,"visibleTauMassL2","visible #tau_{2} mass","GeV"   ,outputDir,40,0,2,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MEtMVA","E_{T}^{miss} MVA","GeV"                        ,outputDir,40,0,100,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MEtMVAPhi","E_{T}^{miss} MVA #phi","units"              ,outputDir,32,-3.2,3.2,   5.0,0,1.5);
  //plotTauTau(125,0,"inclusive",""   ,"MtLeg1MVA","M_{T}(#tau_{1}#nu) MVA","GeV" ,                  outputDir,40,0,160,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MtLeg2MVA","M_{T}(#tau_{2}#nu) MVA","GeV" ,                  outputDir,40,0,160,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MEt","E_{T}^{miss} ","GeV"                        ,outputDir,40,0,100,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MEtPhi","E_{T}^{miss}  #phi","units"              ,outputDir,32,-3.2,3.2,   5.0,0,1.5);
  //plotTauTau(125,0,"inclusive",""   ,"MtLeg1","M_{T}(#tau_{1}#nu) ","GeV" ,                  outputDir,40,0,160,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"MtLeg2","M_{T}(#tau_{2}#nu) ","GeV" ,                  outputDir,40,0,160,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"diTauVisMass","visible mass","GeV"      ,outputDir,50,0,200,5.0,0,1.2); 
  ////plotTauTau(125,0,"inclusive",""   ,"diTauNSVfitMass","SVfit mass","GeV"     ,outputDir,60,0,360,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"etaL1","#tau_{1} #eta", "units"              ,outputDir,25,-2.5, 2.5,5.0,0,2.);
  //plotTauTau(125,0,"inclusive",""   ,"ptL1","#tau_{1} p_{T}", "GeV"                ,outputDir,20,40, 140,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"ptL2","#tau_{2} p_{T}","GeV"           ,outputDir,20,40, 140,5.0,0,1.2);
  //plotTauTau(125,0,"inclusive",""   ,"etaL2","#tau_{2} #eta","units"         ,outputDir,25,-2.5, 2.5,5.0,0,2.);
  //plotTauTau(125,0,"inclusive",""   ,"numPV","reconstructed vertexes","units"             ,outputDir,30,0,30,5.0,0,1.5);
  //plotTauTau(125,0,"inclusive",""   ,"nJets30","jet multiplicity","units"                 ,outputDir,10,0, 10,5.0,1,10);
  //plotTauTau(125,0,"oneJet",""        ,"ptj1", "leading jet p_{T}","GeV"       ,outputDir,50,30, 330,5.0,1,100);
  //plotTauTau(125,0,"oneJet",""        ,"etaj1","leading jet #eta","units"      ,outputDir,21,-5, 5,5.0,0,2.);
  plotTauTau(125,0,"bTag",""        ,"ptB1", "leading b-tagged jet p_{T}","GeV"       ,outputDir,50,30, 330,5.0,1,10);
  //plotTauTau(125,0,"bTag",""        ,"etaB1","leading b-tagged jet #eta","units"      ,outputDir,21,-5, 5,5.0,0,2.);
  
  return;
  /*
  for(unsigned int i = 0 ; i < variables.size(); i++){
    for(unsigned j = 0; j < mH.size(); j++){

      plotTauTau(mH[j],useEmbedded,"inclusive",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"inclusive","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"inclusive","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);

      plotTauTau(mH[j],useEmbedded,"oneJetLow",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"oneJetLow","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"oneJetLow","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      
      plotTauTau(mH[j],useEmbedded,"oneJetHigh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"oneJetHigh","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      plotTauTau(mH[j],useEmbedded,"oneJetHigh","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2);
      
      plotTauTau(mH[j],useEmbedded,"vbf",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2); 
      plotTauTau(mH[j],useEmbedded,"vbf","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2); 
      plotTauTau(mH[j],useEmbedded,"vbf","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,1.2); 
    }
  }
  
  return;
  */
}

int main(int argc, const char* argv[])
{

  std::cout << "plotTauTau()" << std::endl;
  gROOT->SetBatch(true);

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  int mH, nBins, logy, useEmb; 
  float magnify, hltEff, xMin, xMax, maxY;
  float antiWsgn, antiWsdb;
  string category, analysis, variable, xtitle, unity, outputDir, version;

  if(argc==1) plotTauTauAll();
  //   else if(argc==20) { 
  else if(argc>14) { 

    mH         =  (int)atof(argv[1]); 
    category   =  argv[2]; 
    variable   =  argv[3]; 
    xtitle     =  argv[4]; 
    unity      =  argv[5]; 
    nBins      =  (int)atof(argv[6]); 
    xMin       =  atof(argv[7]); 
    xMax       =  atof(argv[8]); 
    magnify    =  atof(argv[9]); 
    logy       =  (int)atof(argv[10]); 
    maxY       =  atof(argv[11]) ;
    outputDir  =  argv[12]; 
    version    =  argv[13];
    useEmb     =  (int)atof(argv[14]);
    analysis   =  argc>15 ? argv[15] : "";
    //     antiWsgn   =  atof(argv[17]);
    //     antiWsdb   =  atof(argv[18]);
    cout << endl << " VERSION : " << version << " | ANALYSIS : " << analysis << endl << endl;

    //     plotTauTau(mH,useEmb,category,analysis,variable,xtitle,unity,outputDir,nBins,xMin,xMax,magnify,hltEff,logy,maxY, version,antiWsgn,antiWsdb);
    plotTauTau(mH,useEmb,category,analysis,variable,xtitle,unity,outputDir,nBins,xMin,xMax,magnify,logy,maxY, version);
  }
  else { cout << "Please put at least 14 arguments" << endl; return 1;}

  cout << "DONE" << endl;
  return 0;
}

