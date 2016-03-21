#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TSystem.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#define RESCALETO1PB false
#define DOSUSY false
#define OldCat false

using namespace std;

void checkValidity(TH1F* h){
  int nBins = h->GetNbinsX();
  for(int i=1 ; i<=nBins+1 ; i++) {
    //f << h->GetBinContent(i) << " | " ;
    if(h->GetBinContent(i) < 0) {
      //h->SetBinContent(i,0);
      h->SetBinContent(i,1e-9);
    }
    if(h->GetBinContent(i) == 0) h->SetBinContent(i,1e-9);
  }
}

float higgsXsection(int mH = 125, string process = "ggH"){

  float xsection = 1.0;

  if(process.find("ggH")!=string::npos){
    if(mH == 125) xsection = 43.92*0.0632 ;
  }
  if(process.find("qqH")!=string::npos){
    if(mH == 125) xsection = 3.748*0.0632; 
  }
  if(process.find("VH")!=string::npos){
    if(mH == 125) xsection = 0.5085*0.0632; // only ttH
  }

  return xsection;
}

void produce(
	     int mH_=125,
	     string variable_ = "diTauVisMass",
	     string analysis_ = "",
	     string bin_ = "inclusive",
	     TString outputDir = "DiTauV1",
	     int useEmb = 1,
	     TString location = "/nfs/dust/cms/user/anayak/CMS/OnSLC6/CMSSW_7414_htt/src/DesyTauAnalyses/NTupleMaker/bin/results/"
	     ){
  
  cout << "Now doing mass mH=" << mH_ << ", for variable " << variable_ << " analysis " << analysis_ << " and bin " << bin_ << endl;
  TFile* fin = 0;
  string analysisFile = analysis_;

  fin = new TFile(Form(location+"/%s/histograms/tauTau_mH%d_%s_%s_%s.root", outputDir.Data(), 125, bin_.c_str() , analysisFile.c_str(), variable_.c_str()), "READ");
  TString nameZtt="hDataEmb";
  if(!useEmb) nameZtt="hZtt";
  
  //TFile* fin_tUp = new TFile(Form(location+"/%s/histograms/tauTau_mH%d_%s_TauUp_%s.root", outputDir.Data(), 125, bin_.c_str() , variable_.c_str()), "READ");
  //TFile* fin_tDown = new TFile(Form(location+"/%s/histograms/tauTau_mH%d_%s_TauDown_%s.root", outputDir.Data(), 125, bin_.c_str() , variable_.c_str()), "READ");
  //TFile* fin_nominal = new TFile(Form(location+"/%s/histograms/tauTau_mH%d_%s__%s.root", outputDir.Data(), 125, bin_.c_str() , variable_.c_str()), "READ");

  float rescaleggH = RESCALETO1PB ? higgsXsection(mH_,"ggH") : 1.0;
  float rescaleqqH = RESCALETO1PB ? higgsXsection(mH_,"qqH") : 1.0;
  float rescaleVH = RESCALETO1PB ? higgsXsection(mH_,"VH") : 1.0;

  string binNameSpace = "";
  if(bin_.find("inclusive")!=string::npos)
    binNameSpace = "inclusive";

  if(bin_.find("boostedZ") !=string::npos)
    binNameSpace = "boostedZ";

  string theory = !DOSUSY ? "sm" : "mssm" ;

  TFile* fTemplOut = new TFile(Form(location+"/%s/datacards/htt_tt.inputs-%s-13TeV.root",outputDir.Data(), theory.c_str()),"UPDATE");
  
  string suffix = "";
  if(analysis_.find("TauUp")!=string::npos)
    suffix = "_CMS_scale_t_tautau_13TeVUp";
  else if(analysis_.find("TauDown")!=string::npos)
    suffix = "_CMS_scale_t_tautau_13TeVDown";

  cout << "Adding histos with suffix " << suffix << endl;
  TString dirName( Form("tt_%s",binNameSpace.c_str()) );
  
  if(! (fTemplOut->cd( dirName.Data() ) ) ){

    cout << "Editing the directory for bin and variable " << binNameSpace << ", " << variable_ << endl;
    TDirectory* dir = fTemplOut->mkdir( dirName.Data() );
    dir->cd();

    TH1F* hData = ((TH1F*)fin->Get("hData"));
    hData->SetName("data_obs");
    hData->Write("data_obs");

    if(!DOSUSY){

      TH1F* hSgn2 = (TH1F*)fin->Get(Form("hGGFH%d",mH_));
      hSgn2->Scale(1./rescaleggH);
      hSgn2->SetName(Form("ggH%d%s" ,mH_,suffix.c_str()));
      hSgn2->Write(Form("ggH%d%s" ,mH_,suffix.c_str()));

      /*
      TH1F* hSgn1 = (TH1F*)fin->Get(Form("hVBFH%d",mH_));
      hSgn1->Scale(1./rescaleqqH);
      hSgn1->SetName(Form("qqH%d%s" ,mH_,suffix.c_str()));
      hSgn1->Write(Form("qqH%d%s" ,mH_,suffix.c_str()));
      
      TH1F* hSgn3 = (TH1F*)fin->Get(Form("hVH%d",mH_));
      hSgn3->Scale(1./rescaleVH);
      hSgn3->SetName(Form("VH%d%s" ,mH_,suffix.c_str()));
      hSgn3->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
      */

    }
    else{
      
      TH1F* hSgn1 = (TH1F*)fin->Get(Form("hGGH%d",mH_));
      hSgn1->SetName(Form("ggH%d%s" ,mH_,suffix.c_str()));
      hSgn1->Write(Form("ggH%d%s" ,mH_,suffix.c_str()));

      TH1F* hSgn2 = (TH1F*)fin->Get(Form("hBBH%d",mH_));
      hSgn2->SetName(Form("bbH%d%s" ,mH_,suffix.c_str()));
      hSgn2->Write(Form("bbH%d%s" ,mH_,suffix.c_str()));

      TH1F* hSMSgn2 = (TH1F*)fin->Get(Form("hGGFH%d",125));
      hSMSgn2->SetName(Form("ggH_SM%d%s" ,125,suffix.c_str()));
      hSMSgn2->Write(Form("ggH_SM%d%s" ,125,suffix.c_str()));

      TH1F* hSMSgn1 = (TH1F*)fin->Get(Form("hVBFH%d",125));
      hSMSgn1->SetName(Form("qqH_SM%d%s" ,125,suffix.c_str()));
      hSMSgn1->Write(Form("qqH_SM%d%s" ,125,suffix.c_str()));
      
      TH1F* hSMSgn3 = (TH1F*)fin->Get(Form("hVH%d",125));
      hSMSgn3->SetName(Form("VH_SM%d%s" ,125,suffix.c_str()));
      hSMSgn3->Write(Form("VH_SM%d%s" ,125,suffix.c_str()));
    }

    if(bin_.find("boostedZ")!=string::npos){

      // ZTT: embedded sample
      TH1F* hDataEmb = ((TH1F*)fin->Get(nameZtt));
      hDataEmb->SetName(Form("ZTT%s",suffix.c_str()));
      hDataEmb->Write(Form("ZTT%s",suffix.c_str()));

      // ----- QCD ------
      TH1F *hQCD = ((TH1F*)fin->Get("hQCD"));
      hQCD->SetName(Form("QCD%s" ,suffix.c_str()));
      checkValidity(hQCD);
      hQCD->Write(Form("QCD%s" ,suffix.c_str()));
      
      // ----- W ------
      TH1F* hW = ((TH1F*)fin->Get("hW"));
      hW->SetName(Form("W%s" ,suffix.c_str()));
      hW->Write(Form("W%s" ,suffix.c_str()));

      // ----- Fakes ------
      TH1F* hZj = ((TH1F*)fin->Get("hZj"));
      hZj->SetName(Form("ZJ%s" ,suffix.c_str()));
      hZj->Write(Form("ZJ%s" ,suffix.c_str()));

      TH1F* hZl = ((TH1F*)fin->Get("hZl"));
      hZl->SetName(Form("ZL%s" ,suffix.c_str()));
      hZl->Write(Form("ZL%s" ,suffix.c_str()));

      TH1F* hZfakes = ((TH1F*)fin->Get("hZfakes"));
      hZfakes->SetName(Form("ZLL%s" ,suffix.c_str()));
      hZfakes->Write(Form("ZLL%s" ,suffix.c_str()));

      TH1F* hTTb = ((TH1F*)fin->Get("hTTb"));
      hTTb->SetName(Form("TT%s" ,suffix.c_str()));
      hTTb->Write(Form("TT%s" ,suffix.c_str()));

      TH1F* hVV = ((TH1F*)fin->Get("hVV"));
      hVV->SetName(Form("VV%s" ,suffix.c_str()));
      hVV->Write(Form("VV%s" ,suffix.c_str()));
    }
    else{
      
      // ZTT: embedded sample
      TH1F* hDataEmb = ((TH1F*)fin->Get(nameZtt));
      hDataEmb->SetName(Form("ZTT%s",suffix.c_str()));
      checkValidity(hDataEmb);
      hDataEmb->Write(Form("ZTT%s",suffix.c_str()));
      
      // ----- QCD ------
      TH1F *hQCD = ((TH1F*)fin->Get("hQCD"));
      hQCD->SetName(Form("QCD%s" ,suffix.c_str()));
      checkValidity(hQCD);
      hQCD->Write(Form("QCD%s" ,suffix.c_str()));

      // ----- W ------
      TH1F* hW = ((TH1F*)fin->Get("hW"));
      hW->SetName(Form("W%s" ,suffix.c_str()));
      checkValidity(hW);
      hW->Write(Form("W%s" ,suffix.c_str()));

      // ----- Fakes ------
      TH1F* hZj = ((TH1F*)fin->Get("hZj"));
      hZj->SetName(Form("ZJ%s" ,suffix.c_str()));
      checkValidity(hZj);
      hZj->Write(Form("ZJ%s" ,suffix.c_str()));

      TH1F* hZl = ((TH1F*)fin->Get("hZl"));
      hZl->SetName(Form("ZL%s" ,suffix.c_str()));
      checkValidity(hZl);
      hZl->Write(Form("ZL%s" ,suffix.c_str()));

      TH1F* hZfakes = ((TH1F*)fin->Get("hZfakes"));
      hZfakes->SetName(Form("ZLL%s" ,suffix.c_str()));
      checkValidity(hZfakes);
      hZfakes->Write(Form("ZLL%s" ,suffix.c_str()));

      TH1F* hTTb = ((TH1F*)fin->Get("hTTb"));
      hTTb->SetName(Form("TT%s" ,suffix.c_str()));
      checkValidity(hTTb);
      hTTb->Write(Form("TT%s" ,suffix.c_str()));

      TH1F* hVV = ((TH1F*)fin->Get("hVV"));
      hVV->SetName(Form("VV%s" ,suffix.c_str()));
      checkValidity(hVV);
      hVV->Write(Form("VV%s" ,suffix.c_str()));

    }

  }
  else{

    cout << "Directory is there, filling new histos..." << endl;
    TDirectory* dir = (TDirectory*)fTemplOut->Get( dirName.Data() );
    fTemplOut->cd( dirName.Data() );

    if(!DOSUSY){

      if(dir->FindObjectAny(Form("ggH%d%s" , mH_,suffix.c_str()))==0 ){
        TH1F* hSgn1 = (TH1F*)fin->Get(Form("hGGFH%d",mH_));
	hSgn1->Scale(1./rescaleggH);
        hSgn1->SetName(Form("ggH%d%s" , mH_,suffix.c_str()));
        hSgn1->Write(Form("ggH%d%s" , mH_,suffix.c_str()));
      }
      /*
      if(dir->FindObjectAny(Form("qqH%d%s" ,mH_,suffix.c_str()))==0 ){
	TH1F* hSgn2 = (TH1F*)fin->Get(Form("hVBFH%d",mH_));
	hSgn2->Scale(1./rescaleqqH);
	hSgn2->SetName(Form("qqH%d%s" ,mH_,suffix.c_str()));
	hSgn2->Write(Form("qqH%d%s" ,mH_,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("VH%d%s" , mH_,suffix.c_str()))==0 ){
	TH1F* hSgn3 = (TH1F*)fin->Get(Form("hVH%d",mH_));
	hSgn3->Scale(1./rescaleVH);
	hSgn3->SetName(Form("VH%d%s" ,mH_,suffix.c_str()));
	hSgn3->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
      }
      */

    }
    else{

      if(dir->FindObjectAny(Form("ggH%d%s" ,mH_,suffix.c_str()))==0 ){
	TH1F* hSgn1 = (TH1F*)fin->Get(Form("hGGH%d",mH_));
	hSgn1->SetName(Form("ggH%d%s" ,mH_,suffix.c_str()));
	hSgn1->Write(Form("ggH%d%s" ,mH_,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("bbH%d%s" , mH_,suffix.c_str()))==0 ){
	TH1F* hSgn2 = (TH1F*)fin->Get(Form("hBBH%d",mH_));
	hSgn2->SetName(Form("bbH%d%s" , mH_,suffix.c_str()));
	hSgn2->Write(Form("bbH%d%s" , mH_,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("ggH_SM%d%s" ,mH_,suffix.c_str()))==0 ){
	TH1F* hSMSgn2 = (TH1F*)fin->Get(Form("hGGFH%d",125));
	hSMSgn2->SetName(Form("ggH_SM%d%s" ,125,suffix.c_str()));
	hSMSgn2->Write(Form("ggH_SM%d%s" ,125,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("qqH_SM%d%s" ,mH_,suffix.c_str()))==0 ){
	TH1F* hSMSgn1 = (TH1F*)fin->Get(Form("hVBFH%d",125));
	hSMSgn1->SetName(Form("qqH_SM%d%s" ,125,suffix.c_str()));
	hSMSgn1->Write(Form("qqH_SM%d%s" ,125,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("VH_SM%d%s" ,mH_,suffix.c_str()))==0 ){
	TH1F* hSMSgn3 = (TH1F*)fin->Get(Form("hVH%d",125));
	hSMSgn3->SetName(Form("VH_SM%d%s" ,125,suffix.c_str()));
	hSMSgn3->Write(Form("VH_SM%d%s" ,125,suffix.c_str()));
      }
    }

    if(bin_.find("boostedZ")!=string::npos){

      // ZTT: embedded sample
      if(dir->FindObjectAny(Form("ZTT%s" ,suffix.c_str()))==0 ){
	TH1F* hDataEmb = ((TH1F*)fin->Get(nameZtt));
	hDataEmb->SetName(Form("ZTT%s",suffix.c_str()));
	hDataEmb->Write(Form("ZTT%s",suffix.c_str()));
      }

      // ----- QCD ------
      if(dir->FindObjectAny(Form("QCD%s" ,suffix.c_str()))==0 ){
	TH1F *hQCD = ((TH1F*)fin->Get("hQCD"));
	hQCD->SetName(Form("QCD%s" ,suffix.c_str()));
	checkValidity(hQCD);
	hQCD->Write(Form("QCD%s" ,suffix.c_str()));
      }
      
      // ----- W ------ 
      if(dir->FindObjectAny(Form("W%s" ,suffix.c_str()))==0 ){
      TH1F* hW = ((TH1F*)fin->Get("hW"));
      hW->SetName(Form("W%s" ,suffix.c_str()));
      hW->Write(Form("W%s" ,suffix.c_str()));
      }
      
      // ----- Fakes ------
      if(dir->FindObjectAny(Form("ZJ%s" ,suffix.c_str()))==0 ){
	TH1F* hZj = ((TH1F*)fin->Get("hZj"));
	hZj->SetName(Form("ZJ%s" ,suffix.c_str()));
	hZj->Write(Form("ZJ%s" ,suffix.c_str()));
      }
      
      if(dir->FindObjectAny(Form("ZL%s" ,suffix.c_str()))==0 ) {
	TH1F* hZl = ((TH1F*)fin->Get("hZl"));
	hZl->SetName(Form("ZL%s" ,suffix.c_str()));
	hZl->Write(Form("ZL%s" ,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("ZLL%s" ,suffix.c_str()))==0 ) {
	TH1F* hZfakes = ((TH1F*)fin->Get("hZfakes"));
	hZfakes->SetName(Form("ZLL%s" ,suffix.c_str()));
	hZfakes->Write(Form("ZLL%s" ,suffix.c_str()));
      }
      
      if(dir->FindObjectAny(Form("TT%s" ,suffix.c_str()))==0 ){
	TH1F* hTTb = ((TH1F*)fin->Get("hTTb"));
	hTTb->SetName(Form("TT%s" ,suffix.c_str()));
	hTTb->Write(Form("TT%s" ,suffix.c_str()));
      }
      
      if(dir->FindObjectAny(Form("VV%s" ,suffix.c_str()))==0 ){
	TH1F* hVV = ((TH1F*)fin->Get("hVV"));
	hVV->SetName(Form("VV%s" ,suffix.c_str()));
	hVV->Write(Form("VV%s" ,suffix.c_str()));
      }

    }
    else{
      
      // ZTT: embedded sample 
      if(dir->FindObjectAny(Form("ZTT%s" ,suffix.c_str()))==0 ){
	TH1F* hDataEmb = ((TH1F*)fin->Get(nameZtt));
        hDataEmb->SetName(Form("ZTT%s",suffix.c_str()));
        hDataEmb->Write(Form("ZTT%s",suffix.c_str()));
      }

      // ----- QCD ------  
      if(dir->FindObjectAny(Form("QCD%s" ,suffix.c_str()))==0 ){
        TH1F *hQCD = ((TH1F*)fin->Get("hQCD"));
        hQCD->SetName(Form("QCD%s" ,suffix.c_str()));
        checkValidity(hQCD);
        hQCD->Write(Form("QCD%s" ,suffix.c_str()));
      }

      // ----- W ------ 
      if(dir->FindObjectAny(Form("W%s" ,suffix.c_str()))==0 ){
	TH1F* hW = ((TH1F*)fin->Get("hW"));
	hW->SetName(Form("W%s" ,suffix.c_str()));
	hW->Write(Form("W%s" ,suffix.c_str()));
      }

      // ----- Fakes ------
      if(dir->FindObjectAny(Form("ZJ%s" ,suffix.c_str()))==0 ){
	TH1F* hZj = ((TH1F*)fin->Get("hZj"));
        hZj->SetName(Form("ZJ%s" ,suffix.c_str()));
	hZj->Write(Form("ZJ%s" ,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("ZL%s" ,suffix.c_str()))==0 ) {
        TH1F* hZl = ((TH1F*)fin->Get("hZl"));
        hZl->SetName(Form("ZL%s" ,suffix.c_str()));
        hZl->Write(Form("ZL%s" ,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("ZLL%s" ,suffix.c_str()))==0 ) {
        TH1F* hZfakes = ((TH1F*)fin->Get("hZfakes"));
        hZfakes->SetName(Form("ZLL%s" ,suffix.c_str()));
        hZfakes->Write(Form("ZLL%s" ,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("TT%s" ,suffix.c_str()))==0 ){
        TH1F* hTTb = ((TH1F*)fin->Get("hTTb"));
        hTTb->SetName(Form("TT%s" ,suffix.c_str()));
        hTTb->Write(Form("TT%s" ,suffix.c_str()));
      }

      if(dir->FindObjectAny(Form("VV%s" ,suffix.c_str()))==0 ){
        TH1F* hVV = ((TH1F*)fin->Get("hVV"));
        hVV->SetName(Form("VV%s" ,suffix.c_str()));
        hVV->Write(Form("VV%s" ,suffix.c_str()));
      }

    }
  }

  fTemplOut->Close();
  delete fTemplOut;

  // edit the datacards only for the nominal analysis
  if(analysis_.find("Up")!=string::npos || analysis_.find("Down")!=string::npos){
    fin->Close(); delete fin;
    //fin_nominal->Close(); delete fin_nominal;
    //fin_tUp->Close(); delete fin_tUp;
    //fin_tDown->Close(); delete fin_tDown;
    return;
  }
  
}


void produceAll(TString outputDir="DiTauV1", int useEmb=0){

  const int nVar=1;
  const int nM=1; // 3 for sm, 21 for mssm
  const int nCat=2; //9 sm; //6 before including "inclusive", 3 for MSSM
  const int nAn=3;

  //string variables[nVar]={"diTauNSVfitMass"};
  string variables[nVar]={"diTauVisMass"};

  int mH[nM]={125};
  //int mH[nM]={80,90,100,110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000};

  string analysis[nAn]={"","TauUp","TauDown"};
  string category[nCat]={"inclusive","boostedZ"};
  //string category[nCat]={"inclusive","novbfLow","novbfMedium","novbfHigh","boostMedium","boostHighhighhiggs","boostHighlowhiggs","vbf","vbfTight"};
  //string category[nCat]={"inclusive","bTag", "nobTag"}; //for MSSM

  for(int iVar = 0 ; iVar < nVar; iVar++)
    for(int iM = 0; iM < nM; iM++)
      for(int iAn=0 ; iAn<nAn ; iAn++)
	for(int iCat=0 ; iCat<nCat ; iCat++)
	  produce(mH[iM],variables[iVar], analysis[iAn], category[iCat], outputDir, useEmb);

}

int main(int argc, const char* argv[])
{
  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  TString outputDir="DiTauV1";
  int useEmb=0;
  if(argc>=2) outputDir = argv[1];
  if(argc>=3) useEmb = (int)atof(argv[2]);
  produceAll(outputDir, useEmb);

}
