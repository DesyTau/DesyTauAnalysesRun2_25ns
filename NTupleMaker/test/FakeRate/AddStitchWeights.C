#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <math.h>

using namespace std;

void AddStitchWeights(TString filename = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",const size_t nSamples=4){
    Float_t refEvents[nSamples+1]; 
    Float_t refXSec[nSamples+1];
    if(filename.Contains("DY")){
	refXSec[0] = 6225.42;
	refXSec[1] = 1.165*877.8;
	refXSec[2] = 1.165*304.4;
	refXSec[3] = 1.165*111.5;
	refXSec[4] = 1.165*44.03;
    }
    else if(filename.Contains("WJets")){
	refXSec[0] = 61526.7;
	refXSec[1] = 1.1622*8104.0;
	refXSec[2] = 1.1622*2793.0;
	refXSec[3] = 1.1622*992.5;
	refXSec[4] = 1.1622*544.3;
    }else{ 
      cout << "ERROR Stiching used on unknown samples" <<endl;
      return;
    }
    TString refSamples[nSamples+1];
    refSamples[0]=filename;
    for(size_t i=1;i<(nSamples+1);i++){
      TString name=filename;
      TString njets=TString::Itoa(i,10);
      njets+="Jets";
      refSamples[i]=name.ReplaceAll("Jets",njets);
    }
    for(size_t i=0; i<(nSamples+1); i++) cout << refSamples[i] <<endl;
    for (size_t i=0;i<(nSamples+1);i++) {
      TFile * file = new TFile("./"+refSamples[i]+".root");
      TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
      refEvents[i] = histWeightsH->GetSumOfWeights();
    }
    Double_t StitchWeight[nSamples+1];
    StitchWeight[0]=refXSec[0]/refEvents[0];
    for(size_t i=1;i<(nSamples+1);i++) {
      StitchWeight[i]=1./(refEvents[0]/refXSec[0]+refEvents[i]/refXSec[i]);
    }
    for(size_t i=0;i<(nSamples+1);i++) cout << StitchWeight[i] << endl;
    //return;

    for(size_t iSample=0;iSample<(nSamples+1);iSample++) {
      TFile *f = new TFile("./"+refSamples[iSample]+".root","update");
      TTree *T = (TTree*)f->Get("ETauFR");
      Double_t stitchweight = 1;
      UInt_t npartons = 0; 
      T->SetBranchAddress("npartons",&npartons);
      TBranch *b_stitchweight=T->Branch("stitchweight",&stitchweight,"stitchweight/D");
      size_t nentries = T->GetEntries();
      cout << refSamples[iSample] << ": " << nentries << "events" <<endl;

      for (size_t i=0;i<nentries;i++) {
	T->GetEntry(i);
	if(iSample==0){
	  if(npartons==0||npartons>4)stitchweight=StitchWeight[0];
	  else stitchweight=StitchWeight[npartons];
	}else stitchweight=StitchWeight[iSample];
	b_stitchweight->Fill();
	if(i%100000==0)cout << i << " events processed" << endl;
      }
      T->Print();
      T->Write("",TObject::kOverwrite);
      delete f;
       
    }
}
