
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Yiwen Wen, yiwen.wen@desy.de
//
//    This code runs on events ntuples (MC), to produce pile-up weight files.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"


#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

int main(int argc, char * argv[]){

	// first argument - file list to be analyzed
    // second argument - config file
  // third argument - optional index of the first file to be analyzed
  // fourth argument - optional index of the first file to be analyzed

	using namespace std;

  gErrorIgnoreLevel = kFatal;

	string cmsswBase = (getenv ("CMSSW_BASE"));

  //config file reading
    Config cfg(argv[1]);
    const string sampleName = cfg.get<string>("sampleName");


  const string infiles = argv[2];
  const string sample = argv[2];


  //file list reading 
  int ifile = 0;
  int jfile = -1;

  if (argc > 3)
    ifile = atoi(argv[3]);
  if (argc > 4)
    jfile = atoi(argv[4]);

  // create input files list
  std::vector<std::string> fileList;  
  if (infiles.find(".root") != std::string::npos){
    ifile = 0;
    jfile = 1;

    fileList.push_back(infiles);
  }
  else{
    ifstream input;
    std::string infile;
    
    input.open(infiles);

    while(true){
      input>>infile;
      if(!input.eof()){
	if (infile.length() > 0)
	  fileList.push_back(infile);
      }
      else
	break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }  

  //output inizialization
  TString rootFileName(sample);
  rootFileName += "_PU.root";
  
  std::string ntupleName("makeroottree/AC1B");

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");
  TString puHistName(sampleName);
  puHistName += "_pileup";
  TH1D * puH = new TH1D(puHistName,"",80,0,80);//80 or 800 bins? it is a question?

  int nTotalFiles = 0;
  int nEvents = 0;
  int nFiles = 0;


  // int Run, Event, Lumi;

  for (int iF=ifile; iF<jfile; ++iF) {  //FILEs LOOP

    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue; 
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");

    if (histoInputEvents==NULL) continue;
    

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++){ //EVENTs LOOP

      analysisTree.GetEntry(iEntry);


      nEvents++;

      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl;
        
      puH->Fill(analysisTree.numtruepileupinteractions);

    }
    
    //closing files
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  cout<<endl;
  //std::cout << "-- after end of loop --" << std::endl;
  file->cd("");
  //std::cout << "-- after file->cd() --" << std::endl;
  file->Write();
  //std::cout << "-- after file->write() --" << std::endl;
  file->Close();
  //std::cout << "-- after file->close() --" << std::endl;
  delete file;
  //std::cout << "-- after delete file --" << std::endl;

  return 0;
}






