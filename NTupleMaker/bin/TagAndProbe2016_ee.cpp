
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//
//    This code runs on events ntuples (data or MC), finds the T&P ee pairs, evaluates 
//    id and trigger criterias, computes isolation and produces the root file with the T&P otree.
//
//    This is code is not commented, one can find comments in the twin file for muons.
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

#include "DesyTauAnalyses/NTupleMaker/interface/TagProbeTree.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#define pi 	3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 		   0.10565837
#define tauMass 		   1.77682
#define pionMass 		   0.1396

typedef std::vector<std::pair<int,int> > lumi_json;

struct compare_lumi { //accepts two pairs, return 1 if left.first < right.first or left.first = right.first e left.second < right.second
  bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
    if (left.first < right.first)
      return 1;
    else if (left.first > right.first)
      return 0;
    else
      return left.second < right.second;
  }
};

int read_json(std::string filename, lumi_json& json);
bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json);
bool isGoodLumi(int run, int lumi, const lumi_json& json);
float abs_Iso(int Index, TString lep, const AC1B * analysisTree, float dRCone);
float rel_Iso(int Index, TString lep, const AC1B * analysisTree, float dRCone); //dr cone 0.3 OR 0.4 

int main(int argc, char * argv[]){

	// first argument - config file for analysis
  // second argument - config file for process
  // third argumen - first file 
  //fourth argument - last file

	using namespace std;

  gErrorIgnoreLevel = kFatal;

	string cmsswBase = (getenv ("CMSSW_BASE"));

  //put config here
  Config cfg(argv[1]);
  const string infiles = argv[2];

  const bool isData = cfg.get<bool>("isData");

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  const float ptProbeElectronCut = cfg.get<float>("ptProbeElectronCut");
  const float etaProbeElectronCut = cfg.get<float>("etaProbeElectronCut");
  const float dxyProbeElectronCut = cfg.get<float>("dxyProbeElectronCut");
  const float dzProbeElectronCut = cfg.get<float>("dzProbeElectronCut");


  const float isoElectronCut = cfg.get<float>("isoElectronCut");
  const float dRPairCut = cfg.get<float>("dRPairCut");


  const float ptElectronCut = cfg.get<float>("ptElectronCut");
  const float etaElectronCut = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut = cfg.get<float>("dzElectronCut");
  const float dRCone = cfg.get<float>("dRCone");


  const float massCut = cfg.get<float>("massCut");

  const float dxyPassingCut = cfg.get<float>("dxyPassingCut");
  const float dzPassingCut = cfg.get<float>("dzPassingCut");

  const float deltaRTrigMatch = cfg.get<float>("deltaRTrigMatch");

  const unsigned int nhlt_check =  cfg.get<unsigned int>("nhlt_check");

  const unsigned int nRunMin = cfg.get<unsigned int>("nRunMin");
  const unsigned int nRunMax = cfg.get<unsigned int>("nRunMax");
  const bool checkRun = cfg.get<bool>("checkRun");


  const bool debug = cfg.get<bool>("debug");

  const bool ApplyTrigger = cfg.get<bool>("ApplyTrigger");
  const float ptTrigObjCut = cfg.get<float>("ptTrigObjCut");

  //pileup distrib
  const bool ApplyPUweight = cfg.get<bool>("ApplyPUweight"); 
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");

  string isoLeg;
  if (isData || ApplyTrigger){
    isoLeg = cfg.get<string>("isoLeg");
  }


  vector<string> hlt;

  for(unsigned int i=1; i<=nhlt_check; i++){
    hlt.push_back(cfg.get<string>("hlt_" + std::to_string(i)));
  }


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
  const string sample = argv[2];
  TString rootFileName(sample);
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_TP_Electron.root";

  std::string ntupleName("makeroottree/AC1B");

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  int p_total=0, p_pass_1=0;
  int *p_pass_hlt = new int[nhlt_check];
  int *p_tot_hlt = new int[nhlt_check];
  for(unsigned int i2=0; i2<nhlt_check; ++i2) {
    p_pass_hlt[i2]=0;
    p_tot_hlt[i2]=0;
  }

  // PU reweighting - initialization
  PileUp * PUofficial = new PileUp();
  if(ApplyPUweight){
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(pileUpInDataFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(pileUpInMCFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  

  file->cd("");
  TTree * tree = new TTree("TagProbe","TagProbe");
  TagProbeTree *otree = new TagProbeTree(tree);

  int nTotalFiles = 0;

  int nEvents = 0;
  int nFiles = 0;


  int Run, Event, Lumi;

  for (int iF=ifile; iF<jfile; ++iF) {  //FILE LOOP

    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  

    if (_tree==NULL) continue; 
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");

    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());

    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();

  

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++){ //EVENT LOOP

      analysisTree.GetEntry(iEntry);

      if (checkRun){
        if ((analysisTree.event_run<nRunMin) || (analysisTree.event_run>nRunMax)){
          continue;
        }
      }

      nEvents++;

      //filters
      unsigned int nfilters = analysisTree.run_hltfilters->size();

      //check isoleg
      unsigned int nIsoLeg = 0;
      bool checkIsoLeg = false;
      if(isData || ApplyTrigger){  
      
            for (unsigned int i=0; i<nfilters; ++i) {
              TString HLTFilter(analysisTree.run_hltfilters->at(i));
              if (HLTFilter==isoLeg) {
                nIsoLeg = i;
                checkIsoLeg = true;
              }
            }
            if (!checkIsoLeg) {
              std::cout << "HLT filter " << isoLeg << " not found" << std::endl;
              exit(-1);
            }
      }

      //hlt filters indices
      int *nHLT = new int[nhlt_check];

      for(unsigned int i2=0; i2<nhlt_check; ++i2) {nHLT[i2] = -1;}

      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
        for(unsigned int i2=0; i2<nhlt_check; ++i2){
          if (HLTFilter==hlt[i2]) nHLT[i2] = i;
        }
      }


      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      otree->run = int(analysisTree.event_run);
      otree->lumi = int(analysisTree.event_luminosityblock);
      otree->evt = int(analysisTree.event_nr); 
      

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;

      otree->pu_weight=1;

      if(ApplyPUweight) {
        otree->pu_weight = 0;
        otree->pu_weight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
      } 

      //tag selection
      unsigned int ip=0, it=0;
      vector<int> probes; probes.clear();
      vector<int> tags; tags.clear();
      TLorentzVector tagLV, probeLV, pairLV;

      for (it = 0; it<analysisTree.electron_count; it++){

        bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[it];

        if (analysisTree.electron_pt[it]<=ptElectronCut) continue;
        if (fabs(analysisTree.electron_eta[it])>=etaElectronCut) continue;
        if (fabs(analysisTree.electron_dxy[it])>=dxyElectronCut) continue;
        if (fabs(analysisTree.electron_dz[it])>=dzElectronCut) continue;

        if (rel_Iso(it, "e", &analysisTree, dRCone)>=isoElectronCut) continue;
        if (!electronMvaId) continue;
        if (analysisTree.electron_nmissinginnerhits[it]>1) continue;
        if (!analysisTree.electron_pass_conversion[it]) continue;

        //trigger match
        bool isSingleLepTrig = false;

        if(isData || ApplyTrigger){  
          for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
            float dRtrig = deltaR(analysisTree.electron_eta[it], analysisTree.electron_phi[it], 
                                  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
  
            if (dRtrig < deltaRTrigMatch){
              if (analysisTree.trigobject_filters[iT][nIsoLeg] && ( isData || analysisTree.trigobject_pt[iT] > ptTrigObjCut)) // Ele23 Leg
                isSingleLepTrig = true;
            }
          }
          
          if (!isSingleLepTrig) continue;
        }
        
        otree->pt_tag = analysisTree.electron_pt[it]; 
        otree->eta_tag = analysisTree.electron_eta[it];
        otree->phi_tag = analysisTree.electron_phi[it];


        //probe selection
        for (ip = 0; ip<analysisTree.electron_count; ip++){
          if (it == ip) continue;

          if (analysisTree.electron_pt[ip]<=ptProbeElectronCut) continue;
          if (fabs(analysisTree.electron_eta[ip])>=etaProbeElectronCut) continue;
          if (fabs(analysisTree.electron_dxy[ip])>=dxyProbeElectronCut) continue;
          if (fabs(analysisTree.electron_dz[ip])>=dzProbeElectronCut) continue;

          if ((analysisTree.electron_charge[it] * analysisTree.electron_charge[ip]) > 0.) continue;
          float dR = deltaR(analysisTree.electron_eta[it], analysisTree.electron_phi[it],
                analysisTree.electron_eta[ip], analysisTree.electron_phi[ip]);
          if (dR < dRPairCut) continue;

          tagLV.SetXYZM(analysisTree.electron_px[it], analysisTree.electron_py[it], analysisTree.electron_pz[it], electronMass);
          probeLV.SetXYZM(analysisTree.electron_px[ip], analysisTree.electron_py[ip], analysisTree.electron_pz[ip], electronMass);
          pairLV = tagLV + probeLV;

          float m_vis = pairLV.M();

          if (m_vis < massCut) continue;
          

          p_total++;

          tags.push_back(it);
          probes.push_back(ip);

          otree->pt_probe = analysisTree.electron_pt[ip];
          otree->eta_probe = analysisTree.electron_eta[ip];
          otree->phi_probe = analysisTree.electron_phi[ip];
          otree->m_vis = m_vis;
          otree->mcweight = 1.;
          if(!isData)
            otree->mcweight = analysisTree.genweight;


          bool id_probe= false;

          if ((analysisTree.electron_nmissinginnerhits[ip]<=1) && analysisTree.electron_pass_conversion[ip] && analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ip]) {
            if (analysisTree.electron_dxy[ip]<dxyPassingCut){
              if (analysisTree.electron_dz[ip]<dxyPassingCut){
                id_probe = true;
              }
            }
          }

          float iso_probe = rel_Iso(ip, "e", &analysisTree, dRCone);

          otree->id_probe = id_probe;
          otree->iso_probe = iso_probe;

          otree->hlt_1_probe = -1;
          otree->hlt_2_probe = -1;
          otree->hlt_3_probe = -1;
          otree->hlt_4_probe = -1;
          otree->hlt_5_probe = -1;
          otree->hlt_6_probe = -1;

          //TRIGGER
          p_pass_1++;

          int *hlt_probe = new int[nhlt_check];
          for(unsigned int i=0; i<nhlt_check; ++i) {hlt_probe[i] = 0;}

          for (unsigned int iTr=0; iTr<analysisTree.trigobject_count; ++iTr){

            float dRtrig = deltaR(analysisTree.electron_eta[ip],analysisTree.electron_phi[ip],
                                  analysisTree.trigobject_eta[iTr],analysisTree.trigobject_phi[iTr]);

            if (dRtrig < deltaRTrigMatch){
              for(unsigned int i=0; i<nhlt_check; ++i){
                if(nHLT[i] == -1){
                  hlt_probe[i] = -1;
                } else if (analysisTree.trigobject_filters[iTr][nHLT[i]]){
                  hlt_probe[i] = 1;
                }
              }
            }
          }
          for(unsigned int i=0; i<nhlt_check; ++i) {if(hlt_probe[i] == 1) p_pass_hlt[i]++;}
          for(unsigned int i=0; i<nhlt_check; ++i) {
            if((hlt_probe[i] != 1)&&(hlt_probe[i] != 0)) {
              //cout<<"!!!! hlt_probe["<<i<<"]="<<hlt_probe[i]<<endl;
            }
          }
          otree->hlt_1_probe = hlt_probe[0];
          otree->hlt_2_probe = hlt_probe[1];
          otree->hlt_3_probe = hlt_probe[2];
          otree->hlt_4_probe = hlt_probe[3];
          otree->hlt_5_probe = hlt_probe[4];
          otree->hlt_6_probe = hlt_probe[5];

          otree->Fill();
        }
      }
    }
  
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
//  cout<<endl<<"Total P: "<<p_total<<endl<<"Passing Id&Iso P: "<<p_pass_1<<" eff: "<<(float)p_pass_1/(float)p_total<<endl;
//  for(unsigned int i=0; i<nhlt_check; ++i) {cout<<"Passing P for hlt: "<<hlt[i]<<" : "<<p_pass_hlt[i]<<" eff: "<<(float)p_pass_hlt[i]/(float)p_pass_1<<endl;}
  cout<<endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
}

///////////////////////////////////////////////
//////////////FUNCTION DEFINITION//////////////
///////////////////////////////////////////////

float rel_Iso(int Index, TString lep, const AC1B * analysisTree, float dRCone){
  if(lep=="m")  return(abs_Iso(Index, lep, analysisTree, dRCone) / analysisTree->muon_pt[Index] );
  else if(lep="e")   return(abs_Iso(Index, lep, analysisTree, dRCone) / analysisTree->electron_pt[Index] );
    else return(-1.);
}

float abs_Iso (int Index, TString lep, const AC1B * analysisTree, float dRCone){
  float neutralHadIso, photonIso, chargedHadIso, puIso;

  if(lep=="m"){
    neutralHadIso = analysisTree->muon_neutralHadIso[Index];
    photonIso =     analysisTree->muon_photonIso[Index];
    chargedHadIso = analysisTree->muon_chargedHadIso[Index];
    puIso =         analysisTree->muon_puIso[Index];
    if (dRCone == 0.3) {
      neutralHadIso =     analysisTree->muon_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r03_sumPUPt[Index];
    }
    if (dRCone == 0.4) {
      neutralHadIso =     analysisTree->muon_r04_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r04_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r04_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r04_sumPUPt[Index];
    }
  }
  if(lep=="e"){
    neutralHadIso = analysisTree->electron_neutralHadIso[Index];
    photonIso =     analysisTree->electron_photonIso[Index];
    chargedHadIso = analysisTree->electron_chargedHadIso[Index];
    puIso =         analysisTree->electron_puIso[Index];
    if (dRCone == 0.3) {
      neutralHadIso =     analysisTree->electron_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->electron_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->electron_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->electron_r03_sumPUPt[Index];
    }
  }

  float neutralIso = neutralHadIso + photonIso -0.5*puIso;
  neutralIso = TMath::Max(float(0), neutralIso);
  return(chargedHadIso + neutralIso);
}

int read_json(std::string filename, lumi_json& json){

  std::pair <int,int> lumi;

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);

  BOOST_FOREACH(boost::property_tree::ptree::value_type &json_run, pt.get_child("")){
    int irun = atoi(json_run.first.data());
    BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_ranges, json_run.second.get_child("")){
      int ilumi[2] = {};

      int count = 0;
      BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_boundaries, lumi_ranges.second.get_child("")){
  ilumi[count] = atoi(lumi_boundaries.second.data().data());
  count++;
      }
      
      for (;ilumi[0] <= ilumi[1]; ilumi[0]++){
  lumi = std::make_pair(irun, ilumi[0]);
  json.push_back(lumi);
      }
    }
  }

  sort( json.begin(), json.end(),  compare_lumi());
  json.erase( unique( json.begin(), json.end() ), json.end() );
  
  return 0;
}


bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json){
  static compare_lumi compare;
  static std::pair<int,int> oldlumi = lumi;
  static bool old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  
  if(lumi.first != oldlumi.first || lumi.second != oldlumi.second){
    oldlumi = lumi;
    old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  }

  return old;
}

//accepts run number, lumi and json, make pair of (run,lumi) and starts isGoodLumi 
bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}