
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

#include "DesyTauAnalyses/NTupleMaker/interface/diMuon2016Tree.h"

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
float abs_Iso(int Index, TString lep, const AC1B * analysisTree, bool isIsoR03);
float rel_Iso(int Index, TString lep, const AC1B * analysisTree, bool isIsoR03);

int main(int argc, char * argv[]){

	// first argument - config file for analysis
  // second argument - config file for process

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


  const float ptMuonLowCut = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut = cfg.get<float>("dzMuonCut");
  const float isoMuonCut = cfg.get<float>("isoMuonCut");

  const float dRPairCut = cfg.get<float>("dRPairCut");

  const unsigned int nRunMin = cfg.get<unsigned int>("nRunMin");
  const unsigned int nRunMax = cfg.get<unsigned int>("nRunMax");

  const bool checkRun = cfg.get<bool>("checkRun");

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
  std::string ntupleName("makeroottree/AC1B");

  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += ".root";

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);

  TH1D * pt_0 = new TH1D("pt_0","",200,0,200);
  TH1D * pt_1 = new TH1D("pt_1","",200,0,200);
  TH1D * eta_0 = new TH1D("eta_0","",50,-2.5,2.5);
  TH1D * eta_1 = new TH1D("eta_1","",50,-2.5,2.5);
  TH1D * m_visDY = new TH1D("m_visDY","",60,60,120);
  TH1D * m_vis = new TH1D("m_vis","",200,0,200);
  TH1D * m_vislow = new TH1D("m_vislow","",500,0,20);
  TH1D * npv = new TH1D("npv","",50,0,50);

  TTree * tree = new TTree("diMuon","diMuon");
  diMuon2016Tree *otree = new diMuon2016Tree(tree);

  int nTotalFiles = 0;

  int nEvents = 0;
  int nFiles = 0;

//  vector<int> runList; runList.clear();
//  vector<int> eventList; eventList.clear();

  int Run, Event, Lumi;

  for (int iF=ifile; iF<jfile; ++iF) {  //FILE LOOP

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


      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      otree->run = int(analysisTree.event_run);
      otree->lumi = int(analysisTree.event_luminosityblock);
      otree->evt = int(analysisTree.event_nr); 

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;

      //first muon selection
      vector<int> muon_0; muon_0.clear();
      vector<int> muon_1; muon_1.clear();
      vector<int> used_muons; used_muons.clear();
      TLorentzVector LV_0, LV_2, LV_pair;

      for (unsigned int im = 0; im<analysisTree.muon_count; im++){

        if(std::find(used_muons.begin(), used_muons.end(), im) != used_muons.end()) continue;

        bool TagmuonMediumId = analysisTree.muon_isMedium[im]; 
        if (analysisTree.muon_pt[im]<=ptMuonHighCut) continue;
        if (rel_Iso(im, "m", &analysisTree, true)>=isoMuonCut) continue;
        if (fabs(analysisTree.muon_eta[im])>=etaMuonCut) continue;
        if (fabs(analysisTree.muon_dxy[im])>=dxyMuonCut) continue;
        if (fabs(analysisTree.muon_dz[im])>=dzMuonCut) continue;
        if (!TagmuonMediumId) continue;

        otree->pt_0 = analysisTree.muon_pt[im];
        otree->eta_0 = analysisTree.muon_eta[im];
        otree->phi_0 = analysisTree.muon_phi[im];
        
        //second muon selection
        for (unsigned im2 = 0; im2<analysisTree.muon_count; im2++){
          if(std::find(used_muons.begin(), used_muons.end(), im2) != used_muons.end()) continue;
          if (im == im2) continue;
          bool TagmuonMediumId2 = analysisTree.muon_isMedium[im2]; 
          if (analysisTree.muon_pt[im2]<=ptMuonLowCut) continue;
          if (rel_Iso(im2, "m", &analysisTree, true)>=isoMuonCut) continue;
          if (fabs(analysisTree.muon_eta[im2])>=etaMuonCut) continue;
          if (fabs(analysisTree.muon_dxy[im2])>=dxyMuonCut) continue;
          if (fabs(analysisTree.muon_dz[im2])>=dzMuonCut) continue;
          if (!TagmuonMediumId2) continue;;

          if ((analysisTree.muon_charge[im] * analysisTree.muon_charge[im2]) > 0.) continue;
          float dR = deltaR(analysisTree.muon_eta[im], analysisTree.muon_phi[im],
                analysisTree.muon_eta[im2], analysisTree.muon_phi[im2]);
          if (dR < dRPairCut) continue;

          otree->pt_1 = analysisTree.muon_pt[im2];
          otree->eta_1 = analysisTree.muon_eta[im2];
          otree->phi_1 = analysisTree.muon_phi[im2];

          LV_0.SetXYZM(analysisTree.muon_px[im], analysisTree.muon_py[im], analysisTree.muon_pz[im], muonMass);
          LV_2.SetXYZM(analysisTree.muon_px[im2], analysisTree.muon_py[im2], analysisTree.muon_pz[im2], muonMass);
          LV_pair = LV_0 + LV_2;

          float mass = LV_pair.M();

          muon_0.push_back(im);
          muon_1.push_back(im2);

          used_muons.push_back(im);
          used_muons.push_back(im2);

          otree->m_vis = mass;
          otree->npv = analysisTree.primvertex_count;

          otree->Fill();

          pt_0->Fill(analysisTree.muon_pt[im]);
          pt_1->Fill(analysisTree.muon_pt[im2]);
          eta_0->Fill(analysisTree.muon_eta[im]);
          eta_1->Fill(analysisTree.muon_eta[im2]);
          m_visDY->Fill(mass);
          m_vis->Fill(mass);
          m_vislow->Fill(mass);
          npv->Fill(analysisTree.primvertex_count);
          break;
        }
      }
    }
  
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }

  file->cd("");
  file->Write();
  file->Close();
  delete file;
}

///////////////////////////////////////////////
//////////////FUNCTION DEFINITION//////////////
///////////////////////////////////////////////


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

  float abs_Iso (int Index, TString lep, const AC1B * analysisTree, bool isIsoR03){
  float neutralHadIso, photonIso, chargedHadIso, puIso;

  if(lep=="m"){
    neutralHadIso = analysisTree->muon_neutralHadIso[Index];
    photonIso =     analysisTree->muon_photonIso[Index];
    chargedHadIso = analysisTree->muon_chargedHadIso[Index];
    puIso =         analysisTree->muon_puIso[Index];
    if (isIsoR03) {
      neutralHadIso =     analysisTree->muon_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r03_sumPUPt[Index];
    }
  }
  if(lep=="e"){
    neutralHadIso = analysisTree->electron_neutralHadIso[Index];
    photonIso =     analysisTree->electron_photonIso[Index];
    chargedHadIso = analysisTree->electron_chargedHadIso[Index];
    puIso =         analysisTree->electron_puIso[Index];
    if (isIsoR03) {
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

//compute the relative isolation for a given lepton labeled by Index in channel ch
float rel_Iso(int Index, TString lep, const AC1B * analysisTree, bool isIsoR03){
  if(lep=="m")  return(abs_Iso(Index, lep, analysisTree, isIsoR03) / analysisTree->muon_pt[Index] );
  else if(lep="e")   return(abs_Iso(Index, lep, analysisTree, isIsoR03) / analysisTree->electron_pt[Index] );
    else return(-1.);
}