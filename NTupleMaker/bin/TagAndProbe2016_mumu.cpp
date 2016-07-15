
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//
//    This code runs on events ntuples (data or MC), finds the T&P mumu pairs, evaluates 
//    id and trigger criterias, computes isolation and produces the root file with the T&P otree.
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
float rel_Iso(int Index, TString lep, const AC1B * analysisTree, float dRCone);//dr cone 0.3 OR 0.4 
bool isICHEPmed(int Index, const AC1B * analysisTree); 

int main(int argc, char * argv[]){

	// first argument - config file for analysis
  // second argument - file list to be analyzed
  // third argument - optional index of the first file to be analyzed
  // fourth argument - optional index of the first file to be analyzed

	using namespace std;

  gErrorIgnoreLevel = kFatal;

	string cmsswBase = (getenv ("CMSSW_BASE"));

  //config file reading

  Config cfg(argv[1]);
  const string infiles = argv[2];

  const bool isData = cfg.get<bool>("isData");

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }


  const float ptProbeMuonCut = cfg.get<float>("ptProbeMuonCut");
  const float etaProbeMuonCut = cfg.get<float>("etaProbeMuonCut");
  const float dxyProbeMuonCut = cfg.get<float>("dxyProbeMuonCut");
  const float dzProbeMuonCut = cfg.get<float>("dzProbeMuonCut");

  const float isoMuonCut = cfg.get<float>("isoMuonCut");
  const float dRPairCut = cfg.get<float>("dRPairCut");

  const float ptMuonCut = cfg.get<float>("ptMuonCut");
  const float etaMuonCut = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut = cfg.get<float>("dzMuonCut");
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

  //pileup distrib config
  const bool ApplyPUweight = cfg.get<bool>("ApplyPUweight"); 
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");

  string isoLeg;
  if (isData || ApplyTrigger){
    isoLeg = cfg.get<string>("isoLeg");
  }

  //hlt filters to be evaluated inizialization
  vector<string> hlt;

  for(unsigned int i=1; i<=nhlt_check; i++){
    hlt.push_back(cfg.get<string>("hlt_" + std::to_string(i)));
  }

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
  const string sample = argv[2];
  TString rootFileName(sample);
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_TP_Muon.root";
  
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

  for (int iF=ifile; iF<jfile; ++iF) {  //FILEs LOOP

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

  

    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++){ //EVENTs LOOP

      analysisTree.GetEntry(iEntry);

      // if checkRun flag is true only use runs between nRunMin and nRunMax
      if (checkRun){
        if ((analysisTree.event_run<nRunMin) || (analysisTree.event_run>nRunMax)){
          continue;
        }
      }

      nEvents++;

      //filters
      unsigned int nfilters = analysisTree.run_hltfilters->size();

      //find the trigger legs to be matchd by the tags
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

      //hlt filters to be evaluated indices finding
      int *nHLT = new int[nhlt_check];

      for(unsigned int i2=0; i2<nhlt_check; ++i2) {nHLT[i2] = -1;}

      for (unsigned int i=0; i<nfilters; ++i) {
      	TString HLTFilter(analysisTree.run_hltfilters->at(i));
        for(unsigned int i2=0; i2<nhlt_check; ++i2){
          if (HLTFilter==hlt[i2]) nHLT[i2] = i;
//          if((i2==11) && (nHLT[i2]!=-1)) cout<<"debug "<<hlt[i2]<<": "<<nHLT[i2]<<endl; //OK
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

      for (it = 0; it<analysisTree.muon_count; it++){

//        bool TagmuonMediumId = analysisTree.muon_isMedium[it]; 

        bool TagmuonMediumId = isICHEPmed(it, &analysisTree);

        if (analysisTree.muon_pt[it]<=ptMuonCut) continue;
        if (rel_Iso(it, "m", &analysisTree, dRCone)>=isoMuonCut) continue;
        if (fabs(analysisTree.muon_eta[it])>=etaMuonCut) continue;
        if (fabs(analysisTree.muon_dxy[it])>=dxyMuonCut) continue;
        if (fabs(analysisTree.muon_dz[it])>=dzMuonCut) continue;
        if (!TagmuonMediumId) continue;

        //trigger match
        bool isSingleLepTrig = false;

        
        if(isData || ApplyTrigger){
          if(debug) cout<<"analysisTree.trigobject_count: "<<analysisTree.trigobject_count<<endl;
          for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
            float dRtrig = deltaR(analysisTree.muon_eta[it], analysisTree.muon_phi[it], 
                                  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
  
            if (dRtrig < deltaRTrigMatch){
              if(debug) cout<<"[iT][nIsoLeg] = "<<iT<<" - "<<nIsoLeg<<" = "<<analysisTree.trigobject_filters[iT][nIsoLeg]<<endl;
              if (analysisTree.trigobject_filters[iT][nIsoLeg] && ( isData || analysisTree.trigobject_pt[iT] > ptTrigObjCut)) // Ele23 Leg
                isSingleLepTrig = true;
            }
          }
          
          if (!isSingleLepTrig) {
            if(debug) {cout<<"debug: tag trigger match failed"<<endl;}
            continue;}
          if(debug) {cout<<"debug: tag trigger match OK"<<endl;}
        }

        
        
        otree->pt_tag = analysisTree.muon_pt[it]; 
        otree->eta_tag = analysisTree.muon_eta[it];
        otree->phi_tag = analysisTree.muon_phi[it];


        //probe selection
        for (ip = 0; ip<analysisTree.muon_count; ip++){
          if (it == ip) continue;

          if (analysisTree.muon_pt[ip]<=ptProbeMuonCut) continue;
          if (fabs(analysisTree.muon_eta[ip])>=etaProbeMuonCut) continue;
          if (fabs(analysisTree.muon_dxy[ip])>=dxyProbeMuonCut) continue;
          if (fabs(analysisTree.muon_dz[ip])>=dzProbeMuonCut) continue;

          if ((analysisTree.muon_charge[it] * analysisTree.muon_charge[ip]) > 0.) continue;
          float dR = deltaR(analysisTree.muon_eta[it], analysisTree.muon_phi[it],
                analysisTree.muon_eta[ip], analysisTree.muon_phi[ip]);
          if (dR < dRPairCut) continue;

          tagLV.SetXYZM(analysisTree.muon_px[it], analysisTree.muon_py[it], analysisTree.muon_pz[it], muonMass);
          probeLV.SetXYZM(analysisTree.muon_px[ip], analysisTree.muon_py[ip], analysisTree.muon_pz[ip], muonMass);
          pairLV = tagLV + probeLV;

          float m_vis = pairLV.M();

          if (m_vis < massCut) continue;
          
          p_total++;

          tags.push_back(it);
          probes.push_back(ip);

          otree->pt_probe = analysisTree.muon_pt[ip];
          otree->eta_probe = analysisTree.muon_eta[ip];
          otree->phi_probe = analysisTree.muon_phi[ip];
          otree->m_vis = m_vis;
          otree->mcweight = 1.;
          if(!isData)
            otree->mcweight = analysisTree.genweight;
          //cout<<analysisTree.genweight<<endl;

          //id evaluating for the probes
          bool id_probe= false;

          if (isICHEPmed(ip, &analysisTree)) {
            if (analysisTree.muon_dxy[ip]<dxyPassingCut){
              if (analysisTree.muon_dz[ip]<dxyPassingCut){
                id_probe = true;
              }
            }
          }

          //computation of the isolation
          float iso_probe = rel_Iso(ip, "m", &analysisTree, dRCone);

          otree->id_probe = id_probe;
          otree->iso_probe = iso_probe;

          //TRIGGER evaluating for the probes

          otree->hlt_1_probe = -1;
          otree->hlt_2_probe = -1;
          otree->hlt_3_probe = -1;
          otree->hlt_4_probe = -1;
          otree->hlt_5_probe = -1;
          otree->hlt_6_probe = -1;
          otree->hlt_7_probe = -1;
          otree->hlt_8_probe = -1;
          otree->hlt_9_probe = -1;
          otree->hlt_10_probe = -1;
          otree->hlt_11_probe = -1;
          otree->hlt_12_probe = -1;
          otree->hlt_13_probe = -1;
          otree->hlt_14_probe = -1;
          otree->hlt_15_probe = -1;
          otree->hlt_16_probe = -1;
          otree->hlt_17_probe = -1;
          otree->hlt_18_probe = -1;
          otree->hlt_19_probe = -1;
          otree->hlt_20_probe = -1;

          
          p_pass_1++;

          int *hlt_probe = new int[nhlt_check];
          for(unsigned int i=0; i<nhlt_check; ++i) {
            hlt_probe[i] = 0;
            if(nHLT[i] == -1) hlt_probe[i] = -1;
          }

          for (unsigned int iTr=0; iTr<analysisTree.trigobject_count; ++iTr){

            float dRtrig = deltaR(analysisTree.muon_eta[ip],analysisTree.muon_phi[ip],
                                  analysisTree.trigobject_eta[iTr],analysisTree.trigobject_phi[iTr]);

            if (dRtrig < deltaRTrigMatch){
              for(unsigned int i=0; i<nhlt_check; ++i){
                if(nHLT[i] == -1){
                  hlt_probe[i] = -1;
                } else{
                  if (analysisTree.trigobject_filters[iTr][nHLT[i]]){
                    hlt_probe[i] = 1;
                  }
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
          otree->hlt_7_probe = hlt_probe[6];
          otree->hlt_8_probe = hlt_probe[7];
          otree->hlt_9_probe = hlt_probe[8];
          otree->hlt_10_probe = hlt_probe[9];
          otree->hlt_11_probe = hlt_probe[10];
          otree->hlt_12_probe = hlt_probe[11];
          otree->hlt_13_probe = hlt_probe[12];
          otree->hlt_14_probe = hlt_probe[13];
          otree->hlt_15_probe = hlt_probe[14];
          otree->hlt_16_probe = hlt_probe[15];
          otree->hlt_17_probe = hlt_probe[16];
          otree->hlt_18_probe = hlt_probe[17];
          otree->hlt_19_probe = hlt_probe[18];
          otree->hlt_20_probe = hlt_probe[19];

          otree->Fill();
        }
      }
    }
    
    //closing files
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


//compute the relative isolation for a given lepton labeled by Index in channel ch with dR 0.3 or 0.4
float rel_Iso(int Index, TString lep, const AC1B * analysisTree, float dRCone){
  if(lep=="m")  return(abs_Iso(Index, lep, analysisTree, dRCone) / analysisTree->muon_pt[Index] );
  else if(lep=="e")   return(abs_Iso(Index, lep, analysisTree, dRCone) / analysisTree->electron_pt[Index] );
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

bool isICHEPmed(int Index, const AC1B * analysisTree) {
        bool goodGlob = analysisTree->muon_isGlobal[Index] && analysisTree->muon_normChi2[Index] < 3 && analysisTree->muon_combQ_chi2LocalPosition[Index] < 12
                                   && analysisTree->muon_combQ_trkKink[Index] < 20;

        bool isICHEPmedium  = analysisTree->muon_isLoose[Index] &&
                                          analysisTree->muon_validFraction[Index] >0.49 &&
                                          analysisTree->muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
        return isICHEPmedium;
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

bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}



