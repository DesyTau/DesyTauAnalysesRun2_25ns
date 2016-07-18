
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

#include "DesyTauAnalyses/NTupleMaker/interface/Spring15Tree.h"

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
void fill_weight(const AC1B * analysisTree, Spring15Tree *otree, PileUp *PUofficial, bool isData);
float abs_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso);
float rel_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso);
void FillMuTau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso);
void FillETau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso);
void FillTau(const AC1B * analysisTree, Spring15Tree *otree, int tauIndex);
bool dilepton_veto_mt(const Config *cfg, const AC1B *analysisTree);
bool dilepton_veto_et(const Config *cfg, const AC1B *analysisTree);
bool extra_electron_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree);
bool extra_muon_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree);
void fillMET(TString ch, int leptonIndex, int tauIndex, const AC1B * analysisTree, Spring15Tree *otree);
void mt_calculation(Spring15Tree *otree);
void counting_jets(const AC1B *analysisTree, Spring15Tree *otree, const Config *cfg);
bool isICHEPmed(int Index, const AC1B * analysisTree);



int main(int argc, char * argv[]){

  // first argument - config file for analysis
  // second argument - file list (MUST BE IN THE SAME DIRECTORY OF THE EXECUTABLE)
  // third argument - channel ("et" or "mt")
  // third argument - index of first file to run on (optional, ignored if only one file is used)
  // fourth argument - index of last file to run on (optional, ignored if only one file is used)

  using namespace std;

  gErrorIgnoreLevel = kFatal;

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // **** configuration analysis

  Config cfg(argv[1]);

  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];

  TString ch = argv[3];
  std::string lep;

  if ((ch != "mt" ) && (ch != "et")) exit(0);

  if (ch=="mt") lep = "Muon"; 
    else if (ch=="et") lep = "Electron";

  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  const bool ApplyPUweight = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyTrigger = cfg.get<bool>("ApplyTrigger"); 

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");

  //lep eff
  const string idIsoEffFile = cfg.get<string>("idIsoEffFile");
  const string trigEffFile = cfg.get<string>("trigEffFile");

  //svfit
  const string svFitPtResFile = cfg.get<string>("svFitPtResFile");

  // MET Recoil Corrections
  const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  const bool isDY = infiles.find("DY") == infiles.rfind("/")+1;
  const bool isWJets = infiles.find("WJets") == infiles.rfind("/")+1;
  const bool isMG = infiles.find("madgraph") != string::npos;
  //const bool applyRecoilCorrections = isDY || isWJets;

  RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
  RecoilCorrector* recoilPuppiMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPuppiMetCorrector));
  RecoilCorrector* recoilMvaMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilMvaMetCorrector));

  if(!isData && applyRecoilCorrections){
    TString RecoilDir("HTT-utilities/RecoilCorrections/data/");

    TString RecoilFileName = RecoilDir; RecoilFileName += "recoilPFMEt_76X";
    if (isMG)
      RecoilFileName += "_MG5";
    RecoilFileName += ".root";
    std::cout<<RecoilFileName<<std::endl;
    recoilPFMetCorrector = new RecoilCorrector( RecoilFileName);
        
    RecoilFileName = RecoilDir; RecoilFileName += "recoilPuppiMet.root";
    std::cout<<RecoilFileName<<std::endl;
    recoilPuppiMetCorrector = new RecoilCorrector( RecoilFileName);

    RecoilFileName = RecoilDir; RecoilFileName += "recoilMvaMEt_76X_newTraining";
    if (isMG)
      RecoilFileName += "_MG5";
    RecoilFileName += ".root";
    std::cout<<RecoilFileName<<std::endl;
    recoilMvaMetCorrector = new RecoilCorrector( RecoilFileName);
  }
  
  // HLT filters
  string isoLeg;
  if (isData){
    isoLeg = cfg.get<string>("isoLegData");
  }
  else {if(ApplyTrigger) isoLeg = cfg.get<string>("isoLegMC");}

  const float ptTrigObjCut  = cfg.get<float>("ptTrigObjCut");

  // tau cuts
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut        = cfg.get<float>("dzTauCut");
  const bool applyTauId      = cfg.get<bool>("ApplyTauId");

  // pair selection
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // extra electron veto
  const float ptVetoElectronCut  = cfg.get<float>("ptVetoElectronCut");  
  const float etaVetoElectronCut = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut = cfg.get<float>("dxyVetoElectronCut");  
  const float dzVetoElectronCut  = cfg.get<float>("dzVetoElectronCut"); 
  const bool applyVetoElectronId = cfg.get<bool>("applyVetoElectronId");
  const float isoVetoElectronCut = cfg.get<float>("isoVetoElectronCut");  
  
  // extra muon veto
  const float ptVetoMuonCut  = cfg.get<float>("ptVetoMuonCut");  
  const float etaVetoMuonCut = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut = cfg.get<float>("dxyVetoMuonCut");  
  const float dzVetoMuonCut  = cfg.get<float>("dzVetoMuonCut"); 
  const bool applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
  const float isoVetoMuonCut = cfg.get<float>("isoVetoMuonCut");

 // lepton cuts
  const float ptLeptonLowCut   = cfg.get<float>("pt"+lep+"LowCut");
  const float ptLeptonHighCut  = cfg.get<float>("pt"+lep+"HighCut");
  const float etaLeptonCut     = cfg.get<float>("eta"+lep+"Cut");
  const float dxyLeptonCut     = cfg.get<float>("dxy"+lep+"Cut");
  const float dzLeptonCut      = cfg.get<float>("dz"+lep+"Cut");

  const bool  applyLeptonId    = cfg.get<bool>("Apply"+lep+"Id");

  //dilepton veto
  const float ptDiLeptonVeto     = cfg.get<float>("ptDi"+lep+"Veto");  
  const float etaDiLeptonVeto    = cfg.get<float>("etaDi"+lep+"Veto");
  const float dxyDiLeptonVeto    = cfg.get<float>("dxyDi"+lep+"Veto");  
  const float dzDiLeptonVeto     = cfg.get<float>("dzDi"+lep+"Veto"); 
  const bool applyDiLeptonVetoId = cfg.get<bool>("applyDi"+lep+"VetoId");
  const bool applyDiLeptonOS     = cfg.get<bool>("applyDi"+lep+"OS");
  const float isoDiLeptonVeto    = cfg.get<float>("isoDi"+lep+"Veto");
  const float drDiLeptonVeto     = cfg.get<float>("drDi"+lep+"Veto"); 

  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
  const float dRiso = cfg.get<float>("dRiso");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");

  // check overlap
  const bool checkOverlap = cfg.get<bool>("CheckOverlap");
  const bool debug = cfg.get<bool>("debug");
  
  // **** end of configuration analysis

  //file list creation

  int ifile = 0;
  int jfile = -1;

  if (argc > 4)
    ifile = atoi(argv[4]);
  if (argc > 5)
    jfile = atoi(argv[5]);
  
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

  
  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");

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

  // Lepton Scale Factors
  // Lepton Id+Iso scale factor
  ScaleFactor * SF_lepIdIso = new ScaleFactor();
  ScaleFactor * SF_lepTrigger = new ScaleFactor();

  if(ApplyLepSF){
    SF_lepIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(idIsoEffFile));
  
    // Electron SingleElectron trigger scale factor
    SF_lepTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(trigEffFile));
  }

  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_" + ch + "_Sync.root";

  std::cout <<rootFileName <<std::endl;  

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5,0.5);
  
  TTree * tree = new TTree("TauCheck","TauCheck");

  Spring15Tree *otree = new Spring15Tree(tree);

  int nTotalFiles = 0;

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  vector<int> runList; runList.clear();
  vector<int> eventList; eventList.clear();

  int nonOverlap = 0;

  std::ifstream fileEvents("overlap.txt");
  int Run, Event, Lumi;
  if (checkOverlap) {
    std::cout << "Non-overlapping events ->" << std::endl;
    while (fileEvents >> Run >> Event >> Lumi) {
      runList.push_back(Run);
      eventList.push_back(Event);
      std::cout << Run << ":" << Event << std::endl;
    }
    std::cout << std::endl;
  }
  std::ofstream fileOutput("overlap.out");


  //svFit
  TH1::AddDirectory(false);  
  //TFile* inputFile_visPtResolution = new TFile(svFitPtResFile.data());


  ///////////////FILE LOOP///////////////

  for (int iF=ifile; iF<jfile; ++iF) {

    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {       
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (isData)
	       nWeightedEventsH->Fill(0., 1.);
      else
	       nWeightedEventsH->Fill(0., analysisTree.genweight);

      unsigned int nIsoLeg = 0;
      bool checkIsoLeg = false;
      if(isData || ApplyTrigger){
            
      
            unsigned int nfilters = analysisTree.run_hltfilters->size();
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


      if (nEvents%10000==0) 
      	cout << "      processed " << nEvents << " events" << endl; 

      otree->run = int(analysisTree.event_run);
      otree->lumi = int(analysisTree.event_luminosityblock);
      otree->evt = int(analysisTree.event_nr);
      
      bool overlapEvent = true;
      for (unsigned int iEvent=0; iEvent<runList.size(); ++iEvent) {
      	if (runList.at(iEvent)==otree->run && eventList.at(iEvent)==otree->evt) {
      	  overlapEvent = false;	  
      	}
      }

      if (overlapEvent&&checkOverlap) continue;
      nonOverlap++;

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
      	continue;
      
      // weights
      if(ApplyPUweight) fill_weight(&analysisTree, otree, PUofficial, isData);
      
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;

      // tau selection

      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
        
        if (analysisTree.tau_pt[it]<=ptTauLowCut) continue;
        if (fabs(analysisTree.tau_eta[it])>=etaTauCut) continue;
        if (fabs(fabs(analysisTree.tau_charge[it])-1)>0.001) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>=dzTauCut) continue;
        if (applyTauId &&
            //analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) 
            analysisTree.tau_decayModeFinding[it] < 0.5) continue;
        
        taus.push_back(it);
      }


      //selecting leptons 
      vector<int> leptons; leptons.clear();
      if(ch == "et"){
        
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {

          bool electronMvaId = analysisTree.electron_mva_wp80_nontrig_Spring15_v1[ie];
    
          if (analysisTree.electron_pt[ie]<=ptLeptonLowCut) continue;
          if (fabs(analysisTree.electron_eta[ie])>=etaLeptonCut) continue;
          if (fabs(analysisTree.electron_dxy[ie])>=dxyLeptonCut) continue;
          if (fabs(analysisTree.electron_dz[ie])>=dzLeptonCut) continue;
          if (!electronMvaId  &&  applyLeptonId) continue;
          if (!analysisTree.electron_pass_conversion[ie]  && applyLeptonId) continue;
          if (analysisTree.electron_nmissinginnerhits[ie]>1 && applyLeptonId) continue;

          leptons.push_back(ie);
        }
      }

      if(ch == "mt"){
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

//          bool muonMediumId = analysisTree.muon_isMedium[im];
          bool muonMediumId = isICHEPmed(im, &analysisTree);
  
          if (analysisTree.muon_pt[im]<=ptLeptonLowCut) continue;
          if (fabs(analysisTree.muon_eta[im])>=etaLeptonCut) continue;
          if (fabs(analysisTree.muon_dxy[im])>=dxyLeptonCut) continue;
          if (fabs(analysisTree.muon_dz[im])>=dzLeptonCut) continue;
          if (!muonMediumId &&  applyLeptonId) continue;

          leptons.push_back(im);
        }
      }

      if (leptons.size()==0) continue;
      if (taus.size()==0) continue;


      // selecting electron and tau pair (OS or SS);
      int leptonIndex = -1;
      int tauIndex = -1;
      
      float isoLepMin = 1e+10;
      float isoTauMin = 1e+10;      


      for (unsigned int il=0; il<leptons.size(); ++il) {
        unsigned int lIndex  = leptons.at(il);

        float relIsoLep = rel_Iso(lIndex, ch, &analysisTree, dRiso);

        bool isSingleLepTrig = false;

        float lep_pt, lep_pt_max, lep_eta, lep_phi;
        if(ch=="mt"){
          lep_pt =      analysisTree.muon_pt[lIndex]; 
          lep_eta =     analysisTree.muon_eta[lIndex]; 
          lep_phi =     analysisTree.muon_phi[lIndex];}
        if(ch=="et"){         
          lep_pt =      analysisTree.electron_pt[lIndex];
          lep_eta =     analysisTree.electron_eta[lIndex]; 
          lep_phi =     analysisTree.electron_phi[lIndex];}

        if(isData || ApplyTrigger){  
                for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
                  float dRtrig = deltaR(lep_eta, lep_phi, analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
        
                  if (dRtrig < deltaRTrigMatch){
                    if (analysisTree.trigobject_filters[iT][nIsoLeg] && ( isData || analysisTree.trigobject_pt[iT] > ptTrigObjCut)) // Ele23 Leg
                      isSingleLepTrig = true;
                  }
                }
                
                if (!isSingleLepTrig) continue;
        }
        
        for (unsigned int it=0; it<taus.size(); ++it) {
          unsigned int tIndex = taus.at(it);

          float absIsoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
          float relIsoTau = absIsoTau / analysisTree.tau_pt[tIndex];

          if(ch=="mt"){lep_pt_max =  analysisTree.muon_pt[leptonIndex];}
          if(ch=="et"){lep_pt_max =  analysisTree.electron_pt[leptonIndex];}

          float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
                lep_eta, lep_phi);

          if (dR<dRleptonsCut) continue;

          // kinematic match
          if (lep_pt<=ptLeptonHighCut) continue;
          
          // change pair
//        cout<<iEntry<<". "<<lIndex<<", "<<tauIndex<<". "<<lep_pt<<" - "<<lep_pt_max;
          bool changePair =  false;
          if (relIsoLep<isoLepMin) {
//          cout<<"ele iso - ";

            changePair = true;
          }
          else if (fabs(relIsoLep - isoLepMin) < 1.e-5) {
            if (lep_pt > lep_pt_max) {
//            cout<<"ele pt - "; 

              changePair = true;
            }     
            else if (fabs(lep_pt - lep_pt_max) < 1.e-5) {
              if (absIsoTau > isoTauMin) {    // -->>> changed here, before was "<" 
//                cout<<"tau iso - "; 

                changePair = true;
              }
              else if ((absIsoTau - isoTauMin) < 1.e-5){
                if (analysisTree.tau_pt[tIndex] > analysisTree.tau_pt[tauIndex]) {
//                cout<<"tau pt"; 

                  changePair = true;
                }
              }
            }
          }
//        cout<<endl;

          if (changePair){

            isoLepMin  = relIsoLep;
            leptonIndex = lIndex;
            isoTauMin = absIsoTau;
            tauIndex = tIndex;
          }
        }
      }
      
      if (leptonIndex<0) continue;
      if (tauIndex<0) continue;
//      cout << "Selected Pair (e,tau) = "<<leptonIndex<<", "<<tauIndex<<std::endl;

      //filling variables
      TLorentzVector leptonLV;

      if(ch=="mt") {
      	FillMuTau(&analysisTree, otree, leptonIndex, dRiso);
      	
        leptonLV.SetXYZM(analysisTree.muon_px[leptonIndex],
					    analysisTree.muon_py[leptonIndex],
					    analysisTree.muon_pz[leptonIndex],
					    muonMass);

        if (!isData && ApplyLepSF) {
              // Scale Factor SingleEle trigger SF_eleTrigger
          otree->trigweight_1 = (SF_lepTrigger->get_ScaleFactor(double(analysisTree.muon_pt[leptonIndex]),double(analysisTree.muon_eta[leptonIndex])));
              // Scale Factor Id+Iso SF_eleIdIso
          otree->idisoweight_1 = (SF_lepIdIso->get_EfficiencyData(double(analysisTree.muon_pt[leptonIndex]),double(analysisTree.muon_eta[leptonIndex])));
          otree->effweight = (otree->trigweight_1)*(otree->idisoweight_1)*(otree->trigweight_2)*(otree->idisoweight_2);
        }

      } else if(ch=="et"){
	      FillETau(&analysisTree, otree, leptonIndex, dRiso);
	      
        leptonLV.SetXYZM(analysisTree.electron_px[leptonIndex],
						    analysisTree.electron_py[leptonIndex],
						    analysisTree.electron_pz[leptonIndex],
						    electronMass);

        if (!isData && ApplyLepSF) {
                // Scale Factor SingleEle trigger SF_eleTrigger
          otree->trigweight_1 = (SF_lepTrigger->get_ScaleFactor(double(analysisTree.electron_pt[leptonIndex]),double(analysisTree.electron_eta[leptonIndex])));
                // Scale Factor Id+Iso SF_eleIdIso
          otree->idisoweight_1 = (SF_lepIdIso->get_ScaleFactor(double(analysisTree.electron_pt[leptonIndex]),double(analysisTree.electron_eta[leptonIndex])));
          otree->effweight = (otree->trigweight_1)*(otree->idisoweight_1)*(otree->trigweight_2)*(otree->idisoweight_2);
        }
	    }

      FillTau(&analysisTree, otree, tauIndex);

      //ditau sytem
      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
					  analysisTree.tau_py[tauIndex],
					  analysisTree.tau_pz[tauIndex],
					  tauMass);

      TLorentzVector metLV; metLV.SetXYZT(analysisTree.pfmet_ex,
					  analysisTree.pfmet_ey,
					  0,
					  TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey));

      TLorentzVector dileptonLV = leptonLV + tauLV;

      // visible mass
      otree->m_vis = dileptonLV.M();
      // visible ditau pt 
      otree->pt_tt = (dileptonLV+metLV).Pt();

      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;

      // dilepton veto
      if(ch=="mt") otree->dilepton_veto = dilepton_veto_mt(&cfg, &analysisTree);
			if(ch=="et") otree->dilepton_veto = dilepton_veto_et(&cfg, &analysisTree);

			//extra letpn veto
			otree->extraelec_veto = extra_electron_veto(leptonIndex, ch, &cfg, &analysisTree);
      otree->extramuon_veto = extra_muon_veto(leptonIndex, ch, &cfg, &analysisTree);

      //MET

      fillMET(ch, leptonIndex, tauIndex, &analysisTree, otree);

      // define MET covariance
/*      TMatrixD covMET(2, 2);
      covMET[0][0] = otree->mvacov00;
      covMET[1][0] = otree->mvacov10;
      covMET[0][1] = otree->mvacov01;
      covMET[1][1] = otree->mvacov11;*/

      // mt calculation and filling
      mt_calculation(otree);

      // bisector of lepton and tau transverse momenta

      float leptonUnitX = leptonLV.Px()/leptonLV.Pt();
      float leptonUnitY = leptonLV.Py()/leptonLV.Pt();

      float tauUnitX = tauLV.Px()/tauLV.Pt();
      float tauUnitY = tauLV.Py()/tauLV.Pt();

      float zetaX = leptonUnitX + tauUnitX;
      float zetaY = leptonUnitY + tauUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorVisX = leptonLV.Px() + tauLV.Px();
      float vectorVisY = leptonLV.Py() + tauLV.Py();

      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;
      otree->pzetamiss = otree->mvamet*TMath::Cos(otree->mvametphi)*zetaX + otree->mvamet*TMath::Sin(otree->mvametphi)*zetaY;
      otree->pfpzetamiss = analysisTree.pfmet_ex*zetaX + analysisTree.pfmet_ey*zetaY;      
      otree->puppipzetamiss = analysisTree.puppimet_ex*zetaX + analysisTree.puppimet_ey*zetaY;

      metLV.SetXYZT(otree->mvamet*TMath::Cos(otree->mvametphi),
                otree->mvamet*TMath::Sin(otree->mvametphi),
                0,
                TMath::Sqrt( otree->mvamet*TMath::Sin(otree->mvametphi)*otree->mvamet*TMath::Sin(otree->mvametphi) +
                         otree->mvamet*TMath::Cos(otree->mvametphi)*otree->mvamet*TMath::Cos(otree->mvametphi)));
      otree->pt_tt = (dileptonLV+metLV).Pt();   

      //counting jet
      counting_jets(&analysisTree, otree, &cfg);

      ////////////////////////////////////////////////////////////
      // MET Recoil Corrections
      ////////////////////////////////////////////////////////////

      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);

      otree->njetshad = otree->njets;
      if (!isData && applyRecoilCorrections){
				genV = genTools::genV(analysisTree);
				genL = genTools::genL(analysisTree);
				otree->njetshad = genTools::nJetsHad(analysisTree);
      }

      // MVA MET      
			// // njetshad, genVis, mean-resolution correction
			genTools::RecoilCorrections( *recoilMvaMetCorrector, (!isData && applyRecoilCorrections) * genTools::MeanResolution,
			                     otree->mvamet, otree->mvametphi,
			                     genV.Px(), genV.Py(),
			                     genL.Px(), genL.Py(),
			                     otree->njetshad,
			                     otree->mvamet_rcmr, otree->mvametphi_rcmr
			                     );
			             
			otree->mt_rcmr_1 = calc::mt(otree->pt_1, otree->phi_1, otree->mvamet_rcmr, otree->mvametphi_rcmr);
			otree->mt_rcmr_2 = calc::mt(otree->pt_2, otree->phi_2, otree->mvamet_rcmr, otree->mvametphi_rcmr);     
			otree->pzetamiss_rcmr = calc::pzetamiss( zetaX, zetaY, otree->mvamet_rcmr, otree->mvametphi_rcmr);

			//end MET Recoil Corrections

      otree->Fill();
      selEvents++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
   }

  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();
  file->Close();
  delete file;
    


}




///////////////////////////////////////////////
//////////////FUNCTION DEFINITION//////////////
///////////////////////////////////////////////

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

//accepts run number, lumi and json, make pair of (run,lumi) and starts isGoodLumi 
bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}

//fill the otree with the weights
void fill_weight(const AC1B * analysisTree, Spring15Tree *otree, PileUp *PUofficial, bool isData){
  float xs = -1.;
  
  otree->mcweight = 1.;
  otree->pu_weight = 0;

  if(!isData)
    otree->mcweight = analysisTree->genweight;

  otree->pu_weight = float(PUofficial->get_PUweight(double(analysisTree->numtruepileupinteractions)));

  otree->xs = xs;
  otree->trigweight_1 = 0;
  otree->trigweight_2 = 1;
  otree->idisoweight_1 = 0;
  otree->idisoweight_2 = 1;

  otree->effweight = 0;
  otree->fakeweight = 0;
  otree->embeddedWeight = 0;
  otree->signalWeight = 0;
  otree->weight = 1;
  otree->lheHt = analysisTree->genparticles_lheHt;
  otree->gen_noutgoing = analysisTree->genparticles_noutgoing;
}

//compute the absolute isolation for a given lepton labeled by Index in channel ch
float abs_Iso (int Index, TString ch, const AC1B * analysisTree, float dRiso){
  float neutralHadIso, photonIso, chargedHadIso, puIso;

  if(ch=="mt"){
    neutralHadIso = analysisTree->muon_neutralHadIso[Index];
    photonIso =     analysisTree->muon_photonIso[Index];
    chargedHadIso = analysisTree->muon_chargedHadIso[Index];
    puIso =         analysisTree->muon_puIso[Index];
    if (dRiso>0.29) {
      neutralHadIso =     analysisTree->muon_r03_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r03_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r03_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r03_sumPUPt[Index];
    }
    if (dRiso>0.39) {
      neutralHadIso =     analysisTree->muon_r04_sumNeutralHadronEt[Index];
      photonIso =         analysisTree->muon_r04_sumPhotonEt[Index];
      chargedHadIso =     analysisTree->muon_r04_sumChargedHadronPt[Index];
      puIso =             analysisTree->muon_r04_sumPUPt[Index];
    }
  }
  if(ch=="et"){
    neutralHadIso = analysisTree->electron_neutralHadIso[Index];
    photonIso =     analysisTree->electron_photonIso[Index];
    chargedHadIso = analysisTree->electron_chargedHadIso[Index];
    puIso =         analysisTree->electron_puIso[Index];
    if (dRiso>0.29) {
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
float rel_Iso(int Index, TString ch, const AC1B * analysisTree, float dRiso){
  if(ch=="mt")  return(abs_Iso(Index, ch, analysisTree, dRiso) / analysisTree->muon_pt[Index] );
  else if(ch=="et")   return(abs_Iso(Index, ch, analysisTree, dRiso) / analysisTree->electron_pt[Index] );
    else return(-1.);
}


////FILLING FUNCTIONS//////

//fill the otree with the muon variables in channel mutau
void FillMuTau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->muon_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->muon_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->muon_phi[leptonIndex];
  otree->m_1 = 		muonMass;
  otree->q_1 = -1;
  if (analysisTree->muon_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->muon_genmatch[leptonIndex];

	otree->iso_1 = abs_Iso(leptonIndex, "mt", analysisTree, dRiso) / analysisTree->muon_pt[leptonIndex];

	otree->d0_1 = analysisTree->muon_dxy[leptonIndex];
	otree->dZ_1 = analysisTree->muon_dz[leptonIndex];

	otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = 0;
	otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->againstElectronLooseMVA5_1 = 0;
	otree->againstElectronMediumMVA5_1 = 0;
	otree->againstElectronTightMVA5_1 = 0;
	otree->againstElectronVLooseMVA5_1 = 0;
	otree->againstElectronVTightMVA5_1 = 0;
	otree->againstElectronVLooseMVA6_2 = 0;
	otree->againstElectronTightMVA6_2 = 0;
	otree->againstMuonLoose3_1 = 0;
	otree->againstMuonTight3_1 = 0;

}

//fill the otree with the electron variables in channel etau
void FillETau(const AC1B * analysisTree, Spring15Tree *otree, int leptonIndex, float dRiso){
	otree->pt_1 = 	analysisTree->electron_pt[leptonIndex];
  otree->eta_1 = 	analysisTree->electron_eta[leptonIndex];
  otree->phi_1 = 	analysisTree->electron_phi[leptonIndex];
  otree->m_1 = 		electronMass;
    otree->q_1 = -1;
  if (analysisTree->electron_charge[leptonIndex]>0)
    otree->q_1 = 1;
  otree->gen_match_1 = analysisTree->electron_genmatch[leptonIndex];

  otree->iso_1 =  abs_Iso(leptonIndex, "et", analysisTree, dRiso) /analysisTree->electron_pt[leptonIndex];
  otree->mva_1 = analysisTree->electron_mva_value_nontrig_Spring15_v1[leptonIndex];

	otree->d0_1 = analysisTree->electron_dxy[leptonIndex];
	otree->dZ_1 = analysisTree->electron_dz[leptonIndex];

	otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = 0;
	otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
	otree->againstElectronLooseMVA5_1 = 0;
	otree->againstElectronMediumMVA5_1 = 0;
	otree->againstElectronTightMVA5_1 = 0;
	otree->againstElectronVLooseMVA5_1 = 0;
	otree->againstElectronVTightMVA5_1 = 0;

	otree->againstMuonLoose3_1 = 0;
	otree->againstMuonTight3_1 = 0;

}

//fill the otree with the tau variables 
void FillTau(const AC1B * analysisTree, Spring15Tree *otree, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->eta_2 = analysisTree->tau_eta[tauIndex];
  otree->phi_2 = analysisTree->tau_phi[tauIndex];
  otree->q_2 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_2 = analysisTree->tau_genmatch[tauIndex];
  //if (analysisTree->tau_charge[tauIndex]>0)
  //otree->q_2 = 1;
  //otree->mva_2 = log(0);
  otree->mva_2 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex]; //before, tau_byIsolationMVArun2v1DBoldDMwLTraw
  otree->d0_2 = analysisTree->tau_leadchargedhadrcand_dxy[tauIndex];
  otree->dZ_2 = analysisTree->tau_leadchargedhadrcand_dz[tauIndex];      
  otree->iso_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->m_2 = analysisTree->tau_mass[tauIndex];
  otree->tau_decay_mode_2 = analysisTree->tau_decayMode[tauIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstElectronLooseMVA5_2 = analysisTree->tau_againstElectronLooseMVA5[tauIndex];
  otree->againstElectronMediumMVA5_2 = analysisTree->tau_againstElectronMediumMVA5[tauIndex];
  otree->againstElectronTightMVA5_2 = analysisTree->tau_againstElectronTightMVA5[tauIndex];
  otree->againstElectronVLooseMVA5_2 = analysisTree->tau_againstElectronVLooseMVA5[tauIndex];
  otree->againstElectronVTightMVA5_2 = analysisTree->tau_againstElectronVTightMVA5[tauIndex];
  otree->againstMuonLoose3_2 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_2 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronVLooseMVA6_2 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronTightMVA6_2 = analysisTree->tau_againstElectronTightMVA6[tauIndex];

  otree->byTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree->tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];

}

//////DILEPTON FUNCTIONS

//returns the dilepton veto for mt channel
bool dilepton_veto_mt(const Config *cfg,const  AC1B *analysisTree){

  for (unsigned int im = 0; im<analysisTree->muon_count; ++im) {
		if (analysisTree->muon_pt[im]<=cfg->get<float>("ptDiMuonVeto")) continue;
		if (fabs(analysisTree->muon_eta[im])>=cfg->get<float>("etaDiMuonVeto")) continue;	
		
		if (fabs(analysisTree->muon_dxy[im])>=cfg->get<float>("dxyDiMuonVeto")) continue;
		if (fabs(analysisTree->muon_dz[im])>=cfg->get<float>("dzDiMuonVeto")) continue;

		float absIsoMu =  abs_Iso(im, "mt", analysisTree, cfg->get<float>("dRiso"));
		float relIsoMu = 	rel_Iso(im, "mt", analysisTree, cfg->get<float>("dRiso"));

		if(relIsoMu >= cfg->get<float>("isoDiMuonVeto")) continue;
		
		//bool passedVetoId =  analysisTree->muon_isMedium[im]; 
    bool passedVetoId = isICHEPmed(im, analysisTree);
		if (!passedVetoId && cfg->get<bool>("applyDiMuonVetoId")) continue;
		
		for (unsigned int je = im+1; je<analysisTree->muon_count; ++je) {
		  if (analysisTree->muon_pt[je]<=cfg->get<float>("ptDiMuonVeto")) continue;
		  if (fabs(analysisTree->muon_eta[je])>=cfg->get<float>("etaDiMuonVeto")) continue;	
		  
		  if (fabs(analysisTree->muon_dxy[je])>=cfg->get<float>("dxyDiMuonVeto")) continue;
		  if (fabs(analysisTree->muon_dz[je])>=cfg->get<float>("dzDiMuonVeto")) continue;
		  
		  if (analysisTree->muon_charge[im] * analysisTree->muon_charge[je] > 0. && cfg->get<bool>("applyDiMuonOS")) continue;

		  float absIsoMu =  abs_Iso(je, "mt", analysisTree, cfg->get<float>("dRiso"));
			float relIsoMu = 	rel_Iso(je, "mt", analysisTree, cfg->get<float>("dRiso"));

		  if(relIsoMu >= cfg->get<float>("isoDiMuonVeto")) continue;	

		  //passedVetoId =  analysisTree->muon_isMedium[je];
      passedVetoId = isICHEPmed(je, analysisTree);

		  if (!passedVetoId && cfg->get<bool>("applyDiMuonVetoId")) continue;
		  
		  float dr = deltaR(analysisTree->muon_eta[im],analysisTree->muon_phi[im],
				    analysisTree->muon_eta[je],analysisTree->muon_phi[je]);

		  if(dr<=cfg->get<float>("drDiMuonVeto")) continue;

		  return(1);
		}
  }
  return(0);
}


//returns the dilepton veto fot et channel
bool dilepton_veto_et(const Config *cfg,const  AC1B *analysisTree){

  for (unsigned int ie = 0; ie<analysisTree->electron_count; ++ie) {
		if (analysisTree->electron_pt[ie]<=cfg->get<float>("ptDiElectronVeto")) continue;
		if (fabs(analysisTree->electron_eta[ie])>=cfg->get<float>("etaDiElectronVeto")) continue;	
		
		if (fabs(analysisTree->electron_dxy[ie])>=cfg->get<float>("dxyDiElectronVeto")) continue;
		if (fabs(analysisTree->electron_dz[ie])>=cfg->get<float>("dzDiElectronVeto")) continue;

		float absIsoEle =   abs_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));
		float relIsoEle =   rel_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));


		if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;
		
		bool passedVetoId =  analysisTree->electron_cutId_veto_Spring15[ie];
		if (!passedVetoId && cfg->get<bool>("applyDiElectronVetoId")) continue;
		
		for (unsigned int je = ie+1; je<analysisTree->electron_count; ++je) {
		  if (analysisTree->electron_pt[je]<=cfg->get<float>("ptDiElectronVeto")) continue;
		  if (fabs(analysisTree->electron_eta[je])>=cfg->get<float>("etaDiElectronVeto")) continue;	
		  
		  if (fabs(analysisTree->electron_dxy[je])>=cfg->get<float>("dxyDiElectronVeto")) continue;
		  if (fabs(analysisTree->electron_dz[je])>=cfg->get<float>("dzDiElectronVeto")) continue;
		  
		  if (analysisTree->electron_charge[ie] * analysisTree->electron_charge[je] > 0. && cfg->get<bool>("applyDiElectronOS")) continue;

			float absIsoEle =  abs_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));
			float relIsoEle = 	rel_Iso(ie, "et", analysisTree, cfg->get<float>("dRiso"));

		  if(relIsoEle >= cfg->get<float>("isoDiElectronVeto")) continue;	

		  passedVetoId =  analysisTree->electron_cutId_veto_Spring15[je];
		  if (!passedVetoId && cfg->get<bool>("applyDiElectronVetoId")) continue;
		  
		  float dr = deltaR(analysisTree->electron_eta[ie],analysisTree->electron_phi[ie],
				    analysisTree->electron_eta[je],analysisTree->electron_phi[je]);

		  if(dr<=cfg->get<float>("drDiElectronVeto")) continue;

		  return(1);
		}
  }
  return(0);
}

//////EXTRA LEPTON VETO FUNCTIONS

//returns the extra electron veto
bool extra_electron_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree){
	for (unsigned int ie = 0; ie<analysisTree->electron_count; ++ie) {
		if (ch=="et") {if (int(ie)==leptonIndex) continue;}
		if (analysisTree->electron_pt[ie]<=cfg->get<float>("ptVetoElectronCut")) continue;
		if (fabs(analysisTree->electron_eta[ie])>=cfg->get<float>("etaVetoElectronCut")) continue;
		if (fabs(analysisTree->electron_dxy[ie])>=cfg->get<float>("dxyVetoElectronCut")) continue;
		if (fabs(analysisTree->electron_dz[ie])>=cfg->get<float>("dzVetoElectronCut")) continue;

		bool electronMvaId = analysisTree->electron_mva_wp90_nontrig_Spring15_v1[ie];
		if (!electronMvaId && cfg->get<bool>("applyVetoElectronId")) continue;
		if (!analysisTree->electron_pass_conversion[ie] && cfg->get<bool>("applyVetoElectronId")) continue;
		if (analysisTree->electron_nmissinginnerhits[ie]>1 && cfg->get<bool>("applyVetoElectronId")) continue;

		float relIsoEle = rel_Iso(ie, ch, analysisTree, cfg->get<float>("dRiso"));

		if (relIsoEle>=cfg->get<float>("isoVetoElectronCut")) continue;

		return(1);		
	}
	return(0);
}			

//returns the extra muon veto
bool extra_muon_veto(int leptonIndex, TString ch, const Config *cfg, const AC1B *analysisTree){
	for (unsigned int im = 0; im<analysisTree->muon_count; ++im) {
		if (ch=="mt") {if (int(im)==leptonIndex) continue;}
		if (analysisTree->muon_pt[im]<cfg->get<float>("ptVetoMuonCut")) continue;
		if (fabs(analysisTree->muon_eta[im])>cfg->get<float>("etaVetoMuonCut")) continue;
		if (fabs(analysisTree->muon_dxy[im])>cfg->get<float>("dxyVetoMuonCut")) continue;
		if (fabs(analysisTree->muon_dz[im])>cfg->get<float>("dzVetoMuonCut")) continue;
		//if (cfg->get<bool>("applyVetoMuonId") && !analysisTree->muon_isMedium[im]) continue;
    if (cfg->get<bool>("applyVetoMuonId") && !(isICHEPmed(im, analysisTree))) continue;

		float relIsoMu = rel_Iso(im, ch, analysisTree, cfg->get<float>("dRiso"));
		if (relIsoMu>cfg->get<float>("isoVetoMuonCut")) continue;

		return(1);
  }
  return(0);
}


//////MET FUNCTIONS

//fill the otree with the met variables
void fillMET(TString ch, int leptonIndex, int tauIndex, const AC1B * analysisTree, Spring15Tree *otree){

  // svfit variables
  otree->m_sv = -9999;
  otree->pt_sv = -9999;
  otree->eta_sv = -9999;
  otree->phi_sv = -9999;

   // pfmet variables
  otree->met = TMath::Sqrt(analysisTree->pfmet_ex*analysisTree->pfmet_ex + analysisTree->pfmet_ey*analysisTree->pfmet_ey);
  otree->metphi = TMath::ATan2(analysisTree->pfmet_ey,analysisTree->pfmet_ex);
  otree->metcov00 = analysisTree->pfmet_sigxx;
  otree->metcov01 = analysisTree->pfmet_sigxy;
  otree->metcov10 = analysisTree->pfmet_sigyx;
  otree->metcov11 = analysisTree->pfmet_sigyy;

  float met_x = analysisTree->pfmet_ex;
  float met_y = analysisTree->pfmet_ey;
  float met_x2 = met_x * met_x;
  float met_y2 = met_y * met_y;

  // puppimet variables
  otree->puppimet = TMath::Sqrt(analysisTree->puppimet_ex*analysisTree->puppimet_ex +
		    analysisTree->puppimet_ey*analysisTree->puppimet_ey);
  otree->puppimetphi = TMath::ATan2(analysisTree->puppimet_ey,analysisTree->puppimet_ex);

  // choosing mva met
/*  int mva_break;
  if (ch=="mt") mva_break = 3;
  	else if (ch=="et") mva_break = 2;

  unsigned int iMet = 0;
  float mvamet_x = 0;
  float mvamet_y = 0;
  otree->mvacov00 = 0.;
  otree->mvacov01 = 0.;
  otree->mvacov10 = 0.;
  otree->mvacov11 = 0.;

  for (; iMet<analysisTree->mvamet_count; ++iMet) {
		if ((int)analysisTree->mvamet_channel[iMet]==mva_break) break;
      }

  if (iMet>=analysisTree->mvamet_count){		
		otree->mvamet = 	log(0);
		otree->mvametphi = log(0);
	 	otree->mvacov00 = log(0);
		otree->mvacov01 = log(0);
		otree->mvacov10 = log(0);
		otree->mvacov11 = log(0);
  } else {
		// choosing mva met
		unsigned int iMet = 0;
		for (; iMet<analysisTree->mvamet_count; ++iMet) {
		  if (analysisTree->mvamet_channel[iMet]==mva_break){
		    if( ((int)analysisTree->mvamet_lep1[iMet])==tauIndex && ((int)analysisTree->mvamet_lep2[iMet])==leptonIndex)
		      break;
		  }
		}

		float mvamet_x = 0;
		float mvamet_y = 0;
		otree->mvacov00 = 0.;
		otree->mvacov01 = 0.;
		otree->mvacov10 = 0.;
		otree->mvacov11 = 0.;
		
		if(iMet < analysisTree->mvamet_count){
		  mvamet_x = analysisTree->mvamet_ex[iMet];
		  mvamet_y = analysisTree->mvamet_ey[iMet];
		  otree->mvacov00 = analysisTree->mvamet_sigxx[iMet];
		  otree->mvacov01 = analysisTree->mvamet_sigxy[iMet];
		  otree->mvacov10 = analysisTree->mvamet_sigyx[iMet];
		  otree->mvacov11 = analysisTree->mvamet_sigyy[iMet];
		}
		else{
		  std::cout<<"MVA MET not found! :: "<< analysisTree->mvamet_count << std::endl;
		  iMet = 0;
		  //std::cout<<"tau = "<<tauIndex<<" mu = "<<muonIndex<<std::endl;
		  //for (; iMet<analysisTree->mvamet_count; ++iMet) {
		    //std::cout<<"imet = "<<analysisTree->mvamet_channel[iMet]<<" itau = "<<analysisTree->mvamet_lep1[iMet]<<" imu = "<<analysisTree->mvamet_lep2[iMet]<<std::endl;
		  //}
		  
		}
		
		float mvamet_x2 = mvamet_x * mvamet_x;
		float mvamet_y2 = mvamet_y * mvamet_y;
		otree->mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
		otree->mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
  }*/

}


//caltulates and fill the otree with the MT variables
void mt_calculation(Spring15Tree *otree){
  otree->mt_1 = sqrt(2*otree->pt_1*otree->mvamet*(1.-cos(otree->phi_1-otree->mvametphi)));
  otree->mt_2 = sqrt(2*otree->pt_2*otree->mvamet*(1.-cos(otree->phi_2-otree->mvametphi)));
  
  otree->pfmt_1 = sqrt(2*otree->pt_1*otree->met*(1.-cos(otree->phi_1-otree->metphi)));
  otree->pfmt_2 = sqrt(2*otree->pt_2*otree->met*(1.-cos(otree->phi_2-otree->metphi)));

  otree->puppimt_1 = sqrt(2*otree->pt_1*otree->puppimet*(1.-cos(otree->phi_1-otree->puppimetphi)));
  otree->puppimt_2 = sqrt(2*otree->pt_2*otree->puppimet*(1.-cos(otree->phi_2-otree->puppimetphi)));
}

//filling otree with bisector variables


void counting_jets(const AC1B *analysisTree, Spring15Tree *otree, const Config *cfg){

  vector<unsigned int> jets; jets.clear();
  vector<unsigned int> jetspt20; jetspt20.clear();
  vector<unsigned int> bjets; bjets.clear();

  int indexLeadingJet = -1;
  float ptLeadingJet = -1;

  int indexSubLeadingJet = -1;
  float ptSubLeadingJet = -1;

  int indexLeadingBJet = -1;
  float ptLeadingBJet = -1;

  for (unsigned int jet=0; jet<analysisTree->pfjet_count; ++jet) {

    float absJetEta = fabs(analysisTree->pfjet_eta[jet]);
    if (absJetEta>=cfg->get<float>("JetEtaCut")) continue;

    float jetPt = analysisTree->pfjet_pt[jet];
    if (jetPt<=cfg->get<float>("JetPtLowCut")) continue;

    float dR1 = deltaR(analysisTree->pfjet_eta[jet],analysisTree->pfjet_phi[jet],
                             otree->eta_1,otree->phi_1);
    if (dR1<=cfg->get<float>("dRJetLeptonCut")) continue;

    float dR2 = deltaR(analysisTree->pfjet_eta[jet],analysisTree->pfjet_phi[jet],
                             otree->eta_2,otree->phi_2);
    if (dR2<=cfg->get<float>("dRJetLeptonCut")) continue;

    // jetId
    float energy = analysisTree->pfjet_e[jet];
          energy *= analysisTree->pfjet_energycorr[jet];
          float chf = analysisTree->pfjet_chargedhadronicenergy[jet]/energy;
          float nhf = analysisTree->pfjet_neutralhadronicenergy[jet]/energy;
          float phf = analysisTree->pfjet_neutralemenergy[jet]/energy;
          float elf = analysisTree->pfjet_chargedemenergy[jet]/energy;
    float muf = analysisTree->pfjet_muonenergy[jet]/energy;
          float chm = analysisTree->pfjet_chargedmulti[jet];
    float nm = analysisTree->pfjet_neutralmulti[jet];
          float npr = analysisTree->pfjet_chargedmulti[jet] + analysisTree->pfjet_neutralmulti[jet];
    //bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
    bool isPFJetId = false;
    if (absJetEta<=3.0)
      isPFJetId = (nhf < 0.99 && phf < 0.99 && npr > 1) && (absJetEta>2.4 || (chf>0 && chm > 0 && elf < 0.99));
    else
      isPFJetId = phf < 0.9 && nm > 10;
    //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99 && muf < 0.8) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
    //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
    
    if (!isPFJetId) continue;

    jetspt20.push_back(jet);

    if (absJetEta<cfg->get<float>("bJetEtaCut") && analysisTree->pfjet_btag[jet][0]>cfg->get<float>("btagCut")) { // b-jet
      bjets.push_back(jet);
      if (jetPt>ptLeadingBJet) {
        ptLeadingBJet = jetPt;
        indexLeadingBJet = jet;
      }
    } 

    if (indexLeadingJet>=0) {
      if (jetPt<ptLeadingJet && jetPt>ptSubLeadingJet) {
        indexSubLeadingJet = jet;
        ptSubLeadingJet = jetPt;
      }
    }

    if (jetPt>ptLeadingJet) {
      indexLeadingJet = jet;
      ptLeadingJet = jetPt;
    }

    if (jetPt<cfg->get<float>("JetPtHighCut")) continue;
    jets.push_back(jet);
  }

  otree->njets = jets.size();
  otree->njetspt20 = jetspt20.size();
  otree->nbtag = bjets.size();
  
  otree->bpt = -9999;
  otree->beta = -9999;
  otree->bphi = -9999;
  
  if (indexLeadingBJet>=0) {
    otree->bpt = analysisTree->pfjet_pt[indexLeadingBJet];
    otree->beta = analysisTree->pfjet_eta[indexLeadingBJet];
    otree->bphi = analysisTree->pfjet_phi[indexLeadingBJet];
  }

  otree->jpt_1 = -9999;
  otree->jeta_1 = -9999;
  otree->jphi_1 = -9999;
  otree->jptraw_1 = -9999;
  otree->jptunc_1 = -9999;
  otree->jmva_1 = -9999;
  otree->jlrm_1 = -9999;
  otree->jctm_1 = -9999;

  if ( indexLeadingJet>=0 && indexSubLeadingJet>=0 && indexLeadingJet==indexSubLeadingJet )
cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;

  if (indexLeadingJet>=0) {
    otree->jpt_1 = analysisTree->pfjet_pt[indexLeadingJet];
    otree->jeta_1 = analysisTree->pfjet_eta[indexLeadingJet];
    otree->jphi_1 = analysisTree->pfjet_phi[indexLeadingJet];
    otree->jptraw_1 = analysisTree->pfjet_pt[indexLeadingJet]*analysisTree->pfjet_energycorr[indexLeadingJet];
    otree->jmva_1 = analysisTree->pfjet_pu_jet_full_mva[indexLeadingJet];
  }

  otree->jpt_2 = -9999;
  otree->jeta_2 = -9999;
  otree->jphi_2 = -9999;
  otree->jptraw_2 = -9999;
  otree->jptunc_2 = -9999;
  otree->jmva_2 = -9999;
  otree->jlrm_2 = -9999;
  otree->jctm_2 = -9999;

  if (indexSubLeadingJet>=0) {
    otree->jpt_2 = analysisTree->pfjet_pt[indexSubLeadingJet];
    otree->jeta_2 = analysisTree->pfjet_eta[indexSubLeadingJet];
    otree->jphi_2 = analysisTree->pfjet_phi[indexSubLeadingJet];
    otree->jptraw_2 = analysisTree->pfjet_pt[indexSubLeadingJet]*analysisTree->pfjet_energycorr[indexSubLeadingJet];
    otree->jmva_2 = analysisTree->pfjet_pu_jet_full_mva[indexSubLeadingJet];
  }

  otree->mjj =  -9999;
  otree->jdeta =  -9999;
  otree->njetingap = -1;

  if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
    otree->njetingap = 0;
    TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree->pfjet_px[indexLeadingJet],
                 analysisTree->pfjet_py[indexLeadingJet],
                 analysisTree->pfjet_pz[indexLeadingJet],
                 analysisTree->pfjet_e[indexLeadingJet]);

    TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree->pfjet_px[indexSubLeadingJet],
                 analysisTree->pfjet_py[indexSubLeadingJet],
                 analysisTree->pfjet_pz[indexSubLeadingJet],
                 analysisTree->pfjet_e[indexSubLeadingJet]);

    otree->mjj = (jet1+jet2).M();
    otree->jdeta = abs(analysisTree->pfjet_eta[indexLeadingJet]-
          analysisTree->pfjet_eta[indexSubLeadingJet]);
   
    float etamax = analysisTree->pfjet_eta[indexLeadingJet];
    float etamin = analysisTree->pfjet_eta[indexSubLeadingJet];
    if (etamax<etamin) {
      float tmp = etamax;
      etamax = etamin;
      etamin = tmp;
    }
    for (unsigned int jet=0; jet<jets.size(); ++jet) {
      int index = jets.at(jet);
      float etaX = analysisTree->pfjet_eta[index];
      if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax) 
        otree->njetingap++;
    }


  }

}
