#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

    const bool isData = cfg.get<bool>("IsData");
    const bool isDY = cfg.get<bool>("IsDY");
    const bool isWJets = cfg.get<bool>("IsWJets");
    const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

    const string jsonFile = cfg.get<string>("jsonFile");
    const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
    const string recoilFileName   = cfg.get<string>("RecoilFileName");
    const string dataPUFile = cfg.get<string>("DataPUFile");
    const string mcPUFile = cfg.get<string>("MCPUFile");
    const string sampleName = cfg.get<string>("sampleName");
    // pile up reweighting
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
    
    // which Higgs final states: TT, MT, ET, EM, MM, EE
    const string HiggsFinalState = cfg.get<string>("HiggsFinalState");
    
    // kinematic cuts on electrons
    const float ptEleCut = cfg.get<float>("ptEleCut");
    const float etaEleCut = cfg.get<float>("etaEleCut");
    const float dxyEleCut = cfg.get<float>("dxyEleCut");
    const float dzEleCut = cfg.get<float>("dzEleCut");

    // kinematic cuts on muons
    const float ptMuonCut  = cfg.get<float>("ptMuonCut");
    const float etaMuonCut = cfg.get<float>("etaMuonCut");
    const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
    const float dzMuonCut      = cfg.get<float>("dzMuonCut");

    //kinematic cuts on taus
    const float ptTauCut = cfg.get<float>("ptTauCut");
    const float etaTauCut = cfg.get<float>("etaTauCut");
    const float dzTauCut = cfg.get<float>("dzTauCut");
    
    // jets
    const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
    const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
    const float jetEtaCut = cfg.get<float>("JetEtaCut");
    const float btagCut = cfg.get<float>("btagCut");
    const string bTagDiscriminator = cfg.get<string>("BTagDiscriminator");
    const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
    
    // topological cuts
    const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
    const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
    const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");

    // trigger
    const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
    const string muonTriggerName  = cfg.get<string>("MuonTriggerName");
    const string muonFilterName   = cfg.get<string>("MuonFilterName");

    const string singleMuonFilterName = cfg.get<string>("SingleMuonFilterName");
    const float singleMuonTriggerPtCut = cfg.get<float>("SingleMuonTriggerPtCut");
    const float singleMuonTriggerEtaCut = cfg.get<float>("SingleMuonTriggerEtaCut");

    TString MuonTriggerName(muonTriggerName);
    TString MuonFilterName(muonFilterName);

    TString SingleMuonFilterName(singleMuonFilterName);

    // vertex cuts
    const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
    const float zVertexCut     = cfg.get<float>("ZVertexCut");
    const float dVertexCut     = cfg.get<float>("DVertexCut");
    
    //scale factor
    const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoFile");
    const string MuonTriggerFile = cfg.get<string>("MuonTriggerFile");
    
    // Z pt weight
    const string ZptweightFile = cfg.get<string>("ZptweightFile");


    // **** end of configuration

    string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

    // Run-lumi selector
    std::vector<Period> periods;
    if (isData)
    { // read the good runs
	  std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
  	  if (inputFileStream.fail() )
      {
            std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
            std::cout << "please check" << std::endl;
            std::cout << "quitting program" << std::endl;
	    exit(-1);
      }
  
        for(std::string s; std::getline(inputFileStream, s); )
        {
           periods.push_back(Period());
           std::stringstream ss(s);
           ss >> periods.back();
        }
    }

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  TH1::SetDefaultSumw2(true);
    
    RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
    recoilPFMetCorrector = new RecoilCorrector(recoilFileName);

  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
    
    //declare the tree
    TTree * SynTree = new TTree("SynTree","SynTree");
    
    TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);

    Float_t         puweight;
    
    //Event ID variable
    ULong64_t run;
    UInt_t lumi;
    ULong64_t evt;
    
    //Pile up variable
    Int_t           npu;
    Int_t           npv;
    Float_t         rho;
    
    //Jet related variables
    Int_t   nbtag;
    Int_t   njets;
    Int_t   njetspt20;
    Float_t jpt_1;
    Float_t jeta_1;
    Float_t jphi_1;
    Float_t jpt_2;
    Float_t jeta_2;
    Float_t jphi_2;
    Float_t bpt_1;
    Float_t beta_1;
    Float_t bphi_1;
    Float_t bhadronFlavor_1;
    Float_t bpt_2;
    Float_t beta_2;
    Float_t bphi_2;
    Float_t bhadronFlavor_2;
    
    //Leg 1 (highest pt elec/mu from Z candidate)
    Float_t pt_1;
    Float_t phi_1;
    Float_t eta_1;
    Float_t iso_1;
    Int_t gen_match_1;
    Float_t IdRawMva_1;
    Float_t pfmt_1;
    Float_t muonSF_1;
    Float_t electronSF_1;

    //Leg 2 (subleading pt elec/mu from Z candidate)
    Float_t pt_2;
    Float_t phi_2;
    Float_t eta_2;
    Float_t iso_2;
    Int_t gen_match_2;
    Float_t IdRawMva_2;
    Float_t pfmt_2;
    Float_t muonSF_2;
    Float_t electronSF_2;
    
    //Leg 3 (Based off of Higgs legs: ET -> E, MT -> M, EM -> E, TT -> highest pt tau)
    Float_t pt_3;
    Float_t phi_3;
    Float_t eta_3;
    Float_t iso_3;
    Int_t gen_match_3;
    Float_t IdRawMva_3;
    Float_t pfmt_3;
    Float_t muonSF_3;
    Float_t electronSF_3;
    Float_t tauSF_3;
    
    //Leg 4 (Based off of Higgs legs: ET -> T, MT -> T, EM -> M, TT -> subleading pt tau)
    Float_t pt_4;
    Float_t phi_4;
    Float_t eta_4;
    Float_t iso_4;
    Int_t gen_match_4;
    Float_t IdRawMva_4;
    Float_t pfmt_4;
    Float_t muonSF_4;
    Float_t electronSF_4;
    Float_t tauSF_4;
    Float_t met;
    Float_t metphi;
    Float_t metcov00;
    Float_t metcov01;
    Float_t metcov10;
    Float_t metcov11;
    
    Float_t DR_12;
    Float_t DR_13;
    Float_t DR_14;
    
    Float_t DR_23;
    Float_t DR_24;
    Float_t DR_34;

    Float_t mll;
    Float_t Z_Pt;
    Float_t Z_DR;
    Bool_t  Z_SS;
    
    Float_t m_vis;
    Float_t pt_tt;
    Float_t H_DR;
    Bool_t  H_SS;
    
    Float_t Mass;

    UInt_t  npartons;

    SynTree->Branch("run",&run,"run/l");
    SynTree->Branch("lumi",&lumi,"lumi/i");
    SynTree->Branch("evt",&evt,"evt/l");
    SynTree->Branch("npv",&npv,"npv/I");
    SynTree->Branch("npu",&npu,"npu/I");
    SynTree->Branch("rho",&rho,"rho/F");
    SynTree->Branch("puweight", &puweight, "puweight/F");
    SynTree->Branch("nbtag",&nbtag,"nbtag/I");
    SynTree->Branch("njets",&njets,"njets/I");
    SynTree->Branch("njetspt20",&njetspt20,"njetspt20/I");
    SynTree->Branch("jpt_1",&jpt_1,"jpt_1/F");
    SynTree->Branch("jeta_1",&jeta_1,"jeta_1/F");
    SynTree->Branch("jphi_1",&jphi_1,"jphi_1/F");
    SynTree->Branch("jpt_2",&jpt_2,"jpt_2/F");
    SynTree->Branch("jeta_2",&jeta_2,"jeta_2/F");
    SynTree->Branch("jphi_2",&jphi_2,"jphi_2/F");
    SynTree->Branch("bpt_1",&bpt_1,"bpt_1/F");
    SynTree->Branch("beta_1",&beta_1,"beta_1/F");
    SynTree->Branch("bphi_1",&bphi_1,"bphi_1/F");
    SynTree->Branch("bhadronFlavor_1",&bhadronFlavor_1,"bhadronFlavor_1/F");
    SynTree->Branch("bpt_2",&bpt_2,"bpt_2/F");
    SynTree->Branch("beta_2",&beta_2,"beta_2/F");
    SynTree->Branch("bphi_2",&bphi_2,"bphi_2/F");
    SynTree->Branch("bhadronFlavor_2",&bhadronFlavor_2,"bhadronFlavor_2/F");
    
    SynTree->Branch("pt_1",&pt_1,"pt_1/F");
    SynTree->Branch("phi_1",&phi_1,"phi_1/F");
    SynTree->Branch("eta_1",&eta_1,"eta_1/F");
    SynTree->Branch("iso_1",&iso_1,"iso_1/F");
    SynTree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");
    SynTree->Branch("IdRawMva_1",&IdRawMva_1,"IdRawMva_1/F");
    SynTree->Branch("pfmt_1",&pfmt_1,"pfmt_1/F");
    SynTree->Branch("muonSF_1",&muonSF_1,"muonSF_1/F");
    SynTree->Branch("electronSF_1",&electronSF_1,"electronSF_1/F");
    
    SynTree->Branch("pt_2",&pt_2,"pt_2/F");
    SynTree->Branch("phi_2",&phi_2,"phi_2/F");
    SynTree->Branch("eta_2",&eta_2,"eta_2/F");
    SynTree->Branch("iso_2",&iso_2,"iso_2/F");
    SynTree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");
    SynTree->Branch("IdRawMva_2",&IdRawMva_2,"IdRawMva_2/F");
    SynTree->Branch("pfmt_2",&pfmt_2,"pfmt_2/F");
    SynTree->Branch("muonSF_2",&muonSF_2,"muonSF_2/F");
    SynTree->Branch("electronSF_2",&electronSF_2,"electronSF_2/F");
    
    SynTree->Branch("pt_3",&pt_3,"pt_3/F");
    SynTree->Branch("phi_3",&phi_3,"phi_3/F");
    SynTree->Branch("eta_3",&eta_3,"eta_3/F");
    SynTree->Branch("iso_3",&iso_3,"iso_3/F");
    SynTree->Branch("gen_match_3",&gen_match_3,"gen_match_3/I");
    SynTree->Branch("IdRawMva_3",&IdRawMva_3,"IdRawMva_3/F");
    SynTree->Branch("pfmt_3",&pfmt_3,"pfmt_3/F");
    SynTree->Branch("muonSF_3",&muonSF_3,"muonSF_3/F");
    SynTree->Branch("electronSF_3",&electronSF_3,"electronSF_3/F");
    SynTree->Branch("tauSF_3",&tauSF_3,"tauSF_3/F");
    
    SynTree->Branch("pt_4",&pt_4,"pt_4/F");
    SynTree->Branch("phi_4",&phi_4,"phi_4/F");
    SynTree->Branch("eta_4",&eta_4,"eta_4/F");
    SynTree->Branch("iso_4",&iso_4,"iso_4/F");
    SynTree->Branch("gen_match_4",&gen_match_4,"gen_match_4/I");
    SynTree->Branch("IdRawMva_4",&IdRawMva_4,"IdRawMva_4/F");
    SynTree->Branch("pfmt_4",&pfmt_4,"pfmt_4/F");
    SynTree->Branch("muonSF_4",&muonSF_4,"muonSF_4/F");
    SynTree->Branch("electronSF_4",&electronSF_4,"electronSF_4/F");
    SynTree->Branch("tauSF_4",&tauSF_4,"tauSF_4/F");
    SynTree->Branch("met",&met,"met/F");
    SynTree->Branch("metphi",&metphi,"metphi/F");
    SynTree->Branch("metcov00",&metcov00,"metcov00/F");
    SynTree->Branch("metcov01",&metcov01,"metcov01/F");
    SynTree->Branch("metcov10",&metcov10,"metcov10/F");
    SynTree->Branch("metcov11",&metcov11,"metcov11/F");
    SynTree->Branch("DR_12",&DR_12,"DR_12/F");
    SynTree->Branch("DR_13",&DR_13,"DR_13/F");
    SynTree->Branch("DR_14",&DR_14,"DR_14/F");
    SynTree->Branch("DR_23",&DR_23,"DR_23/F");
    SynTree->Branch("DR_24",&DR_24,"DR_24/F");
    SynTree->Branch("DR_34",&DR_34,"DR_34/F");

    SynTree->Branch("mll",&mll,"mll/F");
    SynTree->Branch("Z_Pt",&Z_Pt,"Z_Pt/F");
    SynTree->Branch("Z_DR",&Z_DR,"Z_DR/F");
    SynTree->Branch("Z_SS",&Z_SS,"Z_SS/O");

    SynTree->Branch("m_vis",&m_vis,"m_vis/F");
    SynTree->Branch("pt_tt",&pt_tt,"pt_tt/F");
    SynTree->Branch("H_DR",&H_DR,"H_DR/F");
    SynTree->Branch("H_SS",&H_SS,"H_SS/O");

    SynTree->Branch("Mass",&Mass,"Mass/F");
    SynTree->Branch("npartons",&npartons,"npartons/i");

    TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
    TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

    // initialize pile up object
    PileUp * PUofficial = new PileUp();
    
    TString puHistName(sampleName);//for 2017 only
    puHistName += "_pileup";
  
    if (applyPUreweighting)
    {
        TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(dataPUFile),"read");
        TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(mcPUFile), "read");
        TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
        TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
        //TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(puHistName); //for 2017 only
        PUofficial->set_h_data(PU_data);
        PUofficial->set_h_MC(PU_mc);
    }
    // Muon scale factors
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerFile));
    
    // Zpt reweighting for LO DY samples
    TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+ZptweightFile,"read");
    TH2D * h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

    
  int nFiles = 0;
  int nEvents = 0;
  int selEventsIsoMuons = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();

  for (int iF=0; iF<nTotalFiles; ++iF) {
  
    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
      
      TTree * _inittree = NULL;
      _inittree = (TTree*)file_->Get(TString(initNtupleName));
      if (_inittree!=NULL) {
          Float_t genweight;
          if (!isData)
              _inittree->SetBranchAddress("genweight",&genweight);
          Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
          std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
          for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
              _inittree->GetEntry(iEntry);
              if (isData)
                  histWeightsH->Fill(0.,1.);
              else
                  histWeightsH->Fill(0.,genweight);
          }
      }
      
      TTree * _tree = NULL;
      _tree = (TTree*)file_->Get(TString(ntupleName));
      if (_tree==NULL) continue;
      Long64_t numberOfEntries = _tree->GetEntries();
      std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
      AC1B analysisTree(_tree);
    
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
      
    int NE = int(histoInputEvents->GetEntries());
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
    std::cout << "      number of input events         = " << NE << std::endl;

    // EVENT LOOP //
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
          cout << "      processed " << nEvents << " events" << endl;

        float weight = 1;
        
        // Pileup variable
        npv = analysisTree.primvertex_count;
        npu = analysisTree.numtruepileupinteractions;
        rho = analysisTree.rho;
        
        if (!isData)
          weight *=analysisTree.genweight;
        npartons = analysisTree.genparticles_noutgoing;

        
        if (!isData) {
            //	cout << analysisTree.numtruepileupinteractions << endl;
            if (applyPUreweighting) {
                nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
                double Ninteractions = analysisTree.numtruepileupinteractions;
                double PUweight = PUofficial->get_PUweight(Ninteractions);
                weight *= float(PUweight);
                PUweightsOfficialH->Fill(PUweight);
                //	  cout << PUweight << endl;
            }
        }

      if (isData && applyGoodRunSelection){
          bool lumi = false;
          int n=analysisTree.event_run;
          int lum = analysisTree.event_luminosityblock;
	
          std::string num = std::to_string(n);
          std::string lnum = std::to_string(lum);
          for(const auto& a : periods)
          {
	    
              if ( num.c_str() ==  a.name ) {
                  //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
                  //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      
                  for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
                      //	cout<<b->lower<<"  "<<b->bigger<<endl;
                      if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
                  }
                  auto last = std::prev(a.ranges.end());
                  //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
                  if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
              }
          }
    
          if (!lumi) continue;
          //if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
          //std::remove("myinputfile");
      }     

      if (analysisTree.event_run<RunMin)
          RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
          RunMax = analysisTree.event_run;
      
      //std::cout << " Run : " << analysisTree.event_run << std::endl;
      
      bool isNewRun = true;
      if (allRuns.size()>0) {
          for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
              if (analysisTree.event_run==allRuns.at(iR)) {
                  isNewRun = false;
                  break;
              }
          }
      }
        
        run = analysisTree.event_run;
        lumi = analysisTree.event_luminosityblock;
        evt = analysisTree.event_nr;
        
        //pile up weight variable
        if (!isData) {
            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
        }
        
      if (isNewRun) 
          allRuns.push_back(analysisTree.event_run);
      
      bool isTriggerMuon = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
          TString trigName(it->first);
          if (trigName.Contains(MuonTriggerName)) {
              //  std::cout << it->first << " : " << it->second << std::endl;
              if (it->second==1)
                  isTriggerMuon = true;
          }
      }

      if (applyTrigger && !isTriggerMuon) continue;
    
      unsigned int nMuonFilter = 0;
      bool isMuonFilter = false;

      unsigned int nSingleMuonFilter = 0;
      bool isSingleMuonFilter = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MuonFilterName) {
	  nMuonFilter = i;
	  isMuonFilter = true;
	}
	if (HLTFilter==SingleMuonFilterName) {
          nSingleMuonFilter = i;
          isSingleMuonFilter = true;
        }
      }
      if (!isMuonFilter && isData) {
	cout << "Filter " << MuonFilterName << " not found " << endl;
	exit(-1);
      }
      if (!isSingleMuonFilter && isData) {
	cout << "Filter " << SingleMuonFilterName << " not found " << endl;
        exit(-1);
      }
      
        if(!((HiggsFinalState== "MT")||(HiggsFinalState== "ET")||(HiggsFinalState== "EM")||(HiggsFinalState== "EE")||(HiggsFinalState== "MM")||(HiggsFinalState== "TT")))
        {
            std::cout << "you are not considering one of the final states of Higgs decay"<< std::endl;
            exit(-1);
        }
        
        // vertex cuts
        if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
        if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
        float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);
        if (dVertex>dVertexCut) continue;
        
        //initiation of btag
        unsigned int nBTagDiscriminant = 0;
        for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag)
        {
            TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
            if (discr.Contains(bTagDiscriminator))
                nBTagDiscriminant = iBTag;
        }
        
        // electron selection
        vector<unsigned int> idEles; idEles.clear(); //identified electrons
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
            if (analysisTree.electron_pt[ie]<ptEleCut) continue;
            if (fabs(analysisTree.electron_eta[ie])>etaEleCut) continue;
            if (fabs(analysisTree.electron_dxy[ie])>dxyEleCut) continue;
            if (fabs(analysisTree.electron_dz[ie])>dzEleCut) continue;
            if (!analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie]) continue;
            if (!analysisTree.electron_pass_conversion[ie]) continue;
            if (analysisTree.electron_nmissinghits[ie]>=2) continue;
            idEles.push_back(ie);
            //    cout << "pt:" << analysisTree.electron_pt[ie] << "  passed:" << elePassed << endl;
        }
        
        // muon selection
        vector<unsigned int> idMuons; idMuons.clear(); //identified muons
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
            if (fabs(analysisTree.muon_pt[im])<ptMuonCut) continue;
            if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
            if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
            if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
            if (!(analysisTree.muon_isLoose[im])) continue;
            if (!(analysisTree.muon_isGlobal[im]||analysisTree.muon_isTracker[im])) continue;
            idMuons.push_back(im);
            //	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
        }
        
        //tau selection
        vector<unsigned int> idTaus; idTaus.clear(); //identified taus
        for(unsigned int it=0;it<analysisTree.tau_count;++it)
        {
            if (analysisTree.tau_pt[it]<ptTauCut) continue;
            if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
            if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
            if (analysisTree.tau_decayModeFindingNewDMs[it]<=0.5) continue;
            if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
            //(no needed for Synchronization)
            //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[it] <= 0.5) continue;
            if (analysisTree.tau_byLooseDeepTau2017v2p1VSmu[it] <= 0.5) continue;
            if (analysisTree.tau_byVLooseDeepTau2017v2p1VSe[it] <= 0.5) continue;
            idTaus.push_back(it);
        }
        
        //Final state selection: EE
        if(idEles.size() < 2) continue;
        if(HiggsFinalState == "EE")
        {
            if(idEles.size() != 4) continue;
            if(idMuons.size() != 0) continue;
        }
        if(HiggsFinalState == "EM")
        {
            if(idEles.size() != 3) continue;
            if(idMuons.size() != 1) continue;
        }
        if(HiggsFinalState == "ET")
        {
            if(idEles.size() != 3) continue;
            if(idMuons.size() != 0) continue;
        }
        if(HiggsFinalState == "MT")
        {
            if(idEles.size() != 2) continue;
            if(idMuons.size() != 1) continue;
        }
        if(HiggsFinalState == "MM")
        {
            if(idEles.size() != 2) continue;
            if(idMuons.size() != 2) continue;
        }
        if(HiggsFinalState == "TT")
        {
            if(idEles.size() != 2) continue;
            if(idMuons.size() != 0) continue;
        }
        //std::cout << "id electrons:" << idEles.size() << std::endl;
        //std::cout << "id muons:" << idMuons.size() << std::endl;
        //std::cout << "id tau:" << idTaus.size() << std::endl;
        
        //Check sepration among the leptons----------
        //First check electrons
        vector<unsigned int> goodEles; goodEles.clear(); //good electrons
        for(unsigned int ie=0;ie<idEles.size();++ie)
        {
            int indexElectron = idEles[ie];
            bool foundNearElectron = false;
            bool foundNearMuon = false;
            bool foundNearTau = false;
            //checking whether there are other electrons overlapped with the electron
            for(unsigned int ie2=0;ie2<idEles.size();++ie2)
            {
                int indexElectron2 = idEles[ie2];
                if (indexElectron == indexElectron2) continue;
                float dRElectrons = deltaR(analysisTree.electron_eta[indexElectron],analysisTree.electron_phi[indexElectron],analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2]);
                if (dRElectrons<0.3) foundNearElectron = true;
            }
            //checking whether there are other muons overlapped with the electron
            for(unsigned int im=0;im<idMuons.size();++im)
            {
                int indexMuon = idMuons[im];
                float dRElectronMuon = deltaR(analysisTree.electron_eta[indexElectron],analysisTree.electron_phi[indexElectron],analysisTree.muon_eta[indexMuon],analysisTree.muon_phi[indexMuon]);
                if (dRElectronMuon<0.3) foundNearMuon = true;
            }
            //checking whether there are other taus overlapped with the electron
            for(unsigned int it=0;it<idTaus.size();++it)
            {
                int indexTau = idTaus[it];
                float dRElectronTau = deltaR(analysisTree.electron_eta[indexElectron],analysisTree.electron_phi[indexElectron],analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau]);
                if (dRElectronTau<0.3) foundNearTau = true;
            }
            if ((!foundNearElectron)&&(!foundNearMuon)&&(!foundNearTau))
                goodEles.push_back(indexElectron);
        }
        
        //Second check muons
        vector<unsigned int> goodMuons; goodMuons.clear(); //good muons
        for(unsigned int im=0;im<idMuons.size();++im)
        {
            int indexMuon = idMuons[im];
            bool foundNearElectron = false;
            bool foundNearMuon = false;
            bool foundNearTau = false;
            //checking whether there are other electrons overlapped with the muon
            for(unsigned int ie=0;ie<idEles.size();++ie)
            {
                int indexElectron = idEles[ie];
                float dRMuonElectron = deltaR(analysisTree.muon_eta[indexMuon],analysisTree.muon_phi[indexMuon],analysisTree.electron_eta[indexElectron],analysisTree.electron_phi[indexElectron]);
                if (dRMuonElectron<0.3) foundNearElectron = true;
            }
            //checking whether there are other muons overlapped with the muon
            for(unsigned int im2=0;im2<idMuons.size();++im2)
            {
                int indexMuon2 = idMuons[im2];
                if (indexMuon == indexMuon2) continue;
                float dRMuons = deltaR(analysisTree.muon_eta[indexMuon],analysisTree.muon_phi[indexMuon],analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2]);
                if (dRMuons<0.3) foundNearMuon = true;
            }
            //checking whether there are other taus overlapped with the electron
            for(unsigned int it=0;it<idTaus.size();++it)
            {
                int indexTau = idTaus[it];
                float dRMuonTau = deltaR(analysisTree.muon_eta[indexMuon],analysisTree.muon_phi[indexMuon],analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau]);
                if (dRMuonTau<0.3) foundNearTau = true;
            }
            if ((!foundNearElectron)&&(!foundNearMuon)&&(!foundNearTau))
                goodMuons.push_back(indexMuon);
        }
        
        //Third check taus
        vector<unsigned int> goodTaus; goodTaus.clear(); //good taus
        for(unsigned int it=0;it<idTaus.size();++it)
        {
            int indexTau = idTaus[it];
            bool foundNearElectron = false;
            bool foundNearMuon = false;
            bool foundNearTau = false;
            //checking whether there are other electrons overlapped with the tau
            for(unsigned int ie=0;ie<idEles.size();++ie)
            {
                int indexElectron = idEles[ie];
                float dRTauElectron = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau],analysisTree.electron_eta[indexElectron],analysisTree.electron_phi[indexElectron]);
                if (dRTauElectron<0.5) foundNearElectron = true;
            }
            //checking whether there are other muons overlapped with the tau
            for(unsigned int im=0;im<idMuons.size();++im)
            {
                int indexMuon = idMuons[im];
                float dRTauMuon = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau],analysisTree.muon_eta[indexMuon],analysisTree.muon_phi[indexMuon]);
                if (dRTauMuon<0.5) foundNearMuon = true;
            }
            //checking whether there are other taus overlapped with the tau
            for(unsigned int it2=0;it2<idTaus.size();++it2)
            {
                int indexTau2 = idTaus[it2];
                if (indexTau == indexTau2) continue;
                float dRTaus = deltaR(analysisTree.tau_eta[indexTau],analysisTree.tau_phi[indexTau],analysisTree.tau_eta[indexTau2],analysisTree.tau_phi[indexTau2]);
                if (dRTaus<0.5) foundNearTau = true;
            }
            if ((!foundNearElectron)&&(!foundNearMuon)&&(!foundNearTau))
                goodTaus.push_back(indexTau);
        }
        //end checking seprations for leptons ------------
        
        //Final state selection: EE
        if(goodEles.size() < 2) continue;
        //FInal state selection: Higgs candidates
        if(HiggsFinalState == "EE")
        {
            if(goodEles.size() != 4) continue;
            if(goodMuons.size() != 0) continue;
        }
        if(HiggsFinalState == "EM")
        {
            if(goodEles.size() != 3) continue;
            if(goodMuons.size() != 1) continue;
        }
        if(HiggsFinalState == "ET")
        {
            if(goodEles.size() != 3) continue;
            if(goodMuons.size() != 0) continue;
            if(goodTaus.size() < 1) continue;
        }
        if(HiggsFinalState == "MT")
        {
            if(goodEles.size() != 2) continue;
            if(goodMuons.size() != 1) continue;
            if(goodTaus.size() < 1) continue;
        }
        if(HiggsFinalState == "TT")
        {
            if(goodEles.size() != 2) continue;
            if(goodMuons.size() != 0) continue;
            if(goodTaus.size() <= 2) continue;
        }
        
        //std::cout << "good electrons:" << goodEles.size() << std::endl;
        //std::cout << "good muons:" << goodMuons.size() << std::endl;
        //std::cout << "good tau:" << goodTaus.size() << std::endl;
        
        //under construction
        //Good Z candidates selection
        int indexElectron1 = -1;
        int indexElectron2 = -1;
        bool foundZOS = false;
        TLorentzVector Ele_1, Ele_2;
        float m_EE = 0;
        float m_Dist = 9999;
        for(unsigned int ie1=0;ie1<goodEles.size();++ie1)
        {
            int indexTemp1 = goodEles[ie1];
            float q1 = analysisTree.electron_charge[indexTemp1];
            for(unsigned int ie2=ie1+1;ie2<(goodEles.size()-ie1);++ie2)
            {
                int indexTemp2 = goodEles[ie2];
                if (indexTemp1 == indexTemp2) continue;
                float q2 = analysisTree.electron_charge[indexTemp2];
                if (q1*q2 > 0) continue;
                foundZOS = true;
                Ele_1.SetXYZM(analysisTree.electron_px[indexTemp1],analysisTree.electron_py[indexTemp1],analysisTree.electron_pz[indexTemp1],electronMass);
                Ele_2.SetXYZM(analysisTree.electron_px[indexTemp2],analysisTree.electron_py[indexTemp2],analysisTree.electron_pz[indexTemp2],electronMass);
                if (fabs((Ele_1+Ele_2).M() - 91.1876) < m_Dist)
                {
                    m_EE = (Ele_1+Ele_2).M();
                    m_Dist = fabs(m_EE - 91.1876);
                    if (analysisTree.electron_pt[indexTemp1]>analysisTree.electron_pt[indexTemp2])
                    {
                        indexElectron1 = indexTemp1;
                        indexElectron2 = indexTemp2;
                    }
                    else
                    {
                        indexElectron1 = indexTemp2;
                        indexElectron2 = indexTemp1;
                    }
                }
            }
        }
        
        if(foundZOS==false) continue;
        
        //std::cout << "index electron 1: " << indexElectron1 << std::endl;
        //std::cout << "index electron 2: " << indexElectron2 << std::endl;
        //std::cout << "distance between true Z mass and reco mass: " << m_Dist << std::endl;
        //std::cout << "reco Z mass " << m_EE << std::endl;
        //std::cout << "__________________________________________________________ " << std::endl;

        //Higgs boson candidates selection
        int indexLeg3 = -1;
        int indexLeg4 = -1;
        float Higgs_LT = -9999;
        if(HiggsFinalState =="EM")
        for(unsigned int ie=0;ie<goodEles.size();++ie)
        {
                int indexTemp1 = goodEles[ie];
                if ((indexTemp1 == indexElectron1)||(indexTemp1 == indexElectron2))
                    continue;
                for(unsigned int im=0;im<goodMuons.size();++im)
                {
                    int indexTemp2 = goodMuons[im];
                    if (analysisTree.electron_pt[indexTemp1]+analysisTree.muon_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + analysisTree.muon_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
        }
        if(HiggsFinalState =="ET")
        {
            for(unsigned int ie=0;ie<goodEles.size();++ie)
            {
                int indexTemp1 = goodEles[ie];
                if ((indexTemp1 == indexElectron1)||(indexTemp1 == indexElectron2))
                    continue;
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    if (analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTemp2]<0.5) continue;
                    //(no needed for Synchronization)
                    //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTemp2]<0.5) continue;
                    if (analysisTree.electron_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }
        
        if(HiggsFinalState =="MT")
        {
            for(unsigned int im=0;im<goodMuons.size();++im)
            {
                int indexTemp1 = goodMuons[im];
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    if (analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTemp2]<0.5) continue;
                    //(no needed for Synchronization)
                    //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTemp2]<0.5) continue;
                    if (analysisTree.muon_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.muon_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }
        
        if(HiggsFinalState =="TT")
        {
            for(unsigned int it=0;it<goodTaus.size();++it)
            {
                int indexTemp1 = goodTaus[it];
                float q1 = analysisTree.tau_charge[indexTemp1];
                for(unsigned int it2=it+1;it2<(goodTaus.size()-it);++it2)
                {
                    int indexTemp2 = goodTaus[it2];
                    //if (analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTemp2]<0.5) continue;
                    //(no needed for Synchronization)
                    //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTemp2]<0.5) continue;
                    if (analysisTree.tau_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.tau_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }
        //(No need for Synchronization)
        //if(Higgs_LT<=40) continue;
        //if(HiggsFinalState =="TT")
        //{
        //    if(Higgs_LT < 60) continue;
        //}
                
        //std::cout << "index muon 1: " << indexLeg3 << std::endl;
        //std::cout << "index tau 1: " << indexLeg4 << std::endl;
        //std::cout << "Higgs LT (pT_muon + pT_tau): " << Higgs_LT << std::endl;
        //std::cout << "__________________________________________________________ " << std::endl;
        //Filling up lepton kinematics
        pt_1 = analysisTree.electron_pt[indexElectron1];
        phi_1 = analysisTree.electron_phi[indexElectron1];
        eta_1 = analysisTree.electron_eta[indexElectron1];
        
        float absIso_1 = analysisTree.electron_chargedHadIso[indexElectron1];
        float neutralIso_1 = analysisTree.electron_neutralHadIso[indexElectron1] + analysisTree.electron_photonIso[indexElectron1] - 0.5*analysisTree.electron_puIso[indexElectron1];
        neutralIso_1 = TMath::Max(float(0),neutralIso_1);
        absIso_1 += neutralIso_1;
        iso_1 = absIso_1/analysisTree.electron_pt[indexElectron1];
        
        IdRawMva_1 = analysisTree.electron_mva_value_noIso_Fall17_v2[indexElectron1];
        
        pt_2 = analysisTree.electron_pt[indexElectron2];
        phi_2 = analysisTree.electron_phi[indexElectron2];
        eta_2 = analysisTree.electron_eta[indexElectron2];
        
        float absIso_2 = analysisTree.electron_chargedHadIso[indexElectron2];
        float neutralIso_2 = analysisTree.electron_neutralHadIso[indexElectron2] + analysisTree.electron_photonIso[indexElectron2] - 0.5*analysisTree.electron_puIso[indexElectron2];
        neutralIso_2 = TMath::Max(float(0),neutralIso_2);
        absIso_2 += neutralIso_2;
        iso_2 = absIso_2/analysisTree.electron_pt[indexElectron2];
        IdRawMva_2 = analysisTree.electron_mva_value_noIso_Fall17_v2[indexElectron2];
        
        Z_DR = deltaR(analysisTree.electron_eta[indexElectron1],analysisTree.electron_phi[indexElectron1],analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2]);
        Z_SS = false;
        if (analysisTree.electron_charge[indexElectron1]*analysisTree.electron_charge[indexElectron2]>0)
            Z_SS = true;
        
        TLorentzVector Leg3;
        TLorentzVector Leg4;
        if(HiggsFinalState == "EM")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);

            pt_4 = analysisTree.muon_pt[indexLeg4];
            phi_4 = analysisTree.muon_phi[indexLeg4];
            eta_4 = analysisTree.muon_eta[indexLeg4];
            float absIso_4 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg4];
            float neutralIso_4 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg4] + analysisTree.muon_r04_sumPhotonEt[indexLeg4] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg4];
            neutralIso_4 = TMath::Max(float(0),neutralIso_4);
            absIso_4 += neutralIso_4;
            iso_4 = absIso_4/analysisTree.muon_pt[indexLeg4];
            Leg4.SetXYZM(analysisTree.muon_px[indexLeg4],analysisTree.muon_py[indexLeg4],analysisTree.muon_pz[indexLeg4],muonMass);
        }
        if(HiggsFinalState == "ET")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            
            float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);

            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);
            H_SS = false;
            if (analysisTree.electron_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
        }
        if(HiggsFinalState == "MT")
        {
            pt_3 = analysisTree.muon_pt[indexLeg3];
            phi_3 = analysisTree.muon_phi[indexLeg3];
            eta_3 = analysisTree.muon_eta[indexLeg3];
            
            float absIso_3 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg3] + analysisTree.muon_r04_sumPhotonEt[indexLeg3] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.muon_pt[indexLeg3];
            Leg3.SetXYZM(analysisTree.muon_px[indexLeg3],analysisTree.muon_py[indexLeg3],analysisTree.muon_pz[indexLeg3],muonMass);
        
            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);

            H_SS = false;
            if (analysisTree.muon_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
            //std::cout << "H_SS:  " << H_SS << std::endl;
        }
        if(HiggsFinalState == "TT")
        {
            pt_3 = analysisTree.tau_pt[indexLeg3];
            phi_3 = analysisTree.tau_phi[indexLeg3];
            eta_3 = analysisTree.tau_eta[indexLeg3];
            iso_3 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg3];
            Leg3.SetXYZM(analysisTree.tau_px[indexLeg3],analysisTree.tau_py[indexLeg3],analysisTree.tau_pz[indexLeg3],analysisTree.tau_mass[indexLeg3]);

            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);

            H_SS = false;
            if (analysisTree.tau_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
            //std::cout << "H_SS:  " << H_SS << std::endl;
        }
    
        //new MC matching piece of code
        gen_match_1 = 6;
        gen_match_2 = 6;
        gen_match_3 = 6;
        gen_match_4 = 6;
        float minDR_1 = 0.2;
        float minDR_2 = 0.2;
        float minDR_3 = 0.2;
        float minDR_4 = 0.2;
        if (!isData)
        {
            for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen)
            {
                float ptGen = PtoPt(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen]);
                bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
                bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
                bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
                bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
                bool isAnyType = type1 || type2 || type3 || type4;
                if (isAnyType)
                {
                    float etaGen = PtoEta(analysisTree.genparticles_px[igen],analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen]);
                    float phiGen = PtoPhi(analysisTree.genparticles_px[igen],analysisTree.genparticles_py[igen]);
                    float deltaR_1 = deltaR(analysisTree.electron_eta[indexElectron1],analysisTree.electron_phi[indexElectron1],etaGen,phiGen);
                    if (deltaR_1<minDR_1)
                    {
                        minDR_1 = deltaR_1;
                        if (type1) gen_match_1 = 1;
                        else if (type2) gen_match_1 = 2;
                        else if (type3) gen_match_1 = 3;
                        else if (type4) gen_match_1 = 4;
                    }
                    float deltaR_2 = deltaR(analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2],etaGen,phiGen);
                    if (deltaR_2<minDR_2)
                    {
                        minDR_2 = deltaR_2;
                        if (type1) gen_match_2 = 1;
                        else if (type2) gen_match_2 = 2;
                        else if (type3) gen_match_2 = 3;
                        else if (type4) gen_match_2 = 4;
                    }

                    float deltaR_3 = deltaR(Leg3.Eta(),Leg3.Phi(),etaGen,phiGen);
                    if (deltaR_3<minDR_3)
                    {
                        minDR_3 = deltaR_3;
                        if (type1) gen_match_3 = 1;
                        else if (type2) gen_match_3 = 2;
                        else if (type3) gen_match_3 = 3;
                        else if (type4) gen_match_3 = 4;
                    }
                    float deltaR_4 = deltaR(Leg4.Eta(),Leg4.Phi(),etaGen,phiGen);
                    if (deltaR_4<minDR_4)
                    {
                        minDR_4 = deltaR_4;
                        if (type1) gen_match_4 = 1;
                        else if (type2) gen_match_4 = 2;
                        else if (type3) gen_match_4 = 3;
                        else if (type4) gen_match_4 = 4;
                    }
                }
            }
            for (unsigned int igen=0; igen < analysisTree.gentau_count; ++igen)
            {
                if (analysisTree.gentau_visibleNoLep_pt[igen]>15.)
                {
                    float deltaR_1 = deltaR(analysisTree.electron_eta[indexElectron1],analysisTree.electron_phi[indexElectron1],analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_1<minDR_1)
                    {
                        minDR_1 = deltaR_1;
                        gen_match_1 = 5;
                    }
                    float deltaR_2 = deltaR(analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2],analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_2<minDR_2)
                    {
                        minDR_2 = deltaR_2;
                        gen_match_2 = 5;
                    }

                    float deltaR_3 = deltaR(Leg3.Eta(),Leg3.Phi(),analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_3<minDR_3)
                    {
                        minDR_3 = deltaR_3;
                        gen_match_3 = 5;
                    }
                    float deltaR_4 = deltaR(Leg4.Eta(),Leg4.Phi(),analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_4<minDR_4)
                    {
                        minDR_4 = deltaR_4;
                        gen_match_4 = 5;
                    }
                }
            }
        }

        //uncorrected met
        met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
        metphi = analysisTree.pfmetcorr_phi;
        metcov00 = analysisTree.pfmetcorr_sigxx;
        metcov01 = analysisTree.pfmetcorr_sigxy;
        metcov10 = analysisTree.pfmetcorr_sigyx;
        metcov11 = analysisTree.pfmetcorr_sigyy;

        TLorentzVector metLV;
        metLV.SetXYZT(met*TMath::Cos(metphi),met*TMath::Sin(metphi),0,TMath::Sqrt(met*TMath::Sin(metphi)*met*TMath::Sin(metphi) + met*TMath::Cos(metphi)*met*TMath::Cos(metphi)));
        
        TLorentzVector electron1;
        electron1.SetXYZM(analysisTree.electron_px[indexElectron1],analysisTree.electron_py[indexElectron1],analysisTree.electron_pz[indexElectron1],electronMass);
        pfmt_1=mT(electron1,metLV);
        
        TLorentzVector electron2;
        electron2.SetXYZM(analysisTree.electron_px[indexElectron2],analysisTree.electron_py[indexElectron2],analysisTree.electron_pz[indexElectron2],electronMass);
        pfmt_2=mT(electron2,metLV);
        
        pfmt_3=mT(Leg3,metLV);
        pfmt_4=mT(Leg4,metLV);

        m_vis = (Leg3+Leg4).M();
        pt_tt = (Leg3+Leg4).Pt();

        mll = (electron1+electron2).M();
        Z_Pt = (electron1+electron2).Pt();
        
        if ( mll < 60 || mll > 120)
            continue;
        
        Mass = (electron1+electron2+Leg3+Leg4).M();
        H_DR = deltaR(Leg3.Eta(),Leg3.Phi(),Leg4.Eta(),Leg4.Phi());
        DR_12 = Z_DR;
        DR_13 = deltaR(analysisTree.electron_eta[indexElectron1],analysisTree.electron_phi[indexElectron1],Leg3.Eta(),Leg3.Phi());
        DR_14 = deltaR(analysisTree.electron_eta[indexElectron1],analysisTree.electron_phi[indexElectron1],Leg4.Eta(),Leg4.Phi());
        DR_23 = deltaR(analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2],Leg3.Eta(),Leg3.Phi());
        DR_24 = deltaR(analysisTree.electron_eta[indexElectron2],analysisTree.electron_phi[indexElectron2],Leg4.Eta(),Leg4.Phi());
        DR_34 = H_DR;

        // Jets selection
        vector<unsigned int> jets; jets.clear();
        vector<unsigned int> jetspt20; jetspt20.clear();
        vector<unsigned int> bjets; bjets.clear();
        vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
        vector<unsigned int> bjetsRaw; bjetsRaw.clear();
            
        int indexLeadingJet = -1;
        float ptLeadingJet = -1;
            
        int indexSubLeadingJet = -1;
        float ptSubLeadingJet = -1;
            
        int indexLeadingBJet = -1;
        float ptLeadingBJet = -1;
        
        int indexSubLeadingBJet = -1;
        float ptSubLeadingBJet = -1;
        for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
        {
            float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
            float jetEta = analysisTree.pfjet_eta[jet];
            if (absJetEta>jetEtaCut) continue;
                
            float jetPt = analysisTree.pfjet_pt[jet];
            //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
            if (jetPt<jetPtLowCut) continue;
                
            bool isPFJetId = looseJetiD(analysisTree,int(jet));
            if (!isPFJetId) continue;
            // non-cleaned jet at the moment
            
            // jetId
            if (jetPt>jetPtLowCut)
                jetspt20.push_back(jet);
            
            if (absJetEta<bJetEtaCut)
            {   // jet within b-tagging acceptance
                    
                bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                //missing b-tag SFs

                if (tagged)
                {
                    bjets.push_back(jet);
                    if (jetPt<ptLeadingBJet&&jetPt>ptSubLeadingBJet)
                    {
                        indexSubLeadingBJet = jet;
                        ptSubLeadingBJet = jetPt;
                    }
                    if (jetPt>ptLeadingBJet)
                    {
                        ptLeadingBJet = jetPt;
                        indexLeadingBJet = jet;
                    }
                }
            }

            if (jetPt>jetPtHighCut)
                jets.push_back(jet);
                
            if (indexLeadingJet>=0)
            {
                if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet)
                {
                    indexSubLeadingJet = jet;
                    ptSubLeadingJet = jetPt;
                }
            }
                
            if (jetPt>ptLeadingJet)
            {
                indexSubLeadingJet = indexLeadingJet;
                ptSubLeadingJet = ptLeadingJet;
                indexLeadingJet = jet;
                ptLeadingJet = jetPt;
            }
        }//end jet selection
        
        njets = jets.size();
        
        int njetsMax = njets;

        njetspt20 = jetspt20.size();
        nbtag = bjets.size();

        //missing b-tag weight
        
        bpt_1 = -9999;
        beta_1 = -9999;
        bphi_1 = -9999;
        bhadronFlavor_1 = -1;
        
        bpt_2 = -9999;
        beta_2 = -9999;
        bphi_2 = -9999;
        bhadronFlavor_2 = -1;

        if (indexLeadingBJet>=0)
        {
            bpt_1 = analysisTree.pfjet_pt[indexLeadingBJet];
            beta_1 = analysisTree.pfjet_eta[indexLeadingBJet];
            bphi_1 = analysisTree.pfjet_phi[indexLeadingBJet];
        }
        if(!isData)
            bhadronFlavor_1 = abs(analysisTree.pfjet_flavour[indexLeadingBJet]);

        if (indexSubLeadingBJet>=0)
        {
            bpt_2 = analysisTree.pfjet_pt[indexSubLeadingBJet];
            beta_2 = analysisTree.pfjet_eta[indexSubLeadingBJet];
            bphi_2 = analysisTree.pfjet_phi[indexSubLeadingBJet];
        }
        if(!isData)
            bhadronFlavor_2 = abs(analysisTree.pfjet_flavour[indexSubLeadingJet]);

        jpt_1 = -9999;
        jeta_1 = -9999;
        jphi_1 = -9999;
            
        if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
            cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
            
        if (indexLeadingJet>=0)
        {
            jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
            jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
            jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
        //        cout << "Leading jet pt = " << jpt_1 << "   eta = " << jeta_1 << endl;
        //        for (auto const& Name : uncertNames) {
        //          cout << "    " << Name << " : " << jecUncertainties->getUncertainty(Name,jpt_1,jeta_1) << endl;
        //        }
        }
            
        jpt_2 = -9999;
        jeta_2 = -9999;
        jphi_2 = -9999;
            
        if (indexSubLeadingJet>=0)
        {
            jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
            jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
            jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
        }
        
        SynTree->Fill();
        selEventsIsoMuons++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                     = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  //cout << "weight used:" << weight << std::endl;

  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs   :";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  
  file->Write();
  file->Close();
  delete file;
  
  }
