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
#include "TRFIOFile.h"
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
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

    const bool isData = cfg.get<bool>("IsData");
    const bool isData2016BCDEF = cfg.get<bool>("IsData2016BCDEF");
    const bool isData2016GH = cfg.get<bool>("IsData2016GH");
    const bool isDY = cfg.get<bool>("IsDY");
    const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

    const string jsonFile = cfg.get<string>("jsonFile");
    const string dataPUFile = cfg.get<string>("DataPUFile");
    const string mcPUFile = cfg.get<string>("MCPUFile");

    // pile up reweighting
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");

    // kinematic cuts on muons
    const float ptMuonTagCut  = cfg.get<float>("ptMuonTagCut");
    const float etaMuonTagCut = cfg.get<float>("etaMuonTagCut");
    const float dxyMuonTagCut     = cfg.get<float>("dxyMuonTagCut");
    const float dzMuonTagCut      = cfg.get<float>("dzMuonTagCut");
    const float isoMuonTagCut = cfg.get<float>("isoMuonTagCut");

    //kinematic cuts on taus
    const float ptTauCut = cfg.get<float>("ptTauCut");
    const float etaTauCut = cfg.get<float>("etaTauCut");
    const float dzTauCut = cfg.get<float>("dzTauCut");
    
    // topological cuts
    const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
    const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
    const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
    const bool isoDR03         = cfg.get<bool>("IsoDR03");

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
    
    //Tag Muon energy scale
    const float tagmuonScale = cfg.get<float>("TagMuonScale");
    
    //Probe Muon energy scale
    const float probemuonScale = cfg.get<float>("ProbeMuonScale");
    
    //Probe Tau energy scale
    const float probetauScale = cfg.get<float>("ProbeTauScale");
    
    //visble resolution scale
    const float resoScale = cfg.get<float>("ResoScale");

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

  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
    
    //declare the tree
    TTree * muonTree = new TTree("MuTauFR","MuTauFR");

    Bool_t isZLL;
    Bool_t isZMM;
    Bool_t isZEE;
    Bool_t isZTT;
    
    TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
    TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
    TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
    TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
    TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
    
    Float_t         mcweight;
    Float_t         puweight;
    Float_t         trigweight_1;
    Float_t         idweight_1;
    Float_t         isoweight_1;
    Float_t         effweight;
    Float_t         iso_1;
    Int_t           TPmatching_status;
    
    Float_t m_vis;
    Float_t m_vis_gen;
    Float_t m_vis_reso_scale_up;
    Float_t m_vis_reso_scale_down;
    
    //Float_t m_vis_tagmuon_scale_up;
    //Float_t m_vis_tagmuon_scale_down;
    
    //Float_t m_vis_probemuon_scale_up;
    //Float_t m_vis_probemuon_scale_down;
    
    
    Float_t m_vis_probetau_scale_up;
    Float_t m_vis_probetau_scale_down;
    
    Float_t mt_1;
    Bool_t os;
    Bool_t tauagainstMuonLoose;
    Bool_t tauagainstMuonTight;
    Bool_t taubyLooseCombinedIsolationDeltaBetaCorr3Hits;
    
    Float_t PtTag;
    Float_t EtaTag;
    Float_t PtProbe;
    Float_t EtaProbe;
    Float_t met;
    
    muonTree->Branch("isZLL",&isZLL,"isZLL/O");
    muonTree->Branch("isZEE",&isZEE,"isZEE/O");
    muonTree->Branch("isZMM",&isZMM,"isZMM/O");
    muonTree->Branch("isZTT",&isZTT,"isZTT/O");
    
    muonTree->Branch("mcweight",&mcweight,"mcweight/F");
    muonTree->Branch("puweight", &puweight, "puweight/F");
    muonTree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
    muonTree->Branch("idweight_1", &idweight_1, "idweight_1/F");
    muonTree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
    muonTree->Branch("effweight", &effweight, "effweight/F");
    muonTree->Branch("iso_1",&iso_1,"iso_1/F");
    muonTree->Branch("TPmatching_status",&TPmatching_status,"TPmatching_status/I");
    
    muonTree->Branch("m_vis",&m_vis,"m_vis/F");
    muonTree->Branch("m_vis_gen",&m_vis_gen,"m_vis_gen/F");
    muonTree->Branch("m_vis_reso_scale_up",&m_vis_reso_scale_up,"m_vis_reso_scale_up/F");
    muonTree->Branch("m_vis_reso_scale_down",&m_vis_reso_scale_down,"m_vis_reso_scale_down/F");
    
    //muonTree->Branch("m_vis_tagmuon_scale_up",&m_vis_tagmuon_scale_up,"m_vis_tagmuon_scale_up/F");
    //muonTree->Branch("m_vis_tagmuon_scale_down",&m_vis_tagmuon_scale_down,"m_vis_tagmuon_scale_down/F");
    //muonTree->Branch("m_vis_probemuon_scale_up",&m_vis_probemuon_scale_up,"m_vis_probemuon_scale_up/F");
    //muonTree->Branch("m_vis_probemuon_scale_down",&m_vis_probemuon_scale_down,"m_vis_probemuon_scale_down/F");
    muonTree->Branch("m_vis_probetau_scale_up",&m_vis_probetau_scale_up,"m_vis_probetau_scale_up/F");
    muonTree->Branch("m_vis_probetau_scale_down",&m_vis_probetau_scale_down,"m_vis_probetau_scale_down/F");
    
    muonTree->Branch("mt_1",&mt_1,"mt_1/F");
    muonTree->Branch("os",&os,"os/O");
    muonTree->Branch("tauagainstMuonLoose",&tauagainstMuonLoose,"tauagainstMuonLoose/O");
    muonTree->Branch("tauagainstMuonTight",&tauagainstMuonTight,"tauagainstMuonTight/O");
    muonTree->Branch("taubyLooseCombinedIsolationDeltaBetaCorr3Hits",&taubyLooseCombinedIsolationDeltaBetaCorr3Hits,"taubyLooseCombinedIsolationDeltaBetaCorr3Hits/O");
    muonTree->Branch("PtTag",&PtTag,"PtTag/F");
    muonTree->Branch("EtaTag",&EtaTag,"EtaTag/F");
    muonTree->Branch("PtProbe",&PtProbe,"PtProbe/F");
    muonTree->Branch("EtaProbe",&EtaProbe,"EtaProbe/F");
    muonTree->Branch("met",&met,"met/F");


  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
  TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

  // reweighting official recipe 
  // initialize pile up object
  PileUp * PUofficial = new PileUp();
  
  if (applyPUreweighting) {
    //Temprory DY MC Summer16, other MCs are Spring16, need to switch mannually, won't need this when all the Summer16 samples come
    //TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_runBCDEFGH_Rereco_xsec69p2mb_bin50.root","read");
    //TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring16_PU.root", "read");

    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(dataPUFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(mcPUFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUofficial->set_h_data(PU_data); 
    PUofficial->set_h_MC(PU_mc);
  }
    // Muon scale factors
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerFile));

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
        
        
        isZLL = false;
        isZEE = false;
        isZMM = false;
        isZTT = false;
        
        // weights
        puweight = 1;
        trigweight_1 = 1;
        idweight_1 = 1;
        isoweight_1 = 1;
        effweight = 1;
        mcweight =1 ;

      //------------------------------------------------

        if (!isData)
          weight *=analysisTree.genweight;
        
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
      
        // vertex cuts
        if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
        if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
        float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);
        if (dVertex>dVertexCut) continue;
      
      // muon selection

      vector<unsigned int> allMuons; allMuons.clear();
      vector<unsigned int> idMuons; idMuons.clear();
      vector<unsigned int> isoMuons; isoMuons.clear();
      vector<float> isoMuonsValue; isoMuonsValue.clear();
      vector<float> allMuonsIso; allMuonsIso.clear();
      vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
      vector<bool> isMuonMatchedSingleMuFilter; isMuonMatchedSingleMuFilter.clear();
        
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	bool muPassed    = true;
	bool muSingleMatched = false;
	if (analysisTree.muon_pt[im]<ptMuonTagCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonTagCut) continue;
	allMuons.push_back(im);
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonTagCut) muPassed = false;
	if (fabs(analysisTree.muon_dz[im])>dzMuonTagCut) muPassed = false;
    // Muon POG suggestions for Muon ID for Moriond17
    if(isData && isData2016BCDEF)
    {
        if (!analysisTree.muon_isICHEP[im]) muPassed = false;
    }
    if(isData && isData2016GH)
    {
        if (!analysisTree.muon_isMedium[im]) muPassed = false;
    }
    if(!isData)
    {
        if (!analysisTree.muon_isMedium[im]) muPassed = false;
    }
    if (muPassed) idMuons.push_back(im);
	float absIso = 0;
	if (isoDR03) { 
        absIso = analysisTree.muon_r03_sumChargedHadronPt[im];
	  float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] +
	    analysisTree.muon_r03_sumPhotonEt[im] -
	    0.5*analysisTree.muon_r03_sumPUPt[im];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  absIso += neutralIso;
	}
	else {
	  absIso = analysisTree.muon_chargedHadIso[im];
          float neutralIso = analysisTree.muon_neutralHadIso[im] +
            analysisTree.muon_photonIso[im] -
            0.5*analysisTree.muon_puIso[im];
          neutralIso = TMath::Max(float(0),neutralIso);
          absIso += neutralIso;
	}
	float relIso = absIso/analysisTree.muon_pt[im];
	allMuonsIso.push_back(relIso);
	if (muPassed && relIso<isoMuonTagCut) {
	  isoMuons.push_back(im);
	  isoMuonsValue.push_back(relIso);
	}
	if (relIso>isoMuonTagCut) muPassed = false;
	//	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
	isMuonPassedIdIso.push_back(muPassed);
      }
        
        //tau selection
        vector<unsigned int> goodTaus; goodTaus.clear();
        for(unsigned int it=0;it<analysisTree.tau_count;++it)
        {
            if (analysisTree.tau_decayModeFinding[it]<=0.5) continue;
            if (analysisTree.tau_pt[it]<ptTauCut) continue;
            if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
            //if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
            if (fabs(analysisTree.tau_charge[it])<0.5||fabs(analysisTree.tau_charge[it])>1.5) continue;
            if (fabs(analysisTree.primvertex_z-analysisTree.tau_vertexz[it])>0.05) continue;
            //if (analysisTree.tau_againstElectronVLooseMVA6[it] < 0.5) continue;
            goodTaus.push_back(it);
        }
        
        
        vector<unsigned int> genMuons; genMuons.clear();
        if(!isData)
        {
            for(unsigned int igenMuon=0;igenMuon<analysisTree.genparticles_count;++igenMuon)
            {
                //cout <<"pdgid:"<<analysisTree.genparticles_pdgid[igenMuon]<< endl;
                if (fabs(analysisTree.genparticles_pdgid[igenMuon])==13 && analysisTree.genparticles_status[igenMuon]==1 &&analysisTree.genparticles_fromHardProcess[igenMuon])
                {
                    genMuons.push_back(igenMuon);
                }
            }
        }
        //cout<<"size of gen muons:"<< genMuons.size()<<endl;
        
        std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
        std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
        std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
        std::vector<TLorentzVector> promptMuons; promptMuons.clear();
        
        if (!isData)
        {
            for (unsigned int igentau=0; igentau < analysisTree.gentau_count; ++igentau) {
                TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
                                                    analysisTree.gentau_py[igentau],
                                                    analysisTree.gentau_pz[igentau],
                                                    analysisTree.gentau_e[igentau]);
                TLorentzVector tauVisLV; tauVisLV.SetXYZT(analysisTree.gentau_visible_px[igentau],
                                                          analysisTree.gentau_visible_py[igentau],
                                                          analysisTree.gentau_visible_pz[igentau],
                                                          analysisTree.gentau_visible_e[igentau]);
                if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
                    promptTausFirstCopy.push_back(tauLV);
                }
                if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {
                    promptTausLastCopy.push_back(tauVisLV);
                }
            }
            
            for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
                
                TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                    analysisTree.genparticles_py[igen],
                                                    analysisTree.genparticles_pz[igen],
                                                    analysisTree.genparticles_e[igen]);
                
                if (fabs(analysisTree.genparticles_pdgid[igen])==11) {
                    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                        promptElectrons.push_back(genLV);
                    }
                }
                
                if (fabs(analysisTree.genparticles_pdgid[igen])==13) {
                    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                        promptMuons.push_back(genLV);
                    }
                }
            }
            
            if (isDY) {
                if (promptTausFirstCopy.size()==2) {
                    isZTT = true; isZMM = false; isZEE = false;
                }
                else if (promptMuons.size()==2) {
                    isZTT = false; isZMM = true; isZEE = false;
                }
                else {
                    isZTT = false; isZMM = false; isZEE = true;
                }
            }
            
            ALL->Fill(0.0);
            if (isZMM) ZMM->Fill(0.);
            if (isZEE) ZEE->Fill(0.);
            if (isZTT) ZTT->Fill(0.);
            isZLL = isZMM || isZEE;
        }
        //cout<< "size of gen muons:"<<genMuons.size()<<endl;
        
        ////////////////////////////////////////////////////////////////////////////////////////
        bool isIsoMuonsPair = false;
        float isoMin = 9999;
        if (isoMuons.size()>0) {
            for (unsigned int im1=0; im1<isoMuons.size(); ++im1) {
                unsigned int index1 = isoMuons[im1];
                bool isMu1matched = false;
                //cout << "ICHEP ID:"<< analysisTree.muon_isICHEP[index1]<< endl;
                if(isData)
                {
                    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
                    {
                        float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                        if (dRtrig>DRTrigMatch) continue;
                        if (analysisTree.trigobject_filters[iT][nMuonFilter] && analysisTree.muon_pt[index1] > ptMuonTagCut && fabs(analysisTree.muon_eta[index1]) < etaMuonTagCut)
                        isMu1matched = true;
                    }
                }
                bool genRecoTagMuonMatched = false;
                TLorentzVector genMuon1;
                float indexGen1 = -9999;
                if(!isData)
                {
                    for(unsigned int iMuonGen=0; iMuonGen<genMuons.size();++iMuonGen)
                    {
                        unsigned int indexGenMuon = genMuons[iMuonGen];
                        float genMuonEta = PtoEta(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon],analysisTree.genparticles_pz[indexGenMuon]);
                        float genMuonPhi = PtoPhi(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon]);
                        float dRRecoGenTag = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],genMuonEta,genMuonPhi);
                        if(dRRecoGenTag < 0.1)
                        {
                            genRecoTagMuonMatched = true;
                            genMuon1.SetXYZM(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon],analysisTree.genparticles_pz[indexGenMuon],muonMass);
                            //cout << "gen muon 1 px:" << analysisTree.genparticles_px[indexGenMuon]<< endl;
                            indexGen1 = indexGenMuon;
                        }
                    }
                }
                                
                if (isMu1matched || (!isData)) {
                    for (unsigned int iTau=0; iTau<goodTaus.size(); ++iTau) {
                        unsigned int indexProbe = goodTaus[iTau];
                        float q1 = analysisTree.muon_charge[index1];
                        float q2 = analysisTree.tau_charge[indexProbe];
                        float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe]);
                        if (dR<dRleptonsCut) continue;
                        
                        bool genRecoProbeMuonMatched = false;
                        TLorentzVector genMuon2;

                        for(unsigned int iMuonGen=0; iMuonGen<genMuons.size();++iMuonGen)
                        {
                            unsigned int indexGenMuon = genMuons[iMuonGen];
                            if(indexGenMuon!=indexGen1)
                            {
                                float genMuonEta = PtoEta(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon],analysisTree.genparticles_pz[indexGenMuon]);
                                float genMuonPhi = PtoPhi(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon]);
                                float dRRecoGenProbe = deltaR(analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe],genMuonEta,genMuonPhi);
                                if(dRRecoGenProbe<0.1)
                                {
                                    genRecoProbeMuonMatched = true;
                                    genMuon2.SetXYZM(analysisTree.genparticles_px[indexGenMuon],analysisTree.genparticles_py[indexGenMuon],analysisTree.genparticles_pz[indexGenMuon],muonMass);
                                    //cout << "gen muon 2 px:" << analysisTree.genparticles_px[indexGenMuon]<< endl;
                                }
                            }
                        }
                        
                        TLorentzVector genTau2;
                        bool genRecoProbeTauMatched = false;
                        if(isZTT)//to match with GEN Tau object, building template for tau energy shapes systematics
                        {
                            for(unsigned int igentau = 0; igentau < analysisTree.gentau_count; igentau++)
                            {
                                float genTauEta = PtoEta(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau],analysisTree.gentau_pz[igentau]);
                                float genTauPhi = PtoPhi(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau]);
                                float dRRecoGenTauProbe = deltaR(analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe],genTauEta,genTauPhi);
                                if(dRRecoGenTauProbe < 0.1)
                                {
                                    genTau2.SetXYZM(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau],analysisTree.gentau_pz[igentau],analysisTree.tau_mass[igentau]);
                                    genRecoProbeTauMatched = true;
                                }
                            }
                            
                        }
                        
                        //setting tag and probe matching status
                        TPmatching_status = 3;
                        if(genRecoTagMuonMatched && !genRecoProbeMuonMatched && !genRecoProbeTauMatched)
                            TPmatching_status = 0;
                        if(genRecoTagMuonMatched && genRecoProbeMuonMatched)
                            TPmatching_status = 1;//MC truth
                        if(genRecoProbeTauMatched)
                            TPmatching_status = 2;
                        
                        
                        //cout << "TP mathing status:" << TPmatching_status << endl;
                        // choosing mva met
                        unsigned int metMuTau = 0;
                        bool mvaMetFound = false;
                        for (unsigned int iMet=0; iMet<analysisTree.mvamet_count; ++iMet)
                        {
                            //cout << iMet << endl;
                            if (analysisTree.mvamet_channel[iMet]==3)
                            {
                                if (analysisTree.mvamet_lep1[iMet]==indexProbe &&analysisTree.mvamet_lep2[iMet]==index1)
                                {
                                    metMuTau = iMet;
                                    mvaMetFound = true;
                                }
                            }
                        }
                        float mvamet_ex = analysisTree.mvamet_ex[metMuTau];
                        float mvamet_ey = analysisTree.mvamet_ey[metMuTau];
                        float mvamet = TMath::Sqrt(mvamet_ex*mvamet_ex+mvamet_ey*mvamet_ey);
                        
                        met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
                        float dPhiMETMuon = dPhiFrom2P(analysisTree.muon_px[index1],analysisTree.muon_py[index1],analysisTree.pfmetcorr_ex,analysisTree.pfmetcorr_ey);
                        
                        if(!isData)//only temporary for these round's MC
                        {
                            met = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
                            dPhiMETMuon = dPhiFrom2P(analysisTree.muon_px[index1],analysisTree.muon_py[index1],analysisTree.pfmet_ex,analysisTree.pfmet_ey);
                        }
        
                        
                        mt_1 = TMath::Sqrt(2*met*analysisTree.muon_pt[index1]*(1-TMath::Cos(dPhiMETMuon)));
            
                        TLorentzVector muon1;
                        muon1.SetXYZM(analysisTree.muon_px[index1],analysisTree.muon_py[index1],analysisTree.muon_pz[index1],muonMass);
                        //TLorentzVector tagmuon1_Up;
                        //tagmuon1_Up.SetXYZM((1+tagmuonScale)*analysisTree.muon_px[index1],(1+tagmuonScale)*analysisTree.muon_py[index1],(1+tagmuonScale)*analysisTree.muon_pz[index1],muonMass);
                        //TLorentzVector tagmuon1_Down;
                        //tagmuon1_Down.SetXYZM((1-tagmuonScale)*analysisTree.muon_px[index1],(1-tagmuonScale)*analysisTree.muon_py[index1],(1-tagmuonScale)*analysisTree.muon_pz[index1],muonMass);
                        
                        TLorentzVector tau2;
                        tau2.SetXYZM(analysisTree.tau_px[indexProbe],analysisTree.tau_py[indexProbe],analysisTree.tau_pz[indexProbe],analysisTree.tau_mass[indexProbe]);
                        //TLorentzVector probemuon2_Up;
                        //probemuon2_Up.SetXYZM((1+probemuonScale)*analysisTree.tau_px[indexProbe],(1+probemuonScale)*analysisTree.tau_py[indexProbe],(1+probemuonScale)*analysisTree.tau_pz[indexProbe],(1+probemuonScale)*analysisTree.tau_mass[indexProbe]);
                        //TLorentzVector probemuon2_Down;
                        //probemuon2_Down.SetXYZM((1-probemuonScale)*analysisTree.tau_px[indexProbe],(1-probemuonScale)*analysisTree.tau_py[indexProbe],(1-probemuonScale)*analysisTree.tau_pz[indexProbe],(1-probemuonScale)*analysisTree.tau_mass[indexProbe]);
                        
                        TLorentzVector tau2_Up;
                        tau2_Up.SetXYZM((1+probetauScale)*analysisTree.tau_px[indexProbe],(1+probetauScale)*analysisTree.tau_py[indexProbe],(1+probetauScale)*analysisTree.tau_pz[indexProbe],(1+probetauScale)*analysisTree.tau_mass[indexProbe]);
                        TLorentzVector tau2_Down;
                        tau2_Down.SetXYZM((1-probetauScale)*analysisTree.tau_px[indexProbe],(1-probetauScale)*analysisTree.tau_py[indexProbe],(1-probetauScale)*analysisTree.tau_pz[indexProbe],(1-probetauScale)*analysisTree.tau_mass[indexProbe]);
                        
                        //filling the tree
                        
                        //MC (pile up& scale factor)weight variable
                        if (!isData)
                        {
                            mcweight = analysisTree.genweight;
                            isoweight_1 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[index1]),double(analysisTree.muon_eta[index1]));
                            //FIX ME next time for all Summer16 MC
                            trigweight_1 = (float)SF_muonTrig->get_EfficiencyData(double(analysisTree.muon_pt[index1]),double(analysisTree.muon_eta[index1]));
                            if(isDY)
                            {
                                trigweight_1 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[index1]),double(analysisTree.muon_eta[index1]));
   
                            }
                            //      cout << "isoweight_1 = " << isoweight_1 << endl;
                            effweight = isoweight_1*trigweight_1;
                            puweight = float (PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
                        }
                        //cout << "Pile up weight: " << puweight<< endl;
                        float absIso = analysisTree.muon_r04_sumChargedHadronPt[index1];
                        float neutralIso = analysisTree.muon_r04_sumNeutralHadronEt[index1] + analysisTree.muon_r04_sumPhotonEt[index1] - 0.5*analysisTree.muon_r04_sumPUPt[index1];
                        neutralIso = TMath::Max(float(0),neutralIso);
                        absIso += neutralIso;
                        
                        iso_1 = absIso/analysisTree.muon_pt[index1];
                        //os or ss
                        float paircharge = q1*q2;
                        if(paircharge == -1) os = true;
                        if(paircharge == 1) os = false;
                        
                        m_vis = (muon1+tau2).M();
                        //m_vis_tagmuon_scale_up = (tagmuon1_Up+tau2).M();
                        //m_vis_tagmuon_scale_down = (tagmuon1_Down+tau2).M();
                        //m_vis_probemuon_scale_up = (muon1+probemuon2_Up).M();
                        //m_vis_probemuon_scale_down = (muon1+probemuon2_Down).M();
                        if(TPmatching_status==2)
                        {
                            m_vis_probetau_scale_up = (muon1+tau2_Up).M();
                            m_vis_probetau_scale_down = (muon1+tau2_Down).M();
                        }
                        else
                        {
                            m_vis_probetau_scale_up = m_vis;
                            m_vis_probetau_scale_down = m_vis;
                        }
                        if(TPmatching_status==1)
                        {
                            m_vis_gen = (genMuon1+genMuon2).M();
                            //cout <<"m_vis_gen:"<< m_vis_gen<< endl;
                            
                            //cout <<"----------" << endl;
                            m_vis_reso_scale_up = m_vis_gen + (1+resoScale)*(m_vis-m_vis_gen);
                            m_vis_reso_scale_down = m_vis_gen + (1-resoScale)*(m_vis-m_vis_gen);
                            //cout <<"m_vis_reso_up:"<< m_vis_reso_scale_up<< endl;
                            //cout <<"m_vis_reso_down:"<< m_vis_reso_scale_down<< endl;
                            //cout <<"--------------------------------------------" << endl;

                        }
                        else
                        {
                            m_vis_reso_scale_up = m_vis;
                            m_vis_reso_scale_down = m_vis;
                        }
                        
                        //cout <<"m_vis_gen:"<< m_vis_gen<< endl;
                        
                        //cout <<"----------" << endl;
                        
                        tauagainstMuonLoose = analysisTree.tau_againstMuonLoose3[indexProbe];
                        tauagainstMuonTight = analysisTree.tau_againstMuonTight3[indexProbe];
                        taubyLooseCombinedIsolationDeltaBetaCorr3Hits = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexProbe];
                        PtTag = analysisTree.muon_pt[index1];
                        EtaTag = analysisTree.muon_eta[index1];
                        PtProbe = analysisTree.tau_pt[indexProbe];
                        EtaProbe = analysisTree.tau_eta[indexProbe];
                        
                        muonTree->Fill();
                        selEventsIsoMuons++;
                    }
                   
                }
            }
      }

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



