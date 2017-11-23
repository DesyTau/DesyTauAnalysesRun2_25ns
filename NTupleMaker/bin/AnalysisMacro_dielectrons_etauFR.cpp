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

//A macro to perform e->tau fake rate measurement

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

    const bool isData = cfg.get<bool>("IsData");
    const bool isDY = cfg.get<bool>("IsDY");
    const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

    const string jsonFile = cfg.get<string>("jsonFile");
    
    // pile up reweighting
    const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");
    const string dataPUFile = cfg.get<string>("DataPUFile");
    const string mcPUFile = cfg.get<string>("MCPUFile");


    // kinematic cuts on electrons
    const float ptEleTagCut  = cfg.get<float>("ptEleTagCut");
    const float etaEleTagCut = cfg.get<float>("etaEleTagCut");
    const float dxyEleTagCut     = cfg.get<float>("dxyEleTagCut");
    const float dzEleTagCut      = cfg.get<float>("dzEleTagCut");
    const float isoEleTagCut = cfg.get<float>("isoEleTagCut");

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
    const string eleTriggerName  = cfg.get<string>("EleTriggerName");
    const string eleFilterName   = cfg.get<string>("EleFilterName");
    const string singleEleFilterName = cfg.get<string>("SingleEleFilterName");
    const float singleEleTriggerPtCut = cfg.get<float>("SingleEleTriggerPtCut");
    const float singleEleTriggerEtaCut = cfg.get<float>("SingleEleTriggerEtaCut");
  

    TString EleTriggerName(eleTriggerName);
    TString EleFilterName(eleFilterName);

    TString SingleEleFilterName(singleEleFilterName);

    // vertex cuts
    const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
    const float zVertexCut     = cfg.get<float>("ZVertexCut");
    const float dVertexCut     = cfg.get<float>("DVertexCut");
    
    //Perform eles Gen-Level matching
    //const bool genMatching = cfg.get<bool>("genMatching");
    
    //scale factor
    const string EleIdIsoFile = cfg.get<string>("EleIdIsoFile");
    const string EleTriggerFile = cfg.get<string>("EleTriggerFile");
    
    //Tag ele energy scale
    const float tageleScaleBarrel = cfg.get<float>("TagEleScaleBarrel");
    const float tageleScaleEndcap = cfg.get<float>("TagEleScaleEndcap");
    
    //Probe ele energy scale
    const float probeeleScale = cfg.get<float>("ProbeEleScale");
    
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
    TTree * eleTree = new TTree("ETauFR","ETauFR");

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
    
    Int_t           TPmatching_status;
    
    Float_t m_vis;
    Float_t m_vis_gen;
    Float_t m_vis_reso_scale_up;
    Float_t m_vis_reso_scale_down;
    
    Float_t m_vis_tagele_scale_up;
    Float_t m_vis_tagele_scale_down;
    
    Float_t m_vis_probeele_scale_up;
    Float_t m_vis_probeele_scale_down;
    
    Float_t m_vis_probetau_scale_up;
    Float_t m_vis_probetau_scale_down;
    
    Float_t m_vis_Up[50];
    Float_t m_vis_Down[50];
    
    Float_t mt_1;
    Bool_t os;
    Float_t iso_1;
    
    Bool_t tauagainstEleVLoose;
    Bool_t tauagainstEleLoose;
    Bool_t tauagainstEleMedium;
    Bool_t tauagainstEleTight;
    Bool_t tauagainstEleVTight;

    Bool_t taubyLooseCombinedIsolationDeltaBetaCorr3Hits;
    Bool_t taubyLooseIsolationMVArun2v1DBoldDMwLT;
    Bool_t taubyTightIsolationMVArun2v1DBoldDMwLT;

    
    Float_t PtTag;
    Float_t EtaTag;
    Float_t PtProbe;
    Float_t EtaProbe;
    Int_t   DecayModeProbe;
    
    Float_t met;
    UInt_t  npartons;
    
    
    eleTree->Branch("isZLL",&isZLL,"isZLL/O");
    eleTree->Branch("isZEE",&isZEE,"isZEE/O");
    eleTree->Branch("isZMM",&isZMM,"isZMM/O");
    eleTree->Branch("isZTT",&isZTT,"isZTT/O");
    
    eleTree->Branch("mcweight",&mcweight,"mcweight/F");
    eleTree->Branch("puweight", &puweight, "puweight/F");
    eleTree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
    eleTree->Branch("idweight_1", &idweight_1, "idweight_1/F");
    eleTree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
    eleTree->Branch("effweight", &effweight, "effweight/F");
    
    eleTree->Branch("TPmatching_status",&TPmatching_status,"TPmatching_status/I");
    
    eleTree->Branch("m_vis",&m_vis,"m_vis/F");
    eleTree->Branch("m_vis_gen",&m_vis_gen,"m_vis_gen/F");
    eleTree->Branch("m_vis_reso_scale_up",&m_vis_reso_scale_up,"m_vis_reso_scale_up/F");
    eleTree->Branch("m_vis_reso_scale_down",&m_vis_reso_scale_down,"m_vis_reso_scale_down/F");
    
    eleTree->Branch("m_vis_tagele_scale_up",&m_vis_tagele_scale_up,"m_vis_tagele_scale_up/F");
    eleTree->Branch("m_vis_tagele_scale_down",&m_vis_tagele_scale_down,"m_vis_tagele_scale_down/F");
    eleTree->Branch("m_vis_probeele_scale_up",&m_vis_probeele_scale_up,"m_vis_probeele_scale_up/F");
    eleTree->Branch("m_vis_probeele_scale_down",&m_vis_probeele_scale_down,"m_vis_probeele_scale_down/F");
    eleTree->Branch("m_vis_probetau_scale_up",&m_vis_probetau_scale_up,"m_vis_probetau_scale_up/F");
    eleTree->Branch("m_vis_probetau_scale_down",&m_vis_probetau_scale_down,"m_vis_probetau_scale_down/F");
    eleTree->Branch("m_vis_Up",&m_vis_Up,"m_vis_Up[50]/F");
    eleTree->Branch("m_vis_Down",&m_vis_Down,"m_vis_Down[50]/F");

    
    eleTree->Branch("mt_1",&mt_1,"mt_1/F");
    eleTree->Branch("os",&os,"os/O");
    eleTree->Branch("iso_1",&iso_1,"iso_1/F");
    
    eleTree->Branch("tauagainstEleVLoose",&tauagainstEleVLoose,"tauagainstEleVLoose/O");
    eleTree->Branch("tauagainstEleLoose",&tauagainstEleLoose,"tauagainstEleLoose/O");
    eleTree->Branch("tauagainstEleMedium",&tauagainstEleMedium,"tauagainstEleMedium/O");
    eleTree->Branch("tauagainstEleTight",&tauagainstEleTight,"tauagainstEleTight/O");
    eleTree->Branch("tauagainstEleVTight",&tauagainstEleVTight,"tauagainstEleVTight/O");

    eleTree->Branch("taubyLooseCombinedIsolationDeltaBetaCorr3Hits",&taubyLooseCombinedIsolationDeltaBetaCorr3Hits,"taubyLooseCombinedIsolationDeltaBetaCorr3Hits/O");
    eleTree->Branch("taubyLooseIsolationMVArun2v1DBoldDMwLT",&taubyLooseIsolationMVArun2v1DBoldDMwLT,"taubyLooseIsolationMVArun2v1DBoldDMwLT/O");
    eleTree->Branch("taubyTightIsolationMVArun2v1DBoldDMwLT",&taubyTightIsolationMVArun2v1DBoldDMwLT,"taubyTightIsolationMVArun2v1DBoldDMwLT/O");

    eleTree->Branch("PtTag",&PtTag,"PtTag/F");
    eleTree->Branch("EtaTag",&EtaTag,"EtaTag/F");
    eleTree->Branch("PtProbe",&PtProbe,"PtProbe/F");
    eleTree->Branch("EtaProbe",&EtaProbe,"EtaProbe/F");
    eleTree->Branch("DecayModeProbe",&DecayModeProbe,"DecayModeProbe/I");

    eleTree->Branch("met",&met,"met/F");
    eleTree->Branch("npartons",&npartons,"npartons/i");
    
  TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
  TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

  // PILE UP REWEIGHTING - OPTIONS


  // reweighting official recipe 
  // initialize pile up object
  PileUp * PUofficial = new PileUp();
  
  if (applyPUreweighting_official) {
    TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(dataPUFile),"read");
    TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(mcPUFile), "read");
    TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
    TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
    PUofficial->set_h_data(PU_data); 
    PUofficial->set_h_MC(PU_mc);
  }
    // ele scale factors
    ScaleFactor * SF_eleIdIso = new ScaleFactor();
    SF_eleIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(EleIdIsoFile));
    ScaleFactor * SF_eleTrig = new ScaleFactor();
    SF_eleTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(EleTriggerFile));

  int nFiles = 0;
  int nEvents = 0;
  int selEventsIsoEles = 0;

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

        npartons = analysisTree.genparticles_noutgoing;

        
        
      if (!isData) {

          //	cout << analysisTree.numtruepileupinteractions << endl;


        if (applyPUreweighting_official) {
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
      
      bool isTriggerEle = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(EleTriggerName)) {
	  //  std::cout << it->first << " : " << it->second << std::endl;
	  if (it->second==1)
	    isTriggerEle = true;
	}
      }

      if (applyTrigger && !isTriggerEle) continue;
    
      unsigned int nEleFilter = 0;
      bool isEleFilter = false;

      unsigned int nSingleEleFilter = 0;
      bool isSingleEleFilter = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==EleFilterName) {
	  nEleFilter = i;
	  isEleFilter = true;
	}
	if (HLTFilter==SingleEleFilterName) {
          nSingleEleFilter = i;
          isSingleEleFilter = true;
        }
      }
      if (!isEleFilter && isData) {
	cout << "Filter " << EleFilterName << " not found " << endl;
	exit(-1);
      }
      if (!isSingleEleFilter && isData) {
	cout << "Filter " << SingleEleFilterName << " not found " << endl;
        exit(-1);
      }

        // vertex cuts
        if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
        if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
        float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);
        if (dVertex>dVertexCut) continue;
      
      // eletron selection

      vector<unsigned int> allEles; allEles.clear();
      vector<unsigned int> idEles; idEles.clear();
      vector<unsigned int> isoEles; isoEles.clear();
      vector<float> isoElesValue; isoElesValue.clear();
      vector<float> allElesIso; allElesIso.clear();
      vector<bool> isElePassedIdIso; isElePassedIdIso.clear();
      vector<bool> isEleMatchedSingleEleFilter; isEleMatchedSingleEleFilter.clear();
        
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	bool elePassed    = true;
	bool eleSingleMatched = false;
	if (analysisTree.electron_pt[ie]<ptEleTagCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaEleTagCut) continue;
	allEles.push_back(ie);
	if (fabs(analysisTree.electron_dxy[ie])>dxyEleTagCut) elePassed = false;
	if (fabs(analysisTree.electron_dz[ie])>dzEleTagCut) elePassed = false;
	if (!analysisTree.electron_mva_wp80_general_Spring16_v1[ie]) elePassed = false;
    if (!analysisTree.electron_pass_conversion[ie]) elePassed = false;
	if (elePassed) idEles.push_back(ie);
	float absIso = 0;
	if (isoDR03) { 
        absIso = analysisTree.electron_r03_sumChargedHadronPt[ie];
	  float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[ie] +
	    analysisTree.electron_r03_sumPhotonEt[ie] -
	    0.5*analysisTree.electron_r03_sumPUPt[ie];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  absIso += neutralIso;
	}
	else {
	  absIso = analysisTree.electron_chargedHadIso[ie];
          float neutralIso = analysisTree.electron_neutralHadIso[ie] +
            analysisTree.electron_photonIso[ie] -
            0.5*analysisTree.electron_puIso[ie];
          neutralIso = TMath::Max(float(0),neutralIso);
          absIso += neutralIso;
	}
	float relIso = absIso/analysisTree.electron_pt[ie];
	allElesIso.push_back(relIso);
	if (elePassed && relIso<isoEleTagCut) {
	  isoEles.push_back(ie);
	  isoElesValue.push_back(relIso);
	}
	if (relIso>isoEleTagCut) elePassed = false;
	//	cout << "pt:" << analysisTree.electron_pt[ie] << "  passed:" << elePassed << endl;
	isElePassedIdIso.push_back(elePassed);
      }
        
        //tau selection
        vector<unsigned int> goodTaus; goodTaus.clear();
        for(unsigned int it=0;it<analysisTree.tau_count;++it)
        {
            if (analysisTree.tau_decayModeFinding[it]<=0.5) continue;
            if (analysisTree.tau_pt[it]<ptTauCut) continue;
            if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
            if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
            if (fabs(analysisTree.tau_charge[it])<0.5||fabs(analysisTree.tau_charge[it])>1.5) continue;
            if (analysisTree.tau_againstMuonLoose3[it] < 0.5) continue;
            goodTaus.push_back(it);
        }
        
        vector<unsigned int> genEles; genEles.clear();
        if(!isData)
        {
            for(unsigned int igenEle=0;igenEle<analysisTree.genparticles_count;++igenEle)
            {
                //cout <<"pdgid:"<<analysisTree.genparticles_pdgid[igenEle]<< endl;
                if (fabs(analysisTree.genparticles_pdgid[igenEle])==11 && analysisTree.genparticles_status[igenEle]==1 && analysisTree.genparticles_fromHardProcess[igenEle])
                {
                    genEles.push_back(igenEle);
                }
            }
        }
        //cout<<"size of gen electrons:"<< genEles.size()<<endl;

        
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
        //cout<< "size of gen eles:"<<genEles.size()<<endl;
        
        ////////////////////////////////////////////////////////////////////////////////////////
        bool isIsoElesPair = false;
        float isoMin = 9999;
        if (isoEles.size()>0) {
            for (unsigned int ie1=0; ie1<isoEles.size(); ++ie1) {
                unsigned int index1 = isoEles[ie1];
                bool isEle1matched = false;
                //if(isData)
                //{
                    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
                    {
                        float dRtrig = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                        if (dRtrig>DRTrigMatch) continue;
                        if (analysisTree.trigobject_filters[iT][nEleFilter] && analysisTree.electron_pt[index1] > ptEleTagCut && fabs(analysisTree.electron_eta[index1]) < etaEleTagCut)
                        isEle1matched = true;
                    }
                //}
                bool genRecoTagEleMatched = false;
                TLorentzVector genEle1;
                float indexGen1 = -9999;
                if(!isData)
                {
                    for(unsigned int iEleGen=0; iEleGen<genEles.size();++iEleGen)
                    {
                        unsigned int indexGenEle = genEles[iEleGen];
                        float genEleEta = PtoEta(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle],analysisTree.genparticles_pz[indexGenEle]);
                        float genElePhi = PtoPhi(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle]);
                        float dRRecoGenTag = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],genEleEta,genElePhi);
                        if(dRRecoGenTag<0.1)
                        {
                            genRecoTagEleMatched = true;
                            genEle1.SetXYZM(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle],analysisTree.genparticles_pz[indexGenEle],electronMass);
                            //cout << "gen Ele 1 px:" << analysisTree.genparticles_px[indexGenEle]<< endl;
                            indexGen1 = indexGenEle;
                        }
                    }
                }
                //if (isEle1matched || (!isData)) {
                if (isEle1matched) {
                    for (unsigned int iTau=0; iTau<goodTaus.size(); ++iTau) {
                        unsigned int indexProbe = goodTaus[iTau];
                        float q1 = analysisTree.electron_charge[index1];
                        float q2 = analysisTree.tau_charge[indexProbe];
                        float dR = deltaR(analysisTree.electron_eta[index1],analysisTree.electron_phi[index1],analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe]);
                        if (dR<dRleptonsCut) continue;

                        bool genRecoProbeEleMatched = false;
                        TLorentzVector genEle2;
                        for(unsigned int iEleGen=0; iEleGen<genEles.size();++iEleGen)
                        {
                            unsigned int indexGenEle = genEles[iEleGen];
                            if(indexGenEle!=indexGen1)
                            {
                                
                                float genEleEta = PtoEta(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle],analysisTree.genparticles_pz[indexGenEle]);
                                float genElePhi = PtoPhi(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle]);
                                float dRRecoGenProbe = deltaR(analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe],genEleEta,genElePhi);
                                if(dRRecoGenProbe<0.1)
                                {
                                    genRecoProbeEleMatched = true;
                                    genEle2.SetXYZM(analysisTree.genparticles_px[indexGenEle],analysisTree.genparticles_py[indexGenEle],analysisTree.genparticles_pz[indexGenEle],electronMass);
                                        //cout << "gen Ele 2 px:" << analysisTree.genparticles_px[indexGenEle]<< endl;
                                }
                            }
                        }
                        
                        TLorentzVector genTau2;
                        bool genRecoProbeTauMatched = false;
                        if(isZTT)//to match with GEN Tau object, building template for tau energy shapes systematics
                        {
                            //cout << "analysisTree.gentau_count:"<< analysisTree.gentau_count <<endl;
                            for(unsigned int igentau = 0; igentau < analysisTree.gentau_count; igentau++)
                            {
                                float genTauEta = PtoEta(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau],analysisTree.gentau_pz[igentau]);
                                float genTauPhi = PtoPhi(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau]);
                                float dRRecoGenTauProbe = deltaR(analysisTree.tau_eta[indexProbe],analysisTree.tau_phi[indexProbe],genTauEta,genTauPhi);
                                if(dRRecoGenTauProbe < 0.1)
                                {
                                    genRecoProbeTauMatched = true;
                                    genTau2.SetXYZM(analysisTree.gentau_px[igentau],analysisTree.gentau_py[igentau],analysisTree.gentau_pz[igentau],analysisTree.tau_mass[igentau]);
                                }
                            }
                            //cout << "genRecoProbeTauMatched:" << genRecoProbeTauMatched <<endl;
                            
                        }
                        
                        //setting tag and probe matching status
                        TPmatching_status = 3;
                        if(genRecoTagEleMatched && !genRecoProbeEleMatched && !genRecoProbeTauMatched)
                            TPmatching_status = 0;
                        
                        if(genRecoTagEleMatched && genRecoProbeEleMatched)
                            TPmatching_status = 1;
                        
                        if(genRecoProbeTauMatched)
                            TPmatching_status = 2;
                        
                       //cout << "TP matching status:" << TPmatching_status << endl;
                        
                        // choosing mva met
                        unsigned int metEleTau = 0;
                        bool mvaMetFound = false;
                        for (unsigned int iMet=0; iMet<analysisTree.mvamet_count; ++iMet)
                        {
                            //cout << iMet << endl;
                            if (analysisTree.mvamet_channel[iMet]==3)
                            {
                                if (analysisTree.mvamet_lep1[iMet]==indexProbe &&analysisTree.mvamet_lep2[iMet]==index1)
                                {
                                    metEleTau = iMet;
                                    mvaMetFound = true;
                                }
                            }
                        }
                        
                        //using uncorrected MET for now in 2017
                        
                        met = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
                        float dPhiMETEle = dPhiFrom2P(analysisTree.electron_px[index1],analysisTree.electron_py[index1],analysisTree.pfmet_ex,analysisTree.pfmet_ey);

                        mt_1 = TMath::Sqrt(2*met*analysisTree.electron_pt[index1]*(1-TMath::Cos(dPhiMETEle)));
            
                        TLorentzVector ele1;
                        ele1.SetXYZM(analysisTree.electron_px[index1],analysisTree.electron_py[index1],analysisTree.electron_pz[index1],electronMass);
                        
                        TLorentzVector tagele1_Up;
                        TLorentzVector tagele1_Down;
                        if(fabs(analysisTree.electron_eta[index1])<1.460)
                        {
                            tagele1_Up.SetXYZM((1+tageleScaleBarrel)*analysisTree.electron_px[index1],(1+tageleScaleBarrel)*analysisTree.electron_py[index1],(1+tageleScaleBarrel)*analysisTree.electron_pz[index1],electronMass);
                            tagele1_Down.SetXYZM((1-tageleScaleBarrel)*analysisTree.electron_px[index1],(1-tageleScaleBarrel)*analysisTree.electron_py[index1],(1-tageleScaleBarrel)*analysisTree.electron_pz[index1],electronMass);
                        }
                        if(fabs(analysisTree.electron_eta[index1])>1.558)
                        {
                            tagele1_Up.SetXYZM((1+tageleScaleEndcap)*analysisTree.electron_px[index1],(1+tageleScaleEndcap)*analysisTree.electron_py[index1],(1+tageleScaleEndcap)*analysisTree.electron_pz[index1],electronMass);
                            tagele1_Down.SetXYZM((1-tageleScaleEndcap)*analysisTree.electron_px[index1],(1-tageleScaleEndcap)*analysisTree.electron_py[index1],(1-tageleScaleEndcap)*analysisTree.electron_pz[index1],electronMass);
                        }
                        if((fabs(analysisTree.electron_eta[index1])<1.558) && (fabs(analysisTree.electron_eta[index1])>1.460))
                        {
                            tagele1_Up.SetXYZM(analysisTree.electron_px[index1],analysisTree.electron_py[index1],analysisTree.electron_pz[index1],electronMass);
                            tagele1_Down.SetXYZM(analysisTree.electron_px[index1],analysisTree.electron_py[index1],analysisTree.electron_pz[index1],electronMass);
                        }
                        
                        
                        TLorentzVector tau2;
                        tau2.SetXYZM(analysisTree.tau_px[indexProbe],analysisTree.tau_py[indexProbe],analysisTree.tau_pz[indexProbe],analysisTree.tau_mass[indexProbe]);
                        TLorentzVector probeele2_Up;//electron fake tau energy scale
                        TLorentzVector probeele2_Down;

                        TLorentzVector EleFakeTau_Up;// electron fake tau energy scale study
                        TLorentzVector EleFakeTau_Down;// electron fake tau energy scale study

                        for(int i=0;i<50;++i)
                        {
                            EleFakeTau_Up.SetXYZM((1.0+probeeleScale*(i/50.0))*analysisTree.tau_px[indexProbe],(1.0+probeeleScale*(i/50.0))*analysisTree.tau_py[indexProbe],(1.0+probeeleScale*(i/50.0))*analysisTree.tau_pz[indexProbe],(1.0+probeeleScale*(i/50.0))*analysisTree.tau_mass[indexProbe]);
                            EleFakeTau_Down.SetXYZM((1.0-probeeleScale*(i/50.0))*analysisTree.tau_px[indexProbe],(1.0-probeeleScale*(i/50.0))*analysisTree.tau_py[indexProbe],(1.0-probeeleScale*(i/50.0))*analysisTree.tau_pz[indexProbe],(1.0-probeeleScale*(i/50.0))*analysisTree.tau_mass[indexProbe]);
                            m_vis_Up[i] = (ele1+EleFakeTau_Up).M();
                            m_vis_Down[i] = (ele1+EleFakeTau_Down).M();

                        }

                        probeele2_Up.SetXYZM((1+probeeleScale)*analysisTree.tau_px[indexProbe],(1+probeeleScale)*analysisTree.tau_py[indexProbe],(1+probeeleScale)*analysisTree.tau_pz[indexProbe],(1+probeeleScale)*analysisTree.tau_mass[indexProbe]);
                        probeele2_Down.SetXYZM((1-probeeleScale)*analysisTree.tau_px[indexProbe],(1-probeeleScale)*analysisTree.tau_py[indexProbe],(1-probeeleScale)*analysisTree.tau_pz[indexProbe],(1-probeeleScale)*analysisTree.tau_mass[indexProbe]);
                        
                        TLorentzVector tau2inZTT_Up;
                        tau2inZTT_Up.SetXYZM((1+probetauScale)*analysisTree.tau_px[indexProbe],(1+probetauScale)*analysisTree.tau_py[indexProbe],(1+probetauScale)*analysisTree.tau_pz[indexProbe],(1+probetauScale)*analysisTree.tau_mass[indexProbe]);
                        TLorentzVector tau2inZTT_Down;
                        tau2inZTT_Down.SetXYZM((1-probetauScale)*analysisTree.tau_px[indexProbe],(1-probetauScale)*analysisTree.tau_py[indexProbe],(1-probetauScale)*analysisTree.tau_pz[indexProbe],(1-probetauScale)*analysisTree.tau_mass[indexProbe]);
                        
                        //filling the tree
                        
                        //MC (pile up& scale factor)weight variable
                        if (!isData)
                        {
                            mcweight = analysisTree.genweight;
                            isoweight_1 = (float)SF_eleIdIso->get_ScaleFactor(double(analysisTree.electron_pt[index1]),double(analysisTree.electron_eta[index1]));
                            //trigweight_1 = (float)SF_eleTrig->get_EfficiencyData(double(analysisTree.electron_pt[index1]),double(analysisTree.electron_eta[index1]));
                            //      cout << "isoweight_1 = " << isoweight_1 << endl;

                            trigweight_1 = (float)SF_eleTrig->get_ScaleFactor(double(analysisTree.electron_pt[index1]),double(analysisTree.electron_eta[index1]));
                            effweight = isoweight_1*trigweight_1;
                            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
                        }
                       
                        float absIso = analysisTree.electron_r03_sumChargedHadronPt[index1];
                        float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[index1] + analysisTree.electron_r03_sumPhotonEt[index1] - 0.5*analysisTree.electron_r03_sumPUPt[index1];
                        neutralIso = TMath::Max(float(0),neutralIso);
                        absIso += neutralIso;
                        
                        iso_1 = absIso/analysisTree.electron_pt[index1];
                        //os or ss
                        float paircharge = q1*q2;
                        if(paircharge == -1) os = true;
                        if(paircharge == 1) os = false;
                        
                        m_vis = (ele1+tau2).M();
                        
                        if(TPmatching_status==1 || TPmatching_status==2)
                        {
                            m_vis_tagele_scale_up = (tagele1_Up+tau2).M();
                            m_vis_tagele_scale_down = (tagele1_Down+tau2).M();

                        }
                        else
                        {
                            m_vis_tagele_scale_up = m_vis;
                            m_vis_tagele_scale_down = m_vis;
                        }
                        
                        if(TPmatching_status==1)
                        {
                            m_vis_probeele_scale_up = (ele1+probeele2_Up).M();
                            m_vis_probeele_scale_down = (ele1+probeele2_Down).M();

                        }
                        else
                        {
                            m_vis_probeele_scale_up = m_vis;
                            m_vis_probeele_scale_down = m_vis;

                        }
                        
                        
                        if(TPmatching_status==1)
                        {
                            m_vis_gen = (genEle1+genEle2).M();
                            m_vis_reso_scale_up = 90.99 + (1+resoScale)*(m_vis-90.99);
                            m_vis_reso_scale_down = 90.99 + (1-resoScale)*(m_vis-90.99);
                            //m_vis_reso_scale_up = m_vis_gen + (1+resoScale)*(m_vis-m_vis_gen);
                            //m_vis_reso_scale_down = m_vis_gen + (1-resoScale)*(m_vis-m_vis_gen);
                        }
                        else
                        {
                            m_vis_reso_scale_up = m_vis;
                            m_vis_reso_scale_down = m_vis;
                        }
                        
                        
                        if(TPmatching_status==2)
                        {
                            m_vis_probetau_scale_up = (ele1+tau2inZTT_Up).M();
                            m_vis_probetau_scale_down = (ele1+tau2inZTT_Down).M();
                        }
                        else
                        {
                            m_vis_probetau_scale_up = m_vis;
                            m_vis_probetau_scale_down = m_vis;
                        }
                        
                        //cout <<"m_vis_gen:"<< m_vis_gen<< endl;
                        
                        //cout <<"----------" << endl;
                        
                        tauagainstEleVLoose = analysisTree.tau_againstElectronVLooseMVA6[indexProbe];
                        tauagainstEleLoose = analysisTree.tau_againstElectronLooseMVA6[indexProbe];
                        tauagainstEleMedium = analysisTree.tau_againstElectronMediumMVA6[indexProbe];
                        tauagainstEleTight = analysisTree.tau_againstElectronTightMVA6[indexProbe];
                        tauagainstEleVTight = analysisTree.tau_againstElectronVTightMVA6[indexProbe];
                        taubyLooseCombinedIsolationDeltaBetaCorr3Hits = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexProbe];
                        taubyLooseIsolationMVArun2v1DBoldDMwLT = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexProbe];
                        taubyTightIsolationMVArun2v1DBoldDMwLT = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexProbe];

                        PtTag = analysisTree.electron_pt[index1];
                        EtaTag = analysisTree.electron_eta[index1];
                        PtProbe = analysisTree.electron_pt[indexProbe];
                        EtaProbe = analysisTree.electron_eta[indexProbe];
                        DecayModeProbe = analysisTree.tau_decayMode[indexProbe];

                        
                        eleTree->Fill();
                        selEventsIsoEles++;
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
  std::cout << "Total number of selected events (iso ele pairs) = " << selEventsIsoEles << std::endl;
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



