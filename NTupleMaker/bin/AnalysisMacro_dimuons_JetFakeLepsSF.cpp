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



int main(int argc, char * argv[])
{
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
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
    
    // kinematic cuts on muons
    const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
    const float ptMuonCut   = cfg.get<float>("ptMuonCut");
    const float etaMuonCut = cfg.get<float>("etaMuonCut");
    const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
    const float dzMuonCut      = cfg.get<float>("dzMuonCut");
    const float isoMuonCut     = cfg.get<float>("isoMuonCut");
    
    // kinematic cuts on electrons
    const float ptEleCut  = cfg.get<float>("ptEleCut");
    const float etaEleCut = cfg.get<float>("etaEleCut");
    const float dxyEleCut     = cfg.get<float>("dxyEleCut");
    const float dzEleCut      = cfg.get<float>("dzEleCut");
    const float isoEleCut     = cfg.get<float>("isoEleCut");
    
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
    
    // jet related cuts
    const float jetEtaCut      = cfg.get<float>("JetEtaCut");
    const float jetPtCut   = cfg.get<float>("JetPtCut");
    const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
    
    const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
    TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
    
    const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
    TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
    
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
    TTree * ZMMTree = new TTree("JetFakeLep","JetFakeLep");
    
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
    Float_t         zptmassweight;
    Float_t         iso_1;
    Float_t         iso_2;
    Int_t           JetMatchingMuonStatus;
    Int_t           JetMatchingEleStatus;
    
    Bool_t isBadMuonEvents;
    
    Float_t m_vis;

    Bool_t os;
    
    Float_t pt_1;
    Float_t eta_1;
    Float_t phi_1;
    Float_t mt_1;

    Float_t pt_2;
    Float_t eta_2;
    Float_t phi_2;
    Float_t mt_2;
    
    Float_t dilepton_pt;
    Float_t dilepton_eta;
    Float_t dilepton_phi;

    Float_t met;
    UInt_t  npartons;
    
    UInt_t  njets;
    Float_t jet_pt;//leading jets
    Float_t jet_eta;
    Float_t jet_phi;
    
    Float_t         nuPx;  // neutrinos from t -> W(lv) + b
    Float_t         nuPy;  // or from tau->l+v+v events
    Float_t         nuPz;  // (x,y,z) components
    Float_t         nuPt;  // pT
    Float_t         nuPhi; // phi
    
    Float_t         lepPx;
    Float_t         lepPy;
    Float_t         lepPz;
    
    Float_t         bosonPx;
    Float_t         bosonPy;
    Float_t         bosonPz;
    Float_t         bosonPt;
    Float_t         bosonEta;
    Float_t         bosonMass;
    
    ZMMTree->Branch("isZLL",&isZLL,"isZLL/O");
    ZMMTree->Branch("isZEE",&isZEE,"isZEE/O");
    ZMMTree->Branch("isZMM",&isZMM,"isZMM/O");
    ZMMTree->Branch("isZTT",&isZTT,"isZTT/O");
    
    ZMMTree->Branch("mcweight",&mcweight,"mcweight/F");
    ZMMTree->Branch("puweight", &puweight, "puweight/F");
    ZMMTree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");
    ZMMTree->Branch("idweight_1", &idweight_1, "idweight_1/F");
    ZMMTree->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
    ZMMTree->Branch("effweight", &effweight, "effweight/F");
    ZMMTree->Branch("zptmassweight",&zptmassweight,"zptmassweight/F");

    ZMMTree->Branch("iso_1",&iso_1,"iso_1/F");
    ZMMTree->Branch("iso_2",&iso_2,"iso_2/F");
    ZMMTree->Branch("JetMatchingMuonStatus",&JetMatchingMuonStatus,"JetMatchingMuonStatus/I");
    ZMMTree->Branch("JetMatchingEleStatus",&JetMatchingEleStatus,"JetMatchingEleStatus/I");

    ZMMTree->Branch("m_vis",&m_vis,"m_vis/F");
    
    ZMMTree->Branch("isBadMuonEvents",&isBadMuonEvents,"isBadMuonEvents/O");// badMuon events

    ZMMTree->Branch("os",&os,"os/O");

    ZMMTree->Branch("pt_1",&pt_1,"pt_1/F");
    ZMMTree->Branch("eta_1",&eta_1,"eta_1/F");
    ZMMTree->Branch("phi_1",&phi_1,"phi_1/F");

    ZMMTree->Branch("mt_1",&mt_1,"mt_1/F");

    ZMMTree->Branch("pt_2",&pt_2,"pt_2/F");
    ZMMTree->Branch("eta_2",&eta_2,"eta_2/F");
    ZMMTree->Branch("phi_2",&phi_2,"phi_2/F");

    ZMMTree->Branch("mt_2",&mt_2,"mt_2/F");
    
    ZMMTree->Branch("dilepton_pt",&dilepton_pt,"dilepton_pt/F");
    ZMMTree->Branch("dilepton_eta",&dilepton_eta,"dilepton_eta/F");
    ZMMTree->Branch("dilepton_phi",&dilepton_phi,"dilepton_phi/F");
    
    ZMMTree->Branch("met",&met,"met/F");
    ZMMTree->Branch("npartons",&npartons,"npartons/i");
    
    ZMMTree->Branch("njets",&njets,"njets/i");
    ZMMTree->Branch("jet_pt",&jet_pt,"jet_pt/F");
    ZMMTree->Branch("jet_eta",&jet_eta,"jet_eta/F");
    ZMMTree->Branch("jet_phi",&jet_phi,"jet_phi/F");
    
    ZMMTree->Branch("nuPx",&nuPx,"nuPx/F");
    ZMMTree->Branch("nuPy",&nuPy,"nuPy/F");
    ZMMTree->Branch("nuPz",&nuPz,"nuPz/F");
    ZMMTree->Branch("nuPt",&nuPt,"nuPt/F");
    ZMMTree->Branch("nuPhi",&nuPhi,"nuPhi/F");
    
    ZMMTree->Branch("lepPx",&lepPx,"lepPx/F");
    ZMMTree->Branch("lepPy",&lepPy,"lepPy/F");
    ZMMTree->Branch("lepPz",&lepPz,"lepPz/F");
    
    ZMMTree->Branch("bosonPx",&bosonPx,"bosonPx/F");
    ZMMTree->Branch("bosonPy",&bosonPy,"bosonPy/F");
    ZMMTree->Branch("bosonPz",&bosonPz,"bosonPz/F");
    ZMMTree->Branch("bosonPt",&bosonPt,"bosonPt/F");
    ZMMTree->Branch("bosonMass",&bosonMass,"bosonMass/F");
    
    
    TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
    TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);
    
    // reweighting official recipe
    // initialize pile up object
    PileUp * PUofficial = new PileUp();
    
    if (applyPUreweighting) {
        
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
    
    // Z pt mass weights
    TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName);
    if (fileZMassPtWeights->IsZombie()) {
        std::cout << "File " << TString(cmsswBase) << "/src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
        exit(-1);
    }
    TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName);
    if (histZMassPtWeights==NULL) {
        std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName
        << std::endl;
        exit(-1);
    }
    
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
    
    for (int iF=0; iF<nTotalFiles; ++iF)
    {
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
        for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++)
        {
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
            mcweight = 1;
            zptmassweight = 1;
            
            //------------------------------------------------
            if (!isData)
            {
                weight *=analysisTree.genweight;
                npartons = analysisTree.genparticles_noutgoing;
            }
            
            if (!isData)
            {
                //	cout << analysisTree.numtruepileupinteractions << endl;
                if (applyPUreweighting)
                {
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
            
            bool isNewRun = true;
            if (allRuns.size()>0)
            {
                for (unsigned int iR=0; iR<allRuns.size(); ++iR)
                {
                    if (analysisTree.event_run==allRuns.at(iR))
                    {
                        isNewRun = false;
                        break;
                    }
                }
            }
            
            //pile up and gen weight
            if (!isData)
            {
                puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
                mcweight = analysisTree.genweight;
            }
            
            if (isNewRun)
            allRuns.push_back(analysisTree.event_run);

            bool isTriggerMuon = false;
            for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it)
            {
                TString trigName(it->first);
                if (trigName.Contains(MuonTriggerName))
                {
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
            
            for (unsigned int i=0; i<nfilters; ++i)
            {
                TString HLTFilter(analysisTree.run_hltfilters->at(i));
                if (HLTFilter==MuonFilterName)
                {
                    nMuonFilter = i;
                    isMuonFilter = true;
                }
                if (HLTFilter==SingleMuonFilterName)
                {
                    nSingleMuonFilter = i;
                    isSingleMuonFilter = true;
                }
            }
            
            if (!isMuonFilter && isData)
            {
                cout << "Filter " << MuonFilterName << " not found " << endl;
                exit(-1);
            }
            if (!isSingleMuonFilter && isData)
            {
                cout << "Filter " << SingleMuonFilterName << " not found " << endl;
                exit(-1);
            }
            
            
            // vertex cuts
            if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
            if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
            float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);
            if (dVertex>dVertexCut) continue;
            
            
            // bad muon selection
            vector<int> badmuons; badmuons.clear();
            for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
                if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
                if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
                if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
                if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
                bool muonbad = analysisTree.muon_isBad[im];
                bool muonduplicate = analysisTree.muon_isDuplicate[im];
                if (muonbad || muonduplicate)
                    badmuons.push_back(im);
            }
            
            if (badmuons.size()>0)
            {
                isBadMuonEvents = 1;
                //cout << "bad muon number: " << badmuons.size() << endl;
            }
            else
            {
                isBadMuonEvents = 0;
            }
            
            //muon selection
            vector<unsigned int> allMuons; allMuons.clear();
            vector<unsigned int> idMuons; idMuons.clear();
            vector<unsigned int> isoMuons; isoMuons.clear();
            vector<float> isoMuonsValue; isoMuonsValue.clear();
            vector<float> allMuonsIso; allMuonsIso.clear();
            vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
            vector<bool> isMuonMatchedSingleMuFilter; isMuonMatchedSingleMuFilter.clear();
            for (unsigned int im = 0; im<analysisTree.muon_count; ++im)
            {
                bool muPassed    = true;
                bool muSingleMatched = false;
                if (analysisTree.muon_pt[im]<ptMuonLowCut) continue; //low cut 10 GeV
                if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
                allMuons.push_back(im);
                if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) muPassed = false;
                if (fabs(analysisTree.muon_dz[im])>dzMuonCut) muPassed = false;
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
                if (isoDR03)
                {
                    absIso = analysisTree.muon_r03_sumChargedHadronPt[im];
                    float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] + analysisTree.muon_r03_sumPhotonEt[im] - 0.5*analysisTree.muon_r03_sumPUPt[im];
                    neutralIso = TMath::Max(float(0),neutralIso);
                    absIso += neutralIso;
                }
                else
                {
                    absIso = analysisTree.muon_chargedHadIso[im];
                    float neutralIso = analysisTree.muon_neutralHadIso[im] + analysisTree.muon_photonIso[im] - 0.5*analysisTree.muon_puIso[im];
                    neutralIso = TMath::Max(float(0),neutralIso);
                    absIso += neutralIso;
                }
                
                float relIso = absIso/analysisTree.muon_pt[im];
                allMuonsIso.push_back(relIso);
                if (relIso>isoMuonCut) muPassed = false;
                if (muPassed)
                {
                    isoMuons.push_back(im);
                    isoMuonsValue.push_back(relIso);
                }
                isMuonPassedIdIso.push_back(muPassed);
                for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
                {
                    float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                    if (dRtrig>DRTrigMatch) continue;
                    if (analysisTree.trigobject_filters[iT][nSingleMuonFilter] && analysisTree.trigobject_pt[iT]>singleMuonTriggerPtCut && fabs(analysisTree.trigobject_eta[iT])<singleMuonTriggerEtaCut)
                        muSingleMatched = true;
                    if (!applyTrigger)
                    {
                        muSingleMatched = true;
                    }
                }
                isMuonMatchedSingleMuFilter.push_back(muSingleMatched);

            }
            
            //if(isoMuons.size()>2)
               // cout << "iso Muon numbers: " << isoMuons.size() << endl;
            
            //electron selection
            vector<unsigned int> allEles; allEles.clear();
            vector<unsigned int> idEles; idEles.clear();
            vector<unsigned int> isoEles; isoEles.clear();
            vector<float> isoElesValue; isoElesValue.clear();
            vector<float> allElesIso; allElesIso.clear();
            vector<bool> isElePassedIdIso; isElePassedIdIso.clear();
            
            for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie)
            {
                bool elePassed = true;
                if (analysisTree.electron_pt[ie]<ptEleCut) continue;
                if (fabs(analysisTree.electron_eta[ie])>etaEleCut) continue;
                allEles.push_back(ie);
                if (fabs(analysisTree.electron_dxy[ie])>dxyEleCut) elePassed = false;
                if (fabs(analysisTree.electron_dz[ie])>dzEleCut) elePassed = false;
                if (!analysisTree.electron_mva_wp80_general_Spring16_v1[ie]) elePassed = false;
                if (!analysisTree.electron_pass_conversion[ie]) elePassed = false;
                if (elePassed) idEles.push_back(ie);
                float absIso = 0;
                if (isoDR03)
                {
                    absIso = analysisTree.electron_r03_sumChargedHadronPt[ie];
                    float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[ie] + analysisTree.electron_r03_sumPhotonEt[ie] - 0.5*analysisTree.electron_r03_sumPUPt[ie];
                    neutralIso = TMath::Max(float(0),neutralIso);
                    absIso += neutralIso;
                }
                else
                {
                    absIso = analysisTree.electron_chargedHadIso[ie];
                    float neutralIso = analysisTree.electron_neutralHadIso[ie] + analysisTree.electron_photonIso[ie] - 0.5*analysisTree.electron_puIso[ie];
                    neutralIso = TMath::Max(float(0),neutralIso);
                    absIso += neutralIso;
                }
                float relIso = absIso/analysisTree.electron_pt[ie];
                allElesIso.push_back(relIso);
                if (elePassed && relIso<isoEleCut)
                {
                    isoEles.push_back(ie);
                    isoElesValue.push_back(relIso);
                }
                if (relIso>isoEleCut) elePassed = false;
                isElePassedIdIso.push_back(elePassed);

            }
            //if(isoEles.size()>0)
                //cout << "iso ele number:" << isoEles.size()<< endl;
            
            //cout << "-----------------------------------" << endl;
            
            //determining isZMM isZTT or isZEE and implementing zptmasweight
            std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
            std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
            std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
            std::vector<TLorentzVector> promptMuons; promptMuons.clear();
            std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
            std::vector<TLorentzVector> tauNeutrinos; tauNeutrinos.clear();
            
            TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0.001,0.001,0,0);
            TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
            TLorentzVector tauNeutrinosLV;  tauNeutrinosLV.SetXYZT(0,0,0,0);
            TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
            TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
            TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);
            
            TLorentzVector genBosonLV; genBosonLV.SetXYZT(0,0,0,0);
            TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);
            
            if(!isData)
            {
                // computing boson 4-vector
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
                        promptTausLV += tauLV;
                        wDecayProductsLV += tauLV;
                    }
                    if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {
                        promptTausLastCopy.push_back(tauVisLV);
                        promptVisTausLV += tauVisLV;
                    }
                }
                
                for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
                    
                    TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                        analysisTree.genparticles_py[igen],
                                                        analysisTree.genparticles_pz[igen],
                                                        analysisTree.genparticles_e[igen]);
                    
                    bool fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1;
                    bool isMuon = false;
                    bool isElectron = false;
                    bool isNeutrino = false;
                    bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==11) {
                        isElectron = true;
                        if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                            promptElectrons.push_back(genLV);
                            promptElectronsLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                    }
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==13) {
                        isMuon = true;
                        if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                            promptMuons.push_back(genLV);
                            promptMuonsLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                    }
                    
                    if (abs(analysisTree.genparticles_pdgid[igen])==12||
                        abs(analysisTree.genparticles_pdgid[igen])==14||
                        abs(analysisTree.genparticles_pdgid[igen])==16)  {
                        isNeutrino = true;
                        if ((analysisTree.genparticles_fromHardProcess[igen]||analysisTree.genparticles_isPrompt[igen])&&
                            !analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
                            analysisTree.genparticles_status[igen]==1) {
                            promptNeutrinos.push_back(genLV);
                            promptNeutrinosLV += genLV;
                            wDecayProductsLV += genLV;
                        }
                        if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
                            analysisTree.genparticles_status[igen]==1) {
                            tauNeutrinos.push_back(genLV);
                            tauNeutrinosLV += genLV;
                        }
                    }
                    
                    bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
                    bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);
                    
                    if (isBoson)
                        genBosonLV += genLV;
                    if (isVisibleBoson)
                        genVisBosonLV += genLV;
                    
                }
                
                if (isDY)
                {
                    if (promptTausFirstCopy.size()==2)
                    {
                        isZTT = true; isZMM = false; isZEE = false;
                        bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
                        bosonMass = promptTausLV.M();
                        bosonEta  = promptTausLV.Eta();
                        lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
                    }
                    else if (promptMuons.size()==2)
                    {
                        isZTT = false; isZMM = true; isZEE = false;
                        bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
                        bosonMass = promptMuonsLV.M();
                        bosonEta = promptMuonsLV.Eta();
                        lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
                    }
                    else
                    {
                        isZTT = false; isZMM = false; isZEE = true;
                        bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
                        bosonMass = promptElectronsLV.M();
                        bosonEta = promptElectronsLV.Eta();
                        lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
                    }
                    nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
                }
                
                bosonPx = genBosonLV.Px();
                bosonPy = genBosonLV.Py();
                bosonPz = genBosonLV.Pz();
                bosonPt = genBosonLV.Pt();
                bosonMass = genBosonLV.M();
                
                lepPx = genVisBosonLV.Px();
                lepPy = genVisBosonLV.Py();
                lepPz = genVisBosonLV.Pz();
                
                nuPt = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
                nuPhi = TMath::ATan2(nuPy,nuPx);
                
                if (isDY)// applying Z pt mass weights
                {
                    zptmassweight = 1;
                    if (bosonMass>50.0)
                    {
                        float bosonMassX = bosonMass;
                        float bosonPtX = bosonPt;
                        if (bosonMassX>1000.) bosonMassX = 1000.;
                        if (bosonPtX<1.)      bosonPtX = 1.;
                        if (bosonPtX>1000.)   bosonPtX = 1000.;
                        zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
                                                                          histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
                    }
                }
                
                ALL->Fill(0.0);
                if (isZMM) ZMM->Fill(0.);
                if (isZEE) ZEE->Fill(0.);
                if (isZTT) ZTT->Fill(0.);
                isZLL = isZMM || isZEE;
            }
            
            //select muon pair
            unsigned int indx1 = 0;
            unsigned int indx2 = 0;
            bool isIsoMuonsPair = false;
            bool firstTrigger = true;
            float isoMin = 9999;
            if (isoMuons.size()>0)
            {
                for (unsigned int im1=0; im1<isoMuons.size(); ++im1)
                {
                    unsigned int index1 = isoMuons[im1];
                    if(analysisTree.muon_pt[index1] < ptMuonCut) continue;
                    bool isMu1matched = false;
                    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
                    {
                        float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                        if (dRtrig>DRTrigMatch) continue;
                        if (analysisTree.trigobject_filters[iT][nMuonFilter] && analysisTree.muon_pt[index1] > ptMuonCut && fabs(analysisTree.muon_eta[index1]) < etaMuonCut)
                            isMu1matched = true;
                        if (!applyTrigger && analysisTree.muon_pt[index1] > ptMuonCut && fabs(analysisTree.muon_eta[index1]) < etaMuonCut)
                            isMu1matched = true;
                        
                    }
                    for (unsigned int im2=0; im2<isoMuons.size(); ++im2)
                    {
                        unsigned int index2 = isoMuons[im2];
                        if (index2 == index1) continue;
                        if(analysisTree.muon_pt[index1] < ptMuonCut) continue;
                        bool isMu2matched = false;
                        for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
                        {
                            float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
                                                  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                            if (dRtrig>DRTrigMatch) continue;
                            if (analysisTree.trigobject_filters[iT][nMuonFilter] &&
                                analysisTree.muon_pt[index2] > ptMuonCut &&
                                fabs(analysisTree.muon_eta[index2]) < etaMuonCut) 
                            isMu2matched = true;
                        }
                        if (!applyTrigger && analysisTree.muon_pt[index2] > ptMuonCut && fabs(analysisTree.muon_eta[index2]) < etaMuonCut)
                            isMu2matched = true;
                        
                        bool isTriggerMatch = false;
                        if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2])
                        {
                            isTriggerMatch = isMu1matched;
                            firstTrigger = true;
                        }
                        else
                        {
                            isTriggerMatch = isMu2matched;
                            firstTrigger = false;
                        }
                        float dRmumu = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
                        if (isTriggerMatch && dRmumu>dRleptonsCut)
                        {
                            bool sumIso = isoMuonsValue[im1]+isoMuonsValue[im2];
                            if (sumIso<isoMin)
                            {
                                isIsoMuonsPair = true;
                                isoMin = sumIso;
                                if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2])
                                {
                                    indx1 = index1;
                                    indx2 = index2;
                                }
                                else
                                {
                                    indx2 = index1;
                                    indx1 = index2;
                                }
                            }
                        }
                    }//end looping second isolated muons
                }//end looping first isolated muons
            }//end selecting muon pair
            
            if (!isIsoMuonsPair) continue;
            
            //fill the trigger weight if first lepton matched to single trigger
            if(!isData && isIsoMuonsPair && firstTrigger)
            {
                isoweight_1 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indx1]),double(analysisTree.muon_eta[indx1]));
                //cout << "isoweight_1 = " << isoweight_1 << endl;
                trigweight_1 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indx1]),double(analysisTree.muon_eta[indx1]));
                //cout << "trigweight_1 = " << trigweight_1 << endl;
                effweight = isoweight_1*trigweight_1;
            }

            //fill the trigger weight if second lepton matched to single trigger
            if(!isData && isIsoMuonsPair && (!firstTrigger))
            {
                isoweight_1 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indx2]),double(analysisTree.muon_eta[indx2]));
                trigweight_1 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indx2]),double(analysisTree.muon_eta[indx2]));
                effweight = isoweight_1*trigweight_1;
            }
            
            ///////////////////////
            //input dimuon
            TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[indx1],
                                                analysisTree.muon_py[indx1],
                                                analysisTree.muon_pz[indx1],
                                                muonMass);
                
            TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[indx2],
                                                analysisTree.muon_py[indx2],
                                                analysisTree.muon_pz[indx2],
                                                muonMass);
                
            TLorentzVector dimuon = mu1 + mu2;
                
            m_vis = dimuon.M();
            float dilepton_px = dimuon.Px();
            float dilepton_py = dimuon.Py();
            dilepton_pt = dimuon.Pt();
            dilepton_eta = dimuon.Eta();
            dilepton_phi = dimuon.Phi();
                
            JetMatchingMuonStatus = 0;
            JetMatchingEleStatus = 0;
            
            
            njets = 0;
            jet_pt = -9999;
            jet_eta = -9999;
            jet_phi = -9999;
            
            int jet_indx = -1;
            float dPhiJetZ = -1.0;
            
            for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
            {
                float jetEta = analysisTree.pfjet_eta[jet];
                float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                float jetPt = analysisTree.pfjet_pt[jet];
                
                if (analysisTree.pfjet_pt[jet]<jetPtCut) continue;
                if (absJetEta>jetEtaCut) continue;
                    
                float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1]);
                if (dR1<dRJetLeptonCut) continue;
                
                float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
                if (dR2<dRJetLeptonCut) continue;
                    
                // pfJetId
                //bool isPFJetId = looseJetiD(analysisTree,int(jet));
                //if (!isPFJetId) continue;
                
                // jet is required to be separated delta phi > 2.0
                
                float dPhiJetZ_temp = dPhiFrom2P(dilepton_px,dilepton_py,analysisTree.pfjet_px[jet],analysisTree.pfjet_py[jet]);
                
                
                if(dPhiJetZ_temp < 2.0) continue;
                
                njets++;
                
                if(dPhiJetZ < dPhiJetZ_temp)
                {
                    dPhiJetZ = dPhiJetZ_temp;
                    jet_indx = jet;
                }
                    
            }//end looping pf jet
            
            //filling up jet variables
        
            if(njets==0) continue;
            
            jet_pt = analysisTree.pfjet_pt[jet_indx];
            jet_eta = analysisTree.pfjet_eta[jet_indx];
            jet_phi = analysisTree.pfjet_phi[jet_indx];
            
            
            //define jet matching status
            for(unsigned int im=0; im<isoMuons.size(); ++im)
            {
                if((isoMuons[im] == indx1) || (isoMuons[im] == indx2)) continue;
                float dRJetMuon = deltaR(analysisTree.pfjet_eta[jet_indx],analysisTree.pfjet_phi[jet_indx],analysisTree.muon_eta[isoMuons[im]],analysisTree.muon_phi[isoMuons[im]]);
                if (dRJetMuon<0.5)
                    JetMatchingMuonStatus = 1;
            }
            
            for(unsigned int ie=0; ie<isoEles.size(); ++ie)
            {
                float dRJetEle = deltaR(analysisTree.pfjet_eta[jet_indx],analysisTree.pfjet_phi[jet_indx],analysisTree.electron_eta[isoEles[ie]],analysisTree.electron_phi[isoEles[ie]]);
                if(dRJetEle<0.5)
                    JetMatchingEleStatus = 1;
            }
            
            //filling isolation of muon variable
            float absIso_1 = analysisTree.muon_r04_sumChargedHadronPt[indx1];
            float neutralIso_1 = analysisTree.muon_r04_sumNeutralHadronEt[indx1] + analysisTree.muon_r04_sumPhotonEt[indx1] - 0.5*analysisTree.muon_r04_sumPUPt[indx1];
            neutralIso_1 = TMath::Max(float(0),neutralIso_1);
            absIso_1 += neutralIso_1;
            iso_1 = absIso_1/analysisTree.muon_pt[indx1];
            
            float absIso_2 = analysisTree.muon_r04_sumChargedHadronPt[indx2];
            float neutralIso_2 = analysisTree.muon_r04_sumNeutralHadronEt[indx2] + analysisTree.muon_r04_sumPhotonEt[indx2] - 0.5*analysisTree.muon_r04_sumPUPt[indx2];
            neutralIso_2 = TMath::Max(float(0),neutralIso_2);
            absIso_2 += neutralIso_2;
            iso_2 = absIso_2/analysisTree.muon_pt[indx2];
                
            float q1 = analysisTree.muon_charge[indx1];
            float q2 = analysisTree.muon_charge[indx2];
            
            if(q1*q2>0)
                os = false;
            if(q1*q2<0)
                os = true;
                
            pt_1 = analysisTree.muon_pt[indx1];
            eta_1 = analysisTree.muon_eta[indx1];
            phi_1 = analysisTree.muon_phi[indx1];
            
            pt_2 = analysisTree.muon_pt[indx2];
            eta_2 = analysisTree.muon_eta[indx2];
            phi_2 = analysisTree.muon_phi[indx2];
            
            met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
            float dPhiMETMuon = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],analysisTree.pfmetcorr_ex,analysisTree.pfmetcorr_ey);
            mt_1 = TMath::Sqrt(2*met*analysisTree.muon_pt[indx1]*(1-TMath::Cos(dPhiMETMuon)));
            
            ZMMTree->Fill();
            selEventsIsoMuons++;

        }//End event loop
        nFiles++;
        delete _tree;
        file_->Close();
        delete file_;
    }// End of file processing (loop over events in one file)
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



