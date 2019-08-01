#include "NtupleMaker_TauID_all_eras_definitions.h"
#include "NtupleMaker_TauID_all_eras_functions.h"

int main(int argc, char * argv[]) {

   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   
   // ------------------- read config file --------------------
   // **** configuration
   Config cfg(argv[1]);
   
   // general settings
   const string era = cfg.get<string>("Era");
   const bool debug = cfg.get<bool>("Debug");
   const bool isData = cfg.get<bool>("IsData");
   const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
   const string jsonFile = cfg.get<string>("jsonFile");
   
   // trigger information
   const string metHTLName        = cfg.get<string>("MetHLTName");
   const string singleMuonHLTName = cfg.get<string>("SingleMuonHLTName");
   const string singleMuonHLTFilterName = cfg.get<string>("SingleMuonHLTFilterName");
   string singleMuonHLTName1_ = "";
   string singleMuonHLTFilterName1_ = "";
   try{
     singleMuonHLTName1_ = cfg.get<string>("SingleMuonHLTName1");
     singleMuonHLTFilterName1_ = cfg.get<string>("SingleMuonHLTFilterName1");
   }
   catch(...){
     singleMuonHLTName1_ = cfg.get<string>("SingleMuonHLTName");
     singleMuonHLTFilterName1_ = cfg.get<string>("SingleMuonHLTFilterName");
   }
   const string singleMuonHLTName1       = singleMuonHLTName1_;
   const string singleMuonHLTFilterName1 = singleMuonHLTFilterName1_;
   const string pfJet60HLTFilterName = cfg.get<string>("PFJet60HLTFilterName"); 
   const string pfJet80HLTFilterName = cfg.get<string>("PFJet80HLTFilterName"); 
   const string pfJet140HLTFilterName = cfg.get<string>("PFJet140HLTFilterName"); 
   const string pfJet200HLTFilterName = cfg.get<string>("PFJet200HLTFilterName"); 
   const string pfJet260HLTFilterName = cfg.get<string>("PFJet260HLTFilterName"); 
   const string pfJet320HLTFilterName = cfg.get<string>("PFJet320HLTFilterName"); 
   const string pfJet400HLTFilterName = cfg.get<string>("PFJet400HLTFilterName"); 
   const string pfJet450HLTFilterName = cfg.get<string>("PFJet450HLTFilterName"); 
   
   // const string singlePFTau180Trk50Name = cfg.get<string>("SinglePFTau180Trk50Name");
   // const string singlePFTau180Trk50oneprongName = cfg.get<string>("SinglePFTau180Trk50oneprongName");
   
   // TString SinglePFTau180Trk50Name(singlePFTau180Trk50Name);
   // TString SinglePFTau180Trk50oneprongName(singlePFTau180Trk50oneprongName);
   
   TString MetHLTName(metHTLName);
   TString SingleMuonHLTName(singleMuonHLTName);
   TString SingleMuonHLTFilterName(singleMuonHLTFilterName);
   TString SingleMuonHLTName1(singleMuonHLTName1);
   TString SingleMuonHLTFilterName1(singleMuonHLTFilterName1);
   TString PFJet60HLTFilterName(pfJet60HLTFilterName);
   TString PFJet80HLTFilterName(pfJet80HLTFilterName);
   TString PFJet140HLTFilterName(pfJet140HLTFilterName);
   TString PFJet200HLTFilterName(pfJet200HLTFilterName);
   TString PFJet260HLTFilterName(pfJet260HLTFilterName);
   TString PFJet320HLTFilterName(pfJet320HLTFilterName);
   TString PFJet400HLTFilterName(pfJet400HLTFilterName);
   TString PFJet450HLTFilterName(pfJet450HLTFilterName);
   
   // tau cuts
   const float ptTauCut  = cfg.get<float>("PtTauCut");
   const float etaTauCut = cfg.get<float>("EtaTauCut");
   
   // muon selection
   const float ptMuCut       = cfg.get<float>("PtMuCut");
   const float etaMuCut      = cfg.get<float>("EtaMuCut");
   const float isoMuCut      = cfg.get<float>("IsoMuCut");
   const float dxyMuCut      = cfg.get<float>("dxyMuCut");  
   const float dzMuCut       = cfg.get<float>("dzMuCut");
   const float isoSelMuCut   = cfg.get<float>("IsoSelMuCut");
   const float ptSelMuCut    = cfg.get<float>("PtSelMuCut");
   const float ptTrigMuCut   = cfg.get<float>("PtTrigMuCut");
   const bool  isDRIso03 = cfg.get<bool>("IsDRIso03");
   
   // electron selection
   const float ptEleCut   = cfg.get<float>("PtEleCut");
   const float etaEleCut  = cfg.get<float>("EtaEleCut");
   const float isoEleCut  = cfg.get<float>("IsoEleCut");
   const float dxyEleCut  = cfg.get<float>("dxyEleCut");
   const float dzEleCut   = cfg.get<float>("dzEleCut");
   
   // topological cuts (W*->tau+v)
   const float metCut_WTauNu                 = cfg.get<float>("MetCut_WTauNu");
   const float ptTauCut_WTauNu               = cfg.get<float>("PtTauCut_WTauNu");
   
   // topological cuts (W*->mu+v)
   const float ptMuCut_WMuNu               = cfg.get<float>("PtMuCut_WMuNu");
   const float metCut_WMuNu                = cfg.get<float>("MetCut_WMuNu");
   const float etaMuCut_WMuNu              = cfg.get<float>("EtaMuCut_WMuNu");

   // topological cuts (W->muv+Jet)
   const float ptMuCut_WJet              = cfg.get<float>("PtMuCut_WJet");
   const float ptTauCut_WJet             = cfg.get<float>("PtTauCut_WJet");
   const float metCut_WJet               = cfg.get<float>("MetCut_WJet");
   const float etaMuCut_WJet              = cfg.get<float>("EtaMuCut_WJet");

   // topological cuts (dijet)
   const float ptJetCut_DiJet              = cfg.get<float>("PtJetCut_DiJet");
   const float ptTauCut_DiJet              = cfg.get<float>("PtTauCut_DiJet");
   const float etaJetCut_DiJet             = cfg.get<float>("EtaJetCut_DiJet");
   const float deltaPhiTauJetCut_DiJet     = cfg.get<float>("DeltaPhiTauJetCut_DiJet");
   
   // topological cuts (trigger study)
   const float mtCut_Trig = cfg.get<float>("MtCut_Trig");
   
   // trigger eff filename
   const string trigEffFileName = cfg.get<string>("TrigEffFileName");
   
   // trigger weight cuts
   const float mhtNoMu_Trig              = cfg.get<float>("mhtNoMu_Trig");
   const float metNoMu_Trig              = cfg.get<float>("metNoMu_Trig");
   
   // momentum scales
   const float tauMomScale = cfg.get<float>("TauMomScale");
   const string tauDecayMode  = cfg.get<string>("TauDecayMode");
   const float muonMomScale = cfg.get<float>("MuonMomScale");
   const float eleMomScale = cfg.get<float>("EleMomScale");
   const int unclusteredES = cfg.get<int>("UnclusteredES");
   const int jetES = cfg.get<int>("JetES");
   
   const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
   const string MuonTrigFile  = cfg.get<string>("MuonTrigEff");
   
   const string puDataFile = cfg.get<string>("PileUpDataFile");
   const string puMCFile = cfg.get<string>("PileUpMCFile");
   const string samplenameForPUHist = cfg.get<string>("SampleNameForPUHist");
   
   TString PUDataFile(puDataFile);
   TString PUMCFile(puMCFile);
   // **** end of configuration
   
   // ------------------- open files and trees --------------------------------
   std::string rootFileName(argv[2]);
   std::ifstream fileList(argv[2]);
   std::ifstream fileList0(argv[2]);
   std::string ntupleName("makeroottree/AC1B");
   std::string eventHistoName("eventCount/EventCount");
   std::string eventHistoNameData("makeroottree/nEvents");
   std::string weightsHistoName("eventCount/EventWeights");
   
   TString TStrName(rootFileName);
   std::cout <<TStrName <<std::endl;  

   TFile * file = new TFile(TStrName+TString(".root"),"recreate");
   file->cd("");
   
   ntuple_ = new TTree("NTuple","NTuple");
   trigNTuple_ = new TTree("TriggerNTuple","TriggerNTuple");
   
   TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
   TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
   
   TH1D * WMassH     = new TH1D("WMassH","",300,0,3000);
   TH1D * WPtH       = new TH1D("WPtH",  "",100,0,1000);
   TH1D * WDecayH    = new TH1D("WTauDecayH","",5,-1.5,3.5);
   TH1D * WTauDecayH = new TH1D("WDecayH","",11,-1.5,9.5);
   TH1D * dRtauCentralJetH = new TH1D("dRtauCentralJetH","",50,0.,5.0);
   TH1D * dRtauForwardJetH = new TH1D("dRtauForwardJetH","",50,0.,5.0);

   SetupTrees();

  // ------------------- set project directory ----------------------------
  string cmsswBase = (getenv ("CMSSW_BASE"));
  

  // ------------------- initialize good run selection --------------------
  std::vector<Period> periods;
  string fullPathToJsonFile = cmsswBase + "/src/" + jsonFile;
  if (isData) ReadJson(periods,fullPathToJsonFile);
  
  
  // ------------------- initialize PU reweighting ------------------------
  PileUp * PUofficial = new PileUp(); 
  TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/"+PUDataFile,"read");
  TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/"+PUMCFile, "read");
  if (!isData) iniPU(PUofficial,filePUOfficial_data,filePUOfficial_MC, samplenameForPUHist);

  
  // ------------------- initialize SF and efficiencies --------------------
  ScaleFactor * SF_muonIdIso = new ScaleFactor();
  ScaleFactor * SF_muonTrig = new ScaleFactor();
  SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
  SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
  
  map<int,TGraphAsymmErrors*>  map_trigEffData;
  map<int,TGraphAsymmErrors*>  map_trigEffMC;
  TFile * trigEffFile = new TFile(TString(cmsswBase)+"/src/"+trigEffFileName,"read");
  iniTriggerEfficiencies(trigEffFile,map_trigEffData,map_trigEffMC);
  
  // ------------------- set MET filters ---------------------------------
  std::vector<TString> metFlags; metFlags.clear();
  if (era == "2017"){                                      //FIXME: Recommendations changed, has to be updated
     metFlags.push_back("Flag_HBHENoiseFilter");
     metFlags.push_back("Flag_HBHENoiseIsoFilter");
     metFlags.push_back("Flag_globalTightHalo2016Filter");
     metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
     metFlags.push_back("Flag_goodVertices"); 
     if (isData)
        metFlags.push_back("Flag_eeBadScFilter");
     //metFlags.push_back("Flag_muonBadTrackFilter");
     //metFlags.push_back("Flag_chargedHadronTrackResolutionFilter");
     metFlags.push_back("Flag_BadChargedCandidateFilter");
     metFlags.push_back("Flag_BadPFMuonFilter");
     metFlags.push_back("Flag_ecalBadCalibFilter");
  }
  else if (era == "2018"){
     metFlags.push_back("Flag_goodVertices");
     metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
     metFlags.push_back("Flag_HBHENoiseFilter");
     metFlags.push_back("Flag_HBHENoiseIsoFilter");
     metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
     metFlags.push_back("Flag_BadPFMuonFilter");
     metFlags.push_back("Flag_BadChargedCandidateFilter");
     if (isData)
        metFlags.push_back("Flag_eeBadScFilter");
     metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
  }
  else if (era == "2016"){
    metFlags.push_back("Flag_goodVertices");
    if(isData) metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
    metFlags.push_back("Flag_HBHENoiseFilter");
    metFlags.push_back("Flag_HBHENoiseIsoFilter");
    metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
    metFlags.push_back("Flag_BadPFMuonFilter");
  }
  else {
     std::cout << "MET filters not defined for era "<<era<<std::endl;
     exit(-1);
  }
  // ----------------------------------------- ---------------------------------
  // ------------------- open files --------------------------------------------
  // ----------------------------------------- ---------------------------------
  
  // count number of files --->
  std::string dummy;
  while (fileList0 >> dummy) nTotalFiles++;
  
  for (int iF=0; iF<nTotalFiles; ++iF) { // loop over files
     
     std::string filen;
     fileList >> filen;
     
     // opening file
     std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
     TFile * file_ = TFile::Open(TString(filen));
     if (file_->IsZombie()) {
        cout << "Problems opening file : quitting program" << endl;
        exit(-1);
     }
     
     // accessing tree
     TTree * tree_ = (TTree*)file_->Get(TString(ntupleName));
     if (tree_==NULL) { 
        cout << "NTuple " << ntupleName << " is not found in file : quitting program" << endl;
        exit(-1);
     }
     AC1B analysisTree(tree_);
     Long64_t numberOfEntries = analysisTree.GetEntries();
     std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
     
     // ------------------ print era specific information -------------------------
     if( era == "2017") std::cout<<"EE noise cleaning applied!"<<std::endl;
     else std::cout<<"EE noise cleaning NOT applied!"<<std::endl;

     // ----------------------------------------- ---------------------------------
     // ------------------- event loop --------------------------------------------
     // ----------------------------------------- ---------------------------------
     for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
           
        analysisTree.GetEntry(iEntry);
        nEvents++;
        if (nEvents%10000==0) cout << "Processed " << nEvents << endl;
       
        // ------------------- initialize variables --------------------------------------------
        SetDefaultValues();
        
        run_ = analysisTree.event_run;
        lumi_ = analysisTree.event_luminosityblock;
        event_ = analysisTree.event_nr;
        nVert_ = analysisTree.primvertex_count;
        
        if (event_<=0) {
           std::cout << "Event : " << event_ << std::endl;
           std::cout << "From NTuple = " << analysisTree.event_nr << std::endl;
        }
        
        if (debug) {
           cout << "Run                       = " << analysisTree.event_nr           << endl;
           cout << "Event                     = " << analysisTree.event_run          << endl; 
           cout << "Number of gen particles   = " << analysisTree.genparticles_count << endl;
           cout << "Number of taus            = " << analysisTree.gentau_count       << endl;
           cout << "Number of gen taus        = " << analysisTree.tau_count          << endl;
           cout << "Number of jets            = " << analysisTree.pfjet_count        << endl;
           cout << "Number of muons           = " << analysisTree.muon_count         << endl;
           cout << "Number of electrons       = " << analysisTree.electron_count     << endl;
           cout << "Number of tracks          = " << analysisTree.track_count        << endl;
           cout << "Number of PF jets         = " << analysisTree.pfjet_count        << endl;
           cout << "Number of dimuons         = " << analysisTree.dimuon_count       << endl;
           cout << "Number of l1 taus         = " << analysisTree.l1tau_count        << endl;
           cout << "Number of l1 muons        = " << analysisTree.l1muon_count       << endl;
           cout << "Number of l1 egamma       = " << analysisTree.l1egamma_count     << endl;
           cout << "Number of trigger objects = " << analysisTree.l1egamma_count     << endl;
        }	
        
        // Set MC relevant variables
        if (!isData) {
           if (analysisTree.genweight<0) genWeight_ = -1;
           else                          genWeight_ = 1;
           
           weight_ *= genWeight_;
           
           genHt_ = analysisTree.genparticles_lheHt;
           npartons_ = analysisTree.genparticles_noutgoing;
           npartonsNLO_ = analysisTree.genparticles_noutgoing_NLO;
           lheWPt_   = analysisTree.genparticles_lheWPt;
        }
        histWeightsH->Fill(double(0.),double(genWeight_));
      
        
        // **********************************
        // *** Analysis of generator info ***
        // **********************************
        int indexW  = -1;
        int indexMu = -1; // muon from W
        int indexE  = -1; // elec from W
        int indexTau = -1; // tau from W
        vector<TLorentzVector> gentauLV; gentauLV.clear();
        vector<int> gentauDecay; gentauDecay.clear();
        vector<TLorentzVector> genmuonLV; genmuonLV.clear();
        vector<TLorentzVector> genelecLV; genelecLV.clear();
        vector<TLorentzVector> gentaumuonLV; gentaumuonLV.clear();
        vector<TLorentzVector> gentauelecLV; gentauelecLV.clear();
        TLorentzVector wmuonLV; wmuonLV.SetXYZT(0,0,0,0);
        TLorentzVector welecLV; welecLV.SetXYZT(0,0,0,0);
        TLorentzVector wgenvistauLV;  wgenvistauLV.SetXYZT(0,0,0,0);
        TLorentzVector wgentauLV;  wgentauLV.SetXYZT(0,0,0,0);
        TLorentzVector wnuLV;   wnuLV.SetXYZT(0,0,0,0);
        TLorentzVector wallnuLV; wallnuLV.SetXYZT(0,0,0,0);
        
        if (!isData) {
           for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
              
              float pxGen = analysisTree.genparticles_px[igen];
              float pyGen = analysisTree.genparticles_py[igen];
              float pzGen = analysisTree.genparticles_pz[igen];
              float etaGen = PtoEta(pxGen,pyGen,pzGen);
              float ptGen  = PtoPt(pxGen,pyGen);
              
              TLorentzVector genPartLV; genPartLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                          analysisTree.genparticles_py[igen],
                                                          analysisTree.genparticles_pz[igen],
                                                          analysisTree.genparticles_e[igen]);
              
              if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
                 indexW = igen;
              
              if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==12 || 
                  TMath::Abs(analysisTree.genparticles_pdgid[igen])==14 ||
                  TMath::Abs(analysisTree.genparticles_pdgid[igen])==16) { 
                 
                 if (analysisTree.genparticles_info[igen]==(1<<1)) {
                    wnuLV = genPartLV;
                 }
                 if (analysisTree.genparticles_info[igen]==(1<<1) ||
                     analysisTree.genparticles_info[igen]==((1<<1)|(1<<2))) {
                    wallnuLV += genPartLV;
                 }
              }
              
              if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==13) {
                 if ( analysisTree.genparticles_info[igen]==(1<<1) ||
                      analysisTree.genparticles_info[igen]==(1<<0) ) // W/Z->mu
                    genmuonLV.push_back(genPartLV);
                 if ( analysisTree.genparticles_info[igen]==((1<<0)|(1<<2)) ||
                      analysisTree.genparticles_info[igen]==((1<<1)|(1<<2))) // W/Z -> tau -> mu
                    gentaumuonLV.push_back(genPartLV);
                 if ( analysisTree.genparticles_info[igen]==(1<<1) ) { // W -> muv  
                    indexMu = igen;
                    wmuonLV = genPartLV;
                 }
                 //if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->mu
                 //  indexTauMu = igen;
              }
              
              if (TMath::Abs(analysisTree.genparticles_pdgid[igen])==11) { // electron
                 if ( analysisTree.genparticles_info[igen]==(1<<1) || 
                      analysisTree.genparticles_info[igen]==(1<<2) ) // W/Z->e 
                    genelecLV.push_back(genPartLV);
                 if ( analysisTree.genparticles_info[igen]==((1<<0)|(1<<2)) || 
                      analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W/Z -> tau -> e 
                    gentauelecLV.push_back(genPartLV);
                 if ( analysisTree.genparticles_info[igen]==(1<<1) ) { // W->ev
                    indexE = igen;
                    welecLV = genPartLV;
                 }
                 //if ( analysisTree.genparticles_info[igen]==((1<<1)|(1<<2)) ) // W->tau->e
                 //  indexTauE = igen;
              }
           }
           
           for (unsigned int igentau=0; igentau<analysisTree.gentau_count;++igentau) {
              TLorentzVector GenVisTau; GenVisTau.SetXYZT(analysisTree.gentau_visible_px[igentau],
                                                          analysisTree.gentau_visible_py[igentau],
                                                          analysisTree.gentau_visible_pz[igentau],
                                                          analysisTree.gentau_visible_e[igentau]);
              TLorentzVector GenTau; GenTau.SetXYZT(analysisTree.gentau_px[igentau],
                                                    analysisTree.gentau_py[igentau],
                                                    analysisTree.gentau_pz[igentau],
                                                    analysisTree.gentau_e[igentau]);
              
              if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau] ) { // W/Z->tau  
                 gentauLV.push_back(GenVisTau);
                 gentauDecay.push_back(analysisTree.gentau_decayMode[igentau]);
                 indexTau = igentau;
                 wgenvistauLV = GenVisTau;
                 wgentauLV = GenTau;
              }
           }
        }
        
        TLorentzVector lorentzVectorGenW; lorentzVectorGenW.SetXYZT(0,0,0,0);
        if (indexW>=0) 
           lorentzVectorGenW.SetXYZT(analysisTree.genparticles_px[indexW],
                                     analysisTree.genparticles_py[indexW],
                                     analysisTree.genparticles_pz[indexW],
                                     analysisTree.genparticles_e[indexW]);
        
        
        
        if (indexW>=0) {
           wMass_ = lorentzVectorGenW.M();
           wPt_   = lorentzVectorGenW.Pt();
           wEta_  = lorentzVectorGenW.Eta();
           wPhi_  = lorentzVectorGenW.Phi();
           nuWPt_ = wnuLV.Pt();
           nuWEta_ = wnuLV.Eta();
           nuWPhi_ = wnuLV.Phi();
           wDecay_ = 0;
           wTauDecay_ = -1;
           wCharge_ = double(analysisTree.genparticles_pdgid[indexW])/TMath::Abs(double(analysisTree.genparticles_pdgid[indexW]));
           if (indexMu>=0) {
              wDecay_ = 2;
              lepWPt_  = wmuonLV.Pt();
              lepWEta_ = wmuonLV.Eta(); 
              lepWPhi_ = wmuonLV.Phi();
              lepWE_   = wmuonLV.E();
           }
           else if (indexE>=0) {
              wDecay_ = 1;
              lepWPt_  = welecLV.Pt();
              lepWEta_ = welecLV.Eta();
              lepWPhi_ = welecLV.Phi();
              lepWE_   = welecLV.E();
           }
           else if (indexTau>=0) {
              lepWPt_  = wgenvistauLV.Pt();
              lepWEta_ = wgenvistauLV.Eta();
              lepWPhi_ = wgenvistauLV.Phi();
              lepWE_   = wgenvistauLV.E();
              
              genTauWPt_  = wgentauLV.Pt();
              genTauWEta_ = wgentauLV.Eta();
              genTauWPhi_ = wgentauLV.Phi();
              genTauWE_   = wgentauLV.E();
              wDecay_ = 3;
              wTauDecay_ = analysisTree.gentau_decayMode[indexTau];
              if (wTauDecay_<0) wTauDecay_ = -1;
           }
           WMassH->Fill(wMass_,weight_);
           WPtH->Fill(wPt_,weight_);
           WDecayH->Fill(wDecay_,weight_);
           WTauDecayH->Fill(wTauDecay_,weight_);
        }
        
        // *****************************************
        // ***** primary vertex selection, *********
        // ***** PU weights, good run selection ****
        // *****************************************
        
        if (analysisTree.primvertex_count==0) continue; // at least one good primary vertex
        
        if (!isData) { 
           puWeight_ =  float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
           weight_ *= puWeight_; 
           if (debug)
              cout << "nPU = " << analysisTree.numtruepileupinteractions << " --> puweight = " << puWeight_ << endl;
        }
        
        if (isData && applyGoodRunSelection){
           int n=analysisTree.event_run;
           int lum = analysisTree.event_luminosityblock;      
           if (!GoodRunSelection(n,lum,periods)) continue;
        }
        
        // *********************************
        // ***** accessing trigger info ****
        // *********************************
        
        bool isMetHLT = false;
        for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
           TString trigName(it->first);
           if (trigName.Contains(MetHLTName)) {
              if (it->second==1) isMetHLT = true;
           }
        }
        trigger_ = isMetHLT;
        trig_ = isMetHLT;
        
        isSingleMuonHLTFilter1= AccessTriggerInfo(analysisTree,SingleMuonHLTFilterName,nSingleMuonHLTFilter);
        isSingleMuonHLTFilter2= AccessTriggerInfo(analysisTree,SingleMuonHLTFilterName1,nSingleMuonHLTFilter1);
	isSingleMuonHLTFilter = isSingleMuonHLTFilter1 || isSingleMuonHLTFilter2;
        if (!isSingleMuonHLTFilter) {
	  std::cout << "Single Muon HLT filter not found" << std::endl;
	  exit(-1);
	}
	isPFJet60HLTFilter  = AccessTriggerInfo(analysisTree,PFJet60HLTFilterName,nPFJet60HLTFilter);
        isPFJet80HLTFilter  = AccessTriggerInfo(analysisTree,PFJet80HLTFilterName,nPFJet80HLTFilter);
        isPFJet140HLTFilter = AccessTriggerInfo(analysisTree,PFJet140HLTFilterName,nPFJet140HLTFilter);
        isPFJet200HLTFilter = AccessTriggerInfo(analysisTree,PFJet200HLTFilterName,nPFJet200HLTFilter);
        isPFJet260HLTFilter = AccessTriggerInfo(analysisTree,PFJet260HLTFilterName,nPFJet260HLTFilter);
        isPFJet320HLTFilter = AccessTriggerInfo(analysisTree,PFJet320HLTFilterName,nPFJet320HLTFilter);
        isPFJet400HLTFilter = AccessTriggerInfo(analysisTree,PFJet400HLTFilterName,nPFJet400HLTFilter);
        isPFJet450HLTFilter = AccessTriggerInfo(analysisTree,PFJet450HLTFilterName,nPFJet450HLTFilter);
	if (!isPFJet60HLTFilter || !isPFJet80HLTFilter || !isPFJet140HLTFilter || !isPFJet200HLTFilter || !isPFJet260HLTFilter || !isPFJet320HLTFilter || !isPFJet400HLTFilter || !isPFJet450HLTFilter ) {
	  std::cout << "PFJet HLT filter not found" << std::endl;
	  exit(-1);
	}
        //AccessTriggerInfo(analysisTree,SinglePFTau180Trk50Name,nSinglePFTau180Trk50Filter,isSinglePFTau180Trk50Filter);
        //AccessTriggerInfo(analysisTree,SinglePFTau180Trk50oneprongName,nSinglePFTau180Trk50oneprongFilter,isSinglePFTau180Trk50oneprongFilter);
        
        
        // ***************************************************
        // accessing PF MET and changing momentum scale of met
        // ***************************************************
        float pfmetcorr_ex = analysisTree.pfmetcorr_ex;
        float pfmetcorr_ey = analysisTree.pfmetcorr_ey;

        if (!isData) { 
           if (jetES<0) {
              pfmetcorr_ex = analysisTree.pfmetcorr_ex_JetEnDown;
              pfmetcorr_ey = analysisTree.pfmetcorr_ey_JetEnDown;
           }
           else if (jetES>0) {
              pfmetcorr_ex = analysisTree.pfmetcorr_ex_JetEnUp;
              pfmetcorr_ey = analysisTree.pfmetcorr_ey_JetEnUp;
           }
           else if (unclusteredES<0) {
              pfmetcorr_ex = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
              pfmetcorr_ey = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
           }
           else if (unclusteredES>0) {
              pfmetcorr_ex = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
              pfmetcorr_ey = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
           }
           else {
              pfmetcorr_ex = analysisTree.pfmetcorr_ex;
              pfmetcorr_ey = analysisTree.pfmetcorr_ey;
           }
        }
        
        met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
        if (met_<1e-4) met_ = 1e-4;
        metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
        TLorentzVector lorentzVectorMet; lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);
        
      
        // *************************
        // **** accessing muons ****
        // *************************
        TLorentzVector lorentzVectorAllMuons; lorentzVectorAllMuons.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorAllSelMuons; lorentzVectorAllSelMuons.SetXYZT(0,0,0,0);
        std::vector<unsigned int> muonIndexes; muonIndexes.clear();
        std::vector<unsigned int> selMuonIndexes; selMuonIndexes.clear();
        std::vector<TLorentzVector>  lorentzVectorMuons; lorentzVectorMuons.clear();
        std::vector<TLorentzVector>  lorentzVectorSelMuons; lorentzVectorSelMuons.clear();
        int indexTriggerMu = -1;
        float ptTriggerMu  = -1;
        float etaTriggerMu = -1; 
        float muonHt = 0;
        
        // ----------------------- store muons that pass selection ----------------------- 
        for (unsigned int imuon=0; imuon<analysisTree.muon_count; ++imuon) {
           analysisTree.muon_px[imuon] *= muonMomScale;
           analysisTree.muon_py[imuon] *= muonMomScale;
           analysisTree.muon_pz[imuon] *= muonMomScale;
           analysisTree.muon_pt[imuon] *= muonMomScale;
           
           float relIso=10;
         // muon selection
           if (!PassesMuonSelection(analysisTree,imuon,ptMuCut,etaMuCut,dxyMuCut,dzMuCut,isDRIso03,isoMuCut, relIso, era)) continue;
         muonIndexes.push_back(imuon);  
         
         TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[imuon],
                                               analysisTree.muon_py[imuon],
                                               analysisTree.muon_pz[imuon],
                                               muonMass);
         lorentzVectorAllMuons += muonLV;
         lorentzVectorMuons.push_back(muonLV);
         muonHt += analysisTree.muon_pt[imuon];
         
         // selected muons -->
         bool passedSelIso = relIso < isoSelMuCut;
         if (analysisTree.muon_pt[imuon]>ptSelMuCut && passedSelIso) {
            selMuonIndexes.push_back(imuon);
            TLorentzVector muonSelLV; muonSelLV.SetXYZM(analysisTree.muon_px[imuon],
                                                        analysisTree.muon_py[imuon],
                                                        analysisTree.muon_pz[imuon],
                                                        muonMass);
            lorentzVectorAllSelMuons += muonSelLV;
            lorentzVectorSelMuons.push_back(muonSelLV);
         }
         // triggering muon -->
         if (analysisTree.muon_pt[imuon]>ptTrigMuCut && 
             passedSelIso) {
            bool trigMatch1 = TriggerMatching(analysisTree, analysisTree.muon_eta[imuon],analysisTree.muon_phi[imuon],nSingleMuonHLTFilter);
            bool trigMatch2 = TriggerMatching(analysisTree, analysisTree.muon_eta[imuon],analysisTree.muon_phi[imuon],nSingleMuonHLTFilter1);
	    bool trigMatch = trigMatch1 || trigMatch2;
            if (trigMatch&&analysisTree.muon_pt[imuon]>ptTriggerMu) {
               ptTriggerMu = analysisTree.muon_pt[imuon];
               etaTriggerMu = analysisTree.muon_eta[imuon];
               indexTriggerMu = int(imuon);
            }
         }
        }
        metNoMu_    =  (lorentzVectorMet+lorentzVectorAllMuons).Pt();
        metNoSelMu_ =  (lorentzVectorMet+lorentzVectorAllSelMuons).Pt();
        
        nMuon_ = muonIndexes.size();
        nMuonTrig_ = nMuon_;
        nSelMuon_ = selMuonIndexes.size();
        nSelMuonTrig_ = nSelMuon_;
        
        float ptSecondMu  = -1;
        int indexSecondMu = -1;
        TLorentzVector lorentzVectorTriggerMu; lorentzVectorTriggerMu.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorSecondMu;  lorentzVectorSecondMu.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorZ; lorentzVectorZ.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorW; lorentzVectorW.SetXYZT(0,0,0,0);
       
        // ----------------------- set properties of muons ----------------------- 
        if (indexTriggerMu>=0) {
           lorentzVectorTriggerMu.SetXYZM(analysisTree.muon_px[indexTriggerMu],
                                          analysisTree.muon_py[indexTriggerMu],
                                          analysisTree.muon_pz[indexTriggerMu],
                                          muonMass);
           muonPt_  = lorentzVectorTriggerMu.Pt();
           muonPz_  = lorentzVectorTriggerMu.Pz();
           muonEta_ = lorentzVectorTriggerMu.Eta();
           muonPhi_ = lorentzVectorTriggerMu.Phi();
           muonQ_   = int(analysisTree.muon_charge[indexTriggerMu]);
           pfmetcorr_ex = pfmetcorr_ex + lorentzVectorTriggerMu.Px() - lorentzVectorTriggerMu.Px(); 
           pfmetcorr_ey = pfmetcorr_ey + lorentzVectorTriggerMu.Py() - lorentzVectorTriggerMu.Py();
           met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
           metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
           lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);
           mtmuon_  = mT(lorentzVectorTriggerMu,lorentzVectorMet);
           dPhiMetMuon_ = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet); 
           lorentzVectorW = lorentzVectorTriggerMu + lorentzVectorMet;
           //select second muon
           for (unsigned int iMu = 0; iMu < selMuonIndexes.size(); ++iMu) {
              int indexMu = int(selMuonIndexes.at(iMu));
              if (indexMu==indexTriggerMu) continue;
              float netcharge = analysisTree.muon_charge[indexTriggerMu]*analysisTree.muon_charge[indexMu];
              if (netcharge>0) continue;
              if (analysisTree.muon_pt[indexMu]>ptSecondMu) {
                 ptSecondMu = analysisTree.muon_pt[indexMu];
                 //etaSecondMu = analysisTree.muon_eta[indexMu];
                 indexSecondMu = int(indexMu);
              }
           }
           if (indexSecondMu>=0) {
              lorentzVectorSecondMu.SetXYZM(analysisTree.muon_px[indexSecondMu],
                                            analysisTree.muon_py[indexSecondMu],
                                            analysisTree.muon_pz[indexSecondMu],
                                            muonMass);
              lorentzVectorZ = lorentzVectorTriggerMu + lorentzVectorSecondMu;
              muon2Pt_  = lorentzVectorSecondMu.Pt();
              muon2Eta_ = lorentzVectorSecondMu.Eta();
              muon2Phi_ = lorentzVectorSecondMu.Phi();
              muon2Q_   = int(analysisTree.muon_charge[indexSecondMu]);
           }
           isWTrig_ = mtmuon_ >  mtCut_Trig;
        }
        float ptLeadingMu = ptTriggerMu;
        float ptTrailingMu = ptSecondMu;
        
        if (ptTrailingMu>ptLeadingMu) {
           ptLeadingMu = ptSecondMu;
           ptTrailingMu = ptTriggerMu;
        }
     
        if (debug) std::cout << "end accessing muons " << std::endl;
        
     
        // *****************************
        // **** accessing electrons ****
        // *****************************
        
        TLorentzVector lorentzVectorAllElectrons; lorentzVectorAllElectrons.SetXYZT(0,0,0,0);
        std::vector<unsigned int> eleIndexes; eleIndexes.clear();
        std::vector<TLorentzVector>  lorentzVectorElectrons; lorentzVectorElectrons.clear();
        float elecHt = 0;
        for (unsigned int ielec=0; ielec<analysisTree.electron_count; ++ielec) {
           analysisTree.electron_px[ielec] *= eleMomScale;
           analysisTree.electron_py[ielec] *= eleMomScale;
           analysisTree.electron_pz[ielec] *= eleMomScale;
           analysisTree.electron_pt[ielec] *= eleMomScale;
           
           if (PassesElectronSelection(analysisTree, ielec, dxyEleCut, dzEleCut, isoEleCut, ptEleCut, etaEleCut, era) ){
              TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[ielec],
                                                            analysisTree.electron_py[ielec],
                                                            analysisTree.electron_pz[ielec],
                                                            electronMass);
              lorentzVectorAllElectrons += electronLV; 
              eleIndexes.push_back(ielec);
              elecHt = analysisTree.electron_pt[ielec];
           }
        }
        nElec_ = eleIndexes.size();
        
        if (debug)	std::cout << "end accessing electrons " << std::endl;
        
        
        // ************************
        // **** accessing jets ****
        // ************************
        
        TLorentzVector lorentzVectorAllJetsForMht; lorentzVectorAllJetsForMht.SetXYZT(0,0,0,0);
        std::vector<unsigned int> centralJets20Indexes; centralJets20Indexes.clear();
        std::vector<unsigned int> forwardJets20Indexes; forwardJets20Indexes.clear();
        std::vector<unsigned int> centralJets30Indexes; centralJets30Indexes.clear();
        std::vector<unsigned int> forwardJets30Indexes; forwardJets30Indexes.clear();
        float htCentral20 = 0;
        float htCentral30 = 0;
        float htForward20 = 0;
        float htForward30 = 0;
        std::vector<unsigned int> triggerJetsIndexes; triggerJetsIndexes.clear();
        std::vector<bool> jets60trigger; jets60trigger.clear();
        std::vector<bool> jets80trigger; jets80trigger.clear();
        std::vector<bool> jets140trigger; jets140trigger.clear();
        std::vector<bool> jets200trigger; jets200trigger.clear();
        std::vector<bool> jets260trigger; jets260trigger.clear();
        std::vector<bool> jets320trigger; jets320trigger.clear();
        std::vector<bool> jets400trigger; jets400trigger.clear();
        std::vector<bool> jets450trigger; jets450trigger.clear();
        
        for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
           
           // ----------------------- clean for EEnoise jets ----------------------- 
           if( era == "2017" && analysisTree.pfjet_pt[ijet] < 50 && fabs(analysisTree.pfjet_eta[ijet]) < 3.139 && fabs(analysisTree.pfjet_eta[ijet]) > 2.65) continue;
           
           // Scale for sys uncertainties
           float scaleJ = 1;
           
           if (jetES<0)      scaleJ = 1.0 - analysisTree.pfjet_jecUncertainty[ijet];
           else if (jetES>0) scaleJ = 1.0 + analysisTree.pfjet_jecUncertainty[ijet];
           else 	          scaleJ = 1.0;
           
           analysisTree.pfjet_px[ijet] *= scaleJ;
           analysisTree.pfjet_py[ijet] *= scaleJ;
           analysisTree.pfjet_pz[ijet] *= scaleJ;
           analysisTree.pfjet_pt[ijet] *= scaleJ;
           analysisTree.pfjet_e[ijet]  *= scaleJ;
           
           float absJetEta = fabs(analysisTree.pfjet_eta[ijet]);
           if (absJetEta>5.2) continue;
           if (analysisTree.pfjet_pt[ijet]<20.0) continue; 
           
           // ----------------------- jet ID ----------------------- 
	   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
           bool isPFJetId = false;
           if (era == "2017") isPFJetId = tightJetiD_2017(analysisTree,int(ijet));
           else if (era == "2018") isPFJetId = tightJetiD_2018(analysisTree,int(ijet));
	   else if (era == "2016") isPFJetId = looseJetiD_2016(analysisTree,int(ijet)); // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
           else {
              std::cout<<"No Jet Id specified for era "<<era<<std::endl;
              exit(-1);
           }
           if (!isPFJetId) continue;
           
           // jet four-vector
           TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[ijet],
                                               analysisTree.pfjet_py[ijet],
                                               analysisTree.pfjet_pz[ijet],
                                               analysisTree.pfjet_e[ijet]);
           
           // counting jets for Mht
           lorentzVectorAllJetsForMht += jetLV;
           
           if (OverlapWithLepton(analysisTree, ijet, muonIndexes, eleIndexes)) continue;
           
           // pt > 20 GeV
           if (analysisTree.pfjet_pt[ijet]>20.0) {
              if (absJetEta<2.4) {
                 centralJets20Indexes.push_back(ijet);
                 htCentral20 += analysisTree.pfjet_pt[ijet];
              }
              else if (absJetEta<4.7) {
                 forwardJets20Indexes.push_back(ijet);
                 htForward20 += analysisTree.pfjet_pt[ijet]; 
              }
              if (nJets20_<10) {
                 jet20Pt_[nJets20_]  = analysisTree.pfjet_pt[ijet];
                 jet20Eta_[nJets20_] = analysisTree.pfjet_eta[ijet];
                 jet20Phi_[nJets20_] = analysisTree.pfjet_phi[ijet];
                 nJets20_++;
              }
           }
           
           // pt > 30 GeV
           if (analysisTree.pfjet_pt[ijet]>30) {
              if (absJetEta<2.4) {
                 centralJets30Indexes.push_back(ijet);
                 htCentral30 += analysisTree.pfjet_pt[ijet];
              }
              else if (absJetEta<4.7) {
                 forwardJets30Indexes.push_back(ijet);
                 htForward30 += analysisTree.pfjet_pt[ijet];
              }
           }
         
           // triggering jets
           if (analysisTree.pfjet_pt[ijet]<ptJetCut_DiJet) continue;
           if (absJetEta>etaJetCut_DiJet) continue;
           bool trigMatch60 = false;
           bool trigMatch80 = false;
           bool trigMatch140 = false;
           bool trigMatch200 = false;
           bool trigMatch260 = false;
           bool trigMatch320 = false;
           bool trigMatch400 = false;
           bool trigMatch450 = false;
           for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
              float dRtrig = deltaR(analysisTree.pfjet_eta[ijet],analysisTree.pfjet_phi[ijet],
                                    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
              if (dRtrig>0.5) continue;
              
              if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   trigMatch60 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   trigMatch80 = true; 
              if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) trigMatch140 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) trigMatch200 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet260HLTFilter]&&isPFJet260HLTFilter)) trigMatch260 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet320HLTFilter]&&isPFJet320HLTFilter)) trigMatch320 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet400HLTFilter]&&isPFJet400HLTFilter)) trigMatch400 = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet450HLTFilter]&&isPFJet450HLTFilter)) trigMatch450 = true;
           }
           if (!isData) {
              trigMatch60 = true;
              trigMatch80 = true;
              trigMatch140 = true;
              trigMatch200 = true;
              trigMatch260 = true;
              trigMatch320 = true;
              trigMatch400 = true;
              trigMatch450 = true;
         }
           triggerJetsIndexes.push_back(ijet);
           jets60trigger.push_back(trigMatch60);
           jets80trigger.push_back(trigMatch80);
           jets140trigger.push_back(trigMatch140);
           jets200trigger.push_back(trigMatch200);
           jets260trigger.push_back(trigMatch260);
           jets320trigger.push_back(trigMatch320);
           jets400trigger.push_back(trigMatch400);
           jets450trigger.push_back(trigMatch450);
      }
        nJetsCentral20_ = centralJets20Indexes.size();
        nJetsCentral30_ = centralJets30Indexes.size();
        nJetsForward20_ = forwardJets20Indexes.size();
        nJetsForward30_ = forwardJets30Indexes.size();
        JetHt_     = htCentral30 + htForward30;
        SoftJetHt_ = htCentral20 + htForward30;
        Ht_        = JetHt_     + muonHt + elecHt;
        SoftHt_    = SoftJetHt_ + muonHt + elecHt;
        mht_        = (lorentzVectorAllJetsForMht).Pt();
        mhtNoMu_    = (lorentzVectorAllJetsForMht - lorentzVectorAllMuons).Pt();
        mhtNoSelMu_ = (lorentzVectorAllJetsForMht - lorentzVectorAllSelMuons).Pt();
        TLorentzVector lorentzVectorJet; lorentzVectorJet.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorJet2; lorentzVectorJet2.SetXYZT(0,0,0,0);
        if (nJetsCentral30_>0) {
           unsigned int indexJet0 = centralJets30Indexes.at(0);
           jetPt_ = analysisTree.pfjet_pt[indexJet0];
           jetEta_ = analysisTree.pfjet_eta[indexJet0];
           jetPhi_ = analysisTree.pfjet_phi[indexJet0];
           jetBtag_ = analysisTree.pfjet_btag[indexJet0][0];
           jetChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet0];
           jetNeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet0];
           jetChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet0];
           jetNeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];
           jetNeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet0]/analysisTree.pfjet_e[indexJet0];  
           pfJet60_ = false;
           pfJet80_ = false;
           pfJet140_ = false;
           pfJet200_ = false;
           pfJet260_ = false;
           pfJet320_ = false;
           pfJet400_ = false;
           pfJet450_ = false;
           lorentzVectorJet.SetXYZT(analysisTree.pfjet_px[indexJet0],
                                    analysisTree.pfjet_py[indexJet0],
                                    analysisTree.pfjet_pz[indexJet0],
                                    analysisTree.pfjet_e[indexJet0]);
           for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
              float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet0],analysisTree.pfjet_phi[indexJet0],
                                    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
              if (dRtrig>0.5) continue;
              if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pfJet60_  = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pfJet80_  = true; 
              if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pfJet140_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) pfJet200_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet260HLTFilter]&&isPFJet260HLTFilter)) pfJet260_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet320HLTFilter]&&isPFJet320HLTFilter)) pfJet320_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet400HLTFilter]&&isPFJet400HLTFilter)) pfJet400_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet450HLTFilter]&&isPFJet450HLTFilter)) pfJet450_ = true;
           }
        }
        if (nJetsCentral30_>1) {
           unsigned int indexJet1 = centralJets30Indexes.at(1);
           jet2Pt_ = analysisTree.pfjet_pt[indexJet1];
           jet2Eta_ = analysisTree.pfjet_eta[indexJet1];
           jet2Phi_ = analysisTree.pfjet_phi[indexJet1];
           jet2Btag_ = analysisTree.pfjet_btag[indexJet1][0];
           jet2ChargedMult_ = analysisTree.pfjet_chargedmulti[indexJet1];
           jet2NeutralMult_ = analysisTree.pfjet_neutralmulti[indexJet1];
           jet2ChargedHadMult_ = analysisTree.pfjet_chargedhadronmulti[indexJet1];
           jet2NeutralEMEnergyFraction_  = analysisTree.pfjet_neutralemenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];
           jet2NeutralHadEnergyFraction_ = analysisTree.pfjet_neutralhadronicenergy[indexJet1]/analysisTree.pfjet_e[indexJet1];  
           pf2Jet60_ = false;
           pf2Jet80_ = false;
           pf2Jet140_ = false;
           pf2Jet200_ = false;
           lorentzVectorJet2.SetXYZT(analysisTree.pfjet_px[indexJet1],
                                     analysisTree.pfjet_py[indexJet1],
                                     analysisTree.pfjet_pz[indexJet1],
                                     analysisTree.pfjet_e[indexJet1]);
           for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
              float dRtrig = deltaR(analysisTree.pfjet_eta[indexJet1],analysisTree.pfjet_phi[indexJet1],
                                    analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
              if (dRtrig>0.5) continue;
              if ((analysisTree.trigobject_filters[iT][nPFJet60HLTFilter]&&isPFJet60HLTFilter))   pf2Jet60_  = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet80HLTFilter]&&isPFJet80HLTFilter))   pf2Jet80_  = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet140HLTFilter]&&isPFJet140HLTFilter)) pf2Jet140_ = true;
              if ((analysisTree.trigobject_filters[iT][nPFJet200HLTFilter]&&isPFJet200HLTFilter)) pf2Jet200_ = true;
           }
        }
        
        if (debug) std::cout << "end accessing jets" << std::endl;
        
        // ************************
        // **** accessing taus ****
        // ************************
        std::vector<unsigned int> tauIndexes; tauIndexes.clear();
        std::vector<unsigned int> tau20Indexes; tau20Indexes.clear();
        std::vector<unsigned int> tau30Indexes; tau30Indexes.clear();
        std::vector<int> tauGenMatchDecay; tauGenMatchDecay.clear();
        
        for (unsigned int itau=0; itau<analysisTree.tau_count; ++itau) { // loop over taus
           
           if( (tauDecayMode == "1prong0pizeros"     && analysisTree.tau_decayMode[itau]==0) ||
               (tauDecayMode == "1prongUpTo4pizeros" && analysisTree.tau_decayMode[itau]>=1 && analysisTree.tau_decayMode[itau]<=4) ||
               (tauDecayMode == "3prong0pizeros"     && (analysisTree.tau_decayMode[itau]==10 || analysisTree.tau_decayMode[itau]==11)) ){
              
              analysisTree.tau_px[itau]   *= tauMomScale;
              analysisTree.tau_py[itau]   *= tauMomScale;
              analysisTree.tau_pz[itau]   *= tauMomScale;
              analysisTree.tau_pt[itau]   *= tauMomScale;
              analysisTree.tau_e[itau]    *= tauMomScale;
              analysisTree.tau_mass[itau] *= tauMomScale;
           }
           if (!LooseTauSelection(analysisTree,itau)) continue;
           // finding matching mu and e-->
           int matchedMuIndex = MatchingMuon(analysisTree, itau, muonIndexes);
           int matchedEleIndex = MatchingElectron(analysisTree, itau, eleIndexes);
           
           if (analysisTree.tau_pt[itau]>20.&&matchedMuIndex<0&&matchedEleIndex<0) {
              tau20Indexes.push_back(itau);
           }
           
           if (analysisTree.tau_pt[itau]>30.&&matchedMuIndex<0&&matchedEleIndex<0) {
              tau30Indexes.push_back(itau);
           }
           
           if (analysisTree.tau_pt[itau]>ptTauCut&&fabs(analysisTree.tau_eta[itau])<etaTauCut&&matchedMuIndex<0&&matchedEleIndex<0) { 
              tauIndexes.push_back(itau);
              int genMatchDecay = -1;
              float dRmin = 0.2;
              for (unsigned int igentau=0; igentau<gentauLV.size(); ++igentau) {
                 TLorentzVector genTauLV = gentauLV.at(igentau);
                 int decaymode = gentauDecay.at(igentau);
                 float dR = deltaR(analysisTree.tau_eta[itau],analysisTree.tau_phi[itau],
                                   genTauLV.Eta(),genTauLV.Phi());
                 if (dR<dRmin) {
                    dRmin = dR;
                    genMatchDecay = decaymode;
                 }
              }
              tauGenMatchDecay.push_back(genMatchDecay);
           }
        } // end loop over taus
        
        nTaus20_ = tau20Indexes.size();
        nTaus30_ = tau30Indexes.size();
        nSelTaus_ = tauIndexes.size();
        TLorentzVector lorentzVectorTau; lorentzVectorTau.SetXYZT(0,0,0,0);
        TLorentzVector lorentzVectorTauJet; lorentzVectorTauJet.SetXYZT(0,0,0,0);
        
        // -----------------------  set tau properties ----------------------- 
        if (nSelTaus_>0) {
           
           unsigned int indexTau = tauIndexes.at(0);
           lorentzVectorTau.SetXYZM(analysisTree.tau_px[indexTau],
                                    analysisTree.tau_py[indexTau],
                                    analysisTree.tau_pz[indexTau],
                                    analysisTree.tau_mass[indexTau]);
           
           pfmetcorr_ex = pfmetcorr_ex + lorentzVectorTau.Px() - lorentzVectorTau.Px();
           pfmetcorr_ey = pfmetcorr_ey + lorentzVectorTau.Py() - lorentzVectorTau.Py();
           met_ = TMath::Sqrt(pfmetcorr_ex*pfmetcorr_ex+pfmetcorr_ey*pfmetcorr_ey);
           metphi_ = TMath::ATan2(pfmetcorr_ey,pfmetcorr_ex);
           lorentzVectorMet.SetXYZT(pfmetcorr_ex,pfmetcorr_ey,0,met_);
           mttau_ = mT(lorentzVectorTau,lorentzVectorMet);
           mtgen_ = mT(wgentauLV,wnuLV);
           tauPt_ = analysisTree.tau_pt[indexTau];
           tauPz_ = analysisTree.tau_pz[indexTau];
           tauEta_ = analysisTree.tau_eta[indexTau];
           tauPhi_ = analysisTree.tau_phi[indexTau];
           tauMass_ = analysisTree.tau_mass[indexTau];
           tauQ_ = int(analysisTree.tau_charge[indexTau]);
           tauNtrk1_ = analysisTree.tau_ntracks_pt1[indexTau];
           tauNtrk05_ = analysisTree.tau_ntracks_pt05[indexTau];
           tauNtrk08_ = analysisTree.tau_ntracks_pt08[indexTau];
           
           tauLeadingTrackPt_ = PtoPt(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                      analysisTree.tau_leadchargedhadrcand_py[indexTau]);
           
           tauLeadingTrackEta_ = PtoEta(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                        analysisTree.tau_leadchargedhadrcand_py[indexTau],
                                        analysisTree.tau_leadchargedhadrcand_pz[indexTau]);
           
           tauLeadingTrackPhi_ = PtoPhi(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                        analysisTree.tau_leadchargedhadrcand_py[indexTau]);
           
           tauLeadingTrackDz_  = analysisTree.tau_leadchargedhadrcand_dz[indexTau];
           tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];
         
           tauDecay_ = analysisTree.tau_decayMode[indexTau];
           tauGenDecay_ = analysisTree.tau_genDecayMode[indexTau];
           tauGenMatchDecay_ = tauGenMatchDecay.at(0);
           
           if (tauDecay_<0) tauDecay_ = -1;
           if (tauGenDecay_<0) tauGenDecay_ = -1;
           if (tauGenMatchDecay_<0) tauGenMatchDecay_ = -1;
           
           tauGenMatch_ = 6;
           if (tauGenMatchDecay_>=0) tauGenMatch_ = 5;
           if (!isData) tauGenMatch_ = FindTauGenMatch(analysisTree, tauEta_, tauPhi_);
           
           tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
           tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;
           
           tauIso_ = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[indexTau];
           tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
           tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
           tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
           
           tauVLooseMvaIso_ = analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           tauVVTightMvaIso_ = analysisTree.tau_byVVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
           
           tauVVLooseMva2017v2Iso_ = analysisTree.tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauVLooseMva2017v2Iso_ = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauLooseMva2017v2Iso_ = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauMediumMva2017v2Iso_ = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauTightMva2017v2Iso_ = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauVTightMva2017v2Iso_ = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           tauVVTightMva2017v2Iso_ = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
           
           tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
           tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;
           
           tauAntiElectronVLooseMVA6_ = analysisTree.tau_againstElectronVLooseMVA6[indexTau] > 0.5;
           tauAntiElectronLooseMVA6_  = analysisTree.tau_againstElectronLooseMVA6[indexTau] > 0.5;
           tauAntiElectronMediumMVA6_  = analysisTree.tau_againstElectronMediumMVA6[indexTau] > 0.5;
           tauAntiElectronTightMVA6_  = analysisTree.tau_againstElectronTightMVA6[indexTau] > 0.5;
           tauAntiElectronVTightMVA6_ = analysisTree.tau_againstElectronVTightMVA6[indexTau] > 0.5;
           
           taubyDeepTau2017v2VSeraw_ = analysisTree.tau_byDeepTau2017v2VSeraw[indexTau] > 0.5;
           taubyDeepTau2017v2VSjetraw_ = analysisTree.tau_byDeepTau2017v2VSjetraw[indexTau] > 0.5;
           taubyDeepTau2017v2VSmuraw_ = analysisTree.tau_byDeepTau2017v2VSmuraw[indexTau] > 0.5;
           taubyLooseDeepTau2017v2VSe_ = analysisTree.tau_byLooseDeepTau2017v2VSe[indexTau] > 0.5;
           taubyLooseDeepTau2017v2VSjet_ = analysisTree.tau_byLooseDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyLooseDeepTau2017v2VSmu_ = analysisTree.tau_byLooseDeepTau2017v2VSmu[indexTau] > 0.5;
           taubyMediumDeepTau2017v2VSe_ = analysisTree.tau_byMediumDeepTau2017v2VSe[indexTau] > 0.5;
           taubyMediumDeepTau2017v2VSjet_ = analysisTree.tau_byMediumDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyMediumDeepTau2017v2VSmu_ = analysisTree.tau_byMediumDeepTau2017v2VSmu[indexTau] > 0.5;
           taubyTightDeepTau2017v2VSe_ = analysisTree.tau_byTightDeepTau2017v2VSe[indexTau] > 0.5;
           taubyTightDeepTau2017v2VSjet_ = analysisTree.tau_byTightDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyTightDeepTau2017v2VSmu_ = analysisTree.tau_byTightDeepTau2017v2VSmu[indexTau] > 0.5;
           taubyVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVLooseDeepTau2017v2VSe[indexTau] > 0.5;
           taubyVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVLooseDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyVLooseDeepTau2017v2VSmu_ = analysisTree.tau_byVLooseDeepTau2017v2VSmu[indexTau] > 0.5;
           taubyVTightDeepTau2017v2VSe_ = analysisTree.tau_byVTightDeepTau2017v2VSe[indexTau] > 0.5;
           taubyVTightDeepTau2017v2VSjet_ = analysisTree.tau_byVTightDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyVVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVVLooseDeepTau2017v2VSe[indexTau] > 0.5;
           taubyVVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVVLooseDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyVVTightDeepTau2017v2VSe_ = analysisTree.tau_byVVTightDeepTau2017v2VSe[indexTau] > 0.5;
           taubyVVTightDeepTau2017v2VSjet_ = analysisTree.tau_byVVTightDeepTau2017v2VSjet[indexTau] > 0.5;
           taubyVVVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVVVLooseDeepTau2017v2VSe[indexTau] > 0.5;
           taubyVVVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVVVLooseDeepTau2017v2VSjet[indexTau] > 0.5;

           // bool isSingleTau = false;
           // bool isSingleTauOneProng = false;
           // for (unsigned int iT=0; iT<analysisTree.trigobject_count;++iT) {
           //   double dR = deltaR(lorentzVectorTau.Eta(),lorentzVectorTau.Phi(),
           // 		     analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
           //   if (dR>0.5) continue;
           //   if (isSinglePFTau180Trk50Filter) {
           //     if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50Filter]) isSingleTau = true;
           //   }
           //   if (isSinglePFTau180Trk50oneprongFilter) {
           //     if (analysisTree.trigobject_filters[iT][nSinglePFTau180Trk50oneprongFilter]) isSingleTauOneProng = true;
           //   }
           // }
           // tauSinglePFTau180Trk50_ = isSingleTau;
           // tauSinglePFTau180Trk50oneprong_ = isSingleTauOneProng;
           
           // finding matching jet
           bool jetFound = false;
           int indexMatchingJet = MatchingJet(analysisTree, lorentzVectorTau, jetFound, lorentzVectorTauJet);
           
           tauJetFlavor_ = 0;
           tauJetTightId_ = true;
           if (jetFound) {
              tauJetFlavor_ = analysisTree.pfjet_flavour[indexMatchingJet];
              if (era == "2017") tauJetTightId_ = tightJetiD_2017(analysisTree,indexMatchingJet);
              else if (era == "2018") tauJetTightId_ = tightJetiD_2018(analysisTree,indexMatchingJet);
	      else if (era == "2016") tauJetTightId_ = tightJetiD_2016(analysisTree,indexMatchingJet); // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
              else {
                 std::cout<<"Jet Id not set for jets faking taus in era "<<era<<std::endl;
                 exit(-1);
              }
           }
           else {
              lorentzVectorTauJet = lorentzVectorTau;
           }
           
           tauJetPt_  = lorentzVectorTauJet.Pt();
           tauJetEta_ = lorentzVectorTauJet.Eta();
           tauJetPhi_ = lorentzVectorTauJet.Phi();
           
        }
        
        if (debug) std::cout << "end of accessing taus" << std::endl;
        
        // ****************************
        // ****** trigger weight ******
        // ****************************
        float trigEffData = 1.0;
        float trigEffMC   = 1.0;
        
        trigWeight_ = 1;
        
        for (auto const& it : map_trigEffMC)
           {
              if( mhtNoMu_ < it.first)
                 { 
                    trigEffMC    =  it.second->Eval(metNoMu_);
                    trigEffData  =  map_trigEffData[it.first]->Eval(metNoMu_);
                    break;
                 }
           }
        
        if(trigEffMC !=0 ) trigWeight_ = trigEffData / trigEffMC;
        if(trigWeight_ < 0) trigWeight_=0;
        if(mhtNoMu_<mhtNoMu_Trig || metNoMu_< metNoMu_Trig) trigWeight_=0;
        
        if (debug) {
           cout << "MetNoMu  = " << metNoMu_
                << "  MhtNoMu  = " << mhtNoMu_
                << "  trigWeight = " << trigWeight_ << endl;
        }
        weight_ *= trigWeight_;
        
        // setting met filters
        metFilters_ = metFiltersPasses(analysisTree,metFlags,isData);
        
        // ********************************
        // **** filling trigger ntuple ****
        // ********************************
        if (ptTriggerMu>ptTrigMuCut && metNoMu_>metNoMu_Trig && mhtNoMu_>mhtNoMu_Trig) {
           trigNTuple_->Fill();
           TrigEvents++;
        }
        
        // ***************************
        // ******* WJet selection ****
        // ***************************
        if (lorentzVectorW.Pt()>1e-4) {
           recoilDPhi_  = dPhiFromLV(lorentzVectorW,lorentzVectorTau);
           recoilJetRatio_ = lorentzVectorTauJet.Pt()/lorentzVectorW.Pt();
           recoilJetDPhi_ = dPhiFromLV(lorentzVectorW,lorentzVectorTauJet);
           isWJet = ptTriggerMu>ptMuCut_WJet; 
           isWJet = isWJet && met_>metCut_WJet;
           isWJet = isWJet && tauPt_>ptTauCut_WJet;
           isWJet = isWJet && nMuon_ == 1;
           isWJet = isWJet && nElec_ == 0;
           isWJet = isWJet && nSelTaus_ == 1;
           isWJet = isWJet && nJetsCentral30_ == 1;
           isWJet = isWJet && abs(muonEta_)<etaMuCut_WJet;
           isWJet = isWJet && (tauDM_ || tauNewDM_);
           
           
           if (isWJet) {
              if (!isData) {
                 mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
                 mutrigweight = SF_muonTrig->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
              }
              HtNoRecoil_     = Ht_     - ptTriggerMu;
              SoftHtNoRecoil_ = SoftHt_ - ptTriggerMu;
              recoilM_   = lorentzVectorW.M();
              recoilPt_  = lorentzVectorW.Pt();
              recoilEta_ = lorentzVectorW.Eta();
              recoilPhi_ = lorentzVectorW.Phi();
              dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
              selection_ = 1;
              ntuple_->Fill();
              WJetEvents++;
           }
        }
        
        // ********************************
        // ******* W*->MuNu selection *****
        // ********************************
        if (lorentzVectorMet.Pt()>1e-4) {
           recoilDPhi_  = dPhiFromLV(lorentzVectorTriggerMu,lorentzVectorMet);
           recoilJetRatio_ = -1;
           recoilJetDPhi_  = 0;
           isWMuNu = ptTriggerMu>ptMuCut_WMuNu;
           isWMuNu = isWMuNu && met_>metCut_WMuNu;
           isWMuNu = isWMuNu && nMuon_ == 1;
           isWMuNu = isWMuNu && nElec_ == 0;
           isWMuNu = isWMuNu && nSelTaus_ == 0;
           isWMuNu = isWMuNu && nJetsCentral30_ == 0;
           isWMuNu = isWMuNu && abs(muonEta_)<etaMuCut_WMuNu;
           
           if (isWMuNu) {
              if (!isData) {
                 mueffweight  = SF_muonIdIso->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
                 mutrigweight = SF_muonTrig->get_ScaleFactor(ptTriggerMu, etaTriggerMu);
            }
              HtNoRecoil_     = Ht_;
              SoftHtNoRecoil_ = SoftHt_;
              recoilM_   = lorentzVectorMet.M();
              recoilPt_  = lorentzVectorMet.Pt();
              recoilEta_ = lorentzVectorMet.Eta();
              recoilPhi_ = lorentzVectorMet.Phi();
              selection_ = 2;
              ntuple_->Fill();
              WMuNuEvents++;
           }
        }
        
        // ********************************
        // ****** W*->TauNu selection *****
        // ******************************** 
        if (lorentzVectorMet.Pt()>1e-4) {
           recoilDPhi_  = dPhiFromLV(lorentzVectorTau,lorentzVectorMet);
           recoilJetRatio_ = lorentzVectorTauJet.Pt()/lorentzVectorMet.Pt();
           recoilJetDPhi_ = dPhiFromLV(lorentzVectorMet,lorentzVectorTauJet);
           isWTauNu = met_>metCut_WTauNu;
           isWTauNu = isWTauNu && tauPt_>ptTauCut_WTauNu;
           isWTauNu = isWTauNu && nSelTaus_ >= 1;
           isWTauNu = isWTauNu && nJetsCentral30_==1;
           isWTauNu = isWTauNu && nMuon_ == 0;
           isWTauNu = isWTauNu && nElec_ == 0;
           isWTauNu = isWTauNu && (tauDM_ || tauNewDM_);
           isWTauNu = isWTauNu && trigger_;
         
           if (isWTauNu) {
              HtNoRecoil_     = Ht_;
              SoftHtNoRecoil_ = SoftHt_;
              recoilM_   = lorentzVectorMet.M();
              recoilPt_  = lorentzVectorMet.Pt();
              recoilEta_ = lorentzVectorMet.Eta();
              recoilPhi_ = lorentzVectorMet.Phi();
              dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
              selection_ = 3;
              ntuple_->Fill();
              WTauNuEvents++;
           }
        }

        
        // *********************************
        // ****** Jet+Tau selection ********
        // *********************************
        isDiJet = nSelTaus_>0 && triggerJetsIndexes.size()>0;
        bool foundJetTauPair = false;
        if (isDiJet) {
           for (unsigned int iTau=0; iTau<tauIndexes.size(); ++iTau) { // loop over taus
              
              unsigned int indexTau = tauIndexes.at(iTau);
            
              if(nMuon_!=0)          continue;
              if(nElec_!=0)          continue;
              if(nSelTaus_!=1)       continue;
              if(nJetsCentral30_!=2) continue;
              if(analysisTree.tau_pt[indexTau]<ptTauCut_DiJet) continue;
              if(!pfJet60_ && !pfJet80_ && !pfJet140_ && !pfJet200_ && !pfJet260_ && !pfJet320_ && !pfJet400_ && !pfJet450_ ) continue;
              
              TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[indexTau],
                                                  analysisTree.tau_py[indexTau],
                                                  analysisTree.tau_pz[indexTau],
                                                  analysisTree.tau_mass[indexTau]);
              // finding recoiling jet 
              int acceptedJetIndex = -1;
              int acceptedJetDirIndex = -1;
              TLorentzVector recoilJetLV; recoilJetLV.SetXYZT(0,0,0,0);
              float dPhiTauJetMax = deltaPhiTauJetCut_DiJet;
              
              for (unsigned int iJet=0; iJet<triggerJetsIndexes.size(); ++iJet) { // loop over jets
                 unsigned int indexJet = triggerJetsIndexes.at(iJet);
                 TLorentzVector jetLV; jetLV.SetXYZT(analysisTree.pfjet_px[indexJet],
                                                     analysisTree.pfjet_py[indexJet],
                                                     analysisTree.pfjet_pz[indexJet],
                                                     analysisTree.pfjet_e[indexJet]);
                 if (jetLV.Pt()<ptJetCut_DiJet) continue;
                 if (fabs(jetLV.Eta())>etaJetCut_DiJet) continue;
                 float ptTauJetRatio = tauLV.Pt() / jetLV.Pt();
                 //if (ptTauJetRatio<ptTauJetRatioLowerCut_DiJet) continue;
                 //if (ptTauJetRatio>ptTauJetRatioUpperCut_DiJet) continue;
                 float dPhiTauJet = dPhiFromLV(tauLV,jetLV);
                 if (dPhiTauJet>dPhiTauJetMax) {
                    acceptedJetIndex = int(indexJet);
                    acceptedJetDirIndex = int(iJet);
                    dPhiTauJetMax = dPhiTauJet;
                    recoilJetLV = jetLV;
                 }
              } // end loop over jet
              if (acceptedJetIndex>=0) { // recoil jet found
                 
                 foundJetTauPair =  true;
                 
                 mttau_ = mT(tauLV,lorentzVectorMet);
                 mtgen_ = mT(wgentauLV,wnuLV);
                 tauPt_ = analysisTree.tau_pt[indexTau];
                 tauPz_ = analysisTree.tau_pz[indexTau];
                 tauEta_ = analysisTree.tau_eta[indexTau];
                 tauPhi_ = analysisTree.tau_phi[indexTau];
                 tauMass_ = analysisTree.tau_mass[indexTau];
                 tauQ_ = int(analysisTree.tau_charge[indexTau]);
                 tauNtrk1_ = analysisTree.tau_ntracks_pt1[indexTau];
                 tauNtrk05_ = analysisTree.tau_ntracks_pt05[indexTau];
                 tauNtrk08_ = analysisTree.tau_ntracks_pt08[indexTau];
                 
                 // finding matching jet
                 bool jetFound = false;
                 float dRmin = 0.4;
                 int indexMatchingJet = -1;
                 for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
                    TLorentzVector lorentzVectorJ; lorentzVectorJ.SetXYZT(analysisTree.pfjet_px[ijet],
                                                                          analysisTree.pfjet_py[ijet],
                                                                          analysisTree.pfjet_pz[ijet],
                                                                          analysisTree.pfjet_e[ijet]);
                    float drJetTau = deltaR(lorentzVectorJ.Eta(),lorentzVectorJ.Phi(),
                                            tauEta_,tauPhi_);
                    
                    if (drJetTau<dRmin) {
                       dRmin = drJetTau;
                       jetFound = true;
                       indexMatchingJet = ijet;
                       lorentzVectorTauJet = lorentzVectorJ;
                    }
                    
                 }
                 
                 
                 if (!jetFound) {
                    lorentzVectorTauJet = lorentzVectorTau;
                    continue;
                 }
                 
                 tauJetPt_  = lorentzVectorTauJet.Pt();
                 tauJetEta_ = lorentzVectorTauJet.Eta();
                 tauJetPhi_ = lorentzVectorTauJet.Phi();
                 if (era == "2017") tauJetTightId_ = tightJetiD_2017(analysisTree,indexMatchingJet);
                 else if (era == "2018")  tauJetTightId_ = tightJetiD_2018(analysisTree,indexMatchingJet);
		 else if (era == "2016") tauJetTightId_ = tightJetiD_2016(analysisTree,indexMatchingJet);
                 else {
                    std::cout<<"Jet Id for tau faking jets not set (in jet + tau selection) in era "<<era<<std::endl;
                    exit(-1);
                 }
                 
                 tauLeadingTrackPt_ = PtoPt(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                            analysisTree.tau_leadchargedhadrcand_py[indexTau]);
                 
                 tauLeadingTrackEta_ = PtoEta(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                              analysisTree.tau_leadchargedhadrcand_py[indexTau],
                                              analysisTree.tau_leadchargedhadrcand_pz[indexTau]);
                 
                 tauLeadingTrackPhi_ = PtoPhi(analysisTree.tau_leadchargedhadrcand_px[indexTau],
                                              analysisTree.tau_leadchargedhadrcand_py[indexTau]);
                 
                 tauLeadingTrackDz_  = analysisTree.tau_leadchargedhadrcand_dz[indexTau];
                 tauLeadingTrackDxy_ = analysisTree.tau_leadchargedhadrcand_dxy[indexTau];
                 
                 tauDecay_ = analysisTree.tau_decayMode[indexTau];
                 tauGenDecay_ = analysisTree.tau_genDecayMode[indexTau];
                 tauGenMatchDecay_ = tauGenMatchDecay.at(iTau);
                 
                 if (tauDecay_<0) tauDecay_ = -1;
                 if (tauGenDecay_<0) tauGenDecay_ = -1;
                 if (tauGenMatchDecay_<0) tauGenMatchDecay_ = -1;
                 
                 pfJet60_ = jets60trigger.at(acceptedJetDirIndex);
                 pfJet80_ = jets80trigger.at(acceptedJetDirIndex);
                 pfJet140_ = jets140trigger.at(acceptedJetDirIndex);
                 pfJet200_ = jets200trigger.at(acceptedJetDirIndex);
                 pfJet260_ = jets260trigger.at(acceptedJetDirIndex);
                 pfJet320_ = jets320trigger.at(acceptedJetDirIndex);
                 pfJet400_ = jets400trigger.at(acceptedJetDirIndex);
                 pfJet450_ = jets450trigger.at(acceptedJetDirIndex);
                 
                 tauDM_ = analysisTree.tau_decayModeFinding[indexTau] > 0.5;
                 tauNewDM_ = analysisTree.tau_decayModeFindingNewDMs[indexTau] > 0.5;
                 
                 tauIso_ =  analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[indexTau];
                 tauLooseIso_ = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
                 tauMediumIso_ = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
                 tauTightIso_ = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[indexTau] > 0.5;
                 
                 tauVLooseMvaIso_ = analysisTree.tau_byVLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 tauLooseMvaIso_ = analysisTree.tau_byLooseIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 tauMediumMvaIso_ = analysisTree.tau_byMediumIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 tauTightMvaIso_ = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 tauVTightMvaIso_ = analysisTree.tau_byVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 tauVVTightMvaIso_ = analysisTree.tau_byVVTightIsolationMVArun2v1DBoldDMwLT[indexTau] > 0.5;
                 
                 tauVVLooseMva2017v2Iso_ = analysisTree.tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauVLooseMva2017v2Iso_ = analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauLooseMva2017v2Iso_ = analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauMediumMva2017v2Iso_ = analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauTightMva2017v2Iso_ = analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauVTightMva2017v2Iso_ = analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 tauVVTightMva2017v2Iso_ = analysisTree.tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[indexTau] > 0.5;
                 
                 tauAntiMuonLoose3_ = analysisTree.tau_againstMuonLoose3[indexTau] > 0.5;
                 tauAntiMuonTight3_ = analysisTree.tau_againstMuonTight3[indexTau] > 0.5;

                 taubyDeepTau2017v2VSeraw_ = analysisTree.tau_byDeepTau2017v2VSeraw[indexTau] > 0.5;
                 taubyDeepTau2017v2VSjetraw_ = analysisTree.tau_byDeepTau2017v2VSjetraw[indexTau] > 0.5;
                 taubyDeepTau2017v2VSmuraw_ = analysisTree.tau_byDeepTau2017v2VSmuraw[indexTau] > 0.5;
                 taubyLooseDeepTau2017v2VSe_ = analysisTree.tau_byLooseDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyLooseDeepTau2017v2VSjet_ = analysisTree.tau_byLooseDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyLooseDeepTau2017v2VSmu_ = analysisTree.tau_byLooseDeepTau2017v2VSmu[indexTau] > 0.5;
                 taubyMediumDeepTau2017v2VSe_ = analysisTree.tau_byMediumDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyMediumDeepTau2017v2VSjet_ = analysisTree.tau_byMediumDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyMediumDeepTau2017v2VSmu_ = analysisTree.tau_byMediumDeepTau2017v2VSmu[indexTau] > 0.5;
                 taubyTightDeepTau2017v2VSe_ = analysisTree.tau_byTightDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyTightDeepTau2017v2VSjet_ = analysisTree.tau_byTightDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyTightDeepTau2017v2VSmu_ = analysisTree.tau_byTightDeepTau2017v2VSmu[indexTau] > 0.5;
                 taubyVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVLooseDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVLooseDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyVLooseDeepTau2017v2VSmu_ = analysisTree.tau_byVLooseDeepTau2017v2VSmu[indexTau] > 0.5;
                 taubyVTightDeepTau2017v2VSe_ = analysisTree.tau_byVTightDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyVTightDeepTau2017v2VSjet_ = analysisTree.tau_byVTightDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyVVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVVLooseDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyVVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVVLooseDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyVVTightDeepTau2017v2VSe_ = analysisTree.tau_byVVTightDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyVVTightDeepTau2017v2VSjet_ = analysisTree.tau_byVVTightDeepTau2017v2VSjet[indexTau] > 0.5;
                 taubyVVVLooseDeepTau2017v2VSe_ = analysisTree.tau_byVVVLooseDeepTau2017v2VSe[indexTau] > 0.5;
                 taubyVVVLooseDeepTau2017v2VSjet_ = analysisTree.tau_byVVVLooseDeepTau2017v2VSjet[indexTau] > 0.5;

                 recoilDPhi_ = dPhiFromLV(tauLV,recoilJetLV);
                 recoilJetRatio_ = lorentzVectorTauJet.Pt()/recoilJetLV.Pt();
                 recoilJetDPhi_ = dPhiFromLV(lorentzVectorTauJet,recoilJetLV);
                 recoilM_ = recoilJetLV.M();
                 recoilPt_ = recoilJetLV.Pt();
                 recoilEta_ = recoilJetLV.Eta();
                 recoilPhi_ = recoilJetLV.Phi();
                 HtNoRecoil_     = Ht_     - recoilJetLV.Pt();
                 SoftHtNoRecoil_ = SoftHt_ - recoilJetLV.Pt();
                 dPhiMetTau_= dPhiFromLV(lorentzVectorMet,lorentzVectorTau);
                 selection_ = 4;
                 if(!tauDM_ && !tauNewDM_) continue;
                 ntuple_->Fill();
              }
           } // END: loop over tau
        }
        if (foundJetTauPair) JetTauEvents++;
        
     } // end of file processing (loop over events in one file)
     nFiles++;
     delete tree_;
     file_->Close();
     delete file_;
  }
  
  int allEvents   = int(inputEventsH->GetSumOfWeights());
  double sumOfWeights = histWeightsH->GetSumOfWeights();
  std::cout << "Total number of input events      = " << allEvents << std::endl;
  std::cout << "Total weight sum                  = " << sumOfWeights << std::endl;
  std::cout << "Total number of events in Tree    = " << nEvents << std::endl;
  std::cout << "Total number of trigger events    = " << TrigEvents << std::endl;
  std::cout << "Total number of W+Jet events      = " << WJetEvents << std::endl;
  std::cout << "Total number of W->muv events     = " << WMuNuEvents << std::endl;
  std::cout << "Total number of W->tauv events    = " << WTauNuEvents << std::endl;
  std::cout << "Total number of single jet events = " << SingleJetEvents << std::endl;
  std::cout << "Total number of jet+tau events    = " << JetTauEvents << std::endl;
  std::cout << "Total number of dijet events      = " << DiJetEvents << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}



