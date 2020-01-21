#include "SynchNTupleProducer_em_allEras_definitions.h"
#include "NtupleMaker_TauID_all_eras_functions.h"

int main(int argc, char * argv[]) {

    // first argument - config file
    // second argument - filelist
   bool sync = false;  
   
   // read config file =======================================================================================================================================================
   Config cfg(argv[1]);

   const string era = cfg.get<string>("Era");
   
   const bool computeSVFitMass = cfg.get<bool>("ComputeSVFitMass");
   const bool computeFastMTTMass = cfg.get<bool>("ComputeFastMTTMass");
   const bool removeGammaStar = cfg.get<bool>("RemoveGammaStar");
   const bool usePuppiMet = cfg.get<bool>("UsePuppiMET");

   const bool isData = cfg.get<bool>("IsData");
   const bool isDY   = cfg.get<bool>("IsDY");
   const bool isW    = cfg.get<bool>("IsW");
   const bool isEmbedded = cfg.get<bool>("IsEmbedded");
   const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
   const string jsonFile = cfg.get<string>("jsonFile");
   const bool applySimpleRecoilCorrections = cfg.get<bool>("ApplySimpleRecoilCorrections");
   const bool isSignal = cfg.get<bool>("IsSignal");
   
   // kinematic cuts on electrons
   const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
   const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
   const float etaElectronCut     = cfg.get<float>("etaElectronCut");
   const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
   const float dzElectronCut      = cfg.get<float>("dzElectronCut");
   const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
   const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
   const string lowPtLegElectron  = cfg.get<string>("LowPtLegElectron");
   const string highPtLegElectron = cfg.get<string>("HighPtLegElectron");
   
   // veto electrons
   const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
   const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
   const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
   const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
   const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
   const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");
   
   // kinematic cuts on muons
   const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
   const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
   const float etaMuonCut     = cfg.get<float>("etaMuonCut");
   const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
   const float dzMuonCut      = cfg.get<float>("dzMuonCut");
   const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
   const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
   const string lowPtLegMuon  = cfg.get<string>("LowPtLegMuon");
   const string highPtLegMuon = cfg.get<string>("HighPtLegMuon");
   
   // veto muons
   const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
   const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
   const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
   const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
   const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
   const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");
    
   // vertex cuts
   const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
   const float zVertexCut     = cfg.get<float>("ZVertexCut");
   const float dVertexCut     = cfg.get<float>("DVertexCut");
   
   // topological cuts
   const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
   const bool isMuonIsoR03 = cfg.get<bool>("IsMuonIsoR03");
   const bool isElectronIsoR03 = cfg.get<bool>("IsElectronIsoR03");
   const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
   const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
   const bool applyDzFilterMatch = cfg.get<bool>("ApplyDzFilterMatch");
   const string mu23ele12DzFilter = cfg.get<string>("Mu23Ele12DzFilter");
   const string mu8ele23DzFilter = cfg.get<string>("Mu8Ele23DzFilter");
    
    // jets
   const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
   const string bTagFile = cfg.get<string>("BTagFile");
   const string bTagEffFile = cfg.get<string>("BTagEffFile");
   const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
   const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
   const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");
   const float jetEtaCut = cfg.get<float>("JetEtaCut");
   const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
   const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
   const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
   const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
   const float btagCut = cfg.get<float>("btagCut");
   const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
   const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");
   
   const string jec_UncertaintySources = cfg.get<string>("JEC_UncertaintySources");
   
   TString LowPtLegElectron(lowPtLegElectron);
   TString HighPtLegElectron(highPtLegElectron);
   
   TString LowPtLegMuon(lowPtLegMuon);
   TString HighPtLegMuon(highPtLegMuon);
   
   TString Mu23Ele12DzFilter(mu23ele12DzFilter);
   TString Mu8Ele23DzFilter(mu8ele23DzFilter);
   
   TString BTagDiscriminator1(bTagDiscriminator1);
   TString BTagDiscriminator2(bTagDiscriminator2);
   TString BTagDiscriminator3(bTagDiscriminator3);
   
   const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
   const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEff");
   
   const string Muon23TriggerFile = cfg.get<string>("Muon23TriggerEff");
   const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");
   
   const string Electron23TriggerFile = cfg.get<string>("Electron23TriggerEff");
   const string Electron12TriggerFile = cfg.get<string>("Electron12TriggerEff");
   
   const float muonScale = cfg.get<float>("MuonScale");
   const float eleScaleBarrel = cfg.get<float>("EleScaleBarrel");
   const float eleScaleEndcap = cfg.get<float>("EleScaleEndcap");
   
   const string recoilFileName   = cfg.get<string>("RecoilFileName");
   TString RecoilFileName(recoilFileName);
   const string recoilFileNamePuppi   = cfg.get<string>("RecoilFileNamePuppi");
   TString RecoilFileNamePuppi(recoilFileNamePuppi);
   
   const string metSysFileName   = cfg.get<string>("MetSysFileName");
   TString MetSysFileName(metSysFileName);
   const string metSysFileNamePuppi   = cfg.get<string>("MetSysFileNamePuppi");
   TString MetSysFileNamePuppi(metSysFileNamePuppi);
   
   const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
   TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
   
   const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
   TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
   
   const string pileUpDataFile = cfg.get<string>("PileUpDataFile");
   const string pileUpMCFile = cfg.get<string>("PileUpMCFile");
   const string samplenameForPUHist = cfg.get<string>("SampleNameForPUHist");
   TString PileUpDataFile(pileUpDataFile);
   TString PileUpMCFile(pileUpMCFile);
   
   const string correctionWSFile = cfg.get<string>("CorrectionWSFile");

   const bool apply_ggh_reweighting = cfg.get<bool>("ApplygghReweighting");
   const bool apply_ggh_uncertainties = cfg.get<bool>("ApplygghUncertainties");
   const bool apply_vbf_uncertainties = cfg.get<bool>("ApplyVBFUncertainties");
 
   const bool  isSampleForRecoilCorrection = (isW||isDY||isSignal)&&!isData;
   
   
   // Open root files and setup trees  =================================================================================================================================
   std::string rootFileName(argv[2]);
   std::ifstream fileList(argv[2]);
   std::ifstream fileList0(argv[2]);
   std::string ntupleName("makeroottree/AC1B");
   std::string initNtupleName("initroottree/AC1B");
   
   TString TStrName(rootFileName);
   std::cout <<TStrName <<std::endl;
   
   // output fileName with histograms ==================================================================================================================================
   TFile * file = new TFile(TStrName+TString(".root"),"recreate");
   file->cd("");

   TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
   TH1D * ZMM = new TH1D("ZMM","",1,-0.5,0.5);
   TH1D * ZEE = new TH1D("ZEE","",1,-0.5,0.5);
   TH1D * ZTT = new TH1D("ZTT","",1,-0.5,0.5);
   TH1D * ALL = new TH1D("ALL","",1,-0.5,0.5);
   TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
       
   tree = new TTree("TauCheck","TauCheck");
   SetupTree();
    
   // set project directory ============================================================================================================================================
   string cmsswBase = (getenv ("CMSSW_BASE"));

   // prepare uncertainties (NEW) ======================================================================================================================================
   for(auto &uncert : uncertainty_map){
      int count_unc = 0;
      for(auto &unc_var : unc_vars){
         tree->Branch(unc_var+"_"+uncert.first,
                      &uncert.second.container[count_unc],
                      unc_var+"_"+uncert.first+"/F");
         count_unc += 1;
      }
   }
   
   // JEC uncertainties (NEW) ==========================================================================================================================================
   for (int isrc = 0; isrc < nsrc_Eta0To5; isrc++) {
      const char *name = srcnames_Eta0To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To5[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_Eta0To3; isrc++) {
      const char *name = srcnames_Eta0To3[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To3[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_Eta3To5; isrc++) {
      const char *name = srcnames_Eta3To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta3To5[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_RelativeBal; isrc++) {
      const char *name = srcnames_RelativeBal[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_RelativeBal[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_RelativeSampleYear; isrc++) {
      const char *name = srcnames_RelativeSampleYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_RelativeSampleYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_EC2; isrc++) {
      const char *name = srcnames_EC2[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_EC2[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_FlavorQCD; isrc++) {
      const char *name = srcnames_FlavorQCD[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_FlavorQCD[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_AbsoluteYear; isrc++) {
      const char *name = srcnames_AbsoluteYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_AbsoluteYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_HFYear; isrc++) {
      const char *name = srcnames_HFYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_HFYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_EC2Year; isrc++) {
      const char *name = srcnames_EC2Year[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_EC2Year[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_BBEC1Year; isrc++) {
      const char *name = srcnames_BBEC1Year[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_BBEC1Year[isrc] = unc;
   }

   map< TString , vector<JetCorrectionUncertainty*> > jec_unc_map = {
      { "jecUncEta0To5"     , vsrc_Eta0To5 },
      { "jecUncEta0To3"     , vsrc_Eta0To3 },
      { "jecUncEta3To5"     , vsrc_Eta3To5 },
      { "jecUncRelativeBal" , vsrc_RelativeBal }, 
      { "jecUncRelativeSampleYear", vsrc_RelativeSampleYear },
      { "jecUncEC2"         , vsrc_EC2 }, 
      { "jecUncFlavorQCD"   , vsrc_FlavorQCD }, 
      { "jecUncAbsoluteYear"   , vsrc_AbsoluteYear }, 
      { "jecUncHFYear"   , vsrc_HFYear }, 
      { "jecUncEC2Year"   , vsrc_EC2Year }, 
      { "jecUncBBEC1Year"   , vsrc_BBEC1Year }, 
   };
      
   // initialize good run selection ====================================================================================================================================
   std::vector<Period> periods;
   string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
   if (isData || isEmbedded ) ReadJson(periods,fullPathToJsonFile);
   
   // Load CrystalBallEfficiency class =================================================================================================================================
   TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
   LoadCrystalBallEfficiencyClass(pathToCrystalLib);

   // PU reweighting  ==================================================================================================================================================
   PileUp * PUofficial = new PileUp();
   TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read");
   TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read");
   if (!isData && !isEmbedded) iniPU(PUofficial,filePUdistribution_data,filePUdistribution_MC, samplenameForPUHist);
       
   // load Lepton Scale Factors ========================================================================================================================================
   ScaleFactor * SF_muonIdIso = new ScaleFactor();
   SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
   ScaleFactor * SF_electronIdIso = new ScaleFactor();
   SF_electronIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
   ScaleFactor * SF_muon23 = new ScaleFactor();
   SF_muon23->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon23TriggerFile));
   ScaleFactor * SF_muon8 = new ScaleFactor();
   SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Muon8TriggerFile));
   ScaleFactor * SF_electron23 = new ScaleFactor();
   SF_electron23->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron23TriggerFile));
   ScaleFactor * SF_electron12 = new ScaleFactor();
   SF_electron12->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(Electron12TriggerFile));
    
   // set MET filters ==================================================================================================================================================
   // for recommendations see : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
   std::vector<TString> metFlags; metFlags.clear();
   if (era=="2016"){
      metFlags.push_back("Flag_HBHENoiseFilter");
      metFlags.push_back("Flag_HBHENoiseIsoFilter");
      metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
      metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
      metFlags.push_back("Flag_goodVertices");
      if (isData || isEmbedded)
         metFlags.push_back("Flag_eeBadScFilter");
      metFlags.push_back("Flag_BadPFMuonFilter");
      // metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
   }
   else if (era=="2017"){
      metFlags.push_back("Flag_HBHENoiseFilter");
      metFlags.push_back("Flag_HBHENoiseIsoFilter");
      metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
      metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
      metFlags.push_back("Flag_goodVertices");
      if (isData || isEmbedded)
         metFlags.push_back("Flag_eeBadScFilter");
      metFlags.push_back("Flag_BadPFMuonFilter");
      //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
      metFlags.push_back("ecalBadCalibReducedMINIAODFilter");//version we currently use is outdated, will be updated with next ntuple production campaign
   }
   else if (era == "2018"){
      metFlags.push_back("Flag_goodVertices");
      metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
      metFlags.push_back("Flag_HBHENoiseFilter");
      metFlags.push_back("Flag_HBHENoiseIsoFilter");
      metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
      metFlags.push_back("Flag_BadPFMuonFilter");
      //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
      if (isData || isEmbedded)
         metFlags.push_back("Flag_eeBadScFilter");
      metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
   }
   else {
      std::cout << "MET filters not defined for era "<<era<<std::endl;
      exit(-1);
   }
   std::vector<TString> badChargedCandidateFlag; badChargedCandidateFlag.clear();
   badChargedCandidateFlag.push_back("Flag_BadChargedCandidateFilter");//do not apply? not applied for now, check!
   std::vector<TString> badPFMuonFlag; badPFMuonFlag.clear();
   badPFMuonFlag.push_back("Flag_BadPFMuonFilter");//do not apply? not applied for now, check!
   std::vector<TString> badMuonFlag; badMuonFlag.clear();
   badMuonFlag.push_back("Flag_badMuons");//do not apply? not applied for now, check! 

   // initialize recoil corrections ====================================================================================================================================
   kit::MEtSys metSys(MetSysFileName);
   kit::RecoilCorrector recoilMetCorrector(RecoilFileName);
   kit::MEtSys metSysPuppi(MetSysFileNamePuppi);
   kit::RecoilCorrector recoilMetCorrectorPuppi(RecoilFileNamePuppi);

   // initialize calibration of impact parameters ======================================================================================================================
   //CalibrationOfImpactParameters calibrateIP;
   
   // SV fit mass   // FIXME: not needed for ClassicSV FIT, can be removed =============================================================================================
   edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
   TH1::AddDirectory(false);
   TFile * inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
       
   // initialize BTag scale factors ====================================================================================================================================
   BTagCalibration calib(bTagAlgorithm, cmsswBase+"/src/"+bTagFile);
   BTagCalibrationReader reader_BTAG(BTagEntry::OP_MEDIUM,"central",{"up","down"});
   TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/"+bTagEffFile));
      
   reader_BTAG.load(calib,BTagEntry::FLAV_B,"comb");
   reader_BTAG.load(calib,BTagEntry::FLAV_C,"comb");
   reader_BTAG.load(calib,BTagEntry::FLAV_UDSG,"incl");
   
   float etaBTAG[2] = {0.5,2.1};
   float ptBTAG[5] = {25.,35.,50.,100.,200.};
   
   std::cout << std::endl;
   for (int iEta=0; iEta<2; ++iEta) {
      for (int iPt=0; iPt<5; ++iPt) {
         float sfB = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
         float sfC = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
         float sfLight = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
         printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
      }
   }
   std::cout << std::endl;
   
   tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
   tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
   tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
       
   // initialize ggh reweighting =======================================================================================================================================
   TFile *file_ggh_reweighting = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/NNLOPS_reweight.root"));
   TGraph * gr_NNLOPSratio_pt_powheg_0jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_powheg_0jet");
   TGraph * gr_NNLOPSratio_pt_powheg_1jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_powheg_1jet");
   TGraph * gr_NNLOPSratio_pt_powheg_2jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_powheg_2jet");
   TGraph * gr_NNLOPSratio_pt_powheg_3jet  = (TGraph*) file_ggh_reweighting->Get("gr_NNLOPSratio_pt_powheg_3jet");
   
   // load Z pt mass weights ===========================================================================================================================================
   TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName);
   TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName);
   LoadZpTMassWeights(fileZMassPtWeights, ZMassPtWeightsFileName, histZMassPtWeights, ZMassPtWeightsHistName);
      
   // load correction workspaces =======================================================================================================================================
   TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/"+correctionWSFile;
   TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
   RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

   // store whether FastMTT or SVFit is used ===========================================================================================================================
   isSVFitUsed = computeSVFitMass;
   isFastMTTUsed = computeFastMTTMass;

   // store whether PuppiMET is used ===========================================================================================================================
   isPuppiMETUsed = usePuppiMet;
   
   // setup MELA =======================================================================================================================================================
   const int erg_tev = 13;
   const float mPOLE = 125.6;
   TVar::VerbosityLevel verbosity = TVar::SILENT;
   Mela mela(erg_tev, mPOLE, verbosity);

   // ==================================================================================================================================================================
   // open files =======================================================================================================================================================
   // ==================================================================================================================================================================
   int nTotalFiles = 0;
   std::string dummy;
   // count number of files --->
   while (fileList0 >> dummy) nTotalFiles++;
   int nEvents = 0;
   int selEvents = 0;
   int nFiles = 0;
   
   for (int iF=0; iF<nTotalFiles; ++iF) {
      
      std::string filen;
      fileList >> filen;
      
      std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
      TFile * file_ = TFile::Open(TString(filen));
        
      TTree * _inittree = NULL;
      _inittree = (TTree*)file_->Get(TString(initNtupleName));
      FillHistGenWeights( _inittree, histWeightsH, isData);
              
      TTree * _tree = NULL;
      _tree = (TTree*)file_->Get(TString(ntupleName));
      if (_tree==NULL) continue;
            
      if (!FillHistInputEvents(file_, inputEventsH)) continue;
             
      AC1B analysisTree(_tree);
      Long64_t numberOfEntries = analysisTree.GetEntries();
      std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
      
      // event loop ====================================================================================================================================================
      for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {

         analysisTree.GetEntry(iEntry);
         nEvents++;
         
         SetDefaultValues();
         
         if (analysisTree.event_run>maxRun)
            maxRun = analysisTree.event_run;
            
         if (analysisTree.event_run<minRun)
            minRun = analysisTree.event_run;
            
         if (nEvents%10000==0)
            cout << "      processed " << nEvents << " events" << endl;

         // set weights  ===============================================================================================================================================
         if (isEmbedded) mcweight = analysisTree.genweight;
         else {
            if (analysisTree.genweight < 0){
               mcweight = -1.;
            }
            else mcweight = 1.;
         }

         if (sync && mcweight<0) continue; // FIXME: needed for sync?

         weightScale1 = analysisTree.weightScale1;
         weightScale2 = analysisTree.weightScale2;
         weightScale3 = analysisTree.weightScale3;
         weightScale4 = analysisTree.weightScale4;
         weightScale5 = analysisTree.weightScale5;
         weightScale6 = analysisTree.weightScale6;
         weightScale7 = analysisTree.weightScale7;
         weightScale8 = analysisTree.weightScale8;
         
         weightPDFup   = analysisTree.weightPDFup;
         weightPDFdown = analysisTree.weightPDFdown;

	 if(era=="2016" || era=="2017"){
	   prefiringweight     = analysisTree.prefiringweight;
	   prefiringweightup   = analysisTree.prefiringweightup;
	   prefiringweightdown = analysisTree.prefiringweightdown;
	 }

         // store gen-info  ============================================================================================================================================
         if (!isData) {
                
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
               
            }  // end of loop over gen-taus

            for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
               
               TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                   analysisTree.genparticles_py[igen],
                                                   analysisTree.genparticles_pz[igen],
                                                   analysisTree.genparticles_e[igen]);
                    
               if (analysisTree.genparticles_pdgid[igen]==6)
                  topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
               
               if (analysisTree.genparticles_pdgid[igen]==-6)
                  antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                          analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
               
               if (analysisTree.genparticles_pdgid[igen]==22 && analysisTree.genparticles_status[igen]==44)
                  isGSfound = true;
               
               if (analysisTree.genparticles_pdgid[igen]==23) zBosonLV = genLV;
                              
               if (analysisTree.genparticles_pdgid[igen]==25||
                   analysisTree.genparticles_pdgid[igen]==35||
                   analysisTree.genparticles_pdgid[igen]==36) hBosonLV = genLV;
               
               if (abs(analysisTree.genparticles_pdgid[igen])==24) wBosonLV = genLV;
                           
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
               
            } // end of loop over gen-particles
            
            if (isGSfound) {
               //	  std::cout << "gamma* found : " << std::endl;
               if (removeGammaStar) continue;
            }

            if (isDY||isEmbedded) {
 
               if (promptTausFirstCopy.size()==2) {
                  isZTT = true; isZMM = false; isZEE = false;
                  bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
                  bosonMass = promptTausLV.M();
                  lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
                  mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);

                  // set embeddedWeight  ===============================================================================================================================
                  embeddedWeight = GetEmbeddedWeight( promptTausFirstCopy, correctionWS, era);

               }
               else if (promptMuons.size()==2) {
                  isZTT = false; isZMM = true; isZEE = false;
                  bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
                  bosonMass = promptMuonsLV.M();
                  lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
                  mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
               }
               else {
                  isZTT = false; isZMM = false; isZEE = true;
                  bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
                  bosonMass = promptElectronsLV.M();
                  lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
                  if (promptElectrons.size()==2)
                     mtBoson_gen = mT(promptElectrons[0],promptElectrons[1]);
               }
               nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
            } // if DY or embedded
            else if (isW) {
               bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
               bosonMass = wDecayProductsLV.M();
               if (promptTausLastCopy.size()==1) {
                  lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
               }
               else if (promptMuons.size()==1) {
                  lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
               }
               else {
                  lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
               }
               nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
            }
            else {
               TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
               bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
               TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
               lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
               nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
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

            // set zptmassweight  ======================================================================================================================================             
            if (isDY) zptmassweight = GetZptMassWeight(bosonMass, bosonPt, histZMassPtWeights);
 
            ALL->Fill(0.0);
            if (isZMM) ZMM->Fill(0.);
            if (isZEE) ZEE->Fill(0.);
            if (isZTT) ZTT->Fill(0.);
            isZLL = isZMM || isZEE;
         } // end of if !isData

         // ggh reweighting   && ggH and VBF uncertainties================================================================================================================
         if (isSignal){
            njets_HTXS = analysisTree.htxs_njets30;
            higgspt_HTXS = analysisTree.htxs_higgsPt;
            htxs_stage0cat = analysisTree.htxs_stage0cat;
            htxs_stage1cat = analysisTree.htxs_stage1p1cat;
            htxs_stage1p1cat = analysisTree.htxs_stage1p1cat;

            if (apply_ggh_reweighting)
               {
                  if      (njets_HTXS==0) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_powheg_0jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)125.0));
                  else if (njets_HTXS==1) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_powheg_1jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)625.0));
                  else if (njets_HTXS==2) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_powheg_2jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)800.0));
                  else if (njets_HTXS>=3) weight_ggh_NNLOPS = gr_NNLOPSratio_pt_powheg_3jet->Eval(TMath::Min(higgspt_HTXS,(Float_t)925.0));
                  else weight_ggh_NNLOPS = 1.0;
               }
            if (apply_ggh_uncertainties){
               std::vector<double> ggF_unc = qcd_ggF_uncertSF_2017(njets_HTXS, higgspt_HTXS, htxs_stage1cat, 1.0);
               THU_ggH_Mu = ggF_unc[0];
               THU_ggH_Res = ggF_unc[1];
               THU_ggH_Mig01 = ggF_unc[2];
               THU_ggH_Mig12 = ggF_unc[3];
               THU_ggH_VBF2j = ggF_unc[4];
               THU_ggH_VBF3j = ggF_unc[5];
               THU_ggH_PT60 = ggF_unc[6];
               THU_ggH_PT120 = ggF_unc[7];
               THU_ggH_qmtop = ggF_unc[8];
            }
            
            if (apply_vbf_uncertainties){
               THU_qqH_TOT = vbf_uncert_stage_1_1(0, htxs_stage1p1cat, 1.0);
               THU_qqH_PTH200 = vbf_uncert_stage_1_1(1, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj60 = vbf_uncert_stage_1_1(2, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj120 = vbf_uncert_stage_1_1(3, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj350 = vbf_uncert_stage_1_1(4, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj700 = vbf_uncert_stage_1_1(5, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj1000 =vbf_uncert_stage_1_1(6, htxs_stage1p1cat, 1.0);
               THU_qqH_Mjj1500 = vbf_uncert_stage_1_1(7, htxs_stage1p1cat, 1.0);
               THU_qqH_25 = vbf_uncert_stage_1_1(8, htxs_stage1p1cat, 1.0);
               THU_qqH_JET01 = vbf_uncert_stage_1_1(9, htxs_stage1p1cat, 1.0);
            }
                  
         }
         // store event variables  =====================================================================================================================================
         run = int(analysisTree.event_run);
         lumi = int(analysisTree.event_luminosityblock);
         evt = analysisTree.event_nr;
         npv = analysisTree.primvertex_count;
         npu = analysisTree.numtruepileupinteractions;
         rho = analysisTree.rho;
         npartons = analysisTree.genparticles_noutgoing;
         
         // apply good run selection  ==================================================================================================================================

         if ((isData || isEmbedded) && applyGoodRunSelection){

            int n=analysisTree.event_run;
            int lum = analysisTree.event_luminosityblock;      
            if (!GoodRunSelection(n,lum,periods)) continue;

         }

         // pileup and top pt re-weighting weight ======================================================================================================================
         if (!isData && !isEmbedded) {

            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
            if (topPt>0&&antitopPt>0) {
               topptweight = topPtWeight(topPt,antitopPt,true);
               topptweightRun2 = topPtWeight(topPt,antitopPt,false);
            }
         }

         // accessing trigger info =====================================================================================================================================
         isLowPtLegElectron = AccessTriggerInfo(analysisTree,LowPtLegElectron,nLowPtLegElectron);
         isHighPtLegElectron = AccessTriggerInfo(analysisTree,HighPtLegElectron,nHighPtLegElectron);
         isLowPtLegMuon = AccessTriggerInfo(analysisTree,LowPtLegMuon,nLowPtLegMuon);
         isHighPtLegMuon = AccessTriggerInfo(analysisTree,HighPtLegMuon,nHighPtLegMuon);
         if (applyDzFilterMatch){
         isMu23Ele12DzFilter = AccessTriggerInfo(analysisTree,Mu23Ele12DzFilter,nMu23Ele12DzFilter);
         isMu8Ele23DzFilter = AccessTriggerInfo(analysisTree,Mu8Ele23DzFilter,nMu8Ele23DzFilter);
         }
         if (!isLowPtLegElectron || !isHighPtLegElectron || !isLowPtLegMuon || !isHighPtLegMuon ) {
            std::cout << "PFJet HLT filter not found" << std::endl;
            exit(-1);
         }
         if (applyDzFilterMatch && (!isMu23Ele12DzFilter || !isMu8Ele23DzFilter) ) {
            std::cout << "PFJet HLT filter not found" << std::endl;
            exit(-1);
         }
         // searching for btagging discriminant ========================================================================================================================
         unsigned int nBTagDiscriminant1 = 0;
         unsigned int nBTagDiscriminant2 = 0;
         unsigned int nBTagDiscriminant3 = 0;
         SearchForBtagDiscriminant(analysisTree, BTagDiscriminator1, BTagDiscriminator2, BTagDiscriminator3, nBTagDiscriminant1, nBTagDiscriminant2, nBTagDiscriminant3, era);
         
         // electron selection =========================================================================================================================================

         vector<int> electrons; electrons.clear();
         for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
            if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
            if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
            if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
            if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
            bool electronMvaId = true;
            electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie]>0.5;
            if (!electronMvaId) continue;
            if (!analysisTree.electron_pass_conversion[ie]) continue;
            if (analysisTree.electron_nmissinginnerhits[ie]>1) continue;
            electrons.push_back(ie);
         }

         // muon selection =============================================================================================================================================

         vector<int> muons; muons.clear();
         for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
            if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;          
            if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
            if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
            if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
            if (analysisTree.muon_isBad[im]) badMuonFilter_ = false;                     
            if (analysisTree.muon_isDuplicate[im]) duplicateMuonFilter_ = false;                   
            bool muonId = analysisTree.muon_isMedium[im]; 
            if (!muonId) continue;
            muons.push_back(im);
         }
         // select only events with at least one electron and one muon =================================================================================================
         if (electrons.size()==0) continue; 
         if (muons.size()==0) continue;          


         // selecting muon and electron pair (OS or SS) ================================================================================================================

         int electronIndex = -1;
         int muonIndex = -1;
         float isoMuMin = 1e+10;
         float isoEleMin = 1e+10;
         SelectMuonElePair(analysisTree, muons, electrons, isMuonIsoR03, isElectronIsoR03, dRleptonsCut, ptMuonHighCut, ptElectronHighCut, electronIndex, muonIndex, isoMuMin, isoEleMin, era);
         if (electronIndex<0) continue;
         if (muonIndex<0) continue;

         // MET Filters ================================================================================================================================================             
         metFilters_ = metFiltersPasses(analysisTree,metFlags);
         badChargedCandidateFilter_ = metFiltersPasses(analysisTree,badChargedCandidateFlag);
         badPFMuonFilter_ = metFiltersPasses(analysisTree,badPFMuonFlag);
         badMuonFilter_ = metFiltersPasses(analysisTree,badMuonFlag);
         if (era=="2016") metFilters_ = metFilters_;// && badMuonFilter_ && duplicateMuonFilter_ &&  badChargedCandidateFilter_ && badPFMuonFilter_; not applied for now, check!              
         else metFilters_ = metFilters_;
         //if (era=="2016" && badMuonFilter_ < 0.5) continue;
         //if (era=="2016" && duplicateMuonFilter_ < 0.5) continue; 

         // triggered?  ================================================================================================================================================    
         isMu23 = TriggerMatching(analysisTree, analysisTree.muon_eta[muonIndex], analysisTree.muon_phi[muonIndex], nHighPtLegMuon, deltaRTrigMatch) && analysisTree.muon_pt[muonIndex]>ptMuonHighCut;
         isMu8 = TriggerMatching(analysisTree, analysisTree.muon_eta[muonIndex], analysisTree.muon_phi[muonIndex], nLowPtLegMuon, deltaRTrigMatch) && analysisTree.muon_pt[muonIndex]>ptMuonLowCut;
         if (applyDzFilterMatch){
            isMu23dz = TriggerMatching(analysisTree, analysisTree.muon_eta[muonIndex], analysisTree.muon_phi[muonIndex], nMu23Ele12DzFilter, deltaRTrigMatch);
            isMu8dz = TriggerMatching(analysisTree, analysisTree.muon_eta[muonIndex], analysisTree.muon_phi[muonIndex], nMu8Ele23DzFilter, deltaRTrigMatch);
            isMu23 = isMu23 && isMu23dz;
            isMu8 = isMu8 && isMu8dz;
         }
         isEle23 = TriggerMatching(analysisTree, analysisTree.electron_eta[electronIndex], analysisTree.electron_phi[electronIndex], nHighPtLegElectron, deltaRTrigMatch) && analysisTree.electron_pt[electronIndex]>ptElectronHighCut;
         isEle12 = TriggerMatching(analysisTree, analysisTree.electron_eta[electronIndex], analysisTree.electron_phi[electronIndex], nLowPtLegElectron, deltaRTrigMatch) && analysisTree.electron_pt[electronIndex]>ptElectronLowCut; 
         if (applyDzFilterMatch) {
         isEle12dz = TriggerMatching(analysisTree, analysisTree.electron_eta[electronIndex], analysisTree.electron_phi[electronIndex], nMu23Ele12DzFilter, deltaRTrigMatch);
         isEle23dz = TriggerMatching(analysisTree, analysisTree.electron_eta[electronIndex], analysisTree.electron_phi[electronIndex], nMu8Ele23DzFilter, deltaRTrigMatch);
         isEle23 = isEle23 && isEle23dz;
         isEle12 = isEle12 && isEle12dz;
         }		      
         trg_muonelectron = 
            (isMu23&&isEle12&&analysisTree.muon_pt[muonIndex]>ptMuonHighCut) ||
            (isMu8&&isEle23&&analysisTree.electron_pt[electronIndex]>ptElectronHighCut);
         if (!sync && applyTriggerMatch&&trg_muonelectron<0.5) continue;

         // check charge of electron & muon pair =======================================================================================================================
         os = (analysisTree.muon_charge[muonIndex]*analysisTree.electron_charge[electronIndex]) < 0;

         // looking for extra leptons  =================================================================================================================================
         bool foundExtraElectron = false;
         bool foundExtraMuon = false;
         foundExtraElectron = ElectronVeto(analysisTree, electronIndex, ptVetoElectronCut, etaVetoElectronCut, dxyVetoElectronCut, dzVetoElectronCut, era, applyVetoElectronId, isElectronIsoR03, isoVetoElectronCut);
         foundExtraMuon = MuonVeto(analysisTree, muonIndex, ptVetoMuonCut, etaVetoMuonCut, dxyVetoMuonCut, dzVetoMuonCut, applyVetoMuonId, isMuonIsoR03, isoVetoMuonCut);
         extraelec_veto = foundExtraElectron;
         extramuon_veto = foundExtraMuon;

         // filling muon variables   ===================================================================================================================================
         pt_2 = analysisTree.muon_pt[muonIndex];
         eta_2 = analysisTree.muon_eta[muonIndex];
         phi_2 = analysisTree.muon_phi[muonIndex];
         q_2 = -1;
         if (analysisTree.muon_charge[muonIndex]>0)
            q_2 = 1;
         mva_2 = -10;
         d0_2 = analysisTree.muon_dxy[muonIndex];
         dZ_2 = analysisTree.muon_dz[muonIndex];
         iso_2 = isoMuMin;
         m_2 =  classic_svFit::muonMass;

         // filling electron variables    ==============================================================================================================================
         pt_1 = analysisTree.electron_pt[electronIndex];
         eta_1 = analysisTree.electron_eta[electronIndex];
         phi_1 = analysisTree.electron_phi[electronIndex];
         q_1 = -1;
         if (analysisTree.electron_charge[electronIndex]>0)
            q_1 = 1;
         //mva_1 = analysisTree.electron_mva_id_nontrigPhys14[electronIndex];
         d0_1 = analysisTree.electron_dxy[electronIndex];
         dZ_1 = analysisTree.electron_dz[electronIndex];
         iso_1 = isoEleMin;
         m_1 =  classic_svFit::electronMass;

         // set iso, id and trigger weights   ==========================================================================================================================
         if (!isData || isEmbedded) {

            correctionWS->var("e_pt")->setVal(pt_1);
            correctionWS->var("e_eta")->setVal(eta_1);
            correctionWS->var("e_iso")->setVal(iso_1);
            correctionWS->var("m_pt")->setVal(pt_2);
            correctionWS->var("m_eta")->setVal(eta_2);
            correctionWS->var("m_iso")->setVal(iso_2);
            // scale factors
            if (era=="2016"){
               if (isEmbedded) {
                  isoweight_1 = correctionWS->function("e_idiso_ratio_emb")->getVal();
                  isoweight_2 = correctionWS->function("m_idiso_ratio_emb")->getVal();
               }
               else {
                  isoweight_1 = correctionWS->function("e_idiso_ratio")->getVal();
                  isoweight_2 = correctionWS->function("m_idiso_ratio")->getVal();
               }
            }
            else{
               if (isEmbedded) {
                  isoweight_1 = correctionWS->function("e_id90iso_binned_embed_kit_ratio")->getVal();
                  isoweight_2 = correctionWS->function("m_idiso_binned_embed_kit_ratio")->getVal();
               }
               else {
                  isoweight_1 = correctionWS->function("e_id90iso_binned_kit_ratio")->getVal();
                  isoweight_2 = correctionWS->function("m_idiso_binned_kit_ratio")->getVal();
               }
            }
            correctionWS->var("e_pt")->setVal(pt_1);
            correctionWS->var("e_eta")->setVal(eta_1);
            if (era=="2017") idweight_1 = correctionWS->function("e_trk_ratio")->getVal();
            else idweight_1 =1.0;
            correctionWS->var("m_eta")->setVal(eta_2);
            correctionWS->var("m_pt")->setVal(pt_2);
            if (era=="2016") idweight_2 = correctionWS->function("m_trk_ratio")->getVal();
            else idweight_2 =1.0;
            
            isoweight_1 *= idweight_1;
            isoweight_2 *= idweight_2;

            //		cout << "isoweight_1 = " << isoweight_1
            //		     << "isoweight_2 = " << isoweight_2 << endl;
            
            correctionWS->var("e_pt")->setVal(pt_1);
            correctionWS->var("e_eta")->setVal(eta_1);
            correctionWS->var("m_pt")->setVal(pt_2);
            correctionWS->var("m_eta")->setVal(eta_2);
            correctionWS->var("m_iso")->setVal(iso_2);
            correctionWS->var("e_iso")->setVal(iso_1);
            
            Ele23EffData = correctionWS->function("e_trg_23_binned_ic_data")->getVal();
            Ele12EffData = correctionWS->function("e_trg_12_binned_ic_data")->getVal(); 
            Mu23EffData  = correctionWS->function("m_trg_23_binned_ic_data")->getVal();
            Mu8EffData   = correctionWS->function("m_trg_8_binned_ic_data")->getVal();
            
            if (isEmbedded){
               Ele23EffData = correctionWS->function("e_trg_23_binned_ic_data")->getVal();
               Ele12EffData = correctionWS->function("e_trg_12_binned_ic_data")->getVal(); 
               Mu23EffData  = correctionWS->function("m_trg_23_binned_ic_data")->getVal();
               Mu8EffData   = correctionWS->function("m_trg_8_binned_ic_data")->getVal();
            }
            float trigWeightData = Mu23EffData*Ele12EffData + Mu8EffData*Ele23EffData - Mu23EffData*Ele23EffData;

            if (applyTriggerMatch && !isData) {
               if (isEmbedded) {
                  Ele23EffMC = correctionWS->function("e_trg_23_binned_ic_embed")->getVal();
                  Ele12EffMC = correctionWS->function("e_trg_12_binned_ic_embed")->getVal(); 
                  Mu23EffMC  = correctionWS->function("m_trg_23_binned_ic_embed")->getVal();
                  Mu8EffMC   = correctionWS->function("m_trg_8_binned_ic_embed" )->getVal();
               }
               else {
                  Ele23EffMC = correctionWS->function("e_trg_23_binned_ic_mc")->getVal();
                  Ele12EffMC = correctionWS->function("e_trg_12_binned_ic_mc")->getVal();
                  Mu23EffMC  = correctionWS->function("m_trg_23_binned_ic_mc")->getVal();
                  Mu8EffMC   = correctionWS->function("m_trg_8_binned_ic_mc" )->getVal();
               }
               
               float trigWeightMC   = Mu23EffMC*Ele12EffMC     + Mu8EffMC*Ele23EffMC     - Mu23EffMC*Ele23EffMC;
               if (trigWeightMC>1e-6)
                  trigweight = trigWeightData / trigWeightMC;
               
            }
            else {
               trigweight = trigWeightData;
            }
            
            effweight = trigweight;
         }
         // set properties of dilepton system ==========================================================================================================================
         TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
                                               analysisTree.muon_py[muonIndex],
                                               analysisTree.muon_pz[muonIndex],
                                               classic_svFit::muonMass);
         
         TLorentzVector muonUpLV; muonUpLV.SetXYZM((1.0+muonScale)*analysisTree.muon_px[muonIndex],
                                                   (1.0+muonScale)*analysisTree.muon_py[muonIndex],
                                                   (1.0+muonScale)*analysisTree.muon_pz[muonIndex],
                                                   classic_svFit::muonMass);
         
         TLorentzVector muonDownLV; muonDownLV.SetXYZM((1.0-muonScale)*analysisTree.muon_px[muonIndex],
                                                       (1.0-muonScale)*analysisTree.muon_py[muonIndex],
                                                       (1.0-muonScale)*analysisTree.muon_pz[muonIndex],
                                                       classic_svFit::muonMass);
         
         TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
                                                       analysisTree.electron_py[electronIndex],
                                                       analysisTree.electron_pz[electronIndex],
                                                       classic_svFit::electronMass);
         
         TLorentzVector electronUpLV; electronUpLV.SetXYZM(analysisTree.electron_px_energyscale_up[electronIndex],
                                                           analysisTree.electron_py_energyscale_up[electronIndex],
                                                           analysisTree.electron_pz_energyscale_up[electronIndex],
                                                           classic_svFit::electronMass);
         
         TLorentzVector electronDownLV; electronDownLV.SetXYZM(analysisTree.electron_px_energyscale_down[electronIndex],
                                                               analysisTree.electron_py_energyscale_down[electronIndex],
                                                               analysisTree.electron_pz_energyscale_down[electronIndex],
                                                               classic_svFit::electronMass);

         TLorentzVector electronResoUpLV; electronResoUpLV.SetXYZM(analysisTree.electron_px_energysigma_up[electronIndex],
                                                           analysisTree.electron_py_energysigma_up[electronIndex],
                                                           analysisTree.electron_pz_energysigma_up[electronIndex],
                                                           classic_svFit::electronMass);
         
         TLorentzVector electronResoDownLV; electronResoDownLV.SetXYZM(analysisTree.electron_px_energysigma_down[electronIndex],
                                                               analysisTree.electron_py_energysigma_down[electronIndex],
                                                               analysisTree.electron_pz_energysigma_down[electronIndex],
                                                               classic_svFit::electronMass);
         
         TLorentzVector dileptonLV = muonLV + electronLV;
            
         m_vis = dileptonLV.M();
         pt_vis = dileptonLV.Pt();
         dphi_tt = dPhiFrom2P(muonLV.Px(),muonLV.Py(),electronLV.Px(),electronLV.Py());   
         dr_tt = deltaR(muonLV.Eta(),muonLV.Phi(),electronLV.Eta(),electronLV.Phi());

         // counting jets ==============================================================================================================================================
         TLorentzVector metLV;
         float met_x;
         float met_y;
         
         if (!usePuppiMet){
            met_x = analysisTree.pfmetcorr_ex;
            met_y = analysisTree.pfmetcorr_ey;
         }
         else{
            met_x = analysisTree.puppimet_ex;
            met_y = analysisTree.puppimet_ey;
         }
         met = TMath::Sqrt(met_x*met_x + met_y*met_y);
         metLV.SetXYZT(met_x,met_y,0.,met);

        for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
            
            double absJetEta = fabs(analysisTree.pfjet_eta[jet]);
            double jetEta = analysisTree.pfjet_eta[jet];
            if (absJetEta>jetEtaCut) continue;
            
            float jetPt = analysisTree.pfjet_pt[jet];
            
            jetLV.SetPxPyPzE(analysisTree.pfjet_px[jet],
                             analysisTree.pfjet_py[jet],
                             analysisTree.pfjet_pz[jet],
                             analysisTree.pfjet_e[jet]);
            
            map<TString,TLorentzVector> jetLV_jecUnc;
            
            // Include variations for jec uncertainties
            for (auto uncer_split : jec_unc_map) {
               float sum_unc   = 0;
               for (auto single_jec_unc : uncer_split.second){
                  JetCorrectionUncertainty *unc = single_jec_unc;
                  unc->setJetPt(jetPt);
                  unc->setJetEta(jetEta);
                  double unc_ = unc->getUncertainty(true);
                  sum_unc  += pow(unc_,2);
               }
               float unc_total  = TMath::Sqrt(sum_unc);
               jetLV_jecUnc[uncer_split.first + "Up"]   = jetLV * ( 1 + unc_total);
               jetLV_jecUnc[uncer_split.first + "Down"] = jetLV * ( 1 - unc_total);
               // Propagate jec uncertainties to met
               if( metLV_jecUnc.find(uncer_split.first+"Up") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Up"] = metLV;
               if( metLV_jecUnc.find(uncer_split.first+"Down") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Down"] = metLV;
               if (!isSampleForRecoilCorrection) {
                  metLV_jecUnc[uncer_split.first + "Up"]   -= jetLV* unc_total;
                  metLV_jecUnc[uncer_split.first + "Down"] += jetLV* unc_total;
               }
            }

            float jetPt_tocheck = jetPt;
            if (sync) jetPt_tocheck = jetPt;
            else{
               float jetPtmax = jetPt;
               if(jetLV_jecUnc.at("jecUncEta0To5Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To5Down").Pt();
               if(jetLV_jecUnc.at("jecUncEta0To5Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To5Up").Pt();
               if(jetLV_jecUnc.at("jecUncEta0To3Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To3Down").Pt();
               if(jetLV_jecUnc.at("jecUncEta0To3Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To3Up").Pt();
               if(jetLV_jecUnc.at("jecUncEta3To5Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta3To5Down").Pt();
               if(jetLV_jecUnc.at("jecUncEta3To5Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta3To5Up").Pt();
               if(jetLV_jecUnc.at("jecUncRelativeBalDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeBalDown").Pt();
               if(jetLV_jecUnc.at("jecUncRelativeBalUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeBalUp").Pt();
               if(jetLV_jecUnc.at("jecUncEC2Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2Down").Pt();
               if(jetLV_jecUnc.at("jecUncEC2Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2Up").Pt();
               if(jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt();
               if(jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt();
               if(jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt();
               if(jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt();
               if(jetLV_jecUnc.at("jecUncHFYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncHFYearDown").Pt();
               if(jetLV_jecUnc.at("jecUncHFYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncHFYearUp").Pt();
               if(jetLV_jecUnc.at("jecUncEC2YearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2YearDown").Pt();
               if(jetLV_jecUnc.at("jecUncEC2YearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2YearUp").Pt();
               if(jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt();
               if(jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt();
               if(jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt();
               if(jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt();
               jetPt_tocheck = jetPtmax;
            }
            if (jetPt_tocheck<jetPtLowCut) continue;   

            bool isPFJetId =false;
            if (era=="2016") isPFJetId= looseJetiD_2016(analysisTree,int(jet));
            else if (era=="2017") isPFJetId = tightJetiD_2017(analysisTree,int(jet));
            else if (era=="2018") isPFJetId = tightJetiD_2018(analysisTree,int(jet));
            if (!isPFJetId) continue;

            if (era=="2017" && jetPt < 50 && absJetEta > 2.65 && absJetEta < 3.139) continue;

            bool cleanedJet = true;
            float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                               eta_1,phi_1);
            if (dR1<dRJetLeptonCut) cleanedJet = false;
            float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                               eta_2,phi_2);
            if (dR2<dRJetLeptonCut) cleanedJet = false;
            if (!cleanedJet) continue;

            if (jetPt>jetPtLowCut) jetspt20.push_back(jet);
            
            if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

               bool tagged = false;
               bool tagged_mistagUp = false;
               bool tagged_mistagDown = false;
               bool tagged_btagUp = false;
               bool tagged_btagDown = false;
               tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3]>btagCut;
               tagged_mistagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
               tagged_mistagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
               tagged_btagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
               tagged_btagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;          
               bool taggedRaw = tagged;

               if (!isData) {
                  int flavor = abs(analysisTree.pfjet_flavour[jet]);
                  
                  double jet_scalefactor      = 1;
                  double jet_scalefactor_up   = 1;
                  double jet_scalefactor_down = 1;
                  double JetPtForBTag = jetPt;
                  double tageff = 1;

                  if (flavor==5) {
                     if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                     if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                     jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                     jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                     jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                     tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                  }
                  else if (flavor==4) {
                     if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                     if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                     jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                     jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                     jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                     tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                  }
                  else {
                     if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                     if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                     jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                     jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                     jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                     tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                  }
                  if (tageff<1e-5)      tageff = 1e-5;
                  if (tageff>0.99999)   tageff = 0.99999;
                  r.SetSeed((int)((jetEta+5)*100000));
                  double rannum = r.Rndm();
                  
                  if (tagged) { // demote
                     if(jet_scalefactor<1){
                        double fraction = 1-jet_scalefactor;
                        if (rannum<fraction) tagged = false;
                     }
                     if(jet_scalefactor_up<1){
                        double fraction_up = 1-jet_scalefactor_up;
                        if (rannum<fraction_up) tagged_mistagUp = false;
                     }
                     if(jet_scalefactor_down<1){
                        double fraction_down = 1-jet_scalefactor_down;
                        if (rannum<fraction_down) tagged_mistagDown = false;
                     }
                     tagged_btagUp   = tagged;
                     tagged_btagDown = tagged;
                  }
                     else if (!tagged) { // promote
                        if(jet_scalefactor>1){
                           double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                           if (rannum<fraction) tagged = true;
                        }
                        if(jet_scalefactor_up>1){
                           double fraction_up = (jet_scalefactor_up-1.0)/(1.0/tageff-1.0);
                           if (rannum<fraction_up) tagged_btagUp = true;
                        }
                        if(jet_scalefactor_down>1){
                           double fraction_down = (jet_scalefactor_down-1.0)/(1.0/tageff-1.0);
                           if (rannum<fraction_down) tagged_btagDown = true;
                        }
                        tagged_mistagUp   = tagged;
                        tagged_mistagDown = tagged;
                     }
               }
               
               if (taggedRaw) bjetsRaw.push_back(jet);
               
               if (tagged) {
                  bjets.push_back(jet);
                     if (jetPt>ptLeadingBJet) {
                        ptLeadingBJet = jetPt;
                        indexLeadingBJet = jet;
                     }
               }
               
                  if(tagged_mistagUp)   bjets_mistagUp.push_back(jet);
                  if(tagged_mistagDown) bjets_mistagDown.push_back(jet);
                  if(tagged_btagUp)     bjets_btagUp.push_back(jet);
                  if(tagged_btagDown)   bjets_btagDown.push_back(jet);
            }
            if (jetLV_jecUnc.at("jecUncEta0To5Up").Pt()>jetPtHighCut) njets_jecUncEta0To5Up += 1;
            if (jetLV_jecUnc.at("jecUncEta0To5Down").Pt()>jetPtHighCut) njets_jecUncEta0To5Down += 1;
            if (jetLV_jecUnc.at("jecUncEta0To3Up").Pt()>jetPtHighCut) njets_jecUncEta0To3Up += 1;
            if (jetLV_jecUnc.at("jecUncEta0To3Down").Pt()>jetPtHighCut) njets_jecUncEta0To3Down += 1;
            if (jetLV_jecUnc.at("jecUncEta3To5Up").Pt()>jetPtHighCut) njets_jecUncEta3To5Up += 1;
            if (jetLV_jecUnc.at("jecUncEta3To5Down").Pt()>jetPtHighCut) njets_jecUncEta3To5Down += 1;
            if (jetLV_jecUnc.at("jecUncRelativeBalUp").Pt()>jetPtHighCut) njets_jecUncRelativeBalUp += 1;
            if (jetLV_jecUnc.at("jecUncRelativeBalDown").Pt()>jetPtHighCut) njets_jecUncRelativeBalDown += 1;
            if (jetLV_jecUnc.at("jecUncEC2Up").Pt()>jetPtHighCut) njets_jecUncEC2Up += 1;
            if (jetLV_jecUnc.at("jecUncEC2Down").Pt()>jetPtHighCut) njets_jecUncEC2Down += 1;
            if (jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt()>jetPtHighCut) njets_jecUncFlavorQCDUp += 1;
            if (jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt()>jetPtHighCut) njets_jecUncFlavorQCDDown += 1;
            if (jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt()>jetPtHighCut) njets_jecUncRelativeSampleYearUp += 1;
            if (jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt()>jetPtHighCut) njets_jecUncRelativeSampleYearDown += 1;
            if (jetLV_jecUnc.at("jecUncEC2YearUp").Pt()>jetPtHighCut) njets_jecUncEC2YearUp += 1;
            if (jetLV_jecUnc.at("jecUncEC2YearDown").Pt()>jetPtHighCut) njets_jecUncEC2YearDown += 1;
            if (jetLV_jecUnc.at("jecUncHFYearUp").Pt()>jetPtHighCut) njets_jecUncHFYearUp += 1;
            if (jetLV_jecUnc.at("jecUncHFYearDown").Pt()>jetPtHighCut) njets_jecUncHFYearDown += 1;
            if (jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt()>jetPtHighCut) njets_jecUncAbsoluteYearUp += 1;
            if (jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt()>jetPtHighCut) njets_jecUncAbsoluteYearDown += 1;
            if (jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt()>jetPtHighCut) njets_jecUncBBEC1YearUp += 1;
            if (jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt()>jetPtHighCut) njets_jecUncBBEC1YearDown += 1;
            if (jetPt>jetPtHighCut)
               jets.push_back(jet); 
            
            //if (jetPt_tocheck<jetPtHighCut) continue; cut is done later to make sure that the uncertainties are treated correcly

            if (indexLeadingJet>=0) {
               if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
                  indexSubLeadingJet = jet;
                  ptSubLeadingJet = jetPt;
               }
            }
            if (jetPt>ptLeadingJet) {
               indexSubLeadingJet = indexLeadingJet;
               ptSubLeadingJet = ptLeadingJet;
               indexLeadingJet = jet;
               ptLeadingJet = jetPt;
            }
         }
         njets = jets.size();
         int njetsMax = njets;
         
         njetspt20 = jetspt20.size();
         nbtag = bjets.size();
         nbtag_mistagUp   = bjets_mistagUp.size();
         nbtag_mistagDown = bjets_mistagDown.size();
         nbtag_btagUp   = bjets_btagUp.size();
         nbtag_btagDown = bjets_btagDown.size();
         nbtag_noSF = bjetsRaw.size();
         
         if (indexLeadingBJet>=0) {
            bpt = analysisTree.pfjet_pt[indexLeadingBJet];
            beta_1 = analysisTree.pfjet_eta[indexLeadingBJet];
            bphi = analysisTree.pfjet_phi[indexLeadingBJet];
         }
         
         if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
            cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
        
         if (indexLeadingJet>=0) {
            jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
            jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
            jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
            jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
            jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
                            analysisTree.pfjet_py[indexLeadingJet],
                            analysisTree.pfjet_pz[indexLeadingJet],
                            analysisTree.pfjet_e[indexLeadingJet]);

            for (auto uncer_split : jec_unc_map) {
               float sum_unc   = 0;
               for (auto single_jec_unc : uncer_split.second){
                  JetCorrectionUncertainty *unc = single_jec_unc;
                  unc->setJetPt(jet1.Pt());
                  unc->setJetEta(jet1.Eta());
                  double unc_ = unc->getUncertainty(true);
                  sum_unc  += pow(unc_,2);
               }
               float unc_total = TMath::Sqrt(sum_unc);
               jet1LV_jecUnc[uncer_split.first+"Up"]   = jet1*(1+unc_total);
               jet1LV_jecUnc[uncer_split.first+"Down"] = jet1*(1-unc_total);
            }
         }

         if (indexSubLeadingJet>=0) {
            jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
            jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
            jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
            jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];

            jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
                            analysisTree.pfjet_py[indexSubLeadingJet],
                            analysisTree.pfjet_pz[indexSubLeadingJet],
                            analysisTree.pfjet_e[indexSubLeadingJet]);
            
            for (auto uncer_split : jec_unc_map) {
               float sum_unc   = 0;
               for (auto single_jec_unc : uncer_split.second){
                  JetCorrectionUncertainty *unc = single_jec_unc;
                  unc->setJetPt(jet2.Pt());
                  unc->setJetEta(jet2.Eta());
                  double unc_ = unc->getUncertainty(true);
                  sum_unc  += pow(unc_,2);
               }
               float unc_total = TMath::Sqrt(sum_unc);
               jet2LV_jecUnc[uncer_split.first+"Up"]   = jet2*(1+unc_total);
               jet2LV_jecUnc[uncer_split.first+"Down"] = jet2*(1-unc_total);
            }
         }

         if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
            
            mjj      = (jet1+jet2).M();
            dijetpt  = (jet1+jet2).Pt();
            dijetphi = (jet1+jet2).Phi();
            jdeta = fabs(analysisTree.pfjet_eta[indexLeadingJet]-analysisTree.pfjet_eta[indexSubLeadingJet]);
         }

         // METs =======================================================================================================================================================
         float met_x_recoilscaleUp;
         float met_x_recoilscaleDown;
         float met_y_recoilscaleUp;
         float met_y_recoilscaleDown;
         float met_x_recoilresoUp;
         float met_x_recoilresoDown;
         float met_y_recoilresoUp;
         float met_y_recoilresoDown;
         
         float met_unclMetUp_x;
         float met_unclMetUp_y;
         float met_unclMetDown_x;
         float met_unclMetDown_y;
        
         if (!usePuppiMet){
            met_x_recoilscaleUp = analysisTree.pfmetcorr_ex;
            met_x_recoilscaleDown = analysisTree.pfmetcorr_ex;
            met_y_recoilscaleUp = analysisTree.pfmetcorr_ey;
            met_y_recoilscaleDown = analysisTree.pfmetcorr_ey;
            met_x_recoilresoUp = analysisTree.pfmetcorr_ex;
            met_x_recoilresoDown = analysisTree.pfmetcorr_ex;
            met_y_recoilresoUp = analysisTree.pfmetcorr_ey;
            met_y_recoilresoDown = analysisTree.pfmetcorr_ey;
         
            met_unclMetUp_x    = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            met_unclMetUp_y    = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            met_unclMetDown_x  = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            met_unclMetDown_y  = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
         }
         else{
            met_x_recoilscaleUp = analysisTree.puppimet_ex;
            met_x_recoilscaleDown = analysisTree.puppimet_ex;
            met_y_recoilscaleUp = analysisTree.puppimet_ey;
            met_y_recoilscaleDown = analysisTree.puppimet_ey;
            met_x_recoilresoUp = analysisTree.puppimet_ex;
            met_x_recoilresoDown = analysisTree.puppimet_ex;
            met_y_recoilresoUp = analysisTree.puppimet_ey;
            met_y_recoilresoDown = analysisTree.puppimet_ey;
         
            met_unclMetUp_x    = analysisTree.puppimet_ex_UnclusteredEnUp;
            met_unclMetUp_y    = analysisTree.puppimet_ey_UnclusteredEnUp;
            met_unclMetDown_x  = analysisTree.puppimet_ex_UnclusteredEnDown;
            met_unclMetDown_y  = analysisTree.puppimet_ey_UnclusteredEnDown;
         }
         met = TMath::Sqrt(met_x*met_x + met_y*met_y);
         metphi = TMath::ATan2(met_y,met_x);

         if (!usePuppiMet){
            metcov00 = analysisTree.pfmetcorr_sigxx;
            metcov01 = analysisTree.pfmetcorr_sigxy;
            metcov10 = analysisTree.pfmetcorr_sigyx;
            metcov11 = analysisTree.pfmetcorr_sigyy;
         }
         else{
            metcov00 = analysisTree.puppimet_sigxx;
            metcov01 = analysisTree.puppimet_sigxy;
            metcov10 = analysisTree.puppimet_sigyx;
            metcov11 = analysisTree.puppimet_sigyy;
         }
         // recoil corrections =========================================================================================================================================
         int njetsforrecoil = njets;
         if (isW) njetsforrecoil = njets + 1;
         
         met_uncorr = met;
         metphi_uncorr = metphi;
         
         float pfmet_corr_x = met_x;
         float pfmet_corr_y = met_y;
         if ((isW||isDY||isSignal)&&!isData) {
            if (applySimpleRecoilCorrections) {
               if (!usePuppiMet){
                  recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
               }
               else{
                  recoilMetCorrectorPuppi.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
               }
            }
            else {
               if (!usePuppiMet)recoilMetCorrector.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
               else recoilMetCorrectorPuppi.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
            }
         }
         
         met_x = pfmet_corr_x;
         met_y = pfmet_corr_y;
         met = TMath::Sqrt(met_x*met_x+met_y*met_y);
         metphi = TMath::ATan2(met_y,met_x);
         // update metLV (regognize changes if present)
         metLV.SetXYZT(met_x,met_y,0.,met);
         met_recoilscaleUp   = TMath::Sqrt(met_x_recoilscaleUp*met_x_recoilscaleUp+met_y_recoilscaleUp*met_y_recoilscaleUp);
         met_recoilscaleDown = TMath::Sqrt(met_x_recoilscaleDown*met_x_recoilscaleDown+met_y_recoilscaleDown*met_y_recoilscaleDown);
         met_recoilresoUp    = TMath::Sqrt(met_x_recoilresoUp*met_x_recoilresoUp+met_y_recoilresoUp*met_y_recoilresoUp);
         met_recoilresoDown  = TMath::Sqrt(met_x_recoilresoDown*met_x_recoilresoDown+met_y_recoilresoDown*met_y_recoilresoDown);
         TLorentzVector metLV_recoilscaleUp; metLV_recoilscaleUp.SetXYZT(met_x_recoilscaleUp, met_y_recoilscaleUp, 0., met_recoilscaleUp);
         TLorentzVector metLV_recoilscaleDown; metLV_recoilscaleDown.SetXYZT(met_x_recoilscaleDown, met_y_recoilscaleDown, 0., met_recoilscaleDown);
         TLorentzVector metLV_recoilresoUp; metLV_recoilresoUp.SetXYZT(met_x_recoilresoUp, met_y_recoilresoUp, 0., met_recoilresoUp);
         TLorentzVector metLV_recoilresoDown; metLV_recoilresoDown.SetXYZT(met_x_recoilresoDown, met_y_recoilresoDown, 0., met_recoilresoDown);
         
         met_unclMetUp = TMath::Sqrt(met_unclMetUp_x*met_unclMetUp_x+met_unclMetUp_y*met_unclMetUp_y);
         metphi_unclMetUp = TMath::ATan2(met_unclMetUp_y,met_unclMetUp_x);
         met_unclMetDown = TMath::Sqrt(met_unclMetDown_x*met_unclMetDown_x+met_unclMetDown_y*met_unclMetDown_y);
         metphi_unclMetDown = TMath::ATan2(met_unclMetDown_y,met_unclMetDown_x);
         
         if(isSampleForRecoilCorrection){
            met_unclMetUp_x   = met_x;
            met_unclMetDown_x = met_x;
            met_unclMetUp_y   = met_y;
            met_unclMetDown_y = met_y;
            met_unclMetUp     = met;
            met_unclMetDown   = met;
            }
         
         float genmet_ex = analysisTree.genmet_ex;
         float genmet_ey = analysisTree.genmet_ey;
         
         genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
         genmetphi = TMath::ATan2(genmet_ey,genmet_ex);
         
         // bisector of electron and muon transverse momenta
         float electronUnitX = electronLV.Px()/electronLV.Pt();
         float electronUnitY = electronLV.Py()/electronLV.Pt();
         float muonUnitX = muonLV.Px()/muonLV.Pt();
         float muonUnitY = muonLV.Py()/muonLV.Pt();
         float zetaX = electronUnitX + muonUnitX;
         float zetaY = electronUnitY + muonUnitY;
         float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
         
         zetaX = zetaX/normZeta;
         zetaY = zetaY/normZeta;
         
         float vectorX = met_x + muonLV.Px() + electronLV.Px();
         float vectorY = met_y + muonLV.Py() + electronLV.Py();
         float vectorVisX = muonLV.Px() + electronLV.Px();
         float vectorVisY = muonLV.Py() + electronLV.Py();
         pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
         
         double px_escaleUp = analysisTree.electron_px_energyscale_up[electronIndex];
         double py_escaleUp = analysisTree.electron_py_energyscale_up[electronIndex];
         double px_eresoUp = analysisTree.electron_px_energysigma_up[electronIndex];
         double py_eresoUp = analysisTree.electron_py_energysigma_up[electronIndex];

         double px_e = pt_1 * TMath::Cos(phi_1);
         double py_e = pt_1 * TMath::Sin(phi_1);
         
         double px_escaleDown = analysisTree.electron_px_energyscale_down[electronIndex];
         double py_escaleDown = analysisTree.electron_py_energyscale_down[electronIndex];
         double px_eresoDown = analysisTree.electron_px_energysigma_down[electronIndex];
         double py_eresoDown = analysisTree.electron_py_energysigma_down[electronIndex];
         
         double metx_escaleUp = met_x + px_e - px_escaleUp;
         double mety_escaleUp = met_y + py_e - py_escaleUp;
         double metx_escaleDown = met_x + px_e - px_escaleDown;
         double mety_escaleDown = met_y + py_e - py_escaleDown;

         double metx_eresoUp = met_x + px_e - px_eresoUp;
         double mety_eresoUp = met_y + py_e - py_eresoUp;
         double metx_eresoDown = met_x + px_e - px_eresoDown;
         double mety_eresoDown = met_y + py_e - py_eresoDown;
         
         // computation of DZeta variable
         // pfmet
         computeDzeta(met_x,met_y,
                      zetaX,zetaY,pzetavis,pzetamiss,dzeta);
         
         // genmet
         computeDzeta(analysisTree.genmet_ex,analysisTree.genmet_ey,
                      zetaX,zetaY,pzetavis,pzetamiss_genmet,dzeta_genmet);
         
         float met_escaleUp = TMath::Sqrt(metx_escaleUp*metx_escaleUp+mety_escaleUp*mety_escaleUp);
         float met_escaleDown = TMath::Sqrt(metx_escaleDown*metx_escaleDown+mety_escaleDown*mety_escaleDown);
         
         TLorentzVector metLV_escaleUp; metLV_escaleUp.SetXYZT(metx_escaleUp,mety_escaleUp,0.,met_escaleUp);
         TLorentzVector metLV_escaleDown; metLV_escaleDown.SetXYZT(metx_escaleDown,mety_escaleDown,0.,met_escaleDown);

         float met_eresoUp = TMath::Sqrt(metx_eresoUp*metx_eresoUp+mety_eresoUp*mety_eresoUp);
         float met_eresoDown = TMath::Sqrt(metx_eresoDown*metx_eresoDown+mety_eresoDown*mety_eresoDown);
         
         TLorentzVector metLV_eresoUp; metLV_eresoUp.SetXYZT(metx_eresoUp,mety_eresoUp,0.,met_eresoUp);
         TLorentzVector metLV_eresoDown; metLV_eresoDown.SetXYZT(metx_eresoDown,mety_eresoDown,0.,met_eresoDown);

         mt_1 = mT(electronLV,metLV);
         mt_2 = mT(muonLV,metLV);
         mtmax = TMath::Max(float(mt_1),float(mt_2));
         
         mCDF = (muonLV+electronLV+metLV).M();
         
         TLorentzVector metLV_unclMetUp; metLV_unclMetUp.SetXYZT(met_unclMetUp_x,met_unclMetUp_y,0.,met_unclMetUp);
         TLorentzVector metLV_unclMetDown; metLV_unclMetDown.SetXYZT(met_unclMetDown_x,met_unclMetDown_y,0.,met_unclMetDown);
         // computing total transverse mass
         mTtot = totalTransverseMass        ( muonLV ,     electronLV , metLV);
         mTemu   = mT(electronLV,muonLV);
         dphi_mumet = dPhiFrom2P(muonLV.Px(),muonLV.Py(),
                                 metLV.Px(),metLV.Py());
         dphi_emet  = dPhiFrom2P(electronLV.Px(),electronLV.Py(),
                                 metLV.Px(),metLV.Py());
         
         
         
         mTdileptonMET = mT(dileptonLV,metLV);
         
         pt_ttjj = -10;
         if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
            pt_ttjj = (muonLV+electronLV+metLV+jet1+jet2).Pt();
         }
         pt_tt = (muonLV+electronLV+metLV).Pt();
     
         checkSV = false;
         checkFastMTT = false;
         passesPreSel = false;
         passesFastMTTPreSel = false;

         if (sync) {
            checkSV = computeSVFitMass;
            checkFastMTT = computeFastMTTMass;
         }
         else {
            passesPreSel = iso_1<0.15 && iso_2<0.2 && trg_muonelectron > 0.5 && extraelec_veto<0.5 && extramuon_veto<0.5 && (nbtag==0||nbtag_mistagUp==0||nbtag_mistagDown==0||nbtag_btagUp==0||nbtag_btagDown==0);
            passesFastMTTPreSel = iso_1<0.15 && iso_2<0.2 && trg_muonelectron > 0.5 && extraelec_veto<0.5 && extramuon_veto<0.5;
            checkSV = computeSVFitMass && passesPreSel;
            checkFastMTT = computeFastMTTMass && passesFastMTTPreSel;
         }

         TMatrixD covMET(2, 2);

         if (checkSV || checkFastMTT) {

            if (checkSV && checkFastMTT) {
               std::cout<<"Both SVFit and FastMTT are switched on."<<std::endl;
               exit(-1);
            }
            // covariance matrix MET
            covMET[0][0] =  metcov00;
            covMET[1][0] =  metcov10;
            covMET[0][1] =  metcov01;
            covMET[1][1] =  metcov11;
            
            // define electron 4-vector
            classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt_1, eta_1, phi_1, 0.51100e-3);
            // define muon 4-vector
            classic_svFit::MeasuredTauLepton svFitMu(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3); 
            
            if (checkSV){
               ClassicSVfit algo = SVFitMassComputation(svFitEle, svFitMu, met_x, met_y, covMET, inputFile_visPtResolution);

               if ( !algo.isValidSolution() )
                  std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;

               m_sv =static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass();
               mt_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass(); // return value of transverse svfit mass is in units of GeV
               pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt();
               eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
               phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
            }

            if (checkFastMTT){

               LorentzVector ttP4 = FastMTTComputation(svFitEle, svFitMu, met_x, met_y, covMET, tau1P4, tau2P4);

               m_sv = ttP4.M();
               mt_sv = mT_LorentzVector(tau1P4,tau2P4);

               pt_sv = ttP4.Pt();
               eta_sv = ttP4.Eta();
               phi_sv = ttP4.Phi();
            }
         }
         gen_match_1 = 6;
         gen_match_2 = 6;
         
         float minDR_1 = 0.2;
         float minDR_2 = 0.2;
         unsigned int gen_1 = 0;
         bool bgen_1 = false;
         unsigned int gen_2 = 0;
         bool bgen_2 = false;

         if (!isData) {
            for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
               TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
                                                   analysisTree.genparticles_py[igen],
                                                   analysisTree.genparticles_pz[igen],
                                                   analysisTree.genparticles_e[igen]);
               float ptGen = genLV.Pt();
               bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
               bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
               bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
               bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
               bool isAnyType = type1 || type2 || type3 || type4;
               if (isAnyType && analysisTree.genparticles_status[igen]==1) {
                  float etaGen = genLV.Eta();
                  float phiGen = genLV.Phi();
                  float deltaR_1 = deltaR(eta_1,phi_1,
                                          etaGen,phiGen);
                  if (deltaR_1<minDR_1) {
                     minDR_1 = deltaR_1;
                     gen_1 = igen;
                     bgen_1 = true;
                     if (type1) gen_match_1 = 1;
                     else if (type2) gen_match_1 = 2;
                     else if (type3) gen_match_1 = 3;
                     else if (type4) gen_match_1 = 4;
                  }
                  
                  float deltaR_2 = deltaR(eta_2,phi_2,
                                          etaGen,phiGen);
                  if (deltaR_2<minDR_2) {
                     minDR_2 = deltaR_2;
                     gen_2 = igen;
                     bgen_2 = true;
                     if (type1) gen_match_2 = 1;
                     else if (type2) gen_match_2 = 2;
                     else if (type3) gen_match_2 = 3;
                     else if (type4) gen_match_2 = 4;
                  }
               }
            }
            TLorentzVector genElectronLV; genElectronLV.SetPtEtaPhiM(pt_1,eta_1,phi_1,0.51100e-3);
            TLorentzVector genMuonLV; genMuonLV.SetPtEtaPhiM(pt_2,eta_2,phi_2,105.658e-3);
            if (bgen_1) genElectronLV.SetXYZT(analysisTree.genparticles_px[gen_1],
                                              analysisTree.genparticles_py[gen_1],
                                              analysisTree.genparticles_pz[gen_1],
                                              analysisTree.genparticles_e[gen_1]);
            if (bgen_2) genMuonLV.SetXYZT(analysisTree.genparticles_px[gen_2],
                                          analysisTree.genparticles_py[gen_2],
                                          analysisTree.genparticles_pz[gen_2],
                                          analysisTree.genparticles_e[gen_2]);
            
            //	TLorentzVector genDileptonLV = genElectronLV + genMuonLV;
            TLorentzVector genMetLV; genMetLV.SetXYZM(genmet_ex,genmet_ey,
                                                      0,0);
            
            mTemu_gen = mT(genElectronLV,genMuonLV);
            mTemet_gen = mT(genElectronLV,genMetLV);
            mTmumet_gen = mT(genMuonLV,genMetLV);
            mTtot_gen = totalTransverseMass(genElectronLV,genMuonLV,genMetLV);
            
            if (gen_match_1>2&&gen_match_2>3) { 
               isZTT = true;
               isZLL = false;
            }
            else {
               isZTT = false;
               isZLL = true;
            }
            if (gen_match_1 == 3 && gen_match_2 ==4) veto_embedded = true;
            if (gen_match_1 == 3 && gen_match_2 ==4) isZTTEM = true;
            if (gen_match_1 == 4 && gen_match_2 ==4) isZTTMM = true;
            if (gen_match_1 == 3 && gen_match_2 ==3) isZTTEE = true;
            if (gen_match_1 == 1 && gen_match_2 ==1) isPromptZEE = true;
            if (gen_match_1 == 2 && gen_match_2 ==2) isPromptZMM = true;
            }
         
         // QCD estimation =======================================================================================================================================================

         correctionWS->var("e_pt")->setVal(pt_1);
         correctionWS->var("m_pt")->setVal(pt_2);
         correctionWS->var("njets")->setVal(njets);
         correctionWS->var("dR")->setVal(dr_tt);
         double_t em_qcd_osss_binned = correctionWS->function("em_qcd_osss")->getVal();
         double_t em_qcd_osss_binned_0jet_rate_up = correctionWS->function("em_qcd_osss_stat_0jet_unc1_up")->getVal();
         double_t em_qcd_osss_binned_0jet_rate_down = correctionWS->function("em_qcd_osss_stat_0jet_unc1_down")->getVal();
         double_t em_qcd_osss_binned_1jet_rate_up = correctionWS->function("em_qcd_osss_stat_1jet_unc1_up")->getVal();
         double_t em_qcd_osss_binned_1jet_rate_down = correctionWS->function("em_qcd_osss_stat_1jet_unc1_down")->getVal();
         double_t em_qcd_osss_binned_2jet_rate_up = correctionWS->function("em_qcd_osss_stat_2jet_unc1_up")->getVal();
         double_t em_qcd_osss_binned_2jet_rate_down = correctionWS->function("em_qcd_osss_stat_2jet_unc1_down")->getVal();
         double_t em_qcd_osss_binned_0jet_shape_up = correctionWS->function("em_qcd_osss_stat_0jet_unc2_up")->getVal();
         double_t em_qcd_osss_binned_0jet_shape_down = correctionWS->function("em_qcd_osss_stat_0jet_unc2_down")->getVal();
         double_t em_qcd_osss_binned_1jet_shape_up = correctionWS->function("em_qcd_osss_stat_1jet_unc2_up")->getVal();
         double_t em_qcd_osss_binned_1jet_shape_down = correctionWS->function("em_qcd_osss_stat_1jet_unc2_down")->getVal();
         double_t em_qcd_osss_binned_2jet_shape_up = correctionWS->function("em_qcd_osss_stat_2jet_unc2_up")->getVal();
         double_t em_qcd_osss_binned_2jet_shape_down = correctionWS->function("em_qcd_osss_stat_2jet_unc2_down")->getVal();
         double_t em_qcd_extrap_up = correctionWS->function("em_qcd_osss_extrap_up")->getVal();
         double_t em_qcd_extrap_down = correctionWS->function("em_qcd_osss_extrap_down")->getVal();
        
         qcdweight =  em_qcd_osss_binned;
         qcdweight_0jet_rate_up =  em_qcd_osss_binned_0jet_rate_up;
         qcdweight_0jet_rate_down =  em_qcd_osss_binned_0jet_rate_down;   
         qcdweight_1jet_rate_up =  em_qcd_osss_binned_1jet_rate_up;
         qcdweight_1jet_rate_down =  em_qcd_osss_binned_1jet_rate_down;
         qcdweight_2jet_rate_up =  em_qcd_osss_binned_2jet_rate_up;
         qcdweight_2jet_rate_down =  em_qcd_osss_binned_2jet_rate_down;
         qcdweight_0jet_shape_up =  em_qcd_osss_binned_0jet_shape_up;
         qcdweight_0jet_shape_down =  em_qcd_osss_binned_0jet_shape_down;   
         qcdweight_1jet_shape_up =  em_qcd_osss_binned_1jet_shape_up;
         qcdweight_1jet_shape_down =  em_qcd_osss_binned_1jet_shape_down;
         qcdweight_2jet_shape_up =  em_qcd_osss_binned_2jet_shape_up;
         qcdweight_2jet_shape_down = em_qcd_osss_binned_2jet_shape_down;
         
         qcdweight_iso_up = em_qcd_extrap_up;
         qcdweight_iso_down = em_qcd_extrap_down;

         if (!isData) effweight = effweight*isoweight_1*isoweight_2;
         // if (!isData){
         //    if (era=="2016")
         //       {
         //          double weightE = 1;
         //          double weightEUp = 1;
         //          double weightEDown = 1;
         //          if (gen_match_1==6) {
         //             double dRmin = 0.5;
         //             double ptJet = pt_1;
         //             for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
         //                double etaJet = analysisTree.pfjet_eta[ijet];
         //                double phiJet = analysisTree.pfjet_phi[ijet];
         //                double dRjetE = deltaR(etaJet,phiJet,eta_1,phi_1);
         //                if (dRjetE<dRmin) {
         //                   dRmin = dRjetE;
         //                   ptJet = analysisTree.pfjet_pt[ijet];
         //                }
         //             }
         //             weightE     = a_jetEle + b_jetEle*ptJet;
         //             weightEUp   = a_jetEleUp + b_jetEleUp*ptJet;
         //             weightEDown = a_jetEleDown + b_jetEleDown*ptJet;
         //             fakeweight *= weightE;
         //          }
         //          else {
         //             weightE = isoweight_1;
         //             weightEUp = isoweight_1;
         //             weightEDown = isoweight_1;
         //          }
         //          double weightMu = 1;
         //          double weightMuUp = 1;
         //          double weightMuDown = 1;
         //          if (gen_match_2==6) {
         //             double dRmin = 0.5;
         //             double ptJet = pt_2;
         //             for (unsigned int ijet=0; ijet<analysisTree.pfjet_count; ++ijet) {
         //                double etaJet = analysisTree.pfjet_eta[ijet];
         //                double phiJet = analysisTree.pfjet_phi[ijet];
         //                double dRjetM = deltaR(etaJet,phiJet,eta_2,phi_2);
         //                if (dRjetM<dRmin) {
         //                   dRmin = dRjetM;
         //                   ptJet = analysisTree.pfjet_pt[ijet];
         //                }
         //             }
         //             weightMu     = a_jetMu     + b_jetMu*ptJet;
         //             weightMuUp   = a_jetMuUp   + b_jetMuUp*ptJet;
         //             weightMuDown = a_jetMuDown + b_jetMuDown*ptJet;
         //             fakeweight *= weightMu;
         //          }
         //          else {
         //             weightMu = isoweight_2;
         //             weightMuUp = isoweight_2;
         //             weightMuDown = isoweight_2;
         //          }
         //          float effweight0 = effweight;
         //          effweight = effweight0 * weightE * weightMu;
         //          effweight_jetMuUp   = effweight0 * weightE     * weightMuUp;
         //          effweight_jetMuDown = effweight0 * weightE     * weightMuDown; 
         //          effweight_jetEUp    = effweight0 * weightEUp   * weightMu;
         //          effweight_jetEDown  = effweight0 * weightEDown * weightMu; 
                  
         //          if (sync) weight = effweight * puweight * 0.979;
         //       }
         //    else {
         //       effweight = effweight*isoweight_1*isoweight_2;
         //    }
         // }

         // First set the calibrated variables to the uncalibrated versions -> necessary for data
         d0_1_cal=d0_1;
         dZ_1_cal=dZ_1;
         d0_2_cal=d0_2;
         dZ_2_cal=dZ_2;
         
         // if (!isData){               
         //    //check if particle is prompt or non-prompt
         //         if (gen_match_1 == 3) {  
         //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_1,eta_1,"d0_1",d0_1, d0_1_cal);
         //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_1,eta_1,"dZ_1",dZ_1, dZ_1_cal);
         //            //d0_1_cal = doquantileshift_d01_np.shift(d0_1, 0.55);
         //            //dZ_1_cal = doquantileshift_dZ1_np.shift(dZ_1, 0.55);
                  
         //         }
         //         else if (gen_match_1 == 1){
         //            calibrateIP.DoCalibrationForPromptElectrons(pt_1,eta_1,"d0_1",d0_1, d0_1_cal);
         //            calibrateIP.DoCalibrationForPromptElectrons(pt_1,eta_1,"dZ_1",dZ_1, dZ_1_cal);
         //            //d0_1_cal = doquantileshift_d01_pele.shift(d0_1, 0.55);
         //            //dZ_1_cal = doquantileshift_dZ1_pele.shift(dZ_1, 0.55);
         //         }
         //         else {
         //            d0_1_cal = d0_1;
         //            dZ_1_cal = dZ_1;
         //         }
         //         if (gen_match_2 == 4) {  
         //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_2,eta_2,"d0_2",d0_2,d0_2_cal);
         //            calibrateIP.DoCalibrationForNonPromptLeptons(pt_2,eta_2,"dZ_2",dZ_2,dZ_2_cal);
         //            //d0_2_cal = doquantileshift_d02_np.shift(d0_2, 0.55);
         //            //dZ_2_cal = doquantileshift_dZ2_np.shift(dZ_2, 0.55);
         //         }
         //         else if (gen_match_2 == 2){
         //            calibrateIP.DoCalibrationForPromptMuons(pt_2,eta_2,"d0_2",d0_2, d0_2_cal);
         //            calibrateIP.DoCalibrationForPromptMuons(pt_2,eta_2,"dZ_2",dZ_2, dZ_2_cal);
         //            //d0_2_cal = doquantileshift_d02_pmu.shift(d0_2, 0.55);
         //            //dZ_2_cal = doquantileshift_dZ2_pmu.shift(dZ_2, 0.55);
         //         }
         //         else{
         //            d0_2_cal =d0_2;
         //            dZ_2_cal =dZ_2;
         //         }
         //      }       
         
         // Add for all relevant variables the met uncertainty
         for(auto &uncert : uncertainty_map){
            uncert.second.electronLV = electronLV;
            uncert.second.muonLV     = muonLV;
            uncert.second.metLV      = metLV;
            uncert.second.jet1LV     = jet1;
            uncert.second.jet2LV     = jet2;
         }
         
         uncertainty_map.at("unclMetUp").metLV = metLV_unclMetUp;
         uncertainty_map.at("unclMetDown").metLV = metLV_unclMetDown;
         uncertainty_map.at("escaleUp").metLV = metLV_escaleUp;
         uncertainty_map.at("escaleUp").electronLV = electronUpLV;
         uncertainty_map.at("escaleDown").metLV = metLV_escaleDown;
         uncertainty_map.at("escaleDown").electronLV = electronDownLV;
         uncertainty_map.at("eresoUp").metLV = metLV_eresoUp;
         uncertainty_map.at("eresoUp").electronLV = electronResoUpLV;
         uncertainty_map.at("eresoDown").metLV = metLV_eresoDown;
         uncertainty_map.at("eresoDown").electronLV = electronResoDownLV;
         //uncertainty_map.at("mscaleUp").muonLV = muonUpLV;
         //uncertainty_map.at("mscaleDown").muonLV = muonDownLV;
         uncertainty_map.at("recoilscaleUp").metLV = metLV_recoilscaleUp;
         uncertainty_map.at("recoilscaleDown").metLV = metLV_recoilscaleDown;
         uncertainty_map.at("recoilresoUp").metLV = metLV_recoilresoUp;
         uncertainty_map.at("recoilresoDown").metLV = metLV_recoilresoDown;
  
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Up").metLV = metLV_jecUnc.at("jecUncEta0To5Up");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Up").jet1LV = jet1LV_jecUnc.at("jecUncEta0To5Up");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To5Up").jet2LV = jet2LV_jecUnc.at("jecUncEta0To5Up");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Down").metLV = metLV_jecUnc.at("jecUncEta0To5Down");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To5Down").jet1LV = jet1LV_jecUnc.at("jecUncEta0To5Down");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To5Down").jet2LV = jet2LV_jecUnc.at("jecUncEta0To5Down");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Up").metLV = metLV_jecUnc.at("jecUncEta0To3Up");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Up").jet1LV = jet1LV_jecUnc.at("jecUncEta0To3Up");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To3Up").jet2LV = jet2LV_jecUnc.at("jecUncEta0To3Up");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Down").metLV = metLV_jecUnc.at("jecUncEta0To3Down");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta0To3Down").jet1LV = jet1LV_jecUnc.at("jecUncEta0To3Down");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta0To3Down").jet2LV = jet2LV_jecUnc.at("jecUncEta0To3Down");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Up").metLV = metLV_jecUnc.at("jecUncEta3To5Up");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Up").jet1LV = jet1LV_jecUnc.at("jecUncEta3To5Up");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta3To5Up").jet2LV = jet2LV_jecUnc.at("jecUncEta3To5Up");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Down").metLV = metLV_jecUnc.at("jecUncEta3To5Down");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEta3To5Down").jet1LV = jet1LV_jecUnc.at("jecUncEta3To5Down");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEta3To5Down").jet2LV = jet2LV_jecUnc.at("jecUncEta3To5Down");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalUp").metLV = metLV_jecUnc.at("jecUncRelativeBalUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalUp").jet1LV = jet1LV_jecUnc.at("jecUncRelativeBalUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeBalUp").jet2LV = jet2LV_jecUnc.at("jecUncRelativeBalUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalDown").metLV = metLV_jecUnc.at("jecUncRelativeBalDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeBalDown").jet1LV = jet1LV_jecUnc.at("jecUncRelativeBalDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeBalDown").jet2LV = jet2LV_jecUnc.at("jecUncRelativeBalDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEC2Up").metLV = metLV_jecUnc.at("jecUncEC2Up");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEC2Up").jet1LV = jet1LV_jecUnc.at("jecUncEC2Up");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEC2Up").jet2LV = jet2LV_jecUnc.at("jecUncEC2Up");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEC2Down").metLV = metLV_jecUnc.at("jecUncEC2Down");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEC2Down").jet1LV = jet1LV_jecUnc.at("jecUncEC2Down");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEC2Down").jet2LV = jet2LV_jecUnc.at("jecUncEC2Down");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncFlavorQCDUp").metLV = metLV_jecUnc.at("jecUncFlavorQCDUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncFlavorQCDUp").jet1LV = jet1LV_jecUnc.at("jecUncFlavorQCDUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncFlavorQCDUp").jet2LV = jet2LV_jecUnc.at("jecUncFlavorQCDUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncFlavorQCDDown").metLV = metLV_jecUnc.at("jecUncFlavorQCDDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncFlavorQCDDown").jet1LV = jet1LV_jecUnc.at("jecUncFlavorQCDDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncFlavorQCDDown").jet2LV = jet2LV_jecUnc.at("jecUncFlavorQCDDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEC2YearUp").metLV = metLV_jecUnc.at("jecUncEC2YearUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEC2YearUp").jet1LV = jet1LV_jecUnc.at("jecUncEC2YearUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEC2YearUp").jet2LV = jet2LV_jecUnc.at("jecUncEC2YearUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncEC2YearDown").metLV = metLV_jecUnc.at("jecUncEC2YearDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncEC2YearDown").jet1LV = jet1LV_jecUnc.at("jecUncEC2YearDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncEC2YearDown").jet2LV = jet2LV_jecUnc.at("jecUncEC2YearDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncHFYearUp").metLV = metLV_jecUnc.at("jecUncHFYearUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncHFYearUp").jet1LV = jet1LV_jecUnc.at("jecUncHFYearUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncHFYearUp").jet2LV = jet2LV_jecUnc.at("jecUncHFYearUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncHFYearDown").metLV = metLV_jecUnc.at("jecUncHFYearDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncHFYearDown").jet1LV = jet1LV_jecUnc.at("jecUncHFYearDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncHFYearDown").jet2LV = jet2LV_jecUnc.at("jecUncHFYearDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncAbsoluteYearUp").metLV = metLV_jecUnc.at("jecUncAbsoluteYearUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncAbsoluteYearUp").jet1LV = jet1LV_jecUnc.at("jecUncAbsoluteYearUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncAbsoluteYearUp").jet2LV = jet2LV_jecUnc.at("jecUncAbsoluteYearUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncAbsoluteYearDown").metLV = metLV_jecUnc.at("jecUncAbsoluteYearDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncAbsoluteYearDown").jet1LV = jet1LV_jecUnc.at("jecUncAbsoluteYearDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncAbsoluteYearDown").jet2LV = jet2LV_jecUnc.at("jecUncAbsoluteYearDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncBBEC1YearUp").metLV = metLV_jecUnc.at("jecUncBBEC1YearUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncBBEC1YearUp").jet1LV = jet1LV_jecUnc.at("jecUncBBEC1YearUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncBBEC1YearUp").jet2LV = jet2LV_jecUnc.at("jecUncBBEC1YearUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncBBEC1YearDown").metLV = metLV_jecUnc.at("jecUncBBEC1YearDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncBBEC1YearDown").jet1LV = jet1LV_jecUnc.at("jecUncBBEC1YearDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncBBEC1YearDown").jet2LV = jet2LV_jecUnc.at("jecUncBBEC1YearDown");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearUp").metLV = metLV_jecUnc.at("jecUncRelativeSampleYearUp");
         if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearUp").jet1LV = jet1LV_jecUnc.at("jecUncRelativeSampleYearUp");
         if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearUp").jet2LV = jet2LV_jecUnc.at("jecUncRelativeSampleYearUp");
         if (!isSampleForRecoilCorrection && jet1.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearDown").metLV = metLV_jecUnc.at("jecUncRelativeSampleYearDown");
         if(jet1.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearDown").jet1LV = jet1LV_jecUnc.at("jecUncRelativeSampleYearDown");
         if(jet2.E() != 0) uncertainty_map.at("jecUncRelativeSampleYearDown").jet2LV = jet2LV_jecUnc.at("jecUncRelativeSampleYearDown");
         
         for(auto &uncert : uncertainty_map){
            bool is_data_or_embedded = isData || (isEmbedded && !uncert.first.Contains("escale") && !uncert.first.Contains("ereso"));
            
            propagate_uncertainty( uncert.first,
                                   uncert.second.metLV, covMET, inputFile_visPtResolution,
                                   uncert.second.muonLV,
                                   uncert.second.electronLV,
                                   uncert.second.jet1LV,
                                   uncert.second.jet2LV,
                                   uncert.second.container,
                                   is_data_or_embedded, checkSV, checkFastMTT);
            
         }

         //set MELA variables
         if (njets>1){

            TLorentzVector tau1, tau2;
            tau1.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_1);
            tau2.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_2);

            // FIXME: TODO: Why do we not use the jet mass here? (comment from KIT)
            TLorentzVector jet1, jet2;
            jet1.SetPtEtaPhiM(jpt_1, jeta_1, jphi_1, 0);
            jet2.SetPtEtaPhiM(jpt_2, jeta_2, jphi_2, 0);

            // Run MELA
            SimpleParticleCollection_t daughters;
            if (q_1 * q_2 < 0){
               daughters.push_back(SimpleParticle_t(15 * q_1, tau1));
               daughters.push_back(SimpleParticle_t(15 * q_2, tau2));
            }
            else { //Sanitize charge for application on same-sign events
               daughters.push_back(SimpleParticle_t(15 * q_1, tau1)); 
               daughters.push_back(SimpleParticle_t(-15 * q_1, tau2));
            }
            SimpleParticleCollection_t associated;
            associated.push_back(SimpleParticle_t(0, jet1));
            associated.push_back(SimpleParticle_t(0, jet2));

            SimpleParticleCollection_t associated2;
            associated2.push_back(SimpleParticle_t(0, jet2));
            associated2.push_back(SimpleParticle_t(0, jet1));

            mela.resetInputEvent();
            mela.setCandidateDecayMode(TVar::CandidateDecay_ff); //decay into fermions
            mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t *)0, false); //set the decay products and the associated jets

            // Hypothesis: SM VBF Higgs
            mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
            mela.computeProdP(ME_vbf, false);
            mela.computeVBFAngles(ME_q2v1, ME_q2v2, ME_costheta1, ME_costheta2, ME_phi, ME_costhetastar, ME_phi1);

            // Hypothesis ggH + 2 jets
            mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
            mela.selfDHggcoupl[0][gHIGGS_GG_2][0] = 1;
            mela.computeProdP(ME_ggh, false);
            
            // Hypothesis: Z + 2 jets
            // Compute the Hypothesis with flipped jets and sum them up for the discriminator.
            mela.setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
            mela.computeProdP(ME_z2j_1, false);

            mela.resetInputEvent();
            mela.setInputEvent(&daughters, &associated2, (SimpleParticleCollection_t *)0, false);
            mela.computeProdP(ME_z2j_2, false);

            // Compute discriminator for VBF vs Z
            if ((ME_vbf + ME_z2j_1 + ME_z2j_2) != 0.0)
               {
                  ME_vbf_vs_Z = ME_vbf / (ME_vbf + ME_z2j_1 + ME_z2j_2);
               }
            else
               {
                  std::cout << "WARNING: ME_vbf_vs_Z = X / 0. Setting it to default " << -10 << std::endl;
                  ME_vbf_vs_Z = -10;
               }

            // Compute discriminator for ggH vs Z
            if ((ME_ggh + ME_z2j_1 + ME_z2j_2) != 0.0)
               {
                  ME_ggh_vs_Z = ME_ggh / (ME_ggh + ME_z2j_1 + ME_z2j_2);
               }
            else
               {
                  std::cout << "WARNING: ME_ggh_vs_Z = X / 0. Setting it to default " << -10 << std::endl;
                  ME_ggh_vs_Z = -10;
               }

            // Compute discriminator for VBF vs ggH
            if ((ME_vbf + ME_ggh) != 0.0)
               {
                  ME_vbf_vs_ggh = ME_vbf / (ME_vbf + ME_ggh);
               }
            else
               {
                  std::cout << "WARNING: ME_vbf_vs_ggh = X / 0. Setting it to default " << -10 << std::endl;
                  ME_vbf_vs_ggh = -10;
               }
         }

         tree->Fill();
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

