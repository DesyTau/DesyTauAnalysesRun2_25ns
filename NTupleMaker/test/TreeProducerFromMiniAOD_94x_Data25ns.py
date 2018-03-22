import FWCore.ParameterSet.Config as cms

isData = True
is25ns = True
year = 2017
period = '2017'

#configurable options =======================================================================
runOnData=isData #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file) /// OUTDATED for 25ns
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false  == existing as slimmedMETsNoHF
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
#===================================================================

### External JECs =====================================================================================================



# Define the CMSSW process
process = cms.Process("TreeProducer")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(100)
)

# Define the input source
process.source = cms.Source("PoolSource", 
  fileNames = cms.untracked.vstring(
        '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/0C144772-16E5-E711-B272-001E673971CA.root',
        '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/10A7ED13-90EB-E711-8BE5-484D7E8DF114.root',
                                    '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/2E4A57CA-71E7-E711-A7A1-A4BF0112BC3E.root',
                                    '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/3E5F02AC-33E7-E711-AE42-A0369FC5FBA4.root',
                                    '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/3E7110F5-07E6-E711-A856-44A842CF057E.root',
                                    '/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/424CDDA0-E7E5-E711-B839-44A842CFD5F2.root'
        ),
  skipEvents = cms.untracked.uint32(0)
)


# Electron ID ==========================================================================================

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',#new! Fall17ID
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

### END Electron ID ====================================================================================


# NTuple Maker =======================================================================

process.initroottree = cms.EDAnalyzer("InitAnalyzer",
IsData = cms.untracked.bool(isData),
#IsData = cms.untracked.bool(False),
GenParticles = cms.untracked.bool(not isData),
GenJets = cms.untracked.bool(not isData)
)

process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(not isData),
GenJets = cms.untracked.bool(not isData),
SusyInfo = cms.untracked.bool(False),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),
RecPFMetCorr = cms.untracked.bool(False),
RecPuppiMet = cms.untracked.bool(False),
RecMvaMet = cms.untracked.bool(False),                                      
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(True),
RecTau = cms.untracked.bool(True),
L1Objects = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
ElectronCollectionTag = cms.InputTag("slimmedElectrons"),
#######new in 9.4.0
eleMvanoIsoWP90Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
eleMvanoIsoWP80Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
eleMvanoIsoWPLooseFall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
eleMvaIsoWP90Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"),
eleMvaIsoWP80Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
eleMvaIsoWPLooseFall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
mvaValuesIsoFall17Map = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
mvaValuesnoIsoFall17Map = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
#######new in 8.0.25
eleMvaWP90GeneralMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
eleMvaWP80GeneralMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
mvaValuesMapSpring16     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
mvaCategoriesMapSpring16 = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
eleVetoIdSummer16Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
eleLooseIdSummer16Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
eleMediumIdSummer16Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
eleTightIdSummer16Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
###############
eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
eleMvaNonTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
eleMvaNonTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
eleMvaTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
eleMvaTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
mvaNonTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
mvaNonTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
mvaTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
mvaTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
TauCollectionTag = cms.InputTag("slimmedTaus"),
L1MuonCollectionTag = cms.InputTag("gmtStage2Digis:Muon"),
L1EGammaCollectionTag = cms.InputTag("caloStage2Digis:EGamma"),
L1TauCollectionTag = cms.InputTag("caloStage2Digis:Tau"),
L1JetCollectionTag = cms.InputTag("caloStage2Digis:Jet"),
#JetCollectionTag = cms.InputTag("patJetsReapplyJEC::TreeProducer"),
JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs"),
MetCorrCollectionTag = cms.InputTag("slimmedMETs"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi"),
MvaMetCollectionsTag = cms.VInputTag(cms.InputTag("MVAMET","MVAMET","TreeProducer")),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
GenJetCollectionTag = cms.InputTag("slimmedGenJets"),
TriggerObjectCollectionTag = cms.InputTag("slimmedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
LHEEventProductTag = cms.InputTag("externalLHEProducer"),
SusyMotherMassTag = cms.InputTag("susyInfo","SusyMotherMass"),
SusyLSPMassTag = cms.InputTag("susyInfo","SusyLSPMass"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
#SingleMuon
'HLT_IsoMu20_v',
'HLT_IsoMu24_v',
'HLT_IsoMu24_eta2p1_v',
'HLT_IsoMu27_v',
'HLT_IsoMu30_v',
'HLT_Mu50_v',
# Muon-Tau triggers
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v',
# SingleElectron
'HLT_Ele27_WPTight_Gsf_v',
'HLT_Ele32_WPTight_Gsf_v',
'HLT_Ele35_WPTight_Gsf_v',
'HLT_Ele38_WPTight_Gsf_v',
# Electron-Tau triggers
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v',
# Dilepton triggers
'HLT_DoubleEle24_eta2p1_WPTight_Gsf_v',
'HLT_DoubleIsoMu20_eta2p1_v',
'HLT_DoubleIsoMu24_eta2p1_v',
# Muon+Electron triggers
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
# MET Triggers
'HLT_PFMET110_PFMHT110_IDTight_v',
'HLT_PFMET120_PFMHT120_IDTight_v',
'HLT_PFMET130_PFMHT130_IDTight_v',
'HLT_PFMET140_PFMHT140_IDTight_v',
'HLT_PFMET110NoMu_PFMHT110NoMu_IDTight_v',
'HLT_PFMET120NoMu_PFMHT120NoMu_IDTight_v',
# Jet Triggers
'HLT_PFJet60_v',
'HLT_PFJet80_v',
'HLT_PFJet140_v',
'HLT_PFJet200_v',
'HLT_PFJet260_v',
'HLT_PFJet320_v',
'HLT_PFJet400_v',
'HLT_PFJet450_v',
'HLT_PFJet500_v',
'HLT_PFJet550_v'
),
TriggerProcess = cms.untracked.string("HLT"),
Flags = cms.untracked.vstring(
  'Flag_HBHENoiseFilter',
  'Flag_HBHENoiseIsoFilter',
  'Flag_CSCTightHalo2015Filter',
  'Flag_EcalDeadCellTriggerPrimitiveFilter',
  'Flag_goodVertices',
  'Flag_eeBadScFilter',
  'Flag_chargedHadronTrackResolutionFilter',
  'Flag_muonBadTrackFilter',
  'Flag_globalTightHalo2016Filter',
  'Flag_METFilters',
  'allMetFilterPaths'
),
FlagsProcesses = cms.untracked.vstring("RECO","PAT"),
BadChargedCandidateFilter =  cms.InputTag("BadChargedCandidateFilter"),
BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
BadGlobalMuons    = cms.InputTag("badGlobalMuonTagger","bad","TreeProducer"),
BadDuplicateMuons = cms.InputTag("cloneGlobalMuonTagger","bad","TreeProducer"),
# tracks
RecTrackPtMin = cms.untracked.double(0.5),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackDxyMax = cms.untracked.double(1.0),
RecTrackDzMax = cms.untracked.double(1.0),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(5.),
RecMuonEtaMax = cms.untracked.double(2.5),
RecMuonHLTriggerMatching = cms.untracked.vstring(
#SingleMuon
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07',
'HLT_IsoMu24_v.*:hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu24_eta2p1_v.*:hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu27_v.*:hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07',
'HLT_IsoMu30_v.*:hltL3crIsoL1sMu22Or25L1f0L2f10QL3f30QL3trkIsoFiltered0p07',
'HLT_Mu50_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL1sMu18erTau24erIorMu20erTau24er',
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07',
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltL1sSingleMu22er',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoMu24LooseChargedIsoPFTau20',
'HLT_DoubleIsoMu20_eta2p1_v.*:hltL3crIsoL1sDoubleMu18erL1f0L2f10QL3f20QL3trkIsoFiltered0p07',
'HLT_DoubleIsoMu24_eta2p1_v.*:hltL3crIsoL1sDoubleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter'
),
RecMuonNum = cms.untracked.int32(0),
# photons
RecPhotonPtMin = cms.untracked.double(20.),
RecPhotonEtaMax = cms.untracked.double(10000.),
RecPhotonHLTriggerMatching = cms.untracked.vstring(),
RecPhotonNum = cms.untracked.int32(0),
# electrons
RecElectronPtMin = cms.untracked.double(8.),
RecElectronEtaMax = cms.untracked.double(2.6),
RecElectronHLTriggerMatching = cms.untracked.vstring(
#SingleElectron
'HLT_Ele27_WPTight_Gsf_v.*:hltEle27WPTightGsfTrackIsoFilter',
'HLT_Ele32_WPTight_Gsf_v.*:hltEle32WPTightGsfTrackIsoFilter',
'HLT_Ele35_WPTight_Gsf_v.*:hltEle35noerWPTightGsfTrackIsoFilter',
'HLT_Ele38_WPTight_Gsf_v.*:hltEle38noerWPTightGsfTrackIsoFilter',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3,hltL1sIsoEG22erIsoTau26erdEtaMin0p2',
#'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltEle24erWPTightGsfTrackIsoFilterForTau',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30',
'HLT_DoubleEle24_eta2p1_WPTight_Gsf_v.*:hltDoubleEle24erWPTightGsfTrackIsoFilterForTau',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(15),
RecTauEtaMax = cms.untracked.double(2.5),
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL1sMu18erTau24erIorMu20erTau24er',
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltPFTau27TrackLooseChargedIsoAgainstMuon',
'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseChargedIsoAgainstMuon',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoMu24LooseChargedIsoPFTau20',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3,hltL1sIsoEG22erIsoTau26erdEtaMin0p2',
#'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltPFTau30TrackLooseChargedIso',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30'
),
RecTauFloatDiscriminators = cms.untracked.vstring(),
RecTauBinaryDiscriminators = cms.untracked.vstring(),
RecTauNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(18.),
RecJetEtaMax = cms.untracked.double(5.2),
RecJetHLTriggerMatching = cms.untracked.vstring(
'HLT_PFJet60_v.*:hltSinglePFJet60',
'HLT_PFJet80_v.*:hltSinglePFJet80',
'HLT_PFJet140_v.*:hltSinglePFJet140',
'HLT_PFJet200_v.*:hltSinglePFJet200',
'HLT_PFJet260_v.*:hltSinglePFJet260',
'HLT_PFJet320_v.*:hltSinglePFJet320',
'HLT_PFJet400_v.*:hltSinglePFJet400',
'HLT_PFJet450_v.*:hltSinglePFJet450',
'HLT_PFJet500_v.*:hltSinglePFJet500',
'HLT_PFJet550_v.*:hltSinglePFJet550'
),
RecJetBtagDiscriminators = cms.untracked.vstring(
'pfCombinedInclusiveSecondaryVertexV2BJetTags',
'pfJetProbabilityBJetTags'
),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("Data") 
)
#process.patJets.addBTagInfo = cms.bool(True)

process.p = cms.Path(
  process.initroottree*
#  process.BadChargedCandidateFilter *
#  process.BadPFMuonFilter *
#  process.BadGlobalMuonFilter *
#  process.pileupJetIdUpdated * 
#  process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
#  process.fullPatMetSequence * 
#  process.egmPhotonIDSequence *
#  process.puppiMETSequence *
#  process.fullPatMetSequencePuppi *
  process.egmGsfElectronIDSequence * 
#  process.rerunMvaIsolation2SeqRun2 * 
  #process.mvaMetSequence *
  #process.HBHENoiseFilterResultProducer* #produces HBHE bools baseline
  #process.ApplyBaselineHBHENoiseFilter*  #reject events based 
  #process.ApplyBaselineHBHEISONoiseFilter*  #reject events based -- disable the module, performance is being investigated fu
  process.makeroottree
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output_DATA.root")
                                 )

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string('output_particles_DATA.root'),
                                  outputCommands = cms.untracked.vstring(
                                    'keep *_*_bad_TreeProducer'#,
                                    #'drop patJets*_*_*_*'
                                    #'keep *_slimmedMuons_*_*',
                                    #'drop *_selectedPatJetsForMetT1T2Corr_*_*',
                                    #'drop patJets_*_*_*'
                                  ),        
                                  SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p'))
)

#process.end = cms.EndPath(process.output)

#processDumpFile = open('MyRootMaker.dump', 'w')
#print >> processDumpFile, process.dumpPython()

def miniAOD_customizeMETFiltersFastSim(process):
    """Replace some MET filters that don't work in FastSim with trivial bools"""
    for X in 'CSCTightHaloFilter', 'CSCTightHaloTrkMuUnvetoFilter','CSCTightHalo2015Filter','globalTightHalo2016Filter','globalSuperTightHalo2016Filter','HcalStripHaloFilter':
        process.globalReplace(X, cms.EDFilter("HLTBool", result=cms.bool(True)))    
    for X in 'manystripclus53X', 'toomanystripclus53X', 'logErrorTooManyClusters':
        process.globalReplace(X, cms.EDFilter("HLTBool", result=cms.bool(False)))
    return process



def customise_for_gc(process):
	import FWCore.ParameterSet.Config as cms
	from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper

	try:
		maxevents = __MAX_EVENTS__
		process.maxEvents = cms.untracked.PSet(
			input = cms.untracked.int32(max(-1, maxevents))
		)
	except:
		pass

	# Dataset related setup
	try:
		primaryFiles = [__FILE_NAMES__]
		process.source = cms.Source('PoolSource',
			skipEvents = cms.untracked.uint32(__SKIP_EVENTS__),
			fileNames = cms.untracked.vstring(primaryFiles)
		)
		try:
			secondaryFiles = [__FILE_NAMES2__]
			process.source.secondaryFileNames = cms.untracked.vstring(secondaryFiles)
		except:
			pass
		try:
			lumirange = [__LUMI_RANGE__]
			if len(lumirange) > 0:
				process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(lumirange)
				process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
		except:
			pass
	except:
		pass

	if hasattr(process, 'RandomNumberGeneratorService'):
		randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
		randSvc.populate()

	process.AdaptorConfig = cms.Service('AdaptorConfig',
		enable = cms.untracked.bool(True),
		stats = cms.untracked.bool(True),
	)

	# Generator related setup
	try:
		if hasattr(process, 'generator') and process.source.type_() != 'PoolSource':
			process.source.firstLuminosityBlock = cms.untracked.uint32(1 + __MY_JOBID__)
			print 'Generator random seed:', process.RandomNumberGeneratorService.generator.initialSeed
	except:
		pass

	return (process)

process = customise_for_gc(process)

# grid-control: https://ekptrac.physik.uni-karlsruhe.de/trac/grid-control
