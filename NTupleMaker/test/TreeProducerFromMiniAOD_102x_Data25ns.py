import FWCore.ParameterSet.Config as cms

isData = True
is25ns = True
year = 2018
period = '2018'

#configurable options =======================================================================
runOnData=isData #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file) /// OUTDATED for 25ns
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false  == existing as slimmedMETsNoHF
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
#===================================================================

# Define the CMSSW process
process = cms.Process("TreeProducer")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '102X_dataRun2_v10'

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
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
        '/store/data/Run2018A/EGamma/MINIAOD/17Sep2018-v2/60000/FFD234BD-747F-9242-9EC3-4D3BC8E564B0.root'
        ),
  skipEvents = cms.untracked.uint32(0)
)

### JECs ==============================================================================================
# From :  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  labelName = 'UpdatedJEC',
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

### END JECs ==========================================================================================

### MET ===============================================================================================
# from : https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#We_have_a_tool_to_help_you_to_ap

# PFMET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
                           isData=runOnData,
                           #fixEE2017 = False,
                           #fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
                           postfix = "ModifiedMET"
                           )

# PuppiMET
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True );

# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
                           isData=runOnData,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )

process.puppiNoLep.useExistingWeights = False
process.puppi.useExistingWeights = False

# MET filter for 2018 (from : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM)
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
baddetEcallist = cms.vuint32(
                             [872439604,872422825,872420274,872423218,
                              872423215,872416066,872435036,872439336,
                              872420273,872436907,872420147,872439731,
                              872436657,872420397,872439732,872439339,
                              872439603,872422436,872439861,872437051,
                              872437052,872420649,872422436,872421950,
                              872437185,872422564,872421566,872421695,
                              872421955,872421567,872437184,872421951,
                              872421694,872437056,872437057,872437313])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
                                                        "EcalBadCalibFilter",
                                                        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
                                                        ecalMinEt        = cms.double(50.),
                                                        baddetEcal    = baddetEcallist, 
                                                        taggingMode = cms.bool(True),
                                                        debug = cms.bool(False)
                                                        )

### END MET ===========================================================================================

### Electron scale and smearing =======================================================================
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       applyEnergyCorrections=False,
#                       applyVIDOnCorrectedEgamma=False,
#                       isMiniAOD=True,
#                       era='2017-Nov17ReReco')
### END Electron scale and smearing ====================================================================

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
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',#new! Fall17ID
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

### END Electron ID ====================================================================================

# Tau ID ===============================================================================================
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Running_of_the_DNN_based_tau_ID

updatedTauName = "NewTauIDsEmbedded" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
                                          updatedTauName = updatedTauName,
                                          toKeep = [ "2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1", "deepTau2017v1", "DPFTau_2016_v0", "DPFTau_2016_v1" ]
                                          )

tauIdEmbedder.runTauID()
# END Tau ID ===========================================================================================

# NTuple Maker =========================================================================================

process.initroottree = cms.EDAnalyzer("InitAnalyzer",
IsData = cms.untracked.bool(isData),
#IsData = cms.untracked.bool(False),
GenParticles = cms.untracked.bool(not isData),
GenJets = cms.untracked.bool(not isData)
)

process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
IsEmbedded = cms.untracked.bool(False),
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
RecTrack = cms.untracked.bool(True),
RecPFMet = cms.untracked.bool(True),
RecPFMetCorr = cms.untracked.bool(True),
RecPuppiMet = cms.untracked.bool(True),
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
applyElectronESShift = cms.untracked.bool(False),
#######new in 9.4.0
eleMvanoIsoWP90Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
eleMvanoIsoWP80Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
eleMvanoIsoWPLooseFall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
eleMvaIsoWP90Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"),
eleMvaIsoWP80Fall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
eleMvaIsoWPLooseFall17Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
mvaValuesIsoFall17Map = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
mvaValuesnoIsoFall17Map = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
eleVetoIdFall17Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
eleLooseIdFall17Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
eleMediumIdFall17Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
eleTightIdFall17Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),
eleVetoIdFall17V2Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
eleLooseIdFall17V2Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
eleMediumIdFall17V2Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
eleTightIdFall17V2Map = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
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
TauCollectionTag = cms.InputTag("NewTauIDsEmbedded"),
L1MuonCollectionTag = cms.InputTag("gmtStage2Digis:Muon"),
L1EGammaCollectionTag = cms.InputTag("caloStage2Digis:EGamma"),
L1TauCollectionTag = cms.InputTag("caloStage2Digis:Tau"),
L1JetCollectionTag = cms.InputTag("caloStage2Digis:Jet"),
#JetCollectionTag = cms.InputTag("slimmedJets"),
JetCollectionTag = cms.InputTag("updatedPatJetsUpdatedJEC::TreeProducer"),
MetCollectionTag = cms.InputTag("slimmedMETs::@skipCurrentProcess"),
MetCorrCollectionTag = cms.InputTag("slimmedMETsModifiedMET::TreeProducer"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi::TreeProducer"),
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
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v',
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
'HLT_Mu18_Mu9_v',
'HLT_Mu18_Mu9_DZ_v',
'HLT_Mu18_Mu9_SameSign_v',
'HLT_Mu18_Mu9_SameSign_DZ_v',
'HLT_Mu20_Mu10_v',
'HLT_Mu20_Mu10_DZ_v',
'HLT_Mu20_Mu10_SameSign_v',
'HLT_Mu20_Mu10_SameSign_DZ_v',
'HLT_Mu37_TkMu27_v',
# Triple muon
'HLT_TripleMu_12_10_5_v',
'HLT_TripleMu_10_5_5_DZ_v',
# Muon+Electron triggers
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
# Ditau triggers
'HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v',
'HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v',
'HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v',
# Single tau triggers
'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v',
'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v',
# MET Triggers
'HLT_PFMET110_PFMHT110_IDTight_v',
'HLT_PFMET120_PFMHT120_IDTight_v',
'HLT_PFMET130_PFMHT130_IDTight_v',
'HLT_PFMET140_PFMHT140_IDTight_v',
'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v',
'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v',
'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v',
# Single-Jet Triggers
'HLT_PFJet40_v',
'HLT_PFJet60_v',
'HLT_PFJet80_v',
'HLT_PFJet140_v',
'HLT_PFJet200_v',
'HLT_PFJet260_v',
'HLT_PFJet320_v',
'HLT_PFJet400_v',
'HLT_PFJet450_v',
'HLT_PFJet500_v',
'HLT_PFJet550_v',
# Di-Jet Triggers
'HLT_DiPFJetAve40_',
'HLT_DiPFJetAve60_',
'HLT_DiPFJetAve80_',
'HLT_DiPFJetAve140_',
'HLT_DiPFJetAve200_',
'HLT_DiPFJetAve260_',
'HLT_DiPFJetAve320_',
'HLT_DiPFJetAve400_',
'HLT_DiPFJetAve500_'
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
RecTrackPtMin = cms.untracked.double(1.0),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackDxyMax = cms.untracked.double(1.0),
RecTrackDzMax = cms.untracked.double(1.0),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(2.),
RecMuonEtaMax = cms.untracked.double(2.5),
RecMuonHLTriggerMatching = cms.untracked.vstring(
#SingleMuon
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07',
'HLT_IsoMu24_v.*:hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu24_eta2p1_v.*:hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu27_v.*:hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07',
'HLT_IsoMu30_v.*:hltL3crIsoL1sMu22Or25L1f0L2f10QL3f30QL3trkIsoFiltered0p07',
'HLT_Mu50_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL1sMu18erTau24erIorMu20erTau24er',
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07',
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltL1sSingleMu22er',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoMu24LooseChargedIsoPFTau20',
'HLT_DoubleIsoMu20_eta2p1_v.*:hltL3crIsoL1sDoubleMu18erL1f0L2f10QL3f20QL3trkIsoFiltered0p07',
'HLT_DoubleIsoMu24_eta2p1_v.*:hltL3crIsoL1sDoubleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8,hltL3fL1sMu7EG23f0Filtered8',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu18_Mu9_SameSign_DZ_v.*:hltL1sDoubleMu125to157',
'HLT_Mu18_Mu9_SameSign_DZ_v.*:hltL3fL1DoubleMu157fFiltered9',
'HLT_Mu18_Mu9_SameSign_DZ_v.*:hltL3fL1DoubleMu157fFiltered18',
'HLT_Mu18_Mu9_SameSign_DZ_v.*:hltDiMuon189SameSignFiltered',
'HLT_Mu18_Mu9_SameSign_DZ_v.*:hltDiMuon189SameSignDzFiltered0p2',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltL1sDoubleMu125to157',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltL3fL1DoubleMu157fFiltered10',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltL3fL1DoubleMu157fFiltered20',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltDiMuon2010SameSignFiltered',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltDiMuon2010SameSignDzFiltered0p2',
'HLT_Mu37_TkMu27_v.*:hltDiMuonGlbFiltered37TrkFiltered27',
'HLT_Mu37_TkMu27_v.*:hltPreMu37TkMu27',
'HLT_Mu37_TkMu27_v.*:hltL3fL1sMu16orMu25L1f0L2f25L3Filtered37',
'HLT_Mu37_TkMu27_v.*:hltDiMuonGlb37Trk27DzFiltered0p2',
'HLT_Mu37_TkMu27_v.*:hltL1sSingleMu22IorSingleMu25',
'HLT_TripleMu_12_10_5_v.*:hltL1sTripleMu0IorTripleMu553',
'HLT_TripleMu_12_10_5_v.*:hltL3fL1TripleMu553f0PreFiltered555',
'HLT_TripleMu_12_10_5_v.*:hltL3fL1TripleMu553f0Filtered10105',
'HLT_TripleMu_12_10_5_v.*:hltL3fL1TripleMu553f0Filtered12105',
'HLT_TripleMu_10_5_5_DZ_v.*:hltL1sTripleMu0IorTripleMu553',
'HLT_TripleMu_10_5_5_DZ_v.*:hltL3fL1TripleMu553f0Filtered1055',
'HLT_TripleMu_10_5_5_DZ_v.*:hltL3fL1TripleMu553f0PreFiltered555',
'HLT_TripleMu_10_5_5_DZ_v.*:hltPreTripleMu1055DZ',
'HLT_TripleMu_10_5_5_DZ_v.*:hltTripleMu555TripleDZ0p2'
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
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltEle24erWPTightGsfTrackIsoFilterForTau',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30',
'HLT_DoubleEle24_eta2p1_WPTight_Gsf_v.*:hltDoubleEle24erWPTightGsfTrackIsoFilterForTau',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(15),
RecTauEtaMax = cms.untracked.double(2.5),
RecTauHLTriggerMatching = cms.untracked.vstring(
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltL1sMu18erTau24erIorMu20erTau24er',
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltPFTau27TrackLooseChargedIsoAgainstMuon',
#'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v.*:hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseChargedIsoAgainstMuon',
'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoMu24LooseChargedIsoPFTau20',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3,hltL1sIsoEG22erIsoTau26erdEtaMin0p2',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltPFTau30TrackLooseChargedIso',
'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v.*:hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30',
'HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v.*:hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg',
'HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v.*:hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg',
'HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v.*:hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg',
'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v.*:hltSelectedPFTau180MediumChargedIsolationL1HLTMatched',
'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v.*:hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong'
),
RecTauFloatDiscriminators = cms.untracked.vstring(),
RecTauBinaryDiscriminators = cms.untracked.vstring(),
RecTauNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(18.),
RecJetEtaMax = cms.untracked.double(5.2),
RecJetHLTriggerMatching = cms.untracked.vstring(
'HLT_PFJet40_v.*:hltSinglePFJet40',
'HLT_PFJet60_v.*:hltSinglePFJet60',
'HLT_PFJet80_v.*:hltSinglePFJet80',
'HLT_PFJet140_v.*:hltSinglePFJet140',
'HLT_PFJet200_v.*:hltSinglePFJet200',
'HLT_PFJet260_v.*:hltSinglePFJet260',
'HLT_PFJet320_v.*:hltSinglePFJet320',
'HLT_PFJet400_v.*:hltSinglePFJet400',
'HLT_PFJet450_v.*:hltSinglePFJet450',
'HLT_PFJet500_v.*:hltSinglePFJet500',
'HLT_PFJet550_v.*:hltSinglePFJet550',
'HLT_DiPFJetAve40_v.*:hltDiPFJetAve40',
'HLT_DiPFJetAve60_v.*:hltDiPFJetAve60',
'HLT_DiPFJetAve80_v.*:hltDiPFJetAve80',
'HLT_DiPFJetAve140_v.*:hltDiPFJetAve140',
'HLT_DiPFJetAve200_v.*:hltDiPFJetAve200',
'HLT_DiPFJetAve260_v.*:hltDiPFJetAve260',
'HLT_DiPFJetAve320_v.*:hltDiPFJetAve320',
'HLT_DiPFJetAve400_v.*:hltDiPFJetAve400',
'HLT_DiPFJetAve500_v.*:hltDiPFJetAve500'
),
RecJetBtagDiscriminators = cms.untracked.vstring(
'pfCombinedInclusiveSecondaryVertexV2BJetTags',
'pfDeepCSVJetTags:probb',
'pfDeepCSVJetTags:probbb'
),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("Data") 
)
#process.patJets.addBTagInfo = cms.bool(True)

process.p = cms.Path(
  process.initroottree *
  process.jecSequence *  # New JECs
  process.egmPhotonIDSequence * # Puppi MET
  process.puppiMETSequence *  # Puppi MET
  process.fullPatMetSequencePuppi *  # Re-correcting Puppi MET
  process.fullPatMetSequenceModifiedMET *  # Re-correcting PFMET
  process.ecalBadCalibReducedMINIAODFilter *  # MET filter 2018
  process.egmGsfElectronIDSequence *
  process.rerunMvaIsolationSequence *  # Tau IDs
  getattr(process,updatedTauName) *  # Tau IDs
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

