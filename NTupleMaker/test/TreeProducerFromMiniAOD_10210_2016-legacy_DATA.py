import FWCore.ParameterSet.Config as cms

isData = True
is25ns = True
year = 2016
period = '2016'

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
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '94X_dataRun2_v10'

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True), 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(100)
)

# Define the input source
process.source = cms.Source("PoolSource", 
  fileNames = cms.untracked.vstring(
        "/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver1-v1/80000/F2494388-1D8C-E811-BB8E-0242AC1C0502.root"
        ),
  skipEvents = cms.untracked.uint32(0)
)

### JECs =====================================================================================================

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 
            'L2Relative', 
            'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

if runOnData:
    process.patJetCorrFactorsReapplyJEC.levels.append("L2L3Residual")

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJEC = updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
    )

#PFMET
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct and get the proper uncertainties
#runMetCorAndUncFromMiniAOD(process,
#                           isData=runOnData,
#                           fixEE2017 = True,
#                           fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
#                           postfix = "ModifiedMET"
#                           )

#PuppiMET
#from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
#makePuppiesFromMiniAOD( process, True );

#MET filter for 2018
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

# If you only want to re-cluster and get the proper uncertainties
#runMetCorAndUncFromMiniAOD(process,
#                           isData=runOnData,
#                           metType="Puppi",
#                           pfCandColl=cms.InputTag("puppiForMET"),
#                           recoMetFromPFCs=True,
#                           jetFlavor="AK4PFPuppi",
#                           postfix="Puppi"
#                           )

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
from DesyTauAnalyses.NTupleMaker.runTauIdMVA import *
na = TauIDEmbedder(process, cms, # pass tour process object
    debug=True,
#    toKeep = ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1", "deepTau2017v1", "DPFTau_2016_v0","DPFTau_2016_v1"]
    toKeep = ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
                   )
na.runTauID()

tauSrc = cms.InputTag('NewTauIDsEmbedded')
# END Tau ID ===========================================================================================

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
IsEmbedded = cms.untracked.bool(False),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(not isData),
GenJets = cms.untracked.bool(not isData),
SusyInfo = cms.untracked.bool(True),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(True),
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
JetCollectionTag = cms.InputTag("patJetsReapplyJEC::TreeProducer"),
MetCollectionTag = cms.InputTag("slimmedMETs::@skipCurrentProcess"),
#MetCorrCollectionTag = cms.InputTag("slimmedMETs::@skipCurrentProcess"),
#PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi::@skipCurrentProcess"),
#MetCollectionTag = cms.InputTag("slimmedMETs::TreeProducer"),
MetCorrCollectionTag = cms.InputTag("slimmedMETsModifiedMET::TreeProducer"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi::TreeProducer"),
MvaMetCollectionsTag = cms.VInputTag(cms.InputTag("MVAMET","MVAMET","TreeProducer")),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
GenJetCollectionTag = cms.InputTag("slimmedGenJets"),
TriggerObjectCollectionTag = cms.InputTag("slimmedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
PVwithBSCollectionTag =  cms.InputTag("MiniAODRefitVertexBSProducer"),
RefittedPVCollectionTag =  cms.InputTag("AdvancedRefitVertexNoBSProducer"),
RefittedwithBSPVCollectionTag =  cms.InputTag("AdvancedRefitVertexBSProducer"),
LHEEventProductTag = cms.InputTag("externalLHEProducer"),
SusyMotherMassTag = cms.InputTag("susyInfo","SusyMotherMass"),
SusyLSPMassTag = cms.InputTag("susyInfo","SusyLSPMass"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
#SingleMuon
'HLT_IsoMu18_v',
'HLT_IsoMu20_v',
'HLT_IsoMu22_v',
'HLT_IsoMu22_eta2p1_v',
'HLT_IsoMu24_v',
'HLT_IsoMu27_v',
'HLT_IsoTkMu18_v',
'HLT_IsoTkMu20_v',
'HLT_IsoTkMu22_v',
'HLT_IsoTkMu22_eta2p1_v',
'HLT_IsoTkMu24_v',
'HLT_IsoTkMu27_v',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v',
'HLT_Mu20_v',
'HLT_Mu24_eta2p1_v',
'HLT_Mu27_v',
'HLT_Mu45_eta2p1_v',
'HLT_Mu50_v',
'HLT_Mu17_TrkIsoVVL_v',
'HLT_Mu8_TrkIsoVVL_v', 
#SingleElectron
'HLT_Ele22_eta2p1_WPLoose_Gsf_v',
'HLT_Ele23_WPLoose_Gsf_v',
'HLT_Ele24_eta2p1_WPLoose_Gsf_v',
'HLT_Ele25_WPTight_Gsf_v',
'HLT_Ele25_eta2p1_WPLoose_Gsf_v',
'HLT_Ele25_eta2p1_WPTight_Gsf_v',
'HLT_Ele27_WPLoose_Gsf_v',
'HLT_Ele27_WPTight_Gsf_v',
'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
'HLT_Ele27_eta2p1_WPTight_Gsf_v',
'HLT_Ele32_eta2p1_WPTight_Gsf_v',
'HLT_Ele35_WPLoose_Gsf_v',
'HLT_Ele45_WPLoose_Gsf_v',  
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v',
'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v', # for RunG-H
'HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_v',
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
#DoubleMuon
'HLT_Mu17_Mu8_DZ_v',
'HLT_Mu17_Mu8_SameSign_DZ_v',
'HLT_Mu17_Mu8_SameSign_v',
'HLT_Mu17_Mu8_v',
'HLT_Mu17_TkMu8_DZ_v',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
'HLT_DoubleMu18NoFiltersNoVtx_v',
'HLT_DoubleMu23NoFiltersNoVtxDisplaced_v',
'HLT_DoubleMu28NoFiltersNoVtxDisplaced_v',
'HLT_DoubleMu33NoFiltersNoVtx_v',
'HLT_DoubleMu38NoFiltersNoVtx_v',
'HLT_TripleMu_12_10_5_v',
'HLT_TripleMu_5_3_3_v',
'HLT_DoubleIsoMu17_eta2p1_noDzCut_v',
'HLT_DoubleIsoMu17_eta2p1_v',
#DoubleElectron
'HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v',
'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v',
'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
'HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
#Tau
'HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v',
'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v',
'HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v',
#MET
'HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v',
'HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v',
'HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v',
'HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v',
'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v',
'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v',
'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v',
'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
#JetHT
'HLT_PFJet60_v',
'HLT_PFJet80_v',
'HLT_PFJet140_v',
'HLT_PFJet200_v',
'HLT_PFJet260_v',
'HLT_PFJet320_v',
'HLT_PFJet400_v',
'HLT_PFJet450_v',
'HLT_PFJet500_v',
'HLT_DiPFJetAve40_v',
'HLT_DiPFJetAve60_v',
'HLT_DiPFJetAve80_v',
'HLT_DiPFJetAve140_v',
'HLT_DiPFJetAve200_v',
'HLT_DiPFJetAve260_v',
'HLT_DiPFJetAve320_v',
'HLT_DiPFJetAve400_v',
'HLT_DiPFJetAve500_v',
'HLT_HT200_v',
'HLT_HT275_v',
'HLT_HT325_v',
'HLT_HT425_v',
'HLT_HT575_v'
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
'HLT_IsoMu18_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09',
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09',
'HLT_IsoMu22_v.*:hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09',
'HLT_IsoMu22_eta2p1_v.*:hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09',
'HLT_IsoMu24_v.*:hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09',
'HLT_IsoMu27_v.*:hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09',
'HLT_IsoTkMu18_v.*:hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09',
'HLT_IsoTkMu20_v.*:hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09',
'HLT_IsoTkMu22_v.*:hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09',
'HLT_IsoTkMu22_eta2p1_v.*:hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09',
'HLT_IsoTkMu24_v.*:hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09',
'HLT_IsoTkMu27_v.*:hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleMu16er',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu17LooseIsoPFTau20',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltL1sMu16erTau20er',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleMu18erIorSingleMu20er',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu19LooseIsoPFTau20',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v.*:hltL1sMu18erTau20er',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu19LooseIsoPFTau20',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleMu20erIorSingleMu22er',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu21LooseIsoPFTau20',
'HLT_Mu20_v.*:hltL3fL1sMu18L1f0L2f10QL3Filtered20Q',
'HLT_Mu24_eta2p1_v.*:hltL3fL1sMu20Eta2p1Or22L1f0L2f10QL3Filtered24Q',
'HLT_Mu27_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q',
'HLT_Mu45_eta2p1_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered45e2p1Q',
'HLT_Mu50_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
'HLT_Mu17_TrkIsoVVL_v.*:hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
'HLT_Mu17_TrkIsoVVL_v.*:hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
'HLT_Mu8_TrkIsoVVL_v.*:hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4',
'HLT_Mu8_TrkIsoVVL_v.*:hltL3fL1sMu5L1f0L2f5L3Filtered8',
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu12EG10',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu20EG10,hltL1sMu20EG10IorMu23EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sSingleMu20erIorSingleMu22IorSingleMu25,hltL1sSingleMu20erlorSingleMu22lorSingleMu25,hltL1sSingleMu20erIorSingleMu22IorSingleMu25IorMu20IsoEG6',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG15',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG20IorMu5IsoEG18,hltL1sMu5EG20IorMu5IsoEG18IorMu5IsoEG20IorMu5EG23',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8', 
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltL1sMu20EG10IorMu23EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltL1sMu5EG20IorMu5IsoEG18IorMu5IsoEG20IorMu5EG23',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
#DoubleMuon
'HLT_Mu17_Mu8_DZ_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2SameSign',
'HLT_Mu17_Mu8_SameSign_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_SameSign_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_SameSign_v.*:hltDiMuonGlb17Glb8SameSign',
'HLT_Mu17_Mu8_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_TkMu8_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
'HLT_Mu17_TkMu8_DZ_v.*:hltDiMuonGlbFiltered17TrkFiltered8',
'HLT_Mu17_TkMu8_DZ_v.*:hltDiMuonGlb17Trk8DzFiltered0p2',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.*:hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlbFiltered17TrkFiltered8',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v.*:hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v.*:hltDiMuonGlbFiltered17TrkFiltered8',
'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v.*:hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
'HLT_DoubleMu18NoFiltersNoVtx_v.*:hltL3fDimuonL1f0L2NVf10L3NoFiltersNoVtxFiltered18',
'HLT_DoubleMu23NoFiltersNoVtxDisplaced_v.*:hltL3fDimuonL1f0L2NVf10L3NoFiltersNoVtxDisplacedFiltered23',
'HLT_DoubleMu28NoFiltersNoVtxDisplaced_v.*:hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxDisplacedFiltered28',
'HLT_DoubleMu33NoFiltersNoVtx_v.*:hltL3fDimuonL1f0L2NVf10L3NoFiltersNoVtxFiltered33',
'HLT_DoubleMu38NoFiltersNoVtx_v.*:hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxFiltered38',  
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5',
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105',
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105',
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered533,hltL1TripleMu500L2TriMuFiltered0L3TriMuFiltered533',
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered3,hltL1TripleMu500L2TriMuFiltered0L3TriMuFiltered3',
'HLT_DoubleIsoMu17_eta2p1_noDzCut_v.*:hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3L3crIsoRhoFiltered0p15IterTrk02',
'HLT_DoubleIsoMu17_eta2p1_v.*:hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02'
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
'HLT_Ele22_eta2p1_WPLoose_Gsf_v.*:hltSingleEle22WPLooseGsfTrackIsoFilter',
'HLT_Ele23_WPLoose_Gsf_v.*:hltEle23WPLooseGsfTrackIsoFilter',
'HLT_Ele24_eta2p1_WPLoose_Gsf_v.*:hltSingleEle24WPLooseGsfTrackIsoFilter',
'HLT_Ele25_WPTight_Gsf_v.*:hltEle25WPTightGsfTrackIsoFilter',
'HLT_Ele25_eta2p1_WPLoose_Gsf_v.*:hltEle25erWPLooseGsfTrackIsoFilter',
'HLT_Ele25_eta2p1_WPTight_Gsf_v.*:hltEle25erWPTightGsfTrackIsoFilter',
'HLT_Ele27_WPLoose_Gsf_v.*:hltEle27noerWPLooseGsfTrackIsoFilter',
'HLT_Ele27_WPTight_Gsf_v.*:hltEle27WPTightGsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WPLoose_Gsf_v.*:hltEle27erWPLooseGsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WPTight_Gsf_v.*:hltEle27erWPTightGsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WPTight_Gsf_v.*:hltEle32WPTightGsfTrackIsoFilter',
'HLT_Ele35_WPLoose_Gsf_v.*:hltEle35WPLooseGsfTrackIsoFilter',
'HLT_Ele45_WPLoose_Gsf_v.*:hltEle45WPLooseGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleIsoEG20erIorSingleIsoEG22erIorSingleEG40',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleIsoEG22erIorSingleIsoEG24erIorSingleEG40',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20',  
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltL1sIsoEG22erTau20erdEtaMin0p2',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24er,hltL1sSingleEGor',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltEle27erWPLooseGsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoEle27WPLooseGsfLooseIsoPFTau20',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24er,hltL1sSingleEGor',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltEle32WPLooseGsfTrackIsoFilter,hltEle32erWPLooseGsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoEle32WPLooseGsfLooseIsoPFTau20',
'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter',
'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter',
'HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_v.*:hltEle45WPLooseGsfTrackIsoL1TauJetSeededFilter',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v.*:hltEle24WPLooseL1IsoEG22erIsoTau26erGsfTrackIsoFilter',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v.*:hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau30', # for RunG-H
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu12EG10',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu20EG10,hltL1sMu20EG10IorMu23EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sSingleMu20erIorSingleMu22IorSingleMu25,hltL1sSingleMu20erlorSingleMu22lorSingleMu25,hltL1sSingleMu20erIorSingleMu22IorSingleMu25IorMu20IsoEG6',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG15',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG20IorMu5IsoEG18,hltL1sMu5EG20IorMu5IsoEG18IorMu5IsoEG20IorMu5EG23',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter', 
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltL1sMu20EG10IorMu23EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltL1sMu5EG20IorMu5IsoEG18IorMu5IsoEG20IorMu5EG23',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter',
#DoubleEG
'HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v.*:hltEle24Ele22WPLooseGsfleg1TrackIsoFilter',
'HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v.*:hltEle24Ele22WPLooseGsfleg2TrackIsoFilter',
'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v.*:hltDiEle33CaloIdLGsfTrkIdVLMWPMS2UnseededFilter',
'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v.*:hltDiEle33CaloIdLPixelMatchUnseededFilter',
'HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v.*:hltEle37CaloIdLGsfTrkIdVLDPhiUnseededFilter',
'HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v.*:hltDiEle27CaloIdLGsfTrkIdVLDPhiUnseededFilter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(15),
RecTauEtaMax = cms.untracked.double(2.5),
RecTauHLTriggerMatching = cms.untracked.vstring(
#SingleMuon
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu17LooseIsoPFTau20',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu19LooseIsoPFTau20',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu19LooseIsoPFTau20',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoMu21LooseIsoPFTau20',
#SingleElectron
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoEle27WPLooseGsfLooseIsoPFTau20',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoEle32WPLooseGsfLooseIsoPFTau20',
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v.*:hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau30', # for RunG-H
'HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v.*:hltPFTau30TrackLooseIso',# for RunG-H
#Tau
'HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v.*:hltDoublePFTau32TrackPt1MediumIsolationDz02Reg',
'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v.*:hltDoublePFTau35TrackPt1MediumIsolationDz02Reg',
'HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v.*:hltDoublePFTau40TrackPt1MediumIsolationDz02Reg',
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
  process.initroottree*
#  process.BadChargedCandidateFilter *
#  process.BadPFMuonFilter *
#  process.BadGlobalMuonFilter *
#  process.pileupJetIdUpdated * 
#  process.rerunMvaIsolationSequence*

  process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
#  process.egmPhotonIDSequence *
#  process.puppiMETSequence *
#  process.fullPatMetSequencePuppi *
#  process.fullPatMetSequenceModifiedMET *
  process.ecalBadCalibReducedMINIAODFilter*
  process.egmGsfElectronIDSequence *
  process.rerunMvaIsolationSequence *      # add new tau ids
  process.NewTauIDsEmbedded *              # add new tau ids
  #process.rerunMvaIsolation2SeqRun2 *
  #process.mvaMetSequence *
  #process.HBHENoiseFilterResultProducer* #produces HBHE bools baseline
  #process.ApplyBaselineHBHENoiseFilter*  #reject events based 
  #process.ApplyBaselineHBHEISONoiseFilter*  #reject events based -- disable the module, performance is being investigated
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
