import FWCore.ParameterSet.Config as cms

isData = True
is25ns = True
year = 2016
period = 'Spring16'

#if "@SKIM@".lower()=="0":
#    skim = 0
#configurable options =======================================================================
runOnData=isData #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file) /// OUTDATED for 25ns
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false  == existing as slimmedMETsNoHF
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
#===================================================================


process = cms.Process("TreeProducer")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load('Configuration/StandardSequences/Services_cff')
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")


# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 500
# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if runOnData:
  process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'
else:
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'

print 'The conditions are =======>',process.GlobalTag.globaltag


if usePrivateSQlite:
    from CondCore.CondDB.CondDB_cfi import *
    CondDBSetup = CondDB.clone()
    CondDBSetup.__delattr__('connect')
    import os
    if runOnData:
      era="Spring16_25nsV3_DATA"
    else:
      era="Spring16_25nsV3_MC"
    
    dBFile = os.path.expandvars("$CMSSW_BASE/src/DesyTauAnalyses/NTupleMaker/data/JEC/"+era+".db")
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string( "sqlite_file://"+dBFile ),
                               toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
                ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
                ),
            )
                               )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

### =====================================================================================================


### ReRun JEC ===========================================================================================

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

### END ReRun JEC ======================================================================================

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
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



### END Electron ID ====================================================================================


## HBHE noise Filter

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")




fnames = []
if runOnData:
  fnames.append('/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/448/00000/CECFFCBE-CE1C-E611-8660-02163E011A4E.root')
else:
  fnames.append('/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00F0B3DC-211B-E611-A6A0-001E67248A39.root')
    
# Define the input source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring( fnames ),
                            skipEvents = cms.untracked.uint32(0) #33333
)



# filters ==============================================================================================

from RecoMET.METFilters.metFilters_cff import HBHENoiseFilterResultProducer, HBHENoiseFilter, HBHENoiseIsoFilter, hcalLaserEventFilter
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter, eeBadScFilter, ecalLaserCorrFilter, EcalDeadCellBoundaryEnergyFilter
from RecoMET.METFilters.metFilters_cff import primaryVertexFilter, CSCTightHaloFilter, CSCTightHaloTrkMuUnvetoFilter, CSCTightHalo2015Filter, globalTightHalo2016Filter, globalSuperTightHalo2016Filter, HcalStripHaloFilter
from RecoMET.METFilters.metFilters_cff import goodVertices, trackingFailureFilter, trkPOGFilters, manystripclus53X, toomanystripclus53X, logErrorTooManyClusters
from RecoMET.METFilters.metFilters_cff import chargedHadronTrackResolutionFilter, muonBadTrackFilter
from RecoMET.METFilters.metFilters_cff import metFilters

# individual filters
Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseFilter)
Flag_HBHENoiseIsoFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseIsoFilter)
Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(CSCTightHaloTrkMuUnvetoFilter)
Flag_CSCTightHalo2015Filter = cms.Path(CSCTightHalo2015Filter)
Flag_globalTightHalo2016Filter = cms.Path(globalTightHalo2016Filter)
Flag_globalSuperTightHalo2016Filter = cms.Path(globalSuperTightHalo2016Filter)
Flag_HcalStripHaloFilter = cms.Path(HcalStripHaloFilter)
Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(EcalDeadCellBoundaryEnergyFilter)
Flag_goodVertices = cms.Path(primaryVertexFilter)
Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
Flag_eeBadScFilter = cms.Path(eeBadScFilter)
Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
Flag_trkPOGFilters = cms.Path(trkPOGFilters)
Flag_chargedHadronTrackResolutionFilter = cms.Path(chargedHadronTrackResolutionFilter)
Flag_muonBadTrackFilter = cms.Path(muonBadTrackFilter)
# and the sub-filters
Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)

Flag_METFilters = cms.Path(metFilters)

allMetFilterPaths=['HBHENoiseFilter','HBHENoiseIsoFilter','CSCTightHaloFilter','CSCTightHaloTrkMuUnvetoFilter','CSCTightHalo2015Filter','globalTightHalo2016Filter','globalSuperTightHalo2016Filter','HcalStripHaloFilter','hcalLaserEventFilter','EcalDeadCellTriggerPrimitiveFilter','EcalDeadCellBoundaryEnergyFilter','goodVertices','eeBadScFilter',
'ecalLaserCorrFilter','trkPOGFilters','chargedHadronTrackResolutionFilter','muonBadTrackFilter','trkPOG_manystripclus53X','trkPOG_toomanystripclus53X','trkPOG_logErrorTooManyClusters','METFilters']


###########################################
######   NTUPLE PRODUCER    ###############
###########################################
process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(True),
SusyInfo = cms.untracked.bool(False),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(True),
RecPFMet = cms.untracked.bool(True),
RecMvaMet = cms.untracked.bool(False),                                      
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(True),
RecTau = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
JECfile = cms.untracked.string("DesyTauAnalyses/NTupleMaker/data/Fall15_25nsV2/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt"),
# collections
L1MuonCollectionTag = cms.InputTag("l1extraParticles"),
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
ElectronCollectionTag = cms.InputTag("slimmedElectrons"),
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
IsoTauCollectionTag = cms.InputTag("l1extraParticles:IsoTau:@skipCurrentProcess"),
JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs"),
MetCorrCollectionTag = cms.InputTag("slimmedMETs"),
MetCovMatrixTag = cms.InputTag("METSignificance:METCovariance:TreeProducer"),
MetSigTag = cms.InputTag("METSignificance:METSignificance:TreeProducer"),
MetCorrCovMatrixTag = cms.InputTag("METCorrSignificance:METCovariance:TreeProducer"),
MetCorrSigTag = cms.InputTag("METCorrSignificance:METSignificance:TreeProducer"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi"),
MvaMetCollectionsTag = cms.VInputTag("patMvaMetDiTau", "patMvaMetTauMu", "patMvaMetTauEle", "patMvaMetMuEle"),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
SusyMotherMassTag = cms.InputTag("susyInfo","SusyMotherMass"),
SusyLSPMassTag = cms.InputTag("susyInfo","SusyLSPMass"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_IsoMu18_v',
'HLT_IsoMu20_v',
'HLT_IsoMu22_v', #######no prescaled
'HLT_Mu27_v',
'HLT_Mu45_eta2p1_v',
'HLT_Mu17_Mu8_v',
'HLT_Mu17_Mu8_DZ_v',
'HLT_Mu17_Mu8_SameSign_DZ_v',
#'HLT_Mu20_Mu10_v',
#'HLT_Mu20_Mu10_DZ_v',
#'HLT_Mu20_Mu10_SameSign_DZ_v',
'HLT_DoubleMu23NoFiltersNoVtxDisplaced_v',
'HLT_DoubleMu28NoFiltersNoVtxDisplaced_v',
'HLT_DoubleMu33NoFiltersNoVtx_v',
'HLT_DoubleMu38NoFiltersNoVtx_v',
'HLT_TripleMu_5_3_3_v',
'HLT_TripleMu_12_10_5_v',
'HLT_Ele23_WPLoose_Gsf_v',
'HLT_Ele25_eta2p1_WPTight_Gsf_v',
'HLT_Ele27_WPTight_Gsf_v', ##### no prescaled for 1e34

'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v', ##########filters
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', ########## filters no prescales

'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v',
'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v',
'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v',
'HLT_PFJet60_v',
'HLT_PFJet80_v',
'HLT_PFJet140_v'
),
TriggerProcess = cms.untracked.string("HLT"),
Flags = cms.untracked.vstring(
#  'Flag_HBHENoiseFilter',
#  'Flag_HBHENoiseIsoFilter',
#  'Flag_CSCTightHalo2015Filter',
#  'Flag_EcalDeadCellTriggerPrimitiveFilter',
#  'Flag_goodVertices',
#  'Flag_eeBadScFilter',
#  'Flag_chargedHadronTrackResolutionFilter',
#  'Flag_muonBadTrackFilter'
   'allMetFilterPaths'
),
FlagsProcess = cms.untracked.string("RECO"),
# tracks
RecTrackPtMin = cms.untracked.double(1.0),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackDxyMax = cms.untracked.double(1.0),
RecTrackDzMax = cms.untracked.double(1.0),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(5.),
RecMuonEtaMax = cms.untracked.double(2.5),
RecMuonHLTriggerMatching = cms.untracked.vstring(
# isolated muons
'HLT_IsoMu18_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09',
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09',
'HLT_IsoMu22_v.*:hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09',
# muons
'HLT_Mu27_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q',
'HLT_Mu45_eta2p1_v.*:hltL3fL1sMu22Or25L1f0L2f10QL3Filtered45e2p1Q',
# Mu17_Mu8
'HLT_Mu17_Mu8_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_DZ_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2SameSign',
# Mu20_Mu10
'HLT_Mu20_Mu10_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered20,hltL3fL1sDoubleMu114ORDoubleMu125L1f0L2f10OneMuL3Filtered20',
'HLT_Mu20_Mu10_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered10,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered10',
'HLT_Mu20_Mu10_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered20,hltL3fL1sDoubleMu114ORDoubleMu125L1f0L2f10OneMuL3Filtered20',
'HLT_Mu20_Mu10_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered10,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered10',
'HLT_Mu20_Mu10_DZ_v.*:hltDiMuonGlb20Glb10DzFiltered0p2',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered20,hltL3fL1sDoubleMu114ORDoubleMu125L1f0L2f10OneMuL3Filtered20',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered10,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered10',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltDiMuonGlb20Glb10DzFiltered0p2',
'HLT_Mu20_Mu10_SameSign_DZ_v.*:hltDiMuonGlb20Glb10DzFiltered0p2SameSign',
# DoubleMuons (large IP)
'HLT_DoubleMu23NoFiltersNoVtxDisplaced_v.*:hltL3fDimuonL1f0L2NVf10L3NoFiltersNoVtxDisplacedFiltered23',
'HLT_DoubleMu28NoFiltersNoVtxDisplaced_v.*:hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxDisplacedFiltered28',
'HLT_DoubleMu33NoFiltersNoVtx_v.*:hltL3fDimuonL1f0L2NVf10L3NoFiltersNoVtxFiltered33',
'HLT_DoubleMu38NoFiltersNoVtx_v.*:hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxFiltered38',
# TripleMuons
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5',
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105',
'HLT_TripleMu_12_10_5_v.*:hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105',
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered533',
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered3',
# MuonElectron
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17',
# DoubleMuon
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8,hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
#                                        hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8 
#                                        hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
#                                        hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2'
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
# SingleElectron
'HLT_Ele23_WPLoose_Gsf_v.*:hltEle23WPLooseGsfTrackIsoFilter',
'HLT_Ele25_eta2p1_WPTight_Gsf_v*:hltEle25erWPTightGsfTrackIsoFilter',
'HLT_Ele27_WPTight_Gsf_v.*:hltEle27WPTightGsfTrackIsoFilter',
# MuonElectron
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
# DoubleElectron
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(18),
RecTauEtaMax = cms.untracked.double(2.4),                                      
RecTauHLTriggerMatching = cms.untracked.vstring(
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
'HLT_PFJet140_v.*:hltSinglePFJet140'
),
RecJetBtagDiscriminators = cms.untracked.vstring(
'pfCombinedInclusiveSecondaryVertexV2BJetTags',
'pfJetProbabilityBJetTags'
),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("Data") 
)

#process.patJets.addBTagInfo = cms.bool(True)
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.p = cms.Path(
#                     process.particleFlowPtrs +
#                     process.makePatElectrons +
#                     process.makePatMuons + 
#                     process.makePatJets +
#	              process.selectPrimaryVertex *
#                     process.patDefaultSequence * 
#                     process.egmGsfElectronIDSequence *
#                     process.mvaTrigV025nsCSA14 * 
#                     process.mvaNonTrigV025nsCSA14 * 
    #process.ak4PFJets*
#    process.mvaMetSequence*
#    process.puJetIdSequence*
    process.METSignificance *
    process.egmGsfElectronIDSequence *	
    process.makeroottree
    )






process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root")
				   )


#process.end = cms.EndPath(process.Out*process.TFileService)

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



