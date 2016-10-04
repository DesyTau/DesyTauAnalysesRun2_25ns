import FWCore.ParameterSet.Config as cms

isData = False
isFastSim = False
is25ns = True
year = 2016
period = 'Spring16'

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

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.cerr.default = cms.untracked.PSet (
  limit = cms.untracked.int32(10),
  timespan = cms.untracked.int32(60)
)    


# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(1000)
)

### External JECs =====================================================================================================

#from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'
else:
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'

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


# Define the input source

fnames = []
if runOnData:
  fnames.append('/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/448/00000/CECFFCBE-CE1C-E611-8660-02163E011A4E.root')
else:
  fnames.append('/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/40000/0002BEE4-D55B-E611-B35D-0017A4770C08.root')

# Define the input source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring( fnames ),
                            skipEvents = cms.untracked.uint32(0) #33333
)

#####################################################
  
# Pairwise MVA MET ================================================================================= 

## PreSelection for pairwise MVA MEt
process.muonMVAMET = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 8. && abs(eta) < 2.5 && isPFMuon')
)

process.electronMVAMET = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 8. && abs(eta) < 2.5')
)

process.tauMVAMET = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && eta < 2.5 && eta > -2.5')
)

process.leptonPreSelectionSequence = cms.Sequence(process.muonMVAMET+
                                                  process.electronMVAMET+
                                                  process.tauMVAMET)

# mva MET
from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
runMVAMET( process, jetCollectionPF = "patJetsReapplyJEC"  )
process.MVAMET.srcLeptons  = cms.VInputTag("muonMVAMET", "electronMVAMET", "tauMVAMET")
process.MVAMET.requireOS = cms.bool(False)

process.mvaMetSequence  = cms.Sequence(process.leptonPreSelectionSequence +
                                       process.MVAMET)
# END Pairwise MVA MET ==============================================================

########### HBHE



process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)


# MET Filters ==============================================================================================

process.load('Configuration.StandardSequences.Services_cff')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

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


#ICHEPMetFilterPaths=['HBHENoiseFilter','HBHENoiseIsoFilter','globalTightHalo2016Filter','EcalDeadCellTriggerPrimitiveFilter','goodVertices','eeBadScFilter']

#Flag_allMETFilters = cms.Path(allMetFilterPaths)

#Flag_ICHEPMETFilters = cms.Path(ICHEPMetFilterPaths)

#####################################################################################

# NTuple Maker =======================================================================

process.initroottree = cms.EDAnalyzer("InitAnalyzer",
IsData = cms.untracked.bool(isData),
#IsData = cms.untracked.bool(False),
GenParticles = cms.untracked.bool(not isData)
)

JECfile = "DesyTauAnalyses/NTupleMaker/data/JEC/Spring16_25nsV6/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt"
if isData:
  JECfile = "DesyTauAnalyses/NTupleMaker/data/JEC/Spring16_25nsV6/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"

process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
JECfile = cms.untracked.string(JECfile),
# switches of collections
GenParticles = cms.untracked.bool(not isData),
SusyInfo = cms.untracked.bool(True),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),
RecPFMetCorr = cms.untracked.bool(True),
RecPuppiMet = cms.untracked.bool(True),
RecMvaMet = cms.untracked.bool(True),                                      
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(True),
RecTau = cms.untracked.bool(True),
L1Objects = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
ElectronCollectionTag = cms.InputTag("slimmedElectrons"),
#eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
#eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
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
JetCollectionTag = cms.InputTag("patJetsReapplyJEC::TreeProducer"),
#JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs::@skipCurrentProcess"),
MetCovMatrixTag = cms.InputTag("METSignificance:METCovariance:TreeProducer"),
MetSigTag = cms.InputTag("METSignificance:METSignificance:TreeProducer"),
MetCorrCovMatrixTag = cms.InputTag("METCorrSignificance:METCovariance:TreeProducer"),
MetCorrSigTag = cms.InputTag("METCorrSignificance:METSignificance:TreeProducer"),
MetCorrCollectionTag = cms.InputTag("patpfMETT1::TreeProducer"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi"),
MvaMetCollectionsTag = cms.VInputTag(cms.InputTag("MVAMET","MVAMET","TreeProducer")),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
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
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 
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
#JetHT
'HLT_PFJet60_v',
'HLT_PFJet80_v',
'HLT_PFJet140_v'
),
TriggerProcess = cms.untracked.string("HLT2"),
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
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu12EG10',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu20EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sSingleMu20erIorSingleMu22IorSingleMu25,hltL1sSingleMu20erlorSingleMu22lorSingleMu25',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG15',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG20IorMu5IsoEG18',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8', 
#MuonEG
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
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered533',
'HLT_TripleMu_5_3_3_v.*:hltL1TripleMu0L2TriMuFiltered0L3TriMuFiltered3',
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
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltEle32WPLooseGsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterIsoEle32WPLooseGsfLooseIsoPFTau20',
#MuonEG
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu12EG10',
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu20EG10',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sSingleMu20erIorSingleMu22IorSingleMu25,hltL1sSingleMu20erlorSingleMu22lorSingleMu25',
'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG15',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltL1sMu5EG20IorMu5IsoEG18',
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter', 
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
RecJetEtaMin = cms.untracked.double(4.9),
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

process.METCorrSignificance = process.METSignificance.clone(
  srcPfJets = cms.InputTag('patJetsReapplyJEC::TreeProducer'),
  srcMet = cms.InputTag('patpfMETT1::TreeProducer')
)

process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.debug = cms.bool(False)
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.BadPFMuonFilter.debug = cms.bool(False)

process.p = cms.Path(
  process.initroottree*
  process.BadChargedCandidateFilter *
  process.BadPFMuonFilter *
  process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
  process.egmGsfElectronIDSequence * 
  process.mvaMetSequence *
  process.METSignificance * process.METCorrSignificance *
  #process.HBHENoiseFilterResultProducer* #produces HBHE bools baseline
  #process.ApplyBaselineHBHENoiseFilter*  #reject events based 
  #process.ApplyBaselineHBHEISONoiseFilter*  #reject events based -- disable the module, performance is being investigated fu
  process.makeroottree
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output_MC.root")
                                 )

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string('output_particles_MC.root'),
                                  outputCommands = cms.untracked.vstring(
                                    'keep *_slimmedMETs_*_*',
				    'keep *_MVAMET_*_*',
                                    'keep *_patpfMETT1_*_*',
                                    'keep *_*MET*_*_*'
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
