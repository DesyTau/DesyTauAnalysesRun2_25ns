import FWCore.ParameterSet.Config as cms

isData = True
is25ns = True
isPRv4 = True
isRepr05Oct = False
#skim = 0
year = 2015
period = 'Run2015D'

usePUJetID5X = True

#sampleName = 'MonteCarlo'
# sampleName = 'TTJets', "QCD", "DYJetsToLL_M50"

#if "@SKIM@".lower()=="0":
#    skim = 0

#sampleName = "@SAMPLE_NAME@"


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
# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(10000)
)


#configurable options =======================================================================
runOnData=isData #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file)
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false  == existing as slimmedMETsNoHF
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
#===================================================================


### External JECs =====================================================================================================

#from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = '74X_dataRun2_v2'
  if isPRv4:
    #process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v4'
    process.GlobalTag.globaltag = '74X_dataRun2_v5'
  if isRepr05Oct:
    #process.GlobalTag.globaltag = '74X_dataRun2_reMiniAOD_v0'
    process.GlobalTag.globaltag = '74X_dataRun2_v5'
else:
  if is25ns:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'
  else:
    process.GlobalTag.globaltag = '74X_mcRun2_startup_v2'

print 'The conditions are =======>',process.GlobalTag.globaltag
    
if usePrivateSQlite:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    if runOnData:
      era="Summer15_25nsV6_DATA"
    else:
      era="Summer15_25nsV6_MC"
    dBFile = os.path.expandvars("$CMSSW_BASE/src/DesyTauAnalyses/NTupleMaker/data/"+era+".db")
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

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 
            'L2Relative', 
            'L3Absolute',
            'L2L3Residual'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)


### END ReRun JEC ======================================================================================

### PFMET Corrections ==================================================================================

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
jecUncertaintyFile=""
if runOnData:
  if is25ns:
    jecUncertaintyFile="DesyTauAnalyses/NTupleMaker/data/Summer15_25nsV6/Summer15_25nsV6_DATA_UncertaintySources_AK4PFchs.txt"
    print ' JECUNCERTAINTY  ',jecUncertaintyFile
  else:
    jecUncertaintyFile="DesyTauAnalyses/NTupleMaker/data/Summer15_50nsV6/Summer15_50nsV6_DATA_UncertaintySources_AK4PFchs.txt"    
else:
  if is25ns:
    jecUncertaintyFile="DesyTauAnalyses/NTupleMaker/data/Summer15_25nsV6/Summer15_25nsV6_MC_UncertaintySources_AK4PFchs.txt"
    print ' JECUNCERTAINTY  ',jecUncertaintyFile
  else:
    jecUncertaintyFile="DesyTauAnalyses/NTupleMaker/data/Summer15_50nsV6/Summer15_50nsV6_MC_UncertaintySources_AK4PFchs.txt"

### ---------------------------------------------------------------------------
### Removing the HF from the MET computation
### ---------------------------------------------------------------------------
if not useHFCandidates:
    process.noHFCands = cms.EDFilter("CandPtrSelector",
                                     src=cms.InputTag("packedPFCandidates"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )

#jets are rebuilt from those candidates by the tools, no need to do anything else
### =================================================================================

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
runMetCorAndUncFromMiniAOD(process,
                           isData=runOnData,
                           )
process.slimmedMETs.t01Variation = cms.InputTag("slimmedMETs","","PAT") 
if (runOnData and isPRv4): 
  process.slimmedMETs.t01Variation = cms.InputTag("slimmedMETs","","RECO")
    
if not useHFCandidates:
    runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               pfCandColl=cms.InputTag("noHFCands"),
                               postfix="NoHF"
                               )
    process.slimmedMETsNoHF.t01Variation = cms.InputTag("slimmedMETsNoHF","","PAT") 
    if (runOnData and isPRv4): 
      process.slimmedMETsNoHF.t01Variation = cms.InputTag("slimmedMETsNoHF","","RECO")
          
### -------------------------------------------------------------------
### the lines below remove the L2L3 residual corrections when processing data
### -------------------------------------------------------------------
if not applyResiduals:
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")

    if not useHFCandidates:
          process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
          process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
          process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
          process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
          process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
          process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------

### END PFMET CORRECTIONS ==============================================================================


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

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

##################################################################################


# Define the input source
fnames = []
if runOnData:
  if is25ns:
      if isPRv4:
          fnames.append('/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/6CA1C627-246C-E511-8A6A-02163E014147.root')
      if isRepr05Oct:
          fnames.append('/store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/10000/021FD3F0-876F-E511-99D2-0025905A6060.root')
else:
  if is25ns:
    fnames.append('/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2_ext3-v1/10000/020B5100-426E-E511-888A-0026189437F9.root')
    #fnames.append('file:/nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_4/src/SUS-RunIISpring15FSPremix-00070.root')
    #fnames.append('file:/nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_4_14/src/miniAOD/Output/RootFiles/stau_stau400_LSP300_run55863_unwgt_GEN-SIM_chunk4.root')
  else:
    fnames.append('/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/AsymptNoPU_MCRUN2_74_V9A-v2/00000/02AD5DBB-1C0C-E511-8C41-00A0D1EE8E64.root')
    
# Define the input source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring( fnames ),
                            skipEvents = cms.untracked.uint32(0)
)

#####################################################
  

#--------------------------------------------------------------------------------
# produce PAT-tuple
####process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
#jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(
#    process,
#    jetSource = cms.InputTag('ak4PFJetsCHS'),
#    jetCorrections = ( 'AK4PFchs', jetCorrections, "" ),
#    outputModules = []
#)

#process.patJets.addTagInfos = cms.bool(True)
#process.patJets.tagInfoSources = cms.VInputTag("impactParameterTagInfosAOD","secondaryVertexTagInfosAOD","softMuonTagInfosAOD")
#process.patJets.tagInfoSources = cms.VInputTag("secondaryVertexTagInfosEI")


#--------------------------------------------------------------------------------
# switch to HPS PFTaus (and disable all "cleaning" cuts)
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process)

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
#process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

#--------------------------------------------------------------------------------
# electron Id (CSA14) ->
# input tags for electron id are hardcoded
# in class MyRootMaker    
#process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_CSA14_cfi")
#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff")

# specify correct sources for MVA electron Id.
#process.electronIDValueMapProducer.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
#process.electronIDValueMapProducer.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
#process.electronIDValueMapProducer.esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES")

#----------------------------------------------------------------------------------
# produce PU Jet Ids & MVA MET
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# produce PU Jet Ids & MVA MET
#----------------------------------------------------------------------------------


#####MVAMet

# Pairwise MVA MET ================================================================================= 
# as from jan 

## PreSelection for pairwise MVA MEt
## ## DiTau
process.tauPreSelectionDiTau = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 35. && abs(eta) < 2.5 && tauID("decayModeFindingNewDMs") > 0.5')
)
 
## ## TauEle
process.tauPreSelectionTauEle = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFindingNewDMs") > 0.5')
)
process.electronPreSelectionTauEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 18. && abs(eta) < 2.5')
)

## ## TauMu
process.tauPreSelectionTauMu = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFindingNewDMs") > 0.5')
)
process.muonPreSelectionTauMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 18. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)

## ## MuEle
process.muonPreSelectionMuEle = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 8. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)
process.electronPreSelectionMuEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 8. && abs(eta) < 2.5')
)
 
process.leptonPreSelectionSequence = cms.Sequence(process.tauPreSelectionDiTau+
                                                  process.tauPreSelectionTauEle+process.electronPreSelectionTauEle+
                                                  process.tauPreSelectionTauMu+process.muonPreSelectionTauMu+
                                                  process.muonPreSelectionMuEle+process.electronPreSelectionMuEle
                                                  )
# mva MET
from RecoMET.METPUSubtraction.mvaPFMET_cff import pfMVAMEt

mvaMETTauMu = cms.EDProducer('PFMETProducerMVATauTau', 
                        **pfMVAMEt.parameters_())#pfMVAMEt.clone()
if usePUJetID5X:
  mvaMETTauMu.srcCorrJets = cms.InputTag('calibratedAK4PFJetsForPFMVAMEt')
  mvaMETTauMu.srcUncorrJets = cms.InputTag('ak4PFJets')
else:
  mvaMETTauMu.srcCorrJets = cms.InputTag('calibratedAK4PFJetsCHSForPFMVAMEt')
  mvaMETTauMu.srcUncorrJets = cms.InputTag('ak4PFJetsCHS')

mvaMETTauMu.srcPFCandidates = cms.InputTag("packedPFCandidates")
mvaMETTauMu.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
mvaMETTauMu.srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauMu", "", ""),
                                               cms.InputTag("muonPreSelectionTauMu", "", ""))
mvaMETTauMu.permuteLeptons = cms.bool(True)

mvaMETTauMu.inputFileNames = cms.PSet(U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_25NS_July2015.root'),
                                      DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_25NS_July2015.root'),
                                      CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_25NS_July2015.root'),
                                      CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_25NS_July2015.root')
                                      )

process.mvaMETTauMu = mvaMETTauMu

process.mvaMETDiTau = cms.EDProducer('PFMETProducerMVATauTau',
                                      **mvaMETTauMu.parameters_())
process.mvaMETDiTau.srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionDiTau", "", ""),                                      
                                               cms.InputTag("tauPreSelectionDiTau", "", ""))


process.mvaMETMuEle = cms.EDProducer('PFMETProducerMVATauTau',
                                     **mvaMETTauMu.parameters_())
process.mvaMETMuEle.srcLeptons = cms.VInputTag(cms.InputTag("muonPreSelectionMuEle", "", ""),
                                               cms.InputTag("electronPreSelectionMuEle", "", ""))

process.mvaMETTauEle = cms.EDProducer('PFMETProducerMVATauTau',
                                      **mvaMETTauMu.parameters_())
process.mvaMETTauEle.srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauEle", "", ""),
                                                cms.InputTag("electronPreSelectionTauEle", "", ""))

process.mvaMETMuMu = cms.EDProducer('PFMETProducerMVATauTau',
                                     **mvaMETTauMu.parameters_())
process.mvaMETMuMu.srcLeptons = cms.VInputTag(cms.InputTag("muonPreSelectionMuEle", "", ""),
                                              cms.InputTag("muonPreSelectionMuEle", "", ""))

process.mvaMETEleEle = cms.EDProducer('PFMETProducerMVATauTau',
                                      **mvaMETTauMu.parameters_())
process.mvaMETEleEle.srcLeptons = cms.VInputTag(cms.InputTag("electronPreSelectionMuEle", "", ""),
                                                cms.InputTag("electronPreSelectionMuEle", "", ""))

if usePUJetID5X:
  process.ak4PFJets = cms.EDProducer("FastjetJetProducer",
                                     Active_Area_Repeats = cms.int32(1),
                                     doAreaFastjet = cms.bool(True),
                                     voronoiRfact = cms.double(-0.9),
                                     maxBadHcalCells = cms.uint32(9999999),
                                     doAreaDiskApprox = cms.bool(False),
                                     maxRecoveredEcalCells = cms.uint32(9999999),
                                     jetType = cms.string('PFJet'),
                                     minSeed = cms.uint32(14327),
                                     Ghost_EtaMax = cms.double(5.0),
                                     doRhoFastjet = cms.bool(False),
                                     jetAlgorithm = cms.string('AntiKt'),
                                     nSigmaPU = cms.double(1.0),
                                     GhostArea = cms.double(0.01),
                                     Rho_EtaMax = cms.double(4.4),
                                     maxBadEcalCells = cms.uint32(9999999),
                                     useDeterministicSeed = cms.bool(True),
                                     doPVCorrection = cms.bool(False),
                                     maxRecoveredHcalCells = cms.uint32(9999999),
                                     rParam = cms.double(0.4),
                                     maxProblematicHcalCells = cms.uint32(9999999),
                                     doOutputJets = cms.bool(True),
                                     src = cms.InputTag("packedPFCandidates"),
                                     inputEtMin = cms.double(0.0),
                                     srcPVs = cms.InputTag(""),
                                     jetPtMin = cms.double(3.0),
                                     radiusPU = cms.double(0.5),
                                     maxProblematicEcalCells = cms.uint32(9999999),
                                     doPUOffsetCorr = cms.bool(False),
                                     inputEMin = cms.double(0.0)
                                   )
  

  process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

  if runOnData:
    process.calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer("PFJetCorrectionProducer",
                                                            src = cms.InputTag("ak4PFJets"),
                                                            correctors = cms.vstring('ak4PFL1FastL2L3Residual')
                                                          )
  else:
    process.calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer("PFJetCorrectionProducer",
                                                            src = cms.InputTag("ak4PFJets"),
                                                            correctors = cms.vstring('ak4PFL1FastL2L3')
                                                          )  
  process.puJetIdForPFMVAMEt = cms.EDProducer("PileupJetIdProducer",
                                              algos = cms.VPSet(cms.PSet(
                                                tmvaVariables = cms.vstring('nvtx', 
                                                                            'jetPt', 
                                                                            'jetEta', 
                                                                            'jetPhi', 
                                                                            'dZ', 
                                                                            'beta', 
                                                                            'betaStar', 
                                                                            'nCharged', 
                                                                            'nNeutrals', 
                                                                            'dR2Mean', 
                                                                            'ptD', 
                                                                            'frac01', 
                                                                            'frac02', 
                                                                            'frac03', 
                                                                            'frac04', 
                                                                            'frac05'),
                                                etaBinnedWeights=cms.bool(False),
                                                tmvaMethod = cms.string('JetID'),
                                                cutBased = cms.bool(False),
                                                tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml.gz'),
                                                tmvaSpectators = cms.vstring(),
                                                label = cms.string('full'),
                                                version = cms.int32(-1),
                                                JetIdParams = cms.PSet(
                                                  Pt2030_Tight = cms.vdouble(0.3, 0.4, 0.7, 0.8),
                                                  Pt2030_Loose = cms.vdouble(0.0, 0.0, 0.2, 0.6),
                                                  Pt3050_Medium = cms.vdouble(0.3, 0.2, 0.7, 0.8),
                                                  Pt1020_Tight = cms.vdouble(-0.2, 0.2, 0.2, 0.6),
                                                  Pt2030_Medium = cms.vdouble(0.2, 0.2, 0.5, 0.7),
                                                  Pt010_Tight = cms.vdouble(0.5, 0.6, 0.6, 0.9),
                                                  Pt1020_Loose = cms.vdouble(-0.4, -0.4, -0.4, 0.4),
                                                  Pt010_Medium = cms.vdouble(0.2, 0.4, 0.2, 0.6),
                                                  Pt1020_Medium = cms.vdouble(-0.3, 0.0, 0.0, 0.5),
                                                  Pt010_Loose = cms.vdouble(0.0, 0.0, 0.0, 0.2),
                                                  Pt3050_Loose = cms.vdouble(0.0, 0.0, 0.6, 0.2),
                                                  Pt3050_Tight = cms.vdouble(0.5, 0.4, 0.8, 0.9)
                                                ),
                                                impactParTkThreshold = cms.double(0.0)
                                              )),
                                            inputIsCorrected = cms.bool(True),
                                            vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            produceJetIds = cms.bool(True),
                                            jec = cms.string('AK4PF'),
                                            residualsFromTxt = cms.bool(False),
                                            applyJec = cms.bool(True),
                                            jetids = cms.InputTag("pileupJetIdCalculator"),
                                            rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                            jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
                                            runMvas = cms.bool(True)
                                          )
  process.mvaMetSequence  = cms.Sequence(process.leptonPreSelectionSequence +
                                         process.ak4PFJets + process.calibratedAK4PFJetsForPFMVAMEt +
                                         process.puJetIdForPFMVAMEt +
                                         process.mvaMETDiTau + process.mvaMETTauMu + process.mvaMETTauEle + process.mvaMETMuEle + 
                                         process.mvaMETMuMu + process.mvaMETEleEle)
  
else:
  from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
  process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
  process.ak4PFJetsCHS = ak4PFJets.clone(src = 'chs', doAreaFastjet = True)

  process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

  if runOnData:
    process.calibratedAK4PFJetsCHSForPFMVAMEt = cms.EDProducer("PFJetCorrectionProducer",
                                                               src = cms.InputTag("ak4PFJetsCHS"),
                                                               correctors = cms.vstring('ak4PFCHSL1FastL2L3Residual')
                                                             )
  else:
    process.calibratedAK4PFJetsCHSForPFMVAMEt = cms.EDProducer("PFJetCorrectionProducer",
                                                               src = cms.InputTag("ak4PFJetsCHS"),
                                                               correctors = cms.vstring('ak4PFCHSL1FastL2L3')
                                                             )

  from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId
  process.puJetIdForPFMVAMEt = cms.EDProducer("PileupJetIdProducer",
                                              **pileupJetId.parameters_())
  process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
  process.puJetIdForPFMVAMEt.jets = cms.InputTag("calibratedAK4PFJetsCHSForPFMVAMEt")

  process.mvaMetSequence  = cms.Sequence(process.leptonPreSelectionSequence +
                                         process.chs + process.ak4PFJetsCHS + process.calibratedAK4PFJetsCHSForPFMVAMEt +
                                         process.puJetIdForPFMVAMEt +
                                         process.mvaMETDiTau + process.mvaMETTauMu + process.mvaMETTauEle + process.mvaMETMuEle + 
                                         process.mvaMETMuMu + process.mvaMETEleEle)


# END Pairwise MVA MET ==============================================================



#####################################################################################

# NTuple Maker =======================================================================

process.initroottree = cms.EDAnalyzer("InitAnalyzer",
IsData = cms.untracked.bool(isData),
#IsData = cms.untracked.bool(False),
GenParticles = cms.untracked.bool(True)
)

process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(not isData),
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
JetCollectionTag = cms.InputTag("patJetsReapplyJEC::TreeProducer"),
#JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs::RECO" if isPRv4 else "slimmedMETs::PAT"),
MetCovMatrixTag = cms.InputTag("METSignificance:METCovariance:TreeProducer"),
MetSigTag = cms.InputTag("METSignificance:METSignificance:TreeProducer"),
MetCorrCovMatrixTag = cms.InputTag("METCorrSignificance:METCovariance:TreeProducer"),
MetCorrSigTag = cms.InputTag("METCorrSignificance:METSignificance:TreeProducer"),
MetCorrCollectionTag = cms.InputTag("slimmedMETs::TreeProducer"),
#MetCorrCollectionTag = cms.InputTag("slimmedMETsNoHF"),
PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi"),
MvaMetCollectionsTag = cms.VInputTag("mvaMETDiTau", "mvaMETTauMu", "mvaMETTauEle", "mvaMETMuEle", "mvaMETMuMu", "mvaMETEleEle"),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_IsoMu17_eta2p1_v',
'HLT_IsoMu18_v', 
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v',
'HLT_IsoMu22_v', #new
#eltau
'HLT_Ele22_eta2p1_WP75_Gsf_v',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v',

'HLT_Ele32_eta2p1_WPTight_Gsf_v',
'HLT_Ele23_WPLoose_Gsf_v',


#tautau
'HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v',
'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v',

#emu
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v', #new


'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',  #exact hlt filters not available
'HLT_IsoMu17_eta2p1_v',
'HLT_IsoMu24_eta2p1_v',
'HLT_IsoMu18_v',  #new
'HLT_IsoMu20_v', #new
'HLT_IsoMu27_v',
#'HLT_IsoMuTk24_eta2p1_v', no existing in 25ns
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v', #new
'HLT_Ele22_eta2p1_WPLoose_Gsf_v', ##exact hlt filters not available
'HLT_Ele22_eta2p1_WPTight_Gsf_v', 
'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
'HLT_Ele27_eta2p1_WP75_Gsf_v'
		),
TriggerProcess = cms.untracked.string("HLT"),
# tracks
RecTrackPtMin = cms.untracked.double(0.5),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(8.),
RecMuonEtaMax = cms.untracked.double(2.5),
RecMuonHLTriggerMatching = cms.untracked.vstring(


'HLT_IsoMu17_eta2p1_v.*:hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu18_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09,' ,
#mu leg
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',

'HLT_IsoMu22_v.*:hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09',

#mu
'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',

'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4', #new
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2',
'HLT_IsoMu18_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09',
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09',
'HLT_IsoMu27_v.*:hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09', 
'HLT_IsoMu17_eta2p1_v.*:hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu24_eta2p1_v.*:hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09'
),
RecMuonNum = cms.untracked.int32(0),
# photons
RecPhotonPtMin = cms.untracked.double(20.),
RecPhotonEtaMax = cms.untracked.double(10000.),
RecPhotonHLTriggerMatching = cms.untracked.vstring(),
RecPhotonNum = cms.untracked.int32(0),
# electrons
RecElectronPtMin = cms.untracked.double(8.),
RecElectronEtaMax = cms.untracked.double(2.5),
RecElectronHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20', 
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltEle22WPLooseL1IsoEG20erTau20erGsfTrackIsoFilter', 

'HLT_Ele32_eta2p1_WPTight_Gsf_v.*:hltEle32WPTightGsfTrackIsoFilter',
'HLT_Ele23_WPLoose_Gsf_v.*:hltEle23WPLooseGsfTrackIsoFilter',

'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter', 

'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',

'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20',


'HLT_Ele22_eta2p1_WPLoose_Gsf_v.*:hltSingleEle22WPLooseGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WPTight_Gsf_v.*:hltSingleEle22WPTightGsfTrackIsoFilter',

'HLT_Ele22_eta2p1_WP75_Gsf_v.*:hltSingleEle22WP75GsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WP75_Gsf_v.*:hltEle27WP75GsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WPLoose_Gsf_v.*:hltEle27WPLooseGsfTrackIsoFilter'
#'HLT_Ele32_eta2p1_WPLoose_Gsf_v.*:hltEle32WPLooseGsfTrackIsoFilter', 
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(15),
RecTauEtaMax = cms.untracked.double(2.5),                                      
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20',

#tau leg
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIsoAgainstMuon',

'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v.*:hltAK4PFJetsForTaus', #need to check
'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v.*:hltDoublePFTau35TrackPt1MediumIsolationDz02Reg',
'HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v.*:hltDoublePFTau40TrackPt1MediumIsolationDz02Reg'
),


RecTauFloatDiscriminators = cms.untracked.vstring(
#'againstElectronLoose',
#'againstElectronLooseMVA5',
#'againstElectronMVA5category',
#'againstElectronMVA5raw',
#'againstElectronMedium',
#'againstElectronMediumMVA5',
#'againstElectronTight',
#'againstElectronTightMVA5',
#'againstElectronVLooseMVA5',
#'againstElectronVTightMVA5',
#'againstMuonLoose',
#'againstMuonLoose2',
#'againstMuonLoose3',
#'againstMuonLooseMVA',
#'againstMuonMVAraw',
#'againstMuonMedium',
#'againstMuonMedium2',
#'againstMuonMediumMVA',
#'againstMuonTight',
#'againstMuonTight2',
#'againstMuonTight3',
#'againstMuonTightMVA',
#'byCombinedIsolationDeltaBetaCorrRaw3Hits',
#'byIsolationMVA3newDMwLTraw',
#'byIsolationMVA3newDMwoLTraw',
#'byIsolationMVA3oldDMwLTraw',
#'byIsolationMVA3oldDMwoLTraw',
#'byLooseCombinedIsolationDeltaBetaCorr3Hits',
#'byLooseIsolationMVA3newDMwLT',
#'byLooseIsolationMVA3newDMwoLT',
#'byLooseIsolationMVA3oldDMwLT',
#'byLooseIsolationMVA3oldDMwoLT',
#'byMediumCombinedIsolationDeltaBetaCorr3Hits',
#'byMediumIsolationMVA3newDMwLT',
#'byMediumIsolationMVA3newDMwoLT',
#'byMediumIsolationMVA3oldDMwLT',
#'byMediumIsolationMVA3oldDMwoLT',
#'byTightCombinedIsolationDeltaBetaCorr3Hits',
#'byTightIsolationMVA3newDMwLT',
#'byTightIsolationMVA3newDMwoLT',
#'byTightIsolationMVA3oldDMwLT',
#'byTightIsolationMVA3oldDMwoLT',
#'byVLooseIsolationMVA3newDMwLT',
#'byVLooseIsolationMVA3newDMwoLT',
#'byVLooseIsolationMVA3oldDMwLT',
#'byVLooseIsolationMVA3oldDMwoLT',
#'byVTightIsolationMVA3newDMwLT',
#'byVTightIsolationMVA3newDMwoLT',
#'byVTightIsolationMVA3oldDMwLT',
#'byVTightIsolationMVA3oldDMwoLT',
#'byVVTightIsolationMVA3newDMwLT',
#'byVVTightIsolationMVA3newDMwoLT',
#'byVVTightIsolationMVA3oldDMwLT',
#'byVVTightIsolationMVA3oldDMwoLT',
#'chargedIsoPtSum',
#'decayModeFinding',
#'decayModeFindingNewDMs',
#'neutralIsoPtSum',
#'puCorrPtSum'
),
RecTauBinaryDiscriminators = cms.untracked.vstring(),
RecTauNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(18.),
RecJetEtaMax = cms.untracked.double(5.2),
RecJetHLTriggerMatching = cms.untracked.vstring(),
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
  srcMet = cms.InputTag('slimmedMETs::TreeProducer')
)

#################################################################################

process.p = cms.Path(
  process.initroottree*
  process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
  process.METSignificance * process.METCorrSignificance *
  process.mvaMetSequence *
  process.egmGsfElectronIDSequence * 
  #process.HBHENoiseFilterResultProducer* #produces HBHE bools baseline
  #process.ApplyBaselineHBHENoiseFilter*  #reject events based 
  #process.ApplyBaselineHBHEISONoiseFilter*  #reject events based -- disable the module, performance is being investigated fu
  process.makeroottree
)

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("/nfs/dust/cms/user/alkaloge/TauAnalysis/new/CMSSW_7_4_6/src/DesyTauAnalyses/NTupleMaker/test/Ntuple74.root"),
                                   #fileName = cms.string("/nfs/dust/cms/user/alkaloge/TauAnalysis/new/CMSSW_7_2_3_patch1/src/DesyTauAnalyses/NTupleMaker/test/Staus/${1}_NTuple.root"),
                                   fileName = cms.string("output.root")
                               	)

#process.end = cms.EndPath(process.Out*process.TFileService)

#processDumpFile = open('MyRootMaker.dump', 'w')
#print >> processDumpFile, process.dumpPython()

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
