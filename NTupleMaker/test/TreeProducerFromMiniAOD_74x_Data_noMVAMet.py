import FWCore.ParameterSet.Config as cms

isData = True
#skim = 0
year = 2015
period = 'Run2015B'
#sampleName = 'MonteCarlo'
# sampleName = 'TTJets', "QCD", "DYJetsToLL_M50"

#if "@SKIM@".lower()=="0":
#    skim = 0

#sampleName = "@SAMPLE_NAME@"

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

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V7C', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v0', '')

#process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9::Al')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
		# Single Muon Run2015B
#		'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/164/00000/544E58CD-A226-E511-834E-02163E0134D6.root',
#		'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/167/00000/EE9594B8-A826-E511-A18B-02163E011A7D.root',
#		'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/168/00000/4E8E390B-EA26-E511-9EDA-02163E013567.root',
#		'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/168/00000/60FF8405-EA26-E511-A892-02163E01387D.root'
		# Single Electron Run2015B
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/096/00000/22D22D7F-5626-E511-BDE3-02163E011FAB.root',
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/161/00000/7019DC27-9C26-E511-84FF-02163E011CC2.root',
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/163/00000/9C435096-9F26-E511-A1D7-02163E012AB6.root',
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/164/00000/4633CC68-A326-E511-95D0-02163E0124EA.root',
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/167/00000/E661FBF2-A726-E511-8B8B-02163E01207C.root',
#		'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/168/00000/FCB6CB61-BC26-E511-8858-02163E01375B.root'
                # Double EG Run2015B
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/096/00000/8A2D533C-5626-E511-AF3C-02163E011FAB.root',
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/161/00000/AC857A3B-9C26-E511-B32E-02163E012704.root',
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/163/00000/9CB965BF-9F26-E511-8FB1-02163E012040.root',
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/164/00000/2CA2349A-A526-E511-A371-02163E0138D0.root',
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/167/00000/C446EC67-A726-E511-81C1-02163E0119E7.root',
		'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/168/00000/747F782A-BB26-E511-BA24-02163E011EE9.root'
		# MuonEG Run2015B
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/162/00000/427789C7-4127-E511-979F-02163E012096.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/164/00000/0270BD14-A426-E511-9857-02163E01182A.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/167/00000/F668ECD1-A726-E511-901F-02163E011D57.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/168/00000/F8F3DF28-BC26-E511-B98D-02163E0135AC.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/244/00000/10C5D0A8-7527-E511-BFE2-02163E0140E1.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/251/00000/88361DDE-7827-E511-94C3-02163E0127DF.root',
#		'/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1/000/251/252/00000/E09EA71E-9D27-E511-B478-02163E014558.root'
        )		    
)
 
#####################################################
  

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")
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

"""
process.skimOutputModule = cms.OutputModule("PoolOutputModule",                                 
    ##goldenZmumuEventContent,
     outputCommands = cms.untracked.vstring(                                    
        'drop *',
        'keep *_l1extraParticles_*_*'
    )  
)
"""


process.tauPreSelectionDiTau = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 40. && abs(eta) < 2.5 && tauID("decayModeFindingNewDMs") > 0.5')
)
 
 
process.tauPreSelectionTauEle = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFinding") > 0.5')
)
 
 
process.tauPreSelectionTauMu = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFinding") > 0.5')
)

process.muonPreSelectionMuEle = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 9. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)
 
 
process.muonPreSelectionTauMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 20. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)

process.electronPreSelectionMuEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 13. && abs(eta) < 2.5')
)
 
 
process.electronPreSelectionTauEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 20. && abs(eta) < 2.5')
)

process.leptonSkimSequence = cms.Sequence(process.electronPreSelectionMuEle + process.electronPreSelectionTauEle + process.muonPreSelectionTauMu + process.muonPreSelectionMuEle + process.tauPreSelectionTauMu + process.tauPreSelectionTauEle + process.tauPreSelectionDiTau)

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

process.mvaMETDiTau = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionDiTau"), cms.InputTag("tauPreSelectionDiTau")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETTauMu = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauMu"), cms.InputTag("muonPreSelectionTauMu")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETTauEle = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauEle"), cms.InputTag("electronPreSelectionTauEle")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETMuEle = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("muonPreSelectionMuEle"), cms.InputTag("electronPreSelectionMuEle")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.patMvaMetDiTau = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETDiTau'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetTauMu = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETTauMu'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetTauEle = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETTauEle'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetMuEle = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETMuEle'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )

process.mvaMetSequence  = cms.Sequence(process.leptonSkimSequence* process.ak4PFJets * process.calibratedAK4PFJetsForPFMVAMEt * 
                                       process.puJetIdForPFMVAMEt * process.mvaMETDiTau * process.mvaMETTauMu * 
                                       process.mvaMETTauEle * process.mvaMETMuEle * process.patMvaMetDiTau * 
                                       process.patMvaMetTauMu * process.patMvaMetTauEle * process.patMvaMetMuEle)

process.pileupJetIdFull = cms.EDProducer("PileupJetIdProducer",
    produceJetIds = cms.bool(True),
    runMvas = cms.bool(True),
    inputIsCorrected = cms.bool(False),
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    residualsTxt = cms.FileInPath('RecoJets/JetProducers/data/download.url'),
    jec = cms.string('AK4PF'),
    residualsFromTxt = cms.bool(False),
    applyJec = cms.bool(True),
    jetids = cms.InputTag(""),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    jets = cms.InputTag("ak4PFJets"),
    algos = cms.VPSet(cms.PSet(
        tmvaVariables = cms.vstring('nvtx', 
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
        tmvaMethod = cms.string('JetIDMVAHighPt'),
        cutBased = cms.bool(False),
        tmvaWeights = cms.string('CondFormats/JetMETObjects/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),
        #tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml.gz'),
        tmvaSpectators = cms.vstring('jetPt', 
            'jetEta', 
            'jetPhi'),
        label = cms.string('full53x'),
        version = cms.int32(-1),
        JetIdParams = cms.PSet(
            Pt2030_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42),
            Pt2030_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
            Pt3050_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
            Pt1020_MET = cms.vdouble(0.3, -0.2, -0.4, -0.4),
            Pt2030_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
            Pt010_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
            Pt1020_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
            Pt3050_MET = cms.vdouble(0.0, 0.0, -0.1, -0.2),
            Pt010_MET = cms.vdouble(0.0, -0.6, -0.4, -0.4),
            Pt1020_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
            Pt010_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
            Pt1020_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
            Pt2030_MET = cms.vdouble(0.0, 0.0, 0.0, 0.0),
            Pt010_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
            Pt3050_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
            Pt3050_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42)
        ),
        impactParTkThreshold = cms.double(1.0)
    ))
)
process.puJetIdSequence = cms.Sequence(process.pileupJetIdFull)
'''
process.load("PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff")
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
massSearchReplaceAnyInputTag(process.makePatJets,
                             "ak4PFJetsCHS",
                             "ak4PFJets",
                             verbose=True)
massSearchReplaceAnyInputTag(process.makePatJets,
                             'offlinePrimaryVertices',
                             'offlineSlimmedPrimaryVertices',
                             verbose=True)
process.patJetCorrFactors.payload = cms.string('AK4PF')
process.puJetIdSequence += process.makePatJets 
'''
##################################################
# Main
process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(True),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),
RecMvaMet = cms.untracked.bool(True),                                      
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(True),
RecTau = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
ElectronCollectionTag = cms.InputTag("slimmedElectrons"),
TauCollectionTag = cms.InputTag("slimmedTaus"),
JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs"),
MvaMetCollectionsTag = cms.VInputTag("patMvaMetDiTau", "patMvaMetTauMu", "patMvaMetTauEle", "patMvaMetMuEle"),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
'HLT_IsoMu24_eta2p1_v',
'HLT_IsoMu27_v',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v',
'HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v',
'HLT_Ele32_eta2p1_WP75_Gsf_v',
'HLT_Ele22_eta2p1_WPLoose_Gsf_v',
'HLT_Ele22_eta2p1_WPTight_Gsf_v',
'HLT_Ele23_WPLoose_Gsf_v',
'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
'HLT_Ele32_eta2p1_WPLoose_Gsf_v',
'HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v'
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
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8', 
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v.*:hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',
'HLT_IsoMu24_eta2p1_v.*:hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09', 
'HLT_IsoMu27_v.*:hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09' 
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
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter', 
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter',
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*:hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter',
'HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v.*:hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltSingleEle22WPLooseGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20',
'HLT_Ele22_eta2p1_WPLoose_Gsf_v.*:hltSingleEle22WPLooseGsfTrackIsoFilter',
'HLT_Ele22_eta2p1_WPTight_Gsf_v.*:hltSingleEle22WPTightGsfTrackIsoFilter',
'HLT_Ele23_WPLoose_Gsf_v.*:hltEle23WPLooseGsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WP75_Gsf_v.*:hltEle32WP75GsfTrackIsoFilter',
'HLT_Ele27_eta2p1_WPLoose_Gsf_v.*:hltEle27WPLooseGsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WPLoose_Gsf_v.*:hltEle32WPLooseGsfTrackIsoFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(18),
RecTauEtaMax = cms.untracked.double(2.4),                                      
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20', 
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIso',
'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltPFTau20TrackLooseIsoAgainstMuon',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',
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
'jetBProbabilityBJetTags',
'jetProbabilityBJetTags',
'trackCountingHighPurBJetTags',
'trackCountingHighEffBJetTags',
'simpleSecondaryVertexHighEffBJetTags',
'simpleSecondaryVertexHighPurBJetTags',
'combinedInclusiveSecondaryVertexV2BJetTags',
'pfCombinedSecondaryVertexBJetTags',
'combinedMVABJetTags'
),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("Data") 
)

#process.patJets.addBTagInfo = cms.bool(True)


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
	process.ak4PFJets*
	process.mvaMetSequence*
	process.puJetIdSequence*
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



