import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeProducer")

runOnMCSample = True

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        #'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/125A6B71-C56A-E411-9D2B-0025907609BE.root'
        #'file:copyEvent/pickevents_merged.root'
    ),
##    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
)
 
#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# compute event weights for pile-up reweighting
# (Summer'12 MC to 2012 run A data)

##from TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi import vertexMultiplicityReweight
##process.vertexMultiplicityReweight3d2012RunABC = vertexMultiplicityReweight.clone(
##    inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs190456to208686_Mu17_Mu8.root"),
##    type = cms.string("gen3d"),
##    mcPeriod = cms.string("Summer12_S10")
##)
##process.producePFTauIdEffNtuple2Sequence += process.vertexMultiplicityReweight3d2012RunABC
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
#process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
#jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(
#    process,
#    jetSource = cms.InputTag('ak5PFJets'),
#    jetCorrections = ( 'AK5PF', jetCorrections, "" ),
#    outputModules = []
#)

from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = ""

if runOnMCSample:
	usePF2PAT(process,runPF2PAT=True,
		jetAlgo='AK4', runOnMC=True, postfix=postfix,
		jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
		pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
		typeIMetCorrections=True
	)
else:
	usePF2PAT(process,runPF2PAT=True,
		jetAlgo='AK4', runOnMC=False, postfix=postfix,
		jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']),
		pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
#		typeIMetCorrections=True
	)

getattr(process, 'pfPileUp'+postfix).checkClosestZVertex = False

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('patTuple.root'),
        # save only events passing the full path
        #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
        # save PAT Layer 1 output; you need a '*' to
        # unpack the list of commands 'patEventContent'
        outputCommands = cms.untracked.vstring('keep *')
)

# switch on PAT trigger                                                                                                                      
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )

process.makePatTrigger = cms.Sequence(process.patTrigger*process.patTriggerEvent)

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
	"PrimaryVertexObjectFilter",
	filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
	src=cms.InputTag('offlinePrimaryVertices')
)

##################################################
# Main
process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(False),
Year = cms.untracked.uint32(2015),
Period = cms.untracked.string("PHYS14"),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(True),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(False),
RecTau = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("patMuons"), 
ElectronCollectionTag = cms.InputTag("patElectrons"),
TauCollectionTag = cms.InputTag("patTaus"),
JetCollectionTag = cms.InputTag("patJets"),
MetCollectionTag = cms.InputTag("patMETs"),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("genParticles"),
TriggerObjectCollectionTag = cms.InputTag("patTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlinePrimaryVertices"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_IsoMu24_eta2p1_IterTrk02_v',
'HLT_IsoTkMu24_eta2p1_IterTrk02_v',
'HLT_Ele27_eta2p1_WP85_Gsf_v',
'HLT_Ele32_eta2p1_WP85_Gsf_v',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v',
'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v',
'HLT_PFMET170_NoiseCleaned_v'
),
TriggerProcess = cms.untracked.string("HLT"),
# tracks
RecTrackPtMin = cms.untracked.double(0.5),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(10.),
RecMuonEtaMax = cms.untracked.double(2.4),
RecMuonHLTriggerMatching = cms.untracked.vstring(
'HLT_IsoMu24_eta2p1_IterTrk02_v.*:hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02',
'HLT_IsoTkMu24_eta2p1_IterTrk02_v.*:hltL3fL1sMu20L1Eta2p1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltL1Mu12EG7L3IsoMuFiltered23',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8'
),
RecMuonNum = cms.untracked.int32(0),
# photons
RecPhotonPtMin = cms.untracked.double(20.),
RecPhotonEtaMax = cms.untracked.double(10000.),
RecPhotonHLTriggerMatching = cms.untracked.vstring(),
RecPhotonNum = cms.untracked.int32(0),
# electrons
RecElectronPtMin = cms.untracked.double(10.),
RecElectronEta = cms.untracked.double(2.4),
RecElectronHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele27_eta2p1_WP85_Gsf_v.*:hltEle27WP85GsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WP85_Gsf_v.*:hltEle32WP85GsfTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(20),
RecTauEtaMax = cms.untracked.double(2.3),                                      
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v.*:hltPFTau50TrackPt30LooseAbsOrRelIso'
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
RecJetPtMin = cms.untracked.double(30.),
RecJetEtaMax = cms.untracked.double(2.5),
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
SampleName = cms.untracked.string("MC") 
)

# 3. Add them all to the sequence
process.patseq = cms.Sequence(
	process.goodOfflinePrimaryVertices 
#	getattr(process,"patPF2PATSequence"+postfix)
)

process.options = cms.untracked.PSet()
process.options.allowUnscheduled = cms.untracked.bool(True)

process.p = cms.Path(#process.producePFTauIdEffNtuple2Sequence + 
#                     process.particleFlowPtrs +
#                     process.makePatElectrons +
#                     process.makePatMuons + 
#                     process.makePatJets +
#                    process.patDefaultSequence + 
                     process.patseq *  
		     process.makePatTrigger *
                     process.makeroottree
                     )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root")
)

#processDumpFile = open('producePFTauIdEffNtuple2.dump', 'w')
#print >> processDumpFile, process.dumpPython()




