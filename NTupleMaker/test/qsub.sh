#!/bin/sh 
# $1 - filename
# $2 - jobId
rm stau100_LSP0_$1.py
cat > stau100_LSP0_$1.py <<EOF1
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# The following three lines reduce the clutter of repeated printouts
# of the same exception message.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.destinations = ['cerr']
process.MessageLogger.statistics = []
process.MessageLogger.fwkJobReports = []
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(
          'file:/nfs/dust/cms/user/alkaloge/new_HtoTauTau/CMSSW_5_3_9/src/H2to2H1to4Taus/GenInfo/test/LHE/$1'
	)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('output_$2.root')
)
process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('runs Z2* Pythia6'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/GenProduction/python/EightTeV/Hadronizer_TuneZ2star_8TeV_generic_LHE_pythia_tauola_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Hadronizer_TuneZ2star_8TeV_generic_LHE_pythia_tauola_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START52_V9::All'

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            UseTauolaPolarization = cms.bool(True),
            InputCards = cms.PSet(
                mdtau = cms.int32(0),
                pjak2 = cms.int32(0),
                pjak1 = cms.int32(0)
            )
        ),
        parameterSets = cms.vstring('Tauola')
    ),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(8000.0),
    UseExternalGenerators = cms.untracked.bool(True),
    jetMatching = cms.untracked.PSet(
        MEMAIN_showerkt = cms.double(0),
        MEMAIN_maxjets = cms.int32(1),
        MEMAIN_minjets = cms.int32(0),
        MEMAIN_qcut = cms.double(10),
        MEMAIN_excres = cms.string(''),
        MEMAIN_etaclmax = cms.double(5),
        MEMAIN_nqmatch = cms.int32(3),
        outTree_flag = cms.int32(1),
        scheme = cms.string('Madgraph'),
        mode = cms.string('auto')
    ),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
	    'MDCY(C130,1)=1   ! decay k0-longs',
            'MDCY(C211,1)=1   ! decay pions',
            'MDCY(C321,1)=1   ! decay kaons',
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0          ! User defined processes', 
            'PMAS(5,1)=4.4   ! b quark mass', 
            'PMAS(6,1)=172.4 ! t quark mass'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                     maxEventsToPrint = cms.untracked.int32(-1),
                                     printVertex = cms.untracked.bool(False),
                                     src = cms.InputTag("genParticles")
                                   )
	
process.doublemufilter = cms.EDFilter("MCParticlePairFilter",
                                      Status = cms.untracked.vint32(1, 1),
                                      MinDeltaPhi = cms.untracked.double(0.0),
                                      MaxDeltaPhi = cms.untracked.double(6.29),
                                      MinPt = cms.untracked.vdouble(15.0, 8.0),
                                      MinP = cms.untracked.vdouble(15.0, 8.0),
                                      MaxEta = cms.untracked.vdouble(2.4, 2.4),
                                      MinEta = cms.untracked.vdouble(-2.4, -2.4),
                                      ParticleCharge = cms.untracked.int32(1),
                                      MaxInvMass = cms.untracked.double(1000.0),
                                      MinInvMass = cms.untracked.double(14.0),
                                      ParticleID1 = cms.untracked.vint32(13),
                                      ParticleID2 = cms.untracked.vint32(13)
)
process.analysis = cms.EDAnalyzer("GenInfo",
	GenParticleSource = cms.InputTag("genParticles"),
	generatorInfoSource =cms.InputTag("generator")

) 
process.eventCountPostFilter  = cms.EDAnalyzer("EventCount",
					       DoMC = cms.bool(False),
					       GenParticleSource = cms.InputTag("genParticles"))

process.eventCountPreFilter  = cms.EDAnalyzer("EventCount",
					      DoMC = cms.bool(False),
					      GenParticleSource = cms.InputTag("genParticles"))



# Path and EndPath definitions
#process.generation_step = cms.Path(process.pgen*process.analysis)
process.generation_step = cms.Path(process.eventCountPreFilter*process.pgen*process.doublemufilter*process.eventCountPostFilter*process.analysis)
#process.generation_step = cms.Path(process.pgen*process.analysis)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.RAWSIMoutput_step = cms.EndPath()

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step)
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
#process.schedule = cms.Schedule(process.generation_step,process.endjob_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

EOF1
rm stau100_LSP0_$1.zsh
cat > stau100_LSP0_$1.zsh <<EOF
#!/bin/zsh
#
#(make sure the right shell will be used)
#$ -S /bin/zsh
#
#(the cpu time for this job)
#$ -l h_cpu=0:29:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=2000M
#
#(use hh site)
#$ -l site=hh 
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL5
#$ -l os=sld5
#
# use current dir and current environment
#$ -cwd
#$ -V
#
#$ -o stau100_LSP0_$1.out
#
#$ -e stau100_LSP0_$1.err
mkdir $1
cmsRun stau100_LSP0_$1.py

EOF

rm stau100_LSP0_$1.out
chmod u+x stau100_LSP0_$1.zsh
qsub stau100_LSP0_$1.zsh