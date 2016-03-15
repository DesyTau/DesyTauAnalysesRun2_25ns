#!/usr/bin/env python

import os
import re
import shlex
import string
import subprocess

#WorkdirLoc = '/nfs/dust/cms/user/rasp/ntuples/'
#WorkdirData = '/nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet/'
#WorkdirLoc = '/nfs/dust/cms/user/fcost/store/DataCards_MVAMEt_5X/' #'/nfs/dust/cms/user/fcost/store/DataCards/'
WorkdirData = '/nfs/dust/cms/user/rasp/ntuples/Data2015D_MVAMEt_5X/'
#WorkdirMC = '/nfs/dust/cms/user/fcost/store/DataCards_JECv6_SynchedMVAMEt/'
#WorkdirMC = '/nfs/dust/cms/user/fcost/store/76x_MiniAODv2/'
WorkdirMC = '/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/'
#WorkdirLoc = '/nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_Spring15_50ns_v1/'
OutDir     = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples76_25ns_v2/'    

options = {
    ###tau+tau samples

    ##DATA
#    'Tau_Run2015D_05Oct2015' : {
#    'inputFilePath'  : WorkdirData+'Tau_Run2015D_05Oct2015/',
#    'outputFileName' : OutDir+'nTupleRun2015D-05Oct2015-Data.root',
#    'sample'         : 'Run2015D-05Oct2015-Data',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 0,
#    'skimEff'        : 0,
#    'iJson'          : 0,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'Tau_Run2015D_PRv4' : {
#    'inputFilePath'  : WorkdirData+'Tau_Run2015D_PRv4/',
#    'outputFileName' : OutDir+'nTupleRun2015D-PRv4-Data.root',
#    'sample'         : 'Run2015D-PRv4-Data',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 0,
#    'skimEff'        : 0,
#    'iJson'          : 0,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    ##Bkg MC
#    'DYJets_TauTau' : {
#    'inputFilePath'  : WorkdirMC+'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/',
#    'outputFileName' : OutDir+'nTupleDYJets_TauTau.root',
#    'sample'         : 'DYJets_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 6025.2,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WJetsToLNu' : {
#    'inputFilePath'  : WorkdirMC+'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/',    
#    'outputFileName' : OutDir+'nTupleWJets_TauTau.root',
#    'sample'         : 'WJets_TauTau',
#    'xSection'       : 61526.7,
#    'isMCAtNLO'      : 0,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'TTJets' : {
#    'inputFilePath'  : WorkdirMC+'TT_TuneCUETP8M1_13TeV-powheg-pythia8/',
#    'outputFileName' : OutDir+'nTupleTTJets_TauTau.root',
#    'sample'         : 'TTJets_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 831.76,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'SingleTop_t' : {
#    'inputFilePath'  : WorkdirMC+'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/',
#    'outputFileName' : OutDir+'nTupleSTopT_TauTau.root',
#    'sample'         : 'STopT_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 136.05*0.108*3,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'SingleAntiTop_t' : {
#    'inputFilePath'  : WorkdirMC+'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/',
#    'outputFileName' : OutDir+'nTupleSAntiTopT_TauTau.root',
#    'sample'         : 'SAntiTopT_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 80.95*0.108*3,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'SingleTop_tW' : {
#    'inputFilePath'  : WorkdirMC+'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/',
#    'outputFileName' : OutDir+'nTupleSTopTW_TauTau.root',
#    'sample'         : 'STopTW_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 35.6,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'SingleAntiTop_tW' : {
#    'inputFilePath'  : WorkdirMC+'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/',
#    'outputFileName' : OutDir+'nTupleSAntiTopTW_TauTau.root',
#    'sample'         : 'SAntiTopTW_TauTau',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 35.6,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'VVTo2L2Nu' : {
#    'inputFilePath'  : WorkdirMC+'VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleVVTo2L2Nu_TauTau.root',
#    'sample'         : 'VVTo2L2Nu_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 11.95,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WZJets' : {
#    'inputFilePath'  : WorkdirMC+'WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/',
#    'outputFileName' : OutDir+'nTupleWZJets_TauTau.root',
#    'sample'         : 'WZJets_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 5.26,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WWToLNuQQ' : {
#    'inputFilePath'  : WorkdirMC+'WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleWWToLNuQQ_TauTau.root',
#    'sample'         : 'WWToLNuQQ_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 49.997,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WZTo1L1Nu2Q' : {
#    'inputFilePath'  : WorkdirMC+'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleWZTo1L1Nu2Q_TauTau.root',
#    'sample'         : 'WZTo1L1Nu2Q_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 10.71,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WZTo1L3Nu' : {
#    'inputFilePath'  : WorkdirMC+'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleWZTo1L3Nu_TauTau.root',
#    'sample'         : 'WZTo1L3Nu_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 3.05,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'WZTo2L2Q' : {
#    'inputFilePath'  : WorkdirMC+'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleWZTo2L2Q_TauTau.root',
#    'sample'         : 'WZTo2L2Q_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 5.595,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'ZZTo2L2Q' : {
#    'inputFilePath'  : WorkdirMC+'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/',
#    'outputFileName' : OutDir+'nTupleZZTo2L2Q_TauTau.root',
#    'sample'         : 'ZZTo2L2Q_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 3.22,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'ZZTo4L' : {
#    'inputFilePath'  : WorkdirMC+'ZZTo4L_13TeV-amcatnloFXFX-pythia8/',
#    'outputFileName' : OutDir+'nTupleZZTo4L_TauTau.root',
#    'sample'         : 'ZZTo4L_TauTau',
#    'isMCAtNLO'      : 1,
#    'xSection'       : 1.212,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    ##Higgs MC
#    'GGFH125' : {
#    'inputFilePath'  : WorkdirMC+'GluGluHToTauTau_M125_13TeV_powheg_pythia8/',
#    'outputFileName' : OutDir+'nTupleGGFH125_TauTau.root',
#    'sample'         : 'GGFH125',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 43.92*0.0632,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,                                                      
#    'nDiv'           : 1
#    },
#    #'VBFH125' : {
#    #'inputFilePath'  : WorkdirMC+'VBF_HToTauTau_M-125_13TeV-powheg-pythia6',
#    #'outputFileName' : OutDir+'nTupleVBFH125_TauTau.root',
#    #'sample'         : 'VBFH125',
#    #'isMCAtNLO'      : 0,
#    #'xSection'       : 3.748*0.0632,
#    #'skimEff'        : 1.0,
#    #'iJson'          : -1,
#    #'iDiv'           : 0,
#    #'nDiv'           : 1
#    #},
    ######MSSM
    'SUSYGGH160' : {
    'inputFilePath'  : WorkdirMC+'SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/',
    'outputFileName' : OutDir+'nTupleSUSYGGH160_TauTau.root',
    'sample'         : 'SUSYGGH160',
    'isMCAtNLO'      : 0,
    'xSection'       : 1.0,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
#    'SUSYBBH160' : {
#    'inputFilePath'  : WorkdirMC+'SUSYGluGluToBBHToTauTau_M-160/',
#    'outputFileName' : OutDir+'nTupleSUSYBBH160_TauTau.root',
#    'sample'         : 'SUSYGGH160',
#    'isMCAtNLO'      : 0,
#    'xSection'       : 1.0,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
}

config_template = string.Template('''
import FWCore.ParameterSet.Config as cms

import os

process = cms.PSet()

process.fwliteInput = cms.PSet(
fileNames = cms.vstring(),
maxEvents = cms.int32(-1),
#outputEvery = cms.uint32(1000)
)

inputFilePath = '$inputFilePath'
inputFiles = os.listdir(inputFilePath)
process.fwliteInput.fileNames = cms.vstring([ os.path.join(inputFilePath, inputFile) for inputFile in inputFiles ])

process.fwliteOutput = cms.PSet(
fileName = cms.string('$outputFileName')
)

process.preAnalyzerTauTau = cms.PSet(
sample = cms.string('$sample'),
ismcatnlo = cms.int32($isMCAtNLO), 
analysis = cms.string('$analysis'),
xSection = cms.double($xSection),
skimEff = cms.double($skimEff),
iJson = cms.int32($iJson),
iDiv = cms.int32($iDiv),
nDiv = cms.int32($nDiv)
)
''')

currentDirectory    = os.getcwd()
submissionDirectory = os.path.join(currentDirectory, "Configs")

for sample, option in options.items():
    for analysis in [ 'nominal', 'TauUp', 'TauDown' ]:
        if  re.search("Data",sample)!=None and analysis != 'nominal':
            continue
        configOptions = option.copy()
        configOptions['analysis'] = analysis
        if  re.search("Data",sample)==None :
            configOptions['outputFileName'] = configOptions['outputFileName'].replace('.root', '_%s.root' % analysis)
        configFileName = "preAnalyzerTauTau_Summer15_%s_%s_cfg.py" % (sample,analysis)
        configFileName_full = os.path.join(submissionDirectory, configFileName)
        configFile = open(configFileName_full, 'w')
        configConfig = config_template.substitute(configOptions)
        configFile.write(configConfig)
        configFile.close()
