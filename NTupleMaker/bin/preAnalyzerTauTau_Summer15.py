#!/usr/bin/env python

import os
import re
import shlex
import string
import subprocess

#WorkdirLoc = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/Sync2015/'
#WorkdirLoc = '/nfs/dust/cms/user/rasp/ntuples/'
#WorkdirLoc = '/nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_Spring15_25ns_v1/'
#OutDir     = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74/'
WorkdirData = '/nfs/dust/cms/group/susy-desy/Run2/Stau/Data/50ns/v2/'
WorkdirLoc = '/nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_Spring15_50ns_v1/'
OutDir     = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74_50ns/'    

options = {
    ###tau+tau samples

    ##DATA
    'Run2015B-Data_TauTau' : {
    'inputFilePath'  : WorkdirData+'Tau/',
    'outputFileName' : OutDir+'nTupleRun2015B-Data.root',
    'sample'         : 'Run2015B-Data',
    'xSection'       : 0,
    'skimEff'        : 0,
    'iJson'          : 0,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    ##Bkg MC
    'DYJets_TauTau' : {
    'inputFilePath'  : WorkdirLoc+'DYJetsToLL_M-50/',
    'outputFileName' : OutDir+'nTupleDYJets_TauTau.root',
    'sample'         : 'DYJets_TauTau',
    'xSection'       : 6025,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    'WJetsToLNu' : {
    #'inputFilePath'  : WorkdirLoc+'WJetsToLNu/',
    'inputFilePath'  : WorkdirLoc+'WJets/',    
    'outputFileName' : OutDir+'nTupleWJets_TauTau.root',
    'sample'         : 'WJets_TauTau',
    'xSection'       : 61526.7,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    'TTJets' : {
    #'inputFilePath'  : WorkdirLoc+'TTJets_amatnlo/',
    'inputFilePath'  : WorkdirLoc+'TTJets/',
    'outputFileName' : OutDir+'nTupleTTJets_TauTau.root',
    'sample'         : 'TTJets_TauTau',
    'xSection'       : 831.76,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    'SingleTop_t' : {
    'inputFilePath'  : WorkdirLoc+'ST_t-channel/',
    'outputFileName' : OutDir+'nTupleSTopT_TauTau.root',
    'sample'         : 'STopT_TauTau',
    'xSection'       : 136.05*0.108*3,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    #'SingleAntiTop_t' : {
    #'inputFilePath'  : WorkdirLoc+'ST_t-channel_antitop/',
    #'outputFileName' : OutDir+'nTupleSAntiTopT_TauTau.root',
    #'sample'         : 'SAntiTopT_TauTau',
    #'xSection'       : 80.97*0.108*3,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    'SingleTop_tW' : {
    'inputFilePath'  : WorkdirLoc+'ST_tW_top/',
    'outputFileName' : OutDir+'nTupleSTopTW_TauTau.root',
    'sample'         : 'STopTW_TauTau',
    'xSection'       : 35.6,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    'SingleAntiTop_tW' : {
    'inputFilePath'  : WorkdirLoc+'ST_tW_antitop/',
    'outputFileName' : OutDir+'nTupleSAntiTopTW_TauTau.root',
    'sample'         : 'SAntiTopTW_TauTau',
    'xSection'       : 35.6,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    #'WWTo2L2Nu' : {
    #'inputFilePath'  : WorkdirLoc+'WWTo2L2Nu/',
    #'outputFileName' : OutDir+'nTupleWWTo2L2Nu_TauTau.root',
    #'sample'         : 'WWTo2L2Nu_TauTau',
    #'xSection'       : 10.481,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    #'WWTo4Q' : {
    #'inputFilePath'  : WorkdirLoc+'WWTo4Q/',
    #'outputFileName' : OutDir+'nTupleWWTo4Q_TauTau.root',
    #'sample'         : 'WWTo4Q_TauTau',
    #'xSection'       : 45.20,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    #'WWToLNuQQ' : {
    #'inputFilePath'  : WorkdirLoc+'WWToLNuQQ/',
    #'outputFileName' : OutDir+'nTupleWWToLNuQQ_TauTau.root',
    #'sample'         : 'WWToLNuQQ_TauTau',
    #'xSection'       : 43.53,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    #'WZTo1L1Nu2Q' : {
    #'inputFilePath'  : WorkdirLoc+'WZTo1L1Nu2Q/',
    #'outputFileName' : OutDir+'nTupleWZTo1L1Nu2Q_TauTau.root',
    #'sample'         : 'WZTo1L1Nu2Q_TauTau',
    #'xSection'       : 10.96,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    #'WZTo3LNu' : {
    #'inputFilePath'  : WorkdirLoc+'WZTo3LNu/',
    #'outputFileName' : OutDir+'nTupleWZTo3LNu_TauTau.root',
    #'sample'         : 'WZTo3LNu_TauTau',
    #'xSection'       : 4.42965,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
    #'ZZ' : {
    #'inputFilePath'  : WorkdirLoc+'ZZ/',
    #'outputFileName' : OutDir+'nTupleZZ_TauTau.root',
    #'sample'         : 'ZZ_TauTau',
    #'xSection'       : 10.32,
    #'skimEff'        : 1.0,
    #'iJson'          : -1,
    #'iDiv'           : 0,
    #'nDiv'           : 1
    #},
#    ##Higgs MC
#    'GGFH125' : {
#    'inputFilePath'  : WorkdirLoc+'HiggsSM/GluGluToHToTauTau_M-125_MC_TauTau_v2',
#    'outputFileName' : OutDir+'nTupleGGFH125_TauTau.root',
#    'sample'         : 'GGFH125',
#    'xSection'       : 36.80 ,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,                                                      
#    'nDiv'           : 1
#    },
#    'VBFH125' : {
#        'inputFilePath'  : WorkdirLoc+'HiggsSM/VBFHToTauTau_M-125_MC_TauTau_v4',
#        #'inputFilePath'  : WorkdirLoc+'VBF_HToTauTau_M-125_13TeV-powheg-pythia6',
#        'outputFileName' : OutDir+'nTupleVBFH125_TauTau.root',
#        'sample'         : 'VBFH125',
#        'xSection'       : 36.80 ,
#        'skimEff'        : 1.0,
#        'iJson'          : -1,
#        'iDiv'           : 0,
#        'nDiv'           : 1
#        },
    ######MSSM
#    'SUSYGGH160' : {
#    #'inputFilePath'  : WorkdirLoc+'MC_Spring15_v1/SUSYGluGluToHToTauTau_M-160_PY8_25ns/',
#    'inputFilePath'  : WorkdirLoc+'SUSYGluGluToHToTauTau_M-160/',    
#    'outputFileName' : OutDir+'nTupleSUSYGGH160_TauTau.root',
#    'sample'         : 'SUSYGGH160',
#    'xSection'       : 1.0,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,
#    'nDiv'           : 1
#    },
#    'SUSYBBH160' : {
#    'inputFilePath'  : WorkdirLoc+'SUSYGluGluToBBHToTauTau_M-160/',
#    'outputFileName' : OutDir+'nTupleSUSYBBH160_TauTau.root',
#    'sample'         : 'SUSYGGH160',
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
    for analysis in [ 'nominal' ]: #, 'TauUp', 'TauDown' ]:
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
