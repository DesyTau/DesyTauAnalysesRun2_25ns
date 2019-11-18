#!/bin/sh
dirDataSingleMuon=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/SingleMuon
dirMC=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc

ls $dirMC/DYJetsToLL_M-50/*root > DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50/*root > DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50/*root > DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50/*root > DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50/*root > DY4JetsToLL_M-50 
ls $dirMC/DYJetsToLL_M-10to50/*root > DYJetsToLL_M-10to50

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W4JetsToLNu

ls $dirMC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic

ls $dirMC/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f

ls $dirMC/WW_TuneCP5_13TeV-pythia8/*root > WW
ls $dirMC/WZ_TuneCP5_13TeV-pythia8/*root > WZ
ls $dirMC/ZZ_TuneCP5_13TeV-pythia8/*root > ZZ

ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > GluGluHToTauTau_M125
ls $dirMC/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125

ls $dirDataSingleMuon/SingleMuon_Run2018A-17Sep2018-v2/*root > SingleMuon_Run2018A
ls $dirDataSingleMuon/SingleMuon_Run2018B-17Sep2018-v1/*root > SingleMuon_Run2018B
ls $dirDataSingleMuon/SingleMuon_Run2018C-17Sep2018-v1/*root > SingleMuon_Run2018C
ls $dirDataSingleMuon/SingleMuon_Run2018D-22Jan2019-v2/*root > SingleMuon_Run2018D
