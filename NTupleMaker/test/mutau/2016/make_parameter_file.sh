#!/bin/bash

echo "CONFIGFILE,FILELIST" > parameters.txt

## Test for 2016
#data
#./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleMuon_Run2017B 2
#./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleMuon_Run2017C 2
#./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleMuon_Run2017D 2
#./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleMuon_Run2017E 2
#./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleMuon_Run2017F 2

# Signals
#./split_filelist.sh analysisMacroSynch_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8.conf GluGluHToTauTau_M125 2
#./split_filelist.sh analysisMacroSynch_mt_VBFHToTauTau_M125.conf VBFHToTauTau_M125 2
##  commented since now switching to TauSpinner weights
# #./split_filelist.sh analysisMacroSynch_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf SUSYGluGluToHToTauTau 4

# Tau Spinner
#./split_filelist.sh analysisMacroSynch_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8.conf GluGluHToTauTauUncorrDecays_M125 2
#./split_filelist.sh analysisMacroSynch_mt_VBFHToTauTau_M125.conf VBFHToTauTauUncorrDecays_M125 2

# DY
# TODO: add low mass DY (need to switch to the other PU file)
#./split_filelist.sh analysisMacroSynch_mt_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DYJetsToLL_M-50 3
#./split_filelist.sh analysisMacroSynch_mt_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY1JetsToLL_M-50 3
#./split_filelist.sh analysisMacroSynch_mt_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY2JetsToLL_M-50 3
#./split_filelist.sh analysisMacroSynch_mt_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY3JetsToLL_M-50 3
#./split_filelist.sh analysisMacroSynch_mt_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY4JetsToLL_M-50 3

# Embedded
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedMuTau_Run2017B 1
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedMuTau_Run2017C 1
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedMuTau_Run2017D 1
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedMuTau_Run2017E 1
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedMuTau_Run2017F 1

# W+jets
#./split_filelist.sh analysisMacroSynch_mt_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf WJetsToLNu 5
#./split_filelist.sh analysisMacroSynch_mt_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W1JetsToLNu 5
#./split_filelist.sh analysisMacroSynch_mt_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W2JetsToLNu 5
#./split_filelist.sh analysisMacroSynch_mt_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W3JetsToLNu 5
#./split_filelist.sh analysisMacroSynch_mt_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W4JetsToLNu 5

# VV
#./split_filelist.sh analysisMacroSynch_mt_WW_TuneCP5_13TeV-pythia8.conf WW 5
#./split_filelist.sh analysisMacroSynch_mt_WZ_TuneCP5_13TeV-pythia8.conf WZ 5
#./split_filelist.sh analysisMacroSynch_mt_ZZ_TuneCP5_13TeV-pythia8.conf ZZ 5

# TT
#./split_filelist.sh analysisMacroSynch_mt_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTTo2L2Nu 2
#./split_filelist.sh analysisMacroSynch_mt_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToHadronic 2
#./split_filelist.sh analysisMacroSynch_mt_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToSemiLeptonic 2

# Single Top
#./split_filelist.sh analysisMacroSynch_mt_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf ST_t-channel_antitop_4f 5
#./split_filelist.sh analysisMacroSynch_mt_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf ST_t-channel_top_4f 5
#./split_filelist.sh analysisMacroSynch_mt_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf ST_tW_antitop_5f 5
#./split_filelist.sh analysisMacroSynch_mt_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf ST_tW_top_5f 5
