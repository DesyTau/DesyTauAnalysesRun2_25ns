#!/bin/bash

echo "CONFIGFILE,FILELIST" > parameters.txt

## Test for 2017
#data
./split_filelist.sh analysisMacroSynch_mt_17_data.conf SingleMuon_Run2017B 4
./split_filelist.sh analysisMacroSynch_mt_17_data.conf SingleMuon_Run2017C 4
./split_filelist.sh analysisMacroSynch_mt_17_data.conf SingleMuon_Run2017D 4
./split_filelist.sh analysisMacroSynch_mt_17_data.conf SingleMuon_Run2017E 4
./split_filelist.sh analysisMacroSynch_mt_17_data.conf SingleMuon_Run2017F 4

# Signals
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf GluGluHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf VBFHToTauTauUncorrDecays_M125 4

# DY
# TODO: add low mass DY (need to switch to the other PU file)
./split_filelist.sh analysisMacroSynch_mt_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DYJetsToLL_M-50 5
./split_filelist.sh analysisMacroSynch_mt_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY1JetsToLL_M-50 5
./split_filelist.sh analysisMacroSynch_mt_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY2JetsToLL_M-50 5
./split_filelist.sh analysisMacroSynch_mt_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY3JetsToLL_M-50 5
./split_filelist.sh analysisMacroSynch_mt_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY4JetsToLL_M-50 5

# Embedded
./split_filelist.sh analysisMacroSynch_mt_17_embedded.conf EmbeddedMuTau_Run2017B 3
./split_filelist.sh analysisMacroSynch_mt_17_embedded.conf EmbeddedMuTau_Run2017C 3
./split_filelist.sh analysisMacroSynch_mt_17_embedded.conf EmbeddedMuTau_Run2017D 3
./split_filelist.sh analysisMacroSynch_mt_17_embedded.conf EmbeddedMuTau_Run2017E 3
./split_filelist.sh analysisMacroSynch_mt_17_embedded.conf EmbeddedMuTau_Run2017F 3

# W+jets
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf WJetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf W1JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf W2JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W3JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W4JetsToLNu 10

# VV
./split_filelist.sh analysisMacroSynch_mt_WW_TuneCP5_13TeV-pythia8.conf WW 10
./split_filelist.sh analysisMacroSynch_mt_WZ_TuneCP5_13TeV-pythia8.conf WZ 10
./split_filelist.sh analysisMacroSynch_mt_ZZ_TuneCP5_13TeV-pythia8.conf ZZ 10

# TT
./split_filelist.sh analysisMacroSynch_mt_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTTo2L2Nu 3
./split_filelist.sh analysisMacroSynch_mt_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToHadronic 3
./split_filelist.sh analysisMacroSynch_mt_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToSemiLeptonic 3

# Single Top
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf ST_t-channel_antitop_4f 10
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf ST_t-channel_top_4f 10
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf ST_tW_antitop_5f 10
./split_filelist.sh analysisMacroSynch_mt_17_MC.conf ST_tW_top_5f 10
