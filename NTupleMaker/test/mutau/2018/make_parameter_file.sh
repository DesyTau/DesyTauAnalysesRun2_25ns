#!/bin/bash

echo "CONFIGFILE,FILELIST" > parameters.txt

## Test for 2018
#data
./split_filelist.sh analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018A 4
./split_filelist.sh analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018B 4
./split_filelist.sh analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018C 4
./split_filelist.sh analysisMacroSynch_mt_18_data.conf SingleMuon_Run2018D 4

# Signals
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf GluGluHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf VBFHToTauTau_M125 4
##  commented since now switching to TauSpinner weights
# ./split_filelist.sh analysisMacroSynch_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf SUSYGluGluToHToTauTau 4

# Tau Spinner
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf GluGluHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf VBFHToTauTauUncorrDecays_M125 4

# DY
# TODO: add low mass DY (need to switch to the other PU file)
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY4JetsToLL_M-50 10

# Embedded
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedMuTau_Run2018A 3
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedMuTau_Run2018B 3
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedMuTau_Run2018C 3
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedMuTau_Run2018D 3

# W+jets
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WJetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W1JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W2JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W3JetsToLNu 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W4JetsToLNu 10

# VV
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WW 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WZ 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ZZ 10

# TT
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTTo2L2Nu 3
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTToHadronic 3
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTToSemiLeptonic 3

# Single Top
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_t-channel_antitop_4f 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_t-channel_top_4f 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_tW_antitop_5f 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_tW_top_5f 10
