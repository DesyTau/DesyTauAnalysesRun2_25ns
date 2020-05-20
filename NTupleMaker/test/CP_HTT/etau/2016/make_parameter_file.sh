#!/bin/bash

echo "CONFIGFILE,FILELIST" > parameters.txt

#data
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016B 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016C 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016D 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016E 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016F 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016G 10
./split_filelist.sh analysisMacroSynch_mt_16_data.conf SingleElectron_Run2016H 10

# Tau Spinner
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf GluGluHToTauTauUncorrDecays_M125 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf VBFHToTauTauUncorrDecays_M125 10
# ./split_filelist.sh analysisMacroSynch_mt_16_MC.conf WminusHToTauTauUncorrDecays_M125 10
# ./split_filelist.sh analysisMacroSynch_mt_16_MC.conf WplusHToTauTauUncorrDecays_M125 10
# ./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ZHToTauTauUncorrDecays_M125 10

# DY
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DYJetsToLL_M-10to50 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf DY4JetsToLL_M-50 10

# Embedded
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016B 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016C 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016D 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016E 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016F 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016G 10
./split_filelist.sh analysisMacroSynch_mt_16_embedded.conf EmbeddedElTau_Run2016H 10

# W+jets
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf W4JetsToLNu 20

# VV
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf WW 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf WZ 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ZZ 20

# TT
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf TT 10

# Single Top
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_mt_16_MC.conf ST_tW_top_5f 20
