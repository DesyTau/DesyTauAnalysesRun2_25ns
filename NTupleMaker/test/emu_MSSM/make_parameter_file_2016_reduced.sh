#!/bin/bash

CHANNEL=em
echo "CONFIGFILE,FILELIST" > parameters.txt

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016D 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016E 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016F 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016G 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleMuon_Run2016H 10
  
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016D 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016E 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016F 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016G 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf SingleElectron_Run2016H 10

./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016D 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016E 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016F 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016G 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016H 10
    
# Embedded
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016D 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016E 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016F 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016G 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016H 10

# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DYJetsToLL_M-10to50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf DY4JetsToLL_M-50 10

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf W4JetsToLNu 20

# EWK
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf EWKWPlus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf EWKWMinus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf EWKZ2Jets_ZToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf EWKZ2Jets_ZToNuNu 10

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf TT 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_tW_top_5f 20

# SUSY_ggH
for j in $(less list_SUSY_ggH_2016);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf $j 5
done

# SUSY_bbH
for j in $(less list_SUSY_bbH_2016);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf $j 5
done
