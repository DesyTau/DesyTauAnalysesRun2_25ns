#!/bin/bash

CHANNEL=em
echo "CONFIGFILE,FILELIST" > parameters.txt

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018A 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018D 10

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018A 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018D 10

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018A 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018D 10

# Embedded
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018A 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018B 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018C 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018D 10
    
# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DYJetsToLL_M-10to50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY4JetsToLL_M-50 10

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WJetsToLNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W1JetsToLNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W2JetsToLNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W3JetsToLNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W4JetsToLNu 10

# EWK
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKWPlus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKWMinus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKZ2Jets_ZToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKZ2Jets_ZToNuNu 10

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WWTo1L1Nu2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToSemiLeptonic 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_antitop_4f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_top_4f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_antitop_5f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_top_5f 10

# SUSY_ggH
for j in $(less list_SUSY_ggH_2018);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf $j 5
done

# SUSY_bbH
for j in $(less list_SUSY_bbH_2018);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf $j 5
done
