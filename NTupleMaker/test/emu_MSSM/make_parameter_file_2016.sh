#!/bin/bash

CHANNEL=em
echo "CONFIGFILE,FILELIST" > parameters.txt

# data

./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016B 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016C 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016D 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016E 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_data.conf MuonEG_Run2016F 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_dataGH.conf MuonEG_Run2016G 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_dataGH.conf MuonEG_Run2016H 150
    
# Embedded
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016B 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016C 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016D 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016E 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016F 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016G 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_embedded.conf EmbeddedElMu_Run2016H 4

# DY
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

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf TTToSemiLeptonic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf TTToHadronic 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ST_tW_top_5f 20

# H->WW
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf GluGluHToWWTo2L2Nu_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf VBFHToWWTo2L2Nu_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf HWminusJ_HToWW_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf HWplusJ_HToWW_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZHJ_HToWW_M125 4

# H->tautau
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf GluGluHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf VBFHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WplusHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf WminusHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ZHToTauTau_M125_13TeV 4


# SUSY_ggH
for j in $(less list_SUSY_ggH_2016);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ${j}_pythia 5
done

# SUSY_bbH
for j in $(less list_SUSY_bbH_2016);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_16_MC.conf ${j}_amcatnlo 5
done
