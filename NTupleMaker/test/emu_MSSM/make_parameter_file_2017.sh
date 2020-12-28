#!/bin/bash

CHANNEL=em
echo "CONFIGFILE,FILELIST" > parameters.txt

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017B 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017C 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017D 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017E 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017F 100

./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017B 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017C 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017D 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017E 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017F 100

./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017B 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017C 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017D 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017E 100
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017F 100
    
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017B 2
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017C 2
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017D 2
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017E 2
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017F 2

# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY3JetsToLL_M-50_ext1 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY4JetsToLL_M-50 10

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W4JetsToLNu 20

# EWK
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKWPlus2Jets_WToLNu_M-50 10
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKWMinus2Jets_WToLNu_M-50 10
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKZ2Jets_ZToLL_M-50 10
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKZ2Jets_ZToNuNu 10

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTTo2L2Nu_PSweights 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToSemiLeptonic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTToSemiLeptonic_PSweights 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_top_5f 20

# SUSY
# SUSY_ggH
#for j in $(less list_SUSY_ggH_2017);
#do
#    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_${j}.conf $j 5
#done

# SUSY_bbH
#for j in $(less list_SUSY_bbH_2017);
#do
#    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf $j 5
#done
