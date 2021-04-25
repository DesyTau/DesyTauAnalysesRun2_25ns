#!/bin/bash

CHANNEL=em
echo "CONFIGFILE,FILELIST" > parameters.txt

# data
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_dataB.conf MuonEG_Run2017B 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017C 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017D 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017E 150
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf MuonEG_Run2017F 150
    
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017B 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017C 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017D 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017E 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElMu_Run2017F 4

# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DYJetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY1JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY2JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY3JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY3JetsToLL_M-50_ext1 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY4JetsToLL_M-50 20

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W4JetsToLNu 20

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTToSemiLeptonic 10

# H->WW
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf GluGluHToWWTo2L2Nu_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VBFHToWWTo2L2Nu_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf HWminusJ_HToWW_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf HWplusJ_HToWW_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZHJ_HToWW_M125 4

# H->tautau
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf GluGluHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VBFHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WplusHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WminusHToTauTau_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZHToTauTau_M125_13TeV 4

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_top_5f 20

# SUSY
# SUSY_ggH
for j in $(less list_SUSY_ggH_2017_powheg);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_${j}.conf ${j} 2
done

# SUSY_bbH
for j in $(less list_SUSY_bbH_2017_powheg);
do
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ${j} 2
done
