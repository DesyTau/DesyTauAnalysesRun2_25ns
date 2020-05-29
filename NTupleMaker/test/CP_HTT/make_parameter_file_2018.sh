#!/bin/bash

CHANNEL=$1
echo "CONFIGFILE,FILELIST" > parameters.txt

if [[ $CHANNEL == "mt" ]]; then
  # data
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018A 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018B 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018C 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf SingleMuon_Run2018D 4

  # Embedded
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedMuTau_Run2018A 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedMuTau_Run2018B 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedMuTau_Run2018C 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedMuTau_Run2018D 4
else 
  if [[ $CHANNEL == "et" ]]; then
    # data
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018A 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf EGamma_Run2018D 4
    
    # Embedded
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElTau_Run2018A 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElTau_Run2018B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElTau_Run2018C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElTau_Run2018D 4
  else 
    echo "channel is not mt or et - exiting make_parameter_file_2018.sh"
    exit
  fi
fi

# Tau Spinner
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VBFHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WminusHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WplusHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZHToTauTauUncorrDecays_M125 4

# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DYJetsToLL_M-10to50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY4JetsToLL_M-50 10

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W4JetsToLNu 20

# EWK
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKWPlus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKWMinus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKZ2Jets_ZToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf EWKZ2Jets_ZToNuNu 10

# VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WW 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZ 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZ 20

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VVTo2L2Nu 10
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WWTo1L1Nu2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo2L2Q 10
#./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo4L_TuneCP5_13TeV_powheg_pythia8 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo4L 10

# H->WW
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VBFHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluZH_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf HWminusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf HWplusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf HZJ_HToWW_M125 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToSemiLeptonic 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_top_5f 20
