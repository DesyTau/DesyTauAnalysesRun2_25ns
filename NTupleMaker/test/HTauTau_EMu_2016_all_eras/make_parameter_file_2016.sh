#!/bin/sh
echo "CONFIGFILE,FILELIST" > parameters.txt

../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2016B 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2016C 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2016D 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2016E 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2016F 5
../split_filelist.sh analysisMacroSynch_em_GH.conf MuonEG_Run2016G  5
../split_filelist.sh analysisMacroSynch_em_GH.conf MuonEG_Run2016H  5


#DY Samples
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DY1JetsToLL_M-50    5
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DY2JetsToLL_M-50    5
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DY3JetsToLL_M-50    5
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DY4JetsToLL_M-50   5
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DYJetsToLL_M-50   5
../split_filelist.sh  analysisMacroSynch_em_DY.conf  DYJetsToLL_M-10to50  5

#WJets
../split_filelist.sh  analysisMacroSynch_em_W.conf  W1JetsToLNu   5
../split_filelist.sh  analysisMacroSynch_em_W.conf  W2JetsToLNu   5
../split_filelist.sh  analysisMacroSynch_em_W.conf  W3JetsToLNu   5
../split_filelist.sh  analysisMacroSynch_em_W.conf  W4JetsToLNu   5
../split_filelist.sh  analysisMacroSynch_em_W.conf   WJetsToLNu   5


#SingleTop 
../split_filelist.sh  analysisMacroSynch_em_mc.conf  ST_t-channel_antitop  4
../split_filelist.sh  analysisMacroSynch_em_mc.conf  ST_t-channel_top  4
../split_filelist.sh  analysisMacroSynch_em_mc.conf  ST_tW_top   4
../split_filelist.sh  analysisMacroSynch_em_mc.conf  ST_tW_antitop   4

#TTbar
../split_filelist.sh  analysisMacroSynch_em_mc.conf  TTbar   4

#Diboson 
../split_filelist.sh  analysisMacroSynch_em_mc.conf  ZZ  1
../split_filelist.sh  analysisMacroSynch_em_mc.conf  WW  1
../split_filelist.sh  analysisMacroSynch_em_mc.conf  WZ   1


#Signal
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf  GluGluHToTauTau_M125  4
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf  GluGluHToTauTau_M125_ext1 40 
../split_filelist.sh analysisMacroSynch_em_Signal_VBF.conf VBFHToTauTau_M125 4
../split_filelist.sh analysisMacroSynch_em_Signal_VBF.conf VBFHToTauTau_M125_ext1 40

../split_filelist.sh  analysisMacroSynch_em_Signal.conf  WminusHToTauTau_M125 1
../split_filelist.sh  analysisMacroSynch_em_Signal.conf  WplusHToTauTau_M125 1
../split_filelist.sh  analysisMacroSynch_em_Signal.conf  ZHToTauTau_M125 1
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf  GluGluHToWWTo2L2Nu_M125 1
../split_filelist.sh  analysisMacroSynch_em_Signal_VBF.conf  VBFHToWWTo2L2Nu_M125 1
../split_filelist.sh  analysisMacroSynch_em_Signal.conf  ttHJetToTT_M125 1

#WG and WGstar 
../split_filelist.sh  analysisMacroSynch_em_WG.conf  WGToLNuG   4
../split_filelist.sh  analysisMacroSynch_em_W.conf  WGstarToLNuEE   4
../split_filelist.sh  analysisMacroSynch_em_W.conf   WGstarToLNuMuMu  4

#EWK
../split_filelist.sh  analysisMacroSynch_em_W.conf EWKWMinus2Jet 5
../split_filelist.sh  analysisMacroSynch_em_W.conf EWKWPlus2Jets 5
../split_filelist.sh  analysisMacroSynch_em_DY.conf EWKZ2Jets_ZToLL  5
../split_filelist.sh  analysisMacroSynch_em_DY.conf EWKZ2Jets_ZToNuNu  5

../split_filelist.sh  analysisMacroSynch_em_Embedded.conf Embedding_Run2016  1

../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8 3
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8 3
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8 3