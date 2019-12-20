#!/bin/sh

echo "CONFIGFILE,FILELIST" > parameters.txt

#data
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2018A 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2018B 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2018C 5
../split_filelist.sh analysisMacroSynch_em_DATA.conf MuonEG_Run2018D 5

#DY Samples
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_DY.conf  DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8 5

#WJets
../split_filelist.sh   analysisMacroSynch_em_W.conf  WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_W.conf  W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_W.conf  W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_W.conf  W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_W.conf  W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 5


#SingleTop 
../split_filelist.sh   analysisMacroSynch_em_MC.conf  ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf  ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf  ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf  ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 5

#TTbar
../split_filelist.sh   analysisMacroSynch_em_MC.conf  TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf  TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf  TTToHadronic_TuneCP5_13TeV_powheg_pythia8 5

#Diboson 
../split_filelist.sh   analysisMacroSynch_em_MC.conf WW_TuneCP5_13TeV-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf WZ_TuneCP5_13TeV-pythia8 5
../split_filelist.sh   analysisMacroSynch_em_MC.conf ZZ_TuneCP5_13TeV-pythia8 5


#Embedded
../split_filelist.sh   analysisMacroSynch_em_Embedded.conf EmbeddingRun2018A 1
../split_filelist.sh   analysisMacroSynch_em_Embedded.conf EmbeddingRun2018B 1
../split_filelist.sh   analysisMacroSynch_em_Embedded.conf EmbeddingRun2018C 1
../split_filelist.sh   analysisMacroSynch_em_Embedded.conf EmbeddingRun2018D 1


#EWKZ
../split_filelist.sh  analysisMacroSynch_em_W.conf EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8 5
../split_filelist.sh  analysisMacroSynch_em_W.conf EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8 5
../split_filelist.sh  analysisMacroSynch_em_DY.conf EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8 5
../split_filelist.sh  analysisMacroSynch_em_DY.conf EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8 5

#WGamma
../split_filelist.sh  analysisMacroSynch_em_W.conf WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8 5

#Signals
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_M125_13TeV_powheg_pythia8  5
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal_vbf.conf VBFHToTauTau_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal_vbf.conf VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf WminusHToTauTau_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf WplusHToTauTau_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ZHToTauTau_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 5

../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125 5 
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_vbf.conf VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_ggh.conf GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_vbf.conf VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal_vbf.conf VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8 5
../split_filelist.sh  analysisMacroSynch_em_Signal.conf ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8 5