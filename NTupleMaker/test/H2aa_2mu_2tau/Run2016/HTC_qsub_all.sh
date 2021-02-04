#!/bin/sh 
#Submit data samples:
for i in B C D E F G H
do
./HTC_qsub_seq.sh analysis_macro analysisMacro.conf SingleMuon_Run2016$i 20
done

#Submit signal MC samples:
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
do
./HTC_qsub_seq.sh analysis_macro analysisMacro_ggH.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i 3
done

#Submit background MC samples:
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf DYJetsToLL_M-10to50_13TeV-madgraphMLM 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf DYJetsToLL_M-50_13TeV-madgraphMLM 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf ST_tW_top_5f_inclusiveDecays_13TeV-powheg 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf WW_13TeV-pythia8 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf WZ_13TeV-pythia8 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf ZZ_13TeV-pythia8 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf TT_13TeV-powheg-pythia8 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf WJetsToLNu_13TeV-madgraphMLM 30
./HTC_qsub_seq.sh analysis_macro analysisMacro_mc.conf QCD_Pt-20toInf_MuEnrichedPt15_13TeV 30
