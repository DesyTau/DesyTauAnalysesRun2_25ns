#!/bin/sh 
#Join root files corresponding to data samples:
for i in B C D E F G H
do
./hadd.sh SingleMuon_Run2016$i
done

rm SingleMuon_Run2016.root 

hadd SingleMuon_Run2016.root SingleMuon_Run2016B.root SingleMuon_Run2016C.root SingleMuon_Run2016D.root SingleMuon_Run2016E.root SingleMuon_Run2016F.root SingleMuon_Run2016G.root SingleMuon_Run2016H.root

#Join root files corresponding to signal MC samples:
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
do
./hadd.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i
done

#Join root files corresponding to background MC samples:
./hadd.sh DYJetsToLL_M-10to50_13TeV-madgraphMLM
./hadd.sh DYJetsToLL_M-50_13TeV-madgraphMLM
./hadd.sh WW_13TeV-pythia8
./hadd.sh WZ_13TeV-pythia8
./hadd.sh ZZ_13TeV-pythia8
./hadd.sh ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg
./hadd.sh ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg
./hadd.sh ST_tW_top_5f_inclusiveDecays_13TeV-powheg
./hadd.sh ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg
./hadd.sh TT_13TeV-powheg-pythia8
./hadd.sh WJetsToLNu_13TeV-madgraphMLM
./hadd.sh QCD_Pt-20toInf_MuEnrichedPt15_13TeV

