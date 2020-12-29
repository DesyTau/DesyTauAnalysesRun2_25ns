#!/bin/bash
for j in A B C D
do
    ./hadd.sh DoubleMuon_Run2018${j}
done

samples=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
ST_t-channel_top
ST_t-channel_antitop
ST_tW_top
ST_tW_antitop
TTTo2L2Nu
TTToHadronic
TTToSemiLeptonic
WJetsToLNu
QCD_Pt-20to30_MuEnrichedPt5
QCD_Pt-30to50_MuEnrichedPt5
QCD_Pt-50to80_MuEnrichedPt5
QCD_Pt-80to120_MuEnrichedPt5
QCD_Pt-120to170_MuEnrichedPt5
QCD_Pt-170to300_MuEnrichedPt5
QCD_Pt-300to470_MuEnrichedPt5
QCD_Pt-470to600_MuEnrichedPt5
QCD_Pt-600to800_MuEnrichedPt5
QCD_Pt-800to1000_MuEnrichedPt5
QCD_Pt-1000toInf_MuEnrichedPt5
)

j=0
while [ $j -lt ${#samples[@]} ] 
do
    ./hadd.sh ${samples[$j]}
    j=`expr $j + 1` 
done

# change to {4..21} to run on all signal samples
for i in {8..9}
do
    ./hadd.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}
    ./hadd.sh SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done



   
