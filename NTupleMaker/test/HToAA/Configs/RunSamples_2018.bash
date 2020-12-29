#!/bin/bash

for i in A B C D
do
    echo submitting jobs for sample DoubleMuon_Run2018${i}
    ./HTC_submit_seq.sh analysis_macro analysisMacro_2018.conf DoubleMuon_Run2018${i} 200
done

echo submitting jobs for HToAA_AToMuMu_AToTauTau samples
for i in {4..21}
do
    ./HTC_submit_seq.sh analysis_macro analysisMacro_ggH_2018.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i} 5
done

echo submitting jobs for HToAA_To4Tau samples
for i in {4..21}
do
    ./HTC_submit_seq.sh analysis_macro analysisMacro_ggH_2018.conf SUSYGluGluToHToAA_AToTauTau_M-125_M-${i} 5
    ./HTC_submit_seq.sh analysis_macro analysisMacro_VBF_2018.conf SUSYVBFToHToAA_AToTauTau_M-125_M-${i} 5
    ./HTC_submit_seq.sh analysis_macro analysisMacro_VH_2018.conf SUSYVH_HToAA_AToTauTau_M-125_M-${i} 5
    ./HTC_submit_seq.sh analysis_macro analysisMacro_ttH_2018.conf SUSYttH_HToAA_AToTauTau_M-125_M-${i} 5
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
    echo "Submitting jobs on sample " ${samples[$j]} 
    ./HTC_submit_seq.sh analysis_macro analysisMacro_VH_2018.conf ${samples[$j]} 100
    j=`expr $j + 1` 
done

