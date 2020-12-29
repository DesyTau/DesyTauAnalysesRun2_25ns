#!/bin/bash

#File lists for 2018 legacy data:
for i in A B C D
do
    echo creating file list for data sample DoubleMuon_Run2018${i}
    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/data/DoubleMuon_Run2018${i}/*root > DoubleMuon_Run2018${i}
done

#File lists for 2018 MC ntuples v2:
echo creating file lists for HToAA_AToMuMu_AToTauTau samples
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
do
    ls /nfs/dust/cms/user/consuegs/ntuples/NMSSM_2018_v2/'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-'$i''/*root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i
done

echo creating file lists for HToAA_To4Tau samples
for i in {4..21}
do
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}/*root > SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYVBFToHToAA_AToTauTau_M-125_M-${i}/*root > SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYVH_HToAA_AToTauTau_M-125_M-${i} > SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYttH_HToAA_AToTauTau_M-125_M-${i} > SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done


#File lists for background MC samples
samples=(DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
WW_TuneCP5_13TeV-pythia8
WZ_TuneCP5_13TeV-pythia8
ZZ_TuneCP5_13TeV-pythia8
ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
TTToHadronic_TuneCP5_13TeV-powheg-pythia8
TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
)

names=(DYJetsToLL_M-10to50
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
)

names_QCD=(QCD_Pt-20to30_MuEnrichedPt5  
QCD_Pt-30to50_MuEnrichedPt5
QCD_Pt-50to80_MuEnrichedPt5
QCD_Pt-80to120_MuEnrichedPt5
QCD_Pt-120to170_MuEnrichedPt5
QCD_Pt-120to170_MuEnrichedPt5-ext1-v2
QCD_Pt-170to300_MuEnrichedPt5
QCD_Pt-300to470_MuEnrichedPt5
QCD_Pt-300to470_MuEnrichedPt5-ext3-v1
QCD_Pt-470to600_MuEnrichedPt5
QCD_Pt-470to600_MuEnrichedPt5-ext1-v2
QCD_Pt-600to800_MuEnrichedPt5
QCD_Pt-800to1000_MuEnrichedPt5
QCD_Pt-1000toInf_MuEnrichedPt5
)

i=0
while [ $i -lt ${#samples[@]} ] 
do
    echo "Creating file list for sample" ${samples[$i]} 

    #ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc/${samples[$i]}*/*root > ${names[$i]}
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/mc_v2/${samples[$i]}/*root > ${names[$i]}
      
    i=`expr $i + 1` 
done

j=0
while [ $j -lt ${#names_QCD[@]} ] 
do
    echo "Creating file list for sample" ${names_QCD[$j]} 
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/${names_QCD[$j]}/*root > ${names_QCD[$j]}
      
    j=`expr $j + 1` 
done

ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/QCD_Pt-80to120_MuEnrichedPt5-ext1-v2/QCD_Pt-80to120_MuEnrichedPt5-ext1-v2/*root > QCD_Pt-80to120_MuEnrichedPt5-ext1-v2

cat QCD_Pt-80to120_MuEnrichedPt5-ext1-v2 >> QCD_Pt-80to120_MuEnrichedPt5
cat QCD_Pt-120to170_MuEnrichedPt5-ext1-v2 >> QCD_Pt-120to170_MuEnrichedPt5
cat QCD_Pt-300to470_MuEnrichedPt5-ext3-v1 >> QCD_Pt-300to470_MuEnrichedPt5
cat QCD_Pt-470to600_MuEnrichedPt5-ext1-v2 >> QCD_Pt-470to600_MuEnrichedPt5
rm QCD_Pt-80to120_MuEnrichedPt5-ext1-v2
rm QCD_Pt-120to170_MuEnrichedPt5-ext1-v2
rm QCD_Pt-300to470_MuEnrichedPt5-ext3-v1
rm QCD_Pt-470to600_MuEnrichedPt5-ext1-v2
