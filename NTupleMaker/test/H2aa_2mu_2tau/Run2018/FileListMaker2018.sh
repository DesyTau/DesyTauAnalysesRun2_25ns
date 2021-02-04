#File lists for 2018 legacy data:
ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/SingleMuon_v2/SingleMuon_Run2018A-17Sep2018-v2/*root > SingleMuon_Run2018A

for i in B C
do
ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/SingleMuon_v2/SingleMuon_Run2018$i-17Sep2018-v1/*root > SingleMuon_Run2018$i
done

ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/SingleMuon_v2/SingleMuon_Run2018D-22Jan2019-v2/*root > SingleMuon_Run2018D


#File lists for 2018 MC ntuples v2:
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
do
ls /nfs/dust/cms/user/consuegs/ntuples/NMSSM_2018_v2/'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-'$i''/*root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i
done

#File lists for background MC samples
samples=(DYJetsToLL_M-10to50
DYJetsToLL_M-50   
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

samples_QCD=(QCD_Pt-20to30
QCD_Pt-30to50
QCD_Pt-50to80
QCD_Pt-80to120
QCD_Pt-120to170
QCD_Pt-170to300
QCD_Pt-300to470
QCD_Pt-470to600
QCD_Pt-600to800
QCD_Pt-800to1000
QCD_Pt-1000toInf
)

names=(DYJetsToLL_M-10to50_13TeV-madgraphMLM
DYJetsToLL_M-50_13TeV-madgraphMLM
WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg
ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg
ST_tW_top_5f_inclusiveDecays_13TeV-powheg
ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg
TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
TTToHadronic_TuneCP5_13TeV-powheg-pythia8
TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
WJetsToLNu_13TeV-madgraphMLM
)

names_QCD=(QCD_Pt-20to30_MuEnrichedPt5_TuneCP5  
QCD_Pt-30to50_MuEnrichedPt5_TuneCP5
QCD_Pt-50to80_MuEnrichedPt5_TuneCP5
QCD_Pt-80to120_MuEnrichedPt5_TuneCP5
QCD_Pt-120to170_MuEnrichedPt5_TuneCP5
QCD_Pt-170to300_MuEnrichedPt5_TuneCP5
QCD_Pt-300to470_MuEnrichedPt5_TuneCP5
QCD_Pt-470to600_MuEnrichedPt5_TuneCP5
QCD_Pt-600to800_MuEnrichedPt5_TuneCP5
QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5
QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5
)

j=1
i=1  

while [ $i -le ${#samples[@]} ] 
do
    echo "Creating file list for sample" ${samples[$i]} 

    #ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc/${samples[$i]}*/*root > ${names[$i]}
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/mc_v2/${samples[$i]}*/*root > ${names[$i]}
      
    i=`expr $i + 1` 
done

while [ $j -le ${#samples_QCD[@]} ] 
do
    echo "Creating file list for sample" ${samples_QCD[$j]} 
    ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc/${samples_QCD[$j]}*/*root > ${names_QCD[$j]}
      
    j=`expr $j + 1` 
done
