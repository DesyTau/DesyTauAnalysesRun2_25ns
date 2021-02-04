#File lists for 2016 legacy data:
ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/data_v2/SingleMuon/SingleMuon_Run2016B-17Jul2018_ver2/*root > SingleMuon_Run2016B

for i in C D E F G H
do
ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/data_v2/SingleMuon/SingleMuon_Run2016$i/*root > SingleMuon_Run2016$i
done

#File lists for 2016 MC ntuples v2:
for i in 3p6 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
do
ls /nfs/dust/cms/user/consuegs/ntuples/NMSSM_2016_v2/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-${i}/*root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i
done

ls /nfs/dust/cms/user/consuegs/ntuples/NMSSM_2016/SUSYGGH_HToAA_AToTauTau_M-4_madgraph/*root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-4

#File lists for background MC samples
samples=(DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
WW_TuneCUETP8M1_13TeV-pythia8
WZ_TuneCUETP8M1_13TeV-pythia8
ZZ_TuneCUETP8M1_13TeV
ST_t-channel_top_4f
ST_t-channel_antitop_4f
ST_tW_top_5f
ST_tW_antitop_5f
TT_TuneCUETP8M2T4_13TeV-powheg-pythia8
WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8)

names=(DYJetsToLL_M-10to50_13TeV-madgraphMLM
DYJetsToLL_M-50_13TeV-madgraphMLM
WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg
ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg
ST_tW_top_5f_inclusiveDecays_13TeV-powheg
ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg
TT_13TeV-powheg-pythia8
WJetsToLNu_13TeV-madgraphMLM)
i=1

while [ $i -le ${#samples[@]} ] 
do
    echo "Creating file list for sample" ${samples[$i]} 

    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/mc/${samples[$i]}*/*root > ${names[$i]}
       
    i=`expr $i + 1` 
done

samples_DY_LM=(DYJetsToLL_M-1to5_HT-70to100
DYJetsToLL_M-1to5_HT-100to200
DYJetsToLL_M-1to5_HT-200to400
DYJetsToLL_M-1to5_HT-400to600
DYJetsToLL_M-1to5_HT-600toInf
DYJetsToLL_M-5to50_HT-70to100
DYJetsToLL_M-5to50_HT-100to200_ext
DYJetsToLL_M-5to50_HT-200to400_ext
DYJetsToLL_M-5to50_HT-400to600_ext
DYJetsToLL_M-5to50_HT-600toInf_ext)

names_DY_LM=(DYJetsToLL_M-1to5_HT-70to100
DYJetsToLL_M-1to5_HT-100to200
DYJetsToLL_M-1to5_HT-200to400
DYJetsToLL_M-1to5_HT-400to600
DYJetsToLL_M-1to5_HT-600toInf
DYJetsToLL_M-5to50_HT-70to100
DYJetsToLL_M-5to50_HT-100to200
DYJetsToLL_M-5to50_HT-200to400
DYJetsToLL_M-5to50_HT-400to600
DYJetsToLL_M-5to50_HT-600toInf)

i=1
while [ $i -le ${#samples_DY_LM[@]} ]
do
    echo "Creating file list for sample" ${samples_DY_LM[$i]}

    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/NMSSM_2016_v3/${samples_DY_LM[$i]}*/*root > ${names_DY_LM[$i]}

    i=`expr $i + 1`
done


### --- 15to30 30to50 50to80 80to120 120to170 170to300 300to470 470to600 600to800 800to1000 --- ###
ls /pnfs/desy.de/cms/tier2/store/user/dperezad/NTuple_Production/2016/QCD_MuEnriched/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/190702_123614/0000/*.root > QCD_Pt-20toInf_MuEnrichedPt15_13TeV
