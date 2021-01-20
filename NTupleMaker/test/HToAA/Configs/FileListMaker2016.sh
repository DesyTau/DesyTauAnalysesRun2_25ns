##########################################################################################################################################################
##################******************************* List of NTuples for H->aa->4tau Analysis ************************************###########################
##########################################################################################################################################################

for i in 4 5 6 7 8 9 10 11 12 13 14 15 17 19 21
do

######## H->aa->4tau #########
#ggh
ls /nfs/dust/cms/user/perezdan/Store/NTuples/2016/H2aa_4tau/ggH/SUSYGluGluToHToAA_AToTauTau_M-${i}_13TeV_pythia8/*.root > SUSYGluGluToHToAA_AToTauTau_M-${i}

#vbf
ls /nfs/dust/cms/user/perezdan/Store/NTuples/2016/H2aa_4tau/VBF/SUSYVBFToHToAA_AToTauTau_M-${i}_13TeV_pythia8/*.root > SUSYVBFToHToAA_AToTauTau_M-${i}

#vh
ls /nfs/dust/cms/user/perezdan/Store/NTuples/2016/H2aa_4tau/VH/SUSYVH_HToAA_AToTauTau_M-${i}_13TeV_pythia8/*.root > SUSYVH_HToAA_AToTauTau_M-${i}

#tth
ls /nfs/dust/cms/user/perezdan/Store/NTuples/2016/H2aa_4tau/ttH/SUSYttH_HToAA_AToTauTau_M-${i}_13TeV_pythia8/*.root > SUSYttH_HToAA_AToTauTau_M-${i}

######## H->aa->2mu2tau #########
#ggH
ls /nfs/dust/cms/user/perezdan/Store/NTuples/2016/H2aa_2mu2tau/ggH/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}_13TeV_madgraph/*.root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}
 
done

for i in 16 18 20
do
#ggh
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2016/mc/SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}/*.root > SUSYGluGluToHToAA_AToTauTau_M-${i}

#vbf
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2016/mc/SUSYVBFToHToAA_AToTauTau_M-125_M-${i}/*.root > SUSYVBFToHToAA_AToTauTau_M-${i}

#vh
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2016/mc/SUSYVH_HToAA_AToTauTau_M-125_M-${i}/*.root > SUSYVH_HToAA_AToTauTau_M-${i}

#tth
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2016/mc/SUSYttH_HToAA_AToTauTau_M-125_M-${i}/*.root > SUSYttH_HToAA_AToTauTau_M-${i}  

done

###############################################################################################################################
##################******************************* MC Background ************************************###########################

ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/WZ_TuneCUETP8M1_13TeV-pythia8*/*.root > WZ_13TeV-pythia8
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/ZZ_TuneCUETP8M1_13TeV-pythia8*/*.root > ZZ_13TeV-pythia8
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/WW_TuneCUETP8M1_13TeV-pythia8*/*.root > WW_13TeV-pythia8
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*.root > WJetsToLNu_13TeV-madgraphMLM
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*.root > ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*.root > ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8*/*.root > ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8*/*.root > ST_tW_top_5f_inclusiveDecays_13TeV-powheg
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*.root > TT_13TeV-powheg-pythia8
ls /pnfs/desy.de/cms/tier2/store/user/alkaloge/2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ntuples_2016_metv2/190307_211223/0000/*.root > DYJetsToLL_M-10to50_13TeV-madgraphMLM
ls /pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/2016/2016-legacy/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*.root > DYJetsToLL_M-50_13TeV-madgraphMLM

### --- 15to30 30to50 50to80 80to120 120to170 170to300 300to470 470to600 600to800 800to1000 --- ###
ls /pnfs/desy.de/cms/tier2/store/user/dperezad/NTuple_Production/2016/QCD_MuEnriched/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/190702_123614/0000/*.root > QCD_Pt-20toInf_MuEnrichedPt15_13TeV



###############################################################################################################################                           
##################******************************* Data ************************************###########################
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016B-17Jul2018_ver2-v1/*.root > Run2016B-17Jul2018_ver2-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016C-17Jul2018-v1/*.root > Run2016C-17Jul2018-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016D-17Jul2018-v1/*.root > Run2016D-17Jul2018-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016E-17Jul2018-v1/*.root > Run2016E-17Jul2018-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016F-17Jul2018-v1/*.root > Run2016F-17Jul2018-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016G-17Jul2018-v1/*.root > Run2016G-17Jul2018-v1
ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data/DoubleMuon/Run2016H-17Jul2018-v1/*.root > Run2016H-17Jul2018-v1 


