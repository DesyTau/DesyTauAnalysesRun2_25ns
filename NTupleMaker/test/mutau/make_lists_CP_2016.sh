#!/bin/sh
dirDataSingleMuon=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/SingleMuon
dirMC=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc/
dirMC_v2=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v2/

ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY4JetsToLL_M-50

ls $dirMC/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*root > TT

ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > W4JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> W4JetsToLNu

ls $dirMC/ST_t-channel_antitop_4f/*root > ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f/*root > ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f/*root > ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f/*root > ST_tW_top_5f

ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8/*root > WW
ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> WW
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8/*root > WZ
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> WZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV/*root > ZZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV_ext1/*root >> ZZ


ls $dirMC_v2/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v1/*root > GluGluHToTauTau_M125
ls $dirMC_v2/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v3/*root >> GluGluHToTauTau_M125
ls $dirMC_v2/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125

ls $dirDataSingleMuon/SingleMuon_Run2016B-17Jul2018_ver1-v1/*root > SingleMuon_Run2016B
ls $dirDataSingleMuon/SingleMuon_Run2016B-17Jul2018_ver2-v1/*root >> SingleMuon_Run2016B
ls $dirDataSingleMuon/SingleMuon_Run2016C-17Jul2018-v1/*root > SingleMuon_Run2016C
ls $dirDataSingleMuon/SingleMuon_Run2016D-17Jul2018-v1/*root > SingleMuon_Run2016D
ls $dirDataSingleMuon/SingleMuon_Run2016E-17Jul2018-v1/*root > SingleMuon_Run2016E
ls $dirDataSingleMuon/SingleMuon_Run2016F-17Jul2018-v1/*root > SingleMuon_Run2016F
ls $dirDataSingleMuon/SingleMuon_Run2016G-17Jul2018-v1/*root > SingleMuon_Run2016G
ls $dirDataSingleMuon/SingleMuon_Run2016H-17Jul2018-v1/*root > SingleMuon_Run2016H
