#!/bin/sh
dirDataSingleMuon=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data/SingleMuon
dirMC=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirMC_v4=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v4
dirTauSpinner=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v3/

ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY2JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY3JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY4JetsToLL_M-50 
ls $dirMC/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-10to50

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W4JetsToLNu

ls $dirMC_v4/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu
ls $dirMC_v4/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic
ls $dirMC_v4/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic

ls $dirMC_v4/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_top_4f
ls $dirMC_v4/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_antitop_4f
ls $dirMC_v4/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f
ls $dirMC_v4/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f

ls $dirMC/WW_TuneCP5_13TeV-pythia8/*root > WW
ls $dirMC/WZ_TuneCP5_13TeV-pythia8/*root > WZ
ls $dirMC/ZZ_TuneCP5_13TeV-pythia8/*root > ZZ

# ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > GluGluHToTauTau_M125
# ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8_ext1/*root >> GluGluHToTauTau_M125
ls $dirMC_v4/GluGluHToTauTau_M125_13TeV_powheg_pythia8_ext1/*root > GluGluHToTauTau_M125
ls $dirMC_v4/GluGluHToTauTau_M125_13TeV_powheg_pythia8_new_pmx/*root >> GluGluHToTauTau_M125
ls $dirMC/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125

ls $dirMC_v4/SUSYGluGluToHToTauTau_M-130_TuneCUETP8M1_13TeV-pythia8/*root > SUSYGluGluToHToTauTau
ls $dirMC_v4/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau
ls $dirMC_v4/SUSYGluGluToBBHToTauTau_M-130_TuneCUETP8M1_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau
ls $dirMC_v4/SUSYGluGluToBBHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau

ls $dirDataSingleMuon/SingleMuon_Run2017B-31Mar2018-v1/*root > SingleMuon_Run2017B
ls $dirDataSingleMuon/SingleMuon_Run2017C-31Mar2018-v1/*root > SingleMuon_Run2017C
ls $dirDataSingleMuon/SingleMuon_Run2017D-31Mar2018-v1/*root > SingleMuon_Run2017D
ls $dirDataSingleMuon/SingleMuon_Run2017E-31Mar2018-v1/*root > SingleMuon_Run2017E
ls $dirDataSingleMuon/SingleMuon_Run2017F-31Mar2018-v1/*root > SingleMuon_Run2017F

ls $dirTauSpinner/GluGluToHToTauTauNoSpin_M125/*root > GluGluHToTauTauUncorrDecays_M125
ls $dirTauSpinner/VBFHToTauTauNoSpin_M125/*root > GluGluHToTauTauUncorrDecays_M125

ls /pnfs/desy.de/cms/tier2/store/user/mvandekl/2017/mc/GluGluToHToTauTauNoSpin_Unfiltered_Rev1/*.root > GluGluToHToTauTauNoSpinCorr_Unfiltered
ls /pnfs/desy.de/cms/tier2/store/user/mvandekl/2017/mc/VBFHToTauTauNoSpin_Unfiltered_Rev1/*.root > VBFHToTauTauNoSpinCorr_Unfiltered
