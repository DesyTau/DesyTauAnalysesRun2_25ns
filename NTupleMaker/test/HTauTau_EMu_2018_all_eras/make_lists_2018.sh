#!/bin/sh
dirDataSingleMuon=/nfs/dust/cms/user/ywen/Storage/SingleMuon
dirMuonEG=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/2018Data_17Sep2018/
dirDataEGamma=/nfs/dust/cms/user/cardinia/gridjobs/2018/NTuples/DATA/Run2018-17Sep2018
dirDYandTT=/nfs/dust/cms/user/ywen/Storage/MC
dirWVVandST=/nfs/dust/cms/user/cardinia/gridjobs/2018/NTuples/MC/RunIIAutumn18

ls $dirMuonEG/MuonEG_Run2018A/*root > MuonEG_Run2018A
ls $dirMuonEG/MuonEG_Run2018B/*root > MuonEG_Run2018B
ls $dirMuonEG/MuonEG_Run2018C/*root > MuonEG_Run2018C
ls $dirMuonEG/MuonEG_Run2018D/*root > MuonEG_Run2018D

ls $dirDataSingleMuon/SingleMuon_Run2018A_17Sep2018v2/*root > SingleMuon_Run2018A
ls $dirDataSingleMuon/SingleMuon_Run2018B_17Sep2018v1/*root > SingleMuon_Run2018B
ls $dirDataSingleMuon/SingleMuon_Run2018C_17Sep2018v1/*root > SingleMuon_Run2018C
ls $dirDataSingleMuon/SingleMuon_Run2018D_PromptRecov2/*root > SingleMuon_Run2018D


ls $dirDataEGamma/EGamma_Run2018A-17Sep2018-v2/*root > EGamma_Run2018A
ls $dirDataEGamma/EGamma_Run2018B-17Sep2018-v1/*root > EGamma_Run2018B
ls $dirDataEGamma/EGamma_Run2018C-17Sep2018-v1/*root > EGamma_Run2018C
ls $dirDataEGamma/EGamma_Run2018D-PromptReco-v2/*root > EGamma_Run2018D
ls $dirDataEGamma/EGamma_Run2018D-PromptReco-v1/*root >> EGamma_Run2018D


ls $dirDYandTT/DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8/*root > DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8
ls $dirDYandTT/DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8/*root > DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8
ls $dirDYandTT/DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8/*root > DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8
ls $dirDYandTT/DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8/*root > DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8
ls $dirDYandTT/DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8/*root > DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8


ls $dirDYandTT/TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/*root > TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8
ls $dirDYandTT/TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8/*root > TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8
ls $dirDYandTT/TTToHadronic_TuneCP5_13TeV_powheg_pythia8/*root > TTToHadronic_TuneCP5_13TeV_powheg_pythia8


ls $dirWVVandST/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirWVVandST/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirWVVandST/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirWVVandST/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirWVVandST/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8

ls $dirWVVandST/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/*root > ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ls $dirWVVandST/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/*root > ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8
ls $dirWVVandST/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/*root > ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ls $dirWVVandST/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/*root > ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8


ls $dirWVVandST/WW_TuneCP5_13TeV-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > WW_TuneCP5_13TeV-pythia8
ls $dirWVVandST/WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/*root > WZ_TuneCP5_13TeV-pythia8
ls $dirWVVandST/ZZ_TuneCP5_13TeV-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/*root > ZZ_TuneCP5_13TeV-pythia8

