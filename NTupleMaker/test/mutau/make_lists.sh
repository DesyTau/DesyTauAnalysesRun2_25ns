#!/bin/bash

MCdir_VJets=/pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/MC/RunIIFall17MiniAOD-94X_mc2017_realistic/
MCdir_DY=/nfs/dust/cms/user/rasp/ntuples/MC_2017_v3
MCdir_TTVV=/nfs/dust/cms/user/rasp/ntuples/MC_2017_v2
MCdir_Signal=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/MCFall17/Signals
MCdir_ST=/nfs/dust/cms/user/rasp/ntuples/MC_2017_v2

# Signal

ls $MCdir_Signal/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*.root > ggH_125
ls $MCdir_Signal/VBFHToTauTau_M125_13TeV_powheg_pythia8/*.root > VBF_125
ls $MCdir_Signal/SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8/*.root > SUSYGluGluHTauTau_120


#DY

ls $MCdir_DY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DYJetsToLL
ls $MCdir_DY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*.root > DYJetsToLL_ext1
ls $MCdir_DY/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DY1JetsToLL
ls $MCdir_DY/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*.root > DY1JetsToLL_ext1
ls $MCdir_DY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DY2JetsToLL
ls $MCdir_DY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*.root > DY2JetsToLL_ext1
ls $MCdir_DY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DY3JetsToLL
ls $MCdir_DY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*.root > DY3JetsToLL_ext1
ls $MCdir_DY/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DY4JetsToLL
ls $MCdir_DY/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > DYJetsToLL_M-10to50


#Wjets

ls $MCdir_VJets/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > WJetsToLNu
ls $MCdir_VJets/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > W1JetsToLNu
ls $MCdir_VJets/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > W2JetsToLNu
ls $MCdir_VJets/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > W3JetsToLNu
ls $MCdir_VJets/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*.root > W4JetsToLNu


#Top

ls $MCdir_TTVV/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*.root > TTToSemiLeptonic
ls $MCdir_TTVV/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/*.root > TTTo2L2Nu
ls $MCdir_TTVV/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/*.root > TTToHadronic

# ST

ls $MCdir_ST/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root > ST_tW_top
ls $MCdir_ST/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root > ST_tW_antitop
ls $MCdir_ST/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*.root > ST_t_top
ls $MCdir_ST/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*.root > ST_t_antitop


#VV

ls $MCdir_TTVV/WW_TuneCP5_13TeV-pythia8/*.root > WW
ls $MCdir_TTVV/WZ_TuneCP5_13TeV-pythia8/*.root > WZ
ls $MCdir_TTVV/ZZ_TuneCP5_13TeV-pythia8/*.root > ZZ


#DATA

Data17dirMu=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/Run2017_Nov17ReReco/SingleMuon/

ls $Data17dirMu/SingleMuon_Run2017B*/*.root > DATA_MuB 
ls $Data17dirMu/SingleMuon_Run2017C*/*.root > DATA_MuC
ls $Data17dirMu/SingleMuon_Run2017D*/*.root > DATA_MuD
ls $Data17dirMu/SingleMuon_Run2017E*/*.root > DATA_MuE
ls $Data17dirMu/SingleMuon_Run2017F*/*.root > DATA_MuF

cat DATA_MuB DATA_MuC DATA_MuD DATA_MuE DATA_MuF > DATA_SingleMuon

#cat DATA_El* > DATA_SingleElectron

