#!/bin/sh
dirDataSingleMuon=/nfs/dust/cms/user/cardinia/gridjobs/NTuples/2017/DATA/Run2017-31Mar2018
dirDY=/nfs/dust/cms/user/klundert/HiggsCPNTuples/2017/MC
dirTT=/nfs/dust/cms/user/klundert/HiggsCPNTuples/2017/MC
dirST=/nfs/dust/cms/user/klundert/HiggsCPNTuples/2017/MC
dirW=/nfs/dust/cms/user/cardinia/gridjobs/NTuples/2017/MC/RunIIFall17
dirVV=/nfs/dust/cms/user/cardinia/gridjobs/NTuples/2017/MC/RunIIFall17
dirCPsignals=/nfs/dust/cms/user/cardinia/gridjobs/NTuples/2017/MC/RunIIFall17

ls $dirDY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
#ls $dirDY/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root > DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8  
ls $dirDY/DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8

ls $dirW/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8

ls $dirTT/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8

ls $dirST/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ls $dirST/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8

ls $dirVV/WW_TuneCP5_13TeV-pythia8/*root > WW_TuneCP5_13TeV-pythia8
ls $dirVV/WZ_TuneCP5_13TeV-pythia8/*root > WZ_TuneCP5_13TeV-pythia8
ls $dirVV/ZZ_TuneCP5_13TeV-pythia8/*root > ZZ_TuneCP5_13TeV-pythia8

ls $dirCPsignals/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > GluGluHToTauTau_M125_13TeV_powheg_pythia8
ls $dirCPsignals/GluGluHToTauTau_M125_13TeV_powheg_pythia8_ext1/*root >> GluGluHToTauTau_M125_13TeV_powheg_pythia8
ls $dirCPsignals/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125_13TeV_powheg_pythia8
ls $dirCPsignals/SUSYGluGluToHToTauTau_M-130_TuneCP5_13TeV-pythia8/*root > SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8
ls $dirCPsignals/SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8
ls $dirCPsignals/SUSYGluGluToBBHToTauTau_M-130_TuneCP5_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8
ls $dirCPsignals/SUSYGluGluToBBHToTauTau_M-120_TuneCP5_13TeV-pythia8/*root >> SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8

ls $dirDataSingleMuon/SingleMuon_Run2017B-31Mar2018-v1/*root > SingleMuon_Run2017B
ls $dirDataSingleMuon/SingleMuon_Run2017C-31Mar2018-v1/*root > SingleMuon_Run2017C
ls $dirDataSingleMuon/SingleMuon_Run2017D-31Mar2018-v1/*root > SingleMuon_Run2017D
ls $dirDataSingleMuon/SingleMuon_Run2017E-31Mar2018-v1/*root > SingleMuon_Run2017E
ls $dirDataSingleMuon/SingleMuon_Run2017F-31Mar2018-v1/*root > SingleMuon_Run2017F
