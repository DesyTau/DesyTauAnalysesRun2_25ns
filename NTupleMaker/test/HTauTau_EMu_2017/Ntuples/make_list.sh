#!/bin/sh
dirData=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2/MuonEG
dirDY=/nfs/dust/cms/user/tlenz/13TeV/2017/NTuples/MC/12Apr2018_PU2017_METrecipe_v2
dirW=/nfs/dust/cms/user/tlenz/13TeV/2017/NTuples/MC/12Apr2018_PU2017_METrecipe_v2
dirTT=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2
dirST=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2
dirVV=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2
dirSignal=/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2/Signals
dirEm=/nfs/dust/cms/user/rasp/ntuples/Embedding_2017

ls $dirDY/DYJetsToLL_M-50_13TeV-12Apr2018/*root > DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DYJetsToLL_M-50_13TeV-12Apr2018_ext1/*root >> DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY1JetsToLL_M-50_13TeV-12Apr2018_ext1/*root > DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_13TeV-12Apr2018/*root > DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_13TeV-12Apr2018_ext1/*root >> DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY3JetsToLL_M-50_13TeV-12Apr2018/*root > DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY3JetsToLL_M-50_13TeV-12Apr2018_ext1/*root >> DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY4JetsToLL_M-50_13TeV-12Apr2018/*root > DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8  

ls $dirW/WJetsToLNu_TuneCP5_13TeV-12Apr2018/*root > WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/WJetsToLNu_TuneCP5_13TeV-12Apr2018-ext1/*root >> WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W1JetsToLNu_TuneCP5_13TeV-12Apr2018/*root > W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W2JetsToLNu_TuneCP5_13TeV-12Apr2018/*root > W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W3JetsToLNu_TuneCP5_13TeV-12Apr2018/*root > W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W4JetsToLNu_TuneCP5_13TeV-12Apr2018/*root > W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8

ls $dirTT/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8

ls $dirST/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ls $dirST/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8

ls $dirVV/WW_TuneCP5_13TeV-pythia8/*root > WW_TuneCP5_13TeV-pythia8
ls $dirVV/WZ_TuneCP5_13TeV-pythia8/*root > WZ_TuneCP5_13TeV-pythia8
ls $dirVV/ZZ_TuneCP5_13TeV-pythia8/*root > ZZ_TuneCP5_13TeV-pythia8

ls $dirSignal/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > GluGluHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125_13TeV_powheg_pythia8

ls $dirData/MuonEG_Run2017B/*root > MuonEG_Run2017
ls $dirData/MuonEG_Run2017C/*root >> MuonEG_Run2017
ls $dirData/MuonEG_Run2017D/*root >> MuonEG_Run2017
ls $dirData/MuonEG_Run2017E/*root >> MuonEG_Run2017
ls $dirData/MuonEG_Run2017F/*root >> MuonEG_Run2017

ls $dirEm/Embedding_Run2017B/*root > Embedding_Run2017
ls $dirEm/Embedding_Run2017C/*root >> Embedding_Run2017
ls $dirEm/Embedding_Run2017D/*root >> Embedding_Run2017
ls $dirEm/Embedding_Run2017E/*root >> Embedding_Run2017
ls $dirEm/Embedding_Run2017F/*root >> Embedding_Run2017
