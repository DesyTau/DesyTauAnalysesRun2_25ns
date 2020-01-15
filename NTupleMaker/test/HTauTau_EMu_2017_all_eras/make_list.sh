#!/bin/sh
dirData=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data/MuonEG
dirDY=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirW=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirTT=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v4
dirTT2=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirST=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v4
dirVV=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirSignal=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirggh=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v4
dirEm=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/embedded/Embedding_emu_v2
dirWGammaStar=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3
dirEWKZ=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3

ls $dirDY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root > DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext
ls $dirDY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirDY/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root > DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext
ls $dirDY/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8  
ls $dirDY/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-10to50_13TeV-12Apr2018

ls $dirW/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirW/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8

ls $dirTT/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
ls $dirTT2/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
ls $dirTT2/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8
ls $dirTT/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > TTToHadronic_TuneCP5_13TeV-powheg-pythia8
ls $dirTT2/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8

ls $dirST/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
ls $dirST/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8
ls $dirST/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8

ls $dirVV/WW_TuneCP5_13TeV-pythia8/*root > WW_TuneCP5_13TeV-pythia8
ls $dirVV/WZ_TuneCP5_13TeV-pythia8/*root > WZ_TuneCP5_13TeV-pythia8
ls $dirVV/ZZ_TuneCP5_13TeV-pythia8/*root > ZZ_TuneCP5_13TeV-pythia8

ls $dirggh/GluGluHToTauTau_M125_13TeV_powheg_pythia8_new_pmx/*root > GluGluHToTauTau_M125_13TeV_powheg_pythia8
ls $dirggh/GluGluHToTauTau_M125_13TeV_powheg_pythia8_ext1/*root >> GluGluHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/*root > GluGluHToWWTo2L2Nu
ls $dirSignal/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/*root > VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8
ls $dirSignal/WminusHToTauTau_M125_13TeV_powheg_pythia8/*root > WminusHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/WplusHToTauTau_M125_13TeV_powheg_pythia8/*root > WplusHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/ZHToTauTau_M125_13TeV_powheg_pythia8/*root > ZHToTauTau_M125_13TeV_powheg_pythia8
ls $dirSignal/ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/*root > ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
ls $dirSignal/ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8/*root > ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8
ls $dirSignal/ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8/*root > ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8
ls $dirSignal/ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8/*root > ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8

ls $dirggh/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125/*root >  GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125
ls $dirggh/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125/*root >  GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125
ls $dirggh/VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125/*root >  VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125
ls $dirggh/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125/*root >  GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125
ls $dirggh/VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125/*root >  VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125
ls $dirggh/VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125/*root >  VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125
ls $dirggh/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125/*root >  GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125
ls $dirggh/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125/*root >  GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125

ls $dirData/MuonEG_Run2017B-31Mar2018-v1/*root > MuonEG_Run2017B
ls $dirData/MuonEG_Run2017C-31Mar2018-v1/*root > MuonEG_Run2017CtoF
ls $dirData/MuonEG_Run2017D-31Mar2018-v1/*root >> MuonEG_Run2017CtoF
ls $dirData/MuonEG_Run2017E-31Mar2018-v1/*root >> MuonEG_Run2017CtoF
ls $dirData/MuonEG_Run2017F-31Mar2018-v1/*root >> MuonEG_Run2017CtoF

ls $dirEm/EmbeddingRun2017B_ElMuFinalState-v2/*root > Embedding_Run2017
ls $dirEm/EmbeddingRun2017C_ElMuFinalState-v2/*root >> Embedding_Run2017
ls $dirEm/EmbeddingRun2017D_ElMuFinalState-v2/*root >> Embedding_Run2017
ls $dirEm/EmbeddingRun2017E_ElMuFinalState-v2/*root >> Embedding_Run2017
ls $dirEm/EmbeddingRun2017F_ElMuFinalState-v2/*root >> Embedding_Run2017

ls $dirWGammaStar/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/*root > WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8
ls $dirEWKZ/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8/*root > EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8
ls $dirEWKZ/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8/*root > EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8
ls $dirEWKZ/EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8/*root > EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8
ls $dirEWKZ/EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8/*root > EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8