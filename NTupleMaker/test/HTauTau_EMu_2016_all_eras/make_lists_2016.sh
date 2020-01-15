#!/bin/sh
dirData=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/Run2016-17Jul2018/MuonEG
dirMC=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc/
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/embedded/Embedding_emu_v2/
dirMC_v2=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v2/
dirMC_v3=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v3/
dirMC_v4=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v4/

ls $dirData/MuonEG_Run2016B-17Jul2018_ver2-v1/*root > MuonEG_Run2016B
ls $dirData/MuonEG_Run2016C-17Jul2018-v1/*root > MuonEG_Run2016C
ls $dirData/MuonEG_Run2016D-17Jul2018-v1/*root > MuonEG_Run2016D
ls $dirData/MuonEG_Run2016E-17Jul2018-v2/*root > MuonEG_Run2016E
ls $dirData/MuonEG_Run2016F-17Jul2018-v1/*root > MuonEG_Run2016F
ls $dirData/MuonEG_Run2016G-17Jul2018-v1/*root > MuonEG_Run2016G
ls $dirData/MuonEG_Run2016H-17Jul2018-v1/*root > MuonEG_Run2016H

ls $dirMC/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DYJetsToLL_M-10to50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > DY4JetsToLL_M-50

ls $dirMC/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*root > TTbar

ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > W4JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> W4JetsToLNu

ls $dirMC/ST_t-channel_antitop_4f/*root > ST_t-channel_antitop
ls $dirMC/ST_t-channel_top_4f/*root > ST_t-channel_top
ls $dirMC/ST_tW_antitop_5f/*root > ST_tW_antitop
ls $dirMC/ST_tW_top_5f/*root > ST_tW_top

ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8/*root > WW
ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> WW
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8/*root > WZ
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> WZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV/*root > ZZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV_ext1/*root >> ZZ

ls $dirEmbedded/EmbeddingRun2016B_ElMuFinalState-v1/*root > Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016C_ElMuFinalState-v1/*root >> Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016D_ElMuFinalState-v1/*root >> Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016E_ElMuFinalState-v1/*root >> Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016F_ElMuFinalState-v1/*root >> Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016G_ElMuFinalState-v1/*root >> Embedding_Run2016
ls $dirEmbedded/EmbeddingRun2016H_ElMuFinalState-v1/*root >> Embedding_Run2016

ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/*root > EWKWMinus2Jet
ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext1/*root >> EWKWMinus2Jet
ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext2/*root >> EWKWMinus2Jet
ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/*root > EWKWPlus2Jets
ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext1/*root >> EWKWPlus2Jets
ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext2/*root >> EWKWPlus2Jets
ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8//*root > EWKZ2Jets_ZToLL
ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_ext1/*root >> EWKZ2Jets_ZToLL
ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_ext2/*root >> EWKZ2Jets_ZToLL
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8/*root > EWKZ2Jets_ZToNuNu
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_ext1/*root >> EWKZ2Jets_ZToNuNu
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_ext2/*root >> EWKZ2Jets_ZToNuNu

ls $dirMC/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1/*root > WGToLNuG
ls $dirMC/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext2/*root >> WGToLNuG
ls $dirMC/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext3/*root >> WGToLNuG
ls $dirMC/WGstarToLNuEE_012Jets_13TeV-madgraph/*root > WGstarToLNuEE
ls $dirMC/WGstarToLNuMuMu_012Jets_13TeV-madgraph/*root > WGstarToLNuMuMu

ls $dirMC_v2/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v3/*root > GluGluHToTauTau_M125
find $dirMC_v4/GluGluHToTauTau_M125_13TeV_powheg_pythia8_ext1/ -type f -name "*.root" >GluGluHToTauTau_M125_ext1
ls $dirMC_v2/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > GluGluHToWWTo2L2Nu_M125
ls $dirMC_v2/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > VBFHToWWTo2L2Nu_M125
ls $dirMC_v2/WminusHToTauTau_M125_13TeV_powheg_pythia8/*root > WminusHToTauTau_M125
ls $dirMC_v2/WplusHToTauTau_M125_13TeV_powheg_pythia8/*root > WplusHToTauTau_M125
ls $dirMC_v2/ZHToTauTau_M125_13TeV_powheg_pythia8/*root > ZHToTauTau_M125
ls $dirMC_v2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/*root > ttHJetToTT_M125
ls $dirMC_v4/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > VBFHToTauTau_M125
find $dirMC_v4/VBFHToTauTau_M125_13TeV_powheg_pythia8_ext1/ -type f -name "*.root" >VBFHToTauTau_M125_ext1

ls $dirMC_v3/ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8/*root > ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8
ls $dirMC_v3/ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8/*root >  ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8
ls $dirMC_v3/ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8/*root > ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8