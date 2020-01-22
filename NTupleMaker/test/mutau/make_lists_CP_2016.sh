#!/bin/sh
dirDataSingleMuon=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/SingleMuon
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau
dirMC=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc
dirMC_v2=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v2
dirMC_v4=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc_v4
dirTauSpinner=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/mc
OUTDIR=./2016

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50
ls $dirMC/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-10to50

ls $dirMC/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*root > $OUTDIR/TT

ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > $OUTDIR/W4JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/W4JetsToLNu

ls $dirMC/ST_t-channel_antitop_4f/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8/*root > $OUTDIR/WW
ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> $OUTDIR/WW
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8/*root > $OUTDIR/WZ
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> $OUTDIR/WZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV/*root > $OUTDIR/ZZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV_ext1/*root >> $OUTDIR/ZZ


ls $dirMC_v2/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v1/*root > $OUTDIR/GluGluHToTauTau_M125
ls $dirMC_v2/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v3/*root >> $OUTDIR/GluGluHToTauTau_M125
ls $dirMC_v4/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/VBFHToTauTau_M125

ls $dirDataSingleMuon/SingleMuon_Run2016B-17Jul2018_ver1-v1/*root > $OUTDIR/SingleMuon_Run2016B
ls $dirDataSingleMuon/SingleMuon_Run2016B-17Jul2018_ver2-v1/*root >> $OUTDIR/SingleMuon_Run2016B
ls $dirDataSingleMuon/SingleMuon_Run2016C-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016C
ls $dirDataSingleMuon/SingleMuon_Run2016D-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016D
ls $dirDataSingleMuon/SingleMuon_Run2016E-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016E
ls $dirDataSingleMuon/SingleMuon_Run2016F-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016F
ls $dirDataSingleMuon/SingleMuon_Run2016G-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016G
ls $dirDataSingleMuon/SingleMuon_Run2016H-17Jul2018-v1/*root > $OUTDIR/SingleMuon_Run2016H

ls $dirEmbedded/EmbeddingRun2016B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016B
ls $dirEmbedded/EmbeddingRun2016C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016C
ls $dirEmbedded/EmbeddingRun2016D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016D
ls $dirEmbedded/EmbeddingRun2016E_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016E
ls $dirEmbedded/EmbeddingRun2016F_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016F
ls $dirEmbedded/EmbeddingRun2016G_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016G
ls $dirEmbedded/EmbeddingRun2016H_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016H

ls $dirTauSpinner/GluGluHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirTauSpinner/VBFHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125
