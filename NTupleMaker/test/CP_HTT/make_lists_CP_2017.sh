#!/bin/sh
dirDataSingleMuon=/nfs/dust/cms/user/rasp/ntuples/2017/data/SingleMuon
dirMC=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2017/mc
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2017/embedded/Embedding_mutau
OUTDIR=./2017

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50 
ls $dirMC/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-10to50

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_ex1/*root >> $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu

ls $dirMC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic
ls $dirMC/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/*root >> $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root >> $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root >> $OUTDIR/TTToHadronic

ls $dirMC/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_top_5f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_antitop_5f

ls $dirMC/WW_TuneCP5_13TeV-pythia8/*root > $OUTDIR/WW
ls $dirMC/WZ_TuneCP5_13TeV-pythia8/*root > $OUTDIR/WZ
ls $dirMC/ZZ_TuneCP5_13TeV-pythia8/*root > $OUTDIR/ZZ


ls $dirDataSingleMuon/SingleMuon_Run2017B/*root > $OUTDIR/SingleMuon_Run2017B
ls $dirDataSingleMuon/SingleMuon_Run2017C/*root > $OUTDIR/SingleMuon_Run2017C
ls $dirDataSingleMuon/SingleMuon_Run2017D/*root > $OUTDIR/SingleMuon_Run2017D
ls $dirDataSingleMuon/SingleMuon_Run2017E/*root > $OUTDIR/SingleMuon_Run2017E
ls $dirDataSingleMuon/SingleMuon_Run2017F/*root > $OUTDIR/SingleMuon_Run2017F

ls $dirMC/GluGluHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirMC/VBFHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125

# ls /pnfs/desy.de/cms/tier2/store/user/mvandekl/2017/mc/GluGluToHToTauTauNoSpin_Unfiltered_Rev1/*.root > $OUTDIR/GluGluToHToTauTauNoSpinCorr_Unfiltered
# ls /pnfs/desy.de/cms/tier2/store/user/mvandekl/2017/mc/VBFHToTauTauNoSpin_Unfiltered_Rev1/*.root > $OUTDIR/VBFHToTauTauNoSpinCorr_Unfiltered

ls $dirEmbedded/EmbeddingRun2017B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2017B
ls $dirEmbedded/EmbeddingRun2017C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2017C
ls $dirEmbedded/EmbeddingRun2017D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2017D
ls $dirEmbedded/EmbeddingRun2017E_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2017E
ls $dirEmbedded/EmbeddingRun2017F_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2017F
