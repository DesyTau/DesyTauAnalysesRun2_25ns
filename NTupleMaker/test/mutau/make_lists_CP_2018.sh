#!/bin/sh
dirDataSingleMuon=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/data
dirTauSpinner=/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc
dirMC=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/mc
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau
OUTDIR=./2018

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50 
ls $dirMC/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-10to50

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu

ls $dirMC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic

ls $dirMC/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC/WW_TuneCP5_13TeV-pythia8/*root > $OUTDIR/WW
ls $dirMC/WZ_TuneCP5_13TeV-pythia8/*root > $OUTDIR/WZ
ls $dirMC/ZZ_TuneCP5_13TeV-pythia8/*root > $OUTDIR/ZZ

ls $dirTauSpinner/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/GluGluHToTauTau_M125
ls $dirTauSpinner/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/VBFHToTauTau_M125

ls $dirDataSingleMuon/SingleMuon_Run2018A_17Sep2018_v2/*root > $OUTDIR/SingleMuon_Run2018A
ls $dirDataSingleMuon/SingleMuon_Run2018B_17Sep2018_v1/*root > $OUTDIR/SingleMuon_Run2018B
ls $dirDataSingleMuon/SingleMuon_Run2018C_17Sep2018_v1/*root > $OUTDIR/SingleMuon_Run2018C
ls $dirDataSingleMuon/SingleMuon_Run2018D_22Jan2019_v2/*root > $OUTDIR/SingleMuon_Run2018D

ls $dirEmbedded/EmbeddingRun2018A_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018A
ls $dirEmbedded/EmbeddingRun2018B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018B
ls $dirEmbedded/EmbeddingRun2018C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018C
ls $dirEmbedded/EmbeddingRun2018D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018D

ls $dirTauSpinner/GluGluToHToTauUncorrelatedDecay_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirTauSpinner/VBFHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125
