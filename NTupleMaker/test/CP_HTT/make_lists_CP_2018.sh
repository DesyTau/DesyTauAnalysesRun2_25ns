#!/bin/sh

CHANNEL=$1

dirMC=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/mc_v2
if [[ CHANNEL == "mt" ]]; then
  dirData=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/data_v2/SingleMuon
  dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau
  OUTDIR=./mutau/2018
else
  if [[ CHANNEL == "et" ]]; then
    dirData=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/data_v2/SingleMuon
    dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2018/embedded/Embedding_mutau
    OUTDIR=./etau/2018
  else 
    echo
    echo "To produce file lists for a specific channel this script is to be run with a command:"
    echo
    echo "  ./make_lists_CP_2018.sh <channel={mt,et}>"
    echo
    echo "channel is not mt or et - exiting"
    exit
  fi
fi

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

ls $dirMC/GluGluToHToTauUncorrelatedDecay_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirMC/VBFHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125

if [[ CHANNEL == "mt" ]]; then
  ls $dirData/SingleMuon_Run2018A/*root > $OUTDIR/SingleMuon_Run2018A
  ls $dirData/SingleMuon_Run2018B/*root > $OUTDIR/SingleMuon_Run2018B
  ls $dirData/SingleMuon_Run2018C/*root > $OUTDIR/SingleMuon_Run2018C
  ls $dirData/SingleMuon_Run2018D/*root > $OUTDIR/SingleMuon_Run2018D
fi

if [[ CHANNEL == "mt" ]]; then
  ls $dirEmbedded/EmbeddingRun2018A_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018A
  ls $dirEmbedded/EmbeddingRun2018B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018B
  ls $dirEmbedded/EmbeddingRun2018C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018C
  ls $dirEmbedded/EmbeddingRun2018D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018D
fi
