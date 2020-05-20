#!/bin/sh

CHANNEL=$1

dirMC=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2018/mc
dirData=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2018/data
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2018/embedded

if [[ -z "$CMSSW_BASE" ]]; then
  echo "Nah, set up your CMSSW first!"
  exit
fi

if [[ $CHANNEL == "mt" ]]; then
  OUTDIR=$CMSSW_BASE/src/DesyTauAnalyses/NTupleMaker/test/CP_HTT/mutau/2018
else
  if [[ $CHANNEL == "et" ]]; then
    OUTDIR=$CMSSW_BASE/src/DesyTauAnalyses/NTupleMaker/test/CP_HTT/etau/2018
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
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
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
ls $dirMC/WminusHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/WminusHToTauTauUncorrDecays_M125
ls $dirMC/WplusHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/WplusHToTauTauUncorrDecays_M125
ls $dirMC/ZHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/ZHToTauTauUncorrDecays_M125

if [[ $CHANNEL == "mt" ]]; then
  ls $dirData/SingleMuon/SingleMuon_Run2018A/*root > $OUTDIR/SingleMuon_Run2018A
  ls $dirData/SingleMuon/SingleMuon_Run2018B/*root > $OUTDIR/SingleMuon_Run2018B
  ls $dirData/SingleMuon/SingleMuon_Run2018C/*root > $OUTDIR/SingleMuon_Run2018C
  ls $dirData/SingleMuon/SingleMuon_Run2018D/*root > $OUTDIR/SingleMuon_Run2018D

  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2018A_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018A
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2018B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018B
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2018C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018C
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2018D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2018D
fi

if [[ $CHANNEL == "et" ]]; then
  ls $dirData/EGamma/EGamma_Run2018A/*root > $OUTDIR/EGamma_Run2018A
  ls $dirData/EGamma/EGamma_Run2018B/*root > $OUTDIR/EGamma_Run2018B
  ls $dirData/EGamma/EGamma_Run2018C/*root > $OUTDIR/EGamma_Run2018C
  ls $dirData/EGamma/EGamma_Run2018D/*root > $OUTDIR/EGamma_Run2018D
  
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2018A_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2018A
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2018B_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2018B
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2018C_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2018C
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2018D_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2018D
fi
