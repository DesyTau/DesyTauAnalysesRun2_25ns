#!/bin/sh

CHANNEL=$1

dirMC=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/mc
if [[ $CHANNEL == "mt" ]]; then
  dirData=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/data_v2/SingleMuon
  dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau
  OUTDIR=./mutau/2016
else
  if [[ $CHANNEL == "et" ]]; then
    dirData=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/data_v2/SingleMuon
    dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/2016/embedded/Embedding_mutau
    OUTDIR=./etau/2016
  else 
    echo
    echo "To produce file lists for a specific channel this script is to be run with a command:"
    echo
    echo "  ./make_lists_CP_2016.sh <channel={mt,et}>"
    echo
    echo "channel is not mt or et - exiting"
    exit
  fi
fi

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-50
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
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/W4JetsToLNu

ls $dirMC/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8/*root > $OUTDIR/WW
ls $dirMC/WW_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> $OUTDIR/WW
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8/*root > $OUTDIR/WZ
ls $dirMC/WZ_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> $OUTDIR/WZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV-pythia8/*root > $OUTDIR/ZZ
ls $dirMC/ZZ_TuneCUETP8M1_13TeV-pythia8_ext1/*root >> $OUTDIR/ZZ

ls $dirMC/GluGluHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirMC/VBFHToTauTauUncorrelatedDecay_Filtered_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125

if [[ $CHANNEL == "mt" ]]; then
  ls $dirData/SingleMuon_Run2016B-17Jul2018_ver2/*root > $OUTDIR/SingleMuon_Run2016B
  ls $dirData/SingleMuon_Run2016C/*root > $OUTDIR/SingleMuon_Run2016C
  ls $dirData/SingleMuon_Run2016D/*root > $OUTDIR/SingleMuon_Run2016D
  ls $dirData/SingleMuon_Run2016E/*root > $OUTDIR/SingleMuon_Run2016E
  ls $dirData/SingleMuon_Run2016F/*root > $OUTDIR/SingleMuon_Run2016F
  ls $dirData/SingleMuon_Run2016G/*root > $OUTDIR/SingleMuon_Run2016G
  ls $dirData/SingleMuon_Run2016H/*root > $OUTDIR/SingleMuon_Run2016H
fi

if [[ $CHANNEL == "mt" ]]; then
  ls $dirEmbedded/EmbeddingRun2016B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016B
  ls $dirEmbedded/EmbeddingRun2016C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016C
  ls $dirEmbedded/EmbeddingRun2016D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016D
  ls $dirEmbedded/EmbeddingRun2016E_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016E
  ls $dirEmbedded/EmbeddingRun2016F_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016F
  ls $dirEmbedded/EmbeddingRun2016G_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016G
  ls $dirEmbedded/EmbeddingRun2016H_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016H
fi
