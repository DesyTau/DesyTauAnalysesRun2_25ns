#!/bin/sh

CHANNEL=$1

dirMC=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2016/mc
dirData=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2016/data
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/ywen/ntuples_Apr2020/2016/embedded

if [[ $CHANNEL == "mt" ]]; then
    OUTDIR=./mutau/2016
else   
    if [[ $CHANNEL == "et" ]]; then
	OUTDIR=./etau/2016
    else
	echo "ERROR: please run the script with ./make_lists_CP_2016.sh <channel={mt,et}>"
	exit
    fi
fi 

if [ ! -d "$OUTDIR" ]; then
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
fi

ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50
ls $dirMC/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-10to50

ls $dirMC/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic
ls $dirMC/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*root > $OUTDIR/TT

ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/*root > $OUTDIR/EWKWPlus2Jets_WToLNu_M-50
ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext1/*root >> $OUTDIR/EWKWPlus2Jets_WToLNu_M-50
ls $dirMC/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext2/*root >> $OUTDIR/EWKWPlus2Jets_WToLNu_M-50
ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/*root > $OUTDIR/EWKWMinus2Jets_WToLNu_M-50
ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext1/*root >> $OUTDIR/EWKWMinus2Jets_WToLNu_M-50
ls $dirMC/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_ext2/*root >> $OUTDIR/EWKWMinus2Jets_WToLNu_M-50

ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8/*root > $OUTDIR/EWKZ2Jets_ZToLL_M-50
ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_ext1/*root >> $OUTDIR/EWKZ2Jets_ZToLL_M-50
ls $dirMC/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_ext2/*root >> $OUTDIR/EWKZ2Jets_ZToLL_M-50
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8/*root > $OUTDIR/EWKZ2Jets_ZToNuNu
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_ext1/*root >> $OUTDIR/EWKZ2Jets_ZToNuNu
ls $dirMC/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_ext2/*root >> $OUTDIR/EWKZ2Jets_ZToNuNu

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

ls $dirMC/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/VVTo2L2Nu
ls $dirMC/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/WZTo2L2Q
ls $dirMC/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/WZTo3LNu
ls $dirMC/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/ZZTo2L2Q
ls $dirMC/ZZTo4L_13TeV_powheg_pythia8/*root > $OUTDIR/ZZTo4L

ls $dirMC/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > $OUTDIR/GluGluHToWWTo2L2Nu_M125
ls $dirMC/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > $OUTDIR/VBFHToWWTo2L2Nu_M125
ls $dirMC/GluGluZH_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/GluGluZH_HToWW_M125
ls $dirMC/HWminusJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/HWminusJ_HToWW_M125
ls $dirMC/HWplusJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/HWplusJ_HToWW_M125
ls $dirMC/HZJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/HZJ_HToWW_M125

ls $dirMC/GluGluHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/GluGluHToTauTauUncorrDecays_M125
ls $dirMC/VBFHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/VBFHToTauTauUncorrDecays_M125
ls $dirMC/WminusHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/WminusHToTauTauUncorrDecays_M125
ls $dirMC/WplusHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/WplusHToTauTauUncorrDecays_M125
ls $dirMC/ZHToTauTauUncorrelatedDecay_M125/*root > $OUTDIR/ZHToTauTauUncorrDecays_M125

if [[ $CHANNEL == "mt" ]]; then
  ls $dirData/SingleMuon/SingleMuon_Run2016B_ver2/*root > $OUTDIR/SingleMuon_Run2016B
  ls $dirData/SingleMuon/SingleMuon_Run2016C/*root > $OUTDIR/SingleMuon_Run2016C
  ls $dirData/SingleMuon/SingleMuon_Run2016D/*root > $OUTDIR/SingleMuon_Run2016D
  ls $dirData/SingleMuon/SingleMuon_Run2016E/*root > $OUTDIR/SingleMuon_Run2016E
  ls $dirData/SingleMuon/SingleMuon_Run2016F/*root > $OUTDIR/SingleMuon_Run2016F
  ls $dirData/SingleMuon/SingleMuon_Run2016G/*root > $OUTDIR/SingleMuon_Run2016G
  ls $dirData/SingleMuon/SingleMuon_Run2016H/*root > $OUTDIR/SingleMuon_Run2016H

  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016B_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016B
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016C_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016C
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016D_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016D
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016E_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016E
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016F_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016F
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016G_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016G
  ls $dirEmbedded/Embedding_mutau/EmbeddingRun2016H_MuTau/*root > $OUTDIR/EmbeddedMuTau_Run2016H
fi

if [[ $CHANNEL == "et" ]]; then
  ls $dirData/SingleElectron_Run2016B-17Jul2018_ver2-v1/*root > $OUTDIR/SingleElectron_Run2016B
  ls $dirData/SingleElectron_Run2016C-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016C
  ls $dirData/SingleElectron_Run2016D-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016D
  ls $dirData/SingleElectron_Run2016E-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016E
  ls $dirData/SingleElectron_Run2016F-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016F
  ls $dirData/SingleElectron_Run2016G-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016G
  ls $dirData/SingleElectron_Run2016H-17Jul2018-v1/*root > $OUTDIR/SingleElectron_Run2016H

  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016B_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016B
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016C_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016C
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016D_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016D
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016E_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016E
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016F_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016F
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016G_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016G
  ls $dirEmbedded/Embedding_eltau/EmbeddingRun2016H_ElTau/*root > $OUTDIR/EmbeddedElTau_Run2016H
fi
