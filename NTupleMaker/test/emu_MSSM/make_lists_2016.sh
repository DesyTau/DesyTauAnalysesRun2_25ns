#!/bin/sh

CHANNEL=em

dirMC=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2016/mc
dirData=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2016/data
dirSingleMuon=/pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/Dec2020/2016/data
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2016/emb

OUTDIR=./2016

if [ ! -d "$OUTDIR" ]; then
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
fi

ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50

ls $dirMC/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic

ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/*root > $OUTDIR/W4JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2/*root >> $OUTDIR/W4JetsToLNu

ls $dirMC/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/VVTo2L2Nu
ls $dirMC/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/WZTo2L2Q
ls $dirMC/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/WZTo3LNu
ls $dirMC/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/ZZTo2L2Q
ls $dirMC/ZZTo4L_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/ZZTo4L

ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v1/*root > $OUTDIR/GluGluHToTauTau_M125
ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8-v3/*root >> $OUTDIR/GluGluHToTauTau_M125
ls $dirMC/VBFHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/VBFHToTauTau_M125
ls $dirMC/WminusHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/WminusHToTauTau_M125
ls $dirMC/WplusHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/WplusHToTauTau_M125
ls $dirMC/ZHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/ZHToTauTau_M125_13TeV

ls $dirMC/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > $OUTDIR/GluGluHToWWTo2L2Nu_M125
ls $dirMC/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/*root > $OUTDIR/VBFHToWWTo2L2Nu_M125
ls $dirMC/HWminusJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/HWminusJ_HToWW_M125
ls $dirMC/HWplusJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/HWplusJ_HToWW_M125
ls $dirMC/HZJ_HToWW_M125_13TeV_powheg_pythia8/*root > $OUTDIR/ZHJ_HToWW_M125

ls $dirSingleMuon/SingleMuon_Run2016B_ver2/*root > $OUTDIR/SingleMuon_Run2016B
ls $dirSingleMuon/SingleMuon_Run2016C/*root > $OUTDIR/SingleMuon_Run2016C
ls $dirSingleMuon/SingleMuon_Run2016D/*root > $OUTDIR/SingleMuon_Run2016D
ls $dirSingleMuon/SingleMuon_Run2016E/*root > $OUTDIR/SingleMuon_Run2016E
ls $dirSingleMuon/SingleMuon_Run2016F/*root > $OUTDIR/SingleMuon_Run2016F
ls $dirSingleMuon/SingleMuon_Run2016G/*root > $OUTDIR/SingleMuon_Run2016G
ls $dirSingleMuon/SingleMuon_Run2016H/*root > $OUTDIR/SingleMuon_Run2016H

ls $dirData/SingleElectron_Run2016B_ver2/*root > $OUTDIR/SingleElectron_Run2016B
ls $dirData/SingleElectron_Run2016C/*root > $OUTDIR/SingleElectron_Run2016C
ls $dirData/SingleElectron_Run2016D/*root > $OUTDIR/SingleElectron_Run2016D
ls $dirData/SingleElectron_Run2016E/*root > $OUTDIR/SingleElectron_Run2016E
ls $dirData/SingleElectron_Run2016F/*root > $OUTDIR/SingleElectron_Run2016F
ls $dirData/SingleElectron_Run2016G/*root > $OUTDIR/SingleElectron_Run2016G
ls $dirData/SingleElectron_Run2016H/*root > $OUTDIR/SingleElectron_Run2016H

ls $dirData/MuonEG_Run2016B_ver2/*root > $OUTDIR/MuonEG_Run2016B
ls $dirData/MuonEG_Run2016C/*root > $OUTDIR/MuonEG_Run2016C
ls $dirData/MuonEG_Run2016D/*root > $OUTDIR/MuonEG_Run2016D
ls $dirData/MuonEG_Run2016E/*root > $OUTDIR/MuonEG_Run2016E
ls $dirData/MuonEG_Run2016F/*root > $OUTDIR/MuonEG_Run2016F
ls $dirData/MuonEG_Run2016G/*root > $OUTDIR/MuonEG_Run2016G
ls $dirData/MuonEG_Run2016H/*root > $OUTDIR/MuonEG_Run2016H

ls $dirEmbedded/EmbeddingRun2016B_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016B
ls $dirEmbedded/EmbeddingRun2016C_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016C
ls $dirEmbedded/EmbeddingRun2016D_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016D
ls $dirEmbedded/EmbeddingRun2016E_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016E
ls $dirEmbedded/EmbeddingRun2016F_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016F
ls $dirEmbedded/EmbeddingRun2016G_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016G
ls $dirEmbedded/EmbeddingRun2016H_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2016H

for j in $(less list_SUSY_ggH_2016);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_pythia
done

for j in $(less list_SUSY_bbH_2016);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_amcatnlo
done

