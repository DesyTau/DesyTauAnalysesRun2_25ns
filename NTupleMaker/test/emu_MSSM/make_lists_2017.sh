#!/bin/sh

dirMC=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/mc
dirData=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/data
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2017/emb

OUTDIR=./2017

if [ ! -d "$OUTDIR" ]; then
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
fi

ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> $OUTDIR/DYJetsToLL_M-50
ls $dirMC/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root >> $OUTDIR/DY2JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1/*root > $OUTDIR/DY3JetsToLL_M-50_ext1
ls $dirMC/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50 

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_ex1/*root >> $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu

ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic
ls $dirMC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v1/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic

ls $dirMC/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_top_5f
ls $dirMC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_antitop_5f

ls $dirMC/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/VVTo2L2Nu
ls $dirMC/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/WZTo2L2Q
ls $dirMC/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/WZTo3LNu
ls $dirMC/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/ZZTo2L2Q
ls $dirMC/ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/ZZTo4L

ls $dirMC/GluGluHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/GluGluHToTauTau_M125
ls $dirMC/VBFHToTauTau_M125_13TeV_powheg_pythia8/*.root > $OUTDIR/VBFHToTauTau_M125
ls $dirMC/WplusHToTauTau_M125_13TeV_powheg_pythia8/*.root > $OUTDIR/WplusHToTauTau_M125
ls $dirMC/WminusHToTauTau_M125_13TeV_powheg_pythia8/*.root > $OUTDIR/WminusHToTauTau_M125
ls $dirMC/ZHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/ZHToTauTau_M125_13TeV

ls $dirMC/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/*root > $OUTDIR/GluGluHToWWTo2L2Nu_M125
ls $dirMC/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/*root > $OUTDIR/VBFHToWWTo2L2Nu_M125
ls $dirMC/HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5/*root > $OUTDIR/HWminusJ_HToWW_M125
ls $dirMC/HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5/*root > $OUTDIR/HWplusJ_HToWW_M125
ls $dirMC/HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5/*root > $OUTDIR/ZHJ_HToWW_M125

ls $dirData/SingleMuon_Run2017B/*root > $OUTDIR/SingleMuon_Run2017B
ls $dirData/SingleMuon/Run2017C/SingleMuon/*root > $OUTDIR/SingleMuon_Run2017C
ls $dirData/SingleMuon/Run2017D/SingleMuon/*root > $OUTDIR/SingleMuon_Run2017D
ls $dirData/SingleMuon/Run2017E/SingleMuon/*root > $OUTDIR/SingleMuon_Run2017E
ls $dirData/SingleMuon/Run2017F/SingleMuon/*root > $OUTDIR/SingleMuon_Run2017F

ls $dirData/SingleElectron_Run2017B/*root > $OUTDIR/SingleElectron_Run2017B
ls $dirData/SingleElectron_Run2017C/*root > $OUTDIR/SingleElectron_Run2017C
ls $dirData/SingleElectron_Run2017D/*root > $OUTDIR/SingleElectron_Run2017D
ls $dirData/SingleElectron_Run2017E/*root > $OUTDIR/SingleElectron_Run2017E
ls $dirData/SingleElectron_Run2017F/*root > $OUTDIR/SingleElectron_Run2017F

ls $dirData/MuonEG_Run2017B/*.root > $OUTDIR/MuonEG_Run2017B
ls $dirData/MuonEG_Run2017C/*.root > $OUTDIR/MuonEG_Run2017C
ls $dirData/MuonEG_Run2017D/*.root > $OUTDIR/MuonEG_Run2017D
ls $dirData/MuonEG_Run2017E/*.root > $OUTDIR/MuonEG_Run2017E
ls $dirData/MuonEG_Run2017F/*.root > $OUTDIR/MuonEG_Run2017F

ls $dirEmbedded/EmbeddingRun2017B_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2017B
ls $dirEmbedded/EmbeddingRun2017C_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2017C
ls $dirEmbedded/EmbeddingRun2017D_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2017D
ls $dirEmbedded/EmbeddingRun2017E_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2017E
ls $dirEmbedded/EmbeddingRun2017F_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2017F

for j in $(less list_SUSY_ggH_2017);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_pythia
done

for j in $(less list_SUSY_bbH_2017);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_amcatnlo
done

