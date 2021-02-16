#!/bin/sh
dirData=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/data/
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/emb
dirMC=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/mc
dirMC2=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/mc_2
OUTDIR=./2018
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

ls $dirMC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu

ls $dirMC2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu
ls $dirMC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToSemiLeptonic
ls $dirMC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic

ls $dirMC2/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC2/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC2/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/VVTo2L2Nu
ls $dirMC2/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/WZTo2L2Q
ls $dirMC2/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/WZTo3LNu
ls $dirMC2/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/ZZTo2L2Q
ls $dirMC2/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/*root > $OUTDIR/ZZTo4L

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

ls $dirData/MuonEG_Run2018A/*.root > $OUTDIR/MuonEG_Run2018A
ls $dirData/MuonEG_Run2018B/*.root > $OUTDIR/MuonEG_Run2018B
ls $dirData/MuonEG_Run2018C/*.root > $OUTDIR/MuonEG_Run2018C
ls $dirData/MuonEG_Run2018D/*.root > $OUTDIR/MuonEG_Run2018D

ls $dirEmbedded/EmbeddingRun2018A_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018A
ls $dirEmbedded/EmbeddingRun2018B_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018B
ls $dirEmbedded/EmbeddingRun2018C_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018C
ls $dirEmbedded/EmbeddingRun2018D_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018D

for j in $(less list_SUSY_ggH_2018);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_pythia
done

for j in $(less list_SUSY_bbH_2018);
do
    ls $dirMC/${j}/*.root > $OUTDIR/${j}_amcatnlo
done
