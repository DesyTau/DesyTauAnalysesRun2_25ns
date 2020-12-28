#!/bin/sh
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/emb
dirSUSY=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Oct2020/2018/mc
OUTDIR=./2018
if [ ! -d "$OUTDIR" ]; then
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
fi

ls $dirEmbedded/EmbeddingRun2018D_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018D
ls $dirSUSY/SUSYGluGluToHToTauTau_M-1200/*.root > $OUTDIR/SUSYGluGluToHToTauTau_M-1200
