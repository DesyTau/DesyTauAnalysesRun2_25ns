#!/bin/sh
#Submit signal MC samples:
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
do
./HTC_qsub_seq.sh analysis_macro analysisMacro_ggH.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i 3
done

