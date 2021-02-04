#!/bin/sh 
#Submit data samples:
for i in B C D E F G H
do
./HTC_qsub_seq.sh analysis_macro analysisMacro.conf SingleMuon_Run2016$i 20
done
