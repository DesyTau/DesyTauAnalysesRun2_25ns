#!/bin/sh 
#let "n = 0"
#for i in `ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_2_3/src/Submit/miniAOD/stau_stau100_LSP0/`
#for i in `ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_2_3/src/Submit/miniAOD/Output/WJetsToLNu_mlv_400toInf/`
for i in `ls /nfs/dust/cms/user/alkaloge/ACD/NAFtools-RunOnProcessed/CMSSW_7_2_3/src/Submit/miniAOD/Output/stau500_LSP400/`
 do
 echo submitting job for file $i
 ./qsub_flat.sh $i 
# let "n = n + 1"
done
#echo Total $n jobs submitted
