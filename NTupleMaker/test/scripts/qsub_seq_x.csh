#!/bin/csh 
# $1 - code name
# $2 - analysis macro
# $3 - base directory
# /nfs/dust/cms/group/susy-desy/Run2/Stau/MC/25ns/cmssw7414v1_noMVAmet_v2
# /nfs/dust/cms/user/rasp/ntuples/Sync_2015_v1
# $4 - sample name
# $5 - starting index
# $6 - last index
#rm -rf $4_files
#mkdir $4_files
cp $2 $4_files/
cp qsub.sh $4_files/
cd $4_files
set n = $5
while ( $n <= $6 )
    ls $3/$4/$4_${n}.root > $4_$2_$n
    ./qsub.sh $1 $2 $4_$2_$n
    @ n++
end
cd ../
