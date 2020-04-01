#!/bin/sh
# $1 - executable
# $2 - analysis macro
# $3 - filelist
# $4 - channel
# $5 - files per job

let "n = 0"
rm -rf $3_files
mkdir $3_files
cp HTC_submit.sh $3_files/
cp $2 $3_files/
Splitter $3 $5
cd $3_files/
for i in `ls -1v $3_*`
 do
 echo submitting job $n for file $i from list $3
  cat > $3_$n.sh <<EOF
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd ${CMSSW_BASE}/src
cmsenv
cd -
$1 $2 $i $4
EOF
  chmod u+x $3_$n.sh
  ./HTC_submit.sh $3_$n.sh $3_$n
  let "n = n + 1"
done
echo Total $n jobs submitted
cd ../