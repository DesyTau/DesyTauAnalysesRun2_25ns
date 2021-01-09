#!/bin/sh 
era=$1
trigger=$2
dir=/nfs/dust/cms/user/rasp/Run/emu_MSSM/Jan1/

echo "base directory" $dir

folder=${dir}/datacards_${trigger}/${era}
echo "merging file in dir " ${folder}
cd ${folder}
rm htt_em_mssm.root
hadd htt_em_mssm.root *.root	
cd -
