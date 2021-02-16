#!/bin/sh 
era=$1
dir=/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/

echo "base directory" $dir

folder=${dir}/datacards/${era}
echo "merging file in dir " ${folder}
cd ${folder}
rm htt_em_mssm.root
hadd htt_em_mssm.root *.root	
cd -
