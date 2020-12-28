#!/bin/bash

YEAR=20${1}
CHANNEL=em

cp run_synchntyples.sh ./${YEAR}
cp split_filelist.sh ./${YEAR}
cp gc_synch${1}.conf ./${YEAR} 
cp make_parameter_file_${YEAR}.sh ./${YEAR}

./make_config_Run2.sh $1 MC 
./make_config_Run2.sh $1 data 
./make_config_Run2.sh $1 embedded 
./make_lists_2018.sh

cd ./${YEAR}
rm parameters.txt
./make_parameter_file_${YEAR}.sh
cd -

