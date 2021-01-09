#!/bin/bash

YEAR=20${1}
CHANNEL=em

if [ ! -d "$YEAR" ]; then
  mkdir ${YEAR}
fi

cp run_synchntyples.sh ./${YEAR}
cp split_filelist.sh ./${YEAR}
cp gc_synch.conf ./${YEAR} 
cp make_parameter_file_${YEAR}.sh ./${YEAR}
cp list_SUSY_ggH_${YEAR} ./${YEAR}
cp list_SUSY_bbH_${YEAR} ./${YEAR}

./make_config_Run2.sh $1 MC 
./make_config_Run2.sh $1 data 
./make_config_Run2.sh $1 embedded 
./make_lists_${YEAR}.sh

cd ./${YEAR}
rm parameters.txt
./make_parameter_file_${YEAR}.sh
cd -

echo "-----------------------------------------"
echo "DONT FORGET TO EDIT gc_synch.conf file!!!"
echo "-----------------------------------------"
