#!/bin/bash

### the script is to be run with "./gridcontrol_setup_mt_Run2.sh <year={16,17,18}> <channel={mt,et}>"
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh

YEAR=$1
CHANNEL=$2
if [[ $CHANNEL == "mt" ]]; then
    OUTDIR=./mutau/20$YEAR
else  
    if [[ $CHANNEL == "et" ]]; then
	OUTDIR=./etau/20$YEAR
    else
	echo "ERROR: please run the script with ./gridcontrol_setup_mt_Run2.sh <year={16,17,18}> <channel={mt,et}>"
	exit
    fi
fi

./make_config_Run2.sh $YEAR MC $CHANNEL
./make_config_Run2.sh $YEAR data $CHANNEL
./make_config_Run2.sh $YEAR embedded $CHANNEL

./make_lists_CP_20$YEAR.sh $CHANNEL

cp ./make_parameter_file_20${YEAR}.sh $OUTDIR
cd $OUTDIR
./make_parameter_file_20${YEAR}.sh $CHANNEL
cd -

###################################################################################################
#
# Change make_parameter_file.sh in the ./20$YEAR directory before running this script
# After running this script setup the gridcontrol config file (./2017/gc_synch17.conf for example)
# Change project area, and storage area (se path)
# Then run grid-control
#
###################################################################################################
