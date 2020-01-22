#!/bin/bash

### the script is to be run with "./gridcontrol_setup_mt_Run2.sh <year={16,17,18}> "

YEAR=$1
OUTDIR=./20$YEAR

./make_config_Run2.sh $YEAR MC
./make_config_Run2.sh $YEAR data
./make_config_Run2.sh $YEAR embedded

./make_lists_CP_20$YEAR.sh

cd $OUTDIR
./make_parameter_file.sh
cd ..

###################################################################################################
#
# After running this script setup the gridcontrol config file (./2017/gc_synch17.conf for example)
# Change project area, and storage area (se path)
#
###################################################################################################