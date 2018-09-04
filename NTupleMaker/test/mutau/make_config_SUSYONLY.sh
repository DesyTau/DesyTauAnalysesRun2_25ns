#!/bin/bash


cp analysisMacroSynch_lept_mt_MC17.conf analysisMacroSynch_lept_mt_backup.conf



for file in $(ls *.conf)
do
    echo "**********************************"
    echo $file
    echo "**********************************"

    diff analysisMacroSynch_lept_mt_MC17.conf $file
done
