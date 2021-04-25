#!/bin/bash

YEAR=$1
SAMPLE=$2
OUTDIR=20$YEAR
GENERATOR=$3
cp list_${SAMPLE}_20${YEAR} ./${OUTDIR}
cp hadd.sh ./${OUTDIR}
cd ./${OUTDIR}

# Data
for j in $(less list_${SAMPLE}_20${YEAR});
do
    ./hadd.sh ${j}_${GENERATOR}
done

