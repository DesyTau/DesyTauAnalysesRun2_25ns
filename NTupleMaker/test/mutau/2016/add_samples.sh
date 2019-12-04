#!/bin/bash

for dir in $(ls | grep _files) 
do  
    echo $dir
    cd $dir
    sample=${dir%_*}
    echo $sample
    export samplename="$sample.root"
    echo $samplename
    hadd $samplename *.root
    mv $samplename ../$samplename
    cd ..
done    
