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

hadd DATA_SingleMuon.root SingleMuon_Run2017B.root SingleMuon_Run2017C.root SingleMuon_Run2017D.root SingleMuon_Run2017E.root SingleMuon_Run2017F.root

