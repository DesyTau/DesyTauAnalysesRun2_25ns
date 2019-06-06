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

hadd DY_buf.root DYJetsToLL.root DYJetsToLL_ext1.root
mv DY_buf.root DYJetsToLL.root
hadd DY_buf.root  DY1JetsToLL.root DY1JetsToLL_ext1.root
mv DY_buf.root DY1JetsToLL.root
hadd DY_buf.root  DY2JetsToLL.root DY2JetsToLL_ext1.root
mv DY_buf.root DY2JetsToLL.root
hadd DY_buf.root  DY3JetsToLL.root DY3JetsToLL_ext1.root
mv DY_buf.root DY3JetsToLL.root
hadd DATA_buf.root DATA_MuB.root DATA_MuC.root DATA_MuD.root DATA_MuE.root
mv DATA_buf.root DATA_RunsBCDE.root
hadd DATA_buf.root DATA_Mu*.root
mv DATA_buf.root DATA_SingleMuon.root