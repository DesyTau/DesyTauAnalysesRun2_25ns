#!/bin/bash

for dir in $(find -maxdepth 1 -type d -name "*_files")
do
    filename=${dir:2:${#dir}-8}
    echo $filename
    if [ "$filename" == "SingleMuon_Run2017" ]; then
	     for i in {0..9}
	     do
	         hadd -f ${filename}_${i}.root ${filename}_files/*${i}.root
	     done
    else
	     hadd -f ${filename}.root ${filename}_files/*.root
    fi
done

hadd -f MuonEG.root MuonEG_Run2017CtoF.root MuonEG_Run2017B.root
hadd -f  DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root  DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root  DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext.root
hadd -f  DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root  DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root  DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext.root

hadd -f TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root
hadd -f TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root 
hadd -f TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root
