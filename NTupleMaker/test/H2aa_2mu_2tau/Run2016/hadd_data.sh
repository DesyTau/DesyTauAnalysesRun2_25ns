#!/bin/sh 
#Join root files corresponding to data samples:
for i in B C D E F G H
do
./hadd.sh SingleMuon_Run2016$i
done

rm SingleMuon_Run2016.root 

hadd SingleMuon_Run2016.root SingleMuon_Run2016B.root SingleMuon_Run2016C.root SingleMuon_Run2016D.root SingleMuon_Run2016E.root SingleMuon_Run2016F.root SingleMuon_Run2016G.root SingleMuon_Run2016H.root
