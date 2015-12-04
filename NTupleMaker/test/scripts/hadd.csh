#!/bin/csh
cd $1_files
rm $1.root 
hadd $1.root *.root
mv $1.root ../
cd ../
