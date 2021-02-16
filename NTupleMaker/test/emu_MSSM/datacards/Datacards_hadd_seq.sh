#!/bin/sh 

era=$1
trigger=$2

echo "base directory" $dir
for era in 2016 2017 2018
do
    ./Datacards_hadd.bash ${era}
done
