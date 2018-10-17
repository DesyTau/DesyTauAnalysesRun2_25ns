#!/bin/bash

for dir in $(find -maxdepth 1 -type d -name "*_files")
do
    filename=${dir:2:${#dir}-8}
    echo $filename
    hadd -f ${filename}.root ${filename}_files/*.root
done

root -l -b -q createDNNinput_2016.C+

