#!/bin/bash                                                                                          
# $1 - executable
# $2 - analysis macro
# $3 - filelist
# $4 - channel
let "n = 0"
rm -rf $3_files
mkdir $3_files
cp qsub_leptau.sh $3_files/
cp $2 $3_files/
Splitter $3 5
cd $3_files/
for i in $(ls $3_*)
 do
 echo submitting job $n for file $i from list $3
 ./qsub_leptau.sh $1 $2 $i $4
 let "n = n + 1"
done
echo Total $n jobs submitted
cd ../
