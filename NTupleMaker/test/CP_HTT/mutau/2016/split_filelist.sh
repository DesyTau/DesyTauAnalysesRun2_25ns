#!/bin/sh
# $1 - config file
# $2 - filelist (before 3)
# $3 - files per job (before 4)

echo ""
echo "$2 : "
let "n = 0"
rm -rf $2_files
mkdir $2_files
Splitter $2 $3

# Now add all lines to paramters.txt by looping over all files in the directory
for filename in $PWD/$2_files/*; do
    [ -e "$filename" ] || continue
    echo $1,$filename >> parameters.txt
done