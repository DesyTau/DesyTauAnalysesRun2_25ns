
dir=MET_filters

while read line
do

d=`echo $line | awk -F "list_" '{print $2}'  | awk -F ".txt" '{print $1}'`

echo $d

mkdir $dir/$d

cp met/MET_nof $dir/$d/.

done<$1
