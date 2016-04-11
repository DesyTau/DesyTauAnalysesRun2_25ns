

dir=MET_filters
prefix=zpeak
prefix=dijet
while read line
do

d=`echo $line | awk -F "list_" '{print $2}'  | awk -F ".txt" '{print $1}'`

echo $d

cd $dir/$d

hadd -f ${prefix}_${d}.root MET_*.root
#rm MET_*
cp ${prefix}_${d}.root ../RootFiles/.

cd -

done<$1
