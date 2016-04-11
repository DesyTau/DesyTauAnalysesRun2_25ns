while read line
do
	while read run
	do
count=`grep -c "$run" $line`
if [ $count -ne "0" ] ; then

echo $line has a good run $run
echo $run  >> found_runs
fi

done<$1
done<$2
