while read line
do
	run=`echo $line | awk -F ":" '{print $1}'`
	echo for $run
	grep "$run" jobSin*
done<temp
