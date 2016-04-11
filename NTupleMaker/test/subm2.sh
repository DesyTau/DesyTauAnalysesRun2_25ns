flag=$2

while read line
do

	
lt=`echo $line`
echo $line > $lt

	echo submitting  $1 $line
	#qsub run_mc2.sh $line
	qsub $1 $line
done<$2

