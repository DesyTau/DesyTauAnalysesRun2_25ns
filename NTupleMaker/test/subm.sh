flag=$2

while read line
do

	
lt=`echo $line`
echo $line > $lt


	if [ $flag == "eltau" ]

	then

		echo $line > dt
	
	echo submitting  run_mc.sh $line
	qsub run_mc_eltau.sh $line

	fi

	if [ $flag == "mutau" ]
	then

		echo $line > dt
	
	echo submitting  run_mc.sh $line
	qsub run_mc2.sh $line

fi
	if [ $flag == "signal" ]
	then

		echo $line > dt
	
	echo submitting  run_stau.sh $line
	qsub run_stau.sh $line

fi
done<$1

