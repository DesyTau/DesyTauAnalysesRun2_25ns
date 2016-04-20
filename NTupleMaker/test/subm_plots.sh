flag=$2

while read line
do

	
lt=`echo $line | cut -d '/' -f2`


	echo $lt > list_$lt
	
	#echo submitting  run_mc.sh $lt

	#if [[ ! -z "$2" ]] ;then
	if [[ $2 == *"W"* ]] ;then
		echo w template
		qsub run_plots_Wtemplate.sh list_$lt
		#qsub run_plots_WtemplateQCD.sh list_$lt
	fi
	
	if [[ $2 == *"MET"* ]] ;then
	echo met inverted
	 	qsub run_plotsB.sh list_$lt
	fi
	
	if [[ $2 == *"Inv"* ]] ;then
		echo inv region
	 	qsub run_plots_InvTemplate.sh list_$lt

	fi

	if [[ $2 == *"new"* ]] ;then
		echo inv region
	 	qsub run_plots_new.sh list_$lt

	fi


done<$1

