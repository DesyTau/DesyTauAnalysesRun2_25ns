#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:29:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=5000M
#
#(use hh site)
#$ -l site=hh
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL5
#$ -l os=sld6
#
# use current dir and current environment
#$ -cwd
#$ -V
#



cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;
flag=$2

channel=$3

	if [[  -z "$3" ]] ;then

		echo you must provide a channel....
		return 1
	fi

while read line
do


	
lt=`echo $line | cut -d '/' -f2`


	echo $lt > list_$lt
	
	#echo submitting  run_mc.sh $lt

	#if [[ ! -z "$2" ]] ;then
	if [[ $2 == *"W"* ]] ;then
		echo w template
		qsub run_plots_Wtemplate.sh list_$lt $3
		#qsub run_plots_WtemplateQCD.sh list_$lt
	fi
	
	if [[ $2 == *"MET"* ]] ;then
	echo met inverted
	 	qsub run_plotsB.sh list_$lt $3
	fi
	
	if [[ $2 == *"Inv"* ]] ;then
		echo inv region
	 	qsub run_plots_InvTemplate.sh list_$lt $3

	fi

	if [[ $2 == *"new"* ]] ;then
		echo  plots for new workflow 
	 	qsub -N p$3 -l h_rt=1:30:00 -l h_cpu=2000M run_plots_new.sh list_$lt $3
	 	#qsub -N pA$3 run_plots_A.sh list_$lt $3
	 	#qsub -N pB$3 run_plots_B.sh list_$lt $3
	 	#qsub -N pC$3 run_plots_C.sh list_$lt $3
	 	#qsub -N pD$3 run_plots_D.sh list_$lt $3

	fi

	if [[ $2 == *"Ttemplate"* ]] ;then
		echo inv region
	 	qsub run_plots_Ttemplate.sh list_$lt $3

	fi

done<$1

