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

channel=$2

	if [[  -z "$2" ]] ;then

		echo you must provide a channel....
		return 1
	fi

while read line
do


	
lt=`echo $line | cut -d '/' -f2`


	echo $lt > list_$lt
	
	

		echo  plots for channel $3 
	 	qsub -N p$2 -l h_rt=1:30:00 -l h_cpu=2000M run_plots_new.sh list_$lt $2 
	 	#qsub -N pA$3 run_plots_A.sh list_$lt $3
	 	#qsub -N pB$3 run_plots_B.sh list_$lt $3
	 	#qsub -N pC$3 run_plots_C.sh list_$lt $3
	 	#qsub -N pD$3 run_plots_D.sh list_$lt $3



done<$1

