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

systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown"
systematics="Nominal JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown UnclEnUp UnclEnDown TopPtUp TopPtDown ZPtUp ZPtDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test;eval `scramv1 runtime -sh` ;

channel=$2

	if [[  -z "$2" ]] ;then

		echo you must provide a channel....
		return 1
	fi

while read line
do

if [[  -z "$3" || $3 == "Nominal" ]] ;then
systematics="Nominal"
fi

if [[  $3 == "list" ]] ;then
systematics="list"
fi


lt=`echo $line | cut -d '/' -f2`


	echo $lt > list_$lt
	
	for syst in $systematics
	do
	

		echo  plots for channel $2 and syst $syst and $lt 
	 	qsub -N p$2 -l h_rt=1:30:00 -l h_cpu=2000M run_plots_new.sh list_$lt $2 $syst


	done
done<$1

