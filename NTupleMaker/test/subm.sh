flag=$2


systematics="Nominal"

#systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

#systematics="JetEnUp"

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test
while read line
do


for syst in $systematics
do
	
lt=`echo $line`
echo $line > $lt



	#echo $line > dt
	
	echo submitting  run_mc.sh $line for $2 channel and systematic = $syst
	qsub  -l h_vmem=2000M -l h_cpu=2:59:00 run_mc2.sh $line $2 $syst  

done
done<$1

