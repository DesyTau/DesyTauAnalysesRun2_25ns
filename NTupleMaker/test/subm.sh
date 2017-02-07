flag=$2


rm ex
systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

#systematics="JetEnUp"
systematics="Nominal"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown"
systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="MuEnUp MuEnDown"
#systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test
while read line
do

lt=`echo $line`
echo $line > $lt

for syst in $systematics
do
	



	#echo $line > dt
	
	echo submitting  run_mc.sh $line for $2 channel and systematic = $syst
	if [[ $3 == "local" ]] ; then
	       	echo . run_mc2.sh $line $2 $syst  >> ex
	
		else 
			qsub  -l h_vmem=2000M -l h_cpu=2:59:00 run_mc2.sh $line $2 $syst  
		
	fi

done
done<$1
. ex
