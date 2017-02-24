

systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

systematics="Nominal"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown"
systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="MuEnUp MuEnDown"
systematics="Nominal JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

unset systematics

systematics="$3"


if [[ $3 == "list" ]] ; then

#systematics="UnclEnUp UnclEnDown"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
fi


cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test 
while read line
do

lt=`echo $line`
echo $line > $lt

for syst in $systematics
do
	



	#echo $line > dt
	
	echo submitting  run_mc.sh $line for $2 channel and systematic = $syst
	
			qsub  -l h_vmem=2000M -l h_cpu=2:59:00 run_mc2.sh $line $2 $syst  

done
done<$1
