flag=$2

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test
while read line
do

	
lt=`echo $line`
echo $line > $lt




	#echo $line > dt
	
	echo submitting  run_mc.sh $line for $2 channel and point sparticle = $3 lsp = $4
	qsub run_mc2.sh $line $2 $3 1 1 1 1

done<$1

