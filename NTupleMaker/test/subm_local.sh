
rm ex.sh
systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

#systematics="JetEnUp"
systematics="Nominal"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown"
systematics="TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#systematics="MuEnUp MuEnDown"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"

cd /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test
while read line
do

lt=`echo $line`
echo $line > $lt

for syst in $systematics
do
	
	
	. run_mc2.sh $line $2 $syst  >>ex.sh
		

done
done<$1
