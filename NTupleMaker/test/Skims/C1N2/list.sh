dir=/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/New8025/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test
systematics="Nominal _JetEnUp _JetEnDown _UnclEnUp _UnclEnDown _TauEnUp _TauEnDown _ElEnUp _ElEnDown _MuEnUp _MuEnDown _BTagUp _BTagDown"
systematics="Nominal"


dirm=${dir}/mutau
dire=${dir}/eltau
dirme=${dir}/muel

rm $dir/$1

for syst in $systematics
do

	
if [[ $syst != "Nominal" ]] ;then

dirm=${dir}/mutau${syst}
dire=${dir}/eltau${syst}
dirme=${dir}/muel${syst}
fi

echo $dirm
for i in `ls $1_[1-9]*.root`
do
	f=`echo $i | awk -F ".root" '{print $1}'`
	ls `pwd`/$i > $dirm/$f
	ls `pwd`/$i > $dire/$f
	ls `pwd`/$i > $dirme/$f


	#echo $i | awk -F ".root" '{print $1}' > $dir/$f

if [[ $syst == "Nominal" ]] ;then
	echo $f >> $dir/$1
fi

done
done
