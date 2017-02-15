
dir=$1


systematics="Nominal JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown"


alias ls='ls'


for syst in $systematics
do

unset dir
dir=${1}_${syst}

if [[ $syst == "Nominal" ]] ; then

dir=${1}

fi


cd $dir
cp /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/mutau/m .
cp /nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/StauAnalysis/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/mutau/d .
. m
#. d
cd .. 

done

