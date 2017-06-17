
dir=$1

rm datasets${dir}

#sources="/nfs/dust/cms/user/alkaloge/MoriondMC_v1_wWeights/ /nfs/dust/cms/user/rasp/ntuples/Moriond17_MC_v1/"
sources="/nfs/dust/cms/user/alkaloge/MoriondMC_v1_wWeights/  /nfs/dust/cms/user/rasp/ntuples/DYJets/"
source1="/nfs/dust/cms/group/higgs-kit/SingleElectron_Run2016-03Feb2017"
source2="/nfs/dust/cms/group/higgs-kit/MuonEG_Run2016-03Feb2017/"
source3="/nfs/dust/cms/user/fcost/store/Moriond17_03Feb2017/"

systematics="Nominal JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown BTagUp BTagDown"
#systematics="Nominal JetEnUp JetEnDown UnclEnUp UnclEnDown  MuEnUp MuEnDown BTagUp BTagDown"
#systematics="MuEnUp MuEnDown"
systematics="Nominal JetEnUp JetEnDown UnclEnUp UnclEnDown MuEnUp MuEnDown TauEnUp TauEnDown"


alias ls='ls'

#for i in `ls $source/`

for syst in $systematics
do

unset dir
dir=${1}_${syst}

	if  [[ $syst == "Nominal" ]] ; then

		dir=${1}
	fi

	if [[ ! -d ${dir}_${syst} ]] && [[ $syst != "Nominal" ]] ; then
		mkdir ${1}_${syst}
	fi


	echo Putting listing in $dir

	cp C1N2_names/* ${dir}/.
	for source in $sources
		do
			for i in `ls $source/`
				do


	ls $source/$i/*.root > ${dir}/$i
	echo $i > $i
	echo $i >> datasets${dir}
	
			done
	done

	#echo hadd $i.root ${i}_*.root >> ${dir}/merg.sh
	#echo rm ${i}_* >> ${dir}/merg.sh
	#echo "" >> ${dir}/merg.sh

ls $source1/SingleElectron*Run2016*/*.root > ${dir}/SingleElectron
ls $source3/SingleMuon*Run2016*/*.root > ${dir}/SingleMuon
ls $source2/MuonEG*/*.root > ${dir}/MuonEG



echo SingleMuon  > SingleMuon
echo SingleElectron > SingleElectron
echo MuonEG > MuonEG

rm GC*
rm ${dir}/GC*

rm *Glu*
rm $dir/*Glu*

rm AToZh*
rm $dir/AToZh*

rm Charged*
rm $dir/Charged*
done

sed -i '/AToZh/d' datasets$dir
sed -i '/Charged/d' datasets$dir
sed -i '/Glu/d' datasets$dir
sed -i '/GC/d' datasets$dir
