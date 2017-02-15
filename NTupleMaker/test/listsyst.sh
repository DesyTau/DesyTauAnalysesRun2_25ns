
dir=$1

rm datasets${dir}
#sources="/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/ /nfs/dust/cms/group/susy-desy/Run2/Stau/MC/25ns/76x_JECv2_MVAMET0p6/"
#sources="/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/"

#sources="/nfs/dust/cms/group/higgs-kit/80x_v1 /nfs/dust/cms/group/higgs-kit/80x_v2 /nfs/dust/cms/group/higgs-kit/80x_SUSY_v1 80x_wMETFilters_v1"
sources="/nfs/dust/cms/group/higgs-kit/80x_wMETFilters_v1 /nfs/dust/cms/user/alkaloge/80x_wMETFilters_v1/"

source2="/nfs/dust/cms/group/higgs-kit/8020_23SeptReReco_v3"
source3="/nfs/dust/cms/group/higgs-kit/8020_23SeptReReco_v3"
#source3="/nfs/dust/cms/group/susy-desy/Run2/Stau/Data/25ns/8020_23SeptReReco"

systematics="JetEnUp JetEnDown UnclEnUp UnclEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"
#syst _JetEnUp _JetEnDown _TauEnUp _TauEnDown _ElEnUp _ElEnDown _MuEnUp _MuEnDown
systematics="UnclEnUp UnclEnDown"
systematics="JetEnUp JetEnDown TauEnUp TauEnDown ElEnUp ElEnDown MuEnUp MuEnDown"


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


ls $source2/SingleMuon*Run2016*/*.root > ${dir}/SingleMuon
ls $source3/SingleElectron*Run201*/*.root > ${dir}/SingleElectron
ls $source3/MuonEG*Run2016*/*.root > ${dir}/MuonEG

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
