
dir=$1

rm datasets${dir}
#sources="/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/ /nfs/dust/cms/group/susy-desy/Run2/Stau/MC/25ns/76x_JECv2_MVAMET0p6/"
sources="/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/"

source2=/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6_DatawFilters/

alias ls='ls'

#for i in `ls $source/`

for source in $sources
do
for i in `ls $source/`
do
	ls $source/$i/*.root > ${dir}/$i
	echo $i > $i
	echo $i >> datasets${dir}
	

	echo hadd $i.root ${i}_*.root >> ${dir}/merg.sh
	echo rm ${i}_* >> ${dir}/merg.sh
	echo "" >> ${dir}/merg.sh
done
done


ls $source2/SingleMuon*Run2015D*/*.root > ${dir}/SingleMuon
ls $source2/SingleEl*Run2015D*/*.root > ${dir}/SingleElectron
ls $source2/MuonEG*Run2015D*/*.root > ${dir}/MuonEG
#ls $source/METw*/*.root > ${dir}/MET
#echo MET > MET
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


sed -i '/AToZh/d' datasets$dir
sed -i '/Charged/d' datasets$dir
sed -i '/Glu/d' datasets$dir
sed -i '/GC/d' datasets$dir
