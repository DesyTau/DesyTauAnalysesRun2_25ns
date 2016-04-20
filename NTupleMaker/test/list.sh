rm datasets${dir}

dir=$1
alias ls='ls'

for i in `ls /nfs/dust/cms/group/susy-desy//Run2/Stau/MC/25ns/cmssw7414v1_noMVAmet_v2/`
do
	ls /nfs/dust/cms/group/susy-desy//Run2/Stau/MC/25ns/cmssw7414v1_noMVAmet_v2/$i/*.root > ${dir}/$i
	echo $i > $i
	echo $i >> datasets${dir}
	

	echo hadd $i.root ${i}_*.root >> ${dir}/merg.sh
	echo rm ${i}_* >> ${dir}/merg.sh
	echo "" >> ${dir}/merg.sh
done

ls /nfs/dust/cms/group/susy-desy//Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v3/SingleMuon*/*.root > ${dir}/SingleMuon
ls /nfs/dust/cms/group/susy-desy//Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v3/SingleEl*/*.root > ${dir}/SingleElectron
ls /nfs/dust/cms/group/susy-desy//Run2/Stau/Data/25ns/cmssw7414v1_noMVAmet_v3/Jet*/*.root > ${dir}/JetsHT
echo JetsHT > JetsHT
echo SingleMuon  > SingleMuon
echo SingleElectron > SingleElectron

rm GC*
rm ${dir}/GC*

