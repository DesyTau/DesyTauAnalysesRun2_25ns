rm datasets25ns

for i in `ls /nfs/dust/cms/group/susy-desy//Run2/Stau/MC/25ns/cmssw7414v1_noMVAmet_v2/`
do
	ls /nfs/dust/cms/group/susy-desy//Run2/Stau/MC/25ns/cmssw7414v1_noMVAmet_v2/$i/*.root > 25ns/$i
	echo $i > $i
	echo $i >> datasets25ns


	echo hadd $i.root ${i}_*.root >> 25ns/merg.sh
	echo rm ${i}_* >> 25ns/merg.sh
done
