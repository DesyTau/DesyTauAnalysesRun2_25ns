for i in `ls /nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_PHYS14_V4/`
do
	ls /nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_PHYS14_V4/$i/*.root > $i
done
