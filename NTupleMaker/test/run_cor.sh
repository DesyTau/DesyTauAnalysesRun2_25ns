
dir=25ns/InvMET
dir=plots
while read line
do

	if [ -f $dir/${line}_B_np.root ] ;then
	RunCor  all_bkg ${line}_B_np $dir

fi
done<$1
