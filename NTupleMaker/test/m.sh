rm 25ns/merg.sh

#25ns/WJetsToLNu_MG_402_A_SS.root  25ns/WJetsToLNu_MG_402_B_OS.root  25ns/WJetsToLNu_MG_402_InvMET__C_OS.root  25ns/WJetsToLNu_MG_402_InvMET__D_SS.root
sel=InvMET
unset line
while read line 
do






	echo rm ${line}_${sel}_A.root >> 25ns/merg.sh
	echo hadd ${line}_${sel}_A.root ${line}_*_A_SS.root >> 25ns/merg.sh
	echo rm ${line}_*_A_SS.root >> 25ns/merg.sh
	echo "" >> 25ns/merg.sh
done<$1

