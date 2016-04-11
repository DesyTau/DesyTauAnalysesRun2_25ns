dir=$1
rm ${dir}/merg.sh

#${dir}/WJetsToLNu_MG_402_A_SS.root  ${dir}/WJetsToLNu_MG_402_B_OS.root  ${dir}/WJetsToLNu_MG_402_InvMET__C_OS.root  ${dir}/WJetsToLNu_MG_402_InvMET__D_SS.root
sel=InvMET
unset line
while read line 
do






	echo rm ${line}_${sel}_A.root >> ${dir}/merg.sh
	echo hadd ${line}_${sel}_A.root ${line}_*_A_SS.root >> ${dir}/merg.sh
	echo rm ${line}_*_A_SS.root >> ${dir}/merg.sh
	echo "" >> ${dir}/merg.sh
	echo rm ${line}_${sel}.root >> ${dir}/merg.sh
	echo hadd ${line}_${sel}.root ${line}_*_B_OS.root >> ${dir}/merg.sh
	echo rm ${line}_*_B_OS.root >> ${dir}/merg.sh

	echo "" >> ${dir}/merg.sh
	echo rm ${line}_${sel}_C.root >> ${dir}/merg.sh
	echo hadd ${line}_${sel}_C.root ${line}_*_${sel}__C_OS.root >> ${dir}/merg.sh
	echo rm ${line}_*_${sel}__C_OS.root >> ${dir}/merg.sh

	echo "" >> ${dir}/merg.sh
	echo rm ${line}_${sel}_D.root >> ${dir}/merg.sh
	echo hadd ${line}_${sel}_D.root ${line}_*_${sel}__D_SS.root >> ${dir}/merg.sh
	echo rm ${line}_*_${sel}__D_SS.root >> ${dir}/merg.sh
	echo "" >> ${dir}/merg.sh
done<$1

