
dir=$2

rm ${dir}/merg.sh
while read i
do

	echo hadd $i.root ${i}_*.root >> ${dir}/merg.sh
	echo rm ${i}_* >> ${dir}/merg.sh
	echo "" >> ${dir}/merg.sh
done<$1


