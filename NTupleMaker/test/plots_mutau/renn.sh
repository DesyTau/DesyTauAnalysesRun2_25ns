#C1N2_600_LSP125_Nominal_B.root
while read line 
do

	p1=`echo $line | awk -F "_Nom" '{print $1}'`
	p2=`echo $line | awk -F "_Nom" '{print $1}'`

	mv $line ${p1}_B.root

done<$1
