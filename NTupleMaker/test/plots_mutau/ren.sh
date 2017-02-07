ls *Nominal*root > n
while read line
do

	p1=`echo $line | awk -F 'Nominal_' '{print $1}'`
	p2=`echo $line | awk -F 'Nominal_' '{print $2}'`

	cp $line $p1$p2
done<n

