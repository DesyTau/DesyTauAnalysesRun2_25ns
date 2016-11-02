for i in `ls Templates_*.root`
do

ln=`echo $i | cut -d "." -f1`

mv $i ${ln}_mt.root
done
