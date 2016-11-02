
model=$1
ls $1*.root | awk -F "_B.root" '{print $1}' > model



#for i in `ls Templ*.root`
#for i in `ls Templates_MTtot_MT2lester*_${model}*.root`
#for i in `ls Templates_MTtot_MT2lester_*_${model}*.root`
for i in `ls Templates*_${model}*.root`
do

ln=`echo $i  | cut -d '.' -f1`

echo the file is $ln
#. cards.sh $ln stau
. cards.sh $ln model

done
