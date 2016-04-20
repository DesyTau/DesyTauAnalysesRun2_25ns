
while read line
do 
unset f

unset ff
ff=`echo $line | awk -F ".root" '{print $1}' | awk -F "/" '{print $2}'`

cat bss > sub.sh


echo SUSYmutau analysisMacroSUSY.conf $ff 25ns  >> sub.sh

qsub sub.sh

done<25ns/scan

