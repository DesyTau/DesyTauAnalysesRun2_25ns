
#for i in `ls Templ*root`
for i in `ls Templ*Jet*_*17_*root`
do

file=`ls $i | cut -d "_" -f2-4`

unset name
name=`echo $file | cut -d "_" -f3`

if [[ $name -eq "0" ]] ; then
name=all_bins

else name=bin_${name}

fi

cp PlotsForUnrolledDistr_C PlotsForUnrolledDistr.C
sed -i 's/VARSHERE/'$file'/g' PlotsForUnrolledDistr.C
sed -i 's/BINHERE/'$name'/g' PlotsForUnrolledDistr.C
root -l -q -b PlotsForUnrolledDistr.C

done
