
while read line
do
	
unset file
file=`echo $line | cut -d '/' -f2`

fileB=`echo $file | awk -F "_B.root" '{print $1}'`
fileA=`echo $file | awk -F "_A.root" '{print $1}'`

if [ ! -f plots/$file ] ; then
cp analyzer_h analyzer.h
cp analyzer_C analyzer.C

sed -i 's/FILEIN/'$file'/g' analyzer.h
sed -i 's/FILEIN/'$file'/g' analyzer.C

rm plots.root
root -l -q -b runme.C 
mv plots.root plots/$file
fi

done<$1
