
for i in `ls *_B.root.txt`

do

f=`echo $i | awk -F "_B.root.txt" '{print $1}'`

echo $i $f
mv $i $f.root.txt

done

for i in `ls *_B_OS.root.txt`

do

f=`echo $i | awk -F "_B_OS.root.txt" '{print $1}'`

echo $i $f
mv $i $f.root.txt

done
