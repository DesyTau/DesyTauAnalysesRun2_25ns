
dir=$1

cd $1

for i in `ls *_19.pdf`
do

file=`echo $i | awk -F "_19" '{print $1}'`

mv $i ${file}_$2.pdf


done
cd ..
