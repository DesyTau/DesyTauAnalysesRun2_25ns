for i in `ls *_B_OS.root`
do
file=`echo $i | awk -F "_B_OS" '{print $1}'`

mv $i ${file}_B.root
done
