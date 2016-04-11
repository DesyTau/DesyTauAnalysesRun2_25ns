while read line
do
lt=`echo $line`
echo $line > $lt
done<$1
