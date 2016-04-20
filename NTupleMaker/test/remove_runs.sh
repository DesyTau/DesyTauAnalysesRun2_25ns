while read line
do
echo	sed -i "'"/$line/d"'" temp
done<$1
