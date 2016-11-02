while read line 
do

n=`echo $line | awk -F "_LSP" '{print $1}'`
xs=`echo $line  | awk -F "_LSP" '{print $2}' | awk -F " " '{print $2}'`
xss=`echo $line  | awk -F "_B.root" '{print $1}'` 
xs2=`grep " $xss " xsecs | awk -F " " '{print $2}'`

echo " "$xss $xs2

done<$1
