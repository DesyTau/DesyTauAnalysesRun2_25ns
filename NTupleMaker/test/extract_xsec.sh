while read line
do

#line2=`echo $line | awk -F " " '{print $1}' | awk -F "_LSP" '{print $1}'`
line2=`echo $line | awk -F " " '{print $1}'`

line3=`grep "$line2 " xsecs`
line4=`grep "$line" xsecs | awk -F "_LSP" '{print $1}'`
#echo line 2 $line2

dt=`echo $line3 | awk -F " " '{print $1}'`
dt2=`echo $line4 | awk -F "_LSP" '{print $1}'`
lsp=`echo $line4 | awk -F "_LSP" '{print $2}'`
xsec=`echo $line3 | awk -F " " '{print $2}'`
xsec2=`echo $line3 | awk -F " " '{print $3}'`

if [[ $2 == "qcd" ]] ;then
echo ${dt}_A $xsec $xsec2
echo ${dt}_C $xsec $xsec2
echo ${dt}_D $xsec $xsec2
fi

if [[ $2 != "qcd" ]] ;then

echo " "${dt} $xsec $xsec2


fi


done<$1

