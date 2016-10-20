

listy=`cat $1| awk -F " " '{print $3}' | tr '\n\'`
listx=`cat $1| awk -F " " '{print $2}' |  tr '\n' ','`



echo listx={ $listx}
echo listy={ $listy}

echo $listy $listx

while read line
do
        echo will work for $line
        echo $line > ds

root -l -b -q GetEff.C
done<$1

