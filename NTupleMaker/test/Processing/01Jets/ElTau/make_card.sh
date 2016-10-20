
while read line
do

python datacard.py  18 $line  > mutau_$line.txt
done <$1
