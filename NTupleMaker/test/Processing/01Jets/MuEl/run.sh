
index="19"
titles="hDZeta"
rm  SignalTemplates_$titles.root 
while read line
do

echo for $line..

cp GetHisto_C GetHisto.C

sed -i '' 's/HISTOHERE/'$titles'/g' GetHisto.C

sed -i '' 's/FILEIN/'$line'/g' GetHisto.C

sed -i '' 's/INDEX/'$index'/g' GetHisto.C

root -l -q -b GetHisto.C

done<$1
