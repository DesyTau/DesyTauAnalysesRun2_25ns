#variable=met_MTsum_14

#file=Templates_met_MTsum_14_2320invfb
#lumi="2320"
## syntax = cards.sh FILE_WITH_2D.root text_with_signal
while read line 
do

file=$1.root

channel=et

if [[ ! -f $file ]]; then

echo the $1 is not a valid input...
return 0
fi

variable=`echo $1  | cut -d "_" -f2-4`

lumi=`echo $1 | cut -d "_" -f5 | cut -d "i" -f1`



if [[ ! -f cards_${channel}/${line}_${variable}_${channel}_${lumi}invfb.txt ]] ;then

echo Will create the card for channel $channel, signal $line, variable $variable and lumi $lumi...

cp CreateDatacards_C CreateDatacards.C

sed  -i  's/CHANNEL/'$channel'/g' CreateDatacards.C
sed  -i  's/FILE/'$file'/g' CreateDatacards.C
sed  -i  's/VARIABLE/'$variable'/g' CreateDatacards.C
sed  -i  's/SIGNAL/'$line'/g' CreateDatacards.C
sed  -i  's/LUMI/'$lumi'/g' CreateDatacards.C



root -q -b -l CreateDatacards.C
fi


done<$2
