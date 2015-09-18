channel=$2
while read line
do

cat bss > job$line$channel.sh
echo 	SUSY$channel analysisMacroSUSY.conf $line >> job$line$channel.sh

chmod u+x job$line$channel.sh
qsub job$line$channel.sh

##rm job$line.sh

done<$1
