while read line
do
	#rm $line.root
	#if [ ! -f $line.root ] 
	#then
		SUSYmuel analysisMacroSUSY.conf $line
	#fi
done<$1
