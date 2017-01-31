channel=$2 
dir=$2




##type MC or Data
type=Data



cp *.conf Jobs/.

while read line
do

unset xsec
xsec=`grep " ${line} " xsecs | cut -d " " -f3-4`	
#	cp $dir/${line} input$line
#xsec=1
echo FOUND XSEC for ${line} to be $xsec
unset f
while read f
	do

echo $f
unset bas
bas=`basename $f | awk -F ".root" '{print $1}'`

if [ ! -f $dir/$bas.root ] 
then
echo $f > $dir/$bas

#	echo " "$bas $xsec >> xsecs

	if [ -f Jobs/job${line}$channel$dir${bas}_B.sh ] ; then
rm Jobs/job${line}$channel$dir${bas}_B.sh
fi


cat bss > Jobs/job${line}$channel$dir${bas}_B.sh

echo SUSY$channel analysisMacroSUSY_${type}_B.conf ${bas} $dir 1 1 1 1>> Jobs/job${line}$channel$dir${bas}_B.sh
#echo SUSY${channel} analysisMacroSUSY_${type}_B.conf ${bas} $dir>> Jobs/job${line}$channel$dir${bas}_B.sh



#echo $bas $xsec >> xsecs
echo $bas $xsec 

if [ ! -f $dir/${bas}_B_OS.root ] ;then

chmod u+x Jobs/job${line}$channel$dir${bas}_B.sh
echo . Jobs/job${line}$channel$dir${bas}_B.sh 
rm Jobs/job${line}$channel$dir${bas}_B.sh 

fi


fi

done<$dir/${line}
rm input${line}
done<$1



