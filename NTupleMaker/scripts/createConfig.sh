#!/bin/zsh

if [ ! -s haddAll.sh ]
then
    echo "#!/bin/zsh" > haddAll.sh
    echo "" >> haddAll.sh
fi

if [ ! -s runAll.sh ]
then
    echo "#!/bin/zsh" > runAll.sh
    echo "" >> runAll.sh
fi


for file in `ls $1 | grep -v tar.gz`
do 
   ls $1/$file/* > $file 
   echo hadd.sh $file '&' >> haddAll.sh
   echo qsub_seq.sh -e='$MACRO' -c='$CONFIG'.conf -l=$file --nfiles='$NFILES' --ch='$CHANNEL' >> runAll.sh
done

chmod +x haddAll.sh
chmod +x runAll.sh
