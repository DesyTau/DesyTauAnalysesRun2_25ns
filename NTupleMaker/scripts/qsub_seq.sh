#!/bin/bash
# $1 - executable
# $2 - analysis macro
# $3 - filelist
# $4 - nfiles [optional]
# $5 - channel [optional]

EXECUTABLE=""
CONFIGFILE=""
LIST=""
CHANNEL=""
NFILES="20"
SCRIPT=""

for i in "$@"
do
case $i in
    -e=*|--exe=*)
    EXECUTABLE="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--config=*)
    CONFIGFILE="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|--list=*)
    LIST="${i#*=}"
    shift # past argument=value
    ;;
    --ch=*)
    CHANNEL="${i#*=}"
    shift # past argument=value
    ;;
    --nfiles=*)
    NFILES="${i#*=}"
    shift # past argument=value
    ;;
    --script=*)
    SCRIPT="${i#*=}"
    shift # past argument=value
    ;;    
    *)
            # unknown option
    ;;
esac
done

if [ -z "$EXECUTABLE" ]
then
    echo "No executable provided. Use -e option. Exiting.."
    exit
fi

if [ -z "$CONFIGFILE" ]
then
    echo "No config file provided. Use -c option. Exiting.."
    exit
fi

if [ -z "$LIST" ]
then
    echo "No file list provided. Use -l option. Exiting.."
    exit
fi

SAMPLE=${LIST##*/}

let "n = 0"
rm -rf ${SAMPLE}_files
mkdir ${SAMPLE}_files

if [ -z "$SCRIPT" ]
then
    echo "No script template provided. Using default."
    cp $CMSSW_BASE/src/DesyTauAnalyses/NTupleMaker/scripts/qsub.sh ${SAMPLE}_files/
else
    cp $SCRIPT ${SAMPLE}_files/
fi

cp $CONFIGFILE ${SAMPLE}_files/
cp $LIST .
Splitter $SAMPLE $NFILES
rm $SAMPLE

cd ${SAMPLE}_files/
for i in `ls $SAMPLE_*`
 do
 echo submitting job $n for file $i from list $LIST
 ./qsub.sh $EXECUTABLE $CONFIGFILE $i $CHANNEL
 let "n = n + 1"
done
echo Total $n jobs submitted
cd ../
