#!/bin/bash

### the script is to be run with "./make_config_Run2.sh <year={16,17,18}> <data_type={data, MC, embedded}>"

YEAR=$1
DATA_TYPE=$2
CHANNEL=em
if [ $YEAR -lt 16 ] || [ $YEAR -gt 18 ]; then
    echo
    echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
    echo
    echo "  ./make_config_Run2.sh <year={16,17,18}> <data_type={data, MC, embedded}> <channel={mt,et,em}>"
    echo
    echo "year is not 16, 17 or 18 - exiting"
    exit
fi

if [ $DATA_TYPE != 'data' ] && [ $DATA_TYPE != 'MC' ] && [ $DATA_TYPE != 'embedded' ]; then
    echo
    echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
    echo
    echo "  ./make_config_Run2.sh <year={16,17,18}> <data_type={data, MC, embedded}>"
    echo
    echo "data_type is not {data, MC, embedded} - exiting"
    exit

fi
OUTDIR=20$YEAR
if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

if [[ $YEAR -eq 16 ]]; then
  NOT_YEAR=(17 18)
else 
  if [[ $YEAR -eq 17 ]]; then
    NOT_YEAR=(16 18)
  else 
    if [[ $YEAR -eq 18 ]]; then
      NOT_YEAR=(16 17)
    else 
      echo
      echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
      echo
      echo "  ./make_config_Run2.sh <year={16,17,18}> <data_type={data, MC, embedded}> <channel={mt,et}>"
      echo
      echo "year is not 16, 17 or 18 - exiting"
      exit
    fi
  fi
fi

TEMPLATE_CFG_PREFIX="analysisMacroSynch"
TEMPLATE_CFG_NAME=${TEMPLATE_CFG_PREFIX}_${CHANNEL}_${YEAR}_${DATA_TYPE}
cp ${TEMPLATE_CFG_PREFIX}_template.conf ${TEMPLATE_CFG_NAME}_tmp.conf

# remove all the lines for years which is not the one specified 
NOT_YEAR_LEN=${#NOT_YEAR[@]}
for (( i = 0; i < NOT_YEAR_LEN; i++ )); do
  sed -i "/${NOT_YEAR[i]}: /d" ${TEMPLATE_CFG_NAME}_tmp.conf
done

# remove just the strings "YEAR: " leaving the rest of the line intact 
sed -i "s/${YEAR}: //" ${TEMPLATE_CFG_NAME}_tmp.conf

# define parameters which are different between MC and data configs
KEY_LIST=(isData ApplyPUweight ApplyLepSF ApplyRecoilCorrections ApplyBTagScaling)
VALUE_LIST_MC=(false true true true true)
VALUE_LIST_DATA=(true false false false false)
VALUE_LIST_EMBEDDED=(true false true false false)

if [[ $DATA_TYPE == "data" ]]; then
  VALUE_LIST=("${VALUE_LIST_DATA[@]}")
  NOT_DATA_TYPE=("MC" "embedded")
else
  if [[ $DATA_TYPE == "MC" ]]; then
    VALUE_LIST=("${VALUE_LIST_MC[@]}")
    NOT_DATA_TYPE=("data" "embedded")
  else    
    if [[ $DATA_TYPE == "embedded" ]]; then          
      VALUE_LIST=("${VALUE_LIST_EMBEDDED[@]}")
      NOT_DATA_TYPE=("MC" "data")
    else
      echo
      echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
      echo
      echo "  ./make_config_Run2.sh <year={16,17,18}> <data_type={data, MC, embedded}> <channel={mt,et}>"
      echo
      echo "data_type is neither data nor MC - exiting"
      exit
    fi  
  fi
fi

# these parameters are year dependant for MC, so leave them as they are in the config and set to 0 only if it is data config
# also redefine list of the parameters according to the input data type

# add the parameters to the config
KEY_LEN=${#KEY_LIST[@]}
for (( i = 0; i < $KEY_LEN; i++ )); do
        printf '/%s/c\%s\n' "${KEY_LIST[i]} =*" "${KEY_LIST[i]} = ${VALUE_LIST[i]}"
done | sed -r -i -f- ${TEMPLATE_CFG_NAME}_tmp.conf

# add year
sed -i "s/era = /era = 20$YEAR/" ${TEMPLATE_CFG_NAME}_tmp.conf

# remove all the lines which starts with "NOT_DATA_TYPE: " 
NOT_DATA_TYPE_LEN=${#NOT_DATA_TYPE[@]}
for (( i = 0; i < NOT_DATA_TYPE_LEN; i++ )); do
  sed -i "/${NOT_DATA_TYPE[i]}: /d" ${TEMPLATE_CFG_NAME}_tmp.conf
done

# remove just the strings "DATA_TYPE: " leaving the rest of the line intact 
sed -i "s/${DATA_TYPE}: //" ${TEMPLATE_CFG_NAME}_tmp.conf

# lists with the MC samples' names
MC_SAMPLES_LIST=(DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-10to50_13TeV-12Apr2018)
MC_SAMPLES_LIST+=(DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8)
MC_SAMPLES_LIST+=(W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8)
MC_SAMPLES_LIST+=(TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8 TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8)
MC_SAMPLES_LIST+=(WW_TuneCP5_13TeV-pythia8 WZ_TuneCP5_13TeV-pythia8 ZZ_TuneCP5_13TeV-pythia8)

for susyname in $(less list_SUSY_ggH_2017);
do
    MC_SAMPLES_LIST+=(${susyname} )
done

MC_SAMPLES_LEN=${#MC_SAMPLES_LIST[@]}

if [[ $DATA_TYPE == "MC" ]]; then
  if [[ $YEAR -eq 17 ]]; then # for 17 the path in the root file to PU histograms is sample-dependent, pick it from the list
    for (( i = 0; i < $MC_SAMPLES_LEN; i++ )); do
        PU_STR=${MC_SAMPLES_LIST[i]}_pileup
        sed "s/pileUpforMC =/pileUpforMC = ${PU_STR}/" ${TEMPLATE_CFG_NAME}_tmp.conf > $OUTDIR/${TEMPLATE_CFG_NAME}_${MC_SAMPLES_LIST[i]}.conf
    done
    sed "s/pileUpforMC =/pileUpforMC = MC_PU2017_pileup/" ${TEMPLATE_CFG_NAME}_tmp.conf > $OUTDIR/${TEMPLATE_CFG_NAME}.conf
  else
    sed "s/pileUpforMC =/pileUpforMC = pileup/" ${TEMPLATE_CFG_NAME}_tmp.conf > $OUTDIR/${TEMPLATE_CFG_NAME}.conf
  fi # path in the root file to PU histograms for 16 and 18 data; 
else
  cp  ${TEMPLATE_CFG_NAME}_tmp.conf $OUTDIR/${TEMPLATE_CFG_NAME}.conf
fi
rm ${TEMPLATE_CFG_NAME}_tmp.conf

if [[ $DATA_TYPE == "data" ]]; then
    if [[ $YEAR -eq 16 ]]; then
	sed "s/ApplyDzFilterMatch = false/ApplyDzFilterMatch = true/" $OUTDIR/${TEMPLATE_CFG_NAME}.conf > $OUTDIR/${TEMPLATE_CFG_NAME}GH.conf
    fi
    if [[ $YEAR -eq 17 ]]; then
	sed "s/LowPtLegMuon =  hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8/LowPtLegMuon =  hltL3fL1sMu7EG23f0Filtered8/" $OUTDIR/${TEMPLATE_CFG_NAME}.conf > $OUTDIR/${TEMPLATE_CFG_NAME}B.conf
    fi	
fi
