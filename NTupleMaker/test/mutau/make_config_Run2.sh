#!/bin/bash

### Important:
### the script is to be run with "bash make_config_Run2.sh <year={16,17,18}> <data_type={data, MC}>"

YEAR=$1
DATA_TYPE=$2

if [[ $YEAR -ne 16 && $YEAR -ne 17 && $YEAR -ne 18 ]]; then
  echo
  echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
  echo
  echo "  bash make_config_Run2.sh <year={16,17,18}> <data_type={data, MC}>"
  echo
  echo "year is not 16, 17 or 18 - exiting"
  exit
fi

# pick the template for the corresponding year
TEMPLATE_CFG_NAME="analysisMacroSynch_mt_${YEAR}"

# define parameters which are different between MC and data configs
KEY_LIST=(isData ApplyPUweight ApplyLepSF)
VALUE_LIST_MC=(false true true)
VALUE_LIST_DATA=(true false false)

# redefine list of the parameters according to the input data type
if [[ $DATA_TYPE == "data" ]]; then
  VALUE_LIST=("${VALUE_LIST_DATA[@]}")
  NOT_DATA_TYPE="MC"
else
  if [[ $DATA_TYPE == "MC" ]]; then
    VALUE_LIST=("${VALUE_LIST_MC[@]}")
    NOT_DATA_TYPE="data"
  else
    echo
    echo "To produce the scripts for a specific year and either data or MC this script is to be run with a command:"
    echo
    echo "  bash make_config_Run2.sh <year={16,17,18}> <data_type={data, MC}>"
    echo
    echo "data_type is neither data nor MC - exiting"
    exit
  fi
fi

# these parameters are year dependant for MC, so leave them as they are in the config and set to 0 only if it is data config
if [[ $DATA_TYPE == "data" ]]; then
  KEY_LIST+=(TauEnergyScaleShift_OneProng TauEnergyScaleShift_OneProngOnePi0 TauEnergyScaleShift_ThreeProng)
  VALUE_LIST_DATA+=(0.0 0.0 0.0)
  
  KEY_LIST+=(TauEnergyScaleShift_LepFake_OneProng TauEnergyScaleShift_LepFake_OneProngOnePi0 TauEnergyScaleShift_LepFake_ThreeProng)
  VALUE_LIST_DATA+=(0.0 0.0 0.0)
fi

# add the parameters to the config
KEY_LEN=${#KEY_LIST[@]}
for (( i = 0; i < $KEY_LEN; i++ )); do
        printf '/%s/c\%s\n' "${KEY_LIST[i]} =*" "${KEY_LIST[i]} = ${VALUE_LIST[i]}"
done | sed -r -f- $TEMPLATE_CFG_NAME.conf > ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf

# remove all the lines which starts with "NOT_DATA_TYPE: " 
sed -i "/${NOT_DATA_TYPE}: /d" ./${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf

# remove just the strings "DATA_TYPE: " leaving the rest of the line intact 
sed -i "s/${DATA_TYPE}: //" ./${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf

# lists with the MC samples' names
MC_SAMPLES_LIST=(DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8)
MC_SAMPLES_LIST+=(DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8)

MC_SAMPLES_LIST+=(ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8 ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8)
MC_SAMPLES_LIST+=(ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8)

MC_SAMPLES_LIST+=(WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8)
MC_SAMPLES_LIST+=(W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8)

MC_SAMPLES_LIST+=(TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8 TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8 TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8)
MC_SAMPLES_LIST+=(WW_TuneCP5_13TeV-pythia8 WZ_TuneCP5_13TeV-pythia8 ZZ_TuneCP5_13TeV-pythia8)
MC_SAMPLES_LIST+=(GluGluHToTauTau_M125_13TeV_powheg_pythia8 VBFHToTauTau_M125)
MC_SAMPLES_LEN=${#MC_SAMPLES_LIST[@]}

# Path in the root file to PU histograms for 16 and 18 data; 
PU_STR=pileup

if [[ $DATA_TYPE == "MC" ]]; then
  for (( i = 0; i < $MC_SAMPLES_LEN; i++ )); do
      if [[ $YEAR == 17 ]]; then
        # for 17 it is sample dependent, pick it from the list
        PU_STR=${MC_SAMPLES_LIST[i]}_pileup
      fi
      sed "s/pileUpforMC =/pileUpforMC = ${PU_STR}/" ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_${MC_SAMPLES_LIST[i]}.conf
  done
  sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf
fi
