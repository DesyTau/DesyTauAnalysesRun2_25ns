#!/bin/bash

### Important:
### the script is to be run with "bash make_config_Run2.sh"

YEAR=$1
DATA_TYPE=$2

if [[ $YEAR -ne 16 && $YEAR -ne 17 && $YEAR -ne 18 ]]; then
  echo "year is not 16, 17 or 18 - exiting"
  exit
fi

# pick the corresponding template for the year
TEMPLATE_CFG_NAME="analysisMacroSynch_mt_${YEAR}"

# define parameters which are different between MC and data configs
KEY_LIST=(isData ApplyPUweight ApplyLepSF)
VALUE_LIST_MC=(false true true)
VALUE_LIST_DATA=(true false false)

# these parameters are year dependant for MC, so leave them as they are in the config and set to 0 only if it is data config
if [[ $DATA_TYPE == "data" ]]; then
  KEY_LIST+=(TauEnergyScaleShift_OneProng TauEnergyScaleShift_OneProngOnePi0 TauEnergyScaleShift_ThreeProng)
  VALUE_LIST_DATA+=(0.0 0.0 0.0)
  
  KEY_LIST+=(TauEnergyScaleShift_LepFake_OneProng TauEnergyScaleShift_LepFake_OneProngOnePi0 TauEnergyScaleShift_LepFake_ThreeProng)
  VALUE_LIST_DATA+=(0.0 0.0 0.0)
fi

# redefine list of the parameters according to the input data type
if [[ $DATA_TYPE == "data" ]]; then
  VALUE_LIST=("${VALUE_LIST_DATA[@]}")
  NOT_DATA_TYPE="MC"
else
  if [[ $DATA_TYPE == "MC" ]]; then
    VALUE_LIST=("${VALUE_LIST_MC[@]}")
    NOT_DATA_TYPE="data"
  else
    echo "data_type is neither data or MC - exiting"
    exit
  fi
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

cp analysisMacroSynch_lept_mt_MC17.conf analysisMacroSynch_lept_mt_backup.conf

# sed 's/pileUpforMC =/pileUpforMC = DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
sed 's/pileUpforMC =/pileUpforMC = TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = MC_PU2017_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_DYJetsToLL_M-5to50_13TeV-12Apr2018.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = WW_TuneCP5_13TeV-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_WW_TuneCP5_13TeV-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = WZ_TuneCP5_13TeV-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_WZ_TuneCP5_13TeV-pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = ZZ_TuneCP5_13TeV-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_ZZ_TuneCP5_13TeV-pythia8.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_GluGluHToTauTau_M125_13TeV_powheg_pythia8.conf
# sed 's/pileUpforMC =/pileUpforMC = GluGluHToTauTau_M125_13TeV_powheg_pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf
# 
# #Merijn: add a line for VBF:
# sed 's/pileUpforMC =/pileUpforMC = VBFHToTauTau_M125_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_VBF125.conf
# 
# 
# sed 's/pileUpforMC =/pileUpforMC = MC_PU2017_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_WGToLNuG.conf
# 
# sed 's/pileUpforMC =/pileUpforMC = EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_EWKWMinus.conf
# sed 's/pileUpforMC =/pileUpforMC = EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_EWKWPlus.conf
# sed 's/pileUpforMC =/pileUpforMC =  EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_EWKZ2Jets_ZToLL.conf
# sed 's/pileUpforMC =/pileUpforMC =  EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_pileup/' ${TEMPLATE_CFG_NAME}_${DATA_TYPE}.conf > analysisMacroSynch_lept_mt_EWKZ2Jets_ZToNuNu.conf
