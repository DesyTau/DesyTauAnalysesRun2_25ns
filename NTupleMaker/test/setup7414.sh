export USERNAME=YOURUSERNAME

export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USERNAME}/.cmsgit-cache

cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
cmsenv

git cms-init

git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses
cd ${CMSSW_BASE}/src/DesyTauAnalyses/. 
cd ${CMSSW_BASE}/src/

git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout svFit_2015Apr03
cd ${CMSSW_BASE}/src/

## Electron id
git cms-merge-topic ikrav:egm_id_7.4.12_v1

## MVA MEt
git cms-addpkg RecoMET/METPUSubtraction/
cd ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/
git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
cd ${CMSSW_BASE}/src/
## ## adding the pairwise MVA MEt calculation
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/MVAMEt/PFMETProducerMVATauTau* ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/plugins/.

## MEt corrections
##git cms-merge-topic -u cms-met:METCorUnc74X
scram b -j 32


