# Updated instructions to synchronize your area with GitHub -102x starting from 10_2_14 for 2018 ntuple production.
# Before starting, it is required that you already have a github account. Current stable version 10_2_14.

export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700

cmsrel CMSSW_10_2_14
cd CMSSW_10_2_14/src
cmsenv

git cms-init

### EGamma energy corrections (from here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2018_Preliminary_Energy_Correcti)
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone https://github.com/cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
scram b -j 8

### DeepTauId v2  (from here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Running_of_the_DNN_based_tau_ID)
#git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2 
#Add 2017v2 training file by using "git clone" or wget 
#git clone -b DeepTau2017v2 https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles.git RecoTauTag/TrainingFiles/data 
#wget https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles/raw/DeepTau2017v2/DeepTauId/deepTau_2017v2p6_e6_core.pb -P RecoTauTag/TrainingFiles/data/DeepTauId 
#wget https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles/raw/DeepTau2017v2/DeepTauId/deepTau_2017v2p6_e6_inner.pb -P RecoTauTag/TrainingFiles/data/DeepTauId 
#wget https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles/raw/DeepTau2017v2/DeepTauId/deepTau_2017v2p6_e6_outer.pb -P RecoTauTag/TrainingFiles/data/DeepTauId

### SVfit (needed for backward compatibility)
git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd ${CMSSW_BASE}/src/

### SVfit-updated
git clone https://github.com/svfit/ClassicSVfit TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

### Lepton efficiencies interface and ROOT files from CMS-HTT
git clone https://github.com/CMS-HTT/LeptonEff-interface HTT-utilities
cd HTT-utilities/LepEffInterface/
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data
cd ${CMSSW_BASE}/src

### Corrections workspace from CMS-HTT
git clone https://github.com/CMS-HTT/CorrectionsWorkspace HTT-utilities/CorrectionsWorkspace
cd ${CMSSW_BASE}/src/HTT-utilities/CorrectionsWorkspace
root -l -q CrystalBallEfficiency.cxx++
cd ${CMSSW_BASE}/src

### RecoilCorrections
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections

### QCD background model in emu channel
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/QCDModelingEMu.git HTT-utilities/QCDModelingEMu

### Tau trigger SF (clarify who is needing this)
cd ${CMSSW_BASE}/src
git clone -b tauTriggers2017_MCv2_PreReMiniaod https://github.com/truggles/TauTriggerSFs.git HTT-utilities/TauTriggerSFs2017

### DesyTau
cd ${CMSSW_BASE}/src
git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses

### patch for SVFitStandalone (for backward compatibility)
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.h TauAnalysis/SVfitStandalone/interface/
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.cc TauAnalysis/SVfitStandalone/src
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/testSVfitStandalone.cc TauAnalysis/SVfitStandalone/bin
rm TauAnalysis/SVfitStandalone/interface/SVfitStandaloneQuantities.h
rm TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc

scram b -j 32
scram b -j 32
