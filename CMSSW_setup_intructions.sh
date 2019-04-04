export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USERNAME}/.cmsgit-cache

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700

cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv

git cms-init

### SVfit (needed for backward compatibility)
git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd ${CMSSW_BASE}/src/

### SVfit-updated
git clone https://github.com/svfit/ClassicSVfit TauAnalysis/ClassicSVfit  -b release_2018Mar20
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
cd ${CMSSW_BASE}/src 

### QCD background model in emu channel 
cd ${CMSSW_BASE}/src 
git clone https://github.com/CMS-HTT/QCDModelingEMu.git HTT-utilities/QCDModelingEMu 
cd ${CMSSW_BASE}/src 

### Tau trigger SF (clarify who is needing this)
cd ${CMSSW_BASE}/src 
git clone -b tauTriggers2017_MCv2_PreReMiniaod https://github.com/truggles/TauTriggerSFs.git HTT-utilities/TauTriggerSFs2017
cd ${CMSSW_BASE}/src 

### DesyTau 
cd ${CMSSW_BASE}/src
git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses -b 2016-legacy
cd ${CMSSW_BASE}/src

### Refitted vertex (needed for CP analysis)
cd ${CMSSW_BASE}/src
git clone https://github.com/aknayak/TauRefit.git VertexRefit/TauRefit
cd ${CMSSW_BASE}/src

### patch for SVFitStandalone (for backward compatibility)
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.h TauAnalysis/SVfitStandalone/interface/
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.cc TauAnalysis/SVfitStandalone/src
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/testSVfitStandalone.cc TauAnalysis/SVfitStandalone/bin
rm TauAnalysis/SVfitStandalone/interface/SVfitStandaloneQuantities.h
rm TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc

scram b -j 32
scram b -j 32