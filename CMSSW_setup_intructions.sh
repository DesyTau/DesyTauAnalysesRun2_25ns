export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USERNAME}/.cmsgit-cache

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv

git cms-init

## MVA MET
git cms-addpkg RecoMET/METPUSubtraction
git cms-addpkg DataFormats/METReco
git remote add -f mvamet https://github.com/rfriese/cmssw.git
git checkout mvamet/mvamet8020 -b mvamet
mkdir RecoMET/METPUSubtraction/data
cd RecoMET/METPUSubtraction/data
wget https://github.com/rfriese/cmssw/raw/MVAMET2_beta_0.6/RecoMET/METPUSubtraction/data/weightfile.root
cd $CMSSW_BASE/src

## MET filters
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git cms-merge-topic cms-met:METRecipe_8020

### Pileup Jet Id
cd ${CMSSW_BASE}/src
git cms-addpkg RecoJets/JetProducers
cd ${CMSSW_BASE}/src/RecoJets/JetProducers/data
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta0to2p5_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta2p5to2p75_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta2p75to3_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta3to5_BDT.weights.xml.gz
cd ${CMSSW_BASE}/src

## SVfit
cd $CMSSW_BASE/src
git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd ${CMSSW_BASE}/src/

## Lepton efficiencies interface and ROOT files from CMS-HTT
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/LeptonEff-interface  HTT-utilities
cd HTT-utilities/LepEffInterface/
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data

## Corrections workspace from CMS-HTT
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/CorrectionsWorkspace HTT-utilities/CorrectionsWorkspace
cd ${CMSSW_BASE}/src/HTT-utilities/CorrectionsWorkspace
root -l -q CrystalBallEfficiency.cxx++
cd ${CMSSW_BASE}/src

## RecoilCorrections
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections 

### QCD background model in emu channel
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/QCDModelingEMu.git HTT-utilities/QCDModelingEMu

### DesyTau
cd ${CMSSW_BASE}/src
git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses
cd ${CMSSW_BASE}/src/DesyTauAnalyses
#git checkout NTuple_80X_MiniAODv2
cd ${CMSSW_BASE}/src/

#### patch
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/MVAMEt/MVAMET.cc ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/plugins/.

scram b -j 32
scram b -j 32
