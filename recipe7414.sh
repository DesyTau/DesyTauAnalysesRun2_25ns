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

## PileUp Jet Id patch
## ref. https://github.com/cms-sw/cmssw/pull/12398
git-cms-addpkg DataFormats/JetReco
git-cms-addpkg RecoJets/JetProducers

rm -rf ${CMSSW_BASE}/src/RecoJets/JetProducers/data
git clone https://github.com/cms-data/RecoJets-JetProducers.git ${CMSSW_BASE}/src/RecoJets/JetProducers/data

cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/PileUpJetID/DataFormats/JetReco/interface/PileupJetIdentifier.h ${CMSSW_BASE}/src/DataFormats/JetReco/interface/PileupJetIdentifier.h
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/PileUpJetID/RecoJets/JetProducers/plugins/PileupJetIdProducer.cc ${CMSSW_BASE}/src/RecoJets/JetProducers/plugins/PileupJetIdProducer.cc
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/PileUpJetID/RecoJets/JetProducers/python/PileupJetIDParams_cfi.py ${CMSSW_BASE}/src/RecoJets/JetProducers/python/PileupJetIDParams_cfi.py
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/PileUpJetID/RecoJets/JetProducers/src/PileupJetIdAlgo.cc ${CMSSW_BASE}/src/RecoJets/JetProducers/src/PileupJetIdAlgo.cc

## MVA MEt
git cms-addpkg RecoMET/METPUSubtraction/
cd ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/
git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
cd ${CMSSW_BASE}/src/
## ## adding the pairwise MVA MEt calculation
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/MVAMEt/PFMETProducerMVATauTau* ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/plugins/.

## MVA MEt covariance matrix fix as from Raphael's talk
findline='mvaOutputCovU1_ = GetResponse(mvaReaderCovU1_, varForCovU1_)\* mvaOutputU_ \* var_\["particleFlow_U"\]\;'
replaceline='mvaOutputCovU1_ = std::pow(GetResponse(mvaReaderCovU1_, varForCovU1_)\* mvaOutputU_ \* var_\["particleFlow_U"\], 2)\;'
sed -i -e 's/'$findline'/'$replaceline'/' ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/src/PFMETAlgorithmMVA.cc

findline='mvaOutputCovU2_ = GetResponse(mvaReaderCovU2_, varForCovU2_)\* mvaOutputU_ \* var_\["particleFlow_U"\]\;'
replaceline='mvaOutputCovU2_ = std::pow(GetResponse(mvaReaderCovU2_, varForCovU2_)\* mvaOutputU_ \* var_\["particleFlow_U"\], 2)\;'
sed -i -e 's/'$findline'/'$replaceline'/' ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/src/PFMETAlgorithmMVA.cc

## Lepton efficiencies interface and ROOT files from CMS-HTT 
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/LeptonEff-interface.git HTT-utilities 
cd HTT-utilities/LepEffInterface/
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data 

cd ${CMSSW_BASE}/src

scram b -j 32
scram b -j 32

