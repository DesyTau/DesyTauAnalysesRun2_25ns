export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USERNAME}/.cmsgit-cache

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src
cmsenv

git cms-init

### MET re-correction (https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription)
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic ahinzmann:METRecipe_8020_Moriond17

### MET filters (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2)
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

### Bad Muon remedy (https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2786.html)
git cms-merge-topic gpetruc:badMuonFilters_80X_v2

### Pileup Jet Id (https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Information_for_13_TeV_data_anal)
cd ${CMSSW_BASE}/src
git cms-addpkg RecoJets/JetProducers
cd ${CMSSW_BASE}/src/RecoJets/JetProducers/data
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta0to2p5_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta2p5to2p75_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta2p75to3_BDT.weights.xml.gz
wget https://github.com/cms-data/RecoJets-JetProducers/raw/master/pileupJetId_80XvarFix_Eta3to5_BDT.weights.xml.gz
cd $CMSSW_BASE/src

### Electron Id (https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recipes_for_regular_users_common)
git cms-merge-topic ikrav:egm_id_80X_v2
scram b -j 32
cd $CMSSW_BASE/external
cd slc6_amd64_gcc530/
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_BASE/src

### Photon ID (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recipe_for_regular_users_for_8_0)
git cms-merge-topic ikrav:egm_id_80X_v3_photons
cd $CMSSW_BASE/src

### SVfit
git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd ${CMSSW_BASE}/src/

### Lepton efficiencies interface and ROOT files from CMS-HTT
git clone https://github.com/CMS-HTT/LeptonEff-interface  HTT-utilities
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

### DesyTau
cd ${CMSSW_BASE}/src
git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses

#### patch for MVAMET
# cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/MVAMEt/MVAMET.cc ${CMSSW_BASE}/src/RecoMET/METPUSubtraction/plugins/.

scram b -j 32
scram b -j 32
