export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USERNAME}/.cmsgit-cache

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_0_29
cd CMSSW_8_0_29/src
cmsenv

git cms-init

### MET re-correction (https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription)
git cms-merge-topic cms-met:METRecipe_8020_for80Xintegration

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

### Photon ID (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recipe_for_regular_users_for_8_0)
git cms-merge-topic ikrav:egm_id_80X_v3_photons
cd $CMSSW_BASE/src

### Tau MVA ID (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_MiniA)
git cms-merge-topic -u cms-tau-pog:CMSSW_8_0_X_tau-pog_tauIDOnMiniAOD-legacy-backport-81Xv2

### SVfit 
git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone 
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone 
git checkout HIG-16-006 
cd ${CMSSW_BASE}/src/ 

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

### DesyTau 
cd ${CMSSW_BASE}/src 
git clone https://github.com/DesyTau/DesyTauAnalysesRun2_25ns.git ${CMSSW_BASE}/src/DesyTauAnalyses 
cd ${CMSSW_BASE}/src/DesyTauAnalyses
git checkout NTuple_80X_MiniAODv2
cd ${CMSSW_BASE}/src 

### patch for SVFit 
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.h TauAnalysis/SVfitStandalone/interface/ 
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.cc TauAnalysis/SVfitStandalone/src 
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/testSVfitStandalone.cc TauAnalysis/SVfitStandalone/bin 
rm TauAnalysis/SVfitStandalone/interface/SVfitStandaloneQuantities.h 
rm TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc 

scram b -j 32 
scram b -j 32