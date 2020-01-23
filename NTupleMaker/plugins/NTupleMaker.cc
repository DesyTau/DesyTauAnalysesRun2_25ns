#include "DesyTauAnalyses/NTupleMaker/plugins/NTupleMaker.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisByIntegration.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>
#include <DataFormats/METReco/interface/GenMET.h>
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
//#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DesyTauAnalyses/CandidateTools/interface/candidateAuxFunctions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/idAlgos.h"

//#include "ICTauSpinnerProducer.hh"


#include <TString.h>
#include <Compression.h>

using namespace reco;
using namespace pat;
//typedef std::vector<NSVfitEventHypothesisByIntegration> NSVfitEventHypothesisByIntegrationCollection;
typedef ROOT::Math::XYZVector Vector;

// http://root.cern.ch/root/html/ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___.html#ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___:_SMatrix_double_3_3_ROOT__Math__MatRepSym_double_3___
// http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SVectorDoc.html
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;  // Standard Matrix representation for a general 3x3 matrix of type double.
typedef ROOT::Math::SVector<double, 3> SVector3; //// SVector: vector of size 3

// static const unsigned int SKIM_MUTAUTAU   = (1 << 0);     //1  : mu+tau+tau
// static const unsigned int SKIM_ETAUTAU    = (1 << 1);     //2  : e+tau+tau
// static const unsigned int SKIM_MUMU       = (1 << 2);     //4  : mu+mu
// static const unsigned int SKIM_EE         = (1 << 3);     //8  : e+e
// static const unsigned int SKIM_ETAU       = (1 << 4);     //16 : e+tau
// static const unsigned int SKIM_ALL        = (1 << 5);     //32 : all
// static const unsigned int SKIM_EMU        = (1 << 6);     //64 : e+mu
// static const unsigned int SKIM_TAUTAU     = (1 << 7);     //128: tau+tau

#include "DesyTauAnalyses/NTupleMaker/interface/genMatch.h"

void NTupleMaker::computePCA(double * pv, double * refPoint, double * mom, double * pca) {

  double diff[3];
  double norm = 0;
  for (int i=0; i<3; ++i) {
    norm += mom[i]*mom[i];
  }

  norm = TMath::Sqrt(norm);
  double n[3];
  double time = 0;
  for (int i=0; i<3; ++i) {
    n[i] = mom[i]/norm;
    time += n[i]*(pv[i]-refPoint[i]);
  }
  for (int i=0; i<3; ++i) {
    pca[i] = refPoint[i]+time*n[i];
  }

}

Int_t NTupleMaker::find_lep(const Int_t nlep, const Float_t px[], const Float_t py[], const Float_t pz[], const reco::Candidate::LorentzVector& refp4){
  Int_t ilep = -1;
  Float_t err = 1.e-4;

  for (int il= 0; il < nlep; il++){
    if( fabs(1 - refp4.px() / px[il]) > err) continue;
    if( fabs(1 - refp4.py() / py[il]) > err) continue;
    if( fabs(1 - refp4.pz() / pz[il]) > err) continue;

    if (ilep < 0)
      ilep=il;
    else
      std::cout<<"WARNING!!: Pairwise MVAMEt: duplicated lepton found!"<<std::endl;
  }

  //if(ilep < 0)
  //  std::cout<<"WARNING!!: Pairwise MVAMEt: Lep not found, pt = "<<refp4.pt()<<std::endl;

  return ilep;
}


//to set the values from parameter
NTupleMaker::NTupleMaker(const edm::ParameterSet& iConfig) :
  // data, year, period, skim
  cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
  cembedded(iConfig.getUntrackedParameter<bool>("IsEmbedded", false)),
  cFastSim(iConfig.getUntrackedParameter<bool>("IsFastSim", false)),
  cYear(iConfig.getUntrackedParameter<unsigned int>("Year")),
  cPeriod(iConfig.getUntrackedParameter<std::string>("Period")),
  cSkim(iConfig.getUntrackedParameter<unsigned int>("Skim")),
  // switches (collections)
  cgen(iConfig.getUntrackedParameter<bool>("GenParticles", false)),
  csusyinfo(iConfig.getUntrackedParameter<bool>("SusyInfo", false)),
  ctrigger(iConfig.getUntrackedParameter<bool>("Trigger", false)),
  cbeamspot(iConfig.getUntrackedParameter<bool>("RecBeamSpot", false)),
  crectrack(iConfig.getUntrackedParameter<bool>("RecTrack", false)),
  crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
  crecprimvertexwithbs(iConfig.getUntrackedParameter<bool>("RecPrimVertexWithBS", false)),
  crefittedvertex(iConfig.getUntrackedParameter<bool>("RefittedVertex", false)),
  crefittedvertexwithbs(iConfig.getUntrackedParameter<bool>("RefittedVertexWithBS", false)),
  crecTauSpinner(iConfig.getUntrackedParameter<bool>("ApplyTauSpinner", false)),
  crecmuon(iConfig.getUntrackedParameter<bool>("RecMuon", false)),
  crecelectron(iConfig.getUntrackedParameter<bool>("RecElectron", false)),
  crectau(iConfig.getUntrackedParameter<bool>("RecTau", false)),
  cl1objects(iConfig.getUntrackedParameter<bool>("L1Objects", false)),
  crecphoton(iConfig.getUntrackedParameter<bool>("RecPhoton", false)),
  crecpfjet(iConfig.getUntrackedParameter<bool>("RecJet", false)),
  crecpfpuppijet(iConfig.getUntrackedParameter<bool>("RecJetPuppi", false)),
  crecpfmet(iConfig.getUntrackedParameter<bool>("RecPFMet", false)),
  crecpfmetcorr(iConfig.getUntrackedParameter<bool>("RecPFMetCorr", false)),
  crecpuppimet(iConfig.getUntrackedParameter<bool>("RecPuppiMet", false)),
  crecmvamet(iConfig.getUntrackedParameter<bool>("RecMvaMet", false)),
  crecstxs(iConfig.getUntrackedParameter<bool>("RecHTXS", false)),
  // triggers
  cHLTriggerPaths(iConfig.getUntrackedParameter<vector<string> >("HLTriggerPaths")),
  cTriggerProcess(iConfig.getUntrackedParameter<string>("TriggerProcess", "HLT")),
  //TauSpinner
  cTauSpinAngles(iConfig.getParameter<std::string>("TauSpinnerAngles")),
  //cTSinput(iConfig.getUntrackedParameter<edm::InputTag>("TauSpinnerinput")),

  //Flags
  cFlags(iConfig.getUntrackedParameter<vector<string> >("Flags")),
  cFlagsProcesses(iConfig.getUntrackedParameter<vector<string> >("FlagsProcesses")),
  BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter"))),
  ecalBadCalibFilterUpdate_token(consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"))),
  //ecalBadCalibFilterUpdate_token(iConfig.getUntrackedParameter<bool>("ecalBadCalibReducedMINIAODFilterTag")),
  //ecalBadCalibFilterUpdate_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilterTag"))),
  // muons
  cMuPtMin(iConfig.getUntrackedParameter<double>("RecMuonPtMin", 10.)),
  cMuEtaMax(iConfig.getUntrackedParameter<double>("RecMuonEtaMax", 2.5)),
  cMuHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecMuonHLTriggerMatching")),
  cMuNum(iConfig.getUntrackedParameter<int>("RecMuonNum", 0)),
  // electrons
  cElPtMin(iConfig.getUntrackedParameter<double>("RecElectronPtMin", 10.)),
  cElEtaMax(iConfig.getUntrackedParameter<double>("RecElectronEtaMax", 2.5)),
  cElHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecElectronHLTriggerMatching")),
  cElNum(iConfig.getUntrackedParameter<int>("RecElectronNum", 0)),
  // taus
  cTauPtMin(iConfig.getUntrackedParameter<double>("RecTauPtMin", 20.)),
  cTauEtaMax(iConfig.getUntrackedParameter<double>("RecTauEtaMax", 2.5)),
  cTauHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecTauHLTriggerMatching")),
  cTauFloatDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecTauFloatDiscriminators")),
  cTauBinaryDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecTauBinaryDiscriminators")),
  cTauNum(iConfig.getUntrackedParameter<int>("RecTauNum", 0)),
  // tracks
  cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.5)),
  cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 2.5)),
  cTrackDxyMax(iConfig.getUntrackedParameter<double>("RecTrackDxyMax", 2.0)),
  cTrackDzMax(iConfig.getUntrackedParameter<double>("RecTrackDzMax", 2.0)),
  cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
  // photons
  cPhotonPtMin(iConfig.getUntrackedParameter<double>("RecPhotonPtMin", 10.)),
  cPhotonEtaMax(iConfig.getUntrackedParameter<double>("RecPhotonEtaMax", 2.5)),
  cPhotonHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecPhotonHLTriggerMatching")),
  cPhotonNum(iConfig.getUntrackedParameter<int>("RecPhotonNum", 0)),
  // jets
  cJetPtMin(iConfig.getUntrackedParameter<double>("RecJetPtMin", 10.)),
  cJetEtaMax(iConfig.getUntrackedParameter<double>("RecJetEtaMax", 5.0)),
  cBtagDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecJetBtagDiscriminators")),
  cJetHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecJetHLTriggerMatching")),
  cJetNum(iConfig.getUntrackedParameter<int>("RecJetNum", 0)),
  // pre-firing weights
  prefweight_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),
  prefweightup_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
  prefweightdown_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
  // collections
  MuonCollectionToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonCollectionTag"))),
  BadGlobalMuonsToken_(consumes<edm::PtrVector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("BadGlobalMuons"))),
  BadDuplicateMuonsToken_(consumes<edm::PtrVector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("BadDuplicateMuons"))),
  ElectronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("ElectronCollectionTag"))),
  applyElectronESShift_(iConfig.getUntrackedParameter<bool>("applyElectronESShift")),

  TauCollectionToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("TauCollectionTag"))),
  TauMVAIsolationRawToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1raw","","TreeProducer"))),
  TauMVAIsolationVLooseToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1VLoose","","TreeProducer"))),
  TauMVAIsolationLooseToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1Loose","","TreeProducer"))),
  TauMVAIsolationMediumToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1Medium","","TreeProducer"))),
  TauMVAIsolationTightToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1Tight","","TreeProducer"))),
  TauMVAIsolationVTightToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1VTight","","TreeProducer"))),
  TauMVAIsolationVVTightToken_(consumes<pat::PATTauDiscriminator>(edm::InputTag("rerunDiscriminationByIsolationMVArun2v1VVTight","","TreeProducer"))),
  JetCollectionToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag"))),
  PuppiJetCollectionToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("PuppiJetCollectionTag"))),
  MetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCollectionTag"))),
  MetCorrCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCorrCollectionTag"))),
  PuppiMetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("PuppiMetCollectionTag"))),
  MvaMetCollectionsTag_(iConfig.getParameter<std::vector<edm::InputTag> >("MvaMetCollectionsTag")),

  GenParticleCollectionToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag"))),
  GenJetCollectionToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("GenJetCollectionTag"))),
  L1ExtraJetCollectionToken_(consumes<l1extra::L1JetParticleCollection>(edm::InputTag("l1extraParticles","Central"))),
  L1ExtraTauCollectionToken_(consumes<l1extra::L1JetParticleCollection>(edm::InputTag("l1extraParticles","Tau"))),
  PackedCantidateCollectionToken_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
  TriggerObjectCollectionToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjectCollectionTag"))),
  BeamSpotToken_(consumes<BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag"))),
  PVToken_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag"))),
  PVwithBSToken_(consumes<RefitVertexCollection>(iConfig.getParameter<edm::InputTag>("PVwithBSCollectionTag"))),
  RefittedPVToken_(consumes<RefitVertexCollection>(iConfig.getParameter<edm::InputTag>("RefittedPVCollectionTag"))),
  RefittedwithBSPVToken_(consumes<RefitVertexCollection>(iConfig.getParameter<edm::InputTag>("RefittedwithBSPVCollectionTag"))),
  LHEToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductTag"))),
  SusyMotherMassToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("SusyMotherMassTag"))),
  SusyLSPMassToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("SusyLSPMassTag"))),
  htxsToken_(consumes<HTXS::HiggsClassification>(iConfig.getParameter<edm::InputTag>("htxsInfo"))),
  sampleName(iConfig.getUntrackedParameter<std::string>("SampleName", "Higgs")),
  propagatorWithMaterial(0)
{
  setTauBranches = true;

  //  propagatorWithMaterial = NULL;
  if(cYear != 2015 && cYear != 2016 && cYear != 2017 && cYear != 2018)
    throw cms::Exception("NTupleMaker") << "Invalid Year : 2015, 2016 2017 and 2018 are allowed!" << endl << "Why is this check even done?! You my dear shall PhD answer that when running on Run 3 data!";
  //if(cPeriod != "Summer11" && cPeriod != "Fall11" && cPeriod != "Summer12" && cPeriod != "PHYS14" && cPeriod != "Spring15" && cPeriod != "Run2015B" && cPeriod != "Run2015C" && cPeriod != "Run2015D")
  //  throw cms::Exception("NTupleMaker") << "Invalid period, only Summer11, Fall11, Summer12, PHYS14, Spring15, Run2015B, Run2015C and Run2015D are allowed!";

  double barrelRadius = 129.;  //p81, p50, ECAL TDR
  double endcapZ      = 320.5; // fig 3.26, p81, ECAL TDR
  Surface::RotationType rot;
  ecalBarrel         = Cylinder::build(Surface::PositionType(0, 0, 0), rot, barrelRadius);
  ecalNegativeEtaEndcap = Plane::build(Surface::PositionType(0, 0, -endcapZ), rot);
  ecalPositiveEtaEndcap = Plane::build(Surface::PositionType(0, 0, endcapZ), rot);

  const char* cmssw_base = getenv("CMSSW_BASE");
  if(!cmssw_base) throw cms::Exception("No CMSSW runtime environment found");
  std::string prefix = std::string(cmssw_base) + "/src/";
  /*
  //recoilCorrector
  if(sampleName.find("ZJets") !=std::string::npos) {
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else if(sampleName.find("WJets")!=std::string::npos){
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_wjets53X_20pv_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else if(sampleName.find("Higgs")!=std::string::npos){
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else
    corrector_ = 0;
  */


  string cmsswBase = (getenv ("CMSSW_BASE"));

  //  if(cgen && !cdata) {
  //    consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  //    consumes<LHEEventProduct>(edm::InputTag("source"));
  //    consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));
  //    consumes<LHERunInfoProduct, edm::InRun>({"source"});
  //  }

  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", cTriggerProcess));
  consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));

  for(std::vector<string>::iterator it = cFlagsProcesses.begin();
      it != cFlagsProcesses.end(); it++){
    consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", it->data()));
  }

  consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumes<vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));

  consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));
  consumes<pat::JetCollection>(edm::InputTag("slimmedJetsPuppi"));

  for(std::vector<edm::InputTag>::iterator mit = MvaMetCollectionsTag_.begin();
      mit != MvaMetCollectionsTag_.end(); mit++){
    MvaMetCollectionsToken_.push_back(consumes<pat::METCollection>(*mit));
  }

  if(cl1objects){
    L1MuonCollectionToken_   = consumes<BXVector<l1t::Muon> >(iConfig.getParameter<edm::InputTag>("L1MuonCollectionTag"));
    L1EGammaCollectionToken_ = consumes<BXVector<l1t::EGamma> >(iConfig.getParameter<edm::InputTag>("L1EGammaCollectionTag"));
    L1TauCollectionToken_    = consumes<BXVector<l1t::Tau> >(iConfig.getParameter<edm::InputTag>("L1TauCollectionTag"));
    L1JetCollectionToken_    = consumes<BXVector<l1t::Jet> >(iConfig.getParameter<edm::InputTag>("L1JetCollectionTag"));

  }

  HLTPrescaleConfig = new HLTPrescaleProvider(iConfig, consumesCollector(), *this);
}

//destructor
NTupleMaker::~NTupleMaker(){
  if(propagatorWithMaterial != 0){ delete propagatorWithMaterial;}

}



void NTupleMaker::beginJob(){
  edm::Service<TFileService> FS;
  FS->file().SetCompressionAlgorithm(ROOT::kLZMA);
  FS->file().SetCompressionLevel(9);
  tree = FS->make<TTree>("AC1B", "AC1B", 1);
  tree->SetMaxVirtualSize(300000000);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);

  tree->Branch("errors", &errors, "errors/i");
  tree->Branch("event_nr", &event_nr, "event_nr/l");
  tree->Branch("event_run", &event_run, "event_run/i");
  tree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
  tree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
  tree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
  tree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
  tree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
  tree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");
  tree->Branch("_passecalBadCalibFilterUpdate",&_passecalBadCalibFilterUpdate,"_passecalBadCalibFilterUpdate/b");
  // beam spot
  if (cbeamspot) {
    tree->Branch("beamspot_x", &beamspot_x, "beamspot_x/F");
    tree->Branch("beamspot_y", &beamspot_y, "beamspot_y/F");
    tree->Branch("beamspot_z", &beamspot_z, "beamspot_z/F");
    tree->Branch("beamspot_xwidth", &beamspot_xwidth, "beamspot_xwidth/F");
    tree->Branch("beamspot_ywidth", &beamspot_ywidth, "beamspot_ywidth/F");
    tree->Branch("beamspot_zsigma", &beamspot_zsigma, "beamspot_zsigma/F");
    tree->Branch("beamspot_cov", &beamspot_cov, "beamspot_cov[6]/F");
  }

  tree->Branch("rho",&rhoNeutral,"rho/F");

  //tree->Branch("track_count", &track_count, "track_count/i");
  //tree->Branch("track_px", track_px, "track_px[track_count]/F");
  //tree->Branch("track_py", track_py, "track_py[track_count]/F");
  //tree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
  //tree->Branch("track_outerx", track_outerx, "track_outerx[track_count]/F");
  //tree->Branch("track_outery", track_outery, "track_outery[track_count]/F");
  //tree->Branch("track_outerz", track_outerz, "track_outerz[track_count]/F");
  //tree->Branch("track_closestpointx", track_closestpointx, "track_closestpointx[track_count]/F");
  //tree->Branch("track_closestpointy", track_closestpointy, "track_closestpointy[track_count]/F");
  //tree->Branch("track_closestpointz", track_closestpointz, "track_closestpointz[track_count]/F");
  //tree->Branch("track_chi2", track_chi2, "track_chi2[track_count]/F");
  //tree->Branch("track_ndof", track_ndof, "track_ndof[track_count]/F");
  //tree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
  //tree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
  //tree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
  //tree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
  //tree->Branch("track_dedxharmonic2", track_dedxharmonic2, "track_dedxharmonic2[track_count]/F");
  //tree->Branch("track_charge", track_charge, "track_charge[track_count]/I");
  //tree->Branch("track_nhits", track_nhits, "track_nhits[track_count]/b");
  //tree->Branch("track_npixelhits", track_npixelhits, "track_npixelhits[track_count]/b");
  //tree->Branch("track_nmissinghits", track_nmissinghits, "track_nmissinghits[track_count]/b");
  //tree->Branch("track_npixellayers", track_npixellayers, "track_npixellayers[track_count]/b");
  //tree->Branch("track_nstriplayers", track_nstriplayers, "track_nstriplayers[track_count]/b");

  // primary vertex
  if (crecprimvertex) {
    tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i");
    tree->Branch("goodprimvertex_count", &goodprimvertex_count, "goodprimvertex_count/i");
    tree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
    tree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
    tree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
    tree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
    tree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
    tree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_pdf/F");
    tree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
    tree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");
    tree->Branch("primvertex_mindz", &primvertex_mindz, "primvertex_mindz/F");
  }
  // primary vertex with BS
  if (crecprimvertexwithbs) {
    tree->Branch("primvertexwithbs_x", &primvertexwithbs_x, "primvertexwithbs_x/F");
    tree->Branch("primvertexwithbs_y", &primvertexwithbs_y, "primvertexwithbs_y/F");
    tree->Branch("primvertexwithbs_z", &primvertexwithbs_z, "primvertexwithbs_z/F");
    tree->Branch("primvertexwithbs_chi2", &primvertexwithbs_chi2, "primvertexwithbs_chi2/F");
    tree->Branch("primvertexwithbs_ndof", &primvertexwithbs_ndof, "primvertexwithbs_ndof/F");
    tree->Branch("primvertexwithbs_ntracks", &primvertexwithbs_ntracks, "primvertexwithbs_ntracks/I");
    tree->Branch("primvertexwithbs_cov", primvertexwithbs_cov, "primvertexwithbs_cov[6]/F");
  }
  // refitted primary vertex
  if(crefittedvertex) {
    tree->Branch("refitvertex_count", &refitvertex_count, "refitvertex_count/i");
    tree->Branch("refitvertex_x", refitvertex_x, "refitvertex_x[refitvertex_count]/F");
    tree->Branch("refitvertex_y", refitvertex_y, "refitvertex_y[refitvertex_count]/F");
    tree->Branch("refitvertex_z", refitvertex_z, "refitvertex_z[refitvertex_count]/F");
    tree->Branch("refitvertex_chi2", refitvertex_chi2, "refitvertex_chi2[refitvertex_count]/F");
    tree->Branch("refitvertex_ndof", refitvertex_ndof, "refitvertex_ndof[refitvertex_count]/F");
    tree->Branch("refitvertex_ntracks", refitvertex_ntracks, "refitvertex_ntracks[refitvertex_count]/I");
    tree->Branch("refitvertex_cov", refitvertex_cov, "refitvertex_cov[refitvertex_count][6]/F");
    tree->Branch("refitvertex_eleIndex", refitvertex_eleIndex, "refitvertex_eleIndex[refitvertex_count][2]/I");
    tree->Branch("refitvertex_muIndex", refitvertex_muIndex, "refitvertex_muIndex[refitvertex_count][2]/I");
    tree->Branch("refitvertex_tauIndex", refitvertex_tauIndex, "refitvertex_tauIndex[refitvertex_count][2]/I");
  }
  // refitted primary vertex with bs
  if(crefittedvertexwithbs) {
    tree->Branch("refitvertexwithbs_count", &refitvertexwithbs_count, "refitvertexwithbs_count/i");
    tree->Branch("refitvertexwithbs_x", refitvertexwithbs_x, "refitvertexwithbs_x[refitvertexwithbs_count]/F");
    tree->Branch("refitvertexwithbs_y", refitvertexwithbs_y, "refitvertexwithbs_y[refitvertexwithbs_count]/F");
    tree->Branch("refitvertexwithbs_z", refitvertexwithbs_z, "refitvertexwithbs_z[refitvertexwithbs_count]/F");
    tree->Branch("refitvertexwithbs_chi2", refitvertexwithbs_chi2, "refitvertexwithbs_chi2[refitvertexwithbs_count]/F");
    tree->Branch("refitvertexwithbs_ndof", refitvertexwithbs_ndof, "refitvertexwithbs_ndof[refitvertexwithbs_count]/F");
    tree->Branch("refitvertexwithbs_ntracks", refitvertexwithbs_ntracks, "refitvertexwithbs_ntracks[refitvertexwithbs_count]/I");
    tree->Branch("refitvertexwithbs_cov", refitvertexwithbs_cov, "refitvertexwithbs_cov[refitvertexwithbs_count][6]/F");
    tree->Branch("refitvertexwithbs_eleIndex", refitvertexwithbs_eleIndex, "refitvertexwithbs_eleIndex[refitvertexwithbs_count][2]/I");
    tree->Branch("refitvertexwithbs_muIndex", refitvertexwithbs_muIndex, "refitvertexwithbs_muIndex[refitvertexwithbs_count][2]/I");
    tree->Branch("refitvertexwithbs_tauIndex", refitvertexwithbs_tauIndex, "refitvertexwithbs_tauIndex[refitvertexwithbs_count][2]/I");
  }

  // muons
  if (crecmuon) {

    tree->Branch("muon_count", &muon_count, "muon_count/i");

    tree->Branch("muon_helixparameters", muon_helixparameters, "muon_helixparameters[muon_count][5]/F");
    tree->Branch("muon_helixparameters_covar", muon_helixparameters_covar,"muon_helixparameters_covar[muon_count][5][5]/F");
    tree->Branch("muon_referencePoint", muon_referencePoint,"muon_referencePoint[muon_count][3]/F");
    tree->Branch("muon_Bfield", muon_Bfield, "muon_Bfield[muon_count]/F");

    tree->Branch("muon_px", muon_px, "muon_px[muon_count]/F");
    tree->Branch("muon_py", muon_py, "muon_py[muon_count]/F");
    tree->Branch("muon_pz", muon_pz, "muon_pz[muon_count]/F");
    tree->Branch("muon_pt", muon_pt, "muon_pt[muon_count]/F");
    tree->Branch("muon_eta", muon_eta, "muon_eta[muon_count]/F");
    tree->Branch("muon_phi", muon_phi, "muon_phi[muon_count]/F");
    tree->Branch("muon_pterror", muon_pterror, "muon_pterror[muon_count]/F");
    tree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_count]/F");
    tree->Branch("muon_normChi2", muon_normChi2, "muon_normChi2[muon_count]/F");
    tree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_count]/F");
    tree->Branch("muon_charge", muon_charge, "muon_charge[muon_count]/F");
    tree->Branch("muon_miniISO", muon_miniISO, "muon_miniISO[muon_count]/F");
    tree->Branch("muon_combQ_chi2LocalPosition", muon_combQ_chi2LocalPosition, "muon_combQ_chi2LocalPosition[muon_count]/F");
    tree->Branch("muon_combQ_trkKink", muon_combQ_trkKink, "muon_combQ_trkKink[muon_count]/F");
    tree->Branch("muon_validFraction", muon_validFraction, "muon_validFraction[muon_count]/F");
    tree->Branch("muon_segmentComp", muon_segmentComp, "muon_segmentComp[muon_count]/F");

    tree->Branch("muon_nMuonStations", muon_nMuonStations,"muon_nMuonStations[muon_count]/i");
    tree->Branch("muon_nMuonHits", muon_nMuonHits,"muon_nMuonHits[muon_count]/i");
    tree->Branch("muon_nPixelHits", muon_nPixelHits,"muon_nPixelHits[muon_count]/i");
    tree->Branch("muon_nTrackerHits", muon_nTrackerHits,"muon_nTrackerHits[muon_count]/i");
    tree->Branch("muon_dxy",muon_dxy,"muon_dxy[muon_count]/F");
    tree->Branch("muon_dxyerr",muon_dxyerr,"muon_dxyerr[muon_count]/F");
    tree->Branch("muon_dz",muon_dz,"muon_dz[muon_count]/F");
    tree->Branch("muon_dzerr",muon_dzerr,"muon_dzerr[muon_count]/F");
    tree->Branch("muon_vx",muon_vx,"muon_vx[muon_count]/F");
    tree->Branch("muon_vy",muon_vy,"muon_vy[muon_count]/F");
    tree->Branch("muon_vz",muon_vz,"muon_vz[muon_count]/F");
    tree->Branch("muon_chargedHadIso",muon_chargedHadIso,"muon_chargedHadIso[muon_count]/F");
    tree->Branch("muon_neutralHadIso",muon_neutralHadIso,"muon_neutralHadIso[muon_count]/F");
    tree->Branch("muon_photonIso",muon_photonIso,"muon_photonIso[muon_count]/F");
    tree->Branch("muon_puIso",muon_puIso,"muon_puIso[muon_count]/F");

    tree->Branch("muon_r03_sumChargedHadronPt",muon_r03_sumChargedHadronPt,"muon_r03_sumChargedHadronPt[muon_count]/F");
    tree->Branch("muon_r03_sumChargedParticlePt",muon_r03_sumChargedParticlePt,"muon_r03_sumChargedParticlePt[muon_count]/F");
    tree->Branch("muon_r03_sumNeutralHadronEt",muon_r03_sumNeutralHadronEt,"muon_r03_sumNeutralHadronEt[muon_count]/F");
    tree->Branch("muon_r03_sumPhotonEt",muon_r03_sumPhotonEt,"muon_r03_sumPhotonEt[muon_count]/F");
    tree->Branch("muon_r03_sumNeutralHadronEtHighThreshold",muon_r03_sumNeutralHadronEtHighThreshold,"muon_r03_sumNeutralHadronEtHighThreshold[muon_count]/F");
    tree->Branch("muon_r03_sumPhotonEtHighThreshold",muon_r03_sumPhotonEtHighThreshold,"muon_r03_sumPhotonEtHighThreshold[muon_count]/F");
    tree->Branch("muon_r03_sumPUPt",muon_r03_sumPUPt,"muon_r03_sumPUPt[muon_count]/F");

    tree->Branch("muon_r04_sumChargedHadronPt",muon_r04_sumChargedHadronPt,"muon_r04_sumChargedHadronPt[muon_count]/F");
    tree->Branch("muon_r04_sumChargedParticlePt",muon_r04_sumChargedParticlePt,"muon_r04_sumChargedParticlePt[muon_count]/F");
    tree->Branch("muon_r04_sumNeutralHadronEt",muon_r04_sumNeutralHadronEt,"muon_r04_sumNeutralHadronEt[muon_count]/F");
    tree->Branch("muon_r04_sumPhotonEt",muon_r04_sumPhotonEt,"muon_r04_sumPhotonEt[muon_count]/F");
    tree->Branch("muon_r04_sumNeutralHadronEtHighThreshold",muon_r04_sumNeutralHadronEtHighThreshold,"muon_r04_sumNeutralHadronEtHighThreshold[muon_count]/F");
    tree->Branch("muon_r04_sumPhotonEtHighThreshold",muon_r04_sumPhotonEtHighThreshold,"muon_r04_sumPhotonEtHighThreshold[muon_count]/F");
    tree->Branch("muon_r04_sumPUPt",muon_r04_sumPUPt,"muon_r04_sumPUPt[muon_count]/F");

    tree->Branch("muon_isPF",muon_isPF,"muon_isPF[muon_count]/O");
    tree->Branch("muon_isGlobal",muon_isGlobal,"muon_isGlobal[muon_count]/O");
    tree->Branch("muon_isTracker",muon_isTracker,"muon_isTracker[muon_count]/O");
    tree->Branch("muon_isTight",muon_isTight,"muon_isTight[muon_count]/O");
    tree->Branch("muon_isLoose",muon_isLoose,"muon_isLoose[muon_count]/O");
    tree->Branch("muon_isMedium",muon_isMedium,"muon_isMedium[muon_count]/O");
    tree->Branch("muon_isICHEP",muon_isICHEP,"muon_isICHEP[muon_count]/O");
    tree->Branch("muon_genmatch", muon_genmatch, "muon_genmatch[muon_count]/I");
    tree->Branch("muon_isDuplicate",muon_isDuplicate,"muon_isDuplicate[muon_count]/O");
    tree->Branch("muon_isBad",muon_isBad,"muon_isBad[muon_count]/O");

    /*tree->Branch("dimuon_count", &dimuon_count, "dimuon_count/i");
    tree->Branch("dimuon_leading", dimuon_leading, "dimuon_leading[dimuon_count]/i");
    tree->Branch("dimuon_trailing", dimuon_trailing, "dimuon_trailing[dimuon_count]/i");
    tree->Branch("dimuon_dist2D", dimuon_dist2D, "dimuon_dist2D[dimuon_count]/F");
    tree->Branch("dimuon_dist2DE", dimuon_dist2DE, "dimuon_dist2DE[dimuon_count]/F");
    tree->Branch("dimuon_dist3D", dimuon_dist3D, "dimuon_dist3D[dimuon_count]/F");
    tree->Branch("dimuon_dist3DE", dimuon_dist3DE, "dimuon_dist3DE[dimuon_count]/F");
    */

}

  // pf jets
  if (crecpfjet) {
    tree->Branch("pfjet_count", &pfjet_count, "pfjet_count/i");
    tree->Branch("pfjet_e", pfjet_e, "pfjet_e[pfjet_count]/F");
    tree->Branch("pfjet_px", pfjet_px, "pfjet_px[pfjet_count]/F");
    tree->Branch("pfjet_py", pfjet_py, "pfjet_py[pfjet_count]/F");
    tree->Branch("pfjet_pz", pfjet_pz, "pfjet_pz[pfjet_count]/F");
    tree->Branch("pfjet_pt", pfjet_pt, "pfjet_pt[pfjet_count]/F");
    tree->Branch("pfjet_eta", pfjet_eta, "pfjet_eta[pfjet_count]/F");
    tree->Branch("pfjet_phi", pfjet_phi, "pfjet_phi[pfjet_count]/F");
    tree->Branch("pfjet_neutralhadronicenergy", pfjet_neutralhadronicenergy, "pfjet_neutralhadronicenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedhadronicenergy", pfjet_chargedhadronicenergy, "pfjet_chargedhadronicenergy[pfjet_count]/F");
    tree->Branch("pfjet_neutralemenergy", pfjet_neutralemenergy, "pfjet_neutralemenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedemenergy", pfjet_chargedemenergy, "pfjet_chargedemenergy[pfjet_count]/F");
    tree->Branch("pfjet_muonenergy", pfjet_muonenergy, "pfjet_muonenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedmuonenergy", pfjet_chargedmuonenergy, "pfjet_chargedmuonenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedmulti", pfjet_chargedmulti, "pfjet_chargedmulti[pfjet_count]/i");
    tree->Branch("pfjet_neutralmulti", pfjet_neutralmulti, "pfjet_neutralmulti[pfjet_count]/i");
    tree->Branch("pfjet_chargedhadronmulti", pfjet_chargedhadronmulti, "pfjet_chargedhadronmulti[pfjet_count]/i");
    tree->Branch("pfjet_energycorr", pfjet_energycorr, "pfjet_energycorr[pfjet_count]/F");
    tree->Branch("pfjet_energycorr_l1fastjet", pfjet_energycorr_l1fastjet, "pfjet_energycorr_l1fastjet[pfjet_count]/F");
    tree->Branch("pfjet_energycorr_l2relative", pfjet_energycorr_l2relative, "pfjet_energycorr_l2relative[pfjet_count]/F");
    tree->Branch("pfjet_energycorr_l3absolute", pfjet_energycorr_l3absolute, "pfjet_energycorr_l3absolute[pfjet_count]/F");
    tree->Branch("pfjet_energycorr_l2l3residual", pfjet_energycorr_l2l3residual, "pfjet_energycorr_l2l3residual[pfjet_count]/F");
    tree->Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour[pfjet_count]/I");
    tree->Branch("pfjet_btag", pfjet_btag,"pfjet_btag[pfjet_count][10]/F");
    tree->Branch("pfjet_jecUncertainty",pfjet_jecUncertainty,"pfjet_jecUncertainty[pfjet_count]/F");
    tree->Branch("pfjet_pu_jet_fullId_loose", pfjet_pu_jet_fullId_loose, "pfjet_pu_jet_fullId_loose[pfjet_count]/O");
    tree->Branch("pfjet_pu_jet_fullId_medium", pfjet_pu_jet_fullId_medium, "pfjet_pu_jet_fullId_medium[pfjet_count]/O");
    tree->Branch("pfjet_pu_jet_fullId_tight", pfjet_pu_jet_fullId_tight, "pfjet_pu_jet_fullId_tight[pfjet_count]/O");
    tree->Branch("pfjet_pu_jet_fullDisc_mva", pfjet_pu_jet_fullDisc_mva, "pfjet_pu_jet_fullDisc_mva[pfjet_count]/F");
  }

  if (crecpfpuppijet) {
    tree->Branch("pfjetpuppi_count", &pfjetpuppi_count, "pfjetpuppi_count/i");
    tree->Branch("pfjetpuppi_e", pfjetpuppi_e, "pfjetpuppi_e[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_px", pfjetpuppi_px, "pfjetpuppi_px[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_py", pfjetpuppi_py, "pfjetpuppi_py[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_pz", pfjetpuppi_pz, "pfjetpuppi_pz[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_pt", pfjetpuppi_pt, "pfjetpuppi_pt[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_eta", pfjetpuppi_eta, "pfjetpuppi_eta[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_phi", pfjetpuppi_phi, "pfjetpuppi_phi[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_neutralhadronicenergy", pfjetpuppi_neutralhadronicenergy, "pfjetpuppi_neutralhadronicenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_chargedhadronicenergy", pfjetpuppi_chargedhadronicenergy, "pfjetpuppi_chargedhadronicenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_neutralemenergy", pfjetpuppi_neutralemenergy, "pfjetpuppi_neutralemenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_chargedemenergy", pfjetpuppi_chargedemenergy, "pfjetpuppi_chargedemenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_muonenergy", pfjetpuppi_muonenergy, "pfjetpuppi_muonenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_chargedmuonenergy", pfjetpuppi_chargedmuonenergy, "pfjetpuppi_chargedmuonenergy[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_chargedmulti", pfjetpuppi_chargedmulti, "pfjetpuppi_chargedmulti[pfjetpuppi_count]/i");
    tree->Branch("pfjetpuppi_neutralmulti", pfjetpuppi_neutralmulti, "pfjetpuppi_neutralmulti[pfjetpuppi_count]/i");
    tree->Branch("pfjetpuppi_chargedhadronmulti", pfjetpuppi_chargedhadronmulti, "pfjetpuppi_chargedhadronmulti[pfjetpuppi_count]/i");
    tree->Branch("pfjetpuppi_energycorr", pfjetpuppi_energycorr, "pfjetpuppi_energycorr[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_energycorr_l1fastjet", pfjetpuppi_energycorr_l1fastjet, "pfjetpuppi_energycorr_l1fastjet[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_energycorr_l2relative", pfjetpuppi_energycorr_l2relative, "pfjetpuppi_energycorr_l2relative[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_energycorr_l3absolute", pfjetpuppi_energycorr_l3absolute, "pfjetpuppi_energycorr_l3absolute[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_energycorr_l2l3residual", pfjetpuppi_energycorr_l2l3residual, "pfjetpuppi_energycorr_l2l3residual[pfjetpuppi_count]/F");
    tree->Branch("pfjetpuppi_flavour", pfjetpuppi_flavour, "pfjetpuppi_flavour[pfjetpuppi_count]/I");
    tree->Branch("pfjetpuppi_btag", pfjetpuppi_btag,"pfjetpuppi_btag[pfjetpuppi_count][10]/F");
    tree->Branch("pfjetpuppi_jecUncertainty",pfjetpuppi_jecUncertainty,"pfjetpuppi_jecUncertainty[pfjetpuppi_count]/F");
  }
  // electrons
  if (crecelectron) {
    tree->Branch("electron_count", &electron_count, "electron_count/i");
    tree->Branch("electron_px", electron_px, "electron_px[electron_count]/F");
    tree->Branch("electron_py", electron_py, "electron_py[electron_count]/F");
    tree->Branch("electron_pz", electron_pz, "electron_pz[electron_count]/F");
    tree->Branch("electron_pt", electron_pt, "electron_pt[electron_count]/F");
    tree->Branch("electron_eta", electron_eta, "electron_eta[electron_count]/F");
    tree->Branch("electron_phi", electron_phi, "electron_phi[electron_count]/F");
    tree->Branch("electron_px_energyscale_up", electron_px_energyscale_up, "electron_px_energyscale_up[electron_count]/F");
    tree->Branch("electron_px_energyscale_down", electron_px_energyscale_down, "electron_px_energyscale_down[electron_count]/F");
    tree->Branch("electron_py_energyscale_up", electron_py_energyscale_up, "electron_py_energyscale_up[electron_count]/F");
    tree->Branch("electron_py_energyscale_down", electron_py_energyscale_down, "electron_py_energyscale_down[electron_count]/F");
    tree->Branch("electron_pz_energyscale_up", electron_pz_energyscale_up, "electron_pz_energyscale_up[electron_count]/F");
    tree->Branch("electron_pz_energyscale_down", electron_pz_energyscale_down, "electron_pz_energyscale_down[electron_count]/F");
    tree->Branch("electron_pt_energyscale_up", electron_pt_energyscale_up, "electron_pt_energyscale_up[electron_count]/F");
    tree->Branch("electron_pt_energyscale_down", electron_pt_energyscale_down, "electron_pt_energyscale_down[electron_count]/F");
    tree->Branch("electron_px_energysigma_up", electron_px_energysigma_up, "electron_px_energysigma_up[electron_count]/F");
    tree->Branch("electron_px_energysigma_down", electron_px_energysigma_down, "electron_px_energysigma_down[electron_count]/F");
    tree->Branch("electron_py_energysigma_up", electron_py_energysigma_up, "electron_py_energysigma_up[electron_count]/F");
    tree->Branch("electron_py_energysigma_down", electron_py_energysigma_down, "electron_py_energysigma_down[electron_count]/F");
    tree->Branch("electron_pz_energysigma_up", electron_pz_energysigma_up, "electron_pz_energysigma_up[electron_count]/F");
    tree->Branch("electron_pz_energysigma_down", electron_pz_energysigma_down, "electron_pz_energysigma_down[electron_count]/F");
    tree->Branch("electron_pt_energysigma_up", electron_pt_energysigma_up, "electron_pt_energysigma_up[electron_count]/F");
    tree->Branch("electron_pt_energysigma_down", electron_pt_energysigma_down, "electron_pt_energysigma_down[electron_count]/F");
    tree->Branch("electron_trackchi2", electron_trackchi2, "electron_trackchi2[electron_count]/F");
    tree->Branch("electron_trackndof", electron_trackndof, "electron_trackndof[electron_count]/F");
    tree->Branch("electron_outerx", electron_outerx, "electron_outerx[electron_count]/F");
    tree->Branch("electron_outery", electron_outery, "electron_outery[electron_count]/F");
    tree->Branch("electron_outerz", electron_outerz, "electron_outerz[electron_count]/F");
    tree->Branch("electron_vx", electron_vx, "electron_vx[electron_count]/F");
    tree->Branch("electron_vy", electron_vy, "electron_vy[electron_count]/F");
    tree->Branch("electron_vz", electron_vz, "electron_vz[electron_count]/F");
    tree->Branch("electron_esuperclusterovertrack", electron_esuperclusterovertrack, "electron_esuperclusterovertrack[electron_count]/F");
    tree->Branch("electron_eseedclusterovertrack", electron_eseedclusterovertrack, "electron_eseedclusterovertrack[electron_count]/F");
    tree->Branch("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack, "electron_deltaetasuperclustertrack[electron_count]/F");
    tree->Branch("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack, "electron_deltaphisuperclustertrack[electron_count]/F");
    tree->Branch("electron_e1x5", electron_e1x5, "electron_e1x5[electron_count]/F");
    tree->Branch("electron_e2x5", electron_e2x5, "electron_e2x5[electron_count]/F");
    tree->Branch("electron_e5x5", electron_e5x5, "electron_e5x5[electron_count]/F");
    tree->Branch("electron_sigmaetaeta", electron_sigmaetaeta, "electron_sigmaetaeta[electron_count]/F");
    tree->Branch("electron_sigmaietaieta", electron_sigmaietaieta, "electron_sigmaietaieta[electron_count]/F");
    tree->Branch("electron_ehcaloverecal", electron_ehcaloverecal, "electron_ehcaloverecal[electron_count]/F");
    tree->Branch("electron_ehcaloverecaldepth1", electron_ehcaloverecaldepth1, "electron_ehcaloverecaldepth1[electron_count]/F");
    tree->Branch("electron_ehcaloverecaldepth2", electron_ehcaloverecaldepth2, "electron_ehcaloverecaldepth2[electron_count]/F");
    tree->Branch("electron_full5x5_sigmaietaieta", electron_full5x5_sigmaietaieta, "electron_full5x5_sigmaietaieta[electron_count]/F");
    tree->Branch("electron_ooemoop", electron_ooemoop, "electron_ooemoop[electron_count]/F");
    tree->Branch("electron_miniISO", electron_miniISO, "electron_miniISO[electron_count]/F");

    tree->Branch("electron_superclusterEta", electron_superClusterEta, "electron_superclusterEta[electron_count]/F");
    tree->Branch("electron_superclusterPhi", electron_superClusterPhi, "electron_superclusterPhi[electron_count]/F");
    tree->Branch("electron_superclusterX", electron_superClusterX, "electron_superclusterX[electron_count]/F");
    tree->Branch("electron_superclusterY", electron_superClusterY, "electron_superclusterY[electron_count]/F");
    tree->Branch("electron_superclusterZ", electron_superClusterZ, "electron_superclusterZ[electron_count]/F");

    tree->Branch("electron_detaInSeed", electron_detaInSeed, "electron_detaInSeed[electron_count]/F");
    tree->Branch("electron_he", electron_he, "electron_he[electron_count]/F");
    tree->Branch("electron_eaIsolation", electron_eaIsolation, "electron_eaIsolation[electron_count]/F");

    tree->Branch("electron_chargedHadIso", electron_chargedHadIso,"electron_chargedHadIso[electron_count]/F");
    tree->Branch("electron_neutralHadIso", electron_neutralHadIso,"electron_neutralHadIso[electron_count]/F");
    tree->Branch("electron_photonIso",     electron_photonIso,    "electron_photonIso[electron_count]/F");
    tree->Branch("electron_puIso",         electron_puIso,        "electron_puIso[electron_count]/F");

    tree->Branch("electron_r03_sumChargedHadronPt",electron_r03_sumChargedHadronPt,"electron_r03_sumChargedHadronPt[electron_count]/F");
    tree->Branch("electron_r03_sumChargedParticlePt",electron_r03_sumChargedParticlePt,"electron_r03_sumChargedParticlePt[electron_count]/F");
    tree->Branch("electron_r03_sumNeutralHadronEt",electron_r03_sumNeutralHadronEt,"electron_r03_sumNeutralHadronEt[electron_count]/F");
    tree->Branch("electron_r03_sumPhotonEt",electron_r03_sumPhotonEt,"electron_r03_sumPhotonEt[electron_count]/F");
    tree->Branch("electron_r03_sumNeutralHadronEtHighThreshold",electron_r03_sumNeutralHadronEtHighThreshold,"electron_r03_sumNeutralHadronEtHighThreshold[electron_count]/F");
    tree->Branch("electron_r03_sumPhotonEtHighThreshold",electron_r03_sumPhotonEtHighThreshold,"electron_r03_sumPhotonEtHighThreshold[electron_count]/F");
    tree->Branch("electron_r03_sumPUPt",electron_r03_sumPUPt,"electron_r03_sumPUPt[electron_count]/F");

    tree->Branch("electron_nhits", electron_nhits, "electron_nhits[electron_count]/b");
    tree->Branch("electron_npixelhits", electron_npixelhits, "electron_npixelhits[electron_count]/b");
    tree->Branch("electron_nmissinghits", electron_nmissinghits, "electron_nmissinghits[electron_count]/b");
    tree->Branch("electron_nmissinginnerhits", electron_nmissinginnerhits, "electron_nmissinginnerhits[electron_count]/b");
    tree->Branch("electron_npixellayers", electron_npixellayers, "electron_npixellayers[electron_count]/b");
    tree->Branch("electron_nstriplayers", electron_nstriplayers, "electron_nstriplayers[electron_count]/b");
    tree->Branch("electron_dxy", electron_dxy, "electron_dxy[electron_count]/F");
    tree->Branch("electron_dxyerr", electron_dxyerr, "electron_dxyerr[electron_count]/F");
    tree->Branch("electron_dz", electron_dz, "electron_dz[electron_count]/F");
    tree->Branch("electron_dzerr", electron_dzerr, "electron_dzerr[electron_count]/F");
    tree->Branch("electron_convdist", electron_convdist, "electron_convdist[electron_count]/F");

    tree->Branch("electron_gapinfo", electron_gapinfo, "electron_gapinfo[electron_count]/i");
    tree->Branch("electron_chargeinfo", electron_chargeinfo, "electron_chargeinfo[electron_count]/i");
    tree->Branch("electron_fbrems", electron_fbrems, "electron_fbrems[electron_count]/F");
    tree->Branch("electron_numbrems", electron_numbrems, "electron_numbrems[electron_count]/I");
    tree->Branch("electron_charge", electron_charge, "electron_charge[electron_count]/F");
    tree->Branch("electron_superclusterindex", electron_superclusterindex, "electron_superclusterindex[electron_count]/I");
    tree->Branch("electron_info", electron_info, "electron_info[electron_count]/b");

    tree->Branch("electron_cutId_veto_Summer16", electron_cutId_veto_Summer16, "electron_cutId_veto_Summer16[electron_count]/O");
    //tree->Branch("electron_cutId_loose_Summer16", electron_cutId_loose_Summer16, "electron_cutId_loose_Summer16[electron_count]/O");
    //tree->Branch("electron_cutId_medium_Summer16", electron_cutId_medium_Summer16, "electron_cutId_medium_Summer16[electron_count]/O");
    //tree->Branch("electron_cutId_tight_Summer16", electron_cutId_tight_Summer16, "electron_cutId_tight_Summer16[electron_count]/O");
    tree->Branch("electron_mva_value_Spring16_v1", electron_mva_value_Spring16_v1, "electron_mva_value_Spring16_v1[electron_count]/F");
    tree->Branch("electron_mva_wp90_general_Spring16_v1", electron_mva_wp90_general_Spring16_v1, "electron_mva_wp90_general_Spring16_v1[electron_count]/F");
    tree->Branch("electron_mva_wp80_general_Spring16_v1", electron_mva_wp80_general_Spring16_v1, "electron_mva_wp80_general_Spring16_v1[electron_count]/F");

    //new for 9.4.0
    tree->Branch("electron_mva_value_Iso_Fall17_v1", electron_mva_value_Iso_Fall17_v1, "electron_mva_value_Iso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_value_noIso_Fall17_v1", electron_mva_value_noIso_Fall17_v1, "electron_mva_value_noIso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_wp90_Iso_Fall17_v1", electron_mva_wp90_Iso_Fall17_v1, "electron_mva_wp90_Iso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_wp80_Iso_Fall17_v1", electron_mva_wp80_Iso_Fall17_v1, "electron_mva_wp80_Iso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_Loose_Iso_Fall17_v1", electron_mva_Loose_Iso_Fall17_v1, "electron_mva_Loose_Iso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_wp90_noIso_Fall17_v1", electron_mva_wp90_noIso_Fall17_v1, "electron_mva_wp90_noIso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_wp80_noIso_Fall17_v1", electron_mva_wp80_noIso_Fall17_v1, "electron_mva_wp80_noIso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_Loose_noIso_Fall17_v1", electron_mva_Loose_noIso_Fall17_v1, "electron_mva_Loose_noIso_Fall17_v1[electron_count]/F");
    tree->Branch("electron_mva_value_Iso_Fall17_v2", electron_mva_value_Iso_Fall17_v2, "electron_mva_value_Iso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_value_noIso_Fall17_v2", electron_mva_value_noIso_Fall17_v2, "electron_mva_value_noIso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_wp90_Iso_Fall17_v2", electron_mva_wp90_Iso_Fall17_v2, "electron_mva_wp90_Iso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_wp80_Iso_Fall17_v2", electron_mva_wp80_Iso_Fall17_v2, "electron_mva_wp80_Iso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_Loose_Iso_Fall17_v2", electron_mva_Loose_Iso_Fall17_v2, "electron_mva_Loose_Iso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_wp90_noIso_Fall17_v2", electron_mva_wp90_noIso_Fall17_v2, "electron_mva_wp90_noIso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_wp80_noIso_Fall17_v2", electron_mva_wp80_noIso_Fall17_v2, "electron_mva_wp80_noIso_Fall17_v2[electron_count]/F");
    tree->Branch("electron_mva_Loose_noIso_Fall17_v2", electron_mva_Loose_noIso_Fall17_v2, "electron_mva_Loose_noIso_Fall17_v2[electron_count]/F");
    //cut-based Fall17
    tree->Branch("electron_cutId_veto_Fall17", electron_cutId_veto_Fall17, "electron_cutId_veto_Fall17[electron_count]/O");
    //tree->Branch("electron_cutId_loose_Fall17", electron_cutId_loose_Fall17, "electron_cutId_loose_Fall17[electron_count]/O");
    //tree->Branch("electron_cutId_medium_Fall17", electron_cutId_medium_Fall17, "electron_cutId_medium_Fall17[electron_count]/O");
    //tree->Branch("electron_cutId_tight_Fall17", electron_cutId_tight_Fall17, "electron_cutId_tight_Fall17[electron_count]/O");
    //cut-based Fall17V2
    tree->Branch("electron_cutId_veto_Fall17V2", electron_cutId_veto_Fall17V2, "electron_cutId_veto_Fall17V2[electron_count]/O");
    //tree->Branch("electron_cutId_loose_Fall17V2", electron_cutId_loose_Fall17V2, "electron_cutId_loose_Fall17V2[electron_count]/O");
    //tree->Branch("electron_cutId_medium_Fall17V2", electron_cutId_medium_Fall17V2, "electron_cutId_medium_Fall17V2[electron_count]/O");
    //tree->Branch("electron_cutId_tight_Fall17V2", electron_cutId_tight_Fall17V2, "electron_cutId_tight_Fall17V2[electron_count]/O");
    //end of new

    tree->Branch("electron_pass_conversion", electron_pass_conversion, "electron_pass_conversion[electron_count]/O");

    tree->Branch("electron_genmatch", electron_genmatch, "electron_genmatch[electron_count]/I");

  }

  //tree->Branch("photon_count", &photon_count, "photon_count/i");
  //tree->Branch("photon_px", photon_px, "photon_px[photon_count]/F");
  //tree->Branch("photon_py", photon_py, "photon_py[photon_count]/F");
  //tree->Branch("photon_pz", photon_pz, "photon_pz[photon_count]/F");
  //tree->Branch("photon_e1x5", photon_e1x5, "photon_e1x5[photon_count]/F");
  //tree->Branch("photon_e2x5", photon_e2x5, "photon_e2x5[photon_count]/F");
  //tree->Branch("photon_e3x3", photon_e3x3, "photon_e3x3[photon_count]/F");
  //tree->Branch("photon_e5x5", photon_e5x5, "photon_e5x5[photon_count]/F");
  //tree->Branch("photon_sigmaetaeta", photon_sigmaetaeta, "photon_sigmaetaeta[photon_count]/F");
  //tree->Branch("photon_sigmaietaieta", photon_sigmaietaieta, "photon_sigmaietaieta[photon_count]/F");
  //tree->Branch("photon_ehcaloverecal", photon_ehcaloverecal, "photon_ehcaloverecal[photon_count]/F");
  //tree->Branch("photon_ehcaloverecaldepth1", photon_ehcaloverecaldepth1, "photon_ehcaloverecaldepth1[photon_count]/F");
  //tree->Branch("photon_ehcaloverecaldepth2", photon_ehcaloverecaldepth2, "photon_ehcaloverecaldepth2[photon_count]/F");
  //tree->Branch("photon_maxenergyxtal", photon_maxenergyxtal, "photon_maxenergyxtal[photon_count]/F");
  //tree->Branch("photon_isolationr3track", photon_isolationr3track, "photon_isolationr3track[photon_count]/F");
  //tree->Branch("photon_isolationr3trackhollow", photon_isolationr3trackhollow, "photon_isolationr3trackhollow[photon_count]/F");
  //tree->Branch("photon_isolationr3ntrack", photon_isolationr3ntrack, "photon_isolationr3ntrack[photon_count]/i");
  //tree->Branch("photon_isolationr3ntrackhollow", photon_isolationr3ntrackhollow, "photon_isolationr3ntrackhollow[photon_count]/i");
  //tree->Branch("photon_isolationr3ecal", photon_isolationr3ecal, "photon_isolationr3ecal[photon_count]/F");
  //tree->Branch("photon_isolationr3hcal", photon_isolationr3hcal, "photon_isolationr3hcal[photon_count]/F");
  //tree->Branch("photon_isolationr4track", photon_isolationr4track, "photon_isolationr4track[photon_count]/F");
  //tree->Branch("photon_isolationr4trackhollow", photon_isolationr4trackhollow, "photon_isolationr4trackhollow[photon_count]/F");
  //tree->Branch("photon_isolationr4ntrack", photon_isolationr4ntrack, "photon_isolationr4ntrack[photon_count]/i");
  //tree->Branch("photon_isolationr4ntrackhollow", photon_isolationr4ntrackhollow, "photon_isolationr4ntrackhollow[photon_count]/i");
  //tree->Branch("photon_isolationr4ecal", photon_isolationr4ecal, "photon_isolationr4ecal[photon_count]/F");
  //tree->Branch("photon_isolationr4hcal", photon_isolationr4hcal, "photon_isolationr4hcal[photon_count]/F");
  //tree->Branch("photon_superclusterindex", photon_superclusterindex, "photon_superclusterindex[photon_count]/I");
  //tree->Branch("photon_info", photon_info, "photon_info[photon_count]/b");
  //tree->Branch("photon_gapinfo", photon_gapinfo, "photon_gapinfo[photon_count]/i");
  //tree->Branch("photon_trigger", photon_trigger, "photon_trigger[photon_count]/i");
  //tree->Branch("photon_conversionbegin", photon_conversionbegin, "photon_conversionbegin[photon_count]/i");

  // taus
  if (crectau) {
 
    tree->Branch("tau_count", &tau_count, "tau_count/i");

    tree->Branch("tau_helixparameters", tau_helixparameters, "tau_helixparameters[tau_count][5]/F");
    tree->Branch("tau_helixparameters_covar", tau_helixparameters_covar,"tau_helixparameters_covar[tau_count][5][5]/F");
    tree->Branch("tau_referencePoint", tau_referencePoint,"tau_referencePoint[tau_count][3]/F");
    tree->Branch("tau_Bfield", tau_Bfield, "tau_Bfield[tau_count]/F");

    tree->Branch("tau_e", tau_e, "tau_e[tau_count]/F");
    tree->Branch("tau_px", tau_px, "tau_px[tau_count]/F");
    tree->Branch("tau_py", tau_py, "tau_py[tau_count]/F");
    tree->Branch("tau_pz", tau_pz, "tau_pz[tau_count]/F");
    tree->Branch("tau_mass", tau_mass, "tau_mass[tau_count]/F");
    tree->Branch("tau_eta", tau_eta, "tau_eta[tau_count]/F");
    tree->Branch("tau_phi", tau_phi, "tau_phi[tau_count]/F");
    tree->Branch("tau_pt", tau_pt, "tau_pt[tau_count]/F");

    tree->Branch("tau_vertexx", tau_vertexx, "tau_vertexx[tau_count]/F");
    tree->Branch("tau_vertexy", tau_vertexy, "tau_vertexy[tau_count]/F");
    tree->Branch("tau_vertexz", tau_vertexz, "tau_vertexz[tau_count]/F");

    tree->Branch("tau_pca2D_x", tau_pca2D_x, "tau_pca2D_x[tau_count]/F");
    tree->Branch("tau_pca2D_y", tau_pca2D_y, "tau_pca2D_y[tau_count]/F");
    tree->Branch("tau_pca2D_z", tau_pca2D_z, "tau_pca2D_z[tau_count]/F");

    tree->Branch("tau_pca3D_x", tau_pca3D_x, "tau_pca3D_x[tau_count]/F");
    tree->Branch("tau_pca3D_y", tau_pca3D_y, "tau_pca3D_y[tau_count]/F");
    tree->Branch("tau_pca3D_z", tau_pca3D_z, "tau_pca3D_z[tau_count]/F");


    tree->Branch("tau_dxy", tau_dxy, "tau_dxy[tau_count]/F");
    tree->Branch("tau_dxySig", tau_dxySig, "tau_dxySig[tau_count]/F");
    tree->Branch("tau_dz", tau_dz, "tau_dz[tau_count]/F");
    tree->Branch("tau_ip3d", tau_ip3d, "tau_ip3d[tau_count]/F");
    tree->Branch("tau_ip3dSig", tau_ip3dSig, "tau_ip3dSig[tau_count]/F");
    tree->Branch("tau_charge", tau_charge, "tau_charge[tau_count]/F");

    tree->Branch("tau_flightLength", tau_flightLength, "tau_flightLength[tau_count]/F");
    tree->Branch("tau_flightLengthSig", tau_flightLengthSig, "tau_flightLengthSig[tau_count]/F");
    tree->Branch("tau_SV_x", tau_SV_x, "tau_SV_x[tau_count]/F");
    tree->Branch("tau_SV_y", tau_SV_y, "tau_SV_y[tau_count]/F");
    tree->Branch("tau_SV_z", tau_SV_z, "tau_SV_z[tau_count]/F");
    tree->Branch("tau_SV_cov", tau_SV_cov, "tau_SV_cov[tau_count][6]/F");

    tree->Branch("tau_genjet_px", tau_genjet_px, "tau_genjet_px[tau_count]/F");
    tree->Branch("tau_genjet_py", tau_genjet_py, "tau_genjet_py[tau_count]/F");
    tree->Branch("tau_genjet_pz", tau_genjet_pz, "tau_genjet_pz[tau_count]/F");
    tree->Branch("tau_genjet_e", tau_genjet_e, "tau_genjet_e[tau_count]/F");
    tree->Branch("tau_genmatch", tau_genmatch, "tau_genmatch[tau_count]/I");

    tree->Branch("tau_leadchargedhadrcand_px",  tau_leadchargedhadrcand_px,  "tau_leadchargedhadrcand_px[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_py",  tau_leadchargedhadrcand_py,  "tau_leadchargedhadrcand_py[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_pz",  tau_leadchargedhadrcand_pz,  "tau_leadchargedhadrcand_pz[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_mass",tau_leadchargedhadrcand_mass,"tau_leadchargedhadrcand_mass[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_id",  tau_leadchargedhadrcand_id,  "tau_leadchargedhadrcand_id[tau_count]/I");
    tree->Branch("tau_leadchargedhadrcand_dxy", tau_leadchargedhadrcand_dxy, "tau_leadchargedhadrcand_dxy[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_dz",  tau_leadchargedhadrcand_dz,  "tau_leadchargedhadrcand_dz[tau_count]/F");
    tree->Branch("tau_leadchargedhadrcand_lostPixelHits", tau_leadchargedhadrcand_lostPixelHits, "tau_leadchargedhadrcand_lostPixelHits[tau_count]/I");
    tree->Branch("tau_leadchargedhadrcand_pvAssocQ", tau_leadchargedhadrcand_pvAssocQ, "tau_leadchargedhadrcand_pvAssocQ[tau_count]/I");

    tree->Branch("tau_L1trigger_match", tau_L1trigger_match, "tau_L1trigger_match[tau_count]/O");

    tree->Branch("tau_signalChargedHadrCands_size", tau_signalChargedHadrCands_size, "tau_signalChargedHadrCands_size[tau_count]/i");
    tree->Branch("tau_signalNeutralHadrCands_size", tau_signalNeutralHadrCands_size, "tau_signalNeutralHadrCands_size[tau_count]/i");
    tree->Branch("tau_signalGammaCands_size", tau_signalGammaCands_size, "tau_signalGammaCands_size[tau_count]/i");

    tree->Branch("tau_isolationChargedHadrCands_size", tau_isolationChargedHadrCands_size, "tau_isolationChargedHadrCands_size[tau_count]/i");
    tree->Branch("tau_isolationNeutralHadrCands_size", tau_isolationNeutralHadrCands_size, "tau_isolationNeutralHadrCands_size[tau_count]/i");
    tree->Branch("tau_isolationGammaCands_size", tau_isolationGammaCands_size, "tau_isolationGammaCands_size[tau_count]/i");

    tree->Branch("tau_genDecayMode_name", tau_genDecayMode_name, "tau_genDecayMode_name[tau_count]/C");
    tree->Branch("tau_genDecayMode", tau_genDecayMode, "tau_genDecayMode[tau_count]/I");
    tree->Branch("tau_decayMode_name", tau_decayMode_name, "tau_decayMode_name[tau_count]/C");
    tree->Branch("tau_decayMode", tau_decayMode, "tau_decayMode[tau_count]/I");

    tree->Branch("tau_constituents_count", tau_constituents_count, "tau_constituents_count[tau_count]/i");
    tree->Branch("tau_constituents_px", tau_constituents_px, "tau_constituents_px[tau_count][50]/F");
    tree->Branch("tau_constituents_py", tau_constituents_py, "tau_constituents_py[tau_count][50]/F");
    tree->Branch("tau_constituents_pz", tau_constituents_pz, "tau_constituents_pz[tau_count][50]/F");
    tree->Branch("tau_constituents_e", tau_constituents_e, "tau_constituents_e[tau_count][50]/F");
    tree->Branch("tau_constituents_mass", tau_constituents_mass, "tau_constituents_mass[tau_count][50]/F");
    tree->Branch("tau_constituents_charge", tau_constituents_charge, "tau_constituents_charge[tau_count][50]/I");
    tree->Branch("tau_constituents_vx", tau_constituents_vx, "tau_constituents_vx[tau_count][50]/F");
    tree->Branch("tau_constituents_vy", tau_constituents_vy, "tau_constituents_vy[tau_count][50]/F");
    tree->Branch("tau_constituents_vz", tau_constituents_vz, "tau_constituents_vz[tau_count][50]/F");
    tree->Branch("tau_constituents_pdgId", tau_constituents_pdgId, "tau_constituents_pdgId[tau_count][50]/I");
    tree->Branch("tau_constituents_lostPixelHits", tau_constituents_lostPixelHits, "tau_constituents_lostPixelHits[tau_count][50]/I");
  }

  if (crectrack) {
    tree->Branch("track_count", &track_count, "track_count/i");
    tree->Branch("track_px", track_px, "track_px[track_count]/F");
    tree->Branch("track_py", track_py, "track_py[track_count]/F");
    tree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
    tree->Branch("track_pt", track_pt, "track_pt[track_count]/F");
    tree->Branch("track_eta", track_eta, "track_eta[track_count]/F");
    tree->Branch("track_phi", track_phi, "track_phi[track_count]/F");
    tree->Branch("track_charge", track_charge, "track_charge[track_count]/F");
    tree->Branch("track_mass", track_mass, "track_mass[track_count]/F");
    tree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
    tree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
    tree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
    tree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
    tree->Branch("track_vx", track_vx, "track_vx[track_count]/F");
    tree->Branch("track_vy", track_vy, "track_vy[track_count]/F");
    tree->Branch("track_vz", track_vz, "track_vz[track_count]/F");
    tree->Branch("track_ID", track_ID, "track_ID[track_count]/I");
    tree->Branch("track_highPurity", track_highPurity, "track_highPurity[track_count]/O");
  }

  // Met
  if (crecpfmet) {
    tree->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
    tree->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
    tree->Branch("pfmet_ez", &pfmet_ey, "pfmet_ez/F");
    tree->Branch("pfmet_pt", &pfmet_pt, "pfmet_pt/F");
    tree->Branch("pfmet_phi", &pfmet_phi, "pfmet_phi/F");
    tree->Branch("pfmet_sigxx", &pfmet_sigxx, "pfmet_sigxx/F");
    tree->Branch("pfmet_sigxy", &pfmet_sigxy, "pfmet_sigxy/F");
    tree->Branch("pfmet_sigyx", &pfmet_sigyx, "pfmet_sigyx/F");
    tree->Branch("pfmet_sigyy", &pfmet_sigyy, "pfmet_sigyy/F");
    tree->Branch("pfmet_sig", &pfmet_sig, "pfmet_sig/F");

    tree->Branch("genmet_ex", &genmet_ex, "genmet_ex/F");
    tree->Branch("genmet_ey", &genmet_ey, "genmet_ey/F");

    tree->Branch("pfmet_ex_JetEnUp", &pfmet_ex_JetEnUp, "pfmet_ex_JetEnUp/F");
    tree->Branch("pfmet_ey_JetEnUp", &pfmet_ey_JetEnUp, "pfmet_ey_JetEnUp/F");

    tree->Branch("pfmet_ex_JetEnDown", &pfmet_ex_JetEnDown, "pfmet_ex_JetEnDown/F");
    tree->Branch("pfmet_ey_JetEnDown", &pfmet_ey_JetEnDown, "pfmet_ey_JetEnDown/F");

    tree->Branch("pfmet_ex_UnclusteredEnUp", &pfmet_ex_UnclusteredEnUp, "pfmet_ex_UnclusteredEnUp/F");
    tree->Branch("pfmet_ey_UnclusteredEnUp", &pfmet_ey_UnclusteredEnUp, "pfmet_ey_UnclusteredEnUp/F");

    tree->Branch("pfmet_ex_UnclusteredEnDown", &pfmet_ex_UnclusteredEnDown, "pfmet_ex_UnclusteredEnDown/F");
    tree->Branch("pfmet_ey_UnclusteredEnDown", &pfmet_ey_UnclusteredEnDown, "pfmet_ey_UnclusteredEnDown/F");

  }

  if (crecpfmetcorr) {
    tree->Branch("pfmetcorr_ex", &pfmetcorr_ex, "pfmetcorr_ex/F");
    tree->Branch("pfmetcorr_ey", &pfmetcorr_ey, "pfmetcorr_ey/F");
    tree->Branch("pfmetcorr_ez", &pfmetcorr_ey, "pfmetcorr_ez/F");
    tree->Branch("pfmetcorr_pt", &pfmetcorr_pt, "pfmetcorr_pt/F");
    tree->Branch("pfmetcorr_phi", &pfmetcorr_phi, "pfmetcorr_phi/F");
    tree->Branch("pfmetcorr_sigxx", &pfmetcorr_sigxx, "pfmetcorr_sigxx/F");
    tree->Branch("pfmetcorr_sigxy", &pfmetcorr_sigxy, "pfmetcorr_sigxy/F");
    tree->Branch("pfmetcorr_sigyx", &pfmetcorr_sigyx, "pfmetcorr_sigyx/F");
    tree->Branch("pfmetcorr_sigyy", &pfmetcorr_sigyy, "pfmetcorr_sigyy/F");
    tree->Branch("pfmetcorr_sig", &pfmetcorr_sig, "pfmetcorr_sig/F");

    tree->Branch("pfmetcorr_ex_JetEnUp", &pfmetcorr_ex_JetEnUp, "pfmetcorr_ex_JetEnUp/F");
    tree->Branch("pfmetcorr_ey_JetEnUp", &pfmetcorr_ey_JetEnUp, "pfmetcorr_ey_JetEnUp/F");

    tree->Branch("pfmetcorr_ex_JetEnDown", &pfmetcorr_ex_JetEnDown, "pfmetcorr_ex_JetEnDown/F");
    tree->Branch("pfmetcorr_ey_JetEnDown", &pfmetcorr_ey_JetEnDown, "pfmetcorr_ey_JetEnDown/F");

    tree->Branch("pfmetcorr_ex_UnclusteredEnUp", &pfmetcorr_ex_UnclusteredEnUp, "pfmetcorr_ex_UnclusteredEnUp/F");
    tree->Branch("pfmetcorr_ey_UnclusteredEnUp", &pfmetcorr_ey_UnclusteredEnUp, "pfmetcorr_ey_UnclusteredEnUp/F");

    tree->Branch("pfmetcorr_ex_UnclusteredEnDown", &pfmetcorr_ex_UnclusteredEnDown, "pfmetcorr_ex_UnclusteredEnDown/F");
    tree->Branch("pfmetcorr_ey_UnclusteredEnDown", &pfmetcorr_ey_UnclusteredEnDown, "pfmetcorr_ey_UnclusteredEnDown/F");

    tree->Branch("pfmetcorr_ex_JetResUp", &pfmetcorr_ex_JetResUp, "pfmetcorr_ex_JetResUp/F");
    tree->Branch("pfmetcorr_ey_JetResUp", &pfmetcorr_ey_JetResUp, "pfmetcorr_ey_JetResUp/F");

    tree->Branch("pfmetcorr_ex_JetResDown", &pfmetcorr_ex_JetResDown, "pfmetcorr_ex_JetResDown/F");
    tree->Branch("pfmetcorr_ey_JetResDown", &pfmetcorr_ey_JetResDown, "pfmetcorr_ey_JetResDown/F");

/*    tree->Branch("pfmetcorr_ex_smeared", &pfmetcorr_ex_smeared, "pfmetcorr_ex_smeared/F");
    tree->Branch("pfmetcorr_ey_smeared", &pfmetcorr_ey_smeared, "pfmetcorr_ey_smeared/F");
    tree->Branch("pfmetcorr_ez_smeared", &pfmetcorr_ey_smeared, "pfmetcorr_ez_smeared/F");
    tree->Branch("pfmetcorr_pt_smeared", &pfmetcorr_pt_smeared, "pfmetcorr_pt_smeared/F");
    tree->Branch("pfmetcorr_phi_smeared", &pfmetcorr_phi_smeared, "pfmetcorr_phi_smeared/F");
    tree->Branch("pfmetcorr_ex_JetEnUp_smeared", &pfmetcorr_ex_JetEnUp_smeared, "pfmetcorr_ex_JetEnUp_smeared/F");
    tree->Branch("pfmetcorr_ey_JetEnUp_smeared", &pfmetcorr_ey_JetEnUp_smeared, "pfmetcorr_ey_JetEnUp_smeared/F");

    tree->Branch("pfmetcorr_ex_JetEnDown_smeared", &pfmetcorr_ex_JetEnDown_smeared, "pfmetcorr_ex_JetEnDown_smeared/F");
    tree->Branch("pfmetcorr_ey_JetEnDown_smeared", &pfmetcorr_ey_JetEnDown_smeared, "pfmetcorr_ey_JetEnDown_smeared/F");

    tree->Branch("pfmetcorr_ex_UnclusteredEnUp_smeared", &pfmetcorr_ex_UnclusteredEnUp_smeared, "pfmetcorr_ex_UnclusteredEnUp_smeared/F");
    tree->Branch("pfmetcorr_ey_UnclusteredEnUp_smeared", &pfmetcorr_ey_UnclusteredEnUp_smeared, "pfmetcorr_ey_UnclusteredEnUp_smeared/F");

    tree->Branch("pfmetcorr_ex_UnclusteredEnDown_smeared", &pfmetcorr_ex_UnclusteredEnDown_smeared, "pfmetcorr_ex_UnclusteredEnDown_smeared/F");
    tree->Branch("pfmetcorr_ey_UnclusteredEnDown_smeared", &pfmetcorr_ey_UnclusteredEnDown_smeared, "pfmetcorr_ey_UnclusteredEnDown_smeared/F");

    tree->Branch("pfmetcorr_ex_JetResUp_smeared", &pfmetcorr_ex_JetResUp_smeared, "pfmetcorr_ex_JetResUp_smeared/F");
    tree->Branch("pfmetcorr_ey_JetResUp_smeared", &pfmetcorr_ey_JetResUp_smeared, "pfmetcorr_ey_JetResUp_smeared/F");

    tree->Branch("pfmetcorr_ex_JetResDown_smeared", &pfmetcorr_ex_JetResDown_smeared, "pfmetcorr_ex_JetResDown_smeared/F");
    tree->Branch("pfmetcorr_ey_JetResDown_smeared", &pfmetcorr_ey_JetResDown_smeared, "pfmetcorr_ey_JetResDown_smeared/F");
*/
  }


  if (crecpuppimet) {
    tree->Branch("puppimet_ex", &puppimet_ex, "puppimet_ex/F");
    tree->Branch("puppimet_ey", &puppimet_ey, "puppimet_ey/F");
    tree->Branch("puppimet_ez", &puppimet_ey, "puppimet_ez/F");
    tree->Branch("puppimet_pt", &puppimet_pt, "puppimet_pt/F");
    tree->Branch("puppimet_phi", &puppimet_phi, "puppimet_phi/F");
    tree->Branch("puppimet_sigxx", &puppimet_sigxx, "puppimet_sigxx/F");
    tree->Branch("puppimet_sigxy", &puppimet_sigxy, "puppimet_sigxy/F");
    tree->Branch("puppimet_sigyx", &puppimet_sigyx, "puppimet_sigyx/F");
    tree->Branch("puppimet_sigyy", &puppimet_sigyy, "puppimet_sigyy/F");

    tree->Branch("puppimet_ex_JetEnUp", &puppimet_ex_JetEnUp, "puppimet_ex_JetEnUp/F");
    tree->Branch("puppimet_ey_JetEnUp", &puppimet_ey_JetEnUp, "puppimet_ey_JetEnUp/F");

    tree->Branch("puppimet_ex_JetEnDown", &puppimet_ex_JetEnDown, "puppimet_ex_JetEnDown/F");
    tree->Branch("puppimet_ey_JetEnDown", &puppimet_ey_JetEnDown, "puppimet_ey_JetEnDown/F");

    tree->Branch("puppimet_ex_UnclusteredEnUp", &puppimet_ex_UnclusteredEnUp, "puppimet_ex_UnclusteredEnUp/F");
    tree->Branch("puppimet_ey_UnclusteredEnUp", &puppimet_ey_UnclusteredEnUp, "puppimet_ey_UnclusteredEnUp/F");

    tree->Branch("puppimet_ex_UnclusteredEnDown", &puppimet_ex_UnclusteredEnDown, "puppimet_ex_UnclusteredEnDown/F");
    tree->Branch("puppimet_ey_UnclusteredEnDown", &puppimet_ey_UnclusteredEnDown, "puppimet_ey_UnclusteredEnDown/F");


    tree->Branch("puppimet_ex_JetResUp", &puppimet_ex_JetResUp, "puppimet_ex_JetResUp/F");
    tree->Branch("puppimet_ey_JetResUp", &puppimet_ey_JetResUp, "puppimet_ey_JetResUp/F");

    tree->Branch("puppimet_ex_JetResDown", &puppimet_ex_JetResDown, "puppimet_ex_JetResDown/F");
    tree->Branch("puppimet_ey_JetResDown", &puppimet_ey_JetResDown, "puppimet_ey_JetResDown/F");

  }



  if (crecmvamet) {
    tree->Branch("mvamet_count", &mvamet_count, "mvamet_count/i");
    tree->Branch("mvamet_ex", mvamet_ex, "mvamet_ex[mvamet_count]/F");
    tree->Branch("mvamet_ey", mvamet_ey, "mvamet_ey[mvamet_count]/F");
    tree->Branch("mvamet_sigxx", mvamet_sigxx, "mvamet_sigxx[mvamet_count]/F");
    tree->Branch("mvamet_sigxy", mvamet_sigxy, "mvamet_sigxy[mvamet_count]/F");
    tree->Branch("mvamet_sigyx", mvamet_sigyx, "mvamet_sigyx[mvamet_count]/F");
    tree->Branch("mvamet_sigyy", mvamet_sigyy, "mvamet_sigyy[mvamet_count]/F");
    tree->Branch("mvamet_channel", mvamet_channel, "mvamet_channel[mvamet_count]/b");
    tree->Branch("mvamet_lep1", mvamet_lep1, "mvamet_lep1[mvamet_count]/i");
    tree->Branch("mvamet_lep2", mvamet_lep2, "mvamet_lep2[mvamet_count]/i");
    tree->Branch("mvamet_lep1_pt", mvamet_lep1_pt, "mvamet_lep1_pt[mvamet_count]/F");
    tree->Branch("mvamet_lep2_pt", mvamet_lep2_pt, "mvamet_lep2_pt[mvamet_count]/F");
  }

  // generator info

  if(crecstxs){
    tree->Branch("htxs_stage0cat",&htxs_stage0cat,"htxs_stage0cat/I");
    tree->Branch("htxs_stage1p1cat",&htxs_stage1p1cat,"htxs_stage1p1cat/I");
    tree->Branch("htxs_higgsPt",&htxs_higgsPt,"htxs_higgsPt/F");
    tree->Branch("htxs_njets30",&htxs_njets30,"htxs_njets30/I");
  }

  if (cgen) {

    if (!cdata) {
      tree->Branch("genid1", &genid1, "genid1/F");
      tree->Branch("genx1", &genx1, "genx1/F");
      tree->Branch("genid2", &genid2, "genid2/F");
      tree->Branch("genx2", &genx2, "genx2/F");
      tree->Branch("genScale", &genScale, "genScale/F");

      for (int iScale=0; iScale<9; ++iScale) {
	char number[4];
	sprintf(number,"%1i",iScale);
	TString Number(number);
	tree->Branch("weightScale"+Number,&weightScale[iScale],"weightScale"+Number+"/F");
      }

      tree->Branch("weightPDFmax",&weightPDFmax,"weightPDFmax/F");
      tree->Branch("weightPDFmin",&weightPDFmin,"weightPDFmin/F");
      tree->Branch("weightPDFmean",&weightPDFmean,"weightPDFmean/F");
      tree->Branch("weightPDFup",&weightPDFup,"weightPDFup/F");
      tree->Branch("weightPDFdown",&weightPDFdown,"weightPDFdown/F");
      tree->Branch("weightPDFvar",&weightPDFvar,"weightPDFvar/F");

      tree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
      tree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
      tree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
      tree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");

    }

    if (!cdata || cembedded) {

      tree->Branch("genweight", &genweight, "genweight/F");

      // generated taus
      tree->Branch("gentau_count", &gentau_count, "gentau_count/i");
      tree->Branch("gentau_e",  gentau_e,  "tau_e[gentau_count]/F");
      tree->Branch("gentau_charge",  gentau_charge,  "tau_charge[gentau_count]/F");
      tree->Branch("gentau_px", gentau_px, "tau_px[gentau_count]/F");
      tree->Branch("gentau_py", gentau_py, "tau_py[gentau_count]/F");
      tree->Branch("gentau_pz", gentau_pz, "tau_pz[gentau_count]/F");

      tree->Branch("gentau_visible_e",  gentau_visible_e,  "tau_visible_e[gentau_count]/F");
      tree->Branch("gentau_visible_px", gentau_visible_px, "tau_visible_px[gentau_count]/F");
      tree->Branch("gentau_visible_py", gentau_visible_py, "tau_visible_py[gentau_count]/F");
      tree->Branch("gentau_visible_pz", gentau_visible_pz, "tau_visible_pz[gentau_count]/F");

      tree->Branch("gentau_visible_pt",  gentau_visible_pt,  "tau_visible_pt[gentau_count]/F");
      tree->Branch("gentau_visible_eta", gentau_visible_eta, "tau_visible_eta[gentau_count]/F");
      tree->Branch("gentau_visible_phi", gentau_visible_phi, "tau_visible_phi[gentau_count]/F");
      tree->Branch("gentau_visible_mass", gentau_visible_mass, "tau_visible_mass[gentau_count]/F");

      tree->Branch("gentau_visibleNoLep_e",  gentau_visibleNoLep_e,  "tau_visibleNoLep_e[gentau_count]/F");
      tree->Branch("gentau_visibleNoLep_pt",  gentau_visibleNoLep_pt,  "tau_visibleNoLep_pt[gentau_count]/F");
      tree->Branch("gentau_visibleNoLep_eta", gentau_visibleNoLep_eta, "tau_visibleNoLep_eta[gentau_count]/F");
      tree->Branch("gentau_visibleNoLep_phi", gentau_visibleNoLep_phi, "tau_visibleNoLep_phi[gentau_count]/F");
      tree->Branch("gentau_visibleNoLep_mass", gentau_visibleNoLep_mass, "tau_visibleNoLep_mass[gentau_count]/F");

      tree->Branch("gentau_status", gentau_status, "gentau_status[gentau_count]/I");
      tree->Branch("gentau_fromHardProcess", gentau_fromHardProcess, "gentau_fromHardProcess[gentau_count]/I");
      tree->Branch("gentau_fromHardProcessBeforeFSR", gentau_fromHardProcessBeforeFSR, "gentau_fromHardProcessBeforeFSR[gentau_count]/I");
      tree->Branch("gentau_isDecayedLeptonHadron", gentau_isDecayedLeptonHadron, "gentau_isDecayedLeptonHadron[gentau_count]/I");
      tree->Branch("gentau_isDirectHadronDecayProduct", gentau_isDirectHadronDecayProduct, "gentau_isDirectHadronDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isDirectHardProcessTauDecayProduct", gentau_isDirectHardProcessTauDecayProduct, "gentau_isDirectHardProcessTauDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isDirectPromptTauDecayProduct", gentau_isDirectPromptTauDecayProduct, "gentau_isDirectPromptTauDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isDirectTauDecayProduct", gentau_isDirectTauDecayProduct, "gentau_isDirectTauDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isFirstCopy", gentau_isFirstCopy, "gentau_isFirstCopy[gentau_count]/I");
      tree->Branch("gentau_isHardProcess", gentau_isHardProcess, "gentau_isHardProcess[gentau_count]/I");
      tree->Branch("gentau_isHardProcessTauDecayProduct", gentau_isHardProcessTauDecayProduct, "gentau_isHardProcessTauDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isLastCopy", gentau_isLastCopy, "gentau_isLastCopy[gentau_count]/I");
      tree->Branch("gentau_isLastCopyBeforeFSR", gentau_isLastCopyBeforeFSR, "gentau_isLastCopyBeforeFSR[gentau_count]/I");
      tree->Branch("gentau_isPrompt", gentau_isPrompt, "gentau_isPrompt[gentau_count]/I");
      tree->Branch("gentau_isPromptTauDecayProduct", gentau_isPromptTauDecayProduct, "gentau_isPromptTauDecayProduct[gentau_count]/I");
      tree->Branch("gentau_isTauDecayProduct", gentau_isTauDecayProduct, "gentau_isTauDecayProduct[gentau_count]/I");

      tree->Branch("gentau_decayMode",  gentau_decayMode,  "tau_decayMode[gentau_count]/I");
      tree->Branch("gentau_decayMode_name",  gentau_decayMode_name,  "tau_decayMode_name[gentau_count]/C");
      tree->Branch("gentau_mother",gentau_mother,"gentau_mother[gentau_count]/b");

      // generated particles
      tree->Branch("genparticles_lheHt", &genparticles_lheHt, "genparticles_lheHt/F");
      tree->Branch("genparticles_lheWPt", &genparticles_lheWPt, "genparticles_lheWPt/F");
      tree->Branch("genparticles_noutgoing", &genparticles_noutgoing, "genparticles_noutgoing/i");
      tree->Branch("genparticles_noutgoing_NLO", &genparticles_noutgoing_NLO, "genparticles_noutgoing_NLO/I");
      tree->Branch("genparticles_count", &genparticles_count, "genparticles_count/i");
      tree->Branch("genparticles_e", genparticles_e, "genparticles_e[genparticles_count]/F");
      tree->Branch("genparticles_px", genparticles_px, "genparticles_px[genparticles_count]/F");
      tree->Branch("genparticles_py", genparticles_py, "genparticles_py[genparticles_count]/F");
      tree->Branch("genparticles_pz", genparticles_pz, "genparticles_pz[genparticles_count]/F");
      tree->Branch("genparticles_vx", genparticles_vx, "genparticles_vx[genparticles_count]/F");
      tree->Branch("genparticles_vy", genparticles_vy, "genparticles_vy[genparticles_count]/F");
      tree->Branch("genparticles_vz", genparticles_vz, "genparticles_vz[genparticles_count]/F");
      tree->Branch("genparticles_pdgid", genparticles_pdgid, "genparticles_pdgid[genparticles_count]/I");
      tree->Branch("genparticles_status", genparticles_status, "genparticles_status[genparticles_count]/I");
      tree->Branch("genparticles_info", genparticles_info, "genparticles_info[genparticles_count]/i");

      tree->Branch("genparticles_fromHardProcess", genparticles_fromHardProcess, "genparticles_fromHardProcess[genparticles_count]/I");
      tree->Branch("genparticles_fromHardProcessBeforeFSR", genparticles_fromHardProcessBeforeFSR, "genparticles_fromHardProcessBeforeFSR[genparticles_count]/I");
      tree->Branch("genparticles_isDecayedLeptonHadron", genparticles_isDecayedLeptonHadron, "genparticles_isDecayedLeptonHadron[genparticles_count]/I");
      tree->Branch("genparticles_isDirectHadronDecayProduct", genparticles_isDirectHadronDecayProduct, "genparticles_isDirectHadronDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isDirectHardProcessTauDecayProduct", genparticles_isDirectHardProcessTauDecayProduct, "genparticles_isDirectHardProcessTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isDirectPromptTauDecayProduct", genparticles_isDirectPromptTauDecayProduct, "genparticles_isDirectPromptTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isDirectTauDecayProduct", genparticles_isDirectTauDecayProduct, "genparticles_isDirectTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isFirstCopy", genparticles_isFirstCopy, "genparticles_isFirstCopy[genparticles_count]/I");
      tree->Branch("genparticles_isHardProcess", genparticles_isHardProcess, "genparticles_isHardProcess[genparticles_count]/I");
      tree->Branch("genparticles_isHardProcessTauDecayProduct", genparticles_isHardProcessTauDecayProduct, "genparticles_isHardProcessTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isLastCopy", genparticles_isLastCopy, "genparticles_isLastCopy[genparticles_count]/I");
      tree->Branch("genparticles_isLastCopyBeforeFSR", genparticles_isLastCopyBeforeFSR, "genparticles_isLastCopyBeforeFSR[genparticles_count]/I");
      tree->Branch("genparticles_isPrompt", genparticles_isPrompt, "genparticles_isPrompt[genparticles_count]/I");
      tree->Branch("genparticles_isPromptTauDecayProduct", genparticles_isPromptTauDecayProduct, "genparticles_isPromptTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_isTauDecayProduct", genparticles_isTauDecayProduct, "genparticles_isTauDecayProduct[genparticles_count]/I");
      tree->Branch("genparticles_mother", genparticles_mother, "genparticles_mother[genparticles_count]/b");

      tree->Branch("genjets_count", &genjets_count, "genjets_count/i");
      tree->Branch("genjets_e", genjets_e, "genjets_e[genjets_count]/F");
      tree->Branch("genjets_px", genjets_px, "genjets_px[genjets_count]/F");
      tree->Branch("genjets_py", genjets_py, "genjets_py[genjets_count]/F");
      tree->Branch("genjets_pz", genjets_pz, "genjets_pz[genjets_count]/F");
      tree->Branch("genjets_pt", genjets_pt, "genjets_pt[genjets_count]/F");
      tree->Branch("genjets_eta", genjets_eta, "genjets_eta[genjets_count]/F");
      tree->Branch("genjets_phi", genjets_phi, "genjets_phi[genjets_count]/F");
      tree->Branch("genjets_pdgid", genjets_pdgid, "genjets_pdgid[genjets_count]/I");
      tree->Branch("genjets_status", genjets_status, "genjets_status[genjets_count]/I");
      tree->Branch("genjets_em_energy", genjets_em_energy, "genjets_em_energy[genjets_count]/F");
      tree->Branch("genjets_had_energy", genjets_had_energy, "genjets_had_energy[genjets_count]/F");
      tree->Branch("genjets_invisible_energy", genjets_invisible_energy, "genjets_invisible_energy[genjets_count]/F");
      tree->Branch("genjets_auxiliary_energy", genjets_auxiliary_energy, "genjets_auxiliary_energy[genjets_count]/F");
    }
  }
  // TauSpinner branches and definitions
  if (crecTauSpinner) {
    std::string TauSpinnerSettingsPDF="NNPDF30_nlo_as_0118";
    bool Ipp=true;
    int Ipol=0;
    int nonSM2=0;
    int nonSMN=0;
    double CMSENE=13000.0;
    
    Tauolapp::Tauola::initialize();
    LHAPDF::initPDFSetByName(TauSpinnerSettingsPDF);
    TauSpinner::initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);
    
    tree->Branch("TauSpinAngles_count", &TauSpinAngles_count, "TauSpinAngles_count/i");
    tree->Branch("TauSpinnerWeight", TauSpinnerWeight, "TauSpinnerWeight[TauSpinAngles_count]/D");

  }

  // SUSY Info
  if (csusyinfo) {
    tree->Branch("SusyMotherMass",&SusyMotherMass,"SusyMotherMass/F");
    tree->Branch("SusyLSPMass",&SusyLSPMass,"SusyLSPMass/F");
  }

  // L1 objects
  tree->Branch("l1muon_count",       &l1muon_count,       "l1muon_count/i");
  tree->Branch("l1muon_px",           l1muon_px,          "l1muon_px[l1muon_count]/F");
  tree->Branch("l1muon_py",           l1muon_py,          "l1muon_py[l1muon_count]F");
  tree->Branch("l1muon_pz",           l1muon_pz,          "l1muon_pz[l1muon_count]/F");
  tree->Branch("l1muon_pt",           l1muon_pt,          "l1muon_pt[l1muon_count]/F");
//  tree->Branch("l1muon_ipt",          l1muon_ipt,         "l1muon_ipt[l1muon_count]/I");
  tree->Branch("l1muon_eta",          l1muon_eta,         "l1muon_eta[l1muon_count]/I");
  tree->Branch("l1muon_phi",          l1muon_phi,         "l1muon_phi[l1muon_count]/I");
//  tree->Branch("l1muon_qual",         l1muon_qual,        "l1muon_qual[l1muon_count]/I");
  tree->Branch("l1muon_iso",          l1muon_iso,         "l1muon_iso[l1muon_count]/I");
  tree->Branch("l1muon_charge",       l1muon_charge,      "l1muon_charge[l1muon_count]/I");
//  tree->Branch("l1muon_chargeValid",  l1muon_chargeValid, "l1muon_chargeValid[l1muon_count]/I");
  tree->Branch("l1muon_muonIndex",    l1muon_muonIndex,   "l1muon_muonIndex[l1muon_count]/I");
//  tree->Branch("l1muon_tag",          l1muon_tag,         "l1muon_tag[l1muon_count]/I");
//  tree->Branch("l1muon_isoSum",       l1muon_isoSum,      "l1muon_isoSum[l1muon_count]/I");
//  tree->Branch("l1muon_dPhiExtra",    l1muon_dPhiExtra,   "l1muon_dPhiExtra[l1muon_count]/I");
//  tree->Branch("l1muon_dEtaExtra",    l1muon_dEtaExtra,   "l1muon_dEtaExtra[l1muon_count]/I");
//  tree->Branch("l1muon_rank",         l1muon_rank,        "l1muon_rank[l1muon_count]/I");
//  tree->Branch("l1muon_bx",           l1muon_bx,          "l1muon_bx[l1muon_count]/I");

  tree->Branch("l1egamma_count",       &l1egamma_count,       "l1egamma_count/i");
  tree->Branch("l1egamma_px",           l1egamma_px,          "l1egamma_px[l1egamma_count]/F");
  tree->Branch("l1egamma_py",           l1egamma_py,          "l1egamma_py[l1egamma_count]F");
  tree->Branch("l1egamma_pz",           l1egamma_pz,          "l1egamma_pz[l1egamma_count]/F");
  tree->Branch("l1egamma_pt",           l1egamma_pt,          "l1egamma_pt[l1egamma_count]/F");
//  tree->Branch("l1egamma_ipt",          l1egamma_ipt,         "l1egamma_ipt[l1egamma_count]/I");
  tree->Branch("l1egamma_eta",          l1egamma_eta,         "l1egamma_eta[l1egamma_count]/I");
  tree->Branch("l1egamma_phi",          l1egamma_phi,         "l1egamma_phi[l1egamma_count]/I");
//  tree->Branch("l1egamma_qual",         l1egamma_qual,        "l1egamma_qual[l1egamma_count]/I");
  tree->Branch("l1egamma_iso",          l1egamma_iso,         "l1egamma_iso[l1egamma_count]/I");
//  tree->Branch("l1egamma_towerIEta",    l1egamma_towerIEta,   "l1egamma_towerIEta[l1egamma_count]/I");
//  tree->Branch("l1egamma_towerIPhi",    l1egamma_towerIPhi,   "l1egamma_towerIPhi[l1egamma_count]/I");
//  tree->Branch("l1egamma_rawEt",        l1egamma_rawEt,       "l1egamma_rawEt[l1egamma_count]/I");
//  tree->Branch("l1egamma_isoEt",        l1egamma_isoEt,       "l1egamma_isoEt[l1egamma_count]/I");
//  tree->Branch("l1egamma_footprintEt",  l1egamma_footprintEt, "l1egamma_footprintEt[l1egamma_count]/I");
//  tree->Branch("l1egamma_nTT",          l1egamma_nTT,         "l1egamma_nTT[l1egamma_count]/I");
//  tree->Branch("l1egamma_shape",        l1egamma_shape,       "l1egamma_shape[l1egamma_count]/I");
//  tree->Branch("l1egamma_bx",           l1egamma_bx,          "l1egamma_bx[l1egamma_count]/I");

  tree->Branch("l1tau_count",       &l1tau_count,       "l1tau_count/i");
  tree->Branch("l1tau_px",           l1tau_px,          "l1tau_px[l1tau_count]/F");
  tree->Branch("l1tau_py",           l1tau_py,          "l1tau_py[l1tau_count]F");
  tree->Branch("l1tau_pz",           l1tau_pz,          "l1tau_pz[l1tau_count]/F");
  tree->Branch("l1tau_pt",           l1tau_pt,          "l1tau_pt[l1tau_count]/F");
  tree->Branch("l1tau_ipt",          l1tau_ipt,         "l1tau_ipt[l1tau_count]/I");
  tree->Branch("l1tau_eta",          l1tau_eta,         "l1tau_eta[l1tau_count]/I");
  tree->Branch("l1tau_phi",          l1tau_phi,         "l1tau_phi[l1tau_count]/I");
  tree->Branch("l1tau_qual",         l1tau_qual,        "l1tau_qual[l1tau_count]/I");
  tree->Branch("l1tau_iso",          l1tau_iso,         "l1tau_iso[l1tau_count]/I");
  tree->Branch("l1tau_towerIEta",    l1tau_towerIEta,   "l1tau_towerIEta[l1tau_count]/I");
  tree->Branch("l1tau_towerIPhi",    l1tau_towerIPhi,   "l1tau_towerIPhi[l1tau_count]/I");
  tree->Branch("l1tau_rawEt",        l1tau_rawEt,       "l1tau_rawEt[l1tau_count]/I");
  tree->Branch("l1tau_isoEt",        l1tau_isoEt,       "l1tau_isoEt[l1tau_count]/I");
  tree->Branch("l1tau_nTT",          l1tau_nTT,         "l1tau_nTT[l1tau_count]/I");
  tree->Branch("l1tau_hasEM",        l1tau_hasEM,       "l1tau_hasEM[l1tau_count]/I");
  tree->Branch("l1tau_isMerged",     l1tau_isMerged,    "l1tau_isMerged[l1tau_count]/I");
  tree->Branch("l1tau_bx",           l1tau_bx,          "l1tau_bx[l1tau_count]/I");

  tree->Branch("l1isotau_count", &l1isotau_count, "l1isotau_count/i");
  tree->Branch("l1isotau_e", l1isotau_e, "l1isotau_e[l1isotau_count]/F");
  tree->Branch("l1isotau_px", l1isotau_px, "l1isotau_px[l1isotau_count]/F");
  tree->Branch("l1isotau_py", l1isotau_py, "l1isotau_py[l1isotau_count]/F");
  tree->Branch("l1isotau_pz", l1isotau_pz, "l1isotau_pz[l1isotau_count]/F");
  tree->Branch("l1isotau_mass", l1isotau_mass, "l1isotau_mass[l1isotau_count]/F");
  tree->Branch("l1isotau_eta", l1isotau_eta, "l1isotau_eta[l1isotau_count]/F");
  tree->Branch("l1isotau_phi", l1isotau_phi, "l1isotau_phi[l1isotau_count]/F");
  tree->Branch("l1isotau_pt", l1isotau_pt, "l1isotau_pt[l1isotau_count]/F");
  tree->Branch("l1isotau_charge", l1isotau_charge, "l1isotau_charge[l1isotau_count]/F");
  tree->Branch("l1isotau_iso",    l1isotau_iso,    "l1isotau_iso[l1isotau_count]/I");

  // trigger objects
  if (ctrigger) {
    tree->Branch("trigobject_count",&trigobject_count,"trigobject_count/i");
    tree->Branch("trigobject_px",trigobject_px,"trigobject_px[trigobject_count]/F");
    tree->Branch("trigobject_py",trigobject_py,"trigobject_py[trigobject_count]/F");
    tree->Branch("trigobject_pz",trigobject_pz,"trigobject_pz[trigobject_count]/F");
    tree->Branch("trigobject_pt",trigobject_pt,"trigobject_pt[trigobject_count]/F");
    tree->Branch("trigobject_eta",trigobject_eta,"trigobject_eta[trigobject_count]/F");
    tree->Branch("trigobject_phi",trigobject_phi,"trigobject_phi[trigobject_count]/F");
    tree->Branch("trigobject_filters",trigobject_filters,TString("trigobject_filters[trigobject_count][")+(Long64_t)M_hltfiltersmax+"]/O");
    tree->Branch("trigobject_isMuon",trigobject_isMuon,"trigobject_isMuon[trigobject_count]/O");
    tree->Branch("trigobject_isElectron",trigobject_isElectron,"trigobject_isElectron[trigobject_count]/O");
    tree->Branch("trigobject_isTau",trigobject_isTau,"trigobject_isTau[trigobject_count]/O");
    tree->Branch("trigobject_isJet",trigobject_isJet,"trigobject_isJet[trigobject_count]/O");
    tree->Branch("trigobject_isMET",trigobject_isMET,"trigobject_isMET[trigobject_count]/O");
  }

  //add these branches to main tree as well
  tree->Branch("run_hltnames", "std::vector<std::string>", &run_hltnames);
  tree->Branch("run_hltfilters", "std::vector<std::string>",&run_hltfilters);
  tree->Branch("run_hltmufilters", "std::vector<std::string>", &run_hltmufilters);
  tree->Branch("run_hltelectronfilters", "std::vector<std::string>", &run_hltelectronfilters);
  tree->Branch("run_hlttaufilters", "std::vector<std::string>", &run_hlttaufilters);
  tree->Branch("run_hltphotonfilters", "std::vector<std::string>", &run_hltphotonfilters);
  tree->Branch("run_hltjetfilters", "std::vector<std::string>", &run_hltjetfilters);
  tree->Branch("run_floattaudiscriminators", "std::vector<std::string>", &run_floattaudiscriminators);
  tree->Branch("run_binarytaudiscriminators", "std::vector<std::string>", &run_binarytaudiscriminators);
  tree->Branch("run_btagdiscriminators", "std::vector<std::string>", &run_btagdiscriminators);
  hltriggerresults_ = new std::map<std::string, int>();
  hltriggerprescales_ = new std::map<std::string, int>();
  tree->Branch("hltriggerresults", "std::map<std::string, int>", &hltriggerresults_);
  tree->Branch("hltriggerprescales", "std::map<std::string, int>", &hltriggerprescales_);
  tree->Branch("hltriggerresultsV", "std::vector<std::string>", &hltriggerresultsV_);
  // tree->Branch("embeddingWeight", &embeddingWeight_, "embeddingWeight/F");

  // add flags
  flags_ = new std::map<std::string, int>();
  tree->Branch("flags", "std::map<std::string, int>", &flags_);

  // add pre-firing weights
  tree->Branch("prefiringweight", &prefiringweight, "prefiringweight/F");
  tree->Branch("prefiringweightup", &prefiringweightup, "prefiringweightup/F");
  tree->Branch("prefiringweightdown", &prefiringweightdown, "prefiringweightdown/F");

  lumitree = FS->make<TTree>("AC1Blumi", "AC1Blumi", 1);
  lumitree->Branch("lumi_run", &lumi_run, "lumi_run/i");
  lumitree->Branch("lumi_block", &lumi_block, "lumi_block/i");

  if(cdata)
    {
      lumitree->Branch("lumi_value", &lumi_value, "lumi_value/F");
      lumitree->Branch("lumi_valueerr", &lumi_valueerr, "lumi_valueerr/F");
      lumitree->Branch("lumi_livefrac", &lumi_livefrac, "lumi_livefrac/F");
      lumitree->Branch("lumi_deadfrac", &lumi_deadfrac, "lumi_deadfrac/F");
      lumitree->Branch("lumi_quality", &lumi_quality, "lumi_quality/i");
      lumitree->Branch("lumi_eventsprocessed", &lumi_eventsprocessed, "lumi_eventsprocessed/i");
      lumitree->Branch("lumi_eventsfiltered", &lumi_eventsfiltered, "lumi_eventsfiltered/i");
      lumitree->Branch("lumi_hltprescaletable", &lumi_hltprescaletable, "lumi_hltprescaletable/i");
      lumitree->Branch("lumi_l1algoprescaletable", &lumi_l1algoprescaletable, "lumi_l1algoprescaletable/i");
      lumitree->Branch("lumi_l1techprescaletable", &lumi_l1techprescaletable, "lumi_l1techprescaletable/i");
    }

  runtree = FS->make<TTree>("AC1Brun", "AC1Brun", 1);
  runtree->Branch("run_number", &run_number, "run_number/i");
  runtree->Branch("run_hltcount", &run_hltcount, "run_hltcount/i");
  runtree->Branch("run_hltnames", "std::vector<std::string>", &run_hltnames);
  runtree->Branch("run_hltfilters", "std::vector<std::string>",&run_hltfilters);
  runtree->Branch("run_hltmufilters", "std::vector<std::string>", &run_hltmufilters);
  runtree->Branch("run_hltelectronfilters", "std::vector<std::string>", &run_hltelectronfilters);
  runtree->Branch("run_hlttaufilters", "std::vector<std::string>", &run_hlttaufilters);
  runtree->Branch("run_hltphotonfilters", "std::vector<std::string>", &run_hltphotonfilters);
  runtree->Branch("run_hltjetfilters", "std::vector<std::string>", &run_hltjetfilters);
  runtree->Branch("run_floattaudiscriminators", "std::vector<std::string>", &run_floattaudiscriminators);
  runtree->Branch("run_binarytaudiscriminators", "std::vector<std::string>", &run_binarytaudiscriminators);
  runtree->Branch("run_btagdiscriminators", "std::vector<std::string>", &run_btagdiscriminators);
  runtree->Branch("run_hltprescaletablescount", &run_hltprescaletablescount, "run_hltprescaletablescount/i");
  runtree->Branch("run_hltprescaletables", run_hltprescaletables, "run_hltprescaletables[run_hltprescaletablescount]/i");
  runtree->Branch("run_l1algocount", &run_l1algocount, "run_l1algocount/i");
  runtree->Branch("run_l1algoprescaletablescount", &run_l1algoprescaletablescount, "run_l1algoprescaletablescount/i");
  runtree->Branch("run_l1algoprescaletables", run_l1algoprescaletables, "run_l1algoprescaletables[run_l1algoprescaletablescount]/i");
  runtree->Branch("run_l1techcount", &run_l1techcount, "run_l1techcount/i");
  runtree->Branch("run_l1techprescaletablescount", &run_l1techprescaletablescount, "run_l1techprescaletablescount/i");
  runtree->Branch("run_l1techprescaletables", run_l1techprescaletables, "run_l1techprescaletables[run_l1techprescaletablescount]/i");




} //void NTupleMaker::beginJob()



//http://www.boost.org/doc/libs/1_55_0/libs/regex/doc/html/boost_regex/introduction_and_overview.html
void AddTriggerList(const edm::Run& run,
		    const HLTConfigProvider& hltConfig,
                    const std::vector<std::string> run_hltnames,
                    const std::vector<std::pair<boost::regex, std::string> >& regexes,
                    std::vector<std::pair<std::string, std::string> >& foundTriggers,
                    std::string& triggerNames   )
{
  for(std::size_t i = 0; i < hltConfig.size(); ++i) {

    // In some early 2011 runs saveTagsModules() does not give the
    // modules which actually save tags, but an empty list. So check all
    // the modules instead.

    //const std::vector<std::string> saveTagsModules = hltConfig.saveTagsModules(i);
    const std::vector<std::string> saveTagsModules = hltConfig.moduleLabels(i);

    for(std::size_t j = 0; j < regexes.size(); ++j) {
      boost::cmatch what;

      if(boost::regex_match(hltConfig.triggerName(i), regexes[j].first))
	{
	  // NEW: Continue if trigger was not specified in the HLTriggerPaths in the TreeProducer
	  // cout << "regexes[j].first = " << regexes[j].first << endl;
	  // cout<<"regexes[j].second = "<<regexes[j].second<<endl;
	  bool trigger_found = false;
	  for(auto const& s:run_hltnames){
	    if( boost::regex_match(s, regexes[j].first) ){
              trigger_found = true;
	    }
	  }
	  if(!trigger_found){
	    // cout<< "TRIGGER not found in HLTlist in TreeProducer: " << regexes[j].first << endl;
	    continue;
	  }

	  // Check for filter
	  if(regexes[j].second.size() != 0)
	      {
		// First, check for leg
		std::string::size_type legpos = regexes[j].second.find('|');
		std::string leg, filters;
		if(legpos != std::string::npos)
		  {
		    leg     = regexes[j].second.substr(0, legpos);
		    filters = regexes[j].second.substr(legpos + 1);
		  }
		else
		  {
		    filters = regexes[j].second;
		  }

		std::vector<std::string> strs;
		boost::split(strs, filters, boost::is_any_of(","));
		bool foundFilter = false;
		std::string filter;

		for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		  {
		    for(std::size_t l = 0; l < strs.size() && !foundFilter; ++l)
		      {
			if(saveTagsModules[k] == strs[l])
			  {
			    filter = saveTagsModules[k];
			    foundFilter = true;
			  }
		      }
		  }//for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)

		if(!foundFilter)
		  {
		    for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		      std::cout << saveTagsModules[k] << std::endl;

		    std::cout<<std::endl;
		    for(std::size_t l = 0; l < strs.size() && !foundFilter; ++l)
		      std::cout << strs[l] << std::endl;

		    throw cms::Exception("NTupleMaker") << "Did not find filter for trigger " << hltConfig.triggerName(i) << " in run " << run.run() << std::endl;
		  }

		if(leg.empty())
		  triggerNames += hltConfig.triggerName(i) + string(":") + filter + string(" ");
		else
		  triggerNames += hltConfig.triggerName(i) + string("|") + leg + string(":") + filter + string(" ");

		foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), filter));
	      }//if(regexes[j].second.size() != 0)

	  else
	    {
	      triggerNames += hltConfig.triggerName(i) + string(" ");
	      foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), saveTagsModules.back()));
	    }
	}
    }//for(std::size_t j = 0; j < regexes.size(); ++j)
  }//for(std::size_t i = 0; i < hltConfig.size(); ++i)
}// std::string& triggerNames


void NTupleMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  if(propagatorWithMaterial != 0){ delete propagatorWithMaterial;}
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  propagatorWithMaterial = new PropagatorWithMaterial(alongMomentum, 0.10566, &(*magneticField));
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", TTrackBuilder);

  run_number = iRun.run();
  //L1 prescales
  edm::ESHandle<L1GtPrescaleFactors> l1GtPfAlgo;
  iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().get(l1GtPfAlgo);

  unsigned numl1algo = (l1GtPfAlgo.product()->gtPrescaleFactors())[0].size();
  unsigned numl1algotables = (l1GtPfAlgo.product()->gtPrescaleFactors()).size();

  run_l1algoprescaletablescount = numl1algo*numl1algotables;
  run_l1algocount = numl1algo;
  if(l1GtPfAlgo.isValid())
    {
      for(unsigned i = 0 ; i < numl1algotables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1algo ; j++)
	    {
	      run_l1algoprescaletables[j+numl1algo*i] = (l1GtPfAlgo.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }

  edm::ESHandle<L1GtPrescaleFactors> l1GtPfTech;
  iSetup.get<L1GtPrescaleFactorsTechTrigRcd>().get(l1GtPfTech);

  unsigned numl1tech = (l1GtPfTech.product()->gtPrescaleFactors())[0].size();
  unsigned numl1techtables = (l1GtPfTech.product()->gtPrescaleFactors()).size();

  run_l1techprescaletablescount = numl1tech*numl1techtables;
  run_l1techcount = numl1tech;
  if(l1GtPfTech.isValid())
    {
      for(unsigned i = 0 ; i < numl1techtables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1tech ; j++)
	    {
	      run_l1techprescaletables[j+numl1tech*i] = (l1GtPfTech.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }

  //HLT names and prescales
  muontriggers.clear();
  electrontriggers.clear();
  tautriggers.clear();
  photontriggers.clear();
  jettriggers.clear();

  bool changed = true;
  HLTConfiguration.init(iRun, iSetup, cTriggerProcess, changed);
  HLTPrescaleConfig->init(iRun, iSetup, cTriggerProcess, changed);

  for(std::size_t i = 0; i < HLTConfiguration.size(); ++i) {
    TString TrigName(HLTConfiguration.triggerName(i));
    for (unsigned j = 0; j<cHLTriggerPaths.size(); ++j) {
      TString TrigNameConf(cHLTriggerPaths[j]);
      if (TrigName.Contains(TrigNameConf)) {
	//	std::cout << TrigName << " --> " << std::endl;
	const std::vector<std::string> saveTagsModules = HLTConfiguration.moduleLabels(i);
	for (unsigned k = 0; k<cHLTriggerPaths.size(); ++k) {
	  //	  std::cout << "    " << k << " " << saveTagsModules[k] << std::endl;
	}
      }
    }
  }

  vector<pair<boost::regex, string> > muonregexes;
  for(unsigned i = 0 ; i < cMuHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cMuHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      muonregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > electronregexes;
  for(unsigned i = 0 ; i < cElHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cElHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      electronregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > tauregexes;
  for(unsigned i = 0 ; i < cTauHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cTauHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      tauregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > photonregexes;
  for(unsigned i = 0 ; i < cPhotonHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cPhotonHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      photonregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > jetregexes;
  for(unsigned i = 0 ; i < cJetHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cJetHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      jetregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }

  run_hltcount = HLTConfiguration.size();
  string allnames;
  string allmuonnames;
  string allelectronnames;
  string alltaunames;
  string allphotonnames;
  string alljetnames;

  AddTriggerList(iRun, HLTConfiguration, cHLTriggerPaths, muonregexes,     muontriggers,     allmuonnames);
  AddTriggerList(iRun, HLTConfiguration, cHLTriggerPaths, electronregexes, electrontriggers, allelectronnames);
  AddTriggerList(iRun, HLTConfiguration, cHLTriggerPaths, tauregexes,      tautriggers,      alltaunames);
  AddTriggerList(iRun, HLTConfiguration, cHLTriggerPaths, photonregexes,   photontriggers,   allphotonnames);
  AddTriggerList(iRun, HLTConfiguration, cHLTriggerPaths, jetregexes,      jettriggers,      alljetnames);


  run_hltnames.clear();
  for (unsigned int i=0; i<cHLTriggerPaths.size(); ++i)
    run_hltnames.push_back(cHLTriggerPaths.at(i));

  if(muontriggers.size() > 100) throw cms::Exception("NTupleMaker") << "Too many muon triggers!" << std::endl;
  if(electrontriggers.size() > 100) throw cms::Exception("NTupleMaker") << "Too many electron triggers!" << std::endl;
  if(tautriggers.size() > 100) throw cms::Exception("NTupleMaker") << "Too many tau triggers!" << std::endl;
  if(photontriggers.size() > 100) throw cms::Exception("NTupleMaker") << "Too many photon triggers!" << std::endl;
  if(jettriggers.size() > 100) throw cms::Exception("NTupleMaker") << "Too many jet triggers!" << std::endl;

  // adding all filters
  run_hltfilters.clear();
  run_hltmufilters.clear();
  run_hltelectronfilters.clear();
  run_hlttaufilters.clear();
  run_hltphotonfilters.clear();
  run_hltjetfilters.clear();
  for (unsigned int i = 0; i<muontriggers.size(); ++i) {
    run_hltfilters.push_back(muontriggers.at(i).second);
    run_hltmufilters.push_back(muontriggers.at(i).second);
  }
  for (unsigned int i = 0; i<electrontriggers.size(); ++i) {
    run_hltfilters.push_back(electrontriggers.at(i).second);
    run_hltelectronfilters.push_back(electrontriggers.at(i).second);
  }
  for (unsigned int i = 0; i<tautriggers.size(); ++i) {
    run_hltfilters.push_back(tautriggers.at(i).second);
    run_hlttaufilters.push_back(tautriggers.at(i).second);
  }
  for (unsigned int i = 0; i<photontriggers.size(); ++i) {
    run_hltfilters.push_back(photontriggers.at(i).second);
    run_hltphotonfilters.push_back(photontriggers.at(i).second);
  }
  for (unsigned int i = 0; i<jettriggers.size(); ++i) {
    run_hltfilters.push_back(jettriggers.at(i).second);
    run_hltjetfilters.push_back(jettriggers.at(i).second);
  }

  if (run_hltfilters.size()>200) throw cms::Exception("NTupleMaker") << "Too many HLT filters!" << std::endl;

  run_hltprescaletablescount = HLTConfiguration.prescaleSize()*HLTConfiguration.size();

  // adding btag discriminators
  run_btagdiscriminators.clear();
  for (unsigned int i=0; i<cBtagDiscriminators.size(); ++i)
    run_btagdiscriminators.push_back(cBtagDiscriminators.at(i));

  if (run_btagdiscriminators.size()>10) throw cms::Exception("NTupleMaker") << "Too many btag discriminators!" << std::endl;

  // adding tau  discriminators
  run_floattaudiscriminators.clear();
  for (unsigned int i=0; i<cTauFloatDiscriminators.size(); ++i)
    run_floattaudiscriminators.push_back(cTauFloatDiscriminators.at(i));

  run_binarytaudiscriminators.clear();
  for (unsigned int i=0; i<cTauBinaryDiscriminators.size(); ++i)
    run_binarytaudiscriminators.push_back(cTauBinaryDiscriminators.at(i));

  if (run_floattaudiscriminators.size()>10) throw cms::Exception("NTupleMaker") << "Too many float tau discriminators!" << std::endl;
  if (run_binarytaudiscriminators.size()>50) throw cms::Exception("NTupleMaker") << "Too many binary tau discriminators!" << std::endl;

  for(unsigned j = 0 ; j < HLTConfiguration.prescaleSize() ; j++)
    {
      for(unsigned i = 0 ; i < HLTConfiguration.size() ; i++)
	{
	  run_hltprescaletables[i+HLTConfiguration.size()*j] = HLTConfiguration.prescaleValue(j, HLTConfiguration.triggerName(i));
	}
    }
  runtree->Fill();

  // JEC
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  jecUnc = new JetCorrectionUncertainty(JetCorPar);

  /*
  std::cout << "Begin of run ----> " << std::endl;

  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  iRun.getByLabel( "externalLHEProducer", run );
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());

  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
  */

}

void NTupleMaker::endRun(const edm::Run& iRun)
{
  delete jecUnc;

  std::cout << "End of run --->" << std::endl;

  // LHERunInfoProduct

  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  iRun.getByLabel( "externalLHEProducer", run );
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());

  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }



}

void NTupleMaker::beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumi_run = iLumiBlock.run();
  lumi_block = iLumiBlock.luminosityBlock();

  if(cdata)
    {
      //edm::Handle<LumiSummary> lumiSummary;
      //iLumiBlock.getByLabel(edm::InputTag("lumiProducer"), lumiSummary);
      lumi_value = 0;//lumiSummary->avgInsDelLumi();
      lumi_valueerr = 0;//lumiSummary->avgInsDelLumiErr();
      lumi_livefrac = 0;//lumiSummary->lumiSecQual();
      lumi_deadfrac = 0;//lumiSummary->deadFrac();
      lumi_quality = 0;//lumiSummary->liveFrac();
      lumi_eventsprocessed = 0;
      lumi_eventsfiltered = 0;
    }
}

void NTupleMaker::endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumitree->Fill();
}


void NTupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug)  cout<<"inside the analyze function"<< endl;

  track_count = 0;
  goodprimvertex_count = 0;
  primvertex_count = 0;
  refitvertex_count = 0;
  refitvertexwithbs_count = 0;
  muon_count = 0;
  //dimuon_count = 0;
  tau_count = 0;
  l1muon_count = 0;
  l1egamma_count = 0;
  l1tau_count = 0;
  l1isotau_count = 0;
  gentau_count = 0;
  pfjet_count = 0;
  pfjetpuppi_count = 0;
  electron_count = 0;
  photon_count = 0;
  genparticles_count = 0;
  genjets_count = 0;
  errors = 0;
  trigobject_count = 0;
  mvamet_count = 0;

  bool takeevent = true;

  nEvents->Fill(0);
  pv_position = math::XYZPoint(0.,0.,0.);

  lumi_eventsprocessed++;

  event_nr      = iEvent.id().event();
  event_run      = iEvent.id().run();
  event_timeunix = iEvent.time().unixTime();
  event_timemicrosec = iEvent.time().microsecondOffset();
  event_luminosityblock = iEvent.getLuminosityBlock().luminosityBlock();

  //  cout << "NTupleMaker.cc : run = " << event_run << "  event = " << event_nr << std::endl;

  // L1TriggerBits
  // https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h?v=CMSSW_6_2_0_SLHC2#04
  edm::Handle<L1GlobalTriggerReadoutRecord> L1trigger;
  iEvent.getByLabel(edm::InputTag("gtDigis"), L1trigger);
  if (L1trigger.isValid()){

    const TechnicalTriggerWord& L1triggerbits = L1trigger->technicalTriggerWord();
    for(int i  = 0  ; i < 8 ; i++) trigger_level1bits[i] = 0;

    for(unsigned i = 0 ; i < min(unsigned(L1triggerbits.size()), unsigned(64)) ; i++)
      trigger_level1bits[i/8] |= (Byte_t)L1triggerbits[i] << (i % 8);  // bitwise OR -> |

    //trigger_level1bits[i/8] = trigger_level1bits[i/8]  | (Byte_t) L1triggerbits[i] << (i % 8);  // bitwise OR -> |


    //L1TriggerAlgos
    const DecisionWord& L1triggeralgos = L1trigger->decisionWord();
    for(int i = 0  ; i < 128 ; i++){trigger_level1[i] = 0;}
    for(unsigned i = 0 ; i < min(unsigned(L1triggeralgos.size()), unsigned(1024)) ; i++)
      {
	trigger_level1[i/8] |= (Byte_t)L1triggeralgos[i] << (i%8);
      }
    lumi_l1techprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexTech();
    lumi_l1algoprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexAlgo();
  }
  lumi_hltprescaletable = -1;

  //HLTriggerResults
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", cTriggerProcess), HLTrigger);
  hltriggerresults_->clear();
  hltriggerprescales_->clear();
  hltriggerresultsV_.clear();
  if(HLTrigger.isValid()){
    //    std::cout << "Valid HLT" << std::endl;
    for(int i = 0  ; i < 128 ; i++){trigger_HLT[i] = 0;}
    //store trigger bits for selected trigger paths
    const edm::TriggerNames& TrigNames_ = iEvent.triggerNames(*HLTrigger);
    for(unsigned i = 0 ; i < HLTrigger->size(); i++)
      {
	if(!HLTrigger->wasrun(i) )continue;
	std::string trigName=TrigNames_.triggerName(i);
	if(cHLTriggerPaths.size() > 0){
	  for(size_t ip = 0; ip < cHLTriggerPaths.size(); ip++){
	    if(trigName.find(cHLTriggerPaths[ip]) != string::npos){

	      //hltriggerprescales_->insert(std::pair<string, int>(trigName, 1.));// FIXME HLTPrescaleConfig->prescaleValue(iEvent,iSetup,trigName)));
	      hltriggerprescales_->insert(std::pair<string, int>(trigName,
								 int(HLTPrescaleConfig->prescaleValue(iEvent,iSetup,trigName))));
	      hltriggerresults_->insert(std::pair<string, int>(trigName, HLTrigger->accept(i)));
	      TString TriggerName(trigName);
	      //	      std::cout << trigName << " : "
	      //			<< HLTrigger->accept(i)
	      //			<< " ; prescale : "  << HLTConfiguration.prescaleValue(iEvent,iSetup,trigName)
	      //		<< std::endl;
	      if(HLTrigger->accept(i)) hltriggerresultsV_.push_back(trigName);
	    }
	  }
	}
      }
  }

  if (!cFastSim){
    flags_->clear();
    for(std::vector<string>::iterator it = cFlagsProcesses.begin();
	it != cFlagsProcesses.end(); it++){
      //std::cout<<it->data()<<std::endl;
      AddFlags(iEvent,"TriggerResults", "", it->data());
    }

    /*
    edm::Handle<bool> ifilterbadChCand;
    iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
    flags_->insert(std::pair<string, int>("Flag_BadChargedCandidateFilter", *ifilterbadChCand));

    edm::Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
    flags_->insert(std::pair<string, int>("Flag_BadPFMuonFilter", *ifilterbadPFMuon));

    // Bad global muons
    bool ifilterBadGlobalMuon;
    edm::Handle<edm::PtrVector<reco::Muon>> BadGlobalMuons;
    iEvent.getByToken(BadGlobalMuonsToken_, BadGlobalMuons);
    if(BadGlobalMuons->size() == 0 ) ifilterBadGlobalMuon = true;
    else                             ifilterBadGlobalMuon = false;
    flags_->insert(std::pair<string, int>("Flag_BadGlobalMuonFilter", ifilterBadGlobalMuon));
    */
  }
  //int _passecalBadCalibFilterUpdate;// = 0;

  if(cYear == 2018){
    edm::Handle< bool > passecalBadCalibFilterUpdate ;
    iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
    _passecalBadCalibFilterUpdate = *passecalBadCalibFilterUpdate;
  }

  if(cbeamspot)
    {
      edm::Handle<BeamSpot> TheBeamSpot;
      iEvent.getByToken(BeamSpotToken_, TheBeamSpot);
      if(TheBeamSpot.isValid())
	{
	  beamspot_x = TheBeamSpot->x0();
	  beamspot_y = TheBeamSpot->y0();
	  beamspot_z = TheBeamSpot->z0();
	  beamspot_xwidth = TheBeamSpot->BeamWidthX();
	  beamspot_ywidth = TheBeamSpot->BeamWidthY();
	  beamspot_zsigma = TheBeamSpot->sigmaZ();
	  beamspot_cov[0] = TheBeamSpot->covariance(0,0);
	  beamspot_cov[1] = TheBeamSpot->covariance(0,1);
	  beamspot_cov[2] = TheBeamSpot->covariance(0,2);
	  beamspot_cov[3] = TheBeamSpot->covariance(1,1);
	  beamspot_cov[4] = TheBeamSpot->covariance(1,2);
	  beamspot_cov[5] = TheBeamSpot->covariance(2,2);
	  pv_position = math::XYZPoint(TheBeamSpot->x0(), TheBeamSpot->y0(), TheBeamSpot->z0());
	}
      else
	{
	  beamspot_x = 0.;
	  beamspot_y = 0.;
	  beamspot_z = 0.;
	  beamspot_xwidth = 0.;
	  beamspot_ywidth = 0.;
	  beamspot_zsigma = 0.;
	  beamspot_cov[0] = 0.;
	  beamspot_cov[1] = 0.;
	  beamspot_cov[2] = 0.;
	  beamspot_cov[3] = 0.;
	  beamspot_cov[4] = 0.;
	  beamspot_cov[5] = 0.;
	}
    }

  if(crecprimvertex)
    {
      edm::Handle<VertexCollection> Vertex;
      iEvent.getByToken(PVToken_, Vertex);
      if(Vertex.isValid()) {
	primvertex_mindz = 999;
	for(unsigned i = 0 ; i < Vertex->size(); i++) {
	  primvertex_count++;
	  if(i == 0) {
	    primvertex_x = (*Vertex)[i].x();
	    primvertex_y = (*Vertex)[i].y();
	    primvertex_z = (*Vertex)[i].z();
	    primvertex_chi2 = (*Vertex)[i].chi2();
	    primvertex_ndof = (*Vertex)[i].ndof();
	    primvertex_ntracks = (*Vertex)[i].tracksSize();
	    primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
	    primvertex_cov[1] = (*Vertex)[i].covariance(0,1);
	    primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
	    primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
	    primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
	    primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()
	    Float_t ptq = 0.;
	    for(Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin() ; it != (*Vertex)[i].tracks_end() ; ++it)
	      {
		ptq += (*it)->pt() * (*it)->pt();
	      }
	    primvertex_ptq = ptq;

	    pv_position = (*Vertex)[i].position();
	    primvertex = (*Vertex)[i];
	  } else {
	    if(std::abs((*Vertex)[i].z()-(*Vertex)[0].z()) < primvertex_mindz)
	      primvertex_mindz = std::abs((*Vertex)[i].z()-(*Vertex)[0].z()); //minimal longitudinal distance between the PV and other vertex
	  }
	  if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake() && (*Vertex)[i].ndof() >= 4 && (*Vertex)[i].z() > -24 && (*Vertex)[i].z() < 24 && (*Vertex)[i].position().Rho() < 2.)
	    goodprimvertex_count++;
	}
      }
    }

  if(crecprimvertexwithbs)
    {
      edm::Handle<RefitVertexCollection> VertexWithBS;
      iEvent.getByToken(PVwithBSToken_, VertexWithBS);
      if(VertexWithBS.isValid()) {
        for(unsigned i = 0 ; i < VertexWithBS->size(); i++) {
          if(i == 0) {
            primvertexwithbs_x = (*VertexWithBS)[i].x();
            primvertexwithbs_y = (*VertexWithBS)[i].y();
            primvertexwithbs_z = (*VertexWithBS)[i].z();
            primvertexwithbs_chi2 = (*VertexWithBS)[i].chi2();
            primvertexwithbs_ndof = (*VertexWithBS)[i].ndof();
            primvertexwithbs_ntracks = (*VertexWithBS)[i].tracksSize();
            primvertexwithbs_cov[0] = (*VertexWithBS)[i].covariance(0,0); // xError()
            primvertexwithbs_cov[1] = (*VertexWithBS)[i].covariance(0,1);
            primvertexwithbs_cov[2] = (*VertexWithBS)[i].covariance(0,2);
            primvertexwithbs_cov[3] = (*VertexWithBS)[i].covariance(1,1); // yError()
            primvertexwithbs_cov[4] = (*VertexWithBS)[i].covariance(1,2);
            primvertexwithbs_cov[5] = (*VertexWithBS)[i].covariance(2,2); // zError()
	  }
	}
      }
    }
  if (csusyinfo) {
    if(doDebug)  cout<<"add SUSY info"<< endl;
    AddSusyInfo(iEvent);
    //    std::cout << "Run = " << event_run << "   Lumi = " << event_luminosityblock << "   Event = " << event_nr << std::endl;
    //    std::cout << "SUSY Mother Mass = " << SusyMotherMass << std::endl;
    //    std::cout << "SUSY LSP Mass = " << SusyLSPMass << std::endl;
    //    std::cout << std::endl;
  }


  if (crecmuon)
    {
      if(doDebug)  cout<<"add muons"<< endl;
      int numberOfMuons = int(AddMuons(iEvent, iSetup));
    } // crecmuon

  if (crecelectron)
    {
      if(doDebug)  cout<<"add electrons"<< endl;
      int numberOfElectrons = int(AddElectrons(iEvent,iSetup));
    } // crecelectron

  if(crectau)
    {
      if (doDebug) cout<<"add taus"<< endl;
      int numberOfTaus = int(AddTaus(iEvent, iSetup));
      // if (cSkim>0) {
      // 	bool goodTaus = false;
      // 	for (unsigned int i=0; i<tau_count; ++i) {
      // 	  if (tau_pt[i]>50. &&
      // 	      fabs(tau_eta[i])<2.3 &&
      // 	      ( tau_decayModeFinding[i]>0.5 || tau_decayModeFindingNewDMs[i]>0.5) &&
      // 	      tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[i]>0.5) {
      // 	    goodTaus = true;
	    // std::cout << "Good tau : pt = " << tau_pt[i]
	    //  	      << "  eta = " << tau_eta[i]
	    // 	      << "  phi = " << tau_phi[i]
	    //  	      << "  decay = " <<  tau_decayMode_name[i] << std::endl;
      // 	    break;
      // 	  }
      // 	}
      // 	if (!goodTaus) return;
      // }
    }

  //add refitted vertex with tracks removed from pair of leptons
  if (crefittedvertex)
    {
      edm::Handle<RefitVertexCollection> RFVertex;
      iEvent.getByToken(RefittedPVToken_, RFVertex);
      if(RFVertex.isValid()) {
        for(unsigned i = 0 ; i < RFVertex->size(); i++) {

          //Find pair of leptons
          const std::vector<std::string> leptonPairNames = (*RFVertex)[i].userCandNames();
          if(leptonPairNames.size() < 2) continue;
          int nMu = 0, nEle = 0, nTau = 0;
          int muon1 = -1, muon2 = -1, ele1 = -1, ele2 = -1, tau1 = -1, tau2 = -1;
          for(size_t il = 0; il < leptonPairNames.size(); il++){
            const std::string leptonName = leptonPairNames[il];
	    edm::Ptr<reco::Candidate> leptonCand = (*RFVertex)[i].userCand(leptonName);
            if(std::abs(leptonCand->pdgId())==11 ){
              for(size_t ie = 0; ie < electron_count; ie++){
                if(TMath::Abs(leptonCand->pt() - electron_pt[ie]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - electron_eta[ie]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - electron_phi[ie]) < 0.001){
                  nEle++;
                  if(nEle == 1)ele1 = ie;
                  else if(nEle == 2) ele2 = ie;
                }
              }
            }
            else if(std::abs(leptonCand->pdgId())==13 ){
              for(size_t im = 0; im < muon_count; im++){
                if(TMath::Abs(leptonCand->pt() - muon_pt[im]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - muon_eta[im]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - muon_phi[im]) < 0.001){
                  nMu++;
                  if(nMu == 1)muon1 = im;
                  else if(nMu == 2) muon2 = im;
                }
              }
            }
            else { //for tau
              for(size_t it = 0; it < tau_count; it++){
                if(TMath::Abs(leptonCand->pt() - tau_pt[it]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - tau_eta[it]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - tau_phi[it]) < 0.001){
                  nTau++;
                  if(nTau == 1)tau1 = it;
                  else if(nTau == 2) tau2 = it;
                }
              }
            }
          } //end of lepton pair
          if((nEle + nMu + nTau) < 2) continue;  //should have a pair of selected leptons

          refitvertex_eleIndex[refitvertex_count][0] = ele1;
          refitvertex_eleIndex[refitvertex_count][1] = ele2;
          refitvertex_muIndex[refitvertex_count][0] = muon1;
          refitvertex_muIndex[refitvertex_count][1] = muon2;
          refitvertex_tauIndex[refitvertex_count][0] = tau1;
          refitvertex_tauIndex[refitvertex_count][1] = tau2;
          refitvertex_x[refitvertex_count] = (*RFVertex)[i].x();
          refitvertex_y[refitvertex_count] = (*RFVertex)[i].y();
          refitvertex_z[refitvertex_count] = (*RFVertex)[i].z();
          refitvertex_chi2[refitvertex_count] = (*RFVertex)[i].chi2();
          refitvertex_ndof[refitvertex_count] = (*RFVertex)[i].ndof();
          refitvertex_ntracks[refitvertex_count] = (*RFVertex)[i].tracksSize();
          refitvertex_cov[refitvertex_count][0] = (*RFVertex)[i].covariance(0,0); // xError()
          refitvertex_cov[refitvertex_count][1] = (*RFVertex)[i].covariance(0,1);
          refitvertex_cov[refitvertex_count][2] = (*RFVertex)[i].covariance(0,2);
          refitvertex_cov[refitvertex_count][3] = (*RFVertex)[i].covariance(1,1); // yError()
          refitvertex_cov[refitvertex_count][4] = (*RFVertex)[i].covariance(1,2);
          refitvertex_cov[refitvertex_count][5] = (*RFVertex)[i].covariance(2,2); // zError()
          refitvertex_count++;
        }
      }
    }

  //add refitted vertex with tracks removed from pair of leptons along with bs constraint
  if (crefittedvertexwithbs)
    {
      edm::Handle<RefitVertexCollection> RFVertexwithbs;
      iEvent.getByToken(RefittedwithBSPVToken_, RFVertexwithbs);
      if(RFVertexwithbs.isValid()) {
        for(unsigned i = 0 ; i < RFVertexwithbs->size(); i++) {

          //Find pair of leptons
          const std::vector<std::string> leptonPairNames = (*RFVertexwithbs)[i].userCandNames();
          if(leptonPairNames.size() < 2) continue;
          int nMu = 0, nEle = 0, nTau = 0;
          int muon1 = -1, muon2 = -1, ele1 = -1, ele2 = -1, tau1 = -1, tau2 = -1;
          for(size_t il = 0; il < leptonPairNames.size(); il++){
            const std::string leptonName = leptonPairNames[il];
	    edm::Ptr<reco::Candidate> leptonCand = (*RFVertexwithbs)[i].userCand(leptonName);
            if(std::abs(leptonCand->pdgId())==11 ){
              for(size_t ie = 0; ie < electron_count; ie++){
                if(TMath::Abs(leptonCand->pt() - electron_pt[ie]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - electron_eta[ie]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - electron_phi[ie]) < 0.001){
                  nEle++;
                  if(nEle == 1)ele1 = ie;
                  else if(nEle == 2) ele2 = ie;
                }
              }
            }
            else if(std::abs(leptonCand->pdgId())==13 ){
              for(size_t im = 0; im < muon_count; im++){
                if(TMath::Abs(leptonCand->pt() - muon_pt[im]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - muon_eta[im]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - muon_phi[im]) < 0.001){
                  nMu++;
                  if(nMu == 1)muon1 = im;
                  else if(nMu == 2) muon2 = im;
                }
              }
            }
            else { //for tau
              for(size_t it = 0; it < tau_count; it++){
                if(TMath::Abs(leptonCand->pt() - tau_pt[it]) < 0.01 &&
                   TMath::Abs(leptonCand->eta() - tau_eta[it]) < 0.001 &&
                   TMath::Abs(leptonCand->phi() - tau_phi[it]) < 0.001){
                  nTau++;
                  if(nTau == 1)tau1 = it;
                  else if(nTau == 2) tau2 = it;
                }
              }
            }
          } //end of lepton pair
          if((nEle + nMu + nTau) < 2) continue;  //should have a pair of selected leptons

          refitvertexwithbs_eleIndex[refitvertexwithbs_count][0] = ele1;
          refitvertexwithbs_eleIndex[refitvertexwithbs_count][1] = ele2;
          refitvertexwithbs_muIndex[refitvertexwithbs_count][0] = muon1;
          refitvertexwithbs_muIndex[refitvertexwithbs_count][1] = muon2;
          refitvertexwithbs_tauIndex[refitvertexwithbs_count][0] = tau1;
          refitvertexwithbs_tauIndex[refitvertexwithbs_count][1] = tau2;
          refitvertexwithbs_x[refitvertexwithbs_count] = (*RFVertexwithbs)[i].x();
          refitvertexwithbs_y[refitvertexwithbs_count] = (*RFVertexwithbs)[i].y();
          refitvertexwithbs_z[refitvertexwithbs_count] = (*RFVertexwithbs)[i].z();
          refitvertexwithbs_chi2[refitvertexwithbs_count] = (*RFVertexwithbs)[i].chi2();
          refitvertexwithbs_ndof[refitvertexwithbs_count] = (*RFVertexwithbs)[i].ndof();
          refitvertexwithbs_ntracks[refitvertexwithbs_count] = (*RFVertexwithbs)[i].tracksSize();
          refitvertexwithbs_cov[refitvertexwithbs_count][0] = (*RFVertexwithbs)[i].covariance(0,0); // xError()
          refitvertexwithbs_cov[refitvertexwithbs_count][1] = (*RFVertexwithbs)[i].covariance(0,1);
          refitvertexwithbs_cov[refitvertexwithbs_count][2] = (*RFVertexwithbs)[i].covariance(0,2);
          refitvertexwithbs_cov[refitvertexwithbs_count][3] = (*RFVertexwithbs)[i].covariance(1,1); // yError()
          refitvertexwithbs_cov[refitvertexwithbs_count][4] = (*RFVertexwithbs)[i].covariance(1,2);
          refitvertexwithbs_cov[refitvertexwithbs_count][5] = (*RFVertexwithbs)[i].covariance(2,2); // zError()
          refitvertexwithbs_count++;
        }
      }
    }


  if (crectrack) AddPFCand(iEvent, iSetup);

  edm::Handle<BXVector<l1t::Muon> > l1muons;
  iEvent.getByToken( L1MuonCollectionToken_, l1muons);

  edm::Handle<BXVector<l1t::EGamma> > l1egammas;
  iEvent.getByToken( L1EGammaCollectionToken_, l1egammas);

  edm::Handle<BXVector<l1t::Tau> > l1taus;
  iEvent.getByToken( L1TauCollectionToken_, l1taus);

  edm::Handle<BXVector<l1t::Tau> > l1isotaus;
  iEvent.getByToken( L1TauCollectionToken_, l1isotaus);

  if (cl1objects && l1muons.isValid() && l1egammas.isValid() && l1taus.isValid() && l1isotaus.isValid()){
    if (doDebug) cout<<"add L1 Objects"<< endl;

    for(int ibx = l1muons->getFirstBX() ; ibx <= l1muons->getLastBX() ; ++ibx) {
      for(BXVector<l1t::Muon>::const_iterator it=l1muons->begin(ibx); it!=l1muons->end(ibx); it++) {
	if(l1muon_count == M_muonmaxcount) {
	  cerr << "number of L1 Muons > M_muonmaxcount. They are missing." << endl;
	  errors |= 1<<3;
	  break;
	}

	if(it->pt() < cMuPtMin ) continue;

	l1muon_px[l1muon_count]       = it->px();
	l1muon_py[l1muon_count]       = it->py();
	l1muon_pz[l1muon_count]       = it->pz();
	l1muon_pt[l1muon_count]       = it->pt();
//	l1muon_ipt[l1muon_count]      = it->hwPt();
	l1muon_eta[l1muon_count]      = it->hwEta();
	l1muon_phi[l1muon_count]      = it->hwPhi();
	l1muon_iso[l1muon_count]      = it->hwIso();
//	l1muon_qual[l1muon_count]     = it->hwQual();

	l1muon_charge[l1muon_count]      = it->hwCharge();
//	l1muon_chargeValid[l1muon_count] = it->hwChargeValid();
	l1muon_muonIndex[l1muon_count]   = it->tfMuonIndex();
//	l1muon_tag[l1muon_count]         = it->hwTag();
//	l1muon_isoSum[l1muon_count]      = it->hwIsoSum();
//	l1muon_dPhiExtra[l1muon_count]   = it->hwDPhiExtra();
//	l1muon_dEtaExtra[l1muon_count]   = it->hwDEtaExtra();
//	l1muon_rank[l1muon_count]        = it->hwRank();
//	l1muon_bx[l1muon_count]          = ibx;

	l1muon_count++;
      }
    }


    for(int ibx = l1egammas->getFirstBX() ; ibx <= l1egammas->getLastBX() ; ibx++) {
      for(BXVector<l1t::EGamma>::const_iterator it=l1egammas->begin(ibx); it!=l1egammas->end(ibx); it++) {
	if(l1egamma_count == M_electronmaxcount) {
	  cerr << "number of L1 Electrons > M_electronmaxcount. They are missing." << endl;
	  errors |= 1<<3;
	  break;
	}

	if(it->pt() < cElPtMin ) continue;

	l1egamma_px[l1egamma_count]       = it->px();
	l1egamma_py[l1egamma_count]       = it->py();
	l1egamma_pz[l1egamma_count]       = it->pz();
	l1egamma_pt[l1egamma_count]       = it->pt();
//	l1egamma_ipt[l1egamma_count]      = it->hwPt();
	l1egamma_eta[l1egamma_count]      = it->hwEta();
	l1egamma_phi[l1egamma_count]      = it->hwPhi();
	l1egamma_iso[l1egamma_count]      = it->hwIso();
//	l1egamma_qual[l1egamma_count]     = it->hwQual();
//	l1egamma_bx[l1egamma_count]       = ibx;

	//l1egamma_towerIEta[l1egamma_count]   = it->towerIEta();
	//l1egamma_towerIPhi[l1egamma_count]   = it->towerIPhi();
	//l1egamma_rawEt[l1egamma_count]       = it->rawEt();
	//l1egamma_isoEt[l1egamma_count]       = it->isoEt();
	//l1egamma_footprintEt[l1egamma_count] = it->footprintEt();
	//l1egamma_nTT[l1egamma_count]         = it->nTT();
	//l1egamma_shape[l1egamma_count]       = it->shape();

	  l1egamma_count++;
      }
    }

    for (int ibx = l1taus->getFirstBX(); ibx<=l1taus->getLastBX(); ++ibx) {
      for(BXVector<l1t::Tau>::const_iterator it=l1taus->begin(ibx); it!=l1taus->end(ibx); it++) {
	if(l1tau_count == M_taumaxcount) {
	  cerr << "number of L1 Taus > M_taumaxcount. They are missing." << endl;
	  errors |= 1<<3;
	  break;
	}

	if(it->pt() < cTauPtMin ) continue;

	l1tau_px[l1tau_count]       = it->px();
	l1tau_py[l1tau_count]       = it->py();
	l1tau_pz[l1tau_count]       = it->pz();
	l1tau_pt[l1tau_count]       = it->pt();
	l1tau_ipt[l1tau_count]      = it->hwPt();
	l1tau_eta[l1tau_count]      = it->hwEta();
	l1tau_phi[l1tau_count]      = it->hwPhi();
	l1tau_iso[l1tau_count]      = it->hwIso();
	l1tau_qual[l1tau_count]     = it->hwQual();

	l1tau_towerIEta[l1tau_count]  = it->towerIEta();
	l1tau_towerIPhi[l1tau_count]  = it->towerIPhi();
	l1tau_rawEt[l1tau_count]      = it->rawEt();
	l1tau_isoEt[l1tau_count]      = it->isoEt();
	l1tau_nTT[l1tau_count]        = it->nTT();
	l1tau_hasEM[l1tau_count]      = it->hasEM();
	l1tau_isMerged[l1tau_count]   = it->isMerged();
	l1tau_bx[l1tau_count]         = ibx;

	//	std::cout << "l1 tau : pT = " << l1tau_pt[l1tau_count] << "  BX = " << ibx << std::endl;

	l1tau_count++;
      }
    }
    //    std::cout << std::endl;

    /*
    for(unsigned itau = 0 ; itau < l1isotaus->size() ; itau++) {
      if(l1isotau_count == M_taumaxcount) {
	cerr << "number of iso taus > M_taumaxcount. They are missing." << endl;
	errors |= 1<<3;
	break;
      }
      if((*l1isotaus)[itau].pt() < cTauPtMin ) continue;
      l1isotau_e[l1isotau_count]        = (*l1isotaus)[itau].energy();
      l1isotau_px[l1isotau_count]       = (*l1isotaus)[itau].px();
      l1isotau_py[l1isotau_count]       = (*l1isotaus)[itau].py();
      l1isotau_pz[l1isotau_count]       = (*l1isotaus)[itau].pz();
      l1isotau_pt[l1isotau_count]       = (*l1isotaus)[itau].pt();
      l1isotau_phi[l1isotau_count]      = (*l1isotaus)[itau].phi();
      l1isotau_eta[l1isotau_count]      = (*l1isotaus)[itau].eta();
      l1isotau_mass[l1isotau_count]     = (*l1isotaus)[itau].mass();
      l1isotau_charge[l1isotau_count]   = (*l1isotaus)[itau].charge();
      l1isotau_iso[l1isotau_count]         = (*l1isotaus)[itau].hwIso();

      l1isotau_count++;
    }
    */
  }

  //l1EGammas, l1EGammaLabel  = Handle("BXVector<l1t::EGamma>"), "caloStage2Digis:EGamma";

  if (crecpfjet)
    {
      if(doDebug)  cout<<"add PF jets"<< endl;
      int numberOfJets = int(AddPFJets(iEvent,iSetup));
      // if (cSkim) {
      // 	bool goodJets = false;
      // 	for (unsigned int i=0; i<pfjet_count; ++i) {
      // 	  if (pfjet_pt[i]>50. && fabs(pfjet_eta[i])<2.3) {
      // 	    goodJets = true;
      // 	    // std::cout << "Good jet : pt = " << pfjet_pt[i]
      // 	    //  	      << "  eta = " << pfjet_eta[i]
      // 	    // 	      << "  phi = " << pfjet_phi[i] << std::endl;
      // 	    break;
      // 	  }
      // 	}
      // 	if (!goodJets) return;
      // }
    } // crecpfjet

    if (crecpfpuppijet)
      {
        if(doDebug)  cout<<"add PF Puppi jets"<< endl;
        int numberOfJetsPuppi = int(AddPFPuppiJets(iEvent,iSetup));
      }


  if(crecpfmet)
    {
      if(doDebug)  cout<<"add PF MET"<< endl;
      edm::Handle<pat::METCollection> patMet;
      iEvent.getByToken(MetCollectionToken_, patMet);

      assert(patMet->size() > 0);
      pfmet_ex = (*patMet)[0].px();
      pfmet_ey = (*patMet)[0].py();
      pfmet_ez = (*patMet)[0].pz();
      pfmet_pt = (*patMet)[0].pt();
      pfmet_phi = (*patMet)[0].phi();

      pfmet_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
      pfmet_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
      pfmet_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
      pfmet_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);
      pfmet_sig   = (*patMet)[0].significance();

      pfmet_ex_JetEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmet_ey_JetEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmet_ex_JetEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmet_ey_JetEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);

      pfmet_ex_UnclusteredEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmet_ey_UnclusteredEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmet_ex_UnclusteredEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmet_ey_UnclusteredEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);

      if (!cdata) {
	const reco::GenMET * genMET = (*patMet)[0].genMET();
	genmet_ex = genMET->px();
	genmet_ey = genMET->py();
      }
    } // crecpfmet

  if(crecpfmetcorr)
    {
      if(doDebug)  cout<<"add Corrected PF MET"<< endl;
      edm::Handle<pat::METCollection> patMet;
      iEvent.getByToken(MetCorrCollectionToken_, patMet);

      assert(patMet->size() > 0);
      pfmetcorr_ex = (*patMet)[0].px();
      pfmetcorr_ey = (*patMet)[0].py();
      pfmetcorr_ez = (*patMet)[0].pz();
      pfmetcorr_pt = (*patMet)[0].pt();
      pfmetcorr_phi = (*patMet)[0].phi();

      pfmetcorr_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
      pfmetcorr_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
      pfmetcorr_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
      pfmetcorr_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);
      pfmetcorr_sig   = (*patMet)[0].significance();

      pfmetcorr_ex_JetEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_JetEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_UnclusteredEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_UnclusteredEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_UnclusteredEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_UnclusteredEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_JetResUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetResUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_JetResDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetResDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);


  /*    pfmetcorr_ex_smeared = (*patMet)[0].corPx(pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_smeared = (*patMet)[0].corPy(pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_pt_smeared = (*patMet)[0].corPt(pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_phi_smeared = (*patMet)[0].corPhi(pat::MET::METCorrectionLevel::Type1Smear);


      pfmetcorr_ex_JetEnUp_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_JetEnUp_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1Smear);

      pfmetcorr_ex_JetEnDown_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_JetEnDown_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1Smear);

      pfmetcorr_ex_UnclusteredEnUp_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_UnclusteredEnUp_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1Smear);

      pfmetcorr_ex_UnclusteredEnDown_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_UnclusteredEnDown_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1Smear);

      pfmetcorr_ex_JetResUp_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_JetResUp_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1Smear);

      pfmetcorr_ex_JetResDown_smeared = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1Smear);
      pfmetcorr_ey_JetResDown_smeared = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1Smear); */
    } // crecpfmetcorr

  if(crecpuppimet)
    {
      if(doDebug)  cout<<"add Puppi MET"<< endl;
      edm::Handle<pat::METCollection> patMet;
      iEvent.getByToken(PuppiMetCollectionToken_, patMet);

      assert(patMet->size() > 0);
      puppimet_ex = (*patMet)[0].px();
      puppimet_ey = (*patMet)[0].py();
      puppimet_ez = (*patMet)[0].pz();
      puppimet_pt = (*patMet)[0].pt();
      puppimet_phi = (*patMet)[0].phi();

      puppimet_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
      puppimet_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
      puppimet_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
      puppimet_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);

      puppimet_ex_JetEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_JetEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);

      puppimet_ex_JetEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_JetEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);

      puppimet_ex_UnclusteredEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_UnclusteredEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);

      puppimet_ex_UnclusteredEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_UnclusteredEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);

      puppimet_ex_JetResUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_JetResUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);

      puppimet_ex_JetResDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);
      puppimet_ey_JetResDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);
    } // crecpuppimet

  if(doDebug)  cout<<"add MVA MET"<< endl;
  if(crecmvamet)
    {
      for(std::vector<edm::InputTag>::iterator mit = MvaMetCollectionsTag_.begin();
	    mit != MvaMetCollectionsTag_.end(); mit++){

	//collect MVA Mets
	edm::Handle<pat::METCollection> imets;
	iEvent.getByLabel(*mit, imets);

	if(!imets.isValid()) continue;
	if(imets->size() == 0)continue;

	for(std::vector<pat::MET>::const_iterator met = imets->begin(); met != imets->end(); met++){
	  if (mvamet_count==M_mvametmaxcount) {
	    cerr << "number of mvamet > M_mvametmaxcount. They are missing." << endl; errors |= 1<<1;
	    break;
	  }

	  mvamet_ex[mvamet_count] = met->px();
	  mvamet_ey[mvamet_count] = met->py();

	  mvamet_sigxx[mvamet_count] = met->getSignificanceMatrix()(0,0);
	  mvamet_sigxy[mvamet_count] = met->getSignificanceMatrix()(0,1);
	  mvamet_sigyx[mvamet_count] = met->getSignificanceMatrix()(1,0);
	  mvamet_sigyy[mvamet_count] = met->getSignificanceMatrix()(1,1);

	  if(met->userCand("lepton0")->isMuon() && (met->userCand("lepton1")->isMuon())){
	    mvamet_channel[mvamet_count] = MUMU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton1")->p4() );
	  }
	  else if (met->userCand("lepton0")->isMuon() && (met->userCand("lepton1")->isElectron())){
	    mvamet_channel[mvamet_count] = EMU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton1")->p4() );
	  }
	  else if (met->userCand("lepton0")->isMuon() && !(met->userCand("lepton1")->isPhoton())){
	    mvamet_channel[mvamet_count] = MUTAU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton0")->p4() );
	  }
	  else if (met->userCand("lepton0")->isElectron() && (met->userCand("lepton1")->isMuon())){
	    mvamet_channel[mvamet_count] = EMU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton0")->p4() );
	  }
	  else if (met->userCand("lepton0")->isElectron() && (met->userCand("lepton1")->isElectron())){
	    mvamet_channel[mvamet_count] = EE;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton1")->p4() );
	  }
	  else if (met->userCand("lepton0")->isElectron() && !(met->userCand("lepton1")->isPhoton())){
	    mvamet_channel[mvamet_count] = ETAU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton0")->p4() );
	  }
	  else if (!met->userCand("lepton0")->isPhoton() && (met->userCand("lepton1")->isMuon())){
	    mvamet_channel[mvamet_count] = MUTAU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton1")->p4() );
	  }
	  else if (!met->userCand("lepton0")->isPhoton() && (met->userCand("lepton1")->isElectron())){
	    mvamet_channel[mvamet_count] = ETAU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton1")->p4() );
	  }
	  else if (!met->userCand("lepton0")->isPhoton() && !(met->userCand("lepton1")->isPhoton())){
	    mvamet_channel[mvamet_count] = TAUTAU;
	    mvamet_lep1_pt[mvamet_count] = met->userCand("lepton0")->pt();
	    mvamet_lep2_pt[mvamet_count] = met->userCand("lepton1")->pt();
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton0")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	  }
	  else {
	    mvamet_channel[mvamet_count] = UNKNOWN;
	    mvamet_lep1[mvamet_count] = -1;
	    mvamet_lep2[mvamet_count] = -1;
	  }

	  mvamet_count++;
	}
      }
    }// crecmvamet

  if(doDebug)  cout<<"add rho"<< endl;
  // rho neutral
  edm::Handle<double> rho;
  iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAll"), rho);
  assert(rho.isValid());
  rhoNeutral = *rho;

  //  std::cout << "rhoNeutral = " << rhoNeutral << std::endl;

  genweight = 1.;
  numpileupinteractionsminus = -1;
  numpileupinteractions      = -1;
  numpileupinteractionsplus  = -1;
  numtruepileupinteractions  = -1.0f;
  hepNUP_ = -1;

  genparticles_lheHt = 0.;
  genparticles_noutgoing = 0;
  genparticles_noutgoing_NLO = 0;

  // generator info and generated particles
  if(doDebug)  cout<<"add gen info"<< endl;
  if(cembedded) {
    bool haveGenParticles = AddGenParticles(iEvent);
    edm::Handle<GenEventInfoProduct> GenEventInfo;
    iEvent.getByLabel(edm::InputTag("generator"), GenEventInfo);
    genweight = GenEventInfo->weight();
  }
  if(cgen && !cdata)
    {
      AddLHEInformation(iEvent);

      bool haveGenParticles = AddGenParticles(iEvent);
      bool haveGenJets      = AddGenJets(iEvent);

      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByLabel(edm::InputTag("generator"), HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  //	  cout << "Event weight from HEPMC : " << genweight << endl;
	  genid1 = HEPMC->pdf()->id.first;
	  genx1 = HEPMC->pdf()->x.first;
	  genid2 = HEPMC->pdf()->id.first;
	  genx2 = HEPMC->pdf()->x.second;
	  genScale = HEPMC->qScale();
	}

      edm::Handle<vector<PileupSummaryInfo> > PUInfo;
      iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PUInfo);
      if(PUInfo.isValid())
	{
	  for(vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	    {
	      int BX = PVI->getBunchCrossing();
	      if(BX == -1)
		{
		  numpileupinteractionsminus = PVI->getPU_NumInteractions();
		}
	      else if(BX == 0)
		{
		  numpileupinteractions = PVI->getPU_NumInteractions();
		}
	      else if(BX == 1)
		{
		  numpileupinteractionsplus = PVI->getPU_NumInteractions();
		}

	      numtruepileupinteractions = PVI->getTrueNumInteractions();
	    }
	}
    } // cgen

  if(crecTauSpinner&&cgen&&!cdata){

    for(int i=0;i<M_tauspinneranglesmaxcount;i++)TauSpinnerWeight[i]=1.;
    if(doDebug)  cout<<"add TauSpinner weights"<< endl;
    //ICTauSpinnerProducer * TauSpinner = new ICTauSpinnerProducer(iSetup.get<>);
    //TauSpinner->produce(iEvent,iSetup);
    TauSpinAngles_count= GetTauSpinnerweights(iEvent,iSetup);
  }


  if (ctrigger)
    {
      if(doDebug)  cout<<"add trigger info"<< endl;
      int numberOfTriggerObjects = int(AddTriggerObjects(iEvent,TriggerObjectCollectionToken_,*HLTrigger));
      //      std::cout << std::endl;
    } // ctrigger

  prefiringweight = 1.;
  prefiringweightup = 1.;
  prefiringweightdown = 1.;

  if( cYear != 2018) {
    edm::Handle< double > theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight ) ;
    prefiringweight =(*theprefweight);

    edm::Handle< double > theprefweightup;
    iEvent.getByToken(prefweightup_token, theprefweightup ) ;
    prefiringweightup =(*theprefweightup);

    edm::Handle< double > theprefweightdown;
    iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
    prefiringweightdown =(*theprefweightdown);
  }

  if(crecstxs)
    {
      // Get STXS infos
      edm::Handle<HTXS::HiggsClassification> htxs;
      iEvent.getByToken(htxsToken_, htxs);
      htxs_stage0cat = htxs->stage0_cat;
      htxs_stage1p1cat = htxs->stage1_1_cat_pTjet30GeV;
      htxs_higgsPt = htxs->higgs.Pt();
      htxs_njets30 = htxs->jets30.size();
    }


  tree->Fill();

  //Store Tau embedding information
  // edm::Handle<GenFilterInfo> embeddingWeightHandle;
  // iEvent.getByLabel(edm::InputTag("generator","minVisPtFilter",""), embeddingWeightHandle);
  // embeddingWeight_ = embeddingWeightHandle.isValid() ? embeddingWeightHandle->filterEfficiency() : 1.0;
  //cout << "--- EMBEDDING WEIGHT : " << embeddingWeight_ << endl;
  /*
    if (isRhEmb_){
    iEvent.getByLabel("genZdecayToTausForEmbeddingKineReweight", genDiTauHandle);
    iEvent.getByLabel(edm::InputTag("TauSpinnerReco","TauSpinnerWT"), TauSpinnerHandle);
    iEvent.getByLabel(edm::InputTag("ZmumuEvtSelEffCorrWeightProducer","weight"), ZmumuEffHandle);
    if(isMC_){
    iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genDiTauMassVsGenDiTauPt"), diTauMassVSdiTauPtHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genTau2EtaVsGenTau1Eta"), tau2EtaVStau1EtaHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genTau2PtVsGenTau1Pt"), tau2PtVStau1PtHandle);
    }
    else{
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genDiTauMassVsGenDiTauPt"), diTauMassVSdiTauPtHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genTau2EtaVsGenTau1Eta"), tau2EtaVStau1EtaHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genTau2PtVsGenTau1Pt"), tau2PtVStau1PtHandle);
    }
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weight"), muonRadiationHandle);
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weightDown"), muonRadiationDownHandle);
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weightUp"), muonRadiationUpHandle);

    TauSpinnerWeight = TauSpinnerHandle.isValid() ? (*TauSpinnerHandle) : 1.0;
    ZmumuEffWeight   = ZmumuEffHandle.isValid() ? (*ZmumuEffHandle) : 1.0;
    diTauMassVSdiTauPtWeight = diTauMassVSdiTauPtHandle.isValid() ? (*diTauMassVSdiTauPtHandle) : 1.0;
    tau2EtaVStau1EtaWeight = tau2EtaVStau1EtaHandle.isValid() ? (*tau2EtaVStau1EtaHandle) : 1.0;
    tau2PtVStau1PtWeight = tau2PtVStau1PtHandle.isValid() ? (*tau2PtVStau1PtHandle) : 1.0;
    muonRadiationWeight = muonRadiationHandle.isValid() ? (*muonRadiationHandle) : 1.0;
    muonRadiationDownWeight = muonRadiationDownHandle.isValid() ? (*muonRadiationDownHandle) : 1.0;
    muonRadiationUpWeight = muonRadiationUpHandle.isValid() ? (*muonRadiationUpHandle) : 1.0;
    genDiTauMass_ = genDiTauHandle.isValid() && genDiTauHandle->size()>0 ? genDiTauHandle->at(0).mass() : 9999;
  }
  embeddingWeights_->push_back(TauSpinnerWeight);
  embeddingWeights_->push_back(ZmumuEffWeight);
  embeddingWeights_->push_back(diTauMassVSdiTauPtWeight);
  embeddingWeights_->push_back(tau2EtaVStau1EtaWeight);
  embeddingWeights_->push_back(tau2PtVStau1PtWeight);
  embeddingWeights_->push_back(muonRadiationWeight);
  embeddingWeights_->push_back(muonRadiationDownWeight);
  embeddingWeights_->push_back(muonRadiationUpWeight);
  */



} //void NTupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


Int_t NTupleMaker::HasAnyMother(const GenParticle* particle, int id)
{
  vector<unsigned> bknummother;
  vector<const GenParticle*> bkparticle;
  bknummother.reserve(10);
  bkparticle.reserve(10);
  int level = 0;
  bkparticle.push_back(particle);
  bknummother.push_back(0);

  unsigned j = 0;
  while(true)
    {
      if(j == bkparticle[level]->numberOfMothers())
	{
	  level--;
	  if(level == -1){return(0);}
	  j = bknummother[level];
	  bkparticle.resize(level+1);
	  bknummother.resize(level+1);
	  continue;
	}

      if(bkparticle[level]->mother(j)->pdgId() == id) return(2);
      if(abs(bkparticle[level]->mother(j)->pdgId()) == abs(id)) return(1);

      if(bkparticle[level]->mother(j)->numberOfMothers() > 0)
	{
	  bknummother[level] = j+1;
	  bkparticle.push_back(dynamic_cast<const GenParticle*>(bkparticle[level]->mother(j)));
	  bknummother.push_back(0);
	  j = 0;
	  level++;
	  continue;
	}
      j++;
    }
  return(0);
} // Int_t NTupleMaker::HasAnyMother(const GenParticle* particle, int id)

void NTupleMaker::endJob()
{
}

math::XYZPoint NTupleMaker::PositionOnECalSurface(TransientTrack& trTrack)
{
	math::XYZPoint ecalPosition(0.,0.,0.);
	const FreeTrajectoryState myTSOS = trTrack.initialFreeState();
	TrajectoryStateOnSurface stateAtECal = propagatorWithMaterial->propagate(myTSOS, *ecalBarrel);

	if(stateAtECal.isValid() && stateAtECal.globalPosition().eta() > 1.479)
	{
		stateAtECal= propagatorWithMaterial->propagate(myTSOS, *ecalPositiveEtaEndcap);
	}

	if(stateAtECal.isValid() && stateAtECal.globalPosition().eta() < -1.479)
	{
		stateAtECal= propagatorWithMaterial->propagate(myTSOS, *ecalNegativeEtaEndcap);
	}

	if(stateAtECal.isValid())
	{
		ecalPosition = stateAtECal.globalPosition();
	}
	return(ecalPosition);
}

bool NTupleMaker::AddSusyInfo(const edm::Event& iEvent) {

  SusyMotherMass = -1;
  SusyLSPMass = -1;

  bool success = true;


    edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByToken(GenParticleCollectionToken_, GenParticles);


    if(GenParticles.isValid())
    {
      //bool count_partons = false;
      //passed = true;
      for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	{

    if(abs((*GenParticles)[i].pdgId()) == 1000015 ) SusyMotherMass =  (*GenParticles)[i].mass();
    if(abs((*GenParticles)[i].pdgId()) == 1000022 ) SusyLSPMass =  (*GenParticles)[i].mass();
	}
	}
    return success;

  /*
  edm::Handle<double> susyMotherMass;
  iEvent.getByToken(SusyMotherMassToken_,susyMotherMass);
  if (susyMotherMass.isValid()) {
    SusyMotherMass = *susyMotherMass;
  }
  else {
    success = false;
  }
  edm::Handle<double> susyLSPMass;
  iEvent.getByToken(SusyLSPMassToken_,susyLSPMass);
  if (susyLSPMass.isValid()) {
    SusyLSPMass = *susyLSPMass;
  }
  else {
    success = false;
  }

  return success;*/

}

bool NTupleMaker::AddFlags(const edm::Event& iEvent, const char* module, const char* label, const char* process) {
  iEvent.getByLabel(edm::InputTag( module, label, process), Flags);
  if (!Flags.isValid())
    return false;

  const edm::TriggerNames& FlagNames_ = iEvent.triggerNames(*Flags);
  for(unsigned i = 0 ; i < Flags->size(); i++){
    if(!Flags->wasrun(i) )continue;
    std::string flagName=FlagNames_.triggerName(i);
    if(cFlags.size() > 0){
      for(size_t ip = 0; ip < cFlags.size(); ip++){
	if(flagName.find(cFlags[ip]) != string::npos){

	  flags_->insert(std::pair<string, int>(flagName, Flags->accept(i)));
	  TString TriggerName(flagName);
	  //std::cout << flagName << " : " << Flags->accept(i) << std::endl;
	}
      }
    }
  }

  return true;
}

bool NTupleMaker::AddLHEInformation(const edm::Event& iEvent) {
  genparticles_lheHt = 0.;
  genparticles_noutgoing = 0;
  genparticles_noutgoing_NLO = 0;
  genparticles_lheWPt = 0.;

  edm::Handle<LHEEventProduct> lheEventProduct;
  iEvent.getByToken(LHEToken_, lheEventProduct);

  if(!lheEventProduct.isValid())
    return false;

  const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

  size_t numParticles = lheParticles.size();

  LorentzVector lepLV(0,0,0,0);
  LorentzVector nuLV(0,0,0,0);

  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
   int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
   int status = lheEvent.ISTUP[idxParticle];
   if ( absPdgId == 24 ){
     genparticles_lheWPt = TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
   }
   if ( absPdgId == 11 || absPdgId == 13 || absPdgId == 15 ){
     if( status == 1) lepLV.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); // see:  SimDataFormats/GeneratorProducts/interface/LesHouches.h
   }
   if ( absPdgId == 12 || absPdgId == 14 || absPdgId == 16 ){
     if( status == 1) nuLV.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); // see:  SimDataFormats/GeneratorProducts/interface/LesHouches.h
   }
   if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
       genparticles_lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
       ++genparticles_noutgoing;
   }
  }

  // Calculate w-boson pt from decay products
  LorentzVector wLV = lepLV + nuLV;
  if(genparticles_lheWPt == 0)  genparticles_lheWPt = wLV.Pt();

  genparticles_noutgoing_NLO = lheEventProduct->npNLO();

  weightPDFmax = 0.1;
  weightPDFmin = 10;

  float nPDFWeights = 0;
  float nPDFWeightsUp = 0;
  float nPDFWeightsDown = 0;
  float weightPDF2 = 0;

  weightPDFmean = 0;
  weightPDFup = 0;
  weightPDFdown = 0;

  for (unsigned int iW=0; iW<lheEventProduct->weights().size(); iW++) {
    //    std::cout << lheEventProduct->weights()[iW].id << " : " << lheEventProduct->weights()[iW].wgt/lheEventProduct->originalXWGTUP() << std::endl;

    float weightGen = lheEventProduct->weights()[iW].wgt/lheEventProduct->originalXWGTUP();
    if (iW<9)
      weightScale[iW] = weightGen;
    else {
      weightPDFmean += weightGen;
      weightPDF2 += weightGen*weightGen;
      nPDFWeights += 1.0;
      if (weightGen>weightPDFmax) weightPDFmax = weightGen;
      if (weightGen<weightPDFmin) weightPDFmin = weightGen;
      if (weightGen>1.0) {
	weightPDFup += weightGen;
	nPDFWeightsUp += 1.0;
      }
      else {
	weightPDFdown += weightGen;
	nPDFWeightsDown += 1.0;
      }
    }
  }
  //  std::cout << std::endl;
  weightPDFmean = weightPDFmean/nPDFWeights;
  weightPDF2 = weightPDF2/nPDFWeights;
  weightPDFvar = TMath::Sqrt(weightPDF2 - weightPDFmean*weightPDFmean);
  weightPDFup = weightPDFup/nPDFWeightsUp;
  weightPDFdown = weightPDFdown/nPDFWeightsDown;

  return true;

}

bool NTupleMaker::AddGenParticles(const edm::Event& iEvent) {

  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByToken(GenParticleCollectionToken_, GenParticles);

  bool passed = false;

  if(GenParticles.isValid())
    {
      bool count_partons = false;
      passed = true;
      for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	{
	  bool fill = false;
	  UInt_t info = 0;
	  UInt_t mother = 100;
	  const GenStatusFlags statusFlags = (*GenParticles)[i].statusFlags();

	  if(abs((*GenParticles)[i].pdgId()) == 13)// /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
		 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenMuon : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta()
	      //			<< "   phi = " << (*GenParticles)[i].phi()
	      //			<< "   status = " << (*GenParticles)[i].status()
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 11)// /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
                 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenElectron : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta()
	      //			<< "   phi = " << (*GenParticles)[i].phi()
	      //			<< "   status = " << (*GenParticles)[i].status()
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 15 /*&& (*GenParticles)[i].pt() > 10.*/)
	    {
	      fill = false;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
                 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenTau : "
	      //	       		<< "   pt = " << (*GenParticles)[i].pt()
	      //	       		<< "   eta = " << (*GenParticles)[i].eta()
	      //	       		<< "   phi = " << (*GenParticles)[i].phi()
	      //	       		<< "   status = " << (*GenParticles)[i].status()
	      //	       		<< "   mother = " << mother << std::endl;
	      reco::Candidate::LorentzVector tau_visible_p4 = getVisMomentum(&(*GenParticles)[i]);
	      reco::Candidate::LorentzVector tau_visibleNoLep_p4 = utils_genMatch::getVisMomentumNoLep(&(*GenParticles)[i]);
	      // std::cout << "   visible pt = " << tau_visible_p4.pt()
	      // 		<< "   eta = " << tau_visible_p4.eta()
	      // 		<< "   phi = " << tau_visible_p4.phi()
	      // 		<< "   mode = " << getGenTauDecayMode(&(*GenParticles)[i]) << std::endl;

	      std::string genTauDecayMode = getGenTauDecayMode(&(*GenParticles)[i]);

	      gentau_px[gentau_count] = (*GenParticles)[i].px();
	      gentau_py[gentau_count] = (*GenParticles)[i].py();
	      gentau_pz[gentau_count] = (*GenParticles)[i].pz();
	      gentau_e[gentau_count]  = (*GenParticles)[i].energy();
	      gentau_charge[gentau_count]  = (*GenParticles)[i].charge();
	      gentau_status[gentau_count] = (*GenParticles)[i].status();

	      gentau_visible_px[gentau_count] = tau_visible_p4.px();
	      gentau_visible_py[gentau_count] = tau_visible_p4.py();
	      gentau_visible_pz[gentau_count] = tau_visible_p4.pz();
	      gentau_visible_e[gentau_count]  = tau_visible_p4.energy();

	      gentau_visible_pt[gentau_count]   = tau_visible_p4.pt();
	      gentau_visible_eta[gentau_count]  = tau_visible_p4.eta();
	      gentau_visible_phi[gentau_count]  = tau_visible_p4.phi();
	      gentau_visible_mass[gentau_count] = tau_visible_p4.mass();

	      gentau_visibleNoLep_e[gentau_count]  = tau_visibleNoLep_p4.energy();
	      gentau_visibleNoLep_pt[gentau_count]   = tau_visibleNoLep_p4.pt();
	      gentau_visibleNoLep_eta[gentau_count]  = tau_visibleNoLep_p4.eta();
	      gentau_visibleNoLep_phi[gentau_count]  = tau_visibleNoLep_p4.phi();
	      gentau_visibleNoLep_mass[gentau_count] = tau_visibleNoLep_p4.mass();

	      gentau_fromHardProcess[gentau_count] = statusFlags.fromHardProcess();
	      gentau_fromHardProcessBeforeFSR[gentau_count] = statusFlags.fromHardProcessBeforeFSR();
	      gentau_isDecayedLeptonHadron[gentau_count] = statusFlags.isDecayedLeptonHadron();
	      gentau_isDirectHadronDecayProduct[gentau_count] = statusFlags.isDirectHadronDecayProduct();
	      gentau_isDirectHardProcessTauDecayProduct[gentau_count] = statusFlags.isDirectHardProcessTauDecayProduct();
	      gentau_isDirectPromptTauDecayProduct[gentau_count] = statusFlags.isDirectPromptTauDecayProduct();
	      gentau_isDirectTauDecayProduct[gentau_count] = statusFlags.isDirectTauDecayProduct();
	      gentau_isFirstCopy[gentau_count] = statusFlags.isFirstCopy();
	      gentau_isHardProcess[gentau_count] = statusFlags.isHardProcess();
	      gentau_isHardProcessTauDecayProduct[gentau_count] = statusFlags.isHardProcessTauDecayProduct();
	      gentau_isLastCopy[gentau_count] = statusFlags.isLastCopy();
	      gentau_isLastCopyBeforeFSR[gentau_count] = statusFlags.isLastCopyBeforeFSR();
	      gentau_isPrompt[gentau_count] = statusFlags.isPrompt();
	      gentau_isPromptTauDecayProduct[gentau_count] = statusFlags.isPromptTauDecayProduct();
	      gentau_isTauDecayProduct[gentau_count] = statusFlags.isTauDecayProduct();

	      gentau_decayMode_name[gentau_count] = genTauDecayMode;

	      if( genTauDecayMode.find("oneProng0Pi0")!=string::npos )          gentau_decayMode[gentau_count] = 0;
	      else if( genTauDecayMode.find("oneProng1Pi0")!=string::npos )     gentau_decayMode[gentau_count] = 1;
	      else if( genTauDecayMode.find("oneProng2Pi0")!=string::npos )     gentau_decayMode[gentau_count] = 2;
	      else if( genTauDecayMode.find("oneProngOther")!=string::npos )    gentau_decayMode[gentau_count] = 3;
	      else if( genTauDecayMode.find("threeProng0Pi0")!=string::npos )   gentau_decayMode[gentau_count] = 4;
	      else if( genTauDecayMode.find("threeProng1Pi0")!=string::npos )   gentau_decayMode[gentau_count] = 5;
	      else if( genTauDecayMode.find("threeProngOther")!=string::npos )  gentau_decayMode[gentau_count] = 6;
	      else if( genTauDecayMode.find("rare")!=string::npos )             gentau_decayMode[gentau_count] = 7;
	      else if( genTauDecayMode.find("muon")!=string::npos )             gentau_decayMode[gentau_count] = 8;
	      else if( genTauDecayMode.find("electron")!=string::npos )         gentau_decayMode[gentau_count] = 9;
	      else   tau_genDecayMode[gentau_count] = -99;

	      gentau_mother[gentau_count] = mother;

	      gentau_count++;

	    }
	  else if( (abs((*GenParticles)[i].pdgId()) == 16 || abs((*GenParticles)[i].pdgId()) == 14 || abs((*GenParticles)[i].pdgId()) == 12))
	    {
	      if ((*GenParticles)[i].status()==1) {
		fill = true;
		if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;mother=ZBOSON;}
		if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;mother=WBOSON;}
		if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2;mother=TAU;}
		//		std::cout << "GenNeutrino : " << (*GenParticles)[i].pdgId()
		//		 	  << "   pt = " << (*GenParticles)[i].pt()
		//		 	  << "   eta = " << (*GenParticles)[i].eta()
		//		 	  << "   phi = " << (*GenParticles)[i].phi()
		//			  << "   status = " << (*GenParticles)[i].status()
		//		 	  << "   mother =  " << mother << std::endl;

	      }
	    }
	  // Save partons (quarks)
	  else if(abs((*GenParticles)[i].pdgId()) < 6)
	    {
	      if ((*GenParticles)[i].status()==3 && count_partons)
		fill = true;
	    }
	  // Save all tops from Madgraph
	  else if ( abs((*GenParticles)[i].pdgId()) == 6)// && (*GenParticles)[i].status()==62 )
	    {
	      fill = true;
	      //	      	std::cout << "GenTop : " << (*GenParticles)[i].pdgId()
	      //			  << "   pt = " << (*GenParticles)[i].pt()
	      //		 	  << "   eta = " << (*GenParticles)[i].eta()
	      //		 	  << "   phi = " << (*GenParticles)[i].phi()
	      //		 	  << "   status = " << (*GenParticles)[i].status() << std::endl;
	    }
	  // Save partons (gluons)
  	  else if(abs((*GenParticles)[i].pdgId()) == 21)
	    {
	      if ((*GenParticles)[i].status()==3 && count_partons)
		fill = true;
	    }
	  // Save photons
	  else if(abs((*GenParticles)[i].pdgId()) == 22)
	    {
	      if (statusFlags.isPrompt())
		fill = true;
	    }

	  // Save all W/Z bosons from Madgraph
	  else if(abs((*GenParticles)[i].pdgId()) == 23 || abs((*GenParticles)[i].pdgId()) == 24 )
	    {
	      count_partons = true;
	      fill = true;
	      int nDaughters = (*GenParticles)[i].numberOfDaughters();
	      bool posElectronFound = false;
	      bool negElectronFound = false;
	      bool posMuonFound = false;
	      bool negMuonFound = false;
	      bool posTauFound = false;
	      bool negTauFound = false;
	      for (int iD=0 ; iD<nDaughters; ++iD) {
		const reco::Candidate * kid = (*GenParticles)[i].daughter(iD);
		int pdgId = kid->pdgId();
		if (pdgId==11) negElectronFound = true;
		if (pdgId==-11) posElectronFound = true;
		if (pdgId==13) negMuonFound = true;
		if (pdgId==-13) posMuonFound = true;
		if (pdgId==15) negTauFound = true;
		if (pdgId==-15) posTauFound = true;
	      }
	      if (posElectronFound&&negElectronFound) info = 11;
	      else if (posMuonFound&&negMuonFound) info = 13;
	      else if (posTauFound&&negTauFound) info = 15;
	      //	      std::cout << "GenBoson : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta()
	      //			<< "   phi = " << (*GenParticles)[i].phi()
	      //			<< "   status = " << (*GenParticles)[i].status()
	      //			<< "   info = " << info << std::endl;
	    }
	  //Save Higgs bosons
	  else if(abs((*GenParticles)[i].pdgId()) == 25 || abs((*GenParticles)[i].pdgId()) == 35 ||
		  abs((*GenParticles)[i].pdgId()) == 36){
	    fill = true;
	  }

	  else if( ((*GenParticles)[i].status()==1 ||
		    statusFlags.fromHardProcess() ||
		    statusFlags.isDirectHardProcessTauDecayProduct()) &&
		   abs((*GenParticles)[i].pdgId())!=15) {
	    if(HasAnyMother(&(*GenParticles)[i], 23) > 0) {info |= 1<<0;mother=ZBOSON;}
	    if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;mother=WBOSON;}
	    if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2;mother=TAU;}
	    if(HasAnyMother(&(*GenParticles)[i], 25) > 0||
	       HasAnyMother(&(*GenParticles)[i], 35) > 0||
	       HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3;mother=HIGGS;}
	    fill = true;
	  }

	  if(fill)
	    {
	      genparticles_e[genparticles_count] = (*GenParticles)[i].energy();
	      genparticles_px[genparticles_count] = (*GenParticles)[i].px();
	      genparticles_py[genparticles_count] = (*GenParticles)[i].py();
	      genparticles_pz[genparticles_count] = (*GenParticles)[i].pz();
	      genparticles_vx[genparticles_count] = (*GenParticles)[i].vx();
	      genparticles_vy[genparticles_count] = (*GenParticles)[i].vy();
	      genparticles_vz[genparticles_count] = (*GenParticles)[i].vz();
	      genparticles_pdgid[genparticles_count] = (*GenParticles)[i].pdgId();
	      genparticles_status[genparticles_count] = (*GenParticles)[i].status();
	      genparticles_info[genparticles_count] = info;
	      genparticles_mother[genparticles_count] = mother;

	      genparticles_fromHardProcess[genparticles_count] = statusFlags.fromHardProcess();
	      genparticles_fromHardProcessBeforeFSR[genparticles_count] = statusFlags.fromHardProcessBeforeFSR();
	      genparticles_isDecayedLeptonHadron[genparticles_count] = statusFlags.isDecayedLeptonHadron();
	      genparticles_isDirectHadronDecayProduct[genparticles_count] = statusFlags.isDirectHadronDecayProduct();
	      genparticles_isDirectHardProcessTauDecayProduct[genparticles_count] = statusFlags.isDirectHardProcessTauDecayProduct();
	      genparticles_isDirectPromptTauDecayProduct[genparticles_count] = statusFlags.isDirectPromptTauDecayProduct();
	      genparticles_isDirectTauDecayProduct[genparticles_count] = statusFlags.isDirectTauDecayProduct();
	      genparticles_isFirstCopy[genparticles_count] = statusFlags.isFirstCopy();
	      genparticles_isHardProcess[genparticles_count] = statusFlags.isHardProcess();
	      genparticles_isHardProcessTauDecayProduct[genparticles_count] = statusFlags.isHardProcessTauDecayProduct();
	      genparticles_isLastCopy[genparticles_count] = statusFlags.isLastCopy();
	      genparticles_isLastCopyBeforeFSR[genparticles_count] = statusFlags.isLastCopyBeforeFSR();
	      genparticles_isPrompt[genparticles_count] = statusFlags.isPrompt();
	      genparticles_isPromptTauDecayProduct[genparticles_count] = statusFlags.isPromptTauDecayProduct();
	      genparticles_isTauDecayProduct[genparticles_count] = statusFlags.isTauDecayProduct();

	      genparticles_count++;
	    }
	} // for(unsigned i = 0 ; i < GenParticles->size() ; i++)
    } // if(GenParticles.isValid())

  return passed;

} // bool NTupleMaker::AddGenParticles(const edm::Event& iEvent)


bool NTupleMaker::AddGenJets(const edm::Event& iEvent)
{
  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByToken(GenJetCollectionToken_, genjets);

  bool passed = false;

  if(genjets.isValid())
    {
      bool passed = true;

      for(unsigned i = 0 ; i < genjets->size() ; i++)
	{
	  if(genjets_count == M_genjetsmaxcount){
	    cerr << "number of genjets_count > M_genjetsmaxcount. They are missing." << endl;
	    errors |= 1<<4;
	    break;
	  }

	  if(fabs((*genjets)[i].eta()) > 5.2) continue;

	  genjets_e[genjets_count]  = (*genjets)[i].energy();
	  genjets_px[genjets_count] = (*genjets)[i].px();
	  genjets_py[genjets_count] = (*genjets)[i].py();
	  genjets_pz[genjets_count] = (*genjets)[i].pz();
	  genjets_pt[genjets_count] = (*genjets)[i].pt();
	  genjets_eta[genjets_count] = (*genjets)[i].eta();
	  genjets_phi[genjets_count] = (*genjets)[i].phi();
	  genjets_pdgid[genjets_count]  = (*genjets)[i].pdgId();
	  genjets_status[genjets_count] = (*genjets)[i].status();

	  genjets_em_energy[genjets_count]        = (*genjets)[i].emEnergy();
	  genjets_had_energy[genjets_count]       = (*genjets)[i].hadEnergy();
	  genjets_invisible_energy[genjets_count] = (*genjets)[i].invisibleEnergy();
	  genjets_auxiliary_energy[genjets_count] = (*genjets)[i].auxiliaryEnergy();
	  genjets_count++;
	}
    }
  return  passed;
}


unsigned int NTupleMaker::AddPFCand(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<pat::PackedCandidateCollection> Tracks;
  iEvent.getByToken( PackedCantidateCollectionToken_, Tracks);

  if(Tracks.isValid())
    {
      for(unsigned i = 0 ; i < Tracks->size() ; i++){
        if ((*Tracks)[i].pt() < cTrackPtMin) continue;
        if (fabs((*Tracks)[i].eta()) > cTrackEtaMax) continue;
        if (fabs((*Tracks)[i].charge()) < 0.5) continue;
	if (fabs((*Tracks)[i].dxy()) > cTrackDxyMax) continue;
	if (fabs((*Tracks)[i].dz()) > cTrackDzMax) continue;
        track_px[track_count] = (*Tracks)[i].px();
        track_py[track_count] = (*Tracks)[i].py();
        track_pz[track_count] = (*Tracks)[i].pz();
        track_pt[track_count] = (*Tracks)[i].pt();
        track_eta[track_count] = (*Tracks)[i].eta();
        track_phi[track_count] = (*Tracks)[i].phi();
        track_charge[track_count] = (*Tracks)[i].charge();
        track_mass[track_count] = (*Tracks)[i].mass();
        track_dxy[track_count] = (*Tracks)[i].dxy();
        track_dz[track_count] = (*Tracks)[i].dz();
	if((*Tracks)[i].hasTrackDetails()){
	  track_dxyerr[track_count] = (*Tracks)[i].dxyError();
	  track_dzerr[track_count] = (*Tracks)[i].dzError();
	}
	else{
	  track_dxyerr[track_count] = -9999;
	  track_dzerr[track_count] = -9999;
	}
	track_vx[track_count] = (*Tracks)[i].vertex().x();
	track_vy[track_count] = (*Tracks)[i].vertex().y();
	track_vz[track_count] = (*Tracks)[i].vertex().z();
        track_ID[track_count] = (*Tracks)[i].pdgId();
	const reco::Track * trkRef = (*Tracks)[i].bestTrack();
	track_highPurity[track_count] = false;
	if (trkRef != NULL) {
	  track_highPurity[track_count] = trkRef->quality(reco::Track::highPurity);
	}
        track_count++;

        if (track_count==M_trackmaxcount) {
          cerr << "number of tracks > M_trackmaxcount. They are missing." << endl; errors |= 1<<1;
          break;
        }
      }

    }

  return track_count;

}

unsigned int NTupleMaker::AddMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);
  const TransientTrackBuilder * transientTrackBuilder = builder.product();

  edm::Handle<pat::MuonCollection> Muons;
  iEvent.getByToken(MuonCollectionToken_, Muons);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken( PackedCantidateCollectionToken_, pfcands);
  /*
  edm::Handle<edm::PtrVector<reco::Muon>> BadDuplicateMuons;
  iEvent.getByToken(BadDuplicateMuonsToken_, BadDuplicateMuons);

  edm::Handle<edm::PtrVector<reco::Muon>> BadGlobalMuons;
  iEvent.getByToken(BadGlobalMuonsToken_, BadGlobalMuons);
  */


  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){

	if (muon_count==M_muonmaxcount) {
	  cerr << "number of muons > M_muonmaxcount. They are missing." << endl; errors |= 1<<1;
	  break;
	}

	if ((*Muons)[i].pt() < cMuPtMin) continue;
	if (fabs(((*Muons)[i].eta()))>cMuEtaMax) continue;


//code below is to store the track param vec + covariances, a references pt on track, and B field in ref pt. For CP measurement

	reco::TrackRef leadTrk = (*Muons)[i].innerTrack();
	if (leadTrk.isNonnull()) {
	TrackBase::ParameterVector ParamVecMu=leadTrk->parameters();
	TrackBase::CovarianceMatrix CVMTrack=leadTrk->covariance();

	for(int index=0; index<ParamVecMu.kSize;index++){
	  muon_helixparameters[muon_count][index]=ParamVecMu[index];
	  for(int index2=0; index2<ParamVecMu.kSize;index2++){
	    muon_helixparameters_covar[muon_count][index][index2]=CVMTrack[index][index2];
	  }
	}

	TrackBase::Point RFptMu=leadTrk->referencePoint();

	muon_referencePoint[muon_count][0]=RFptMu.X();
	muon_referencePoint[muon_count][1]=RFptMu.Y();
	muon_referencePoint[muon_count][2]=RFptMu.Z();
	
	double	magneticField = (TTrackBuilder.product() ? TTrackBuilder.product()->field()->inInverseGeV(GlobalPoint(RFptMu.X(), RFptMu.Y(), RFptMu.Z())).z() : 0.0);
	muon_Bfield[muon_count]=magneticField;
	}
	else{
	for(int index=0; index<5;index++){
	  muon_helixparameters[muon_count][index]=-999;
	  for(int index2=0; index2<5;index2++){
	    muon_helixparameters_covar[muon_count][index][index2]=-999;
	  }
	}

	muon_referencePoint[muon_count][0]=-999;
	muon_referencePoint[muon_count][1]=-999;
	muon_referencePoint[muon_count][2]=-999;	  
	muon_Bfield[muon_count]=-999;
}

	muon_px[muon_count] = (*Muons)[i].px();
	muon_py[muon_count] = (*Muons)[i].py();
	muon_pz[muon_count] = (*Muons)[i].pz();
	muon_pt[muon_count] = (*Muons)[i].pt();
	muon_eta[muon_count] = (*Muons)[i].eta();
	muon_phi[muon_count] = (*Muons)[i].phi();
	muon_charge[muon_count] = (*Muons)[i].charge();
	muon_vx[muon_count] = (*Muons)[i].vx(); // gives the same as (*Muons)[i].muonBestTrack()->vx()
	muon_vy[muon_count] = (*Muons)[i].vy();
	muon_vz[muon_count] = (*Muons)[i].vz();


	const pat::Muon &lep = (*Muons)[i];
	muon_miniISO[muon_count]=getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false);

	if((*Muons)[i].globalTrack().isNonnull())
	  {
	    muon_globalTrack[muon_count] = true;
	    muon_pterror[muon_count] = (*Muons)[i].globalTrack()->ptError();
	    muon_chi2[muon_count] = (*Muons)[i].globalTrack()->chi2();
	    muon_ndof[muon_count] = (*Muons)[i].globalTrack()->ndof();
	    muon_nMuonHits[muon_count] = (*Muons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
	    muon_normChi2[muon_count] = (*Muons)[i].globalTrack()->normalizedChi2();
	  }
	else
	  {
	    muon_globalTrack[muon_count] = false;
	    muon_pterror[muon_count] = -1.;
	    muon_chi2[muon_count] = -1.;
	    muon_ndof[muon_count] = 0;
	    muon_nMuonHits[muon_count] = 0;
	    muon_normChi2[muon_count] = -1;
	  }

	//	std::cout << "  chi2 = " << muon_chi2[muon_count] << "  ndof = " << muon_ndof[muon_count] << std::endl;

	muon_nMuonStations[muon_count] = (*Muons)[i].numberOfMatchedStations();

	muon_isTracker[muon_count] = (*Muons)[i].isTrackerMuon();
	muon_isPF[muon_count] = (*Muons)[i].isPFMuon();
	muon_isTight[muon_count] = (*Muons)[i].isTightMuon(primvertex);
	muon_isLoose[muon_count] = (*Muons)[i].isLooseMuon();
	muon_isGlobal[muon_count] = (*Muons)[i].isGlobalMuon();
	muon_isMedium[muon_count] = (*Muons)[i].isMediumMuon();
	muon_isICHEP[muon_count] = idAlgos::isICHEPMuon((*Muons)[i]);

	muon_chargedHadIso[muon_count] = (*Muons)[i].chargedHadronIso();
	muon_neutralHadIso[muon_count] = (*Muons)[i].neutralHadronIso();
	muon_photonIso[muon_count] = (*Muons)[i].photonIso();
	muon_puIso[muon_count] = (*Muons)[i].puChargedHadronIso();

	muon_r03_sumChargedHadronPt[muon_count] = (*Muons)[i].pfIsolationR03().sumChargedHadronPt;
	muon_r03_sumChargedParticlePt[muon_count] = (*Muons)[i].pfIsolationR03().sumChargedParticlePt;
	muon_r03_sumNeutralHadronEt[muon_count] = (*Muons)[i].pfIsolationR03().sumNeutralHadronEt;
	muon_r03_sumPhotonEt[muon_count] = (*Muons)[i].pfIsolationR03().sumPhotonEt;
	muon_r03_sumNeutralHadronEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR03().sumNeutralHadronEtHighThreshold;
	muon_r03_sumPhotonEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR03().sumPhotonEtHighThreshold;
	muon_r03_sumPUPt[muon_count] = (*Muons)[i].pfIsolationR03().sumPUPt;

        muon_r04_sumChargedHadronPt[muon_count] = (*Muons)[i].pfIsolationR04().sumChargedHadronPt;
        muon_r04_sumChargedParticlePt[muon_count] = (*Muons)[i].pfIsolationR04().sumChargedParticlePt;
        muon_r04_sumNeutralHadronEt[muon_count] = (*Muons)[i].pfIsolationR04().sumNeutralHadronEt;
        muon_r04_sumPhotonEt[muon_count] = (*Muons)[i].pfIsolationR04().sumPhotonEt;
        muon_r04_sumNeutralHadronEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR04().sumNeutralHadronEtHighThreshold;
        muon_r04_sumPhotonEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR04().sumPhotonEtHighThreshold;
        muon_r04_sumPUPt[muon_count] = (*Muons)[i].pfIsolationR04().sumPUPt;

	TrackRef innertrack = (*Muons)[i].innerTrack();
	TrackRef bestTrack  = (*Muons)[i].muonBestTrack();

	muon_combQ_trkKink[muon_count] = (*Muons)[i].combinedQuality().trkKink;
	muon_combQ_chi2LocalPosition[muon_count] = (*Muons)[i].combinedQuality().chi2LocalPosition;
	muon_segmentComp[muon_count] = (*Muons)[i].segmentCompatibility();

	if (bestTrack.isNonnull()) {
	  muon_dxy[muon_count]    = bestTrack->dxy(pv_position);
	  muon_dz[muon_count]     = bestTrack->dz(pv_position);
	  muon_dxyerr[muon_count]    = bestTrack->dxyError();
	  muon_dzerr[muon_count]     = bestTrack->dzError();
	}
	else {
	  muon_dxy[muon_count]    = -9999;
	  muon_dz[muon_count]     = -9999;
	  muon_dxyerr[muon_count]    = -9999;
	  muon_dzerr[muon_count]     = -9999;
	}

	if(innertrack.isNonnull())
	  {
	    muon_innerTrack[muon_count] = true;
	    muon_nPixelHits[muon_count] = innertrack->hitPattern().numberOfValidPixelHits();
	    muon_nTrackerHits[muon_count] = innertrack->hitPattern().trackerLayersWithMeasurement();
	    muon_validFraction[muon_count] = innertrack->validFraction();
	  }
	else
	  {
	    muon_innerTrack[muon_count] = false;
	    muon_nPixelHits[muon_count] = 0;
	    muon_nTrackerHits[muon_count] = 0;
	    muon_validFraction[muon_count] = 0;
	  }

	//	bool goodGlb = muon_isGlobal[muon_count] && muon_normChi2[muon_count]  < 3 && muon_normChi2[muon_count] > 0
	//	 && muon_combQ_chi2LocalPosition[muon_count] < 12 && muon_combQ_trkKink[muon_count] < 20;
	//	muon_isMedium[muon_count] =  muon_isLoose[muon_count] && muon_validFraction[muon_count] > 0.8 && muon_segmentComp[muon_count] > (goodGlb ? 0.303 : 0.451);

	muon_genmatch[muon_count] = 0;
	if(cgen && (!cdata || cembedded)){
	  edm::Handle<reco::GenParticleCollection> GenParticles;
	  iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	  if(GenParticles.isValid())
	    muon_genmatch[muon_count] = utils_genMatch::genMatch( (*Muons)[i].p4(), *GenParticles);
	}

	// Duplicate & bad muons
	muon_isDuplicate[muon_count] = false;
	//	for(unsigned j = 0; j < BadDuplicateMuons->size(); j++){
	//	  if((*BadDuplicateMuons)[j]->pt() == muon_pt[muon_count]) muon_isDuplicate[muon_count] = true;
	//	}
	muon_isBad[muon_count] = false;
	//	for(unsigned j = 0; j < BadGlobalMuons->size(); j++){
	//	  if((*BadGlobalMuons)[j]->pt() == muon_pt[muon_count]) muon_isBad[muon_count] = true;
	//	}

	// Dimuons
	/*if( !(*Muons)[i].innerTrack().isNull()){
	  for(unsigned j = i+1 ; j < Muons->size() ; j++){
	    if (dimuon_count==M_muonmaxcount*(M_muonmaxcount-1)/2) {
	      cerr << "number of dimuons > M_muonmaxcount*(M_muonmaxcount-1)/2. They are missing." << endl; errors |= 1<<1;
	      break;
	    }

	    if ((*Muons)[j].pt() < cMuPtMin) continue;
	    if (fabs(((*Muons)[j].eta()))>cMuEtaMax) continue;
	    if( (*Muons)[j].innerTrack().isNull()) continue;

	    dimuon_leading[dimuon_count] = i;
	    dimuon_trailing[dimuon_count] = j;
	    if ((*Muons)[i].pt() < (*Muons)[j].pt()){
	      dimuon_leading[dimuon_count] = j;
	      dimuon_trailing[dimuon_count] = i;
	    }

	    if (fabs( (*Muons)[i].pt() - (*Muons)[j].pt()) < 1.e-4){
	      std::cout<<"WTF!!"<<std::endl;
	    }

	    dimuon_dist2D[dimuon_count] = -1.;
	    dimuon_dist2DE[dimuon_count] = -1.;
	    dimuon_dist3D[dimuon_count] = -1.;
	    dimuon_dist3DE[dimuon_count] = -1.;

	    const pat::Muon* leading = &(*Muons)[dimuon_leading[dimuon_count]];
	    const pat::Muon* trailing = &(*Muons)[dimuon_trailing[dimuon_count]];

	    reco::TrackRef leadTrk = leading->innerTrack();
	    reco::TrackRef trailTrk = trailing->innerTrack();

	    TransientTrack trLeadTrk = transientTrackBuilder->build(*leadTrk);
	    TransientTrack trTrailTrk = transientTrackBuilder->build(*trailTrk);

	    FreeTrajectoryState leadState = trLeadTrk.impactPointTSCP().theState();
	    FreeTrajectoryState trailState = trTrailTrk.impactPointTSCP().theState();

	    if (trLeadTrk.impactPointTSCP().isValid() && trTrailTrk.impactPointTSCP().isValid()) {
	      TwoTrackMinimumDistance minDist;

	      typedef ROOT::Math::SVector<double, 3> SVector3;
	      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;

	      minDist.calculate(leadState,trailState);
	      if (minDist.status()) {

		//float dist3D = minDist.distance();
		std::pair<GlobalPoint,GlobalPoint> pcaMuons = minDist.points();
		//GlobalPoint posPCA = pcaMuons.first;
		//GlobalPoint negPCA = pcaMuons.second;

		ParticleMass muon_mass = 0.105658;
		float muon_sigma = muon_mass*1.e-6;

		//Creating a KinematicParticleFactory
		KinematicParticleFactoryFromTransientTrack pFactory;

		//initial chi2 and ndf before kinematic fits.
		float chi = 0.;
		float ndf = 0.;
		RefCountedKinematicParticle leadPart = pFactory.particle(trLeadTrk,muon_mass,chi,ndf,muon_sigma);
		RefCountedKinematicParticle trailPart = pFactory.particle(trTrailTrk,muon_mass,chi,ndf,muon_sigma);

		SVector3 distanceVector(pcaMuons.first.x()-pcaMuons.second.x(),
					pcaMuons.first.y()-pcaMuons.second.y(),
					pcaMuons.first.z()-pcaMuons.second.z());

		dimuon_dist3D[dimuon_count] = ROOT::Math::Mag(distanceVector);

		std::vector<float> vvv(6);
		vvv[0] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(0,0);
		vvv[1] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(0,1);
		vvv[2] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(1,1);
		vvv[3] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(0,2);
		vvv[4] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(1,2);
		vvv[5] = leadPart->stateAtPoint(pcaMuons.first).kinematicParametersError().matrix()(2,2);
		SMatrixSym3D leadPCACov(vvv.begin(),vvv.end());

		vvv[0] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(0,0);
		vvv[1] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(0,1);
		vvv[2] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(1,1);
		vvv[3] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(0,2);
		vvv[4] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(1,2);
		vvv[5] = trailPart->stateAtPoint(pcaMuons.second).kinematicParametersError().matrix()(2,2);
		SMatrixSym3D trailPCACov(vvv.begin(),vvv.end());


		SMatrixSym3D totCov = leadPCACov + trailPCACov;

		dimuon_dist3DE[dimuon_count] = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/dimuon_dist3D[dimuon_count];

		distanceVector(2) = 0.0;
		dimuon_dist2D[dimuon_count] = ROOT::Math::Mag(distanceVector);
		dimuon_dist2DE[dimuon_count] = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/dimuon_dist2D[dimuon_count];

	      }
	    }

	    dimuon_count++;
	  }
	} */

	muon_count++;

      }
    }
  return muon_count;
}

bool NTupleMaker::GetL1ExtraTriggerMatch(const l1extra::L1JetParticleCollection* l1jets,
					 const l1extra::L1JetParticleCollection* l1taus,
					 const LeafCandidate& leg2)
{
  bool matched = false;
  //check matching to l1tau 44 or l1jet 64
  if(l1taus)
    {
      //check matching with l1tau Pt>44 |eta|<2.172
      matched = false;
      for(unsigned int i=0; i<l1taus->size(); ++i)
	{
	  if( (*l1taus)[i].pt() < 44 || fabs((*l1taus)[i].eta() ) > 2.172 ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR( (*l1taus)[i].p4(), leg2.p4() )  < 0.5 )
	    {
	      matched = true;
	      break;
	    }
	}// for(unsigned int i=0; i<l1taus->size(); ++i)

      if(!matched)
	{
	  if(l1jets){//check matching with l1jet Pt>64 |eta|<2.172
	    for(unsigned int i=0; i < l1jets->size(); ++i)
	      {
		if( (*l1jets)[i].pt() < 64 || fabs((*l1jets)[i].eta() ) > 2.172 ) continue;
		if( ROOT::Math::VectorUtil::DeltaR((*l1jets)[i].p4(), leg2.p4() ) < 0.5 ) {
		  matched = true;
		  break;
		}
	      }//for(unsigned int i=0; i<l1jets->size(); ++i)
	  }
	}
    } //if(l1taus)
  return matched;
}

bool NTupleMaker::GetL1ExtraTriggerMatch(const BXVector<l1t::Jet>* l1jets,
					 const BXVector<l1t::Tau>* l1taus,
					 const LeafCandidate& leg2)
{
  bool matched = false;
  //check matching to l1tau 44 or l1jet 64
  if(l1taus)
    {
      //check matching with l1tau Pt>44 |eta|<2.172
      matched = false;
      for(unsigned int i=0; i<l1taus->size(); ++i)
	{
	  if( (*l1taus)[i].pt() < 44 || fabs((*l1taus)[i].eta() ) > 2.172 ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR( (*l1taus)[i].p4(), leg2.p4() )  < 0.5 )
	    {
	      matched = true;
	      break;
	    }
	}// for(unsigned int i=0; i<l1taus->size(); ++i)

      if(!matched)
	{
	  if(l1jets){//check matching with l1jet Pt>64 |eta|<2.172
	    for(unsigned int i=0; i < l1jets->size(); ++i)
	      {
		if( (*l1jets)[i].pt() < 64 || fabs((*l1jets)[i].eta() ) > 2.172 ) continue;
		if( ROOT::Math::VectorUtil::DeltaR((*l1jets)[i].p4(), leg2.p4() ) < 0.5 ) {
		  matched = true;
		  break;
		}
	      }//for(unsigned int i=0; i<l1jets->size(); ++i)
	  }
	}
    } //if(l1taus)
  return matched;
}

unsigned int NTupleMaker::AddTriggerObjects(const edm::Event& iEvent,
					    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken,
					    const edm::TriggerResults & triggerResults) {

  // trigger objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(TriggerObjectCollectionToken, triggerObjects);

  assert(triggerObjects.isValid());

  /*
  std::cout << "Filters : " << std::endl;
  for (unsigned int n=0; n < run_hltfilters.size(); ++n) {
    std::cout << "  " << run_hltfilters.at(n) << std::endl;
  }
  */

  for (unsigned int iTO=0; iTO<triggerObjects->size(); ++iTO) {
    if (trigobject_count==M_trigobjectmaxcount) {
      cerr << "number of trigger objects > M_trigobjectmaxcount. They are missing." << endl;
      errors |= 1<<5;
      break;
    }

    (*triggerObjects)[iTO].unpackNamesAndLabels( iEvent, triggerResults  );
    //    std::vector<std::string> filterLabels = (*triggerObjects)[iTO].pathNames();
    std::vector<std::string> filterLabels = (*triggerObjects)[iTO].filterLabels();
    bool matchFound = false;
    //    std::cout << "Object " << iTO << " filters : " <<  filterLabels.size() << std::endl;
    //    for (unsigned int i=0; i < filterLabels.size(); ++i) {
    //      TString FilterLabel(filterLabels.at(i));
    //      if (FilterLabel.Contains("hltL1sMu18erTau24erIorMu20erTau24er"))
    //	std::cout << "    " << filterLabels.at(i) << std::endl;
    //    }
    std::vector<bool> passedFilters; passedFilters.clear();
    //    std::vector<std::string> matchedFilters; matchedFilters.clear();
    for (unsigned int n=0; n < run_hltfilters.size(); ++n) {
      TString HltFilter(run_hltfilters.at(n));
      bool thisMatched = false;
      for (unsigned int i=0; i < filterLabels.size(); ++i) {
	TString FilterName(filterLabels.at(i));
	if (HltFilter==FilterName) {
	  matchFound = true;
	  thisMatched = true;
	  //	  matchedFilters.push_back(filterLabels.at(i));
	  break;
	}
      }
      passedFilters.push_back(thisMatched);
    }
    if (matchFound) {
      //      std::cout << "   trigger object " << iTO
      //		<< "   pt = " << (*triggerObjects)[iTO].pt()
      //		<< "   eta = " << (*triggerObjects)[iTO].eta()
      //		<< "   phi = " << (*triggerObjects)[iTO].phi() << std::endl;
      //      for (unsigned int ifilter=0; ifilter<matchedFilters.size(); ++ifilter)
      //	std::cout << "    " << matchedFilters[ifilter] << std::endl;
      for (unsigned int n=0; n < M_hltfiltersmax; ++n) {
	if (n<passedFilters.size())
	  trigobject_filters[trigobject_count][n] = passedFilters.at(n);
	else
	  trigobject_filters[trigobject_count][n] = false;
      }
      trigobject_px[trigobject_count] = (*triggerObjects)[iTO].px();
      trigobject_py[trigobject_count] = (*triggerObjects)[iTO].py();
      trigobject_pz[trigobject_count] = (*triggerObjects)[iTO].pz();
      trigobject_pt[trigobject_count] = (*triggerObjects)[iTO].pt();
      trigobject_eta[trigobject_count] = (*triggerObjects)[iTO].eta();
      trigobject_phi[trigobject_count] = (*triggerObjects)[iTO].phi();
      trigobject_isMuon[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerMuon);
      trigobject_isElectron[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerElectron);
      trigobject_isTau[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerTau);
      trigobject_isJet[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerJet);
      trigobject_isMET[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerMET);
      trigobject_count++;
    }
  }
  return trigobject_count;
}

//bool NTupleMaker::AddPhotons(const edm::Event& iEvent)
//{
	// int NumGood = 0;
	// //edm::Handle<pat::PhotonCollection> Photons;
	// //iEvent.getByLabel(edm::InputTag("patPhotons"), Photons);
	// edm::Handle<PhotonCollection> Photons;
	// iEvent.getByLabel(edm::InputTag("photons"), Photons);

	// if(Photons.isValid())
	// {
	// 	for(unsigned i = 0 ; i < Photons->size() ; i++)
	// 	{
	// 		photon_px[i] = (*Photons)[i].px();
	// 		photon_py[i] = (*Photons)[i].py();
	// 		photon_pz[i] = (*Photons)[i].pz();
	// 		photon_e1x5[i] = (*Photons)[i].e1x5();
	// 		photon_e2x5[i] = (*Photons)[i].e2x5();
	// 		photon_e3x3[i] = (*Photons)[i].e3x3();
	// 		photon_e5x5[i] = (*Photons)[i].e5x5();
	// 		photon_maxenergyxtal[i] = (*Photons)[i].maxEnergyXtal();
	// 		photon_sigmaetaeta[i] = (*Photons)[i].sigmaEtaEta();
	// 		photon_sigmaietaieta[i] = (*Photons)[i].sigmaIetaIeta();
	// 		photon_ehcaloverecal[i] = (*Photons)[i].hadronicOverEm();
	// 		photon_ehcaloverecaldepth1[i] = (*Photons)[i].hadronicDepth1OverEm();
	// 		photon_ehcaloverecaldepth2[i] = (*Photons)[i].hadronicDepth2OverEm();
	// 		photon_isolationr3track[i] = (*Photons)[i].trkSumPtSolidConeDR03();
	// 		photon_isolationr3trackhollow[i] = (*Photons)[i].trkSumPtHollowConeDR03();
	// 		photon_isolationr3ecal[i] = (*Photons)[i].ecalRecHitSumEtConeDR03();
	// 		photon_isolationr3hcal[i] = (*Photons)[i].hcalTowerSumEtConeDR03();
	// 		photon_isolationr3ntrack[i] = (*Photons)[i].nTrkSolidConeDR03();
	// 		photon_isolationr3ntrackhollow[i] = (*Photons)[i].nTrkHollowConeDR03();
	// 		photon_isolationr4track[i] = (*Photons)[i].trkSumPtSolidConeDR04();
	// 		photon_isolationr4trackhollow[i] = (*Photons)[i].trkSumPtHollowConeDR04();
	// 		photon_isolationr4ecal[i] = (*Photons)[i].ecalRecHitSumEtConeDR04();
	// 		photon_isolationr4hcal[i] = (*Photons)[i].hcalTowerSumEtConeDR04();
	// 		photon_isolationr4ntrack[i] = (*Photons)[i].nTrkSolidConeDR04();
	// 		photon_isolationr4ntrackhollow[i] = (*Photons)[i].nTrkHollowConeDR04();
	// 		//			photon_superclusterindex[i] = getSuperCluster((*Photons)[i].superCluster()->energy(), (*Photons)[i].superCluster()->x(), (*Photons)[i].superCluster()->y(), (*Photons)[i].superCluster()->z());
	// 		photon_info[i] = 0;
	// 		photon_info[i] |= (*Photons)[i].isPhoton() << 0;
	// 		photon_info[i] |= (*Photons)[i].hasConversionTracks() << 1;
	// 		photon_info[i] |= (*Photons)[i].hasPixelSeed() << 2;
	// 		photon_gapinfo[i] = 0;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEB() << 0;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEE() << 1;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEBGap() << 2;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEBEtaGap() << 3;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEBPhiGap() << 4;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEEGap() << 5;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEERingGap() << 6;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEEDeeGap() << 7;
	// 		photon_gapinfo[i] |= (*Photons)[i].isEBEEGap() << 8;
	// 		photon_conversionbegin[i] = conversion_count;
	// 		photon_trigger[i] = GetTriggerMatch((*Photons)[i], photontriggers);

	// 		ConversionRefVector conversions = (*Photons)[i].conversions();
	// 		for(unsigned j = 0 ; j < conversions.size() ; j++)
	// 		{
	// 			ConversionRef currconv = conversions[j];

	// 			conversion_info[conversion_count] = 0;
	// 			conversion_info[conversion_count] |= currconv->isConverted() << 0;
	// 			conversion_info[conversion_count] |= currconv->conversionVertex().isValid() << 1;
	// 			conversion_vx[conversion_count] = currconv->conversionVertex().x();
	// 			conversion_vy[conversion_count] = currconv->conversionVertex().y();
	// 			conversion_vz[conversion_count] = currconv->conversionVertex().z();
	// 			conversion_chi2[conversion_count] = currconv->conversionVertex().chi2();
	// 			conversion_ndof[conversion_count] = currconv->conversionVertex().ndof();
	// 			conversion_cov[conversion_count][0] = currconv->conversionVertex().covariance(0,0);
	// 			conversion_cov[conversion_count][1] = currconv->conversionVertex().covariance(0,1);
	// 			conversion_cov[conversion_count][2] = currconv->conversionVertex().covariance(0,2);
	// 			conversion_cov[conversion_count][3] = currconv->conversionVertex().covariance(1,1);
	// 			conversion_cov[conversion_count][4] = currconv->conversionVertex().covariance(1,2);
	// 			conversion_cov[conversion_count][5] = currconv->conversionVertex().covariance(2,2);
	// 			conversion_mvaout[conversion_count] = currconv->MVAout();

	// 			conversion_trackndof[conversion_count][0] = -1.;
	// 			conversion_trackndof[conversion_count][1] = -1.;
	// 			if(currconv->nTracks() == 2)
	// 			{
	// 				edm::RefToBase<Track> trA  = currconv->tracks()[0];
	// 				TransientTrack SVTTrackA = TTrackBuilder->build(*trA);
	// 				TrajectoryStateClosestToPoint TTrackStateA = SVTTrackA.trajectoryStateClosestToPoint(GlobalPoint(currconv->conversionVertex().x(), currconv->conversionVertex().y(), currconv->conversionVertex().z()));
	// 				edm::RefToBase<Track> trB  = currconv->tracks()[1];
	// 				TransientTrack SVTTrackB = TTrackBuilder->build(*trB);
	// 				TrajectoryStateClosestToPoint TTrackStateB = SVTTrackB.trajectoryStateClosestToPoint(GlobalPoint(currconv->conversionVertex().x(), currconv->conversionVertex().y(), currconv->conversionVertex().z()));

	// 				if(TTrackStateB.isValid() && TTrackStateA.isValid())
	// 				{

	// 					conversion_trackecalpointx[conversion_count][0] = currconv->ecalImpactPosition()[0].X();
	// 					conversion_trackecalpointy[conversion_count][0] = currconv->ecalImpactPosition()[0].Y();
	// 					conversion_trackecalpointz[conversion_count][0] = currconv->ecalImpactPosition()[0].Z();
	// 					conversion_trackpx[conversion_count][0] = TTrackStateA.momentum().x();
	// 					conversion_trackpy[conversion_count][0] = TTrackStateA.momentum().y();
	// 					conversion_trackpz[conversion_count][0] = TTrackStateA.momentum().z();
	// 					conversion_trackclosestpointx[conversion_count][0] =  TTrackStateA.position().x();
	// 					conversion_trackclosestpointy[conversion_count][0] =  TTrackStateA.position().y();
	// 					conversion_trackclosestpointz[conversion_count][0] =  TTrackStateA.position().z();
	// 					conversion_trackchi2[conversion_count][0] = currconv->tracks()[0]->chi2();
	// 					conversion_trackndof[conversion_count][0] = currconv->tracks()[0]->ndof();
	// 					conversion_trackdxy[conversion_count][0] = TTrackStateA.perigeeParameters().transverseImpactParameter();
	// 					conversion_trackdxyerr[conversion_count][0] = TTrackStateA.perigeeError().transverseImpactParameterError();
	// 					conversion_trackdz[conversion_count][0] = TTrackStateA.perigeeParameters().longitudinalImpactParameter();
	// 					conversion_trackdzerr[conversion_count][0] = TTrackStateA.perigeeError().longitudinalImpactParameterError();
	// 					conversion_trackcharge[conversion_count][0] = currconv->tracks()[0]->charge();
	// 					conversion_tracknhits[conversion_count][0] = currconv->tracks()[0]->numberOfValidHits();
	// 					conversion_tracknmissinghits[conversion_count][0] = currconv->tracks()[0]->numberOfLostHits();
	// 					conversion_tracknpixelhits[conversion_count][0] = currconv->tracks()[0]->hitPattern().numberOfValidPixelHits();
	// 					conversion_tracknpixellayers[conversion_count][0] = currconv->tracks()[0]->hitPattern().pixelLayersWithMeasurement();
	// 					conversion_tracknstriplayers[conversion_count][0] = currconv->tracks()[0]->hitPattern().stripLayersWithMeasurement();
	// 					conversion_trackecalpointx[conversion_count][1] = currconv->ecalImpactPosition()[1].X();
	// 					conversion_trackecalpointy[conversion_count][1] = currconv->ecalImpactPosition()[1].Y();
	// 					conversion_trackecalpointz[conversion_count][1] = currconv->ecalImpactPosition()[1].Z();
	// 					conversion_trackpx[conversion_count][1] = TTrackStateB.momentum().x();
	// 					conversion_trackpy[conversion_count][1] = TTrackStateB.momentum().y();
	// 					conversion_trackpz[conversion_count][1] = TTrackStateB.momentum().z();
	// 					conversion_trackclosestpointx[conversion_count][1] =  TTrackStateB.position().x();
	// 					conversion_trackclosestpointy[conversion_count][1] =  TTrackStateB.position().y();
	// 					conversion_trackclosestpointz[conversion_count][1] =  TTrackStateB.position().z();
	// 					conversion_trackchi2[conversion_count][1] = currconv->tracks()[1]->chi2();
	// 					conversion_trackndof[conversion_count][1] = currconv->tracks()[1]->ndof();
	// 					conversion_trackdxy[conversion_count][1] = TTrackStateB.perigeeParameters().transverseImpactParameter();
	// 					conversion_trackdxyerr[conversion_count][1] = TTrackStateB.perigeeError().transverseImpactParameterError();
	// 					conversion_trackdz[conversion_count][1] = TTrackStateB.perigeeParameters().longitudinalImpactParameter();
	// 					conversion_trackdzerr[conversion_count][1] = TTrackStateB.perigeeError().longitudinalImpactParameterError();
	// 					conversion_trackcharge[conversion_count][1] = currconv->tracks()[1]->charge();
	// 					conversion_tracknhits[conversion_count][1] = currconv->tracks()[1]->numberOfValidHits();
	// 					conversion_tracknmissinghits[conversion_count][1] = currconv->tracks()[1]->numberOfLostHits();
	// 					conversion_tracknpixelhits[conversion_count][1] = currconv->tracks()[1]->hitPattern().numberOfValidPixelHits();
	// 					conversion_tracknpixellayers[conversion_count][1] = currconv->tracks()[1]->hitPattern().pixelLayersWithMeasurement();
	// 					conversion_tracknstriplayers[conversion_count][1] = currconv->tracks()[1]->hitPattern().stripLayersWithMeasurement();
	// 				}
	// 			}
	// 			conversion_count++;
	// 			if(conversion_count == M_conversionmaxcount){cerr << "number of conversions > M_conversionmaxcount. They are missing." << endl; errors |= 1<<4; break;}
	// 		}
	// 		photon_count++;
	// 		if(photon_count == M_photonmaxcount || conversion_count == M_conversionmaxcount){cerr << "number of photon > M_photonmaxcount. They are missing." << endl; errors |= 1<<3; break;}
	// 		if((*Photons)[i].pt() >= cPhotonPtMin && fabs((*Photons)[i].eta()) <= cPhotonEtaMax) NumGood++;
	// 	}
	// }

	// if(NumGood >= cPhotonNum) return(true);
//	return 0;
//}


unsigned int NTupleMaker::AddTaus(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug) cout<<"inside the AddTaus()"<< endl;
  edm::Handle<BXVector<l1t::Jet> > l1jetsHandle;
  const BXVector<l1t::Jet>* l1jets = 0;

  edm::Handle<BXVector<l1t::Tau> > l1tausHandle;
  const BXVector<l1t::Tau>* l1taus = 0;

  edm::Handle<l1extra::L1JetParticleCollection> l1extrajetsHandle;
  const l1extra::L1JetParticleCollection* l1extrajets = 0;

  edm::Handle<l1extra::L1JetParticleCollection> l1extratausHandle;
  const l1extra::L1JetParticleCollection* l1extrataus = 0;

  //  edm::Handle<pat::PackedCandidateCollection> packedPFCandidates;

  // Obtain Collections
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);

  iEvent.getByToken( L1JetCollectionToken_, l1jetsHandle);
  if( !l1jetsHandle.isValid() ) {
    iEvent.getByToken( L1ExtraJetCollectionToken_, l1extrajetsHandle);
    if ( !l1extrajetsHandle.isValid() )
      edm::LogError("DataNotAvailable") << "No L1CentralJets collection available \n"<<std::endl;
    else
      l1extrajets = l1extrajetsHandle.product();
  }
  else
    l1jets = l1jetsHandle.product();

  iEvent.getByToken( L1TauCollectionToken_, l1tausHandle);
  if( !l1tausHandle.isValid() ){
    iEvent.getByToken( L1ExtraTauCollectionToken_, l1extratausHandle);
    if( !l1extratausHandle.isValid() )
      edm::LogError("DataNotAvailable")  << "No L1TauJets collection available \n";
    else
      l1extrataus = l1extratausHandle.product();
  }
  else
    l1taus = l1tausHandle.product();	

  //  iEvent.getByLabel(edm::InputTag("packedPFCandidates"),packedPFCandidates);
  //  assert(packedPFCandidates.isValid());

  //tau collection
  edm::Handle<pat::TauCollection> Taus;
  iEvent.getByToken(TauCollectionToken_, Taus);
  /*
  edm::Handle<pat::PATTauDiscriminator> mvaIsoRaw;
  iEvent.getByToken(TauMVAIsolationRawToken_,mvaIsoRaw);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoVLoose;
  iEvent.getByToken(TauMVAIsolationVLooseToken_,mvaIsoVLoose);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoLoose;
  iEvent.getByToken(TauMVAIsolationLooseToken_,mvaIsoLoose);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoMedium;
  iEvent.getByToken(TauMVAIsolationMediumToken_,mvaIsoMedium);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoTight;
  iEvent.getByToken(TauMVAIsolationTightToken_,mvaIsoTight);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoVTight;
  iEvent.getByToken(TauMVAIsolationVTightToken_,mvaIsoVTight);
  edm::Handle<pat::PATTauDiscriminator> mvaIsoVVTight;
  iEvent.getByToken(TauMVAIsolationVVTightToken_,mvaIsoVVTight);
  */
  //  std::cout << "PV (x,y,z)=("
  //	    << primvertex_x << ","
  //	    << primvertex_y << ","
  //	    << primvertex_z << ")"
  //	    << std::endl;
  if(Taus.isValid())
    {

      for(unsigned i = 0 ; i < Taus->size() ; i++)
	{

	  if(tau_count == M_taumaxcount) {
	    cerr << "number of taus > M_taumaxcount. They are missing." << endl;
	    errors |= 1<<3;
	    break;
	  }

	  if( (*Taus)[i].pt() < cTauPtMin) continue;
	  if(fabs((*Taus)[i].eta()) > cTauEtaMax ) continue;
	  if((*Taus)[i].tauID("decayModeFinding") < 0.5
	     && (*Taus)[i].tauID("decayModeFindingNewDMs") < 0.5 ) continue; //remove this cut from here OR apply new DMF cuts

	  if(doDebug) cout << "Skimmed events..."<< endl;

//calculation for helix parameters for tau. First get PackedCandidate, from this get reco::Track.
	pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>((*Taus)[i].leadChargedHadrCand().get());
const reco::Track* leadTrk=packedLeadTauCand->bestTrack();

if (leadTrk != nullptr){ 
	  TrackBase::ParameterVector ParamVecTau=leadTrk->parameters();
	  TrackBase::CovarianceMatrix CVMTrackTau=leadTrk->covariance();

	  //storage them
	  for(int index=0; index<ParamVecTau.kSize;index++){
	    tau_helixparameters[tau_count][index]=ParamVecTau[index];
	    for(int index2=0; index2<ParamVecTau.kSize;index2++){
	      tau_helixparameters_covar[tau_count][index][index2]=CVMTrackTau[index][index2];
	    }}

	  //get reference point in track
	  TrackBase::Point RFptTau=leadTrk->referencePoint();
	  tau_referencePoint[tau_count][0]=RFptTau.X();
	  tau_referencePoint[tau_count][1]=RFptTau.Y();
	  tau_referencePoint[tau_count][2]=RFptTau.Z();

	  double magneticField = (TTrackBuilder.product() ? TTrackBuilder.product()->field()->inInverseGeV(GlobalPoint(RFptTau.X(), RFptTau.Y(), RFptTau.Z())).z() : 0.0);
	  tau_Bfield[tau_count]=magneticField;
}
else{
	  for(int index=0; index<5;index++){
	    tau_helixparameters[tau_count][index]=-999;
	    for(int index2=0; index2<5;index2++){
	      tau_helixparameters_covar[tau_count][index][index2]=-999;
	    }}

	  tau_referencePoint[tau_count][0]=-999;
	  tau_referencePoint[tau_count][1]=-999;
	  tau_referencePoint[tau_count][2]=-999;
	  tau_Bfield[tau_count]=-999;
}

	  tau_e[tau_count]                                        = (*Taus)[i].energy();
	  tau_px[tau_count]                                       = (*Taus)[i].px();
	  tau_py[tau_count]                                       = (*Taus)[i].py();
	  tau_pz[tau_count]                                       = (*Taus)[i].pz();
	  tau_pt[tau_count]                                       = (*Taus)[i].pt();
	  tau_phi[tau_count]                                      = (*Taus)[i].phi();
	  tau_eta[tau_count]                                      = (*Taus)[i].eta();
	  tau_mass[tau_count]                                     = (*Taus)[i].mass();
	  tau_charge[tau_count]                                   = (*Taus)[i].charge();

	  tau_signalChargedHadrCands_size[tau_count]              = (*Taus)[i].signalChargedHadrCands().size();
	  tau_signalNeutralHadrCands_size[tau_count]              = (*Taus)[i].signalNeutrHadrCands().size();
	  tau_signalGammaCands_size[tau_count]                    = (*Taus)[i].signalGammaCands().size();

	  tau_isolationChargedHadrCands_size[tau_count]           = (*Taus)[i].isolationChargedHadrCands().size();
	  tau_isolationNeutralHadrCands_size[tau_count]           = (*Taus)[i].isolationNeutrHadrCands().size();
          tau_isolationGammaCands_size[tau_count]                 = (*Taus)[i].isolationGammaCands().size();

	  // tau constituents
	  // Set all entries first -999
	  for(unsigned int m=0; m<50; m++){
	    tau_constituents_px[tau_count][m]     = -9999;
	    tau_constituents_py[tau_count][m]     = -9999;
	    tau_constituents_pz[tau_count][m]     = -9999;
	    tau_constituents_e[tau_count][m]      = -9999;
	    tau_constituents_mass[tau_count][m]   = -9999;
	    tau_constituents_charge[tau_count][m] = -9999;
	    tau_constituents_vx[tau_count][m]     = -9999;
	    tau_constituents_vy[tau_count][m]     = -9999;
	    tau_constituents_vz[tau_count][m]     = -9999;
	    tau_constituents_pdgId[tau_count][m]  = -9999;
	    tau_constituents_lostPixelHits[tau_count][m]  = -9999;
	  }
	  UInt_t tau_constituents_count_ = 0;
	  double momHardestTrack[3];
	  double momMax = -1;
	  for(unsigned int m=0; m<(*Taus)[i].signalCands().size(); m++){
	    tau_constituents_px[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->px();
	    tau_constituents_py[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->py();
	    tau_constituents_pz[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->pz();
	    tau_constituents_e[tau_count][tau_constituents_count_]      = (*Taus)[i].signalCands()[m]->energy();
	    tau_constituents_mass[tau_count][tau_constituents_count_]   = (*Taus)[i].signalCands()[m]->mass();
	    tau_constituents_charge[tau_count][tau_constituents_count_] = (*Taus)[i].signalCands()[m]->charge();
	    tau_constituents_vx[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->vx();
	    tau_constituents_vy[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->vy();
	    tau_constituents_vz[tau_count][tau_constituents_count_]     = (*Taus)[i].signalCands()[m]->vz();
	    tau_constituents_pdgId[tau_count][tau_constituents_count_]  = (*Taus)[i].signalCands()[m]->pdgId();

	    if((*Taus)[i].signalCands()[m]->charge()!=0){
	      if ((*Taus)[i].signalCands()[m]->pt()>momMax) {
		momMax = (*Taus)[i].signalCands()[m]->pt();
		momHardestTrack[0] = (*Taus)[i].signalCands()[m]->px();
		momHardestTrack[1] = (*Taus)[i].signalCands()[m]->py();
		momHardestTrack[2] = (*Taus)[i].signalCands()[m]->pz();
	      }
	      //Pixel his information:
	      // -1: valid hit in 1st pixel barrel layer,
	      //  0: noLostInnerHits - no hit in 1st pixel barrel layer, but not expected there e.g. due to geometry,
	      //  1: one lost hit, 2: two or more lost hits
	      const pat::PackedCandidate* pCand = dynamic_cast<const pat::PackedCandidate*>((*Taus)[i].signalCands()[m].get());
	      if(pCand!=nullptr)
		tau_constituents_lostPixelHits[tau_count][tau_constituents_count_] = pCand->lostInnerHits();
	    }
	    /*
	    std::cout << tau_count << ":" << tau_constituents_count_
		      << "   (vx,vy,vz)=("
		      << tau_constituents_vx[tau_count][tau_constituents_count_] << ","
		      << tau_constituents_vy[tau_count][tau_constituents_count_] << ","
		      << tau_constituents_vz[tau_count][tau_constituents_count_] << ")" << std::endl;
	    */
	    tau_constituents_count_ ++;
	  }

	  // tau vertex
	  tau_vertexx[tau_count] = (*Taus)[i].vertex().x();
	  tau_vertexy[tau_count] = (*Taus)[i].vertex().y();
	  tau_vertexz[tau_count] = (*Taus)[i].vertex().z();

	  double pv[3];
	  pv[0] = primvertex_x;
	  pv[1] = primvertex_y;
	  pv[2] = primvertex_z;
	  double refPoint[3];
	  refPoint[0] = (*Taus)[i].dxy_PCA().X();
	  refPoint[1] = (*Taus)[i].dxy_PCA().Y();
	  refPoint[2] = (*Taus)[i].dxy_PCA().Z();
	  double pca[3];
	  computePCA(pv,refPoint,momHardestTrack,pca);
	  tau_pca2D_x[tau_count] = (*Taus)[i].dxy_PCA().X();
	  tau_pca2D_y[tau_count] = (*Taus)[i].dxy_PCA().Y();
	  tau_pca2D_z[tau_count] = (*Taus)[i].dxy_PCA().Z();
	  tau_pca3D_x[tau_count] = pca[0];
	  tau_pca3D_y[tau_count] = pca[1];
	  tau_pca3D_z[tau_count] = pca[2];
	  /*
	  std::cout << "tau-lepton " << tau_count << " -> PCA(xy) (x,y,z)=("
		    << (*Taus)[i].dxy_PCA().X() << ","
		    << (*Taus)[i].dxy_PCA().Y() << ","
		    << (*Taus)[i].dxy_PCA().Z() << ") : PCA(3D) (x,y,z)=("
		    << pca[0] << ","
		    << pca[1] << ","
		    << pca[2] << ") : vertex (x,y,z)=("
		    << tau_vertexx[tau_count] << ","
		    << tau_vertexy[tau_count] << ","
		    << tau_vertexz[tau_count] << ")" << std::endl;
	  */
	  tau_constituents_count[tau_count] = tau_constituents_count_;

	  // discriminators
	  if(setTauBranches){
	    std::vector<std::pair<std::string, float> > idpairs = (*Taus)[i].tauIDs();
	    for (unsigned int id = 0; id < idpairs.size(); id++){

	      TString name1 = "tau_"; name1+=idpairs[id].first;
	      TString name2 = name1; name2+="[tau_count]/F";
	      TBranch* nb = tree->Branch( name1, tau_ids[id], name2);

	      for(Long64_t ient = 0; ient < tree->GetEntries(); ient++)
		nb->Fill();

	      tauIdIndx.push_back(std::make_pair( idpairs[id].first, id));
	    }

	    setTauBranches = 0;
	  }

	  for(unsigned int id = 0; id < tauIdIndx.size(); id++)
	    tau_ids[tauIdIndx[id].second][tau_count]=(*Taus)[i].tauID(tauIdIndx[id].first);

	  if( ((*Taus)[i].genJet())) {
	    std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(*((*Taus)[i].genJet()));
	    if( genTauDecayMode.find("oneProng0Pi0")!=string::npos )          tau_genDecayMode[tau_count] = 0;
	    else if( genTauDecayMode.find("oneProng1Pi0")!=string::npos )     tau_genDecayMode[tau_count] = 1;
	    else if( genTauDecayMode.find("oneProng2Pi0")!=string::npos )     tau_genDecayMode[tau_count] = 2;
	    else if( genTauDecayMode.find("oneProngOther")!=string::npos )    tau_genDecayMode[tau_count] = 3;
	    else if( genTauDecayMode.find("threeProng0Pi0")!=string::npos )   tau_genDecayMode[tau_count] = 4;
	    else if( genTauDecayMode.find("threeProng1Pi0")!=string::npos )   tau_genDecayMode[tau_count] = 5;
	    else if( genTauDecayMode.find("threeProngOther")!=string::npos )  tau_genDecayMode[tau_count] = 6;
	    else if( genTauDecayMode.find("rare")!=string::npos )             tau_genDecayMode[tau_count] = 7;
	    else   tau_genDecayMode[tau_count] = -99;

	    tau_genDecayMode_name[tau_count] = genTauDecayMode;

	    tau_genjet_px[tau_count]        = (*Taus)[i].genJet()->px();
	    tau_genjet_py[tau_count]        = (*Taus)[i].genJet()->py();
	    tau_genjet_pz[tau_count]        = (*Taus)[i].genJet()->pz();
	    tau_genjet_e[tau_count]         = (*Taus)[i].genJet()->energy();
	  }
	  else {
	    tau_genDecayMode[tau_count]              = -99;
	    tau_genDecayMode_name[tau_count]         = string("");
	    tau_genjet_px[tau_count]                 = 0;
            tau_genjet_py[tau_count]                 = 0;
            tau_genjet_pz[tau_count]                 = 0;
            tau_genjet_e[tau_count]                  = 0;
	  }

	  //add reco decay mode
	  tau_decayMode[tau_count] = (*Taus)[i].decayMode();
	  if ( tau_decayMode[tau_count]==0 )  tau_decayMode_name[tau_count] = string("oneProng0Pi0");
	  else if ( tau_decayMode[tau_count]==1 ) tau_decayMode_name[tau_count] = string("oneProng1Pi0");
	  else if ( tau_decayMode[tau_count]==2 ) tau_decayMode_name[tau_count] = string("oneProng2Pi0");
	  else if ( tau_decayMode[tau_count]<10 ) tau_decayMode_name[tau_count] = string("oneProngOther");
	  else if ( tau_decayMode[tau_count]==10) tau_decayMode_name[tau_count] = string("threeProng0Pi0");
	  else if ( tau_decayMode[tau_count]==11) tau_decayMode_name[tau_count] = string("threeProng1Pi0");
	  else if ( tau_decayMode[tau_count]>11)  tau_decayMode_name[tau_count] = string("threeProngOther");
	  else tau_decayMode_name[tau_count] = string("other");

	  SignedImpactParameter3D signed_ip3D;


	  if((*Taus)[i].leadChargedHadrCand().isNonnull())
	    {
	      tau_leadchargedhadrcand_px[tau_count]   = (*Taus)[i].leadChargedHadrCand()->px();
	      tau_leadchargedhadrcand_py[tau_count]   = (*Taus)[i].leadChargedHadrCand()->py();
	      tau_leadchargedhadrcand_pz[tau_count]   = (*Taus)[i].leadChargedHadrCand()->pz();
	      tau_leadchargedhadrcand_mass[tau_count] = (*Taus)[i].leadChargedHadrCand()->mass();
	      tau_leadchargedhadrcand_id[tau_count]   = (*Taus)[i].leadChargedHadrCand()->pdgId();

	      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>((*Taus)[i].leadChargedHadrCand().get());
	      if(packedLeadTauCand!=nullptr){
		tau_leadchargedhadrcand_dxy[tau_count]   = packedLeadTauCand->dxy();
		tau_leadchargedhadrcand_dz[tau_count]    = packedLeadTauCand->dz();
		tau_leadchargedhadrcand_lostPixelHits[tau_count] = packedLeadTauCand->lostInnerHits();
		if(packedLeadTauCand->vertexRef().key()==0){//the PV
		  //documented at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#PV_Assignment
		  //and DataFormats/PatCandidates/interface/PackedCandidate.h
		  tau_leadchargedhadrcand_pvAssocQ[tau_count] = packedLeadTauCand->pvAssociationQuality();
		} else {
		  tau_leadchargedhadrcand_pvAssocQ[tau_count] = -1;
		}
	      }
	    }
	  else
	    {
	      tau_leadchargedhadrcand_px[tau_count]   = -999;
	      tau_leadchargedhadrcand_py[tau_count]   = -999;
	      tau_leadchargedhadrcand_pz[tau_count]   = -999;
	      tau_leadchargedhadrcand_mass[tau_count] = -999;
	      tau_leadchargedhadrcand_id[tau_count]   = -999;
	      tau_leadchargedhadrcand_dxy[tau_count]  = -999;
	      tau_leadchargedhadrcand_dz[tau_count]   = -999;
	      tau_leadchargedhadrcand_lostPixelHits[tau_count] = -999;
	      tau_leadchargedhadrcand_pvAssocQ[tau_count] = -999;
	    }

	  tau_dxy[tau_count]     = -100.0f;
	  tau_dz[tau_count]      = -100.0f;
	  tau_ip3d[tau_count]    = -1.0f;
	  tau_ip3dSig[tau_count] = -1.0f;
	  tau_dxy[tau_count] = (*Taus)[i].dxy();
	  tau_dxySig[tau_count] = (*Taus)[i].dxy_Sig();
	  tau_ip3d[tau_count] = (*Taus)[i].ip3d();
	  tau_ip3dSig[tau_count] = (*Taus)[i].ip3d_Sig();

	  // tau secondary vertex
	  if( (*Taus)[i].hasSecondaryVertex() ){
	    tau_flightLength[tau_count] = sqrt( (*Taus)[i].flightLength().mag2() );
	    tau_flightLengthSig[tau_count] = (*Taus)[i].flightLengthSig();
	    // (*Taus)[i].secondaryVertex() is a reference (edm::Ref) to
	    //  a vertex collection which is not stored in MiniAOD,
	    //  so the SV of Tau should be refit:
	    //  1) Get tracks form Tau signal charged candidates
	    std::vector<reco::TransientTrack> transTrk;
	    TransientVertex transVtx;
	    const reco::CandidatePtrVector cands = (*Taus)[i].signalChargedHadrCands();
	    for(const auto& cand : cands) {
	      if(cand.isNull()) continue;
	      const pat::PackedCandidate* pCand = dynamic_cast<const pat::PackedCandidate*>(cand.get());
	      if (pCand != nullptr && pCand->hasTrackDetails())
		transTrk.push_back(transTrackBuilder->build(&pCand->pseudoTrack()));
	    }
	    // 2) Fit the secondary vertex
	    bool fitOK(true);
	    KalmanVertexFitter kvf(true);
	    if(transTrk.size() > 1){
	      try{
		transVtx = kvf.vertex(transTrk); //KalmanVertexFitter
	      } catch(...){
		fitOK = false;
	      }
	    } else{
	      fitOK = false;
	    }
	    if(!transVtx.hasRefittedTracks()) fitOK = false;
	    if(transVtx.refittedTracks().size()!=transTrk.size()) fitOK = false;
	    if(fitOK){
	      tau_SV_x[tau_count] = transVtx.position().x();
	      tau_SV_y[tau_count] = transVtx.position().y();
	      tau_SV_z[tau_count] = transVtx.position().z();
	      tau_SV_cov[tau_count][0] = transVtx.positionError().cxx(); // xError()
	      tau_SV_cov[tau_count][1] = transVtx.positionError().cyx();
	      tau_SV_cov[tau_count][2] = transVtx.positionError().czx();
	      tau_SV_cov[tau_count][3] = transVtx.positionError().cyy(); // yError()
	      tau_SV_cov[tau_count][4] = transVtx.positionError().czy();
	      tau_SV_cov[tau_count][5] = transVtx.positionError().czz(); // zError()
	    }
	    else{
	      tau_SV_x[tau_count] = tau_vertexx[tau_count];
	      tau_SV_y[tau_count] = tau_vertexy[tau_count];
	      tau_SV_z[tau_count] = tau_vertexz[tau_count];
	      tau_SV_cov[tau_count][0] = 0.;
	      tau_SV_cov[tau_count][1] = 0.;
	      tau_SV_cov[tau_count][2] = 0.;
	      tau_SV_cov[tau_count][3] = 0.;
	      tau_SV_cov[tau_count][4] = 0.;
	      tau_SV_cov[tau_count][5] = 0.;
	    }
	  }
	  else{
	    tau_flightLength[tau_count] = -1.0f;
	    tau_flightLengthSig[tau_count] = -1.0f;
	    tau_SV_x[tau_count] = tau_vertexx[tau_count];
	    tau_SV_y[tau_count] = tau_vertexy[tau_count];
	    tau_SV_z[tau_count] = tau_vertexz[tau_count];
	    tau_SV_cov[tau_count][0] = 0.;
	    tau_SV_cov[tau_count][1] = 0.;
	    tau_SV_cov[tau_count][2] = 0.;
	    tau_SV_cov[tau_count][3] = 0.;
	    tau_SV_cov[tau_count][4] = 0.;
	    tau_SV_cov[tau_count][5] = 0.;
	  }

	  // l1 match
	  if(l1jets!=0 && l1taus!=0)
	    tau_L1trigger_match[tau_count] = GetL1ExtraTriggerMatch(l1jets, l1taus,  (*Taus)[i] );
	  else
	    tau_L1trigger_match[tau_count] = GetL1ExtraTriggerMatch(l1extrajets, l1extrataus,  (*Taus)[i] );

	  // number of tracks in dR cone

	  tau_genmatch[tau_count] = 0;
	  if(cgen && (!cdata || cembedded)){
	    edm::Handle<reco::GenParticleCollection> GenParticles;
	    iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	    if(GenParticles.isValid()) tau_genmatch[tau_count] = utils_genMatch::genMatch( (*Taus)[i].p4(), *GenParticles);
	  }

	  tau_count++;

	} // for(unsigned i = 0 ; i < Taus->size() ; i++)
    }
  //  std::cout << std::endl;
  return tau_count;
}


// #if 0
// template<typename TCollection>
// const NSVfitEventHypothesisBase* matchSVfitHypothesis(const NSVfitEventHypothesisBase& hypo, const TCollection& coll) //const NSVfitEventHypothesisBaseCollection& coll)
// {
//   for(unsigned int i = 0; i < coll.size(); ++i)
//     {
//       const NSVfitEventHypothesisBase& cand = coll[i];
//       assert(hypo.numResonances() == cand.numResonances());

//       bool particle_match = true;
//       for(unsigned int j = 0; j < hypo.numResonances(); ++j)
// 	{
// 	  const NSVfitResonanceHypothesisBase* hypo_res = hypo.resonance(j);
// 	  const NSVfitResonanceHypothesisBase* cand_res = cand.resonance(j);
// 	  assert(hypo_res && cand_res);

// 	  assert(hypo_res->numDaughters() == cand_res->numDaughters());
// 	  for(unsigned int k = 0; k < hypo_res->numDaughters(); ++k)
// 	    {
// 	      if(hypo_res->daughter(k)->particle() != cand_res->daughter(k)->particle())
// 		particle_match = false;
// 	    }
// 	}

//       if(particle_match) return &cand;
//     }

//   assert(false);
//   return NULL;
// }
// #endif

NTupleMaker::DCA NTupleMaker::calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2)
{
        // TODO: Use the reconstructed decay mode, and then only the mass of the decay products?
	float tauMass = 1.777;
	// TODO: What is this sigma?
	float tauSigma = tauMass*1e-6;

	DCA dca = { -1.0, -1.0f, -1.0f, -1.0f };
	reco::TrackRef track1 = tau1.leadTrack();
	reco::TrackRef track2 = tau2.leadTrack();
	if(track1.isNonnull() && track2.isNonnull())
	{
		reco::TransientTrack transientTrack1 = TTrackBuilder->build(*track1);
		reco::TransientTrack transientTrack2 = TTrackBuilder->build(*track2);
		if(transientTrack1.impactPointTSCP().isValid() && transientTrack2.impactPointTSCP().isValid())
		{
			FreeTrajectoryState state1 = transientTrack1.impactPointTSCP().theState();
			FreeTrajectoryState state2 = transientTrack2.impactPointTSCP().theState();
			TwoTrackMinimumDistance minDist;
			minDist.calculate(state1, state2);
			if(minDist.status())
			{
				const float dist3D = minDist.distance();
				std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
				GlobalPoint pca1 = pcas.first;
				GlobalPoint pca2 = pcas.second;

				ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
				const float twoTauDist3D = ROOT::Math::Mag(distanceVector);
				assert(fabs(dist3D - twoTauDist3D) < 1e-3);

				float chi2 = 0.0f, ndf = 0.0f;
				KinematicParticleFactoryFromTransientTrack pFactory;
				RefCountedKinematicParticle tau1Particle = pFactory.particle(transientTrack1, tauMass, chi2, ndf, tauSigma);
				RefCountedKinematicParticle tau2Particle = pFactory.particle(transientTrack2, tauMass, chi2, ndf, tauSigma);

				float sig[6];
				sig[0] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
				sig[1] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
				sig[2] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
				sig[3] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
				sig[4] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
				sig[5] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);

				sig[0] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
				sig[1] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
				sig[2] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
				sig[3] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
				sig[4] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
				sig[5] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);

				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
				const float twoTauDist3DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist3D;

				distanceVector(2) = 0.0;
				const float twoTauDist2D = ROOT::Math::Mag(distanceVector);
				const float twoTauDist2DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist2D;

				dca.dca2d = twoTauDist2D;
				dca.dca2dErr = twoTauDist2DErr;
				dca.dca3d = twoTauDist3D;
				dca.dca3dErr = twoTauDist3DErr;
			}
		}
	}

	return dca;
}


LorentzVector NTupleMaker::GetRescaledTau(const pat::Tau& tau, double shift)
{
  double scale =  1 + shift;
  if((tau.signalPFChargedHadrCands()).size()==1 && (tau.signalPFGammaCands()).size() <= 0) //1 prong e/mu
    {
      scale = sqrt( tau.energy()*(1+shift) * tau.energy()*(1+shift) - tau.mass()*tau.mass() )/tau.p();
    }
  LorentzVector ShiftedTau( tau.px()*scale , tau.py()*scale, tau.pz()*scale, tau.energy()*(1+shift) );

  if(doDebug)cout<<"scaled tau : "<< ShiftedTau.px() <<", "<< ShiftedTau.py()<< ", "<< ShiftedTau.pz()<<"," << ShiftedTau.E()<< endl;
  if(doDebug) std::cout<<" Rescaling Pat::Tau of DecayMode "<<tau.decayMode()<<" by "<<shift<<"% ==> scale = "<<scale<<std::endl;

  return ShiftedTau;
}


unsigned int NTupleMaker::AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<pat::JetCollection> pfjets;
  iEvent.getByToken(JetCollectionToken_, pfjets);

  edm::Handle<edm::ValueMap<float> > puJetIdMVAFull;
  iEvent.getByLabel(edm::InputTag("pileupJetIdFull","full53xDiscriminant"), puJetIdMVAFull);

  edm::Handle<reco::PFJetCollection> ak4jets;
  iEvent.getByLabel(edm::InputTag("slimmedJets"), ak4jets);

  if(pfjets.isValid())
    {

      for(unsigned i = 0 ; i < pfjets->size() ; i++)
	{

	  if(pfjet_count == M_jetmaxcount){
	    cerr << "number of pfjet > M_jetmaxcount. They are missing." << endl;
	    errors |= 1<<4;
	    break;
	  }

	  if((*pfjets)[i].pt() < cJetPtMin) continue;
	  if(fabs((*pfjets)[i].eta()) > cJetEtaMax) continue;

	  pfjet_e[pfjet_count] = (*pfjets)[i].energy();
	  pfjet_px[pfjet_count] = (*pfjets)[i].px();
	  pfjet_py[pfjet_count] = (*pfjets)[i].py();
	  pfjet_pz[pfjet_count] = (*pfjets)[i].pz();
          pfjet_pt[pfjet_count] = (*pfjets)[i].pt();
          pfjet_eta[pfjet_count] = (*pfjets)[i].eta();
          pfjet_phi[pfjet_count] = (*pfjets)[i].phi();
	  pfjet_neutralhadronicenergy[pfjet_count] = (*pfjets)[i].neutralHadronEnergy();
	  pfjet_chargedhadronicenergy[pfjet_count] = (*pfjets)[i].chargedHadronEnergy();
	  pfjet_neutralemenergy[pfjet_count] = (*pfjets)[i].neutralEmEnergy();
	  pfjet_chargedemenergy[pfjet_count] = (*pfjets)[i].chargedEmEnergy();
	  pfjet_muonenergy[pfjet_count] = (*pfjets)[i].muonEnergy();
	  pfjet_chargedmuonenergy[pfjet_count] = (*pfjets)[i].chargedMuEnergy();
	  pfjet_chargedmulti[pfjet_count] = (*pfjets)[i].chargedMultiplicity();
	  pfjet_neutralmulti[pfjet_count] = (*pfjets)[i].neutralMultiplicity();
	  pfjet_chargedhadronmulti[pfjet_count] = (*pfjets)[i].chargedHadronMultiplicity();

	  pfjet_energycorr[pfjet_count] = -1.;
	  pfjet_energycorr_l1fastjet[pfjet_count] = -1.;
	  pfjet_energycorr_l2relative[pfjet_count] = -1.;
	  pfjet_energycorr_l3absolute[pfjet_count] = -1.;
	  pfjet_energycorr_l2l3residual[pfjet_count] = -1.;

	  if((*pfjets)[i].jecSetsAvailable())
	    {
	      pfjet_energycorr[pfjet_count] = (*pfjets)[i].jecFactor("Uncorrected");
	      pfjet_energycorr_l1fastjet[pfjet_count] = (*pfjets)[i].jecFactor("L1FastJet");
	      pfjet_energycorr_l2relative[pfjet_count] = (*pfjets)[i].jecFactor("L2Relative");
	      pfjet_energycorr_l3absolute[pfjet_count] = (*pfjets)[i].jecFactor("L3Absolute");
	      if (cdata) pfjet_energycorr_l2l3residual[pfjet_count] = (*pfjets)[i].jecFactor("L2L3Residual");
	    }

	  jecUnc->setJetEta(pfjet_eta[pfjet_count]);
	  jecUnc->setJetPt(pfjet_pt[pfjet_count]);
	  pfjet_jecUncertainty[pfjet_count] = jecUnc->getUncertainty(true);
	  pfjet_flavour[pfjet_count] = (*pfjets)[i].partonFlavour();

	  //pileup jet id

	  pfjet_pu_jet_fullId_loose[pfjet_count]  = false;
	  pfjet_pu_jet_fullId_medium[pfjet_count] = false;
	  pfjet_pu_jet_fullId_tight[pfjet_count]  = false;
	  /*
	  pfjet_pu_jet_fullDisc_mva[pfjet_count]  = (*pfjets)[i].userFloat("pileupJetIdUpdated:fullDiscriminant");
	  pfjet_pu_jet_fullId_loose[pfjet_count]  = ( (*pfjets)[i].userInt("pileupJetIdUpdated:fullId") & (1<<2) );
	  pfjet_pu_jet_fullId_medium[pfjet_count] = ( (*pfjets)[i].userInt("pileupJetIdUpdated:fullId") & (1<<1) );
	  pfjet_pu_jet_fullId_tight[pfjet_count]  = ( (*pfjets)[i].userInt("pileupJetIdUpdated:fullId") & (1<<0) );
	  */

	  for(unsigned n = 0 ; n < cBtagDiscriminators.size() ; n++)
	    {
	      pfjet_btag[pfjet_count][n] = -1000;
	      if(cBtagDiscriminators[n] != "F"){
		//		std::cout << " " << cBtagDiscriminators.at(n) << "  : " <<  (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) << std::endl;
		pfjet_btag[pfjet_count][n] = (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) ;
	      }
	    }
	  pfjet_count++;
	}
    }
  return  pfjet_count;
}

unsigned int NTupleMaker::AddPFPuppiJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::JetCollection> pfjetspuppi;
  iEvent.getByToken(PuppiJetCollectionToken_, pfjetspuppi);

  if(pfjetspuppi.isValid())
    {
      for(unsigned i = 0 ; i < pfjetspuppi->size() ; i++)
	     {

         if(pfjetpuppi_count == M_jetmaxcount){
	          cerr << "number of pfjetspuppi > M_jetmaxcount. They are missing." << endl;
	          errors |= 1<<4;
	          break;
	       }
        if((*pfjetspuppi)[i].pt() < cJetPtMin) continue;
	      if(fabs((*pfjetspuppi)[i].eta()) > cJetEtaMax) continue;

	       pfjetpuppi_e[pfjetpuppi_count] = (*pfjetspuppi)[i].energy();
	       pfjetpuppi_px[pfjetpuppi_count] = (*pfjetspuppi)[i].px();
	      pfjetpuppi_py[pfjetpuppi_count] = (*pfjetspuppi)[i].py();
	      pfjetpuppi_pz[pfjetpuppi_count] = (*pfjetspuppi)[i].pz();
        pfjetpuppi_pt[pfjetpuppi_count] = (*pfjetspuppi)[i].pt();
          pfjetpuppi_eta[pfjetpuppi_count] = (*pfjetspuppi)[i].eta();
          pfjetpuppi_phi[pfjetpuppi_count] = (*pfjetspuppi)[i].phi();
	  pfjetpuppi_neutralhadronicenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].neutralHadronEnergy();
	  pfjetpuppi_chargedhadronicenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].chargedHadronEnergy();
	  pfjetpuppi_neutralemenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].neutralEmEnergy();
	  pfjetpuppi_chargedemenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].chargedEmEnergy();
	  pfjetpuppi_muonenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].muonEnergy();
	  pfjetpuppi_chargedmuonenergy[pfjetpuppi_count] = (*pfjetspuppi)[i].chargedMuEnergy();
	  pfjetpuppi_chargedmulti[pfjetpuppi_count] = (*pfjetspuppi)[i].chargedMultiplicity();
	  pfjetpuppi_neutralmulti[pfjetpuppi_count] = (*pfjetspuppi)[i].neutralMultiplicity();
	  pfjetpuppi_chargedhadronmulti[pfjetpuppi_count] = (*pfjetspuppi)[i].chargedHadronMultiplicity();

	  pfjetpuppi_energycorr[pfjetpuppi_count] = -1.;
	  pfjetpuppi_energycorr_l1fastjet[pfjetpuppi_count] = -1.;
	  pfjetpuppi_energycorr_l2relative[pfjetpuppi_count] = -1.;
	  pfjetpuppi_energycorr_l3absolute[pfjetpuppi_count] = -1.;
	  pfjetpuppi_energycorr_l2l3residual[pfjetpuppi_count] = -1.;

	  if((*pfjetspuppi)[i].jecSetsAvailable())
	    {
	      pfjetpuppi_energycorr[pfjetpuppi_count] = (*pfjetspuppi)[i].jecFactor("Uncorrected");
	      pfjetpuppi_energycorr_l1fastjet[pfjetpuppi_count] = (*pfjetspuppi)[i].jecFactor("L1FastJet");
	      pfjetpuppi_energycorr_l2relative[pfjetpuppi_count] = (*pfjetspuppi)[i].jecFactor("L2Relative");
	      pfjetpuppi_energycorr_l3absolute[pfjetpuppi_count] = (*pfjetspuppi)[i].jecFactor("L3Absolute");
	      if (cdata) pfjetpuppi_energycorr_l2l3residual[pfjetpuppi_count] = (*pfjetspuppi)[i].jecFactor("L2L3Residual");
	    }

	  jecUnc->setJetEta(pfjetpuppi_eta[pfjetpuppi_count]);
	  jecUnc->setJetPt(pfjetpuppi_pt[pfjetpuppi_count]);
	  pfjetpuppi_jecUncertainty[pfjetpuppi_count] = jecUnc->getUncertainty(true);
	  pfjetpuppi_flavour[pfjetpuppi_count] = (*pfjetspuppi)[i].partonFlavour();

	  for(unsigned n = 0 ; n < cBtagDiscriminators.size() ; n++)
	    {
	      pfjetpuppi_btag[pfjetpuppi_count][n] = -1000;
	      if(cBtagDiscriminators[n] != "F"){
		//		std::cout << " " << cBtagDiscriminators.at(n) << "  : " <<  (*pfjetspuppi)[i].bDiscriminator(cBtagDiscriminators[n]) << std::endl;
		pfjetpuppi_btag[pfjetpuppi_count][n] = (*pfjetspuppi)[i].bDiscriminator(cBtagDiscriminators[n]) ;
	      }
	    }
	  pfjetpuppi_count++;
	}
    }
  return  pfjetpuppi_count;
}

unsigned int NTupleMaker::AddElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Electron> > Electrons;
  //edm::Handle<edm:View<pat::Electron> > Electrons;
  iEvent.getByToken(ElectronCollectionToken_, Electrons);
        edm::Handle<pat::PackedCandidateCollection> pfcands;
        iEvent.getByToken( PackedCantidateCollectionToken_, pfcands);

	/*if(crecelectrontrigger)
	{
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "NonIsolated"), L1Electrons);
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "Isolated"), L1ElectronsIso);
	}*/

	assert(Electrons.isValid());

	for(unsigned i = 0 ; i < Electrons->size() ; i++){

	  if(electron_count == M_electronmaxcount)
	    {
	      cerr << "number of electron > M_electronmaxcount. They are missing." << endl;
	      errors |= 1<<2;
	      break;
	    }

	  const auto el = Electrons->ptrAt(i);

	  if (el->pt()<cElPtMin) continue;
	  if (fabs(el->eta())>cElEtaMax) continue;

	  // Electron scale and smearing corrections
	  if (applyElectronESShift_) {
	    auto corrP4  = el->p4() * el->userFloat("ecalTrkEnergyPostCorr") / el->energy();
	    electron_px[electron_count] = corrP4.Px();
	    electron_py[electron_count] = corrP4.Py();
	    electron_pz[electron_count] = corrP4.Pz();
	    electron_pt[electron_count] = corrP4.Pt();
      auto corrP4_energyscale_up  = el->p4() * el->userFloat("energyScaleUp") / el->energy();
      auto corrP4_energyscale_down  = el->p4() * el->userFloat("energyScaleDown") / el->energy();
      electron_px_energyscale_up[electron_count] = corrP4_energyscale_up.Px();
      electron_px_energyscale_down[electron_count] = corrP4_energyscale_down.Px();
      electron_py_energyscale_up[electron_count] = corrP4_energyscale_up.Py();
      electron_py_energyscale_down[electron_count] = corrP4_energyscale_down.Py();
      electron_pz_energyscale_up[electron_count] = corrP4_energyscale_up.Pz();
      electron_pz_energyscale_down[electron_count] = corrP4_energyscale_down.Pz();
      electron_pt_energyscale_up[electron_count] = corrP4_energyscale_up.Pt();
      electron_pt_energyscale_down[electron_count] = corrP4_energyscale_down.Pt();
      auto corrP4_energysigma_up  = el->p4() * el->userFloat("energySigmaUp") / el->energy();
      auto corrP4_energysigma_down = el->p4() * el->userFloat("energySigmaDown") / el->energy();
      electron_px_energysigma_up[electron_count] = corrP4_energysigma_up.Px();
      electron_px_energysigma_down[electron_count] = corrP4_energysigma_down.Px();
      electron_py_energysigma_up[electron_count] = corrP4_energysigma_up.Py();
      electron_py_energysigma_down[electron_count] = corrP4_energysigma_down.Py();
      electron_pz_energysigma_up[electron_count] = corrP4_energysigma_up.Pz();
      electron_pz_energysigma_down[electron_count] = corrP4_energysigma_down.Pz();
      electron_pt_energysigma_up[electron_count] = corrP4_energysigma_up.Pt();
      electron_pt_energysigma_down[electron_count] = corrP4_energysigma_down.Pt();
	  }
	  else {
	    electron_px[electron_count] = el->px();
	    electron_py[electron_count] = el->py();
	    electron_pz[electron_count] = el->pz();
	    electron_pt[electron_count] = el->pt();
	  }
	  electron_eta[electron_count] = el->eta();
	  electron_phi[electron_count] = el->phi();
	  electron_charge[electron_count] = el->charge();

	  const pat::Electron &lep = (*Electrons)[i];
          electron_miniISO[electron_count]=getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false);

	  electron_esuperclusterovertrack[electron_count] = el->eSuperClusterOverP();
	  electron_eseedclusterovertrack[electron_count] = el->eSeedClusterOverP();
	  electron_deltaetasuperclustertrack[electron_count] = el->deltaEtaSuperClusterTrackAtVtx();
	  electron_deltaphisuperclustertrack[electron_count] = el->deltaPhiSuperClusterTrackAtVtx();
	  electron_e1x5[electron_count] = el->e1x5();
	  electron_e2x5[electron_count] = el->e2x5Max();
	  electron_e5x5[electron_count] = el->e5x5();
	  electron_sigmaetaeta[electron_count] = el->sigmaEtaEta();
	  electron_sigmaietaieta[electron_count] = el->sigmaIetaIeta();
	  electron_full5x5_sigmaietaieta[electron_count] = el->full5x5_sigmaIetaIeta();
	  electron_ehcaloverecal[electron_count] = el->hcalOverEcal();
	  electron_ehcaloverecaldepth1[electron_count] = el->hcalDepth1OverEcal();
	  electron_ehcaloverecaldepth2[electron_count] = el->hcalDepth2OverEcal();
	  electron_info[electron_count] = 0;
	  electron_info[electron_count] |= el->isElectron();
	  electron_ooemoop[electron_count] = (1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy());

	  electron_superClusterEta[electron_count] = el->superCluster()->eta();
	  electron_superClusterPhi[electron_count] = el->superCluster()->phi();
	  electron_superClusterX[electron_count] = el->superCluster()->x();
	  electron_superClusterY[electron_count] = el->superCluster()->y();
	  electron_superClusterZ[electron_count] = el->superCluster()->z();

	  electron_detaInSeed[electron_count] = el->deltaEtaSuperClusterTrackAtVtx()
	    - el->superCluster()->eta()
	    + el->superCluster()->seed()->eta();

	  electron_he[electron_count] = el->hadronicOverEm();

	  electron_chargedHadIso[electron_count] = el->chargedHadronIso();
	  electron_neutralHadIso[electron_count] = el->neutralHadronIso();
	  electron_photonIso[electron_count] = el->photonIso();
	  electron_puIso[electron_count] = el->puChargedHadronIso();

	  electron_r03_sumChargedHadronPt[electron_count] = el->pfIsolationVariables().sumChargedHadronPt;
	  electron_r03_sumChargedParticlePt[electron_count] = el->pfIsolationVariables().sumChargedParticlePt;
	  electron_r03_sumNeutralHadronEt[electron_count] = el->pfIsolationVariables().sumNeutralHadronEt;
	  electron_r03_sumPhotonEt[electron_count] = el->pfIsolationVariables().sumPhotonEt;
	  electron_r03_sumNeutralHadronEtHighThreshold[electron_count] = el->pfIsolationVariables().sumNeutralHadronEtHighThreshold;
	  electron_r03_sumPhotonEtHighThreshold[electron_count] = el->pfIsolationVariables().sumPhotonEtHighThreshold;
	  electron_r03_sumPUPt[electron_count] = el->pfIsolationVariables().sumPUPt;

	  float  eA = getEffectiveArea( fabs(electron_superClusterEta[electron_count]) );
	  electron_eaIsolation[electron_count] = electron_r03_sumChargedHadronPt[electron_count] +
	    TMath::Max(0.0f,electron_r03_sumNeutralHadronEt[electron_count]+electron_r03_sumPhotonEt[electron_count]-eA*rhoNeutral);

	  electron_gapinfo[electron_count] = 0;
	  electron_gapinfo[electron_count] |= el->isEB() << 0;
	  electron_gapinfo[electron_count] |= el->isEE() << 1;
	  electron_gapinfo[electron_count] |= el->isEBGap() << 2;
	  electron_gapinfo[electron_count] |= el->isEBEtaGap() << 3;
	  electron_gapinfo[electron_count] |= el->isEBPhiGap() << 4;
	  electron_gapinfo[electron_count] |= el->isEEGap() << 5;
	  electron_gapinfo[electron_count] |= el->isEERingGap() << 6;
	  electron_gapinfo[electron_count] |= el->isEEDeeGap() << 7;
	  electron_gapinfo[electron_count] |= el->isEBEEGap() << 8;

	  electron_chargeinfo[electron_count] = 0;
	  if(el->isGsfCtfChargeConsistent()) electron_chargeinfo[electron_count] |= (1 << 0);
	  if(el->isGsfCtfScPixChargeConsistent()) electron_chargeinfo[electron_count] |= (1 << 1);
	  if(el->isGsfScPixChargeConsistent()) electron_chargeinfo[electron_count] |= (1 << 2);

	  electron_fbrems[electron_count] = el->fbrem();
	  electron_numbrems[electron_count] = el->numberOfBrems();

	  GsfTrackRef gsfTr_e = el->gsfTrack();
	  TransientTrack TTrack = TTrackBuilder->build(gsfTr_e);
	  math::XYZPoint ecalPos = PositionOnECalSurface(TTrack);
	  electron_outerx[electron_count] = ecalPos.x();
	  electron_outery[electron_count] = ecalPos.y();
	  electron_outerz[electron_count] = ecalPos.z();
	  //TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));

	  electron_trackchi2[electron_count] = gsfTr_e->chi2();
	  electron_trackndof[electron_count] = gsfTr_e->ndof();
	  electron_vx[electron_count] = gsfTr_e->vx();
	  electron_vy[electron_count] = gsfTr_e->vy();
	  electron_vz[electron_count] = gsfTr_e->vz();

	  electron_nhits[electron_count]        = gsfTr_e->numberOfValidHits();
	  electron_nmissinghits[electron_count] = gsfTr_e->numberOfLostHits();
	  electron_npixelhits[electron_count]   = (gsfTr_e->hitPattern()).numberOfValidPixelHits();
	  electron_npixellayers[electron_count]   = (gsfTr_e->hitPattern()).pixelLayersWithMeasurement();
	  electron_nstriplayers[electron_count]   = (gsfTr_e->hitPattern()).stripLayersWithMeasurement();
	  //electron_nmissinginnerhits[electron_count] = gsfTr_e->trackerExpectedHitsInner().numberOfHits();
	  electron_nmissinginnerhits[electron_count] = gsfTr_e->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

	  electron_dxy[electron_count]          = gsfTr_e->dxy(pv_position);
	  electron_dxyerr[electron_count]       = gsfTr_e->dxyError();
	  electron_dz[electron_count]           = gsfTr_e->dz(pv_position);
	  electron_dzerr[electron_count]        = gsfTr_e->dzError();

	  //	  std::cout << "   dxy = " << electron_dxy[electron_count] << "   dz = " << electron_dz[electron_count] << std::endl;

	  // Electron Ids
     electron_cutId_veto_Summer16[electron_count] = el ->electronID("cutBasedElectronID-Summer16-80X-V1-veto");
     //electron_cutId_loose_Summer16[electron_count] =  el ->electronID("cutBasedElectronID-Summer16-80X-V1-loose");
     //electron_cutId_medium_Summer16[electron_count] =  el ->electronID("cutBasedElectronID-Summer16-80X-V1-medium");
     //electron_cutId_tight_Summer16[electron_count]=  el ->electronID("cutBasedElectronID-Summer16-80X-V1-tight");

     electron_mva_value_Spring16_v1[electron_count] = el ->userFloat("ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values");
     electron_mva_wp90_general_Spring16_v1[electron_count] = el ->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp90");
	  electron_mva_wp80_general_Spring16_v1[electron_count] = el ->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp80");

	  //new for 9.4.0 Fall17 Electron id
	  electron_mva_value_Iso_Fall17_v1[electron_count] = el ->userFloat("ElectronMVAEstimatorRun2Fall17IsoV1Values");
	  electron_mva_value_noIso_Fall17_v1[electron_count] = el ->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values");

	  electron_mva_wp90_Iso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V1-wp90");
	  electron_mva_wp80_Iso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V1-wp80");
	  electron_mva_Loose_Iso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V1-wpLoose");

	  electron_mva_wp90_noIso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V1-wp90");
	  electron_mva_wp80_noIso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V1-wp80");
	  electron_mva_Loose_noIso_Fall17_v1[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V1-wpLoose");

     electron_mva_value_Iso_Fall17_v2[electron_count] = el ->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
     electron_mva_value_noIso_Fall17_v2[electron_count] = el ->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");

	  electron_mva_wp90_Iso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V2-wp90");
	  electron_mva_wp80_Iso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V2-wp80");
	  electron_mva_Loose_Iso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-iso-V2-wpLoose");

	  electron_mva_wp90_noIso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V2-wp90");
	  electron_mva_wp80_noIso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V2-wp80");
	  electron_mva_Loose_noIso_Fall17_v2[electron_count] = el ->electronID("mvaEleID-Fall17-noIso-V2-wpLoose");

	  electron_cutId_veto_Fall17[electron_count] =  el ->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
	  //electron_cutId_loose_Fall17[electron_count] =  el ->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
	  //electron_cutId_medium_Fall17[electron_count] =  el ->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
	  //electron_cutId_tight_Fall17[electron_count] =  el ->electronID("cutBasedElectronID-Fall17-94X-V1-tight");

	  electron_cutId_veto_Fall17V2[electron_count] = el ->electronID("cutBasedElectronID-Fall17-94X-V2-veto");
	  //electron_cutId_loose_Fall17V2[electron_count] = el ->electronID("cutBasedElectronID-Fall17-94X-V2-loose");
	  //electron_cutId_medium_Fall17V2[electron_count] = el ->electronID("cutBasedElectronID-Fall17-94X-V2-medium");
	  //electron_cutId_tight_Fall17V2[electron_count] = el ->electronID("cutBasedElectronID-Fall17-94X-V2-tight");
	  //ending for 9.4.0 electron id

	  electron_pass_conversion[electron_count] = (*Electrons)[i].passConversionVeto();

	  //	  std::cout << "  passed conversion veto = " << electron_pass_conversion[electron_count] << std::endl;

	  electron_genmatch[electron_count] = 0;
	  if(cgen && (!cdata || cembedded)){
	    edm::Handle<reco::GenParticleCollection> GenParticles;
	    iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	    if(GenParticles.isValid())
	      electron_genmatch[electron_count] = utils_genMatch::genMatch(  (*Electrons)[i].p4(), *GenParticles);
	  }

	  electron_count++;
	}

	return electron_count;
}

double NTupleMaker::ComputeDiTauMass(LorentzVector leg1, LorentzVector leg2, LorentzVector met, TMatrixD cov)
{
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  NSVfitStandalone::Vector measuredMET( met.px(), met.py(), 0);
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, leg1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, leg2));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, cov, 0);
  algo.addLogM(false);
  //algo.integrateMarkovChain();
  algo.integrateVEGAS();
  double diTauNSVfitMass_ = algo.getMass();
  return diTauNSVfitMass_;
}


//TLorentzVector NTupleMaker::RecoilCorrectedMET(LorentzVector pfMet_, LorentzVector Leg1p4_, LorentzVector Leg2p4_, const reco::GenParticle *boson_, string sampleName_, int nJets_)
//{
  // double newPfMetPt_ = pfMet_.pt();
  // double newPfMetPhi_ = pfMet_.phi();
  // LorentzVector genVisLeg1_(0, 0, 0, 0);
  // LorentzVector genVisLeg2_(0, 0, 0, 0);
  // LorentzVector diLepton_(0, 0, 0, 0);
  // LorentzVector finalLepLeg1_(0, 0, 0, 0);
  // LorentzVector finalLepLeg2_(0, 0, 0, 0);

  // double u1 = 0; double u2 = 0;

  // //case-I
  // if(abs(boson_->pdgId()) == 23 ||  abs(boson_->pdgId()) ==25 || abs(boson_->pdgId()) ==35 || abs(boson_->pdgId()) == 36 ) {  // zjets
  //   for(unsigned j = 0 ; j < boson_->numberOfDaughters() ; j++) {

  //     const reco::Candidate *lepton = boson_->daughter(j);

  //     //checking daugther
  //     if(lepton->pdgId() == -15) { //for tau-
  // 	genVisLeg1_ = lepton->p4();

  // 	//checking invisilbe decay
  // 	for(int k = 0; k < fabs(lepton->numberOfDaughters()); ++k) {
  // 	  const reco::Candidate *daughter2 = lepton->daughter(k);
  // 	  int daughter2Id     = daughter2->pdgId();
  // 	  int daughter2Status = daughter2->status();

  // 	  if((fabs(daughter2Id)==12 || fabs(daughter2Id)==14 || fabs(daughter2Id)==16) && daughter2Status==1 )
  // 	    genVisLeg1_ -= daughter2->p4();
  // 	} // k-loop
  //     }// tau-

  //     if(lepton->pdgId() == +15) { //for tau+
  // 	genVisLeg1_ = lepton->p4();

  // 	//checking invisilbe decay
  // 	for(int k = 0; k < fabs(lepton->numberOfDaughters()); ++k) {
  // 	  const reco::Candidate *daughter2 = lepton->daughter(k);
  // 	  int daughter2Id     = daughter2->pdgId();
  // 	  int daughter2Status = daughter2->status();

  // 	  if((fabs(daughter2Id)==12 || fabs(daughter2Id)==14 || fabs(daughter2Id)==16) && daughter2Status==1 )
  // 	    genVisLeg1_ -= daughter2->p4();
  // 	} // k-loop
  //     }// tau+
  //   } //#daughter

  //   // finalLeg1, finalLeg2
  //   if(genVisLeg1_.Pt() > 0) {
  //     if(ROOT::Math::VectorUtil::DeltaR(genVisLeg1_, Leg1p4_) < 0.3)
  // 	finalLepLeg1_ = Leg1p4_;
  //     else if(ROOT::Math::VectorUtil::DeltaR(genVisLeg1_, Leg2p4_) < 0.3)
  // 	finalLepLeg1_ = Leg2p4_;
  //     else
  // 	finalLepLeg1_ = genVisLeg1_;
  //   }

  //   if(genVisLeg2_.Pt() > 0) {
  //     if(ROOT::Math::VectorUtil::DeltaR(genVisLeg2_, Leg1p4_) < 0.3)
  // 	finalLepLeg2_ = Leg1p4_;
  //     else if(ROOT::Math::VectorUtil::DeltaR(genVisLeg2_, Leg2p4_) < 0.3)
  // 	finalLepLeg2_ = Leg2p4_;
  //     else
  // 	finalLepLeg2_ = Leg2p4_;
  //   }

  //   if(finalLepLeg1_.Pt() > 0 && finalLepLeg2_.Pt() > 0)
  //     diLepton_ = finalLepLeg1_ + finalLepLeg2_;
  //   else diLepton_ = Leg1p4_ + Leg2p4_;
  // } //zjets

  // //case-II
  // if(abs(boson_->pdgId()) == 24) {
  //   for(unsigned j = 0 ; j < boson_->numberOfDaughters() ; j++) {
  //     const reco::Candidate *lepton = boson_->daughter(j);

  //     //checking daugther
  //     if(abs(lepton->pdgId()) == 15) { //for tau-/tau+
  // 	genVisLeg1_ = lepton->p4();

  // 	//checking invisilbe decay
  // 	for(int k = 0; k < fabs(lepton->numberOfDaughters()); ++k) {
  // 	  const reco::Candidate *daughter2 = lepton->daughter(k);
  // 	  int daughter2Id     = daughter2->pdgId();
  // 	  int daughter2Status = daughter2->status();

  // 	  if((fabs(daughter2Id)==12 || fabs(daughter2Id)==14 || fabs(daughter2Id)==16) && daughter2Status==1 )
  // 	    genVisLeg1_ -= daughter2->p4();
  // 	} // k-loop
  //     }// tau-/tau+

  //     else  if(abs(lepton->pdgId()) == 11 || abs(lepton->pdgId()) == 13) {  // ele, mu
  // 	genVisLeg1_ = lepton->p4();
  //     }
  //   } //#daughter

  //   // finalLeg1
  //   if(genVisLeg1_.Pt() > 0){
  //     if(ROOT::Math::VectorUtil::DeltaR(genVisLeg1_, Leg1p4_ ) < 0.3)
  // 	finalLepLeg1_ = Leg1p4_;
  //     else if(ROOT::Math::VectorUtil::DeltaR(genVisLeg1_, Leg2p4_) < 0.3)
  // 	finalLepLeg1_ = Leg2p4_;
  //     else
  // 	finalLepLeg1_ = Leg1p4_;
  //   }

  //   if(finalLepLeg1_.Pt() > 0)
  //     diLepton_ = finalLepLeg1_;
  //   else diLepton_ = Leg1p4_;
  // } //Wjets

  // // TypeI Correction
  // TLorentzVector corMET_;
  // // cout<<"correctir _ : "<< corrector_ << endl;
  // //corrector_->CorrectType1(newPfMetPt_,newPfMetPhi_, boson_->pt() , boson_->phi(), diLepton_.Pt(), diLepton_.Phi(), u1, u2 , 0 , 0, TMath::Min(nJets_,2) );
  // corMET_.SetPtEtaPhiM(newPfMetPt_, 0, newPfMetPhi_, 0);
  //  return corMET_;
  //}


LorentzVector NTupleMaker::GetShiftedMomentum(const pat::Tau& tau, double shift) {
  // shift in momentum by 1%
  cout<<"=========== inside GetShiftedMomentum, scaling by: "<<shift << endl;
  double shiftP    = 1.;
  //  double shiftMass = 1.;

  LorentzVector tauP4 = tau.p4();

  if ( tau.genJet() && deltaR(tauP4, tau.genJet()->p4()) < 0.5 && tau.genJet()->pt() > 8. ) {
    if((tau.signalPFChargedHadrCands()).size()==1 && (tau.signalPFGammaCands()).size() > 0){    // 1-prong with EM decay
      cout<<"case-I -> 1 prong EM decay"<< endl;
      shiftP = 1 + shift; //New correction for Winter2013
      //      shiftMass = 1 + shift;
    }
    else if((tau.signalPFChargedHadrCands()).size()==1 && (tau.signalPFGammaCands()).size()==0){ // 1-prong with hadronic pions
      cout<<"case-II -> 1 prong nonEM decay"<< endl;
      shiftP = 1 + shift; //New correction for Winter2013
      //      shiftMass = 1.;
    }
    else if((tau.signalPFChargedHadrCands()).size()==3) { //3-prong
      cout<<"case-III -> 3 prong"<< endl;
      shiftP = 1 + shift; //New correction for Winter2013
      //      shiftMass = 1 + shift;
    }
  }

  double scaledPx = tau.px()*shiftP;
  double scaledPy = tau.py()*shiftP;
  double scaledPz = tau.pz()*shiftP;
  double scaledM  = tau.mass()*shiftP;
  double scaledE  = TMath::Sqrt(scaledPx*scaledPx + scaledPy*scaledPy + scaledPz*scaledPz + scaledM*scaledM);


  LorentzVector ShiftedTau(scaledPx, scaledPy, scaledPz, scaledE);
  cout<<"unscaled Tau : "<< tauP4.px() <<", "<< tauP4.py()<< ", "<< tauP4.pz()<<", " << tauP4.E()<< endl;
  cout<<"scaled tau :   "<< ShiftedTau.px() <<", "<< ShiftedTau.py()<< ", "<< ShiftedTau.pz()<<", " << ShiftedTau.E()<< endl;

  return ShiftedTau;
}

reco::GenParticle NTupleMaker::getLastCopy(reco::GenParticle part, edm::Handle<reco::GenParticleCollection> parts_handle) {
  reco::GenParticle last_part;
  reco::GenStatusFlags  flags = part.statusFlags();
  bool isLastCopy = flags.isLastCopy();

  if (!isLastCopy) {
    int pdgid = part.pdgId();
    for (unsigned i = 0; i < part.daughterRefVector().size(); ++i) {
      int daughter_index = static_cast<int>(part.daughterRefVector().at(i).key());
      reco::GenParticle daughter = parts_handle->at(daughter_index);
      int daughter_pdgid = daughter.pdgId();
      if (daughter_pdgid==pdgid) last_part = getLastCopy(daughter, parts_handle);
    }
  }
  else last_part =  part;
  return last_part;
}

void NTupleMaker::getTauDaughters(std::vector<reco::GenParticle> &tau_daughters, unsigned &type, reco::GenParticle tau, edm::Handle<reco::GenParticleCollection> parts_handle){
  for (unsigned i = 0; i < tau.daughterRefVector().size(); ++i) {
    int daughter_index = static_cast<int>(tau.daughterRefVector().at(i).key());
    reco::GenParticle daughter = parts_handle->at(daughter_index);
    int daughter_pdgid = fabs(daughter.pdgId());
    if(daughter_pdgid == 22 || daughter_pdgid == 111 || daughter_pdgid == 211 || daughter_pdgid == 321 || daughter_pdgid == 130 || daughter_pdgid == 310 || daughter_pdgid == 11 || daughter_pdgid == 12 || daughter_pdgid == 13 || daughter_pdgid == 14 || daughter_pdgid == 16){
      if(daughter_pdgid == 11) type = 1;
      if(daughter_pdgid == 13) type = 2;
      tau_daughters.push_back(daughter);
    }
    else getTauDaughters(tau_daughters, type, daughter, parts_handle);
  }
}

TauSpinner::SimpleParticle NTupleMaker::ConvertToSimplePart(reco::GenParticle input_part){
  return TauSpinner::SimpleParticle(input_part.px(), input_part.py(), input_part.pz(), input_part.energy(), input_part.pdgId());
}

unsigned int  NTupleMaker::GetTauSpinnerweights(const edm::Event& event, const edm::EventSetup& setup) {


  //std::vector<std::pair<std::string,double>>theta_vec_ = ICTauSpinnerProducer::SplitString(cTauSpinAngles);
  std::vector<std::pair<std::string,double>> theta_vec_;
  std::stringstream ss(cTauSpinAngles);   
  std::string splitstring;  
  while(std::getline(ss, splitstring, ',')){
    double val = std::stod(splitstring);
    if(splitstring.find(".") != std::string::npos) splitstring.replace(splitstring.find("."),1,"p");
    if(splitstring.find("-") != std::string::npos) splitstring.replace(splitstring.find("-"),1,"minus");
    std::string weight_name = "wt_cp_"+splitstring;    
    theta_vec_.push_back(std::make_pair(weight_name,val)); 
  }

  edm::Handle<reco::GenParticleCollection> parts_handle;
  event.getByToken(GenParticleCollectionToken_, parts_handle);
  
  //reco::GenParticle boson = ICTauSpinnerProducer::getBoson(parts_handle);
  reco::GenParticle boson;
  bool foundBoson = false;
  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::GenParticle const& part = parts_handle->at(i);
    bool isLastCopy = part.isLastCopy();
    int part_pdgid = fabs(part.pdgId());
    if (part_pdgid==25&&isLastCopy){
        boson = part;
        foundBoson = true;
        break;
    }
  }
  if(!foundBoson) {std::cout << "ICTauSpinnerProducer: Gen boson not found. Throwing exception." << std::endl; throw;}

  //std::vector<reco::GenParticle> taus = ICTauSpinnerProducer::getTaus(boson,parts_handle);
  std::vector<reco::GenParticle> taus;
  unsigned nTaus=0;
  for (unsigned i = 0; i < boson.daughterRefVector().size(); ++i) {
    int daughter_index = static_cast<int>(boson.daughterRefVector().at(i).key());
    reco::GenParticle daughter = parts_handle->at(daughter_index);
    int daughter_pdgid = fabs(daughter.pdgId());
    if (daughter_pdgid != 15) continue;
    reco::GenParticle tau = NTupleMaker::getLastCopy(daughter, parts_handle);
    taus.push_back(tau);
    nTaus++;
  }  

  if(nTaus!=2) {std::cout << "ICTauSpinnerProducer: Found " << nTaus << " taus, expected 2. Throwing exception." << std::endl; throw;}

  std::vector<reco::GenParticle> tau1_daughters;
  std::vector<reco::GenParticle> tau2_daughters;
  unsigned type1 = 0; 
  unsigned type2 = 0;
  NTupleMaker::getTauDaughters(tau1_daughters,type1,taus[0],parts_handle);  
  NTupleMaker::getTauDaughters(tau2_daughters,type2,taus[1],parts_handle); 

  TauSpinner::SimpleParticle simple_boson = NTupleMaker::ConvertToSimplePart(boson);
  TauSpinner::SimpleParticle simple_tau1 = NTupleMaker::ConvertToSimplePart(taus[0]);
  TauSpinner::SimpleParticle simple_tau2 = NTupleMaker::ConvertToSimplePart(taus[1]);
  std::vector<TauSpinner::SimpleParticle> simple_tau1_daughters;
  std::vector<TauSpinner::SimpleParticle> simple_tau2_daughters;
  
  for(unsigned i=0; i<tau1_daughters.size(); ++i) simple_tau1_daughters.push_back(NTupleMaker::ConvertToSimplePart(tau1_daughters[i]));
  for(unsigned i=0; i<tau2_daughters.size(); ++i) simple_tau2_daughters.push_back(NTupleMaker::ConvertToSimplePart(tau2_daughters[i]));

  for(unsigned i=0; i<theta_vec_.size(); ++i){
    double theta_val_ = theta_vec_[i].second;
    // Can make this more general by having boson pdgid as input or have option for set boson type
    TauSpinner::setHiggsParametersTR(-cos(2*M_PI*theta_val_),cos(2*M_PI*theta_val_),-sin(2*M_PI*theta_val_),-sin(2*M_PI*theta_val_));
    Double_t weight_ = TauSpinner::calculateWeightFromParticlesH(simple_boson,simple_tau1,simple_tau2,simple_tau1_daughters,simple_tau2_daughters); 
   
    TauSpinnerWeight[i]=weight_;	

  }
  return theta_vec_.size();
}

