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
//#include "AnalysisDataFormats/TauAnalysis/interface/PFMEtSignCovMatrix.h"
//#include "RecoJets/JetProducers/interface/PileupJetIdentifier.h"
//#include "DataFormats/METReco/interface/PFMEtSignCovMatrix.h"
//#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
//#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DesyTauAnalyses/CandidateTools/interface/candidateAuxFunctions.h"

#include <TString.h>

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

  if(ilep < 0)
    std::cout<<"WARNING!!: Pairwise MVAMEt: Lep not found, pt = "<<refp4.pt()<<std::endl;
  
  return ilep;
}	


//to set the values from parameter
NTupleMaker::NTupleMaker(const edm::ParameterSet& iConfig) :  
  // data, year, period, skim
  cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
  cYear(iConfig.getUntrackedParameter<unsigned int>("Year")),
  cPeriod(iConfig.getUntrackedParameter<std::string>("Period")),
  cSkim(iConfig.getUntrackedParameter<unsigned int>("Skim")),
  // switches (collections)
  cgen(iConfig.getUntrackedParameter<bool>("GenParticles", false)),
  ctrigger(iConfig.getUntrackedParameter<bool>("Trigger", false)),
  cbeamspot(iConfig.getUntrackedParameter<bool>("BeamSpot", false)),
  crectrack(iConfig.getUntrackedParameter<bool>("RecTrack", false)),
  crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
  crecmuon(iConfig.getUntrackedParameter<bool>("RecMuon", false)),
  crecelectron(iConfig.getUntrackedParameter<bool>("RecElectron", false)),
  crectau(iConfig.getUntrackedParameter<bool>("RecTau", false)),
  crecphoton(iConfig.getUntrackedParameter<bool>("RecPhoton", false)),
  crecpfjet(iConfig.getUntrackedParameter<bool>("RecJet", false)),
  crecpfmet(iConfig.getUntrackedParameter<bool>("RecPFMet", false)),
  crecpfmetcorr(iConfig.getUntrackedParameter<bool>("RecPFMetCorr", false)),
  crecpuppimet(iConfig.getUntrackedParameter<bool>("RecPuppiMet", false)),
  crecmvamet(iConfig.getUntrackedParameter<bool>("RecMvaMet", false)),
  // triggers
  cHLTriggerPaths(iConfig.getUntrackedParameter<vector<string> >("HLTriggerPaths")),
  cTriggerProcess(iConfig.getUntrackedParameter<string>("TriggerProcess", "HLT")),
  // muons
  cMuPtMin(iConfig.getUntrackedParameter<double>("RecMuonPtMin", 10.)),
  cMuEtaMax(iConfig.getUntrackedParameter<double>("RecMuonEtaMax", 2.4)),
  cMuHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecMuonHLTriggerMatching")),
  cMuNum(iConfig.getUntrackedParameter<int>("RecMuonNum", 0)),
  // electrons
  cElPtMin(iConfig.getUntrackedParameter<double>("RecElectronPtMin", 10.)),
  cElEtaMax(iConfig.getUntrackedParameter<double>("RecElectronEtaMax", 2.4)),
  cElHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecElectronHLTriggerMatching")),
  cElNum(iConfig.getUntrackedParameter<int>("RecElectronNum", 0)),
  // taus
  cTauPtMin(iConfig.getUntrackedParameter<double>("RecTauPtMin", 20.)),
  cTauEtaMax(iConfig.getUntrackedParameter<double>("RecTauEtaMax", 2.3)),
  cTauHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecTauHLTriggerMatching")),
  cTauFloatDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecTauFloatDiscriminators")),
  cTauBinaryDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecTauBinaryDiscriminators")),
  cTauNum(iConfig.getUntrackedParameter<int>("RecTauNum", 0)),
  // tracks
  cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.5)),
  cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 2.4)),
  cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
  // photons
  cPhotonPtMin(iConfig.getUntrackedParameter<double>("RecPhotonPtMin", 20.)),
  cPhotonEtaMax(iConfig.getUntrackedParameter<double>("RecPhotonEtaMax", 1000000.)),
  cPhotonHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecPhotonHLTriggerMatching")),
  cPhotonNum(iConfig.getUntrackedParameter<int>("RecPhotonNum", 0)),
  // jets
  cJetPtMin(iConfig.getUntrackedParameter<double>("RecJetPtMin", 30.)),
  cJetEtaMax(iConfig.getUntrackedParameter<double>("RecJetEtaMax", 4.5)),
  cBtagDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecJetBtagDiscriminators")),
  cJetHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecJetHLTriggerMatching")),
  cJetNum(iConfig.getUntrackedParameter<int>("RecJetNum", 0)),
  // collections
  MuonCollectionToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonCollectionTag"))),
  ElectronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("ElectronCollectionTag"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleMvaNonTrigWP80MapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaNonTrigIdWP80Map"))),
  eleMvaNonTrigWP90MapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaNonTrigIdWP90Map"))),
  eleMvaTrigWP80MapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaTrigIdWP80Map"))),
  eleMvaTrigWP90MapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaTrigIdWP90Map"))),
  mvaNonTrigValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaNonTrigValuesMap"))),
  mvaNonTrigCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaNonTrigCategoriesMap"))),
  mvaTrigValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaTrigValuesMap"))),
  mvaTrigCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaTrigCategoriesMap"))),
  TauCollectionToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("TauCollectionTag"))),
  JetCollectionToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag"))),
  MetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCollectionTag"))),
  MetCovMatrixToken_(consumes<CovMatrix2D>(iConfig.getParameter<edm::InputTag>("MetCovMatrixTag"))),
  MetSigToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("MetSigTag"))),
  MetCorrCovMatrixToken_(consumes<CovMatrix2D>(iConfig.getParameter<edm::InputTag>("MetCorrCovMatrixTag"))),
  MetCorrSigToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("MetCorrSigTag"))),  
  MetCorrCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCorrCollectionTag"))),
  PuppiMetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("PuppiMetCollectionTag"))),
  MvaMetCollectionsTag_(iConfig.getParameter<std::vector<edm::InputTag> >("MvaMetCollectionsTag")),

  GenParticleCollectionToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag"))),
  L1JetCollectionToken_(consumes<l1extra::L1JetParticleCollection>(edm::InputTag("l1extraParticles","Central"))),
  L1TauCollectionToken_(consumes<l1extra::L1JetParticleCollection>(edm::InputTag("l1extraParticles","Tau"))),
  PackedCantidateCollectionToken_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
  TriggerObjectCollectionToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjectCollectionTag"))),
  BeamSpotToken_(consumes<BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag"))),
  PVToken_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag"))),
  sampleName(iConfig.getUntrackedParameter<std::string>("SampleName", "Higgs")),
  propagatorWithMaterial(0)
{
  setTauBranches = true;
  
  //  propagatorWithMaterial = NULL;
  if(cYear != 2011 && cYear != 2012 && cYear != 2015)
    throw cms::Exception("NTupleMaker") << "Invalid Year, only 2011, 2012 and 2015  are allowed!";
  if(cPeriod != "Summer11" && cPeriod != "Fall11" && cPeriod != "Summer12" && cPeriod != "PHYS14" && cPeriod != "Spring15" && cPeriod != "Run2015B" && cPeriod != "Run2015C" && cPeriod != "Run2015D")
    throw cms::Exception("NTupleMaker") << "Invalid period, only Summer11, Fall11, Summer12, PHYS14, Spring15, Run2015B, Run2015C and Run2015D are allowed!";
  
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

  std::vector<std::string> myManualCatWeigths;
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml");
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml");
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml");
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml");
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml");
  myManualCatWeigths.push_back("DesyTauAnalyses/NTupleMaker/data/ElectronID/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml");
  
  vector<string> myManualCatWeigthsTrig;
  string the_path;
  for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
    the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
    myManualCatWeigthsTrig.push_back(the_path);
  }
  
  myMVAnonTrigPhys14 = new EGammaMvaEleEstimatorCSA14();
  myMVAnonTrigPhys14->initialize("BDT",
				 EGammaMvaEleEstimatorCSA14::kNonTrigPhys14,
				 true,
				 myManualCatWeigthsTrig);

  string cmsswBase = (getenv ("CMSSW_BASE"));
  jecUnc = new JetCorrectionUncertainty(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JEC/Summer15_25nsV5_MC_Uncertainty_AK4PFchs.txt");  


  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", cTriggerProcess));
  consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));

  consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumes<vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));

  consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));

  for(std::vector<edm::InputTag>::iterator mit = MvaMetCollectionsTag_.begin();
      mit != MvaMetCollectionsTag_.end(); mit++){
    MvaMetCollectionsToken_.push_back(consumes<pat::METCollection>(*mit));
  }
}

//destructor
NTupleMaker::~NTupleMaker(){
  if(propagatorWithMaterial != 0){ delete propagatorWithMaterial;}
  
  delete myMVAnonTrigPhys14;

}


void NTupleMaker::beginJob(){
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("AC1B", "AC1B", 1);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);
  
  tree->Branch("errors", &errors, "errors/i");
  tree->Branch("event_nr", &event_nr, "event_nr/i");
  tree->Branch("event_run", &event_run, "event_run/i");
  tree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
  tree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
  tree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
  tree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
  tree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
  tree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");
  
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
  }  

  // muons
  if (crecmuon) {
    tree->Branch("muon_count", &muon_count, "muon_count/i");
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
    
    tree->Branch("muon_genmatch", muon_genmatch, "muon_genmatch[muon_count]/I");


    tree->Branch("dimuon_count", &dimuon_count, "dimuon_count/i");
    tree->Branch("dimuon_leading", dimuon_leading, "dimuon_leading[dimuon_count]/i");
    tree->Branch("dimuon_trailing", dimuon_trailing, "dimuon_trailing[dimuon_count]/i");   
    tree->Branch("dimuon_dist2D", dimuon_dist2D, "dimuon_dist2D[dimuon_count]/F");
    tree->Branch("dimuon_dist2DE", dimuon_dist2DE, "dimuon_dist2DE[dimuon_count]/F");
    tree->Branch("dimuon_dist3D", dimuon_dist3D, "dimuon_dist3D[dimuon_count]/F");
    tree->Branch("dimuon_dist3DE", dimuon_dist3DE, "dimuon_dist3DE[dimuon_count]/F");
    
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
    // tree->Branch("pfjet_energycorr_l2l3residual", pfjet_energycorr_l2l3residual, "pfjet_energycorr_l2l3residual[pfjet_count]/F");
    // tree->Branch("pfjet_pu_jet_cut_loose", pfjet_pu_jet_cut_loose, "pfjet_pu_jet_cut_loose[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_cut_medium", pfjet_pu_jet_cut_medium, "pfjet_pu_jet_cut_medium[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_cut_tight", pfjet_pu_jet_cut_tight, "pfjet_pu_jet_cut_tight[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_cut_mva", pfjet_pu_jet_cut_mva, "pfjet_pu_jet_cut_mva[pfjet_count]/F");
    // tree->Branch("pfjet_pu_jet_simple_loose", pfjet_pu_jet_simple_loose, "pfjet_pu_jet_simple_loose[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_simple_medium", pfjet_pu_jet_simple_medium, "pfjet_pu_jet_simple_medium[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_simple_tight", pfjet_pu_jet_simple_tight, "pfjet_pu_jet_simple_tight[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_simple_mva", pfjet_pu_jet_simple_mva, "pfjet_pu_jet_simple_mva[pfjet_count]/F");
    // tree->Branch("pfjet_pu_jet_full_loose", pfjet_pu_jet_full_loose, "pfjet_pu_jet_full_loose[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_full_medium", pfjet_pu_jet_full_medium, "pfjet_pu_jet_full_medium[pfjet_count]/O");
    // tree->Branch("pfjet_pu_jet_full_tight", pfjet_pu_jet_full_tight, "pfjet_pu_jet_full_tight[pfjet_count]/O");
    tree->Branch("pfjet_pu_jet_full_mva", pfjet_pu_jet_full_mva, "pfjet_pu_jet_full_mva[pfjet_count]/F");
    tree->Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour[pfjet_count]/I");
    tree->Branch("pfjet_btag", pfjet_btag,"pfjet_btag[pfjet_count][10]/F");
    tree->Branch("pfjet_jecUncertainty",pfjet_jecUncertainty,"pfjet_jecUncertainty[pfjet_count]/F");
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
    tree->Branch("electron_trackchi2", electron_trackchi2, "electron_trackchi2[electron_count]/F");
    tree->Branch("electron_trackndof", electron_trackndof, "electron_trackndof[electron_count]/F");
    tree->Branch("electron_outerx", electron_outerx, "electron_outerx[electron_count]/F");
    tree->Branch("electron_outery", electron_outery, "electron_outery[electron_count]/F");
    tree->Branch("electron_outerz", electron_outerz, "electron_outerz[electron_count]/F");
    tree->Branch("electron_closestpointx", electron_closestpointx, "electron_closestpointx[electron_count]/F");
    tree->Branch("electron_closestpointy", electron_closestpointy, "electron_closestpointy[electron_count]/F");
    tree->Branch("electron_closestpointz", electron_closestpointz, "electron_closestpointz[electron_count]/F");
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

    tree->Branch("electron_mva_id_nontrigPhys14", electron_mva_id_nontrigPhys14, "electron_mva_id_nontrigPhys14[electron_count]/F");
    tree->Branch("electron_mva_value_nontrig_Spring15_v1", electron_mva_value_nontrig_Spring15_v1, "electron_mva_value_nontrig_Spring15_v1[electron_count]/F");
    tree->Branch("electron_mva_value_trig_Spring15_v1", electron_mva_value_trig_Spring15_v1, "electron_mva_value_trig_Spring15_v1[electron_count]/F");
    tree->Branch("electron_mva_category_nontrig_Spring15_v1", electron_mva_category_nontrig_Spring15_v1, "electron_mva_category_nontrig_Spring15_v1[electron_count]/I");
    tree->Branch("electron_mva_category_trig_Spring15_v1", electron_mva_category_trig_Spring15_v1, "electron_mva_category_trig_Spring15_v1[electron_count]/I");

    tree->Branch("electron_mva_wp80_nontrig_Spring15_v1", electron_mva_wp80_nontrig_Spring15_v1, "electron_mva_wp80_nontrig_Spring15_v1[electron_count]/O");
    tree->Branch("electron_mva_wp90_nontrig_Spring15_v1", electron_mva_wp90_nontrig_Spring15_v1, "electron_mva_wp90_nontrig_Spring15_v1[electron_count]/O");
    tree->Branch("electron_mva_wp80_trig_Spring15_v1", electron_mva_wp80_trig_Spring15_v1, "electron_mva_wp80_trig_Spring15_v1[electron_count]/O");
    tree->Branch("electron_mva_wp90_trig_Spring15_v1", electron_mva_wp90_trig_Spring15_v1, "electron_mva_wp90_trig_Spring15_v1[electron_count]/O");

    tree->Branch("electron_cutId_veto_Spring15", electron_cutId_veto_Spring15, "electron_cutId_veto_Spring15[electron_count]/O");
    tree->Branch("electron_cutId_loose_Spring15", electron_cutId_loose_Spring15, "electron_cutId_loose_Spring15[electron_count]/O");
    tree->Branch("electron_cutId_medium_Spring15", electron_cutId_medium_Spring15, "electron_cutId_medium_Spring15[electron_count]/O");
    tree->Branch("electron_cutId_tight_Spring15", electron_cutId_tight_Spring15, "electron_cutId_tight_Spring15[electron_count]/O");

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

    tree->Branch("tau_dxy", tau_dxy, "tau_dxy[tau_count]/F");
    tree->Branch("tau_dz", tau_dz, "tau_dz[tau_count]/F");
    tree->Branch("tau_ip3d", tau_ip3d, "tau_ip3d[tau_count]/F");
    tree->Branch("tau_ip3dSig", tau_ip3dSig, "tau_ip3dSig[tau_count]/F");
    tree->Branch("tau_charge", tau_charge, "tau_charge[tau_count]/F");
    
    tree->Branch("tau_genjet_px", tau_genjet_px, "tau_genjet_px[tau_count]/F");
    tree->Branch("tau_genjet_py", tau_genjet_py, "tau_genjet_py[tau_count]/F");
    tree->Branch("tau_genjet_pz", tau_genjet_pz, "tau_genjet_pz[tau_count]/F");
    tree->Branch("tau_genjet_e", tau_genjet_e, "tau_genjet_e[tau_count]/F");
    tree->Branch("tau_genmatch", tau_genmatch, "tau_genmatch[tau_count]/I");

    tree->Branch("tau_ntracks_pt05", tau_ntracks_pt05, "tau_ntracks_pt05[tau_count]/i");
    tree->Branch("tau_ntracks_pt08", tau_ntracks_pt05, "tau_ntracks_pt05[tau_count]/i");
    tree->Branch("tau_ntracks_pt1",  tau_ntracks_pt1,  "tau_ntracks_pt1[tau_count]/i");
    
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
  if (cgen && !cdata) {
    tree->Branch("genweight", &genweight, "genweight/F");
    tree->Branch("genid1", &genid1, "genid1/F");
    tree->Branch("genx1", &genx1, "genx1/F");
    tree->Branch("genid2", &genid2, "genid2/F");
    tree->Branch("genx2", &genx2, "genx2/F");
    tree->Branch("genScale", &genScale, "genScale/F");
    
    tree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
    tree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
    tree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
    tree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");

    // generated taus
    tree->Branch("gentau_count", &gentau_count, "gentau_count/i");
    tree->Branch("gentau_e",  gentau_e,  "tau_e[gentau_count]/F");
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
    tree->Branch("gentau_visibleNoLep_px", gentau_visibleNoLep_px, "tau_visibleNoLep_px[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_py", gentau_visibleNoLep_py, "tau_visibleNoLep_py[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_pz", gentau_visibleNoLep_pz, "tau_visibleNoLep_pz[gentau_count]/F");

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
  }    

  // trigger objects
  if (ctrigger) {
    tree->Branch("trigobject_count",&trigobject_count,"trigobject_count/i");
    tree->Branch("trigobject_px",trigobject_px,"trigobject_px[trigobject_count]/F");
    tree->Branch("trigobject_py",trigobject_py,"trigobject_py[trigobject_count]/F");
    tree->Branch("trigobject_pz",trigobject_pz,"trigobject_pz[trigobject_count]/F");
    tree->Branch("trigobject_pt",trigobject_pt,"trigobject_pt[trigobject_count]/F");
    tree->Branch("trigobject_eta",trigobject_eta,"trigobject_eta[trigobject_count]/F");
    tree->Branch("trigobject_phi",trigobject_phi,"trigobject_phi[trigobject_count]/F");
    tree->Branch("trigobject_filters",trigobject_filters,"trigobject_filters[trigobject_count][50]/O");
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
	//if(boost::regex_match(HLTConfiguration.triggerName(i).c_str(), what, muonregexes[j].first) && muontriggers.size() < 32)
	{
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
  //HLTPrescaleConfig->init(iRun, iSetup, cTriggerProcess, changed);

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
  
  AddTriggerList(iRun, HLTConfiguration, muonregexes,     muontriggers,     allmuonnames);
  AddTriggerList(iRun, HLTConfiguration, electronregexes, electrontriggers, allelectronnames);
  AddTriggerList(iRun, HLTConfiguration, tauregexes,      tautriggers,      alltaunames);
  AddTriggerList(iRun, HLTConfiguration, photonregexes,   photontriggers,   allphotonnames);
  AddTriggerList(iRun, HLTConfiguration, jetregexes,      jettriggers,      alljetnames);
  
  run_hltnames.clear();
  for (unsigned int i=0; i<cHLTriggerPaths.size(); ++i)
    run_hltnames.push_back(cHLTriggerPaths.at(i));
	
  if(muontriggers.size() > 32) throw cms::Exception("NTupleMaker") << "Too many muon triggers!" << std::endl;
  if(electrontriggers.size() > 32) throw cms::Exception("NTupleMaker") << "Too many electron triggers!" << std::endl;
  if(tautriggers.size() > 32) throw cms::Exception("NTupleMaker") << "Too many tau triggers!" << std::endl;
  if(photontriggers.size() > 32) throw cms::Exception("NTupleMaker") << "Too many photon triggers!" << std::endl;
  if(jettriggers.size() > 32) throw cms::Exception("NTupleMaker") << "Too many jet triggers!" << std::endl;

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
  
  if (run_hltfilters.size()>50) throw cms::Exception("NTupleMaker") << "Too many HLT filters!" << std::endl;

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
  muon_count = 0;
  dimuon_count = 0;
  tau_count = 0;
  gentau_count = 0;
  pfjet_count = 0;
  electron_count = 0;
  photon_count = 0;
  genparticles_count = 0;
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
  assert(L1trigger.isValid());
  
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
  lumi_hltprescaletable = -1;
  
  //HLTriggerResults
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", cTriggerProcess), HLTrigger);
  assert(HLTrigger.isValid());
  for(int i = 0  ; i < 128 ; i++){trigger_HLT[i] = 0;}
  

  //store trigger bits for selected trigger paths
  hltriggerresults_->clear();
  hltriggerprescales_->clear();
  hltriggerresultsV_.clear();
  const edm::TriggerNames& TrigNames_ = iEvent.triggerNames(*HLTrigger);
  for(unsigned i = 0 ; i < HLTrigger->size(); i++)
    {
      if(!HLTrigger->wasrun(i) )continue;
      std::string trigName=TrigNames_.triggerName(i);
      if(cHLTriggerPaths.size() > 0){
  	for(size_t ip = 0; ip < cHLTriggerPaths.size(); ip++){
  	  if(trigName.find(cHLTriggerPaths[ip]) != string::npos){
	    
	    hltriggerprescales_->insert(std::pair<string, int>(trigName, 1.));// FIXME HLTPrescaleConfig->prescaleValue(iEvent,iSetup,trigName)));
  	    hltriggerresults_->insert(std::pair<string, int>(trigName, HLTrigger->accept(i)));
	    TString TriggerName(trigName);
	    //	    std::cout << trigName << " : " 
	    //		      << HLTrigger->accept(i) << " ; prescale : " 
	    //		      << HLTConfiguration.prescaleValue(iEvent,iSetup,trigName) << std::endl;
  	    if(HLTrigger->accept(i)) hltriggerresultsV_.push_back(trigName);
  	  }
  	}
      }
    }
  //  std::cout << std::endl;

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
	  }
	  if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake() && (*Vertex)[i].ndof() >= 4 && (*Vertex)[i].z() > -24 && (*Vertex)[i].z() < 24 && (*Vertex)[i].position().Rho() < 2.)
	    goodprimvertex_count++;
	}
      }
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

      // if (cSkim>0) {
      // 	float pfmet_et = sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      // 	if (pfmet_et<100.) return;
      // 	// else cout << "  PFMet = " << pfmet_et << std::endl;
      // }
      
      edm::Handle<CovMatrix2D> metcov;
      iEvent.getByToken( MetCovMatrixToken_, metcov);
      pfmet_sigxx = (*metcov)(0,0);
      pfmet_sigxy = (*metcov)(0,1);
      pfmet_sigyx = (*metcov)(1,0);
      pfmet_sigyy = (*metcov)(1,1);

      edm::Handle<double> metsig;
      iEvent.getByToken( MetSigToken_, metsig);
      assert(metsig.isValid());
      pfmet_sig = *metsig;

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
      
      edm::Handle<CovMatrix2D> metcov;
      iEvent.getByToken( MetCorrCovMatrixToken_, metcov);
      pfmetcorr_sigxx = (*metcov)(0,0);
      pfmetcorr_sigxy = (*metcov)(0,1);
      pfmetcorr_sigyx = (*metcov)(1,0);
      pfmetcorr_sigyy = (*metcov)(1,1);

      edm::Handle<double> metsig;
      iEvent.getByToken( MetCorrSigToken_, metsig);
      assert(metsig.isValid());
      pfmetcorr_sig = *metsig;

      /*pfmetcorr_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
      pfmetcorr_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
      pfmetcorr_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
      pfmetcorr_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);*/

      pfmetcorr_ex_JetEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_JetEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_JetEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_UnclusteredEnUp = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_UnclusteredEnUp = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);

      pfmetcorr_ex_UnclusteredEnDown = (*patMet)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);
      pfmetcorr_ey_UnclusteredEnDown = (*patMet)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);

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

    } // crecpuppimet
  
  if(doDebug)  cout<<"add MVA MET"<< endl; 
  if(crecmvamet)
    {
      for(std::vector<edm::InputTag>::iterator mit = MvaMetCollectionsTag_.begin();
	
         mit != MvaMetCollectionsTag_.end(); mit++){

	//collect MVA Mets
	edm::Handle<pat::METCollection> imets;
	//edm::Handle<reco::PFMETCollection> imets;
	iEvent.getByLabel(*mit, imets);
        //iEvent.getByToken(*mit, imets);

	if(!imets.isValid()) continue;
	if(imets->size() == 0)continue;
	
	for(std::vector<pat::MET>::const_iterator met = imets->begin(); met != imets->end(); met++){
	  //for(std::vector<reco::PFMET>::const_iterator met = imets->begin(); met != imets->end(); met++){
	  mvamet_ex[mvamet_count] = met->px();
	  mvamet_ey[mvamet_count] = met->py();

	  mvamet_sigxx[mvamet_count] = met->getSignificanceMatrix()(0,0);
	  mvamet_sigxy[mvamet_count] = met->getSignificanceMatrix()(0,1);
	  mvamet_sigyx[mvamet_count] = met->getSignificanceMatrix()(1,0);
	  mvamet_sigyy[mvamet_count] = met->getSignificanceMatrix()(1,1);

	  mvamet_lep1_pt[mvamet_count] = met->userCand("lepton1")->pt();
	  mvamet_lep2_pt[mvamet_count] = met->userCand("lepton2")->pt();
	  
	  if(mit->label().find("MuEle") != std::string::npos){
	    mvamet_channel[mvamet_count] = EMU;
	    mvamet_lep1[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton2")->p4() );
	  }
	  else if (mit->label().find("TauEle") != std::string::npos){
	    mvamet_channel[mvamet_count] = ETAU;
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton2")->p4() );
	  }
	  else if(mit->label().find("TauMu") != std::string::npos){
	    mvamet_channel[mvamet_count] = MUTAU;
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton2")->p4() );
	  }
	  else if(mit->label().find("DiTau") != std::string::npos){
	    mvamet_channel[mvamet_count] = TAUTAU;
	    mvamet_lep1[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(tau_count, tau_px, tau_py, tau_pz, met->userCand("lepton2")->p4() );
	  }
	  else if(mit->label().find("MuMu") != std::string::npos){
	    mvamet_channel[mvamet_count] = MUMU;
	    mvamet_lep1[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(muon_count, muon_px, muon_py, muon_pz, met->userCand("lepton2")->p4() );
	  }
	  else if(mit->label().find("EleEle") != std::string::npos){
	    mvamet_channel[mvamet_count] = EE;
	    mvamet_lep1[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton1")->p4() );
	    mvamet_lep2[mvamet_count] = find_lep(electron_count, electron_px, electron_py, electron_pz, met->userCand("lepton2")->p4() );
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

  // generator info and generated particles 
  if(doDebug)  cout<<"add gen info"<< endl; 
  if(cgen && !cdata)
    {
      AddGenHt(iEvent);

      bool haveGenParticles = AddGenParticles(iEvent);

      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByLabel(edm::InputTag("generator"), HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  //	  cout << "Event weight from HEPMC : " << genweight << endl;
	  genid1 = HEPMC->pdf()->id.first;
	  genx1 = HEPMC->pdf()->x.second;
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
  if (ctrigger) 
    {
      if(doDebug)  cout<<"add trigger info"<< endl; 
      int numberOfTriggerObjects = int(AddTriggerObjects(iEvent));
    } // ctrigger

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

bool NTupleMaker::AddGenHt(const edm::Event& iEvent) {

  genparticles_lheHt = 0.;

  edm::Handle<LHEEventProduct> lheEventProduct;  
  iEvent.getByLabel( "externalLHEProducer", lheEventProduct);
  
  if(!lheEventProduct.isValid())
    return false;

  const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

  size_t numParticles = lheParticles.size();
  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
   int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
   int status = lheEvent.ISTUP[idxParticle];
   if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
       genparticles_lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
   } 
  }

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
	  if(abs((*GenParticles)[i].pdgId()) == 13 /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
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
	  else if(abs((*GenParticles)[i].pdgId()) == 11 /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
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
	      gentau_status[gentau_count] = (*GenParticles)[i].status();

	      gentau_visible_px[gentau_count] = tau_visible_p4.px();
	      gentau_visible_py[gentau_count] = tau_visible_p4.py();
	      gentau_visible_pz[gentau_count] = tau_visible_p4.pz();
	      gentau_visible_e[gentau_count]  = tau_visible_p4.energy();
	      
	      gentau_visible_pt[gentau_count]   = tau_visible_p4.pt();
	      gentau_visible_eta[gentau_count]  = tau_visible_p4.eta();
	      gentau_visible_phi[gentau_count]  = tau_visible_p4.phi();
	      gentau_visible_mass[gentau_count] = tau_visible_p4.mass();

	      gentau_visibleNoLep_px[gentau_count] = tau_visibleNoLep_p4.px();
	      gentau_visibleNoLep_py[gentau_count] = tau_visibleNoLep_p4.py();
	      gentau_visibleNoLep_pz[gentau_count] = tau_visibleNoLep_p4.pz();
	      gentau_visibleNoLep_e[gentau_count]  = tau_visibleNoLep_p4.energy();
	      
	      gentau_visibleNoLep_pt[gentau_count]   = tau_visibleNoLep_p4.pt();
	      gentau_visibleNoLep_eta[gentau_count]  = tau_visibleNoLep_p4.eta();
	      gentau_visibleNoLep_phi[gentau_count]  = tau_visibleNoLep_p4.phi();
	      gentau_visibleNoLep_mass[gentau_count] = tau_visibleNoLep_p4.mass();
	      
	      const GenStatusFlags statusFlags = (*GenParticles)[i].statusFlags();
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
	  else if ( abs((*GenParticles)[i].pdgId()) == 6 && (*GenParticles)[i].status()==62 ) 
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

	      const GenStatusFlags statusFlags = (*GenParticles)[i].statusFlags();
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

unsigned int NTupleMaker::AddMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);  
  const TransientTrackBuilder * transientTrackBuilder = builder.product();

  edm::Handle<pat::MuonCollection> Muons;
  iEvent.getByToken(MuonCollectionToken_, Muons);
  
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken( PackedCantidateCollectionToken_, pfcands);

  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){
	
	if ((*Muons)[i].pt() < cMuPtMin) continue;
	if (fabs(((*Muons)[i].eta()))>cMuEtaMax) continue;

	//	std::cout << "Selected pat::Muon " << i << std::endl;
	
	muon_px[muon_count] = (*Muons)[i].px();
	muon_py[muon_count] = (*Muons)[i].py();
	muon_pz[muon_count] = (*Muons)[i].pz();
	muon_pt[muon_count] = (*Muons)[i].pt();
	muon_eta[muon_count] = (*Muons)[i].eta();
	muon_phi[muon_count] = (*Muons)[i].phi();
	muon_charge[muon_count] = (*Muons)[i].charge();

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
	if(cgen && !cdata){
	  edm::Handle<reco::GenParticleCollection> GenParticles;
	  iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	  if(GenParticles.isValid())
	    muon_genmatch[muon_count] = utils_genMatch::genMatch( (*Muons)[i].p4(), *GenParticles);
	}

	muon_count++;
	
	if( !(*Muons)[i].innerTrack().isNull()){
	  for(unsigned j = i+1 ; j < Muons->size() ; j++){
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
	}

	if (muon_count==M_muonmaxcount) {
	  cerr << "number of muons > M_muonmaxcount. They are missing." << endl; errors |= 1<<1; 
	  break;
	}

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

unsigned int NTupleMaker::AddTriggerObjects(const edm::Event& iEvent) {

  // trigger objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(TriggerObjectCollectionToken_, triggerObjects);
  assert(triggerObjects.isValid());
  
  for (unsigned int iTO=0; iTO<triggerObjects->size(); ++iTO) {
    vector<string> filterLabels = (*triggerObjects)[iTO].filterLabels();
    bool matchFound = false;
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
      for (unsigned int n=0; n < 50; ++n) {
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
      if (trigobject_count==M_trigobjectmaxcount) {
	 cerr << "number of trigger objects > M_trigobjectmaxcount. They are missing." << endl; 
	 errors |= 1<<5; 
	 break;
      }
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
  edm::Handle<l1extra::L1JetParticleCollection> l1jetsHandle;
  const l1extra::L1JetParticleCollection* l1jets = 0;
  
  edm::Handle<l1extra::L1JetParticleCollection> l1tausHandle;
  const l1extra::L1JetParticleCollection* l1taus = 0;

  //  edm::Handle<pat::PackedCandidateCollection> packedPFCandidates;

  // Obtain Collections 
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);

  iEvent.getByToken( L1JetCollectionToken_, l1jetsHandle);
  if( !l1jetsHandle.isValid() )  edm::LogError("DataNotAvailable")  << "No L1CentralJets collection available \n";
  else  l1jets = l1jetsHandle.product();
  
  iEvent.getByToken( L1TauCollectionToken_, l1tausHandle);
  if( !l1jetsHandle.isValid() )  edm::LogError("DataNotAvailable")  << "No L1TauJets collection available \n";
  else  l1taus = l1tausHandle.product(); 
 
  //  iEvent.getByLabel(edm::InputTag("packedPFCandidates"),packedPFCandidates);
  //  assert(packedPFCandidates.isValid());

  //tau collection
  edm::Handle<pat::TauCollection> Taus;
  iEvent.getByToken(TauCollectionToken_, Taus);

  if(Taus.isValid())
    {
      
      // std::cout << "size of tau collection : " << Taus->size() << std::endl;

      // if (Taus->size()>0) {

      // 	std::vector<pat::Tau::IdPair> tauDiscriminatorPairs = (*Taus)[0].tauIDs();
      
      // 	unsigned tauIdSize = tauDiscriminatorPairs.size();
      
      // 	for (unsigned nId = 0; nId < tauIdSize; ++nId) {
      // 	  std::cout << tauDiscriminatorPairs[nId].first << "  :  " << tauDiscriminatorPairs[nId].second << std::endl;
      // 	}

      // }

      for(unsigned i = 0 ; i < Taus->size() ; i++)
	{
	  // std::cout << i << " :  decayModeFinding = " << (*Taus)[i].tauID("decayModeFinding") 
	  //  	    << "   decayModeFindingNewDMs = " << (*Taus)[i].tauID("decayModeFindingNewDMs") 
	  //  	    << "   pt =  " << (*Taus)[i].pt() 
	  // 	    << "  eta = " << (*Taus)[i].eta() 
	  // 	    << "  phi = " << (*Taus)[i].phi() << std::endl;

	  if( (*Taus)[i].pt() < cTauPtMin) continue;
	  if(fabs((*Taus)[i].eta()) > cTauEtaMax ) continue;
	  if((*Taus)[i].tauID("decayModeFinding") < 0.5 
	     && (*Taus)[i].tauID("decayModeFindingNewDMs") < 0.5 ) continue; //remove this cut from here OR apply new DMF cuts

	  if(doDebug) cout << "Skimmed events..."<< endl;

	  tau_e[tau_count]                                        = (*Taus)[i].energy();
	  tau_px[tau_count]                                       = (*Taus)[i].px();
	  tau_py[tau_count]                                       = (*Taus)[i].py();
	  tau_pz[tau_count]                                       = (*Taus)[i].pz();
	  tau_pt[tau_count]                                       = (*Taus)[i].pt();
	  tau_phi[tau_count]                                      = (*Taus)[i].phi();
	  tau_eta[tau_count]                                      = (*Taus)[i].eta();
	  tau_mass[tau_count]                                     = (*Taus)[i].mass();
	  tau_charge[tau_count]                                   = (*Taus)[i].charge();

	  // std::cout << "Tau " << i 
	  // 	    << "   pt = "  << tau_pt[tau_count]
	  //  	    << "   eta = " << tau_eta[tau_count] 
	  //  	    << "   phi = " << tau_phi[tau_count];

	  tau_signalChargedHadrCands_size[tau_count]              = (*Taus)[i].signalChargedHadrCands().size();   
	  tau_signalNeutralHadrCands_size[tau_count]              = (*Taus)[i].signalNeutrHadrCands().size();   
	  tau_signalGammaCands_size[tau_count]                    = (*Taus)[i].signalGammaCands().size();

	  tau_isolationChargedHadrCands_size[tau_count]           = (*Taus)[i].isolationChargedHadrCands().size();
	  tau_isolationNeutralHadrCands_size[tau_count]           = (*Taus)[i].isolationNeutrHadrCands().size();
          tau_isolationGammaCands_size[tau_count]                 = (*Taus)[i].isolationGammaCands().size();

	  // main discriminators
	  if(setTauBranches){
	    std::vector<std::pair<std::string, float> > idpairs = (*Taus)[i].tauIDs();
	    for (unsigned int id = 0; id < idpairs.size(); id++){

	      TString name1 = "tau_"; name1+=idpairs[id].first;
	      TString name2 = name1; name2+="[tau_count]/F";
	      tree->Branch( name1, tau_ids[id], name2);

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

	    //	    std::cout << "   GenDecayMode = " << genTauDecayMode;

	  }
	  else {
	    tau_genDecayMode[tau_count]              = -99;
	    tau_genDecayMode_name[tau_count]         = string("");    
	    tau_genjet_px[tau_count]                 = 0;
            tau_genjet_py[tau_count]                 = 0;
            tau_genjet_pz[tau_count]                 = 0;
            tau_genjet_e[tau_count]                  = 0;

	    //	    std::cout << "   Missing genTau " << std::endl;

	  }

	  //	  std::cout << std::endl;

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

	      //	      std::cout << "leading charged hadron : dxy = " << (*Taus)[i].leadChargedHadrCand()->dxy() 
	      //       		<< "   dz = " << (*Taus)[i].leadChargedHadrCand()->dz() << std::endl;
	      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>((*Taus)[i].leadChargedHadrCand().get());
	      //	      std::cout << "leading charged hadron : dxy = " << packedLeadTauCand->dxy()
	      //			<< "   dz = " << packedLeadTauCand->dz() << std::endl;

	      tau_leadchargedhadrcand_dxy[tau_count]   = packedLeadTauCand->dxy();
	      tau_leadchargedhadrcand_dz[tau_count]   = packedLeadTauCand->dz();


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

	    }
	  
	  // TrackRef track = (*Taus)[i].leadTrack();

	  // std::cout << "Tau dxy = " << (*Taus)[i].dxy() << std::endl;
	  // std::cout << "PV : x = " << pv_position.x() 
	  // 	    << "  y = " << pv_position.y() 
	  // 	    << "  z = " << pv_position.z() << std::endl;
	  // std::cout << "Tau vertex : x = " <<  (*Taus)[i].vertex().x() 
	  //  	    << "  y = " << (*Taus)[i].vertex().y()
	  //  	    << "  z = " << (*Taus)[i].vertex().z() << std::endl;
	  // std::cout << "Leading ch.cand. : x = " << (*Taus)[i].leadChargedHadrCand()->vx()
	  // 	    << "  y = " << (*Taus)[i].leadChargedHadrCand()->vy()
	  // 	    << "  z = " << (*Taus)[i].leadChargedHadrCand()->vz() << std::endl;

	  // if(track.isNonnull())
	  //   {
	  //     const reco::TransientTrack transientTrack = TTrackBuilder->build(*track);
	  //     const Measurement1D meas = signed_ip3D.apply(transientTrack, GlobalVector(track->px(), track->py(), track->pz()), primvertex).second;
	      
	  //     tau_dxy[tau_count] = track->dxy(pv_position);
	  //     tau_dz[tau_count] = track->dz(pv_position);
	  //     tau_ip3d[tau_count] = meas.value();
	  //     tau_ip3dSig[tau_count] = meas.significance();

	  //     //	      std::cout << "Track is present !" << std::endl;

	  //   }
	  // else
	  //   {
	  //     tau_dxy[tau_count]     = -100.0f;
	  //     tau_dz[tau_count]      = -100.0f;
	  //     tau_ip3d[tau_count]    = -1.0f;
	  //     tau_ip3dSig[tau_count] = -1.0f;
	  //     tau_dxy[tau_count] = (*Taus)[i].dxy();
	  //   }

	  tau_dxy[tau_count]     = -100.0f;
	  tau_dz[tau_count]      = -100.0f;
	  tau_ip3d[tau_count]    = -1.0f;
	  tau_ip3dSig[tau_count] = -1.0f;
	  tau_dxy[tau_count] = (*Taus)[i].dxy();

	  // tau vertex
	  tau_vertexx[tau_count] = (*Taus)[i].vertex().x();
	  tau_vertexy[tau_count] = (*Taus)[i].vertex().y();
	  tau_vertexz[tau_count] = (*Taus)[i].vertex().z();

	  // l1 match
	  tau_L1trigger_match[tau_count] = GetL1ExtraTriggerMatch(l1jets, l1taus,  (*Taus)[i] );
	  
	  // number of tracks in dR cone
	  tau_ntracks_pt05[tau_count] = 0;
	  tau_ntracks_pt08[tau_count] = 0;
	  tau_ntracks_pt1[tau_count]  = 0;
	  // for (unsigned int ipf=0; ipf<packedPFCandidates->size(); ++ipf) {
	  //   if (fabs((*packedPFCandidates)[ipf].dz((*Taus)[i].vertex()))<0.2 && 
	  // 	fabs((*packedPFCandidates)[ipf].dxy((*Taus)[i].vertex()))<0.05 ) {
	  //     if(ROOT::Math::VectorUtil::DeltaR((*Taus)[i].p4(),(*packedPFCandidates)[ipf].p4()) < 0.5) {
	  // 	if ((*packedPFCandidates)[ipf].pt()>0.5) tau_ntracks_pt05[tau_count]++; 
	  // 	if ((*packedPFCandidates)[ipf].pt()>0.8) tau_ntracks_pt08[tau_count]++; 	
	  // 	if ((*packedPFCandidates)[ipf].pt()>1.0) tau_ntracks_pt1[tau_count]++; 	
	  //     }
	  //   }
	  // }
	  //	  std::cout << " tracks around :  pt>0.5 = " << tau_ntracks_pt05[tau_count]
	  //		    << "  pt>0.8 = " << tau_ntracks_pt08[tau_count]
	  //		    << "  pt>1.0 = " << tau_ntracks_pt1[tau_count] << std::endl;

	  tau_genmatch[tau_count] = 0;
	  if(cgen && !cdata){
	    edm::Handle<reco::GenParticleCollection> GenParticles;
	    iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	    if(GenParticles.isValid())
	      tau_genmatch[tau_count] = utils_genMatch::genMatch( (*Taus)[i].p4(), *GenParticles);
	  }
	  
	  tau_count++;
	  if(tau_count == M_taumaxcount) {
	    cerr << "number of taus > M_taumaxcount. They are missing." << endl; 
	    errors |= 1<<3; 
	    break;
	  }
	} // for(unsigned i = 0 ; i < Taus->size() ; i++)
    }
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
  
  //	edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > svInfos;
  //	iEvent.getByLabel(edm::InputTag("secondaryVertexTagInfosEI"), svInfos);
  //	assert(svInfos.isValid());
  
  edm::Handle<edm::ValueMap<float> > puJetIdMVAFull;
  //iEvent.getByLabel(edm::InputTag("puJetIdForPFMVAMEt","fullDiscriminant"), puJetIdMVAFull);
  iEvent.getByLabel(edm::InputTag("pileupJetIdFull","full53xDiscriminant"), puJetIdMVAFull);

  edm::Handle<reco::PFJetCollection> ak4jets;
  //iEvent.getByLabel(edm::InputTag("calibratedAK4PFJetsForPFMVAMEt"), ak4jets);
  //iEvent.getByLabel(edm::InputTag("ak4PFJets"), ak4jets);
  //iEvent.getByLabel(edm::InputTag("AK4PFCHS"), ak4jets);
  iEvent.getByLabel(edm::InputTag("slimmedJets"), ak4jets);
  
  //	edm::Handle<edm::ValueMap<int> > puJetIdFlagFull;
  //	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullId"), puJetIdFlagFull);

  if(pfjets.isValid())
    {
      // if (pfjets->size()>0) {
      // 	const std::vector< std::pair< std::string, float > > pairDiscriVector = (*pfjets)[0].getPairDiscri();
      // 	int nDiscri = pairDiscriVector.size();
      // 	std::cout << "Number of discriminators = " << nDiscri << std::endl;
      // 	for (int iD=0;iD<nDiscri;++iD) {
      // 	  std::pair<std::string, float> pairDiscri = pairDiscriVector[iD];
      // 	  std::cout << "Dicsr = " << pairDiscriVector[iD].first << std::endl;
      // 	}
      // }
      // std::cout << std::endl;
      
      
      for(unsigned i = 0 ; i < pfjets->size() ; i++)
	{
	  if((*pfjets)[i].pt() < cJetPtMin) continue;
	  if(fabs((*pfjets)[i].eta()) > cJetEtaMax) continue;
	  //	  std::cout << "Jet  " << i <<  ", pT=" <<  (*pfjets)[i].pt() << std::endl;
	  
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
		    
	  // std::cout << "Jet Energy corrections : " << std::endl;
	  // std::cout << "    L1FastJet    = " << pfjet_energycorr_l1fastjet[pfjet_count] << std::endl;
	  // std::cout << "    L2Relative   = " << pfjet_energycorr_l2relative[pfjet_count] << std::endl;
	  // std::cout << "    L3Absolute   = " << pfjet_energycorr_l3absolute[pfjet_count] << std::endl;
	  // std::cout << "    L2L3Residual = " << pfjet_energycorr_l2l3residual[pfjet_count] << std::endl;
	  // std::cout << "    Total (Uncor)= " << pfjet_energycorr[pfjet_count] << std::endl;
	  
	  
	  //pfjet_pu_jet_simple_loose[pfjet_count] = false;
	  //pfjet_pu_jet_simple_medium[pfjet_count] = false;
	  //pfjet_pu_jet_simple_tight[pfjet_count] = false;
	  //pfjet_pu_jet_simple_mva[pfjet_count] = (*puJetIdMVAFull)[(*pfjets)[i].originalObjectRef()];
	  
	  
	  //		    pfjet_pu_jet_full_loose[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kLoose);
	  //		    pfjet_pu_jet_full_medium[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kMedium);
	  //		    pfjet_pu_jet_full_tight[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kTight);
	  //                pfjet_pu_jet_full_mva[pfjet_count] = (*puJetIdMVAFull)[(*pfjets)[i].originalObjectRef()];
	  
	  //get MVA Id
          //for(reco::PFJetCollection::const_iterator iak4jets = ak4jets->begin(); iak4jets != ak4jets->end(); iak4jets++){
	  //	  pfjet_pu_jet_full_mva[pfjet_count] = -9999;
	  //	  if (puJetIdMVAFull.isValid()&&ak4jets.isValid()) {
	  //	    for(size_t ij = 0; ij < ak4jets->size(); ij++){
	  //	      reco::PFJetRef jetRef (ak4jets, ij);
	  //	      if(deltaR((*pfjets)[i].p4(), jetRef->p4()) < 0.3){
	  //		if(deltaR((*pfjets)[i].p4(), jetRef->p4()) > 0.1)
	  //		  std::cout<<"original jet pt "<<(*pfjets)[i].pt()<<" re-recoed jet pt "<<jetRef->pt()<<" pu mva value "<<(*puJetIdMVAFull)[jetRef]<<std::endl;
	  //		pfjet_pu_jet_full_mva[pfjet_count] = (*puJetIdMVAFull)[jetRef];
	  //	      }
	  //	    }
	  //	  }	  
	  pfjet_pu_jet_full_mva[pfjet_count] = (*pfjets)[i].userFloat("pileupJetId:fullDiscriminant");
	  jecUnc->setJetEta(pfjet_eta[pfjet_count]);
	  jecUnc->setJetPt(pfjet_pt[pfjet_count]);
	  pfjet_jecUncertainty[pfjet_count] = jecUnc->getUncertainty(true);
	  pfjet_flavour[pfjet_count] = (*pfjets)[i].partonFlavour();
		
	  for(unsigned n = 0 ; n < cBtagDiscriminators.size() ; n++)
	    {
	      pfjet_btag[pfjet_count][n] = -1000;
	      if(cBtagDiscriminators[n] != "F"){
		//		std::cout << " " << cBtagDiscriminators.at(n) << "  : " <<  (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) << std::endl;
		pfjet_btag[pfjet_count][n] = (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) ;
	      }
	    }
	  pfjet_count++;
	  if(pfjet_count == M_jetmaxcount){
	    cerr << "number of pfjet > M_jetmaxcount. They are missing." << endl; 
	    errors |= 1<<4; 
	    break;
	  }
	}
    }
  
  return  pfjet_count;
}


unsigned int NTupleMaker::AddElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Electron> > Electrons;
  //edm::Handle<edm:View<pat::Electron> > Electrons;
  iEvent.getByToken(ElectronCollectionToken_, Electrons);
        edm::Handle<pat::PackedCandidateCollection> pfcands;
        iEvent.getByToken( PackedCantidateCollectionToken_, pfcands);

	// cut based
	edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
        iEvent.getByToken(eleVetoIdMapToken_,veto_id_decisions);
        iEvent.getByToken(eleLooseIdMapToken_,loose_id_decisions);
        iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
        iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

	// mva
	edm::Handle<edm::ValueMap<bool> > nontrig_wp80_decisions;
	edm::Handle<edm::ValueMap<bool> > nontrig_wp90_decisions;
	edm::Handle<edm::ValueMap<bool> > trig_wp80_decisions;
	edm::Handle<edm::ValueMap<bool> > trig_wp90_decisions;
        iEvent.getByToken(eleMvaNonTrigWP80MapToken_,nontrig_wp80_decisions);
        iEvent.getByToken(eleMvaNonTrigWP90MapToken_,nontrig_wp90_decisions);
        iEvent.getByToken(eleMvaTrigWP80MapToken_,trig_wp80_decisions);
        iEvent.getByToken(eleMvaTrigWP90MapToken_,trig_wp90_decisions);

	edm::Handle<edm::ValueMap<float> > mvaNonTrigValues;
	edm::Handle<edm::ValueMap<int> > mvaNonTrigCategories;
        iEvent.getByToken(mvaNonTrigValuesMapToken_,mvaNonTrigValues);
        iEvent.getByToken(mvaNonTrigCategoriesMapToken_,mvaNonTrigCategories);

	edm::Handle<edm::ValueMap<float> > mvaTrigValues;
	edm::Handle<edm::ValueMap<int> > mvaTrigCategories;
        iEvent.getByToken(mvaTrigValuesMapToken_,mvaTrigValues);
        iEvent.getByToken(mvaTrigCategoriesMapToken_,mvaTrigCategories);

	/*if(crecelectrontrigger)
	{
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "NonIsolated"), L1Electrons);
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "Isolated"), L1ElectronsIso);
	}*/

	assert(Electrons.isValid());

	for(unsigned i = 0 ; i < Electrons->size() ; i++){

	  const auto el = Electrons->ptrAt(i);

	  if (el->pt()<cElPtMin) continue;
	  if (fabs(el->eta())>cElEtaMax) continue;

	  electron_px[electron_count] = el->px();
	  electron_py[electron_count] = el->py();
	  electron_pz[electron_count] = el->pz();
	  electron_pt[electron_count] = el->pt();
	  electron_eta[electron_count] = el->eta();
	  electron_phi[electron_count] = el->phi(); 
	  electron_charge[electron_count] = el->charge();
	  
	  const pat::Electron &lep = (*Electrons)[i];
          electron_miniISO[electron_count]=getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false);

	  electron_esuperclusterovertrack[electron_count] = el->eSuperClusterOverP();
	  electron_eseedclusterovertrack[electron_count] = el->eSeedClusterOverP();
	  electron_deltaetasuperclustertrack[electron_count] = el->deltaEtaSuperClusterTrackAtVtx();
	  electron_deltaphisuperclustertrack[electron_count] = el->deltaEtaSuperClusterTrackAtVtx();
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
	  electron_closestpointx[electron_count] = gsfTr_e->vx();
	  electron_closestpointy[electron_count] = gsfTr_e->vy();
	  electron_closestpointz[electron_count] = gsfTr_e->vz();
	  
	  electron_nhits[electron_count]        = gsfTr_e->numberOfValidHits();
	  electron_nmissinghits[electron_count] = gsfTr_e->numberOfLostHits();
	  electron_npixelhits[electron_count]   = (gsfTr_e->hitPattern()).numberOfValidPixelHits();
	  electron_npixellayers[electron_count]   = (gsfTr_e->hitPattern()).pixelLayersWithMeasurement();
	  electron_nstriplayers[electron_count]   = (gsfTr_e->hitPattern()).stripLayersWithMeasurement();
	  //electron_nmissinginnerhits[electron_count] = gsfTr_e->trackerExpectedHitsInner().numberOfHits();
	  electron_nmissinginnerhits[electron_count] = gsfTr_e->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

	  electron_dxy[electron_count]          = gsfTr_e->dxy(pv_position);
	  electron_dxyerr[electron_count]       = gsfTr_e->dxyError();
	  electron_dz[electron_count]           = gsfTr_e->dz(pv_position);
	  electron_dzerr[electron_count]        = gsfTr_e->dzError();
	  
	  //	  std::cout << "   dxy = " << electron_dxy[electron_count] << "   dz = " << electron_dz[electron_count] << std::endl;

	  electron_mva_id_nontrigPhys14[electron_count] = myMVAnonTrigPhys14->mvaValue(Electrons->at(i),false);

	  electron_mva_value_nontrig_Spring15_v1[electron_count] = (*mvaNonTrigValues)[el];
	  electron_mva_category_nontrig_Spring15_v1[electron_count] = (*mvaNonTrigCategories)[el];
	  electron_mva_value_trig_Spring15_v1[electron_count] = (*mvaTrigValues)[el];
	  electron_mva_category_trig_Spring15_v1[electron_count] = (*mvaTrigCategories)[el];

          electron_cutId_veto_Spring15[electron_count] = (*veto_id_decisions)[el];
          electron_cutId_loose_Spring15[electron_count] = (*loose_id_decisions)[el];
          electron_cutId_medium_Spring15[electron_count] = (*medium_id_decisions)[el];
          electron_cutId_tight_Spring15[electron_count] = (*tight_id_decisions)[el];

	  electron_mva_wp80_nontrig_Spring15_v1[electron_count] = (*nontrig_wp80_decisions)[el];
	  electron_mva_wp90_nontrig_Spring15_v1[electron_count] = (*nontrig_wp90_decisions)[el];
	  electron_mva_wp80_trig_Spring15_v1[electron_count] = (*trig_wp80_decisions)[el];
	  electron_mva_wp90_trig_Spring15_v1[electron_count] = (*trig_wp90_decisions)[el];
	  
	  electron_pass_conversion[electron_count] = (*Electrons)[i].passConversionVeto();
	  
	  //	  std::cout << "  passed conversion veto = " << electron_pass_conversion[electron_count] << std::endl;

	  electron_genmatch[electron_count] = 0;
	  if(cgen && !cdata){
	    edm::Handle<reco::GenParticleCollection> GenParticles;
	    iEvent.getByToken(GenParticleCollectionToken_, GenParticles);
	    if(GenParticles.isValid())
	      electron_genmatch[electron_count] = utils_genMatch::genMatch(  (*Electrons)[i].p4(), *GenParticles);
	  }
	  
	  electron_count++;

	  if(electron_count == M_electronmaxcount)
	    {
	      cerr << "number of electron > M_electronmaxcount. They are missing." << endl;
	      errors |= 1<<2;
	      break;
	    }
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



