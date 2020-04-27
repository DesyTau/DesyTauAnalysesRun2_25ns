#ifndef NTupleTauIdEfficiencyW_h
#define NTupleTauIdEfficiencyW_h
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <memory>

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

//#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
//#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
//#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HiggsCPinTauDecays/TauRefit/interface/RefitVertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DesyTauAnalyses/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

//includes for Helix parameter calculations
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

//needed for tauspinner
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"

using namespace std;
using namespace reco;

#define M_Nweights 20 //number of GenEventInfo weights to be stored
#define M_trackmaxcount 1000
#define M_superclustermaxcount 1000
#define M_superclustermembermaxcount 1000
#define M_superclusterhitmaxcount 5000
#define M_primvertexmaxcount 1000
#define M_muonmaxcount 1000
#define M_taumaxcount 1000
#define M_electronmaxcount 1000
#define M_photonmaxcount 1000
#define M_conversionmaxcount 1000
#define M_jetmaxcount 1000
#define M_mvametmaxcount 2000
#define M_genparticlesmaxcount 1000
#define M_genjetsmaxcount 1000
#define M_trigobjectmaxcount 1000
#define M_hltfiltersmax 200
#define M_refitvtxmaxcount 1000
#define M_tauspinneranglesmaxcount 10
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double, 2> > CovMatrix2D;

bool doDebug = false;
class NTupleMaker : public edm::EDAnalyzer{
 public:
  explicit NTupleMaker( const edm::ParameterSet& iConfig );
  ~NTupleMaker();

  float getEffectiveArea(float eta) {
    float effArea = 0.2393;
    float absEta = fabs(eta);
    if (absEta<1.0) effArea = 0.1703;
    else if (absEta < 1.4790) effArea = 0.1715;
    else if (absEta < 2.0) effArea = 0.1213;
    else if (absEta < 2.2) effArea = 0.1230;
    else if (absEta < 2.3) effArea = 0.1635;
    else if (absEta < 2.4) effArea = 0.1937;
    return effArea;

  }

 double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,
                        double r_iso_min, double r_iso_max, double kt_scale,
                        bool charged_only) {

    if (ptcl->pt()<5.) return 99999.;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;

      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    double iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
  }


 private:
  enum MotherNames{HIGGS=1, WBOSON, ZBOSON, TAU};
  enum MvaMetChannel{EMU=1, ETAU, MUTAU, TAUTAU, MUMU, EE, UNKNOWN};

  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void endRun(const edm::Run& iRun);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup );

  void computePCA(double * pv, double * refPoint, double * mom, double * pca);

  bool AddLHEInformation(const edm::Event& iEvent);
  bool AddGenParticles(const edm::Event& iEvent);
  bool AddGenJets(const edm::Event& iEvent);
  unsigned int AddElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  //  unsigned int AddPhotons(const edm::Event& iEvent);
  unsigned int AddTaus(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddPFCand(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddPFPuppiJets(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddTriggerObjects(const edm::Event& iEvent, edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken, const edm::TriggerResults & trigRes);
  unsigned int GetTauSpinnerweights(const edm::Event&, const edm::EventSetup&);
  reco::GenParticle getLastCopy(reco::GenParticle part, edm::Handle<reco::GenParticleCollection> parts_handle);
  void getTauDaughters(std::vector<reco::GenParticle> &tau_daughters, unsigned &type, reco::GenParticle tau, edm::Handle<reco::GenParticleCollection> parts_handle);
  TauSpinner::SimpleParticle ConvertToSimplePart(reco::GenParticle input_part);
  bool foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB);
  bool AddSusyInfo(const edm::Event& iEvent);
  bool AddFlags(const edm::Event& iEvent, const char* module, const char* label, const char* process);

  UInt_t GenParticleInfo(const GenParticle* particle);
  bool GetL1ExtraTriggerMatch(const l1extra::L1JetParticleCollection* l1jets,  const l1extra::L1JetParticleCollection* l1taus, const LeafCandidate& leg2);
  bool GetL1ExtraTriggerMatch(const BXVector<l1t::Jet>* l1jets, const BXVector<l1t::Tau>* l1taus, const LeafCandidate& leg2);
  Int_t HasAnyMother(const GenParticle* particle, int id);
  math::XYZPoint PositionOnECalSurface(reco::TransientTrack&);
  double ComputeDiTauMass(LorentzVector leg1, LorentzVector leg2, LorentzVector met, TMatrixD cov);
  //  TLorentzVector RecoilCorrectedMET(LorentzVector pfMet_, LorentzVector Leg1p4_, LorentzVector Leg2p4_, const reco::GenParticle *boson_, string sampleName_, int nJets_);
  LorentzVector GetShiftedMomentum(const pat::Tau& tau, double shift);
  LorentzVector GetRescaledTau(const pat::Tau& tau, double shift);

  Int_t find_lep(const Int_t nlep, const Float_t px[], const Float_t py[], const Float_t pz[], const reco::Candidate::LorentzVector& refp4);

  struct DCA {
    float dca2d;
    float dca2dErr;
    float dca3d;
    float dca3dErr;
  };
  DCA calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2);

  TTree* tree;
  TTree* lumitree;
  TTree* runtree;
  TH1D*  nEvents;

  //Configuration (steering cards)

  bool cdata;
  bool cembedded;
  bool cFastSim;
  unsigned int cYear;
  std::string cPeriod;
  unsigned int cSkim;
  std::string cJECfile;

  bool cgen;
  bool csusyinfo;
  bool ctrigger;
  bool cbeamspot;
  bool crectrack;
  bool crecprimvertex;
  bool crecprimvertexwithbs;
  bool crefittedvertex;
  bool crefittedvertexwithbs;
  bool crecTauSpinner;
  bool crecmuon;
  bool crecelectron;
  bool crectau;
  bool cl1isotau;
  bool cl1objects;
  bool crecphoton;
  bool crecpfjet;
  bool crecpfpuppijet;
  bool crecpfmet;
  bool crecpfmetcorr;
  bool crecpuppimet;
  bool crecmvamet;
  bool crecstxs;

  vector<string> cHLTriggerPaths;
  string cTriggerProcess;
  string cTauSpinAngles;

  string cMyTriggerProcess;
  vector<string> cFlags;
  vector<string> cFlagsProcesses;
  edm::EDGetTokenT<bool> BadChCandFilterToken_;
  edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdate_token;

  double cMuPtMin;
  double cMuEtaMax;
  vector<string> cMuHLTriggerMatching;
  int cMuNum;

  double cElPtMin;
  double cElEtaMax;
  vector<string> cElHLTriggerMatching;
  int cElNum;

  double cTauPtMin;
  double cTauEtaMax;
  vector<string> cTauHLTriggerMatching;
  vector<string> cTauFloatDiscriminators;
  vector<string> cTauBinaryDiscriminators;
  int cTauNum;

  double cTrackPtMin;
  double cTrackEtaMax;
  double cTrackDxyMax;
  double cTrackDzMax;
  int cTrackNum;

  double cPhotonPtMin;
  double cPhotonEtaMax;
  vector<string> cPhotonHLTriggerMatching;
  int cPhotonNum;

  double cJetPtMin;
  double cJetEtaMax;
  vector<string> cBtagDiscriminators;
  vector<string> cJetHLTriggerMatching;
  int cJetNum;

  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
  edm::EDGetTokenT<edm::PtrVector<reco::Muon>> BadGlobalMuonsToken_;
  edm::EDGetTokenT<edm::PtrVector<reco::Muon>> BadDuplicateMuonsToken_;

  // Electron Configuration
  edm::EDGetTokenT<edm::View<pat::Electron> > ElectronCollectionToken_;
  // Apply electron energy scale shift
  bool applyElectronESShift_;

  //// New for Spring16
  edm::EDGetTokenT<pat::TauCollection> TauCollectionToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationRawToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationVLooseToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationLooseToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationMediumToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationTightToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationVTightToken_;
  edm::EDGetTokenT<pat::PATTauDiscriminator> TauMVAIsolationVVTightToken_;

  edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
  edm::EDGetTokenT<pat::JetCollection> PuppiJetCollectionToken_;
  edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
  edm::EDGetTokenT<CovMatrix2D> MetCovMatrixToken_;
  edm::EDGetTokenT<double> MetSigToken_;
  edm::EDGetTokenT<CovMatrix2D> MetCorrCovMatrixToken_;
  edm::EDGetTokenT<double> MetCorrSigToken_;
  edm::EDGetTokenT<pat::METCollection> MetCorrCollectionToken_;
  edm::EDGetTokenT<pat::METCollection> PuppiMetCollectionToken_;
  std::vector<edm::InputTag> MvaMetCollectionsTag_;
  std::vector<edm::EDGetTokenT<pat::METCollection> > MvaMetCollectionsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> GenParticleCollectionToken_;
  edm::EDGetTokenT<reco::GenJetCollection> GenJetCollectionToken_;
  edm::EDGetTokenT<l1extra::L1JetParticleCollection> L1ExtraJetCollectionToken_;
  edm::EDGetTokenT<l1extra::L1JetParticleCollection> L1ExtraTauCollectionToken_;
  edm::EDGetTokenT<BXVector<l1t::Jet> > L1JetCollectionToken_;
  edm::EDGetTokenT<BXVector<l1t::Tau> > L1IsoTauCollectionToken_;
  edm::EDGetTokenT<BXVector<l1t::Muon> > L1MuonCollectionToken_;
  edm::EDGetTokenT<BXVector<l1t::EGamma> > L1EGammaCollectionToken_;
  edm::EDGetTokenT<BXVector<l1t::Tau> > L1TauCollectionToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> PackedCantidateCollectionToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> myTriggerObjectCollectionToken_;
  edm::EDGetTokenT<BeamSpot> BeamSpotToken_;
  edm::EDGetTokenT<VertexCollection> PVToken_;
  edm::EDGetTokenT<RefitVertexCollection> PVwithBSToken_;
  edm::EDGetTokenT<RefitVertexCollection>RefittedPVToken_;
  edm::EDGetTokenT<RefitVertexCollection> RefittedwithBSPVToken_;
  edm::EDGetTokenT<LHEEventProduct> LHEToken_;
  edm::EDGetTokenT<double> SusyMotherMassToken_;
  edm::EDGetTokenT<double> SusyLSPMassToken_;
  edm::EDGetTokenT<HTXS::HiggsClassification> htxsToken_;
  std::string sampleName;

  PropagatorWithMaterial*               propagatorWithMaterial;

  //Variables
  edm::ESHandle<TransientTrackBuilder>  TTrackBuilder         ;
  edm::ESHandle<MagneticField>          magneticField         ;
  Cylinder::ConstCylinderPointer        ecalBarrel            ;
  Plane::ConstPlanePointer              ecalNegativeEtaEndcap ;
  Plane::ConstPlanePointer              ecalPositiveEtaEndcap ;
  //  MEtRecoilCorrection *                 metRecCorr;
  //  RecoilCorrector *                     corrector_ ;

  HLTConfigProvider HLTConfiguration;
  HLTPrescaleProvider* HLTPrescaleConfig;
  edm::Handle<edm::TriggerResults> HLTrigger;
  edm::Handle<edm::TriggerResults> MyHLTrigger;
  edm::Handle<trigger::TriggerEvent> HLTriggerEvent;
  //edm::Handle<l1extra::L1MuonParticleCollection> L1Muons;
  //edm::Handle<l1extra::L1EmParticleCollection> L1Electrons;
  //edm::Handle<l1extra::L1EmParticleCollection> L1ElectronsIso;
  vector<std::pair<string,string> > muontriggers;
  vector<std::pair<string,string> > electrontriggers;
  vector<std::pair<string,string> > tautriggers;
  vector<std::pair<string,string> > photontriggers;
  vector<std::pair<string,string> > jettriggers;

  vector<int> HLTriggerIndexSelection;
  vector<int> tauIndexSelection;
  vector<int> DiTauIndex;

  edm::Handle<edm::TriggerResults> Flags;
  math::XYZPoint pv_position;
  Vertex primvertex;

  //Data

  UInt_t errors;
  ULong64_t event_nr;
  UInt_t event_luminosityblock;
  UInt_t event_run;
  UInt_t event_timeunix;
  UInt_t event_timemicrosec;
  UChar_t trigger_level1bits[8];
  UChar_t trigger_level1[128];
  UChar_t trigger_HLT[128];
  Bool_t _passecalBadCalibFilterUpdate;
  // beam spot
  Float_t beamspot_x;
  Float_t beamspot_y;
  Float_t beamspot_z;
  Float_t beamspot_xwidth;
  Float_t beamspot_ywidth;
  Float_t beamspot_zsigma;
  Float_t beamspot_cov[6];

  // primary vertex
  UInt_t  primvertex_count;
  UInt_t  goodprimvertex_count;
  Float_t primvertex_x;
  Float_t primvertex_y;
  Float_t primvertex_z;
  Float_t primvertex_chi2;
  Float_t primvertex_ndof;
  Float_t primvertex_ptq;
  Int_t   primvertex_ntracks;
  Float_t primvertex_cov[6];
  Float_t primvertex_mindz;

  // primary vertex
  Float_t primvertexwithbs_x;
  Float_t primvertexwithbs_y;
  Float_t primvertexwithbs_z;
  Float_t primvertexwithbs_chi2;
  Float_t primvertexwithbs_ndof;
  Int_t   primvertexwithbs_ntracks;
  Float_t primvertexwithbs_cov[6];

  // re-fitted vertex
  UInt_t  refitvertex_count;
  Float_t refitvertex_x[M_refitvtxmaxcount];
  Float_t refitvertex_y[M_refitvtxmaxcount];
  Float_t refitvertex_z[M_refitvtxmaxcount];
  Float_t refitvertex_chi2[M_refitvtxmaxcount];
  Float_t refitvertex_ndof[M_refitvtxmaxcount];
  Int_t   refitvertex_ntracks[M_refitvtxmaxcount];
  Float_t refitvertex_cov[M_refitvtxmaxcount][6];
  Float_t refitvertex_mindz[M_refitvtxmaxcount];
  Int_t   refitvertex_eleIndex[M_refitvtxmaxcount][2];
  Int_t   refitvertex_muIndex[M_refitvtxmaxcount][2];
  Int_t   refitvertex_tauIndex[M_refitvtxmaxcount][2];

  // re-fitted vertex with bs
  UInt_t  refitvertexwithbs_count;
  Float_t refitvertexwithbs_x[M_refitvtxmaxcount];
  Float_t refitvertexwithbs_y[M_refitvtxmaxcount];
  Float_t refitvertexwithbs_z[M_refitvtxmaxcount];
  Float_t refitvertexwithbs_chi2[M_refitvtxmaxcount];
  Float_t refitvertexwithbs_ndof[M_refitvtxmaxcount];
  Int_t   refitvertexwithbs_ntracks[M_refitvtxmaxcount];
  Float_t refitvertexwithbs_cov[M_refitvtxmaxcount][6];
  Float_t refitvertexwithbs_mindz[M_refitvtxmaxcount];
  Int_t   refitvertexwithbs_eleIndex[M_refitvtxmaxcount][2];
  Int_t   refitvertexwithbs_muIndex[M_refitvtxmaxcount][2];
  Int_t   refitvertexwithbs_tauIndex[M_refitvtxmaxcount][2];

  // tracks
  UInt_t track_count;
  Float_t track_px[M_trackmaxcount];
  Float_t track_py[M_trackmaxcount];
  Float_t track_pz[M_trackmaxcount];
  Float_t track_pt[M_trackmaxcount];
  Float_t track_eta[M_trackmaxcount];
  Float_t track_phi[M_trackmaxcount];
  Float_t track_mass[M_trackmaxcount];
  Float_t track_charge[M_trackmaxcount];
  Float_t track_outerx[M_trackmaxcount];
  Float_t track_outery[M_trackmaxcount];
  Float_t track_outerz[M_trackmaxcount];
  Float_t track_closestpointx[M_trackmaxcount];
  Float_t track_closestpointy[M_trackmaxcount];
  Float_t track_closestpointz[M_trackmaxcount];
  Float_t track_chi2[M_trackmaxcount];
  Float_t track_ndof[M_trackmaxcount];
  Float_t track_dxy[M_trackmaxcount];
  Float_t track_dxyerr[M_trackmaxcount];
  Float_t track_dz[M_trackmaxcount];
  Float_t track_dzerr[M_trackmaxcount];
  Float_t track_dedxharmonic2[M_trackmaxcount];
  UChar_t track_nhits[M_trackmaxcount];
  UChar_t track_nmissinghits[M_trackmaxcount];
  UChar_t track_npixelhits[M_trackmaxcount];
  UChar_t track_npixellayers[M_trackmaxcount];
  UChar_t track_nstriplayers[M_trackmaxcount];
  Int_t track_ID[M_trackmaxcount];
  Bool_t track_highPurity[M_trackmaxcount];
  Float_t track_vx[M_trackmaxcount];
  Float_t track_vy[M_trackmaxcount];
  Float_t track_vz[M_trackmaxcount];


  // pat muons
  UInt_t muon_count;
  Float_t muon_helixparameters[M_muonmaxcount][5];
  Float_t muon_helixparameters_covar[M_muonmaxcount][5][5];
  Float_t muon_referencePoint[M_muonmaxcount][3];
  Float_t muon_Bfield[M_muonmaxcount];
  Float_t muon_px[M_muonmaxcount];
  Float_t muon_py[M_muonmaxcount];
  Float_t muon_pz[M_muonmaxcount];
  Float_t muon_pt[M_muonmaxcount];
  Float_t muon_eta[M_muonmaxcount];
  Float_t muon_phi[M_muonmaxcount];
  Float_t muon_pterror[M_muonmaxcount];
  Float_t muon_chi2[M_muonmaxcount];
  Float_t muon_normChi2[M_muonmaxcount];
  Float_t muon_ndof[M_muonmaxcount];
  Float_t muon_charge[M_muonmaxcount];
  Float_t muon_miniISO[M_muonmaxcount];
  // needed for medium muon Id
  Float_t muon_combQ_chi2LocalPosition[M_muonmaxcount];
  Float_t muon_combQ_trkKink[M_muonmaxcount];
  Float_t muon_validFraction[M_muonmaxcount];
  Float_t muon_segmentComp[M_muonmaxcount];

  UInt_t muon_nMuonStations[M_muonmaxcount];
  UInt_t muon_nMuonHits[M_muonmaxcount];
  UInt_t muon_nPixelHits[M_muonmaxcount];
  UInt_t muon_nTrackerHits[M_muonmaxcount];

  Float_t muon_dz[M_muonmaxcount];
  Float_t muon_dzerr[M_muonmaxcount];
  Float_t muon_dxy[M_muonmaxcount];
  Float_t muon_dxyerr[M_muonmaxcount];
  Float_t muon_vx[M_muonmaxcount];
  Float_t muon_vy[M_muonmaxcount];
  Float_t muon_vz[M_muonmaxcount];

  Float_t muon_chargedHadIso[M_muonmaxcount];
  Float_t muon_neutralHadIso[M_muonmaxcount];
  Float_t muon_photonIso[M_muonmaxcount];
  Float_t muon_puIso[M_muonmaxcount];

  Float_t muon_r03_sumChargedHadronPt[M_muonmaxcount];
  Float_t muon_r03_sumChargedParticlePt[M_muonmaxcount];
  Float_t muon_r03_sumNeutralHadronEt[M_muonmaxcount];
  Float_t muon_r03_sumPhotonEt[M_muonmaxcount];
  Float_t muon_r03_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
  Float_t muon_r03_sumPhotonEtHighThreshold[M_muonmaxcount];
  Float_t muon_r03_sumPUPt[M_muonmaxcount];

  Float_t muon_r04_sumChargedHadronPt[M_muonmaxcount];
  Float_t muon_r04_sumChargedParticlePt[M_muonmaxcount];
  Float_t muon_r04_sumNeutralHadronEt[M_muonmaxcount];
  Float_t muon_r04_sumPhotonEt[M_muonmaxcount];
  Float_t muon_r04_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
  Float_t muon_r04_sumPhotonEtHighThreshold[M_muonmaxcount];
  Float_t muon_r04_sumPUPt[M_muonmaxcount];

  Bool_t muon_isPF[M_muonmaxcount];
  Bool_t muon_isGlobal[M_muonmaxcount];
  Bool_t muon_isTracker[M_muonmaxcount];
  Bool_t muon_isTight[M_muonmaxcount];
  Bool_t muon_isLoose[M_muonmaxcount];
  Bool_t muon_isMedium[M_muonmaxcount];
  Bool_t muon_isICHEP[M_muonmaxcount];
  Bool_t muon_isDuplicate[M_muonmaxcount];
  Bool_t muon_isBad[M_muonmaxcount];

  Bool_t muon_globalTrack[M_muonmaxcount];
  Bool_t muon_innerTrack[M_muonmaxcount];

  Int_t muon_genmatch[M_muonmaxcount];

  /*UInt_t dimuon_count;
  UInt_t dimuon_leading[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  UInt_t dimuon_trailing[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist2D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist2DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist3D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist3DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  */

  // pat jets
  UInt_t pfjet_count;
  Float_t pfjet_e[M_jetmaxcount];
  Float_t pfjet_px[M_jetmaxcount];
  Float_t pfjet_py[M_jetmaxcount];
  Float_t pfjet_pz[M_jetmaxcount];
  Float_t pfjet_pt[M_jetmaxcount];
  Float_t pfjet_eta[M_jetmaxcount];
  Float_t pfjet_phi[M_jetmaxcount];

  Float_t pfjet_neutralhadronicenergy[M_jetmaxcount];
  Float_t pfjet_chargedhadronicenergy[M_jetmaxcount];
  Float_t pfjet_neutralemenergy[M_jetmaxcount];
  Float_t pfjet_chargedemenergy[M_jetmaxcount];
  Float_t pfjet_muonenergy[M_jetmaxcount];
  Float_t pfjet_chargedmuonenergy[M_jetmaxcount];

  UInt_t pfjet_chargedmulti[M_jetmaxcount];
  UInt_t pfjet_neutralmulti[M_jetmaxcount];
  UInt_t pfjet_chargedhadronmulti[M_jetmaxcount];
  Float_t pfjet_energycorr[M_jetmaxcount];
  Float_t pfjet_energycorr_l1fastjet[M_jetmaxcount];
  Float_t pfjet_energycorr_l2relative[M_jetmaxcount];
  Float_t pfjet_energycorr_l3absolute[M_jetmaxcount];
  Float_t pfjet_energycorr_l2l3residual[M_jetmaxcount];
  Int_t pfjet_flavour[M_jetmaxcount];
  Float_t pfjet_btag[M_jetmaxcount][10];
  Float_t pfjet_jecUncertainty[M_jetmaxcount];
  Bool_t pfjet_pu_jet_fullId_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_fullId_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_fullId_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_fullDisc_mva[M_jetmaxcount];

  // pat jets puppi
  UInt_t pfjetpuppi_count;
  Float_t pfjetpuppi_e[M_jetmaxcount];
  Float_t pfjetpuppi_px[M_jetmaxcount];
  Float_t pfjetpuppi_py[M_jetmaxcount];
  Float_t pfjetpuppi_pz[M_jetmaxcount];
  Float_t pfjetpuppi_pt[M_jetmaxcount];
  Float_t pfjetpuppi_eta[M_jetmaxcount];
  Float_t pfjetpuppi_phi[M_jetmaxcount];

  Float_t pfjetpuppi_neutralhadronicenergy[M_jetmaxcount];
  Float_t pfjetpuppi_chargedhadronicenergy[M_jetmaxcount];
  Float_t pfjetpuppi_neutralemenergy[M_jetmaxcount];
  Float_t pfjetpuppi_chargedemenergy[M_jetmaxcount];
  Float_t pfjetpuppi_muonenergy[M_jetmaxcount];
  Float_t pfjetpuppi_chargedmuonenergy[M_jetmaxcount];

  UInt_t pfjetpuppi_chargedmulti[M_jetmaxcount];
  UInt_t pfjetpuppi_neutralmulti[M_jetmaxcount];
  UInt_t pfjetpuppi_chargedhadronmulti[M_jetmaxcount];
  Float_t pfjetpuppi_energycorr[M_jetmaxcount];
  Float_t pfjetpuppi_energycorr_l1fastjet[M_jetmaxcount];
  Float_t pfjetpuppi_energycorr_l2relative[M_jetmaxcount];
  Float_t pfjetpuppi_energycorr_l3absolute[M_jetmaxcount];
  Float_t pfjetpuppi_energycorr_l2l3residual[M_jetmaxcount];
  Int_t pfjetpuppi_flavour[M_jetmaxcount];
  Float_t pfjetpuppi_btag[M_jetmaxcount][10];
  Float_t pfjetpuppi_jecUncertainty[M_jetmaxcount];

  // pat electrons
  UInt_t electron_count;
  Float_t electron_helixparameters[M_electronmaxcount][5];
  Float_t electron_helixparameters_covar[M_electronmaxcount][5][5];
  Float_t electron_referencePoint[M_electronmaxcount][3];
  Float_t electron_Bfield[M_electronmaxcount];
  Float_t electron_px[M_electronmaxcount];
  Float_t electron_py[M_electronmaxcount];
  Float_t electron_pz[M_electronmaxcount];
  Float_t electron_pt[M_electronmaxcount];
  Float_t electron_eta[M_electronmaxcount];
  Float_t electron_px_energyscale_down[M_electronmaxcount];
  Float_t electron_px_energyscale_up[M_electronmaxcount];
  Float_t electron_py_energyscale_down[M_electronmaxcount];
  Float_t electron_py_energyscale_up[M_electronmaxcount];
  Float_t electron_pz_energyscale_down[M_electronmaxcount];
  Float_t electron_pz_energyscale_up[M_electronmaxcount];
  Float_t electron_pt_energyscale_down[M_electronmaxcount];
  Float_t electron_pt_energyscale_up[M_electronmaxcount];
  Float_t electron_px_energysigma_down[M_electronmaxcount];
  Float_t electron_px_energysigma_up[M_electronmaxcount];
  Float_t electron_py_energysigma_down[M_electronmaxcount];
  Float_t electron_py_energysigma_up[M_electronmaxcount];
  Float_t electron_pz_energysigma_down[M_electronmaxcount];
  Float_t electron_pz_energysigma_up[M_electronmaxcount];
  Float_t electron_pt_energysigma_down[M_electronmaxcount];
  Float_t electron_pt_energysigma_up[M_electronmaxcount];
  Float_t electron_phi[M_electronmaxcount];
  Float_t electron_trackchi2[M_electronmaxcount];
  Float_t electron_trackndof[M_electronmaxcount];
  Float_t electron_outerx[M_electronmaxcount];
  Float_t electron_outery[M_electronmaxcount];
  Float_t electron_outerz[M_electronmaxcount];
  Float_t electron_vx[M_electronmaxcount];
  Float_t electron_vy[M_electronmaxcount];
  Float_t electron_vz[M_electronmaxcount];
  Float_t electron_esuperclusterovertrack[M_electronmaxcount];
  Float_t electron_eseedclusterovertrack[M_electronmaxcount];
  Float_t electron_deltaetasuperclustertrack[M_electronmaxcount];
  Float_t electron_deltaphisuperclustertrack[M_electronmaxcount];
  Float_t electron_e1x5[M_electronmaxcount];
  Float_t electron_e2x5[M_electronmaxcount];
  Float_t electron_e5x5[M_electronmaxcount];
  Float_t electron_sigmaetaeta[M_electronmaxcount];
  Float_t electron_sigmaietaieta[M_electronmaxcount];
  Float_t electron_ehcaloverecal[M_electronmaxcount];
  Float_t electron_ehcaloverecaldepth1[M_electronmaxcount];
  Float_t electron_ehcaloverecaldepth2[M_electronmaxcount];
  Float_t electron_full5x5_sigmaietaieta[M_electronmaxcount];
  Float_t electron_ooemoop[M_electronmaxcount];

  Float_t electron_detaInSeed[M_electronmaxcount];
  Float_t electron_he[M_electronmaxcount];
  Float_t electron_eaIsolation[M_electronmaxcount];

  Float_t electron_superClusterEta[M_electronmaxcount];
  Float_t electron_superClusterPhi[M_electronmaxcount];

  Float_t electron_superClusterX[M_electronmaxcount];
  Float_t electron_superClusterY[M_electronmaxcount];
  Float_t electron_superClusterZ[M_electronmaxcount];

  Float_t electron_chargedHadIso[M_electronmaxcount];
  Float_t electron_neutralHadIso[M_electronmaxcount];
  Float_t electron_photonIso[M_electronmaxcount];
  Float_t electron_puIso[M_electronmaxcount];
  Float_t electron_miniISO[M_electronmaxcount];

  Float_t electron_r03_sumChargedHadronPt[M_electronmaxcount];
  Float_t electron_r03_sumChargedParticlePt[M_electronmaxcount];
  Float_t electron_r03_sumNeutralHadronEt[M_electronmaxcount];
  Float_t electron_r03_sumPhotonEt[M_electronmaxcount];
  Float_t electron_r03_sumNeutralHadronEtHighThreshold[M_electronmaxcount];
  Float_t electron_r03_sumPhotonEtHighThreshold[M_electronmaxcount];
  Float_t electron_r03_sumPUPt[M_electronmaxcount];

  Float_t electron_charge[M_electronmaxcount];
  UChar_t electron_nhits[M_electronmaxcount];
  UChar_t electron_nmissinghits[M_electronmaxcount];
  UChar_t electron_npixelhits[M_electronmaxcount];
  UChar_t electron_nmissinginnerhits[M_electronmaxcount];
  UChar_t electron_npixellayers[M_electronmaxcount];
  UChar_t electron_nstriplayers[M_electronmaxcount];
  Float_t electron_dxy[M_electronmaxcount];
  Float_t electron_dxyerr[M_electronmaxcount];
  Float_t electron_dz[M_electronmaxcount];
  Float_t electron_dzerr[M_electronmaxcount];
  Float_t electron_convdist[M_electronmaxcount];
  Float_t electron_convdcot[M_electronmaxcount];
  Float_t electron_convradius[M_electronmaxcount];
  UInt_t electron_gapinfo[M_electronmaxcount];
  UInt_t electron_chargeinfo[M_electronmaxcount];
  Float_t electron_fbrems[M_electronmaxcount];
  Int_t electron_numbrems[M_electronmaxcount];
  Int_t electron_superclusterindex[M_electronmaxcount];
  UChar_t electron_info[M_electronmaxcount];

  Float_t electron_mva_value_Spring16_v1[M_electronmaxcount];
  Float_t electron_mva_wp80_general_Spring16_v1[M_electronmaxcount];
  Float_t electron_mva_wp90_general_Spring16_v1[M_electronmaxcount];

    //new for 9.4.0 Fall17
  Float_t electron_mva_value_Iso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_value_noIso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_wp90_Iso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_wp80_Iso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_Loose_Iso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_wp90_noIso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_wp80_noIso_Fall17_v1[M_electronmaxcount];
  Float_t electron_mva_Loose_noIso_Fall17_v1[M_electronmaxcount];

  Float_t electron_mva_value_Iso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_value_noIso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_wp90_Iso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_wp80_Iso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_Loose_Iso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_wp90_noIso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_wp80_noIso_Fall17_v2[M_electronmaxcount];
  Float_t electron_mva_Loose_noIso_Fall17_v2[M_electronmaxcount];

  Bool_t electron_cutId_veto_Summer16[M_electronmaxcount];
  //Bool_t electron_cutId_loose_Summer16[M_electronmaxcount];
  //Bool_t electron_cutId_medium_Summer16[M_electronmaxcount];
  //Bool_t electron_cutId_tight_Summer16[M_electronmaxcount];

  Bool_t electron_cutId_veto_Fall17[M_electronmaxcount];
  //Bool_t electron_cutId_loose_Fall17[M_electronmaxcount];
  //Bool_t electron_cutId_medium_Fall17[M_electronmaxcount];
  //Bool_t electron_cutId_tight_Fall17[M_electronmaxcount];

  Bool_t electron_cutId_veto_Fall17V2[M_electronmaxcount];
  //Bool_t electron_cutId_loose_Fall17V2[M_electronmaxcount];
  //Bool_t electron_cutId_medium_Fall17V2[M_electronmaxcount];
  //Bool_t electron_cutId_tight_Fall17V2[M_electronmaxcount];

  Bool_t electron_pass_conversion[M_electronmaxcount];

  Int_t electron_genmatch[M_electronmaxcount];

  UInt_t photon_count;
  Float_t photon_px[M_photonmaxcount];
  Float_t photon_py[M_photonmaxcount];
  Float_t photon_pz[M_photonmaxcount];
  Float_t photon_e1x5[M_photonmaxcount];
  Float_t photon_e2x5[M_photonmaxcount];
  Float_t photon_e3x3[M_photonmaxcount];
  Float_t photon_e5x5[M_photonmaxcount];
  Float_t photon_sigmaetaeta[M_photonmaxcount];
  Float_t photon_sigmaietaieta[M_photonmaxcount];
  Float_t photon_ehcaloverecal[M_photonmaxcount];
  Float_t photon_ehcaloverecaldepth1[M_photonmaxcount];
  Float_t photon_ehcaloverecaldepth2[M_photonmaxcount];
  Float_t photon_maxenergyxtal[M_photonmaxcount];
  Float_t photon_isolationr3track[M_photonmaxcount];
  Float_t photon_isolationr3trackhollow[M_photonmaxcount];
  UInt_t photon_isolationr3ntrack[M_photonmaxcount];
  UInt_t photon_isolationr3ntrackhollow[M_photonmaxcount];
  Float_t photon_isolationr3ecal[M_photonmaxcount];
  Float_t photon_isolationr3hcal[M_photonmaxcount];
  Float_t photon_isolationr4track[M_photonmaxcount];
  Float_t photon_isolationr4trackhollow[M_photonmaxcount];
  UInt_t photon_isolationr4ntrack[M_photonmaxcount];
  UInt_t photon_isolationr4ntrackhollow[M_photonmaxcount];
  Float_t photon_isolationr4ecal[M_photonmaxcount];
  Float_t photon_isolationr4hcal[M_photonmaxcount];
  Int_t photon_superclusterindex[M_photonmaxcount];
  UChar_t photon_info[M_photonmaxcount];
  UInt_t photon_gapinfo[M_photonmaxcount];
  UInt_t photon_conversionbegin[M_photonmaxcount];

  // taus
  UInt_t tau_count;
  Float_t tau_helixparameters[M_taumaxcount][5];
  Float_t tau_helixparameters_covar[M_taumaxcount][5][5];
  Float_t tau_referencePoint[M_taumaxcount][3];
  Float_t tau_Bfield[M_muonmaxcount];

  Float_t tau_e[M_taumaxcount];
  Float_t tau_px[M_taumaxcount];
  Float_t tau_py[M_taumaxcount];
  Float_t tau_pz[M_taumaxcount];
  Float_t tau_pt[M_taumaxcount];
  Float_t tau_eta[M_taumaxcount];
  Float_t tau_phi[M_taumaxcount];
  Float_t tau_mass[M_taumaxcount];

  Float_t tau_leadchargedhadrcand_px[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_py[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_pz[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_mass[M_taumaxcount];
  Int_t   tau_leadchargedhadrcand_id[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_dxy[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_dz[M_taumaxcount];
  Int_t   tau_leadchargedhadrcand_lostPixelHits[M_taumaxcount];
  Int_t   tau_leadchargedhadrcand_pvAssocQ[M_taumaxcount];

  Float_t tau_vertexx[M_taumaxcount];
  Float_t tau_vertexy[M_taumaxcount];
  Float_t tau_vertexz[M_taumaxcount];

  Float_t tau_pca2D_x[M_taumaxcount];
  Float_t tau_pca2D_y[M_taumaxcount];
  Float_t tau_pca2D_z[M_taumaxcount];

  Float_t tau_pca3D_x[M_taumaxcount];
  Float_t tau_pca3D_y[M_taumaxcount];
  Float_t tau_pca3D_z[M_taumaxcount];

  Float_t tau_dxy[M_taumaxcount];
  Float_t tau_dxySig[M_taumaxcount];
  Float_t tau_dz[M_taumaxcount];
  Float_t tau_ip3d[M_taumaxcount];
  Float_t tau_ip3dSig[M_taumaxcount];
  Float_t tau_flightLength[M_taumaxcount];
  Float_t tau_flightLengthSig[M_taumaxcount];
  Float_t tau_SV_x[M_taumaxcount];
  Float_t tau_SV_y[M_taumaxcount];
  Float_t tau_SV_z[M_taumaxcount];
  Float_t tau_SV_cov[M_taumaxcount][6];
  Float_t tau_SV_cov2[M_taumaxcount][3];

  Float_t tau_charge[M_taumaxcount];
  Float_t tau_genjet_e[M_taumaxcount];
  Float_t tau_genjet_px[M_taumaxcount];
  Float_t tau_genjet_py[M_taumaxcount];
  Float_t tau_genjet_pz[M_taumaxcount];
  Int_t tau_genmatch[M_taumaxcount];

  // main tau discriminators
  bool setTauBranches;
  std::vector<std::pair<std::string, unsigned int> >tauIdIndx;
  Float_t tau_ids[200][M_taumaxcount];

  // l1 match
  Bool_t  tau_L1trigger_match[M_taumaxcount];

  UInt_t tau_signalChargedHadrCands_size[M_taumaxcount];
  UInt_t tau_signalNeutralHadrCands_size[M_taumaxcount];
  UInt_t tau_signalGammaCands_size[M_taumaxcount];
  UInt_t tau_isolationChargedHadrCands_size[M_taumaxcount];
  UInt_t tau_isolationNeutralHadrCands_size[M_taumaxcount];
  UInt_t tau_isolationGammaCands_size[M_taumaxcount];

  string tau_genDecayMode_name[M_taumaxcount];
  Int_t tau_genDecayMode[M_taumaxcount];

  string tau_decayMode_name[M_taumaxcount];
  Int_t tau_decayMode[M_taumaxcount];

  // tau consituents
  UInt_t tau_constituents_count[M_taumaxcount];
  Float_t tau_constituents_px[M_taumaxcount][50];
  Float_t tau_constituents_py[M_taumaxcount][50];
  Float_t tau_constituents_pz[M_taumaxcount][50];
  Float_t tau_constituents_e[M_taumaxcount][50];
  Float_t tau_constituents_mass[M_taumaxcount][50];
  Int_t tau_constituents_charge[M_taumaxcount][50];
  Float_t tau_constituents_vx[M_taumaxcount][50];
  Float_t tau_constituents_vy[M_taumaxcount][50];
  Float_t tau_constituents_vz[M_taumaxcount][50];
  Int_t tau_constituents_pdgId[M_taumaxcount][50];
  Int_t tau_constituents_lostPixelHits[M_taumaxcount][50];

  // generated tau
  UInt_t gentau_count;
  Float_t gentau_px[M_taumaxcount];
  Float_t gentau_py[M_taumaxcount];
  Float_t gentau_pz[M_taumaxcount];
  Float_t gentau_e[M_taumaxcount];
  Float_t gentau_charge[M_taumaxcount];

  Float_t gentau_visible_px[M_taumaxcount];
  Float_t gentau_visible_py[M_taumaxcount];
  Float_t gentau_visible_pz[M_taumaxcount];
  Float_t gentau_visible_e[M_taumaxcount];

  Float_t gentau_visible_pt[M_taumaxcount];
  Float_t gentau_visible_eta[M_taumaxcount];
  Float_t gentau_visible_phi[M_taumaxcount];
  Float_t gentau_visible_mass[M_taumaxcount];
  Float_t gentau_visibleNoLep_e[M_taumaxcount];

  Float_t gentau_visibleNoLep_pt[M_taumaxcount];
  Float_t gentau_visibleNoLep_eta[M_taumaxcount];
  Float_t gentau_visibleNoLep_phi[M_taumaxcount];
  Float_t gentau_visibleNoLep_mass[M_taumaxcount];

  Int_t gentau_status[M_taumaxcount];
  Int_t gentau_fromHardProcess[M_taumaxcount];
  Int_t gentau_fromHardProcessBeforeFSR[M_taumaxcount];
  Int_t gentau_isDecayedLeptonHadron[M_taumaxcount];
  Int_t gentau_isDirectHadronDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectHardProcessTauDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectPromptTauDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectTauDecayProduct[M_taumaxcount];
  Int_t gentau_isFirstCopy[M_taumaxcount];
  Int_t gentau_isHardProcess[M_taumaxcount];
  Int_t gentau_isHardProcessTauDecayProduct[M_taumaxcount];
  Int_t gentau_isLastCopy[M_taumaxcount];
  Int_t gentau_isLastCopyBeforeFSR[M_taumaxcount];
  Int_t gentau_isPrompt[M_taumaxcount];
  Int_t gentau_isPromptTauDecayProduct[M_taumaxcount];
  Int_t gentau_isTauDecayProduct[M_taumaxcount];

  Int_t gentau_decayMode[M_taumaxcount];
  string gentau_decayMode_name[M_taumaxcount];
  UChar_t gentau_mother[M_taumaxcount];

  // L1 Objects
  UInt_t  l1muon_count;
  Float_t l1muon_px[M_muonmaxcount];
  Float_t l1muon_py[M_muonmaxcount];
  Float_t l1muon_pz[M_muonmaxcount];
  Float_t l1muon_pt[M_muonmaxcount];
//  Int_t   l1muon_ipt[M_muonmaxcount];
  Int_t   l1muon_eta[M_muonmaxcount];
  Int_t   l1muon_phi[M_muonmaxcount];
  Int_t   l1muon_iso[M_muonmaxcount];
//  Int_t   l1muon_qual[M_muonmaxcount];
  Int_t   l1muon_charge[M_muonmaxcount];
//  Int_t   l1muon_chargeValid[M_muonmaxcount];
  Int_t   l1muon_muonIndex[M_muonmaxcount];
//  Int_t   l1muon_tag[M_muonmaxcount];
//  Int_t   l1muon_isoSum[M_muonmaxcount];
//  Int_t   l1muon_dPhiExtra[M_muonmaxcount];
//  Int_t   l1muon_dEtaExtra[M_muonmaxcount];
//  Int_t   l1muon_rank[M_muonmaxcount];
//  Int_t   l1muon_bx[M_muonmaxcount];

  UInt_t  l1egamma_count;
  Float_t l1egamma_px[M_electronmaxcount];
  Float_t l1egamma_py[M_electronmaxcount];
  Float_t l1egamma_pz[M_electronmaxcount];
  Float_t l1egamma_pt[M_electronmaxcount];
//  Int_t   l1egamma_ipt[M_electronmaxcount];
  Int_t   l1egamma_eta[M_electronmaxcount];
  Int_t   l1egamma_phi[M_electronmaxcount];
  Int_t   l1egamma_iso[M_electronmaxcount];
//  Int_t   l1egamma_qual[M_electronmaxcount];
//  Int_t   l1egamma_towerIEta[M_electronmaxcount];
//  Int_t   l1egamma_towerIPhi[M_electronmaxcount];
//  Int_t   l1egamma_rawEt[M_electronmaxcount];
//  Int_t   l1egamma_isoEt[M_electronmaxcount];
//  Int_t   l1egamma_footprintEt[M_electronmaxcount];
//  Int_t   l1egamma_nTT[M_electronmaxcount];
//  Int_t   l1egamma_shape[M_electronmaxcount];
//  Int_t   l1egamma_bx[M_electronmaxcount];

  UInt_t  l1tau_count;
  Float_t l1tau_px[M_taumaxcount];
  Float_t l1tau_py[M_taumaxcount];
  Float_t l1tau_pz[M_taumaxcount];
  Float_t l1tau_pt[M_taumaxcount];
  Int_t   l1tau_ipt[M_taumaxcount];
  Int_t l1tau_eta[M_taumaxcount];
  Int_t l1tau_phi[M_taumaxcount];
  Int_t   l1tau_iso[M_taumaxcount];
  Int_t   l1tau_qual[M_taumaxcount];
  Int_t   l1tau_towerIEta[M_taumaxcount];
  Int_t   l1tau_towerIPhi[M_taumaxcount];
  Int_t   l1tau_rawEt[M_taumaxcount];
  Int_t   l1tau_isoEt[M_taumaxcount];
  Int_t   l1tau_nTT[M_taumaxcount];
  Int_t   l1tau_hasEM[M_taumaxcount];
  Int_t   l1tau_isMerged[M_taumaxcount];
  Int_t l1tau_bx[M_taumaxcount];

  UInt_t l1isotau_count;
  Float_t l1isotau_e[M_taumaxcount];
  Float_t l1isotau_px[M_taumaxcount];
  Float_t l1isotau_py[M_taumaxcount];
  Float_t l1isotau_pz[M_taumaxcount];
  Float_t l1isotau_pt[M_taumaxcount];
  Float_t l1isotau_eta[M_taumaxcount];
  Float_t l1isotau_phi[M_taumaxcount];
  Float_t l1isotau_mass[M_taumaxcount];
  Float_t l1isotau_charge[M_taumaxcount];
  Int_t   l1isotau_iso[M_taumaxcount];


  // rho neutral
  Float_t rhoNeutral;

  // met
  Float_t pfmet_ex;
  Float_t pfmet_ey;
  Float_t pfmet_ez;
  Float_t pfmet_pt;
  Float_t pfmet_phi;
  Float_t pfmet_sumet;

  Float_t pfmet_sig;
  Float_t pfmet_sigxx;
  Float_t pfmet_sigxy;
  Float_t pfmet_sigyx;
  Float_t pfmet_sigyy;

  Float_t pfmet_ex_JetEnUp;
  Float_t pfmet_ey_JetEnUp;

  Float_t pfmet_ex_JetEnDown;
  Float_t pfmet_ey_JetEnDown;

  Float_t pfmet_ex_UnclusteredEnUp;
  Float_t pfmet_ey_UnclusteredEnUp;

  Float_t pfmet_ex_UnclusteredEnDown;
  Float_t pfmet_ey_UnclusteredEnDown;

  Float_t pfmetcorr_ex;
  Float_t pfmetcorr_ey;
  Float_t pfmetcorr_ez;
  Float_t pfmetcorr_pt;
  Float_t pfmetcorr_phi;
  Float_t pfmetcorr_sumet;

  Float_t pfmetcorr_sig;
  Float_t pfmetcorr_sigxx;
  Float_t pfmetcorr_sigxy;
  Float_t pfmetcorr_sigyx;
  Float_t pfmetcorr_sigyy;

  Float_t pfmetcorr_ex_JetEnUp;
  Float_t pfmetcorr_ey_JetEnUp;

  Float_t pfmetcorr_ex_JetEnDown;
  Float_t pfmetcorr_ey_JetEnDown;

  Float_t pfmetcorr_ex_UnclusteredEnUp;
  Float_t pfmetcorr_ey_UnclusteredEnUp;

  Float_t pfmetcorr_ex_UnclusteredEnDown;
  Float_t pfmetcorr_ey_UnclusteredEnDown;

  Float_t      pfmetcorr_ex_JetResDown;
  Float_t      pfmetcorr_ey_JetResDown;
  Float_t      pfmetcorr_ex_JetResUp;
  Float_t      pfmetcorr_ey_JetResUp;

/*  Float_t      pfmetcorr_ex_smeared;
  Float_t      pfmetcorr_ey_smeared;
  Float_t      pfmetcorr_ez_smeared;
  Float_t      pfmetcorr_pt_smeared;
  Float_t      pfmetcorr_phi_smeared;


   Float_t     pfmetcorr_ex_JetEnUp_smeared;
   Float_t     pfmetcorr_ey_JetEnUp_smeared;

  Float_t      pfmetcorr_ex_JetEnDown_smeared;
  Float_t      pfmetcorr_ey_JetEnDown_smeared;

  Float_t      pfmetcorr_ex_UnclusteredEnUp_smeared;
  Float_t      pfmetcorr_ey_UnclusteredEnUp_smeared;

  Float_t      pfmetcorr_ex_UnclusteredEnDown_smeared;
  Float_t      pfmetcorr_ey_UnclusteredEnDown_smeared;

  Float_t      pfmetcorr_ex_JetResUp_smeared;
  Float_t      pfmetcorr_ey_JetResUp_smeared;

  Float_t      pfmetcorr_ex_JetResDown_smeared;
  Float_t      pfmetcorr_ey_JetResDown_smeared; */

  Float_t puppimet_ex;
  Float_t puppimet_ey;
  Float_t puppimet_ez;
  Float_t puppimet_pt;
  Float_t puppimet_phi;
  Float_t puppimet_sumet;

  Float_t puppimet_ex_JetEnUp;
  Float_t puppimet_ey_JetEnUp;

  Float_t puppimet_ex_JetEnDown;
  Float_t puppimet_ey_JetEnDown;

  Float_t puppimet_ex_UnclusteredEnUp;
  Float_t puppimet_ey_UnclusteredEnUp;

  Float_t puppimet_ex_UnclusteredEnDown;
  Float_t puppimet_ey_UnclusteredEnDown;

  Float_t puppimet_ex_JetResUp;
  Float_t puppimet_ey_JetResUp;

  Float_t puppimet_ex_JetResDown;
  Float_t puppimet_ey_JetResDown;

  Float_t puppimet_sigxx;
  Float_t puppimet_sigxy;
  Float_t puppimet_sigyx;
  Float_t puppimet_sigyy;

  UInt_t mvamet_count;
  Float_t mvamet_ex[M_mvametmaxcount];
  Float_t mvamet_ey[M_mvametmaxcount];

  Float_t mvamet_sigxx[M_mvametmaxcount];
  Float_t mvamet_sigxy[M_mvametmaxcount];
  Float_t mvamet_sigyx[M_mvametmaxcount];
  Float_t mvamet_sigyy[M_mvametmaxcount];
  UChar_t mvamet_channel[M_mvametmaxcount];
  UInt_t mvamet_lep1[M_mvametmaxcount];
  UInt_t mvamet_lep2[M_mvametmaxcount];
  Float_t mvamet_lep1_pt[M_mvametmaxcount];
  Float_t mvamet_lep2_pt[M_mvametmaxcount];

  Float_t genmet_ex;
  Float_t genmet_ey;

  //Generator Information
  Int_t htxs_stage0cat;
  Int_t htxs_stage1p1cat;
  Int_t htxs_stage1p1finecat;
  Float_t htxs_higgsPt;
  Int_t htxs_njets30;

  Float_t genweight;
  Float_t genid1;
  Float_t genx1;
  Float_t genid2;
  Float_t genx2;
  Float_t genScale;
  Float_t gen_pythiaweights[M_Nweights];

  Float_t weightScale[9];
  Float_t weightPDFmax;
  Float_t weightPDFmin;
  Float_t weightPDFmean;
  Float_t weightPDFup;
  Float_t weightPDFdown;
  Float_t weightPDFvar;

  Int_t numpileupinteractionsminus;
  Int_t numpileupinteractions;
  Int_t numpileupinteractionsplus;
  Float_t numtruepileupinteractions;
  Int_t hepNUP_;

  Float_t genparticles_lheHt;
  Float_t genparticles_lheWPt;
  UInt_t genparticles_noutgoing;
  Int_t genparticles_noutgoing_NLO;
  UInt_t genparticles_count;
  Float_t genparticles_e[M_genparticlesmaxcount];
  Float_t genparticles_px[M_genparticlesmaxcount];
  Float_t genparticles_py[M_genparticlesmaxcount];
  Float_t genparticles_pz[M_genparticlesmaxcount];
  Float_t genparticles_vx[M_genparticlesmaxcount];
  Float_t genparticles_vy[M_genparticlesmaxcount];
  Float_t genparticles_vz[M_genparticlesmaxcount];
  Int_t genparticles_pdgid[M_genparticlesmaxcount];
  Int_t genparticles_status[M_genparticlesmaxcount];
  UInt_t genparticles_info[M_genparticlesmaxcount];
  UChar_t genparticles_mother[M_genparticlesmaxcount];

  Int_t genparticles_fromHardProcess[M_genparticlesmaxcount];
  Int_t genparticles_fromHardProcessBeforeFSR[M_genparticlesmaxcount];
  Int_t genparticles_isDecayedLeptonHadron[M_genparticlesmaxcount];
  Int_t genparticles_isDirectHadronDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectHardProcessTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectPromptTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isFirstCopy[M_genparticlesmaxcount];
  Int_t genparticles_isHardProcess[M_genparticlesmaxcount];
  Int_t genparticles_isHardProcessTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isLastCopy[M_genparticlesmaxcount];
  Int_t genparticles_isLastCopyBeforeFSR[M_genparticlesmaxcount];
  Int_t genparticles_isPrompt[M_genparticlesmaxcount];
  Int_t genparticles_isPromptTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isTauDecayProduct[M_genparticlesmaxcount];

  // GenJets
  UInt_t genjets_count;
  Float_t genjets_e[M_genjetsmaxcount];
  Float_t genjets_px[M_genjetsmaxcount];
  Float_t genjets_py[M_genjetsmaxcount];
  Float_t genjets_pz[M_genjetsmaxcount];
  Float_t genjets_pt[M_genjetsmaxcount];
  Float_t genjets_eta[M_genjetsmaxcount];
  Float_t genjets_phi[M_genjetsmaxcount];
  Int_t genjets_pdgid[M_genjetsmaxcount];
  Int_t genjets_status[M_genjetsmaxcount];

  Float_t genjets_em_energy[M_genjetsmaxcount];
  Float_t genjets_had_energy[M_genjetsmaxcount];
  Float_t genjets_invisible_energy[M_genjetsmaxcount];
  Float_t genjets_auxiliary_energy[M_genjetsmaxcount];

  // trigger objects
  UInt_t trigobject_count;
  Float_t trigobject_px[M_trigobjectmaxcount];
  Float_t trigobject_py[M_trigobjectmaxcount];
  Float_t trigobject_pz[M_trigobjectmaxcount];
  Float_t trigobject_pt[M_trigobjectmaxcount];
  Float_t trigobject_eta[M_trigobjectmaxcount];
  Float_t trigobject_phi[M_trigobjectmaxcount];
  Bool_t  trigobject_filters[M_trigobjectmaxcount][M_hltfiltersmax];
  Bool_t trigobject_isMuon[M_trigobjectmaxcount];
  Bool_t trigobject_isElectron[M_trigobjectmaxcount];
  Bool_t trigobject_isTau[M_trigobjectmaxcount];
  Bool_t trigobject_isJet[M_trigobjectmaxcount];
  Bool_t trigobject_isMET[M_trigobjectmaxcount];

  //lumitree
  UInt_t lumi_run;
  UInt_t lumi_block;
  Float_t lumi_value;
  Float_t lumi_valueerr;
  Float_t lumi_livefrac;
  Float_t lumi_deadfrac;
  UInt_t lumi_quality;
  UInt_t lumi_eventsprocessed;
  UInt_t lumi_eventsfiltered;
  UInt_t lumi_hltprescaletable;
  UInt_t lumi_l1algoprescaletable;
  UInt_t lumi_l1techprescaletable;

  //runtree
  UInt_t run_number;
  UInt_t run_hltcount;
  vector<string> run_hltnames;
  vector<string> run_hltfilters;
  vector<string> run_hltmufilters;
  vector<string> run_hltelectronfilters;
  vector<string> run_hlttaufilters;
  vector<string> run_hltphotonfilters;
  vector<string> run_hltjetfilters;
  vector<string> run_floattaudiscriminators;
  vector<string> run_binarytaudiscriminators;
  vector<string> run_btagdiscriminators;
  UInt_t run_hltprescaletablescount;
  UInt_t run_hltprescaletables[10000];
  UInt_t run_l1algocount;
  UInt_t run_l1algoprescaletablescount;
  UInt_t run_l1algoprescaletables[10000];
  UInt_t run_l1techcount;
  UInt_t run_l1techprescaletablescount;
  UInt_t run_l1techprescaletables[10000];

  // susy info
  Float_t SusyMotherMass;
  Float_t SusyLSPMass;

  std::map<std::string, int>* hltriggerresults_;
  std::map<std::string, int>* hltriggerprescales_;
  std::vector<std::string>hltriggerresultsV_;
  std::map<std::string, int>* flags_;
  float embeddingWeight_;

  Float_t prefiringweight;
  Float_t prefiringweightup;
  Float_t prefiringweightdown;
  //std::vector< double > embeddingWeights_; //for RhEmb
  Double_t TauSpinnerWeight[M_tauspinneranglesmaxcount];
  UInt_t TauSpinAngles_count;
  JetCorrectionUncertainty *jecUnc;

};

DEFINE_FWK_MODULE(NTupleMaker);

#endif
