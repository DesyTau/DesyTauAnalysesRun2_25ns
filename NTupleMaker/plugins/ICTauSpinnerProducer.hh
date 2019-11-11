//Merijn: I think that this is the right way to do it, not sure..
#ifndef DesyTauAnalyses_ICHiggsTauTau_ICTauSpinnerProducer_h
#define DesyTauAnalyses_ICHiggsTauTau_ICTauSpinnerProducer_h  //DesyTauAnalyses/NTupleMaker/plugins
#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

//Merijn add:
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
//#include "DesyTauAnalyses/NTupleMaker/plugins/NTupleMaker.h" currently will generate an issue with multiple definitions 
#define NrAnglestoStore   5

enum pdgId {
  Gamma = 22,
  PiZero = 111,
  PiPlus = 211,
  KPlus = 321,
  KLong = 130,
  KShort = 310,
  KZero=311,
  Electron = 11,
  NuE = 12,
  Muon = 13,
  NuMu = 14,
  NuTau = 16
};

/**
 * @brief See documentation [here](\ref objs-candidate)
 */
class ICTauSpinnerProducer : public edm::EDProducer {
 public:
  explicit ICTauSpinnerProducer(const edm::ParameterSet &);
  ~ICTauSpinnerProducer();

 private:
  virtual void beginJob();
  virtual void initialize();
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  reco::GenParticle getBoson(edm::Handle<edm::View<reco::GenParticle> > parts_handle);
  reco::GenParticle getLastCopy(reco::GenParticle part, edm::Handle<edm::View<reco::GenParticle> > parts_handle);
  std::vector<reco::GenParticle> getTaus(reco::GenParticle boson, edm::Handle<edm::View<reco::GenParticle> > parts_handle);
  void getTauDaughters(std::vector<reco::GenParticle> &tau_daughters, unsigned &type, reco::GenParticle tau, edm::Handle<edm::View<reco::GenParticle> > parts_handle);
  void removeGammas(std::vector<TauSpinner::SimpleParticle> &tau_daughters);
  void removeSecondGamma(std::vector<TauSpinner::SimpleParticle> &tau_daughters);
  TauSpinner::SimpleParticle ConvertToSimplePart(reco::GenParticle input_part);
  bool channelMatch(std::vector<reco::GenParticle> parts, std::vector<int> matches, bool signed_);
  std::vector<std::pair<std::string,double>> SplitString(std::string instring);
 // ic::EventInfo *info_;
  std::vector<std::pair<std::string,double>> theta_vec_;

  edm::InputTag input_;
  std::string branch_;
  std::string theta_;
  int bosonPdgId_;
  
  std::string TauSpinnerSettingsPDF;
  bool Ipp;
  int Ipol;
  int nonSM2;
  int nonSMN;
  double CMSENE;


  //Merijn add:
  TTree* tree;
  TTree* treeAngles;//this we'll fill only once!
  TH1D*  nEvents;//for test purpose..
  int NThetaAngles;
  double WeightsPtr[NrAnglestoStore];//in case we like to store 10 angles


};

#endif
