#ifndef TauAnalysis_CandidateTools_candidateAuxFunctions_h
#define TauAnalysis_CandidateTools_candidateAuxFunctions_h

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include <TVector2.h>

#include <vector>

/// Find the best match to the input four vector out of the GenParticleCollection
/// If pdgIds is non null (and pdgIdStrict is true), look for specific pdgIds.  If
/// pdgIdStrict is false, take the best candidate, but always prefer those which match
/// the pdgId list
const reco::GenParticle* findGenParticle(const reco::Candidate::LorentzVector& toMatch,
					 const reco::GenParticleCollection& genParticles,
                                         double dRMax = 0.5, int status = -1,
					 const std::vector<int>* pdgIds = 0, bool pdgIdStrict = true);

void findDaughters(const reco::GenParticle*, std::vector<const reco::GenParticle*>&, int = -1);

/// Find the effective secondary vertex of a generator level particle
reco::Candidate::Point getDecayVertex(const reco::GenParticle*);

bool isNeutrino(const reco::GenParticle*);

reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>&, int = 1);
reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle*);
reco::Candidate::LorentzVector getInvisMomentum(const std::vector<const reco::GenParticle*>&, int = 1);
reco::Candidate::LorentzVector getInvisMomentum(const reco::GenParticle*);

void compX1X2byCollinearApprox(double&, double&, double, double, double, double, double, double);
double getPhysX(double x, bool& isWithinPhysRange);

reco::Candidate::LorentzVector boostToRestFrame(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);

std::string getGenTauDecayMode(const reco::GenParticle*);

TVector2 getDiTauBisectorDirection(const reco::Candidate::LorentzVector& leg1P4, const reco::Candidate::LorentzVector& leg2P4);

std::string getTauDecayModeName(int);

const reco::Candidate* getDistPion(const pat::Tau&);
const reco::Candidate* getDistPion(const reco::GenJet&);

std::pair<double, double> compMEtProjU(const reco::Candidate::LorentzVector&, double, double, int&, bool = true);
std::pair<double, double> compMEtProjU(double, double, double, double, int&, bool = true);

std::vector<double> compTrackPtSums(const reco::VertexCollection&);
size_t getNumVerticesPtGtThreshold(const std::vector<double>&, double);

#endif
