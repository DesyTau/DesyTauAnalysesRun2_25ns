#ifndef NTupleMakerGenMatch_h
#define NTupleMakerGenMatch_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

namespace utils_genMatch{

  int genMatch( const reco::Candidate::LorentzVector& p4, const std::vector<reco::GenParticle>& genPart);
    
  reco::Candidate::LorentzVector getVisMomentumNoLep(const std::vector<const reco::GenParticle*>&, int = 1);
  reco::Candidate::LorentzVector getVisMomentumNoLep(const reco::GenParticle*);
}

#endif
