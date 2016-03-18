#include "DesyTauAnalyses/NTupleMaker/interface/genMatch.h"

#include "DesyTauAnalyses/CandidateTools/interface/candidateAuxFunctions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/genMatch.h"
#include "Math/GenVector/VectorUtil.h"

int utils_genMatch::genMatch( const reco::Candidate::LorentzVector& p4, const std::vector<reco::GenParticle>& genParticles){
  int match = 6;
  float dr = 0.2;
  float dr_new = 0.2;
  
  for ( std::vector<reco::GenParticle>::const_iterator genp = genParticles.begin();
	genp != genParticles.end(); ++genp ) {
    dr_new = ROOT::Math::VectorUtil::DeltaR( p4, genp->p4());

    if (dr_new < dr){
      if (TMath::Abs(genp->pdgId()) == 11) {
	if (genp->pt() > 8) {
	  if (genp->statusFlags().isPrompt()) {
	    match = 1;
	    dr = dr_new;
	  }
	  else if (genp->statusFlags().isDirectPromptTauDecayProduct()) {
	    match = 3;
	    dr = dr_new;
	  }
	}
      }
      else if (TMath::Abs(genp->pdgId()) == 13) {
	if (genp->pt() > 8) {
	  if (genp->isPromptFinalState()){
	    match = 2;
	    dr = dr_new;
	  }
	  if (genp->isDirectPromptTauDecayProductFinalState()){
	    match = 4;
	    dr = dr_new;
	  }
	}
      }
      else if (TMath::Abs(genp->pdgId()) == 15) {
	if (genp->statusFlags().isPrompt()) {
	  reco::Candidate::LorentzVector tau_p4 = utils_genMatch::getVisMomentumNoLep(&(*genp));
	  if (tau_p4.pt() > 15) {
	    match = 5;
	    dr = dr_new;
	  }
	}
      }
    }
  }

  return match;
}

reco::Candidate::LorentzVector utils_genMatch::getVisMomentumNoLep(const std::vector<const reco::GenParticle*>& daughters, int status)
{
  reco::Candidate::LorentzVector p4Vis(0,0,0,0);

  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if (status != -1 && (*daughter)->status() != status) continue;
    if (isNeutrino(*daughter)) continue;
    if (TMath::Abs((*daughter)->pdgId()) == 11) continue;
    if (TMath::Abs((*daughter)->pdgId()) == 13) continue;
    //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
    //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
    p4Vis += (*daughter)->p4();
  }

  //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;

  return p4Vis;
}

reco::Candidate::LorentzVector utils_genMatch::getVisMomentumNoLep(const reco::GenParticle* genLeg)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(genLeg, stableDaughters, 1);

  reco::Candidate::LorentzVector p4Vis = utils_genMatch::getVisMomentumNoLep(stableDaughters);

  return p4Vis;
}
