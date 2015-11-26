#include "DesyTauAnalyses/NTupleMaker/interface/genMatch.h"

#include "DesyTauAnalyses/CandidateTools/interface/candidateAuxFunctions.h"

int utils_genMatch::genMatch( const reco::Candidate::LorentzVector& p4, const std::vector<const reco::GenParticle*>& genPart){


  //if (!(*daughter)->statusFlags().isPrompt()) continue;
  return 0;
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
