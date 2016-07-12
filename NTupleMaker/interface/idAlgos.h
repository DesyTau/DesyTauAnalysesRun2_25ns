#ifndef NTupleMakerIdAlgos_h
#define NTupleMakerIdAlgos_h

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

namespace idAlgos{

  bool isICHEPMuon(const reco::Muon & recoMu) 
   {
      bool goodGlob = recoMu.isGlobalMuon() && 
                      recoMu.globalTrack()->normalizedChi2() < 3 && 
                      recoMu.combinedQuality().chi2LocalPosition < 12 && 
                      recoMu.combinedQuality().trkKink < 20; 
      bool isMedium = muon::isLooseMuon(recoMu) && 
                      recoMu.innerTrack()->validFraction() > 0.49 && 
                      muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
      return isMedium; 
   }
}

#endif
