// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MC_330_CALCULATORS_H
#define MC_330_CALCULATORS_H

#include "Mt2Vectors.h"

namespace Mt2 {
/**
   Class which knows how to calculate Tovey's M_C -- the Contralinear-Boost Invariant Mass.
   (Tovey,http://arxiv.org/abs/0802.2879 )
 
   @author  Chris Lester
*/

  double mc_330   (const LorentzVector& visibleA,  // 3 d.o.f. 
		   const LorentzVector& visibleB); // 3 d.o.f.
  double mc_330_Sq(const LorentzVector& visibleA,  // 3 d.o.f. 
		   const LorentzVector& visibleB); // 3 d.o.f.
  
   
}

#endif 
