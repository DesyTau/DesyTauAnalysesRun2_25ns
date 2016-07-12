// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MCT_330_CALCULATORS_H
#define MCT_330_CALCULATORS_H

#include "Mt2Vectors.h"

namespace Mt2 {
/**
   Class which knows how to calculate Tovey's M_CT -- the Contransverse Mass.
   (Tovey,http://arxiv.org/abs/0802.2879 )
 
   @author  Chris Lester
*/

  double mct_330   (const LorentzTransverseVector& visibleA,  // 3 d.o.f. 
		    const LorentzTransverseVector& visibleB); // 3 d.o.f.
  double mct_330_Sq(const LorentzTransverseVector& visibleA,  // 3 d.o.f. 
		    const LorentzTransverseVector& visibleB); // 3 d.o.f.
  
   
}

#endif 
