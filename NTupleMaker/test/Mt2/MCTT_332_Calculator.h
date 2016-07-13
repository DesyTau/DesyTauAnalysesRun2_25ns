// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MCTT_332_CALCULATORS_H
#define MCTT_332_CALCULATORS_H

#include "Mt2Vectors.h"

namespace Mt2 {
/**
   This is "MCT Perp" from eq (11) from KM's arXiv:0910.1584v1
   @author  Chris Lester
*/

  double mctt_332   (const LorentzTransverseVector& visibleA,  // 3 d.o.f. 
		     const LorentzTransverseVector& visibleB,  // 3 d.o.f.
                     const TwoVector& utm);

  double mctt_332_Sq(const LorentzTransverseVector& visibleA,  // 3 d.o.f. 
		     const LorentzTransverseVector& visibleB,  // 3 d.o.f.
		     const TwoVector& utm);
   
}

#endif 
