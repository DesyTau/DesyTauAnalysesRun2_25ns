// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Advanced_Mt2_332_Calculator_H
#define Advanced_Mt2_332_Calculator_H

#include "Mt2/Mt2_332_Calculator.h"
#include <vector>

namespace Mt2 {

  class Advanced_Mt2_332_Calculator : public Mt2_332_Calculator {
  public:
    double mt2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		   const LorentzTransverseVector& visibleB, // 3 d.o.f.
		   const TwoVector& ptmiss,                 // 2 d.o.f.
		   const double mInvisible);
    double mt2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		      const LorentzTransverseVector& visibleB, // 3 d.o.f.
		      const TwoVector& ptmiss,                 // 2 d.o.f.
		      const double mInvisible);

     // You can enable the fancyAssistor if you want to play with it, but the fancy assistor is still a bit naughty about how it deals with the limit of massless visible particles, so the fancy assistor is disabled by default. See more comments inside the implementation of the  mt2_332_Sq method. 

    Advanced_Mt2_332_Calculator(const bool useFancyAssistor = false) : Mt2_332_Calculator("Advanced_Mt2_332"), m_useFancyAssistor(useFancyAssistor) {} 
  private:
    const bool m_useFancyAssistor;
  };

} //end namespace Mt2

#endif
