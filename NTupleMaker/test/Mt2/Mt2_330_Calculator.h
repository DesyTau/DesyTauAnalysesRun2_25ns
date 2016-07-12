// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_330_CALCULATORS_H
#define MT2_330_CALCULATORS_H

#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate M_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2_330_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_330
	
	Calculate M_T2 for events in which PTMISS is exactly
        balanced by the two visible momenta fed into MT2 ...
        such as may be the case in MTGEN type events.
	
    */
    virtual double mt2_330(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			   const LorentzTransverseVector& visibleB, // 3 d.o.f.
			                    // 0 d.o.f for ptmiss ... it's calculable from minus the transverse cpts of visiableA and visibleB in this case.
			   const double mInvisible) = 0;
    virtual double mt2_330_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			      const LorentzTransverseVector& visibleB, // 3 d.o.f.
			                      // 0 d.o.f for ptmiss ... it's calculable from minus the transverse cpts of visiableA and visibleB in this case.
			      const double mInvisible) = 0;
    Mt2_330_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
 
}

#endif 
