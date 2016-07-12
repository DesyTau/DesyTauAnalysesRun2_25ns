// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_302_CALCULATOR_H
#define MT2_302_CALCULATOR_H

#include "Mt2_300_Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate MT2 for the special case in which ONE
   of the visible particles has ZERO lorentz-1+2 momentum.

   This not something one would be able to measure in a real detector, but is
   an idealisation which is useful for other MT2 calculators internally.

   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2_302_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_302
    */
    double mt2_302(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			  // 0 d.o.f. : Momentum of visibleB is taken to be zero!
			  const TwoVector & ptmiss,                // 2 d.o.f.   ... does not actually get used!
			  double mInvisible) {
      // There is an analytic answer for this, so ww could write it in here!
      // See for example (6) on 4/7/2007 in lester lab book 6.
      // (The answer would change with a "root S" constraint for small visibleA mass and large mChi)
      // However since the answer doesn't even use ptmiss, and I don't think it is good to type the
      // same expression in too many placed independently, we will grab the result from the simplest
      // implementation it is equivalent to:
      static Mt2_300_Calculator calc;
      const double ans = calc.mt2_300(visibleA, mInvisible);
      m_lastSolutionType = calc.lastSolutionType();
      return ans;
    }
    double mt2_302_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			     // 0 d.o.f. : Momentum of visibleB is taken to be zero!
			     const TwoVector& ptmiss,                 // 2 d.o.f. ... does not actually get used!
			     double mInvisible) {
      // There is an analytic answer for this, so ww could write it in here!
      // See for example (6) on 4/7/2007 in lester lab book 6.
      // (The answer would change with a "root S" constraint for small visibleA mass and large mChi)
      // However since the answer doesn't even use ptmiss, and I don't think it is good to type the
      // same expression in too many placed independently, we will grab the result from the simplest
      // implementation it is equivalent to:
      static Mt2_300_Calculator calc;
      const double ans = calc.mt2_300_Sq(visibleA, mInvisible);
      m_lastSolutionType = calc.lastSolutionType();
      return ans;
    }
    Mt2_302_Calculator(const std::string & algName="Mt2_302_Calculator") : Mt2Calculator(algName) {};
  };
 
}

#endif
