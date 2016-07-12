// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_300_CALCULATOR_H
#define MT2_300_CALCULATOR_H

#include "Mt2Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate MT2 for the special case in which ONE
   of the visible particles has ZERO lorentz-1+2 momentum, and which the PTMiss
   is ENTIRELY OPPOSITE to the visible paticle.  Actually the answer is in fact
   also general and would be the same EVEN IF the ptmiss were to point somewhere
   else.

   This not something one would be able to measure in a real detector, but is
   an idealisation which is useful for other MT2 calculators internally.

   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2_300_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_300
    */
    double mt2_300(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			  // 0 d.o.f. : Momentum of visibleB is taken to be zero!
			  // 0 d.o.f. : Momentum of ptmiss is taken to be opposite to visibleA !
			  double mInvisible) {
      // There is an analytic answer for this, so we may as well write it in here!
      // See for example (6) on 4/7/2007 in lester lab book 6.
      // (The answer would change with a "root S" constraint for small visibleA mass and large mChi)
      m_lastSolutionType = SolutionType(ASideGlobalMin);
      return visibleA.mass() + mInvisible;
    }
    double mt2_300_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			     // 0 d.o.f. : Momentum of visibleB is taken to be zero!
			     // 0 d.o.f. : Momentum of ptmiss is taken to be opposite to visibleA !
			     double mInvisible) {
      const double m = mt2_300(visibleA, mInvisible);
      m_lastSolutionType = SolutionType(ASideGlobalMin);
      return m*m;
    }
    Mt2_300_Calculator(const std::string & algName="Mt2_300_Calculator") : Mt2Calculator(algName) {};
  };
 
}

#endif
