// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Analytic_Mt2_330_Calculator_H
#define Analytic_Mt2_330_Calculator_H

#include "SolutionType.h"
#include "Mt2_330_Calculator.h"

namespace Mt2 {

  class Analytic_Mt2_330_Calculator : public Mt2_330_Calculator {
  public:
    double mt2_330(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		   const LorentzTransverseVector& visibleB, // 3 d.o.f.
		                    // 0 d.o.f. for ptmiss ... it must be opposite to whavever visible stuff was seen!
		   const double mInvisible);
    double mt2_330_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		      const LorentzTransverseVector& visibleB, // 3 d.o.f.
		                      // 0 d.o.f. for ptmiss ... it must be opposite to whavever visible stuff was seen!
		     const double mInvisible);
    Analytic_Mt2_330_Calculator() : Mt2_330_Calculator("Analytic_Mt2_330") {};
    struct Notes {
      SolutionType soln;
      double rootS;
      LorentzTransverseVector gamma;
      Notes() { clear(); }
      void clear() {
	soln = SolutionType(NotSpecified);
        rootS=0;
	gamma = LorentzTransverseVector();
      }	
    };
    Notes notes() const {
	return m_notes;
    }
  private:
    Notes m_notes;
  };

} //end namespace Mt2

#endif
