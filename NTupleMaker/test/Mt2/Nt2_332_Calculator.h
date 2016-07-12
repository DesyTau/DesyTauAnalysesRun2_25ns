// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef NT2_332_CALCULATOR_H
#define NT2_332_CALCULATOR_H
#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate N_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Nt2_332_Calculator : public Mt2Calculator {
  public:
    /** 
	nt2_332
	
	Calculate N_T2 knowing ptmiss not pvis-transverse-lorentz vec
	
    */
    virtual double nt2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			   const LorentzTransverseVector& visibleB, // 3 d.o.f.
			   const TwoVector& ptmiss,                 // 2 d.o.f.
			   const double mHeavy1, const double mHeavy2) = 0;
    virtual double nt2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			      const LorentzTransverseVector& visibleB, // 3 d.o.f.
			      const TwoVector& ptmiss,                 // 2 d.o.f.
			      const double mHeavy1, const double mHeavy2) = 0;
    Nt2_332_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
 
}

#endif
