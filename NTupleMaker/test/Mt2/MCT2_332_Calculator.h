// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MCT2_332_CALCULATOR_H
#define MCT2_332_CALCULATOR_H
#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

namespace Mt2 {
/**
   Base class for classes which know how to calculate M_CT2 of http://arxiv.org/abs/0912.2354
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date April 2010 and onwards
*/

  class MCT2_332_Calculator : public Mt2Calculator {
  public:
    /** 
	mct2_332
    */
    virtual double mct2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			    const LorentzTransverseVector& visibleB, // 3 d.o.f.
			    const TwoVector& ptmiss,                 // 2 d.o.f.
			    double mInvisible) = 0;
    virtual double mct2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			       const LorentzTransverseVector& visibleB, // 3 d.o.f.
			       const TwoVector& ptmiss,                 // 2 d.o.f.
			       double mInvisible) = 0;
    MCT2_332_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
 
}

#endif
