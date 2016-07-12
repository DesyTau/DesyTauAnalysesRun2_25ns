// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_Asymmetric332_CALCULATOR_H
#define MT2_Asymmetric332_CALCULATOR_H
#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate the specialisation of M_T2 to the case of disimilar parent masses (but with a fixed ratio of parent masses).
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2_Asymmetric332_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_Asymmetric332 assumes that the heavy parent particles have different masses .... whose ratio is fixed in terms of the paramter theta according to: tan(theta) = mHeavyB/mHeavyA.
	
    */
    virtual double mt2_Asymmetric332(const double theta,  // tan(theta) is mHeavyB/mHeavyA. 0<=theta<=Pi/2
				const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			   const LorentzTransverseVector& visibleB, // 3 d.o.f.
			   const TwoVector& ptmiss,                 // 2 d.o.f.
			   double mInvisible) = 0;
    virtual double mt2_Asymmetric332_Sq(const double theta, // tan(theta) is mHeavyB/mHeavyA. 0<=theta<=Pi/2
				const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			      const LorentzTransverseVector& visibleB, // 3 d.o.f.
			      const TwoVector& ptmiss,                 // 2 d.o.f.
			      double mInvisible) = 0;
    Mt2_Asymmetric332_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
 
}

#endif
