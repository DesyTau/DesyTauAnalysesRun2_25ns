// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MtGen_330_Calculator_H
#define MtGen_330_Calculator_H

#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

#include <vector>

namespace Mt2 {
/**
   Class which knows how to calculate M_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class MtGen_330_Calculator : public Mt2Calculator {
  public:
    virtual double mtGen_330(const std::vector<LorentzTransverseVector> & interestingParticles,
			     const double mChi) = 0;
    virtual double mtGen_330(const std::vector<LorentzVector> & interestingParticles,
			     const double mChi) = 0;
    MtGen_330_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };

}

#endif
