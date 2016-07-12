// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MPairProd_Calculator_H
#define MPairProd_Calculator_H

#include "Mt2Vectors.h"
#include "Mt2Calculator.h"

#include <vector>

namespace Mt2 {
/**
   Class which knows how to calculate a quantity like MT2 but without the missing particles and using all four dimensions.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class MPairProd_Calculator : public Mt2Calculator {
  public:
    virtual double mPairProd(const std::vector<LorentzVector> & interestingParticles) const = 0;
    MPairProd_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };

}

#endif
