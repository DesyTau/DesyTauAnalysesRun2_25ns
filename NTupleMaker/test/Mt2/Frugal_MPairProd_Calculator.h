// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Frugal_MPairProd_Calculator_H
#define Frugal_MPairProd_Calculator_H

#include "MPairProd_Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate a quantity like MT2 but without the missing particles and using all four dimensions.  Supposed to be equivalent to Basic_MPairProd_Calculator, but using fewer operations.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Frugal_MPairProd_Calculator : public MPairProd_Calculator {
  private:
    static unsigned long s_partitions;
  public:
    virtual double mPairProd(const std::vector<LorentzVector> & interestingParticles) const {
      const unsigned int n=interestingParticles.size();
      if (debug()) {
	std::cerr << __FUNCTION__ << " called in " << this->algorithmName() << " with " << n << " interesting particles." << std::endl;
      }
      if (n==0) { return 0; } // not really needed ... but will save some time
      s_partitions = 0;

      //Find total mom
      const LorentzVector & totalMom = total(interestingParticles);

      //create a variable to store the best MPP value so far:
      double mppSq = totalMom.m2(); // which you get by assigning nothing to side A, and everything to side B.  i.e.:
      const LorentzVector presentSideA;

      // Now make a note to try adding things to side a starting at position 0 in the vector, and see if anything better is found.
      tryAlsoFrom(0, interestingParticles, presentSideA, totalMom, mppSq);
      if (debug()) {
	std::cerr << "Probed " << s_partitions << " partitions" << std::endl;
      }
      return sqrt(mppSq);
    }
    Frugal_MPairProd_Calculator(const std::string & algName = "Frugal_MPairProd_Calculator") : MPairProd_Calculator(algName) {};
  private:
    static void tryAlsoFrom(const unsigned int startHere,
			    const std::vector<LorentzVector> & interestingParticles,
			    const LorentzVector & previousSideA,
			    const LorentzVector & totalMom, 
			    double & bestMppSqSoFar) {
      for (unsigned int index /*to add to side A*/ = startHere;
	   index != interestingParticles.size();
	   ++index) {
	++s_partitions;
	const LorentzVector & currentSideA = previousSideA + interestingParticles[index];
	const double mA2 = currentSideA.m2();
	if (mA2 >= bestMppSqSoFar) {
	  //This partition makes A too heavy, so forget it.  Futhermore, don't bother adding anything extra to A as this will only make A heavier.  So just continue to an entirely different side A
	  continue;
	}
	const LorentzVector & currentSideB = totalMom - currentSideA;
	const double mB2 = currentSideB.m2();
	const bool aIsLargest = (mA2 >= mB2);
	if (mB2 >= bestMppSqSoFar) {
	  // THIS partition is not a winner, but we may yet get one by adding more things to side A,
	  // so unlike the A case above, we DON'T "continue" here !
	} else {
	  // we have a new leader
	  bestMppSqSoFar = (aIsLargest ? mA2 : mB2);
	}
	if (aIsLargest) {
	  // no point in trying to add anything further to the current A as this will only increase the mass of side A and we are trying to find a MINIMUM invariant mass over all splittings.  However, we can try different things with the PREVIOUS side A.
	  continue;
	}
	// If we are here, then b was largest.  We should try adding further things to the current side A as this will decrease the invariant mass of side B, assuming there is anything left to 
	tryAlsoFrom(index+1, // when there are no particles left to try (index+1 == interestingParticles.size(), the function call will return straight away
		    interestingParticles,
		    currentSideA,
		    totalMom,
		    bestMppSqSoFar);
      }
    }		    
    static LorentzVector total(const std::vector<LorentzVector> & particles) {
      LorentzVector v;
      for (unsigned int i=0; i<particles.size(); ++i) {
	v += particles[i];
      }
      return v;
    }
  };

  unsigned long Frugal_MPairProd_Calculator::s_partitions;
}

#endif
