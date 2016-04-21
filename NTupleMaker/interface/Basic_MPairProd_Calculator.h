// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_MPairProd_Calculator_H
#define Basic_MPairProd_Calculator_H

#include "MPairProd_Calculator.h"
#include <vector>

namespace Mt2 {
/**
   Class which knows how to calculate a quantity like MT2 but without the missing particles and using all four dimensions.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Basic_MPairProd_Calculator : public MPairProd_Calculator {
  private:
    static void bszero(std::vector<bool> & bs) {
      for (unsigned int i=0; i<bs.size(); ++i) {
	bs[i]=false;
      }
    }
    static void bsincr(std::vector<bool> & bs) {
      for(unsigned long bit=0; bit<bs.size(); ++bit) {
	bs[bit] = (!(bs[bit]));
	if (bs[bit]) {
	  // we incremented 0 to 1, so we are done
	  return;
	}
	// we incremented 1 to 0, so we had better add the carry bit...
      } 
      // we overflowed, but do not complain
    }
  public:
    virtual double mPairProd(const std::vector<LorentzVector> & interestingParticles) const {
      const unsigned int n=interestingParticles.size();
      if (debug()) {
	std::cerr << __FUNCTION__ << " called in " << this->algorithmName() << " with " << n << " interesting particles." << std::endl;
      }
      if (n<=1) { return 0; }
      {
	static bool firstWarning = true;
	if (firstWarning && n>20) {
	  firstWarning = false;
	  std::cerr << "Warning, calling " << __FUNCTION__ << " with more than 20 particles may take a long time on many computers!  Perhaps you wanted to use a different MPairProd evaluator than [" << this->algorithmName() << "]" <<  std::endl;
	}
      }
      std::vector<bool> combinationsPlusTwoOverTwo(n);
      bszero(combinationsPlusTwoOverTwo);
      combinationsPlusTwoOverTwo[n-1]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may      
      double mPairProd2=0;
      bool first=true;
      std::vector<bool> comb(n);
      bszero(comb);
      comb[0]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may want to change that!!!
      for (;
	   comb != combinationsPlusTwoOverTwo; 
	   bsincr(comb)) {
	LorentzVector a;
	LorentzVector b;
	for (unsigned int bit=0; bit<n; ++bit) {
	  if (comb[bit]) {
	    a = a + interestingParticles[bit];
	    //if (debug()) {
	    //std::cout << "A";
	    //} 
	  } else {
	    b = b + interestingParticles[bit];
	    //if (debug()) {
	    //		std::cout << "b";
	    //} 
	  }
	}
	const double mA2 = a.m2();
	const double mB2 = b.m2();
	const double mMax2 = (mA2>mB2) ? mA2 : mB2; 
	//if(debug()) { std::cout << " " << sqrt(mMax2) << std::endl; }
	
	if (first || (mMax2 >=0 && mMax2<mPairProd2)) {
	  first = false;
	  mPairProd2 = mMax2;
	}
      }
      if (debug()) {
	std::cerr << "Used approx " << pow(2.0,(n-1.0)) << " iterations" << std::endl;
      }
      return sqrt(mPairProd2);
    }
    Basic_MPairProd_Calculator(const std::string & algName = "Basic_MPairProd_Calculator") : MPairProd_Calculator(algName) {};
  };

}

#endif
