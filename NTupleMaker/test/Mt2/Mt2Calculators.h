// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_CALCULATORS_H
#define MT2_CALCULATORS_H
#include "Mt2Vectors.h"
#include "Mt2Calculator.h"
#include "Mt2Minimiser.h"
#include "Mt2_332_Calculator.h"

#include <string>
#include <vector>

namespace Mt2 {
/**
   Class which knows how to calculate M_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class MtGen_PossiblySomeJunkInEvent_Calculator : public Mt2Calculator {
  public:
    virtual double mtGen(const std::vector<LorentzTransverseVector> & interestingParticles,
		         const TwoVector & ptmiss,
                         const double mChi) const = 0;
    MtGen_PossiblySomeJunkInEvent_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
   
  class Mt2_4441_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_4441
	
	Calculate M_T2 knowing visible energy - in principle better than
	the method above, as more information is known.  In practice there
	is very little difference between the two methods, however.
    */
    virtual double mt2_4441(const Mt2::LorentzVector& visibleA,      // 4 d.o.f. 
			    const Mt2::LorentzVector& visibleB,      // 4 d.o.f. 
			    const Mt2::LorentzVector& totalVisible,  // 4 d.o.f.
			    const double rootSMax, /* eg 14 TeV */ // 1 d.o.f
			    double mInvisible) = 0;
    Mt2_4441_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
  
  class Mt2_222_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_222
	
	Calculate M_T2 knowing ptmiss not pvis-transverse-lorentz vec,
	but assuming visible particles to be massless (or of
	negligible mass)
    */
    virtual double mt2_222(const TwoVector& visibleA, // 2 d.o.f. 
			   const TwoVector& visibleB, // 2 d.o.f.
			   const TwoVector& ptmiss,   // 2 d.o.f.
			   double mInvisible) = 0;
    Mt2_222_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
  
  class Mt2_2220_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_2220
	
	Calculate M_T2 knowing ptmiss not pvis-transverse-lorentz vec,
	but assuming visible particles to be massless (or of
	negligible mass) and assuming chi=0,
    */
    virtual double mt2_2220(const TwoVector& visibleA, // 2 d.o.f. 
			    const TwoVector& visibleB, // 2 d.o.f.
			    const TwoVector& ptmiss   // 2 d.o.f.
			    // 0 d.o.f in chi!
			    ) = 0;
    virtual double mt2_2220_Sq(const TwoVector& visibleA, // 2 d.o.f. 
			       const TwoVector& visibleB, // 2 d.o.f.
			       const TwoVector& ptmiss   // 2 d.o.f.
			       // 0 d.o.f in chi!
			       ) = 0;
    Mt2_2220_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
  
  class Suspect_Mt2_332_Calculator : public Mt2_332_Calculator {
  public:
    double mt2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		   const LorentzTransverseVector& visibleB, // 3 d.o.f.
		   const TwoVector& ptmiss,                 // 2 d.o.f.
		   const double mInvisible);
    
    Suspect_Mt2_332_Calculator(const std::string & algName) : Mt2_332_Calculator(algName) {};
  };

  class Basic_Mt2_4441_Calculator : public Mt2_4441_Calculator {
  public:
    double mt2_4441(const Mt2::LorentzVector& visibleA,      // 4 d.o.f. 
		    const Mt2::LorentzVector& visibleB,      // 4 d.o.f. 
		    const Mt2::LorentzVector& otherVisible,  // 4 d.o.f.
		    const double rootSMax, /* eg 14 TeV */ // 1 d.o.f
		    const double mInvisible);
    Basic_Mt2_4441_Calculator() : Mt2_4441_Calculator("Basic_Mt2_4441") {};

  private:
    double mt2_4441(const Mt2::LorentzVector& visibleA,      // 4 d.o.f. 
		    const Mt2::LorentzVector& visibleB,      // 4 d.o.f. 
		    const Mt2::LorentzVector& otherVisible,  // 4 d.o.f.
		    const double rootSMax, /* eg 14 TeV */ // 1 d.o.f
		    const double mInvisible,
		    const double sign);
  };

  /*
  // Data and methods below this point should only be used by 
  // the other classes in the package

  double m_tolerance;
  static bool m_debug;

  // input variables
  LorentzTransverseVector alpha;     /// first visible particle
  LorentzTransverseVector beta;      /// second visible particle
  LorentzTransverseVector bigSigma;  /// total visible
  double mlspsq;
  
  // derrived quantities
  LorentzTransverseVector sigma;     /// sum of two visibles
  LorentzTransverseVector delta;     /// difference of two visibles
  LorentzTransverseVector lambda;    /// lab frame
  
  // dot products of input variables
  double lambdaDotBigSigma;
  double lambdaDotSigma;
  double lambdaDotDelta;
  double deltaDotBigSigma;
  double deltaDotSigma;
  double sigmaDotBigSigma;
  
  double bigSigmaSq;
  double sigmaSq;
  double deltaSq;
  double alphaSqPlusBetaSq;

  /// minimum CMF energy;
  double rootSMin;
  double rootSMax;

  //// 
  //// calculate Mt2 at fixed value of rootS = rootSMin + par^2
  ////
  double trialMt2(double contribution) const;
  
  void printVariables() const;

  
  mutable unsigned m_calls;
  mutable unsigned m_evaluations;
  mutable unsigned m_unbalancedEvaluations;

 private:

  Mt2Minimiser* m_minimiser; /// abstract class which does the minimisation
  static bool m_initialised;
private:
  static bool m_initialised;
  static bool m_debug;
};
  */
 
}

/*
std::ostream& operator << (std::ostream& os , const Mt2::LorentzTransverseVector& v);

std::ostream& operator << (std::ostream& os , const Mt2::TwoVector& v);
*/

#endif
