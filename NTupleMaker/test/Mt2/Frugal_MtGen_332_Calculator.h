// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Frugal_MtGen_332_Calculator_H
#define Frugal_MtGen_332_Calculator_H

#include "MtGen_332_Calculator.h"
#include "Mt2_332_Calculator.h"
#include "Mt2_302_Calculator.h"

namespace Mt2 {
/**
   Class which knows how to calculate a quantity like MT2 but without the missing particles and using all four dimensions.  Supposed to be equivalent to Basic_MtGen_332_Calculator, but using fewer operations.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Frugal_MtGen_332_Calculator : public MtGen_332_Calculator {
  private:
    mutable unsigned long m_partitions;
    mutable unsigned long m_mt2evaluations;

  public:
    virtual double mtGen_332(const std::vector<LorentzVector> & interestingParticles,
                             const TwoVector & ptMiss,
			     const double mChi) {
      const unsigned int n=interestingParticles.size();
      if (debug()) {
	std::cerr << __FUNCTION__ << " called in " << this->algorithmName() << " with " << n << " interesting particles." << std::endl;
      }
      if (n==0) {
        m_lastSolutionType = SolutionType(BalancedSolution);
	return mChi;
      }

      //Find total mom
      const LorentzVector & totalMom = total(interestingParticles);

      // First interesting particle will always be on SOME side.
      // WLOG we choose to name the side that IT is on as "side B".
      // One example of this is ALL ON SIDE B and NOTHING ON SIDE A.
      /// We will use this example to set things going.
      const LorentzVector   zeroMom;
      const LorentzVector & presentSideA = zeroMom;
      const LorentzVector & presentSideB = totalMom;

      m_partitions = 1;
      m_mt2evaluations = 1;
      static Mt2_302_Calculator calc;
      double mtGenSq = calc.mt2_302_Sq(LorentzTransverseVector(presentSideB), 
						      /* presentSideA = zero */
						      ptMiss,
						      mChi);
      m_lastSolutionType =  calc.lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << " " << m_lastSolutionType << std::endl;



      // Now make a note to try adding things to side A to see if anything better is found.
      tryAlsoFrom(1, // no need to start at 0 .. we NEVER want the 0th particle on side A by definition.
		  interestingParticles,
		  ptMiss,
		  presentSideA,
		  totalMom,
		  mChi,
		  mtGenSq);
      if (debug()) {
	std::cerr << "Probed " << m_partitions << " partitions " << std::endl;
	std::cerr << "Made   " << m_mt2evaluations << " mt2 evaluations " << std::endl;
      }
      return sqrt(mtGenSq);
    }
    virtual double mtGen_332(const std::vector<LorentzTransverseVector> & interestingParticles,
                             const TwoVector & ptMiss,
			     const double mChi) {
      const unsigned int n=interestingParticles.size();
      if (debug()) {
	std::cerr << __FUNCTION__ << " called in " << this->algorithmName() << " with " << n << " interesting particles." << std::endl;
      }
      if (n==0) {
        m_lastSolutionType = SolutionType(BalancedSolution);
	return mChi;
      }

      //Find total mom
      const LorentzTransverseVector & totalMom = total(interestingParticles);

      // First interesting particle will always be on SOME side.
      // WLOG we choose to name the side that IT is on as "side B".
      // One example of this is ALL ON SIDE B and NOTHING ON SIDE A.
      /// We will use this example to set things going.
      const LorentzTransverseVector   zeroMom;
      const LorentzTransverseVector & presentSideA = zeroMom;
      const LorentzTransverseVector & presentSideB = totalMom;

      m_partitions = 1;
      m_mt2evaluations = 1;
      static Mt2_302_Calculator calc;
      double mtGenSq = calc.mt2_302_Sq(presentSideB, 
						      /* presentSideA = zero */
						      ptMiss,
						      mChi);
      m_lastSolutionType =  calc.lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << " " << m_lastSolutionType << std::endl;



      // Now make a note to try adding things to side A to see if anything better is found.
      tryAlsoFrom(1, // no need to start at 0 .. we NEVER want the 0th particle on side A by definition.
		  interestingParticles,
		  ptMiss,
		  presentSideA,
		  totalMom,
		  mChi,
		  mtGenSq);
      if (debug()) {
	std::cerr << "Probed " << m_partitions << " partitions " << std::endl;
	std::cerr << "Made   " << m_mt2evaluations << " mt2 evaluations " << std::endl;
      }
      return sqrt(mtGenSq);
    }
    Frugal_MtGen_332_Calculator(Mt2_332_Calculator & mt2_332_calc,
					  const std::string & algName = "Frugal_MtGen_332_Calculator") : 
      MtGen_332_Calculator(algName),
      m_mt2_332_Calculator(mt2_332_calc) {
    }
  private:
    Mt2_332_Calculator & m_mt2_332_Calculator;
  private:
    static inline double qAdd(const double aSq, const double mChi) {
      return aSq + 2.*sqrt(aSq)*mChi + mChi*mChi; // (a+b)^2
    }
    void tryAlsoFrom(const unsigned int startHere,
		     const std::vector<LorentzVector> & interestingParticles,
		     const TwoVector ptMiss,
		     const LorentzVector & previousSideA,
		     const LorentzVector & totalMom,
		     const double mChi,
		     double & bestMt2SqSoFar) {
      for (unsigned int index /*to add to side A*/ = startHere;
	   index != interestingParticles.size();
	   ++index) {
	++m_partitions;
	const LorentzVector & currentSideA = previousSideA + interestingParticles[index];
	const double mT2Sq_AMin = qAdd(currentSideA.masssq(), mChi);
	if (mT2Sq_AMin >= bestMt2SqSoFar) {
	  //This partition makes A too heavy, so forget it.  Futhermore, don't bother adding anything extra to A as this will only make A heavier.  So just continue to an entirely different side A
	  continue;
	}
	const LorentzVector & currentSideB = totalMom - currentSideA;
	const double mT2Sq_BMin = qAdd(currentSideB.masssq(), mChi);
	//const bool aIsLargest = (mA2 >= mB2);
	if (mT2Sq_BMin >= bestMt2SqSoFar) {
	  // THIS partition is not a winner, but we may yet get one by adding more things to side A,
	  // so unlike the A case above, we DON'T "continue" here !
	} else {
	  // We had better actually calculate MT2 now and see if we get a new winner.
	  const double mt2Sq = m_mt2_332_Calculator.mt2_332_Sq(LorentzTransverseVector(currentSideA), LorentzTransverseVector(currentSideB), ptMiss, mChi);
	  ++m_mt2evaluations;
	  if (mt2Sq < bestMt2SqSoFar) {
	    bestMt2SqSoFar = mt2Sq;
            m_lastSolutionType = m_mt2_332_Calculator.lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << " " << m_lastSolutionType << std::endl;
	  }
	  // might still get a better winner by adding more things to A ... due to the magic freedoms in the behaiour of the missing momenta splittings ... so don't continue!
	}

	tryAlsoFrom(index+1, // when there are no particles left to try (index+1 == interestingParticles.size(), the function call will return straight away
		    interestingParticles,
		    ptMiss,
		    currentSideA,
		    totalMom,
		    mChi,
		    bestMt2SqSoFar);
      }
    }		    
    void tryAlsoFrom(const unsigned int startHere,
		     const std::vector<LorentzTransverseVector> & interestingParticles,
		     const TwoVector ptMiss,
		     const LorentzTransverseVector & previousSideA,
		     const LorentzTransverseVector & totalMom,
		     const double mChi,
		     double & bestMt2SqSoFar) {
      for (unsigned int index /*to add to side A*/ = startHere;
	   index != interestingParticles.size();
	   ++index) {
	++m_partitions;
	const LorentzTransverseVector & currentSideA = previousSideA + interestingParticles[index];
	const double mT2Sq_AMin = qAdd(currentSideA.masssq(), mChi);
	if (mT2Sq_AMin >= bestMt2SqSoFar) {
	  //This partition makes A too heavy, so forget it.  Futhermore, don't bother adding anything extra to A as this will only make A heavier.  So just continue to an entirely different side A
	  continue;
	}
	const LorentzTransverseVector & currentSideB = totalMom - currentSideA;
	const double mT2Sq_BMin = qAdd(currentSideB.masssq(), mChi);
	//const bool aIsLargest = (mA2 >= mB2);
	if (mT2Sq_BMin >= bestMt2SqSoFar) {
	  // THIS partition is not a winner, but we may yet get one by adding more things to side A,
	  // so unlike the A case above, we DON'T "continue" here !
	} else {
	  // We had better actually calculate MT2 now and see if we get a new winner.
	  const double mt2Sq = m_mt2_332_Calculator.mt2_332_Sq(currentSideA, currentSideB, ptMiss, mChi);
	  ++m_mt2evaluations;
	  if (mt2Sq < bestMt2SqSoFar) {
	    bestMt2SqSoFar = mt2Sq;
            m_lastSolutionType = m_mt2_332_Calculator.lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << " " << m_lastSolutionType << std::endl;
	  }
	  // might still get a better winner by adding more things to A ... due to the magic freedoms in the behaiour of the missing momenta splittings ... so don't continue!
	}

	tryAlsoFrom(index+1, // when there are no particles left to try (index+1 == interestingParticles.size(), the function call will return straight away
		    interestingParticles,
		    ptMiss,
		    currentSideA,
		    totalMom,
		    mChi,
		    bestMt2SqSoFar);
      }
    }		    
    static LorentzVector total(const std::vector<LorentzVector> & particles) {
      LorentzVector v;
      for (unsigned int i=0; i<particles.size(); ++i) {
	v += particles[i];
      }
      return v;
    }
    static LorentzTransverseVector total(const std::vector<LorentzTransverseVector> & particles) {
      LorentzTransverseVector v;
      for (unsigned int i=0; i<particles.size(); ++i) {
	v += particles[i];
      }
      return v;
    }
  };
}

#endif
