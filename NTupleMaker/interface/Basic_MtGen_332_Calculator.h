// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_MtGen_332_Calculator_H
#define Basic_MtGen_332_Calculator_H

#include "MtGen_332_Calculator.h"
#include "Mt2_332_Calculator.h"
#include "Mt2_302_Calculator.h"

#include <vector>
#include <iostream>

namespace Mt2 {

  class Basic_MtGen_332_Calculator : public MtGen_332_Calculator {
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
       if (n==1) {
         static Mt2_302_Calculator calc;
	 const double ans = calc.mt2_302(LorentzTransverseVector(interestingParticles[0]),
					    /* zero, */
					     ptMiss,
					     mChi);
         m_lastSolutionType = calc.lastSolutionType();
         return ans;
       }
       
       // Now carry on assuming n>=2 ..
       {
	 static bool firstWarning = true;
	 if (firstWarning && n>20) {
	   firstWarning = false;
	   std::cerr << "Warning, calling " << __FUNCTION__ << " with more than 20 particles may take a long time on many computers!  Perhaps you wanted to use a different MTGEN evaluator than [" << this->algorithmName() << "]" <<  std::endl;
	 }
       }
       unsigned long mt2evaluations=0;
       std::vector<bool> combinationsPlusTwoOverTwo(n);
       bszero(combinationsPlusTwoOverTwo);
       combinationsPlusTwoOverTwo[n-1]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may want to change that!!!
       double mtGen=0;
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
		//  std::cout << "A";
                //} 
            } else {
                b = b + interestingParticles[bit];
		//if (debug()) {
		//  std::cout << "b";
                //} 
            }
	  }
          //if(debug()) { std::cout << std::endl; }
	  const double mt2 = m_mt2_332_Calculator.mt2_332(LorentzTransverseVector(a),LorentzTransverseVector(b),ptMiss,mChi);
	  ++mt2evaluations;
          if (first || mt2<mtGen) {
	    first=false;
            mtGen = mt2;
            m_lastSolutionType = m_mt2_332_Calculator.lastSolutionType();
          }
       }
       if (debug()) {
	 std::cerr << "Used " << mt2evaluations << " mt2 evaluations" << std::endl;
       }
       return mtGen;
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
       if (n==1) {
         static Mt2_302_Calculator calc;
	 const double ans = calc.mt2_302(interestingParticles[0],
					    /* zero, */
					     ptMiss,
					     mChi);
         m_lastSolutionType = calc.lastSolutionType();
         return ans;
       }
       
       // Now carry on assuming n>=2 ..
       {
	 static bool firstWarning = true;
	 if (firstWarning && n>20) {
	   firstWarning = false;
	   std::cerr << "Warning, calling " << __FUNCTION__ << " with more than 20 particles may take a long time on many computers!  Perhaps you wanted to use a different MTGEN evaluator than [" << this->algorithmName() << "]" <<  std::endl;
	 }
       }
       unsigned long mt2evaluations=0;
       std::vector<bool> combinationsPlusTwoOverTwo(n);
       bszero(combinationsPlusTwoOverTwo);
       combinationsPlusTwoOverTwo[n-1]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may want to change that!!!
       double mtGen=0;
       bool first=true;
       std::vector<bool> comb(n);
       bszero(comb);
       comb[0]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may want to change that!!!
       for (;
	    comb != combinationsPlusTwoOverTwo;
	    bsincr(comb)) {
	  LorentzTransverseVector a;
	  LorentzTransverseVector b;
          for (unsigned int bit=0; bit<n; ++bit) {
            if (comb[bit]) {
		a = a + interestingParticles[bit];
		//if (debug()) {
		//  std::cout << "A";
                //} 
            } else {
                b = b + interestingParticles[bit];
		//if (debug()) {
		//  std::cout << "b";
                //} 
            }
	  }
          //if(debug()) { std::cout << std::endl; }
	  const double mt2 = m_mt2_332_Calculator.mt2_332(a,b,ptMiss,mChi);
	  ++mt2evaluations;
          if (first || mt2<mtGen) {
	    first=false;
            mtGen = mt2;
            m_lastSolutionType = m_mt2_332_Calculator.lastSolutionType();
          }
       }
       if (debug()) {
	 std::cerr << "Used " << mt2evaluations << " mt2 evaluations" << std::endl;
       }
       return mtGen;
    }
    Basic_MtGen_332_Calculator(Mt2_332_Calculator & mt2_332_calc) : 
      MtGen_332_Calculator("Basic_MtGen_332_Calculator"),
      m_mt2_332_Calculator(mt2_332_calc) {
    }
  private:
    Mt2_332_Calculator & m_mt2_332_Calculator;
  };

} //end namespace Mt2

#endif
