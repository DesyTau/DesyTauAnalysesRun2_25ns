// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_MtGen_330_Calculator_H
#define Basic_MtGen_330_Calculator_H

#include "Mt2/MtGen_330_Calculator.h"
#include "Mt2/Mt2_330_Calculator.h"
#include "Basic_MtGen_332_Calculator.h"
#include "Mt2ApproximatingAdapter_332_from_330.h"

#include <iostream>

namespace Mt2 {

  class Basic_MtGen_330_Calculator : public MtGen_330_Calculator {
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
    // Want this:
    virtual double mtGen_330(const std::vector<LorentzVector> & interestingParticles, const double mChi) {
      //warning! The next line LOOKS like it is an error (it looks like we oughtto be setting ptMiss to minus the pt sum of the interesting particles) but in actual fact there is no error here.  The key is to understand that this class uses an "Mt2ApproximatingAdapter_332_from_330" to calculate MT2 values, and this sort of caclulator deliberately throws away any ptmiss information it is fed.  Why then do we construct a variable, only to throw it away later?  Really this is so that we can piggy back on the implementation of the basic_MtGen_332_Calculator to avoid duplication of code.
      static const TwoVector ptMiss_ignored;
      const double ans = m_basic_MtGen_332_Calculator->mtGen_332(interestingParticles, ptMiss_ignored, mChi);
      m_lastSolutionType = m_basic_MtGen_332_Calculator->lastSolutionType();
      return ans;
    }
    // Want to get rid of this:
    virtual double mtGen_330(const std::vector<LorentzTransverseVector> & interestingParticles, const double mChi) {
      //warning! The next line LOOKS like it is an error (it looks like we oughtto be setting ptMiss to minus the pt sum of the interesting particles) but in actual fact there is no error here.  The key is to understand that this class uses an "Mt2ApproximatingAdapter_332_from_330" to calculate MT2 values, and this sort of caclulator deliberately throws away any ptmiss information it is fed.  Why then do we construct a variable, only to throw it away later?  Really this is so that we can piggy back on the implementation of the basic_MtGen_332_Calculator to avoid duplication of code.
      static const TwoVector ptMiss_ignored;
      const double ans = m_basic_MtGen_332_Calculator->mtGen_332(interestingParticles, ptMiss_ignored, mChi);
      m_lastSolutionType = m_basic_MtGen_332_Calculator->lastSolutionType();
      return ans;
    }
    Basic_MtGen_330_Calculator(Mt2_330_Calculator & mt2_330_calc) : 
      MtGen_330_Calculator("Basic_MtGen_330_Calculator"),
      m_mt2_332_Calculator(0),
      m_basic_MtGen_332_Calculator(0) {
      m_mt2_332_Calculator = new Mt2ApproximatingAdapter_332_from_330(mt2_330_calc);
      if (!m_mt2_332_Calculator) {
	throw __FILE__ ;
      }
      m_basic_MtGen_332_Calculator = new Basic_MtGen_332_Calculator(*m_mt2_332_Calculator);
      if (!m_basic_MtGen_332_Calculator) {
	delete m_mt2_332_Calculator;
	throw __FILE__ ;
      }
    }
  public:
    void setDebug(const bool val) {
      (static_cast<Mt2Calculator*>(this))->setDebug(val);
      m_basic_MtGen_332_Calculator->setDebug(val);
      m_mt2_332_Calculator->setDebug(val);
    }
  public:
    virtual ~Basic_MtGen_330_Calculator() {
      delete m_basic_MtGen_332_Calculator;
      delete m_mt2_332_Calculator;
    }
  private:
    Mt2_332_Calculator         *         m_mt2_332_Calculator;
    Basic_MtGen_332_Calculator * m_basic_MtGen_332_Calculator;
  };

} //end namespace Mt2

#endif
