// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Frugal_MtGen_330_Calculator_H
#define Frugal_MtGen_330_Calculator_H

#include "Mt2_330_Calculator.h"
#include "MtGen_330_Calculator.h"
#include "Frugal_MtGen_332_Calculator.h"
#include "Mt2ApproximatingAdapter_332_from_330.h"

namespace Mt2 {
/**
   Class which knows how to calculate a quantity like MT2 but without the missing particles and using all four dimensions.  Supposed to be equivalent to Basic_MtGen_330_Calculator, but using fewer operations.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Frugal_MtGen_330_Calculator : public MtGen_330_Calculator {
  public:
    virtual double mtGen_330(const std::vector<LorentzVector> & interestingParticles,
			     const double mChi) {
      //warning! The next line LOOKS like it is an error (it looks like we oughtto be setting ptMiss to minus the pt sum of the interesting particles) but in actual fact there is no error here.  The key is to understand that this class uses an "Mt2ApproximatingAdapter_332_from_330" to calculate MT2 values, and this sort of caclulator deliberately throws away any ptmiss information it is fed.  Why then do we construct a variable, only to throw it away later?  Really this is so that we can piggy back on the implementation of the frugal_MtGen_332_Calculator to avoid duplication of code.
      static const TwoVector ptMiss_ignored;
      const double ans = m_frugal_MtGen_332_Calculator->mtGen_332(interestingParticles, ptMiss_ignored, mChi);
      m_lastSolutionType = m_frugal_MtGen_332_Calculator->lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << m_lastSolutionType << std::endl;
      return ans;
    }
    virtual double mtGen_330(const std::vector<LorentzTransverseVector> & interestingParticles,
			     const double mChi) {
      //warning! The next line LOOKS like it is an error (it looks like we oughtto be setting ptMiss to minus the pt sum of the interesting particles) but in actual fact there is no error here.  The key is to understand that this class uses an "Mt2ApproximatingAdapter_332_from_330" to calculate MT2 values, and this sort of caclulator deliberately throws away any ptmiss information it is fed.  Why then do we construct a variable, only to throw it away later?  Really this is so that we can piggy back on the implementation of the frugal_MtGen_332_Calculator to avoid duplication of code.
      static const TwoVector ptMiss_ignored;
      const double ans = m_frugal_MtGen_332_Calculator->mtGen_332(interestingParticles, ptMiss_ignored, mChi);
      m_lastSolutionType = m_frugal_MtGen_332_Calculator->lastSolutionType();
      //std::cout << "NNN" << __FILE__ << " " << __LINE__ << m_lastSolutionType << std::endl;
      return ans;
    }
    Frugal_MtGen_330_Calculator(Mt2_330_Calculator & mt2_330_calc,
				const std::string & algName = "Frugal_MtGen_330_Calculator") : 
      MtGen_330_Calculator("Frugal_MtGen_330_Calculator"),
      m_mt2_332_Calculator(0),
      m_frugal_MtGen_332_Calculator(0) {
      m_mt2_332_Calculator = new Mt2ApproximatingAdapter_332_from_330(mt2_330_calc);
      if (!m_mt2_332_Calculator) {
	throw __FILE__ ;
      }
      m_frugal_MtGen_332_Calculator = new Frugal_MtGen_332_Calculator(*m_mt2_332_Calculator);
      if (!m_frugal_MtGen_332_Calculator) {
	delete m_mt2_332_Calculator;
	throw __FILE__ ;
      }
    }
  public:
    void setDebug(const bool val) {
      (static_cast<Mt2Calculator*>(this))->setDebug(val);
      m_frugal_MtGen_332_Calculator->setDebug(val);
      m_mt2_332_Calculator->setDebug(val);
    }
  public:
    virtual ~Frugal_MtGen_330_Calculator() {
      delete m_frugal_MtGen_332_Calculator;
      delete m_mt2_332_Calculator;
    }
  private:
    Mt2_332_Calculator          *          m_mt2_332_Calculator;
    Frugal_MtGen_332_Calculator * m_frugal_MtGen_332_Calculator;
  };
}

#endif
