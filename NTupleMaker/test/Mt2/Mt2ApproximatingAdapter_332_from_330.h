// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2APPROXIMATINGADAPTER_332_FOMR_330_H
#define MT2APPROXIMATINGADAPTER_332_FOMR_330_H

#include "Mt2_330_Calculator.h"
#include "Mt2_332_Calculator.h"

#include <sstream>

namespace Mt2 {

  // This class is used to PRETEND that a 330 mt2 alg is actually a 332 alg.
  // The interface that pretends to be 332, just throws away the information in ptmiss.
  // This approximation is valid in the no ISR case only.

  class Mt2ApproximatingAdapter_332_from_330 : public Mt2_332_Calculator {
  private:
    static std::string incorporateName(const std::string & n) {
      std::ostringstream os;
      os << "Mt2ApproximatingAdapter_332_from_330_using_" << n;
      return os.str();
    }
  public:
    Mt2ApproximatingAdapter_332_from_330(Mt2_330_Calculator & calc_330) :
      Mt2_332_Calculator(incorporateName(calc_330.algorithmName())), 
      m_calc_330(calc_330) {
    }
  public:
    virtual double mt2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			   const LorentzTransverseVector& visibleB, // 3 d.o.f.
			   const TwoVector& ptmiss,                 // 2 d.o.f. -- will be thrown away!
			   double mInvisible) {
      
      const double ans = m_calc_330.mt2_330(visibleA, visibleB, mInvisible);
      m_lastSolutionType = m_calc_330.lastSolutionType();
      return ans;
    }
    virtual double mt2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
			      const LorentzTransverseVector& visibleB, // 3 d.o.f.
			      const TwoVector& ptmiss,                 // 2 d.o.f. -- will be thrown away!
			      double mInvisible) {
      
      const double ans = m_calc_330.mt2_330_Sq(visibleA, visibleB, mInvisible);
      m_lastSolutionType = m_calc_330.lastSolutionType();
      return ans;
    }
  private:
    Mt2_330_Calculator & m_calc_330;
  };

}

#endif
