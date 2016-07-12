// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Analytic_Mt2_332_Assistor_H
#define Analytic_Mt2_332_Assistor_H

// Sometimes the momentum configuration that leads to the minimised Mt2 value is an "unbalanced" configuration in which MTa != MTb.  Those configurations can always be determined analytically, even if the balanced solutions are not always calculable analytically.

// The purpose of this file is to answer the question: "Does the minimised Mt2 value occur at an (easy!) unbalanced solution, and if it does, then what is that Mt2 value?"

#include "Mt2/Mt2Vectors.h"
#include "Mt2/SolutionType.h"

namespace Mt2 {

  class Analytic_Mt2_332_Assistor {
    bool m_debug;
  public:
    Analytic_Mt2_332_Assistor() : m_debug(false) {}
  public:
    bool debug() const { return m_debug; }
  public:
    // Is it worth integrating the following struct 
    // with Analytic_Mt2_330_Calculator::Notes?
    struct Ans {
      SolutionType solutionType;
      double       mt2;
      TwoVector    ptInvisP;
      TwoVector    ptInvisQ;
    };
  public:
    Ans assist(const LorentzTransverseVector& visA, // 3 d.o.f. 
	       const LorentzTransverseVector& visB, // 3 d.o.f.
	       const TwoVector& ptMiss,             // 2 d.o.f.
	       const double mChi) {

      Ans ans;

      const double mASq = visA.masssq();
      const double mBSq = visB.masssq();
      const double mA = sqrt(mASq);
      const double mB = sqrt(mBSq);

      const double mtgma = mA + mChi; // mass for A side if soln is an unconstrained min
      const double mtgmb = mB + mChi; // mass for B side if soln is an unconstrained min
      
#ifndef __CINT__
#warning .... should deal with mA = 0 (massless cases) separately:
#endif
      //std::cout << "mChi " << mChi << " mA " << mA << " mB " << mB << std::endl;
      const TwoVector pIfASideGlobMin = visA.transverse()*mChi/mA;
      const TwoVector qIfBSideGlobMin = visB.transverse()*mChi/mB;
      
      const LorentzTransverseVector sigma = visA + visB;
      //const LorentzTransverseVector delta = visA - visB;

      const TwoVector qIfASideGlobMin = -pIfASideGlobMin + ptMiss; // OLD BUG: - (sigma.transverse()); -- NB "Sigma" (if we had it) rather than "sigma" would have been correct
      const TwoVector pIfBSideGlobMin = -qIfBSideGlobMin + ptMiss; // OLD BUG: - (sigma.transverse()); -- NB "Sigma" (if we had it) rather than "sigma" would have been correct
      
      const LorentzTransverseVector sideBCorrspToASideGlobMin = visB + LorentzTransverseVector(qIfASideGlobMin, mChi);
      const LorentzTransverseVector sideACorrspToBSideGlobMin = visA + LorentzTransverseVector(pIfBSideGlobMin, mChi);
      
      if (mtgma*mtgma >= sideBCorrspToASideGlobMin.masssq()) {
	if(debug()) {std::cout << "Found unbalanced solution at A side unconstrained min." << std::endl;}
	ans.solutionType = SolutionType(ASideGlobalMin);
	ans.mt2 = mtgma;
	ans.ptInvisP = pIfASideGlobMin;
	ans.ptInvisQ = qIfASideGlobMin;
	return ans;
      }

      if (mtgmb*mtgmb >= sideACorrspToBSideGlobMin.masssq()) {
        if(debug()){std::cout << "Found unbalanced solution at B side unconstrained min." << std::endl;}
	ans.solutionType = SolutionType(BSideGlobalMin);
	ans.mt2 = mtgmb;
	ans.ptInvisP = pIfBSideGlobMin;
	ans.ptInvisQ = qIfBSideGlobMin;
	return ans;
      }

      ans.solutionType = SolutionType(NotSpecified);
      return ans;
    }
  };

} //end namespace Mt2

#endif
