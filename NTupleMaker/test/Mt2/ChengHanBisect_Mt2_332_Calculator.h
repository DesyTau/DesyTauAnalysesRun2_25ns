// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef ChengHanBisect_Mt2_332_Calculator_H
#define ChengHanBisect_Mt2_332_Calculator_H

#include "Mt2/Mt2_332_Calculator.h"
#include <vector>

/** Adapted from http://daneel.physics.ucdavis.edu/~zhenyuhan/mt2.html
    which is described in 
    Minimal Kinematic Constraints and MT2
    arXiv:0810.5178
    Hsin-Chia Cheng, Zhenyu Han 
*/

namespace ZhenyuHan {
  class mt2;
}

namespace Mt2 {

  class ChengHanBisect_Mt2_332_Calculator : public Mt2_332_Calculator {
  public:
    double mt2_332(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		   const LorentzTransverseVector& visibleB, // 3 d.o.f.
		   const TwoVector& ptmiss,                 // 2 d.o.f.
		   const double mInvisible);
    double mt2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		      const LorentzTransverseVector& visibleB, // 3 d.o.f.
		      const TwoVector& ptmiss,                 // 2 d.o.f.
		      const double mInvisible);
    ChengHanBisect_Mt2_332_Calculator();
    
  private:
    
    ZhenyuHan::mt2 * p_ZhenyuHan_mt2_pointer;
  };

} //end namespace Mt2

#endif
