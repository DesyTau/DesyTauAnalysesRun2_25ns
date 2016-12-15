#ifndef LESTER_MTTRUE_EXPORT_H
#define LESTER_MTTRUE_EXPORT_H

/**

  Copyright Christopher Lester.   Please contact Alan Barr or Christopher Lester if you wish to incorporate it in your code.  We will grant permission (!) but just want to have an idea of whether the code is un use.

  Kinematic variable providing lower bound on the mass of a resonance R undergoing the decay decay:
	
	R -> S T I1 I2

  where S and T are visible particles
  and I1 and I2 are invisible particles which are only constrained by knowledge of the sum of their masses and by their total transverse momentum (ptmiss).

  The mass of particle I1 plus the mass of particle I2 is denoted by the parameter "mmass".
  Note this is not the same as the invariant (or transverse) mass of the pair of invisible paricles!

  If you use this code, please cite the paper where we introduced this variable in the context of H -> WW -> lnulnu decays:

	http://arxiv.org/abs/0902.4864

  or the "squirrel paper" in which it is described as M1T (eq 111)

	http://arxiv.org/abs/1105.2977
	
   Use the mTTrue method defined below.
*/

#include <cmath>

namespace Lester {

  double mTTrue(const double se, const double sx, const double sy, const double sz,
                const double te, const double tx, const double ty, const double tz,
                const double pmissx, const double pmissy,
		const double mmass) {
    
    const double ve = se + te;
    const double vx = sx + tx;
    const double vy = sy + ty;
    const double vz = sz + tz;
    
    const double mVisSq = fabs(ve*ve - vx*vx - vy*vy - vz*vz);
    const double mInvisSq = mmass*mmass;
    
    const double etVis = sqrt(mVisSq + vx*vx + vy*vy);
    const double etInvis = sqrt(mInvisSq + pmissx*pmissx + pmissy*pmissy);

    return sqrt(
		fabs(
		     mVisSq + 
		     mInvisSq +
		     2.*(
			 etVis*etInvis - vx*pmissx - vy*pmissy
			 )
		     )
		);
  }
  
}

#endif 


