// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

/** 
    Tools for calculation of the transverse mass.
    If you use this package please cite: 
    (1) C.G.Lester, D.J.Summers; Phys.Lett.B.463:99-103 (1999) hep-ph/9906349
    (2) A.J.Barr, C.G.Lester, P.Stephens; J.Phys.G 29:2343-2363 (2003) hep-ph/0304226    

   @author Alan Barr
   @date 9 Feb 2006
*/

#ifndef MT2_LORENTZ_VECTOR_H
#define MT2_LORENTZ_VECTOR_H

// fwd dec:
namespace Mt2 {
	class LorentzVector;
}

#include <iostream>
#include <exception>
#include <cmath>

#include "Mt2LorentzTransverseVector.h"

namespace Mt2 {

class LorentzVector {
public:
	void setEPxPyPz(const double te,
			const double tpx,
			const double tpy,
			const double tpz) {
		e  = te ;
		px = tpx;
		py = tpy;
		pz = tpz;
	}
	void setVectM(const double tpx,
                        const double tpy,
                        const double tpz,
			const double m) {
		px=tpx;
		py=tpy;
		pz=tpz;
		e=sqrt(m*m + p2());
	}
	LorentzVector operator+(const LorentzVector & other) const {
		return LorentzVector(	e +other.e ,
					px+other.px,
					py+other.py,
					pz+other.pz);
	}
	void operator+=(const LorentzVector & other) {
	  e +=other.e;
	  px+=other.px;
	  py+=other.py;
	  pz+=other.pz;
	}
	LorentzVector operator-(const LorentzVector & other) const {
		return LorentzVector(	e -other.e ,
					px-other.px,
					py-other.py,
					pz-other.pz);
	}
	LorentzVector operator*(const double d) const {
		return LorentzVector(d*e,d*px,d*py,d*pz);
	}
	double masssq() const {
		return m2();
	}
	double m2() const {
		return e*e - p2();
	}
	double m() const {
	  const double mm = m2();
	  if (mm<0) {
	    return -sqrt(-mm);
	  } else {
	    return sqrt(mm);
	  }
	}
	double ET() const {
	  return sqrt(ET2());
	}
	double ET2() const {
	  return e*e - pz*pz;
	}
	double p2() const {
		return px*px + py*py + pz*pz;
	}
	double pT2() const {
		return px*px + py*py;
	}
	double cosineOfSpatialSeparationAngleFrom(const LorentzVector & other) const {
		return (px*other.px + py*other.py + pz*other.pz)/sqrt(p2()*other.p2());
	}
        /// Lorentz invariant 4-dot product
        double dot(const LorentzVector & other) const {
	  return e*other.e - (px*other.px + py*other.py + pz*other.pz);
        }
        /// Lorentz contralinear boost invariant 4-dot product, used by MCT and MC etc
        double contralinearDot(const LorentzVector & other) const {
	  return e*other.e + (px*other.px + py*other.py + pz*other.pz);
        }


	LorentzTransverseVector getLorentzTransverseVector() const {
	  return LorentzTransverseVector(TwoVector(px,py), m());
	}

	TwoVector transverse() const {
          return TwoVector(px,py);
        }

	// boost by a beta vector:
	LorentzVector boostBy(double bx, double by, double bz) const {
   double b2 = bx*bx + by*by + bz*bz;
   double gamma = 1.0 / sqrt(1.0 - b2);
   double bp = bx*px + by*py + bz*pz;
   double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
 
	LorentzVector ans; // Order of the lines below is important!
   ans.px=(px + gamma2*bp*bx + gamma*bx*e);
   ans.py=(py + gamma2*bp*by + gamma*by*e);
   ans.pz=(pz + gamma2*bp*bz + gamma*bz*e);
   ans.e=(gamma*(e + bp));
   return ans;
	}
public:
	double e;
	double px;
	double py;
	double pz;
	struct InitEPxPyPz {
		InitEPxPyPz(const double e, const double px, const double py, const double pz) : e(e), px(px), py(py), pz(pz) {}
                template <class T> InitEPxPyPz(T t) : e(t.e()), px(t.px()), py(t.py), pz(t.pz()) {}
		double e;
        	double px;
        	double py;
        	double pz;
	};
public:
	LorentzVector() : e(0), px(0), py(0), pz(0) {}
	LorentzVector(const InitEPxPyPz & i) : e(i.e), px(i.px), py(i.py), pz(i.pz) {}
        /** Set from class implementing blah.e(), blah.px(), blah.py(), blah.pz() accessors, eg CLHEP LorentzVectors */
        template <class EPxPyPz>
        void setFromEPxPyPz(const EPxPyPz & other) {
    e  = other.e ();
    px = other.px();
    py = other.py();
    pz = other.pz();
  }

private:
	/// This constructor is private because the order of its parameters is unchecked.
	/// I would rather users called setEPxPyPz(...);
	LorentzVector(const double e, const double px, const double py, const double pz) : e(e), px(px), py(py), pz(pz) {}
};

}

inline Mt2::LorentzVector operator*(const double d, const Mt2::LorentzVector & v) {
	return v*d;
}

std::ostream& operator << (std::ostream& os , const Mt2::LorentzTransverseVector& v);
std::ostream& operator << (std::ostream& os , const Mt2::TwoVector& v);
inline std::ostream& operator << (std::ostream& os , const Mt2::LorentzVector& v) {
	return os << "Mt2::LorentzVector[e=" << v.e << ", px=" << v.px << ", py="<<v.py<<", pz="<<v.pz<<"]";
}

/** 
     constructor from a lorentz vector 
*/
inline Mt2::LorentzTransverseVector::LorentzTransverseVector(const Mt2::LorentzVector & lv) : 
m_px(lv.px), 
m_py(lv.py),
m_et(lv.ET()) {
} 


#endif
