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

#ifndef MT2_LORENTZ_TRANSVERSE_VECTOR_H
#define MT2_LORENTZ_TRANSVERSE_VECTOR_H

//fwd dec
namespace Mt2 {
	class LorentzTransverseVector;
}

#include <iostream>
#include <exception>
#include <cmath>

#include "Mt2TwoVector.h"

namespace Mt2 {

/**
   A lorentz vector with 2 spatial and 1 time dimensions
   Has data members px, py, mass, and calculates
   Et, pt and mass on request.
   Compiler-default copy, assignment operators and destructors are used.

   @author Alan Barr
   @date 9 Feb 2006
*/

//fwd dec:
class LorentzVector;

class LorentzTransverseVector {
 public:
  /**
     vector with all px, py and mass set to zero.
   */
  LorentzTransverseVector();
  /** 
      constructor for a massless particle;
  */
  LorentzTransverseVector(double px, double py);
  /** 
      constructor for a massive particle;
  */  
  LorentzTransverseVector(double et, double px, double py);
  /** 
      constructor for a massless particle
  */
  LorentzTransverseVector(const TwoVector& momentum);
  /** 
      constructor for a massive particle;
  */  
  LorentzTransverseVector(double Et, const TwoVector& momentum);
  /** 
      constructor for a particle, from momentum and mass
  */
  LorentzTransverseVector(const TwoVector& momentum, double mass);
  /** 
      constructor from a lorentz vector 
  */
  explicit LorentzTransverseVector(const LorentzVector & lv);

  /** return the x compoment */
  virtual double px() const;
  /** return the y compoment */
  virtual double py() const;
  /** return the transeverse momentum */
  virtual double pt() const;
  /** return the transeverse momentum squared*/
  virtual double ptsq() const;
  /** return of the transverse energy */
  virtual double Et() const;
  /** return of the transverse energy squared*/
  virtual double Etsq() const;
  /** impliment return the mass */
  virtual double mass() const;  
  /** impliment return the mass squared */
  virtual double masssq() const;  
  /** dot product operator */
  virtual double dot(const LorentzTransverseVector& second) const;
  /** dot product operator invariant under contralinear boosts (used by MCT) */
  virtual double contralinearDot(const LorentzTransverseVector& second) const;
  /** get the spatial TwoVector */
  virtual TwoVector transverse() const;
  /** binary addition operator */
  virtual LorentzTransverseVector operator+(const LorentzTransverseVector& second) const;
  /** += operator */
  virtual void operator+=(const LorentzTransverseVector& second);
  /** unary minus operator */
  //  virtual LorentzTransverseVector operator-() const;
  /** binary subtraction operator */
  virtual LorentzTransverseVector operator-(const LorentzTransverseVector& second) const;
  /** multiplication by a factor */
  virtual LorentzTransverseVector operator*(double factor) const;
  /** division by a factor */
  virtual LorentzTransverseVector operator/(double factor) const;
  /** print function */
  void print(std::ostream& os) const;
  /** Set from class implementing blah.m(), blah.x(), blah.y() accessors, eg CLHEP LorentzVectors */
  template <class MXY>
  void setFromMXY(const MXY & other) {
    m_px = other.x();
    m_py = other.y();
    const double m = other.m();
    m_et = sqrt(m_px*m_px + m_py*m_py + m*m); 
  }
 protected:
  double m_px; // x-component
  double m_py; // y-component
  double m_et; // mass component
};

inline LorentzTransverseVector::LorentzTransverseVector(double px, double py)
  : m_px(px), m_py(py)
{
  m_et = transverse().pt();
}

inline LorentzTransverseVector::LorentzTransverseVector()
  : m_px(0.), m_py(0.), m_et(0.)
{
}

inline LorentzTransverseVector::LorentzTransverseVector(const TwoVector& v)
  : m_px(v.px()), m_py(v.py()), m_et(v.pt())
{
}

inline LorentzTransverseVector::LorentzTransverseVector(double et, double px, double py)
  : m_px(px), m_py(py), m_et(et)
{
}

inline LorentzTransverseVector::LorentzTransverseVector(double et, const TwoVector& v)
  : m_px(v.px()), m_py(v.py()), m_et(et)
{
}

inline LorentzTransverseVector::LorentzTransverseVector(const TwoVector& v, double mass)
  : m_px(v.px()), m_py(v.py())
{
  m_et = sqrt( Util::square(mass) + v.ptsq() );
}


inline double LorentzTransverseVector::px() const {
  return m_px;
}

inline double LorentzTransverseVector::py() const {
  return m_py;
}

inline TwoVector LorentzTransverseVector::transverse() const {
  return TwoVector(px(), py());
}

inline double LorentzTransverseVector::mass() const {
  double mass_sq = masssq();
  if (mass_sq<0) return -sqrt(-mass_sq);
  return sqrt(mass_sq);
}

inline double LorentzTransverseVector::masssq() const {
  return Util::square(m_et) - ptsq();
}

inline double LorentzTransverseVector::pt() const {
  return sqrt( ptsq() );
}

inline double LorentzTransverseVector::ptsq() const {
  return Util::square(px()) + Util::square(py());
}

inline double LorentzTransverseVector::Et() const {
  return m_et;
}

inline double LorentzTransverseVector::Etsq() const {
  return Util::square(m_et);
}

inline LorentzTransverseVector LorentzTransverseVector::operator*(double factor) const {  
  return LorentzTransverseVector(Et()*factor, px()*factor, py()*factor);
}

inline LorentzTransverseVector LorentzTransverseVector::operator/(double factor) const {  
  if (factor==0.) throw Mt2Exception("divide by zero error");
  return (*this) * (1. / factor);
}

inline void LorentzTransverseVector::operator+=(const LorentzTransverseVector& second) {  
  m_px += second.m_px;
  m_py += second.m_py;
  m_et += second.m_et;
}

inline LorentzTransverseVector LorentzTransverseVector::operator+(const LorentzTransverseVector& second) const {  
  double pxtot = px() + second.px();
  double pytot = py() + second.py();
  double etot = Et() + second.Et();
  return LorentzTransverseVector(etot, pxtot, pytot);
}

inline LorentzTransverseVector LorentzTransverseVector::operator-(const LorentzTransverseVector& second) const {  
  double pxtot = px() - second.px();
  double pytot = py() - second.py();
  double etot = Et() - second.Et();
  return LorentzTransverseVector(etot, pxtot, pytot);
}

inline double LorentzTransverseVector::dot(const LorentzTransverseVector& second) const {  
  return Et()*second.Et() - px()*second.px() - py()*second.py();
}
inline double LorentzTransverseVector::contralinearDot(const LorentzTransverseVector& second) const {  
  return Et()*second.Et() + px()*second.px() + py()*second.py();
}

/// The following follows the definitions of 8, 9 in http://arxiv.org/pdf/0910.1584v1
  class ResolvedLTV {
  public:
    ResolvedLTV(const LorentzTransverseVector & a,
		const TwoVector & u) {
      const TwoVector & aTwoVec = a.transverse();
      const TwoVector paraBit = u*(aTwoVec.dot(u)/u.ptsq());
      const TwoVector perpBit = aTwoVec - paraBit;
      m_para = LorentzTransverseVector(paraBit, a.mass());
      m_perp = LorentzTransverseVector(perpBit, a.mass());
    } 
    const LorentzTransverseVector & para() const {
      return m_para;
    }
    const LorentzTransverseVector & perp() const {
      return m_perp;
    }
  private:
    LorentzTransverseVector m_para;
    LorentzTransverseVector m_perp;
  };


}

#endif
