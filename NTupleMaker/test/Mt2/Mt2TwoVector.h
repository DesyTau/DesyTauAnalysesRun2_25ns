// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_TWO_VECTORS_H
#define MT2_TWO_VECTORS_H

#include <iostream>
#include <exception>
#include <cmath>

#include "Mt2Util.h"

/** 
    Tools for calculation of the transverse mass.
    If you use this package please cite: 
    (1) C.G.Lester, D.J.Summers; Phys.Lett.B.463:99-103 (1999) hep-ph/9906349
    (2) A.J.Barr, C.G.Lester, P.Stephens; J.Phys.G 29:2343-2363 (2003) hep-ph/0304226    

   @author Alan Barr
   @date 9 Feb 2006
*/
namespace Mt2 {

class TwoVector{
 public:
  /** constructor initialised with x and y components */
  TwoVector(double px, double py);
  /** constructor initialised with x and y components = 0*/
  TwoVector();
  /** return x component */
  double px() const;
  /** return y component */
  double py() const;
  /** calculate magnitude */
  double pt() const;
  /** calculate magnitude squared */
  double ptsq() const;
  /** print function */
  void print(std::ostream& os) const;
  /** dot product with another TwoVector */
  double dot(const TwoVector& other) const;
  /** unary negative operator*/
  TwoVector operator-() const;
  /** binary negative operator*/
  TwoVector operator-(const TwoVector& other) const;
  /** binary self addition operator*/
  void operator+=(const TwoVector& other) ;
  /** binary self subtraction operator*/
  void operator-=(const TwoVector& other) ;
  /** binary sum operator*/
  TwoVector operator+(const TwoVector& other) const;
  /** multiplication by scale factor */
  TwoVector operator*(const double factor) const;
  /** division by scale factor */
  TwoVector operator/(const double factor) const;
  /** set contents from any class implementing blah.x() and blah.y() accessors -- eg CLHEP LorentzVectors */ 
  template<class XY>
  void setFromXY(const XY & other) { m_px = other.x(); m_py = other.y(); }
 protected:
  double m_px;
  double m_py;
};

inline TwoVector::TwoVector(double px, double py)
  : m_px (px), m_py(py) {
}

inline TwoVector::TwoVector()
  : m_px (0.), m_py(0.) {
}

inline double TwoVector::px() const {
  return m_px;
} 

inline double TwoVector::py() const {
  return m_py;
} 

inline double TwoVector::pt() const {
  return sqrt( ptsq() );
}

inline double TwoVector::ptsq() const {
  return Util::square(px()) +  Util::square(py());
}

inline double TwoVector::dot(const TwoVector& other) const{
  return px()*other.px() + py()*other.py();
}

inline TwoVector TwoVector::operator-(const TwoVector& other) const {
  return (*this) + (-other);
}

inline TwoVector TwoVector::operator+(const TwoVector& other) const {
  return TwoVector(px() + other.px(), 
		   py() + other.py());
}

inline TwoVector TwoVector::operator-() const {
  return TwoVector(-px(), -py());
}

inline void TwoVector::operator+=(const TwoVector& other) {
  m_px += other.m_px;
  m_py += other.m_py;
}

inline void TwoVector::operator-=(const TwoVector& other) {
  m_px -= other.m_px;
  m_py -= other.m_py;
}

inline TwoVector TwoVector::operator*(const double f) const {
  return TwoVector(px()*f, py()*f);
}

inline TwoVector TwoVector::operator/(const double f) const {
  return TwoVector(px()/f, py()/f);
}

}

#endif
