// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef mT2Fcn_4441_a_H
#define mT2Fcn_4441_a_H

#include "Minuit2/FCNBase.h"
#include "Mt2/Mt2Vectors.h"
#include <vector>

namespace Mt2 {

class mT2Fcn_4441_a : public ROOT::Minuit2::FCNBase {

 public:

  mT2Fcn_4441_a(const Mt2::LorentzVector& pxpypzeAlpha_in,      // 4 d.o.f. 
		const Mt2::LorentzVector& pxpypzeBeta_in,       // 4 d.o.f. 
		const Mt2::LorentzVector& pxpypzeG_in,          // 4 d.o.f.
		const double rootS_in,    /* eg 14 TeV */     // 1 d.o.f
		const double mChi_in,
		const double sign_in,
		const bool debug = false) : 
    theErrorDef(1.),
    pxpypzeAlpha(pxpypzeAlpha_in),     
    pxpypzeBeta(pxpypzeBeta_in),       
    pxpypzeG(pxpypzeG_in),         
    rootS(rootS_in),   
    mChi(mChi_in),
    sign(sign_in),
    m_debug(debug) {
  }

  ~mT2Fcn_4441_a() {}

  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const std::vector<double>&) const;

  void setErrorDef(double def) {theErrorDef = def;}
  
private:
  static double trip(const double a,
		     const double b,
		     const double c);
  static Mt2::LorentzVector boostFirstToRFOfSecond(const Mt2::LorentzVector & a,
						 const Mt2::LorentzVector & b);
private:
  double theErrorDef;
  const Mt2::LorentzVector& pxpypzeAlpha;     
  const Mt2::LorentzVector& pxpypzeBeta;       
  const Mt2::LorentzVector& pxpypzeG;         
  const double rootS;
  const double mChi;
  const double sign;
  bool m_debug;
};

}

#endif //mT2Fcn_4441_a_H
  
