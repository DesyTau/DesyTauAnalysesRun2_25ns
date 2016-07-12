// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#ifndef SUSYPhys_Mt2_222_Calculator_H
#define SUSYPhys_Mt2_222_Calculator_H

#include "Mt2/Mt2Calculators.h"
#include "Minuit2/FCNBase.h"
#include <vector>

namespace Mt2 {

  class SUSYPhys_Mt2_222_Calculator : public Mt2_222_Calculator {
  public:
    double mt2_222(const TwoVector& visibleA, // 2 d.o.f. 
		   const TwoVector& visibleB, // 2 d.o.f.
		   const TwoVector& ptmiss,   // 2 d.o.f.
		   const double mInvisible);
    SUSYPhys_Mt2_222_Calculator() : Mt2_222_Calculator("SUSYPhys_Mt2_222") {};
  private:
    class mT2Fcn : public ROOT::Minuit2::FCNBase {
    public:
      mT2Fcn(const double exmiss, 
	     const double eymiss, 
	     const double mchi,
	     const double pT1x,
	     const double pT1y,
	     const double pT2x,
	     const double pT2y) : 
	theExmiss(exmiss),
	theEymiss(eymiss),
	theMchi(mchi),
	thePT1x(pT1x),
	thePT1y(pT1y),
	thePT2x(pT2x),
	thePT2y(pT2y),
	theErrorDef(1.) {}
      
      ~mT2Fcn() {}
      
      virtual double Up() const {return theErrorDef;}
      virtual double operator()(const std::vector<double>&) const;
      
      const double exMiss() const {return theExmiss;}
      const double eyMiss() const {return theEymiss;}
      const double mChi() const {return theMchi;}
      const double pT1x() const {return thePT1x;}
      const double pT1y() const {return thePT1y;}
      const double pT2x() const {return thePT2x;}
      const double pT2y() const {return thePT2y;}
      
      void setErrorDef(double def) {theErrorDef = def;}
      
    private:
      
      double theExmiss;
      double theEymiss;
      double theMchi;
      
      double thePT1x;
      double thePT1y;
      
      double thePT2x;
      double thePT2y;
      
      double theErrorDef;
    };
    
  };  // end class SUSYPhys_Mt2_222_Calculator

} // end Mt2 namespace

#endif
