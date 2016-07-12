// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_Mt2_Asymmetric332_Calculator_H
#define Basic_Mt2_Asymmetric332_Calculator_H

#include "Mt2/Mt2_Asymmetric332_Calculator.h"
#include "Minuit2/FCNBase.h"
#include <vector>

namespace Mt2 {

  class Basic_Mt2_Asymmetric332_Calculator : public Mt2_Asymmetric332_Calculator {
  public:
    double mt2_Asymmetric332(const double theta,//tan(theta)=mHeavyA/mHeavyB
                   const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		   const LorentzTransverseVector& visibleB, // 3 d.o.f.
		   const TwoVector& ptmiss,                 // 2 d.o.f.
		   const double mInvisible);
    double mt2_Asymmetric332_Sq(const double theta,//tan(theta)=mHeavyA/mHeavyB
                       const LorentzTransverseVector& visibleA, // 3 d.o.f. 
		      const LorentzTransverseVector& visibleB, // 3 d.o.f.
		      const TwoVector& ptmiss,                 // 2 d.o.f.
		      const double mInvisible);
    Basic_Mt2_Asymmetric332_Calculator() : Mt2_Asymmetric332_Calculator("Basic_Mt2_Asymmetric332") {};
    
  private:
    class mT2Fcn : public ROOT::Minuit2::FCNBase {
    public:
      mT2Fcn(const double theta,//tan(theta)=mHeavyA/mHeavyB
             const double exmiss, 
	     const double eymiss, 
	     const double mchi,
	     const double pT1x,
	     const double pT1y,
	     const double mass1,
	     const double pT2x,
	     const double pT2y,
	     const double mass2) : 
	theTheta(theta),
	theMHeavyAOnMHeavyB(tan(theta)),
        theExmiss(exmiss),
	theEymiss(eymiss),
	theMchi(mchi),
	thePT1x(pT1x),
	thePT1y(pT1y),
	theMass1(mass1),
	thePT2x(pT2x),
	thePT2y(pT2y),
	theMass2(mass2),
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
      
      const double theTheta;
      const double theMHeavyAOnMHeavyB;

      double theExmiss;
      double theEymiss;
      double theMchi;
      
      double thePT1x;
      double thePT1y;
      double theMass1;
      
      double thePT2x;
      double thePT2y;
      double theMass2;
      
      double theErrorDef;
    };
  };

} //end namespace Mt2

#endif
