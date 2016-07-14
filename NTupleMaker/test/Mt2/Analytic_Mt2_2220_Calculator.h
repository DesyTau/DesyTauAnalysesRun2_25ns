// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#ifndef Analytic_Mt2_2220_Calculator_H
#define Analytic_Mt2_2220_Calculator_H

#include "Mt2/Mt2Calculators.h"
#include <vector>
#include <complex>

namespace Mt2 {
  
  class Analytic_Mt2_2220_Calculator : public Mt2_2220_Calculator {
  public:
    double mt2_2220(const TwoVector& visibleA, // 2 d.o.f. 
		   const TwoVector& visibleB, // 2 d.o.f.
		   const TwoVector& ptmiss);
    double mt2_2220_Sq(const TwoVector& visibleA, // 2 d.o.f. 
		       const TwoVector& visibleB, // 2 d.o.f.
		       const TwoVector& ptmiss);
    typedef enum { FERRARI, R8POLY4 } RootFindingMethod;
    Analytic_Mt2_2220_Calculator(const bool debugMode=false, RootFindingMethod meth = R8POLY4) : Mt2_2220_Calculator("Analytic_Mt2_2220"), m_debugMode(debugMode), m_rootFindingMethod(meth) {};
  public:
    typedef enum enum_error_codes { no_mt2_value_found=-100 } EnumErrorCode;
    struct Info { /*info about last evaluation*/
      double Kss,Kcc,Ks,Kc,Kcs,K1;
      double A,B,C,D,E;
      std::complex<double> r1,r2,r3,r4; // quartic roots for sintheta -- only filled if you use the numerical root filling routine.
      Info() : Kss(0), Kcc(0), Ks(0), Kc(0), Kcs(0), K1(0),
	       A(0),B(0),C(0),D(0),E(0),
	       r1(0),r2(0),r3(0),r4(0) {
      }
    };
    Info info() const {
      return m_info;
    }
  private:
    void insertIfOK(const std::complex<double> &r1,
		    std::vector<double> & goodSinValues);
  private:
    double square(const double x) { return x*x; }
    double cube(const double x) { return x*x*x; }
    const bool m_debugMode;
    Info m_info; // info about last evaluation
    const RootFindingMethod m_rootFindingMethod;
  };  // end class Analytic_Mt2_222_Calculator

} // end Mt2 namespace

std::ostream & operator<<(std::ostream & os, const Mt2::Analytic_Mt2_2220_Calculator::Info & info) {
  return os 
    << " Kss= " << info.Kss
    << " Kcc= " << info.Kcc
    << " Kcs= " << info.Kcs
    << " Ks= " << info.Ks
    << " Kc= " << info.Kc
    << " K1= " << info.K1
    << " A= " << info.A
    << " B= " << info.B
    << " C= " << info.C
    << " D= " << info.D
    << " E= " << info.E
    << " r1= " << info.r1
    << " r2= " << info.r2
    << " r3= " << info.r3
    << " r4= " << info.r4;
}

#endif
