// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_MINIMISERGSL_H
#define MT2_MINIMISERGSL_H
#include "Mt2Minimiser.h"

namespace Mt2 {
  /**
     minimise MT2 using the GSL (gnu scientific library)
     one dimensional minimisation
     routines.  See http://www.gnu.org/software/gsl
  */
  class Mt2MinimiserGsl : public Mt2Minimiser {
  public:
    /// constructor sets default tolerance
    Mt2MinimiserGsl();
    /// constructor sets default tolerance
    virtual ~Mt2MinimiserGsl();
    /// over-ride base class function
    virtual double minimise(const Mt2Calculator& calc);
    /// static function for gsl to use. Calulator is pointer to by the comm struct.
    static double static_funct(double x, void* params);

    /** 
	set the machine tolerance. Default it at 1.E-8
    */
    virtual void setTolerance(double tol);
    /**
       set the maximum number of evaluations allowed in trying to find the min
     */
    virtual void setMaxEvaluations(int i);
  private:
    struct Params { const Mt2Calculator * calc; };
    double m_tolerance;
    int m_max_fun;
    static const int m_debug = 0;
  };
}

#endif
