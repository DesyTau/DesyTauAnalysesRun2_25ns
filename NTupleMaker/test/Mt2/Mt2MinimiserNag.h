// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_MINIMISERNAG_H
#define MT2_MINIMISERNAG_H
#include "Mt2Minimiser.h"
#include "nag.h"

namespace Mt2 {
  /**
     minimise MT2 using the nag c-library function 
     nag_opt_one_var_no_deriv (e04abc)
     http://www.nag.co.uk/numeric/CL/nagdoc_cl08/xhtml/E04/e04abc.xml
  */
  class Mt2MinimiserNag : public Mt2Minimiser {
  public:
    /// constructor sets default tolerance
    Mt2MinimiserNag();
    /// constructor sets default tolerance
    virtual ~Mt2MinimiserNag();
    /// over-ride base class function
    virtual double minimise(const Mt2Calculator& calc);
    /// static function for nag to use. Calulator is pointer to by the comm struct.
    static void static_funct(double x, double* f, Nag_Comm* comm);

    /** 
	set the machine tolerance. Default it at 1.E-8
    */
    virtual void setTolerance(double tol);
    /**
       set the maximum number of evaluations allowed in trying to find the min
     */
    virtual void setMaxEvaluations(int i);
  private:
    double m_tolerance;
    int m_max_fun;
  };
}

#endif
