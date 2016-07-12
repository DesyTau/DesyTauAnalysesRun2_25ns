// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_MINIMISER_H
#define MT2_MINIMISER_H

namespace Mt2 {
  class Mt2Calculator;
  /*
    Abstract base class defines interface for concrete minimisers.
    The job of this class is to take a Mt2Calculator, and 
    to minimise MT2 over sqrt(s) to within tolerance of 
    its minimum.
  */
  class Mt2Minimiser{
  public:
    /// virtual destructor
    virtual ~Mt2Minimiser(){};
    /// set the maximum number of evaluations of the function allowed for 
    virtual void setMaxEvaluations(int n) = 0;
    /// set the tolerance (relative uncertainty in the function)
    /// default value 1E-8
    virtual void setTolerance(double tol) = 0;
    /// do the minimisation
    virtual double minimise(const Mt2Calculator& calc) = 0;
  };
}

#endif
