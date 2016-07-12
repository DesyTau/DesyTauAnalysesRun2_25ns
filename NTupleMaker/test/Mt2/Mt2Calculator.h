// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_CALCULATOR_H
#define MT2_CALCULATOR_H

#include "SolutionType.h"
#include <string>

namespace Mt2 {
/**
   Class which knows how to calculate M_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2Calculator {
  private:
    static bool m_initialised;
    bool m_debug;
    std::string m_algName;
  public:
    Mt2Calculator(const std::string & algName);
    void setDebug(const bool val) { m_debug=val; }
    bool debug() const { return m_debug; }
    const std::string & algorithmName() const { return m_algName; }
    virtual ~Mt2Calculator() {}
    SolutionType lastSolutionType() const { return m_lastSolutionType; }
  protected:
    SolutionType m_lastSolutionType;
  };

}

#endif
