// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef SOLUTIONTYPE_H
#define SOLUTIONTYPE_H

#include <iostream>

namespace Mt2 {

  typedef enum { 
    ASideGlobalMin=0,
    BSideGlobalMin=1,
    BalancedSolution=2,
    NotSpecified=3
  } SolutionTypeEnum;
  
  struct SolutionType {
    SolutionType() : type(NotSpecified) {}
    explicit SolutionType(SolutionTypeEnum t) : type(t) {}
    SolutionTypeEnum type;
  };
  
}

inline std::ostream & operator<< (std::ostream & os, const Mt2::SolutionType & st) {
  os << st.type << "@";
  switch (st.type) {
  case Mt2::ASideGlobalMin:
    os << "ASide";
    break;
  case Mt2::BSideGlobalMin:
    os << "BSide";
    break;
  case Mt2::BalancedSolution:
    os << "Balanced";
    break;
  case Mt2::NotSpecified:
    os << "NS";
    break;
  default:
    os << "Error in " << __FILE__ << " at line " << __LINE__;
    break;
  }
  return os;
}

#endif
