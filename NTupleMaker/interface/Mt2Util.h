// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2_UTIL_H
#define MT2_UTIL_H

#include <iostream>
#include <exception>
#include <cmath>

/** 
    Tools for calculation of the transverse mass.
    If you use this package please cite: 
    (1) C.G.Lester, D.J.Summers; Phys.Lett.B.463:99-103 (1999) hep-ph/9906349
    (2) A.J.Barr, C.G.Lester, P.Stephens; J.Phys.G 29:2343-2363 (2003) hep-ph/0304226    

   @author Alan Barr
   @date 9 Feb 2006
*/
namespace Mt2 {

/// The following should be a simple inline function, but CINT cannot cope with utility function in a namespace it seems
struct Util {
  static inline double square(const double x) {
    return x*x;
  }
};

/** internal class to allow construction from CLHEP without library dependency */
struct CLHEP_HepLorentzVector_Type_Specifier {};


/**
   internal class defines an exception which may be thrown if 

   @author Alan Barr
   @date 9 Feb 2006
*/
class Mt2Exception : public std::exception {
 public:
  /** constructor with reason. does not throw */
  Mt2Exception(const std::string & reason) throw();
  /** virtual constructor with reason. does not throw */
  virtual ~Mt2Exception() throw();
  /** override std::exception what() method */
  virtual const char* what() const throw ();
 private:
  std::string m_reason;
};

/// An exception class to indicate that the exception was caused by the user's stupidity, not a failure in the program:
class UserWasStupidException : public Mt2Exception {
public:
  UserWasStupidException(const std::string & reason) :
    Mt2Exception(reason) {}
};

/// An exception class to indicate that the program failed to find a value of MT2 due to its own incompetence, rather than because a user abused it:
class ProgramWasNotCleverEnoughException : public Mt2Exception {
public:
  ProgramWasNotCleverEnoughException(const std::string & reason) :
    Mt2Exception(reason) {}
};

class EventIsIncompatibleWithRootS_Exception : public UserWasStupidException {
public:
  EventIsIncompatibleWithRootS_Exception() :
    UserWasStupidException("This event is not compatible with stated maximum value of sqrt(s), and so no mT2 value can be assigned to it.") {}
};
class EventIsAlmostIncompatibleWithRootS_Exception : public UserWasStupidException {
public:
  EventIsAlmostIncompatibleWithRootS_Exception() :
    UserWasStupidException("The kinematics of the event are so close to the limit of being incompatible with the stated maximum value of sqrt(s) that it is not safe or meaningful to calculate mT2 for this event with that value of sqrt(s).") {}
};

}


#endif
