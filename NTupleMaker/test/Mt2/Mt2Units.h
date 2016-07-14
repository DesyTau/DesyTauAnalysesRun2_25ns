// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef MT2UNITS_H
#define MT2UNITS_H

/* This header defines some VERY BASIC units for the Mt2 package.  The idea is to
   have the same units as in CLHEP but without needing am explicit dependancy on
   CLHEP */

namespace Mt2 {
	const double MeV = 1;
	const double GeV = 1000*MeV;
}

#endif
