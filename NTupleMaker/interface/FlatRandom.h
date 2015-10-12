#ifndef LESTER_FLATRANDOM_H
#define LESTER_FLATRANDOM_H


/* Alter the implementation of the function FlatRandom() below so that
   you get your random numbers from your preferred source. */

//#include "CLHEP/Random/RandFlat.h"
#include "TMath.h"
#include <TRandom3.h>
#include "TROOT.h"

namespace Lester {

	double FlatRandom() {
		TRandom *rnd       = new TRandom();
		//return CLHEP::RandFlat::shoot();
		//return TRandom::Uniform(1);
		return rnd->Uniform();
        }

}

#endif




