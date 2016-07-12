// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_M2C_332_Calculator_H
#define Basic_M2C_332_Calculator_H

#include "Mt2/Mt2_332_Calculator.h"
#include <vector>

// M2C part of the Code written by Mario A. Serna Jr
// Last updated on 08 Feb 2008
// The variable M2C was introduced in
//   G.Ross and M.Serna, "Mass Determination of New States at Hadron Colliders", 
//   arXiv:0712.0943 [hep-ph] 06 Dec 2007

// M2C_332Calculator
//   Applies to new state decaying to light-stable state via a three body decay.
//   Needs the edge of the m12 distribution in the three-body decay to get a measure on m2-m1.
//   Now for each event, given the 2+1 momenta on branch A & B and the missing pT 
//   and the mass difference, this will find the smallest mass of the parent (second to lightest new state)
//   consistant with these states 
// Assumes momenta are all given in unist of GeV.

double M2C_332Calculator(Mt2::LorentzTransverseVector& visA, // 2+1 momenta from branch A  (eo, px, py)
					     Mt2::LorentzTransverseVector& visB, // 2+1 momenta from branch B  (eo, px, py)
					     Mt2::TwoVector& ptmiss, // the 2 momenta , the missing transverse momenta (px, py)
					     double deltam,  // The value of m2-m1 read off from the three-body decay m12 end point 
                                             Mt2::Mt2_332_Calculator * mt2CalculatorToUse=0
                    );


#endif
