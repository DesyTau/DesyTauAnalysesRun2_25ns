// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_Nt2_332_Calculator_H
#define Basic_Nt2_332_Calculator_H

#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/Mt2/Nt2_332_Calculator.h"
#include "/nfs/dust/cms/user/alkaloge/TauAnalysis/new/new/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/Mt2/Mt2Vectors.h"
#include "Minuit2/FCNBase.h"
#include <vector>
namespace Mt2 {


  double Basic_Nt2_332_Calculator::nt2_332(const LorentzTransverseVector& visA, 
					   const LorentzTransverseVector& visB,
					   const TwoVector& ptmiss, 
					   const double massParent1,
					   const double massParent2) {
    const double mSq = nt2_332_Sq(visA, visB, ptmiss, massParent1, massParent2);

    if (mSq>=0) {
      return sqrt(mSq);
    } else {
      // something went wrong -- maybe just rounding errors;
      return sqrt(fabs(mSq));
    }
  }


  double Basic_Nt2_332_Calculator::nt2_332_Sq(const LorentzTransverseVector& visA,
					      const LorentzTransverseVector& visB,
					      const TwoVector& ptmiss,
					      const double massParent1,
					      const double massParent2){


    mT2Fcn theFCN(ptmiss.px(),
		  ptmiss.py(),
                  massParent1,
		  massParent2,
		  visA.px(),
		  visA.py(),
		  visA.mass(),
		  visB.px(),
		  visB.py(),
		  visB.mass());

    const double massScale = (
			      ptmiss.pt() +
			      massParent1 +
			      massParent2 +
			      visA.pt()  +
			      visA.mass() +
			      visB.pt() +
			      visB.mass()
			      )/7.0;
    // DANG! Try to get rid of Minuit output:
    //std::ofstream    DANG_log("/dev/null");
    //std::streambuf * DANG_save = std::cerr.rdbuf();
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_log.rdbuf());
    // }

    double guessx = 0.5*(ptmiss.px()+visA.px()+visB.px());
    double guessy = 0.5*(ptmiss.py()+visA.py()+visB.py());
    
    MnUserParameters upar;
    upar.Add("etx", guessx, 0.02*massScale); 
    upar.Add("ety", guessy, 0.02*massScale);
    const int highQuality=2;    

    // Usually migrad produces the best minumum.
    // But when the minimum is in a fold, migrad can fail badly.
    // On the fold, simplex does well.  We therefore do both separately
    // and record the answer of the one that did best.
    
    // Further to the above notes, it now seems that by choosing the massScale sensibly, and by making the "tolerance" (the second argument to the call that extracts the FunctionMinimum below) much smaller, the simplex algorithm seems to work so well that we don't need the migrad algorithm.  This is good news, as the migrad algorithm produces lots of error output that we can't get rid of.

    MnSimplex simplex(theFCN, upar, highQuality);
    FunctionMinimum minS = simplex(0,massScale*0.000001);
    //const double etxAtMin = minS.UserState().Value("etx");
    //const double etyAtMin = minS.UserState().Value("ety");
    //MnMigrad migrad(theFCN, upar, highQuality);
    //FunctionMinimum minM = migrad(0,massScale*0.000001);
    //const double best = fmin(minS.Fval(), minM.Fval());
    const double best = minS.Fval();

    // DANG! Undoing our attempt to get rid of Minuit output:
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_save);
    //}

    return -best;    
  }
  
  double  Basic_Nt2_332_Calculator::mT2Fcn::operator()(const std::vector<double>& par) const {
    
    double qT1x = par[0];
    double qT1y = par[1];
    
    double qT2x = theExmiss + thePT1x + thePT2x - qT1x;
    double qT2y = theEymiss + thePT1y + thePT2y - qT1y;
    
    double ETj1 = sqrt(theMass1*theMass1 +  thePT1x*thePT1x + thePT1y*thePT1y);    // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
    double ETj2 = sqrt(theMass2*theMass2 +  thePT2x*thePT2x + thePT2y*thePT2y);     // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless! 
    
    double ETparent1 = sqrt(qT1x*qT1x + qT1y*qT1y + theParent1*theParent1);
    double ETparent2 = sqrt(qT2x*qT2x + qT2y*qT2y + theParent2*theParent2);
    
    double mTsq1 
      = theMass1*theMass1 // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theParent1*theParent1 - 2.0*(ETj1*ETparent1 - (thePT1x*qT1x + thePT1y*qT1y));
    double mTsq2 
      = theMass2*theMass2  // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theParent2*theParent2 - 2.0*(ETj2*ETparent2 - (thePT2x*qT2x + thePT2y*qT2y));
   

    // std::cout << "MOO " << sqrt(fmax(mTsq1,mTsq2)) << "\t" << qT1x << "\t" << qT1y << "\t" << 0.5*(qT1x-qT2x) << "\t" << 0.5*(qT1y-qT2y ) << std::endl;
 
    return -fmin(mTsq1,mTsq2);
    
  }

} // end of Mt2 Namespace

