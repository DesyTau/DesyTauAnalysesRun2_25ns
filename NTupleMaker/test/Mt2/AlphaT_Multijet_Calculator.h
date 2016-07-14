// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#ifndef MT2_ALPHAT_MULTIJET_CALCULATOR_h
#define MT2_ALPHAT_MULTIJET_CALCULATOR_h

#include <vector>
#include "Mt2/Mt2LorentzTransverseVector.h"
#include <cmath>

namespace Mt2 {

  class AlphaT_Multijet_Calculator {
  public:
    AlphaT_Multijet_Calculator() {
      clear();
    }
    // The "clear()" method resets the calculator, clearing its store of
    // jet momenta.  In order to calculate alphaT, first call "clear()"
    // before pushing_back jet momenta.
    void clear() {
      m_jets.clear();
      m_tot = Mt2::LorentzTransverseVector(0,0,0);
      m_sumPt = 0;
    }
    void push_back(const Mt2::LorentzTransverseVector & jet) {
      m_jets.push_back(jet);
      m_tot = m_tot + jet;
      m_sumPt += jet.pt();
    }
    // The CMS paper http://cdsweb.cern.ch/record/1194509/files/SUS-09-001-pas.pdf defines H_T in one place to be the scalar some of the pts of the input jets, but elsewhere defines it to be the scalar some of the ETs of the input jets.  These two definitions don't agree in the case of massive jets.  I am waiting for clarification from someone in CMS as to the nature of the exact definition.  For the moment, I am assuming that it is the scalar some of the ETs.
    double H_T() const {
      // any sensible person would imagine that this would do:
      // return m_tot.Et();
      // however the CMS definition is
      return m_sumPt;// see  http://cms-physics.web.cern.ch/cms-physics/public/SUS-10-001-pas.pdf	
    }
    double H_T_Sq() const {
      const double HT = H_T();
      return HT*HT;
    }
    double MHT_Sq() /* missing HT magnitude squared */ const {
      return m_tot.ptsq(); // see  http://cms-physics.web.cern.ch/cms-physics/public/SUS-10-001-pas.pdf
    }
    // The CMS paper http://cdsweb.cern.ch/record/1194509/files/SUS-09-001-pas.pdf defines H_T^miss_vec as the two vector sum of the negatives of the transverse two momenta of the jets.
    double alphaT() const {
      if (m_jets.size()<=1) {
	return -1; /// This flags some stupid user behavour!
      };
      //assert(m_jets.size()>=2);
      const unsigned int n = m_jets.size();
      unsigned long deltahtevaluations=0;
      std::vector<bool> combinationsPlusTwoOverTwo(n);
      bszero(combinationsPlusTwoOverTwo);
      combinationsPlusTwoOverTwo[n-1]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle
      bool first=true;
      double smallestDeltaHTSoFar=0;
      //double HTForPseudoJets=0;
      std::vector<bool> comb(n);
      bszero(comb);
      comb[0]=true; // Here we enforce the "Each side has AT LEAST ONE visible particle ... you may want to change that!!!
      for (;
	   comb != combinationsPlusTwoOverTwo;
	   bsincr(comb)) {
	LorentzTransverseVector eta(0,0,0);
	LorentzTransverseVector etb(0,0,0);
	double hta=0;
	double htb=0;
	for (unsigned int bit=0; bit<n; ++bit) {
	  if (comb[bit]) {
	    eta = eta + m_jets[bit];
	    hta += (m_jets[bit]).pt();
	    //if (debug()) {
	    //  std::cout << "A";
	    //}
	  } else {
	    etb = etb + m_jets[bit];
	    htb += (m_jets[bit]).pt();
	    //if (debug()) {
	    //  std::cout << "b";
	    //}
	  }
	}
	//if(debug()) { std::cout << std::endl; }
	// OLD pre  Sarah Williams alphaT 25 definition:
	//const double deltaht = fabs(eta.pt()-etb.pt());
	// NEW post Sarah Williams alphaT 25 definition:
	const double deltaht = fabs(hta-htb);
	++deltahtevaluations;
	if (first || deltaht<smallestDeltaHTSoFar) {
	  first=false;
	  smallestDeltaHTSoFar = deltaht;
	  //HTForPseudoJets = eta.pt() + etb.pt(); // There is a better description of why this is pt rather than et in http://cms-physics.web.cern.ch/cms-physics/public/SUS-10-001-pas.pdf     What follows is an older comment. LOOK AT THIS CRAZY DEFINITION ... ONE WOULD IMAGINE CMS WANTED ET+ET, BUT IF WE DO THAT ALPHAT CAN BE BIGGER THAN 0.5.  This means in the Dijet alphaT definition H_T really means p_T whereas when it is ued in equation (2) of http://cdsweb.cern.ch/record/1194509/files/SUS-09-001-pas.pdf it instead means E_T ... YUK YUK YUK.

	}
      }
      if (debug()) {
	std::cerr << "Used " << deltahtevaluations << " mt2 evaluations" << std::endl;
      }
      return 0.5*(H_T()-smallestDeltaHTSoFar)/sqrt(H_T_Sq() - MHT_Sq());
      
      
    }
  private:
    bool debug() const { return false; }
    // There is overlap between code in the two methods below and code in Basic_MtGen_332_Calculator.h and Basic_MtGen_330_Calculator.h that should be remove.  FIXME
    static void bszero(std::vector<bool> & bs) {
      for (unsigned int i=0; i<bs.size(); ++i) {
	bs[i]=false;
      }
    }
    static void bsincr(std::vector<bool> & bs) {
      for(unsigned long bit=0; bit<bs.size(); ++bit) {
	bs[bit] = (!(bs[bit]));
	if (bs[bit]) { 
	  // we incremented 0 to 1, so we are done
	  return;
	} 
	// we incremented 1 to 0, so we had better add the carry bit...
      }
      // we overflowed, but do not complain
    }
    
  private:
    std::vector<Mt2::LorentzTransverseVector> m_jets;
    Mt2::LorentzTransverseVector m_tot;
    double m_sumPt; 
  };
  
}

#endif
