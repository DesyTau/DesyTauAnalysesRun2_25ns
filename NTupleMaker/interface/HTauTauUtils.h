#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <cmath>
#include <string>

#ifndef TauTauUtils
#define TauTauUtils
class HTauTauUtils {

	public:

		HTauTauUtils();
		~HTauTauUtils();
		
		void CalculateU1U2FromMet(float MetPx,
		                          float MetPy,
		                          float genZPx,
		                          float genZPy,
		                          float diLepPx,
		                          float diLepPy,
		                          float & U1,
		                          float & U2,
		                          float & metU1,
		                          float & metU2);

		void CalculateMetFromU1U2(float U1,
		                          float U2,
		                          float genZPx,
		                          float genZPy,
		                          float diLepPx,
		                          float diLepPy,
		                          float & metPx,
		                          float & metPy);

		float cosRestFrame(TLorentzVector boost, TLorentzVector vect); 

		float cosProdPlane(TLorentzVector boson, TLorentzVector lepton, TLorentzVector beam);

		float QToEta(float Q);

		float EtaToQ(float Eta);

		float PtoEta(float Px, float Py, float Pz);

		float dPhiFrom2P(float Px1, float Py1, float Px2, float Py2);

		float dPhiFrom2Phi(float phi1, float phi2);

		float deltaEta(float Px1, float Py1, float Pz1, float Px2, float Py2, float Pz2);

		float deltaPhi(float Phi1, float Phi2);

		float deltaR(float Eta1, float Phi1, float Eta2, float Phi2);

		float PtEtaToP(float Pt, float Eta);

		float InvMassFromP(float px1,float py1, float pz1, float px2, float py2,float pz2);


};
#endif

