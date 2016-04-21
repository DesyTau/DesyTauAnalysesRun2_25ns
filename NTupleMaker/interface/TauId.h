#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"


float leadchargedhadrcand_dz = 0.2;
float leadchargedhadrcand_dxy = 0;
float decayModeFinding  = 0.5;
float decayModeFindingNewDMs  = 0.5;
float againstElectronVLooseMVA5  = 0.5;
float againstMuonTight3  = 0.5;
float againstMuonLoose3  = 0.5;
float vertexz =  0.2;
float byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;



bool tightTauID(AC1B &tree_, int &tau){
      bool isTight = false;

	if (tree_.tau_decayModeFindingNewDMs[tau]<decayModeFindingNewDMs) continue;
	if ( fabs(tree_.tau_leadchargedhadrcand_dz[tau])> leadchargedhadrcand_dz) continue;
	if (tree_.tau_againstElectronVLooseMVA5[tau]<againstElectronVLooseMVA5) continue;
	if (tree_.tau_againstMuonTight3[tau]<againstMuonTight3) continue;
	
	if (tree_.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau] > byCombinedIsolationDeltaBetaCorrRaw3Hits ) continue;

	double  tauIso = tree_.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau];


      	isTight = true;
		return isTight;
     }	

bool looseTauID(AC1B &tree_, int &tau){
	bool isLoose = false;

	if (tree_.tau_decayModeFindingNewDMs[tau]<decayModeFindingNewDMs) continue;
	if ( fabs(tree_.tau_leadchargedhadrcand_dz[tau])> leadchargedhadrcand_dz) continue;
	if (tree_.tau_againstElectronVLooseMVA5[tau]<againstElectronVLooseMVA5) continue;
	if (tree_.tau_againstMuonTight3[tau]<againstMuonTight3) continue;
	
	if (tree_.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau] < byCombinedIsolationDeltaBetaCorrRaw3Hits ) continue;
        if (tauIso > 2*byCombinedIsolationDeltaBetaCorrRaw3Hits ) continue; 


      	isLoose = true;
 		return isLoose;     
     }	
