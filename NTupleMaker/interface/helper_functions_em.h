#include <iostream>
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

void computeDzeta(float metX,  float metY, float zetaX, float zetaY, float pzetavis,
                  float & pzetamiss,
                  float & dzeta) {

    pzetamiss = metX*zetaX + metY*zetaY;
    dzeta = pzetamiss - 0.85*pzetavis;
}

float totalTransverseMass(TLorentzVector l1, TLorentzVector l2, TLorentzVector l3) {

    TLorentzVector totalLV = l1 + l2 + l3;
    float totalET = l1.Pt() +  l2.Pt() + l3.Pt();
    float mTtot = TMath::Sqrt(totalET*totalET-totalLV.Pt()*totalLV.Pt());
    return  mTtot;
}


float topPtWeight(float pt1, float pt2, bool run1) {

  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter

  if (run1) {
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
    a = 0.156;    // Run1 a parameter
    b = -0.00137;  // Run1 b parameter
  }
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);

  return TMath::Sqrt(w1*w2);

}

ClassicSVfit SVFitMassComputation(classic_svFit::MeasuredTauLepton svFitEle,
                                              classic_svFit::MeasuredTauLepton svFitMu,
                                              double measuredMVAMETx,
                                              double measuredMVAMETy,
                                              TMatrixD covMVAMET,
                                              TFile * inputFile_visPtResolution
                                              ) {
    
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitEle);
    measuredTauLeptons.push_back(svFitMu);
    
    int verbosity = 0;
    ClassicSVfit svFitAlgo(verbosity);
    double kappa = 3.; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
    svFitAlgo.addLogM_fixed(true, kappa);
    svFitAlgo.integrate(measuredTauLeptons, measuredMVAMETx, measuredMVAMETy, covMVAMET);
    
    return svFitAlgo;
    
}



bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags) {
    
    bool passed = true;
    unsigned int nFlags = metFlags.size();
    //  std::cout << "MEt filters : " << std::endl;
    for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
        TString flagName(it->first);
        //    std::cout << it->first << " : " << it->second << std::endl;
        for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
            if (flagName.Contains(metFlags[iFilter])) {
                if (it->second==0) {
                    passed = false;
                    break;
                }
            }
        }
    }
    //  std::cout << "Passed : " << passed << std::endl;
    return passed;
    
}

//class MetPropagatedToVariable :

void propagate_uncertainty(TString uncertainty_name,
			   TLorentzVector metLV, TMatrixD covMET, TFile* inputFile_visPtResolution,
			   TLorentzVector muonLV,
			   TLorentzVector electronLV,
			   TLorentzVector jet1LV,
			   TLorentzVector jet2LV,
			   float uncertainty_container[],
			   bool isData,
			   bool svfit_on
			   )
{

  TLorentzVector dileptonLV = muonLV + electronLV;
  float electronUnitX = electronLV.Px()/electronLV.Pt();
  float electronUnitY = electronLV.Py()/electronLV.Pt();
  float muonUnitX = muonLV.Px()/muonLV.Pt();
  float muonUnitY = muonLV.Py()/muonLV.Pt();
  float zetaX = electronUnitX + muonUnitX;
  float zetaY = electronUnitY + muonUnitY;
  float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
  zetaX = zetaX/normZeta;
  zetaY = zetaY/normZeta;
  float vectorVisX = muonLV.Px()+electronLV.Px();
  float vectorVisY = muonLV.Py()+electronLV.Py();
  float pzetavis = vectorVisX*zetaX+vectorVisY*zetaY;

  // initliaze all values to -10
  for (int i = 0; i < 50; ++i) uncertainty_container[i] = -10;

  // met // 0
  uncertainty_container[0] = metLV.E();
  // metphi // 1
  uncertainty_container[1] = metLV.Phi();
  // mTtot // 2
  uncertainty_container[2] = totalTransverseMass( muonLV, electronLV, metLV);
  // mTdileptonMET // 3
  uncertainty_container[3] = mT(dileptonLV,metLV);
  // pt_tt // 4
  uncertainty_container[4] = (muonLV+electronLV+metLV).Pt();
  // pt_ttjj // 5
  if(jet1LV.E()!=0 && jet2LV.E()!=0) uncertainty_container[5] = (muonLV+electronLV+metLV+jet1LV+jet2LV).Pt();
  // pzetamiss // 6 + dzeta // 7
  computeDzeta(metLV.Px(), metLV.Py(), zetaX, zetaY, pzetavis, uncertainty_container[6], uncertainty_container[7]);
  // mt_1 // 8
  uncertainty_container[8] = mT(electronLV,metLV);
  // mt_2 // 9
  uncertainty_container[9] = mT(muonLV,metLV);
  // mtmax // 10
  uncertainty_container[10] = TMath::Max(float(mT(electronLV,metLV)),float(mT(muonLV,metLV)));
  // dphi_emet // 11
  uncertainty_container[11] = dPhiFrom2P(electronLV.Px(), electronLV.Py(), metLV.Px(), metLV.Py());
  // dphi_mumet // 12
  uncertainty_container[12] = dPhiFrom2P(muonLV.Px(), muonLV.Py(), metLV.Px(), metLV.Py());
  // pzetavis // 13
  uncertainty_container[13] = pzetavis;
  // m_vis // 14
  uncertainty_container[14] = dileptonLV.M();
  // pt_vis // 15
  uncertainty_container[15] = dileptonLV.Pt();
  // pt_1 // 16
  uncertainty_container[16] = electronLV.Pt();
  // pt_2 // 17
  uncertainty_container[17] = muonLV.Pt();
  // jpt_1 // 18
  if(jet1LV.E()!=0) uncertainty_container[18] = jet1LV.Pt();
  // jpt_2 // 19
  if(jet2LV.E()!=0) uncertainty_container[19] = jet2LV.Pt();
  // mjj // 20
  if(jet1LV.E()!=0 && jet2LV.E()!=0) uncertainty_container[20] = (jet1LV+jet2LV).M();;
  // dijetphi // 21
  if(jet1LV.E()!=0 && jet2LV.E()!=0) uncertainty_container[21] = (jet1LV+jet2LV).Phi();;
  // dijetpt // 22
  if(jet1LV.E()!=0 && jet2LV.E()!=0) uncertainty_container[22] = (jet1LV+jet2LV).Pt();;
  // m_sv // 23, pt_sv // 24, eta_sv // 25, phi_sv // 26, mt_sv // 27
  bool checkSV = svfit_on && uncertainty_container[7]>-35 && uncertainty_container[3]<60;
  if(!isData && checkSV){
    classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, electronLV.Pt(), electronLV.Eta(), electronLV.Phi(), 0.51100e-3);
    classic_svFit::MeasuredTauLepton svFitMu(classic_svFit::MeasuredTauLepton::kTauToMuDecay, muonLV.Pt(), muonLV.Eta(), muonLV.Phi(), 105.658e-3);
    ClassicSVfit algo = SVFitMassComputation(svFitEle, svFitMu, metLV.Px(), metLV.Py(), covMET, inputFile_visPtResolution);
    uncertainty_container[23] = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass();
    uncertainty_container[24] = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt(); 
    uncertainty_container[25]= static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
    uncertainty_container[26]= static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
    uncertainty_container[27]= static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass();
  }
  // mTemu // 28
  uncertainty_container[28] = mT(electronLV,muonLV);
}
