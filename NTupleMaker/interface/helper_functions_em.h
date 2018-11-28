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
			       TTree* out_tree,
			       float met_x , float met_y,
			       TLorentzVector muonLV,
			       TLorentzVector electronLV,
			       TLorentzVector jet1,
			       TLorentzVector jet2,
			       float uncertainty_container[]
			       )
{


  TLorentzVector metLV; metLV.SetXYZT(met_x, met_y, 0., TMath::Sqrt( met_x*met_x + met_y*met_y ));
  /* TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex], */
  /* 					analysisTree.muon_py[muonIndex], */
  /* 					analysisTree.muon_pz[muonIndex], */
  /* 					classic_svFit::muonMass); */
  /* TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex], */
  /* 						analysisTree.electron_py[electronIndex], */
  /* 						analysisTree.electron_pz[electronIndex], */
  /* 						classic_svFit::electronMass); */
  /* TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet], */
  /* 						   analysisTree.pfjet_py[indexLeadingJet], */
  /* 						   analysisTree.pfjet_pz[indexLeadingJet], */
  /* 						   analysisTree.pfjet_e[indexLeadingJet]); */

  /* TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet], */
  /* 						   analysisTree.pfjet_py[indexSubLeadingJet], */
  /* 						   analysisTree.pfjet_pz[indexSubLeadingJet], */
  /* 						   analysisTree.pfjet_e[indexSubLeadingJet]); */

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
  for (int i = 0; i < 10; ++i) uncertainty_container[i] = -10;

  // met
  out_tree->Branch("met_"+uncertainty_name, &uncertainty_container[0], "met_"+uncertainty_name+"/F");
  uncertainty_container[0] = metLV.E();

  // metphi
  out_tree->Branch("metphi_"+uncertainty_name, &uncertainty_container[1], "metphi_"+uncertainty_name+"/F");
  uncertainty_container[1] = metLV.Phi();
  
  // mTtot
  out_tree->Branch("mTtot_"+uncertainty_name, &uncertainty_container[2], "mTtot_" + uncertainty_name+"/F");
  uncertainty_container[2] = totalTransverseMass( muonLV, electronLV, metLV);

  // mTdileptonMET
  out_tree->Branch("mTdileptonMET_"+uncertainty_name, &uncertainty_container[3], "mTdileptonMET_"+uncertainty_name+"/F");
  uncertainty_container[3] = mT(dileptonLV,metLV);

  // pt_tt
  out_tree->Branch("pt_tt_"+uncertainty_name, &uncertainty_container[4], "pt_tt_"+uncertainty_name+"/F");
  uncertainty_container[4] = (muonLV+electronLV+metLV).Pt();

  // pt_ttjj
  out_tree->Branch("pt_ttjj_"+uncertainty_name, &uncertainty_container[5], "pt_ttjj_"+uncertainty_name+"/F");
  if(jet1.E()!=0 && jet2.E()!=0) uncertainty_container[5] = (muonLV+electronLV+metLV+jet1+jet2).Pt();

  // pzetamiss + dzeta
  out_tree->Branch("pzetamiss_"+uncertainty_name, &uncertainty_container[6], "pzetamiss_"+uncertainty_name+"/F");
  out_tree->Branch("dzeta_"+uncertainty_name, &uncertainty_container[7], "dzeta_"+uncertainty_name+"/F");
  computeDzeta(met_x, met_y, zetaX, zetaY, pzetavis, uncertainty_container[6], uncertainty_container[7]);

  // mt_1
  out_tree->Branch("mt_1_"+uncertainty_name, &uncertainty_container[8], "mt_1_"+uncertainty_name+"/F");
  uncertainty_container[8] = mT(electronLV,metLV);

  // mt_2
  out_tree->Branch("mt_2_"+uncertainty_name, &uncertainty_container[9], "mt_2_"+uncertainty_name+"/F");
  uncertainty_container[9] = mT(muonLV,metLV);

  // mtmax
  out_tree->Branch("mtmax_"+uncertainty_name, &uncertainty_container[10], "mtmax_"+uncertainty_name+"/F");
  uncertainty_container[10] = TMath::Max(float(uncertainty_container[8]),float(uncertainty_container[9]));

  // dphi_emet
  out_tree->Branch("dphi_emet_"+uncertainty_name, &uncertainty_container[11], "dphi_emet_"+uncertainty_name+"/F");
  uncertainty_container[11] = dPhiFrom2P(electronLV.Px(), electronLV.Py(), met_x, met_y);

  // dphi_mumet
  out_tree->Branch("dphi_mumet_"+uncertainty_name, &uncertainty_container[12], "dphi_mumet_"+uncertainty_name+"/F");
  uncertainty_container[12] = dPhiFrom2P(muonLV.Px(), muonLV.Py(), met_x, met_y);

  // pzetavis
  out_tree->Branch("pzetavis_"+uncertainty_name, &uncertainty_container[13], "pzetavis_"+uncertainty_name+"/F");
  uncertainty_container[13] = pzetavis;

  // m_vis
  out_tree->Branch("m_vis_"+uncertainty_name, &uncertainty_container[14], "m_vis_"+uncertainty_name+"/F");
  uncertainty_container[14] = dileptonLV.M();

  // pt_vis
  out_tree->Branch("pt_vis_"+uncertainty_name, &uncertainty_container[15], "pt_vis_"+uncertainty_name+"/F");
  uncertainty_container[15] = dileptonLV.Pt();

  // pt_1
  out_tree->Branch("pt_1_"+uncertainty_name, &uncertainty_container[15], "pt_1_"+uncertainty_name+"/F");
  uncertainty_container[15] = electronLV.Pt();

  // pt_2
  out_tree->Branch("pt_2_"+uncertainty_name, &uncertainty_container[16], "pt_2_"+uncertainty_name+"/F");
  uncertainty_container[16] = muonLV.Pt();
}
