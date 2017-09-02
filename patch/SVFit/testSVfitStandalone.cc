
/**
   \class testSVfitStandalone testSVfitStandalone.cc "TauAnalysis/SVfitStandalone/bin/testSVfitStandalone.cc"
   \brief Basic example of the use of the standalone version of SVfit

   This is an example executable to show the use of the standalone version of SVfit 
   from a flat n-tuple or single event.
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

void singleEvent()
{

  /* 
     This is a single event for testing in the integration mode.
  */
  // define MET

  /*
  float met = 84.32;
  float metphi = -1.3005;

  float cov00 = 622.057;
  float cov01 = 26.5763;
  float cov11 = 710.285;

  float pt_1 = 72.2457;
  float eta_1 = 0.321129;
  float phi_1 = -1.88906;

  float pt_2 = 95.3255;
  float eta_2 = 1.36582;
  float phi_2 = 0.28302;



  float met = 14.3559;
  float metphi = -2.32674;

  float cov00 = 558.561;
  float cov01 = -29.2334;
  float cov11 = 511.677;

  float pt_1 = 38.0104;
  float eta_1 = 1.1496;
  float phi_1 = -2.93753;

  float pt_2 = 17.5338;
  float eta_2 = 2.10357;
  float phi_2 = 2.43002;


  float met = 38.1073;
  float metphi = -2.7665;

  float cov00 = 1133.65;
  float cov01 = 173.975;
  float cov11 = 922.363;

  float pt_1 = 17.4679;
  float eta_1 = 0.86439;
  float phi_1 = -2.54987;

  float pt_2 = 34.9185;
  float eta_2 = -0.9838;
  float phi_2 = 0.541699;
  */

  float met = 24.2576;
  float metphi = 0.718649;

  float cov00 = 757.51;
  float cov01 = 50.324;
  float cov11 = 738.451;

  float pt_1 = 27.9314;
  float eta_1 = -1.56379;
  float phi_1 = 0.168914;

  float pt_2 = 10.8825;
  float eta_2 = -2.05159;
  float phi_2 = 0.994835;

  float measuredMETx =  met*TMath::Cos(metphi);
  float measuredMETy =  met*TMath::Sin(metphi); 
  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] =  cov00;
  covMET[1][0] =  cov01;
  covMET[0][1] =  cov01;
  covMET[1][1] =  cov11;
  // define lepton four vectors
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, pt_1, eta_1, phi_1, 0.51100e-3)); // tau -> electron decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay,  pt_2, eta_2, phi_2, 105.658e-3)); // tau -> muon (Pt, eta, phi, mass)
  // define algorithm (set the debug level to 3 for testing)
  unsigned verbosity = 2;
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, verbosity);
  algo.addLogM(false);  
  //  algo.addLogM(true, 1.);
  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  algo.shiftVisPt(true, inputFile_visPtResolution);
  //algo.shiftVisPt2(true);
  /* 
     the following lines show how to use the different methods on a single event
  */
  // minuit fit method
  //algo.fit();
  // integration by VEGAS (same as function algo.integrate() that has been in use when markov chain integration had not yet been implemented)
  //algo.integrateVEGAS();
  // integration by markov chain MC
  algo.integrateMarkovChain();

  float mass = algo.getMass(); // full mass of tau lepton pair in units of GeV
  float transverseMass = algo.transverseMass(); // transverse mass of tau lepton pair in units of GeV
  if ( algo.isValidSolution() ) {
    std::cout << "found valid solution: mass = " << mass << " (expected value = 124.646), transverse mass = " << transverseMass << " (expected value = 123.026)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }

  delete inputFile_visPtResolution;
}

void eventsFromTree(int argc, char* argv[]) 
{
  // parse arguments
  if ( argc < 3 ) {
    std::cout << "Usage : " << argv[0] << " [inputfile.root] [tree_name]" << std::endl;
    return;
  }
  // get intput directory up to one before mass points
  TFile* file = new TFile(argv[1]); 
  // access tree in file
  TTree* tree = (TTree*) file->Get(argv[2]);
  // input variables
  float met, metPhi;
  float covMet11, covMet12; 
  float covMet21, covMet22;
  float l1Pt, l1Eta, l1Phi, l1Mass;
  float l2Pt, l2Eta, l2Phi, l2Mass;
  float mTrue;
  // branch adresses
  tree->SetBranchAddress("met", &met);
  tree->SetBranchAddress("metphi", &metPhi);
  tree->SetBranchAddress("metcov00", &covMet11);
  tree->SetBranchAddress("metcov01", &covMet12);
  tree->SetBranchAddress("metcov10", &covMet21);
  tree->SetBranchAddress("metcov11", &covMet22);
  tree->SetBranchAddress("pt_1", &l1Pt);
  tree->SetBranchAddress("eta_1", &l1Eta);
  tree->SetBranchAddress("phi_1", &l1Phi);
  tree->SetBranchAddress("m_1", &l1Mass);
  tree->SetBranchAddress("pt_2", &l2Pt);
  tree->SetBranchAddress("eta_2", &l2Eta);
  tree->SetBranchAddress("phi_2", &l2Phi);
  tree->SetBranchAddress("m_2", &l2Mass);
  tree->SetBranchAddress("m_sv", &mTrue);
  int nevent = tree->GetEntries();
  for ( int i = 0; i < nevent; ++i ) {
    tree->GetEvent(i);
    std::cout << "event " << (i + 1) << std::endl;

    if (mTrue<0) continue;

    // setup MET input vector
    double measuredMETx = met*TMath::Cos(metPhi);
    double measuredMETy = met*TMath::Sin(metPhi);
    // setup the MET significance
    TMatrixD covMET(2,2);
    covMET[0][0] = covMet11;
    covMET[0][1] = covMet12;
    covMET[1][0] = covMet21;
    covMET[1][1] = covMet22;
    // setup measure tau lepton vectors 
    
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, l1Pt, l1Eta, l1Phi, 0.51100e-3));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, l2Pt, l2Eta, l2Phi, 105.658e-3));
    // construct the class object from the minimal necesarry information
    SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
    // apply customized configurations if wanted (examples are given below)
    //    algo.maxObjFunctionCalls(5000);
    edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
    TH1::AddDirectory(false);  
    TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
    algo.shiftVisPt(true, inputFile_visPtResolution);
    algo.addLogM(false);
    //algo.metPower(0.5)
    // minuit fit method
    //algo.fit();
    // integration by VEGAS (default)
    //    algo.integrateVEGAS();
    // integration by markov chain MC
    algo.integrateMarkovChain();
    // retrieve the results upon success
    std::cout << "... m svfit from tree : " << mTrue << std::endl;
    if ( algo.isValidSolution() ) {
      std::cout << "... m svfit computed  : " << algo.mass() << std::endl; // return value is in units of GeV
    } else {
      std::cout << "... m svfit : ---" << std::endl;
    }
  }
  return;
}

int main(int argc, char* argv[]) 
{
  eventsFromTree(argc, argv);
  //  singleEvent();
  return 0;
}
