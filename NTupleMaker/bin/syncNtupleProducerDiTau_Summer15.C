#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TCut.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TH1F.h"

#define DEBUG true

//Preselection 
//
//Requirement Objects must pass lepton ids and default kinematic selection

void synchNtuple(string sample = "GGFH125", string stream = "MuTau", bool incl=false) {

  TFile *lOFile = new TFile(("SYNC_"+sample+".root").c_str(),"RECREATE");
  TTree *lOTree = new TTree("tree","tree");

   //Bookeeping
   int   lRun         = 0; lOTree->Branch("run"        ,&lRun           ,"lRun/I"     );//Run
   int   lLumi        = 0; lOTree->Branch("lumi"       ,&lLumi          ,"lLumi/I"    );//Lumi
   int   lEvt         = 0; lOTree->Branch("evt"        ,&lEvt           ,"lEvt/I"     );//Evt

   //Event Variables
   int   lNPV         = 0; lOTree->Branch("npv"        ,&lNPV           ,"lNPV/I"     );//NPV
   float lNPU         = 0; lOTree->Branch("npu"        ,&lNPU           ,"lNPU/F"     );//NPU
   float lRho         = 0; lOTree->Branch("rho"        ,&lRho           ,"lRho/F"     );//Rho
   
   //Event Weights
   float lMCWeight    = 0; lOTree->Branch("mcweight"   ,&lMCWeight      ,"lMCWeight/F");//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
   float lPUWeight    = 0; lOTree->Branch("puweight"   ,&lPUWeight      ,"lPUWeight/F");//Pielup Weight
   float lEffWeight   = 0; lOTree->Branch("effweight"  ,&lEffWeight     ,"lEffWeight/F");//Effieiency Scale factor (all components multiplied in)
   float lWeight      = 0; lOTree->Branch("weight"     ,&lWeight        ,"lWeight/F"  );//mcweight*puweight*effweight

   float lTrigweight_1 = 0; lOTree->Branch("trigweight_1"     ,&lTrigweight_1        ,"lTrigweight_1/F"  );
   float lIdweight_1 = 0; lOTree->Branch("idweight_1"     ,&lIdweight_1        ,"lIdweight_1/F"  );
   float lIsoweight_1 = 0; lOTree->Branch("isoweight_1"     ,&lIsoweight_1        ,"lIsoweight_1/F"  );
   float lTrigweight_2 = 0; lOTree->Branch("trigweight_2"     ,&lTrigweight_2        ,"lTrigweight_2/F"  ); 
   float lIdweight_2 = 0; lOTree->Branch("idweight_2"     ,&lIdweight_2        ,"lIdweight_2/F"  ); 
   float lIsoweight_2 = 0; lOTree->Branch("isoweight_2"     ,&lIsoweight_2        ,"lIsoweight_2/F"  );

   //SV Fit variables
   float lMSV         = 0; lOTree->Branch("m_sv"       ,&lMSV           ,"lMSV/F"     );//SV Fit using integration method
   float lMSVUp       = 0; lOTree->Branch("m_sv_Up"    ,&lMSVUp         ,"lMSVUp/F"   );//High Energy scale shape
   float lMSVDown     = 0; lOTree->Branch("m_sv_Down"  ,&lMSVDown       ,"lMSVDown/F" );//Low Energy Scale Shape
 
   float lPtSV         = 0; lOTree->Branch("pt_sv"       ,&lPtSV           ,"lPtSV/F"     );//SV Fit using integration method
   float lEtaSV        = 0; lOTree->Branch("eta_sv"      ,&lEtaSV          ,"lEtaSV/F"     );
   float lPhiSV        = 0; lOTree->Branch("phi_sv"      ,&lPhiSV          ,"lPhiSV/F"     );
   float lMetSV        = 0; lOTree->Branch("met_sv"      ,&lMetSV          ,"lMetSV/F"     );//using MarkovChain MC integration svfit.fittedMET().Rho()
   float lPtH          = 0; lOTree->Branch("pth"         ,&lPtH            ,"lPtH/F"     );//

   ///First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
   float lPt1         = 0; lOTree->Branch("pt_1"       ,&lPt1           ,"lPt1/F"     ); //pT 
   float lPhi1        = 0; lOTree->Branch("phi_1"      ,&lPhi1          ,"lPhi1/F"    ); //Phi 
   float lEta1        = 0; lOTree->Branch("eta_1"      ,&lEta1          ,"lEta1/F"    ); //Eta 
   float lM1          = 0; lOTree->Branch("m_1"        ,&lM1            ,"lM1/F"      ); //Mass 
   float lIso1        = 0; lOTree->Branch("iso_1"      ,&lIso1          ,"lIso1/F"    ); //Delta Beta iso value 
   //float lMVA1        = 0; lOTree->Branch("mva_1"      ,&lMVA1          ,"lMVA1/F"   );//MVA id (when using electron) 0 otherwise
   float lD01         = 0; lOTree->Branch("d0_1"       ,&lD01           ,"lD01/F"      );//d0 with respect to primary vertex
   float lDZ1         = 0; lOTree->Branch("dZ_1"       ,&lDZ1           ,"lDZ1/F"      );//dZ with respect to primary vertex
   //bool  lPassId1     = 0; lOTree->Branch("passid_1"   ,&lPassId1       ,"lPassId1/B" );//Whether it passes id  (not necessarily iso)
   //bool  lPassIso1    = 0; lOTree->Branch("passiso_1"  ,&lPassIso1      ,"lPassIso1/B");//Whether it passes iso (not necessarily id)
   float lMt1         = 0; lOTree->Branch("mt_1"       ,&lMt1           ,"lMt1/F"     );//mT of  first lepton wrt to pf met
   float lMtMVA1      = 0; lOTree->Branch("mtmva_1"    ,&lMtMVA1        ,"lMtMVA1/F"     );//mT of  first lepton wrt to mva met    

   float lQ1          = 0; lOTree->Branch("q_1"        ,&lQ1            ,"lQ1/F"      ); //Charge
   float lagainstElectronLooseMVA5_1  = 0; lOTree->Branch("againstElectronLooseMVA5_1",   &lagainstElectronLooseMVA5_1, "lagainstElectronLooseMVA5_1/F") ;
   float lagainstElectronMediumMVA5_1 = 0; lOTree->Branch("againstElectronMediumMVA5_1",   &lagainstElectronMediumMVA5_1, "lagainstElectronMediumMVA5_1/F");
   float lagainstElectronTightMVA5_1  = 0; lOTree->Branch("againstElectronTightMVA5_1",   &lagainstElectronTightMVA5_1, "lagainstElectronTightMVA5_1/F");
   float lagainstElectronVLooseMVA5_1 = 0; lOTree->Branch("againstElectronVLooseMVA5_1",   &lagainstElectronVLooseMVA5_1, "lagainstElectronVLooseMVA5_1/F");
   float lagainstMuonLoose3_1         = 0; lOTree->Branch("againstMuonLoose3_1",     &lagainstMuonLoose3_1,   "lagainstMuonLoose3_1/F");
   float lagainstMuonTight3_1         = 0; lOTree->Branch("againstMuonTight3_1",     &lagainstMuonTight3_1,   "lagainstMuonTight3_1/F");
   float lbyCombinedIsolationDeltaBetaCorrRaw3Hits_1 = 0; lOTree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",   &lbyCombinedIsolationDeltaBetaCorrRaw3Hits_1,   "lbyCombinedIsolationDeltaBetaCorrRaw3Hits_1/F");
   float ldecayModeFinding_1 = 0;  lOTree->Branch("decayModeFinding_1",  &ldecayModeFinding_1,   "ldecayModeFinding_1/F");
   float ldecayModeFindingNewDMs_1 = 0; lOTree->Branch("decayModeFindingNewDMs_1",    &ldecayModeFindingNewDMs_1,     "ldecayModeFindingNewDMs_1/F");

   ///Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
   float lPt2         = 0; lOTree->Branch("pt_2"       ,&lPt2           ,"lPt2/F"     );//pT
   float lPhi2        = 0; lOTree->Branch("phi_2"      ,&lPhi2          ,"lPhi2/F"    );//Phi
   float lEta2        = 0; lOTree->Branch("eta_2"      ,&lEta2          ,"lEta2/F"    );//Eta
   float lM2          = 0; lOTree->Branch("m_2"        ,&lM2            ,"lM2/F"      );//Mass (visible mass for hadronic Tau)
   float lIso2        = 0; lOTree->Branch("iso_2"      ,&lIso2          ,"lIso2/F"    );//MVA iso for hadronic Tau, Delta Beta for muon
   float lD02         = 0; lOTree->Branch("d0_2"       ,&lD02           ,"lD02/F"      );//d0 with respect to primary vertex
   float lDZ2         = 0; lOTree->Branch("dZ_2"       ,&lDZ2           ,"lDZ2/F"      );//dZ with respect to primary vertex 
   //float lMVA2        = 0; lOTree->Branch("mva_2"      ,&lMVA2          ,"lMVA2/F"   );//MVA id (for anti electron id)
   //bool  lPassId2     = 0; lOTree->Branch("passid_2"   ,&lPassId2       ,"lPassId2/B" );//Whether it passes id  (not necessarily iso)
   //bool  lPassIso2    = 0; lOTree->Branch("passiso_2"  ,&lPassIso2      ,"lPassIso2/B");//Whether it passes iso (not necessarily id)
   float lMt2         = 0; lOTree->Branch("mt_2"       ,&lMt2           ,"lMt2/F"     );//mT of 2nd lepton wrt to pf met
   float lMtMVA2      = 0; lOTree->Branch("mtmva_2"    ,&lMtMVA2        ,"lMtMVA2/F"     );//mT of  first lepton wrt to mva met
   float lQ2          = 0; lOTree->Branch("q_2"        ,&lQ2            ,"lQ2/F"      ); //Charge
   float lagainstElectronLooseMVA5_2  = 0; lOTree->Branch("againstElectronLooseMVA5_2",   &lagainstElectronLooseMVA5_2, "lagainstElectronLooseMVA5_2/F");
   float lagainstElectronMediumMVA5_2 = 0; lOTree->Branch("againstElectronMediumMVA5_2",   &lagainstElectronMediumMVA5_2, "lagainstElectronMediumMVA5_2/F");
   float lagainstElectronTightMVA5_2  = 0; lOTree->Branch("againstElectronTightMVA5_2",   &lagainstElectronTightMVA5_2, "lagainstElectronTightMVA5_2/F");
   float lagainstElectronVLooseMVA5_2 = 0; lOTree->Branch("againstElectronVLooseMVA5_2",   &lagainstElectronVLooseMVA5_2, "lagainstElectronVLooseMVA5_2/F");
   float lagainstMuonLoose3_2         = 0; lOTree->Branch("againstMuonLoose3_2",     &lagainstMuonLoose3_2,   "lagainstMuonLoose3_2/F");
   float lagainstMuonTight3_2         = 0; lOTree->Branch("againstMuonTight3_2",     &lagainstMuonTight3_2,   "lagainstMuonTight3_2/F");
   float lbyCombinedIsolationDeltaBetaCorrRaw3Hits_2 = 0; lOTree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",   &lbyCombinedIsolationDeltaBetaCorrRaw3Hits_2,   "lbyCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
   float ldecayModeFinding_2 = 0;  lOTree->Branch("decayModeFinding_2",  &ldecayModeFinding_2,   "ldecayModeFinding_2/F");
   float ldecayModeFindingNewDMs_2 = 0; lOTree->Branch("decayModeFindingNewDMs_2",    &ldecayModeFindingNewDMs_2,     "ldecayModeFindingNewDMs_2/F");

   //veto lepton
   int lNVetoMuon     ; lOTree->Branch("nVetoMuon",     &lNVetoMuon );
   int lNVetoElectron ; lOTree->Branch("nVetoElectron",      &lNVetoElectron );
   bool lExtraelec_veto; lOTree->Branch("extraelec_veto",    &lExtraelec_veto);
   bool lExtramuon_veto; lOTree->Branch("extramuon_veto",    &lExtramuon_veto);

   //Met related variables
   float lMet         = 0; lOTree->Branch("met"        ,&lMet           ,"lMet/F"      ); //pfmet
   float lMetPhi      = 0; lOTree->Branch("metphi"     ,&lMetPhi        ,"lMetPhi/F"   ); //pfmet Phi
   float lMVAMet      = 0; lOTree->Branch("mvamet"     ,&lMVAMet        ,"lMet/F"      ); //mvamet
   float lMVAMetPhi   = 0; lOTree->Branch("mvametphi"  ,&lMVAMetPhi     ,"lMetPhi/F"   ); //mvamet Phi
   //float lPZetaVis    = 0; lOTree->Branch("pzetavis"   ,&lPZetaVis      ,"lPZetaVis/F" ); //pZeta Visible
   //float lPZetaMiss   = 0; lOTree->Branch("pzetamiss"  ,&lPZetaMiss     ,"lPZetaMiss/F"); //pZeta Missing
   //MET covariance matrices
   float lMetCov00   = 0; lOTree->Branch("metcov00"   ,&lMetCov00      ,"lMetCov00/F"); //pf met covariance matrix 00 
   float lMetCov01    = 0; lOTree->Branch("metcov01"   ,&lMetCov01      ,"lMetCov01/F"); //pf met covariance matrix 01 
   float lMetCov10    = 0; lOTree->Branch("metcov10"   ,&lMetCov10      ,"lMetCov10/F"); //pf met covariance matrix 10 
   float lMetCov11    = 0; lOTree->Branch("metcov11"   ,&lMetCov11      ,"lMetCov11/F"); //pf met covariance matrix 11 
   //MVAMet covariance matrices
   float lMVACov00    = 0; lOTree->Branch("mvacov00"   ,&lMVACov00      ,"lMVACov00/F"); //mva met covariance matrix 00 
   float lMVACov01    = 0; lOTree->Branch("mvacov01"   ,&lMVACov01      ,"lMVACov01/F"); //mva met covariance matrix 01 
   float lMVACov10    = 0; lOTree->Branch("mvacov10"   ,&lMVACov10      ,"lMVACov10/F"); //mva met covariance matrix 10 
   float lMVACov11    = 0; lOTree->Branch("mvacov11"   ,&lMVACov11      ,"lMVACov11/F"); //mva met covariance matrix 11 
   

   //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
   float lJPt1       = 0; lOTree->Branch("jpt_1"      ,&lJPt1          ,"lJPt1/F"     );//Jet Pt after corrections
   float lJEta1      = 0; lOTree->Branch("jeta_1"     ,&lJEta1         ,"lJEta1/F"    );//Jet Eta
   float lJPhi1      = 0; lOTree->Branch("jphi_1"     ,&lJPhi1         ,"lJPhi1/F"    );//Jet Phi     
   float lJPtRaw1    = 0; lOTree->Branch("jptraw_1"   ,&lJPtRaw1       ,"lJPtRaw1/F"  );//Jet Raw Pt (before corrections)
   float lJPtUnc1    = 0; lOTree->Branch("jptunc_1"   ,&lJPtUnc1       ,"lJPtUnc1/F"  );//Jet Unc (relative to Jet corrected pT)
   float lJMVA1      = 0; lOTree->Branch("jmva_1"     ,&lJMVA1         ,"lJMVA1/F"    );//Jet MVA id value
   bool  lJPass1     = 0; lOTree->Branch("jpass_1"    ,&lJPass1        ,"lJPass1/B"   );//Whether Jet pass PU Id Loose WP
   float lJCsv1      = 0; lOTree->Branch("jcsv_1"     ,&lJCsv1         ,"lJCsv1/F"     );//CSV discriminator

   //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
   float lJPt2       = 0; lOTree->Branch("jpt_2"      ,&lJPt2          ,"lJPt2/F"     );//Jet Pt after corrections
   float lJEta2      = 0; lOTree->Branch("jeta_2"     ,&lJEta2         ,"lJEta2/F"    );//Jet Eta
   float lJPhi2      = 0; lOTree->Branch("jphi_2"     ,&lJPhi2         ,"lJPhi2/F"    );//Jet Phi
   float lJPtRaw2    = 0; lOTree->Branch("jptraw_2"   ,&lJPtRaw2       ,"lJPtRaw2/F"  );//Jet Raw Pt (before corrections)
   float lJPtUnc2    = 0; lOTree->Branch("jptunc_2"   ,&lJPtUnc2       ,"lJPtUnc2/F"  );//Jet Unc (relative to Jet corrected pT)
   float lJMVA2      = 0; lOTree->Branch("jmva_2"     ,&lJMVA2         ,"lJMVA2/F"    );//Jet MVA id value
   bool  lJPass2     = 0; lOTree->Branch("jpass_2"    ,&lJPass2        ,"lJPass2/B"   );//Whether jet passes PU Id Loose WP 
   float lJCsv2      = 0; lOTree->Branch("jcsv_2"     ,&lJCsv2         ,"lJCsv2/F"     );//CSV discriminator

   //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
   float lBTagPt     = 0; lOTree->Branch("bpt_1"        ,&lBTagPt        ,"lBTagPt/F"   );//Corrected BTag Pt
   float lBTagEta    = 0; lOTree->Branch("beta_1"       ,&lBTagEta       ,"lBTagEta/F"  );//Btag Eta
   float lBTagPhi    = 0; lOTree->Branch("bphi_1"       ,&lBTagPhi       ,"lBTagPhi/F"  );//Btag Phi
 
   //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets 
   float lMJJ        = 0; lOTree->Branch("mjj"        ,&lMJJ           ,"lMJJ/F"      );//Mass Di Jet system  
   float lJDEta      = 0; lOTree->Branch("jdeta"      ,&lJDEta         ,"lJDEta/F"    );//|jeta_1-jeta_2| 
   int   lNJetInGap  = 0; lOTree->Branch("njetingap"  ,&lNJetInGap     ,"lNJetInGap/I");//# of Jets between two jets
   //float lMVA        = 0; lOTree->Branch("mva"        ,&lMVA           ,"lMVA/F"      );//VBF MVA value
   
   //Variables that go into the VBF MVA
   float lJDPhi      = 0; lOTree->Branch("jdphi"      ,&lJDPhi         ,"lJDPhi/F"    );//Delta Phi between two leading jets
   float lDiJetPt    = 0; lOTree->Branch("dijetpt"    ,&lDiJetPt       ,"lDiJetPt/F"  );//Pt of the di jet system
   float lDiJetPhi   = 0; lOTree->Branch("dijetphi"   ,&lDiJetPhi      ,"lDiJetPhi/F" );//Phi of the di jet system
   float lHDJetPhi   = 0; lOTree->Branch("hdijetphi"  ,&lHDJetPhi      ,"lHDJetPhi/F" );//Phi of the di jet system - Higgs system phi
   float lVisJetEta  = 0; lOTree->Branch("visjeteta"  ,&lVisJetEta     ,"lVisJetEta/F");//TMath::Min(eta_vis - jeta,eta_vis,jeta2);
   float lPtVis      = 0; lOTree->Branch("ptvis"      ,&lPtVis         ,"lPtVis/F"    );//Pt Vis
  
   //number of btags passing btag id ( pt > 20 )
   int   lNBTag      = 0; lOTree->Branch("nbtag"      ,&lNBTag         ,"lNBTag/I");

   //number of jets passing jet id ( pt > 30 )
   int   lNJets      = 0; lOTree->Branch("njets"      ,&lNJets         ,"lNJets/I");
   int   lNJetsPt20  = 0; lOTree->Branch("njetspt20"      ,&lNJetsPt20         ,"lNJetsPt20/I");

   int   lMuFlag     = 0; 
   if(stream.find("MuTau")!=string::npos)
     lOTree->Branch("muFlag"      ,&lMuFlag         ,"muFlag/I");
   else if(stream.find("EleTau")!=string::npos)
     lOTree->Branch("elecFlag"    ,&lMuFlag         ,"elecFlag/I");
   int lTriLepVeto = 0;
   lOTree->Branch("triLepVeto"    ,&lTriLepVeto         ,"triLepVeto/I");

   int   lPairIndex   = 0; lOTree->Branch("pairIndex"   ,&lPairIndex      ,"pairIndex/I");
   
   int   lHLTx        = 0; lOTree->Branch("HLTx"        ,&lHLTx ,        "HLTx/I");
   int   lHLTmatch1    = 0; lOTree->Branch("HLTmatch1"    ,&lHLTmatch1 ,     "HLTmatch1/I");
   int   lHLTmatch2    = 0; lOTree->Branch("HLTmatch2"    ,&lHLTmatch2 ,     "HLTmatch2/I");
   float ldiTauCharge = 0; lOTree->Branch("diTauCharge"      ,&ldiTauCharge, "diTauCharge/F");


   /////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////


   //TFile *file = new TFile(("./batch/nTuple"+sample+"_Open_"+stream+"Stream.root").c_str(),"READ");
   //TTree *tree_orig = (TTree*)file->Get("outTreePtOrd")->Clone("tree");
   TChain *tree = new TChain("outTree");
   tree->Add(("/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74/nTuple"+sample+".root").c_str());

   //TCut sbinPair("etaL1<999.");
   //if( stream.find("TauTau")!=string::npos)
   //sbinPair = TCut("ptL1>45 && ptL2>45 && tightestHPSMVA3oldDMwLTWPL1>2 && tightestHPSMVA3oldDMwLTWPL2>2 && genVMass>(0.7*130) && genVMass<(1.3*130) && HLTx"); // && HLTmatchL1 && HLTmatchL2 && HLTx");

   int iPairIndex      ; tree->SetBranchAddress("pairIndex"     ,&iPairIndex);
   int  iHLTx          ; tree->SetBranchAddress("HLTx"        ,&iHLTx   );
   int  iHLTmatch1      ; tree->SetBranchAddress("HLTmatchL1"    ,&iHLTmatch1);
   int  iHLTmatch2      ; tree->SetBranchAddress("HLTmatchL2"    ,&iHLTmatch2);

   //Bookeeping
   ULong64_t   iRun         ; tree->SetBranchAddress("run"        ,&iRun   );//Run
   ULong64_t   iLumi        ; tree->SetBranchAddress("lumi"       ,&iLumi    );//Lumi
   ULong64_t   iEvt         ; tree->SetBranchAddress("event"      ,&iEvt      );//Evt

   //Event Variables
   float iNPV         ; tree->SetBranchAddress("numPV"        ,&iNPV       );//NPV
   float iNPU         ; tree->SetBranchAddress("npu"  ,&iNPU       );//NPU
   float iRho         ; tree->SetBranchAddress("rho"   ,&iRho         );//Rho
   
   //Event Weights
   //float iMCWeight1    ; tree->SetBranchAddress("sampleWeight"   ,&iMCWeight1    );//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
   //float iMCWeight2    ; tree->SetBranchAddress("HqTWeight"      ,&iMCWeight2    );//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
   //float iPUWeight    ; tree->SetBranchAddress("puWeight"        ,&iPUWeight    );//Pielup Weight
   ////SV Fit variables
   float iMSV         ; tree->SetBranchAddress("diTauNSVfitMass"       ,&iMSV            );//SV Fit using integration method
   float iMSVUp       ; tree->SetBranchAddress("diTauNSVfitMassErrUp"    ,&iMSVUp         );//High Energy scale shape
   float iMSVDown     ; tree->SetBranchAddress("diTauNSVfitMassErrDown"  ,&iMSVDown       );//Low Energy Scale Shape
   float iPtSV        ; tree->SetBranchAddress("diTauNSVfitPt"        ,&iPtSV);
   float iEtaSV       ; tree->SetBranchAddress("diTauNSVfitEta"       ,&iEtaSV);
   float iPhiSV       ; tree->SetBranchAddress("diTauNSVfitPhi"       ,&iPhiSV);
   float iPtH         ; tree->SetBranchAddress("diTauRecoPt"         ,&iPtH);
   float iPtVis       ; tree->SetBranchAddress("diTauVisPt"         ,&iPtVis);
   
   ///First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
   float iPt1         ; tree->SetBranchAddress("ptL1"       ,&iPt1            ); //pT 
   float iPhi1        ; tree->SetBranchAddress("phiL1"      ,&iPhi1           ); //Phi 
   float iEta1        ; tree->SetBranchAddress("etaL1"      ,&iEta1           ); //Eta 
   float iDZ1         ; tree->SetBranchAddress("dzL1"       ,&iDZ1            );//dZ with respect to primary vertex 
   float iD01         ; tree->SetBranchAddress("dxyL1"      ,&iD01            );//d0 with respect to primary vertex 
   float iMt1         ; tree->SetBranchAddress("MtLeg1"     ,&iMt1            );//mT of  first lepton wrt to PF met
   float iMtMVA1      ; tree->SetBranchAddress("MtLeg1MVA"  ,&iMtMVA1            );//mT of  first lepton wrt to MVA met      
   float iM1          ; tree->SetBranchAddress("visibleTauMassL1"        ,&iM1 );
   float iQ1          ; tree->SetBranchAddress("chargeL1"        ,&iQ1 );
   int idecayModeFinding1   ; tree->SetBranchAddress("decayModeFindingL1",   &idecayModeFinding1 );
   int idecayModeFindingNewDM1  ; tree->SetBranchAddress("decayModeFindingNewDML1",  &idecayModeFindingNewDM1 );
   int itightestAntiEMVA5WP1   ; tree->SetBranchAddress("tightestAntiEMVA5WPL1",  &itightestAntiEMVA5WP1 );
   int itightestAntiMu3WP1     ; tree->SetBranchAddress("tightestAntiMu3WPL1",   &itightestAntiMu3WP1 );
   float ihpsDB3H1 ; tree->SetBranchAddress("hpsDB3HL1",    &ihpsDB3H1 );
   //int  iPassId1      ; tree->SetBranchAddress("tightestHPSMVA3oldDMwLTWPL1"   ,&iPassId1      ); 

   //float iscEta1;
   //if(stream.find("ElecTau")!=string::npos )
   //tree->SetBranchAddress("scEtaL1"      ,&iscEta1           ); //supercluster Eta

   ///Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
   float iPt2         ; tree->SetBranchAddress("ptL2"       ,&iPt2            );//pT
   float iPhi2        ; tree->SetBranchAddress("phiL2"      ,&iPhi2           );//Phi
   float iEta2        ; tree->SetBranchAddress("etaL2"      ,&iEta2           );//Eta
   float iDZ2         ; tree->SetBranchAddress("dzL2"       ,&iDZ2            );//dZ with respect to primary vertex
   float iD02         ; tree->SetBranchAddress("dxyL2"      ,&iD02            );//d0 with respect to primary vertex
   float iMt2         ; tree->SetBranchAddress("MtLeg2"     ,&iMt2       );//mT of 2nd lepton wrt to PF met
   float iMtMVA2      ; tree->SetBranchAddress("MtLeg2MVA"  ,&iMtMVA2            );//mT of  2nd lepton wrt to MVA met  
   float iM2          ; tree->SetBranchAddress("visibleTauMassL2"        ,&iM2 );//Mass (visible mass for hadronic Tau)
   int idecayModeFinding2   ; tree->SetBranchAddress("decayModeFindingL2",   &idecayModeFinding2 );
   int idecayModeFindingNewDM2  ; tree->SetBranchAddress("decayModeFindingNewDML2",  &idecayModeFindingNewDM2);
   int itightestAntiEMVA5WP2  ; tree->SetBranchAddress("tightestAntiEMVA5WPL2",  &itightestAntiEMVA5WP2 );
   int itightestAntiMu3WP2   ; tree->SetBranchAddress("tightestAntiMu3WPL2",   &itightestAntiMu3WP2 );
   float ihpsDB3H2 ; tree->SetBranchAddress("hpsDB3HL2",    &ihpsDB3H2 );
   //int  iPassId2      ; tree->SetBranchAddress("tightestHPSMVA3oldDMwLTWPL2"   ,&iPassId2      );//Whether it passes id  (not necessarily iso)
   
   //veto leptons
   int iNVetoMuon     ; tree->SetBranchAddress("nVetoMuon",     &iNVetoMuon );
   int iNVetoElectron ; tree->SetBranchAddress("nVetoElectron",      &iNVetoElectron );

   float iMet         ; tree->SetBranchAddress("MEt"     ,&iMet              ); //pfmet
   float iMetPhi      ; tree->SetBranchAddress("MEtPhi"  ,&iMetPhi           ); //pfmet phi
   float iMVAMet      ; tree->SetBranchAddress("MEtMVA"     ,&iMVAMet        ); //mvamet
   float iMVAMetPhi   ; tree->SetBranchAddress("MEtMVAPhi"  ,&iMVAMetPhi     ); //mvamet Phi
   //MET covariance matrices
   float iMetCov00    ; tree->SetBranchAddress("MEtCov00"   ,&iMetCov00     ); //pf met covariance matrix 00 
   float iMetCov01    ; tree->SetBranchAddress("MEtCov01"   ,&iMetCov01     ); //pf met covariance matrix 01 
   float iMetCov10    ; tree->SetBranchAddress("MEtCov10"   ,&iMetCov10     ); //pf met covariance matrix 10 
   float iMetCov11    ; tree->SetBranchAddress("MEtCov11"   ,&iMetCov11     ); //pf met covariance matrix 11 
   float iMetMvaCov00    ; tree->SetBranchAddress("MEtMVACov00"   ,&iMetMvaCov00     ); //mva met covariance matrix 00 
   float iMetMvaCov01    ; tree->SetBranchAddress("MEtMVACov01"   ,&iMetMvaCov01     ); //mva met covariance matrix 01
   float iMetMvaCov10    ; tree->SetBranchAddress("MEtMVACov10"   ,&iMetMvaCov10     ); //mva met covariance matrix 10 
   float iMetMvaCov11    ; tree->SetBranchAddress("MEtMVACov11"   ,&iMetMvaCov11     ); //mva met covariance matrix 11

   //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
   float iJPt1       ; tree->SetBranchAddress("ptj1"      ,&iJPt1            );//Jet Pt after corrections
   float iJEta1      ; tree->SetBranchAddress("etaj1"     ,&iJEta1           );//Jet Eta
   float iJPhi1      ; tree->SetBranchAddress("phij1"     ,&iJPhi1           );//Jet Phi     
   float iJMVA1      ; tree->SetBranchAddress("pumvaj1"     ,&iJMVA1           );//mva pu jet id 
   float iJCsv1      ; tree->SetBranchAddress("csvj1"     ,&iJCsv1           );//csv discriminator
   float iJPtRaw1    ; tree->SetBranchAddress("ptrawj1"   ,&iJPtRaw1         );//raw jet pt

   //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
   float iJPt2       ; tree->SetBranchAddress("ptj2"      ,&iJPt2            );//Jet Pt after corrections
   float iJEta2      ; tree->SetBranchAddress("etaj2"     ,&iJEta2           );//Jet Eta
   float iJPhi2      ; tree->SetBranchAddress("phij2"     ,&iJPhi2           );//Jet Phi
   float iJMVA2      ; tree->SetBranchAddress("pumvaj2"     ,&iJMVA2           );//mva pu jet id
   float iJCsv2      ; tree->SetBranchAddress("csvj2"     ,&iJCsv2           );//csv discriminator
   float iJPtRaw2    ; tree->SetBranchAddress("ptrawj2"   ,&iJPtRaw2         );//raw jet pt

   //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
   float iBTagPt     ; tree->SetBranchAddress("ptB1"           ,&iBTagPt         );//Corrected BTag Pt
   float iBTagEta    ; tree->SetBranchAddress("etaB1"       ,&iBTagEta        );//Btag Eta
   float iBTagPhi    ; tree->SetBranchAddress("phiB1"       ,&iBTagPhi         );//Btag Phi
 
   //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets 
   float iMJJ        ; tree->SetBranchAddress("Mjj"           ,&iMJJ               );//Mass Di Jet system  
   float iJDEta      ; tree->SetBranchAddress("Detajj"          ,&iJDEta            );//|jeta_1-jeta_2| 
   int   iNJetInGap  ; tree->SetBranchAddress("nVetoJets"  ,&iNJetInGap    );//# of Jets between two jets
   float iptVeto     ; tree->SetBranchAddress("ptVeto"        ,&iptVeto       ); //pt of veto jets
   
   //Variables that go into the VBF MVA
   float iJDPhi      ; tree->SetBranchAddress("Dphijj"      ,&iJDPhi          );//Delta Phi between two leading jets
   float iDiJetPt    ; tree->SetBranchAddress("diJetPt"    ,&iDiJetPt       );//Pt of the di jet system
   float iDiJetPhi   ; tree->SetBranchAddress("diJetPhi"   ,&iDiJetPhi      );//Phi of the di jet system
  
   //number of btags passing btag id ( pt > 20 )
   int   iNBTag      ; tree->SetBranchAddress("nJets20BTagged"      ,&iNBTag         );

   //number of jets passing jet id ( pt > 30 )
   int   iNJets      ; tree->SetBranchAddress("nJets30"      ,&iNJets         );
   int   iNJets20    ; tree->SetBranchAddress("nJets20"      ,&iNJets20       );
  
   float iDiTauCharge; tree->SetBranchAddress("diTauCharge"      ,&iDiTauCharge         );
   
   //float igenVMass = 0; tree->SetBranchAddress("genVMass",    &igenVMass);
   /*
   tree->Draw(">>+skim", sbinPair,"entrylist");   
   
   TEntryList *skim = (TEntryList*)gDirectory->Get("skim");   
   int nEntries = skim->GetN();   
   tree->SetEntryList(skim);   
   if(DEBUG) cout << "-- produced skim : " << skim->GetN() << "entries" << endl;

   // Variables for the loop   
   int nReject=0, nAccept=0;   
   int treenum, iEntry, chainEntry, lastRun, lastEvent;   
   treenum = iEntry = chainEntry = lastRun = lastEvent = -1;  
   */
   for (Long64_t i0=0; i0<tree->GetEntries();i0++) {
     //for (Long64_t i0=0; i0<nEntries; i0++) {
     if (i0 % 1000 == 0) std::cout << "--- ... Processing event: " << double(i0) << std::endl;

     tree->GetEntry(i0);
     //iEntry     = skim->GetEntryAndTree(i0, treenum);   
     //chainEntry = iEntry + (tree->GetTreeOffset())[treenum];   
     //tree->GetEntry(chainEntry);
     //
     //if(iRun==lastRun && iEvt==lastEvent) continue;
     //lastRun  = iRun;    
     //lastEvent= iEvt;
     
     //int jPairIndex = iPairIndex; //iPairIndex[5];
     //bool passDiMuMass = (sample.find("Emb") != string::npos) ? igenDiTauMass > 50.0 : true;
     //bool passDiMuTrig = (sample.find("Emb") != string::npos) ? iHLTxMu17Mu8 > 0 : true;

     //if( stream.find("TauTau")!=string::npos   &&   !(iHLTx>0.5 && /*iHLTmatch1>0.5 && iHLTmatch2>0.5 && */ igenVMass>(0.7*130) && igenVMass<(1.3*130) && iPt1>45 && iPt2>45 && iPassId1>2 && iPassId2>2 && iDZ1 < 0.2 && iDZ2<0.2))
     if(stream.find("TauTau")!=string::npos   &&   !(iHLTx>0.5 && iHLTmatch1>0.5 && iHLTmatch2>0.5 && iPt1>45 && iPt2>45))
       continue;
     if(iPairIndex > 0) continue;
     if(incl && !(iDiTauCharge < 0)) continue; //inclusive selection

     //cout<<iRun<<"  "<<iLumi<<"  "<<iEvt<<endl;
     //cout << iMVA2 << endl;

     lPairIndex  = iPairIndex; //jPairIndex;
     //lMuFlag     = iMuFlag;
     //lTriLepVeto = ivetoEvent;

     lHLTx       = int(iHLTx);
     lHLTmatch1   = int(iHLTmatch1);
     lHLTmatch2   = int(iHLTmatch2);

     //Bookeeping
     lRun         = iRun ;
     lLumi        = iLumi;
     lEvt         = iEvt;

     //Event Variables
     lNPV         = int(iNPV);
     lNPU         = iNPU;
     lRho         = iRho;
   
     //Event Weights
     //lMCWeight    = iMCWeight1*iMCWeight2 ;
     //lPUWeight    = iPUWeight;
     /*if(sample.find("Emb") != string::npos){
       lEffWeight   = stream.find("MuTau")!=string::npos ? iHLTTau*iHLTMu*iSFMuID : iHLTTau*iHLTElec*iSFElecID;
     }
     else
       lEffWeight   = stream.find("MuTau")!=string::npos ? iHLTweightTau*iHLTweightMu*iSFTau*iSFMu : iHLTweightTau*iHLTweightElec*iSFTau*iSFElec;
     lWeight      = lMCWeight*lPUWeight*lEffWeight;
     lTrigweight_1 = stream.find("MuTau")!=string::npos ? iHLTweightMu : iHLTweightElec;
     lTrigweight_2 = iHLTweightTau;
     lIdweight_1  = stream.find("MuTau")!=string::npos ? iSFMuID : iSFElecID;
     lIsoweight_1  = stream.find("MuTau")!=string::npos ? iSFMuIso : iSFElecIso;
     lIdweight_2  = iSFTau;
     lIsoweight_2 = iSFTau;
     */

     //SV Fit variables
     lMSV         = iMSV;
     lMSVUp       = iMSVUp;
     lMSVDown     = iMSVDown;
     lPtSV        = iPtSV;
     lEtaSV       = iEtaSV;
     lPhiSV       = iPhiSV;
     
     ///First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
     lPt1         = iPt1;
     lPhi1        = iPhi1;
     lEta1        = iEta1;
     lM1          = iM1;
     lD01         = iD01;
     lDZ1         = iDZ1;
     lMt1         = iMt1; 
     lMtMVA1      = iMtMVA1;
     lQ1          = iQ1;
     //lPassId1     = stream.find("MuTau")!=string::npos ? (iPassId1L && iPassId1T) : (iPassId1L && ((TMath::Abs(iEta1)<0.80 && iPassId1T_F>0.925) || (TMath::Abs(iEta1)<1.479 && TMath::Abs(iEta1)>0.80 && iPassId1T_F>0.975) || (TMath::Abs(iEta1)>1.479 && iPassId1T_F>0.985)));
     //lPassId1     = iPassId1>2;
     //lPassIso1    = iPassId1>2;
     lagainstElectronLooseMVA5_1  = (itightestAntiEMVA5WP1 >= 2) ? 1 : 0;
     lagainstElectronMediumMVA5_1 = (itightestAntiEMVA5WP1 >= 3) ? 1 : 0;
     lagainstElectronTightMVA5_1  = (itightestAntiEMVA5WP1 >= 4) ? 1 : 0;
     lagainstElectronVLooseMVA5_1 = (itightestAntiEMVA5WP1 >= 1) ? 1 : 0;
     lagainstMuonLoose3_1 = (itightestAntiMu3WP1 >= 1) ? 1 : 0;
     lagainstMuonTight3_1 = (itightestAntiMu3WP1 >= 2) ? 1 : 0;
     lbyCombinedIsolationDeltaBetaCorrRaw3Hits_1 = ihpsDB3H1;
     ldecayModeFinding_1 = idecayModeFinding1;
     ldecayModeFindingNewDMs_1 = idecayModeFindingNewDM1;
     
     ///Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
     lPt2         = iPt2;
     lPhi2        = iPhi2;
     lEta2        = iEta2;
     lM2          = iM2;
     lMtMVA2      = iMtMVA2;
     lD02         = iD02;
     lDZ2         = iDZ2;
     lMt2         = iMt2; 
     lQ2          = iDiTauCharge/iQ1;
     //lIso2        = iIso2;
     //lMVA2        = iMVA2;
     //lPassId2     = iPassId2>2;
     //lPassIso2    = iPassId2>2;
     lagainstElectronLooseMVA5_2  = (itightestAntiEMVA5WP2 >= 2)? 1 : 0;
     lagainstElectronMediumMVA5_2 = (itightestAntiEMVA5WP2 >= 3) ? 1 : 0;
     lagainstElectronTightMVA5_2 = (itightestAntiEMVA5WP2 >= 4) ? 1 : 0;
     lagainstElectronVLooseMVA5_2 = (itightestAntiEMVA5WP2 >= 1) ? 1 : 0;
     lagainstMuonLoose3_2 = (itightestAntiMu3WP2 >= 1) ? 1 : 0;
     lagainstMuonTight3_2 = (itightestAntiMu3WP2 >= 2) ? 1 : 0;
     lbyCombinedIsolationDeltaBetaCorrRaw3Hits_2 = ihpsDB3H2;
     ldecayModeFinding_2 = idecayModeFinding2;
     ldecayModeFindingNewDMs_2 = idecayModeFindingNewDM2;

     ldiTauCharge = iDiTauCharge;

     //veto lepton
     lNVetoMuon = iNVetoMuon;
     lNVetoElectron = iNVetoElectron;
     lExtraelec_veto = (iNVetoElectron > 0);
     lExtramuon_veto = (iNVetoMuon > 0);
     
     //Met related variables
     lMet         = iMet;
     lMetPhi      = iMetPhi;
     lMVAMet      = iMVAMet;
     lMVAMetPhi   = iMVAMetPhi;
     //lPZetaVis    = iPZetaVis;
     //lPZetaMiss   = iPZetaMiss-iPZetaVis;
     //Met covariance matrices
     lMetCov00   = iMetCov00;
     lMetCov01   = iMetCov01;
     lMetCov10   = iMetCov10;
     lMetCov11   = iMetCov11;
     lMVACov00   = iMetMvaCov00;
     lMVACov01   = iMetMvaCov01;
     lMVACov10   = iMetMvaCov10;
     lMVACov11   = iMetMvaCov11;


     //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
     lJPt1       = iJPt1;
     lJEta1      = iJEta1;
     lJPhi1      = iJPhi1;
     lJPtRaw1    = iJPtRaw1;
     lJPtUnc1    = -99;
     lJMVA1      = iJMVA1;
     lJPass1     = true; //iJPass1>0.5;
     lJCsv1      = iJCsv1;

     //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
     lJPt2       = iJPt2;
     lJEta2      = iJEta2;
     lJPhi2      = iJPhi2;
     lJPtRaw2    = iJPtRaw2;
     lJPtUnc2    = -99;
     lJMVA2      = iJMVA2;
     lJPass2     = true; //iJPass2>0.5;
     lJCsv2      = iJCsv2;

     //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
     lBTagPt     = iBTagPt;
     lBTagEta    = iBTagEta;
     lBTagPhi    = iBTagPhi;
 
     //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets 
     lMJJ        = iMJJ;
     lJDEta      = iJDEta;
     lNJetInGap  = (iptVeto > 30) ? iNJetInGap : 0;
     //lMVA        = iMVA;
   
     //Variables that go to the VBF MVA
     lJDPhi      = iJDPhi;
     lDiJetPt    = iDiJetPt;
     lDiJetPhi   = iDiJetPhi;
     //lHDJetPhi   = iHDJetPhi;
     //lVisJetEta  = iVisJetEta;
     //lPtVis      = iPtVis;
  
     //number of btags passing btag id ( pt > 20 )
     lNBTag      = iNBTag;

     //number of jets passing jet id ( pt > 30 )
     lNJets      = iNJets;
     lNJetsPt20    = iNJets20;

     /*
       luPerp    = iuPerp;
       luParl    = iuParl;
       lmetParl    = imetParl;
       lmetPerp    = imetPerp;
       lmetSigmaParl    = imetSigmaParl;
       lmetSigmaPerp    = imetSigmaPerp;
       lmetPullParl    = (imetSigmaParl > 0) ? imetParl/imetSigmaParl : -999.; 
       lmetPullPerp    = (imetSigmaPerp > 0) ? imetPerp/imetSigmaPerp : -999.; 

       lembeddedWeight = iembeddingWeight;
     */
     lOTree->Fill();
     if (i0 % 10000 == 0) std::cout << "--- ... Filling event: " << double(i0) << std::endl;
   }  

   lOFile->cd();
   lOTree->Write();
   lOFile->Close();

   //file->Close();
   //skim->Reset();
   delete tree;

}

void synchNtupleAll()
{
  //synchNtuple("GGFH125_TauTau_nominal", "TauTau", false);
  //synchNtuple("VBFH125_TauTau_nominal", "TauTau", false);
  synchNtuple("SUSYGGH160_TauTau_nominal", "TauTau", false); 
}

