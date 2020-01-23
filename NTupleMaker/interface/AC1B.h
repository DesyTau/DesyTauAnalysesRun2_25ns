//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 27 12:51:14 2018 by ROOT version 6.10/09
// from TTree AC1B/AC1B
// found on file: output_MC.root
//////////////////////////////////////////////////////////

#ifndef AC1B_h
#define AC1B_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "map"

class AC1B {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          errors;
   ULong64_t       event_nr;
   UInt_t          event_run;
   UInt_t          event_timeunix;
   UInt_t          event_timemicrosec;
   UInt_t          event_luminosityblock;
   UChar_t         trigger_level1bits[8];
   UChar_t         trigger_level1[128];
   UChar_t         trigger_HLT[128];
   Float_t         rho;
   UInt_t          primvertex_count;
   UInt_t          goodprimvertex_count;
   Float_t         primvertex_x;
   Float_t         primvertex_y;
   Float_t         primvertex_z;
   Float_t         primvertex_chi2;
   Float_t         primvertex_ndof;
   Float_t         primvertex_ptq;
   Int_t           primvertex_ntracks;
   Float_t         primvertex_cov[6];
   
   UInt_t          primvertexwithbs_count;
   UInt_t          goodprimvertexwithbs_count;
   Float_t         primvertexwithbs_x;
   Float_t         primvertexwithbs_y;
   Float_t         primvertexwithbs_z;
   Float_t         primvertexwithbs_chi2;
   Float_t         primvertexwithbs_ndof;
   Float_t         primvertexwithbs_ptq;
   Int_t           primvertexwithbs_ntracks;
   Float_t         primvertexwithbs_cov[6];   
   

   //Declaration of refitted vertices
   UInt_t          refitvertex_count;
   //UInt_t          goodrefitvertex_count;
   Float_t         refitvertex_x[100];
   Float_t         refitvertex_y[100];
   Float_t         refitvertex_z[100];
   Float_t         refitvertex_chi2[100];
   Float_t         refitvertex_ndof[100];
   Float_t         refitvertex_ptq[100];
   Int_t           refitvertex_ntracks[100];
   Float_t         refitvertex_cov[100][6];
   Int_t   refitvertex_eleIndex[1000][2];
   Int_t   refitvertex_muIndex[1000][2];
   Int_t   refitvertex_tauIndex[1000][2];
   
   UInt_t          refitvertexwithbs_count;
   Float_t         refitvertexwithbs_x[100];
   Float_t         refitvertexwithbs_y[100];
   Float_t         refitvertexwithbs_z[100];
   Float_t         refitvertexwithbs_chi2[100];
   Float_t         refitvertexwithbs_ndof[100];
   Float_t         refitvertexwithbs_ptq[100];
   Int_t           refitvertexwithbs_ntracks[100];
   Float_t         refitvertexwithbs_cov[100][6];
   Int_t   refitvertexwithbs_eleIndex[1000][2];
   Int_t   refitvertexwithbs_muIndex[1000][2];
   Int_t   refitvertexwithbs_tauIndex[1000][2];
   
   //...................................   
   UInt_t          muon_count;
   Float_t         muon_px[100];   //[muon_count]
   Float_t         muon_py[100];   //[muon_count]
   Float_t         muon_pz[100];   //[muon_count]
   Float_t         muon_pt[100];   //[muon_count]
   Float_t         muon_eta[100];   //[muon_count]
   Float_t         muon_phi[100];   //[muon_count]
   Float_t         muon_pterror[100];   //[muon_count]
   Float_t         muon_chi2[100];   //[muon_count]
   Float_t         muon_normChi2[100];   //[muon_count]
   Float_t         muon_ndof[100];   //[muon_count]
   Float_t         muon_charge[100];   //[muon_count]
   Float_t         muon_miniISO[100];   //[muon_count]
   Float_t         muon_combQ_chi2LocalPosition[100];   //[muon_count]
   Float_t         muon_combQ_trkKink[100];   //[muon_count]
   Float_t         muon_validFraction[100];   //[muon_count]
   Float_t         muon_segmentComp[100];   //[muon_count]
   UInt_t          muon_nMuonStations[100];   //[muon_count]
   UInt_t          muon_nMuonHits[100];   //[muon_count]
   UInt_t          muon_nPixelHits[100];   //[muon_count]
   UInt_t          muon_nTrackerHits[100];   //[muon_count]
   Float_t         muon_dxy[100];   //[muon_count]
   Float_t         muon_dxyerr[100];   //[muon_count]
   Float_t         muon_dz[100];   //[muon_count]
   Float_t         muon_dzerr[100];   //[muon_count]
   Float_t         muon_vx[100];   //[muon_count]
   Float_t         muon_vy[100];   //[muon_count]
   Float_t         muon_vz[100];   //[muon_count]
   Float_t         muon_chargedHadIso[100];   //[muon_count]
   Float_t         muon_neutralHadIso[100];   //[muon_count]
   Float_t         muon_photonIso[100];   //[muon_count]
   Float_t         muon_puIso[100];   //[muon_count]
   Float_t         muon_r03_sumChargedHadronPt[100];   //[muon_count]
   Float_t         muon_r03_sumChargedParticlePt[100];   //[muon_count]
   Float_t         muon_r03_sumNeutralHadronEt[100];   //[muon_count]
   Float_t         muon_r03_sumPhotonEt[100];   //[muon_count]
   Float_t         muon_r03_sumNeutralHadronEtHighThreshold[100];   //[muon_count]
   Float_t         muon_r03_sumPhotonEtHighThreshold[100];   //[muon_count]
   Float_t         muon_r03_sumPUPt[100];   //[muon_count]
   Float_t         muon_r04_sumChargedHadronPt[100];   //[muon_count]
   Float_t         muon_r04_sumChargedParticlePt[100];   //[muon_count]
   Float_t         muon_r04_sumNeutralHadronEt[100];   //[muon_count]
   Float_t         muon_r04_sumPhotonEt[100];   //[muon_count]
   Float_t         muon_r04_sumNeutralHadronEtHighThreshold[100];   //[muon_count]
   Float_t         muon_r04_sumPhotonEtHighThreshold[100];   //[muon_count]
   Float_t         muon_r04_sumPUPt[100];   //[muon_count]
   Bool_t          muon_isPF[100];   //[muon_count]
   Bool_t          muon_isGlobal[100];   //[muon_count]
   Bool_t          muon_isTracker[100];   //[muon_count]
   Bool_t          muon_isTight[100];   //[muon_count]
   Bool_t          muon_isLoose[100];   //[muon_count]
   Bool_t          muon_isMedium[100];   //[muon_count]
   Bool_t          muon_isICHEP[100];   //[muon_count]
   Int_t           muon_genmatch[100];   //[muon_count]
   Bool_t          muon_isDuplicate[100];   //[muon_count]
   Bool_t          muon_isBad[100];   //[muon_count]
   Float_t         muon_helixparameters[100][5];
   Float_t         muon_helixparameters_covar[100][5][5];
   Float_t         muon_referencePoint[100][3];
   Float_t         muon_Bfield[100];
   UInt_t          dimuon_count;
   UInt_t          dimuon_leading[100*99/2];   //[dimuon_count]
   UInt_t          dimuon_trailing[100*99/2];   //[dimuon_count]
   Float_t         dimuon_dist2D[100*99/2];   //[dimuon_count]
   Float_t         dimuon_dist2DE[100*99/2];   //[dimuon_count]
   Float_t         dimuon_dist3D[100*99/2];   //[dimuon_count]
   Float_t         dimuon_dist3DE[100*99/2];   //[dimuon_count]
   UInt_t          pfjet_count;
   Float_t         pfjet_e[200];   //[pfjet_count]
   Float_t         pfjet_px[200];   //[pfjet_count]
   Float_t         pfjet_py[200];   //[pfjet_count]
   Float_t         pfjet_pz[200];   //[pfjet_count]
   Float_t         pfjet_pt[200];   //[pfjet_count]
   Float_t         pfjet_eta[200];   //[pfjet_count]
   Float_t         pfjet_phi[200];   //[pfjet_count]
   Float_t         pfjet_neutralhadronicenergy[200];   //[pfjet_count]
   Float_t         pfjet_chargedhadronicenergy[200];   //[pfjet_count]
   Float_t         pfjet_neutralemenergy[200];   //[pfjet_count]
   Float_t         pfjet_chargedemenergy[200];   //[pfjet_count]
   Float_t         pfjet_muonenergy[200];   //[pfjet_count]
   Float_t         pfjet_chargedmuonenergy[200];   //[pfjet_count]
   UInt_t          pfjet_chargedmulti[200];   //[pfjet_count]
   UInt_t          pfjet_neutralmulti[200];   //[pfjet_count]
   UInt_t          pfjet_chargedhadronmulti[200];   //[pfjet_count]
   Float_t         pfjet_energycorr[200];   //[pfjet_count]
   Float_t         pfjet_energycorr_l1fastjet[200];   //[pfjet_count]
   Float_t         pfjet_energycorr_l2relative[200];   //[pfjet_count]
   Float_t         pfjet_energycorr_l3absolute[200];   //[pfjet_count]
   Float_t         pfjet_energycorr_l2l3residual[200];   //[pfjet_count]
   Int_t           pfjet_flavour[200];   //[pfjet_count]
   Float_t         pfjet_btag[200][10];   //[pfjet_count]
   Float_t         pfjet_jecUncertainty[200];   //[pfjet_count]
   Bool_t          pfjet_pu_jet_fullId_loose[200];   //[pfjet_count]
   Bool_t          pfjet_pu_jet_fullId_medium[200];   //[pfjet_count]
   Bool_t          pfjet_pu_jet_fullId_tight[200];   //[pfjet_count]
   Float_t         pfjet_pu_jet_fullDisc_mva[200];   //[pfjet_count]
   UInt_t          electron_count;
   Float_t         electron_px[100];   //[electron_count]
   Float_t         electron_py[100];   //[electron_count]
   Float_t         electron_pz[100];   //[electron_count]
   Float_t         electron_pt[100];   //[electron_count]
   Float_t         electron_eta[100];   //[electron_count]
   Float_t         electron_phi[100];   //[electron_count]
   Float_t         electron_trackchi2[100];   //[electron_count]
   Float_t         electron_trackndof[100];   //[electron_count]
   Float_t         electron_outerx[100];   //[electron_count]
   Float_t         electron_outery[100];   //[electron_count]
   Float_t         electron_outerz[100];   //[electron_count]
   Float_t         electron_vx[100];   //[electron_count]
   Float_t         electron_vy[100];   //[electron_count]
   Float_t         electron_vz[100];   //[electron_count]
   Float_t         electron_esuperclusterovertrack[100];   //[electron_count]
   Float_t         electron_eseedclusterovertrack[100];   //[electron_count]
   Float_t         electron_deltaetasuperclustertrack[100];   //[electron_count]
   Float_t         electron_deltaphisuperclustertrack[100];   //[electron_count]
   Float_t         electron_e1x5[100];   //[electron_count]
   Float_t         electron_e2x5[100];   //[electron_count]
   Float_t         electron_e5x5[100];   //[electron_count]
   Float_t         electron_sigmaetaeta[100];   //[electron_count]
   Float_t         electron_sigmaietaieta[100];   //[electron_count]
   Float_t         electron_ehcaloverecal[100];   //[electron_count]
   Float_t         electron_ehcaloverecaldepth1[100];   //[electron_count]
   Float_t         electron_ehcaloverecaldepth2[100];   //[electron_count]
   Float_t         electron_full5x5_sigmaietaieta[100];   //[electron_count]
   Float_t         electron_ooemoop[100];   //[electron_count]
   Float_t         electron_miniISO[100];   //[electron_count]
   Float_t         electron_superclusterEta[100];   //[electron_count]
   Float_t         electron_superclusterPhi[100];   //[electron_count]
   Float_t         electron_superclusterX[100];   //[electron_count]
   Float_t         electron_superclusterY[100];   //[electron_count]
   Float_t         electron_superclusterZ[100];   //[electron_count]
   Float_t         electron_detaInSeed[100];   //[electron_count]
   Float_t         electron_he[100];   //[electron_count]
   Float_t         electron_eaIsolation[100];   //[electron_count]
   Float_t         electron_chargedHadIso[100];   //[electron_count]
   Float_t         electron_neutralHadIso[100];   //[electron_count]
   Float_t         electron_photonIso[100];   //[electron_count]
   Float_t         electron_puIso[100];   //[electron_count]
   Float_t         electron_r03_sumChargedHadronPt[100];   //[electron_count]
   Float_t         electron_r03_sumChargedParticlePt[100];   //[electron_count]
   Float_t         electron_r03_sumNeutralHadronEt[100];   //[electron_count]
   Float_t         electron_r03_sumPhotonEt[100];   //[electron_count]
   Float_t         electron_r03_sumNeutralHadronEtHighThreshold[100];   //[electron_count]
   Float_t         electron_r03_sumPhotonEtHighThreshold[100];   //[electron_count]
   Float_t         electron_r03_sumPUPt[100];   //[electron_count]
   UChar_t         electron_nhits[100];   //[electron_count]
   UChar_t         electron_npixelhits[100];   //[electron_count]
   UChar_t         electron_nmissinghits[100];   //[electron_count]
   UChar_t         electron_nmissinginnerhits[100];   //[electron_count]
   UChar_t         electron_npixellayers[100];   //[electron_count]
   UChar_t         electron_nstriplayers[100];   //[electron_count]
   Float_t         electron_dxy[100];   //[electron_count]
   Float_t         electron_dxyerr[100];   //[electron_count]
   Float_t         electron_dz[100];   //[electron_count]
   Float_t         electron_dzerr[100];   //[electron_count]
   Float_t         electron_convdist[100];   //[electron_count]
   UInt_t          electron_gapinfo[100];   //[electron_count]
   UInt_t          electron_chargeinfo[100];   //[electron_count]
   Float_t         electron_fbrems[100];   //[electron_count]
   Int_t           electron_numbrems[100];   //[electron_count]
   Float_t         electron_charge[100];   //[electron_count]
   Int_t           electron_superclusterindex[100];   //[electron_count]
   UChar_t         electron_info[100];   //[electron_count]
   Float_t         electron_mva_value_nontrig_Spring15_v1[100];   //[electron_count]
   Float_t         electron_mva_value_trig_Spring15_v1[100];   //[electron_count]
   Int_t           electron_mva_category_nontrig_Spring15_v1[100];   //[electron_count]
   Int_t           electron_mva_category_trig_Spring15_v1[100];   //[electron_count]
   Bool_t          electron_mva_wp80_nontrig_Spring15_v1[100];   //[electron_count]
   Bool_t          electron_mva_wp90_nontrig_Spring15_v1[100];   //[electron_count]
   Bool_t          electron_mva_wp80_trig_Spring15_v1[100];   //[electron_count]
   Bool_t          electron_mva_wp90_trig_Spring15_v1[100];   //[electron_count]
   Bool_t          electron_cutId_veto_Spring15[100];   //[electron_count]
   Bool_t          electron_cutId_loose_Spring15[100];   //[electron_count]
   Bool_t          electron_cutId_medium_Spring15[100];   //[electron_count]
   Bool_t          electron_cutId_tight_Spring15[100];   //[electron_count]
   Bool_t          electron_cutId_veto_Summer16[100];   //[electron_count]
   Bool_t          electron_cutId_loose_Summer16[100];   //[electron_count]
   Bool_t          electron_cutId_medium_Summer16[100];   //[electron_count]
   Bool_t          electron_cutId_tight_Summer16[100];   //[electron_count]
   Float_t         electron_mva_value_Spring16_v1[100];   //[electron_count]
   Int_t           electron_mva_category_Spring16_v1[100];   //[electron_count]
   Float_t         electron_mva_wp90_general_Spring16_v1[100];   //[electron_count]
   Float_t         electron_mva_wp80_general_Spring16_v1[100];   //[electron_count]
   Float_t         electron_mva_value_Iso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_value_noIso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_wp90_Iso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_wp80_Iso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_Loose_Iso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_wp90_noIso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_wp90_noIso_Fall17_v2[100];   //[electron_count]
   Float_t         electron_mva_wp80_noIso_Fall17_v1[100];   //[electron_count]
   Float_t         electron_mva_Loose_noIso_Fall17_v1[100];   //[electron_count]
   Bool_t          electron_cutId_veto_Fall17[100];   //[electron_count]
   Bool_t          electron_cutId_loose_Fall17[100];   //[electron_count]
   Bool_t          electron_cutId_medium_Fall17[100];   //[electron_count]
   Bool_t          electron_cutId_tight_Fall17[100];   //[electron_count]
   Bool_t          electron_cutId_veto_Fall17V2[100];   //[electron_count]
   Bool_t          electron_cutId_loose_Fall17V2[100];   //[electron_count]
   Bool_t          electron_cutId_medium_Fall17V2[100];   //[electron_count]
   Bool_t          electron_cutId_tight_Fall17V2[100];   //[electron_count]
   Bool_t          electron_pass_conversion[100];   //[electron_count]
   Int_t           electron_genmatch[100];   //[electron_count]
   Float_t         electron_px_energyscale_up[100]; //[electron_count]
   Float_t         electron_px_energyscale_down[100]; //[electron_count]
   Float_t         electron_py_energyscale_up[100]; //[electron_count]
   Float_t         electron_py_energyscale_down[100]; //[electron_count]
   Float_t         electron_pz_energyscale_up[100]; //[electron_count]
   Float_t         electron_pz_energyscale_down[100]; //[electron_count]
   Float_t         electron_pt_energyscale_up[100]; //[electron_count]
   Float_t         electron_pt_energyscale_down[100]; //[electron_count]
   Float_t         electron_px_energysigma_up[100]; //[electron_count]
   Float_t         electron_px_energysigma_down[100]; //[electron_count]
   Float_t         electron_py_energysigma_up[100]; //[electron_count]
   Float_t         electron_py_energysigma_down[100]; //[electron_count]
   Float_t         electron_pz_energysigma_up[100]; //[electron_count]
   Float_t         electron_pz_energysigma_down[100]; //[electron_count]
   Float_t         electron_pt_energysigma_up[100]; //[electron_count]
   Float_t         electron_pt_energysigma_down[100]; //[electron_count]
   UInt_t          tau_count;
   Float_t         tau_e[100];   //[tau_count]
   Float_t         tau_px[100];   //[tau_count]
   Float_t         tau_py[100];   //[tau_count]
   Float_t         tau_pz[100];   //[tau_count]
   Float_t         tau_mass[100];   //[tau_count]
   Float_t         tau_eta[100];   //[tau_count]
   Float_t         tau_phi[100];   //[tau_count]
   Float_t         tau_pt[100];   //[tau_count]
   Float_t         tau_vertexx[100];   //[tau_count]
   Float_t         tau_vertexy[100];   //[tau_count]
   Float_t         tau_vertexz[100];   //[tau_count]
   Float_t         tau_pca2D_x[100];   //[tau_count]
   Float_t         tau_pca2D_y[100];   //[tau_count]
   Float_t         tau_pca2D_z[100];   //[tau_count]
   Float_t         tau_pca3D_x[100];   //[tau_count]
   Float_t         tau_pca3D_y[100];   //[tau_count]
   Float_t         tau_pca3D_z[100];   //[tau_count]
   Float_t         tau_SV_x[100];   //[tau_count]
   Float_t         tau_SV_y[100];   //[tau_count]
   Float_t         tau_SV_z[100];   //[tau_count]
   Float_t         tau_SV_cov[100][6]; //[tau_count]
   Float_t         tau_dxy[100];   //[tau_count]
   Float_t         tau_dz[100];   //[tau_count]
   Float_t         tau_ip3d[100];   //[tau_count]
   Float_t         tau_ip3dSig[100];   //[tau_count]
   Float_t         tau_charge[100];   //[tau_count]
   Float_t         tau_genjet_px[100];   //[tau_count]
   Float_t         tau_genjet_py[100];   //[tau_count]
   Float_t         tau_genjet_pz[100];   //[tau_count]
   Float_t         tau_genjet_e[100];   //[tau_count]
   Int_t           tau_genmatch[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_px[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_py[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_pz[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_mass[100];   //[tau_count]
   Int_t           tau_leadchargedhadrcand_id[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_dxy[100];   //[tau_count]
   Float_t         tau_leadchargedhadrcand_dz[100];   //[tau_count]
   UInt_t          tau_ntracks_pt05[100];   //[tau_count]
   UInt_t          tau_ntracks_pt08[100];   //[tau_count]
   UInt_t          tau_ntracks_pt1[100];   //[tau_count]
   Bool_t          tau_L1trigger_match[100];   //[tau_count]
   UInt_t          tau_signalChargedHadrCands_size[100];   //[tau_count]
   UInt_t          tau_signalNeutralHadrCands_size[100];   //[tau_count]
   UInt_t          tau_signalGammaCands_size[100];   //[tau_count]
   UInt_t          tau_isolationChargedHadrCands_size[100];   //[tau_count]
   UInt_t          tau_isolationNeutralHadrCands_size[100];   //[tau_count]
   UInt_t          tau_isolationGammaCands_size[100];   //[tau_count]
   Char_t          tau_genDecayMode_name[100];   //[tau_count]
   Int_t           tau_genDecayMode[100];   //[tau_count]
   Char_t          tau_decayMode_name[100];   //[tau_count]
   Int_t           tau_decayMode[100];   //[tau_count]
   UInt_t          tau_constituents_count[100];   //[tau_count]
   Float_t         tau_constituents_px[100][50];   //[tau_count]
   Float_t         tau_constituents_py[100][50];   //[tau_count]
   Float_t         tau_constituents_pz[100][50];   //[tau_count]
   Float_t         tau_constituents_e[100][50];   //[tau_count]
   Float_t         tau_constituents_mass[100][50];   //[tau_count]
   Int_t           tau_constituents_charge[100][50];   //[tau_count]
   Float_t         tau_constituents_vx[100][50];   //[tau_count]
   Float_t         tau_constituents_vy[100][50];   //[tau_count]
   Float_t         tau_constituents_vz[100][50];   //[tau_count]
   Int_t           tau_constituents_pdgId[100][50];   //[tau_count]
   Float_t         tau_helixparameters[100][5];
   Float_t         tau_helixparameters_covar[100][5][5];
   Float_t         tau_referencePoint[100][3];
   Float_t         tau_Bfield[100];
   UInt_t          track_count;
   Float_t         track_px[1000];   //[track_count]
   Float_t         track_py[1000];   //[track_count]
   Float_t         track_pz[1000];   //[track_count]
   Float_t         track_pt[1000];   //[track_count]
   Float_t         track_eta[1000];   //[track_count]
   Float_t         track_phi[1000];   //[track_count]
   Float_t         track_charge[1000];   //[track_count]
   Float_t         track_mass[1000];   //[track_count]
   Float_t         track_dxy[1000];   //[track_count]
   Float_t         track_dxyerr[1000];   //[track_count]
   Float_t         track_dz[1000];   //[track_count]
   Float_t         track_dzerr[1000];   //[track_count]
   Float_t         track_vx[1000];   //[track_count]
   Float_t         track_vy[1000];   //[track_count]
   Float_t         track_vz[1000];   //[track_count]
   Int_t           track_ID[1000];   //[track_count]
   Bool_t          track_highPurity[1000];   //[track_count]
   Float_t         pfmet_ex;
   Float_t         pfmet_ey;
   Float_t         pfmet_ez;
   Float_t         pfmet_pt;
   Float_t         pfmet_phi;
   Float_t         pfmet_sigxx;
   Float_t         pfmet_sigxy;
   Float_t         pfmet_sigyx;
   Float_t         pfmet_sigyy;
   Float_t         pfmet_sig;
   Float_t         genmet_ex;
   Float_t         genmet_ey;
   Float_t         pfmet_ex_JetEnUp;
   Float_t         pfmet_ey_JetEnUp;
   Float_t         pfmet_ex_JetEnDown;
   Float_t         pfmet_ey_JetEnDown;
   Float_t         pfmet_ex_UnclusteredEnUp;
   Float_t         pfmet_ey_UnclusteredEnUp;
   Float_t         pfmet_ex_UnclusteredEnDown;
   Float_t         pfmet_ey_UnclusteredEnDown;
   Float_t         pfmetcorr_ex;
   Float_t         pfmetcorr_ey;
   Float_t         pfmetcorr_ez;
   Float_t         pfmetcorr_pt;
   Float_t         pfmetcorr_phi;
   Float_t         pfmetcorr_sigxx;
   Float_t         pfmetcorr_sigxy;
   Float_t         pfmetcorr_sigyx;
   Float_t         pfmetcorr_sigyy;
   Float_t         pfmetcorr_sig;
   Float_t         pfmetcorr_ex_JetEnUp;
   Float_t         pfmetcorr_ey_JetEnUp;
   Float_t         pfmetcorr_ex_JetEnDown;
   Float_t         pfmetcorr_ey_JetEnDown;
   Float_t         pfmetcorr_ex_UnclusteredEnUp;
   Float_t         pfmetcorr_ey_UnclusteredEnUp;
   Float_t         pfmetcorr_ex_UnclusteredEnDown;
   Float_t         pfmetcorr_ey_UnclusteredEnDown;
   Float_t         pfmetcorr_ex_JetResUp;
   Float_t         pfmetcorr_ey_JetResUp;
   Float_t         pfmetcorr_ex_JetResDown;
   Float_t         pfmetcorr_ey_JetResDown;
   Float_t         pfmetcorr_ex_smeared;
   Float_t         pfmetcorr_ey_smeared;
   Float_t         pfmetcorr_ez_smeared;
   Float_t         pfmetcorr_pt_smeared;
   Float_t         pfmetcorr_phi_smeared;
   Float_t         pfmetcorr_ex_JetEnUp_smeared;
   Float_t         pfmetcorr_ey_JetEnUp_smeared;
   Float_t         pfmetcorr_ex_JetEnDown_smeared;
   Float_t         pfmetcorr_ey_JetEnDown_smeared;
   Float_t         pfmetcorr_ex_UnclusteredEnUp_smeared;
   Float_t         pfmetcorr_ey_UnclusteredEnUp_smeared;
   Float_t         pfmetcorr_ex_UnclusteredEnDown_smeared;
   Float_t         pfmetcorr_ey_UnclusteredEnDown_smeared;
   Float_t         pfmetcorr_ex_JetResUp_smeared;
   Float_t         pfmetcorr_ey_JetResUp_smeared;
   Float_t         pfmetcorr_ex_JetResDown_smeared;
   Float_t         pfmetcorr_ey_JetResDown_smeared;
   Float_t         puppimet_ex;
   Float_t         puppimet_ey;
   Float_t         puppimet_ez;
   Float_t         puppimet_pt;
   Float_t         puppimet_phi;
   Float_t         puppimet_sigxx;
   Float_t         puppimet_sigxy;
   Float_t         puppimet_sigyx;
   Float_t         puppimet_sigyy;
   Float_t         puppimet_ex_JetEnUp;
   Float_t         puppimet_ey_JetEnUp;
   Float_t         puppimet_ex_JetEnDown;
   Float_t         puppimet_ey_JetEnDown;
   Float_t         puppimet_ex_UnclusteredEnUp;
   Float_t         puppimet_ey_UnclusteredEnUp;
   Float_t         puppimet_ex_UnclusteredEnDown;
   Float_t         puppimet_ey_UnclusteredEnDown;
   Float_t         puppimet_ex_JetResUp;
   Float_t         puppimet_ey_JetResUp;
   Float_t         puppimet_ex_JetResDown;
   Float_t         puppimet_ey_JetResDown;
   Float_t         genweight;
   Float_t         genid1;
   Float_t         genx1;
   Float_t         genid2;
   Float_t         genx2;
   Float_t         genScale;
   Float_t         weightScale0;
   Float_t         weightScale1;
   Float_t         weightScale2;
   Float_t         weightScale3;
   Float_t         weightScale4;
   Float_t         weightScale5;
   Float_t         weightScale6;
   Float_t         weightScale7;
   Float_t         weightScale8;
   Float_t         weightPDFmax;
   Float_t         weightPDFmin;
   Float_t         weightPDFmean;
   Float_t         weightPDFup;
   Float_t         weightPDFdown;
   Float_t         weightPDFvar;
   Float_t         prefiringweight;
   Float_t         prefiringweightup;
   Float_t         prefiringweightdown;
   Int_t           numpileupinteractionsminus;
   Int_t           numpileupinteractions;
   Int_t           numpileupinteractionsplus;
   Float_t         numtruepileupinteractions;
   UInt_t          gentau_count;
   Float_t         gentau_e[100];   //[gentau_count]
   Float_t         gentau_charge[100];   //[gentau_count]
   Float_t         gentau_px[100];   //[gentau_count]
   Float_t         gentau_py[100];   //[gentau_count]
   Float_t         gentau_pz[100];   //[gentau_count]
   Float_t         gentau_visible_e[100];   //[gentau_count]
   Float_t         gentau_visible_px[100];   //[gentau_count]
   Float_t         gentau_visible_py[100];   //[gentau_count]
   Float_t         gentau_visible_pz[100];   //[gentau_count]
   Float_t         gentau_visible_pt[100];   //[gentau_count]
   Float_t         gentau_visible_eta[100];   //[gentau_count]
   Float_t         gentau_visible_phi[100];   //[gentau_count]
   Float_t         gentau_visible_mass[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_e[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_px[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_py[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_pz[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_pt[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_eta[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_phi[100];   //[gentau_count]
   Float_t         gentau_visibleNoLep_mass[100];   //[gentau_count]
   Int_t           gentau_status[100];   //[gentau_count]
   Int_t           gentau_fromHardProcess[100];   //[gentau_count]
   Int_t           gentau_fromHardProcessBeforeFSR[100];   //[gentau_count]
   Int_t           gentau_isDecayedLeptonHadron[100];   //[gentau_count]
   Int_t           gentau_isDirectHadronDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isDirectHardProcessTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isDirectPromptTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isDirectTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isFirstCopy[100];   //[gentau_count]
   Int_t           gentau_isHardProcess[100];   //[gentau_count]
   Int_t           gentau_isHardProcessTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isLastCopy[100];   //[gentau_count]
   Int_t           gentau_isLastCopyBeforeFSR[100];   //[gentau_count]
   Int_t           gentau_isPrompt[100];   //[gentau_count]
   Int_t           gentau_isPromptTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_isTauDecayProduct[100];   //[gentau_count]
   Int_t           gentau_decayMode[100];   //[gentau_count]
   Char_t          gentau_decayMode_name[100];   //[gentau_count]
   UChar_t         gentau_mother[100];   //[gentau_count]
   Float_t         genparticles_lheHt;
   Float_t         genparticles_lheWPt;
   UInt_t          genparticles_noutgoing;
   Int_t           genparticles_noutgoing_NLO;
   UInt_t          genparticles_count;
   Float_t         genparticles_e[200];   //[genparticles_count]
   Float_t         genparticles_px[200];   //[genparticles_count]
   Float_t         genparticles_py[200];   //[genparticles_count]
   Float_t         genparticles_pz[200];   //[genparticles_count]
   Float_t         genparticles_vx[200];   //[genparticles_count]
   Float_t         genparticles_vy[200];   //[genparticles_count]
   Float_t         genparticles_vz[200];   //[genparticles_count]
   Int_t           genparticles_pdgid[200];   //[genparticles_count]
   Int_t           genparticles_status[200];   //[genparticles_count]
   UInt_t          genparticles_info[200];   //[genparticles_count]
   Int_t           genparticles_fromHardProcess[200];   //[genparticles_count]
   Int_t           genparticles_fromHardProcessBeforeFSR[200];   //[genparticles_count]
   Int_t           genparticles_isDecayedLeptonHadron[200];   //[genparticles_count]
   Int_t           genparticles_isDirectHadronDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isDirectHardProcessTauDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isDirectPromptTauDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isDirectTauDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isFirstCopy[200];   //[genparticles_count]
   Int_t           genparticles_isHardProcess[200];   //[genparticles_count]
   Int_t           genparticles_isHardProcessTauDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isLastCopy[200];   //[genparticles_count]
   Int_t           genparticles_isLastCopyBeforeFSR[200];   //[genparticles_count]
   Int_t           genparticles_isPrompt[200];   //[genparticles_count]
   Int_t           genparticles_isPromptTauDecayProduct[200];   //[genparticles_count]
   Int_t           genparticles_isTauDecayProduct[200];   //[genparticles_count]
   UChar_t         genparticles_mother[200];   //[genparticles_count]
   UInt_t          genjets_count;
   Float_t         genjets_e[100];   //[genjets_count]
   Float_t         genjets_px[100];   //[genjets_count]
   Float_t         genjets_py[100];   //[genjets_count]
   Float_t         genjets_pz[100];   //[genjets_count]
   Float_t         genjets_pt[100];   //[genjets_count]
   Float_t         genjets_eta[100];   //[genjets_count]
   Float_t         genjets_phi[100];   //[genjets_count]
   Int_t           genjets_pdgid[100];   //[genjets_count]
   Int_t           genjets_status[100];   //[genjets_count]
   Float_t         genjets_em_energy[100];   //[genjets_count]
   Float_t         genjets_had_energy[100];   //[genjets_count]
   Float_t         genjets_invisible_energy[100];   //[genjets_count]
   Float_t         genjets_auxiliary_energy[100];   //[genjets_count]
   UInt_t          l1muon_count;
   Float_t         l1muon_px[100];   //[l1muon_count]
   Float_t         l1muon_py[100];   //[l1muon_count]
   Float_t         l1muon_pz[100];   //[l1muon_count]
   Float_t         l1muon_pt[100];   //[l1muon_count]
   Int_t           l1muon_ipt[100];   //[l1muon_count]
   Int_t           l1muon_eta[100];   //[l1muon_count]
   Int_t           l1muon_phi[100];   //[l1muon_count]
   Int_t           l1muon_qual[100];   //[l1muon_count]
   Int_t           l1muon_iso[100];   //[l1muon_count]
   Int_t           l1muon_charge[100];   //[l1muon_count]
   Int_t           l1muon_chargeValid[100];   //[l1muon_count]
   Int_t           l1muon_muonIndex[100];   //[l1muon_count]
   Int_t           l1muon_tag[100];   //[l1muon_count]
   Int_t           l1muon_isoSum[100];   //[l1muon_count]
   Int_t           l1muon_dPhiExtra[100];   //[l1muon_count]
   Int_t           l1muon_dEtaExtra[100];   //[l1muon_count]
   Int_t           l1muon_rank[100];   //[l1muon_count]
   Int_t           l1muon_bx[100];   //[l1muon_count]
   UInt_t          l1egamma_count;
   Float_t         l1egamma_px[100];   //[l1egamma_count]
   Float_t         l1egamma_py[100];   //[l1egamma_count]
   Float_t         l1egamma_pz[100];   //[l1egamma_count]
   Float_t         l1egamma_pt[100];   //[l1egamma_count]
   Int_t           l1egamma_ipt[100];   //[l1egamma_count]
   Int_t           l1egamma_eta[100];   //[l1egamma_count]
   Int_t           l1egamma_phi[100];   //[l1egamma_count]
   Int_t           l1egamma_qual[100];   //[l1egamma_count]
   Int_t           l1egamma_iso[100];   //[l1egamma_count]
   Int_t           l1egamma_towerIEta[100];   //[l1egamma_count]
   Int_t           l1egamma_towerIPhi[100];   //[l1egamma_count]
   Int_t           l1egamma_rawEt[100];   //[l1egamma_count]
   Int_t           l1egamma_isoEt[100];   //[l1egamma_count]
   Int_t           l1egamma_footprintEt[100];   //[l1egamma_count]
   Int_t           l1egamma_nTT[100];   //[l1egamma_count]
   Int_t           l1egamma_shape[100];   //[l1egamma_count]
   Int_t           l1egamma_bx[100];   //[l1egamma_count]
   UInt_t          l1tau_count;
   Float_t         l1tau_px[100];   //[l1tau_count]
   Float_t         l1tau_py[100];   //[l1tau_count]
   Float_t         l1tau_pz[100];   //[l1tau_count]
   Float_t         l1tau_pt[100];   //[l1tau_count]
   Int_t           l1tau_ipt[100];   //[l1tau_count]
   Int_t           l1tau_eta[100];   //[l1tau_count]
   Int_t           l1tau_phi[100];   //[l1tau_count]
   Int_t           l1tau_qual[100];   //[l1tau_count]
   Int_t           l1tau_iso[100];   //[l1tau_count]
   Int_t           l1tau_towerIEta[100];   //[l1tau_count]
   Int_t           l1tau_towerIPhi[100];   //[l1tau_count]
   Int_t           l1tau_rawEt[100];   //[l1tau_count]
   Int_t           l1tau_isoEt[100];   //[l1tau_count]
   Int_t           l1tau_nTT[100];   //[l1tau_count]
   Int_t           l1tau_hasEM[100];   //[l1tau_count]
   Int_t           l1tau_isMerged[100];   //[l1tau_count]
   Int_t           l1tau_bx[100];   //[l1tau_count]
   UInt_t          l1isotau_count;
   Float_t         l1isotau_e[50];   //[l1isotau_count]
   Float_t         l1isotau_px[50];   //[l1isotau_count]
   Float_t         l1isotau_py[50];   //[l1isotau_count]
   Float_t         l1isotau_pz[50];   //[l1isotau_count]
   Float_t         l1isotau_mass[50];   //[l1isotau_count]
   Float_t         l1isotau_eta[50];   //[l1isotau_count]
   Float_t         l1isotau_phi[50];   //[l1isotau_count]
   Float_t         l1isotau_pt[50];   //[l1isotau_count]
   Float_t         l1isotau_charge[50];   //[l1isotau_count]
   Int_t           l1isotau_iso[50];   //[l1isotau_count]
   UInt_t          trigobject_count;
   Float_t         trigobject_px[1000];   //[trigobject_count]
   Float_t         trigobject_py[1000];   //[trigobject_count]
   Float_t         trigobject_pz[1000];   //[trigobject_count]
   Float_t         trigobject_pt[1000];   //[trigobject_count]
   Float_t         trigobject_eta[1000];   //[trigobject_count]
   Float_t         trigobject_phi[1000];   //[trigobject_count]
   Bool_t          trigobject_filters[1000][200];   //[trigobject_count]
   Bool_t          trigobject_isMuon[1000];   //[trigobject_count]
   Bool_t          trigobject_isElectron[1000];   //[trigobject_count]
   Bool_t          trigobject_isTau[1000];   //[trigobject_count]
   Bool_t          trigobject_isJet[1000];   //[trigobject_count]
   Bool_t          trigobject_isMET[1000];   //[trigobject_count]
   std::vector<std::string>  *run_hltnames;
   std::vector<std::string>  *run_hltfilters;
   std::vector<std::string>  *run_hltmufilters;
   std::vector<std::string>  *run_hltelectronfilters;
   std::vector<std::string>  *run_hlttaufilters;
   std::vector<std::string>  *run_hltphotonfilters;
   std::vector<std::string>  *run_hltjetfilters;
   std::vector<std::string>  *run_floattaudiscriminators;
   std::vector<std::string>  *run_binarytaudiscriminators;
   std::vector<std::string>  *run_btagdiscriminators;
   std::map<std::string,int> *hltriggerresults;
   std::map<std::string,int> *hltriggerprescales;
   std::vector<std::string>  *hltriggerresultsV;
   std::map<std::string,int> *flags;
   Float_t         tau_againstElectronLooseMVA6[100];   //[tau_count]
   Float_t         tau_againstElectronMVA6Raw[100];   //[tau_count]
   Float_t         tau_againstElectronMVA6category[100];   //[tau_count]
   Float_t         tau_againstElectronMediumMVA6[100];   //[tau_count]
   Float_t         tau_againstElectronTightMVA6[100];   //[tau_count]
   Float_t         tau_againstElectronVLooseMVA6[100];   //[tau_count]
   Float_t         tau_againstElectronVTightMVA6[100];   //[tau_count]
   Float_t         tau_againstMuonLoose3[100];   //[tau_count]
   Float_t         tau_againstMuonTight3[100];   //[tau_count]
   Float_t         tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1DBdR03oldDMwLTraw[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1DBnewDMwLTraw[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1DBoldDMwLTraw[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1PWdR03oldDMwLTraw[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1PWnewDMwLTraw[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1PWoldDMwLTraw[100];   //[tau_count]
   Float_t         tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_byPhotonPtSumOutsideSignalCone[100];   //[tau_count]
   Float_t         tau_byTightCombinedIsolationDeltaBetaCorr3Hits[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1DBnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1DBoldDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1PWnewDMwLT[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1PWoldDMwLT[100];   //[tau_count]
   Float_t         tau_chargedIsoPtSum[100];   //[tau_count]
   Float_t         tau_chargedIsoPtSumdR03[100];   //[tau_count]
   Float_t         tau_decayModeFinding[100];   //[tau_count]
   Float_t         tau_decayModeFindingNewDMs[100];   //[tau_count]
   Float_t         tau_footprintCorrection[100];   //[tau_count]
   Float_t         tau_footprintCorrectiondR03[100];   //[tau_count]
   Float_t         tau_neutralIsoPtSum[100];   //[tau_count]
   Float_t         tau_neutralIsoPtSumWeight[100];   //[tau_count]
   Float_t         tau_neutralIsoPtSumWeightdR03[100];   //[tau_count]
   Float_t         tau_neutralIsoPtSumdR03[100];   //[tau_count]
   Float_t         tau_photonPtSumOutsideSignalCone[100];   //[tau_count]
   Float_t         tau_photonPtSumOutsideSignalConedR03[100];   //[tau_count]
   Float_t         tau_puCorrPtSum[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2017v1DBoldDMwLTraw2017[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2017v2DBnewDMwLTraw2017[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1DBnewDMwLTraw2016[100];   //[tau_count]
   Float_t         tau_byIsolationMVArun2v1DBoldDMwLTraw2016[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byLooseIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byMediumIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byTightIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVTightIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016[100];   //[tau_count]
   Float_t         tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016[100];   //[tau_count]
   Float_t	   tau_byDeepTau2017v2p1VSeraw[100];	//[tau_count]
   Float_t	   tau_byDeepTau2017v2p1VSjetraw[100];	//[tau_count]
   Float_t	   tau_byDeepTau2017v2p1VSmuraw[100];	//[tau_count]
   Float_t	   tau_byLooseDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byLooseDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byLooseDeepTau2017v2p1VSmu[100];	//[tau_count]
   Float_t	   tau_byMediumDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byMediumDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byMediumDeepTau2017v2p1VSmu[100];	//[tau_count]
   Float_t	   tau_byTightDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byTightDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byTightDeepTau2017v2p1VSmu[100];	//[tau_count]
   Float_t	   tau_byVLooseDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byVLooseDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byVLooseDeepTau2017v2p1VSmu[100];	//[tau_count]
   Float_t	   tau_byVTightDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byVTightDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byVVLooseDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byVVLooseDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byVVTightDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byVVTightDeepTau2017v2p1VSjet[100];	//[tau_count]
   Float_t	   tau_byVVVLooseDeepTau2017v2p1VSe[100];	//[tau_count]
   Float_t	   tau_byVVVLooseDeepTau2017v2p1VSjet[100];	//[tau_count]
   
   Float_t         tau_MVADM2017v1[100];
   Float_t         tau_MVADM2017v1DM0raw[100];
   Float_t         tau_MVADM2017v1DM10raw[100];
   Float_t         tau_MVADM2017v1DM11raw[100];
   Float_t         tau_MVADM2017v1DM1raw[100];
   Float_t         tau_MVADM2017v1DM2raw[100];
   Float_t         tau_MVADM2017v1DMotherraw[100];
   
   UInt_t          TauSpinAngles_count;
   Double_t        TauSpinnerWeight[20];//[TauSpinAngles_count]

   Int_t           htxs_stage0cat;
   Int_t           htxs_stage1cat;
   Int_t           htxs_stage1p1cat;
   Float_t         htxs_higgsPt;
   Int_t           htxs_njets30;

   // List of branches
   TBranch        *b_errors;   //!
   TBranch        *b_event_nr;   //!
   TBranch        *b_event_run;   //!
   TBranch        *b_event_timeunix;   //!
   TBranch        *b_event_timemicrosec;   //!
   TBranch        *b_event_luminosityblock;   //!
   TBranch        *b_trigger_level1bits;   //!
   TBranch        *b_trigger_level1;   //!
   TBranch        *b_trigger_HLT;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_primvertex_count;   //!
   TBranch        *b_goodprimvertex_count;   //!
   TBranch        *b_primvertex_x;   //!
   TBranch        *b_primvertex_y;   //!
   TBranch        *b_primvertex_z;   //!
   TBranch        *b_primvertex_chi2;   //!
   TBranch        *b_primvertex_ndof;   //!
   TBranch        *b_primvertex_pdf;   //!
   TBranch        *b_primvertex_ntracks;   //!
   TBranch        *b_primvertex_cov;   //!
   
   TBranch        *b_primvertexwithbs_count;   //!
   TBranch        *b_goodprimvertexwithbs_count;   //!
   TBranch        *b_primvertexwithbs_x;   //!
   TBranch        *b_primvertexwithbs_y;   //!
   TBranch        *b_primvertexwithbs_z;   //!
   TBranch        *b_primvertexwithbs_chi2;   //!
   TBranch        *b_primvertexwithbs_ndof;   //!
   TBranch        *b_primvertexwithbs_pdf;   //!
   TBranch        *b_primvertexwithbs_ntracks;   //!
   TBranch        *b_primvertexwithbs_cov;   //!
   
   //refitvertix
   TBranch        *b_refitvertex_count;   //!
   TBranch        *b_refitvertex_x;   //!
   TBranch        *b_refitvertex_y;   //!
   TBranch        *b_refitvertex_z;   //!
   TBranch        *b_refitvertex_chi2;   //!
   TBranch        *b_refitvertex_ndof;   //!
   TBranch        *b_refitvertex_pdf;   //!
   TBranch        *b_refitvertex_ntracks;   //!
   //TBranch        *b_refitvertex_cov;   //!
   TBranch        *b_refitvertex_eleIndex;
   TBranch        *b_refitvertex_muIndex;
   TBranch        *b_refitvertex_tauIndex;
   
   TBranch        *b_refitvertexwithbs_count;   //!
   TBranch        *b_refitvertexwithbs_x;   //!
   TBranch        *b_refitvertexwithbs_y;   //!
   TBranch        *b_refitvertexwithbs_z;   //!
   TBranch        *b_refitvertexwithbs_chi2;   //!
   TBranch        *b_refitvertexwithbs_ndof;   //!
   TBranch        *b_refitvertexwithbs_pdf;   //!
   TBranch        *b_refitvertexwithbs_ntracks;   //!
   //TBranch        *b_refitvertexwithbs_cov;   //!
   TBranch        *b_refitvertexwithbs_eleIndex;
   TBranch        *b_refitvertexwithbs_muIndex;
   TBranch        *b_refitvertexwithbs_tauIndex;
   
   TBranch        *b_muon_count;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_pterror;   //!
   TBranch        *b_muon_chi2;   //!
   TBranch        *b_muon_normChi2;   //!
   TBranch        *b_muon_ndof;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_miniISO;   //!
   TBranch        *b_muon_combQ_chi2LocalPosition;   //!
   TBranch        *b_muon_combQ_trkKink;   //!
   TBranch        *b_muon_validFraction;   //!
   TBranch        *b_muon_segmentComp;   //!
   TBranch        *b_muon_nMuonStations;   //!
   TBranch        *b_muon_nMuonHits;   //!
   TBranch        *b_muon_nPixelHits;   //!
   TBranch        *b_muon_nTrackerHits;   //!
   TBranch        *b_muon_dxy;   //!
   TBranch        *b_muon_dxyerr;   //!
   TBranch        *b_muon_dz;   //!
   TBranch        *b_muon_dzerr;   //!
   TBranch        *b_muon_vx;   //!
   TBranch        *b_muon_vy;   //!
   TBranch        *b_muon_vz;   //!
   TBranch        *b_muon_chargedHadIso;   //!
   TBranch        *b_muon_neutralHadIso;   //!
   TBranch        *b_muon_photonIso;   //!
   TBranch        *b_muon_puIso;   //!
   TBranch        *b_muon_r03_sumChargedHadronPt;   //!
   TBranch        *b_muon_r03_sumChargedParticlePt;   //!
   TBranch        *b_muon_r03_sumNeutralHadronEt;   //!
   TBranch        *b_muon_r03_sumPhotonEt;   //!
   TBranch        *b_muon_r03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_muon_r03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_muon_r03_sumPUPt;   //!
   TBranch        *b_muon_r04_sumChargedHadronPt;   //!
   TBranch        *b_muon_r04_sumChargedParticlePt;   //!
   TBranch        *b_muon_r04_sumNeutralHadronEt;   //!
   TBranch        *b_muon_r04_sumPhotonEt;   //!
   TBranch        *b_muon_r04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_muon_r04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_muon_r04_sumPUPt;   //!
   TBranch        *b_muon_isPF;   //!
   TBranch        *b_muon_isGlobal;   //!
   TBranch        *b_muon_isTracker;   //!
   TBranch        *b_muon_isTight;   //!
   TBranch        *b_muon_isLoose;   //!
   TBranch        *b_muon_isMedium;   //!
   TBranch        *b_muon_isICHEP;   //!
   TBranch        *b_muon_genmatch;   //!
   TBranch        *b_muon_isDuplicate;   //!
   TBranch        *b_muon_isBad;   //!
   TBranch        *b_muon_helixparameters;
   TBranch        *b_muon_helixparameters_covar;
   TBranch        *b_muon_referencePoint;
   TBranch        *b_muon_Bfield;
   TBranch        *b_dimuon_count;   //!
   TBranch        *b_dimuon_leading;   //!
   TBranch        *b_dimuon_trailing;   //!
   TBranch        *b_dimuon_dist2D;   //!
   TBranch        *b_dimuon_dist2DE;   //!
   TBranch        *b_dimuon_dist3D;   //!
   TBranch        *b_dimuon_dist3DE;   //!
   TBranch        *b_pfjet_count;   //!
   TBranch        *b_pfjet_e;   //!
   TBranch        *b_pfjet_px;   //!
   TBranch        *b_pfjet_py;   //!
   TBranch        *b_pfjet_pz;   //!
   TBranch        *b_pfjet_pt;   //!
   TBranch        *b_pfjet_eta;   //!
   TBranch        *b_pfjet_phi;   //!
   TBranch        *b_pfjet_neutralhadronicenergy;   //!
   TBranch        *b_pfjet_chargedhadronicenergy;   //!
   TBranch        *b_pfjet_neutralemenergy;   //!
   TBranch        *b_pfjet_chargedemenergy;   //!
   TBranch        *b_pfjet_muonenergy;   //!
   TBranch        *b_pfjet_chargedmuonenergy;   //!
   TBranch        *b_pfjet_chargedmulti;   //!
   TBranch        *b_pfjet_neutralmulti;   //!
   TBranch        *b_pfjet_chargedhadronmulti;   //!
   TBranch        *b_pfjet_energycorr;   //!
   TBranch        *b_pfjet_energycorr_l1fastjet;   //!
   TBranch        *b_pfjet_energycorr_l2relative;   //!
   TBranch        *b_pfjet_energycorr_l3absolute;   //!
   TBranch        *b_pfjet_energycorr_l2l3residual;   //!
   TBranch        *b_pfjet_flavour;   //!
   TBranch        *b_pfjet_btag;   //!
   TBranch        *b_pfjet_jecUncertainty;   //!
   TBranch        *b_pfjet_pu_jet_fullId_loose;   //!
   TBranch        *b_pfjet_pu_jet_fullId_medium;   //!
   TBranch        *b_pfjet_pu_jet_fullId_tight;   //!
   TBranch        *b_pfjet_pu_jet_fullDisc_mva;   //!
   TBranch        *b_electron_count;   //!
   TBranch        *b_electron_px;   //!
   TBranch        *b_electron_py;   //!
   TBranch        *b_electron_pz;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_trackchi2;   //!
   TBranch        *b_electron_trackndof;   //!
   TBranch        *b_electron_outerx;   //!
   TBranch        *b_electron_outery;   //!
   TBranch        *b_electron_outerz;   //!
   TBranch        *b_electron_vx;   //!
   TBranch        *b_electron_vy;   //!
   TBranch        *b_electron_vz;   //!
   TBranch        *b_electron_esuperclusterovertrack;   //!
   TBranch        *b_electron_eseedclusterovertrack;   //!
   TBranch        *b_electron_deltaetasuperclustertrack;   //!
   TBranch        *b_electron_deltaphisuperclustertrack;   //!
   TBranch        *b_electron_e1x5;   //!
   TBranch        *b_electron_e2x5;   //!
   TBranch        *b_electron_e5x5;   //!
   TBranch        *b_electron_sigmaetaeta;   //!
   TBranch        *b_electron_sigmaietaieta;   //!
   TBranch        *b_electron_ehcaloverecal;   //!
   TBranch        *b_electron_ehcaloverecaldepth1;   //!
   TBranch        *b_electron_ehcaloverecaldepth2;   //!
   TBranch        *b_electron_full5x5_sigmaietaieta;   //!
   TBranch        *b_electron_ooemoop;   //!
   TBranch        *b_electron_miniISO;   //!
   TBranch        *b_electron_superclusterEta;   //!
   TBranch        *b_electron_superclusterPhi;   //!
   TBranch        *b_electron_superclusterX;   //!
   TBranch        *b_electron_superclusterY;   //!
   TBranch        *b_electron_superclusterZ;   //!
   TBranch        *b_electron_detaInSeed;   //!
   TBranch        *b_electron_he;   //!
   TBranch        *b_electron_eaIsolation;   //!
   TBranch        *b_electron_chargedHadIso;   //!
   TBranch        *b_electron_neutralHadIso;   //!
   TBranch        *b_electron_photonIso;   //!
   TBranch        *b_electron_puIso;   //!
   TBranch        *b_electron_r03_sumChargedHadronPt;   //!
   TBranch        *b_electron_r03_sumChargedParticlePt;   //!
   TBranch        *b_electron_r03_sumNeutralHadronEt;   //!
   TBranch        *b_electron_r03_sumPhotonEt;   //!
   TBranch        *b_electron_r03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_electron_r03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_electron_r03_sumPUPt;   //!
   TBranch        *b_electron_nhits;   //!
   TBranch        *b_electron_npixelhits;   //!
   TBranch        *b_electron_nmissinghits;   //!
   TBranch        *b_electron_nmissinginnerhits;   //!
   TBranch        *b_electron_npixellayers;   //!
   TBranch        *b_electron_nstriplayers;   //!
   TBranch        *b_electron_dxy;   //!
   TBranch        *b_electron_dxyerr;   //!
   TBranch        *b_electron_dz;   //!
   TBranch        *b_electron_dzerr;   //!
   TBranch        *b_electron_convdist;   //!
   TBranch        *b_electron_gapinfo;   //!
   TBranch        *b_electron_chargeinfo;   //!
   TBranch        *b_electron_fbrems;   //!
   TBranch        *b_electron_numbrems;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_electron_superclusterindex;   //!
   TBranch        *b_electron_info;   //!
   TBranch        *b_electron_mva_value_nontrig_Spring15_v1;   //!
   TBranch        *b_electron_mva_value_trig_Spring15_v1;   //!
   TBranch        *b_electron_mva_category_nontrig_Spring15_v1;   //!
   TBranch        *b_electron_mva_category_trig_Spring15_v1;   //!
   TBranch        *b_electron_mva_wp80_nontrig_Spring15_v1;   //!
   TBranch        *b_electron_mva_wp90_nontrig_Spring15_v1;   //!
   TBranch        *b_electron_mva_wp80_trig_Spring15_v1;   //!
   TBranch        *b_electron_mva_wp90_trig_Spring15_v1;   //!
   TBranch        *b_electron_cutId_veto_Spring15;   //!
   TBranch        *b_electron_cutId_loose_Spring15;   //!
   TBranch        *b_electron_cutId_medium_Spring15;   //!
   TBranch        *b_electron_cutId_tight_Spring15;   //!
   TBranch        *b_electron_cutId_veto_Summer16;   //!
   TBranch        *b_electron_cutId_loose_Summer16;   //!
   TBranch        *b_electron_cutId_medium_Summer16;   //!
   TBranch        *b_electron_cutId_tight_Summer16;   //!
   TBranch        *b_electron_mva_value_Spring16_v1;   //!
   TBranch        *b_electron_mva_category_Spring16_v1;   //!
   TBranch        *b_electron_mva_wp90_general_Spring16_v1;   //!
   TBranch        *b_electron_mva_wp80_general_Spring16_v1;   //!
   TBranch        *b_electron_mva_value_Iso_Fall17_v1;   //!
   TBranch        *b_electron_mva_value_noIso_Fall17_v1;   //!
   TBranch        *b_electron_mva_wp90_Iso_Fall17_v1;   //!
   TBranch        *b_electron_mva_wp80_Iso_Fall17_v1;   //!
   TBranch        *b_electron_mva_Loose_Iso_Fall17_v1;   //!
   TBranch        *b_electron_mva_wp90_noIso_Fall17_v1;   //!
   TBranch        *b_electron_mva_wp90_noIso_Fall17_v2;   //!
   TBranch        *b_electron_mva_wp80_noIso_Fall17_v1;   //!
   TBranch        *b_electron_mva_Loose_noIso_Fall17_v1;   //!
   TBranch        *b_electron_cutId_veto_Fall17;   //!
   TBranch        *b_electron_cutId_loose_Fall17;   //!
   TBranch        *b_electron_cutId_medium_Fall17;   //!
   TBranch        *b_electron_cutId_tight_Fall17;   //!
   TBranch        *b_electron_cutId_veto_Fall17V2;   //!
   TBranch        *b_electron_cutId_loose_Fall17V2;   //!
   TBranch        *b_electron_cutId_medium_Fall17V2;   //!
   TBranch        *b_electron_cutId_tight_Fall17V2;   //!
   TBranch        *b_electron_pass_conversion;   //!
   TBranch        *b_electron_genmatch;   //!;   //!
   TBranch        *b_electron_px_energyscale_up;   //!
   TBranch        *b_electron_px_energyscale_down;   //!
   TBranch        *b_electron_py_energyscale_up;   //!
   TBranch        *b_electron_py_energyscale_down;   //!
   TBranch        *b_electron_pz_energyscale_up;   //!
   TBranch        *b_electron_pz_energyscale_down;   //!
   TBranch        *b_electron_pt_energyscale_up;   //!
   TBranch        *b_electron_pt_energyscale_down;   //!
   TBranch        *b_electron_px_energysigma_up;   //!
   TBranch        *b_electron_px_energysigma_down;   //!
   TBranch        *b_electron_py_energysigma_up;   //!
   TBranch        *b_electron_py_energysigma_down;   //!
   TBranch        *b_electron_pz_energysigma_up;   //!
   TBranch        *b_electron_pz_energysigma_down;   //!
   TBranch        *b_electron_pt_energysigma_up;   //!
   TBranch        *b_electron_pt_energysigma_down;   //!
   TBranch        *b_tau_count;   //!
   TBranch        *b_tau_e;   //!
   TBranch        *b_tau_px;   //!
   TBranch        *b_tau_py;   //!
   TBranch        *b_tau_pz;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_vertexx;   //!
   TBranch        *b_tau_vertexy;   //!
   TBranch        *b_tau_vertexz;   //!
   TBranch        *b_tau_pca2D_x;   //!
   TBranch        *b_tau_pca2D_y;   //!
   TBranch        *b_tau_pca2D_z;   //!
   TBranch        *b_tau_pca3D_x;   //!
   TBranch        *b_tau_pca3D_y;   //!
   TBranch        *b_tau_pca3D_z;   //!
   TBranch        *b_tau_SV_x;   //!
   TBranch        *b_tau_SV_y;   //!
   TBranch        *b_tau_SV_z;   //!
   TBranch        *b_tau_SV_cov;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_tau_ip3d;   //!
   TBranch        *b_tau_ip3dSig;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_genjet_px;   //!
   TBranch        *b_tau_genjet_py;   //!
   TBranch        *b_tau_genjet_pz;   //!
   TBranch        *b_tau_genjet_e;   //!
   TBranch        *b_tau_genmatch;   //!
   TBranch        *b_tau_leadchargedhadrcand_px;   //!
   TBranch        *b_tau_leadchargedhadrcand_py;   //!
   TBranch        *b_tau_leadchargedhadrcand_pz;   //!
   TBranch        *b_tau_leadchargedhadrcand_mass;   //!
   TBranch        *b_tau_leadchargedhadrcand_id;   //!
   TBranch        *b_tau_leadchargedhadrcand_dxy;   //!
   TBranch        *b_tau_leadchargedhadrcand_dz;   //!
   TBranch        *b_tau_ntracks_pt05;   //!
   TBranch        *b_tau_ntracks_pt08;   //!
   TBranch        *b_tau_ntracks_pt1;   //!
   TBranch        *b_tau_L1trigger_match;   //!
   TBranch        *b_tau_signalChargedHadrCands_size;   //!
   TBranch        *b_tau_signalNeutralHadrCands_size;   //!
   TBranch        *b_tau_signalGammaCands_size;   //!
   TBranch        *b_tau_isolationChargedHadrCands_size;   //!
   TBranch        *b_tau_isolationNeutralHadrCands_size;   //!
   TBranch        *b_tau_isolationGammaCands_size;   //!
   TBranch        *b_tau_genDecayMode_name;   //!
   TBranch        *b_tau_genDecayMode;   //!
   TBranch        *b_tau_decayMode_name;   //!
   TBranch        *b_tau_decayMode;   //!
   TBranch        *b_tau_constituents_count;   //!
   TBranch        *b_tau_constituents_px;   //!
   TBranch        *b_tau_constituents_py;   //!
   TBranch        *b_tau_constituents_pz;   //!
   TBranch        *b_tau_constituents_e;   //!
   TBranch        *b_tau_constituents_mass;   //!
   TBranch        *b_tau_constituents_charge;   //!
   TBranch        *b_tau_constituents_vx;   //!
   TBranch        *b_tau_constituents_vy;   //!
   TBranch        *b_tau_constituents_vz;   //!
   TBranch        *b_tau_constituents_pdgId;   //!
   TBranch        *b_tau_helixparameters;
   TBranch        *b_tau_helixparameters_covar;
   TBranch        *b_tau_referencePoint;
   TBranch        *b_tau_Bfield;
   TBranch        *b_track_count;   //!
   TBranch        *b_track_px;   //!
   TBranch        *b_track_py;   //!
   TBranch        *b_track_pz;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_mass;   //!
   TBranch        *b_track_dxy;   //!
   TBranch        *b_track_dxyerr;   //!
   TBranch        *b_track_dz;   //!
   TBranch        *b_track_dzerr;   //!
   TBranch        *b_track_vx;   //!
   TBranch        *b_track_vy;   //!
   TBranch        *b_track_vz;   //!
   TBranch        *b_track_ID;   //!
   TBranch        *b_track_highPurity;   //!
   TBranch        *b_pfmet_ex;   //!
   TBranch        *b_pfmet_ey;   //!
   TBranch        *b_pfmet_ez;   //!
   TBranch        *b_pfmet_pt;   //!
   TBranch        *b_pfmet_phi;   //!
   TBranch        *b_pfmet_sigxx;   //!
   TBranch        *b_pfmet_sigxy;   //!
   TBranch        *b_pfmet_sigyx;   //!
   TBranch        *b_pfmet_sigyy;   //!
   TBranch        *b_pfmet_sig;   //!
   TBranch        *b_genmet_ex;   //!
   TBranch        *b_genmet_ey;   //!
   TBranch        *b_pfmet_ex_JetEnUp;   //!
   TBranch        *b_pfmet_ey_JetEnUp;   //!
   TBranch        *b_pfmet_ex_JetEnDown;   //!
   TBranch        *b_pfmet_ey_JetEnDown;   //!
   TBranch        *b_pfmet_ex_UnclusteredEnUp;   //!
   TBranch        *b_pfmet_ey_UnclusteredEnUp;   //!
   TBranch        *b_pfmet_ex_UnclusteredEnDown;   //!
   TBranch        *b_pfmet_ey_UnclusteredEnDown;   //!
   TBranch        *b_pfmetcorr_ex;   //!
   TBranch        *b_pfmetcorr_ey;   //!
   TBranch        *b_pfmetcorr_ez;   //!
   TBranch        *b_pfmetcorr_pt;   //!
   TBranch        *b_pfmetcorr_phi;   //!
   TBranch        *b_pfmetcorr_sigxx;   //!
   TBranch        *b_pfmetcorr_sigxy;   //!
   TBranch        *b_pfmetcorr_sigyx;   //!
   TBranch        *b_pfmetcorr_sigyy;   //!
   TBranch        *b_pfmetcorr_sig;   //!
   TBranch        *b_pfmetcorr_ex_JetEnUp;   //!
   TBranch        *b_pfmetcorr_ey_JetEnUp;   //!
   TBranch        *b_pfmetcorr_ex_JetEnDown;   //!
   TBranch        *b_pfmetcorr_ey_JetEnDown;   //!
   TBranch        *b_pfmetcorr_ex_UnclusteredEnUp;   //!
   TBranch        *b_pfmetcorr_ey_UnclusteredEnUp;   //!
   TBranch        *b_pfmetcorr_ex_UnclusteredEnDown;   //!
   TBranch        *b_pfmetcorr_ey_UnclusteredEnDown;   //!
   TBranch        *b_pfmetcorr_ex_JetResUp;   //!
   TBranch        *b_pfmetcorr_ey_JetResUp;   //!
   TBranch        *b_pfmetcorr_ex_JetResDown;   //!
   TBranch        *b_pfmetcorr_ey_JetResDown;   //!
   TBranch        *b_pfmetcorr_ex_smeared;   //!
   TBranch        *b_pfmetcorr_ey_smeared;   //!
   TBranch        *b_pfmetcorr_ez_smeared;   //!
   TBranch        *b_pfmetcorr_pt_smeared;   //!
   TBranch        *b_pfmetcorr_phi_smeared;   //!
   TBranch        *b_pfmetcorr_ex_JetEnUp_smeared;   //!
   TBranch        *b_pfmetcorr_ey_JetEnUp_smeared;   //!
   TBranch        *b_pfmetcorr_ex_JetEnDown_smeared;   //!
   TBranch        *b_pfmetcorr_ey_JetEnDown_smeared;   //!
   TBranch        *b_pfmetcorr_ex_UnclusteredEnUp_smeared;   //!
   TBranch        *b_pfmetcorr_ey_UnclusteredEnUp_smeared;   //!
   TBranch        *b_pfmetcorr_ex_UnclusteredEnDown_smeared;   //!
   TBranch        *b_pfmetcorr_ey_UnclusteredEnDown_smeared;   //!
   TBranch        *b_pfmetcorr_ex_JetResUp_smeared;   //!
   TBranch        *b_pfmetcorr_ey_JetResUp_smeared;   //!
   TBranch        *b_pfmetcorr_ex_JetResDown_smeared;   //!
   TBranch        *b_pfmetcorr_ey_JetResDown_smeared;   //!
   TBranch        *b_puppimet_ex;   //!
   TBranch        *b_puppimet_ey;   //!
   TBranch        *b_puppimet_ez;   //!
   TBranch        *b_puppimet_pt;   //!
   TBranch        *b_puppimet_phi;   //!
   TBranch        *b_puppimet_sigxx;   //!
   TBranch        *b_puppimet_sigxy;   //!
   TBranch        *b_puppimet_sigyx;   //!
   TBranch        *b_puppimet_sigyy;   //!
   TBranch        *b_puppimet_ex_JetEnUp;   //!
   TBranch        *b_puppimet_ey_JetEnUp;   //!
   TBranch        *b_puppimet_ex_JetEnDown;   //!
   TBranch        *b_puppimet_ey_JetEnDown;   //!
   TBranch        *b_puppimet_ex_UnclusteredEnUp;   //!
   TBranch        *b_puppimet_ey_UnclusteredEnUp;   //!
   TBranch        *b_puppimet_ex_UnclusteredEnDown;   //!
   TBranch        *b_puppimet_ey_UnclusteredEnDown;   //!
   TBranch        *b_puppimet_ex_JetResUp;   //!
   TBranch        *b_puppimet_ey_JetResUp;   //!
   TBranch        *b_puppimet_ex_JetResDown;   //!
   TBranch        *b_puppimet_ey_JetResDown;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_genid1;   //!
   TBranch        *b_genx1;   //!
   TBranch        *b_genid2;   //!
   TBranch        *b_genx2;   //!
   TBranch        *b_genScale;   //!
   TBranch        *b_weightScale0;   //!
   TBranch        *b_weightScale1;   //!
   TBranch        *b_weightScale2;   //!
   TBranch        *b_weightScale3;   //!
   TBranch        *b_weightScale4;   //!
   TBranch        *b_weightScale5;   //!
   TBranch        *b_weightScale6;   //!
   TBranch        *b_weightScale7;   //!
   TBranch        *b_weightScale8;   //!
   TBranch        *b_weightPDFmax;   //!
   TBranch        *b_weightPDFmin;   //!
   TBranch        *b_weightPDFmean;   //!
   TBranch        *b_weightPDFup;   //!
   TBranch        *b_weightPDFdown;   //!
   TBranch        *b_weightPDFvar;   //!
   TBranch        *b_prefiringweight;
   TBranch        *b_prefiringweightup;
   TBranch        *b_prefiringweightdown;
   TBranch        *b_numpileupinteractionsminus;   //!
   TBranch        *b_numpileupinteractions;   //!
   TBranch        *b_numpileupinteractionsplus;   //!
   TBranch        *b_numtruepileupinteractions;   //!
   TBranch        *b_gentau_count;   //!
   TBranch        *b_gentau_e;   //!
   TBranch        *b_gentau_charge;   //!
   TBranch        *b_gentau_px;   //!
   TBranch        *b_gentau_py;   //!
   TBranch        *b_gentau_pz;   //!
   TBranch        *b_gentau_visible_e;   //!
   TBranch        *b_gentau_visible_px;   //!
   TBranch        *b_gentau_visible_py;   //!
   TBranch        *b_gentau_visible_pz;   //!
   TBranch        *b_gentau_visible_pt;   //!
   TBranch        *b_gentau_visible_eta;   //!
   TBranch        *b_gentau_visible_phi;   //!
   TBranch        *b_gentau_visible_mass;   //!
   TBranch        *b_gentau_visibleNoLep_e;   //!
   TBranch        *b_gentau_visibleNoLep_px;   //!
   TBranch        *b_gentau_visibleNoLep_py;   //!
   TBranch        *b_gentau_visibleNoLep_pz;   //!
   TBranch        *b_gentau_visibleNoLep_pt;   //!
   TBranch        *b_gentau_visibleNoLep_eta;   //!
   TBranch        *b_gentau_visibleNoLep_phi;   //!
   TBranch        *b_gentau_visibleNoLep_mass;   //!
   TBranch        *b_gentau_status;   //!
   TBranch        *b_gentau_fromHardProcess;   //!
   TBranch        *b_gentau_fromHardProcessBeforeFSR;   //!
   TBranch        *b_gentau_isDecayedLeptonHadron;   //!
   TBranch        *b_gentau_isDirectHadronDecayProduct;   //!
   TBranch        *b_gentau_isDirectHardProcessTauDecayProduct;   //!
   TBranch        *b_gentau_isDirectPromptTauDecayProduct;   //!
   TBranch        *b_gentau_isDirectTauDecayProduct;   //!
   TBranch        *b_gentau_isFirstCopy;   //!
   TBranch        *b_gentau_isHardProcess;   //!
   TBranch        *b_gentau_isHardProcessTauDecayProduct;   //!
   TBranch        *b_gentau_isLastCopy;   //!
   TBranch        *b_gentau_isLastCopyBeforeFSR;   //!
   TBranch        *b_gentau_isPrompt;   //!
   TBranch        *b_gentau_isPromptTauDecayProduct;   //!
   TBranch        *b_gentau_isTauDecayProduct;   //!
   TBranch        *b_gentau_decayMode;   //!
   TBranch        *b_gentau_decayMode_name;   //!
   TBranch        *b_gentau_mother;   //!
   TBranch        *b_genparticles_lheHt;   //!
   TBranch        *b_genparticles_lheWPt;   //!
   TBranch        *b_genparticles_noutgoing;   //!
   TBranch        *b_genparticles_noutgoing_NLO;   //!
   TBranch        *b_genparticles_count;   //!
   TBranch        *b_genparticles_e;   //!
   TBranch        *b_genparticles_px;   //!
   TBranch        *b_genparticles_py;   //!
   TBranch        *b_genparticles_pz;   //!
   TBranch        *b_genparticles_vx;   //!
   TBranch        *b_genparticles_vy;   //!
   TBranch        *b_genparticles_vz;   //!
   TBranch        *b_genparticles_pdgid;   //!
   TBranch        *b_genparticles_status;   //!
   TBranch        *b_genparticles_info;   //!
   TBranch        *b_genparticles_fromHardProcess;   //!
   TBranch        *b_genparticles_fromHardProcessBeforeFSR;   //!
   TBranch        *b_genparticles_isDecayedLeptonHadron;   //!
   TBranch        *b_genparticles_isDirectHadronDecayProduct;   //!
   TBranch        *b_genparticles_isDirectHardProcessTauDecayProduct;   //!
   TBranch        *b_genparticles_isDirectPromptTauDecayProduct;   //!
   TBranch        *b_genparticles_isDirectTauDecayProduct;   //!
   TBranch        *b_genparticles_isFirstCopy;   //!
   TBranch        *b_genparticles_isHardProcess;   //!
   TBranch        *b_genparticles_isHardProcessTauDecayProduct;   //!
   TBranch        *b_genparticles_isLastCopy;   //!
   TBranch        *b_genparticles_isLastCopyBeforeFSR;   //!
   TBranch        *b_genparticles_isPrompt;   //!
   TBranch        *b_genparticles_isPromptTauDecayProduct;   //!
   TBranch        *b_genparticles_isTauDecayProduct;   //!
   TBranch        *b_genparticles_mother;   //!
   TBranch        *b_genjets_count;   //!
   TBranch        *b_genjets_e;   //!
   TBranch        *b_genjets_px;   //!
   TBranch        *b_genjets_py;   //!
   TBranch        *b_genjets_pz;   //!
   TBranch        *b_genjets_pt;   //!
   TBranch        *b_genjets_eta;   //!
   TBranch        *b_genjets_phi;   //!
   TBranch        *b_genjets_pdgid;   //!
   TBranch        *b_genjets_status;   //!
   TBranch        *b_genjets_em_energy;   //!
   TBranch        *b_genjets_had_energy;   //!
   TBranch        *b_genjets_invisible_energy;   //!
   TBranch        *b_genjets_auxiliary_energy;   //!
   TBranch        *b_l1muon_count;   //!
   TBranch        *b_l1muon_px;   //!
   TBranch        *b_l1muon_py;   //!
   TBranch        *b_l1muon_pz;   //!
   TBranch        *b_l1muon_pt;   //!
   TBranch        *b_l1muon_ipt;   //!
   TBranch        *b_l1muon_eta;   //!
   TBranch        *b_l1muon_phi;   //!
   TBranch        *b_l1muon_qual;   //!
   TBranch        *b_l1muon_iso;   //!
   TBranch        *b_l1muon_charge;   //!
   TBranch        *b_l1muon_chargeValid;   //!
   TBranch        *b_l1muon_muonIndex;   //!
   TBranch        *b_l1muon_tag;   //!
   TBranch        *b_l1muon_isoSum;   //!
   TBranch        *b_l1muon_dPhiExtra;   //!
   TBranch        *b_l1muon_dEtaExtra;   //!
   TBranch        *b_l1muon_rank;   //!
   TBranch        *b_l1muon_bx;   //!
   TBranch        *b_l1egamma_count;   //!
   TBranch        *b_l1egamma_px;   //!
   TBranch        *b_l1egamma_py;   //!
   TBranch        *b_l1egamma_pz;   //!
   TBranch        *b_l1egamma_pt;   //!
   TBranch        *b_l1egamma_ipt;   //!
   TBranch        *b_l1egamma_eta;   //!
   TBranch        *b_l1egamma_phi;   //!
   TBranch        *b_l1egamma_qual;   //!
   TBranch        *b_l1egamma_iso;   //!
   TBranch        *b_l1egamma_towerIEta;   //!
   TBranch        *b_l1egamma_towerIPhi;   //!
   TBranch        *b_l1egamma_rawEt;   //!
   TBranch        *b_l1egamma_isoEt;   //!
   TBranch        *b_l1egamma_footprintEt;   //!
   TBranch        *b_l1egamma_nTT;   //!
   TBranch        *b_l1egamma_shape;   //!
   TBranch        *b_l1egamma_bx;   //!
   TBranch        *b_l1tau_count;   //!
   TBranch        *b_l1tau_px;   //!
   TBranch        *b_l1tau_py;   //!
   TBranch        *b_l1tau_pz;   //!
   TBranch        *b_l1tau_pt;   //!
   TBranch        *b_l1tau_ipt;   //!
   TBranch        *b_l1tau_eta;   //!
   TBranch        *b_l1tau_phi;   //!
   TBranch        *b_l1tau_qual;   //!
   TBranch        *b_l1tau_iso;   //!
   TBranch        *b_l1tau_towerIEta;   //!
   TBranch        *b_l1tau_towerIPhi;   //!
   TBranch        *b_l1tau_rawEt;   //!
   TBranch        *b_l1tau_isoEt;   //!
   TBranch        *b_l1tau_nTT;   //!
   TBranch        *b_l1tau_hasEM;   //!
   TBranch        *b_l1tau_isMerged;   //!
   TBranch        *b_l1tau_bx;   //!
   TBranch        *b_l1isotau_count;   //!
   TBranch        *b_l1isotau_e;   //!
   TBranch        *b_l1isotau_px;   //!
   TBranch        *b_l1isotau_py;   //!
   TBranch        *b_l1isotau_pz;   //!
   TBranch        *b_l1isotau_mass;   //!
   TBranch        *b_l1isotau_eta;   //!
   TBranch        *b_l1isotau_phi;   //!
   TBranch        *b_l1isotau_pt;   //!
   TBranch        *b_l1isotau_charge;   //!
   TBranch        *b_l1isotau_iso;   //!
   TBranch        *b_trigobject_count;   //!
   TBranch        *b_trigobject_px;   //!
   TBranch        *b_trigobject_py;   //!
   TBranch        *b_trigobject_pz;   //!
   TBranch        *b_trigobject_pt;   //!
   TBranch        *b_trigobject_eta;   //!
   TBranch        *b_trigobject_phi;   //!
   TBranch        *b_trigobject_filters;   //!
   TBranch        *b_trigobject_isMuon;   //!
   TBranch        *b_trigobject_isElectron;   //!
   TBranch        *b_trigobject_isTau;   //!
   TBranch        *b_trigobject_isJet;   //!
   TBranch        *b_trigobject_isMET;   //!
   TBranch        *b_run_hltnames;   //!
   TBranch        *b_run_hltfilters;   //!
   TBranch        *b_run_hltmufilters;   //!
   TBranch        *b_run_hltelectronfilters;   //!
   TBranch        *b_run_hlttaufilters;   //!
   TBranch        *b_run_hltphotonfilters;   //!
   TBranch        *b_run_hltjetfilters;   //!
   TBranch        *b_run_floattaudiscriminators;   //!
   TBranch        *b_run_binarytaudiscriminators;   //!
   TBranch        *b_run_btagdiscriminators;   //!
   TBranch        *b_hltriggerresults;   //!
   TBranch        *b_hltriggerprescales;   //!
   TBranch        *b_hltriggerresultsV;   //!
   TBranch        *b_flags;   //!
   TBranch        *b_tau_againstElectronLooseMVA6;   //!
   TBranch        *b_tau_againstElectronMVA6Raw;   //!
   TBranch        *b_tau_againstElectronMVA6category;   //!
   TBranch        *b_tau_againstElectronMediumMVA6;   //!
   TBranch        *b_tau_againstElectronTightMVA6;   //!
   TBranch        *b_tau_againstElectronVLooseMVA6;   //!
   TBranch        *b_tau_againstElectronVTightMVA6;   //!
   TBranch        *b_tau_againstMuonLoose3;   //!
   TBranch        *b_tau_againstMuonTight3;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWdR03oldDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWnewDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWoldDMwLTraw;   //!
   TBranch        *b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byPhotonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_tau_chargedIsoPtSumdR03;   //!
   TBranch        *b_tau_decayModeFinding;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_footprintCorrection;   //!
   TBranch        *b_tau_footprintCorrectiondR03;   //!
   TBranch        *b_tau_neutralIsoPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSumWeight;   //!
   TBranch        *b_tau_neutralIsoPtSumWeightdR03;   //!
   TBranch        *b_tau_neutralIsoPtSumdR03;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalConedR03;   //!
   TBranch        *b_tau_puCorrPtSum;   //!
   TBranch        *b_tau_byIsolationMVArun2017v1DBoldDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBnewDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBnewDMwLTraw2016;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBoldDMwLTraw2016;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch	  *b_tau_byDeepTau2017v2p1VSeraw;	//!
   TBranch	  *b_tau_byDeepTau2017v2p1VSjetraw;	//!
   TBranch	  *b_tau_byDeepTau2017v2p1VSmuraw;	//!
   TBranch	  *b_tau_byLooseDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byLooseDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byLooseDeepTau2017v2p1VSmu;	//!
   TBranch	  *b_tau_byMediumDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byMediumDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byMediumDeepTau2017v2p1VSmu;	//!
   TBranch	  *b_tau_byTightDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byTightDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byTightDeepTau2017v2p1VSmu;	//!
   TBranch	  *b_tau_byVLooseDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byVLooseDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byVLooseDeepTau2017v2p1VSmu;	//!
   TBranch	  *b_tau_byVTightDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byVTightDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byVVLooseDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byVVLooseDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byVVTightDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byVVTightDeepTau2017v2p1VSjet;	//!
   TBranch	  *b_tau_byVVVLooseDeepTau2017v2p1VSe;	//!
   TBranch	  *b_tau_byVVVLooseDeepTau2017v2p1VSjet;	//!
   
   TBranch        *b_tau_MVADM2017v1;
   TBranch        *b_tau_MVADM2017v1DM0raw;
   TBranch        *b_tau_MVADM2017v1DM10raw;
   TBranch        *b_tau_MVADM2017v1DM11raw;
   TBranch        *b_tau_MVADM2017v1DM1raw;
   TBranch        *b_tau_MVADM2017v1DM2raw;
   TBranch        *b_tau_MVADM2017v1DMotherraw;

   TBranch        *b_TauSpinAngles_count;
   TBranch        *b_TauSpinnerWeight;

   TBranch        *b_htxs_stage0cat;   //!
   TBranch        *b_htxs_stage1cat;   //!
   TBranch        *b_htxs_stage1p1cat;
   TBranch        *b_htxs_higgsPt;   //!
   TBranch        *b_htxs_njets30;   //!
    
   AC1B(TTree *tree=0, bool isData=false);
   virtual ~AC1B();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool isData);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Long64_t GetEntries();
};

#endif

#ifdef AC1B_cxx
AC1B::AC1B(TTree *tree, bool isData) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_MC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_MC.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("output_MC.root:/makeroottree");
      dir->GetObject("AC1B",tree);

   }
   Init(tree, isData);
}

AC1B::~AC1B()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AC1B::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AC1B::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AC1B::Init(TTree *tree, bool isData)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   tree->SetMaxVirtualSize(3000000);


   // Set object pointer
   run_hltnames = 0;
   run_hltfilters = 0;
   run_hltmufilters = 0;
   run_hltelectronfilters = 0;
   run_hlttaufilters = 0;
   run_hltphotonfilters = 0;
   run_hltjetfilters = 0;
   run_floattaudiscriminators = 0;
   run_binarytaudiscriminators = 0;
   run_btagdiscriminators = 0;
   hltriggerresults = 0;
   hltriggerprescales = 0;
   hltriggerresultsV = 0;
   flags = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("errors", &errors, &b_errors);
   fChain->SetBranchAddress("event_nr", &event_nr, &b_event_nr);
   fChain->SetBranchAddress("event_run", &event_run, &b_event_run);
   fChain->SetBranchAddress("event_timeunix", &event_timeunix, &b_event_timeunix);
   fChain->SetBranchAddress("event_timemicrosec", &event_timemicrosec, &b_event_timemicrosec);
   fChain->SetBranchAddress("event_luminosityblock", &event_luminosityblock, &b_event_luminosityblock);
   fChain->SetBranchAddress("trigger_level1bits", trigger_level1bits, &b_trigger_level1bits);
   fChain->SetBranchAddress("trigger_level1", trigger_level1, &b_trigger_level1);
   fChain->SetBranchAddress("trigger_HLT", trigger_HLT, &b_trigger_HLT);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("primvertex_count", &primvertex_count, &b_primvertex_count);
   fChain->SetBranchAddress("goodprimvertex_count", &goodprimvertex_count, &b_goodprimvertex_count);
   fChain->SetBranchAddress("primvertex_x", &primvertex_x, &b_primvertex_x);
   fChain->SetBranchAddress("primvertex_y", &primvertex_y, &b_primvertex_y);
   fChain->SetBranchAddress("primvertex_z", &primvertex_z, &b_primvertex_z);
   fChain->SetBranchAddress("primvertex_chi2", &primvertex_chi2, &b_primvertex_chi2);
   fChain->SetBranchAddress("primvertex_ndof", &primvertex_ndof, &b_primvertex_ndof);
   fChain->SetBranchAddress("primvertex_ptq", &primvertex_ptq, &b_primvertex_pdf);
   fChain->SetBranchAddress("primvertex_ntracks", &primvertex_ntracks, &b_primvertex_ntracks);
   fChain->SetBranchAddress("primvertex_cov", primvertex_cov, &b_primvertex_cov);
   
   fChain->SetBranchAddress("primvertexwithbs_count", &primvertexwithbs_count, &b_primvertexwithbs_count);
   fChain->SetBranchAddress("goodprimvertexwithbs_count", &goodprimvertexwithbs_count, &b_goodprimvertexwithbs_count);
   fChain->SetBranchAddress("primvertexwithbs_x", &primvertexwithbs_x, &b_primvertexwithbs_x);
   fChain->SetBranchAddress("primvertexwithbs_y", &primvertexwithbs_y, &b_primvertexwithbs_y);
   fChain->SetBranchAddress("primvertexwithbs_z", &primvertexwithbs_z, &b_primvertexwithbs_z);
   fChain->SetBranchAddress("primvertexwithbs_chi2", &primvertexwithbs_chi2, &b_primvertexwithbs_chi2);
   fChain->SetBranchAddress("primvertexwithbs_ndof", &primvertexwithbs_ndof, &b_primvertexwithbs_ndof);
   fChain->SetBranchAddress("primvertexwithbs_ptq", &primvertexwithbs_ptq, &b_primvertexwithbs_pdf);
   fChain->SetBranchAddress("primvertexwithbs_ntracks", &primvertexwithbs_ntracks, &b_primvertexwithbs_ntracks);
   fChain->SetBranchAddress("primvertexwithbs_cov", primvertexwithbs_cov, &b_primvertexwithbs_cov);
   
   //refit vertices
   fChain->SetBranchAddress("refitvertex_count", &refitvertex_count, &b_refitvertex_count);
   //fChain->SetBranchAddress("goodrefitvertex_count", &goodrefitvertex_count, &b_goodrefitvertex_count);
   fChain->SetBranchAddress("refitvertex_x", refitvertex_x, &b_refitvertex_x);
   fChain->SetBranchAddress("refitvertex_y", refitvertex_y, &b_refitvertex_y);
   fChain->SetBranchAddress("refitvertex_z", refitvertex_z, &b_refitvertex_z);
   fChain->SetBranchAddress("refitvertex_chi2", refitvertex_chi2, &b_refitvertex_chi2);
   fChain->SetBranchAddress("refitvertex_ndof", refitvertex_ndof, &b_refitvertex_ndof);
   fChain->SetBranchAddress("refitvertex_ptq", refitvertex_ptq, &b_refitvertex_pdf);
   fChain->SetBranchAddress("refitvertex_ntracks",refitvertex_ntracks, &b_refitvertex_ntracks);
   //fChain->SetBranchAddress("refitvertex_cov",refitvertex_cov, &b_refitvertex_cov);
   fChain->SetBranchAddress("refitvertex_eleIndex", refitvertex_eleIndex, &b_refitvertex_eleIndex);
   fChain->SetBranchAddress("refitvertex_muIndex", refitvertex_muIndex, &b_refitvertex_muIndex);
   fChain->SetBranchAddress("refitvertex_tauIndex", refitvertex_tauIndex, &b_refitvertex_tauIndex);
   
   fChain->SetBranchAddress("refitvertexwithbs_count", &refitvertexwithbs_count, &b_refitvertexwithbs_count);
   //fChain->SetBranchAddress("goodrefitvertex_count", &goodrefitvertex_count, &b_goodrefitvertex_count);
   fChain->SetBranchAddress("refitvertexwithbs_x", refitvertexwithbs_x, &b_refitvertexwithbs_x);
   fChain->SetBranchAddress("refitvertexwithbs_y", refitvertexwithbs_y, &b_refitvertexwithbs_y);
   fChain->SetBranchAddress("refitvertexwithbs_z", refitvertexwithbs_z, &b_refitvertexwithbs_z);
   fChain->SetBranchAddress("refitvertexwithbs_chi2", refitvertexwithbs_chi2, &b_refitvertexwithbs_chi2);
   fChain->SetBranchAddress("refitvertexwithbs_ndof", refitvertexwithbs_ndof, &b_refitvertexwithbs_ndof);
   fChain->SetBranchAddress("refitvertexwithbs_ptq", refitvertexwithbs_ptq, &b_refitvertexwithbs_pdf);
   fChain->SetBranchAddress("refitvertexwithbs_ntracks",refitvertexwithbs_ntracks, &b_refitvertexwithbs_ntracks);
   //fChain->SetBranchAddress("refitvertexwithbs_cov",refitvertexwithbs_cov, &b_refitvertexwithbs_cov);
   fChain->SetBranchAddress("refitvertexwithbs_eleIndex", refitvertexwithbs_eleIndex, &b_refitvertexwithbs_eleIndex);
   fChain->SetBranchAddress("refitvertexwithbs_muIndex", refitvertexwithbs_muIndex, &b_refitvertexwithbs_muIndex);
   fChain->SetBranchAddress("refitvertexwithbs_tauIndex", refitvertexwithbs_tauIndex, &b_refitvertexwithbs_tauIndex);
   
   fChain->SetBranchAddress("muon_count", &muon_count, &b_muon_count);
   fChain->SetBranchAddress("muon_px", muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon_py", muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon_pz", muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_pterror", muon_pterror, &b_muon_pterror);
   fChain->SetBranchAddress("muon_chi2", muon_chi2, &b_muon_chi2);
   fChain->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon_ndof", muon_ndof, &b_muon_ndof);
   fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_miniISO", muon_miniISO, &b_muon_miniISO);
   fChain->SetBranchAddress("muon_combQ_chi2LocalPosition", muon_combQ_chi2LocalPosition, &b_muon_combQ_chi2LocalPosition);
   fChain->SetBranchAddress("muon_combQ_trkKink", muon_combQ_trkKink, &b_muon_combQ_trkKink);
   fChain->SetBranchAddress("muon_validFraction", muon_validFraction, &b_muon_validFraction);
   fChain->SetBranchAddress("muon_segmentComp", muon_segmentComp, &b_muon_segmentComp);
   fChain->SetBranchAddress("muon_nMuonStations", muon_nMuonStations, &b_muon_nMuonStations);
   fChain->SetBranchAddress("muon_nMuonHits", muon_nMuonHits, &b_muon_nMuonHits);
   fChain->SetBranchAddress("muon_nPixelHits", muon_nPixelHits, &b_muon_nPixelHits);
   fChain->SetBranchAddress("muon_nTrackerHits", muon_nTrackerHits, &b_muon_nTrackerHits);
   fChain->SetBranchAddress("muon_dxy", muon_dxy, &b_muon_dxy);
   fChain->SetBranchAddress("muon_dxyerr", muon_dxyerr, &b_muon_dxyerr);
   fChain->SetBranchAddress("muon_dz", muon_dz, &b_muon_dz);
   fChain->SetBranchAddress("muon_dzerr", muon_dzerr, &b_muon_dzerr);
   fChain->SetBranchAddress("muon_vx", muon_vx, &b_muon_vx);
   fChain->SetBranchAddress("muon_vy", muon_vy, &b_muon_vy);
   fChain->SetBranchAddress("muon_vz", muon_vz, &b_muon_vz);
   fChain->SetBranchAddress("muon_chargedHadIso", muon_chargedHadIso, &b_muon_chargedHadIso);
   fChain->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso, &b_muon_neutralHadIso);
   fChain->SetBranchAddress("muon_photonIso", muon_photonIso, &b_muon_photonIso);
   fChain->SetBranchAddress("muon_puIso", muon_puIso, &b_muon_puIso);
   fChain->SetBranchAddress("muon_r03_sumChargedHadronPt", muon_r03_sumChargedHadronPt, &b_muon_r03_sumChargedHadronPt);
   fChain->SetBranchAddress("muon_r03_sumChargedParticlePt", muon_r03_sumChargedParticlePt, &b_muon_r03_sumChargedParticlePt);
   fChain->SetBranchAddress("muon_r03_sumNeutralHadronEt", muon_r03_sumNeutralHadronEt, &b_muon_r03_sumNeutralHadronEt);
   fChain->SetBranchAddress("muon_r03_sumPhotonEt", muon_r03_sumPhotonEt, &b_muon_r03_sumPhotonEt);
   fChain->SetBranchAddress("muon_r03_sumNeutralHadronEtHighThreshold", muon_r03_sumNeutralHadronEtHighThreshold, &b_muon_r03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("muon_r03_sumPhotonEtHighThreshold", muon_r03_sumPhotonEtHighThreshold, &b_muon_r03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("muon_r03_sumPUPt", muon_r03_sumPUPt, &b_muon_r03_sumPUPt);
   fChain->SetBranchAddress("muon_r04_sumChargedHadronPt", muon_r04_sumChargedHadronPt, &b_muon_r04_sumChargedHadronPt);
   fChain->SetBranchAddress("muon_r04_sumChargedParticlePt", muon_r04_sumChargedParticlePt, &b_muon_r04_sumChargedParticlePt);
   fChain->SetBranchAddress("muon_r04_sumNeutralHadronEt", muon_r04_sumNeutralHadronEt, &b_muon_r04_sumNeutralHadronEt);
   fChain->SetBranchAddress("muon_r04_sumPhotonEt", muon_r04_sumPhotonEt, &b_muon_r04_sumPhotonEt);
   fChain->SetBranchAddress("muon_r04_sumNeutralHadronEtHighThreshold", muon_r04_sumNeutralHadronEtHighThreshold, &b_muon_r04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("muon_r04_sumPhotonEtHighThreshold", muon_r04_sumPhotonEtHighThreshold, &b_muon_r04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("muon_r04_sumPUPt", muon_r04_sumPUPt, &b_muon_r04_sumPUPt);
   fChain->SetBranchAddress("muon_isPF", muon_isPF, &b_muon_isPF);
   fChain->SetBranchAddress("muon_isGlobal", muon_isGlobal, &b_muon_isGlobal);
   fChain->SetBranchAddress("muon_isTracker", muon_isTracker, &b_muon_isTracker);
   fChain->SetBranchAddress("muon_isTight", muon_isTight, &b_muon_isTight);
   fChain->SetBranchAddress("muon_isLoose", muon_isLoose, &b_muon_isLoose);
   fChain->SetBranchAddress("muon_isMedium", muon_isMedium, &b_muon_isMedium);
   fChain->SetBranchAddress("muon_isICHEP", muon_isICHEP, &b_muon_isICHEP);
   fChain->SetBranchAddress("muon_genmatch", muon_genmatch, &b_muon_genmatch);
   fChain->SetBranchAddress("muon_isDuplicate", muon_isDuplicate, &b_muon_isDuplicate);
   fChain->SetBranchAddress("muon_isBad", muon_isBad, &b_muon_isBad);
   fChain->SetBranchAddress("muon_helixparameters", muon_helixparameters, &b_muon_helixparameters);
   fChain->SetBranchAddress("muon_helixparameters_covar", muon_helixparameters_covar, &b_muon_helixparameters_covar);
   fChain->SetBranchAddress("muon_referencePoint", muon_referencePoint, &b_muon_referencePoint);
   fChain->SetBranchAddress("muon_Bfield", muon_Bfield, &b_muon_Bfield);
   fChain->SetBranchAddress("dimuon_count", &dimuon_count, &b_dimuon_count);
   fChain->SetBranchAddress("dimuon_leading", dimuon_leading, &b_dimuon_leading);
   fChain->SetBranchAddress("dimuon_trailing", dimuon_trailing, &b_dimuon_trailing);
   fChain->SetBranchAddress("dimuon_dist2D", dimuon_dist2D, &b_dimuon_dist2D);
   fChain->SetBranchAddress("dimuon_dist2DE", dimuon_dist2DE, &b_dimuon_dist2DE);
   fChain->SetBranchAddress("dimuon_dist3D", dimuon_dist3D, &b_dimuon_dist3D);
   fChain->SetBranchAddress("dimuon_dist3DE", dimuon_dist3DE, &b_dimuon_dist3DE);
   fChain->SetBranchAddress("pfjet_count", &pfjet_count, &b_pfjet_count);
   fChain->SetBranchAddress("pfjet_e", pfjet_e, &b_pfjet_e);
   fChain->SetBranchAddress("pfjet_px", pfjet_px, &b_pfjet_px);
   fChain->SetBranchAddress("pfjet_py", pfjet_py, &b_pfjet_py);
   fChain->SetBranchAddress("pfjet_pz", pfjet_pz, &b_pfjet_pz);
   fChain->SetBranchAddress("pfjet_pt", pfjet_pt, &b_pfjet_pt);
   fChain->SetBranchAddress("pfjet_eta", pfjet_eta, &b_pfjet_eta);
   fChain->SetBranchAddress("pfjet_phi", pfjet_phi, &b_pfjet_phi);
   fChain->SetBranchAddress("pfjet_neutralhadronicenergy", pfjet_neutralhadronicenergy, &b_pfjet_neutralhadronicenergy);
   fChain->SetBranchAddress("pfjet_chargedhadronicenergy", pfjet_chargedhadronicenergy, &b_pfjet_chargedhadronicenergy);
   fChain->SetBranchAddress("pfjet_neutralemenergy", pfjet_neutralemenergy, &b_pfjet_neutralemenergy);
   fChain->SetBranchAddress("pfjet_chargedemenergy", pfjet_chargedemenergy, &b_pfjet_chargedemenergy);
   fChain->SetBranchAddress("pfjet_muonenergy", pfjet_muonenergy, &b_pfjet_muonenergy);
   fChain->SetBranchAddress("pfjet_chargedmuonenergy", pfjet_chargedmuonenergy, &b_pfjet_chargedmuonenergy);
   fChain->SetBranchAddress("pfjet_chargedmulti", pfjet_chargedmulti, &b_pfjet_chargedmulti);
   fChain->SetBranchAddress("pfjet_neutralmulti", pfjet_neutralmulti, &b_pfjet_neutralmulti);
   fChain->SetBranchAddress("pfjet_chargedhadronmulti", pfjet_chargedhadronmulti, &b_pfjet_chargedhadronmulti);
   fChain->SetBranchAddress("pfjet_energycorr", pfjet_energycorr, &b_pfjet_energycorr);
   fChain->SetBranchAddress("pfjet_energycorr_l1fastjet", pfjet_energycorr_l1fastjet, &b_pfjet_energycorr_l1fastjet);
   fChain->SetBranchAddress("pfjet_energycorr_l2relative", pfjet_energycorr_l2relative, &b_pfjet_energycorr_l2relative);
   fChain->SetBranchAddress("pfjet_energycorr_l3absolute", pfjet_energycorr_l3absolute, &b_pfjet_energycorr_l3absolute);
   fChain->SetBranchAddress("pfjet_energycorr_l2l3residual", pfjet_energycorr_l2l3residual, &b_pfjet_energycorr_l2l3residual);
   fChain->SetBranchAddress("pfjet_flavour", pfjet_flavour, &b_pfjet_flavour);
   fChain->SetBranchAddress("pfjet_btag", pfjet_btag, &b_pfjet_btag);
   fChain->SetBranchAddress("pfjet_jecUncertainty", pfjet_jecUncertainty, &b_pfjet_jecUncertainty);
   fChain->SetBranchAddress("pfjet_pu_jet_fullId_loose", pfjet_pu_jet_fullId_loose, &b_pfjet_pu_jet_fullId_loose);
   fChain->SetBranchAddress("pfjet_pu_jet_fullId_medium", pfjet_pu_jet_fullId_medium, &b_pfjet_pu_jet_fullId_medium);
   fChain->SetBranchAddress("pfjet_pu_jet_fullId_tight", pfjet_pu_jet_fullId_tight, &b_pfjet_pu_jet_fullId_tight);
   fChain->SetBranchAddress("pfjet_pu_jet_fullDisc_mva", pfjet_pu_jet_fullDisc_mva, &b_pfjet_pu_jet_fullDisc_mva);
   fChain->SetBranchAddress("electron_count", &electron_count, &b_electron_count);
   fChain->SetBranchAddress("electron_px", electron_px, &b_electron_px);
   fChain->SetBranchAddress("electron_py", electron_py, &b_electron_py);
   fChain->SetBranchAddress("electron_pz", electron_pz, &b_electron_pz);
   fChain->SetBranchAddress("electron_pt", electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_eta", electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_phi", electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_trackchi2", electron_trackchi2, &b_electron_trackchi2);
   fChain->SetBranchAddress("electron_trackndof", electron_trackndof, &b_electron_trackndof);
   fChain->SetBranchAddress("electron_outerx", electron_outerx, &b_electron_outerx);
   fChain->SetBranchAddress("electron_outery", electron_outery, &b_electron_outery);
   fChain->SetBranchAddress("electron_outerz", electron_outerz, &b_electron_outerz);
   fChain->SetBranchAddress("electron_vx", electron_vx, &b_electron_vx);
   fChain->SetBranchAddress("electron_vy", electron_vy, &b_electron_vy);
   fChain->SetBranchAddress("electron_vz", electron_vz, &b_electron_vz);
   fChain->SetBranchAddress("electron_esuperclusterovertrack", electron_esuperclusterovertrack, &b_electron_esuperclusterovertrack);
   fChain->SetBranchAddress("electron_eseedclusterovertrack", electron_eseedclusterovertrack, &b_electron_eseedclusterovertrack);
   fChain->SetBranchAddress("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack, &b_electron_deltaetasuperclustertrack);
   fChain->SetBranchAddress("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack, &b_electron_deltaphisuperclustertrack);
   fChain->SetBranchAddress("electron_e1x5", electron_e1x5, &b_electron_e1x5);
   fChain->SetBranchAddress("electron_e2x5", electron_e2x5, &b_electron_e2x5);
   fChain->SetBranchAddress("electron_e5x5", electron_e5x5, &b_electron_e5x5);
   fChain->SetBranchAddress("electron_sigmaetaeta", electron_sigmaetaeta, &b_electron_sigmaetaeta);
   fChain->SetBranchAddress("electron_sigmaietaieta", electron_sigmaietaieta, &b_electron_sigmaietaieta);
   fChain->SetBranchAddress("electron_ehcaloverecal", electron_ehcaloverecal, &b_electron_ehcaloverecal);
   fChain->SetBranchAddress("electron_ehcaloverecaldepth1", electron_ehcaloverecaldepth1, &b_electron_ehcaloverecaldepth1);
   fChain->SetBranchAddress("electron_ehcaloverecaldepth2", electron_ehcaloverecaldepth2, &b_electron_ehcaloverecaldepth2);
   fChain->SetBranchAddress("electron_full5x5_sigmaietaieta", electron_full5x5_sigmaietaieta, &b_electron_full5x5_sigmaietaieta);
   fChain->SetBranchAddress("electron_ooemoop", electron_ooemoop, &b_electron_ooemoop);
   fChain->SetBranchAddress("electron_miniISO", electron_miniISO, &b_electron_miniISO);
   fChain->SetBranchAddress("electron_superclusterEta", electron_superclusterEta, &b_electron_superclusterEta);
   fChain->SetBranchAddress("electron_superclusterPhi", electron_superclusterPhi, &b_electron_superclusterPhi);
   fChain->SetBranchAddress("electron_superclusterX", electron_superclusterX, &b_electron_superclusterX);
   fChain->SetBranchAddress("electron_superclusterY", electron_superclusterY, &b_electron_superclusterY);
   fChain->SetBranchAddress("electron_superclusterZ", electron_superclusterZ, &b_electron_superclusterZ);
   fChain->SetBranchAddress("electron_detaInSeed", electron_detaInSeed, &b_electron_detaInSeed);
   fChain->SetBranchAddress("electron_he", electron_he, &b_electron_he);
   fChain->SetBranchAddress("electron_eaIsolation", electron_eaIsolation, &b_electron_eaIsolation);
   fChain->SetBranchAddress("electron_chargedHadIso", electron_chargedHadIso, &b_electron_chargedHadIso);
   fChain->SetBranchAddress("electron_neutralHadIso", electron_neutralHadIso, &b_electron_neutralHadIso);
   fChain->SetBranchAddress("electron_photonIso", electron_photonIso, &b_electron_photonIso);
   fChain->SetBranchAddress("electron_puIso", electron_puIso, &b_electron_puIso);
   fChain->SetBranchAddress("electron_r03_sumChargedHadronPt", electron_r03_sumChargedHadronPt, &b_electron_r03_sumChargedHadronPt);
   fChain->SetBranchAddress("electron_r03_sumChargedParticlePt", electron_r03_sumChargedParticlePt, &b_electron_r03_sumChargedParticlePt);
   fChain->SetBranchAddress("electron_r03_sumNeutralHadronEt", electron_r03_sumNeutralHadronEt, &b_electron_r03_sumNeutralHadronEt);
   fChain->SetBranchAddress("electron_r03_sumPhotonEt", electron_r03_sumPhotonEt, &b_electron_r03_sumPhotonEt);
   fChain->SetBranchAddress("electron_r03_sumNeutralHadronEtHighThreshold", electron_r03_sumNeutralHadronEtHighThreshold, &b_electron_r03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("electron_r03_sumPhotonEtHighThreshold", electron_r03_sumPhotonEtHighThreshold, &b_electron_r03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("electron_r03_sumPUPt", electron_r03_sumPUPt, &b_electron_r03_sumPUPt);
   fChain->SetBranchAddress("electron_nhits", electron_nhits, &b_electron_nhits);
   fChain->SetBranchAddress("electron_npixelhits", electron_npixelhits, &b_electron_npixelhits);
   fChain->SetBranchAddress("electron_nmissinghits", electron_nmissinghits, &b_electron_nmissinghits);
   fChain->SetBranchAddress("electron_nmissinginnerhits", electron_nmissinginnerhits, &b_electron_nmissinginnerhits);
   fChain->SetBranchAddress("electron_npixellayers", electron_npixellayers, &b_electron_npixellayers);
   fChain->SetBranchAddress("electron_nstriplayers", electron_nstriplayers, &b_electron_nstriplayers);
   fChain->SetBranchAddress("electron_dxy", electron_dxy, &b_electron_dxy);
   fChain->SetBranchAddress("electron_dxyerr", electron_dxyerr, &b_electron_dxyerr);
   fChain->SetBranchAddress("electron_dz", electron_dz, &b_electron_dz);
   fChain->SetBranchAddress("electron_dzerr", electron_dzerr, &b_electron_dzerr);
   fChain->SetBranchAddress("electron_convdist", electron_convdist, &b_electron_convdist);
   fChain->SetBranchAddress("electron_gapinfo", electron_gapinfo, &b_electron_gapinfo);
   fChain->SetBranchAddress("electron_chargeinfo", electron_chargeinfo, &b_electron_chargeinfo);
   fChain->SetBranchAddress("electron_fbrems", electron_fbrems, &b_electron_fbrems);
   fChain->SetBranchAddress("electron_numbrems", electron_numbrems, &b_electron_numbrems);
   fChain->SetBranchAddress("electron_charge", electron_charge, &b_electron_charge);
   fChain->SetBranchAddress("electron_superclusterindex", electron_superclusterindex, &b_electron_superclusterindex);
   fChain->SetBranchAddress("electron_info", electron_info, &b_electron_info);
   fChain->SetBranchAddress("electron_mva_value_nontrig_Spring15_v1", electron_mva_value_nontrig_Spring15_v1, &b_electron_mva_value_nontrig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_value_trig_Spring15_v1", electron_mva_value_trig_Spring15_v1, &b_electron_mva_value_trig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_category_nontrig_Spring15_v1", electron_mva_category_nontrig_Spring15_v1, &b_electron_mva_category_nontrig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_category_trig_Spring15_v1", electron_mva_category_trig_Spring15_v1, &b_electron_mva_category_trig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_wp80_nontrig_Spring15_v1", electron_mva_wp80_nontrig_Spring15_v1, &b_electron_mva_wp80_nontrig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_wp90_nontrig_Spring15_v1", electron_mva_wp90_nontrig_Spring15_v1, &b_electron_mva_wp90_nontrig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_wp80_trig_Spring15_v1", electron_mva_wp80_trig_Spring15_v1, &b_electron_mva_wp80_trig_Spring15_v1);
   fChain->SetBranchAddress("electron_mva_wp90_trig_Spring15_v1", electron_mva_wp90_trig_Spring15_v1, &b_electron_mva_wp90_trig_Spring15_v1);
   fChain->SetBranchAddress("electron_cutId_veto_Spring15", electron_cutId_veto_Spring15, &b_electron_cutId_veto_Spring15);
   fChain->SetBranchAddress("electron_cutId_loose_Spring15", electron_cutId_loose_Spring15, &b_electron_cutId_loose_Spring15);
   fChain->SetBranchAddress("electron_cutId_medium_Spring15", electron_cutId_medium_Spring15, &b_electron_cutId_medium_Spring15);
   fChain->SetBranchAddress("electron_cutId_tight_Spring15", electron_cutId_tight_Spring15, &b_electron_cutId_tight_Spring15);
   fChain->SetBranchAddress("electron_cutId_veto_Summer16", electron_cutId_veto_Summer16, &b_electron_cutId_veto_Summer16);
   fChain->SetBranchAddress("electron_cutId_loose_Summer16", electron_cutId_loose_Summer16, &b_electron_cutId_loose_Summer16);
   fChain->SetBranchAddress("electron_cutId_medium_Summer16", electron_cutId_medium_Summer16, &b_electron_cutId_medium_Summer16);
   fChain->SetBranchAddress("electron_cutId_tight_Summer16", electron_cutId_tight_Summer16, &b_electron_cutId_tight_Summer16);
   fChain->SetBranchAddress("electron_mva_value_Spring16_v1", electron_mva_value_Spring16_v1, &b_electron_mva_value_Spring16_v1);
   fChain->SetBranchAddress("electron_mva_category_Spring16_v1", electron_mva_category_Spring16_v1, &b_electron_mva_category_Spring16_v1);
   fChain->SetBranchAddress("electron_mva_wp90_general_Spring16_v1", electron_mva_wp90_general_Spring16_v1, &b_electron_mva_wp90_general_Spring16_v1);
   fChain->SetBranchAddress("electron_mva_wp80_general_Spring16_v1", electron_mva_wp80_general_Spring16_v1, &b_electron_mva_wp80_general_Spring16_v1);
   fChain->SetBranchAddress("electron_mva_value_Iso_Fall17_v1", electron_mva_value_Iso_Fall17_v1, &b_electron_mva_value_Iso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_value_noIso_Fall17_v1", electron_mva_value_noIso_Fall17_v1, &b_electron_mva_value_noIso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_wp90_Iso_Fall17_v1", electron_mva_wp90_Iso_Fall17_v1, &b_electron_mva_wp90_Iso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_wp80_Iso_Fall17_v1", electron_mva_wp80_Iso_Fall17_v1, &b_electron_mva_wp80_Iso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_Loose_Iso_Fall17_v1", electron_mva_Loose_Iso_Fall17_v1, &b_electron_mva_Loose_Iso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_wp90_noIso_Fall17_v1", electron_mva_wp90_noIso_Fall17_v1, &b_electron_mva_wp90_noIso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_wp90_noIso_Fall17_v2", electron_mva_wp90_noIso_Fall17_v2, &b_electron_mva_wp90_noIso_Fall17_v2);
   fChain->SetBranchAddress("electron_mva_wp80_noIso_Fall17_v1", electron_mva_wp80_noIso_Fall17_v1, &b_electron_mva_wp80_noIso_Fall17_v1);
   fChain->SetBranchAddress("electron_mva_Loose_noIso_Fall17_v1", electron_mva_Loose_noIso_Fall17_v1, &b_electron_mva_Loose_noIso_Fall17_v1);
   fChain->SetBranchAddress("electron_cutId_veto_Fall17", electron_cutId_veto_Fall17, &b_electron_cutId_veto_Fall17);
   fChain->SetBranchAddress("electron_cutId_loose_Fall17", electron_cutId_loose_Fall17, &b_electron_cutId_loose_Fall17);
   fChain->SetBranchAddress("electron_cutId_medium_Fall17", electron_cutId_medium_Fall17, &b_electron_cutId_medium_Fall17);
   fChain->SetBranchAddress("electron_cutId_tight_Fall17", electron_cutId_tight_Fall17, &b_electron_cutId_tight_Fall17);
   fChain->SetBranchAddress("electron_cutId_veto_Fall17V2", electron_cutId_veto_Fall17V2, &b_electron_cutId_veto_Fall17V2);
   fChain->SetBranchAddress("electron_cutId_loose_Fall17V2", electron_cutId_loose_Fall17V2, &b_electron_cutId_loose_Fall17V2);
   fChain->SetBranchAddress("electron_cutId_medium_Fall17V2", electron_cutId_medium_Fall17V2, &b_electron_cutId_medium_Fall17V2);
   fChain->SetBranchAddress("electron_cutId_tight_Fall17V2", electron_cutId_tight_Fall17V2, &b_electron_cutId_tight_Fall17V2);
   fChain->SetBranchAddress("electron_pass_conversion", electron_pass_conversion, &b_electron_pass_conversion);
   fChain->SetBranchAddress("electron_genmatch", electron_genmatch, &b_electron_genmatch);
   fChain->SetBranchAddress("electron_px_energyscale_up", electron_px_energyscale_up, &b_electron_px_energyscale_up);
   fChain->SetBranchAddress("electron_px_energyscale_down", electron_px_energyscale_down, &b_electron_px_energyscale_down);
   fChain->SetBranchAddress("electron_py_energyscale_up", electron_py_energyscale_up, &b_electron_py_energyscale_up);
   fChain->SetBranchAddress("electron_py_energyscale_down", electron_py_energyscale_down, &b_electron_py_energyscale_down);
   fChain->SetBranchAddress("electron_pz_energyscale_up", electron_pz_energyscale_up, &b_electron_pz_energyscale_up);
   fChain->SetBranchAddress("electron_pz_energyscale_down", electron_pz_energyscale_down, &b_electron_pz_energyscale_down);
   fChain->SetBranchAddress("electron_pt_energyscale_up", electron_pt_energyscale_up, &b_electron_pt_energyscale_up);
   fChain->SetBranchAddress("electron_pt_energyscale_down", electron_pt_energyscale_down, &b_electron_pt_energyscale_down);
   fChain->SetBranchAddress("electron_px_energysigma_up", electron_px_energysigma_up, &b_electron_px_energysigma_up);
   fChain->SetBranchAddress("electron_px_energysigma_down", electron_px_energysigma_down, &b_electron_px_energysigma_down);
   fChain->SetBranchAddress("electron_py_energysigma_up", electron_py_energysigma_up, &b_electron_py_energysigma_up);
   fChain->SetBranchAddress("electron_py_energysigma_down", electron_py_energysigma_down, &b_electron_py_energysigma_down);
   fChain->SetBranchAddress("electron_pz_energysigma_up", electron_pz_energysigma_up, &b_electron_pz_energysigma_up);
   fChain->SetBranchAddress("electron_pz_energysigma_down", electron_pz_energysigma_down, &b_electron_pz_energysigma_down);
   fChain->SetBranchAddress("electron_pt_energysigma_up", electron_pt_energysigma_up, &b_electron_pt_energysigma_up);
   fChain->SetBranchAddress("electron_pt_energysigma_down", electron_pt_energysigma_down, &b_electron_pt_energysigma_down);
   fChain->SetBranchAddress("tau_count", &tau_count, &b_tau_count);
   fChain->SetBranchAddress("tau_e", tau_e, &b_tau_e);
   fChain->SetBranchAddress("tau_px", tau_px, &b_tau_px);
   fChain->SetBranchAddress("tau_py", tau_py, &b_tau_py);
   fChain->SetBranchAddress("tau_pz", tau_pz, &b_tau_pz);
   fChain->SetBranchAddress("tau_mass", tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_eta", tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_vertexx", tau_vertexx, &b_tau_vertexx);
   fChain->SetBranchAddress("tau_vertexy", tau_vertexy, &b_tau_vertexy);
   fChain->SetBranchAddress("tau_vertexz", tau_vertexz, &b_tau_vertexz);
   fChain->SetBranchAddress("tau_pca2D_x", tau_pca2D_x, &b_tau_pca2D_x);
   fChain->SetBranchAddress("tau_pca2D_y", tau_pca2D_y, &b_tau_pca2D_y);
   fChain->SetBranchAddress("tau_pca2D_z", tau_pca2D_z, &b_tau_pca2D_z);
   fChain->SetBranchAddress("tau_pca3D_x", tau_pca3D_x, &b_tau_pca3D_x);
   fChain->SetBranchAddress("tau_pca3D_y", tau_pca3D_y, &b_tau_pca3D_y);
   fChain->SetBranchAddress("tau_pca3D_z", tau_pca3D_z, &b_tau_pca3D_z);
   fChain->SetBranchAddress("tau_SV_x", tau_SV_x, &b_tau_SV_x);
   fChain->SetBranchAddress("tau_SV_y", tau_SV_y, &b_tau_SV_y);
   fChain->SetBranchAddress("tau_SV_z", tau_SV_z, &b_tau_SV_z);
   fChain->SetBranchAddress("tau_SV_cov", tau_SV_cov, &b_tau_SV_cov);
   fChain->SetBranchAddress("tau_dxy", tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("tau_dz", tau_dz, &b_tau_dz);
   fChain->SetBranchAddress("tau_ip3d", tau_ip3d, &b_tau_ip3d);
   fChain->SetBranchAddress("tau_ip3dSig", tau_ip3dSig, &b_tau_ip3dSig);
   fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_genjet_px", tau_genjet_px, &b_tau_genjet_px);
   fChain->SetBranchAddress("tau_genjet_py", tau_genjet_py, &b_tau_genjet_py);
   fChain->SetBranchAddress("tau_genjet_pz", tau_genjet_pz, &b_tau_genjet_pz);
   fChain->SetBranchAddress("tau_genjet_e", tau_genjet_e, &b_tau_genjet_e);
   fChain->SetBranchAddress("tau_genmatch", tau_genmatch, &b_tau_genmatch);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_px", tau_leadchargedhadrcand_px, &b_tau_leadchargedhadrcand_px);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_py", tau_leadchargedhadrcand_py, &b_tau_leadchargedhadrcand_py);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_pz", tau_leadchargedhadrcand_pz, &b_tau_leadchargedhadrcand_pz);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_mass", tau_leadchargedhadrcand_mass, &b_tau_leadchargedhadrcand_mass);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_id", tau_leadchargedhadrcand_id, &b_tau_leadchargedhadrcand_id);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_dxy", tau_leadchargedhadrcand_dxy, &b_tau_leadchargedhadrcand_dxy);
   fChain->SetBranchAddress("tau_leadchargedhadrcand_dz", tau_leadchargedhadrcand_dz, &b_tau_leadchargedhadrcand_dz);
   fChain->SetBranchAddress("tau_ntracks_pt05", tau_ntracks_pt05, &b_tau_ntracks_pt05);
   fChain->SetBranchAddress("tau_ntracks_pt08", tau_ntracks_pt08, &b_tau_ntracks_pt08);
   fChain->SetBranchAddress("tau_ntracks_pt1", tau_ntracks_pt1, &b_tau_ntracks_pt1);
   fChain->SetBranchAddress("tau_L1trigger_match", tau_L1trigger_match, &b_tau_L1trigger_match);
   fChain->SetBranchAddress("tau_signalChargedHadrCands_size", tau_signalChargedHadrCands_size, &b_tau_signalChargedHadrCands_size);
   fChain->SetBranchAddress("tau_signalNeutralHadrCands_size", tau_signalNeutralHadrCands_size, &b_tau_signalNeutralHadrCands_size);
   fChain->SetBranchAddress("tau_signalGammaCands_size", tau_signalGammaCands_size, &b_tau_signalGammaCands_size);
   fChain->SetBranchAddress("tau_isolationChargedHadrCands_size", tau_isolationChargedHadrCands_size, &b_tau_isolationChargedHadrCands_size);
   fChain->SetBranchAddress("tau_isolationNeutralHadrCands_size", tau_isolationNeutralHadrCands_size, &b_tau_isolationNeutralHadrCands_size);
   fChain->SetBranchAddress("tau_isolationGammaCands_size", tau_isolationGammaCands_size, &b_tau_isolationGammaCands_size);
   fChain->SetBranchAddress("tau_genDecayMode_name", tau_genDecayMode_name, &b_tau_genDecayMode_name);
   fChain->SetBranchAddress("tau_genDecayMode", tau_genDecayMode, &b_tau_genDecayMode);
   fChain->SetBranchAddress("tau_decayMode_name", tau_decayMode_name, &b_tau_decayMode_name);
   fChain->SetBranchAddress("tau_decayMode", tau_decayMode, &b_tau_decayMode);
   fChain->SetBranchAddress("tau_constituents_count", tau_constituents_count, &b_tau_constituents_count);
   fChain->SetBranchAddress("tau_constituents_px", tau_constituents_px, &b_tau_constituents_px);
   fChain->SetBranchAddress("tau_constituents_py", tau_constituents_py, &b_tau_constituents_py);
   fChain->SetBranchAddress("tau_constituents_pz", tau_constituents_pz, &b_tau_constituents_pz);
   fChain->SetBranchAddress("tau_constituents_e", tau_constituents_e, &b_tau_constituents_e);
   fChain->SetBranchAddress("tau_constituents_mass", tau_constituents_mass, &b_tau_constituents_mass);
   fChain->SetBranchAddress("tau_constituents_charge", tau_constituents_charge, &b_tau_constituents_charge);
   fChain->SetBranchAddress("tau_constituents_vx", tau_constituents_vx, &b_tau_constituents_vx);
   fChain->SetBranchAddress("tau_constituents_vy", tau_constituents_vy, &b_tau_constituents_vy);
   fChain->SetBranchAddress("tau_constituents_vz", tau_constituents_vz, &b_tau_constituents_vz);
   fChain->SetBranchAddress("tau_constituents_pdgId", tau_constituents_pdgId, &b_tau_constituents_pdgId);
   fChain->SetBranchAddress("tau_helixparameters", tau_helixparameters, &b_tau_helixparameters);
   fChain->SetBranchAddress("tau_helixparameters_covar", tau_helixparameters_covar, &b_tau_helixparameters_covar);
   fChain->SetBranchAddress("tau_referencePoint", tau_referencePoint, &b_tau_referencePoint);
   fChain->SetBranchAddress("tau_Bfield", tau_Bfield, &b_tau_Bfield);
   fChain->SetBranchAddress("track_count", &track_count, &b_track_count);
   fChain->SetBranchAddress("track_px", track_px, &b_track_px);
   fChain->SetBranchAddress("track_py", track_py, &b_track_py);
   fChain->SetBranchAddress("track_pz", track_pz, &b_track_pz);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_mass", track_mass, &b_track_mass);
   fChain->SetBranchAddress("track_dxy", track_dxy, &b_track_dxy);
   fChain->SetBranchAddress("track_dxyerr", track_dxyerr, &b_track_dxyerr);
   fChain->SetBranchAddress("track_dz", track_dz, &b_track_dz);
   fChain->SetBranchAddress("track_dzerr", track_dzerr, &b_track_dzerr);
   fChain->SetBranchAddress("track_vx", track_vx, &b_track_vx);
   fChain->SetBranchAddress("track_vy", track_vy, &b_track_vy);
   fChain->SetBranchAddress("track_vz", track_vz, &b_track_vz);
   fChain->SetBranchAddress("track_ID", track_ID, &b_track_ID);
   fChain->SetBranchAddress("track_highPurity", track_highPurity, &b_track_highPurity);
   fChain->SetBranchAddress("pfmet_ex", &pfmet_ex, &b_pfmet_ex);
   fChain->SetBranchAddress("pfmet_ey", &pfmet_ey, &b_pfmet_ey);
   fChain->SetBranchAddress("pfmet_ez", &pfmet_ez, &b_pfmet_ez);
   fChain->SetBranchAddress("pfmet_pt", &pfmet_pt, &b_pfmet_pt);
   fChain->SetBranchAddress("pfmet_phi", &pfmet_phi, &b_pfmet_phi);
   fChain->SetBranchAddress("pfmet_sigxx", &pfmet_sigxx, &b_pfmet_sigxx);
   fChain->SetBranchAddress("pfmet_sigxy", &pfmet_sigxy, &b_pfmet_sigxy);
   fChain->SetBranchAddress("pfmet_sigyx", &pfmet_sigyx, &b_pfmet_sigyx);
   fChain->SetBranchAddress("pfmet_sigyy", &pfmet_sigyy, &b_pfmet_sigyy);
   fChain->SetBranchAddress("pfmet_sig", &pfmet_sig, &b_pfmet_sig);
   fChain->SetBranchAddress("genmet_ex", &genmet_ex, &b_genmet_ex);
   fChain->SetBranchAddress("genmet_ey", &genmet_ey, &b_genmet_ey);
   fChain->SetBranchAddress("pfmet_ex_JetEnUp", &pfmet_ex_JetEnUp, &b_pfmet_ex_JetEnUp);
   fChain->SetBranchAddress("pfmet_ey_JetEnUp", &pfmet_ey_JetEnUp, &b_pfmet_ey_JetEnUp);
   fChain->SetBranchAddress("pfmet_ex_JetEnDown", &pfmet_ex_JetEnDown, &b_pfmet_ex_JetEnDown);
   fChain->SetBranchAddress("pfmet_ey_JetEnDown", &pfmet_ey_JetEnDown, &b_pfmet_ey_JetEnDown);
   fChain->SetBranchAddress("pfmet_ex_UnclusteredEnUp", &pfmet_ex_UnclusteredEnUp, &b_pfmet_ex_UnclusteredEnUp);
   fChain->SetBranchAddress("pfmet_ey_UnclusteredEnUp", &pfmet_ey_UnclusteredEnUp, &b_pfmet_ey_UnclusteredEnUp);
   fChain->SetBranchAddress("pfmet_ex_UnclusteredEnDown", &pfmet_ex_UnclusteredEnDown, &b_pfmet_ex_UnclusteredEnDown);
   fChain->SetBranchAddress("pfmet_ey_UnclusteredEnDown", &pfmet_ey_UnclusteredEnDown, &b_pfmet_ey_UnclusteredEnDown);
   fChain->SetBranchAddress("pfmetcorr_ex", &pfmetcorr_ex, &b_pfmetcorr_ex);
   fChain->SetBranchAddress("pfmetcorr_ey", &pfmetcorr_ey, &b_pfmetcorr_ey);
   fChain->SetBranchAddress("pfmetcorr_ez", &pfmetcorr_ez, &b_pfmetcorr_ez);
   fChain->SetBranchAddress("pfmetcorr_pt", &pfmetcorr_pt, &b_pfmetcorr_pt);
   fChain->SetBranchAddress("pfmetcorr_phi", &pfmetcorr_phi, &b_pfmetcorr_phi);
   fChain->SetBranchAddress("pfmetcorr_sigxx", &pfmetcorr_sigxx, &b_pfmetcorr_sigxx);
   fChain->SetBranchAddress("pfmetcorr_sigxy", &pfmetcorr_sigxy, &b_pfmetcorr_sigxy);
   fChain->SetBranchAddress("pfmetcorr_sigyx", &pfmetcorr_sigyx, &b_pfmetcorr_sigyx);
   fChain->SetBranchAddress("pfmetcorr_sigyy", &pfmetcorr_sigyy, &b_pfmetcorr_sigyy);
   fChain->SetBranchAddress("pfmetcorr_sig", &pfmetcorr_sig, &b_pfmetcorr_sig);
   fChain->SetBranchAddress("pfmetcorr_ex_JetEnUp", &pfmetcorr_ex_JetEnUp, &b_pfmetcorr_ex_JetEnUp);
   fChain->SetBranchAddress("pfmetcorr_ey_JetEnUp", &pfmetcorr_ey_JetEnUp, &b_pfmetcorr_ey_JetEnUp);
   fChain->SetBranchAddress("pfmetcorr_ex_JetEnDown", &pfmetcorr_ex_JetEnDown, &b_pfmetcorr_ex_JetEnDown);
   fChain->SetBranchAddress("pfmetcorr_ey_JetEnDown", &pfmetcorr_ey_JetEnDown, &b_pfmetcorr_ey_JetEnDown);
   fChain->SetBranchAddress("pfmetcorr_ex_UnclusteredEnUp", &pfmetcorr_ex_UnclusteredEnUp, &b_pfmetcorr_ex_UnclusteredEnUp);
   fChain->SetBranchAddress("pfmetcorr_ey_UnclusteredEnUp", &pfmetcorr_ey_UnclusteredEnUp, &b_pfmetcorr_ey_UnclusteredEnUp);
   fChain->SetBranchAddress("pfmetcorr_ex_UnclusteredEnDown", &pfmetcorr_ex_UnclusteredEnDown, &b_pfmetcorr_ex_UnclusteredEnDown);
   fChain->SetBranchAddress("pfmetcorr_ey_UnclusteredEnDown", &pfmetcorr_ey_UnclusteredEnDown, &b_pfmetcorr_ey_UnclusteredEnDown);
   fChain->SetBranchAddress("pfmetcorr_ex_JetResUp", &pfmetcorr_ex_JetResUp, &b_pfmetcorr_ex_JetResUp);
   fChain->SetBranchAddress("pfmetcorr_ey_JetResUp", &pfmetcorr_ey_JetResUp, &b_pfmetcorr_ey_JetResUp);
   fChain->SetBranchAddress("pfmetcorr_ex_JetResDown", &pfmetcorr_ex_JetResDown, &b_pfmetcorr_ex_JetResDown);
   fChain->SetBranchAddress("pfmetcorr_ey_JetResDown", &pfmetcorr_ey_JetResDown, &b_pfmetcorr_ey_JetResDown);
   fChain->SetBranchAddress("pfmetcorr_ex_smeared", &pfmetcorr_ex_smeared, &b_pfmetcorr_ex_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_smeared", &pfmetcorr_ey_smeared, &b_pfmetcorr_ey_smeared);
   fChain->SetBranchAddress("pfmetcorr_ez_smeared", &pfmetcorr_ez_smeared, &b_pfmetcorr_ez_smeared);
   fChain->SetBranchAddress("pfmetcorr_pt_smeared", &pfmetcorr_pt_smeared, &b_pfmetcorr_pt_smeared);
   fChain->SetBranchAddress("pfmetcorr_phi_smeared", &pfmetcorr_phi_smeared, &b_pfmetcorr_phi_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_JetEnUp_smeared", &pfmetcorr_ex_JetEnUp_smeared, &b_pfmetcorr_ex_JetEnUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_JetEnUp_smeared", &pfmetcorr_ey_JetEnUp_smeared, &b_pfmetcorr_ey_JetEnUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_JetEnDown_smeared", &pfmetcorr_ex_JetEnDown_smeared, &b_pfmetcorr_ex_JetEnDown_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_JetEnDown_smeared", &pfmetcorr_ey_JetEnDown_smeared, &b_pfmetcorr_ey_JetEnDown_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_UnclusteredEnUp_smeared", &pfmetcorr_ex_UnclusteredEnUp_smeared, &b_pfmetcorr_ex_UnclusteredEnUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_UnclusteredEnUp_smeared", &pfmetcorr_ey_UnclusteredEnUp_smeared, &b_pfmetcorr_ey_UnclusteredEnUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_UnclusteredEnDown_smeared", &pfmetcorr_ex_UnclusteredEnDown_smeared, &b_pfmetcorr_ex_UnclusteredEnDown_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_UnclusteredEnDown_smeared", &pfmetcorr_ey_UnclusteredEnDown_smeared, &b_pfmetcorr_ey_UnclusteredEnDown_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_JetResUp_smeared", &pfmetcorr_ex_JetResUp_smeared, &b_pfmetcorr_ex_JetResUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_JetResUp_smeared", &pfmetcorr_ey_JetResUp_smeared, &b_pfmetcorr_ey_JetResUp_smeared);
   fChain->SetBranchAddress("pfmetcorr_ex_JetResDown_smeared", &pfmetcorr_ex_JetResDown_smeared, &b_pfmetcorr_ex_JetResDown_smeared);
   fChain->SetBranchAddress("pfmetcorr_ey_JetResDown_smeared", &pfmetcorr_ey_JetResDown_smeared, &b_pfmetcorr_ey_JetResDown_smeared);
   fChain->SetBranchAddress("puppimet_ex", &puppimet_ex, &b_puppimet_ex);
   fChain->SetBranchAddress("puppimet_ey", &puppimet_ey, &b_puppimet_ey);
   fChain->SetBranchAddress("puppimet_ez", &puppimet_ez, &b_puppimet_ez);
   fChain->SetBranchAddress("puppimet_pt", &puppimet_pt, &b_puppimet_pt);
   fChain->SetBranchAddress("puppimet_phi", &puppimet_phi, &b_puppimet_phi);
   fChain->SetBranchAddress("puppimet_sigxx", &puppimet_sigxx, &b_puppimet_sigxx);
   fChain->SetBranchAddress("puppimet_sigxy", &puppimet_sigxy, &b_puppimet_sigxy);
   fChain->SetBranchAddress("puppimet_sigyx", &puppimet_sigyx, &b_puppimet_sigyx);
   fChain->SetBranchAddress("puppimet_sigyy", &puppimet_sigyy, &b_puppimet_sigyy);
   fChain->SetBranchAddress("puppimet_ex_JetEnUp", &puppimet_ex_JetEnUp, &b_puppimet_ex_JetEnUp);
   fChain->SetBranchAddress("puppimet_ey_JetEnUp", &puppimet_ey_JetEnUp, &b_puppimet_ey_JetEnUp);
   fChain->SetBranchAddress("puppimet_ex_JetEnDown", &puppimet_ex_JetEnDown, &b_puppimet_ex_JetEnDown);
   fChain->SetBranchAddress("puppimet_ey_JetEnDown", &puppimet_ey_JetEnDown, &b_puppimet_ey_JetEnDown);
   fChain->SetBranchAddress("puppimet_ex_UnclusteredEnUp", &puppimet_ex_UnclusteredEnUp, &b_puppimet_ex_UnclusteredEnUp);
   fChain->SetBranchAddress("puppimet_ey_UnclusteredEnUp", &puppimet_ey_UnclusteredEnUp, &b_puppimet_ey_UnclusteredEnUp);
   fChain->SetBranchAddress("puppimet_ex_UnclusteredEnDown", &puppimet_ex_UnclusteredEnDown, &b_puppimet_ex_UnclusteredEnDown);
   fChain->SetBranchAddress("puppimet_ey_UnclusteredEnDown", &puppimet_ey_UnclusteredEnDown, &b_puppimet_ey_UnclusteredEnDown);
   fChain->SetBranchAddress("puppimet_ex_JetResUp", &puppimet_ex_JetResUp, &b_puppimet_ex_JetResUp);
   fChain->SetBranchAddress("puppimet_ey_JetResUp", &puppimet_ey_JetResUp, &b_puppimet_ey_JetResUp);
   fChain->SetBranchAddress("puppimet_ex_JetResDown", &puppimet_ex_JetResDown, &b_puppimet_ex_JetResDown);
   fChain->SetBranchAddress("puppimet_ey_JetResDown", &puppimet_ey_JetResDown, &b_puppimet_ey_JetResDown);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("genid1", &genid1, &b_genid1);
   fChain->SetBranchAddress("genx1", &genx1, &b_genx1);
   fChain->SetBranchAddress("genid2", &genid2, &b_genid2);
   fChain->SetBranchAddress("genx2", &genx2, &b_genx2);
   fChain->SetBranchAddress("genScale", &genScale, &b_genScale);
   fChain->SetBranchAddress("weightScale0", &weightScale0, &b_weightScale0);
   fChain->SetBranchAddress("weightScale1", &weightScale1, &b_weightScale1);
   fChain->SetBranchAddress("weightScale2", &weightScale2, &b_weightScale2);
   fChain->SetBranchAddress("weightScale3", &weightScale3, &b_weightScale3);
   fChain->SetBranchAddress("weightScale4", &weightScale4, &b_weightScale4);
   fChain->SetBranchAddress("weightScale5", &weightScale5, &b_weightScale5);
   fChain->SetBranchAddress("weightScale6", &weightScale6, &b_weightScale6);
   fChain->SetBranchAddress("weightScale7", &weightScale7, &b_weightScale7);
   fChain->SetBranchAddress("weightScale8", &weightScale8, &b_weightScale8);
   fChain->SetBranchAddress("weightPDFmax", &weightPDFmax, &b_weightPDFmax);
   fChain->SetBranchAddress("weightPDFmin", &weightPDFmin, &b_weightPDFmin);
   fChain->SetBranchAddress("weightPDFmean", &weightPDFmean, &b_weightPDFmean);
   fChain->SetBranchAddress("weightPDFup", &weightPDFup, &b_weightPDFup);
   fChain->SetBranchAddress("weightPDFdown", &weightPDFdown, &b_weightPDFdown);
   fChain->SetBranchAddress("weightPDFvar", &weightPDFvar, &b_weightPDFvar);
   fChain->SetBranchAddress("prefiringweight", &prefiringweight, &b_prefiringweight);
   fChain->SetBranchAddress("prefiringweightup", &prefiringweightup, &b_prefiringweightup);
   fChain->SetBranchAddress("prefiringweightdown", &prefiringweightdown, &b_prefiringweightdown);
   fChain->SetBranchAddress("numpileupinteractionsminus", &numpileupinteractionsminus, &b_numpileupinteractionsminus);
   fChain->SetBranchAddress("numpileupinteractions", &numpileupinteractions, &b_numpileupinteractions);
   fChain->SetBranchAddress("numpileupinteractionsplus", &numpileupinteractionsplus, &b_numpileupinteractionsplus);
   fChain->SetBranchAddress("numtruepileupinteractions", &numtruepileupinteractions, &b_numtruepileupinteractions);
   fChain->SetBranchAddress("gentau_count", &gentau_count, &b_gentau_count);
   fChain->SetBranchAddress("gentau_e", gentau_e, &b_gentau_e);
   fChain->SetBranchAddress("gentau_charge", gentau_charge, &b_gentau_charge);
   fChain->SetBranchAddress("gentau_px", gentau_px, &b_gentau_px);
   fChain->SetBranchAddress("gentau_py", gentau_py, &b_gentau_py);
   fChain->SetBranchAddress("gentau_pz", gentau_pz, &b_gentau_pz);
   fChain->SetBranchAddress("gentau_visible_e", gentau_visible_e, &b_gentau_visible_e);
   fChain->SetBranchAddress("gentau_visible_px", gentau_visible_px, &b_gentau_visible_px);
   fChain->SetBranchAddress("gentau_visible_py", gentau_visible_py, &b_gentau_visible_py);
   fChain->SetBranchAddress("gentau_visible_pz", gentau_visible_pz, &b_gentau_visible_pz);
   fChain->SetBranchAddress("gentau_visible_pt", gentau_visible_pt, &b_gentau_visible_pt);
   fChain->SetBranchAddress("gentau_visible_eta", gentau_visible_eta, &b_gentau_visible_eta);
   fChain->SetBranchAddress("gentau_visible_phi", gentau_visible_phi, &b_gentau_visible_phi);
   fChain->SetBranchAddress("gentau_visible_mass", gentau_visible_mass, &b_gentau_visible_mass);
   fChain->SetBranchAddress("gentau_visibleNoLep_e", gentau_visibleNoLep_e, &b_gentau_visibleNoLep_e);
   fChain->SetBranchAddress("gentau_visibleNoLep_px", gentau_visibleNoLep_px, &b_gentau_visibleNoLep_px);
   fChain->SetBranchAddress("gentau_visibleNoLep_py", gentau_visibleNoLep_py, &b_gentau_visibleNoLep_py);
   fChain->SetBranchAddress("gentau_visibleNoLep_pz", gentau_visibleNoLep_pz, &b_gentau_visibleNoLep_pz);
   fChain->SetBranchAddress("gentau_visibleNoLep_pt", gentau_visibleNoLep_pt, &b_gentau_visibleNoLep_pt);
   fChain->SetBranchAddress("gentau_visibleNoLep_eta", gentau_visibleNoLep_eta, &b_gentau_visibleNoLep_eta);
   fChain->SetBranchAddress("gentau_visibleNoLep_phi", gentau_visibleNoLep_phi, &b_gentau_visibleNoLep_phi);
   fChain->SetBranchAddress("gentau_visibleNoLep_mass", gentau_visibleNoLep_mass, &b_gentau_visibleNoLep_mass);
   fChain->SetBranchAddress("gentau_status", gentau_status, &b_gentau_status);
   fChain->SetBranchAddress("gentau_fromHardProcess", gentau_fromHardProcess, &b_gentau_fromHardProcess);
   fChain->SetBranchAddress("gentau_fromHardProcessBeforeFSR", gentau_fromHardProcessBeforeFSR, &b_gentau_fromHardProcessBeforeFSR);
   fChain->SetBranchAddress("gentau_isDecayedLeptonHadron", gentau_isDecayedLeptonHadron, &b_gentau_isDecayedLeptonHadron);
   fChain->SetBranchAddress("gentau_isDirectHadronDecayProduct", gentau_isDirectHadronDecayProduct, &b_gentau_isDirectHadronDecayProduct);
   fChain->SetBranchAddress("gentau_isDirectHardProcessTauDecayProduct", gentau_isDirectHardProcessTauDecayProduct, &b_gentau_isDirectHardProcessTauDecayProduct);
   fChain->SetBranchAddress("gentau_isDirectPromptTauDecayProduct", gentau_isDirectPromptTauDecayProduct, &b_gentau_isDirectPromptTauDecayProduct);
   fChain->SetBranchAddress("gentau_isDirectTauDecayProduct", gentau_isDirectTauDecayProduct, &b_gentau_isDirectTauDecayProduct);
   fChain->SetBranchAddress("gentau_isFirstCopy", gentau_isFirstCopy, &b_gentau_isFirstCopy);
   fChain->SetBranchAddress("gentau_isHardProcess", gentau_isHardProcess, &b_gentau_isHardProcess);
   fChain->SetBranchAddress("gentau_isHardProcessTauDecayProduct", gentau_isHardProcessTauDecayProduct, &b_gentau_isHardProcessTauDecayProduct);
   fChain->SetBranchAddress("gentau_isLastCopy", gentau_isLastCopy, &b_gentau_isLastCopy);
   fChain->SetBranchAddress("gentau_isLastCopyBeforeFSR", gentau_isLastCopyBeforeFSR, &b_gentau_isLastCopyBeforeFSR);
   fChain->SetBranchAddress("gentau_isPrompt", gentau_isPrompt, &b_gentau_isPrompt);
   fChain->SetBranchAddress("gentau_isPromptTauDecayProduct", gentau_isPromptTauDecayProduct, &b_gentau_isPromptTauDecayProduct);
   fChain->SetBranchAddress("gentau_isTauDecayProduct", gentau_isTauDecayProduct, &b_gentau_isTauDecayProduct);
   fChain->SetBranchAddress("gentau_decayMode", gentau_decayMode, &b_gentau_decayMode);
   fChain->SetBranchAddress("gentau_decayMode_name", gentau_decayMode_name, &b_gentau_decayMode_name);
   fChain->SetBranchAddress("gentau_mother", gentau_mother, &b_gentau_mother);
   fChain->SetBranchAddress("genparticles_lheHt", &genparticles_lheHt, &b_genparticles_lheHt);
   fChain->SetBranchAddress("genparticles_lheWPt", &genparticles_lheWPt, &b_genparticles_lheWPt);
   fChain->SetBranchAddress("genparticles_noutgoing", &genparticles_noutgoing, &b_genparticles_noutgoing);
   fChain->SetBranchAddress("genparticles_noutgoing_NLO", &genparticles_noutgoing_NLO, &b_genparticles_noutgoing_NLO);
   fChain->SetBranchAddress("genparticles_count", &genparticles_count, &b_genparticles_count);
   fChain->SetBranchAddress("genparticles_e", genparticles_e, &b_genparticles_e);
   fChain->SetBranchAddress("genparticles_px", genparticles_px, &b_genparticles_px);
   fChain->SetBranchAddress("genparticles_py", genparticles_py, &b_genparticles_py);
   fChain->SetBranchAddress("genparticles_pz", genparticles_pz, &b_genparticles_pz);
   fChain->SetBranchAddress("genparticles_vx", genparticles_vx, &b_genparticles_vx);
   fChain->SetBranchAddress("genparticles_vy", genparticles_vy, &b_genparticles_vy);
   fChain->SetBranchAddress("genparticles_vz", genparticles_vz, &b_genparticles_vz);
   fChain->SetBranchAddress("genparticles_pdgid", genparticles_pdgid, &b_genparticles_pdgid);
   fChain->SetBranchAddress("genparticles_status", genparticles_status, &b_genparticles_status);
   fChain->SetBranchAddress("genparticles_info", genparticles_info, &b_genparticles_info);
   fChain->SetBranchAddress("genparticles_fromHardProcess", genparticles_fromHardProcess, &b_genparticles_fromHardProcess);
   fChain->SetBranchAddress("genparticles_fromHardProcessBeforeFSR", genparticles_fromHardProcessBeforeFSR, &b_genparticles_fromHardProcessBeforeFSR);
   fChain->SetBranchAddress("genparticles_isDecayedLeptonHadron", genparticles_isDecayedLeptonHadron, &b_genparticles_isDecayedLeptonHadron);
   fChain->SetBranchAddress("genparticles_isDirectHadronDecayProduct", genparticles_isDirectHadronDecayProduct, &b_genparticles_isDirectHadronDecayProduct);
   fChain->SetBranchAddress("genparticles_isDirectHardProcessTauDecayProduct", genparticles_isDirectHardProcessTauDecayProduct, &b_genparticles_isDirectHardProcessTauDecayProduct);
   fChain->SetBranchAddress("genparticles_isDirectPromptTauDecayProduct", genparticles_isDirectPromptTauDecayProduct, &b_genparticles_isDirectPromptTauDecayProduct);
   fChain->SetBranchAddress("genparticles_isDirectTauDecayProduct", genparticles_isDirectTauDecayProduct, &b_genparticles_isDirectTauDecayProduct);
   fChain->SetBranchAddress("genparticles_isFirstCopy", genparticles_isFirstCopy, &b_genparticles_isFirstCopy);
   fChain->SetBranchAddress("genparticles_isHardProcess", genparticles_isHardProcess, &b_genparticles_isHardProcess);
   fChain->SetBranchAddress("genparticles_isHardProcessTauDecayProduct", genparticles_isHardProcessTauDecayProduct, &b_genparticles_isHardProcessTauDecayProduct);
   fChain->SetBranchAddress("genparticles_isLastCopy", genparticles_isLastCopy, &b_genparticles_isLastCopy);
   fChain->SetBranchAddress("genparticles_isLastCopyBeforeFSR", genparticles_isLastCopyBeforeFSR, &b_genparticles_isLastCopyBeforeFSR);
   fChain->SetBranchAddress("genparticles_isPrompt", genparticles_isPrompt, &b_genparticles_isPrompt);
   fChain->SetBranchAddress("genparticles_isPromptTauDecayProduct", genparticles_isPromptTauDecayProduct, &b_genparticles_isPromptTauDecayProduct);
   fChain->SetBranchAddress("genparticles_isTauDecayProduct", genparticles_isTauDecayProduct, &b_genparticles_isTauDecayProduct);
   fChain->SetBranchAddress("genparticles_mother", genparticles_mother, &b_genparticles_mother);
   fChain->SetBranchAddress("genjets_count", &genjets_count, &b_genjets_count);
   fChain->SetBranchAddress("genjets_e", genjets_e, &b_genjets_e);
   fChain->SetBranchAddress("genjets_px", genjets_px, &b_genjets_px);
   fChain->SetBranchAddress("genjets_py", genjets_py, &b_genjets_py);
   fChain->SetBranchAddress("genjets_pz", genjets_pz, &b_genjets_pz);
   fChain->SetBranchAddress("genjets_pt", genjets_pt, &b_genjets_pt);
   fChain->SetBranchAddress("genjets_eta", genjets_eta, &b_genjets_eta);
   fChain->SetBranchAddress("genjets_phi", genjets_phi, &b_genjets_phi);
   fChain->SetBranchAddress("genjets_pdgid", genjets_pdgid, &b_genjets_pdgid);
   fChain->SetBranchAddress("genjets_status", genjets_status, &b_genjets_status);
   fChain->SetBranchAddress("genjets_em_energy", genjets_em_energy, &b_genjets_em_energy);
   fChain->SetBranchAddress("genjets_had_energy", genjets_had_energy, &b_genjets_had_energy);
   fChain->SetBranchAddress("genjets_invisible_energy", genjets_invisible_energy, &b_genjets_invisible_energy);
   fChain->SetBranchAddress("genjets_auxiliary_energy", genjets_auxiliary_energy, &b_genjets_auxiliary_energy);
   fChain->SetBranchAddress("l1muon_count", &l1muon_count, &b_l1muon_count);
   fChain->SetBranchAddress("l1muon_px", l1muon_px, &b_l1muon_px);
   fChain->SetBranchAddress("l1muon_py", l1muon_py, &b_l1muon_py);
   fChain->SetBranchAddress("l1muon_pz", l1muon_pz, &b_l1muon_pz);
   fChain->SetBranchAddress("l1muon_pt", l1muon_pt, &b_l1muon_pt);
   fChain->SetBranchAddress("l1muon_ipt", l1muon_ipt, &b_l1muon_ipt);
   fChain->SetBranchAddress("l1muon_eta", l1muon_eta, &b_l1muon_eta);
   fChain->SetBranchAddress("l1muon_phi", l1muon_phi, &b_l1muon_phi);
   fChain->SetBranchAddress("l1muon_qual", l1muon_qual, &b_l1muon_qual);
   fChain->SetBranchAddress("l1muon_iso", l1muon_iso, &b_l1muon_iso);
   fChain->SetBranchAddress("l1muon_charge", l1muon_charge, &b_l1muon_charge);
   fChain->SetBranchAddress("l1muon_chargeValid", l1muon_chargeValid, &b_l1muon_chargeValid);
   fChain->SetBranchAddress("l1muon_muonIndex", l1muon_muonIndex, &b_l1muon_muonIndex);
   fChain->SetBranchAddress("l1muon_tag", l1muon_tag, &b_l1muon_tag);
   fChain->SetBranchAddress("l1muon_isoSum", l1muon_isoSum, &b_l1muon_isoSum);
   fChain->SetBranchAddress("l1muon_dPhiExtra", l1muon_dPhiExtra, &b_l1muon_dPhiExtra);
   fChain->SetBranchAddress("l1muon_dEtaExtra", l1muon_dEtaExtra, &b_l1muon_dEtaExtra);
   fChain->SetBranchAddress("l1muon_rank", l1muon_rank, &b_l1muon_rank);
   fChain->SetBranchAddress("l1muon_bx", l1muon_bx, &b_l1muon_bx);
   fChain->SetBranchAddress("l1egamma_count", &l1egamma_count, &b_l1egamma_count);
   fChain->SetBranchAddress("l1egamma_px", l1egamma_px, &b_l1egamma_px);
   fChain->SetBranchAddress("l1egamma_py", l1egamma_py, &b_l1egamma_py);
   fChain->SetBranchAddress("l1egamma_pz", l1egamma_pz, &b_l1egamma_pz);
   fChain->SetBranchAddress("l1egamma_pt", l1egamma_pt, &b_l1egamma_pt);
   fChain->SetBranchAddress("l1egamma_ipt", l1egamma_ipt, &b_l1egamma_ipt);
   fChain->SetBranchAddress("l1egamma_eta", l1egamma_eta, &b_l1egamma_eta);
   fChain->SetBranchAddress("l1egamma_phi", l1egamma_phi, &b_l1egamma_phi);
   fChain->SetBranchAddress("l1egamma_qual", l1egamma_qual, &b_l1egamma_qual);
   fChain->SetBranchAddress("l1egamma_iso", l1egamma_iso, &b_l1egamma_iso);
   fChain->SetBranchAddress("l1egamma_towerIEta", l1egamma_towerIEta, &b_l1egamma_towerIEta);
   fChain->SetBranchAddress("l1egamma_towerIPhi", l1egamma_towerIPhi, &b_l1egamma_towerIPhi);
   fChain->SetBranchAddress("l1egamma_rawEt", l1egamma_rawEt, &b_l1egamma_rawEt);
   fChain->SetBranchAddress("l1egamma_isoEt", l1egamma_isoEt, &b_l1egamma_isoEt);
   fChain->SetBranchAddress("l1egamma_footprintEt", l1egamma_footprintEt, &b_l1egamma_footprintEt);
   fChain->SetBranchAddress("l1egamma_nTT", l1egamma_nTT, &b_l1egamma_nTT);
   fChain->SetBranchAddress("l1egamma_shape", l1egamma_shape, &b_l1egamma_shape);
   fChain->SetBranchAddress("l1egamma_bx", l1egamma_bx, &b_l1egamma_bx);
   fChain->SetBranchAddress("l1tau_count", &l1tau_count, &b_l1tau_count);
   fChain->SetBranchAddress("l1tau_px", l1tau_px, &b_l1tau_px);
   fChain->SetBranchAddress("l1tau_py", l1tau_py, &b_l1tau_py);
   fChain->SetBranchAddress("l1tau_pz", l1tau_pz, &b_l1tau_pz);
   fChain->SetBranchAddress("l1tau_pt", l1tau_pt, &b_l1tau_pt);
   fChain->SetBranchAddress("l1tau_ipt", l1tau_ipt, &b_l1tau_ipt);
   fChain->SetBranchAddress("l1tau_eta", l1tau_eta, &b_l1tau_eta);
   fChain->SetBranchAddress("l1tau_phi", l1tau_phi, &b_l1tau_phi);
   fChain->SetBranchAddress("l1tau_qual", l1tau_qual, &b_l1tau_qual);
   fChain->SetBranchAddress("l1tau_iso", l1tau_iso, &b_l1tau_iso);
   fChain->SetBranchAddress("l1tau_towerIEta", l1tau_towerIEta, &b_l1tau_towerIEta);
   fChain->SetBranchAddress("l1tau_towerIPhi", l1tau_towerIPhi, &b_l1tau_towerIPhi);
   fChain->SetBranchAddress("l1tau_rawEt", l1tau_rawEt, &b_l1tau_rawEt);
   fChain->SetBranchAddress("l1tau_isoEt", l1tau_isoEt, &b_l1tau_isoEt);
   fChain->SetBranchAddress("l1tau_nTT", l1tau_nTT, &b_l1tau_nTT);
   fChain->SetBranchAddress("l1tau_hasEM", l1tau_hasEM, &b_l1tau_hasEM);
   fChain->SetBranchAddress("l1tau_isMerged", l1tau_isMerged, &b_l1tau_isMerged);
   fChain->SetBranchAddress("l1tau_bx", l1tau_bx, &b_l1tau_bx);
   fChain->SetBranchAddress("l1isotau_count", &l1isotau_count, &b_l1isotau_count);
   fChain->SetBranchAddress("l1isotau_e", &l1isotau_e, &b_l1isotau_e);
   fChain->SetBranchAddress("l1isotau_px", &l1isotau_px, &b_l1isotau_px);
   fChain->SetBranchAddress("l1isotau_py", &l1isotau_py, &b_l1isotau_py);
   fChain->SetBranchAddress("l1isotau_pz", &l1isotau_pz, &b_l1isotau_pz);
   fChain->SetBranchAddress("l1isotau_mass", &l1isotau_mass, &b_l1isotau_mass);
   fChain->SetBranchAddress("l1isotau_eta", &l1isotau_eta, &b_l1isotau_eta);
   fChain->SetBranchAddress("l1isotau_phi", &l1isotau_phi, &b_l1isotau_phi);
   fChain->SetBranchAddress("l1isotau_pt", &l1isotau_pt, &b_l1isotau_pt);
   fChain->SetBranchAddress("l1isotau_charge", &l1isotau_charge, &b_l1isotau_charge);
   fChain->SetBranchAddress("l1isotau_iso", &l1isotau_iso, &b_l1isotau_iso);
   fChain->SetBranchAddress("trigobject_count", &trigobject_count, &b_trigobject_count);
   fChain->SetBranchAddress("trigobject_px", trigobject_px, &b_trigobject_px);
   fChain->SetBranchAddress("trigobject_py", trigobject_py, &b_trigobject_py);
   fChain->SetBranchAddress("trigobject_pz", trigobject_pz, &b_trigobject_pz);
   fChain->SetBranchAddress("trigobject_pt", trigobject_pt, &b_trigobject_pt);
   fChain->SetBranchAddress("trigobject_eta", trigobject_eta, &b_trigobject_eta);
   fChain->SetBranchAddress("trigobject_phi", trigobject_phi, &b_trigobject_phi);
   fChain->SetBranchAddress("trigobject_filters", trigobject_filters, &b_trigobject_filters);
   fChain->SetBranchAddress("trigobject_isMuon", trigobject_isMuon, &b_trigobject_isMuon);
   fChain->SetBranchAddress("trigobject_isElectron", trigobject_isElectron, &b_trigobject_isElectron);
   fChain->SetBranchAddress("trigobject_isTau", trigobject_isTau, &b_trigobject_isTau);
   fChain->SetBranchAddress("trigobject_isJet", trigobject_isJet, &b_trigobject_isJet);
   fChain->SetBranchAddress("trigobject_isMET", trigobject_isMET, &b_trigobject_isMET);
   fChain->SetBranchAddress("run_hltnames", &run_hltnames, &b_run_hltnames);
   fChain->SetBranchAddress("run_hltfilters", &run_hltfilters, &b_run_hltfilters);
   fChain->SetBranchAddress("run_hltmufilters", &run_hltmufilters, &b_run_hltmufilters);
   fChain->SetBranchAddress("run_hltelectronfilters", &run_hltelectronfilters, &b_run_hltelectronfilters);
   fChain->SetBranchAddress("run_hlttaufilters", &run_hlttaufilters, &b_run_hlttaufilters);
   fChain->SetBranchAddress("run_hltphotonfilters", &run_hltphotonfilters, &b_run_hltphotonfilters);
   fChain->SetBranchAddress("run_hltjetfilters", &run_hltjetfilters, &b_run_hltjetfilters);
   fChain->SetBranchAddress("run_floattaudiscriminators", &run_floattaudiscriminators, &b_run_floattaudiscriminators);
   fChain->SetBranchAddress("run_binarytaudiscriminators", &run_binarytaudiscriminators, &b_run_binarytaudiscriminators);
   fChain->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators, &b_run_btagdiscriminators);
   fChain->SetBranchAddress("hltriggerresults", &hltriggerresults, &b_hltriggerresults);
   fChain->SetBranchAddress("hltriggerprescales", &hltriggerprescales, &b_hltriggerprescales);
   fChain->SetBranchAddress("hltriggerresultsV", &hltriggerresultsV, &b_hltriggerresultsV);
   fChain->SetBranchAddress("flags", &flags, &b_flags);
   fChain->SetBranchAddress("tau_againstElectronLooseMVA6", tau_againstElectronLooseMVA6, &b_tau_againstElectronLooseMVA6);
   fChain->SetBranchAddress("tau_againstElectronMVA6Raw", tau_againstElectronMVA6Raw, &b_tau_againstElectronMVA6Raw);
   fChain->SetBranchAddress("tau_againstElectronMVA6category", tau_againstElectronMVA6category, &b_tau_againstElectronMVA6category);
   fChain->SetBranchAddress("tau_againstElectronMediumMVA6", tau_againstElectronMediumMVA6, &b_tau_againstElectronMediumMVA6);
   fChain->SetBranchAddress("tau_againstElectronTightMVA6", tau_againstElectronTightMVA6, &b_tau_againstElectronTightMVA6);
   fChain->SetBranchAddress("tau_againstElectronVLooseMVA6", tau_againstElectronVLooseMVA6, &b_tau_againstElectronVLooseMVA6);
   fChain->SetBranchAddress("tau_againstElectronVTightMVA6", tau_againstElectronVTightMVA6, &b_tau_againstElectronVTightMVA6);
   fChain->SetBranchAddress("tau_againstMuonLoose3", tau_againstMuonLoose3, &b_tau_againstMuonLoose3);
   fChain->SetBranchAddress("tau_againstMuonTight3", tau_againstMuonTight3, &b_tau_againstMuonTight3);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBdR03oldDMwLTraw", tau_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw", tau_byIsolationMVArun2v1DBnewDMwLTraw, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw", tau_byIsolationMVArun2v1DBoldDMwLTraw, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWdR03oldDMwLTraw", tau_byIsolationMVArun2v1PWdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1PWdR03oldDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWnewDMwLTraw", tau_byIsolationMVArun2v1PWnewDMwLTraw, &b_tau_byIsolationMVArun2v1PWnewDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWoldDMwLTraw", tau_byIsolationMVArun2v1PWoldDMwLTraw, &b_tau_byIsolationMVArun2v1PWoldDMwLTraw);
   fChain->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT", tau_byLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT", tau_byLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT", tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWnewDMwLT", tau_byLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byLooseIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWoldDMwLT", tau_byLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT", tau_byMediumIsolationMVArun2v1DBnewDMwLT, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT", tau_byMediumIsolationMVArun2v1DBoldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT", tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWnewDMwLT", tau_byMediumIsolationMVArun2v1PWnewDMwLT, &b_tau_byMediumIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWoldDMwLT", tau_byMediumIsolationMVArun2v1PWoldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byPhotonPtSumOutsideSignalCone", tau_byPhotonPtSumOutsideSignalCone, &b_tau_byPhotonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT", tau_byTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT", tau_byTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT", tau_byTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWnewDMwLT", tau_byTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWoldDMwLT", tau_byTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT", tau_byVLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT", tau_byVLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT", tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWnewDMwLT", tau_byVLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWoldDMwLT", tau_byVLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT", tau_byVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT", tau_byVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT", tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWnewDMwLT", tau_byVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWoldDMwLT", tau_byVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT", tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT", tau_byVVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBoldDMwLT", tau_byVVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT", tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWnewDMwLT", tau_byVVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWoldDMwLT", tau_byVVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_chargedIsoPtSumdR03", tau_chargedIsoPtSumdR03, &b_tau_chargedIsoPtSumdR03);
   fChain->SetBranchAddress("tau_decayModeFinding", tau_decayModeFinding, &b_tau_decayModeFinding);
   fChain->SetBranchAddress("tau_decayModeFindingNewDMs", tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("tau_footprintCorrection", tau_footprintCorrection, &b_tau_footprintCorrection);
   fChain->SetBranchAddress("tau_footprintCorrectiondR03", tau_footprintCorrectiondR03, &b_tau_footprintCorrectiondR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeight", tau_neutralIsoPtSumWeight, &b_tau_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeightdR03", tau_neutralIsoPtSumWeightdR03, &b_tau_neutralIsoPtSumWeightdR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSumdR03", tau_neutralIsoPtSumdR03, &b_tau_neutralIsoPtSumdR03);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalCone", tau_photonPtSumOutsideSignalCone, &b_tau_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalConedR03", tau_photonPtSumOutsideSignalConedR03, &b_tau_photonPtSumOutsideSignalConedR03);
   fChain->SetBranchAddress("tau_puCorrPtSum", tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017", tau_byIsolationMVArun2017v1DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v1DBoldDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBnewDMwLTraw2017", tau_byIsolationMVArun2017v2DBnewDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBnewDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017", tau_byIsolationMVArun2017v2DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw2016", tau_byIsolationMVArun2v1DBnewDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw2016);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw2016", tau_byIsolationMVArun2v1DBoldDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw2016);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017", tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017", tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017", tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT2016", tau_byLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT2016", tau_byLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017", tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017", tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017", tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT2016", tau_byMediumIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT2016", tau_byMediumIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017", tau_byTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBnewDMwLT2017", tau_byTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017", tau_byTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT2016", tau_byTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT2016", tau_byTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017", tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017", tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017", tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016", tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016", tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017", tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017", tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017", tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT2016", tau_byVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT2016", tau_byVTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017", tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017", tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017", tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017", tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017", tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016", tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSeraw", tau_byDeepTau2017v2p1VSeraw, &b_tau_byDeepTau2017v2p1VSeraw);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSjetraw", tau_byDeepTau2017v2p1VSjetraw, &b_tau_byDeepTau2017v2p1VSjetraw);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSmuraw", tau_byDeepTau2017v2p1VSmuraw, &b_tau_byDeepTau2017v2p1VSmuraw);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSe", tau_byLooseDeepTau2017v2p1VSe, &b_tau_byLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSjet", tau_byLooseDeepTau2017v2p1VSjet, &b_tau_byLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSmu", tau_byLooseDeepTau2017v2p1VSmu, &b_tau_byLooseDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSe", tau_byMediumDeepTau2017v2p1VSe, &b_tau_byMediumDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSjet", tau_byMediumDeepTau2017v2p1VSjet, &b_tau_byMediumDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSmu", tau_byMediumDeepTau2017v2p1VSmu, &b_tau_byMediumDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSe", tau_byTightDeepTau2017v2p1VSe, &b_tau_byTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSjet", tau_byTightDeepTau2017v2p1VSjet, &b_tau_byTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSmu", tau_byTightDeepTau2017v2p1VSmu, &b_tau_byTightDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSe", tau_byVLooseDeepTau2017v2p1VSe, &b_tau_byVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSjet", tau_byVLooseDeepTau2017v2p1VSjet, &b_tau_byVLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSmu", tau_byVLooseDeepTau2017v2p1VSmu, &b_tau_byVLooseDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byVTightDeepTau2017v2p1VSe", tau_byVTightDeepTau2017v2p1VSe, &b_tau_byVTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVTightDeepTau2017v2p1VSjet", tau_byVTightDeepTau2017v2p1VSjet, &b_tau_byVTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVLooseDeepTau2017v2p1VSe", tau_byVVLooseDeepTau2017v2p1VSe, &b_tau_byVVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVVLooseDeepTau2017v2p1VSjet", tau_byVVLooseDeepTau2017v2p1VSjet, &b_tau_byVVLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVTightDeepTau2017v2p1VSe", tau_byVVTightDeepTau2017v2p1VSe, &b_tau_byVVTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVVTightDeepTau2017v2p1VSjet", tau_byVVTightDeepTau2017v2p1VSjet, &b_tau_byVVTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVVLooseDeepTau2017v2p1VSe", tau_byVVVLooseDeepTau2017v2p1VSe, &b_tau_byVVVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVVVLooseDeepTau2017v2p1VSjet", tau_byVVVLooseDeepTau2017v2p1VSjet, &b_tau_byVVVLooseDeepTau2017v2p1VSjet);

   fChain->SetBranchAddress("tau_MVADM2017v1", tau_MVADM2017v1, &b_tau_MVADM2017v1);
   fChain->SetBranchAddress("tau_MVADM2017v1DM0raw", tau_MVADM2017v1DM0raw, &b_tau_MVADM2017v1DM0raw);
   fChain->SetBranchAddress("tau_MVADM2017v1DM10raw", tau_MVADM2017v1DM10raw, &b_tau_MVADM2017v1DM10raw);
   fChain->SetBranchAddress("tau_MVADM2017v1DM11raw", tau_MVADM2017v1DM11raw, &b_tau_MVADM2017v1DM11raw);
   fChain->SetBranchAddress("tau_MVADM2017v1DM1raw", tau_MVADM2017v1DM1raw, &b_tau_MVADM2017v1DM1raw);
   fChain->SetBranchAddress("tau_MVADM2017v1DM2raw", tau_MVADM2017v1DM2raw, &b_tau_MVADM2017v1DM2raw);
   fChain->SetBranchAddress("tau_MVADM2017v1DMotherraw", tau_MVADM2017v1DMotherraw, &b_tau_MVADM2017v1DMotherraw);
   
   fChain->SetBranchAddress("TauSpinAngles_count", &TauSpinAngles_count, &b_TauSpinAngles_count);
   fChain->SetBranchAddress("TauSpinnerWeight", TauSpinnerWeight, &b_TauSpinnerWeight);

   fChain->SetBranchAddress("htxs_stage0cat",&htxs_stage0cat, &b_htxs_stage0cat);
   fChain->SetBranchAddress("htxs_stage1cat",&htxs_stage1cat , &b_htxs_stage1cat);
   fChain->SetBranchAddress("htxs_stage1p1cat",&htxs_stage1p1cat , &b_htxs_stage1p1cat);
   fChain->SetBranchAddress("htxs_higgsPt",&htxs_higgsPt , &b_htxs_higgsPt);
   fChain->SetBranchAddress("htxs_njets30", &htxs_njets30, &b_htxs_njets30);
   

   Notify();
}

Long64_t AC1B::GetEntries()
{
   if (!fChain) return 0;
   return fChain->GetEntries();
}

Bool_t AC1B::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AC1B::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AC1B::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AC1B_cxx
