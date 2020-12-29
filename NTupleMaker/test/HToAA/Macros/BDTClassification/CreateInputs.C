#include "HttStylesNew.cc"
#include "HtoH.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////           It contains 4 void functions:                                                                                      ////
////                                                                                                                              ////                               
////            1-) BDTClassificationSignal -> Classification Output of BDT for all Signal Samples                                ////
////            2-) BDTClassificationData   -> Classification Output of BDT for Background Model, Background Tests and Data in SR //// 
////            3-) InputsDataCards         -> Creation of DataCards and corresponding input files for limit computation          ////
////            4-) RunAllInputs            -> It runs for all samples the classifications and creates the datacards              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//********************************************************//
//************** Switch of input variables ***************//
//********************************************************//
bool MuLMuT_Mass_isIN            = false;
bool MuLMuT_DR_isIN              = true;
bool MuLMuT_DPhi_isIN            = false;
bool MuLTrk_Mass_isIN            = true;
bool MuLTrk_Pt_isIN              = true;
bool MuLTrk_DR_isIN              = true;
bool MuLTrkMET_DPhi_isIN         = false;
bool MuTTrk_Mass_isIN            = true;
bool MuTTrk_Pt_isIN              = true;
bool MuTTrk_DR_isIN              = true;
bool MuTTrkMET_DPhi_isIN         = true;
bool MuLTrkMuTTrk_Mass_isIN      = true;
bool MuLTrkMuTTrk_Pt_isIN        = false;
bool MuLTrkMuTTrkMET_Mass_isIN   = false;
bool MET_Pt_isIN                 = true;



void BDTClassificationSignal(TString mass="5", TString process="GGH"){

  float MuLMuT_Mass, MuLMuT_DR, MuLMuT_DPhi, MuLTrk_Mass, MuLTrk_Pt, // 5
        MuLTrk_DR, MuLTrkMET_DPhi, MuTTrk_Mass, MuTTrk_Pt, MuTTrk_DR, // 5
        MuTTrkMET_DPhi, MuLTrkMuTTrk_Mass, MuLTrkMuTTrk_Pt, MET_Pt, MuLTrkMuTTrkMET_Mass; // 5
  float Eventweight; // 1
  float Eventweight_TrkIso_Up,Eventweight_TrkIso_Down;
  float mvavalue;
  
  ///////////////////////////////// Region Strings ////////////////////////////////////////////
  int nRegions = 6;
  
  TString regions[6] = {"SemiIso", "LooseIso", "LooseSemiIso", "LeadingSemiIso", "LeadingLooseIso", "Sel"};

  
  TString sample_name = "SUSYGluGluToHToAA_AToTauTau_M-";
  if (process=="VBF") sample_name = "SUSYVBFToHToAA_AToTauTau_M-";
  if (process=="VH") sample_name = "SUSYVH_HToAA_AToTauTau_M-";
  if (process=="ttH") sample_name = "SUSYttH_HToAA_AToTauTau_M-";
  if (process=="MMTT") sample_name = "SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-";

  std::cout << std::endl<< std::endl;
  std::cout << "                    Start BDT Classification" << std::endl;
  std::cout << std::endl<< std::endl;
  std::cout << "                    ****     ***     *******                    " << std::endl;
  std::cout << "                    *   *    *   *      *                       " << std::endl;
  std::cout << "                    * **     *   *      *                       " << std::endl;
  std::cout << "                    *   *    *   *      *                       " << std::endl;
  std::cout << "                    ****     ***        *                       " << std::endl;
  std::cout << std::endl<< std::endl;
  
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << "***          Starting loop over the Signal Samples           ***" << std::endl;
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << std::endl<< std::endl;
  
  
  //********** Loading TMVA Reader ***********//
  
  TMVA::Reader * reader = new TMVA::Reader( "!Color:!Silent" );
  
  if (MuLMuT_Mass_isIN)
  reader->AddVariable("MuLMuT_Mass",&MuLMuT_Mass);
  if (MuLMuT_DR_isIN)
  reader->AddVariable("MuLMuT_DR",&MuLMuT_DR);
  if (MuLMuT_DPhi_isIN)
  reader->AddVariable("MuLMuT_DPhi",&MuLMuT_DPhi);
  if (MuLTrk_Mass_isIN)
  reader->AddVariable("MuLTrk_Mass",&MuLTrk_Mass);
  if (MuLTrk_Pt_isIN)
  reader->AddVariable("MuLTrk_Pt",&MuLTrk_Pt);
  if (MuLTrk_DR_isIN)
  reader->AddVariable("MuLTrk_DR",&MuLTrk_DR);
  if (MuLTrkMET_DPhi_isIN)
  reader->AddVariable("MuLTrkMET_DPhi",&MuLTrkMET_DPhi);
  if (MuTTrk_Mass_isIN)
  reader->AddVariable("MuTTrk_Mass",&MuLTrk_Mass);
  if (MuTTrk_Pt_isIN)
  reader->AddVariable("MuTTrk_Pt",&MuLTrk_Pt);
  if (MuTTrk_DR_isIN)
  reader->AddVariable("MuTTrk_DR",&MuLTrk_DR);
  if (MuTTrkMET_DPhi_isIN)
  reader->AddVariable("MuTTrkMET_DPhi",&MuLTrkMET_DPhi);
  if (MuLTrkMuTTrk_Mass_isIN)
  reader->AddVariable("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass);
  if (MuLTrkMuTTrk_Pt_isIN)
  reader->AddVariable("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt);
  if (MuLTrkMuTTrkMET_Mass_isIN)
  reader->AddVariable("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass);
  if (MET_Pt_isIN)
  reader->AddVariable("MET_Pt",&MET_Pt);
  
  if(atof(mass)<22.5) reader->BookMVA("BDT", "/nfs/dust/cms/user/perezdan/Run/Run2016/H2aa/H2aa_4tau/H2aa_4tau_MVA/MVA_BDT/ma"+mass+"/weights/TMVAClassification_BDT.weights.xml");
  //else reader->BookMVA("BDT", "/nfs/dust/cms/user/perezdan/Run/Run2016/H2aa/H2aa_4tau/H2aa_4tau_MVA/MVA_BDT/ma8/weights/TMVAClassification_BDT.weights.xml");
  
  //********** Input & Output Files ***********//
  
  TFile * fileIn  = new TFile("../"+sample_name+mass+".root");
  TFile * fileOut = new TFile("SUSY"+process+"_BDTOutput_M-"+mass+".root","recreate");
  
  std::cout<<std::endl;
  std::cout << "--- BDT Classification    : Using input file: " <<fileIn->GetName()<< std::endl;

   //********** Keep these Histos ***********//
   
   TH1D * histWeights = (TH1D*)fileIn->Get("histWeightsH");
   histWeights->Write("histWeightsH");
  
   TH1D * counter_FinalEventsH = (TH1D*)fileIn->Get("counter_FinalEventsH");
   counter_FinalEventsH->Write("counter_FinalEventsH");
   
   
   ////////////////////////////////////////////////////////
   //********** BDT Classification in Regions ***********//
   ////////////////////////////////////////////////////////
   
   for (int iregion=0;iregion<nRegions;iregion++) // Loop over regions
   {
		   
		////// Corresponding Tree //////
		TTree * tree = (TTree*)fileIn->Get("tree_"+regions[iregion]);
        tree->SetBranchAddress("MuLMuT_Mass",&MuLMuT_Mass);
        tree->SetBranchAddress("MuLMuT_DR",&MuLMuT_DR);
        tree->SetBranchAddress("MuLMuT_DPhi",&MuLMuT_DPhi);
        tree->SetBranchAddress("MuLTrk_Mass",&MuLTrk_Mass);
        tree->SetBranchAddress("MuLTrk_Pt",&MuLTrk_Pt);
        tree->SetBranchAddress("MuLTrk_DR",&MuLTrk_DR);
        tree->SetBranchAddress("MuLTrkMET_DPhi",&MuLTrkMET_DPhi);
        tree->SetBranchAddress("MuTTrk_Mass",&MuTTrk_Mass);
        tree->SetBranchAddress("MuTTrk_Pt",&MuTTrk_Pt);
        tree->SetBranchAddress("MuTTrk_DR",&MuTTrk_DR);
        tree->SetBranchAddress("MuTTrkMET_DPhi",&MuTTrkMET_DPhi);
        tree->SetBranchAddress("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass);
        tree->SetBranchAddress("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt);
        tree->SetBranchAddress("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass);
        tree->SetBranchAddress("MET_Pt",&MET_Pt);
        tree->SetBranchAddress("Eventweight",&Eventweight);
        if(regions[iregion]=="Sel")
        {
			tree->SetBranchAddress("Eventweight_TrkIso_Up",&Eventweight_TrkIso_Up);
            tree->SetBranchAddress("Eventweight_TrkIso_Down",&Eventweight_TrkIso_Down);
		}
        
        
        TH1D * Histo = new TH1D("MVABDTOutput"+regions[iregion]+"H","",2000,-1.,1.);
        TH1D * Histo_TrkIso_Up   = new TH1D("MVABDTOutput"+regions[iregion]+"_TrkIso_UpH","",2000,-1.,1.);
        TH1D * Histo_TrkIso_Down = new TH1D("MVABDTOutput"+regions[iregion]+"_TrkIso_DownH","",2000,-1.,1.);
               
        
        for(int ievt=0; ievt<tree->GetEntries();ievt++) // Loop over tree entries
        {
            tree->GetEntry(ievt); 
       
            //**************** Cuts on input variables ***************//
            bool VariableCut = MuLTrk_Mass       < 22. &&
                               MuTTrk_Mass       < 22. &&
                               MuLTrkMuTTrk_Mass < 125. ;

            if(!VariableCut) continue;
       
            mvavalue = reader->EvaluateMVA("BDT");

            Histo->Fill(mvavalue,Eventweight);
            
            if(regions[iregion]=="Sel")
            {
            Histo_TrkIso_Up->Fill(mvavalue,Eventweight_TrkIso_Up);
            Histo_TrkIso_Down->Fill(mvavalue,Eventweight_TrkIso_Down);
		    }
      
        } // End loop over tree entries
        
        Histo->Write("MVABDTOutput"+regions[iregion]+"H");
        
        if(regions[iregion]=="Sel")
        {
            Histo_TrkIso_Up->Write("MVABDTOutput"+regions[iregion]+"_TrkIso_UpH");
            Histo_TrkIso_Down->Write("MVABDTOutput"+regions[iregion]+"_TrkIso_DownH");
		}
        

   } // End loop over regions


  //********** End of Classification ***********//
  
  std::cout << "--- Created root file: \""<<fileOut->GetName()<<"\" containing the BDT output histograms" << std::endl;
  std::cout << "==> BDT Classification is done!" << std::endl << std::endl;

  fileOut->Close();

  std::cout << std::endl;
  std::cout << "End loop over the Signal Samples" << std::endl;
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << std::endl<< std::endl;
  

  
  std::cout << "                     End BDT Classification" << std::endl;
  
}


void BDTClassificationData(TString mass="5"){

  float MuLMuT_Mass, MuLMuT_DR, MuLMuT_DPhi, MuLTrk_Mass, MuLTrk_Pt, // 5
        MuLTrk_DR, MuLTrkMET_DPhi, MuTTrk_Mass, MuTTrk_Pt, MuTTrk_DR, // 5
        MuTTrkMET_DPhi, MuLTrkMuTTrk_Mass, MuLTrkMuTTrk_Pt, MET_Pt, MuLTrkMuTTrkMET_Mass; // 5
  float Eventweight; // 1
  float mvavalue;
  
  ///////////////////////////////// Region Strings ////////////////////////////////////////////
  int nRegions = 6;
  
  TString regions[6] = {"SemiIso", "LooseIso", "LooseSemiIso", "LeadingSemiIso", "LeadingLooseIso", "Sel"};

  std::cout << std::endl<< std::endl;
  std::cout << "                    Start BDT Classification" << std::endl;
  std::cout << std::endl<< std::endl;
  std::cout << "                    ****     ***     *******                    " << std::endl;
  std::cout << "                    *   *    *   *      *                       " << std::endl;
  std::cout << "                    * **     *   *      *                       " << std::endl;
  std::cout << "                    *   *    *   *      *                       " << std::endl;
  std::cout << "                    ****     ***        *                       " << std::endl;
  std::cout << std::endl<< std::endl;
  
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << "***          Starting loop over the Data           ***" << std::endl;
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << std::endl<< std::endl;
  
  
  //********** Loading TMVA Reader ***********//
  
  TMVA::Reader * reader = new TMVA::Reader( "!Color:!Silent" );
  
  if (MuLMuT_Mass_isIN)
  reader->AddVariable("MuLMuT_Mass",&MuLMuT_Mass);
  if (MuLMuT_DR_isIN)
  reader->AddVariable("MuLMuT_DR",&MuLMuT_DR);
  if (MuLMuT_DPhi_isIN)
  reader->AddVariable("MuLMuT_DPhi",&MuLMuT_DPhi);
  if (MuLTrk_Mass_isIN)
  reader->AddVariable("MuLTrk_Mass",&MuLTrk_Mass);
  if (MuLTrk_Pt_isIN)
  reader->AddVariable("MuLTrk_Pt",&MuLTrk_Pt);
  if (MuLTrk_DR_isIN)
  reader->AddVariable("MuLTrk_DR",&MuLTrk_DR);
  if (MuLTrkMET_DPhi_isIN)
  reader->AddVariable("MuLTrkMET_DPhi",&MuLTrkMET_DPhi);
  if (MuTTrk_Mass_isIN)
  reader->AddVariable("MuTTrk_Mass",&MuLTrk_Mass);
  if (MuTTrk_Pt_isIN)
  reader->AddVariable("MuTTrk_Pt",&MuLTrk_Pt);
  if (MuTTrk_DR_isIN)
  reader->AddVariable("MuTTrk_DR",&MuLTrk_DR);
  if (MuTTrkMET_DPhi_isIN)
  reader->AddVariable("MuTTrkMET_DPhi",&MuLTrkMET_DPhi);
  if (MuLTrkMuTTrk_Mass_isIN)
  reader->AddVariable("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass);
  if (MuLTrkMuTTrk_Pt_isIN)
  reader->AddVariable("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt);
  if (MuLTrkMuTTrkMET_Mass_isIN)
  reader->AddVariable("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass);
  if (MET_Pt_isIN)
  reader->AddVariable("MET_Pt",&MET_Pt);
  
  if(atof(mass)<22.5) reader->BookMVA("BDT", "/nfs/dust/cms/user/perezdan/Run/Run2016/H2aa/H2aa_4tau/H2aa_4tau_MVA/MVA_BDT/ma"+mass+"/weights/TMVAClassification_BDT.weights.xml");
  //else reader->BookMVA("BDT", "/nfs/dust/cms/user/perezdan/Run/Run2016/H2aa/H2aa_4tau/H2aa_4tau_MVA/MVA_BDT/ma8/weights/TMVAClassification_BDT.weights.xml");
  

  //********** Input & Output Files ***********//
  
  TFile * fileIn  = new TFile("../DoubleMuon_Run2016.root");
  TFile * fileOut = new TFile("DoubleMuon_BDTOutput_M-"+mass+".root","recreate");
  
  std::cout<<std::endl;
  std::cout << "--- BDT Classification    : Using input file: " <<fileIn->GetName()<< std::endl;

   //********** Keep these Histos ***********//
   
   TH1D * histWeightsData = (TH1D*)fileIn->Get("histWeightsH");
   histWeightsData->Write("histWeightsH");
  
   TH1D * counter_FinalEventsH = (TH1D*)fileIn->Get("counter_FinalEventsH");
   counter_FinalEventsH->Write("counter_FinalEventsH");


   ////////////////////////////////////////////////////////
   //********** BDT Classification in Regions ***********//
   ////////////////////////////////////////////////////////
   
   for (int iregion=0;iregion<nRegions;iregion++) // Loop over regions
   {
		   
		////// Corresponding Tree //////
		TTree * tree = (TTree*)fileIn->Get("tree_"+regions[iregion]);
        tree->SetBranchAddress("MuLMuT_Mass",&MuLMuT_Mass);
        tree->SetBranchAddress("MuLMuT_DR",&MuLMuT_DR);
        tree->SetBranchAddress("MuLMuT_DPhi",&MuLMuT_DPhi);
        tree->SetBranchAddress("MuLTrk_Mass",&MuLTrk_Mass);
        tree->SetBranchAddress("MuLTrk_Pt",&MuLTrk_Pt);
        tree->SetBranchAddress("MuLTrk_DR",&MuLTrk_DR);
        tree->SetBranchAddress("MuLTrkMET_DPhi",&MuLTrkMET_DPhi);
        tree->SetBranchAddress("MuTTrk_Mass",&MuTTrk_Mass);
        tree->SetBranchAddress("MuTTrk_Pt",&MuTTrk_Pt);
        tree->SetBranchAddress("MuTTrk_DR",&MuTTrk_DR);
        tree->SetBranchAddress("MuTTrkMET_DPhi",&MuTTrkMET_DPhi);
        tree->SetBranchAddress("MuLTrkMuTTrk_Mass",&MuLTrkMuTTrk_Mass);
        tree->SetBranchAddress("MuLTrkMuTTrk_Pt",&MuLTrkMuTTrk_Pt);
        tree->SetBranchAddress("MuLTrkMuTTrkMET_Mass",&MuLTrkMuTTrkMET_Mass);
        tree->SetBranchAddress("MET_Pt",&MET_Pt);
        tree->SetBranchAddress("Eventweight",&Eventweight);
    
        
        TH1D * Histo = new TH1D("MVABDTOutput"+regions[iregion]+"H","",2000,-1.,1.);
        
        for(int ievt=0; ievt<tree->GetEntries();ievt++) // Loop over tree entries
        {
            tree->GetEntry(ievt); 
       
            //**************** Cuts on input variables ***************//
            bool VariableCut = MuLTrk_Mass       < 22. &&
                               MuTTrk_Mass       < 22. &&
                               MuLTrkMuTTrk_Mass < 125. ;

            if(!VariableCut) continue;
       
            mvavalue = reader->EvaluateMVA("BDT");

            Histo->Fill(mvavalue);
      
        } // End loop over tree entries
        
        Histo->Write("MVABDTOutput"+regions[iregion]+"H");

        

   } // End loop over regions
   
  
  ////////////////////////////////////////////
  //********** Combining MC&DATA ***********//
  ////////////////////////////////////////////
  TFile * fileMC = new TFile("MC_BDTOutput_M-"+mass+".root");
  
  TH1D * EWKH = (TH1D*)fileMC->Get("MVABDTOutput_EWK_SelH");
  TH1D * TTH =  (TH1D*)fileMC->Get("MVABDTOutput_TT_SelH");
  TH1D * QCDH = (TH1D*)fileMC->Get("MVABDTOutput_QCD_SelH");
  TH1D * DYH =  (TH1D*)fileMC->Get("MVABDTOutput_DY_SelH");
  
  TH1D * LooseIsoH = (TH1D*)fileOut->Get("MVABDTOutputLooseIsoH")->Clone("LooseIsoH");
  
  TH1D * Histo = new TH1D("MVABDTOutputMixingMCandDATAH","",2000,-1.,1.);
 
  float TotBkgd = EWKH->GetSumOfWeights()+TTH->GetSumOfWeights()+QCDH->GetSumOfWeights()+DYH->GetSumOfWeights();
  
  EWKH->Scale(LooseIsoH->GetSumOfWeights()/TotBkgd);
  TTH->Scale(LooseIsoH->GetSumOfWeights()/TotBkgd);
  DYH->Scale(LooseIsoH->GetSumOfWeights()/TotBkgd);
  LooseIsoH->Scale(QCDH->GetSumOfWeights()/TotBkgd);
  
  Histo->Add(Histo,EWKH,1.,1.);
  Histo->Add(Histo,TTH,1.,1.);
  Histo->Add(Histo,DYH,1.,1.);
  Histo->Add(Histo,LooseIsoH,1.,1.);
  
  fileOut->cd();
  Histo->Write("MVABDTOutputMixingMCandDATAH");
  
  //********** End of Classification ***********//
  
  std::cout << "--- Created root file: \""<<fileOut->GetName()<<"\" containing the BDT output histograms" << std::endl;
  std::cout << "==> BDT Classification is done!" << std::endl << std::endl;
  fileOut->Close();

  
  std::cout << std::endl;
  std::cout << "End loop over the Data" << std::endl;
  std::cout << std::endl<< std::endl;
  std::cout << "                     End BDT Classification" << std::endl;
  
}

   

void InputsDataCards(TString mass="15") {
  
  double massD = atof( mass );

  SetStyle();

  std::cout << std::endl;
  std::cout << "Mass = " << mass << std::endl;
  std::cout << std::endl;
  
  double massTau = 1.777;
  double massMu  = 0.106;
  double massRatio = (massMu*massMu)/(massTau*massTau);
  double aF = 2*massTau/massD;
  double SF = 2*massRatio/TMath::Sqrt(1-aF*aF);

  double lumi = 35890;
  double xsecGGH = 48.52; ////48.52;43.92
  double xsecVBF = 3.779;
  double xsecVH  = 1.369 + 0.8824;
  double xsecTTH = 0.5065;
  double xsecMMTT = (xsecGGH+xsecVBF+xsecVH+xsecTTH) * SF;
  
    
  TFile * file = new TFile("DoubleMuon_BDTOutput_M-"+mass+".root");
  TFile * fileGGH = new TFile("SUSYGGH_BDTOutput_M-"+mass+".root");
  TFile * fileVBF = new TFile("SUSYVBF_BDTOutput_M-"+mass+".root");
  TFile * fileVH = new TFile("SUSYVH_BDTOutput_M-"+mass+".root");
  TFile * filettH = new TFile("SUSYttH_BDTOutput_M-"+mass+".root");
  TFile * fileMMTT = new TFile("SUSYMMTT_BDTOutput_M-"+mass+".root");
   

  //Shape (Nominal)
  TH1D * histBkgd_Old    = (TH1D*)file->Get("MVABDTOutputSemiIsoH");
  
  TH1D * histGGH_Old    = (TH1D*)fileGGH->Get("MVABDTOutputSelH");
  TH1D * histVBF_Old    = (TH1D*)fileVBF->Get("MVABDTOutputSelH");
  TH1D * histVH_Old    = (TH1D*)fileVH->Get("MVABDTOutputSelH");
  TH1D * histttH_Old    = (TH1D*)filettH->Get("MVABDTOutputSelH");
  TH1D * histMMTT_Old    = (TH1D*)fileMMTT->Get("MVABDTOutputSelH");
  
  //Shape (Uncertainty)
  TH1D * histBkgdUnc1_Old    = (TH1D*)file->Get("MVABDTOutputLeadingLooseIsoH");
  TH1D * histBkgdUnc2_Old    = (TH1D*)file->Get("MVABDTOutputLooseSemiIsoH");
  
  TH1D * histGGH_TrkIso_Up_Old    = (TH1D*)fileGGH->Get("MVABDTOutputSel_TrkIso_UpH");
  TH1D * histVBF_TrkIso_Up_Old    = (TH1D*)fileVBF->Get("MVABDTOutputSel_TrkIso_UpH");
  TH1D * histVH_TrkIso_Up_Old    = (TH1D*)fileVH->Get("MVABDTOutputSel_TrkIso_UpH");
  TH1D * histttH_TrkIso_Up_Old    = (TH1D*)filettH->Get("MVABDTOutputSel_TrkIso_UpH");
  TH1D * histMMTT_TrkIso_Up_Old    = (TH1D*)fileMMTT->Get("MVABDTOutputSel_TrkIso_UpH");
  
  TH1D * histGGH_TrkIso_Down_Old    = (TH1D*)fileGGH->Get("MVABDTOutputSel_TrkIso_DownH");
  TH1D * histVBF_TrkIso_Down_Old    = (TH1D*)fileVBF->Get("MVABDTOutputSel_TrkIso_DownH");
  TH1D * histVH_TrkIso_Down_Old    = (TH1D*)fileVH->Get("MVABDTOutputSel_TrkIso_DownH");
  TH1D * histttH_TrkIso_Down_Old    = (TH1D*)filettH->Get("MVABDTOutputSel_TrkIso_DownH");
  TH1D * histMMTT_TrkIso_Down_Old    = (TH1D*)fileMMTT->Get("MVABDTOutputSel_TrkIso_DownH"); 
  
  //Data Observed
  TH1D * hist_obs_Old    = (TH1D*)file->Get("MVABDTOutputSelH");
  
  
  ///////////// Finding Optimal Binning ////////////////
  histBkgd_Old->Rebin(100);
  histBkgdUnc1_Old->Rebin(100);
  histBkgdUnc2_Old->Rebin(100);
  
  int nBins = histBkgd_Old->GetNbinsX();
  float iniBin =-1.;
  float endBin = 1.;
  float bins[nBins] ; bins[0]=iniBin;
  for(int i=1;i<=nBins;i++) bins[i] = bins[0]+i*(endBin-iniBin)/nBins;
  
  for(int iBins=1;iBins<histBkgd_Old->GetNbinsX();iBins++)
  {
   if(histBkgd_Old->GetBinContent(iBins)==0. || histBkgdUnc1_Old->GetBinContent(iBins)==0. || histBkgdUnc2_Old->GetBinContent(iBins)==0.)
   {
    float temp=bins[iBins-(histBkgd_Old->GetNbinsX()-nBins)];
    for(int jArr=iBins-(histBkgd_Old->GetNbinsX()-nBins);jArr<nBins;jArr++)
    {
     bins[jArr]=bins[jArr+1];
    }
    bins[nBins]=temp;
    nBins=nBins-1;
   }
  }
  if(histBkgd_Old->GetBinContent(histBkgd_Old->GetNbinsX())==0. || histBkgdUnc1_Old->GetBinContent(histBkgd_Old->GetNbinsX())==0. || histBkgdUnc2_Old->GetBinContent(histBkgd_Old->GetNbinsX())==0. ) 
  {
	  bins[nBins-1]=bins[nBins];
	  nBins=nBins-1;
  }
  //// Merging First and Second Bin ////
  for(int i=1;i<=nBins;i++) bins[i] = bins[i+1];
  nBins=nBins-1;
   
  
  ///////////// Getting Normalizations ////////////////
  TH1D * histWeightsGGH = (TH1D*)fileGGH->Get("histWeightsH");
  double nGenGGH = histWeightsGGH->GetSumOfWeights();
  TH1D * histWeightsVBF = (TH1D*)fileVBF->Get("histWeightsH");
  double nGenVBF = histWeightsVBF->GetSumOfWeights();
  TH1D * histWeightsVH = (TH1D*)fileVH->Get("histWeightsH");
  double nGenVH = histWeightsVH->GetSumOfWeights();
  TH1D * histWeightsttH = (TH1D*)filettH->Get("histWeightsH");
  double nGenttH = histWeightsttH->GetSumOfWeights();
  TH1D * histWeightsMMTT = (TH1D*)fileMMTT->Get("histWeightsH");
  double nGenMMTT = histWeightsMMTT->GetSumOfWeights();
  
  
  double gghNorm = xsecGGH*lumi/nGenGGH;
  double vbfNorm = xsecVBF*lumi/nGenVBF;
  double vhNorm = xsecVH*lumi/nGenVH;
  double tthNorm = xsecTTH*lumi/nGenttH;
  double mmttNorm = xsecMMTT*lumi/nGenMMTT;
  
  double bkgNorm = hist_obs_Old->GetSumOfWeights();
  
  
  ///////////// Re-binning Histos ////////////////
  TH1D * histBkgd = (TH1D*)TH1DtoTH1D(histBkgd_Old,nBins,bins,true,"_new");
  TH1D * hist_obs = (TH1D*)TH1DtoTH1D(hist_obs_Old,nBins,bins,true,"_new_Data");
  TH1D * histGGH = (TH1D*)TH1DtoTH1D(histGGH_Old,nBins,bins,true,"_new_ggh");
  TH1D * histVBF = (TH1D*)TH1DtoTH1D(histVBF_Old,nBins,bins,true,"_new_vbf");
  TH1D * histVH = (TH1D*)TH1DtoTH1D(histVH_Old,nBins,bins,true,"_new_vh");
  TH1D * histttH = (TH1D*)TH1DtoTH1D(histttH_Old,nBins,bins,true,"_new_tth");
  TH1D * histMMTT = (TH1D*)TH1DtoTH1D(histMMTT_Old,nBins,bins,true,"_new_mmtt");
  
  
  TH1D * histBkgdUnc1 = (TH1D*)TH1DtoTH1D(histBkgdUnc1_Old,nBins,bins,true,"_new");
  TH1D * histBkgdUnc2 = (TH1D*)TH1DtoTH1D(histBkgdUnc2_Old,nBins,bins,true,"_new");
  
  TH1D * histGGH_TrkIso_Up = (TH1D*)TH1DtoTH1D(histGGH_TrkIso_Up_Old,nBins,bins,true,"_new_ggh");
  TH1D * histVBF_TrkIso_Up = (TH1D*)TH1DtoTH1D(histVBF_TrkIso_Up_Old,nBins,bins,true,"_new_vbf");
  TH1D * histVH_TrkIso_Up = (TH1D*)TH1DtoTH1D(histVH_TrkIso_Up_Old,nBins,bins,true,"_new_vh");
  TH1D * histttH_TrkIso_Up = (TH1D*)TH1DtoTH1D(histttH_TrkIso_Up_Old,nBins,bins,true,"_new_tth");
  TH1D * histMMTT_TrkIso_Up = (TH1D*)TH1DtoTH1D(histMMTT_TrkIso_Up_Old,nBins,bins,true,"_new_mmtt");
  
  TH1D * histGGH_TrkIso_Down = (TH1D*)TH1DtoTH1D(histGGH_TrkIso_Down_Old,nBins,bins,true,"_new_ggh");
  TH1D * histVBF_TrkIso_Down = (TH1D*)TH1DtoTH1D(histVBF_TrkIso_Down_Old,nBins,bins,true,"_new_vbf");
  TH1D * histVH_TrkIso_Down = (TH1D*)TH1DtoTH1D(histVH_TrkIso_Down_Old,nBins,bins,true,"_new_vh");
  TH1D * histttH_TrkIso_Down = (TH1D*)TH1DtoTH1D(histttH_TrkIso_Down_Old,nBins,bins,true,"_new_tth");
  TH1D * histMMTT_TrkIso_Down = (TH1D*)TH1DtoTH1D(histMMTT_TrkIso_Down_Old,nBins,bins,true,"_new_mmtt");
  
  
  ///////////// Normalization for Histos ////////////////
  histBkgd->Scale(1/histBkgd->GetSumOfWeights()*bkgNorm);
  
  histGGH->Scale(gghNorm);
  histVBF->Scale(vbfNorm);
  histVH->Scale(vhNorm);
  histttH->Scale(tthNorm);
  histMMTT->Scale(mmttNorm);
  
  histBkgdUnc1->Scale(1/histBkgdUnc1->GetSumOfWeights()*bkgNorm);
  histBkgdUnc2->Scale(1/histBkgdUnc2->GetSumOfWeights()*bkgNorm); 
  
  histGGH_TrkIso_Up->Scale(gghNorm);
  histVBF_TrkIso_Up->Scale(vbfNorm);
  histVH_TrkIso_Up->Scale(vhNorm);
  histttH_TrkIso_Up->Scale(tthNorm);
  histMMTT_TrkIso_Up->Scale(mmttNorm);
  
  histGGH_TrkIso_Down->Scale(gghNorm);
  histVBF_TrkIso_Down->Scale(vbfNorm);
  histVH_TrkIso_Down->Scale(vhNorm);
  histttH_TrkIso_Down->Scale(tthNorm);
  histMMTT_TrkIso_Down->Scale(mmttNorm);
  
  
  ///////////// Backgound Shape Uncertainty ////////////////
  TH1D * histBkgdUnc1Up = (TH1D*)histBkgdUnc1->Clone("histBkgdUnc1Up");
  TH1D * histBkgdUnc1Down = (TH1D*)histBkgdUnc1->Clone("histBkgdUnc1Down");
  TH1D * histBkgdUnc2Up = (TH1D*)histBkgdUnc2->Clone("histBkgdUnc2Up");
  TH1D * histBkgdUnc2Down = (TH1D*)histBkgdUnc2->Clone("histBkgdUnc2Down");
  
  for (int iB=1; iB<=nBins; ++iB) {
	  
	   if (histBkgdUnc1->GetBinContent(iB) !=0) histBkgdUnc1Up->SetBinContent(iB,histBkgd->GetBinContent(iB)*histBkgd->GetBinContent(iB)/histBkgdUnc1->GetBinContent(iB));
	   else { histBkgdUnc1Up->SetBinContent(iB,0.);}
	   if (histBkgdUnc2->GetBinContent(iB) !=0) histBkgdUnc2Up->SetBinContent(iB,histBkgd->GetBinContent(iB)*histBkgd->GetBinContent(iB)/histBkgdUnc2->GetBinContent(iB));
	   else { histBkgdUnc2Up->SetBinContent(iB,0.);}
  }
  histBkgdUnc1Up->Scale(1/histBkgdUnc1Up->GetSumOfWeights()*bkgNorm);
  histBkgdUnc1Down->Scale(1/histBkgdUnc1Down->GetSumOfWeights()*bkgNorm);
  histBkgdUnc2Up->Scale(1/histBkgdUnc2Up->GetSumOfWeights()*bkgNorm);
  histBkgdUnc2Down->Scale(1/histBkgdUnc2Down->GetSumOfWeights()*bkgNorm);
  
  ///////////// Printing Yields ////////////////
  std::cout << "Bkg  Norm = " << bkgNorm << std::endl;
  std::cout << "ggH  Norm = " << histGGH->GetSumOfWeights() << std::endl;
  std::cout << "VBF  Norm = " << histVBF->GetSumOfWeights() << std::endl;
  std::cout << "VH  Norm = " << histVH->GetSumOfWeights() << std::endl;
  std::cout << "ttH  Norm = " << histttH->GetSumOfWeights() << std::endl;
  std::cout << "MMTT  Norm = " << histMMTT->GetSumOfWeights() << std::endl;
  std::cout << "Data      : " << histBkgd->GetSumOfWeights() << std::endl;
  
  
  ///////////// Creating Datacards & Root Files for limit extraction ////////////////
  TString DirName = "DataCards/";
  TString BaseName = "haa-13TeV_2016_ma" + mass;
  TString rootFileName = BaseName+".root";
  
  TFile * fileInputs = new TFile(DirName+rootFileName,"recreate");
  hist_obs->Write("data_obs");
  histBkgd->Write("bkgd");
  histGGH->Write("ggh");
  histVBF->Write("vbf");
  histVH->Write("vh");
  histttH->Write("tth");
  histMMTT->Write("mmtt");
  histBkgdUnc1Up->Write("bkgd_CMS_bdtbkg_unc1Up");
  histBkgdUnc1Down->Write("bkgd_CMS_bdtbkg_unc1Down");
  histBkgdUnc2Up->Write("bkgd_CMS_bdtbkg_unc2Up");
  histBkgdUnc2Down->Write("bkgd_CMS_bdtbkg_unc2Down");
  histGGH_TrkIso_Up->Write("ggh_CMS_TrkIso_uncUp");
  histVBF_TrkIso_Up->Write("vbf_CMS_TrkIso_uncUp");
  histVH_TrkIso_Up->Write("vh_CMS_TrkIso_uncUp");
  histttH_TrkIso_Up->Write("tth_CMS_TrkIso_uncUp");
  histMMTT_TrkIso_Down->Write("mmtt_CMS_TrkIso_uncUp");
  histGGH_TrkIso_Down->Write("ggh_CMS_TrkIso_uncDown");
  histVBF_TrkIso_Down->Write("vbf_CMS_TrkIso_uncDown");
  histVH_TrkIso_Down->Write("vh_CMS_TrkIso_uncDown");
  histttH_TrkIso_Down->Write("tth_CMS_TrkIso_uncDown");
  histMMTT_TrkIso_Down->Write("mmtt_CMS_TrkIso_uncDown");
  
  fileInputs->Close();

  ostringstream str;
  str << DirName+BaseName << ".txt";
  string nn = str.str();
  const char * p = nn.c_str();

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "observation " << hist_obs->GetSumOfWeights() << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * " << rootFileName << "  $PROCESS    $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin";
  textFile << "             haa       haa        haa        haa        haa        haa  "   << std::endl;
  textFile << "process      mmtt      tth	      vh	    vbf	       ggh       bkgd" << std::endl;
  textFile << "process       -4       -3         -2         -1          0         1" << std::endl;
  textFile << "rate      " 
	   << histMMTT->GetSumOfWeights() << "      " 
	   << histttH->GetSumOfWeights() << "      "
	   << histVH->GetSumOfWeights() << "      " 
	   << histVBF->GetSumOfWeights() << "      "
	   << histGGH->GetSumOfWeights() << "      " 
	   << bkgNorm << std::endl;
  textFile << "-----------------------------" << std::endl;
  textFile << "CMS_lumi                    lnN   1.025  1.025   1.025   1.025   1.025      -" << std::endl;
  textFile << "CMS_eff_m                   lnN   1.04    1.04    1.04    1.04    1.04      -" << std::endl;
  textFile << "CMS_TrkIso_unc             shape  1.00    1.00    1.00    1.00    1.00      -" << std::endl;
  textFile << "CMS_bdtbkg_unc1             shape    -       -       -       -       -    1.00" << std::endl;
  //textFile << "CMS_bdtbkg_unc2             shape    -       -       -       -       -    1.00" << std::endl;
  
  textFile << "QCDScale_ggH                lnN   1.046/0.933   -       -       -  1.046/0.933  -" << std::endl;
  textFile << "QCDScale_vbf                lnN      -          -       -  1.004/0.997  -       -" << std::endl;
  textFile << "QCDScale_vh                 lnN      -          - 1.018/0.983   -       -       -" << std::endl;
  textFile << "QCDScale_ttH                lnN      -    1.058/0.908   -    -       -       -" << std::endl;

  textFile << "PDF_ggh                     lnN   1.032      -       -       -    1.032      -" << std::endl;
  textFile << "PDF_vbf                     lnN      -       -       -    1.021      -       -" << std::endl;
  textFile << "PDF_vh                      lnN      -       -    1.018      -       -       -" << std::endl;
  textFile << "PDF_tth                     lnN      -    1.036      -       -       -       -" << std::endl;

  textFile << "bkgNorm  rateParam  haa  bkgd  1  [0.5,1.5]" << std::endl;
  textFile << std::endl;

  

}

void RunAllInputs(){
  
  int ngenmass = 15;
  TString genmassString[] = {"4","5","6","7","8","9","10","11","12","13","14","15","17","19","21"};
  int nsigprocesses = 5;
  TString processes[]={"GGH","VBF","VH","ttH","MMTT"};

  //********************************************************//
  //**********Loop over the generated mass points***********//
  //********************************************************//

  for (int igenMass=0;igenMass<ngenmass;igenMass++) // Loop over the generated mass points
  {
	 for (int ipro=0;ipro<nsigprocesses;ipro++) // Loop over signal processes
     { 
       BDTClassificationSignal(genmassString[igenMass],processes[ipro]);
     }
     
  BDTClassificationData(genmassString[igenMass]);
  InputsDataCards(genmassString[igenMass]);
	 
  }



}


