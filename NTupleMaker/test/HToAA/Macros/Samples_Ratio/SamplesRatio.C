#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"
#include "Plotting.h"
#include "Plotting_Style.h"


void CreateHistos()
{
	
  
  float MuLMuT_Mass, MuLMuT_DR, MuLMuT_DPhi, MuLTrk_Mass, MuLTrk_Pt, // 5
        MuLTrk_DR, MuLTrkMET_DPhi, MuTTrk_Mass, MuTTrk_Pt, MuTTrk_DR, // 5
        MuTTrkMET_DPhi, MuLTrkMuTTrk_Mass, MuLTrkMuTTrk_Pt, MET_Pt, MuLTrkMuTTrkMET_Mass; // 5
  float Eventweight; // 1

  
  float	xLow[]     = {0., 0, 0., 0., 0.,
                      0., 0., 0., 0., 0.,
                      0., 0., 0., 0., 0.};

  float xUp[]     = {400., 10., 4., 22., 200.,
                     1.5, 4, 22., 200., 1.5,
                     4., 300., 300., 300. , 300.};
                     
  
  ///////////////////////////////// Sample Strings //////////////////////////////////////////// 
  int nSamples = 14;
                     
  TString samples[14] = {"WW_13TeV-pythia8",                                             // (0)
			             "WZ_13TeV-pythia8",                                             // (1) 
			             "ZZ_13TeV-pythia8",                                             // (2)
			             "WJetsToLNu_13TeV-madgraphMLM",                                 // (3)
			             "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg",         // (4)
			             "ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg",             // (5)
			             "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg",                // (6)
			             "ST_tW_top_5f_inclusiveDecays_13TeV-powheg",                    // (7)
			             "TT_13TeV-powheg-pythia8",                                      // (8)
			             "DYJetsToLL_M-10to50_13TeV-madgraphMLM",                        // (9)
			             "DYJetsToLL_M-50_13TeV-madgraphMLM",                            // (10)
			             "QCD_Pt-20toInf_MuEnrichedPt15_13TeV",                          // (11)
			             "SUSYGluGluToHToAA_AToTauTau_M-10",                             // (12)
			             "DoubleMuon_Run2016"                                            // (13)
  };
  
  
  
  ///////////////////////////////// Region Strings ////////////////////////////////////////////
  int nRegions = 6;

  TString regions[6] = {"SemiIso", "LooseIso", "LooseSemiIso", "LeadingSemiIso", "LeadingLooseIso", "Sel"};
  
  
  
  ///////////////////////////////// Variable Strings ////////////////////////////////////////////
  int nVariables = 15;
  
  TString variables[15]= {"MuLMuT_Mass", "MuLMuT_DR", "MuLMuT_DPhi", "MuLTrk_Mass", "MuLTrk_Pt",                        // 5
                          "MuLTrk_DR", "MuLTrkMET_DPhi", "MuTTrk_Mass", "MuTTrk_Pt", "MuTTrk_DR",                       // 10
                          "MuTTrkMET_DPhi", "MuLTrkMuTTrk_Mass", "MuLTrkMuTTrk_Pt", "MET_Pt", "MuLTrkMuTTrkMET_Mass"};  // 15
                          
                      
  ///////////////////////////////// Directory String ////////////////////////////////////////////                    
  TString dir = "/nfs/dust/cms/user/perezdan/Run/Run2016/H2aa/H2aa_4tau/H2aa_4tau_MVA/";
  
  
  std::cout << std::endl<< std::endl;
  std::cout << "Creating Histograms is starting ... " << std::endl;
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << std::endl<< std::endl;
  
  
  for (int isample=0;isample<nSamples;isample++) // Loop over samples
  { 
	   //// Input and Outout files ////
       TFile * fileIn  = new TFile(dir+samples[isample]+".root");
       TFile * fileOut = new TFile(samples[isample]+".root","recreate");
       
       //// Save this histos ////
       TH1D * histWeights = (TH1D*)fileIn->Get("histWeightsH");
       histWeights->Write("histWeightsH");
  
       TH1D * counter_FinalEventsH = (TH1D*)fileIn->Get("counter_FinalEventsH");
       counter_FinalEventsH->Write("counter_FinalEventsH");
       
       std::cout<<std::endl;
       std::cout << "-----> Using input file: " <<fileIn->GetName()<< std::endl;
       
       for (int iregion=0;iregion<nRegions;iregion++) // Loop over regions
       {
		   
		   ////// Corresponding Tree //////
		   TTree* tree = (TTree*)fileIn->Get("tree_"+regions[iregion]);
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
           
           
           TH1D * Histos[nVariables];
           for (int ivariable=0;ivariable<nVariables;ivariable++) // Loop over variables
           {
			
			   Histos[ivariable] = new TH1D(variables[ivariable]+"_"+regions[iregion]+"H","",40,xLow[ivariable],xUp[ivariable]);
			   
			        
		   } // End loop over variables
		   
		   
		   for(int ievt=0; ievt<tree->GetEntries();ievt++) // Loop over tree entries
           {
              
              tree->GetEntry(ievt); 
              
              //**************** Cuts on input variables ***************//
              bool VariableCut = MuLTrk_Mass       < 22. &&
                                 MuTTrk_Mass       < 22. &&
                                 MuLTrkMuTTrk_Mass < 125. ;
              
              if (!VariableCut) continue;

              Histos[0]->Fill(MuLMuT_Mass,Eventweight);
              Histos[1]->Fill(MuLMuT_DR,Eventweight);
              Histos[2]->Fill(MuLMuT_DPhi,Eventweight);
              Histos[3]->Fill(MuLTrk_Mass,Eventweight);
              Histos[4]->Fill(MuLTrk_Pt,Eventweight);
              Histos[5]->Fill(MuLTrk_DR,Eventweight);
              Histos[6]->Fill(MuLTrkMET_DPhi,Eventweight);
              Histos[7]->Fill(MuTTrk_Mass,Eventweight);
              Histos[8]->Fill(MuTTrk_Pt,Eventweight);
              Histos[9]->Fill(MuTTrk_DR,Eventweight);
              Histos[10]->Fill(MuTTrkMET_DPhi,Eventweight);
              Histos[11]->Fill(MuLTrkMuTTrk_Mass,Eventweight);
              Histos[12]->Fill(MuLTrkMuTTrk_Pt,Eventweight);
              Histos[13]->Fill(MET_Pt,Eventweight);
              Histos[14]->Fill(MuLTrkMuTTrkMET_Mass,Eventweight);
             
      
           } // End loop over tree entries
           
           for (int ivariable=0;ivariable<nVariables;ivariable++) Histos[ivariable]->Write(variables[ivariable]+"_"+regions[iregion]+"H");
           
           
	   } // End loop over regions
	   
	   std::cout << std::endl;
       std::cout << "-----> Created root file: "<<fileOut->GetName()<<"  containing the output histograms" << std::endl;
       
  } // End loop over samples
  
  
  std::cout << std::endl<< std::endl;
  std::cout << "Creating Histograms is finished !!" << std::endl;
  std::cout << "--- ------------------------------------------------------------" << std::endl;
  std::cout << std::endl<< std::endl;
  
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////*********************************** Histogram Names and Labeling ******************************************//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TString VarName[]= {"MuLMuT_Mass", "MuLMuT_DR", "MuLMuT_DPhi", "MuLTrk_Mass", "MuLTrk_Pt", // 5
//                      "MuLTrk_DR", "MuLTrkMET_DPhi", "MuTTrk_Mass", "MuTTrk_Pt", "MuTTrk_DR", // 5
//                      "MuTTrkMET_DPhi", "MuLTrkMuTTrk_Mass", "MuLTrkMuTTrk_Pt", "MET_Pt", "MuLTrkMuTTrkMET_Mass"}; 
//
//  TString xTitle[] = {"m_{#mu_{L}, #mu_{T}} [GeV]", "#DeltaR_{#mu_{L}, #mu_{T}}", "#Delta#phi_{#mu_{L}, #mu_{T}}", "m_{#mu_{L}, trk_{L}} [GeV]", "p_{T}^{#mu_{L}, trk_{L}} [GeV]", // 5
//                      "#DeltaR_{#mu_{L}, trk_{L}}", "#Delta#phi_{#mu_{L}trk_{L}, E^{miss}_{T}}", "m_{#mu_{T}, trk_{T}} [GeV]", "p_{T}^{#mu_{T}, trk_{T}} [GeV]", // 5
//                      "#DeltaR_{#mu_{T}, trk_{T}}", "#Delta#phi_{#mu_{T}trk_{T}, E^{miss}_{T}}", "m_{#mu_{L}, trk_{L}, #mu_{T}, trk_{T}} [GeV]", "p_{T}^{#mu_{L}, trk_{L}, #mu_{T}, trk_{T}} [GeV]", "E^{miss}_{T} [GeV]", "m_{#mu_{L}, trk_{L}, #mu_{T}, trk_{T}, E^{miss}_{T}} [GeV]"};
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//TString regions[6] = {"SemiIso", "LooseIso", "LooseSemiIso", "LeadingSemiIso", "LeadingLooseIso", "Sel"};

void PlotSamplesRatio(TString RegionName = "Sel",                      // region name
                      TString VariableName = "MuLTrk_DR",                  // variable name
 		              TString xTitle = "#DeltaR_{#mu_{L}, trk_{L}}",          // title of x axis
		              TString ytitle = "Events / Bin",                       // title of y axis		             
		              float xLower = 0.0,                                    // lower boundary of x axis 
		              float xUpper = 1.5,                                   // upper boundary of x axis
		              float yLower = 10,                                    // lower boundary of y axis (in case when plotting y axis in log scale)
		              bool logY = true,                                      // log or linear scale of Y axis
		              bool drawLeg = true)                                   // draw or not legend
{ 

  gROOT->SetBatch();



  float qcdScale = 0.5;
  float xsecGGH = 48.52 * 0.20 * 100;
  float lumi = 35890;

  TFile * file_data = new TFile("DoubleMuon_Run2016.root");
  
  
  int nMCSamples = 13;
  
  TString mc_samples[13] = {"WW_13TeV-pythia8",                                             // (0)
	  		                "WZ_13TeV-pythia8",                                             // (1) 
			                "ZZ_13TeV-pythia8",                                             // (2)
			                "WJetsToLNu_13TeV-madgraphMLM",                                 // (3)
			                "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powheg",         // (4)
			                "ST_t-channel_top_4f_inclusiveDecays_13TeV-powheg",             // (5)
			                "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg",                // (6)
			                "ST_tW_top_5f_inclusiveDecays_13TeV-powheg",                    // (7)
			                "TT_13TeV-powheg-pythia8",                                      // (8)
			                "DYJetsToLL_M-10to50_13TeV-madgraphMLM",                        // (9)
			                "DYJetsToLL_M-50_13TeV-madgraphMLM",                            // (10)
			                "QCD_Pt-20toInf_MuEnrichedPt15_13TeV",                          // (11)
			                "SUSYGluGluToHToAA_AToTauTau_M-10"                              // (12)
  };

  
  float xsec[13] = {12.14,                             // WW (0)
                    27.60,                             // WZ (1)
                    12.14,                             // ZZ (2)
                    61526.7,                           // WJets                (3)
                    80.95,                             // ST_t-channel_antitop (4)
                    136.02,                            // ST_t-channel_top     (5)
                    35.6,                              // ST_tW_antitop        (6)
                    35.6,                              // ST_tW_top            (7)
                    831.8,                             // TTbar                (8)
                    18610,                             // DYJetsToLL_M-10to50  (9)
                    6077.22,                           // DYJetsToLL_M-50      (10)
                    720648000 * 0.00042,               // QCD                  (11)
                    xsecGGH                            //ma=10                 (12)
  };

  xsec[11]=xsec[11]*qcdScale;
  
  
  TH1D * DataOldH = (TH1D*)file_data->Get(VariableName+"_"+RegionName+"H");

  ////// **** Rebining **** ////////
  int nBins = DataOldH->GetNbinsX();
  float xMin = DataOldH->GetBinLowEdge(1);
  float xMax = DataOldH->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << DataOldH->GetName() << " : " << "nbins = " << nBins
	        << " , min = " << xMin
	        << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  float bins[40];
  int nBinsNew = nBins;
  std::cout << "New number of bins : ";
  std::cin >> nBinsNew;

  if (nBins%nBinsNew!=0) { 
    std::cout << "new number of bins = " << nBinsNew 
	          << "  not multiple of " << nBins << std::endl;
    return;
  }
  
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB) bins[iB] = xMin + float(iB)*binWidth;


  ////// **** Histograms **** ////////
  TH1D * DataH = TH1DtoTH1D(DataOldH,nBinsNew,bins,true,"_Data_new");

  TH1D * EWKH = new TH1D("EWKH","",nBinsNew,bins);
  TH1D * TTH  = new TH1D("TTH", "",nBinsNew,bins);
  TH1D * QCDH = new TH1D("QCDH","",nBinsNew,bins);
  TH1D * DYH   = new TH1D("DYH",  "",nBinsNew,bins);
  TH1D * SUSYH   = new TH1D("SUSYH",  "",nBinsNew,bins);


  ////// **** MC Histograms **** ////////
  for (int isample=0; isample<nMCSamples; isample++) 
  {
	  
    TFile * fileMC = new TFile(mc_samples[isample]+".root");
    
    TH1D * HistOld = (TH1D*)fileMC->Get(VariableName+"_"+RegionName+"H");
    
    TH1D * Hist = TH1DtoTH1D(HistOld,nBinsNew,bins,true,"_new_"+mc_samples[isample]);
    
    TH1D * eventCount = (TH1D*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[isample]*lumi/nGen;
    
    
    TH1D * tempHist = EWKH;
    
    if (isample>3&&isample<9) tempHist = TTH;
    
    if (isample==9||isample==10) tempHist = DYH;
    
    if (isample==11) tempHist = QCDH;
    
    if (isample==12) tempHist = SUSYH;
    
    tempHist->Add(tempHist,Hist,1.,norm);

  }


  ////// **** Total Yields **** ////////
  float EWKTot = 0;
  float TTTot = 0;
  float QCDTot = 0;
  float DYTot = 0;
  float SUSYTot = 0;

  ////// **** Uncertainties **** ////////
  ///// Systematics /////
  float lumiSys = 0.03;
  float lepSys  = 0.04;
  
  float EWKSys = 0.20;
  float TTSys  = 0.15;
  float QCDSys = 0.30;
  float DYSys   = 0.05;
  
  ///// Statistical /////
  float EWKE2 = 0;
  float TTE2 = 0;
  float QCDE2 = 0;
  float DYE2 = 0;
  float SUSYE2 = 0;
  
  for (int iB=1; iB<=nBinsNew; ++iB) 
  {
    
    float ewkX = EWKH->GetBinContent(iB);
    float ewkE = EWKH->GetBinError(iB);
    
    float ttX  = TTH->GetBinContent(iB);
    float ttE  = TTH->GetBinError(iB);
    
    float qcdX = QCDH->GetBinContent(iB);
    float qcdE = QCDH->GetBinError(iB);
    
    float dyX  = DYH->GetBinContent(iB);
    float dyE  = DYH->GetBinError(iB);
    
    float susyX  = SUSYH->GetBinContent(iB);
    float susyE  = SUSYH->GetBinError(iB);
    
    EWKE2 += ewkE*ewkE;
    TTE2 += ttE*ttE;
    QCDE2 += qcdE*qcdE;
    DYE2 += dyE*dyE;
    SUSYE2 += susyE*susyE;

    EWKTot += ewkX;
    TTTot += ttX;
    QCDTot += qcdX;
    DYTot += dyX;
    SUSYTot += susyX;
    

    if (dyX<0) dyX = 0;
    
    float ewkErr  = ewkX*EWKSys;
    float ttErr   = ttX*TTSys;
    float qcdErr  = qcdX*QCDSys;
    float dyErr    = dyX*DYSys;
    
    ttX  += ewkX;
    dyX   += ttX;
    qcdX += dyX;

    float lumiErr = qcdX*lumiSys;
    float lepErr  = qcdX*lepSys;

    float totErr = TMath::Sqrt(lumiErr*lumiErr+
			                   lepErr*lepErr+
			                   ttErr*ttErr+
			                   ewkErr*ewkErr+
			                   qcdErr*qcdErr+
			                   ewkE*ewkE+
			                   dyE*dyE+
			                   ttE*ttE+
			                   qcdE*qcdE);
    
    EWKH->SetBinContent(iB,ewkX);
    TTH->SetBinContent(iB,ttX);
    DYH->SetBinContent(iB,dyX);
    QCDH->SetBinContent(iB,qcdX);
    QCDH->SetBinError(iB,totErr);
    SUSYH->SetBinError(iB,0);
  }

  float EWKE = TMath::Sqrt(EWKE2);
  float TTE = TMath::Sqrt(TTE2);
  float QCDE = TMath::Sqrt(QCDE2);
  float DYE = TMath::Sqrt(DYE2);
  float SUSYE = TMath::Sqrt(SUSYE2);
 
  
  std::cout << std::endl;
  std::cout << "QCD : " << QCDTot << "+/-" << QCDE << std::endl;
  std::cout << "EWK : " << EWKTot << "+/-" << EWKE << std::endl;
  std::cout << "TTJ : " << TTTot  << "+/-" << TTE << std::endl;
  std::cout << "DY  : " << DYTot  << "+/-" << DYE << std::endl;
  std::cout << std::endl;
  
  float BkgTot = QCDTot + EWKTot + TTTot + DYTot;
  std::cout << std::endl;
  std::cout << "BKG : " << BkgTot << std::endl;
  std::cout << "DAT : " << DataH->GetSumOfWeights() << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "ma(10) : " << SUSYTot << "+/-" << SUSYE << std::endl;
  std::cout << std::endl;
  
  std::cout << "Background Composition ( " << RegionName << " )" << std::endl;
  std::cout << "*************************************************************************************************"<< std::endl;
  std::cout << std::setprecision(3);
  std::cout << "QCD : " << 100*QCDTot/BkgTot << "+/-" 
  << 100*QCDTot/BkgTot*TMath::Sqrt(1/BkgTot+(QCDE/QCDTot)*(QCDE/QCDTot)) << " % | ";
  std::cout << "EWK : " << 100*EWKTot/BkgTot << "+/-" 
  << 100*EWKTot/BkgTot*TMath::Sqrt(1/BkgTot+(EWKE/EWKTot)*(EWKE/EWKTot)) << " % | ";
  std::cout << "TTJ : " << 100*TTTot/BkgTot  << "+/-" 
  << 100*TTTot/BkgTot*TMath::Sqrt(1/BkgTot+(TTE/TTTot)*(TTE/TTTot))      << " % | ";
  std::cout << "DY  : " << 100*DYTot/BkgTot  << "+/-" 
  << 100*DYTot/BkgTot*TMath::Sqrt(1/BkgTot+(DYE/DYTot)*(DYE/DYTot))      << " % | ";
  std::cout << std::endl;
  std::cout << "*************************************************************************************************"<< std::endl;
  std::cout << std::endl;
 
 
  /////// Ratio Histograms ////////
  TH1D * RatioH = (TH1D*)QCDH->Clone("RatioH");
  TH1D * UnitErrorH = (TH1D*)QCDH->Clone("UnitErrorH");
  
  for (int iB=1; iB<=nBinsNew; iB++) 
  {
    
    double y2 = QCDH->GetBinContent(iB);
    double y2_err = QCDH->GetBinError(iB);
    
    double y1     = DataH->GetBinContent(iB);
    double y1_err = DataH->GetBinError(iB);
    
    double r;
    double r_err;
    
    double u;
    double u_err;
    
    if(y2<0.001) 
    {
		r=1000.;
		r_err=0.;
		u = 1.;
		u_err=0.;
    }
	else
	{
		r = y1/y2;
        r_err = y1_err/y2;
        u = 1.;   
        u_err = y2_err/y2;
	}
    
    RatioH->SetBinContent(iB,r);
    RatioH->SetBinError(iB,r_err);
    
    UnitErrorH->SetBinContent(iB,u);
    UnitErrorH->SetBinError(iB,u_err);
    
    QCDH->SetBinError(iB,0.);
    EWKH->SetBinError(iB,0);
    TTH->SetBinError(iB,0);
    DYH->SetBinError(iB,0);
    SUSYH->SetBinError(iB,0);
 
  }
  
  
  ////// **** Plotting **** ////////
  TCanvas *c = new TCanvas("c","c",900,900) ;
  c->cd();
  
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  
  upper->cd();
  
  QCDH->Draw();
  DYH->Draw("sameh");
  TTH->Draw("sameh");
  EWKH->Draw("sameh");
  SUSYH->Draw("same");
  DataH->Draw("e1psame");
  
  
  QCDH->SetStats(0);
  QCDH->GetXaxis()->SetRangeUser(xLower,xUpper); 
  QCDH->GetXaxis()->SetLabelSize(0.);
  QCDH->GetXaxis()->SetLabelOffset(0.);
  QCDH->GetYaxis()->SetTitleOffset(1.2);
  QCDH->GetYaxis()->SetTitle(ytitle);
  QCDH->GetYaxis()->SetTitleSize(0.045);
  QCDH->GetYaxis()->SetLabelSize(0.04);
  QCDH->GetYaxis()->SetTickLength(0.055);
  QCDH->GetYaxis()->SetTickSize(0.013);
  QCDH->GetYaxis()->SetRangeUser(0,2*QCDH->GetMaximum());
  if (logY) QCDH->GetYaxis()->SetRangeUser(yLower,100000*QCDH->GetMaximum());
  QCDH->SetFillColorAlpha(kOrange-4,0.99);
  QCDH->SetLineColor(kBlack);
  
  DYH->SetStats(0);
  DYH->SetLineColor(kBlack);
  DYH->SetFillColorAlpha(kBlue,0.9);
  
  TTH->SetStats(0);
  TTH->SetLineColor(kBlack);
  TTH->SetFillColorAlpha(kRed,0.9);
  
  EWKH->SetStats(0);
  EWKH->SetLineColor(kBlack);
  EWKH->SetFillColorAlpha(kGreen+2,0.9);
  
  SUSYH->SetStats(0);
  SUSYH->SetLineColor(kTeal);
  SUSYH->SetLineWidth(2);
  SUSYH->SetLineStyle(2);
  
  DataH->SetStats(0);
  DataH->SetLineColor(kBlack);
  DataH->SetLineWidth(1);
  DataH->SetMarkerSize(1);
  DataH->SetMarkerStyle(kFullCircle);

  upper->Modified();
  upper->SetTicks();
  upper->SetLeftMargin(0.13);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);


  TLegend * leg = new TLegend(0.5,0.45,0.8,0.75);
  SetLegendStyle(leg);
  leg->SetTextSize(0.04);
  leg->AddEntry(DataH,"Data","lep");
  leg->AddEntry(QCDH,"QCD #mu-enriched","f");
  leg->AddEntry(DYH,"Drell-Yan","f");
  leg->AddEntry(TTH,"t#bar{t}+Single top","f");
  leg->AddEntry(EWKH,"Electroweak","f");
  leg->AddEntry(SUSYH,"ggH (4#tau) (x100), m_{a_{1}}=10 GeV","l");
  if (drawLeg) leg->Draw();
  
  
  writeExtraText = true;
  extraText   = "Preliminary";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);
    
  upper->RedrawAxis();
  upper->Update();
  
  
  lower->cd();
  
  TLine *line = new TLine(RatioH->GetBinLowEdge(RatioH->FindBin(xLower)),1.,RatioH->GetBinLowEdge(RatioH->FindBin(xUpper))+RatioH->GetBinWidth(RatioH->FindBin(xUpper)),1.);
  line->SetLineColorAlpha(kBlack,0.75);
  line->SetLineStyle(2);
  
  RatioH->Draw("e1p");
  UnitErrorH->Draw("e2same");
  RatioH->Draw("e1psame");

  line->Draw("same");
 
  RatioH->SetStats(0);
  RatioH->GetXaxis()->SetRangeUser(xLower,xUpper); 
  RatioH->GetXaxis()->SetLabelFont(42);
  RatioH->GetXaxis()->SetLabelOffset(0.03);
  RatioH->GetXaxis()->SetLabelSize(0.14);
  RatioH->GetXaxis()->SetTitleSize(0.13);
  RatioH->GetXaxis()->SetTitleOffset(1.35);
  RatioH->GetXaxis()->SetTitle(xTitle);
  RatioH->GetXaxis()->SetTickLength(0.025);
  RatioH->GetXaxis()->SetTickSize(0.08);
  RatioH->GetYaxis()->SetRangeUser(0.,1.99);
  RatioH->GetYaxis()->SetLabelOffset(0.008);
  RatioH->GetYaxis()->SetLabelSize(0.08);
  RatioH->GetYaxis()->SetTitleSize(0.1);
  RatioH->GetYaxis()->SetNdivisions(6);
  RatioH->GetYaxis()->SetTitleOffset(0.45);
  RatioH->GetYaxis()->SetTitle("Ratio    ");
  RatioH->GetYaxis()->SetTickLength(0.025);
  RatioH->GetYaxis()->SetTickSize(0.02);
  RatioH->SetLineColor(kAzure+2);
  RatioH->SetLineWidth(1);
  RatioH->SetLineStyle(1);
  RatioH->SetMarkerStyle(kFullCircle);
  RatioH->SetMarkerColor(kAzure+2);
  RatioH->SetMarkerSize(1);
  
  UnitErrorH->SetStats(0);
  UnitErrorH->SetFillColorAlpha(kGray,0.8);
  UnitErrorH->SetLineColor(kGray);
  UnitErrorH->SetLineWidth(2);
  UnitErrorH->SetLineStyle(1);
  UnitErrorH->SetFillStyle(1001);
  
  
  lower->Modified();
  lower->SetTicks();
  lower->SetLeftMargin(0.13);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.45);
  
  
  lower->RedrawAxis();
  lower->Update();

  
  c->Print(VariableName+"_"+RegionName+"H.pdf","Portrait pdf");

}
