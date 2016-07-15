
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//
//    This code uses the root files produced by the TagAndProbe code, both Data and MC at the same time, and computes
//		the efficiencies fitting the Z peak with the fitting tool.
//		In this code the eta and pt bins are defined. Root files with scale factor and efficiencies will be produced. 
//		The plots with the various fits will be stored in dedicated directories.
//		This code is for muon.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "DesyTauAnalyses/NTupleMaker/test/TP_2016/FitPassAndFail.C"

void SF_mu(TString fileName_Data = "SingleMuon__Run2016-PromptReco-v2.root", //data otree
					TString fileName_MC ="DYJetsToLL_TP_Muon.root", //MC otree
					TString what = "IdIso", // or Mu8, IsoMu, ...
					float iso = 0.15, //isolation cut to be used
					float norm = 1 // luminosity normalization factor (1 for data) 
					) 
{

	gErrorIgnoreLevel = kFatal;

  // output inizialization 
  TString lepton = "Muon";
	TString OutFileName = lepton + "_" + what + "_IsoLt" + Form("%.2f", iso) + "_eff_Spring16";
	TFile * outputFile = new TFile(OutFileName+".root","recreate");

	//open input file
	TFile * file1 = new TFile(fileName_Data);
	TFile * file2 = new TFile(fileName_MC);
  file1->cd();
  TTree *t1 = (TTree*)(file1->Get("TagProbe"));
  file2->cd();
  TTree *t2 = (TTree*)(file2->Get("TagProbe"));
  file1->cd();

	//names of final graphs - suffix
	TString SampleName = "_Data";

	// Title of axis in plots  
	TString yTitle = "Efficiency";
	TString xTitle = lepton+"  p_{T}[GeV]";
	TString xtit; 
	xtit = "m_{#mu#mu}[GeV]"; 

  //binning inizialization

  int nEtaBins = 3;
	float etaBins[4] = {0,0.9,1.2,2.4};

	TString EtaBins[3] = {"EtaLt0p9",
				"Eta0p9to1p2",
				"EtaGt1p2"};

	float ptBins_def[8] = {10,15,20,25,30,40,60,1000};

	TString PtBins_def[7] = {"Pt10to15",
       "Pt15to20",
       "Pt20to25",
       "Pt25to30",
       "Pt30to40",
       "Pt40to60",
       "PtGt60"};

  float ptBinsTrig_def[17] = {
  			10,
			  13,
			  16,
			  19,
			  22,
			  25,
			  28,
			  31,
			  34,
			  37,
			  40,
			  45,
			  50,
			  60,
			  70,
			  100,
				1000};

	TString PtBinsTrig_def[16] = {"Pt10to13",
		    "Pt13to16",
		    "Pt16to19",
		    "Pt19to22",
		    "Pt22to25",
		    "Pt25to28",
		    "Pt28to31",
		    "Pt31to34",
		    "Pt34to37",
		    "Pt37to40",
		    "Pt40to45",
		    "Pt45to50",
		    "Pt50to60",
		    "Pt60to70",
		    "Pt70to100",
		    "PtGt100"};

	int nPtBins = 16; if(what == "IdIso") nPtBins = 7;
	float * ptBins = new float[nPtBins+1];
	TString * PtBins = new TString[nPtBins];

	if(what == "IdIso"){
		for(int i=0; i<nPtBins; ++i){
			ptBins[i] = ptBins_def[i];
			PtBins[i] = PtBins_def [i];
		}
		ptBins[nPtBins] = ptBins_def[nPtBins];
	} else {
		for(int i=0; i<nPtBins; ++i){
			ptBins[i] = ptBinsTrig_def[i];
			PtBins[i] = PtBinsTrig_def[i];
		}
		ptBins[nPtBins] = ptBinsTrig_def[nPtBins];
	}


  //	create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //

	TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
  etaBinsH->Draw();
  etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
  for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}
  etaBinsH->Draw();


	//	create pt histogram_s with pt ranges associated to their names (eg. Pt10to13, ..)   ***** //

	TH1D * ptBinsH =  new TH1D("ptBinsH", "ptBinsH", nPtBins, ptBins);
  ptBinsH->Draw();
  ptBinsH->GetXaxis()->Set(nPtBins, ptBins);
  for (int i=0; i<nPtBins; i++){ ptBinsH->GetXaxis()->SetBinLabel(i+1, PtBins[i]);}
  ptBinsH->Draw();

  float ptBins_edges[nPtBins+1];

	for (int i=0; i<nPtBins; i++) { ptBins_edges[i]=ptBinsH->GetBinLowEdge(i+1); }
	ptBins_edges[nPtBins]= ptBinsH->GetBinLowEdge(nPtBins+1); 

  
	// define if in the fit of failing probes
  // the FSR component will be used in the 
  // signal function 
  bool fitWithFSR[nPtBins];

  for (int i=0; i<nPtBins; i++)  fitWithFSR[i] = true;

  if(what == "IdIso"){
   	fitWithFSR[0]=false;
    fitWithFSR[1]=false;
    fitWithFSR[5]=false;
    fitWithFSR[6]=false;
  } else{ for (int i=0; i<nPtBins; i++)  fitWithFSR[i] = false; }

	// building the histogram base name
  TString prefix = "ZMass";
  TString which = what; if (what == "IdIso") which = "";
  TString histBaseName; 

  TCut cut_flag_idiso_pass, cut_flag_hlt_pass, cut_flag_hlt_fail, cut_pt, cut_eta;

  if (what == "IdIso") {
  	cut_flag_idiso_pass = Form("id_probe == 1 && iso_probe < %f", iso);
  } else{
  		if(what == "hlt_1") {cut_flag_hlt_pass = "hlt_1_probe == 1"; cut_flag_hlt_fail = "hlt_1_probe == 0"; }
  		if(what == "hlt_2") {cut_flag_hlt_pass = "hlt_2_probe == 1"; cut_flag_hlt_fail = "hlt_2_probe == 0"; }
  		if(what == "hlt_3") {cut_flag_hlt_pass = "hlt_3_probe == 1"; cut_flag_hlt_fail = "hlt_3_probe == 0"; }
  		if(what == "hlt_4") {cut_flag_hlt_pass = "hlt_4_probe == 1"; cut_flag_hlt_fail = "hlt_4_probe == 0"; }
  		if(what == "hlt_5") {cut_flag_hlt_pass = "hlt_5_probe == 1"; cut_flag_hlt_fail = "hlt_5_probe == 0"; }
  		if(what == "hlt_6") {cut_flag_hlt_pass = "hlt_6_probe == 1"; cut_flag_hlt_fail = "hlt_6_probe == 0"; }
  		if(what == "hlt_7") {cut_flag_hlt_pass = "hlt_7_probe == 1"; cut_flag_hlt_fail = "hlt_7_probe == 0"; }
  		if(what == "hlt_8") {cut_flag_hlt_pass = "hlt_8_probe == 1"; cut_flag_hlt_fail = "hlt_8_probe == 0"; }
  		if(what == "hlt_9") {cut_flag_hlt_pass = "hlt_9_probe == 1"; cut_flag_hlt_fail = "hlt_9_probe == 0"; }
  		if(what == "hlt_10") {cut_flag_hlt_pass = "hlt_10_probe == 1"; cut_flag_hlt_fail = "hlt_10_probe == 0"; }
  		if(what == "hlt_11") {cut_flag_hlt_pass = "hlt_11_probe == 1"; cut_flag_hlt_fail = "hlt_11_probe == 0"; }
  		if(what == "hlt_12") {cut_flag_hlt_pass = "hlt_12_probe == 1"; cut_flag_hlt_fail = "hlt_12_probe == 0"; }
  		if(what == "hlt_13") {cut_flag_hlt_pass = "hlt_13_probe == 1"; cut_flag_hlt_fail = "hlt_13_probe == 0"; }
  		if(what == "hlt_14") {cut_flag_hlt_pass = "hlt_14_probe == 1"; cut_flag_hlt_fail = "hlt_14_probe == 0"; }
  		if(what == "hlt_15") {cut_flag_hlt_pass = "hlt_15_probe == 1"; cut_flag_hlt_fail = "hlt_15_probe == 0"; }
  		if(what == "hlt_16") {cut_flag_hlt_pass = "hlt_16_probe == 1"; cut_flag_hlt_fail = "hlt_16_probe == 0"; }
  		if(what == "hlt_17") {cut_flag_hlt_pass = "hlt_17_probe == 1"; cut_flag_hlt_fail = "hlt_17_probe == 0"; }
  		if(what == "hlt_18") {cut_flag_hlt_pass = "hlt_18_probe == 1"; cut_flag_hlt_fail = "hlt_18_probe == 0"; }
  		if(what == "hlt_19") {cut_flag_hlt_pass = "hlt_19_probe == 1"; cut_flag_hlt_fail = "hlt_19_probe == 0"; }
  		if(what == "hlt_20") {cut_flag_hlt_pass = "hlt_20_probe == 1"; cut_flag_hlt_fail = "hlt_20_probe == 0"; }
  	}



	TString dir_name1 = "Muon_";
	TString dir_name2 = "Muon_";
	dir_name1 += what;
	dir_name2 += what;
	dir_name1 += Form("%.2f", iso);
	dir_name2 += Form("%.2f", iso);
	dir_name2 += "_MC";
	dir_name1 += "_eff";
	dir_name2 += "_eff";
	gSystem->mkdir(dir_name1, kTRUE);
	gSystem->mkdir(dir_name2, kTRUE);

////////////////////DATA

	for (int iEta = 0; iEta < nEtaBins; iEta++) {

		histBaseName = prefix+which+EtaBins[iEta];

		cut_eta = Form("abs(eta_probe)>= %f && abs(eta_probe)< %f", etaBins[iEta], etaBins[iEta+1]);

		TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins_edges);
		TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins_edges);
	
	  for (int iPt=0; iPt<nPtBins; ++iPt) {

	  	cut_pt = Form("pt_probe > %f && pt_probe < %f", ptBins[iPt], ptBins[iPt+1]);

	  	TH1F * histPassOld = new TH1F("histPassOld","",250,50,300);
	  	TH1F * histFailOld = new TH1F("histFailOld","",250,50,300);
	  	
	  	//Drawing histogram of passing and failing probes
	  	if (what == "IdIso") {
		  	t1->Draw("m_vis>>histPassOld", (cut_eta && cut_pt && cut_flag_idiso_pass));
		  	t1->Draw("m_vis>>histFailOld", (cut_eta && cut_pt && !cut_flag_idiso_pass));
		  }else{
		  	t1->Draw("m_vis>>histPassOld", (cut_eta && cut_pt && cut_flag_hlt_pass && cut_flag_idiso_pass));
		  	t1->Draw("m_vis>>histFailOld", (cut_eta && cut_pt && cut_flag_hlt_fail && cut_flag_idiso_pass));
		  }

	  	int nBinsX = histPassOld->GetNbinsX();

	    for (int iB=1;iB<=nBinsX;++iB) {
	      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
	      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
	      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
	      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
	    }

	    float output[2];
	    TCanvas * c1 = new TCanvas("c1","",700,600);
	    TCanvas * c2 = new TCanvas("c2","",700,600);
	    bool fitPass = true; 
	    bool fitFail = true;
	    bool rebinPass = false;
	    bool rebinFail = false;
	    if(what != "IdIso") {fitPass= false; fitFail = false;}


	    FitPassAndFail(fileName_Data,
	    	histBaseName+PtBins[iPt],
	    	xtit,
	    	histPassOld,
	    	histFailOld,
	    	fitPass,
	    	fitFail,
	    	fitWithFSR[iPt],
	    	rebinPass,
	    	rebinFail,
	    	c1,
	    	c2,
	    	output,
	    	dir_name1);


	    c1->cd();
	    c1->Update();
	    c2->cd();
	    c2->Update();
	    numeratorH->SetBinContent(iPt+1,output[0]);
	    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);

	  }

	  outputFile->cd();

	  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
	  eff->Divide(numeratorH,denominatorH);

	  eff->GetXaxis()->SetTitle(xTitle);
	  //  eff->GetXaxis()->SetRangeUser(10.01,59.99);
	  eff->GetYaxis()->SetRangeUser(0,1.0);
	  eff->GetXaxis()->SetRangeUser(0,99.99);
	  eff->GetYaxis()->SetTitle(yTitle);
	  eff->GetXaxis()->SetTitleOffset(1.1);
	  eff->GetXaxis()->SetNdivisions(510);
	  eff->GetYaxis()->SetTitleOffset(1.1);
	  eff->SetMarkerStyle(21);
	  eff->SetMarkerSize(1);
	  eff->SetMarkerColor(kBlue);
	  eff->SetLineWidth(2);
	  eff->SetLineColor(kBlue);


	  TCanvas * canv = new TCanvas("canv","",700,600);
	  eff->Draw("APE");
	  canv->SetGridx();
	  canv->SetGridy();
	  canv->Update();

	  canv->SaveAs(dir_name1 + "/" + fileName_Data+"_" + histBaseName + ".png");
	  eff->Write(histBaseName+SampleName);

	}

//////////////////////////////////////////////MC

	SampleName = "_MC";

	for (int iEta = 0; iEta < nEtaBins; iEta++) {

		histBaseName = prefix+which+EtaBins[iEta];

		cut_eta = Form("abs(eta_probe)>= %f && abs(eta_probe)< %f", etaBins[iEta], etaBins[iEta+1]);

		TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins_edges);
		TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins_edges);
	
	  for (int iPt=0; iPt<nPtBins; ++iPt) {

	  	cut_pt = Form("pt_probe > %f && pt_probe < %f", ptBins[iPt], ptBins[iPt+1]);

	  	TH1F * histPassOld = new TH1F("histPassOld","",250,50,300);
	  	TH1F * histFailOld = new TH1F("histFailOld","",250,50,300);

	  	
	  	//Drawing histogram of passing and failing probes
	  	if (what == "IdIso") {
		  	t1->Draw("m_vis>>histPassOld", "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_idiso_pass));
		  	t1->Draw("m_vis>>histFailOld", "pu_weight*mcweight" + (cut_eta && cut_pt && !cut_flag_idiso_pass));
		  }else{
		  	t1->Draw("m_vis>>histPassOld", "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_hlt_pass && cut_flag_idiso_pass));
		  	t1->Draw("m_vis>>histFailOld", "pu_weight*mcweight" + (cut_eta && cut_pt && cut_flag_hlt_fail && cut_flag_idiso_pass));
		  }


	  	int nBinsX = histPassOld->GetNbinsX();

	    for (int iB=1;iB<=nBinsX;++iB) {
	      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
	      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
	      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
	      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
	    }

	    float output[2];
	    TCanvas * c1 = new TCanvas("c1","",700,600);
	    TCanvas * c2 = new TCanvas("c2","",700,600);
	    bool fitPass = true; 
	    bool fitFail = true;
	    bool rebinPass = false;
	    bool rebinFail = false;
	    if(what != "IdIso") {fitPass= false; fitFail = false;}


	    FitPassAndFail(fileName_MC,
	    	histBaseName+PtBins[iPt],
	    	xtit,
	    	histPassOld,
	    	histFailOld,
	    	fitPass,
	    	fitFail,
	    	fitWithFSR[iPt],
	    	rebinPass,
	    	rebinFail,
	    	c1,
	    	c2,
	    	output,
	    	dir_name2);


	    c1->cd();
	    c1->Update();
	    c2->cd();
	    c2->Update();
	    numeratorH->SetBinContent(iPt+1,output[0]);
	    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);

	  }

	  outputFile->cd();

	  TGraphAsymmErrors * eff = new TGraphAsymmErrors();
	  eff->Divide(numeratorH,denominatorH);

	  eff->GetXaxis()->SetTitle(xTitle);
	  //  eff->GetXaxis()->SetRangeUser(10.01,59.99);
	  eff->GetYaxis()->SetRangeUser(0,1.0);
	  eff->GetXaxis()->SetRangeUser(0,99.99);
	  eff->GetYaxis()->SetTitle(yTitle);
	  eff->GetXaxis()->SetTitleOffset(1.1);
	  eff->GetXaxis()->SetNdivisions(510);
	  eff->GetYaxis()->SetTitleOffset(1.1);
	  eff->SetMarkerStyle(21);
	  eff->SetMarkerSize(1);
	  eff->SetLineWidth(2);
  	eff->SetMarkerColor(kRed);
  	eff->SetLineColor(kRed);



	  TCanvas * canv = new TCanvas("canv","",700,600);
	  eff->Draw("APE");
	  canv->SetGridx();
	  canv->SetGridy();
	  canv->Update();

	  canv->SaveAs(dir_name2 + "/" + fileName_MC+"_" + histBaseName + ".png");
	  eff->Write(histBaseName+SampleName);

	}

	outputFile->cd(); 
  etaBinsH->Write();
  outputFile->Close();

}