
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Author: Alberto Bragagnolo alberto.bragagnolo.3@studenti.unipd.it
//
//    This code uses the root files produced by the TagAndProbe code and computes the efficiencies fitting the Z peak 
//		with the plotting tool. In this code the eta and pt bins are defined. Root files with scale factor (for IdIso) or
//		efficiencies (for triggers) will be produces. The plots with the various fits will be stored in dedicated directories.
//		This code is for muon.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "DesyTauAnalyses/NTupleMaker/test/TP_2016/FitPassAndFail.C"

void TP_eff_mu(TString fileName = "SingleMuon_Run2016B_TP", // RooT file with tag-&-probe histograms 
 	   TString what = "IdIso", //what do you want to evaluated
 	   float iso = 0.1, //isolation cut to be used
	   float norm = 1 // luminosity normalization factor (1 for data) 
	   ) 
{

	gErrorIgnoreLevel = kFatal;

  // output inizialization 
  TString lepton = "Muon";
	TString OutFileName = fileName + "_" + lepton + "_" + what + "_IsoLt" + Form("%.2f", iso) + "_eff_Spring16";
	TFile * outputFile = new TFile(OutFileName+".root","recreate");

	// Title of axis in plots  
	TString yTitle = "Efficiency";
	TString xTitle = lepton+"  p_{T}[GeV]";
	TString xtit; 
	xtit = "m_{#mu#mu}[GeV]"; 

	//names of final graphs - suffix
	bool isData=false;
	if (fileName.Contains("SingleMuon")) isData = true;
  TString SampleName("_MC");
  if (isData) SampleName = "_Data";

  //open input file
	TFile * file = new TFile(fileName+".root");
  file->cd();
  TTree *t = (TTree*)(file->Get("TagProbe"));


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

  float ptBinsTrig_def[17] = {10,
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
  TString which = what; 
  if (what == "IdIso") which = "";
  which = "";
  TString histBaseName; 

  //definition of the passing and failing criterias

  TCut cut_flag_idiso_pass, cut_flag_idiso_fail, cut_flag_hlt_pass, cut_flag_hlt_fail, cut_pt, cut_eta;

  if (what == "IdIso") {
  	cut_flag_idiso_pass = Form("id_probe == 1 && iso_probe < %f", iso);
  	cut_flag__idiso_fail = Form("id_probe == 0 || iso_probe >= %f", iso);
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
  	}

  //Definition of the output directory names and creation of it
	TString dir_name = "Muon_";
	dir_name += what;
	dir_name += Form("%.2f", iso);
	if (!isData) dir_name += "_MC";
	dir_name += "_eff";
	gSystem->mkdir(dir_name, kTRUE);

	for (int iEta = 0; iEta < nEtaBins; iEta++) { //loop on eta bins

		histBaseName = prefix+which+EtaBins[iEta];

		//eta cuts
		cut_eta = Form("abs(eta_probe)>= %f && (eta_probe)< %f", etaBins[iEta], etaBins[iEta+1]);

		TH1F * numeratorH   = new TH1F("numeratorH","",nPtBins,ptBins_edges);
		TH1F * denominatorH = new TH1F("denominatorH","",nPtBins,ptBins_edges);
	
	  for (int iPt=0; iPt<nPtBins; ++iPt) { //loop on pt bins

	  	//pt cuts
	  	cut_pt = Form("pt_probe > %f && pt_probe < %f", ptBins[iPt], ptBins[iPt+1]);

	  	TH1F * histPassOld = new TH1F("histPassOld","",250,50,300);
	  	TH1F * histFailOld = new TH1F("histFailOld","",250,50,300);
	  	
	  	//Drawing histogram of passing and failing probes
	  	if (what == "IdIso") {
		  	t->Draw("m_vis>>histPassOld", "pu_weight"*"mcweight"*(cut_eta && cut_pt && cut_flag_idiso_pass));
		  	t->Draw("m_vis>>histFailOld", "pu_weight"*"mcweight"*(cut_eta && cut_pt && cut_flag_idiso_fail));
		  }else{
		  	t->Draw("m_vis>>histPassOld", "pu_weight"*"mcweight"*(cut_eta && cut_pt && cut_flag_hlt_pass && cut_flag_idiso_pass));
		  	t->Draw("m_vis>>histFailOld", "pu_weight"*"mcweight"*(cut_eta && cut_pt && cut_flag_hlt_fail && cut_flag_idiso_pass));
		  }

	  	int nBinsX = histPassOld->GetNbinsX();

	  	//lumi renormalization
	    for (int iB=1;iB<=nBinsX;++iB) {
	      histPassOld->SetBinContent(iB,norm*histPassOld->GetBinContent(iB));
	      histPassOld->SetBinError(iB,norm*histPassOld->GetBinError(iB));
	      histFailOld->SetBinContent(iB,norm*histFailOld->GetBinContent(iB));
	      histFailOld->SetBinError(iB,norm*histFailOld->GetBinError(iB));
	    }

	    float output[2];
	    TCanvas * c1 = new TCanvas("c1","",700,600);
	    TCanvas * c2 = new TCanvas("c2","",700,600);


	    //defining fit options
	    bool fitPass = true; 
	    bool fitFail = true;
	    bool rebinPass = false;
	    bool rebinFail = false;
	    if(what != "IdIso") {fitPass= false; fitFail = false;}

	    //fitting
	    FitPassAndFail(fileName,
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
	    	dir_name);


	    c1->cd();
	    c1->Update();
	    c2->cd();
	    c2->Update();
	    numeratorH->SetBinContent(iPt+1,output[0]);
	    denominatorH->SetBinContent(iPt+1,output[0]+output[1]);

	  }

	  outputFile->cd();

	  //produce efficiencies plot
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
	  if(!isData){
	  	eff->SetMarkerColor(kRed);
	  	eff->SetLineColor(kRed);
	  }


	  TCanvas * canv = new TCanvas("canv","",700,600);
	  eff->Draw("APE");
	  canv->SetGridx();
	  canv->SetGridy();
	  canv->Update();

	  canv->SaveAs(dir_name + "/" + fileName+"_" + histBaseName + ".png");
	  eff->Write(histBaseName+SampleName);

/*	  for(int ip=0; ip<nPtBins; ++ip){
		cout<<"PtBins "<<ip<<" content: "<<numeratorH->GetBinContent(ip)/ denominatorH->GetBinContent(ip)<<endl;
	}*/

	}

	//closing
	outputFile->cd(); 
	etaBinsH->Write();
	outputFile->Close();

}