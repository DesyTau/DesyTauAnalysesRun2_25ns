#include <TH1.h>

void plot_eff(TString fileName1 = "SingleMuon__Run2016-PromptReco-v2_Muon_IdIso_IsoLt0.15_eff_Spring16.root",
							TString fileName2 = "DYJetsToLL_TP_Muon_Muon_IdIso_IsoLt0.15_eff_Spring16.root"
	) {
	TFile * file1 = new TFile(fileName1);
	TFile * file2 = new TFile(fileName2);

	TString lep;
	bool IdIso=false;
	
	if (fileName1.Contains("SingleElectron")) lep = "Electron";
	if (fileName1.Contains("SingleMuon")) lep = "Muon";
	if (fileName1.Contains("IdIso")) IdIso=true;

	int nEtaBins=3;
	if(lep == "Electron") nEtaBins=2;


	TString *names = new TString[nEtaBins];

	if(lep == "Muon"){
		names[0]="ZMassEtaLt0p9";
		names[1]="ZMassEta0p9to1p2";
		names[2]="ZMassEtaGt1p2";
	}

	if(lep == "Electron"){
		names[0]="ZMassEtaLt1p48";
		names[1]="ZMassEta1p48to2p5";
	}

	for(int i=0; i<nEtaBins; ++i){
		file1->cd();
		TGraphAsymmErrors *gr_data = (TGraphAsymmErrors*)gDirectory->Get(names[i]+"_Data");

		TCanvas * c1 = new TCanvas("c1", "", 700, 800);

		TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
		upper->Draw();
		upper->cd();

		upper->SetFillColor(0);
    upper->SetBorderMode(0);
    upper->SetBorderSize(10);
    upper->SetTickx(1);
    upper->SetTicky(1); 
    upper->SetFrameFillStyle(0);
    upper->SetFrameLineStyle(0);
    upper->SetFrameLineWidth(2);
    upper->SetFrameBorderMode(0);
    upper->SetFrameBorderSize(10);
    upper->SetFrameFillStyle(0);
    upper->SetFrameLineStyle(0);
    upper->SetFrameLineWidth(2);
    upper->SetFrameBorderMode(0);
    upper->SetFrameBorderSize(10);


		gr_data->Draw("APE");
		gr_data->GetXaxis()->SetTitleOffset(0.85);
		gr_data->GetYaxis()->SetTitleOffset(0.85);
		gr_data->GetXaxis()->SetRangeUser(0,59.99);
		upper->SetGridx();
		upper->SetGridy();
		
		file2->cd();
		TGraphAsymmErrors *gr_MC = (TGraphAsymmErrors*)gDirectory->Get(names[i]+"_MC");
		gr_MC->Draw("PE");

		upper->Draw("SAME");
    upper->RedrawAxis();
    upper->Modified();
    upper->Update();
    c1->cd();

		//RATIO
		int N = gr_data->GetN();
		TGraphAsymmErrors * ratio = new TGraphAsymmErrors(N);

		for(int in=0; in<N; ++in){
			Double_t x; Double_t y1; Double_t y2;
			gr_data->GetPoint(in, x, y1);
			gr_MC->GetPoint(in, x, y2);
			ratio->SetPoint(in, x, y1/y2);
			ratio->SetPointEXhigh(in, gr_data->GetErrorXhigh(in));
			ratio->SetPointEXlow(in, gr_data->GetErrorXlow(in));
			cout<<"in: "<<in<<", x: "<<x<<", y: "<<y1/y2<<endl;
		}

		ratio->SetTitle("Scale Factor"); 		
   	ratio->GetYaxis()->SetRangeUser(0.8,1.2);
   	if(lep == "Muon") ratio->GetYaxis()->SetRangeUser(0.9,1.1);
		ratio->GetYaxis()->SetNdivisions(505);
		ratio->GetYaxis()->SetDecimals();
		ratio->GetYaxis()->SetLabelSize(0.1);
		ratio->GetXaxis()->SetLabelSize(0);
		ratio->GetXaxis()->SetRangeUser(0,59.99);
		ratio->GetXaxis()->SetLabelSize(0);
		ratio->SetMarkerStyle(21);
	  ratio->SetMarkerSize(1);

  	TPad *lower = new TPad("lower", "pad",0,0,1,0.30);
  	lower->Draw();
  	lower->cd();
  	lower->SetFillColor(0);
  	lower->SetBorderMode(0);
  	lower->SetBorderSize(10);
  	lower->SetGridy();
  	lower->SetGridx();
  	lower->SetTickx(1);
  	lower->SetTicky(1);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);
  	ratio->Draw("APE");
  	lower->Update();

  	TLine *l=new TLine(lower->GetUxmin(),1.0,lower->GetUxmax(),1.0);
  	l->SetLineWidth(3);
  	l->Draw();

  	TPaveText *title = (TPaveText*)lower->GetPrimitive("title");
   	title->SetTextSize(0.09);

  	lower->Modified();
  	lower->RedrawAxis();

  	c1->cd();
  	c1->Modified();
  	c1->cd();
  	c1->SetSelected(c1);

		c1->SaveAs(lep + "_" + names[i] + "_eff.png");
	}


}