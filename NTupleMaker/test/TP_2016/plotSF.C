#include <TH1.h>
 
void plot_eff(TString fileName = "Muon_IdIso_eff_Spring16.root"
    ) {
    TFile * file = new TFile(fileName);
    file->cd();
    TString lep; TString what;
    bool IdIso=false;
     
    if (fileName.Contains("Electron")) lep = "Electron";
    if (fileName.Contains("Muon")) lep = "Muon";
    if (fileName.Contains("IdIso")) {IdIso=true; what="IdIso";}
    TString lep_latex;
    if(lep == "Muon") lep_latex = "#mu";
    if(lep == "Electron") lep_latex = "e";
 
    int nEtaBins=3;
    if(lep == "Electron") nEtaBins=2;
 
 
    TString *names = new TString[nEtaBins];
 
    if(lep == "Muon"){
        names[0]="EtaLt0p9";
        names[1]="Eta0p9to1p2";
        names[2]="EtaGt1p2";
    }
 
    if(lep == "Electron"){
        names[0]="EtaLt1p48";
        names[1]="EtaGt1p48";
    }
 
    for(int i=0; i<nEtaBins; ++i){
 
        TGraphAsymmErrors *gr_data = (TGraphAsymmErrors*)gDirectory->Get("ZMass"+names[i]+"_Data");
        TGraphAsymmErrors *gr_MC = (TGraphAsymmErrors*)gDirectory->Get("ZMass"+names[i]+"_MC");
 
        int N = gr_data->GetN();
        Double_t x; Double_t y1; Double_t y2;
        gr_data->GetPoint(N-1, x, y1);
        gr_data->SetPoint(N-1, 70, y1);
        gr_data->SetPointEXhigh(N-1, 10);
        gr_data->SetPointEXlow(N-1, 10);
 
        gr_MC->GetPoint(N-1, x, y1);
        gr_MC->SetPoint(N-1, 70, y1);
        gr_MC->SetPointEXhigh(N-1, 10);
        gr_MC->SetPointEXlow(N-1, 10);
 
        TCanvas * c1 = new TCanvas("c1", "", 700, 800);
 
        TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
        upper->Draw();
        upper->cd();
 
        upper->SetFillColor(0);
    upper->SetBorderMode(0);
    upper->SetBorderSize(10);
    upper->SetTickx(1);
    upper->SetTicky(1); 
 
    upper->SetLeftMargin(0.15);
    upper->SetRightMargin(0.05);
    upper->SetBottomMargin(0.02);
    upper->SetTopMargin(0.05);
 
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
        gr_data->GetXaxis()->SetTitleSize(0);
        gr_data->GetXaxis()->SetLabelSize(0);
//      gr_data->GetXaxis()->SetRangeUser(0,59.99);
 
        gr_data->GetYaxis()->SetTitle("Efficiency");
        gr_data->GetYaxis()->SetTitleOffset(1);
    gr_data->GetYaxis()->SetTitleSize(0.06);
    gr_data->GetYaxis()->SetLabelSize(0.05);
        gr_data->GetYaxis()->SetRangeUser(0,1);
        gr_data->GetYaxis()->SetDecimals();
 
        gr_data->Draw("APE");
        upper->Update();
        upper->SetGridx();
        upper->SetGridy();
         
        gr_MC->Draw("PE");
 
        //LEGEND
 
        TLegend * leg = new TLegend(0.6,0.1,0.92,0.3);
        TString eta_string;
        leg->SetTextSize(0.05);
        leg->AddEntry(gr_data, "Data", "P");
        leg->AddEntry(gr_MC, "MC", "P");
        if(names[i].Contains("Lt0p9"))      eta_string = "|#eta|<0.9";
        if(names[i].Contains("0p9to1p2")) eta_string = "0.9<|#eta|<1.2";
        if(names[i].Contains("Gt1p2"))      eta_string = "|#eta|>1.2";
        if(names[i].Contains("Lt1p48"))     eta_string = "|#eta|<1.48";
        if(names[i].Contains("Gt1p48"))     eta_string = "|#eta|>1.48";
        leg->SetHeader(lep_latex+" "+what+ " "+ eta_string);
 
 
        leg->Draw();
 
        upper->Draw("SAME");
    upper->RedrawAxis();
    upper->Modified();
    upper->Update();
    c1->cd();
 
        //RATIO
 
        TGraphAsymmErrors * ratio = new TGraphAsymmErrors(N);
 
         
        for(int in=0; in<N; ++in){
            gr_data->GetPoint(in, x, y1);
            gr_MC->GetPoint(in, x, y2);
            if(in != N-1){
                ratio->SetPoint(in, x, y1/y2);
                ratio->SetPointEXhigh(in, gr_data->GetErrorXhigh(in));
                ratio->SetPointEXlow(in, gr_data->GetErrorXlow(in));  
            }else{
                ratio->SetPoint(in, 70, y1/y2);
                ratio->SetPointEXhigh(in, 10);
                ratio->SetPointEXlow(in, 10);
            }       
        }
 
        ratio->SetTitle("");         
        ratio->GetYaxis()->SetTitle("Scale Factor");
 
        ratio->GetYaxis()->SetLabelFont(42);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetDecimals();
        ratio->GetYaxis()->SetLabelSize(0.08);        
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.6);       
    ratio->GetYaxis()->SetRangeUser(0.8,1.2);
        ratio->GetYaxis()->SetLabelSize(0.1);
    if(lep == "Muon") ratio->GetYaxis()->SetRangeUser(0.9,1.1);
 
    ratio->GetXaxis()->SetTitle(gr_data->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetLabelFont(42);
    ratio->GetXaxis()->SetLabelOffset(0.04);
    ratio->GetXaxis()->SetLabelSize(0.14);
    ratio->GetXaxis()->SetTitleSize(0.13);
    ratio->GetXaxis()->SetTitleOffset(1.2);
    ratio->GetXaxis()->SetLabelFont(42);
//      ratio->GetXaxis()->SetRangeUser(0,59.99);
        ratio->GetXaxis()->SetTickLength(0.07);
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
 
        lower->SetLeftMargin(0.15);
        lower->SetRightMargin(0.05);
        lower->SetTopMargin(0.05);
        lower->SetBottomMargin(0.35);
 
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
 
//      TPaveText *title = (TPaveText*)lower->GetPrimitive("title");
//      title->SetTextSize(0.09);
 
    lower->Modified();
    lower->RedrawAxis();
 
    c1->cd();
    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
        c1->SaveAs(lep + "_" + "ZMass" + names[i] + "_eff.png");
    }
 
 
}