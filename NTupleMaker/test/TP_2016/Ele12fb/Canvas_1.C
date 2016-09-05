void Canvas_1()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Thu Jul 28 21:29:31 2016) by ROOT version6.06/01
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,74,970,577);
   Canvas_1->Range(-27.4975,-0.125,247.4775,1.125);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   Double_t Graph_from_ratioH_fx3001[17] = {
   11.5,
   14.5,
   17.5,
   20.5,
   23.5,
   26.5,
   29.5,
   32.5,
   35.5,
   38.5,
   42.5,
   47.5,
   55,
   65,
   85,
   150,
   5100};
   Double_t Graph_from_ratioH_fy3001[17] = {
   1,
   0,
   0,
   0,
   0.01771524,
   0.2943746,
   0.4625348,
   0.5249608,
   0.5746665,
   0.613711,
   0.6641073,
   0.7019153,
   0.71673,
   0.7435814,
   0.7597724,
   0.7770616,
   0.7787418};
   Double_t Graph_from_ratioH_felx3001[17] = {
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   2.5,
   2.5,
   5,
   5,
   15,
   50,
   4900};
   Double_t Graph_from_ratioH_fely3001[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph_from_ratioH_fehx3001[17] = {
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   1.5,
   2.5,
   2.5,
   5,
   5,
   15,
   50,
   4900};
   Double_t Graph_from_ratioH_fehy3001[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(17,Graph_from_ratioH_fx3001,Graph_from_ratioH_fy3001,Graph_from_ratioH_felx3001,Graph_from_ratioH_fehx3001,Graph_from_ratioH_fely3001,Graph_from_ratioH_fehy3001);
   grae->SetName("Graph_from_ratioH");
   grae->SetTitle("");
   grae->SetFillColor(2);
   grae->SetFillStyle(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineStyle(0);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   
   TH1F *Graph_Graph_from_ratioH3001 = new TH1F("Graph_Graph_from_ratioH3001","",100,0,10999);
   Graph_Graph_from_ratioH3001->SetMinimum(0);
   Graph_Graph_from_ratioH3001->SetMaximum(1);
   Graph_Graph_from_ratioH3001->SetDirectory(0);
   Graph_Graph_from_ratioH3001->SetStats(0);
   Graph_Graph_from_ratioH3001->SetFillColor(2);
   Graph_Graph_from_ratioH3001->SetFillStyle(0);
   Graph_Graph_from_ratioH3001->SetLineStyle(0);
   Graph_Graph_from_ratioH3001->SetLineWidth(3);
   Graph_Graph_from_ratioH3001->SetMarkerStyle(20);
   Graph_Graph_from_ratioH3001->SetMarkerSize(1.4);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetTitle("Electron  p_{T}[GeV]");
   Graph_Graph_from_ratioH3001->GetXaxis()->SetRange(1,2);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetLabelFont(42);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetTitleSize(0.055);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph_from_ratioH3001->GetXaxis()->SetTitleFont(42);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetTitle("Efficiency");
   Graph_Graph_from_ratioH3001->GetYaxis()->SetLabelFont(42);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetTitleSize(0.055);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph_from_ratioH3001->GetYaxis()->SetTitleFont(42);
   Graph_Graph_from_ratioH3001->GetZaxis()->SetLabelFont(42);
   Graph_Graph_from_ratioH3001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_from_ratioH3001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_from_ratioH3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_from_ratioH3001);
   
   grae->Draw("alp");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
