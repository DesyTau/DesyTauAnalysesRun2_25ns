#include "DesyTauAnalyses/NTupleMaker/interface/CalibrationOfImpactParameters.h"


CalibrationOfImpactParameters::CalibrationOfImpactParameters()
{
   // TString NameFile_eta = "/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/ImpactParameterCalibration/"+IP_Name+"_binned_eta.root";
   // f_eta_np = new TFile(NameFile_eta, "READ");
   // TString NameFile_pt = "/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/ImpactParameterCalibration/"+IP_Name+"_binned_pT.root";
   // f_pt_np = new TFile(NameFile_pt, "READ");

   
   TString NameFile = "/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/ImpactParameterCalibration/IP_NonPrompt.root";
   f_np= new TFile(NameFile, "READ");
   
   if (f_np->IsZombie()) {
      std::cout << "File is not found." << std::endl;
      exit(-1);
   }

   f_p_ele= new TFile( "/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/ImpactParameterCalibration/IP_PromptElectron.root", "READ");
   
   if (f_p_ele->IsZombie()) {
      std::cout << "File dO_1.root is not found." << std::endl;
      exit(-1);
   }
   
   f_p_muon = new TFile( "/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/ImpactParameterCalibration/IP_PromptMuon.root", "READ");

   if (f_p_muon->IsZombie()) {
      std::cout << "File dO_1.root is not found." << std::endl;
      exit(-1);
   }
   

}

CalibrationOfImpactParameters::~CalibrationOfImpactParameters()
{
   //f_eta_np->Close();
   //f_pt_np->Close();
  
   f_np->Close();
   f_p_ele->Close();
   f_p_muon->Close();
}

void CalibrationOfImpactParameters::DoCalibrationForNonPromptLeptons(Float_t pt,Float_t eta,TString IP_Name, Float_t IP, Float_t &IP_cal){

   bool debug =false;
   if (debug) std::cout<<"In Calibration for electrons"<<std::endl;
   
   static TH1D *Data=new TH1D("","",50,-0.5,0.5);; 
   static TH1D *MC=new TH1D("","",50,-0.5,0.5);;
   

   TString NameData;
   TString NameMC;
   // if (10 <= pt && pt < 15) {
   //    NameData = IP_Name+"_pT_10to15_Data";
   //    NameMC = IP_Name+"_pT_10to15_MC";
   // }
   // if (15 <= pt && pt < 20) {
   //    NameData = IP_Name+"_pT_15to20_Data";
   //    NameMC = IP_Name+"_pT_15to20_MC";
   // }
   // if (20 <= pt && pt < 25) {
   //    NameData = IP_Name+"_pT_20to25_Data";
   //    NameMC = IP_Name+"_pT_20to25_MC";
   // }
   // if (25 <= pt && pt < 30) {
   //    NameData = IP_Name+"_pT_25to30_Data";
   //    NameMC = IP_Name+"_pT_25to30_MC";
   // }
   // if (30 <= pt && pt < 50) {
   //    NameData = IP_Name+"_pT_30to50_Data";
   //    NameMC = IP_Name+"_pT_30to50_MC";
   // }
   // if (50 <= pt) {
   //    NameData = IP_Name+"_pT_50_Data";
   //    NameMC = IP_Name+"_pT_50_MC";
   // }
   // Data = (TH1D*)f_pt_np->Get(NameData);
   // MC = (TH1D*)f_pt_np->Get(NameMC);

   // if (IP_Name =="d0_2" || IP_Name == "dZ_2"){
   //    if (fabs(eta) <= 0.9) {
   //       NameData = IP_Name+"_eta_0to09_Data";
   //       NameMC = IP_Name+"_eta_0to09_MC";
   //    }
   //    if (fabs(eta) > 0.9 && fabs(eta) <= 1.2) {
   //       NameData = IP_Name+"_eta_09to12_Data";
   //       NameMC = IP_Name+"_eta_09to12_MC";
   //    }
   //    if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) {
   //       NameData = IP_Name+"_eta_12to21_Data";
   //       NameMC = IP_Name+"_eta_12to21_MC";
   //    }
   //    if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) {
   //       NameData = IP_Name+"_eta_21to24_Data";
   //       NameMC = IP_Name+"_eta_21to24_MC";
   //    }
   // }

   // if (IP_Name =="d0_1" || IP_Name == "dZ_1"){
   //    if (fabs(eta) <= 0.9) {
   //       NameData = IP_Name+"_eta_0to09_Data";
   //       NameMC = IP_Name+"_eta_0to09_MC";
   //    }
   //    if (fabs(eta) > 0.9 && fabs(eta) <= 1.2) {
   //       NameData = IP_Name+"_eta_09to12_Data";
   //       NameMC = IP_Name+"_eta_09to12_MC";
   //    }
   //    if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) {
   //       NameData = IP_Name+"_eta_12to21_Data";
   //       NameMC = IP_Name+"_eta_12to21_MC";
   //    }
   //    if (fabs(eta) > 2.1 && fabs(eta) <= 2.5) {
   //       NameData = IP_Name+"_eta_21to25_Data";
   //       NameMC = IP_Name+"_eta_21to25_MC";
   //    }
   // }
   

   // std::cout<<"Name Data: "<<NameData<<std::endl;
   // std::cout<<"Name MC: "<<NameMC<<std::endl;
   // Data = (TH1D*)f_eta_np->Get(NameData);
   // MC = (TH1D*)f_eta_np->Get(NameMC);
   
   
   NameData = IP_Name+"_Data";
   NameMC = IP_Name+"_MC";
   
   Data = (TH1D*)f_np->Get(NameData);
   MC = (TH1D*)f_np->Get(NameMC);
      
   
   //To DO: Add check that Integral is One!
   if (Data ==NULL){
      std::cout << "Data histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
   if (MC ==NULL){
      std::cout << "MC histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
      
   double_t integral = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral != ((double)1.)){
      std::cout << "Data histogram not normalized to unity."<< std::endl;
      exit(-1);
   }
   double_t integral_mc = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral_mc !=1.){
      std::cout << "MC histogram not normalized to unity."<< std::endl;
      exit(-1);
   }

   IP_cal = Calibrate(Data, MC, IP);
   if (debug) std::cout<<"Calibrated Variable is : "<<IP_cal<<std::endl;
   
   
   return;
}


void CalibrationOfImpactParameters::DoCalibrationForPromptElectrons(Float_t pt,Float_t eta,TString IP_Name, Float_t IP, Float_t &IP_cal){

   bool debug =false;
   if (debug) std::cout<<"In Calibration for prompt electrons"<<std::endl;
   
  
   
   static TH1D *Data=new TH1D("","",50,-0.5,0.5); 
   static TH1D *MC=new TH1D("","",50,-0.5,0.5);;
   
   TString NameData;
   TString NameMC;
      
   NameData = IP_Name+"_Data";
   NameMC = IP_Name+"_MC";
   
   Data = (TH1D*)f_p_ele->Get(NameData);
   MC = (TH1D*)f_p_ele->Get(NameMC);

   // static TH1D *Data_2=new TH1D("","",50,-0.5,0.5); 
   // static TH1D *MC_2=new TH1D("","",50,-0.5,0.5);

   // if (IP_Name=="d0_1") {
   //    Data_2 = (TH1D*)f_p_ele->Get("d0_3_Data");
   //    MC_2= (TH1D*)f_p_ele->Get("d0_3_MC");
   // }
   // if (IP_Name=="dZ_1")
   // {
   //    Data_2 = (TH1D*)f_p_ele->Get("dZ_3_Data");
   //    MC_2= (TH1D*)f_p_ele->Get("dZ_3_MC");
   // }
   // Data->Add(Data_2);
   // MC->Add(MC_2);

   MC->Scale(1./MC->Integral());
   Data->Scale(1./Data->Integral());
   
   if (Data ==NULL){
      std::cout << "Data histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
   if (MC ==NULL){
      std::cout << "MC histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
      
   double_t integral = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral != ((double)1.)){
      std::cout << "Data histogram not normalized to unity."<< std::endl;
      exit(-1);
   }
   double_t integral_mc = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral_mc !=1.){
      std::cout << "MC histogram not normalized to unity."<< std::endl;
      exit(-1);
   }

   IP_cal = Calibrate(Data, MC, IP);
   if (debug) std::cout<<"Calibrated Variable is : "<<IP_cal<<std::endl;

   return;
}

void CalibrationOfImpactParameters::DoCalibrationForPromptMuons(Float_t pt,Float_t eta,TString IP_Name, Float_t IP, Float_t &IP_cal){

   bool debug =false;
   if (debug) std::cout<<"In Calibration for prompt muons"<<std::endl;
   
   static TH1D *Data=new TH1D("","",50,-0.5,0.5); 
   static TH1D *MC=new TH1D("","",50,-0.5,0.5);;
   
   TString NameData;
   TString NameMC;
      
   NameData = IP_Name+"_Data";
   NameMC = IP_Name+"_MC";
   
   Data = (TH1D*)f_p_muon->Get(NameData);
   MC = (TH1D*)f_p_muon->Get(NameMC);

   // static TH1D *Data_2=new TH1D("","",50,-0.5,0.5); 
   // static TH1D *MC_2=new TH1D("","",50,-0.5,0.5);

   // if (IP_Name=="d0_2") {
   //    Data_2 = (TH1D*)f_p_muon->Get("d0_4_Data");
   //    MC_2= (TH1D*)f_p_muon->Get("d0_4_MC");
   // }
   // if (IP_Name=="dZ_2") {
   //    Data_2 = (TH1D*)f_p_muon->Get("dZ_4_Data");
   //    MC_2= (TH1D*)f_p_muon->Get("dZ_4_MC");
   // }

   // Data->Add(Data_2);
   // MC->Add(MC_2);

   MC->Scale(1./MC->Integral());
   Data->Scale(1./Data->Integral());
   
   if (Data ==NULL){
      std::cout << "Data histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
   if (MC ==NULL){
      std::cout << "MC histogram used for quatile mapping not found."<< std::endl;
      exit(-1);
   }
      
   double_t integral = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral != ((double)1.)){
      std::cout << "Data histogram not normalized to unity."<< std::endl;
      exit(-1);
   }
   double_t integral_mc = round(Data->Integral() * 10000.0 ) / 10000.0;
   if (integral_mc !=1.){
      std::cout << "MC histogram not normalized to unity."<< std::endl;
      exit(-1);
   }

   IP_cal = Calibrate(Data, MC, IP);
   if (debug) std::cout<<"Calibrated Variable is : "<<IP_cal<<std::endl;

   return;
}




Float_t CalibrationOfImpactParameters::Calibrate(TH1D * Data, TH1D* MC, Float_t IP){
   
   bool debug =false;
   if (debug) std::cout<<"Do the Calibration!"<<std::endl;
   
   //calculated cumulative distribution in Data and MC, compute quantile that gives you corrections
   int i = MC ->FindBin(IP);
   
   if (debug) {
      std::cout<<"Bin is "<<i<<std::endl;
      std::cout<<"Integral is "<< MC->Integral(1,i-1)<<std::endl;
      std::cout<<"Bin content is "<< MC->GetBinContent(i)<<std::endl;
      std::cout<<"Binlow edge is "<<MC->GetXaxis()->GetBinLowEdge(i) <<std::endl;
      std::cout<<"Bin width is "<<MC->GetXaxis()->GetBinWidth(i) <<std::endl;
   }
   
   //double F_ip_sim = MC->Integral(1,i-1) + MC->GetBinContent(i) * (IP -MC->GetXaxis()->GetBinLowEdge(i))/MC->GetXaxis()->GetBinWidth(i) ; 
   double F_ip_sim = MC->Integral(1,i-1) + MC->GetBinContent(i) * (IP -MC->GetXaxis()->GetBinLowEdge(i))/MC->GetXaxis()->GetBinWidth(i) ; 
   
   double a = (MC->GetBinContent(i+1)-MC->GetBinContent(i-1))/MC->GetXaxis()->GetBinWidth(i);
   double b = MC->GetBinContent(i-1) - a * MC->GetXaxis()->GetBinLowEdge(i);
   double y_IP = a*IP + b;
   double F_Tr = ((IP-MC->GetXaxis()->GetBinLowEdge(i)) * (y_IP - MC->GetBinContent(i-1))) /2;
   
   F_ip_sim = F_ip_sim - F_Tr;
                                                                                               
   


   if (debug) std::cout<<"F_ip_sim "<<  F_ip_sim <<std::endl;
   
   // if (debug){
   // double_t sum =0;
   // int foundbin = 0;
   // for (unsigned int j =1; j<Data->GetNbinsX()+1; j++){
   //    sum=0;
   //    sum = Data->Integral(1,j);
   //    std::cout<<"Sum is : "<<sum<<" at bin "<<j<<std::endl;
   //    if (sum>F_ip_sim) {
   //       foundbin = j-1;
   //       std::cout<<"Sum above is : "<<sum<<std::endl;
   //       break;
   //    }
   // }
   // std::cout<<"Foundbin is: "<<foundbin<<std::endl;
   // double corr_IP = ((F_ip_sim -Data->Integral(1,foundbin-1)) * Data->GetXaxis()->GetBinWidth(foundbin)/Data->GetBinContent(foundbin))  + Data->GetXaxis()->GetBinLowEdge(foundbin);
   // std::cout<<"Cross check: corrected IP: "<<corr_IP<<std::endl;
   // }
   
   Int_t nprobSum = 1;
   double_t probSum[1] = {F_ip_sim};
   Double_t q[1];
   
   Data->GetQuantiles(nprobSum, q, probSum);
   return q[0];

}

//input: Histograms binned in pt and eta, TH1F for prompt / non-prompt electrons/ muons, data histograms: background subtracted from MC











