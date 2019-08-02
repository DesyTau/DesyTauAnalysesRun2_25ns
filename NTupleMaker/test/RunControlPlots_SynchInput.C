#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
/*
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"*/

/*
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
*/

//xsec_lumi_weight

//Macro to create the control plots serially:
/*
 pt_tt           = 160.873
 m_vis           = 31.947
 mt_tot          = 19.3558
 m_sv            = -9999
 pt_sv           = -9999
 eta_sv          = -9999
 phi_sv          = -9999
 met_sv          = -9999
 mt_sv           = -9999
 mjj             = 142.255
 jdeta           = 1.5279
 dijetpt         = 101.305
 jdphi           = 1.53945
*/

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"


void RunControlPlots_SynchInput(){
//Plot_lept_mutau_DNN

//think we need to specify observable and axis range

TString aap="mt_tot";
string aap2="mt_tot";

//gROOT->ProcessLine(".L Plot_lept_mutau_DNN.C++(%s)");
//gROOT->ProcessLine(".L Plot_lept_mutau_DNN.C++");
cout<<"here 1"<<endl;

TString cmd="Plot_lept_mutau_DNN(";

//Plot_lept_mutau_DNN(aap);
//gROOT->ProcessLine("Plot_lept_mutau_DNN()",TString("mt_tot"));


//gROOT->ProcessLine(Form(".x Plot_lept_mutau_DNN.C++","mt_tot"));
//Plot_lept_mutau_DNN.C++

//gROOT->ProcessLine(Form(".x Plot_lept_mutau_DNN.C++(%i)",1));

gROOT->ProcessLine(".L Plot_lept_mutau_Updated_Seq.C++");

for(int i=7; i<19;i++) gROOT->ProcessLine(Form("Plot_lept_mutau_Updated_Seq(%i)",i));


//gROOT->ProcessLine(Form(".x Plot_lept_mutau_DNN.C++(%s)",aap2)); 

//root 'Plot_lept_mutau_DNN.C("mt_tot")'

  return;
}

/*
void Plot_lept_mutau_DNN(TString Variable = "m_vis",
			     TString xtitle = "m_{vis} [GeV]",
			     int nBins  =   30,
			     float xmin =    0,
			     float xmax =  300,
			     TString Weight = "puweight*effweight*mcweight*",
			     TString ytitle = "Events",
			     TString directory = "./Inputs/NTuples_mt_Klundert_V2_2017/",//"./mutau_2019_4_8/",
			     TString outputDir = "./Plots/",
			     int year=2017,
			     bool logY = false
			     ){
*/
