
#include <string.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "TList.h"
#include <string>
#include "TObject.h"
#include "TBranch.h"
#include <functional>
#include "TAxis.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Ranker.h"

using namespace std;




int main (int argc, char *argv[])
{

 TTree *tSignal, *tBkg; TFile *hfileS, *hfileB;
      bool over_table;bool correl_table;bool over_cleaned_sorted_table,first_loop,DoCors;
      DoCors=false;
      over_table=true;
      correl_table=true;
      over_cleaned_sorted_table= true;
      first_loop=false;
      string merged_file_ ;
      string sm_file_ ;
      string np_file_;
      char sgnl_file_char_[100];
      char bkg_file_char_[100];

vector <float> Correl_Cut_;
  Correl_Cut_.clear();
  /*for (unsigned int j=0;j<4;j++){
    float v_cut=0.5 + 0.1*j;
    Correl_Cut_.push_back(v_cut);
    }*/
  Correl_Cut_.push_back(0.8);


      first_loop=true;
      TString era=argv[3];

      //merged_file_ = "Merged_Plots_" + datasets[0]->Title () + "_" + datasets[d]->Title () + ".root";
      merged_file_ = era+"/"+"Merged_Plots_"+argv[2]+".root";
      //string target_file_ = datasets[0]->Title () + postfix+"_tree.root" ;
      string target_file_ = "SM_tree.root" ;
      TFile *merged_file= TFile::Open(merged_file_.c_str(),"read");
      TFile *target_file = TFile::Open(target_file_.c_str(),"read");
 /*     if (!merged_file)
        cout<<" Merged_file "<<merged_file_<<"  DOES NOT EXIST..... ========================== Will be glad to make it for you ;-)"<<endl;
      if (merged_file || !target_file)
        {
          cout<<" Merged File "<<merged_file_<<"  DOES EXIST..... =============================== Will go to the next dataset"<<endl;
        }
*/
      sm_file_ = era+"/"+argv[1];
      np_file_ = era+"/"+argv[2];
      cout << " Will merge and compute overlap for " << sm_file_ << "  and  " << np_file_ << "  file will be " << merged_file_.c_str () << endl;
      sprintf (bkg_file_char_, era +"/"+ argv[1] + ".root");
      sprintf (sgnl_file_char_, era + "/" + argv[2] + ".root");
      hfileS = new TFile (sgnl_file_char_);
      hfileB = new TFile (bkg_file_char_);
      cout << "Filename for signal: " << sgnl_file_char_ << "  - Filename for background: " << bkg_file_char_ << endl;
      //tSignal = (TTree *) hfileS->Get ("T");
      //tBkg = (TTree *) hfileB->Get ("T");
      //tSignal->Print();
       
      Ranker (tSignal, tBkg, merged_file_, np_file_, sm_file_, over_table, correl_table, over_cleaned_sorted_table, DoCors, Correl_Cut_);

}
