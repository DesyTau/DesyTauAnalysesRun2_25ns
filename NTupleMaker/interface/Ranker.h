//#include "../interface/Ranker.h"
//#include "../interface/Observables.h"
#include "TCut.h"
#include "math.h"
#include <algorithm>
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
#include "TDirectory.h"
#include "TBranch.h"
#include <functional>
#include "TAxis.h"


using namespace std;


  void ComputeOverlap(string &fout, string &merged_file,string &np_file);
  void ModifyHist(TH1D* &h, Color_t lcolor);
  void Impose( TDirectory *ttarget, TList *ssourcelist, string &np_legend );
  void PrintOverlapTable(bool & table_overlap_);
  void PrintCorrelationTable(bool &table_correlation_ , float & Cor_Cut);
  void DoOnlyCorrelation(bool &table_correlation_ );
  void CorrelationBasedOnOverlap(bool &table_correlation_ );
  void CorrelationBasedOnOverlap(bool &table_correlation_ ,float & Cor_Cut);
  void DoCorrelation(bool &table_correlation_ ,  float &Cor_Cut);
  void PrintOnlyCorrelation(bool &table_correlation_ , float & Cor_Cut);
  void PrintOverlapSortedCleanedTable(bool &table_overlap_sorted_cleaned_, string & fin, string & fout);
  void PrintOverlapSortedCleanedTable(bool &table_overlap_sorted_cleaned_);
  void SaveBestObservablesList(string &fout, float & Cor);
  void SaveOverlapList(string &fout);
  void SaveCorrelanceObservablesList(string &fout);
  void Read(string & fout, float & Cor);
  void ReadOverlap(string & fout);


struct sort_pred {
    bool operator()(const std::pair<string,float> &left, const std::pair<string,float> &right) {
	    //if(left.second!=right.second)
	//	return left.second < right.second;
	//        return left.first < right.first;
        return left.second < right.second;
    }
};


struct sort_string_int {
    bool operator()(const std::pair<string,int> &left, const std::pair<string,int> &right) {
        return left.second < right.second;
    }
};


struct sort_two_int {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }


bool cmp(const pair<string, long> &p1, const pair<string, long> &p2)
{
	    if(p1.second!=p2.second)
		            return p1.second < p2.second;
	        return p1.first < p2.first;
}


};



 float GetRankedVariableByPos(int position);

 float GetRankedVariableByName(const string);

 string GetRankedVariableName(int position);

 string GetVariableName(int position);

 string GetVarsForPseudo(int position);

 int GetRankingofVariable(const string );

 float GetCorrelation(const string, const string);



  Double_t CS;
  Double_t bincontS;
  Double_t bincontB;

   //TTrees
   TTree* tree1_;
   TTree* tree2_;
   TTree * tree_Signal;
   TTree * tree_Bkg;
   //overlap
   vector<pair<string , float> > overlaps_;
   vector<pair<string , float> > overlaps_sorted_;
   vector<pair<string, float> >  overlaps_sorted_cleaned_; 
   vector<pair<string, float> >  overlaps_sorted_pro_cleaned_; 
   vector <string> overlaps_sorted_cleaned_name_;
   vector <float> overlaps_sorted_cleaned_value_;
   //correlations
   vector<vector<float> >   correlations_;
   vector<vector<float> >   correlations_sorted_;
   vector<string> variables_;  
   vector<float> correlation_;
   vector<float> correlation_sorted_;
   vector<pair<string, float> > correlation_with_name_;
   vector<pair<string, float> > correlation_with_name_sorted_;
   vector<pair<string, float> > all_correlation_table_;
   
   vector<pair<string, int> > high_correlation_position_;
   vector<pair<string,string> > high_correlation_matrix_;
   vector<float> correlation_sorted_cleaned_;
   vector<vector<float> > correlations_sorted_cleaned_;

void Ranker (string & sm_file_, bool &tbl_correlation,  bool & DoCorrs)
{
  if (DoCorrs){
    DoOnlyCorrelation(tbl_correlation);
    SaveCorrelanceObservablesList(sm_file_);
  }
}


void Ranker (  bool & DoCorrs)
{
 cout<<" ===================================== END  "<<endl;

}

/*
void ~Ranker ()
{
  delete tree_Signal;;
  delete tree_Bkg;
}
*/


void Ranker (TTree * treeSignal, TTree * treeBkg, string & merged_file_, string & np_file_, string & sm_file_, bool &tbl_overlap, bool &tbl_correlation, bool &tbl_overlap_sorted_cleaned , bool & DoCorrs,vector <float> cor_cut_)
{
  /*
    treeBkg->SetBranchStatus("*Had*",0);
    treeBkg->SetBranchStatus("*Lep*",0);
    treeBkg->SetBranchStatus("*BTag*",0);
    treeBkg->SetBranchStatus("*Boosted*",0);
    treeBkg->SetBranchStatus("*Chi2*",0);
    treeBkg->SetBranchStatus("*Ttbar*",0);
 
    treeSignal->SetBranchStatus("*Had*",0);
    treeSignal->SetBranchStatus("*Lep*",0);
    treeSignal->SetBranchStatus("*BTag*",0);
    treeSignal->SetBranchStatus("*Boosted*",0);
    treeSignal->SetBranchStatus("*Chi2*",0);
    treeSignal->SetBranchStatus("*Ttbar*",0);
  */

  ComputeOverlap (merged_file_, np_file_, sm_file_);
  tree_Signal = treeSignal;
  tree_Bkg = treeBkg;
  if (DoCorrs){

	 
    DoOnlyCorrelation(tbl_correlation);
    CorrelationBasedOnOverlap(tbl_correlation);
    SaveCorrelanceObservablesList(sm_file_);
  }

  for (unsigned int i=0;i<cor_cut_.size();i++){

    //correlations
    correlations_.clear();
    correlations_sorted_.clear();
    correlation_.clear();
    correlation_sorted_.clear();
    correlation_with_name_.clear();
    correlation_with_name_sorted_.clear();
    all_correlation_table_.clear();
    high_correlation_position_.clear();
    high_correlation_matrix_.clear();
    correlation_sorted_cleaned_.clear();
    correlations_sorted_cleaned_.clear();

   // DoCorrelation(tbl_correlation,cor_cut_[i]);
    PrintOverlapTable (tbl_overlap);
    //PrintCorrelationTable (tbl_correlation, cor_cut_[i]);
   // PrintOverlapSortedCleanedTable(tbl_overlap_sorted_cleaned);
    //SaveBestObservablesList(merged_file_,cor_cut_[i]);
    SaveOverlapList(merged_file_);
    ReadOverlap(merged_file_);
  }
}



void Ranker (string & merged_file_, string & np_file_, string & sm_file_, bool &tbl_overlap, bool &tbl_correlation, bool &tbl_overlap_sorted_cleaned, float & cor_cut)
{


  ComputeOverlap (merged_file_, np_file_, sm_file_);
   
  //for (int i=0;i<cor_cut_.size();i++){
  DoCorrelation(tbl_correlation,cor_cut);
  PrintOverlapTable (tbl_overlap);
  PrintCorrelationTable (tbl_correlation, cor_cut);
  PrintOverlapSortedCleanedTable(tbl_overlap_sorted_cleaned);
  //}
}






void
SaveCorrelanceObservablesList(string & fout){
  TFile * target;
  string file = fout+"_tree.root";
  target = TFile::Open(file.c_str(),"update");
  target->cd();
  
  TTree *Tree_Cor = new TTree("CorrelanceVector","List of Correlance for Observables");
  
  Tree_Cor->Branch("all_correlation_table_",&all_correlation_table_);
  //  for (unsigned int i=0;i<all_correlation_table_.size();i++){
  //cout<<" all------------ "<<all_correlation_table_[i].first<<"   "<<fout<<endl;}
  
  Tree_Cor->Fill();
  target->Write();
  target->Close();
  delete target;
}

void
SaveBestObservablesList(string & fout, float & cor){
  TFile * target;
  target = TFile::Open(fout.c_str(),"update");
  target->cd();
  char tree_name[50];
  sprintf(tree_name,"ListVector with Cor of %f",cor);
  TTree *t = new TTree(tree_name,"List of Sorted and Cleaned Observables");
  t->Branch("overlaps_sorted_cleaned_",&overlaps_sorted_cleaned_);
  t->Fill();
  target->Write();
  target->Close();
  delete target;
}

void
SaveOverlapList(string & fout){
  TFile * target;
  target = TFile::Open(fout.c_str(),"update");
  target->cd();
  char tree_name[50];
  TTree *t = new TTree("Overlap","Overlap List");
  t->Branch("overlaps_sorted_",&overlaps_sorted_);
  t->Fill();
  target->Write();
  target->Close();
  delete target;
}

void
ReadOverlap(string &fout){
  TFile *f = TFile::Open(fout.c_str(),"read");
  TTree *t;
  f->GetObject("Overlap",t);
  TBranch *br_list=0;
  vector<pair<string, float> >* myListP=&overlaps_sorted_;
  t->SetBranchAddress("overlaps_sorted_",&myListP,&br_list);
  Long64_t tentry = t->LoadTree(0);
  br_list->GetEntry(tentry);
  for (unsigned int j=0;j<overlaps_sorted_.size();++j){
    cout<<" NOW READ OUT THE VECTOR  "<<overlaps_sorted_[j].first<<"    "<<overlaps_sorted_[j].second<<endl;
  }
  delete t;
  f->Close();
  delete f;


}
	    
void
Read(string &fout,float & cor){
  TFile *f = TFile::Open(fout.c_str(),"read");
  TTree *t;
  char tree_name[50];
  sprintf(tree_name,"ListVector with Cor of %f",cor);
  f->GetObject(tree_name,t);
  TBranch *br_list=0;
  vector<pair<string, float> >* myListP=&overlaps_sorted_cleaned_;
  t->SetBranchAddress("overlaps_sorted_cleaned_",&myListP,&br_list);
  Long64_t tentry = t->LoadTree(0);
  br_list->GetEntry(tentry);
  for (unsigned int j=0;j<overlaps_sorted_cleaned_.size();++j){
    cout<<" NOW READ OUT THE VECTOR  "<<overlaps_sorted_cleaned_[j].first<<"    "<<overlaps_sorted_cleaned_[j].second<<endl;
  }
  delete t;
  f->Close();
  delete f;


}

void
ComputeOverlap (string & fout, string & np_file, string & sm_file)
{
  gROOT->SetStyle ("Plain");
  gStyle->SetPalette (1);

  TList *FileList;
  TFile *Target;
  string fileout = fout;
  string file_np = np_file + ".root";
  string file_sm = sm_file + ".root";
  FileList = new TList ();
  FileList->Add (TFile::Open (file_np.c_str ()));
  FileList->Add (TFile::Open (file_sm.c_str ()));

  //TFile *Target = new TFile (fileout.c_str (), "RECREATE");
  Target = TFile::Open (fout.c_str (), "RECREATE");

  // more than 2 source files: not tested
  string np_title = np_file;
  Impose (Target, FileList, np_title);
  delete FileList;
  delete Target;


}


void
Impose (TDirectory * target, TList * sourcelist, string & np_title_)
{
  cout << "	" << "========================================================" << endl;
  cout << "	" << "This is a macro to superimpose plots of different root files." << endl;
  cout << "	" << "Only TH1D objects are superimposed." << endl;
  cout << "	" << "Target path: " << target->GetPath () << endl;
  TString path ((char *) strstr (target->GetPath (), ":"));
  path.Remove (0, 2);
  std::vector < pair < string, float > > BinContent_;


  overlaps_.clear();
  overlaps_sorted_.clear();

  TFile *first_source = (TFile *) sourcelist->First ();
  first_source->cd ();
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus ();
  TH1::AddDirectory (kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey (current_sourcedir->GetListOfKeys ());
  TKey *key, *oldkey = 0;
  TCanvas *c1 = new TCanvas ("c1", "c1", 500, 500);
  while ((key = (TKey *) nextkey ())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp (oldkey->GetName (), key->GetName ()))
      continue;

    // read object from first source file and create a canvas
    first_source->cd ();

    TObject *obj = key->ReadObj ();
    string nn   = obj->GetName ();

    if (std::string::npos == nn.find("_14")) continue;
    c1->SetName(obj->GetName ());

    if (!obj->IsA ()->InheritsFrom ("TH1D") )  continue;
    if (obj->IsA ()->InheritsFrom ("TH1D") ) {
      // descendant of TH1D -> prepare the histograms to be superimposed

  	 
        int counter = 0;
      TH1D *h1 = (TH1D *) obj;
      ModifyHist (h1, kBlue);
 //           cout << "Modifying histogram " << obj->GetName() << endl;      

      // loop over all source files and modify the correspondant
      // histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile *) sourcelist->After (first_source);
      while (nextsource) {
        counter=0;
	// make sure we are at the correct directory level by cd'ing to path
	nextsource->cd ();
	TKey *key2 = (TKey *) gDirectory->GetListOfKeys ()->FindObject (h1->GetName ());
	if (key2) {
		counter++;
	//	if (counter>1) break;
	  TH1D *h2 = (TH1D *) key2->ReadObj ();
	  ModifyHist (h2, kRed);
           // cout << "Modifying histogram for background  " << obj->GetName() << endl;      
	  TH1D* htest;
	  double maxh1;
	  double maxh2;

	  cout<<"        ============================= "<<h1->GetTitle() <<"    "<< h2->GetTitle()<<"    "<<h1->GetName()<<"   "<<h2->GetName()<<endl;

	  CS = 0.;
	  if (h1->GetNbinsX() != h2->GetNbinsX() ) break;

	 // if (h1->Integral() == 0 ) continue;
	 // if (h2->Integral() == 0 ) continue;
	  string nstr = h1->GetName();
	  nstr.erase(nstr.size()-2);
	  h1->SetName(nstr.c_str());

	  nstr = h2->GetName();
	  nstr.erase(nstr.size()-2);
	//  h2->SetName(nstr.c_str());

	  //cout <<" new name "<<nstr.c_str()<<"  "<<h1->GetName()<<"  "<<h2->GetName()<<"  "<<h1->GetNbinsX()<<endl;
          h1->SetStats(00000);
          h2->SetStats(00000);
	  htest=h2;
	   //overlaps_.clear();overlaps_sorted_.clear();
	  for (Int_t b = 0; b < htest->GetNbinsX ()+1 ; b++) {
	    bincontS = float(h2->GetBinContent (b));
	    bincontB = float(h1->GetBinContent (b));
	    //if (bincontS ==0 || bincontB==0) CS=1;
            if (bincontS >0 && bincontB>0){
	    
	      if (bincontS < bincontB)
		CS += bincontS;
	      if (bincontS >= bincontB)
		CS += bincontB;
	      //cout<<" ==========================    BIN CONTENT h1 name "<<h1->GetName()<<" h2 name "<<h2->GetName()<<" cont_S "<<bincontS<<" cont_B  "<<bincontB<<" BIN "<<b<<" BINS_H1 "<<h1->GetNbinsX()<<" BINS_H2  "<<h2->GetNbinsX()<<" CS  "<<CS<<endl; 	  	
	    }
	    //else CS=1;
	  }
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);
	h1->SetLineStyle(2);
	h2->SetLineStyle(2);
	h2->SetMinimum(0.001);
	h1->SetMinimum(0.001);
	
	//h2->Add(h1,-1);
	//CS = h1->Integral();
	  //if (CS>0. && CS<1.){
	  //if (CS<1.){
	  overlaps_.push_back (pair < string, float >((htest->GetName ()), CS));
	  overlaps_.erase( unique( overlaps_.begin(), overlaps_.end() ), overlaps_.end() );
	  //if (CS>0)
	  overlaps_sorted_.push_back (pair < string, float >((htest->GetName ()), CS));
	  overlaps_sorted_.erase( unique( overlaps_sorted_.begin(), overlaps_sorted_.end() ), overlaps_sorted_.end() );
	  // overlaps_sorted_cleaned_.push_back(pair< string, float> ((h1->GetName()), CS));
	  maxh1 = h1->GetMaximum (10000.);
	  maxh2 = h2->GetMaximum (10000.);
	  //maxh1 = h1->GetBinCenter(h1->GetNbinsX()+1);
          //maxh2 = h2->GetBinCenter(h2->GetNbinsX()+1);
	  if (maxh1 > maxh2 ) {
	    h1->Draw ("hist text");
	    h2->Draw ("histsames");
	    h1->SetMaximum (1.5*h1->GetMaximum());
	    h2->SetMaximum (1.5*h1->GetMaximum());
	  }
	  else {
	    h2->Draw ("hist text");
	    h1->Draw ("histsames");
	    h1->SetMaximum (1.5*h2->GetMaximum());
	    h2->SetMaximum (1.5*h2->GetMaximum());

	  }
	//h1->Draw();
	  TLegend *legend_c1 = new TLegend (0.55, 0.80, 0.9, 0.70);
	  legend_c1->SetTextFont (22);
	  legend_c1->SetTextSize (0.03);

	  //string new_np_title = np_title_.substr(0, np_title_.find("Skimmed/"));
	  //string new_np_title = np_title_.substr(0, np_title_.find("Skimmed/"));
	  np_title_.erase(0,8);//= np_title_.substr(0, np_title_.find("Skimmed/"));
	  
	  legend_c1->AddEntry (h1, "NP", "L");	// change this aBranch_Bkg_ording to first source file
	  legend_c1->AddEntry (h2, "SM", "L");	// change this aBranch_Bkg_ording to second source file
	  legend_c1->Draw ("SAMES");
		//}//CS>0


	}
	nextsource = (TFile *) sourcelist->After (nextsource);
      }
      // while ( nextsource )
    }

    /*
    else if (obj->IsA ()->InheritsFrom ("TTree")) {	// not tested

      // loop over all source files create a chain of Trees "globChain"
      const char *obj_name = obj->GetName ();

      globChain = new TChain (obj_name);
      globChain->Add (first_source->GetName ());
      TFile *nextsource = (TFile *) sourcelist->After (first_source);
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
      while (nextsource) {

	globChain->Add (nextsource->GetName ());
	nextsource = (TFile *) sourcelist->After (nextsource);
      }

    }
    */

    /*else if (obj->IsA ()->InheritsFrom ("TDirectory")) {	// not tested
    // it's a subdirectory

    cout << "Found subdirectory " << obj->GetName () << endl;

    // create a new subdir of same name and title in the target file
    target->cd ();
    ////     TDirectory *newdir = target->mkdir (obj->GetName (), obj->GetTitle ());

    // newdir is now the starting point of another round of superimposing
    // newdir still knows its depth within the target file via
    // GetPath(), so we can still figure out where we are in the recursion
    // Impose( newdir, sourcelist );

    }
    else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " << obj->GetName () << ", object type: " << obj->ClassName () << endl;
    }
    */
    // now draw and write the superimposed histograms to the target file
    // note that this will just store the canvas c1 in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if (obj) {
      target->cd ();
      
      //!!if the object is a tree, it is stored in globChain...     
      if (obj->IsA ()->InheritsFrom ("TTree"))	// not tested
	globChain->Merge (target->GetFile (), 0, "keep");
      else
	c1->Write (key->GetName ());
    }

    // delete first_source;
    //   delete obj;
    //    delete legend_c1;
    //    delete c1;
 
  }				// while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf (kTRUE);
  TH1D::AddDirectory (status);
  cout << "	" << "========================================================" << endl;
  cout<< " Ended SuperImpose of files.... " <<endl;

  // for (int i=0;i<SValueVariables_.size();i++){
  //cout<<SValueVariables_[i].second<<endl;}


  /*   TObject *obj = key->ReadObj ();
       TCanvas *c1 = new TCanvas ("c1", obj->GetName (), 500, 500);
       TLegend *legend_c1 = new TLegend (0.65, 0.80, 0.89, 0.70);
  */



}


void
ModifyHist (TH1D * &h, Color_t lcolor)
{
  double temp_integral;

  h->SetLineColor (lcolor);
  h->SetFillColor (lcolor);
  h->SetMarkerColor (lcolor);
  h->SetMarkerSize (1.);
  if (lcolor== kBlue)  h->SetFillStyle (3004);
  if (lcolor== kRed)  h->SetFillStyle (3005);
  int nbins=h->GetNbinsX();
  int nn=1;

 

  float over_ = h->GetBinContent(nbins+1);
  float contlast = h->GetBinContent(nbins);
  h->SetBinContent(nbins,contlast+over_);

  if (  nbins==50)  nn=2.5;
    if (  nbins==100)  nn=4;
    if (  nbins==150)  nn=4;
    if (  nbins>150)  nn=5;
    if (  nbins>200)  nn=5;
    if (  nbins<30)  nn=1;
    if (  nbins==70)  nn=3.5;
    if (  nbins==80)  nn=4;
//    if (  nbins>249)  nn=5;

/*
  if( nbins%10==0)   nn=10;
  if( nbins%20==0)   nn=20;
  if( nbins%30==0)   nn=30;
*/
    //  }

	  h->Rebin(nn);
  temp_integral = h->Integral ();
  //cout << temp_integral << endl;
  //float lumi_scale=5000;
  // h->Scale (pow (lumi_scale, -1));
  if (temp_integral != 0)
   //h->Scale (pow (temp_integral, -1));
   h->Scale (1/h->Integral());
}


void
PrintOverlapTable (bool & cout_table_over)
{

 
  //sort (overlaps_.begin (), overlaps_.end ());
  sort (overlaps_sorted_.begin (), overlaps_sorted_.end (), sort_pred ());
  // sort (overlaps_sorted_cleaned_.begin (), overlaps_sorted_cleaned_.end (), sort_pred ());
  cout << " OVERLAPS AND CORRELATION SORTED OVERLAPS " << overlaps_sorted_.size () << "  " << correlations_.size () << endl;
  if (cout_table_over){
    cout << "\\begin{table}" << endl;
    cout << "\\centering" << endl;
    cout << "\\begin{tabular}{l|";
    for (unsigned int j = 0; j < overlaps_sorted_.size (); j++)


      cout << "c";
    cout << "}" << endl;
    cout << "\\hline" << endl;
    //fill first line of the table
    cout << "Variables & \t";
    //loop over overlaps
  }


  if (cout_table_over) cout << " Now will print the Sorted Table according to Overlap ascending ...." << endl;
  for (unsigned int i = 0; i < overlaps_sorted_.size (); i++) {
    if (cout_table_over)  cout << overlaps_sorted_[i].first << " &\t";


    for (unsigned int j = 0; j < overlaps_.size (); j++) {
      if (overlaps_[j].first == overlaps_sorted_[i].first) {
	//cout<< " FOUND THE PAIR , it was "<<"\t"<<overlaps_[j].first<<" and now is   "<<"\t"<<overlaps_sorted_[i].first<<" was in place  "<<"\t"<<j<<" and now is in place "<<"\t"<<i<<"\t"<<" with value  "<<overlaps_[j].second<<"  "<<overlaps_sorted_[i].second<<endl;
	//	correspondant_.push_back (pair < int, int >(j, i));
	//	correspondant_cleaned_.push_back (pair < int, int >(j, i));
      }
    }
  }
  if (cout_table_over){
    cout << " \\\\" << endl;
    cout << "\\hline" << endl;
    cout << "Overlap & \t";
    //loop over overlaps
    for (unsigned int i = 0; i < overlaps_sorted_.size (); i++) {

      cout << overlaps_sorted_[i].second << " &\t";
    }
    cout << " \\\\" << endl;
    cout << "\\hline" << endl;

    cout << "  Finished with overlap table.... " << endl;

  }

}


void
PrintOverlapSortedCleanedTable (bool & cout_table_srt_cln)
{
 
  //sort (overlaps_.begin (), overlaps_.end ());
  sort (overlaps_sorted_cleaned_.begin (), overlaps_sorted_cleaned_.end (), sort_pred ());
  sort (overlaps_sorted_cleaned_.begin (), overlaps_sorted_cleaned_.end (), sort_pred ());
  //  cout << " OVERLAPS AND CORRELATION SORTED AND CLEANED OVERLAPS " << overlaps_sorted_cleaned_.size () << "  " << correlations_.size () << endl;


  if (cout_table_srt_cln){
    cout << "\\begin{table}" << endl;
    cout << "\\centering" << endl;
    cout << "\\begin{tabular}{l|";
    for (unsigned int j = 0; j < overlaps_sorted_cleaned_.size (); j++)


      cout << "c";
    cout << "}" << endl;
    cout << "\\hline" << endl;
    //fill first line of the table
    cout << "Variables & \t";
    //loop over overlaps



    cout << " Now will print the Sorted and Cleaned Table according to Overlap ascending ...." << endl;
    for (unsigned int i = 0; i < overlaps_sorted_cleaned_.size (); i++) {
      cout << overlaps_sorted_cleaned_[i].first << " &\t";
  
    }
    cout << " \\\\" << endl;
    cout << "\\hline" << endl;
    cout << "Overlap & \t";
    //loop over overlaps
    for (unsigned int i = 0; i < overlaps_sorted_cleaned_.size (); i++) {

      cout << overlaps_sorted_cleaned_[i].second << " &\t";
    }
    cout << " \\\\" << endl;
    cout << "\\hline" << endl;

    cout << "  Finished with overlap sorted-cleaned table.... " << endl;

  }

}

void
DoOnlyCorrelation ( bool &cout_table_corr_srt_cln )
{



  TObjArray *Branches_Bkg = tree_Bkg->GetListOfBranches ();
 

  Int_t Entries_Bkg = Branches_Bkg->GetEntries ();
  

  for (Int_t i = 0; i < Entries_Bkg; i++) {

    TBranch *Branch_Sgnl_ = (TBranch *) Branches_Bkg->At (i);


    char Name_Bkg_[100];

    sprintf (Name_Bkg_, "%s", Branch_Sgnl_->GetName ());


   
    for (Int_t j = 0; j < Entries_Bkg; j++) {


      TBranch *Branch_Bkg_ = (TBranch *) Branches_Bkg->At (j);


      char Name_Bkg2_[100];
      sprintf (Name_Bkg2_, "%s", Branch_Bkg_->GetName ());
      tree_Signal->Draw (Form ("%s.Draw():%s.Draw()>>hcorr", Name_Bkg_, Name_Bkg2_), "goff");
      //tree_Signal->Draw (Form ("%s:%s>>hcorr", Name_Bkg_, Name_Bkg2_),"", "goff");
      TH2 *hcorr = (TH2 *) gDirectory->Get ("hcorr");
      cout << "The correlation from TTJets " << Name_Bkg_ << " with " << Name_Bkg2_ << "   correl factor =   " << '\t' << hcorr->GetCorrelationFactor () << "  " << endl;

      float hCorr = hcorr->GetCorrelationFactor ();
      //   correlation_sorted_.push_back (hCorr);
      string name =  Name_Bkg_; 

      all_correlation_table_.push_back(pair <string, float > (Name_Bkg_,hCorr));

    }
  }

}

void
CorrelationBasedOnOverlap ( bool &cout_table_corr_srt_cln )
{


  for (unsigned int h=0;h<overlaps_sorted_.size();h++){
    if (overlaps_sorted_[h].first=="weight")overlaps_sorted_.erase (overlaps_sorted_.begin()+h);}

  sort (overlaps_sorted_.begin (), overlaps_sorted_.end (), sort_pred ());
//  for (unsigned int h=0;h<overlaps_sorted_.size();h++){

//    cout<<"  Here we are... "<<h<<"  "<<overlaps_sorted_[h].first<<"   "<<overlaps_sorted_[h].second<<endl;
//  }

  /*
    TObjArray *Branches_Bkg = tree_Bkg->GetListOfBranches ();
    TObjArray *Branches_Sgn = tree_Signal->GetListOfBranches ();
 

    //Int_t Entries_Bkg = Branches_Bkg->GetEntries ();
    Int_t Entries_Bkg = overlaps_sorted_.size();
  

    //  for (Int_t i = 0; i < Entries_Bkg; i++) {
    for ( unsigned int i=0;i<overlaps_sorted_.size();i++){
    //for ( unsigned int i=0;i<5;i++){

    TBranch *Branch_Sgnl_ = (TBranch *) Branches_Sgn->At (i);


    char Name_Bkg_[100];

    sprintf (Name_Bkg_, "%s", Branch_Sgnl_->GetName ());

   
    //	     for (Int_t j = 0; j < Entries_Bkg; j++) {
    for ( unsigned int j=0;j<overlaps_sorted_.size();j++){
    //for ( unsigned int j=0;j<5;j++){


    TBranch *Branch_Bkg_ = (TBranch *) Branches_Sgn->At (j);


    char Name_Bkg2_[100];


    sprintf (Name_Bkg2_, "%s", Branch_Bkg_->GetName ());
	


    tree_Signal->Draw (Form ("%s:%s>>hcorr", Name_Bkg_, Name_Bkg2_, "", "goff"));
    TH2 *hcorr = (TH2 *) gDirectory->Get ("hcorr");
    //   cout << "The correlation from TTJets " << Name_Bkg_ << " with " << Name_Bkg2_ << "   correl factor =   " << '\t' << hcorr->GetCorrelationFactor () << "  " << endl;

    float hCorr = hcorr->GetCorrelationFactor ();
    //   correlation_sorted_.push_back (hCorr);
    string name =  Name_Bkg_; 

    all_correlation_table_.push_back(pair <string, float > (Name_Bkg_,hCorr));
    //cout<<"  all correl  "<<all_correlation_table_[all_correlation_table_.size()-1].first<<"   "<<all_correlation_table_[all_correlation_table_.size()-1].second<<endl;
    }
    }
  */
}


void
DoCorrelation ( bool &cout_table_corr_srt_cln , float & cor_cut)
{


  //  high_correlation_matrix_.clear();

  sort (overlaps_sorted_.begin (), overlaps_sorted_.end (), sort_pred ());

  TObjArray *Branches_Sgn = tree_Signal->GetListOfBranches ();
   
  //Int_t Entries_Sgn = Branches_Sgn->GetEntries ();
  Int_t Entries_Sgn = overlaps_sorted_.size();

    char Name_Sgnl_[100];
    char Name_Sgnl2_[100];
    char cors[100];
   char nn[100]; 

   cout<< " Entries  "<<Entries_Sgn<<"  "<<endl;
  for (int i = 0; i < Entries_Sgn; i++) {
		
    TBranch *Branch_Sgnl_ = (TBranch *) Branches_Sgn->At (i);
		
		
    sprintf (Name_Sgnl_, "%s", Branch_Sgnl_->GetName ());
		
    for (int j = 0; j < Entries_Sgn; j++) {
			
      TBranch *Branch_Sgnl2_ = (TBranch *) Branches_Sgn->At (j);

      sprintf (Name_Sgnl2_, "%s", Branch_Sgnl2_->GetName ());
       //tree_Signal->Draw (Form ("%s:%s>>hcorr", Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName (), "", "goff"));
      //tree_Signal->Draw (Form ("%s:%s>>hcorr_%i%i", Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName(),i,j,"","goff"));

      cout<<" Names  "<<Branch_Sgnl_->GetName ()<<endl;

        tree_Signal->Draw (  Form ("%s.Draw():%s.Draw()>>hcorr_%i%i", Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName(),i,j),"","goff");
    //  tree_Signal->Draw (Form ("%s:%s>>hcorr_%i%i", Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName(),i,j),"goff");
     
    //  sprintf(cors,"%s:%s>>%s",Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName ());

      //tree_Signal->Draw(cors);
      sprintf(nn,"hcorr_%i%i",i,j);
      //TH2 *hcorr_ = (TH2 *) gDirectory->Get ("hcorr");
      TH2 *hcorr_ = (TH2 *) gDirectory->Get (nn);
      
      float hCorrs = -100;hCorrs= hcorr_->GetCorrelationFactor ();
      cout << "The correlation is..." << Name_Sgnl_ << " with " << Name_Sgnl2_ << "   correl factor =   " << '\t' << hcorr_->GetCorrelationFactor () <<"  "<<nn<< endl;

      all_correlation_table_.push_back(pair <string, float > (Name_Sgnl_,float(hCorrs)));

      if (fabs (hCorrs) > cor_cut && Name_Sgnl_ != Name_Sgnl2_ && j > i) {
	high_correlation_matrix_.push_back (pair < string, string > (Name_Sgnl_,Name_Sgnl2_));
      }
     // delete hcorr_;
    }
  }
}



void
PrintCorrelationTable ( bool &cout_table_corr_srt_cln , float & cor_cut)
{



  high_correlation_position_.clear();
  overlaps_sorted_pro_cleaned_.clear();
  
  sort (overlaps_sorted_.begin (), overlaps_sorted_.end (), sort_pred ());

  correlation_with_name_sorted_.clear ();

  char Name_Bkg_[100];
  char Name_Bkg2_[100];
  
  for (unsigned int i = 0; i < high_correlation_matrix_.size(); i++) {


    sprintf (Name_Bkg_, "%s", high_correlation_matrix_[i].first.c_str());

    sprintf (Name_Bkg2_, "%s", high_correlation_matrix_[i].second.c_str());


    //for (unsigned int ii = 0; ii < high_correlation_matrix_.size(); ii++) {
    //if (Name_Bkg_==Name_Bkg2_) cout <<"  WARNINGGGGGGGGGGG NOT UNIQUE RECORD  "<<Name_Bkg_<<"   "<<Name_Bkg2_<<endl;}

    for (unsigned int kk=0;kk<overlaps_sorted_.size();kk++){

      if (overlaps_sorted_[kk].first == Name_Bkg_){ // && high_correlation_matrix_[i].first.c_str() !=high_correlation_matrix_[j].second.c_str()){ 

	for (unsigned int ll=kk;ll<overlaps_sorted_.size();ll++){

	  if (overlaps_sorted_[ll].first == Name_Bkg2_){
	      
	    high_correlation_position_.push_back (pair < string, int >(Name_Bkg2_, ll));

	    //  cout << " Will now remove the " << Name_Bkg2_ <<'\t'<< "  as has high Cor with  " << Name_Bkg_ <<'\t'<<" Var to be kept.. " <<'\t'<< overlaps_sorted_[kk].first <<'\t'<< " with  " <<'\t'<< overlaps_sorted_[kk].second <<" var to be removed..."<< '\t'<< overlaps_sorted_[ll].first <<'\t'<< " with  " <<'\t'<< overlaps_sorted_[ll].second <<" equal  "<<Name_Bkg2_<<endl;
	  

	  }
	}
      }
	
      else if (overlaps_sorted_[kk].first == Name_Bkg2_){ // && high_correlation_matrix_[i].first.c_str() !=high_correlation_matrix_[j].second.c_str()){ 

	for (unsigned int ll=kk;ll<overlaps_sorted_.size();ll++){

	  if (overlaps_sorted_[ll].first == Name_Bkg_){
	      
	    high_correlation_position_.push_back (pair < string, int >(Name_Bkg_, ll));
	    // cout << " Will now remove the " << Name_Bkg_ <<'\t'<< "  as has high Cor with  " << Name_Bkg2_ <<'\t'<<" Var to be kept.. " <<'\t'<< overlaps_sorted_[kk].first <<'\t'<< " with  " <<'\t'<< overlaps_sorted_[kk].second <<" var to be removed..."<< '\t'<< overlaps_sorted_[ll].first <<'\t'<< " with  " <<'\t'<< overlaps_sorted_[ll].second <<" equal  "<<Name_Bkg2_<<endl;

	  }
	}
	
	
      }      
    }
  }
  sort (high_correlation_position_.begin (), high_correlation_position_.end (), sort_string_int ());

  
  high_correlation_position_.erase (unique (high_correlation_position_.begin (), high_correlation_position_.end ()), high_correlation_position_.end ());

  overlaps_sorted_pro_cleaned_ = overlaps_sorted_;


  //   for (unsigned int l = 0; l < overlaps_sorted_.size (); l++) {
  //cout<<" Overlaps  "<<overlaps_sorted_[l].first<<"   "<<overlaps_sorted_pro_cleaned_[l].first<<"   "<<overlaps_sorted_[l].second<<"   "<<overlaps_sorted_pro_cleaned_[l].second<<"   "<<endl;   }

  for (unsigned int l = 0; l < high_correlation_position_.size (); l++) {
   
	   
    for (unsigned int m = 0; m < overlaps_sorted_pro_cleaned_.size (); m++) {

      if (high_correlation_position_[l].first==overlaps_sorted_pro_cleaned_[m].first){
	cout<<"  Will have to move from the list .... "<<"   "<<high_correlation_position_[l].first<<"  "<<overlaps_sorted_pro_cleaned_[m].first<<endl;
	overlaps_sorted_pro_cleaned_.erase(overlaps_sorted_pro_cleaned_.begin() + m);}
    }

  }


  /* for (unsigned int m = 0; m < overlaps_sorted_pro_cleaned_.size (); m++) {

  cout<<"  should be ok now... "<<overlaps_sorted_pro_cleaned_[m].first<<endl;
  }
  */


  char names[100];
  char namet[100];
  char nn[100];
	  
  for (unsigned int j = 0; j < overlaps_sorted_pro_cleaned_.size (); j++) {

    sprintf (names, "%s", overlaps_sorted_pro_cleaned_[j].first.c_str());

    for (unsigned int jj = 0; jj < overlaps_sorted_pro_cleaned_.size(); jj++) {

      sprintf (namet, "%s", overlaps_sorted_pro_cleaned_[jj].first.c_str());

      //tree_Signal->Draw (Form ("%s:%s>>hcorr", overlaps_sorted_pro_cleaned_[j].first.c_str(), overlaps_sorted_pro_cleaned_[jj].first.c_str(), "", "goff"));
      //tree_Signal->Draw (Form ("%s:%s>>hcorr", names,namet, "", "goff"));
      //TH2 *hcorrs = (TH2 *) gDirectory->Get ("hcorr");
      //float hCorr=0 ; hCorr= hcorrs->GetCorrelationFactor ();


       //tree_Signal->Draw (Form ("%s:%s>>hcorr", Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName (), "", "goff"));
        tree_Signal->Draw (Form ("%s.Draw():%s.Draw()>>hcorrt_%i%i", names, namet,j,jj),"","goff");
     
    //  sprintf(cors,"%s:%s>>%s",Branch_Sgnl_->GetName (), Branch_Sgnl2_->GetName ());

      //tree_Signal->Draw(cors);
      sprintf(nn,"hcorrt_%i%i",j,jj);
      //TH2 *hcorr_ = (TH2 *) gDirectory->Get ("hcorr");
      TH2 *hcorr_ = (TH2 *) gDirectory->Get (nn);
      
      float hCorrs = -100;hCorrs= hcorr_->GetCorrelationFactor ();

      correlation_sorted_cleaned_.push_back (hCorrs);

      if (hCorrs > cor_cut && names!=namet && j!=jj ){ cout<<" ---------------->>>>>>>>>>>>>>>>  THE Correlance Factor is "<<hCorrs<< " and should be lower that that..... CHECK!!!!!!!!!!!!!!!! " <<names<<"   and   "<<namet<<"   "<<correlation_sorted_cleaned_.size()<<endl; 
      }
      delete hcorr_;
    }
  }



  /*  for (unsigned int i = 0; i < correlation_sorted_.size (); i++) {
      correlations_sorted_.push_back (correlation_sorted_);


      }

  */
  for (unsigned int i = 0; i < correlation_sorted_cleaned_.size (); i++) {
    correlations_sorted_cleaned_.push_back (correlation_sorted_cleaned_);


  }


  overlaps_sorted_cleaned_.clear();
  overlaps_sorted_cleaned_ = overlaps_sorted_pro_cleaned_;

  if (cout_table_corr_srt_cln){

    cout << "" << endl;
    cout << "  NOW WILL PRINT SORTED AND CLEANED TABLE  FOR A CORRELATION OF " <<cor_cut<< endl;
    cout << "" << endl;





    for (unsigned int ii = 0; ii < overlaps_sorted_cleaned_.size (); ii++) {


      cout << overlaps_sorted_cleaned_[ii].first << "  " << "&\t";
    }
    cout << "\\\\" << endl;

    int ovv = overlaps_sorted_cleaned_.size ();
    for (unsigned int i = 0; i < overlaps_sorted_cleaned_.size (); i++) {
      cout << overlaps_sorted_cleaned_[i].first << "&\t";


      for (unsigned int k = 0; k < overlaps_sorted_cleaned_.size (); k++) {

	cout << correlation_sorted_cleaned_[k + ovv * i] << "&\t";
      }

      cout << " \\\\" << endl;
    }
    cout << "\\hline" << endl;

  }
  ////

}




float
GetRankedVariableByName (const string label)
{


  //overlaps_sorted_cleaned_name_

  //cout <<" size is " << overlaps_sorted_cleaned_.size ()<<endl;
  for (unsigned int i = 0; i < overlaps_sorted_cleaned_.size (); i++) {
    if (label == overlaps_sorted_cleaned_[i].first) {
      return overlaps_sorted_cleaned_[i].second;
    }
  }
  return -66666.;
}



float
GetRankedVariableByPos (const int position)
{

  return overlaps_sorted_cleaned_[position - 1].second;


  //return -66666.;
}


string GetRankedVariableName (const int position)
{

  return overlaps_sorted_cleaned_[position - 1].first;


  //return -66666.;
}



int
GetRankingofVariable (string label)
{
  for (unsigned int i = 0; i < overlaps_sorted_cleaned_.size (); i++) {
    if (label == overlaps_sorted_cleaned_[i].first) {

      return i + 1;
    }

  }
  return -9999;
}



float
GetCorrelation (string var1_, string var2_)
{




  for (unsigned int i = 0; i < overlaps_sorted_cleaned_name_.size (); i++) {

    int ovv = overlaps_sorted_cleaned_name_.size ();

    if (var1_ == overlaps_sorted_cleaned_name_[i])
      {

	for (unsigned int k = 0; k < overlaps_sorted_cleaned_name_.size (); k++) {


	  if (var2_ == overlaps_sorted_cleaned_name_[k]) {

	    // cout<<"  ADASDASSD "<<endl;
	    return correlation_sorted_cleaned_[k + ovv * i];
	  }
	}
      }
  }

  return -9999;
}
