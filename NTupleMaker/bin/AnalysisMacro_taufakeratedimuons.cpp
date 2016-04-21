#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"


#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "TGraphAsymmErrors.h"



int main(int argc, char * argv[]) {

	// first argument - config file 
	// second argument - filelist

	using namespace std;

	// **** configuration
	Config cfg(argv[1]);

	const bool isData = cfg.get<bool>("IsData");
	const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

	// pile up reweighting
	const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");

	//const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
	const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

	// kinematic cuts on muons
	const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
	const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
	const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
	const float etaMuonLowCut = cfg.get<float>("etaMuonLowCut");
	const double etaMuonCut     = cfg.get<double>("etaMuonCut");
	const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
	const float dzMuonCut      = cfg.get<float>("dzMuonCut");
	const float isoMuonCut     = cfg.get<float>("isoMuonCut");
	//const bool  applyTauTauSelection = cfg.get<bool>("ApplyTauTauSelection");
	//const bool  selectZToTauTauMuMu = cfg.get<bool>("SelectZToTauTauMuMu");

	// topological cuts

	// trigger
	const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
	const string muonTriggerName  = cfg.get<string>("MuonTriggerName");
	const string muonFilterName   = cfg.get<string>("MuonFilterName");
	const string muon17FilterName = cfg.get<string>("Muon17FilterName"); 
	const string muon8FilterName = cfg.get<string>("Muon8FilterName"); 
	const string singleMuonFilterName = cfg.get<string>("SingleMuonFilterName");
	const float singleMuonTriggerPtCut = cfg.get<float>("SingleMuonTriggerPtCut");
	const float singleMuonTriggerEtaCut = cfg.get<float>("SingleMuonTriggerEtaCut");


	TString MuonTriggerName(muonTriggerName);
	TString MuonFilterName(muonFilterName);

	TString Muon17FilterName(muon17FilterName);
	TString Muon8FilterName(muon8FilterName);
	TString SingleMuonFilterName(singleMuonFilterName);

	const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
	const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");

	const double etaJetCut   = cfg.get<double>("etaJetCut");
	const double ptJetCut   = cfg.get<double>("ptJetCut");


	// topological cuts
	const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
	const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
	const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
	const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
	const double dZetaCut       = cfg.get<double>("dZetaCut");
	const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
	const bool oppositeSign    = cfg.get<bool>("oppositeSign");
	const bool isIsoR03 = cfg.get<bool>("IsIsoR03");

	// tau
	const double taupt    = cfg.get<double>("taupt");
	const double taueta    = cfg.get<double>("taueta");
	const double decayModeFinding    = cfg.get<double>("decayModeFinding");
	const double   decayModeFindingNewDMs  = cfg.get<double>("decayModeFindingNewDMs");
	const double   againstElectronVLooseMVA5  = cfg.get<double>("againstElectronVLooseMVA5");
	const double   againstMuonTight3  = cfg.get<double>("againstMuonTight3");
	const double   vertexz =  cfg.get<double>("vertexz");
	const double   byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");




	// vertex cuts
	const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
	const float zVertexCut     = cfg.get<float>("ZVertexCut");
	const float dVertexCut     = cfg.get<float>("DVertexCut");

	// jet related cuts
	//const float jetEtaCut      = cfg.get<float>("JetEtaCut");
	//const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");
	//const float jetPtHighCut   = cfg.get<float>("JetPtHighCut");
	//const float jetPtLowCut    = cfg.get<float>("JetPtLowCut");
	//const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");

	// Run range
	const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
	const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");

	//
	const string dataBaseDir = cfg.get<string>("DataBaseDir");

	// vertex distributions filenames and histname
	const string vertDataFileName = cfg.get<string>("VertexDataFileName");
	const string vertMcFileName   = cfg.get<string>("VertexMcFileName");
	const string vertHistName     = cfg.get<string>("VertexHistName");

	// lepton scale factors
	const string muonSfDataBarrel = cfg.get<string>("MuonSfDataBarrel");
	const string muonSfDataEndcap = cfg.get<string>("MuonSfDataEndcap");
	const string muonSfMcBarrel = cfg.get<string>("MuonSfMcBarrel");
	const string muonSfMcEndcap = cfg.get<string>("MuonSfMcEndcap");

	const string jsonFile = cfg.get<string>("jsonFile");
	string TrigLeg  ;
	if (!isData) TrigLeg  = cfg.get<string>("Mu17LegMC");
	if (isData) TrigLeg  = cfg.get<string>("Mu18LegData");
  	const string Region  = cfg.get<string>("Region");
  	const string Sign  = cfg.get<string>("Sign");
	
	TString MainTrigger(TrigLeg);
	// **** end of configuration

	string cmsswBase = (getenv ("CMSSW_BASE"));
	string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

	// Run-lumi selector
	std::vector<Period> periods;  
	if (isData) { // read the good runs 
		std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
		if (inputFileStream.fail() ) {
			std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
			std::cout << "please check" << std::endl;
			std::cout << "quitting program" << std::endl;
			exit(-1);
		}

		for(std::string s; std::getline(inputFileStream, s); ) {
			periods.push_back(Period());
			std::stringstream ss(s);
			ss >> periods.back();
		}
	}
	char ff[100];

	sprintf(ff,"%s/%s",argv[3],argv[2]);
	// file name and tree name
	std::string rootFileName(argv[2]);
	std::ifstream fileList(ff);
	std::ifstream fileList0(ff);
	std::string ntupleName("makeroottree/AC1B");

	TString era=argv[3];
  	TString TStrName(rootFileName+"_"+Region+"_"+Sign);
	std::cout <<TStrName <<std::endl;  
  	datasetName = rootFileName.c_str();

	TFile * file ;//= new TFile(era+"/"+TStrName+TString(".root"),"update");
  	if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  	if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");
	file->cd("");
	/*
	// file name and tree name
	std::string rootFileName(argv[2]);
	std::ifstream fileList(argv[2]);
	std::ifstream fileList0(argv[2]);
	std::string ntupleName("makeroottree/AC1B");

	TString TStrName(rootFileName);
	TString era=argv[3];
	std::cout <<TStrName <<std::endl;  

	TFile * file = new TFile(TStrName+TString(".root"),"recreate");
	file->cd("");
	*/
	std::string initNtupleName("initroottree/AC1B");
	TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
	TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
	TH1D * histWeightsSkimmedH = new TH1D("histWeightsSkimmedH","",1,-0.5,0.5);

	// Histograms after selecting unique dimuon pair
	TH1D * massSelH = new TH1D("massSelH","",200,0,200);
	TH1D * metSelH  = new TH1D("metSelH","",200,0,400);

	TH1D * hDiJetmet = new TH1D("hDiJetmet","",200,0,400);
	TH1D * hDiJetmass = new TH1D("hDiJetmass","",500,0,1000);
	TH1D * hHT_ = new TH1D("hHT_","",800,0,1600);
	TH1D * metAll  = new TH1D("metAll","",200,0,400);

	TH1D * muon1PtH = new TH1D("muon1PtH","",200,0,400);
	TH1D * muon2PtH = new TH1D("muon2PtH","",200,0,400);

	TH1F * NumberOfVerticesH = new TH1F("NumberOfVerticesH","",51,-0.5,50.5);


	//*****  create eta histogram with eta ranges associated to their names (eg. endcap, barrel)   ***** //

	TFile * fileDataNVert = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertDataFileName);
  	TFile * fileMcNVert   = new TFile(TString(cmsswBase)+"/src/"+dataBaseDir+"/"+vertMcFileName);


	TH1D * vertexDataH = (TH1D*)fileDataNVert->Get(TString(vertHistName));
	TH1D * vertexMcH   = (TH1D*)fileMcNVert->Get(TString(vertHistName));

	float normVertexData = vertexDataH->GetSumOfWeights();
	float normVertexMc   = vertexMcH->GetSumOfWeights();

	vertexDataH->Scale(1/normVertexData);
	vertexMcH->Scale(1/normVertexMc);


	PileUp * PUofficial = new PileUp();

	TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/Data_Pileup_2015D_Nov17.root","read");
	TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/MC_Spring15_PU25_Startup.root", "read");
	TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
	TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
	PUofficial->set_h_data(PU_data);
	PUofficial->set_h_MC(PU_mc);


	TFile *f10= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataBarrel);  // mu SF barrel data
	TFile *f11 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfDataEndcap); // mu SF endcap data
	TFile *f12= new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcBarrel);  // mu SF barrel MC
	TFile *f13 = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/"+muonSfMcEndcap); // mu SF endcap MC 

	TGraphAsymmErrors *hEffBarrelData = (TGraphAsymmErrors*)f10->Get("ZMassBarrel");
	TGraphAsymmErrors *hEffEndcapData = (TGraphAsymmErrors*)f11->Get("ZMassEndcap");
	TGraphAsymmErrors *hEffBarrelMC = (TGraphAsymmErrors*)f12->Get("ZMassBarrel");
	TGraphAsymmErrors *hEffEndcapMC = (TGraphAsymmErrors*)f13->Get("ZMassEndcap");

	double * dataEffBarrel = new double[10];
	double * dataEffEndcap = new double[10];
	double * mcEffBarrel = new double[10];
	double * mcEffEndcap = new double[10];

	dataEffBarrel = hEffBarrelData->GetY();
	dataEffEndcap = hEffEndcapData->GetY();
	mcEffBarrel = hEffBarrelMC->GetY();
	mcEffEndcap = hEffEndcapMC->GetY();
	double Weight=0;
	int nTotalFiles = 0;
	int nominator=-1;int denominator = -1;


	file->cd();
	TH1D * hMT = new TH1D("hMT","",20,0,200);
	TH1D * hMass = new TH1D("hMass","",20,0,200);
	TH1D * hRatioSum = new TH1D("hRatioSum","",10,0,1);
	TH1D * hDPhi = new TH1D("hDPhi","",70,0,3.5);

	int nPtBins = 3;
	//float ptBins[6] = {19,25,30,40,60,1000};
	float ptBins[4] = {19,25,40,1000};


/*
	TString PtBins[5] = {"Pt19to25",
		"Pt25to30",
		"Pt30to40",
		"Pt40to60",
		"PtGt60"};//,		       "Pt100to150",		       "Pt150to200",		        "PtGt200"};
*/
	TString PtBins[3] = {"Pt19to25",
		"Pt25to40",
		"PtGt40"};//,		       "Pt100to150",		       "Pt150to200",		        "PtGt200"};


//int nEtaBins = 3;
int nEtaBins = 1;
int nCuts = 4;
//float etaBins[4] = {0,0.9,1.2,2.4}; 
//float etaBins[3] = {0,1.48,2.4}; 
float etaBins[2] = {0,2.4}; 

TString EtaBins[1] = {
	//"EtaLt0p9",
	//"Eta0p9to1p2",
	//"EtaGt1p2"};
	"EtaLt2p4"};

TString Cuts[4] = {"Ratio","mT","DPhi","All"};
/////first stands for the Eta bin, second array for the cut 
TH1D * FakeRatePtIncLoose[1][4];
TH1D * FakeRatePtIncTight[2][4];


//TH1D * FakeRatePt[3][7];
//TH1D * FakeRateNV[3][7];
//TH1D * FakeRateEta[3][7];


TH1D * etaBinsH = new TH1D("etaBinsH", "etaBinsH", nEtaBins, etaBins);
etaBinsH->Draw();
etaBinsH->GetXaxis()->Set(nEtaBins, etaBins);
for (int i=0; i<nEtaBins; i++){ etaBinsH->GetXaxis()->SetBinLabel(i+1, EtaBins[i]);}


for (int iEta=0; iEta<nEtaBins; ++iEta) {
	for (int iCut=0; iCut<nCuts; ++iCut) {
		FakeRatePtIncLoose[iEta][iCut] = new TH1D("FakeRatePtIncLoose"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);
		FakeRatePtIncTight[iEta][iCut] = new TH1D("FakeRatePtIncTight"+EtaBins[iEta]+Cuts[iCut],"",nPtBins,ptBins);

	}

	//for (int iPt=0; iPt<nPtBins; ++iPt) {
	//  FakeRatePt[iEta][iPt] = new TH1D("FakeRatePt"+EtaBins[iEta]+PtBins[iPt],"",100,0,1000);
	//  FakeRateNV[iEta][iPt] = new TH1D("FakeRateNV"+EtaBins[iEta]+PtBins[iPt],"",50,0,50);
	//  FakeRateEta[iEta][iPt] = new TH1D("FakeRateEta"+EtaBins[iEta]+PtBins[iPt],"",80,-4,4);
	// }

}





int nFiles = 0;
int nEvents = 0;
int selEventsAllMuons = 0;
int selEventsIdMuons = 0;
int selEventsIsoMuons = 0;


std::string dummy;
// count number of files --->
while (fileList0 >> dummy) nTotalFiles++;

unsigned int RunMin = 9999999;
unsigned int RunMax = 0;
  ifstream ifs("xsecs");
  string line;

  while(std::getline(ifs, line)) // read one line from ifs
    {

      fact=fact2=1;
      istringstream iss(line); // access line as a stream

      // we only need the first two columns
      string dt,st1,st2;st1="stau2_1";st2="stau5_2";
      iss >> dt >> xs >> fact >> fact2;
      //datasetName = dt.c_str();
      //ifs >> dt >> xs; // no need to read further
      //cout<< " "<<dt<<"  "<<endl;
      //cout<< "For sample ========================"<<dt<<" xsecs is "<<xs<<" XSec "<<XSec<<"  "<<fact<<"  "<<fact2<<endl;
      //if (dt==argv[2]) {
      //if (std::string::npos != dt.find(argv[2])) {
      if (  dt == argv[2]) {
	XSec= xs*fact*fact2;
	cout<<" Found the correct cross section "<<xs<<" for Dataset "<<dt<<" XSec "<<XSec<<endl;
      }
      /*
	if ( argv[2] == st1) {ChiMass=100;mIntermediate=200;}
	else if (argv[2] == st2) {ChiMass=200;mIntermediate=500;}
      */
      if (isData) XSec=1.;
      ChiMass=0.0;
    }

  if (XSec<0&& !isData) {cout<<" Something probably wrong with the xsecs...please check  - the input was "<<argv[2]<<endl;return 0;}
	
  xsecs=XSec;

std::vector<unsigned int> allRuns; allRuns.clear();
vector <unsigned int> run_;
vector <unsigned int> lumi_;
vector <unsigned int> event_;
run_.clear();
lumi_.clear();
event_.clear();

std::vector<Event> EventList;
std::ifstream EventsFile;
TString file_events=argv[4]; //eventlist_csc2015.txt
//EventsFile.open("MET_filters/eventlist_"+file_events+".txt");
EventsFile.open(era+"/"+file_events+".txt");
//cout<<"  limits  int -> "<<std::numeric_limits<int>::max()<<"  long int -> "<<std::numeric_limits<long int>::max()<<"  unsigned int -> "<<std::numeric_limits<unsigned int>::max()<<endl;
cout<<" The file that will be used will be "<<era<<"/"<<file_events<<".txt"<<endl;


while (getline(EventsFile, line))
{
	std::vector<std::string> columns = split(line,':');
	run_.push_back(std::stoi(columns[0]));
	lumi_.push_back(std::stoi(columns[1]));
	event_.push_back(std::stoul(columns[2]));
	/*   Event events_;

	     events_.name     = "Test";
	     events_.run     = std::stoi(columns[0]);
	     events_.lumi   = std::stoi(columns[1]);
	     events_.eventrn = std::stof(columns[2]);

	     EventList.push_back(events_);
	     */
}

cout<<" In total there are "<<run_.size()<< " entries for "<<file_events<<" filter "<<endl;
EventsFile.close();
//----Attention----//
//if(XSec!=1) nTotalFiles=20;
//nTotalFiles=5;

if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);

cout<<" There are in total "<<nTotalFiles<<endl;
for (int iF=0; iF<nTotalFiles; ++iF) {

	std::string filen;
	fileList >> filen;

	std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
	TFile * file_ = TFile::Open(TString(filen));

	TH1D * histoInputEvents = NULL;
	histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
	if (histoInputEvents==NULL) continue;
	int NE = int(histoInputEvents->GetEntries());
	for (int iE=0;iE<NE;++iE)
		inputEventsH->Fill(0.);
	std::cout << "      number of input events         = " << NE << std::endl;

	TTree * _inittree = NULL;
	_inittree = (TTree*)file_->Get(TString(initNtupleName));
	if (_inittree==NULL) continue;
	Float_t genweight;
	if (!isData)
		_inittree->SetBranchAddress("genweight",&genweight);
	Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
	std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
	for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
		_inittree->GetEntry(iEntry);
		if (isData)
			histWeightsH->Fill(0.,1.);
		else
			histWeightsH->Fill(0.,genweight);
	}

	TTree * _tree = NULL;
	_tree = (TTree*)file_->Get(TString(ntupleName));
	if (_tree==NULL) continue;
	Long64_t numberOfEntries = _tree->GetEntries();
	std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
	AC1B analysisTree(_tree);


	// EVENT LOOP //
	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 

		analysisTree.GetEntry(iEntry);
		nEvents++;

		if (nEvents%50000==0) 
			cout << "      processed " << nEvents << " events" << endl; 

		float weight = 1;

		//------------------------------------------------

		if (!isData) 
			weight *=analysisTree.genweight;

		histWeightsSkimmedH->Fill(float(0),weight);

		if (!isData) {/*
			if (applyPUreweighting_vertices) {
				int binNvert = vertexDataH->FindBin(analysisTree.primvertex_count);
				float_t dataNvert = vertexDataH->GetBinContent(binNvert);
				float_t mcNvert = vertexMcH->GetBinContent(binNvert);
				if (mcNvert < 1e-10){mcNvert=1e-10;}
				float_t vertWeight = dataNvert/mcNvert;
				weight *= vertWeight;
				//	  cout << "NVert = " << analysisTree.primvertex_count << "   weight = " << vertWeight << endl;
			}
*/
			if (applyPUreweighting) {

				double Ninteractions = analysisTree.numtruepileupinteractions;
				double PUweight = PUofficial->get_PUweight(Ninteractions);
				weight *= PUweight;
				//PUweightsOfficialH->Fill(PUweight);

			}
/*
			if (applyTauTauSelection) {
				unsigned int nTaus = 0;
				if (analysisTree.gentau_count>0) {
					//	  cout << "Generated taus present" << endl;
					for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++itau) {
						//	    cout << itau << "  : pt = " 
						//		 << analysisTree.gentau_visible_pt[itau] 
						//		 << "   eta = " <<  analysisTree.gentau_visible_eta[itau]
						//		 << "   mother = " << int(analysisTree.gentau_mother[itau]) << endl;
						if (int(analysisTree.gentau_mother[itau])==3) nTaus++;

					}
				}
				bool notTauTau = nTaus < 2;
				//	  std::cout << "nTaus = " << nTaus << std::endl;

				if (selectZToTauTauMuMu&&notTauTau) { 
					//	    std::cout << "Skipping event..." << std::endl;
					//	    cout << endl;
					continue;
				}
				if (!selectZToTauTauMuMu&&!notTauTau) { 
					//	    std::cout << "Skipping event..." << std::endl;
					//	    cout << endl;
					continue;
				}
				//	  cout << endl;
			}*/
		}

		if (isData && applyGoodRunSelection){


			bool lumi = false;
			int n=analysisTree.event_run;
			int lum = analysisTree.event_luminosityblock;
			int nr = analysisTree.event_nr;

			std::string num = std::to_string(n);
			std::string lnum = std::to_string(lum);
			for(const auto& a : periods)
			{

				if ( num.c_str() ==  a.name ) {
					//std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
					//     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

					for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

						//	cout<<b->lower<<"  "<<b->bigger<<endl;
						if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
					}
					auto last = std::prev(a.ranges.end());
					//    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
					if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


				}

			}

			if (!lumi) continue;
			bool runbool=false;
			bool lumibool=false;
			bool eventbool=false;
			runbool= std::find(run_.begin(), run_.end(), n) != run_.end();
			lumibool= std::find(lumi_.begin(), lumi_.end(), lum) != lumi_.end();
			eventbool= std::find(event_.begin(), event_.end(), nr) != event_.end();

			if (runbool && lumibool && eventbool) continue;
			//if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
			//std::remove("myinputfile");
		}     
		float metall = sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex+analysisTree.pfmet_ey*analysisTree.pfmet_ey);
		metAll->Fill(metall,weight);

		if (analysisTree.event_run<RunMin)
			RunMin = analysisTree.event_run;

		if (analysisTree.event_run>RunMax)
			RunMax = analysisTree.event_run;

		//std::cout << " Run : " << analysisTree.event_run << std::endl;

		bool isNewRun = true;
		if (allRuns.size()>0) {
			for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
				if (analysisTree.event_run==allRuns.at(iR)) {
					isNewRun = false;
					break;
				}
			}
		}

		if (isNewRun) 
			allRuns.push_back(analysisTree.event_run);
			unsigned int nMainTrigger = 0;
			bool isMainTrigger = false;

		unsigned int nMuonFilter = 0;

			unsigned int nfilters = analysisTree.run_hltfilters->size();
			//  std::cout << "nfiltres = " << nfilters << std::endl;
			for (unsigned int i=0; i<nfilters; ++i) {
				//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
				TString HLTFilter(analysisTree.run_hltfilters->at(i));
				if (HLTFilter==MainTrigger) {
					nMainTrigger = i;
					isMainTrigger = true;
					nMuonFilter = i;
				}



			}

			if (!isMainTrigger) {
				std::cout << "Fail on main HLT Filter " << MainTrigger << " not found" << std::endl;
				return(-1);
			}



/*

		bool isTriggerMuon = false;
		for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
			TString trigName(it->first);
			if (trigName.Contains(MuonTriggerName)) {
				//	  std::cout << it->first << " : " << it->second << std::endl;
				if (it->second==1)
					isTriggerMuon = true;
			}
		}

		if (applyTrigger && !isTriggerMuon) continue;

		unsigned int nMuonFilter = 0;
		bool isMuonFilter = false;

		unsigned int nMuon17Filter = 0;
		bool isMuon17Filter = false;

		unsigned int nMuon8Filter = 0;
		bool isMuon8Filter = false;

		unsigned int nSingleMuonFilter = 0;
		bool isSingleMuonFilter = false;

		unsigned int nfilters = analysisTree.run_hltfilters->size();
		for (unsigned int i=0; i<nfilters; ++i) {
			TString HLTFilter(analysisTree.run_hltfilters->at(i));
			if (HLTFilter==MuonFilterName) {
				nMuonFilter = i;
				isMuonFilter = true;
			}
			if (HLTFilter==Muon17FilterName) {
				nMuon17Filter = i;
				isMuon17Filter = true;
			}
			if (HLTFilter==Muon8FilterName) {
				nMuon8Filter = i;
				isMuon8Filter = true;
			}
			if (HLTFilter==SingleMuonFilterName) {
				nSingleMuonFilter = i;
				isSingleMuonFilter = true;
			}
		}
		if (!isMuonFilter) {
			cout << "Filter " << MuonFilterName << " not found " << endl;
			exit(-1);
		}
		if (!isMuon17Filter) {
			cout << "Filter " << Muon17FilterName << " not found " << endl;
			exit(-1);
		}
		if (!isMuon8Filter) {
			cout << "Filter " << Muon8FilterName << " not found " << endl;
			exit(-1);
		}
		if (!isSingleMuonFilter) {
			cout << "Filter " << SingleMuonFilterName << " not found " << endl;
			exit(-1);
		}
*/
		// vertex cuts

		if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
		if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
		float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
				analysisTree.primvertex_y*analysisTree.primvertex_y);
		if (dVertex>dVertexCut) continue;

		// muon selection

		vector<unsigned int> allMuons; allMuons.clear();
		vector<unsigned int> idMuons; idMuons.clear();
		vector<unsigned int> isoMuons; isoMuons.clear();
		vector<float> isoMuonsValue; isoMuonsValue.clear();
		vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
		vector<bool> isMuonMatched23Filter; isMuonMatched23Filter.clear();
		vector<bool> isMuonMatched17Filter; isMuonMatched17Filter.clear();
		vector<bool> isMuonMatched8Filter; isMuonMatched8Filter.clear();
		vector<bool> isMuonMatchedSingleMuFilter; isMuonMatchedSingleMuFilter.clear();

		for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
			allMuons.push_back(im);
			bool muPassed    = true;
			bool mu23Matched = false;
			bool mu17Matched = false;
			bool mu8Matched  = false;
			bool muSingleMatched = false;
			if (analysisTree.muon_pt[im]<5) muPassed = false;
			if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) muPassed = false;
			if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) muPassed = false;
			if (fabs(analysisTree.muon_dz[im])>dzMuonCut) muPassed = false;
			if (!analysisTree.muon_isMedium[im]) muPassed = false;
			if (muPassed) idMuons.push_back(im);
			float absIso = analysisTree.muon_r03_sumChargedHadronPt[im];
			float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] + 
				analysisTree.muon_r03_sumPhotonEt[im] - 
				0.5*analysisTree.muon_r03_sumPUPt[im];
			neutralIso = TMath::Max(float(0),neutralIso); 
			absIso += neutralIso;
			float relIso = absIso/analysisTree.muon_pt[im];
			if (relIso>isoMuonCut) muPassed = false;
			if (muPassed && analysisTree.muon_pt[im]>ptMuonLowCut) { 
				isoMuons.push_back(im);
				isoMuonsValue.push_back(relIso);
			}

			isMuonPassedIdIso.push_back(muPassed);


			for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
				if (analysisTree.trigobject_filters[iT][nMainTrigger]) { // Mu17 Leg
					double dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
							analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
					if (dRtrig>deltaRTrigMatch) continue;
					//if (!isData && analysisTree.trigobject_filters[iT][nMainTrigger] && analysisTree.trigobject_pt[iT]>18) 
					isMainTrigger = true;

				}
			}

			if (!isMainTrigger) continue;


/*
			//	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
			for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
				float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
						analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
				if (dRtrig>DRTrigMatch) continue;
				if (analysisTree.trigobject_filters[iT][nMuon17Filter] && analysisTree.trigobject_pt[iT]>23.0) 
					mu23Matched = true;
				if (analysisTree.trigobject_filters[iT][nMuon17Filter] && analysisTree.trigobject_pt[iT]>17.0)
					mu17Matched = true;
				if (analysisTree.trigobject_filters[iT][nSingleMuonFilter] && analysisTree.trigobject_pt[iT]>singleMuonTriggerPtCut)
					muSingleMatched = true;
				if (analysisTree.trigobject_filters[iT][nMuon8Filter] && analysisTree.trigobject_pt[iT]>8.0) 
					mu8Matched = true;

			}
			isMuonMatched23Filter.push_back(mu23Matched);
			isMuonMatched17Filter.push_back(mu17Matched);
			isMuonMatched8Filter.push_back(mu8Matched);
			isMuonMatchedSingleMuFilter.push_back(muSingleMatched);*/
		}


		unsigned int indx1 = 0;
		unsigned int indx2 = 0;
		bool isIsoMuonsPair = false;
		float isoMin = 9999;
		if (isoMuons.size()>0) {
			for (unsigned int im1=0; im1<isoMuons.size(); ++im1) {
				unsigned int index1 = isoMuons[im1];
				bool isMu1matched = false;
				for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
					float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
							analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
					if (dRtrig>DRTrigMatch) continue;
					if (analysisTree.trigobject_filters[iT][nMainTrigger] && 
							analysisTree.muon_pt[index1] > ptMuonHighCut &&
							fabs(analysisTree.muon_eta[index1]) < etaMuonHighCut) 
						isMu1matched = true;
				}
				if (isMu1matched) {
					for (unsigned int iMu=0; iMu<allMuons.size(); ++iMu) {
						unsigned int indexProbe = allMuons[iMu];
						if (index1==indexProbe) continue;
						float q1 = analysisTree.muon_charge[index1];
						float q2 = analysisTree.muon_charge[indexProbe];
						if (q1*q2>0) continue;
						float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
								analysisTree.muon_eta[indexProbe],analysisTree.muon_phi[indexProbe]);
						if (dR<dRleptonsCut) continue; 
						float dPhi = dPhiFrom2P(analysisTree.muon_px[index1],analysisTree.muon_py[index1],
								analysisTree.muon_px[indexProbe],analysisTree.muon_py[indexProbe]);
						if (dPhi>dPhileptonsCut) continue;
						TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
								analysisTree.muon_py[index1],
								analysisTree.muon_pz[index1],
								muonMass);
						TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[indexProbe],
								analysisTree.muon_py[indexProbe],
								analysisTree.muon_pz[indexProbe],
								muonMass);


						float mass = (muon1+muon2).M();
					}
				}
				for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
					unsigned int index2 = isoMuons[im2];
					float q1 = analysisTree.muon_charge[index1];
					float q2 = analysisTree.muon_charge[index2];
					bool isMu2matched = false;
					for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
						float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
								analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
						if (dRtrig>DRTrigMatch) continue;
						if (analysisTree.trigobject_filters[iT][nMainTrigger] && 
								analysisTree.muon_pt[index2] > ptMuonHighCut &&
								fabs(analysisTree.muon_eta[index2]) < etaMuonHighCut) 
							isMu2matched = true;
					}
					bool isPairSelected = q1*q2 > 0;
					if (oppositeSign) isPairSelected = q1*q2 < 0;
					bool isTriggerMatch = (isMu1matched || isMu2matched);
					float dRmumu = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
							analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
					if (isTriggerMatch && isPairSelected && dRmumu>dRleptonsCut) {
						bool sumIso = isoMuonsValue[im1]+isoMuonsValue[im2];
						if (sumIso<isoMin) {
							isIsoMuonsPair = true;
							isoMin = sumIso;
							if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]) {
								indx1 = index1;
								indx2 = index2;
							}
							else {
								indx2 = index1;
								indx1 = index2;
							}
						}
					}
				}
			}
		}
		if (isIsoMuonsPair) {      
			//match to genparticles

			double isoTauMin = 999;
			bool tau_iso = false;
			bool isTight = false;
			bool isLoose = false;
			unsigned tau_loose=-1;
			vector<int> tau; tau.clear();


			for (unsigned  int it = 0; it<analysisTree.tau_count; ++it) {

				if (analysisTree.tau_pt[it] < ptMuonLowCut || fabs(analysisTree.tau_eta[it])> etaMuonCut) continue;
				if (analysisTree.tau_decayModeFindingNewDMs[it]<decayModeFindingNewDMs) continue;
				if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
				double  tauIso = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[it];

				isLoose  = true;
				tau_loose = int(it);

				if (analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits [it]> 0.5  && analysisTree.tau_againstElectronVLooseMVA5[it]>againstElectronVLooseMVA5 
						&& analysisTree.tau_againstMuonTight3[it]>againstMuonTight3) {isTight = true;	  tau_tight = int(it);}




			}


			if (!isLoose) continue;
			TLorentzVector JetsV;

			JetsMV.clear();
			int countjets=0;
			for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
				float absJetEta = fabs(analysisTree.pfjet_eta[jet]);

				if (analysisTree.pfjet_pt[jet] < 19) continue;
				if (absJetEta > etaJetCut) continue;
				//if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
				bool  looseJetID = looseJetiD(analysisTree,jet);
				if (!looseJetID) continue;	  

				double dRmuJet1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
						analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1]);

				if (dRmuJet1 < 0.5) continue;
				double dRmuJet2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
						analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);

				if (dRmuJet2 < 0.5) continue;

				JetsV.SetPxPyPzE(0.,0.,0.,0.);
				JetsV.SetPxPyPzE(analysisTree.pfjet_px[jet], analysisTree.pfjet_py[jet], analysisTree.pfjet_pz[jet], analysisTree.pfjet_e[jet]);
				JetsMV.push_back(JetsV);	
				countjets++;
			}
			sort(JetsMV.begin(), JetsMV.end(),ComparePt);

			if (countjets ==0) continue;



			double dPhi=-1;double MT=-1 ; double RatioSums=-1;


			double met = sqrt ( analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
			// w = mu+MET
			// ptW - ptJ/ptW+ptJ      
			//
			//

			TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[indx1],
					analysisTree.muon_py[indx1],
					analysisTree.muon_pz[indx1],
					muonMass);
			TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[indx2],
					analysisTree.muon_py[indx2],
					analysisTree.muon_pz[indx2],
					muonMass);
			TLorentzVector DiM = muon1+muon2;

			RatioSums = analysisTree.tau_pt[tau_loose]/DiM.Pt();


			TLorentzVector MetV; 
			MetV.SetPx(analysisTree.pfmet_ex); 
			MetV.SetPy(analysisTree.pfmet_ey);

			TLorentzVector tauV;  tauV.SetPtEtaPhiM(analysisTree.tau_pt[tau_loose], analysisTree.tau_eta[tau_loose], analysisTree.tau_phi[tau_loose], tauMass);

			TLorentzVector DiL = DiM  + tauV;

			//dPhi = dPhiFrom2P( DiM.Px(), DiM.Py(), MetV.Px(),  MetV.Py() );
			dPhi = dPhiFrom2P( DiM.Px(), DiM.Py(), analysisTree.tau_px[tau_loose],analysisTree.tau_py[tau_loose]);
			MT = TMath::Sqrt(2*DiM.Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));


			hRatioSum->Fill(RatioSums,weight);
			hMT->Fill(MT,weight);
			hMass->Fill(DiM.M(),weight);
			hDPhi->Fill(dPhi, weight);


			//if (tau.size()==0 || !tau_iso ) continue;

			//cout<<"  "<<endl;
			if (isLoose) denominator++;
			if (isTight) nominator++;

			if (isLoose){
				float ptProbe = TMath::Min(float(analysisTree.tau_pt[tau_loose]),float(ptBins[nPtBins]-0.1));
				float absEtaProbe = fabs(analysisTree.tau_eta[tau_loose]);
				int ptBin = binNumber(ptProbe,nPtBins,ptBins);
				if (ptBin<0) continue;
				int etaBin = binNumber(absEtaProbe,nEtaBins,etaBins);
				if (etaBin<0) continue;

				//cout<< "filling here  "<<analysisTree.tau_pt[tau_loose]<<"  "<<ptBin<<"  "<<etaBin<<"  "<<weight<<endl;
				//FakeRatePt[etaBin][ptBin]->Fill(analysisTree.tau_pt[tau_loose],weight);

				bool bRatio = (RatioSums < 1.2  && RatioSums > 0.8);
				bool bMass = (DiM.M() < 120 && DiM.M() > 60);
				if (isLoose){
					if (bRatio) FakeRatePtIncLoose[etaBin][0]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (bMass && bRatio) FakeRatePtIncLoose[etaBin][1]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (dPhi > 2.5  && bRatio) FakeRatePtIncLoose[etaBin][2]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (bMass  && dPhi > 2.5 && bRatio) FakeRatePtIncLoose[etaBin][3]->Fill(analysisTree.tau_pt[tau_loose],weight);

				}

				if (isTight) 

				{	
					if (bRatio) FakeRatePtIncTight[etaBin][0]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (bMass  && bRatio) FakeRatePtIncTight[etaBin][1]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (dPhi > 2.5 && bRatio) FakeRatePtIncTight[etaBin][2]->Fill(analysisTree.tau_pt[tau_loose],weight);
					if (bMass  && dPhi > 2.5 && bRatio) FakeRatePtIncTight[etaBin][3]->Fill(analysisTree.tau_pt[tau_loose],weight);

				}

				//FakeRateEta[etaBin][ptBin]->Fill(analysisTree.tau_eta[tau_loose],weight);
				//FakeRateNV[etaBin][ptBin]->Fill(analysisTree.tau_vertexz[tau_loose],weight);

			}
		}
	} // end of file processing (loop over events in one file)
	nFiles++;
	delete _tree;
	file_->Close();
	delete file_;
}
std::cout << std::endl;
int allEvents = int(inputEventsH->GetEntries());
std::cout << "Total number of input events                     = " << allEvents << std::endl;
std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
std::cout << std::endl;
std::cout << "RunMin = " << RunMin << std::endl;
std::cout << "RunMax = " << RunMax << std::endl;

//cout << "weight used:" << weight << std::endl;

// using object as comp
std::sort (allRuns.begin(), allRuns.end(), myobject);
std::cout << "Runs   :";
for (unsigned int iR=0; iR<allRuns.size(); ++iR)
std::cout << " " << allRuns.at(iR);
std::cout << std::endl;

file->cd();
hxsec->Fill(XSec);
hxsec->Write();
inputEventsH->Write();
histWeightsH->Write();
file->Write();
file->Close();
delete file;

}



