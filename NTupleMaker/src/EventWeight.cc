#include "DesyTauAnalyses/NTupleMaker/interface/EventWeight.h"
#include "DesyTauAnalyses/NTupleMaker/interface/HTauTauUtils.h"

Double_t FitFuncGauss(Double_t * x, Double_t * par) {

  Double_t aLeft  = (x[0]-par[1])/par[2];
  Double_t bLeft  = (x[0]-par[1])/par[3];
  Double_t aRight = (x[0]-par[1])/par[4];
  Double_t bRight = (x[0]-par[1])/par[5];
  
  Double_t result = 1.0;
  
  if (x[0]<par[1]) 
    result = par[6]*TMath::Exp(-0.5*aLeft*aLeft)+(1-par[6])*TMath::Exp(-0.5*bLeft*bLeft);
  else 
    result = par[7]*TMath::Exp(-0.5*aRight*aRight)+(1-par[7])*TMath::Exp(-0.5*bRight*bRight);

  
  return par[0]*result;

}

TF1 * EventWeight::getFuncRecoil(TF1 * initFunc, bool left) {

  double xminD;
  double xmaxD;

  initFunc->GetRange(xminD,xmaxD);

  float xmin = float(xminD);
  float xmax = float(xmaxD);

  TF1 * func = new TF1("func",FitFuncGauss,xmin,xmax,8);

  func->SetParameter(0,initFunc->GetParameter(0));
  func->SetParameter(1,initFunc->GetParameter(1));
  if (left) {
    func->SetParameter(2,initFunc->GetParameter(2));
    func->SetParameter(3,initFunc->GetParameter(3));
    func->SetParameter(4,initFunc->GetParameter(2));
    func->SetParameter(5,initFunc->GetParameter(3));
    func->SetParameter(6,initFunc->GetParameter(6));
    func->SetParameter(7,initFunc->GetParameter(6));
  }
  else {
    func->SetParameter(2,initFunc->GetParameter(4));
    func->SetParameter(3,initFunc->GetParameter(5));
    func->SetParameter(4,initFunc->GetParameter(4));
    func->SetParameter(5,initFunc->GetParameter(5));
    func->SetParameter(6,initFunc->GetParameter(7));
    func->SetParameter(7,initFunc->GetParameter(7));
  }

  return func;

}


EventWeight::EventWeight(TString baseDir) {

  _baseDir = baseDir;
  _epsrel = 1e-3;
  _epsabs = 1e-3;


}

EventWeight::~EventWeight() {

}

void EventWeight::InitMEtWeights(TString fileName,
				 TString  _perpZStr,
				 TString  _paralZStr,
				 int nZPtBins,
				 float * ZPtBins,
				 TString * _ZPtStr,
				 int nJetsBins,
				 TString * _nJetsStr) {

	std::vector<float> newZPtBins;
	std::vector<std::string> newZPtStr;
	std::vector<std::string> newNJetsStr;

	std::string newPerpZStr  = std::string(_perpZStr);
	std::string newParalZStr = std::string(_paralZStr);

	for (int idx=0; idx<nZPtBins+1; ++idx) 
	  newZPtBins.push_back(ZPtBins[idx]);
	for (int idx=0; idx<nZPtBins; ++idx ) 
	  newZPtStr.push_back(std::string(_ZPtStr[idx]));

	for (int idx=0; idx<nJetsBins; ++idx)
	  newNJetsStr.push_back(std::string(_nJetsStr[idx]));

	InitMEtWeights(fileName,
		       newZPtBins,
		       newPerpZStr,
		       newParalZStr,
		       newZPtStr,
		       newNJetsStr);
}


void EventWeight::InitMEtWeights(TString fileName,
				 const std::vector<float>& ZPtBins,
				 const std::string _perpZStr,
				 const std::string _paralZStr,
				 const std::vector<std::string>& _ZPtStr,
				 const std::vector<std::string>& _nJetsStr)
{

  TFile * _fileMet = new TFile(_baseDir+"/"+fileName);

  // checking files
  if (_fileMet->IsZombie()) {
    std::cout << "File " << fileName << " is not found in directory " << _baseDir << std::endl;
    std::cout << "quitting program..." << std::endl;
    exit(-1);
  }

  _nZPtBins = ZPtBins.size()-1; // the -1 is on purpose!
  _nJetsBins = _nJetsStr.size();
  
  _ZPtBins = ZPtBins;

  for (int ZPtBin=0; ZPtBin<_nZPtBins; ++ZPtBin) {
    for (int jetBin=0; jetBin<_nJetsBins; ++jetBin) {

      TString binStrPerpData  = _perpZStr  + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_data";
      TString binStrParalData = _paralZStr + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_data";
      TString binStrPerpMC    = _perpZStr  + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc";
      TString binStrParalMC   = _paralZStr + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc";

      _metZParalData[ZPtBin][jetBin] = (TF1*)_fileMet->Get(binStrParalData);
      _metZPerpData[ZPtBin][jetBin]  = (TF1*)_fileMet->Get(binStrPerpData);
      _metZParalMC[ZPtBin][jetBin]   = (TF1*)_fileMet->Get(binStrParalMC);
      _metZPerpMC[ZPtBin][jetBin]    = (TF1*)_fileMet->Get(binStrPerpMC);


      // checking functions
      if (_metZParalData[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrParalData
		  << " is not found in file " << fileName << "... quitting program..." << std::endl;
	exit(-1);

      }
      if (_metZPerpData[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrPerpData
		  << " is not found in file " << fileName << "... quitting program..." << std::endl;
	exit(-1);
	
      }

      if (_metZParalMC[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrParalMC
		  << " is not found in file " << fileName << "... quitting program..." << std::endl;
	exit(-1);

      }
      if (_metZPerpMC[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrPerpMC
		  << " is not found in file " << fileName << "... quitting program..." << std::endl;
	exit(-1);
	
      }

      // Met Paral Data

      std::cout << _ZPtStr[ZPtBin] << " : " << _nJetsStr[jetBin] << std::endl;
      
      double xminD,xmaxD;

      _metZParalData[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZParalData[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZParalData[ZPtBin][jetBin] = float(xmaxD);

      TF1 * func = getFuncRecoil(_metZParalData[ZPtBin][jetBin],true);
      _rmsLeftMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZParalData[ZPtBin][jetBin],false);
      _rmsRightMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));
      delete func;
      
      _meanMetZParalData[ZPtBin][jetBin] = _metZParalData[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(_metZParalData[ZPtBin][jetBin]->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));

      float integral =  _metZParalData[ZPtBin][jetBin]->IntegralOneDim(_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin],_epsrel,_epsabs,_error);
      std::cout << "   Data Paral ----> " << std::endl;
      std::cout << "   X0 = " << _meanMetZParalData[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZParalData[ZPtBin][jetBin] << std::endl; 
      std::cout << "   Integral [" << _xminMetZParalData[ZPtBin][jetBin] << "," << _xmaxMetZParalData[ZPtBin][jetBin] << "] = " << integral << std::endl;

      // Met Perp Data

      _metZPerpData[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZPerpData[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZPerpData[ZPtBin][jetBin] = float(xmaxD);


      func = getFuncRecoil(_metZPerpData[ZPtBin][jetBin],true);
      _rmsLeftMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZPerpData[ZPtBin][jetBin],false);
      _rmsRightMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));
      delete func;
      
      _meanMetZPerpData[ZPtBin][jetBin] = _metZPerpData[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(_metZPerpData[ZPtBin][jetBin]->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));

      integral =  _metZPerpData[ZPtBin][jetBin]->IntegralOneDim(_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin],_epsrel,_epsabs,_error);
      std::cout << "   Data Perp -----> " << std::endl;
      std::cout << "   X0 = " << _meanMetZPerpData[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZPerpData[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZPerpData[ZPtBin][jetBin] << "," << _xmaxMetZPerpData[ZPtBin][jetBin] << "] = " << integral << std::endl;

     
      // Met Paral MC
      
      _metZParalMC[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZParalMC[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZParalMC[ZPtBin][jetBin] = float(xmaxD);

      func = getFuncRecoil(_metZParalMC[ZPtBin][jetBin],true);
      _rmsLeftMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZParalMC[ZPtBin][jetBin],false);
      _rmsRightMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));
      delete func;

      _meanMetZParalMC[ZPtBin][jetBin] = _metZParalMC[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(_metZParalMC[ZPtBin][jetBin]->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));

      integral =  _metZParalMC[ZPtBin][jetBin]->IntegralOneDim(_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin],_epsrel,_epsabs,_error);
      std::cout << "   MC Paral ------> " << std::endl;
      std::cout << "   X0 = " << _meanMetZParalMC[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZParalMC[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZParalMC[ZPtBin][jetBin] << "," << _xmaxMetZParalMC[ZPtBin][jetBin] << "] = " << integral << std::endl;
     

      // Met Perp MC

      _metZPerpMC[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZPerpMC[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZPerpMC[ZPtBin][jetBin] = float(xmaxD);

      func = getFuncRecoil(_metZPerpMC[ZPtBin][jetBin],true);
      _rmsLeftMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZPerpMC[ZPtBin][jetBin],false);
      _rmsRightMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));
      delete func;

      _meanMetZPerpMC[ZPtBin][jetBin] = _metZPerpMC[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(_metZPerpMC[ZPtBin][jetBin]->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));

      integral =  _metZPerpMC[ZPtBin][jetBin]->IntegralOneDim(_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin],_epsrel,_epsabs,_error);
      std::cout << "   MC Perp -------> " << std::endl;
      std::cout << "   X0 = " << _meanMetZPerpMC[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZPerpMC[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZPerpMC[ZPtBin][jetBin] << "," << _xmaxMetZPerpMC[ZPtBin][jetBin] << "] = " << integral << std::endl;

      _xminMetZParal[ZPtBin][jetBin] = TMath::Max(_xminMetZParalData[ZPtBin][jetBin],_xminMetZParalMC[ZPtBin][jetBin]);
      _xmaxMetZParal[ZPtBin][jetBin] = TMath::Min(_xmaxMetZParalData[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]);

      _xminMetZPerp[ZPtBin][jetBin] = TMath::Max(_xminMetZPerpData[ZPtBin][jetBin],_xminMetZPerpMC[ZPtBin][jetBin]);
      _xmaxMetZPerp[ZPtBin][jetBin] = TMath::Min(_xmaxMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]);
      
      std::cout << "   Perp  : min = " <<  _xminMetZPerp[ZPtBin][jetBin] << "    xmax = " <<  _xmaxMetZPerp[ZPtBin][jetBin] << std::endl;
      std::cout << "   Paral : min = " <<  _xminMetZParal[ZPtBin][jetBin] << "    xmax = " <<  _xmaxMetZParal[ZPtBin][jetBin] << std::endl;
      std::cout << std::endl;

      // specifying ranges ---->
      //       float range = 55;
      //       if (jetBin>0)
      // 	range = 50;
      
      //       _xminMetZPerp[ZPtBin][jetBin] = -range;
      //       _xmaxMetZPerp[ZPtBin][jetBin] = range;
      
      //       _xminMetZParal[ZPtBin][jetBin] = -range;
      //       _xmaxMetZParal[ZPtBin][jetBin] = range;
      

    }
  }

  std::cout << "Met functions downloaded...." << std::endl;
  //  exit(-1);

}

void EventWeight::RecoilCorrected(float & MetPx,
				  float & MetPy,
				  float genZPx, 
				  float genZPy,
				  float diLepPx,
				  float diLepPy,
				  int njets,
				  int method) {
  
  // input parameters
  // MetPx, MetPy - missing transverse momentum 
  //                ( corrected replaces uncorrected )
  // genZPx, genZPy - generated transverse momentum of Z
  // diLepPx, diLepPy - dilepton transverse momentum 
  // njets - number of jets 
  // method : 1 - corrections by sampling 
  //              ( calculations of quantiles )
  // method : 2 - corrections by width w(MC)=w(Data)/w(MC) w(Process)

  float Zpt = TMath::Sqrt(genZPx*genZPx + genZPy*genZPy);

  HTauTauUtils utils;

  float U1 = 0;
  float U2 = 0;
  float metU1 = 0;
  float metU2 = 0;

  utils.CalculateU1U2FromMet(MetPx,
			     MetPy,
			     genZPx,
			     genZPy,
			     diLepPx,
			     diLepPy,
			     U1,
			     U2,
			     metU1,
			     metU2);
  if (Zpt>1000)
    Zpt = 999;

  if (njets>=_nJetsBins)
    njets = _nJetsBins - 1;

  int ZptBin = binNumber(Zpt, _ZPtBins);

  if (method==1) {

    TF1 * metZParalData = _metZParalData[ZptBin][njets];
    TF1 * metZPerpData  = _metZPerpData[ZptBin][njets];
  
    TF1 * metZParalMC   = _metZParalMC[ZptBin][njets];
    TF1 * metZPerpMC     = _metZPerpMC[ZptBin][njets];

    if (U1>_xminMetZParal[ZptBin][njets]&&U1<_xmaxMetZParal[ZptBin][njets]) {

      int nSumProb = 1;
      double q[1];
      double sumProb[1];

      sumProb[0] = metZParalMC->IntegralOneDim(_xminMetZParalMC[ZptBin][njets],U1,_epsrel,_epsabs,_error);

      if (sumProb[0]<0) {
	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 1e-10;
      }
       if (sumProb[0]>1) {
	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 0.9999999;
      }
     
      
      metZParalData->GetQuantiles(nSumProb,q,sumProb);

      float U1reco = float(q[0]);

      if (U1reco>_xminMetZParal[ZptBin][njets]&&U1reco<_xmaxMetZParal[ZptBin][njets]) 
	U1 = U1reco;

    }

    if (U2>_xminMetZPerp[ZptBin][njets]&&U2<_xmaxMetZPerp[ZptBin][njets]) {

      int nSumProb = 1;
      double q[1];
      double sumProb[1];
      
      sumProb[0] = metZPerpMC->IntegralOneDim(_xminMetZPerpMC[ZptBin][njets],U2,_epsrel,_epsabs,_error);
      
      if (sumProb[0]<0) {
	//	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 1e-10;
      }
       if (sumProb[0]>1) {
	 //	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 0.9999999;
      }

      metZPerpData->GetQuantiles(nSumProb,q,sumProb);

      float U2reco = float(q[0]);
 
      if (U2reco>_xminMetZPerp[ZptBin][njets]&&U2reco<_xmaxMetZPerp[ZptBin][njets]) 
	U2 = U2reco;
      

    }

  }
  else 
    U1U2CorrectionsByWidth(U1,U2,ZptBin,njets);
  
  utils. CalculateMetFromU1U2(U1,U2,genZPx,genZPy,diLepPx,diLepPy,MetPx,MetPy);


}

float EventWeight::CorrectionsBySampling(float x, TF1 * funcMC, TF1 * funcData) {

  int nSumProb = 1;
  double q[1];
  double sumProb[1];
  
  double xD = double(x);

  double xminD = 0;
  double xmaxD = 0;

  funcMC->GetRange(xminD,xmaxD);

  float xmin = float(xminD);

  sumProb[0] = funcMC->IntegralOneDim(xmin,xD,_epsrel,_epsabs,_error);
  
  funcData->GetQuantiles(nSumProb,q,sumProb);

  float output = float(q[0]);

  return output;

}

void EventWeight::U1U2CorrectionsByWidth(float & U1, 
					 float & U2,
					 int ZptBin,
					 int njets) {

  if (njets>=_nJetsBins)
    njets = _nJetsBins - 1;
    
  // ********* U1 *************

  float width = U1 - _meanMetZParalMC[ZptBin][njets];

  if (width<0)
    width *= _rmsLeftMetZParalData[ZptBin][njets]/_rmsLeftMetZParalMC[ZptBin][njets];
  else
    width *= _rmsRightMetZParalData[ZptBin][njets]/_rmsRightMetZParalMC[ZptBin][njets];

  U1 = _meanMetZParalData[ZptBin][njets] + width;

  // ********* U2 *************

  width = U2 - _meanMetZPerpMC[ZptBin][njets];


  if (width<0)
    width *= _rmsLeftMetZPerpData[ZptBin][njets]/_rmsLeftMetZPerpMC[ZptBin][njets];
  else 
    width *= _rmsRightMetZPerpData[ZptBin][njets]/_rmsRightMetZPerpMC[ZptBin][njets];

  U2 = _meanMetZPerpData[ZptBin][njets] + width;

}

