//----------------------Version 1.0-------------------------//
//Plotting Macro for Mu -> Tau Fake Rates study
//Author: Yiwen Wen
//DESY
//-----------------------------------------------------------//
#include "PlotShapes.C"

void PlotEveryShapes(TString eta = "Gt1p558", TString wp = "VTight",int nbins =11, float xmin = 60,float xmax = 120)
{

    PlotShapes("m_vis","m_{vis} [GeV]","dN/dm_{vis}[GeV^{-1}]",eta,wp,xmin,xmax,nbins,false,true,false,false);
    PlotShapes("m_vis","m_{vis} [GeV]","dN/dm_{vis}[GeV^{-1}]",eta,wp,xmin,xmax,1,false,false,false,false);

    PlotShapes("m_vis","m_{vis} [GeV]","dN/dm_{vis}[GeV^{-1}]",eta,wp,xmin,xmax,nbins,true,true,false,false);
    //PlotShapes("m_vis","m_{vis} [GeV]","Events / bin",eta,wp,60,120,1,true,false,false,false);
    
    //float fakeratePostfit = fakeratePrefit*SF;
    //float errfakeratePostfit = errfakeratePrefit*SF;
    
    //cout << "fake rate post-fit:" << fakeratePostfit << endl;
    //cout <<"fake rate post-fit error:" << errfakeratePostfit << endl;

}

