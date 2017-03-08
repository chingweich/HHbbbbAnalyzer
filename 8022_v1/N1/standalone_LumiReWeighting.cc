/**
   \class    standalone_LumiReWeighting standalone_LumiReWeighting.h "PhysicsTools/Utilities/interface/standalone_LumiReWeighting.h"
   \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data

   This class will trivially take two histograms:
   1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
   2. A histogram generated from the "estimatePileup" macro here:

   https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup

   and produce weights to convert the input distribution (1) to the latter (2).

   \author Shin-Shan Eiko Yu, Salvatore Rappoccio, modified by Mike Hildreth
  
*/
#ifndef standalone_LumiReWeighting_cxx
#define standalone_LumiReWeighting_cxx
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "standalone_LumiReWeighting.h"

//=======================================================
// For 2012 Data and MC
//=======================================================

Double_t Summer2012_S10[60] = {
  1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 
		};



double Data2012[60]={
238797,
837543,
2.30843e+06,
3.12475e+06,
4.47619e+06,
5.99591e+06,
7.0009e+06,
1.28916e+07,
3.52617e+07,
7.87011e+07,
1.76945e+08,
3.60088e+08,
6.0276e+08,
8.76504e+08,
1.17444e+09,
1.48901e+09,
1.75929e+09,
1.94385e+09,
2.04911e+09,
2.10153e+09,
2.13274e+09,
2.14905e+09,
2.12894e+09,
2.06261e+09,
1.96285e+09,
1.84185e+09,
1.70411e+09,
1.5545e+09,
1.39946e+09,
1.24349e+09,
1.08878e+09,
9.37258e+08,
7.91998e+08,
6.56675e+08,
5.34429e+08,
4.27096e+08,
3.35081e+08,
2.57707e+08,
1.93739e+08,
1.41823e+08,
1.00667e+08,
6.90111e+07,
4.55386e+07,
2.88467e+07,
1.75059e+07,
1.01625e+07,
5.63771e+06,
2.98725e+06,
1.51199e+06,
731842,
339821,
152545,


};


double Data2012Up[60]={
232683,
659470,
2.18387e+06,
2.74582e+06,
4.07155e+06,
5.39991e+06,
6.38596e+06,
9.0417e+06,
2.37826e+07,
5.40018e+07,
1.16012e+08,
2.46045e+08,
4.43373e+08,
6.80092e+08,
9.3784e+08,
1.21877e+09,
1.50164e+09,
1.7304e+09,
1.88258e+09,
1.96863e+09,
2.01262e+09,
2.04033e+09,
2.05469e+09,
2.03687e+09,
1.97828e+09,
1.88941e+09,
1.78123e+09,
1.65793e+09,
1.52312e+09,
1.38222e+09,
1.23964e+09,
1.09764e+09,
9.57667e+08,
8.2189e+08,
6.93274e+08,
5.74787e+08,
4.68605e+08,
3.75745e+08,
2.96156e+08,
2.29062e+08,
1.73372e+08,
1.27958e+08,
9.17531e+07,
6.37043e+07,
4.27027e+07,
2.75724e+07,
1.71182e+07,
1.02056e+07,
5.83745e+06,
3.2018e+06,
1.68401e+06,
849964,




};

double Data2012Down[60]={
247409,
1.06921e+06,
2.42825e+06,
3.56688e+06,
4.99168e+06,
6.59299e+06,
8.09691e+06,
1.99609e+07,
5.191e+07,
1.19767e+08,
2.72745e+08,
5.13122e+08,
8.02267e+08,
1.11867e+09,
1.46287e+09,
1.78008e+09,
2.00506e+09,
2.13496e+09,
2.19839e+09,
2.23388e+09,
2.25248e+09,
2.22965e+09,
2.15416e+09,
2.04166e+09,
1.9056e+09,
1.75109e+09,
1.5846e+09,
1.4136e+09,
1.24253e+09,
1.07363e+09,
9.0977e+08,
7.55103e+08,
6.13891e+08,
4.8913e+08,
3.81991e+08,
2.92068e+08,
2.18024e+08,
1.58233e+08,
1.11114e+08,
7.51472e+07,
4.8753e+07,
3.02455e+07,
1.79006e+07,
1.00897e+07,
5.41014e+06,
2.75826e+06,
1.33743e+06,
617795,
273138,
116942,
49859.9,
22465.6,



};



standalone_LumiReWeighting::standalone_LumiReWeighting(int mode) {

  std::cout << "=======================================================================" << std::endl;
  
  std::vector<double> MC_distr;
  std::vector<double> Lumi_distr;

  MC_distr.clear();
  Lumi_distr.clear();
  switch (mode)
    {
    case 0:
      std::cout << "Using central value " << std::endl;
      break;
    case 1:
      std::cout << "Using +1 sigma 5% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 5% value " << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
    } // end of switch

  Int_t NBins = 60;
  
  for( int i=0; i< NBins; ++i) {
    switch (mode){
    case 0:
      Lumi_distr.push_back(Data2012[i]);
      break;
    case 1:
      Lumi_distr.push_back(Data2012Up[i]);
      break;
    case -1:
      Lumi_distr.push_back(Data2012Down[i]);
      break;
    default:
      Lumi_distr.push_back(Data2012[i]);
      break;
    } // end of switch

    MC_distr.push_back(Summer2012_S10[i]);
  } // end of loop over bins

  // no histograms for input: use vectors
  
  // now, make histograms out of them:

  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr <<"ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";

  }


  weights_ = new TH1D(Form("luminumer_%d",mode),
 		      Form("luminumer_%d",mode),
 		      NBins,0.0, double(NBins));

  TH1D* den = new TH1D(Form("lumidenom_%d",mode),
 		       Form("lumidenom_%d",mode),
 		       NBins,0.0, double(NBins));


  
  for(int ibin = 1; ibin<NBins+1; ++ibin ) {
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    den->SetBinContent(ibin,MC_distr[ibin-1]);
  }
/*
  std::cout << "Data Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }
  std::cout << "MC Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << den->GetBinContent(ibin) << std::endl;
  }
*/
  // check integrals, make sure things are normalized

  double deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }
  double deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

  std::cout << "Reweighting: Computed Weights per In-Time Nint " << std::endl;

/*
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }

  std::cout << "=======================================================================" << std::endl;
*/
}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}



double standalone_LumiReWeighting::weight( double npv ) {
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
}


#endif	
