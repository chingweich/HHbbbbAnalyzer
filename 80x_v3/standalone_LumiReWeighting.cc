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
 0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
 		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
 		0.0,
 		0.0,
		0.0 
		};



double Data2012[60]={
5045.95,
240723,
783161,
1.73514e+06,
2.3678e+06,
3.41347e+06,
6.1202e+06,
2.43051e+07,
6.78253e+07,
1.44538e+08,
2.56896e+08,
4.06022e+08,
5.62421e+08,
7.05479e+08,
8.40326e+08,
9.53264e+08,
1.02925e+09,
1.0643e+09,
1.06023e+09,
1.01887e+09,
9.45857e+08,
8.50802e+08,
7.40484e+08,
6.18853e+08,
4.92509e+08,
3.71891e+08,
2.6666e+08,
1.81724e+08,
1.17459e+08,
7.17341e+07,
4.12926e+07,
2.24099e+07,
1.14838e+07,
5.56604e+06,
2.55802e+06,
1.11961e+06,
470081,
191679,
77758.7,
33015.2,
16073.9,
9871.76,
7673.3,
6917.35,
6659.5,
6558.48,
6487.19,
6398.03,
6277.76,
6119.3,
0,


};


double Data2012Up[60]={
3262.51,
194494,
691041,
1.47463e+06,
2.11508e+06,
3.00073e+06,
4.31454e+06,
1.41709e+07,
4.45189e+07,
1.00256e+08,
1.86453e+08,
3.0538e+08,
4.48272e+08,
5.85501e+08,
7.12692e+08,
8.31392e+08,
9.26756e+08,
9.88637e+08,
1.01522e+09,
1.00822e+09,
9.69127e+08,
9.02831e+08,
8.17342e+08,
7.18417e+08,
6.09255e+08,
4.94735e+08,
3.8311e+08,
2.82962e+08,
1.9957e+08,
1.34291e+08,
8.59511e+07,
5.21703e+07,
3.00064e+07,
1.63699e+07,
8.48254e+06,
4.18166e+06,
1.96629e+06,
885801,
385002,
163382,
69275,
30754,
15498.3,
9638,
7451.85,
6658.89,
6375.28,
6265.24,
6200.28,
6126.6,
0,


};

double Data2012Down[60]={
7532.07,
297915,
906101,
1.99923e+06,
2.72768e+06,
3.93039e+06,
1.02854e+07,
4.03937e+07,
1.03019e+08,
2.0642e+08,
3.53703e+08,
5.28528e+08,
6.92137e+08,
8.44773e+08,
9.78252e+08,
1.07173e+09,
1.11768e+09,
1.11765e+09,
1.07395e+09,
9.93159e+08,
8.86869e+08,
7.63139e+08,
6.26992e+08,
4.87417e+08,
3.57483e+08,
2.47705e+08,
1.62168e+08,
9.99646e+07,
5.77806e+07,
3.12758e+07,
1.5874e+07,
7.56895e+06,
3.39869e+06,
1.44339e+06,
584073,
228089,
88154.5,
35624.4,
16695.1,
10124,
7922.96,
7209.47,
6977.38,
6882.3,
6799.25,
6689.13,
6538.52,
6342.06,
6100.03,
5816.08,
0,



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