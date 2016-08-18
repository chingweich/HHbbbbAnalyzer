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
2089.61,
159956,
614014,
1.24534e+06,
1.9424e+06,
2.63405e+06,
3.53975e+06,
8.46087e+06,
2.92458e+07,
6.98017e+07,
1.36508e+08,
2.29961e+08,
3.52183e+08,
4.83336e+08,
6.04198e+08,
7.18704e+08,
8.22509e+08,
9.02549e+08,
9.52373e+08,
9.71555e+08,
9.61587e+08,
9.23768e+08,
8.62489e+08,
7.84529e+08,
6.94739e+08,
5.95764e+08,
4.91364e+08,
3.88146e+08,
2.9356e+08,
2.12821e+08,
1.47878e+08,
9.82612e+07,
6.22552e+07,
3.75496e+07,
2.15681e+07,
1.18113e+07,
6.17454e+06,
3.08639e+06,
1.47928e+06,
682895,
305781,
134410,
59350.3,
27539.8,
14454.7,
9219.55,
7180.32,
6406.24,
6116.86,
6002.58,
5940.78,
5877.42,
0,



};


double Data2012Up[60]={
1291.27,
131956,
540114,
1.04264e+06,
1.79842e+06,
2.29331e+06,
3.10406e+06,
5.3571e+06,
1.83328e+07,
4.77505e+07,
9.78151e+07,
1.71086e+08,
2.69498e+08,
3.88238e+08,
5.05185e+08,
6.12933e+08,
7.15979e+08,
8.06414e+08,
8.74111e+08,
9.14958e+08,
9.29211e+08,
9.18092e+08,
8.82715e+08,
8.26869e+08,
7.56284e+08,
6.75156e+08,
5.85725e+08,
4.90804e+08,
3.95536e+08,
3.06272e+08,
2.28088e+08,
1.63446e+08,
1.12544e+08,
7.42736e+07,
4.68804e+07,
2.82877e+07,
1.63301e+07,
9.02883e+06,
4.78652e+06,
2.43724e+06,
1.19533e+06,
567137,
262074,
119311,
54682.9,
26287,
14145.2,
9083.61,
7024.56,
6206.71,
5887.44,
5759.36,
0,



};

double Data2012Down[60]={
3325.05,
196217,
694602,
1.48527e+06,
2.12398e+06,
3.0174e+06,
4.36373e+06,
1.45022e+07,
4.53251e+07,
1.01841e+08,
1.88995e+08,
3.0916e+08,
4.52808e+08,
5.90239e+08,
7.17781e+08,
8.36487e+08,
9.31253e+08,
9.92161e+08,
1.0176e+09,
1.0093e+09,
9.68818e+08,
9.01277e+08,
8.14753e+08,
7.14885e+08,
6.04881e+08,
4.89828e+08,
3.78163e+08,
2.78431e+08,
1.95725e+08,
1.31226e+08,
8.36508e+07,
5.05525e+07,
2.89435e+07,
1.57162e+07,
8.10478e+06,
3.97605e+06,
1.86065e+06,
834398,
361206,
152850,
64797.9,
28921.3,
14776.4,
9366.98,
7358.13,
6632.61,
6373.35,
6270.91,
6206.36,
6130.47,
6030.04,
5897.58,
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