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
 0.000108643,
                0.000388957,
                0.000332882,
                0.00038397,
                0.000549167,
                0.00105412,
                0.00459007,
                0.0210314,
                0.0573688,
                0.103986,
                0.142369,
                0.157729,
                0.147685,
                0.121027,
                0.08855,
                0.0582866,
                0.0348526,
                0.019457,
                0.0107907,
                0.00654313,
                0.00463195,
                0.00370927,
                0.0031137,
                0.00261141,
                0.00215499,
                0.00174491,
                0.00138268,
                0.00106731,
                0.000798828,
                0.00057785,
                0.00040336,
                0.00027161,
                0.000176535,
                0.00011092,
                6.75502e-05,
                4.00323e-05,
                2.32123e-05,
                1.32585e-05,
                7.51611e-06,
                4.25902e-06,
                2.42513e-06,
                1.39077e-06,
                8.02452e-07,
                4.64159e-07,
                2.67845e-07,
                1.5344e-07,
                8.68966e-08,
                4.84931e-08,
                2.6606e-08,
                1.433e-08};



double Data2012[60]={
 147883,
666514,
886594,
1.31856e+06,
2.24527e+06,
5.14763e+06,
1.67634e+07,
6.8125e+07,
2.00021e+08,
3.57497e+08,
4.51859e+08,
4.60518e+08,
3.93021e+08,
2.75288e+08,
1.56235e+08,
7.26515e+07,
2.93166e+07,
1.19476e+07,
5.79312e+06,
3.0476e+06,
1.40504e+06,
512249,
146792,
35303.4,
8271.33,
2235.61,
721.38,
258.85,
97.2709,
36.8716,
13.7275,
4.93171,
1.69241,
0.551894,
0.170601,
0.0499358,
0.0138339,
0.00362662,
0.000899606,
0.000211148,
4.68928e-05,
9.85392e-06,
1.95929e-06,
3.6862e-07,
6.56248e-08,
1.10534e-08,
1.76248e-09,
2.61497e-10,
4.768e-11,
};


double Data2012Up[60]={
 124589,
611524,
813988,
1.17491e+06,
1.8654e+06,
3.78638e+06,
1.06087e+07,
3.92691e+07,
1.31161e+08,
2.76067e+08,
3.96212e+08,
4.43964e+08,
4.2086e+08,
3.3767e+08,
2.24967e+08,
1.23818e+08,
5.73059e+07,
2.38044e+07,
1.02984e+07,
5.2607e+06,
2.8523e+06,
1.3637e+06,
526640,
163061,
42592.2,
10584.7,
2923.75,
955.685,
351.3,
136.772,
54.3136,
21.4253,
8.23874,
3.05183,
1.08184,
0.36577,
0.117761,
0.0360763,
0.0105133,
0.00291401,
0.000768174,
0.000192592,
4.59223e-05,
1.04141e-05,
2.2461e-06,
4.60737e-07,
8.98872e-08,
1.66763e-08,
2.94239e-09,
4.96039e-10,
9.24658e-11,
1.00203e-12,
};

double Data2012Down[60]={
 173216,
727331,
978923,
1.50592e+06,
2.793e+06,
7.41786e+06,
2.86966e+07,
1.1731e+08,
2.87814e+08,
4.37113e+08,
4.90842e+08,
4.50019e+08,
3.3499e+08,
1.98131e+08,
9.33802e+07,
3.68064e+07,
1.4096e+07,
6.43336e+06,
3.27032e+06,
1.44844e+06,
495063,
130080,
28609.7,
6333.36,
1680.79,
533.046,
185.383,
66.6646,
23.8662,
8.28654,
2.74673,
0.862062,
0.255191,
0.0711318,
0.0186567,
0.00460324,
0.00106833,
0.000233212,
4.78844e-05,
9.24779e-06,
1.67991e-06,
2.8704e-07,
4.6131e-08,
6.97509e-09,
9.94095e-10,
1.33153e-10,
1.59604e-11,
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