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
1782.75,
26901.5,
178266,
471412,
760778,
1.02496e+06,
1.47743e+06,
7.35145e+06,
2.29578e+07,
3.75378e+07,
6.0116e+07,
9.31793e+07,
1.40615e+08,
2.08958e+08,
2.87676e+08,
3.53149e+08,
3.93263e+08,
4.08686e+08,
3.99528e+08,
3.6893e+08,
3.23124e+08,
2.68661e+08,
2.11647e+08,
1.57058e+08,
1.08779e+08,
6.9628e+07,
4.08606e+07,
2.18857e+07,
1.06993e+07,
4.79698e+06,
1.99105e+06,
775953,
289879,
107460,
42190.1,
19509.5,
11629.1,
8731.49,
7497.94,
6849.6,
6443.46,
6164.29,
5962.51,
5806.68,
5667.29,
5528.68,
5378.08,
5206.21,
5008.45,
4783.92,
4534.42,
4263.6,
0,
0,
0,

};


double Data2012Up[60]={
1107.85,
22613.4,
152712,
375694,
685952,
916877,
1.12966e+06,
3.68767e+06,
1.55721e+07,
2.93011e+07,
4.50509e+07,
7.0282e+07,
1.0477e+08,
1.54871e+08,
2.22135e+08,
2.92434e+08,
3.46681e+08,
3.7869e+08,
3.89495e+08,
3.79114e+08,
3.50522e+08,
3.08928e+08,
2.59727e+08,
2.07986e+08,
1.57894e+08,
1.1276e+08,
7.51012e+07,
4.6301e+07,
2.62952e+07,
1.37381e+07,
6.62017e+06,
2.96151e+06,
1.24254e+06,
496097,
192638,
75547.7,
32051.7,
16219.6,
10402.1,
8119.24,
7082.61,
6511.4,
6143.72,
5886.6,
5698.46,
5553.03,
5425.45,
5300.65,
5168.3,
5019.53,
4849.38,
4656.05,
0,
0,

};

double Data2012Down[60]={
2620.96,
33782.7,
211551,
576943,
861247,
1.1334e+06,
2.48634e+06,
1.37439e+07,
3.09392e+07,
4.99028e+07,
8.12365e+07,
1.25924e+08,
1.93167e+08,
2.7913e+08,
3.57593e+08,
4.08055e+08,
4.29427e+08,
4.22108e+08,
3.89358e+08,
3.38669e+08,
2.78074e+08,
2.14996e+08,
1.55391e+08,
1.03841e+08,
6.34446e+07,
3.51459e+07,
1.75881e+07,
7.96713e+06,
3.29246e+06,
1.25811e+06,
453482,
159186,
57793,
24222.2,
13227.8,
9459.54,
7967.6,
7225.06,
6774.02,
6469.84,
6252.81,
6084.77,
5931.42,
5775.53,
5601.99,
5401.3,
5169.62,
4907.36,
4617.92,
4306.56,
3979.57,
3643.61,
0,
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