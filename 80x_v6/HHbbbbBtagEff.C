//Histogram 
#include <TH1.h>
//#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
//#include <TGraph.h>
//#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
//#include <fstream>
//#include <sstream>

//root feature
//#include <TLegend.h>
//#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
//#include <TCanvas.h>
//#include "TSystem.h"
//#include "TStyle.h"
#include <TClonesArray.h>
#include <TF1.h>
//math 
#include <cmath>
//#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "untuplizer.h"
//#include "jetEnergyScale.h"

//#include "getPUPPIweight.h"
#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"
//#include "HHbbbbBtagMakeEff_76.C"
#include "HH4bBtagEffBase_80_v2.C"

using namespace std;

void HHbbbbBtagEff(int a,int b){

	string st1[50]={
		/*0-11*/
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B1000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B1200.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B1400.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B1600.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B1800.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B2000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B2500.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B3000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B3500.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B4000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/B/B4500.root",
		/*12-21*/
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R1000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R1200.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R1400.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R1600.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R1800.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R2000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R2500.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R3000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R3500.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R4000.root",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/R/R4500.root",
		/*22-32*/
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/700_2/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/1000_2/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/1500_2/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/2000_1/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/QCD/2000_2/NCUGlobalTuples_",	
		
		
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/B/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/C/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/D/NCUGlobalTuples_",
		
	};
	string  fileName[50]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2",
	
	"dataB","dataC","dataD",
	};
	
	HH4bBtagEffBase_80_v2(b,b+1,st1[a],fileName[a]);	
	/*
	if(a==38)HH4bBtagEffBase_80_v2(1,400,dataPathB,"data1");	
	else if(a==39)HH4bBtagEffBase_80_v2(400,800,dataPathB,"data2");	
	else if(a==40)HH4bBtagEffBase_80_v2(800,1048,dataPathB,"data3");	
	else if(a==41)HH4bBtagEffBase_80_v2(1,348,dataPathC,"data4");	
	else if(a==42)HH4bBtagEffBase_80_v2(1,300,dataPathD,"data5");	
	else if(a==43)HH4bBtagEffBase_80_v2(300,584,dataPathD,"data6");
	else if(a>21)HH4bBtagEffBase_80_v2(1,aa[a-22],st1[a],fileName[a],"");
	else HH4bBtagEffBase_80_v2(1,2,st1[a],fileName[a],"");
	*/
}