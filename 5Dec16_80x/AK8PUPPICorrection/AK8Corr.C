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

//math 
#include <cmath>
//#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../../untuplizer.h"
//#include "jetEnergyScale.h"

//#include "standalone_LumiReWeighting.cc"
//#include "standalone_LumiReWeighting.h"

//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "AK8CorrBase.C"

using namespace std;

void AK8Corr(){

	string st1[50]={
		/*0-11*/
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B400.root",
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B500.root",
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B600.root",
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B700.root",
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B800.root",
		"/data7/syu/NCUGlobalTuples/80X_MiniAODv2/80X_puppi/a64c3c6/BulkGravTohhTo4b/B900.root",
		
		
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1200.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1400.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1600.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1800.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B2000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B2500.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B3000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B3500.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B4000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B4500.root",
		/*12-21*/
		"/data7/chchen/80x_dsv/R/R1000.root",
		"/data7/chchen/80x_dsv/R/R1200.root",
		"/data7/chchen/80x_dsv/R/R1400.root",
		"/data7/chchen/80x_dsv/R/R1600.root",
		"/data7/chchen/80x_dsv/R/R1800.root",
		"/data7/chchen/80x_dsv/R/R2000.root",
		"/data7/chchen/80x_dsv/R/R2500.root",
		"/data7/chchen/80x_dsv/R/R3000.root",
		"/data7/chchen/80x_dsv/R/R3500.root",
		"/data7/chchen/80x_dsv/R/R4000.root",
		"/data7/chchen/80x_dsv/R/R4500.root",
		/*22-32*/
		"/data7/chchen/80x_dsv/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/2000_1/NCUGlobalTuples_",	
		
		"/data7/chchen/80x_dsv/QCDHTbGen/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/2000/NCUGlobalTuples_",
		
		"/data7/chchen/80x_dsv/QCDHTbEnriched/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/2000/NCUGlobalTuples_",
	};
	string  fileName[50]={
	"B400","B500","B600","B700","B800","B900","B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2",
	"bGen700","bGen1000","bGen1500","bGen2000",
	"bEnriched700","bEnriched1000","bEnriched1500","bEnriched2000",
	};
	int aa[40]={119,225,40,80,31,62,17,33,27,12,5,3,9,3,2,3};
	string dataPathB="/data7/chchen/80x_dsv/JetHT_runB/NCUGlobalTuples_";
	string dataPathC="/data7/chchen/80x_dsv/JetHT_runC/NCUGlobalTuples_";
	string dataPathD="/data7/chchen/80x_dsv/JetHT_runD/NCUGlobalTuples_";
	
	for(int i=0;i<17;i++)AK8CorrBase(1,2,st1[i],fileName[i]);
		/*
	a=0;
	if(a==38)HH4bCategoryBase_80(1,400,dataPathB,"data1");	
	else if(a==39)HH4bCategoryBase_80(400,800,dataPathB,"data2");	
	else if(a==40)HH4bCategoryBase_80(800,1048,dataPathB,"data3");	
	else if(a==41)HH4bCategoryBase_80(1,348,dataPathC,"data4");	
	else if(a==42)HH4bCategoryBase_80(1,300,dataPathD,"data5");	
	else if(a==43)HH4bCategoryBase_80(300,584,dataPathD,"data6");
	else if (a>21)HH4bCategoryBase_80(1,aa[a-22],st1[a],fileName[a],"");
	else HH4bCategoryBase_80(1,2,st1[a],fileName[a],"");
	*/
}