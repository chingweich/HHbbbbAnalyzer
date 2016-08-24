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
#include "untuplizer.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "HH4bBtagMakeFastBase_80.C"

using namespace std;

void HH4bBtagMakeFast_80(int a){

	string st1[40]={
		/*0-11*/
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*12-21*/
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
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
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2",
	"bGen700","bGen1000","bGen1500","bGen2000",
	"bEnriched700","bEnriched1000","bEnriched1500","bEnriched2000",
	};
	string dataPathB="/data7/chchen/80x_dsv/JetHT_runB/NCUGlobalTuples_";
	string dataPathC="/data7/chchen/80x_dsv/JetHT_runC/NCUGlobalTuples_";
	string dataPathD="/data7/chchen/80x_dsv/JetHT_runD/NCUGlobalTuples_";
	int aa[40]={119,225,40,80,31,62,17,33,27,12,5,3,9,3,2,3};
	if(a==21)for(int j=0;j<11;j++)HH4bBtagMakeFastBase_80(1,2,st1[j],fileName[j],"");
	else if(a==38)HH4bBtagMakeFastBase_80(1,400,dataPathB,"data1");	
	else if(a==39)HH4bBtagMakeFastBase_80(400,800,dataPathB,"data2");	
	else if(a==40)HH4bBtagMakeFastBase_80(800,1048,dataPathB,"data3");	
	else if(a==41)HH4bBtagMakeFastBase_80(1,348,dataPathC,"data4");	
	else if(a==42)HH4bBtagMakeFastBase_80(1,300,dataPathD,"data5");	
	else if(a==43)HH4bBtagMakeFastBase_80(300,584,dataPathD,"data6");
	else if (a>21)HH4bBtagMakeFastBase_80(1,aa[a-22],st1[a],fileName[a],"");
}