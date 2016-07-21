//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//root feature
#include <TLegend.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <TSystem.h>

//math 
#include <cmath>
#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "untuplizer.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
//#include "HHbbbbBtagMakeEff_76.C"
#include "HH4bBtagEffBase_80.C"
//#include "HHbbbbAnalyzerBaseC_76.C"
//#include "testBtag.C"

using namespace std;

void HHbbbbBtagEff(int a){

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
		"/data7/chchen/80x/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/2000_1/NCUGlobalTuples_",
		
		
	};
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2"
	};
	
	int aa[40]={
		355,673,115,235,90,181,47,94
	};
	string dataPath="/data7/chchen/80x/JetHT/NCUGlobalTuples_";
	//cout<<"a"<<endl;
	//testBtag();
	if(a==21){
		for(int j=0;j<11;j++){
			//HHbbbbBtagMakeEff_76(1,2,st1[j],fileName[j],"");
			HH4bBtagEffBase_80(1,2,st1[j],fileName[j],"");
		}
	}
	else if(a==30){
		//HHbbbbBtagMakeEff_76(1,500,dataPath,"data1");	
		HH4bBtagEffBase_80(1,500,dataPath,"data1");	
	}
	else if(a==31){
		//HHbbbbBtagMakeEff_76(500,999,dataPath,"data2");	
		HH4bBtagEffBase_80(500,999,dataPath,"data2");	
	}
	else if(a==32){
		//HHbbbbBtagMakeEff_76(999,1498,dataPath,"data3");	
		HH4bBtagEffBase_80(999,1498,dataPath,"data3");	
	}
	else if(a==33){
		//HHbbbbBtagMakeEff_76(1498,1997,dataPath,"data4");	
		HH4bBtagEffBase_80(1498,1997,dataPath,"data4");	
	}
	else if(a==34){
		//HHbbbbBtagMakeEff_76(1997,2026,dataPath,"data5");	
		HH4bBtagEffBase_80(1997,2026,dataPath,"data5");	
	}
	else if (a>21){
		//HHbbbbBtagMakeEff_76(1,aa[a-22],st1[a],fileName[a],"");
		HH4bBtagEffBase_80(1,aa[a-22],st1[a],fileName[a],"");
	}

}